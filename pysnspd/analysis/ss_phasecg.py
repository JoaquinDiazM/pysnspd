"""Analysis products for the corrected Allmaras phase-drive SS diagnostics."""

from __future__ import annotations

from pathlib import Path
from typing import Any, Mapping

import numpy as np

from pysnspd.analysis.ss_run import SSRunData, build_ss_plot_dataset
from pysnspd.mesh.operators import (
    build_fv_operators,
    edge_scalar_to_node_vector_least_squares,
)
from pysnspd.solver.targets import (
    dynamic_stationarity_diagnostics,
    stationarity_diagnostics,
)
from pysnspd.thermal.evolution import thermal_stationarity_diagnostics


# Current production photon-readiness policy.  Old SS summaries are not trusted
# for this visualization because several predate the dynamic-attractor and
# field-tail thermal criteria.  Reanalysis uses only persisted fields.
_REANALYSIS_EVALUATION_INTERVAL_PS = 0.5
_REANALYSIS_TAIL_PS = 5.0
_REANALYSIS_DYNAMIC_PROFILE_REL = 1.0e-2
_REANALYSIS_DYNAMIC_VOLTAGE_REL = 2.0e-2
_REANALYSIS_THERMAL_REL = 3.0e-3
_REANALYSIS_THERMAL_P99_K = 3.0e-3
_REANALYSIS_THERMAL_PROJECTION_PS = 20.0
_REANALYSIS_THERMAL_PROJECTION_REL = 1.0e-2
_REANALYSIS_CONSECUTIVE = 5


def build_phasecg_diagnostic_dataset(
    run: SSRunData,
    *,
    thickness_m: float,
    center_width_m: float = 100.0e-9,
    measured_wall_time_s: float | None = None,
) -> dict[str, Any]:
    """Build physical, snapshot, and numerical diagnostics for one SS run.

    The existing solver does not store wall-clock time for every accepted step.
    When a measured total wall time is supplied, the per-step curve is therefore
    an explicitly labelled estimate proportional to the number of nonlinear
    solve attempts in that accepted step.  Its cumulative integral is scaled to
    the measured total.
    """

    base = build_ss_plot_dataset(run)
    history = run.history
    snapshots = _snapshot_source(history, run.raw_ss / "stationary_snapshots.npz")
    summary = run.summary
    solver = _mapping(summary.get("solver"))

    nodes_m = np.asarray(run.mesh.nodes, dtype=float)[:, :2]
    x_m = nodes_m[:, 0]
    n_nodes = x_m.size
    triangles = np.asarray(run.mesh.triangles, dtype=np.int64)
    ops = build_fv_operators(run.mesh, run.edge_data)

    snapshot_t_s = _first_array(
        snapshots,
        ("snapshot_t_s", "delta_snapshot_t_s", "phi_snapshot_t_s"),
    )
    snapshot_t_ps = np.asarray(snapshot_t_s, dtype=float).reshape(-1) / 1.0e-12
    n_snapshots = snapshot_t_ps.size
    if n_snapshots == 0:
        raise ValueError("The SS run does not contain stored snapshots.")

    delta0_meV = _scalar(snapshots.get("delta0_meV"), default=base.get("delta0_meV", np.nan))
    javg_A_m2 = abs(_scalar(snapshots.get("javg_A_m2"), default=base.get("javg_A_m2", np.nan)))
    if not np.isfinite(delta0_meV) or delta0_meV <= 0.0:
        raise ValueError("A positive delta0_meV is required for phase-CG diagnostics.")
    if not np.isfinite(javg_A_m2) or javg_A_m2 <= 0.0:
        raise ValueError("A positive javg_A_m2 is required for current normalization.")

    width_m = float(getattr(run.mesh, "width_m", np.ptp(nodes_m[:, 1])))
    thickness = float(thickness_m)
    if not np.isfinite(thickness) or thickness <= 0.0:
        raise ValueError("thickness_m must be positive and finite.")
    cross_section_m2 = width_m * thickness

    target_current_A = float(
        summary.get("target_current_uA", np.nan) * 1.0e-6
        if "target_current_uA" in summary
        else solver.get("target_current_A", np.nan)
    )
    if not np.isfinite(target_current_A):
        target_current_A = javg_A_m2 * cross_section_m2

    delta_meV = _snapshot_matrix(snapshots, ("delta_snapshot_meV",), n_snapshots, n_nodes)
    delta_over_delta0 = delta_meV / delta0_meV
    phi_V = _snapshot_matrix(snapshots, ("phi_snapshot_V",), n_snapshots, n_nodes)
    Te_K = _snapshot_matrix(snapshots, ("Te_snapshot_K",), n_snapshots, n_nodes)
    Tph_K = _snapshot_matrix(snapshots, ("Tph_snapshot_K",), n_snapshots, n_nodes)

    jtot_x = _snapshot_matrix(
        snapshots,
        ("jtot_snapshot_x_A_m2", "current_density_snapshot_x_A_m2"),
        n_snapshots,
        n_nodes,
    )
    jtot_y = _snapshot_matrix(
        snapshots,
        ("jtot_snapshot_y_A_m2", "current_density_snapshot_y_A_m2"),
        n_snapshots,
        n_nodes,
    )
    js_x = _snapshot_matrix(
        snapshots,
        ("js_us_snapshot_x_A_m2", "supercurrent_Usadel_density_snapshot_x_A_m2"),
        n_snapshots,
        n_nodes,
    )
    js_y = _snapshot_matrix(
        snapshots,
        ("js_us_snapshot_y_A_m2", "supercurrent_Usadel_density_snapshot_y_A_m2"),
        n_snapshots,
        n_nodes,
    )
    jn_x = _snapshot_matrix(
        snapshots,
        ("jn_snapshot_x_A_m2", "normal_current_density_snapshot_x_A_m2"),
        n_snapshots,
        n_nodes,
    )
    jn_y = _snapshot_matrix(
        snapshots,
        ("jn_snapshot_y_A_m2", "normal_current_density_snapshot_y_A_m2"),
        n_snapshots,
        n_nodes,
    )

    jtot_mag_over_javg = np.hypot(jtot_x, jtot_y) / javg_A_m2
    js_mag_over_javg = np.hypot(js_x, js_y) / javg_A_m2
    jn_mag_over_javg = np.hypot(jn_x, jn_y) / javg_A_m2

    edge_i = np.asarray(snapshots.get("edge_i", []), dtype=np.int64).reshape(-1)
    edge_j = np.asarray(snapshots.get("edge_j", []), dtype=np.int64).reshape(-1)
    edge_q = np.asarray(
        _first_array(snapshots, ("edge_Q_snapshot_m_inv", "edge_phase_gradient_snapshot_m_inv")),
        dtype=float,
    )
    if edge_q.shape[0] != n_snapshots:
        edge_q = np.resize(edge_q, (n_snapshots, edge_i.size))
    q_node_m_inv = np.zeros((n_snapshots, n_nodes), dtype=float)
    for snapshot_index, q_projection in enumerate(edge_q):
        qx, qy = edge_scalar_to_node_vector_least_squares(q_projection, ops)
        q_node_m_inv[snapshot_index] = np.hypot(qx, qy)

    xi_m = _scalar(snapshots.get("stationarity_xi_m"), default=np.nan)
    if not np.isfinite(xi_m) or xi_m <= 0.0:
        xi_m = float(_mapping(solver.get("contact_recovery")).get("xi_m", np.nan))
    if not np.isfinite(xi_m) or xi_m <= 0.0:
        raise ValueError("The run does not provide a positive coherence length.")

    div_j_A_m3 = _snapshot_matrix(
        snapshots,
        ("div_jtot_snapshot_A_m3", "divergence_snapshot_A_m3"),
        n_snapshots,
        n_nodes,
    )
    div_j_normalized = xi_m * div_j_A_m3 / javg_A_m2

    phase_drive = _snapshot_matrix(
        snapshots,
        ("allmaras_phase_drive_abs_over_delta0_snapshot",),
        n_snapshots,
        n_nodes,
    )

    xmin = float(np.nanmin(x_m))
    xmax = float(np.nanmax(x_m))
    x_center = 0.5 * (xmin + xmax)
    half_width = 0.5 * float(center_width_m)
    center_mask = np.abs(x_m - x_center) <= half_width
    if not np.any(center_mask):
        center_mask = np.ones(n_nodes, dtype=bool)

    current_sign = _current_orientation(jtot_x[:, center_mask], target_current_A)
    jtot_x_over_javg = current_sign * jtot_x / javg_A_m2
    js_x_over_javg = current_sign * js_x / javg_A_m2
    jn_x_over_javg = current_sign * jn_x / javg_A_m2

    terminal_left = _nearest_x_column_mask(x_m, xmin)
    terminal_right = _nearest_x_column_mask(x_m, xmax)
    probe_left = _nearest_x_column_mask(x_m, x_center - half_width)
    probe_right = _nearest_x_column_mask(x_m, x_center + half_width)

    center_weights = np.asarray(ops.node_area_m2, dtype=float)[center_mask]
    current_total_A = current_sign * _weighted_rows(jtot_x[:, center_mask], center_weights) * cross_section_m2
    current_super_A = current_sign * _weighted_rows(js_x[:, center_mask], center_weights) * cross_section_m2
    current_normal_A = current_sign * _weighted_rows(jn_x[:, center_mask], center_weights) * cross_section_m2

    voltage_terminal_V = np.abs(
        np.nanmean(phi_V[:, terminal_right], axis=1)
        - np.nanmean(phi_V[:, terminal_left], axis=1)
    )
    voltage_center_V = np.abs(
        np.nanmean(phi_V[:, probe_right], axis=1)
        - np.nanmean(phi_V[:, probe_left], axis=1)
    )

    central_delta = delta_over_delta0[:, center_mask]
    delta_center_min = np.nanmin(central_delta, axis=1)
    delta_center_mean = _weighted_rows(central_delta, center_weights)
    delta_center_max = np.nanmax(central_delta, axis=1)
    normal_fraction = np.abs(current_normal_A) / np.maximum(np.abs(current_total_A), 1.0e-300)

    bulk_mask = np.asarray(history.get("allmaras_bulk_node_mask", center_mask), dtype=bool).reshape(-1)
    if bulk_mask.size != n_nodes:
        bulk_mask = np.resize(bulk_mask, n_nodes).astype(bool)
    div_bulk = np.abs(div_j_normalized[:, bulk_mask])
    div_normalized_max = np.nanmax(div_bulk, axis=1)
    div_normalized_rms = np.sqrt(np.nanmean(div_bulk**2, axis=1))

    history_t_ps = np.asarray(base.get("t_ps", []), dtype=float)
    rejected_per_step = _history_series(history, "adaptive_rejected_attempts", history_t_ps.size)
    solve_attempts_per_step = 1.0 + np.maximum(rejected_per_step, 0.0)
    cumulative_rejected_attempts = np.cumsum(np.maximum(rejected_per_step, 0.0))
    cumulative_solve_attempts = np.cumsum(solve_attempts_per_step)

    wall_total = None if measured_wall_time_s is None else float(measured_wall_time_s)
    estimated_wall_step_s = np.array([], dtype=float)
    estimated_wall_cumulative_s = np.array([], dtype=float)
    if wall_total is not None:
        if not np.isfinite(wall_total) or wall_total <= 0.0:
            raise ValueError("measured_wall_time_s must be positive and finite when supplied.")
        estimated_wall_step_s = wall_total * solve_attempts_per_step / np.sum(solve_attempts_per_step)
        estimated_wall_cumulative_s = np.cumsum(estimated_wall_step_s)

    out = dict(base)
    out.update(
        {
            "snapshot_t_ps": snapshot_t_ps,
            "nodes_x_nm": x_m * 1.0e9,
            "nodes_y_nm": nodes_m[:, 1] * 1.0e9,
            "triangles": triangles,
            "delta0_meV": delta0_meV,
            "delta_snapshot_over_delta0": delta_over_delta0,
            "phi_snapshot_mV": phi_V * 1.0e3,
            "Te_snapshot_K": Te_K,
            "Tph_snapshot_K": Tph_K,
            "qxi_snapshot": q_node_m_inv * xi_m,
            "jtot_snapshot_over_javg": jtot_mag_over_javg,
            "js_snapshot_over_javg": js_mag_over_javg,
            "jn_snapshot_over_javg": jn_mag_over_javg,
            "jtot_x_snapshot_over_javg": jtot_x_over_javg,
            "js_x_snapshot_over_javg": js_x_over_javg,
            "jn_x_snapshot_over_javg": jn_x_over_javg,
            "jtot_y_snapshot_over_javg": current_sign * jtot_y / javg_A_m2,
            "js_y_snapshot_over_javg": current_sign * js_y / javg_A_m2,
            "jn_y_snapshot_over_javg": current_sign * jn_y / javg_A_m2,
            "node_area_m2": np.asarray(ops.node_area_m2, dtype=float),
            "div_j_snapshot_normalized": div_j_normalized,
            "phase_drive_snapshot_over_delta0": phase_drive,
            "target_current_uA": target_current_A * 1.0e6,
            "current_total_snapshot_uA": current_total_A * 1.0e6,
            "current_super_snapshot_uA": current_super_A * 1.0e6,
            "current_normal_snapshot_uA": current_normal_A * 1.0e6,
            "voltage_terminal_snapshot_mV": voltage_terminal_V * 1.0e3,
            "voltage_center_snapshot_mV": voltage_center_V * 1.0e3,
            "delta_center_min": delta_center_min,
            "delta_center_mean": delta_center_mean,
            "delta_center_max": delta_center_max,
            "normal_current_fraction_snapshot": normal_fraction,
            "div_j_normalized_max_snapshot": div_normalized_max,
            "div_j_normalized_rms_snapshot": div_normalized_rms,
            "xi_m": xi_m,
            "cross_section_area_m2": cross_section_m2,
            "center_width_nm": float(center_width_m) * 1.0e9,
            "solve_attempts_per_step": solve_attempts_per_step,
            "cumulative_rejected_attempts": cumulative_rejected_attempts,
            "cumulative_solve_attempts": cumulative_solve_attempts,
            "measured_wall_time_s": wall_total,
            "estimated_wall_step_s": estimated_wall_step_s,
            "estimated_wall_cumulative_s": estimated_wall_cumulative_s,
            "allmaras_phase_convergence_converged": _history_series(
                history,
                "allmaras_phase_convergence_converged",
                history_t_ps.size,
            ).astype(bool),
            "allmaras_phase_convergence_iterations": _history_series(
                history,
                "allmaras_phase_convergence_iterations",
                history_t_ps.size,
            ),
            "allmaras_phase_convergence_residual_rel": _history_series(
                history,
                "allmaras_phase_convergence_residual_rel",
                history_t_ps.size,
            ),
            "allmaras_phase_continued_node_count": _history_series(
                history,
                "allmaras_phase_continued_node_count",
                history_t_ps.size,
            ),
            "allmaras_phase_direct_node_count": _history_series(
                history,
                "allmaras_phase_direct_node_count",
                history_t_ps.size,
            ),
            "allmaras_phase_zero_amplitude_node_count": _history_series(
                history,
                "allmaras_phase_zero_amplitude_node_count",
                history_t_ps.size,
            ),
            "poisson_residual_rel": _history_series(
                history,
                "pytdgl_like_poisson_residual_rel",
                history_t_ps.size,
            ),
            "allmaras_update_forcing_max_abs": _history_series(
                history,
                "allmaras_update_forcing_max_abs",
                history_t_ps.size,
            ),
            "allmaras_phase_drive_rms_snapshot": np.asarray(
                history.get("allmaras_phase_drive_rms_over_delta0", []),
                dtype=float,
            ),
            "allmaras_phase_drive_max_snapshot": np.asarray(
                history.get("allmaras_phase_drive_max_over_delta0", []),
                dtype=float,
            ),
            "usadel_vs_gl_relative_l2_snapshot": np.asarray(
                history.get("usadel_vs_gl_edge_relative_l2", []),
                dtype=float,
            ),
            "phase_convergence_tolerance": float(
                _mapping(solver.get("allmaras_phase_continuation")).get("tolerance", np.nan)
            ),
            "poisson_tolerance": float(
                _mapping(solver.get("continuity")).get("tolerance_poisson", np.nan)
            ),
            "stationarity_passes": bool(_mapping(solver.get("stationarity")).get("passes", False)),
            "dynamic_stationarity_passes": bool(
                _mapping(solver.get("dynamic_stationarity")).get("passes", False)
            ),
            "thermal_stationarity_passes": bool(
                _mapping(solver.get("thermal_stationarity")).get("passes", False)
            ),
            "thermal_enabled": bool(
                np.any(_history_series(history, "thermal_enabled", history_t_ps.size) > 0.5)
            ),
            "continuity_passes": bool(_mapping(solver.get("continuity")).get("passes", False)),
            "stored_photon_ready": solver.get("photon_ready"),
            "stored_stationarity_summary": dict(_mapping(solver.get("stationarity"))),
            "stored_dynamic_stationarity_summary": dict(
                _mapping(solver.get("dynamic_stationarity"))
            ),
            "stored_thermal_stationarity_summary": dict(
                _mapping(solver.get("thermal_stationarity"))
            ),
            "stored_circuit_stationarity_summary": dict(
                _mapping(solver.get("circuit_stationarity"))
            ),
        }
    )
    out.update(
        _build_stationarity_reanalysis(
            snapshots=snapshots,
            history=history,
            solver=solver,
            nodes_m=nodes_m,
            delta0_meV=delta0_meV,
            snapshot_t_ps=snapshot_t_ps,
            delta_over_delta0=delta_over_delta0,
            phi_V=phi_V,
            Te_K=Te_K,
            Tph_K=Tph_K,
            terminal_voltage_V=voltage_terminal_V,
            div_normalized_max=div_normalized_max,
            div_normalized_rms=div_normalized_rms,
            history_t_ps=history_t_ps,
        )
    )
    return out


def _snapshot_source(
    history: Mapping[str, Any],
    snapshot_path: str | Path,
) -> Mapping[str, Any]:
    """Reuse snapshot arrays already present in relaxation history.

    Production histories contain the complete snapshot bundle.  Loading the
    separate multi-gigabyte convenience NPZ as well would double peak memory.
    """

    required = (
        "snapshot_t_s",
        "delta_snapshot_meV",
        "phi_snapshot_V",
        "Te_snapshot_K",
        "Tph_snapshot_K",
    )
    if all(key in history for key in required):
        return history
    return _load_npz(snapshot_path)


def _build_stationarity_reanalysis(
    *,
    snapshots: Mapping[str, Any],
    history: Mapping[str, Any],
    solver: Mapping[str, Any],
    nodes_m: np.ndarray,
    delta0_meV: float,
    snapshot_t_ps: np.ndarray,
    delta_over_delta0: np.ndarray,
    phi_V: np.ndarray,
    Te_K: np.ndarray,
    Tph_K: np.ndarray,
    terminal_voltage_V: np.ndarray,
    div_normalized_max: np.ndarray,
    div_normalized_rms: np.ndarray,
    history_t_ps: np.ndarray,
) -> dict[str, Any]:
    """Re-evaluate current photon-readiness criteria from persisted fields."""

    eval_indices = _evaluation_snapshot_indices(
        snapshot_t_ps,
        interval_ps=_REANALYSIS_EVALUATION_INTERVAL_PS,
    )
    eval_t = np.asarray(snapshot_t_ps, dtype=float)[eval_indices]
    n_eval = eval_t.size
    nan_series = lambda: np.full(n_eval, np.nan, dtype=float)

    strict_q_margin = nan_series()
    strict_phi_margin = nan_series()
    strict_pass = np.zeros(n_eval, dtype=bool)
    stationarity_summary = _mapping(solver.get("stationarity"))
    q_rel_tol = float(stationarity_summary.get("tolerance_phase_gradient_rel", 3.0e-1))
    phi_rel_tol = float(stationarity_summary.get("tolerance_phi_gradient_rel", 2.5e-1))
    q_abs_tol = float(stationarity_summary.get("tolerance_phase_gradient_abs_m_inv", 6.0e6))
    phi_abs_tol = float(stationarity_summary.get("tolerance_phi_gradient_abs_V_m", 2.0e3))
    active_threshold = float(stationarity_summary.get("edge_active_threshold_over_bulk", 0.05))
    bulk_exclusion_xi = float(stationarity_summary.get("bulk_exclusion_xi", 4.0))

    edge_q = np.asarray(
        _first_array(
            snapshots,
            ("edge_phase_gradient_snapshot_m_inv", "edge_Q_snapshot_m_inv"),
        ),
        dtype=float,
    )
    edge_delta = np.asarray(
        snapshots.get("edge_delta_amp_over_delta0_snapshot", []),
        dtype=float,
    )
    edge_i = np.asarray(snapshots.get("edge_i", history.get("edge_i", [])), dtype=np.int64)
    edge_j = np.asarray(snapshots.get("edge_j", history.get("edge_j", [])), dtype=np.int64)
    edge_length = np.asarray(
        snapshots.get("edge_length_m", history.get("edge_length_m", [])),
        dtype=float,
    )
    strict_available = bool(
        edge_q.ndim == 2
        and edge_q.shape[0] == snapshot_t_ps.size
        and edge_i.size
        and edge_j.size == edge_i.size
        and edge_length.size == edge_i.size
    )
    static_stationarity = {
        "edge_i": edge_i,
        "edge_j": edge_j,
        "edge_length_m": edge_length,
        "edge_distance_from_contact_m": np.asarray(
            history.get("edge_distance_from_contact_m", []), dtype=float
        ),
        "stationarity_xi_m": np.asarray(
            snapshots.get("stationarity_xi_m", history.get("stationarity_xi_m", [])),
            dtype=float,
        ),
        "normal_terminal_edge_mask": np.asarray(
            history.get("normal_terminal_edge_mask", []), dtype=bool
        ),
    }
    for position in range(1, n_eval):
        if not strict_available:
            break
        pair = eval_indices[position - 1 : position + 1]
        strict_history = dict(static_stationarity)
        strict_history.update(
            {
                "edge_phase_gradient_snapshot_m_inv": edge_q[pair],
                "phi_snapshot_V": phi_V[pair],
                "eta_R": np.array([np.nan, np.nan]),
            }
        )
        if edge_delta.ndim == 2 and edge_delta.shape[0] == snapshot_t_ps.size:
            strict_history["edge_delta_amp_over_delta0_snapshot"] = edge_delta[pair]
        diagnostic = stationarity_diagnostics(
            history=strict_history,
            material=None,
            phase_gradient_rel_tol=q_rel_tol,
            phi_gradient_rel_tol=phi_rel_tol,
            phase_gradient_abs_tol_m_inv=q_abs_tol,
            phi_gradient_abs_tol_V_m=phi_abs_tol,
            edge_active_threshold=active_threshold,
            bulk_exclusion_xi=bulk_exclusion_xi,
        )
        strict_q_margin[position] = min(
            diagnostic.phase_gradient_rel_change / max(q_rel_tol, 1.0e-300),
            diagnostic.phase_gradient_abs_change_m_inv / max(q_abs_tol, 1.0e-300),
        )
        strict_phi_margin[position] = min(
            diagnostic.phi_gradient_rel_change / max(phi_rel_tol, 1.0e-300),
            diagnostic.phi_gradient_abs_change_V_m / max(phi_abs_tol, 1.0e-300),
        )
        strict_pass[position] = bool(diagnostic.passes)

    eval_delta = np.asarray(delta_over_delta0, dtype=float)[eval_indices]
    eval_Te = np.asarray(Te_K, dtype=float)[eval_indices]
    eval_Tph = np.asarray(Tph_K, dtype=float)[eval_indices]
    eval_voltage = np.asarray(terminal_voltage_V, dtype=float)[eval_indices]
    xi_values = np.asarray(
        snapshots.get("stationarity_xi_m", history.get("stationarity_xi_m", [])),
        dtype=float,
    ).reshape(-1)
    xi_m = float(xi_values[0]) if xi_values.size else float("nan")

    dynamic_profile_margin = nan_series()
    dynamic_voltage_margin = nan_series()
    dynamic_pass = np.zeros(n_eval, dtype=bool)
    thermal_relative_margin = nan_series()
    thermal_p99_margin = nan_series()
    thermal_projected_margin = nan_series()
    thermal_pass = np.zeros(n_eval, dtype=bool)
    thermal_runtime = _mapping(solver.get("thermal_runtime"))
    thermal_enabled = bool(thermal_runtime.get("enabled", False))
    thermal_start_ps = float(thermal_runtime.get("start_time_ps", 0.0))
    thermal_bath_K = float(
        thermal_runtime.get(
            "bath_K",
            np.nanmedian(eval_Tph[0]) if eval_Tph.size else 0.0,
        )
    )
    for position in range(n_eval):
        tail_start = int(
            np.searchsorted(
                eval_t,
                float(eval_t[position]) - _REANALYSIS_TAIL_PS,
                side="right",
            )
        ) - 1
        tail_start = max(0, tail_start)
        tail = slice(tail_start, position + 1)
        dynamic_history = {
            "snapshot_t_s": eval_t[tail] * 1.0e-12,
            "delta_snapshot_meV": eval_delta[tail] * float(delta0_meV),
            "stationarity_xi_m": np.asarray([xi_m]),
            "t_s": eval_t[tail] * 1.0e-12,
            "terminal_voltage_V": eval_voltage[tail],
        }
        dynamic = dynamic_stationarity_diagnostics(
            history=dynamic_history,
            nodes_m=nodes_m,
            delta0_J=float(delta0_meV) * 1.602176634e-22,
            minimum_tail_duration_ps=_REANALYSIS_TAIL_PS,
            profile_relative_tolerance=_REANALYSIS_DYNAMIC_PROFILE_REL,
            voltage_relative_tolerance=_REANALYSIS_DYNAMIC_VOLTAGE_REL,
            psl_threshold_over_delta0=0.75,
            bulk_exclusion_xi=bulk_exclusion_xi,
        )
        dynamic_profile_margin[position] = max(
            dynamic.profile_relative_fluctuation,
            dynamic.profile_relative_drift,
        ) / _REANALYSIS_DYNAMIC_PROFILE_REL
        dynamic_voltage_margin[position] = max(
            dynamic.voltage_relative_span,
            dynamic.voltage_relative_drift,
        ) / _REANALYSIS_DYNAMIC_VOLTAGE_REL
        dynamic_pass[position] = bool(dynamic.passes)

        if thermal_enabled:
            thermal = thermal_stationarity_diagnostics(
                {
                    "snapshot_t_s": eval_t[tail] * 1.0e-12,
                    "Te_snapshot_K": eval_Te[tail],
                    "Tph_snapshot_K": eval_Tph[tail],
                },
                enabled=True,
                start_time_s=thermal_start_ps * 1.0e-12,
                bath_K=thermal_bath_K,
                minimum_tail_duration_ps=_REANALYSIS_TAIL_PS,
                relative_tolerance=_REANALYSIS_THERMAL_REL,
                p99_tolerance_K=_REANALYSIS_THERMAL_P99_K,
                projection_horizon_ps=_REANALYSIS_THERMAL_PROJECTION_PS,
                projection_relative_tolerance=_REANALYSIS_THERMAL_PROJECTION_REL,
            )
            thermal_relative_margin[position] = max(
                float(thermal.get("relative_rms_drift", np.inf)),
                float(thermal.get("relative_rms_fluctuation", np.inf)),
            ) / _REANALYSIS_THERMAL_REL
            thermal_p99_margin[position] = float(
                thermal.get("p99_abs_drift_K", np.inf)
            ) / _REANALYSIS_THERMAL_P99_K
            thermal_projected_margin[position] = float(
                thermal.get("projected_relative_rms_drift", np.inf)
            ) / _REANALYSIS_THERMAL_PROJECTION_REL
            thermal_pass[position] = bool(thermal.get("passes", False))
        else:
            thermal_relative_margin[position] = 0.0
            thermal_p99_margin[position] = 0.0
            thermal_projected_margin[position] = 0.0
            thermal_pass[position] = True

    continuity_summary = _mapping(solver.get("continuity"))
    continuity_rms_tol = float(continuity_summary.get("tolerance_rms", 1.0e-6))
    continuity_max_tol = float(continuity_summary.get("tolerance_max", 1.0e-3))
    poisson_tol = float(continuity_summary.get("tolerance_poisson", 1.0e-9))
    continuity_rms_margin = np.asarray(div_normalized_rms, dtype=float)[eval_indices] / max(
        continuity_rms_tol, 1.0e-300
    )
    continuity_max_margin = np.asarray(div_normalized_max, dtype=float)[eval_indices] / max(
        continuity_max_tol, 1.0e-300
    )
    poisson_values = _sample_previous(
        history_t_ps,
        np.asarray(history.get("pytdgl_like_poisson_residual_rel", []), dtype=float),
        eval_t,
    )
    poisson_margin = poisson_values / max(poisson_tol, 1.0e-300)
    continuity_pass = (
        (continuity_rms_margin <= 1.0)
        & (continuity_max_margin <= 1.0)
        & (poisson_margin <= 1.0)
    )

    phase_summary = _mapping(solver.get("allmaras_phase_continuation"))
    phase_tol = float(phase_summary.get("tolerance", 3.0e-3))
    phase_residual = _sample_previous(
        history_t_ps,
        np.asarray(history.get("allmaras_phase_convergence_residual_rel", []), dtype=float),
        eval_t,
    )
    phase_converged = _sample_previous(
        history_t_ps,
        np.asarray(history.get("allmaras_phase_convergence_converged", []), dtype=float),
        eval_t,
    ) > 0.5
    phase_margin = phase_residual / max(phase_tol, 1.0e-300)
    phase_pass = phase_converged & np.isfinite(phase_margin) & (phase_margin <= 1.0)

    circuit_value_margin, circuit_rhs_margin, circuit_pass, circuit_start_ps, circuit_hold_ps = (
        _circuit_stationarity_series(history=history, solver=solver, eval_t_ps=eval_t)
    )
    contact_pass_value = bool(_mapping(solver.get("contact_recovery")).get("passes", False))
    contact_pass = np.full(n_eval, contact_pass_value, dtype=bool)
    mesoscopic_pass = strict_pass | dynamic_pass
    minimum_gate_time = max(circuit_start_ps, thermal_start_ps) + max(
        circuit_hold_ps,
        _REANALYSIS_TAIL_PS,
        _REANALYSIS_EVALUATION_INTERVAL_PS,
    )
    eligible = eval_t >= minimum_gate_time - 1.0e-9
    instantaneous_ready = (
        eligible
        & mesoscopic_pass
        & contact_pass
        & continuity_pass
        & thermal_pass
        & circuit_pass
        & phase_pass
    )
    photon_ready = _persistent_true(instantaneous_ready, _REANALYSIS_CONSECUTIVE)
    ready_times = eval_t[photon_ready]
    first_ready = float(ready_times[0]) if ready_times.size else None

    return {
        "stationarity_eval_t_ps": eval_t,
        "strict_q_tolerance_margin": strict_q_margin,
        "strict_phi_tolerance_margin": strict_phi_margin,
        "strict_stationarity_pass_history": strict_pass,
        "dynamic_profile_tolerance_margin": dynamic_profile_margin,
        "dynamic_voltage_tolerance_margin": dynamic_voltage_margin,
        "dynamic_stationarity_pass_history": dynamic_pass,
        "mesoscopic_stationarity_pass_history": mesoscopic_pass,
        "continuity_rms_tolerance_margin": continuity_rms_margin,
        "continuity_max_tolerance_margin": continuity_max_margin,
        "poisson_tolerance_margin": poisson_margin,
        "continuity_pass_history": continuity_pass,
        "phase_cg_tolerance_margin": phase_margin,
        "phase_cg_pass_history": phase_pass,
        "thermal_relative_tolerance_margin": thermal_relative_margin,
        "thermal_p99_tolerance_margin": thermal_p99_margin,
        "thermal_projected_tolerance_margin": thermal_projected_margin,
        "thermal_stationarity_pass_history": thermal_pass,
        "circuit_value_tolerance_margin": circuit_value_margin,
        "circuit_rhs_tolerance_margin": circuit_rhs_margin,
        "circuit_stationarity_pass_history": circuit_pass,
        "contact_recovery_pass_history": contact_pass,
        "instantaneous_photon_ready_history": instantaneous_ready,
        "photon_ready_history": photon_ready,
        "photon_ready_reanalysis_summary": {
            "policy": "current_photon_ready_from_persisted_fields_v1",
            "passes": bool(photon_ready[-1]) if photon_ready.size else False,
            "first_ready_time_ps": first_ready,
            "required_consecutive_evaluations": _REANALYSIS_CONSECUTIVE,
            "evaluation_interval_ps": _REANALYSIS_EVALUATION_INTERVAL_PS,
            "minimum_tail_duration_ps": _REANALYSIS_TAIL_PS,
            "stored_photon_ready": solver.get("photon_ready"),
        },
    }


def _evaluation_snapshot_indices(times_ps: np.ndarray, *, interval_ps: float) -> np.ndarray:
    times = np.asarray(times_ps, dtype=float).reshape(-1)
    if times.size == 0:
        return np.array([], dtype=np.int64)
    targets = np.arange(float(times[0]), float(times[-1]) + 0.5 * interval_ps, interval_ps)
    right = np.searchsorted(times, targets, side="left")
    right = np.clip(right, 0, times.size - 1)
    left = np.clip(right - 1, 0, times.size - 1)
    choose_left = np.abs(times[left] - targets) <= np.abs(times[right] - targets)
    indices = np.where(choose_left, left, right).astype(np.int64)
    indices = np.unique(np.r_[indices, times.size - 1])
    return indices


def _sample_previous(times_ps: np.ndarray, values: np.ndarray, targets_ps: np.ndarray) -> np.ndarray:
    times = np.asarray(times_ps, dtype=float).reshape(-1)
    series = np.asarray(values, dtype=float).reshape(-1)
    targets = np.asarray(targets_ps, dtype=float).reshape(-1)
    n = min(times.size, series.size)
    if n == 0:
        return np.full(targets.size, np.nan, dtype=float)
    times = times[:n]
    series = series[:n]
    indices = np.searchsorted(times, targets, side="right") - 1
    indices = np.clip(indices, 0, n - 1)
    return series[indices]


def _circuit_stationarity_series(
    *,
    history: Mapping[str, Any],
    solver: Mapping[str, Any],
    eval_t_ps: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, float, float]:
    runtime = _mapping(solver.get("circuit_runtime"))
    final_diag = _mapping(solver.get("circuit_stationarity"))
    enabled = bool(runtime.get("enabled", False))
    n_eval = np.asarray(eval_t_ps).size
    if not enabled:
        return (
            np.zeros(n_eval, dtype=float),
            np.zeros(n_eval, dtype=float),
            np.ones(n_eval, dtype=bool),
            0.0,
            0.0,
        )

    start_ps = float(runtime.get("start_time_ps", 5.0))
    hold_ps = float(final_diag.get("hold_time_ps", 5.0))
    time = np.asarray(history.get("circuit_time_s", history.get("t_s", [])), dtype=float).reshape(-1)
    time_ps = time / 1.0e-12
    active = np.asarray(history.get("circuit_active", np.ones(time.size)), dtype=bool).reshape(-1)
    component_keys = (
        "circuit_I_b_A",
        "circuit_I_s_A",
        "circuit_I_rf_A",
        "circuit_V_out_V",
        "circuit_v_c_V",
        "circuit_V_tdgl_center_V",
    )
    rhs_keys = ("circuit_dI_b_A_s", "circuit_dI_s_A_s", "circuit_dv_c_V_s")
    stored_component_tolerances = _mapping(final_diag.get("component_tolerances"))
    stored_rhs_tolerances = _mapping(final_diag.get("rhs_tolerances"))
    initial_current = abs(float(np.asarray(history.get("circuit_I_s_A", [0.0])).reshape(-1)[0]))
    initial_voltage = abs(
        float(np.asarray(history.get("circuit_V_tdgl_center_V", [0.0])).reshape(-1)[0])
    )
    current_tol = max(0.05e-6, 1.0e-2 * initial_current)
    voltage_tol = max(10.0e-6, 1.0e-2 * initial_voltage)
    component_tolerances = {
        key: float(
            stored_component_tolerances.get(
                key,
                current_tol if "_I_" in key else voltage_tol,
            )
        )
        for key in component_keys
    }
    hold_s = max(hold_ps * 1.0e-12, 1.0e-300)
    rhs_tolerances = {
        key: float(
            stored_rhs_tolerances.get(
                key,
                (current_tol if "dI_" in key else voltage_tol) / hold_s,
            )
        )
        for key in rhs_keys
    }
    value_margin = np.full(n_eval, np.inf, dtype=float)
    rhs_margin = np.full(n_eval, np.inf, dtype=float)
    passes = np.zeros(n_eval, dtype=bool)
    for position, evaluation_time in enumerate(np.asarray(eval_t_ps, dtype=float)):
        stop = int(np.searchsorted(time_ps, evaluation_time, side="right"))
        start = int(np.searchsorted(time_ps, evaluation_time - hold_ps, side="left"))
        if stop - start < 2 or evaluation_time < start_ps + hold_ps:
            continue
        tail_active = active[start:stop]
        if not np.any(tail_active):
            continue
        duration = time_ps[start:stop][tail_active]
        if duration.size < 2 or duration[-1] - duration[0] < 0.999 * hold_ps:
            continue
        ratios = []
        for key in component_keys:
            values = np.asarray(history.get(key, []), dtype=float).reshape(-1)
            if values.size < stop:
                ratios.append(np.inf)
                continue
            tail_values = values[start:stop][tail_active]
            span = float(np.nanmax(tail_values) - np.nanmin(tail_values))
            ratios.append(span / max(component_tolerances[key], 1.0e-300))
        rhs_ratios = []
        for key in rhs_keys:
            values = np.asarray(history.get(key, []), dtype=float).reshape(-1)
            if values.size < stop:
                rhs_ratios.append(np.inf)
            else:
                rhs_ratios.append(abs(float(values[stop - 1])) / max(rhs_tolerances[key], 1.0e-300))
        value_margin[position] = max(ratios)
        rhs_margin[position] = max(rhs_ratios)
        passes[position] = bool(value_margin[position] <= 1.0 and rhs_margin[position] <= 1.0)
    return value_margin, rhs_margin, passes, start_ps, hold_ps


def _persistent_true(values: np.ndarray, required: int) -> np.ndarray:
    source = np.asarray(values, dtype=bool).reshape(-1)
    out = np.zeros(source.size, dtype=bool)
    count = 0
    latched = False
    for index, value in enumerate(source):
        count = count + 1 if value else 0
        latched = bool(latched or count >= max(1, int(required)))
        out[index] = latched
    return out


def _load_npz(path: str | Path) -> dict[str, np.ndarray]:
    source = Path(path)
    if not source.exists():
        raise FileNotFoundError(f"Missing SS snapshots: {source}")
    with np.load(source, allow_pickle=True) as data:
        return {key: np.asarray(data[key]) for key in data.files}


def _mapping(value: Any) -> Mapping[str, Any]:
    return value if isinstance(value, Mapping) else {}


def _scalar(value: Any, *, default: Any = np.nan) -> float:
    try:
        array = np.asarray(value if value is not None else default, dtype=float).reshape(-1)
        return float(array[-1]) if array.size else float(default)
    except Exception:
        return float(default)


def _first_array(data: Mapping[str, Any], keys: tuple[str, ...]) -> np.ndarray:
    for key in keys:
        if key in data:
            return np.asarray(data[key])
    return np.array([], dtype=float)


def _snapshot_matrix(
    data: Mapping[str, Any],
    keys: tuple[str, ...],
    n_snapshots: int,
    n_nodes: int,
) -> np.ndarray:
    values = np.asarray(_first_array(data, keys), dtype=float)
    if values.size == 0:
        return np.zeros((n_snapshots, n_nodes), dtype=float)
    if values.shape != (n_snapshots, n_nodes):
        values = np.resize(values, (n_snapshots, n_nodes))
    return values


def _weighted_rows(values: np.ndarray, weights: np.ndarray) -> np.ndarray:
    matrix = np.asarray(values, dtype=float)
    weight = np.asarray(weights, dtype=float).reshape(-1)
    if matrix.ndim != 2 or matrix.shape[1] != weight.size:
        raise ValueError("Weighted row average requires shape (n_rows, n_weights).")
    finite_weight = np.where(np.isfinite(weight) & (weight > 0.0), weight, 0.0)
    finite_values = np.isfinite(matrix)
    numerator = np.sum(np.where(finite_values, matrix, 0.0) * finite_weight[None, :], axis=1)
    denominator = np.sum(finite_values * finite_weight[None, :], axis=1)
    return numerator / np.maximum(denominator, 1.0e-300)


def _nearest_x_column_mask(x_m: np.ndarray, target_m: float) -> np.ndarray:
    x = np.asarray(x_m, dtype=float)
    distance = np.abs(x - float(target_m))
    minimum = float(np.nanmin(distance))
    tolerance = max(1.0e-15, 1.0e-9 * max(float(np.nanmax(np.abs(x))), 1.0e-12))
    return distance <= minimum + tolerance


def _current_orientation(jx_center: np.ndarray, target_current_A: float) -> float:
    median = float(np.nanmedian(jx_center))
    target_sign = 1.0 if target_current_A >= 0.0 else -1.0
    if not np.isfinite(median) or median == 0.0:
        return target_sign
    return target_sign * np.sign(median)


def _history_series(history: Mapping[str, Any], key: str, n: int) -> np.ndarray:
    values = np.asarray(history.get(key, []), dtype=float).reshape(-1)
    if values.size == 0:
        return np.zeros(n, dtype=float)
    return values if values.size == n else np.resize(values, n)
