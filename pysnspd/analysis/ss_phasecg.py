"""Analysis products for the corrected Allmaras phase-drive SS diagnostics."""

from __future__ import annotations

from pathlib import Path
from typing import Any, Mapping

import numpy as np

from pysnspd.analysis.ss_run import SSRunData, build_ss_plot_dataset
from pysnspd.gtdgl.operators import (
    build_fv_operators,
    edge_scalar_to_node_vector_least_squares,
)


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
    snapshots = _load_npz(run.raw_ss / "stationary_snapshots.npz")
    history = run.history
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

    terminal_left = _nearest_x_column_mask(x_m, xmin)
    terminal_right = _nearest_x_column_mask(x_m, xmax)
    probe_left = _nearest_x_column_mask(x_m, x_center - half_width)
    probe_right = _nearest_x_column_mask(x_m, x_center + half_width)

    current_sign = _current_orientation(jtot_x[:, center_mask], target_current_A)
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
            "delta_snapshot_over_delta0": delta_over_delta0,
            "phi_snapshot_mV": phi_V * 1.0e3,
            "qxi_snapshot": q_node_m_inv * xi_m,
            "jtot_snapshot_over_javg": jtot_mag_over_javg,
            "js_snapshot_over_javg": js_mag_over_javg,
            "jn_snapshot_over_javg": jn_mag_over_javg,
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
            "continuity_passes": bool(_mapping(solver.get("continuity")).get("passes", False)),
        }
    )
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
