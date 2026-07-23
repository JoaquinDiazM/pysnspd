"""Adapters between pySNSPD OE7 data and the pyTDGL-like solver core."""
from __future__ import annotations

from types import SimpleNamespace
from typing import Any

import numpy as np

from pysnspd.gtdgl.material import GTDGLMaterial
from pysnspd.mesh.operators import FVOperators, terminal_voltage, edge_scalar_to_node_vector_least_squares
from pysnspd.gtdgl.state import GTDGLStationaryState, RelaxationResult
from pysnspd.solver.diagnostics import (
    current_residual,
    current_density_maxima_A_m2,
    seed_target_current_A,
    target_current_density_A_m2,
)
from pysnspd.gtdgl.currents import native_edge_currents_to_current_fields, native_current_scale_A_m2
from pysnspd.gtdgl.usadel_current import compute_usadel_supercurrent_diagnostic
from pysnspd.gtdgl.allmaras import (
    PhaseDriveContinuationSolver,
    allmaras_coefficients,
    compute_allmaras_appendix_b_diagnostic,
    compute_allmaras_forcing_dimensionless,
    rms as _allmaras_rms,
    max_abs as _allmaras_max_abs,
)
from pysnspd.mesh.device import build_pytdgl_like_device
from .options import SolverOptions, SparseSolver
from .core import TDGLSolver
from pysnspd.thermal.evolution import ThermalRuntimeConfig, ThermalRuntimeController, thermal_stationarity_diagnostics
from .targets import (
    apply_terminal_proximity_seed,
    contact_recovery_diagnostics,
    continuity_diagnostics,
    dynamic_stationarity_diagnostics,
    stationarity_diagnostics,
)

MEV_J = 1.602176634e-22

from pysnspd.solver.callbacks import (
    _terminal_edge_mask_from_device,
    _terminal_site_mask_from_device,
)

def _build_history(
    *,
    solution,
    mesh,
    edge_data,
    ops: FVOperators,
    material: GTDGLMaterial,
    Te: np.ndarray,
    Tph: np.ndarray,
    psi_final_J: np.ndarray,
    phi_final_V: np.ndarray,
    currents,
    target_current_A: float,
    usadel_catalog: Any | None,
    n_snapshots: int,
    phase_drive_continuation: PhaseDriveContinuationSolver,
    stationarity_bulk_exclusion_xi: float = 4.0,
) -> dict[str, np.ndarray]:
    raw = solution.history
    terminal_site_mask = _terminal_site_mask_from_device(solution.device, mesh.n_nodes)
    blocked_edge_mask = _terminal_edge_mask_from_device(solution.device, ops)
    tau0 = float(material.tau0_GL_s)
    dt_dimless = np.asarray(raw.get("dt", []), dtype=float)
    t_s = np.cumsum(dt_dimless) * tau0 if dt_dimless.size else np.array([0.0])
    max_d = np.asarray(raw.get("max_d_abs_sq_psi", np.zeros_like(t_s)), dtype=float)
    if max_d.size != t_s.size:
        max_d = np.resize(max_d, t_s.size)
    mu_ptp = np.asarray(raw.get("mu_ptp", np.zeros_like(t_s)), dtype=float) * solution.device.voltage_scale_V
    if mu_ptp.size != t_s.size:
        mu_ptp = np.resize(mu_ptp, t_s.size)

    def _raw_series(name: str, default_value: float = 0.0) -> np.ndarray:
        arr = np.asarray(raw.get(name, np.full(t_s.shape, default_value)), dtype=float)
        if arr.size != t_s.size:
            arr = np.resize(arr, t_s.size)
        return arr

    dt_attempt_dimless = _raw_series("dt_attempt", np.nan)
    dt_accepted_dimless = _raw_series("dt_accepted", np.nan)
    dt_next_dimless = _raw_series("dt_next", np.nan)
    adaptive_target_dt_dimless = _raw_series("adaptive_target_dt", np.nan)
    adaptive_retries = _raw_series("adaptive_retries", 0.0)
    adaptive_rejected_attempts = _raw_series("adaptive_rejected_attempts", 0.0)
    adaptive_window_mean_d_abs_sq = _raw_series("adaptive_window_mean_d_abs_sq", np.nan)

    delta_abs = np.abs(psi_final_J)
    javg = abs(target_current_density_A_m2(material, target_current_A))
    jn_max_A_m2, jt_max_A_m2 = current_density_maxima_A_m2(currents)
    residual = float(current_residual(currents, mesh, material, target_current_A))

    hist: dict[str, np.ndarray] = {
        "t_s": t_s,
        "dt_s": dt_dimless * tau0,
        "dt_attempt_s": dt_attempt_dimless * tau0,
        "dt_accepted_s": dt_accepted_dimless * tau0,
        "dt_next_s": dt_next_dimless * tau0,
        "adaptive_target_dt_s": adaptive_target_dt_dimless * tau0,
        "adaptive_retries": adaptive_retries.astype(np.int64, copy=False),
        "adaptive_rejected_attempts": adaptive_rejected_attempts.astype(np.int64, copy=False),
        "adaptive_window_mean_d_abs_sq": adaptive_window_mean_d_abs_sq,
        "adaptive_enabled": np.array([bool(raw.get("adaptive_enabled", np.array([False], dtype=bool))[0])], dtype=bool),
        "adaptive_window": np.array([int(raw.get("adaptive_window", np.array([0], dtype=int))[0])], dtype=np.int64),
        "adaptive_time_step_multiplier": np.array([float(raw.get("adaptive_time_step_multiplier", np.array([np.nan]))[0])], dtype=float),
        "adaptive_growth_factor": np.array([float(raw.get("adaptive_growth_factor", np.array([np.nan]))[0])], dtype=float),
        "adaptive_total_rejected_attempts": np.array([int(raw.get("total_rejected_attempts", np.array([0], dtype=int))[0])], dtype=np.int64),
        "adaptive_max_retries_per_step": np.array([int(raw.get("max_adaptive_retries_per_step", np.array([0], dtype=int))[0])], dtype=np.int64),
        "eta_R": max_d,
        "max_amp2_change_rel": max_d,
        "current_residual": np.full(t_s.shape, residual),
        "pairbreaking_max": np.full(t_s.shape, float(np.nanmax(currents.node_pairbreaking_ratio))),
        "terminal_voltage_V": mu_ptp,
        "delta_min_over_delta0": np.full(t_s.shape, float(np.min(delta_abs) / material.delta0_J)),
        "delta_max_over_delta0": np.full(t_s.shape, float(np.max(delta_abs) / material.delta0_J)),
        "normal_current_max_A_m2": np.full(t_s.shape, float(jn_max_A_m2)),
        "total_current_max_A_m2": np.full(t_s.shape, float(jt_max_A_m2)),
        "delta0_meV": np.array([float(material.delta0_J / MEV_J)]),
        "javg_A_m2": np.array([javg]),
        "pytdgl_like_poisson_rhs_norm": np.asarray(raw.get("poisson_rhs_norm", np.zeros_like(t_s)), dtype=float),
        "pytdgl_like_poisson_residual_norm": np.asarray(raw.get("poisson_residual_norm", np.zeros_like(t_s)), dtype=float),
        "pytdgl_like_poisson_residual_rel": np.asarray(raw.get("poisson_residual_rel", np.zeros_like(t_s)), dtype=float),
        "pytdgl_like_poisson_residual_max_abs": np.asarray(raw.get("poisson_residual_max_abs", np.zeros_like(t_s)), dtype=float),
        "pytdgl_like_div_supercurrent_norm": np.asarray(raw.get("div_supercurrent_norm", np.zeros_like(t_s)), dtype=float),
        "pytdgl_like_boundary_rhs_norm": np.asarray(raw.get("boundary_rhs_norm", np.zeros_like(t_s)), dtype=float),
        "pytdgl_like_mu_boundary_max_abs": np.asarray(raw.get("mu_boundary_max_abs", np.zeros_like(t_s)), dtype=float),
        "allmaras_update_forcing_max_abs": np.asarray(raw.get("allmaras_update_forcing_max_abs", np.zeros_like(t_s)), dtype=float),
        "allmaras_phase_convergence_converged": _raw_series("allmaras_phase_convergence_converged", 0.0).astype(bool),
        "allmaras_phase_convergence_iterations": _raw_series("allmaras_phase_convergence_iterations", 0.0).astype(np.int64),
        "allmaras_phase_convergence_residual_rel": _raw_series("allmaras_phase_convergence_residual_rel", np.nan),
        "allmaras_phase_continued_node_count": _raw_series("allmaras_phase_continued_node_count", 0.0).astype(np.int64),
        "allmaras_phase_direct_node_count": _raw_series("allmaras_phase_direct_node_count", 0.0).astype(np.int64),
        "allmaras_phase_zero_amplitude_node_count": _raw_series("allmaras_phase_zero_amplitude_node_count", 0.0).astype(np.int64),
        "normal_terminal_node_mask": terminal_site_mask.astype(bool),
        "normal_terminal_edge_mask": blocked_edge_mask.astype(bool),
    }
    for key in (
        "thermal_enabled",
        "thermal_active",
        "thermal_active_n_nodes",
        "thermal_substeps",
        "thermal_max_Te_K",
        "thermal_mean_Te_K",
        "thermal_max_Tph_K",
        "thermal_mean_Tph_K",
        "thermal_max_abs_dTe_K",
        "thermal_max_abs_dTph_K",
        "thermal_max_rate_K_per_ps",
        "thermal_max_P_J_W_m3",
        "thermal_max_P_ep_W_m3",
        "thermal_max_P_esc_W_m3",
        "thermal_max_P_diff_W_m3",
    ):
        hist[key] = _raw_series(key, 0.0)

    # Store actual lightweight trajectory snapshots captured by TDGLSolver.
    # If the solver did not capture frames for any reason, fall back to the
    # final state so plotting remains robust.
    raw_snap_t = np.asarray(raw.get("snapshot_t", []), dtype=float).reshape(-1)
    psi_snap = np.asarray(raw.get("psi_snapshot", []), dtype=np.complex128)
    mu_snap = np.asarray(raw.get("mu_snapshot", []), dtype=float)
    native_super_snap = np.asarray(raw.get("supercurrent_snapshot", []), dtype=float)
    native_gl_super_snap = np.asarray(raw.get("gl_supercurrent_snapshot", []), dtype=float)
    native_normal_snap = np.asarray(raw.get("normal_current_snapshot", []), dtype=float)
    native_rhs_snap = np.asarray(raw.get("poisson_rhs_snapshot", []), dtype=float)
    native_lhs_snap = np.asarray(raw.get("poisson_lhs_snapshot", []), dtype=float)
    native_res_snap = np.asarray(raw.get("poisson_residual_snapshot", []), dtype=float)
    native_divs_snap = np.asarray(raw.get("div_supercurrent_snapshot", []), dtype=float)
    native_brhs_snap = np.asarray(raw.get("boundary_rhs_snapshot", []), dtype=float)
    Te_snap_raw = np.asarray(raw.get("Te_snapshot_K", []), dtype=float)
    Tph_snap_raw = np.asarray(raw.get("Tph_snapshot_K", []), dtype=float)

    if psi_snap.ndim != 2 or psi_snap.shape[1] != mesh.n_nodes or raw_snap_t.size == 0:
        ns = max(2, int(n_snapshots))
        snap_t = np.linspace(0.0, float(t_s[-1]) if t_s.size else 0.0, ns)
        psi_snap_J = np.tile(psi_final_J, (ns, 1))
        phi_snap_V = np.tile(phi_final_V, (ns, 1))
        Te_frames = np.tile(np.asarray(Te, dtype=float), (ns, 1))
        Tph_frames = np.tile(np.asarray(Tph, dtype=float), (ns, 1))
    else:
        ns = min(int(n_snapshots), raw_snap_t.size, psi_snap.shape[0])
        if ns < 2:
            ns = min(raw_snap_t.size, psi_snap.shape[0])
        snap_t = raw_snap_t[:ns] * tau0
        psi_snap_J = psi_snap[:ns] * material.delta0_J
        if mu_snap.ndim == 2 and mu_snap.shape[0] >= ns and mu_snap.shape[1] == mesh.n_nodes:
            phi_snap_V = mu_snap[:ns] * solution.device.voltage_scale_V
            phi_snap_V = phi_snap_V - np.mean(phi_snap_V, axis=1, keepdims=True)
        else:
            phi_snap_V = np.tile(phi_final_V, (ns, 1))
        if Te_snap_raw.ndim == 2 and Te_snap_raw.shape[0] >= ns and Te_snap_raw.shape[1] == mesh.n_nodes:
            Te_frames = Te_snap_raw[:ns].copy()
        else:
            Te_frames = np.tile(np.asarray(Te, dtype=float), (ns, 1))
        if Tph_snap_raw.ndim == 2 and Tph_snap_raw.shape[0] >= ns and Tph_snap_raw.shape[1] == mesh.n_nodes:
            Tph_frames = Tph_snap_raw[:ns].copy()
        else:
            Tph_frames = np.tile(np.asarray(Tph, dtype=float), (ns, 1))

    delta_mev = np.abs(psi_snap_J) / MEV_J
    def _native_snap(arr: np.ndarray, ncols: int) -> np.ndarray:
        if arr.ndim == 2 and arr.shape[0] >= ns and arr.shape[1] == ncols:
            return arr[:ns]
        return np.zeros((ns, ncols), dtype=float)

    native_super_si = _native_snap(native_super_snap, ops.n_edges)
    native_gl_si = _native_snap(native_gl_super_snap, ops.n_edges)
    if not np.any(native_gl_si) and np.any(native_super_si):
        native_gl_si = native_super_si.copy()
    native_normal_si = _native_snap(native_normal_snap, ops.n_edges)
    current_frames = []
    native_diags = []
    usadel_diags = []
    allmaras_diags = []
    for k in range(psi_snap_J.shape[0]):
        psi_dim = psi_snap_J[k] / material.delta0_J
        c, d = native_edge_currents_to_current_fields(
            psi_dimensionless=psi_dim,
            native_supercurrent=native_super_si[k],
            native_normal_current=native_normal_si[k],
            device=solution.device,
            mesh=mesh,
            edge_data=edge_data,
            ops=ops,
            material=material,
            Te_K=Te_frames[k],
            target_current_A=target_current_A,
        )
        udiag = compute_usadel_supercurrent_diagnostic(
            usadel_catalog=usadel_catalog,
            psi_dimensionless=psi_dim,
            material=material,
            Te_K=Te_frames[k],
            ops=ops,
            blocked_edge_mask=blocked_edge_mask,
        )
        adiag = compute_allmaras_appendix_b_diagnostic(
            psi_dimensionless=psi_dim,
            material=material,
            Te_K=Te_frames[k],
            ops=ops,
            terminal_node_mask=terminal_site_mask,
            blocked_edge_mask=blocked_edge_mask,
            edge_js_usadel_A_m2=(udiag.edge_js_usadel_A_m2 if udiag.available else None),
            phase_drive_continuation=phase_drive_continuation,
        )
        current_frames.append(c)
        native_diags.append(d)
        usadel_diags.append(udiag)
        allmaras_diags.append(adiag)

    jtot_x = np.vstack([c.node_jtot_x_A_m2 for c in current_frames])
    jtot_y = np.vstack([c.node_jtot_y_A_m2 for c in current_frames])
    scale = native_current_scale_A_m2(solution.device)
    edge_js_gl_from_psi = scale * native_gl_si
    gl_node_vectors = [edge_scalar_to_node_vector_least_squares(edge_js_gl_from_psi[k], ops) for k in range(edge_js_gl_from_psi.shape[0])]
    js_gl_x = np.vstack([v[0] for v in gl_node_vectors])
    js_gl_y = np.vstack([v[1] for v in gl_node_vectors])
    js_us_x = np.vstack([d.node_js_usadel_x_A_m2 for d in usadel_diags])
    js_us_y = np.vstack([d.node_js_usadel_y_A_m2 for d in usadel_diags])
    jn_x = np.vstack([c.node_jn_x_A_m2 for c in current_frames])
    jn_y = np.vstack([c.node_jn_y_A_m2 for c in current_frames])

    jtot_mag = np.sqrt(jtot_x**2 + jtot_y**2)
    js_gl_mag = np.sqrt(js_gl_x**2 + js_gl_y**2)
    js_us_mag = np.sqrt(js_us_x**2 + js_us_y**2)
    jn_mag = np.sqrt(jn_x**2 + jn_y**2)

    div = np.vstack([c.node_div_jtot_A_m3 for c in current_frames])
    pairbreaking = np.vstack([c.node_pairbreaking_ratio for c in current_frames])
    edge_q = np.vstack([c.edge_Q_m_inv for c in current_frames])
    edge_js_actual = np.vstack([c.edge_js_us_A_m2 for c in current_frames])
    edge_js_gl = edge_js_gl_from_psi
    edge_js_us = np.vstack([d.edge_js_usadel_A_m2 for d in usadel_diags])
    div_js_us = np.vstack([d.node_div_js_usadel_A_m3 for d in usadel_diags])
    edge_jn = np.vstack([c.edge_jn_A_m2 for c in current_frames])
    edge_jtot = np.vstack([c.edge_jtot_A_m2 for c in current_frames])

    allmaras_bulk_mask = allmaras_diags[-1].bulk_node_mask if allmaras_diags else np.ones(mesh.n_nodes, dtype=bool)
    allmaras_edge_js_us = np.vstack([d.edge_js_us_allmaras_A_m2 for d in allmaras_diags])
    allmaras_edge_js_gl = np.vstack([d.edge_js_gl_allmaras_A_m2 for d in allmaras_diags])
    allmaras_div_us = np.vstack([d.node_div_js_us_allmaras_A_m3 for d in allmaras_diags])
    allmaras_div_gl = np.vstack([d.node_div_js_gl_allmaras_A_m3 for d in allmaras_diags])
    allmaras_mismatch_div = np.vstack([d.node_mismatch_divergence_A_m3 for d in allmaras_diags])
    allmaras_phase_drive_abs_J = np.vstack([d.node_phase_drive_abs_J for d in allmaras_diags])
    allmaras_phase_drive_abs_over_delta0 = np.vstack([d.node_phase_drive_abs_over_delta0 for d in allmaras_diags])
    allmaras_delta_mod_over_delta0 = np.vstack([d.coefficients.delta_mod_over_delta0 for d in allmaras_diags])
    allmaras_rho_kwt = np.vstack([d.coefficients.rho_kwt for d in allmaras_diags])
    allmaras_xi_mod2_m2 = np.vstack([d.coefficients.xi_mod2_m2 for d in allmaras_diags])
    allmaras_C = np.vstack([d.coefficients.correction_C_J2_m3_A for d in allmaras_diags])

    time_aliases = {
        "snapshot_t_s": snap_t,
        "psi_snapshot_t_s": snap_t,
        "delta_snapshot_t_s": snap_t,
        "phase_snapshot_t_s": snap_t,
        "phi_snapshot_t_s": snap_t,
        "current_snapshot_t_s": snap_t,
        "jtot_snapshot_t_s": snap_t,
        "supercurrent_snapshot_t_s": snap_t,
        "supercurrent_GL_snapshot_t_s": snap_t,
        "supercurrent_Usadel_snapshot_t_s": snap_t,
        "js_gl_snapshot_t_s": snap_t,
        "js_us_snapshot_t_s": snap_t,
        "normal_current_snapshot_t_s": snap_t,
        "divergence_snapshot_t_s": snap_t,
        "div_jtot_snapshot_t_s": snap_t,
        "pairbreaking_snapshot_t_s": snap_t,
        "edge_snapshot_t_s": snap_t,
        "allmaras_snapshot_t_s": snap_t,
        "allmaras_mismatch_snapshot_t_s": snap_t,
        "allmaras_phase_drive_snapshot_t_s": snap_t,
    }
    hist.update({key: np.asarray(value) for key, value in time_aliases.items()})

    edge_i = np.asarray(ops.edge_i, dtype=np.int64)
    edge_j = np.asarray(ops.edge_j, dtype=np.int64)
    nodes_xy = np.asarray(mesh.nodes, dtype=float)[:, :2]
    edge_length_m = np.linalg.norm(
        nodes_xy[edge_j] - nodes_xy[edge_i],
        axis=1,
    )
    edge_length_m = np.maximum(edge_length_m, 1.0e-300)
    edge_center_x_m = 0.5 * (nodes_xy[edge_i, 0] + nodes_xy[edge_j, 0])
    xmin = float(np.nanmin(nodes_xy[:, 0]))
    xmax = float(np.nanmax(nodes_xy[:, 0]))
    edge_distance_from_contact_m = np.minimum(
        np.maximum(edge_center_x_m - xmin, 0.0),
        np.maximum(xmax - edge_center_x_m, 0.0),
    )
    xi2 = np.asarray(material.xi_mod_squared_m2(Te), dtype=float)
    stationarity_xi_m = float(np.sqrt(np.nanmedian(np.maximum(xi2, 1.0e-300))))
    bulk_exclusion_m = max(0.0, float(stationarity_bulk_exclusion_xi)) * max(stationarity_xi_m, 1.0e-300)
    stationarity_bulk_edge_mask = edge_distance_from_contact_m >= bulk_exclusion_m
    edge_phi_grad = (phi_snap_V[:, edge_j] - phi_snap_V[:, edge_i]) / edge_length_m[None, :]
    edge_amp_over_delta0 = 0.5 * (
        np.abs(psi_snap_J[:, edge_i]) + np.abs(psi_snap_J[:, edge_j])
    ) / max(float(material.delta0_J), 1.0e-300)
    final_edge_amp = edge_amp_over_delta0[-1]
    finite_edge_amp = final_edge_amp[np.isfinite(final_edge_amp)]
    if finite_edge_amp.size:
        bulk_edge_amp = float(np.nanpercentile(finite_edge_amp, 90.0))
        stationarity_active_edge_mask = final_edge_amp >= 0.05 * max(bulk_edge_amp, 1.0e-300)
    else:
        stationarity_active_edge_mask = np.ones(ops.n_edges, dtype=bool)
    stationarity_active_edge_mask = (
        stationarity_active_edge_mask
        & stationarity_bulk_edge_mask
        & ~blocked_edge_mask
    )

    for key, arr in {
        "psi_snapshot_real_J": np.real(psi_snap_J),
        "psi_snapshot_imag_J": np.imag(psi_snap_J),
        "delta_snapshot_meV": delta_mev,
        "phi_snapshot_V": phi_snap_V,
        "Te_snapshot_K": Te_frames,
        "Tph_snapshot_K": Tph_frames,
        "current_density_snapshot_A_m2": jtot_mag,
        "jtot_snapshot_mag_A_m2": jtot_mag,
        "current_density_snapshot_x_A_m2": jtot_x,
        "current_density_snapshot_y_A_m2": jtot_y,
        "jtot_snapshot_x_A_m2": jtot_x,
        "jtot_snapshot_y_A_m2": jtot_y,
        "supercurrent_density_snapshot_A_m2": js_gl_mag,
        "supercurrent_GL_density_snapshot_A_m2": js_gl_mag,
        "supercurrent_density_snapshot_x_A_m2": js_gl_x,
        "supercurrent_density_snapshot_y_A_m2": js_gl_y,
        "supercurrent_GL_density_snapshot_x_A_m2": js_gl_x,
        "supercurrent_GL_density_snapshot_y_A_m2": js_gl_y,
        "js_gl_snapshot_mag_A_m2": js_gl_mag,
        "js_gl_snapshot_x_A_m2": js_gl_x,
        "js_gl_snapshot_y_A_m2": js_gl_y,
        "supercurrent_Usadel_density_snapshot_A_m2": js_us_mag,
        "supercurrent_Usadel_density_snapshot_x_A_m2": js_us_x,
        "supercurrent_Usadel_density_snapshot_y_A_m2": js_us_y,
        "js_us_snapshot_mag_A_m2": js_us_mag,
        "js_us_snapshot_x_A_m2": js_us_x,
        "js_us_snapshot_y_A_m2": js_us_y,
        "normal_current_density_snapshot_A_m2": jn_mag,
        "jn_snapshot_mag_A_m2": jn_mag,
        "normal_current_density_snapshot_x_A_m2": jn_x,
        "normal_current_density_snapshot_y_A_m2": jn_y,
        "jn_snapshot_x_A_m2": jn_x,
        "jn_snapshot_y_A_m2": jn_y,
        "divergence_snapshot_A_m3": div,
        "div_jtot_snapshot_A_m3": div,
        "pairbreaking_ratio_snapshot": pairbreaking,
        "edge_Q_snapshot_m_inv": edge_q,
        "edge_phase_gradient_snapshot_m_inv": edge_q,
        "edge_phi_gradient_snapshot_V_m": edge_phi_grad,
        "edge_length_m": edge_length_m,
        "edge_center_x_m": edge_center_x_m,
        "edge_distance_from_contact_m": edge_distance_from_contact_m,
        "stationarity_xi_m": np.array([stationarity_xi_m], dtype=float),
        "stationarity_bulk_exclusion_xi": np.array([float(stationarity_bulk_exclusion_xi)], dtype=float),
        "stationarity_bulk_exclusion_m": np.array([bulk_exclusion_m], dtype=float),
        "stationarity_bulk_edge_mask": stationarity_bulk_edge_mask.astype(bool),
        "edge_delta_amp_over_delta0_snapshot": edge_amp_over_delta0,
        "stationarity_active_edge_mask": stationarity_active_edge_mask,
        "edge_js_actual_snapshot_A_m2": edge_js_actual,
        "edge_js_gl_snapshot_A_m2": edge_js_gl,
        "edge_js_us_snapshot_A_m2": edge_js_us,
        "edge_js_usadel_snapshot_A_m2": edge_js_us,
        "div_js_usadel_snapshot_A_m3": div_js_us,
        "edge_jn_snapshot_A_m2": edge_jn,
        "edge_jtot_snapshot_A_m2": edge_jtot,
        "allmaras_edge_js_us_snapshot_A_m2": allmaras_edge_js_us,
        "allmaras_edge_js_gl_snapshot_A_m2": allmaras_edge_js_gl,
        "allmaras_div_js_us_snapshot_A_m3": allmaras_div_us,
        "allmaras_div_js_gl_snapshot_A_m3": allmaras_div_gl,
        "allmaras_mismatch_divergence_snapshot_A_m3": allmaras_mismatch_div,
        "allmaras_phase_drive_abs_snapshot_J": allmaras_phase_drive_abs_J,
        "allmaras_phase_drive_abs_over_delta0_snapshot": allmaras_phase_drive_abs_over_delta0,
        "allmaras_delta_mod_over_delta0_snapshot": allmaras_delta_mod_over_delta0,
        "allmaras_rho_kwt_snapshot": allmaras_rho_kwt,
        "allmaras_xi_mod2_snapshot_m2": allmaras_xi_mod2_m2,
        "allmaras_correction_C_snapshot_J2_m3_A": allmaras_C,
        "allmaras_bulk_node_mask": allmaras_bulk_mask.astype(bool),
        "allmaras_bulk_mask_policy": np.array(["terminal_nodes_excluded_only"]),
        "edge_i": edge_i,
        "edge_j": edge_j,
    }.items():
        hist[key] = np.asarray(arr)

    # Native pyTDGL-like solver snapshots, kept in the internal operator units.
    # These are for debugging the sparse Poisson system and should not be
    # interpreted as SI currents without an explicit adapter conversion.
    hist["pytdgl_like_native_supercurrent_snapshot"] = native_super_si
    hist["pytdgl_like_native_normal_current_snapshot"] = native_normal_si
    hist["pytdgl_like_native_total_current_snapshot"] = (
        hist["pytdgl_like_native_supercurrent_snapshot"]
        + hist["pytdgl_like_native_normal_current_snapshot"]
    )
    hist["pytdgl_like_poisson_rhs_snapshot"] = _native_snap(native_rhs_snap, mesh.n_nodes)
    hist["pytdgl_like_poisson_lhs_snapshot"] = _native_snap(native_lhs_snap, mesh.n_nodes)
    hist["pytdgl_like_poisson_residual_snapshot"] = _native_snap(native_res_snap, mesh.n_nodes)
    hist["pytdgl_like_div_supercurrent_snapshot"] = _native_snap(native_divs_snap, mesh.n_nodes)
    hist["pytdgl_like_boundary_rhs_snapshot"] = _native_snap(native_brhs_snap, mesh.n_nodes)
    hist["pytdgl_like_snapshot_t_s"] = snap_t
    hist["pytdgl_like_native_si_current_scale_A_m2"] = np.array(
        [float(native_diags[-1].current_scale_A_m2) if native_diags else 0.0]
    )
    hist["pytdgl_like_native_si_residual_plus_boundary_rms_A_m3"] = np.asarray(
        [d.residual_plus_boundary_rms_A_m3 for d in native_diags], dtype=float
    )
    hist["pytdgl_like_native_si_residual_minus_boundary_rms_A_m3"] = np.asarray(
        [d.residual_minus_boundary_rms_A_m3 for d in native_diags], dtype=float
    )
    hist["pytdgl_like_native_si_residual_no_boundary_rms_A_m3"] = np.asarray(
        [d.residual_no_boundary_rms_A_m3 for d in native_diags], dtype=float
    )

    us_available = bool(usadel_diags and all(d.available for d in usadel_diags))
    hist["usadel_current_available"] = np.array([us_available], dtype=bool)
    hist["usadel_current_backend"] = np.array([usadel_diags[-1].backend if usadel_diags else "unavailable"], dtype=object)
    hist["usadel_current_reason"] = np.array([usadel_diags[-1].reason if usadel_diags else "not computed"], dtype=object)
    gl_norm = np.sqrt(np.nansum(edge_js_gl * edge_js_gl, axis=1))
    if us_available:
        diff = edge_js_us - edge_js_gl
        diff_norm = np.sqrt(np.nansum(diff * diff, axis=1))
        us_norm = np.sqrt(np.nansum(edge_js_us * edge_js_us, axis=1))
        hist["usadel_vs_gl_edge_relative_l2"] = diff_norm / np.maximum(us_norm, 1.0e-300)
        hist["usadel_vs_gl_edge_max_abs_diff_A_m2"] = np.nanmax(np.abs(diff), axis=1)
        hist["usadel_supercurrent_max_A_m2"] = np.nanmax(np.abs(edge_js_us), axis=1)
        hist["usadel_vs_gl_edge_usadel_norm_A_m2"] = us_norm
    else:
        hist["usadel_vs_gl_edge_relative_l2"] = np.full(edge_js_gl.shape[0], np.nan)
        hist["usadel_vs_gl_edge_max_abs_diff_A_m2"] = np.full(edge_js_gl.shape[0], np.nan)
        hist["usadel_supercurrent_max_A_m2"] = np.full(edge_js_gl.shape[0], np.nan)
        hist["usadel_vs_gl_edge_usadel_norm_A_m2"] = np.full(edge_js_gl.shape[0], np.nan)
    hist["gl_supercurrent_max_A_m2"] = np.nanmax(np.abs(edge_js_gl), axis=1)
    hist["usadel_vs_gl_edge_gl_norm_A_m2"] = gl_norm

    bulk_mask = np.asarray(hist["allmaras_bulk_node_mask"], dtype=bool)
    hist["allmaras_mismatch_div_rms_A_m3"] = np.asarray(
        [_allmaras_rms(row) for row in allmaras_mismatch_div],
        dtype=float,
    )
    hist["allmaras_mismatch_div_bulk_rms_A_m3"] = np.asarray(
        [_allmaras_rms(row, bulk_mask) for row in allmaras_mismatch_div],
        dtype=float,
    )
    hist["allmaras_mismatch_div_max_abs_A_m3"] = np.asarray(
        [_allmaras_max_abs(row) for row in allmaras_mismatch_div],
        dtype=float,
    )
    hist["allmaras_mismatch_div_bulk_max_abs_A_m3"] = np.asarray(
        [_allmaras_max_abs(row, bulk_mask) for row in allmaras_mismatch_div],
        dtype=float,
    )
    hist["allmaras_phase_drive_rms_over_delta0"] = np.asarray(
        [_allmaras_rms(row) for row in allmaras_phase_drive_abs_over_delta0],
        dtype=float,
    )
    hist["allmaras_phase_drive_bulk_rms_over_delta0"] = np.asarray(
        [_allmaras_rms(row, bulk_mask) for row in allmaras_phase_drive_abs_over_delta0],
        dtype=float,
    )
    hist["allmaras_phase_drive_max_over_delta0"] = np.asarray(
        [_allmaras_max_abs(row) for row in allmaras_phase_drive_abs_over_delta0],
        dtype=float,
    )
    hist["allmaras_phase_drive_bulk_max_over_delta0"] = np.asarray(
        [_allmaras_max_abs(row, bulk_mask) for row in allmaras_phase_drive_abs_over_delta0],
        dtype=float,
    )
    hist["allmaras_delta_mod_min_over_delta0"] = np.nanmin(allmaras_delta_mod_over_delta0, axis=1)
    hist["allmaras_delta_mod_median_over_delta0"] = np.nanmedian(allmaras_delta_mod_over_delta0, axis=1)
    hist["allmaras_delta_mod_max_over_delta0"] = np.nanmax(allmaras_delta_mod_over_delta0, axis=1)
    hist["allmaras_rho_kwt_min"] = np.nanmin(allmaras_rho_kwt, axis=1)
    hist["allmaras_rho_kwt_median"] = np.nanmedian(allmaras_rho_kwt, axis=1)
    hist["allmaras_rho_kwt_max"] = np.nanmax(allmaras_rho_kwt, axis=1)
    return hist
