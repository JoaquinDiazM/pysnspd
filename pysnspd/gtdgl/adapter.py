"""Adapters between pySNSPD OE7 data and the pyTDGL-like solver core."""
from __future__ import annotations

from types import SimpleNamespace
from typing import Any

import numpy as np

from pysnspd.gtdgl.material import GTDGLMaterial
from pysnspd.gtdgl.operators import FVOperators, terminal_voltage, edge_scalar_to_node_vector_least_squares
from pysnspd.gtdgl.state import GTDGLStationaryState, RelaxationResult
from pysnspd.gtdgl.diagnostics import (
    current_residual,
    current_density_maxima_A_m2,
    seed_target_current_A,
    target_current_density_A_m2,
)
from .currents import native_edge_currents_to_current_fields, native_current_scale_A_m2
from .usadel_current import compute_usadel_supercurrent_diagnostic
from .allmaras import (
    allmaras_coefficients,
    compute_allmaras_appendix_b_diagnostic,
    compute_allmaras_forcing_dimensionless,
    rms as _allmaras_rms,
    max_abs as _allmaras_max_abs,
)
from .device import build_pytdgl_like_device
from .options import SolverOptions, SparseSolver
from .solver import TDGLSolver
from .ss_targets import (
    apply_terminal_proximity_seed,
    contact_recovery_diagnostics,
    continuity_diagnostics,
    stationarity_diagnostics,
)

MEV_J = 1.602176634e-22


def solve_stationary_pytdgl_like(
    *,
    mesh,
    edge_data,
    seed,
    material: GTDGLMaterial,
    ops: FVOperators,
    steps: int | None = None,
    total_time_s: float | None = None,
    dt_s: float = 1.0e-17,
    target_current_A: float | None = None,
    usadel_catalog: Any | None = None,
    terminal_psi: complex | float | None = 0.0,
    adaptive: bool = True,
    adaptive_window: int = 6,
    max_solve_retries: int = 8,
    adaptive_time_step_multiplier: float = 0.5,
    adaptive_growth_factor: float = 1.5,
    dt_max_factor: float = 6.0,
    n_snapshots: int = 6,
    progress: bool = False,
    supercurrent_law: str = "gl",
    terminal_healing_xi: float | None = None,
    terminal_healing_fraction: float = 0.95,
    stationarity_eta: float = 1.0e-5,
    stationarity_phase_gradient_rel: float | None = None,
    stationarity_phi_gradient_rel: float | None = None,
    stationarity_q_abs_m_inv: float = 1.0e3,
    stationarity_phi_gradient_abs_V_m: float = 1.0e2,
    stationarity_edge_active_threshold: float = 0.05,
    stationarity_bulk_exclusion_xi: float = 4.0,
    # Deprecated aliases from the pre-gauge-gradient diagnostic.  If the new
    # tolerances are not provided, these values are used as compatibility
    # aliases for phase-gradient and phi-gradient relative tolerances.
    stationarity_delta_rel: float | None = None,
    stationarity_phi_rel: float | None = None,
    convergence_min_steps: int = 50,
    stop_on_convergence: bool = False,
    continuity_rms_tol: float = 1.0e-6,
    continuity_max_tol: float = 1.0e-3,
    continuity_poisson_tol: float = 1.0e-9,
    recovery_min_xi: float = 1.5,
    recovery_max_xi: float = 4.0,
    allmaras_contact_guard_layers: int = 2,
) -> RelaxationResult:
    """Run the essential pyTDGL-like stationary solver on a pySNSPD seed.

    The public adapter keeps pySNSPD inputs and outputs in SI units.  The
    pyTDGL-shaped core may use scaled coordinates/potentials internally, but
    terminal currents are passed in amperes and converted only inside the
    Poisson boundary operator.
    """

    if dt_s <= 0:
        raise ValueError("dt_s must be positive.")
    if total_time_s is None:
        if steps is None:
            steps = 2000
        if int(steps) <= 0:
            raise ValueError("steps must be positive when total_time_s is not provided.")
        total_time_s = int(steps) * float(dt_s)
    if float(total_time_s) <= 0:
        raise ValueError("total_time_s must be positive.")
    if target_current_A is None:
        target_current_A = seed_target_current_A(seed)
    target_current_A = float(target_current_A)
    supercurrent_law = _normalize_supercurrent_law(supercurrent_law)

    psi0_J = (
        np.asarray(seed.node_psi_real_J, dtype=float)
        + 1j * np.asarray(seed.node_psi_imag_J, dtype=float)
    )
    phi0_V = np.asarray(seed.node_phi_electric_V, dtype=float)
    Te = np.asarray(seed.node_Te_K, dtype=float).copy()
    Tph = np.asarray(seed.node_Tph_K, dtype=float).copy()

    device = build_pytdgl_like_device(
        mesh=mesh,
        edge_data=edge_data,
        material=material,
        ops=ops,
        Te_K=Te,
        target_current_A=target_current_A,
    )

    tau0 = float(material.tau0_GL_s)
    dt_dimless = float(dt_s) / tau0
    solve_time = float(total_time_s) / tau0
    requested_steps_equivalent = int(np.ceil(float(total_time_s) / float(dt_s)))

    options = SolverOptions(
        solve_time=solve_time,
        dt_init=dt_dimless,
        dt_max=dt_dimless if not adaptive else max(dt_dimless, float(dt_max_factor) * dt_dimless),
        adaptive=bool(adaptive),
        adaptive_window=int(adaptive_window),
        max_solve_retries=int(max_solve_retries),
        adaptive_time_step_multiplier=float(adaptive_time_step_multiplier),
        terminal_psi=terminal_psi,
        sparse_solver=SparseSolver.SUPERLU,
        include_screening=False,
    )
    # Additional pySNSPD adaptive controller knob.  SolverOptions intentionally
    # mirrors pyTDGL, so this is attached dynamically and consumed with getattr
    # inside TDGLSolver.
    options.adaptive_growth_factor = float(adaptive_growth_factor)

    # Appendix-B Allmaras/KWT coefficients inside the unchanged pyTDGL
    # algebraic update.  The local solver still uses its existing
    # ``X_i + z_i |X_i|^2 = w_i`` form; only the coefficient fields are now
    # constructed from Appendix B:
    #
    #   rho_KWT = sqrt(1 + gamma_i^2 |psi_i|^2),
    #   gamma_i = 2 Delta0 tau_sc(Te_i) / hbar,
    #   epsilon_i = Delta_mod^2(Te_i) / Delta0^2.
    #
    # Setting u=1 makes the existing ``dt/u * rho_KWT`` factor match the
    # Appendix-B tau0-scaled form.  The missing non-unit cubic denominator
    # Delta0^2/Delta_mod^2 is intentionally not introduced here because that
    # would change the local w/z algebra.
    psi0 = psi0_J / material.delta0_J
    allmaras0 = allmaras_coefficients(
        psi_dimensionless=psi0,
        material=material,
        Te_K=Te,
    )
    epsilon = np.clip(allmaras0.solver_epsilon, 0.0, 1.0)
    gamma = np.asarray(allmaras0.gamma_kwt_dimensionless, dtype=float)
    gamma = np.where(np.isfinite(gamma) & (gamma >= 0.0), gamma, 0.0)
    u = 1.0
    device.layer.u = u
    device.layer.gamma = gamma
    gamma_report = float(np.nanmedian(gamma)) if gamma.size else 0.0

    # pySNSPD policy: terminal currents are always SI amperes.  The
    # pyTDGL-like solver converts these physical currents to the internal
    # Neumann value needed by its Poisson operator using sigma_n, film
    # thickness and the Josephson voltage scale.
    terminal_currents_A = {"left": -target_current_A, "right": target_current_A}

    psi0 = psi0_J / material.delta0_J
    psi0, proximity_seed_diag = apply_terminal_proximity_seed(
        psi0,
        nodes_m=np.asarray(mesh.nodes, dtype=float)[:, :2],
        material=material,
        Te_K=Te,
        healing_target_xi=terminal_healing_xi,
        target_bulk_fraction=float(terminal_healing_fraction),
        terminal_value=0.0 if terminal_psi is None else terminal_psi,
    )
    mu0 = (phi0_V - float(np.mean(phi0_V))) / device.voltage_scale_V
    seed_solution = {
        "psi": psi0,
        "mu": mu0,
        "supercurrent": np.zeros(ops.n_edges, dtype=float),
        "normal_current": np.zeros(ops.n_edges, dtype=float),
    }

    def disorder_epsilon(r, *, t=None, vectorized=True):
        """Return the pySNSPD epsilon field using pyTDGL's vectorized callable API.

        pyTDGL detects vectorized disorder callables through the keyword-only
        default ``vectorized=True``.  Without this marker the solver evaluates
        the callable one site at a time and expects a scalar, which is not what
        the pySNSPD adapter wants here.
        """
        del t
        arr = np.asarray(r)
        if vectorized and arr.ndim >= 2:
            if arr.shape[0] != epsilon.size:
                return np.full(arr.shape[0], float(np.nanmedian(epsilon)), dtype=float)
            return np.asarray(epsilon, dtype=float).copy()
        return float(np.nanmedian(epsilon))

    supercurrent_override = None
    if supercurrent_law == "usadel_poisson":
        supercurrent_override = _build_usadel_poisson_supercurrent_override(
            usadel_catalog=usadel_catalog,
            device=device,
            material=material,
            Te_K=Te,
            ops=ops,
        )

    allmaras_forcing_callback = _build_allmaras_forcing_callback(
        usadel_catalog=usadel_catalog,
        device=device,
        material=material,
        Te_K=Te,
        ops=ops,
        blocked_edge_mask=_terminal_edge_mask_from_device(device, ops),
        require_usadel=(supercurrent_law == "usadel_poisson"),
        contact_guard_layers=int(allmaras_contact_guard_layers),
    )

    solver = TDGLSolver(
        device=device,
        options=options,
        applied_vector_potential=0.0,
        terminal_currents=terminal_currents_A,
        disorder_epsilon=disorder_epsilon,
        seed_solution=seed_solution,
        progress=progress,
        supercurrent_override=supercurrent_override,
        supercurrent_law=supercurrent_law,
        allmaras_forcing_callback=allmaras_forcing_callback,
        stop_eta=float(stationarity_eta) if stationarity_eta is not None and stationarity_eta > 0.0 else None,
        stop_min_steps=int(convergence_min_steps),
        stop_on_convergence=bool(stop_on_convergence),
    )
    solver.snapshot_count = max(2, int(n_snapshots))
    solution = solver.solve()
    if solution is None:
        raise RuntimeError("pytdgl_like solver returned None.")

    psi_final = solution.tdgl_data.psi
    mu_final = solution.tdgl_data.mu
    terminal_site_mask = _terminal_site_mask_from_device(device, mesh.n_nodes)
    psi_final_J = psi_final * material.delta0_J
    phi_final_V = mu_final * device.voltage_scale_V
    phi_final_V = phi_final_V - float(np.mean(phi_final_V))

    currents, native_diag = native_edge_currents_to_current_fields(
        psi_dimensionless=psi_final,
        native_supercurrent=solution.tdgl_data.supercurrent,
        native_normal_current=solution.tdgl_data.normal_current,
        device=device,
        mesh=mesh,
        edge_data=edge_data,
        ops=ops,
        material=material,
        Te_K=Te,
        target_current_A=target_current_A,
    )

    history = _build_history(
        solution=solution,
        mesh=mesh,
        edge_data=edge_data,
        ops=ops,
        material=material,
        Te=Te,
        Tph=Tph,
        psi_final_J=psi_final_J,
        phi_final_V=phi_final_V,
        currents=currents,
        target_current_A=target_current_A,
        usadel_catalog=usadel_catalog,
        n_snapshots=n_snapshots,
        stationarity_bulk_exclusion_xi=float(stationarity_bulk_exclusion_xi),
    )

    phase_gradient_rel_tol = (
        float(stationarity_phase_gradient_rel)
        if stationarity_phase_gradient_rel is not None
        else (float(stationarity_delta_rel) if stationarity_delta_rel is not None else 1.0e-4)
    )
    phi_gradient_rel_tol = (
        float(stationarity_phi_gradient_rel)
        if stationarity_phi_gradient_rel is not None
        else (float(stationarity_phi_rel) if stationarity_phi_rel is not None else 1.0e-4)
    )
    stationarity_diag = stationarity_diagnostics(
        history=history,
        material=material,
        phase_gradient_rel_tol=phase_gradient_rel_tol,
        phi_gradient_rel_tol=phi_gradient_rel_tol,
        phase_gradient_abs_tol_m_inv=float(stationarity_q_abs_m_inv),
        phi_gradient_abs_tol_V_m=float(stationarity_phi_gradient_abs_V_m),
        edge_active_threshold=float(stationarity_edge_active_threshold),
        bulk_exclusion_xi=float(stationarity_bulk_exclusion_xi),
        eta_tol=float(stationarity_eta),
    )
    recovery_diag = contact_recovery_diagnostics(
        psi_dimensionless=psi_final,
        nodes_m=np.asarray(mesh.nodes, dtype=float)[:, :2],
        material=material,
        Te_K=Te,
        threshold_fraction=float(terminal_healing_fraction),
        min_allowed_xi=float(recovery_min_xi),
        max_allowed_xi=float(recovery_max_xi),
        bin_width_m=float(getattr(mesh, "target_spacing_m", 0.0) or 0.0),
    )
    continuity_diag = continuity_diagnostics(
        currents=currents,
        mesh=mesh,
        material=material,
        target_current_A=target_current_A,
        history=history,
        rms_tol=float(continuity_rms_tol),
        max_tol=float(continuity_max_tol),
        poisson_tol=float(continuity_poisson_tol),
    )
    magic_ready = bool(stationarity_diag.passes and recovery_diag.passes and continuity_diag.passes)

    jn_max_A_m2, jt_max_A_m2 = current_density_maxima_A_m2(currents)
    summary: dict[str, Any] = {
        "backend": "pytdgl_like_minimal_no_screening",
        "supercurrent_law": supercurrent_law,
        "converged": bool(solution.history.get("converged", np.array([False], dtype=bool))[0]),
        "convergence_reason": str(solution.history.get("convergence_reason", np.array(["not_reported"], dtype=object))[0]),
        "stop_reason": str(solution.history.get("stop_reason", np.array(["not_reported"], dtype=object))[0]),
        "stop_on_convergence": bool(solution.history.get("stop_on_convergence", np.array([False], dtype=bool))[0]),
        "eta_converged": bool(solution.history.get("eta_converged", np.array([False], dtype=bool))[0]),
        "eta_convergence_step": int(solution.history.get("eta_convergence_step", np.array([-1], dtype=int))[0]),
        "eta_convergence_time_ps": float(solution.history.get("eta_convergence_time", np.array([float("nan")]))[0] * tau0 / 1.0e-12),
        "first_magic_ready": magic_ready,
        "accepted_steps": int(solution.history.get("final_step", np.array([0]))[0]),
        "rejected_steps": int(solution.history.get("total_rejected_attempts", np.array([0]))[0]),
        "requested_time_ps": float(total_time_s) / 1.0e-12,
        "requested_steps_equivalent_at_dt_init": int(requested_steps_equivalent),
        "final_time_ps": float(solution.history.get("final_time", np.array([0.0]))[0] * tau0 / 1.0e-12),
        "requested_time_reached": bool(
            float(solution.history.get("final_time", np.array([0.0]))[0] * tau0)
            >= 0.999999 * float(total_time_s)
        ),
        "dt_init_s": float(dt_s),
        "adaptive_enabled": bool(solution.history.get("adaptive_enabled", np.array([False], dtype=bool))[0]),
        "adaptive_window": int(solution.history.get("adaptive_window", np.array([0], dtype=int))[0]),
        "adaptive_time_step_multiplier": float(solution.history.get("adaptive_time_step_multiplier", np.array([float("nan")]))[0]),
        "adaptive_growth_factor": float(solution.history.get("adaptive_growth_factor", np.array([float("nan")]))[0]),
        "adaptive_dt_max_factor": float(dt_max_factor),
        "adaptive_max_retries_per_step": int(solution.history.get("max_adaptive_retries_per_step", np.array([0], dtype=int))[0]),
        "adaptive_total_rejected_attempts": int(solution.history.get("total_rejected_attempts", np.array([0], dtype=int))[0]),
        "dt_final_s": float(history.get("dt_s", np.array([float("nan")]))[-1]) if history.get("dt_s", np.array([])).size else float("nan"),
        "dt_min_s": float(np.nanmin(history.get("dt_s", np.array([float("nan")])))),
        "dt_max_used_s": float(np.nanmax(history.get("dt_s", np.array([float("nan")])))),
        "dt_mean_s": float(np.nanmean(history.get("dt_s", np.array([float("nan")])))),
        "tau0_GL_s": tau0,
        "material_tau_ee_Tc_s": float(material.tau_ee_Tc_s),
        "material_tau_ep_Tc_s": float(material.tau_ep_Tc_s),
        "material_tau_ee_Tc_ps": float(material.tau_ee_Tc_s / 1.0e-12),
        "material_tau_ep_Tc_ps": float(material.tau_ep_Tc_s / 1.0e-12),
        "material_tau_ee_median_s": float(np.nanmedian(material.tau_ee_s(Te))),
        "material_tau_ep_median_s": float(np.nanmedian(material.tau_ep_s(Te))),
        "material_tau_sc_median_s": float(np.nanmedian(material.tau_sc_s(Te))),
        "pytdgl_u": u,
        "pytdgl_gamma": gamma_report,
        "allmaras_coefficients_backend": "appendix_b_allmaras_wz_update_v1",
        "allmaras_update_backend": "appendix_b_explicit_forcing_rho_kwt_wz_v1_contact_guarded",
        "allmaras_contact_correction_guard_layers": int(allmaras_contact_guard_layers),
        "stationarity_bulk_exclusion_xi": float(stationarity_bulk_exclusion_xi),
        "allmaras_solver_u": float(u),
        "allmaras_gamma_min": float(np.nanmin(gamma)) if gamma.size else float("nan"),
        "allmaras_gamma_median": gamma_report,
        "allmaras_gamma_max": float(np.nanmax(gamma)) if gamma.size else float("nan"),
        "allmaras_rho_seed_min": float(np.nanmin(allmaras0.rho_kwt)),
        "allmaras_rho_seed_median": float(np.nanmedian(allmaras0.rho_kwt)),
        "allmaras_rho_seed_max": float(np.nanmax(allmaras0.rho_kwt)),
        "allmaras_delta_mod_seed_min_over_delta0": float(np.nanmin(allmaras0.delta_mod_over_delta0)),
        "allmaras_delta_mod_seed_median_over_delta0": float(np.nanmedian(allmaras0.delta_mod_over_delta0)),
        "allmaras_delta_mod_seed_max_over_delta0": float(np.nanmax(allmaras0.delta_mod_over_delta0)),
        "proximity_seed": proximity_seed_diag.as_dict(),
        "stationarity": stationarity_diag.as_dict(),
        "contact_recovery": recovery_diag.as_dict(),
        "continuity": continuity_diag.as_dict(),
        "terminal_psi": None if terminal_psi is None else str(terminal_psi),
        "normal_terminal_enforced": bool(terminal_psi is not None and np.any(terminal_site_mask)),
        "normal_terminal_n_nodes": int(np.count_nonzero(terminal_site_mask)),
        "normal_terminal_delta_max_over_delta0": (
            float(np.max(np.abs(psi_final[terminal_site_mask]))) if np.any(terminal_site_mask) else float("nan")
        ),
        "normal_terminal_delta_rms_over_delta0": (
            float(np.sqrt(np.mean(np.abs(psi_final[terminal_site_mask]) ** 2))) if np.any(terminal_site_mask) else float("nan")
        ),
        "target_current_A": target_current_A,
        "terminal_voltage_V": terminal_voltage(mesh.nodes, phi_final_V, length_m=mesh.length_m),
        "current_residual": float(current_residual(currents, mesh, material, target_current_A)),
        "native_si_current_scale_A_m2": float(native_diag.current_scale_A_m2),
        "native_si_residual_no_boundary_rms_A_m3": float(native_diag.residual_no_boundary_rms_A_m3),
        "native_si_residual_plus_boundary_rms_A_m3": float(native_diag.residual_plus_boundary_rms_A_m3),
        "native_si_residual_minus_boundary_rms_A_m3": float(native_diag.residual_minus_boundary_rms_A_m3),
        "native_si_selected_boundary_sign": float(native_diag.selected_boundary_sign),
        "native_si_boundary_currents_from_total_A": native_diag.boundary_currents_from_total_A,
        "usadel_current_available": bool(history.get("usadel_current_available", np.array([False], dtype=bool))[0]),
        "usadel_current_backend": str(history.get("usadel_current_backend", np.array(["unavailable"], dtype=object))[0]),
        "usadel_current_reason": str(history.get("usadel_current_reason", np.array(["not computed"], dtype=object))[0]),
        "usadel_vs_gl_edge_relative_l2_final": float(history.get("usadel_vs_gl_edge_relative_l2", np.array([float("nan")]))[-1]),
        "usadel_vs_gl_edge_max_abs_diff_A_m2_final": float(history.get("usadel_vs_gl_edge_max_abs_diff_A_m2", np.array([float("nan")]))[-1]),
        "usadel_supercurrent_max_A_m2_final": float(history.get("usadel_supercurrent_max_A_m2", np.array([float("nan")]))[-1]),
        "gl_supercurrent_max_A_m2_final": float(history.get("gl_supercurrent_max_A_m2", np.array([float("nan")]))[-1]),
        "eta_R_final": float(history["eta_R"][-1]) if history["eta_R"].size else float("nan"),
        "min_delta_over_delta0": float(np.min(np.abs(psi_final_J)) / material.delta0_J),
        "mean_delta_over_delta0": float(np.mean(np.abs(psi_final_J)) / material.delta0_J),
        "max_pairbreaking_ratio": float(np.nanmax(currents.node_pairbreaking_ratio)),
        "normal_current_max_A_m2": float(jn_max_A_m2),
        "total_current_max_A_m2": float(jt_max_A_m2),
        "normal_current_fraction_max": float(jn_max_A_m2 / max(jt_max_A_m2, 1.0e-300)),
        "delta0_meV": float(material.delta0_J / MEV_J),
        "boundary_currents_A": {
            "left_A": float(terminal_currents_A.get("left", 0.0)),
            "right_A": float(terminal_currents_A.get("right", 0.0)),
            "net_A": float(sum(terminal_currents_A.values())),
        },
        "terminal_neumann_current_unit_A": float(device.terminal_neumann_current_unit_A),
        "native_poisson_residual_rel_final": float(history.get("pytdgl_like_poisson_residual_rel", np.array([float("nan")]))[-1]),
        "native_poisson_residual_norm_final": float(history.get("pytdgl_like_poisson_residual_norm", np.array([float("nan")]))[-1]),
        "native_poisson_rhs_norm_final": float(history.get("pytdgl_like_poisson_rhs_norm", np.array([float("nan")]))[-1]),
        "native_boundary_rhs_norm_final": float(history.get("pytdgl_like_boundary_rhs_norm", np.array([float("nan")]))[-1]),
        "native_mu_boundary_max_abs_final": float(history.get("pytdgl_like_mu_boundary_max_abs", np.array([float("nan")]))[-1]),
        "allmaras_bulk_guard_layers": int(history.get("allmaras_bulk_guard_layers", np.array([1]))[0]),
        "allmaras_bulk_n_nodes": int(np.count_nonzero(history.get("allmaras_bulk_node_mask", np.zeros(mesh.n_nodes, dtype=bool)))),
        "allmaras_mismatch_div_rms_A_m3_final": float(history.get("allmaras_mismatch_div_rms_A_m3", np.array([float("nan")]))[-1]),
        "allmaras_mismatch_div_bulk_rms_A_m3_final": float(history.get("allmaras_mismatch_div_bulk_rms_A_m3", np.array([float("nan")]))[-1]),
        "allmaras_mismatch_div_max_abs_A_m3_final": float(history.get("allmaras_mismatch_div_max_abs_A_m3", np.array([float("nan")]))[-1]),
        "allmaras_mismatch_div_bulk_max_abs_A_m3_final": float(history.get("allmaras_mismatch_div_bulk_max_abs_A_m3", np.array([float("nan")]))[-1]),
        "allmaras_phase_drive_rms_over_delta0_final": float(history.get("allmaras_phase_drive_rms_over_delta0", np.array([float("nan")]))[-1]),
        "allmaras_phase_drive_bulk_rms_over_delta0_final": float(history.get("allmaras_phase_drive_bulk_rms_over_delta0", np.array([float("nan")]))[-1]),
        "allmaras_phase_drive_max_over_delta0_final": float(history.get("allmaras_phase_drive_max_over_delta0", np.array([float("nan")]))[-1]),
        "allmaras_phase_drive_bulk_max_over_delta0_final": float(history.get("allmaras_phase_drive_bulk_max_over_delta0", np.array([float("nan")]))[-1]),
        "allmaras_update_forcing_max_abs_final": float(history.get("allmaras_update_forcing_max_abs", np.array([float("nan")]))[-1]),
    }

    state = GTDGLStationaryState(
        psi_J=psi_final_J,
        phi_V=phi_final_V,
        Te_K=Te,
        Tph_K=Tph,
        currents=currents,
        metadata={
            "backend": "pytdgl_like_minimal_no_screening",
            "supercurrent_law": supercurrent_law,
            "length_scale_m": float(device.length_scale_m),
            "voltage_scale_V": float(device.voltage_scale_V),
            "terminal_neumann_current_unit_A": float(device.terminal_neumann_current_unit_A),
            "current_scale_A": float(device.current_scale_A),
            "native_si_current_scale_A_m2": float(native_diag.current_scale_A_m2),
            "allmaras_coefficients_backend": "appendix_b_allmaras_wz_update_v1",
            "allmaras_update_backend": "appendix_b_explicit_forcing_rho_kwt_wz_v1_contact_guarded",
            "allmaras_contact_correction_guard_layers": int(allmaras_contact_guard_layers),
            "pytdgl_reference": "loganbvh/py-tdgl solver/operator structure, MIT license",
            "first_magic_ready": magic_ready,
            "stationarity": stationarity_diag.as_dict(),
            "contact_recovery": recovery_diag.as_dict(),
            "continuity": continuity_diag.as_dict(),
        },
    )
    return RelaxationResult(state=state, history=history, summary=summary)



def _terminal_site_mask_from_device(device, n_nodes: int) -> np.ndarray:
    """Return a boolean mask for metallic normal-terminal sites."""
    mask = np.zeros(int(n_nodes), dtype=bool)
    try:
        terminal_info = device.terminal_info()
    except Exception:
        terminal_info = []
    for terminal in terminal_info:
        idx = np.asarray(getattr(terminal, "site_indices", []), dtype=np.int64)
        idx = idx[(idx >= 0) & (idx < mask.size)]
        if idx.size:
            mask[idx] = True
    return mask


def _terminal_edge_mask_from_device(device, ops: FVOperators) -> np.ndarray:
    """Return edges incident on normal-terminal sites.

    The GL current is automatically zero on these edges when terminal psi is
    clamped to zero.  Usadel-Poisson uses an external constitutive table, so we
    explicitly block the same contact edges to keep the metallic terminal
    condition consistent.
    """
    node_mask = _terminal_site_mask_from_device(device, ops.n_nodes)
    return node_mask[np.asarray(ops.edge_i, dtype=np.int64)] | node_mask[np.asarray(ops.edge_j, dtype=np.int64)]



def _expand_node_mask_by_edges(mask: np.ndarray, *, ops: FVOperators, layers: int) -> np.ndarray:
    """Expand a node mask by graph-neighbor layers on the FV edge graph."""
    out = np.asarray(mask, dtype=bool).reshape(-1).copy()
    if out.size != ops.n_nodes:
        raise ValueError(f"mask has length {out.size}, expected {ops.n_nodes}.")
    edge_i = np.asarray(ops.edge_i, dtype=np.int64)
    edge_j = np.asarray(ops.edge_j, dtype=np.int64)
    for _ in range(max(0, int(layers))):
        touch = out[edge_i] | out[edge_j]
        if not np.any(touch):
            break
        out[edge_i[touch]] = True
        out[edge_j[touch]] = True
    return out

def _normalize_supercurrent_law(value: str) -> str:
    law = str(value).strip().lower().replace("-", "_")
    aliases = {
        "gl": "gl",
        "pytdgl": "gl",
        "native_gl": "gl",
        "usadel": "usadel_poisson",
        "usadel_poisson": "usadel_poisson",
        "poisson_usadel": "usadel_poisson",
    }
    if law not in aliases:
        raise ValueError(
            "supercurrent_law must be one of gl or usadel_poisson "
            f"(got {value!r})."
        )
    return aliases[law]


def _build_usadel_poisson_supercurrent_override(
    *,
    usadel_catalog: Any | None,
    device,
    material: GTDGLMaterial,
    Te_K: np.ndarray,
    ops: FVOperators,
):
    scale = native_current_scale_A_m2(device)
    blocked_edge_mask = _terminal_edge_mask_from_device(device, ops)

    def usadel_poisson_supercurrent(psi_dimensionless: np.ndarray, gl_supercurrent_native: np.ndarray) -> np.ndarray:
        del gl_supercurrent_native
        diag = compute_usadel_supercurrent_diagnostic(
            usadel_catalog=usadel_catalog,
            psi_dimensionless=psi_dimensionless,
            material=material,
            Te_K=Te_K,
            ops=ops,
            blocked_edge_mask=blocked_edge_mask,
        )
        if not diag.available:
            raise RuntimeError(
                "--ss-supercurrent-law usadel-poisson requires a PRE Usadel "
                f"supercurrent table. Diagnostic reason: {diag.reason}"
            )
        return np.asarray(diag.edge_js_usadel_A_m2, dtype=float) / max(scale, 1.0e-300)

    return usadel_poisson_supercurrent


def _build_allmaras_forcing_callback(
    *,
    usadel_catalog: Any | None,
    device,
    material: GTDGLMaterial,
    Te_K: np.ndarray,
    ops: FVOperators,
    blocked_edge_mask: np.ndarray,
    require_usadel: bool,
    contact_guard_layers: int = 2,
):
    """Build the Appendix-B forcing callback used by ``TDGLSolver``.

    The callback is explicit in the current order parameter.  For the official
    ``usadel_poisson`` path it uses the PRE Matsubara/Usadel supercurrent table
    for both Poisson and the Allmaras current-divergence correction.
    """

    Te = np.asarray(Te_K, dtype=float).copy()
    L0 = float(device.length_scale_m)
    terminal_node_mask = _terminal_site_mask_from_device(device, ops.n_nodes)
    contact_guard_node_mask = _expand_node_mask_by_edges(
        terminal_node_mask,
        ops=ops,
        layers=max(0, int(contact_guard_layers)),
    )

    def callback(psi_dimensionless: np.ndarray, psi_laplacian) -> np.ndarray:
        psi = np.asarray(psi_dimensionless, dtype=np.complex128)
        edge_js = None
        if require_usadel:
            diag = compute_usadel_supercurrent_diagnostic(
                usadel_catalog=usadel_catalog,
                psi_dimensionless=psi,
                material=material,
                Te_K=Te,
                ops=ops,
                blocked_edge_mask=blocked_edge_mask,
            )
            if not diag.available:
                raise RuntimeError(
                    "Appendix-B Allmaras update with usadel_poisson requires a PRE "
                    f"Matsubara supercurrent table. Diagnostic reason: {diag.reason}"
                )
            edge_js = diag.edge_js_usadel_A_m2

        forcing = compute_allmaras_forcing_dimensionless(
            psi_dimensionless=psi,
            psi_laplacian_dimensionless=psi_laplacian @ psi,
            material=material,
            Te_K=Te,
            ops=ops,
            length_scale_m=L0,
            edge_js_usadel_A_m2=edge_js,
            blocked_edge_mask=blocked_edge_mask,
        )
        out = np.asarray(forcing.forcing_dimensionless, dtype=np.complex128).copy()
        # Do not assume the Allmaras current-divergence correction is valid in
        # the normal-metal contact conversion region.  Keep diffusion/reaction
        # active but remove only the imaginary mismatch drive on the guarded
        # contact nodes.
        if np.any(contact_guard_node_mask):
            out[contact_guard_node_mask] = (
                np.asarray(forcing.diffusion_dimensionless, dtype=np.complex128)[contact_guard_node_mask]
                + np.asarray(forcing.reaction_dimensionless, dtype=np.complex128)[contact_guard_node_mask]
            )
        return out

    return callback


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
        "normal_terminal_node_mask": terminal_site_mask.astype(bool),
        "normal_terminal_edge_mask": blocked_edge_mask.astype(bool),
    }

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

    if psi_snap.ndim != 2 or psi_snap.shape[1] != mesh.n_nodes or raw_snap_t.size == 0:
        ns = max(2, int(n_snapshots))
        snap_t = np.linspace(0.0, float(t_s[-1]) if t_s.size else 0.0, ns)
        psi_snap_J = np.tile(psi_final_J, (ns, 1))
        phi_snap_V = np.tile(phi_final_V, (ns, 1))
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
            Te_K=Te,
            target_current_A=target_current_A,
        )
        udiag = compute_usadel_supercurrent_diagnostic(
            usadel_catalog=usadel_catalog,
            psi_dimensionless=psi_dim,
            material=material,
            Te_K=Te,
            ops=ops,
            blocked_edge_mask=blocked_edge_mask,
        )
        adiag = compute_allmaras_appendix_b_diagnostic(
            psi_dimensionless=psi_dim,
            material=material,
            Te_K=Te,
            ops=ops,
            terminal_node_mask=terminal_site_mask,
            blocked_edge_mask=blocked_edge_mask,
            edge_js_usadel_A_m2=(udiag.edge_js_usadel_A_m2 if udiag.available else None),
            bulk_guard_layers=1,
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
    allmaras_C = np.vstack([d.coefficients.correction_C_J_m3_A for d in allmaras_diags])

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
        "allmaras_correction_C_snapshot_J_m3_A": allmaras_C,
        "allmaras_bulk_node_mask": allmaras_bulk_mask.astype(bool),
        "allmaras_bulk_guard_layers": np.array([1], dtype=np.int64),
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
