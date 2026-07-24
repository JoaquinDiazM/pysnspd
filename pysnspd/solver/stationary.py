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
from .steady_gate import build_total_stationarity_callback
from pysnspd.thermal.evolution import ThermalRuntimeConfig, ThermalRuntimeController, thermal_stationarity_diagnostics
from pysnspd.circuit.readout import (
    CircuitParams,
    CircuitRuntimeConfig,
    CircuitRuntimeController,
    central_tdgl_voltage_V,
    circuit_stationarity_diagnostics,
)
from .targets import (
    apply_terminal_proximity_seed,
    contact_recovery_diagnostics,
    continuity_diagnostics,
    dynamic_stationarity_diagnostics,
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
    dynamic_stationarity_tail_snapshots: int = 4,
    dynamic_stationarity_minimum_tail_ps: float = 2.0,
    dynamic_stationarity_profile_rel: float = 5.0e-2,
    dynamic_stationarity_voltage_rel: float = 5.0e-2,
    dynamic_stationarity_psl_threshold: float = 0.75,
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
    allmaras_phase_direct_amplitude_fraction: float = 1.0e-2,
    allmaras_phase_convergence_tol: float = 1.0e-3,
    allmaras_phase_convergence_max_iterations: int = 64,
    thermal_enabled: bool = False,
    thermal_power_table_npz: str | None = None,
    thermal_window_m: float = 100.0e-9,
    thermal_start_time_s: float = 2.0e-12,
    thermal_bath_K: float | None = None,
    thermal_min_K: float | None = None,
    thermal_max_K: float | None = None,
    thermal_max_step_K: float = 0.05,
    thermal_max_substeps: int = 64,
    thermal_stationarity_rate_K_per_ps: float = 1.0e-2,
    circuit_enabled: bool = False,
    circuit_params: CircuitParams | None = None,
    circuit_runtime_config: CircuitRuntimeConfig | None = None,
    center_voltage_width_m: float = 100.0e-9,
    center_voltage_probe_band_m: float | None = None,
    stop_on_total_stationarity: bool = False,
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

    psi0_J = (
        np.asarray(seed.node_psi_real_J, dtype=float)
        + 1j * np.asarray(seed.node_psi_imag_J, dtype=float)
    )
    phi0_V = np.asarray(seed.node_phi_electric_V, dtype=float)
    Te = np.asarray(seed.node_Te_K, dtype=float).copy()
    Tph = np.asarray(seed.node_Tph_K, dtype=float).copy()
    bath_K = float(thermal_bath_K) if thermal_bath_K is not None else float(np.nanmedian(Tph))
    thermal_controller = None
    thermal_config = ThermalRuntimeConfig(
        enabled=bool(thermal_enabled and thermal_power_table_npz),
        window_m=float(thermal_window_m),
        start_time_s=float(thermal_start_time_s),
        bath_K=bath_K,
        min_K=thermal_min_K,
        max_K=thermal_max_K,
        max_step_K=float(thermal_max_step_K),
        max_substeps=int(thermal_max_substeps),
        stationarity_rate_K_per_ps=float(thermal_stationarity_rate_K_per_ps),
    )
    if thermal_config.enabled:
        thermal_controller = ThermalRuntimeController(
            nodes_m=np.asarray(mesh.nodes, dtype=float)[:, :2],
            ops=ops,
            material=material,
            Te_K=Te,
            Tph_K=Tph,
            power_table_npz=str(thermal_power_table_npz),
            config=thermal_config,
        )

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

    circuit_controller = None
    circuit_runtime = (
        circuit_runtime_config or CircuitRuntimeConfig(start_time_s=thermal_start_time_s)
    ).validated()
    if circuit_enabled:
        initial_center_voltage = central_tdgl_voltage_V(
            nodes_m=np.asarray(mesh.nodes, dtype=float)[:, :2],
            phi_V=phi0_V,
            center_width_m=float(center_voltage_width_m),
            probe_band_m=center_voltage_probe_band_m,
        )
        circuit_controller = CircuitRuntimeController(
            I_ss_A=target_current_A,
            V_tdgl_ss_V=initial_center_voltage,
            params=circuit_params or CircuitParams(),
            config=circuit_runtime,
        )

    # pySNSPD policy: terminal currents are always SI amperes.  The
    # pyTDGL-like solver converts these physical currents to the internal
    # Neumann value needed by its Poisson operator using sigma_n, film
    # thickness and the Josephson voltage scale.
    if circuit_controller is None:
        terminal_currents_A = {"left": -target_current_A, "right": target_current_A}
    else:
        def terminal_currents_A(time_dimensionless: float) -> dict[str, float]:
            current_A = circuit_controller.terminal_current_A(
                float(time_dimensionless) * tau0
            )
            return {"left": -current_A, "right": current_A}

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

        The field is rebuilt from the mutable runtime ``Te`` array.  This is
        what makes the thermal coupling feed back into the local Allmaras/KWT
        coefficients on the next gTDGL step.
        """
        del t
        arr = np.asarray(r)
        coeff = allmaras_coefficients(
            psi_dimensionless=np.ones_like(Te, dtype=np.complex128),
            material=material,
            Te_K=Te,
        )
        eps_now = np.clip(np.asarray(coeff.solver_epsilon, dtype=float), 0.0, 1.0)
        gamma_now = np.asarray(coeff.gamma_kwt_dimensionless, dtype=float)
        gamma_now = np.where(np.isfinite(gamma_now) & (gamma_now >= 0.0), gamma_now, 0.0)
        device.layer.gamma = gamma_now
        if vectorized and arr.ndim >= 2:
            if arr.shape[0] != eps_now.size:
                return np.full(arr.shape[0], float(np.nanmedian(eps_now)), dtype=float)
            return eps_now.copy()
        return float(np.nanmedian(eps_now))

    supercurrent_override = None
    if supercurrent_law == "usadel_poisson":
        supercurrent_override = _build_usadel_poisson_supercurrent_override(
            usadel_catalog=usadel_catalog,
            device=device,
            material=material,
            Te_K=Te,
            ops=ops,
        )

    phase_drive_continuation = PhaseDriveContinuationSolver.from_operators(
        ops,
        direct_amplitude_fraction=float(allmaras_phase_direct_amplitude_fraction),
        tolerance=float(allmaras_phase_convergence_tol),
        max_iterations=int(allmaras_phase_convergence_max_iterations),
    )
    allmaras_forcing_callback = _build_allmaras_forcing_callback(
        usadel_catalog=usadel_catalog,
        device=device,
        material=material,
        Te_K=Te,
        ops=ops,
        blocked_edge_mask=_terminal_edge_mask_from_device(device, ops),
        require_usadel=(supercurrent_law == "usadel_poisson"),
        phase_drive_continuation=phase_drive_continuation,
    )

    thermal_step_callback = None
    thermal_snapshot_callback = None
    if thermal_controller is not None:
        current_scale_A_m2 = native_current_scale_A_m2(device)

        def thermal_step_callback(**kwargs):
            return thermal_controller.step(
                time_s=float(kwargs["time"]) * tau0,
                dt_s=float(kwargs["dt"]) * tau0,
                psi_dimensionless=np.asarray(kwargs["psi"], dtype=np.complex128),
                native_normal_current=np.asarray(kwargs["normal_current"], dtype=float),
                current_scale_A_m2=float(current_scale_A_m2),
            )

        thermal_snapshot_callback = thermal_controller.snapshot_payload

    circuit_step_callback = None
    circuit_snapshot_callback = None
    if circuit_controller is not None:
        def circuit_step_callback(**kwargs):
            phi_V = np.asarray(kwargs["mu"], dtype=float) * device.voltage_scale_V
            phi_V = phi_V - float(np.mean(phi_V))
            V_tdgl = central_tdgl_voltage_V(
                nodes_m=np.asarray(mesh.nodes, dtype=float)[:, :2],
                phi_V=phi_V,
                center_width_m=float(center_voltage_width_m),
                probe_band_m=center_voltage_probe_band_m,
            )
            return circuit_controller.step(
                time_s=float(kwargs["time"]) * tau0,
                dt_s=float(kwargs["dt"]) * tau0,
                V_tdgl_V=V_tdgl,
            )

        circuit_snapshot_callback = circuit_controller.snapshot_payload

    early_stop_state: dict[str, Any] = {}
    early_stop_callback = None
    if stop_on_total_stationarity:
        early_stop_callback = build_total_stationarity_callback(
            mesh=mesh,
            edge_data=edge_data,
            ops=ops,
            material=material,
            device=device,
            target_current_A=target_current_A,
            circuit_controller=circuit_controller,
            circuit_runtime=circuit_runtime,
            thermal_config=thermal_config,
            Te_K=Te,
            Tph_K=Tph,
            terminal_healing_fraction=terminal_healing_fraction,
            recovery_min_xi=recovery_min_xi,
            recovery_max_xi=recovery_max_xi,
            stationarity_phase_gradient_rel=phase_gradient_rel_tol,
            stationarity_phi_gradient_rel=phi_gradient_rel_tol,
            stationarity_q_abs_m_inv=stationarity_q_abs_m_inv,
            stationarity_phi_gradient_abs_V_m=stationarity_phi_gradient_abs_V_m,
            stationarity_edge_active_threshold=stationarity_edge_active_threshold,
            stationarity_bulk_exclusion_xi=stationarity_bulk_exclusion_xi,
            dynamic_stationarity_tail_snapshots=dynamic_stationarity_tail_snapshots,
            dynamic_stationarity_minimum_tail_ps=dynamic_stationarity_minimum_tail_ps,
            dynamic_stationarity_profile_rel=dynamic_stationarity_profile_rel,
            dynamic_stationarity_voltage_rel=dynamic_stationarity_voltage_rel,
            dynamic_stationarity_psl_threshold=dynamic_stationarity_psl_threshold,
            continuity_rms_tol=continuity_rms_tol,
            continuity_max_tol=continuity_max_tol,
            continuity_poisson_tol=continuity_poisson_tol,
            thermal_stationarity_rate_K_per_ps=thermal_stationarity_rate_K_per_ps,
            requested_total_time_s=float(total_time_s),
            tau0=tau0,
            state=early_stop_state,
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
        thermal_step_callback=thermal_step_callback,
        thermal_snapshot_callback=thermal_snapshot_callback,
        circuit_step_callback=circuit_step_callback,
        circuit_snapshot_callback=circuit_snapshot_callback,
        early_stop_callback=early_stop_callback,
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

    final_target_current_A = (
        float(circuit_controller.state.I_s_A)
        if circuit_controller is not None
        else target_current_A
    )
    final_terminal_currents_A = {
        "left": -final_target_current_A,
        "right": final_target_current_A,
    }
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
        target_current_A=final_target_current_A,
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
        target_current_A=final_target_current_A,
        usadel_catalog=usadel_catalog,
        n_snapshots=n_snapshots,
        stationarity_bulk_exclusion_xi=float(stationarity_bulk_exclusion_xi),
        phase_drive_continuation=phase_drive_continuation,
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
    dynamic_stationarity_diag = dynamic_stationarity_diagnostics(
        history=history,
        nodes_m=np.asarray(mesh.nodes, dtype=float)[:, :2],
        delta0_J=float(material.delta0_J),
        tail_snapshots=int(dynamic_stationarity_tail_snapshots),
        minimum_tail_duration_ps=float(dynamic_stationarity_minimum_tail_ps),
        profile_relative_tolerance=float(dynamic_stationarity_profile_rel),
        voltage_relative_tolerance=float(dynamic_stationarity_voltage_rel),
        psl_threshold_over_delta0=float(dynamic_stationarity_psl_threshold),
        bulk_exclusion_xi=float(stationarity_bulk_exclusion_xi),
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
        target_current_A=final_target_current_A,
        history=history,
        rms_tol=float(continuity_rms_tol),
        max_tol=float(continuity_max_tol),
        poisson_tol=float(continuity_poisson_tol),
    )
    thermal_stationarity_diag = thermal_stationarity_diagnostics(
        history,
        enabled=bool(thermal_config.enabled),
        start_time_s=float(thermal_config.start_time_s),
        requested_total_time_s=float(total_time_s),
        rate_tol_K_per_ps=float(thermal_config.stationarity_rate_K_per_ps),
    )
    circuit_stationarity_diag = circuit_stationarity_diagnostics(
        history,
        config=circuit_runtime,
    ) if circuit_controller is not None else {
        "enabled": False,
        "passes": True,
        "reason": "circuit coupling disabled",
    }
    magic_ready = bool(
        stationarity_diag.passes
        and dynamic_stationarity_diag.passes
        and recovery_diag.passes
        and continuity_diag.passes
        and bool(thermal_stationarity_diag.get("passes", False))
        and bool(circuit_stationarity_diag.get("passes", False))
    )

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
        "dynamic_stationarity_passes": bool(dynamic_stationarity_diag.passes),
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
        "material_D_base_m2_s": float(material.D_base_m2_s if material.D_base_m2_s is not None else material.D_m2_s),
        "material_D_effective_factor": float(material.D_effective_factor),
        "material_D_m2_s": float(material.D_m2_s),
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
        "allmaras_update_backend": "appendix_b_normalized_phase_drive_harmonic_continuation_v2",
        "allmaras_phase_continuation": {
            "method": "jacobi_preconditioned_cg_harmonic_continuation",
            "direct_amplitude_fraction": float(allmaras_phase_direct_amplitude_fraction),
            "tolerance": float(allmaras_phase_convergence_tol),
            "max_iterations": int(allmaras_phase_convergence_max_iterations),
            "final_converged": bool(history.get("allmaras_phase_convergence_converged", np.array([False]))[-1]),
            "final_iterations": int(history.get("allmaras_phase_convergence_iterations", np.array([0]))[-1]),
            "final_residual_rel": float(history.get("allmaras_phase_convergence_residual_rel", np.array([float("nan")]))[-1]),
        },
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
        "dynamic_stationarity": dynamic_stationarity_diag.as_dict(),
        "contact_recovery": recovery_diag.as_dict(),
        "continuity": continuity_diag.as_dict(),
        "thermal_stationarity": thermal_stationarity_diag,
        "circuit_stationarity": circuit_stationarity_diag,
        "circuit_runtime": {
            "enabled": circuit_controller is not None,
            "start_time_ps": float(circuit_runtime.start_time_s / 1.0e-12),
            "params": (
                circuit_controller.params.as_dict()
                if circuit_controller is not None
                else None
            ),
            "final_state": (
                circuit_controller.state.as_dict()
                if circuit_controller is not None
                else None
            ),
            "early_stop_diagnostics": dict(early_stop_state),
        },
        "thermal_runtime": {
            "enabled": bool(thermal_config.enabled),
            "power_table_npz": None if thermal_power_table_npz is None else str(thermal_power_table_npz),
            "window_nm": float(thermal_config.window_m / 1.0e-9),
            "start_time_ps": float(thermal_config.start_time_s / 1.0e-12),
            "bath_K": float(thermal_config.bath_K),
            "min_K": None if thermal_config.min_K is None else float(thermal_config.min_K),
            "max_K": None if thermal_config.max_K is None else float(thermal_config.max_K),
            "max_step_K": float(thermal_config.max_step_K),
            "max_substeps": int(thermal_config.max_substeps),
            "active_n_nodes": int(np.count_nonzero(thermal_controller.mask)) if thermal_controller is not None else 0,
        },
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
        "current_residual": float(current_residual(currents, mesh, material, final_target_current_A)),
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
            "left_A": float(final_terminal_currents_A["left"]),
            "right_A": float(final_terminal_currents_A["right"]),
            "net_A": float(sum(final_terminal_currents_A.values())),
        },
        "terminal_neumann_current_unit_A": float(device.terminal_neumann_current_unit_A),
        "native_poisson_residual_rel_final": float(history.get("pytdgl_like_poisson_residual_rel", np.array([float("nan")]))[-1]),
        "native_poisson_residual_norm_final": float(history.get("pytdgl_like_poisson_residual_norm", np.array([float("nan")]))[-1]),
        "native_poisson_rhs_norm_final": float(history.get("pytdgl_like_poisson_rhs_norm", np.array([float("nan")]))[-1]),
        "native_boundary_rhs_norm_final": float(history.get("pytdgl_like_boundary_rhs_norm", np.array([float("nan")]))[-1]),
        "native_mu_boundary_max_abs_final": float(history.get("pytdgl_like_mu_boundary_max_abs", np.array([float("nan")]))[-1]),
        "allmaras_bulk_mask_policy": str(history.get("allmaras_bulk_mask_policy", np.array(["terminal_nodes_excluded_only"]))[0]),
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
            "allmaras_update_backend": "appendix_b_normalized_phase_drive_harmonic_continuation_v2",
            "allmaras_phase_continuation_method": "jacobi_preconditioned_cg_harmonic_continuation",
            "pytdgl_reference": "loganbvh/py-tdgl solver/operator structure, MIT license",
            "first_magic_ready": magic_ready,
            "dynamic_stationarity_passes": bool(dynamic_stationarity_diag.passes),
            "stationarity": stationarity_diag.as_dict(),
            "dynamic_stationarity": dynamic_stationarity_diag.as_dict(),
            "contact_recovery": recovery_diag.as_dict(),
            "continuity": continuity_diag.as_dict(),
            "thermal_stationarity": thermal_stationarity_diag,
            "circuit_stationarity": circuit_stationarity_diag,
            "circuit_runtime": {
                "enabled": circuit_controller is not None,
                "params": (
                    circuit_controller.params.as_dict()
                    if circuit_controller is not None
                    else None
                ),
                "final_state": (
                    circuit_controller.state.as_dict()
                    if circuit_controller is not None
                    else None
                ),
            },
        },
    )
    return RelaxationResult(state=state, history=history, summary=summary)


from pysnspd.solver.callbacks import (
    _build_allmaras_forcing_callback,
    _build_usadel_poisson_supercurrent_override,
    _normalize_supercurrent_law,
    _terminal_edge_mask_from_device,
    _terminal_site_mask_from_device,
)
from pysnspd.solver.history import _build_history
