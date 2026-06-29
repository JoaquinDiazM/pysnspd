"""Adapters between pySNSPD OE7 data and the pyTDGL-like solver core."""
from __future__ import annotations

from types import SimpleNamespace
from typing import Any

import numpy as np

from pysnspd.gtdgl.material import GTDGLMaterial
from pysnspd.gtdgl.operators import FVOperators, terminal_voltage
from pysnspd.gtdgl.state import GTDGLStationaryState, RelaxationResult
from pysnspd.gtdgl.fields import compute_current_fields
from pysnspd.gtdgl.diagnostics import (
    current_residual,
    current_density_maxima_A_m2,
    seed_target_current_A,
    target_current_density_A_m2,
)
from .device import build_pytdgl_like_device
from .options import SolverOptions, SparseSolver
from .solver import TDGLSolver

MEV_J = 1.602176634e-22


def solve_stationary_pytdgl_like(
    *,
    mesh,
    edge_data,
    seed,
    material: GTDGLMaterial,
    ops: FVOperators,
    steps: int = 2000,
    dt_s: float = 1.0e-17,
    target_current_A: float | None = None,
    terminal_psi: complex | float | None = 0.0,
    adaptive: bool = True,
    adaptive_window: int = 10,
    max_solve_retries: int = 10,
    adaptive_time_step_multiplier: float = 0.25,
    n_snapshots: int = 6,
) -> RelaxationResult:
    """Run the essential pyTDGL-like stationary solver on a pySNSPD seed.

    The state evolved internally is dimensionless, as in pyTDGL.  The returned
    ``RelaxationResult`` is converted back to pySNSPD physical units so existing
    diagnostics and plotting code can consume it.
    """

    if steps <= 0:
        raise ValueError("steps must be positive.")
    if dt_s <= 0:
        raise ValueError("dt_s must be positive.")
    if target_current_A is None:
        target_current_A = seed_target_current_A(seed)
    target_current_A = float(target_current_A)

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
    solve_time = max(1, int(steps)) * dt_dimless

    options = SolverOptions(
        solve_time=solve_time,
        dt_init=dt_dimless,
        dt_max=dt_dimless if not adaptive else max(dt_dimless, 100.0 * dt_dimless),
        adaptive=bool(adaptive),
        adaptive_window=int(adaptive_window),
        max_solve_retries=int(max_solve_retries),
        adaptive_time_step_multiplier=float(adaptive_time_step_multiplier),
        terminal_psi=terminal_psi,
        sparse_solver=SparseSolver.SUPERLU,
        include_screening=False,
    )

    # pySNSPD modified coefficients inside the pyTDGL algebra.
    delta_mod2 = material.delta_mod_squared_J2(Te)
    epsilon = np.clip(delta_mod2 / (material.delta0_J**2), 0.0, 1.0)
    tau_sc = material.tau_sc_s(Te)
    gamma = float(np.nanmedian(2.0 * material.delta0_J * tau_sc / 1.054571817e-34))
    if not np.isfinite(gamma) or gamma < 0:
        gamma = 0.0
    # Keep u scalar to preserve pyTDGL's method signature.
    rho0 = material.rho_kwt(Te, np.maximum(np.abs(psi0_J), 1.0e-12 * material.delta0_J))
    u = float(np.nanmedian(np.maximum(rho0, 1.0e-12)))
    if not np.isfinite(u) or u <= 0:
        u = 1.0
    device.layer.u = u
    device.layer.gamma = gamma

    I_norm = 0.0 if target_current_A == 0 else 1.0
    terminal_currents = {"left": -I_norm, "right": I_norm}

    psi0 = psi0_J / material.delta0_J
    mu0 = (phi0_V - float(np.mean(phi0_V))) / device.voltage_scale_V
    seed_solution = {
        "psi": psi0,
        "mu": mu0,
        "supercurrent": np.zeros(ops.n_edges, dtype=float),
        "normal_current": np.zeros(ops.n_edges, dtype=float),
    }

    solver = TDGLSolver(
        device=device,
        options=options,
        applied_vector_potential=0.0,
        terminal_currents=terminal_currents,
        disorder_epsilon=lambda r: epsilon,
        seed_solution=seed_solution,
    )
    solution = solver.solve()
    if solution is None:
        raise RuntimeError("pytdgl_like solver returned None.")

    psi_final = solution.tdgl_data.psi
    mu_final = solution.tdgl_data.mu
    psi_final_J = psi_final * material.delta0_J
    phi_final_V = mu_final * device.voltage_scale_V
    phi_final_V = phi_final_V - float(np.mean(phi_final_V))

    currents = compute_current_fields(
        psi_J=psi_final_J,
        phi_V=phi_final_V,
        Te_K=Te,
        material=material,
        ops=ops,
    )

    history = _build_history(
        solution=solution,
        mesh=mesh,
        ops=ops,
        material=material,
        Te=Te,
        Tph=Tph,
        psi_final_J=psi_final_J,
        phi_final_V=phi_final_V,
        currents=currents,
        target_current_A=target_current_A,
        n_snapshots=n_snapshots,
    )

    jmax = current_density_maxima_A_m2(currents)
    summary: dict[str, Any] = {
        "backend": "pytdgl_like_minimal_no_screening",
        "converged": False,
        "accepted_steps": int(solution.history.get("final_step", np.array([0]))[0]),
        "rejected_steps": 0,
        "final_time_ps": float(solution.history.get("final_time", np.array([0.0]))[0] * tau0 / 1.0e-12),
        "dt_init_s": float(dt_s),
        "tau0_GL_s": tau0,
        "pytdgl_u": u,
        "pytdgl_gamma": gamma,
        "terminal_psi": None if terminal_psi is None else str(terminal_psi),
        "target_current_A": target_current_A,
        "terminal_voltage_V": terminal_voltage(mesh.nodes, phi_final_V, length_m=mesh.length_m),
        "current_residual": float(current_residual(currents, mesh)),
        "eta_R_final": float(history["eta_R"][-1]) if history["eta_R"].size else float("nan"),
        "min_delta_over_delta0": float(np.min(np.abs(psi_final_J)) / material.delta0_J),
        "mean_delta_over_delta0": float(np.mean(np.abs(psi_final_J)) / material.delta0_J),
        "max_pairbreaking_ratio": float(np.nanmax(currents.node_pairbreaking_ratio)),
        "normal_current_max_A_m2": jmax["normal_current_max_A_m2"],
        "total_current_max_A_m2": jmax["total_current_max_A_m2"],
        "normal_current_fraction_max": float(jmax["normal_current_max_A_m2"] / max(jmax["total_current_max_A_m2"], 1.0e-300)),
        "delta0_meV": float(material.delta0_J / MEV_J),
        "boundary_currents_A": {},
    }

    state = GTDGLStationaryState(
        psi_J=psi_final_J,
        phi_V=phi_final_V,
        Te_K=Te,
        Tph_K=Tph,
        currents=currents,
        metadata={
            "backend": "pytdgl_like_minimal_no_screening",
            "length_scale_m": float(device.length_scale_m),
            "voltage_scale_V": float(device.voltage_scale_V),
            "current_scale_A": float(device.current_scale_A),
            "pytdgl_reference": "loganbvh/py-tdgl solver/operator structure, MIT license",
        },
    )
    return RelaxationResult(state=state, history=history, summary=summary)


def _build_history(
    *,
    solution,
    mesh,
    ops: FVOperators,
    material: GTDGLMaterial,
    Te: np.ndarray,
    Tph: np.ndarray,
    psi_final_J: np.ndarray,
    phi_final_V: np.ndarray,
    currents,
    target_current_A: float,
    n_snapshots: int,
) -> dict[str, np.ndarray]:
    raw = solution.history
    tau0 = float(material.tau0_GL_s)
    dt_dimless = np.asarray(raw.get("dt", []), dtype=float)
    t_s = np.cumsum(dt_dimless) * tau0 if dt_dimless.size else np.array([0.0])
    max_d = np.asarray(raw.get("max_d_abs_sq_psi", np.zeros_like(t_s)), dtype=float)
    if max_d.size != t_s.size:
        max_d = np.resize(max_d, t_s.size)
    mu_ptp = np.asarray(raw.get("mu_ptp", np.zeros_like(t_s)), dtype=float) * solution.device.voltage_scale_V
    if mu_ptp.size != t_s.size:
        mu_ptp = np.resize(mu_ptp, t_s.size)

    delta_abs = np.abs(psi_final_J)
    javg = abs(target_current_density_A_m2(material, target_current_A))
    jmax = current_density_maxima_A_m2(currents)
    residual = float(current_residual(currents, mesh))

    hist: dict[str, np.ndarray] = {
        "t_s": t_s,
        "dt_s": dt_dimless * tau0,
        "eta_R": max_d,
        "max_amp2_change_rel": max_d,
        "current_residual": np.full(t_s.shape, residual),
        "pairbreaking_max": np.full(t_s.shape, float(np.nanmax(currents.node_pairbreaking_ratio))),
        "terminal_voltage_V": mu_ptp,
        "delta_min_over_delta0": np.full(t_s.shape, float(np.min(delta_abs) / material.delta0_J)),
        "delta_max_over_delta0": np.full(t_s.shape, float(np.max(delta_abs) / material.delta0_J)),
        "normal_current_max_A_m2": np.full(t_s.shape, jmax["normal_current_max_A_m2"]),
        "total_current_max_A_m2": np.full(t_s.shape, jmax["total_current_max_A_m2"]),
        "delta0_meV": np.array([float(material.delta0_J / MEV_J)]),
        "javg_A_m2": np.array([javg]),
    }

    # Store final-state snapshots repeatedly at representative times.  This is
    # enough for existing plotting utilities and avoids pretending we have a full
    # HDF5 pyTDGL trajectory in this first comparison backend.
    ns = max(2, int(n_snapshots))
    snap_t = np.linspace(0.0, float(t_s[-1]) if t_s.size else 0.0, ns)
    delta_mev = np.abs(psi_final_J) / MEV_J
    jtot_mag = np.sqrt(currents.node_jtot_x_A_m2**2 + currents.node_jtot_y_A_m2**2)
    js_mag = np.sqrt(currents.node_js_us_x_A_m2**2 + currents.node_js_us_y_A_m2**2)
    jn_mag = np.sqrt(currents.node_jn_x_A_m2**2 + currents.node_jn_y_A_m2**2)
    for key, arr in {
        "snapshot_t_s": snap_t,
        "psi_snapshot_real_J": np.tile(np.real(psi_final_J), (ns, 1)),
        "psi_snapshot_imag_J": np.tile(np.imag(psi_final_J), (ns, 1)),
        "delta_snapshot_meV": np.tile(delta_mev, (ns, 1)),
        "phi_snapshot_V": np.tile(phi_final_V, (ns, 1)),
        "current_density_snapshot_A_m2": np.tile(jtot_mag, (ns, 1)),
        "current_density_snapshot_x_A_m2": np.tile(currents.node_jtot_x_A_m2, (ns, 1)),
        "current_density_snapshot_y_A_m2": np.tile(currents.node_jtot_y_A_m2, (ns, 1)),
        "supercurrent_density_snapshot_A_m2": np.tile(js_mag, (ns, 1)),
        "supercurrent_density_snapshot_x_A_m2": np.tile(currents.node_js_us_x_A_m2, (ns, 1)),
        "supercurrent_density_snapshot_y_A_m2": np.tile(currents.node_js_us_y_A_m2, (ns, 1)),
        "normal_current_density_snapshot_A_m2": np.tile(jn_mag, (ns, 1)),
        "normal_current_density_snapshot_x_A_m2": np.tile(currents.node_jn_x_A_m2, (ns, 1)),
        "normal_current_density_snapshot_y_A_m2": np.tile(currents.node_jn_y_A_m2, (ns, 1)),
        "divergence_snapshot_A_m3": np.tile(currents.node_div_jtot_A_m3, (ns, 1)),
        "pairbreaking_ratio_snapshot": np.tile(currents.node_pairbreaking_ratio, (ns, 1)),
        "edge_Q_snapshot_m_inv": np.tile(currents.edge_Q_m_inv, (ns, 1)),
        "edge_js_us_snapshot_A_m2": np.tile(currents.edge_js_us_A_m2, (ns, 1)),
        "edge_jn_snapshot_A_m2": np.tile(currents.edge_jn_A_m2, (ns, 1)),
        "edge_jtot_snapshot_A_m2": np.tile(currents.edge_jtot_A_m2, (ns, 1)),
    }.items():
        hist[key] = np.asarray(arr)
    return hist
