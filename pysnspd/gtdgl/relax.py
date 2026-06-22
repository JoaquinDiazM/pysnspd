"""Stationary gTDGL/Poisson relaxation for pySNSPD OE7.

The OE7 solver starts from the analytic OE6 seed, keeps Te and Tph frozen,
does not activate the external circuit, and advances only the mesoscopic
gTDGL/Poisson sector.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any
import math

import numpy as np

from pysnspd.gtdgl.material import (
    E_CHARGE_C,
    HBAR_J_S,
    K_B_J_K,
    GTDGLMaterial,
)
from pysnspd.gtdgl.operators import (
    FVOperators,
    boundary_currents_from_node_vectors,
    divergence_from_edge_scalar,
    edge_average,
    edge_phase_gradient_from_psi,
    edge_scalar_gradient,
    edge_scalar_to_node_vector,
    laplacian,
    terminal_voltage,
)

try:
    from scipy.sparse import coo_matrix
    from scipy.sparse.linalg import spsolve
except Exception:
    coo_matrix = None
    spsolve = None


@dataclass(frozen=True)
class CurrentFields:
    """Edge and node current diagnostics."""

    edge_Q_m_inv: np.ndarray
    edge_js_us_A_m2: np.ndarray
    edge_js_gl_A_m2: np.ndarray
    edge_jn_A_m2: np.ndarray
    edge_jtot_A_m2: np.ndarray
    node_div_js_us_A_m3: np.ndarray
    node_div_js_gl_A_m3: np.ndarray
    node_div_jtot_A_m3: np.ndarray
    node_js_us_x_A_m2: np.ndarray
    node_js_us_y_A_m2: np.ndarray
    node_jn_x_A_m2: np.ndarray
    node_jn_y_A_m2: np.ndarray
    node_jtot_x_A_m2: np.ndarray
    node_jtot_y_A_m2: np.ndarray


@dataclass(frozen=True)
class GTDGLStationaryState:
    """Node-based stationary gTDGL state."""

    psi_J: np.ndarray
    phi_V: np.ndarray
    Te_K: np.ndarray
    Tph_K: np.ndarray
    currents: CurrentFields
    metadata: dict[str, Any]


@dataclass(frozen=True)
class RelaxationResult:
    """Final state and compact history for one stationary relaxation run."""

    state: GTDGLStationaryState
    history: dict[str, np.ndarray]
    summary: dict[str, Any]


def compute_current_fields(
    *,
    psi_J: np.ndarray,
    phi_V: np.ndarray,
    Te_K: np.ndarray,
    material: GTDGLMaterial,
    ops: FVOperators,
) -> CurrentFields:
    """Evaluate Usadel-like, GL, normal and total currents on edges and nodes."""
    psi = np.asarray(psi_J, dtype=np.complex128)
    phi = np.asarray(phi_V, dtype=float)
    Te = np.asarray(Te_K, dtype=float)

    R_node = np.abs(psi)
    R_edge = edge_average(R_node, ops)
    Te_edge = np.maximum(edge_average(Te, ops), 1.0e-12)

    Q_edge = edge_phase_gradient_from_psi(psi, ops)

    coeff_us = (
        math.pi
        * material.sigma_n_S_m
        / (2.0 * E_CHARGE_C)
        * R_edge
        * np.tanh(R_edge / (2.0 * K_B_J_K * Te_edge))
    )
    edge_js_us = coeff_us * Q_edge

    coeff_gl = (
        math.pi
        * material.sigma_n_S_m
        * R_edge**2
        / (4.0 * E_CHARGE_C * K_B_J_K * material.Tc_K)
    )
    edge_js_gl = coeff_gl * Q_edge

    grad_phi = edge_scalar_gradient(phi, ops)
    edge_jn = -material.sigma_n_S_m * grad_phi
    edge_jtot = edge_js_us + edge_jn

    div_us = divergence_from_edge_scalar(edge_js_us, ops)
    div_gl = divergence_from_edge_scalar(edge_js_gl, ops)
    div_total = divergence_from_edge_scalar(edge_jtot, ops)

    js_x, js_y = edge_scalar_to_node_vector(edge_js_us, ops)
    jn_x, jn_y = edge_scalar_to_node_vector(edge_jn, ops)
    jt_x, jt_y = edge_scalar_to_node_vector(edge_jtot, ops)

    return CurrentFields(
        edge_Q_m_inv=Q_edge,
        edge_js_us_A_m2=edge_js_us,
        edge_js_gl_A_m2=edge_js_gl,
        edge_jn_A_m2=edge_jn,
        edge_jtot_A_m2=edge_jtot,
        node_div_js_us_A_m3=div_us,
        node_div_js_gl_A_m3=div_gl,
        node_div_jtot_A_m3=div_total,
        node_js_us_x_A_m2=js_x,
        node_js_us_y_A_m2=js_y,
        node_jn_x_A_m2=jn_x,
        node_jn_y_A_m2=jn_y,
        node_jtot_x_A_m2=jt_x,
        node_jtot_y_A_m2=jt_y,
    )


def solve_poisson_potential(
    *,
    edge_js_us_A_m2: np.ndarray,
    material: GTDGLMaterial,
    ops: FVOperators,
) -> np.ndarray:
    """Solve sigma_n Laplacian(phi) = div(j_s^Us) with mean-zero gauge."""
    js = np.asarray(edge_js_us_A_m2, dtype=float)
    conductance = material.sigma_n_S_m * ops.dual_face_length_m / ops.edge_length_m
    i = ops.edge_i
    j = ops.edge_j

    rhs = np.zeros(ops.n_nodes, dtype=float)
    edge_flux = ops.dual_face_length_m * js
    np.add.at(rhs, i, -edge_flux)
    np.add.at(rhs, j, edge_flux)

    if coo_matrix is not None and spsolve is not None:
        rows = np.concatenate([i, i, j, j])
        cols = np.concatenate([i, j, j, i])
        data = np.concatenate([conductance, -conductance, conductance, -conductance])
        matrix = coo_matrix((data, (rows, cols)), shape=(ops.n_nodes, ops.n_nodes)).tolil()
        matrix[0, :] = 0.0
        matrix[0, 0] = 1.0
        rhs[0] = 0.0
        phi = np.asarray(spsolve(matrix.tocsr(), rhs), dtype=float)
    else:
        matrix = np.zeros((ops.n_nodes, ops.n_nodes), dtype=float)
        for a, b, g in zip(i, j, conductance):
            matrix[a, a] += g
            matrix[a, b] -= g
            matrix[b, b] += g
            matrix[b, a] -= g
        matrix[0, :] = 0.0
        matrix[0, 0] = 1.0
        rhs[0] = 0.0
        phi = np.linalg.solve(matrix, rhs)

    phi -= float(np.mean(phi))
    return phi


def gtdgl_forcing(
    *,
    psi_J: np.ndarray,
    Te_K: np.ndarray,
    currents: CurrentFields,
    material: GTDGLMaterial,
    ops: FVOperators,
    regularization_fraction: float = 1.0e-9,
) -> np.ndarray:
    """Evaluate the explicit nonlinear forcing F[Delta]."""
    psi = np.asarray(psi_J, dtype=np.complex128)
    Te = np.asarray(Te_K, dtype=float)
    R2 = np.abs(psi) ** 2

    xi2 = material.xi_mod_squared_m2(Te)
    delta_mod2 = material.delta_mod_squared_J2(Te)
    local = 1.0 - Te / material.Tc_K - R2 / delta_mod2

    eps2 = (regularization_fraction * material.delta0_J) ** 2
    R2_safe = np.maximum(R2, eps2)
    div_diff = currents.node_div_js_us_A_m3 - currents.node_div_js_gl_A_m3
    correction = 1j * material.allmaras_C(Te) * div_diff * psi / R2_safe

    return xi2 * laplacian(psi, ops) + local * psi + correction


def kwt_local_update(
    *,
    psi_J: np.ndarray,
    phi_V: np.ndarray,
    Te_K: np.ndarray,
    forcing_J: np.ndarray,
    dt_s: float,
    material: GTDGLMaterial,
    discriminant_tol: float = 1.0e-12,
) -> tuple[np.ndarray, bool, float]:
    """Advance Delta by one local semi-implicit KWT step."""
    psi = np.asarray(psi_J, dtype=np.complex128)
    phi = np.asarray(phi_V, dtype=float)
    Te = np.asarray(Te_K, dtype=float)
    F = np.asarray(forcing_J, dtype=np.complex128)

    R2_old = np.abs(psi) ** 2
    rho = material.rho_kwt(Te, np.sqrt(R2_old))
    alpha = material.alpha_kwt_J_inv2(Te)

    U = np.exp(-1j * (2.0 * E_CHARGE_C / HBAR_J_S) * phi * dt_s)
    U_conj = np.conjugate(U)

    z = alpha * U_conj * psi
    w = U_conj * (psi + alpha * psi * R2_old + dt_s * F / (material.tau0_GL_s * rho))

    abs_z2 = np.abs(z) ** 2
    abs_w2 = np.abs(w) ** 2
    B = 1.0 + 2.0 * np.real(w * np.conjugate(z))
    disc = B**2 - 4.0 * abs_z2 * abs_w2
    min_disc = float(np.min(disc))

    scale = np.maximum(B**2 + 4.0 * abs_z2 * abs_w2, 1.0)
    bad = disc < -discriminant_tol * scale
    if np.any(bad):
        return psi.copy(), False, min_disc

    disc = np.maximum(disc, 0.0)
    denom = B + np.sqrt(disc)
    tiny = np.abs(denom) < 1.0e-300
    r = np.empty_like(abs_w2, dtype=float)
    r[~tiny] = 2.0 * abs_w2[~tiny] / denom[~tiny]
    r[tiny] = abs_w2[tiny]

    X = w - z * r
    psi_new = U * X
    return psi_new, True, min_disc


def terminal_node_mask(mesh) -> np.ndarray:
    """Return boolean mask for left and right terminal nodes."""
    nodes = np.asarray(mesh.nodes, dtype=float)
    x = nodes[:, 0]
    tol = max(1.0e-15, 1.0e-9 * float(mesh.length_m))
    return (np.abs(x - np.min(x)) <= tol) | (np.abs(x - np.max(x)) <= tol)


def apply_terminal_dirichlet(
    *,
    psi_trial_J: np.ndarray,
    psi_reference_J: np.ndarray,
    mesh,
    enabled: bool = True,
) -> np.ndarray:
    """Pin the complex order parameter on the longitudinal terminals."""
    psi = np.array(psi_trial_J, dtype=np.complex128, copy=True)
    if enabled:
        mask = terminal_node_mask(mesh)
        psi[mask] = np.asarray(psi_reference_J, dtype=np.complex128)[mask]
    return psi


def relax_stationary_gtdgl(
    *,
    mesh,
    edge_data,
    seed,
    material: GTDGLMaterial,
    ops: FVOperators,
    steps: int = 2000,
    dt_s: float = 2.5e-16,
    min_steps: int = 10,
    tolerance_eta: float = 1.0e-9,
    tolerance_current_residual: float = 1.0e-6,
    eta_reject: float = 5.0e-2,
    adapt_dt: bool = True,
    dt_min_s: float = 1.0e-18,
    dt_max_s: float = 2.0e-15,
    lock_terminals: bool = True,
) -> RelaxationResult:
    """Relax the OE6 seed with frozen temperatures and gTDGL/Poisson active."""
    if steps <= 0:
        raise ValueError("steps must be positive.")
    if dt_s <= 0.0:
        raise ValueError("dt_s must be positive.")

    psi0 = np.asarray(seed.node_psi_real_J, dtype=float) + 1j * np.asarray(seed.node_psi_imag_J, dtype=float)
    psi = psi0.copy()
    phi = np.asarray(seed.node_phi_electric_V, dtype=float).copy()
    Te = np.asarray(seed.node_Te_K, dtype=float).copy()
    Tph = np.asarray(seed.node_Tph_K, dtype=float).copy()

    t_s = 0.0
    accepted = 0
    rejected = 0
    converged = False

    hist_t: list[float] = []
    hist_dt: list[float] = []
    hist_eta: list[float] = []
    hist_res: list[float] = []
    hist_v: list[float] = []
    hist_ir: list[float] = []
    hist_il: list[float] = []

    currents = compute_current_fields(psi_J=psi, phi_V=phi, Te_K=Te, material=material, ops=ops)

    for _ in range(int(steps)):
        forcing = gtdgl_forcing(
            psi_J=psi,
            Te_K=Te,
            currents=currents,
            material=material,
            ops=ops,
        )
        psi_trial, ok, min_disc = kwt_local_update(
            psi_J=psi,
            phi_V=phi,
            Te_K=Te,
            forcing_J=forcing,
            dt_s=dt_s,
            material=material,
        )

        if not ok:
            rejected += 1
            if adapt_dt and dt_s > dt_min_s:
                dt_s = max(dt_min_s, 0.5 * dt_s)
                continue
            raise FloatingPointError(f"KWT update failed; min discriminant={min_disc:.6e}")

        psi_trial = apply_terminal_dirichlet(
            psi_trial_J=psi_trial,
            psi_reference_J=psi0,
            mesh=mesh,
            enabled=lock_terminals,
        )

        trial_currents_no_phi = compute_current_fields(
            psi_J=psi_trial,
            phi_V=np.zeros_like(phi),
            Te_K=Te,
            material=material,
            ops=ops,
        )
        phi_trial = solve_poisson_potential(
            edge_js_us_A_m2=trial_currents_no_phi.edge_js_us_A_m2,
            material=material,
            ops=ops,
        )
        trial_currents = compute_current_fields(
            psi_J=psi_trial,
            phi_V=phi_trial,
            Te_K=Te,
            material=material,
            ops=ops,
        )

        eta = float(np.max(np.abs(np.abs(psi_trial) ** 2 - np.abs(psi) ** 2)) / material.delta0_J**2)
        if eta > eta_reject and adapt_dt and dt_s > dt_min_s:
            rejected += 1
            dt_s = max(dt_min_s, 0.5 * dt_s)
            continue

        psi = psi_trial
        phi = phi_trial
        currents = trial_currents
        t_s += dt_s
        accepted += 1

        residual = current_residual(currents, mesh)
        voltage = terminal_voltage(np.asarray(mesh.nodes, dtype=float), phi, length_m=float(mesh.length_m))
        boundary = boundary_currents_from_node_vectors(
            mesh=mesh,
            edge_data=edge_data,
            jx_A_m2=currents.node_jtot_x_A_m2,
            jy_A_m2=currents.node_jtot_y_A_m2,
            thickness_m=material.thickness_m,
        )

        hist_t.append(t_s)
        hist_dt.append(dt_s)
        hist_eta.append(eta)
        hist_res.append(residual)
        hist_v.append(voltage)
        hist_ir.append(boundary["right_A"])
        hist_il.append(boundary["left_A"])

        if accepted >= min_steps and eta < tolerance_eta and residual < tolerance_current_residual:
            converged = True
            break

        if adapt_dt and eta < 0.1 * tolerance_eta:
            dt_s = min(dt_max_s, 1.2 * dt_s)

    metadata = {
        "backend": "oe7_stationary_gtdgl_poisson_v1",
        "description": (
            "Frozen-temperature stationary gTDGL/Poisson relaxation from the OE6 "
            "analytic seed. External circuit and thermal evolution are inactive."
        ),
        "accepted_steps": int(accepted),
        "rejected_steps": int(rejected),
        "requested_steps": int(steps),
        "converged": bool(converged),
        "final_time_s": float(t_s),
        "tau_scale": float(material.tau_scale),
        "tau_ee_Tc_effective_s": float(material.tau_scale * material.tau_ee_Tc_s),
        "tau_ep_Tc_effective_s": float(material.tau_scale * material.tau_ep_Tc_s),
        "lock_terminals": bool(lock_terminals),
        "thermal_policy": "frozen_Te_Tph",
        "circuit_policy": "inactive",
        "poisson_gauge": "mean_zero",
    }

    state = GTDGLStationaryState(
        psi_J=psi,
        phi_V=phi,
        Te_K=Te,
        Tph_K=Tph,
        currents=currents,
        metadata=metadata,
    )
    history = {
        "t_s": np.asarray(hist_t, dtype=float),
        "dt_s": np.asarray(hist_dt, dtype=float),
        "eta_R": np.asarray(hist_eta, dtype=float),
        "current_residual": np.asarray(hist_res, dtype=float),
        "terminal_voltage_V": np.asarray(hist_v, dtype=float),
        "integrated_right_current_A": np.asarray(hist_ir, dtype=float),
        "integrated_left_current_A": np.asarray(hist_il, dtype=float),
    }
    summary = stationary_summary(
        mesh=mesh,
        edge_data=edge_data,
        state=state,
        material=material,
        history=history,
    )
    return RelaxationResult(state=state, history=history, summary=summary)


def current_residual(currents: CurrentFields, mesh) -> float:
    """Dimensionless current-continuity residual from Appendix B Eq. (177)."""
    div = np.asarray(currents.node_div_jtot_A_m3, dtype=float)
    jmag = np.sqrt(currents.node_jtot_x_A_m2**2 + currents.node_jtot_y_A_m2**2)
    javg = max(float(np.mean(jmag)), 1.0e-300)
    xi_mesh = max(float(getattr(mesh, "target_spacing_m", 1.0)), 1.0e-300)
    return float(np.sqrt(np.mean(div**2)) / (javg / xi_mesh))


def stationary_summary(
    *,
    mesh,
    edge_data,
    state: GTDGLStationaryState,
    material: GTDGLMaterial,
    history: dict[str, np.ndarray],
) -> dict[str, Any]:
    """Build a compact YAML/console summary for the stationary state."""
    R = np.abs(state.psi_J)
    theta = np.unwrap(np.angle(state.psi_J))
    div = state.currents.node_div_jtot_A_m3
    voltage = terminal_voltage(np.asarray(mesh.nodes, dtype=float), state.phi_V, length_m=float(mesh.length_m))
    boundary = boundary_currents_from_node_vectors(
        mesh=mesh,
        edge_data=edge_data,
        jx_A_m2=state.currents.node_jtot_x_A_m2,
        jy_A_m2=state.currents.node_jtot_y_A_m2,
        thickness_m=material.thickness_m,
    )

    final_eta = float(history["eta_R"][-1]) if history["eta_R"].size else float("nan")
    final_res = float(history["current_residual"][-1]) if history["current_residual"].size else current_residual(state.currents, mesh)

    return {
        "backend": state.metadata["backend"],
        "converged": bool(state.metadata["converged"]),
        "accepted_steps": int(state.metadata["accepted_steps"]),
        "rejected_steps": int(state.metadata["rejected_steps"]),
        "final_time_ps": float(state.metadata["final_time_s"] / 1.0e-12),
        "tau_scale": float(material.tau_scale),
        "tau_ee_Tc_original_ps": float(material.tau_ee_Tc_s / 1.0e-12),
        "tau_ep_Tc_original_ps": float(material.tau_ep_Tc_s / 1.0e-12),
        "tau_ee_Tc_effective_ps": float(material.tau_scale * material.tau_ee_Tc_s / 1.0e-12),
        "tau_ep_Tc_effective_ps": float(material.tau_scale * material.tau_ep_Tc_s / 1.0e-12),
        "terminal_voltage_V": float(voltage),
        "phi_min_V": float(np.min(state.phi_V)),
        "phi_max_V": float(np.max(state.phi_V)),
        "delta_min_meV": float(np.min(R) / 1.602176634e-22),
        "delta_max_meV": float(np.max(R) / 1.602176634e-22),
        "theta_min_rad": float(np.min(theta)),
        "theta_max_rad": float(np.max(theta)),
        "jx_mean_A_m2": float(np.mean(state.currents.node_jtot_x_A_m2)),
        "jy_mean_A_m2": float(np.mean(state.currents.node_jtot_y_A_m2)),
        "divergence_rms_A_m3": float(np.sqrt(np.mean(div**2))),
        "current_residual": float(final_res),
        "eta_R_final": float(final_eta),
        "boundary_currents_A": boundary,
        "thermal_policy": state.metadata["thermal_policy"],
        "circuit_policy": state.metadata["circuit_policy"],
    }


def save_stationary_state_npz(state: GTDGLStationaryState, path: str | Path) -> Path:
    """Save final stationary state to NPZ."""
    output = Path(path)
    output.parent.mkdir(parents=True, exist_ok=True)
    c = state.currents
    np.savez_compressed(
        output,
        psi_real_J=np.real(state.psi_J),
        psi_imag_J=np.imag(state.psi_J),
        delta_J=np.abs(state.psi_J),
        theta_wrapped_rad=np.angle(state.psi_J),
        phi_V=state.phi_V,
        Te_K=state.Te_K,
        Tph_K=state.Tph_K,
        edge_Q_m_inv=c.edge_Q_m_inv,
        edge_js_us_A_m2=c.edge_js_us_A_m2,
        edge_js_gl_A_m2=c.edge_js_gl_A_m2,
        edge_jn_A_m2=c.edge_jn_A_m2,
        edge_jtot_A_m2=c.edge_jtot_A_m2,
        node_div_js_us_A_m3=c.node_div_js_us_A_m3,
        node_div_js_gl_A_m3=c.node_div_js_gl_A_m3,
        node_div_jtot_A_m3=c.node_div_jtot_A_m3,
        node_js_us_x_A_m2=c.node_js_us_x_A_m2,
        node_js_us_y_A_m2=c.node_js_us_y_A_m2,
        node_jn_x_A_m2=c.node_jn_x_A_m2,
        node_jn_y_A_m2=c.node_jn_y_A_m2,
        node_jtot_x_A_m2=c.node_jtot_x_A_m2,
        node_jtot_y_A_m2=c.node_jtot_y_A_m2,
        metadata=np.array(state.metadata, dtype=object),
    )
    return output


def save_relaxation_history_npz(history: dict[str, np.ndarray], path: str | Path) -> Path:
    """Save stationary relaxation history to NPZ."""
    output = Path(path)
    output.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(output, **history)
    return output