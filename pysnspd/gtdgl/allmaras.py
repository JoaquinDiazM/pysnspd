"""Appendix-B Allmaras diagnostics for the pyTDGL-like backend.

This module implements Appendix-B coefficient fields, the
current-divergence correction, and the dimensionless forcing injected into the
promoted flat gTDGL solver.  The physical Usadel/Matsubara current catalogue is
owned by ``usadel_current.py``; when available, its edge current is supplied to
this module so the Allmaras correction uses the same supercurrent closure as the
Poisson projection.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from pysnspd.gtdgl.material import GTDGLMaterial, K_B_J_K
from pysnspd.gtdgl.operators import (
    FVOperators,
    divergence_from_edge_scalar,
    edge_average,
    edge_phase_gradient_from_psi,
    edge_scalar_to_node_vector_least_squares,
)

HBAR_J_S = 1.054571817e-34
E_ABS_C = 1.602176634e-19
MEV_J = 1.602176634e-22


@dataclass(frozen=True)
class AllmarasCoefficientFields:
    """Node-local Appendix-B material coefficients."""

    tau_sc_s: np.ndarray
    rho_kwt: np.ndarray
    gamma_kwt_dimensionless: np.ndarray
    alpha_kwt_inv_J2: np.ndarray
    xi_mod2_m2: np.ndarray
    delta_mod2_J2: np.ndarray
    delta_mod_over_delta0: np.ndarray
    solver_epsilon: np.ndarray
    correction_C_J_m3_A: np.ndarray


@dataclass(frozen=True)
class AllmarasDiagnosticFields:
    """Snapshot diagnostic fields for the Appendix-B correction."""

    coefficients: AllmarasCoefficientFields
    edge_Q_m_inv: np.ndarray
    edge_R_J: np.ndarray
    edge_Te_K: np.ndarray
    edge_js_us_allmaras_A_m2: np.ndarray
    edge_js_gl_allmaras_A_m2: np.ndarray
    node_js_us_allmaras_x_A_m2: np.ndarray
    node_js_us_allmaras_y_A_m2: np.ndarray
    node_js_gl_allmaras_x_A_m2: np.ndarray
    node_js_gl_allmaras_y_A_m2: np.ndarray
    node_div_js_us_allmaras_A_m3: np.ndarray
    node_div_js_gl_allmaras_A_m3: np.ndarray
    node_mismatch_divergence_A_m3: np.ndarray
    node_phase_drive_complex_J: np.ndarray
    node_phase_drive_abs_J: np.ndarray
    node_phase_drive_abs_over_delta0: np.ndarray
    bulk_node_mask: np.ndarray



@dataclass(frozen=True)
class AllmarasForcingFields:
    """Appendix-B forcing fields used by the local KWT update.

    ``forcing_dimensionless`` is the complete bracket in Eq. (139), divided by
    ``Delta0`` and written in the dimensionless coordinates of the promoted
    pyTDGL-like solver.  The local solver then multiplies it by
    ``rho_KWT * dt/tau0`` through the existing ``w_i, z_i`` algebra.
    """

    coefficients: AllmarasCoefficientFields
    forcing_dimensionless: np.ndarray
    diffusion_dimensionless: np.ndarray
    reaction_dimensionless: np.ndarray
    phase_drive_dimensionless: np.ndarray
    edge_Q_m_inv: np.ndarray
    edge_R_J: np.ndarray
    edge_Te_K: np.ndarray
    edge_js_us_A_m2: np.ndarray
    edge_js_gl_A_m2: np.ndarray
    node_div_js_us_A_m3: np.ndarray
    node_div_js_gl_A_m3: np.ndarray
    node_mismatch_divergence_A_m3: np.ndarray


def allmaras_coefficients(
    *,
    psi_dimensionless: np.ndarray,
    material: GTDGLMaterial,
    Te_K: np.ndarray,
) -> AllmarasCoefficientFields:
    """Return Appendix-B local coefficient fields in SI units.

    The solver currently advances ``psi = Delta / Delta0``.  Therefore the KWT
    gamma used by the pyTDGL-like local algebra is

        gamma_i = 2 Delta0 tau_sc(Te_i) / hbar,

    so that ``sqrt(1 + gamma_i^2 |psi_i|^2)`` equals
    ``rho_KWT(Te_i, |Delta_i|)``.
    """

    psi = np.asarray(psi_dimensionless, dtype=np.complex128)
    Te = np.asarray(Te_K, dtype=float).reshape(-1)
    if psi.shape != Te.shape:
        raise ValueError(f"psi and Te must have the same node shape, got {psi.shape} and {Te.shape}.")
    Te_safe = np.maximum(Te, 1.0e-12)
    R_J = np.abs(psi) * float(material.delta0_J)

    tau_sc = np.asarray(material.tau_sc_s(Te_safe), dtype=float)
    rho = np.asarray(material.rho_kwt(Te_safe, np.maximum(R_J, 0.0)), dtype=float)
    gamma = 2.0 * float(material.delta0_J) * tau_sc / HBAR_J_S
    alpha = 2.0 * tau_sc * tau_sc / (HBAR_J_S * HBAR_J_S)

    delta_mod2 = np.asarray(material.delta_mod_squared_J2(Te_safe), dtype=float)
    delta0 = float(material.delta0_J)
    solver_epsilon = np.clip(delta_mod2 / max(delta0 * delta0, 1.0e-300), 0.0, 1.0)
    delta_mod_over_delta0 = np.sqrt(np.maximum(delta_mod2, 0.0)) / max(delta0, 1.0e-300)

    xi_mod2 = (
        np.pi
        * HBAR_J_S
        * float(material.D_m2_s)
        / (
            4.0
            * np.sqrt(2.0)
            * K_B_J_K
            * float(material.Tc_K)
            * np.sqrt(np.maximum(1.0 + Te_safe / float(material.Tc_K), 1.0e-300))
        )
    )

    correction_C = (
        HBAR_J_S
        * E_ABS_C
        * float(material.D_m2_s)
        / (
            float(material.sigma_n_S_m)
            * np.sqrt(2.0)
            * np.sqrt(np.maximum(1.0 + Te_safe / float(material.Tc_K), 1.0e-300))
        )
    )

    return AllmarasCoefficientFields(
        tau_sc_s=tau_sc,
        rho_kwt=rho,
        gamma_kwt_dimensionless=gamma,
        alpha_kwt_inv_J2=alpha,
        xi_mod2_m2=np.asarray(xi_mod2, dtype=float),
        delta_mod2_J2=delta_mod2,
        delta_mod_over_delta0=delta_mod_over_delta0,
        solver_epsilon=solver_epsilon,
        correction_C_J_m3_A=np.asarray(correction_C, dtype=float),
    )


def compute_allmaras_appendix_b_diagnostic(
    *,
    psi_dimensionless: np.ndarray,
    material: GTDGLMaterial,
    Te_K: np.ndarray,
    ops: FVOperators,
    terminal_node_mask: np.ndarray | None = None,
    blocked_edge_mask: np.ndarray | None = None,
    edge_js_usadel_A_m2: np.ndarray | None = None,
    bulk_guard_layers: int = 1,
    r_epsilon_fraction: float = 1.0e-6,
) -> AllmarasDiagnosticFields:
    """Compute the Appendix-B current-divergence correction as a diagnostic.

    This evaluates Eqs. (137), (159), and (160) from Appendix B on the current
    snapshot.  The result is diagnostic only: it is not yet injected into the
    ``w_i, z_i`` local update.
    """

    psi = np.asarray(psi_dimensionless, dtype=np.complex128)
    Te = np.asarray(Te_K, dtype=float).reshape(-1)
    coeff = allmaras_coefficients(psi_dimensionless=psi, material=material, Te_K=Te)

    Q = edge_phase_gradient_from_psi(psi, ops)
    R_node_J = np.abs(psi) * float(material.delta0_J)
    R_edge_J = edge_average(R_node_J, ops)
    Te_edge = np.maximum(edge_average(Te, ops), 1.0e-12)

    if edge_js_usadel_A_m2 is None:
        # Analytic Allmaras fallback retained for diagnostics and tests.  The
        # production usadel_poisson path supplies the Matsubara-table current
        # through ``edge_js_usadel_A_m2``.
        js_us = analytic_allmaras_usadel_current_edges(
            edge_Q_m_inv=Q,
            edge_R_J=R_edge_J,
            edge_Te_K=Te_edge,
            material=material,
        )
    else:
        js_us = np.asarray(edge_js_usadel_A_m2, dtype=float).reshape(-1)
        if js_us.shape != (ops.n_edges,):
            raise ValueError(f"edge_js_usadel_A_m2 must have shape ({ops.n_edges},), got {js_us.shape}.")
    js_gl = gl_supercurrent_edges(
        edge_Q_m_inv=Q,
        edge_R_J=R_edge_J,
        material=material,
    )

    if blocked_edge_mask is not None:
        mask = np.asarray(blocked_edge_mask, dtype=bool).reshape(-1)
        if mask.size != ops.n_edges:
            raise ValueError(f"blocked_edge_mask has length {mask.size}, expected {ops.n_edges}.")
        if np.any(mask):
            js_us = js_us.copy()
            js_gl = js_gl.copy()
            js_us[mask] = 0.0
            js_gl[mask] = 0.0

    js_us_x, js_us_y = edge_scalar_to_node_vector_least_squares(js_us, ops)
    js_gl_x, js_gl_y = edge_scalar_to_node_vector_least_squares(js_gl, ops)
    div_us = divergence_from_edge_scalar(js_us, ops)
    div_gl = divergence_from_edge_scalar(js_gl, ops)
    mismatch = div_us - div_gl

    terminal_mask = (
        np.zeros(ops.n_nodes, dtype=bool)
        if terminal_node_mask is None
        else np.asarray(terminal_node_mask, dtype=bool).reshape(-1).copy()
    )
    if terminal_mask.size != ops.n_nodes:
        raise ValueError(f"terminal_node_mask has length {terminal_mask.size}, expected {ops.n_nodes}.")
    bulk_mask = bulk_node_mask_from_terminal_mask(ops, terminal_mask, guard_layers=bulk_guard_layers)

    Delta_J = psi * float(material.delta0_J)
    r_eps2 = np.maximum(
        np.abs(Delta_J) ** 2,
        (float(r_epsilon_fraction) * float(material.delta0_J)) ** 2,
    )
    drive_complex = 1j * coeff.correction_C_J_m3_A * mismatch * Delta_J / r_eps2
    drive_abs = np.abs(drive_complex)

    return AllmarasDiagnosticFields(
        coefficients=coeff,
        edge_Q_m_inv=Q,
        edge_R_J=R_edge_J,
        edge_Te_K=Te_edge,
        edge_js_us_allmaras_A_m2=np.asarray(js_us, dtype=float),
        edge_js_gl_allmaras_A_m2=np.asarray(js_gl, dtype=float),
        node_js_us_allmaras_x_A_m2=np.asarray(js_us_x, dtype=float),
        node_js_us_allmaras_y_A_m2=np.asarray(js_us_y, dtype=float),
        node_js_gl_allmaras_x_A_m2=np.asarray(js_gl_x, dtype=float),
        node_js_gl_allmaras_y_A_m2=np.asarray(js_gl_y, dtype=float),
        node_div_js_us_allmaras_A_m3=np.asarray(div_us, dtype=float),
        node_div_js_gl_allmaras_A_m3=np.asarray(div_gl, dtype=float),
        node_mismatch_divergence_A_m3=np.asarray(mismatch, dtype=float),
        node_phase_drive_complex_J=np.asarray(drive_complex, dtype=np.complex128),
        node_phase_drive_abs_J=np.asarray(drive_abs, dtype=float),
        node_phase_drive_abs_over_delta0=np.asarray(drive_abs / max(float(material.delta0_J), 1.0e-300), dtype=float),
        bulk_node_mask=bulk_mask,
    )



def analytic_allmaras_usadel_current_edges(
    *,
    edge_Q_m_inv: np.ndarray,
    edge_R_J: np.ndarray,
    edge_Te_K: np.ndarray,
    material: GTDGLMaterial,
) -> np.ndarray:
    """Analytic Allmaras current closure on edges.

    This is kept as a fallback diagnostic.  The production path for the present
    thesis uses the PRE Matsubara/Usadel table instead.
    """

    Q = np.asarray(edge_Q_m_inv, dtype=float)
    R = np.asarray(edge_R_J, dtype=float)
    Te = np.maximum(np.asarray(edge_Te_K, dtype=float), 1.0e-12)
    return (
        np.pi
        * float(material.sigma_n_S_m)
        / (2.0 * E_ABS_C)
        * R
        * np.tanh(R / (2.0 * K_B_J_K * Te))
        * Q
    )


def gl_supercurrent_edges(
    *,
    edge_Q_m_inv: np.ndarray,
    edge_R_J: np.ndarray,
    material: GTDGLMaterial,
) -> np.ndarray:
    """Auxiliary GL supercurrent used only in the Allmaras mismatch."""

    Q = np.asarray(edge_Q_m_inv, dtype=float)
    R = np.asarray(edge_R_J, dtype=float)
    return (
        np.pi
        * float(material.sigma_n_S_m)
        * R
        * R
        / (4.0 * E_ABS_C * K_B_J_K * float(material.Tc_K))
        * Q
    )


def compute_allmaras_forcing_dimensionless(
    *,
    psi_dimensionless: np.ndarray,
    psi_laplacian_dimensionless: np.ndarray,
    material: GTDGLMaterial,
    Te_K: np.ndarray,
    ops: FVOperators,
    length_scale_m: float,
    edge_js_usadel_A_m2: np.ndarray | None = None,
    blocked_edge_mask: np.ndarray | None = None,
    r_epsilon_fraction: float = 1.0e-6,
) -> AllmarasForcingFields:
    """Return the Appendix-B forcing injected into the local KWT update.

    The continuous equation used here is

        tau0 [D_t Delta + alpha_KWT d_t |Delta|^2 Delta]
        = rho_KWT F[Delta],

    with

        F/Delta0 = (xi_mod^2/L0^2) L' psi
                 + [1 - Te/Tc - |Delta|^2/Delta_mod^2] psi
                 + i C (div j_s^Us - div j_s^GL) Delta/(R_eps^2 Delta0).

    ``L'`` is the promoted pyTDGL-like Laplacian in dimensionless coordinates
    x' = x/L0.  The returned array is dimensionless; the solver multiplies it by
    ``rho_KWT * dt`` because its time variable is already scaled by ``tau0``.
    """

    psi = np.asarray(psi_dimensionless, dtype=np.complex128).reshape(-1)
    lap = np.asarray(psi_laplacian_dimensionless, dtype=np.complex128).reshape(-1)
    Te = np.asarray(Te_K, dtype=float).reshape(-1)
    if psi.shape != (ops.n_nodes,):
        raise ValueError(f"psi_dimensionless must have shape ({ops.n_nodes},), got {psi.shape}.")
    if lap.shape != (ops.n_nodes,):
        raise ValueError(f"psi_laplacian_dimensionless must have shape ({ops.n_nodes},), got {lap.shape}.")
    if Te.shape != (ops.n_nodes,):
        raise ValueError(f"Te_K must have shape ({ops.n_nodes},), got {Te.shape}.")

    L0 = float(length_scale_m)
    if not np.isfinite(L0) or L0 <= 0.0:
        raise ValueError("length_scale_m must be positive and finite.")

    coeff = allmaras_coefficients(psi_dimensionless=psi, material=material, Te_K=Te)
    Q = edge_phase_gradient_from_psi(psi, ops)
    Delta_J = psi * float(material.delta0_J)
    R_node_J = np.abs(Delta_J)
    R_edge_J = edge_average(R_node_J, ops)
    Te_edge = np.maximum(edge_average(Te, ops), 1.0e-12)

    if edge_js_usadel_A_m2 is None:
        js_us = analytic_allmaras_usadel_current_edges(
            edge_Q_m_inv=Q,
            edge_R_J=R_edge_J,
            edge_Te_K=Te_edge,
            material=material,
        )
    else:
        js_us = np.asarray(edge_js_usadel_A_m2, dtype=float).reshape(-1)
        if js_us.shape != (ops.n_edges,):
            raise ValueError(f"edge_js_usadel_A_m2 must have shape ({ops.n_edges},), got {js_us.shape}.")
    js_gl = gl_supercurrent_edges(edge_Q_m_inv=Q, edge_R_J=R_edge_J, material=material)

    if blocked_edge_mask is not None:
        mask = np.asarray(blocked_edge_mask, dtype=bool).reshape(-1)
        if mask.size != ops.n_edges:
            raise ValueError(f"blocked_edge_mask has length {mask.size}, expected {ops.n_edges}.")
        if np.any(mask):
            js_us = js_us.copy()
            js_gl = js_gl.copy()
            js_us[mask] = 0.0
            js_gl[mask] = 0.0

    div_us = divergence_from_edge_scalar(js_us, ops)
    div_gl = divergence_from_edge_scalar(js_gl, ops)
    mismatch = div_us - div_gl

    diffusion = (coeff.xi_mod2_m2 / (L0 * L0)) * lap

    delta_mod2 = np.asarray(coeff.delta_mod2_J2, dtype=float)
    delta0 = float(material.delta0_J)
    delta_mod_floor2 = (float(r_epsilon_fraction) * delta0) ** 2
    denom = np.maximum(delta_mod2, delta_mod_floor2)
    reaction = (1.0 - Te / float(material.Tc_K) - (np.abs(Delta_J) ** 2) / denom) * psi

    r_eps2 = np.maximum(np.abs(Delta_J) ** 2, delta_mod_floor2)
    phase_drive = 1j * coeff.correction_C_J_m3_A * mismatch * Delta_J / r_eps2

    forcing = diffusion + reaction + phase_drive
    forcing = np.asarray(forcing, dtype=np.complex128)
    forcing[~np.isfinite(forcing)] = 0.0

    return AllmarasForcingFields(
        coefficients=coeff,
        forcing_dimensionless=forcing,
        diffusion_dimensionless=np.asarray(diffusion, dtype=np.complex128),
        reaction_dimensionless=np.asarray(reaction, dtype=np.complex128),
        phase_drive_dimensionless=np.asarray(phase_drive, dtype=np.complex128),
        edge_Q_m_inv=Q,
        edge_R_J=R_edge_J,
        edge_Te_K=Te_edge,
        edge_js_us_A_m2=np.asarray(js_us, dtype=float),
        edge_js_gl_A_m2=np.asarray(js_gl, dtype=float),
        node_div_js_us_A_m3=np.asarray(div_us, dtype=float),
        node_div_js_gl_A_m3=np.asarray(div_gl, dtype=float),
        node_mismatch_divergence_A_m3=np.asarray(mismatch, dtype=float),
    )


def bulk_node_mask_from_terminal_mask(
    ops: FVOperators,
    terminal_node_mask: np.ndarray,
    *,
    guard_layers: int = 1,
) -> np.ndarray:
    """Return a bulk diagnostic mask excluding terminals and guard layers.

    For production meshes, one or more graph guard layers remove the nodes
    adjacent to the normal-metal contacts, preventing contact singularities from
    dominating bulk diagnostics.  Very small unit-test meshes may not contain an
    interior after one guard expansion.  In that case the physically meaningful
    fallback is the terminal-only exclusion, not an empty bulk mask.
    """

    terminal = np.asarray(terminal_node_mask, dtype=bool).reshape(-1).copy()
    if terminal.size != ops.n_nodes:
        raise ValueError(f"terminal_node_mask has length {terminal.size}, expected {ops.n_nodes}.")
    base_bulk = ~terminal
    if not np.any(base_bulk):
        return np.ones(ops.n_nodes, dtype=bool)

    excluded = terminal.copy()
    for _ in range(max(0, int(guard_layers))):
        touch = excluded[np.asarray(ops.edge_i, dtype=np.int64)] | excluded[np.asarray(ops.edge_j, dtype=np.int64)]
        if not np.any(touch):
            continue
        new_excluded = excluded.copy()
        new_excluded[np.asarray(ops.edge_i, dtype=np.int64)[touch]] = True
        new_excluded[np.asarray(ops.edge_j, dtype=np.int64)[touch]] = True
        excluded = new_excluded

    guarded_bulk = ~excluded
    if not np.any(guarded_bulk):
        return base_bulk
    return guarded_bulk


def rms(values: np.ndarray, mask: np.ndarray | None = None) -> float:
    """Finite-safe RMS helper for summaries."""

    arr = np.asarray(values, dtype=float)
    if mask is not None:
        m = np.asarray(mask, dtype=bool)
        arr = arr[m]
    arr = arr[np.isfinite(arr)]
    if arr.size == 0:
        return float("nan")
    return float(np.sqrt(np.mean(arr * arr)))


def max_abs(values: np.ndarray, mask: np.ndarray | None = None) -> float:
    """Finite-safe max-absolute helper for summaries."""

    arr = np.asarray(values, dtype=float)
    if mask is not None:
        m = np.asarray(mask, dtype=bool)
        arr = arr[m]
    arr = arr[np.isfinite(arr)]
    if arr.size == 0:
        return float("nan")
    return float(np.max(np.abs(arr)))
