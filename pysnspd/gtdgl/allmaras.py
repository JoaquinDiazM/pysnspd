"""Appendix-B Allmaras diagnostics for the pyTDGL-like backend.

This module implements Appendix-B coefficient fields and a diagnostic analytic
Allmaras current-divergence reference.  The active solver path is the
pyTDGL-like backend; the microscopic Usadel/Matsubara current catalogue is
handled by ``usadel_current.py`` and by the adapter-level diagnostics.  This
module still does not alter the ``X_i + z_i |X_i|^2 = w_i`` local algebraic
update.
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

    js_us = (
        np.pi
        * float(material.sigma_n_S_m)
        / (2.0 * E_ABS_C)
        * R_edge_J
        * np.tanh(R_edge_J / (2.0 * K_B_J_K * Te_edge))
        * Q
    )
    js_gl = (
        np.pi
        * float(material.sigma_n_S_m)
        * R_edge_J
        * R_edge_J
        / (4.0 * E_ABS_C * K_B_J_K * float(material.Tc_K))
        * Q
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
