"""Formula fields and current diagnostics for OE7 gTDGL."""
from __future__ import annotations

import math

import numpy as np

from pysnspd.gtdgl.material import E_CHARGE_C, K_B_J_K, HBAR_J_S, GTDGLMaterial
from pysnspd.gtdgl.operators import (
    FVOperators,
    divergence_from_edge_scalar,
    edge_average,
    edge_phase_gradient_from_psi,
    edge_scalar_gradient,
    edge_scalar_to_node_vector_least_squares,
    laplacian,
)
from pysnspd.gtdgl.state import CurrentFields, FormulaFields


def compute_current_fields(
    *,
    psi_J: np.ndarray,
    phi_V: np.ndarray,
    Te_K: np.ndarray,
    material: GTDGLMaterial,
    ops: FVOperators,
    boundary_accum_A_m: np.ndarray | None = None,
) -> CurrentFields:
    """Evaluate Usadel-like, GL, normal and total current fields.

    Edge currents follow the notebook formulas. Node vectors are diagnostics
    reconstructed from edge projections by the notebook LS method.
    """
    defs = compute_formula_fields(
        psi_J=psi_J,
        Te_K=Te_K,
        material=material,
        ops=ops,
    )
    phi = np.asarray(phi_V, dtype=float)
    grad_phi = edge_scalar_gradient(phi, ops)
    edge_jn = -material.sigma_n_S_m * grad_phi
    edge_jtot = defs.edge_js_us_A_m2 + edge_jn

    div_jn = divergence_from_edge_scalar(edge_jn, ops)
    div_total = divergence_from_edge_scalar(
        edge_jtot,
        ops,
        boundary_accum_A_m=boundary_accum_A_m,
    )

    js_x, js_y = edge_scalar_to_node_vector_least_squares(defs.edge_js_us_A_m2, ops)
    jn_x, jn_y = edge_scalar_to_node_vector_least_squares(edge_jn, ops)
    jt_x, jt_y = edge_scalar_to_node_vector_least_squares(edge_jtot, ops)

    edge_pb = pairbreaking_ratio_edges(
        Q_edge_m_inv=defs.edge_Q_m_inv,
        Te_edge_K=edge_average(Te_K, ops),
        material=material,
    )
    node_pb = edge_to_node_weighted_average(edge_pb, ops)

    return CurrentFields(
        edge_Q_m_inv=defs.edge_Q_m_inv,
        edge_js_us_A_m2=defs.edge_js_us_A_m2,
        edge_js_gl_A_m2=defs.edge_js_gl_A_m2,
        edge_jn_A_m2=edge_jn,
        edge_jtot_A_m2=edge_jtot,
        node_div_js_us_A_m3=defs.node_div_js_us_A_m3,
        node_div_js_gl_A_m3=defs.node_div_js_gl_A_m3,
        node_div_jtot_A_m3=div_total,
        node_js_us_x_A_m2=js_x,
        node_js_us_y_A_m2=js_y,
        node_jn_x_A_m2=jn_x,
        node_jn_y_A_m2=jn_y,
        node_jtot_x_A_m2=jt_x,
        node_jtot_y_A_m2=jt_y,
        edge_pairbreaking_ratio=edge_pb,
        node_pairbreaking_ratio=node_pb,
    )


def compute_formula_fields(
    *,
    psi_J: np.ndarray,
    Te_K: np.ndarray,
    material: GTDGLMaterial,
    ops: FVOperators,
) -> FormulaFields:
    """Compute the notebook ``defs`` object from the current Delta field."""
    psi = np.asarray(psi_J, dtype=np.complex128)
    Te = np.asarray(Te_K, dtype=float)
    R = np.abs(psi)
    R_safe = safe_abs_delta(material, R)

    tau_sc = material.tau_sc_s(Te)
    rho = material.rho_kwt(Te, R_safe)
    alpha = 2.0 * tau_sc**2 / HBAR_J_S**2
    xi2 = material.xi_mod_squared_m2(Te)
    xi = np.sqrt(np.maximum(xi2, 0.0))
    delta_mod2 = material.delta_mod_squared_J2(Te)
    delta_mod = np.sqrt(np.maximum(delta_mod2, 0.0))

    lap_delta = laplacian(psi, ops)
    js_us, Q_edge = edge_supercurrent_usadel(
        psi_J=psi,
        Te_K=Te,
        material=material,
        ops=ops,
    )
    js_gl = edge_supercurrent_gl(
        psi_J=psi,
        material=material,
        ops=ops,
    )

    div_us = divergence_from_edge_scalar(js_us, ops)
    div_gl = divergence_from_edge_scalar(js_gl, ops)
    div_corr = div_us - div_gl

    C = material.allmaras_C(Te)
    reaction = (1.0 - Te / material.Tc_K - R**2 / delta_mod2) * psi
    correction = 1j * C * div_corr * psi / (R_safe**2)
    forcing = xi2 * lap_delta + reaction + correction

    return FormulaFields(
        Te_K=Te,
        rho=rho,
        tau_sc_s=tau_sc,
        alpha_kwt_J_inv2=alpha,
        xi_mod_m=xi,
        delta_mod_J=delta_mod,
        delta_abs_J=R,
        edge_Q_m_inv=Q_edge,
        edge_js_us_A_m2=js_us,
        edge_js_gl_A_m2=js_gl,
        node_div_js_us_A_m3=div_us,
        node_div_js_gl_A_m3=div_gl,
        node_div_correction_A_m3=div_corr,
        node_lap_delta_J_m2=lap_delta,
        forcing_J=forcing,
    )


def edge_supercurrent_usadel(
    *,
    psi_J: np.ndarray,
    Te_K: np.ndarray,
    material: GTDGLMaterial,
    ops: FVOperators,
) -> tuple[np.ndarray, np.ndarray]:
    """Notebook Usadel-like edge supercurrent and Q."""
    R = np.abs(np.asarray(psi_J, dtype=np.complex128))
    R_e = edge_average(R, ops)
    Te_e = np.maximum(edge_average(Te_K, ops), 1.0e-12)
    Q_e = edge_phase_gradient_from_psi(psi_J, ops)
    pref = math.pi * material.sigma_n_S_m / (2.0 * E_CHARGE_C)
    js = pref * R_e * np.tanh(R_e / (2.0 * K_B_J_K * Te_e)) * Q_e
    return js, Q_e


def edge_supercurrent_gl(
    *,
    psi_J: np.ndarray,
    material: GTDGLMaterial,
    ops: FVOperators,
) -> np.ndarray:
    """Notebook GL edge supercurrent used only in the correction term."""
    R = np.abs(np.asarray(psi_J, dtype=np.complex128))
    R_e = edge_average(R, ops)
    Q_e = edge_phase_gradient_from_psi(psi_J, ops)
    pref = math.pi * material.sigma_n_S_m / (4.0 * E_CHARGE_C * K_B_J_K * material.Tc_K)
    return pref * R_e**2 * Q_e


def safe_abs_delta(material: GTDGLMaterial, R_J: np.ndarray) -> np.ndarray:
    """Notebook small floor for denominators only."""
    return np.maximum(np.asarray(R_J, dtype=float), 1.0e-10 * material.delta0_J)


def pairbreaking_ratio_edges(
    *,
    Q_edge_m_inv: np.ndarray,
    Te_edge_K: np.ndarray,
    material: GTDGLMaterial,
) -> np.ndarray:
    """Return xi^2 Q^2/(1 - T/Tc)."""
    Q = np.asarray(Q_edge_m_inv, dtype=float)
    Te = np.asarray(Te_edge_K, dtype=float)
    a = np.maximum(1.0 - Te / material.Tc_K, 1.0e-30)
    xi2 = material.xi_mod_squared_m2(Te)
    return np.asarray(xi2 * Q * Q / a, dtype=float)


def edge_to_node_weighted_average(edge_values: np.ndarray, ops: FVOperators) -> np.ndarray:
    """Average edge quantities to nodes using notebook dual/length weights."""
    values = np.asarray(edge_values, dtype=float)
    if values.shape != (ops.n_edges,):
        raise ValueError(f"edge_values must have shape ({ops.n_edges},), got {values.shape}.")
    weights = ops.dual_face_length_m / np.maximum(ops.edge_length_m, 1.0e-300)
    weights = np.maximum(weights, 1.0e-300)
    out = np.zeros(ops.n_nodes, dtype=float)
    wsum = np.zeros(ops.n_nodes, dtype=float)
    np.add.at(out, ops.edge_i, weights * values)
    np.add.at(out, ops.edge_j, weights * values)
    np.add.at(wsum, ops.edge_i, weights)
    np.add.at(wsum, ops.edge_j, weights)
    return out / np.maximum(wsum, 1.0e-300)

