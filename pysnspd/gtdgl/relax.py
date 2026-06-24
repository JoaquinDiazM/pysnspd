"""Stationary gTDGL/Poisson relaxation using the notebook OE7 ordering.

This file intentionally ports the solver order from ``gTDGL_model.ipynb`` while
keeping the pySNSPD package API used by ``pipelines/02_ss_run_template.py``:

1. apply stationary current-Neumann boundary constraints to Delta;
2. compute formula fields from the constrained state;
3. perform the local semi-implicit KWT quadratic update with temporal gauge link;
4. re-apply boundary constraints;
5. solve the conservative Poisson projection for varphi;
6. recompute fields and diagnostics;
7. adapt dt from max change in |Delta|^2, with KWT discriminant retries.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any
import json
import math

import numpy as np

try:
    from tqdm.auto import trange
except Exception:  # pragma: no cover
    trange = None

try:
    from scipy.sparse import coo_matrix, csr_matrix, bmat
    from scipy.sparse.linalg import splu, spsolve
except Exception:  # pragma: no cover
    coo_matrix = None
    csr_matrix = None
    bmat = None
    splu = None
    spsolve = None

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
    edge_flux_accumulator_A_m,
    edge_phase_gradient_from_psi,
    edge_scalar_gradient,
    edge_scalar_to_node_vector_least_squares,
    laplacian,
    terminal_boundary_accum_A_m,
    terminal_voltage,
)

MEV_J = 1.602176634e-22


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
    edge_pairbreaking_ratio: np.ndarray
    node_pairbreaking_ratio: np.ndarray


@dataclass(frozen=True)
class FormulaFields:
    """Notebook formula-field bundle evaluated from Delta, varphi and Te."""

    Te_K: np.ndarray
    rho: np.ndarray
    tau_sc_s: np.ndarray
    alpha_kwt_J_inv2: np.ndarray
    xi_mod_m: np.ndarray
    delta_mod_J: np.ndarray
    delta_abs_J: np.ndarray
    edge_Q_m_inv: np.ndarray
    edge_js_us_A_m2: np.ndarray
    edge_js_gl_A_m2: np.ndarray
    node_div_js_us_A_m3: np.ndarray
    node_div_js_gl_A_m3: np.ndarray
    node_div_correction_A_m3: np.ndarray
    node_lap_delta_J_m2: np.ndarray
    forcing_J: np.ndarray


@dataclass(frozen=True)
class PoissonResult:
    """Poisson projection result, following the notebook naming."""

    phi_V: np.ndarray
    edge_jn_A_m2: np.ndarray
    edge_jtot_A_m2: np.ndarray
    node_div_jtot_A_m3: np.ndarray
    node_div_js_A_m3: np.ndarray
    node_div_jn_A_m3: np.ndarray
    lambda_mean: float


@dataclass(frozen=True)
class StepInfo:
    """One accepted KWT/Poisson step diagnostic."""

    dt_eff_s: float
    retries: int
    discr_min: float
    max_amp2_change_rel: float
    max_Q_m_inv: float
    p95_Q_m_inv: float
    median_Q_m_inv: float
    max_js_A_m2: float
    max_j_A_m2: float
    rms_divj_rel: float
    max_divj_rel: float
    delta_min_norm: float
    delta_max_norm: float
    phi_ptp_V: float


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


@dataclass(frozen=True)
class _PoissonOperator:
    A_aug: Any
    solver: Any


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


def build_poisson_operator(*, material: GTDGLMaterial, ops: FVOperators) -> _PoissonOperator:
    """Build the notebook Poisson operator with mean-zero gauge row."""
    sigma = material.sigma_n_S_m
    g = sigma * ops.dual_face_length_m / ops.edge_length_m
    i = ops.edge_i.astype(np.int64)
    j = ops.edge_j.astype(np.int64)
    n = int(ops.n_nodes)
    rows = np.concatenate([i, i, j, j])
    cols = np.concatenate([i, j, j, i])
    data = np.concatenate([g, -g, g, -g])

    if coo_matrix is not None and csr_matrix is not None and bmat is not None:
        A = coo_matrix((data, (rows, cols)), shape=(n, n)).tocsr()
        ones_col = csr_matrix(np.ones((n, 1), dtype=float))
        ones_row = csr_matrix(np.ones((1, n), dtype=float))
        zero = csr_matrix((1, 1), dtype=float)
        A_aug = bmat([[A, ones_col], [ones_row, zero]], format="csc")
        solver = splu(A_aug) if splu is not None else None
        return _PoissonOperator(A_aug=A_aug, solver=solver)

    A = np.zeros((n, n), dtype=float)  # pragma: no cover
    for a, b, gg in zip(i, j, g):
        A[a, a] += gg
        A[a, b] -= gg
        A[b, b] += gg
        A[b, a] -= gg
    A_aug = np.zeros((n + 1, n + 1), dtype=float)
    A_aug[:n, :n] = A
    A_aug[:n, n] = 1.0
    A_aug[n, :n] = 1.0
    return _PoissonOperator(A_aug=A_aug, solver=None)


def solve_varphi_poisson(
    *,
    edge_js_us_A_m2: np.ndarray,
    material: GTDGLMaterial,
    ops: FVOperators,
    poisson_op: _PoissonOperator | None = None,
    boundary_accum_A_m: np.ndarray | None = None,
) -> PoissonResult:
    """Notebook Poisson projection for varphi and normal current."""
    js = np.asarray(edge_js_us_A_m2, dtype=float)
    if js.shape != (ops.n_edges,):
        raise ValueError(f"edge_js_us_A_m2 must have shape ({ops.n_edges},).")
    if poisson_op is None:
        poisson_op = build_poisson_operator(material=material, ops=ops)
    if boundary_accum_A_m is None:
        boundary = np.zeros(ops.n_nodes, dtype=float)
    else:
        boundary = np.asarray(boundary_accum_A_m, dtype=float)

    # Notebook RHS: b_i += -s_ij js_ij, b_j += +s_ij js_ij, plus
    # b_boundary = - outward_boundary_accumulator.
    b = -edge_flux_accumulator_A_m(js, ops) - boundary

    rhs_aug = np.concatenate([b, [0.0]])
    if poisson_op.solver is not None:
        sol = np.asarray(poisson_op.solver.solve(rhs_aug), dtype=float)
    elif spsolve is not None:  # pragma: no cover
        sol = np.asarray(spsolve(poisson_op.A_aug, rhs_aug), dtype=float)
    else:  # pragma: no cover
        sol = np.asarray(np.linalg.solve(poisson_op.A_aug, rhs_aug), dtype=float)

    phi = np.asarray(sol[:-1], dtype=float)
    phi -= float(np.mean(phi))
    lam = float(sol[-1])

    edge_jn = -material.sigma_n_S_m * edge_scalar_gradient(phi, ops)
    edge_jtot = js + edge_jn
    div_js = divergence_from_edge_scalar(js, ops, boundary_accum_A_m=boundary)
    div_jn = divergence_from_edge_scalar(edge_jn, ops)
    div_j = divergence_from_edge_scalar(edge_jtot, ops, boundary_accum_A_m=boundary)

    return PoissonResult(
        phi_V=phi,
        edge_jn_A_m2=edge_jn,
        edge_jtot_A_m2=edge_jtot,
        node_div_jtot_A_m3=div_j,
        node_div_js_A_m3=div_js,
        node_div_jn_A_m3=div_jn,
        lambda_mean=lam,
    )


def solve_poisson_potential(
    *,
    edge_js_us_A_m2: np.ndarray,
    material: GTDGLMaterial,
    ops: FVOperators,
    boundary_accum_A_m: np.ndarray | None = None,
) -> np.ndarray:
    """Backward-compatible wrapper returning only the mean-zero potential.

    The notebook-port solver uses :func:`solve_varphi_poisson`, which returns
    the full Poisson projection bundle. Older tests and scripts imported
    ``solve_poisson_potential`` and expected just ``phi``.
    """
    return solve_varphi_poisson(
        edge_js_us_A_m2=edge_js_us_A_m2,
        material=material,
        ops=ops,
        boundary_accum_A_m=boundary_accum_A_m,
    ).phi_V


def target_terminal_boundary_accum_A_m(
    *,
    edge_data,
    ops: FVOperators,
    material: GTDGLMaterial,
    target_current_A: float,
) -> np.ndarray:
    """Return fixed left/right terminal boundary accumulator."""
    return terminal_boundary_accum_A_m(
        edge_data,
        n_nodes=ops.n_nodes,
        target_current_A=float(target_current_A),
        thickness_m=material.thickness_m,
    )


def kwt_delta_update_attempt(
    *,
    psi_J: np.ndarray,
    phi_V: np.ndarray,
    defs: FormulaFields,
    dt_s: float,
    material: GTDGLMaterial,
) -> tuple[np.ndarray | None, float]:
    """Notebook local semi-implicit KWT update attempt.

    Solves locally
        Delta^{n+1} + z |Delta^{n+1}|^2 = w
    with temporal gauge link U = exp(i 2e varphi dt / hbar).
    """
    tau0 = material.tau0_GL_s
    Delta_n = np.asarray(psi_J, dtype=np.complex128)
    phi = np.asarray(phi_V, dtype=float)
    amp2_n = np.abs(Delta_n) ** 2

    U = np.exp(1j * (2.0 * E_CHARGE_C / HBAR_J_S) * phi * float(dt_s))
    Uinv = np.conjugate(U)

    alpha = defs.alpha_kwt_J_inv2
    z = alpha * Uinv * Delta_n
    w = Uinv * (
        Delta_n
        + alpha * Delta_n * amp2_n
        + (float(dt_s) / tau0) * defs.rho * defs.forcing_J
    )

    ccoef = np.real(w * np.conjugate(z))
    absz2 = np.abs(z) ** 2
    absw2 = np.abs(w) ** 2
    B = 1.0 + 2.0 * ccoef
    discr = B**2 - 4.0 * absz2 * absw2
    discr_min = float(np.nanmin(discr))
    if discr_min < -1.0e-14:
        return None, discr_min

    discr = np.maximum(discr, 0.0)
    denom = B + np.sqrt(discr)
    denom = np.where(np.abs(denom) < 1.0e-300, np.inf, denom)
    amp2_new = 2.0 * absw2 / denom
    amp2_new = np.maximum(amp2_new, 0.0)
    Delta_new = w - z * amp2_new

    if not (
        np.all(np.isfinite(np.real(Delta_new)))
        and np.all(np.isfinite(np.imag(Delta_new)))
    ):
        return None, discr_min
    return Delta_new, discr_min


def kwt_local_update(
    *,
    psi_J: np.ndarray,
    phi_V: np.ndarray,
    Te_K: np.ndarray,
    forcing_J: np.ndarray,
    dt_s: float,
    material: GTDGLMaterial,
    max_phase_step_rad: float = 0.25,
    use_phi_phase: bool = False,
) -> tuple[np.ndarray, bool, float]:
    """Backward-compatible local KWT update used by legacy tests.

    The production OE7 notebook-port path uses ``kwt_delta_update_attempt``
    because it needs the full notebook formula-field bundle. This wrapper keeps
    the old public API available for tests and scripts that only pass an
    externally computed forcing. With zero forcing and ``use_phi_phase=False``
    it preserves the input state exactly, matching the old smoke-test contract.
    """
    psi = np.asarray(psi_J, dtype=np.complex128)
    phi = np.asarray(phi_V, dtype=float)
    Te = np.asarray(Te_K, dtype=float)
    forcing = np.asarray(forcing_J, dtype=np.complex128)

    if dt_s <= 0.0:
        raise ValueError("dt_s must be positive.")
    if psi.shape != forcing.shape:
        raise ValueError("psi_J and forcing_J must have the same shape.")

    rho = np.maximum(material.rho_kwt(Te, np.abs(psi)), 1.0e-30)
    psi_new = psi + (float(dt_s) / material.tau0_GL_s) * forcing / rho

    max_abs_phase = 0.0
    if use_phi_phase:
        phase_step = (2.0 * E_CHARGE_C / HBAR_J_S) * phi * float(dt_s)
        max_abs_phase = float(np.max(np.abs(phase_step))) if phase_step.size else 0.0
        if max_abs_phase > max_phase_step_rad:
            return psi.copy(), False, max_abs_phase
        psi_new = np.exp(1j * phase_step) * psi_new

    if not (
        np.all(np.isfinite(np.real(psi_new)))
        and np.all(np.isfinite(np.imag(psi_new)))
    ):
        return psi.copy(), False, max_abs_phase

    return psi_new, True, max_abs_phase


def apply_stationary_boundary_conditions(
    *,
    psi_trial_J: np.ndarray,
    mesh,
    seed,
    q_bias_m_inv: float,
    material: GTDGLMaterial,
    enabled: bool = True,
) -> np.ndarray:
    """Notebook current-Neumann/vacuum boundary constraints for Delta."""
    psi = np.asarray(psi_trial_J, dtype=np.complex128)
    out = np.nan_to_num(psi, nan=0.0, posinf=0.0, neginf=0.0).astype(np.complex128, copy=True)
    out = clip_gap_amplitude(out, material)
    if not enabled:
        return out

    # Vacuum Neumann top/bottom first.
    for side in ("bottom", "top"):
        boundary, inner = nearest_inward_boundary_pairs(mesh, side)
        out[boundary] = out[inner]

    nodes = np.asarray(mesh.nodes, dtype=float)
    R_bc = seed_delta_bias_J(seed, fallback=float(np.nanmedian(np.abs(out))))
    pairs = terminal_inner_node_pairs(mesh)

    left_boundary, left_inner = pairs["left"]
    dx_left = np.abs(nodes[left_inner, 0] - nodes[left_boundary, 0])
    theta_left = np.angle(out[left_inner]) - float(q_bias_m_inv) * dx_left
    out[left_boundary] = R_bc * np.exp(1j * theta_left)

    right_boundary, right_inner = pairs["right"]
    dx_right = np.abs(nodes[right_boundary, 0] - nodes[right_inner, 0])
    theta_right = np.angle(out[right_inner]) + float(q_bias_m_inv) * dx_right
    out[right_boundary] = R_bc * np.exp(1j * theta_right)

    return clip_gap_amplitude(out, material)


def clip_gap_amplitude(psi_J: np.ndarray, material: GTDGLMaterial) -> np.ndarray:
    """Clip only nonphysical overshoots; do not enforce a positive floor."""
    psi = np.asarray(psi_J, dtype=np.complex128)
    amp = np.abs(psi)
    phase = np.exp(1j * np.angle(psi))
    amp = np.clip(amp, 0.0, 2.0 * material.delta0_J)
    return amp * phase


def boundary_node_masks(mesh) -> dict[str, np.ndarray]:
    """Return rectangular-boundary masks for left/right/bottom/top."""
    nodes = np.asarray(mesh.nodes, dtype=float)
    x = nodes[:, 0]
    y = nodes[:, 1]
    length_m = float(getattr(mesh, "length_m", np.ptp(x)))
    width_m = float(getattr(mesh, "width_m", np.ptp(y)))
    tol = max(1.0e-15, 1.0e-9 * max(length_m, width_m))
    xmin = float(np.min(x))
    xmax = float(np.max(x))
    ymin = float(np.min(y))
    ymax = float(np.max(y))
    return {
        "left": np.abs(x - xmin) <= tol,
        "right": np.abs(x - xmax) <= tol,
        "bottom": np.abs(y - ymin) <= tol,
        "top": np.abs(y - ymax) <= tol,
    }


def terminal_inner_node_pairs(mesh) -> dict[str, tuple[np.ndarray, np.ndarray]]:
    """Pair terminal boundary nodes with their nearest inward neighbours."""
    nodes = np.asarray(mesh.nodes, dtype=float)
    x = nodes[:, 0]
    y = nodes[:, 1]
    masks = boundary_node_masks(mesh)
    h = float(getattr(mesh, "target_spacing_m", 1.0e-9))
    tol = max(1.0e-15, 1.0e-9 * float(getattr(mesh, "length_m", np.ptp(x))))

    def pair(side: str) -> tuple[np.ndarray, np.ndarray]:
        boundary = np.where(masks[side])[0]
        inner = np.empty_like(boundary, dtype=np.int64)
        for k, b in enumerate(boundary):
            if side == "left":
                candidates = np.where(x > x[b] + tol)[0]
                dx = x[candidates] - x[b]
            elif side == "right":
                candidates = np.where(x < x[b] - tol)[0]
                dx = x[b] - x[candidates]
            else:
                raise ValueError(f"Invalid terminal side {side!r}.")
            if candidates.size == 0:
                raise ValueError(f"No inward candidates found for {side} node {b}.")
            dy = y[candidates] - y[b]
            score = (dy / max(h, 1.0e-300)) ** 2 + 0.05 * (dx / max(h, 1.0e-300)) ** 2
            inner[k] = int(candidates[int(np.argmin(score))])
        return boundary.astype(np.int64), inner

    return {"left": pair("left"), "right": pair("right")}


def nearest_inward_boundary_pairs(mesh, side: str) -> tuple[np.ndarray, np.ndarray]:
    """Pair top/bottom boundary nodes with their nearest inward neighbours."""
    if side not in {"bottom", "top"}:
        raise ValueError(f"side must be 'bottom' or 'top', got {side!r}.")
    nodes = np.asarray(mesh.nodes, dtype=float)
    x = nodes[:, 0]
    y = nodes[:, 1]
    masks = boundary_node_masks(mesh)
    h = float(getattr(mesh, "target_spacing_m", 1.0e-9))
    tol = max(
        1.0e-15,
        1.0e-9
        * max(
            float(getattr(mesh, "length_m", np.ptp(x))),
            float(getattr(mesh, "width_m", np.ptp(y))),
        ),
    )
    boundary = np.where(masks[side])[0]
    inner = np.empty_like(boundary, dtype=np.int64)
    for k, b in enumerate(boundary):
        if side == "bottom":
            candidates = np.where(y > y[b] + tol)[0]
            dy = y[candidates] - y[b]
        else:
            candidates = np.where(y < y[b] - tol)[0]
            dy = y[b] - y[candidates]
        if candidates.size == 0:
            raise ValueError(f"No inward candidates found for {side} node {b}.")
        dx = x[candidates] - x[b]
        score = (dx / max(h, 1.0e-300)) ** 2 + 0.05 * (dy / max(h, 1.0e-300)) ** 2
        inner[k] = int(candidates[int(np.argmin(score))])
    return boundary.astype(np.int64), inner


def current_residual(
    currents: CurrentFields,
    mesh,
    material: GTDGLMaterial | None = None,
    target_current_A: float | None = None,
) -> float:
    """Dimensionless RMS residual of div(j_tot).

    New notebook-order runs pass ``material`` and ``target_current_A`` so the
    scale is the imposed average current density divided by the mesh spacing.
    Older smoke tests called this with only ``(currents, mesh)``; in that case
    we fall back to the reconstructed total-current RMS scale.
    """
    div = np.asarray(currents.node_div_jtot_A_m3, dtype=float)
    h = float(getattr(mesh, "target_spacing_m", getattr(mesh, "xi_mesh_m", 1.0e-9)))
    if material is not None and target_current_A is not None:
        jscale = abs(target_current_density_A_m2(material, float(target_current_A)))
    else:
        jscale = float(
            np.sqrt(
                np.nanmean(
                    currents.node_jtot_x_A_m2**2 + currents.node_jtot_y_A_m2**2
                )
            )
        )
    scale = max(jscale / max(h, 1.0e-300), 1.0)
    return float(np.sqrt(np.nanmean(div * div)) / scale)


def max_current_residual(
    currents: CurrentFields,
    mesh,
    material: GTDGLMaterial | None = None,
    target_current_A: float | None = None,
) -> float:
    div = np.asarray(currents.node_div_jtot_A_m3, dtype=float)
    h = float(getattr(mesh, "target_spacing_m", getattr(mesh, "xi_mesh_m", 1.0e-9)))
    if material is not None and target_current_A is not None:
        jscale = abs(target_current_density_A_m2(material, float(target_current_A)))
    else:
        jscale = float(
            np.sqrt(
                np.nanmean(
                    currents.node_jtot_x_A_m2**2 + currents.node_jtot_y_A_m2**2
                )
            )
        )
    scale = max(jscale / max(h, 1.0e-300), 1.0)
    return float(np.nanmax(np.abs(div)) / scale)


def normal_current_fraction_rms(currents: CurrentFields) -> float:
    """RMS normal-current fraction relative to total-current edge scale."""
    num = float(np.sqrt(np.nanmean(currents.edge_jn_A_m2**2))) if currents.edge_jn_A_m2.size else 0.0
    den = float(np.sqrt(np.nanmean(currents.edge_jtot_A_m2**2))) if currents.edge_jtot_A_m2.size else 0.0
    return num / max(den, 1.0e-300)


def current_density_maxima_A_m2(currents: CurrentFields) -> tuple[float, float]:
    """Return max |j_n| and max |j_tot| from edge fields."""
    jn_max = float(np.nanmax(np.abs(currents.edge_jn_A_m2))) if currents.edge_jn_A_m2.size else 0.0
    jt_max = float(np.nanmax(np.abs(currents.edge_jtot_A_m2))) if currents.edge_jtot_A_m2.size else 0.0
    return jn_max, jt_max


def normal_current_fraction_max(currents: CurrentFields) -> float:
    jn_max, jt_max = current_density_maxima_A_m2(currents)
    return jn_max / max(jt_max, 1.0e-300)


def seed_target_current_A(seed) -> float:
    """Extract imposed transport current from an OE6 seed-like object."""
    for name in ("I_bias_A", "target_current_A", "current_A"):
        if hasattr(seed, name):
            value = float(getattr(seed, name))
            if np.isfinite(value):
                return value
    metadata = getattr(seed, "metadata", None)
    if isinstance(metadata, dict):
        for name in ("I_bias_A", "target_current_A", "current_A"):
            if name in metadata:
                value = float(metadata[name])
                if np.isfinite(value):
                    return value
    return 0.0


def seed_q_bias_m_inv(seed, *, target_current_A: float | None = None) -> float:
    """Extract seed phase-gradient q."""
    for name in ("q_bias_m_inv", "target_q_m_inv", "q_m_inv"):
        if hasattr(seed, name):
            value = float(getattr(seed, name))
            if np.isfinite(value):
                return value
    metadata = getattr(seed, "metadata", None)
    if isinstance(metadata, dict):
        for name in ("q_bias_m_inv", "target_q_m_inv", "q_m_inv"):
            if name in metadata:
                value = float(metadata[name])
                if np.isfinite(value):
                    return value
    if target_current_A is not None and abs(float(target_current_A)) <= 0.0:
        return 0.0
    return 0.0


def seed_delta_bias_J(seed, *, fallback: float) -> float:
    """Extract stationary terminal amplitude from the OE6 seed if available."""
    for name in ("delta_bias_J", "target_delta_J", "Delta_bias_J"):
        if hasattr(seed, name):
            value = float(getattr(seed, name))
            if np.isfinite(value) and value >= 0.0:
                return value
    metadata = getattr(seed, "metadata", None)
    if isinstance(metadata, dict):
        for name in ("delta_bias_J", "target_delta_J", "Delta_bias_J"):
            if name in metadata:
                value = float(metadata[name])
                if np.isfinite(value) and value >= 0.0:
                    return value
    if hasattr(seed, "node_delta_J"):
        arr = np.asarray(seed.node_delta_J, dtype=float)
        finite = arr[np.isfinite(arr)]
        if finite.size:
            return float(np.nanmedian(finite))
    return float(fallback)


def target_current_density_A_m2(material: GTDGLMaterial, target_current_A: float) -> float:
    return float(target_current_A) / max(material.width_m * material.thickness_m, 1.0e-300)


def suggest_next_dt(
    *,
    dt_s: float,
    max_amp2_change_rel: float,
    retries: int,
    adaptive: bool,
    target: float,
    shrink_factor: float,
    grow_factor: float,
    dt_min_s: float,
    dt_max_s: float,
) -> float:
    """Notebook adaptive dt rule."""
    if not adaptive:
        return float(dt_s)
    if retries > 0 or max_amp2_change_rel > 2.0 * target:
        return max(float(dt_s) * float(shrink_factor), float(dt_min_s))
    if max_amp2_change_rel < 0.25 * target:
        return min(float(dt_s) * float(grow_factor), float(dt_max_s))
    return min(float(dt_s), float(dt_max_s))


def relax_stationary_gtdgl(
    *,
    mesh,
    edge_data,
    seed,
    material: GTDGLMaterial,
    ops: FVOperators,
    steps: int = 2000,
    dt_s: float = 1.0e-17,
    min_steps: int = 10,
    tolerance_eta: float = 1.0e-9,
    tolerance_current_residual: float = 1.0e-6,
    eta_reject: float = 5.0e-4,
    adapt_dt: bool = True,
    dt_min_s: float = 1.0e-22,
    dt_max_s: float = 1.0e-13,
    lock_terminals: bool = True,
    target_current_A: float | None = None,
    progress: bool = False,
    n_phi_snapshots: int = 6,
    use_phi_phase: bool = True,
) -> RelaxationResult:
    """Relax the OE6 seed with frozen temperatures and notebook solver ordering."""
    del use_phi_phase  # Notebook always uses the temporal gauge link in KWT.

    if steps <= 0:
        raise ValueError("steps must be positive.")
    if dt_s <= 0.0:
        raise ValueError("dt_s must be positive.")
    if min_steps < 0:
        raise ValueError("min_steps must be non-negative.")

    if target_current_A is None:
        target_current_A = seed_target_current_A(seed)
    target_current_A = float(target_current_A)
    q_bias = seed_q_bias_m_inv(seed, target_current_A=target_current_A)
    javg = target_current_density_A_m2(material, target_current_A)
    q_ref = abs(q_bias) if abs(q_bias) > 0.0 else 1.0

    psi0 = (
        np.asarray(seed.node_psi_real_J, dtype=float)
        + 1j * np.asarray(seed.node_psi_imag_J, dtype=float)
    )
    phi0 = np.asarray(seed.node_phi_electric_V, dtype=float).copy()
    Te = np.asarray(seed.node_Te_K, dtype=float).copy()
    Tph = np.asarray(seed.node_Tph_K, dtype=float).copy()

    psi = apply_stationary_boundary_conditions(
        psi_trial_J=psi0,
        mesh=mesh,
        seed=seed,
        q_bias_m_inv=q_bias,
        material=material,
        enabled=lock_terminals,
    )
    phi = phi0 - float(np.mean(phi0))

    boundary_accum = target_terminal_boundary_accum_A_m(
        edge_data=edge_data,
        ops=ops,
        material=material,
        target_current_A=target_current_A,
    )
    poisson_op = build_poisson_operator(material=material, ops=ops)

    # Notebook initial projection: compute defs, then Poisson, then recompute fields.
    defs0 = compute_formula_fields(psi_J=psi, Te_K=Te, material=material, ops=ops)
    poisson0 = solve_varphi_poisson(
        edge_js_us_A_m2=defs0.edge_js_us_A_m2,
        material=material,
        ops=ops,
        poisson_op=poisson_op,
        boundary_accum_A_m=boundary_accum,
    )
    phi = poisson0.phi_V
    currents = compute_current_fields(
        psi_J=psi,
        phi_V=phi,
        Te_K=Te,
        material=material,
        ops=ops,
        boundary_accum_A_m=boundary_accum,
    )

    t_s = 0.0
    accepted = 0
    rejected = 0
    converged = False

    hist_keys = [
        "t_s",
        "dt_s",
        "retries",
        "discr_min",
        "eta_R",
        "max_amp2_change_rel",
        "current_residual",
        "current_residual_max",
        "terminal_voltage_V",
        "pairbreaking_max",
        "delta_min_over_delta0",
        "delta_max_over_delta0",
        "normal_current_fraction_rms",
        "normal_current_fraction_max",
        "normal_current_max_A_m2",
        "total_current_max_A_m2",
        "median_Q_m_inv",
        "p95_Q_m_inv",
        "max_Q_m_inv",
        "max_js_A_m2",
        "max_j_A_m2",
    ]
    hist: dict[str, list[float]] = {key: [] for key in hist_keys}

    n_phi_snapshots = max(2, int(n_phi_snapshots))
    snapshot_steps = set(np.linspace(0, int(steps) - 1, n_phi_snapshots, dtype=int).tolist())
    snapshots: dict[str, list[np.ndarray] | list[float]] = {
        "snapshot_t_s": [],
        "psi_snapshot_real_J": [],
        "psi_snapshot_imag_J": [],
        "delta_snapshot_meV": [],
        "phi_snapshot_V": [],
        "current_density_snapshot_A_m2": [],
        "current_density_snapshot_x_A_m2": [],
        "current_density_snapshot_y_A_m2": [],
        "supercurrent_density_snapshot_A_m2": [],
        "supercurrent_density_snapshot_x_A_m2": [],
        "supercurrent_density_snapshot_y_A_m2": [],
        "normal_current_density_snapshot_A_m2": [],
        "normal_current_density_snapshot_x_A_m2": [],
        "normal_current_density_snapshot_y_A_m2": [],
        "divergence_snapshot_A_m3": [],
        "pairbreaking_ratio_snapshot": [],
    }

    def append_snapshot() -> None:
        jtot_mag = np.sqrt(currents.node_jtot_x_A_m2**2 + currents.node_jtot_y_A_m2**2)
        js_mag = np.sqrt(currents.node_js_us_x_A_m2**2 + currents.node_js_us_y_A_m2**2)
        jn_mag = np.sqrt(currents.node_jn_x_A_m2**2 + currents.node_jn_y_A_m2**2)
        snapshots["snapshot_t_s"].append(float(t_s))
        snapshots["psi_snapshot_real_J"].append(np.real(psi).copy())
        snapshots["psi_snapshot_imag_J"].append(np.imag(psi).copy())
        snapshots["delta_snapshot_meV"].append(np.abs(psi).copy() / MEV_J)
        snapshots["phi_snapshot_V"].append(phi.copy())
        snapshots["current_density_snapshot_A_m2"].append(jtot_mag.copy())
        snapshots["current_density_snapshot_x_A_m2"].append(currents.node_jtot_x_A_m2.copy())
        snapshots["current_density_snapshot_y_A_m2"].append(currents.node_jtot_y_A_m2.copy())
        snapshots["supercurrent_density_snapshot_A_m2"].append(js_mag.copy())
        snapshots["supercurrent_density_snapshot_x_A_m2"].append(currents.node_js_us_x_A_m2.copy())
        snapshots["supercurrent_density_snapshot_y_A_m2"].append(currents.node_js_us_y_A_m2.copy())
        snapshots["normal_current_density_snapshot_A_m2"].append(jn_mag.copy())
        snapshots["normal_current_density_snapshot_x_A_m2"].append(currents.node_jn_x_A_m2.copy())
        snapshots["normal_current_density_snapshot_y_A_m2"].append(currents.node_jn_y_A_m2.copy())
        snapshots["divergence_snapshot_A_m3"].append(currents.node_div_jtot_A_m3.copy())
        snapshots["pairbreaking_ratio_snapshot"].append(currents.node_pairbreaking_ratio.copy())

    append_snapshot()

    iterator = range(int(steps))
    if progress and trange is not None:
        iterator = trange(int(steps), desc="OE7 notebook KWT", leave=True)

    for n in iterator:
        retries = 0
        dt_eff = float(dt_s)
        while True:
            defs_n = compute_formula_fields(psi_J=psi, Te_K=Te, material=material, ops=ops)
            psi_new, discr_min = kwt_delta_update_attempt(
                psi_J=psi,
                phi_V=phi,
                defs=defs_n,
                dt_s=dt_eff,
                material=material,
            )
            if psi_new is not None:
                break
            retries += 1
            rejected += 1
            if retries > 30:
                raise RuntimeError(
                    "Failed KWT update: negative discriminant after "
                    f"{retries} retries. Last min={discr_min:.3e}"
                )
            dt_eff = max(dt_min_s, 0.5 * dt_eff)

        psi_trial = apply_stationary_boundary_conditions(
            psi_trial_J=psi_new,
            mesh=mesh,
            seed=seed,
            q_bias_m_inv=q_bias,
            material=material,
            enabled=lock_terminals,
        )
        defs_pre = compute_formula_fields(psi_J=psi_trial, Te_K=Te, material=material, ops=ops)
        poisson = solve_varphi_poisson(
            edge_js_us_A_m2=defs_pre.edge_js_us_A_m2,
            material=material,
            ops=ops,
            poisson_op=poisson_op,
            boundary_accum_A_m=boundary_accum,
        )
        phi_trial = poisson.phi_V
        trial_currents = compute_current_fields(
            psi_J=psi_trial,
            phi_V=phi_trial,
            Te_K=Te,
            material=material,
            ops=ops,
            boundary_accum_A_m=boundary_accum,
        )

        amp2_change_rel = float(
            np.nanmax(np.abs(np.abs(psi_trial) ** 2 - np.abs(psi) ** 2))
            / material.delta0_J**2
        )

        psi = psi_trial
        phi = phi_trial
        currents = trial_currents
        t_s += dt_eff
        accepted += 1

        residual = current_residual(currents, mesh, material, target_current_A)
        residual_max = max_current_residual(currents, mesh, material, target_current_A)
        voltage = terminal_voltage(np.asarray(mesh.nodes, dtype=float), phi, length_m=float(mesh.length_m))
        pb_max = float(np.nanmax(currents.node_pairbreaking_ratio))
        delta_min_ratio = float(np.nanmin(np.abs(psi)) / material.delta0_J)
        delta_max_ratio = float(np.nanmax(np.abs(psi)) / material.delta0_J)
        normal_frac = normal_current_fraction_rms(currents)
        normal_max_frac = normal_current_fraction_max(currents)
        jn_max_A_m2, jt_max_A_m2 = current_density_maxima_A_m2(currents)
        Qabs = np.abs(currents.edge_Q_m_inv)

        values = {
            "t_s": t_s,
            "dt_s": dt_eff,
            "retries": float(retries),
            "discr_min": float(discr_min),
            "eta_R": amp2_change_rel,
            "max_amp2_change_rel": amp2_change_rel,
            "current_residual": residual,
            "current_residual_max": residual_max,
            "terminal_voltage_V": voltage,
            "pairbreaking_max": pb_max,
            "delta_min_over_delta0": delta_min_ratio,
            "delta_max_over_delta0": delta_max_ratio,
            "normal_current_fraction_rms": normal_frac,
            "normal_current_fraction_max": normal_max_frac,
            "normal_current_max_A_m2": jn_max_A_m2,
            "total_current_max_A_m2": jt_max_A_m2,
            "median_Q_m_inv": float(np.nanmedian(Qabs)),
            "p95_Q_m_inv": float(np.nanpercentile(Qabs, 95.0)),
            "max_Q_m_inv": float(np.nanmax(Qabs)),
            "max_js_A_m2": float(np.nanmax(np.abs(currents.edge_js_us_A_m2))),
            "max_j_A_m2": float(np.nanmax(np.abs(currents.edge_jtot_A_m2))),
        }
        for key in hist_keys:
            hist[key].append(values[key])

        if n in snapshot_steps:
            append_snapshot()

        if progress and hasattr(iterator, "set_postfix") and accepted % 10 == 0:
            iterator.set_postfix(
                dA2=f"{amp2_change_rel:.2e}",
                eps=f"{residual:.2e}",
                V=f"{voltage:.2e}",
                chi=f"{pb_max:.3g}",
                dt_fs=f"{dt_eff / 1.0e-15:.3g}",
            )

        if accepted >= min_steps and amp2_change_rel < tolerance_eta and residual < tolerance_current_residual:
            converged = True
            break

        dt_s = suggest_next_dt(
            dt_s=dt_eff,
            max_amp2_change_rel=amp2_change_rel,
            retries=retries,
            adaptive=adapt_dt,
            target=float(eta_reject),
            shrink_factor=0.7,
            grow_factor=1.05,
            dt_min_s=dt_min_s,
            dt_max_s=dt_max_s,
        )

        if not np.all(np.isfinite(psi)) or not np.all(np.isfinite(phi)):
            raise FloatingPointError(f"Stopped: non-finite state at accepted step {accepted}.")

    if len(snapshots["snapshot_t_s"]) == 0 or snapshots["snapshot_t_s"][-1] != t_s:
        append_snapshot()

    # Keep exactly n_phi_snapshots, preserving first and final snapshots.
    n_snap = len(snapshots["snapshot_t_s"])
    if n_snap > n_phi_snapshots:
        keep = np.unique(np.rint(np.linspace(0, n_snap - 1, n_phi_snapshots)).astype(int))
        if keep[-1] != n_snap - 1:
            keep[-1] = n_snap - 1
        for key, seq in list(snapshots.items()):
            snapshots[key] = [seq[int(i)] for i in keep]

    boundary = boundary_currents_from_node_vectors(
        mesh=mesh,
        edge_data=edge_data,
        jx_A_m2=currents.node_jtot_x_A_m2,
        jy_A_m2=currents.node_jtot_y_A_m2,
        thickness_m=material.thickness_m,
    )
    voltage = terminal_voltage(np.asarray(mesh.nodes, dtype=float), phi, length_m=float(mesh.length_m))
    normal_ohmic_voltage = (
        float(target_current_A)
        * float(mesh.length_m)
        / (material.sigma_n_S_m * material.width_m * material.thickness_m)
    )
    normal_max_A_m2, total_max_A_m2 = current_density_maxima_A_m2(currents)

    summary = {
        "backend": "oe7_notebook_order_kwt_poisson_v1",
        "gauge_policy": "notebook_temporal_gauge_link_in_kwt",
        "converged": bool(converged),
        "accepted_steps": int(accepted),
        "rejected_steps": int(rejected),
        "final_time_ps": float(t_s / 1.0e-12),
        "tau_scale": float(material.tau_scale),
        "tau_ee_Tc_effective_ps": float(material.tau_scale * material.tau_ee_Tc_s / 1.0e-12),
        "tau_ep_Tc_effective_ps": float(material.tau_scale * material.tau_ep_Tc_s / 1.0e-12),
        "target_current_A": float(target_current_A),
        "target_q_bias_m_inv": float(q_bias),
        "target_j_bias_A_m2": float(javg),
        "terminal_voltage_V": float(voltage),
        "normal_ohmic_voltage_V": float(normal_ohmic_voltage),
        "terminal_voltage_over_normal": float(
            voltage / normal_ohmic_voltage if normal_ohmic_voltage != 0.0 else float("nan")
        ),
        "normal_current_fraction_rms": float(normal_current_fraction_rms(currents)),
        "normal_current_fraction_max": float(normal_max_A_m2 / max(total_max_A_m2, 1.0e-300)),
        "normal_current_max_A_m2": float(normal_max_A_m2),
        "total_current_max_A_m2": float(total_max_A_m2),
        "current_residual": float(current_residual(currents, mesh, material, target_current_A)),
        "eta_R_final": float(hist["eta_R"][-1]) if hist["eta_R"] else float("nan"),
        "divergence_rms_A_m3": float(np.sqrt(np.nanmean(currents.node_div_jtot_A_m3**2))),
        "min_delta_over_delta0": float(np.nanmin(np.abs(psi)) / material.delta0_J),
        "mean_delta_over_delta0": float(np.nanmean(np.abs(psi)) / material.delta0_J),
        "max_pairbreaking_ratio": float(np.nanmax(currents.node_pairbreaking_ratio)),
        "p99_pairbreaking_ratio": float(np.nanpercentile(currents.node_pairbreaking_ratio, 99.0)),
        "edge_Q_max_m_inv": float(np.nanmax(np.abs(currents.edge_Q_m_inv))),
        "boundary_currents_A": boundary,
    }

    metadata = {
        "backend": summary["backend"],
        "description": "Notebook-order frozen-temperature gTDGL/Poisson relaxation.",
        "thermal_policy": "frozen_Te_Tph",
        "circuit_policy": "inactive",
        "boundary_policy": "current_neumann_from_seed_q_and_seed_delta",
        "poisson_policy": "notebook_conservative_FV_mean_zero_gauge",
        "pairbreaking_ratio": "xi^2 Q^2 / (1 - T/Tc)",
    }

    state = GTDGLStationaryState(
        psi_J=psi,
        phi_V=phi,
        Te_K=Te,
        Tph_K=Tph,
        currents=currents,
        metadata=metadata,
    )

    history: dict[str, np.ndarray] = {key: np.asarray(val, dtype=float) for key, val in hist.items()}
    snapshot_t_s = np.asarray(snapshots["snapshot_t_s"], dtype=float)
    history.update(
        {
            "delta0_meV": np.asarray([material.delta0_J / MEV_J], dtype=float),
            "javg_A_m2": np.asarray([javg], dtype=float),
            "qref_m_inv": np.asarray([q_ref], dtype=float),
            "snapshot_t_s": snapshot_t_s,
            "phi_snapshot_t_s": snapshot_t_s,
            "phi_snapshot_V": np.asarray(snapshots["phi_snapshot_V"], dtype=float),
            "psi_snapshot_t_s": snapshot_t_s,
            "psi_snapshot_real_J": np.asarray(snapshots["psi_snapshot_real_J"], dtype=float),
            "psi_snapshot_imag_J": np.asarray(snapshots["psi_snapshot_imag_J"], dtype=float),
            "delta_snapshot_t_s": snapshot_t_s,
            "delta_snapshot_meV": np.asarray(snapshots["delta_snapshot_meV"], dtype=float),
            "current_snapshot_t_s": snapshot_t_s,
            "current_density_snapshot_A_m2": np.asarray(snapshots["current_density_snapshot_A_m2"], dtype=float),
            "current_density_snapshot_x_A_m2": np.asarray(snapshots["current_density_snapshot_x_A_m2"], dtype=float),
            "current_density_snapshot_y_A_m2": np.asarray(snapshots["current_density_snapshot_y_A_m2"], dtype=float),
            "jtot_snapshot_t_s": snapshot_t_s,
            "jtot_snapshot_mag_A_m2": np.asarray(snapshots["current_density_snapshot_A_m2"], dtype=float),
            "jtot_snapshot_x_A_m2": np.asarray(snapshots["current_density_snapshot_x_A_m2"], dtype=float),
            "jtot_snapshot_y_A_m2": np.asarray(snapshots["current_density_snapshot_y_A_m2"], dtype=float),
            "supercurrent_snapshot_t_s": snapshot_t_s,
            "supercurrent_density_snapshot_A_m2": np.asarray(snapshots["supercurrent_density_snapshot_A_m2"], dtype=float),
            "supercurrent_density_snapshot_x_A_m2": np.asarray(snapshots["supercurrent_density_snapshot_x_A_m2"], dtype=float),
            "supercurrent_density_snapshot_y_A_m2": np.asarray(snapshots["supercurrent_density_snapshot_y_A_m2"], dtype=float),
            "js_us_snapshot_t_s": snapshot_t_s,
            "js_us_snapshot_mag_A_m2": np.asarray(snapshots["supercurrent_density_snapshot_A_m2"], dtype=float),
            "js_us_snapshot_x_A_m2": np.asarray(snapshots["supercurrent_density_snapshot_x_A_m2"], dtype=float),
            "js_us_snapshot_y_A_m2": np.asarray(snapshots["supercurrent_density_snapshot_y_A_m2"], dtype=float),
            "normal_current_snapshot_t_s": snapshot_t_s,
            "normal_current_density_snapshot_A_m2": np.asarray(snapshots["normal_current_density_snapshot_A_m2"], dtype=float),
            "normal_current_density_snapshot_x_A_m2": np.asarray(snapshots["normal_current_density_snapshot_x_A_m2"], dtype=float),
            "normal_current_density_snapshot_y_A_m2": np.asarray(snapshots["normal_current_density_snapshot_y_A_m2"], dtype=float),
            "jn_snapshot_t_s": snapshot_t_s,
            "jn_snapshot_mag_A_m2": np.asarray(snapshots["normal_current_density_snapshot_A_m2"], dtype=float),
            "jn_snapshot_x_A_m2": np.asarray(snapshots["normal_current_density_snapshot_x_A_m2"], dtype=float),
            "jn_snapshot_y_A_m2": np.asarray(snapshots["normal_current_density_snapshot_y_A_m2"], dtype=float),
            "divergence_snapshot_t_s": snapshot_t_s,
            "divergence_snapshot_A_m3": np.asarray(snapshots["divergence_snapshot_A_m3"], dtype=float),
            "pairbreaking_snapshot_t_s": snapshot_t_s,
            "pairbreaking_ratio_snapshot": np.asarray(snapshots["pairbreaking_ratio_snapshot"], dtype=float),
        }
    )

    return RelaxationResult(state=state, history=history, summary=summary)


def save_stationary_state_npz(state: GTDGLStationaryState, output_path: str | Path) -> Path:
    """Save a relaxed stationary state to NPZ."""
    output = Path(output_path)
    output.parent.mkdir(parents=True, exist_ok=True)
    c = state.currents
    np.savez_compressed(
        output,
        psi_real_J=np.real(state.psi_J),
        psi_imag_J=np.imag(state.psi_J),
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
        edge_pairbreaking_ratio=c.edge_pairbreaking_ratio,
        node_pairbreaking_ratio=c.node_pairbreaking_ratio,
        metadata_json=json.dumps(state.metadata, sort_keys=True),
    )
    return output


def save_relaxation_history_npz(history: dict[str, np.ndarray], output_path: str | Path) -> Path:
    """Save compact relaxation history to NPZ."""
    output = Path(output_path)
    output.parent.mkdir(parents=True, exist_ok=True)
    arrays = {key: np.asarray(value) for key, value in history.items()}
    np.savez_compressed(output, **arrays)
    return output
