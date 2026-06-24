"""Stationary gTDGL/Poisson relaxation for pySNSPD OE7.

OE7 starts from the analytic OE6 stationary seed, keeps Te and Tph frozen,
keeps the external circuit inactive, and relaxes only the mesoscopic
gTDGL/Poisson sector.

Mathematical policy used here
-----------------------------
1. The longitudinal terminals impose the superconducting phase gradient from
   the OE6 seed,

       theta_boundary = theta_inner +/- q_bias * dx,

   and use a Neumann condition for the gap amplitude,

       |Delta_boundary| = |Delta_inner|.

   This is the closest discrete analogue of the notebook-style node-to-node
   continuity condition.

2. The insulating top/bottom boundaries impose vacuum Neumann conditions,

       Delta_boundary = Delta_inner.

3. Poisson is retained. It solves the conservative FV equation

       div_h(j_s + j_n) + b_h = 0,
       j_n = -sigma_n grad(phi),

   where b_h is only the prescribed terminal flux accumulator. We do not
   reconstruct an additional terminal current from the evolving field and we
   do not impose the transport current twice.

4. The electrostatic potential is solved with a mean-zero gauge constraint.
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
    from scipy.sparse.linalg import spsolve
except Exception:  # pragma: no cover
    coo_matrix = None
    csr_matrix = None
    bmat = None
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
    edge_phase_gradient_from_psi,
    edge_scalar_gradient,
    edge_scalar_to_node_vector,
    laplacian,
    terminal_boundary_accum_A_m,
    terminal_voltage,
)


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
    boundary_accum_A_m: np.ndarray | None = None,
) -> CurrentFields:
    """Evaluate Usadel-like, GL, normal and total current fields.

    Currents live on oriented edges. Node vectors are diagnostics obtained by
    direct edge-to-node averaging, not by a least-squares reconstruction.
    """

    psi = np.asarray(psi_J, dtype=np.complex128)
    phi = np.asarray(phi_V, dtype=float)
    Te = np.asarray(Te_K, dtype=float)

    if psi.shape != (ops.n_nodes,):
        raise ValueError(f"psi_J must have shape ({ops.n_nodes},), got {psi.shape}.")
    if phi.shape != (ops.n_nodes,):
        raise ValueError(f"phi_V must have shape ({ops.n_nodes},), got {phi.shape}.")
    if Te.shape != (ops.n_nodes,):
        raise ValueError(f"Te_K must have shape ({ops.n_nodes},), got {Te.shape}.")

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
    div_total = divergence_from_edge_scalar(
        edge_jtot,
        ops,
        boundary_accum_A_m=boundary_accum_A_m,
    )

    js_x, js_y = edge_scalar_to_node_vector(edge_js_us, ops)
    jn_x, jn_y = edge_scalar_to_node_vector(edge_jn, ops)
    jt_x, jt_y = edge_scalar_to_node_vector(edge_jtot, ops)

    edge_pb = pairbreaking_ratio_edges(
        Q_edge_m_inv=Q_edge,
        Te_edge_K=Te_edge,
        material=material,
    )
    node_pb = edge_to_node_weighted_average(edge_pb, ops)

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
        edge_pairbreaking_ratio=edge_pb,
        node_pairbreaking_ratio=node_pb,
    )


def pairbreaking_ratio_edges(
    *,
    Q_edge_m_inv: np.ndarray,
    Te_edge_K: np.ndarray,
    material: GTDGLMaterial,
) -> np.ndarray:
    """Return the local GL pairbreaking ratio xi^2 Q^2 / (1 - T/Tc).

    The stationary GL amplitude satisfies

        |Delta|^2 / Delta_mod^2 = 1 - T/Tc - xi^2 Q^2.

    Therefore chi_pb = 1 is the local threshold where the GL stationary
    amplitude vanishes.
    """

    Q = np.asarray(Q_edge_m_inv, dtype=float)
    Te = np.asarray(Te_edge_K, dtype=float)

    a = np.maximum(1.0 - Te / material.Tc_K, 1.0e-30)
    xi2 = material.xi_mod_squared_m2(Te)

    return np.asarray(xi2 * Q * Q / a, dtype=float)


def edge_to_node_weighted_average(edge_values: np.ndarray, ops: FVOperators) -> np.ndarray:
    """Average edge quantities to nodes using dual-face weights."""

    values = np.asarray(edge_values, dtype=float)
    if values.shape != (ops.n_edges,):
        raise ValueError(
            f"edge_values must have shape ({ops.n_edges},), got {values.shape}."
        )

    weights = np.maximum(np.asarray(ops.dual_face_length_m, dtype=float), 1.0e-300)
    out = np.zeros(ops.n_nodes, dtype=float)
    wsum = np.zeros(ops.n_nodes, dtype=float)

    np.add.at(out, ops.edge_i, weights * values)
    np.add.at(out, ops.edge_j, weights * values)
    np.add.at(wsum, ops.edge_i, weights)
    np.add.at(wsum, ops.edge_j, weights)

    return out / np.maximum(wsum, 1.0e-300)


def solve_poisson_potential(
    *,
    edge_js_us_A_m2: np.ndarray,
    material: GTDGLMaterial,
    ops: FVOperators,
    boundary_accum_A_m: np.ndarray | None = None,
) -> np.ndarray:
    """Solve the conservative FV Poisson equation for phi.

    We solve

        div_h(j_s - sigma grad phi) + b_h = 0.

    Multiplying by the nodal control volume gives

        A_phi phi = -(accum_s + b),

    where A_phi is the positive graph operator associated with
    -sigma grad(phi). The gauge is fixed by sum_i phi_i = 0.
    """

    js = np.asarray(edge_js_us_A_m2, dtype=float)
    if js.shape != (ops.n_edges,):
        raise ValueError(f"edge_js_us_A_m2 must have shape ({ops.n_edges},).")

    if boundary_accum_A_m is None:
        boundary = np.zeros(ops.n_nodes, dtype=float)
    else:
        boundary = np.asarray(boundary_accum_A_m, dtype=float)
        if boundary.shape != (ops.n_nodes,):
            raise ValueError(
                f"boundary_accum_A_m must have shape ({ops.n_nodes},), "
                f"got {boundary.shape}."
            )

    source = edge_flux_accumulator_A_m(js, ops) + boundary
    source = project_source_to_zero_sum(source, ops)

    g = material.sigma_n_S_m * ops.dual_face_length_m / ops.edge_length_m
    i = np.asarray(ops.edge_i, dtype=np.int64)
    j = np.asarray(ops.edge_j, dtype=np.int64)
    n = int(ops.n_nodes)

    rows = np.concatenate([i, i, j, j])
    cols = np.concatenate([i, j, j, i])
    data = np.concatenate([g, -g, g, -g])

    rhs = -source

    if coo_matrix is not None and csr_matrix is not None and bmat is not None and spsolve is not None:
        A = coo_matrix((data, (rows, cols)), shape=(n, n)).tocsr()
        ones_col = csr_matrix(np.ones((n, 1), dtype=float))
        ones_row = csr_matrix(np.ones((1, n), dtype=float))
        zero = csr_matrix((1, 1), dtype=float)

        A_aug = bmat(
            [
                [A, ones_col],
                [ones_row, zero],
            ],
            format="csr",
        )
        rhs_aug = np.concatenate([rhs, [0.0]])
        sol = np.asarray(spsolve(A_aug, rhs_aug), dtype=float)
        phi = sol[:n]
    else:  # pragma: no cover
        A = np.zeros((n, n), dtype=float)
        for a, b, gg in zip(i, j, g):
            A[a, a] += gg
            A[a, b] -= gg
            A[b, b] += gg
            A[b, a] -= gg

        A_aug = np.zeros((n + 1, n + 1), dtype=float)
        A_aug[:n, :n] = A
        A_aug[:n, n] = 1.0
        A_aug[n, :n] = 1.0

        rhs_aug = np.concatenate([rhs, [0.0]])
        sol = np.linalg.solve(A_aug, rhs_aug)
        phi = np.asarray(sol[:n], dtype=float)

    phi -= float(np.mean(phi))
    return phi


def edge_flux_accumulator_A_m(
    edge_current_i_to_j_A_m2: np.ndarray,
    ops: FVOperators,
) -> np.ndarray:
    """Return conservative edge-flux accumulator before division by node area."""

    current = np.asarray(edge_current_i_to_j_A_m2, dtype=float)
    if current.shape != (ops.n_edges,):
        raise ValueError(f"edge current must have shape ({ops.n_edges},).")

    flux = ops.dual_face_length_m * current
    out = np.zeros(ops.n_nodes, dtype=float)

    np.add.at(out, ops.edge_i, flux)
    np.add.at(out, ops.edge_j, -flux)

    return out


def project_source_to_zero_sum(source_A_m: np.ndarray, ops: FVOperators) -> np.ndarray:
    """Project a Neumann source onto the range of the graph Laplacian."""

    out = np.asarray(source_A_m, dtype=float).copy()
    total = float(np.sum(out))
    norm = float(np.sum(np.abs(out)))

    if not np.isfinite(total):
        raise FloatingPointError("Poisson source has non-finite total sum.")

    if norm == 0.0 or abs(total) <= 1.0e-14 * max(norm, 1.0):
        return out

    weights = np.asarray(ops.node_area_m2, dtype=float)
    weights = np.maximum(weights, 1.0e-300)

    out -= total * weights / float(np.sum(weights))
    return out


def target_terminal_boundary_accum_A_m(
    *,
    edge_data,
    ops: FVOperators,
    material: GTDGLMaterial,
    target_current_A: float,
) -> np.ndarray:
    """Return the fixed left/right terminal boundary accumulator.

    Positive current flows from left to right. The outward terminal flux is
    negative at the left contact and positive at the right contact.
    """

    return terminal_boundary_accum_A_m(
        edge_data,
        n_nodes=ops.n_nodes,
        target_current_A=float(target_current_A),
        thickness_m=material.thickness_m,
    )


def gtdgl_forcing(
    *,
    psi_J: np.ndarray,
    Te_K: np.ndarray,
    currents: CurrentFields,
    material: GTDGLMaterial,
    ops: FVOperators,
    regularization_fraction: float = 1.0e-9,
) -> np.ndarray:
    """Evaluate the explicit nonlinear gTDGL forcing F[Delta]."""

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
    max_phase_step_rad: float = 0.25,
    use_phi_phase: bool = False,
) -> tuple[np.ndarray, bool, float]:
    """Advance Delta by one explicit KWT relaxation step.

    OE7 is a stationary fixed-bias relaxation problem, not a voltage-driven
    time-domain detector simulation. Poisson is retained to compute the
    normal-current correction

        j_n = -sigma_n grad(phi),

    and to enforce the conservative FV condition

        div_h(j_s + j_n) + b_h = 0.

    However, the Poisson potential is not used by default as an additional
    Josephson phase rotator during OE7. The superconducting momentum Q is
    reconstructed from the phase of psi itself:

        Q_ij = Arg(psi_j psi_i*) / l_ij.

    Therefore applying exp[-i 2e phi dt / hbar] here would feed the Poisson
    projection back into Q and can artificially drive pairbreaking in a
    nominally stationary below-Ic state.

    Set use_phi_phase=True only for a genuinely voltage-driven dynamic run,
    not for the OE7 stationary branch.
    """
    psi = np.asarray(psi_J, dtype=np.complex128)
    phi = np.asarray(phi_V, dtype=float)
    Te = np.asarray(Te_K, dtype=float)
    forcing = np.asarray(forcing_J, dtype=np.complex128)

    if dt_s <= 0.0:
        raise ValueError("dt_s must be positive.")

    R = np.abs(psi)
    rho = material.rho_kwt(Te, R)
    rho = np.maximum(rho, 1.0e-30)

    psi_euler = psi + (dt_s / material.tau0_GL_s) * forcing / rho

    if use_phi_phase:
        phase_step = (2.0 * E_CHARGE_C / HBAR_J_S) * phi * dt_s
        max_abs_phase = float(np.max(np.abs(phase_step))) if phase_step.size else 0.0

        if max_abs_phase > max_phase_step_rad:
            return psi.copy(), False, max_abs_phase

        psi_new = np.exp(-1j * phase_step) * psi_euler
    else:
        max_abs_phase = 0.0
        psi_new = psi_euler

    if (
        not np.all(np.isfinite(np.real(psi_new)))
        or not np.all(np.isfinite(np.imag(psi_new)))
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
    """Apply notebook-style OE7 stationary boundary conditions.

    Insulating edges:
        Delta_boundary = Delta_inner.

    Terminals:
        |Delta_boundary| = |Delta_inner|,
        theta_boundary = theta_inner +/- q_bias dx.

    Terminal constraints are applied after top/bottom constraints, so terminal
    physics wins at corner nodes.
    """

    psi = np.asarray(psi_trial_J, dtype=np.complex128)
    out = np.array(psi, dtype=np.complex128, copy=True)

    out = (np.nan_to_num(np.real(out), nan=0.0, posinf=0.0, neginf=0.0) + 1j * np.nan_to_num(np.imag(out), nan=0.0, posinf=0.0, neginf=0.0))
    out = clip_gap_amplitude(out, material)

    if not enabled:
        return out

    for side in ("bottom", "top"):
        boundary, inner = nearest_inward_boundary_pairs(mesh, side)
        out[boundary] = out[inner]

    nodes = np.asarray(mesh.nodes, dtype=float)

    pairs = terminal_inner_node_pairs(mesh)

    left_boundary, left_inner = pairs["left"]
    dx_left = np.abs(nodes[left_inner, 0] - nodes[left_boundary, 0])
    amp_left = np.abs(out[left_inner])
    theta_left = np.angle(out[left_inner]) - float(q_bias_m_inv) * dx_left
    out[left_boundary] = amp_left * np.exp(1j * theta_left)

    right_boundary, right_inner = pairs["right"]
    dx_right = np.abs(nodes[right_boundary, 0] - nodes[right_inner, 0])
    amp_right = np.abs(out[right_inner])
    theta_right = np.angle(out[right_inner]) + float(q_bias_m_inv) * dx_right
    out[right_boundary] = amp_right * np.exp(1j * theta_right)

    out = clip_gap_amplitude(out, material)
    return out


def clip_gap_amplitude(psi_J: np.ndarray, material: GTDGLMaterial) -> np.ndarray:
    """Clip only nonphysical overshoots, never enforce a positive floor."""

    psi = np.asarray(psi_J, dtype=np.complex128)
    amp = np.abs(psi)
    phase = np.exp(1j * np.angle(psi))

    amp = np.clip(amp, 0.0, 1.5 * material.delta0_J)
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


def terminal_node_mask(mesh) -> np.ndarray:
    """Return boolean mask for left and right terminal nodes."""

    masks = boundary_node_masks(mesh)
    return masks["left"] | masks["right"]


def terminal_inner_node_pairs(mesh) -> dict[str, tuple[np.ndarray, np.ndarray]]:
    """Pair left/right terminal nodes with their immediate inward neighbours."""

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
            score = (dy / max(h, 1.0e-300)) ** 2 + 0.05 * (
                dx / max(h, 1.0e-300)
            ) ** 2
            inner[k] = int(candidates[int(np.argmin(score))])

        return boundary.astype(np.int64), inner

    return {
        "left": pair("left"),
        "right": pair("right"),
    }


def nearest_inward_boundary_pairs(mesh, side: str) -> tuple[np.ndarray, np.ndarray]:
    """Pair top/bottom boundary nodes with their immediate inward neighbours."""

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
        score = (dx / max(h, 1.0e-300)) ** 2 + 0.05 * (
            dy / max(h, 1.0e-300)
        ) ** 2
        inner[k] = int(candidates[int(np.argmin(score))])

    return boundary.astype(np.int64), inner


def current_residual(currents: CurrentFields, mesh) -> float:
    """Dimensionless RMS residual of div(j_tot)."""

    div = np.asarray(currents.node_div_jtot_A_m3, dtype=float)
    jmag = np.sqrt(
        currents.node_jtot_x_A_m2**2
        + currents.node_jtot_y_A_m2**2
    )

    length_scale = max(
        float(getattr(mesh, "target_spacing_m", 1.0e-9)),
        1.0e-30,
    )
    scale = max(float(np.nanmean(jmag)) / length_scale, 1.0)

    return float(np.sqrt(np.nanmean(div * div)) / scale)


def normal_current_fraction_rms(currents: CurrentFields) -> float:
    """RMS normal-current fraction relative to total-current scale."""

    jn2 = currents.edge_jn_A_m2**2
    jt2 = currents.edge_jtot_A_m2**2

    num = float(np.sqrt(np.nanmean(jn2))) if jn2.size else 0.0
    den = float(np.sqrt(np.nanmean(jt2))) if jt2.size else 0.0

    return num / max(den, 1.0e-300)


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
    """Extract the seed phase-gradient q.

    For zero-current unit tests, missing q is interpreted as zero.
    """

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
    target_current_A: float | None = None,
    progress: bool = False,
    n_phi_snapshots: int = 6,
    use_phi_phase: bool = False,
) -> RelaxationResult:
    """Relax the OE6 seed with frozen temperatures and active Poisson."""

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

    psi0 = (
        np.asarray(seed.node_psi_real_J, dtype=float)
        + 1j * np.asarray(seed.node_psi_imag_J, dtype=float)
    )
    phi = np.asarray(seed.node_phi_electric_V, dtype=float).copy()
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

    boundary_accum = target_terminal_boundary_accum_A_m(
        edge_data=edge_data,
        ops=ops,
        material=material,
        target_current_A=target_current_A,
    )

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

    hist_t: list[float] = []
    hist_dt: list[float] = []
    hist_eta: list[float] = []
    hist_res: list[float] = []
    hist_v: list[float] = []
    hist_pairbreaking_max: list[float] = []
    hist_delta_min_ratio: list[float] = []
    hist_normal_fraction: list[float] = []

    n_phi_snapshots = max(2, int(n_phi_snapshots))
    phi_snapshot_t_s: list[float] = [0.0]
    phi_snapshot_V: list[np.ndarray] = [phi.copy()]
    phi_snapshot_steps = set(
        np.unique(
            np.rint(np.linspace(1, int(steps), n_phi_snapshots - 1)).astype(int)
        ).tolist()
    )

    iterator = range(int(steps))
    if progress and trange is not None:
        iterator = trange(int(steps), desc="OE7 SS gTDGL", leave=True)

    for _ in iterator:
        forcing = gtdgl_forcing(
            psi_J=psi,
            Te_K=Te,
            currents=currents,
            material=material,
            ops=ops,
        )

        psi_trial, ok, update_measure = kwt_local_update(
            psi_J=psi,
            phi_V=phi,
            Te_K=Te,
            forcing_J=forcing,
            dt_s=dt_s,
            material=material,
            use_phi_phase=use_phi_phase,
        )

        if not ok:
            rejected += 1
            if adapt_dt and dt_s > dt_min_s:
                dt_s = max(dt_min_s, 0.5 * dt_s)
                continue
            raise FloatingPointError(
                "KWT update failed; "
                f"max gauge phase step = {update_measure:.6e} rad."
            )

        psi_trial = apply_stationary_boundary_conditions(
            psi_trial_J=psi_trial,
            mesh=mesh,
            seed=seed,
            q_bias_m_inv=q_bias,
            material=material,
            enabled=lock_terminals,
        )

        trial_currents_no_phi = compute_current_fields(
            psi_J=psi_trial,
            phi_V=np.zeros_like(phi),
            Te_K=Te,
            material=material,
            ops=ops,
            boundary_accum_A_m=None,
        )

        phi_trial = solve_poisson_potential(
            edge_js_us_A_m2=trial_currents_no_phi.edge_js_us_A_m2,
            material=material,
            ops=ops,
            boundary_accum_A_m=boundary_accum,
        )

        trial_currents = compute_current_fields(
            psi_J=psi_trial,
            phi_V=phi_trial,
            Te_K=Te,
            material=material,
            ops=ops,
            boundary_accum_A_m=boundary_accum,
        )

        eta = float(
            np.max(np.abs(np.abs(psi_trial) ** 2 - np.abs(psi) ** 2))
            / material.delta0_J**2
        )

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
        voltage = terminal_voltage(
            np.asarray(mesh.nodes, dtype=float),
            phi,
            length_m=float(mesh.length_m),
        )
        pb_max = float(np.nanmax(currents.node_pairbreaking_ratio))
        delta_min_ratio = float(np.nanmin(np.abs(psi)) / material.delta0_J)
        normal_frac = normal_current_fraction_rms(currents)

        hist_t.append(t_s)
        hist_dt.append(dt_s)
        hist_eta.append(eta)
        hist_res.append(residual)
        hist_v.append(voltage)
        hist_pairbreaking_max.append(pb_max)
        hist_delta_min_ratio.append(delta_min_ratio)
        hist_normal_fraction.append(normal_frac)

        if accepted in phi_snapshot_steps:
            phi_snapshot_t_s.append(t_s)
            phi_snapshot_V.append(phi.copy())

        if progress and hasattr(iterator, "set_postfix") and accepted % 10 == 0:
            iterator.set_postfix(
                eta=f"{eta:.2e}",
                eps=f"{residual:.2e}",
                V=f"{voltage:.2e}",
                chi=f"{pb_max:.3g}",
                dt_fs=f"{dt_s / 1.0e-15:.3g}",
            )

        if accepted >= min_steps and eta < tolerance_eta and residual < tolerance_current_residual:
            converged = True
            break

        if adapt_dt and eta < 0.1 * tolerance_eta:
            dt_s = min(dt_max_s, 1.2 * dt_s)

    if len(phi_snapshot_t_s) == 0 or phi_snapshot_t_s[-1] != t_s:
        phi_snapshot_t_s.append(t_s)
        phi_snapshot_V.append(phi.copy())

    if len(phi_snapshot_t_s) > n_phi_snapshots:
        keep = np.unique(
            np.rint(np.linspace(0, len(phi_snapshot_t_s) - 1, n_phi_snapshots)).astype(
                int
            )
        )
        if keep[-1] != len(phi_snapshot_t_s) - 1:
            keep[-1] = len(phi_snapshot_t_s) - 1

        phi_snapshot_t_s = [phi_snapshot_t_s[int(i)] for i in keep]
        phi_snapshot_V = [phi_snapshot_V[int(i)] for i in keep]

    boundary = boundary_currents_from_node_vectors(
        mesh=mesh,
        edge_data=edge_data,
        jx_A_m2=currents.node_jtot_x_A_m2,
        jy_A_m2=currents.node_jtot_y_A_m2,
        thickness_m=material.thickness_m,
    )

    voltage = terminal_voltage(
        np.asarray(mesh.nodes, dtype=float),
        phi,
        length_m=float(mesh.length_m),
    )

    j_bias = (
        float(target_current_A)
        / (material.width_m * material.thickness_m)
        if material.width_m > 0.0 and material.thickness_m > 0.0
        else float("nan")
    )
    normal_ohmic_voltage = (
        float(target_current_A)
        * float(mesh.length_m)
        / (
            material.sigma_n_S_m
            * material.width_m
            * material.thickness_m
        )
    )

    summary = {
        "backend": "oe7_stationary_gtdgl_poisson_fixed_phase_gauge_v2",
        "gauge_policy": "poisson_retained_but_phi_not_used_as_stationary_phase_rotator",
        "use_phi_phase": bool(use_phi_phase),
        "converged": bool(converged),
        "accepted_steps": int(accepted),
        "rejected_steps": int(rejected),
        "final_time_ps": float(t_s / 1.0e-12),
        "tau_scale": float(material.tau_scale),
        "tau_ee_Tc_effective_ps": float(material.tau_scale * material.tau_ee_Tc_s / 1.0e-12),
        "tau_ep_Tc_effective_ps": float(material.tau_scale * material.tau_ep_Tc_s / 1.0e-12),
        "target_current_A": float(target_current_A),
        "target_q_bias_m_inv": float(q_bias),
        "target_j_bias_A_m2": float(j_bias),
        "terminal_voltage_V": float(voltage),
        "normal_ohmic_voltage_V": float(normal_ohmic_voltage),
        "terminal_voltage_over_normal": float(
            voltage / normal_ohmic_voltage
            if normal_ohmic_voltage != 0.0
            else float("nan")
        ),
        "normal_current_fraction_rms": float(normal_current_fraction_rms(currents)),
        "current_residual": float(current_residual(currents, mesh)),
        "eta_R_final": float(hist_eta[-1]) if hist_eta else float("nan"),
        "divergence_rms_A_m3": float(
            np.sqrt(np.nanmean(currents.node_div_jtot_A_m3**2))
        ),
        "min_delta_over_delta0": float(np.nanmin(np.abs(psi)) / material.delta0_J),
        "mean_delta_over_delta0": float(np.nanmean(np.abs(psi)) / material.delta0_J),
        "max_pairbreaking_ratio": float(np.nanmax(currents.node_pairbreaking_ratio)),
        "p99_pairbreaking_ratio": float(
            np.nanpercentile(currents.node_pairbreaking_ratio, 99.0)
        ),
        "edge_Q_max_m_inv": float(np.nanmax(np.abs(currents.edge_Q_m_inv))),
        "boundary_currents_A": boundary,
    }

    metadata = {
        "backend": summary["backend"],
        "description": (
            "Frozen-temperature stationary gTDGL/Poisson relaxation from the OE6 "
            "analytic seed. Poisson is retained for current conservation and "
            "normal-current diagnostics, but the stationary OE7 branch does not "
            "apply phi as an additional Josephson phase rotator."
        ),
        "thermal_policy": "frozen_Te_Tph",
        "circuit_policy": "inactive",
        "boundary_policy": "notebook_style_q_bias_phase_continuation",
        "poisson_policy": "conservative_FV_mean_zero_gauge",
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

    history = {
        "t_s": np.asarray(hist_t, dtype=float),
        "dt_s": np.asarray(hist_dt, dtype=float),
        "eta_R": np.asarray(hist_eta, dtype=float),
        "current_residual": np.asarray(hist_res, dtype=float),
        "terminal_voltage_V": np.asarray(hist_v, dtype=float),
        "pairbreaking_max": np.asarray(hist_pairbreaking_max, dtype=float),
        "delta_min_over_delta0": np.asarray(hist_delta_min_ratio, dtype=float),
        "normal_current_fraction_rms": np.asarray(hist_normal_fraction, dtype=float),
        "phi_snapshot_t_s": np.asarray(phi_snapshot_t_s, dtype=float),
        "phi_snapshot_V": np.asarray(phi_snapshot_V, dtype=float),
    }

    return RelaxationResult(state=state, history=history, summary=summary)


def save_stationary_state_npz(
    state: GTDGLStationaryState,
    output_path: str | Path,
) -> Path:
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


def save_relaxation_history_npz(
    history: dict[str, np.ndarray],
    output_path: str | Path,
) -> Path:
    """Save compact relaxation history to NPZ."""

    output = Path(output_path)
    output.parent.mkdir(parents=True, exist_ok=True)

    arrays = {key: np.asarray(value) for key, value in history.items()}
    np.savez_compressed(output, **arrays)

    return output