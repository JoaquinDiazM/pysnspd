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
try:
    from tqdm.auto import trange
except Exception:  # pragma: no cover
    trange = None

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
    unwrap_phase_graph,
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
    boundary_accum_A_m: np.ndarray | None = None,
) -> CurrentFields:
    """Evaluate Usadel-like, GL, normal and total currents.

    The superconducting and normal currents live on edges. The total divergence
    can optionally include prescribed terminal fluxes through
    ``boundary_accum_A_m``.
    """
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
    div_total = divergence_from_edge_scalar(
        edge_jtot,
        ops,
        boundary_accum_A_m=boundary_accum_A_m,
    )

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
    edge_data=None,
    target_current_A: float | None = None,
    boundary_accum_A_m: np.ndarray | None = None,
) -> np.ndarray:
    """Solve Poisson with conservative terminal-current boundary terms.

    The discrete condition is

        div_edges(j_s + j_n) + div_boundary(j_terminal) = 0,

    with j_n = -sigma_n grad(phi). This prevents Poisson from cancelling the
    imposed transport current.
    """
    js = np.asarray(edge_js_us_A_m2, dtype=float)
    conductance = material.sigma_n_S_m * ops.dual_face_length_m / ops.edge_length_m
    i = ops.edge_i
    j = ops.edge_j

    rhs = np.zeros(ops.n_nodes, dtype=float)

    edge_flux = ops.dual_face_length_m * js
    np.add.at(rhs, i, -edge_flux)
    np.add.at(rhs, j, edge_flux)

    if boundary_accum_A_m is not None:
        boundary_accum = np.asarray(boundary_accum_A_m, dtype=float)
        if boundary_accum.shape != rhs.shape:
            raise ValueError(
                "boundary_accum_A_m must have shape "
                f"{rhs.shape}, got {boundary_accum.shape}."
            )
        rhs -= boundary_accum
    elif target_current_A is not None:
        if edge_data is None:
            raise ValueError("edge_data is required when target_current_A is used.")
        boundary_accum = terminal_boundary_accum_A_m(
            edge_data,
            n_nodes=ops.n_nodes,
            target_current_A=float(target_current_A),
            thickness_m=material.thickness_m,
        )
        rhs -= boundary_accum

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

def terminal_inner_node_pairs(mesh) -> dict[str, tuple[np.ndarray, np.ndarray]]:
    """Pair terminal nodes with their immediate inward neighbours.

    The current mesh has protected, non-jittered boundary layers, so this gives
    a clean terminal normal direction while still being robust to small jitter.
    """
    nodes = np.asarray(mesh.nodes, dtype=float)
    x = nodes[:, 0]
    y = nodes[:, 1]

    xmin = float(np.min(x))
    xmax = float(np.max(x))
    h = float(getattr(mesh, "target_spacing_m", 1.0e-9))
    tol = max(1.0e-15, 1.0e-9 * float(mesh.length_m))

    left = np.where(np.abs(x - xmin) <= tol)[0]
    right = np.where(np.abs(x - xmax) <= tol)[0]

    def pair_one_side(boundary: np.ndarray, *, side: str) -> tuple[np.ndarray, np.ndarray]:
        inner = np.empty_like(boundary, dtype=np.int64)

        for k, b in enumerate(boundary):
            if side == "left":
                candidates = np.where(x > x[b] + tol)[0]
                dx = x[candidates] - x[b]
            elif side == "right":
                candidates = np.where(x < x[b] - tol)[0]
                dx = x[b] - x[candidates]
            else:
                raise ValueError(f"Unknown side {side!r}.")

            if candidates.size == 0:
                raise ValueError(f"No inward candidates found for {side} terminal node {b}.")

            dy = y[candidates] - y[b]

            # Strongly prefer same-y first inward neighbours; fall back to nearest.
            score = (dy / max(h, 1.0e-300)) ** 2 + 0.05 * (dx / max(h, 1.0e-300)) ** 2
            inner[k] = int(candidates[int(np.argmin(score))])

        return boundary.astype(np.int64), inner

    return {
        "left": pair_one_side(left, side="left"),
        "right": pair_one_side(right, side="right"),
    }


def apply_terminal_supercurrent_bc(
    *,
    psi_trial_J: np.ndarray,
    Te_K: np.ndarray,
    mesh,
    material: GTDGLMaterial,
    target_current_A: float,
    enabled: bool = True,
) -> np.ndarray:
    """Apply superconducting terminal boundary condition.

    Terminal amplitude is copied from the immediate inward neighbour,

        R_boundary = R_inner,

    while the terminal phase is chosen so that

        j_s^Us(R_inner, T, q_b) = I_bias/(w d).

    This replaces the older hard Dirichlet condition Psi_terminal=Psi_seed.
    """
    psi = np.array(psi_trial_J, dtype=np.complex128, copy=True)
    if not enabled:
        return psi

    j_target = float(target_current_A) / (material.width_m * material.thickness_m)
    pairs = terminal_inner_node_pairs(mesh)

    for side, sign in (("left", -1.0), ("right", +1.0)):
        boundary, inner = pairs[side]

        R_inner = np.abs(psi[inner])
        theta_inner = np.angle(psi[inner])
        Te_inner = np.maximum(np.asarray(Te_K, dtype=float)[inner], 1.0e-12)

        coeff = (
            np.pi
            * material.sigma_n_S_m
            / (2.0 * E_CHARGE_C)
            * R_inner
            * np.tanh(R_inner / (2.0 * K_B_J_K * Te_inner))
        )

        if np.any(coeff <= 0.0) or not np.all(np.isfinite(coeff)):
            raise FloatingPointError(
                "Cannot impose superconducting terminal current because "
                "the local Usadel current coefficient is non-positive or non-finite."
            )

        q_b = j_target / coeff

        dx = np.abs(
            np.asarray(mesh.nodes, dtype=float)[boundary, 0]
            - np.asarray(mesh.nodes, dtype=float)[inner, 0]
        )

        # left: theta_boundary = theta_inner - q dx
        # right: theta_boundary = theta_inner + q dx
        theta_boundary = theta_inner + sign * q_b * dx

        psi[boundary] = R_inner * np.exp(1j * theta_boundary)

    return psi

def _edge_flux_accumulator_A_m(
    edge_current_i_to_j: np.ndarray,
    ops: FVOperators,
) -> np.ndarray:
    """Return conservative edge-flux accumulator before division by node area.

    Units are [A/m], i.e. current per film thickness. This is the numerator of
    the finite-volume divergence used by divergence_from_edge_scalar().
    """
    current = np.asarray(edge_current_i_to_j, dtype=float)
    if current.shape != (ops.n_edges,):
        raise ValueError(f"edge_current_i_to_j must have shape ({ops.n_edges},).")

    out = np.zeros(ops.n_nodes, dtype=float)
    flux = ops.dual_face_length_m * current
    np.add.at(out, ops.edge_i, flux)
    np.add.at(out, ops.edge_j, -flux)
    return out


def _supercurrent_terminal_boundary_accum_A_m(
    *,
    edge_js_us_A_m2: np.ndarray,
    mesh,
    ops: FVOperators,
) -> np.ndarray:
    """Build Poisson boundary accumulator compatible with imposed j_s.

    Option A: the transport current is imposed through the superconducting
    phase-gradient condition on Psi. Poisson must not impose I_bias again.
    Therefore the terminal boundary accumulator cancels only the discrete
    edge-flux imbalance of the already-imposed supercurrent at left/right
    terminal nodes.

    Interior nodes are left untouched, so Poisson can still correct genuine
    interior continuity errors through j_n = -sigma_n grad(phi).
    """
    raw = _edge_flux_accumulator_A_m(edge_js_us_A_m2, ops)
    terminal = terminal_node_mask(mesh)

    out = np.zeros(ops.n_nodes, dtype=float)
    out[terminal] = -raw[terminal]
    return out

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
    target_current_A: float | None = None,
    progress: bool = False,
    n_phi_snapshots: int = 6,
) -> RelaxationResult:
    """Relax the OE6 seed with frozen temperatures and active gTDGL/Poisson.

    Option A:
    - The terminal transport current is imposed through the superconducting
      phase-gradient boundary condition on Psi.
    - Poisson does not impose I_bias a second time.
    - Instead, Poisson uses a terminal accumulator compatible with the actual
      discrete supercurrent after the boundary condition has been applied.
    """
    if steps <= 0:
        raise ValueError("steps must be positive.")
    if dt_s <= 0.0:
        raise ValueError("dt_s must be positive.")
    if min_steps < 0:
        raise ValueError("min_steps must be non-negative.")

    if target_current_A is None:
        target_current_A = _seed_target_current_A(seed)
    target_current_A = float(target_current_A)

    psi0 = (
        np.asarray(seed.node_psi_real_J, dtype=float)
        + 1j * np.asarray(seed.node_psi_imag_J, dtype=float)
    )
    psi = psi0.copy()
    phi = np.asarray(seed.node_phi_electric_V, dtype=float).copy()
    Te = np.asarray(seed.node_Te_K, dtype=float).copy()
    Tph = np.asarray(seed.node_Tph_K, dtype=float).copy()

    # Initial compatible boundary accumulator from the initial supercurrent.
    seed_currents_no_phi = compute_current_fields(
        psi_J=psi,
        phi_V=np.zeros_like(phi),
        Te_K=Te,
        material=material,
        ops=ops,
        boundary_accum_A_m=None,
    )
    boundary_accum = _supercurrent_terminal_boundary_accum_A_m(
        edge_js_us_A_m2=seed_currents_no_phi.edge_js_us_A_m2,
        mesh=mesh,
        ops=ops,
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
    hist_ir: list[float] = []
    hist_il: list[float] = []

    n_phi_snapshots = max(2, int(n_phi_snapshots))
    phi_snapshot_t_s: list[float] = [0.0]
    phi_snapshot_V: list[np.ndarray] = [phi.copy()]
    phi_snapshot_steps = set(
        np.unique(
            np.rint(np.linspace(1, int(steps), n_phi_snapshots - 1)).astype(int)
        ).tolist()
    )

    currents = compute_current_fields(
        psi_J=psi,
        phi_V=phi,
        Te_K=Te,
        material=material,
        ops=ops,
        boundary_accum_A_m=boundary_accum,
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
            raise FloatingPointError(
                f"KWT update failed; min discriminant={min_disc:.6e}"
            )

        # Terminal BC: Neumann amplitude + imposed superconducting q_b.
        psi_trial = apply_terminal_supercurrent_bc(
            psi_trial_J=psi_trial,
            Te_K=Te,
            mesh=mesh,
            material=material,
            target_current_A=target_current_A,
            enabled=lock_terminals,
        )

        # First evaluate j_s from the trial order parameter.
        trial_currents_no_phi = compute_current_fields(
            psi_J=psi_trial,
            phi_V=np.zeros_like(phi),
            Te_K=Te,
            material=material,
            ops=ops,
            boundary_accum_A_m=None,
        )

        # Option A: boundary accumulator compatible with the already-imposed
        # superconducting boundary current, not an independent I_bias source.
        trial_boundary_accum = _supercurrent_terminal_boundary_accum_A_m(
            edge_js_us_A_m2=trial_currents_no_phi.edge_js_us_A_m2,
            mesh=mesh,
            ops=ops,
        )

        phi_trial = solve_poisson_potential(
            edge_js_us_A_m2=trial_currents_no_phi.edge_js_us_A_m2,
            material=material,
            ops=ops,
            boundary_accum_A_m=trial_boundary_accum,
        )

        trial_currents = compute_current_fields(
            psi_J=psi_trial,
            phi_V=phi_trial,
            Te_K=Te,
            material=material,
            ops=ops,
            boundary_accum_A_m=trial_boundary_accum,
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
        boundary_accum = trial_boundary_accum
        t_s += dt_s
        accepted += 1

        residual = current_residual(currents, mesh)
        voltage = terminal_voltage(
            np.asarray(mesh.nodes, dtype=float),
            phi,
            length_m=float(mesh.length_m),
        )
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

        if accepted in phi_snapshot_steps:
            phi_snapshot_t_s.append(t_s)
            phi_snapshot_V.append(phi.copy())

        if progress and hasattr(iterator, "set_postfix") and accepted % 10 == 0:
            iterator.set_postfix(
                eta=f"{eta:.2e}",
                eps=f"{residual:.2e}",
                V=f"{voltage:.2e}",
                dt_fs=f"{dt_s / 1.0e-15:.3g}",
            )

        if (
            accepted >= min_steps
            and eta < tolerance_eta
            and residual < tolerance_current_residual
        ):
            converged = True
            break

        if adapt_dt and eta < 0.1 * tolerance_eta:
            dt_s = min(dt_max_s, 1.2 * dt_s)

    if not phi_snapshot_t_s or phi_snapshot_t_s[-1] != t_s:
        phi_snapshot_t_s.append(t_s)
        phi_snapshot_V.append(phi.copy())

    if len(phi_snapshot_t_s) > n_phi_snapshots:
        keep = np.unique(
            np.rint(np.linspace(0, len(phi_snapshot_t_s) - 1, n_phi_snapshots)).astype(int)
        )
        if keep[-1] != len(phi_snapshot_t_s) - 1:
            keep[-1] = len(phi_snapshot_t_s) - 1

        phi_snapshot_t_s = [phi_snapshot_t_s[int(i)] for i in keep]
        phi_snapshot_V = [phi_snapshot_V[int(i)] for i in keep]

    metadata = {
        "backend": "oe7_stationary_gtdgl_poisson_v3_supercurrent_terminal_bc",
        "description": (
            "Frozen-temperature stationary gTDGL/Poisson relaxation from the OE6 "
            "analytic seed. The transport current is imposed through the "
            "superconducting terminal phase-gradient boundary condition. Poisson "
            "uses a terminal accumulator compatible with that discrete supercurrent, "
            "so I_bias is not imposed a second time. External circuit and thermal "
            "evolution are inactive."
        ),
        "accepted_steps": int(accepted),
        "rejected_steps": int(rejected),
        "requested_steps": int(steps),
        "converged": bool(converged),
        "final_time_s": float(t_s),
        "target_current_A": float(target_current_A),
        "tau_scale": float(material.tau_scale),
        "tau_ee_Tc_effective_s": float(material.tau_scale * material.tau_ee_Tc_s),
        "tau_ep_Tc_effective_s": float(material.tau_scale * material.tau_ep_Tc_s),
        "lock_terminals": bool(lock_terminals),
        "thermal_policy": "frozen_Te_Tph",
        "circuit_policy": "inactive",
        "poisson_gauge": "mean_zero",
        "poisson_boundary_policy": "supercurrent_compatible_terminal_flux",
        "terminal_order_parameter_policy": (
            "Neumann amplitude from immediate inward neighbour plus "
            "Usadel supercurrent phase-gradient condition"
        ),
        "n_phi_snapshots": int(len(phi_snapshot_t_s)),
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
        "phi_snapshot_t_s": np.asarray(phi_snapshot_t_s, dtype=float),
        "phi_snapshot_V": np.vstack(phi_snapshot_V),
    }

    summary = stationary_summary(
        mesh=mesh,
        edge_data=edge_data,
        state=state,
        material=material,
        history=history,
    )

    return RelaxationResult(
        state=state,
        history=history,
        summary=summary,
    )
def _seed_target_current_A(seed) -> float:
    """Extract the imposed transport current from an OE6 seed object."""
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

    raise AttributeError(
        "Could not infer target_current_A from seed. Pass target_current_A "
        "explicitly or ensure the seed contains I_bias_A."
    )

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
    theta = unwrap_phase_graph(state.psi_J, edge_data.edges)
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

    jtot_rms = float(
        np.sqrt(
            np.mean(
                state.currents.node_jtot_x_A_m2**2
                + state.currents.node_jtot_y_A_m2**2
            )
        )
    )
    jn_rms = float(
        np.sqrt(
            np.mean(
                state.currents.node_jn_x_A_m2**2
                + state.currents.node_jn_y_A_m2**2
            )
        )
    )
    js_rms = float(
        np.sqrt(
            np.mean(
                state.currents.node_js_us_x_A_m2**2
                + state.currents.node_js_us_y_A_m2**2
            )
        )
    )
    normal_current_fraction_rms = jn_rms / max(jtot_rms, 1.0e-300)

    target_current_A = float(state.metadata.get("target_current_A", np.nan))
    normal_resistance_ohm = (
        float(mesh.length_m)
        / (material.sigma_n_S_m * material.width_m * material.thickness_m)
    )
    normal_ohmic_voltage_V = target_current_A * normal_resistance_ohm
    terminal_voltage_over_normal = (
        float(voltage) / normal_ohmic_voltage_V
        if np.isfinite(normal_ohmic_voltage_V) and abs(normal_ohmic_voltage_V) > 0.0
        else float("nan")
    )

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
        "normal_resistance_ohm": float(normal_resistance_ohm),
        "normal_ohmic_voltage_V": float(normal_ohmic_voltage_V),
        "terminal_voltage_over_normal": float(terminal_voltage_over_normal),
        "phi_min_V": float(np.min(state.phi_V)),
        "phi_max_V": float(np.max(state.phi_V)),
        "delta_min_meV": float(np.min(R) / 1.602176634e-22),
        "delta_max_meV": float(np.max(R) / 1.602176634e-22),
        "theta_min_rad": float(np.min(theta)),
        "theta_max_rad": float(np.max(theta)),
        "jx_mean_A_m2": float(np.mean(state.currents.node_jtot_x_A_m2)),
        "jy_mean_A_m2": float(np.mean(state.currents.node_jtot_y_A_m2)),
        "js_us_rms_A_m2": float(js_rms),
        "jn_rms_A_m2": float(jn_rms),
        "jtot_rms_A_m2": float(jtot_rms),
        "normal_current_fraction_rms": float(normal_current_fraction_rms),
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

def _seed_compatible_boundary_accum_A_m(
    *,
    seed,
    material: GTDGLMaterial,
    ops: FVOperators,
) -> np.ndarray:
    """Build a discrete terminal accumulator compatible with the OE6 seed.

    This makes the analytic seed satisfy the same finite-volume continuity
    equation used by OE7 with phi=0. It avoids generating a spurious linear
    electrostatic potential just to compensate small mesh/operator mismatches.
    """
    psi0 = (
        np.asarray(seed.node_psi_real_J, dtype=float)
        + 1j * np.asarray(seed.node_psi_imag_J, dtype=float)
    )
    phi0 = np.zeros(ops.n_nodes, dtype=float)
    Te0 = np.asarray(seed.node_Te_K, dtype=float)

    seed_currents = compute_current_fields(
        psi_J=psi0,
        phi_V=phi0,
        Te_K=Te0,
        material=material,
        ops=ops,
        boundary_accum_A_m=None,
    )

    div_seed_no_boundary = divergence_from_edge_scalar(
        seed_currents.edge_js_us_A_m2,
        ops,
        boundary_accum_A_m=None,
    )

    return -div_seed_no_boundary * ops.node_area_m2