"""Poisson projection and terminal-flux helpers for OE7."""
from __future__ import annotations

import numpy as np

try:
    from scipy.sparse import coo_matrix, csr_matrix, bmat
    from scipy.sparse.linalg import splu, spsolve
except Exception:  # pragma: no cover
    coo_matrix = None
    csr_matrix = None
    bmat = None
    splu = None
    spsolve = None

from pysnspd.gtdgl.material import GTDGLMaterial
from pysnspd.gtdgl.operators import (
    FVOperators,
    divergence_from_edge_scalar,
    edge_flux_accumulator_A_m,
    edge_scalar_gradient,
    terminal_boundary_accum_A_m,
)
from pysnspd.gtdgl.state import PhiBoundaryConditions, PoissonResult, _PoissonOperator
from pysnspd.gtdgl.stationary_boundary import (
    boundary_node_masks,
    nearest_inward_boundary_pairs,
    terminal_inner_node_pairs,
)
from pysnspd.gtdgl.diagnostics import target_current_density_A_m2

VALID_POISSON_TERMINAL_POLICIES = {
    "target_flux",
    "zero_flux",
}


def build_poisson_operator(
    *,
    material: GTDGLMaterial,
    ops: FVOperators,
    phi_bc: PhiBoundaryConditions | None = None,
) -> _PoissonOperator:
    """Build the pyTDGL-style sparse-LU Poisson operator.

    This is the physical-units analogue of pyTDGL Eq. (17):

        sum_j s_ij * (mu_j - mu_i) / e_ij = sum_j s_ij * J_s,ij

    with ``mu`` replaced by the electric potential ``phi`` in volts and
    ``J_n,ij = -sigma_n * (phi_j - phi_i) / e_ij``.  Equivalently,

        sum_j sigma_n * s_ij/e_ij * (phi_i - phi_j)
        = - sum_j s_ij * j_s,ij - F_i^ext.

    The operator is pure Neumann and therefore singular up to a constant; we
    append one gauge row/column enforcing mean(phi)=0.  The optional ``phi_bc``
    argument is accepted only for backward compatibility with the previous
    experimental branch; it is deliberately ignored here.
    """
    del phi_bc

    sigma = material.sigma_n_S_m
    g = sigma * ops.dual_face_length_m / np.maximum(ops.edge_length_m, 1.0e-300)
    i = ops.edge_i.astype(np.int64)
    j = ops.edge_j.astype(np.int64)
    n = int(ops.n_nodes)

    # Row i gets +g*(phi_i - phi_j); row j gets +g*(phi_j - phi_i).
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


def _build_constrained_poisson_operator(
    *,
    material: GTDGLMaterial,
    ops: FVOperators,
    phi_bc: PhiBoundaryConditions,
) -> _PoissonOperator:
    """Deprecated experimental constrained-Poisson builder.

    Retained only so old imports do not break.  The OE7 production path now
    follows pyTDGL Eq. (17) and no longer calls this function.
    """
    del phi_bc
    return build_poisson_operator(material=material, ops=ops, phi_bc=None)


def solve_varphi_poisson(
    *,
    edge_js_us_A_m2: np.ndarray,
    material: GTDGLMaterial,
    ops: FVOperators,
    poisson_op: _PoissonOperator | None = None,
    boundary_accum_A_m: np.ndarray | None = None,
    phi_bc: PhiBoundaryConditions | None = None,
) -> PoissonResult:
    """Solve the pyTDGL Eq. (17) Poisson projection for ``phi``.

    Given the already-updated supercurrent ``j_s^{n+1}``, solve the finite-
    volume equation

        sum_j sigma_n * s_ij/e_ij * (phi_i - phi_j)
        = - sum_j s_ij * j_s,ij - F_i^ext,

    by sparse LU factorization of the mean-zero Neumann system.  The normal
    current is then evaluated only after the solve as

        j_n,ij = -sigma_n * (phi_j - phi_i) / e_ij.

    This intentionally does not impose any first-edge electric constraint.
    ``phi_bc`` is accepted for backward compatibility and ignored.
    """
    del phi_bc

    js = np.asarray(edge_js_us_A_m2, dtype=float)
    if js.shape != (ops.n_edges,):
        raise ValueError(f"edge_js_us_A_m2 must have shape ({ops.n_edges},).")

    if poisson_op is None:
        poisson_op = build_poisson_operator(material=material, ops=ops, phi_bc=None)

    if boundary_accum_A_m is None:
        boundary = np.zeros(ops.n_nodes, dtype=float)
    else:
        boundary = np.asarray(boundary_accum_A_m, dtype=float)
        if boundary.shape != (ops.n_nodes,):
            raise ValueError(
                f"boundary_accum_A_m must have shape ({ops.n_nodes},), got {boundary.shape}."
            )

    # Eq. (17) in physical units.  ``edge_flux_accumulator_A_m(js, ops)``
    # returns sum_j s_ij*j_s,ij with the edge orientation convention.
    # ``boundary`` is the prescribed outward flux through true terminals.
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
    phi_bc: PhiBoundaryConditions | None = None,
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
        phi_bc=phi_bc,
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


def build_phi_boundary_conditions(
    *,
    mesh,
    ops: FVOperators,
    material: GTDGLMaterial,
    seed,
    target_current_A: float,
    enabled: bool = True,
) -> PhiBoundaryConditions | None:
    """Build first-normal-edge electric BCs for varphi.

    The constraints enforce the total current on the same inward normal edges
    used by the Delta boundary condition:

    * left/right: j_tot,bk = j_avg * (dx/ell), laminar longitudinal injection;
    * top/bottom: j_tot,bk = 0, insulating boundary.

    Corners are assigned to the longitudinal terminals first; top/bottom
    constraints exclude left/right nodes to avoid over-constraining a corner.
    """
    if not enabled:
        return None

    nodes = np.asarray(mesh.nodes, dtype=float)
    n_nodes = int(nodes.shape[0])
    masks = boundary_node_masks(mesh)
    terminal_mask = masks["left"] | masks["right"]
    javg = target_current_density_A_m2(material, float(target_current_A))

    boundary_parts: list[np.ndarray] = []
    inner_parts: list[np.ndarray] = []
    target_parts: list[np.ndarray] = []

    terminal_pairs = terminal_inner_node_pairs(mesh, ops=ops)
    for side in ("left", "right"):
        b, k = terminal_pairs[side]
        dx = nodes[k, 0] - nodes[b, 0]
        dy = nodes[k, 1] - nodes[b, 1]
        ell = np.maximum(np.sqrt(dx * dx + dy * dy), 1.0e-300)
        boundary_parts.append(np.asarray(b, dtype=np.int64))
        inner_parts.append(np.asarray(k, dtype=np.int64))
        target_parts.append(javg * dx / ell)

    for side in ("bottom", "top"):
        b, k = nearest_inward_boundary_pairs(mesh, side, ops=ops)
        keep = ~terminal_mask[b]
        if np.any(keep):
            boundary_parts.append(np.asarray(b[keep], dtype=np.int64))
            inner_parts.append(np.asarray(k[keep], dtype=np.int64))
            target_parts.append(np.zeros(int(np.count_nonzero(keep)), dtype=float))

    if not boundary_parts:
        return None

    boundary = np.concatenate(boundary_parts).astype(np.int64)
    inner = np.concatenate(inner_parts).astype(np.int64)
    target = np.concatenate(target_parts).astype(float)

    # Keep the first constraint for each boundary node.  This preserves the
    # terminal priority at corners and avoids duplicate matrix rows.
    seen: set[int] = set()
    keep_idx: list[int] = []
    for idx, b in enumerate(boundary.tolist()):
        if b in seen:
            continue
        seen.add(b)
        keep_idx.append(idx)
    keep_arr = np.asarray(keep_idx, dtype=np.int64)
    boundary = boundary[keep_arr]
    inner = inner[keep_arr]
    target = target[keep_arr]

    edge_index, edge_sign = _edge_indices_and_signs_for_pairs(
        ops=ops,
        boundary=boundary,
        inner=inner,
    )
    edge_length = np.asarray(ops.edge_length_m, dtype=float)[edge_index]
    boundary_mask = np.zeros(n_nodes, dtype=bool)
    boundary_mask[boundary] = True

    return PhiBoundaryConditions(
        boundary_nodes=boundary,
        inner_nodes=inner,
        edge_index=edge_index,
        edge_sign_b_to_inner=edge_sign,
        target_edge_A_m2=target,
        edge_length_m=edge_length,
        boundary_mask=boundary_mask,
    )


def _edge_indices_and_signs_for_pairs(
    *,
    ops: FVOperators,
    boundary: np.ndarray,
    inner: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Return edge index and orientation sign for boundary -> inner pairs."""
    edge_i = np.asarray(ops.edge_i, dtype=np.int64)
    edge_j = np.asarray(ops.edge_j, dtype=np.int64)
    lookup: dict[tuple[int, int], tuple[int, float]] = {}
    for idx, (a, b) in enumerate(zip(edge_i.tolist(), edge_j.tolist())):
        lookup[(int(a), int(b))] = (idx, +1.0)
        lookup[(int(b), int(a))] = (idx, -1.0)

    edge_index = np.empty(boundary.size, dtype=np.int64)
    edge_sign = np.empty(boundary.size, dtype=float)
    for m, (b, k) in enumerate(zip(boundary.tolist(), inner.tolist())):
        key = (int(b), int(k))
        if key not in lookup:
            raise ValueError(
                "Phi boundary pair is not an FV edge: "
                f"boundary={b}, inner={k}."
            )
        idx, sign = lookup[key]
        edge_index[m] = int(idx)
        edge_sign[m] = float(sign)
    return edge_index, edge_sign

