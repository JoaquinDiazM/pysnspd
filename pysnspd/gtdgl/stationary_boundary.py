"""Stationary Delta boundary conditions for OE7 gTDGL."""
from __future__ import annotations

import numpy as np

from pysnspd.gtdgl.material import E_CHARGE_C, K_B_J_K, GTDGLMaterial
from pysnspd.gtdgl.operators import FVOperators
from pysnspd.gtdgl.diagnostics import (
    seed_delta_bias_J,
    seed_target_current_A,
    target_current_density_A_m2,
)

VALID_DELTA_BOUNDARY_POLICIES = {
    "current_inversion",
    "vacuum_only",
    "normal_terminal",
    "none",
}



def apply_delta_boundary_policy(
    *,
    psi_trial_J: np.ndarray,
    mesh,
    seed,
    q_bias_m_inv: float,
    material: GTDGLMaterial,
    ops: FVOperators | None = None,
    Te_K: np.ndarray | None = None,
    target_current_A: float | None = None,
    policy: str = "current_inversion",
) -> np.ndarray:
    """Apply selectable Delta boundary policy for OE7 diagnostics.

    Policies
    --------
    current_inversion
        Current OE7 policy: top/bottom vacuum plus left/right Usadel-current
        inversion on the first inward terminal edge.

    vacuum_only
        Keep only the superconducting-vacuum Neumann condition on top/bottom.
        Do not force terminal phase/current through Delta.

    normal_terminal
        pyTDGL-like terminal diagnostic: top/bottom vacuum and psi=0 on
        left/right terminal nodes. The transport current must then enter only
        through the Poisson terminal flux.

    none
        No Delta boundary projection beyond amplitude clipping.
    """
    policy = str(policy)

    if policy == "current_inversion":
        return apply_stationary_boundary_conditions(
            psi_trial_J=psi_trial_J,
            mesh=mesh,
            seed=seed,
            q_bias_m_inv=q_bias_m_inv,
            material=material,
            ops=ops,
            Te_K=Te_K,
            target_current_A=target_current_A,
            enabled=True,
        )

    out = np.asarray(psi_trial_J, dtype=np.complex128).copy()
    out = clip_gap_amplitude(out, material)

    if policy == "none":
        return out

    if policy not in {"vacuum_only", "normal_terminal"}:
        raise ValueError(
            "Unknown Delta boundary policy "
            f"{policy!r}. Expected current_inversion, vacuum_only, "
            "normal_terminal, or none."
        )

    for side in ("bottom", "top"):
        boundary, inner = nearest_inward_boundary_pairs(mesh, side, ops=ops)
        out[boundary] = out[inner]

    if policy == "normal_terminal":
        pairs = terminal_inner_node_pairs(mesh, ops=ops)
        for side in ("left", "right"):
            boundary, _inner = pairs[side]
            out[boundary] = 0.0 + 0.0j

    return clip_gap_amplitude(out, material)


def apply_stationary_boundary_conditions(
    *,
    psi_trial_J: np.ndarray,
    mesh,
    seed,
    q_bias_m_inv: float,
    material: GTDGLMaterial,
    ops: FVOperators | None = None,
    Te_K: np.ndarray | None = None,
    target_current_A: float | None = None,
    enabled: bool = True,
) -> np.ndarray:
    """Current-inverted boundary constraints for Delta.

    Top/bottom keep the vacuum condition Delta_boundary = Delta_inner on the
    first inward normal edge.

    Left/right do not impose an absolute phase and do not impose q_bias
    directly.  Instead, they invert the exact Usadel-like edge-current formula
    used later by ``compute_formula_fields``:

        j_s,e = K_e Q_e,
        K_e = pi*sigma_n/(2e) * R_e * tanh(R_e/(2 kB Te_e)),
        Q_e = (theta_inner - theta_boundary)/ell_e.

    For the first inward normal edge b -> k, impose the laminar target
    projection j_s,e = j_avg * d_x/ell_e.  Therefore

        theta_boundary = theta_inner - (j_avg/K_e) * d_x.

    The terminal amplitude is the predicted current-biased gap magnitude
    R_bias from the OE6/Usadel seed.  The global phase remains free because the
    boundary phase follows theta_inner.
    """
    psi = np.asarray(psi_trial_J, dtype=np.complex128)
    out = np.nan_to_num(psi, nan=0.0, posinf=0.0, neginf=0.0).astype(
        np.complex128,
        copy=True,
    )
    out = clip_gap_amplitude(out, material)
    if not enabled:
        return out

    # Vacuum Neumann top/bottom first, using actual inward mesh edges.
    for side in ("bottom", "top"):
        boundary, inner = nearest_inward_boundary_pairs(mesh, side, ops=ops)
        out[boundary] = out[inner]

    nodes = np.asarray(mesh.nodes, dtype=float)
    R_bc = seed_delta_bias_J(seed, fallback=float(np.nanmedian(np.abs(out))))
    R_bc = float(np.clip(R_bc, 0.0, 2.0 * material.delta0_J))

    if target_current_A is None:
        target_current_A = seed_target_current_A(seed)
    javg = target_current_density_A_m2(material, float(target_current_A))

    Te_nodes = _boundary_temperature_nodes(seed=seed, material=material, Te_K=Te_K, n_nodes=out.size)
    pairs = terminal_inner_node_pairs(mesh, ops=ops)

    for side in ("left", "right"):
        boundary, inner = pairs[side]
        dx = nodes[inner, 0] - nodes[boundary, 0]
        dy = nodes[inner, 1] - nodes[boundary, 1]
        ell = np.maximum(np.sqrt(dx * dx + dy * dy), 1.0e-300)

        R_inner = np.abs(out[inner])
        R_edge = np.maximum(0.5 * (R_bc + R_inner), 1.0e-30 * material.delta0_J)
        Te_edge = np.maximum(0.5 * (Te_nodes[boundary] + Te_nodes[inner]), 1.0e-12)

        K_edge = (
            np.pi
            * material.sigma_n_S_m
            / (2.0 * E_CHARGE_C)
            * R_edge
            * np.tanh(R_edge / (2.0 * K_B_J_K * Te_edge))
        )
        K_edge = np.maximum(K_edge, 1.0e-300)

        # Reverse-engineered from the plotted edge current:
        # j_s,e = K_e * (theta_inner - theta_boundary)/ell
        # and j_s,e_target = javg * dx/ell.
        theta_boundary = np.angle(out[inner]) - (javg / K_edge) * dx
        out[boundary] = R_bc * np.exp(1j * theta_boundary)

    return clip_gap_amplitude(out, material)


def _boundary_temperature_nodes(
    *,
    seed,
    material: GTDGLMaterial,
    Te_K: np.ndarray | None,
    n_nodes: int,
) -> np.ndarray:
    """Return node temperatures for boundary-current inversion."""
    if Te_K is not None:
        arr = np.asarray(Te_K, dtype=float).reshape(-1)
        if arr.size == n_nodes:
            return np.maximum(arr, 1.0e-12)

    if hasattr(seed, "node_Te_K"):
        arr = np.asarray(seed.node_Te_K, dtype=float).reshape(-1)
        if arr.size == n_nodes:
            return np.maximum(arr, 1.0e-12)

    metadata = getattr(seed, "metadata", None)
    if isinstance(metadata, dict):
        for name in ("T_bias_K", "Te_K", "temperature_K"):
            if name in metadata:
                value = float(metadata[name])
                if np.isfinite(value) and value > 0.0:
                    return np.full(n_nodes, value, dtype=float)

    return np.full(n_nodes, max(1.0e-12, 0.1 * material.Tc_K), dtype=float)


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


def terminal_inner_node_pairs(
    mesh,
    ops: FVOperators | None = None,
) -> dict[str, tuple[np.ndarray, np.ndarray]]:
    """Pair left/right terminal nodes with their first inward normal-edge node.

    The protected mesh has, for every terminal node, an actual edge that leaves
    the terminal approximately orthogonally.  We use that edge instead of a
    generic nearest-node search.  Diagonal and tangent terminal edges are still
    part of the FV operator, but they no longer define the imposed boundary
    phase drop.
    """
    return {
        "left": _edge_aware_boundary_pairs(mesh, "left", ops=ops),
        "right": _edge_aware_boundary_pairs(mesh, "right", ops=ops),
    }


def nearest_inward_boundary_pairs(
    mesh,
    side: str,
    ops: FVOperators | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    """Pair top/bottom boundary nodes with their first inward normal-edge node."""
    if side not in {"bottom", "top"}:
        raise ValueError(f"side must be 'bottom' or 'top', got {side!r}.")
    return _edge_aware_boundary_pairs(mesh, side, ops=ops)


def _edge_aware_boundary_pairs(
    mesh,
    side: str,
    ops: FVOperators | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    """Return boundary -> inward pairs, preferring actual orthogonal edges.

    If ``ops`` is available, candidates are restricted to true mesh neighbors
    of the boundary node.  Among those, we choose the inward edge with the
    smallest tangential displacement and shortest normal step.  This makes the
    longitudinal current-Neumann condition act on the first normal edge rather
    than on diagonal or tangent edges.

    The fallback branch preserves the previous nearest-node behavior for unit
    tests or legacy callers that do not pass ``ops``.
    """
    if side not in {"left", "right", "bottom", "top"}:
        raise ValueError(f"Invalid boundary side {side!r}.")

    nodes = np.asarray(mesh.nodes, dtype=float)
    x = nodes[:, 0]
    y = nodes[:, 1]
    masks = boundary_node_masks(mesh)
    boundary = np.where(masks[side])[0].astype(np.int64)
    inner = np.empty_like(boundary, dtype=np.int64)

    h = float(getattr(mesh, "target_spacing_m", 1.0e-9))
    length_m = float(getattr(mesh, "length_m", np.ptp(x)))
    width_m = float(getattr(mesh, "width_m", np.ptp(y)))
    tol = max(1.0e-15, 1.0e-9 * max(length_m, width_m))

    adjacency = _node_adjacency_from_ops_or_mesh(mesh, ops=ops)

    for k, b in enumerate(boundary):
        if adjacency is not None:
            candidates = np.asarray(sorted(adjacency.get(int(b), [])), dtype=np.int64)
        else:
            candidates = np.arange(nodes.shape[0], dtype=np.int64)

        if side == "left":
            candidates = candidates[x[candidates] > x[b] + tol]
            normal = x[candidates] - x[b]
            tangent = y[candidates] - y[b]
        elif side == "right":
            candidates = candidates[x[candidates] < x[b] - tol]
            normal = x[b] - x[candidates]
            tangent = y[candidates] - y[b]
        elif side == "bottom":
            candidates = candidates[y[candidates] > y[b] + tol]
            normal = y[candidates] - y[b]
            tangent = x[candidates] - x[b]
        else:  # top
            candidates = candidates[y[candidates] < y[b] - tol]
            normal = y[b] - y[candidates]
            tangent = x[candidates] - x[b]

        if candidates.size == 0 and adjacency is not None:
            # Robust fallback: if an externally supplied mesh lacks the expected
            # normal edge in adjacency, retry with all nodes rather than failing.
            candidates = np.arange(nodes.shape[0], dtype=np.int64)
            if side == "left":
                candidates = candidates[x[candidates] > x[b] + tol]
                normal = x[candidates] - x[b]
                tangent = y[candidates] - y[b]
            elif side == "right":
                candidates = candidates[x[candidates] < x[b] - tol]
                normal = x[b] - x[candidates]
                tangent = y[candidates] - y[b]
            elif side == "bottom":
                candidates = candidates[y[candidates] > y[b] + tol]
                normal = y[candidates] - y[b]
                tangent = x[candidates] - x[b]
            else:
                candidates = candidates[y[candidates] < y[b] - tol]
                normal = y[b] - y[candidates]
                tangent = x[candidates] - x[b]

        if candidates.size == 0:
            raise ValueError(f"No inward candidates found for {side} node {b}.")

        # Lexicographic intent in one scalar: enforce near-zero tangential
        # displacement first, then choose the nearest inward layer.
        score = (
            (tangent / max(h, 1.0e-300)) ** 2
            + 1.0e-3 * (normal / max(h, 1.0e-300)) ** 2
        )
        inner[k] = int(candidates[int(np.argmin(score))])

    return boundary, inner


def _node_adjacency_from_ops_or_mesh(
    mesh,
    ops: FVOperators | None = None,
) -> dict[int, set[int]] | None:
    """Build node adjacency from FV operators or mesh triangles."""
    n_nodes = int(np.asarray(mesh.nodes).shape[0])
    adjacency: dict[int, set[int]] = {i: set() for i in range(n_nodes)}

    if ops is not None:
        for a, b in zip(np.asarray(ops.edge_i, dtype=np.int64), np.asarray(ops.edge_j, dtype=np.int64)):
            adjacency[int(a)].add(int(b))
            adjacency[int(b)].add(int(a))
        return adjacency

    triangles = getattr(mesh, "triangles", None)
    if triangles is None:
        return None

    tri = np.asarray(triangles, dtype=np.int64)
    if tri.ndim != 2 or tri.shape[1] != 3:
        return None

    for a, b, c in tri:
        adjacency[int(a)].update((int(b), int(c)))
        adjacency[int(b)].update((int(a), int(c)))
        adjacency[int(c)].update((int(a), int(b)))
    return adjacency

