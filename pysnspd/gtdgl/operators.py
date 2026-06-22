"""Edge-based finite-volume operators for pySNSPD gTDGL.

Scalar fields live on nodes. Phase gradients and edge currents live on mesh
edges. Divergences are accumulated back to nodes in conservative form.
"""
from __future__ import annotations

from dataclasses import dataclass
import numpy as np


@dataclass(frozen=True)
class FVOperators:
    """Geometry factors for edge-based finite-volume operations."""

    edges: np.ndarray
    edge_i: np.ndarray
    edge_j: np.ndarray
    edge_vec_m: np.ndarray
    edge_length_m: np.ndarray
    edge_unit: np.ndarray
    dual_face_length_m: np.ndarray
    node_area_m2: np.ndarray

    @property
    def n_nodes(self) -> int:
        return int(self.node_area_m2.size)

    @property
    def n_edges(self) -> int:
        return int(self.edges.shape[0])


def build_fv_operators(mesh, edge_data) -> FVOperators:
    """Construct positive barycentric-dual FV geometry factors."""
    nodes = np.asarray(mesh.nodes, dtype=float)
    triangles = np.asarray(mesh.triangles, dtype=np.int64)
    edges = np.asarray(edge_data.edges, dtype=np.int64)
    edge_triangles = np.asarray(edge_data.edge_triangles, dtype=np.int64)

    n_nodes = int(nodes.shape[0])
    tri_area = triangle_areas(nodes, triangles)

    node_area = np.zeros(n_nodes, dtype=float)
    for local in range(3):
        np.add.at(node_area, triangles[:, local], tri_area / 3.0)

    edge_i = edges[:, 0].astype(np.int64, copy=True)
    edge_j = edges[:, 1].astype(np.int64, copy=True)
    edge_vec = nodes[edge_j] - nodes[edge_i]
    edge_length = np.linalg.norm(edge_vec, axis=1)
    if np.any(edge_length <= 0.0):
        raise ValueError("All mesh edges must have positive length.")

    edge_unit = edge_vec / edge_length[:, None]

    adjacent_area = np.zeros(edges.shape[0], dtype=float)
    for k in range(edges.shape[0]):
        tri_ids = edge_triangles[k]
        valid = tri_ids[tri_ids >= 0]
        adjacent_area[k] = float(np.sum(tri_area[valid]))

    dual_face = adjacent_area / (3.0 * edge_length)
    positive = dual_face > 0.0
    if not np.all(positive):
        replacement = float(np.min(dual_face[positive])) if np.any(positive) else 1.0
        dual_face = np.where(positive, dual_face, replacement)

    if np.any(node_area <= 0.0):
        raise ValueError("All mesh nodes must have positive dual area.")

    return FVOperators(
        edges=edges,
        edge_i=edge_i,
        edge_j=edge_j,
        edge_vec_m=edge_vec,
        edge_length_m=edge_length,
        edge_unit=edge_unit,
        dual_face_length_m=dual_face,
        node_area_m2=node_area,
    )


def triangle_areas(nodes: np.ndarray, triangles: np.ndarray) -> np.ndarray:
    """Return triangle areas in square meters."""
    p0 = nodes[triangles[:, 0]]
    p1 = nodes[triangles[:, 1]]
    p2 = nodes[triangles[:, 2]]
    return 0.5 * np.abs(
        (p1[:, 0] - p0[:, 0]) * (p2[:, 1] - p0[:, 1])
        - (p1[:, 1] - p0[:, 1]) * (p2[:, 0] - p0[:, 0])
    )


def edge_average(values: np.ndarray, ops: FVOperators) -> np.ndarray:
    """Arithmetic average of a node scalar on each edge."""
    v = np.asarray(values)
    return 0.5 * (v[ops.edge_i] + v[ops.edge_j])


def edge_scalar_gradient(values: np.ndarray, ops: FVOperators) -> np.ndarray:
    """Directional scalar gradient from edge_i to edge_j."""
    v = np.asarray(values)
    return (v[ops.edge_j] - v[ops.edge_i]) / ops.edge_length_m


def edge_phase_gradient_from_psi(psi: np.ndarray, ops: FVOperators) -> np.ndarray:
    """Gauge-safe phase gradient from the complex order parameter."""
    z = np.asarray(psi, dtype=np.complex128)
    dtheta = np.angle(z[ops.edge_j] * np.conjugate(z[ops.edge_i]))
    return dtheta / ops.edge_length_m


def laplacian(values: np.ndarray, ops: FVOperators) -> np.ndarray:
    """Graph finite-volume Laplacian of a node scalar or complex field."""
    v = np.asarray(values)
    out = np.zeros(ops.n_nodes, dtype=np.result_type(v, np.complex128))
    flux = ops.dual_face_length_m * (v[ops.edge_j] - v[ops.edge_i]) / ops.edge_length_m
    np.add.at(out, ops.edge_i, flux)
    np.add.at(out, ops.edge_j, -flux)
    return out / ops.node_area_m2


def divergence_from_edge_scalar(
    edge_current_i_to_j: np.ndarray,
    ops: FVOperators,
    *,
    boundary_accum_A_m: np.ndarray | None = None,
) -> np.ndarray:
    """Finite-volume divergence of an edge current scalar.

    The edge scalar is positive when current flows from edge_i to edge_j.

    Parameters
    ----------
    edge_current_i_to_j:
        Edge current density [A/m^2] projected along the oriented edge.
    ops:
        Finite-volume edge geometry.
    boundary_accum_A_m:
        Optional nodal boundary flux accumulator [A/m]. This is where
        prescribed terminal currents enter the conservative divergence.
    """
    current = np.asarray(edge_current_i_to_j, dtype=float)
    out = np.zeros(ops.n_nodes, dtype=float)

    flux = ops.dual_face_length_m * current
    np.add.at(out, ops.edge_i, flux)
    np.add.at(out, ops.edge_j, -flux)

    if boundary_accum_A_m is not None:
        boundary = np.asarray(boundary_accum_A_m, dtype=float)
        if boundary.shape != out.shape:
            raise ValueError(
                "boundary_accum_A_m must have shape "
                f"{out.shape}, got {boundary.shape}."
            )
        out += boundary

    return out / ops.node_area_m2


def boundary_node_measure_m(
    edge_data,
    *,
    n_nodes: int,
    tag: str,
) -> np.ndarray:
    """Lump boundary-edge lengths onto boundary nodes.

    The returned array has units of meters and sums to the geometric length of
    the tagged boundary. Each boundary edge contributes half of its length to
    each endpoint.
    """
    tags = np.asarray(edge_data.tags).astype(str)
    edges = np.asarray(edge_data.edges, dtype=np.int64)
    lengths = np.asarray(edge_data.lengths, dtype=float)

    measure = np.zeros(int(n_nodes), dtype=float)
    mask = tags == str(tag)
    if not np.any(mask):
        return measure

    e = edges[mask]
    ell = lengths[mask]
    np.add.at(measure, e[:, 0], 0.5 * ell)
    np.add.at(measure, e[:, 1], 0.5 * ell)
    return measure


def terminal_outward_flux_densities_A_m2(
    edge_data,
    *,
    n_nodes: int,
    target_current_A: float,
    thickness_m: float,
) -> dict[str, float]:
    """Return prescribed outward current densities at left/right terminals.

    Positive current flows from left to right. Therefore the outward flux is
    negative at the left terminal and positive at the right terminal.
    """
    left_measure = boundary_node_measure_m(edge_data, n_nodes=n_nodes, tag="left")
    right_measure = boundary_node_measure_m(edge_data, n_nodes=n_nodes, tag="right")

    width_left_m = float(np.sum(left_measure))
    width_right_m = float(np.sum(right_measure))
    if width_left_m <= 0.0 or width_right_m <= 0.0:
        raise ValueError("Left/right terminal boundary measures must be positive.")
    if thickness_m <= 0.0:
        raise ValueError("thickness_m must be positive.")

    I = float(target_current_A)
    return {
        "left_A_m2": -I / (float(thickness_m) * width_left_m),
        "right_A_m2": I / (float(thickness_m) * width_right_m),
        "width_left_m": width_left_m,
        "width_right_m": width_right_m,
    }


def terminal_boundary_accum_A_m(
    edge_data,
    *,
    n_nodes: int,
    target_current_A: float,
    thickness_m: float,
) -> np.ndarray:
    """Build nodal terminal-flux accumulator for div(j)=0.

    The accumulator has units [A/m], i.e. current per film thickness. It must be
    added to the edge-flux accumulator before division by node area.
    """
    left_measure = boundary_node_measure_m(edge_data, n_nodes=n_nodes, tag="left")
    right_measure = boundary_node_measure_m(edge_data, n_nodes=n_nodes, tag="right")
    flux = terminal_outward_flux_densities_A_m2(
        edge_data,
        n_nodes=n_nodes,
        target_current_A=target_current_A,
        thickness_m=thickness_m,
    )

    out = np.zeros(int(n_nodes), dtype=float)
    out += left_measure * flux["left_A_m2"]
    out += right_measure * flux["right_A_m2"]
    return out


def unwrap_phase_graph(
    psi: np.ndarray,
    edges: np.ndarray,
    *,
    seed_index: int | None = None,
    subtract_mean: bool = False,
) -> np.ndarray:
    """Unwrap phase by walking the mesh graph.

    This avoids the artificial 2*pi accumulation that appears when np.unwrap is
    applied to a flattened 2D mesh.
    """
    z = np.asarray(psi, dtype=np.complex128).reshape(-1)
    edges = np.asarray(edges, dtype=np.int64)
    n = int(z.size)

    if n == 0:
        return np.array([], dtype=float)

    if seed_index is None:
        seed_index = 0
    seed_index = int(seed_index)

    adj: list[list[tuple[int, float]]] = [[] for _ in range(n)]
    dtheta = np.angle(z[edges[:, 1]] * np.conjugate(z[edges[:, 0]]))
    for (i, j), dth in zip(edges, dtheta):
        ii = int(i)
        jj = int(j)
        adj[ii].append((jj, float(dth)))
        adj[jj].append((ii, -float(dth)))

    theta = np.full(n, np.nan, dtype=float)
    visited = np.zeros(n, dtype=bool)

    starts = [seed_index] + [i for i in range(n) if i != seed_index]
    for start in starts:
        if visited[start]:
            continue
        theta[start] = float(np.angle(z[start]))
        visited[start] = True
        stack = [start]

        while stack:
            i = stack.pop()
            for j, dth in adj[i]:
                if not visited[j]:
                    theta[j] = theta[i] + dth
                    visited[j] = True
                    stack.append(j)

    if subtract_mean:
        theta -= float(np.nanmean(theta))

    return theta

def edge_scalar_to_node_vector(edge_current_i_to_j: np.ndarray, ops: FVOperators) -> tuple[np.ndarray, np.ndarray]:
    """Average edge-oriented current scalars into node vector components."""
    current = np.asarray(edge_current_i_to_j, dtype=float)
    vec = current[:, None] * ops.edge_unit
    wx = ops.dual_face_length_m
    sum_w = np.zeros(ops.n_nodes, dtype=float)
    out = np.zeros((ops.n_nodes, 2), dtype=float)

    for comp in range(2):
        contrib = wx * vec[:, comp]
        np.add.at(out[:, comp], ops.edge_i, contrib)
        np.add.at(out[:, comp], ops.edge_j, contrib)
    np.add.at(sum_w, ops.edge_i, wx)
    np.add.at(sum_w, ops.edge_j, wx)

    safe = np.maximum(sum_w, 1.0e-300)
    out[:, 0] /= safe
    out[:, 1] /= safe
    return out[:, 0], out[:, 1]


def terminal_voltage(nodes: np.ndarray, phi_V: np.ndarray, *, length_m: float) -> float:
    """Return <phi>_right - <phi>_left."""
    nodes = np.asarray(nodes, dtype=float)
    phi = np.asarray(phi_V, dtype=float)
    x = nodes[:, 0]
    tol = max(1.0e-15, 1.0e-9 * float(length_m))
    left = np.abs(x - np.min(x)) <= tol
    right = np.abs(x - np.max(x)) <= tol
    if not np.any(left) or not np.any(right):
        return float("nan")
    return float(np.mean(phi[right]) - np.mean(phi[left]))


def boundary_currents_from_node_vectors(
    *,
    mesh,
    edge_data,
    jx_A_m2: np.ndarray,
    jy_A_m2: np.ndarray,
    thickness_m: float,
) -> dict[str, float]:
    """Integrate node-vector current through tagged geometric boundaries."""
    tags = np.asarray(edge_data.tags).astype(str)
    edges = np.asarray(edge_data.edges, dtype=np.int64)
    lengths = np.asarray(edge_data.lengths, dtype=float)

    jx = np.asarray(jx_A_m2, dtype=float)
    jy = np.asarray(jy_A_m2, dtype=float)

    normals = {
        "left": (-1.0, 0.0),
        "right": (1.0, 0.0),
        "bottom": (0.0, -1.0),
        "top": (0.0, 1.0),
    }

    out: dict[str, float] = {}
    for tag, (nx, ny) in normals.items():
        mask = tags == tag
        if not np.any(mask):
            out[f"{tag}_A"] = 0.0
            continue
        e = edges[mask]
        jx_edge = 0.5 * (jx[e[:, 0]] + jx[e[:, 1]])
        jy_edge = 0.5 * (jy[e[:, 0]] + jy[e[:, 1]])
        flux_density = jx_edge * nx + jy_edge * ny
        out[f"{tag}_A"] = float(thickness_m * np.sum(lengths[mask] * flux_density))

    out["net_boundary_current_A"] = float(sum(out.values()))
    return out