"""Notebook-style edge finite-volume operators for pySNSPD gTDGL.

This module keeps the public OE7 operator API, but changes the core geometry to
match the notebook implementation more closely:

* node control volumes are barycentric triangle areas;
* edge dual lengths are circumcenter/Voronoi lengths ``edge_s``;
* terminal boundary fluxes use lumped boundary-edge lengths;
* edge scalar projections can be reconstructed to node vectors by local LS.

Scalar fields live on nodes. Directed edge scalars are oriented edge_i -> edge_j.
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

    # Optional notebook-style auxiliary data. Tests that instantiate the old
    # dataclass signature remain valid because these fields have defaults.
    triangle_area_m2: np.ndarray | None = None
    triangle_circumcenter_m: np.ndarray | None = None
    edge_triangles: np.ndarray | None = None
    left_nodes: np.ndarray | None = None
    right_nodes: np.ndarray | None = None
    left_boundary_measure_m: np.ndarray | None = None
    right_boundary_measure_m: np.ndarray | None = None
    xi_mesh_m: float | None = None

    @property
    def n_nodes(self) -> int:
        return int(self.node_area_m2.size)

    @property
    def n_edges(self) -> int:
        return int(self.edges.shape[0])


def build_fv_operators(mesh, edge_data) -> FVOperators:
    """Construct notebook-style Delaunay/Voronoi finite-volume factors.

    The previous OE7 implementation used a positive barycentric dual
    ``adjacent_area/(3*l_ij)``. The notebook uses a circumcenter dual length:

    * interior edge: distance between the two adjacent triangle circumcenters;
    * boundary edge: distance from the triangle circumcenter to the edge midpoint.

    The fallback remains positive for degenerate/ill-conditioned edges.
    """
    nodes = np.asarray(mesh.nodes, dtype=float)
    triangles = np.asarray(mesh.triangles, dtype=np.int64)
    edges = np.asarray(edge_data.edges, dtype=np.int64)

    if nodes.ndim != 2 or nodes.shape[1] < 2:
        raise ValueError("mesh.nodes must have shape (n_nodes, >=2).")
    if triangles.ndim != 2 or triangles.shape[1] != 3:
        raise ValueError("mesh.triangles must have shape (n_triangles, 3).")
    if edges.ndim != 2 or edges.shape[1] != 2:
        raise ValueError("edge_data.edges must have shape (n_edges, 2).")

    n_nodes = int(nodes.shape[0])
    tri_area = triangle_areas(nodes, triangles)
    node_area = np.zeros(n_nodes, dtype=float)
    for local in range(3):
        np.add.at(node_area, triangles[:, local], tri_area / 3.0)

    edge_i = edges[:, 0].astype(np.int64, copy=True)
    edge_j = edges[:, 1].astype(np.int64, copy=True)
    edge_vec = nodes[edge_j, :2] - nodes[edge_i, :2]
    edge_length = np.linalg.norm(edge_vec, axis=1)
    if np.any(edge_length <= 0.0):
        raise ValueError("All mesh edges must have positive length.")
    edge_unit = edge_vec / edge_length[:, None]

    circum = triangle_circumcenters(nodes, triangles)
    edge_triangles = _edge_triangles_from_edge_data_or_mesh(edge_data, triangles, edges)
    dual_face = _circumcenter_dual_lengths(
        nodes=nodes,
        edges=edges,
        edge_length=edge_length,
        circumcenters=circum,
        edge_triangles=edge_triangles,
    )

    if np.any(node_area <= 0.0):
        raise ValueError("All mesh nodes must have positive dual area.")

    left_measure_full = boundary_node_measure_m(edge_data, n_nodes=n_nodes, tag="left")
    right_measure_full = boundary_node_measure_m(edge_data, n_nodes=n_nodes, tag="right")
    left_nodes = np.flatnonzero(left_measure_full > 0.0).astype(np.int64)
    right_nodes = np.flatnonzero(right_measure_full > 0.0).astype(np.int64)

    xi_mesh = getattr(mesh, "target_spacing_m", None)
    if xi_mesh is None:
        xi_mesh = float(np.sqrt(np.nanmedian(node_area)))

    return FVOperators(
        edges=edges,
        edge_i=edge_i,
        edge_j=edge_j,
        edge_vec_m=edge_vec,
        edge_length_m=edge_length,
        edge_unit=edge_unit,
        dual_face_length_m=dual_face,
        node_area_m2=node_area,
        triangle_area_m2=tri_area,
        triangle_circumcenter_m=circum,
        edge_triangles=edge_triangles,
        left_nodes=left_nodes,
        right_nodes=right_nodes,
        left_boundary_measure_m=left_measure_full[left_nodes],
        right_boundary_measure_m=right_measure_full[right_nodes],
        xi_mesh_m=float(xi_mesh),
    )


def triangle_areas(nodes: np.ndarray, triangles: np.ndarray) -> np.ndarray:
    """Return triangle areas in square meters."""
    p0 = nodes[triangles[:, 0], :2]
    p1 = nodes[triangles[:, 1], :2]
    p2 = nodes[triangles[:, 2], :2]
    return 0.5 * np.abs(
        (p1[:, 0] - p0[:, 0]) * (p2[:, 1] - p0[:, 1])
        - (p1[:, 1] - p0[:, 1]) * (p2[:, 0] - p0[:, 0])
    )


def triangle_circumcenters(nodes: np.ndarray, triangles: np.ndarray) -> np.ndarray:
    """Return circumcenters for all triangles."""
    pts = np.asarray(nodes, dtype=float)[:, :2]
    circum = np.zeros((triangles.shape[0], 2), dtype=float)
    for k, tri in enumerate(np.asarray(triangles, dtype=np.int64)):
        circum[k] = _circumcenter(pts[tri])
    return circum


def _circumcenter(p: np.ndarray) -> np.ndarray:
    A, B, C = p[0], p[1], p[2]
    D = 2.0 * (
        A[0] * (B[1] - C[1])
        + B[0] * (C[1] - A[1])
        + C[0] * (A[1] - B[1])
    )
    if abs(D) < 1.0e-300:
        return np.mean(p, axis=0)
    a2 = float(np.dot(A, A))
    b2 = float(np.dot(B, B))
    c2 = float(np.dot(C, C))
    ux = (a2 * (B[1] - C[1]) + b2 * (C[1] - A[1]) + c2 * (A[1] - B[1])) / D
    uy = (a2 * (C[0] - B[0]) + b2 * (A[0] - C[0]) + c2 * (B[0] - A[0])) / D
    return np.array([ux, uy], dtype=float)


def _edge_triangles_from_edge_data_or_mesh(edge_data, triangles: np.ndarray, edges: np.ndarray) -> np.ndarray:
    if hasattr(edge_data, "edge_triangles"):
        raw = np.asarray(edge_data.edge_triangles, dtype=np.int64)
        if raw.shape[0] == edges.shape[0]:
            if raw.ndim == 1:
                out = -np.ones((edges.shape[0], 2), dtype=np.int64)
                out[:, 0] = raw
                return out
            if raw.ndim == 2:
                out = -np.ones((edges.shape[0], 2), dtype=np.int64)
                cols = min(2, raw.shape[1])
                out[:, :cols] = raw[:, :cols]
                return out

    edge_lookup = {tuple(sorted(map(int, e))): k for k, e in enumerate(edges)}
    adj: list[list[int]] = [[] for _ in range(edges.shape[0])]
    for t, tri in enumerate(np.asarray(triangles, dtype=np.int64)):
        for u, v in ((tri[0], tri[1]), (tri[1], tri[2]), (tri[2], tri[0])):
            key = tuple(sorted((int(u), int(v))))
            k = edge_lookup.get(key)
            if k is not None:
                adj[k].append(int(t))

    out = -np.ones((edges.shape[0], 2), dtype=np.int64)
    for k, ts in enumerate(adj):
        if ts:
            out[k, : min(2, len(ts))] = ts[:2]
    return out


def _circumcenter_dual_lengths(
    *,
    nodes: np.ndarray,
    edges: np.ndarray,
    edge_length: np.ndarray,
    circumcenters: np.ndarray,
    edge_triangles: np.ndarray,
) -> np.ndarray:
    dual = np.zeros(edges.shape[0], dtype=float)
    pts = np.asarray(nodes, dtype=float)[:, :2]
    for k, (u, v) in enumerate(np.asarray(edges, dtype=np.int64)):
        ts = np.asarray(edge_triangles[k], dtype=np.int64)
        ts = ts[ts >= 0]
        mid = 0.5 * (pts[int(u)] + pts[int(v)])
        if ts.size == 1:
            sij = float(np.linalg.norm(circumcenters[int(ts[0])] - mid))
        elif ts.size >= 2:
            sij = float(np.linalg.norm(circumcenters[int(ts[0])] - circumcenters[int(ts[1])]))
        else:
            sij = 0.5 * float(edge_length[k])
        if not np.isfinite(sij) or sij <= 0.0:
            sij = 0.5 * float(edge_length[k])
        dual[k] = sij
    return dual


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
    """Notebook-style graph FV Laplacian of a node scalar/complex field."""
    v = np.asarray(values)
    out = np.zeros(ops.n_nodes, dtype=np.result_type(v, np.complex128))
    flux = ops.dual_face_length_m * (v[ops.edge_j] - v[ops.edge_i]) / ops.edge_length_m
    np.add.at(out, ops.edge_i, flux)
    np.add.at(out, ops.edge_j, -flux)
    return out / ops.node_area_m2


def edge_flux_accumulator_A_m(edge_current_i_to_j_A_m2: np.ndarray, ops: FVOperators) -> np.ndarray:
    """Return conservative edge-flux accumulator before division by area."""
    current = np.asarray(edge_current_i_to_j_A_m2, dtype=float)
    if current.shape != (ops.n_edges,):
        raise ValueError(f"edge current must have shape ({ops.n_edges},), got {current.shape}.")
    flux = ops.dual_face_length_m * current
    out = np.zeros(ops.n_nodes, dtype=float)
    np.add.at(out, ops.edge_i, flux)
    np.add.at(out, ops.edge_j, -flux)
    return out


def divergence_from_edge_scalar(
    edge_current_i_to_j: np.ndarray,
    ops: FVOperators,
    *,
    boundary_accum_A_m: np.ndarray | None = None,
) -> np.ndarray:
    """Finite-volume divergence of an oriented edge-current density."""
    out = edge_flux_accumulator_A_m(edge_current_i_to_j, ops)
    if boundary_accum_A_m is not None:
        boundary = np.asarray(boundary_accum_A_m, dtype=float)
        if boundary.shape != out.shape:
            raise ValueError(f"boundary_accum_A_m must have shape {out.shape}, got {boundary.shape}.")
        out += boundary
    return out / ops.node_area_m2


def boundary_node_measure_m(edge_data, *, n_nodes: int, tag: str) -> np.ndarray:
    """Lump tagged boundary-edge lengths onto boundary nodes."""
    tags = np.asarray(edge_data.tags).astype(str)
    edges = np.asarray(edge_data.edges, dtype=np.int64)
    if hasattr(edge_data, "lengths"):
        lengths = np.asarray(edge_data.lengths, dtype=float)
    else:
        raise ValueError("edge_data.lengths is required for boundary measures.")
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
    """Return prescribed outward current densities at left/right terminals."""
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
    """Build nodal terminal-flux accumulator for div(j)=0."""
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
    """Unwrap phase by walking the mesh graph."""
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
    """Simple edge-to-node vector average, retained for backward compatibility."""
    current = np.asarray(edge_current_i_to_j, dtype=float)
    vec = current[:, None] * ops.edge_unit
    cnt = np.zeros(ops.n_nodes, dtype=float)
    out = np.zeros((ops.n_nodes, 2), dtype=float)
    for comp in range(2):
        np.add.at(out[:, comp], ops.edge_i, vec[:, comp])
        np.add.at(out[:, comp], ops.edge_j, vec[:, comp])
    np.add.at(cnt, ops.edge_i, 1.0)
    np.add.at(cnt, ops.edge_j, 1.0)
    out[:, 0] /= np.maximum(cnt, 1.0)
    out[:, 1] /= np.maximum(cnt, 1.0)
    return out[:, 0], out[:, 1]


def edge_scalar_to_node_vector_least_squares(
    edge_current_i_to_j: np.ndarray,
    ops: FVOperators,
    *,
    ridge: float = 1.0e-30,
) -> tuple[np.ndarray, np.ndarray]:
    """Notebook-style local LS reconstruction from edge projections.

    The notebook uses weights ``edge_s/edge_len``. This is important for
    unstructured Delaunay meshes; using raw ``edge_s`` changes the local LS
    balance and can distort the reconstructed magnitude.
    """
    vals = np.asarray(edge_current_i_to_j, dtype=float)
    if vals.shape != (ops.n_edges,):
        raise ValueError(f"edge_current_i_to_j must have shape ({ops.n_edges},).")
    ex = ops.edge_unit[:, 0]
    ey = ops.edge_unit[:, 1]
    w = ops.dual_face_length_m / np.maximum(ops.edge_length_m, 1.0e-300)
    w = np.maximum(w, 1.0e-300)

    Axx = np.zeros(ops.n_nodes, dtype=float)
    Axy = np.zeros(ops.n_nodes, dtype=float)
    Ayy = np.zeros(ops.n_nodes, dtype=float)
    bx = np.zeros(ops.n_nodes, dtype=float)
    by = np.zeros(ops.n_nodes, dtype=float)

    cxx = w * ex * ex
    cxy = w * ex * ey
    cyy = w * ey * ey
    rx = w * vals * ex
    ry = w * vals * ey
    for nodes in (ops.edge_i, ops.edge_j):
        np.add.at(Axx, nodes, cxx)
        np.add.at(Axy, nodes, cxy)
        np.add.at(Ayy, nodes, cyy)
        np.add.at(bx, nodes, rx)
        np.add.at(by, nodes, ry)

    trace = Axx + Ayy
    reg = ridge * np.maximum(trace, 1.0)
    Axx = Axx + reg
    Ayy = Ayy + reg
    det = Axx * Ayy - Axy * Axy
    good = np.abs(det) > 1.0e-300

    vx = np.zeros(ops.n_nodes, dtype=float)
    vy = np.zeros(ops.n_nodes, dtype=float)
    vx[good] = (Ayy[good] * bx[good] - Axy[good] * by[good]) / det[good]
    vy[good] = (-Axy[good] * bx[good] + Axx[good] * by[good]) / det[good]

    if not np.all(good):
        fx, fy = edge_scalar_to_node_vector(vals, ops)
        vx[~good] = fx[~good]
        vy[~good] = fy[~good]
    return vx, vy


def boundary_currents_from_edge_scalar_least_squares(
    *,
    mesh,
    edge_data,
    ops: FVOperators,
    edge_current_i_to_j: np.ndarray,
    thickness_m: float,
) -> dict[str, float]:
    """Boundary-current diagnostic using LS node-vector reconstruction."""
    jx, jy = edge_scalar_to_node_vector_least_squares(edge_current_i_to_j, ops)
    return boundary_currents_from_node_vectors(
        mesh=mesh,
        edge_data=edge_data,
        jx_A_m2=jx,
        jy_A_m2=jy,
        thickness_m=thickness_m,
    )


def strip_transport_current_profile_from_node_vectors(
    *,
    mesh,
    jx_A_m2: np.ndarray,
    thickness_m: float,
    n_bins: int = 41,
) -> tuple[np.ndarray, np.ndarray]:
    """Estimate I(x)=d*w*<jx>_strip as a diagnostic profile."""
    nodes = np.asarray(mesh.nodes, dtype=float)
    jx = np.asarray(jx_A_m2, dtype=float)
    x = nodes[:, 0]
    xmin = float(np.min(x))
    xmax = float(np.max(x))
    if n_bins < 2:
        raise ValueError("n_bins must be at least 2.")
    edges_x = np.linspace(xmin, xmax, int(n_bins) + 1)
    centers = 0.5 * (edges_x[:-1] + edges_x[1:])
    current_A = np.full(int(n_bins), np.nan, dtype=float)
    for k in range(int(n_bins)):
        if k == int(n_bins) - 1:
            mask = (x >= edges_x[k]) & (x <= edges_x[k + 1])
        else:
            mask = (x >= edges_x[k]) & (x < edges_x[k + 1])
        if np.any(mask):
            current_A[k] = float(thickness_m * mesh.width_m * np.mean(jx[mask]))
    return centers, current_A


def terminal_voltage(nodes: np.ndarray, phi_V: np.ndarray, *, length_m: float) -> float:
    """Return mean(phi_right)-mean(phi_left)."""
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
