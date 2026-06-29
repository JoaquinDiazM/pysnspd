"""pyTDGL-like rectangular Delaunay mesh generation for pySNSPD.

This module intentionally follows the finite-volume mesh idea used by
``pyTDGL``: a triangular Delaunay mesh whose vertices carry scalar fields and
whose dual Voronoi cells define the finite-volume control volumes.  pySNSPD
continues to store coordinates in SI meters and to use its existing PRE-run
file formats.

The generator below replaces the previous protected structured grid for the
pyTDGL-like branch.  It keeps the rectangle and current terminals explicit, but
uses an unstructured Delaunay point cloud plus optional Laplacian smoothing of
interior sites, analogous to ``tdgl.finite_volume.Mesh.smooth``.
"""
from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Any, Mapping

import numpy as np
from scipy.spatial import Delaunay

from pysnspd.mesh.delaunay import MeshData, orient_triangles_counterclockwise


@dataclass(frozen=True)
class PyTDGLLikeMeshParameters:
    """Parameters used by the pyTDGL-like rectangular mesh generator."""

    length_m: float
    width_m: float
    target_spacing_m: float
    seed: int
    boundary_spacing_factor: float = 1.0
    interior_jitter_fraction: float = 0.20
    smoothing_iterations: int = 2
    terminal_contact_mode: str = "normal_left_right"


def parameters_from_config(
    config: Mapping[str, Any],
    *,
    jitter_fraction: float = 0.20,
    boundary_guard_layers: int = 1,
    boundary_spacing_factor: float = 1.0,
    smoothing_iterations: int | None = None,
) -> PyTDGLLikeMeshParameters:
    """Resolve pyTDGL-like mesh parameters from the pySNSPD config.

    ``boundary_guard_layers`` is accepted for compatibility with the old
    PRE-run CLI.  In this unstructured backend it controls the default number
    of Laplacian smoothing iterations when ``smoothing_iterations`` is not
    supplied: more protected layers in the old mesh correspond to a smoother
    interior mesh here.
    """

    # Mesh generation must not validate the full pySNSPD runtime config.
    # Unit tests and small mesh experiments often provide only the sections
    # needed here.  The full PRE-run pipeline is still responsible for loading
    # and validating the complete config before it reaches this helper.
    cfg = dict(config)
    if "material" not in cfg:
        raise KeyError("config must contain a 'material' section")
    if "mesh" not in cfg:
        raise KeyError("config must contain a 'mesh' section")

    material = cfg["material"]
    mesh_cfg = cfg["mesh"]

    width_m = float(material["width_m"])
    spacing_m = float(mesh_cfg["target_spacing_m"])
    seed = int(mesh_cfg.get("seed", 0))
    if "length_m" in mesh_cfg:
        length_m = float(mesh_cfg["length_m"])
    elif "geometry" in cfg and isinstance(cfg["geometry"], Mapping) and "length_m" in cfg["geometry"]:
        length_m = float(cfg["geometry"]["length_m"])
    else:
        length_m = 2.0 * width_m

    if smoothing_iterations is None:
        smoothing_iterations = max(0, int(boundary_guard_layers))

    return PyTDGLLikeMeshParameters(
        length_m=length_m,
        width_m=width_m,
        target_spacing_m=spacing_m,
        seed=seed,
        boundary_spacing_factor=float(boundary_spacing_factor),
        interior_jitter_fraction=float(jitter_fraction),
        smoothing_iterations=int(smoothing_iterations),
    )


def generate_rectangular_pytdgl_like_mesh(
    config: Mapping[str, Any],
    *,
    jitter_fraction: float = 0.20,
    boundary_guard_layers: int = 1,
    boundary_spacing_factor: float = 1.0,
    smoothing_iterations: int | None = None,
) -> MeshData:
    """Generate an unstructured pyTDGL-like Delaunay mesh in SI meters.

    The rectangle is ``x in [0, L]`` and ``y in [-W/2, W/2]``.  Boundary
    vertices are placed exactly on the four physical boundaries so the normal
    contacts at ``left`` and ``right`` remain explicit.  Interior vertices are
    generated from a slightly jittered triangular lattice and triangulated with
    SciPy/Qhull Delaunay.
    """

    params = parameters_from_config(
        config,
        jitter_fraction=jitter_fraction,
        boundary_guard_layers=boundary_guard_layers,
        boundary_spacing_factor=boundary_spacing_factor,
        smoothing_iterations=smoothing_iterations,
    )
    return generate_rectangular_pytdgl_like_mesh_from_parameters(params)


def generate_rectangular_pytdgl_like_mesh_from_parameters(
    params: PyTDGLLikeMeshParameters,
) -> MeshData:
    """Generate a pyTDGL-like mesh from resolved parameters."""

    _validate_parameters(params)
    rng = np.random.default_rng(int(params.seed))

    boundary = _rectangular_boundary_points(
        params.length_m,
        params.width_m,
        params.target_spacing_m * params.boundary_spacing_factor,
    )
    interior = _triangular_lattice_interior_points(
        params.length_m,
        params.width_m,
        params.target_spacing_m,
        jitter_fraction=params.interior_jitter_fraction,
        rng=rng,
    )
    nodes = np.vstack([boundary, interior]) if interior.size else boundary.copy()
    nodes = _unique_rows_with_tolerance(nodes, tol=1.0e-15)
    triangles = _delaunay_triangles_inside_rectangle(nodes, params.length_m, params.width_m)
    nodes, triangles = _compact_used_nodes(nodes, triangles)
    triangles = orient_triangles_counterclockwise(nodes, triangles)

    boundary_mask = _boundary_node_mask(nodes, params.length_m, params.width_m)
    if params.smoothing_iterations > 0:
        nodes = laplacian_smooth_interior_nodes(
            nodes,
            triangles,
            boundary_mask=boundary_mask,
            iterations=params.smoothing_iterations,
        )
        nodes = _snap_rectangular_boundary(nodes, params.length_m, params.width_m)
        triangles = orient_triangles_counterclockwise(nodes, triangles)

    return MeshData(
        nodes=np.asarray(nodes, dtype=float),
        triangles=np.asarray(triangles, dtype=np.int64),
        length_m=float(params.length_m),
        width_m=float(params.width_m),
        target_spacing_m=float(params.target_spacing_m),
        seed=int(params.seed),
        triangulation_method="pytdgl_like_unstructured_delaunay_v1",
        boundary_guard_layers=int(params.smoothing_iterations),
    )


def laplacian_smooth_interior_nodes(
    nodes: np.ndarray,
    triangles: np.ndarray,
    *,
    boundary_mask: np.ndarray,
    iterations: int,
) -> np.ndarray:
    """Laplacian-smooth interior nodes while holding boundary nodes fixed.

    This mirrors pyTDGL's ``Mesh.smooth`` behavior: each interior site moves to
    the arithmetic average of its neighboring sites, then boundary sites are
    restored exactly.
    """

    out = np.asarray(nodes, dtype=float).copy()
    boundary_mask = np.asarray(boundary_mask, dtype=bool)
    edges = _unique_triangle_edges(np.asarray(triangles, dtype=np.int64))
    n = out.shape[0]
    boundary_values = out[boundary_mask].copy()

    for _ in range(max(0, int(iterations))):
        num_neighbors = np.bincount(edges.ravel(), minlength=n).astype(float)
        new_sites = np.zeros_like(out)
        for dim in (0, 1):
            vals = out[:, dim]
            accum = np.bincount(edges[:, 0], weights=vals[edges[:, 1]], minlength=n)
            accum += np.bincount(edges[:, 1], weights=vals[edges[:, 0]], minlength=n)
            new_sites[:, dim] = accum / np.maximum(num_neighbors, 1.0)
        out[~boundary_mask] = new_sites[~boundary_mask]
        out[boundary_mask] = boundary_values
    return out


def _validate_parameters(params: PyTDGLLikeMeshParameters) -> None:
    if params.length_m <= 0.0:
        raise ValueError("length_m must be positive.")
    if params.width_m <= 0.0:
        raise ValueError("width_m must be positive.")
    if params.target_spacing_m <= 0.0:
        raise ValueError("target_spacing_m must be positive.")
    if params.boundary_spacing_factor <= 0.0:
        raise ValueError("boundary_spacing_factor must be positive.")
    if not (0.0 <= params.interior_jitter_fraction < 0.5):
        raise ValueError("interior_jitter_fraction must satisfy 0 <= f < 0.5.")


def _rectangular_boundary_points(length_m: float, width_m: float, spacing_m: float) -> np.ndarray:
    half_w = 0.5 * width_m
    nx = max(2, int(math.ceil(length_m / spacing_m)))
    ny = max(2, int(math.ceil(width_m / spacing_m)))

    xs = np.linspace(0.0, length_m, nx + 1)
    ys = np.linspace(-half_w, half_w, ny + 1)

    bottom = np.column_stack([xs, np.full_like(xs, -half_w)])
    right = np.column_stack([np.full_like(ys[1:-1], length_m), ys[1:-1]])
    top = np.column_stack([xs[::-1], np.full_like(xs, half_w)])
    left = np.column_stack([np.full_like(ys[-2:0:-1], 0.0), ys[-2:0:-1]])
    return np.vstack([bottom, right, top, left])


def _triangular_lattice_interior_points(
    length_m: float,
    width_m: float,
    spacing_m: float,
    *,
    jitter_fraction: float,
    rng: np.random.Generator,
) -> np.ndarray:
    half_w = 0.5 * width_m
    dy = 0.5 * math.sqrt(3.0) * spacing_m
    xs_base = np.arange(spacing_m, length_m, spacing_m)
    ys = np.arange(-half_w + dy, half_w, dy)
    pts: list[np.ndarray] = []
    margin = 0.35 * spacing_m
    jitter = jitter_fraction * spacing_m

    for row, y in enumerate(ys):
        offset = 0.5 * spacing_m if row % 2 else 0.0
        xs = xs_base + offset
        xs = xs[(xs > margin) & (xs < length_m - margin)]
        if xs.size == 0:
            continue
        row_pts = np.column_stack([xs, np.full_like(xs, y)])
        if jitter > 0:
            row_pts += rng.uniform(-jitter, jitter, size=row_pts.shape)
        keep = (
            (row_pts[:, 0] > margin)
            & (row_pts[:, 0] < length_m - margin)
            & (row_pts[:, 1] > -half_w + margin)
            & (row_pts[:, 1] < half_w - margin)
        )
        if np.any(keep):
            pts.append(row_pts[keep])
    if not pts:
        return np.empty((0, 2), dtype=float)
    return np.vstack(pts)


def _delaunay_triangles_inside_rectangle(nodes: np.ndarray, length_m: float, width_m: float) -> np.ndarray:
    tri = Delaunay(nodes)
    triangles = np.asarray(tri.simplices, dtype=np.int64)
    centroids = nodes[triangles].mean(axis=1)
    half_w = 0.5 * width_m
    keep = (
        (centroids[:, 0] >= -1.0e-15)
        & (centroids[:, 0] <= length_m + 1.0e-15)
        & (centroids[:, 1] >= -half_w - 1.0e-15)
        & (centroids[:, 1] <= half_w + 1.0e-15)
    )
    triangles = triangles[keep]
    if triangles.size == 0:
        raise RuntimeError("Delaunay triangulation produced no in-domain triangles.")
    return triangles


def _compact_used_nodes(nodes: np.ndarray, triangles: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    used = np.unique(triangles.ravel())
    remap = -np.ones(nodes.shape[0], dtype=np.int64)
    remap[used] = np.arange(used.size, dtype=np.int64)
    return nodes[used], remap[triangles]


def _boundary_node_mask(nodes: np.ndarray, length_m: float, width_m: float) -> np.ndarray:
    half_w = 0.5 * width_m
    tol = max(1.0e-15, 1.0e-8 * min(length_m, width_m))
    return (
        np.isclose(nodes[:, 0], 0.0, atol=tol)
        | np.isclose(nodes[:, 0], length_m, atol=tol)
        | np.isclose(nodes[:, 1], -half_w, atol=tol)
        | np.isclose(nodes[:, 1], half_w, atol=tol)
    )


def _snap_rectangular_boundary(nodes: np.ndarray, length_m: float, width_m: float) -> np.ndarray:
    out = np.asarray(nodes, dtype=float).copy()
    half_w = 0.5 * width_m
    tol = max(1.0e-15, 1.0e-8 * min(length_m, width_m))
    out[np.isclose(out[:, 0], 0.0, atol=tol), 0] = 0.0
    out[np.isclose(out[:, 0], length_m, atol=tol), 0] = length_m
    out[np.isclose(out[:, 1], -half_w, atol=tol), 1] = -half_w
    out[np.isclose(out[:, 1], half_w, atol=tol), 1] = half_w
    return out


def _unique_triangle_edges(triangles: np.ndarray) -> np.ndarray:
    edges = np.vstack(
        [
            triangles[:, [0, 1]],
            triangles[:, [1, 2]],
            triangles[:, [2, 0]],
        ]
    )
    edges = np.sort(edges, axis=1)
    return np.unique(edges, axis=0)


def _unique_rows_with_tolerance(points: np.ndarray, *, tol: float) -> np.ndarray:
    if points.size == 0:
        return points.reshape(0, 2)
    scaled = np.round(points / tol).astype(np.int64)
    _, idx = np.unique(scaled, axis=0, return_index=True)
    return points[np.sort(idx)]
