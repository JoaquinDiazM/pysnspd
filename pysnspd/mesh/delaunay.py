"""
Protected rectangular triangulation for pySNSPD.

OE2 scope:
- Build a reproducible rectangular nanowire mesh.
- Preserve all boundary nodes.
- Avoid Qhull/Delaunay coplanar-node loss at straight boundaries.
- Keep the geometry simple and explicit.

The public function keeps the historical name
``generate_rectangular_delaunay_mesh`` because later modules already call it,
but the implementation now uses a protected rectangular grid with a local
diagonal choice that maximizes triangle quality inside each cell.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any, Mapping
import math

import numpy as np

from pysnspd.config import validate_config


@dataclass(frozen=True)
class MeshData:
    """
    Container for a 2D nanowire mesh.

    Coordinates are in meters. The rectangle is represented as

        x in [0, length_m],
        y in [-width_m/2, width_m/2].
    """
    nodes: np.ndarray
    triangles: np.ndarray
    length_m: float
    width_m: float
    target_spacing_m: float
    seed: int
    triangulation_method: str = "protected_structured_local_delaunay"
    boundary_guard_layers: int = 1

    @property
    def n_nodes(self) -> int:
        """Return the number of mesh nodes."""
        return int(self.nodes.shape[0])

    @property
    def n_triangles(self) -> int:
        """Return the number of triangles."""
        return int(self.triangles.shape[0])

    @property
    def extent_m(self) -> tuple[float, float, float, float]:
        """Return mesh extent as ``(xmin, xmax, ymin, ymax)``."""
        return (0.0, self.length_m, -0.5 * self.width_m, 0.5 * self.width_m)


def geometry_from_config(config: Mapping[str, Any]) -> dict[str, float]:
    """
    Resolve the rectangular nanowire geometry from the config.

    Priority for length:

    1. ``mesh.length_m`` if present.
    2. ``geometry.length_m`` if present.
    3. fallback ``2 * material.width_m``.

    The fallback follows the OE2 rule: keep the mesh short enough to avoid
    wasting nodes in longitudinal regions that are not yet dynamically useful.
    """
    cfg = validate_config(config, require_big_data_root_exists=False)

    width_m = float(cfg["material"]["width_m"])
    spacing_m = float(cfg["mesh"]["target_spacing_m"])
    seed = int(cfg["mesh"]["seed"])

    if "length_m" in cfg["mesh"]:
        length_m = float(cfg["mesh"]["length_m"])
    elif "geometry" in cfg and isinstance(cfg["geometry"], Mapping) and "length_m" in cfg["geometry"]:
        length_m = float(cfg["geometry"]["length_m"])
    else:
        length_m = 2.0 * width_m

    if length_m <= 0.0:
        raise ValueError("Nanowire length_m must be positive.")
    if width_m <= 0.0:
        raise ValueError("Nanowire width_m must be positive.")
    if spacing_m <= 0.0:
        raise ValueError("mesh.target_spacing_m must be positive.")

    return {
        "length_m": length_m,
        "width_m": width_m,
        "target_spacing_m": spacing_m,
        "seed": seed,
    }


def generate_rectangular_delaunay_mesh(
    config: Mapping[str, Any],
    *,
    jitter_fraction: float = 0.20,
    boundary_guard_layers: int = 1,
) -> MeshData:
    """
    Generate a protected rectangular triangulation for a nanowire.

    The old free-Delaunay approach can drop collinear boundary nodes. That is
    dangerous for later current-conservation and boundary-condition operators.
    This implementation instead starts from a rectangular point lattice, jitters
    only sufficiently interior nodes, and triangulates each cell using the
    diagonal that maximizes the minimum angle of the two local triangles.

    Parameters
    ----------
    config:
        Valid pySNSPD configuration dictionary.
    jitter_fraction:
        Interior random displacement as a fraction of the nominal spacing.
        Must satisfy ``0 <= jitter_fraction < 0.5``.
    boundary_guard_layers:
        Number of grid layers near each boundary that are not jittered.
        ``1`` means the boundary layer and the first interior layer are kept
        regular.

    Returns
    -------
    MeshData
        Mesh container with nodes and triangles.
    """
    if not (0.0 <= jitter_fraction < 0.5):
        raise ValueError("jitter_fraction must satisfy 0 <= jitter_fraction < 0.5.")

    if not isinstance(boundary_guard_layers, int):
        raise ValueError("boundary_guard_layers must be an integer.")

    if boundary_guard_layers < 0:
        raise ValueError("boundary_guard_layers must be nonnegative.")

    geom = geometry_from_config(config)
    length_m = geom["length_m"]
    width_m = geom["width_m"]
    spacing_m = geom["target_spacing_m"]
    seed = int(geom["seed"])

    nx = max(3, int(math.ceil(length_m / spacing_m)) + 1)
    ny = max(3, int(math.ceil(width_m / spacing_m)) + 1)

    xs = np.linspace(0.0, length_m, nx)
    ys = np.linspace(-0.5 * width_m, 0.5 * width_m, ny)

    xx, yy = np.meshgrid(xs, ys, indexing="xy")
    nodes = np.column_stack([xx.ravel(), yy.ravel()])

    node_id = np.arange(nx * ny, dtype=np.int64).reshape(ny, nx)

    jj, ii = np.meshgrid(np.arange(ny), np.arange(nx), indexing="ij")
    protected = (
        (ii <= boundary_guard_layers)
        | (ii >= nx - 1 - boundary_guard_layers)
        | (jj <= boundary_guard_layers)
        | (jj >= ny - 1 - boundary_guard_layers)
    ).ravel()

    rng = np.random.default_rng(seed)

    dx = xs[1] - xs[0] if len(xs) > 1 else spacing_m
    dy = ys[1] - ys[0] if len(ys) > 1 else spacing_m
    nominal_step = min(dx, dy)
    jitter = jitter_fraction * nominal_step

    if jitter > 0.0:
        movable = ~protected
        nodes[movable, :] += rng.uniform(
            low=-jitter,
            high=jitter,
            size=(int(np.count_nonzero(movable)), 2),
        )

        nodes[:, 0] = np.clip(nodes[:, 0], 0.0, length_m)
        nodes[:, 1] = np.clip(nodes[:, 1], -0.5 * width_m, 0.5 * width_m)

        left = node_id[:, 0]
        right = node_id[:, -1]
        bottom = node_id[0, :]
        top = node_id[-1, :]

        nodes[left, 0] = 0.0
        nodes[right, 0] = length_m
        nodes[bottom, 1] = -0.5 * width_m
        nodes[top, 1] = 0.5 * width_m

    triangles = _triangulate_structured_grid_locally(nodes, node_id)
    triangles = orient_triangles_counterclockwise(nodes, triangles)

    return MeshData(
        nodes=np.asarray(nodes, dtype=float),
        triangles=np.asarray(triangles, dtype=np.int64),
        length_m=float(length_m),
        width_m=float(width_m),
        target_spacing_m=float(spacing_m),
        seed=int(seed),
        triangulation_method="protected_structured_local_delaunay",
        boundary_guard_layers=int(boundary_guard_layers),
    )


def _triangulate_structured_grid_locally(
    nodes: np.ndarray,
    node_id: np.ndarray,
) -> np.ndarray:
    """
    Triangulate each structured cell using the better local diagonal.

    For each quadrilateral cell, two diagonal choices are possible. We choose
    the one with the larger minimum triangle angle. This removes the global
    Qhull dependency while preserving a Delaunay-like local quality criterion.
    """
    ny, nx = node_id.shape
    triangles: list[tuple[int, int, int]] = []

    for j in range(ny - 1):
        for i in range(nx - 1):
            a = int(node_id[j, i])
            b = int(node_id[j, i + 1])
            c = int(node_id[j + 1, i])
            d = int(node_id[j + 1, i + 1])

            candidate_1 = [(a, b, d), (a, d, c)]
            candidate_2 = [(a, b, c), (b, d, c)]

            score_1 = _two_triangle_quality(nodes, candidate_1)
            score_2 = _two_triangle_quality(nodes, candidate_2)

            if score_2 > score_1:
                triangles.extend(candidate_2)
            elif score_1 > score_2:
                triangles.extend(candidate_1)
            else:
                if (i + j) % 2 == 0:
                    triangles.extend(candidate_1)
                else:
                    triangles.extend(candidate_2)

    return np.asarray(triangles, dtype=np.int64)


def _two_triangle_quality(
    nodes: np.ndarray,
    triangles: list[tuple[int, int, int]],
) -> float:
    """
    Quality score for two triangles.

    The score is the smallest internal angle among both triangles. Maximizing
    this quantity avoids unnecessarily skinny triangles.
    """
    return min(_triangle_min_angle_rad(nodes[list(tri)]) for tri in triangles)


def _triangle_min_angle_rad(points: np.ndarray) -> float:
    """
    Return the minimum internal angle of one triangle in radians.
    """
    p0, p1, p2 = points

    vectors = [
        (p1 - p0, p2 - p0),
        (p0 - p1, p2 - p1),
        (p0 - p2, p1 - p2),
    ]

    angles = []
    for u, v in vectors:
        nu = np.linalg.norm(u)
        nv = np.linalg.norm(v)
        if nu <= 0.0 or nv <= 0.0:
            return 0.0
        cosang = float(np.dot(u, v) / (nu * nv))
        cosang = max(-1.0, min(1.0, cosang))
        angles.append(math.acos(cosang))

    return float(min(angles))


def orient_triangles_counterclockwise(
    nodes: np.ndarray,
    triangles: np.ndarray,
) -> np.ndarray:
    """
    Return triangles with counter-clockwise orientation.
    """
    tri = np.array(triangles, dtype=np.int64, copy=True)

    p0 = nodes[tri[:, 0]]
    p1 = nodes[tri[:, 1]]
    p2 = nodes[tri[:, 2]]

    signed_area_2 = (
        (p1[:, 0] - p0[:, 0]) * (p2[:, 1] - p0[:, 1])
        -
        (p1[:, 1] - p0[:, 1]) * (p2[:, 0] - p0[:, 0])
    )

    flip = signed_area_2 < 0.0
    old_1 = tri[flip, 1].copy()
    tri[flip, 1] = tri[flip, 2]
    tri[flip, 2] = old_1

    return tri


def triangle_areas(nodes: np.ndarray, triangles: np.ndarray) -> np.ndarray:
    """
    Compute triangle areas in square meters.
    """
    p0 = nodes[triangles[:, 0]]
    p1 = nodes[triangles[:, 1]]
    p2 = nodes[triangles[:, 2]]

    return 0.5 * np.abs(
        (p1[:, 0] - p0[:, 0]) * (p2[:, 1] - p0[:, 1])
        -
        (p1[:, 1] - p0[:, 1]) * (p2[:, 0] - p0[:, 0])
    )


def mesh_summary(mesh: MeshData) -> dict[str, Any]:
    """
    Build a compact summary dictionary for diagnostics and manifests.
    """
    areas = triangle_areas(mesh.nodes, mesh.triangles)
    used_nodes = np.unique(mesh.triangles.reshape(-1))
    area_rectangle = float(mesh.length_m * mesh.width_m)
    area_total = float(np.sum(areas))

    relative_area_error = (
        abs(area_total - area_rectangle) / area_rectangle
        if area_rectangle > 0.0
        else float("nan")
    )

    return {
        "n_nodes": mesh.n_nodes,
        "n_used_nodes": int(used_nodes.size),
        "n_unused_nodes": int(mesh.n_nodes - used_nodes.size),
        "n_triangles": mesh.n_triangles,
        "length_m": mesh.length_m,
        "width_m": mesh.width_m,
        "target_spacing_m": mesh.target_spacing_m,
        "seed": mesh.seed,
        "triangulation_method": mesh.triangulation_method,
        "boundary_guard_layers": mesh.boundary_guard_layers,
        "area_total_from_triangles_m2": area_total,
        "area_rectangle_m2": area_rectangle,
        "area_relative_error": float(relative_area_error),
        "triangle_area_min_m2": float(np.min(areas)),
        "triangle_area_max_m2": float(np.max(areas)),
        "triangle_area_mean_m2": float(np.mean(areas)),
    }


def save_mesh_npz(mesh: MeshData, path: str | Path) -> Path:
    """
    Save mesh arrays to a compressed ``.npz`` file.
    """
    output = Path(path)
    output.parent.mkdir(parents=True, exist_ok=True)

    np.savez_compressed(
        output,
        nodes=mesh.nodes,
        triangles=mesh.triangles,
        length_m=np.array(mesh.length_m),
        width_m=np.array(mesh.width_m),
        target_spacing_m=np.array(mesh.target_spacing_m),
        seed=np.array(mesh.seed, dtype=np.int64),
        triangulation_method=np.array(mesh.triangulation_method),
        boundary_guard_layers=np.array(mesh.boundary_guard_layers, dtype=np.int64),
    )

    return output


def load_mesh_npz(path: str | Path) -> MeshData:
    """
    Load a mesh saved by :func:`save_mesh_npz`.
    """
    source = Path(path)

    with np.load(source) as data:
        if "triangulation_method" in data:
            method = str(data["triangulation_method"].item())
        else:
            method = "legacy"

        if "boundary_guard_layers" in data:
            guard_layers = int(data["boundary_guard_layers"])
        else:
            guard_layers = 0

        return MeshData(
            nodes=np.asarray(data["nodes"], dtype=float),
            triangles=np.asarray(data["triangles"], dtype=np.int64),
            length_m=float(data["length_m"]),
            width_m=float(data["width_m"]),
            target_spacing_m=float(data["target_spacing_m"]),
            seed=int(data["seed"]),
            triangulation_method=method,
            boundary_guard_layers=guard_layers,
        )