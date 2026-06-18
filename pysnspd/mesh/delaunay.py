"""
Delaunay mesh generation for pySNSPD.

OE2 scope:
- Build a reproducible rectangular nanowire mesh.
- Keep the geometry simple and explicit.
- Save/load mesh arrays in a format that later modules can reuse.

No superconducting physics is implemented here.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any, Mapping
import math

import numpy as np

try:
    from scipy.spatial import Delaunay
except ImportError as exc:  # pragma: no cover
    raise ImportError(
        "pySNSPD mesh generation requires scipy. "
        "Install it with: python -m pip install scipy"
    ) from exc

from pysnspd.config import validate_config


@dataclass(frozen=True)
class MeshData:
    """
    Container for a 2D nanowire Delaunay mesh.

    Coordinates are in meters. The rectangle is represented as

        x in [0, length_m],
        y in [-width_m/2, width_m/2].

    Attributes
    ----------
    nodes:
        Array of shape ``(n_nodes, 2)`` with columns ``x`` and ``y``.
    triangles:
        Array of shape ``(n_triangles, 3)`` with node indices.
    length_m:
        Nanowire length in meters.
    width_m:
        Nanowire width in meters.
    target_spacing_m:
        Requested nominal point spacing in meters.
    seed:
        Random seed used for interior jitter.
    """
    nodes: np.ndarray
    triangles: np.ndarray
    length_m: float
    width_m: float
    target_spacing_m: float
    seed: int

    @property
    def n_nodes(self) -> int:
        """Return the number of mesh nodes."""
        return int(self.nodes.shape[0])

    @property
    def n_triangles(self) -> int:
        """Return the number of Delaunay triangles."""
        return int(self.triangles.shape[0])

    @property
    def extent_m(self) -> tuple[float, float, float, float]:
        """Return mesh extent as ``(xmin, xmax, ymin, ymax)``."""
        return (0.0, self.length_m, -0.5 * self.width_m, 0.5 * self.width_m)


def geometry_from_config(config: Mapping[str, Any]) -> dict[str, float]:
    """
    Resolve the rectangular nanowire geometry from the config.

    Current template configs do not yet contain an explicit nanowire length.
    Therefore, OE2 uses the following priority:

    1. ``mesh.length_m`` if present.
    2. ``geometry.length_m`` if present.
    3. fallback ``6 * material.width_m``.

    The fallback is only a development default for mesh testing.
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
        length_m = 6.0 * width_m

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
) -> MeshData:
    """
    Generate a reproducible Delaunay mesh for a rectangular nanowire.

    Boundary nodes are kept exactly on the rectangle boundary. Only interior
    nodes are jittered. This makes boundary tagging robust while avoiding a
    perfectly structured triangulation in the interior.

    Parameters
    ----------
    config:
        Valid pySNSPD configuration dictionary.
    jitter_fraction:
        Interior random displacement as a fraction of the nominal spacing.
        Must satisfy ``0 <= jitter_fraction < 0.5``.

    Returns
    -------
    MeshData
        Mesh container with nodes and triangles.
    """
    if not (0.0 <= jitter_fraction < 0.5):
        raise ValueError("jitter_fraction must satisfy 0 <= jitter_fraction < 0.5.")

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

    is_left = np.isclose(nodes[:, 0], 0.0)
    is_right = np.isclose(nodes[:, 0], length_m)
    is_bottom = np.isclose(nodes[:, 1], -0.5 * width_m)
    is_top = np.isclose(nodes[:, 1], 0.5 * width_m)
    is_boundary = is_left | is_right | is_bottom | is_top

    rng = np.random.default_rng(seed)
    nominal_step = min(
        xs[1] - xs[0] if len(xs) > 1 else spacing_m,
        ys[1] - ys[0] if len(ys) > 1 else spacing_m,
    )
    jitter = jitter_fraction * nominal_step

    if jitter > 0.0:
        interior = ~is_boundary
        nodes[interior, :] += rng.uniform(
            low=-jitter,
            high=jitter,
            size=(int(np.count_nonzero(interior)), 2),
        )

        nodes[:, 0] = np.clip(nodes[:, 0], 0.0, length_m)
        nodes[:, 1] = np.clip(nodes[:, 1], -0.5 * width_m, 0.5 * width_m)

        nodes[is_left, 0] = 0.0
        nodes[is_right, 0] = length_m
        nodes[is_bottom, 1] = -0.5 * width_m
        nodes[is_top, 1] = 0.5 * width_m

    tri = Delaunay(nodes)
    triangles = np.asarray(tri.simplices, dtype=np.int64)
    triangles = orient_triangles_counterclockwise(nodes, triangles)

    return MeshData(
        nodes=np.asarray(nodes, dtype=float),
        triangles=triangles,
        length_m=float(length_m),
        width_m=float(width_m),
        target_spacing_m=float(spacing_m),
        seed=int(seed),
    )


def orient_triangles_counterclockwise(
    nodes: np.ndarray,
    triangles: np.ndarray,
) -> np.ndarray:
    """
    Return triangles with counter-clockwise orientation.

    This is useful later for finite-volume and edge-based operators.
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
    tri[flip, 1], tri[flip, 2] = tri[flip, 2], tri[flip, 1].copy()

    return tri


def triangle_areas(nodes: np.ndarray, triangles: np.ndarray) -> np.ndarray:
    """
    Compute triangle areas in square meters.
    """
    p0 = nodes[triangles[:, 0]]
    p1 = nodes[triangles[:, 1]]
    p2 = nodes[triangles[:, 2]]

    area = 0.5 * np.abs(
        (p1[:, 0] - p0[:, 0]) * (p2[:, 1] - p0[:, 1])
        -
        (p1[:, 1] - p0[:, 1]) * (p2[:, 0] - p0[:, 0])
    )
    return area


def mesh_summary(mesh: MeshData) -> dict[str, Any]:
    """
    Build a compact summary dictionary for diagnostics and manifests.
    """
    areas = triangle_areas(mesh.nodes, mesh.triangles)

    return {
        "n_nodes": mesh.n_nodes,
        "n_triangles": mesh.n_triangles,
        "length_m": mesh.length_m,
        "width_m": mesh.width_m,
        "target_spacing_m": mesh.target_spacing_m,
        "seed": mesh.seed,
        "area_total_from_triangles_m2": float(np.sum(areas)),
        "area_rectangle_m2": float(mesh.length_m * mesh.width_m),
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
    )

    return output


def load_mesh_npz(path: str | Path) -> MeshData:
    """
    Load a mesh saved by :func:`save_mesh_npz`.
    """
    source = Path(path)

    with np.load(source) as data:
        return MeshData(
            nodes=np.asarray(data["nodes"], dtype=float),
            triangles=np.asarray(data["triangles"], dtype=np.int64),
            length_m=float(data["length_m"]),
            width_m=float(data["width_m"]),
            target_spacing_m=float(data["target_spacing_m"]),
            seed=int(data["seed"]),
        )