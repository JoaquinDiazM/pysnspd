"""pyTDGL-like rectangular triangulation compatibility layer for pySNSPD."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any, Mapping

import numpy as np


@dataclass(frozen=True)
class MeshData:
    """Container for a 2D nanowire mesh in SI meters."""

    nodes: np.ndarray
    triangles: np.ndarray
    length_m: float
    width_m: float
    target_spacing_m: float
    seed: int
    triangulation_method: str = "pytdgl_generate_mesh_meshpy_triangle_v1"
    boundary_guard_layers: int = 0

    @property
    def n_nodes(self) -> int:
        return int(self.nodes.shape[0])

    @property
    def n_triangles(self) -> int:
        return int(self.triangles.shape[0])

    @property
    def extent_m(self) -> tuple[float, float, float, float]:
        return (0.0, self.length_m, -0.5 * self.width_m, 0.5 * self.width_m)


def geometry_from_config(config: Mapping[str, Any]) -> dict[str, float]:
    """Resolve the rectangular nanowire geometry from a full or minimal config."""

    material = config.get("material", {})
    mesh_cfg = config.get("mesh", {})
    geometry = config.get("geometry", {}) if isinstance(config.get("geometry", {}), Mapping) else {}

    width_m = float(material.get("width_m", mesh_cfg.get("width_m", 0.0)))
    spacing_m = float(mesh_cfg.get("target_spacing_m", material.get("target_spacing_m", 0.0)))
    seed = int(mesh_cfg.get("seed", 12345))

    if "length_m" in mesh_cfg:
        length_m = float(mesh_cfg["length_m"])
    elif "length_m" in geometry:
        length_m = float(geometry["length_m"])
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
) -> MeshData:
    """Generate the official pyTDGL-style rectangular mesh.

    The historical function name is kept so older PRE-stage code can keep using
    ``generate_rectangular_delaunay_mesh``. The implementation is no longer the
    old jittered Delaunay point cloud. It delegates to the pyTDGL-style meshpy
    finite-volume path and does not expose a ``jitter_fraction`` argument.
    """

    from pysnspd.mesh.pytdgl_like import generate_rectangular_pytdgl_like_mesh

    return generate_rectangular_pytdgl_like_mesh(config)


def orient_triangles_counterclockwise(nodes: np.ndarray, triangles: np.ndarray) -> np.ndarray:
    tri = np.array(triangles, dtype=np.int64, copy=True)
    p0 = nodes[tri[:, 0]]
    p1 = nodes[tri[:, 1]]
    p2 = nodes[tri[:, 2]]
    signed_area_2 = (
        (p1[:, 0] - p0[:, 0]) * (p2[:, 1] - p0[:, 1])
        - (p1[:, 1] - p0[:, 1]) * (p2[:, 0] - p0[:, 0])
    )
    flip = signed_area_2 < 0.0
    old_1 = tri[flip, 1].copy()
    tri[flip, 1] = tri[flip, 2]
    tri[flip, 2] = old_1
    return tri


def triangle_areas(nodes: np.ndarray, triangles: np.ndarray) -> np.ndarray:
    p0 = nodes[triangles[:, 0]]
    p1 = nodes[triangles[:, 1]]
    p2 = nodes[triangles[:, 2]]
    return 0.5 * np.abs(
        (p1[:, 0] - p0[:, 0]) * (p2[:, 1] - p0[:, 1])
        - (p1[:, 1] - p0[:, 1]) * (p2[:, 0] - p0[:, 0])
    )


def mesh_summary(mesh: MeshData) -> dict[str, Any]:
    areas = triangle_areas(mesh.nodes, mesh.triangles)
    used_nodes = np.unique(mesh.triangles.reshape(-1))
    area_rectangle = float(mesh.length_m * mesh.width_m)
    area_total = float(np.sum(areas))
    relative_area_error = abs(area_total - area_rectangle) / area_rectangle if area_rectangle > 0 else float("nan")

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
    source = Path(path)
    with np.load(source) as data:
        method = str(data["triangulation_method"].item()) if "triangulation_method" in data else "legacy"
        guard_layers = int(data["boundary_guard_layers"]) if "boundary_guard_layers" in data else 0
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
