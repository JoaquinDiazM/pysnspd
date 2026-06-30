"""pyTDGL-faithful rectangular PRE mesh generation in SI units.

This module replaces the earlier jittered/staggered rectangular point cloud with
the same meshing path used by pyTDGL:

    polygon boundary
    -> meshpy.triangle.build
    -> Mesh.from_triangulation
    -> optional Mesh.smooth
    -> EdgeMesh + Voronoi control volumes.

The only deliberate difference from pyTDGL is units: every coordinate remains in
meters. No pyTDGL coherence-length scaling is performed here.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any, Mapping

import numpy as np

from pysnspd.gtdgl.finite_volume import Mesh
from pysnspd.gtdgl.finite_volume.meshing import generate_mesh
from pysnspd.mesh.delaunay import MeshData, orient_triangles_counterclockwise


@dataclass(frozen=True)
class PyTDGLLikeMeshParameters:
    """Parameters for pyTDGL-style rectangular meshing in SI meters."""

    length_m: float
    width_m: float
    target_spacing_m: float
    seed: int
    max_edge_length_m: float
    min_angle_deg: float = 32.5
    smooth: int = 100
    min_points: int | None = None
    terminal_contact_mode: str = "normal_left_right"


def parameters_from_config(
    config: Mapping[str, Any],
    *,
    jitter_fraction: float = 0.0,
    boundary_guard_layers: int = 1,
    max_edge_length_m: float | None = None,
    min_angle_deg: float | None = None,
    smooth: int | None = None,
    min_points: int | None = None,
) -> PyTDGLLikeMeshParameters:
    """Resolve pyTDGL-like meshing parameters from a full or minimal config.

    This follows pyTDGL's public ``Device.make_mesh`` contract as closely as the
    rectangular pySNSPD PRE stage allows:

    * ``max_edge_length`` controls Triangle refinement.
    * ``min_angle`` is passed to ``meshpy.triangle.build``.
    * ``smooth`` is an explicit Laplacian-smoothing iteration count.
    * the old pySNSPD ``jitter_fraction`` and ``boundary_guard_layers`` are
      accepted for CLI compatibility, but they do not silently change the
      pyTDGL meshing controls.

    Coordinates remain in SI meters.
    """

    del jitter_fraction, boundary_guard_layers

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

    if max_edge_length_m is None:
        max_edge_length_m = mesh_cfg.get(
            "pytdgl_max_edge_length_m",
            mesh_cfg.get("max_edge_length_m", spacing_m),
        )

    if min_angle_deg is None:
        min_angle_deg = mesh_cfg.get(
            "pytdgl_min_angle_deg",
            mesh_cfg.get("min_angle_deg", 32.5),
        )

    if smooth is None:
        smooth = mesh_cfg.get("pytdgl_smooth", mesh_cfg.get("smooth", 0))

    if min_points is None and "pytdgl_min_points" in mesh_cfg:
        value = mesh_cfg.get("pytdgl_min_points")
        min_points = None if value is None else int(value)

    params = PyTDGLLikeMeshParameters(
        length_m=length_m,
        width_m=width_m,
        target_spacing_m=spacing_m,
        seed=seed,
        max_edge_length_m=float(max_edge_length_m),
        min_angle_deg=float(min_angle_deg),
        smooth=int(smooth),
        min_points=min_points,
    )
    _validate_parameters(params)
    return params


def generate_rectangular_pytdgl_like_mesh(
    config: Mapping[str, Any],
    *,
    jitter_fraction: float = 0.0,
    boundary_guard_layers: int = 1,
    max_edge_length_m: float | None = None,
    min_angle_deg: float | None = None,
    smooth: int | None = None,
    min_points: int | None = None,
) -> MeshData:
    """Generate a rectangular pyTDGL-style mesh and return pySNSPD MeshData."""

    params = parameters_from_config(
        config,
        jitter_fraction=jitter_fraction,
        boundary_guard_layers=boundary_guard_layers,
        max_edge_length_m=max_edge_length_m,
        min_angle_deg=min_angle_deg,
        smooth=smooth,
        min_points=min_points,
    )
    mesh = generate_rectangular_pytdgl_fvm_mesh_from_parameters(params)

    return MeshData(
        nodes=np.asarray(mesh.sites, dtype=float),
        triangles=orient_triangles_counterclockwise(mesh.sites, mesh.elements),
        length_m=float(params.length_m),
        width_m=float(params.width_m),
        target_spacing_m=float(params.target_spacing_m),
        seed=int(params.seed),
        triangulation_method="pytdgl_generate_mesh_meshpy_triangle_v1",
        boundary_guard_layers=int(params.smooth),
    )


def generate_rectangular_pytdgl_fvm_mesh_from_parameters(
    params: PyTDGLLikeMeshParameters,
) -> Mesh:
    """Generate the full pyTDGL-like finite-volume Mesh for a rectangle."""

    _validate_parameters(params)

    poly_coords = rectangular_boundary_points(
        params.length_m,
        params.width_m,
        params.target_spacing_m,
    )

    points, triangles = generate_mesh(
        poly_coords=poly_coords,
        hole_coords=None,
        min_points=params.min_points,
        max_edge_length=params.max_edge_length_m,
        convex_hull=False,
        boundary=poly_coords,
        min_angle=params.min_angle_deg,
    )

    primary_mesh = Mesh.from_triangulation(
        points,
        triangles,
        create_submesh=False,
    )

    if params.smooth > 0:
        primary_mesh = primary_mesh.smooth(
            params.smooth,
            create_submesh=False,
        )

    return Mesh.from_triangulation(
        primary_mesh.sites,
        primary_mesh.elements,
        create_submesh=True,
    )


def rectangular_boundary_points(
    length_m: float,
    width_m: float,
    spacing_m: float,
) -> np.ndarray:
    """Return the rectangular film polygon vertices in meters.

    The polygon supplies only its geometric vertices, while ``meshpy.triangle``
    is responsible for inserting boundary and interior mesh sites during
    refinement. A dense, pre-resampled boundary can over-constrain Triangle and
    produce malformed Voronoi control cells near the first interior row.

    The curve is open: the first point is not repeated at the end.
    """

    del spacing_m

    length = float(length_m)
    half_w = 0.5 * float(width_m)

    return np.array(
        [
            [0.0, -half_w],
            [length, -half_w],
            [length, half_w],
            [0.0, half_w],
        ],
        dtype=float,
    )


def save_pytdgl_like_mesh_npz(mesh: Mesh, path: str | Path) -> Path:
    """Save full pyTDGL-like finite-volume mesh arrays to ``.npz``."""

    output = Path(path)
    output.parent.mkdir(parents=True, exist_ok=True)

    polygons = mesh.voronoi_polygons or []
    if polygons:
        split_indices = np.cumsum([len(p) for p in polygons[:-1]])
        polygons_flat = np.concatenate(polygons, axis=0)
    else:
        split_indices = np.array([], dtype=np.int64)
        polygons_flat = np.empty((0, 2), dtype=float)

    edge_mesh = mesh.edge_mesh
    if edge_mesh is None:
        raise ValueError("Mesh must have edge_mesh to save pyTDGL-like full mesh.")

    np.savez_compressed(
        output,
        sites=mesh.sites,
        elements=mesh.elements,
        boundary_indices=mesh.boundary_indices,
        areas=mesh.areas,
        dual_sites=mesh.dual_sites,
        voronoi_polygons_flat=polygons_flat,
        voronoi_split_indices=split_indices,
        edge_centers=edge_mesh.centers,
        edge_edges=edge_mesh.edges,
        edge_boundary_edge_indices=edge_mesh.boundary_edge_indices,
        edge_directions=edge_mesh.directions,
        edge_lengths=edge_mesh.edge_lengths,
        edge_dual_edge_lengths=edge_mesh.dual_edge_lengths,
    )

    return output


def build_pytdgl_like_mesh_summary(mesh: Mesh) -> dict[str, Any]:
    """Summary for the full pyTDGL-like finite-volume mesh."""

    edge_mesh = mesh.edge_mesh
    if edge_mesh is None:
        raise ValueError("Mesh must include an EdgeMesh.")

    return {
        "backend": "pytdgl_like_meshpy_triangle_fvm_v1",
        "n_sites": int(len(mesh.sites)),
        "n_elements": int(len(mesh.elements)),
        "n_boundary_sites": int(len(mesh.boundary_indices)),
        "n_edges": int(len(edge_mesh.edges)),
        "n_boundary_edges": int(len(edge_mesh.boundary_edge_indices)),
        "area_sum_m2": float(np.sum(mesh.areas)),
        "area_min_m2": float(np.min(mesh.areas)),
        "area_max_m2": float(np.max(mesh.areas)),
        "edge_length_min_m": float(np.min(edge_mesh.edge_lengths)),
        "edge_length_max_m": float(np.max(edge_mesh.edge_lengths)),
        "dual_edge_length_min_m": float(np.min(edge_mesh.dual_edge_lengths)),
        "dual_edge_length_max_m": float(np.max(edge_mesh.dual_edge_lengths)),
    }


def _validate_parameters(params: PyTDGLLikeMeshParameters) -> None:
    if params.length_m <= 0.0:
        raise ValueError("length_m must be positive.")
    if params.width_m <= 0.0:
        raise ValueError("width_m must be positive.")
    if params.target_spacing_m <= 0.0:
        raise ValueError("target_spacing_m must be positive.")
    if params.max_edge_length_m <= 0.0:
        raise ValueError("max_edge_length_m must be positive.")
    if params.min_angle_deg <= 0.0:
        raise ValueError("min_angle_deg must be positive.")
    if params.smooth < 0:
        raise ValueError("smooth must be non-negative.")