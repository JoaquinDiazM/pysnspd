"""pyTDGL-identical rectangular PRE mesh generation in SI units.

This module builds the rectangular nanowire through the same geometry and
mesh-generation path used by pyTDGL:

    tdgl.geometry.box
    -> tdgl.device.Polygon
    -> tdgl.device.Device.make_mesh
    -> tdgl.device.meshing.generate_mesh
    -> Mesh.from_triangulation

Only the unit convention differs: pySNSPD keeps coordinates in SI meters. The
compatibility Device therefore uses ``coherence_length = 1.0`` so the
pyTDGL-style dimensionless mesh is numerically equal to the SI mesh.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any, Mapping

import numpy as np

from pysnspd.gtdgl.finite_volume import Mesh
from pysnspd.gtdgl.geometry import box
from pysnspd.gtdgl.tdgl_compat import Device, Layer, Polygon
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
    smooth: int = 0
    min_points: int | None = None
    boundary_points: int = 101
    terminal_contact_mode: str = "normal_left_right"


def parameters_from_config(
    config: Mapping[str, Any],
    *,
    max_edge_length_m: float | None = None,
    min_angle_deg: float | None = None,
    smooth: int | None = None,
    min_points: int | None = None,
    boundary_points: int | None = None,
) -> PyTDGLLikeMeshParameters:
    """Resolve pyTDGL meshing parameters from a full or minimal config."""

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
        # Match pyTDGL's Device.make_mesh default. Smoothing should be explicit
        # for this backend via mesh.pytdgl_smooth.
        smooth = mesh_cfg.get("pytdgl_smooth", 0)

    if min_points is None and "pytdgl_min_points" in mesh_cfg:
        value = mesh_cfg.get("pytdgl_min_points")
        min_points = None if value is None else int(value)

    if boundary_points is None:
        boundary_points = int(mesh_cfg.get("pytdgl_boundary_points", 101))

    params = PyTDGLLikeMeshParameters(
        length_m=float(length_m),
        width_m=float(width_m),
        target_spacing_m=float(spacing_m),
        seed=int(seed),
        max_edge_length_m=float(max_edge_length_m),
        min_angle_deg=float(min_angle_deg),
        smooth=int(smooth),
        min_points=min_points,
        boundary_points=int(boundary_points),
    )
    _validate_parameters(params)
    return params


def generate_rectangular_pytdgl_like_mesh(
    config: Mapping[str, Any],
    *,
    max_edge_length_m: float | None = None,
    min_angle_deg: float | None = None,
    smooth: int | None = None,
    min_points: int | None = None,
    boundary_points: int | None = None,
) -> MeshData:
    """Generate a rectangular pyTDGL-style mesh and return pySNSPD MeshData."""

    params = parameters_from_config(
        config,
        max_edge_length_m=max_edge_length_m,
        min_angle_deg=min_angle_deg,
        smooth=smooth,
        min_points=min_points,
        boundary_points=boundary_points,
    )
    mesh = generate_rectangular_pytdgl_fvm_mesh_from_parameters(params)

    return MeshData(
        nodes=np.asarray(mesh.sites, dtype=float),
        triangles=orient_triangles_counterclockwise(mesh.sites, mesh.elements),
        length_m=float(params.length_m),
        width_m=float(params.width_m),
        target_spacing_m=float(params.target_spacing_m),
        seed=int(params.seed),
        triangulation_method="pytdgl_device_make_mesh_box_generate_mesh_exact_fv_v1",
        boundary_guard_layers=int(params.smooth),
    )


def generate_rectangular_pytdgl_fvm_mesh_from_parameters(
    params: PyTDGLLikeMeshParameters,
) -> Mesh:
    """Generate the full finite-volume mesh through pyTDGL Device.make_mesh."""

    _validate_parameters(params)

    film = Polygon(
        "film",
        points=rectangular_boundary_points(
            params.length_m,
            params.width_m,
            points=params.boundary_points,
        ),
    )

    # coherence_length=1 keeps Device._create_dimensionless_mesh numerically in
    # meters, while preserving the exact pyTDGL Device.make_mesh sequence.
    layer = Layer(
        london_lambda=1.0,
        coherence_length=1.0,
        thickness=1.0,
        conductivity=None,
        u=5.79,
        gamma=10.0,
        z0=0.0,
    )
    device = Device(
        "rectangular_nanowire",
        layer=layer,
        film=film,
        holes=None,
        terminals=None,
        length_units="m",
    )
    device.make_mesh(
        max_edge_length=params.max_edge_length_m,
        min_points=params.min_points,
        smooth=params.smooth,
        min_angle=params.min_angle_deg,
    )
    if device.mesh is None:
        raise RuntimeError("pyTDGL Device.make_mesh did not create a mesh.")
    return device.mesh


def rectangular_boundary_points(
    length_m: float,
    width_m: float,
    *,
    points: int = 101,
) -> np.ndarray:
    """Return the pyTDGL ``box`` boundary for the rectangular film."""

    length = float(length_m)
    width = float(width_m)
    if length <= 0.0:
        raise ValueError("length_m must be positive.")
    if width <= 0.0:
        raise ValueError("width_m must be positive.")
    if int(points) < 4:
        raise ValueError("points must be at least 4.")

    return box(
        width=length,
        height=width,
        points=int(points),
        center=(0.5 * length, 0.0),
        angle=0,
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
        "backend": "pytdgl_device_make_mesh_box_generate_mesh_exact_fv_v1",
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
    if params.boundary_points < 4:
        raise ValueError("boundary_points must be at least 4.")
