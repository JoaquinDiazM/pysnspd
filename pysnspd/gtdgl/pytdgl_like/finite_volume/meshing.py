"""pyTDGL mesh generator wrapper in SI units.

The public function ``generate_mesh`` intentionally keeps pyTDGL's signature.
It uses ``meshpy.triangle.build`` just as pyTDGL does.  pySNSPD supplies
coordinates in meters and receives coordinates in meters; no nondimensional
length conversion is performed here.
"""
from __future__ import annotations

import logging
from typing import List, Tuple, Union

import numpy as np
from scipy import spatial
from shapely.geometry.polygon import Polygon

from .util import ensure_unique, get_max_edge_length

logger = logging.getLogger("pysnspd.gtdgl.pytdgl_like.device")


def generate_mesh(
    poly_coords: np.ndarray,
    hole_coords: Union[List[np.ndarray], None] = None,
    min_points: Union[int, None] = None,
    max_edge_length: Union[float, None] = None,
    convex_hull: bool = False,
    boundary: Union[np.ndarray, None] = None,
    min_angle: float = 32.5,
    **kwargs,
) -> Tuple[np.ndarray, np.ndarray]:
    """Generates a Delaunay mesh for a given polygon.

    This follows ``tdgl.device.meshing.generate_mesh``.  Additional keyword
    arguments are passed to ``meshpy.triangle.build``.
    """
    try:
        from meshpy import triangle
    except ImportError as exc:
        raise ImportError(
            "pyTDGL-like meshing requires meshpy. Install it in the snspd "
            "environment, e.g. `python -m pip install meshpy` or "
            "`conda install -c conda-forge meshpy`."
        ) from exc

    poly_coords = ensure_unique(poly_coords)
    if hole_coords is None:
        hole_coords = []
    hole_coords = [ensure_unique(coords) for coords in hole_coords]

    coords = np.concatenate([poly_coords] + hole_coords, axis=0)
    xmin = coords[:, 0].min()
    dx = np.ptp(coords[:, 0])
    ymin = coords[:, 1].min()
    dy = np.ptp(coords[:, 1])
    r0 = np.array([[xmin, ymin]]) + np.array([[dx, dy]]) / 2

    coords = coords - r0
    indices = np.arange(len(poly_coords), dtype=int)
    if convex_hull:
        if boundary is not None:
            raise ValueError("Cannot have both boundary is not None and convex_hull = True.")
        facets = spatial.ConvexHull(coords).simplices
    else:
        if boundary is not None:
            boundary = list(map(tuple, ensure_unique(boundary - r0)))
            indices = np.array([i for i in indices if tuple(coords[i]) in boundary], dtype=int)
        facets = np.array([indices, np.roll(indices, -1)]).T

    for hole in hole_coords:
        hole_indices = np.arange(indices[-1] + 1, indices[-1] + 1 + len(hole), dtype=int)
        hole_facets = np.array([hole_indices, np.roll(hole_indices, -1)]).T
        indices = np.concatenate([indices, hole_indices], axis=0)
        facets = np.concatenate([facets, hole_facets], axis=0)

    mesh_info = triangle.MeshInfo()
    mesh_info.set_points(coords)
    mesh_info.set_facets(facets)
    if hole_coords:
        holes = [np.array(Polygon(hole).centroid.coords[0]) - r0.squeeze() for hole in hole_coords]
        mesh_info.set_holes(holes)

    kwargs = kwargs.copy()
    kwargs["min_angle"] = min_angle
    mesh = triangle.build(mesh_info=mesh_info, **kwargs)
    points = np.array(mesh.points) + r0
    triangles = np.array(mesh.elements, dtype=np.int64)
    if min_points is None and (max_edge_length is None or max_edge_length <= 0):
        return points, triangles

    kwargs["max_volume"] = dx * dy / 100
    i = 1
    if min_points is None:
        min_points = 0
    if max_edge_length is None or max_edge_length <= 0:
        max_edge_length = np.inf
    max_length = get_max_edge_length(points, triangles)
    while (len(points) < min_points) or (max_length > max_edge_length):
        mesh = triangle.build(mesh_info=mesh_info, **kwargs)
        points = np.array(mesh.points) + r0
        triangles = np.array(mesh.elements, dtype=np.int64)
        max_length = get_max_edge_length(points, triangles)
        logger.info(
            f"Iteration {i}: {len(points)} points, {len(triangles)} triangles, "
            f"max_edge_length: {max_length:.2e} (target: {max_edge_length:.2e})."
        )
        if np.isfinite(max_edge_length):
            kwargs["max_volume"] *= min(0.98, np.sqrt(max_edge_length / max_length))
        else:
            kwargs["max_volume"] *= 0.98
        i += 1
    return points, triangles
