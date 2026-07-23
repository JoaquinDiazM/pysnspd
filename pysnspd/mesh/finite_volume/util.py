"""Finite-volume utilities following pyTDGL's implementation.

This file is a pySNSPD-local copy of the finite-volume utility subset used by
``tdgl.finite_volume.util`` in pyTDGL, with only import-path and packaging
adaptations. Coordinates are still whatever the caller provides; in pySNSPD PRE
runs they are SI meters.

Source compatibility target:
    loganbvh/py-tdgl, ``tdgl/finite_volume/util.py``
    MIT License, Copyright (c) 2022-2026 Logan Bishop-Van Horn.
"""

from __future__ import annotations

import logging
from collections import defaultdict
from typing import List, Tuple

import numpy as np
import scipy.sparse as sp
from scipy.spatial import Delaunay
from shapely.geometry import MultiLineString
from shapely.ops import orient, polygonize
from tqdm import tqdm

from pysnspd.mesh.geometry import (
    get_convex_polygon_area,
    orient_convex_polygon,
)

logger = logging.getLogger("tdgl.finite_volume")


def get_edges(elements: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """Finds the edges from a list of triangle indices.

    Args:
        elements: The triangle indices, shape ``(n, 3)``.

    Returns:
        A tuple containing an integer array of edges and a boolean array
        indicating whether each edge is on the boundary.
    """

    edges = np.concatenate([elements[:, e] for e in [(0, 1), (1, 2), (2, 0)]])
    edges = np.sort(edges, axis=1)
    edges, counts = np.unique(edges, return_counts=True, axis=0)
    return edges, counts == 1


def get_max_edge_length(points: np.ndarray, elements: np.ndarray) -> float:
    """Returns the maximum edge length in a triangulation."""

    edges = np.concatenate([elements[:, e] for e in [(0, 1), (1, 2), (2, 0)]])
    return np.linalg.norm(np.diff(points[edges], axis=1), axis=2).max()


def get_dual_edge_lengths(
    edge_centers: np.ndarray,
    elements: np.ndarray,
    dual_sites: np.ndarray,
    edges: np.ndarray,
    num_sites: float,
) -> np.ndarray:
    """Compute the lengths of the dual edges."""

    adj = make_adj_directed_tri_indices(elements, num_sites)
    edge_to_element = defaultdict(list)
    for i, j, v in zip(*sp.find(adj)):
        edge_to_element[frozenset((i, j))].append(v - 1)
    edge_to_element = dict(edge_to_element)

    dual_lengths = np.zeros(len(edge_centers), dtype=float)
    for i, edge in enumerate(edges):
        indices = edge_to_element[frozenset(edge)]
        if len(indices) == 1:
            # Boundary edge.
            dual_lengths[i] = np.linalg.norm(dual_sites[indices[0]] - edge_centers[i])
        else:
            # Inner edge.
            dual_lengths[i] = np.linalg.norm(
                dual_sites[indices[0]] - dual_sites[indices[1]]
            )
    return dual_lengths


def generate_voronoi_vertices(
    sites: np.ndarray, elements: np.ndarray
) -> Tuple[np.ndarray, np.ndarray]:
    """Compute the vertices of the Voronoi lattice.

    This is done by computing the circumcenters of the triangles in the
    tessellation.
    """

    A = sites[elements[:, 0]]
    B = sites[elements[:, 1]] - A
    C = sites[elements[:, 2]] - A

    D = 2 * B[:, 0] * C[:, 1] - 2 * B[:, 1] * C[:, 0]
    Ux = (C[:, 1] * (B**2).sum(axis=1) - B[:, 1] * (C**2).sum(axis=1)) / D
    Uy = (B[:, 0] * (C**2).sum(axis=1) - C[:, 0] * (B**2).sum(axis=1)) / D

    return np.array([Ux, Uy]).T + A


def make_adj_directed_tri_indices(elements: np.ndarray, num_sites: int) -> sp.csc_array:
    """Construct the directed adjacency matrix.

    Each element ``(i, j)`` represents an edge in the mesh, and the value at
    ``(i, j)`` is ``1 +`` the index of a triangle containing that edge.
    """

    t0 = elements[:, 0]
    t1 = elements[:, 1]
    t2 = elements[:, 2]
    i = np.column_stack([t0, t1, t2]).ravel()
    j = np.column_stack([t1, t2, t0]).ravel()
    data = np.repeat(np.arange(1, elements.shape[0] + 1), 3)
    return sp.csc_array((data, (i, j)), shape=(num_sites, num_sites))


def get_voronoi_polygon_indices(elements: np.ndarray, num_sites: int) -> List[np.ndarray]:
    """Find the polygons surrounding each site."""

    adj = make_adj_directed_tri_indices(elements, num_sites).tolil()
    return [np.array(tri) - 1 for tri in adj.data]


def compute_voronoi_polygon_areas(
    sites: np.ndarray,
    dual_sites: np.ndarray,
    boundary: np.ndarray,
    edges: np.ndarray,
    boundary_edge_indices: np.ndarray,
    polygons: List[np.ndarray],
) -> Tuple[np.ndarray, List[np.ndarray]]:
    """Compute the areas of the surrounding polygons.

    Areas of boundary points are handled by adding additional points on the
    boundary to make a convex polygon.
    """

    boundary_set = set(boundary)
    boundary_edges = edges[boundary_edge_indices]
    areas = np.zeros(len(polygons), dtype=float)
    voronoi_sites = []
    warning_str = (
        "Malformed Voronoi cell surrounding boundary site {site}."
        " Try changing the number of boundary mesh sites using"
        " Polygon.resample() or Polygon.buffer(eps) where eps"
        " is 0 or a small positive float."
    )

    for site, polygon in enumerate(
        tqdm(polygons, desc="Constructing Voronoi polygons")
    ):
        poly = dual_sites[polygon]

        if site not in boundary_set:
            areas[site], is_convex = get_convex_polygon_area(poly)
            if not is_convex:
                raise ValueError(warning_str.format(site=site))
            voronoi_sites.append(orient_convex_polygon(poly))
            continue

        connected_boundary_edges = boundary_edges[(boundary_edges == site).any(axis=1)]
        midpoints = sites[connected_boundary_edges].mean(axis=1)
        coords = orient_convex_polygon(np.concatenate([poly, midpoints], axis=0))
        coords = [tuple(xy) for xy in coords]
        indices = sorted([coords.index(tuple(mid)) for mid in midpoints])
        if indices[1] == indices[0] + 1:
            coords.insert(indices[1], sites[site])
        else:
            if indices[0] != 0:
                logger.warning(warning_str.format(site=site))
            coords.append(sites[site])
        poly = np.array(coords)
        areas[site], is_convex = get_convex_polygon_area(poly)
        if not is_convex:
            triangle_area, is_convex = get_convex_polygon_area(
                np.concatenate([midpoints, [sites[site]]], axis=0)
            )
            assert is_convex
            areas[site] -= triangle_area
        voronoi_sites.append(poly)
    return areas, voronoi_sites
