"""Finite-volume utilities copied/adapted from pyTDGL.

This module intentionally mirrors ``tdgl.finite_volume.util`` as closely as
possible while keeping all coordinates in pySNSPD's SI units.  The only
adaptations are packaging/import paths and small numerical guards/documentation.

Reference: loganbvh/py-tdgl, MIT license.
"""
from __future__ import annotations

import logging
from collections import defaultdict
from typing import List, Tuple

import numpy as np
import scipy.sparse as sp
from scipy.spatial import ConvexHull, Delaunay, QhullError
from shapely.geometry import MultiLineString
from shapely.ops import orient, polygonize

logger = logging.getLogger("pysnspd.gtdgl.pytdgl_like.finite_volume")


def close_curve(points: np.ndarray) -> np.ndarray:
    """Return a closed version of an ``(n, 2)`` coordinate array."""
    points = np.asarray(points, dtype=float)
    if len(points) == 0:
        return points
    if np.allclose(points[0], points[-1]):
        return points
    return np.vstack([points, points[0]])


def ensure_unique(points: np.ndarray, *, decimals: int = 14) -> np.ndarray:
    """Remove repeated coordinate rows while preserving first occurrence order."""
    points = np.asarray(points, dtype=float)
    if points.size == 0:
        return points.reshape(0, 2)
    rounded = np.round(points, decimals=decimals)
    _, indices = np.unique(rounded, axis=0, return_index=True)
    indices = np.sort(indices)
    return points[indices]


def get_edges(elements: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """Finds the edges from a list of triangle indices.

    Args:
        elements: The triangle indices, shape ``(n, 3)``.

    Returns:
        A tuple containing an integer array of edges and a boolean array
        indicating whether each edge is on the boundary.
    """
    elements = np.asarray(elements, dtype=np.int64)
    edges = np.concatenate([elements[:, e] for e in [(0, 1), (1, 2), (2, 0)]])
    edges = np.sort(edges, axis=1)
    edges, counts = np.unique(edges, return_counts=True, axis=0)
    return edges, counts == 1


def get_edge_lengths(points: np.ndarray, elements: np.ndarray) -> np.ndarray:
    """Returns the lengths of all unique edges in a triangulation."""
    edges, _ = get_edges(elements)
    return np.linalg.norm(np.diff(points[edges], axis=1), axis=2).squeeze()


def get_max_edge_length(points: np.ndarray, elements: np.ndarray) -> float:
    """Returns the maximum edge length in a triangulation."""
    elements = np.asarray(elements, dtype=np.int64)
    edges = np.concatenate([elements[:, e] for e in [(0, 1), (1, 2), (2, 0)]])
    return float(np.linalg.norm(np.diff(points[edges], axis=1), axis=2).max())


def get_dual_edge_lengths(
    edge_centers: np.ndarray,
    elements: np.ndarray,
    dual_sites: np.ndarray,
    edges: np.ndarray,
    num_sites: int,
) -> np.ndarray:
    """Compute the lengths of the Voronoi dual edges."""
    adj = make_adj_directed_tri_indices(elements, int(num_sites))
    edge_to_element: dict[frozenset[int], list[int]] = defaultdict(list)
    for i, j, v in zip(*sp.find(adj)):
        edge_to_element[frozenset((int(i), int(j)))].append(int(v) - 1)
    dual_lengths = np.zeros(len(edge_centers), dtype=float)
    for i, edge in enumerate(edges):
        indices = edge_to_element[frozenset(map(int, edge))]
        if len(indices) == 1:
            dual_lengths[i] = np.linalg.norm(dual_sites[indices[0]] - edge_centers[i])
        else:
            dual_lengths[i] = np.linalg.norm(dual_sites[indices[0]] - dual_sites[indices[1]])
    return dual_lengths


def generate_voronoi_vertices(sites: np.ndarray, elements: np.ndarray) -> np.ndarray:
    """Compute circumcenters of all triangles in the tessellation."""
    sites = np.asarray(sites, dtype=float)
    elements = np.asarray(elements, dtype=np.int64)
    A = sites[elements[:, 0]]
    B = sites[elements[:, 1]] - A
    C = sites[elements[:, 2]] - A
    D = 2 * B[:, 0] * C[:, 1] - 2 * B[:, 1] * C[:, 0]
    tiny = np.finfo(float).tiny
    D = np.where(np.abs(D) < tiny, np.sign(D) * tiny + (D == 0) * tiny, D)
    Ux = (C[:, 1] * (B**2).sum(axis=1) - B[:, 1] * (C**2).sum(axis=1)) / D
    Uy = (B[:, 0] * (C**2).sum(axis=1) - C[:, 0] * (B**2).sum(axis=1)) / D
    return np.array([Ux, Uy]).T + A


def make_adj_directed_tri_indices(elements: np.ndarray, num_sites: int) -> sp.csc_array:
    """Construct the directed adjacency matrix storing triangle indices + 1."""
    elements = np.asarray(elements, dtype=np.int64)
    t0 = elements[:, 0]
    t1 = elements[:, 1]
    t2 = elements[:, 2]
    i = np.column_stack([t0, t1, t2]).ravel()
    j = np.column_stack([t1, t2, t0]).ravel()
    data = np.repeat(np.arange(1, elements.shape[0] + 1), 3)
    return sp.csc_array((data, (i, j)), shape=(int(num_sites), int(num_sites)))


def get_voronoi_polygon_indices(elements: np.ndarray, num_sites: int) -> List[np.ndarray]:
    """Find Voronoi-vertex indices surrounding each mesh site."""
    adj = make_adj_directed_tri_indices(elements, int(num_sites)).tolil()
    return [np.array(tri, dtype=np.int64) - 1 for tri in adj.data]


def compute_voronoi_polygon_areas(
    sites: np.ndarray,
    dual_sites: np.ndarray,
    boundary: np.ndarray,
    edges: np.ndarray,
    boundary_edge_indices: np.ndarray,
    polygons: List[np.ndarray],
) -> Tuple[np.ndarray, List[np.ndarray]]:
    """Compute Voronoi cell areas and counter-clockwise polygon vertices."""
    boundary_set = set(np.asarray(boundary, dtype=np.int64))
    boundary_edges = np.asarray(edges, dtype=np.int64)[np.asarray(boundary_edge_indices, dtype=np.int64)]
    areas = np.zeros(len(polygons), dtype=float)
    voronoi_sites: list[np.ndarray] = []
    warning_str = (
        "Malformed Voronoi cell surrounding boundary site {site}. Try changing "
        "the number of boundary mesh sites using boundary resampling."
    )
    for site, polygon in enumerate(polygons):
        poly = np.asarray(dual_sites[np.asarray(polygon, dtype=np.int64)], dtype=float)
        if site not in boundary_set:
            areas[site], is_convex = get_convex_polygon_area(poly)
            if not is_convex:
                raise ValueError(warning_str.format(site=site))
            voronoi_sites.append(orient_convex_polygon(poly))
            continue

        connected_boundary_edges = boundary_edges[(boundary_edges == site).any(axis=1)]
        if len(connected_boundary_edges) < 2:
            # Fallback for pathological boundary points: use the site and all
            # available Voronoi vertices.  This mirrors pyTDGL's preference for
            # exposing malformed cells while keeping the PRE-run inspectable.
            coords = orient_convex_polygon(np.vstack([poly, sites[[site]]]))
            areas[site], _ = get_convex_polygon_area(coords)
            voronoi_sites.append(coords)
            continue
        midpoints = sites[connected_boundary_edges].mean(axis=1)
        coords = orient_convex_polygon(np.concatenate([poly, midpoints], axis=0))
        coords_list = [tuple(xy) for xy in coords]
        mid_list = [tuple(xy) for xy in midpoints]
        try:
            indices = sorted([coords_list.index(mid) for mid in mid_list])
            coords_mut = list(coords_list)
            if indices[1] == indices[0] + 1:
                coords_mut.insert(indices[1], tuple(sites[site]))
            else:
                if indices[0] != 0:
                    logger.warning(warning_str.format(site=site))
                coords_mut.append(tuple(sites[site]))
            poly_full = np.array(coords_mut, dtype=float)
        except ValueError:
            logger.warning(warning_str.format(site=site))
            poly_full = np.vstack([coords, sites[[site]]])
        areas[site], is_convex = get_convex_polygon_area(poly_full)
        if not is_convex and len(midpoints) >= 2:
            triangle_area, is_tri_convex = get_convex_polygon_area(
                np.concatenate([midpoints[:2], [sites[site]]], axis=0)
            )
            assert is_tri_convex
            areas[site] -= triangle_area
        voronoi_sites.append(poly_full)
    return areas, voronoi_sites


def get_convex_polygon_area(coords: np.ndarray) -> Tuple[float, bool]:
    """Compute area of a convex polygon or convex hull."""
    coords = np.asarray(coords, dtype=float)
    try:
        hull = ConvexHull(coords)
    except QhullError:
        return 0.0, True
    is_convex = len(hull.vertices) == len(coords)
    return float(hull.volume), bool(is_convex)


def triangle_areas(points: np.ndarray, triangles: np.ndarray) -> np.ndarray:
    """Calculate signed areas of triangular elements."""
    xy = np.asarray(points, dtype=float)[np.asarray(triangles, dtype=np.int64)]
    s = xy[:, [2, 0]] - xy[:, [1, 2]]
    a = np.linalg.det(s)
    return a * 0.5


def orient_convex_polygon(vertices: np.ndarray) -> np.ndarray:
    """Return vertices sorted counterclockwise."""
    vertices = np.asarray(vertices, dtype=float)
    if len(vertices) <= 2:
        return vertices
    diffs = vertices - vertices.mean(axis=0)
    return vertices[np.argsort(np.arctan2(diffs[:, 1], diffs[:, 0]))]


def convex_polygon_centroid(points: np.ndarray) -> Tuple[float, float]:
    """Calculate the centroid of a convex polygon."""
    points = np.asarray(points, dtype=float)
    triangles = Delaunay(points).simplices
    areas = triangle_areas(points, triangles)
    centroids = points[triangles].mean(axis=1)
    return tuple(np.average(centroids, weights=areas, axis=0))


def get_oriented_boundary(points: np.ndarray, boundary_edges: np.ndarray) -> List[np.ndarray]:
    """Return arrays of boundary vertex indices ordered counterclockwise."""
    points = np.asarray(points, dtype=float)
    points_list = [tuple(xy) for xy in points]
    edges = MultiLineString([points[np.asarray(edge, dtype=np.int64), :] for edge in boundary_edges])
    polygons = list(polygonize(edges))
    polygon_indices = []
    for p in polygons:
        polygon = orient(p)
        indices = np.array([points_list.index(tuple(xy)) for xy in polygon.exterior.coords])
        polygon_indices.append(indices[:-1])
    return polygon_indices
