"""pyTDGL-compatible geometry helpers used by pySNSPD meshing.

This module mirrors the geometry subset used by pyTDGL for constructing
polygonal device boundaries. Coordinates remain in SI meters in pySNSPD.

Source compatibility target:
    loganbvh/py-tdgl, ``tdgl/geometry.py``
    MIT License, Copyright (c) 2022-2026 Logan Bishop-Van Horn.
"""

from __future__ import annotations

from typing import Optional, Tuple

import numpy as np
from scipy.spatial import ConvexHull, QhullError


def rotation_matrix(angle_radians: float) -> np.ndarray:
    """Returns a 2D rotation matrix."""

    c = np.cos(angle_radians)
    s = np.sin(angle_radians)
    return np.array([[c, -s], [s, c]])


def rotate(coords: np.ndarray, angle_degrees: float) -> np.ndarray:
    """Rotates an array of ``(x, y)`` coordinates counterclockwise."""

    coords = np.asarray(coords)
    assert coords.ndim == 2
    assert coords.shape[1] == 2
    R = rotation_matrix(np.radians(angle_degrees))
    return (R @ coords.T).T


def box(
    width: float,
    height: Optional[float] = None,
    points: int = 101,
    center: Tuple[float, float] = (0, 0),
    angle: float = 0,
) -> np.ndarray:
    """Returns the coordinates for a rectangle with a given width and height.

    This follows pyTDGL's ``tdgl.geometry.box`` convention. Depending on the
    point allocation and endpoint handling, the first coordinate can also appear
    at the end of the returned boundary. Downstream pyTDGL meshing removes
    repeated coordinates with ``ensure_unique`` before constructing Triangle
    facets.
    """

    width = abs(width)
    if height is None:
        height = width
    height = abs(height)
    x0, y0 = center
    perimeter = 2 * (width + height)
    x_points = round(points * width / perimeter)
    y_points = round(points * height / perimeter)

    xs = np.concatenate(
        [
            width / 2 * np.ones(y_points),
            np.linspace(width / 2, -width / 2, x_points),
            -width / 2 * np.ones(y_points),
            np.linspace(-width / 2, width / 2, x_points),
        ]
    )
    ys = np.concatenate(
        [
            np.linspace(-height / 2, height / 2, y_points),
            height / 2 * np.ones(x_points),
            np.linspace(height / 2, -height / 2, y_points),
            -height / 2 * np.ones(x_points),
        ]
    )
    coords = np.array([xs, ys]).T + np.array([[x0, y0]])
    if angle:
        coords = rotate(coords, angle)
    return coords


def close_curve(points: np.ndarray) -> np.ndarray:
    """Close a curve if it is not already closed."""

    if not np.allclose(points[0], points[-1]):
        points = np.concatenate([points, points[:1]], axis=0)
    return points


def ensure_unique(coords: np.ndarray) -> np.ndarray:
    """Remove duplicate coordinates while preserving order."""

    coords = np.asarray(coords)
    _, ix = np.unique(coords, return_index=True, axis=0)
    coords = coords[np.sort(ix)]
    return coords


def get_convex_polygon_area(coords: np.ndarray) -> tuple[float, bool]:
    """Compute a convex polygon's area, or the area of its convex hull."""
    try:
        hull = ConvexHull(coords)
    except QhullError:
        return 0.0, True
    return float(hull.volume), len(hull.vertices) == len(coords)


def orient_convex_polygon(vertices: np.ndarray) -> np.ndarray:
    """Return counterclockwise-oriented vertices for a convex polygon."""
    diffs = vertices - vertices.mean(axis=0)
    return vertices[np.argsort(np.arctan2(diffs[:, 1], diffs[:, 0]))]
