"""Mesh plotting utilities inspired by pyTDGL finite-volume views.

The active pySNSPD PRE-run stores meshes as ``pysnspd.mesh.delaunay.MeshData``
objects with node coordinates and triangular elements. This module keeps the
visual spirit of pyTDGL's finite-volume mesh plot while remaining compatible
with the lightweight pySNSPD mesh objects currently written by the PRE-run.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Sequence

import matplotlib

matplotlib.use("Agg", force=True)
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.lines import Line2D
from matplotlib.patches import Circle
import numpy as np

from pysnspd.plotting.style import THESIS_WIDTH_IN, apply_thesis_style


Color = str | Sequence[float] | None


def _mesh_sites_elements(mesh: Any) -> tuple[np.ndarray, np.ndarray]:
    if hasattr(mesh, "sites") and hasattr(mesh, "elements"):
        sites = np.asarray(getattr(mesh, "sites"), dtype=float)
        elements = np.asarray(getattr(mesh, "elements"), dtype=np.int64)
    elif hasattr(mesh, "nodes") and hasattr(mesh, "triangles"):
        sites = np.asarray(getattr(mesh, "nodes"), dtype=float)
        elements = np.asarray(getattr(mesh, "triangles"), dtype=np.int64)
    else:
        raise TypeError("mesh must expose either (sites, elements) or (nodes, triangles).")

    sites = np.asarray(sites).squeeze()
    elements = np.asarray(elements).squeeze()
    if sites.ndim != 2 or sites.shape[1] != 2:
        raise ValueError(f"Mesh sites/nodes must have shape (n, 2), got {sites.shape!r}.")
    if elements.ndim != 2 or elements.shape[1] != 3:
        raise ValueError(
            f"Mesh elements/triangles must have shape (m, 3), got {elements.shape!r}."
        )
    return sites, elements.astype(np.int64, copy=False)



def _edge_segments(sites_plot: np.ndarray, elements: np.ndarray) -> np.ndarray:
    edges, _ = _unique_edges_with_neighbors(elements)
    if edges.size == 0:
        return np.empty((0, 2, 2), dtype=float)
    return sites_plot[edges]



def _dual_segments_from_triangulation(sites_plot: np.ndarray, elements: np.ndarray) -> np.ndarray:
    if elements.size == 0:
        return np.empty((0, 2, 2), dtype=float)

    centers = _triangle_circumcenters(sites_plot, elements)
    edges, neighbors = _unique_edges_with_neighbors(elements)
    p0 = sites_plot[edges[:, 0]]
    p1 = sites_plot[edges[:, 1]]
    midpoints = 0.5 * (p0 + p1)

    segments: list[np.ndarray] = []
    for idx in range(edges.shape[0]):
        tri_a = int(neighbors[idx, 0])
        tri_b = int(neighbors[idx, 1])
        if tri_a >= 0 and tri_b >= 0:
            segments.append(np.stack([centers[tri_a], centers[tri_b]], axis=0))
        elif tri_a >= 0:
            segments.append(np.stack([centers[tri_a], midpoints[idx]], axis=0))
        elif tri_b >= 0:
            segments.append(np.stack([centers[tri_b], midpoints[idx]], axis=0))
    if not segments:
        return np.empty((0, 2, 2), dtype=float)
    return np.stack(segments, axis=0)



def _unique_edges_with_neighbors(elements: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    if elements.size == 0:
        return np.empty((0, 2), dtype=np.int64), np.empty((0, 2), dtype=np.int64)

    edge_map: dict[tuple[int, int], list[int]] = {}
    for t_idx, tri in enumerate(np.asarray(elements, dtype=np.int64)):
        local_edges = [(int(tri[0]), int(tri[1])), (int(tri[1]), int(tri[2])), (int(tri[2]), int(tri[0]))]
        for a, b in local_edges:
            key = (a, b) if a < b else (b, a)
            edge_map.setdefault(key, []).append(int(t_idx))

    edge_keys = sorted(edge_map.keys())
    edges = np.asarray(edge_keys, dtype=np.int64)
    neighbors = -np.ones((edges.shape[0], 2), dtype=np.int64)
    for i, key in enumerate(edge_keys):
        local = edge_map[key]
        neighbors[i, : min(2, len(local))] = local[:2]
    return edges, neighbors



def _triangle_circumcenters(sites_plot: np.ndarray, elements: np.ndarray) -> np.ndarray:
    p0 = sites_plot[elements[:, 0]]
    p1 = sites_plot[elements[:, 1]]
    p2 = sites_plot[elements[:, 2]]

    x0, y0 = p0[:, 0], p0[:, 1]
    x1, y1 = p1[:, 0], p1[:, 1]
    x2, y2 = p2[:, 0], p2[:, 1]

    d = 2.0 * (x0 * (y1 - y2) + x1 * (y2 - y0) + x2 * (y0 - y1))
    centers = np.empty((elements.shape[0], 2), dtype=float)

    safe = np.abs(d) > 1.0e-18
    centers[:] = (p0 + p1 + p2) / 3.0
    if np.any(safe):
        u_x = ((x0**2 + y0**2) * (y1 - y2) + (x1**2 + y1**2) * (y2 - y0) + (x2**2 + y2**2) * (y0 - y1)) / d
        u_y = ((x0**2 + y0**2) * (x2 - x1) + (x1**2 + y1**2) * (x0 - x2) + (x2**2 + y2**2) * (x1 - x0)) / d
        centers[safe, 0] = u_x[safe]
        centers[safe, 1] = u_y[safe]
    return centers



def _mesh_voronoi_polygons(mesh: Any, *, coordinate_scale: float) -> list[np.ndarray]:
    raw = getattr(mesh, "voronoi_polygons", None)
    if raw is None:
        return []
    polygons: list[np.ndarray] = []
    for polygon in raw:
        arr = np.asarray(polygon, dtype=float).squeeze()
        if arr.ndim == 2 and arr.shape[1] == 2 and arr.shape[0] >= 3:
            polygons.append(coordinate_scale * arr)
    return polygons



def _polygon_segments(polygons: list[np.ndarray]) -> np.ndarray:
    segments: list[np.ndarray] = []
    for polygon in polygons:
        closed = _close_curve(polygon)
        if closed.shape[0] >= 2:
            segments.append(np.stack([closed[:-1], closed[1:]], axis=1))
    if not segments:
        return np.empty((0, 2, 2), dtype=float)
    return np.concatenate(segments, axis=0)



def _close_curve(points: np.ndarray) -> np.ndarray:
    pts = np.asarray(points, dtype=float)
    if pts.size == 0:
        return pts.reshape(0, 2)
    if np.allclose(pts[0], pts[-1]):
        return pts
    return np.vstack([pts, pts[0]])



def _convex_polygon_centroid(points: np.ndarray) -> np.ndarray:
    pts = np.asarray(points, dtype=float)
    if pts.shape[0] < 3:
        return np.nanmean(pts, axis=0)
    closed = _close_curve(pts)
    x = closed[:, 0]
    y = closed[:, 1]
    cross = x[:-1] * y[1:] - x[1:] * y[:-1]
    area2 = np.sum(cross)
    if abs(area2) <= 1.0e-300:
        return np.nanmean(pts, axis=0)
    cx = np.sum((x[:-1] + x[1:]) * cross) / (3.0 * area2)
    cy = np.sum((y[:-1] + y[1:]) * cross) / (3.0 * area2)
    return np.asarray([cx, cy], dtype=float)



def _classify_nodes(
    mesh: Any,
    sites_plot: np.ndarray,
    *,
    coordinate_scale: float,
) -> dict[str, np.ndarray]:
    x = sites_plot[:, 0]
    y = sites_plot[:, 1]
    xmin = float(np.nanmin(x))
    xmax = float(np.nanmax(x))
    ymin = float(np.nanmin(y))
    ymax = float(np.nanmax(y))

    spacing = _maybe_scaled_attr(mesh, "target_spacing_m", coordinate_scale)
    span = max(float(xmax - xmin), float(ymax - ymin), 1.0)
    if np.isfinite(spacing) and spacing > 0.0:
        tol = max(0.35 * spacing, 1.0e-6)
    else:
        tol = 1.0e-6 * span

    current = (np.abs(x - xmin) <= tol) | (np.abs(x - xmax) <= tol)
    insulating = (~current) & ((np.abs(y - ymin) <= tol) | (np.abs(y - ymax) <= tol))
    interior = ~(current | insulating)

    return {
        "current": current,
        "insulating": insulating,
        "interior": interior,
    }



def _bulk_center(sites_plot: np.ndarray) -> tuple[float, float]:
    return (float(np.nanmedian(sites_plot[:, 0])), float(np.nanmedian(sites_plot[:, 1])))



def _resolve_inset_radius(
    mesh: Any,
    sites_plot: np.ndarray,
    *,
    coordinate_scale: float,
    requested_radius: float | None,
) -> float:
    if requested_radius is not None:
        radius = float(requested_radius)
        if radius <= 0.0:
            raise ValueError("inset_radius must be positive when provided.")
        return radius

    spans = np.nanmax(sites_plot, axis=0) - np.nanmin(sites_plot, axis=0)
    min_span = float(max(np.nanmin(spans), 1.0))
    spacing = _maybe_scaled_attr(mesh, "target_spacing_m", coordinate_scale)
    if np.isfinite(spacing) and spacing > 0.0:
        return float(max(2.0 * spacing, 0.055 * min_span))
    return float(0.07 * min_span)



def _maybe_scaled_attr(mesh: Any, name: str, coordinate_scale: float) -> float:
    if not hasattr(mesh, name):
        return float("nan")
    try:
        return float(getattr(mesh, name)) * coordinate_scale
    except Exception:
        return float("nan")



def _inset_bounds(location: str) -> list[float]:
    loc = location.lower().strip()
    mapping = {
        "upper right": [0.772, 0.555, 0.205, 0.385],
        "upper left": [0.025, 0.555, 0.205, 0.385],
        "lower right": [0.772, 0.060, 0.205, 0.385],
        "lower left": [0.025, 0.060, 0.205, 0.385],
    }
    if loc not in mapping:
        raise ValueError(
            "inset_location must be one of: upper right, upper left, lower right, lower left."
        )
    return mapping[loc]



def _opposite_location(location: str) -> str:
    loc = location.lower().strip()
    mapping = {
        "upper right": "upper left",
        "upper left": "upper right",
        "lower right": "upper left",
        "lower left": "upper right",
    }
    return mapping.get(loc, "upper left")



def _text_anchor(location: str) -> tuple[float, float, str, str]:
    loc = location.lower().strip()
    mapping = {
        "upper left": (0.020, 0.965, "left", "top"),
        "upper right": (0.980, 0.965, "right", "top"),
        "lower left": (0.020, 0.035, "left", "bottom"),
        "lower right": (0.980, 0.035, "right", "bottom"),
    }
    return mapping.get(loc, mapping["upper left"])



def _fmt_number(value: float) -> str:
    val = float(value)
    if not np.isfinite(val):
        return "nan"
    if abs(val) >= 100:
        return f"{val:.0f}"
    if abs(val) >= 10:
        return f"{val:.1f}"
    if abs(val) >= 1:
        return f"{val:.2f}"
    return f"{val:.3g}"



def _prepare_output(output_path: str | Path) -> Path:
    output = Path(output_path)
    output.parent.mkdir(parents=True, exist_ok=True)
    return output
