"""Mesh plotting utilities inspired by pyTDGL's finite-volume mesh plot.

The active pySNSPD PRE-run stores meshes as :class:`pysnspd.mesh.delaunay.MeshData`,
which contains only primal mesh data: node coordinates and triangular elements. The
pyTDGL ``Mesh.plot`` method can additionally draw Voronoi/dual edges when those are
stored on the mesh object. This module keeps that same visual interface where possible
while remaining compatible with pySNSPD's lightweight ``MeshData`` objects.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Sequence

import matplotlib.pyplot as plt
import numpy as np


def plot_mesh_pytdgl_style(
    mesh: Any,
    output_path: str | Path,
    *,
    dpi: int = 480,
    show_sites: bool = True,
    show_edges: bool = True,
    show_dual_edges: bool = True,
    show_voronoi_centroids: bool = False,
    site_color: str | Sequence[float] | None = None,
    edge_color: str | Sequence[float] | None = "k",
    centroid_color: str | Sequence[float] | None = None,
    dual_edge_color: str | Sequence[float] | None = "k",
    linewidth: float = 0.55,
    linestyle: str = "-",
    marker: str = ".",
    coordinate_scale: float = 1.0e9,
    coordinate_unit: str = "nm",
    title: str | None = None,
) -> Path:
    """Save a pyTDGL-style mesh figure.

    Parameters
    ----------
    mesh:
        Mesh-like object. Supported active pySNSPD objects expose ``nodes`` and
        ``triangles``. Full finite-volume mesh objects exposing ``sites``,
        ``elements`` and optionally ``voronoi_polygons`` are also supported.
    output_path:
        Destination PNG/PDF path.
    dpi:
        Figure resolution.
    show_sites:
        Whether to show mesh vertices.
    show_edges:
        Whether to show triangular primal edges.
    show_dual_edges:
        Whether to show Voronoi/dual edges when available on ``mesh``.
    show_voronoi_centroids:
        Whether to show Voronoi-cell centroids when dual polygons are available.
    site_color, edge_color, centroid_color, dual_edge_color:
        Matplotlib colors. Defaults intentionally mirror pyTDGL's ``Mesh.plot``
        style as closely as possible while using Matplotlib defaults for sites.
    linewidth, linestyle, marker:
        Line and marker style controls.
    coordinate_scale:
        Factor applied to stored coordinates before plotting. pySNSPD stores SI
        coordinates, so the default converts meters to nanometers.
    coordinate_unit:
        Label used on the axes after applying ``coordinate_scale``.
    title:
        Optional title. If omitted, a compact default is used.

    Returns
    -------
    pathlib.Path
        Path to the saved figure.
    """

    output = _prepare_output(output_path)
    sites, elements = _mesh_sites_elements(mesh)
    sites_plot = coordinate_scale * sites
    polygons = _mesh_voronoi_polygons(mesh, coordinate_scale=coordinate_scale)

    fig, ax = plt.subplots(figsize=(7.2, 4.2))
    _draw_mesh_on_axes(
        ax,
        sites_plot=sites_plot,
        elements=elements,
        voronoi_polygons=polygons,
        show_sites=show_sites,
        show_edges=show_edges,
        show_dual_edges=show_dual_edges,
        show_voronoi_centroids=show_voronoi_centroids,
        site_color=site_color,
        edge_color=edge_color,
        centroid_color=centroid_color,
        dual_edge_color=dual_edge_color,
        linewidth=linewidth,
        linestyle=linestyle,
        marker=marker,
    )

    ax.set_title(title or "Mesh: pyTDGL-style finite-volume view")
    ax.set_xlabel(f"x [{coordinate_unit}]")
    ax.set_ylabel(f"y [{coordinate_unit}]")
    ax.grid(False)
    fig.tight_layout()
    fig.savefig(output, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    return output


def draw_mesh_pytdgl_style(
    mesh: Any,
    ax: plt.Axes | None = None,
    *,
    show_sites: bool = True,
    show_edges: bool = True,
    show_dual_edges: bool = True,
    show_voronoi_centroids: bool = False,
    site_color: str | Sequence[float] | None = None,
    edge_color: str | Sequence[float] | None = "k",
    centroid_color: str | Sequence[float] | None = None,
    dual_edge_color: str | Sequence[float] | None = "k",
    linewidth: float = 0.55,
    linestyle: str = "-",
    marker: str = ".",
    coordinate_scale: float = 1.0e9,
) -> plt.Axes:
    """Draw a pyTDGL-style mesh on an existing Matplotlib axes.

    This mirrors the public behavior of pyTDGL's ``Mesh.plot`` and returns the
    axes object, which is useful for composition in future plotting functions.
    """

    if ax is None:
        _, ax = plt.subplots()

    sites, elements = _mesh_sites_elements(mesh)
    sites_plot = coordinate_scale * sites
    polygons = _mesh_voronoi_polygons(mesh, coordinate_scale=coordinate_scale)

    return _draw_mesh_on_axes(
        ax,
        sites_plot=sites_plot,
        elements=elements,
        voronoi_polygons=polygons,
        show_sites=show_sites,
        show_edges=show_edges,
        show_dual_edges=show_dual_edges,
        show_voronoi_centroids=show_voronoi_centroids,
        site_color=site_color,
        edge_color=edge_color,
        centroid_color=centroid_color,
        dual_edge_color=dual_edge_color,
        linewidth=linewidth,
        linestyle=linestyle,
        marker=marker,
    )


def _draw_mesh_on_axes(
    ax: plt.Axes,
    *,
    sites_plot: np.ndarray,
    elements: np.ndarray,
    voronoi_polygons: list[np.ndarray],
    show_sites: bool,
    show_edges: bool,
    show_dual_edges: bool,
    show_voronoi_centroids: bool,
    site_color: str | Sequence[float] | None,
    edge_color: str | Sequence[float] | None,
    centroid_color: str | Sequence[float] | None,
    dual_edge_color: str | Sequence[float] | None,
    linewidth: float,
    linestyle: str,
    marker: str,
) -> plt.Axes:
    ax.set_aspect("equal")

    x, y = sites_plot.T
    if show_edges:
        ax.triplot(
            x,
            y,
            elements,
            color=edge_color,
            linestyle=linestyle,
            linewidth=linewidth,
        )

    if show_dual_edges and voronoi_polygons:
        for polygon in voronoi_polygons:
            closed = _close_curve(polygon)
            ax.plot(
                closed[:, 0],
                closed[:, 1],
                color=dual_edge_color,
                linestyle=linestyle,
                linewidth=linewidth,
            )

    if show_sites:
        ax.plot(x, y, marker=marker, linestyle="", color=site_color)

    if show_voronoi_centroids and voronoi_polygons:
        centroids = np.asarray([_convex_polygon_centroid(p) for p in voronoi_polygons])
        if centroids.size:
            ax.plot(
                centroids[:, 0],
                centroids[:, 1],
                marker=marker,
                linestyle="",
                color=centroid_color,
            )

    return ax


def _mesh_sites_elements(mesh: Any) -> tuple[np.ndarray, np.ndarray]:
    if hasattr(mesh, "sites") and hasattr(mesh, "elements"):
        sites = np.asarray(getattr(mesh, "sites"), dtype=float)
        elements = np.asarray(getattr(mesh, "elements"), dtype=np.int64)
    elif hasattr(mesh, "nodes") and hasattr(mesh, "triangles"):
        sites = np.asarray(getattr(mesh, "nodes"), dtype=float)
        elements = np.asarray(getattr(mesh, "triangles"), dtype=np.int64)
    else:
        raise TypeError(
            "mesh must expose either (sites, elements) or (nodes, triangles)."
        )

    sites = np.asarray(sites).squeeze()
    elements = np.asarray(elements).squeeze()
    if sites.ndim != 2 or sites.shape[1] != 2:
        raise ValueError(f"Mesh sites/nodes must have shape (n, 2), got {sites.shape!r}.")
    if elements.ndim != 2 or elements.shape[1] != 3:
        raise ValueError(
            f"Mesh elements/triangles must have shape (m, 3), got {elements.shape!r}."
        )
    return sites, elements.astype(np.int64, copy=False)


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


def _prepare_output(output_path: str | Path) -> Path:
    output = Path(output_path)
    output.parent.mkdir(parents=True, exist_ok=True)
    return output


__all__ = ["draw_mesh_pytdgl_style", "plot_mesh_pytdgl_style"]
