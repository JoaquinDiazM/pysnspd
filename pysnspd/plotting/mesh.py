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


Color = str | Sequence[float] | None


def plot_mesh_pytdgl_style(
    mesh: Any,
    output_path: str | Path,
    *,
    dpi: int = 480,
    show_sites: bool = True,
    show_edges: bool = True,
    show_dual_edges: bool = True,
    show_voronoi_centroids: bool = False,
    show_center_inset: bool = True,
    show_stats_box: bool = True,
    show_node_class_legend: bool = True,
    interior_color: Color = "C0",
    current_color: Color = "C3",
    insulating_color: Color = "C2",
    edge_color: Color = "0.18",
    centroid_color: Color = "C1",
    dual_edge_color: Color = "0.38",
    linewidth: float = 0.15,
    inset_linewidth: float = 0.34,
    dual_linewidth: float = 0.30,
    site_size: float = 0.95,
    boundary_site_size: float = 2.30,
    inset_site_size: float = 5.7,
    inset_boundary_site_size: float = 8.0,
    site_alpha: float = 0.90,
    edge_alpha: float = 0.40,
    dual_edge_alpha: float = 0.55,
    coordinate_scale: float = 1.0e9,
    coordinate_unit: str = "nm",
    title: str | None = None,
    inset_location: str = "lower right",
    inset_radius: float | None = None,
) -> Path:
    """Save a compact pyTDGL-style mesh figure.

    Parameters
    ----------
    mesh:
        Mesh-like object exposing either ``(nodes, triangles)`` or
        ``(sites, elements)``.
    output_path:
        Destination image path.
    dpi:
        Figure resolution.
    show_sites, show_edges:
        Toggle primal mesh vertices and triangular edges.
    show_dual_edges, show_voronoi_centroids:
        Toggle optional finite-volume / Voronoi information.
    show_center_inset:
        Add a circular zoom inset centered in the wire bulk.
    show_stats_box:
        Add a compact mesh-statistics box.
    show_node_class_legend:
        Add a compact legend for current / insulating / interior nodes.
    coordinate_scale:
        Factor applied to stored coordinates. pySNSPD stores SI coordinates, so
        the default converts meters to nanometers.
    coordinate_unit:
        Axis unit label after scaling.
    title:
        Optional figure title.
    inset_location:
        One of ``"upper right"``, ``"upper left"``, ``"lower right"`` or
        ``"lower left"``.
    inset_radius:
        Circular inset radius in plotted coordinates. If omitted, a radius is
        inferred from the mesh spacing and strip width.
    """
    output = _prepare_output(output_path)
    sites, elements = _mesh_sites_elements(mesh)
    sites_plot = coordinate_scale * sites
    edge_segments = _edge_segments(sites_plot, elements)
    dual_segments = _dual_segments_from_triangulation(sites_plot, elements)
    polygons = _mesh_voronoi_polygons(mesh, coordinate_scale=coordinate_scale)
    node_classes = _classify_nodes(mesh, sites_plot, coordinate_scale=coordinate_scale)

    fig, ax = plt.subplots(figsize=(7.45, 3.95), constrained_layout=False)
    fig.subplots_adjust(left=0.086, right=0.985, bottom=0.150, top=0.875)

    # Main panel: primal mesh + colored node classes. The dual is reserved for
    # the inset so the global view remains readable.
    draw_mesh_pytdgl_style(
        mesh,
        ax=ax,
        show_sites=show_sites,
        show_edges=show_edges,
        show_dual_edges=False,
        show_voronoi_centroids=show_voronoi_centroids,
        interior_color=interior_color,
        current_color=current_color,
        insulating_color=insulating_color,
        edge_color=edge_color,
        centroid_color=centroid_color,
        dual_edge_color=dual_edge_color,
        linewidth=linewidth,
        dual_linewidth=dual_linewidth,
        site_size=site_size,
        boundary_site_size=boundary_site_size,
        site_alpha=site_alpha,
        edge_alpha=edge_alpha,
        dual_edge_alpha=dual_edge_alpha,
        coordinate_scale=coordinate_scale,
        node_classes=node_classes,
    )

    _format_mesh_axes(
        ax,
        sites_plot,
        coordinate_unit=coordinate_unit,
        title=title or "Mesh: pyTDGL-style finite-volume view",
    )

    radius = _resolve_inset_radius(
        mesh,
        sites_plot,
        coordinate_scale=coordinate_scale,
        requested_radius=inset_radius,
    )
    center = _bulk_center(sites_plot)

    if show_center_inset:
        _draw_center_zoom_inset(
            ax,
            sites_plot=sites_plot,
            edge_segments=edge_segments,
            dual_segments=dual_segments,
            voronoi_polygons=polygons,
            center=center,
            radius=radius,
            location=inset_location,
            show_sites=show_sites,
            show_edges=show_edges,
            show_dual_edges=show_dual_edges,
            site_alpha=site_alpha,
            edge_alpha=0.72,
            dual_edge_alpha=dual_edge_alpha,
            linewidth=inset_linewidth,
            dual_linewidth=max(inset_linewidth, dual_linewidth),
            interior_color=interior_color,
            current_color=current_color,
            insulating_color=insulating_color,
            edge_color=edge_color,
            dual_edge_color=dual_edge_color,
            site_size=inset_site_size,
            boundary_site_size=inset_boundary_site_size,
            node_classes=node_classes,
        )

    opposite = _opposite_location(inset_location)
    if show_stats_box:
        _draw_stats_box(
            ax,
            mesh=mesh,
            n_edges=edge_segments.shape[0],
            coordinate_scale=coordinate_scale,
            coordinate_unit=coordinate_unit,
            location=opposite,
        )
    if show_node_class_legend:
        _draw_node_class_legend(
            ax,
            node_classes=node_classes,
            interior_color=interior_color,
            current_color=current_color,
            insulating_color=insulating_color,
            location=opposite,
        )

    fig.savefig(output, dpi=dpi, bbox_inches="tight", pad_inches=0.05)
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
    interior_color: Color = "C0",
    current_color: Color = "C3",
    insulating_color: Color = "C2",
    edge_color: Color = "0.18",
    centroid_color: Color = "C1",
    dual_edge_color: Color = "0.38",
    linewidth: float = 0.15,
    dual_linewidth: float = 0.30,
    site_size: float = 0.95,
    boundary_site_size: float = 2.30,
    site_alpha: float = 0.90,
    edge_alpha: float = 0.40,
    dual_edge_alpha: float = 0.55,
    coordinate_scale: float = 1.0e9,
    node_classes: dict[str, np.ndarray] | None = None,
) -> plt.Axes:
    """Draw a pyTDGL-style mesh on an existing Matplotlib axes."""
    if ax is None:
        _, ax = plt.subplots()

    sites, elements = _mesh_sites_elements(mesh)
    sites_plot = coordinate_scale * sites
    edge_segments = _edge_segments(sites_plot, elements)
    polygons = _mesh_voronoi_polygons(mesh, coordinate_scale=coordinate_scale)
    dual_segments = _dual_segments_from_triangulation(sites_plot, elements)

    if node_classes is None:
        node_classes = _classify_nodes(mesh, sites_plot, coordinate_scale=coordinate_scale)

    ax.set_aspect("equal", adjustable="box")

    if show_edges and edge_segments.size:
        collection = LineCollection(
            edge_segments,
            colors=edge_color,
            linewidths=linewidth,
            alpha=edge_alpha,
            zorder=1,
        )
        ax.add_collection(collection)

    if show_dual_edges:
        if polygons:
            dual_segments_to_draw = _polygon_segments(polygons)
        else:
            dual_segments_to_draw = dual_segments
        if dual_segments_to_draw.size:
            collection = LineCollection(
                dual_segments_to_draw,
                colors=dual_edge_color,
                linewidths=max(dual_linewidth, 0.14),
                alpha=dual_edge_alpha,
                zorder=2,
            )
            ax.add_collection(collection)

    if show_sites:
        _scatter_node_classes(
            ax,
            sites_plot,
            node_classes=node_classes,
            interior_color=interior_color,
            current_color=current_color,
            insulating_color=insulating_color,
            site_size=site_size,
            boundary_site_size=boundary_site_size,
            alpha=site_alpha,
        )

    if show_voronoi_centroids and polygons:
        centroids = np.asarray([_convex_polygon_centroid(p) for p in polygons])
        if centroids.size:
            ax.scatter(
                centroids[:, 0],
                centroids[:, 1],
                s=max(site_size, 2.0),
                marker="o",
                color=centroid_color,
                linewidths=0.0,
                alpha=site_alpha,
                zorder=4,
            )

    return ax


def _scatter_node_classes(
    ax: plt.Axes,
    sites_plot: np.ndarray,
    *,
    node_classes: dict[str, np.ndarray],
    interior_color: Color,
    current_color: Color,
    insulating_color: Color,
    site_size: float,
    boundary_site_size: float,
    alpha: float,
) -> None:
    classes = [
        ("interior", interior_color, site_size, 3),
        ("insulating", insulating_color, boundary_site_size, 4),
        ("current", current_color, boundary_site_size, 5),
    ]
    for name, color, size, zorder in classes:
        mask = node_classes[name]
        if not np.any(mask):
            continue
        ax.scatter(
            sites_plot[mask, 0],
            sites_plot[mask, 1],
            s=size,
            marker="o",
            color=color,
            linewidths=0.0,
            alpha=alpha,
            zorder=zorder,
        )



def _draw_center_zoom_inset(
    ax: plt.Axes,
    *,
    sites_plot: np.ndarray,
    edge_segments: np.ndarray,
    dual_segments: np.ndarray,
    voronoi_polygons: list[np.ndarray],
    center: tuple[float, float],
    radius: float,
    location: str,
    show_sites: bool,
    show_edges: bool,
    show_dual_edges: bool,
    interior_color: Color,
    current_color: Color,
    insulating_color: Color,
    edge_color: Color,
    dual_edge_color: Color,
    node_classes: dict[str, np.ndarray],
    site_size: float,
    boundary_site_size: float,
    site_alpha: float,
    edge_alpha: float,
    dual_edge_alpha: float,
    linewidth: float,
    dual_linewidth: float,
) -> None:
    bounds = _inset_bounds(location)
    inset = ax.inset_axes(bounds)
    inset.set_aspect("equal", adjustable="box")
    inset.set_facecolor("white")

    cx, cy = center
    circle_clip = Circle((cx, cy), radius, transform=inset.transData)

    if show_edges and edge_segments.size:
        collection = LineCollection(
            edge_segments,
            colors=edge_color,
            linewidths=linewidth,
            alpha=edge_alpha,
            zorder=1,
        )
        collection.set_clip_path(circle_clip)
        inset.add_collection(collection)

    if show_dual_edges:
        if voronoi_polygons:
            dual_segments_to_draw = _polygon_segments(voronoi_polygons)
        else:
            dual_segments_to_draw = dual_segments
        if dual_segments_to_draw.size:
            collection = LineCollection(
                dual_segments_to_draw,
                colors=dual_edge_color,
                linewidths=max(dual_linewidth, 0.16),
                alpha=dual_edge_alpha,
                zorder=2,
            )
            collection.set_clip_path(circle_clip)
            inset.add_collection(collection)

    if show_sites:
        _scatter_node_classes(
            inset,
            sites_plot,
            node_classes=node_classes,
            interior_color=interior_color,
            current_color=current_color,
            insulating_color=insulating_color,
            site_size=site_size,
            boundary_site_size=boundary_site_size,
            alpha=site_alpha,
        )
        for coll in inset.collections[-3:]:
            try:
                coll.set_clip_path(circle_clip)
            except Exception:
                pass

    inset.set_xlim(cx - radius, cx + radius)
    inset.set_ylim(cy - radius, cy + radius)
    inset.set_xticks([])
    inset.set_yticks([])
    inset.grid(False)

    border = Circle(
        (cx, cy),
        radius,
        transform=inset.transData,
        fill=False,
        edgecolor="0.10",
        linewidth=0.80,
        zorder=10,
    )
    inset.add_patch(border)

    for spine in inset.spines.values():
        spine.set_visible(False)

    main_circle = Circle(
        (cx, cy),
        radius,
        transform=ax.transData,
        fill=False,
        edgecolor="0.18",
        linestyle=":",
        linewidth=0.78,
        alpha=0.72,
        zorder=8,
    )
    ax.add_patch(main_circle)
    inset.text(
        0.5,
        0.04,
        "center zoom",
        transform=inset.transAxes,
        ha="center",
        va="bottom",
        fontsize=6.5,
        color="0.20",
        bbox={"boxstyle": "round,pad=0.14", "facecolor": "white", "edgecolor": "none", "alpha": 0.82},
        zorder=11,
    )



def _draw_stats_box(
    ax: plt.Axes,
    *,
    mesh: Any,
    n_edges: int,
    coordinate_scale: float,
    coordinate_unit: str,
    location: str,
) -> None:
    sites, elements = _mesh_sites_elements(mesh)
    length = _maybe_scaled_attr(mesh, "length_m", coordinate_scale)
    width = _maybe_scaled_attr(mesh, "width_m", coordinate_scale)
    spacing = _maybe_scaled_attr(mesh, "target_spacing_m", coordinate_scale)

    lines = [
        "mesh stats",
        f"nodes: {sites.shape[0]}",
        f"triangles: {elements.shape[0]}",
        f"edges: {int(n_edges)}",
    ]
    if np.isfinite(length) and np.isfinite(width):
        lines.append(f"L×W: {_fmt_number(length)}×{_fmt_number(width)} {coordinate_unit}")
    if np.isfinite(spacing):
        lines.append(f"h: {_fmt_number(spacing)} {coordinate_unit}")

    x, y, ha, va = _text_anchor(location)
    y_stats = y if va == "top" else min(0.28, y + 0.23)
    ax.text(
        x,
        y_stats,
        "\n".join(lines),
        transform=ax.transAxes,
        ha=ha,
        va=va,
        fontsize=7.0,
        linespacing=1.18,
        color="0.08",
        bbox={
            "boxstyle": "round,pad=0.28",
            "facecolor": "white",
            "edgecolor": "0.25",
            "linewidth": 0.55,
            "alpha": 0.88,
        },
        zorder=20,
    )



def _draw_node_class_legend(
    ax: plt.Axes,
    *,
    node_classes: dict[str, np.ndarray],
    interior_color: Color,
    current_color: Color,
    insulating_color: Color,
    location: str,
) -> None:
    labels = [
        ("current", current_color, int(np.count_nonzero(node_classes["current"]))),
        ("insulating", insulating_color, int(np.count_nonzero(node_classes["insulating"]))),
        ("interior", interior_color, int(np.count_nonzero(node_classes["interior"]))),
    ]
    handles = [
        Line2D([0], [0], marker="o", linestyle="None", markersize=4.2, markerfacecolor=color, markeredgewidth=0.0)
        for _, color, _ in labels
    ]
    legend_labels = [f"{name} (N={count})" for name, _, count in labels]

    loc = location.lower().strip()
    if loc == "upper left":
        bbox = (0.015, 0.67)
        legend_loc = "upper left"
    elif loc == "upper right":
        bbox = (0.985, 0.67)
        legend_loc = "upper right"
    elif loc == "lower left":
        bbox = (0.015, 0.33)
        legend_loc = "lower left"
    else:
        bbox = (0.985, 0.33)
        legend_loc = "lower right"

    legend = ax.legend(
        handles,
        legend_labels,
        title="node classes",
        loc=legend_loc,
        bbox_to_anchor=bbox,
        fontsize=6.7,
        title_fontsize=6.9,
        frameon=True,
        facecolor="white",
        edgecolor="0.25",
        framealpha=0.88,
        borderpad=0.28,
        labelspacing=0.26,
        handlelength=0.9,
        handletextpad=0.45,
        borderaxespad=0.0,
    )
    legend.set_zorder(21)



def _format_mesh_axes(
    ax: plt.Axes,
    sites_plot: np.ndarray,
    *,
    coordinate_unit: str,
    title: str,
) -> None:
    #ax.set_title(title, fontsize=13) # For now commenting out the title to keep the figure compact
    ax.set_xlabel(f"x [{coordinate_unit}]")
    ax.set_ylabel(f"y [{coordinate_unit}]")
    ax.grid(False)
    ax.set_aspect("equal", adjustable="box")

    xmin, ymin = np.nanmin(sites_plot, axis=0)
    xmax, ymax = np.nanmax(sites_plot, axis=0)
    dx = max(float(xmax - xmin), 1.0)
    dy = max(float(ymax - ymin), 1.0)
    ax.set_xlim(float(xmin - 0.018 * dx), float(xmax + 0.018 * dx))
    ax.set_ylim(float(ymin - 0.050 * dy), float(ymax + 0.050 * dy))



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


__all__ = ["draw_mesh_pytdgl_style", "plot_mesh_pytdgl_style"]
