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

from pysnspd.plotting.style import THESIS_DPI, THESIS_WIDTH_IN, apply_thesis_style

apply_thesis_style()


Color = str | Sequence[float] | None


def plot_mesh_pytdgl_style(
    mesh: Any,
    output_path: str | Path,
    *,
    dpi: int = THESIS_DPI,
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
    apply_thesis_style()
    output = _prepare_output(output_path)
    sites, elements = _mesh_sites_elements(mesh)
    sites_plot = coordinate_scale * sites
    edge_segments = _edge_segments(sites_plot, elements)
    dual_segments = _dual_segments_from_triangulation(sites_plot, elements)
    polygons = _mesh_voronoi_polygons(mesh, coordinate_scale=coordinate_scale)
    node_classes = _classify_nodes(mesh, sites_plot, coordinate_scale=coordinate_scale)

    fig, ax = plt.subplots(figsize=(THESIS_WIDTH_IN, 3.15), constrained_layout=False)
    fig.subplots_adjust(left=0.095, right=0.985, bottom=0.145, top=0.965)

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
            "Center zoom",
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
        "Mesh statistics",
        f"Nodes: {sites.shape[0]}",
        f"Triangles: {elements.shape[0]}",
        f"Edges: {int(n_edges)}",
    ]
    if np.isfinite(length) and np.isfinite(width):
        lines.append(f"L x W: {_fmt_number(length)} x {_fmt_number(width)} [{coordinate_unit}]")
    if np.isfinite(spacing):
        lines.append(f"Target h: {_fmt_number(spacing)} [{coordinate_unit}]")

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
        ("Current contacts", current_color, int(np.count_nonzero(node_classes["current"]))),
        ("Insulating boundary", insulating_color, int(np.count_nonzero(node_classes["insulating"]))),
        ("Interior", interior_color, int(np.count_nonzero(node_classes["interior"]))),
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
        title="Node classes",
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

from pysnspd.plotting.mesh_geometry import (
    Color,
    _bulk_center,
    _classify_nodes,
    _convex_polygon_centroid,
    _dual_segments_from_triangulation,
    _edge_segments,
    _fmt_number,
    _inset_bounds,
    _maybe_scaled_attr,
    _mesh_sites_elements,
    _mesh_voronoi_polygons,
    _opposite_location,
    _polygon_segments,
    _prepare_output,
    _resolve_inset_radius,
    _text_anchor,
)
