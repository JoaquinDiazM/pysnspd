"""Mesh plotting utilities inspired by pyTDGL finite-volume views.

The active pySNSPD PRE-run stores meshes as ``pysnspd.mesh.delaunay.MeshData``
objects with node coordinates and triangular elements.  pyTDGL's finite-volume
``Mesh.plot`` can also draw dual/Voronoi data when available.  This module keeps
that visual spirit while remaining compatible with the lightweight pySNSPD mesh
objects currently written by the PRE-run.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Sequence

import matplotlib

matplotlib.use("Agg", force=True)
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
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
    site_color: Color = "C0",
    edge_color: Color = "0.10",
    centroid_color: Color = "C3",
    dual_edge_color: Color = "0.55",
    linewidth: float = 0.16,
    inset_linewidth: float = 0.34,
    site_size: float = 1.25,
    inset_site_size: float = 7.0,
    site_alpha: float = 0.90,
    edge_alpha: float = 0.42,
    dual_edge_alpha: float = 0.50,
    coordinate_scale: float = 1.0e9,
    coordinate_unit: str = "nm",
    title: str | None = None,
    inset_location: str = "upper right",
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
        Toggle optional dual/Voronoi information if present on ``mesh``.
    show_center_inset:
        Add a circular zoom inset centered in the wire bulk.
    show_stats_box:
        Add a compact mesh-statistics box in the corner opposite the inset.
    site_color, edge_color, centroid_color, dual_edge_color:
        Matplotlib-compatible colors.
    linewidth, inset_linewidth:
        Main and inset edge linewidths.
    site_size, inset_site_size:
        Main and inset scatter marker areas in pt^2.
    coordinate_scale:
        Factor applied to stored coordinates.  pySNSPD stores SI coordinates,
        so the default converts meters to nanometers.
    coordinate_unit:
        Axis unit label after scaling.
    title:
        Optional title.
    inset_location:
        One of ``"upper right"``, ``"upper left"``, ``"lower right"`` or
        ``"lower left"``.
    inset_radius:
        Circular inset radius in plotted coordinates.  If omitted, a radius is
        inferred from the target spacing and strip width.
    """
    output = _prepare_output(output_path)
    sites, elements = _mesh_sites_elements(mesh)
    sites_plot = coordinate_scale * sites
    edge_segments = _edge_segments(sites_plot, elements)
    polygons = _mesh_voronoi_polygons(mesh, coordinate_scale=coordinate_scale)

    fig, ax = plt.subplots(figsize=(7.4, 3.95), constrained_layout=False)
    fig.subplots_adjust(left=0.085, right=0.985, bottom=0.150, top=0.875)

    draw_mesh_pytdgl_style(
        mesh,
        ax=ax,
        show_sites=show_sites,
        show_edges=show_edges,
        show_dual_edges=show_dual_edges,
        show_voronoi_centroids=show_voronoi_centroids,
        site_color=site_color,
        edge_color=edge_color,
        centroid_color=centroid_color,
        dual_edge_color=dual_edge_color,
        linewidth=linewidth,
        site_size=site_size,
        site_alpha=site_alpha,
        edge_alpha=edge_alpha,
        dual_edge_alpha=dual_edge_alpha,
        coordinate_scale=coordinate_scale,
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
            voronoi_polygons=polygons,
            center=center,
            radius=radius,
            location=inset_location,
            show_sites=show_sites,
            show_edges=show_edges,
            show_dual_edges=show_dual_edges,
            site_color=site_color,
            edge_color=edge_color,
            dual_edge_color=dual_edge_color,
            site_size=inset_site_size,
            site_alpha=site_alpha,
            edge_alpha=0.72,
            dual_edge_alpha=dual_edge_alpha,
            linewidth=inset_linewidth,
        )

    if show_stats_box:
        _draw_stats_box(
            ax,
            mesh=mesh,
            n_edges=edge_segments.shape[0],
            coordinate_scale=coordinate_scale,
            coordinate_unit=coordinate_unit,
            location=_opposite_location(inset_location),
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
    site_color: Color = "C0",
    edge_color: Color = "0.10",
    centroid_color: Color = "C3",
    dual_edge_color: Color = "0.55",
    linewidth: float = 0.16,
    site_size: float = 1.25,
    site_alpha: float = 0.90,
    edge_alpha: float = 0.42,
    dual_edge_alpha: float = 0.50,
    coordinate_scale: float = 1.0e9,
) -> plt.Axes:
    """Draw a pyTDGL-style mesh on an existing Matplotlib axes."""
    if ax is None:
        _, ax = plt.subplots()

    sites, elements = _mesh_sites_elements(mesh)
    sites_plot = coordinate_scale * sites
    edge_segments = _edge_segments(sites_plot, elements)
    polygons = _mesh_voronoi_polygons(mesh, coordinate_scale=coordinate_scale)

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

    if show_dual_edges and polygons:
        dual_segments = _polygon_segments(polygons)
        if dual_segments.size:
            collection = LineCollection(
                dual_segments,
                colors=dual_edge_color,
                linewidths=max(linewidth, 0.12),
                alpha=dual_edge_alpha,
                zorder=2,
            )
            ax.add_collection(collection)

    if show_sites:
        ax.scatter(
            sites_plot[:, 0],
            sites_plot[:, 1],
            s=site_size,
            marker="o",
            color=site_color,
            linewidths=0.0,
            alpha=site_alpha,
            zorder=3,
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


def _draw_center_zoom_inset(
    ax: plt.Axes,
    *,
    sites_plot: np.ndarray,
    edge_segments: np.ndarray,
    voronoi_polygons: list[np.ndarray],
    center: tuple[float, float],
    radius: float,
    location: str,
    show_sites: bool,
    show_edges: bool,
    show_dual_edges: bool,
    site_color: Color,
    edge_color: Color,
    dual_edge_color: Color,
    site_size: float,
    site_alpha: float,
    edge_alpha: float,
    dual_edge_alpha: float,
    linewidth: float,
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

    if show_dual_edges and voronoi_polygons:
        dual_segments = _polygon_segments(voronoi_polygons)
        if dual_segments.size:
            collection = LineCollection(
                dual_segments,
                colors=dual_edge_color,
                linewidths=linewidth,
                alpha=dual_edge_alpha,
                zorder=2,
            )
            collection.set_clip_path(circle_clip)
            inset.add_collection(collection)

    if show_sites:
        sc = inset.scatter(
            sites_plot[:, 0],
            sites_plot[:, 1],
            s=site_size,
            marker="o",
            color=site_color,
            linewidths=0.0,
            alpha=site_alpha,
            zorder=3,
        )
        sc.set_clip_path(circle_clip)

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
        edgecolor="0.10",
        linestyle=":" ,
        linewidth=0.80,
        alpha=0.75,
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
    guard_layers = getattr(mesh, "boundary_guard_layers", None)
    method = str(getattr(mesh, "triangulation_method", "mesh"))

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
    if guard_layers is not None:
        lines.append(f"guard: {int(guard_layers)}")
    if method:
        lines.append(_compact_method_label(method))

    x, y, ha, va = _text_anchor(location)
    ax.text(
        x,
        y,
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


def _format_mesh_axes(
    ax: plt.Axes,
    sites_plot: np.ndarray,
    *,
    coordinate_unit: str,
    title: str,
) -> None:
    ax.set_title(title, fontsize=13)
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
    if elements.size == 0:
        return np.empty((0, 2, 2), dtype=float)
    edges = np.vstack(
        [
            elements[:, [0, 1]],
            elements[:, [1, 2]],
            elements[:, [2, 0]],
        ]
    )
    edges = np.sort(edges, axis=1)
    edges = np.unique(edges, axis=0)
    return sites_plot[edges]


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
        return float(max(4.0 * spacing, 0.10 * min_span))
    return float(0.13 * min_span)


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
        "lower right": [0.772, 0.065, 0.205, 0.385],
        "lower left": [0.025, 0.065, 0.205, 0.385],
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
        "lower right": "lower left",
        "lower left": "lower right",
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


def _compact_method_label(method: str) -> str:
    label = method.replace("pytdgl_generate_mesh_", "")
    label = label.replace("protected_structured_local_delaunay", "structured+delaunay")
    label = label.replace("_", " ")
    if len(label) > 31:
        label = label[:28].rstrip() + "..."
    return label


def _prepare_output(output_path: str | Path) -> Path:
    output = Path(output_path)
    output.parent.mkdir(parents=True, exist_ok=True)
    return output


__all__ = ["draw_mesh_pytdgl_style", "plot_mesh_pytdgl_style"]
