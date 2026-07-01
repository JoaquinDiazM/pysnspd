"""PRE-run diagnostic plots for mesh, edges and Usadel calibration.

These plots are intentionally cheap and deterministic. They are meant to be
looked at before starting the SS/PHOTON stages so geometry, boundary tags and
Matsubara supercurrent tables can be checked independently of the dynamic
solver.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import numpy as np

from pysnspd.mesh.delaunay import MeshData, triangle_areas
from pysnspd.mesh.edges import EdgeData


def write_pre_diagnostic_plots(
    *,
    mesh: MeshData,
    edge_data: EdgeData,
    usadel_catalog: Any,
    output_dir: str | Path,
    dpi: int = 480,
) -> dict[str, str]:
    """Write the standard PRE diagnostic plots and return their paths.

    Parameters
    ----------
    mesh:
        Rectangular pyTDGL-style finite-volume triangulation stored in SI units.
    edge_data:
        Edge connectivity and boundary tags derived from ``mesh``.
    usadel_catalog:
        Dirty-limit Usadel catalogue returned by
        :func:`pysnspd.usadel.catalog.build_usadel_catalog_from_config`.
    output_dir:
        Directory where ``.png`` files will be written.
    dpi:
        Output resolution for all plots.
    """

    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)

    paths = {
        "mesh_boundary_tags_png": plot_mesh_boundary_tags(
            mesh,
            edge_data,
            out / "mesh_boundary_tags.png",
            dpi=dpi,
        ),
        "mesh_triangle_area_hist_png": plot_triangle_area_histogram(
            mesh,
            out / "mesh_triangle_area_hist.png",
            dpi=dpi,
        ),
        "mesh_edge_length_hist_png": plot_edge_length_histogram(
            edge_data,
            out / "mesh_edge_length_hist.png",
            dpi=dpi,
        ),
        "usadel_supercurrent_curve_png": plot_usadel_supercurrent_curve(
            usadel_catalog,
            out / "usadel_supercurrent_curve.png",
            dpi=dpi,
        ),
    }

    return {key: str(value) for key, value in paths.items()}


def plot_mesh_boundary_tags(
    mesh: MeshData,
    edge_data: EdgeData,
    output_path: str | Path,
    *,
    dpi: int = 480,
) -> Path:
    """Plot the triangulation and overlay boundary-edge tags."""

    output = _prepare_output(output_path)
    nodes_nm = np.asarray(mesh.nodes, dtype=float) * 1.0e9

    fig, ax = plt.subplots(figsize=(8.0, 4.2))
    ax.set_aspect("equal")
    ax.triplot(
        nodes_nm[:, 0],
        nodes_nm[:, 1],
        np.asarray(mesh.triangles, dtype=np.int64),
        linewidth=0.25,
        alpha=0.35,
    )

    boundary = np.asarray(edge_data.is_boundary, dtype=bool)
    tags = np.asarray(edge_data.tags).astype(str)
    edges = np.asarray(edge_data.edges, dtype=np.int64)

    ordered_tags = ["left", "right", "bottom", "top", "boundary_unknown"]
    for tag in ordered_tags:
        mask = boundary & (tags == tag)
        if not np.any(mask):
            continue
        _plot_edge_segments(
            ax,
            nodes_nm,
            edges[mask],
            label=f"{tag} ({int(np.count_nonzero(mask))})",
            linewidth=1.0,
        )

    unknown_boundary = boundary & ~np.isin(tags, ordered_tags)
    if np.any(unknown_boundary):
        _plot_edge_segments(
            ax,
            nodes_nm,
            edges[unknown_boundary],
            label=f"other boundary ({int(np.count_nonzero(unknown_boundary))})",
            linewidth=1.0,
        )

    ax.set_title("PRE mesh: Delaunay triangulation and boundary tags")
    ax.set_xlabel("x [nm]")
    ax.set_ylabel("y [nm]")
    ax.legend(loc="best", fontsize=8)
    ax.grid(False)
    fig.tight_layout()
    fig.savefig(output, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    return output


def plot_triangle_area_histogram(
    mesh: MeshData,
    output_path: str | Path,
    *,
    dpi: int = 480,
) -> Path:
    """Plot the distribution of triangle areas in nm^2."""

    output = _prepare_output(output_path)
    areas_nm2 = triangle_areas(mesh.nodes, mesh.triangles) * 1.0e18

    fig, ax = plt.subplots(figsize=(6.2, 4.0))
    ax.hist(areas_nm2, bins=_safe_histogram_bins(areas_nm2, max_bins=60))
    ax.axvline(float(np.mean(areas_nm2)), linestyle="--", linewidth=1.0, label="mean")
    ax.set_title("PRE mesh: triangle area distribution")
    ax.set_xlabel(r"triangle area [nm$^2$]")
    ax.set_ylabel("count")
    ax.legend(loc="best", fontsize=8)
    ax.grid(False)
    fig.tight_layout()
    fig.savefig(output, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    return output


def plot_edge_length_histogram(
    edge_data: EdgeData,
    output_path: str | Path,
    *,
    dpi: int = 480,
) -> Path:
    """Plot the interior/boundary edge-length distributions in nm."""

    output = _prepare_output(output_path)
    lengths_nm = np.asarray(edge_data.lengths, dtype=float) * 1.0e9
    boundary = np.asarray(edge_data.is_boundary, dtype=bool)

    fig, ax = plt.subplots(figsize=(6.2, 4.0))
    if np.any(~boundary):
        interior_lengths = lengths_nm[~boundary]
        ax.hist(
            interior_lengths,
            bins=_safe_histogram_bins(interior_lengths, max_bins=60),
            alpha=0.65,
            label="interior",
        )
    if np.any(boundary):
        boundary_lengths = lengths_nm[boundary]
        ax.hist(
            boundary_lengths,
            bins=_safe_histogram_bins(boundary_lengths, max_bins=60),
            alpha=0.65,
            label="boundary",
        )
    ax.axvline(float(np.mean(lengths_nm)), linestyle="--", linewidth=1.0, label="mean")
    ax.set_title("PRE mesh: edge-length distribution")
    ax.set_xlabel("edge length [nm]")
    ax.set_ylabel("count")
    ax.legend(loc="best", fontsize=8)
    ax.grid(False)
    fig.tight_layout()
    fig.savefig(output, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    return output


def plot_usadel_supercurrent_curve(
    usadel_catalog: Any,
    output_path: str | Path,
    *,
    dpi: int = 480,
) -> Path:
    """Plot the Matsubara Usadel calibration current used by SS runs."""

    output = _prepare_output(output_path)

    q = np.asarray(usadel_catalog.calibration_q_values_m_inv, dtype=float)
    current_uA = 1.0e6 * np.asarray(
        usadel_catalog.calibration_current_values_A,
        dtype=float,
    )
    q_1e7 = q / 1.0e7

    finite = np.isfinite(q_1e7) & np.isfinite(current_uA)
    if not np.any(finite):
        raise ValueError("Usadel calibration current table has no finite points.")

    fig, ax = plt.subplots(figsize=(6.4, 4.0))
    ax.plot(q_1e7[finite], current_uA[finite], marker=".", linewidth=1.0)

    i_max = int(np.nanargmax(current_uA[finite]))
    q_plot = q_1e7[finite]
    i_plot = current_uA[finite]
    ax.plot(q_plot[i_max], i_plot[i_max], marker="o", markersize=5.0, label="model Ic")

    metadata = getattr(usadel_catalog, "metadata", {})
    target = metadata.get("Ic_target_A")
    if target is not None:
        ax.axhline(1.0e6 * float(target), linestyle="--", linewidth=1.0, label="target Ic")

    ax.set_title("Usadel/Matsubara supercurrent calibration")
    ax.set_xlabel(r"q [$10^7$ m$^{-1}$]")
    ax.set_ylabel(r"I$_s$ [$\mu$A]")
    ax.legend(loc="best", fontsize=8)
    ax.grid(False)
    fig.tight_layout()
    fig.savefig(output, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    return output


def _safe_histogram_bins(values: np.ndarray, *, max_bins: int = 60) -> np.ndarray:
    """Return finite histogram edges, including for constant large-valued data.

    NumPy/Matplotlib can fail when all values are identical and very large,
    because automatic range expansion is smaller than floating-point spacing at
    that magnitude. PRE smoke tests intentionally use simple synthetic meshes,
    so equal-area/equal-length distributions must still produce a diagnostic
    figure instead of failing.
    """

    finite = np.asarray(values, dtype=float)
    finite = finite[np.isfinite(finite)]
    if finite.size == 0:
        raise ValueError("Cannot plot histogram: no finite values were provided.")

    vmin = float(np.min(finite))
    vmax = float(np.max(finite))
    scale = max(abs(vmin), abs(vmax), 1.0)
    span = vmax - vmin

    if (not np.isfinite(span)) or span <= 16.0 * np.finfo(float).eps * scale:
        center = float(np.mean(finite))
        pad = max(1.0e-6 * max(abs(center), 1.0), 1.0e-12)
        return np.array([center - pad, center + pad], dtype=float)

    n_bins = int(min(max_bins, max(1, finite.size)))
    return np.linspace(vmin, vmax, n_bins + 1, dtype=float)


def _plot_edge_segments(
    ax: Any,
    nodes_nm: np.ndarray,
    edges: np.ndarray,
    *,
    label: str,
    linewidth: float,
) -> None:
    """Plot many edge segments using the current Matplotlib property cycle."""

    first = True
    for i, j in np.asarray(edges, dtype=np.int64):
        p = nodes_nm[[int(i), int(j)]]
        ax.plot(
            p[:, 0],
            p[:, 1],
            linewidth=linewidth,
            label=label if first else None,
        )
        first = False


def _prepare_output(path: str | Path) -> Path:
    output = Path(path)
    output.parent.mkdir(parents=True, exist_ok=True)
    return output
