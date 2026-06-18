"""
Basic plotting utilities for pySNSPD.

OE2 scope:
- Mesh geometry plot.
- Boundary tag diagnostic plot.

These functions save figures to disk and do not display interactive windows.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np

import matplotlib
matplotlib.use("Agg", force=True)

import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection

from pysnspd.mesh.delaunay import MeshData
from pysnspd.mesh.edges import EdgeData


def plot_mesh_geometry(
    mesh: MeshData,
    edge_data: EdgeData,
    output_path: str | Path,
    *,
    dpi: int = 180,
) -> Path:
    """
    Plot nodes, triangles and boundary edges.
    """
    output = Path(output_path)
    output.parent.mkdir(parents=True, exist_ok=True)

    nodes_nm = mesh.nodes * 1.0e9

    fig, ax = plt.subplots(figsize=(9.0, 3.0))

    ax.triplot(
        nodes_nm[:, 0],
        nodes_nm[:, 1],
        mesh.triangles,
        linewidth=0.35,
        alpha=0.8,
    )

    ax.scatter(
        nodes_nm[:, 0],
        nodes_nm[:, 1],
        s=3,
        alpha=0.9,
    )

    boundary_edges = edge_data.edges[edge_data.is_boundary]
    if boundary_edges.size:
        segments = nodes_nm[boundary_edges]
        collection = LineCollection(segments, linewidths=1.2)
        ax.add_collection(collection)

    ax.set_title("Delaunay nanowire mesh")
    ax.set_xlabel("x [nm]")
    ax.set_ylabel("y [nm]")
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlim(-0.03 * mesh.length_m * 1.0e9, 1.03 * mesh.length_m * 1.0e9)
    ax.set_ylim(-0.65 * mesh.width_m * 1.0e9, 0.65 * mesh.width_m * 1.0e9)
    ax.grid(True, linewidth=0.25, alpha=0.4)

    fig.tight_layout()
    fig.savefig(output, dpi=dpi)
    plt.close(fig)

    return output


def plot_boundary_tags(
    mesh: MeshData,
    edge_data: EdgeData,
    output_path: str | Path,
    *,
    dpi: int = 180,
) -> Path:
    """
    Plot boundary edges colored by tag.
    """
    output = Path(output_path)
    output.parent.mkdir(parents=True, exist_ok=True)

    nodes_nm = mesh.nodes * 1.0e9

    fig, ax = plt.subplots(figsize=(9.0, 3.0))

    ax.triplot(
        nodes_nm[:, 0],
        nodes_nm[:, 1],
        mesh.triangles,
        linewidth=0.25,
        alpha=0.25,
    )

    tag_styles = {
        "left": {"linewidth": 2.5},
        "right": {"linewidth": 2.5},
        "top": {"linewidth": 1.5},
        "bottom": {"linewidth": 1.5},
        "boundary_unknown": {"linewidth": 2.0},
    }

    for tag, style in tag_styles.items():
        mask = edge_data.tags == tag
        if not np.any(mask):
            continue

        segments = nodes_nm[edge_data.edges[mask]]
        collection = LineCollection(
            segments,
            linewidths=style["linewidth"],
            label=f"{tag} ({int(np.count_nonzero(mask))})",
        )
        ax.add_collection(collection)

    ax.set_title("Boundary tags")
    ax.set_xlabel("x [nm]")
    ax.set_ylabel("y [nm]")
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlim(-0.03 * mesh.length_m * 1.0e9, 1.03 * mesh.length_m * 1.0e9)
    ax.set_ylim(-0.65 * mesh.width_m * 1.0e9, 0.65 * mesh.width_m * 1.0e9)
    ax.grid(True, linewidth=0.25, alpha=0.4)
    ax.legend(loc="upper center", ncol=5, fontsize=7, frameon=True)

    fig.tight_layout()
    fig.savefig(output, dpi=dpi)
    plt.close(fig)

    return output


def plot_stationary_state(config, run_name):
    """Plot stationary gTDGL fields and current-conservation diagnostics."""
    return 0


def plot_photon_transient(config, run_name):
    """Plot photon-run snapshots, histories, and circuit observables."""
    return 0


def plot_catalog_diagnostics(config, run_name):
    """Plot Usadel and phase-space catalog sanity checks."""
    return 0
