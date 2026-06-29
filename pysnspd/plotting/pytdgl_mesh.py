"""pyTDGL-style mesh plots for pySNSPD PRE-runs."""
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from pysnspd.gtdgl.pytdgl_like.finite_volume.mesh import Mesh


def plot_pytdgl_fvm_mesh(mesh: Mesh, output_path: str | Path, *, dpi: int = 480) -> Path:
    """Plot Delaunay and Voronoi meshes using pyTDGL-style options."""
    output = Path(output_path)
    output.parent.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(1, 3, figsize=(12, 3.5), gridspec_kw={"width_ratios": [1.4, 1.0, 1.0]})
    ax, bx, cx = axes
    for a in axes:
        a.set_aspect("equal")
        a.set_xlabel("x [nm]")
        a.set_ylabel("y [nm]")
    x_nm = mesh.sites[:, 0] * 1e9
    y_nm = mesh.sites[:, 1] * 1e9
    ax.triplot(x_nm, y_nm, mesh.elements, lw=0.35, alpha=0.65)
    ax.plot(x_nm, y_nm, ".", ms=1.6, alpha=0.8)
    ax.set_title("sites + Delaunay")

    # Draw dual polygons in nm on the center panel.
    bx.triplot(x_nm, y_nm, mesh.elements, lw=0.25, alpha=0.35)
    if mesh.voronoi_polygons is not None:
        for poly in mesh.voronoi_polygons:
            p = np.asarray(poly) * 1e9
            if len(p):
                q = np.vstack([p, p[0]])
                bx.plot(q[:, 0], q[:, 1], lw=0.35, alpha=0.7)
    bx.plot(x_nm, y_nm, ".", ms=1.2)
    bx.set_title("Delaunay + Voronoi")

    # Zoom near an interior node closest to center, like pyTDGL's py-mesh notebook.
    i0 = mesh.closest_site((float(np.mean(mesh.x)), float(np.mean(mesh.y))))
    d_nm = 4.0 * np.median(mesh.edge_mesh.edge_lengths) * 1e9 if mesh.edge_mesh is not None else 20.0
    cx.triplot(x_nm, y_nm, [tri for tri in mesh.elements if i0 in tri], lw=1.1)
    if mesh.voronoi_polygons is not None:
        p = np.asarray(mesh.voronoi_polygons[i0]) * 1e9
        if len(p):
            q = np.vstack([p, p[0]])
            cx.fill(q[:, 0], q[:, 1], alpha=0.3)
            cx.plot(q[:, 0], q[:, 1], lw=1.1)
    cx.plot(x_nm, y_nm, "ko", ms=2.0)
    cx.plot(x_nm[i0], y_nm[i0], "o", ms=5.0)
    cx.set_xlim(x_nm[i0] - d_nm, x_nm[i0] + d_nm)
    cx.set_ylim(y_nm[i0] - d_nm, y_nm[i0] + d_nm)
    cx.set_title("local Voronoi cell")

    fig.suptitle("pyTDGL-like finite-volume mesh")
    fig.tight_layout()
    fig.savefig(output, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    return output
