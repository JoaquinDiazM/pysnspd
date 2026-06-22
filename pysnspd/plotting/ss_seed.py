"""Plots for OE6 stationary analytic seed diagnostics."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg", force=True)
import matplotlib.pyplot as plt
import matplotlib.tri as mtri


MEV_J = 1.602176634e-22


def plot_ss_seed_delta(mesh, seed, output_path: str | Path, *, dpi: int = 480) -> Path:
    """Plot the analytic seed order-parameter amplitude."""
    delta_meV = seed.node_delta_J / MEV_J
    vmax = max(float(np.nanmax(delta_meV)), 1.0e-30)

    return _plot_node_scalar(
        mesh,
        delta_meV,
        output_path,
        title="OE6 seed: Δ",
        label="Δ [meV]",
        vmin=0.0,
        vmax=vmax,
        dpi=dpi,
    )


def plot_ss_seed_phase(mesh, seed, output_path: str | Path, *, dpi: int = 480) -> Path:
    """Plot the unwrapped seed phase."""
    return _plot_node_scalar(
        mesh,
        seed.node_theta_rad,
        output_path,
        title=r"OE6 seed: unwrapped phase $\theta$",
        label=r"$\theta$ [rad]",
        dpi=dpi,
    )


def plot_ss_seed_current_density(
    mesh,
    seed,
    output_path: str | Path,
    *,
    dpi: int = 480,
) -> Path:
    """Plot current-density magnitude and sparse arrows."""
    output = Path(output_path)
    output.parent.mkdir(parents=True, exist_ok=True)

    nodes = np.asarray(mesh.nodes, dtype=float)
    x_nm = nodes[:, 0] * 1.0e9
    y_nm = nodes[:, 1] * 1.0e9

    tri = mtri.Triangulation(x_nm, y_nm, np.asarray(mesh.triangles, dtype=np.int64))

    jx = np.asarray(seed.node_jtot_x_A_m2, dtype=float)
    jy = np.asarray(seed.node_jtot_y_A_m2, dtype=float)
    jmag = np.sqrt(jx * jx + jy * jy)

    fig, ax = plt.subplots(figsize=(8.0, 3.2))
    vmax = max(float(np.nanmax(jmag)), 1.0e-30)
    im = ax.tripcolor(tri, jmag, shading="gouraud", vmin=0.0, vmax=vmax)
    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label(r"$|\mathbf{j}|$ [A m$^{-2}$]")

    step = max(1, nodes.shape[0] // 150)
    ax.quiver(
        x_nm[::step],
        y_nm[::step],
        jx[::step],
        jy[::step],
        angles="xy",
        scale_units="xy",
        scale=None,
        width=0.002,
    )

    ax.set_title(r"OE6 seed: current density")
    ax.set_xlabel("x [nm]")
    ax.set_ylabel("y [nm]")
    ax.set_aspect("equal", adjustable="box")
    ax.grid(False)

    fig.tight_layout()
    fig.savefig(output, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    return output


def plot_ss_seed_divergence(mesh, seed, output_path: str | Path, *, dpi: int = 480) -> Path:
    """Plot the analytic seed divergence diagnostic."""
    return _plot_node_scalar(
        mesh,
        seed.node_div_j_A_m3,
        output_path,
        title=r"OE6 seed: analytic $\nabla\cdot\mathbf{j}$ diagnostic",
        label=r"$\nabla\cdot\mathbf{j}$ [$\mathrm{A\,m^{-3}}$]",
        dpi=dpi,
    )


def plot_ss_seed_boundary_currents(
    seed,
    output_path: str | Path,
    *,
    dpi: int = 480,
) -> Path:
    """Plot integrated boundary currents."""
    output = Path(output_path)
    output.parent.mkdir(parents=True, exist_ok=True)

    currents = seed.metadata["boundary_currents_A"]
    labels = ["left", "right", "bottom", "top"]
    values = [currents[f"{label}_A"] for label in labels]

    fig, ax = plt.subplots(figsize=(6.0, 3.6))
    ax.bar(labels, values)
    ax.axhline(0.0, linewidth=0.8)
    ax.set_title("OE6 seed: integrated boundary currents")
    ax.set_ylabel("current [A]")
    ax.grid(True, axis="y", linewidth=0.25, alpha=0.35)

    fig.tight_layout()
    fig.savefig(output, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    return output


def _plot_node_scalar(
    mesh,
    values,
    output_path: str | Path,
    *,
    title: str,
    label: str,
    vmin=None,
    vmax=None,
    dpi: int = 480,
) -> Path:
    output = Path(output_path)
    output.parent.mkdir(parents=True, exist_ok=True)

    nodes = np.asarray(mesh.nodes, dtype=float)
    x_nm = nodes[:, 0] * 1.0e9
    y_nm = nodes[:, 1] * 1.0e9

    tri = mtri.Triangulation(
        x_nm,
        y_nm,
        np.asarray(mesh.triangles, dtype=np.int64),
    )

    z = np.asarray(values, dtype=float).reshape(-1)

    fig, ax = plt.subplots(figsize=(8.0, 3.2))
    im = ax.tripcolor(
        tri,
        z,
        shading="gouraud",
        vmin=vmin,
        vmax=vmax,
    )

    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label(label)

    ax.set_title(title)
    ax.set_xlabel("x [nm]")
    ax.set_ylabel("y [nm]")
    ax.set_aspect("equal", adjustable="box")
    ax.grid(False)

    fig.tight_layout()
    fig.savefig(output, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    return output

##