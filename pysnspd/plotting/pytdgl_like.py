"""Native diagnostic plots for the OE7 pyTDGL-like backend.

These plots intentionally inspect the pyTDGL-like sparse system directly,
without forcing the arrays through the legacy OE7 SI-current plotting adapter.
They are meant to answer one question first: does the pyTDGL-like linear system
close internally, or is the mismatch introduced by the SI adapter/diagnostics?
"""
from __future__ import annotations

from pathlib import Path

import numpy as np

try:  # pragma: no cover - plotting is optional in headless tests.
    import matplotlib.pyplot as plt
    import matplotlib.tri as mtri
    from matplotlib.collections import LineCollection
except Exception:  # pragma: no cover
    plt = None
    mtri = None
    LineCollection = None


def plot_pytdgl_like_native_history(history: dict, output_path: str | Path, *, dpi: int = 480) -> Path:
    """Plot native pyTDGL-like scalar diagnostics stored by the solver."""

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    if plt is None:
        output_path.touch()
        return output_path

    t_ps = np.asarray(history.get("t_s", []), dtype=float) / 1.0e-12
    if t_ps.size == 0:
        t_ps = np.arange(1, dtype=float)

    def y(key: str) -> np.ndarray:
        arr = np.asarray(history.get(key, np.zeros_like(t_ps)), dtype=float).reshape(-1)
        if arr.size != t_ps.size:
            arr = np.resize(arr, t_ps.size)
        return arr

    fig, axes = plt.subplots(3, 1, figsize=(12, 10), constrained_layout=True)
    ax = axes[0]
    ax.semilogy(t_ps, np.maximum(y("pytdgl_like_poisson_residual_rel"), 1.0e-300), label="||Lμ-rhs|| / ||rhs||")
    ax.semilogy(t_ps, np.maximum(y("pytdgl_like_poisson_residual_max_abs"), 1.0e-300), label="max |Lμ-rhs|")
    ax.semilogy(t_ps, np.maximum(y("pytdgl_like_mu_boundary_max_abs"), 1.0e-300), label="max |μ boundary|")
    ax.set_title("pyTDGL-like native Poisson diagnostics")
    ax.set_xlabel("t [ps]")
    ax.set_ylabel("native value")
    ax.legend(loc="best")

    ax = axes[1]
    ax.semilogy(t_ps, np.maximum(y("pytdgl_like_div_supercurrent_norm"), 1.0e-300), label="||div J_s||")
    ax.semilogy(t_ps, np.maximum(y("pytdgl_like_boundary_rhs_norm"), 1.0e-300), label="||boundary rhs||")
    ax.semilogy(t_ps, np.maximum(y("pytdgl_like_poisson_rhs_norm"), 1.0e-300), label="||rhs||")
    ax.set_title("native RHS balance")
    ax.set_xlabel("t [ps]")
    ax.set_ylabel("native norm")
    ax.legend(loc="best")

    ax = axes[2]
    ax.semilogy(t_ps, np.maximum(y("pytdgl_like_native_si_residual_no_boundary_rms_A_m3"), 1.0e-300), label="SI rms div(j), no boundary")
    ax.semilogy(t_ps, np.maximum(y("pytdgl_like_native_si_residual_plus_boundary_rms_A_m3"), 1.0e-300), label="SI rms div(j)+boundary")
    ax.semilogy(t_ps, np.maximum(y("pytdgl_like_native_si_residual_minus_boundary_rms_A_m3"), 1.0e-300), label="SI rms div(j)-boundary")
    ax.set_title("native currents converted back to SI")
    ax.set_xlabel("t [ps]")
    ax.set_ylabel("rms residual [A m$^{-3}$]")
    ax.legend(loc="best")

    fig.savefig(output_path, dpi=dpi)
    plt.close(fig)
    return output_path


def plot_pytdgl_like_poisson_snapshots(mesh, history: dict, output_path: str | Path, *, dpi: int = 480, ncols: int = 3) -> Path:
    """Plot node-wise native Poisson terms at stored trajectory snapshots."""

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    if plt is None or mtri is None:
        output_path.touch()
        return output_path

    fields = [
        ("pytdgl_like_div_supercurrent_snapshot", "div J_s"),
        ("pytdgl_like_boundary_rhs_snapshot", "boundary rhs"),
        ("pytdgl_like_poisson_residual_snapshot", "Lμ - rhs"),
    ]
    t_s = np.asarray(history.get("pytdgl_like_snapshot_t_s", history.get("snapshot_t_s", [])), dtype=float)
    nodes = np.asarray(mesh.nodes, dtype=float)
    triangles = np.asarray(mesh.triangles, dtype=int)
    tri = mtri.Triangulation(nodes[:, 0] * 1.0e9, nodes[:, 1] * 1.0e9, triangles)

    ns = int(max([np.asarray(history.get(key, [])).shape[0] if np.asarray(history.get(key, [])).ndim == 2 else 0 for key, _ in fields] + [0]))
    if ns == 0:
        output_path.touch()
        return output_path
    cols = max(1, int(ncols))
    rows = len(fields)
    fig, axes = plt.subplots(rows, cols, figsize=(4.8 * cols, 3.4 * rows), constrained_layout=True)
    axes = np.asarray(axes).reshape(rows, cols)

    if t_s.size < cols:
        t_s = np.resize(t_s, cols)
    # Use evenly spaced stored frames if there are more frames than columns.
    idxs = np.unique(np.linspace(0, ns - 1, cols).astype(int))
    if idxs.size < cols:
        idxs = np.resize(idxs, cols)

    for r, (key, label) in enumerate(fields):
        arr = np.asarray(history.get(key), dtype=float)
        vmax = float(np.nanmax(np.abs(arr))) if arr.size else 1.0
        vmax = max(vmax, 1.0e-300)
        for c, idx in enumerate(idxs):
            ax = axes[r, c]
            values = arr[min(int(idx), arr.shape[0] - 1)]
            im = ax.tripcolor(tri, values, shading="gouraud", vmin=-vmax, vmax=vmax, cmap="coolwarm")
            ax.set_aspect("equal", adjustable="box")
            ax.set_xlabel("x [nm]")
            ax.set_ylabel("y [nm]")
            if r == 0:
                ax.set_title(f"t = {t_s[min(int(idx), t_s.size - 1)] / 1e-12:.4g} ps")
            if c == 0:
                ax.text(-0.18, 0.5, label, rotation=90, va="center", ha="center", transform=ax.transAxes)
            fig.colorbar(im, ax=ax, fraction=0.046, pad=0.02)

    fig.suptitle("pyTDGL-like native Poisson terms")
    fig.savefig(output_path, dpi=dpi)
    plt.close(fig)
    return output_path


def plot_pytdgl_like_native_edge_currents(mesh, history: dict, output_path: str | Path, *, dpi: int = 480) -> Path:
    """Plot final native edge super/normal/total currents on the mesh."""

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    if plt is None or LineCollection is None:
        output_path.touch()
        return output_path

    nodes = np.asarray(mesh.nodes, dtype=float)[:, :2] * 1.0e9
    edges = np.column_stack([
        np.asarray(history.get("edge_i"), dtype=int),
        np.asarray(history.get("edge_j"), dtype=int),
    ])
    if edges.size == 0:
        output_path.touch()
        return output_path
    segments = np.stack([nodes[edges[:, 0]], nodes[edges[:, 1]]], axis=1)
    datasets = [
        ("pytdgl_like_native_supercurrent_snapshot", "native J_s"),
        ("pytdgl_like_native_normal_current_snapshot", "native J_n"),
        ("pytdgl_like_native_total_current_snapshot", "native J_total"),
    ]
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.5), constrained_layout=True)
    for ax, (key, title) in zip(axes, datasets):
        arr = np.asarray(history.get(key, np.zeros((1, edges.shape[0]))), dtype=float)
        vals = arr[-1] if arr.ndim == 2 and arr.shape[1] == edges.shape[0] else np.zeros(edges.shape[0])
        vmax = max(float(np.nanmax(np.abs(vals))), 1.0e-300)
        lc = LineCollection(segments, array=vals, cmap="coolwarm", linewidths=1.0)
        lc.set_clim(-vmax, vmax)
        ax.add_collection(lc)
        ax.autoscale()
        ax.set_aspect("equal", adjustable="box")
        ax.set_title(title)
        ax.set_xlabel("x [nm]")
        ax.set_ylabel("y [nm]")
        fig.colorbar(lc, ax=ax, fraction=0.046, pad=0.02)
    fig.suptitle("pyTDGL-like native final edge currents")
    fig.savefig(output_path, dpi=dpi)
    plt.close(fig)
    return output_path


def plot_pytdgl_like_usadel_gl_comparison(history: dict, output_path: str | Path, *, dpi: int = 480) -> Path:
    """Plot diagnostic Usadel-vs-GL supercurrent mismatch histories."""

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    if plt is None:
        output_path.touch()
        return output_path

    available = np.asarray(history.get("usadel_current_available", [False]), dtype=bool).reshape(-1)
    if available.size and not bool(available[0]):
        output_path.touch()
        return output_path

    t_ps = np.asarray(history.get("snapshot_t_s", history.get("t_s", [])), dtype=float) / 1.0e-12
    rel = np.asarray(history.get("usadel_vs_gl_edge_relative_l2", []), dtype=float).reshape(-1)
    diff = np.asarray(history.get("usadel_vs_gl_edge_max_abs_diff_A_m2", []), dtype=float).reshape(-1)
    usmax = np.asarray(history.get("usadel_supercurrent_max_A_m2", []), dtype=float).reshape(-1)
    glmax = np.asarray(history.get("gl_supercurrent_max_A_m2", []), dtype=float).reshape(-1)
    if rel.size == 0:
        output_path.touch()
        return output_path
    if t_ps.size != rel.size:
        t_ps = np.linspace(0.0, float(t_ps[-1]) if t_ps.size else 0.0, rel.size)

    fig, axes = plt.subplots(2, 1, figsize=(10, 7), constrained_layout=True)
    axes[0].semilogy(t_ps, np.maximum(rel, 1.0e-300), label=r"$||j_s^{Us}-j_s^{GL}||/||j_s^{Us}||$")
    axes[0].set_title("Usadel vs GL supercurrent diagnostic")
    axes[0].set_xlabel("t [ps]")
    axes[0].set_ylabel("relative value")
    axes[0].legend(loc="best")

    axes[1].semilogy(t_ps, np.maximum(usmax, 1.0e-300), label=r"max $|j_s^{Usadel}|$")
    axes[1].semilogy(t_ps, np.maximum(glmax, 1.0e-300), label=r"max $|j_s^{GL}|$")
    axes[1].semilogy(t_ps, np.maximum(diff, 1.0e-300), label=r"max $|j_s^{Usadel}-j_s^{GL}|$")
    axes[1].set_xlabel("t [ps]")
    axes[1].set_ylabel(r"current density [A m$^{-2}$]")
    axes[1].legend(loc="best")

    fig.savefig(output_path, dpi=dpi)
    plt.close(fig)
    return output_path
