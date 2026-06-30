"""Plots for Appendix-B Allmaras diagnostics."""
from __future__ import annotations

from pathlib import Path

import numpy as np

import matplotlib
matplotlib.use("Agg", force=True)

import matplotlib.pyplot as plt
import matplotlib.tri as mtri


def plot_allmaras_appendix_b_diagnostics(mesh, history: dict, output_dir: str | Path, *, dpi: int = 480, ncols: int = 3) -> dict[str, Path]:
    """Plot the Allmaras Appendix-B diagnostics if present."""
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)
    saved: dict[str, Path] = {}
    specs = [
        (
            "mismatch_divergence",
            "allmaras_mismatch_divergence_snapshot_A_m3",
            "allmaras_mismatch_snapshot_t_s",
            out / "allmaras_mismatch_divergence_snapshots.png",
            r"$\nabla\cdot(j_s^{Us}-j_s^{GL})$ [A m$^{-3}$]",
            True,
        ),
        (
            "phase_drive",
            "allmaras_phase_drive_abs_over_delta0_snapshot",
            "allmaras_phase_drive_snapshot_t_s",
            out / "allmaras_phase_drive_snapshots.png",
            r"$|S_{All}|/\Delta_0$",
            False,
        ),
        (
            "delta_mod",
            "allmaras_delta_mod_over_delta0_snapshot",
            "allmaras_snapshot_t_s",
            out / "allmaras_delta_mod_snapshots.png",
            r"$\Delta_{mod}(T_e)/\Delta_0$",
            False,
        ),
        (
            "rho_kwt",
            "allmaras_rho_kwt_snapshot",
            "allmaras_snapshot_t_s",
            out / "allmaras_rho_kwt_snapshots.png",
            r"$\rho_{KWT}$",
            False,
        ),
    ]
    for name, key, tkey, path, label, symmetric in specs:
        if key not in history:
            continue
        arr = np.asarray(history[key], dtype=float)
        if arr.ndim != 2:
            continue
        t_s = np.asarray(history.get(tkey, history.get("snapshot_t_s", np.arange(arr.shape[0]))), dtype=float)
        saved[name] = _plot_snapshot_grid(
            mesh,
            arr,
            t_s,
            path,
            title=f"Appendix-B Allmaras diagnostic: {name.replace('_', ' ')}",
            label=label,
            symmetric=symmetric,
            dpi=dpi,
            ncols=ncols,
        )
    return saved


def _plot_snapshot_grid(mesh, values, t_s, output_path, *, title, label, symmetric=False, dpi=480, ncols=3):
    nodes = np.asarray(mesh.nodes, dtype=float)
    tri = np.asarray(mesh.triangles, dtype=np.int64)
    z = np.asarray(values, dtype=float)
    if z.ndim != 2 or z.shape[1] != nodes.shape[0]:
        raise ValueError(f"values must have shape (n_snapshots, n_nodes), got {z.shape}.")
    n = z.shape[0]
    ncols = max(1, int(ncols))
    nrows = int(np.ceil(n / ncols))
    triang = mtri.Triangulation(nodes[:, 0] * 1e9, nodes[:, 1] * 1e9, tri)

    finite = z[np.isfinite(z)]
    if finite.size:
        if symmetric:
            vmax = max(float(np.nanpercentile(np.abs(finite), 99.0)), 1.0e-300)
            vmin = -vmax
        else:
            vmin = 0.0 if np.nanmin(finite) >= 0 else float(np.nanpercentile(finite, 1.0))
            vmax = max(float(np.nanpercentile(finite, 99.0)), vmin + 1.0e-300)
    else:
        vmin, vmax = (-1.0, 1.0) if symmetric else (0.0, 1.0)

    fig, axes = plt.subplots(nrows, ncols, figsize=(4.4 * ncols + 0.8, 3.3 * nrows), squeeze=False)
    mappable = None
    for k, ax in enumerate(axes.ravel()):
        if k >= n:
            ax.axis("off")
            continue
        mappable = ax.tripcolor(triang, z[k], shading="gouraud", vmin=vmin, vmax=vmax)
        ax.set_aspect("equal", adjustable="box")
        ax.set_xlabel("x [nm]")
        ax.set_ylabel("y [nm]")
        tt = float(t_s[k]) / 1e-12 if k < len(t_s) else float(k)
        ax.set_title(f"t = {tt:.5g} ps")
    if mappable is not None:
        cbar = fig.colorbar(mappable, ax=axes.ravel().tolist(), shrink=0.86, pad=0.01)
        cbar.set_label(label)
    fig.suptitle(title)
    output = Path(output_path)
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=dpi)
    plt.close(fig)
    return output
