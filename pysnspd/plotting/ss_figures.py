"""Presentation-quality figures for stationary SS runs.

The functions in this module expect already-extracted arrays from
``pysnspd.analysis.ss_run``.  They should not open raw ``.npz`` files or
interpret solver metadata; that keeps numerical analysis separate from the
visual methods.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Mapping

import numpy as np
import matplotlib

matplotlib.use("Agg", force=True)
import matplotlib.pyplot as plt
import matplotlib.tri as mtri

from pysnspd.plotting.style import (
    THESIS_DOUBLE_FIGSIZE,
    THESIS_DPI,
    THESIS_WIDTH_IN,
    apply_thesis_style,
)

apply_thesis_style()


def make_ss_run_figures(
    *,
    mesh: Any,
    dataset: Mapping[str, Any],
    output_dir: str | Path,
    dpi: int = THESIS_DPI,
) -> dict[str, Path]:
    """Create the standard SS figure set under ``plots/<run>/figures``."""

    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)
    saved: dict[str, Path] = {}
    saved["overview"] = plot_ss_final_overview(mesh, dataset, out / "ss_final_overview.png", dpi=dpi)
    saved["relaxation"] = plot_ss_relaxation_monitors(dataset, out / "ss_relaxation_monitors.png", dpi=dpi)
    saved["adaptive"] = plot_ss_adaptive_summary(dataset, out / "ss_adaptive_summary.png", dpi=dpi)
    saved["masks"] = plot_ss_region_masks(mesh, dataset, out / "ss_region_masks.png", dpi=dpi)
    return saved

def _legend_if_labels(ax, *, frameon: bool = False, loc: str = "best") -> None:
    """Draw legend only when the axis has visible labeled artists."""
    handles, labels = ax.get_legend_handles_labels()
    filtered = [
        (handle, label)
        for handle, label in zip(handles, labels)
        if label and not label.startswith("_")
    ]
    if filtered:
        handles, labels = zip(*filtered)
        ax.legend(handles, labels, frameon=frameon, loc=loc)

        
def plot_ss_final_overview(
    mesh: Any,
    dataset: Mapping[str, Any],
    output_path: str | Path,
    *,
    dpi: int = THESIS_DPI,
) -> Path:
    """Plot final SS fields in a compact six-panel overview."""

    output = Path(output_path)
    output.parent.mkdir(parents=True, exist_ok=True)
    tri = _triangulation(mesh, dataset)
    jscale = float(dataset.get("javg_A_m2", 1.0))
    jscale = abs(jscale) if np.isfinite(jscale) and abs(jscale) > 0.0 else 1.0

    panels = [
        ("delta_over_delta0", r"$|\Delta|/\Delta_0$", r"$|\Delta|/\Delta_0$", False, 0.0),
        ("phi_mV", r"$\phi$ [mV]", r"$\phi$ [mV]", True, None),
        ("jtot_mag_A_m2", r"$|j|/j_{avg}$", r"$|j|/j_{avg}$", False, 0.0),
        ("js_us_mag_A_m2", r"$|j_s^{Us}|/j_{avg}$", r"$|j_s^{Us}|/j_{avg}$", False, 0.0),
        ("jn_mag_A_m2", r"$|j_n|/j_{avg}$", r"$|j_n|/j_{avg}$", False, 0.0),
        ("pairbreaking_ratio", r"pairbreaking ratio", r"$\chi_{pb}$", False, 0.0),
    ]

    fig, axes = plt.subplots(
        2,
        3,
        figsize=(THESIS_WIDTH_IN, 4.5),
        constrained_layout=False,
    )
    fig.subplots_adjust(
        left=0.065,
        right=0.975,
        bottom=0.085,
        top=0.900,
        wspace=0.58,
        hspace=0.38,
    )
    fig.suptitle(f"SS final state: {dataset.get('run_name', '')}", y=0.975)

    for ax, (key, title, label, symmetric, vmin) in zip(axes.ravel(), panels):
        z = np.asarray(dataset.get(key, []), dtype=float)
        if key.endswith("_A_m2"):
            z = z / jscale
        _draw_node_scalar(
            ax,
            tri,
            z,
            title=title,
            label=label,
            symmetric=symmetric,
            vmin=vmin,
            compact=True,
        )

    fig.savefig(output, dpi=dpi, bbox_inches="tight", pad_inches=0.08)
    plt.close(fig)
    return output


def plot_ss_relaxation_monitors(
    dataset: Mapping[str, Any],
    output_path: str | Path,
    *,
    dpi: int = THESIS_DPI,
) -> Path:
    """Plot physical and numerical monitors versus physical time."""

    output = Path(output_path)
    output.parent.mkdir(parents=True, exist_ok=True)
    t = np.asarray(dataset.get("t_ps", []), dtype=float)
    fig, axes = plt.subplots(
        3,
        1,
        figsize=(THESIS_WIDTH_IN, 6.0),
        sharex=True,
        constrained_layout=False,
    )
    fig.subplots_adjust(left=0.110, right=0.970, bottom=0.085, top=0.930, hspace=0.32)
    fig.suptitle(f"SS relaxation monitors: {dataset.get('run_name', '')}", y=0.985)

    _plot_positive_curve(axes[0], t, dataset.get("eta_R"), r"$\eta_R$", log=True)
    _plot_positive_curve(axes[0], t, dataset.get("pairbreaking_max_history"), r"$\max\chi_{pb}$", log=True)
    axes[0].set_ylabel("stiff monitors")
    axes[0].legend(frameon=False)
    axes[0].grid(False)

    v_tdgl_t = np.asarray(dataset.get("tdgl_probe_voltage_t_ps", []), dtype=float)
    v_tdgl = np.asarray(dataset.get("tdgl_probe_voltage_mV", []), dtype=float)
    if v_tdgl_t.size == 0 or v_tdgl.size == 0:
        v_tdgl_t = t
        v_tdgl = np.asarray(dataset.get("terminal_voltage_mV", []), dtype=float)
    _plot_curve(axes[1], v_tdgl_t, v_tdgl, r"$V_{TDGL}^{\pm 50\,nm}$ [mV]")
    _plot_curve(axes[1], t, dataset.get("normal_current_fraction"), r"$\max|j_n|/\max|j|$")
    axes[1].set_ylabel("physical monitors")
    axes[1].legend(frameon=False)
    axes[1].grid(False)

    _plot_positive_curve(axes[2], t, dataset.get("dt_accepted_fs"), r"accepted $\Delta t$ [fs]", log=True)
    _plot_positive_curve(axes[2], t, dataset.get("dt_next_fs"), r"next tentative $\Delta t$ [fs]", log=True)
    axes[2].set_xlabel("t [ps]")
    axes[2].set_ylabel(r"$\Delta t$ [fs]")
    axes[2].legend(frameon=False)
    axes[2].grid(False)

    fig.savefig(output, dpi=dpi, bbox_inches="tight", pad_inches=0.08)
    plt.close(fig)
    return output



def plot_ss_adaptive_summary(
    dataset: Mapping[str, Any],
    output_path: str | Path,
    *,
    dpi: int = THESIS_DPI,
) -> Path:
    """Plot adaptive Euler information from relaxation history."""

    output = Path(output_path)
    output.parent.mkdir(parents=True, exist_ok=True)
    t = np.asarray(dataset.get("t_ps", []), dtype=float)
    fig, axes = plt.subplots(
        2,
        1,
        figsize=(THESIS_WIDTH_IN, 4.8),
        sharex=True,
        constrained_layout=False,
    )
    fig.subplots_adjust(left=0.110, right=0.970, bottom=0.095, top=0.920, hspace=0.34)
    fig.suptitle(f"SS adaptive Euler summary: {dataset.get('run_name', '')}", y=0.985)

    _plot_positive_curve(axes[0], t, dataset.get("dt_attempt_fs"), r"attempted $\Delta t$", log=True)
    _plot_positive_curve(axes[0], t, dataset.get("dt_accepted_fs"), r"accepted $\Delta t$", log=True)
    _plot_positive_curve(axes[0], t, dataset.get("dt_next_fs"), r"next tentative $\Delta t$", log=True)
    _plot_positive_curve(axes[0], t, dataset.get("adaptive_target_dt_fs"), r"window target $\Delta t$", log=True)
    axes[0].set_ylabel(r"$\Delta t$ [fs]")
    axes[0].legend(frameon=False)
    axes[0].grid(False)

    _plot_curve(axes[1], t, dataset.get("adaptive_retries"), "retries")
    _plot_curve(axes[1], t, dataset.get("adaptive_rejected_attempts"), "cumulative rejected attempts")
    axes2 = axes[1].twinx()
    _plot_positive_curve(axes2, t, dataset.get("adaptive_window_mean_d_abs_sq"), r"window mean $\Delta|\psi|^2$", log=True)
    axes[1].set_xlabel("t [ps]")
    axes[1].set_ylabel("retry count")
    axes2.set_ylabel(r"window mean $\Delta|\psi|^2$")
    _legend_if_labels(axes[1], frameon=False, loc="upper left")
    _legend_if_labels(axes2, frameon=False, loc="upper right")
    axes[1].grid(False)

    fig.savefig(output, dpi=dpi, bbox_inches="tight", pad_inches=0.08)
    plt.close(fig)
    return output


def plot_ss_region_masks(
    mesh: Any,
    dataset: Mapping[str, Any],
    output_path: str | Path,
    *,
    dpi: int = THESIS_DPI,
) -> Path:
    """Visualize terminal and bulk masks used by post-run analysis."""

    output = Path(output_path)
    output.parent.mkdir(parents=True, exist_ok=True)
    tri = _triangulation(mesh, dataset)
    terminal = np.asarray(dataset.get("normal_terminal_node_mask", []), dtype=bool)
    bulk = np.asarray(dataset.get("bulk_node_mask", []), dtype=bool)
    if terminal.size != bulk.size:
        n = np.asarray(dataset.get("x_nm", [])).size
        terminal = np.resize(terminal, n).astype(bool)
        bulk = np.resize(bulk, n).astype(bool)
    z = np.zeros(bulk.size, dtype=float)
    z[bulk] = 1.0
    z[terminal] = 2.0

    fig, ax = plt.subplots(figsize=THESIS_DOUBLE_FIGSIZE)
    _draw_node_scalar(ax, tri, z, title="analysis masks: bulk and metallic terminals", label="mask id", vmin=0.0)
    fig.savefig(output, dpi=dpi, bbox_inches="tight", pad_inches=0.08)
    plt.close(fig)
    return output


def _triangulation(mesh: Any, dataset: Mapping[str, Any]) -> mtri.Triangulation:
    x = np.asarray(dataset.get("x_nm", np.asarray(mesh.nodes)[:, 0] * 1.0e9), dtype=float)
    y = np.asarray(dataset.get("y_nm", np.asarray(mesh.nodes)[:, 1] * 1.0e9), dtype=float)
    triangles = np.asarray(dataset.get("triangles", mesh.triangles), dtype=np.int64)
    return mtri.Triangulation(x, y, triangles)


def _draw_node_scalar(
    ax,
    tri: mtri.Triangulation,
    values: np.ndarray,
    *,
    title: str,
    label: str,
    symmetric: bool = False,
    vmin: float | None = None,
    compact: bool = False,
):
    z = np.asarray(values, dtype=float)
    if z.size != tri.x.size:
        z = np.resize(z, tri.x.size)
    finite = z[np.isfinite(z)]
    if finite.size == 0:
        finite = np.array([0.0])
    if symmetric:
        vmax = float(np.nanpercentile(np.abs(finite), 99.5))
        vmax = max(vmax, 1.0e-30)
        vmin = -vmax
    else:
        vmax = float(np.nanpercentile(finite, 99.5))
        if vmin is None:
            vmin = float(np.nanpercentile(finite, 0.5))
        if not np.isfinite(vmax) or vmax <= vmin:
            vmax = vmin + 1.0
    mappable = ax.tripcolor(tri, z, shading="gouraud", vmin=vmin, vmax=vmax)
    if compact:
        ax.set_title(title, fontsize=8.0)
        ax.set_xlabel("x [nm]", fontsize=8.0)
        ax.set_ylabel("y [nm]", fontsize=8.0)
        ax.tick_params(labelsize=7.0)
    else:
        ax.set_title(title)
        ax.set_xlabel("x [nm]")
        ax.set_ylabel("y [nm]")
    ax.set_aspect("equal", adjustable="box")
    ax.grid(False)
    cb = ax.figure.colorbar(mappable, ax=ax, shrink=0.86)
    if compact:
        cb.set_label(label, fontsize=8.0)
        cb.ax.tick_params(labelsize=7.0)
    else:
        cb.set_label(label)


def _plot_curve(ax, x, y, label: str):
    x_arr = np.asarray(x, dtype=float)
    y_arr = np.asarray(y if y is not None else [], dtype=float)
    if x_arr.size == 0 or y_arr.size == 0:
        return
    if y_arr.size != x_arr.size:
        y_arr = np.resize(y_arr, x_arr.size)
    ax.plot(x_arr, y_arr, label=label)


def _plot_positive_curve(ax, x, y, label: str, *, log: bool = False):
    x_arr = np.asarray(x, dtype=float)
    y_arr = np.asarray(y if y is not None else [], dtype=float)
    if x_arr.size == 0 or y_arr.size == 0:
        return
    if y_arr.size != x_arr.size:
        y_arr = np.resize(y_arr, x_arr.size)
    mask = np.isfinite(x_arr) & np.isfinite(y_arr)
    if log:
        mask &= y_arr > 0.0
    if not np.any(mask):
        return
    ax.plot(x_arr[mask], y_arr[mask], label=label)
    if log:
        ax.set_yscale("log")
