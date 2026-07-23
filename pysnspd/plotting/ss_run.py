"""Adaptive time-step diagnostics for stationary solver runs."""
from __future__ import annotations

from pathlib import Path

import numpy as np

import matplotlib
matplotlib.use("Agg", force=True)

import matplotlib.pyplot as plt

from pysnspd.plotting.style import THESIS_DOUBLE_FIGSIZE, THESIS_DPI, apply_thesis_style

apply_thesis_style()


def plot_ss_adaptive_timestep_history(
    history: dict,
    output_path: str | Path,
    *,
    dpi: int = THESIS_DPI,
) -> Path:
    """Plot adaptive-Euler step-size diagnostics.

    The upper panel shows the actual accepted ``dt`` together with the
    tentative step attempted at the start of each update and the next tentative
    step selected from the moving ``max_d_abs_sq_psi`` window.  The lower panel
    shows how often the pyTDGL-style algebraic |psi|^2 solve had to shrink the
    time step before accepting the update.
    """
    output = Path(output_path)
    output.parent.mkdir(parents=True, exist_ok=True)
    t_ps = np.asarray(history.get("t_s", []), dtype=float) / 1.0e-12

    fig, axes = plt.subplots(
        2,
        1,
        figsize=THESIS_DOUBLE_FIGSIZE,
        sharex=True,
        constrained_layout=False,
    )
    fig.subplots_adjust(left=0.105, right=0.970, bottom=0.090, top=0.925, hspace=0.42)
    ax_dt, ax_retry = axes

    if t_ps.size:
        for key, label in (
            ("dt_attempt_s", r"attempted $\Delta t$"),
            ("dt_s", r"accepted $\Delta t$"),
            ("dt_next_s", r"next tentative $\Delta t$"),
            ("adaptive_target_dt_s", r"window target $\Delta t$"),
        ):
            arr = np.asarray(history.get(key, []), dtype=float)
            if arr.size != t_ps.size:
                continue
            mask = np.isfinite(arr) & (arr > 0.0)
            if np.any(mask):
                ax_dt.plot(t_ps[mask], arr[mask] / 1.0e-15, linewidth=1.2, label=label)

        retries = np.asarray(history.get("adaptive_retries", []), dtype=float)
        rejected = np.asarray(history.get("adaptive_rejected_attempts", []), dtype=float)
        if retries.size == t_ps.size:
            ax_retry.step(t_ps, retries, where="post", linewidth=1.0, label="retries before accepted step")
        if rejected.size == t_ps.size and np.any(rejected):
            cumulative = np.cumsum(np.maximum(rejected, 0.0))
            denom = max(float(np.nanmax(cumulative)), 1.0)
            ax_retry.plot(t_ps, cumulative / denom, linewidth=1.0, label="normalized cumulative rejections")

        mean_d = np.asarray(history.get("adaptive_window_mean_d_abs_sq", []), dtype=float)
        if mean_d.size == t_ps.size and np.any(np.isfinite(mean_d)):
            ax_extra = ax_retry.twinx()
            mask = np.isfinite(mean_d) & (mean_d > 0.0)
            if np.any(mask):
                ax_extra.semilogy(t_ps[mask], mean_d[mask], linewidth=1.0, alpha=0.75, label=r"window mean $\Delta |\psi|^2$")
                ax_extra.set_ylabel(r"window mean $\Delta |\psi|^2$")
                lines_1, labels_1 = ax_retry.get_legend_handles_labels()
                lines_2, labels_2 = ax_extra.get_legend_handles_labels()
                ax_retry.legend(lines_1 + lines_2, labels_1 + labels_2, frameon=False, loc="best")
            else:
                ax_retry.legend(frameon=False, loc="best")
        else:
            ax_retry.legend(frameon=False, loc="best")

    ax_dt.set_title("SS adaptive Euler time-step evolution")
    ax_dt.set_ylabel(r"$\Delta t$ [fs]")
    ax_dt.set_yscale("log")
    ax_dt.grid(False)
    ax_dt.legend(frameon=False, loc="best")

    ax_retry.set_title("pyTDGL-style retry/shrink diagnostics")
    ax_retry.set_xlabel("t [ps]")
    ax_retry.set_ylabel("retries")
    ax_retry.grid(False)

    fig.savefig(output, dpi=dpi, bbox_inches="tight", pad_inches=0.08)
    plt.close(fig)
    return output
