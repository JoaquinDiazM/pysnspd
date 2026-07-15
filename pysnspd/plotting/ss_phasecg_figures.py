"""Focused figures for the normalized Allmaras phase-drive SS runs."""

from __future__ import annotations

from pathlib import Path
from typing import Any, Mapping

import matplotlib

matplotlib.use("Agg", force=True)
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import numpy as np


def make_phasecg_ss_figures(
    *,
    dataset: Mapping[str, Any],
    output_dir: str | Path,
    dpi: int = 240,
) -> dict[str, Path]:
    """Create the three E2 diagnostics requested for one completed SS run."""

    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)
    return {
        "snapshot_fields": plot_phasecg_snapshot_fields(
            dataset,
            out / "E2_phasecg_snapshot_fields.png",
            dpi=dpi,
        ),
        "physical_evolution": plot_phasecg_physical_evolution(
            dataset,
            out / "E2_phasecg_physical_evolution.png",
            dpi=dpi,
        ),
        "numerical_diagnostics": plot_phasecg_numerical_diagnostics(
            dataset,
            out / "E2_phasecg_numerical_diagnostics.png",
            dpi=dpi,
        ),
    }


def plot_phasecg_snapshot_fields(
    dataset: Mapping[str, Any],
    output_path: str | Path,
    *,
    dpi: int = 240,
) -> Path:
    """Plot the fundamental mesoscopic fields at every stored snapshot."""

    output = _prepare_output(output_path)
    times = np.asarray(dataset.get("snapshot_t_ps", []), dtype=float)
    if times.size == 0:
        raise ValueError("No snapshot times are available for the E2 field figure.")

    tri = mtri.Triangulation(
        np.asarray(dataset["nodes_x_nm"], dtype=float),
        np.asarray(dataset["nodes_y_nm"], dtype=float),
        np.asarray(dataset["triangles"], dtype=np.int64),
    )
    fields = [
        (
            np.asarray(dataset["delta_snapshot_over_delta0"], dtype=float),
            r"$|\Delta|/\Delta_{\mathrm{BCS}}(0)$",
            "viridis",
            False,
            0.0,
            1.0,
        ),
        (
            np.asarray(dataset["phi_snapshot_mV"], dtype=float),
            r"$\phi$ [mV]",
            "coolwarm",
            True,
            None,
            None,
        ),
        (
            np.asarray(dataset["qxi_snapshot"], dtype=float),
            r"$|\mathbf{q}|\,\xi$",
            "magma",
            False,
            0.0,
            None,
        ),
        (
            np.asarray(dataset["js_snapshot_over_javg"], dtype=float),
            r"$|\mathbf{j}_s^{\mathrm{Us}}|/j_{\mathrm{avg}}$",
            "viridis",
            False,
            0.0,
            None,
        ),
        (
            np.asarray(dataset["jn_snapshot_over_javg"], dtype=float),
            r"$|\mathbf{j}_n|/j_{\mathrm{avg}}$",
            "plasma",
            False,
            0.0,
            None,
        ),
        (
            np.asarray(dataset["jtot_snapshot_over_javg"], dtype=float),
            r"$|\mathbf{j}|/j_{\mathrm{avg}}$",
            "viridis",
            False,
            0.0,
            None,
        ),
        (
            np.asarray(dataset["div_j_snapshot_normalized"], dtype=float),
            r"$\xi\,\nabla\!\cdot\!\mathbf{j}/j_{\mathrm{avg}}$",
            "coolwarm",
            True,
            None,
            None,
        ),
    ]

    n_rows = times.size
    n_cols = len(fields)
    fig, axes = plt.subplots(
        n_rows,
        n_cols,
        figsize=(18.2, max(5.0, 1.15 * n_rows + 1.7)),
        squeeze=False,
    )
    fig.subplots_adjust(left=0.045, right=0.965, bottom=0.045, top=0.925, wspace=0.10, hspace=0.12)
    fig.suptitle(
        f"Current-driven SS snapshots, thermal dynamics disabled\n{dataset.get('run_name', '')}",
        y=0.988,
        fontsize=14,
    )

    column_mappables = []
    for col, (values, label, cmap, symmetric, forced_min, forced_max) in enumerate(fields):
        vmin, vmax = _global_limits(
            values,
            symmetric=symmetric,
            forced_min=forced_min,
            forced_max=forced_max,
        )
        for row in range(n_rows):
            ax = axes[row, col]
            mappable = ax.tripcolor(
                tri,
                values[row],
                shading="gouraud",
                cmap=cmap,
                vmin=vmin,
                vmax=vmax,
                rasterized=True,
            )
            ax.set_aspect("equal", adjustable="box")
            ax.set_xlim(float(np.nanmin(tri.x)), float(np.nanmax(tri.x)))
            ax.set_ylim(float(np.nanmin(tri.y)), float(np.nanmax(tri.y)))
            ax.grid(False)
            if row < n_rows - 1:
                ax.set_xticklabels([])
            else:
                ax.set_xlabel("x [nm]")
            if col == 0:
                ax.set_ylabel("y [nm]")
            else:
                ax.set_yticklabels([])
            if col == n_cols - 1:
                ax.text(
                    1.04,
                    0.5,
                    rf"$t={times[row]:.3g}$ ps",
                    transform=ax.transAxes,
                    rotation=-90,
                    va="center",
                    ha="left",
                    fontsize=9,
                )
        column_mappables.append((mappable, label))

    for col, (mappable, label) in enumerate(column_mappables):
        position = axes[0, col].get_position()
        colorbar_axis = fig.add_axes([position.x0, position.y1 + 0.010, position.width, 0.009])
        colorbar = fig.colorbar(mappable, cax=colorbar_axis, orientation="horizontal")
        colorbar_axis.xaxis.set_ticks_position("top")
        colorbar_axis.xaxis.set_label_position("top")
        colorbar.set_label(label, labelpad=2)
        colorbar_axis.tick_params(labelsize=8, pad=1)

    fig.savefig(output, dpi=dpi, bbox_inches="tight", pad_inches=0.08)
    plt.close(fig)
    return output


def plot_phasecg_physical_evolution(
    dataset: Mapping[str, Any],
    output_path: str | Path,
    *,
    dpi: int = 240,
) -> Path:
    """Plot current partition, voltage, and condensate response."""

    output = _prepare_output(output_path)
    snapshot_t = np.asarray(dataset.get("snapshot_t_ps", []), dtype=float)
    history_t = np.asarray(dataset.get("t_ps", []), dtype=float)
    fig, axes = plt.subplots(3, 1, figsize=(10.4, 9.2), sharex=True)
    fig.subplots_adjust(left=0.105, right=0.885, bottom=0.080, top=0.915, hspace=0.25)
    target_current = float(dataset.get("target_current_uA", np.nan))
    fig.suptitle(
        f"Physical response of the {target_current:.3g} \N{MICRO SIGN}A thermal-off SS run\n"
        f"{dataset.get('run_name', '')}",
        y=0.982,
        fontsize=14,
    )

    axes[0].axhline(
        float(dataset.get("target_current_uA", np.nan)),
        color="black",
        linestyle="--",
        linewidth=1.5,
        label=r"target $I$",
    )
    axes[0].plot(snapshot_t, dataset.get("current_total_snapshot_uA"), "o-", label=r"$I_{\mathrm{tot}}$")
    axes[0].plot(snapshot_t, dataset.get("current_super_snapshot_uA"), "o-", label=r"$I_s^{\mathrm{Us}}$")
    axes[0].plot(snapshot_t, dataset.get("current_normal_snapshot_uA"), "o-", label=r"$I_n$")
    axes[0].set_ylabel("current [\N{MICRO SIGN}A]")
    axes[0].legend(frameon=False, ncol=4, loc="best")

    terminal_t, terminal_v = _decimate(history_t, dataset.get("terminal_voltage_mV"), max_points=7000)
    axes[1].plot(terminal_t, terminal_v, color="tab:blue", linewidth=1.1, label=r"terminal $V$")
    axes[1].plot(
        snapshot_t,
        dataset.get("voltage_center_snapshot_mV"),
        "o-",
        color="tab:purple",
        label=r"central 100 nm $V$",
    )
    axes[1].plot(
        snapshot_t,
        dataset.get("voltage_terminal_snapshot_mV"),
        "o",
        color="tab:cyan",
        label="terminal snapshots",
    )
    axes[1].set_ylabel("voltage [mV]")
    axes[1].legend(frameon=False, ncol=3, loc="best")

    axes[2].plot(
        snapshot_t,
        dataset.get("delta_center_min"),
        "o-",
        label=r"central min $|\Delta|/\Delta_{\mathrm{BCS}}(0)$",
    )
    axes[2].plot(
        snapshot_t,
        dataset.get("delta_center_mean"),
        "o-",
        label=r"central mean $|\Delta|/\Delta_{\mathrm{BCS}}(0)$",
    )
    axes[2].plot(
        snapshot_t,
        dataset.get("delta_center_max"),
        "o-",
        label=r"central max $|\Delta|/\Delta_{\mathrm{BCS}}(0)$",
    )
    axes[2].set_ylabel(r"$|\Delta|/\Delta_{\mathrm{BCS}}(0)$")
    axes[2].set_xlabel("t [ps]")
    fraction_axis = axes[2].twinx()
    fraction_axis.plot(
        snapshot_t,
        dataset.get("normal_current_fraction_snapshot"),
        "s--",
        color="tab:red",
        label=r"$|I_n/I_{\mathrm{tot}}|$",
    )
    fraction_axis.set_ylabel("normal-current fraction", color="tab:red")
    fraction_axis.tick_params(axis="y", colors="tab:red")
    _combined_legend(axes[2], fraction_axis, ncol=2)

    for ax in axes:
        ax.grid(True, alpha=0.22)
        ax.set_xlim(left=0.0)

    fig.savefig(output, dpi=dpi, bbox_inches="tight", pad_inches=0.08)
    plt.close(fig)
    return output


def plot_phasecg_numerical_diagnostics(
    dataset: Mapping[str, Any],
    output_path: str | Path,
    *,
    dpi: int = 240,
) -> Path:
    """Plot adaptive, continuity, and Allmaras-CG implementation metrics."""

    output = _prepare_output(output_path)
    t = np.asarray(dataset.get("t_ps", []), dtype=float)
    snap_t = np.asarray(dataset.get("snapshot_t_ps", []), dtype=float)
    fig, axes = plt.subplots(3, 3, figsize=(16.5, 11.0))
    fig.subplots_adjust(
        left=0.065,
        right=0.950,
        bottom=0.075,
        top=0.900,
        wspace=0.42,
        hspace=0.34,
    )
    fig.suptitle(
        f"Numerical diagnostics: normalized phase drive and adaptive solve\n{dataset.get('run_name', '')}",
        y=0.978,
        fontsize=14,
    )

    ax = axes[0, 0]
    _plot_decimated(ax, t, dataset.get("dt_attempt_fs"), "attempted", positive=True)
    _plot_decimated(ax, t, dataset.get("dt_accepted_fs"), "accepted", positive=True)
    _plot_decimated(ax, t, dataset.get("dt_next_fs"), "next", positive=True)
    ax.set_yscale("log")
    ax.set_ylabel(r"$\Delta t$ [fs]")
    ax.set_title("adaptive time step")
    ax.legend(frameon=False, fontsize=8)

    ax = axes[0, 1]
    _plot_decimated(ax, t, dataset.get("solve_attempts_per_step"), "solve attempts / accepted step")
    _plot_decimated(ax, t, dataset.get("adaptive_retries"), "retries")
    rejected_axis = ax.twinx()
    _plot_decimated(
        rejected_axis,
        t,
        dataset.get("cumulative_rejected_attempts"),
        "cumulative rejected",
        color="tab:red",
    )
    ax.set_ylabel("count per accepted step")
    rejected_axis.set_ylabel("cumulative rejected", color="tab:red")
    rejected_axis.tick_params(axis="y", colors="tab:red")
    ax.set_title("nonlinear solve effort")
    _combined_legend(ax, rejected_axis, fontsize=8)

    ax = axes[0, 2]
    estimated = np.asarray(dataset.get("estimated_wall_step_s", []), dtype=float)
    cumulative = np.asarray(dataset.get("estimated_wall_cumulative_s", []), dtype=float)
    if estimated.size:
        _plot_decimated(ax, t, 1.0e3 * estimated, "estimated step wall time", color="tab:green")
        wall_axis = ax.twinx()
        _plot_decimated(wall_axis, t, cumulative / 3600.0, "integrated measured wall time", color="tab:purple")
        ax.set_ylabel("estimated step time [ms]")
        wall_axis.set_ylabel("cumulative wall time [h]", color="tab:purple")
        wall_axis.tick_params(axis="y", colors="tab:purple")
        _combined_legend(ax, wall_axis, fontsize=8)
    else:
        _plot_decimated(ax, t, dataset.get("solve_attempts_per_step"), "attempt-unit proxy")
        _plot_decimated(ax, t, dataset.get("cumulative_solve_attempts"), "cumulative attempts")
        ax.set_ylabel("solve-attempt units")
        ax.legend(frameon=False, fontsize=8)
    ax.set_title("wall-time accounting proxy")

    ax = axes[1, 0]
    _plot_decimated(ax, t, dataset.get("eta_R"), r"$\eta_R$", positive=True)
    _plot_decimated(
        ax,
        t,
        dataset.get("allmaras_update_forcing_max_abs"),
        r"max $|F|/\Delta_{\mathrm{BCS}}(0)$",
        positive=True,
    )
    ax.set_yscale("log")
    ax.set_title("local-update stiffness")
    ax.legend(frameon=False, fontsize=8)

    ax = axes[1, 1]
    residual = np.asarray(dataset.get("poisson_residual_rel", []), dtype=float)
    _plot_decimated(ax, t, residual, "relative Poisson residual", positive=True)
    poisson_tol = float(dataset.get("poisson_tolerance", np.nan))
    if np.isfinite(poisson_tol) and poisson_tol > 0.0:
        ax.axhline(poisson_tol, color="black", linestyle="--", linewidth=1.2, label="tolerance")
    ax.set_yscale("log")
    ax.set_title("Poisson current conservation")
    ax.legend(frameon=False, fontsize=8)

    ax = axes[1, 2]
    ax.plot(snap_t, dataset.get("div_j_normalized_max_snapshot"), "o-", label="bulk max")
    ax.plot(snap_t, dataset.get("div_j_normalized_rms_snapshot"), "o-", label="bulk RMS")
    ax.set_yscale("log")
    ax.set_title(r"$\xi|\nabla\!\cdot\!\mathbf{j}|/j_{\mathrm{avg}}$")
    ax.legend(frameon=False, fontsize=8)

    ax = axes[2, 0]
    _plot_decimated(ax, t, dataset.get("allmaras_phase_convergence_residual_rel"), "CG residual", positive=True)
    cg_tol = float(dataset.get("phase_convergence_tolerance", np.nan))
    if np.isfinite(cg_tol) and cg_tol > 0.0:
        ax.axhline(cg_tol, color="black", linestyle="--", linewidth=1.2, label="CG tolerance")
    iteration_axis = ax.twinx()
    _plot_decimated(
        iteration_axis,
        t,
        dataset.get("allmaras_phase_convergence_iterations"),
        "CG iterations",
        color="tab:orange",
    )
    ax.set_yscale("log")
    ax.set_ylabel("relative residual")
    iteration_axis.set_ylabel("iterations", color="tab:orange")
    iteration_axis.tick_params(axis="y", colors="tab:orange")
    ax.set_title("harmonic continuation convergence")
    _combined_legend(ax, iteration_axis, fontsize=8)

    ax = axes[2, 1]
    _plot_decimated(ax, t, dataset.get("allmaras_phase_direct_node_count"), "direct")
    _plot_decimated(ax, t, dataset.get("allmaras_phase_continued_node_count"), "continued")
    _plot_decimated(ax, t, dataset.get("allmaras_phase_zero_amplitude_node_count"), "exactly zero")
    ax.set_ylabel("nodes")
    ax.set_title("phase-drive domains")
    ax.legend(frameon=False, fontsize=8)

    ax = axes[2, 2]
    phase_rms = np.asarray(dataset.get("allmaras_phase_drive_rms_snapshot", []), dtype=float)
    phase_max = np.asarray(dataset.get("allmaras_phase_drive_max_snapshot", []), dtype=float)
    if phase_rms.size:
        ax.plot(snap_t[: phase_rms.size], phase_rms, "o-", label="phase-drive RMS")
    if phase_max.size:
        ax.plot(snap_t[: phase_max.size], phase_max, "o-", label="phase-drive max")
    ax.set_yscale("log")
    ax.set_ylabel(r"$|F_\phi|/\Delta_{\mathrm{BCS}}(0)$")
    mismatch_axis = ax.twinx()
    relative_l2 = np.asarray(dataset.get("usadel_vs_gl_relative_l2_snapshot", []), dtype=float)
    if relative_l2.size:
        mismatch_axis.plot(
            snap_t[: relative_l2.size],
            relative_l2,
            "s--",
            color="tab:red",
            label=r"$\|j_s^{Us}-j_s^{GL}\|/\|j_s^{Us}\|$",
        )
    mismatch_axis.set_ylabel("Usadel-GL relative mismatch", color="tab:red")
    mismatch_axis.tick_params(axis="y", colors="tab:red")
    ax.set_title("current-law correction")
    _combined_legend(ax, mismatch_axis, fontsize=8)

    for ax in axes.ravel():
        ax.set_xlabel("t [ps]")
        ax.grid(True, alpha=0.22)
        ax.set_xlim(left=0.0)

    converged = np.asarray(dataset.get("allmaras_phase_convergence_converged", []), dtype=bool)
    wall = dataset.get("measured_wall_time_s")
    wall_text = "not supplied"
    if wall is not None:
        wall_text = f"{float(wall) / 3600.0:.3f} h (measured total; per-step allocation is an attempt-count estimate)"
    footer = (
        f"continuity pass={bool(dataset.get('continuity_passes', False))}; "
        f"stationarity pass={bool(dataset.get('stationarity_passes', False))}; "
        f"CG converged on all accepted steps={bool(converged.size and np.all(converged))}; "
        f"wall time={wall_text}"
    )
    fig.text(0.5, 0.018, footer, ha="center", va="bottom", fontsize=9)

    fig.savefig(output, dpi=dpi, bbox_inches="tight", pad_inches=0.10)
    plt.close(fig)
    return output


def _prepare_output(path: str | Path) -> Path:
    output = Path(path)
    output.parent.mkdir(parents=True, exist_ok=True)
    return output


def _global_limits(
    values: np.ndarray,
    *,
    symmetric: bool,
    forced_min: float | None,
    forced_max: float | None,
) -> tuple[float, float]:
    finite = np.asarray(values, dtype=float)
    finite = finite[np.isfinite(finite)]
    if finite.size == 0:
        return -1.0 if symmetric else 0.0, 1.0
    if symmetric:
        limit = float(np.nanpercentile(np.abs(finite), 99.7))
        limit = max(limit, 1.0e-30)
        return -limit, limit
    vmin = float(forced_min) if forced_min is not None else float(np.nanpercentile(finite, 0.3))
    vmax = float(forced_max) if forced_max is not None else float(np.nanpercentile(finite, 99.7))
    if not np.isfinite(vmax) or vmax <= vmin:
        vmax = vmin + 1.0
    return vmin, vmax


def _decimate(x: Any, y: Any, *, max_points: int = 6000) -> tuple[np.ndarray, np.ndarray]:
    x_values = np.asarray(x, dtype=float).reshape(-1)
    y_values = np.asarray(y if y is not None else [], dtype=float).reshape(-1)
    if x_values.size == 0 or y_values.size == 0:
        return np.array([], dtype=float), np.array([], dtype=float)
    if y_values.size != x_values.size:
        y_values = np.resize(y_values, x_values.size)
    stride = max(1, int(np.ceil(x_values.size / max(1, int(max_points)))))
    return x_values[::stride], y_values[::stride]


def _plot_decimated(
    ax,
    x: Any,
    y: Any,
    label: str,
    *,
    positive: bool = False,
    color: str | None = None,
) -> None:
    x_values, y_values = _decimate(x, y)
    mask = np.isfinite(x_values) & np.isfinite(y_values)
    if positive:
        mask &= y_values > 0.0
    if np.any(mask):
        ax.plot(x_values[mask], y_values[mask], label=label, color=color, linewidth=1.0)


def _combined_legend(ax, twin, *, ncol: int = 1, fontsize: int | None = None) -> None:
    handles_a, labels_a = ax.get_legend_handles_labels()
    handles_b, labels_b = twin.get_legend_handles_labels()
    handles = handles_a + handles_b
    labels = labels_a + labels_b
    if handles:
        ax.legend(handles, labels, frameon=False, ncol=ncol, fontsize=fontsize, loc="best")
