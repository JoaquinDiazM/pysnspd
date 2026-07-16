"""Thesis figures for a normalized phase-continuation SS run."""

from __future__ import annotations

from pathlib import Path
from typing import Any, Mapping

import matplotlib

matplotlib.use("Agg", force=True)
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
from matplotlib.colors import LogNorm, Normalize, SymLogNorm
import numpy as np

from pysnspd.plotting.style import THESIS_DPI, THESIS_WIDTH_IN, apply_thesis_style

apply_thesis_style()


def make_phasecg_ss_figures(
    *,
    dataset: Mapping[str, Any],
    output_dir: str | Path,
    dpi: int = THESIS_DPI,
) -> dict[str, Path]:
    """Create the final E2 physical and numerical figures from stored data."""

    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)
    saved = {
        "snapshot_fields": plot_phasecg_snapshot_fields(
            dataset,
            out / "E2_ss_snapshot_fields.pdf",
            dpi=dpi,
        ),
        "physical_evolution": plot_phasecg_physical_evolution(
            dataset,
            out / "E2_ss_physical_evolution.pdf",
            dpi=dpi,
        ),
        "numerical_diagnostics": plot_phasecg_numerical_diagnostics(
            dataset,
            out / "E2_ss_numerical_diagnostics.pdf",
            dpi=dpi,
        ),
    }
    conversion_keys = (
        "jtot_x_snapshot_over_javg",
        "js_x_snapshot_over_javg",
        "jn_x_snapshot_over_javg",
        "node_area_m2",
    )
    if all(np.asarray(dataset.get(key, [])).size for key in conversion_keys):
        saved["current_conversion_profiles"] = plot_current_conversion_profiles(
            dataset,
            out / "E2_ss_current_conversion_profiles.pdf",
            dpi=dpi,
        )
    if np.asarray(dataset.get("joule_snapshot_W_m3", [])).size:
        saved["thermal_balance"] = plot_ss_thermal_balance(
            dataset,
            out / "E2_ss_thermal_balance.pdf",
            dpi=dpi,
        )
    return saved


def plot_phasecg_snapshot_fields(
    dataset: Mapping[str, Any],
    output_path: str | Path,
    *,
    dpi: int = THESIS_DPI,
) -> Path:
    """Plot six fundamental fields at every stored SS snapshot."""

    output = _prepare_output(output_path)
    times = np.asarray(dataset.get("snapshot_t_ps", []), dtype=float)
    if times.size == 0:
        raise ValueError("No snapshot times are available for the E2 field figure.")

    tri = _triangulation(dataset)
    delta_field = np.asarray(dataset["delta_snapshot_over_delta0"], dtype=float)
    Te_field = np.asarray(dataset.get("Te_snapshot_K", np.zeros_like(delta_field)), dtype=float)
    Tc_K = _critical_temperature_K(dataset)
    finite_Te = Te_field[np.isfinite(Te_field)]
    Te_vmax = float(np.nanpercentile(finite_Te, 99.7)) if finite_Te.size else 1.0
    if np.isfinite(Tc_K):
        Te_vmax = max(Te_vmax, Tc_K)
    fields = [
        (
            delta_field,
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
            r"$|\mathbf{q}|\xi$",
            "magma",
            False,
            0.0,
            None,
        ),
        (
            Te_field,
            r"$T_e$ [K]",
            "inferno",
            False,
            None,
            Te_vmax,
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
    ]

    n_rows = times.size
    n_cols = len(fields)
    fig, axes = plt.subplots(
        n_rows,
        n_cols,
        figsize=(THESIS_WIDTH_IN, max(4.7, 0.49 * n_rows + 1.0)),
        squeeze=False,
    )
    fig.subplots_adjust(left=0.064, right=0.948, bottom=0.055, top=0.910, wspace=0.08, hspace=0.10)

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
            _format_strip_axis(ax, tri)
            if row < n_rows - 1:
                ax.tick_params(axis="x", labelbottom=False)
            else:
                ax.set_xlabel(r"$x$ [nm]", labelpad=1.0)
            if col == 0:
                ax.set_ylabel(r"$y$ [nm]", labelpad=1.0)
            else:
                ax.tick_params(axis="y", labelleft=False)
            if col == n_cols - 1:
                ax.text(
                    1.05,
                    0.5,
                    rf"$t={times[row]:.3g}$ ps",
                    transform=ax.transAxes,
                    rotation=-90,
                    va="center",
                    ha="left",
                    fontsize=7.2,
                )
        column_mappables.append((mappable, label))

    fig.canvas.draw()
    for col, (mappable, label) in enumerate(column_mappables):
        position = axes[0, col].get_position()
        cax = fig.add_axes([position.x0, position.y1 + 0.009, position.width, 0.010])
        colorbar = fig.colorbar(mappable, cax=cax, orientation="horizontal")
        cax.xaxis.set_ticks_position("top")
        cax.xaxis.set_label_position("top")
        colorbar.set_label(label, labelpad=1.5, fontsize=8.2)
        if label == r"$T_e$ [K]" and np.isfinite(Tc_K):
            lower = float(mappable.norm.vmin)
            middle = 0.5 * (lower + Tc_K)
            colorbar.set_ticks([lower, middle, Tc_K])
            colorbar.set_ticklabels([f"{lower:.2g}", f"{middle:.2f}", r"$T_c$"])
        cax.tick_params(labelsize=6.8, pad=0.8, length=2.0)
        tick_labels = colorbar.ax.get_xticklabels()
        if len(tick_labels) >= 2:
            tick_labels[0].set_ha("left")
            tick_labels[-1].set_ha("right")

    fig.savefig(output, dpi=dpi)
    plt.close(fig)
    return output


def plot_phasecg_physical_evolution(
    dataset: Mapping[str, Any],
    output_path: str | Path,
    *,
    dpi: int = THESIS_DPI,
) -> Path:
    """Plot current partition, voltage, temperatures and condensate response."""

    output = _prepare_output(output_path)
    snapshot_t = np.asarray(dataset.get("snapshot_t_ps", []), dtype=float)
    history_t = np.asarray(dataset.get("t_ps", []), dtype=float)
    fig, axes = plt.subplots(4, 1, figsize=(THESIS_WIDTH_IN, 6.25), sharex=True)
    fig.subplots_adjust(left=0.115, right=0.870, bottom=0.080, top=0.985, hspace=0.18)

    target_current = float(dataset.get("target_current_uA", np.nan))
    axes[0].axhline(
        target_current,
        color="0.25",
        linestyle="--",
        linewidth=0.9,
        label="Target current",
    )
    axes[0].plot(snapshot_t, dataset.get("current_total_snapshot_uA"), "o-", label=r"$I_{\mathrm{tot}}$")
    axes[0].plot(snapshot_t, dataset.get("current_super_snapshot_uA"), "o-", label=r"$I_s^{\mathrm{Us}}$")
    axes[0].plot(snapshot_t, dataset.get("current_normal_snapshot_uA"), "o-", label=r"$I_n$")
    axes[0].set_ylabel(r"Current [$\mu$A]")
    axes[0].legend(frameon=False, ncol=4, loc="best")

    terminal_t, terminal_v = _decimate(history_t, dataset.get("terminal_voltage_mV"), max_points=7000)
    axes[1].plot(terminal_t, terminal_v, linewidth=0.9, label="Terminal voltage")
    axes[1].plot(
        snapshot_t,
        dataset.get("voltage_center_snapshot_mV"),
        "o-",
        label="Central 100 nm voltage",
    )
    axes[1].plot(
        snapshot_t,
        dataset.get("voltage_terminal_snapshot_mV"),
        "o",
        label="Terminal snapshots",
    )
    axes[1].set_ylabel("Voltage [mV]")
    axes[1].legend(frameon=False, ncol=3, loc="best")

    Te = np.asarray(dataset.get("Te_snapshot_K", []), dtype=float)
    Tph = np.asarray(dataset.get("Tph_snapshot_K", []), dtype=float)
    if Te.ndim == 2 and Te.shape[0] == snapshot_t.size:
        axes[2].plot(snapshot_t, np.nanmax(Te, axis=1), "o-", label=r"Max $T_e$")
        axes[2].plot(snapshot_t, np.nanmean(Te, axis=1), "o-", label=r"Mean $T_e$")
    if Tph.ndim == 2 and Tph.shape[0] == snapshot_t.size:
        axes[2].plot(snapshot_t, np.nanmax(Tph, axis=1), "s--", label=r"Max $T_{ph}$")
    axes[2].set_ylabel("Temperature [K]")
    thermal_handles, thermal_labels = axes[2].get_legend_handles_labels()
    if thermal_handles:
        axes[2].legend(thermal_handles, thermal_labels, frameon=False, ncol=3, loc="best")

    axes[3].plot(snapshot_t, dataset.get("delta_center_min"), "o-", label="Central minimum")
    axes[3].plot(snapshot_t, dataset.get("delta_center_mean"), "o-", label="Central mean")
    axes[3].plot(snapshot_t, dataset.get("delta_center_max"), "o-", label="Central maximum")
    axes[3].set_ylabel(r"$|\Delta|/\Delta_{\mathrm{BCS}}(0)$")
    axes[3].set_xlabel(r"$t$ [ps]")
    fraction_axis = axes[3].twinx()
    fraction_axis.plot(
        snapshot_t,
        dataset.get("normal_current_fraction_snapshot"),
        "s--",
        color="tab:red",
        label=r"$|I_n/I_{\mathrm{tot}}|$",
    )
    fraction_axis.set_ylabel("Normal-current fraction", color="tab:red")
    fraction_axis.tick_params(axis="y", colors="tab:red")
    _combined_legend(axes[3], fraction_axis, ncol=2)

    for ax in axes:
        ax.grid(True)
        ax.set_xlim(left=0.0)

    fig.savefig(output, dpi=dpi)
    plt.close(fig)
    return output


def plot_current_conversion_profiles(
    dataset: Mapping[str, Any],
    output_path: str | Path,
    *,
    dpi: int = THESIS_DPI,
) -> Path:
    """Plot the final longitudinal current conversion and condensate profiles."""

    output = _prepare_output(output_path)
    x = np.asarray(dataset.get("nodes_x_nm", []), dtype=float)
    weights = np.asarray(dataset.get("node_area_m2", np.ones_like(x)), dtype=float)
    if x.size == 0:
        raise ValueError("The E2 conversion figure requires mesh coordinates.")

    profiles = {}
    for key in (
        "jtot_x_snapshot_over_javg",
        "js_x_snapshot_over_javg",
        "jn_x_snapshot_over_javg",
        "delta_snapshot_over_delta0",
        "phi_snapshot_mV",
    ):
        values = np.asarray(dataset.get(key, []), dtype=float)
        if values.ndim != 2 or values.shape[1] != x.size:
            raise ValueError(f"Missing node-resolved snapshot field: {key}")
        profiles[key] = _binned_profile(x, values[-1], weights, n_bins=90)

    x_profile = profiles["jtot_x_snapshot_over_javg"][0]
    jtot = profiles["jtot_x_snapshot_over_javg"][1]
    js = profiles["js_x_snapshot_over_javg"][1]
    jn = profiles["jn_x_snapshot_over_javg"][1]
    delta = profiles["delta_snapshot_over_delta0"][1]
    phi = profiles["phi_snapshot_mV"][1]

    fig, axes = plt.subplots(2, 1, figsize=(THESIS_WIDTH_IN, 4.05), sharex=True)
    fig.subplots_adjust(left=0.115, right=0.870, bottom=0.120, top=0.975, hspace=0.16)

    axes[0].plot(x_profile, jtot, label=r"$j_{\mathrm{tot},x}/j_{\mathrm{avg}}$")
    axes[0].plot(x_profile, js, label=r"$j_{s,x}^{\mathrm{Us}}/j_{\mathrm{avg}}$")
    axes[0].plot(x_profile, jn, label=r"$j_{n,x}/j_{\mathrm{avg}}$")
    for side, color in (("left", "0.25"), ("right", "0.45")):
        fitted = _fit_conversion_exponential(x_profile, jn, side=side)
        if fitted is not None:
            fit_values, length_nm = fitted
            axes[0].plot(
                x_profile,
                fit_values,
                linestyle="--",
                color=color,
                linewidth=0.9,
                label=rf"{side.capitalize()} exponential guide, $\ell_Q={length_nm:.1f}$ nm",
            )
    axes[0].set_ylabel(r"Current density / $j_{\mathrm{avg}}$")
    axes[0].legend(frameon=False, ncol=2, loc="best")

    axes[1].plot(x_profile, delta, color="tab:blue", label=r"$|\Delta|/\Delta_{\mathrm{BCS}}(0)$")
    axes[1].set_ylabel(r"$|\Delta|/\Delta_{\mathrm{BCS}}(0)$", color="tab:blue")
    axes[1].tick_params(axis="y", colors="tab:blue")
    phi_axis = axes[1].twinx()
    phi_axis.plot(x_profile, phi, color="tab:red", label=r"$\phi$")
    phi_axis.set_ylabel(r"$\phi$ [mV]", color="tab:red")
    phi_axis.tick_params(axis="y", colors="tab:red")
    axes[1].set_xlabel(r"$x$ [nm]")
    _combined_legend(axes[1], phi_axis, ncol=2)

    for ax in axes:
        ax.grid(True)
        ax.set_xlim(float(np.nanmin(x_profile)), float(np.nanmax(x_profile)))

    fig.savefig(output, dpi=dpi)
    plt.close(fig)
    return output


def plot_ss_thermal_balance(
    dataset: Mapping[str, Any],
    output_path: str | Path,
    *,
    dpi: int = THESIS_DPI,
) -> Path:
    """Plot the stored Joule, electron-phonon and escape power balance."""

    output = _prepare_output(output_path)
    tri = _triangulation(dataset)
    times = np.asarray(dataset.get("snapshot_t_ps", []), dtype=float)
    joule = np.asarray(dataset.get("joule_snapshot_W_m3", []), dtype=float)
    electron_phonon = np.asarray(dataset.get("P_total_snapshot_W_m3", []), dtype=float)
    escape = np.asarray(dataset.get("P_esc_snapshot_W_m3", []), dtype=float)
    channels = [
        (joule, r"$P_J$", "magma"),
        (electron_phonon, r"$P_{e\text{-}ph}=P_S+P_R$", "coolwarm"),
        (escape, r"$P_{\mathrm{esc}}$", "coolwarm"),
    ]

    fig = plt.figure(figsize=(THESIS_WIDTH_IN, 4.5))
    grid = fig.add_gridspec(2, 3, height_ratios=(1.0, 0.88), left=0.075, right=0.985, bottom=0.105, top=0.965, wspace=0.42, hspace=0.38)
    final_time = float(times[-1]) if times.size else np.nan

    for col, (values, label, cmap) in enumerate(channels):
        ax = fig.add_subplot(grid[0, col])
        final = values[-1]
        norm = _power_norm(final, positive=(col == 0))
        mappable = ax.tripcolor(
            tri,
            final,
            shading="gouraud",
            cmap=cmap,
            norm=norm,
            rasterized=True,
        )
        _format_strip_axis(ax, tri)
        ax.set_xlabel(r"$x$ [nm]", labelpad=1.0)
        if col == 0:
            ax.set_ylabel(r"$y$ [nm]", labelpad=1.0)
        else:
            ax.tick_params(axis="y", labelleft=False)
        ax.set_title(label)
        cbar = fig.colorbar(mappable, ax=ax, orientation="horizontal", fraction=0.075, pad=0.18)
        cbar.set_label(r"W m$^{-3}$", fontsize=8.0, labelpad=1.0)
        cbar.ax.tick_params(labelsize=6.8, pad=1.0)
        max_abs = float(np.nanmax(np.abs(final))) if np.isfinite(final).any() else np.nan
        ax.text(
            0.98,
            0.94,
            rf"$\max|P|={max_abs:.2e}$",
            transform=ax.transAxes,
            ha="right",
            va="top",
            fontsize=7.0,
            color="white" if max_abs != 0.0 else "0.15",
        )

    joule_axis = fig.add_subplot(grid[1, :2])
    joule_max = np.nanmax(np.abs(joule), axis=1)
    joule_axis.plot(times, joule_max, "o-", label=r"Max $|P_J|$")
    joule_axis.set_yscale("log")
    joule_axis.set_xlabel(r"$t$ [ps]")
    joule_axis.set_ylabel(r"$\max |P_J|$ [W m$^{-3}$]")
    joule_axis.grid(True)
    joule_axis.legend(frameon=False, loc="best")
    joule_axis.text(0.99, 0.06, rf"Final snapshot: $t={final_time:.3g}$ ps", transform=joule_axis.transAxes, ha="right", fontsize=8.0)

    exchange_axis = fig.add_subplot(grid[1, 2])
    electron_phonon_max = np.nanmax(np.abs(electron_phonon), axis=1)
    escape_max = np.nanmax(np.abs(escape), axis=1)
    exchange_axis.plot(times, electron_phonon_max, "o-", label=r"Max $|P_{e\text{-}ph}|$")
    if np.any(escape_max > 0.0):
        exchange_axis.plot(times, escape_max, "s--", label=r"Max $|P_{\mathrm{esc}}|$")
    else:
        exchange_axis.text(0.98, 0.08, r"$P_{\mathrm{esc}}=0$", transform=exchange_axis.transAxes, ha="right", fontsize=8.0)
    positive_exchange = electron_phonon_max[electron_phonon_max > 0.0]
    if positive_exchange.size:
        exchange_axis.set_yscale("log")
    exchange_axis.set_xlabel(r"$t$ [ps]")
    exchange_axis.set_ylabel(r"Interaction power [W m$^{-3}$]")
    exchange_axis.grid(True)
    exchange_axis.legend(frameon=False, loc="best")

    fig.savefig(output, dpi=dpi)
    plt.close(fig)
    return output


def plot_phasecg_numerical_diagnostics(
    dataset: Mapping[str, Any],
    output_path: str | Path,
    *,
    dpi: int = THESIS_DPI,
) -> Path:
    """Plot adaptive, thermal, continuity and continuation diagnostics."""

    output = _prepare_output(output_path)
    t = np.asarray(dataset.get("t_ps", []), dtype=float)
    snap_t = np.asarray(dataset.get("snapshot_t_ps", []), dtype=float)
    fig, axes = plt.subplots(3, 3, figsize=(THESIS_WIDTH_IN, 6.25))
    fig.subplots_adjust(left=0.085, right=0.930, bottom=0.105, top=0.975, wspace=0.48, hspace=0.42)

    ax = axes[0, 0]
    _plot_decimated(ax, t, dataset.get("dt_attempt_fs"), "Attempted", positive=True)
    _plot_decimated(ax, t, dataset.get("dt_accepted_fs"), "Accepted", positive=True)
    _plot_decimated(ax, t, dataset.get("dt_next_fs"), "Next", positive=True)
    ax.set_yscale("log")
    ax.set_ylabel(r"$\Delta t$ [fs]")
    ax.set_title("Adaptive time step")
    ax.legend(frameon=False, fontsize=7.0)

    ax = axes[0, 1]
    _plot_decimated(ax, t, dataset.get("solve_attempts_per_step"), "Solve attempts")
    _plot_decimated(ax, t, dataset.get("adaptive_retries"), "Retries")
    rejected_axis = ax.twinx()
    _plot_decimated(rejected_axis, t, dataset.get("cumulative_rejected_attempts"), "Cumulative rejected", color="tab:red")
    ax.set_ylabel("Count / accepted step")
    rejected_axis.tick_params(axis="y", labelright=False, right=False)
    ax.set_title("Nonlinear solve effort")
    _combined_legend(ax, rejected_axis, fontsize=7.0)

    ax = axes[0, 2]
    _plot_decimated(ax, t, dataset.get("thermal_max_Te_K_history"), r"Max $T_e$")
    _plot_decimated(ax, t, dataset.get("thermal_max_Tph_K_history"), r"Max $T_{ph}$")
    rate_axis = ax.twinx()
    _plot_decimated(rate_axis, t, dataset.get("thermal_max_rate_K_per_ps_history"), "Max thermal rate", positive=True, color="tab:red")
    ax.set_ylabel("Temperature [K]")
    rate_axis.tick_params(axis="y", colors="tab:red")
    rate_axis.set_yscale("log")
    ax.set_title("Thermal evolution")
    _combined_legend(ax, rate_axis, fontsize=7.0)

    ax = axes[1, 0]
    _plot_decimated(ax, t, dataset.get("eta_R"), r"$\eta_R$", positive=True)
    _plot_decimated(ax, t, dataset.get("allmaras_update_forcing_max_abs"), r"Max $|F|/\Delta_{\mathrm{BCS}}(0)$", positive=True)
    ax.set_yscale("log")
    ax.set_title("Local-update stiffness")
    ax.legend(frameon=False, fontsize=7.0)

    ax = axes[1, 1]
    _plot_decimated(ax, t, dataset.get("poisson_residual_rel"), "Relative Poisson residual", positive=True)
    poisson_tol = float(dataset.get("poisson_tolerance", np.nan))
    if np.isfinite(poisson_tol) and poisson_tol > 0.0:
        ax.axhline(poisson_tol, color="0.25", linestyle="--", linewidth=0.9, label="Tolerance")
    ax.set_yscale("log")
    ax.set_title("Poisson current conservation")
    ax.legend(frameon=False, fontsize=7.0)

    ax = axes[1, 2]
    ax.plot(snap_t, dataset.get("div_j_normalized_max_snapshot"), "o-", label="Bulk maximum")
    ax.plot(snap_t, dataset.get("div_j_normalized_rms_snapshot"), "o-", label="Bulk RMS")
    ax.set_yscale("log")
    ax.set_title(r"$\xi|\nabla\!\cdot\!\mathbf{j}|/j_{\mathrm{avg}}$")
    ax.legend(frameon=False, fontsize=7.0)

    ax = axes[2, 0]
    _plot_decimated(ax, t, dataset.get("allmaras_phase_convergence_residual_rel"), "CG residual", positive=True)
    cg_tol = float(dataset.get("phase_convergence_tolerance", np.nan))
    if np.isfinite(cg_tol) and cg_tol > 0.0:
        ax.axhline(cg_tol, color="0.25", linestyle="--", linewidth=0.9, label="CG tolerance")
    iteration_axis = ax.twinx()
    _plot_decimated(iteration_axis, t, dataset.get("allmaras_phase_convergence_iterations"), "CG iterations", color="tab:red")
    ax.set_yscale("log")
    ax.set_ylabel("Relative residual")
    iteration_axis.tick_params(axis="y", labelright=False, right=False)
    ax.set_title("Harmonic continuation")
    _combined_legend(ax, iteration_axis, fontsize=7.0)

    ax = axes[2, 1]
    _plot_decimated(ax, t, dataset.get("allmaras_phase_direct_node_count"), "Direct nodes")
    _plot_decimated(ax, t, dataset.get("allmaras_phase_continued_node_count"), "Continued nodes")
    _plot_decimated(ax, t, dataset.get("allmaras_phase_zero_amplitude_node_count"), "Zero-amplitude nodes")
    ax.set_ylabel("Nodes")
    ax.set_title("Continuation domains")
    ax.legend(frameon=False, fontsize=7.0)

    ax = axes[2, 2]
    phase_rms = np.asarray(dataset.get("allmaras_phase_drive_rms_snapshot", []), dtype=float)
    phase_max = np.asarray(dataset.get("allmaras_phase_drive_max_snapshot", []), dtype=float)
    if phase_rms.size:
        ax.plot(snap_t[: phase_rms.size], phase_rms, "o-", label="Phase-drive RMS")
    if phase_max.size:
        ax.plot(snap_t[: phase_max.size], phase_max, "o-", label="Phase-drive maximum")
    ax.set_yscale("log")
    mismatch_axis = ax.twinx()
    mismatch = np.asarray(dataset.get("usadel_vs_gl_relative_l2_snapshot", []), dtype=float)
    if mismatch.size:
        mismatch_axis.plot(snap_t[: mismatch.size], mismatch, "s--", color="tab:red", label="Usadel-GL mismatch")
    mismatch_axis.set_ylabel("Relative mismatch", color="tab:red")
    mismatch_axis.tick_params(axis="y", colors="tab:red")
    ax.set_title("Current-law correction")
    _combined_legend(ax, mismatch_axis, fontsize=7.0)

    for ax in axes.ravel():
        ax.set_xlabel(r"$t$ [ps]")
        ax.grid(True)
        ax.set_xlim(left=0.0)

    status = (
        f"Continuity: {bool(dataset.get('continuity_passes', False))}; "
        f"dynamic stationarity: {bool(dataset.get('dynamic_stationarity_passes', False))}; "
        f"thermal stationarity: {bool(dataset.get('thermal_stationarity_passes', False))}"
    )
    fig.text(0.5, 0.025, status, ha="center", fontsize=7.5)
    fig.savefig(output, dpi=dpi)
    plt.close(fig)
    return output


def _prepare_output(path: str | Path) -> Path:
    output = Path(path)
    if output.suffix.lower() != ".pdf":
        output = output.with_suffix(".pdf")
    output.parent.mkdir(parents=True, exist_ok=True)
    return output


def _critical_temperature_K(dataset: Mapping[str, Any]) -> float:
    try:
        Tc_K = float(dataset.get("Tc_K", np.nan))
    except Exception:
        Tc_K = float("nan")
    if np.isfinite(Tc_K) and Tc_K > 0.0:
        return Tc_K
    try:
        delta0_meV = float(np.asarray(dataset.get("delta0_meV", np.nan), dtype=float).reshape(-1)[-1])
    except Exception:
        return float("nan")
    if not np.isfinite(delta0_meV) or delta0_meV <= 0.0:
        return float("nan")
    return delta0_meV / (1.764 * 8.617333262e-2)


def _triangulation(dataset: Mapping[str, Any]) -> mtri.Triangulation:
    return mtri.Triangulation(
        np.asarray(dataset["nodes_x_nm"], dtype=float),
        np.asarray(dataset["nodes_y_nm"], dtype=float),
        np.asarray(dataset["triangles"], dtype=np.int64),
    )


def _format_strip_axis(ax: plt.Axes, tri: mtri.Triangulation) -> None:
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlim(float(np.nanmin(tri.x)), float(np.nanmax(tri.x)))
    ax.set_ylim(float(np.nanmin(tri.y)), float(np.nanmax(tri.y)))
    ax.grid(False)
    ax.tick_params(labelsize=6.8, length=2.0, pad=1.0)


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
        limit = max(float(np.nanpercentile(np.abs(finite), 99.7)), 1.0e-30)
        return -limit, limit
    vmin = float(forced_min) if forced_min is not None else float(np.nanpercentile(finite, 0.3))
    vmax = float(forced_max) if forced_max is not None else float(np.nanpercentile(finite, 99.7))
    if not np.isfinite(vmax) or vmax <= vmin:
        vmax = vmin + 1.0
    return vmin, vmax


def _binned_profile(
    x_nm: np.ndarray,
    values: np.ndarray,
    weights: np.ndarray,
    *,
    n_bins: int,
) -> tuple[np.ndarray, np.ndarray]:
    x = np.asarray(x_nm, dtype=float)
    y = np.asarray(values, dtype=float)
    w = np.asarray(weights, dtype=float)
    edges = np.linspace(float(np.nanmin(x)), float(np.nanmax(x)), int(n_bins) + 1)
    index = np.clip(np.digitize(x, edges) - 1, 0, int(n_bins) - 1)
    centers = 0.5 * (edges[:-1] + edges[1:])
    averaged = np.full(int(n_bins), np.nan, dtype=float)
    for i in range(int(n_bins)):
        mask = (index == i) & np.isfinite(y) & np.isfinite(w) & (w > 0.0)
        if np.any(mask):
            averaged[i] = float(np.average(y[mask], weights=w[mask]))
    valid = np.isfinite(averaged)
    return centers[valid], averaged[valid]


def _fit_conversion_exponential(
    x_nm: np.ndarray,
    normal_current: np.ndarray,
    *,
    side: str,
) -> tuple[np.ndarray, float] | None:
    x = np.asarray(x_nm, dtype=float)
    current = np.asarray(normal_current, dtype=float)
    length = float(np.nanmax(x) - np.nanmin(x))
    if x.size < 8 or not np.isfinite(length) or length <= 0.0:
        return None
    center = (x >= np.nanmin(x) + 0.40 * length) & (x <= np.nanmin(x) + 0.60 * length)
    baseline = float(np.nanmedian(current[center])) if np.any(center) else float(np.nanmedian(current))
    if side == "left":
        distance = x - float(np.nanmin(x))
    elif side == "right":
        distance = float(np.nanmax(x)) - x
    else:
        raise ValueError(f"Unknown side: {side}")
    excess_signed = current - baseline
    contact = distance <= 0.35 * length
    amplitude = float(np.nanmax(np.abs(excess_signed[contact]))) if np.any(contact) else 0.0
    mask = contact & np.isfinite(excess_signed) & (np.abs(excess_signed) > max(0.05 * amplitude, 1.0e-10))
    if np.count_nonzero(mask) < 5:
        return None
    slope, intercept = np.polyfit(distance[mask], np.log(np.abs(excess_signed[mask])), 1)
    if not np.isfinite(slope) or slope >= 0.0:
        return None
    length_nm = -1.0 / float(slope)
    if not np.isfinite(length_nm) or length_nm <= 0.0 or length_nm > length:
        return None
    sign = float(np.sign(np.nanmedian(excess_signed[mask & (distance <= np.nanpercentile(distance[mask], 40.0))])))
    if sign == 0.0:
        sign = 1.0
    fit = baseline + sign * np.exp(intercept + slope * distance)
    fit[distance > 0.35 * length] = np.nan
    return fit, length_nm


def _power_norm(values: np.ndarray, *, positive: bool):
    finite = np.asarray(values, dtype=float)
    finite = finite[np.isfinite(finite)]
    if finite.size == 0 or float(np.nanmax(np.abs(finite))) == 0.0:
        return Normalize(vmin=-1.0, vmax=1.0)
    if positive and np.all(finite > 0.0):
        vmin = max(float(np.nanpercentile(finite, 1.0)), float(np.nanmax(finite)) * 1.0e-8)
        vmax = max(float(np.nanpercentile(finite, 99.7)), vmin * 10.0)
        return LogNorm(vmin=vmin, vmax=vmax)
    vmax = max(float(np.nanpercentile(np.abs(finite), 99.7)), 1.0)
    return SymLogNorm(linthresh=max(vmax * 1.0e-6, 1.0), vmin=-vmax, vmax=vmax, base=10.0)


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
    ax: plt.Axes,
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
        ax.plot(x_values[mask], y_values[mask], label=label, color=color, linewidth=0.8)


def _combined_legend(
    ax: plt.Axes,
    twin: plt.Axes,
    *,
    ncol: int = 1,
    fontsize: float | None = None,
) -> None:
    handles_a, labels_a = ax.get_legend_handles_labels()
    handles_b, labels_b = twin.get_legend_handles_labels()
    handles = handles_a + handles_b
    labels = labels_a + labels_b
    if handles:
        ax.legend(handles, labels, frameon=False, ncol=ncol, fontsize=fontsize, loc="best")


__all__ = [
    "make_phasecg_ss_figures",
    "plot_current_conversion_profiles",
    "plot_phasecg_numerical_diagnostics",
    "plot_phasecg_physical_evolution",
    "plot_phasecg_snapshot_fields",
    "plot_ss_thermal_balance",
]
