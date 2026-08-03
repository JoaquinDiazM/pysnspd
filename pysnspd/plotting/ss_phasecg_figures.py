"""Thesis figures for a normalized phase-continuation SS run."""

from __future__ import annotations

from pathlib import Path
from typing import Any, Mapping, Sequence

import matplotlib

matplotlib.use("Agg", force=True)
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
from matplotlib.colors import ListedColormap, LogNorm, Normalize, SymLogNorm
from matplotlib.patches import FancyArrowPatch
import matplotlib.patheffects as path_effects
import numpy as np

from pysnspd.plotting.style import THESIS_DPI, THESIS_WIDTH_IN, apply_thesis_style

apply_thesis_style()

DEFAULT_SNAPSHOT_TIMES_PS = (0.0, 0.2, 0.5, 2.0, 20.0, 100.0, 200.0)
_OBSOLETE_E2_OUTPUTS = (
    "E2_ss_current_conversion_profiles.pdf",
    "E2_ss_numerical_diagnostics.pdf",
    "E2_ss_thermal_balance.pdf",
)


def make_phasecg_ss_figures(
    *,
    dataset: Mapping[str, Any],
    output_dir: str | Path,
    dpi: int = THESIS_DPI,
    requested_snapshot_times_ps: Sequence[float] | None = None,
) -> dict[str, Path]:
    """Create the final E2 physical and numerical figures from stored data."""

    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)
    snapshot_times = (
        DEFAULT_SNAPSHOT_TIMES_PS
        if requested_snapshot_times_ps is None
        else requested_snapshot_times_ps
    )
    saved = {
        "snapshot_fields": plot_phasecg_snapshot_fields(
            dataset,
            out / "E2_ss_snapshot_fields.pdf",
            dpi=dpi,
            requested_times_ps=snapshot_times,
        ),
        "physical_evolution": plot_phasecg_physical_evolution(
            dataset,
            out / "E2_ss_physical_evolution.pdf",
            dpi=dpi,
        ),
        "numerical_procedure": plot_phasecg_numerical_procedure(
            dataset,
            out / "E2_ss_numerical_procedure.pdf",
            dpi=dpi,
        ),
        "stationarity_evolution": plot_phasecg_stationarity_evolution(
            dataset,
            out / "E2_ss_stationarity_evolution.pdf",
            dpi=dpi,
        ),
    }
    conversion_keys = (
        "jtot_x_snapshot_over_javg",
        "js_x_snapshot_over_javg",
        "jn_x_snapshot_over_javg",
        "node_area_m2",
        "Te_snapshot_K",
        "Tph_snapshot_K",
    )
    if all(np.asarray(dataset.get(key, [])).size for key in conversion_keys):
        saved["final_longitudinal_profiles"] = plot_final_longitudinal_profiles(
            dataset,
            out / "E2_ss_final_longitudinal_profiles.pdf",
            dpi=dpi,
        )
    if np.asarray(dataset.get("joule_snapshot_W_m3", [])).size:
        saved["power_density_snapshots"] = plot_ss_power_density_snapshots(
            dataset,
            out / "E2_ss_power_density_snapshots.pdf",
            dpi=dpi,
            requested_times_ps=snapshot_times,
        )
    energy_keys = (
        "u_e_snapshot_J_m3",
        "u_ph_snapshot_J_m3",
        "C_e_snapshot_J_m3_K",
        "C_ph_snapshot_J_m3_K",
    )
    if all(np.asarray(dataset.get(key, [])).size for key in energy_keys):
        saved["energy_heat_capacity_snapshots"] = plot_ss_energy_heat_capacity_snapshots(
            dataset,
            out / "E2_ss_energy_heat_capacity_snapshots.pdf",
            dpi=dpi,
            requested_times_ps=snapshot_times,
        )
    for filename in _OBSOLETE_E2_OUTPUTS:
        (out / filename).unlink(missing_ok=True)
    return saved


def plot_phasecg_snapshot_fields(
    dataset: Mapping[str, Any],
    output_path: str | Path,
    *,
    dpi: int = THESIS_DPI,
    requested_times_ps: Sequence[float] | None = None,
) -> Path:
    """Plot six fundamental fields at selected stored SS snapshots."""

    output = _prepare_output(output_path)
    stored_times = np.asarray(dataset.get("snapshot_t_ps", []), dtype=float)
    if stored_times.size == 0:
        raise ValueError("No snapshot times are available for the E2 field figure.")
    indices = _nearest_unique_snapshot_indices(stored_times, requested_times_ps)
    times = stored_times[indices]

    tri = _triangulation(dataset)
    delta_field = np.asarray(dataset["delta_snapshot_over_delta0"], dtype=float)[indices]
    Te_values = np.asarray(
        dataset.get("Te_snapshot_K", np.zeros_like(np.asarray(dataset["delta_snapshot_over_delta0"]))),
        dtype=float,
    )
    Te_field = Te_values[indices]
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
            np.asarray(dataset["phi_snapshot_mV"], dtype=float)[indices],
            r"$\phi$ [mV]",
            "coolwarm",
            True,
            None,
            None,
        ),
        (
            np.asarray(dataset["qxi_snapshot"], dtype=float)[indices],
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
            np.asarray(dataset["js_snapshot_over_javg"], dtype=float)[indices],
            r"$|\mathbf{j}_s^{\mathrm{Us}}|/j_{\mathrm{avg}}$",
            "viridis",
            False,
            0.0,
            None,
        ),
        (
            np.asarray(dataset["jn_snapshot_over_javg"], dtype=float)[indices],
            r"$|\mathbf{j}_n|/j_{\mathrm{avg}}$",
            "plasma",
            False,
            0.0,
            None,
        ),
    ]
    node_weights = np.asarray(
        dataset.get("node_area_m2", np.ones(np.asarray(dataset["nodes_x_nm"]).size)),
        dtype=float,
    )
    current_vectors: dict[int, tuple[np.ndarray, np.ndarray, float]] = {}
    for column, x_key, y_key in (
        (4, "js_x_snapshot_over_javg", "js_y_snapshot_over_javg"),
        (5, "jn_x_snapshot_over_javg", "jn_y_snapshot_over_javg"),
    ):
        x_values = np.asarray(dataset.get(x_key, []), dtype=float)
        y_values = np.asarray(dataset.get(y_key, []), dtype=float)
        if x_values.ndim == 2 and y_values.shape == x_values.shape:
            mean_x = _weighted_rows_for_plot(x_values[indices], node_weights)
            mean_y = _weighted_rows_for_plot(y_values[indices], node_weights)
            max_magnitude = max(float(np.nanmax(np.hypot(mean_x, mean_y))), 1.0e-300)
            current_vectors[column] = (mean_x, mean_y, max_magnitude)

    n_rows = times.size
    n_cols = len(fields)
    fig, axes = plt.subplots(
        n_rows,
        n_cols,
        figsize=(THESIS_WIDTH_IN, max(4.7, 0.49 * n_rows + 1.0)),
        squeeze=False,
    )
    fig.subplots_adjust(left=0.076, right=0.948, bottom=0.072, top=0.910, wspace=0.08, hspace=0.10)

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
            if col in current_vectors:
                mean_x, mean_y, max_magnitude = current_vectors[col]
                _add_current_direction_arrow(
                    ax,
                    tri,
                    mean_x=float(mean_x[row]),
                    mean_y=float(mean_y[row]),
                    reference_magnitude=max_magnitude,
                )
            if row < n_rows - 1:
                ax.tick_params(axis="x", labelbottom=False)
            if col != 0:
                ax.tick_params(axis="y", labelleft=False)
            if col == n_cols - 1:
                ax.text(
                    1.05,
                    0.5,
                    rf"$t={times[row]:.3g}$ [ps]",
                    transform=ax.transAxes,
                    rotation=-90,
                    va="center",
                    ha="left",
                    fontsize=7.2,
                )
        column_mappables.append((mappable, label))

    fig.supxlabel(r"$x$ [nm]", y=0.018, fontsize=8.5)
    fig.supylabel(r"$y$ [nm]", x=0.016, fontsize=8.5)

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


def _nearest_unique_snapshot_indices(
    stored_times_ps: np.ndarray,
    requested_times_ps: Sequence[float] | None,
) -> np.ndarray:
    stored = np.asarray(stored_times_ps, dtype=float).reshape(-1)
    if stored.size == 0:
        raise ValueError("No stored snapshot times are available.")
    if requested_times_ps is None:
        return np.arange(stored.size, dtype=np.int64)

    requested = np.asarray(list(requested_times_ps), dtype=float).reshape(-1)
    requested = requested[np.isfinite(requested)]
    if requested.size == 0:
        raise ValueError("At least one finite snapshot time must be requested.")

    indices: list[int] = []
    seen: set[int] = set()
    for value in requested:
        index = int(np.nanargmin(np.abs(stored - float(value))))
        if index not in seen:
            seen.add(index)
            indices.append(index)
    return np.asarray(indices, dtype=np.int64)


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
    fig.subplots_adjust(left=0.115, right=0.875, bottom=0.080, top=0.985, hspace=0.16)

    target_current = float(dataset.get("target_current_uA", np.nan))
    axes[0].axhline(
        target_current,
        color="0.25",
        linestyle="--",
        linewidth=0.9,
        label="Target current",
    )
    axes[0].plot(snapshot_t, dataset.get("current_total_snapshot_uA"), color="tab:blue", label=r"$I_{\mathrm{tot}}$")
    axes[0].plot(snapshot_t, dataset.get("current_super_snapshot_uA"), color="tab:purple", label=r"$I_s^{\mathrm{Us}}$")
    axes[0].plot(snapshot_t, dataset.get("current_normal_snapshot_uA"), color="tab:orange", label=r"$I_n$")
    axes[0].set_ylabel(r"Current [$\mu$A]")
    axes[0].legend(frameon=False, ncol=4, loc="best")

    terminal_t, terminal_v = _decimate(history_t, dataset.get("terminal_voltage_mV"), max_points=7000)
    axes[1].plot(terminal_t, terminal_v, color="tab:blue", linewidth=0.9, label=r"$V_{\mathrm{terminal}}$")
    axes[1].plot(
        snapshot_t,
        dataset.get("voltage_center_snapshot_mV"),
        color="tab:orange",
        label=r"$V_{\mathrm{TDGL}}$",
    )
    axes[1].set_ylabel("Voltage [mV]")
    axes[1].legend(frameon=False, ncol=2, loc="best")

    Te = np.asarray(dataset.get("Te_snapshot_K", []), dtype=float)
    Tph = np.asarray(dataset.get("Tph_snapshot_K", []), dtype=float)
    phonon_axis = axes[2].twinx()
    if Te.ndim == 2 and Te.shape[0] == snapshot_t.size:
        axes[2].plot(snapshot_t, np.nanmax(Te, axis=1), color="tab:red", label=r"Max $T_e$")
        axes[2].plot(snapshot_t, np.nanmean(Te, axis=1), color="tab:brown", label=r"Mean $T_e$")
    if Tph.ndim == 2 and Tph.shape[0] == snapshot_t.size:
        phonon_axis.plot(snapshot_t, np.nanmax(Tph, axis=1), color="tab:green", label=r"Max $T_{ph}$")
        phonon_axis.plot(snapshot_t, np.nanmean(Tph, axis=1), color="tab:cyan", linestyle="--", label=r"Mean $T_{ph}$")
    axes[2].set_ylabel(r"$T_e$ [K]", color="tab:red")
    axes[2].tick_params(axis="y", colors="tab:red")
    phonon_axis.set_ylabel(r"$T_{ph}$ [K]", color="tab:green")
    phonon_axis.tick_params(axis="y", colors="tab:green")
    _combined_legend(axes[2], phonon_axis, ncol=2, loc="best")

    axes[3].plot(snapshot_t, dataset.get("delta_center_min"), color="tab:blue", label="Central minimum")
    axes[3].plot(snapshot_t, dataset.get("delta_center_mean"), color="tab:purple", label="Central mean")
    axes[3].plot(snapshot_t, dataset.get("delta_center_max"), color="tab:green", label="Central maximum")
    axes[3].set_ylabel(r"$|\Delta|/\Delta_{\mathrm{BCS}}(0)$")
    axes[3].set_xlabel(r"$t$ [ps]")
    fraction_axis = axes[3].twinx()
    fraction_axis.plot(
        snapshot_t,
        dataset.get("normal_current_fraction_snapshot"),
        linestyle="--",
        color="tab:orange",
        label=r"$|I_n/I_{\mathrm{tot}}|$",
    )
    fraction_axis.set_ylabel("Normal-current fraction", color="tab:orange")
    fraction_axis.tick_params(axis="y", colors="tab:orange")
    _combined_legend(axes[3], fraction_axis, ncol=2, loc="best")

    for ax in axes:
        ax.grid(True)
        ax.set_xlim(left=0.0)

    fig.savefig(output, dpi=dpi)
    plt.close(fig)
    return output


def plot_final_longitudinal_profiles(
    dataset: Mapping[str, Any],
    output_path: str | Path,
    *,
    dpi: int = THESIS_DPI,
) -> Path:
    """Plot final longitudinal current, condensate, voltage and temperature profiles."""

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
        "Te_snapshot_K",
        "Tph_snapshot_K",
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
    Te = profiles["Te_snapshot_K"][1]
    Tph = profiles["Tph_snapshot_K"][1]

    fig, axes = plt.subplots(3, 1, figsize=(THESIS_WIDTH_IN, 5.55), sharex=True)
    fig.subplots_adjust(left=0.115, right=0.875, bottom=0.090, top=0.980, hspace=0.15)

    axes[0].plot(x_profile, jtot, color="tab:blue", label=r"$j_{\mathrm{tot},x}/j_{\mathrm{avg}}$")
    axes[0].plot(x_profile, js, color="tab:purple", label=r"$j_{s,x}^{\mathrm{Us}}/j_{\mathrm{avg}}$")
    axes[0].plot(x_profile, jn, color="tab:orange", label=r"$j_{n,x}/j_{\mathrm{avg}}$")
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
                label=rf"{side.capitalize()} exponential fit",
            )
    axes[0].set_ylabel(r"$j / j_{\mathrm{avg}}$")
    axes[0].legend(frameon=False, ncol=2, loc="best")

    axes[1].plot(x_profile, delta, color="tab:green", label=r"$|\Delta|/\Delta_{\mathrm{BCS}}(0)$")
    axes[1].set_ylabel(r"$|\Delta|/\Delta_{\mathrm{BCS}}(0)$", color="tab:green")
    axes[1].tick_params(axis="y", colors="tab:green")
    phi_axis = axes[1].twinx()
    phi_axis.plot(x_profile, phi, color="tab:red", label=r"$\phi$")
    phi_axis.set_ylabel(r"$\phi$ [mV]", color="tab:red")
    phi_axis.tick_params(axis="y", colors="tab:red")
    _combined_legend(axes[1], phi_axis, ncol=2, loc="lower right")

    axes[2].plot(x_profile, Te, color="tab:brown", label=r"$T_e$")
    axes[2].set_ylabel(r"$T_e$ [K]", color="tab:brown")
    axes[2].tick_params(axis="y", colors="tab:brown")
    phonon_axis = axes[2].twinx()
    phonon_axis.plot(x_profile, Tph, color="tab:cyan", label=r"$T_{ph}$")
    phonon_axis.set_ylabel(r"$T_{ph}$ [K]", color="tab:cyan")
    phonon_axis.tick_params(axis="y", colors="tab:cyan")
    axes[2].set_xlabel(r"$x$ [nm]")
    _combined_legend(axes[2], phonon_axis, ncol=2, loc="lower right")

    for ax in axes:
        ax.grid(True)
        ax.set_xlim(float(np.nanmin(x_profile)), float(np.nanmax(x_profile)))

    fig.savefig(output, dpi=dpi)
    plt.close(fig)
    return output


def plot_ss_power_density_snapshots(
    dataset: Mapping[str, Any],
    output_path: str | Path,
    *,
    dpi: int = THESIS_DPI,
    requested_times_ps: Sequence[float] | None = None,
) -> Path:
    """Plot stored SS power-density maps at the common snapshot times."""

    stored_times = np.asarray(dataset.get("snapshot_t_ps", []), dtype=float)
    indices = _nearest_unique_snapshot_indices(stored_times, requested_times_ps)
    joule_all = np.asarray(dataset.get("joule_snapshot_W_m3", []), dtype=float)
    diffusion_all = np.asarray(dataset.get("P_diff_snapshot_W_m3", []), dtype=float)
    if diffusion_all.shape != joule_all.shape:
        diffusion_all = np.zeros_like(joule_all)
    escape_all = np.asarray(dataset.get("P_esc_snapshot_W_m3", []), dtype=float)
    if escape_all.shape != joule_all.shape:
        escape_all = np.zeros_like(joule_all)
    channels = [
        (
            joule_all[indices],
            r"$P_J$ [W m$^{-3}$]",
            "magma",
            True,
        ),
        (
            np.asarray(dataset.get("P_total_snapshot_W_m3", []), dtype=float)[indices],
            r"$P_{e\mathrm{-}ph}=P_S+P_R$ [W m$^{-3}$]",
            "coolwarm",
            False,
        ),
        (
            diffusion_all[indices],
            r"$P_{\mathrm{diff}}$ [W m$^{-3}$]",
            "PuOr_r",
            False,
        ),
        (
            escape_all[indices],
            r"$P_{\mathrm{esc}}$ [W m$^{-3}$]",
            "cividis",
            True,
        ),
    ]
    return _plot_snapshot_scalar_atlas(
        dataset=dataset,
        output_path=output_path,
        times_ps=stored_times[indices],
        channels=channels,
        dpi=dpi,
    )


def plot_ss_energy_heat_capacity_snapshots(
    dataset: Mapping[str, Any],
    output_path: str | Path,
    *,
    dpi: int = THESIS_DPI,
    requested_times_ps: Sequence[float] | None = None,
) -> Path:
    """Plot electronic/phononic energy densities and volumetric heat capacities."""

    stored_times = np.asarray(dataset.get("snapshot_t_ps", []), dtype=float)
    indices = _nearest_unique_snapshot_indices(stored_times, requested_times_ps)
    channels = [
        (
            np.asarray(dataset.get("u_e_snapshot_J_m3", []), dtype=float)[indices],
            r"$u_e$ [J m$^{-3}$]",
            "coolwarm",
            False,
        ),
        (
            np.asarray(dataset.get("u_ph_snapshot_J_m3", []), dtype=float)[indices],
            r"$u_{ph}$ [J m$^{-3}$]",
            "viridis",
            True,
        ),
        (
            np.asarray(dataset.get("C_e_snapshot_J_m3_K", []), dtype=float)[indices],
            r"$C_e$ [J m$^{-3}$ K$^{-1}$]",
            "magma",
            True,
        ),
        (
            np.asarray(dataset.get("C_ph_snapshot_J_m3_K", []), dtype=float)[indices],
            r"$C_{ph}$ [J m$^{-3}$ K$^{-1}$]",
            "cividis",
            True,
        ),
    ]
    return _plot_snapshot_scalar_atlas(
        dataset=dataset,
        output_path=output_path,
        times_ps=stored_times[indices],
        channels=channels,
        dpi=dpi,
    )


def plot_phasecg_numerical_procedure(
    dataset: Mapping[str, Any],
    output_path: str | Path,
    *,
    dpi: int = THESIS_DPI,
) -> Path:
    """Plot strictly computational solver procedure diagnostics."""

    output = _prepare_output(output_path)
    t = np.asarray(dataset.get("t_ps", []), dtype=float)
    fig, axes = plt.subplots(2, 2, figsize=(THESIS_WIDTH_IN, 5.15), sharex=True)
    fig.subplots_adjust(left=0.105, right=0.895, bottom=0.095, top=0.955, wspace=0.34, hspace=0.24)

    ax = axes[0, 0]
    _plot_binned_envelope(ax, t, dataset.get("dt_accepted_fs"), label="Accepted median", color="tab:blue", band_label="Accepted 5-95%")
    _plot_binned_line(ax, t, dataset.get("dt_next_fs"), label="Next median", color="tab:orange")
    _plot_binned_line(ax, t, dataset.get("dt_attempt_fs"), label="Attempted median", color="tab:green")
    ax.set_yscale("log")
    ax.set_ylabel(r"$\Delta t$ [fs]")
    ax.set_title("Adaptive time step")
    ax.legend(frameon=False, fontsize=7.0, loc="best")

    ax = axes[0, 1]
    _plot_binned_line(ax, t, dataset.get("solve_attempts_per_step"), label="Solve attempts", color="tab:blue")
    _plot_binned_line(ax, t, dataset.get("adaptive_retries"), label="Retries", color="tab:orange")
    rejected_axis = ax.twinx()
    _plot_decimated(rejected_axis, t, dataset.get("cumulative_rejected_attempts"), "Cumulative rejected", color="tab:red")
    ax.set_ylabel("Per accepted step")
    rejected_axis.set_ylabel("Cumulative rejected", color="tab:red")
    rejected_axis.tick_params(axis="y", colors="tab:red")
    ax.set_title("Nonlinear solve effort")
    _combined_legend(ax, rejected_axis, fontsize=7.0, loc="upper left")

    ax = axes[1, 0]
    cg_tol = float(dataset.get("phase_convergence_tolerance", np.nan))
    poisson_tol = float(dataset.get("poisson_tolerance", np.nan))
    cg_values = np.asarray(dataset.get("allmaras_phase_convergence_residual_rel", []), dtype=float)
    poisson_values = np.asarray(dataset.get("poisson_residual_rel", []), dtype=float)
    if np.isfinite(cg_tol) and cg_tol > 0.0:
        _plot_binned_envelope(ax, t, cg_values / cg_tol, label="Phase-CG median", color="tab:blue", band_label="Phase-CG 5-95%")
    if np.isfinite(poisson_tol) and poisson_tol > 0.0:
        _plot_binned_line(ax, t, poisson_values / poisson_tol, label="Poisson median", color="tab:orange")
    ax.axhline(1.0, color="0.25", linestyle="--", linewidth=0.9, label="Tolerance")
    iteration_axis = ax.twinx()
    _plot_binned_line(iteration_axis, t, dataset.get("allmaras_phase_convergence_iterations"), label="CG iterations", color="tab:red")
    ax.set_yscale("log")
    ax.set_ylabel("Residual / tolerance")
    iteration_axis.tick_params(axis="y", right=False, labelright=False)
    ax.set_title("Linear-solver convergence")
    _combined_legend(ax, iteration_axis, fontsize=7.0, loc="best")

    ax = axes[1, 1]
    direct = np.asarray(dataset.get("allmaras_phase_direct_node_count", []), dtype=float)
    continued = np.asarray(dataset.get("allmaras_phase_continued_node_count", []), dtype=float)
    zero = np.asarray(dataset.get("allmaras_phase_zero_amplitude_node_count", []), dtype=float)
    total = np.maximum(direct + continued + zero, 1.0)
    _plot_binned_line(ax, t, 100.0 * direct / total, label="Direct", color="tab:blue")
    _plot_binned_line(ax, t, 100.0 * continued / total, label="Continued", color="tab:orange")
    _plot_binned_line(ax, t, 100.0 * zero / total, label="Zero amplitude", color="tab:green")
    ax.set_ylabel("Node fraction [%]")
    ax.set_ylim(-2.0, 102.0)
    ax.set_title("Phase-continuation domain")
    ax.legend(frameon=False, fontsize=7.0, loc="best")

    for row, ax_row in enumerate(axes):
        for axis in ax_row:
            if row == 1:
                axis.set_xlabel(r"$t$ [ps]")
            else:
                axis.tick_params(axis="x", labelbottom=False)
            axis.grid(True)
            axis.set_xlim(left=0.0)
    fig.savefig(output, dpi=dpi)
    plt.close(fig)
    return output


def plot_phasecg_stationarity_evolution(
    dataset: Mapping[str, Any],
    output_path: str | Path,
    *,
    dpi: int = THESIS_DPI,
) -> Path:
    """Plot current photon-readiness tolerances and reconstructed gate history."""

    output = _prepare_output(output_path)
    t = np.asarray(dataset.get("stationarity_eval_t_ps", []), dtype=float)
    if t.size == 0:
        raise ValueError("Stationarity reanalysis did not provide evaluation times.")
    fig, axes = plt.subplots(2, 2, figsize=(THESIS_WIDTH_IN, 5.25), sharex=True)
    fig.subplots_adjust(left=0.125, right=0.875, bottom=0.135, top=0.950, wspace=0.30, hspace=0.25)

    _plot_margin(axes[0, 0], t, dataset.get("strict_q_tolerance_margin"), r"Strict $Q$", "tab:blue")
    _plot_margin(axes[0, 0], t, dataset.get("strict_phi_tolerance_margin"), r"Strict $\nabla\phi$", "tab:orange")
    _plot_margin(axes[0, 0], t, dataset.get("dynamic_profile_tolerance_margin"), "Dynamic profile", "tab:green")
    _plot_margin(axes[0, 0], t, dataset.get("dynamic_voltage_tolerance_margin"), "Dynamic voltage", "tab:purple")
    _format_tolerance_axis(axes[0, 0], title="Mesoscopic stationarity")

    _plot_margin(axes[0, 1], t, dataset.get("continuity_rms_tolerance_margin"), "Continuity RMS", "tab:blue")
    _plot_margin(axes[0, 1], t, dataset.get("continuity_max_tolerance_margin"), "Continuity maximum", "tab:orange")
    _plot_margin(axes[0, 1], t, dataset.get("poisson_tolerance_margin"), "Poisson", "tab:green")
    _plot_margin(axes[0, 1], t, dataset.get("phase_cg_tolerance_margin"), "Phase CG", "tab:red")
    _format_tolerance_axis(axes[0, 1], title="Validity and phase solve")

    _plot_margin(axes[1, 0], t, dataset.get("thermal_relative_tolerance_margin"), "Thermal relative", "tab:blue")
    _plot_margin(axes[1, 0], t, dataset.get("thermal_p99_tolerance_margin"), "Thermal p99", "tab:orange")
    _plot_margin(axes[1, 0], t, dataset.get("thermal_projected_tolerance_margin"), "Thermal projected", "tab:green")
    _plot_margin(axes[1, 0], t, dataset.get("circuit_value_tolerance_margin"), "Circuit values", "tab:purple")
    _plot_margin(axes[1, 0], t, dataset.get("circuit_rhs_tolerance_margin"), "Circuit RHS", "tab:red")
    _format_tolerance_axis(axes[1, 0], title="Thermal and circuit tails")

    gate_labels = (
        "Strict", "Dynamic", "Mesoscopic", "Contact",
        "Continuity", "Thermal", "Circuit", "Phase CG", "Photon ready",
    )
    gate_keys = (
        "strict_stationarity_pass_history", "dynamic_stationarity_pass_history",
        "mesoscopic_stationarity_pass_history", "contact_recovery_pass_history",
        "continuity_pass_history", "thermal_stationarity_pass_history",
        "circuit_stationarity_pass_history", "phase_cg_pass_history", "photon_ready_history",
    )
    gates = np.vstack([np.asarray(dataset.get(key, np.zeros(t.size)), dtype=bool) for key in gate_keys])
    ax = axes[1, 1]
    right = float(t[-1]) if t.size > 1 else float(t[0]) + 1.0
    ax.imshow(gates, aspect="auto", interpolation="nearest", origin="lower", extent=(float(t[0]), right, -0.5, len(gate_labels) - 0.5), cmap=ListedColormap(["0.88", "#2a9d62"]), vmin=0.0, vmax=1.0)
    ax.set_yticks(np.arange(len(gate_labels)))
    ax.set_yticklabels(gate_labels, fontsize=7.0)
    ax.yaxis.tick_right()
    ax.tick_params(axis="y", labelleft=False, labelright=True, pad=2.0)
    ax.set_title("Photon-readiness gates")
    ax.text(0.99, 0.02, "gray: fail   green: pass", transform=ax.transAxes, ha="right", va="bottom", fontsize=6.8)

    summary = dict(dataset.get("photon_ready_reanalysis_summary", {}))
    first_ready = summary.get("first_ready_time_ps")
    if first_ready is not None and np.isfinite(float(first_ready)):
        for axis in axes.ravel():
            axis.axvline(float(first_ready), color="0.2", linestyle=":", linewidth=0.9)
        ax.text(float(first_ready), len(gate_labels) - 0.6, rf" ready at {float(first_ready):.1f} ps", ha="left", va="top", fontsize=7.0)

    for row, ax_row in enumerate(axes):
        for axis in ax_row:
            if row == 1:
                axis.set_xlabel(r"$t$ [ps]")
            else:
                axis.tick_params(axis="x", labelbottom=False)
            axis.set_xlim(left=0.0)
            if axis is not axes[1, 1]:
                axis.grid(True)
    status = "PASS" if bool(summary.get("passes", False)) else "NOT REACHED"
    fig.text(0.5, 0.025, f"Current-policy reanalysis: {status}; stored legacy photon_ready: {summary.get('stored_photon_ready')!s}", ha="center", fontsize=7.2)
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
    ax.set_ylabel("Count per accepted step [count]")
    rejected_axis.tick_params(axis="y", labelright=False, right=False)
    ax.set_title("Nonlinear solve effort")
    _combined_legend(ax, rejected_axis, fontsize=7.0)

    ax = axes[0, 2]
    _plot_decimated(ax, t, dataset.get("thermal_max_Te_K_history"), r"Max $T_e$")
    _plot_decimated(ax, t, dataset.get("thermal_max_Tph_K_history"), r"Max $T_{ph}$")
    rate_axis = ax.twinx()
    _plot_decimated(rate_axis, t, dataset.get("thermal_max_rate_K_per_ps_history"), "Max thermal rate", positive=True, color="tab:red")
    ax.set_ylabel("Temperature [K]")
    rate_axis.set_ylabel(r"Thermal rate [K ps$^{-1}$]", color="tab:red")
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
    ax.set_ylabel("Relative residual [-]")
    iteration_axis.tick_params(axis="y", labelright=False, right=False)
    ax.set_title("Harmonic continuation")
    _combined_legend(ax, iteration_axis, fontsize=7.0)

    ax = axes[2, 1]
    _plot_decimated(ax, t, dataset.get("allmaras_phase_direct_node_count"), "Direct nodes")
    _plot_decimated(ax, t, dataset.get("allmaras_phase_continued_node_count"), "Continued nodes")
    _plot_decimated(ax, t, dataset.get("allmaras_phase_zero_amplitude_node_count"), "Zero-amplitude nodes")
    ax.set_ylabel("Nodes [count]")
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
    mismatch_axis.set_ylabel("Relative mismatch [-]", color="tab:red")
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


def _plot_snapshot_scalar_atlas(
    *,
    dataset: Mapping[str, Any],
    output_path: str | Path,
    times_ps: np.ndarray,
    channels: Sequence[tuple[np.ndarray, str, str, bool]],
    dpi: int,
) -> Path:
    output = _prepare_output(output_path)
    tri = _triangulation(dataset)
    times = np.asarray(times_ps, dtype=float)
    n_rows = times.size
    n_cols = len(channels)
    if n_rows == 0:
        raise ValueError("No stored times are available for the snapshot atlas.")
    fig, axes = plt.subplots(
        n_rows,
        n_cols,
        figsize=(THESIS_WIDTH_IN, max(4.9, 0.50 * n_rows + 1.0)),
        squeeze=False,
    )
    fig.subplots_adjust(left=0.078, right=0.948, bottom=0.070, top=0.900, wspace=0.08, hspace=0.10)
    column_mappables = []
    for column, (values, label, cmap, positive) in enumerate(channels):
        field = np.asarray(values, dtype=float)
        if field.shape != (n_rows, tri.x.size):
            raise ValueError(f"Snapshot channel {label} has shape {field.shape}, expected {(n_rows, tri.x.size)}.")
        norm = _scalar_field_norm(field, positive=positive)
        for row in range(n_rows):
            axis = axes[row, column]
            mappable = axis.tripcolor(
                tri,
                field[row],
                shading="gouraud",
                cmap=cmap,
                norm=norm,
                rasterized=True,
            )
            _format_strip_axis(axis, tri)
            if row < n_rows - 1:
                axis.tick_params(axis="x", labelbottom=False)
            if column != 0:
                axis.tick_params(axis="y", labelleft=False)
            if column == n_cols - 1:
                axis.text(1.05, 0.5, rf"$t={times[row]:.3g}$ [ps]", transform=axis.transAxes, rotation=-90, va="center", ha="left", fontsize=7.2)
        column_mappables.append((mappable, label, field))
    fig.supxlabel(r"$x$ [nm]", y=0.017, fontsize=8.5)
    fig.supylabel(r"$y$ [nm]", x=0.016, fontsize=8.5)
    fig.canvas.draw()
    for column, (mappable, label, field) in enumerate(column_mappables):
        position = axes[0, column].get_position()
        color_axis = fig.add_axes([position.x0, position.y1 + 0.009, position.width, 0.010])
        colorbar = fig.colorbar(mappable, cax=color_axis, orientation="horizontal")
        color_axis.xaxis.set_ticks_position("top")
        color_axis.xaxis.set_label_position("top")
        colorbar.set_label(label, labelpad=1.5, fontsize=7.8)
        _set_compact_colorbar_ticks(colorbar, field)
        color_axis.tick_params(labelsize=6.6, pad=0.8, length=2.0)
        tick_labels = colorbar.ax.get_xticklabels()
        if len(tick_labels) >= 2:
            tick_labels[0].set_ha("left")
            tick_labels[-1].set_ha("right")
    fig.savefig(output, dpi=dpi)
    plt.close(fig)
    return output


def _scalar_field_norm(values: np.ndarray, *, positive: bool):
    finite = np.asarray(values, dtype=float)
    finite = finite[np.isfinite(finite)]
    if finite.size == 0 or float(np.nanmax(finite) - np.nanmin(finite)) == 0.0:
        center = float(finite[0]) if finite.size else 0.0
        width = max(abs(center) * 0.05, 1.0)
        return Normalize(vmin=center - width, vmax=center + width)
    if positive or float(np.nanmin(finite)) >= 0.0:
        positive_values = finite[finite > 0.0]
        vmax = float(np.nanpercentile(finite, 99.7))
        if positive_values.size and vmax / max(float(np.nanpercentile(positive_values, 1.0)), 1.0e-300) >= 30.0:
            vmin = max(float(np.nanpercentile(positive_values, 1.0)), vmax * 1.0e-8)
            return LogNorm(vmin=vmin, vmax=max(vmax, vmin * 10.0))
        return Normalize(vmin=max(0.0, float(np.nanpercentile(finite, 0.3))), vmax=vmax)
    vmax = max(float(np.nanpercentile(np.abs(finite), 99.7)), 1.0e-300)
    return SymLogNorm(linthresh=max(vmax * 1.0e-5, 1.0e-300), vmin=-vmax, vmax=vmax, base=10.0)


def _set_compact_colorbar_ticks(colorbar, values: np.ndarray) -> None:
    """Keep narrow atlas colorbars legible without implying nonexistent ranges."""

    finite = np.asarray(values, dtype=float)
    finite = finite[np.isfinite(finite)]
    if finite.size == 0:
        return
    lo = float(np.nanmin(finite))
    hi = float(np.nanmax(finite))
    if np.isclose(lo, hi, rtol=1.0e-10, atol=max(abs(lo), abs(hi), 1.0) * 1.0e-12):
        ticks = [0.5 * (lo + hi)]
    else:
        norm = colorbar.mappable.norm
        norm_lo = float(norm.vmin)
        norm_hi = float(norm.vmax)
        if isinstance(norm, LogNorm) and norm_lo > 0.0:
            ticks = [norm_lo, float(np.sqrt(norm_lo * norm_hi)), norm_hi]
        elif norm_lo < 0.0 < norm_hi:
            ticks = [norm_lo, 0.0, norm_hi]
        else:
            ticks = [norm_lo, 0.5 * (norm_lo + norm_hi), norm_hi]
    colorbar.set_ticks(ticks)
    colorbar.set_ticklabels([_compact_number(value) for value in ticks])


def _compact_number(value: float) -> str:
    magnitude = abs(float(value))
    if magnitude == 0.0:
        return "0"
    if magnitude >= 1.0e3 or magnitude < 1.0e-2:
        return f"{value:.1e}"
    return f"{value:.3g}"


def _weighted_rows_for_plot(values: np.ndarray, weights: np.ndarray) -> np.ndarray:
    matrix = np.asarray(values, dtype=float)
    weight = np.asarray(weights, dtype=float).reshape(-1)
    if matrix.ndim != 2 or matrix.shape[1] != weight.size:
        return np.nanmean(matrix, axis=1)
    valid_weight = np.where(np.isfinite(weight) & (weight > 0.0), weight, 0.0)
    finite = np.isfinite(matrix)
    numerator = np.sum(np.where(finite, matrix, 0.0) * valid_weight[None, :], axis=1)
    denominator = np.sum(finite * valid_weight[None, :], axis=1)
    return numerator / np.maximum(denominator, 1.0e-300)


def _add_current_direction_arrow(
    ax: plt.Axes,
    tri: mtri.Triangulation,
    *,
    mean_x: float,
    mean_y: float,
    reference_magnitude: float,
) -> None:
    magnitude = float(np.hypot(mean_x, mean_y))
    if not np.isfinite(magnitude) or magnitude <= 0.0:
        return
    ux = mean_x / magnitude
    uy = mean_y / magnitude
    length = float(np.ptp(tri.x))
    width = float(np.ptp(tri.y))
    arrow_length = 0.20 * length * min(magnitude / max(reference_magnitude, 1.0e-300), 1.0)
    if abs(uy) > 1.0e-12:
        arrow_length = min(arrow_length, 0.65 * width / abs(uy))
    center_x = 0.5 * (float(np.nanmin(tri.x)) + float(np.nanmax(tri.x)))
    center_y = 0.5 * (float(np.nanmin(tri.y)) + float(np.nanmax(tri.y)))
    start = (center_x - 0.5 * ux * arrow_length, center_y - 0.5 * uy * arrow_length)
    stop = (center_x + 0.5 * ux * arrow_length, center_y + 0.5 * uy * arrow_length)
    arrow = FancyArrowPatch(start, stop, arrowstyle="-|>", mutation_scale=8.5, linewidth=1.2, color="white", zorder=8)
    arrow.set_path_effects([path_effects.Stroke(linewidth=2.2, foreground="0.1"), path_effects.Normal()])
    ax.add_patch(arrow)


def _binned_statistics(x: Any, y: Any, *, n_bins: int = 320) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    x_values = np.asarray(x, dtype=float).reshape(-1)
    y_values = np.asarray(y if y is not None else [], dtype=float).reshape(-1)
    n = min(x_values.size, y_values.size)
    if n == 0:
        empty = np.array([], dtype=float)
        return empty, empty, empty, empty
    x_values = x_values[:n]
    y_values = y_values[:n]
    finite = np.isfinite(x_values) & np.isfinite(y_values)
    x_values = x_values[finite]
    y_values = y_values[finite]
    if x_values.size == 0:
        empty = np.array([], dtype=float)
        return empty, empty, empty, empty
    bins = min(max(1, int(n_bins)), x_values.size)
    edges = np.linspace(float(np.nanmin(x_values)), float(np.nanmax(x_values)), bins + 1)
    index = np.clip(np.digitize(x_values, edges) - 1, 0, bins - 1)
    centers = 0.5 * (edges[:-1] + edges[1:])
    median = np.full(bins, np.nan)
    low = np.full(bins, np.nan)
    high = np.full(bins, np.nan)
    for bin_index in range(bins):
        selected = y_values[index == bin_index]
        if selected.size:
            median[bin_index] = np.nanmedian(selected)
            low[bin_index], high[bin_index] = np.nanpercentile(selected, [5.0, 95.0])
    valid = np.isfinite(median)
    return centers[valid], median[valid], low[valid], high[valid]


def _plot_binned_line(ax: plt.Axes, x: Any, y: Any, *, label: str, color: str) -> None:
    centers, median, _, _ = _binned_statistics(x, y)
    if centers.size:
        ax.plot(centers, median, color=color, label=label, linewidth=1.0)


def _plot_binned_envelope(
    ax: plt.Axes,
    x: Any,
    y: Any,
    *,
    label: str,
    color: str,
    band_label: str,
) -> None:
    centers, median, low, high = _binned_statistics(x, y)
    if centers.size:
        ax.fill_between(centers, low, high, color=color, alpha=0.18, linewidth=0.0, label=band_label)
        ax.plot(centers, median, color=color, label=label, linewidth=1.0)


def _plot_margin(ax: plt.Axes, x: np.ndarray, values: Any, label: str, color: str) -> None:
    y = np.asarray(values if values is not None else [], dtype=float).reshape(-1)
    n = min(np.asarray(x).size, y.size)
    if n:
        finite = np.isfinite(y[:n]) & (y[:n] > 0.0)
        if np.any(finite):
            ax.plot(np.asarray(x)[:n][finite], y[:n][finite], color=color, label=label, linewidth=1.0)


def _format_tolerance_axis(ax: plt.Axes, *, title: str) -> None:
    ax.axhline(1.0, color="0.25", linestyle="--", linewidth=0.9, label="Tolerance")
    ax.set_yscale("log")
    ax.set_ylabel("Metric / tolerance")
    ax.set_title(title)
    ax.legend(frameon=False, fontsize=6.8, ncol=2, loc="best")


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
    loc: str = "best",
) -> None:
    handles_a, labels_a = ax.get_legend_handles_labels()
    handles_b, labels_b = twin.get_legend_handles_labels()
    handles = handles_a + handles_b
    labels = labels_a + labels_b
    if handles:
        ax.legend(handles, labels, frameon=False, ncol=ncol, fontsize=fontsize, loc=loc)


__all__ = [
    "make_phasecg_ss_figures",
    "plot_final_longitudinal_profiles",
    "plot_phasecg_numerical_procedure",
    "plot_phasecg_stationarity_evolution",
    "plot_phasecg_physical_evolution",
    "plot_phasecg_snapshot_fields",
    "plot_ss_energy_heat_capacity_snapshots",
    "plot_ss_power_density_snapshots",
]
