"""Thesis diagnostics for one completed photon/circuit transient."""

from __future__ import annotations

from pathlib import Path
from typing import Any, Mapping, Sequence

import matplotlib

matplotlib.use("Agg", force=True)
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import numpy as np

from pysnspd.plotting.photon_figures import phase_gradient_q_abs_m_inv
from pysnspd.plotting.ss_phasecg_figures import (
    plot_ss_energy_heat_capacity_snapshots,
    plot_ss_power_density_snapshots,
)
from pysnspd.plotting.style import THESIS_DPI, THESIS_WIDTH_IN, apply_thesis_style

apply_thesis_style()

MEV_J = 1.602176634e-22


def make_photon_run_diagnostic_figures(
    *,
    mesh: Any,
    history: Mapping[str, Any],
    snapshots: Mapping[str, Any],
    summary: Mapping[str, Any],
    delta0_meV: float,
    xi_m: float,
    requested_times_ps: Sequence[float],
    output_dir: str | Path,
    dpi: int = THESIS_DPI,
    timing: Mapping[str, Any] | None = None,
    snapshot_diagnostics: Mapping[str, Any] | None = None,
) -> dict[str, Path]:
    """Create scalar-evolution and selected-field figures for one photon run."""

    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)
    diagnostics = dict(snapshot_diagnostics or {})
    saved = {
        "scalar_evolution": plot_photon_scalar_evolution(
            history=history,
            summary=summary,
            timing=timing or {},
            output_path=out / "E3_photon_scalar_evolution.pdf",
            dpi=dpi,
        ),
        "field_evolution": plot_photon_field_evolution(
            mesh=mesh,
            snapshots=snapshots,
            summary=summary,
            delta0_meV=delta0_meV,
            xi_m=xi_m,
            requested_times_ps=requested_times_ps,
            snapshot_diagnostics=diagnostics,
            output_path=out / "E3_photon_field_evolution.pdf",
            dpi=dpi,
        ),
    }
    if diagnostics:
        saved["power_density_snapshots"] = plot_ss_power_density_snapshots(
            diagnostics,
            out / "E3_photon_power_density_snapshots.pdf",
            dpi=dpi,
        )
        saved["energy_heat_capacity_snapshots"] = plot_ss_energy_heat_capacity_snapshots(
            diagnostics,
            out / "E3_photon_energy_heat_capacity_snapshots.pdf",
            dpi=dpi,
        )

    recovery_path = out / "E3_photon_censored_recovery_diagnostics.pdf"
    if _is_detected_but_unrecovered(timing or {}):
        saved["censored_recovery_diagnostics"] = plot_photon_censored_recovery_diagnostics(
            history=history,
            summary=summary,
            timing=timing or {},
            output_path=recovery_path,
            dpi=dpi,
        )
    elif recovery_path.exists():
        recovery_path.unlink()
    return saved


def plot_photon_scalar_evolution(
    *,
    history: Mapping[str, Any],
    summary: Mapping[str, Any],
    timing: Mapping[str, Any],
    output_path: str | Path,
    dpi: int = THESIS_DPI,
) -> Path:
    """Plot the complete circuit, thermal and condensate histories."""

    output = _prepare_output(output_path)
    time = _history_time_ps(history)
    if time.size < 2:
        raise ValueError("The photon history must contain at least two time samples.")

    current_s = 1.0e6 * _series(history, "I_s_A", time.size)
    current_rf = _series(history, "I_rf_A", time.size)
    if not np.isfinite(current_rf).any():
        current_rf = _series(history, "I_b_A", time.size) - _series(
            history, "I_s_A", time.size
        )
    current_rf = 1.0e6 * current_rf

    voltage_tdgl = 1.0e3 * _series(history, "V_tdgl_center_V", time.size)
    voltage_out = _series(history, "V_out_V", time.size)
    if not np.isfinite(voltage_out).any():
        resistance = _nested_float(
            summary,
            ("circuit", "params", "R_load_ohm"),
            default=50.0,
        )
        voltage_out = resistance * current_rf * 1.0e-6
    voltage_out = 1.0e3 * voltage_out

    fig, axes = plt.subplots(4, 1, figsize=(THESIS_WIDTH_IN, 6.3), sharex=True)
    fig.subplots_adjust(left=0.115, right=0.870, bottom=0.080, top=0.985, hspace=0.18)

    axes[0].plot(time, current_s - float(current_s[0]), label=r"$\Delta I_s$")
    axes[0].plot(time, current_rf, linestyle="--", label=r"$I_{\mathrm{RF}}$")
    axes[0].set_ylabel(r"Current [$\mu$A]")

    axes[1].plot(time, voltage_tdgl, label=r"$V_{\mathrm{TDGL}}$")
    axes[1].plot(time, voltage_out, linestyle="--", label=r"$V_{\mathrm{out}}$")
    axes[1].set_ylabel("Voltage [mV]")

    axes[2].plot(time, _series(history, "max_Te_K", time.size), label=r"Max $T_e$")
    axes[2].plot(
        time,
        _series(history, "max_Tph_K", time.size),
        linestyle="--",
        label=r"Max $T_{ph}$",
    )
    axes[2].set_ylabel("Temperature [K]")

    axes[3].plot(
        time,
        _series(history, "mean_delta_over_delta0", time.size),
        label=r"Mean $|\Delta|/\Delta_0$",
    )
    axes[3].plot(
        time,
        _series(history, "min_delta_over_delta0", time.size),
        linestyle="--",
        label=r"Min $|\Delta|/\Delta_0$",
    )
    pairbreaking = _series(history, "max_pairbreaking_ratio", time.size)
    pair_axis = None
    if np.isfinite(pairbreaking).any():
        pair_axis = axes[3].twinx()
        pair_axis.plot(
            time,
            pairbreaking,
            color="tab:red",
            linestyle=":",
            label="Max pair-breaking ratio",
        )
        pair_axis.set_ylabel("Pair breaking", color="tab:red")
        pair_axis.tick_params(axis="y", colors="tab:red")
        positive_pairbreaking = pairbreaking[
            np.isfinite(pairbreaking) & (pairbreaking > 0.0)
        ]
        if (
            positive_pairbreaking.size
            and float(np.nanmax(positive_pairbreaking))
            / max(float(np.nanmin(positive_pairbreaking)), 1.0e-300)
            > 1.0e3
        ):
            pair_axis.set_yscale("log")
    axes[3].set_ylabel(r"$|\Delta|/\Delta_0$")
    axes[3].set_xlabel(r"$t$ [ps]")

    event_time = _photon_time_ps(history, summary)
    latency = dict(timing.get("latency", {}))
    recovery = dict(dict(timing.get("recovery", {})).get("selected", {}))
    markers = (
        (event_time, "0.25", ":", "Photon arrival"),
        (latency.get("crossing_time_ps"), "tab:blue", "--", "Detection"),
        (recovery.get("entry_time_ps"), "tab:green", "-.", "Recovery"),
    )
    for index, axis in enumerate(axes):
        for value, color, style, label in markers:
            if value is not None and np.isfinite(float(value)):
                axis.axvline(
                    float(value),
                    color=color,
                    linestyle=style,
                    linewidth=0.9,
                    label=label if index == 0 else None,
                )
        axis.grid(True)
        axis.set_xlim(left=0.0)

    for index, axis in enumerate(axes):
        if index == 3 and pair_axis is not None:
            _combined_legend(
                axis,
                pair_axis,
                loc="lower center",
                bbox_to_anchor=(0.48, 0.20),
            )
        else:
            axis.legend(frameon=False, ncol=3 if index == 0 else 2, loc="best")

    timing_text = (
        f"t_lat={_timing_value(latency, 't_lat_ps')}; "
        f"t_rec[{recovery.get('mode', 'electrical')}]={_timing_value(recovery, 't_rec_ps')}"
    )
    axes[3].text(
        0.99,
        0.96,
        timing_text,
        transform=axes[3].transAxes,
        ha="right",
        va="top",
        fontsize=7.0,
        bbox={"facecolor": "white", "edgecolor": "0.75", "alpha": 0.82, "pad": 1.5},
    )

    fig.savefig(output, dpi=dpi)
    plt.close(fig)
    return output


def plot_photon_field_evolution(
    *,
    mesh: Any,
    snapshots: Mapping[str, Any],
    summary: Mapping[str, Any],
    delta0_meV: float,
    xi_m: float,
    requested_times_ps: Sequence[float],
    snapshot_diagnostics: Mapping[str, Any] | None = None,
    output_path: str | Path,
    dpi: int = THESIS_DPI,
) -> Path:
    """Plot selected photon fields, including dimensionless superfluid momentum."""

    output = _prepare_output(output_path)
    if not np.isfinite(delta0_meV) or delta0_meV <= 0.0:
        raise ValueError("delta0_meV must be finite and positive.")
    if not np.isfinite(xi_m) or xi_m <= 0.0:
        raise ValueError("xi_m must be finite and positive.")

    tri = _triangulation(mesh)
    times = _snapshot_times_ps(snapshots)
    indices = nearest_unique_snapshot_indices(times, requested_times_ps)
    selected_times = times[indices]
    n_nodes = tri.x.size

    delta_all = _snapshot_matrix(snapshots, "delta_snapshot_meV", n_nodes)
    phi_all = _snapshot_matrix(snapshots, "phi_snapshot_V", n_nodes)
    Te_all = _snapshot_matrix(snapshots, "Te_snapshot_K", n_nodes)
    Tph_all = _snapshot_matrix(snapshots, "Tph_snapshot_K", n_nodes)
    delta = delta_all[indices] / delta0_meV
    phi = 1.0e3 * phi_all[indices]
    Te = Te_all[indices]
    Tph = Tph_all[indices]
    real = _snapshot_matrix(snapshots, "psi_real_snapshot_J", n_nodes, required=False)
    imag = _snapshot_matrix(snapshots, "psi_imag_snapshot_J", n_nodes, required=False)

    qxi, qxi_limits = _selected_qxi_and_global_limits(
        tri=tri,
        snapshots=snapshots,
        snapshot_diagnostics=snapshot_diagnostics or {},
        selected_indices=indices,
        xi_m=xi_m,
        real=real,
        imag=imag,
        n_nodes=n_nodes,
    )

    fields = (
        (delta, r"$|\Delta|/\Delta_{\mathrm{BCS}}(0)$", "viridis", False, 0.0, None, _finite_limits(delta_all) / delta0_meV),
        (phi, r"$\phi$ [mV]", "coolwarm", True, None, None, 1.0e3 * _finite_limits(phi_all)),
        (qxi, r"$|\mathbf{q}|\xi$", "magma", False, 0.0, None, qxi_limits),
        (Te, r"$T_e$ [K]", "inferno", False, None, None, _finite_limits(Te_all)),
        (Tph, r"$T_{ph}$ [K]", "inferno", False, None, None, _finite_limits(Tph_all)),
    )
    impact = _impact_coordinates_nm(summary)
    n_rows = len(indices)
    fig, axes = plt.subplots(
        n_rows,
        len(fields),
        figsize=(THESIS_WIDTH_IN, max(3.2, 0.55 * n_rows + 1.0)),
        squeeze=False,
    )
    fig.subplots_adjust(left=0.078, right=0.965, bottom=0.072, top=0.905, wspace=0.08, hspace=0.12)

    mappables = []
    for col, (values, label, cmap, symmetric, forced_min, forced_max, global_limits) in enumerate(fields):
        vmin, vmax = _field_limits(
            values,
            symmetric=symmetric,
            forced_min=forced_min,
            forced_max=forced_max,
            global_limits=global_limits,
        )
        for row in range(n_rows):
            axis = axes[row, col]
            mappable = axis.tripcolor(
                tri,
                values[row],
                shading="gouraud",
                cmap=cmap,
                vmin=vmin,
                vmax=vmax,
                rasterized=True,
            )
            _format_strip_axis(axis, tri)
            if np.isfinite(impact).all():
                axis.plot(
                    impact[0],
                    impact[1],
                    marker="x",
                    markersize=3.6,
                    markeredgewidth=0.8,
                    color="white",
                )
            if row < n_rows - 1:
                axis.tick_params(axis="x", labelbottom=False)
            if col != 0:
                axis.tick_params(axis="y", labelleft=False)
            if col == len(fields) - 1:
                axis.text(
                    1.05,
                    0.5,
                    rf"$t={selected_times[row]:.3g}$ [ps]",
                    transform=axis.transAxes,
                    rotation=-90,
                    va="center",
                    ha="left",
                    fontsize=7.2,
                )
        mappables.append((mappable, label))

    fig.supxlabel(r"$x$ [nm]", y=0.018, fontsize=8.5)
    fig.supylabel(r"$y$ [nm]", x=0.016, fontsize=8.5)
    fig.canvas.draw()
    for col, (mappable, label) in enumerate(mappables):
        position = axes[0, col].get_position()
        color_axis = fig.add_axes(
            [position.x0, position.y1 + 0.010, position.width, 0.011]
        )
        colorbar = fig.colorbar(mappable, cax=color_axis, orientation="horizontal")
        color_axis.xaxis.set_ticks_position("top")
        color_axis.xaxis.set_label_position("top")
        colorbar.set_label(label, labelpad=1.5, fontsize=8.4)
        color_axis.tick_params(labelsize=6.8, pad=0.8, length=2.0)

    fig.savefig(output, dpi=dpi)
    plt.close(fig)
    return output


def plot_photon_censored_recovery_diagnostics(
    *,
    history: Mapping[str, Any],
    summary: Mapping[str, Any],
    timing: Mapping[str, Any],
    output_path: str | Path,
    dpi: int = THESIS_DPI,
) -> Path:
    """Diagnose a detected transient whose selected recovery stayed censored."""

    output = _prepare_output(output_path)
    time_ps = _history_time_ps(history)
    event_ps = _photon_time_ps(history, summary)
    if time_ps.size < 2 or not np.isfinite(event_ps):
        raise ValueError("Censored-recovery diagnostics require a photon event and history.")
    post = np.isfinite(time_ps) & (time_ps >= event_ps)
    relative_time = time_ps[post] - event_ps
    baseline = dict(timing.get("baseline", {})).get("values", {})
    recovery = dict(dict(timing.get("recovery", {})).get("selected", {}))
    tolerances = dict(recovery.get("absolute_tolerances", {}))
    if not isinstance(baseline, Mapping) or not tolerances:
        raise ValueError("Electrical recovery baselines and tolerances are unavailable.")

    current_keys = ("I_b_A", "I_s_A", "I_rf_A")
    voltage_keys = ("V_out_V", "v_c_V", "V_tdgl_center_V")
    labels = {
        "I_b_A": r"$I_b$",
        "I_s_A": r"$I_s$",
        "I_rf_A": r"$I_{\mathrm{RF}}$",
        "V_out_V": r"$V_{\mathrm{out}}$",
        "v_c_V": r"$v_c$",
        "V_tdgl_center_V": r"$V_{\mathrm{TDGL}}$",
    }
    ratios: dict[str, np.ndarray] = {}
    for key in current_keys + voltage_keys:
        tolerance = float(tolerances.get(key, np.nan))
        reference = float(baseline.get(key, np.nan))
        values = _series(history, key, time_ps.size)
        ratios[key] = np.abs(values[post] - reference) / max(tolerance, 1.0e-300)

    fig, axes = plt.subplots(2, 2, figsize=(THESIS_WIDTH_IN, 5.35))
    fig.subplots_adjust(left=0.115, right=0.965, bottom=0.105, top=0.900, wspace=0.40, hspace=0.32)
    colors = ("tab:blue", "tab:orange", "tab:green")
    for axis, keys, title in (
        (axes[0, 0], current_keys, "Current recovery margins"),
        (axes[0, 1], voltage_keys, "Voltage recovery margins"),
    ):
        for color, key in zip(colors, keys):
            axis.plot(
                relative_time,
                np.maximum(ratios[key], 1.0e-6),
                color=color,
                linewidth=1.0,
                label=labels[key],
            )
        axis.axhline(1.0, color="0.2", linestyle="--", linewidth=0.9, label="Tolerance")
        axis.set_yscale("log")
        axis.set_xlabel(r"$t-t_\gamma$ [ps]")
        axis.set_ylabel("Residual / tolerance")
        axis.set_title(title)
        axis.legend(frameon=False, fontsize=7.2, ncol=2, loc="upper right")
        axis.grid(True)
        axis.set_xlim(left=0.0)

    final_keys = current_keys + voltage_keys
    final_ratios = np.asarray(
        [float(ratios[key][-1]) if ratios[key].size else np.nan for key in final_keys],
        dtype=float,
    )
    axis = axes[1, 0]
    y = np.arange(len(final_keys))
    bar_values = np.maximum(final_ratios, 1.0e-4)
    bar_colors = ["#d1495b" if value > 1.0 else "#2a9d62" for value in final_ratios]
    axis.barh(y, bar_values, color=bar_colors, alpha=0.88)
    axis.axvline(1.0, color="0.2", linestyle="--", linewidth=0.9)
    axis.set_xscale("log")
    axis.set_yticks(y)
    axis.set_yticklabels([labels[key] for key in final_keys])
    axis.invert_yaxis()
    axis.set_xlabel("Final residual / tolerance")
    axis.set_title("Distance from electrical recovery")
    axis.grid(True, axis="x")
    for row, value in enumerate(final_ratios):
        if np.isfinite(value):
            axis.text(
                max(value, 1.0e-4) * 1.08,
                row,
                f"{value:.2g}x",
                va="center",
                fontsize=7.0,
            )

    axis = axes[1, 1]
    mode_labels, mode_values, mode_colors = _circuit_timescale_bars(
        summary=summary,
        timing=timing,
        post_photon_window_ps=float(relative_time[-1]),
    )
    if mode_values:
        row = np.arange(len(mode_values))
        axis.barh(row, mode_values, color=mode_colors, alpha=0.88)
        axis.set_xscale("log")
        axis.set_yticks(row)
        axis.set_yticklabels(mode_labels, fontsize=7.2)
        axis.invert_yaxis()
        axis.set_xlabel("Timescale [ps]")
        axis.grid(True, axis="x")
        for index, value in enumerate(mode_values):
            axis.text(value * 1.08, index, f"{value:.3g}", va="center", fontsize=6.8)
        decay_values = [
            value for label, value in zip(mode_labels, mode_values) if "decay" in label
        ]
        if decay_values and np.isfinite(post_window := float(relative_time[-1])):
            axis.text(
                0.98,
                0.02,
                f"window / slowest decay = {post_window / max(decay_values):.3f}",
                transform=axis.transAxes,
                ha="right",
                va="bottom",
                fontsize=6.8,
            )
    else:
        axis.text(0.5, 0.5, "Circuit parameters unavailable", transform=axis.transAxes, ha="center", va="center")
        axis.set_xticks([])
        axis.set_yticks([])
    axis.set_title(
        "Overdamped modes vs available window"
        if mode_labels and not any("Osc." in label for label in mode_labels)
        else "Circuit modes vs available window"
    )

    latency = dict(timing.get("latency", {}))
    t_lat = latency.get("t_lat_ps")
    lower_bound = recovery.get("lower_bound_ps")
    title = (
        rf"Detected: $t_{{\mathrm{{lat}}}}={float(t_lat):.3g}$ ps; "
        rf"electrical recovery censored beyond {float(lower_bound):.3g} ps"
        if t_lat is not None and lower_bound is not None
        else "Detected transient with censored electrical recovery"
    )
    fig.suptitle(title, y=0.975, fontsize=10.0)
    fig.savefig(output, dpi=dpi)
    plt.close(fig)
    return output


def _is_detected_but_unrecovered(timing: Mapping[str, Any]) -> bool:
    latency = dict(timing.get("latency", {}))
    recovery = dict(dict(timing.get("recovery", {})).get("selected", {}))
    detected = bool(latency.get("detected", latency.get("t_lat_ps") is not None))
    recovered = bool(recovery.get("recovered", recovery.get("t_rec_ps") is not None))
    return bool(detected and not recovered)


def _circuit_timescale_bars(
    *,
    summary: Mapping[str, Any],
    timing: Mapping[str, Any],
    post_photon_window_ps: float,
) -> tuple[list[str], list[float], list[str]]:
    params = dict(dict(summary.get("circuit", {})).get("params", {}))
    try:
        Rl = float(params["R_load_ohm"])
        Rb = float(params["R_bias_ohm"])
        Lb = float(params["L_bias_H"])
        Lk = float(params["L_k_H"])
        capacitance = float(params["C_couple_F"])
    except (KeyError, TypeError, ValueError):
        return [], [], []
    if not all(np.isfinite(value) and value > 0.0 for value in (Rl, Rb, Lb, Lk, capacitance)):
        return [], [], []

    matrix = np.asarray(
        [
            [-(Rb + Rl) / Lb, Rl / Lb, -1.0 / Lb],
            [Rl / Lk, -Rl / Lk, 1.0 / Lk],
            [1.0 / capacitance, -1.0 / capacitance, 0.0],
        ],
        dtype=float,
    )
    eigenvalues = np.linalg.eigvals(matrix)
    real_modes = [value for value in eigenvalues if abs(float(np.imag(value))) <= 1.0e-8 * max(abs(value), 1.0)]
    complex_modes = [value for value in eigenvalues if float(np.imag(value)) > 0.0]
    labels: list[str] = []
    values: list[float] = []
    colors: list[str] = []
    for index, value in enumerate(sorted(real_modes, key=lambda item: abs(float(np.real(item))))):
        if float(np.real(value)) < 0.0:
            labels.append(f"Mode {index + 1}: decay")
            values.append(-1.0e12 / float(np.real(value)))
            colors.append("tab:blue")
    for index, value in enumerate(complex_modes):
        if float(np.real(value)) < 0.0:
            labels.append(f"Osc. mode {index + 1}: decay")
            values.append(-1.0e12 / float(np.real(value)))
            colors.append("tab:purple")
        labels.append(f"Osc. mode {index + 1}: period")
        values.append(2.0 * np.pi * 1.0e12 / abs(float(np.imag(value))))
        colors.append("tab:orange")
    hold_ps = 1.0e12 * float(dict(timing.get("recovery_criteria", {})).get("hold_s", np.nan))
    if np.isfinite(hold_ps) and hold_ps > 0.0:
        labels.append("Recovery hold")
        values.append(hold_ps)
        colors.append("tab:green")
    if np.isfinite(post_photon_window_ps) and post_photon_window_ps > 0.0:
        labels.append("Post-photon window")
        values.append(post_photon_window_ps)
        colors.append("0.35")
    return labels, values, colors


def nearest_unique_snapshot_indices(
    stored_times_ps: Sequence[float],
    requested_times_ps: Sequence[float],
) -> np.ndarray:
    """Resolve finite requested times to unique nearest stored indices."""

    stored = np.asarray(stored_times_ps, dtype=float).reshape(-1)
    if stored.size == 0:
        raise ValueError("The photon run has no stored snapshot times.")
    requested = np.asarray(list(requested_times_ps), dtype=float).reshape(-1)
    requested = requested[np.isfinite(requested)]
    if requested.size == 0:
        raise ValueError("At least one finite photon snapshot time is required.")
    indices: list[int] = []
    seen: set[int] = set()
    for value in requested:
        index = int(np.nanargmin(np.abs(stored - float(value))))
        if index not in seen:
            seen.add(index)
            indices.append(index)
    return np.asarray(indices, dtype=np.int64)


def _history_time_ps(history: Mapping[str, Any]) -> np.ndarray:
    if "t_ps" in history:
        return np.asarray(history["t_ps"], dtype=float).reshape(-1)
    if "t_s" in history:
        return np.asarray(history["t_s"], dtype=float).reshape(-1) / 1.0e-12
    return np.array([], dtype=float)


def _snapshot_times_ps(snapshots: Mapping[str, Any]) -> np.ndarray:
    if "snapshot_t_ps" in snapshots:
        return np.asarray(snapshots["snapshot_t_ps"], dtype=float).reshape(-1)
    if "snapshot_t_s" in snapshots:
        return np.asarray(snapshots["snapshot_t_s"], dtype=float).reshape(-1) / 1.0e-12
    return np.array([], dtype=float)


def _series(history: Mapping[str, Any], key: str, size: int) -> np.ndarray:
    values = np.asarray(history.get(key, []), dtype=float).reshape(-1)
    if values.size == 0:
        return np.full(size, np.nan, dtype=float)
    if values.size != size:
        values = np.resize(values, size)
    return values


def _snapshot_matrix(
    snapshots: Mapping[str, Any],
    key: str,
    n_nodes: int,
    *,
    required: bool = True,
) -> np.ndarray:
    values = np.asarray(snapshots.get(key, []), dtype=float)
    if values.ndim == 2 and values.shape[1] == n_nodes:
        return values
    if required:
        raise ValueError(f"Photon snapshots lack node-resolved field {key}.")
    return np.empty((0, n_nodes), dtype=float)


def _triangulation(mesh: Any) -> mtri.Triangulation:
    nodes = np.asarray(mesh.nodes, dtype=float)
    triangles = np.asarray(mesh.triangles, dtype=np.int64)
    return mtri.Triangulation(
        1.0e9 * nodes[:, 0],
        1.0e9 * nodes[:, 1],
        triangles,
    )


def _field_limits(
    values: np.ndarray,
    *,
    symmetric: bool,
    forced_min: float | None,
    forced_max: float | None,
    global_limits: np.ndarray | None = None,
) -> tuple[float, float]:
    limits = np.asarray(global_limits if global_limits is not None else [], dtype=float).reshape(-1)
    if limits.size == 2 and np.all(np.isfinite(limits)):
        finite = limits
    else:
        finite = np.asarray(values, dtype=float)
        finite = finite[np.isfinite(finite)]
    if finite.size == 0:
        return (-1.0, 1.0) if symmetric else (0.0, 1.0)
    if symmetric:
        limit = max(float(np.nanmax(np.abs(finite))), 1.0e-30)
        return -limit, limit
    lower = (
        float(forced_min)
        if forced_min is not None
        else float(np.nanmin(finite))
    )
    upper = (
        float(forced_max)
        if forced_max is not None
        else float(np.nanmax(finite))
    )
    if not np.isfinite(upper) or upper <= lower:
        upper = lower + 1.0
    return lower, upper


def _selected_qxi_and_global_limits(
    *,
    tri: mtri.Triangulation,
    snapshots: Mapping[str, Any],
    snapshot_diagnostics: Mapping[str, Any],
    selected_indices: np.ndarray,
    xi_m: float,
    real: np.ndarray,
    imag: np.ndarray,
    n_nodes: int,
) -> tuple[np.ndarray, np.ndarray]:
    diagnostic_indices = np.asarray(
        snapshot_diagnostics.get("selected_indices", []), dtype=np.int64
    ).reshape(-1)
    diagnostic_q = np.asarray(
        snapshot_diagnostics.get("q_abs_snapshot_m_inv", []), dtype=float
    )
    global_limits = snapshot_diagnostics.get("snapshot_global_limits", {})
    if (
        np.array_equal(diagnostic_indices, np.asarray(selected_indices, dtype=np.int64))
        and diagnostic_q.shape == (selected_indices.size, n_nodes)
        and isinstance(global_limits, Mapping)
        and "q_abs_snapshot_m_inv" in global_limits
    ):
        limits = np.asarray(global_limits["q_abs_snapshot_m_inv"], dtype=float).reshape(-1)
        if limits.size == 2 and np.all(np.isfinite(limits)):
            return xi_m * diagnostic_q, xi_m * limits

    q_selected = np.empty((selected_indices.size, n_nodes), dtype=float)
    selected_position = {
        int(index): position for position, index in enumerate(selected_indices)
    }
    extrema = [float("inf"), float("-inf")]
    delta = _snapshot_matrix(snapshots, "delta_snapshot_meV", n_nodes)
    for index in range(delta.shape[0]):
        if real.size and imag.size:
            psi = real[index] + 1j * imag[index]
        else:
            psi = delta[index] * MEV_J + 0.0j
        q_abs = phase_gradient_q_abs_m_inv(
            tri,
            psi,
            x_nm=np.asarray(tri.x, dtype=float),
            y_nm=np.asarray(tri.y, dtype=float),
        )
        _update_pair_extrema(extrema, q_abs)
        if index in selected_position:
            q_selected[selected_position[index]] = q_abs
    limits = np.asarray(extrema, dtype=float)
    if not np.all(np.isfinite(limits)):
        limits = np.asarray([0.0, 1.0], dtype=float)
    return xi_m * q_selected, xi_m * limits


def _finite_limits(values: np.ndarray) -> np.ndarray:
    array = np.asarray(values, dtype=float)
    lo = float("inf")
    hi = float("-inf")
    rows = array.shape[0] if array.ndim else 1
    for start in range(0, rows, 64):
        chunk = array[start : start + 64] if array.ndim else array.reshape(1)
        finite = chunk[np.isfinite(chunk)]
        if finite.size:
            lo = min(lo, float(np.min(finite)))
            hi = max(hi, float(np.max(finite)))
    if not np.isfinite(lo) or not np.isfinite(hi):
        return np.asarray([0.0, 1.0], dtype=float)
    return np.asarray([lo, hi], dtype=float)


def _update_pair_extrema(bounds: list[float], values: np.ndarray) -> None:
    finite = np.asarray(values, dtype=float)
    finite = finite[np.isfinite(finite)]
    if finite.size:
        bounds[0] = min(bounds[0], float(np.min(finite)))
        bounds[1] = max(bounds[1], float(np.max(finite)))


def _format_strip_axis(axis: plt.Axes, tri: mtri.Triangulation) -> None:
    axis.set_aspect("equal", adjustable="box")
    axis.set_xlim(float(np.nanmin(tri.x)), float(np.nanmax(tri.x)))
    axis.set_ylim(float(np.nanmin(tri.y)), float(np.nanmax(tri.y)))
    axis.grid(False)
    axis.tick_params(labelsize=6.8, length=2.0, pad=1.0)


def _photon_time_ps(
    history: Mapping[str, Any],
    summary: Mapping[str, Any],
) -> float:
    value = _nested_float(summary, ("photon", "time_s"), default=np.nan)
    if np.isfinite(value):
        return 1.0e12 * value
    time = _history_time_ps(history)
    applied = np.asarray(history.get("photon_applied", []), dtype=bool).reshape(-1)
    count = min(time.size, applied.size)
    indices = np.flatnonzero(applied[:count])
    return float(time[indices[0]]) if indices.size else np.nan


def _impact_coordinates_nm(summary: Mapping[str, Any]) -> np.ndarray:
    return np.asarray(
        [
            1.0e9 * _nested_float(summary, ("photon", "x_m"), default=np.nan),
            1.0e9 * _nested_float(summary, ("photon", "y_m"), default=np.nan),
        ],
        dtype=float,
    )


def _nested_float(
    mapping: Mapping[str, Any],
    path: tuple[str, ...],
    *,
    default: float,
) -> float:
    value: Any = mapping
    for key in path:
        if not isinstance(value, Mapping) or key not in value:
            return float(default)
        value = value[key]
    try:
        return float(value)
    except (TypeError, ValueError):
        return float(default)


def _timing_value(mapping: Mapping[str, Any], key: str) -> str:
    value = mapping.get(key)
    return f"{float(value):.3g} ps" if value is not None else "censored"


def _combined_legend(
    axis: plt.Axes,
    twin: plt.Axes,
    *,
    loc: str = "best",
    bbox_to_anchor: tuple[float, float] | None = None,
) -> None:
    handles_a, labels_a = axis.get_legend_handles_labels()
    handles_b, labels_b = twin.get_legend_handles_labels()
    axis.legend(
        handles_a + handles_b,
        labels_a + labels_b,
        frameon=False,
        ncol=2,
        loc=loc,
        bbox_to_anchor=bbox_to_anchor,
    )


def _prepare_output(path: str | Path) -> Path:
    output = Path(path)
    if output.suffix.lower() != ".pdf":
        output = output.with_suffix(".pdf")
    output.parent.mkdir(parents=True, exist_ok=True)
    return output


__all__ = [
    "make_photon_run_diagnostic_figures",
    "nearest_unique_snapshot_indices",
    "plot_photon_censored_recovery_diagnostics",
    "plot_photon_field_evolution",
    "plot_photon_scalar_evolution",
]
