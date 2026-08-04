"""Thesis figures comparing two completed photon-impact transients."""

from __future__ import annotations

from pathlib import Path
from typing import Any, Mapping, Sequence

import matplotlib

matplotlib.use("Agg", force=True)
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
from matplotlib.colors import LogNorm, Normalize, SymLogNorm
import numpy as np

from pysnspd.plotting.style import THESIS_DPI, THESIS_WIDTH_IN, apply_thesis_style

apply_thesis_style()


def make_photon_position_figures(
    *,
    mesh: Any,
    center_history: Mapping[str, Any],
    center_snapshots: Mapping[str, Any],
    center_summary: Mapping[str, Any],
    edge_history: Mapping[str, Any],
    edge_snapshots: Mapping[str, Any],
    edge_summary: Mapping[str, Any],
    delta0_meV: float,
    xi_m: float,
    requested_times_ps: Sequence[float],
    output_dir: str | Path,
    dpi: int = THESIS_DPI,
    center_timing: Mapping[str, Any] | None = None,
    edge_timing: Mapping[str, Any] | None = None,
    center_snapshot_diagnostics: Mapping[str, Any] | None = None,
    edge_snapshot_diagnostics: Mapping[str, Any] | None = None,
) -> dict[str, Path]:
    """Create matched field, thermodynamic and circuit comparisons."""

    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)
    center_summary_with_timing = dict(center_summary)
    center_summary_with_timing["_timing"] = dict(center_timing or {})
    edge_summary_with_timing = dict(edge_summary)
    edge_summary_with_timing["_timing"] = dict(edge_timing or {})
    runs = (
        (
            "Center",
            center_history,
            center_snapshots,
            center_summary_with_timing,
            dict(center_snapshot_diagnostics or {}),
        ),
        (
            "Edge",
            edge_history,
            edge_snapshots,
            edge_summary_with_timing,
            dict(edge_snapshot_diagnostics or {}),
        ),
    )
    saved = {
        "field_comparison": plot_photon_position_field_comparison(
            mesh,
            runs,
            delta0_meV=delta0_meV,
            xi_m=xi_m,
            requested_times_ps=requested_times_ps,
            output_path=out / "E3_photon_position_field_comparison.pdf",
            dpi=dpi,
        ),
        "circuit_comparison": plot_photon_position_circuit_comparison(
            runs,
            output_path=out / "E3_photon_position_circuit_comparison.pdf",
            dpi=dpi,
        ),
        "power_density_comparison": plot_photon_position_power_density_comparison(
            mesh,
            runs,
            requested_times_ps=requested_times_ps,
            output_path=out / "E3_photon_position_power_density_comparison.pdf",
            dpi=dpi,
        ),
        "energy_heat_capacity_comparison": (
            plot_photon_position_energy_heat_capacity_comparison(
                mesh,
                runs,
                requested_times_ps=requested_times_ps,
                output_path=(
                    out / "E3_photon_position_energy_heat_capacity_comparison.pdf"
                ),
                dpi=dpi,
            )
        ),
    }
    recovery_path = out / "E3_photon_position_censored_recovery_diagnostics.pdf"
    if any(_is_detected_but_unrecovered(run[3].get("_timing", {})) for run in runs):
        saved["censored_recovery_diagnostics"] = (
            plot_photon_position_censored_recovery_diagnostics(
                runs,
                output_path=recovery_path,
                dpi=dpi,
            )
        )
    elif recovery_path.exists():
        recovery_path.unlink()
    return saved


def plot_photon_position_field_comparison(
    mesh: Any,
    runs: Sequence[
        tuple[
            str,
            Mapping[str, Any],
            Mapping[str, Any],
            Mapping[str, Any],
            Mapping[str, Any],
        ]
    ],
    *,
    delta0_meV: float,
    xi_m: float,
    requested_times_ps: Sequence[float],
    output_path: str | Path,
    dpi: int = THESIS_DPI,
) -> Path:
    """Plot matched condensate, potential, momentum and temperature maps."""

    output = _prepare_output(output_path)
    nodes = np.asarray(mesh.nodes, dtype=float)
    x_nm = 1.0e9 * nodes[:, 0]
    y_nm = 1.0e9 * nodes[:, 1]
    tri = mtri.Triangulation(x_nm, y_nm, np.asarray(mesh.triangles, dtype=np.int64))
    requested = [float(value) for value in requested_times_ps]
    if not requested:
        raise ValueError("At least one E3 comparison time is required.")
    if not np.isfinite(delta0_meV) or delta0_meV <= 0.0:
        raise ValueError("A positive Delta_BCS(0) is required for E3 normalization.")
    if not np.isfinite(xi_m) or xi_m <= 0.0:
        raise ValueError("A positive xi is required for E3 momentum normalization.")

    rows: list[dict[str, Any]] = []
    for requested_time in requested:
        for label, _, snapshots, summary, diagnostics in runs:
            stored_time = np.asarray(snapshots.get("snapshot_t_ps", []), dtype=float)
            if stored_time.size == 0:
                raise ValueError(f"{label} run has no photon snapshots.")
            index = int(np.nanargmin(np.abs(stored_time - requested_time)))
            rows.append(
                {
                    "label": label,
                    "requested_time_ps": requested_time,
                    "stored_time_ps": float(stored_time[index]),
                    "delta": np.asarray(snapshots["delta_snapshot_meV"], dtype=float)[index] / delta0_meV,
                    "phi": 1.0e3 * np.asarray(snapshots["phi_snapshot_V"], dtype=float)[index],
                    "qxi": xi_m * _selected_diagnostic_map(
                        diagnostics,
                        key="q_abs_snapshot_m_inv",
                        snapshot_index=index,
                        label=label,
                    ),
                    "Te": np.asarray(snapshots["Te_snapshot_K"], dtype=float)[index],
                    "Tph": np.asarray(snapshots["Tph_snapshot_K"], dtype=float)[index],
                    "impact": _impact_coordinates_nm(summary),
                }
            )

    field_specs = (
        ("delta", r"$|\Delta|/\Delta_{\mathrm{BCS}}(0)$", "viridis", False, 0.0),
        ("phi", r"$\phi$ [mV]", "coolwarm", True, None),
        ("qxi", r"$|\mathbf{q}|\xi$", "magma", False, 0.0),
        ("Te", r"$T_e$ [K]", "inferno", False, None),
        ("Tph", r"$T_{ph}$ [K]", "inferno", False, None),
    )
    limits: dict[str, tuple[float, float]] = {}
    for key, _, _, symmetric, forced_min in field_specs:
        if key == "qxi":
            limits[key] = _combined_diagnostic_limits(
                runs,
                key="q_abs_snapshot_m_inv",
                scale=xi_m,
                symmetric=False,
                forced_min=0.0,
            )
        else:
            stack = _full_snapshot_field(runs, key=key, delta0_meV=delta0_meV)
            limits[key] = _field_limits(
                stack, symmetric=symmetric, forced_min=forced_min
            )

    n_rows = len(rows)
    fig, axes = plt.subplots(
        n_rows,
        len(field_specs),
        figsize=(THESIS_WIDTH_IN, max(4.8, 0.58 * n_rows + 1.0)),
        squeeze=False,
    )
    fig.subplots_adjust(left=0.078, right=0.985, bottom=0.072, top=0.905, wspace=0.08, hspace=0.12)

    mappables = []
    for col, (key, label, cmap, _, _) in enumerate(field_specs):
        vmin, vmax = limits[key]
        for row_index, row in enumerate(rows):
            ax = axes[row_index, col]
            mappable = ax.tripcolor(
                tri,
                row[key],
                shading="gouraud",
                cmap=cmap,
                vmin=vmin,
                vmax=vmax,
                rasterized=True,
            )
            _format_strip_axis(ax, tri)
            impact_x, impact_y = row["impact"]
            if np.isfinite(impact_x) and np.isfinite(impact_y):
                ax.plot(impact_x, impact_y, marker="x", markersize=3.6, markeredgewidth=0.8, color="white")
            if row_index < n_rows - 1:
                ax.tick_params(axis="x", labelbottom=False)
            if col != 0:
                ax.tick_params(axis="y", labelleft=False)
            if col == len(field_specs) - 1:
                ax.text(
                    0.975,
                    0.88,
                    rf"{row['label']}, $t={row['stored_time_ps']:.3g}$ [ps]",
                    transform=ax.transAxes,
                    va="top",
                    ha="right",
                    fontsize=7.1,
                    color="white",
                    bbox={"facecolor": "0.1", "edgecolor": "none", "alpha": 0.65, "pad": 1.0},
                )
        mappables.append((mappable, label))

    fig.supxlabel(r"$x$ [nm]", y=0.018, fontsize=8.5)
    fig.supylabel(r"$y$ [nm]", x=0.016, fontsize=8.5)

    fig.canvas.draw()
    for col, (mappable, label) in enumerate(mappables):
        position = axes[0, col].get_position()
        cax = fig.add_axes([position.x0, position.y1 + 0.010, position.width, 0.011])
        colorbar = fig.colorbar(mappable, cax=cax, orientation="horizontal")
        cax.xaxis.set_ticks_position("top")
        cax.xaxis.set_label_position("top")
        colorbar.set_label(label, labelpad=1.5, fontsize=8.4)
        cax.tick_params(labelsize=6.8, pad=0.8, length=2.0)

    fig.savefig(output, dpi=dpi)
    plt.close(fig)
    return output


def plot_photon_position_circuit_comparison(
    runs: Sequence[
        tuple[
            str,
            Mapping[str, Any],
            Mapping[str, Any],
            Mapping[str, Any],
            Mapping[str, Any],
        ]
    ],
    *,
    output_path: str | Path,
    dpi: int = THESIS_DPI,
) -> Path:
    """Compare circuit, thermal and condensate histories for both positions."""

    output = _prepare_output(output_path)
    fig, axes = plt.subplots(4, 1, figsize=(THESIS_WIDTH_IN, 6.3), sharex=True)
    fig.subplots_adjust(left=0.115, right=0.975, bottom=0.080, top=0.985, hspace=0.18)
    colors = {"Center": "tab:blue", "Edge": "tab:red"}
    photon_times = []

    timing_rows: list[str] = []
    for label, history, _, summary, _ in runs:
        time = np.asarray(history.get("t_ps", []), dtype=float)
        if time.size == 0:
            raise ValueError(f"{label} run has no photon history.")
        color = colors.get(label, None)
        current_s = 1.0e6 * np.asarray(history.get("I_s_A", []), dtype=float)
        current_rf = 1.0e6 * np.asarray(history.get("I_rf_A", []), dtype=float)
        delta_current_s = current_s - float(current_s[0])
        axes[0].plot(time, delta_current_s, color=color, label=rf"{label}: $\Delta I_s$")
        axes[0].plot(time, current_rf, color=color, linestyle="--", label=rf"{label}: $I_{{\mathrm{{RF}}}}$")

        voltage_tdgl = 1.0e3 * np.asarray(history.get("V_tdgl_center_V", []), dtype=float)
        voltage_out = 1.0e3 * np.asarray(history.get("V_out_V", []), dtype=float)
        axes[1].plot(time, voltage_tdgl, color=color, label=rf"{label}: $V_{{\mathrm{{TDGL}}}}$")
        axes[1].plot(time, voltage_out, color=color, linestyle="--", label=rf"{label}: $V_{{\mathrm{{out}}}}$")

        axes[2].plot(time, history.get("max_Te_K"), color=color, label=rf"{label}: max $T_e$")
        axes[2].plot(time, history.get("max_Tph_K"), color=color, linestyle="--", label=rf"{label}: max $T_{{ph}}$")
        axes[3].plot(time, history.get("mean_delta_over_delta0"), color=color, label=label)
        photon_times.append(_photon_time_ps(history))
        timing = dict(summary.get("_timing", {}))
        latency = dict(timing.get("latency", {}))
        recovery = dict(dict(timing.get("recovery", {})).get("selected", {}))
        crossing = latency.get("crossing_time_ps")
        recovery_entry = recovery.get("entry_time_ps")
        if crossing is not None and np.isfinite(float(crossing)):
            axes[1].axvline(float(crossing), color=color, linestyle=":", linewidth=0.8)
        if recovery_entry is not None and np.isfinite(float(recovery_entry)):
            axes[1].axvline(float(recovery_entry), color=color, linestyle="-.", linewidth=0.8)
        latency_text = (
            f"{float(latency['t_lat_ps']):.3g} ps"
            if latency.get("t_lat_ps") is not None
            else "censored"
        )
        recovery_text = (
            f"{float(recovery['t_rec_ps']):.3g} ps"
            if recovery.get("t_rec_ps") is not None
            else "censored"
        )
        timing_rows.append(
            f"{label}: t_lat={latency_text}; "
            f"t_rec[{recovery.get('mode', 'electrical')}]={recovery_text}"
        )

    axes[0].set_ylabel(r"Current [$\mu$A]")
    axes[1].set_ylabel("Voltage [mV]")
    axes[2].set_ylabel("Temperature [K]")
    axes[3].set_ylabel(r"Mean $|\Delta|/\Delta_{\mathrm{BCS}}(0)$")
    axes[3].set_xlabel(r"$t$ [ps]")

    photon_time = float(np.nanmedian(photon_times)) if np.isfinite(photon_times).any() else np.nan
    for index, ax in enumerate(axes):
        if np.isfinite(photon_time):
            ax.axvline(
                photon_time,
                color="0.25",
                linestyle=":",
                linewidth=0.9,
                label="Photon arrival" if index == 0 else None,
            )
        ax.grid(True)
        ax.set_xlim(left=0.0)
        ax.legend(frameon=False, ncol=3 if index == 0 else 2, loc="best", fontsize=8.0)
    axes[3].text(
        0.01,
        0.04,
        "\n".join(timing_rows),
        transform=axes[3].transAxes,
        ha="left",
        va="bottom",
        fontsize=7.0,
        bbox={"facecolor": "white", "edgecolor": "0.75", "alpha": 0.82, "pad": 1.5},
    )

    fig.savefig(output, dpi=dpi)
    plt.close(fig)
    return output


def plot_photon_position_power_density_comparison(
    mesh: Any,
    runs: Sequence[tuple[str, Mapping[str, Any], Mapping[str, Any], Mapping[str, Any], Mapping[str, Any]]],
    *,
    requested_times_ps: Sequence[float],
    output_path: str | Path,
    dpi: int = THESIS_DPI,
) -> Path:
    """Compare runtime-consistent power-density maps using shared scales."""

    channels = (
        ("joule_snapshot_W_m3", r"$P_J$ [W m$^{-3}$]", "magma", True),
        (
            "P_total_snapshot_W_m3",
            r"$P_{e\mathrm{-}ph}=P_S+P_R$ [W m$^{-3}$]",
            "coolwarm",
            False,
        ),
        ("P_diff_snapshot_W_m3", r"$P_{\mathrm{diff}}$ [W m$^{-3}$]", "PuOr_r", False),
        ("P_esc_snapshot_W_m3", r"$P_{\mathrm{esc}}$ [W m$^{-3}$]", "cividis", True),
    )
    return _plot_photon_position_diagnostic_atlas(
        mesh=mesh,
        runs=runs,
        requested_times_ps=requested_times_ps,
        channels=channels,
        output_path=output_path,
        dpi=dpi,
    )


def plot_photon_position_energy_heat_capacity_comparison(
    mesh: Any,
    runs: Sequence[tuple[str, Mapping[str, Any], Mapping[str, Any], Mapping[str, Any], Mapping[str, Any]]],
    *,
    requested_times_ps: Sequence[float],
    output_path: str | Path,
    dpi: int = THESIS_DPI,
) -> Path:
    """Compare electronic/phononic energy and heat-capacity maps."""

    channels = (
        ("u_e_snapshot_J_m3", r"$u_e$ [J m$^{-3}$]", "coolwarm", False),
        ("u_ph_snapshot_J_m3", r"$u_{ph}$ [J m$^{-3}$]", "viridis", True),
        (
            "C_e_snapshot_J_m3_K",
            r"$C_e$ [J m$^{-3}$ K$^{-1}$]",
            "magma",
            True,
        ),
        (
            "C_ph_snapshot_J_m3_K",
            r"$C_{ph}$ [J m$^{-3}$ K$^{-1}$]",
            "cividis",
            True,
        ),
    )
    return _plot_photon_position_diagnostic_atlas(
        mesh=mesh,
        runs=runs,
        requested_times_ps=requested_times_ps,
        channels=channels,
        output_path=output_path,
        dpi=dpi,
    )


def _plot_photon_position_diagnostic_atlas(
    *,
    mesh: Any,
    runs: Sequence[tuple[str, Mapping[str, Any], Mapping[str, Any], Mapping[str, Any], Mapping[str, Any]]],
    requested_times_ps: Sequence[float],
    channels: Sequence[tuple[str, str, str, bool]],
    output_path: str | Path,
    dpi: int,
) -> Path:
    output = _prepare_output(output_path)
    nodes = np.asarray(mesh.nodes, dtype=float)
    tri = mtri.Triangulation(
        1.0e9 * nodes[:, 0],
        1.0e9 * nodes[:, 1],
        np.asarray(mesh.triangles, dtype=np.int64),
    )
    rows = _diagnostic_rows(runs, requested_times_ps=requested_times_ps)
    if not rows:
        raise ValueError("No matched photon diagnostics are available.")
    norms = {
        key: _diagnostic_norm(runs, key=key, positive=positive)
        for key, _, _, positive in channels
    }
    fig, axes = plt.subplots(
        len(rows),
        len(channels),
        figsize=(THESIS_WIDTH_IN, max(4.8, 0.58 * len(rows) + 1.0)),
        squeeze=False,
    )
    fig.subplots_adjust(
        left=0.078,
        right=0.985,
        bottom=0.072,
        top=0.905,
        wspace=0.08,
        hspace=0.12,
    )
    mappables: list[tuple[Any, str]] = []
    for col, (key, label, cmap, _) in enumerate(channels):
        for row_index, row in enumerate(rows):
            axis = axes[row_index, col]
            mappable = axis.tripcolor(
                tri,
                row["diagnostics"][key],
                shading="gouraud",
                cmap=cmap,
                norm=norms[key],
                rasterized=True,
            )
            _format_strip_axis(axis, tri)
            impact_x, impact_y = row["impact"]
            if np.isfinite(impact_x) and np.isfinite(impact_y):
                axis.plot(
                    impact_x,
                    impact_y,
                    marker="x",
                    markersize=3.6,
                    markeredgewidth=0.8,
                    color="white",
                )
            if row_index < len(rows) - 1:
                axis.tick_params(axis="x", labelbottom=False)
            if col != 0:
                axis.tick_params(axis="y", labelleft=False)
            if col == len(channels) - 1:
                axis.text(
                    0.975,
                    0.88,
                    rf"{row['label']}, $t={row['stored_time_ps']:.3g}$ [ps]",
                    transform=axis.transAxes,
                    va="top",
                    ha="right",
                    fontsize=7.1,
                    color="white",
                    bbox={"facecolor": "0.1", "edgecolor": "none", "alpha": 0.65, "pad": 1.0},
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
        if isinstance(mappable.norm, SymLogNorm):
            symmetric_max = max(
                abs(float(mappable.norm.vmin)),
                abs(float(mappable.norm.vmax)),
            )
            colorbar.set_ticks([-symmetric_max, 0.0, symmetric_max])
            colorbar.set_ticklabels(
                [f"-{symmetric_max:.0e}", "0", f"{symmetric_max:.0e}"]
            )
            tick_labels = colorbar.ax.get_xticklabels()
            if len(tick_labels) >= 2:
                tick_labels[0].set_ha("left")
                tick_labels[-1].set_ha("right")
        color_axis.xaxis.set_ticks_position("top")
        color_axis.xaxis.set_label_position("top")
        colorbar.set_label(label, labelpad=1.5, fontsize=8.4)
        color_axis.tick_params(labelsize=6.8, pad=0.8, length=2.0)
    fig.savefig(output, dpi=dpi)
    plt.close(fig)
    return output


def plot_photon_position_censored_recovery_diagnostics(
    runs: Sequence[tuple[str, Mapping[str, Any], Mapping[str, Any], Mapping[str, Any], Mapping[str, Any]]],
    *,
    output_path: str | Path,
    dpi: int = THESIS_DPI,
) -> Path:
    """Compare recovery margins and slow circuit modes in four panels."""

    output = _prepare_output(output_path)
    payloads = [payload for run in runs if (payload := _recovery_payload(run)) is not None]
    if not payloads:
        raise ValueError("Recovery baselines and tolerances are unavailable for both runs.")
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
    colors = {"Center": "tab:blue", "Edge": "tab:red"}
    styles = ("-", "--", ":")
    fig, axes = plt.subplots(2, 2, figsize=(THESIS_WIDTH_IN, 5.35))
    fig.subplots_adjust(left=0.115, right=0.965, bottom=0.105, top=0.900, wspace=0.40, hspace=0.32)
    for axis, keys, title in (
        (axes[0, 0], current_keys, "Current recovery margins"),
        (axes[0, 1], voltage_keys, "Voltage recovery margins"),
    ):
        for payload in payloads:
            for style, key in zip(styles, keys):
                axis.plot(
                    payload["relative_time_ps"],
                    np.maximum(payload["ratios"][key], 1.0e-6),
                    color=colors.get(payload["label"]),
                    linestyle=style,
                    linewidth=1.0,
                    label=f"{payload['label']}: {labels[key]}",
                )
        axis.axhline(1.0, color="0.2", linestyle="--", linewidth=0.9, label="Tolerance")
        axis.set_yscale("log")
        axis.set_xlabel(r"$t-t_\gamma$ [ps]")
        axis.set_ylabel("Residual / tolerance")
        axis.set_title(title)
        axis.legend(frameon=False, fontsize=6.2, ncol=2, loc="upper right")
        axis.grid(True)
        axis.set_xlim(left=0.0)

    axis = axes[1, 0]
    final_keys = current_keys + voltage_keys
    y = np.arange(len(final_keys), dtype=float)
    width = 0.34
    offsets = np.linspace(-0.5 * width, 0.5 * width, len(payloads)) if len(payloads) > 1 else np.zeros(1)
    for offset, payload in zip(offsets, payloads):
        final = np.asarray([payload["ratios"][key][-1] for key in final_keys], dtype=float)
        axis.barh(
            y + offset,
            np.maximum(final, 1.0e-4),
            height=width,
            color=colors.get(payload["label"]),
            alpha=0.82,
            label=payload["label"],
        )
    axis.axvline(1.0, color="0.2", linestyle="--", linewidth=0.9)
    axis.set_xscale("log")
    axis.set_yticks(y)
    axis.set_yticklabels([labels[key] for key in final_keys])
    axis.invert_yaxis()
    axis.set_xlabel("Final residual / tolerance")
    axis.set_title("Distance from electrical recovery")
    axis.legend(frameon=False, fontsize=7.0)
    axis.grid(True, axis="x")

    axis = axes[1, 1]
    mode_rows = _comparison_circuit_timescales(runs, payloads)
    if mode_rows:
        row = np.arange(len(mode_rows))
        mode_values = [item[1] for item in mode_rows]
        axis.barh(row, mode_values, color=[item[2] for item in mode_rows], alpha=0.88)
        axis.set_xscale("log")
        axis.set_yticks(row)
        axis.set_yticklabels([item[0] for item in mode_rows], fontsize=7.0)
        axis.invert_yaxis()
        axis.set_xlabel("Timescale [ps]")
        axis.grid(True, axis="x")
        largest_mode = max(mode_values)
        for index, (_, value, _) in enumerate(mode_rows):
            inside = value >= 0.8 * largest_mode
            axis.text(
                value / 1.08 if inside else value * 1.08,
                index,
                f"{value:.3g}",
                va="center",
                ha="right" if inside else "left",
                color="white" if inside else "black",
                fontsize=6.6,
            )
        axis.set_xlim(right=largest_mode * 1.25)
    else:
        axis.text(0.5, 0.5, "Circuit parameters unavailable", transform=axis.transAxes, ha="center", va="center")
        axis.set_xticks([])
        axis.set_yticks([])
    axis.set_title("Circuit modes vs available windows")
    censored = [payload["label"] for payload in payloads if payload["censored"]]
    fig.suptitle(
        "Detected transient; recovery censored: " + ", ".join(censored),
        y=0.975,
        fontsize=10.0,
    )
    fig.savefig(output, dpi=dpi)
    plt.close(fig)
    return output


def _prepare_output(path: str | Path) -> Path:
    output = Path(path)
    if output.suffix.lower() != ".pdf":
        output = output.with_suffix(".pdf")
    output.parent.mkdir(parents=True, exist_ok=True)
    return output


def _impact_coordinates_nm(summary: Mapping[str, Any]) -> tuple[float, float]:
    photon = summary.get("photon", {}) if isinstance(summary, Mapping) else {}
    if not isinstance(photon, Mapping):
        return np.nan, np.nan
    return 1.0e9 * float(photon.get("x_m", np.nan)), 1.0e9 * float(photon.get("y_m", np.nan))


def _photon_time_ps(history: Mapping[str, Any]) -> float:
    time = np.asarray(history.get("t_ps", []), dtype=float)
    applied = np.asarray(history.get("photon_applied", []), dtype=bool)
    if time.size == 0 or applied.size == 0:
        return np.nan
    if applied.size != time.size:
        applied = np.resize(applied, time.size)
    indices = np.flatnonzero(applied)
    return float(time[indices[0]]) if indices.size else np.nan


def _field_limits(
    values: np.ndarray,
    *,
    symmetric: bool,
    forced_min: float | None,
) -> tuple[float, float]:
    finite = np.asarray(values, dtype=float)
    finite = finite[np.isfinite(finite)]
    if finite.size == 0:
        return (-1.0, 1.0) if symmetric else (0.0, 1.0)
    if symmetric:
        limit = max(float(np.nanmax(np.abs(finite))), 1.0e-30)
        return -limit, limit
    lower = float(forced_min) if forced_min is not None else float(np.nanmin(finite))
    upper = float(np.nanmax(finite))
    if forced_min == 0.0:
        upper = max(1.0, upper)
    if not np.isfinite(upper) or upper <= lower:
        upper = lower + 1.0
    return lower, upper


def _full_snapshot_field(
    runs: Sequence[
        tuple[
            str,
            Mapping[str, Any],
            Mapping[str, Any],
            Mapping[str, Any],
            Mapping[str, Any],
        ]
    ],
    *,
    key: str,
    delta0_meV: float,
) -> np.ndarray:
    source = {
        "delta": ("delta_snapshot_meV", 1.0 / float(delta0_meV)),
        "phi": ("phi_snapshot_V", 1.0e3),
        "Te": ("Te_snapshot_K", 1.0),
        "Tph": ("Tph_snapshot_K", 1.0),
    }
    if key not in source:
        raise KeyError(f"Unknown photon field: {key}")
    source_key, scale = source[key]
    parts = []
    for label, _, snapshots, _, _ in runs:
        values = np.asarray(snapshots.get(source_key, []), dtype=float)
        if values.size == 0:
            raise ValueError(f"{label} run has no stored field {source_key}.")
        parts.append((scale * values).reshape(-1))
    return np.concatenate(parts)


def _selected_diagnostic_map(
    diagnostics: Mapping[str, Any],
    *,
    key: str,
    snapshot_index: int,
    label: str,
) -> np.ndarray:
    selected = np.asarray(diagnostics.get("selected_indices", []), dtype=np.int64).reshape(-1)
    matches = np.flatnonzero(selected == int(snapshot_index))
    values = np.asarray(diagnostics.get(key, []), dtype=float)
    if matches.size != 1 or values.ndim != 2 or values.shape[0] != selected.size:
        raise ValueError(f"{label} diagnostics do not contain {key} for snapshot {snapshot_index}.")
    return values[int(matches[0])]


def _combined_diagnostic_limits(
    runs: Sequence[tuple[str, Mapping[str, Any], Mapping[str, Any], Mapping[str, Any], Mapping[str, Any]]],
    *,
    key: str,
    scale: float,
    symmetric: bool,
    forced_min: float | None,
) -> tuple[float, float]:
    pairs = []
    for label, _, _, _, diagnostics in runs:
        limits = np.asarray(
            dict(diagnostics.get("snapshot_global_limits", {})).get(key, []),
            dtype=float,
        ).reshape(-1)
        if limits.size != 2 or not np.all(np.isfinite(limits)):
            raise ValueError(f"{label} diagnostics lack global limits for {key}.")
        pairs.append(float(scale) * limits)
    return _field_limits(
        np.concatenate(pairs),
        symmetric=symmetric,
        forced_min=forced_min,
    )


def _diagnostic_rows(
    runs: Sequence[tuple[str, Mapping[str, Any], Mapping[str, Any], Mapping[str, Any], Mapping[str, Any]]],
    *,
    requested_times_ps: Sequence[float],
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for requested in requested_times_ps:
        for label, _, _, summary, diagnostics in runs:
            times = np.asarray(diagnostics.get("snapshot_t_ps", []), dtype=float).reshape(-1)
            if times.size == 0:
                raise ValueError(f"{label} run has no selected photon diagnostics.")
            index = int(np.nanargmin(np.abs(times - float(requested))))
            rows.append(
                {
                    "label": label,
                    "stored_time_ps": float(times[index]),
                    "impact": _impact_coordinates_nm(summary),
                    "diagnostics": {
                        key: np.asarray(value, dtype=float)[index]
                        for key, value in diagnostics.items()
                        if isinstance(value, np.ndarray)
                        and np.asarray(value).ndim == 2
                        and np.asarray(value).shape[0] == times.size
                    },
                }
            )
    return rows


def _diagnostic_norm(
    runs: Sequence[tuple[str, Mapping[str, Any], Mapping[str, Any], Mapping[str, Any], Mapping[str, Any]]],
    *,
    key: str,
    positive: bool,
):
    bounds = []
    selected_values = []
    for label, _, _, _, diagnostics in runs:
        limits = np.asarray(
            dict(diagnostics.get("snapshot_global_limits", {})).get(key, []),
            dtype=float,
        ).reshape(-1)
        if limits.size != 2 or not np.all(np.isfinite(limits)):
            raise ValueError(f"{label} diagnostics lack global limits for {key}.")
        bounds.append(limits)
        selected_values.append(np.asarray(diagnostics.get(key, []), dtype=float).reshape(-1))
    combined = np.concatenate(bounds)
    if positive:
        vmax = max(float(np.nanmax(combined)), 1.0e-300)
        values = np.concatenate(selected_values)
        finite_positive = values[np.isfinite(values) & (values > 0.0)]
        if finite_positive.size and vmax > 0.0:
            vmin = max(float(np.nanmin(finite_positive)), vmax * 1.0e-8)
            if vmax > vmin:
                return LogNorm(vmin=vmin, vmax=vmax)
        return Normalize(vmin=0.0, vmax=max(vmax, 1.0))
    vmax = max(float(np.nanmax(np.abs(combined))), 1.0e-300)
    return SymLogNorm(
        linthresh=max(vmax * 1.0e-6, 1.0e-300),
        vmin=-vmax,
        vmax=vmax,
        base=10.0,
    )


def _is_detected_but_unrecovered(timing: Mapping[str, Any]) -> bool:
    latency = dict(timing.get("latency", {}))
    recovery = dict(dict(timing.get("recovery", {})).get("selected", {}))
    detected = bool(latency.get("detected", latency.get("t_lat_ps") is not None))
    recovered = bool(recovery.get("recovered", recovery.get("t_rec_ps") is not None))
    return detected and not recovered


def _recovery_payload(
    run: tuple[str, Mapping[str, Any], Mapping[str, Any], Mapping[str, Any], Mapping[str, Any]],
) -> dict[str, Any] | None:
    label, history, _, summary, _ = run
    timing = dict(summary.get("_timing", {}))
    time_ps = np.asarray(history.get("t_ps", []), dtype=float).reshape(-1)
    event_ps = _photon_time_ps(history)
    if time_ps.size < 2 or not np.isfinite(event_ps):
        return None
    post = np.isfinite(time_ps) & (time_ps >= event_ps)
    baseline = dict(timing.get("baseline", {})).get("values", {})
    recovery = dict(dict(timing.get("recovery", {})).get("selected", {}))
    tolerances = dict(recovery.get("absolute_tolerances", {}))
    keys = ("I_b_A", "I_s_A", "I_rf_A", "V_out_V", "v_c_V", "V_tdgl_center_V")
    if not isinstance(baseline, Mapping) or not all(key in tolerances for key in keys):
        return None
    ratios: dict[str, np.ndarray] = {}
    for key in keys:
        values = np.asarray(history.get(key, []), dtype=float).reshape(-1)
        if values.size != time_ps.size:
            values = np.resize(values, time_ps.size)
        tolerance = max(float(tolerances[key]), 1.0e-300)
        ratios[key] = np.abs(values[post] - float(baseline.get(key, np.nan))) / tolerance
    if not np.any(post) or any(values.size == 0 for values in ratios.values()):
        return None
    return {
        "label": label,
        "relative_time_ps": time_ps[post] - event_ps,
        "ratios": ratios,
        "timing": timing,
        "censored": _is_detected_but_unrecovered(timing),
    }


def _comparison_circuit_timescales(
    runs: Sequence[tuple[str, Mapping[str, Any], Mapping[str, Any], Mapping[str, Any], Mapping[str, Any]]],
    payloads: Sequence[Mapping[str, Any]],
) -> list[tuple[str, float, str]]:
    params = dict(dict(runs[0][3].get("circuit", {})).get("params", {}))
    try:
        Rl = float(params["R_load_ohm"])
        Rb = float(params["R_bias_ohm"])
        Lb = float(params["L_bias_H"])
        Lk = float(params["L_k_H"])
        capacitance = float(params["C_couple_F"])
    except (KeyError, TypeError, ValueError):
        return []
    if not all(np.isfinite(value) and value > 0.0 for value in (Rl, Rb, Lb, Lk, capacitance)):
        return []
    matrix = np.asarray(
        [
            [-(Rb + Rl) / Lb, Rl / Lb, -1.0 / Lb],
            [Rl / Lk, -Rl / Lk, 1.0 / Lk],
            [1.0 / capacitance, -1.0 / capacitance, 0.0],
        ],
        dtype=float,
    )
    eigenvalues = np.linalg.eigvals(matrix)
    rows: list[tuple[str, float, str]] = []
    real_modes = [value for value in eigenvalues if abs(float(np.imag(value))) <= 1.0e-8 * max(abs(value), 1.0)]
    complex_modes = [value for value in eigenvalues if float(np.imag(value)) > 0.0]
    for index, value in enumerate(sorted(real_modes, key=lambda item: abs(float(np.real(item))))):
        if float(np.real(value)) < 0.0:
            rows.append((f"Mode {index + 1}: decay", -1.0e12 / float(np.real(value)), "tab:blue"))
    for index, value in enumerate(complex_modes):
        if float(np.real(value)) < 0.0:
            rows.append((f"Osc. {index + 1}: decay", -1.0e12 / float(np.real(value)), "tab:purple"))
        rows.append((f"Osc. {index + 1}: period", 2.0 * np.pi * 1.0e12 / abs(float(np.imag(value))), "tab:orange"))
    hold_values = []
    for payload in payloads:
        hold_s = float(dict(payload["timing"].get("recovery_criteria", {})).get("hold_s", np.nan))
        if np.isfinite(hold_s) and hold_s > 0.0:
            hold_values.append(1.0e12 * hold_s)
    if hold_values:
        rows.append(("Recovery hold", float(np.nanmedian(hold_values)), "tab:green"))
    window_colors = {"Center": "tab:blue", "Edge": "tab:red"}
    for payload in payloads:
        relative = np.asarray(payload["relative_time_ps"], dtype=float)
        if relative.size and np.isfinite(relative[-1]) and relative[-1] > 0.0:
            rows.append((f"{payload['label']} window", float(relative[-1]), window_colors.get(str(payload["label"]), "0.35")))
    return rows


def _format_strip_axis(ax: plt.Axes, tri: mtri.Triangulation) -> None:
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlim(float(np.nanmin(tri.x)), float(np.nanmax(tri.x)))
    ax.set_ylim(float(np.nanmin(tri.y)), float(np.nanmax(tri.y)))
    ax.grid(False)
    ax.tick_params(labelsize=6.8, length=2.0, pad=1.0)


__all__ = [
    "make_photon_position_figures",
    "plot_photon_position_censored_recovery_diagnostics",
    "plot_photon_position_circuit_comparison",
    "plot_photon_position_energy_heat_capacity_comparison",
    "plot_photon_position_field_comparison",
    "plot_photon_position_power_density_comparison",
]
