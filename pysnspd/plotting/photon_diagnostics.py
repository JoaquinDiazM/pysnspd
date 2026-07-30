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
) -> dict[str, Path]:
    """Create scalar-evolution and selected-field figures for one photon run."""

    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)
    return {
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
            output_path=out / "E3_photon_field_evolution.pdf",
            dpi=dpi,
        ),
    }


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
        pair_axis.set_ylabel("Pair breaking [-]", color="tab:red")
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
            _combined_legend(axis, pair_axis)
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

    delta = _snapshot_matrix(snapshots, "delta_snapshot_meV", n_nodes)[indices] / delta0_meV
    phi = 1.0e3 * _snapshot_matrix(snapshots, "phi_snapshot_V", n_nodes)[indices]
    Te = _snapshot_matrix(snapshots, "Te_snapshot_K", n_nodes)[indices]
    Tph = _snapshot_matrix(snapshots, "Tph_snapshot_K", n_nodes)[indices]
    real = _snapshot_matrix(snapshots, "psi_real_snapshot_J", n_nodes, required=False)
    imag = _snapshot_matrix(snapshots, "psi_imag_snapshot_J", n_nodes, required=False)

    qxi_rows = []
    for index in indices:
        if real.size and imag.size:
            psi = real[index] + 1j * imag[index]
        else:
            psi = (
                _snapshot_matrix(snapshots, "delta_snapshot_meV", n_nodes)[index]
                * MEV_J
                + 0.0j
            )
        qxi_rows.append(
            xi_m
            * phase_gradient_q_abs_m_inv(
                tri,
                psi,
                x_nm=np.asarray(tri.x, dtype=float),
                y_nm=np.asarray(tri.y, dtype=float),
            )
        )
    qxi = np.asarray(qxi_rows, dtype=float)

    fields = (
        (delta, r"$|\Delta|/\Delta_{\mathrm{BCS}}(0)$", "viridis", False, 0.0, 1.0),
        (phi, r"$\phi$ [mV]", "coolwarm", True, None, None),
        (qxi, r"$|\mathbf{q}|\xi$", "magma", False, 0.0, None),
        (Te, r"$T_e$ [K]", "inferno", False, None, None),
        (Tph, r"$T_{ph}$ [K]", "inferno", False, None, None),
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
    for col, (values, label, cmap, symmetric, forced_min, forced_max) in enumerate(fields):
        vmin, vmax = _field_limits(
            values,
            symmetric=symmetric,
            forced_min=forced_min,
            forced_max=forced_max,
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
) -> tuple[float, float]:
    finite = np.asarray(values, dtype=float)
    finite = finite[np.isfinite(finite)]
    if finite.size == 0:
        return (-1.0, 1.0) if symmetric else (0.0, 1.0)
    if symmetric:
        limit = max(float(np.nanpercentile(np.abs(finite), 99.8)), 1.0e-30)
        return -limit, limit
    lower = (
        float(forced_min)
        if forced_min is not None
        else float(np.nanpercentile(finite, 0.1))
    )
    upper = (
        float(forced_max)
        if forced_max is not None
        else float(np.nanpercentile(finite, 99.9))
    )
    if not np.isfinite(upper) or upper <= lower:
        upper = lower + 1.0
    return lower, upper


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


def _combined_legend(axis: plt.Axes, twin: plt.Axes) -> None:
    handles_a, labels_a = axis.get_legend_handles_labels()
    handles_b, labels_b = twin.get_legend_handles_labels()
    axis.legend(
        handles_a + handles_b,
        labels_a + labels_b,
        frameon=False,
        ncol=2,
        loc="best",
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
    "plot_photon_field_evolution",
    "plot_photon_scalar_evolution",
]
