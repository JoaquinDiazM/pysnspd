"""Thesis figures comparing two completed photon-impact transients."""

from __future__ import annotations

from pathlib import Path
from typing import Any, Mapping, Sequence

import matplotlib

matplotlib.use("Agg", force=True)
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
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
    requested_times_ps: Sequence[float],
    output_dir: str | Path,
    dpi: int = THESIS_DPI,
) -> dict[str, Path]:
    """Create matched field and circuit-response comparisons."""

    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)
    runs = (
        ("Center", center_history, center_snapshots, center_summary),
        ("Edge", edge_history, edge_snapshots, edge_summary),
    )
    return {
        "field_comparison": plot_photon_position_field_comparison(
            mesh,
            runs,
            delta0_meV=delta0_meV,
            requested_times_ps=requested_times_ps,
            output_path=out / "E3_photon_position_field_comparison.pdf",
            dpi=dpi,
        ),
        "circuit_comparison": plot_photon_position_circuit_comparison(
            runs,
            output_path=out / "E3_photon_position_circuit_comparison.pdf",
            dpi=dpi,
        ),
    }


def plot_photon_position_field_comparison(
    mesh: Any,
    runs: Sequence[tuple[str, Mapping[str, Any], Mapping[str, Any], Mapping[str, Any]]],
    *,
    delta0_meV: float,
    requested_times_ps: Sequence[float],
    output_path: str | Path,
    dpi: int = THESIS_DPI,
) -> Path:
    """Plot matched condensate, potential and temperature maps."""

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

    rows: list[dict[str, Any]] = []
    for requested_time in requested:
        for label, _, snapshots, summary in runs:
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
                    "Te": np.asarray(snapshots["Te_snapshot_K"], dtype=float)[index],
                    "Tph": np.asarray(snapshots["Tph_snapshot_K"], dtype=float)[index],
                    "impact": _impact_coordinates_nm(summary),
                }
            )

    field_specs = (
        ("delta", r"$|\Delta|/\Delta_{\mathrm{BCS}}(0)$", "viridis", False, 0.0),
        ("phi", r"$\phi$ [mV]", "coolwarm", True, None),
        ("Te", r"$T_e$ [K]", "inferno", False, None),
        ("Tph", r"$T_{ph}$ [K]", "inferno", False, None),
    )
    limits: dict[str, tuple[float, float]] = {}
    for key, _, _, symmetric, forced_min in field_specs:
        stack = _full_snapshot_field(runs, key=key, delta0_meV=delta0_meV)
        limits[key] = _field_limits(stack, symmetric=symmetric, forced_min=forced_min)
    Tc_K = delta0_meV / (1.764 * 8.617333262e-2)
    for key in ("Te", "Tph"):
        lower, upper = limits[key]
        limits[key] = lower, max(upper, Tc_K)

    n_rows = len(rows)
    fig, axes = plt.subplots(
        n_rows,
        len(field_specs),
        figsize=(THESIS_WIDTH_IN, max(4.8, 0.58 * n_rows + 1.0)),
        squeeze=False,
    )
    fig.subplots_adjust(left=0.066, right=0.985, bottom=0.055, top=0.905, wspace=0.08, hspace=0.12)

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
            else:
                ax.set_xlabel(r"$x$ [nm]", labelpad=1.0)
            if col == 0:
                ax.set_ylabel(r"$y$ [nm]", labelpad=1.0)
            else:
                ax.tick_params(axis="y", labelleft=False)
            if col == len(field_specs) - 1:
                ax.text(
                    0.975,
                    0.88,
                    rf"{row['label']}, {row['stored_time_ps']:.3g} ps",
                    transform=ax.transAxes,
                    va="top",
                    ha="right",
                    fontsize=7.1,
                    color="white",
                    bbox={"facecolor": "0.1", "edgecolor": "none", "alpha": 0.65, "pad": 1.0},
                )
        mappables.append((mappable, label))

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
    runs: Sequence[tuple[str, Mapping[str, Any], Mapping[str, Any], Mapping[str, Any]]],
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

    for label, history, _, _ in runs:
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
        limit = max(float(np.nanpercentile(np.abs(finite), 99.8)), 1.0e-30)
        return -limit, limit
    lower = float(forced_min) if forced_min is not None else float(np.nanpercentile(finite, 0.1))
    upper = float(np.nanpercentile(finite, 99.95))
    if forced_min == 0.0:
        upper = max(1.0, upper)
    if not np.isfinite(upper) or upper <= lower:
        upper = lower + 1.0
    return lower, upper


def _full_snapshot_field(
    runs: Sequence[tuple[str, Mapping[str, Any], Mapping[str, Any], Mapping[str, Any]]],
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
    for label, _, snapshots, _ in runs:
        values = np.asarray(snapshots.get(source_key, []), dtype=float)
        if values.size == 0:
            raise ValueError(f"{label} run has no stored field {source_key}.")
        parts.append((scale * values).reshape(-1))
    return np.concatenate(parts)


def _format_strip_axis(ax: plt.Axes, tri: mtri.Triangulation) -> None:
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlim(float(np.nanmin(tri.x)), float(np.nanmax(tri.x)))
    ax.set_ylim(float(np.nanmin(tri.y)), float(np.nanmax(tri.y)))
    ax.grid(False)
    ax.tick_params(labelsize=6.8, length=2.0, pad=1.0)


__all__ = [
    "make_photon_position_figures",
    "plot_photon_position_circuit_comparison",
    "plot_photon_position_field_comparison",
]
