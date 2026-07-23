"""Current-sweep plotting helpers for Z-series multi-run analysis.

Notes
-----
The raw stationary SS runs currently save only the final-state fields, not a
long dense voltage time series suitable for a temporal average over the PSL
oscillation cycle. Because of that, the IV figure keeps the raw endpoint
samples as points and overlays a monotone best-fit curve instead of connecting
neighboring points directly. The monotone fit is computed with isotonic
regression (nondecreasing least-squares fit), which is a pragmatic way to
represent the expected macroscopic IV trend while acknowledging the residual
phase-of-oscillation ambiguity of the saved endpoint voltage.
"""

from __future__ import annotations

import csv
import re
from pathlib import Path
from typing import Any, Mapping, Sequence

import matplotlib
matplotlib.use("Agg", force=True)
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize
import numpy as np
import yaml

from pysnspd.analysis.ss_run import build_ss_plot_dataset, load_ss_run
MEV_J = 1.602176634e-22

def _add_delta_insets(ax: plt.Axes, delta_insets: Sequence[Mapping[str, Any]]) -> None:
    """Place four |Delta| colormap insets inside the main IV axes in a 2x2 grid, with a shared colorbar."""
    positions = [
        [0.055, 0.545, 0.15, 0.255],
        [0.235, 0.545, 0.15, 0.255],
        [0.055, 0.205, 0.15, 0.255],
        [0.235, 0.205, 0.15, 0.255],
    ]
    delta_fields_meV = [_extract_delta_field_meV(item.get("dataset", {})) for item in delta_insets]
    finite_maxima = [float(np.nanmax(field)) for field in delta_fields_meV if field.size and np.any(np.isfinite(field))]
    vmax_meV = max(finite_maxima) if finite_maxima else 1.0
    vmax_meV = max(float(vmax_meV), 1.0e-6)
    norm = Normalize(vmin=0.0, vmax=vmax_meV)
    cmap = plt.get_cmap("viridis")

    cax = ax.inset_axes([0.055, 0.86, 0.33, 0.024])
    sm = ScalarMappable(norm=norm, cmap=cmap)
    cbar = plt.colorbar(sm, cax=cax, orientation="horizontal")
    cbar.set_label(r"$|\Delta|$ [meV]", labelpad=4.0)
    cbar.ax.xaxis.set_label_position("top")
    cbar.ax.xaxis.set_ticks_position("top")
    ticks = np.linspace(0.0, vmax_meV, 5)
    cbar.set_ticks(ticks)
    cbar.set_ticklabels([f"{tick:.2f}" for tick in ticks])
    cbar.ax.tick_params(length=2.2, pad=1.5)

    for pos, inset, field_meV in zip(positions, delta_insets, delta_fields_meV):
        ax_in = ax.inset_axes(pos)
        _draw_delta_inset(ax_in, inset, field_meV=field_meV, norm=norm, cmap=cmap)


def _add_terminal_delta_panels(
    fig: plt.Figure,
    subplot_spec,
    delta_insets: Sequence[Mapping[str, Any]],
) -> None:
    """Draw three vertically stacked full-strip gap maps beside the terminal IV curve."""
    grid = subplot_spec.subgridspec(
        4,
        1,
        height_ratios=(0.16, 1.0, 1.0, 1.0),
        hspace=0.19,
    )
    delta_fields_meV = [_extract_delta_field_meV(item.get("dataset", {})) for item in delta_insets]
    finite_maxima = [
        float(np.nanmax(field))
        for field in delta_fields_meV
        if field.size and np.any(np.isfinite(field))
    ]
    vmax_meV = max(max(finite_maxima) if finite_maxima else 1.0, 1.0e-6)
    norm = Normalize(vmin=0.0, vmax=vmax_meV)
    cmap = plt.get_cmap("viridis")

    cax = fig.add_subplot(grid[0, 0])
    cbar = fig.colorbar(ScalarMappable(norm=norm, cmap=cmap), cax=cax, orientation="horizontal")
    cbar.set_label(r"$|\Delta|$ [meV]", labelpad=3.0)
    cbar.ax.xaxis.set_label_position("top")
    cbar.ax.xaxis.set_ticks_position("top")
    cbar.ax.tick_params(length=2.2, pad=1.5)

    for row, (inset, field_meV) in enumerate(zip(delta_insets, delta_fields_meV), start=1):
        panel = fig.add_subplot(grid[row, 0])
        _draw_full_delta_panel(
            panel,
            inset,
            field_meV=field_meV,
            norm=norm,
            cmap=cmap,
            show_x_label=row == 3,
            show_y_label=True,
        )


def _draw_full_delta_panel(
    ax: plt.Axes,
    inset: Mapping[str, Any],
    *,
    field_meV: np.ndarray,
    norm: Normalize,
    cmap,
    show_x_label: bool,
    show_y_label: bool,
) -> None:
    dataset = inset.get("dataset", {})
    x_nm = np.asarray(dataset.get("x_nm", []), dtype=float)
    y_nm = np.asarray(dataset.get("y_nm", []), dtype=float)
    triangles = np.asarray(dataset.get("triangles", []), dtype=np.int64)
    if x_nm.size == 0 or y_nm.size == 0 or triangles.size == 0 or field_meV.size != x_nm.size:
        ax.text(0.5, 0.5, r"missing $|\Delta|$ data", ha="center", va="center", transform=ax.transAxes)
        ax.set_axis_off()
        return

    triang = mtri.Triangulation(x_nm, y_nm, triangles)
    ax.tripcolor(
        triang,
        field_meV,
        shading="gouraud",
        cmap=cmap,
        norm=norm,
        rasterized=True,
    )
    ax.set_xlim(float(np.nanmin(x_nm)), float(np.nanmax(x_nm)))
    ax.set_ylim(float(np.nanmin(y_nm)), float(np.nanmax(y_nm)))
    ax.set_aspect("equal", adjustable="box")
    ax.tick_params(axis="both", labelsize="x-small", length=2.5)
    if show_x_label:
        ax.set_xlabel(r"$x$ [nm]", labelpad=1.5)
    else:
        ax.tick_params(axis="x", labelbottom=False)
    if show_y_label:
        ax.set_ylabel(r"$y$ [nm]", labelpad=1.5)

    index = int(inset.get("index", 0))
    current = float(inset.get("actual_current_uA", np.nan))
    requested = float(inset.get("requested_current_uA", np.nan))
    label = rf"#{index}  {current:.0f} [$\mu$A]"
    if np.isfinite(requested) and abs(requested - current) > 0.05:
        label += rf" (requested {requested:.0f} [$\mu$A])"
    ax.text(
        0.985,
        0.93,
        label,
        ha="right",
        va="top",
        transform=ax.transAxes,
        fontsize="x-small",
        color="white",
        bbox={
            "boxstyle": "round,pad=0.18",
            "facecolor": "red",
            "edgecolor": "white",
            "linewidth": 0.6,
            "alpha": 0.84,
        },
        zorder=10,
    )



def _highlight_snapshot_points(
    ax: plt.Axes,
    points: Sequence[Mapping[str, Any]],
    delta_insets: Sequence[Mapping[str, Any]],
    *,
    voltage_key: str = "voltage_mV",
):
    xs = []
    ys = []
    indices = []
    for inset in delta_insets:
        x = float(inset.get("actual_current_uA", np.nan))
        y = _lookup_voltage(points, x, voltage_key=voltage_key)
        if np.isfinite(x) and np.isfinite(y):
            xs.append(x)
            ys.append(y)
            indices.append(int(inset.get("index", 0)))
    if not xs:
        return None

    handle = ax.scatter(
        xs,
        ys,
        s=74.0,
        facecolors="red",
        edgecolors="black",
        linewidths=0.9,
        zorder=5.0,
    )
    for x, y, idx in zip(xs, ys, indices):
        ax.text(
            x,
            y,
            str(idx),
            ha="center",
            va="center",
            fontsize="x-small",
            color="white",
            fontweight="bold",
            zorder=6.0,
        )
    return handle



def _extract_delta_field_meV(dataset: Mapping[str, Any]) -> np.ndarray:
    """Return |Delta| in meV using the best available dataset keys."""
    direct_keys = (
        "delta_meV",
        "delta_abs_meV",
        "delta_magnitude_meV",
        "abs_delta_meV",
        "delta_mod_meV",
        "delta_mag_meV",
    )
    for key in direct_keys:
        if key in dataset:
            arr = np.asarray(dataset.get(key, []), dtype=float)
            if arr.size:
                return np.abs(arr)

    profiles = dataset.get("profiles", {})
    if isinstance(profiles, Mapping):
        for key in direct_keys:
            if key in profiles:
                arr = np.asarray(profiles.get(key, []), dtype=float)
                if arr.size:
                    return np.abs(arr)

    delta_over = np.asarray(dataset.get("delta_over_delta0", []), dtype=float)
    if delta_over.size:
        summary_scalars = dataset.get("summary_scalars", {})
        if not isinstance(summary_scalars, Mapping):
            summary_scalars = {}
        delta0_meV = _find_numeric_recursive(
            summary_scalars,
            keys=("delta0_meV", "delta_eq_meV", "delta_ref_meV", "Delta0_meV", "gap0_meV"),
        )
        if not np.isfinite(delta0_meV):
            delta0_J = _find_numeric_recursive(
                summary_scalars,
                keys=("delta0_J", "delta_eq_J", "Delta0_J", "gap0_J"),
            )
            if np.isfinite(delta0_J):
                delta0_meV = float(delta0_J / MEV_J)
        if np.isfinite(delta0_meV) and delta0_meV > 0.0:
            return np.abs(delta_over) * float(delta0_meV)
        # Fall back to unitless scale if nothing else is available; label will still be meV-scale placeholder.
        return np.abs(delta_over)

    return np.array([], dtype=float)



def _draw_delta_inset(
    ax: plt.Axes,
    inset: Mapping[str, Any],
    *,
    field_meV: np.ndarray,
    norm: Normalize,
    cmap,
) -> None:
    dataset = inset.get("dataset", {})
    x_nm = np.asarray(dataset.get("x_nm", []), dtype=float)
    y_nm = np.asarray(dataset.get("y_nm", []), dtype=float)
    triangles = np.asarray(dataset.get("triangles", []), dtype=np.int64)
    if x_nm.size == 0 or y_nm.size == 0 or triangles.size == 0 or field_meV.size != x_nm.size:
        ax.text(0.5, 0.5, "missing\n|Δ| data", ha="center", va="center", transform=ax.transAxes, fontsize="x-small")
        ax.set_xticks([])
        ax.set_yticks([])
        return

    x_center = 0.5 * (float(np.nanmin(x_nm)) + float(np.nanmax(x_nm)))
    x_left = x_center - 50.0
    x_right = x_center + 50.0
    y_min = float(np.nanmin(y_nm))
    y_max = float(np.nanmax(y_nm))

    triang = mtri.Triangulation(x_nm, y_nm, triangles)
    ax.tripcolor(
        triang,
        field_meV,
        shading="gouraud",
        cmap=cmap,
        norm=norm,
        rasterized=True,
    )
    ax.set_xlim(x_left, x_right)
    ax.set_ylim(y_min, y_max)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xticks([])
    ax.set_yticks([])

    idx = int(inset.get("index", 0))
    current = float(inset.get("actual_current_uA", np.nan))
    requested = float(inset.get("requested_current_uA", np.nan))
    if np.isfinite(requested) and abs(requested - current) > 0.05:
        label = f"#{idx}  {current:.0f} [μA]\n(req {requested:.0f} [μA])"
    else:
        label = f"#{idx}  {current:.0f} [μA]"
    ax.text(
        0.965,
        0.965,
        label,
        ha="right",
        va="top",
        transform=ax.transAxes,
        fontsize="x-small",
        color="white",
        bbox={"boxstyle": "round,pad=0.18", "facecolor": "red", "edgecolor": "white", "linewidth": 0.6, "alpha": 0.82},
        zorder=10,
    )
    for spine in ax.spines.values():
        spine.set_linewidth(0.7)
        spine.set_edgecolor("0.1")



def _lookup_voltage(
    points: Sequence[Mapping[str, Any]],
    current_uA: float,
    *,
    voltage_key: str = "voltage_mV",
) -> float:
    for item in points:
        if np.isfinite(float(item.get("current_uA", np.nan))) and abs(float(item.get("current_uA")) - current_uA) < 1.0e-9:
            return float(item.get(voltage_key, np.nan))
    return np.nan



def _window_or_interp(x: np.ndarray, y: np.ndarray, *, center: float, half_window: float) -> float:
    mask = np.isfinite(x) & np.isfinite(y) & (np.abs(x - center) <= half_window)
    if np.any(mask):
        return float(np.nanmean(y[mask]))
    order = np.argsort(x)
    xs = np.asarray(x[order], dtype=float)
    ys = np.asarray(y[order], dtype=float)
    finite = np.isfinite(xs) & np.isfinite(ys)
    if np.count_nonzero(finite) < 2:
        raise ValueError("insufficient finite x-profile samples for voltage interpolation.")
    return float(np.interp(float(center), xs[finite], ys[finite]))


def _last_finite(value: Any) -> float:
    values = np.asarray(value if value is not None else [], dtype=float).reshape(-1)
    finite = values[np.isfinite(values)]
    return float(finite[-1]) if finite.size else np.nan



def _default_profile_half_window_nm(x: np.ndarray) -> float:
    xs = np.asarray(x, dtype=float)
    diffs = np.diff(np.unique(xs[np.isfinite(xs)]))
    diffs = diffs[diffs > 0.0]
    if diffs.size == 0:
        return 1.0
    return float(max(0.55 * np.nanmedian(diffs), 1.0))



def _infer_bias_current_uA(*, run_name: str, summary: Mapping[str, Any], dataset: Mapping[str, Any]) -> float:
    scalar = _find_first_numeric(
        summary,
        keys=(
            "target_current_A", "current_A", "I_bias_A", "bias_current_A",
            "target_current_uA", "current_uA", "I_bias_uA", "bias_current_uA",
        ),
    )
    if scalar is not None:
        key, value = scalar
        if key.endswith("_uA"):
            return float(value)
        return 1.0e6 * float(value)

    summary_scalars = dataset.get("summary_scalars", {})
    if isinstance(summary_scalars, Mapping) and "target_current_A" in summary_scalars:
        try:
            return 1.0e6 * float(summary_scalars["target_current_A"])
        except Exception:
            pass

    match = re.search(r"(?:^|_)I(?P<i>[-+]?\d+(?:\.\d+)?)uA(?:_|$)", run_name)
    if match:
        return float(match.group("i"))
    match = re.search(r"base(?P<i>[-+]?\d+(?:\.\d+)?)uA", run_name)
    if match:
        return float(match.group("i"))
    raise ValueError(f"Could not infer bias current from run '{run_name}'.")



def _find_first_numeric(obj: Any, *, keys: Sequence[str]) -> tuple[str, float] | None:
    if isinstance(obj, Mapping):
        for key in keys:
            if key in obj:
                try:
                    return str(key), float(obj[key])
                except Exception:
                    pass
        for value in obj.values():
            found = _find_first_numeric(value, keys=keys)
            if found is not None:
                return found
    elif isinstance(obj, (list, tuple)):
        for value in obj:
            found = _find_first_numeric(value, keys=keys)
            if found is not None:
                return found
    return None



def _find_numeric_recursive(obj: Any, *, keys: Sequence[str]) -> float:
    found = _find_first_numeric(obj, keys=keys)
    return float(found[1]) if found is not None else np.nan



def _orient_positive_voltage(
    points: Sequence[dict[str, Any]],
    *,
    voltage_key: str = "voltage_mV",
) -> bool:
    voltages = np.asarray([float(item.get(voltage_key, np.nan)) for item in points], dtype=float)
    currents = np.asarray([float(item.get("current_uA", np.nan)) for item in points], dtype=float)
    mask = np.isfinite(voltages) & np.isfinite(currents) & (currents > 0.0) & (np.abs(voltages) > 0.0)
    if not np.any(mask):
        return False
    median_v = float(np.nanmedian(voltages[mask]))
    if median_v >= 0.0:
        return False
    for item in points:
        try:
            item[voltage_key] = -float(item.get(voltage_key, np.nan))
        except Exception:
            pass
    return True



def _to_builtin(value: Any) -> Any:
    if isinstance(value, Mapping):
        return {str(k): _to_builtin(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [_to_builtin(v) for v in value]
    if isinstance(value, (np.floating, np.integer, np.bool_)):
        return value.item()
    if isinstance(value, np.ndarray):
        return value.tolist()
    return value
