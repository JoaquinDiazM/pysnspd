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
from pysnspd.plotting.style import THESIS_DOUBLE_FIGSIZE, THESIS_DPI, apply_thesis_style

apply_thesis_style()

MEV_J = 1.602176634e-22


def plot_current_sweep_iv(
    points: Sequence[Mapping[str, Any]],
    output_path: str | Path,
    *,
    dpi: int = THESIS_DPI,
    voltage_probe_offset_nm: float = 50.0,
    include_origin: bool = True,
    delta_insets: Sequence[Mapping[str, Any]] | None = None,
) -> Path:
    """Plot current on x-axis and central TDGL voltage on y-axis."""
    output = Path(output_path)
    output.parent.mkdir(parents=True, exist_ok=True)

    x_data, y_data = _extract_iv_arrays(points)
    valid = np.isfinite(x_data) & np.isfinite(y_data)
    x_valid = x_data[valid]
    y_valid = y_data[valid]

    normal_x, normal_y = _extract_normal_curve(points)
    fit_x, fit_y = _monotone_fit_curve(x_valid, y_valid)

    fig, ax = plt.subplots(figsize=THESIS_DOUBLE_FIGSIZE, constrained_layout=False)
    fig.subplots_adjust(left=0.105, right=0.975, bottom=0.125, top=0.965)

    fit_line = None
    normal_line = None
    if fit_x.size:
        fit_line, = ax.plot(
            fit_x,
            fit_y,
            linewidth=1.55,
            color="black",
            label="Monotone fit",
            zorder=2.0,
        )
    if normal_x.size:
        normal_line, = ax.plot(
            normal_x,
            normal_y,
            linestyle=(0, (5.5, 2.6)),
            linewidth=1.9,
            color="tab:orange",
            alpha=0.98,
            label="Ohmic behavior",
            zorder=2.6,
        )
    raw_scatter = ax.scatter(
        x_valid,
        y_valid,
        s=22.0,
        color="tab:blue",
        label="Raw data points",
        zorder=3.0,
    )

    if include_origin:
        ax.scatter([0.0], [0.0], s=26.0, color="tab:orange", zorder=4.0)

    snapshot_handle = None
    if delta_insets:
        snapshot_handle = _highlight_snapshot_points(ax, points, delta_insets)
        _add_delta_insets(ax, delta_insets)

    ax.set_xlabel(r"$I_{\mathrm{bias}}$ [$\mu$A]")
    ax.set_ylabel(r"$V_{\mathrm{TDGL}}$ [mV]")
    ax.tick_params(axis="both", which="major")
    ax.grid(True, linewidth=0.45, alpha=0.33)

    if x_valid.size:
        xmin = float(np.nanmin(x_valid))
        xmax = float(np.nanmax(x_valid))
        dx = max(xmax - xmin, 1.0)
        ax.set_xlim(min(0.0, xmin) - 0.03 * dx, xmax + 0.04 * dx)
    if y_valid.size or normal_y.size or fit_y.size:
        all_y = np.concatenate([arr for arr in (y_valid, normal_y, fit_y) if arr.size])
        ymin = float(np.nanmin(all_y))
        ymax = float(np.nanmax(all_y))
        dy = max(ymax - ymin, 1.0e-6)
        lower = min(0.0, ymin - 0.06 * dy)
        upper = max(0.0, ymax + 0.08 * dy)
        if upper <= lower:
            upper = lower + 1.0
        ax.set_ylim(lower, upper)

    handles = [raw_scatter]
    labels = ["Raw data points"]
    if fit_line is not None:
        handles.append(fit_line)
        labels.append("Monotone fit")
    if normal_line is not None:
        handles.append(normal_line)
        labels.append("Ohmic behavior")
    if snapshot_handle is not None:
        handles.append(snapshot_handle)
        labels.append("Order-parameter snapshots")

    legend = ax.legend(
        handles,
        labels,
        loc="lower right",
        frameon=True,
    )
    legend.get_frame().set_alpha(0.95)

    fig.savefig(output, dpi=dpi, bbox_inches="tight", pad_inches=0.06)
    plt.close(fig)
    return output


def plot_terminal_current_sweep_iv(
    points: Sequence[Mapping[str, Any]],
    output_path: str | Path,
    *,
    dpi: int = THESIS_DPI,
    include_origin: bool = True,
    delta_insets: Sequence[Mapping[str, Any]] | None = None,
) -> Path:
    """Plot terminal voltage and three full-device gap snapshots."""
    output = Path(output_path)
    output.parent.mkdir(parents=True, exist_ok=True)

    x_data, y_data = _extract_iv_arrays(points, voltage_key="terminal_voltage_mV")
    valid = np.isfinite(x_data) & np.isfinite(y_data)
    x_valid = x_data[valid]
    y_valid = y_data[valid]
    normal_x, normal_y = _extract_normal_curve(
        points,
        voltage_key="normal_terminal_voltage_mV",
    )
    fit_x, fit_y = _monotone_fit_curve(x_valid, y_valid)

    fig = plt.figure(figsize=THESIS_DOUBLE_FIGSIZE, constrained_layout=False)
    outer = fig.add_gridspec(
        1,
        2,
        width_ratios=(1.7, 1.0),
        left=0.075,
        right=0.975,
        bottom=0.105,
        top=0.955,
        wspace=0.18,
    )
    ax = fig.add_subplot(outer[0, 0])

    fit_line = None
    normal_line = None
    if fit_x.size:
        fit_line, = ax.plot(
            fit_x,
            fit_y,
            linewidth=1.55,
            color="black",
            label="Monotone fit",
            zorder=2.0,
        )
    if normal_x.size:
        normal_line, = ax.plot(
            normal_x,
            normal_y,
            linestyle=(0, (5.5, 2.6)),
            linewidth=1.9,
            color="tab:orange",
            alpha=0.98,
            label="Ohmic behavior",
            zorder=2.6,
        )
    raw_scatter = ax.scatter(
        x_valid,
        y_valid,
        s=25.0,
        color="tab:blue",
        label="Raw data points",
        zorder=3.0,
    )
    if include_origin:
        ax.scatter([0.0], [0.0], s=29.0, color="tab:orange", zorder=4.0)

    snapshot_handle = None
    if delta_insets:
        snapshot_handle = _highlight_snapshot_points(
            ax,
            points,
            delta_insets,
            voltage_key="terminal_voltage_mV",
        )
        _add_terminal_delta_panels(fig, outer[0, 1], delta_insets)
    else:
        empty_ax = fig.add_subplot(outer[0, 1])
        empty_ax.axis("off")

    ax.set_xlabel(r"$I_{\mathrm{bias}}$ [$\mu$A]")
    ax.set_ylabel(r"$V_{\mathrm{terminal}}$ [mV]")
    ax.tick_params(axis="both", which="major")
    ax.grid(True, linewidth=0.45, alpha=0.33)
    _set_iv_limits(ax, x_valid, y_valid, normal_y, fit_y)

    handles = [raw_scatter]
    labels = ["Raw data points"]
    if fit_line is not None:
        handles.append(fit_line)
        labels.append("Monotone fit")
    if normal_line is not None:
        handles.append(normal_line)
        labels.append("Ohmic behavior")
    if snapshot_handle is not None:
        handles.append(snapshot_handle)
        labels.append("Order-parameter snapshots")
    legend = ax.legend(handles, labels, loc="lower right", frameon=True)
    legend.get_frame().set_alpha(0.95)

    fig.savefig(output, dpi=dpi, bbox_inches="tight", pad_inches=0.06)
    plt.close(fig)
    return output


def _set_iv_limits(
    ax: plt.Axes,
    current_uA: np.ndarray,
    voltage_mV: np.ndarray,
    normal_voltage_mV: np.ndarray,
    fit_voltage_mV: np.ndarray,
) -> None:
    if current_uA.size:
        xmin = float(np.nanmin(current_uA))
        xmax = float(np.nanmax(current_uA))
        dx = max(xmax - xmin, 1.0)
        ax.set_xlim(min(0.0, xmin) - 0.03 * dx, xmax + 0.04 * dx)
    arrays = [arr for arr in (voltage_mV, normal_voltage_mV, fit_voltage_mV) if arr.size]
    if not arrays:
        return
    all_y = np.concatenate(arrays)
    ymin = float(np.nanmin(all_y))
    ymax = float(np.nanmax(all_y))
    dy = max(ymax - ymin, 1.0e-6)
    lower = min(0.0, ymin - 0.06 * dy)
    upper = max(0.0, ymax + 0.08 * dy)
    ax.set_ylim(lower, upper if upper > lower else lower + 1.0)



def write_current_sweep_iv_csv(points: Sequence[Mapping[str, Any]], output_path: str | Path) -> Path:
    path = Path(output_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "run_name", "current_uA", "voltage_mV", "normal_voltage_mV", "normal_resistance_probe_ohm",
        "terminal_voltage_mV", "normal_terminal_voltage_mV", "normal_resistance_terminal_ohm",
        "probe_left_x_nm", "probe_right_x_nm", "probe_left_phi_mV", "probe_right_phi_mV",
        "profile_x_center_nm", "profile_x_min_nm", "profile_x_max_nm",
        "voltage_probe_offset_nm", "voltage_probe_half_window_nm",
        "pre_run_name", "raw_ss", "source",
    ]
    with path.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for item in points:
            writer.writerow({key: item.get(key, "") for key in fieldnames})
    return path



def write_current_sweep_iv_yaml(points: Sequence[Mapping[str, Any]], meta: Mapping[str, Any], output_path: str | Path) -> Path:
    path = Path(output_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    payload = {"meta": dict(meta), "points": [_to_builtin(item) for item in points]}
    with path.open("w", encoding="utf-8") as f:
        yaml.safe_dump(payload, f, sort_keys=False, allow_unicode=True, default_flow_style=False)
    return path



def write_skipped_runs_yaml(skipped: Sequence[Mapping[str, Any]], output_path: str | Path) -> Path:
    path = Path(output_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as f:
        yaml.safe_dump([_to_builtin(item) for item in skipped], f, sort_keys=False, allow_unicode=True, default_flow_style=False)
    return path



def write_iv_insets_yaml(insets: Sequence[Mapping[str, Any]], output_path: str | Path) -> Path:
    path = Path(output_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    serializable = []
    for item in insets:
        serializable.append({k: _to_builtin(v) for k, v in item.items() if k != "dataset"})
    with path.open("w", encoding="utf-8") as f:
        yaml.safe_dump(serializable, f, sort_keys=False, allow_unicode=True, default_flow_style=False)
    return path



def _build_iv_point(
    *,
    run_name: str,
    run: Any,
    dataset: Mapping[str, Any],
    project_config: Mapping[str, Any] | None,
    voltage_probe_offset_nm: float,
    voltage_probe_half_window_nm: float | None,
) -> tuple[dict[str, Any], float, float, float]:
    current_uA = _infer_bias_current_uA(run_name=run_name, summary=getattr(run, "summary", {}), dataset=dataset)
    x_profile = np.asarray(dataset.get("x_profile_nm", []), dtype=float)
    profiles = dataset.get("profiles", {})
    has_profile = (
        x_profile.size > 0
        and isinstance(profiles, Mapping)
        and np.asarray(profiles.get("phi_mV", [])).size > 0
    )
    if has_profile:
        left_phi, right_phi, left_x, right_x, half_window_nm = _extract_profile_probe_values(
            dataset,
            voltage_probe_offset_nm=voltage_probe_offset_nm,
            voltage_probe_half_window_nm=voltage_probe_half_window_nm,
        )
        central_voltage_mV = float(right_phi - left_phi)
        central_source = "x_profile_phi_mV"
    else:
        central_voltage_mV = _last_finite(dataset.get("tdgl_probe_voltage_mV"))
        left_x = float(dataset.get("tdgl_probe_left_x_nm", np.nan))
        right_x = float(dataset.get("tdgl_probe_right_x_nm", np.nan))
        left_phi = np.nan
        right_phi = np.nan
        half_window_nm = (
            float(voltage_probe_half_window_nm)
            if voltage_probe_half_window_nm is not None
            else np.nan
        )
        central_source = "tdgl_probe_voltage_mV_history_final"

    spatial_x_nm = x_profile if x_profile.size else np.asarray(dataset.get("x_nm", []), dtype=float)
    rn_probe_ohm = _infer_probe_normal_resistance_ohm(
        run=run,
        dataset=dataset,
        project_config=project_config,
        voltage_probe_offset_nm=voltage_probe_offset_nm,
    )
    normal_voltage_mV = 1.0e-3 * float(current_uA) * float(rn_probe_ohm) if np.isfinite(rn_probe_ohm) else np.nan
    rn_terminal_ohm = _infer_terminal_normal_resistance_ohm(
        run=run,
        dataset=dataset,
        project_config=project_config,
    )
    terminal_voltage_mV = _last_finite(dataset.get("terminal_voltage_mV"))
    if not np.isfinite(terminal_voltage_mV):
        terminal_voltage_V = _find_numeric_recursive(
            getattr(run, "summary", {}),
            keys=("terminal_voltage_V",),
        )
        terminal_voltage_mV = 1.0e3 * terminal_voltage_V if np.isfinite(terminal_voltage_V) else np.nan
    normal_terminal_voltage_mV = (
        1.0e-3 * float(current_uA) * float(rn_terminal_ohm)
        if np.isfinite(rn_terminal_ohm)
        else np.nan
    )

    point = {
        "run_name": run_name,
        "current_uA": float(current_uA),
        "voltage_mV": float(central_voltage_mV),
        "normal_voltage_mV": float(normal_voltage_mV),
        "normal_resistance_probe_ohm": float(rn_probe_ohm),
        "terminal_voltage_mV": float(terminal_voltage_mV),
        "normal_terminal_voltage_mV": float(normal_terminal_voltage_mV),
        "normal_resistance_terminal_ohm": float(rn_terminal_ohm),
        "probe_left_x_nm": float(left_x),
        "probe_right_x_nm": float(right_x),
        "probe_left_phi_mV": float(left_phi),
        "probe_right_phi_mV": float(right_phi),
        "profile_x_center_nm": float(0.5 * (np.nanmin(spatial_x_nm) + np.nanmax(spatial_x_nm))) if spatial_x_nm.size else np.nan,
        "profile_x_min_nm": float(np.nanmin(spatial_x_nm)) if spatial_x_nm.size else np.nan,
        "profile_x_max_nm": float(np.nanmax(spatial_x_nm)) if spatial_x_nm.size else np.nan,
        "voltage_probe_offset_nm": float(voltage_probe_offset_nm),
        "voltage_probe_half_window_nm": float(half_window_nm),
        "pre_run_name": getattr(run, "pre_run_name", None),
        "raw_ss": str(getattr(run, "raw_ss", "")) if getattr(run, "raw_ss", None) is not None else None,
        "source": central_source + "_and_terminal_voltage_mV_history_final",
    }
    return point, float(half_window_nm), float(rn_probe_ohm), float(rn_terminal_ohm)



def _extract_profile_probe_values(
    dataset: Mapping[str, Any],
    *,
    voltage_probe_offset_nm: float,
    voltage_probe_half_window_nm: float | None,
) -> tuple[float, float, float, float, float]:
    x = np.asarray(dataset.get("x_profile_nm", []), dtype=float)
    profiles = dataset.get("profiles", {})
    if not isinstance(profiles, Mapping):
        profiles = {}
    phi = np.asarray(profiles.get("phi_mV", []), dtype=float)
    if x.size == 0 or phi.size == 0:
        raise ValueError("dataset does not provide x_profile_nm / profiles['phi_mV'].")
    if phi.size != x.size:
        phi = np.resize(phi, x.size)

    xmin = float(np.nanmin(x))
    xmax = float(np.nanmax(x))
    x_center = 0.5 * (xmin + xmax)
    left_x = x_center - float(voltage_probe_offset_nm)
    right_x = x_center + float(voltage_probe_offset_nm)
    half_window_nm = float(voltage_probe_half_window_nm) if voltage_probe_half_window_nm is not None else _default_profile_half_window_nm(x)
    if half_window_nm <= 0.0:
        raise ValueError("voltage_probe_half_window_nm must be positive.")

    left_phi = _window_or_interp(x, phi, center=left_x, half_window=half_window_nm)
    right_phi = _window_or_interp(x, phi, center=right_x, half_window=half_window_nm)
    return float(left_phi), float(right_phi), float(left_x), float(right_x), float(half_window_nm)



def _infer_probe_normal_resistance_ohm(
    *,
    run: Any,
    dataset: Mapping[str, Any],
    project_config: Mapping[str, Any] | None,
    voltage_probe_offset_nm: float,
) -> float:
    summary = getattr(run, "summary", {})
    material_cfg = project_config.get("material", {}) if isinstance(project_config, Mapping) else {}
    if not isinstance(material_cfg, Mapping):
        material_cfg = {}

    width_candidates = [
        getattr(run.mesh, "width_m", np.nan),
        _find_numeric_recursive(summary, keys=("width_m", "wire_width_m", "w_m", "device_width_m")),
        _find_numeric_recursive(dataset.get("summary_scalars", {}), keys=("width_m", "wire_width_m", "w_m", "device_width_m")),
        _find_numeric_recursive(material_cfg, keys=("width_m", "wire_width_m", "w_m")),
    ]
    thickness_candidates = [
        _find_numeric_recursive(summary, keys=("thickness_m", "thickness", "d_m", "film_thickness_m")),
        _find_numeric_recursive(dataset.get("summary_scalars", {}), keys=("thickness_m", "thickness", "d_m", "film_thickness_m")),
        _find_numeric_recursive(material_cfg, keys=("thickness_m", "thickness", "d_m", "film_thickness_m")),
    ]
    sigma_candidates = [
        _find_numeric_recursive(summary, keys=("sigma_n_S_m", "sigma_n_S_per_m", "sigma_n", "normal_conductivity_S_m")),
        _find_numeric_recursive(dataset.get("summary_scalars", {}), keys=("sigma_n_S_m", "sigma_n_S_per_m", "sigma_n", "normal_conductivity_S_m")),
        _find_numeric_recursive(material_cfg, keys=("sigma_n_S_m", "sigma_n_S_per_m", "sigma_n", "normal_conductivity_S_m")),
    ]

    width_m = next((float(v) for v in width_candidates if np.isfinite(v) and float(v) > 0.0), np.nan)
    thickness_m = next((float(v) for v in thickness_candidates if np.isfinite(v) and float(v) > 0.0), np.nan)
    sigma_n = next((float(v) for v in sigma_candidates if np.isfinite(v) and float(v) > 0.0), np.nan)
    if not np.isfinite(width_m) or not np.isfinite(thickness_m) or not np.isfinite(sigma_n):
        return np.nan

    length_m = 2.0 * float(voltage_probe_offset_nm) * 1.0e-9
    return float(length_m / (float(sigma_n) * width_m * float(thickness_m)))


def _infer_terminal_normal_resistance_ohm(
    *,
    run: Any,
    dataset: Mapping[str, Any],
    project_config: Mapping[str, Any] | None,
) -> float:
    """Return the normal resistance for the complete simulated strip."""
    summary = getattr(run, "summary", {})
    material_cfg = project_config.get("material", {}) if isinstance(project_config, Mapping) else {}
    if not isinstance(material_cfg, Mapping):
        material_cfg = {}

    width_candidates = [
        getattr(run.mesh, "width_m", np.nan),
        _find_numeric_recursive(summary, keys=("width_m", "wire_width_m", "w_m", "device_width_m")),
        _find_numeric_recursive(dataset.get("summary_scalars", {}), keys=("width_m", "wire_width_m", "w_m", "device_width_m")),
        _find_numeric_recursive(material_cfg, keys=("width_m", "wire_width_m", "w_m")),
    ]
    thickness_candidates = [
        _find_numeric_recursive(summary, keys=("thickness_m", "thickness", "d_m", "film_thickness_m")),
        _find_numeric_recursive(dataset.get("summary_scalars", {}), keys=("thickness_m", "thickness", "d_m", "film_thickness_m")),
        _find_numeric_recursive(material_cfg, keys=("thickness_m", "thickness", "d_m", "film_thickness_m")),
    ]
    sigma_candidates = [
        _find_numeric_recursive(summary, keys=("sigma_n_S_m", "sigma_n_S_per_m", "sigma_n", "normal_conductivity_S_m")),
        _find_numeric_recursive(dataset.get("summary_scalars", {}), keys=("sigma_n_S_m", "sigma_n_S_per_m", "sigma_n", "normal_conductivity_S_m")),
        _find_numeric_recursive(material_cfg, keys=("sigma_n_S_m", "sigma_n_S_per_m", "sigma_n", "normal_conductivity_S_m")),
    ]

    mesh = getattr(run, "mesh", None)
    length_candidates = [getattr(mesh, "length_m", np.nan)]
    nodes = np.asarray(getattr(mesh, "nodes", []), dtype=float)
    if nodes.ndim == 2 and nodes.shape[0] and nodes.shape[1]:
        length_candidates.append(float(np.nanmax(nodes[:, 0]) - np.nanmin(nodes[:, 0])))
    x_nm = np.asarray(dataset.get("x_nm", []), dtype=float)
    if x_nm.size:
        length_candidates.append(1.0e-9 * float(np.nanmax(x_nm) - np.nanmin(x_nm)))

    width_m = next((float(v) for v in width_candidates if np.isfinite(v) and float(v) > 0.0), np.nan)
    thickness_m = next((float(v) for v in thickness_candidates if np.isfinite(v) and float(v) > 0.0), np.nan)
    sigma_n = next((float(v) for v in sigma_candidates if np.isfinite(v) and float(v) > 0.0), np.nan)
    length_m = next((float(v) for v in length_candidates if np.isfinite(v) and float(v) > 0.0), np.nan)
    if not all(np.isfinite(v) for v in (length_m, width_m, thickness_m, sigma_n)):
        return np.nan
    return float(length_m / (sigma_n * width_m * thickness_m))



def _extract_iv_arrays(
    points: Sequence[Mapping[str, Any]],
    *,
    voltage_key: str = "voltage_mV",
) -> tuple[np.ndarray, np.ndarray]:
    current_uA = np.asarray([float(item.get("current_uA", np.nan)) for item in points], dtype=float)
    voltage_mV = np.asarray([float(item.get(voltage_key, np.nan)) for item in points], dtype=float)
    return current_uA, voltage_mV



def _extract_normal_curve(
    points: Sequence[Mapping[str, Any]],
    *,
    voltage_key: str = "normal_voltage_mV",
) -> tuple[np.ndarray, np.ndarray]:
    current_uA = np.asarray([float(item.get("current_uA", np.nan)) for item in points], dtype=float)
    normal_mV = np.asarray([float(item.get(voltage_key, np.nan)) for item in points], dtype=float)
    mask = np.isfinite(current_uA) & np.isfinite(normal_mV)
    if not np.any(mask):
        return np.array([], dtype=float), np.array([], dtype=float)
    order = np.argsort(current_uA[mask])
    return current_uA[mask][order], normal_mV[mask][order]



def _monotone_fit_curve(x: np.ndarray, y: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    if x.size == 0 or y.size == 0:
        return np.array([], dtype=float), np.array([], dtype=float)
    order = np.argsort(x)
    xs = np.asarray(x[order], dtype=float)
    ys = np.asarray(y[order], dtype=float)
    fitted = _isotonic_regression(ys)
    x_dense = np.linspace(float(xs[0]), float(xs[-1]), max(400, 10 * xs.size))
    y_dense = np.interp(x_dense, xs, fitted)
    return x_dense, y_dense



def _isotonic_regression(y: np.ndarray, w: np.ndarray | None = None) -> np.ndarray:
    y = np.asarray(y, dtype=float).reshape(-1)
    n = y.size
    if n == 0:
        return y.copy()
    if w is None:
        w = np.ones(n, dtype=float)
    else:
        w = np.asarray(w, dtype=float).reshape(-1)
        if w.size != n:
            raise ValueError("weights must have the same length as y.")

    block_start = list(range(n))
    block_end = list(range(n))
    block_value = [float(v) for v in y]
    block_weight = [float(max(ww, 1.0e-15)) for ww in w]

    i = 0
    while i < len(block_value) - 1:
        if block_value[i] <= block_value[i + 1] + 1.0e-15:
            i += 1
            continue
        new_weight = block_weight[i] + block_weight[i + 1]
        new_value = (block_weight[i] * block_value[i] + block_weight[i + 1] * block_value[i + 1]) / new_weight
        block_end[i] = block_end[i + 1]
        block_weight[i] = new_weight
        block_value[i] = new_value
        del block_start[i + 1]
        del block_end[i + 1]
        del block_weight[i + 1]
        del block_value[i + 1]
        if i > 0:
            i -= 1

    fitted = np.empty(n, dtype=float)
    for start, end, value in zip(block_start, block_end, block_value):
        fitted[start : end + 1] = float(value)
    return fitted

from pysnspd.plotting.current_sweep_insets import (
    _add_delta_insets,
    _add_terminal_delta_panels,
    _default_profile_half_window_nm,
    _find_numeric_recursive,
    _highlight_snapshot_points,
    _infer_bias_current_uA,
    _last_finite,
    _orient_positive_voltage,
    _to_builtin,
    _window_or_interp,
)
