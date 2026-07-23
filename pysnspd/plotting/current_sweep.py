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


def _load_project_config(config_path: str | Path) -> dict[str, Any]:
    path = Path(config_path)
    with path.open("r", encoding="utf-8") as f:
        data = yaml.safe_load(f) or {}
    return data if isinstance(data, dict) else {}


def make_current_sweep_figures(
    *,
    config_path: str | Path,
    records: Sequence[Mapping[str, Any]],
    output_dir: str | Path,
    dpi: int = THESIS_DPI,
    voltage_probe_offset_nm: float = 50.0,
    voltage_probe_half_window_nm: float | None = None,
    include_origin: bool = True,
    delta_inset_currents_uA: Sequence[float] | None = None,
    terminal_delta_inset_currents_uA: Sequence[float] | None = None,
) -> dict[str, Any]:
    """Create current-sweep inventory products and figures."""
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)

    config_dict = _load_project_config(config_path)

    points, skipped, meta = collect_current_sweep_iv_points(
        config_path=config_path,
        project_config=config_dict,
        records=records,
        voltage_probe_offset_nm=voltage_probe_offset_nm,
        voltage_probe_half_window_nm=voltage_probe_half_window_nm,
        include_origin=include_origin,
    )
    inset_runs: list[dict[str, Any]] = []
    if delta_inset_currents_uA is not None:
        inset_runs = select_delta_inset_runs(
            config_path=config_path,
            points=points,
            requested_currents_uA=delta_inset_currents_uA,
        )
    terminal_inset_runs: list[dict[str, Any]] = []
    if terminal_delta_inset_currents_uA is not None:
        terminal_inset_runs = select_terminal_delta_inset_runs(
            config_path=config_path,
            points=points,
            requested_currents_uA=terminal_delta_inset_currents_uA,
        )

    saved: dict[str, Any] = {}
    saved["iv_curve"] = plot_current_sweep_iv(
        points,
        out / "Z2_iv_curve.pdf",
        dpi=dpi,
        voltage_probe_offset_nm=voltage_probe_offset_nm,
        include_origin=include_origin,
        delta_insets=inset_runs,
    )
    saved["terminal_iv_curve"] = plot_terminal_current_sweep_iv(
        points,
        out / "Z2_terminal_iv_curve.pdf",
        dpi=dpi,
        include_origin=include_origin,
        delta_insets=terminal_inset_runs,
    )
    saved["iv_points_csv"] = write_current_sweep_iv_csv(points, out / "Z2_iv_points.csv")
    saved["iv_points_yaml"] = write_current_sweep_iv_yaml(points, meta, out / "Z2_iv_points.yaml")
    saved["iv_skipped_yaml"] = write_skipped_runs_yaml(skipped, out / "Z2_iv_skipped.yaml")
    saved["iv_insets_yaml"] = write_iv_insets_yaml(inset_runs, out / "Z2_iv_insets.yaml")
    saved["terminal_iv_insets_yaml"] = write_iv_insets_yaml(
        terminal_inset_runs,
        out / "Z2_terminal_iv_insets.yaml",
    )
    saved["iv_summary"] = {
        "n_points": int(len(points)),
        "n_runs_loaded": int(meta.get("n_runs_loaded", 0)),
        "n_runs_skipped": int(len(skipped)),
        "include_origin": bool(include_origin),
        "voltage_probe_offset_nm": float(voltage_probe_offset_nm),
        "voltage_probe_half_window_nm": float(meta.get("voltage_probe_half_window_nm", np.nan)),
        "voltage_sign_flipped": bool(meta.get("voltage_sign_flipped", False)),
        "terminal_voltage_sign_flipped": bool(
            meta.get("terminal_voltage_sign_flipped", False)
        ),
        "normal_resistance_terminal_ohm": float(
            meta.get("normal_resistance_terminal_ohm", np.nan)
        ),
        "delta_inset_currents_uA": [float(v) for v in delta_inset_currents_uA] if delta_inset_currents_uA is not None else [],
        "delta_inset_resolved_currents_uA": [float(item["actual_current_uA"]) for item in inset_runs],
        "terminal_delta_inset_currents_uA": (
            [float(v) for v in terminal_delta_inset_currents_uA]
            if terminal_delta_inset_currents_uA is not None
            else []
        ),
        "terminal_delta_inset_resolved_currents_uA": [
            float(item["actual_current_uA"]) for item in terminal_inset_runs
        ],
    }
    return saved



def collect_current_sweep_iv_points(
    *,
    config_path: str | Path,
    project_config: Mapping[str, Any] | None,
    records: Sequence[Mapping[str, Any]],
    voltage_probe_offset_nm: float = 50.0,
    voltage_probe_half_window_nm: float | None = None,
    include_origin: bool = True,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]], dict[str, Any]]:
    points: list[dict[str, Any]] = []
    skipped: list[dict[str, Any]] = []
    used_half_window_nm: float | None = None
    normal_resistance_ohm: float | None = None
    normal_resistance_terminal_ohm: float | None = None

    for record in records:
        run_name = str(record.get("run_name", ""))
        stages = record.get("stages", {})
        stage_ss = stages.get("ss", {}) if isinstance(stages, Mapping) else {}
        if not isinstance(stage_ss, Mapping) or not stage_ss.get("exists", False):
            skipped.append({"run_name": run_name, "reason": "ss stage not found"})
            continue
        try:
            run = load_ss_run(config_path=config_path, run_name=run_name)
            dataset = build_ss_plot_dataset(run)
            point, half_window_nm, rn_probe, rn_terminal = _build_iv_point(
                run_name=run_name,
                run=run,
                dataset=dataset,
                project_config=project_config,
                voltage_probe_offset_nm=voltage_probe_offset_nm,
                voltage_probe_half_window_nm=voltage_probe_half_window_nm,
            )
            if used_half_window_nm is None and np.isfinite(half_window_nm):
                used_half_window_nm = float(half_window_nm)
            if normal_resistance_ohm is None and np.isfinite(rn_probe):
                normal_resistance_ohm = float(rn_probe)
            if normal_resistance_terminal_ohm is None and np.isfinite(rn_terminal):
                normal_resistance_terminal_ohm = float(rn_terminal)
            points.append(point)
        except Exception as exc:
            skipped.append({"run_name": run_name, "reason": f"{type(exc).__name__}: {exc}"})

    points.sort(key=lambda item: (float(item.get("current_uA", np.nan)), str(item.get("run_name", ""))))
    sign_flipped = _orient_positive_voltage(points)
    terminal_sign_flipped = _orient_positive_voltage(
        points,
        voltage_key="terminal_voltage_mV",
    )

    if include_origin:
        origin = {
            "run_name": "synthetic_origin",
            "current_uA": 0.0,
            "voltage_mV": 0.0,
            "normal_voltage_mV": 0.0,
            "terminal_voltage_mV": 0.0,
            "normal_terminal_voltage_mV": 0.0,
            "normal_resistance_probe_ohm": float(normal_resistance_ohm if normal_resistance_ohm is not None else np.nan),
            "normal_resistance_terminal_ohm": float(
                normal_resistance_terminal_ohm
                if normal_resistance_terminal_ohm is not None
                else np.nan
            ),
            "probe_left_x_nm": float("nan"),
            "probe_right_x_nm": float("nan"),
            "probe_left_phi_mV": 0.0,
            "probe_right_phi_mV": 0.0,
            "profile_x_center_nm": float("nan"),
            "profile_x_min_nm": float("nan"),
            "profile_x_max_nm": float("nan"),
            "voltage_probe_offset_nm": float(voltage_probe_offset_nm),
            "voltage_probe_half_window_nm": float(used_half_window_nm if used_half_window_nm is not None else np.nan),
            "pre_run_name": None,
            "raw_ss": None,
            "source": "synthetic_origin",
        }
        points = [origin] + points

    meta = {
        "n_runs_loaded": len(points) - (1 if include_origin else 0),
        "voltage_probe_offset_nm": float(voltage_probe_offset_nm),
        "voltage_probe_half_window_nm": float(used_half_window_nm if used_half_window_nm is not None else np.nan),
        "voltage_sign_flipped": bool(sign_flipped),
        "terminal_voltage_sign_flipped": bool(terminal_sign_flipped),
        "normal_resistance_probe_ohm": float(normal_resistance_ohm if normal_resistance_ohm is not None else np.nan),
        "normal_resistance_terminal_ohm": float(
            normal_resistance_terminal_ohm
            if normal_resistance_terminal_ohm is not None
            else np.nan
        ),
    }
    return points, skipped, meta



def select_delta_inset_runs(
    *,
    config_path: str | Path,
    points: Sequence[Mapping[str, Any]],
    requested_currents_uA: Sequence[float],
) -> list[dict[str, Any]]:
    """Resolve exactly four current requests to the nearest available SS runs."""
    requested = [float(v) for v in requested_currents_uA]
    if len(requested) != 4:
        raise ValueError("delta_inset_currents_uA must contain exactly four currents.")

    available = [
        item for item in points
        if str(item.get("run_name", "")) != "synthetic_origin" and np.isfinite(float(item.get("current_uA", np.nan)))
    ]
    if not available:
        return []

    resolved: list[dict[str, Any]] = []
    for idx, req in enumerate(requested, start=1):
        nearest = min(available, key=lambda item: abs(float(item.get("current_uA", np.nan)) - req))
        run_name = str(nearest.get("run_name", ""))
        run = load_ss_run(config_path=config_path, run_name=run_name)
        dataset = build_ss_plot_dataset(run)
        resolved.append(
            {
                "index": int(idx),
                "requested_current_uA": float(req),
                "actual_current_uA": float(nearest.get("current_uA", np.nan)),
                "run_name": run_name,
                "dataset": dataset,
            }
        )
    return resolved


def select_terminal_delta_inset_runs(
    *,
    config_path: str | Path,
    points: Sequence[Mapping[str, Any]],
    requested_currents_uA: Sequence[float],
) -> list[dict[str, Any]]:
    """Resolve exactly three terminal-IV snapshot requests."""
    requested = [float(value) for value in requested_currents_uA]
    if len(requested) != 3:
        raise ValueError("terminal_delta_inset_currents_uA must contain exactly three currents.")
    return _select_nearest_delta_runs(
        config_path=config_path,
        points=points,
        requested_currents_uA=requested,
    )


def _select_nearest_delta_runs(
    *,
    config_path: str | Path,
    points: Sequence[Mapping[str, Any]],
    requested_currents_uA: Sequence[float],
) -> list[dict[str, Any]]:
    available = [
        item
        for item in points
        if str(item.get("run_name", "")) != "synthetic_origin"
        and np.isfinite(float(item.get("current_uA", np.nan)))
    ]
    if not available:
        return []

    resolved: list[dict[str, Any]] = []
    for index, requested in enumerate(requested_currents_uA, start=1):
        nearest = min(
            available,
            key=lambda item: abs(float(item.get("current_uA", np.nan)) - float(requested)),
        )
        run_name = str(nearest.get("run_name", ""))
        run = load_ss_run(config_path=config_path, run_name=run_name)
        resolved.append(
            {
                "index": int(index),
                "requested_current_uA": float(requested),
                "actual_current_uA": float(nearest.get("current_uA", np.nan)),
                "run_name": run_name,
                "dataset": build_ss_plot_dataset(run),
            }
        )
    return resolved



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


__all__ = [
    "collect_current_sweep_iv_points",
    "make_current_sweep_figures",
    "plot_current_sweep_iv",
    "plot_terminal_current_sweep_iv",
    "select_delta_inset_runs",
    "select_terminal_delta_inset_runs",
    "write_current_sweep_iv_csv",
    "write_current_sweep_iv_yaml",
    "write_iv_insets_yaml",
    "write_skipped_runs_yaml",
]
