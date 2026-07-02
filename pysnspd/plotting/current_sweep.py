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
import numpy as np
import yaml

from pysnspd.analysis.ss_run import build_ss_plot_dataset, load_ss_run

MEV_J = 1.602176634e-22


def make_current_sweep_figures(
    *,
    config_path: str | Path,
    records: Sequence[Mapping[str, Any]],
    output_dir: str | Path,
    dpi: int = 480,
    voltage_probe_offset_nm: float = 50.0,
    voltage_probe_half_window_nm: float | None = None,
    include_origin: bool = True,
    delta_inset_currents_uA: Sequence[float] | None = None,
) -> dict[str, Any]:
    """Create current-sweep inventory products and figures."""
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)

    points, skipped, meta = collect_current_sweep_iv_points(
        config_path=config_path,
        records=records,
        voltage_probe_offset_nm=voltage_probe_offset_nm,
        voltage_probe_half_window_nm=voltage_probe_half_window_nm,
        include_origin=include_origin,
    )
    inset_runs = []
    if delta_inset_currents_uA is not None:
        inset_runs = select_delta_inset_runs(
            config_path=config_path,
            points=points,
            requested_currents_uA=delta_inset_currents_uA,
        )

    saved: dict[str, Any] = {}
    saved["iv_curve"] = plot_current_sweep_iv(
        points,
        out / "Z1_iv_curve.png",
        dpi=dpi,
        voltage_probe_offset_nm=voltage_probe_offset_nm,
        include_origin=include_origin,
        delta_insets=inset_runs,
    )
    saved["iv_points_csv"] = write_current_sweep_iv_csv(points, out / "Z1_iv_points.csv")
    saved["iv_points_yaml"] = write_current_sweep_iv_yaml(points, meta, out / "Z1_iv_points.yaml")
    saved["iv_skipped_yaml"] = write_skipped_runs_yaml(skipped, out / "Z1_iv_skipped.yaml")
    saved["iv_insets_yaml"] = write_iv_insets_yaml(inset_runs, out / "Z1_iv_insets.yaml")
    saved["iv_summary"] = {
        "n_points": int(len(points)),
        "n_runs_loaded": int(meta.get("n_runs_loaded", 0)),
        "n_runs_skipped": int(len(skipped)),
        "include_origin": bool(include_origin),
        "voltage_probe_offset_nm": float(voltage_probe_offset_nm),
        "voltage_probe_half_window_nm": float(meta.get("voltage_probe_half_window_nm", np.nan)),
        "voltage_sign_flipped": bool(meta.get("voltage_sign_flipped", False)),
        "delta_inset_currents_uA": [float(v) for v in delta_inset_currents_uA] if delta_inset_currents_uA is not None else [],
        "delta_inset_resolved_currents_uA": [float(item["actual_current_uA"]) for item in inset_runs],
    }
    return saved



def collect_current_sweep_iv_points(
    *,
    config_path: str | Path,
    records: Sequence[Mapping[str, Any]],
    voltage_probe_offset_nm: float = 50.0,
    voltage_probe_half_window_nm: float | None = None,
    include_origin: bool = True,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]], dict[str, Any]]:
    points: list[dict[str, Any]] = []
    skipped: list[dict[str, Any]] = []
    used_half_window_nm: float | None = None
    normal_resistance_ohm: float | None = None

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
            point, half_window_nm, rn_probe = _build_iv_point(
                run_name=run_name,
                run=run,
                dataset=dataset,
                voltage_probe_offset_nm=voltage_probe_offset_nm,
                voltage_probe_half_window_nm=voltage_probe_half_window_nm,
            )
            if used_half_window_nm is None and np.isfinite(half_window_nm):
                used_half_window_nm = float(half_window_nm)
            if normal_resistance_ohm is None and np.isfinite(rn_probe):
                normal_resistance_ohm = float(rn_probe)
            points.append(point)
        except Exception as exc:
            skipped.append({"run_name": run_name, "reason": f"{type(exc).__name__}: {exc}"})

    points.sort(key=lambda item: (float(item.get("current_uA", np.nan)), str(item.get("run_name", ""))))
    sign_flipped = _orient_positive_voltage(points)

    if include_origin:
        origin = {
            "run_name": "synthetic_origin",
            "current_uA": 0.0,
            "voltage_mV": 0.0,
            "normal_voltage_mV": 0.0,
            "normal_resistance_probe_ohm": float(normal_resistance_ohm if normal_resistance_ohm is not None else np.nan),
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
        "normal_resistance_probe_ohm": float(normal_resistance_ohm if normal_resistance_ohm is not None else np.nan),
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



def plot_current_sweep_iv(
    points: Sequence[Mapping[str, Any]],
    output_path: str | Path,
    *,
    dpi: int = 480,
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

    fig, ax = plt.subplots(figsize=(7.8, 4.9), constrained_layout=False)
    fig.subplots_adjust(left=0.105, right=0.975, bottom=0.125, top=0.965)

    normal_line = None
    fit_line = None
    if normal_x.size:
        normal_line, = ax.plot(
            normal_x,
            normal_y,
            linestyle="--",
            linewidth=1.45,
            color="tab:orange",
            label="normal-state Ohmic line",
            zorder=1,
        )
    if fit_x.size:
        fit_line, = ax.plot(
            fit_x,
            fit_y,
            linewidth=1.55,
            color="black",
            label="monotone fit",
            zorder=2,
        )
    raw_scatter = ax.scatter(
        x_valid,
        y_valid,
        s=22.0,
        color="tab:blue",
        label="raw endpoint voltage",
        zorder=3,
    )

    if include_origin:
        ax.scatter([0.0], [0.0], s=26.0, color="tab:orange", zorder=4)

    snapshot_handle = None
    if delta_insets:
        snapshot_handle = _highlight_snapshot_points(ax, points, delta_insets)

    ax.set_xlabel(r"$I_{\mathrm{bias}}$ [$\mu$A]")
    ax.set_ylabel(r"$V_{\mathrm{TDGL}}$ [mV]")
    ax.grid(False)

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

    if delta_insets:
        _add_delta_insets(ax, delta_insets)

    probe_text = (
        rf"probe: $V = \phi(x_c + {voltage_probe_offset_nm:.0f}\,\mathrm{{nm}}) - "
        rf"\phi(x_c - {voltage_probe_offset_nm:.0f}\,\mathrm{{nm}})$"
    )
    handles = [raw_scatter]
    labels = ["raw endpoint voltage"]
    if fit_line is not None:
        handles.append(fit_line)
        labels.append("monotone fit")
    if normal_line is not None:
        handles.append(normal_line)
        labels.append("normal-state Ohmic line")
    if snapshot_handle is not None:
        handles.append(snapshot_handle)
        labels.append("|Δ| snapshot used")

    legend = ax.legend(
        handles,
        labels,
        loc="lower right",
        frameon=True,
        fontsize=8.0,
        title=probe_text,
        title_fontsize=8.0,
    )
    legend.get_frame().set_alpha(0.95)

    fig.savefig(output, dpi=dpi, bbox_inches="tight", pad_inches=0.06)
    plt.close(fig)
    return output


def write_current_sweep_iv_csv(points: Sequence[Mapping[str, Any]], output_path: str | Path) -> Path:
    path = Path(output_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "run_name", "current_uA", "voltage_mV", "normal_voltage_mV", "normal_resistance_probe_ohm",
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
    voltage_probe_offset_nm: float,
    voltage_probe_half_window_nm: float | None,
) -> tuple[dict[str, Any], float, float]:
    current_uA = _infer_bias_current_uA(run_name=run_name, summary=getattr(run, "summary", {}), dataset=dataset)
    left_phi, right_phi, left_x, right_x, half_window_nm = _extract_profile_probe_values(
        dataset,
        voltage_probe_offset_nm=voltage_probe_offset_nm,
        voltage_probe_half_window_nm=voltage_probe_half_window_nm,
    )
    x_profile = np.asarray(dataset.get("x_profile_nm", []), dtype=float)
    rn_probe_ohm = _infer_probe_normal_resistance_ohm(run=run, dataset=dataset, voltage_probe_offset_nm=voltage_probe_offset_nm)
    normal_voltage_mV = 1.0e-3 * float(current_uA) * float(rn_probe_ohm) if np.isfinite(rn_probe_ohm) else np.nan

    point = {
        "run_name": run_name,
        "current_uA": float(current_uA),
        "voltage_mV": float(right_phi - left_phi),
        "normal_voltage_mV": float(normal_voltage_mV),
        "normal_resistance_probe_ohm": float(rn_probe_ohm),
        "probe_left_x_nm": float(left_x),
        "probe_right_x_nm": float(right_x),
        "probe_left_phi_mV": float(left_phi),
        "probe_right_phi_mV": float(right_phi),
        "profile_x_center_nm": float(0.5 * (np.nanmin(x_profile) + np.nanmax(x_profile))) if x_profile.size else np.nan,
        "profile_x_min_nm": float(np.nanmin(x_profile)) if x_profile.size else np.nan,
        "profile_x_max_nm": float(np.nanmax(x_profile)) if x_profile.size else np.nan,
        "voltage_probe_offset_nm": float(voltage_probe_offset_nm),
        "voltage_probe_half_window_nm": float(half_window_nm),
        "pre_run_name": getattr(run, "pre_run_name", None),
        "raw_ss": str(getattr(run, "raw_ss", "")) if getattr(run, "raw_ss", None) is not None else None,
        "source": "x_profile_phi_mV",
    }
    return point, float(half_window_nm), float(rn_probe_ohm)



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



def _infer_probe_normal_resistance_ohm(*, run: Any, dataset: Mapping[str, Any], voltage_probe_offset_nm: float) -> float:
    width_m = float(getattr(run.mesh, "width_m", np.nan))
    thickness_m = _find_numeric_recursive(
        getattr(run, "summary", {}),
        keys=("thickness_m", "thickness", "d_m", "film_thickness_m"),
    )
    sigma_n = _find_numeric_recursive(
        getattr(run, "summary", {}),
        keys=("sigma_n_S_m", "sigma_n_S_per_m", "sigma_n", "normal_conductivity_S_m"),
    )
    if not np.isfinite(width_m) or not np.isfinite(thickness_m) or not np.isfinite(sigma_n):
        return np.nan
    if width_m <= 0 or thickness_m <= 0 or sigma_n <= 0:
        return np.nan
    length_m = 2.0 * float(voltage_probe_offset_nm) * 1.0e-9
    return float(length_m / (float(sigma_n) * width_m * float(thickness_m)))



def _extract_iv_arrays(points: Sequence[Mapping[str, Any]]) -> tuple[np.ndarray, np.ndarray]:
    current_uA = np.asarray([float(item.get("current_uA", np.nan)) for item in points], dtype=float)
    voltage_mV = np.asarray([float(item.get("voltage_mV", np.nan)) for item in points], dtype=float)
    return current_uA, voltage_mV



def _extract_normal_curve(points: Sequence[Mapping[str, Any]]) -> tuple[np.ndarray, np.ndarray]:
    current_uA = np.asarray([float(item.get("current_uA", np.nan)) for item in points], dtype=float)
    normal_mV = np.asarray([float(item.get("normal_voltage_mV", np.nan)) for item in points], dtype=float)
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
    """Place four |Delta| colormap insets inside the main IV axes in a 2x2 grid."""
    positions = [
        [0.055, 0.58, 0.16, 0.29],
        [0.235, 0.58, 0.16, 0.29],
        [0.055, 0.265, 0.16, 0.29],
        [0.235, 0.265, 0.16, 0.29],
    ]
    for pos, inset in zip(positions, delta_insets):
        ax_in = ax.inset_axes(pos)
        _draw_delta_inset(ax_in, inset)



def _highlight_snapshot_points(
    ax: plt.Axes,
    points: Sequence[Mapping[str, Any]],
    delta_insets: Sequence[Mapping[str, Any]],
):
    xs = []
    ys = []
    indices = []
    for inset in delta_insets:
        x = float(inset.get("actual_current_uA", np.nan))
        y = _lookup_voltage(points, x)
        if np.isfinite(x) and np.isfinite(y):
            xs.append(x)
            ys.append(y)
            indices.append(int(inset.get("index", 0)))
    if not xs:
        return None

    handle = ax.scatter(
        xs,
        ys,
        s=72.0,
        facecolors="red",
        edgecolors="black",
        linewidths=0.9,
        zorder=5,
    )
    for x, y, idx in zip(xs, ys, indices):
        ax.text(
            x,
            y,
            str(idx),
            ha="center",
            va="center",
            fontsize=7.0,
            color="white",
            fontweight="bold",
            zorder=6,
        )
    return handle



def _draw_delta_inset(ax: plt.Axes, inset: Mapping[str, Any]) -> None:
    dataset = inset.get("dataset", {})
    x_nm = np.asarray(dataset.get("x_nm", []), dtype=float)
    y_nm = np.asarray(dataset.get("y_nm", []), dtype=float)
    triangles = np.asarray(dataset.get("triangles", []), dtype=np.int64)
    delta = np.asarray(dataset.get("delta_over_delta0", []), dtype=float)
    if x_nm.size == 0 or y_nm.size == 0 or triangles.size == 0 or delta.size != x_nm.size:
        ax.text(0.5, 0.5, "missing\n|Δ| data", ha="center", va="center", transform=ax.transAxes, fontsize=7.0)
        ax.set_xticks([])
        ax.set_yticks([])
        return

    x_center = 0.5 * (float(np.nanmin(x_nm)) + float(np.nanmax(x_nm)))
    x_left = x_center - 50.0
    x_right = x_center + 50.0
    y_min = float(np.nanmin(y_nm))
    y_max = float(np.nanmax(y_nm))

    triang = mtri.Triangulation(x_nm, y_nm, triangles)
    ax.tripcolor(triang, delta, shading="gouraud", vmin=0.0, vmax=1.0)
    ax.set_xlim(x_left, x_right)
    ax.set_ylim(y_min, y_max)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xticks([])
    ax.set_yticks([])

    idx = int(inset.get("index", 0))
    current = float(inset.get("actual_current_uA", np.nan))
    requested = float(inset.get("requested_current_uA", np.nan))
    if np.isfinite(requested) and abs(requested - current) > 0.05:
        label = f"#{idx}  {current:.0f} μA\n(req {requested:.0f})"
    else:
        label = f"#{idx}  {current:.0f} μA"
    ax.text(
        0.965,
        0.965,
        label,
        ha="right",
        va="top",
        transform=ax.transAxes,
        fontsize=7.0,
        color="white",
        bbox={"boxstyle": "round,pad=0.16", "facecolor": "black", "edgecolor": "white", "linewidth": 0.45, "alpha": 0.72},
        zorder=10,
    )
    for spine in ax.spines.values():
        spine.set_linewidth(0.7)
        spine.set_edgecolor("0.1")



def _lookup_voltage(points: Sequence[Mapping[str, Any]], current_uA: float) -> float:
    for item in points:
        if np.isfinite(float(item.get("current_uA", np.nan))) and abs(float(item.get("current_uA")) - current_uA) < 1.0e-9:
            return float(item.get("voltage_mV", np.nan))
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



def _orient_positive_voltage(points: Sequence[dict[str, Any]]) -> bool:
    voltages = np.asarray([float(item.get("voltage_mV", np.nan)) for item in points], dtype=float)
    currents = np.asarray([float(item.get("current_uA", np.nan)) for item in points], dtype=float)
    mask = np.isfinite(voltages) & np.isfinite(currents) & (currents > 0.0) & (np.abs(voltages) > 0.0)
    if not np.any(mask):
        return False
    median_v = float(np.nanmedian(voltages[mask]))
    if median_v >= 0.0:
        return False
    for item in points:
        try:
            item["voltage_mV"] = -float(item.get("voltage_mV", np.nan))
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
    "select_delta_inset_runs",
    "write_current_sweep_iv_csv",
    "write_current_sweep_iv_yaml",
    "write_iv_insets_yaml",
    "write_skipped_runs_yaml",
]
