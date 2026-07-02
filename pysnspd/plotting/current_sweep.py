"""Current-sweep plotting helpers for Z-series multi-run analysis."""

from __future__ import annotations

import csv
import re
from pathlib import Path
from typing import Any, Mapping, Sequence

import matplotlib

matplotlib.use("Agg", force=True)
import matplotlib.pyplot as plt
import numpy as np
import yaml

from pysnspd.analysis.ss_run import build_ss_plot_dataset, load_ss_run


def make_current_sweep_figures(
    *,
    config_path: str | Path,
    records: Sequence[Mapping[str, Any]],
    output_dir: str | Path,
    dpi: int = 480,
    voltage_probe_offset_nm: float = 50.0,
    voltage_probe_half_window_nm: float | None = None,
    include_origin: bool = True,
) -> dict[str, Any]:
    """Create first-pass current-sweep figures and tabulated IV data."""
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)

    points, skipped, meta = collect_current_sweep_iv_points(
        config_path=config_path,
        records=records,
        voltage_probe_offset_nm=voltage_probe_offset_nm,
        voltage_probe_half_window_nm=voltage_probe_half_window_nm,
        include_origin=include_origin,
    )

    saved: dict[str, Any] = {}
    saved["iv_curve"] = plot_current_sweep_iv(
        points,
        out / "Z1_iv_curve.png",
        dpi=dpi,
        voltage_probe_offset_nm=voltage_probe_offset_nm,
        include_origin=include_origin,
    )
    saved["iv_points_csv"] = write_current_sweep_iv_csv(points, out / "Z1_iv_points.csv")
    saved["iv_points_yaml"] = write_current_sweep_iv_yaml(
        points,
        meta,
        out / "Z1_iv_points.yaml",
    )
    saved["iv_skipped_yaml"] = write_skipped_runs_yaml(skipped, out / "Z1_iv_skipped.yaml")
    saved["iv_summary"] = {
        "n_points": int(len(points)),
        "n_runs_loaded": int(meta.get("n_runs_loaded", 0)),
        "n_runs_skipped": int(len(skipped)),
        "include_origin": bool(include_origin),
        "voltage_probe_offset_nm": float(voltage_probe_offset_nm),
        "voltage_probe_half_window_nm": float(meta.get("voltage_probe_half_window_nm", np.nan)),
        "voltage_sign_flipped": bool(meta.get("voltage_sign_flipped", False)),
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
    """Load SS runs and extract IV points from x-profile voltages."""
    points: list[dict[str, Any]] = []
    skipped: list[dict[str, Any]] = []

    used_half_window_nm: float | None = None
    for record in records:
        run_name = str(record.get("run_name", ""))
        stages = record.get("stages", {})
        stage_ss = stages.get("ss", {}) if isinstance(stages, Mapping) else {}
        if not isinstance(stage_ss, Mapping) or not stage_ss.get("exists", False):
            skipped.append(
                {
                    "run_name": run_name,
                    "reason": "ss stage not found",
                }
            )
            continue
        try:
            run = load_ss_run(config_path=config_path, run_name=run_name)
            dataset = build_ss_plot_dataset(run)
            point, half_window_nm = _build_iv_point(
                run_name=run_name,
                run=run,
                dataset=dataset,
                voltage_probe_offset_nm=voltage_probe_offset_nm,
                voltage_probe_half_window_nm=voltage_probe_half_window_nm,
            )
            if used_half_window_nm is None and np.isfinite(half_window_nm):
                used_half_window_nm = float(half_window_nm)
            points.append(point)
        except Exception as exc:
            skipped.append(
                {
                    "run_name": run_name,
                    "reason": f"{type(exc).__name__}: {exc}",
                }
            )

    points.sort(key=lambda item: (float(item.get("current_uA", np.nan)), str(item.get("run_name", ""))))
    sign_flipped = _orient_positive_voltage(points)

    if include_origin:
        origin = {
            "run_name": "synthetic_origin",
            "current_uA": 0.0,
            "voltage_mV": 0.0,
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
    }
    return points, skipped, meta



def plot_current_sweep_iv(
    points: Sequence[Mapping[str, Any]],
    output_path: str | Path,
    *,
    dpi: int = 480,
    voltage_probe_offset_nm: float = 50.0,
    include_origin: bool = True,
) -> Path:
    """Plot current on x-axis and central TDGL voltage on y-axis."""
    output = Path(output_path)
    output.parent.mkdir(parents=True, exist_ok=True)

    current_uA = np.asarray([float(item.get("current_uA", np.nan)) for item in points], dtype=float)
    voltage_mV = np.asarray([float(item.get("voltage_mV", np.nan)) for item in points], dtype=float)
    valid = np.isfinite(current_uA) & np.isfinite(voltage_mV)

    fig, ax = plt.subplots(figsize=(7.0, 4.8), constrained_layout=False)
    fig.subplots_adjust(left=0.120, right=0.970, bottom=0.120, top=0.935)

    if np.any(valid):
        ax.plot(current_uA[valid], voltage_mV[valid], linewidth=1.1, alpha=0.90, zorder=2)
        ax.scatter(current_uA[valid], voltage_mV[valid], s=16.0, alpha=0.95, zorder=3)

    if include_origin:
        ax.scatter([0.0], [0.0], s=22.0, zorder=4)

    ax.set_xlabel(r"$I_{\mathrm{bias}}$ [$\mu$A]")
    ax.set_ylabel(r"$V_{\mathrm{TDGL}}$ [mV]")
    ax.grid(False)

    finite_v = voltage_mV[valid]
    finite_i = current_uA[valid]
    if finite_i.size:
        xmin = float(np.nanmin(finite_i))
        xmax = float(np.nanmax(finite_i))
        dx = max(xmax - xmin, 1.0)
        ax.set_xlim(min(0.0, xmin) - 0.03 * dx, xmax + 0.04 * dx)
    if finite_v.size:
        ymin = float(np.nanmin(finite_v))
        ymax = float(np.nanmax(finite_v))
        dy = max(ymax - ymin, 1.0e-6)
        lower = min(0.0, ymin - 0.06 * dy)
        upper = max(0.0, ymax + 0.08 * dy)
        if upper <= lower:
            upper = lower + 1.0
        ax.set_ylim(lower, upper)

    probe_text = (
        rf"probe: $V = \phi(x_c+{voltage_probe_offset_nm:.0f}\,\mathrm{{nm}}) - "
        rf"\phi(x_c-{voltage_probe_offset_nm:.0f}\,\mathrm{{nm}})$"
    )
    ax.text(
        0.02,
        0.98,
        probe_text,
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=8.0,
        bbox={
            "boxstyle": "round,pad=0.26",
            "facecolor": "white",
            "edgecolor": "0.35",
            "linewidth": 0.55,
            "alpha": 0.88,
        },
        zorder=10,
    )

    fig.savefig(output, dpi=dpi, bbox_inches="tight", pad_inches=0.08)
    plt.close(fig)
    return output



def write_current_sweep_iv_csv(points: Sequence[Mapping[str, Any]], output_path: str | Path) -> Path:
    path = Path(output_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "run_name",
        "current_uA",
        "voltage_mV",
        "probe_left_x_nm",
        "probe_right_x_nm",
        "probe_left_phi_mV",
        "probe_right_phi_mV",
        "profile_x_center_nm",
        "profile_x_min_nm",
        "profile_x_max_nm",
        "voltage_probe_offset_nm",
        "voltage_probe_half_window_nm",
        "pre_run_name",
        "raw_ss",
        "source",
    ]
    with path.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for item in points:
            writer.writerow({key: item.get(key, "") for key in fieldnames})
    return path



def write_current_sweep_iv_yaml(
    points: Sequence[Mapping[str, Any]],
    meta: Mapping[str, Any],
    output_path: str | Path,
) -> Path:
    path = Path(output_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        "meta": dict(meta),
        "points": [_to_builtin(item) for item in points],
    }
    with path.open("w", encoding="utf-8") as f:
        yaml.safe_dump(payload, f, sort_keys=False, allow_unicode=True, default_flow_style=False)
    return path



def write_skipped_runs_yaml(skipped: Sequence[Mapping[str, Any]], output_path: str | Path) -> Path:
    path = Path(output_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as f:
        yaml.safe_dump([_to_builtin(item) for item in skipped], f, sort_keys=False, allow_unicode=True, default_flow_style=False)
    return path



def _build_iv_point(
    *,
    run_name: str,
    run: Any,
    dataset: Mapping[str, Any],
    voltage_probe_offset_nm: float,
    voltage_probe_half_window_nm: float | None,
) -> tuple[dict[str, Any], float]:
    current_uA = _infer_bias_current_uA(run_name=run_name, summary=getattr(run, "summary", {}), dataset=dataset)
    left_phi, right_phi, left_x, right_x, half_window_nm = _extract_profile_probe_values(
        dataset,
        voltage_probe_offset_nm=voltage_probe_offset_nm,
        voltage_probe_half_window_nm=voltage_probe_half_window_nm,
    )
    x_profile = np.asarray(dataset.get("x_profile_nm", []), dtype=float)

    point = {
        "run_name": run_name,
        "current_uA": float(current_uA),
        "voltage_mV": float(right_phi - left_phi),
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
    return point, float(half_window_nm)



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
            "target_current_A",
            "current_A",
            "I_bias_A",
            "bias_current_A",
            "target_current_uA",
            "current_uA",
            "I_bias_uA",
            "bias_current_uA",
        ),
    )
    if scalar is not None:
        key, value = scalar
        if key.endswith("_uA"):
            return float(value)
        return 1.0e6 * float(value)

    summary_scalars = dataset.get("summary_scalars", {})
    if isinstance(summary_scalars, Mapping):
        if "target_current_A" in summary_scalars:
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
    "write_current_sweep_iv_csv",
    "write_current_sweep_iv_yaml",
    "write_skipped_runs_yaml",
]
