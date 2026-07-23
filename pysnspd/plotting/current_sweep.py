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
from pysnspd.plotting.style import THESIS_DPI

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

from pysnspd.plotting.current_sweep_iv import (
    _build_iv_point,
    plot_current_sweep_iv,
    plot_terminal_current_sweep_iv,
    write_current_sweep_iv_csv,
    write_current_sweep_iv_yaml,
    write_iv_insets_yaml,
    write_skipped_runs_yaml,
)
from pysnspd.plotting.current_sweep_insets import _orient_positive_voltage
