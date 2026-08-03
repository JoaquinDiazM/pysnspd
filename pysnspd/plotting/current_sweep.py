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
from pysnspd.plotting.current_sweep_summary import (
    build_current_sweep_regime_summary,
    plot_current_sweep_regime_summary,
    write_current_sweep_regime_summary,
)

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
    ohmic_relative_tolerance: float = 0.10,
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
        ohmic_relative_tolerance=ohmic_relative_tolerance,
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
    regime_summary = build_current_sweep_regime_summary(
        points,
        skipped,
        ohmic_relative_tolerance=ohmic_relative_tolerance,
    )
    saved["regime_summary_curve"] = plot_current_sweep_regime_summary(
        points,
        skipped,
        out / "Z2_sweep_regime_summary.pdf",
        dpi=dpi,
        ohmic_relative_tolerance=ohmic_relative_tolerance,
    )
    saved["regime_summary_yaml"] = write_current_sweep_regime_summary(
        regime_summary,
        out / "Z2_sweep_regime_summary.yaml",
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
        "ohmic_relative_tolerance": float(ohmic_relative_tolerance),
        "regime_counts": dict(regime_summary.get("counts", {})),
        "sampled_ranges": dict(regime_summary.get("sampled_ranges", {})),
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
    ohmic_relative_tolerance: float = 0.10,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]], dict[str, Any]]:
    points: list[dict[str, Any]] = []
    skipped: list[dict[str, Any]] = []
    used_half_window_nm: float | None = None
    normal_resistance_ohm: float | None = None
    normal_resistance_terminal_ohm: float | None = None

    if not np.isfinite(ohmic_relative_tolerance) or not 0.0 < ohmic_relative_tolerance < 1.0:
        raise ValueError("ohmic_relative_tolerance must lie strictly between 0 and 1.")

    for record in records:
        run_name = str(record.get("run_name", ""))
        stages = record.get("stages", {})
        stage_ss = stages.get("ss", {}) if isinstance(stages, Mapping) else {}
        if not isinstance(stage_ss, Mapping) or not stage_ss.get("exists", False):
            skipped.append(_incomplete_case(run_name, "ss stage not found"))
            continue
        missing = _missing_required_endpoint_outputs(stage_ss)
        if missing:
            skipped.append(
                _incomplete_case(
                    run_name,
                    "missing required endpoint outputs: " + ", ".join(missing),
                )
            )
            continue
        try:
            run = load_ss_run(
                config_path=config_path,
                run_name=run_name,
                load_history=False,
            )
            dataset = build_ss_plot_dataset(run, load_snapshots=False)
            point, half_window_nm, rn_probe, rn_terminal = _build_iv_point(
                run_name=run_name,
                run=run,
                dataset=dataset,
                project_config=project_config,
                voltage_probe_offset_nm=voltage_probe_offset_nm,
                voltage_probe_half_window_nm=voltage_probe_half_window_nm,
            )
            point.update(_summary_diagnostics(run.summary))
            if point.get("complete") is not True:
                skipped.append(
                    {
                        **_incomplete_case(run_name, "requested simulation horizon was not reached"),
                        **point,
                    }
                )
                continue
            if used_half_window_nm is None and np.isfinite(half_window_nm):
                used_half_window_nm = float(half_window_nm)
            if normal_resistance_ohm is None and np.isfinite(rn_probe):
                normal_resistance_ohm = float(rn_probe)
            if normal_resistance_terminal_ohm is None and np.isfinite(rn_terminal):
                normal_resistance_terminal_ohm = float(rn_terminal)
            points.append(point)
        except Exception as exc:
            skipped.append(
                _incomplete_case(run_name, f"{type(exc).__name__}: {exc}")
            )

    points.sort(key=lambda item: (float(item.get("current_uA", np.nan)), str(item.get("run_name", ""))))
    sign_flipped = _orient_positive_voltage(points)
    terminal_sign_flipped = _orient_positive_voltage(
        points,
        voltage_key="terminal_voltage_mV",
    )
    _annotate_normal_voltage_ratios(
        points,
        ohmic_relative_tolerance=float(ohmic_relative_tolerance),
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
        "ohmic_relative_tolerance": float(ohmic_relative_tolerance),
    }
    return points, skipped, meta


def _missing_required_endpoint_outputs(stage_ss: Mapping[str, Any]) -> list[str]:
    npz_names = {
        Path(str(item.get("path") or item.get("relative_path") or "")).name
        for item in stage_ss.get("npz_files", []) or []
        if isinstance(item, Mapping)
    }
    summary_names = {
        Path(str(item.get("path") or item.get("relative_path") or "")).name
        for item in stage_ss.get("summary_files", []) or []
        if isinstance(item, Mapping)
    }
    missing: list[str] = []
    if "ss_summary.yaml" not in summary_names:
        missing.append("ss_summary.yaml")
    if "stationary_state.npz" not in npz_names:
        missing.append("stationary_state.npz")
    return missing


def _incomplete_case(run_name: str, reason: str) -> dict[str, Any]:
    return {
        "run_name": str(run_name),
        "current_uA": _current_from_run_name(run_name),
        "complete": False,
        "reason": str(reason),
        "strict_stationarity_passes": None,
        "dynamic_stationarity_passes": None,
        "photon_ready": None,
        "approximately_ohmic": None,
    }


def _current_from_run_name(run_name: str) -> float:
    matches = re.findall(r"(?:^|_)I([-+]?\d+(?:\.\d+)?)uA(?:_|$)", str(run_name))
    if matches:
        return float(matches[-1])
    sweep = re.search(r"(?:^|_)I([-+]?\d+(?:\.\d+)?)to[-+]?\d+(?:\.\d+)?uA", str(run_name))
    return float(sweep.group(1)) if sweep else float("nan")


def _summary_diagnostics(summary: Mapping[str, Any]) -> dict[str, Any]:
    solver = summary.get("solver", {})
    if not isinstance(solver, Mapping):
        solver = {}
    strict = _nested_optional_bool(solver, "stationarity", "passes")
    dynamic = _first_optional_bool(
        solver.get("dynamic_stationarity_passes"),
        _nested_value(solver, "dynamic_stationarity", "passes"),
    )
    contact = _nested_optional_bool(solver, "contact_recovery", "passes")
    continuity = _nested_optional_bool(solver, "continuity", "passes")
    thermal = _nested_optional_bool(solver, "thermal_stationarity", "passes")
    circuit = _nested_optional_bool(solver, "circuit_stationarity", "passes")
    phase = _first_optional_bool(
        solver.get("phase_drive_converged"),
        _nested_value(solver, "allmaras_phase_continuation", "final_converged"),
    )
    stored_ready = _optional_bool(solver.get("photon_ready"))
    if stored_ready is None:
        mesoscopic = _logical_or_optional(strict, dynamic)
        final_ready = _logical_and_optional(
            mesoscopic,
            contact,
            continuity,
            thermal,
            circuit,
            phase,
        )
        ready_source = "reconstructed_final_gate" if final_ready is not None else "unavailable"
    else:
        final_ready = stored_ready
        ready_source = "stored"

    requested_time = _finite_or_nan(solver.get("requested_time_ps"))
    final_time = _finite_or_nan(solver.get("final_time_ps"))
    reached = _optional_bool(solver.get("requested_time_reached"))
    if reached is None and np.isfinite(requested_time) and requested_time > 0.0 and np.isfinite(final_time):
        reached = bool(final_time >= 0.999999 * requested_time)
    accepted = _finite_or_nan(solver.get("accepted_steps"))
    rejected = _finite_or_nan(solver.get("rejected_steps"))
    rejected_ratio = (
        rejected / accepted
        if np.isfinite(accepted) and accepted > 0.0 and np.isfinite(rejected)
        else np.nan
    )
    return {
        "complete": bool(reached) if reached is not None else False,
        "strict_stationarity_passes": strict,
        "dynamic_stationarity_passes": dynamic,
        "photon_ready": final_ready,
        "photon_ready_source": ready_source,
        "contact_recovery_passes": contact,
        "continuity_passes": continuity,
        "thermal_stationarity_passes": thermal,
        "circuit_stationarity_passes": circuit,
        "phase_drive_converged": phase,
        "mean_delta_over_delta0": _finite_or_nan(solver.get("mean_delta_over_delta0")),
        "normal_like_fraction_final": _finite_or_nan(
            _nested_value(solver, "dynamic_stationarity", "normal_like_fraction_final")
        ),
        "accepted_steps": accepted,
        "rejected_steps": rejected,
        "rejected_over_accepted": rejected_ratio,
        "final_time_ps": final_time,
    }


def _annotate_normal_voltage_ratios(
    points: Sequence[dict[str, Any]],
    *,
    ohmic_relative_tolerance: float,
) -> None:
    for point in points:
        central = _safe_ratio(point.get("voltage_mV"), point.get("normal_voltage_mV"))
        terminal = _safe_ratio(
            point.get("terminal_voltage_mV"),
            point.get("normal_terminal_voltage_mV"),
        )
        point["normal_voltage_ratio"] = central
        point["terminal_normal_voltage_ratio"] = terminal
        point["approximately_ohmic"] = (
            bool(
                abs(central - 1.0) <= ohmic_relative_tolerance
                and abs(terminal - 1.0) <= ohmic_relative_tolerance
            )
            if np.isfinite(central) and np.isfinite(terminal)
            else None
        )


def _safe_ratio(value: Any, reference: Any) -> float:
    numerator = _finite_or_nan(value)
    denominator = _finite_or_nan(reference)
    if not np.isfinite(numerator) or not np.isfinite(denominator) or abs(denominator) <= 0.0:
        return float("nan")
    return float(abs(numerator) / abs(denominator))


def _finite_or_nan(value: Any) -> float:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return float("nan")
    return number if np.isfinite(number) else float("nan")


def _nested_value(mapping: Mapping[str, Any], *keys: str) -> Any:
    value: Any = mapping
    for key in keys:
        if not isinstance(value, Mapping) or key not in value:
            return None
        value = value[key]
    return value


def _optional_bool(value: Any) -> bool | None:
    return bool(value) if isinstance(value, (bool, np.bool_)) else None


def _nested_optional_bool(mapping: Mapping[str, Any], *keys: str) -> bool | None:
    return _optional_bool(_nested_value(mapping, *keys))


def _first_optional_bool(*values: Any) -> bool | None:
    for value in values:
        parsed = _optional_bool(value)
        if parsed is not None:
            return parsed
    return None


def _logical_or_optional(left: bool | None, right: bool | None) -> bool | None:
    if left is True or right is True:
        return True
    if left is False and right is False:
        return False
    return None


def _logical_and_optional(*values: bool | None) -> bool | None:
    if any(value is False for value in values):
        return False
    if all(value is True for value in values):
        return True
    return None



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
        run = load_ss_run(
            config_path=config_path,
            run_name=run_name,
            load_history=False,
        )
        dataset = build_ss_plot_dataset(run, load_snapshots=False)
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
        run = load_ss_run(
            config_path=config_path,
            run_name=run_name,
            load_history=False,
        )
        resolved.append(
            {
                "index": int(index),
                "requested_current_uA": float(requested),
                "actual_current_uA": float(nearest.get("current_uA", np.nan)),
                "run_name": run_name,
                "dataset": build_ss_plot_dataset(run, load_snapshots=False),
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
