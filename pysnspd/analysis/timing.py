"""Detection-latency and recovery metrics for photon transients.

The functions in this module are deliberately independent from the solver.
They can be called incrementally while a transient is running or later from
``transient_history.npz`` and ``transient_snapshots.npz``.  Times are stored in
SI internally and reported in picoseconds.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
from typing import Any, Mapping

import numpy as np


@dataclass(frozen=True)
class DetectionCriteria:
    """Operational leading-edge detection criterion."""

    threshold_V: float = 100.0e-6
    polarity: str = "positive"
    baseline_window_s: float = 10.0e-12
    confirmation_s: float = 0.5e-12
    hysteresis_fraction: float = 0.10
    peak_confirmation_s: float = 2.0e-12
    peak_drop_fraction: float = 0.01
    post_peak_safety_s: float = 10.0e-12

    def validated(self) -> "DetectionCriteria":
        if not np.isfinite(self.threshold_V) or self.threshold_V <= 0.0:
            raise ValueError("Detection threshold must be positive and finite.")
        if self.polarity not in {"positive", "negative", "auto"}:
            raise ValueError("Detection polarity must be positive, negative, or auto.")
        for name in (
            "baseline_window_s",
            "confirmation_s",
            "peak_confirmation_s",
            "post_peak_safety_s",
        ):
            if not np.isfinite(getattr(self, name)) or getattr(self, name) < 0.0:
                raise ValueError(f"{name} must be finite and non-negative.")
        if not 0.0 <= self.hysteresis_fraction < 1.0:
            raise ValueError("hysteresis_fraction must lie in [0, 1).")
        if not 0.0 <= self.peak_drop_fraction < 1.0:
            raise ValueError("peak_drop_fraction must lie in [0, 1).")
        return self


@dataclass(frozen=True)
class RecoveryCriteria:
    """Recovery tolerances shared by online and post-processed analysis."""

    mode: str = "electrical"
    hold_s: float = 10.0e-12
    efficiency_fraction: float = 0.90
    current_relative_tolerance: float = 1.0e-2
    current_absolute_tolerance_A: float = 0.05e-6
    voltage_relative_tolerance: float = 1.0e-2
    voltage_absolute_tolerance_V: float = 10.0e-6
    temperature_absolute_tolerance_K: float = 0.05
    condensate_relative_tolerance: float = 2.0e-2
    spatial_quantile: float = 0.995
    spatial_max_guard_factor: float = 4.0

    def validated(self) -> "RecoveryCriteria":
        if self.mode not in {"electrical", "efficiency90", "state"}:
            raise ValueError("Recovery mode must be electrical, efficiency90, or state.")
        if not np.isfinite(self.hold_s) or self.hold_s <= 0.0:
            raise ValueError("Recovery hold time must be positive and finite.")
        if not 0.0 < self.efficiency_fraction <= 1.0:
            raise ValueError("efficiency_fraction must lie in (0, 1].")
        for name in (
            "current_relative_tolerance",
            "current_absolute_tolerance_A",
            "voltage_relative_tolerance",
            "voltage_absolute_tolerance_V",
            "temperature_absolute_tolerance_K",
            "condensate_relative_tolerance",
        ):
            if not np.isfinite(getattr(self, name)) or getattr(self, name) < 0.0:
                raise ValueError(f"{name} must be finite and non-negative.")
        if not 0.5 <= self.spatial_quantile <= 1.0:
            raise ValueError("spatial_quantile must lie in [0.5, 1].")
        if not np.isfinite(self.spatial_max_guard_factor) or self.spatial_max_guard_factor < 1.0:
            raise ValueError("spatial_max_guard_factor must be finite and at least one.")
        return self


def analyze_photon_timing(
    history: Mapping[str, Any],
    *,
    snapshots: Mapping[str, Any] | None = None,
    detection: DetectionCriteria | None = None,
    recovery: RecoveryCriteria | None = None,
) -> dict[str, Any]:
    """Return latency, recovery, and early-stop readiness from saved data."""

    det = (detection or DetectionCriteria()).validated()
    rec = (recovery or RecoveryCriteria()).validated()
    t_s = _history_time_s(history)
    n = t_s.size
    if n == 0:
        return _unavailable_result(det, rec, "history has no time samples")

    event_index = _event_index(history, n)
    if event_index is None:
        return _unavailable_result(det, rec, "photon event is not present in history")
    event_s = float(t_s[event_index])
    baseline_mask = (t_s >= event_s - det.baseline_window_s) & (t_s < event_s)
    if not np.any(baseline_mask):
        baseline_mask = np.arange(n) < event_index
    if not np.any(baseline_mask):
        return _unavailable_result(det, rec, "history has no pre-photon baseline")

    baselines = {
        key: _robust_baseline(history, key, baseline_mask, n)
        for key in (
            "I_b_A",
            "I_s_A",
            "I_rf_A",
            "V_out_V",
            "v_c_V",
            "V_tdgl_center_V",
            "min_delta_over_delta0",
            "mean_delta_over_delta0",
            "max_Te_K",
            "max_Tph_K",
        )
    }
    voltage = _series(history, "V_out_V", n)
    polarity = _resolve_polarity(
        det.polarity,
        voltage[event_index:] - baselines["V_out_V"],
    )
    signal = polarity * (voltage - baselines["V_out_V"])
    detection_result = _detection_result(
        t_s,
        signal,
        event_index=event_index,
        event_s=event_s,
        criteria=det,
        polarity=polarity,
    )

    start_index = int(detection_result.get("peak_index", event_index))
    detected = bool(detection_result["detected"])
    electrical_mask, electrical_components = _electrical_recovery_mask(
        history,
        n=n,
        baselines=baselines,
        criteria=rec,
    )
    electrical = _recovery_result(
        t_s,
        electrical_mask,
        event_s=event_s,
        start_index=start_index,
        hold_s=rec.hold_s,
        enabled=detected,
        label="electrical",
        details=electrical_components,
    )

    efficiency_mask = _efficiency_recovery_mask(
        history,
        n=n,
        baseline_current_A=baselines["I_s_A"],
        fraction=rec.efficiency_fraction,
        start_index=start_index,
    )
    efficiency = _recovery_result(
        t_s,
        efficiency_mask,
        event_s=event_s,
        start_index=start_index,
        hold_s=rec.hold_s,
        enabled=detected,
        label="efficiency90",
        details={
            "definition": "detector-current proxy for relative detection-efficiency recovery",
            "target_fraction": float(rec.efficiency_fraction),
        },
    )

    state = _state_recovery_result(
        history,
        snapshots=snapshots,
        history_t_s=t_s,
        baselines=baselines,
        electrical_mask=electrical_mask,
        event_s=event_s,
        event_index=event_index,
        start_index=start_index,
        enabled=detected,
        criteria=rec,
    )
    recovery_results = {
        "electrical": electrical,
        "efficiency90": efficiency,
        "state": state,
    }
    selected = dict(recovery_results[rec.mode])
    selected["mode"] = rec.mode

    latency_stop_s = detection_result.get("latency_stop_ready_time_s")
    recovery_stop_s = selected.get("confirmed_time_s")
    return {
        "schema_version": 1,
        "event": {
            "photon_time_s": event_s,
            "photon_time_ps": event_s / 1.0e-12,
            "history_index": int(event_index),
        },
        "baseline": {
            "window_s": float(det.baseline_window_s),
            "window_ps": float(det.baseline_window_s / 1.0e-12),
            "sample_count": int(np.count_nonzero(baseline_mask)),
            "values": {key: float(value) for key, value in baselines.items()},
        },
        "detection_criteria": asdict(det),
        "recovery_criteria": asdict(rec),
        "latency": detection_result,
        "recovery": {
            "selected": selected,
            "electrical": electrical,
            "efficiency90": efficiency,
            "state": state,
        },
        "termination": {
            "latency_ready": bool(
                latency_stop_s is not None
                and np.isfinite(float(latency_stop_s))
                and t_s[-1] >= float(latency_stop_s)
            ),
            "latency_ready_time_s": _finite_or_none(latency_stop_s),
            "latency_ready_time_ps": _ps_or_none(latency_stop_s),
            "recovery_ready": bool(selected.get("recovered", False)),
            "recovery_ready_time_s": _finite_or_none(recovery_stop_s),
            "recovery_ready_time_ps": _ps_or_none(recovery_stop_s),
        },
    }


def _detection_result(
    t_s: np.ndarray,
    signal: np.ndarray,
    *,
    event_index: int,
    event_s: float,
    criteria: DetectionCriteria,
    polarity: float,
) -> dict[str, Any]:
    post_t = t_s[event_index:]
    post_y = signal[event_index:]
    crossing_s = _confirmed_crossing_time(
        post_t,
        post_y,
        threshold=float(criteria.threshold_V),
        confirmation_s=float(criteria.confirmation_s),
        hysteresis_fraction=float(criteria.hysteresis_fraction),
    )
    detected = crossing_s is not None
    finite_signal = np.where(np.isfinite(post_y), post_y, -np.inf)
    peak_rel = int(np.argmax(finite_signal)) if finite_signal.size else 0
    peak_index = event_index + peak_rel
    peak_s = float(t_s[peak_index])
    peak_V = float(signal[peak_index])
    enough_tail = float(t_s[-1] - peak_s) >= float(criteria.peak_confirmation_s)
    after_peak = signal[peak_index:]
    no_later_higher_peak = bool(
        after_peak.size
        and np.nanmax(after_peak) <= peak_V * (1.0 + 1.0e-9) + 1.0e-15
    )
    dropped = bool(
        after_peak.size
        and peak_V > 0.0
        and np.nanmin(after_peak) <= (1.0 - criteria.peak_drop_fraction) * peak_V
    )
    peak_confirmed = bool(detected and enough_tail and no_later_higher_peak and dropped)
    fractions: dict[str, float | None] = {}
    if detected and peak_V > 0.0:
        for fraction in (0.10, 0.50):
            value = _first_crossing_interpolated(post_t, post_y, fraction * peak_V)
            fractions[f"{fraction:.2f}"] = (
                None if value is None else float((value - event_s) / 1.0e-12)
            )
    ready_s = (
        peak_s + float(criteria.post_peak_safety_s)
        if peak_confirmed
        else None
    )
    return {
        "detected": detected,
        "censored": not detected,
        "polarity": "positive" if polarity > 0.0 else "negative",
        "threshold_V": float(criteria.threshold_V),
        "threshold_uV": float(criteria.threshold_V * 1.0e6),
        "crossing_time_s": _finite_or_none(crossing_s),
        "crossing_time_ps": _ps_or_none(crossing_s),
        "t_lat_s": None if crossing_s is None else float(crossing_s - event_s),
        "t_lat_ps": None if crossing_s is None else float((crossing_s - event_s) / 1.0e-12),
        "lower_bound_ps": None if detected else float((t_s[-1] - event_s) / 1.0e-12),
        "peak_index": int(peak_index),
        "peak_time_s": peak_s,
        "peak_time_ps": peak_s / 1.0e-12,
        "peak_excursion_V": peak_V,
        "peak_excursion_uV": peak_V * 1.0e6,
        "peak_confirmed": peak_confirmed,
        "constant_fraction_latency_ps": fractions,
        "latency_stop_ready_time_s": _finite_or_none(ready_s),
        "latency_stop_ready_time_ps": _ps_or_none(ready_s),
    }


def _electrical_recovery_mask(
    history: Mapping[str, Any],
    *,
    n: int,
    baselines: Mapping[str, float],
    criteria: RecoveryCriteria,
) -> tuple[np.ndarray, dict[str, Any]]:
    current_keys = ("I_b_A", "I_s_A", "I_rf_A")
    voltage_keys = ("V_out_V", "v_c_V", "V_tdgl_center_V")
    mask = np.ones(n, dtype=bool)
    final_residuals: dict[str, float] = {}
    tolerances: dict[str, float] = {}
    for key in current_keys + voltage_keys:
        values = _series(history, key, n)
        baseline = float(baselines[key])
        if key in current_keys:
            tolerance = max(
                float(criteria.current_absolute_tolerance_A),
                float(criteria.current_relative_tolerance) * abs(baseline),
            )
        else:
            tolerance = max(
                float(criteria.voltage_absolute_tolerance_V),
                float(criteria.voltage_relative_tolerance) * abs(baseline),
            )
        residual = np.abs(values - baseline)
        mask &= np.isfinite(residual) & (residual <= tolerance)
        final_residuals[key] = float(residual[-1])
        tolerances[key] = float(tolerance)
    return mask, {
        "final_absolute_residuals": final_residuals,
        "absolute_tolerances": tolerances,
        "components": [*current_keys, *voltage_keys],
    }


def _efficiency_recovery_mask(
    history: Mapping[str, Any],
    *,
    n: int,
    baseline_current_A: float,
    fraction: float,
    start_index: int,
) -> np.ndarray:
    current = _series(history, "I_s_A", n)
    mask = np.zeros(n, dtype=bool)
    if not np.isfinite(baseline_current_A) or baseline_current_A == 0.0:
        return mask
    target = float(fraction) * abs(float(baseline_current_A))
    mask[start_index:] = np.abs(current[start_index:]) >= target
    return mask


def _state_recovery_result(
    history: Mapping[str, Any],
    *,
    snapshots: Mapping[str, Any] | None,
    history_t_s: np.ndarray,
    baselines: Mapping[str, float],
    electrical_mask: np.ndarray,
    event_s: float,
    event_index: int,
    start_index: int,
    enabled: bool,
    criteria: RecoveryCriteria,
) -> dict[str, Any]:
    n = history_t_s.size
    scalar_mask = electrical_mask.copy()
    scalar_details: dict[str, float] = {}
    for key in ("max_Te_K", "max_Tph_K"):
        values = _series(history, key, n)
        residual = np.abs(values - float(baselines[key]))
        scalar_mask &= residual <= float(criteria.temperature_absolute_tolerance_K)
        scalar_details[f"{key}_final_abs_K"] = float(residual[-1])
    for key in ("min_delta_over_delta0", "mean_delta_over_delta0"):
        values = _series(history, key, n)
        scale = max(abs(float(baselines[key])), 1.0e-12)
        residual = np.abs(values - float(baselines[key])) / scale
        scalar_mask &= residual <= float(criteria.condensate_relative_tolerance)
        scalar_details[f"{key}_final_rel"] = float(residual[-1])

    online_spatial_keys = (
        "state_delta_rms_relative",
        "state_Te_quantile_abs_K",
        "state_Tph_quantile_abs_K",
        "state_delta_max_relative",
        "state_Te_max_abs_K",
        "state_Tph_max_abs_K",
    )
    if all(key in history for key in online_spatial_keys):
        delta_rms = _series(history, online_spatial_keys[0], n)
        te_quantile = _series(history, online_spatial_keys[1], n)
        tph_quantile = _series(history, online_spatial_keys[2], n)
        delta_max = _series(history, online_spatial_keys[3], n)
        te_max = _series(history, online_spatial_keys[4], n)
        tph_max = _series(history, online_spatial_keys[5], n)
        guard = float(criteria.spatial_max_guard_factor)
        scalar_mask &= (
            (delta_rms <= float(criteria.condensate_relative_tolerance))
            & (te_quantile <= float(criteria.temperature_absolute_tolerance_K))
            & (tph_quantile <= float(criteria.temperature_absolute_tolerance_K))
            & (delta_max <= guard * float(criteria.condensate_relative_tolerance))
            & (te_max <= guard * float(criteria.temperature_absolute_tolerance_K))
            & (tph_max <= guard * float(criteria.temperature_absolute_tolerance_K))
        )
        return _recovery_result(
            history_t_s,
            scalar_mask,
            event_s=event_s,
            start_index=start_index,
            hold_s=criteria.hold_s,
            enabled=enabled,
            label="state",
            details={
                "scalar_details": scalar_details,
                "spatial_details": {
                    "available": True,
                    "source": "online history diagnostics",
                    "quantile": float(criteria.spatial_quantile),
                    "final_delta_rms_relative": float(delta_rms[-1]),
                    "final_Te_quantile_abs_K": float(te_quantile[-1]),
                    "final_Tph_quantile_abs_K": float(tph_quantile[-1]),
                    "max_guard_factor": guard,
                    "final_delta_max_relative": float(delta_max[-1]),
                    "final_Te_max_abs_K": float(te_max[-1]),
                    "final_Tph_max_abs_K": float(tph_max[-1]),
                },
            },
        )

    if snapshots is None:
        return {
            "mode": "state",
            "recovered": False,
            "censored": True,
            "reason": "state recovery requires transient snapshots",
            "lower_bound_ps": float((history_t_s[-1] - event_s) / 1.0e-12),
            "scalar_details": scalar_details,
        }

    snap_t = _snapshot_time_s(snapshots)
    if snap_t.size < 2:
        return {
            "mode": "state",
            "recovered": False,
            "censored": True,
            "reason": "snapshot history is insufficient for spatial recovery",
            "lower_bound_ps": float((history_t_s[-1] - event_s) / 1.0e-12),
            "scalar_details": scalar_details,
        }
    baseline_candidates = np.flatnonzero(snap_t < event_s)
    if not baseline_candidates.size:
        return {
            "mode": "state",
            "recovered": False,
            "censored": True,
            "reason": "no pre-photon spatial snapshot is available",
            "lower_bound_ps": float((history_t_s[-1] - event_s) / 1.0e-12),
            "scalar_details": scalar_details,
        }
    baseline_snapshot = int(baseline_candidates[-1])
    spatial_mask, spatial_details = _spatial_recovery_mask(
        snapshots,
        baseline_index=baseline_snapshot,
        criteria=criteria,
    )
    m = min(snap_t.size, spatial_mask.size)
    snap_t = snap_t[:m]
    spatial_mask = spatial_mask[:m]
    nearest_history = np.searchsorted(history_t_s, snap_t, side="left")
    nearest_history = np.clip(nearest_history, 0, n - 1)
    combined = spatial_mask & scalar_mask[nearest_history]
    start_s = float(history_t_s[start_index])
    snap_start = int(np.searchsorted(snap_t, start_s, side="left"))
    result = _recovery_result(
        snap_t,
        combined,
        event_s=event_s,
        start_index=snap_start,
        hold_s=criteria.hold_s,
        enabled=enabled,
        label="state",
        details={
            "scalar_details": scalar_details,
            "spatial_details": spatial_details,
            "baseline_snapshot_index": baseline_snapshot,
        },
    )
    return result


def _spatial_recovery_mask(
    snapshots: Mapping[str, Any],
    *,
    baseline_index: int,
    criteria: RecoveryCriteria,
) -> tuple[np.ndarray, dict[str, Any]]:
    real = np.asarray(snapshots.get("psi_real_snapshot_J", []), dtype=float)
    imag = np.asarray(snapshots.get("psi_imag_snapshot_J", []), dtype=float)
    te = np.asarray(snapshots.get("Te_snapshot_K", []), dtype=float)
    tph = np.asarray(snapshots.get("Tph_snapshot_K", []), dtype=float)
    arrays = [arr for arr in (real, imag, te, tph) if arr.ndim == 2]
    if len(arrays) < 4:
        size = min((arr.shape[0] for arr in arrays), default=0)
        return np.zeros(size, dtype=bool), {"available": False}
    n = min(arr.shape[0] for arr in arrays)
    psi_abs = np.hypot(real[:n], imag[:n])
    base_delta = psi_abs[baseline_index]
    delta_scale = max(float(np.sqrt(np.nanmean(base_delta**2))), 1.0e-300)
    delta_rms = np.sqrt(np.nanmean((psi_abs - base_delta[None, :]) ** 2, axis=1)) / delta_scale
    q = float(criteria.spatial_quantile)
    te_abs = np.nanquantile(np.abs(te[:n] - te[baseline_index][None, :]), q, axis=1)
    tph_abs = np.nanquantile(np.abs(tph[:n] - tph[baseline_index][None, :]), q, axis=1)
    guard = float(criteria.spatial_max_guard_factor)
    delta_max = np.nanmax(np.abs(psi_abs - base_delta[None, :]), axis=1) / delta_scale
    te_max = np.nanmax(np.abs(te[:n] - te[baseline_index][None, :]), axis=1)
    tph_max = np.nanmax(np.abs(tph[:n] - tph[baseline_index][None, :]), axis=1)
    mask = (
        np.isfinite(delta_rms)
        & (delta_rms <= float(criteria.condensate_relative_tolerance))
        & (te_abs <= float(criteria.temperature_absolute_tolerance_K))
        & (tph_abs <= float(criteria.temperature_absolute_tolerance_K))
        & (delta_max <= guard * float(criteria.condensate_relative_tolerance))
        & (te_max <= guard * float(criteria.temperature_absolute_tolerance_K))
        & (tph_max <= guard * float(criteria.temperature_absolute_tolerance_K))
    )
    return mask, {
        "available": True,
        "quantile": q,
        "final_delta_rms_relative": float(delta_rms[-1]),
        "final_Te_quantile_abs_K": float(te_abs[-1]),
        "final_Tph_quantile_abs_K": float(tph_abs[-1]),
        "max_guard_factor": guard,
        "final_delta_max_relative": float(delta_max[-1]),
        "final_Te_max_abs_K": float(te_max[-1]),
        "final_Tph_max_abs_K": float(tph_max[-1]),
    }


def _recovery_result(
    t_s: np.ndarray,
    mask: np.ndarray,
    *,
    event_s: float,
    start_index: int,
    hold_s: float,
    enabled: bool,
    label: str,
    details: Mapping[str, Any],
) -> dict[str, Any]:
    if not enabled:
        return {
            "mode": label,
            "recovered": False,
            "censored": True,
            "reason": "recovery is not classified before a confirmed detection",
            "lower_bound_ps": float((t_s[-1] - event_s) / 1.0e-12),
            **dict(details),
        }
    entry = _first_true_interval(t_s, mask, start_index=start_index, duration_s=hold_s)
    if entry is None:
        return {
            "mode": label,
            "recovered": False,
            "censored": True,
            "reason": "criterion was not maintained for the requested hold window",
            "lower_bound_ps": float((t_s[-1] - event_s) / 1.0e-12),
            **dict(details),
        }
    entry_s, confirmed_s = entry
    return {
        "mode": label,
        "recovered": True,
        "censored": False,
        "reason": "criterion maintained for the requested hold window",
        "entry_time_s": float(entry_s),
        "entry_time_ps": float(entry_s / 1.0e-12),
        "t_rec_s": float(entry_s - event_s),
        "t_rec_ps": float((entry_s - event_s) / 1.0e-12),
        "confirmed_time_s": float(confirmed_s),
        "confirmed_time_ps": float(confirmed_s / 1.0e-12),
        "hold_s": float(hold_s),
        "hold_ps": float(hold_s / 1.0e-12),
        **dict(details),
    }


def _confirmed_crossing_time(
    t_s: np.ndarray,
    y: np.ndarray,
    *,
    threshold: float,
    confirmation_s: float,
    hysteresis_fraction: float,
) -> float | None:
    if t_s.size < 2:
        return None
    candidates = np.flatnonzero(
        np.isfinite(y[:-1])
        & np.isfinite(y[1:])
        & (y[:-1] < threshold)
        & (y[1:] >= threshold)
    )
    for index in candidates:
        crossing = _interpolate(t_s[index], t_s[index + 1], y[index], y[index + 1], threshold)
        end = crossing + confirmation_s
        confirm = (t_s >= crossing) & (t_s <= end + 1.0e-30)
        if t_s[-1] < end - 1.0e-30:
            continue
        floor = threshold * (1.0 - hysteresis_fraction)
        if np.any(confirm) and np.nanmin(y[confirm]) >= floor:
            return float(crossing)
    return None


def _first_crossing_interpolated(
    t_s: np.ndarray,
    y: np.ndarray,
    threshold: float,
) -> float | None:
    candidates = np.flatnonzero(
        np.isfinite(y[:-1])
        & np.isfinite(y[1:])
        & (y[:-1] < threshold)
        & (y[1:] >= threshold)
    )
    if not candidates.size:
        return None
    i = int(candidates[0])
    return _interpolate(t_s[i], t_s[i + 1], y[i], y[i + 1], threshold)


def _first_true_interval(
    t_s: np.ndarray,
    mask: np.ndarray,
    *,
    start_index: int,
    duration_s: float,
) -> tuple[float, float] | None:
    start = max(0, int(start_index))
    active_start: int | None = None
    for index in range(start, min(t_s.size, mask.size)):
        if bool(mask[index]):
            if active_start is None:
                active_start = index
            if float(t_s[index] - t_s[active_start]) >= duration_s - 1.0e-30:
                return float(t_s[active_start]), float(t_s[index])
        else:
            active_start = None
    return None


def _event_index(history: Mapping[str, Any], n: int) -> int | None:
    applied = np.asarray(history.get("photon_applied", []), dtype=bool).reshape(-1)
    if applied.size:
        indices = np.flatnonzero(applied[:n])
        if indices.size:
            return int(indices[0])
    return None


def _resolve_polarity(requested: str, signal: np.ndarray) -> float:
    if requested == "positive":
        return 1.0
    if requested == "negative":
        return -1.0
    positive = float(np.nanmax(signal)) if signal.size else 0.0
    negative = abs(float(np.nanmin(signal))) if signal.size else 0.0
    return 1.0 if positive >= negative else -1.0


def _robust_baseline(
    history: Mapping[str, Any],
    key: str,
    mask: np.ndarray,
    n: int,
) -> float:
    values = _series(history, key, n)
    selected = values[mask & np.isfinite(values)]
    return float(np.nanmedian(selected)) if selected.size else float("nan")


def _series(history: Mapping[str, Any], key: str, n: int) -> np.ndarray:
    values = np.asarray(history.get(key, np.full(n, np.nan)), dtype=float).reshape(-1)
    if values.size != n:
        values = np.resize(values, n)
    return values


def _history_time_s(history: Mapping[str, Any]) -> np.ndarray:
    if "t_s" in history:
        values = np.asarray(history["t_s"], dtype=float).reshape(-1)
        if values.size:
            return values
    if "t_ps" in history:
        return np.asarray(history["t_ps"], dtype=float).reshape(-1) * 1.0e-12
    return np.empty(0, dtype=float)


def _snapshot_time_s(snapshots: Mapping[str, Any]) -> np.ndarray:
    if "snapshot_t_s" in snapshots:
        values = np.asarray(snapshots["snapshot_t_s"], dtype=float).reshape(-1)
        if values.size:
            return values
    if "snapshot_t_ps" in snapshots:
        return np.asarray(snapshots["snapshot_t_ps"], dtype=float).reshape(-1) * 1.0e-12
    return np.empty(0, dtype=float)


def _interpolate(t0: float, t1: float, y0: float, y1: float, target: float) -> float:
    if y1 == y0:
        return float(t1)
    fraction = (target - y0) / (y1 - y0)
    return float(t0 + fraction * (t1 - t0))


def _finite_or_none(value: Any) -> float | None:
    if value is None:
        return None
    number = float(value)
    return number if np.isfinite(number) else None


def _ps_or_none(value: Any) -> float | None:
    number = _finite_or_none(value)
    return None if number is None else float(number / 1.0e-12)


def _unavailable_result(
    detection: DetectionCriteria,
    recovery: RecoveryCriteria,
    reason: str,
) -> dict[str, Any]:
    return {
        "schema_version": 1,
        "available": False,
        "reason": reason,
        "detection_criteria": asdict(detection),
        "recovery_criteria": asdict(recovery),
        "latency": {"detected": False, "censored": True, "reason": reason},
        "recovery": {
            "selected": {
                "mode": recovery.mode,
                "recovered": False,
                "censored": True,
                "reason": reason,
            }
        },
        "termination": {"latency_ready": False, "recovery_ready": False},
    }


__all__ = [
    "DetectionCriteria",
    "RecoveryCriteria",
    "analyze_photon_timing",
]
