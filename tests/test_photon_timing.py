"""Operational latency and recovery metrics from saved photon histories."""

from __future__ import annotations

import numpy as np

from pysnspd.analysis.timing import (
    DetectionCriteria,
    RecoveryCriteria,
    analyze_photon_timing,
)


def _synthetic_detected_history() -> dict[str, np.ndarray]:
    t_ps = np.arange(0.0, 100.5, 0.5)
    excursion_uV = np.zeros_like(t_ps)
    rising = (t_ps >= 20.0) & (t_ps <= 30.0)
    falling = (t_ps > 30.0) & (t_ps <= 60.0)
    excursion_uV[rising] = 1000.0 * (t_ps[rising] - 20.0) / 10.0
    excursion_uV[falling] = 1000.0 * (60.0 - t_ps[falling]) / 30.0
    voltage = excursion_uV * 1.0e-6
    I_b = np.full_like(t_ps, 30.0e-6)
    I_rf = voltage / 50.0
    I_s = I_b - I_rf
    photon = t_ps >= 20.0
    recovery_envelope = np.zeros_like(t_ps)
    recovery_envelope[rising] = (t_ps[rising] - 20.0) / 10.0
    recovery_envelope[falling] = (60.0 - t_ps[falling]) / 30.0
    return {
        "t_s": t_ps * 1.0e-12,
        "photon_applied": photon,
        "I_b_A": I_b,
        "I_s_A": I_s,
        "I_rf_A": I_rf,
        "V_out_V": voltage,
        "v_c_V": voltage,
        "V_tdgl_center_V": np.zeros_like(t_ps),
        "min_delta_over_delta0": 1.0 - 0.2 * recovery_envelope,
        "mean_delta_over_delta0": 1.0 - 0.05 * recovery_envelope,
        "max_Te_K": 1.0 + recovery_envelope,
        "max_Tph_K": 1.0 + 0.5 * recovery_envelope,
        "state_delta_rms_relative": 0.1 * recovery_envelope,
        "state_Te_quantile_abs_K": recovery_envelope,
        "state_Tph_quantile_abs_K": 0.5 * recovery_envelope,
        "state_delta_max_relative": 0.2 * recovery_envelope,
        "state_Te_max_abs_K": 2.0 * recovery_envelope,
        "state_Tph_max_abs_K": recovery_envelope,
    }


def test_latency_uses_confirmed_interpolated_leading_edge():
    result = analyze_photon_timing(_synthetic_detected_history())

    assert result["latency"]["detected"] is True
    assert np.isclose(result["latency"]["t_lat_ps"], 1.0)
    assert result["latency"]["peak_confirmed"] is True
    assert result["termination"]["latency_ready"] is True


def test_electrical_recovery_is_default_and_requires_hold_window():
    result = analyze_photon_timing(_synthetic_detected_history())
    selected = result["recovery"]["selected"]

    assert selected["mode"] == "electrical"
    assert selected["recovered"] is True
    assert selected["confirmed_time_ps"] - selected["entry_time_ps"] >= 10.0
    assert result["termination"]["recovery_ready"] is True


def test_state_and_efficiency_recovery_are_selectable():
    history = _synthetic_detected_history()
    state = analyze_photon_timing(
        history,
        recovery=RecoveryCriteria(mode="state"),
    )
    efficiency = analyze_photon_timing(
        history,
        recovery=RecoveryCriteria(mode="efficiency90"),
    )

    assert state["recovery"]["selected"]["mode"] == "state"
    assert state["recovery"]["selected"]["recovered"] is True
    assert efficiency["recovery"]["selected"]["mode"] == "efficiency90"
    assert efficiency["recovery"]["selected"]["recovered"] is True


def test_no_threshold_crossing_is_reported_as_right_censored():
    history = _synthetic_detected_history()
    history["V_out_V"] = 0.01 * history["V_out_V"]
    result = analyze_photon_timing(
        history,
        detection=DetectionCriteria(threshold_V=100.0e-6),
    )

    assert result["latency"]["detected"] is False
    assert result["latency"]["censored"] is True
    assert result["recovery"]["selected"]["recovered"] is False
    assert result["termination"]["recovery_ready"] is False
