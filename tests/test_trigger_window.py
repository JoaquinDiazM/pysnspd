"""Synthetic horizon controls; no physical threshold or material admission."""
import json

import pytest

from pysnspd.experimental.trigger_window import TriggerWindow, TriggerWindowConfig


def config(**changes):
    # Synthetic SI choices serve algebraic tests only.
    values = dict(event_time_s=0., baseline_V=3., baseline_reference="synthetic accepted DC",
                  threshold_V=1., polarity="positive", confirmation_s=.25,
                  hysteresis_fraction=.1, post_trigger_margin_s=1., no_trigger_time_ceiling_s=4.)
    values.update(changes)
    return TriggerWindowConfig(**values)


def test_rising_signal_stops_without_peak_or_recovery():
    window = TriggerWindow(config())
    window.observe_accepted(0., 3.)
    confirmed = window.observe_accepted(1., 5.)
    assert confirmed["detected"]
    assert confirmed["crossing_time_s"] == pytest.approx(.5)
    assert window.next_required_time_s() == pytest.approx(1.5)
    final = window.observe_accepted(1.5, 6.)
    assert final["margin_complete"]


def test_first_crossing_margin_and_confirmation_are_separate():
    window = TriggerWindow(config())
    window.observe_accepted(0., 3.)
    result = window.observe_accepted(.5, 4.)
    assert result["unconfirmed_crossing_s"] == pytest.approx(.5)
    assert window.next_required_time_s() == pytest.approx(.75)
    confirmed = window.observe_accepted(.75, 4.5)
    assert confirmed["latency_s"] == pytest.approx(.5)
    assert confirmed["confirmation_time_s"] == pytest.approx(.75)
    assert window.next_required_time_s() == pytest.approx(1.5)
    final = window.observe_accepted(1.5, 6.)
    assert final["status"] == "TRIGGER_MARGIN_COMPLETE"
    assert final["stop_overshoot_s"] == 0
    assert not final["requires_peak_confirmation"]
    assert not final["requires_full_recovery"]


def test_coarse_confirmation_and_zero_margin_stop_at_accepted_time():
    window = TriggerWindow(config(post_trigger_margin_s=0.))
    window.observe_accepted(0., 3.)
    final = window.observe_accepted(1., 5.)
    assert final["crossing_time_s"] == pytest.approx(.5)
    assert final["confirmation_deadline_s"] == pytest.approx(.75)
    assert final["confirmation_time_s"] == 1.
    assert final["stop_target_s"] == 1.
    assert final["actual_stop_s"] == 1.


def test_rejected_crossing_then_first_confirmed_recrossing():
    window = TriggerWindow(config(confirmation_s=.5))
    for time, voltage in [(0., 3.), (.5, 4.), (.75, 3.5), (1., 4.), (1.5, 4.2)]:
        result = window.observe_accepted(time, voltage)
    assert result["rejected_candidate_count"] == 1
    assert result["crossing_time_s"] == 1.
    assert result["confirmation_time_s"] == 1.5


def test_confirmation_uses_segment_value_at_nominal_deadline():
    window = TriggerWindow(config(confirmation_s=.1, post_trigger_margin_s=2.))
    window.observe_accepted(0., 3.)
    window.observe_accepted(.5, 4.)
    # Falling below floor after the confirmation deadline does not undo a
    # crossing that already met its declared persistence interval.
    result = window.observe_accepted(.75, 3.8)
    assert result["detected"]
    assert result["crossing_time_s"] == .5


def test_no_trigger_ceiling_does_not_cut_a_confirmed_margin_short():
    window = TriggerWindow(config(no_trigger_time_ceiling_s=1.))
    window.observe_accepted(0., 3.)
    window.observe_accepted(.5, 4.)
    window.observe_accepted(.75, 4.5)
    at_search_ceiling = window.observe_accepted(1., 5.)
    assert at_search_ceiling["status"] == "TRIGGER_CONFIRMED"
    assert window.next_required_time_s() == 1.5
    final = window.observe_accepted(1.5, 6.)
    assert final["status"] == "TRIGGER_MARGIN_COMPLETE"
    assert final["latency_s"] == .5
    assert not final["censored"]
    assert final["margin_complete"]


def test_no_trigger_and_unconfirmed_trigger_are_censored_without_latency():
    for voltage in (3.5, 4.):
        window = TriggerWindow(config(no_trigger_time_ceiling_s=.5))
        window.observe_accepted(0., 3.)
        final = window.observe_accepted(.5, voltage)
        assert final["status"] == "TIME_CENSORED_NO_TRIGGER"
        assert final["censored"]
        assert final["latency_s"] is None
        assert final["censored_lower_bound_s"] == .5


@pytest.mark.parametrize("initial_time,initial_voltage", [(0., 4.), (.1, 3.)])
def test_no_invented_crossing_for_invalid_initial_observation(initial_time, initial_voltage):
    final = TriggerWindow(config()).observe_accepted(initial_time, initial_voltage)
    assert final["status"] == "INVALID_INITIAL_EVENT"
    assert final["latency_s"] is None
    assert not final["censored"]


def test_negative_polarity_and_nonzero_baseline():
    window = TriggerWindow(config(polarity="negative", confirmation_s=0.))
    window.observe_accepted(0., 3.)
    final = window.observe_accepted(2., 1.)
    assert final["latency_s"] == 1.
    assert final["margin_complete"]


def test_numerical_failure_is_separate_and_invalid_sample_is_not_consumed():
    window = TriggerWindow(config())
    window.observe_accepted(0., 3.)
    before = window.checkpoint()
    with pytest.raises(ValueError, match="finite"):
        window.observe_accepted(.5, float("nan"))
    assert before == window.checkpoint()
    final = window.report_numerical_failure("coupled state validation failed")
    assert final["status"] == "NUMERICAL_FAILURE"
    assert not final["censored"]
    assert final["latency_s"] is None


@pytest.mark.parametrize("times", [(0., 0.), (0., -1.), (0., 5.)])
def test_duplicate_rejected_or_over_ceiling_times_do_not_advance(times):
    window = TriggerWindow(config())
    window.observe_accepted(times[0], 3.)
    before = window.checkpoint()
    with pytest.raises(ValueError):
        window.observe_accepted(times[1], 3.5)
    assert before == window.checkpoint()


@pytest.mark.parametrize("prefix", [
    [(0., 3.)],
    [(0., 3.), (.5, 4.)],
    [(0., 3.), (.5, 4.), (.75, 4.5)],
])
def test_controller_checkpoint_roundtrip_preserves_baseline_candidate_and_target(prefix):
    original = TriggerWindow(config())
    for sample in prefix:
        original.observe_accepted(*sample)
    restored = TriggerWindow.from_checkpoint(json.loads(json.dumps(original.checkpoint())))
    assert restored.summary() == original.summary()
    for sample in [(t, v) for t, v in [(.5, 4.), (.75, 4.5), (1.5, 6.)]
                   if t > prefix[-1][0]]:
        assert restored.observe_accepted(*sample) == original.observe_accepted(*sample)
    assert restored.summary()["margin_complete"]
    assert restored.summary()["checkpoint_scope"] == "observation_controller_only"


@pytest.mark.parametrize("field,value", [
    ("threshold_V", 0.), ("threshold_V", float("nan")), ("polarity", "auto"),
    ("baseline_reference", ""), ("confirmation_s", -.1),
    ("post_trigger_margin_s", -.1), ("no_trigger_time_ceiling_s", 0.),
    ("hysteresis_fraction", 1.), ("event_time_s", True),
])
def test_all_operational_choices_must_be_explicit_and_valid(field, value):
    with pytest.raises(ValueError):
        TriggerWindow(config(**{field: value}))


def test_checkpoint_rejects_incomplete_and_inconsistent_state():
    window = TriggerWindow(config())
    window.observe_accepted(0., 3.)
    checkpoint = window.checkpoint()
    del checkpoint["state"]["last_excursion_V"]
    with pytest.raises(ValueError, match="incomplete"):
        TriggerWindow.from_checkpoint(checkpoint)
    checkpoint = window.checkpoint()
    checkpoint["state"]["status"] = "TRIGGER_CONFIRMED"
    with pytest.raises(ValueError, match="requires"):
        TriggerWindow.from_checkpoint(checkpoint)


def test_picosecond_scale_and_nonzero_transfer_clock():
    window = TriggerWindow(config(event_time_s=10e-12, no_trigger_time_ceiling_s=20e-12,
                                  confirmation_s=.25e-12, post_trigger_margin_s=1e-12))
    window.observe_accepted(10e-12, 3.)
    window.observe_accepted(10.5e-12, 4.)
    window.observe_accepted(10.75e-12, 4.5)
    final = window.observe_accepted(11.5e-12, 6.)
    assert final["latency_s"] == pytest.approx(.5e-12, abs=1e-27)
    assert final["margin_complete"]


def test_unclamped_last_step_reports_overshoot():
    window = TriggerWindow(config())
    for sample in [(0., 3.), (.5, 4.), (.75, 4.5), (2., 6.)]:
        final = window.observe_accepted(*sample)
    assert final["stop_target_s"] == 1.5
    assert final["actual_stop_s"] == 2.
    assert final["stop_overshoot_s"] == .5


@pytest.mark.parametrize("trajectory, failure", [
    ([], True),
    ([(0., 3.)], True),
    ([(0., 3.), (.5, 4.), (.75, 4.5), (1.5, 6.)], False),
    ([(0., 4.)], False),
    ([(0., 3.), (4., 3.)], False),
])
def test_terminal_and_failure_checkpoint_roundtrip(trajectory, failure):
    original = TriggerWindow(config())
    for sample in trajectory:
        original.observe_accepted(*sample)
    if failure:
        original.report_numerical_failure("synthetic state failure")
    restored = TriggerWindow.from_checkpoint(json.loads(json.dumps(original.checkpoint())))
    assert restored.summary() == original.summary()


def test_finite_large_voltage_bracket_does_not_overflow_crossing():
    window = TriggerWindow(config(baseline_V=0., confirmation_s=0., post_trigger_margin_s=0.))
    window.observe_accepted(0., -1e308)
    final = window.observe_accepted(2., 1e308)
    assert final["crossing_time_s"] == pytest.approx(1.)
    assert final["detected"]


def test_checkpoint_rejects_inconsistent_actual_stop():
    window = TriggerWindow(config())
    window.observe_accepted(0., 3.)
    checkpoint = window.checkpoint()
    checkpoint["state"]["actual_stop_s"] = 0.
    with pytest.raises(ValueError, match="actual stop"):
        TriggerWindow.from_checkpoint(checkpoint)


def test_checkpoint_rejects_fabricated_partial_trigger_without_observations():
    checkpoint = TriggerWindow(config()).checkpoint()
    checkpoint["state"]["candidate_crossing_s"] = .5
    with pytest.raises(ValueError, match="no accepted observations"):
        TriggerWindow.from_checkpoint(checkpoint)


def test_checkpoint_rejects_altered_crossing_and_stale_active_ceiling():
    window = TriggerWindow(config())
    for sample in [(0., 3.), (.5, 4.), (.75, 4.5)]:
        window.observe_accepted(*sample)
    checkpoint = window.checkpoint()
    checkpoint["state"]["crossing_time_s"] = 0.
    with pytest.raises(ValueError):
        TriggerWindow.from_checkpoint(checkpoint)
    window = TriggerWindow(config())
    window.observe_accepted(0., 3.)
    window.observe_accepted(4., 3.)
    checkpoint = window.checkpoint()
    checkpoint["state"]["status"] = "ACTIVE"
    checkpoint["state"]["actual_stop_s"] = None
    with pytest.raises(ValueError, match="ceiling"):
        TriggerWindow.from_checkpoint(checkpoint)


def test_observation_resolution_is_measured_without_invented_tolerance():
    window = TriggerWindow(config())
    for sample in [(0., 3.), (.5, 4.), (.75, 4.5)]:
        result = window.observe_accepted(*sample)
    assert result["maximum_observation_step_s"] == .5
    assert result["confirmation_interpolation"] == "linear between accepted observations"


def test_late_199ps_trigger_completes_10ps_margin_at_209ps():
    # Synthetic electrical threshold; only the user's 200ps search ceiling is
    # represented. This test does not select a detector threshold or margin.
    window = TriggerWindow(config(baseline_V=0., threshold_V=100e-6,
                                  confirmation_s=.5e-12, post_trigger_margin_s=10e-12,
                                  no_trigger_time_ceiling_s=200e-12))
    window.observe_accepted(0., 0.)
    window.observe_accepted(198.5e-12, 0.)
    confirmed = window.observe_accepted(199.5e-12, 200e-6)
    assert confirmed["crossing_time_s"] == pytest.approx(199e-12, abs=1e-25)
    assert confirmed["confirmation_time_s"] == 199.5e-12
    assert window.next_required_time_s() == pytest.approx(209e-12, abs=1e-25)
    beyond_search_ceiling = window.observe_accepted(201e-12, 250e-6)
    assert beyond_search_ceiling["status"] == "TRIGGER_CONFIRMED"
    restored = TriggerWindow.from_checkpoint(json.loads(json.dumps(window.checkpoint())))
    final = restored.observe_accepted(209e-12, 300e-6)
    assert final["margin_complete"]
    assert not final["censored"]
    assert final["actual_stop_s"] == 209e-12


def test_unconfirmed_late_crossing_remains_censored_at_200ps():
    window = TriggerWindow(config(baseline_V=0., threshold_V=100e-6,
                                  confirmation_s=.5e-12, post_trigger_margin_s=10e-12,
                                  no_trigger_time_ceiling_s=200e-12))
    window.observe_accepted(0., 0.)
    window.observe_accepted(199.6e-12, 0.)
    final = window.observe_accepted(200e-12, 200e-6)
    assert final["status"] == "TIME_CENSORED_NO_TRIGGER"
    assert final["unconfirmed_crossing_s"] == pytest.approx(199.8e-12, abs=1e-25)
    assert final["latency_s"] is None


def test_initial_accepted_clock_roundtrip_preserves_nextafter_time():
    import math
    event_time = 10e-12
    accepted_time = math.nextafter(event_time, 0.)
    original = TriggerWindow(config(event_time_s=event_time,
                                   no_trigger_time_ceiling_s=200e-12))
    original.observe_accepted(accepted_time, 3.)
    assert original.summary()["status"] == "ACTIVE"
    restored = TriggerWindow.from_checkpoint(json.loads(json.dumps(original.checkpoint())))
    assert restored.checkpoint() == original.checkpoint()
    assert restored.last_time_s == accepted_time
    assert restored.last_time_s != event_time


@pytest.mark.parametrize("status", ["ACTIVE", "TRIGGER_CONFIRMED"])
def test_checkpoint_crossing_cannot_be_fabricated_from_one_observation(status):
    window = TriggerWindow(config())
    window.observe_accepted(0., 3.)
    checkpoint = window.checkpoint()
    checkpoint["state"].update(last_time_s=.75, last_excursion_V=1.5,
                               candidate_crossing_s=.5, status=status)
    if status == "TRIGGER_CONFIRMED":
        checkpoint["state"].update(crossing_time_s=.5, confirmation_time_s=.75,
                                   stop_target_s=1.5)
    with pytest.raises(ValueError, match="at least two"):
        TriggerWindow.from_checkpoint(checkpoint)


@pytest.mark.parametrize("field,value", [("last_time_s", .75), ("last_excursion_V", 1.)])
def test_active_unique_observation_requires_valid_transfer_time_and_voltage(field, value):
    window = TriggerWindow(config())
    window.observe_accepted(0., 3.)
    checkpoint = window.checkpoint()
    checkpoint["state"][field] = value
    with pytest.raises(ValueError, match="unique observation"):
        TriggerWindow.from_checkpoint(checkpoint)
