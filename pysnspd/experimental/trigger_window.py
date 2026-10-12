"""Opt-in observation horizon for an explicitly configured V_out trigger.

This module changes no physical equations and supplies no physical defaults.
Call observe_accepted only after a successful coupled integration step, including
an initial observation at the photon-transfer time. Times and voltages use SI.

A checkpoint here stores the observation controller only. A resumable physical
run must also preserve its full coupled state, circuit, clock, deposition state,
mesh, boundaries, parameters and numerical settings. This object is never a
substitute for that physical checkpoint.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass
import math
from typing import Any, Mapping

SCHEMA = "pysnspd.trigger_window.v1"


@dataclass(frozen=True)
class TriggerWindowConfig:
    """All measurement choices are required; polarity cannot depend on the tail."""

    event_time_s: float
    baseline_V: float
    baseline_reference: str
    threshold_V: float
    polarity: str
    confirmation_s: float
    hysteresis_fraction: float
    post_trigger_margin_s: float
    no_trigger_time_ceiling_s: float

    def validated(self) -> "TriggerWindowConfig":
        for name in ("event_time_s", "baseline_V", "threshold_V", "confirmation_s",
                     "hysteresis_fraction", "post_trigger_margin_s", "no_trigger_time_ceiling_s"):
            value = getattr(self, name)
            if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
                raise ValueError(name + " must be an explicit finite number")
        if self.event_time_s < 0:
            raise ValueError("event_time_s must be nonnegative")
        if self.threshold_V <= 0:
            raise ValueError("threshold_V must be positive")
        if self.confirmation_s < 0 or self.post_trigger_margin_s < 0:
            raise ValueError("confirmation and extra margin must be nonnegative")
        if not 0 <= self.hysteresis_fraction < 1:
            raise ValueError("hysteresis_fraction must lie in [0, 1)")
        if self.no_trigger_time_ceiling_s <= self.event_time_s:
            raise ValueError("no_trigger_time_ceiling_s must follow event_time_s")
        if self.polarity not in ("positive", "negative"):
            raise ValueError("an explicit positive or negative polarity is required")
        if not isinstance(self.baseline_reference, str) or not self.baseline_reference.strip():
            raise ValueError("a valid pre-event baseline reference is required")
        return self


def _time_tolerance(a: float, b: float) -> float:
    return 64 * math.ulp(max(abs(a), abs(b), 1e-300))


def _reached(time_s: float, target_s: float) -> bool:
    return time_s >= target_s - _time_tolerance(time_s, target_s)


class TriggerWindow:
    """Track the first accepted upward crossing and a finite post-trigger window.

    The crossing is interpolated between accepted observations. Confirmation
    uses accepted observations and checks the linear segments connecting them
    against the registered hysteresis floor. The caller must establish temporal
    resolution by refinement; sparse samples do not prove unsampled persistence.
    The no-trigger ceiling bounds the search for a confirmed crossing. Once a
    crossing is confirmed, its finite margin target may extend beyond that
    ceiling. There is no implicit absolute cap that cuts this margin short.
    The confirmation time is the accepted time when confirmation is established,
    separately from the crossing and its nominal confirmation deadline.
    """

    _state_keys = (
        "last_time_s", "last_excursion_V", "sample_count", "maximum_observation_step_s",
        "candidate_crossing_s", "crossing_time_s", "confirmation_time_s",
        "stop_target_s", "actual_stop_s", "status", "failure_reason",
        "rejected_candidate_count",
    )
    _terminal_statuses = frozenset((
        "TRIGGER_MARGIN_COMPLETE", "TIME_CENSORED_NO_TRIGGER",
        "INVALID_INITIAL_EVENT",
        "NUMERICAL_FAILURE",
    ))

    def __init__(self, config: TriggerWindowConfig):
        self.config = config.validated()
        self.last_time_s = None
        self.last_excursion_V = None
        self.sample_count = 0
        self.maximum_observation_step_s = 0.
        self.candidate_crossing_s = None
        self.crossing_time_s = None
        self.confirmation_time_s = None
        self.stop_target_s = None
        self.actual_stop_s = None
        self.status = "ACTIVE"
        self.failure_reason = None
        self.rejected_candidate_count = 0

    @property
    def terminal(self) -> bool:
        return self.status in self._terminal_statuses

    def next_required_time_s(self) -> float:
        """Return the next horizon cap; the integrator still chooses its step."""
        if self.terminal:
            if self.last_time_s is None:
                raise RuntimeError("no accepted time is available after failure")
            return float(self.last_time_s)
        if self.stop_target_s is not None:
            return float(self.stop_target_s)
        cap = self.config.no_trigger_time_ceiling_s
        if self.candidate_crossing_s is not None:
            deadline = self.candidate_crossing_s + self.config.confirmation_s
            if self.last_time_s is None or deadline > self.last_time_s:
                cap = min(cap, deadline)
        return float(cap)

    def observe_accepted(self, time_s: float, V_out_V: float) -> dict[str, Any]:
        """Consume an accepted output voltage; never call this for rejected steps."""
        if self.terminal:
            raise RuntimeError("the observation window is already terminal")
        for name, value in (("time_s", time_s), ("V_out_V", V_out_V)):
            if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
                raise ValueError(name + " must be finite; report numerical failure separately")
        t, voltage = float(time_s), float(V_out_V)
        if t < self.config.event_time_s - _time_tolerance(t, self.config.event_time_s):
            raise ValueError("observations must start at or after the transfer time")
        if self.crossing_time_s is None and t > (
            self.config.no_trigger_time_ceiling_s
            + _time_tolerance(t, self.config.no_trigger_time_ceiling_s)
        ):
            raise ValueError("accepted observation exceeds the finite no-trigger ceiling")
        if self.last_time_s is not None and t <= self.last_time_s:
            raise ValueError("accepted observation times must increase strictly")

        sign = 1.0 if self.config.polarity == "positive" else -1.0
        y = sign * (voltage - self.config.baseline_V)
        if not math.isfinite(y):
            raise ValueError("nonfinite excursion; report numerical failure separately")
        previous_t, previous_y = self.last_time_s, self.last_excursion_V
        self.last_time_s, self.last_excursion_V = t, y
        self.sample_count += 1
        if previous_t is not None:
            self.maximum_observation_step_s = max(self.maximum_observation_step_s, t - previous_t)
        threshold = self.config.threshold_V

        if previous_t is None:
            if not _reached(self.config.event_time_s, t):
                self.status = "INVALID_INITIAL_EVENT"
                self.failure_reason = "initial accepted observation is missing at transfer time"
            elif y >= threshold:
                self.status = "INVALID_INITIAL_EVENT"
                self.failure_reason = "initial output is already at or above threshold"
        elif self.crossing_time_s is None:
            floor = threshold * (1 - self.config.hysteresis_fraction)
            if self.candidate_crossing_s is not None:
                deadline = self.candidate_crossing_s + self.config.confirmation_s
                # If the segment extends past the nominal deadline, its value
                # at the deadline also tests the confirmation interval.
                end_y = y
                if previous_t < deadline < t:
                    fraction = (deadline - previous_t) / (t - previous_t)
                    end_y = previous_y * (1 - fraction) + y * fraction
                if min(previous_y, end_y) < floor:
                    self.candidate_crossing_s = None
                    self.rejected_candidate_count += 1
                elif _reached(t, deadline):
                    self._confirm(t)
            if self.candidate_crossing_s is None and self.crossing_time_s is None:
                if previous_y < threshold <= y:
                    scale = max(abs(previous_y), abs(y), threshold)
                    fraction = ((threshold / scale - previous_y / scale) /
                                (y / scale - previous_y / scale))
                    cross = previous_t + fraction * (t - previous_t)
                    self.candidate_crossing_s = cross
                    if _reached(t, cross + self.config.confirmation_s):
                        self._confirm(t)

        if self.status == "INVALID_INITIAL_EVENT":
            self.actual_stop_s = t
        elif self.stop_target_s is not None and _reached(t, self.stop_target_s):
            self.status = "TRIGGER_MARGIN_COMPLETE"
            self.actual_stop_s = t
        elif self.crossing_time_s is None and _reached(t, self.config.no_trigger_time_ceiling_s):
            self.status = "TIME_CENSORED_NO_TRIGGER"
            self.actual_stop_s = t
        return self.summary()

    def _confirm(self, accepted_time_s: float) -> None:
        self.crossing_time_s = self.candidate_crossing_s
        self.confirmation_time_s = accepted_time_s
        self.stop_target_s = max(
            self.crossing_time_s + self.config.post_trigger_margin_s,
            self.confirmation_time_s,
        )
        self.status = "TRIGGER_CONFIRMED"

    def report_numerical_failure(self, reason: str) -> dict[str, Any]:
        """Record failure separately from time-window censoring."""
        if self.terminal:
            raise RuntimeError("the observation window is already terminal")
        if not isinstance(reason, str) or not reason.strip():
            raise ValueError("an explicit failure reason is required")
        self.status = "NUMERICAL_FAILURE"
        self.failure_reason = reason
        self.actual_stop_s = self.last_time_s
        return self.summary()

    def summary(self) -> dict[str, Any]:
        detected = self.crossing_time_s is not None
        censored = self.status == "TIME_CENSORED_NO_TRIGGER"
        elapsed = (None if self.last_time_s is None else
                   self.last_time_s - self.config.event_time_s)
        return dict(
            schema=SCHEMA, config=asdict(self.config), observable="V_out_V",
            baseline_relative=True, status=self.status, terminal=self.terminal,
            detected=detected, censored=censored,
            crossing_time_s=self.crossing_time_s,
            confirmation_deadline_s=(None if self.crossing_time_s is None else
                                     self.crossing_time_s + self.config.confirmation_s),
            confirmation_time_s=self.confirmation_time_s,
            latency_s=(None if not detected else
                       self.crossing_time_s - self.config.event_time_s),
            censored_lower_bound_s=elapsed if censored else None,
            unconfirmed_crossing_s=(self.candidate_crossing_s if not detected else None),
            stop_target_s=self.stop_target_s, actual_stop_s=self.actual_stop_s,
            margin_complete=self.status == "TRIGGER_MARGIN_COMPLETE",
            stop_overshoot_s=(None if self.actual_stop_s is None or self.stop_target_s is None
                              else max(0., self.actual_stop_s - self.stop_target_s)),
            sample_count=self.sample_count,
            maximum_observation_step_s=self.maximum_observation_step_s,
            confirmation_interpolation="linear between accepted observations",
            rejected_candidate_count=self.rejected_candidate_count,
            failure_reason=self.failure_reason,
            requires_peak_confirmation=False, requires_full_recovery=False,
            checkpoint_scope="observation_controller_only",
        )

    def checkpoint(self) -> dict[str, Any]:
        """JSON-ready controller state; pair with a complete physical checkpoint."""
        return dict(schema=SCHEMA, config=asdict(self.config),
                    state={key: getattr(self, key) for key in self._state_keys})

    @classmethod
    def from_checkpoint(cls, checkpoint: Mapping[str, Any]) -> "TriggerWindow":
        """Restore a controller without recomputing its baseline or trigger."""
        if checkpoint.get("schema") != SCHEMA:
            raise ValueError("unsupported trigger-window checkpoint")
        result = cls(TriggerWindowConfig(**checkpoint["config"]))
        state = dict(checkpoint["state"])
        if set(state) != set(cls._state_keys):
            raise ValueError("incomplete or unexpected controller checkpoint state")
        for key in ("sample_count", "rejected_candidate_count"):
            if type(state[key]) is not int or state[key] < 0:
                raise ValueError("invalid checkpoint count: " + key)
        numeric = set(cls._state_keys) - {"sample_count", "rejected_candidate_count", "status", "failure_reason"}
        for key in numeric:
            value = state[key]
            if value is not None and (isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value)):
                raise ValueError("nonfinite checkpoint field: " + key)
        if state["status"] not in cls._terminal_statuses | {"ACTIVE", "TRIGGER_CONFIRMED"}:
            raise ValueError("invalid controller status")
        if state["failure_reason"] is not None and not isinstance(state["failure_reason"], str):
            raise ValueError("invalid checkpoint failure reason")
        if (state["last_time_s"] is None) != (state["sample_count"] == 0):
            raise ValueError("checkpoint accepted time and count disagree")
        if (state["last_time_s"] is None) != (state["last_excursion_V"] is None):
            raise ValueError("checkpoint accepted observation is incomplete")
        if state["last_time_s"] is not None and state["last_time_s"] < (
            result.config.event_time_s - _time_tolerance(state["last_time_s"], result.config.event_time_s)
        ):
            raise ValueError("checkpoint time precedes the observation window")
        if state["last_time_s"] is not None and state["crossing_time_s"] is None and state["last_time_s"] > (
            result.config.no_trigger_time_ceiling_s
            + _time_tolerance(state["last_time_s"], result.config.no_trigger_time_ceiling_s)
        ):
            raise ValueError("checkpoint time exceeds the no-trigger observation window")
        if state["maximum_observation_step_s"] is None or state["maximum_observation_step_s"] < 0:
            raise ValueError("checkpoint observation resolution is invalid")
        if (state["maximum_observation_step_s"] == 0) != (state["sample_count"] < 2):
            raise ValueError("checkpoint observation resolution and count disagree")
        confirmed = state["crossing_time_s"] is not None
        if confirmed and state["last_time_s"] is None:
            raise ValueError("checkpoint confirmed trigger has no accepted observations")
        if state["candidate_crossing_s"] is not None and state["last_time_s"] is None:
            raise ValueError("checkpoint candidate has no accepted observations")
        if state["sample_count"] < 2 and (
            state["candidate_crossing_s"] is not None or state["crossing_time_s"] is not None
        ):
            raise ValueError("checkpoint crossing requires at least two accepted observations")
        if state["sample_count"] == 1 and state["status"] != "INVALID_INITIAL_EVENT":
            if abs(state["last_time_s"] - result.config.event_time_s) > _time_tolerance(
                state["last_time_s"], result.config.event_time_s
            ) or state["last_excursion_V"] >= result.config.threshold_V:
                raise ValueError("checkpoint unique observation is inconsistent with a valid initial event")
        if confirmed != (state["confirmation_time_s"] is not None) or confirmed != (state["stop_target_s"] is not None):
            raise ValueError("checkpoint confirmed trigger is incomplete")
        if confirmed and not (
            result.config.event_time_s <= state["crossing_time_s"] <= state["confirmation_time_s"] <= state["last_time_s"]
            and state["confirmation_time_s"] <= result.config.no_trigger_time_ceiling_s + _time_tolerance(
                state["confirmation_time_s"], result.config.no_trigger_time_ceiling_s)
            and _reached(state["confirmation_time_s"], state["crossing_time_s"] + result.config.confirmation_s)
        ):
            raise ValueError("checkpoint trigger times are inconsistent")
        if confirmed and state["stop_target_s"] != max(
            state["crossing_time_s"] + result.config.post_trigger_margin_s,
            state["confirmation_time_s"],
        ):
            raise ValueError("checkpoint stop target disagrees with the declared margin")
        if state["candidate_crossing_s"] is not None and not (
            result.config.event_time_s <= state["candidate_crossing_s"] <= state["last_time_s"]
        ):
            raise ValueError("checkpoint candidate is outside accepted observations")
        if state["status"] in ("TRIGGER_CONFIRMED", "TRIGGER_MARGIN_COMPLETE") and not confirmed:
            raise ValueError("checkpoint status requires a confirmed trigger")
        if state["status"] == "TIME_CENSORED_NO_TRIGGER" and confirmed:
            raise ValueError("checkpoint censoring conflicts with a confirmed trigger")
        if confirmed and state["candidate_crossing_s"] != state["crossing_time_s"]:
            raise ValueError("checkpoint candidate and confirmed crossing disagree")
        if state["rejected_candidate_count"] > state["sample_count"]:
            raise ValueError("checkpoint rejected count exceeds accepted observations")
        terminal = state["status"] in cls._terminal_statuses
        if (state["actual_stop_s"] is not None) != (terminal and state["last_time_s"] is not None):
            raise ValueError("checkpoint actual stop and terminal status disagree")
        if state["actual_stop_s"] is not None and state["actual_stop_s"] != state["last_time_s"]:
            raise ValueError("checkpoint actual stop differs from its accepted state")
        failed = state["status"] in ("INVALID_INITIAL_EVENT", "NUMERICAL_FAILURE")
        if failed != bool(state["failure_reason"]):
            raise ValueError("checkpoint failure reason and status disagree")
        if state["status"] == "ACTIVE" and state["last_time_s"] is not None and _reached(
            state["last_time_s"], result.config.no_trigger_time_ceiling_s
        ):
            raise ValueError("active checkpoint has already reached its ceiling")
        if state["status"] == "ACTIVE" and confirmed:
            raise ValueError("active checkpoint cannot contain a confirmed trigger")
        if state["status"] == "INVALID_INITIAL_EVENT" and (state["sample_count"] != 1 or confirmed):
            raise ValueError("invalid initial event checkpoint is inconsistent")
        if state["status"] == "TIME_CENSORED_NO_TRIGGER":
            if state["last_time_s"] is None or not _reached(state["last_time_s"], result.config.no_trigger_time_ceiling_s):
                raise ValueError("checkpoint temporal censoring requires the finite ceiling")
        if state["status"] == "TRIGGER_MARGIN_COMPLETE" and not _reached(
            state["last_time_s"], state["stop_target_s"]
        ):
            raise ValueError("checkpoint margin is incomplete")
        if state["status"] == "TRIGGER_CONFIRMED" and _reached(
            state["last_time_s"], state["stop_target_s"]
        ):
            raise ValueError("checkpoint should already have a terminal status")
        if state["status"] == "ACTIVE" and state["candidate_crossing_s"] is not None and _reached(
            state["last_time_s"], state["candidate_crossing_s"] + result.config.confirmation_s
        ):
            raise ValueError("checkpoint pending crossing has already reached confirmation")
        for key, value in state.items():
            setattr(result, key, value)
        return result


__all__ = ["TriggerWindowConfig", "TriggerWindow", "SCHEMA"]
