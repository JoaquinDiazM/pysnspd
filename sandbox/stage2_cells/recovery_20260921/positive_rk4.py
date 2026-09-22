"""RK4 with local interval bisection after explicitly classified support failure.

This is a corrective integration experiment, not an admission certificate.
The supplied validator alone classifies support violations. All other errors
propagate. No state, energy, or population is projected, clipped, or repaired.
"""
from __future__ import annotations

import hashlib
from pathlib import Path
import numpy as np


class SupportViolation(ValueError):
    """Raised by validate(state), with explicit sector/index/value information."""

    def __init__(self, message, details=None):
        super().__init__(message)
        self.details = dict(details or {})


class IntegrationLimitError(RuntimeError):
    """A hard RHS/depth/representability limit stopped an incomplete trajectory."""

    def __init__(self, message, stats):
        super().__init__(message)
        self.stats = stats


class _RejectedStage(Exception):
    def __init__(self, stage, time, violation):
        self.stage = stage
        self.time = float(time)
        self.violation = violation


def _jsonable(value):
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, dict):
        return {str(key): _jsonable(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_jsonable(item) for item in value]
    return value


def integrate(rhs, initial, duration, steps, validate, *, max_depth=20,
              max_rhs_calls=100000, callback=None, stats=None):
    """Return macro output times, complete states, and a separate event ledger.

    ``validate(state)`` must return None or raise SupportViolation. Invalid
    initial states, nonfinite values and arbitrary RHS/validator exceptions
    abort; they are not reasons to silently reduce a step. Only a trial stage
    or endpoint classified by this validator triggers local interval bisection.

    A successful left half is retained when advancing the right half. The
    accepted macro states are written only after both halves have completed.
    Pass an empty mutable ``stats`` dict to retain all diagnostics on failure.
    ``callback`` receives metadata events only, never writable state vectors.
    With no rejection the arithmetic and macro times match rk4_trajectory.
    """
    duration = float(duration)
    if not np.isfinite(duration) or duration <= 0:
        raise ValueError('duration must be finite and positive')
    if not isinstance(steps, (int, np.integer)) or isinstance(steps, bool) or steps < 1:
        raise ValueError('steps must be a positive integer')
    if type(max_depth) is not int or max_depth < 0:
        raise ValueError('max_depth must be a nonnegative integer')
    if type(max_rhs_calls) is not int or max_rhs_calls < 1:
        raise ValueError('max_rhs_calls must be a positive integer')
    if not callable(rhs) or not callable(validate) or (callback is not None and not callable(callback)):
        raise TypeError('rhs, validate and optional callback must be callable')
    initial = np.asarray(initial, float)
    if np.any(~np.isfinite(initial)):
        raise ValueError('initial state must be finite')
    if stats is None:
        stats = {}
    elif not isinstance(stats, dict) or stats:
        raise ValueError('stats must be an empty mutable dictionary')
    times = np.linspace(0, duration, steps+1)
    states = np.empty((steps+1, *initial.shape))
    states[0] = initial
    dt = duration/steps
    stats.update(schema='pysnspd.stage2.positive-rk4-experiment.v1',
        status='RUNNING', method='classical RK4 with local support-triggered bisection',
        source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        duration=duration, macro_steps=int(steps), macro_dt=dt,
        max_depth=max_depth, max_rhs_calls=max_rhs_calls, rhs_calls=0,
        completed_macro_steps=0, accepted_substeps=[], rejected_attempts=[],
        max_depth_reached=0, last_accepted_time=0., minimum_accepted_dt=None,
        maximum_accepted_dt=None, no_clipping_or_projection=True,
        scope='Corrective experiment only; positivity checks are not temporal-error or model admission.')

    def emit(event, **data):
        if callback is not None:
            callback(_jsonable(dict(event=event, rhs_calls=stats['rhs_calls'], **data)))

    def fail_limit(message):
        stats['status'] = 'INCOMPLETE_LIMIT'
        stats['reason'] = message
        raise IntegrationLimitError(message, stats)

    def guard(state, stage, time):
        if state.shape != initial.shape or np.any(~np.isfinite(state)):
            raise FloatingPointError('Trial state is nonfinite or has changed shape at '+stage)
        try:
            result = validate(state)
        except SupportViolation as exc:
            raise _RejectedStage(stage, time, exc) from exc
        if result is not None:
            raise TypeError('validate(state) must return None or raise SupportViolation')

    def evaluate(time, state):
        if stats['rhs_calls'] >= max_rhs_calls:
            fail_limit('Maximum RHS call budget reached before another evaluation')
        stats['rhs_calls'] += 1
        value = np.asarray(rhs(float(time), state), float)
        if value.shape != initial.shape or np.any(~np.isfinite(value)):
            raise FloatingPointError('RHS is nonfinite or has an incorrect shape')
        return value

    def attempt(time, y, h):
        guard(y, 'k1', time)
        k1 = evaluate(time, y)
        y2 = y+h*k1/2
        guard(y2, 'k2', time+h/2)
        k2 = evaluate(time+h/2, y2)
        y3 = y+h*k2/2
        guard(y3, 'k3', time+h/2)
        k3 = evaluate(time+h/2, y3)
        y4 = y+h*k3
        guard(y4, 'k4', time+h)
        k4 = evaluate(time+h, y4)
        candidate = y+h*(k1+2*k2+2*k3+k4)/6
        guard(candidate, 'endpoint', time+h)
        return candidate

    def advance(time, y, h, depth, macro):
        stats['max_depth_reached'] = max(stats['max_depth_reached'], depth)
        try:
            candidate = attempt(time, y, h)
        except _RejectedStage as exc:
            event = _jsonable(dict(macro_index=macro, start_time=float(time), dt=float(h),
                depth=depth, failed_stage=exc.stage, failed_stage_time=exc.time,
                reason=str(exc.violation), violation=exc.violation.details,
                rhs_calls=stats['rhs_calls']))
            stats['rejected_attempts'].append(event)
            emit('REJECT', **{k:v for k,v in event.items() if k != 'rhs_calls'})
            if depth >= max_depth:
                fail_limit('Maximum local bisection depth reached after '+exc.stage+' support violation')
            half = h/2
            if half == 0 or time+half == time or time+h == time+half:
                fail_limit('Bisection cannot represent two distinct time subintervals')
            left = advance(time, y, half, depth+1, macro)
            return advance(time+half, left, half, depth+1, macro)
        event = dict(macro_index=macro, start_time=float(time), end_time=float(time+h),
                     dt=float(h), depth=depth, rhs_calls=stats['rhs_calls'])
        stats['accepted_substeps'].append(event)
        stats['last_accepted_time'] = float(time+h)
        old_min, old_max = stats['minimum_accepted_dt'], stats['maximum_accepted_dt']
        stats['minimum_accepted_dt'] = float(h) if old_min is None else min(old_min, float(h))
        stats['maximum_accepted_dt'] = float(h) if old_max is None else max(old_max, float(h))
        emit('ACCEPT', **{k:v for k,v in event.items() if k != 'rhs_calls'})
        return candidate

    emit('START', duration=duration, macro_steps=int(steps), macro_dt=dt)
    try:
        # Initial support cannot be repaired by reducing a future timestep.
        validate_initial = validate(states[0])
        if validate_initial is not None:
            raise TypeError('validate(state) must return None or raise SupportViolation')
        for i, time in enumerate(times[:-1]):
            states[i+1] = advance(float(time), states[i], dt, 0, i)
            stats['completed_macro_steps'] = i+1
        stats['status'] = 'COMPLETED_NOT_ADJUDICATED'
        emit('STOP', status=stats['status'], completed_macro_steps=steps)
        return times, states, stats
    except BaseException as exc:
        if stats['status'] == 'RUNNING':
            stats['status'] = 'FAILED'
            stats['reason'] = str(exc)
            stats['exception'] = type(exc).__name__
        emit('STOP', status=stats['status'], completed_macro_steps=stats['completed_macro_steps'])
        raise
