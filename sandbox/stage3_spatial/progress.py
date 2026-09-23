"""Checkpoint-driven terminal progress without worker threads or dependencies.

``advance`` confirms completed work; ``checkpoint`` reports an active task without
changing its completion count. Seed ``completed`` with verified reused tasks.
Their historical execution time is not used to estimate this session's speed.

Weights describe relative expected task cost. With ``total_weight`` supplied,
pass each task's weight to ``start_task`` (or explicitly to ``advance``). Without
weights the estimate assumes equally costly units. An explicit task estimate
alone does not estimate the duration of unknown heterogeneous future tasks.
"""
from __future__ import annotations

import math
import sys
import time
from typing import Callable, TextIO


def _nonnegative(value: float, name: str) -> float:
    result = float(value)
    if not math.isfinite(result) or result < 0:
        raise ValueError(f'{name} must be finite and nonnegative')
    return result


def _count(value: int, name: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value < 0:
        raise ValueError(f'{name} must be a nonnegative integer')
    return value


def _duration(seconds: float | None) -> str:
    if seconds is None:
        return 'unknown'
    if seconds < 60:
        return f'{seconds:.1f}s'
    minutes, seconds = divmod(int(seconds), 60)
    if minutes < 60:
        return f'{minutes}m{seconds:02d}s'
    hours, minutes = divmod(minutes, 60)
    return f'{hours}h{minutes:02d}m{seconds:02d}s'


class Progress:
    """Report completed units and honest, approximate time estimates.

    ``callback`` receives a fresh snapshot dictionary for each emitted line.
    Checkpoints and advances are throttled by ``min_interval``; task changes,
    notes and the final outcome are always flushed. Callers choose checkpoints:
    there is no background heartbeat and no per-RHS instrumentation.
    """

    def __init__(self, total: int, label: str, stream: TextIO | None = None,
                 clock: Callable[[], float] = time.monotonic, *,
                 total_weight: float | None = None, completed: int = 0,
                 completed_weight: float = 0, min_interval: float = 1.0,
                 callback: Callable[[dict], None] | None = None):
        self.total = _count(total, 'total')
        self.completed = _count(completed, 'completed')
        if self.completed > self.total:
            raise ValueError('completed cannot exceed total')
        self.label = str(label)
        self.stream = sys.stdout if stream is None else stream
        self.clock = clock
        self.callback = callback
        self.min_interval = _nonnegative(min_interval, 'min_interval')
        self.weighted = total_weight is not None
        self.total_weight = (_nonnegative(total_weight, 'total_weight')
                             if self.weighted else float(self.total))
        if self.total and self.total_weight == 0:
            raise ValueError('nonempty work requires positive total_weight')
        if not self.total and self.total_weight != 0:
            raise ValueError('empty work requires zero total_weight')
        initial_weight = _nonnegative(completed_weight, 'completed_weight')
        if not self.weighted and initial_weight != 0:
            raise ValueError('completed_weight requires total_weight')
        self.completed_weight = initial_weight if self.weighted else float(completed)
        if self.completed_weight > self.total_weight:
            raise ValueError('completed_weight cannot exceed total_weight')
        if self.weighted and ((completed == 0 and initial_weight != 0)
                or (completed > 0 and initial_weight == 0)
                or (completed == total and initial_weight != self.total_weight)
                or (completed < total and initial_weight >= self.total_weight)):
            raise ValueError('resumed task count and weight disagree')
        self.resumed = completed
        self.newly_completed = 0
        self.work_elapsed = 0.0
        self._new_weight = 0.0
        self._task = None
        self._task_started = None
        self._task_weight = 1.0
        self._task_estimate = None
        self._status = 'RUNNING'
        self._last_time = None
        self._started = self._now()
        self._boundary = self._started
        self._last_emitted = None
        self._emit(self._started, force=True,
                   message=f'reused={completed}' if completed else None)

    def _now(self) -> float:
        now = float(self.clock())
        if not math.isfinite(now):
            raise ValueError('clock must return a finite monotonic value')
        if self._last_time is not None and now < self._last_time:
            raise ValueError('clock moved backwards')
        self._last_time = now
        return now

    def _running(self):
        if self._status != 'RUNNING':
            raise RuntimeError('progress has already finished')

    @property
    def current_task(self) -> str | None:
        return self._task

    def start_task(self, label: str, estimated_seconds: float | None = None,
                   *, weight: float = 1.0) -> dict:
        self._running()
        if self.completed == self.total:
            raise ValueError('all registered work is already completed')
        if self._task_started is not None:
            raise ValueError('complete the active task before starting another')
        weight = _nonnegative(weight, 'weight')
        if weight == 0:
            raise ValueError('task weight must be positive')
        if not self.weighted and weight != 1:
            raise ValueError('task weights require total_weight')
        remaining = self.total_weight-self.completed_weight
        if self.weighted and weight > remaining and not math.isclose(weight, remaining, rel_tol=1e-12):
            raise ValueError('task weight exceeds remaining registered weight')
        estimate = (None if estimated_seconds is None else
                    _nonnegative(estimated_seconds, 'estimated_seconds'))
        now = self._now()
        self._task, self._task_started = str(label), now
        self._task_weight, self._task_estimate = weight, estimate
        return self._emit(now, force=True)

    def advance(self, count: int = 1, *, weight: float | None = None) -> dict:
        self._running()
        count = _count(count, 'count')
        if not count or self.completed+count > self.total:
            raise ValueError('advance must complete positive registered work only')
        if not self.weighted and weight is not None:
            raise ValueError('explicit weights require total_weight')
        amount = (self._task_weight*count if weight is None and self.weighted
                  else float(count) if weight is None else _nonnegative(weight, 'weight'))
        new_weight = self.completed_weight+amount
        if amount <= 0 or (new_weight > self.total_weight and not math.isclose(
                new_weight, self.total_weight, rel_tol=1e-12)):
            raise ValueError('advance weight must be positive and within total_weight')
        if self.completed+count == self.total and not math.isclose(
                new_weight, self.total_weight, rel_tol=1e-12, abs_tol=0):
            raise ValueError('final completed task weight does not match total_weight')
        if self.completed+count == self.total:
            new_weight = self.total_weight  # bookkeeping roundoff, not a physical state
        now = self._now()
        interval_start = self._boundary if self._task_started is None else self._task_started
        self.work_elapsed += now-interval_start
        self.newly_completed += count
        self._new_weight += amount
        self.completed += count
        self.completed_weight = new_weight
        completed_task = self._task
        self._task = self._task_started = self._task_estimate = None
        self._task_weight = 1.0
        self._boundary = now
        return self._emit(now, force=self.completed == self.total,
                          message=f'completed: {completed_task}' if completed_task else None)

    def _snapshot(self, now: float, message: str | None = None) -> dict:
        elapsed = now-self._started
        active_elapsed = None if self._task_started is None else now-self._task_started
        remaining_weight = max(0., self.total_weight-self.completed_weight)
        rate = (self.work_elapsed/self._new_weight
                if self._new_weight > 0 and self.work_elapsed > 0 else None)
        task_eta = None
        if active_elapsed is not None:
            estimated_total = self._task_estimate
            if estimated_total is None and rate is not None:
                estimated_total = rate*self._task_weight
            if estimated_total is not None and estimated_total > active_elapsed:
                task_eta = estimated_total-active_elapsed
        eta = None
        if self._status == 'SUCCESS':
            eta = 0.0
        elif self._status == 'RUNNING' and self.completed < self.total:
            if self._task_estimate is not None:
                future_weight = max(0., remaining_weight-self._task_weight)
                if task_eta is not None:
                    if future_weight == 0:
                        eta = task_eta
                    elif rate is not None:
                        eta = task_eta+rate*future_weight
                    elif self.weighted:
                        eta = task_eta+self._task_estimate/self._task_weight*future_weight
            elif rate is not None:
                estimated = rate*remaining_weight-(active_elapsed or 0.)
                eta = estimated if estimated > 0 else None
            elif self.total-self.completed == 1:
                eta = task_eta
        if self._status == 'FAIL':
            percent = None
            eta = task_eta = None
        elif self._status == 'SUCCESS':
            percent = 100.0
        else:
            fraction = self.completed_weight/self.total_weight if self.total_weight else 0.
            percent = min(99.9, 100*fraction)
        return dict(status=self._status, label=self.label, completed=self.completed,
                    total=self.total, percent=percent, elapsed=elapsed,
                    estimated_remaining=eta, current_task=self._task,
                    task_elapsed=active_elapsed, estimated_task_remaining=task_eta,
                    completed_weight=self.completed_weight, total_weight=self.total_weight,
                    resumed=self.resumed, newly_completed=self.newly_completed,
                    work_elapsed=self.work_elapsed, message=message)

    def _emit(self, now: float, *, force: bool = False,
              message: str | None = None) -> dict:
        state = self._snapshot(now, message)
        if not force and self._last_emitted is not None and now-self._last_emitted < self.min_interval:
            return state
        prefix = f'[{state["status"]}] {self.label}'
        if state['percent'] is not None:
            filled = min(20, int(state['percent']/5))
            bar = '#'*filled+'-'*(20-filled)
            prefix += f' [{bar}] {state["percent"]:.1f}%'
        line = f'{prefix} completed={self.completed}/{self.total} elapsed={_duration(state["elapsed"])}'
        if self._status != 'FAIL':
            line += f' ETA={_duration(state["estimated_remaining"])}'
        if self._task is not None:
            line += f' task={self._task} task_elapsed={_duration(state["task_elapsed"])}'
            if self._status != 'FAIL':
                line += f' task_ETA={_duration(state["estimated_task_remaining"])}'
        if message:
            line += f' | {message}'
        print(line, file=self.stream, flush=True)
        self._last_emitted = now
        if self.callback is not None:
            self.callback(dict(state))
        return state

    def checkpoint(self, message: str | None = None, *, force: bool = False) -> dict:
        """Explicit heartbeat from a caller's checkpoint; never marks work done."""
        self._running()
        return self._emit(self._now(), force=force, message=message)

    def note(self, message: str) -> dict:
        self._running()
        return self._emit(self._now(), force=True, message=str(message))

    def finish(self, success: bool = True) -> dict:
        """Finalize once; incomplete success requests emit FAIL and raise ValueError."""
        desired = 'SUCCESS' if success else 'FAIL'
        if self._status != 'RUNNING':
            if self._status != desired:
                raise RuntimeError('cannot change the final outcome')
            return self._snapshot(self._last_time)
        now = self._now()
        incomplete = success and self.completed != self.total
        self._status = 'FAIL' if incomplete else desired
        state = self._emit(now, force=True,
                           message='incomplete work; success refused' if incomplete else None)
        if incomplete:
            raise ValueError('cannot report success before every registered unit completes')
        return state
