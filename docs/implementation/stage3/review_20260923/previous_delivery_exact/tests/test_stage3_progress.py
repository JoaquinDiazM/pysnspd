"""Progress accounting uses fake time/streams; no simulation or sleeps."""
import io
import pytest

from sandbox.stage3_spatial.progress import Progress


class Clock:
    def __init__(self):
        self.value = 0.0

    def __call__(self):
        return self.value

    def tick(self, seconds):
        self.value += seconds


class FlushedStream(io.StringIO):
    def __init__(self):
        super().__init__()
        self.flush_count = 0

    def flush(self):
        self.flush_count += 1


def fixture(total=3, **kwargs):
    clock, stream = Clock(), FlushedStream()
    return Progress(total, 'static batch', stream, clock, **kwargs), clock, stream


def test_unknown_initial_eta_and_explicit_heartbeat_do_not_complete_work():
    progress, clock, stream = fixture()
    assert 'ETA=unknown' in stream.getvalue()
    progress.start_task('energy derivative')
    clock.tick(3)
    state = progress.checkpoint('mesh checkpoint')
    assert state['completed'] == 0
    assert state['percent'] == 0
    assert state['elapsed'] == state['task_elapsed'] == 3
    assert state['estimated_remaining'] is None
    assert 'task=energy derivative' in stream.getvalue()
    assert stream.flush_count == 3


def test_measured_eta_excludes_completed_reused_work_and_between_task_overhead():
    progress, clock, _ = fixture(5, completed=2)
    clock.tick(100)  # setup overhead, not two freshly completed tasks
    progress.start_task('first new task')
    clock.tick(10)
    state = progress.advance()
    assert state['elapsed'] == 110
    assert state['work_elapsed'] == 10
    assert state['estimated_remaining'] == 20
    assert state['resumed'] == 2 and state['newly_completed'] == 1
    clock.tick(50)
    state = progress.start_task('second new task')
    assert state['estimated_remaining'] == 20
    clock.tick(4)
    state = progress.checkpoint()
    assert state['estimated_remaining'] == 16
    assert state['completed'] == 3


def test_weighted_heterogeneous_estimate_and_accounting():
    progress, clock, _ = fixture(3, total_weight=12)
    state = progress.start_task('small mesh', estimated_seconds=10, weight=2)
    assert state['estimated_task_remaining'] == 10
    assert state['estimated_remaining'] == 60
    clock.tick(10)
    state = progress.advance()
    assert state['percent'] == pytest.approx(100/6)
    assert state['estimated_remaining'] == 50
    progress.start_task('large mesh', weight=8)
    clock.tick(30)
    state = progress.advance()
    assert state['estimated_remaining'] == 8
    progress.start_task('last mesh', weight=2)
    clock.tick(8)
    progress.advance()
    assert progress.finish()['percent'] == 100


def test_task_estimate_does_not_invent_duration_of_unknown_future_tasks():
    progress, clock, _ = fixture(2)
    state = progress.start_task('known task', estimated_seconds=5)
    assert state['estimated_task_remaining'] == 5
    assert state['estimated_remaining'] is None
    clock.tick(8)
    state = progress.checkpoint()
    assert state['estimated_task_remaining'] is None
    assert state['estimated_remaining'] is None


def test_decimal_task_weights_do_not_fail_due_to_binary_roundoff():
    progress, clock, _ = fixture(2, total_weight=.3)
    for weight in (.1,.2):
        progress.start_task('weighted', weight=weight)
        clock.tick(1)
        progress.advance()
    assert progress.finish()['completed_weight'] == .3


def test_last_task_estimate_available_but_overrun_does_not_claim_zero_eta():
    progress, clock, _ = fixture(1)
    state = progress.start_task('only task', estimated_seconds=2)
    assert state['estimated_remaining'] == 2
    clock.tick(3)
    assert progress.checkpoint()['estimated_remaining'] is None


def test_explicit_task_estimate_takes_precedence_over_prior_throughput():
    progress, clock, _ = fixture(2)
    progress.start_task('fast task')
    clock.tick(2)
    progress.advance()
    state = progress.start_task('known slow task', estimated_seconds=100)
    assert state['estimated_remaining'] == 100
    clock.tick(101)
    assert progress.checkpoint()['estimated_remaining'] is None


def test_starting_another_task_cannot_hide_unfinished_work():
    progress, _, _ = fixture()
    progress.start_task('unfinished')
    with pytest.raises(ValueError, match='active task'):
        progress.start_task('replacement')
    assert progress.current_task == 'unfinished'


def test_notes_and_redirected_output_flushed_without_carriage_return_spam():
    snapshots = []
    progress, clock, stream = fixture(callback=snapshots.append, min_interval=5)
    progress.start_task('counted operation')
    lines = len(stream.getvalue().splitlines())
    for _ in range(100):
        progress.checkpoint()
    assert len(stream.getvalue().splitlines()) == lines
    progress.note('reading already computed data')
    assert len(stream.getvalue().splitlines()) == lines+1
    clock.tick(5)
    progress.checkpoint()
    assert len(snapshots) == lines+2
    assert stream.flush_count == len(snapshots)
    assert '\r' not in stream.getvalue()


def test_success_only_finalizes_at_finish_and_is_idempotent():
    progress, clock, stream = fixture(1)
    progress.start_task('operation')
    clock.tick(1)
    state = progress.advance()
    assert state['percent'] == 99.9
    assert '100.0%' not in stream.getvalue()
    state = progress.finish()
    assert state['status'] == 'SUCCESS' and state['estimated_remaining'] == 0
    final_text = stream.getvalue()
    assert '100.0%' in final_text.splitlines()[-1]
    progress.finish()
    assert stream.getvalue() == final_text
    with pytest.raises(RuntimeError):
        progress.advance()
    with pytest.raises(RuntimeError):
        progress.finish(False)


@pytest.mark.parametrize('complete_before_failure', [False, True])
def test_failure_never_prints_full_completion_or_a_false_eta(complete_before_failure):
    progress, clock, stream = fixture(1)
    progress.start_task('operation')
    clock.tick(1)
    if complete_before_failure:
        progress.advance()
    state = progress.finish(False)
    final_line = stream.getvalue().splitlines()[-1]
    assert state['status'] == 'FAIL'
    assert state['percent'] is None and state['estimated_remaining'] is None
    assert '[FAIL]' in final_line and 'ETA=' not in final_line
    assert '100.0%' not in stream.getvalue()


def test_incomplete_success_request_reports_failure_then_raises():
    progress, _, stream = fixture()
    with pytest.raises(ValueError, match='every registered unit'):
        progress.finish()
    assert '[FAIL]' in stream.getvalue().splitlines()[-1]
    assert '100.0%' not in stream.getvalue()


def test_fully_reused_tasks_have_no_invented_throughput():
    progress, _, _ = fixture(2, completed=2)
    state = progress.finish()
    assert state['resumed'] == 2 and state['newly_completed'] == 0
    assert state['work_elapsed'] == 0
    assert state['status'] == 'SUCCESS'


def test_empty_batch_can_finish_without_fake_work():
    progress, _, _ = fixture(0)
    assert progress.finish()['status'] == 'SUCCESS'


@pytest.mark.parametrize('kwargs', [
    {'total': -1}, {'total': True}, {'total': 2, 'completed': 3},
    {'total': 2, 'total_weight': 0}, {'total': 0, 'total_weight': 1},
    {'total': 2, 'completed_weight': 1},
    {'total': 2, 'total_weight': 3, 'completed': 1},
    {'total': 2, 'total_weight': 3, 'completed': 1, 'completed_weight': 3},
    {'total': 2, 'total_weight': float('inf')},
    {'total': 2, 'min_interval': -1},
])
def test_invalid_configuration_rejected(kwargs):
    with pytest.raises(ValueError):
        Progress(label='invalid', stream=io.StringIO(), **kwargs)


def test_invalid_advance_does_not_mutate_completion():
    progress, _, _ = fixture(2, total_weight=3)
    with pytest.raises(ValueError):
        progress.advance(0)
    with pytest.raises(ValueError):
        progress.advance(3)
    with pytest.raises(ValueError):
        progress.advance(2, weight=2)
    assert progress.completed == 0 and progress.completed_weight == 0


def test_non_monotonic_clock_is_rejected():
    progress, clock, _ = fixture()
    clock.tick(-1)
    with pytest.raises(ValueError, match='backwards'):
        progress.checkpoint()


def test_default_stream_is_resolved_at_construction(monkeypatch):
    stream = FlushedStream()
    monkeypatch.setattr('sys.stdout', stream)
    progress = Progress(0, 'stdout', clock=Clock())
    progress.finish()
    assert stream.flush_count == 2


def test_callback_errors_propagate_instead_of_being_hidden():
    def broken_callback(_):
        raise RuntimeError('callback failure')
    with pytest.raises(RuntimeError, match='callback failure'):
        fixture(callback=broken_callback)
