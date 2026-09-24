"""Plan-driven parallel subprocess runner; contains no SNSPD physics.

Execution is explicit, output directories are never reused, and independent
cases finish even when a sibling fails. Each allocation includes the child
coordinator plus its numerical workers within one runtime resource budget.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import re
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from sandbox.stage4_core.parallel_runtime import (
    THREAD_VARIABLES, limit_thread_environment, linux_resources, resource_budget,
)

SCHEMA = 'pysnspd.stage4.coupled_campaign.v1'
FINISHED = {'SUCCEEDED', 'FAILED', 'LAUNCH_FAILED', 'MISSING_OUTPUT', 'INTERRUPTED'}


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for chunk in iter(lambda: stream.read(1024*1024), b''):
            digest.update(chunk)
    return digest.hexdigest()


def atomic_json(path, data):
    path = Path(path)
    temporary = path.with_name(path.name+'.tmp')
    with temporary.open('w', encoding='utf8', newline='\n') as stream:
        json.dump(data, stream, indent=2, allow_nan=False)
        stream.write('\n'); stream.flush(); os.fsync(stream.fileno())
    os.replace(temporary, path)
    if hasattr(os, 'O_DIRECTORY'):
        descriptor = os.open(path.parent, os.O_RDONLY | os.O_DIRECTORY)
        try:
            os.fsync(descriptor)
        finally:
            os.close(descriptor)


def _positive_integer(value, name):
    if type(value) is not int or value < 1:
        raise ValueError(name+' must be a positive integer')
    return value


def validate_plan(plan):
    if plan.get('schema') != SCHEMA:
        raise ValueError('Unsupported campaign schema')
    cases = plan.get('cases')
    if not isinstance(cases, list) or not cases:
        raise ValueError('At least one independent case is required')
    names = set()
    for case in cases:
        name = case.get('id', '')
        if not isinstance(name, str) or re.fullmatch(r'[a-z0-9][a-z0-9_-]{0,79}', name) is None:
            raise ValueError('Case ids must be simple lowercase names')
        if name in names:
            raise ValueError('Case ids must be unique')
        names.add(name)
        argv = case.get('argv')
        if (not isinstance(argv, list) or not argv or
                any(not isinstance(item, str) or not item or '\0' in item for item in argv)):
            raise ValueError('Each argv must be a nonempty list of nonempty strings')
        if not any('{workers}' in item for item in argv):
            raise ValueError('Each argv must expose the allocated {workers}')
        if not any('{output}' in item for item in argv):
            raise ValueError('Each argv must expose its distinct {output}')
        estimate = case.get('estimated_seconds')
        if estimate is not None and (isinstance(estimate, bool) or
                not isinstance(estimate, (float, int)) or not math.isfinite(estimate) or estimate <= 0):
            raise ValueError('estimated_seconds must be positive and finite')
        if not isinstance(case.get('expected_outputs', []), list):
            raise ValueError('expected_outputs must be a list')
        for item in case.get('expected_outputs', []):
            if not isinstance(item, str) or not item or Path(item).is_absolute() or '..' in Path(item).parts:
                raise ValueError('Expected outputs must be relative paths within one case')
    _positive_integer(plan.get('maximum_parallel_cases', len(cases)), 'maximum_parallel_cases')
    if 'maximum_workers_per_case' in plan:
        _positive_integer(plan['maximum_workers_per_case'], 'maximum_workers_per_case')
    heartbeat = plan.get('heartbeat_seconds', 5.)
    if isinstance(heartbeat, bool) or not isinstance(heartbeat, (float, int)) or not 0.1 <= heartbeat <= 60:
        raise ValueError('heartbeat_seconds must be between 0.1 and 60')
    return plan


def allocate(plan, resources):
    """Disjoint CPU slots: one campaign coordinator, then case coordinators/workers."""
    validate_plan(plan)
    budget = resource_budget(resources, .9, .9)
    cpus = budget['worker_affinity_cpus']
    parallel = min(len(plan['cases']), plan.get('maximum_parallel_cases', len(plan['cases'])), len(cpus)//2)
    case_coordinator_reserve = 3*1024**3
    while parallel:
        affordable_workers = (budget['memory_limit_bytes']-budget['coordinator_reserve_bytes']-
            parallel*case_coordinator_reserve)//budget['worker_reserve_bytes']
        if affordable_workers >= parallel:
            break
        parallel -= 1
    if parallel < 1:
        raise ValueError('Budget cannot fit a case coordinator and one numerical worker')
    width, extra = divmod(min(len(cpus), parallel+affordable_workers), parallel)
    slots, offset = [], 0
    for index in range(parallel):
        count = width+(index < extra)
        if 'maximum_workers_per_case' in plan:
            count = min(count, plan['maximum_workers_per_case']+1)
        assigned = cpus[offset:offset+count]
        offset += count
        slots.append(dict(slot=index, cpus=assigned, workers=count-1,
                          case_coordinator_cpu=assigned[-1], numerical_worker_cpus=assigned[:-1]))
    return dict(shared_budget=budget, parallel_cases=parallel, slots=slots,
        numerical_worker_ceiling=sum(row['workers'] for row in slots),
        total_process_cpu_ceiling=1+sum(len(row['cpus']) for row in slots),
        case_coordinator_reserve_bytes=case_coordinator_reserve,
        estimated_memory_reservation_bytes=budget['coordinator_reserve_bytes']+
            parallel*case_coordinator_reserve+sum(row['workers'] for row in slots)*budget['worker_reserve_bytes'],
        accounting='Each case slot includes its subprocess coordinator plus workers; BLAS/OpenMP=1. CPU sets are disjoint.')


def exact_argv(case, workers, output):
    substitutions = {'{workers}': str(workers), '{output}': str(Path(output).resolve()),
                     '{python}': sys.executable, '{repo}': str(ROOT)}
    values = []
    for item in case['argv']:
        for key, value in substitutions.items():
            item = item.replace(key, value)
        values.append(item)
    return values


def _launch(argv, cwd, environment, cpus, stream):
    # POSIX preexec is safe here: the coordinator never starts Python threads.
    return subprocess.Popen(argv, cwd=cwd, env=environment, shell=False,
        stdin=subprocess.DEVNULL, stdout=stream, stderr=subprocess.STDOUT,
        start_new_session=True, preexec_fn=lambda: os.sched_setaffinity(0, set(cpus)))


def read_progress(active, final=False):
    """Read bounded increments from the durable console log, retaining JSON telemetry only."""
    stream = active['reader']
    chunk = stream.read(65536)
    if not chunk and not final:
        return
    text = active.get('partial', '')+chunk
    lines = text.split('\n')
    active['partial'] = lines.pop()
    if final and active['partial']:
        lines.append(active.pop('partial'))
        active['partial'] = ''
    if len(active.get('partial', '')) > 1024*1024:
        active['partial'] = ''  # The complete bytes remain in the console log.
    for line in lines:
        try:
            value = json.loads(line)
        except (ValueError, TypeError, RecursionError):
            continue
        if not isinstance(value, dict):
            continue
        progress = active['record'].setdefault('progress', {})
        if isinstance(value.get('event'), str):
            progress['last_event'] = value['event'][:160]
        for key in ('fraction', 'progress_fraction', 'eta_seconds'):
            number = value.get(key)
            if isinstance(number, bool) or not isinstance(number, (int, float)) or not math.isfinite(number):
                continue
            if key == 'eta_seconds' and number >= 0:
                progress['eta_seconds'] = float(number)
            elif key != 'eta_seconds' and 0 <= number <= 1:
                progress['fraction'] = float(number)
        if value.get('horizon_ps', 0) and isinstance(value.get('time_ps'), (int, float)):
            horizon = value['horizon_ps']
            if isinstance(horizon, (int, float)) and math.isfinite(horizon) and horizon > 0:
                fraction = value['time_ps']/horizon
                if math.isfinite(fraction) and 0 <= fraction <= 1:
                    progress['fraction'] = float(fraction)
        progress['read_elapsed_seconds'] = time.monotonic()-active['campaign_started']


def estimated_eta(records, parallel, now):
    remaining, pending = [], []
    for record in records:
        if record['status'] in FINISHED:
            continue
        if record['status'] == 'PENDING':
            if record.get('estimated_seconds') is None:
                return None
            pending.append(record['estimated_seconds']); continue
        progress = record.get('progress', {})
        if 'eta_seconds' in progress:
            since = max(0., now-progress.get('read_elapsed_seconds', now))
            remaining.append(max(0., progress['eta_seconds']-since))
        elif record.get('estimated_seconds') is not None:
            if 'fraction' in progress:
                remaining.append(record['estimated_seconds']*(1-progress['fraction']))
            else:
                remaining.append(max(0., record['estimated_seconds']-(now-record['start_elapsed_seconds'])))
        else:
            return None
    remaining += [0.]*(parallel-len(remaining))
    for estimate in pending:
        index = min(range(len(remaining)), key=remaining.__getitem__)
        remaining[index] += estimate
    return max(remaining, default=0.)


def run_campaign(plan, output_root, *, resources=None):
    """Execute a user-requested campaign; each case is attempted once."""
    validate_plan(plan)
    output_root = Path(output_root).resolve()
    if output_root.exists():
        raise FileExistsError('Campaign output already exists; no reuse, overwrite or implicit resume')
    resources = linux_resources() if resources is None else resources
    allocation = allocate(plan, resources)
    if not hasattr(os, 'sched_setaffinity'):
        raise RuntimeError('Execution requires Linux process affinity')
    cwd = Path(plan.get('working_directory', ROOT)).resolve()
    if not cwd.is_dir():
        raise ValueError('The declared working directory must already exist')
    output_root.mkdir(parents=True)
    (output_root/'logs').mkdir(); (output_root/'case_records').mkdir(); (output_root/'cases').mkdir()
    atomic_json(output_root/'executed_plan.json', plan)
    started = time.monotonic()
    records = [dict(id=case['id'], status='PENDING', attempts=0,
                    estimated_seconds=case.get('estimated_seconds')) for case in plan['cases']]
    manifest = dict(schema=SCHEMA, status='RUNNING', plan_sha256=sha(output_root/'executed_plan.json'),
        resources=resources, allocation=allocation, cases=records,
        numerical_physics='None: this runner only schedules declared argv lists',
        retry_policy='One attempt per case; independent siblings continue after failure')
    environment = os.environ.copy()
    for variable in THREAD_VARIABLES:
        environment[variable] = '1'
    environment['PYTHONUNBUFFERED'] = '1'
    previous_affinity = os.sched_getaffinity(0)
    active, pending = {}, list(range(len(records)))
    event_stream = (output_root/'progress.jsonl').open('x', encoding='utf8', buffering=1)

    def event(kind, **values):
        item = dict(event=kind, elapsed_seconds=time.monotonic()-started, **values)
        line = json.dumps(item, allow_nan=False)
        event_stream.write(line+'\n'); event_stream.flush(); os.fsync(event_stream.fileno())
        print(line, flush=True)

    def persist():
        manifest['elapsed_seconds'] = time.monotonic()-started
        atomic_json(output_root/'campaign.json', manifest)
        for record in records:
            atomic_json(output_root/'case_records'/f'{record["id"]}.json', record)

    def finish(slot_id, returncode):
        item = active.pop(slot_id); record = item['record']
        read_progress(item, final=True)
        item['reader'].close()
        item['log'].flush(); os.fsync(item['log'].fileno()); item['log'].close()
        record.update(returncode=returncode, end_elapsed_seconds=time.monotonic()-started,
                      status='SUCCEEDED' if returncode == 0 else 'FAILED')
        record['duration_seconds'] = record['end_elapsed_seconds']-record['start_elapsed_seconds']
        record['console_sha256'] = sha(output_root/record['console_log'])
        outputs = []
        for relative in item['case'].get('expected_outputs', []):
            path = (Path(record['output'])/relative).resolve()
            case_root = Path(record['output']).resolve()
            if not path.is_relative_to(case_root):
                outputs.append(dict(path=relative, status='OUTSIDE_CASE_DIRECTORY')); continue
            outputs.append(dict(path=relative, status='PRESENT' if path.is_file() else 'MISSING',
                                sha256=sha(path) if path.is_file() else None))
        record['outputs'] = outputs
        if returncode == 0 and any(row['status'] != 'PRESENT' for row in outputs):
            record['status'] = 'MISSING_OUTPUT'
        event('CASE_COMPLETE', id=record['id'], status=record['status'], returncode=returncode,
              duration_seconds=record['duration_seconds'])
        persist()

    try:
        os.sched_setaffinity(0, {allocation['shared_budget']['coordinator_cpu']})
        persist(); event('CAMPAIGN_START', allocation=allocation)
        last_heartbeat = -float('inf')
        while pending or active:
            for slot in allocation['slots']:
                if not pending or slot['slot'] in active:
                    continue
                index = pending.pop(0); case, record = plan['cases'][index], records[index]
                case_output = output_root/'cases'/case['id']
                argv = exact_argv(case, slot['workers'], case_output)
                log_path = output_root/'logs'/f'{case["id"]}.console.log'
                log = log_path.open('xb')
                record.update(status='RUNNING', attempts=1, argv=argv, workers=slot['workers'],
                    affinity_cpus=slot['cpus'], output=str(case_output),
                    console_log=str(log_path.relative_to(output_root)),
                    start_elapsed_seconds=time.monotonic()-started)
                persist()
                try:
                    case_environment = dict(environment)
                    case_environment['PYSNSPD_SHARED_ALLOCATION'] = json.dumps(dict(
                        parent_pid=os.getpid(), case_id=case['id'], cpus=slot['cpus'],
                        workers=slot['workers'], shared_budget=allocation['shared_budget'],
                        total_process_cpu_ceiling=allocation['total_process_cpu_ceiling']))
                    process = _launch(argv, cwd, case_environment, slot['cpus'], log)
                except (OSError, ValueError, subprocess.SubprocessError) as error:
                    log.close()
                    record.update(status='LAUNCH_FAILED', returncode=None, reason=str(error),
                                  end_elapsed_seconds=time.monotonic()-started)
                    event('CASE_LAUNCH_FAILED', id=case['id'], reason=str(error)); persist(); continue
                record['pid'] = process.pid
                active[slot['slot']] = dict(process=process, record=record, case=case,
                    reader=log_path.open('r', encoding='utf8', errors='replace'), log=log,
                    campaign_started=started, partial='')
                event('CASE_START', id=case['id'], argv=argv, workers=slot['workers'], cpus=slot['cpus'])
                persist()
            for slot_id, item in list(active.items()):
                read_progress(item)
                returncode = item['process'].poll()
                if returncode is not None:
                    finish(slot_id, returncode)
            now = time.monotonic()-started
            if now-last_heartbeat >= plan.get('heartbeat_seconds', 5.) or not (pending or active):
                completed = sum(row['status'] in FINISHED for row in records)
                progress = (completed+sum(item['record'].get('progress', {}).get('fraction', 0.)
                                          for item in active.values()))/len(records)
                eta = estimated_eta(records, allocation['parallel_cases'], now)
                event('HEARTBEAT', completed_cases=completed, total_cases=len(records),
                    progress_fraction=progress, eta_seconds=eta,
                    eta_basis='Case-reported ETA or declared estimates; not a guaranteed completion time',
                    active_cases=[item['record']['id'] for item in active.values()],
                    numerical_workers=sum(item['record']['workers'] for item in active.values()),
                    case_coordinators=len(active), campaign_coordinators=1)
                bar = ('#'*int(24*progress)).ljust(24, '-')
                print(f'[{bar}] {completed}/{len(records)} cases | {now:.1f}s | '
                      f'ETA {"unknown" if eta is None else f"{eta:.1f}s"}', flush=True)
                last_heartbeat = now; persist()
            if pending or active:
                time.sleep(.05)
        manifest['status'] = 'SUCCEEDED' if all(row['status'] == 'SUCCEEDED' for row in records) else 'COMPLETED_WITH_FAILURES'
        persist(); event('CAMPAIGN_COMPLETE', status=manifest['status'])
        return manifest
    except BaseException as error:
        manifest.update(status='INTERRUPTED', reason=f'{type(error).__name__}: {error}')
        for item in active.values():
            # Subprocesses have their own sessions; terminate their process groups
            # so nested numerical workers cannot outlive an interrupted campaign.
            try:
                import signal
                os.killpg(item['process'].pid, signal.SIGTERM)
            except ProcessLookupError:
                pass
        for item in active.values():
            try:
                item['process'].wait(timeout=3)
            except subprocess.TimeoutExpired:
                import signal
                os.killpg(item['process'].pid, signal.SIGKILL)
                item['process'].wait(timeout=3)
            item['reader'].close()
            item['log'].flush(); os.fsync(item['log'].fileno()); item['log'].close()
            item['record'].update(status='INTERRUPTED', returncode=item['process'].returncode)
        persist(); event('CAMPAIGN_INTERRUPTED', reason=manifest['reason']); raise
    finally:
        os.sched_setaffinity(0, previous_affinity)
        event_stream.close()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--plan', required=True, type=Path)
    parser.add_argument('--output-root', required=True, type=Path)
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument('--execute', action='store_true')
    mode.add_argument('--dry-run', action='store_true', help='Default: inspect argv and runtime resource allocation only')
    args = parser.parse_args()
    limit_thread_environment()
    plan = validate_plan(json.loads(args.plan.read_text(encoding='utf8')))
    resources = linux_resources()
    if not args.execute:
        allocation = allocate(plan, resources)
        previews = []
        for index, case in enumerate(plan['cases']):
            slot = allocation['slots'][index % allocation['parallel_cases']]
            previews.append(dict(id=case['id'], workers=slot['workers'],
                example_argv=exact_argv(case, slot['workers'], args.output_root/'cases'/case['id'])))
        print(json.dumps(dict(status='DRY_RUN', output_exists=args.output_root.exists(), resources=resources,
            allocation=allocation, cases=previews, note='Slot assignment can change as independent cases finish'), indent=2))
        return 0
    result = run_campaign(plan, args.output_root, resources=resources)
    return 0 if result['status'] == 'SUCCEEDED' else 1


if __name__ == '__main__':
    raise SystemExit(main())
