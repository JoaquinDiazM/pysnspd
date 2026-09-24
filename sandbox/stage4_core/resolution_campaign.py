"""User-run moment-resolution campaign on the existing thermal cores.

One shared-budget spectral pool evaluates each distinct case/eta/energy once.
Frozen kinetic and potential-projection workers then reuse those maps. Nested
energy grids compare integrated moments without duplicate spectral solves.
No thermal gap relaxation, material fitting or detector transient is performed.
Without --execute this program validates inputs/resources and writes nothing.
"""
from __future__ import annotations
import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
import hashlib
import json
import multiprocessing as mp
import os
from pathlib import Path
import shutil
import sys
import time

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(HERE))
from parallel_runtime import initialize_affinity, limit_thread_environment, linux_resources, resource_budget
limit_thread_environment()
sys.path.insert(0, str(ROOT))
import numpy as np
import next_retarded_oracle as spectral_backend
import next_kinetic_response as kinetic_backend
import next_potential_projection as projection_backend

DEFAULT_PLAN = ROOT/'docs/implementation/stage4/self_consistent_review_20260924/resolution_campaign/plan.json'


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write(path, value):
    with Path(path).open('x', encoding='utf8') as stream:
        json.dump(value, stream, indent=2, allow_nan=False); stream.write('\n')


def validate_plan(plan):
    if plan['schema'] != 'pysnspd.stage4.moment_resolution.v1':
        raise ValueError('Unknown moment-resolution campaign')
    grids = plan['energy_grids']
    for name, count in (('base31', 31), ('fine50', 50)):
        values = grids[name]
        if len(values) != count or sorted(set(values)) != values or values[0] != 0 or values[-1] != 5:
            raise ValueError('Incorrect ordered energy grid: '+name)
    if not set(grids['base31']) <= set(grids['fine50']):
        raise ValueError('The base energy grid must be nested in the fine grid')
    cases = {case['id']:case for case in plan['cases']}
    seen = set()
    for group in plan['groups']:
        if group['case_id'] not in cases or group['grid'] not in grids or group['eta_relative'] <= 0:
            raise ValueError('Invalid registered case/eta/grid group')
        for energy in grids[group['grid']]:
            identifier = spectral_backend.job_id(cases[group['case_id']], group['eta_relative'], energy)
            if identifier in seen:
                raise ValueError('Duplicate spectral query or ambiguous formatted identifier')
            seen.add(identifier)
    if len(seen) != 181:
        raise ValueError('The registered campaign must contain exactly181 unique queries')
    for name, digest in {**plan['frozen_sources'], **plan['inputs']}.items():
        if sha(ROOT/name) != digest:
            raise ValueError('Frozen source/input changed: '+name)
    return cases


def spectral_jobs(plan, pilot=False):
    cases = validate_plan(plan)
    jobs = []
    for group in plan['groups']:
        case = cases[group['case_id']]
        for energy in plan['energy_grids'][group['grid']]:
            if pilot and not (case['id'] == 'radial_129_N256' and group['eta_relative'] == .01 and energy == 1.):
                continue
            identifier = spectral_backend.job_id(case, group['eta_relative'], energy)
            jobs.append(dict(id=identifier, plan=plan['spectral'], case=case,
                energy_relative=energy, eta_relative=group['eta_relative']))
    return sorted(jobs, key=lambda job:(job['energy_relative'], job['case']['id'], job['eta_relative']))


def phase_budget(phase, jobs, maximum_workers=None):
    resources = linux_resources()
    budget = resource_budget(resources, .9, .9, min(jobs, maximum_workers or jobs))
    # Projection workers retain all 50 sparse operators. Reserve more than the
    # generic 1GiB pool estimate; this is scheduling headroom, not an OS limit.
    reserve = (2 if phase == 'projection' else 1)*1024**3
    available = budget['memory_limit_bytes']-budget['coordinator_reserve_bytes']
    count = min(budget['workers'], available//reserve)
    if count < 1:
        raise RuntimeError('Insufficient reserved memory for this phase')
    budget.update(workers=int(count), total_processes=int(count)+1,
        worker_affinity_cpus=budget['worker_affinity_cpus'][:count], worker_reserve_bytes=reserve,
        estimated_memory_reservation_bytes=budget['coordinator_reserve_bytes']+count*reserve)
    return dict(resources=resources, budget=budget)


def disk_budget(plan, output, pilot=False):
    case_sizes = {}
    for case in plan['cases']:
        with np.load(ROOT/case['fields_path']) as data:
            case_sizes[case['id']] = (len(data['d']), len(data['edges']))
    total = 0
    for group in plan['groups']:
        n, e = case_sizes[group['case_id']]
        count = len(plan['energy_grids'][group['grid']])
        # Uncompressed spectral + two-probe kinetic arrays, ZIP/header margin.
        total += count*(256*n+80*e+65536)
        total += (2 if group['grid'] == 'fine50' else 1)*(400*n+120*e+65536)
    projected = int(max(1024**3, 2*total))
    if pilot:
        n,e = case_sizes['radial_129_N256']; projected = 2*(128*n+16*e+65536)
    parent = output.resolve().parent
    while not parent.exists():
        parent = parent.parent
    free = shutil.disk_usage(parent).free
    if projected > .9*free:
        raise RuntimeError('Conservative raw/checkpoint projection exceeds90% of available output volume')
    return dict(projected_uncompressed_with_margin_bytes=projected, available_bytes=free,
        volume_checked=str(parent), available_fraction_limit=.9,
        note='All completed jobs are saved once. No full thermal checkpoint or recovery transient is generated.')


def projection_task(job):
    record, fields = projection_backend.worker(job)
    record.update(id=job['id'], grid=job['grid'], group_id=job['group_id'])
    return record, fields


def run_phase(name, jobs, worker, output, maximum_workers, event):
    environment = phase_budget(name, len(jobs), maximum_workers)
    budget = environment['budget']
    output.mkdir(); fields_dir = output/'fields'; fields_dir.mkdir()
    write(output/'resource_snapshot.json', environment)
    start, records = time.monotonic(), []
    original_affinity = os.sched_getaffinity(0)
    try:
        event('PHASE_START', phase=name, jobs=len(jobs), **environment)
        os.sched_setaffinity(0, {budget['coordinator_cpu']})
        context = mp.get_context('spawn'); counter = context.Value('i', 0)
        with ProcessPoolExecutor(max_workers=budget['workers'], mp_context=context,
            initializer=initialize_affinity, initargs=(budget['worker_affinity_cpus'], counter)) as pool:
            futures = [pool.submit(worker, job) for job in jobs]
            try:
                for future in as_completed(futures):
                    record, fields = future.result()
                    path = fields_dir/(record['id']+'.npz')
                    np.savez_compressed(path, **fields)
                    record.update(fields_path=path.relative_to(output).as_posix(), fields_sha256=sha(path))
                    write(path.with_suffix('.json'), record); records.append(record)
                    fraction = len(records)/len(jobs); elapsed = time.monotonic()-start
                    eta = elapsed*(1-fraction)/fraction
                    bar = '#'*int(24*fraction)+'-'*(24-int(24*fraction))
                    print(f'[{bar}] {name} {len(records)}/{len(jobs)} | {elapsed/60:.1f} min | ETA {eta/60:.1f} min', flush=True)
                    event('PHASE_PROGRESS', phase=name, complete=len(records), total=len(jobs), eta_seconds=eta)
            except BaseException:
                for future in futures:
                    future.cancel()
                raise
        result = dict(status='PHASE_COMPLETE', phase=name, runtime_seconds=time.monotonic()-start, records=records)
        write(output/'summary.json', result)
        event('PHASE_COMPLETE', phase=name, seconds=result['runtime_seconds'])
        return result
    finally:
        os.sched_setaffinity(0, original_affinity)


def projection_jobs(plan, spectral, kinetic, output):
    cases = {case['id']:case for case in plan['cases']}
    lookup = {row['id']:row for row in kinetic['records']}
    jobs = []
    for group in plan['groups']:
        grids = ('base31','fine50') if group['grid'] == 'fine50' else ('base31',)
        for grid in grids:
            energies = plan['energy_grids'][grid]
            selected = sorted([row for row in spectral['records'] if row['case_id'] == group['case_id']
                and row['eta_relative'] == group['eta_relative'] and row['energy_relative'] in energies],
                key=lambda row:row['energy_relative'])
            if [row['energy_relative'] for row in selected] != energies:
                raise ValueError('Missing registered energy in downstream projection')
            jobs.append(dict(id=group['id']+'_'+grid, grid=grid, group_id=group['id'],
                case=cases[group['case_id']], plan=plan['kinetic'], selected=selected,
                kinetic_records=lookup, spectra_root=str((output/'spectra').resolve()),
                kinetic_root=str((output/'kinetic').resolve()),
                temperature_ratio=plan['spectral']['T_K']/plan['spectral']['Tc_K']))
    return jobs


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--plan', type=Path, default=DEFAULT_PLAN)
    parser.add_argument('--output-root', type=Path, required=True)
    parser.add_argument('--pilot', action='store_true')
    parser.add_argument('--max-workers', type=int)
    parser.add_argument('--execute', action='store_true')
    args = parser.parse_args()
    plan = json.loads(args.plan.read_text())
    jobs = spectral_jobs(plan, args.pilot)
    if args.output_root.exists():
        raise FileExistsError('Output already exists: use a fresh directory. No overwrite or implicit resume.')
    environment = phase_budget('spectra', len(jobs), args.max_workers)
    sources = dict(plan['frozen_sources'])
    for name in ('resolution_campaign.py','resolution_campaign_analysis.py'):
        path = HERE/name; sources[path.relative_to(ROOT).as_posix()] = sha(path)
    identity = dict(plan_sha256=sha(args.plan), sources=sources, inputs=plan['inputs'],
        spectral_queries=len(jobs), kinetic_queries=0 if args.pilot else len(jobs),
        projection_groups=0 if args.pilot else 7, pilot=args.pilot,
        disk_budget=disk_budget(plan,args.output_root,args.pilot), initial_resources=environment,
        production_changed=False, physical_time_steps=0, thermal_relaxations=0,
        purpose='Resolution of integrated force,current,energy flux and full-vs-potential projection; not a physical transient')
    if not args.execute:
        print(json.dumps(dict(status='DRY_RUN_NO_WRITES', **identity), indent=2)); return
    args.output_root.mkdir(parents=True)
    write(args.output_root/'identity.json', identity)
    shutil.copyfile(args.plan, args.output_root/'executed_plan.json')
    start = time.monotonic()
    log = (args.output_root/'progress.jsonl').open('x', encoding='utf8', buffering=1)
    def event(name, **values):
        line = json.dumps(dict(event=name, elapsed_seconds=time.monotonic()-start, **values))
        print(line, flush=True); log.write(line+'\n')
    try:
        spectral = run_phase('spectra', jobs, spectral_backend.worker, args.output_root/'spectra', args.max_workers, event)
        if args.pilot:
            result = dict(status='HARDEST_QUERY_PILOT_COMPLETE', runtime_seconds=time.monotonic()-start,
                spectral_queries=len(jobs), stage4_complete=False, production_changed=False)
        else:
            cases = {case['id']:case for case in plan['cases']}
            kinetic_jobs = [dict(plan=plan['kinetic'], case=cases[row['case_id']], record=row,
                spectra_root=str((args.output_root/'spectra').resolve())) for row in spectral['records']]
            kinetic = run_phase('kinetic', kinetic_jobs, kinetic_backend.worker, args.output_root/'kinetic', args.max_workers, event)
            projected_jobs = projection_jobs(plan, spectral, kinetic, args.output_root)
            projection = run_phase('projection', projected_jobs, projection_task, args.output_root/'projection', args.max_workers, event)
            from resolution_campaign_analysis import analyze
            analysis = analyze(plan, args.output_root, projection)
            write(args.output_root/'moment_comparisons.json', analysis)
            result = dict(status='MOMENT_RESOLUTION_CAMPAIGN_COMPLETE', runtime_seconds=time.monotonic()-start,
                spectral_queries=len(jobs), kinetic_queries=len(kinetic_jobs), projection_groups=len(projected_jobs),
                stage4_complete=False, detector_state_admitted=False, production_changed=False,
                next='Interpret observed resolution changes before choosing the dynamic charge closure')
        write(args.output_root/'summary.json', result); event('COMPLETE', **result)
    except BaseException as exc:
        event('FAILED', error=type(exc).__name__, reason=str(exc))
        write(args.output_root/'failure.json', dict(error=type(exc).__name__, reason=str(exc),
            policy='Stop on error; preserve completed jobs. No fallback, clipping, overwrite or automatic retry.'))
        raise
    finally:
        log.close()


if __name__ == '__main__':
    main()
