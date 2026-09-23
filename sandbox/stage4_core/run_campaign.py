"""Manual stage4A diagnostics, with immutable task outputs and progress/ETA.

Default is a dry run. --execute is deliberate; no background work, retry or
fallback. A negative D.36 sign is a diagnostic result, never a time trajectory.
"""
from pathlib import Path
import argparse
import hashlib
import json
import os
import platform
import sys
import time
import traceback

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from controls import catalogue, local_case, spatial_case, material_reference


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write(path, value):
    if path.exists():
        raise FileExistsError('Preserve previous result: '+str(path))
    path.write_text(json.dumps(value, indent=2, ensure_ascii=False, allow_nan=False)+'\n', encoding='utf8')


def sources():
    names = ['energy_catalog.py','refined_cells.py','cell_closures.py','spatial_functional.py',
             'spatial_nodal.py','spatial_open.py','mixed_spatial.py','rectangular_spatial.py','spatial_dynamics.py']
    paths = [ROOT/'pysnspd/experimental'/name for name in names]
    paths += [Path(__file__), Path(__file__).with_name('controls.py')]
    return {p.relative_to(ROOT).as_posix(): sha(p) for p in paths}


def seconds(value):
    if value is None:
        return 'estimando'
    return f'{int(value)//3600:02d}:{int(value)%3600//60:02d}:{int(value)%60:02d}'


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--plan', type=Path, required=True)
    parser.add_argument('--output-root', type=Path, required=True)
    parser.add_argument('--phase', choices=['all','local','spatial','pilot'], default='all')
    parser.add_argument('--execute', action='store_true')
    args = parser.parse_args()
    plan = json.loads(args.plan.read_text(encoding='utf8'))
    if plan['schema'] != 'pysnspd.stage4A.controls.v1' or plan['photon'] is not False:
        raise ValueError('Unrecognized non-photon control plan')
    tasks = (plan['pilot_cases'] if args.phase == 'pilot' else
             [c for c in plan['cases'] if args.phase == 'all' or c['phase'] == args.phase])
    if not tasks or len({c['id'] for c in tasks}) != len(tasks):
        raise ValueError('Task IDs must be unique and nonempty')
    identity = dict(plan_sha256=sha(args.plan), sources=sources(), phase=args.phase,
                    task_ids=[c['id'] for c in tasks], material=material_reference(plan))
    if not args.execute:
        print(json.dumps(dict(status='DRY_RUN_NO_PHYSICS', **identity, tasks=len(tasks),
            expected_scope='Static and instantaneous only; no photon, circuit or time integrator'), indent=2))
        return
    if args.output_root.exists():
        raise FileExistsError('Choose a NEW output root; previous runs are never overwritten: '+str(args.output_root))
    args.output_root.mkdir(parents=True)
    write(args.output_root/'identity.json', dict(**identity, python=platform.python_version(),
          platform=platform.platform(), threads={k:os.environ.get(k) for k in
          ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']}))
    start = time.monotonic()
    records = []
    active = None
    log = (args.output_root/'progress.jsonl').open('x', encoding='utf8', buffering=1)
    def event(name, **data):
        value = dict(event=name, elapsed_seconds=time.monotonic()-start, **data)
        log.write(json.dumps(value, ensure_ascii=False, allow_nan=False)+'\n')
        print(json.dumps(value, ensure_ascii=False, allow_nan=False), flush=True)
    event('START', tasks=len(tasks), phase=args.phase)
    try:
        for index, case in enumerate(tasks):
            active = case['id']
            output = args.output_root/active
            output.mkdir()
            task_start = time.monotonic()
            stage_start = task_start
            stage = None
            last_print = 0.
            def progress(done, total, label):
                nonlocal stage, stage_start, last_print
                now = time.monotonic()
                if label != stage:
                    stage, stage_start, last_print = label, now, 0.
                if done != total and now-last_print < 5:
                    return
                last_print = now
                fraction = done/total
                remaining = (now-stage_start)*(total-done)/done if done > 1 else None
                bar = '#'*int(24*fraction)+'-'*(24-int(24*fraction))
                print(f'[{bar}] {index+1}/{len(tasks)} {active} | {label} {done}/{total} | transcurrido {seconds(now-start)} | ETA sección {seconds(remaining)}', flush=True)
                log.write(json.dumps(dict(event='PROGRESS',task=active,stage=label,done=done,total=total,
                    elapsed_seconds=now-start,section_eta_seconds=remaining))+'\n')
            prior_elapsed = time.monotonic()-start
            eta = prior_elapsed/index*(len(tasks)-index) if index else None
            event('TASK_START', index=index+1, total=len(tasks), task=active,
                  completed_average_eta_seconds=eta,
                  eta_note='Prior task average is approximate; local and spatial costs differ')
            cat = catalogue(plan, case['resolution'])
            if case['phase'] == 'local':
                result = local_case(plan, case, cat, progress)
            else:
                result = spatial_case(plan, case, cat, progress, output)
            elapsed = time.monotonic()-task_start
            result.update(case=case, runtime_seconds=elapsed, catalogue_diagnostics=cat.diagnostics())
            write(output/'result.json', result)
            records.append(dict(id=active, runtime_seconds=elapsed,
                result=active+'/result.json', sha256=sha(output/'result.json'),
                files=[dict(path=p.relative_to(args.output_root).as_posix(), sha256=sha(p), bytes=p.stat().st_size)
                       for p in sorted(output.iterdir()) if p.is_file()]))
            event('TASK_COMPLETE', task=active, runtime_seconds=elapsed,
                  diagnostics=cat.diagnostics())
        write(args.output_root/'summary.json', dict(status='COMPLETED_DIAGNOSTICS_NOT_PHYSICAL_ADMISSION',
            cases=records, runtime_seconds=time.monotonic()-start, stage4_complete=False,
            photon=False, time_trajectories=0, production_promotion=False))
        event('COMPLETE', tasks=len(records), summary=str(args.output_root/'summary.json'))
    except BaseException as exc:
        value = dict(status='STOPPED', task=active, reason=str(exc), exception=type(exc).__name__,
                     completed_cases=records, runtime_seconds=time.monotonic()-start,
                     policy='No retry, clipping, fallback, overwrite or automatic continuation')
        write(args.output_root/'failure.json', value)
        event('STOPPED', task=active, reason=str(exc))
        traceback.print_exc()
        raise
    finally:
        log.close()


if __name__ == '__main__':
    main()
