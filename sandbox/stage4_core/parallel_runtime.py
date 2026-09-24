"""Parallel exact-spectrum execution of independent static stage4 controls.

Only the expensive independent spectral queries are parallelized. Request
collection and the final small spatial assemblies run in the coordinator.
Both cases share one work queue, so a completed case releases its capacity.
No interpolation, rounded cache keys, or implicit serial fallback is used.
"""
from __future__ import annotations

import argparse
from collections import deque
from concurrent.futures import ProcessPoolExecutor, FIRST_COMPLETED, wait
import json
import math
import multiprocessing as mp
import os
from pathlib import Path
import sys
import time


ROOT = Path(__file__).resolve().parents[2]
THREAD_VARIABLES = ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS',
                    'BLIS_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'NUMEXPR_NUM_THREADS')
_WORKER_PLAN = None
_WORKER_CATALOGUES = {}


def cgroup_limits(cgroup_file=Path('/proc/self/cgroup'), cgroup_root=Path('/sys/fs/cgroup')):
    """Read every visible ancestor limit; missing files are explicitly reported.

    v2 uses its unified mount. v1 recognizes conventional controller mounts;
    a nonstandard unresolved mount is reported as unavailable, not unlimited.
    """
    metadata = dict(cpu_records=[], memory_records=[], unavailable=[])
    cpu_limits, memory_headroom = [], []
    try:
        lines = Path(cgroup_file).read_text().splitlines()
    except (FileNotFoundError, PermissionError) as exc:
        metadata['unavailable'].append(type(exc).__name__+': '+str(cgroup_file))
        return dict(cpu_quota_units=None, memory_headroom_bytes=None, **metadata)

    def ancestors(root, suffix):
        parts = Path(suffix.lstrip('/')).parts
        if '..' in parts:
            raise ValueError('Invalid parent traversal in cgroup path')
        here = root.joinpath(*parts)
        while True:
            yield here
            if here == root:
                break
            here = here.parent

    def readable(path):
        try:
            return path.read_text().strip()
        except FileNotFoundError:
            return None
        except PermissionError:
            metadata['unavailable'].append('Unreadable: '+str(path))
            return None

    def cpu_record(directory, version):
        path = directory/('cpu.max' if version == 2 else 'cpu.cfs_quota_us')
        value = readable(path)
        if value is None:
            return
        try:
            if version == 2:
                quota, period = value.split()
                period = int(period)
                if period <= 0:
                    raise ValueError('nonpositive period')
                amount = None if quota == 'max' else int(quota)/period
            else:
                quota = int(value)
                if quota == -1:
                    amount = None
                else:
                    period = int(readable(directory/'cpu.cfs_period_us'))
                    if period <= 0:
                        raise ValueError('nonpositive period')
                    amount = quota/period
            if amount is not None and (not math.isfinite(amount) or amount <= 0):
                raise ValueError('nonpositive/nonfinite quota')
        except (ValueError, TypeError) as exc:
            raise ValueError('Invalid readable cgroup CPU limit at '+str(path)) from exc
        metadata['cpu_records'].append(dict(path=str(path), version=version,
            cpu_units=amount, status='unlimited' if amount is None else 'finite'))
        if amount is not None:
            cpu_limits.append(amount)

    def memory_record(directory, version):
        path = directory/('memory.max' if version == 2 else 'memory.limit_in_bytes')
        value = readable(path)
        if value is None:
            return
        if version == 2 and value == 'max':
            metadata['memory_records'].append(dict(path=str(path), version=version, status='unlimited'))
            return
        current_path = directory/('memory.current' if version == 2 else 'memory.usage_in_bytes')
        try:
            maximum, current = int(value), int(readable(current_path))
            if maximum < 0 or current < 0:
                raise ValueError('negative memory value')
        except (ValueError, TypeError) as exc:
            raise ValueError('Invalid/incomplete finite cgroup memory limit at '+str(path)) from exc
        headroom = max(0, maximum-current)
        metadata['memory_records'].append(dict(path=str(path), version=version,
            status='finite', maximum_bytes=maximum, current_bytes=current, headroom_bytes=headroom))
        memory_headroom.append(headroom)

    root = Path(cgroup_root)
    for line in lines:
        try:
            hierarchy, controllers, suffix = line.split(':', 2)
        except ValueError as exc:
            raise ValueError('Invalid readable /proc/self/cgroup record') from exc
        if hierarchy == '0' and not controllers:
            for directory in ancestors(root, suffix):
                cpu_record(directory, 2)
                memory_record(directory, 2)
        else:
            names = controllers.split(',')
            for controller, callback in (('cpu', cpu_record), ('memory', memory_record)):
                if controller not in names:
                    continue
                candidates = [root/controllers, root/controller]
                if controller == 'cpu':
                    candidates += [root/'cpu,cpuacct', root/'cpuacct,cpu']
                mount = next((candidate for candidate in candidates if candidate.is_dir()), None)
                if mount is None:
                    metadata['unavailable'].append('Unresolved v1 controller mount: '+controllers)
                    continue
                for directory in ancestors(mount, suffix):
                    callback(directory, 1)
    for kind in ('cpu', 'memory'):
        if not metadata[kind+'_records']:
            metadata['unavailable'].append('No readable '+kind+' limits in visible cgroup ancestors')
    return dict(cpu_quota_units=min(cpu_limits) if cpu_limits else None,
        memory_headroom_bytes=min(memory_headroom) if memory_headroom else None, **metadata)


def linux_resources():
    """Read the actual affinity and physical topology; no guessed SMT factor."""
    if not hasattr(os, 'sched_getaffinity'):
        raise RuntimeError('This runner requires Linux affinity/topology; use serial runner on this host')
    allowed = set(os.sched_getaffinity(0))
    groups = {}
    for cpu in sorted(allowed):
        base = Path(f'/sys/devices/system/cpu/cpu{cpu}/topology')
        package = int((base/'physical_package_id').read_text())
        core = int((base/'core_id').read_text())
        siblings = set()
        for piece in (base/'thread_siblings_list').read_text().strip().split(','):
            ends = piece.split('-')
            siblings.update(range(int(ends[0]), int(ends[-1])+1))
        # Reserve/allocate only physical cores whose complete sibling set is available.
        if siblings <= allowed:
            groups[(package, core)] = sorted(siblings)
    mem = {line.split(':')[0]: int(line.split()[1])*1024
           for line in Path('/proc/meminfo').read_text().splitlines()
           if line.startswith(('MemTotal:', 'MemAvailable:'))}
    limits = cgroup_limits()
    available = mem['MemAvailable']
    if limits['memory_headroom_bytes'] is not None:
        available = min(available, limits['memory_headroom_bytes'])
    return dict(logical_cpus=sorted(allowed), physical_core_groups=list(groups.values()),
                available_memory_bytes=available, host_available_memory_bytes=mem['MemAvailable'],
                total_memory_bytes=mem['MemTotal'], cpu_quota_units=limits['cpu_quota_units'],
                cgroup_limits=limits)


def resource_budget(resources, cpu_fraction=.9, memory_fraction=.9, max_workers=None):
    """Bound coordinator + workers, reserving whole physical cores and RAM."""
    if not 0 < cpu_fraction <= .9 or not 0 < memory_fraction <= .9:
        raise ValueError('CPU and memory fractions must be in (0,0.9]')
    cores = resources['physical_core_groups']
    core_count = math.floor(cpu_fraction*len(cores))
    effective_cpu_units = len(resources['logical_cpus'])
    if resources.get('cpu_quota_units') is not None:
        effective_cpu_units = min(effective_cpu_units, resources['cpu_quota_units'])
    logical_limit = math.floor(cpu_fraction*effective_cpu_units)
    selected = []
    selected_cores = 0
    for group in cores[:core_count]:
        if len(selected)+len(group) <= logical_limit:
            selected.extend(group)
            selected_cores += 1
    if len(selected) < 2:
        raise ValueError('Budget cannot fit a coordinator and worker on whole physical cores; no silent oversubscription')
    memory_limit = math.floor(memory_fraction*resources['available_memory_bytes'])
    coordinator_reserve = 3*1024**3
    worker_reserve = 1024**3
    workers = min(len(selected)-1, (memory_limit-coordinator_reserve)//worker_reserve)
    if max_workers is not None:
        if type(max_workers) is not int or max_workers < 1:
            raise ValueError('max_workers must be a positive integer')
        workers = min(workers, max_workers)
    if workers < 1:
        raise ValueError('Available RAM cannot fit the conservative coordinator/worker reservations')
    return dict(cpu_fraction=cpu_fraction, memory_fraction=memory_fraction,
        selected_affinity_cpus=selected, coordinator_cpu=selected[-1],
        worker_affinity_cpus=selected[:workers], workers=workers, threads_per_process=1,
        total_processes=workers+1, logical_budget=logical_limit,
        effective_cpu_quota_units=effective_cpu_units,
        physical_cores_used=selected_cores,
        physical_cores_reserved=len(cores)-selected_cores,
        memory_limit_bytes=memory_limit, coordinator_reserve_bytes=coordinator_reserve,
        worker_reserve_bytes=worker_reserve,
        estimated_memory_reservation_bytes=coordinator_reserve+workers*worker_reserve,
        memory_policy='Conservative scheduling reservation, not an OS memory limit; abort budget construction if it exceeds90% available RAM')


def collect_exact_keys(plan, case, cat, progress=None):
    """Trace the actual symbol request path using zero placeholder kernels.

    Stencil locations depend on fields and support, not placeholder values.
    Every numerical result from this trace is discarded. The real kernels
    subsequently recompute energy, forces, temperature and D36 unchanged.
    """
    import numpy as np
    from controls import spatial_field, spatial_profile, population
    from pysnspd.experimental.rectangular_spatial import RectangularSpatialFunctional
    if case.get('phase') != 'spatial' or case.get('finite_differences', True):
        raise ValueError('Parallel runner supports spatial cases with explicit finite_differences=false only')
    model = RectangularSpatialFunctional(cat, plan['geometry']['length_m'],
        plan['material']['width_m'], plan['material']['thickness_m'],
        elements_x=case['elements_x'], elements_y=case['elements_y'],
        degree=case['degree'], delta_regularizer_bar=case['delta_reg'])
    z = spatial_field(model, case['profile'])
    fields = model.sample_fields(z)
    p = population(cat, case['population'])
    keys = {}
    class Recorder:
        vacuum = cat.vacuum
        count_nodes, count_weights = cat.count_nodes, cat.count_weights
        def energy_kernel(self, amplitude, gamma):
            keys[(float(amplitude), float(gamma))] = None
            return (np.zeros_like(self.count_nodes),)*3
    model.catalog = Recorder()
    for i in range(model.cells):
        # Cover the exact base query even at q=0, plus all actual Gamma stencils.
        # Array abs (energy assembly) and scalar abs (temperature inversion)
        # may differ by an ULP; record both without rounding either key.
        model._potential(float(fields.amplitude_quadrature_bar[i]),
                         float(fields.gamma_quadrature_bar[i]), p, {})
        model._potential(float(abs(z[i])), float(fields.gamma_quadrature_bar[i]), p, {})
        model.principal_symbol(z[i], fields.derivative_quadrature_bar[i], p)
        if progress is not None:
            progress(i+1, model.cells, 'collect exact spectral requests')
    center = int(np.argmin(np.sum((model.dof_coordinates_bar-[model.length_bar/2, 0.])**2, axis=1)))
    exact = spatial_profile(model.dof_coordinates_bar, model.length_bar, case['profile'])[1]
    model.principal_symbol(z[center], exact[center], p)
    return list(keys)


def interleaved_jobs(case_requests):
    """Deduplicate exact keys while sharing one fair queue between cases."""
    queues = [deque((row['resolution'], key, row['id']) for key in row['keys']) for row in case_requests]
    jobs = {}
    while any(queues):
        for queue in queues:
            if not queue:
                continue
            resolution, key, owner = queue.popleft()
            item = jobs.setdefault((resolution, key), dict(resolution=resolution, key=key, owners=[]))
            if owner not in item['owners']:
                item['owners'].append(owner)
    return list(jobs.values())


def limit_thread_environment():
    """Call before importing NumPy/SciPy in the coordinator and spawned workers."""
    for name in THREAD_VARIABLES:
        os.environ[name] = '1'


def initialize_affinity(cpus, counter):
    """Reusable ProcessPool initializer: one CPU and one BLAS thread per worker."""
    limit_thread_environment()
    with counter.get_lock():
        index = counter.value
        counter.value += 1
    os.sched_setaffinity(0, {cpus[index % len(cpus)]})


def _init_worker(plan, cpus, counter):
    global _WORKER_PLAN, _WORKER_CATALOGUES
    initialize_affinity(cpus, counter)
    sys.path.insert(0, str(ROOT))
    _WORKER_PLAN, _WORKER_CATALOGUES = plan, {}


def _worker_kernel(resolution, key):
    from controls import catalogue
    if resolution not in _WORKER_CATALOGUES:
        _WORKER_CATALOGUES[resolution] = catalogue(_WORKER_PLAN, resolution)
    start = time.monotonic()
    values = _WORKER_CATALOGUES[resolution].source.energy_kernel(*key)
    return values, time.monotonic()-start, os.getpid()


def prefetch_jobs(jobs, executor, workers, on_result, kernel=_worker_kernel):
    """Keep bounded in-flight work; propagate the first failure without retry."""
    pending, iterator = {}, iter(jobs)
    def submit():
        for _ in range(max(0, 2*workers-len(pending))):
            job = next(iterator, None)
            if job is None:
                break
            future = executor.submit(kernel, job['resolution'], job['key'])
            pending[future] = job
    submit()
    try:
        while pending:
            completed, _ = wait(pending, return_when=FIRST_COMPLETED)
            for future in completed:
                job = pending.pop(future)
                values, elapsed, pid = future.result()
                on_result(job, values, elapsed, pid)
            submit()
    except BaseException:
        for future in pending:
            future.cancel()
        raise


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--plan', type=Path, required=True)
    parser.add_argument('--output-root', type=Path, required=True)
    parser.add_argument('--cpu-fraction', type=float, default=.9)
    parser.add_argument('--memory-fraction', type=float, default=.9)
    parser.add_argument('--max-workers', type=int)
    parser.add_argument('--execute', action='store_true')
    args = parser.parse_args()
    limit_thread_environment()
    sys.path.insert(0, str(ROOT))
    from controls import catalogue, spatial_case, material_reference
    from run_campaign import sha, sources, write
    plan = json.loads(args.plan.read_text(encoding='utf8'))
    if plan['schema'] != 'pysnspd.stage4A.controls.v1' or plan['photon'] is not False:
        raise ValueError('Unrecognized non-photon plan')
    cases = plan['cases']
    if not cases or len({case['id'] for case in cases}) != len(cases):
        raise ValueError('Nonempty unique case IDs required')
    if any(case.get('phase') != 'spatial' or case.get('finite_differences', True) for case in cases):
        raise ValueError('Only spatial cases with finite_differences=false are supported; use serial runner otherwise')
    resources = linux_resources()
    budget = resource_budget(resources, args.cpu_fraction, args.memory_fraction, args.max_workers)
    source_ids = sources()
    source_ids[Path(__file__).relative_to(ROOT).as_posix()] = sha(__file__)
    identity = dict(plan_sha256=sha(args.plan), sources=source_ids, task_ids=[c['id'] for c in cases],
        material=material_reference(plan), resource_snapshot=resources, budget=budget,
        execution='Shared exact-key pool across cases; request collection and final spatial assembly serial',
        thread_environment={name: os.environ[name] for name in THREAD_VARIABLES})
    if not args.execute:
        print(json.dumps(dict(status='DRY_RUN_NO_PHYSICS', **identity), indent=2))
        return
    if args.output_root.exists():
        raise FileExistsError('Use a NEW output root; previous results are preserved')
    args.output_root.mkdir(parents=True)
    write(args.output_root/'identity.json', identity)
    original_affinity = os.sched_getaffinity(0)
    os.sched_setaffinity(0, {budget['coordinator_cpu']})
    start, last_print = time.monotonic(), 0.
    log = (args.output_root/'progress.jsonl').open('x', encoding='utf8', buffering=1)
    records = []
    def event(name, **values):
        record = dict(event=name, elapsed_seconds=time.monotonic()-start, **values)
        text = json.dumps(record, ensure_ascii=False, allow_nan=False)
        log.write(text+'\n')
        print(text, flush=True)
    def progress_for(case_id):
        section_start = time.monotonic()
        previous_section = None
        def progress(done, total, section):
            nonlocal last_print, section_start, previous_section
            now = time.monotonic()
            if section != previous_section:
                section_start, previous_section = now, section
            if done != total and now-last_print < 5:
                return
            last_print = now
            event('PROGRESS', task=case_id, section=section, done=done, total=total,
                eta_seconds=(now-section_start)*(total-done)/done if done > 1 else None)
        return progress
    try:
        event('START', tasks=len(cases), budget=budget)
        catalogues, requests, prefetched = {}, [], {}
        for case in cases:
            cat = catalogue(plan, case['resolution'])
            catalogues[case['id']] = cat
            keys = collect_exact_keys(plan, case, cat, progress_for(case['id']))
            requests.append(dict(id=case['id'], resolution=case['resolution'], keys=keys))
            prefetched[case['id']] = {}
        jobs = interleaved_jobs(requests)
        totals = {row['id']: len(row['keys']) for row in requests}
        completed = {name: 0 for name in totals}
        worker_seconds, worker_pids = 0., set()
        done_jobs = 0
        prefetch_start = time.monotonic()
        event('PREFETCH_START', unique_queries=len(jobs), requests_by_case=totals)
        def received(job, values, elapsed, pid):
            nonlocal worker_seconds, done_jobs, last_print
            worker_seconds += elapsed
            worker_pids.add(pid)
            done_jobs += 1
            for name in job['owners']:
                prefetched[name][job['key']] = values
                completed[name] += 1
            now = time.monotonic()
            if done_jobs == len(jobs) or now-last_print >= 5:
                last_print = now
                fraction = done_jobs/len(jobs)
                bar = '#'*int(24*fraction)+'-'*(24-int(24*fraction))
                remaining = (now-prefetch_start)*(len(jobs)-done_jobs)/done_jobs
                print(f'[{bar}] espectros {done_jobs}/{len(jobs)} | ETA {remaining/60:.1f} min | casos {completed}', flush=True)
                event('SPECTRAL_PROGRESS', done=done_jobs, total=len(jobs), completed_by_case=completed,
                    total_by_case=totals, eta_seconds=remaining, worker_processes_seen=len(worker_pids))
        if jobs:
            context = mp.get_context('spawn')
            counter = context.Value('i', 0)
            with ProcessPoolExecutor(max_workers=budget['workers'], mp_context=context,
                initializer=_init_worker, initargs=(plan, budget['worker_affinity_cpus'], counter)) as executor:
                prefetch_jobs(jobs, executor, budget['workers'], received)
        event('PREFETCH_COMPLETE', runtime_seconds=time.monotonic()-prefetch_start,
            summed_worker_compute_seconds=worker_seconds, unique_queries=done_jobs,
            worker_processes_seen=len(worker_pids))
        for case in cases:
            cat = catalogues[case['id']]
            cat.install_prefetched(prefetched.pop(case['id']))
            output = args.output_root/case['id']
            output.mkdir()
            case_start = time.monotonic()
            event('ASSEMBLY_START', task=case['id'])
            result = spatial_case(plan, case, cat, progress_for(case['id']), output)
            result.update(case=case, runtime_seconds=time.monotonic()-case_start,
                runtime_scope='Serial assembly after shared parallel prefetch; see batch summary for full cost',
                catalogue_diagnostics=cat.diagnostics())
            write(output/'result.json', result)
            records.append(dict(id=case['id'], runtime_seconds=result['runtime_seconds'],
                result=case['id']+'/result.json', sha256=sha(output/'result.json'),
                files=[dict(path=p.relative_to(args.output_root).as_posix(), sha256=sha(p), bytes=p.stat().st_size)
                       for p in sorted(output.iterdir()) if p.is_file()]))
            event('TASK_COMPLETE', task=case['id'], diagnostics=cat.diagnostics())
        write(args.output_root/'summary.json', dict(status='COMPLETED_DIAGNOSTICS_NOT_PHYSICAL_ADMISSION',
            cases=records, runtime_seconds=time.monotonic()-start, budget=budget,
            unique_spectral_queries=done_jobs, summed_worker_compute_seconds=worker_seconds,
            stage4_complete=False, photon=False, time_trajectories=0, production_promotion=False))
        event('COMPLETE', tasks=len(records))
    except BaseException as exc:
        write(args.output_root/'failure.json', dict(status='STOPPED', reason=str(exc),
            exception=type(exc).__name__, completed_cases=records,
            policy='No retry, fallback, clipping or overwrite'))
        event('STOPPED', reason=str(exc))
        raise
    finally:
        log.close()
        os.sched_setaffinity(0, original_affinity)


if __name__ == '__main__':
    main()
