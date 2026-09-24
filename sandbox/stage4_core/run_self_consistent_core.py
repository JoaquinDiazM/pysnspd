"""Finite-sum thermal core minimization by exact alternating blocks.

The iteration index is an optimization index, never physical time. A checkpoint
contains d and its stationary spectral fields; next_d is a separate proposed
gap-block minimum. Every run has a fresh directory, including explicit resumes.
"""
from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import asdict
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
from pysnspd.experimental import thermal_spatial_usadel as thermal
from run_spatial_energy import build_case, reference_field, sha, write_json, REFERENCE

DEFAULT_PLAN = ROOT/'docs/implementation/stage4/spatial_energy_20260924/self_consistent_plan.json'
_GEOMETRY_CACHE = {}


def gap_coefficient(t, count):
    if not np.isfinite(t) or t <= 0 or type(count) is not int or count < 1:
        raise ValueError('Positive t and positive integer finite cutoff required')
    value = float(np.log(t)+np.sum(1/(np.arange(count)+.5)))
    if not np.isfinite(value) or value <= 0:
        raise ValueError('The declared finite-sum gap block is not strictly convex')
    return value


def exact_gap_block(graph, d, f_sum, t, count, fixed_nodes):
    """Exact free-node minimum and two independently evaluated energy drops."""
    coefficient = gap_coefficient(t, count)
    d = np.asarray(d, dtype=complex)
    f_sum = np.asarray(f_sum, dtype=complex)
    if d.shape != (graph.n_nodes,) or f_sum.shape != d.shape:
        raise ValueError('One gap and summed anomalous spectrum per node required')
    if not np.all(np.isfinite(d)) or not np.all(np.isfinite(f_sum)):
        raise ValueError('Nonfinite gap block inputs')
    target = 2*np.pi*t*f_sum/coefficient
    next_d = target.copy()
    next_d[fixed_nodes] = d[fixed_nodes]
    difference = d-next_d
    expected = -coefficient*float(np.dot(graph.area_weights, abs(difference)**2))
    # Evaluate the actual d-dependent quadratic terms at both endpoints. All
    # spectral-only terms cancel exactly because this block holds u fixed.
    density_old = coefficient*abs(d)**2-4*np.pi*t*np.real(np.conj(d)*f_sum)
    density_new = coefficient*abs(next_d)**2-4*np.pi*t*np.real(np.conj(next_d)*f_sum)
    observed = float(np.dot(graph.area_weights, density_new-density_old))
    scale = max(1., abs(expected), abs(observed), float(np.dot(graph.area_weights, abs(density_old))))
    if abs(observed-expected) > 128*np.finfo(float).eps*scale:
        raise RuntimeError('Exact gap-block energy identity failed')
    return next_d, dict(coefficient=coefficient, predicted_energy_change=expected,
        observed_energy_change=observed, identity_absolute_error=abs(observed-expected))


def residual_metrics(graph, radius, d, next_d, gap_reference, core_radius):
    free = np.ones(graph.n_nodes, bool)
    free[graph.boundary_nodes] = False
    core = free & (radius <= core_radius)
    if not np.any(core):
        raise ValueError('Core acceptance region contains no free nodes')
    residual = abs(d-next_d)/gap_reference
    return dict(core_mass_rms_relative=float(np.sqrt(np.dot(graph.area_weights[core], residual[core]**2)
        /np.sum(graph.area_weights[core]))),
        global_free_maximum_relative=float(np.max(residual[free])),
        core_maximum_relative=float(np.max(residual[core])))


def core_observables(graph, xy, radius, d, current, gap_reference, core_radius):
    core = radius <= core_radius
    suppressed = core & (abs(d) <= gap_reference/2)
    halfgap_radius = float(np.sqrt(np.sum(graph.area_weights[suppressed])/np.pi))
    minimum_index = int(np.argmin(abs(d)))
    xs, ys = np.unique(xy[:, 0]), np.unique(xy[:, 1])
    ids = np.arange(len(d)).reshape(len(xs), len(ys))
    lo_x, hi_x = int(np.argmin(abs(xs+2.))), int(np.argmin(abs(xs-2.)))
    lo_y, hi_y = int(np.argmin(abs(ys+2.))), int(np.argmin(abs(ys-2.)))
    loop = np.r_[ids[lo_x:hi_x+1, lo_y], ids[hi_x, lo_y+1:hi_y+1],
        ids[hi_x-1:lo_x-1:-1, hi_y], ids[lo_x, hi_y-1:lo_y:-1]]
    if len(loop) < 4 or np.min(abs(d[loop])) <= 1e-12*gap_reference:
        winding = None
    else:
        winding = float(np.sum(np.angle(np.conj(d[loop])*np.roll(d[loop], -1)))/(2*np.pi))
    edge_length = np.linalg.norm(xy[graph.edges[:, 1]]-xy[graph.edges[:, 0]], axis=1)
    dual_width = graph.conductance*edge_length
    active = dual_width > 0
    return dict(halfgap_area_equivalent_radius_ell0=halfgap_radius,
        halfgap_area_definition='Area within r<=core_radius with |d|<=gap_reference/2; sqrt(area/pi)',
        minimum_gap_relative=float(abs(d[minimum_index])/gap_reference),
        minimum_gap_coordinate_ell0=xy[minimum_index].tolist(),
        winding_on_square_near_two_ell0=winding,
        maximum_absolute_link_current_bar=float(np.max(abs(current))),
        maximum_absolute_link_current_per_dual_width=float(np.max(abs(current[active])/dual_width[active])))


def worker(job):
    start = time.monotonic()
    key = json.dumps(job['case'], sort_keys=True)
    if key not in _GEOMETRY_CACHE:
        _GEOMETRY_CACHE[key] = build_case(job['case'], job['gap'])
    graph, xy, radius, direction, initial_d, bump = _GEOMETRY_CACHE[key]
    reference_u = reference_field(job['n'], radius, direction)
    initial_u = reference_u if job['initial_u'] is None else job['initial_u']
    solution = thermal.solve_frequency(graph, job['d'], 2*np.pi*job['t']*(job['n']+.5),
        fixed_nodes=graph.boundary_nodes, fixed_u=reference_u[graph.boundary_nodes],
        initial_u=initial_u, tol=job['spectral_tolerance'], max_iterations=job['max_newton_iterations'])
    return dict(case_id=job['case']['id'], n=job['n'], solution=solution,
        seconds=time.monotonic()-start, pid=os.getpid())


def checkpoint_storage_budget(plan, output, pilot_sweeps=None, resumed=None):
    """Conservative uncompressed storage for this run, including its last sweep."""
    resumed = {} if resumed is None else resumed
    cadence = plan['checkpoint_cadence']
    cases = {}
    for case in plan['cases']:
        previous = resumed.get(case['id'])
        completed = 0 if previous is None else previous['completed_sweeps']
        accepted = False if previous is None else previous['metrics']['core_criterion_met']
        remaining = max(0, plan['max_outer_iterations']-completed)
        if pilot_sweeps is not None:
            remaining = min(remaining, pilot_sweeps)
        if accepted:
            remaining = 0
        nodes, modes = case['nx']*case['ny'], case['matsubara_count']
        edges = (case['nx']-1)*case['ny']+case['nx']*(case['ny']-1)
        # u/f complex, g real; spatial fields/topology; small ZIP/header margin.
        raw = 40*nodes*modes+72*nodes+32*edges+8*(2*case['nx']+2*case['ny']-4)
        per_checkpoint = int(np.ceil(1.01*raw))+65536
        checkpoints = min(remaining, int(np.ceil(remaining/cadence))+2)
        cases[case['id']] = dict(remaining_maximum_sweeps=remaining,
            projected_full_checkpoints=checkpoints, bytes_per_checkpoint=per_checkpoint,
            projected_bytes=checkpoints*per_checkpoint)
    parent = Path(output).resolve().parent
    while not parent.exists():
        parent = parent.parent
    free = shutil.disk_usage(parent).free
    projected = sum(case['projected_bytes'] for case in cases.values())
    budget = dict(checkpoint_cadence=cadence, maximum_recomputed_sweeps=cadence-1,
        policy='Full fields at first/cadence/accepted/pilot-end/budget-end; metrics every sweep; explicit resume may recompute at most four completed unsaved sweeps',
        projection='Conservative uncompressed fields plus header margin; actual compressed files can be smaller',
        available_disk_bytes=free, available_fraction_limit=.9,
        projected_bytes=projected, cases=cases)
    if projected > .9*free:
        raise RuntimeError(f'Checkpoint storage projection {projected} bytes exceeds 90% of {free} available bytes')
    return budget


def save_checkpoint(output, case_id, sweep, state, observation, next_d, metrics, *, full_fields=True):
    name = f'{case_id}_sweep{sweep:04d}'
    path = output/'checkpoints'/(name+'.npz')
    graph = state['geometry'][0]
    ordered = state['solutions']
    if full_fields:
        np.savez_compressed(path, d=state['d'], next_d=next_d,
            u=np.asarray([value.u for value in ordered]),
            f=np.asarray([value.f for value in ordered]), g=np.asarray([value.g for value in ordered]),
            force=observation.gap_gradient/graph.area_weights, current_bar=observation.current_bar,
            coordinates_bar=state['geometry'][1], area_weights=graph.area_weights,
            edges=graph.edges, conductance=graph.conductance, boundary_nodes=graph.boundary_nodes)
    record = dict(case_id=case_id, completed_sweeps=sweep,
        fields_path=str(path.resolve()) if full_fields else None,
        fields_sha256=sha(path) if full_fields else None, full_checkpoint=full_fields, metrics=metrics,
        consistent_pair='d,u,f,g,force,current are evaluated together before the separately stored next_d block',
        next_d_is_not_a_stationary_spectral_pair=True)
    write_json(output/'checkpoints'/(name+'.json'), record)
    return record


def load_resume(directory, identity):
    directory = Path(directory)
    old = json.loads((directory/'identity.json').read_text(encoding='utf8'))
    for key in ('plan_sha256', 'sources', 'radial_reference_identity_sha256'):
        if old[key] != identity[key]:
            raise ValueError('Resume requires identical plan, sources and reference: '+key)
    indices = sorted(directory.glob('resume_index_*.json'))
    if not indices:
        raise ValueError('No completed sweep checkpoint available for explicit resume')
    index = json.loads(indices[-1].read_text(encoding='utf8'))
    for item in index['cases'].values():
        if sha(item['fields_path']) != item['fields_sha256']:
            raise ValueError('Resume checkpoint hash changed')
    return index['cases']


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--plan', type=Path, default=DEFAULT_PLAN)
    parser.add_argument('--output-root', type=Path, required=True)
    parser.add_argument('--max-workers', type=int)
    parser.add_argument('--pilot-sweeps', type=int, choices=(1, 2))
    parser.add_argument('--resume-from', type=Path)
    parser.add_argument('--execute', action='store_true')
    args = parser.parse_args()
    plan = json.loads(args.plan.read_text(encoding='utf8'))
    if plan['schema'] != 'pysnspd.stage4.self_consistent_thermal_core.v1':
        raise ValueError('Unknown plan schema')
    if plan['checkpoint_cadence'] != 5:
        raise ValueError('This campaign uses the declared five-sweep checkpoint cadence')
    if len({case['id'] for case in plan['cases']}) != len(plan['cases']):
        raise ValueError('Duplicate case identifiers')
    gap = json.loads((REFERENCE/'identity.json').read_text())['gap_kBTc']
    t = plan['T_K']/plan['Tc_K']
    for case in plan['cases']:
        gap_coefficient(t, case['matsubara_count'])
    budget = resource_budget(linux_resources(), .9, .9, args.max_workers)
    sources = [Path(__file__), HERE/'run_spatial_energy.py', HERE/'parallel_runtime.py',
        ROOT/'pysnspd/experimental/thermal_spatial_usadel.py']
    identity = dict(schema=plan['schema'], plan_sha256=sha(args.plan),
        sources={path.relative_to(ROOT).as_posix():sha(path) for path in sources},
        radial_reference_identity_sha256=sha(REFERENCE/'identity.json'), gap_reference=gap,
        budget=budget, pilot_sweeps=args.pilot_sweeps,
        resumed_from=None if args.resume_from is None else str(args.resume_from.resolve()),
        scope='Alternating minimization of one finite-N thermal graph action; optimization index is not physical time',
        production_changed=False, physical_time_steps=0)
    resumed = {} if args.resume_from is None else load_resume(args.resume_from, identity)
    identity['checkpoint_budget'] = checkpoint_storage_budget(plan, args.output_root, args.pilot_sweeps, resumed)
    if not args.execute:
        print(json.dumps(dict(status='DRY_RUN', **identity), indent=2)); return
    if args.output_root.exists():
        raise FileExistsError('Choose a fresh output root, also when resuming')
    args.output_root.mkdir(parents=True)
    (args.output_root/'checkpoints').mkdir()
    write_json(args.output_root/'identity.json', identity)
    states = {}
    for case in plan['cases']:
        geometry = build_case(case, gap)
        state = dict(case=case, geometry=geometry, d=geometry[4].copy(), warm=None,
            completed=0, accepted=False, previous_postgap_energy=None, last=None, last_checkpoint=None)
        if case['id'] in resumed:
            saved = resumed[case['id']]
            with np.load(saved['fields_path']) as arrays:
                state['d'] = arrays['next_d'].copy()
                state['warm'] = arrays['u'].copy()
            state.update(completed=saved['completed_sweeps'], accepted=saved['metrics']['core_criterion_met'],
                previous_postgap_energy=saved['metrics']['energy_after_gap_block'], last=saved, last_checkpoint=saved)
        state['stop_sweep'] = min(plan['max_outer_iterations'], state['completed']+(args.pilot_sweeps or plan['max_outer_iterations']))
        states[case['id']] = state
    started = time.monotonic()
    original_affinity = os.sched_getaffinity(0)
    log = (args.output_root/'progress.jsonl').open('x', encoding='utf8', buffering=1)
    completed_jobs, last_print, index = 0, 0., 0
    def event(name, **fields):
        text = json.dumps(dict(event=name, elapsed_seconds=time.monotonic()-started, **fields), allow_nan=False)
        log.write(text+'\n'); print(text, flush=True)
    try:
        event('START', budget=budget, cases=len(states), max_outer_iterations=plan['max_outer_iterations'])
        os.sched_setaffinity(0, {budget['coordinator_cpu']})
        context = mp.get_context('spawn'); counter = context.Value('i', 0)
        with ProcessPoolExecutor(max_workers=budget['workers'], mp_context=context,
                initializer=initialize_affinity, initargs=(budget['worker_affinity_cpus'], counter)) as pool:
            while True:
                active = [state for state in states.values() if not state['accepted'] and state['completed'] < state['stop_sweep']]
                if not active:
                    break
                jobs = []
                for state in active:
                    state['pending'] = {}; state['mode_records'] = []
                    for n in range(state['case']['matsubara_count']):
                        jobs.append(dict(case=state['case'], d=state['d'], n=n, gap=gap, t=t,
                            initial_u=None if state['warm'] is None else state['warm'][n],
                            spectral_tolerance=plan['spectral_tolerance'], max_newton_iterations=plan['max_newton_iterations']))
                jobs.sort(key=lambda job:(job['n'], job['case']['id']))
                futures = {pool.submit(worker, job):job['case']['id'] for job in jobs}
                wave_done = 0
                try:
                    for future in as_completed(futures):
                        result = future.result(); state = states[result['case_id']]
                        solution = result.pop('solution'); state['pending'][result['n']] = solution
                        state['mode_records'].append(dict(**result, residual=solution.residual,
                            gradient_residual=solution.gradient_residual, newton_iterations=solution.iterations,
                            line_search_steps=list(solution.line_search_steps)))
                        completed_jobs += 1; wave_done += 1
                        now = time.monotonic()
                        if now-last_print > 3 or wave_done == len(jobs):
                            last_print = now
                            remaining = len(jobs)-wave_done+sum((s['stop_sweep']-s['completed']-1)*s['case']['matsubara_count'] for s in active)
                            eta = (now-started)/completed_jobs*remaining
                            fraction = wave_done/len(jobs)
                            bar = '#'*int(24*fraction)+'-'*(24-int(24*fraction))
                            print(f'[{bar}] espectros del barrido {wave_done}/{len(jobs)} | transcurrido {(now-started)/60:.1f} min | estimación para presupuesto restante {eta/60:.1f} min (coste medio observado; no predice convergencia)', flush=True)
                            event('PROGRESS', completed_jobs=completed_jobs, wave_complete=wave_done,
                                wave_total=len(jobs), remaining_maximum_jobs=remaining,
                                estimated_remaining_work_budget_seconds=eta)
                except BaseException:
                    for future in futures:
                        future.cancel()
                    raise
                for state in active:
                    graph, xy, radius, direction, initial_d, bump = state['geometry']
                    count = state['case']['matsubara_count']
                    state['solutions'] = [state['pending'][n] for n in range(count)]
                    observation = thermal.evaluate_thermal(graph, state['d'], t, state['solutions'])
                    f_sum = np.sum([solution.f for solution in state['solutions']], axis=0)
                    next_d, block = exact_gap_block(graph, state['d'], f_sum, t, count, graph.boundary_nodes)
                    residual = residual_metrics(graph, radius, state['d'], next_d, gap, plan['core_radius_ell0'])
                    accepted = residual['core_mass_rms_relative'] <= plan['core_mass_rms_relative_tolerance']
                    previous = state['previous_postgap_energy']
                    spectral_change = None if previous is None else observation.energy-previous
                    if previous is not None and spectral_change > 512*np.finfo(float).eps*max(1., abs(previous), abs(observation.energy)):
                        raise RuntimeError('The spectral block increased the common energy')
                    metrics = dict(**residual, energy=observation.energy,
                        energy_after_gap_block=observation.energy+block['observed_energy_change'], gap_block=block,
                        spectral_block_energy_change=spectral_change, core_criterion_met=accepted,
                        maximum_spectral_residual=observation.maximum_spectral_residual,
                        observables=core_observables(graph, xy, radius, state['d'], observation.current_bar, gap, plan['core_radius_ell0']),
                        modes=state['mode_records'], physical_time_steps=0, no_tail=True)
                    state['completed'] += 1
                    full_fields = (state['completed'] == 1 or state['completed'] % plan['checkpoint_cadence'] == 0
                        or accepted or state['completed'] == state['stop_sweep'])
                    state['last'] = save_checkpoint(args.output_root, state['case']['id'], state['completed'], state,
                        observation, next_d, metrics, full_fields=full_fields)
                    if full_fields:
                        state['last_checkpoint'] = state['last']
                    event('CASE_SWEEP', case=state['case']['id'], sweep=state['completed'],
                        core_mass_rms_relative=residual['core_mass_rms_relative'], energy=observation.energy,
                        core_criterion_met=accepted, gap_energy_change=block['observed_energy_change'])
                    state.update(d=next_d, warm=np.asarray([s.u for s in state['solutions']]), accepted=accepted,
                        previous_postgap_energy=metrics['energy_after_gap_block'])
                    del state['solutions'], state['pending'], state['mode_records']
                index += 1
                write_json(args.output_root/f'resume_index_{index:04d}.json',
                    dict(maximum_recomputed_sweeps=plan['checkpoint_cadence']-1,
                        cases={key:state['last_checkpoint'] for key,state in states.items() if state['last_checkpoint'] is not None}))
        all_accepted = all(state['accepted'] for state in states.values())
        status = ('PILOT_ONLY_NOT_PHYSICAL_ADMISSION' if args.pilot_sweeps is not None else
            'FINITE_SUM_CORE_STATIONARY' if all_accepted else 'INCOMPLETE_MAXIMUM_ITERATIONS')
        summary = dict(status=status, all_core_criteria_met=all_accepted, runtime_seconds=time.monotonic()-started,
            completed_jobs=completed_jobs, cases={key:state['last'] for key,state in states.items()},
            checkpoint_budget=identity['checkpoint_budget'],
            stage4_complete=False, production_changed=False, physical_time_steps=0,
            scope='Finite cutoff stationary core with fixed gap and spectral boundary; not stability, barrier, dynamics or device admission')
        write_json(args.output_root/'summary.json', summary)
        event('COMPLETE', status=status)
    except BaseException as exc:
        write_json(args.output_root/'failure.json', dict(reason=str(exc), exception=type(exc).__name__,
            completed_jobs=completed_jobs, policy='Stop without retry or overwrite; explicit resume only from a complete checkpoint index'))
        event('STOPPED', reason=str(exc)); raise
    finally:
        os.sched_setaffinity(0, original_affinity); log.close()


if __name__ == '__main__':
    main()
