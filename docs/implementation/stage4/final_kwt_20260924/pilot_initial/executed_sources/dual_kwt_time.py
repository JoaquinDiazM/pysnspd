"""Smooth dual-mesh thermal trajectory with the actual inherited KWT Euler step.

The uniform finite-Matsubara reference is analytic at each spectral frequency.
Every nonuniform probe retains full spatial spectral solves. No GL force,
linear-response replacement, ETD integrator, photon or population closure is used.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import multiprocessing as mp
import os
from pathlib import Path
import shutil
import sys
import time
import traceback

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(HERE))
from parallel_runtime import limit_thread_environment, linux_resources, resource_budget
limit_thread_environment()
sys.path.insert(0, str(ROOT))
import numpy as np
from scipy.optimize import brentq
from pysnspd.experimental import thermal_spatial_usadel as thermal
from pysnspd.experimental.thermal_stable_newton import solve_frequency
from pysnspd.experimental.thermal_snapshot import spectral_difference
from pysnspd.experimental.thermal_weak_response import ThermalKWTNormal
from pysnspd.experimental.heredado_kwt_bridge import inherited_kwt_step


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def save_json(path, value):
    with Path(path).open('x', encoding='utf8') as stream:
        json.dump(value, stream, indent=2, allow_nan=False)
        stream.write('\n')


def load_graph(path):
    with np.load(path) as data:
        graph = thermal.ThermalGraph(data['area_weights'], data['edges'],
            data['conductance'], data['coordinates_bar'], data['fixed_nodes'])
        profile = data['smooth_profile'].copy()
    return graph, profile


def uniform_reference(graph, temperature_ratio, count):
    """Exact uniform solution for this finite positive-frequency sum."""
    epsilon = 2*np.pi*temperature_ratio*(np.arange(count)+.5)
    def equation(d):
        return np.log(temperature_ratio)+2*np.pi*temperature_ratio*np.sum(
            1/epsilon-1/np.hypot(epsilon, d))
    amplitude = brentq(equation, 1e-8, 4., xtol=1e-14)
    gap = np.full(graph.n_nodes, amplitude, complex)
    gradient = 2*graph.area_weights*gap*equation(amplitude)
    return gap, gradient, epsilon


def explicit_step_bound(graph, d0, epsilon, plan):
    """Conservative uniform linearized Euler bound, including normal potential.

    lambda bounds eigenvalues of M^-1 L by Gershgorin. Hphase(lambda) is
    the EXACT uniform spectral Schur-complement symbol on this same graph.
    The normal-potential contribution uses its lambda->0 maximum. The small
    nonlinear probe is checked independently by dt versus dt/2 and energy.
    This is a starting bound, not a theorem for arbitrary nonuniform states.
    """
    model = ThermalKWTNormal(graph, d0, **material_options(plan))
    degree = np.zeros(graph.n_nodes)
    np.add.at(degree, graph.edges[:, 0], graph.conductance)
    np.add.at(degree, graph.edges[:, 1], graph.conductance)
    lam = float(np.max(2*degree/graph.area_weights))
    d = float(d0[0].real)
    ratio = plan['T_K']/plan['Tc_K']
    g = epsilon/np.hypot(epsilon, d)
    phase_h = float(4*np.pi*ratio*np.sum(g*g*lam/(epsilon*(epsilon+g*lam))))
    radial_h = float(2*np.log(ratio)+4*np.pi*ratio*np.sum(
        1/epsilon-g**3/(epsilon+g*lam)))
    factor = float(model.denominator[0]/graph.area_weights[0])
    angular = float(model.stretch[0]*phase_h/(factor*model.tD_ps))
    radial = float(radial_h/(factor*model.stretch[0]*model.tD_ps))
    normal = float(.5*d*d*4*np.pi*ratio*np.sum(1/(epsilon*epsilon+d*d))/model.tD_ps)
    rate = max(radial, angular+normal)
    return dict(gershgorin_lambda_max=lam, phase_material_rate_bound_per_ps=angular,
        radial_material_rate_bound_per_ps=radial, normal_rate_bound_per_ps=normal,
        maximum_rate_bound_per_ps=rate, primary_step_ps=float(plan['euler_safety']/rate),
        euler_safety=plan['euler_safety'], scope='Uniform linearization with fixed contacts; nonlinear acceptance uses temporal refinement and integrated dissipation')


def material_options(plan):
    return {key: plan[key] for key in ('T_K', 'Tc_K', 'tau_ee_Tc_ps', 'tau_ep_Tc_ps')}


def spectral_actor(pipe, cpu, graph, tasks, d0, epsilon, plan):
    """Resident guesses; independent (trajectory,frequency) pairs share CPUs."""
    try:
        os.sched_setaffinity(0, {cpu})
        cache = {}
        scale0 = graph.area_weights
        free = np.ones(graph.n_nodes, bool)
        free[graph.boundary_nodes] = False
        weight = 2*np.pi*plan['T_K']/plan['Tc_K']
        pipe.send(dict(event='READY', cpu=cpu, task_count=len(tasks)))
        while True:
            command, fields = pipe.recv()
            if command == 'STOP':
                break
            if command != 'EVALUATE':
                raise ValueError('Unknown spectral worker command')
            output = {key: [np.zeros(graph.n_nodes, complex),
                np.zeros(len(graph.edges)), 0.] for key in fields}
            maximum = 0.
            solved = reused = iterations = 0
            for name, n in tasks:
                if name not in fields:
                    continue
                gap = fields[name]
                e = epsilon[n]
                guess = cache.get((name, n))
                if guess is None:
                    guess = gap/e
                    guess[graph.boundary_nodes] = d0[graph.boundary_nodes]/e
                value = thermal.spectral_energy_gradient(graph, gap, e, guess)
                residual = float(np.max(abs(value.residual[free])/
                    (scale0[free]*np.maximum(1., abs(gap[free])))))
                if residual <= plan['spectral_tolerance']:
                    u = guess
                    reused += 1
                else:
                    solution = solve_frequency(graph, gap, e, initial_u=guess,
                        fixed_nodes=graph.boundary_nodes,
                        fixed_u=d0[graph.boundary_nodes]/e,
                        tol=plan['spectral_tolerance'],
                        max_iterations=plan['maximum_newton_iterations'])
                    u, residual = solution.u, solution.residual
                    iterations += solution.iterations
                    solved += 1
                cache[(name, n)] = u
                diff = spectral_difference(graph, d0, d0/e, gap, u, e)
                output[name][0] += weight*diff['gap_gradient_difference']
                output[name][1] += weight*diff['current_difference']
                output[name][2] += weight*diff['renormalized_energy_difference']
                maximum = max(maximum, residual)
            pipe.send(dict(output=output, maximum_residual=maximum,
                newton_solves=solved, stationary_reuses=reused, newton_iterations=iterations))
    except BaseException as error:
        try:
            pipe.send(dict(event='ERROR', reason=str(error), traceback=traceback.format_exc()))
        except (BrokenPipeError, EOFError):
            pass
    finally:
        pipe.close()


class SpectralPool:
    def __init__(self, graph, names, d0, epsilon, plan, budget, event):
        self.graph, self.d0 = graph, d0
        self.logt = np.log(plan['T_K']/plan['Tc_K'])
        self.processes, self.pipes = [], []
        self.batches = self.newton_solves = self.stationary_reuses = self.iterations = 0
        self.maximum_residual = 0.
        tasks = [(name, n) for n in range(len(epsilon)) for name in names]
        context = mp.get_context('spawn')
        try:
            for index, cpu in enumerate(budget['worker_affinity_cpus']):
                parent, child = context.Pipe()
                process = context.Process(target=spectral_actor, args=(child, cpu,
                    graph, tasks[index::budget['workers']], d0, epsilon, plan))
                process.start(); child.close()
                self.processes.append(process); self.pipes.append(parent)
            for index, pipe in enumerate(self.pipes):
                result = pipe.recv()
                if result.get('event') != 'READY':
                    raise RuntimeError(str(result))
                event('WORKER_READY', worker=index, **{k:v for k,v in result.items() if k != 'event'})
        except BaseException:
            self.close()
            raise

    def evaluate(self, fields):
        for pipe in self.pipes:
            pipe.send(('EVALUATE', fields))
        mass = self.graph.area_weights
        results = {}
        for name, gap in fields.items():
            delta = gap-self.d0
            results[name] = [2*mass*delta*self.logt, np.zeros(len(self.graph.edges)),
                float(np.dot(mass, 2*np.real(np.conj(self.d0)*delta)+abs(delta)**2)*self.logt)]
        for pipe in self.pipes:
            part = pipe.recv()
            if part.get('event') == 'ERROR':
                raise RuntimeError(str(part))
            for name, data in part['output'].items():
                for index in range(3):
                    results[name][index] += data[index]
            self.newton_solves += part['newton_solves']
            self.stationary_reuses += part['stationary_reuses']
            self.iterations += part['newton_iterations']
            self.maximum_residual = max(self.maximum_residual, part['maximum_residual'])
        self.batches += 1
        return results

    def close(self):
        for process, pipe in zip(self.processes, self.pipes):
            if process.is_alive():
                try:
                    pipe.send(('STOP', None))
                except (BrokenPipeError, EOFError):
                    pass
        for process in self.processes:
            process.join(timeout=3)
            if process.is_alive():
                process.terminate(); process.join(timeout=3)
        for pipe in self.pipes:
            pipe.close()


def response(graph, gap, data, G0, plan):
    gradient, current, energy = data
    model = ThermalKWTNormal(graph, gap, **material_options(plan))
    result = model.response(G0+gradient, current)
    result.update(energy=float(energy), gradient=G0+gradient, current=current,
        loss_per_ps=(result['kwt_loss']+result['normal_loss'])/model.tD_ps)
    return model, result


def observables(graph, gap, base, result):
    mass = graph.area_weights
    return dict(gap_difference=gap-base, current=result['current'],
        force_density=result['gradient']/mass,
        phase_torque_density=np.imag(np.conj(gap)*result['gradient'])/mass)


def observation_record(time_ps, gap, d0, result, integral, initial_energy):
    return dict(time_ps=float(time_ps), energy_excess=result['energy'],
        initial_energy_excess=float(initial_energy), integrated_loss=float(integral),
        integrated_balance_residual=float(result['energy']-initial_energy+integral),
        instantaneous_power_residual=result['dissipation_residual'],
        continuity_max=result['continuity_max'], noether_max=result['noether_max'],
        loss_per_ps=result['loss_per_ps'],
        maximum_relative_gap_change=float(np.max(abs(gap-d0))/abs(d0[0])),
        min_gap_kBTc=float(np.min(abs(gap))))


def norm(graph, kind, values):
    if kind == 'current':
        active = graph.conductance > 0
        return float(np.sqrt(np.sum(abs(values[active])**2/graph.conductance[active])))
    return float(np.sqrt(np.dot(graph.area_weights, abs(values)**2)))


def compare_observations(graph, observations, plan):
    records = []
    for index, (time_ps, fields) in enumerate(observations):
        for probe in ('amplitude', 'angular_phase'):
            for kind in ('gap_difference', 'current', 'force_density', 'phase_torque_density'):
                left, right = fields['primary_'+probe][kind], fields['refined_'+probe][kind]
                initial = observations[0][1]
                # The same observable's largest initial probe sets the floor;
                # an identically-zero cross-response never defines a relative tolerance.
                scale = max(norm(graph, kind, initial['refined_'+name][kind])
                            for name in ('amplitude', 'angular_phase'))
                error = norm(graph, kind, left-right)
                limit = plan['refinement_absolute_norm']+plan['refinement_relative_limit']*scale
                records.append(dict(time_ps=time_ps, probe=probe, observable=kind,
                    error_L2=error, initial_signal_scale=scale, admitted_limit=limit,
                    relative_to_initial_signal=None if scale == 0 else error/scale,
                    admitted=bool(error <= limit)))
    return dict(records=records, all_admitted=all(row['admitted'] for row in records),
        norm='Node: sqrt(sum m|x|^2); edge current: sqrt(sum |I|^2/c), inactive edges omitted',
        scale='Largest initial amplitude/phase response of the SAME observable; no normalization by a vanishing late tail')


def execute(plan, graph, profile, d0, G0, epsilon, bound, pool, output, event, horizon,
            *, maximum_steps=None):
    names = ['primary_amplitude', 'primary_angular_phase',
             'refined_amplitude', 'refined_angular_phase']
    amplitude = plan['perturbation_amplitude']
    initial = dict(amplitude=d0*(1+amplitude*profile),
                   angular_phase=d0*np.exp(1j*amplitude*profile))
    fields = {name: initial[name.split('_', 1)[1]].copy() for name in names}
    for gap in fields.values():
        gap[graph.boundary_nodes] = d0[graph.boundary_nodes]
    raw = pool.evaluate(fields)
    values = {name: response(graph, fields[name], raw[name], G0, plan) for name in names}
    initial_energies = {name: pair[1]['energy'] for name, pair in values.items()}
    if min(initial_energies.values()) <= 0:
        raise RuntimeError('Perturbations must have positive resolved excess energy')
    integrals = {name: 0. for name in names}
    maximum_balance = {name: 0. for name in names}
    maximum_energy_increase = {name: 0. for name in names}
    maximum_energy_above_initial = {name: 0. for name in names}
    history, observations = [], []
    targets = sorted(set([0., *[float(t) for t in plan['observation_times_ps'] if t <= horizon], float(horizon)]))
    times = {name: 0. for name in names}
    step_counts = {name: 0 for name in names}
    maximum_dt = bound['primary_step_ps']
    actual_steps = 0
    started = time.monotonic()
    def snapshot(target):
        index = len(observations)
        obs = {name: observables(graph, fields[name], d0, values[name][1]) for name in names}
        observations.append((float(target), obs))
        arrays = {'baseline_gap': d0, 'area_weights': graph.area_weights,
            'coordinates_bar': graph.coordinates_bar, 'edges': graph.edges,
            'conductance': graph.conductance, 'fixed_nodes': graph.boundary_nodes,
            'time_ps': np.array(target)}
        records = {}
        for name in names:
            arrays[name+'_gap'] = fields[name]
            for key, value in obs[name].items():
                arrays[name+'_'+key] = value
            arrays[name+'_potential_v'] = values[name][1]['potential_v']
            records[name] = observation_record(target, fields[name], d0,
                values[name][1], integrals[name], initial_energies[name])
        path = output/f'observation_{index:03d}.npz'
        np.savez_compressed(path, **arrays)
        save_json(output/f'observation_{index:03d}.json', dict(time_ps=target,
            states=records, fields_path=path.name, fields_sha256=sha(path),
            accepted_steps=step_counts.copy(), baseline='Exact uniform stationary finite-sum solution; zero currents and dissipation'))
        event('OBSERVATION', time_ps=target, observations=len(observations),
              states=records, fields_path=path.name)
    snapshot(0.)
    for target in targets[1:]:
        # Use an integer common partition: exact matching horizons, no 1e-16 ps final step.
        start = times[names[0]]
        segments = max(1, int(np.ceil((target-start)/maximum_dt)))
        primary_dt = (target-start)/segments
        for segment in range(segments):
            for substep in range(2):
                active = names if substep == 0 else names[2:]
                next_fields = {}
                for name in active:
                    model, result = values[name]
                    dt = primary_dt if name.startswith('primary') else primary_dt/2
                    update = inherited_kwt_step(model, result['gradient'], result['current'],
                        dt_ps=dt, delta0_over_kBTc=plan['delta0_over_kBTc'])
                    if np.any(~np.isfinite(update.gap)):
                        raise RuntimeError('Nonfinite inherited KWT field')
                    if np.max(abs(update.gap-d0))/abs(d0[0]) > plan['maximum_relative_displacement']:
                        raise RuntimeError('Large departure from the admitted smooth thermal neighborhood')
                    if not np.array_equal(update.gap[graph.boundary_nodes], d0[graph.boundary_nodes]):
                        raise RuntimeError('Fixed contacts changed')
                    next_fields[name] = update.gap
                evaluated = pool.evaluate(next_fields)
                for name in active:
                    dt = primary_dt if name.startswith('primary') else primary_dt/2
                    previous = values[name][1]
                    fields[name] = next_fields[name]
                    values[name] = response(graph, fields[name], evaluated[name], G0, plan)
                    now = values[name][1]
                    integrals[name] += .5*dt*(previous['loss_per_ps']+now['loss_per_ps'])
                    balance_residual = now['energy']-initial_energies[name]+integrals[name]
                    maximum_balance[name] = max(maximum_balance[name], abs(balance_residual))
                    maximum_energy_increase[name] = max(maximum_energy_increase[name],
                        now['energy']-previous['energy'])
                    maximum_energy_above_initial[name] = max(maximum_energy_above_initial[name],
                        now['energy']-initial_energies[name])
                    times[name] = min(target, times[name]+dt)
                    step_counts[name] += 1
                    history.append(dict(trajectory=name, time_ps=times[name], dt_ps=dt,
                        energy_excess=now['energy'], integrated_loss=integrals[name],
                        balance_residual=balance_residual))
                actual_steps += 1
                elapsed = time.monotonic()-started
                covered = min(times.values())
                fraction = covered/horizon
                event('KWT_STEP', time_ps=covered, horizon_ps=horizon,
                    actual_steps=step_counts.copy(), spectral_batches=pool.batches,
                    newton_solves=pool.newton_solves, stationary_reuses=pool.stationary_reuses,
                    eta_seconds=None if not covered else elapsed*(horizon-covered)/covered)
                print(f"[{('#'*int(24*fraction)).ljust(24,'-')}] dual KWT {covered:.6g}/{horizon:g} ps | "
                    f"{elapsed:.1f}s | ETA {(elapsed*(horizon-covered)/covered if covered else 0):.1f}s", flush=True)
                if maximum_steps is not None and actual_steps >= maximum_steps:
                    path = output/'pilot_last_states.npz'
                    np.savez_compressed(path, **fields, times_ps=np.array([times[name] for name in names]))
                    save_json(output/'step_history.json', history)
                    return dict(status='BOUNDED_STEP_PILOT_COMPLETE', stage4_complete=False,
                        covered_ps=covered, pilot_spectral_batches=pool.batches,
                        accepted_steps=step_counts, runtime_seconds=elapsed,
                        projected_full_seconds=None if not covered else elapsed*horizon/covered,
                        initial_energy_excess=initial_energies, step_bound=bound,
                        maximum_balance_residual=maximum_balance,
                        maximum_single_step_energy_increase=maximum_energy_increase,
                        reason='Pilot stopped after registered count; not a completed horizon or accuracy certificate')
            # Checkpoints at a common accepted primary/refined time.
            arrays = dict(fields, time_ps=np.array(start+(segment+1)*primary_dt),
                integrals=np.array([integrals[name] for name in names]),
                initial_energy_excess=np.array([initial_energies[name] for name in names]),
                accepted_steps=np.array([step_counts[name] for name in names]))
            path = output/f'checkpoint_{sum(step_counts.values()):07d}.npz'
            np.savez_compressed(path, **arrays)
        times = {name: target for name in names}
        snapshot(target)
    comparison = compare_observations(graph, observations, plan)
    balances = {}
    for name in names:
        residual = values[name][1]['energy']-initial_energies[name]+integrals[name]
        limit = plan['energy_absolute_floor']+plan['energy_relative_limit']*initial_energies[name]
        balances[name] = dict(residual=float(residual), initial_energy_excess=initial_energies[name],
            relative_to_initial_excess=float(residual/initial_energies[name]), admitted_limit=float(limit),
            maximum_absolute_balance_residual=maximum_balance[name],
            maximum_single_step_energy_increase=maximum_energy_increase[name],
            maximum_energy_above_initial=maximum_energy_above_initial[name],
            admitted=bool(max(maximum_balance[name], maximum_energy_increase[name],
                maximum_energy_above_initial[name]) <= limit))
    save_json(output/'step_history.json', history)
    save_json(output/'refinement.json', comparison)
    passed = comparison['all_admitted'] and all(row['admitted'] for row in balances.values())
    return dict(status='DUAL_KWT_THERMAL_COMPLETE' if passed else 'COMPLETED_WITH_DIAGNOSTIC_MARGINS',
        completed_horizon_ps=horizon, accepted_steps=step_counts, integrated_balance=balances,
        all_practical_criteria_met=passed, refinement_all_admitted=comparison['all_admitted'],
        stage4_complete=False, nonthermal_work_closed=False, photon=False,
        production_changed=False, time_order=1, energy_contract='Fixed-temperature free-energy decrease plus KWT and normal Joule dissipation; not coupled-population internal energy')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('plan', 'mesh', 'output-root'):
        parser.add_argument('--'+name, type=Path, required=True)
    parser.add_argument('--execute', action='store_true')
    parser.add_argument('--pilot-steps', type=int)
    args = parser.parse_args()
    plan = json.loads(args.plan.read_text())
    if args.pilot_steps is not None and args.pilot_steps < 1:
        raise ValueError('Positive pilot step count required')
    graph, profile = load_graph(args.mesh)
    if sha(args.mesh) != plan['mesh_sha256']:
        raise ValueError('Mesh hash does not match plan')
    for path, digest in plan.get('source_sha256', {}).items():
        if sha(ROOT/path) != digest:
            raise ValueError('Registered source changed: '+path)
    d0, G0, epsilon = uniform_reference(graph, plan['T_K']/plan['Tc_K'], plan['matsubara_count'])
    bound = explicit_step_bound(graph, d0, epsilon, plan)
    resources = linux_resources()
    budget = resource_budget(resources, .9, .9, max_workers=plan.get('maximum_workers', 27))
    identity = dict(plan_sha256=sha(args.plan), mesh_sha256=sha(args.mesh), resources=resources,
        budget=budget, step_bound=bound, uniform_gap_kBTc=float(d0[0].real),
        maximum_uniform_force=float(np.max(abs(G0))), baseline_spatial_newton_solves=0,
        time_order=1, method='Actual TDGLSolver.solve_for_psi_squared through inherited_kwt_step',
        boundary_conditions='Fixed equilibrium gap and zero normal potential on end contacts only; natural insulating side walls; alpha=0',
        pilot_steps=args.pilot_steps, units='d=Delta/(kB Tc), epsilon=omega/(kB Tc), X=x/ell0, physical time ps')
    if not args.execute:
        print(json.dumps(dict(status='DRY_RUN', **identity), indent=2)); return
    if args.output_root.exists():
        raise FileExistsError('Preserve previous evidence; choose a fresh output root')
    parent = args.output_root.resolve().parent
    while not parent.exists():
        parent = parent.parent
    if plan.get('output_reserve_bytes', 2*1024**3) > .9*shutil.disk_usage(parent).free:
        raise RuntimeError('Insufficient disk space for checkpoint reserve')
    args.output_root.mkdir(parents=True)
    save_json(args.output_root/'identity.json', identity)
    shutil.copyfile(args.plan, args.output_root/'executed_plan.json')
    log = (args.output_root/'progress.jsonl').open('x', buffering=1)
    started = time.monotonic()
    def event(name, **values):
        line = json.dumps(dict(event=name, elapsed_seconds=time.monotonic()-started, **values))
        print(line, flush=True); log.write(line+'\n')
    affinity = os.sched_getaffinity(0)
    pool = None
    try:
        os.sched_setaffinity(0, {budget['coordinator_cpu']})
        event('START', **identity)
        names = ['primary_amplitude', 'primary_angular_phase', 'refined_amplitude', 'refined_angular_phase']
        pool = SpectralPool(graph, names, d0, epsilon, plan, budget, event)
        result = execute(plan, graph, profile, d0, G0, epsilon, bound, pool,
            args.output_root, event, float(plan['observation_times_ps'][-1]), maximum_steps=args.pilot_steps)
        result.update(runtime_seconds=time.monotonic()-started, newton_solves=pool.newton_solves,
            stationary_reuses=pool.stationary_reuses, maximum_spectral_residual=pool.maximum_residual,
            worker_count=budget['workers'], spectral_batches=pool.batches)
        save_json(args.output_root/'summary.json', result)
        event('COMPLETE', **result)
    except BaseException as error:
        save_json(args.output_root/'failure.json', dict(reason=str(error), type=type(error).__name__,
            completed=False, policy='Keep snapshots; no clipping, hidden fallback or physical change'))
        event('FAILED', reason=str(error)); raise
    finally:
        if pool is not None:
            pool.close()
        os.sched_setaffinity(0, affinity)
        log.close()


if __name__ == '__main__':
    main()
