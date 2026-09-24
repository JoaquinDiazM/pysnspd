"""Stationary current-carrying thermal reference on the admitted dual mesh.

Alternating exact spectral/gap blocks minimize one finite Matsubara action.
Iteration is an optimization index, never physical time. End contacts have a
declared phase difference; lateral edges retain natural boundary conditions.
The reference is for later weak dynamic coupling, not a photon trajectory.
"""
from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor, FIRST_COMPLETED, wait
import json
import multiprocessing as mp
import os
from pathlib import Path
import sys
import time

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from sandbox.stage4_core.parallel_runtime import (
    limit_thread_environment, linux_resources, initialize_affinity, THREAD_VARIABLES,
)
limit_thread_environment()
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import dijkstra
from pysnspd.experimental import thermal_spatial_usadel as thermal
from pysnspd.experimental.thermal_stable_newton import solve_frequency
from sandbox.stage4_core.dual_kwt_time import load_graph, uniform_reference
from sandbox.stage4_core.run_self_consistent_core import exact_gap_block
from sandbox.stage4_core.coupled_campaign import atomic_json, sha
from sandbox.stage4_core.coupled_response import admit_runtime, runtime_telemetry

SCHEMA = 'pysnspd.stage4.biased_strip_reference.v1'
_GRAPH = None
_PLAN = None


def validate_plan(plan):
    if plan.get('schema') != SCHEMA:
        raise ValueError('Unsupported stationary-reference schema')
    for key in ('matsubara_count', 'max_outer_iterations', 'max_newton_iterations', 'checkpoint_cadence'):
        if type(plan.get(key)) is not int or plan[key] < 1:
            raise ValueError(key+' must be a positive integer')
    for key in ('gap_rms_relative_tolerance', 'spectral_tolerance'):
        if not np.isfinite(plan[key]) or plan[key] <= 0:
            raise ValueError(key+' must be positive and finite')
    if not 0 < plan['T_K'] < plan['Tc_K'] or not np.isfinite(plan['phase_bias']):
        raise ValueError('A finite phase bias and 0 < T < Tc are required')
    return plan


def contact_partition(graph):
    x = graph.coordinates_bar[:, 0]
    tolerance = 1e-10*max(1., np.ptp(x))
    left = np.flatnonzero(abs(x-x.min()) <= tolerance)
    right = np.flatnonzero(abs(x-x.max()) <= tolerance)
    fixed = np.sort(np.r_[left, right])
    if not np.array_equal(np.sort(graph.boundary_nodes), fixed):
        raise ValueError('Declared fixed_nodes must be the two end contacts only; sides are natural')
    if len(fixed) == graph.n_nodes:
        raise ValueError('At least one free interior node is required')
    return left, right, fixed


def initial_state(graph, plan):
    left, right, fixed = contact_partition(graph)
    t = plan['T_K']/plan['Tc_K']
    uniform, _, epsilon = uniform_reference(graph, t, plan['matsubara_count'])
    x = graph.coordinates_bar[:, 0]
    theta = plan['phase_bias']*((x-x.min())/np.ptp(x)-.5)
    d = uniform*np.exp(1j*theta)
    return d, float(uniform[0].real), epsilon


def phase_path(graph):
    """A shortest geometric interior path used only to report the phase branch."""
    left,right,_=contact_partition(graph)
    y=graph.coordinates_bar[:,1]
    start=int(left[np.argmin(abs(y[left]-np.mean(y[left])))])
    end=int(right[np.argmin(abs(y[right]-np.mean(y[right])))])
    active=graph.conductance>0
    tail,head=graph.edges[active].T
    length=np.linalg.norm(graph.coordinates_bar[head]-graph.coordinates_bar[tail],axis=1)
    adjacency=coo_matrix((np.r_[length,length],(np.r_[tail,head],np.r_[head,tail])),
        shape=(graph.n_nodes,graph.n_nodes)).tocsr()
    distances,predecessors=dijkstra(adjacency,indices=start,return_predecessors=True)
    if not np.isfinite(distances[end]):
        raise ValueError('End contacts are not connected through positive-conductance links')
    path=[end]
    while path[-1]!=start:
        path.append(int(predecessors[path[-1]]))
    return np.asarray(path[::-1],dtype=int)


def solve_mode(graph, plan, d, epsilon, initial_u=None):
    return solve_frequency(graph, d, epsilon, fixed_nodes=graph.boundary_nodes,
        fixed_u=d[graph.boundary_nodes]/epsilon, initial_u=initial_u,
        tol=plan['spectral_tolerance'], max_iterations=plan['max_newton_iterations'])


def initialize_worker(graph, plan, cpus, counter):
    global _GRAPH, _PLAN
    initialize_affinity(cpus, counter)
    _GRAPH, _PLAN = graph, plan


def worker(job):
    index, d, epsilon, warm = job
    started = time.monotonic()
    solution = solve_mode(_GRAPH, _PLAN, d, epsilon, warm)
    return index, solution, dict(elapsed_seconds=time.monotonic()-started,
        newton_iterations=solution.iterations, residual=solution.residual,
        backtracks=list(solution.line_search_steps), **runtime_telemetry())


def measure(graph, d, next_d, d0, t, solutions, observation, block):
    left, right, fixed = contact_partition(graph)
    free = np.ones(graph.n_nodes, bool); free[fixed] = False
    m = graph.area_weights
    residual = abs(d-next_d)/d0
    divergence = np.zeros(graph.n_nodes)
    np.add.at(divergence, graph.edges[:, 0], observation.current_bar)
    np.add.at(divergence, graph.edges[:, 1], -observation.current_bar)
    spectral_reaction = np.zeros(graph.n_nodes)
    for solution in solutions:
        evaluation = thermal.spectral_energy_gradient(graph, d, solution.epsilon, solution.u)
        spectral_reaction += 2*np.pi*t*np.imag(np.conj(solution.u)*evaluation.gradient_u)
    gap_reaction = np.imag(np.conj(d)*observation.gap_gradient)
    contacts = {}
    for name, ids in (('left',left), ('right',right)):
        contacts[name] = dict(outgoing_link_current_bar=float(np.sum(divergence[ids])),
            gap_phase_reaction_bar=float(np.sum(gap_reaction[ids])),
            spectral_phase_reaction_bar=float(np.sum(spectral_reaction[ids])),
            total_phase_reaction_bar=float(np.sum(gap_reaction[ids]+spectral_reaction[ids])))
    path=phase_path(graph)
    jumps=np.angle(np.conj(d[path[:-1]])*d[path[1:]])
    advance=float(np.sum(jumps))
    principal=float(np.angle(np.conj(d[path[0]])*d[path[-1]]))
    metrics = dict(gap_mass_rms_relative=float(np.sqrt(np.dot(m[free],residual[free]**2)/m[free].sum())),
        gap_max_relative=float(np.max(residual[free])),
        gap_residual_definition='|delta_bar-next_delta_bar|/uniform_gap_bar; mass-weighted RMS over free nodes',
        force_density_mass_rms=float(np.sqrt(np.sum(abs(observation.gap_gradient[free])**2/m[free])/m[free].sum())),
        energy=observation.energy, energy_after_exact_gap_block=observation.energy+block['observed_energy_change'],
        exact_gap_block=block, maximum_spectral_residual=observation.maximum_spectral_residual,
        minimum_gap_relative=float(np.min(abs(d))/d0),
        minimum_matsubara_g=float(min(np.min(solution.g) for solution in solutions)),
        maximum_free_current_divergence=float(np.max(abs(divergence[free]))),
        maximum_noether_residual=float(np.max(abs(observation.noether_residual))),
        maximum_absolute_link_current_bar=float(np.max(abs(observation.current_bar))),
        contact_reactions=contacts,phase_advance_rad=advance,
        phase_branch_integer=int(round((advance-principal)/(2*np.pi))),
        maximum_phase_jump_on_path_rad=float(np.max(abs(jumps))),
        phase_path_minimum_gap_relative=float(np.min(abs(d[path]))/d0),
        phase_branch_definition='Sum of principal link phase differences along the saved left-to-right shortest geometric path, alpha=0; diagnostic, not an imposed constraint')
    return metrics, divergence, gap_reaction, spectral_reaction


def checkpoint(output, sweep, graph, d, next_d, epsilon, solutions, observation, metrics,
               divergence, gap_reaction, spectral_reaction, *, final=False):
    name = f'state_{sweep:04d}.npz'
    target = output/'reference.npz' if final else output/'checkpoints'/name
    temporary = target.with_suffix('.tmp')
    with temporary.open('wb') as stream:
        np.savez_compressed(stream, delta_bar=d, next_delta_bar=next_d, epsilon_bar=epsilon,
            u=np.asarray([s.u for s in solutions]), f=np.asarray([s.f for s in solutions]),
            g=np.asarray([s.g for s in solutions]), force=observation.gap_gradient/graph.area_weights,
            current_bar=observation.current_bar, current_divergence_bar=divergence,
            gap_phase_reaction_bar=gap_reaction, spectral_phase_reaction_bar=spectral_reaction,
            coordinates_bar=graph.coordinates_bar, area_weights=graph.area_weights,
            edges=graph.edges, conductance=graph.conductance, boundary_nodes=graph.boundary_nodes,
            alpha=np.zeros(len(graph.edges)), phase_path_nodes=phase_path(graph))
        stream.flush(); os.fsync(stream.fileno())
    os.replace(temporary,target)
    record = dict(sweep=sweep, path=target.relative_to(output).as_posix(), sha256=sha(target), metrics=metrics,
        pair='delta_bar,u,f,g,force,current are consistent; next_delta_bar is a separate proposed update')
    atomic_json(output/'latest_checkpoint.json',record)
    return record


def solve_reference(graph, plan, evaluate, *, output=None, event=lambda *a,**kw:None):
    """The evaluator returns consecutive stationary positive-frequency spectra."""
    d, d0, epsilon = initial_state(graph, plan)
    t = plan['T_K']/plan['Tc_K']
    warm, previous_postgap = None, None
    history, last_checkpoint = [], None
    for sweep in range(1,plan['max_outer_iterations']+1):
        solutions, mode_records = evaluate(d,epsilon,warm,sweep)
        observation = thermal.evaluate_thermal(graph,d,t,solutions)
        f_sum = np.sum([solution.f for solution in solutions],axis=0)
        next_d, block = exact_gap_block(graph,d,f_sum,t,len(epsilon),graph.boundary_nodes)
        metrics, divergence, gap_reaction, spectral_reaction = measure(
            graph,d,next_d,d0,t,solutions,observation,block)
        metrics['spectral_block_energy_change'] = None if previous_postgap is None else observation.energy-previous_postgap
        # The useful stopping scale is the gap residual. This guard rejects
        # material energy growth, not subtraction noise in extensive energies.
        energy_floor = 1e-10*max(1.,abs(observation.energy),abs(previous_postgap or 0.))
        if previous_postgap is not None and observation.energy-previous_postgap > energy_floor:
            raise RuntimeError('Spectral minimization increased the common action beyond the roundoff guard')
        accepted = metrics['gap_mass_rms_relative'] <= plan['gap_rms_relative_tolerance']
        metrics.update(criterion_met=accepted, sweep=sweep, modes=mode_records, physical_time_steps=0)
        history.append(metrics)
        full = sweep == 1 or sweep % plan['checkpoint_cadence'] == 0 or accepted or sweep == plan['max_outer_iterations']
        if output is not None:
            atomic_json(output/f'sweep_{sweep:04d}.json',metrics)
            if full:
                last_checkpoint = checkpoint(output,sweep,graph,d,next_d,epsilon,solutions,observation,
                    metrics,divergence,gap_reaction,spectral_reaction,
                    final=accepted or sweep == plan['max_outer_iterations'])
        event('SWEEP_COMPLETE',sweep=sweep,maximum_sweeps=plan['max_outer_iterations'],
            gap_mass_rms_relative=metrics['gap_mass_rms_relative'],criterion_met=accepted,
            minimum_gap_relative=metrics['minimum_gap_relative'],
            contact_current_bar=metrics['contact_reactions']['left']['outgoing_link_current_bar'])
        if accepted:
            break
        d, warm = next_d, np.asarray([solution.u for solution in solutions])
        previous_postgap = metrics['energy_after_exact_gap_block']
    return dict(status='FINITE_SUM_STATIONARY' if accepted else 'INCOMPLETE_MAXIMUM_ITERATIONS',
        criterion_met=accepted, metrics=history[-1], completed_sweeps=sweep,
        gap_reference_bar=d0, latest_checkpoint=last_checkpoint,
        scope='Fixed phase-biased contacts, natural lateral boundaries, finite positive Matsubara action; iteration is not time',
        production_changed=False, physical_time_steps=0, stage4_complete=False)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--plan',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--workers',type=int,default=1)
    parser.add_argument('--execute',action='store_true')
    args=parser.parse_args()
    plan=validate_plan(json.loads(args.plan.read_text(encoding='utf8')))
    if args.output.exists():
        raise FileExistsError('A fresh output directory is required; no overwrite or implicit resume')
    mesh=ROOT/plan['mesh']
    if sha(mesh)!=plan['mesh_sha256']:
        raise ValueError('Mesh hash differs from the admitted plan')
    resources=linux_resources()
    inherited=os.environ.get('PYSNSPD_SHARED_ALLOCATION')
    budget=admit_runtime(args.workers,resources,json.loads(inherited) if inherited else None)
    graph,_=load_graph(mesh)
    contact_partition(graph)
    sources=[Path(__file__), ROOT/'sandbox/stage4_core/parallel_runtime.py',
        ROOT/'sandbox/stage4_core/coupled_response.py', ROOT/'sandbox/stage4_core/coupled_campaign.py',
        ROOT/'sandbox/stage4_core/dual_kwt_time.py', ROOT/'sandbox/stage4_core/run_self_consistent_core.py',
        ROOT/'pysnspd/experimental/thermal_spatial_usadel.py', ROOT/'pysnspd/experimental/thermal_stable_newton.py',
        ROOT/'pysnspd/experimental/thermal_snapshot.py']
    manifest=dict(schema=SCHEMA,status='DRY_RUN',plan=plan,plan_sha256=sha(args.plan),
        sources={p.relative_to(ROOT).as_posix():sha(p) for p in sources}, resources=resources,
        budget=budget,nodes=graph.n_nodes,edges=len(graph.edges),
        thread_environment={k:os.environ[k] for k in THREAD_VARIABLES})
    if not args.execute:
        print(json.dumps(manifest,indent=2));return 0
    args.output.mkdir(parents=True);(args.output/'checkpoints').mkdir()
    started=time.monotonic();manifest['status']='RUNNING'
    atomic_json(args.output/'manifest.json',manifest)
    original=os.sched_getaffinity(0)
    log=(args.output/'progress.jsonl').open('x',encoding='utf8',buffering=1)
    def event(name,**values):
        value=dict(event=name,elapsed_seconds=time.monotonic()-started,**values)
        line=json.dumps(value,allow_nan=False);log.write(line+'\n');print(line,flush=True)
    try:
        os.sched_setaffinity(0,{budget['coordinator_cpu']})
        event('START',budget=budget,phase_bias=plan['phase_bias'],nodes=graph.n_nodes)
        context=mp.get_context('spawn');counter=context.Value('i',0)
        with ProcessPoolExecutor(max_workers=args.workers,mp_context=context,
                initializer=initialize_worker,initargs=(graph,plan,budget['worker_affinity_cpus'],counter)) as pool:
            def evaluate(d,epsilon,warm,sweep):
                jobs=[(n,d,float(e),None if warm is None else warm[n]) for n,e in enumerate(epsilon)]
                futures={pool.submit(worker,job):job[0] for job in jobs}
                solutions=[None]*len(jobs); records=[None]*len(jobs); done=0;wave=time.monotonic()
                while futures:
                    ready,_=wait(futures,timeout=5.,return_when=FIRST_COMPLETED)
                    for future in ready:
                        futures.pop(future)
                        index,solution,record=future.result()
                        solutions[index],records[index]=solution,record;done+=1
                    elapsed=time.monotonic()-wave
                    remaining=(len(jobs)-done)+(plan['max_outer_iterations']-sweep)*len(jobs)
                    eta=elapsed*remaining/done if done else None
                    event('SPECTRAL_PROGRESS',sweep=sweep,completed=done,total=len(jobs),
                        fraction=((sweep-1)+done/len(jobs))/plan['max_outer_iterations'],eta_seconds=eta,
                        eta_basis='Remaining maximum iteration budget using this sweep mean; convergence can finish earlier')
                return solutions,records
            summary=solve_reference(graph,plan,evaluate,output=args.output,event=event)
        summary.update(elapsed_seconds=time.monotonic()-started,coordinator=runtime_telemetry())
        atomic_json(args.output/'summary.json',summary)
        manifest.update(status=summary['status'],elapsed_seconds=summary['elapsed_seconds'],
            summary_sha256=sha(args.output/'summary.json'),latest_checkpoint=summary['latest_checkpoint'])
        atomic_json(args.output/'manifest.json',manifest)
        event('COMPLETE',status=summary['status'],fraction=1.,eta_seconds=0.)
        return 0 if summary['criterion_met'] else 2
    except BaseException as error:
        manifest.update(status='FAILED',exception=type(error).__name__,reason=str(error),
            elapsed_seconds=time.monotonic()-started)
        atomic_json(args.output/'manifest.json',manifest);event('FAILED',reason=str(error));raise
    finally:
        os.sched_setaffinity(0,original);log.close()


if __name__=='__main__':
    raise SystemExit(main())
