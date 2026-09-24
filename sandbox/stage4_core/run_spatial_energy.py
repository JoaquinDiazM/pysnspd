"""Shared-budget 2D thermal Usadel reference; no physical time evolution.

The spectral degrees of freedom are solved independently at each positive
Matsubara frequency. All cases share one pool. Finite sums are deliberately
reported as finite sums; no force-only tail or local K0 correction is added.
"""
from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import asdict
import hashlib
import json
import multiprocessing as mp
import os
from pathlib import Path
import sys
import time

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(HERE))
from parallel_runtime import initialize_affinity, limit_thread_environment, linux_resources, resource_budget
limit_thread_environment()
sys.path.insert(0, str(ROOT))
import numpy as np
from scipy.interpolate import CubicHermiteSpline
from pysnspd.experimental import thermal_spatial_usadel as thermal

DATA = ROOT/'docs/implementation/stage4/spatial_energy_20260924'
REFERENCE = ROOT/'docs/implementation/stage4/followup_20260923/raw/stage4_radial_reference_20260923'
_CASE_CACHE = {}


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write_json(path, value):
    with Path(path).open('x', encoding='utf8') as stream:
        json.dump(value, stream, indent=2, ensure_ascii=False, allow_nan=False)
        stream.write('\n')


def build_case(case, gap):
    graph = thermal.rectangular_graph(case['nx'], case['ny'], case['length_ell0'], case['width_ell0'])
    coordinates = np.asarray(graph.coordinates_bar).copy()
    coordinates -= (coordinates.min(axis=0)+coordinates.max(axis=0))/2
    x, y = coordinates.T
    radius = np.hypot(x, y)
    direction = np.divide(x+1j*y, radius, out=np.zeros_like(x, complex), where=radius>0)
    d = gap*np.tanh(radius)*direction
    bump = np.maximum(0., 1-radius**2/16)**3
    if case['profile'] == 'asymmetric_vortex':
        d *= 1+x/4*bump
    elif case['profile'] != 'radial_vortex':
        raise ValueError('Unsupported profile')
    perturbation = case.get('perturbation')
    h = float(case.get('perturbation_step', 0.))
    if perturbation == 'amplitude':
        d += h*direction*np.tanh(radius)*bump
    elif perturbation == 'phase':
        d *= np.exp(1j*h*y/4*bump)
    elif perturbation is not None:
        raise ValueError('Unsupported perturbation')
    return graph, coordinates, radius, direction, d, bump


def reference_field(n, radius, direction):
    path = REFERENCE/'tasks'/f'R12_n{n:04d}.npz'
    with np.load(path) as arrays:
        r, theta, derivative = arrays['radius'], arrays['theta'], arrays['derivative']
        interpolant = CubicHermiteSpline(np.r_[0.,r], np.r_[0.,theta], np.r_[derivative[0],derivative])
        if max(radius) > r[-1]:
            raise ValueError('Matched radial reference does not cover the entire square')
        u = np.tan(interpolant(radius))*direction
    return u


def worker(job):
    started = time.monotonic()
    case = job['case']
    cache_key = json.dumps(case, sort_keys=True)
    if cache_key not in _CASE_CACHE:
        _CASE_CACHE[cache_key] = build_case(case, job['gap'])
    graph, xy, radius, direction, d, bump = _CASE_CACHE[cache_key]
    initial = reference_field(job['n'], radius, direction)
    eps = 2*np.pi*job['t']*(job['n']+.5)
    solution = thermal.solve_frequency(graph, d, eps,
        fixed_nodes=graph.boundary_nodes, fixed_u=initial[graph.boundary_nodes],
        initial_u=initial, tol=job['tolerance'], max_iterations=job['max_iterations'])
    metadata = {key:value for key,value in asdict(solution).items() if not isinstance(value,np.ndarray)}
    return dict(case_id=case['id'], n=job['n'], solution=solution, metadata=metadata,
                worker_seconds=time.monotonic()-started, worker_pid=os.getpid())


def summarize_case(plan, case, gap, stored, output):
    graph, xy, radius, direction, d, bump = build_case(case,gap)
    order = sorted(stored)
    if order != list(range(case['matsubara_count'])):
        raise ValueError('A physical sum needs every declared frequency, not the pilot subset')
    results = [stored[n] for n in order]
    cutoffs = case.get('report_cutoffs',[case['matsubara_count']])
    summaries = {}
    for count in cutoffs:
        observation = thermal.evaluate_thermal(graph,d,plan['T_K']/plan['Tc_K'],results[:count])
        gradient = observation.gap_gradient
        force = gradient/graph.area_weights
        amplitude_force = np.real(np.conj(direction)*force)
        core = (radius <= 4.) & (radius > 0.)
        norm = lambda v: float(np.sqrt(np.sum(graph.area_weights[core]*np.abs(v[core])**2)))
        ref_force = np.zeros_like(d)
        ref_f = []
        for n in order[:count]:
            u = reference_field(n,radius,direction)
            ref_f.append(u/np.sqrt(1+np.abs(u)**2))
        frequencies=2*np.pi*(plan['T_K']/plan['Tc_K'])*(np.arange(count)+.5)
        radial_gap=gap*np.tanh(radius)*direction
        ref_force = 2*radial_gap*np.log(plan['T_K']/plan['Tc_K'])+4*np.pi*(plan['T_K']/plan['Tc_K'])*np.sum(radial_gap[None,:]/frequencies[:,None]-np.asarray(ref_f),axis=0)
        phase_source = np.imag(np.conj(d)*gradient)
        divergence = np.zeros(len(d))
        np.add.at(divergence,graph.edges[:,0],observation.current_bar)
        np.add.at(divergence,graph.edges[:,1],-observation.current_bar)
        interior = np.ones(len(d),bool);interior[graph.boundary_nodes]=False
        ward = phase_source+divergence
        boundary = ~interior
        spectral_phase_mismatch = np.angle(results[0].f[core]*np.conj(d[core]))
        filename = f"{case['id']}_N{count}.npz"
        np.savez_compressed(output/filename,coordinates_bar=xy,area_weights=graph.area_weights,
            edges=graph.edges,conductance=graph.conductance,d=d,force=force,amplitude_force=amplitude_force,
            current_bar=observation.current_bar,reference_radial_force=ref_force,ward_residual=ward,
            f_lowest=results[0].f,g_lowest=results[0].g,core_mask=core,boundary_mask=boundary)
        ref_norm = norm(ref_force)
        summaries[str(count)] = dict(energy=observation.energy,core_force_L2=norm(force),
            radial_force_relative_difference=(norm(force-ref_force)/ref_norm if case['profile']=='radial_vortex' else None),
            reference_comparison_scope='Same prescribed radial vortex, same finite sum, interpolated independent BVP, matched spectral boundary; only radial profile has an oracle',
            maximum_interior_gauge_identity_residual=float(np.max(np.abs(ward[interior]))),
            maximum_lowest_mode_phase_mismatch_rad=float(np.max(np.abs(spectral_phase_mismatch))),
            amplitude_directional_work=float(np.real(np.vdot(gradient,direction*np.tanh(radius)*bump))),
            phase_directional_work=float(np.real(np.vdot(gradient,1j*d*xy[:,1]/4*bump))),
            fields_file=filename,fields_sha256=sha(output/filename),
            derivatives_hold_spectral_boundary_fixed=True,tail_correction_applied=False)
    return dict(case=case,counts=summaries)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--plan',type=Path,default=DATA/'campaign_plan.json')
    parser.add_argument('--output-root',type=Path,required=True)
    parser.add_argument('--max-workers',type=int)
    parser.add_argument('--pilot',action='store_true')
    parser.add_argument('--execute',action='store_true')
    args=parser.parse_args()
    plan=json.loads(args.plan.read_text(encoding='utf8'))
    if plan['schema']!='pysnspd.stage4.thermal_spatial_campaign.v1':raise ValueError('Unknown plan')
    base_identity=json.loads((REFERENCE/'identity.json').read_text())
    gap=base_identity['gap_kBTc']
    jobs=[]
    for case in plan['cases']:
        ns=(plan['pilot_modes'] if args.pilot and case['id'] in plan['pilot_cases']
            else [] if args.pilot else range(case['matsubara_count']))
        for n in ns:
            jobs.append(dict(case=case,n=n,gap=gap,t=plan['T_K']/plan['Tc_K'],
                tolerance=plan['spectral_tolerance'],max_iterations=plan['max_newton_iterations']))
    # Interleave equal-frequency work across geometries, never one process group per case.
    jobs.sort(key=lambda item:(item['n'],item['case']['id']))
    budget=resource_budget(linux_resources(),.9,.9,args.max_workers)
    sources=[Path(__file__),ROOT/'pysnspd/experimental/thermal_spatial_usadel.py',HERE/'parallel_runtime.py']
    identity=dict(plan_sha256=sha(args.plan),sources={p.relative_to(ROOT).as_posix():sha(p) for p in sources},
        radial_reference_identity_sha256=sha(REFERENCE/'identity.json'),gap_kBTc=gap,
        jobs=len(jobs),budget=budget,pilot=args.pilot,scope='Thermal stationary spectral states on prescribed2D fields; not a kinetic trajectory',
        physical_time_steps=0,production_changed=False)
    if not args.execute:
        print(json.dumps(dict(status='DRY_RUN',**identity),indent=2));return
    if args.output_root.exists():raise FileExistsError('Choose a new output root; no overwriting')
    args.output_root.mkdir(parents=True);(args.output_root/'modes').mkdir()
    write_json(args.output_root/'identity.json',identity)
    started=time.monotonic();records=[];stored={case['id']:{} for case in plan['cases']}
    affinity=os.sched_getaffinity(0)
    log=(args.output_root/'progress.jsonl').open('x',encoding='utf8',buffering=1)
    last_print=0.
    def event(name,**fields):
        line=json.dumps(dict(event=name,elapsed_seconds=time.monotonic()-started,**fields),allow_nan=False)
        log.write(line+'\n');print(line,flush=True)
    try:
        event('START',jobs=len(jobs),budget=budget)
        os.sched_setaffinity(0,{budget['coordinator_cpu']})
        context=mp.get_context('spawn');counter=context.Value('i',0)
        with ProcessPoolExecutor(max_workers=budget['workers'],mp_context=context,
            initializer=initialize_affinity,initargs=(budget['worker_affinity_cpus'],counter)) as pool:
            futures={pool.submit(worker,job):job for job in jobs}
            try:
                for future in as_completed(futures):
                    result=future.result();cid,n=result['case_id'],result['n'];solution=result.pop('solution')
                    stored[cid][n]=solution
                    path=args.output_root/'modes'/f'{cid}_n{n:04d}.npz'
                    np.savez_compressed(path,u=solution.u,f=solution.f,g=solution.g)
                    result.update(path=path.relative_to(args.output_root).as_posix(),sha256=sha(path))
                    records.append(result)
                    now=time.monotonic()
                    if now-last_print>3 or len(records)==len(jobs):
                        last_print=now;fraction=len(records)/len(jobs)
                        eta=(now-started)*(1-fraction)/fraction
                        bar='#'*int(24*fraction)+'-'*(24-int(24*fraction))
                        print(f'[{bar}] espectros2D {len(records)}/{len(jobs)} | {(now-started)/60:.1f} min | ETA {eta/60:.1f} min',flush=True)
                        event('PROGRESS',complete=len(records),eta_seconds=eta)
            except BaseException:
                for future in futures:future.cancel()
                raise
        summary=dict(status='PILOT_ONLY_NO_PHYSICAL_SUM' if args.pilot else 'THERMAL_2D_REFERENCE_COMPUTED',
            records=records,stage4_complete=False,production_changed=False,physical_time_steps=0)
        if not args.pilot:
            summary['cases']={case['id']:summarize_case(plan,case,gap,stored[case['id']],args.output_root) for case in plan['cases']}
            comparisons={}
            for pair in plan.get('directional_checks',[]):
                n=str(pair['count']);base=summary['cases'][pair['base']]['counts'][n]
                high=summary['cases'][pair['plus']]['counts'][n];low=summary['cases'][pair['minus']]['counts'][n]
                finite=(high['energy']-low['energy'])/(2*pair['step'])
                predicted=base[pair['kind']+'_directional_work']
                comparisons[pair['kind']]=dict(finite_energy_derivative=finite,force_work=predicted,
                    absolute_difference=abs(finite-predicted),relative_difference=abs(finite-predicted)/max(abs(predicted),1e-15),
                    step=pair['step'],no_universal_acceptance_tolerance=True)
            summary['directional_checks']=comparisons
        summary['runtime_seconds']=time.monotonic()-started
        write_json(args.output_root/'summary.json',summary);event('COMPLETE')
    except BaseException as exc:
        write_json(args.output_root/'failure.json',dict(reason=str(exc),exception=type(exc).__name__,completed=records,
            policy='No retry, clipping, different closure, overwrite or silent continuation'))
        event('STOPPED',reason=str(exc));raise
    finally:
        os.sched_setaffinity(0,affinity);log.close()


if __name__=='__main__':main()
