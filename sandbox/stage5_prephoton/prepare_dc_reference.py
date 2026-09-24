"""Relax the dual mesh against intrinsic current-carrying bulk reservoirs.

The iteration counter is not time. No normal contacts, AC source, background
force subtraction or manufactured volume forcing is used.
"""
from __future__ import annotations
import argparse
from concurrent.futures import as_completed
import json
import multiprocessing as mp
import os
from pathlib import Path
import sys
import time

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from sandbox.stage4_core.parallel_runtime import limit_thread_environment, linux_resources, initialize_affinity
limit_thread_environment()
import numpy as np
from pysnspd.experimental import thermal_spatial_usadel as thermal
from pysnspd.experimental.bulk_current_reference import (
    solve_bulk_current, solve_bulk_reference, intrinsic_boundary_fields,
    differential_inductance_per_length_H_m, physical_scales,
)
from pysnspd.experimental.thermal_stable_newton import solve_frequency
from sandbox.stage4_core.dual_kwt_time import load_graph
from sandbox.stage4_core.biased_strip_reference import contact_partition, measure, checkpoint
from sandbox.stage4_core.run_self_consistent_core import exact_gap_block
from sandbox.stage4_core.coupled_campaign import atomic_json, sha
from sandbox.stage4_core.coupled_response import admit_runtime
from sandbox.stage4_core.biased_coupled_response import fail_fast_pool

SCHEMA = 'pysnspd.stage5.intrinsic_dc_reference.v1'
_G = _P = _BC = None


def initialize(graph, plan, boundary, cpus, counter):
    global _G, _P, _BC
    initialize_affinity(cpus, counter)
    _G, _P, _BC = graph, plan, boundary


def mode(job):
    n, gap, ep, initial = job
    return n, solve_frequency(_G, gap, ep, initial_u=initial,
        fixed_nodes=_G.boundary_nodes, fixed_u=_BC[n],
        tol=_P['spectral_tolerance'], max_iterations=_P['max_newton_iterations'])


def material(plan):
    return {k:plan[k] for k in ('width_m','Tc_K','diffusion_m2_s','sheet_resistance_ohm')}


def bulk_state(graph, plan):
    scale=physical_scales(**{k:plan[k] for k in ('Tc_K','diffusion_m2_s','sheet_resistance_ohm')})
    if graph.coordinates_bar is None:
        raise ValueError('The intrinsic strip requires physical graph coordinates')
    mesh_width=float(np.ptp(graph.coordinates_bar[:,1])*scale['ell0_m'])
    if not np.isclose(mesh_width,plan['width_m'],rtol=1e-8,atol=0.):
        raise ValueError('Physical mesh width differs from the width used to calibrate the bulk current')
    if not np.isfinite(plan['current_A']) or plan['current_A']==0:
        raise ValueError('This current-biased admission requires a finite nonzero DC current')
    reference = solve_bulk_current(plan['T_K']/plan['Tc_K'],plan['matsubara_count'],
        plan['current_A'],**material(plan))
    fields = intrinsic_boundary_fields(graph, reference)
    return reference, fields


def fixed_inductance_partition(graph, reference, plan):
    """Freeze the omitted series wire at the declared DC operating point.

    The active strip inductance partitions additively with length. Increasing
    the resolved window removes exactly that window's equilibrium inductance
    from the external circuit; it does not change the physical device.
    """
    scale=physical_scales(**{k:plan[k] for k in ('Tc_K','diffusion_m2_s','sheet_resistance_ohm')})
    length=float(np.ptp(graph.coordinates_bar[:,0])*scale['ell0_m'])
    active=float(plan['active_length_m']);added=float(plan['added_inductance_H'])
    if not np.isfinite(active) or not np.isfinite(added) or not 0<length<=active or added<0:
        raise ValueError('Invalid partition of fixed external inductance')
    per_length=differential_inductance_per_length_H_m(reference,**material(plan))
    external=added+(active-length)*per_length
    if external<=0:
        raise ValueError('A positive fixed external inductance is required')
    return dict(resolved_length_m=length,
        bulk_differential_inductance_per_length_H_m=per_length,
        fixed_external_inductance_H=external,
        resolved_equilibrium_inductance_H=length*per_length,
        total_equilibrium_inductance_H=added+active*per_length)


def metrics_bulk(graph, gap, observation, reference, plan, current_unit):
    left,right,_=contact_partition(graph)
    free=np.ones(graph.n_nodes,bool);free[graph.boundary_nodes]=False
    mass=graph.area_weights
    modulus=abs(gap)
    mean=float(np.average(modulus[free],weights=mass[free]))
    divergence=np.zeros(graph.n_nodes)
    np.add.at(divergence,graph.edges[:,0],observation.current_bar)
    np.add.at(divergence,graph.edges[:,1],-observation.current_bar)
    il=float(divergence[left].sum()*current_unit)
    ir=float(-divergence[right].sum()*current_unit)
    achieved=.5*(il+ir)
    x=graph.coordinates_bar[:,0]
    values=np.unique(x); cuts=.5*(values[:-1]+values[1:])
    tail,head=graph.edges.T
    cross=np.asarray([np.sum(observation.current_bar*((x[tail]<cut).astype(int)
            -(x[head]<cut).astype(int)))*current_unit for cut in cuts])
    return dict(gap_bulk_bar=reference.gap_bar,
        gap_relative_rms_vs_bulk=float(np.sqrt(np.average((modulus[free]/reference.gap_bar-1)**2,weights=mass[free]))),
        gap_peak_to_peak_relative=float(np.ptp(modulus[free])/reference.gap_bar),
        gap_mean_relative_to_bulk=mean/reference.gap_bar,
        current_left_A=il,current_right_A=ir,current_reference_A=achieved,
        current_target_A=plan['current_A'],current_target_relative_error=abs(achieved/plan['current_A']-1),
        crosscut_max_relative_error=float(np.max(abs(cross-achieved))/abs(plan['current_A'])),
        crosscut_count=len(cuts),
        definition='Gap magnitudes on every free node, area-weighted RMS. Current is oriented sum over every distinct-x interval; crosscut deviations from mean end current are normalized by the nonzero nominal current.')


def run(graph, plan, pool, output, event):
    ref,initial=bulk_state(graph,plan)
    ep=ref.epsilon_bar;gap=initial['delta_bar'];warm=initial['u']
    t=plan['T_K']/plan['Tc_K']
    zero=solve_bulk_reference(t,plan['matsubara_count'],0.)
    scale=physical_scales(**{k:plan[k] for k in ('Tc_K','diffusion_m2_s','sheet_resistance_ohm')})
    previous_post=None; history=[]; accepted=False
    for sweep in range(1,plan['max_outer_iterations']+1):
        start=time.monotonic();solutions=[None]*len(ep)
        jobs=[pool.submit(mode,(n,gap,float(e),warm[n])) for n,e in enumerate(ep)]
        for done,future in enumerate(as_completed(jobs),1):
            n,solution=future.result();solutions[n]=solution
            if done==len(ep) or done%max(1,len(ep)//4)==0:
                event('SPECTRA',sweep=sweep,completed=done,total=len(ep),
                    fraction=((sweep-1)+done/len(ep))/plan['max_outer_iterations'],
                    eta_seconds=(time.monotonic()-start)*(len(ep)-done+(plan['max_outer_iterations']-sweep)*len(ep))/done,
                    eta_basis='Maximum remaining optimization budget; convergence may stop earlier')
        obs=thermal.evaluate_thermal(graph,gap,t,solutions)
        next_gap,block=exact_gap_block(graph,gap,np.sum([s.f for s in solutions],axis=0),t,len(ep),graph.boundary_nodes)
        measures,div,gap_react,spectral_react=measure(graph,gap,next_gap,zero.gap_bar,t,solutions,obs,block)
        bulk=metrics_bulk(graph,gap,obs,ref,plan,scale['current_unit_A'])
        accepted=measures['gap_mass_rms_relative']<=plan['gap_rms_relative_tolerance']
        row=dict(sweep=sweep,gap_residual=measures['gap_mass_rms_relative'],criterion_met=accepted,**bulk)
        history.append(row)
        if previous_post is not None and obs.energy>previous_post+1e-10*max(1.,abs(previous_post)):
            raise RuntimeError('Common action grew in a spectral minimization block')
        event('DC_REFERENCE',**row,fraction=sweep/plan['max_outer_iterations'])
        if accepted or sweep==plan['max_outer_iterations']:
            saved=checkpoint(output,sweep,graph,gap,next_gap,ep,solutions,obs,measures,div,gap_react,spectral_react,final=True)
            break
        gap=next_gap;warm=np.asarray([s.u for s in solutions]);previous_post=obs.energy+block['observed_energy_change']
    inductance=fixed_inductance_partition(graph,ref,plan)
    physical_admission=accepted and all((
        bulk['gap_relative_rms_vs_bulk']<=plan['bulk_gap_relative_limit'],
        bulk['gap_peak_to_peak_relative']<=plan['bulk_gap_relative_limit'],
        bulk['current_target_relative_error']<=plan['current_relative_limit'],
        bulk['crosscut_max_relative_error']<=plan['current_relative_limit']))
    result=dict(schema=SCHEMA,status='DC_REFERENCE_ADMITTED' if physical_admission else 'DC_REFERENCE_REQUIRES_REVIEW',
        stationary=accepted,admitted_for_dark_hold=physical_admission,plan=plan,history=history,
        metrics=measures,bulk=bulk,reference_sha256=saved['sha256'],q_bar=ref.q_bar,
        **inductance,
        scope='Intrinsic DC bulk continuation on both cuts; finite Matsubara equilibrium. Not a photon or a general nonthermal energy certificate.',
        equilibrium_populations='hL=tanh(E/(2 kB Tb)), hT=0; phonons nB(Omega,Tb). No population perturbation.',
        circuit_initialization='Ib=Is=measured mean end current, vc=0, Vb=Rb*Is; record nominal-current discrepancy; Lext frozen once.',
        physical_time_steps=0,production_changed=False)
    atomic_json(output/'summary.json',result)
    return result


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--plan',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--workers',type=int,required=True)
    args=parser.parse_args()
    plan=json.loads(args.plan.read_text())
    if plan.get('schema')!=SCHEMA:raise ValueError('Unsupported DC reference plan')
    mesh=ROOT/plan['mesh']
    if sha(mesh)!=plan['mesh_sha256']:raise ValueError('Mesh changed')
    graph,_=load_graph(mesh);contact_partition(graph)
    inherited=os.environ.get('PYSNSPD_SHARED_ALLOCATION')
    budget=admit_runtime(args.workers,linux_resources(),json.loads(inherited) if inherited else None)
    if args.output.exists():raise FileExistsError('Output exists; no overwrite or implicit resume')
    ref,fields=bulk_state(graph,plan)
    args.output.mkdir(parents=True)
    manifest=dict(schema=SCHEMA,status='RUNNING',plan=plan,budget=budget,
        boundary='Current-carrying intrinsic reservoir: gap AND spectral fields from same bulk state',
        q_bar=ref.q_bar,bulk_gap_bar=ref.gap_bar)
    atomic_json(args.output/'manifest.json',manifest)
    started=time.monotonic();original=os.sched_getaffinity(0)
    with (args.output/'progress.jsonl').open('x',encoding='utf8',buffering=1) as log:
        def event(kind,**values):
            line=json.dumps(dict(event=kind,elapsed_seconds=time.monotonic()-started,**values),allow_nan=False)
            print(line,flush=True);log.write(line+'\n')
        try:
            os.sched_setaffinity(0,{budget['coordinator_cpu']})
            context=mp.get_context('spawn');counter=context.Value('i',0)
            with fail_fast_pool(args.output,manifest,max_workers=args.workers,mp_context=context,initializer=initialize,
                    initargs=(graph,plan,fields['u'][:,graph.boundary_nodes],budget['worker_affinity_cpus'],counter)) as pool:
                result=run(graph,plan,pool,args.output,event)
            manifest.update(status=result['status'],elapsed_seconds=time.monotonic()-started,
                reference_sha256=result['reference_sha256'])
            atomic_json(args.output/'manifest.json',manifest)
            event('COMPLETE',status=result['status'],fraction=1.,eta_seconds=0.)
            return 0 if result['admitted_for_dark_hold'] else 2
        except BaseException as error:
            manifest.update(status='FAILED',reason=str(error),exception=type(error).__name__)
            atomic_json(args.output/'manifest.json',manifest)
            raise
        finally:os.sched_setaffinity(0,original)


if __name__=='__main__':raise SystemExit(main())
