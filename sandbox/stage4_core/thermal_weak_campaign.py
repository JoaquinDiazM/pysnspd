"""Prepare and validate a full-node weak thermal evolution operator.

Independent Matsubara tasks share one capped process pool. Each task solves
the accepted fixed gap once, factors its spectral Jacobian once, applies
multiple full-node directions, and compares against independent nonlinear
perturbed spectral roots. No physical time integration is claimed.
"""
from __future__ import annotations
import argparse
from concurrent.futures import ProcessPoolExecutor,as_completed
import hashlib
import json
import multiprocessing as mp
import os
from pathlib import Path
import shutil
import sys
import time

HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]
sys.path.insert(0,str(HERE))
from parallel_runtime import initialize_affinity,limit_thread_environment,linux_resources,resource_budget
limit_thread_environment()
sys.path.insert(0,str(ROOT))
import numpy as np
from pysnspd.experimental import thermal_spatial_usadel as thermal
from pysnspd.experimental.thermal_weak_response import SpectralTangent,ThermalKWTNormal,divergence
from run_spatial_energy import reference_field

DEFAULT_PLAN=ROOT/'docs/implementation/stage4/moment_review_20260924/thermal_weak/plan.json'
_CACHE={}


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def read(path):return json.loads(Path(path).read_text(encoding='utf8'))
def write(path,value):
    with Path(path).open('x',encoding='utf8') as stream:
        json.dump(value,stream,indent=2,allow_nan=False);stream.write('\n')


def fields(plan):
    key=plan['case']['fields_path']
    if key not in _CACHE:
        with np.load(ROOT/key) as a:
            d,xy=a['d'],a['coordinates_bar']
            graph=thermal.ThermalGraph(a['area_weights'],a['edges'],a['conductance'],xy,a['boundary_nodes'])
        radius=np.linalg.norm(xy,axis=1)
        winding=np.divide(xy[:,0]+1j*xy[:,1],radius,out=np.zeros(len(d),complex),where=radius>0)
        bump=np.maximum(0.,1-radius**2/plan['probe_radius_ell0']**2)**3
        directions=np.asarray([d*bump,1j*d*xy[:,0]/plan['probe_radius_ell0']*bump])
        directions[:,graph.boundary_nodes]=0
        _CACHE[key]=(graph,d,xy,radius,winding,directions)
    return _CACHE[key]


def worker(job):
    started=time.monotonic();plan,n=job['plan'],job['n']
    graph,d,xy,radius,winding,directions=fields(plan)
    t=plan['T_K']/plan['Tc_K'];epsilon=2*np.pi*t*(n+.5)
    initial=reference_field(n,radius,winding)
    fixed=graph.boundary_nodes
    def solved(gap,guess):
        return thermal.solve_frequency(graph,gap,epsilon,fixed_nodes=fixed,fixed_u=initial[fixed],
            initial_u=guess,tol=plan['spectral_tolerance'],max_iterations=plan['maximum_newton_iterations'])
    base=solved(d,initial)
    factor=SpectralTangent(graph,d,base)
    tangent=[factor.apply(v) for v in directions]
    value=thermal.spectral_energy_gradient(graph,d,epsilon,base.u)
    m=graph.area_weights;weight=2*np.pi*t
    shape=(len(directions),len(plan['difference_steps']),2)
    force_states=np.zeros(shape+(len(d),),complex)
    current_states=np.zeros(shape+(len(graph.edges),))
    energies=np.zeros(shape)
    residuals=[];iterations=[]
    for j,(direction,action) in enumerate(zip(directions,tangent)):
        for q,h in enumerate(plan['difference_steps']):
            for p,sign in enumerate((-1,1)):
                gap=d+sign*h*direction
                other=solved(gap,base.u+sign*h*action.du)
                evaluated=thermal.spectral_energy_gradient(graph,gap,epsilon,other.u)
                force_states[j,q,p]=2*weight*m*(gap/epsilon-evaluated.f)
                current_states[j,q,p]=-weight*evaluated.link_derivative_alpha
                energies[j,q,p]=weight*(evaluated.energy+np.dot(m,abs(gap)**2)/epsilon)
                residuals.append(other.residual);iterations.append(other.iterations)
    jac=factor.jacobian
    arrays=dict(u=base.u,f=base.f,g=base.g,
        base_force=2*weight*m*(d/epsilon-base.f),base_current=-weight*value.link_derivative_alpha,
        base_energy=np.array(weight*(value.energy+np.dot(m,abs(d)**2)/epsilon)),
        direction_du=np.array([a.du for a in tangent]),direction_df=np.array([a.df for a in tangent]),
        hessian_actions=np.array([2*weight*m*(v/epsilon-a.df) for v,a in zip(directions,tangent)]),
        current_actions=np.array([weight*a.current_derivative for a in tangent]),
        nonlinear_forces=force_states,nonlinear_currents=current_states,nonlinear_energies=energies,
        jacobian_data=jac.data,jacobian_indices=jac.indices,jacobian_indptr=jac.indptr,
        jacobian_shape=np.array(jac.shape),free_real_components=factor.components)
    record=dict(id=f'n{n:04d}',n=n,epsilon=epsilon,worker_seconds=time.monotonic()-started,
        worker_pid=os.getpid(),spectral_residual=base.residual,spectral_iterations=base.iterations,
        factor_nnz_L=int(factor.factor.L.nnz),factor_nnz_U=int(factor.factor.U.nnz),
        maximum_tangent_equation_residual=max(a.equation_residual for a in tangent),
        maximum_perturbed_spectral_residual=max(residuals),maximum_perturbed_iterations=max(iterations),
        nonlinear_roots=1+len(residuals),factor_reused_for_directions=len(directions),
        scope='Full spatial DOFs; two input directions test this operator, not a two-mode physical model')
    return record,arrays


def aggregate(plan,output,records):
    graph,d,xy,radius,winding,directions=fields(plan)
    if sorted(r['n'] for r in records)!=list(range(plan['matsubara_count'])):
        raise ValueError('A full finite Matsubara sum requires every registered frequency')
    t=plan['T_K']/plan['Tc_K'];m=graph.area_weights;count=len(directions)
    G=2*m*d*np.log(t);I=np.zeros(len(graph.edges));F=float(np.dot(m,abs(d)**2)*np.log(t))
    H=2*m*directions*np.log(t);J=np.zeros((count,len(graph.edges)))
    shape=(count,len(plan['difference_steps']),2)
    gradients=np.zeros(shape+(len(d),),complex);currents=np.zeros(shape+(len(graph.edges),));energies=np.zeros(shape)
    for j,v in enumerate(directions):
        for q,h in enumerate(plan['difference_steps']):
            for p,sign in enumerate((-1,1)):
                dd=d+sign*h*v
                gradients[j,q,p]=2*m*dd*np.log(t)
                energies[j,q,p]=np.dot(m,abs(dd)**2)*np.log(t)
    for row in sorted(records,key=lambda r:r['n']):
        with np.load(output/row['fields_path']) as a:
            G+=a['base_force'];I+=a['base_current'];F+=float(a['base_energy'])
            H+=a['hessian_actions'];J+=a['current_actions']
            gradients+=a['nonlinear_forces'];currents+=a['nonlinear_currents'];energies+=a['nonlinear_energies']
    options=dict(Tc_K=plan['Tc_K'],T_K=plan['T_K'],tau_ee_Tc_ps=plan['tau_ee_Tc_ps'],tau_ep_Tc_ps=plan['tau_ep_Tc_ps'])
    model=ThermalKWTNormal(graph,d,**options);base=model.response(G,I)
    free=model.free;core=free&(radius<=plan['probe_radius_ell0'])
    edge_mid=(xy[graph.edges[:,0]]+xy[graph.edges[:,1]])/2
    edge_core=(graph.conductance>0)&(np.linalg.norm(edge_mid,axis=1)<=plan['probe_radius_ell0'])
    node_norm=lambda v:float(np.sqrt(np.sum(m[core]*abs(v[core])**2)))
    edge_norm=lambda v:float(np.sqrt(np.sum(abs(v[edge_core])**2/graph.conductance[edge_core])))
    rel=lambda numerator,denominator:float(numerator/denominator) if denominator>0 else None
    probe_records=[];rhs_actions=[];corrections=[];v_actions=[];fd_rhs=[]
    for j,(name,v) in enumerate(zip(plan['probes'],directions)):
        tangent=model.rhs_tangent(v,H[j],G,I,J[j])
        predicted=tangent['velocity_direction'];rhs_actions.append(predicted)
        corrections.append(tangent['baseline_mobility_gauge_correction']);v_actions.append(tangent['potential_direction'])
        ward=np.imag(np.conj(v)*G+np.conj(d)*H[j])+divergence(graph,J[j])
        comparisons=[];probe_fd=[]
        for q,h in enumerate(plan['difference_steps']):
            gh=(gradients[j,q,1]-gradients[j,q,0])/(2*h)
            ih=(currents[j,q,1]-currents[j,q,0])/(2*h)
            endpoints=[ThermalKWTNormal(graph,d+sign*h*v,**options).response(gradients[j,q,p],currents[j,q,p])
                for p,sign in enumerate((-1,1))]
            rhs=(endpoints[1]['velocity']-endpoints[0]['velocity'])/(2*h);probe_fd.append(rhs)
            energy_derivative=(energies[j,q,1]-energies[j,q,0])/(2*h)
            comparisons.append(dict(step=h,
                hessian_relative_error=rel(node_norm((gh-H[j])/m),node_norm(H[j]/m)),
                current_relative_error=rel(edge_norm(ih-J[j]),edge_norm(J[j])),
                rhs_relative_error=rel(node_norm(rhs-predicted),node_norm(predicted)),
                force_norm=node_norm(H[j]/m),current_norm=edge_norm(J[j]),rhs_norm=node_norm(predicted),
                energy_directional_difference=float(energy_derivative-np.real(np.vdot(G,v))),
                maximum_endpoint_power_residual=max(abs(z['dissipation_residual']) for z in endpoints),
                maximum_endpoint_continuity=max(z['continuity_max'] for z in endpoints)))
        fd_rhs.append(probe_fd)
        probe_records.append(dict(probe=name,comparisons=comparisons,
            maximum_free_linearized_noether=float(np.max(abs(ward[free]),initial=0.)),
            maximum_free_linearized_continuity=tangent['continuity_direction_max'],
            baseline_correction_norm=node_norm(tangent['baseline_mobility_gauge_correction']),
            rhs_total_norm=node_norm(predicted),
            baseline_correction_relative=rel(node_norm(tangent['baseline_mobility_gauge_correction']),node_norm(predicted))))
    hmatrix=np.array([[np.real(np.vdot(v,H[j])) for j in range(count)] for v in directions])
    np.savez_compressed(output/'full_node_operator_checks.npz',coordinates_bar=xy,d=d,area_weights=m,
        edges=graph.edges,conductance=graph.conductance,boundary_nodes=graph.boundary_nodes,
        gap_gradient=G,current=I,baseline_velocity=base['velocity'],baseline_potential_v=base['potential_v'],
        input_directions=directions,hessian_actions=H,current_actions=J,rhs_actions=np.asarray(rhs_actions),
        rhs_baseline_corrections=np.asarray(corrections),potential_actions=np.asarray(v_actions),
        finite_difference_rhs=np.asarray(fd_rhs),difference_steps=np.asarray(plan['difference_steps']),
        nonlinear_energies=energies)
    smallest=[p['comparisons'][-1] for p in probe_records]
    calculus_flags={key:all(p[key] is not None and p[key]<=plan['calculus_relative_tolerance'] for p in smallest)
        for key in ('hessian_relative_error','current_relative_error','rhs_relative_error')}
    return dict(status='FULL_NODE_THERMAL_OPERATOR_VALIDATION_COMPLETE',
        base_free_energy=F,base_gradient_density_core_norm=node_norm(G/m),
        base_velocity_core_norm=node_norm(base['velocity']),
        baseline={key:value for key,value in base.items() if not isinstance(value,np.ndarray)},
        tD_ps=model.tD_ps,taupsi_ps=model.taupsi_ps,probes=probe_records,
        tested_hessian_bilinear_matrix=hmatrix.tolist(),maximum_bilinear_asymmetry=float(np.max(abs(hmatrix-hmatrix.T))),
        calculus_relative_tolerance=plan['calculus_relative_tolerance'],calculus_flags=calculus_flags,
        all_calculus_flags_met=all(calculus_flags.values()),
        scope='Full-node affine thermal RHS/Jacobian prepared. Finite-sum free energy, fixed T and contacts. No time trajectory, nonequilibrium work law, circuit or photon.',
        production_changed=False,physical_time_steps=0,stage4_complete=False)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--plan',type=Path,default=DEFAULT_PLAN)
    parser.add_argument('--output-root',type=Path,required=True)
    parser.add_argument('--pilot',action='store_true');parser.add_argument('--execute',action='store_true')
    args=parser.parse_args();plan=read(args.plan)
    if plan['schema']!='pysnspd.stage4.thermal_weak_operator.v1':raise ValueError('Unknown plan schema')
    sources={**plan['frozen_sources'],**{name:sha(ROOT/name) for name in plan['new_sources']}}
    for name,digest in {**plan['frozen_sources'],**plan['inputs']}.items():
        if sha(ROOT/name)!=digest:raise ValueError('Source/input changed: '+name)
    modes=plan['pilot_modes'] if args.pilot else list(range(plan['matsubara_count']))
    jobs=[dict(plan=plan,n=n) for n in modes]
    resources=linux_resources();budget=resource_budget(resources,.9,.9,max_workers=len(jobs))
    if args.output_root.exists():raise FileExistsError('Use a fresh output directory; no implicit resume or overwrite')
    parent=args.output_root.resolve().parent
    while not parent.exists():parent=parent.parent
    available=shutil.disk_usage(parent).free
    reserve=max(128*1024**2,len(jobs)*8*1024**2)
    if reserve>.9*available:raise RuntimeError('Output reservation exceeds90% of available storage')
    identity=dict(plan_sha256=sha(args.plan),sources=sources,inputs=plan['inputs'],resources=resources,budget=budget,
        frequency_tasks=len(jobs),nonlinear_roots=9*len(jobs),factorizations=len(jobs),pilot=args.pilot,
        output_reservation_bytes=reserve,available_storage_bytes=available,physical_time_steps=0,
        factor_storage='Sparse Jacobian CSC and solved u saved; LU factors reused in each task and rebuilt from saved matrices by a future temporal driver')
    if not args.execute:print(json.dumps(dict(status='DRY_RUN_NO_WRITES',**identity),indent=2));return
    args.output_root.mkdir(parents=True);(args.output_root/'modes').mkdir()
    write(args.output_root/'identity.json',identity);shutil.copyfile(args.plan,args.output_root/'executed_plan.json')
    start=time.monotonic();records=[];affinity=os.sched_getaffinity(0)
    log=(args.output_root/'progress.jsonl').open('x',encoding='utf8',buffering=1)
    def event(name,**values):
        line=json.dumps(dict(event=name,elapsed_seconds=time.monotonic()-start,**values))
        print(line,flush=True);log.write(line+'\n')
    try:
        event('START',tasks=len(jobs),budget=budget)
        os.sched_setaffinity(0,{budget['coordinator_cpu']});context=mp.get_context('spawn');counter=context.Value('i',0)
        with ProcessPoolExecutor(max_workers=budget['workers'],mp_context=context,initializer=initialize_affinity,
            initargs=(budget['worker_affinity_cpus'],counter)) as pool:
            futures=[pool.submit(worker,job) for job in jobs]
            try:
                for future in as_completed(futures):
                    row,arrays=future.result();path=args.output_root/'modes'/(row['id']+'.npz')
                    np.savez_compressed(path,**arrays);row.update(fields_path=path.relative_to(args.output_root).as_posix(),fields_sha256=sha(path))
                    write(path.with_suffix('.json'),row);records.append(row)
                    elapsed=time.monotonic()-start;fraction=len(records)/len(jobs)
                    print(f"[{'#'*int(24*fraction)}{'-'*(24-int(24*fraction))}] Matsubara {len(records)}/{len(jobs)} | {elapsed:.1f} s | ETA {elapsed*(1-fraction)/fraction:.1f} s",flush=True)
                    event('PROGRESS',complete=len(records),total=len(jobs),mode=row['n'],worker_seconds=row['worker_seconds'],eta_seconds=elapsed*(1-fraction)/fraction)
            except BaseException:
                for future in futures:future.cancel()
                raise
        result=(dict(status='THERMAL_WEAK_COST_PILOT_COMPLETE',scope='Incomplete frequency subset; no total Hessian or physical result')
            if args.pilot else aggregate(plan,args.output_root,records))
        result.update(runtime_seconds=time.monotonic()-start,records=records,physical_time_steps=0,production_changed=False)
        write(args.output_root/'summary.json',result);event('COMPLETE',status=result['status'],runtime_seconds=result['runtime_seconds'])
    except BaseException as exc:
        write(args.output_root/'failure.json',dict(error=type(exc).__name__,reason=str(exc),
            policy='Stop with completed checkpoints preserved; no fallback, clipping, retry or overwrite'))
        event('FAILED',error=type(exc).__name__,reason=str(exc));raise
    finally:
        os.sched_setaffinity(0,affinity);log.close()


if __name__=='__main__':main()
