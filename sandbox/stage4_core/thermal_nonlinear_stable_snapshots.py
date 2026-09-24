"""Nonlinear thermal constitutive checks at saved weak-trajectory states.

Parallelism is across independent Matsubara frequencies. Each frequency reuses
its original spectrum and Jacobian for all registered snapshots. No trajectory
is integrated again and no material, boundary or time scale is changed.
"""
from __future__ import annotations
import argparse
from concurrent.futures import ProcessPoolExecutor,as_completed
import hashlib,json,multiprocessing as mp,os
from pathlib import Path
import shutil,sys,time
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]
sys.path.insert(0,str(HERE))
from parallel_runtime import initialize_affinity,limit_thread_environment,linux_resources,resource_budget
limit_thread_environment();sys.path.insert(0,str(ROOT))
import numpy as np
from pysnspd.experimental import thermal_spatial_usadel as thermal
from pysnspd.experimental import thermal_stable_newton as stable
from pysnspd.experimental.thermal_time_operator import SavedSpectralMode
from pysnspd.experimental.thermal_weak_response import ThermalKWTNormal
from pysnspd.experimental.thermal_snapshot import spectral_difference
_CACHE={}


def read(path):return json.loads(Path(path).read_text())
def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def write(path,value):
    with Path(path).open('x',encoding='utf8') as stream:
        json.dump(value,stream,indent=2,allow_nan=False);stream.write('\n')


def snapshot_fields(plan,operator_root,time_root):
    key=(operator_root,time_root)
    if key not in _CACHE:
        with np.load(Path(operator_root)/'full_node_operator_checks.npz') as a:
            graph=thermal.ThermalGraph(a['area_weights'],a['edges'],a['conductance'],a['coordinates_bar'],a['boundary_nodes'])
            d0,G0,I0=a['d'],a['gap_gradient'],a['current']
        xs=[];hs=[];currents=[];velocities=[];energy_taylor=[]
        summary=read(Path(time_root)/'refined/summary.json')
        for state in plan['states']:
            record=next(r for r in summary['observations'] if r['time_ps']==state['time_ps'])
            if record['fields_path']!=state['snapshot_path']:raise ValueError('Snapshot time/path changed')
            with np.load(Path(time_root)/state['snapshot_path']) as a:
                if not np.array_equal(a['d0'],d0) or not np.array_equal(a['G0'],G0):raise ValueError('Different reference core')
                prefix=state['kind']+'_'
                xs.append(a[prefix+'displacement']);hs.append(a[prefix+'hessian'])
                currents.append(a[prefix+'current_increment']);velocities.append(a[prefix+'velocity'])
            energy_taylor.append(next(r['quadratic_free_energy'] for r in record['states'] if r['name']==state['kind']))
        _CACHE[key]=(graph,d0,G0,I0,np.array(xs),np.array(hs),np.array(currents),np.array(velocities),np.array(energy_taylor))
    return _CACHE[key]


def worker(job):
    started=time.monotonic();plan=job['plan'];row=job['mode']
    graph,d0,G0,I0,xs,_,_,_,_=snapshot_fields(plan,job['operator_root'],job['time_root'])
    t=plan['T_K']/plan['Tc_K'];weight=2*np.pi*t;epsilon=row['epsilon']
    with np.load(Path(job['operator_root'])/row['fields_path']) as a:
        u0=a['u'].copy();mode=SavedSpectralMode(graph,a,epsilon,t)
    n=graph.n_nodes;count=len(xs)
    rhs=np.repeat(graph.area_weights,2)[:,None]*np.stack((xs.real,xs.imag),axis=-1).reshape(count,2*n).T
    du=np.zeros((2*n,count));du[mode.components]=mode.factor.solve(rhs[mode.components])
    guesses=u0[None,:]+(du[::2]+1j*du[1::2]).T
    gradients=np.zeros((count,n),complex);currents=np.zeros((count,len(graph.edges)));energies=np.zeros(count)
    spectral_u=[];checks=[]
    for index,(state,x,guess) in enumerate(zip(plan['states'],xs,guesses)):
        if state['reuse_baseline']:
            if np.any(x!=0):raise ValueError('A reused initial baseline must have zero displacement')
            u=u0;residual=row['spectral_residual'];iterations=0
        else:
            solution=stable.solve_frequency(graph,d0+x,epsilon,fixed_nodes=graph.boundary_nodes,
                fixed_u=u0[graph.boundary_nodes],initial_u=guess,tol=plan['spectral_tolerance'],max_iterations=plan['maximum_newton_iterations'])
            u=solution.u;residual=solution.residual;iterations=solution.iterations
        diff=spectral_difference(graph,d0,u0,d0+x,u,epsilon)
        gradients[index]=weight*diff['gap_gradient_difference']
        currents[index]=weight*diff['current_difference'];energies[index]=weight*diff['renormalized_energy_difference']
        spectral_u.append(u)
        checks.append(dict(state_id=state['id'],spectral_residual=float(residual),newton_iterations=int(iterations),reused=state['reuse_baseline']))
    return dict(id=row['id'],n=row['n'],epsilon=epsilon,worker_seconds=time.monotonic()-started,
        nonlinear_roots=sum(not s['reuse_baseline'] for s in plan['states']),checks=checks),dict(delta_gradients=gradients,
        delta_currents=currents,delta_energies=energies,spectral_u=np.array(spectral_u))


def aggregate(plan,operator_root,time_root,output,records):
    graph,d0,G0,I0,xs,H,I,taylor_rhs,taylor_energies=snapshot_fields(plan,str(operator_root),str(time_root))
    mass=graph.area_weights;t=plan['T_K']/plan['Tc_K'];count=len(xs)
    dG=2*mass[None,:]*np.log(t)*xs;dI=np.zeros((count,len(graph.edges)))
    dF=np.array([np.dot(mass,2*np.real(np.conj(d0)*x)+abs(x)**2)*np.log(t) for x in xs])
    for row in sorted(records,key=lambda r:r['n']):
        with np.load(output/row['fields_path']) as a:
            dG+=a['delta_gradients'];dI+=a['delta_currents'];dF+=a['delta_energies']
    options={key:plan[key] for key in ('T_K','Tc_K','tau_ee_Tc_ps','tau_ep_Tc_ps')}
    nonlinear_rhs=[];potential=[];state_metrics=[]
    core=np.linalg.norm(graph.coordinates_bar,axis=1)<=plan['probe_radius_ell0']
    edge_mid=(graph.coordinates_bar[graph.edges[:,0]]+graph.coordinates_bar[graph.edges[:,1]])/2
    active=(graph.conductance>0)&(np.linalg.norm(edge_mid,axis=1)<=plan['probe_radius_ell0'])
    norm_node=lambda z:float(np.sqrt(np.sum(mass[core]*abs(z[core])**2)))
    norm_edge=lambda z:float(np.sqrt(np.sum(abs(z[active])**2/graph.conductance[active])))
    base_energy=read(operator_root/'summary.json')['base_free_energy']
    for state,x,g,current,energy,grhs in zip(plan['states'],xs,dG,dI,dF,taylor_rhs):
        value=ThermalKWTNormal(graph,d0+x,**options).response(G0+g,I0+current)
        nonlinear_rhs.append(value['velocity']);potential.append(value['potential_v'])
        state_metrics.append(dict(state_id=state['id'],time_ps=state['time_ps'],kind=state['kind'],
            exact_energy_change=float(energy),quadratic_energy_change=float(taylor_energies[len(state_metrics)]-base_energy),
            rhs_difference_L2=norm_node(value['velocity']-grhs),exact_rhs_L2=norm_node(value['velocity']),
            free_energy_rate=value['free_energy_rate'],kwt_loss=value['kwt_loss'],normal_loss=value['normal_loss'],
            power_residual=value['dissipation_residual'],continuity_max=value['continuity_max'],noether_max=value['noether_max']))
    nonlinear_rhs=np.array(nonlinear_rhs);potential=np.array(potential);comparisons=[]
    for index,state in enumerate(plan['states']):
        if state['kind']=='baseline':continue
        base=next(j for j,s in enumerate(plan['states']) if s['kind']=='baseline' and s['time_ps']==state['time_ps'])
        initial=next(j for j,s in enumerate(plan['states']) if s['kind']==state['kind'] and s['time_ps']==0)
        initial_base=next(j for j,s in enumerate(plan['states']) if s['kind']=='baseline' and s['time_ps']==0)
        fields={
            'rhs':(nonlinear_rhs[index]-nonlinear_rhs[base],taylor_rhs[index]-taylor_rhs[base],nonlinear_rhs[initial]-nonlinear_rhs[initial_base],norm_node),
            'force_density':((dG[index]-dG[base])/mass,(H[index]-H[base])/mass,(dG[initial]-dG[initial_base])/mass,norm_node),
            'current':(dI[index]-dI[base],I[index]-I[base],dI[initial]-dI[initial_base],norm_edge),
            'phase_torque_density':(np.imag(np.conj(d0+xs[index])*(G0+dG[index])-np.conj(d0+xs[base])*(G0+dG[base]))/mass,
                np.imag(np.conj(d0)*(H[index]-H[base])+np.conj(xs[index]-xs[base])*G0)/mass,
                np.imag(np.conj(d0+xs[initial])*(G0+dG[initial])-np.conj(d0+xs[initial_base])*(G0+dG[initial_base]))/mass,norm_node)}
        for kind,(exact,linear,initial_value,norm) in fields.items():
            error=norm(exact-linear);remaining=norm(exact);initial_norm=norm(initial_value)
            comparisons.append(dict(state_id=state['id'],time_ps=state['time_ps'],probe=state['kind'],observable=kind,
                difference_norm=error,exact_response_norm=remaining,initial_response_norm=initial_norm,
                relative_to_response=error/remaining if remaining>plan['near_zero_norm'] else None,
                relative_to_initial=error/initial_norm if initial_norm>plan['near_zero_norm'] else None,
                within_registered_weak_margin=bool(error<=plan['near_zero_norm']+plan['weak_response_relative_margin']*max(remaining,initial_norm))))
        energy_exact=dF[index]-dF[base];energy_linear=taylor_energies[index]-taylor_energies[base]
        initial_energy=dF[initial]-dF[initial_base]
        comparisons.append(dict(state_id=state['id'],time_ps=state['time_ps'],probe=state['kind'],observable='energy_change',
            exact_response=float(energy_exact),quadratic_response=float(energy_linear),difference=float(energy_exact-energy_linear),
            relative_to_initial=abs(float((energy_exact-energy_linear)/initial_energy)) if abs(initial_energy)>plan['near_zero_norm'] else None))
    path=output/'nonlinear_comparison_fields.npz'
    np.savez_compressed(path,state_ids=np.array([s['id'] for s in plan['states']]),d0=d0,G0=G0,I0=I0,displacements=xs,
        delta_gradients_exact=dG,delta_currents_exact=dI,delta_energies_exact=dF,
        delta_gradients_taylor=H,delta_currents_taylor=I,rhs_taylor=taylor_rhs,rhs_exact=nonlinear_rhs,potential_exact=potential,
        coordinates_bar=graph.coordinates_bar,area_weights=mass,edges=graph.edges,conductance=graph.conductance,boundary_nodes=graph.boundary_nodes)
    return dict(status='NONLINEAR_THERMAL_SNAPSHOTS_COMPLETE',states=state_metrics,comparisons=comparisons,
        comparison_fields_path=path.name,comparison_fields_sha256=sha(path),
        all_registered_response_margins_met=all(row.get('within_registered_weak_margin',True) for row in comparisons),
        scope='Nonlinear constitutive evaluations at saved trajectory states; no new trajectory. Response comparisons subtract exact same-time baseline.',
        production_changed=False,nonthermal_work_closed=False,detector_validated=False)


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--plan',type=Path,required=True)
    p.add_argument('--operator-root',type=Path,required=True);p.add_argument('--time-root',type=Path,required=True)
    p.add_argument('--output-root',type=Path,required=True);p.add_argument('--execute',action='store_true');args=p.parse_args()
    plan=read(args.plan);modes=read(args.operator_root/'summary.json')['records']
    if sorted(row['n'] for row in modes)!=list(range(plan['matsubara_count'])):raise ValueError('Incomplete operator frequencies')
    for name,digest in plan['frozen_sources'].items():
        if sha(ROOT/name)!=digest:raise ValueError('Frozen source changed: '+name)
    inputs={}
    for base,mapping in ((args.operator_root,plan['operator_inputs']),(args.time_root,plan['time_inputs'])):
        for name,digest in mapping.items():
            if sha(base/name)!=digest:raise ValueError('Frozen admitted input changed: '+name)
            inputs[str(base/name)]=digest
    for row in modes:
        if sha(args.operator_root/row['fields_path'])!=row['fields_sha256']:raise ValueError('Mode checkpoint changed')
        inputs[str(args.operator_root/row['fields_path'])]=row['fields_sha256']
    resources=linux_resources();budget=resource_budget(resources,.9,.9,max_workers=27)
    if args.output_root.exists():raise FileExistsError('Fresh output required; no overwrite or implicit resume')
    parent=args.output_root.resolve().parent
    while not parent.exists():parent=parent.parent
    if plan['disk_reserve_bytes']>.9*shutil.disk_usage(parent).free:raise RuntimeError('Insufficient output volume')
    sources={**plan['frozen_sources'],**{name:sha(ROOT/name) for name in plan['new_sources']}}
    identity=dict(plan_sha256=sha(args.plan),sources=sources,inputs=inputs,resources=resources,budget=budget,
        frequency_tasks=len(modes),new_spectral_roots=len(modes)*sum(not s['reuse_baseline'] for s in plan['states']),
        previous_time_steps_repeated=0,changes_to_physical_parameters=False)
    if not args.execute:print(json.dumps(dict(status='DRY_RUN_NO_WRITES',**identity),indent=2));return
    args.output_root.mkdir(parents=True);(args.output_root/'modes').mkdir();write(args.output_root/'identity.json',identity)
    shutil.copyfile(args.plan,args.output_root/'executed_plan.json');start=time.monotonic();records=[];affinity=os.sched_getaffinity(0)
    log=(args.output_root/'progress.jsonl').open('x',buffering=1)
    def event(name,**values):
        line=json.dumps(dict(event=name,elapsed_seconds=time.monotonic()-start,**values));print(line,flush=True);log.write(line+'\n')
    try:
        os.sched_setaffinity(0,{budget['coordinator_cpu']});context=mp.get_context('spawn');counter=context.Value('i',0)
        with ProcessPoolExecutor(max_workers=budget['workers'],mp_context=context,initializer=initialize_affinity,
            initargs=(budget['worker_affinity_cpus'],counter)) as pool:
            futures=[pool.submit(worker,dict(plan=plan,mode=row,operator_root=str(args.operator_root.resolve()),time_root=str(args.time_root.resolve()))) for row in modes]
            for future in as_completed(futures):
                row,arrays=future.result();path=args.output_root/'modes'/(row['id']+'.npz');np.savez_compressed(path,**arrays)
                row.update(fields_path=path.relative_to(args.output_root).as_posix(),fields_sha256=sha(path));write(path.with_suffix('.json'),row);records.append(row)
                elapsed=time.monotonic()-start;fraction=len(records)/len(modes)
                print(f"[{'#'*int(24*fraction)}{'-'*(24-int(24*fraction))}] {len(records)}/{len(modes)} | {elapsed:.1f}s | ETA {elapsed*(1-fraction)/fraction:.1f}s",flush=True)
                event('PROGRESS',complete=len(records),total=len(modes),mode=row['n'],eta_seconds=elapsed*(1-fraction)/fraction)
        result=aggregate(plan,args.operator_root,args.time_root,args.output_root,records)
        result.update(runtime_seconds=time.monotonic()-start,records=records)
        write(args.output_root/'summary.json',result);event('COMPLETE',runtime_seconds=result['runtime_seconds'],status=result['status'])
    except BaseException as exc:
        write(args.output_root/'failure.json',dict(type=type(exc).__name__,reason=str(exc),policy='Retain completed frequency outputs; no retry, fallback or changed tolerance'))
        raise
    finally:os.sched_setaffinity(0,affinity);log.close()


if __name__=='__main__':main()
