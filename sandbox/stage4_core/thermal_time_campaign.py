"""User-run full-node weak thermal affine evolution from saved spectral factors.

The physical tangent remains full spatial dimension. One resident worker pool
owns disjoint frequency shards and reuses their LU factors for every Krylov
action. A planned refined trajectory measures observable variation; no claim
of exact nonlinear free-energy conservation or detector evolution is made.
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

HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]
sys.path.insert(0,str(HERE))
from parallel_runtime import linux_resources,resource_budget,limit_thread_environment
limit_thread_environment()
sys.path.insert(0,str(ROOT))
import numpy as np
from pysnspd.experimental import thermal_spatial_usadel as thermal
from pysnspd.experimental.thermal_weak_response import ThermalKWTNormal,divergence
from pysnspd.experimental.thermal_time_operator import SavedSpectralMode,FullNodeCoordinates
from pysnspd.experimental.affine_krylov import advance


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def read(path):return json.loads(Path(path).read_text())
def write(path,value):
    with Path(path).open('x',encoding='utf8') as stream:
        json.dump(value,stream,indent=2,allow_nan=False);stream.write('\n')


def actor(pipe,cpu,graph,rows,source,t):
    try:
        os.sched_setaffinity(0,{cpu});started=time.monotonic();modes=[]
        for row in rows:
            with np.load(Path(source)/row['fields_path']) as arrays:
                modes.append(SavedSpectralMode(graph,arrays,row['epsilon'],t))
        pipe.send(dict(event='READY',modes=len(modes),seconds=time.monotonic()-started,cpu=cpu))
        while True:
            command,value=pipe.recv()
            if command=='STOP':break
            if command!='APPLY':raise ValueError('Unknown resident-worker command')
            directions=np.asarray(value,complex)
            H=np.zeros_like(directions);I=np.zeros((len(directions),len(graph.edges)));worst=0.
            for mode in modes:
                for j,direction in enumerate(directions):
                    h,current,residual=mode.apply(direction)
                    H[j]+=h;I[j]+=current;worst=max(worst,residual)
            pipe.send((H,I,worst))
    except BaseException as exc:
        try:pipe.send(dict(event='ERROR',type=type(exc).__name__,reason=str(exc),traceback=traceback.format_exc()))
        except (BrokenPipeError,EOFError):pass
    finally:pipe.close()


class ResidentSpectrum:
    def __init__(self,graph,rows,source,t,budget,event):
        self.graph=graph;self.t=t;self.actions=0;self.maximum_linear_solve_residual=0.;self.processes=[];self.pipes=[]
        context=mp.get_context('spawn');count=budget['workers'];ordered=sorted(rows,key=lambda r:r['n'])
        # Frequency ownership is fixed: no task migration or duplicated factors.
        for index,cpu in enumerate(budget['worker_affinity_cpus']):
            parent,child=context.Pipe();shard=ordered[index::count]
            process=context.Process(target=actor,args=(child,cpu,graph,shard,str(source),t))
            process.start();child.close();self.processes.append(process);self.pipes.append(parent)
        try:
            for index,pipe in enumerate(self.pipes):
                message=pipe.recv()
                if message.get('event')!='READY':raise RuntimeError(str(message))
                event('FACTORS_READY',worker=index,**{k:v for k,v in message.items() if k!='event'})
        except BaseException:
            self.close();raise

    def apply(self,directions):
        values=np.asarray(directions,complex)
        if values.ndim==1:values=values[None,:]
        for pipe in self.pipes:pipe.send(('APPLY',values))
        H=2*self.graph.area_weights[None,:]*np.log(self.t)*values
        I=np.zeros((len(values),len(self.graph.edges)))
        for pipe in self.pipes:
            result=pipe.recv()
            if isinstance(result,dict):raise RuntimeError(str(result))
            h,current,residual=result;H+=h;I+=current
            self.maximum_linear_solve_residual=max(self.maximum_linear_solve_residual,residual)
        self.actions+=1
        return H,I

    def close(self):
        for process,pipe in zip(self.processes,self.pipes):
            if process.is_alive():
                try:pipe.send(('STOP',None))
                except (BrokenPipeError,EOFError):pass
        for process in self.processes:
            process.join(timeout=5)
            if process.is_alive():process.terminate();process.join(timeout=5)
        for pipe in self.pipes:pipe.close()


def metrics(graph,d0,G0,I0,base_energy,model_options,spectrum,states,names):
    H,I=spectrum.apply(states);mass=graph.area_weights
    records=[];arrays={}
    for name,x,h,current in zip(names,states,H,I):
        model=ThermalKWTNormal(graph,d0,**model_options)
        b=model.response(G0,I0)['velocity']
        tangent=model.rhs_tangent(x,h,G0,I0,current)
        velocity=b+tangent['velocity_direction']
        force=G0+h;total_current=I0+current;d=d0+x
        constitutive=ThermalKWTNormal(graph,d,**model_options).response(force,total_current)
        F=base_energy+np.real(np.vdot(G0,x))+.5*np.real(np.vdot(x,h))
        rate=float(np.real(np.vdot(force,velocity)))
        # O(x^2) residual of the thermodynamic/constitutive Taylor truncation.
        truncation=velocity-constitutive['velocity']
        free=model.free
        item=dict(name=name,quadratic_free_energy=float(F),quadratic_free_energy_rate=rate,
            approximate_kwt_loss=constitutive['kwt_loss'],approximate_normal_loss=constitutive['normal_loss'],
            quadratic_rate_plus_approximate_losses=rate+constitutive['kwt_loss']+constitutive['normal_loss'],
            constitutive_taylor_rhs_defect_L2=float(np.sqrt(np.sum(mass*abs(truncation)**2))),
            full_affine_rhs_L2=float(np.sqrt(np.sum(mass*abs(velocity)**2))),
            maximum_linearized_continuity=tangent['continuity_direction_max'],
            approximate_noether_residual=constitutive['noether_max'],
            displacement_L2=float(np.sqrt(np.sum(mass*abs(x)**2))),
            maximum_displacement=float(np.max(abs(x[free]))),
            baseline_mobility_gauge_correction_L2=float(np.sqrt(np.sum(mass*abs(tangent['baseline_mobility_gauge_correction'])**2))))
        records.append(item)
        for key,value in dict(displacement=x,velocity=velocity,hessian=h,current_increment=current,
            potential_increment=tangent['potential_direction'],constitutive_taylor_rhs_defect=truncation).items():arrays[name+'_'+key]=value
    return records,arrays


def run_pass(label,settings,plan,graph,d0,G0,I0,initial_directions,base_energy,spectrum,output,event):
    coordinates=FullNodeCoordinates(graph,plan['gap_reference_kBTc']);size=coordinates.size
    options={key:plan[key] for key in ('Tc_K','T_K','tau_ee_Tc_ps','tau_ep_Tc_ps')}
    model=ThermalKWTNormal(graph,d0,**options);base=model.response(G0,I0)
    deltas=plan['perturbation_amplitude']*initial_directions
    references=np.array([max(np.linalg.norm(coordinates.encode(z)) for z in deltas)]+
        [np.linalg.norm(coordinates.encode(z)) for z in deltas])
    if np.any(references<=0):raise ValueError('Nonzero registered perturbations required')
    initial=np.r_[np.concatenate([np.zeros(size)]+[coordinates.encode(z)/scale for z,scale in zip(deltas,references[1:])]),1.]
    baseline_encoded=coordinates.encode(base['velocity'])/references[0]
    def unpack(value):
        return np.array([coordinates.decode(block*scale) for block,scale in zip(value[:-1].reshape(3,size),references)])
    def action(value):
        states=unpack(value);H,I=spectrum.apply(states);result=[]
        for x,h,current,scale in zip(states,H,I,references):
            direction=model.rhs_tangent(x,h,G0,I0,current)['velocity_direction']
            result.append(coordinates.encode(direction)/scale)
        result[0]+=value[-1]*baseline_encoded
        return np.r_[np.concatenate(result),0.]
    folder=output/label;folder.mkdir();state=initial.copy();history=[];observations=[]
    times=np.array(plan['observation_times_ps'])
    if settings['midpoint_intervals']:
        times=np.sort(np.r_[times,(times[:-1]+times[1:])/2])
    start=time.monotonic();completed=0.
    def guard(value,time_ps,tag):
        blocks=unpack(value);states=np.array([blocks[0],blocks[0]+blocks[1],blocks[0]+blocks[2]])
        maximum=float(np.max(abs(states))/plan['gap_reference_kBTc'])
        if maximum>plan['maximum_relative_displacement']:
            np.savez_compressed(folder/('weak_domain_exit_'+tag+'.npz'),time_ps=time_ps,states=states,affine_state=value)
            raise RuntimeError(f'Declared weak domain exceeded at {time_ps:.6g} ps: {maximum:.6g}; no clipping or continuation')
        return blocks,states
    for index,target in enumerate(times):
        if index:
            left=float(times[index-1]);interval=(target-left)/model.tD_ps
            def progress(row,value):
                current_ps=left+(row['start']+row['step'])*model.tD_ps
                entry=dict(pass_name=label,interval=index,time_ps=current_ps,**row);history.append(entry)
                if value is not None:
                    guard(value,current_ps,str(len(history)))
                    np.savez_compressed(folder/f'accepted_{len(history):05d}.npz',time_ps=current_ps,affine_state=value)
                fraction=current_ps/plan['observation_times_ps'][-1]
                elapsed=time.monotonic()-start
                event('ARNOLDI_STEP',pass_name=label,time_ps=current_ps,accepted=row['accepted'],dimension=row['dimension'],
                    integrated_defect=row['integrated_defect'],tolerance=row['tolerance'],operator_batches=spectrum.actions,
                    eta_seconds=elapsed*(1-fraction)/fraction if fraction>0 else None)
            state,steps=advance(action,state,interval,rtol=settings['defect_rtol'],atol=settings['defect_atol'],
                max_dimension=settings['maximum_krylov_dimension'],maximum_steps=plan['maximum_restarts_per_interval'],
                active_size=len(state)-1,progress=progress)
        blocks,states=guard(state,float(target),'observation_'+str(index))
        if target not in plan['observation_times_ps']:continue
        rec,arrays=metrics(graph,d0,G0,I0,base_energy,options,spectrum,states,['baseline','amplitude','angular_phase'])
        arrays.update(d0=d0,G0=G0,coordinates_bar=graph.coordinates_bar,area_weights=graph.area_weights,edges=graph.edges,
            conductance=graph.conductance,boundary_nodes=graph.boundary_nodes,
            amplitude_difference=blocks[1],angular_phase_difference=blocks[2])
        path=folder/f'observation_{len(observations):03d}.npz';np.savez_compressed(path,**arrays)
        observations.append(dict(time_ps=float(target),states=rec,fields_path=path.relative_to(output).as_posix(),fields_sha256=sha(path)))
        event('OBSERVATION',pass_name=label,time_ps=float(target),count=len(observations),total=len(plan['observation_times_ps']))
    result=dict(runtime_seconds=time.monotonic()-start,observations=observations,history=history,
        scaled_reference_norms=references.tolist(),maximum_accepted_defect=max((r['integrated_defect'] for r in history if r['accepted']),default=0.))
    write(folder/'summary.json',result);return result


def compare_passes(plan,graph,output,first,second):
    records=[];mass=graph.area_weights
    for left,right in zip(first['observations'],second['observations']):
        if left['time_ps']!=right['time_ps']:raise ValueError('Different physical observation horizon')
        with np.load(output/left['fields_path']) as a,np.load(output/right['fields_path']) as b:
            for probe in ('amplitude','angular_phase'):
                def field(data,kind):
                    if kind=='displacement':return data[probe+'_difference']
                    if kind=='current':return data[probe+'_current_increment']-data['baseline_current_increment']
                    h=data[probe+'_hessian']-data['baseline_hessian']
                    if kind=='force_density':return h/mass
                    return np.imag(np.conj(data['d0'])*h+np.conj(data[probe+'_difference'])*data['G0'])/mass
                for kind in ('displacement','current','force_density','phase_torque_density'):
                    active=graph.conductance>0
                    norm=(lambda z:float(np.sqrt(np.sum(abs(z[active])**2/graph.conductance[active])))) if kind=='current' else (lambda z:float(np.sqrt(np.sum(mass*abs(z)**2))))
                    error=norm(field(a,kind)-field(b,kind));reference=norm(field(b,kind))
                    initial=first['observations'][0]
                    with np.load(output/initial['fields_path']) as initial_fields:scale=norm(field(initial_fields,kind))
                    records.append(dict(time_ps=left['time_ps'],probe=probe,observable=kind,error_L2=error,
                        relative_to_response=error/reference if reference>1e-14 else None,
                        relative_to_initial=error/scale if scale>1e-14 else None,
                        admitted=error<=plan['refinement_absolute_norm']+plan['refinement_relative_limit']*max(reference,scale)))
    return dict(records=records,all_comparisons_met=all(r['admitted'] for r in records),
        tolerance=plan['refinement_relative_limit'],meaning='Observed change under planned Krylov/timing refinement, not rigorous forward error or nonlinear detector accuracy')


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--plan',type=Path,required=True)
    p.add_argument('--operator-root',type=Path,required=True);p.add_argument('--output-root',type=Path,required=True)
    p.add_argument('--execute',action='store_true');args=p.parse_args()
    plan=read(args.plan);old=read(args.operator_root/'summary.json');old_identity=read(args.operator_root/'identity.json')
    if plan.get('link_field')!='alpha_zero':raise ValueError('This saved operator admits only the declared alpha=0 core')
    if sha(args.operator_root/'summary.json')!=plan['operator_summary_sha256']:
        raise ValueError('The admitted operator summary changed')
    if sha(args.operator_root/'identity.json')!=plan['operator_identity_sha256']:
        raise ValueError('The admitted operator identity changed')
    original=read(args.operator_root/'executed_plan.json')
    if old.get('status')!='FULL_NODE_THERMAL_OPERATOR_VALIDATION_COMPLETE' or not old.get('all_calculus_flags_met'):
        raise RuntimeError('The full thermal operator must first pass its registered validation')
    if sorted(row['n'] for row in old['records'])!=list(range(plan['matsubara_count'])):
        raise ValueError('Incomplete, duplicated or shifted finite spectral sum')
    for key in ('T_K','Tc_K','matsubara_count','tau_ee_Tc_ps','tau_ep_Tc_ps'):
        if original[key]!=plan[key]:raise ValueError('Physical/operator contract changed: '+key)
    sources={}
    for name,digest in old_identity['sources'].items():
        if sha(ROOT/name)!=digest:raise ValueError('Prepared operator source changed: '+name)
        sources[name]=digest
    for name in plan['new_sources']:sources[name]=sha(ROOT/name)
    inputs={}
    for row in old['records']:
        path=args.operator_root/row['fields_path']
        if sha(path)!=row['fields_sha256']:raise ValueError('Prepared mode changed')
        inputs[row['fields_path']]=row['fields_sha256']
    checks_path=args.operator_root/'full_node_operator_checks.npz';inputs[checks_path.name]=sha(checks_path)
    if inputs[checks_path.name]!=plan['operator_checks_sha256']:raise ValueError('Admitted aggregate operator fields changed')
    with np.load(checks_path) as a:
        graph=thermal.ThermalGraph(a['area_weights'],a['edges'],a['conductance'],a['coordinates_bar'],a['boundary_nodes'])
        d0,G0,I0,initial,Hexpected,Iexpected=[a[key].copy() for key in ('d','gap_gradient','current','input_directions','hessian_actions','current_actions')]
    resources=linux_resources();budget=resource_budget(resources,.9,.9,max_workers=plan['maximum_workers'])
    if args.output_root.exists():raise FileExistsError('Fresh output required; no overwrite or implicit resume')
    parent=args.output_root.resolve().parent
    while not parent.exists():parent=parent.parent
    if plan['output_reserve_bytes']>.9*shutil.disk_usage(parent).free:raise RuntimeError('Insufficient output volume')
    identity=dict(plan_sha256=sha(args.plan),operator_summary_sha256=sha(args.operator_root/'summary.json'),
        sources=sources,inputs=inputs,resources=resources,budget=budget,
        operator_root=str(args.operator_root.resolve()),physical_gap_dofs=2*(graph.n_nodes-len(graph.boundary_nodes)),
        spectral_roots=0,lu_factorizations=len(old['records']),time_contract='Full-node affine thermal system; baseline retained; same material/contacts/times',
        cost_estimate=plan['cost_estimate'])
    if not args.execute:print(json.dumps(dict(status='READY_DRY_RUN_NO_WRITES',**identity),indent=2));return
    args.output_root.mkdir(parents=True);write(args.output_root/'identity.json',identity);shutil.copyfile(args.plan,args.output_root/'executed_plan.json')
    log=(args.output_root/'progress.jsonl').open('x',buffering=1);started=time.monotonic();pool=None;affinity=os.sched_getaffinity(0)
    def event(name,**values):
        row=dict(event=name,elapsed_seconds=time.monotonic()-started,**values)
        line=json.dumps(row);print(line,flush=True);log.write(line+'\n')
    try:
        os.sched_setaffinity(0,{budget['coordinator_cpu']})
        pool=ResidentSpectrum(graph,old['records'],args.operator_root,plan['T_K']/plan['Tc_K'],budget,event)
        H,I=pool.apply(initial)
        herror=np.linalg.norm(H-Hexpected)/np.linalg.norm(Hexpected);ierror=np.linalg.norm(I-Iexpected)/np.linalg.norm(Iexpected)
        if max(herror,ierror)>plan['composition_relative_tolerance']:
            raise RuntimeError('Saved full-node operator composition differs from validated actions')
        event('COMPOSITION_CONFIRMED',hessian_relative_difference=float(herror),current_relative_difference=float(ierror))
        results={}
        for name,settings in plan['passes'].items():
            results[name]=run_pass(name,settings,plan,graph,d0,G0,I0,initial,old['base_free_energy'],pool,args.output_root,event)
        comparison=compare_passes(plan,graph,args.output_root,results['primary'],results['refined']);write(args.output_root/'refinement.json',comparison)
        if not comparison['all_comparisons_met']:raise RuntimeError('Planned temporal refinement does not support interpretation; no additional automatic retries')
        result=dict(status='FULL_NODE_AFFINE_THERMAL_TIME_COMPLETE',runtime_seconds=time.monotonic()-started,
            operator_batches=pool.actions,maximum_spectral_linear_solve_residual=pool.maximum_linear_solve_residual,
            horizon_ps=plan['observation_times_ps'][-1],refinement=comparison,production_changed=False,
            nonlinear_transient_admitted=False,photon=False,circuit_replaced=False,
            scope='Thermal linearized dynamics with affine baseline and exact P6 derivative; quadratic free energy and Taylor constitutive defect reported')
        write(args.output_root/'summary.json',result);event('COMPLETE',**result)
    except BaseException as exc:
        write(args.output_root/'failure.json',dict(type=type(exc).__name__,reason=str(exc),policy='Retain accepted checkpoints; no clipping, fallback, changed constants or automatic repeated campaign'))
        event('FAILED',reason=str(exc));raise
    finally:
        if pool is not None:pool.close()
        os.sched_setaffinity(0,affinity);log.close()


if __name__=='__main__':main()
