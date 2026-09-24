"""Full-node thermal nonlinear evolution; ETD2 and resident spectral workers.

The saved J0 is an exponential integrator's linear part, not a replacement
constitutive model. Every nonlinear RHS solves the same finite Matsubara
action at the current gap and keeps the inherited KWT/normal-potential law.
"""
from __future__ import annotations
import argparse,hashlib,json,multiprocessing as mp,os,shutil,sys,time,traceback
from pathlib import Path
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]
sys.path.insert(0,str(HERE))
from parallel_runtime import limit_thread_environment,linux_resources,resource_budget
limit_thread_environment();sys.path.insert(0,str(ROOT))
import numpy as np
from pysnspd.experimental import thermal_spatial_usadel as thermal
from pysnspd.experimental.thermal_stable_newton import solve_frequency
from pysnspd.experimental.thermal_snapshot import spectral_difference
from pysnspd.experimental.thermal_time_operator import SavedSpectralMode,FullNodeCoordinates
from pysnspd.experimental.thermal_weak_response import ThermalKWTNormal
from pysnspd.experimental.exponential_actions import etd2_trial


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def read(path):return json.loads(Path(path).read_text())
def write(path,value):
    with Path(path).open('x',encoding='utf8') as f:json.dump(value,f,indent=2,allow_nan=False);f.write('\n')


def actor(pipe,cpu,graph,rows,source,d0,plan):
    try:
        os.sched_setaffinity(0,{cpu});modes=[];t=plan['T_K']/plan['Tc_K']
        for row in rows:
            with np.load(Path(source)/row['fields_path']) as a:
                modes.append(SavedSpectralMode(graph,a,row['epsilon'],t))
        pipe.send(dict(event='READY',modes=len(modes),cpu=cpu))
        while True:
            command,values=pipe.recv()
            if command=='STOP':break
            if command not in ('LINEAR','NONLINEAR'):raise ValueError('Unknown worker request')
            values=np.asarray(values,complex);G=np.zeros_like(values)
            I=np.zeros((len(values),len(graph.edges)));E=np.zeros(len(values));maximum=0.;roots=0
            for mode in modes:
                for index,x in enumerate(values):
                    if command=='LINEAR':
                        g,current,residual=mode.apply(x);G[index]+=g;I[index]+=current
                    elif np.any(x):
                        n=graph.n_nodes;rhs=np.repeat(graph.area_weights,2)*np.column_stack((x.real,x.imag)).ravel()
                        du=np.zeros(2*n);du[mode.components]=mode.factor.solve(rhs[mode.components])
                        guess=mode.u+du[::2]+1j*du[1::2]
                        solution=solve_frequency(graph,d0+x,mode.epsilon,fixed_nodes=graph.boundary_nodes,
                            fixed_u=mode.u[graph.boundary_nodes],initial_u=guess,
                            tol=plan['spectral_tolerance'],max_iterations=plan['maximum_newton_iterations'])
                        diff=spectral_difference(graph,d0,mode.u,d0+x,solution.u,mode.epsilon)
                        G[index]+=mode.weight*diff['gap_gradient_difference']
                        I[index]+=mode.weight*diff['current_difference']
                        E[index]+=mode.weight*diff['renormalized_energy_difference']
                        residual=solution.residual;roots+=1
                    else:residual=0.
                    maximum=max(maximum,float(residual))
            pipe.send((G,I,E,maximum,roots))
    except BaseException as exc:
        try:pipe.send(dict(event='ERROR',reason=str(exc),traceback=traceback.format_exc()))
        except (BrokenPipeError,EOFError):pass
    finally:pipe.close()


class ResidentNonlinearSpectrum:
    def __init__(self,graph,rows,source,d0,plan,budget,event):
        self.graph=graph;self.d0=d0;self.logt=np.log(plan['T_K']/plan['Tc_K'])
        self.processes=[];self.pipes=[];self.linear_batches=0;self.nonlinear_batches=0
        self.roots=0;self.maximum_spectral_residual=0.;self.maximum_linear_residual=0.
        context=mp.get_context('spawn');ordered=sorted(rows,key=lambda row:row['n']);count=budget['workers']
        try:
            for index,cpu in enumerate(budget['worker_affinity_cpus']):
                parent,child=context.Pipe()
                process=context.Process(target=actor,args=(child,cpu,graph,ordered[index::count],str(source),d0,plan))
                process.start();child.close();self.processes.append(process);self.pipes.append(parent)
            for index,pipe in enumerate(self.pipes):
                message=pipe.recv()
                if message.get('event')!='READY':raise RuntimeError(str(message))
                event('WORKER_READY',worker=index,**{k:v for k,v in message.items() if k!='event'})
        except BaseException:self.close();raise

    def evaluate(self,command,values):
        values=np.asarray(values,complex)
        for pipe in self.pipes:pipe.send((command,values))
        m=self.graph.area_weights;G=2*m[None,:]*self.logt*values
        I=np.zeros((len(values),len(self.graph.edges)));E=np.zeros(len(values))
        if command=='NONLINEAR':E=np.sum(m[None,:]*(2*np.real(np.conj(self.d0)*values)+abs(values)**2),axis=1)*self.logt
        for pipe in self.pipes:
            result=pipe.recv()
            if isinstance(result,dict):raise RuntimeError(str(result))
            g,current,energy,residual,roots=result;G+=g;I+=current;E+=energy;self.roots+=roots
            if command=='NONLINEAR':self.maximum_spectral_residual=max(self.maximum_spectral_residual,residual)
            else:self.maximum_linear_residual=max(self.maximum_linear_residual,residual)
        if command=='LINEAR':self.linear_batches+=1
        else:self.nonlinear_batches+=1
        return G,I,E

    def close(self):
        for process,pipe in zip(self.processes,self.pipes):
            if process.is_alive():
                try:pipe.send(('STOP',None))
                except (BrokenPipeError,EOFError):pass
        for process in self.processes:
            process.join(timeout=5)
            if process.is_alive():process.terminate();process.join(timeout=5)
        for pipe in self.pipes:pipe.close()


def run_pass(label,settings,plan,graph,d0,G0,I0,initial,spectrum,folder,event,horizon):
    folder.mkdir();coordinates=FullNodeCoordinates(graph,plan['gap_reference_kBTc']);size=coordinates.size
    deltas=plan['perturbation_amplitude']*initial
    scales=np.array([max(np.linalg.norm(coordinates.encode(x)) for x in deltas)]+
        [np.linalg.norm(coordinates.encode(x)) for x in deltas])
    if np.any(scales<=0):raise ValueError('Nonzero response scales required')
    options={key:plan[key] for key in ('T_K','Tc_K','tau_ee_Tc_ps','tau_ep_Tc_ps')}
    model0=ThermalKWTNormal(graph,d0,**options)
    def encode(fields):return np.concatenate([coordinates.encode(x)/scale for x,scale in zip(fields,scales)])
    def decode(vector):return np.array([coordinates.decode(x*scale) for x,scale in zip(vector.reshape(3,size),scales)])
    def guard(fields):
        maximum=float(np.max(abs(fields))/plan['gap_reference_kBTc'])
        if maximum>plan['maximum_relative_displacement']:
            raise RuntimeError('Registered weak neighborhood exceeded; no clipping or extrapolated continuation')
        return maximum
    def linear(vector):
        xs=decode(vector);G,I,_=spectrum.evaluate('LINEAR',xs)
        return encode([model0.rhs_tangent(x,g,G0,I0,current)['velocity_direction'] for x,g,current in zip(xs,G,I)])
    def evaluate(vector):
        xs=decode(vector);guard(xs);G,I,E=spectrum.evaluate('NONLINEAR',xs)
        responses=[ThermalKWTNormal(graph,d0+x,**options).response(G0+g,I0+current) for x,g,current in zip(xs,G,I)]
        return encode([r['velocity'] for r in responses]),(xs,G,I,E,responses)
    def rhs(vector):return evaluate(vector)[0]
    names=['baseline','amplitude','angular_phase'];x=encode(np.array([np.zeros_like(d0),deltas[0],deltas[1]]))
    f,cache=evaluate(x);time_ps=0.;step_ps=settings['initial_step_ps'];started=time.monotonic()
    observations=[];history=[];attempt=0
    targets=[float(t) for t in plan['observation_times_ps'] if t<=horizon]
    if targets[-1]!=horizon:targets.append(float(horizon))
    for target in targets:
        while time_ps<target:
            attempt+=1
            if attempt>plan['maximum_step_attempts']:raise RuntimeError('ETD2 attempt budget exhausted; no automatic rerun')
            hps=min(step_ps,settings['maximum_step_ps'],target-time_ps);h=hps/model0.tD_ps
            candidate,detail=etd2_trial(linear,rhs,x,h,rhs_at_state=f,
                rtol=settings['krylov_rtol'],atol=settings['krylov_atol'],max_dimension=settings['maximum_krylov_dimension'])
            successful=detail['accepted'];indicator=None
            if successful:
                indicator=float(max(np.linalg.norm(v) for v in detail['correction'].reshape(3,size)))
                successful=bool(indicator<=settings['embedded_absolute_response_tolerance'])
            row=dict(attempt=attempt,start_ps=time_ps,step_ps=hps,accepted=successful,
                embedded_correction=indicator,embedded_tolerance=settings['embedded_absolute_response_tolerance'],
                first_dimension=detail['first'].dimension,first_defect=detail['first'].integrated_defect,
                second_dimension=None if detail['second'] is None else detail['second'].dimension,
                second_defect=None if detail['second'] is None else detail['second'].integrated_defect)
            if successful:
                f,cache=evaluate(candidate);x=candidate;time_ps=target if hps==target-time_ps else time_ps+hps
                path=folder/f'accepted_{attempt:05d}.npz';np.savez_compressed(path,time_ps=time_ps,displacements=cache[0])
                row.update(time_ps=time_ps,fields_path=path.name,fields_sha256=sha(path),
                    maximum_relative_displacement=guard(cache[0]))
                factor=2. if not indicator else min(2.,max(.5,.8*np.sqrt(settings['embedded_absolute_response_tolerance']/indicator)))
                step_ps=hps*factor
            else:
                step_ps=hps/2
                if step_ps<plan['minimum_step_ps']:raise RuntimeError('Minimum ETD2 step reached; no changed physics or fallback')
            history.append(row);fraction=time_ps/horizon;elapsed=time.monotonic()-started
            event('ETD2_STEP',pass_name=label,**row,nonlinear_batches=spectrum.nonlinear_batches,
                spectral_roots=spectrum.roots,eta_seconds=None if not fraction else elapsed*(1-fraction)/fraction)
            print(f"[{('#'*int(24*fraction)).ljust(24,'-')}] {label} {time_ps:.6g}/{horizon:g} ps | {elapsed:.1f}s | ETA {elapsed*(1-fraction)/fraction:.1f}s" if fraction else f'{label}: adjusting initial step',flush=True)
        xs,G,I,E,responses=cache;arrays=dict(d0=d0,G0=G0,area_weights=graph.area_weights,coordinates_bar=graph.coordinates_bar,
            edges=graph.edges,conductance=graph.conductance,boundary_nodes=graph.boundary_nodes)
        states=[]
        for name,delta,g,current,energy,r in zip(names,xs,G,I,E,responses):
            for key,value in dict(displacement=delta,gap_gradient_increment=g,current_increment=current,velocity=r['velocity'],potential=r['potential_v']).items():arrays[name+'_'+key]=value
            # Compatibility with the comparison helper: this is the exact
            # nonlinear force difference; NOT a Hessian action.
            arrays[name+'_hessian']=g
            states.append(dict(name=name,free_energy_change=float(energy),free_energy_rate=r['free_energy_rate'],
                kwt_loss=r['kwt_loss'],normal_loss=r['normal_loss'],power_residual=r['dissipation_residual'],
                continuity_max=r['continuity_max'],noether_max=r['noether_max']))
        arrays.update(amplitude_difference=xs[1]-xs[0],angular_phase_difference=xs[2]-xs[0])
        path=folder/f'observation_{len(observations):03d}.npz';np.savez_compressed(path,**arrays)
        observations.append(dict(time_ps=target,states=states,fields_path=path.relative_to(folder.parent).as_posix(),fields_sha256=sha(path)))
        event('OBSERVATION',pass_name=label,time_ps=target,count=len(observations),total=len(targets))
    result=dict(observations=observations,history=history,runtime_seconds=time.monotonic()-started,
        response_scales=scales.tolist(),auxiliary_constant_state=False)
    write(folder/'summary.json',result);return result


def compare_exact(plan,graph,output,first,second):
    """Use exact nonlinear torque instead of the old helper's tangent torque."""
    records=[];mass=graph.area_weights;active=graph.conductance>0
    if len(first['observations'])!=len(second['observations']):
        raise ValueError('Observation counts differ')
    for left,right in zip(first['observations'],second['observations']):
        if left['time_ps']!=right['time_ps']:raise ValueError('Time horizons differ')
        with np.load(output/left['fields_path']) as a,np.load(output/right['fields_path']) as b,np.load(output/second['observations'][0]['fields_path']) as initial:
            def field(data,probe,kind):
                if probe=='baseline':
                    if kind=='displacement':return data['baseline_displacement']
                    if kind=='current':return data['baseline_current_increment']
                    G=data['G0'];d=data['d0'];gb=data['baseline_gap_gradient_increment']
                    if kind=='force_density':return gb/mass
                    return np.imag(np.conj(d+data['baseline_displacement'])*(G+gb)-np.conj(d)*G)/mass
                if kind=='displacement':return data[probe+'_difference']
                if kind=='current':return data[probe+'_current_increment']-data['baseline_current_increment']
                G=data['G0'];d=data['d0'];gp=data[probe+'_gap_gradient_increment'];gb=data['baseline_gap_gradient_increment']
                if kind=='force_density':return (gp-gb)/mass
                return np.imag(np.conj(d+data[probe+'_displacement'])*(G+gp)-np.conj(d+data['baseline_displacement'])*(G+gb))/mass
            for probe in ('baseline','amplitude','angular_phase'):
                for kind in ('displacement','current','force_density','phase_torque_density'):
                    norm=(lambda z:float(np.sqrt(np.sum(abs(z[active])**2/graph.conductance[active])))) if kind=='current' else (lambda z:float(np.sqrt(np.sum(mass*abs(z)**2))))
                    error=norm(field(a,probe,kind)-field(b,probe,kind))
                    scale=(max(norm(field(initial,p,kind)) for p in ('amplitude','angular_phase'))
                        if probe=='baseline' else norm(field(initial,probe,kind)))
                    remaining=norm(field(b,probe,kind))
                    records.append(dict(time_ps=left['time_ps'],probe=probe,observable=kind,error_L2=error,
                        initial_scale_source='maximum initial perturbation response' if probe=='baseline' else 'initial perturbation response',
                        relative_to_initial=error/scale if scale>1e-14 else None,
                        relative_to_remaining=error/remaining if remaining>1e-14 else None,
                        admitted=bool(error<=plan['refinement_absolute_norm']+plan['refinement_relative_limit']*max(scale,remaining))))
    return dict(records=records,all_comparisons_met=all(row['admitted'] for row in records),
        meaning='One planned temporal refinement of baseline and both differences against initial response scales; tiny tails do not define relative accuracy')


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    for name in ('plan','operator-root','output-root'):parser.add_argument('--'+name,type=Path,required=True)
    parser.add_argument('--execute',action='store_true');parser.add_argument('--stop-after-ps',type=float)
    args=parser.parse_args();plan=read(args.plan);old=read(args.operator_root/'summary.json')
    for name,digest in plan['frozen_sources'].items():
        if sha(ROOT/name)!=digest:raise ValueError('Frozen numerical source changed: '+name)
    for name,digest in plan['operator_inputs'].items():
        if sha(args.operator_root/name)!=digest:raise ValueError('Admitted operator input changed: '+name)
    original=read(args.operator_root/'executed_plan.json')
    for key in ('T_K','Tc_K','tau_ee_Tc_ps','tau_ep_Tc_ps','matsubara_count'):
        if plan[key]!=original[key]:raise ValueError('Physical system changed: '+key)
    if sorted(row['n'] for row in old['records'])!=list(range(plan['matsubara_count'])):raise ValueError('Incomplete spectral sum')
    for row in old['records']:
        if sha(args.operator_root/row['fields_path'])!=row['fields_sha256']:raise ValueError('Spectral checkpoint changed')
    if args.output_root.exists():raise FileExistsError('Fresh output required; no overwrite or implicit resume')
    horizon=plan['observation_times_ps'][-1] if args.stop_after_ps is None else args.stop_after_ps
    if not np.isfinite(horizon) or not 0<horizon<=plan['observation_times_ps'][-1]:raise ValueError('Invalid observation horizon')
    with np.load(args.operator_root/'full_node_operator_checks.npz') as a:
        graph=thermal.ThermalGraph(a['area_weights'],a['edges'],a['conductance'],a['coordinates_bar'],a['boundary_nodes'])
        d0,G0,I0,initial=[a[k].copy() for k in ('d','gap_gradient','current','input_directions')]
    resources=linux_resources();budget=resource_budget(resources,.9,.9,max_workers=27)
    parent=args.output_root.resolve().parent
    while not parent.exists():parent=parent.parent
    if plan['output_reserve_bytes']>.9*shutil.disk_usage(parent).free:raise RuntimeError('Insufficient disk space')
    identity=dict(plan_sha256=sha(args.plan),sources={**plan['frozen_sources'],**{p:sha(ROOT/p) for p in plan['new_sources']}},
        operator_inputs=plan['operator_inputs'],resources=resources,budget=budget,operator_root=str(args.operator_root.resolve()),
        horizon_ps=horizon,pilot=args.stop_after_ps is not None,physical_constants_changed=False)
    if not args.execute:print(json.dumps(dict(status='DRY_RUN_NO_WRITES',**identity),indent=2));return
    args.output_root.mkdir(parents=True);write(args.output_root/'identity.json',identity);shutil.copyfile(args.plan,args.output_root/'executed_plan.json')
    log=(args.output_root/'progress.jsonl').open('x',buffering=1);started=time.monotonic();pool=None;affinity=os.sched_getaffinity(0)
    def event(name,**values):
        line=json.dumps(dict(event=name,elapsed_seconds=time.monotonic()-started,**values));print(line,flush=True);log.write(line+'\n')
    try:
        os.sched_setaffinity(0,{budget['coordinator_cpu']})
        pool=ResidentNonlinearSpectrum(graph,old['records'],args.operator_root,d0,plan,budget,event)
        results={name:run_pass(name,settings,plan,graph,d0,G0,I0,initial,pool,args.output_root/name,event,horizon) for name,settings in plan['passes'].items()}
        comparison=compare_exact(plan,graph,args.output_root,results['primary'],results['refined']);write(args.output_root/'refinement.json',comparison)
        result=dict(status='NONLINEAR_THERMAL_TIME_PILOT_COMPLETE' if args.stop_after_ps is not None else 'NONLINEAR_THERMAL_TIME_COMPLETE',
            runtime_seconds=time.monotonic()-started,refinement=comparison,horizon_ps=horizon,
            nonlinear_batches=pool.nonlinear_batches,linear_batches=pool.linear_batches,spectral_roots=pool.roots,
            maximum_spectral_residual=pool.maximum_spectral_residual,maximum_linear_residual=pool.maximum_linear_residual,
            stage4_complete=False,nonthermal_work_closed=False,production_changed=False,photon=False,circuit_replaced=False)
        if not comparison['all_comparisons_met']:result['status']='TEMPORAL_REFINEMENT_NOT_MET'
        write(args.output_root/'summary.json',result)
        if not comparison['all_comparisons_met']:
            raise RuntimeError('Planned temporal refinement was not met; retain both trajectories and revise explicitly')
        event('COMPLETE',**{k:v for k,v in result.items() if k!='refinement'})
    except BaseException as exc:
        write(args.output_root/'failure.json',dict(reason=str(exc),type=type(exc).__name__,policy='Retain evidence; no changed constants, fallback, clipping or automatic rerun'))
        event('FAILED',reason=str(exc));raise
    finally:
        if pool is not None:pool.close()
        os.sched_setaffinity(0,affinity);log.close()


if __name__=='__main__':main()
