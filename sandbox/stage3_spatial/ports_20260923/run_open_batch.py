"""Registered, static open-wire experiment. No time integration or 2D interface.

Default is a read-only description. --pilot or --execute creates a NEW output
directory. Failures are saved; no fallback, overwrite or retry is performed.
"""
from pathlib import Path
import argparse
from dataclasses import asdict
import hashlib
import json
import sys
import time
import traceback

import numpy as np
from scipy.optimize import minimize

ROOT=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(ROOT))
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from progress import Progress
from pysnspd.experimental.energy_catalog import OccupationEnergyCatalog
from pysnspd.experimental.refined_cells import refined_count_catalog
from pysnspd.experimental.spatial_open import OpenSpatialFunctional
from pysnspd.experimental.superconducting_reservoir import SuperconductingReservoir


def plain(value):
    if isinstance(value,np.ndarray):return value.tolist()
    if isinstance(value,np.generic):return value.item()
    if isinstance(value,dict):return {k:plain(v) for k,v in value.items()}
    if isinstance(value,(tuple,list)):return [plain(v) for v in value]
    return value


def save(path,value):
    Path(path).write_text(json.dumps(plain(value),indent=2,allow_nan=False)+'\n',encoding='utf-8')


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def sources(registration):
    names=['pysnspd/experimental/'+n+'.py' for n in
        ('energy_catalog','refined_cells','spatial_functional','spatial_open','superconducting_reservoir')]
    names+=['sandbox/stage3_spatial/progress.py',
            'sandbox/stage3_spatial/ports_20260923/run_open_batch.py',
            'docs/implementation/stage1_r2/catalogs/occupation_catalog.npz',
            str(Path(registration).resolve().relative_to(ROOT)).replace('\\','/')]
    return {n:sha(ROOT/n) for n in names}


def solve_case(catalog,reg,case,progress,output):
    start=time.monotonic();geometry=reg['geometry'];solver=reg['solver'];criteria=reg['criteria']
    model=OpenSpatialFunctional(catalog,case['length_m'],geometry['cross_section_m2'],
                                case['elements'],degree=reg['degree'])
    reservoir=SuperconductingReservoir(catalog,reg['bath_theta'],
        cross_section_m2=geometry['cross_section_m2'],resolved_length_m=case['length_m'])
    progress.checkpoint('resolviendo rama térmica del reservorio',force=True)
    branch=reservoir.solve_at_q(reg['q_bare_bar'],reg['amplitude_bracket'],
                               continuation_qs=reg['continuation_qs'])
    partition=reservoir.inductance_partition(branch)
    save(output/'reservoir.json',dict(branch=asdict(branch),inductance_partition=asdict(partition)))
    n=model.cells;a0=branch.amplitude_bar;j0=branch.current_bar
    # Fixing theta_R=0 removes the null global phase only. theta_L is free;
    # the Gibbs boundary work +J*theta_L imposes the left current flux.
    theta0=branch.q_bare_bar*(model.x_bar-model.length_bar)
    initial=np.r_[np.full(n-2,a0),theta0[:-1]]
    phase_radius=solver['phase_trust_radius_rad']
    bounds=[tuple(reg['amplitude_bracket'])]*(n-2)+[(v-phase_radius,v+phase_radius) for v in theta0[:-1]]
    evaluations=0;last_x=None;last=None;history=[]
    def state(y):
        return np.r_[a0,y[:n-2],a0]*np.exp(1j*np.r_[y[n-2:],0.])
    def objective(y):
        nonlocal evaluations,last_x,last
        if last_x is not None and np.array_equal(y,last_x):return last
        z=state(y)
        result=model.evaluate_thermal(z,reg['bath_theta'],require_stability=False,
            on_node=lambda i,total:progress.checkpoint(f"consulta {evaluations+1}, nodo {i+1}/{total}"))
        g=result.gradient_cartesian_bar[:,0]+1j*result.gradient_cartesian_bar[:,1]
        radial=np.real(np.conj(z)/abs(z)*g)
        phase=np.imag(np.conj(z)*g);phase[0]+=j0
        gradient=np.r_[radial[1:-1],phase[:-1]]
        value=result.free_energy_bar+j0*y[n-2]
        evaluations+=1;last_x=y.copy();last=(float(value),gradient)
        history.append(dict(evaluation=evaluations,seconds=time.monotonic()-start,
            gibbs_bar=value,max_gradient=float(np.max(abs(gradient)))))
        save(output/'iterations.json',history)
        return last
    result=minimize(objective,initial,jac=True,method='L-BFGS-B',bounds=bounds,
        options=dict(maxiter=solver['max_iterations'],maxfun=solver['max_evaluations'],
                     gtol=solver['gradient_tolerance'],ftol=solver['function_tolerance'],maxls=20))
    z=state(result.x)
    progress.checkpoint('comprobando D.36 en la solución abierta',force=True)
    final=model.evaluate_thermal(z,reg['bath_theta'],require_stability=True,
        on_node=lambda i,total:progress.checkpoint(f"D.36 nodo {i+1}/{total}"))
    _,grad=objective(result.x)
    current_relative=float(np.max(abs(final.link_current_bar-j0))/abs(j0))
    endpoint_q=final.q_delta_nodes_bar[[0,-1]]/(a0*a0/(a0*a0+.01))
    q_relative=float(np.max(abs(endpoint_q-branch.q_bare_bar))/abs(branch.q_bare_bar))
    amplitude_absolute=float(np.max(abs(abs(z)-a0)))
    active=any(min(v-lo,hi-v)<1e-6 for v,(lo,hi) in zip(result.x,bounds))
    # Residuals decide acceptance; optimizer termination alone is not evidence.
    accepted=(np.max(abs(grad))<=criteria['stationarity_absolute'] and
              current_relative<=criteria['current_relative'] and
              q_relative<=criteria['endpoint_q_relative'] and
              amplitude_absolute<=criteria['amplitude_absolute'] and not active)
    record=dict(status='PASS_STATIC_OPEN_CASE' if accepted else 'STATIC_OPEN_CASE_NOT_ACCEPTED',
        case=case,runtime_seconds=time.monotonic()-start,ell0_m=model.ell0_m,
        solver=dict(success=bool(result.success),message=str(result.message),iterations=int(result.nit),
                    evaluations=evaluations,active_bounds=active),
        metrics=dict(stationarity_absolute=float(np.max(abs(grad))),current_relative=current_relative,
                     endpoint_q_relative=q_relative,amplitude_absolute=amplitude_absolute,
                     minimum_principal_eigenvalue=min(s.eigenvalues[0] for s in final.principal_symbols),
                     maximum_symbol_uncertainty=max(s.uncertainty for s in final.principal_symbols),
                     global_phase_residual=final.global_phase_residual_bar),
        x_m=model.x_bar*model.ell0_m,amplitude_bar=abs(z),phase_rad=np.r_[result.x[n-2:],0.],
        current_A=final.current_A,q_delta_bar=final.q_delta_nodes_bar,
        free_energy_bar=final.free_energy_bar,internal_energy_bar=final.energy_bar,
        reference_current_A=branch.current_A,reference_amplitude_bar=a0,
        reference_q_bare_bar=branch.q_bare_bar,L_res_H=partition.resolved_differential_H,
        L_ext_H=partition.exterior_fixed_H,stage3_closed=False,production=False)
    save(output/'result.json',record)
    return record


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--registration',default='docs/implementation/stage3/ports_20260923/open_registration.json')
    parser.add_argument('--output-root',default='tmp/stage3b_open_full_20260923')
    group=parser.add_mutually_exclusive_group();group.add_argument('--execute',action='store_true');group.add_argument('--pilot',action='store_true')
    args=parser.parse_args();reg_path=ROOT/args.registration
    reg=json.loads(reg_path.read_text(encoding='utf-8'));case_list=[reg['pilot']] if args.pilot else reg['cases']
    if not (args.execute or args.pilot):
        print(json.dumps(dict(status='READY_NOT_EXECUTED',cases=case_list,criteria=reg['criteria'],sources=sources(reg_path)),indent=2));return
    out=ROOT/args.output_root
    if out.exists():raise FileExistsError('Output directory exists; choose a new path. No overwrite or automatic retry.')
    out.mkdir(parents=True)
    save(out/'manifest.json',dict(sources=sources(reg_path),registration=reg,arguments=vars(args),
        started_utc=__import__('datetime').datetime.now(__import__('datetime').timezone.utc).isoformat()))
    progress=Progress(len(case_list),'Etapa3B abierta',total_weight=sum(c['elements'] for c in case_list),
        callback=lambda event: (out/'progress.jsonl').open('a',encoding='utf-8').write(json.dumps(event)+'\n'))
    start=time.monotonic();records=[]
    try:
        catalog=refined_count_catalog(OccupationEnergyCatalog.load(ROOT/reg['catalogue']))
        for case in case_list:
            progress.start_task(case['id'],weight=case['elements'])
            dest=out/case['id'];dest.mkdir()
            record=solve_case(catalog,reg,case,progress,dest);records.append(record)
            if record['status']!='PASS_STATIC_OPEN_CASE':raise RuntimeError('Registered open-case criteria not met; no fallback or continuation.')
            progress.advance()
        summary=dict(status='PILOT_COMPLETE_SCOPE_LIMITED' if args.pilot else 'OPEN_STATIC_BATCH_COMPLETE_SCOPE_LIMITED',
            runtime_seconds=time.monotonic()-start,cases=records,stage3_closed=False,
            excluded=['2D-1D transparent interface','nonuniform driven boundary perturbation',
                      'coupled kinetic/condensate/circuit transient','production detector prediction'])
        save(out/'summary.json',summary);progress.finish(success=True)
    except BaseException as error:
        save(out/'failure.json',dict(exception=type(error).__name__,reason=str(error),
            runtime_seconds=time.monotonic()-start,completed_cases=len(records),traceback=traceback.format_exc()))
        progress.finish(success=False);raise


if __name__=='__main__':main()
