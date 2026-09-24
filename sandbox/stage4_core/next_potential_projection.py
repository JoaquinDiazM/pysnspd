"""Coarse static projection hT(E,i)=chi(E)*v_i on frozen Usadel spectra.

The integrated charge equation determines a dimensionless electrochemical-shift
coordinate v, not an asserted memory electrostatic potential or time law. The
12-point energy grid is a diagnostic: its failure to integrate chi is retained.
No spectrum, material, collision term or population state is fitted here.
"""
from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
import hashlib
import json
import multiprocessing as mp
import os
from pathlib import Path
import sys
import time
import warnings

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(HERE))
from parallel_runtime import initialize_affinity, limit_thread_environment, linux_resources, resource_budget
limit_thread_environment()
sys.path.insert(0, str(ROOT))
import numpy as np
from scipy.sparse.linalg import MatrixRankWarning, spsolve
from pysnspd.experimental import thermal_spatial_usadel as thermal
from pysnspd.experimental import frozen_kinetic_usadel as kinetic

DATA = ROOT/'docs/implementation/stage4/self_consistent_review_20260924'


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read(path):
    return json.loads(Path(path).read_text(encoding='utf8'))


def write(path, value):
    with Path(path).open('x', encoding='utf8') as stream:
        json.dump(value, stream, indent=2, allow_nan=False); stream.write('\n')


def worker(job):
    started = time.monotonic()
    case, plan, selected = job['case'],job['plan'],job['selected']
    with np.load(ROOT/case['fields_path']) as a:
        d,xy,mass = a['d'],a['coordinates_bar'],a['area_weights']
        graph = thermal.ThermalGraph(mass,a['edges'],a['conductance'],xy,a['boundary_nodes'])
    free=np.ones(len(d),bool);free[graph.boundary_nodes]=False
    core=free & (np.linalg.norm(xy,axis=1)<=plan['probe_radius_ell0'])
    edge_midpoint=(xy[graph.edges[:,0]]+xy[graph.edges[:,1]])/2
    core_edges=(graph.conductance>0)&(np.linalg.norm(edge_midpoint,axis=1)<=plan['probe_radius_ell0'])
    node_norm=lambda v:float(np.sqrt(np.sum(mass[core]*abs(v[core])**2)))
    flux_norm=lambda v:float(np.sqrt(np.sum(abs(v[core_edges])**2/graph.conductance[core_edges])))
    energy=plan['gap_reference_kBTc']*np.array([r['energy_relative'] for r in selected])
    weights=np.zeros(len(energy));weights[:-1]+=np.diff(energy)/2;weights[1:]+=np.diff(energy)/2
    temperature=job['temperature_ratio']
    # Stable derivative of tanh(E/2T). Do not rescale to correct the quadrature.
    exp_minus=np.exp(-energy/temperature)
    chi=2*exp_minus/(temperature*(1+exp_minus)**2)
    operators=[];full_responses=[]
    charge=np.arange(1,2*len(d),2);longitudinal=charge-1
    projected_matrix=None
    for weight,susceptibility,row in zip(weights,chi,selected):
        with np.load(Path(job['spectra_root'])/row['fields_path']) as a:
            if not np.array_equal(d,a['d']):
                raise ValueError('Frozen gap differs')
            op=kinetic.assemble(graph,d,a['g'],a['f'],a['f_tilde'])
        operators.append(op)
        contribution=weight*susceptibility*op.matrix[charge][:,charge]
        projected_matrix=contribution if projected_matrix is None else projected_matrix+contribution
        response=job['kinetic_records'][row['id']]
        with np.load(Path(job['kinetic_root'])/response['fields_path']) as a:
            full_responses.append({name:a[name].copy() for name in a.files})
    fields=dict(coordinates_bar=xy,d=d,edges=graph.edges,area_weights=mass,conductance=graph.conductance,
        boundary_nodes=graph.boundary_nodes,energy_nodes_kBTc=energy,energy_weights_kBTc=weights,chi=chi)
    records={}
    for kind in plan['probes']:
        radius=np.linalg.norm(xy,axis=1)
        bump=np.maximum(0.,1-radius**2/plan['probe_radius_ell0']**2)**3
        spatial=bump if kind=='radial' else xy[:,0]/plan['probe_radius_ell0']*bump
        hLs=[row['energy_relative']*max(0.,1-row['energy_relative']**2/25)**3*spatial for row in selected]
        rhs=-sum(w*(op.matrix[charge][:,longitudinal]@hL) for w,op,hL in zip(weights,operators,hLs))
        v=np.zeros(len(d))
        with warnings.catch_warnings():
            warnings.simplefilter('error',MatrixRankWarning)
            v[free]=spsolve(projected_matrix[free][:,free],rhs[free])
        if np.any(~np.isfinite(v)):
            raise RuntimeError('Nonfinite projected potential, no regularizer added')
        integrated={name:dict(force=np.zeros(len(d),complex),current=np.zeros(len(graph.edges)),
            energy_flux=np.zeros(len(graph.edges)),charge_residual=np.zeros(len(d)),energy_residual=np.zeros(len(d)))
            for name in ('zero_hT','full','projected')}
        per_energy=[]
        for n,(w,E,ch,op,hL,row) in enumerate(zip(weights,energy,chi,operators,hLs,selected)):
            saved=full_responses[n]
            if not np.array_equal(hL,saved[kind+'_hL']):
                raise ValueError('Energy-mode probe differs from frozen response')
            states={'zero_hT':np.zeros(len(d)),'full':saved[kind+'_hT'],'projected':ch*v}
            per_row=dict(energy_kBTc=float(E),chi=float(ch))
            for name,hT in states.items():
                obs=op.evaluate(np.column_stack((hL,hT)),eta=row['z_real'])
                if name=='full':
                    if np.max(abs(obs.gap_force_increment-saved[kind+'_delta_gap_force']))>1e-11:
                        raise ValueError('Saved and recomputed full charge responses differ')
                target=integrated[name]
                target['force']+=w*obs.gap_force_increment
                target['current']+=w*obs.charge_current_increment
                target['energy_flux']+=w*E*obs.edge_flux[:,0]
                target['charge_residual']+=w*obs.residual[:,1]
                target['energy_residual']+=w*E*obs.residual[:,0]
                per_row[name+'_charge_residual_density_L2']=node_norm(obs.residual[:,1]/mass)
                per_row[name+'_hT_core_L2']=node_norm(hT)
            per_energy.append(per_row)
        comparisons={}
        phase=np.divide(d,abs(d),out=np.zeros_like(d),where=abs(d)>1e-10*plan['gap_reference_kBTc'])
        for name,mapping in integrated.items():
            divergence=np.zeros(len(d))
            np.add.at(divergence,graph.edges[:,0],mapping['current']);np.add.at(divergence,graph.edges[:,1],-mapping['current'])
            ward=divergence+np.imag(np.conj(d)*mapping['force'])-2*mapping['charge_residual']
            mapping['ward_residual']=ward
            mapping['amplitude_force_density']=np.real(np.conj(phase)*mapping['force'])/mass
            mapping['phase_torque_density']=np.imag(np.conj(d)*mapping['force'])/mass
            for field,value in mapping.items():fields[kind+'_'+name+'_'+field]=value
            comparisons[name]=dict(charge_current_norm=flux_norm(mapping['current']),
                energy_weighted_energy_flux_norm=flux_norm(mapping['energy_flux']),
                gap_force_density_norm=node_norm(mapping['force']/mass),
                amplitude_force_density_norm=node_norm(mapping['amplitude_force_density']),
                phase_torque_density_norm=node_norm(mapping['phase_torque_density']),
                integrated_charge_residual_density_L2=node_norm(mapping['charge_residual']/mass),
                maximum_free_integrated_charge_residual=float(np.max(abs(mapping['charge_residual'][free]))),
                maximum_joint_ward_residual=float(np.max(abs(ward))))
        reference=integrated['full']
        errors={}
        for label,key,norm in (('charge_current','current',flux_norm),('energy_weighted_energy_flux','energy_flux',flux_norm),
            ('gap_force_density','force',lambda value:node_norm(value/mass)),
            ('amplitude_force_density','amplitude_force_density',node_norm),('phase_torque_density','phase_torque_density',node_norm)):
            denominator=norm(reference[key]);difference=norm(integrated['projected'][key]-reference[key])
            errors[label]=dict(full_norm=denominator,projected_norm=norm(integrated['projected'][key]),
                difference_norm=difference,relative_difference=(difference/denominator if denominator>1e-12 else None))
        fields[kind+'_v']=v
        records[kind]=dict(v_maximum_absolute=float(np.max(abs(v))),v_core_L2=node_norm(v),
            norms=comparisons,projection_errors=errors,per_energy=per_energy,
            projected_equation_maximum_free_residual=float(np.max(abs((projected_matrix@v-rhs)[free]))),
            projection_enforces='Integrated charge balance only; not zero residual at every energy',
            v_convention='Dimensionless shift hdiag approximately h0(E+v),h0(E-v); sign mapping to memory phi not asserted')
    return dict(case_id=case['id'],eta_relative=selected[0]['eta_relative'],worker_seconds=time.monotonic()-started,
        temperature_ratio=temperature,chi_quadrature=float(np.dot(weights,chi)),
        chi_integral_exact_window=float(np.tanh(energy[-1]/(2*temperature))-np.tanh(energy[0]/(2*temperature))),
        chi_not_renormalized=True,probes=records),fields


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--plan',type=Path,default=DATA/'next_kinetic_plan.json')
    parser.add_argument('--spectra-root',type=Path,required=True)
    parser.add_argument('--kinetic-root',type=Path,required=True)
    parser.add_argument('--output-root',type=Path,required=True)
    parser.add_argument('--execute',action='store_true')
    args=parser.parse_args()
    plan=read(args.plan);spectral=read(args.spectra_root/'summary.json');responses=read(args.kinetic_root/'summary.json')
    if sha(args.spectra_root/'summary.json')!=plan['spectral_summary_sha256']:
        raise ValueError('Spectral summary differs from frozen plan')
    kinetic_identity=read(args.kinetic_root/'identity.json')
    if kinetic_identity['plan_sha256']!=sha(args.plan):
        raise ValueError('Kinetic plan changed')
    for name,expected in kinetic_identity['sources'].items():
        if sha(ROOT/name)!=expected:raise ValueError('Frozen kinetic source differs: '+name)
    inputs={}
    for case in plan['cases']:
        if sha(ROOT/case['fields_path'])!=case['fields_sha256']:raise ValueError('Core map changed')
        inputs[case['fields_path']]=case['fields_sha256']
    for prefix,base,rows in (('spectral',args.spectra_root,spectral['records']),('kinetic',args.kinetic_root,responses['records'])):
        for row in rows:
            if sha(base/row['fields_path'])!=row['fields_sha256']:raise ValueError('Frozen '+prefix+' map changed')
            inputs[prefix+':'+row['fields_path']]=row['fields_sha256']
    thermal_plan=read(ROOT/'docs/implementation/stage4/spatial_energy_20260924/self_consistent_plan.json')
    temperature=thermal_plan['T_K']/thermal_plan['Tc_K']
    lookup={r['id']:r for r in responses['records']}
    jobs=[]
    for case in plan['cases']:
        for eta in sorted({r['eta_relative'] for r in spectral['records']}):
            rows=sorted([r for r in spectral['records'] if r['case_id']==case['id'] and r['eta_relative']==eta],key=lambda r:r['energy_relative'])
            if [r['energy_relative'] for r in rows]!=plan['energies_relative']:raise ValueError('Incomplete energy grid')
            jobs.append(dict(case=case,plan=plan,selected=rows,kinetic_records=lookup,
                spectra_root=str(args.spectra_root),kinetic_root=str(args.kinetic_root),temperature_ratio=temperature))
    resources=linux_resources();budget=resource_budget(resources,.9,.9,max_workers=len(jobs))
    sources=[Path(__file__),HERE/'parallel_runtime.py',ROOT/'pysnspd/experimental/frozen_kinetic_usadel.py',ROOT/'pysnspd/experimental/thermal_spatial_usadel.py']
    identity=dict(plan_sha256=sha(args.plan),thermal_plan_sha256=sha(ROOT/'docs/implementation/stage4/spatial_energy_20260924/self_consistent_plan.json'),
        sources={path.relative_to(ROOT).as_posix():sha(path) for path in sources},inputs=inputs,
        spectral_summary_sha256=sha(args.spectra_root/'summary.json'),kinetic_summary_sha256=sha(args.kinetic_root/'summary.json'),
        resources=resources,budget=budget,jobs=len(jobs),new_spectral_solves=0,physical_time_steps=0)
    if not args.execute:print(json.dumps(dict(status='DRY_RUN',**identity),indent=2));return
    if args.output_root.exists():raise FileExistsError('Choose fresh output, no overwrite')
    args.output_root.mkdir(parents=True);(args.output_root/'fields').mkdir();write(args.output_root/'identity.json',identity)
    started=time.monotonic();records=[];affinity=os.sched_getaffinity(0)
    try:
        os.sched_setaffinity(0,{budget['coordinator_cpu']})
        context=mp.get_context('spawn');counter=context.Value('i',0)
        with ProcessPoolExecutor(max_workers=budget['workers'],mp_context=context,initializer=initialize_affinity,
            initargs=(budget['worker_affinity_cpus'],counter)) as pool:
            for future in as_completed([pool.submit(worker,job) for job in jobs]):
                record,fields=future.result();key=f"{record['case_id']}_eta{record['eta_relative']:.4f}".replace('.','p')
                path=args.output_root/'fields'/(key+'.npz');np.savez_compressed(path,**fields)
                record.update(id=key,fields_path=path.relative_to(args.output_root).as_posix(),fields_sha256=sha(path));records.append(record)
                elapsed=time.monotonic()-started
                print(json.dumps(dict(event='PROGRESS',complete=len(records),total=len(jobs),elapsed_seconds=elapsed,
                    eta_seconds=elapsed/len(records)*(len(jobs)-len(records)))),flush=True)
        summary=dict(status='COARSE_POTENTIAL_PROJECTION_DIAGNOSTIC',runtime_seconds=time.monotonic()-started,records=records,
            physical_time_steps=0,new_spectral_solves=0,chi_renormalized=False,production_changed=False,
            stage4_complete=False,exact_memory_ohmic_law_admitted=False,instantaneous_charge_dynamics_admitted=False,
            warning='Coarse 12-node integration poorly resolves thermal chi; compare static moment errors, do not certify a physical potential closure')
        write(args.output_root/'summary.json',summary)
        print(json.dumps(dict(event='COMPLETE',runtime_seconds=summary['runtime_seconds'],cases=len(records))),flush=True)
    finally:os.sched_setaffinity(0,affinity)


if __name__=='__main__':main()
