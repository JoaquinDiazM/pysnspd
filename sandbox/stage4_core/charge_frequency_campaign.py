"""Harmonic charge-response diagnostic using already computed spatial spectra.

No spectral solve, gap relaxation, physical collision, photon or circuit is
introduced. Without --execute validates the immutable inputs and resource budget.
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
import warnings

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
sys.path.insert(0,str(HERE))
from parallel_runtime import initialize_affinity,limit_thread_environment,linux_resources,resource_budget
limit_thread_environment()
sys.path.insert(0,str(ROOT))
import numpy as np
from scipy.constants import hbar,k
from scipy.sparse import diags
from scipy.sparse.linalg import MatrixRankWarning,spsolve
from pysnspd.experimental import thermal_spatial_usadel as thermal
from pysnspd.experimental import frozen_kinetic_usadel as kinetic
from pysnspd.experimental import frozen_charge_response as response


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def read(path):return json.loads(Path(path).read_text(encoding='utf8'))
def write(path,value):
    with Path(path).open('x',encoding='utf8') as stream:
        json.dump(value,stream,indent=2,allow_nan=False);stream.write('\n')


def worker(job):
    start=time.monotonic();plan=job['plan'];nu=job['nu'];rows=job['selected']
    with np.load(ROOT/job['case']['fields_path']) as a:
        d,xy,m=a['d'],a['coordinates_bar'],a['area_weights']
        graph=thermal.ThermalGraph(m,a['edges'],a['conductance'],xy,a['boundary_nodes'])
    free=np.ones(len(d),bool);free[graph.boundary_nodes]=False
    radius=np.linalg.norm(xy,axis=1);core=free&(radius<=plan['probe_radius_ell0'])
    edge_mid=(xy[graph.edges[:,0]]+xy[graph.edges[:,1]])/2
    edge_core=np.linalg.norm(edge_mid,axis=1)<=plan['probe_radius_ell0']
    node_norm=lambda value:float(np.sqrt(np.sum(m[core,None]*abs(value[core])**2))) if value.ndim==2 else float(np.sqrt(np.sum(m[core]*abs(value[core])**2)))
    flux_norm=lambda value:float(np.sqrt(np.sum(abs(value[edge_core])**2/graph.conductance[edge_core])))
    energies=plan['gap_reference_kBTc']*np.array([row['energy_relative'] for row in rows])
    weights=np.zeros(len(energies));weights[:-1]+=np.diff(energies)/2;weights[1:]+=np.diff(energies)/2
    exp_minus=np.exp(-energies/plan['temperature_ratio'])
    chi=2*exp_minus/(plan['temperature_ratio']*(1+exp_minus)**2)
    operators=[];projected_matrix=None;projected_mass=np.zeros(len(d));charge=np.arange(1,2*len(d),2)
    for weight,ch,row in zip(weights,chi,rows):
        with np.load(Path(job['spectra_root'])/row['fields_path']) as a:
            if not np.array_equal(d,a['d']):raise ValueError('Frozen gap changed')
            op=kinetic.assemble(graph,d,a['g'],a['f'],a['f_tilde'])
        operators.append(op)
        item=weight*ch*(-op.matrix[charge][:,charge])
        projected_matrix=item if projected_matrix is None else projected_matrix+item
        projected_mass+=weight*ch*m*op.R[:,0,0].real
    projected_matrix=projected_matrix-1j*nu*diags(projected_mass)
    fields=dict(coordinates_bar=xy,d=d,area_weights=m,edges=graph.edges,conductance=graph.conductance,
        boundary_nodes=graph.boundary_nodes,energies_kBTc=energies,weights_kBTc=weights,chi=chi)
    records={}
    for probe in plan['probes']:
        bump=np.maximum(0.,1-radius**2/plan['probe_radius_ell0']**2)**3
        spatial=bump if probe=='radial' else xy[:,0]/plan['probe_radius_ell0']*bump
        hLs=[row['energy_relative']*max(0.,1-row['energy_relative']**2/25)**3*spatial for row in rows]
        rhs=sum(w*(op.matrix[charge][:,charge-1]@hL) for w,op,hL in zip(weights,operators,hLs))
        v=np.zeros(len(d),complex)
        with warnings.catch_warnings():
            warnings.simplefilter('error',MatrixRankWarning)
            v[free]=spsolve(projected_matrix[free][:,free],rhs[free])
        if np.any(~np.isfinite(v)):raise RuntimeError('Nonfinite one-coordinate response')
        integrated={name:dict(force_cartesian=np.zeros((len(d),2),complex),current=np.zeros(len(graph.edges),complex),
            energy_flux=np.zeros(len(graph.edges),complex),charge_residual=np.zeros(len(d),complex),
            charge_storage=np.zeros(len(d),complex),energy_residual=np.zeros(len(d),complex))
            for name in ('zero_hT','full','projected')}
        per_energy=[]
        for weight,E,ch,op,hL,row in zip(weights,energies,chi,operators,hLs,rows):
            full=response.solve(op,hL,graph.boundary_nodes,nu)
            entry=dict(energy_relative=row['energy_relative'],
                full_equation_free_max=float(np.max(abs(full.equation_residual[free]))),
                maximum_hT=float(np.max(abs(full.directions[:,1]))))
            for name,hT in (('zero_hT',np.zeros(len(d))),('full',full.directions[:,1]),('projected',ch*v)):
                obs=response.observe(op,np.column_stack((hL,hT)),eta=row['z_real'])
                target=integrated[name]
                for key in ('force_cartesian','current','charge_residual'):
                    target[key]+=weight*obs[key]
                for key in ('energy_flux','energy_residual'):
                    target[key]+=weight*E*obs[key]
                target['charge_storage']+=weight*full.storage*hT
                if name=='full':
                    entry['eta_artificial_charge_leakage_density_norm']=node_norm(obs['artificial_eta_leakage'][:,1]/m)
            per_energy.append(entry)
        phase=np.divide(d,abs(d),out=np.zeros_like(d),where=abs(d)>1e-10*plan['gap_reference_kBTc'])
        norms={}
        for name,mapping in integrated.items():
            force=mapping['force_cartesian']
            mapping['amplitude_force_density']=(phase.real*force[:,0]+phase.imag*force[:,1])/m
            mapping['phase_torque_density']=(d.real*force[:,1]-d.imag*force[:,0])/m
            divergence=np.zeros(len(d),complex)
            np.add.at(divergence,graph.edges[:,0],mapping['current'])
            np.add.at(divergence,graph.edges[:,1],-mapping['current'])
            ward=divergence+m*mapping['phase_torque_density']-2*mapping['charge_residual']
            dynamic=mapping['charge_residual']+1j*nu*mapping['charge_storage']
            norms[name]=dict(current=flux_norm(mapping['current']),energy_flux=flux_norm(mapping['energy_flux']),
                force_density=node_norm(force/m[:,None]),amplitude_force_density=node_norm(mapping['amplitude_force_density']),
                phase_torque_density=node_norm(mapping['phase_torque_density']),
                ward_max=float(np.max(abs(ward))),integrated_dynamic_residual_free_max=float(np.max(abs(dynamic[free]))))
            for key,value in mapping.items():fields[probe+'_'+name+'_'+key]=value
        fields[probe+'_projected_coordinate']=v
        records[probe]=dict(norms=norms,per_energy=per_energy,
            projected_equation_free_max=float(np.max(abs((projected_matrix@v-rhs)[free]))))
    eta=rows[0]['eta_relative']*plan['gap_reference_kBTc']
    return dict(id=job['id'],case_id=job['case']['id'],eta_relative=rows[0]['eta_relative'],nu=nu,
        worker_seconds=time.monotonic()-start,energy_nodes=len(rows),probes=records,
        omega_rad_per_ps=nu/plan['tD_ps'],inverse_omega_ps=None if nu==0 else plan['tD_ps']/nu,
        period_ps=None if nu==0 else 2*np.pi*plan['tD_ps']/nu,
        hbar_omega_over_gap=2*nu/plan['gap_reference_kBTc'],hbar_omega_over_kBTb=2*nu/plan['temperature_ratio'],
        hbar_omega_over_eta=2*nu/eta,
        label='static' if nu==0 else 'slow anchor; not uniform AC validation' if nu<=.001 else 'exploratory leading-adiabatic kinetic closure'),fields


def compare(output,records):
    comparisons=[]
    for row in records:
        static=next(r for r in records if r['case_id']==row['case_id'] and r['eta_relative']==row['eta_relative'] and r['nu']==0)
        with np.load(output/row['fields_path']) as a,np.load(output/static['fields_path']) as b:
            m,xy,c=a['area_weights'],a['coordinates_bar'],a['conductance']
            core=np.linalg.norm(xy,axis=1)<=4
            edge_core=np.linalg.norm((xy[a['edges'][:,0]]+xy[a['edges'][:,1]])/2,axis=1)<=4
            for probe in ('radial','angular'):
                for key in ('current','energy_flux','force_cartesian','amplitude_force_density','phase_torque_density'):
                    value=a[probe+'_full_'+key];staticvalue=b[probe+'_full_'+key];potential=a[probe+'_projected_'+key]
                    def norm(v):
                        if key in ('current','energy_flux'):return float(np.sqrt(np.sum(abs(v[edge_core])**2/c[edge_core])))
                        if key=='force_cartesian':return float(np.sqrt(np.sum(abs(v[core])**2/m[core,None])))
                        return float(np.sqrt(np.sum(m[core]*abs(v[core])**2)))
                    ref=norm(staticvalue);dynamic=norm(value-staticvalue);proj=norm(potential-value)
                    comparisons.append(dict(case_id=row['case_id'],eta_relative=row['eta_relative'],nu=row['nu'],probe=probe,moment=key,
                        static_norm=ref,dynamic_change=dynamic,relative_dynamic_change=dynamic/ref if ref>1e-12 else None,
                        dynamic_full_norm=norm(value),potential_error=proj,potential_relative_error=proj/norm(value) if norm(value)>1e-12 else None))
    return comparisons


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--plan',type=Path,required=True);parser.add_argument('--spectra-root',type=Path,required=True)
    parser.add_argument('--output-root',type=Path,required=True);parser.add_argument('--execute',action='store_true')
    parser.add_argument('--pilot',action='store_true');args=parser.parse_args()
    plan=read(args.plan);spectral=read(args.spectra_root/'summary.json')
    if sha(args.spectra_root/'summary.json')!=plan['spectral_summary_sha256']:raise ValueError('Spectral summary changed')
    sources=dict(plan['frozen_sources'])
    for name in ('sandbox/stage4_core/charge_frequency_campaign.py','pysnspd/experimental/frozen_charge_response.py','tests/test_frozen_charge_response.py'):
        sources[name]=sha(ROOT/name)
    for name,digest in plan['frozen_sources'].items():
        if sha(ROOT/name)!=digest:raise ValueError('Frozen source changed: '+name)
    cases={case['id']:case for case in plan['cases']};jobs=[];inputs={}
    for case in cases.values():
        if sha(ROOT/case['fields_path'])!=case['fields_sha256']:raise ValueError('Core input changed')
        inputs[case['fields_path']]=case['fields_sha256']
    for group in plan['groups']:
        rows=sorted([r for r in spectral['records'] if r['case_id']==group['case_id'] and r['eta_relative']==group['eta_relative']],key=lambda r:r['energy_relative'])
        if [r['energy_relative'] for r in rows]!=group['energies_relative']:raise ValueError('Energy grid differs')
        for row in rows:
            if sha(args.spectra_root/row['fields_path'])!=row['fields_sha256']:raise ValueError('Frozen spectral map changed')
            inputs['spectra:'+row['fields_path']]=row['fields_sha256']
        for nu in plan['nu_values']:
            if args.pilot and not(group['case_id']=='radial_129_N256' and nu==.01):continue
            jobs.append(dict(id=group['id']+'_nu'+str(nu).replace('.','p'),case=cases[group['case_id']],plan=plan,
                selected=rows,nu=nu,spectra_root=str(args.spectra_root.resolve())))
    resources=linux_resources();budget=resource_budget(resources,.9,.9,max_workers=len(jobs))
    reserve=2*1024**3
    count=min(budget['workers'],(budget['memory_limit_bytes']-budget['coordinator_reserve_bytes'])//reserve)
    if count<1:raise RuntimeError('Not enough reserved memory')
    budget.update(workers=int(count),total_processes=int(count)+1,worker_affinity_cpus=budget['worker_affinity_cpus'][:count],
        worker_reserve_bytes=reserve,estimated_memory_reservation_bytes=budget['coordinator_reserve_bytes']+count*reserve)
    if args.output_root.exists():raise FileExistsError('Fresh output required; no overwrite or resume')
    parent=args.output_root.resolve().parent
    while not parent.exists():parent=parent.parent
    diskfree=shutil.disk_usage(parent).free
    if 512*1024**2>.9*diskfree:raise RuntimeError('Insufficient output volume reserve')
    identity=dict(plan_sha256=sha(args.plan),sources=sources,inputs=inputs,resources=resources,budget=budget,
        jobs=len(jobs),spectral_solves=0,pilot=args.pilot,output_reserve_bytes=512*1024**2,available_disk_bytes=diskfree)
    if not args.execute:print(json.dumps(dict(status='DRY_RUN',**identity),indent=2));return
    args.output_root.mkdir(parents=True);(args.output_root/'fields').mkdir();write(args.output_root/'identity.json',identity)
    shutil.copyfile(args.plan,args.output_root/'executed_plan.json')
    started=time.monotonic();records=[];affinity=os.sched_getaffinity(0)
    try:
        os.sched_setaffinity(0,{budget['coordinator_cpu']});context=mp.get_context('spawn');counter=context.Value('i',0)
        with ProcessPoolExecutor(max_workers=budget['workers'],mp_context=context,initializer=initialize_affinity,
            initargs=(budget['worker_affinity_cpus'],counter)) as pool:
            futures=[pool.submit(worker,job) for job in jobs]
            for future in as_completed(futures):
                row,fields=future.result();path=args.output_root/'fields'/(row['id']+'.npz')
                np.savez_compressed(path,**fields);row.update(fields_path=path.relative_to(args.output_root).as_posix(),fields_sha256=sha(path))
                write(path.with_suffix('.json'),row);records.append(row)
                elapsed=time.monotonic()-started;fraction=len(records)/len(jobs)
                print(f"[{'#'*int(24*fraction)}{'-'*(24-int(24*fraction))}] {len(records)}/{len(jobs)} | {elapsed:.1f} s | ETA {elapsed*(1-fraction)/fraction:.1f} s",flush=True)
        if not args.pilot:write(args.output_root/'comparisons.json',compare(args.output_root,records))
        write(args.output_root/'summary.json',dict(status='FIXED_SPECTRUM_ADIABATIC_FREQUENCY_DIAGNOSTIC_COMPLETE',runtime_seconds=time.monotonic()-started,
            records=records,spectral_solves=0,production_changed=False,detector_response_admitted=False,
            warning='A low-frequency kinetic closure diagnostic, not exact AC response, physical electrostatics or moving-gap dynamics'))
        print(json.dumps(dict(event='COMPLETE',runtime_seconds=time.monotonic()-started,jobs=len(jobs))),flush=True)
    except BaseException as exc:
        write(args.output_root/'failure.json',dict(type=type(exc).__name__,reason=str(exc),policy='No retry, fallback or overwrite'))
        raise
    finally:os.sched_setaffinity(0,affinity)


if __name__=='__main__':main()
