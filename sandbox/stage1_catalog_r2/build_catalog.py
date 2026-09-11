"""Build one bounded R2 candidate and measure field interpolation errors.

This is one stage of the frozen acceptance contract, not an acceptance decision.
The independent reference, quadrature, cutoff and eta checks must also pass.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import platform
import sys
import time

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT))
import numpy as np
from pysnspd.experimental.energy_catalog import (
    E_CHARGE_C, build_vacuum_catalog, build_occupation_catalog,
    segmented_count_quadrature, energy_at_count_batch, retarded_spectrum_batch,
    vacuum_state,
)

CRITERIA=ROOT/'docs/implementation/stage1_r2/acceptance_criteria.json'


def population(name,x):
    if name=='vacuum': return np.zeros_like(x)
    if name=='low_energy': return .2*np.exp(-x/.25)
    if name=='finite_energy_band': return .16*np.exp(-((x-.95)/.3)**2)
    if name=='near_edge_band': return .2*np.exp(-((x-.035)/.01)**2)
    if name=='high_energy_band': return .12*np.exp(-((x-3.2)/.55)**2)
    raise ValueError(name)


def direct(delta,gamma,x,w,p,eta):
    energy=energy_at_count_batch(x,delta=delta,gamma=gamma,eta=eta)
    c,s=retarded_spectrum_batch(energy,delta=delta,gamma=gamma,eta=eta)
    return np.asarray(vacuum_state(delta,gamma))+4*np.sum(
        np.array([energy,s.imag/c.real,-s.real*s.imag/c.real])*(p*w)[None,:],axis=1)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--delta-nodes',type=int,default=33)
    parser.add_argument('--gamma-low-nodes',type=int,default=49)
    parser.add_argument('--gamma-high-nodes',type=int,default=33)
    parser.add_argument('--count-order',type=int,default=10)
    parser.add_argument('--count-cutoff',type=float,default=12)
    parser.add_argument('--eta',type=float,default=1e-8)
    parser.add_argument('--ratio-coordinate',action='store_true')
    parser.add_argument('--output',type=Path,default=ROOT/'docs/implementation/stage1_r2/candidates/default')
    parser.add_argument('--catalog',type=Path)
    args=parser.parse_args()
    started=time.perf_counter()
    if min(args.delta_nodes,args.gamma_low_nodes,args.gamma_high_nodes)<4:
        raise ValueError('field axes require at least four nodes per segment')
    if args.count_cutoff<8: raise ValueError('candidate cutoff must cover the high-energy profile')
    criteria=json.loads(CRITERIA.read_text(encoding='utf-8'))
    args.output.mkdir(parents=True,exist_ok=True)
    deltas=np.geomspace(.08,1.5,args.delta_nodes)
    gamma_min=.03*args.eta**1.5/np.sqrt(1.5)
    gammas=np.r_[0,np.geomspace(gamma_min,.03,args.gamma_low_nodes),
                 np.linspace(.03,1.2,args.gamma_high_nodes)[1:]]
    ratio_axis=None
    coordinate_scale=args.eta**1.5/np.sqrt(1.5)
    if args.ratio_coordinate:
        coordinate_scale=(args.eta/1.5)**1.5
        distances=np.geomspace(1e-8,.1,args.gamma_high_nodes)
        ratio_axis=np.unique(np.r_[0,np.geomspace(.03*coordinate_scale,.03,args.gamma_low_nodes),
                                  np.linspace(.03,.9,2*args.gamma_high_nodes-1)[1:],
                                  (1-distances[::-1])[1:],1,1+distances,
                                  np.geomspace(1.1,15,args.gamma_high_nodes)[1:]])
        gammas=np.linspace(0,1.2,4)  # physical support of the exact vacuum
    breaks=np.array([0,1e-6,1e-5,3e-5,1e-4,3e-4,1e-3,.003,.01,.02,.05,.1,.25,.5,1,2,4,8])
    if args.count_cutoff>8: breaks=np.r_[breaks,args.count_cutoff]
    count,weights=segmented_count_quadrature(breaks,args.count_order)
    D,sigma=1.581e-4,4.2e5
    vacuum=build_vacuum_catalog(deltas,gammas,Tc_K=8.65,N0_per_J_m3=sigma/(2*E_CHARGE_C**2*D),
                                D_m2_s=D,analytic=True,reference_hashes={
                                    CRITERIA.relative_to(ROOT).as_posix():hashlib.sha256(CRITERIA.read_bytes()).hexdigest()})
    vacuum.metadata.update(iteration='stage1_r2',admission_status='PENDING_ALL_ACCEPTANCE_GATES',
                           quadrature_description='common segmented Gauss-Legendre nodes fixed across all fields',
                           count_breaks=breaks.tolist(),count_order_per_segment=args.count_order,
                           inversion_method='causal batched quartic plus safeguarded vector Newton in state count')
    table=build_occupation_catalog(vacuum,eta=args.eta,count_nodes=count,count_weights=weights,
                                   inverse_method='batch',compensated_gamma=True,direct_energy=True,
                                   gamma_coordinate_scale=coordinate_scale,gamma_ratio_axis=ratio_axis,
                                   gamma_offset_threshold=1e-10,gamma_offset_order=8)
    build_seconds=time.perf_counter()-started
    print(f'Built {table.excitation_energies.shape} in {build_seconds:.3f}s',flush=True)
    catalog_path=args.catalog or args.output/'occupation_catalog.npz'
    table.save(catalog_path)
    vacuum.save(args.output/'vacuum_catalog.npz')
    fields=[(f'fixed_{i:02d}',*point) for i,point in enumerate(criteria['mandatory_field_points'])]
    rng=np.random.default_rng(7431)
    for i in range(12):
        d=float(rng.uniform(.085,1.49))
        g=float(10**rng.uniform(-12,-1.5) if i<6 else rng.uniform(.031,1.19))
        fields.append((f'offnode_{i:02d}',d,g))
    rows=[]
    rejected_fields=[]
    for field_id,delta,gamma in fields:
        try:
            table.energy_kernel(delta,gamma)
        except FloatingPointError as error:
            rejected_fields.append({'field_id':field_id,'delta':delta,'gamma':gamma,'reason':str(error)})
            continue
        for profile in criteria['mandatory_population_profiles']:
            p=population(profile['id'],count)
            measured=np.asarray(table.evaluate(delta,gamma,p))
            reference=direct(delta,gamma,count,weights,p,args.eta)
            difference=np.abs(measured-reference)
            rows.append({'case_id':field_id+'_'+profile['id'],'field_id':field_id,
                         'delta':delta,'gamma':gamma,'profile':profile['id'],
                         'candidate':measured.tolist(),'same_eta_same_quadrature_reference':reference.tolist(),
                         'scaled_error':(difference/np.maximum(1,np.abs(reference))).tolist(),
                         'relative_error':[float(difference[k]/abs(reference[k])) if reference[k]!=0 else None for k in range(3)]})
    worst=[]
    for k,name in enumerate(('energy','amplitude_force','gamma_response')):
        maximum=max(rows,key=lambda row:row['scaled_error'][k])
        nonzero=[r for r in rows if r['relative_error'][k] is not None]
        rel=max(nonzero,key=lambda row:row['relative_error'][k])
        worst.append({'observable':name,'max_scaled_error':maximum['scaled_error'][k],
                      'max_scaled_case':maximum['case_id'],'max_relative_error':rel['relative_error'][k],
                      'max_relative_case':rel['case_id']})
    gates=criteria['numerical_gates']
    failures=[dict(row,failed_gates=['positive_monotone_energy']) for row in rejected_fields]
    for row in rows:
        limits=[gates['energy_max_scaled_error'],gates['amplitude_force_max_scaled_error'],gates['gamma_response_max_scaled_error']]
        failed=[name for i,name in enumerate(('energy','amplitude_force','gamma_response')) if row['scaled_error'][i]>limits[i]]
        failed += [name+'_relative' for i,name in ((1,'amplitude_force'),(2,'gamma_response'))
                   if row['relative_error'][i] is not None and row['relative_error'][i]>gates['nonzero_force_and_gamma_max_relative_error']]
        if failed: failures.append({'case_id':row['case_id'],'failed_gates':failed})
    report={'schema':'pysnspd.stage1_r2.candidate.v1','host':platform.node(),'python':platform.python_version(),
            'numpy':np.__version__,'criteria_id':criteria['criteria_id'],
            'criteria_sha256':hashlib.sha256(CRITERIA.read_bytes()).hexdigest(),
            'source_sha256':hashlib.sha256((ROOT/'pysnspd/experimental/energy_catalog.py').read_bytes()).hexdigest(),
            'runner_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            'catalog_sha256':hashlib.sha256(catalog_path.read_bytes()).hexdigest(),
            'catalog_bytes':catalog_path.stat().st_size,'build_seconds':build_seconds,
            'runtime_seconds':time.perf_counter()-started,
            'eta':args.eta,'shape':list(table.excitation_energies.shape),
            'delta_axis':deltas.tolist(),'gamma_axis':gammas.tolist(),
            'gamma_ratio_axis':None if ratio_axis is None else ratio_axis.tolist(),
            'count_breaks':breaks.tolist(),'count_order_per_segment':args.count_order,
            'vacuum_evaluation':'closed_form','worst_errors':worst,'rows':rows,'failures':failures,
            'rejected_fields':rejected_fields,
            'interpolation_only_status':'PASS_SAMPLED' if not failures else 'FAIL',
            'overall_acceptance_status':'PENDING_INDEPENDENT_REFERENCES_QUADRATURE_CUTOFF_REGULATOR_AND_API',
            'production_connected':False}
    (args.output/'candidate_diagnostics.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    print(json.dumps({k:report[k] for k in ('shape','build_seconds','runtime_seconds','worst_errors','interpolation_only_status')},indent=2))


if __name__=='__main__':
    main()
