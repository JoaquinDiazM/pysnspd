"""Separate count quadrature, cutoff and regulator errors without a field table.

The highest two Gauss orders provide a resolved finite-eta reference check;
the causal-limit values come from the independent review parametrization.
Every fixed field/profile remains in the output, including any solver failure.
"""
from pathlib import Path
import argparse
import hashlib
import json
import sys
import time

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT))
import numpy as np
from scipy.special import expit, xlogy
from pysnspd.experimental.energy_catalog import (
    energy_at_count_batch,retarded_spectrum_batch,segmented_count_quadrature,vacuum_state,
)
from build_catalog import population


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,default=ROOT/'docs/implementation/stage1_r2/quadrature_diagnostics.json')
    args=parser.parse_args()
    started=time.perf_counter()
    criteria_path=ROOT/'docs/implementation/stage1_r2/acceptance_criteria.json'
    review=ROOT/'docs/implementation/stage1_r2/review'
    criteria=json.loads(criteria_path.read_text())
    causal=json.loads((review/'causal_reference.json').read_text())
    thermal=json.loads((review/'thermal_reference.json').read_text())
    references={}
    for row in causal['rows']:
        key=(row['amplitude'],row['gamma'],row['profile'])
        references[key]=row['reference']
    for row in thermal['rows']:
        key=(row['amplitude'],row['gamma'],'FD_'+str(row['temperature']))
        references[key]=row['causal_reference']
    fields={tuple(field):[] for field in criteria['mandatory_field_points']}
    for field in fields:
        fields[field]=[profile['id'] for profile in criteria['mandatory_population_profiles']]
    for row in thermal['rows']:
        fields.setdefault((row['amplitude'],row['gamma']),[]).append('FD_'+str(row['temperature']))
    base_breaks=np.array([0,1e-6,1e-5,3e-5,1e-4,3e-4,1e-3,.003,.01,.02,.05,.1,.25,.5,1,2,4,6,8,12,20])
    configurations=[(eta,order,12.) for eta in (1e-7,3e-8,1e-8) for order in (8,12,16,24,32)]
    configurations += [(1e-8,32,cutoff) for cutoff in (6.,8.,20.)]
    rows=[]; failures=[]; timings=[]
    for eta,order,cutoff in configurations:
        config_started=time.perf_counter()
        x,w=segmented_count_quadrature(base_breaks[base_breaks<=cutoff],order)
        for (amplitude,gamma),profiles in fields.items():
            try:
                energy=energy_at_count_batch(x,delta=amplitude,gamma=gamma,eta=eta)
                c,s=retarded_spectrum_batch(energy,delta=amplitude,gamma=gamma,eta=eta)
                kernels=np.array([energy,s.imag/c.real,-s.real*s.imag/c.real])
            except FloatingPointError as error:
                failures.append(dict(eta=eta,order=order,cutoff=cutoff,amplitude=amplitude,gamma=gamma,reason=str(error)))
                continue
            for profile in profiles:
                if profile.startswith('FD_'):
                    temperature=float(profile[3:])
                    p=expit(-energy/temperature)
                else:
                    temperature=0.
                    p=population(profile,x)
                values=np.asarray(vacuum_state(amplitude,gamma))+4*np.sum(kernels*(p*w)[None,:],axis=1)
                if temperature:
                    values[0]+=4*temperature*np.dot(w,xlogy(p,p)+xlogy(1-p,1-p))
                reference=np.asarray(references[(amplitude,gamma,profile)])
                error=abs(values-reference)
                rows.append(dict(eta=eta,order=order,cutoff=cutoff,count_nodes=len(x),amplitude=amplitude,
                                 gamma=gamma,profile=profile,values=values.tolist(),causal_reference=reference.tolist(),
                                 causal_scaled_error=(error/np.maximum(1,abs(reference))).tolist(),
                                 causal_relative_error=[float(error[k]/abs(reference[k])) if reference[k]!=0 else None for k in range(3)]))
        timings.append(dict(eta=eta,order=order,cutoff=cutoff,seconds=time.perf_counter()-config_started))
        print(json.dumps(timings[-1]),flush=True)
    index={(r['eta'],r['order'],r['cutoff'],r['amplitude'],r['gamma'],r['profile']):r for r in rows}
    for row in rows:
        ref=index.get((row['eta'],32,12.,row['amplitude'],row['gamma'],row['profile']))
        if ref is not None:
            error=abs(np.asarray(row['values'])-ref['values'])
            row['same_eta_order32_cutoff12_scaled_difference']=(error/np.maximum(1,abs(np.asarray(ref['values'])))).tolist()
            row['same_eta_order32_cutoff12_relative_difference']=[float(error[k]/abs(ref['values'][k])) if ref['values'][k]!=0 else None for k in range(3)]
    result=dict(schema='pysnspd.stage1_r2.separated_count_diagnostics.v1',
                criteria_sha256=hashlib.sha256(criteria_path.read_bytes()).hexdigest(),
                source_sha256=hashlib.sha256((ROOT/'pysnspd/experimental/energy_catalog.py').read_bytes()).hexdigest(),
                runner_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                reference_sha256={name:hashlib.sha256((review/name).read_bytes()).hexdigest() for name in ('causal_reference.json','thermal_reference.json')},
                interpretation='Direct causal spectral quadrature only. No field-interpolation acceptance is inferred. Order24-to32 checks finite-eta reference resolution; final table requires its own total-error gates.',
                runtime_seconds=time.perf_counter()-started,timings=timings,failures=failures,rows=rows)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps(dict(runtime_seconds=result['runtime_seconds'],rows=len(rows),failures=len(failures))))


if __name__=='__main__':
    main()
