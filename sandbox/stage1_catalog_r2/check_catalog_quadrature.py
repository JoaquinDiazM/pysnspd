"""Check the exact saved count rule against resolved, same-eta quadrature."""
from pathlib import Path
import argparse
import hashlib
import json
import sys
import time

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT))
import numpy as np
from scipy.special import expit,xlogy
from pysnspd.experimental.energy_catalog import OccupationEnergyCatalog,energy_at_count_batch,retarded_spectrum_batch,vacuum_state
from build_catalog import population


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('catalog',type=Path)
    parser.add_argument('--reference',type=Path,default=ROOT/'docs/implementation/stage1_r2/quadrature_diagnostics.json')
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args(); started=time.perf_counter()
    table=OccupationEnergyCatalog.load(args.catalog)
    diagnostics=json.loads(args.reference.read_text())
    refrows=[r for r in diagnostics['rows'] if r['eta']==table.eta and r['order']==32 and r['cutoff']==12]
    if not refrows:
        raise ValueError('no same-eta resolved reference for this catalogue')
    check24={(r['amplitude'],r['gamma'],r['profile']):r for r in diagnostics['rows']
             if r['eta']==table.eta and r['order']==24 and r['cutoff']==12}
    cache={}; rows=[]
    for ref in refrows:
        a,g=ref['amplitude'],ref['gamma']; key=(a,g)
        if key not in cache:
            energy=energy_at_count_batch(table.count_nodes,delta=a,gamma=g,eta=table.eta)
            c,s=retarded_spectrum_batch(energy,delta=a,gamma=g,eta=table.eta)
            cache[key]=np.array([energy,s.imag/c.real,-s.real*s.imag/c.real])
        kernels=cache[key]; temperature=float(ref['profile'][3:]) if ref['profile'].startswith('FD_') else 0.
        p=expit(-kernels[0]/temperature) if temperature else population(ref['profile'],table.count_nodes)
        value=np.asarray(vacuum_state(a,g))+4*np.sum(kernels*(p*table.count_weights)[None,:],axis=1)
        if temperature:
            value[0]+=4*temperature*np.dot(table.count_weights,xlogy(p,p)+xlogy(1-p,1-p))
        reference=np.asarray(ref['values']); error=abs(value-reference)
        uncertainty=abs(reference-check24[(a,g,ref['profile'])]['values'])
        rows.append(dict(amplitude=a,gamma=g,profile=ref['profile'],value=value.tolist(),reference=reference.tolist(),
                         scaled_error=(error/np.maximum(1,abs(reference))).tolist(),
                         relative_error=[float(error[i]/abs(reference[i])) if reference[i]!=0 else None for i in range(3)],
                         reference_order24_to32_change=uncertainty.tolist()))
    worst=[]
    for i,name in enumerate(('energy','amplitude_force','Gamma_response')):
        scaled=max(rows,key=lambda r:r['scaled_error'][i]); relative=max(rows,key=lambda r:r['relative_error'][i] or 0.)
        worst.append(dict(observable=name,max_scaled_error=scaled['scaled_error'][i],
                          scaled_case=[scaled['amplitude'],scaled['gamma'],scaled['profile']],
                          max_relative_error=relative['relative_error'][i],
                          relative_case=[relative['amplitude'],relative['gamma'],relative['profile']]))
    result=dict(schema='pysnspd.stage1_r2.exact_count_rule_check.v1',runtime_seconds=time.perf_counter()-started,
                catalog_sha256=hashlib.sha256(args.catalog.read_bytes()).hexdigest(),
                reference_sha256=hashlib.sha256(args.reference.read_bytes()).hexdigest(),
                source_sha256=hashlib.sha256((ROOT/'pysnspd/experimental/energy_catalog.py').read_bytes()).hexdigest(),
                eta=table.eta,count_nodes=len(table.count_nodes),count_support=float(np.sum(table.count_weights)),
                interpretation='Exact saved count rule, without field interpolation; same eta in candidate and reference.',
                worst_errors=worst,rows=rows)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps(dict(runtime_seconds=result['runtime_seconds'],worst_errors=worst),indent=2))


if __name__=='__main__':
    main()
