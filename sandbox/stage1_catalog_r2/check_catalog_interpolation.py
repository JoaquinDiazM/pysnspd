"""Measure only field interpolation, using the saved eta and exact count rule."""
from pathlib import Path
import argparse
import hashlib
import json
import sys
import time

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT))
import numpy as np
from pysnspd.experimental.energy_catalog import OccupationEnergyCatalog
from build_catalog import population,direct


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('catalog',type=Path)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();started=time.perf_counter()
    criteria_path=ROOT/'docs/implementation/stage1_r2/acceptance_criteria.json'
    criteria=json.loads(criteria_path.read_text())
    table=OccupationEnergyCatalog.load(args.catalog)
    fields=[(f'fixed_{i:02d}',*point) for i,point in enumerate(criteria['mandatory_field_points'])]
    rng=np.random.default_rng(7431)
    for i in range(12):
        a=float(rng.uniform(.085,1.49));g=float(10**rng.uniform(-12,-1.5) if i<6 else rng.uniform(.031,1.19))
        fields.append((f'offnode_{i:02d}',a,g))
    rows=[];failures=[];limits=[1e-5,1e-3,1e-3]
    for field,a,g in fields:
        for profile in criteria['mandatory_population_profiles']:
            name=profile['id'];p=population(name,table.count_nodes)
            try:
                value=np.asarray(table.evaluate(a,g,p))
            except FloatingPointError as error:
                failures.append(dict(field_id=field,amplitude=a,gamma=g,profile=name,reason=str(error)))
                continue
            reference=direct(a,g,table.count_nodes,table.count_weights,p,table.eta)
            error=abs(value-reference);scaled=error/np.maximum(1,abs(reference))
            relative=[float(error[k]/abs(reference[k])) if reference[k]!=0 else None for k in range(3)]
            row=dict(case_id=field+'_'+name,amplitude=a,gamma=g,profile=name,value=value.tolist(),
                     reference=reference.tolist(),scaled_error=scaled.tolist(),relative_error=relative)
            rows.append(row)
            if np.any(scaled>limits) or any(relative[k] is not None and relative[k]>1e-3 for k in (1,2)):
                failures.append(dict(case_id=row['case_id'],scaled_error=scaled.tolist(),relative_error=relative))
    worst=[]
    for i,name in enumerate(('energy','amplitude_force','Gamma_response')):
        scaled=max(rows,key=lambda r:r['scaled_error'][i]);relative=max(rows,key=lambda r:r['relative_error'][i] or 0.)
        worst.append(dict(observable=name,max_scaled_error=scaled['scaled_error'][i],scaled_case=scaled['case_id'],
                          max_relative_error=relative['relative_error'][i],relative_case=relative['case_id']))
    result=dict(schema='pysnspd.stage1_r2.saved_interpolation_check.v1',runtime_seconds=time.perf_counter()-started,
                catalog_sha256=hashlib.sha256(args.catalog.read_bytes()).hexdigest(),
                source_sha256=hashlib.sha256((ROOT/'pysnspd/experimental/energy_catalog.py').read_bytes()).hexdigest(),
                criteria_sha256=hashlib.sha256(criteria_path.read_bytes()).hexdigest(),
                eta=table.eta,shape=list(table.excitation_energies.shape),status='FAIL' if failures else 'PASS_SAMPLED',
                interpretation='Identical count nodes, weights, occupation and eta; only field interpolation differs.',
                worst_errors=worst,failures=failures,rows=rows)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps(dict(runtime_seconds=result['runtime_seconds'],status=result['status'],worst_errors=worst),indent=2))


if __name__=='__main__':
    main()
