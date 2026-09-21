"""Separate occupation quadrature from the already audited pointwise kernels."""
from pathlib import Path
import json
import sys
import time
import numpy as np
ROOT=Path(__file__).resolve().parents[4]
sys.path[:0]=[str(ROOT),str(ROOT/'docs/implementation/stage1_r2/review')]
from gamma_zero_reference import reference as independent_bcs,population
from pysnspd.experimental.energy_catalog import OccupationEnergyCatalog
from pysnspd.experimental.refined_cells import refined_count_catalog


def sha(path):
    import hashlib
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    started=time.perf_counter();out=Path(__file__).resolve().parent
    path=ROOT/'docs/implementation/stage1_r2/catalogs/occupation_catalog.npz'
    base=OccupationEnergyCatalog.load(path)
    candidates=[refined_count_catalog(base,k) for k in (1,2,4)]
    fields=[(.08,0.),(.17,0.),(.72,0.),(1.4,0.),(.72,1e-12),(.72,1e-8),
            (.72,.2),(.72,.719999),(.72,.720001),(.35,.65)]
    # Independent analytic BCS integrals are the regulator-layer anchor.
    # NonzeroGamma uses independently increased quadrature, while the spectral
    # kernels themselves are separately audited against60-digit w-plane roots.
    references=[refined_count_catalog(base,1,bulk_spacing=.01,order=q,edge_order=4*q)
                for q in (6,8)]
    profiles=('low_energy','near_edge_band')
    rows=[]
    for a,g in fields:
        if g==0:
            ref=[independent_bcs(a,profile,base.eta) for profile in profiles]
            wanted=[np.asarray(v[0]) for v in ref]
            ref_error=[max(v[1]) for v in ref]
        else:
            values=[]
            for catalog in references:
                kernels=np.asarray(catalog.energy_kernel(a,g))
                vacuum=np.asarray(catalog.vacuum.evaluate(a,g))
                values.append([vacuum+4*(kernels*catalog.count_weights)@np.array([population(profile,x) for x in catalog.count_nodes]) for profile in profiles])
            wanted=values[-1]
            ref_error=[float(np.max(abs(np.asarray(f)-c)/np.maximum(1,abs(f)))) for c,f in zip(values[0],values[1])]
        for mesh,level in zip(candidates,(1,2,4)):
            kernels=np.asarray(mesh.energy_kernel(a,g));vacuum=np.asarray(mesh.vacuum.evaluate(a,g))
            for i,profile in enumerate(profiles):
                p=np.array([population(profile,x) for x in mesh.count_nodes])
                actual=vacuum+4*(kernels*mesh.count_weights)@p
                error=abs(actual-wanted[i])/np.maximum(1,abs(wanted[i]))
                rows.append(dict(amplitude=a,gamma=g,profile=profile,refinement=level,
                    electron_states=len(mesh.count_nodes),actual=actual.tolist(),reference=wanted[i].tolist(),
                    scaled_errors=error.tolist(),reference_error_scaled_or_absolute=ref_error[i],
                    relative_Gamma_error=float(abs(actual[2]/wanted[i][2]-1)) if wanted[i][2] else None))
    maxima={str(level):np.max([r['scaled_errors'] for r in rows if r['refinement']==level],axis=0).tolist() for level in (1,2,4)}
    gates=dict(reference=max(r['reference_error_scaled_or_absolute'] for r in rows)<=2e-4,
        energies=all(v[0]<=1e-5 for v in maxima.values()),forces=all(max(v[1:])<=1e-3 for v in maxima.values()))
    data=dict(schema='pysnspd.stage2.complementary_moment_quadrature.v1',status='PASS' if all(gates.values()) else 'FAIL',
        gates=gates,maxima_by_refinement=maxima,rows=rows,
        method='Gamma0: independent analytic finite-eta adaptive integral. NonzeroGamma: separate order6/8 bulk and24/32 edge reference, with point kernels previously checked independently at60digits.',
        catalogue_sha256=sha(path),wrapper_sha256=sha(ROOT/'pysnspd/experimental/refined_cells.py'),
        analytic_reference_sha256=sha(ROOT/'docs/implementation/stage1_r2/review/gamma_zero_reference.py'),
        criteria_sha256=sha(ROOT/'docs/implementation/stage2/acceptance_criteria.json'),
        source_sha256=sha(Path(__file__)),runtime_seconds=time.perf_counter()-started)
    (out/'complementary_moment_quadrature.json').write_text(json.dumps(data,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(dict(status=data['status'],gates=gates,maxima=maxima,runtime_seconds=data['runtime_seconds'])))
    if not all(gates.values()):raise SystemExit(1)


if __name__=='__main__':main()
