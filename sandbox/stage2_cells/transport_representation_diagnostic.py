"""Diagnostic only: separate geometric attenuation from logit interpolation.

The normalized geometric rates below are evaluated ONLY for strictly interior
smooth profiles. They are not an accepted operator: normalization has no unique
continuous extension at native brackets containing both zero and one.
"""
from pathlib import Path
import hashlib
import json
import sys
import time
import numpy as np
from scipy.special import expit, logit

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT))
from transport_diagnostics import FIELDS, population, reference, population_norm
from pysnspd.experimental.cell_validation import ElectronicCell
from pysnspd.experimental.energy_catalog import OccupationEnergyCatalog
from pysnspd.experimental.cell_transport import NativeTransportEvents


def deposit(op, cells, flux):
    out=[]
    for side,cell in enumerate(cells):
        rate=np.bincount(op.indices[side].ravel(),
            weights=(op.barycentric[side]*flux[:,None]).ravel(),minlength=len(cell.energies))
        out.append((2*side-1)*rate/(4*cell.weights))
    return out


def main():
    started=time.perf_counter()
    path=ROOT/'docs/implementation/stage1_r2/catalogs/occupation_catalog.npz'
    catalog=OccupationEnergyCatalog.load(path)
    rows=[]
    for field,parameters in FIELDS.items():
        cells=[ElectronicCell(catalog,*p) for p in parameters]
        op=NativeTransportEvents(*cells,gauss_order=4)
        for profile in ['thermal','nonthermal']:
            native=[population(c.energies,profile,i) for i,c in enumerate(cells)]
            reconstructed=[]; factors=[]
            for side,p in enumerate(native):
                selected=p[op.indices[side]]; b=op.barycentric[side]
                if np.any(selected<=0) or np.any(selected>=1):
                    raise ValueError('This diagnostic requires strictly interior populations')
                la=np.sum(b*np.log(selected),axis=1)
                lb=np.sum(b*np.log1p(-selected),axis=1)
                reconstructed.append(expit(la-lb))
                factors.append(np.exp(la)+np.exp(lb))
            unattenuated=op.conductance*(reconstructed[0]-reconstructed[1])
            exact=op.conductance*(population(op.energies,profile,0)-population(op.energies,profile,1))
            raw=op.rates(native)
            ref,power,_=reference(cells,profile,.03,48)
            versions={}
            for name,flux in [('original',raw),('normalized_interior_only',unattenuated),('analytic_population_same_Gauss4',exact)]:
                versions[name]={**population_norm(cells,deposit(op,cells,flux),ref),
                    'power_relative_error':abs(np.dot(op.energies,flux)-power)/abs(power)}
            activity=np.abs(unattenuated)
            versions['geometric_factor']={
                'minimum_Z_left_times_Z_right':float(np.min(factors[0]*factors[1])),
                'activity_weighted_attenuation':float(np.dot(activity,1-factors[0]*factors[1])/np.sum(activity)),
                'factorization_absolute_residual':float(np.max(np.abs(raw-unattenuated*factors[0]*factors[1])))}
            rows.append({'field':field,'profile':profile,'decomposition':versions})
    source=ROOT/'pysnspd/experimental/cell_transport.py'
    result={'schema':'pysnspd.stage2.transport-representation-diagnostic.v1',
        'scope':'Diagnostic only; normalized interior geometric products have no unique continuous boundary limit at brackets [0,1]. No operator is changed or admitted.',
        'operator_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),
        'runner_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'catalog_sha256':hashlib.sha256(path.read_bytes()).hexdigest(),
        'rows':rows,'runtime_seconds':time.perf_counter()-started}
    out=ROOT/'docs/implementation/stage2/pilot_initial/transport_representation_diagnostic.json'
    out.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    print(json.dumps(result,indent=2))


if __name__=='__main__':main()
