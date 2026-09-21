"""Independent algebra and support checks for the two-cell prototype.

Partition-of-unity identities integrate the reconstructed E(x) directly;
matrix exponentials check the advertised pair solution without its formula.
This is not a convergence result for continuum D.17.
"""
from pathlib import Path
import hashlib
import json
import sys
import time

import numpy as np
from scipy.linalg import expm

ROOT=Path(__file__).resolve().parents[4]
sys.path.insert(0,str(ROOT))
from pysnspd.experimental.energy_catalog import OccupationEnergyCatalog
from pysnspd.experimental.cell_validation import ElectronicCell,EnergyFacePair


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    started=time.perf_counter()
    path=ROOT/'docs/implementation/stage1_r2/catalogs/occupation_catalog.npz'
    cat=OccupationEnergyCatalog.load(path)
    rows=[]
    for fields in [((.72,.2),(1.,.03)),((.35,.65),(.72,.002)),((1.,0.),(.72,0.))]:
        cells=[ElectronicCell(cat,*field) for field in fields]
        emin=min(cell.energies[0] for cell in cells)/2
        emax=max(cell.energies[-1] for cell in cells)
        common_upper=min(cell.energies[-1] for cell in cells)
        axis=np.geomspace(emin,common_upper,257)
        if emax>common_upper:
            axis=np.r_[axis,emax]
        pair=EnergyFacePair(*cells,axis,.13)
        moment_rows=[]
        for i,cell in enumerate(cells):
            xx=np.r_[0.,cat.count_nodes]
            ee=np.r_[0.,cell.energies]
            low,high=max(ee[0],axis[0]),min(ee[-1],axis[-1])
            split=np.r_[low,ee[(ee>low)&(ee<high)],high]
            counts=np.interp(split,ee,xx)
            count_ref=counts[-1]-counts[0]
            energy_ref=float(np.sum((split[:-1]+split[1:])*np.diff(counts)/2))
            moment_rows.append(dict(count_sum_error=abs(float(np.sum(pair.count_moments[i]))-count_ref),
                                    energy_sum_error=abs(float(np.sum(pair.energy_moments[i]))-energy_ref),
                                    linear_reproduction_error=abs(float(np.dot(axis[:pair.common_size],pair.count_moments[i,:pair.common_size])-np.sum(pair.energy_moments[i,:pair.common_size]))),
                                    omitted_lower_count=float(counts[0]),
                                    unextrapolated_last_count=float(counts[-1])))
        f=np.array([.08+.12*np.exp(-axis/.3),.03+.02*np.exp(-((axis-.5)/.3)**2)])
        rhs=pair.rhs(f)
        moment=float(4*np.sum(pair.energy_moments*rhs))
        common=np.array([f[0],f[0]])
        steady=float(np.max(abs(pair.rhs(common))))
        chosen=np.flatnonzero(pair.active)[::max(1,int(np.count_nonzero(pair.active)/25))]
        exact=pair.exact(f,.4)
        errors=[]
        for j in chosen:
            cl,cr=pair.capacities[:,j]
            k=pair.conductance[j]
            generator=np.array([[-k/cl,k/cl],[k/cr,-k/cr]])
            reference=expm(.4*generator)@f[:,j]
            errors.append(float(np.max(abs(exact[:,j]-reference))))
        e0=float(np.sum(pair.reconstructed_energy(f)))
        e1=float(np.sum(pair.reconstructed_energy(exact)))
        private=axis>common_upper
        private_error=float(np.max(abs(exact[:,private]-f[:,private]))) if np.any(private) else 0.
        rows.append(dict(fields=fields,shared_energy_nodes=len(axis),
                         moment_identities=moment_rows,instantaneous_energy_moment=abs(moment),
                         common_distribution_rhs=steady,
                         matrix_exponential_population_abs_error=max(errors),
                         exact_trajectory_energy_drift=abs(e1-e0),
                         private_tail_population_change=private_error,
                         pauli_minimum=float(np.min(exact)),pauli_maximum=float(np.max(exact)),
                         initial_quasiparticle_count=pair.quasiparticle_count(f),
                         final_quasiparticle_count=pair.quasiparticle_count(exact)))
    maxima={key:max(row[key] for row in rows) for key in
            ('instantaneous_energy_moment','common_distribution_rhs','matrix_exponential_population_abs_error',
             'exact_trajectory_energy_drift','private_tail_population_change')}
    maxima['moment_identity_abs_error']=max(abs(value) for row in rows for moment in row['moment_identities']
        for key,value in moment.items() if key.endswith('error'))
    passed=all(value<=1e-10 for value in maxima.values()) and all(
        row['pauli_minimum']>=0 and row['pauli_maximum']<=1 for row in rows)
    result=dict(schema='pysnspd.stage1_closure.independent_transport.v1',status='PASS' if passed else 'FAIL',
                scope='Algebra of reconstructed-energy transport only. Native-R2 bias and continuum/remap convergence require separate measured evidence.',
                criteria_sha256=sha(ROOT/'docs/implementation/stage1_closure/acceptance_criteria.json'),
                catalog_sha256=sha(path),cell_source_sha256=sha(ROOT/'pysnspd/experimental/cell_validation.py'),
                reviewer_source_sha256=sha(Path(__file__)),summary=maxima,rows=rows,
                runtime_seconds=time.perf_counter()-started)
    Path(__file__).with_name('transport_checks.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({key:result[key] for key in ('status','summary','runtime_seconds')},indent=2))


if __name__=='__main__':
    main()
