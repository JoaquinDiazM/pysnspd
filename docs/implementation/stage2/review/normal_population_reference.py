"""Independent normal-state collocation RHS from one-dimensional Omega integrals.

This tests the native pair quadrature at its own electron collocation points.
It is not a Galerkin/hat reference for a barycentrically projected operator.
The latter requires its own weak reference; the distinction is explicit here.
"""
from pathlib import Path
import json
import sys
import time
import numpy as np
from continuous_reactions import ROOT,PROFILES,populations,sha

sys.path.insert(0,str(ROOT))
from pysnspd.experimental.energy_catalog import OccupationEnergyCatalog
from pysnspd.experimental.cell_validation import ElectronicCell
from pysnspd.experimental.kinetic_events import PhononGrid,ElectronPhononEvents


def collocation_rhs(energies,profile,order=64,omega_min=.01,omega_max=4.):
    """Normalized prefactor1, p capacity4 dx, normal N1=1, R2=0."""
    z,w=np.polynomial.legendre.leggauss(order)
    result=np.zeros((2,len(energies)))
    emax=float(energies[-1])
    for i,ei in enumerate(energies):
        pi=populations(profile,ei,1.)[0]
        for kind in ('upper','lower','recombination'):
            low=omega_min;high=omega_max
            if kind=='upper':high=min(high,emax-ei)
            elif kind=='lower':high=min(high,ei)
            else:low=max(low,ei);high=min(high,ei+emax)
            if high<=low:continue
            # Fixed phonon features plus translated narrow electronic packet.
            cuts=sorted(set([low,high]+[v for v in [.03,.1,.25,.5,1.,1.35,1.7,2.,2.35,2.7,3.,
                1.0-ei,1.3-ei,1.6-ei,ei-1.6,ei-1.3,ei-1.,ei+1.,ei+1.3,ei+1.6] if low<v<high]))
            for lo,hi in zip(cuts[:-1],cuts[1:]):
                om=lo+(hi-lo)*(z+1)/2;measure=(hi-lo)*w/2
                ej=ei+om if kind=='upper' else ei-om if kind=='lower' else om-ei
                pj,n=populations(profile,ej,om)
                if kind=='upper':net=(1-pi)*pj*(1+n)-pi*(1-pj)*n;sign=1
                elif kind=='lower':net=(1-pj)*pi*(1+n)-pj*(1-pi)*n;sign=-1
                else:net=pi*pj*(1+n)-(1-pi)*(1-pj)*n;sign=-1
                result[int(kind=='recombination'),i]+=sign*.25*np.dot(measure,.03*(om/4)**2*net)
    return result


def main():
    started=time.perf_counter();path=Path(__file__).resolve().parent
    catalogue=ROOT/'docs/implementation/stage1_r2/catalogs/occupation_catalog.npz'
    cell=ElectronicCell(OccupationEnergyCatalog.load(catalogue),0,0)
    # Only coupling support is used below; exact n(Omega) removes Bose-grid error.
    phonons=PhononGrid(np.array([.01,4.]),np.ones(2),lambda om:.03*(om/4)**2,(.01,4.),'reference-only')
    events=ElectronPhononEvents.from_cell(cell,phonons,rate_prefactor=1.)
    rows=[]
    for profile in PROFILES:
        refinements=[collocation_rhs(cell.energies,profile,order) for order in (32,64,128)]
        p,n=populations(profile,cell.energies,events.omega)
        pi,pj=p[events.index_i],p[events.index_j]
        forward=np.where(events.recombination,pi*pj,(1-pi)*pj)*(1+n)
        reverse=np.where(events.recombination,(1-pi)*(1-pj),pi*(1-pj))*n
        rate=events.coefficients*(forward-reverse)
        channels={}
        for index,kind in enumerate(('scattering','recombination')):
            mask=events.recombination if index else ~events.recombination
            actual=events.rhs_from_rates(np.where(mask,rate,0.))[0]
            reference=refinements[-1][index]
            scale=float(np.dot(4*cell.weights,abs(reference)))
            channels[kind]=dict(reference=reference.tolist(),native_pair_rhs=actual.tolist(),
                reference_capacity_weighted_L1=scale,
                reference_refinement_relative_L1=float(np.dot(4*cell.weights,abs(reference-refinements[-2][index]))/scale),
                native_pair_relative_L1=float(np.dot(4*cell.weights,abs(actual-reference))/scale))
        rows.append(dict(profile=profile,channels=channels))
    data=dict(schema='pysnspd.stage2.independent_normal_collocation.v1',
        status='DIAGNOSTIC_NATIVE_PAIR_COLLOCATION_NOT_PROJECTED_OPERATOR_VERDICT',
        method='Independent Omega integration of B.7 at each E_i, normal DOS, exact analytic populations and Bose factors',
        orders=[32,64,128],support=[float(cell.energies[0]),float(cell.energies[-1])],
        physical_lower_in_reference=0.,rate_prefactor=1.,rows=rows,
        energy=cell.energies.tolist(),electron_capacities=(4*cell.weights).tolist(),
        catalogue_sha256=sha(catalogue),source_sha256=sha(Path(__file__)),
        event_source_sha256=sha(ROOT/'pysnspd/experimental/kinetic_events.py'),
        criteria_sha256=sha(ROOT/'docs/implementation/stage2/acceptance_criteria.json'),
        runtime_seconds=time.perf_counter()-started)
    (path/'normal_population_reference.json').write_text(json.dumps(data,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(dict(runtime_seconds=data['runtime_seconds'],
        errors=[(r['profile'],ch,v['native_pair_relative_L1'],v['reference_refinement_relative_L1']) for r in rows for ch,v in r['channels'].items()])))


if __name__=='__main__':main()
