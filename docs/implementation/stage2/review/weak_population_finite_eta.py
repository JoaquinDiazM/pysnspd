"""Independent finite-eta BCS weak B.7 reference on fixed native energy hats.

This is a separate regulator diagnostic; the published ideal references remain
unchanged. Neither an event enumerator nor an interpolated spectral kernel is
used in the integration. The catalog supplies only the native count coordinates.
"""
from pathlib import Path
import json
import sys
import time
import numpy as np
from continuous_reactions import ROOT, PROFILES, populations, sha
sys.path.insert(0, str(ROOT))
from pysnspd.experimental.energy_catalog import OccupationEnergyCatalog
from pysnspd.experimental.refined_cells import refined_count_catalog


def energy_at_x(x, a, eta):
    return x * np.sqrt(1 + a*a/(x*x + eta*eta))


def count_at_energy(e, a, eta):
    return np.sqrt((np.asarray(e) + 1j*eta)**2 - a*a).real


def ratio_at_x(x, a, eta):
    return a*x/np.sqrt((x*x+eta*eta)*(x*x+eta*eta+a*a))


def point_rhs(xi, a, eta, emax, order):
    """Integrate each partner state in its exact finite-eta count coordinate."""
    energy=energy_at_x(xi,a,eta)
    ratio=ratio_at_x(xi,a,eta)
    z,w=np.polynomial.legendre.leggauss(order)
    result=np.zeros((3,2,len(energy)))
    breaks=np.unique(np.r_[0.,np.geomspace(1e-10,.01,18),
        .025,.05,.1,.2,.4,.65,.8,1.,1.2,1.5,2.,2.5,3.,4.,6.,9.,13.])
    for kind in ('upper','lower','recombination'):
        if kind=='upper': lower=energy+.01;upper=np.minimum(emax,energy+4.)
        elif kind=='lower': lower=np.maximum(0.,energy-4.);upper=energy-.01
        else: lower=np.maximum(0.,.01-energy);upper=np.minimum(emax,4.-energy)
        valid=upper>lower
        xmin=count_at_energy(np.maximum(0.,lower),a,eta)
        xmax=count_at_energy(np.maximum(0.,upper),a,eta)
        for lo,hi in zip(breaks[:-1],breaks[1:]):
            begin=np.maximum(xmin,lo);end=np.minimum(xmax,hi)
            ids=np.flatnonzero(valid&(end>begin))
            if not len(ids): continue
            x=begin[ids,None]+(end-begin)[ids,None]*(z[None,:]+1)/2
            weights=(end-begin)[ids,None]*w[None,:]/2
            ei=energy[ids,None];ej=energy_at_x(x,a,eta)
            omega=ej-ei if kind=='upper' else ei-ej if kind=='lower' else ei+ej
            coherence=1+(1 if kind=='recombination' else -1)*ratio[ids,None]*ratio_at_x(x,a,eta)
            factor=.03*(omega/4.)**2*coherence*weights/4
            for ip,profile in enumerate(PROFILES):
                pi,n=populations(profile,ei,omega);pj,_=populations(profile,ej,omega)
                if kind=='upper': net=(1-pi)*pj*(1+n)-pi*(1-pj)*n;sign=1
                elif kind=='lower': net=(1-pj)*pi*(1+n)-pj*(1-pi)*n;sign=-1
                else: net=pi*pj*(1+n)-(1-pi)*(1-pj)*n;sign=-1
                result[ip,int(kind=='recombination'),ids]+=sign*np.sum(factor*net,axis=1)
    return result


def weak_reference(counts,a,eta,outer,inner):
    z,w=np.polynomial.legendre.leggauss(outer)
    grid=energy_at_x(counts,a,eta)
    x=counts[:-1,None]+np.diff(counts)[:,None]*(z[None,:]+1)/2
    measure=np.diff(counts)[:,None]*w[None,:]/2
    energy=energy_at_x(x,a,eta)
    beta=(energy-grid[:-1,None])/np.diff(grid)[:,None]
    rhs=point_rhs(x.ravel(),a,eta,float(grid[-1]),inner).reshape(3,2,len(counts)-1,outer)
    result=np.zeros((3,2,len(counts)))
    result[:,:,:-1]+=4*np.sum(rhs*measure[None,None,:,:]*(1-beta)[None,None,:,:],axis=-1)
    result[:,:,1:]+=4*np.sum(rhs*measure[None,None,:,:]*beta[None,None,:,:],axis=-1)
    return grid,result


def main():
    start=time.perf_counter();out=Path(__file__).resolve().parent
    catalog_path=ROOT/'docs/implementation/stage1_r2/catalogs/occupation_catalog.npz'
    catalog=refined_count_catalog(OccupationEnergyCatalog.load(catalog_path),1)
    a=.72;eta=1e-8;levels=[];orders=((2,16),(4,32),(8,64))
    for outer,inner in orders:
        tick=time.perf_counter()
        grid,result=weak_reference(catalog.count_nodes,a,eta,outer,inner)
        levels.append(result)
        print(json.dumps(dict(outer=outer,inner=inner,seconds=time.perf_counter()-tick)),flush=True)
    rows=[]
    for ip,profile in enumerate(PROFILES):
        channels={}
        for ic,channel in enumerate(('scattering','recombination')):
            expected=levels[-1][ip,ic];scale=float(np.sum(abs(expected)))
            channels[channel]=dict(number_rhs_reference=expected.tolist(),reference_L1=scale,
                last_refinement_relative_L1=float(np.sum(abs(expected-levels[-2][ip,ic]))/scale),
                earlier_refinement_relative_L1=float(np.sum(abs(expected-levels[0][ip,ic]))/scale))
        rows.append(dict(amplitude=a,gamma=0.,profile=profile,energy_nodes=grid.tolist(),channels=channels))
    worst=max(ch['last_refinement_relative_L1'] for row in rows for ch in row['channels'].values())
    data=dict(schema='pysnspd.stage2.independent_weak_population.finite_eta.v1',
        status='PASS_REFERENCE_BUDGET' if worst<=2e-4 else 'FAIL_REFERENCE_BUDGET',
        definition='Independent B.7 partner and outer count integration; native hats linear in exact finite-eta BCS energy; analytic continuous p(E), n(Omega)',
        limitation='Regulator diagnostic for delta=.72, Gamma=0, 630 states, Omega in [.01,4]. It does not certify the selected .005-cut trajectory or dynamic convergence.',
        eta=eta,rate_prefactor=1.,omega_support=[.01,4.],count_refinement=1,
        electron_states=len(catalog.count_nodes),count_nodes=catalog.count_nodes.tolist(),
        count_weights=catalog.count_weights.tolist(),orders=orders,rows=rows,
        worst_reference_refinement=worst,catalog_sha256=sha(catalog_path),
        wrapper_sha256=sha(ROOT/'pysnspd/experimental/refined_cells.py'),
        source_sha256=sha(Path(__file__)),criteria_sha256=sha(ROOT/'docs/implementation/stage2/acceptance_criteria.json'),
        runtime_seconds=time.perf_counter()-start)
    target=out/'weak_population_reference_finite_eta_bcs_1.json'
    target.write_text(json.dumps(data,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(dict(status=data['status'],output=str(target),runtime_seconds=data['runtime_seconds'],worst_reference_refinement=worst)),flush=True)


if __name__=='__main__':main()
