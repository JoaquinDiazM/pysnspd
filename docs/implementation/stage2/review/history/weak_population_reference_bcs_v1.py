"""Independent B.7 collocation integral, then a native-energy hat weak moment.

No event constructor/enumerator is used. For normal/BCS spectra the partner
DOS is removed analytically by E'=sqrt(x'^2+Delta^2). This computes the physical
continuous-population reference, not a re-evaluation of geometric activities.
"""
from pathlib import Path
import argparse
import json
import sys
import time
import numpy as np
from continuous_reactions import ROOT,PROFILES,populations,sha
sys.path.insert(0,str(ROOT))
from pysnspd.experimental.energy_catalog import OccupationEnergyCatalog
from pysnspd.experimental.refined_cells import refined_count_catalog


def point_rhs(energy,a,emax,order):
    """All three analytic profiles, two channels, arbitrary fixed test energies."""
    z,w=np.polynomial.legendre.leggauss(order)
    result=np.zeros((3,2,len(energy)))
    breaks=np.array([0.,.025,.05,.1,.2,.4,.65,.8,1.,1.2,1.5,2.,2.5,3.,4.,6.,9.,13.])
    for kind in ('upper','lower','recombination'):
        if kind=='upper':lower=energy+.01;upper=np.minimum(emax,energy+4.)
        elif kind=='lower':lower=np.maximum(a,energy-4.);upper=energy-.01
        else:lower=np.maximum(a,.01-energy);upper=np.minimum(emax,4.-energy)
        valid=upper>lower
        xmin=np.sqrt(np.maximum(0.,lower*lower-a*a))
        xmax=np.sqrt(np.maximum(0.,upper*upper-a*a))
        for lo,hi in zip(breaks[:-1],breaks[1:]):
            begin=np.maximum(xmin,lo);end=np.minimum(xmax,hi)
            ids=np.flatnonzero(valid&(end>begin))
            if not len(ids):continue
            x=begin[ids,None]+(end-begin)[ids,None]*(z[None,:]+1)/2
            weights=(end-begin)[ids,None]*w[None,:]/2
            ei=energy[ids,None];ej=np.hypot(a,x)
            omega=ej-ei if kind=='upper' else ei-ej if kind=='lower' else ei+ej
            coherence=1+(1 if kind=='recombination' else -1)*a*a/(ei*ej)
            factor=.03*(omega/4.)**2*coherence*weights/4
            for ip,profile in enumerate(PROFILES):
                pi,n=populations(profile,ei,omega);pj,_=populations(profile,ej,omega)
                if kind=='upper':net=(1-pi)*pj*(1+n)-pi*(1-pj)*n;sign=1
                elif kind=='lower':net=(1-pj)*pi*(1+n)-pj*(1-pi)*n;sign=-1
                else:net=pi*pj*(1+n)-(1-pi)*(1-pj)*n;sign=-1
                result[ip,int(kind=='recombination'),ids]+=sign*np.sum(factor*net,axis=1)
    return result


def weak_reference(grid,a,outer_order,inner_order):
    # Native hats are linear in ENERGY, with exact breakpoints at native E_i.
    # In the ideal BCS reference no real states exist below a, even though the
    # finite-eta native grid includes explicitly regulated subgap states.
    z,w=np.polynomial.legendre.leggauss(outer_order)
    low=np.maximum(grid[:-1],a);high=grid[1:]
    valid=high>low
    xl=np.sqrt(np.maximum(0,low[valid]**2-a*a));xh=np.sqrt(high[valid]**2-a*a)
    x=xl[:,None]+(xh-xl)[:,None]*(z[None,:]+1)/2
    measure=(xh-xl)[:,None]*w[None,:]/2
    e=np.hypot(a,x)
    ids=np.flatnonzero(valid)
    beta=(e-grid[ids,None])/(grid[ids+1,None]-grid[ids,None])
    rhs=point_rhs(e.ravel(),a,float(grid[-1]),inner_order).reshape(3,2,len(ids),outer_order)
    values=np.zeros((3,2,len(grid)))
    values[:,:,ids]+=4*np.sum(rhs*measure[None,None,:,:]*(1-beta)[None,None,:,:],axis=-1)
    values[:,:,ids+1]+=4*np.sum(rhs*measure[None,None,:,:]*beta[None,None,:,:],axis=-1)
    return values


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--refinement',type=int,default=1)
    parser.add_argument('--output',type=Path)
    args=parser.parse_args();start=time.perf_counter();out=Path(__file__).resolve().parent
    basepath=ROOT/'docs/implementation/stage1_r2/catalogs/occupation_catalog.npz'
    cat=refined_count_catalog(OccupationEnergyCatalog.load(basepath),args.refinement)
    rows=[]
    for a in (0.,.72):
        grid=cat.energy_kernel(a,0.)[0]
        levels=[]
        for outer,inner in ((2,16),(4,32),(8,64)):
            tick=time.perf_counter();values=weak_reference(grid,a,outer,inner)
            levels.append(values)
            print(json.dumps(dict(amplitude=a,outer=outer,inner=inner,seconds=time.perf_counter()-tick)),flush=True)
        for ip,profile in enumerate(PROFILES):
            channels={}
            for ic,channel in enumerate(('scattering','recombination')):
                expected=levels[-1][ip,ic];scale=max(float(np.sum(abs(expected))),1e-100)
                channels[channel]=dict(number_rhs_reference=expected.tolist(),reference_L1=scale,
                    last_refinement_relative_L1=float(np.sum(abs(expected-levels[-2][ip,ic]))/scale),
                    earlier_refinement_relative_L1=float(np.sum(abs(expected-levels[0][ip,ic]))/scale))
            rows.append(dict(amplitude=a,gamma=0.,profile=profile,energy_nodes=grid.tolist(),channels=channels))
    data=dict(schema='pysnspd.stage2.independent_weak_population.v1',
        status='IDEAL_CONTINUOUS_NORMAL_BCS_REFERENCE',
        definition='number_rhs_i = integral phi_i(E) *4*dp(E)/dt dx; native hats linear in E; independent B.7 partner integral at fixed E, then outer hat integral',
        limitation='Ideal eta0 reference; finite-eta and interpolation of the analytic population by geometric nodal activities are separate measured defects',
        fields=[[0.,0.],[.72,0.]],count_refinement=args.refinement,electron_states=len(cat.count_nodes),
        count_nodes=cat.count_nodes.tolist(),count_weights=cat.count_weights.tolist(),
        orders=[[2,16],[4,32],[8,64]],rows=rows,rate_prefactor=1.,omega_support=[.01,4.],
        catalog_sha256=sha(basepath),wrapper_sha256=sha(ROOT/'pysnspd/experimental/refined_cells.py'),
        source_sha256=sha(Path(__file__)),criteria_sha256=sha(ROOT/'docs/implementation/stage2/acceptance_criteria.json'),
        runtime_seconds=time.perf_counter()-start)
    target=args.output or out/f'weak_population_reference_{args.refinement}.json'
    target.write_text(json.dumps(data,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(dict(output=str(target),runtime_seconds=data['runtime_seconds'],
        worst_reference_refinement=max(ch['last_refinement_relative_L1'] for row in rows for ch in row['channels'].values()))))


if __name__=='__main__':main()
