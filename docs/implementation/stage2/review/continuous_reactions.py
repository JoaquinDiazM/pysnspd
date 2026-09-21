"""Independent continuous B.5-B.9 reference in (E,Omega), not an event list.

Normal and Gamma=0 BCS spectra are analytic. The gapless spectrum uses the
real eta=0 Usadel-angle parametrization, with original-equation residual checks.
BCS substitutions remove the DOS endpoint singularities analytically. The
reference prefactor is8*pi*Delta0*t_ref/hbar=1; R therefore carries1/2.
No experimental implementation is imported. Finite-eta bias is a separate test.
"""
from pathlib import Path
import argparse
import hashlib
import json
import math
import sys
import time

import numpy as np
from scipy.special import expit

ROOT=Path(__file__).resolve().parents[4]
sys.path.insert(0,str(ROOT/'docs/implementation/stage1_r2/review'))
from causal_reference import parametric_state,t_at_count

OUT=Path(__file__).resolve().parent
PROFILES=('hot_electrons','phonon_bubble','nonthermal')
FIELDS=((0.,0.),(.72,0.),(.35,.65))
OBSERVABLES=('phonon_power_net','phonon_power_gross','phonon_number_rate_net',
             'phonon_number_rate_gross','electronic_second_moment_rate')


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def populations(profile,e,omega):
    if profile=='hot_electrons':
        return expit(-e/.4),1/np.expm1(omega/.15)
    if profile=='phonon_bubble':
        return expit(-e/.12),1/np.expm1(omega/.12)+.3*np.exp(-((omega-2)/.35)**2)
    return (.12*np.exp(-e/.25)+.1*np.exp(-((e-1.3)/.25)**2),
            1/np.expm1(omega/.18)+.1*np.exp(-((omega-1.7)/.4)**2))


def real_angle(t,a,g):
    """Vector form of the real cubic for a gapless, ideal causal spectrum."""
    b=a/g;cap=min(1.,b)
    s=cap*np.exp(-t)
    rhs=s*((b-cap)+cap*(-np.expm1(-t)))
    r=np.minimum(rhs/(s*s+2),np.cbrt(rhs))
    for _ in range(12):
        r-=(r*(r*(r+3)+s*s+2)-rhs)/(3*r*(r+2)+s*s+2)
    y=1+r;h=r*(2+r)
    k=(1-cap*cap)+cap*cap*(-np.expm1(-2*t))
    v=np.sqrt(h)
    energy=g*v*(h+k)/s
    count=g*y*v*k**1.5/s
    rp=s*(2*s*y-b)/(3*y*y+s*s-1)
    derivative=count*(1+rp*(1/y+y/h)+3*s*s/k)/(np.sqrt(k)*y)
    return energy,derivative,s,y,v,k


def gapless_spectrum(energy,a,g):
    energy=np.asarray(energy,float)
    if np.any(energy<=0) or g<=a:
        raise ValueError('positive energy and gapless field required')
    lower=np.zeros_like(energy)
    upper=np.log1p((energy/g)**2)+4
    t=np.log1p((energy/g)**2)
    for _ in range(50):
        obtained,derivative,*_=real_angle(t,a,g)
        relative=abs(obtained-energy)/energy
        if np.max(relative)<3e-13:
            break
        lower=np.where(obtained<energy,t,lower)
        upper=np.where(obtained>=energy,t,upper)
        candidate=t-(obtained-energy)/derivative
        updated=np.where((candidate>lower)&(candidate<upper),candidate,(lower+upper)/2)
        # Keep individually converged coordinates fixed. Re-bisecting an exact
        # Newton solution merely because another array entry is still active
        # can destroy its achieved precision and make a global test oscillate.
        t=np.where(relative<3e-13,t,updated)
    else:
        raise ArithmeticError(f'gapless real-angle inversion failed: {relative.max()}')
    obtained,_,s,y,v,k=real_angle(t,a,g)
    c=np.sqrt(k)*y-1j*s*v
    anomalous=s*y+1j*np.sqrt(k)*v
    residual=a*c-(g*c-1j*energy)*anomalous
    error=float(np.max(abs(residual)/np.maximum(1,energy)))
    if error>1e-11:
        raise ArithmeticError(f'causal reference residual {error}')
    return c.real,anomalous.imag,error


def segments(bounds,order,*,square=False):
    z,w=np.polynomial.legendre.leggauss(order)
    result=[]
    for low,high in zip(bounds[:-1],bounds[1:]):
        if high<=low:
            continue
        if square:
            lo,hi=np.sqrt(low),np.sqrt(high)
            t=lo+(z+1)*(hi-lo)/2
            result.append((t*t,w*(hi-lo)*t))
        else:
            result.append((low+(z+1)*(high-low)/2,w*(high-low)/2))
    return np.concatenate([row[0] for row in result]),np.concatenate([row[1] for row in result])


def project_hats(omega,rate,weights,axis):
    upper=np.searchsorted(axis,omega,side='left')
    lower=upper-1
    beta=(omega-axis[lower])/(axis[upper]-axis[lower])
    return (np.bincount(lower,weights=weights*(1-beta)*rate,minlength=len(axis))+
            np.bincount(upper,weights=weights*beta*rate,minlength=len(axis)))


def evaluate_field(a,g,order,omega_min=.01,xcut=12.):
    emax=math.hypot(xcut,a) if g==0 else parametric_state(t_at_count(xcut,a,g),a,g)[1]
    spectral_error=0.
    result={p:{ch:np.zeros(len(OBSERVABLES)) for ch in ('scattering','recombination')} for p in PROFILES}
    projected={p:{ch:np.zeros(17) for ch in ('scattering','recombination')} for p in PROFILES}
    z,w=np.polynomial.legendre.leggauss(order)
    energy_breaks=np.array([0.,.03,.1,.25,.5,.8,1.1,1.5,2.,3.,5.,8.,14.])
    for channel in ('scattering','recombination'):
        low=max(omega_min,2*a) if channel=='recombination' and g==0 else omega_min
        if low>=4:
            continue
        breaks=sorted(set([low,4.]+[v for v in (.02,.05,.1,.2,.5,1.,1.44,2.,2.5,3.) if low<v<4]))
        omega,wo=segments(breaks,order,square=True)
        inner={p:np.zeros((len(omega),3)) for p in PROFILES} # net,gross,electronic second moment
        if channel=='scattering':
            maximum=np.sqrt((emax-omega)**2-a*a) if g==0 and a>0 else emax-omega
            blocks=[]
            for lo,hi in zip(energy_breaks[:-1],energy_breaks[1:]):
                upper=np.minimum(maximum,hi)
                width=np.maximum(0,upper-lo)
                coordinate=lo+width[:,None]*(z[None,:]+1)/2
                measure=width[:,None]*w[None,:]/2
                blocks.append((coordinate,measure))
        else:
            if g==0 and a>0:
                theta=(z+1)*np.pi/4
                blocks=[(np.broadcast_to(theta,(len(omega),order)),np.broadcast_to(w*np.pi/4,(len(omega),order)))]
            else:
                # E in[0,Omega] since Omega_max4 is below either energy cutoff.
                blocks=[(omega[:,None]*(z[None,:]+1)/2,omega[:,None]*w[None,:]/2)]
        for coordinate,measure in blocks:
            oo=omega[:,None]
            if channel=='scattering':
                ei=np.hypot(a,coordinate) if g==0 and a>0 else coordinate
                ej=ei+oo
                if g==0 and a>0:
                    jacobian=(ej-a*a/ei)/np.sqrt((ej-a)*(ej+a))
                elif a==0:
                    jacobian=np.ones_like(ei)
                else:
                    ni,ri,err1=gapless_spectrum(ei,a,g)
                    nj,rj,err2=gapless_spectrum(ej,a,g)
                    spectral_error=max(spectral_error,err1,err2)
                    jacobian=ni*nj-ri*rj
                electron_delta=ei*ei-ej*ej
                factor=1.
            else:
                if g==0 and a>0:
                    length=oo-2*a
                    ei=a+length*np.sin(coordinate)**2
                    ej=a+length*np.cos(coordinate)**2
                    jacobian=2*(ei*ej+a*a)/np.sqrt((ei+a)*(ej+a))
                else:
                    ei=coordinate;ej=oo-ei
                    if a==0:
                        jacobian=np.ones_like(ei)
                    else:
                        ni,ri,err1=gapless_spectrum(ei,a,g)
                        nj,rj,err2=gapless_spectrum(ej,a,g)
                        spectral_error=max(spectral_error,err1,err2)
                        jacobian=ni*nj+ri*rj
                electron_delta=-ei*ei-ej*ej
                factor=.5 # D.15 relative to scattering prefactor, full E in[0,Omega]
            weighted=measure*jacobian*factor
            for profile in PROFILES:
                pi,n=populations(profile,ei,oo)
                pj,_=populations(profile,ej,oo)
                if channel=='scattering':
                    forward=pj*(1-pi)*(1+n);reverse=pi*(1-pj)*n
                else:
                    forward=pi*pj*(1+n);reverse=(1-pi)*(1-pj)*n
                net=forward-reverse
                inner[profile][:,0]+=np.sum(weighted*net,axis=1)
                inner[profile][:,1]+=np.sum(weighted*(forward+reverse),axis=1)
                inner[profile][:,2]+=np.sum(weighted*net*electron_delta,axis=1)
        coupling=.03*(omega/4)**2
        for profile in PROFILES:
            net,gross,second=(inner[profile][:,i]*coupling for i in range(3))
            result[profile][channel]=np.array([np.dot(wo*omega,net),np.dot(wo*omega,gross),
                                               np.dot(wo,net),np.dot(wo,gross),np.dot(wo,second)])
            projected[profile][channel]=project_hats(omega,net,wo,np.linspace(0,4,17))
    return result,projected,spectral_error,emax


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--orders',default='16,32,64')
    parser.add_argument('--output',type=Path,default=OUT/'continuous_reactions.json')
    args=parser.parse_args();orders=[int(v) for v in args.orders.split(',')]
    started=time.perf_counter();rows=[]
    for a,g in FIELDS:
        refinements=[]
        for order in orders:
            tick=time.perf_counter()
            values,projection,residual,emax=evaluate_field(a,g,order)
            refinements.append((values,projection))
            print(json.dumps(dict(field=[a,g],order=order,seconds=time.perf_counter()-tick,causal_residual=residual)),flush=True)
        for profile in PROFILES:
            channels={}
            for channel in ('scattering','recombination'):
                final=refinements[-1][0][profile][channel]
                previous=refinements[-2][0][profile][channel]
                error=abs(final-previous)
                channels[channel]=dict(reference=final.tolist(),refinement_absolute_change=error.tolist(),
                    refinement_relative_change=(error/np.maximum(abs(final),1e-100)).tolist(),
                    all_orders=[values[profile][channel].tolist() for values,_ in refinements],
                    phonon_rate_projected_17=refinements[-1][1][profile][channel].tolist())
            rows.append(dict(amplitude=a,gamma=g,profile=profile,energy_upper_at_count12=emax,channels=channels))
    output=dict(schema='pysnspd.stage2.independent_continuous_reactions.v1',
                status='REFERENCE_READY_REQUIRES_SEPARATE_FINITE_ETA_CONTROL',
                method='B.5-B.9 in(E,Omega), analytic normal/BCS DOS and singularity-cancelling substitutions; real-angle ideal gapless spectrum. No event enumerator or experimental operator imports.',
                eta=0.,rate_prefactor=1.,omega_support=[.01,4.],alpha2F='.03*(Omega/4)^2 on the declared support',
                count_cutoff=12.,orders=orders,observable_order=OBSERVABLES,
                profile_definitions=dict(hot_electrons='p=FD(E,.4); n=BE(Omega,.15)',
                    phonon_bubble='p=FD(E,.12); n=BE(Omega,.12)+.3 exp[-((Omega-2)/.35)^2]',
                    nonthermal='p=.12 exp(-E/.25)+.1 exp[-((E-1.3)/.25)^2]; n=BE(Omega,.18)+.1 exp[-((Omega-1.7)/.4)^2]'),
                criteria_sha256=sha(ROOT/'docs/implementation/stage2/acceptance_criteria.json'),
                reviewer_source_sha256=sha(Path(__file__)),
                prior_real_angle_source_sha256=sha(ROOT/'docs/implementation/stage1_r2/review/causal_reference.py'),
                rows=rows,runtime_seconds=time.perf_counter()-started)
    args.output.write_text(json.dumps(output,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(dict(output=str(args.output),runtime_seconds=output['runtime_seconds'],cases=len(rows))),flush=True)


if __name__=='__main__':
    main()
