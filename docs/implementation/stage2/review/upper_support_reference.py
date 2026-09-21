"""Independent bound on electron absorption leaving the represented support.

This bounds missing inside-to-outside scattering without assuming the unknown
external hole probability: 1-p_out <= 1. It cannot bound emission from arbitrary
unrepresented electrons. The analytic-profile emission below is a separately
labelled tail assumption, never a universal closure of that missing population.
"""
from pathlib import Path
import json
import time

import numpy as np

from continuous_reactions import (FIELDS, PROFILES, ROOT, gapless_spectrum,
                                  populations, segments, sha)


def upper_absorption_bound(amplitude, gamma, energy_upper, electron_population,
                           phonon_population, *, order=32, omega_min=.01,
                           omega_max=4., coupling=None,
                           external_population=None):
    """Return rates with reference prefactor 8*pi*Delta0*tref/hbar = 1.

    ``electron_population`` is only called INSIDE its represented support.
    Callers must supply a physical population on that interval, including its
    endpoint; this helper does not extrapolate their array. Gamma=0 or gamma>a
    ideal spectra are supported. Finite-eta error remains a separate control.
    """
    a=float(amplitude);g=float(gamma);emax=float(energy_upper)
    gap=a if g==0 else 0.
    if not (a>=0 and g>=0 and emax>omega_max+gap and 0<omega_min<omega_max):
        raise ValueError('unsupported field/cutoff; lower tail must remain above the gap')
    if g!=0 and g<=a:
        raise ValueError('this independent helper only covers normal, BCS and gapless fields')
    cuts=sorted(set([omega_min,omega_max]+[v for v in (.02,.05,.1,.2,.5,1.,1.5,2.,2.5,3.) if omega_min<v<omega_max]))
    om,wo=segments(cuts,order,square=True)
    z,w=np.polynomial.legendre.leggauss(order)
    inside=emax-om[:,None]+om[:,None]*(z[None,:]+1)/2
    outside=inside+om[:,None]
    measure=om[:,None]*w[None,:]/2
    if a==0:
        coherence=np.ones_like(inside)
    elif g==0:
        coherence=(inside*outside-a*a)/np.sqrt((inside*inside-a*a)*(outside*outside-a*a))
    else:
        ni,ri,_=gapless_spectrum(inside,a,g)
        no,ro,_=gapless_spectrum(outside,a,g)
        coherence=ni*no-ri*ro
    p=np.asarray(electron_population(inside),float)
    n=np.asarray(phonon_population(om),float)
    if np.any(~np.isfinite(p)) or np.any((p<0)|(p>1)) or np.any(~np.isfinite(n)) or np.any(n<0):
        raise ValueError('invalid populations in the tail bound')
    alpha=.03*(om/4)**2 if coupling is None else np.asarray(coupling(om),float)
    if np.any(~np.isfinite(alpha)) or np.any(alpha<0):
        raise ValueError('invalid coupling')
    weighted=wo[:,None]*measure*coherence*alpha[:,None]
    reverse_bound=p*n[:,None]
    result={
        'absorption_number_rate_upper_bound':float(np.sum(weighted*reverse_bound)),
        'absorption_phonon_power_upper_bound':float(np.sum(weighted*reverse_bound*om[:,None])),
        'absorption_inside_electron_energy_loss_upper_bound':float(np.sum(weighted*reverse_bound*inside)),
        'absorption_external_electron_energy_gain_upper_bound':float(np.sum(weighted*reverse_bound*outside)),
        'recombination_outside_support':0.,
    }
    if external_population is not None:
        po=np.asarray(external_population(outside),float)
        if np.any(~np.isfinite(po)) or np.any((po<0)|(po>1)):
            raise ValueError('invalid separately assumed external population')
        emission=po*(1-p)*(1+n[:,None])
        result['assumed_profile_external_emission_phonon_power']=float(np.sum(weighted*emission*om[:,None]))
    return result


def main():
    started=time.perf_counter()
    path=Path(__file__).resolve().parent
    reference=json.loads((path/'continuous_reactions.json').read_text())
    rows=[]
    for a,g in FIELDS:
        for profile in PROFILES:
            ref=next(v for v in reference['rows'] if v['amplitude']==a and v['gamma']==g and v['profile']==profile)
            p=lambda e: populations(profile,e,1.)[0]
            n=lambda om: populations(profile,1.,om)[1]
            refinements=[upper_absorption_bound(a,g,ref['energy_upper_at_count12'],p,n,
                                                order=order,external_population=p)
                         for order in (16,32,64)]
            final=refinements[-1]
            power=sum(ch['reference'][1] for ch in ref['channels'].values())
            count=sum(ch['reference'][3] for ch in ref['channels'].values())
            rows.append(dict(amplitude=a,gamma=g,profile=profile,
                energy_upper=ref['energy_upper_at_count12'],all_orders=refinements,
                final=final,
                last_refinement_relative_change={k:abs(v-refinements[-2][k])/max(abs(v),1e-100) for k,v in final.items()},
                absorption_power_bound_over_in_support_gross_power=final['absorption_phonon_power_upper_bound']/power,
                absorption_count_bound_over_in_support_gross_count=final['absorption_number_rate_upper_bound']/count,
                assumed_external_emission_power_over_in_support_gross_power=final['assumed_profile_external_emission_phonon_power']/power))
    data=dict(schema='pysnspd.stage2.independent_upper_support.v1',
        status='REFERENCE_PROFILES_ONLY_ACTUAL_TRAJECTORIES_REQUIRE_THEIR_POPULATIONS',
        method='Independent (E,Omega) integral for E in [Emax-Omega,Emax], external hole <=1. No experimental module imports.',
        limitations=['eta=0 reference; finite-eta bias is separately controlled',
                     'No bound on arbitrary incoming electrons from unrepresented energies',
                     'Lower support and quadrature truncation are separate gates',
                     'External emission columns use the stated analytic continuation of each profile'],
        orders=[16,32,64],rate_prefactor=1.,omega_support=[.01,4.],
        source_sha256=sha(Path(__file__)),
        continuous_reference_sha256=sha(path/'continuous_reactions.json'),
        criteria_sha256=sha(ROOT/'docs/implementation/stage2/acceptance_criteria.json'),
        rows=rows,runtime_seconds=time.perf_counter()-started)
    (path/'upper_support_reference.json').write_text(json.dumps(data,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(dict(runtime_seconds=data['runtime_seconds'],cases=len(rows),
        max_absorption_power_relative=max(r['absorption_power_bound_over_in_support_gross_power'] for r in rows),
        max_external_emission_relative=max(r['assumed_external_emission_power_over_in_support_gross_power'] for r in rows))))


if __name__=='__main__':
    main()
