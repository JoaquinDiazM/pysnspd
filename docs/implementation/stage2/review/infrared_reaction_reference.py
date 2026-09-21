"""Independent omitted [0,Omega_min] reaction integrals, without cancellation."""
from pathlib import Path
import json
import time
import numpy as np
import continuous_reactions as ref


def main():
    started=time.perf_counter();out=Path(__file__).resolve().parent
    mainref=json.loads((out/'continuous_reactions.json').read_text())
    original_segments=ref.segments
    omitted={};convergence={}
    for cutoff in (.02,.01,.005):
        def restricted_segments(bounds,order,**kwargs):
            low=bounds[0]
            if low>=cutoff:return np.array([]),np.array([])
            return original_segments([low,cutoff],order,**kwargs)
        ref.segments=restricted_segments
        for a,g in ref.FIELDS:
            prev=ref.evaluate_field(a,g,32,omega_min=0.)[0]
            final=ref.evaluate_field(a,g,64,omega_min=0.)[0]
            omitted[(cutoff,a,g)]=final
            convergence[(cutoff,a,g)]=prev
    rows=[]
    for a,g in ref.FIELDS:
        for profile in ref.PROFILES:
            old=next(r for r in mainref['rows'] if r['amplitude']==a and r['gamma']==g and r['profile']==profile)
            for cutoff in (.02,.01,.005):
                channels={}
                for channel in ('scattering','recombination'):
                    full=np.asarray(old['channels'][channel]['reference'])+omitted[(.01,a,g)][profile][channel]
                    removal=omitted[(cutoff,a,g)][profile][channel]
                    previous=convergence[(cutoff,a,g)][profile][channel]
                    denominator=np.maximum(abs(full),1e-100)
                    channels[channel]=dict(full_uncut_reference=full.tolist(),omitted=removal.tolist(),
                        retained=(full-removal).tolist(),
                        omitted_relative_to_own_observable=(abs(removal)/denominator).tolist(),
                        omitted_reference_refinement_relative_to_full=(abs(removal-previous)/denominator).tolist(),
                        omitted_net_power_over_full_gross_power=float(abs(removal[0])/max(full[1],1e-100)),
                        omitted_net_count_over_full_gross_count=float(abs(removal[2])/max(full[3],1e-100)))
                rows.append(dict(amplitude=a,gamma=g,profile=profile,infrared_cutoff=cutoff,channels=channels))
    data=dict(schema='pysnspd.stage2.independent_infrared_reactions.v1',
        status='INDEPENDENT_IDEAL_CONTINUUM_CUTOFF_CONTROL',
        method='Integrate omitted0..Omega_min directly in(E,Omega), avoiding subtraction of nearly equal rates. Add0..0.01 to the independently converged historical[0.01,4] reference.',
        eta=0.,rate_prefactor=1.,alpha2F='.03*(Omega/4)^2, extended to Omega0; Debye uppercut4',
        orders=[32,64],observable_order=ref.OBSERVABLES,rows=rows,
        criteria_sha256=ref.sha(ref.ROOT/'docs/implementation/stage2/acceptance_criteria.json'),
        parent_reference_sha256=ref.sha(out/'continuous_reactions.json'),
        source_sha256=ref.sha(Path(__file__)),runtime_seconds=time.perf_counter()-started)
    (out/'infrared_reaction_reference.json').write_text(json.dumps(data,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(dict(runtime_seconds=data['runtime_seconds'],
        worst_reference_refinement=max(max(ch['omitted_reference_refinement_relative_to_full']) for r in rows for ch in r['channels'].values()),
        worst_omitted_by_cutoff={str(c):max(max(ch['omitted_relative_to_own_observable'][:4]) for r in rows if r['infrared_cutoff']==c for ch in r['channels'].values()) for c in (.02,.01,.005)})))


if __name__=='__main__':main()
