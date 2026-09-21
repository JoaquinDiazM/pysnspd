"""Refine the independent phonon measure on 17 and 33 common nodal hats.

Each Omega hat knot is a quadrature break. This removes a projection-specific
kink error that a converged total power alone would not detect.
"""
from pathlib import Path
import json
import time
import numpy as np
import continuous_reactions as reference


def main():
    started=time.perf_counter();out=Path(__file__).resolve().parent
    original_segments=reference.segments;original_project=reference.project_hats
    previous=json.loads((out/'continuous_reactions.json').read_text())
    rows=[]
    for size in (17,33):
        axis=np.linspace(0,4,size)
        def aligned_segments(bounds,order,**kwargs):
            cuts=sorted(set(list(bounds)+[float(v) for v in axis if bounds[0]<v<bounds[-1]]))
            return original_segments(cuts,order,**kwargs)
        reference.segments=aligned_segments
        reference.project_hats=lambda omega,rate,weights,unused:original_project(omega,rate,weights,axis)
        # The historical helper allocates projected vectors dynamically on
        # assignment, so both bases are valid without changing its source.
        for a,g in reference.FIELDS:
            coarse=reference.evaluate_field(a,g,32)
            fine=reference.evaluate_field(a,g,64)
            for profile in reference.PROFILES:
                base=next(v for v in previous['rows'] if v['amplitude']==a and v['gamma']==g and v['profile']==profile)
                channels={}
                for channel in ('scattering','recombination'):
                    v=fine[1][profile][channel];old=coarse[1][profile][channel]
                    scale=max(float(np.sum(abs(v))),1e-100)
                    item=dict(reference=v.tolist(),reference_L1=float(np.sum(abs(v))),
                              last_refinement_relative_L1=float(np.sum(abs(v-old))/scale))
                    if size==17:
                        prev=np.asarray(base['channels'][channel]['phonon_rate_projected_17'])
                        item['previous_unaligned_reference_relative_L1']=float(np.sum(abs(v-prev))/scale)
                    channels[channel]=item
                rows.append(dict(amplitude=a,gamma=g,profile=profile,hat_nodes=size,channels=channels))
    data=dict(schema='pysnspd.stage2.continuous_projection_reference.v1',
        status='IDEAL_ETA_ZERO_REFERENCE',
        method='Same independent (E,Omega) integrals with every common-hat knot an Omega quadrature boundary',
        hat_axes={str(size):np.linspace(0,4,size).tolist() for size in (17,33)},
        orders=[32,64],eta=0.,rows=rows,
        criteria_sha256=reference.sha(reference.ROOT/'docs/implementation/stage2/acceptance_criteria.json'),
        historical_reference_sha256=reference.sha(out/'continuous_reactions.json'),
        source_sha256=reference.sha(Path(__file__)),
        runtime_seconds=time.perf_counter()-started)
    (out/'continuous_projection_check.json').write_text(json.dumps(data,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(dict(runtime_seconds=data['runtime_seconds'],
        worst_refinement=max(ch['last_refinement_relative_L1'] for r in rows for ch in r['channels'].values()),
        worst_historical_change=max(ch.get('previous_unaligned_reference_relative_L1',0) for r in rows for ch in r['channels'].values()))))


if __name__=='__main__':main()
