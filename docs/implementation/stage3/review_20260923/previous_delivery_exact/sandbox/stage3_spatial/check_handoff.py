"""Read-only preflight for the prepared Geminga command; no physics is run."""
from pathlib import Path
import json
from run_static_batch import ROOT, frozen_sources, verify_reuse


def main():
    registration_path=ROOT/'docs/implementation/stage3/iteration_20260923/registration.json'
    registration=json.loads(registration_path.read_text(encoding='utf-8'))
    case=next(c for c in registration['cases'] if c['id']=='weak_phase_thermal')
    sources=frozen_sources(registration_path)
    pilot=verify_reuse(ROOT/'tmp/stage3_start_20260923/pilot_reviewed',
                       'weak_phase_thermal_n8',sources,case,8)
    if pilot is None:
        raise RuntimeError('The reviewed pilot is absent')
    output=ROOT/'tmp/stage3a_static_20260923'
    if output.exists():
        raise RuntimeError('Output already exists; inspect it before starting another run')
    print(json.dumps(dict(status='READY_FOR_USER_EXECUTION',
                         reviewed_pilot_reusable=True,output_directory_is_fresh=True,
                         total_cases=len(registration['cases'])*len(registration['cell_counts']),
                         calculations_launched=False),indent=2))


if __name__=='__main__':
    main()
