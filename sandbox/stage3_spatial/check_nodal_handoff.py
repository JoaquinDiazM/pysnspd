"""Read-only preflight for the user-run nodal campaign."""
import json
from run_nodal_batch import ROOT,frozen_sources,verify_reuse


def main():
    reg=ROOT/'docs/implementation/stage3/nodal_20260923/registration.json'
    value=json.loads(reg.read_text(encoding='utf-8'))
    case=next(c for c in value['cases'] if c['id']=='weak_phase_thermal')
    pilot=verify_reuse(ROOT/'tmp/stage3_nodal_20260923/pilot',
                      'weak_phase_thermal_n8',frozen_sources(reg),case,8)
    if pilot is None:raise RuntimeError('Nodal pilot absent')
    if (ROOT/'tmp/stage3a_nodal_full_20260923').exists():
        raise RuntimeError('Output already exists; inspect before relaunching')
    print(json.dumps(dict(status='READY_FOR_USER_EXECUTION',reusable_pilot=True,
                         output_is_fresh=True,calculations_launched=False)))


if __name__=='__main__':main()
