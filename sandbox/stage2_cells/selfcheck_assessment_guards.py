"""Lightweight negative controls for the resumed stage-2 assessment guards."""
from pathlib import Path
import argparse
import hashlib
import json
import time

import numpy as np

from assess_grid_refinement import self_test as grid_self_test
from assess_time_refinement import invariant_checks
from check_trajectory_fields import require_sample_count as require_field_samples
from check_trajectory_support import require_sample_count as require_support_samples

ROOT=Path(__file__).resolve().parents[2]


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    started=time.perf_counter()
    criteria_path=ROOT/'docs/implementation/stage2/acceptance_criteria.json'
    criteria=json.loads(criteria_path.read_text())
    minimum=criteria['catalogue_field_checks']['minimum_actual_trajectory_samples']
    record=dict(energy_ledger_scaled_max=0.,instantaneous_residual_max=0.,
                minimum_electron=0.,maximum_electron=1.,minimum_phonon=0.)
    checks={'valid_invariants_accepted':all(invariant_checks(record,criteria).values())}
    for key,limit_name in (('energy_ledger_scaled_max','energy_ledger_scaled_max'),
                           ('instantaneous_residual_max','instantaneous_balance_scaled_max')):
        bad=record|{key:2*criteria['coupled_trajectories'][limit_name]}
        checks[key+'_rejected_at_any_resolution']=not all(invariant_checks(bad,criteria).values())
    for key,value in (('minimum_electron',-1e-6),('maximum_electron',1.000001),
                      ('minimum_phonon',-1e-6)):
        checks[key+'_rejected']=not all(invariant_checks(record|{key:value},criteria).values())
    for label,guard in (('fields',require_field_samples),('support',require_support_samples)):
        result=guard(np.linspace(0,1,minimum),np.arange(minimum),2,minimum)
        checks[label+'_minimum_applied_per_cell']=result['distinct_physical_times_per_cell']==minimum
        for suffix,times,indices in (
                ('three_points',np.array([0.,.05,.1]),np.arange(3)),
                ('duplicate_physical_times',np.r_[np.zeros(minimum-1),1.],np.arange(minimum)),
                ('nonfinite_time',np.r_[np.arange(minimum-1,dtype=float),np.nan],np.arange(minimum))):
            try:
                guard(times,indices,2,minimum)
            except ValueError:
                checks[label+'_'+suffix+'_rejected']=True
            else:
                checks[label+'_'+suffix+'_rejected']=False
    grid=grid_self_test()
    checks.update({'grid_'+name:passed for name,passed in grid['checks'].items()})
    sources=[Path(__file__),*(ROOT/'sandbox/stage2_cells'/name for name in (
        'assess_time_refinement.py','assess_grid_refinement.py',
        'check_trajectory_fields.py','check_trajectory_support.py'))]
    result=dict(schema='pysnspd.stage2.resumed-assessment-selfcheck.v1',
        status='PASS' if all(checks.values()) else 'FAIL',checks=checks,
        passed_checks=sum(checks.values()),total_checks=len(checks),
        criteria_sha256=sha(criteria_path),
        sources={p.relative_to(ROOT).as_posix():sha(p) for p in sources},
        scope='Nonphysical fixtures and negative controls test adjudication only; no physical trajectory is admitted.',
        runtime_seconds=time.perf_counter()-started)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({key:result[key] for key in ('status','passed_checks','total_checks','runtime_seconds')}))
    if not all(checks.values()):raise SystemExit(1)


if __name__=='__main__':main()
