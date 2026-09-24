"""Read-only verification of the short ETD2 execution and its saved fields."""
import hashlib,json
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[2]
DATA=ROOT/'docs/implementation/stage4/time_review_20260924/nonlinear_time'
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def read(path):return json.loads(path.read_text(encoding='utf8'))


def certificate():
    raw=DATA/'pilot_raw';identity=read(raw/'identity.json');summary=read(raw/'summary.json')
    assert summary['status']=='NONLINEAR_THERMAL_TIME_PILOT_COMPLETE'
    assert identity['pilot'] and identity['horizon_ps']==summary['horizon_ps']==.001
    assert not summary['stage4_complete'] and not summary['production_changed']
    assert identity['plan_sha256']==sha(raw/'executed_plan.json')==sha(DATA/'plan.json')
    for name,digest in identity['sources'].items():assert sha(ROOT/name)==digest,name
    fields=[];maximum_balance=0.;histories=[]
    for label in ('primary','refined'):
        result=read(raw/label/'summary.json');histories.append(result['history'])
        assert [item['time_ps'] for item in result['observations']]==[0.,.001]
        for record in result['observations']:
            path=raw/record['fields_path'];assert sha(path)==record['fields_sha256']
            with np.load(path) as a:
                assert all(np.all(np.isfinite(a[k])) for k in a.files)
                fixed=a['boundary_nodes']
                for name in ('baseline','amplitude','angular_phase'):
                    assert np.all(a[name+'_displacement'][fixed]==0.)
                    assert np.all(a[name+'_velocity'][fixed]==0.)
            fields.append(dict(path=path.relative_to(raw).as_posix(),sha256=sha(path)))
            maximum_balance=max(maximum_balance,max(abs(row['power_residual']) for row in record['states']))
        for record in result['history']:
            if not record['accepted']:continue
            path=raw/label/record['fields_path'];assert sha(path)==record['fields_sha256']
            fields.append(dict(path=path.relative_to(raw).as_posix(),sha256=sha(path)))
    assert summary['refinement']['all_comparisons_met']
    budget=identity['budget'];assert budget['workers']<=27
    return dict(status='PILOT_INTEGRITY_PASSED',horizon_ps=.001,runtime_seconds=summary['runtime_seconds'],
        roots=summary['spectral_roots'],maximum_spectral_residual=summary['maximum_spectral_residual'],
        maximum_power_residual=maximum_balance,checked_fields=fields,
        same_accepted_step_pattern=[(row['start_ps'],row['step_ps']) for row in histories[0] if row['accepted']]==[(row['start_ps'],row['step_ps']) for row in histories[1] if row['accepted']],
        interpretation='Short execution verifies the assembled pipeline and fixed contacts. Both passes take one identical step; this is not independent temporal convergence of the 1ps trajectory.',
        long_trajectory_executed=False,timeout_seconds=180,exit_code=0,
        source_sha256=sha(Path(__file__)),files={p.relative_to(raw).as_posix():sha(p) for p in raw.rglob('*') if p.is_file()})


if __name__=='__main__':print(json.dumps(certificate(),indent=2))
