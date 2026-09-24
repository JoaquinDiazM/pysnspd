"""Read-only verification of the thermal spatial delivery and its predecessors."""
from pathlib import Path
import hashlib
import json
import runpy
import sys

ROOT=Path(__file__).resolve().parents[2]
DATA=ROOT/'docs/implementation/stage4/spatial_energy_20260924'
ARCHIVE=DATA/'previous_delivery_exact'


def read(path):
    return json.loads(Path(path).read_text(encoding='utf8'))


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def records(items, archive=None):
    for item in items:
        path=ROOT/item['path']
        if archive is not None and (archive/item['path']).is_file():path=archive/item['path']
        if sha(path)!=item['sha256'] or path.stat().st_size!=item['bytes']:
            raise RuntimeError('Changed delivery file: '+str(path))
    return len(items)


def main():
    predecessor=runpy.run_path(str(ROOT/'sandbox/stage4_core/verify_followup.py'))['predecessor']()
    previous=read(ROOT/'docs/implementation/stage4/followup_20260923/delivery_manifest.json')
    previous_count=records(previous['files'],ARCHIVE)
    current=read(DATA/'delivery_manifest.json')
    current_count=records(current['files'])
    raw=DATA/'raw/stage4_spatial_energy_reference_20260924'
    identity=read(raw/'identity.json')
    assert identity['plan_sha256']==sha(DATA/'campaign_plan.json')
    for name,digest in identity['sources'].items():assert sha(ROOT/name)==digest,name
    summary=read(raw/'summary.json');count=0
    for case in summary['cases'].values():
        for values in case['counts'].values():
            assert sha(raw/values['fields_file'])==values['fields_sha256'];count+=1
    certificate=read(DATA/'remote_mode_integrity.json')
    assert certificate['summary_sha256']==sha(raw/'summary.json') and certificate['records']==len(summary['records'])==2048
    pilot=DATA/'self_consistent_pilot'
    pilot_identity=read(pilot/'identity.json')
    assert pilot_identity['plan_sha256']==sha(DATA/'self_consistent_plan.json')
    for name,digest in pilot_identity['sources'].items():assert sha(ROOT/name)==digest,name
    pilot_summary=read(pilot/'summary.json')
    pilot_receipt=read(DATA/'self_consistent_pilot_receipt.json')
    assert pilot_receipt['summary_sha256']==sha(pilot/'summary.json')
    assert pilot_summary['completed_jobs']==1792 and not pilot_summary['all_core_criteria_met']
    resume_certificate=read(DATA/'remote_checkpoint_integrity.json')
    assert resume_certificate['summary_sha256']==sha(pilot/'summary.json')
    assert resume_certificate['checkpoints_verified']==8
    for name,case in pilot_summary['cases'].items():
        assert case['fields_sha256']==pilot_receipt['cases'][name]['checkpoint_sha256']
    tests=read(DATA/'test_receipt.json')
    assert tests['passed'] and tests['tests']==36
    for name,digest in tests['test_sources'].items():assert sha(ROOT/name)==digest,name
    assert not current['production_changed'] and not current['stage4_complete']
    external=Path('/home/jdiaz/GEMINGA_COMMANDS.md')
    if sys.platform.startswith('linux') and external.exists():assert sha(external)==sha(ROOT/'docs/GEMINGA_COMMANDS.md')
    print(json.dumps(dict(verified=True,current_files=current_count,previous_files=previous_count,
        previous=predecessor,spatial_maps=count,remote_mode_records=2048,self_consistent_pilot_jobs=1792,
        remote_pilot_checkpoints=8,focused_tests=36,
        mode_files_scope='Remote integrity receipt; raw modes remain on Geminga, not rehashed locally',
        stage4_complete=False,production_changed=False),indent=2))


if __name__=='__main__':main()
