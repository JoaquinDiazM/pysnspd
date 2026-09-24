"""Verify the published core review, source identities and preserved history."""
from pathlib import Path
import json,hashlib,runpy,sys

ROOT=Path(__file__).resolve().parents[2]
DATA=ROOT/'docs/implementation/stage4/self_consistent_review_20260924'
PREVIOUS=ROOT/'docs/implementation/stage4/spatial_energy_20260924'
def read(p):return json.loads(Path(p).read_text(encoding='utf8'))
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def records(items,archive=None):
    for item in items:
        p=ROOT/item['path']
        if archive is not None and (archive/item['path']).is_file():p=archive/item['path']
        assert sha(p)==item['sha256'] and p.stat().st_size==item['bytes'],str(p)
    return len(items)
def main():
    oldest=runpy.run_path(str(ROOT/'sandbox/stage4_core/verify_followup.py'))['predecessor']()
    followup=records(read(ROOT/'docs/implementation/stage4/followup_20260923/delivery_manifest.json')['files'],PREVIOUS/'previous_delivery_exact')
    previous=records(read(PREVIOUS/'delivery_manifest.json')['files'],DATA/'previous_delivery_exact')
    current=read(DATA/'delivery_manifest.json');count=records(current['files'])
    assert not current['stage4_complete'] and not current['production_changed']
    raw=DATA/'raw';identity=read(raw/'identity.json');summary=read(raw/'summary.json')
    assert summary['all_core_criteria_met'] and summary['status']=='FINITE_SUM_CORE_STATIONARY'
    assert len(summary['cases'])==4 and summary['physical_time_steps']==0
    assert identity['plan_sha256']==sha(raw/'executed_plan.json')==sha(PREVIOUS/'self_consistent_plan.json')
    for name,digest in identity['sources'].items():assert sha(ROOT/name)==digest,name
    receipt=read(raw/'extraction_receipt.json');analysis=read(DATA/'analysis.json')
    assert receipt['summary_sha256']==sha(raw/'summary.json') and receipt['identity_sha256']==sha(raw/'identity.json')
    assert receipt['extraction_script_sha256']==sha(ROOT/'sandbox/stage4_core/extract_self_consistent_results.py')
    assert receipt['new_spectral_solves']==0 and receipt['completed_sweep_records']==64
    assert len(receipt['checkpoint_files_verified'])==23
    assert analysis['provenance']['analysis_script_sha256']==sha(ROOT/'sandbox/stage4_core/analyze_self_consistent_results.py')
    assert analysis['provenance']['extraction_receipt_sha256']==sha(raw/'extraction_receipt.json')
    for name,case in receipt['final_cases'].items():
        assert sha(raw/case['compact_file'])==case['compact_sha256']
        assert case['full_checkpoint_sha256']==summary['cases'][name]['fields_sha256']
    for name,digest in analysis['provenance']['history_sha256'].items():assert sha(raw/'histories'/name)==digest,name
    oracle=DATA/'retarded_oracle';spectral_identity=read(oracle/'identity.json')
    spectral_summary=read(oracle/'summary.json');spectral_review=read(oracle/'analysis.json')
    assert spectral_summary['status']=='RETARDED_ORACLE_COMPLETE' and len(spectral_summary['records'])==72
    assert spectral_review['jobs']==72 and spectral_review['continuation_steps']==1728
    assert spectral_identity['plan_sha256']==sha(DATA/'next_retarded_plan.json')
    for name,digest in {**spectral_identity['sources'],**spectral_identity['inputs']}.items():assert sha(ROOT/name)==digest,name
    tests=read(DATA/'test_receipt.json');assert tests['passed'] and tests['tests']==44
    for name,digest in tests['test_sources'].items():assert sha(ROOT/name)==digest,name
    kinetic=DATA/'kinetic_response';ki=read(kinetic/'identity.json');ks=read(kinetic/'summary.json')
    assert ks['status']=='FROZEN_KINETIC_RESPONSE_COMPLETE' and len(ks['records'])==72
    assert ki['spectral_summary_sha256']==sha(oracle/'summary.json')
    assert ki['plan_sha256']==sha(DATA/'next_kinetic_plan.json')
    for name,digest in {**ki['sources'],**ki['inputs']}.items():assert sha(ROOT/name)==digest,name
    ka=read(kinetic/'analysis.json')
    assert ka['mathematical_probe_responses']==144 and ka['integrated_checks']==12
    assert ka['verified_raw_files']==84 and not ka['stage4_complete']
    tests2=read(DATA/'kinetic_test_receipt.json');assert tests2['passed'] and tests2['tests']==52
    for name,digest in tests2['test_sources'].items():assert sha(ROOT/name)==digest,name
    projection=DATA/'potential_projection';pr=read(projection/'verification_receipt.json')
    ps=read(projection/'summary.json')
    assert pr['identity_sha256']==sha(projection/'identity.json')
    assert pr['summary_sha256']==sha(projection/'summary.json')
    assert len(pr['verified_maps'])==6 and pr['new_spectral_solves']==0
    assert ps['new_spectral_solves']==0 and ps['physical_time_steps']==0 and not ps['chi_renormalized']
    assert not ps['exact_memory_ohmic_law_admitted'] and not ps['instantaneous_charge_dynamics_admitted']
    for name,digest in pr['sources'].items():assert sha(ROOT/name)==digest,name
    for name,digest in pr['verified_maps'].items():assert sha(projection/name)==digest,name
    campaign=DATA/'resolution_campaign';dry=read(campaign/'dry_run.json')
    assert dry['status']=='DRY_RUN_NO_WRITES' and dry['spectral_queries']==181
    assert dry['kinetic_queries']==181 and dry['projection_groups']==7
    assert dry['plan_sha256']==sha(campaign/'plan.json')
    for name,digest in {**dry['sources'],**dry['inputs']}.items():assert sha(ROOT/name)==digest,name
    ct=read(campaign/'test_receipt.json');assert ct['passed'] and ct['unit_tests']==4
    for name,digest in ct['sources'].items():assert sha(ROOT/name)==digest,name
    pilot=read(campaign/'stage4_resolution_eta001_pilot_20260924/summary.json')
    assert pilot['status']=='HARDEST_QUERY_PILOT_COMPLETE' and pilot['spectral_queries']==1
    qa=read(DATA/'report_qa.json')
    assert qa['status']=='PASS' and qa['pages']==4
    assert qa['pdf_sha256']==sha(ROOT/'output/pdf/implementation/Informe_etapa_4_nucleo_autoconsistente.pdf')
    external=Path('/home/jdiaz/GEMINGA_COMMANDS.md')
    if sys.platform.startswith('linux') and external.exists():assert sha(external)==sha(ROOT/'docs/GEMINGA_COMMANDS.md')
    print(json.dumps(dict(verified=True,current_files=count,previous_files=previous,followup_files=followup,
        historical_chain_verified=oldest['verified'],accepted_thermal_cases=4,sweep_records=64,
        remote_checkpoints_certified=23,retarded_queries=72,focused_tests=52,
        frozen_kinetic_responses=144,coarse_potential_projections=6,
        next_campaign_preflight_queries=181,campaign_focused_tests=4,long_campaign_executed=False,
        full_checkpoint_scope='Recorded remote hashes, not full spectral arrays in Git',
        stage4_complete=False,production_changed=False),indent=2))
if __name__=='__main__':main()
