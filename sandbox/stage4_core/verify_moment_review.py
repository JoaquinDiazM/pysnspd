"""Verify the charge-moment review and its preserved delivery chain."""
from pathlib import Path
import json,hashlib,runpy,sys
ROOT=Path(__file__).resolve().parents[2]
DATA=ROOT/'docs/implementation/stage4/moment_review_20260924'
PRIOR=ROOT/'docs/implementation/stage4/self_consistent_review_20260924'
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
    spatial=ROOT/'docs/implementation/stage4/spatial_energy_20260924'
    followup=records(read(ROOT/'docs/implementation/stage4/followup_20260923/delivery_manifest.json')['files'],spatial/'previous_delivery_exact')
    earlier=records(read(spatial/'delivery_manifest.json')['files'],PRIOR/'previous_delivery_exact')
    previous=records(read(PRIOR/'delivery_manifest.json')['files'],DATA/'previous_delivery_exact')
    current=read(DATA/'delivery_manifest.json');count=records(current['files'])
    assert not current['stage4_complete'] and not current['production_changed']
    raw=DATA/'raw';identity=read(raw/'identity.json');summary=read(raw/'summary.json')
    assert summary['status']=='MOMENT_RESOLUTION_CAMPAIGN_COMPLETE'
    assert summary['spectral_queries']==summary['kinetic_queries']==181 and summary['projection_groups']==7
    assert identity['plan_sha256']==sha(raw/'executed_plan.json')==sha(PRIOR/'resolution_campaign/plan.json')
    for name,digest in {**identity['sources'],**identity['inputs']}.items():assert sha(ROOT/name)==digest,name
    receipt=read(raw/'extraction_receipt.json');analysis=read(DATA/'analysis.json')
    assert receipt['status']=='ALL_RAW_MAPS_VERIFIED_NO_SOLVES'
    assert receipt['new_spectral_solves']==receipt['physical_time_steps']==0
    assert len(receipt['maps_verified'])==369 and receipt['compact_projection_maps']==7
    assert receipt['script_sha256']==sha(ROOT/'sandbox/stage4_core/extract_moment_results.py')
    for name,digest in receipt['copied_file_sha256'].items():assert sha(raw/name)==digest,name
    p=analysis['provenance']
    assert p['analysis_source_sha256']==sha(ROOT/'sandbox/stage4_core/analyze_moment_results.py')
    assert p['extraction_receipt_sha256']==sha(raw/'extraction_receipt.json')
    assert p['verified_raw_fields']==369 and p['moment_recomputation_maximum_difference']<1e-12
    frequency=DATA/'charge_frequency'; fid=read(frequency/'identity.json')
    assert fid['plan_sha256']==sha(frequency/'plan.json')
    for name,digest in fid['sources'].items():assert sha(ROOT/name)==digest,name
    for name,digest in fid['inputs'].items():
        if not name.startswith('spectra:'):assert sha(ROOT/name)==digest,name
    fs=read(frequency/'summary.json');fa=read(frequency/'analysis.json')
    assert fs['status']=='FIXED_SPECTRUM_ADIABATIC_FREQUENCY_DIAGNOSTIC_COMPLETE'
    assert len(fs['records'])==fa['jobs']==24
    assert max(item['maximum_absolute_difference'] for item in fa['static_overlap'])<1e-12
    compact=read(frequency/'compact_identity.json')
    assert compact['analysis_source_sha256']==sha(ROOT/'sandbox/stage4_core/charge_frequency_analysis.py')
    assert compact['summary_sha256']==sha(frequency/'summary.json')
    assert compact['comparisons_sha256']==sha(frequency/'comparisons.json')
    assert compact['compact_sha256']==sha(frequency/'representative_angular_phasors.npz')
    tests=read(frequency/'test_receipt.json')
    assert tests['passed'] and tests['tests']==61 and tests['failures']==tests['errors']==0
    for name,digest in tests['test_sources'].items():assert sha(ROOT/name)==digest,name
    independent=read(DATA/'frequency_review.json')['provenance']
    assert independent['source_identity_sha256']==sha(frequency/'identity.json')
    assert independent['postprocess_analysis_sha256']==sha(frequency/'analysis.json')
    assert independent['local_compact_field_recomputation_maximum_difference']<1e-12
    weak=DATA/'thermal_weak';wid=read(weak/'raw/identity.json')
    assert wid['plan_sha256']==sha(weak/'plan.json')
    for name,digest in {**wid['sources'],**wid['inputs']}.items():assert sha(ROOT/name)==digest,name
    ws=read(weak/'raw/summary.json');wr=read(weak/'verification_receipt.json')
    assert ws['status']=='FULL_NODE_THERMAL_OPERATOR_VALIDATION_COMPLETE' and ws['all_calculus_flags_met']
    assert len(ws['records'])==wr['checkpoint_hashes_verified']==256
    assert sum(row['nonlinear_roots'] for row in ws['records'])==wr['nonlinear_roots']==2304
    assert not wr['bad_hashes'] and wr['physical_time_steps']==0
    for name,digest in wr['local_files'].items():assert sha(weak/name)==digest,name
    wt=read(weak/'unit_tests_receipt.json')
    assert wt['status']=='PASSED' and wt['tests']==5 and wt['failures']==wt['errors']==0
    for name,digest in wt['sources'].items():assert sha(ROOT/name)==digest,name
    assert sha(weak/wt['log_path'])==wt['log_sha256']
    temporal=DATA/'thermal_time';tp=read(temporal/'plan.json')
    assert tp['operator_summary_sha256']==sha(weak/'raw/summary.json')
    assert tp['operator_checks_sha256']==sha(weak/'raw/full_node_operator_checks.npz')
    assert tp['operator_identity_sha256']==sha(weak/'raw/identity.json')
    assert tp['observation_times_ps'][-1]==1.0 and tp['maximum_workers']<=27
    assert tp['link_field']=='alpha_zero' and tp['perturbation_amplitude']==.001
    tv=read(temporal/'validation_receipt.json')
    assert tv['status']=='PREPARED_NOT_PHYSICALLY_EXECUTED' and not tv['physical_transient_executed']
    for name,digest in tv['files'].items():assert sha(temporal/name)==digest,name
    tt=read(temporal/'unit_tests_receipt_v2.json')
    assert tt['passed'] and tt['tests']==6 and tt['failures']==tt['errors']==0
    for name,digest in tt['sources'].items():assert sha(ROOT/name)==digest,name
    dry=read(temporal/'dry_run_v2.json')
    assert dry['status']=='READY_DRY_RUN_NO_WRITES' and dry['plan_sha256']==sha(temporal/'plan.json')
    for name,digest in dry['sources'].items():assert sha(ROOT/name)==digest,name
    assert dry['budget']['total_processes']<=dry['budget']['logical_budget']
    qa=read(DATA/'report_qa.json')
    assert qa['status']=='PASS' and qa['pdf_sha256']==sha(ROOT/'output/pdf/implementation/Informe_etapa_4_respuesta_de_carga.pdf')
    external=Path('/home/jdiaz/GEMINGA_COMMANDS.md')
    if sys.platform.startswith('linux') and external.exists():assert sha(external)==sha(ROOT/'docs/GEMINGA_COMMANDS.md')
    print(json.dumps(dict(verified=True,current_files=count,previous_files=previous,
        earlier_files=earlier,followup_files=followup,historical_chain_verified=oldest['verified'],
        completed_spectral_queries=181,completed_kinetic_queries=181,projection_maps=7,
        remote_raw_maps_verified=369,charge_frequency_jobs=24,charge_frequency_tests=61,
        thermal_operator_frequencies=256,thermal_operator_roots=2304,thermal_operator_tests=5,
        thermal_time_horizon_ps=1.,thermal_time_prepared=True,thermal_time_completed=False,
        temporal_driver_tests=6,
        stage4_complete=False,production_changed=False),indent=2))
if __name__=='__main__':main()
