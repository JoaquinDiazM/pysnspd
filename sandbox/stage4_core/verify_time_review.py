"""Read-only verification of the thermal-time delivery and preserved history.

No tests, spectral roots, trajectories or automatic retries are executed.
Remote-only arrays are checked against their extraction certificates locally;
when their original paths are mounted on Geminga their bytes are rehashed too.
"""
from pathlib import Path
import hashlib
import json
import runpy
import subprocess
import sys

ROOT=Path(__file__).resolve().parents[2]
DATA=ROOT/'docs/implementation/stage4/time_review_20260924'
MOMENT=ROOT/'docs/implementation/stage4/moment_review_20260924'
TAG_COMMIT='5ea0cd6a08d2cae414c82944f8230469ca5aa7d6'


def read(path):return json.loads(Path(path).read_text(encoding='utf8'))
def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def relative(name):return Path(name.replace('\\','/'))
def require(condition,message):
    if not condition:raise RuntimeError(message)


def records(items,archive=None):
    names=set()
    for row in items:
        name=relative(row['path']);require(not name.is_absolute() and '..' not in name.parts,'Unsafe manifest path')
        require(name.as_posix() not in names,'Duplicate manifest path: '+str(name));names.add(name.as_posix())
        path=ROOT/name
        if archive is not None and (archive/name).is_file():path=archive/name
        require(path.is_file() and sha(path)==row['sha256'] and path.stat().st_size==row['bytes'],
                'Changed delivery file: '+str(path))
    return len(names)


def source_hashes(mapping):
    for name,digest in mapping.items():
        require(sha(ROOT/relative(name))==digest,'Changed frozen source/input: '+name)


def remote_certificate(rows,root):
    names=set();available=Path(root).is_dir();rehashes=0
    for row in rows:
        name=relative(row['path']);digest=row['sha256']
        require(name.as_posix() not in names,'Duplicate certified checkpoint');names.add(name.as_posix())
        require(len(digest)==64 and all(c in '0123456789abcdef' for c in digest) and row['bytes']>0,
                'Invalid checkpoint hash/size certificate')
        if available:
            path=Path(root)/name
            require(sha(path)==digest and path.stat().st_size==row['bytes'],'Changed remote checkpoint: '+str(path))
            rehashes+=1
    return dict(certified=len(rows),rehashed_here=rehashes,originals_available_here=available)


def production_contract():
    def git(*args):
        return subprocess.run(['git',*args],cwd=ROOT,check=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE).stdout
    tag=git('rev-parse','v1.0.0^{}').decode().strip()
    require(tag==TAG_COMMIT,'Frozen thesis tag v1.0.0 moved')
    tracked=[p for p in git('ls-files','-z','--','pysnspd').decode().split('\0')
             if p and not p.startswith('pysnspd/experimental/')]
    unknown=[p for p in git('ls-files','--others','--exclude-standard','-z','--','pysnspd').decode().split('\0')
             if p and not p.startswith('pysnspd/experimental/')]
    require(not unknown,'Untracked production files: '+str(unknown))
    require(not git('diff','--name-only','HEAD','--',*tracked).strip(),'Production differs from HEAD')
    require(not git('diff','--name-only',tag,'HEAD','--',*tracked).strip(),'Production differs from frozen thesis tag')
    return dict(head=git('rev-parse','HEAD').decode().strip(),preserved_tag_commit=tag,
                tracked_production_files=len(tracked),worktree_equals_HEAD=True,HEAD_equals_thesis_tag=True)


def historical_chain():
    # The predecessor routine only checks hashes/receipts; it does not run tests.
    oldest=runpy.run_path(str(ROOT/'sandbox/stage4_core/verify_followup.py'))['predecessor']()
    stage=ROOT/'docs/implementation/stage4'
    pairs=[('followup_20260923','spatial_energy_20260924'),
           ('spatial_energy_20260924','self_consistent_review_20260924'),
           ('self_consistent_review_20260924','moment_review_20260924'),
           ('moment_review_20260924','time_review_20260924')]
    counts={}
    for old,new in pairs:
        counts[old]=records(read(stage/old/'delivery_manifest.json')['files'],stage/new/'previous_delivery_exact')
    require(oldest['verified'],'Earlier historical delivery verification failed')
    return dict(earlier_chain_verified=True,files_by_delivery=counts)


def verify_contents():
    raw=DATA/'raw';identity=read(raw/'identity.json');summary=read(raw/'summary.json')
    plan=read(raw/'executed_plan.json');receipt=read(raw/'extraction_receipt.json');analysis=read(DATA/'analysis.json')
    require(summary['status']=='FULL_NODE_AFFINE_THERMAL_TIME_COMPLETE','User-run time trajectory incomplete')
    require(summary['horizon_ps']==1. and summary['operator_batches']==246,'Different registered user run')
    require(identity['plan_sha256']==sha(raw/'executed_plan.json')==sha(MOMENT/'thermal_time/plan.json'),'Changed user-run plan')
    source_hashes(identity['sources']);source_hashes(receipt['sources_matching_identity'])
    require(receipt['status']=='ALL_THERMAL_TIME_CHECKPOINTS_VERIFIED_NO_SOLVES','Extraction incomplete')
    require(receipt['new_physics_solves']==0 and receipt['operator_inputs_verified']==257,'Extraction is not the recorded read-only pass')
    require(receipt['script_sha256']==sha(ROOT/'sandbox/stage4_core/extract_time_results.py'),'Extraction source changed')
    for name,digest in receipt['copied_file_sha256'].items():
        require(sha(raw/relative(name))==digest,'Changed copied user-run file: '+name)
    checkpoints=remote_certificate(receipt['checkpoint_certificate'],receipt['raw_directory'])
    require(checkpoints['certified']==37 and len(receipt['accepted_steps'])==21,'Wrong checkpoint/step count')
    require(all(row['arrays_finite'] for row in receipt['checkpoint_certificate']),'Nonfinite checkpoint certificate')
    require(receipt['maximum_accepted_relative_displacement']<=plan['maximum_relative_displacement'],'Weak-domain exit was not admitted')
    require(summary['refinement']['all_comparisons_met'] and analysis['all_registered_refinement_comparisons_met'],
            'Registered temporal refinement did not pass')
    require(analysis['checkpoints_verified']==37 and analysis['operator_inputs_verified']==257,'Analysis has different input counts')
    require(analysis['maximum_independent_summary_metric_difference']<1e-12
            and analysis['maximum_independent_refinement_metric_difference']<1e-12,'Independent compact-field analysis disagrees')
    for name,digest in analysis['provenance'].items():
        path=ROOT/'sandbox/stage4_core/analyze_time_results.py' if name=='script_sha256' else DATA/relative(name)
        require(sha(path)==digest,'Changed analysis provenance: '+name)
    require(not analysis['stage4_complete'] and not analysis['nonlinear_transient_admitted'],'Overstated affine trajectory scope')

    folder=DATA/'nonlinear_snapshots';nid=read(folder/'identity.json');ns=read(folder/'summary.json')
    nc=read(folder/'verification_receipt.json');np=read(folder/'stable_plan.json');failed_plan=read(folder/'plan.json')
    require(nid['plan_sha256']==sha(folder/'stable_plan.json'),'Changed stable nonlinear snapshot plan')
    source_hashes(nid['sources'])
    require(ns['status']=='NONLINEAR_THERMAL_SNAPSHOTS_COMPLETE' and ns['all_registered_response_margins_met'],
            'Nonlinear constitutive snapshots not admitted')
    require(len(ns['states'])==9 and len(ns['records'])==256,'Incomplete nonlinear states/frequencies')
    require(sum(row['nonlinear_roots'] for row in ns['records'])==2048,'Unexpected nonlinear root count')
    require(nc['status']=='NONLINEAR_SNAPSHOTS_CERTIFIED','Missing nonlinear integrity certificate')
    require(not nc['spectral_tolerance_changed'] and not nc['physical_equations_changed'],'Changed model or tolerance in cancellation fix')
    require(np['spectral_tolerance']==failed_plan['spectral_tolerance'],'Snapshot residual tolerance changed')
    require(nc['certificate_source_sha256']==sha(ROOT/'sandbox/stage4_core/certify_nonlinear_snapshots.py'),'Nonlinear certificate source changed')
    for name,digest in nc['compact_files'].items():
        path=folder/('stable_plan.json' if name=='executed_plan.json' else name)
        require(sha(path)==digest,'Changed compact nonlinear result: '+name)
    nonlinear=remote_certificate(nc['frequency_modes'],nc['full_root'])
    require(nonlinear['certified']==256,'Wrong certified nonlinear frequency count')
    require(sorted(row['n'] for row in nc['frequency_modes'])==list(range(256)),'Duplicate/missing nonlinear frequency')
    for row in ns['records']:
        cert=next(c for c in nc['frequency_modes'] if c['n']==row['n'])
        require(cert['path']==row['fields_path'] and cert['sha256']==row['fields_sha256'],'Nonlinear summary/certificate mismatch')
    require(nc['original_failed_root']==np['prior_failed_output'],'Original failed output identity changed')
    require(nc['original_failure']['type']=='RuntimeError' and 'energy line search failed' in nc['original_failure']['reason'],
            'Original nonlinear failure was lost')
    failed=remote_certificate(nc['original_preserved_modes'],nc['original_failed_root'])
    require(failed['certified']==45,'Original completed frequency checkpoints not preserved')
    oldroot=Path(nc['original_failed_root'])
    if oldroot.is_dir():
        require(sha(oldroot/'identity.json')==nc['original_identity_sha256'],'Changed failed-run identity')
        require(read(oldroot/'failure.json')==nc['original_failure'],'Changed failed-run cause')
    require(sha(folder/'roundoff_fixture.npz')==np['roundoff_regression_fixture_sha256'],'Roundoff regression fixture changed')
    require(not ns['production_changed'] and not ns['nonthermal_work_closed'] and not ns['detector_validated'],
            'Nonlinear snapshot scope overstated')

    coupling=DATA/'longitudinal_coupling';test=read(coupling/'unit_tests_receipt.json')
    require(test['status']=='PASSED' and test['tests']==6 and test['failures']==test['errors']==0,'Longitudinal focused tests incomplete')
    source_hashes(test['sources']);require(sha(coupling/test['log_path'])==test['log_sha256'],'Longitudinal test log changed')
    temporal=DATA/'nonlinear_time';pilot_raw=temporal/'pilot_raw'
    pilot_receipt=read(temporal/'pilot_receipt.json')
    certificate_path=ROOT/'sandbox/stage4_core/certify_nonlinear_time_pilot.py'
    require(pilot_receipt['source_sha256']==sha(certificate_path),'Pilot certificate source changed')
    # The certificate only reloads arrays and checks hashes/finiteness/contacts.
    # Calling it performs no tests, roots, trajectories or automatic retries.
    checked_pilot=runpy.run_path(str(certificate_path))['certificate']()
    require(checked_pilot==pilot_receipt,'Pilot fields, provenance or recorded scope differ from certificate')
    require(pilot_receipt['status']=='PILOT_INTEGRITY_PASSED' and pilot_receipt['horizon_ps']==.001
            and pilot_receipt['roots']==4096 and not pilot_receipt['long_trajectory_executed'],
            'Pilot was confused with a completed nonlinear trajectory')
    accepted=[]
    for phase in ('primary','refined'):
        history=read(pilot_raw/phase/'summary.json')['history']
        rows=[r for r in history if r['accepted']]
        require(len(rows)==1 and rows[0]['start_ps']==0. and rows[0]['step_ps']==.001,
                'Different short-pilot integration pattern')
        accepted.append(rows[0])
    require(accepted[0]['fields_sha256']==accepted[1]['fields_sha256'],
            'The registered pilot no longer has two identical accepted states')
    etd=read(temporal/'etd2_tests_receipt.json')
    require(etd['status']=='PASSED' and etd['tests']==9,'Prepared ETD2 tests incomplete')
    source_hashes(etd['sources'])
    require(sha(temporal/'etd2_tests.log')==etd['log_sha256'],'ETD2 test log changed')
    return dict(user_run=checkpoints,nonlinear_success=nonlinear,nonlinear_failed_preserved=failed,
                thermal_time_horizon_ps=1.,nonlinear_snapshots=9,nonlinear_roots=2048,
                longitudinal_tests=6,nonlinear_time_pilot_horizon_ps=.001,nonlinear_time_pilot_roots=4096,
                nonlinear_time_pilot_fields=len(pilot_receipt['checked_fields']),etd2_tests=9,
                pilot_proves_temporal_convergence=False,nonlinear_long_trajectory_executed=False,
                tests_executed_by_verifier=0,physical_solves_executed_by_verifier=0)


def main():
    manifest=read(DATA/'delivery_manifest.json');count=records(manifest['files'])
    require(not manifest['stage4_complete'] and not manifest['production_changed'],'Incorrect delivery scope flags')
    require(manifest['thermal_time_completed'],'Completed thermal time trajectory missing from delivery')
    require(manifest.get('thermal_nonlinear_time_prepared',False),'Prepared nonlinear trajectory driver missing from delivery')
    require(not manifest.get('thermal_nonlinear_time_completed',False),'No completed nonlinear trajectory is registered by this delivery')
    history=historical_chain();content=verify_contents();production=production_contract()
    qa=read(DATA/'report_qa.json')
    require(qa['status']=='PASS' and qa['pdf_sha256']==sha(ROOT/'output/pdf/implementation/Informe_etapa_4_evolucion_termica.pdf'),
            'Report PDF does not match visual QA')
    external=Path('/home/jdiaz/GEMINGA_COMMANDS.md')
    if sys.platform.startswith('linux') and external.exists():
        require(sha(external)==sha(ROOT/'docs/GEMINGA_COMMANDS.md'),'Geminga command notebook differs from delivery')
    print(json.dumps(dict(verified=True,manifest_sha256=sha(DATA/'delivery_manifest.json'),current_files=count,
        historical_chain=history,contents=content,production=production,stage4_complete=False,production_changed=False),indent=2))


if __name__=='__main__':main()
