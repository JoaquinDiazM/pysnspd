"""Verify completed data and delivery; run no tests, spectra or trajectories."""
from pathlib import Path
import hashlib,json,runpy,subprocess,sys
ROOT=Path(__file__).resolve().parents[2]
HERE=ROOT/'docs/implementation/stage4/practical_time_review_20260924'
OLD=ROOT/'docs/implementation/stage4/time_review_20260924'
helpers=runpy.run_path(str(ROOT/'sandbox/stage4_core/verify_time_review.py'))
read,sha,require=helpers['read'],helpers['sha'],helpers['require']


def inherited_sources(mapping):
    """Keep executed hashes; tolerate checkout EOLs only for unchanged production."""
    for name,digest in mapping.items():
        if sha(ROOT/name)==digest:continue
        require(name.startswith('pysnspd/') and not name.startswith('pysnspd/experimental/'),
                'Changed exact adapter source: '+name)
        blob=subprocess.check_output(['git','show','HEAD:'+name],cwd=ROOT)
        lf=blob.replace(b'\r\n',b'\n')
        checkout_hashes={hashlib.sha256(variant).hexdigest() for variant in (blob,lf,lf.replace(b'\n',b'\r\n'))}
        require(digest in checkout_hashes,'Executed production source differs from Git, including checkout EOLs: '+name)
        require((ROOT/name).read_bytes().replace(b'\r\n',b'\n')==blob.replace(b'\r\n',b'\n'),
                'Local production source content differs: '+name)


def contents():
    raw=HERE/'raw';summary=read(raw/'summary.json');identity=read(raw/'identity.json')
    cert=read(raw/'extraction_receipt.json');analysis=read(HERE/'analysis.json')
    require(summary['status']=='TEMPORAL_REFINEMENT_NOT_MET','Original unsuccessful comparison changed')
    require(summary['horizon_ps']==1. and summary['spectral_roots']==85504,'Different trajectory')
    require(cert['status']=='CHECKPOINTS_VERIFIED_ORIGINAL_GATE_RETAINED' and cert['new_physics_solves']==0,'Invalid extraction scope')
    require(cert['script_sha256']==sha(ROOT/'sandbox/stage4_core/extract_nonlinear_time_results.py'),'Extraction script changed')
    helpers['source_hashes'](identity['sources']);helpers['source_hashes'](cert['sources_matching_identity'])
    require(identity['plan_sha256']==sha(raw/'executed_plan.json')==sha(OLD/'nonlinear_time/plan.json'),'Original plan changed')
    for name,digest in cert['copied_file_sha256'].items():require(sha(raw/name)==digest,'Changed copied result: '+name)
    certificate=helpers['remote_certificate'](cert['checkpoint_certificate'],cert['raw_directory'])
    require(certificate['certified']==71 and len(cert['accepted_steps'])==55,'Wrong checkpoint or step count')
    require(summary['refinement']==read(raw/'refinement.json'),'Changed original comparison')
    require(sum(row['admitted'] for row in summary['refinement']['records'])==95,'Changed comparison count')
    require(len(summary['refinement']['records'])==96,'Incomplete original comparison')
    require('Planned temporal refinement' in read(raw/'failure.json')['reason'],'Original failure cause not retained')
    require(analysis['original_failure_preserved'] and analysis['new_physics_solves']==0,'Original failure/scope changed')
    require(analysis['analysis_source_sha256']==sha(ROOT/'sandbox/stage4_core/analyze_nonlinear_time_results.py'),'Analysis source changed')
    require(analysis['raw_extraction_receipt_sha256']==sha(raw/'extraction_receipt.json'),'Extraction certificate changed')
    require(analysis['maximum_registered_metric_reproduction_difference']<1e-12,'Reproduced fields differ from comparison')
    require(not analysis['original_gate']['changed'] and analysis['original_gate']['failed_comparisons']==1,'Original gate was rewritten')
    decision=read(HERE/'decision.json')
    require(decision['execution_complete'] and decision['thermal_benchmark_development_closed'],'Missing explicit limited development decision')
    require(not decision['original_temporal_certificate_passed'] and not decision['stage4_complete'],'Overstated admission')
    require(not decision['long_run_requested'],'Unexpected repeated long campaign')
    mesh=read(HERE/'dual_mesh/receipt.json');tests=read(HERE/'dual_mesh/unit_tests_receipt.json')
    require(mesh['status']=='MESH_PREPARED_NO_SPECTRAL_OR_TIME_SOLVE','Mesh preparation scope changed')
    require(sha(HERE/'dual_mesh/mesh.npz')==mesh['mesh_sha256'],'Changed mesh fixture')
    inherited_sources(mesh['source_sha256']);inherited_sources(tests['source_sha256'])
    require(tests['status']=='PASS' and tests['tests']==3,'Mesh focused tests missing')
    require(sha(HERE/'dual_mesh/unit_tests.log')==tests['log_sha256'],'Mesh test log changed')
    remesh=read(HERE/'dual_mesh/resampled/receipt.json')
    require(remesh['source_mesh_sha256']==mesh['mesh_sha256'],'Original mesh lost from remeshing record')
    require(remesh['mesh_sha256']==sha(HERE/'dual_mesh/resampled/mesh.npz'),'Resampled mesh changed')
    require(remesh['runner_sha256']==sha(ROOT/'sandbox/stage4_core/prepare_dual_mesh_boundary_sampling.py'),'Boundary preparation source changed')
    require(remesh['area_relative_error']<1e-12 and not remesh['mesh_warnings'],'Boundary remeshing has unresolved geometry diagnostics')
    bridge=read(HERE/'kwt_bridge/receipt.json');inherited_sources(bridge['source_sha256'])
    require(bridge['status']=='PASSED' and bridge['direct_inherited_solver_calls']==3,'Inherited KWT update not exercised')
    require(all(row['fixed_contact_change']==0 for row in bridge['refinements']),'KWT contacts moved')
    require(all(1.8<ratio<2.2 for ratio in bridge['first_order_ratios']),'Inherited update conversion inconsistent')
    return dict(checkpoints=certificate,steps=55,completed_physical_horizon_ps=1.,original_gate_passed=False,
        comparisons_passed=95,comparisons_total=96,development_decision='Accept thermal benchmark with explicit late cross-torque limit',
        first_mesh_nodes=mesh['n_nodes'],resampled_dual_mesh_nodes=remesh['n_nodes'],
        inherited_kwt_calls=3,tests_executed_now=0,physics_solves_executed_now=0)


def main():
    manifest=read(HERE/'delivery_manifest.json');count=helpers['records'](manifest['files'])
    require(not manifest['stage4_complete'] and not manifest['production_changed'],'Delivery overstates scope')
    older=helpers['historical_chain']()
    previous=helpers['records'](read(OLD/'delivery_manifest.json')['files'],HERE/'previous_delivery_exact')
    helpers['verify_contents']()
    current=contents();production=helpers['production_contract']()
    qa=read(HERE/'report_qa.json')
    require(qa['status']=='PASS' and qa['pdf_sha256']==sha(ROOT/'output/pdf/implementation/Informe_etapa_4_revision_practica.pdf'),'PDF QA does not match')
    external=Path('/home/jdiaz/GEMINGA_COMMANDS.md')
    if sys.platform.startswith('linux') and external.exists():require(sha(external)==sha(ROOT/'docs/GEMINGA_COMMANDS.md'),'External notebook differs')
    print(json.dumps(dict(verified=True,current_files=count,previous_files=previous,older_chain=older,
        contents=current,production=production,stage4_complete=False),indent=2))


if __name__=='__main__':main()
