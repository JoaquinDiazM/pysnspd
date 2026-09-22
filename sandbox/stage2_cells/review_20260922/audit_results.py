"""Read-only audit of collected stage-2 batches; never evaluates a physical RHS."""
from pathlib import Path
import argparse
import hashlib
import json
import sys
import time

ROOT=Path(__file__).resolve().parents[3]
sys.path[:0]=[str(ROOT),str(ROOT/'sandbox/stage2_cells')]
import numpy as np
from assess_grid_refinement import read as strict_read
from assess_time_refinement import load_run,require_same_problem,shared_indices,checkpoint_errors,invariant_checks

CRITERIA=ROOT/'docs/implementation/stage2/acceptance_criteria.json'
FROZEN_CRITERIA='48a56a76b64b58b81b176325a6a34cc2e6c690bfde3b4eb739a8cd2723c702cf'
FIELDS=('amplitudes','excitation_energy','electron_energy','phonon_energy',
        'escape','input','electron_to_phonon','condensate_heat','transport')


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def read(path):return json.loads(Path(path).read_text(encoding='utf-8-sig'))
def relative(path):return Path(path).resolve().relative_to(ROOT).as_posix()


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--raw',type=Path,default=ROOT/'docs/implementation/stage2/review_20260922/raw')
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();args.raw=args.raw.resolve();started=time.perf_counter()
    criteria=read(CRITERIA);checks={'criteria_exactly_preserved':sha(CRITERIA)==FROZEN_CRITERIA}
    hashes={relative(path):sha(path) for path in sorted(args.raw.rglob('*')) if path.is_file()}
    candidates=[*args.raw.rglob('*'),*(ROOT/'docs/implementation/stage2/recovery_20260921/audit_inputs').glob('*'),
                *(ROOT/'docs/implementation/stage2/resume_20260921').glob('*')]
    candidates=[p for p in candidates if p.is_file()]
    def resolve(remote,digest=None):
        possibilities=[p for p in candidates if p.name==Path(remote).name]
        if remote.startswith('/home/jdiaz/pysnspd/'):
            local=ROOT/remote.removeprefix('/home/jdiaz/pysnspd/')
            if local.is_file():possibilities.insert(0,local)
        for path in possibilities:
            if digest is None or sha(path)==digest:return path
        # The same registered initial state can be archived under a different
        # trajectory name. Only an exact byte hash permits that alias.
        if digest is not None:
            for path in candidates:
                if path.suffix==Path(remote).suffix and sha(path)==digest:return path
        raise ValueError('Cannot resolve registered artifact '+remote)
    checkpoint_path=ROOT/'docs/implementation/stage2/stage2_admission.json'
    checkpoint=read(checkpoint_path)
    inherited_static=[]
    dynamic_ids={'simultaneous_trajectories_and_time_convergence','dynamic_electron_phonon_resolution',
                 'selected_boundary_equilibrium_and_activity','actual_trajectory_fields_and_support'}
    for gate in checkpoint['gates']:
        if gate['id'] not in dynamic_ids:
            inherited_static.append(dict(id=gate['id'],recorded_status=gate['status'],detail=gate['detail'],
                provenance='Retained historical gate from '+relative(checkpoint_path)))
    for name,digest in checkpoint.get('audited_source_sha256',{}).items():
        checks['historical_source:'+name]=sha(ROOT/name)==digest
    for key in ('R2_catalogue','R2_query_source','stage1_closure_certificate','stage1_criteria'):
        checks['immutable:'+key]=sha(ROOT/criteria['immutable_inputs'][key+'_path'])==criteria['immutable_inputs'][key+'_sha256']
    manifests=[];receipts=[];batches=[];trajectories=[]
    for family in ('rk4','ssp'):
        directory=args.raw/family;manifest_path=directory/'batch_manifest.json';manifest=read(manifest_path)
        plan_path=ROOT/('docs/implementation/stage2/resume_20260921/manual_validation_plan.json'
            if family=='rk4' else 'docs/implementation/stage2/recovery_20260921/limited_acceptance_plan.json')
        plan=read(plan_path);registered_tasks={task['id']:task for task in plan['tasks']}
        checks[family+':frozen_plan_hash']=sha(plan_path)==manifest['plan_sha256']
        for name,digest in manifest['source_hashes'].items():checks[family+':manifest:'+name]=sha(ROOT/name)==digest
        manifests.append(dict(path=relative(manifest_path),sha256=sha(manifest_path),
                              plan_sha256=manifest['plan_sha256']))
        events=[json.loads(line) for line in (directory/'batch_events.jsonl').read_text().splitlines() if line.strip()]
        completed=[event['task'] for event in events if event['event']=='TASK_COMPLETE']
        started_tasks=[event['task'] for event in events if event['event']=='TASK_START']
        batches.append(dict(family=family,complete_tasks=completed,completed_task_count=len(completed),
            registered_task_count=len(registered_tasks),
            started_tasks=started_tasks,terminal_event=events[-1],
            batch_complete=any(event['event']=='BATCH_COMPLETE' for event in events)))
        for receipt_path in sorted(directory.glob('*.receipt.json')):
            receipt=read(receipt_path);contract=receipt['contract'];verified={}
            verified['manifest']=contract['batch_manifest_sha256']==sha(manifest_path)
            task=contract['task']
            verified['registered_task']=registered_tasks.get(task['id'])==task
            for remote,digest in receipt['output_hashes'].items():
                local=resolve(remote,digest);verified[remote]=sha(local)==digest
                hashes.setdefault(relative(local),sha(local))
                if local.suffix=='.json' and task['kind'] in ('trajectory','existing_reference'):
                    record=read(local)
                    verified['parameters_except_output']={k:v for k,v in record['parameters'].items() if k!='output'}==task['parameters']
            receipts.append(dict(path=relative(receipt_path),verified=all(verified.values()),checks=verified))
            checks[family+':receipt:'+receipt_path.name]=all(verified.values())
        for path in sorted(directory.glob('*.json')):
            record=read(path)
            if record.get('schema')!='pysnspd.stage2.coupled_run.v1':continue
            row=dict(path=relative(path),sha256=sha(path),status=record['status'],parameters=record['parameters'],
                runtime_seconds=record.get('runtime_seconds'),completed=record['status']=='COMPLETED_NOT_YET_ADJUDICATED')
            checks[family+':sources:'+path.name]=all(sha(ROOT/name)==digest for name,digest in record['source_hashes'].items())
            if row['completed']:
                loaded=strict_read(path,ROOT);states=loaded['arrays']['states'];snapshots=record['snapshots']
                scale=max(1.,abs(snapshots[0]['total']),float(np.sum(snapshots[-1]['input'])))
                measured_ledger=max(abs(s['conserved']-snapshots[0]['conserved']) for s in snapshots)/scale
                checks[family+':recomputed_ledger:'+path.name]=measured_ledger==record['energy_ledger_scaled_max']
                row.update(invariant_checks=invariant_checks(record,criteria),
                    ledger=measured_ledger,ledger_limit=criteria['coupled_trajectories']['energy_ledger_scaled_max'],
                    ledger_times_limit=measured_ledger/criteria['coupled_trajectories']['energy_ledger_scaled_max'],
                    instantaneous_balance=record['instantaneous_residual_max'],
                    initial_amplitudes=record['initial']['amplitudes'],final_amplitudes=record['final']['amplitudes'],
                    minimum_electron=record['minimum_electron'],maximum_electron=record['maximum_electron'],
                    minimum_phonon=record['minimum_phonon'],time_integration=record.get('time_integration'),
                    physical_archive_verified=True)
            else:
                initial=path.with_name(path.stem+'_initial.npz')
                checks[family+':failed_initial:'+path.name]=sha(initial)==record['initial_conditions_sha256']
                row.update(reason=record.get('reason'),exception=record.get('exception'),
                    completed_trajectory_present=path.with_suffix('.npz').exists(),numerical_admission=False)
            trajectories.append(row)
    temporal=[]
    for case in ('one','two'):
        assessment_path=args.raw/'rk4'/(case+'_time_assessment.json');assessment=read(assessment_path)
        reference=load_run(resolve(assessment['reference_path'],assessment['reference_sha256']))
        runs=[load_run(resolve(row['path'],row['sha256'])) for row in assessment['cases']]
        for run in [*runs,reference]:strict_read(run['path'],ROOT)
        for run in runs:require_same_problem(run,reference)
        common=shared_indices(runs,reference);rows=[]
        for index,run in enumerate(runs):
            errors=[]
            for j,indices in common:
                measured=checkpoint_errors(run,indices[index],reference,j,FIELDS)
                for name,value in measured.items():
                    if value['assessment_error'] is None:raise ValueError('Undefined time comparison: '+name)
                    errors.append(dict(observable=name,error=value['assessment_error'],time=float(reference['times'][j])))
            worst=max(errors,key=lambda row:row['error']);steps=run['record']['parameters']['steps']
            original=next(row for row in assessment['cases'] if row['steps']==steps)
            checks[case+':independent_max_error:'+str(steps)]=worst['error']==original['max_error']
            row=dict(steps=steps,worst=worst,invariants=invariant_checks(run['record'],criteria))
            rows.append(row)
        rows.sort(key=lambda row:row['steps'])
        reductions=[a['worst']['error']/b['worst']['error'] for a,b in zip(rows[:-1],rows[1:])]
        limits=criteria['coupled_trajectories']
        passed=(len(rows)>=limits['minimum_time_resolutions'] and
            rows[-1]['worst']['error']<=limits['finest_time_observable_relative_error_max'] and
            all(value>=limits['time_refinement_minimum_error_reduction'] for value in reductions) and
            all(all(row['invariants'].values()) for row in rows) and
            all(invariant_checks(reference['record'],criteria).values()))
        checks[case+':original_status_matches_recomputed']=passed==(assessment['status']=='PASS')
        checks[case+':time_assessor_source']=assessment['reviewer_sha256']==sha(ROOT/'sandbox/stage2_cells/assess_time_refinement.py')
        temporal.append(dict(case=case,status='PASS_CANDIDATE_RK4_TIME_ONLY' if passed else 'FAIL',
            steps=[row['steps'] for row in rows],reference_steps=reference['record']['parameters']['steps'],
            common_checkpoint_count=len(common),rows=rows,finest_error=rows[-1]['worst']['error'],
            tolerance=limits['finest_time_observable_relative_error_max'],error_reductions=reductions,
            scope='Fixed 630-electron/1025-phonon mesh, registered duration2 scenario only; no finer-grid or SSP inference.'))
    diagnostic=args.raw/'diagnostics/selected_continuum.json';selected=read(diagnostic)
    old_selected_path=ROOT/'docs/implementation/stage2/resume_20260921/selected_continuum.json';old_selected=read(old_selected_path)
    differences=[key for key in sorted(set(selected)|set(old_selected)) if selected.get(key)!=old_selected.get(key)]
    checks['selected_continuum_only_runtime_changed']=differences==['runtime_seconds']
    mapping={'criteria':CRITERIA,'measured':ROOT/'docs/implementation/stage2/selected_event_measurements.json',
             'catalog':ROOT/'docs/implementation/stage1_r2/catalogs/occupation_catalog.npz',
             'script':ROOT/'sandbox/stage2_cells/resume_20260921/close_selected_continuum.py',
             'kernel':ROOT/'pysnspd/experimental/kinetic_events.py','factory':ROOT/'pysnspd/experimental/refined_cells.py'}
    for name,digest in selected['hashes'].items():
        path=mapping.get(name,ROOT/'docs/implementation/stage2/review'/name)
        checks['selected_continuum_source:'+name]=sha(path)==digest
    weak=selected['weak_comparisons'];moments=selected['moments_and_hats'];hats=[h for row in moments for h in row['hats'].values()]
    max_error=max([row['relative_L1'] for row in weak]+[e for row in moments for e in row['relative_errors']]+[h['relative_L1'] for h in hats])
    max_budget=max([row['reference_refinement_relative_L1'] for row in weak]+[e for row in moments for e in row['reference_refinement_relative_errors']]+[h['reference_refinement_relative_L1'] for h in hats])
    selected_pass=max_error<=criteria['continuous_consistency']['relative_error_max'] and max_budget<=criteria['continuous_consistency']['reference_refinement_relative_budget_max']
    checks['selected_continuum_thresholds_pass']=selected_pass
    smooth_path=args.raw/'diagnostics/rhs_smoothness_recheck.json';smooth=read(smooth_path)
    checks['smoothness_sources_current']=all(sha(ROOT/name)==digest for name,digest in smooth['source_hashes'].items())
    smooth_initial=resolve(smooth['initial_path'],smooth['initial_sha256'])
    checks['smoothness_initial_hash']=sha(smooth_initial)==smooth['initial_sha256']
    activity_record=read(args.raw/'rk4/two_time_160.json')
    initial=activity_record['initial'];final=activity_record['final']
    activity_scale=activity_record['initial_excitation_energy']+sum(initial['phonon_energy'])
    activity=dict(preregistered_initial_excitation_energy=activity_record['initial_excitation_energy'],
        initial_phonon_energy=sum(initial['phonon_energy']),scale=activity_scale,
        maximum_amplitude_change=max(abs(a-b) for a,b in zip(initial['amplitudes'],final['amplitudes'])),
        eph_fraction=sum(abs(value) for value in final['electron_to_phonon'])/activity_scale,
        face_fraction=abs(final['transport'])/activity_scale,
        scope='Measured from the time-converged candidate RK4 trajectory; dynamic mesh admission remains open.')
    activity['passes']=activity['maximum_amplitude_change']>=.02 and activity['eph_fraction']>=.01 and activity['face_fraction']>=.01
    missing=[
        dict(id='electronic_dynamic_mesh',status='FAILED_REFINEMENT_INCOMPLETE_GATE',reason='2520-state RK4 trajectory failed; no complete three-mesh comparison.'),
        dict(id='phonon_dynamic_mesh',status='NOT_REACHED',reason='2049/4097-node trajectories and three-grid assessment were not reached.'),
        dict(id='SSP_temporal_admission',status='FAILED_COARSE_LEDGER_GATE_NOT_CONVERGED',reason='The only complete SSP run fails the ledger gate; 80/160 runs and time comparison were not reached.'),
        dict(id='SSP_fine_mesh_paired_diagnostic',status='NOT_REACHED',reason='No fine SSP trajectory or paired temporal diagnostic exists in this batch.'),
        dict(id='equilibrium_and_population_boundary_trajectories',status='NOT_REACHED',reason='Registered final equilibrium/vacuum/sparse trajectories were not reached; static inward-RHS tests are distinct.'),
        dict(id='actual_final_trajectory_fields_and_support',status='PENDING',reason='Historical two-cell pilot checks exist; final accepted mesh/trajectory checks were not reached.'),
        dict(id='final_report_and_spatial_transition',status='NOT_ADMITTED',reason='The required numerical gates remain failed or missing.')]
    result=dict(schema='pysnspd.stage2.independent-collected-results-audit.v1',date='2026-09-22',
        status='AUDITED_STAGE2_NOT_CLOSED' if all(checks.values()) else 'PROVENANCE_OR_AUDIT_CHECK_FAILED',
        closure='NO',numerical_admission=False,production_promotion=False,final_report_eligible=False,
        physical_RHS_evaluations=0,new_trajectories=0,thresholds_changed=False,
        criteria_sha256=sha(CRITERIA),checks=checks,checked_items=len(checks),passed_checks=sum(checks.values()),
        batch_summaries=batches,manifests=manifests,receipts=receipts,trajectories=trajectories,
        temporal_candidate_RK4=temporal,activity_candidate_RK4=activity,
        selected_static_continuum=dict(status='PASS_STATIC_ONLY' if selected_pass else 'FAIL',
            changed_top_level_keys=differences,maximum_error=max_error,maximum_reference_budget=max_budget,
            source_path=relative(diagnostic),prior_source_path=relative(old_selected_path)),
        smoothness_recheck=dict(path=relative(smooth_path),scope=smooth['scope'],rhs_calls=smooth['rhs_calls'],
            verified_initial_archive=relative(smooth_initial),
            interpretation='A fixed-state finite-difference diagnostic does not prove global smoothness, exclude all stiffness, or close temporal/mesh gates.'),
        historical_static_gates=inherited_static,failed_or_missing_gates=missing,
        conclusions=[
            'The repeat RK4 batch confirms candidate-mesh temporal convergence but reproduces the incomplete finest electronic mesh.',
            'The SSP40 trajectory is complete and positive, but its integrated ledger exceeds the frozen tolerance; a completed file is not acceptance.',
            'No event limiter activated in SSP40. Its failure therefore does not establish a limiter-induced physical bias.',
            'Energy conservation and source consistency alone cannot replace the missing distribution, mesh, boundary and actual-trajectory validations.',
            'A final closure report is not justified by these outputs. Preserve this audit as a partial-results review; do not launch further long calculations under this read-only request.'
        ],artifact_sha256=hashes,auditor_sha256=sha(__file__),runtime_seconds=time.perf_counter()-started)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    print(json.dumps(dict(status=result['status'],closure=result['closure'],checks_passed=result['passed_checks'],
                         checks_total=result['checked_items'],runtime_seconds=result['runtime_seconds']),indent=2))
    if not all(checks.values()):raise SystemExit(1)


if __name__=='__main__':main()
