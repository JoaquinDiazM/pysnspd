"""Offline audit of the completed one-cell SSP temporal series; no RHS calls."""
from pathlib import Path
import argparse
import hashlib
import json
import sys
import time

ROOT=Path(__file__).resolve().parents[3]
sys.path[:0]=[str(ROOT),str(ROOT/'sandbox/stage2_cells'),
              str(ROOT/'sandbox/stage2_cells/recovery_20260921')]
import numpy as np
from assess_grid_refinement import read as strict_read
from assess_time_refinement import load_run,shared_indices,checkpoint_errors,invariant_checks
from assess_limited_time import compare as verify_same_problem,verify_contract

CRITERIA=ROOT/'docs/implementation/stage2/acceptance_criteria.json'
PLAN=ROOT/'docs/implementation/stage2/review_20260922/manual_time_plan.json'
FROZEN_CRITERIA='48a56a76b64b58b81b176325a6a34cc2e6c690bfde3b4eb739a8cd2723c702cf'
FIELDS=('amplitudes','excitation_energy','electron_energy','phonon_energy',
        'escape','input','electron_to_phonon','condensate_heat','transport')


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def read(path):return json.loads(Path(path).read_text(encoding='utf-8-sig'))
def relative(path):return Path(path).resolve().relative_to(ROOT).as_posix()


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--raw',type=Path,default=ROOT/'docs/implementation/stage2/time_pass_20260922/raw')
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();args.raw=args.raw.resolve();started=time.perf_counter()
    criteria=read(CRITERIA);plan=read(PLAN);checks={};artifacts={}
    limits=criteria['coupled_trajectories'];threshold=limits['finest_time_observable_relative_error_max']
    checks['acceptance_criteria_unchanged']=sha(CRITERIA)==FROZEN_CRITERIA
    manifest_path=args.raw/'batch_manifest.json';manifest=read(manifest_path)
    checks['batch_uses_registered_plan']=manifest['plan_sha256']==sha(PLAN)
    checks['batch_source_contract_equals_plan']=manifest['source_hashes']==plan['source_hashes']
    for name,digest in plan['source_hashes'].items():checks['source:'+name]=sha(ROOT/name)==digest
    events=[json.loads(line) for line in (args.raw/'batch_events.jsonl').read_text().splitlines() if line.strip()]
    completed=[row['task'] for row in events if row['event']=='TASK_COMPLETE']
    checks['all_five_tasks_completed_in_order']=completed==[row['id'] for row in plan['tasks']]
    checks['batch_complete_does_not_claim_stage2_admission']=(events[-1]['event']=='BATCH_COMPLETE'
        and events[-1].get('stage2_admission') is False)
    tasks={row['id']:row for row in plan['tasks']};runs={};summaries=[]
    for task in plan['tasks']:
        identifier=task['id'];path=args.raw/(identifier+'.json');record=read(path)
        receipt_path=args.raw/(identifier+'.receipt.json');receipt=read(receipt_path)
        contract=receipt['contract']
        checks[identifier+':receipt_registered_task']=contract['task']==task
        checks[identifier+':receipt_manifest']=contract['batch_manifest_sha256']==sha(manifest_path)
        for remote,digest in receipt['output_hashes'].items():
            archive=args.raw/Path(remote).name
            checks[identifier+':output:'+archive.name]=sha(archive)==digest
        for name,digest in contract['dependency_hashes'].items():
            checks[identifier+':dependency:'+name]=sha(args.raw/(name+'.json'))==digest
        if task['kind']!='trajectory':continue
        loaded=strict_read(path,ROOT);run=load_run(path);verify_contract(run)
        runs[identifier]=run
        checks[identifier+':parameter_contract']={k:v for k,v in record['parameters'].items() if k!='output'}==task['parameters']
        checks[identifier+':versions_match_manifest']=all(record[key]==manifest[key] for key in ('python','numpy','scipy'))
        snapshots=record['snapshots'];initial=record['initial'];final=record['final']
        scale=max(1.,abs(initial['total']),float(np.sum(final['input'])))
        ledger=max(abs(s['conserved']-initial['conserved']) for s in snapshots)/scale
        checks[identifier+':ledger_recomputed']=ledger==record['energy_ledger_scaled_max']
        states=loaded['arrays']['states'];ne=len(loaded['arrays']['electron_count']);nph=len(loaded['arrays']['phonon_energies'])
        populations=states[:,1:1+ne];phonons=states[:,1+ne:1+ne+nph]
        checks[identifier+':population_extrema_recomputed']=(float(populations.min())==record['minimum_electron']
            and float(populations.max())==record['maximum_electron'] and float(phonons.min())==record['minimum_phonon'])
        stats=record['time_integration']['stats']
        summaries.append(dict(id=identifier,steps=record['parameters']['steps'],
            runtime_seconds=record['runtime_seconds'],ledger=ledger,ledger_limit=limits['energy_ledger_scaled_max'],
            instantaneous_balance=record['instantaneous_residual_max'],
            invariant_checks=invariant_checks(record,criteria),
            minimum_electron=record['minimum_electron'],maximum_electron=record['maximum_electron'],
            minimum_phonon=record['minimum_phonon'],
            limiter=dict(limited_event_evaluations=stats['limited_event_evaluations'],
                event_evaluations=stats['event_evaluations'],minimum_event_factor=stats['minimum_event_factor'],
                integrated_energy_weighted_flux_defect=stats['integrated_energy_weighted_flux_defect']),
            initial_amplitude=initial['amplitudes'][0],final_amplitude=final['amplitudes'][0],
            physical_archive_verified=True))
    assessment_path=args.raw/'one_ssp_time_assessment.json';assessment=read(assessment_path)
    reference=runs['one_ssp_1280'];compared=[runs['one_ssp_'+str(n)] for n in (160,320,640)]
    checks['assessment_reference_hash']=assessment['reference_sha256']==sha(reference['path'])
    checks['assessment_criteria']=assessment['criteria_sha256']==sha(CRITERIA)
    checks['base_assessor_source']=assessment['reviewer_sha256']==sha(ROOT/'sandbox/stage2_cells/assess_time_refinement.py')
    checks['limited_assessor_source']=assessment['integration_contract_comparison']['assessor_sha256']==sha(ROOT/'sandbox/stage2_cells/recovery_20260921/assess_limited_time.py')
    for run in compared:verify_same_problem(run,reference)
    common=shared_indices(compared,reference);rows=[]
    for index,run in enumerate(compared):
        errors=[]
        for j,matches in common:
            values=checkpoint_errors(run,matches[index],reference,j,FIELDS)
            for name,value in values.items():
                if value['assessment_error'] is None:raise ValueError('Undefined observable comparison: '+name)
                errors.append(dict(observable=name,time=float(reference['times'][j]),error=value['assessment_error']))
        worst=max(errors,key=lambda row:row['error']);steps=run['record']['parameters']['steps']
        stored=next(case for case in assessment['cases'] if case['steps']==steps)
        checks['assessment_case_hash:'+str(steps)]=stored['sha256']==sha(run['path'])
        checks['error_recomputed:'+str(steps)]=stored['max_error']==worst['error']
        per_observable={name:max(row['error'] for row in errors if row['observable']==name)
                        for name in [*FIELDS,'electronic_population','phonon_population']}
        invariants=invariant_checks(run['record'],criteria)
        own_precision=worst['error']<=threshold
        rows.append(dict(steps=steps,dt=run['record']['parameters']['duration']/steps,
            worst=worst,maximum_errors_by_observable=per_observable,relative_error_limit=threshold,
            own_precision_pass=own_precision,invariant_checks=invariants,
            own_time_resolution_admissible=own_precision and all(invariants.values()),
            comparison_scope='Registered one-cell duration2, 630 electronic/1025 phonon states, reference SSP1280.'))
    reductions=[a['worst']['error']/b['worst']['error'] for a,b in zip(rows[:-1],rows[1:])]
    convergence=all(value>=limits['time_refinement_minimum_error_reduction'] for value in reductions)
    series_pass=(len(rows)>=limits['minimum_time_resolutions'] and rows[-1]['own_precision_pass'] and convergence
        and all(all(row['invariant_checks'].values()) for row in rows)
        and all(invariant_checks(reference['record'],criteria).values()))
    checks['reported_PASS_equals_recomputed_series']=series_pass==(assessment['status']=='PASS')
    checks['reported_reductions_equal_recomputed']=assessment['successive_max_error_reduction']==reductions
    checks['reported_shared_times_equal_recomputed']=assessment['shared_checkpoint_times']==[float(reference['times'][j]) for j,_ in common]
    checks['reported_finest_error_equal_recomputed']=assessment['finest_max_error']==rows[-1]['worst']['error']
    for path in args.raw.rglob('*'):
        if path.is_file():artifacts[relative(path)]=sha(path)
    eligible=[row['steps'] for row in rows if row['own_time_resolution_admissible']]
    all_inactive=all(row['limiter']['limited_event_evaluations']==0 for row in summaries)
    result=dict(schema='pysnspd.stage2.one-cell-time-offline-audit.v1',date='2026-09-22',
        status='AUDITED_ONE_CELL_TIME_PASS' if all(checks.values()) and series_pass else 'AUDIT_OR_TIME_GATE_FAILED',
        stage2_closure=False,production_admission=False,new_RHS_evaluations=0,new_trajectories=0,
        criteria_sha256=sha(CRITERIA),plan_sha256=sha(PLAN),checks=checks,
        checks_passed=sum(checks.values()),checks_total=len(checks),
        batch_elapsed_seconds=events[-1]['elapsed_seconds'],series_time_gate_pass=series_pass,
        registered_scope=dict(case='one',scenario='driven',duration=2.,electron_states=630,phonon_states=1025,
            gamma=0.,escape='inf',heating=0.,method='ssprk3_common_flux_limited'),
        reference=dict(steps=1280,rule='Separately computed half-step reference to640; allowed by frozen criteria, not an independent integration method.'),
        shared_checkpoint_count=len(common),rows=rows,error_reductions=reductions,trajectories=summaries,
        admissible_measured_time_resolutions=eligible,
        minimum_admissible_measured_steps=min(eligible) if eligible else None,
        all_measured_limiters_inactive=all_inactive,
        conclusions=[
            'A PASS of the temporal sequence applies to its finest measured comparison and does not automatically admit every coarser timestep.',
            'The explicit per-resolution verdicts above apply only to this one-cell physical initial-value problem and candidate population grid.',
            'If limiter activity is zero throughout, this series does not validate the active-limiter fine-grid case that previously motivated the recovery.',
            'Two-cell SSP time accuracy, active-limiter finest-mesh checks, full electronic/phonon mesh convergence, final boundaries/equilibrium and actual-trajectory field/support checks remain separate gates.',
            'Do not issue a global stage2 closure or proceed to stage3 implementation from this one-cell result alone.'
        ],artifact_sha256=artifacts,auditor_sha256=sha(__file__),runtime_seconds=time.perf_counter()-started)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    print(json.dumps({key:result[key] for key in ('status','checks_passed','checks_total','series_time_gate_pass',
        'admissible_measured_time_resolutions','minimum_admissible_measured_steps','all_measured_limiters_inactive',
        'runtime_seconds')},indent=2))
    if not all(checks.values()):raise SystemExit(1)


if __name__=='__main__':main()
