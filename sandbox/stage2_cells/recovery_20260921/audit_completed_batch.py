"""Read-only numerical/provenance audit of the completed pre-recovery outputs."""
from pathlib import Path
import argparse
import hashlib
import json
import sys
import time

ROOT=Path(__file__).resolve().parents[3]
sys.path[:0]=[str(ROOT),str(ROOT/'sandbox/stage2_cells')]
import numpy as np
from assess_time_refinement import (load_run,require_same_problem,shared_indices,
                                   checkpoint_errors,invariant_checks)
from assess_grid_refinement import read as strict_read

FIELDS=('amplitudes','excitation_energy','electron_energy','phonon_energy',
        'escape','input','electron_to_phonon','condensate_heat','transport')


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--inputs',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();args.inputs=args.inputs.resolve();started=time.perf_counter()
    criteria_path=ROOT/'docs/implementation/stage2/acceptance_criteria.json'
    criteria=json.loads(criteria_path.read_text());checks={};hashes={};summaries=[]
    remote_batch='/home/jdiaz/pysnspd/tmp/stage2_validation_20260921_203419/'
    remote_root='/home/jdiaz/pysnspd/'
    def resolve(name):
        if name.startswith(remote_batch):return args.inputs/Path(name).name
        if name.startswith(remote_root):return ROOT/name.removeprefix(remote_root)
        raise ValueError('unexpected remote path '+name)
    manifest_path=args.inputs/'batch_manifest.json'
    manifest=json.loads(manifest_path.read_text());manifest_sha=sha(manifest_path)
    hashes[manifest_path.relative_to(ROOT).as_posix()]=manifest_sha
    for kind in ('one','two'):
        assessment_path=args.inputs/f'{kind}_time_assessment.json'
        assessment=json.loads(assessment_path.read_text())
        hashes[assessment_path.relative_to(ROOT).as_posix()]=sha(assessment_path)
        reference=load_run(resolve(assessment['reference_path']))
        strict_read(reference['path'])
        checks[kind+'_reference_record_hash']=sha(reference['path'])==assessment['reference_sha256']
        runs=[];rows=[]
        for case in assessment['cases']:
            run=load_run(resolve(case['path']));strict_read(run['path'])
            require_same_problem(run,reference);runs.append(run)
            checks[kind+f"_{case['steps']}_record_hash"]=sha(run['path'])==case['sha256']
        common=shared_indices(runs,reference)
        for index,run in enumerate(runs):
            record=run['record'];steps=record['parameters']['steps']
            errors=[]
            for j,matches in common:
                measured=checkpoint_errors(run,matches[index],reference,j,FIELDS)
                for name,value in measured.items():
                    if value['assessment_error'] is None:raise ValueError('undefined comparison '+name)
                    errors.append(dict(observable=name,time=float(reference['times'][j]),
                                       error=value['assessment_error']))
            worst=max(errors,key=lambda v:v['error'])
            original=next(case for case in assessment['cases'] if case['steps']==steps)
            checks[kind+f'_{steps}_recomputed_error_matches']=worst['error']==original['max_error']
            checks[kind+f'_{steps}_all_invariants_pass']=all(invariant_checks(record,criteria).values())
            receipt_path=args.inputs/f'{kind}_time_{steps}.receipt.json'
            receipt=json.loads(receipt_path.read_text())
            checks[kind+f'_{steps}_receipt_manifest']=receipt['contract']['batch_manifest_sha256']==manifest_sha
            actual_parameters={k:v for k,v in record['parameters'].items() if k!='output'}
            checks[kind+f'_{steps}_receipt_parameter_contract']=actual_parameters==receipt['contract']['task']['parameters']
            for name,expected in receipt['output_hashes'].items():
                local=resolve(name);digest=sha(local)
                checks[kind+f'_{steps}_receipt_'+Path(name).suffix.lstrip('.')+'_'+Path(name).stem]=digest==expected
                hashes[local.relative_to(ROOT).as_posix()]=digest
            rows.append(dict(steps=steps,worst=worst,ledger=record['energy_ledger_scaled_max'],
                instantaneous_balance=record['instantaneous_residual_max'],
                physical_populations=dict(minimum_electron=record['minimum_electron'],
                    maximum_electron=record['maximum_electron'],minimum_phonon=record['minimum_phonon']),
                runtime_seconds=record['runtime_seconds'],
                reused=bool(receipt['contract']['task'].get('reuse_from'))))
        ref_record=reference['record']
        checks[kind+'_reference_all_invariants_pass']=all(invariant_checks(ref_record,criteria).values())
        checks[kind+'_criteria_unchanged']=assessment['criteria_sha256']==sha(criteria_path)
        checks[kind+'_reviewer_source_unchanged']=assessment['reviewer_sha256']==sha(ROOT/'sandbox/stage2_cells/assess_time_refinement.py')
        checks[kind+'_declared_PASS']=assessment['status']=='PASS'
        rows.sort(key=lambda row:row['steps'])
        ratios=[a['worst']['error']/b['worst']['error'] for a,b in zip(rows[:-1],rows[1:])]
        checks[kind+'_independent_gate']=rows[-1]['worst']['error']<=1e-4 and all(r>=1.5 for r in ratios)
        summaries.append(dict(case=kind,source_contract_verified=True,
            reference_steps=ref_record['parameters']['steps'],
            reference_ledger=ref_record['energy_ledger_scaled_max'],
            shared_checkpoint_count=len(common),rows=rows,
            successive_error_reductions=ratios,finest_error=rows[-1]['worst']['error']))
        for path in (reference['path'],reference['path'].with_suffix('.npz'),
                     reference['path'].with_name(reference['path'].stem+'_initial.npz')):
            hashes[path.relative_to(ROOT).as_posix()]=sha(path)
    e2_path=args.inputs/'two_electron_2.json';e2=strict_read(e2_path)
    e2receipt=json.loads((args.inputs/'two_electron_2.receipt.json').read_text())
    checks['electron2_complete']=e2['record']['status']=='COMPLETED_NOT_YET_ADJUDICATED'
    checks['electron2_current_physical_sources']=all(sha(ROOT/name)==digest for name,digest in e2['record']['source_hashes'].items())
    checks['electron2_receipt_parameters']={k:v for k,v in e2['record']['parameters'].items() if k!='output'}==e2receipt['contract']['task']['parameters']
    for name,expected in e2receipt['output_hashes'].items():
        path=resolve(name);hashes[path.relative_to(ROOT).as_posix()]=sha(path)
        checks['electron2_receipt_'+Path(name).name]=sha(path)==expected
    e4_path=args.inputs/'two_electron_4.json';e4=json.loads(e4_path.read_text())
    hashes[e4_path.relative_to(ROOT).as_posix()]=sha(e4_path)
    failure=dict(status=e4['status'],reason=e4.get('reason'),runtime_seconds=e4['runtime_seconds'],
                 steps=e4['parameters']['steps'],electron_refinement=e4['parameters']['electron_refinement'],
                 trajectory_archive_present=e4_path.with_suffix('.npz').exists(),
                 admitted=False)
    result=dict(schema='pysnspd.stage2.completed-batch-audit.v1',
        status='PASS_READ_ONLY_AUDIT' if all(checks.values()) else 'FAIL_AUDIT',
        stage2_admission=False,remote_batch=remote_batch,
        criteria_sha256=sha(criteria_path),checks=checks,temporal=summaries,
        completed_electron2=dict(steps=e2['record']['parameters']['steps'],
            electron_states=e2['record']['occupation_mesh']['electron_states'],
            ledger=e2['record']['energy_ledger_scaled_max'],
            instantaneous_balance=e2['record']['instantaneous_residual_max'],
            reuse_at_identical_steps_only=True),
        failed_finest_mesh=failure,artifact_sha256=hashes,
        recovery_rules=[
            'Preserve the failed e4 trajectory registration and initial archive; no population clipping.',
            'A smaller timestep on e4 alone cannot be mixed into the existing same-time grid assessor. All compared electronic meshes must use the same method, steps and tolerance.',
            'Candidate e1 at320 steps already exists and may be reused if320 is selected; e2 at160 cannot stand in for e2 at320.',
            'Positivity of one smaller-step trial proves local stage admissibility only, not full time convergence or a continuous-grid PASS.',
            'If finer grids reveal additional fast dynamics, independently bound their temporal error with a tighter complete trajectory before attributing differences to mesh resolution.',
            'The candidate630/1025 must pass its own population/observable errors and convergence; a successful finer trajectory or conserved energy alone cannot admit it.'
        ],source_sha256=sha(__file__),runtime_seconds=time.perf_counter()-started)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(dict(status=result['status'],checks_passed=sum(checks.values()),checks_total=len(checks),
                         temporal=[dict(case=s['case'],finest_error=s['finest_error'],ratios=s['successive_error_reductions']) for s in summaries]),indent=2))
    if not all(checks.values()):raise SystemExit(1)


if __name__=='__main__':main()
