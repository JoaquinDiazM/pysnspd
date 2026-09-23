"""Supplemental two-resolution check, never a three-level temporal admission."""
from pathlib import Path
import argparse,json
import assess_guarded_time as contract

original=contract.original
FIELDS=('amplitudes','excitation_energy','electron_energy','phonon_energy',
        'escape','input','electron_to_phonon','condensate_heat','transport')

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('runs',nargs=2,type=Path)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    if args.output.exists():raise ValueError('Do not overwrite prior evidence')
    coarse,fine=[original.load_run(p) for p in args.runs]
    contract.compare(coarse,fine)
    if any(r['record']['parameters']['method']!='ssprk3_common_flux_guarded' for r in (coarse,fine)):
        raise ValueError('Pair requires the same guarded SSP integrator')
    if fine['record']['parameters']['steps']!=2*coarse['record']['parameters']['steps']:
        raise ValueError('Fine step count must be twice the coarse count')
    criteria_path=contract.ROOT/'docs/implementation/stage2/acceptance_criteria.json'
    criteria=json.loads(criteria_path.read_text())
    if fine['record']['criteria_sha256']!=contract.sha(criteria_path):
        raise ValueError('Acceptance criteria changed')
    tolerance=criteria['coupled_trajectories']['finest_time_observable_relative_error_max']/4
    checks=[original.invariant_checks(r['record'],criteria) for r in (coarse,fine)]
    rows=[]
    for j,matches in original.shared_indices([coarse],fine):
        rows.append(dict(time=float(fine['times'][j]),
            errors=original.checkpoint_errors(coarse,matches[0],fine,j,FIELDS)))
    defined=all(v['comparison_defined'] for row in rows for v in row['errors'].values())
    maximum=max(v['assessment_error'] for row in rows for v in row['errors'].values()) if defined else None
    passed=defined and maximum<=tolerance and all(all(c.values()) for c in checks)
    result=dict(status='PASS_PAIR_DIAGNOSTIC' if passed else 'FAIL_PAIR_DIAGNOSTIC',
        scope='Supplemental fine-grid pair check only. Two levels cannot establish temporal convergence or close stage 2.',
        final_stage2_admission=False,three_level_time_admission=False,
        reviewer_sha256=contract.sha(__file__),contract_sha256=contract.sha(contract.__file__),
        criteria_sha256=contract.sha(criteria_path),
        runs=[dict(path=str(p),sha256=contract.sha(p)) for p in args.runs],
        invariant_checks=checks,checkpoints=rows,maximum_error=maximum,tolerance=tolerance,
        tolerance_policy='One quarter of the unchanged temporal criterion; an additional stricter diagnostic.',
        comparison_defined=defined)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({k:result[k] for k in ('status','maximum_error','tolerance')}))
    if not passed:raise SystemExit(1)

if __name__=='__main__':main()
