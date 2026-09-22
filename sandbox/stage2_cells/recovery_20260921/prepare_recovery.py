"""Freeze reviewed recovery inputs; preserve the preceding admission checkpoint."""
from pathlib import Path
import json,sys
ROOT=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(Path(__file__).parent))
import assess_limited_time as contract
import run_limited_acceptance_batch as batch
DATA=ROOT/'docs/implementation/stage2/recovery_20260921'
def read(name):return json.loads((DATA/name).read_text())
def main():
    assessment=read('limited_short_time_assessment.json')
    if assessment['status']!='PASS':raise ValueError('Short temporal prerequisite failed')
    if assessment['integration_contract_comparison']['assessor_sha256']!=contract.sha(contract.__file__):
        raise ValueError('Short assessment predates the frozen assessor')
    for steps in (10,20,40):
        contract.verify_contract(contract.original.load_run(DATA/f'one_ssp_short_{steps}.json'))
    probe=contract.original.load_run(DATA/'limited_continuous_first_interval.json')
    contract.verify_contract(probe)
    criteria=json.loads((ROOT/contract.CRITERIA).read_text())
    if not all(contract.original.invariant_checks(probe['record'],criteria).values()):
        raise ValueError('Fine first-interval prerequisite failed')
    regression=read('regression_result.json')
    if regression['exit_code']!=0 or not regression['files_unchanged']:
        raise ValueError('Regression prerequisite failed')
    guards=read('limited_guard_tests.json')
    if guards['exit_code']!=0:raise ValueError('Guard tests failed')
    plan_path=DATA/'limited_acceptance_plan.json'
    plan=json.loads(plan_path.read_text())
    if plan['status']!='DRAFT_PENDING_SOURCE_FREEZE':raise ValueError('Plan already frozen; preserve it')
    plan['status']='READY_FOR_AUTHORIZED_SCREEN_RUN'
    plan['estimates']['wall_time_hours']='3-6, uncertain'
    plan['prerequisite_results']={name:contract.sha(DATA/name) for name in (
        'limited_short_time_assessment.json','limited_continuous_first_interval.json',
        'limited_continuous_first_interval.npz','limited_continuous_first_interval_initial.npz',
        'regression_result.json','limited_guard_tests.json','completed_batch_audit.json')}
    plan['source_hashes']={p:contract.sha(ROOT/p) for p in batch.required_contract_paths(plan)}
    batch.validate_plan(plan,True)
    old=ROOT/'docs/implementation/stage2/stage2_admission.json'
    prior=DATA/'admission_before_recovery.json'
    if prior.exists():raise ValueError('Do not overwrite prior admission snapshot')
    prior.write_bytes(old.read_bytes())
    plan_path.write_text(json.dumps(plan,indent=2)+'\n',encoding='utf-8')
    current=json.loads(old.read_text())
    current.update(status='READY_FOR_AUTHORIZED_SCREEN_RECOVERY',stage2_status='NOT_CLOSED',
        numerical_admission=False,production_promotion=False,
        recommendation='Run the authorized, source-frozen SSP recovery in code_000, detach, then review the resulting time/grid tests before admission.',
        generator_sha256=contract.sha(__file__),
        preceding_checkpoint=dict(path=prior.relative_to(ROOT).as_posix(),sha256=contract.sha(prior)),
        computation_handoff=dict(mode='EXPLICIT_AUTHORIZED_SCREEN_LAUNCH',screen='3040582.code_000',
            output_root='/home/jdiaz/pysnspd/tmp/stage2_recovery_20260921_ssp',
            wait_for_completion=False,polling=False),
        recovery=dict(plan=dict(path=plan_path.relative_to(ROOT).as_posix(),sha256=contract.sha(plan_path)),
            status='PREPARED_NOT_YET_FULLY_VALIDATED',physical_kernel_changes=False,
            temporal_integrator_changed=True,population_clipping=False,posthoc_energy_repair=False,
            old_RK4_time_tests='PASS on candidate mesh; not sufficient for new SSP integrator or fine mesh.',
            RK4_failure='First internal stage on electronic refinement4 leaves the physical population domain.',
            bisection_diagnostic='24 RHS calls, no accepted steps; archived, not selected for the relaunch.',
            short_SSP_time_error=assessment['finest_max_error'],
            fine_first_interval_ledger_error=probe['record']['energy_ledger_scaled_max'],
            final_admission=False))
    # Earlier detailed gates remain historical scoped evidence, not SSP admission.
    current['recovery']['inherited_gates_scope']='Retained prior checkpoint gates; every full SSP time/grid gate remains pending.'
    for name in (*plan['prerequisite_results'],'limited_acceptance_plan.json','command_notebook_update.json'):
        current['evidence'].append(dict(path=(DATA/name).relative_to(ROOT).as_posix(),sha256=contract.sha(DATA/name)))
    old.write_text(json.dumps(current,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(dict(status=current['status'],plan_sha256=contract.sha(plan_path),
        tasks=len(plan['tasks']),frozen_files=len(plan['source_hashes']))))
if __name__=='__main__':main()
