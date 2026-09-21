"""Record resumed evidence while retaining the preceding admission checkpoint."""
from pathlib import Path
import hashlib,json
ROOT=Path(__file__).resolve().parents[3]
DATA=ROOT/'docs/implementation/stage2'
RESUME=DATA/'resume_20260921'
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def read(name):return json.loads((RESUME/name).read_text())

def main():
    path=DATA/'stage2_admission.json'
    previous=RESUME/'admission_before_resume.json'
    if not previous.exists():
        assert sha(path)=='5ff42f809fcaa1593e897364037d49b2c71f15633547702969617dad1ab5b01e'
        previous.write_bytes(path.read_bytes())
    record=json.loads(previous.read_text())
    gates={g['id']:g for g in record['gates']}
    gates['reaction_continuum_IR_001'].update(status='PASS_SEPARATED_ELECTRONIC_REFINEMENT',
        detail='The archived combined errors remain nonmonotone. Independent matched-eta references and point-Bose comparisons separate fixed phonon interpolation, Pauli reconstruction, quadrature and regulator bias: all18 electronic series decrease over630/1260/2520. This historical check uses IR=.01; selected IR=.005 is separately verified at630.',
        separated_maximum_errors=[2.91333e-4,7.28419e-5,2.18741e-5])
    gates['selected_IR_0005_independent_reactions'].update(status='PASS_DECLARED_STATIC_CASES',
        detail='Independent continuous reference at actual IR=.005 and630/1025 for all three registered fields/profiles; full static responses meet1e-3. No dynamic admission follows.',
        maximum_weak_electron_error=2.9703e-4,maximum_power_error=7.3653e-4,
        maximum_phonon_hat_error=3.8811e-4)
    gates['simultaneous_trajectories_and_time_convergence'].update(status='PARTIAL_SHORT_REFERENCE_PASS_FULL_DURATION_PENDING',
        detail='User-run DOP853 complete on[0,.1]. RK4 10/20/40 matches at3 common times, max1.574710e-6 and reductions3.70/4.15. Full-duration one/two-cell time refinement remains pending.')
    gates['actual_trajectory_fields_and_support'].update(status='PASS_TWO_CELL_PILOT_FINAL_TRAJECTORIES_PENDING',
        detail='Actual630/1025 TWO40 trajectory:22 field samples PASS; omitted absorption bound6.67e-10. Accepted time/grid trajectories still need their own checks.')
    regression=read('regression_result.json')
    assert regression['exit_code']==0 and regression['files_unchanged']
    for p,s in regression['hashes_after'].items():assert sha(ROOT/p)==s,p
    gates['final_regression_suite'].update(status='PASS',detail='Fresh independent repository suite:509passed; source hashes unchanged.',
        evidence='docs/implementation/stage2/resume_20260921/regression_result.json')
    record.update(schema='pysnspd.stage2.admission_checkpoint.v2',status='PENDING_MANUAL_BATCH',
        stage2_status='NOT_CLOSED',numerical_admission=False,production_promotion=False,
        recommendation='Run the preregistered foreground RK4 batch. Review time/grid convergence, boundary cases and activity before issuing a final admission.',
        preceding_checkpoint=dict(path=previous.relative_to(ROOT).as_posix(),sha256=sha(previous)),
        resumed_reference=read('manual_reference_import.json'),
        manual_batch=dict(path='docs/implementation/stage2/resume_20260921/manual_validation_plan.json',
            sha256=sha(RESUME/'manual_validation_plan.json'),agent_executed=False,
            final_adjudication_required=True),
        kernel_changes=False,tolerances_relaxed=False)
    record['computation_handoff']['historical_only']=True
    record['computation_handoff']['resolution']='The user-run selected1025-node reference completed; preserved above. New pending work is the full-duration and mesh batch, not a repeat of that reference.'
    for p in ('manual_reference_import.json','manual_reference/one_reference_probe.json',
        'short_time_assessment.json','selected_continuum.json','electronic_three_grid_decomposition.json',
        'rhs_smoothness.json','two_active_40.json','two_active_fields.json','two_active_support.json',
        'assessment_guard_selfcheck.json','batch_guard_packaged_selftest.json','regression_result.json',
        'manual_validation_plan.json','command_notebook_update.json'):
        item=RESUME/p
        record['evidence'].append(dict(path=item.relative_to(ROOT).as_posix(),sha256=sha(item)))
    for item in record['evidence']:assert sha(ROOT/item['path'])==item['sha256'],item['path']
    record['generator_sha256']=sha(__file__)
    path.write_text(json.dumps(record,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(dict(status=record['status'],stage2_closed=False,evidence=len(record['evidence']),
        regression_passed=509,sha256=sha(path)),indent=2))

if __name__=='__main__':main()
