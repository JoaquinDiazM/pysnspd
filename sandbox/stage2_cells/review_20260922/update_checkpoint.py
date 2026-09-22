"""Record the reviewed failures without weakening gates or closing stage 2."""
from pathlib import Path
import json,hashlib
ROOT=Path(__file__).resolve().parents[3]
DATA=ROOT/'docs/implementation/stage2/review_20260922'
sha=lambda path:hashlib.sha256(Path(path).read_bytes()).hexdigest()
def main():
    target=ROOT/'docs/implementation/stage2/stage2_admission.json'
    prior=DATA/'admission_before_review.json'
    if prior.exists():raise ValueError('Prior checkpoint already archived; do not overwrite')
    prior.write_bytes(target.read_bytes())
    record=json.loads(target.read_text(encoding='utf-8'))
    audit=json.loads((DATA/'independent_audit.json').read_text(encoding='utf-8'))
    if audit['closure']!='NO' or audit['passed_checks']!=audit['checked_items']:
        raise ValueError('Unexpected audit contract')
    record.update(status='REVIEWED_NOT_CLOSED_PENDING_MANUAL_TIME',date='2026-09-22',
        stage2_status='NOT_CLOSED',numerical_admission=False,production_promotion=False,
        recommendation='Run only the new manual one-cell SSP temporal diagnostic. Retain failed RK4/SSP batches; do not restart them or proceed to spatial implementation.',
        generator_sha256=sha(__file__),
        preceding_checkpoint=dict(path=prior.relative_to(ROOT).as_posix(),sha256=sha(prior)),
        computation_handoff=dict(mode='USER_FOREGROUND_COMMAND',
            plan='docs/implementation/stage2/review_20260922/manual_time_plan.json',
            output_root='/home/jdiaz/pysnspd/tmp/stage2_ssp_time_20260922',
            expected_minutes=26,threads=1,memory_reservation_GiB=2,
            agent_executed=False,new_trajectories_this_review=0))
    record['previous_scoped_gates']=record['gates']
    record['gates']=[dict(id='retained_static_evidence',status='RETAINED_SCOPED_PASSES',
        detail='Prior static gates remain valid only in their declared scenarios; see previous_scoped_gates.'),
        *[dict(id='RK4_candidate_time_'+r['case'],status=r['status'],
            finest_error=r['finest_error'],scope=r['scope']) for r in audit['temporal_candidate_RK4']],
        dict(id='RK4_candidate_activity',status='PASS' if audit['activity_candidate_RK4']['passes'] else 'FAIL',
            detail=audit['activity_candidate_RK4']),*audit['failed_or_missing_gates']]
    record['manual_batch']['status']='SUPERSEDED_DO_NOT_RESTART'
    record['recovery']['status']='FAILED_SSP40_ENERGY_LEDGER'
    record['recovery']['final_admission']=False
    record['latest_review']=dict(audit='docs/implementation/stage2/review_20260922/independent_audit.json',
        audit_sha256=sha(DATA/'independent_audit.json'),report='output/pdf/implementation/Informe_revision_etapa_2_20260922.pdf',
        report_sha256=sha(ROOT/'output/pdf/implementation/Informe_revision_etapa_2_20260922.pdf'),
        current_plan_sha256=sha(DATA/'manual_time_plan.json'),physical_model_changed=False,thresholds_changed=False,
        final_report_eligible=False,report_scope='Results review; no closure certificate.')
    for name in ('import_inventory.json','independent_audit.json','diagnosed_ssp.json',
                 'manual_time_plan.json','manual_time_offline_tests.json','command_notebook_cleanup.json','QA_report.json'):
        record['evidence'].append(dict(path=(DATA/name).relative_to(ROOT).as_posix(),sha256=sha(DATA/name)))
    target.write_text(json.dumps(record,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(dict(status=record['status'],gates=len(record['gates']),new_physical_trajectories=0)))
if __name__=='__main__':main()
