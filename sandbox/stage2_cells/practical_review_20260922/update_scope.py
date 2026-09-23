"""Record the new workflow scope without changing old acceptance results."""
from pathlib import Path
import hashlib,json
ROOT=Path(__file__).resolve().parents[3];DATA=ROOT/'docs/implementation/stage2/practical_review_20260922'
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
read=lambda p:json.loads(p.read_text(encoding='utf-8'))
write=lambda p,d:p.write_text(json.dumps(d,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
admpath=ROOT/'docs/implementation/stage2/stage2_admission.json'
entrypath=ROOT/'docs/implementation/stage3/entry_contract.json'
for p in (admpath,entrypath):
    archive=DATA/('before_'+p.name)
    if not archive.exists():archive.write_bytes(p.read_bytes())
pair=read(DATA/'raw/two_guarded_ph2049_pair_320_640.json')
assert pair['status']=='FAIL_PAIR_DIAGNOSTIC' and pair['tolerance']==2.5e-5
scope=dict(schema='pysnspd.validation.scope-amendment.v1',date='2026-09-22',
    basis='User requests proportional verification and restoration of the thesis readout circuit.',
    old_numerical_tolerances_changed=False,old_failures_reclassified=False,
    decision='The auxiliary pair budget is retained as a failed diagnostic, but does not block unrelated static spatial implementation/preparation.',
    essential_controls=['source/unit provenance','consistent physical signs and energy bookkeeping','valid state/population domain','applicable stationary/invariant checks'],
    static_stage3a=dict(allowed_scope='Small static spatial energy, Cartesian derivatives and current at fixed thermal/nonthermal populations.',
        prerequisites=['Own preregistration','R2 energy/force domain checks','D.36 evaluated on tested states','Independent static derivative/current checks'],
        requires_complete_dynamic_stage2_mesh_certificate=False,implementation_started=False,
        excludes=['kinetic time integration','photon','circuit transient','detector pulse','production']),
    future_dynamic_accuracy=dict(status='OBSERVABLE_BUDGET_TO_REGISTER_BEFORE_NEW_RESULTS',
        rule='Choose the physical difference to resolve and numerical uncertainty in those observables. Do not choose a new number merely to pass the existing pair.',
        old_certificate='Full original stage2 certificate remains incomplete; no full-dynamic admission is inferred.'),
    computation_policy=dict(new_long_batch_required_now=False,old_guarded_batch='STOPPED_AFTER21_COMPLETE_TASKS_AT_SUPPLEMENTARY_PAIR',
        reuse='Reevaluate compatible saved trajectories. Repeat dynamics only when a model/numerical change or unresolved decision requires it.'),
    failed_pair=dict(path=(DATA/'raw/two_guarded_ph2049_pair_320_640.json').relative_to(ROOT).as_posix(),sha256=sha(DATA/'raw/two_guarded_ph2049_pair_320_640.json'),
        maximum_difference=pair['maximum_error'],original_auxiliary_budget=pair['tolerance'],status=pair['status']))
write(DATA/'scope_amendment.json',scope)
adm=read(admpath);adm.update(status='GUARDED_TIME_VALIDATED_AUXILIARY_PAIR_FAILED_STATIC_WORK_SEPARATELY_ALLOWED',
    numerical_admission=False,production_promotion=False,stage2_status='STRICT_DYNAMIC_MESH_CERTIFICATE_INCOMPLETE',
    recommendation='Do not repeat the35-task batch automatically. Preserve its21 completed tasks and failed auxiliary pair. Prepare/implement the separately registered static stage3A scope; register useful dynamic observable accuracy before further costly certification.')
adm['candidate']['proposed_method_status']='GUARDED_TIME_VALIDATED_CANDIDATE_GRID_MESH_CERTIFICATE_INCOMPLETE'
adm['computation_handoff']=dict(mode='NO_NEW_LONG_RUN_REQUIRED_NOW',previous_batch='tmp/stage2_guarded_acceptance_20260922',
    status='PARTIALLY_COMPLETED_SUPERSEDED_AS_UNIVERSAL_WORKFLOW_GATE',executed_by='user',agent_launched_long_run=False,
    preserve_original_plan=True,active_execute_command=None,
    reason='The stopped assessment is supplementary precision separation; it is not a prerequisite of frozen-population static spatial work.')
for gate in adm['gates']:
    if gate['id']=='guarded_complete_validation_plan':gate.update(status='USER_EXECUTED_PARTIALLY_AUXILIARY_PAIR_FAILED',
        sources_frozen=True,completed_tasks=21,scope='Original full certificate remains incomplete; no blanket block on unrelated static stage3A.')
adm['gates'] += [dict(id='guarded_new_ONE_TWO_time',status='PASS_REGISTERED_CANDIDATE_GRID',
    scope='New guarded160/320/640 vs1280 trajectories on630/1025, no inference to missing fine grids.',
    evidence=['docs/implementation/stage2/practical_review_20260922/raw/one_guarded_time_assessment.json','docs/implementation/stage2/practical_review_20260922/raw/two_guarded_time_assessment.json']),
    dict(id='guarded_ph2049_supplementary_pair',**scope['failed_pair']),
    dict(id='static_stage3a_scope',status='SEPARATE_PREPARATION_ALLOWED_NOT_IMPLEMENTED',scope=scope['static_stage3a'])]
adm['latest_review']=dict(directory=DATA.relative_to(ROOT).as_posix(),status=adm['status'],
    report='output/pdf/implementation/Informe_revision_criterios_y_circuito_20260922.pdf',global_stage2_closed=False,
    scope_amendment='docs/implementation/stage2/practical_review_20260922/scope_amendment.json',
    raw_batch='docs/implementation/stage2/practical_review_20260922/raw')
adm['latest_evidence']=[dict(path=p.relative_to(ROOT).as_posix(),sha256=sha(p)) for p in (
    DATA/'raw/batch_manifest.json',DATA/'raw/batch_events.jsonl',DATA/'raw/two_guarded_ph2049_pair_320_640.json',
    DATA/'scope_amendment.json',DATA/'validation_scope_review.md')]
adm['previous_checkpoint']=dict(path=(DATA/'before_stage2_admission.json').relative_to(ROOT).as_posix(),sha256=sha(DATA/'before_stage2_admission.json'))
write(admpath,adm)
entry=read(entrypath);entry['execution_status']='STATIC_STAGE3A_SCOPE_PREPARED_DYNAMIC_STAGE3_BLOCKED'
entry['schema']='pysnspd.stage3.entry_contract.draft.v2'
entry['purpose']='Prepare a static spatial3A experiment independently from full dynamic mesh certification; later coupled dynamics retains its own prerequisites.'
entry['continuous_contract']['circuit_override']=dict(path='docs/modelo_v0_4/actualizaciones/circuito_memoria_20260922.md',
    equations=['CM.1-CM.11'],supersedes=['C.32-C.34','D.28-D.29','D.33','circuit initialization in D.30','I becomes Is in D.26-D.27'],
    implementation_started=False,topology='Thesis4.16 voltage bias, Rbias/Lbias, detector branch, coupling capacitor and readout load.',
    state=['Ib','Is','vc'])
entry['current_stage2_snapshot'].update(guarded_new_ONE_TWO='PASS_REGISTERED_CANDIDATE_TIME',
    guarded_ph2049_pair='FAIL_SUPPLEMENTARY_PRECISION_BUDGET_NOT_PHYSICAL_SUPPORT',
    guarded_plan_status='EXECUTED_PARTIALLY_NO_AUTOMATIC_RERUN',proposed_guarded_method='TIME_VALIDATED_CANDIDATE_MESH_ONLY',
    global_numerical_admission=False,
    remaining=['Complete dynamic mesh certificate if that strict numerical claim is sought','Own3A static registration and applicable checks before3A results','Observable-based dynamic uncertainty budget before subsequent detector-oriented tests'],
    scope_rule='Historical/current dynamic passes retain their scopes. No dynamic certificate is inferred; static3A has separate dependencies.')
entry['entry_requirements'][0].update(applies_to='Coupled dynamic stage3 and any production promotion; not an independent frozen-population static3A experiment.',
    status='BLOCKED_FOR_COUPLED_DYNAMICS')
for item in entry['entry_requirements'][1:]:
    if item['id'] in ('configuration_and_provenance','inherited_numerical_contract'):
        item['status']='REQUIRED_IN_APPLICABLE_STATIC_OR_DYNAMIC_SCOPE'
entry['first_bounded_work']['status']='PREPARED_STATIC_SCOPE_REQUIRES_OWN_REGISTRATION_NOT_RUN'
entry['static_stage3a_scope']=scope['static_stage3a']
entry['conclusion_limit']='Static3A preparation is allowed with its own relevant checks. No spatial implementation, photon/circuit trajectory or production promotion occurred in this update.'
write(entrypath,entry)
print('Updated scope; historical tolerance and failure preserved, no global numerical admission.')
