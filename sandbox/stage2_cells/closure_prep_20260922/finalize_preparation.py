"""Freeze reviewed inputs and record a pending manual batch, never admission."""
from pathlib import Path
import hashlib,json,sys
ROOT=Path(__file__).resolve().parents[3]
DATA=ROOT/'docs/implementation/stage2/closure_prep_20260922'
sys.path.insert(0,str(Path(__file__).parent))
import run_guarded_acceptance_batch as batch
read=lambda p:json.loads(p.read_text(encoding='utf-8'))
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
write=lambda p,d:p.write_text(json.dumps(d,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
main=read(DATA/'isolated_guarded/isolated_results.json')
supplement=read(DATA/'isolated_guarded_negative_v2/isolated_results.json')
assert main['status']=='PASS_GUARDED_PROTOTYPE_ONLY_NOT_GLOBALLY_ADMITTED'
assert supplement['status'].startswith('PASS')
assert all(main['checks'].values())
source=ROOT/'sandbox/stage2_cells/closure_prep_20260922/limited_ssp_guarded.py'
assert sha(source)=='9b77b91ff7b4e68786ea07a7e915f119d330f1c957ca22d2dd7b07028b00fad4'
planpath=DATA/'guarded_acceptance_plan.json';plan=read(planpath)
draft=DATA/'guarded_acceptance_plan_draft_before_review.json'
if not draft.exists():draft.write_bytes(planpath.read_bytes())
plan['source_hashes']={p:sha(ROOT/p) for p in batch.required_contract_paths(plan)}
plan['pending_source_freeze_paths']=[]
plan['pending_review']=dict(status='REVIEWED_PROTOTYPE_SCOPE_ONLY',items=[],
    evidence=['isolated_guarded/isolated_results.json','isolated_guarded_negative_v2/isolated_results.json','guarded_independent_review.json'],
    missing_long_evidence='All19 new trajectories and their acceptance assessments remain unexecuted.')
plan['status']='READY_FOR_USER_FOREGROUND_RUN'
write(planpath,plan);batch.validate_plan(plan,True)
state=dict(schema='pysnspd.stage2.closure-preparation.v1',date='2026-09-22',
    status='GUARDED_PROTOTYPE_TESTED_FULL_VALIDATION_PENDING',stage2_closed=False,
    numerical_admission=False,production_promotion=False,push_performed=False,
    summary='<b>La etapa 2 todavía no se puede cerrar.</b> Pasan las pruebas temporales recibidas y los controles especiales registrados. El fallo numérico adicional está localizado y corregido en un prototipo; falta acreditar sus trayectorias largas y mallas con una validación propia.',
    report_gate_rows=[
        ['SSP histórico: tiempo','PASS ONE/TWO','85/85 y 98/98 controles'],
        ['Equilibrio y fronteras','PASS histórico','3 trayectorias; energía < 2e-9'],
        ['Campos y soporte','PASS en muestras','Cota equilibrio: 1,44e-19'],
        ['Transporte aislado histórico','FAIL preservado','Redondeo subnormal; no se acepta'],
        ['Corrección aritmética','PASS prototipo','14 controles y reacciones inversas'],
        ['Nueva convergencia y mallas','PENDIENTE','19 trayectorias nuevas; 35 tareas'],
        ['Entrada a etapa 3','BLOQUEADA','Preparación documental disponible']],
    next_action='Ejecutar el nuevo lote completo con el integrador protegido. Repite equilibrio, fronteras y referencias temporales; después compara tres mallas fonónicas y tres electrónicas, con control temporal por pareja. Todas las tolerancias se conservan. El plan anterior de sólo mallas queda bloqueado para evitar validar una versión ya rechazada.',
    remaining_runtime_note='Estimación del lote manual: <b>9 a 15 horas, un hilo y reserva de 12 GiB</b>. El coste de la rama extendida en la malla más fina y el pico de memoria no están medidos; puede tardar más. La orden única queda en /home/jdiaz/GEMINGA_COMMANDS.md. No se lanzó ese cálculo.',
    legacy_two_audit=dict(checks=98,passed=98,scope='Frozen historical limited SSP only'),
    guard_sha256=sha(source),guard_checks=main['checks'],supplement_status=supplement['status'],
    manual_plan=dict(path=planpath.relative_to(ROOT).as_posix(),sha256=sha(planpath),tasks=35,new_trajectories=19,
                     status=plan['status'],executed=False),
    exclusions=['absolute NbN calibration','spatial dynamics','circuit','full SNSPD transient','production'],
    failed_evidence_preserved=['isolated serialization failure','isolated transport support failure',
        'equilibrium quadrature self-convergence failure','initial reverse-reaction fixture precheck failure'],
    criteria_changed=False,physical_kernels_changed=False,population_clipping=False,posthoc_energy_repair=False)
write(DATA/'closure_status.json',state)
admpath=ROOT/'docs/implementation/stage2/stage2_admission.json';adm=read(admpath)
adm['status']='HISTORICAL_TIME_PASS_GUARDED_PROTOTYPE_PASS_GLOBAL_STAGE2_PENDING'
adm['recommendation']='Run the source-frozen guarded validation plan manually. All method-specific long temporal and mesh gates remain pending; historical passes do not admit the new map. Stage3 is prepared, not started.'
for gate in adm['gates']:
    if gate['id']=='guarded_arithmetic_prototype':
        gate.update(status='PASS_PROTOTYPE_SCOPE_ONLY',scope='Fourteen checks and supplementary absorption/creation microcases pass; ordinary5step equivalence and complete isolated transport/escape are measured. No long moving-condensate or mesh admission.',evidence=['docs/implementation/stage2/closure_prep_20260922/isolated_guarded/isolated_results.json','docs/implementation/stage2/closure_prep_20260922/isolated_guarded_negative_v2/isolated_results.json'])
    if gate['id']=='guarded_complete_validation_plan':
        gate.update(status='READY_PENDING_USER_LONG_RUN',plan_status_at_update=plan['status'],sources_frozen=True,
                    scope='Fresh full validation prepared, not executed.19 new trajectories; no historical pass is transferred.')
adm['computation_handoff']=dict(mode='USER_FOREGROUND_COMMAND',plan=planpath.relative_to(ROOT).as_posix(),plan_sha256=sha(planpath),
    plan_status_at_update=plan['status'],sources_frozen=True,output_root='/home/jdiaz/pysnspd/tmp/stage2_guarded_acceptance_20260922',
    expected_hours=[9,15],threads=1,memory_reservation_GiB=12,memory_peak_measured=False,agent_executed=False,new_trajectories_pending=19,
    scope='Fresh guarded special, ONE/TWO time and true energy-mesh trajectories with all listed field/support checks; independent global review remains necessary.')
adm['latest_review'].update(status=state['status'],report='output/pdf/implementation/Informe_validacion_dos_celdas_etapa_2_20260922.pdf',
    report_status='RESULTS_REPORT_PREPARED_REQUIRES_VISUAL_QA',guarded_QA_status='PASS_PROTOTYPE_ONLY',guarded_plan_status=plan['status'])
for record in adm.get('latest_evidence',[]):
    p=ROOT/record['path']
    if p.exists():record['sha256']=sha(p)
for p in (DATA/'closure_status.json',planpath,DATA/'isolated_guarded/isolated_results.json',DATA/'isolated_guarded_negative_v2/isolated_results.json'):
    relative=p.relative_to(ROOT).as_posix()
    if not any(r['path']==relative for r in adm['latest_evidence']):adm['latest_evidence'].append(dict(path=relative,sha256=sha(p)))
write(admpath,adm)
entrypath=ROOT/'docs/implementation/stage3/entry_contract.json';entry=read(entrypath);snap=entry['current_stage2_snapshot']
snap['proposed_guarded_method']='PASS_PROTOTYPE_SCOPE_ONLY_LONG_VALIDATION_PENDING';snap['guarded_plan_status']=plan['status']
snap['remaining']=[s for s in snap['remaining'] if s!='Independent QA and source freeze of the guarded arithmetic map']
write(entrypath,entry)
print(json.dumps(dict(status=state['status'],plan_sha256=sha(planpath),frozen_files=len(plan['source_hashes']),long_calculation_launched=False)))
