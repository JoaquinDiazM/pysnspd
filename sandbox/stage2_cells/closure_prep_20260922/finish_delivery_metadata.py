"""Consolidate reviewed artifact/provenance records without running dynamics."""
from pathlib import Path
import hashlib,json
ROOT=Path(__file__).resolve().parents[3];DATA=ROOT/'docs/implementation/stage2/closure_prep_20260922'
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
read=lambda p:json.loads(p.read_text(encoding='utf-8'))
write=lambda p,d:p.write_text(json.dumps(d,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
inventory=read(DATA/'isolated_guarded_inventory.json')
assert all(sha(ROOT/r['path'])==r['sha256'] for r in inventory['artifacts'])
reviewpath=DATA/'guarded_independent_review.json';review=read(reviewpath)
review.update(status='REVIEWED_SCOPED_FINDINGS_RESOLVED_LONG_VALIDATION_PENDING',resolutions=dict(
    reverse_reaction_coverage='isolated_guarded_negative_v2/isolated_results.json',
    isolated_result_provenance='isolated_guarded_inventory.json'))
write(reviewpath,review)
pdf=ROOT/'output/pdf/implementation/Informe_validacion_dos_celdas_etapa_2_20260922.pdf'
qa=dict(status='PASS_VISUAL_REVIEW',pdf=pdf.relative_to(ROOT).as_posix(),pdf_sha256=sha(pdf),pages=4,
    checked_pages=[1,2,3,4],method='Rendered each page with Poppler and visually inspected all four PNGs.',
    findings=dict(clipped_text=False,overlap=False,unreadable_glyphs=False,ambiguous_superposed_result_curves=False),
    figures=['two_time_convergence.png','two_cell_dynamics.png','guarded_rounding_and_escape.png'],
    scope='Legacy measured temporal/special results, diagnosed numerical defect, tested guarded prototype and explicit pending long validation.')
write(DATA/'report_visual_QA.json',qa)
statepath=DATA/'closure_status.json';state=read(statepath)
state['readiness_checks']=dict(guarded_microtests='PASS',reverse_reaction_microtests='PASS',external_inventory_files_verified=len(inventory['artifacts']),
    frozen_orchestrator_guards=22,real_entrypoint_smoke='PASS_ENTRYPOINT_CONTRACT_ONLY',pdf_visual_QA='PASS')
state['evidence']=[dict(path=p.relative_to(ROOT).as_posix(),sha256=sha(p)) for p in (
    DATA/'two_cell_audit.json',DATA/'isolated_guarded/isolated_results.json',DATA/'isolated_guarded_negative_v2/isolated_results.json',
    DATA/'isolated_guarded_inventory.json',DATA/'entrypoint_smoke/smoke_result.json',DATA/'guarded_frozen_offline_tests.json',
    DATA/'equilibrium_envelope_entrypoint_check.json',DATA/'report_visual_QA.json')]
write(statepath,state)
admpath=ROOT/'docs/implementation/stage2/stage2_admission.json';adm=read(admpath)
adm['latest_review']['report_status']='VISUALLY_VERIFIED_FOUR_PAGES'
for record in adm['latest_evidence']:record['sha256']=sha(ROOT/record['path'])
for p in (DATA/'guarded_independent_review.json',DATA/'isolated_guarded_inventory.json',DATA/'entrypoint_smoke/smoke_result.json',DATA/'report_visual_QA.json'):
    relative=p.relative_to(ROOT).as_posix()
    if not any(r['path']==relative for r in adm['latest_evidence']):adm['latest_evidence'].append(dict(path=relative,sha256=sha(p)))
write(admpath,adm)
print(json.dumps(dict(status='DELIVERY_METADATA_CONSOLIDATED',global_admission=adm['numerical_admission'],verified_isolated_artifacts=len(inventory['artifacts']),report_pages=4)))
