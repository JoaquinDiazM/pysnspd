"""Finish documentary scope/provenance checks; never run a physical model."""
from pathlib import Path
import hashlib,json,subprocess
ROOT=Path(__file__).resolve().parents[3];DATA=ROOT/'docs/implementation/stage2/practical_review_20260922'
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
read=lambda p:json.loads(p.read_text(encoding='utf-8'))
write=lambda p,d:p.write_text(json.dumps(d,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
entrypath=ROOT/'docs/implementation/stage3/entry_contract.json';entry=read(entrypath)
entry['entry_requirements'][0].update(mandatory_for_static_stage3a=False,mandatory_for_coupled_dynamics=True)
entry['draft_spatial_acceptance']['later_stage3_gates_to_register']=[s.replace('Circuit storage and global energy balance D.28-D.33','Thesis three-state circuit, port power and storage balance CM.4/CM.8/CM.9') for s in entry['draft_spatial_acceptance']['later_stage3_gates_to_register']]
entry['continuous_contract']['circuit_override']['sha256']=sha(ROOT/entry['continuous_contract']['circuit_override']['path'])
write(entrypath,entry)
audit=read(DATA/'saved_results_audit.json')
adm_path=ROOT/'docs/implementation/stage2/stage2_admission.json';adm=read(adm_path)
for gate in adm['gates']:
    if gate['id']=='final_report_and_spatial_transition':gate.update(status='FULL_DYNAMIC_ADMISSION_PENDING_STATIC_PREPARATION_SEPARATE',
        reason='Full dynamic mesh certification is incomplete. A frozen-population static3A experiment uses its own relevant prerequisites; see scope_amendment.json.')
adm['gates']=[g for g in adm['gates'] if g['id']!='guarded_new_special_fields_and_support']
adm['gates'].append(dict(id='guarded_new_special_fields_and_support',status='PASS_REGISTERED_SPECIAL_CASES',
    evidence=audit['recorded_field_and_support_assessments'],scope='New equilibrium/vacuum/sparse cases only. Driven final trajectory field/upper-support checks were not reached in this batch.'))
for p in (DATA/'saved_results_audit.json',ROOT/'docs/modelo_v0_4/actualizaciones/circuito_memoria_20260922.md',DATA/'circuit_diagnostic.json'):
    relative=p.relative_to(ROOT).as_posix()
    adm['latest_evidence']=[r for r in adm['latest_evidence'] if r['path']!=relative]
    adm['latest_evidence'].append(dict(path=relative,sha256=sha(p)))
adm['latest_review'].update(report_status='VISUALLY_VERIFIED_FOUR_PAGES',circuits='THESIS_TOPOLOGY_DOCUMENTED_NOT_INTEGRATED')
write(adm_path,adm)
previous=read(ROOT/'tmp/stage2_cells/preexisting_work.json')
preserved=[dict(path=p,unchanged=sha(ROOT/p)==d) for p,d in previous.items()]
assert all(p['unchanged'] for p in preserved)
release=subprocess.check_output(['git','rev-parse','v1.0.0^{}'],cwd=ROOT,text=True).strip()
assert release=='5ea0cd6a08d2cae414c82944f8230469ca5aa7d6'
oldmanifest=read(ROOT/'tmp/stage2_closure_prep_20260922/sync_manifest.json')
old_sources={r['path']:r['sha256'] for r in oldmanifest['files']}
immutable=['docs/modelo_v0_4/C_condensado_corriente_y_senal_v0_4.md','docs/modelo_v0_4/D_sintesis_plan_y_verificaciones_v0_4.md','pysnspd/circuit/readout.py']
source_rows=[]
for p in immutable:
    git_bytes=subprocess.check_output(['git','show','HEAD:'+p],cwd=ROOT)
    assert git_bytes.replace(b'\r\n',b'\n')==(ROOT/p).read_bytes().replace(b'\r\n',b'\n'),p
    assert subprocess.run(['git','diff','--quiet','HEAD','--',p],cwd=ROOT).returncode==0,p
    source_rows.append(dict(path=p,sha256=sha(ROOT/p),matches_HEAD_after_line_ending_normalization=True))
write(DATA/'preservation_check.json',dict(status='PASS',release_commit=release,preexisting_files=preserved,unchanged_physical_sources=source_rows))
pdf=ROOT/'output/pdf/implementation/Informe_revision_criterios_y_circuito_20260922.pdf'
write(DATA/'report_visual_QA.json',dict(status='PASS_VISUAL_REVIEW',pages=4,reviewed_pages=[1,2,3,4],
    pdf_sha256=sha(pdf),pdf=pdf.relative_to(ROOT).as_posix(),
    checks=['No clipped text or tables','Circuit symbols and polarity reviewed','Different curves labelled and distinguishable','Nonnegative difference axis','Historical FAIL and unimplemented scope explicit'],
    source='All pages rendered with Poppler and visually inspected.'))
print(json.dumps(dict(status='METADATA_READY',preserved_user_files=len(preserved),historical_release_unchanged=True,
    new_detector_transients=0,new_long_runs=0,physical_production_changes=0)))
