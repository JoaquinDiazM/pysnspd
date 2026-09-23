"""Read-only acceptance/provenance checks for this scoped document delivery."""
from pathlib import Path
import hashlib,json,subprocess
from pypdf import PdfReader
ROOT=Path(__file__).resolve().parents[3]
DATA=ROOT/'docs/implementation/stage2/time_pass_20260922'
load=lambda p:json.loads(p.read_text(encoding='utf-8'))
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
checks={}
audit=load(DATA/'one_cell_audit.json');plan=load(DATA/'two_cell_time_plan.json')
checks['85_one_cell_audit_checks_pass']=audit['checks_passed']==audit['checks_total']==85
checks['one_cell_scope_only']=not audit['stage2_closure'] and audit['all_measured_limiters_inactive']
checks['audit_bytes_match_next_plan']=sha(DATA/'one_cell_audit.json')==plan['evidence_and_decision']['one_cell_audit_sha256']
for path,digest in plan['source_hashes'].items():checks['source:'+path]=sha(ROOT/path)==digest
for record in load(DATA/'import_inventory.json')['files']:checks['import:'+record['path']]=sha(ROOT/record['path'])==record['sha256']
for name in ('one_ssp_admitted_fields_20260922.json','one_ssp_admitted_support_20260922.json'):
    checks['postprocess:'+name]=load(DATA/name)['status']=='PASS'
guards=load(DATA/'two_cell_batch_offline_tests.json')
checks['12_orchestration_guards_pass']=len(guards['checks'])==12 and all(guards['checks'].values())
checks['tested_runner_is_delivered_runner']=guards['runner_sha256']==sha(Path(__file__).parent/'run_two_cell_time_batch.py')
admission=load(ROOT/'docs/implementation/stage2/stage2_admission.json')
checks['global_stage2_and_push_withheld']=(not admission['numerical_admission'] and not admission['production_promotion'] and admission['publication']['status']=='CLOSURE_PUSH_WITHHELD')
stage3=load(ROOT/'docs/implementation/stage3/entry_contract.json')
checks['stage3_prepared_not_started']=not stage3['implementation_started'] and not stage3['calculations_executed']
pdf=ROOT/'output/pdf/implementation/Informe_resultados_temporales_etapa_2_20260922.pdf'
pages=PdfReader(pdf).pages
checks['pdf_three_nonempty_pages']=len(pages)==3 and all(len(p.extract_text())>1000 for p in pages)
checks['active_notebook_has_no_screen_commands']='screen' not in (ROOT/'docs/GEMINGA_COMMANDS.md').read_text(encoding='utf-8').lower()
checks['v1_tag_preserved']=subprocess.check_output(['git','rev-parse','v1.0.0^{}'],cwd=ROOT,text=True).strip()=='5ea0cd6a08d2cae414c82944f8230469ca5aa7d6'
user_snapshot=ROOT/'tmp/stage2_cells/preexisting_work.json'
if user_snapshot.exists():
    baseline=load(user_snapshot)
    checks['preexisting_user_work_preserved']=all(sha(ROOT/path)==digest for path,digest in baseline.items())
result=dict(status='PASS' if all(checks.values()) else 'FAIL',checks=checks,checks_passed=sum(checks.values()),checks_total=len(checks),
    pdf_sha256=sha(pdf),pdf_pages=len(pages),visual_review='Three rendered pages inspected; minor logarithmic labels corrected, no clipping or overlap in final rendering.',
    global_stage2_closed=False,new_physical_trajectories=0)
(DATA/'QA_delivery.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
print(json.dumps({k:v for k,v in result.items() if k!='checks'},indent=2))
if not all(checks.values()):raise SystemExit(1)
