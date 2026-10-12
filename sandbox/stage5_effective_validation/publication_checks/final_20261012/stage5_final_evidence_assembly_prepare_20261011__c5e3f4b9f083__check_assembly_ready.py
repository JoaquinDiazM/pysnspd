"""Read-only final assembly preflight. Never audits, evolves or creates outputs."""
from pathlib import Path
import argparse,hashlib,json,math
HERE=Path(__file__).resolve().parent
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def read(p):return json.loads(Path(p).read_text())
def check(x,message):
 if not x:raise ValueError(message)
def live_owners(run):
 found=[]
 for p in Path("/proc").iterdir():
  if not p.name.isdigit():continue
  try:
   if p.stat().st_uid!=1011:continue
   argv=(p/"cmdline").read_bytes().split(b"\0")
  except OSError:continue
  if b"--output-dir" in argv and str(run).encode() in argv:found.append(int(p.name))
 return found
def readiness(plan):
 checks={};failures={}
 for name,run in plan["runs"].items():
  root=Path(run);owners=live_owners(root);mp=root/"manifest.json"
  if not mp.is_file():
   checks[name]=False;failures[name]={"reason":"EXPECTED_FINAL_RUN_NOT_CREATED"};continue
  m=read(mp);status=m.get("status");complete=m.get("run_complete") is True
  success=status in {"TRIGGER_MARGIN_COMPLETE","TIME_CENSORED_NO_TRIGGER","DEVELOPMENT_DURATION_REACHED","DC_HOLD_COMPLETED"}
  ok=not owners and complete and success
  if name in {"narrow","wide"}:
   ok=ok and status=="DEVELOPMENT_DURATION_REACHED" and m.get("parameter_identity",{}).get("duration_ps")==20.
  checks[name]=ok
  if not ok:failures[name]={"live_owner_PIDs":owners,"run_complete":m.get("run_complete"),"original_status":status,
   "reason":"CURRENT_RUN_ACTIVE_OR_INCOMPLETE_OR_UNSUCCESSFUL"}
 return checks,failures
def require_ready():
 plan=read(HERE/"assembly_plan.json")
 for path,digest in plan["source_sha256"].items():check(sha(path)==digest,"Approved analysis source changed:"+path)
 checks,failures=readiness(plan)
 check(all(checks.values()),"Final assembly blocked by original run states: "+json.dumps(failures,sort_keys=True))
 review_path=HERE/"root_saved_data_assembly_inputs_review.json"
 check(review_path.is_file(),"Final assembly blocked: actual inputs require root review; template grants no approval")
 review=read(review_path)
 check(review.get("approved_for_saved_data_audits") is True and review.get("stage5_admission") is False,
  "Only explicit saved-data assembly approval is accepted here; no automatic stage5 admission")
 check(review.get("assembly_plan_sha256")==sha(HERE/"assembly_plan.json") and
  review.get("commands_sha256")==sha(HERE/"run_final_saved_data_audits.PENDING.sh"),"Reviewed plan/commands changed")
 for entry in review.get("actual_input_files",[]):
  check(isinstance(entry,dict) and set(entry)=={"path","sha256"} and len(entry["sha256"])==64,
   "Every actual final input requires explicit path/SHA")
  check(sha(entry["path"])==entry["sha256"],"Final input advanced:"+entry["path"])
 required={str(Path(run)/f) for name,run in plan["runs"].items() if name!="dc"
  for f in ["manifest.json","trajectory.npz","snapshots.npz","checkpoint_latest.npz"]}
 actual={e["path"] for e in review.get("actual_input_files",[])}
 check(required<=actual,"Review omits exact completed original run files")
 evidence=Path(plan["expected_new_evidence_root"])
 protocol=evidence/"bounded_guard_protocol.json"
 check(protocol.is_file() and str(protocol) in actual,"Actual four-control guard protocol must be reviewed")
 p=read(protocol)
 check(set(p.get("controls",{}))=={"temporal","mesh","domain","spectral"} and
  p.get("schema")=="pysnspd.bounded_guard_comparison_protocol.v1","Actual four-control transfer proof required")
 declaration=Path(plan["final_common_window_declaration_path"])
 check(declaration.is_file() and str(declaration) in actual,"Final common window declaration must be selected and SHA-reviewed after completed controls")
 d=read(declaration)
 check(d.get("schema")=="pysnspd.comparison_window_declaration.v1" and d.get("start_ps")==0 and
  isinstance(d.get("end_ps"),(int,float)) and not isinstance(d.get("end_ps"),bool) and math.isfinite(d["end_ps"]) and d["end_ps"]>0 and d.get("declared_utc"),"Final common interval is still unassigned")
 check(not any((evidence/name).exists() for name in ["primary_saved_audits","temporal_comparison","mesh_comparison",
  "domain_comparison","spectral_comparison"]),"Existing final audit evidence must be preserved; choose a reviewed fresh evidence root")
 return {"ready_for_saved_data_audits_only":True,"stage5":False,"run_checks":checks}
if __name__=="__main__":
 p=argparse.ArgumentParser(description=__doc__);p.add_argument("--inspect",action="store_true");a=p.parse_args()
 if a.inspect:
  checks,failures=readiness(read(HERE/"assembly_plan.json"))
  print(json.dumps({"all_original_runs_complete_and_inactive":all(checks.values()),"run_checks":checks,
   "blockers":failures,"stage5":False},indent=2))
 else:print(json.dumps(require_ready(),indent=2))
