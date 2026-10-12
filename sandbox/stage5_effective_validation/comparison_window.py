"""Proposed explicit finite-window comparison; no runtime or terminal-state mutation."""
from __future__ import annotations
import copy,hashlib,json,math
from pathlib import Path
import numpy as np
COMPLETION_CHECKS={"complete_declared_window","successful_saved_window_status"}
COMPARISON_COMPLETION_CHECKS={"reference_saved_window_passed","candidate_saved_window_passed"}
SAFE_RUN_STATUSES={"RUNNING","PROVISIONAL_LIVE_CHECKPOINT_VIEW","INTERRUPTED","EXPLICIT_TRANSFER_PREPARED","TRIGGER_MARGIN_COMPLETE","TIME_CENSORED_NO_TRIGGER","DEVELOPMENT_DURATION_REACHED","DC_HOLD_COMPLETED"}
SCHEMA="pysnspd.explicit_comparison_window.v1"
def passed_flag(value):return isinstance(value,(bool,np.bool_)) and bool(value)
def checked_interval(start_ps,end_ps):
 a,b=float(start_ps)*1e-12,float(end_ps)*1e-12
 if not math.isfinite(a+b) or a<0 or b<=a:raise ValueError("Finite positive comparison duration required")
 return a,b
def time_close(a,b):return abs(a-b)<=32*np.finfo(float).eps*max(abs(a),abs(b),1e-12)
def confirmed_event_in_window(event,a,b):
 return bool(event.get("confirmed") and event.get("crossing_time_s") is not None
  and event.get("confirmation_time_s") is not None
  and (event["crossing_time_s"]>=a or time_close(event["crossing_time_s"],a))
  and (event["confirmation_time_s"]<=b or time_close(event["confirmation_time_s"],b)))
def covers(t,a,b):
 t=np.asarray(t,float)
 return bool(t.ndim==1 and len(t)>1 and np.isfinite(t).all() and np.all(np.diff(t)>0)
  and (t[0]<=a or time_close(t[0],a)) and (t[-1]>=b or time_close(t[-1],b)))
def clipped_grid(t,a,b):
 t=np.asarray(t,float)
 if not covers(t,a,b):raise ValueError("Saved samples do not cover the declared comparison window")
 return np.unique(np.r_[a,t[(t>a)&(t<b)],b])
def comparison_view(run,a,b):
 # An analysis-only copy. The original manifest, summary and physical arrays are never changed.
 view=dict(run);summary=copy.deepcopy(run["summary"]);view["summary"]=summary
 tr=run["traces"];t=np.asarray(tr["time_s"],float);grid=clipped_grid(t,a,b)
 view["traces"]={k:np.interp(grid,t,np.asarray(v,float)) for k,v in tr.items()}
 view["traces"]["time_s"]=grid
 ht=np.array([r["time_s"] for r in run["hotbelt"]]);hg=clipped_grid(ht,a,b)
 gap=np.array([r["gap_min_meV"] for r in run["hotbelt"]])
 view["hotbelt"]=[dict(time_s=float(q),gap_min_meV=float(np.interp(q,ht,gap))) for q in hg]
 summary["final_time_ps"]=b*1e12
 if not confirmed_event_in_window(summary["event"],a,b):
  summary["event"]=dict(confirmed=False,latency_s=None,crossing_time_s=None,
   confirmation_time_s=None,reason="Original event absent or not fully covered by this analysis window")
 # Do not set passed=True: original run completion remains false and visible to the comparator.
 return view
def evaluate_window(kind,reference,candidate,*,start_ps,end_ps,compare_function,
                    declaration,critical_coverage_review=None):
 if declaration.get("schema")!="pysnspd.comparison_window_declaration.v1":
  raise ValueError("Explicit comparison-window declaration required")
 if declaration.get("start_ps")!=start_ps or declaration.get("end_ps")!=end_ps:
  raise ValueError("Requested interval differs from the pinned declaration")
 if not declaration.get("declared_utc") or not declaration.get("purpose"):
  raise ValueError("Declaration requires UTC and physical purpose")
 a,b=checked_interval(start_ps,end_ps)
 originals={}
 coverage={}
 for label,run in (("reference",reference),("candidate",candidate)):
  s=run["summary"];checks=s["checks"]
  if not COMPLETION_CHECKS<=set(checks):raise ValueError("Original completion checks must remain available")
  failure=[k for k,v in checks.items() if k not in COMPLETION_CHECKS and not passed_flag(v)]
  originals[label]=dict(run_complete=run["manifest"].get("run_complete"),
   terminal_status=run["manifest"].get("status"),saved_audit_passed=s["passed"],
   completion_checks={k:checks[k] for k in sorted(COMPLETION_CHECKS)},
   noncompletion_failed_checks=failure,original_event=copy.deepcopy(s["event"]))
  coverage[label]=covers(run["traces"]["time_s"],a,b) and covers([r["time_s"] for r in run["hotbelt"]],a,b)
 checks=dict(reference_original_status_safe=originals["reference"]["terminal_status"] in SAFE_RUN_STATUSES,
  candidate_original_status_safe=originals["candidate"]["terminal_status"] in SAFE_RUN_STATUSES,
  reference_physical_saved_checks_passed=not originals["reference"]["noncompletion_failed_checks"],
  candidate_physical_saved_checks_passed=not originals["candidate"]["noncompletion_failed_checks"],
  reference_declared_interval_covered=coverage["reference"],candidate_declared_interval_covered=coverage["candidate"])
 result=None;curves={};gapcurves={}
 if all(coverage.values()):
  rv=comparison_view(reference,a,b);cv=comparison_view(candidate,a,b)
  result,curves,gapcurves=compare_function(kind,rv,cv)
  # Explicit window result never rewrites the ordinary full-run audit result.
  checks.update({"comparison_"+k:v for k,v in result["checks"].items() if k not in COMPARISON_COMPLETION_CHECKS})
 latency_certified=all(confirmed_event_in_window(r["summary"]["event"],a,b) for r in (reference,candidate))
 review=copy.deepcopy(critical_coverage_review)
 # A numerical window pass is separate from the scientific choice that the critical interval is long enough.
 out=dict(schema=SCHEMA,status="EXPLICIT_FINITE_COMPARISON_WINDOW_AUDITED",stage5_admission_inferred=False,
  declared_interval_ps=[start_ps,end_ps],declaration=copy.deepcopy(declaration),
  originals=originals,coverage=coverage,checks=checks,passed=all(passed_flag(v) for v in checks.values()),
  original_comparison_result=result,latency_certified_in_this_window=latency_certified,
  electrical_events_generated=False,critical_coverage_review=review,
  wavelength_nm=reference["manifest"].get("parameter_identity",{}).get("wavelength_nm"),
  scope="Only declared finite interval. No original run_complete/status/controller is changed; no latency or censoring is assigned outside coverage. Extension is required if later critical dynamics/event is needed.",
  interpolation="Saved waveform and nodal-min scalars only for analysis at interval boundaries. Complex fields use the unchanged spatial auditor; no electrical event is interpolated or invented.")
 return out,curves,gapcurves

def audit_saved_window(kind,reference_dir,candidate_dir,declaration_path,output_file):
 from sandbox.stage5_effective_validation import audit_runs as audit
 output=Path(output_file).resolve()
 if not output.is_relative_to(Path("/home/jdiaz/scratch")):raise ValueError("Window output must remain in scratch")
 if output.exists():raise FileExistsError("Preserve existing window evidence")
 dp=Path(declaration_path).resolve();decl=json.loads(dp.read_text())
 files={}
 for path in (Path(reference_dir).resolve(),Path(candidate_dir).resolve()):
  for name in ("manifest.json","trajectory.npz","snapshots.npz","checkpoint_latest.npz"):
   p=path/name;files[str(p)]=audit.sha256(p)
 files[str(dp)]=audit.sha256(dp)
 reference=audit.audit_run("reference",reference_dir);candidate=audit.audit_run("candidate",candidate_dir)
 out,curves,gapcurves=evaluate_window(kind,reference,candidate,start_ps=decl["start_ps"],end_ps=decl["end_ps"],compare_function=audit.compare_runs,declaration=decl)
 if any(audit.sha256(p)!=h for p,h in files.items()):raise ValueError("Saved inputs advanced during window audit; use a coherent atomic capture")
 out.update(inputs_sha256=files,auditor_source_sha256=audit.sha256(audit.__file__),window_source_sha256=audit.sha256(__file__))
 output.parent.mkdir(parents=True,exist_ok=True)
 output.write_text(json.dumps(audit.json_ready(out),indent=2,allow_nan=False)+"\n")
 return out
