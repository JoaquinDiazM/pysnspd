"""Proposed strict bridge from a finite comparison window to final reporting."""
from __future__ import annotations
import copy,math
from pathlib import Path
LIMITS={"temporal":.02,"spectral":.02,"mesh":.05,"domain":.05}
REQUIRED_OBSERVABLES={"temporal":{"Vout","gap_min"},"spectral":{"Vout","gap_min"},
 "mesh":{"Vout","central_gap_L2","central_deficit"},"domain":{"Vout","central_gap_L2","central_deficit"}}
COMPLETION={"reference_saved_window_passed","candidate_saved_window_passed"}
def positive_number(x):
 return not isinstance(x,bool) and isinstance(x,(float,int)) and math.isfinite(x) and x>=0
def validate_window_evidence(window,window_sha256,coverage_reviews,*,verify_evidence,primary_summary):
 if window.get("schema")!="pysnspd.explicit_comparison_window.v1":raise ValueError("Unsupported finite-window evidence schema")
 if window.get("passed") is not True or window.get("stage5_admission_inferred") is not False:
  raise ValueError("A window must pass its own explicit checks and infer no stage5 admission")
 if window.get("electrical_events_generated") is not False:raise ValueError("Window may never generate electrical events")
 inputs=window.get("inputs_sha256")
 if not isinstance(inputs,dict) or not inputs:raise ValueError("Stable SHA-bound saved window inputs required")
 for path,digest in inputs.items():
  if not isinstance(digest,str) or len(digest)!=64:raise ValueError("Window input digest missing")
  verify_evidence({"path":path,"sha256":digest})
 checks=window.get("checks")
 if not isinstance(checks,dict) or not checks or not all(v is True for v in checks.values()):raise ValueError("Finite-window checks failed or missing")
 originals=window.get("originals",{})
 for name in ("reference","candidate"):
  o=originals.get(name,{})
  if not {"run_complete","terminal_status","saved_audit_passed","completion_checks","noncompletion_failed_checks","original_event"}<=set(o):
   raise ValueError("Original full-run state must remain visible")
  if o["noncompletion_failed_checks"]:raise ValueError("A physical original failure cannot support a finite window")
  if set(o["completion_checks"])!={"complete_declared_window","successful_saved_window_status"}:
   raise ValueError("Only the two exact completion checks are separate")
 result=window.get("original_comparison_result")
 if not isinstance(result,dict) or result.get("kind") not in LIMITS:raise ValueError("Original numerical comparison required")
 if result.get("robust_spatial") is not None:
  from sandbox.stage5_effective_validation.guard_protocol_dispatcher import physical_gap_module
  physical_gap_module().validate_physical_gap_correction(result)
 rc=result.get("checks",{})
 required={"same_frozen_model_sources","same_nonrefined_parameters","actual_declared_refinement","voltage_waveform","electrical_latency_if_observed","temporal_refinement_at_least_two"}
 required|=({"common_physical_ROI_available","central_gap_magnitude_L2","central_deficit_over_reference_peak"} if result["kind"] in {"mesh","domain"} else {"gap_min_waveform"})
 later_limits={name:LIMITS[result["kind"]] for name in REQUIRED_OBSERVABLES[result["kind"]]}
 if "spectral_observable_policy" in result:
  if result["kind"]!="spectral":raise ValueError("Spectral policy cannot change another coverage axis")
  from sandbox.stage5_effective_validation.guard_protocol_dispatcher import validate_spectral_policy_application
  policy_limits=validate_spectral_policy_application(result)
  required|={"spectral_common_physical_ROI_available","spectral_central_gap_magnitude_L2","spectral_central_deficit_over_reference_peak"}
  if policy_limits["nodal_minimum_over_DC"] is None:
   required.discard("gap_min_waveform");later_limits.pop("gap_min",None)
  else:later_limits["gap_min"]=policy_limits["nodal_minimum_over_DC"]
  later_limits.update(central_gap_L2=policy_limits["central_field_L2_over_DC"],
   central_deficit=policy_limits["central_deficit_over_actual_reference_peak"])
 if not required<=set(rc):raise ValueError("Original numerical comparison omits mandatory checks")
 if not COMPLETION<=set(rc) or any(v is not True for k,v in rc.items() if k not in COMPLETION):
  raise ValueError("All original physical/observable/source checks remain mandatory")
 interval=window.get("declared_interval_ps")
 if not isinstance(interval,list) or len(interval)!=2 or not all(positive_number(v) for v in interval) or interval[1]<=interval[0]:
  raise ValueError("Invalid declared window")
 declaration=window.get("declaration",{})
 if declaration.get("start_ps")!=interval[0] or declaration.get("end_ps")!=interval[1] or not declaration.get("declared_utc"):
  raise ValueError("Window declaration missing or inconsistent")
 if not isinstance(coverage_reviews,dict) or coverage_reviews.get("schema")!="pysnspd.finite_control_coverage_review.v1" or coverage_reviews.get("manual_scientific_review") is not True:
  raise ValueError("Separate mandatory manual critical-coverage review required")
 matches=[r for r in coverage_reviews.get("windows",[]) if r.get("window_sha256")==window_sha256]
 if len(matches)!=1:raise ValueError("Exactly one SHA-bound manual review required")
 review=matches[0]
 if review.get("declared_interval_ps")!=interval or not review.get("reviewed_utc") or not review.get("reviewed_by") or not review.get("rationale"):
  raise ValueError("Incomplete explicit scientific coverage decision")
 evidence=review.get("evidence")
 if not isinstance(evidence,list) or not evidence:raise ValueError("Hashed physical coverage evidence required")
 for entry in evidence:
  if not isinstance(entry,dict) or not entry.get("path") or not isinstance(entry.get("sha256"),str) or len(entry["sha256"])!=64:raise ValueError("Every coverage evidence entry must be SHA-pinned")
  verify_evidence(entry)
 critical=review.get("critical_interval_end_ps")
 if not positive_number(critical):raise ValueError("Finite critical-region endpoint required")
 confirmation=review.get("nominal_event_confirmation_time_ps")
 pe=primary_summary["event"]
 actual_confirmation=None if not pe["confirmed"] else pe["confirmation_time_s"]*1e12
 if confirmation!=actual_confirmation:raise ValueError("Manual review event differs from final saved primary controller")
 if review.get("nominal_saved_final_time_ps")!=primary_summary["final_time_ps"]:raise ValueError("Coverage review must bind the final saved primary window")
 if critical>primary_summary["final_time_ps"]:raise ValueError("Critical coverage cannot exceed observed primary history")
 if confirmation is not None and not positive_number(confirmation):raise ValueError("Invalid saved nominal confirmation")
 needs_later=critical>interval[1] or (confirmation is not None and confirmation>interval[1])
 if review.get("latency_claimed") is True and window.get("latency_certified_in_this_window") is not True:
  raise ValueError("Latency cannot be certified outside both saved-event coverages")
 decision=review.get("decision")
 if decision=="COVERED":
  if needs_later:raise ValueError("Later critical dynamics or detection requires extension or a quantified bound")
 elif decision=="QUANTIFIED_LATER_BOUND":
  bounds=review.get("later_observable_bounds",{})
  if set(bounds)!=set(later_limits):raise ValueError("Every gated observable requires an explicit later bound")
  for name,bound in bounds.items():
   error=bound.get("maximum_relative_error");limit=bound.get("declared_relative_limit")
   if not positive_number(error) or limit!=later_limits[name] or error>limit or not bound.get("expression") or not bound.get("physical_argument"):
    raise ValueError("Missing or insufficient quantitative later bound: "+name)
   if not bound.get("evidence"):raise ValueError("Each later bound requires hashed evidence")
   for entry in bound["evidence"]:
    if not isinstance(entry,dict) or not entry.get("path") or not isinstance(entry.get("sha256"),str) or len(entry["sha256"])!=64:raise ValueError("Every later bound entry must be SHA-pinned")
    verify_evidence(entry)
 else:raise ValueError("Critical coverage is not admitted")
 # New record explicitly describes window semantics; the failed full-run result remains intact.
 row=copy.deepcopy(result)
 row.update(passed=True,checks=copy.deepcopy(checks),comparison_scope="EXPLICIT_FINITE_COMPARISON_WINDOW",
  declared_comparison_interval_ps=interval,original_full_run_comparison=copy.deepcopy(result),
  original_full_run_states=copy.deepcopy(originals),critical_coverage_review=copy.deepcopy(review),
  latency_certified_in_this_window=window["latency_certified_in_this_window"])
 return row
