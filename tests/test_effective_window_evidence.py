import copy,importlib.util,unittest
from pathlib import Path
from sandbox.stage5_effective_validation import window_evidence as m
HASH="a"*64
ENTRY={"path":"/home/jdiaz/scratch/unit_fixture.json","sha256":"b"*64}
def primary():return {"event":{"confirmed":False,"confirmation_time_s":None},"final_time_ps":30.}
def window():
 old={"reference_saved_window_passed":False,"candidate_saved_window_passed":False}
 old.update({k:True for k in ["same_frozen_model_sources","same_nonrefined_parameters","actual_declared_refinement","voltage_waveform","electrical_latency_if_observed","temporal_refinement_at_least_two","common_physical_ROI_available","central_gap_magnitude_L2","central_deficit_over_reference_peak"]})
 original={"run_complete":False,"terminal_status":"RUNNING","saved_audit_passed":False,"completion_checks":{"complete_declared_window":False,"successful_saved_window_status":False},"noncompletion_failed_checks":[],"original_event":{"confirmed":False}}
 return {"schema":"pysnspd.explicit_comparison_window.v1","passed":True,"stage5_admission_inferred":False,"electrical_events_generated":False,"inputs_sha256":{"/home/jdiaz/scratch/unit_fixture.npz":"b"*64},"checks":{"reference_declared_interval_covered":True,"candidate_declared_interval_covered":True},"originals":{"reference":copy.deepcopy(original),"candidate":copy.deepcopy(original)},"original_comparison_result":{"kind":"mesh","passed":False,"checks":old},"declared_interval_ps":[0,20],"declaration":{"start_ps":0,"end_ps":20,"declared_utc":"2026-10-11T15:42:04Z"},"latency_certified_in_this_window":False}
def review():
 return {"schema":"pysnspd.finite_control_coverage_review.v1","manual_scientific_review":True,"windows":[{"window_sha256":HASH,"declared_interval_ps":[0,20],"decision":"COVERED","reviewed_utc":"PENDING_UNIT_FIXTURE","reviewed_by":"unit fixture","rationale":"Fixture only; no physics admission","evidence":[copy.deepcopy(ENTRY)],"critical_interval_end_ps":19.,"nominal_event_confirmation_time_ps":None,"nominal_saved_final_time_ps":30.,"latency_claimed":False}]}
class TestBridge(unittest.TestCase):
 def validate(self,w=None,r=None,p=None):return m.validate_window_evidence(w or window(),HASH,r or review(),verify_evidence=lambda e:None,primary_summary=p or primary())
 def test_separate_window_record_preserves_original_failed_comparison(self):
  out=self.validate();self.assertTrue(out["passed"]);self.assertFalse(out["original_full_run_comparison"]["passed"]);self.assertFalse(out["original_full_run_states"]["reference"]["run_complete"])
 def test_manual_review_required(self):
  with self.assertRaises(ValueError):self.validate(r={"schema":"wrong"})
 def test_window_hash_bound_review(self):
  r=review();r["windows"][0]["window_sha256"]="c"*64
  with self.assertRaises(ValueError):self.validate(r=r)
 def test_later_critical_region_requires_extension_or_bound(self):
  r=review();r["windows"][0]["critical_interval_end_ps"]=25
  with self.assertRaises(ValueError):self.validate(r=r)
 def test_later_detection_cannot_be_labeled_covered(self):
  r=review();r["windows"][0]["nominal_event_confirmation_time_ps"]=25
  p=primary();p["event"]={"confirmed":True,"confirmation_time_s":25e-12}
  with self.assertRaises(ValueError):self.validate(r=r,p=p)
 def test_latency_outside_coverage_rejected(self):
  r=review();r["windows"][0]["latency_claimed"]=True
  with self.assertRaises(ValueError):self.validate(r=r)
 def test_primary_event_must_match_actual(self):
  p=primary();p["event"]={"confirmed":True,"confirmation_time_s":25e-12}
  with self.assertRaises(ValueError):self.validate(p=p)
 def test_final_primary_window_must_match(self):
  p=primary();p["final_time_ps"]=40
  with self.assertRaises(ValueError):self.validate(p=p)
 def test_physical_failures_never_promoted(self):
  w=window();w["originals"]["reference"]["noncompletion_failed_checks"]=["energy"]
  with self.assertRaises(ValueError):self.validate(w=w)
 def test_source_failure_never_promoted(self):
  w=window();w["original_comparison_result"]["checks"]["same_frozen_model_sources"]=False
  with self.assertRaises(ValueError):self.validate(w=w)
 def test_missing_spatial_check_rejected(self):
  w=window();del w["original_comparison_result"]["checks"]["central_deficit_over_reference_peak"]
  with self.assertRaises(ValueError):self.validate(w=w)
 def test_unhashed_evidence_rejected(self):
  r=review();r["windows"][0]["evidence"]=["a_path"]
  with self.assertRaises(ValueError):self.validate(r=r)
 def test_quantified_bound_requires_every_gate(self):
  r=review();q=r["windows"][0];q.update(decision="QUANTIFIED_LATER_BOUND",critical_interval_end_ps=25,later_observable_bounds={})
  with self.assertRaises(ValueError):self.validate(r=r)
 def test_quantified_bound_does_not_relax_limits(self):
  r=review();q=r["windows"][0];q.update(decision="QUANTIFIED_LATER_BOUND",critical_interval_end_ps=25,later_observable_bounds={k:{"maximum_relative_error":.06,"declared_relative_limit":.05,"expression":"unit fixture","physical_argument":"unit fixture","evidence":[ENTRY]} for k in ["Vout","central_gap_L2","central_deficit"]})
  with self.assertRaises(ValueError):self.validate(r=r)
 def test_explicit_quantified_bound_is_separate_not_latency(self):
  r=review();q=r["windows"][0];q.update(decision="QUANTIFIED_LATER_BOUND",critical_interval_end_ps=25,later_observable_bounds={k:{"maximum_relative_error":.01,"declared_relative_limit":.05,"expression":"unit fixture","physical_argument":"unit fixture","evidence":[ENTRY]} for k in ["Vout","central_gap_L2","central_deficit"]})
  out=self.validate(r=r);self.assertTrue(out["passed"]);self.assertFalse(out["latency_certified_in_this_window"])
if __name__=="__main__":unittest.main()
