import copy,importlib.util,json,unittest
from pathlib import Path
import numpy as np
from sandbox.stage5_effective_validation import comparison_window as w
def run():
 checks={"complete_declared_window":False,"successful_saved_window_status":False,"energy_relative_defect":True,"same_sources":True,"positive_global_currents":True}
 return dict(manifest={"run_complete":False,"status":"RUNNING"},
 summary={"passed":False,"checks":checks,"final_time_ps":30.,"event":{"confirmed":False,"crossing_time_s":None,"confirmation_time_s":None,"latency_s":None}},
 traces={"time_s":np.array([0.,10.,20.,30.])*1e-12,"Vout_relative_V":np.array([0.,1.,2.,3.])*1e-6},
 hotbelt=[{"time_s":t*1e-12,"gap_min_meV":1.-t/100} for t in [0.,10.,20.,30.]])
def declaration(end=20):return dict(schema="pysnspd.comparison_window_declaration.v1",start_ps=0,end_ps=end,declared_utc="2026-10-11T16:00:00+00:00",purpose="Initial physical comparison; extend for critical dynamics")
def comparator(kind,r,c):
 return dict(kind=kind,passed=False,checks={"reference_saved_window_passed":r["summary"]["passed"],"candidate_saved_window_passed":c["summary"]["passed"],"same_frozen_model_sources":True,"voltage_waveform":True}),{},{}
class TestWindow(unittest.TestCase):
 def check(self,r=None,c=None,**kw):
  return w.evaluate_window("mesh",r or run(),c or run(),start_ps=0,end_ps=kw.get("end",20),compare_function=kw.get("compare",comparator),declaration=kw.get("declaration",declaration(kw.get("end",20))))[0]
 def test_partial_original_remains_visible(self):
  r=run();out=self.check(r=r);self.assertTrue(out["passed"]);self.assertFalse(out["originals"]["reference"]["run_complete"]);self.assertFalse(out["originals"]["reference"]["saved_audit_passed"]);self.assertFalse(r["summary"]["passed"]);self.assertEqual(r["manifest"]["status"],"RUNNING")
 def test_no_terminal_or_latency_invented(self):
  out=self.check();self.assertFalse(out["latency_certified_in_this_window"]);self.assertFalse(out["electrical_events_generated"]);self.assertFalse(out["stage5_admission_inferred"])
 def test_short_trajectory_pending(self):
  out=self.check(end=40);self.assertFalse(out["passed"]);self.assertIsNone(out["original_comparison_result"])
 def test_energy_failure_never_ignored(self):
  r=run();r["summary"]["checks"]["energy_relative_defect"]=False;out=self.check(r=r);self.assertFalse(out["passed"]);self.assertIn("energy_relative_defect",out["originals"]["reference"]["noncompletion_failed_checks"])
 def test_other_physical_failure_never_ignored(self):
  r=run();r["summary"]["checks"]["positive_global_currents"]=False;self.assertFalse(self.check(r=r)["passed"])
 def test_source_pin_comparator_failure_kept(self):
  def bad(k,r,c):
   result,_,_=comparator(k,r,c);result["checks"]["same_frozen_model_sources"]=False;return result,{},{}
  self.assertFalse(self.check(compare=bad)["passed"])
 def test_waveform_failure_kept(self):
  def bad(k,r,c):
   result,_,_=comparator(k,r,c);result["checks"]["voltage_waveform"]=False;return result,{},{}
  self.assertFalse(self.check(compare=bad)["passed"])
 def test_clipping_interpolates_analysis_not_event(self):
  r=run();view=w.comparison_view(r,0,17e-12);self.assertEqual(view["traces"]["time_s"][-1],17e-12);self.assertAlmostEqual(view["traces"]["Vout_relative_V"][-1],1.7e-6);self.assertFalse(view["summary"]["event"]["confirmed"]);self.assertEqual(r["summary"]["final_time_ps"],30.)
 def test_later_confirmed_event_outside_window(self):
  r=run();r["summary"]["event"].update(confirmed=True,crossing_time_s=25e-12,confirmation_time_s=26e-12,latency_s=25e-12)
  self.assertFalse(self.check(r=r,c=copy.deepcopy(r))["latency_certified_in_this_window"])
 def test_both_fully_covered_events_certified(self):
  r=run();r["summary"]["event"].update(confirmed=True,crossing_time_s=15e-12,confirmation_time_s=15.5e-12,latency_s=15e-12)
  self.assertTrue(self.check(r=r,c=copy.deepcopy(r))["latency_certified_in_this_window"])
 def test_confirmation_outside_window_not_certified(self):
  r=run();r["summary"]["event"].update(confirmed=True,crossing_time_s=19.8e-12,confirmation_time_s=20.3e-12,latency_s=19.8e-12)
  self.assertFalse(self.check(r=r,c=copy.deepcopy(r))["latency_certified_in_this_window"])
 def test_declaration_must_match(self):
  with self.assertRaises(ValueError):self.check(declaration=declaration(21))
 def test_missing_original_completion_rejected(self):
  r=run();del r["summary"]["checks"]["complete_declared_window"]
  with self.assertRaises(ValueError):self.check(r=r)
 def test_numpy_true_is_valid_boolean(self):
  r=run();r["summary"]["checks"]["energy_relative_defect"]=np.bool_(True);self.assertTrue(self.check(r=r)["passed"])
 def test_numpy_false_is_failure(self):
  r=run();r["summary"]["checks"]["energy_relative_defect"]=np.bool_(False);self.assertFalse(self.check(r=r)["passed"])
 def test_integer_true_is_not_physical_boolean(self):
  r=run();r["summary"]["checks"]["energy_relative_defect"]=1;self.assertFalse(self.check(r=r)["passed"])
 def test_numerical_failure_status_never_hidden(self):
  r=run();r["manifest"]["status"]="INTEGRATION_FAILURE";self.assertFalse(self.check(r=r)["passed"])
 def test_unknown_status_never_hidden(self):
  r=run();r["manifest"]["status"]="UNSPECIFIED";self.assertFalse(self.check(r=r)["passed"])
 def test_invalid_interval_rejected(self):
  for a,b in [(0,0),(-1,5),(0,float("inf"))]:
   with self.assertRaises(ValueError):w.checked_interval(a,b)
if __name__=="__main__":unittest.main()
