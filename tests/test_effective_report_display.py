"""Display-only final-report inset and admission gate controls; no rendering."""
import importlib.util,json,tempfile,unittest
from pathlib import Path
from types import SimpleNamespace
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location("effective_final_report_display",ROOT/"sandbox/stage5_effective_validation/build_final_report.py")
report=importlib.util.module_from_spec(spec);spec.loader.exec_module(report)

class Label:
 def set_size(self,value):self.size=value
class Inset:
 def __init__(self):
  self.lines=[];self.settings={};self.xaxis=SimpleNamespace(label=Label())
  self.yaxis=SimpleNamespace(label=Label());self.title=Label()
 def plot(self,x,y,**kwargs):self.lines.append((np.array(x),np.array(y),kwargs))
 def set(self,**kwargs):self.settings.update(kwargs)
 def tick_params(self,**kwargs):self.tick_settings=kwargs
 def grid(self,**kwargs):self.grid_settings=kwargs
class Axes:
 def __init__(self):self.insets=[]
 def inset_axes(self,bounds):
  inset=Inset();self.insets.append((bounds,inset));return inset

class EffectiveReportDisplayTests(unittest.TestCase):
 def test_repository_receipt_keeps_full_arrays_in_scratch(self):
  scratch=Path("/home/jdiaz/scratch/stage5_pdf_compiler_repair_20261012")
  with tempfile.TemporaryDirectory(dir=scratch) as tmp:
   full=Path(tmp)/"full.json";full.write_text("{}")
   saved=report.repository_receipt({"pdf_sha256":"a"*64,"compiled":True,
    "input_sha256":{"input":"b"*64},"comparisons":{"large":[1]*10000},
    "primary_run_summaries":{"large":[2]*10000}},full)
   self.assertNotIn("comparisons",saved)
   self.assertNotIn("primary_run_summaries",saved)
   self.assertEqual(saved["pdf_sha256"],"a"*64)
   self.assertEqual(saved["full_numerical_receipt"]["sha256"],report.sha(full))
   self.assertLess(len(json.dumps(saved)),1000)
 def test_saved_numpy_scalar_receipt_preserves_types_and_rejects_nonfinite(self):
  payload={"checks":{"pass":np.bool_(True)},"count":np.int64(3),"value":np.float64(1.25)}
  saved=json.loads(json.dumps(payload,allow_nan=False,default=report.json_scalar))
  self.assertIs(saved["checks"]["pass"],True)
  self.assertEqual(saved,{"checks":{"pass":True},"count":3,"value":1.25})
  with self.assertRaises(ValueError):
   json.dumps({"bad":np.float32(np.nan)},allow_nan=False,default=report.json_scalar)
  with self.assertRaises(TypeError):
   json.dumps({"unsupported":object()},default=report.json_scalar)
 def runs(self,pulse_V=(0.,1e-6,2e-6),baseline_V=7e-6):
  return {"775":{"manifest":{"observation":{"baseline_Vout_V":baseline_V}},
   "history":{"time_s":np.array([0.,1e-12,2e-12]),"Vout_V":baseline_V+np.array(pulse_V)}}}
 def test_small_signal_inset_preserves_saved_values_and_subtracts_baseline(self):
  runs=self.runs();original=runs["775"]["history"]["Vout_V"].copy();ax=Axes()
  self.assertTrue(report.add_voltage_overview_inset(ax,runs,.1))
  x,y,_=ax.insets[0][1].lines[0]
  np.testing.assert_allclose(x,[0.,1.,2.],rtol=0,atol=1e-14)
  np.testing.assert_allclose(y,[0.,1.,2.],rtol=0,atol=1e-13)
  np.testing.assert_array_equal(runs["775"]["history"]["Vout_V"],original)
  self.assertIn("µV",ax.insets[0][1].settings["ylabel"])
 def test_signal_near_threshold_needs_no_inset(self):
  ax=Axes();self.assertFalse(report.add_voltage_overview_inset(ax,self.runs((0.,50e-6,90e-6)),.1))
  self.assertEqual(ax.insets,[])
 def test_equilibrium_needs_no_signal_inset(self):
  ax=Axes();self.assertFalse(report.add_voltage_overview_inset(ax,self.runs((0.,0.,0.)),.1))
  self.assertEqual(ax.insets,[])
 def test_positive_saved_threshold_required(self):
  for value in (0.,-1.,np.nan):
   with self.assertRaises(ValueError):report.add_voltage_overview_inset(Axes(),self.runs(),value)
 def test_stage5_false_gate_remains_blocked_without_output(self):
  scratch=Path("/home/jdiaz/scratch/stage5_external_audit_development_20261011")
  scratch.mkdir(exist_ok=True)
  with tempfile.TemporaryDirectory(dir=scratch) as tmp:
   path=Path(tmp)/"not_admitted.json";path.write_text(json.dumps({"stage5":False}))
   with self.assertRaisesRegex(ValueError,"stage5: true"):report.admission_inputs(path)
   with self.assertRaisesRegex(ValueError,"stage5: true"):
    report.main(["--admission",str(path),"--build-root",str(Path(tmp)/"forbidden_build"),"--compile","--cpu","4"])
   self.assertEqual(sorted(p.name for p in Path(tmp).iterdir()),["not_admitted.json"])
 def comparison(self,kind="mesh"):
  return dict(kind=kind,reference="nominal",candidate="control",relative_limit=.05 if kind in {"mesh","domain"} else .02,
   policy_version="spatial_observables_v2",
   Vout={"relative_error":.005},gap_min={"relative_error":.30,"diagnostic_only":kind in {"mesh","domain"},"former_gate_passed":False},
   robust_spatial={"L2_field":{"maximum_relative_error":.04},"deficit":{"relative_error_actual_peak":.01}})
 def test_mesh_uses_both_robust_metrics_preserving_failed_former_nodal_gate(self):
  result=report.comparison_admission_metrics(self.comparison())
  self.assertAlmostEqual(result["primary_gap_margin_ratio"],.8)
  self.assertAlmostEqual(result["deficit_margin_ratio"],.2)
  self.assertAlmostEqual(result["nodal_margin_ratio"],6.)
  self.assertTrue(result["nodal_diagnostic"])
  self.assertFalse(result["nodal_former_gate_passed"])
 def test_temporal_and_spectral_keep_nodal_metric(self):
  for kind in ("temporal","spectral"):
   result=report.comparison_admission_metrics(self.comparison(kind))
   self.assertAlmostEqual(result["primary_gap_margin_ratio"],15.)
   self.assertIsNone(result["deficit_margin_ratio"]);self.assertFalse(result["nodal_diagnostic"])
 def test_legacy_spatial_schema_cannot_be_silently_presented_as_new_policy(self):
  row=self.comparison();row["policy_version"]="legacy"
  with self.assertRaises(ValueError):report.comparison_admission_metrics(row)
  row=self.comparison();row["gap_min"]["diagnostic_only"]=False
  with self.assertRaises(ValueError):report.comparison_admission_metrics(row)
 def test_nonfinite_metric_or_margin_is_rejected(self):
  row=self.comparison();row["relative_limit"]=np.nan
  with self.assertRaises(ValueError):report.comparison_admission_metrics(row)
  row=self.comparison();row["robust_spatial"]["L2_field"]["maximum_relative_error"]=np.inf
  with self.assertRaises(ValueError):report.comparison_admission_metrics(row)
 def test_temporal_display_reads_actual_post_anchor_steps_and_guard_pair(self):
  row=self.comparison("temporal")
  row["declared_numerical_control"]={"reference_local_energy_tolerance":1e-5,"candidate_local_energy_tolerance":1.25e-6}
  row["accepted_step_statistics_after_both_guard_anchors"]=dict(available=True,interval_ps=[12.42,20.],
   reference={"available":True,"median_ps":.005},candidate={"available":True,"median_ps":.0025})
  result=report.temporal_refinement_display_rows([row])[0]
  self.assertEqual(result["reference_guard"],1e-5);self.assertEqual(result["candidate_guard"],1.25e-6)
  self.assertEqual(result["reference_accepted_steps"]["median_ps"],.005)
  self.assertEqual(result["candidate_accepted_steps"]["median_ps"],.0025)
  self.assertEqual(result["post_anchor_interval_ps"],[12.42,20.])
 def test_temporal_display_does_not_fabricate_post_anchor_window(self):
  row=self.comparison("temporal")
  row["accepted_step_statistics_after_both_guard_anchors"]={"available":False}
  result=report.temporal_refinement_display_rows([row])[0]
  self.assertIsNone(result["post_anchor_interval_ps"]);self.assertIsNone(result["reference_guard"])

 def finite_comparison(self,decision="COVERED",certified=False):
  row=self.comparison("domain")
  row.update(comparison_scope="EXPLICIT_FINITE_COMPARISON_WINDOW",
   declared_comparison_interval_ps=[0.,20.],critical_coverage_review={"decision":decision},
   latency_certified_in_this_window=certified,
   original_full_run_states={"reference":{"run_complete":True,"terminal_status":"CENSORED_MAX_TIME"},
    "candidate":{"run_complete":False,"terminal_status":"DEVELOPMENT_DURATION_REACHED"}})
  return row
 def test_finite_window_display_keeps_original_completion_states(self):
  import copy
  row=self.finite_comparison();before=copy.deepcopy(row)
  displayed=report.finite_window_display_rows([row])[0]
  self.assertEqual(displayed["interval_ps"],[0.,20.]);self.assertEqual(displayed["decision"],"COVERED")
  self.assertFalse(displayed["latency_certified"]);self.assertEqual(row,before)
  self.assertTrue(displayed["original_run_states"]["reference"]["run_complete"])
  self.assertFalse(displayed["original_run_states"]["candidate"]["run_complete"])
 def test_quantified_later_bound_display_never_invents_latency(self):
  row=self.finite_comparison("QUANTIFIED_LATER_BOUND")
  text=report.finite_window_scope_latex([row])
  self.assertIn("BOUND & no",text);self.assertIn("0--20",text)
  self.assertIn("EXPLICIT",text);self.assertIn("run\\_complete",text)
  self.assertIn("sin alteración",text)
 def test_finite_window_display_can_show_actual_latency_certification(self):
  self.assertIn("COVERED & sí",report.finite_window_scope_latex([self.finite_comparison(certified=True)]))
 def test_full_run_comparison_has_no_finite_window_table(self):
  self.assertEqual(report.finite_window_display_rows([self.comparison()]),[])
  self.assertEqual(report.finite_window_scope_latex([self.comparison()]),"")
 def test_finite_window_display_rejects_missing_scope_metadata(self):
  for field in ("declared_comparison_interval_ps","critical_coverage_review","original_full_run_states","latency_certified_in_this_window"):
   row=self.finite_comparison();row.pop(field)
   with self.assertRaises(ValueError):report.finite_window_display_rows([row])
  row=self.finite_comparison("PENDING")
  with self.assertRaises(ValueError):report.finite_window_display_rows([row])
  row=self.finite_comparison();row["declared_comparison_interval_ps"]=[0.,float("nan")]
  with self.assertRaises(ValueError):report.finite_window_display_rows([row])
  row=self.finite_comparison();row["latency_certified_in_this_window"]=1
  with self.assertRaises(ValueError):report.finite_window_display_rows([row])

if __name__=="__main__":unittest.main()
