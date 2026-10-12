"""Spatial saved-field and explicit local-energy ledger audit controls."""
import copy,importlib.util,json,math,os,tempfile,unittest
from pathlib import Path
import numpy as np
ROOT=Path(os.environ.get("PYSNSPD_TEST_ROOT",Path(__file__).resolve().parents[1]))
def import_path(name,path):
 spec=importlib.util.spec_from_file_location(name,path)
 module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module);return module
old=import_path("legacy_audit_tests",ROOT/"tests/test_effective_saved_audit.py")
audit=import_path("candidate_audit",Path(os.environ.get("AUDIT_CANDIDATE",ROOT/"sandbox/stage5_effective_validation/audit_runs.py")))
old.audit=audit

class SpatialAuditTests(old.EffectiveSavedAuditTests):
 def fixture(self,path):
  manifest,hist,snap=super().fixture(path)
  manifest["source"]["center_m"]=[80e-9,0.]
  (path/"manifest.json").write_text(json.dumps(manifest))
  return manifest,hist,snap
 def pair(self,path,kind="mesh"):
  self.fixture(path);a=audit.audit_run("ref",path);b=copy.deepcopy(a)
  b["summary"]["name"]="candidate"
  keys={"mesh":("spacing_nm",2.5),"domain":("length_nm",320.),
   "temporal":("dt_ps",.005),"spectral":("xi_order",8)}
  key,value=keys[kind];b["manifest"]["parameter_identity"][key]=value
  return a,b
 def test_mesh_nodal_gate_is_visible_diagnostic_not_false_relabel(self):
  with self.workspace() as tmp:
   a,b=self.pair(Path(tmp))
   b["hotbelt"][-1]["gap_min_meV"]=.6
   result,_,_=audit.compare_runs("mesh",a,b)
   self.assertTrue(result["passed"])
   self.assertFalse(result["legacy_nodal_policy_passed"])
   self.assertFalse(result["gap_min"]["former_gate_passed"])
   self.assertTrue(result["gap_min"]["diagnostic_only"])
   self.assertNotIn("gap_min_waveform",result["checks"])
 def test_same_mesh_temporal_and_spectral_keep_nodal_gate(self):
  with self.workspace() as tmp:
   for kind in ("temporal","spectral"):
    a,b=self.pair(Path(tmp),kind);b["hotbelt"][-1]["gap_min_meV"]=.6
    result,_,_=audit.compare_runs(kind,a,b)
    self.assertFalse(result["passed"]);self.assertFalse(result["checks"]["gap_min_waveform"])
    self.assertFalse(result["gap_min"]["diagnostic_only"])
 def test_large_field_difference_fails_mesh_even_identical_nodal_metrics(self):
  with self.workspace() as tmp:
   a,b=self.pair(Path(tmp))
   b["snap"]["delta_J"][1:]*=.5
   b["spatial"]=audit.spatial_snapshot_diagnostics(b["manifest"],b["snap"])
   result,_,_=audit.compare_runs("mesh",a,b)
   self.assertFalse(result["checks"]["central_gap_magnitude_L2"])
   self.assertFalse(result["passed"])
 def test_deficit_uses_actual_peak_not_DC_area_to_dilute_signal(self):
  with self.workspace() as tmp:
   a,b=self.pair(Path(tmp))
   a["spatial"]["deficit_J_m2"]=np.array([0.,1.,1.])*1e-40
   b["spatial"]["deficit_J_m2"]=np.array([0.,1.1,1.1])*1e-40
   result,_,_=audit.compare_runs("mesh",a,b)
   self.assertAlmostEqual(result["robust_spatial"]["deficit"]["relative_error_actual_peak"],.1)
   self.assertLess(result["robust_spatial"]["deficit"]["relative_error_DCarea"],.05)
   self.assertFalse(result["checks"]["central_deficit_over_reference_peak"])
 def test_zero_deficit_reference_has_no_hidden_normalization_floor(self):
  with self.workspace() as tmp:
   a,b=self.pair(Path(tmp))
   a["spatial"]["deficit_J_m2"][:]=0
   b["spatial"]["deficit_J_m2"][:]=1e-45
   result,_,_=audit.compare_runs("domain",a,b)
   self.assertTrue(math.isinf(result["robust_spatial"]["deficit"]["relative_error_actual_peak"]))
   self.assertFalse(result["checks"]["central_deficit_over_reference_peak"])
 def test_missing_source_center_cannot_fabricate_common_physical_ROI(self):
  with self.workspace() as tmp:
   a,b=self.pair(Path(tmp));b["spatial"]=None
   result,_,_=audit.compare_runs("mesh",a,b)
   self.assertFalse(result["checks"]["common_physical_ROI_available"])
   self.assertFalse(result["passed"])
 def test_strict_sources_and_both_numerical_tolerances_remain_strict(self):
  with self.workspace() as tmp:
   a,b=self.pair(Path(tmp))
   for key in ("step_tolerance","local_energy_tolerance"):
    changed=copy.deepcopy(b);changed["manifest"]["parameter_identity"][key]=1e-4
    result,_,_=audit.compare_runs("mesh",a,changed)
    self.assertFalse(result["checks"]["same_nonrefined_parameters"])
    self.assertIn(key,result["forbidden_parameter_changes"])
   b["manifest"]["source_pins_sha256"]["solver"]="different"
   result,_,_=audit.compare_runs("mesh",a,b)
   self.assertFalse(result["checks"]["same_frozen_model_sources"])
 def test_energy_failure_remains_failure_under_new_spatial_observables(self):
  with self.workspace() as tmp:
   a,b=self.pair(Path(tmp));b["summary"]["passed"]=False
   result,_,_=audit.compare_runs("mesh",a,b)
   self.assertFalse(result["checks"]["candidate_saved_window_passed"])
   self.assertFalse(result["passed"])
 def test_complex_interpolation_precedes_magnitude(self):
  with self.workspace() as tmp:
   a,b=self.pair(Path(tmp))
   b["snap"]["delta_J"][:,0]=-b["snap"]["delta_J"][:,0]
   b["spatial"]=audit.spatial_snapshot_diagnostics(b["manifest"],b["snap"])
   metric,_=audit.compare_robust_spatial(a,b,np.array([0.,1e-12]))
   self.assertGreater(metric["L2_field"]["maximum_relative_error"],0.)
 def test_common_physical_x_translation_leaves_comparison_unchanged(self):
  with self.workspace() as tmp:
   a,b=self.pair(Path(tmp),"domain")
   b["snap"]["coordinates_m"][:,0]+=80e-9;b["manifest"]["source"]["center_m"][0]+=80e-9
   b["spatial"]=audit.spatial_snapshot_diagnostics(b["manifest"],b["snap"])
   metric,_=audit.compare_robust_spatial(a,b,np.array([0.,1e-12]))
   self.assertLess(metric["L2_field"]["maximum_relative_error"],1e-14)
   self.assertLess(metric["deficit"]["relative_error_actual_peak"],1e-14)

 def test_only_declared_temporal_dt_guard_pair_gets_exception(self):
  with self.workspace() as tmp:
   a,b=self.pair(Path(tmp),"temporal")
   a["manifest"]["parameter_identity"].update(dt_ps=.02,local_energy_tolerance=1e-5,step_tolerance=2e-7)
   b["manifest"]["parameter_identity"].update(dt_ps=.01,local_energy_tolerance=1.25e-6,step_tolerance=2e-7)
   result,_,_=audit.compare_runs("temporal",a,b)
   self.assertTrue(result["passed"])
   self.assertIsNotNone(result["declared_numerical_control"])
   self.assertIn("local_energy_tolerance",result["changed_parameters"])
   self.assertNotIn("local_energy_tolerance",result["forbidden_parameter_changes"])
 def test_temporal_exception_does_not_generalize_wrong_guard_cap_NL_or_axis(self):
  with self.workspace() as tmp:
   a,b=self.pair(Path(tmp),"temporal")
   a["manifest"]["parameter_identity"].update(dt_ps=.02,local_energy_tolerance=1e-5,step_tolerance=2e-7)
   b["manifest"]["parameter_identity"].update(dt_ps=.01,local_energy_tolerance=1.25e-6,step_tolerance=2e-7)
   for key,value in (("local_energy_tolerance",2e-6),("dt_ps",.005),("step_tolerance",2e-8)):
    changed=copy.deepcopy(b);changed["manifest"]["parameter_identity"][key]=value
    result,_,_=audit.compare_runs("temporal",a,changed)
    self.assertFalse(result["passed"]);self.assertIsNone(result["declared_numerical_control"])
    self.assertIn("local_energy_tolerance",result["forbidden_parameter_changes"])
   for axis in ("mesh","domain","spectral"):
    result,_,_=audit.compare_runs(axis,a,b)
    self.assertIsNone(result["declared_numerical_control"])
    self.assertFalse(result["checks"]["same_nonrefined_parameters"])
 def test_actual_accepted_steps_reported_without_inventing_refinement(self):
  with self.workspace() as tmp:
   a,b=self.pair(Path(tmp),"temporal")
   a["history"]["accepted_step_ps"]=np.array([0.,.002,.001])
   b["history"]["accepted_step_ps"]=np.array([0.,.002,.001])
   result,_,_=audit.compare_runs("temporal",a,b)
   self.assertEqual(result["accepted_step_statistics"]["reference"]["median_ps"],.0015)
   self.assertEqual(result["accepted_step_statistics"]["candidate"]["median_ps"],.0015)

 def test_post_anchor_steps_do_not_hide_refinement_under_legacy_prefix(self):
  with self.workspace() as tmp:
   a,b=self.pair(Path(tmp),"temporal")
   a["history"]["accepted_step_ps"]=np.array([0.,.1,.002])
   b["history"]["accepted_step_ps"]=np.array([0.,.1,.001])
   for run,anchor in ((a,.5),(b,1.)):
    run["summary"]["local_energy_guard"]=dict(applicable=True,guard_active_since_time_ps=anchor)
   result,_,_=audit.compare_runs("temporal",a,b)
   stats=result["accepted_step_statistics_after_both_guard_anchors"]
   self.assertTrue(stats["available"]);self.assertEqual(stats["interval_ps"],[1.,2.])
   self.assertEqual(stats["reference"]["median_ps"],.002)
   self.assertEqual(stats["candidate"]["median_ps"],.001)

 def test_roundoff_only_post_anchor_window_is_not_available(self):
  with self.workspace() as tmp:
   a,b=self.pair(Path(tmp),"temporal")
   end=float(a["history"]["time_s"][-1])
   for run in (a,b):
    run["summary"]["local_energy_guard"]=dict(applicable=True,
     guard_active_since_time_ps=float(np.nextafter(end,0.))*1e12)
   result=audit.guarded_common_step_statistics(a,b,0.,end)
   self.assertFalse(result["available"])

class EnergyGuardAuditTests(unittest.TestCase):
 def fixture(self,anchor_ps=1.,defects=(2e-22,0.)):
  retained=2e-19;work=np.array([0.,1e-20,2e-20])
  energy=2e-16+work+np.r_[0.,np.cumsum(defects)]
  history=dict(time_s=np.array([0.,1e-12,2e-12]),energy_total_J=energy,external_energy_integral_J=work)
  dw=np.diff(work);defect=np.diff(energy)-dw;den=retained+np.cumsum(abs(dw));fraction=abs(defect)/den
  post=history["time_s"][:-1]>=anchor_ps*1e-12
  meta=dict(absolute_external_energy_integral_J=math.fsum(abs(float(v)) for v in dw),
   sum_abs_local_energy_defect_J=math.fsum(abs(float(v)) for v in defect),
   max_accepted_local_energy_fraction=float(max(fraction)),
   max_guarded_accepted_local_energy_fraction=float(max(fraction[post])) if np.any(post) else 0.,
   guard_active_since_time_ps=anchor_ps,local_energy_reference_energy_J=retained,
   local_energy_rejections=2,energy_guard_ledger_origin="synthetic_coherent_anchor")
  manifest=dict(parameter_identity={"local_energy_tolerance":1e-5},numerical={"local_energy_tolerance":1e-5},
   local_energy_control=copy.deepcopy(meta))
  return manifest,history,meta,retained
 def check(self,values):return audit.audit_local_energy_guard(*values)
 def test_legacy_large_defect_not_retroactively_failed(self):
  result=self.check(self.fixture())
  self.assertTrue(result["passed"])
  self.assertGreater(result["maximum_all_reconstructed_fraction"],1e-5)
  self.assertLess(result["maximum_guarded_reconstructed_fraction"],1e-5)
  self.assertEqual(result["legacy_accepted_steps"],1)
 def test_post_anchor_violation_fails_despite_consistent_metadata(self):
  result=self.check(self.fixture(defects=(0.,2e-22)))
  self.assertFalse(result["checks"]["post_anchor_local_defects_bounded"])
  self.assertFalse(result["passed"])
 def test_falsified_maximum_guarded_is_not_trusted(self):
  values=self.fixture();values[2]["max_guarded_accepted_local_energy_fraction"]+=1e-4
  values[0]["local_energy_control"]=copy.deepcopy(values[2])
  result=self.check(values)
  self.assertFalse(result["checks"]["maximum_guarded_ledger"])
 def test_falsified_absolute_ledger_is_not_trusted(self):
  for key in ("absolute_external_energy_integral_J","sum_abs_local_energy_defect_J"):
   values=self.fixture();values[2][key]+=1e-22;values[0]["local_energy_control"]=copy.deepcopy(values[2])
   result=self.check(values);self.assertFalse(result["passed"])
 def test_anchor_must_be_saved_accepted_time(self):
  result=self.check(self.fixture(anchor_ps=.5))
  self.assertFalse(result["checks"]["anchor_is_saved_accepted_time"])
 def test_guard_requires_strict_identity_and_complete_metadata(self):
  values=self.fixture();del values[0]["parameter_identity"]["local_energy_tolerance"]
  self.assertFalse(self.check(values)["passed"])
  values=self.fixture();del values[2]["guard_active_since_time_ps"]
  self.assertFalse(self.check(values)["passed"])
 def test_guard_reference_cannot_change_photon_normalization(self):
  values=self.fixture();values=list(values);values[-1]=1e-19
  result=self.check(values)
  self.assertFalse(result["checks"]["positive_photon_reference_matches_retained"])
 def test_manifest_cannot_disagree_with_checkpoint_ledger(self):
  values=self.fixture();values[0]["local_energy_control"]["absolute_external_energy_integral_J"]=0.
  self.assertFalse(self.check(values)["checks"]["checkpoint_manifest_guard_ledger_matches"])
 def test_legacy_without_declared_guard_is_not_certified(self):
  result=audit.audit_local_energy_guard({}, {}, {},2e-19)
  self.assertFalse(result["applicable"]);self.assertTrue(result["passed"])
if __name__=="__main__":unittest.main()
