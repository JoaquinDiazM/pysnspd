"""Physical-width sensitivity cannot become numerical or event admission."""
import copy,json,math,tempfile,unittest,hashlib
from pathlib import Path
import numpy as np
from sandbox.stage5_effective_validation import physical_sensitivity as ps
from tests.test_effective_report_display import report

class EffectivePhysicalSensitivityTests(unittest.TestCase):
 def run_case(self,sigma=7.549834435):
  times=np.array([0.,15.3,20.])*1e-12
  parameters=dict(sigma_nm=sigma,wavelength_nm=775.,dt_ps=.02,step_tolerance=2e-7,
   local_energy_tolerance=1e-5,duration_ps=20.,threads=2,length_nm=320.,spacing_nm=5.)
  source=dict(sigma_nm=sigma,wavelength_nm=775.,photon_energy_J=2.563e-19,
   retained_energy_J=1.7096e-19,retained_fraction=.667,center_m=[160e-9,0.],
   time_origin="effective post-cascade",initialization="native energy FD/BE",maximum_initial_temperature_K=20.)
  manifest=dict(parameter_identity=parameters,source=source,physicalparams={"Tc_K":8.65},
   run_complete=True,source_pins_sha256={str(k):str(k)*64 for k in range(13)})
  summary=dict(passed=True,status="DEVELOPMENT_DURATION_REACHED",final_time_ps=20.,
   event=dict(confirmed=False,latency_s=None,confirmation_time_s=None),checks={"physical":True})
  return dict(manifest=manifest,summary=summary,traces=dict(time_s=times,Vout_relative_V=np.array([0.,1.,2.])*1e-6),
   snap=dict(time_s=times,Te_K=np.ones((3,2))*20.),spatial=dict(time_s=times,deficit_J_m2=np.array([0.,1.,2.])))
 def inputs(self):
  ref=self.run_case()
  cases={"narrow":self.run_case(ref["manifest"]["parameter_identity"]["sigma_nm"]/math.sqrt(2)),
   "wide":self.run_case(ref["manifest"]["parameter_identity"]["sigma_nm"]*math.sqrt(2))}
  return dict(reference_wavelength_nm=775,common_interval_ps=[0.,20.],map_time_ps=15.3),ref,cases
 def test_width_sensitivity_is_separate_from_numerical_and_experimental_admission(self):
  cfg,ref,cases=self.inputs();before=copy.deepcopy(cases["narrow"]["summary"])
  result=ps.analyze_sensitivity(cfg,ref,cases)
  self.assertFalse(result["numerical_admission_axis"]);self.assertFalse(result["experimental_calibration"])
  self.assertFalse(result["stage5_admission_inferred"]);self.assertEqual(cases["narrow"]["summary"],before)
  self.assertEqual({row["status"] for row in result["cases"]},{"DEVELOPMENT_DURATION_REACHED"})
  self.assertTrue(all(row["latency_in_this_window_ps"] is None for row in result["cases"]))
 def test_initial_temperature_comes_from_saved_postcascade_snapshot(self):
  cfg,ref,cases=self.inputs();cases["narrow"]["manifest"]["source"]["maximum_initial_temperature_K"]=24.7
  cases["narrow"]["snap"]["Te_K"][0,:]=24.7
  rows={row["label"]:row for row in ps.analyze_sensitivity(cfg,ref,cases)["cases"]}
  self.assertEqual(rows["narrow"]["maximum_initial_temperature_K"],24.7)
  self.assertEqual(rows["nominal"]["maximum_initial_temperature_K"],20.)
 def test_fabricated_initial_temperature_is_rejected(self):
  cfg,ref,cases=self.inputs();cases["wide"]["manifest"]["source"]["maximum_initial_temperature_K"]=99.
  with self.assertRaisesRegex(ValueError,"actual saved"):ps.analyze_sensitivity(cfg,ref,cases)
 def test_own_saved_event_is_used_only_inside_common_window(self):
  cfg,ref,cases=self.inputs()
  cases["narrow"]["summary"]["event"]=dict(confirmed=True,latency_s=14e-12,confirmation_time_s=14.5e-12)
  cases["wide"]["summary"]["event"]=dict(confirmed=True,latency_s=21e-12,confirmation_time_s=21.5e-12)
  rows={r["label"]:r for r in ps.analyze_sensitivity(cfg,ref,cases)["cases"]}
  self.assertAlmostEqual(rows["narrow"]["latency_in_this_window_ps"],14.)
  self.assertIsNone(rows["wide"]["latency_in_this_window_ps"])
 def test_every_physical_or_solver_parameter_except_sigma_is_protected(self):
  cfg,ref,cases=self.inputs()
  for key,value in (("dt_ps",.01),("step_tolerance",2e-9),("local_energy_tolerance",1.25e-6),
   ("wavelength_nm",1550.),("spacing_nm",3.5),("length_nm",480.)):
   candidate=copy.deepcopy(cases["narrow"]);candidate["manifest"]["parameter_identity"][key]=value
   with self.assertRaises(ValueError):ps.validate_case_parameters(ref,candidate,"narrow")
 def test_changed_source_or_material_is_rejected(self):
  cfg,ref,cases=self.inputs()
  for key,value in (("retained_energy_J",1e-19),("retained_fraction",.5),("center_m",[0.,0.]),("time_origin","absorption")):
   candidate=copy.deepcopy(cases["narrow"]);candidate["manifest"]["source"][key]=value
   with self.assertRaises(ValueError):ps.validate_case_parameters(ref,candidate,"narrow")
  candidate=copy.deepcopy(cases["narrow"]);candidate["manifest"]["physicalparams"]["Tc_K"]=9.
  with self.assertRaises(ValueError):ps.validate_case_parameters(ref,candidate,"narrow")
 def test_source_pins_cannot_be_ignored(self):
  cfg,ref,cases=self.inputs()
  for action in ("missing","changed"):
   candidate=copy.deepcopy(cases["narrow"])
   if action=="missing":candidate["manifest"]["source_pins_sha256"].pop("1")
   else:candidate["manifest"]["source_pins_sha256"]["1"]="f"*64
   with self.assertRaises(ValueError):ps.validate_case_parameters(ref,candidate,"narrow")
 def test_exact_common_map_time_no_nearest_or_interpolated_state(self):
  cfg,ref,cases=self.inputs();cfg["map_time_ps"]=15.35
  with self.assertRaisesRegex(ValueError,"actual saved snapshot"):ps.analyze_sensitivity(cfg,ref,cases)
  self.assertEqual(ps.exact_saved_index(np.array([0.,15.3,20.])*1e-12,15.3),1)
 def test_common_window_coverage_cannot_be_faked(self):
  cfg,ref,cases=self.inputs();cases["wide"]["traces"]["time_s"][-1]=19.9e-12
  with self.assertRaisesRegex(ValueError,"cover declared"):ps.analyze_sensitivity(cfg,ref,cases)
 def test_failed_or_incomplete_physical_case_is_rejected(self):
  for key in ("passed","run_complete"):
   cfg,ref,cases=self.inputs()
   target=cases["narrow"]["summary"] if key=="passed" else cases["narrow"]["manifest"]
   target[key]=False
   with self.assertRaises(ValueError):ps.analyze_sensitivity(cfg,ref,cases)
 def test_width_ratios_and_case_set_are_explicit(self):
  cfg,ref,cases=self.inputs();cases["narrow"]["manifest"]["parameter_identity"]["sigma_nm"]*=.9
  with self.assertRaises(ValueError):ps.analyze_sensitivity(cfg,ref,cases)
  cfg,ref,cases=self.inputs();cases.pop("wide")
  with self.assertRaises(ValueError):ps.analyze_sensitivity(cfg,ref,cases)
 def test_planned_page_count_includes_optional_sensitivity_and_vortices(self):
  self.assertEqual(report.expected_report_pages({},False),11)
  self.assertEqual(report.expected_report_pages({},True),12)
  self.assertEqual(report.expected_report_pages({"include_vortex_diagnostics":True},True),14)
 def test_sensitivity_cannot_bypass_stage5_main_gate(self):
  root=Path("/home/jdiaz/scratch/stage5_report_physical_sensitivity_source_20261011")
  with tempfile.TemporaryDirectory(dir=root) as directory:
   p=Path(directory)/"unadmitted.json";p.write_text(json.dumps({"stage5":False,"physical_sensitivity_runs":self.inputs()[0]}))
   with self.assertRaisesRegex(ValueError,"stage5: true"):report.main(["--admission",str(p),"--build-root",str(Path(directory)/"forbidden")])
   self.assertFalse((Path(directory)/"forbidden").exists())

 def loader_fixture(self,directory):
  cfg,ref,cases=self.inputs();cfg["cases"]={}
  for label,run in cases.items():
   root=Path(directory)/label;root.mkdir()
   pins={}
   for filename in ps.FILES:
    path=root/filename;path.write_bytes((label+filename).encode())
    pins[filename]=hashlib.sha256(path.read_bytes()).hexdigest()
   audit=root/"audit.json";summary=copy.deepcopy(run["summary"]);summary["path"]=str(root)
   document=dict(schema="pysnspd.effective_saved_data_audit.v2",all_requested_checks_passed=True,errors={},
    auditor_source_sha256=hashlib.sha256(Path(ps.audit.__file__).read_bytes()).hexdigest(),runs={label:summary})
   audit.write_text(json.dumps(document))
   cfg["cases"][label]=dict(path=str(root),sha256=pins,audit=dict(path=str(audit),sha256=hashlib.sha256(audit.read_bytes()).hexdigest()))
  return cfg,ref,cases
 def test_loader_requires_and_verifies_all_data_and_audit_sha(self):
  root=Path("/home/jdiaz/scratch/stage5_report_physical_sensitivity_source_20261011")
  with tempfile.TemporaryDirectory(dir=root) as directory:
   cfg,ref,cases=self.loader_fixture(directory);hashes={}
   result,runs,records=ps.load_sensitivity(cfg,ref,input_hashes=hashes,verify_entry=report.pinned_entry,load_run=lambda label,path:cases[label])
   self.assertEqual(set(runs),{"nominal","narrow","wide"});self.assertEqual(len(records),2);self.assertEqual(len(hashes),12)
   cfg["cases"]["narrow"]["sha256"]["checkpoint_latest.npz"]="f"*64
   with self.assertRaisesRegex(ValueError,"checksum changed"):
    ps.load_sensitivity(cfg,ref,input_hashes={},verify_entry=report.pinned_entry,load_run=lambda label,path:cases[label])
 def test_missing_file_pin_or_audit_pin_is_rejected(self):
  root=Path("/home/jdiaz/scratch/stage5_report_physical_sensitivity_source_20261011")
  for field in ("sha256","audit"):
   with tempfile.TemporaryDirectory(dir=root) as directory:
    cfg,ref,cases=self.loader_fixture(directory)
    if field=="sha256":cfg["cases"]["narrow"][field].pop("snapshots.npz")
    else:cfg["cases"]["narrow"][field].pop("sha256")
    with self.assertRaises(ValueError):
     ps.load_sensitivity(cfg,ref,input_hashes={},verify_entry=report.pinned_entry,load_run=lambda label,path:cases[label])
 def test_input_arrays_changing_during_analysis_are_rejected(self):
  root=Path("/home/jdiaz/scratch/stage5_report_physical_sensitivity_source_20261011")
  with tempfile.TemporaryDirectory(dir=root) as directory:
   cfg,ref,cases=self.loader_fixture(directory)
   def changing(label,path):
    (path/"trajectory.npz").write_bytes(b"changed after pin verification")
    return cases[label]
   with self.assertRaisesRegex(ValueError,"checksum changed"):
    ps.load_sensitivity(cfg,ref,input_hashes={},verify_entry=report.pinned_entry,load_run=changing)

if __name__=="__main__":unittest.main()
