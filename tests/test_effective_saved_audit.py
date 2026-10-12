"""Saved-data energy, comparison and honest-event controls for effective audits."""
import copy,importlib.util,json,tempfile,unittest
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location("effective_saved_audit",ROOT/"sandbox/stage5_effective_validation/audit_runs.py")
audit=importlib.util.module_from_spec(spec);spec.loader.exec_module(audit)

class EffectiveSavedAuditTests(unittest.TestCase):
 def fixture(self,path):
  xy=np.array([[x,y] for y in (-40,0,40) for x in (0,80,160)],float)*1e-9
  edges=np.array([[0,1],[1,2],[3,4],[4,5],[6,7],[7,8],[0,3],[3,6],[1,4],[4,7],[2,5],[5,8]])
  tri=np.array([[0,1,4],[0,4,3],[1,2,5],[1,5,4],[3,4,7],[3,7,6],[4,5,8],[4,8,7]])
  gap=np.full((3,9),audit.MEV_J,complex);gap[1:, [1,4,7]]*=.1
  snap=dict(time_s=np.array([0.,1e-12,2e-12]),coordinates_m=xy,edges=edges,triangles=tri,
   delta_J=gap,Te_K=np.full((3,9),.9),Tph_K=np.full((3,9),.9))
  areas,_=audit.node_areas({},snap)
  snap["u_excess_J_m3"]=np.zeros((3,9))
  snap["u_excess_J_m3"][0,4]=2e-19/(areas[4]*7e-9)
  E0=2e-16;external=np.array([0.,1e-20,0.]);energy=E0+external+np.array([0.,2e-22,0.])
  hist=dict(time_s=snap["time_s"],Vout_V=np.array([0.,1e-4,2e-4]),Vdev_V=np.zeros(3),
   Ib_A=np.full(3,21e-6),Is_A=np.array([21.,19.,17.])*1e-6,
   energy_total_J=energy,energy_defect_J=energy-E0-external,
   external_energy_integral_J=external,instant_power_residual_W=np.zeros(3),
   escape_power_W=np.zeros(3),reservoir_power_out_W=np.zeros(3))
  manifest=dict(run_complete=True,status="DEVELOPMENT_DURATION_REACHED",
   source={"retained_energy_J":2e-19},physicalparams={"film_thickness_m":7e-9,"circuit":{"R_bias_ohm":1e4,"R_load_ohm":50.,"V_bias_V":.21}},energy_reference={"initial_total_energy_J":E0},
   observation={"mode":"photon","trigger_confirmed":False,"baseline_Vout_V":0.,"threshold_voltage_V":1e-4},
   source_pins_sha256={"solver":"synthetic_test_hash"},
   parameter_identity={"dt_ps":.01,"length_nm":160.,"spacing_nm":5.,"sigma_nm":7.55,"xi_cutoff":32.,"xi_order":6,"phonon_nodes":96,"transport_order":48})
  p=np.full((9,2),.1);n=np.full((9,3),.01);target=np.zeros(9);target[4]=2e-19
  np.savez_compressed(path/"trajectory.npz",**hist);np.savez_compressed(path/"snapshots.npz",**snap)
  np.savez_compressed(path/"source_initial_state.npz",p=p,n=n,node_deposition_J=target,coordinates_m=xy)
  np.savez_compressed(path/"checkpoint_latest.npz",delta_bar=np.ones(9,complex),p=p,n=n,
   circuit_SI=np.array([21e-6,17e-6,0.]),metadata_json=np.array(json.dumps({"time_ps":2.})))
  manifest["output_sha256"]={name:audit.sha256(path/name) for name in
   ("trajectory.npz","snapshots.npz","source_initial_state.npz","checkpoint_latest.npz")}
  (path/"manifest.json").write_text(json.dumps(manifest))
  return manifest,hist,snap
 def workspace(self):
  path=Path("/home/jdiaz/scratch/stage5_external_audit_development_20261011")
  path.mkdir(exist_ok=True);return tempfile.TemporaryDirectory(dir=path)
 def test_energy_denominator_keeps_absolute_external_work_despite_cancellation(self):
  with self.workspace() as tmp:
   path=Path(tmp);self.fixture(path);run=audit.audit_run("test",path)
   self.assertTrue(run["summary"]["passed"])
   self.assertAlmostEqual(run["summary"]["energy"]["absolute_external_work_J"]/2e-20,1.,delta=2e-15)
   self.assertAlmostEqual(run["traces"]["energy_relative_defect"][1],
    abs(run["traces"]["energy_recomputed_defect_J"][1])/(2e-19+1e-20),delta=1e-15)
   self.assertFalse(run["summary"]["stage5_admission_inferred"])
 def test_hotbelt_is_saved_connected_transverse_diagnostic_not_electrical_trigger(self):
  with self.workspace() as tmp:
   path=Path(tmp);self.fixture(path);run=audit.audit_run("test",path)
   first=run["summary"]["hotbelt"]["first_saved_width_spanning"]
   self.assertEqual(first["bracket_s"],[0.,1e-12])
   self.assertTrue(run["hotbelt"][1]["width_spanning"])
   self.assertFalse(run["summary"]["event"]["confirmed"])
   self.assertIsNone(run["summary"]["event"]["latency_s"])
 def test_tampered_endpoint_support_and_checksum_fail(self):
  with self.workspace() as tmp:
   path=Path(tmp);self.fixture(path);data=audit.load_npz(path/"checkpoint_latest.npz")
   data["p"][0,0]=-1e-20;np.savez_compressed(path/"checkpoint_latest.npz",**data)
   run=audit.audit_run("test",path)
   self.assertFalse(run["summary"]["checks"]["checkpoint_support"])
   self.assertFalse(run["summary"]["checks"]["saved_output_checksums"])
   self.assertFalse(run["summary"]["passed"])
 def test_confirmed_event_requires_saved_fields_and_consistent_crossing(self):
  with self.workspace() as tmp:
   path=Path(tmp);manifest,hist,snap=self.fixture(path)
   manifest["observation"].update(trigger_confirmed=True,crossing_time_s=1e-12,
    confirmation_time_s=2e-12,latency_s=1e-12)
   (path/"manifest.json").write_text(json.dumps(manifest));run=audit.audit_run("test",path)
   self.assertTrue(run["summary"]["checks"]["saved_electrical_event_consistent"])
   manifest["observation"]["crossing_time_s"]=.5e-12
   (path/"manifest.json").write_text(json.dumps(manifest));run=audit.audit_run("test",path)
   self.assertFalse(run["summary"]["checks"]["saved_electrical_event_consistent"])
   with self.assertRaises(ValueError):
    audit.event_record({"observation":{"trigger_confirmed":True}},2e-12)
 def test_common_time_comparison_and_temporal_latency_absolute_allowance(self):
  with self.workspace() as tmp:
   path=Path(tmp);self.fixture(path);ref=audit.audit_run("base",path);candidate=copy.deepcopy(ref)
   candidate["summary"]["name"]="half";candidate["manifest"]["parameter_identity"]["dt_ps"]=.005
   result,curves,gap=audit.compare_runs("temporal",ref,candidate)
   self.assertTrue(result["passed"]);self.assertIsNone(result["latency"]["passed"])
   np.testing.assert_array_equal(curves["difference_Vout_mV"],0.)
   for run,latency in ((ref,1e-12),(candidate,1.09e-12)):
    run["summary"]["event"].update(confirmed=True,latency_s=latency,crossing_time_s=latency,confirmation_time_s=2e-12)
   result,_,_=audit.compare_runs("temporal",ref,candidate)
   self.assertTrue(result["latency"]["passed"])
   candidate["summary"]["event"]["latency_s"]=1.11e-12
   result,_,_=audit.compare_runs("temporal",ref,candidate)
   self.assertFalse(result["latency"]["passed"])
 def test_refinement_cannot_silently_change_photon_or_physical_parameters(self):
  with self.workspace() as tmp:
   path=Path(tmp);self.fixture(path);ref=audit.audit_run("base",path);candidate=copy.deepcopy(ref)
   candidate["summary"]["name"]="half";candidate["manifest"]["parameter_identity"].update(dt_ps=.005,sigma_nm=10.)
   result,_,_=audit.compare_runs("temporal",ref,candidate)
   self.assertFalse(result["checks"]["same_nonrefined_parameters"])
   self.assertIn("sigma_nm",result["forbidden_parameter_changes"])
 def test_spatial_source_energy_is_independent_of_unchanged_target_vector(self):
  with self.workspace() as tmp:
   path=Path(tmp);manifest,hist,snap=self.fixture(path)
   snap["u_excess_J_m3"][0]*=.5
   np.savez_compressed(path/"snapshots.npz",**snap)
   manifest["output_sha256"]["snapshots.npz"]=audit.sha256(path/"snapshots.npz")
   (path/"manifest.json").write_text(json.dumps(manifest))
   run=audit.audit_run("test",path)
   self.assertTrue(run["summary"]["checks"]["source_initial_support_and_deposition"])
   self.assertTrue(run["summary"]["checks"]["saved_output_checksums"])
   self.assertFalse(run["summary"]["checks"]["independent_initial_source_energy"])
   self.assertAlmostEqual(run["summary"]["initial_source_density_integral"]["relative_defect"],.5)
 def test_refinement_controls_require_a_real_ordered_change(self):
  with self.workspace() as tmp:
   path=Path(tmp);self.fixture(path);ref=audit.audit_run("base",path)
   for kind,key,value in (("mesh","spacing_nm",2.5),("domain","length_nm",320.),
     ("spectral","xi_order",8)):
    candidate=copy.deepcopy(ref);candidate["summary"]["name"]="refined"
    result,_,_=audit.compare_runs(kind,ref,candidate)
    self.assertFalse(result["checks"]["actual_declared_refinement"])
    candidate["manifest"]["parameter_identity"][key]=value
    result,_,_=audit.compare_runs(kind,ref,candidate)
    self.assertTrue(result["checks"]["actual_declared_refinement"])
   candidate["manifest"]["parameter_identity"]["phonon_nodes"]=48
   result,_,_=audit.compare_runs("spectral",ref,candidate)
   self.assertFalse(result["checks"]["actual_declared_refinement"])
 def test_arbitrary_single_output_hash_does_not_authenticate_used_data(self):
  with self.workspace() as tmp:
   path=Path(tmp);manifest,_,_=self.fixture(path)
   manifest["output_sha256"]={"source_initial_state.npz":audit.sha256(path/"source_initial_state.npz")}
   (path/"manifest.json").write_text(json.dumps(manifest))
   run=audit.audit_run("test",path)
   self.assertFalse(run["summary"]["checks"]["saved_output_checksums"])
 def test_instantaneous_power_defect_fails_despite_exact_integrated_ledger(self):
  with self.workspace() as tmp:
   path=Path(tmp);manifest,hist,snap=self.fixture(path)
   hist["instant_power_residual_W"][1]=1e-12
   np.savez_compressed(path/"trajectory.npz",**hist)
   manifest["output_sha256"]["trajectory.npz"]=audit.sha256(path/"trajectory.npz")
   (path/"manifest.json").write_text(json.dumps(manifest))
   run=audit.audit_run("test",path)
   self.assertTrue(run["summary"]["checks"]["energy_relative_defect"])
   self.assertTrue(run["summary"]["checks"]["saved_output_checksums"])
   self.assertFalse(run["summary"]["checks"]["instantaneous_energy_balance"])
   self.assertGreater(run["summary"]["instantaneous_power"]["maximum_residual_over_tolerance"],1.)
 def test_instantaneous_scale_uses_net_CM_power_and_retains_device_escape_cuts(self):
  with self.workspace() as tmp:
   path=Path(tmp);manifest,hist,snap=self.fixture(path)
   hist["Vdev_V"][:]=1e-4;hist["escape_power_W"][:]=2e-9
   hist["reservoir_power_out_W"][:]=-3e-9
   np.savez_compressed(path/"trajectory.npz",**hist)
   manifest["output_sha256"]["trajectory.npz"]=audit.sha256(path/"trajectory.npz")
   (path/"manifest.json").write_text(json.dumps(manifest))
   run=audit.audit_run("test",path)
   net=.21*hist["Ib_A"]-1e4*hist["Ib_A"]**2-50*(hist["Ib_A"]-hist["Is_A"])**2
   expected=abs(hist["Is_A"]*hist["Vdev_V"])+2e-9+3e-9+abs(net)
   np.testing.assert_allclose(run["traces"]["instant_power_scale_W"],expected,rtol=2e-15,atol=1e-30)
   np.testing.assert_allclose(run["traces"]["instant_power_tolerance_W"],1e-20+1e-10*expected,rtol=2e-15,atol=0.)
 def test_no_trigger_is_not_given_a_latency_even_if_stray_time_exists(self):
  event=audit.event_record({"observation":{"trigger_confirmed":False,"trigger_time_s":1e-12,
   "status":"TIME_CENSORED_NO_TRIGGER","censored":True}},200e-12)
  self.assertIsNone(event["latency_s"]);self.assertFalse(event["confirmed"]);self.assertTrue(event["censored"])

if __name__=="__main__":unittest.main()
