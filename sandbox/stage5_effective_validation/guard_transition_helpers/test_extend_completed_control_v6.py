"""Strict20ps control extension: synthetic transactions, no physical evolution."""
import copy,importlib.util,json,types
from pathlib import Path
from unittest.mock import patch
import numpy as np
import pytest
HERE=Path(__file__).parent
def load(name):
 s=importlib.util.spec_from_file_location(name,HERE/(name+".py"));m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m
w=load("extend_completed_control_v6")
f=load("test_checkpoint_transfer_v3_v5_legacy")
v=load("test_checkpoint_transfer_v5_v6")
def save(run,arrays,manifest):
 np.savez_compressed(run/"checkpoint_latest.npz",**arrays);w.h.writejson(run/"manifest.json",manifest)
def source_initial(run,arrays):
 np.savez_compressed(run/"source_initial_state.npz",delta_bar=arrays["snapshots__delta_J"][0],
  p=arrays["p"],n=arrays["n"],temperature_K=np.ones(len(arrays["delta_bar"])),
  node_deposition_J=np.array([1.]),coordinates_m=np.zeros((len(arrays["delta_bar"]),2)))
def completed(run,arrays,parent,*,fresh):
 m=w.h.metadata(arrays);identity=copy.deepcopy(m["parameter_identity"])
 identity.update(duration_ps=20.,ceiling_ps=200.,local_energy_tolerance=1e-5)
 m.update(parameter_identity=identity,time_ps=20.,stable_accepted_steps=7,consecutive_rejections=2)
 m["controller"]["config"].update(no_trigger_time_ceiling_s=200e-12,threshold_V=1e-4)
 m["controller"]["state"].update(last_time_s=20e-12,rejected_candidate_count=0,failure_reason=None,
  **{k:None for k in w.EVENT_KEYS})
 m["source"].update(old_prepared_photon_state_reused=False)
 if fresh:
  m.update(absolute_external_energy_integral_J=1.,sum_abs_local_energy_defect_J=0.,local_energy_rejections=3,
   max_accepted_local_energy_fraction=0.,max_guarded_accepted_local_energy_fraction=0.,
   guard_active_since_time_ps=0.,local_energy_reference_energy_J=1.,
   energy_guard_ledger_origin="accepted_midpoint_receipts_since_time_zero")
 arrays["delta_bar"]=np.array([1.1+.1j])
 arrays["snapshots__delta_J"][-1]=arrays["delta_bar"]
 arrays["trajectory__time_s"][-1]=20e-12;arrays["snapshots__time_s"][-1]=20e-12
 arrays["metadata_json"]=np.array(json.dumps(m))
 threads=identity["threads"];start=0 if fresh else 12
 resources=dict(actual_processes=1,actual_threads=threads,requested_threads=threads,
  reserved_physical_core_count=threads,actual_affinity_cpus=list(range(start,start+threads)))
 parent.update(parameter_identity=identity,run_complete=True,status="DEVELOPMENT_DURATION_REACHED",
  source=m["source"],resources=resources,physicalparams={"energy_scale_J":1.,"circuit":{"R_load_ohm":50.}},
  observation={"development_duration_s":20e-12,"status":"DEVELOPMENT_DURATION_REACHED"})
 save(run,arrays,parent);source_initial(run,arrays)
 return m
def fresh_fixture(tmp_path):
 run,repo,a,m=f.fixture(tmp_path,counters=True);meta=completed(run,a,m,fresh=True)
 dc=w.h.decode((run/"dc/effective_dc.npz").read_bytes());dc.update(delta_bar=np.array([1.+.1j]),coordinates_m=np.zeros((1,2)))
 np.savez_compressed(run/"dc/effective_dc.npz",**dc)
 w.h.writejson(run/"dc/effective_dc_receipt.json",{"fields_sha256":w.h.filehash(run/"dc/effective_dc.npz")})
 started=dict(status="STARTED",output_dir=str(run),time_ps=0.,nodes=1,electron_states=2,phonon_modes=2,
  source=meta["source"],dt_ps=meta["parameter_identity"]["dt_ps"],**w.recorded_resources(m))
 log=tmp_path/"original.log";log.write_text(json.dumps(started)+"\n"+json.dumps({"status":"RUNNING","time_ps":20.})+"\n")
 return run,a,m,meta,log
def inherited_fixture(tmp_path):
 run,target,certificate,c=v.fixture(tmp_path)
 a=w.h.decode((run/"checkpoint_latest.npz").read_bytes());m=json.loads((run/"manifest.json").read_text())
 meta=completed(run,a,m,fresh=False)
 return run,a,m,meta,None
class FakeRunner:
 EffectiveState=staticmethod(lambda d,p,n,c:types.SimpleNamespace(delta_bar=d,p=p,n=n,circuit_SI=c))
 TriggerWindow=types.SimpleNamespace(from_checkpoint=lambda c:c)
 class Parser:
  def parse_args(self,argv):return types.SimpleNamespace(threads=4)
 parser=staticmethod(lambda:FakeRunner.Parser())
 @staticmethod
 def build_system(a,output,resume=False):
  assert resume
  system=types.SimpleNamespace(
   functional=types.SimpleNamespace(evaluate=lambda d,p:types.SimpleNamespace(energy_J=float(10+10*(d[0].real-1)))),
   kinetics=types.SimpleNamespace(phonon_energy_weights_J=np.zeros((1,2))),
   circuit=types.SimpleNamespace(energy=lambda c:0.),kinetic_residual_energy_scale_J=1.)
  g=types.SimpleNamespace(receipt={"mesh_sha256":"mesh"},coordinates_m=np.zeros((1,2)),
   triangles=np.empty((0,3),int),graph=types.SimpleNamespace(edges=np.empty((0,2),int)))
  return system,None,g,types.SimpleNamespace(manifest={"test":1}),{"energy_scale_J":1.,"circuit":{"R_load_ohm":50.}}
def sandbox(tmp_path):
 return patch.multiple(w.h,SCRATCH=tmp_path,load_runner=lambda p:FakeRunner)
@pytest.mark.parametrize("kind",["fresh","inherited"])
@pytest.mark.parametrize("duration",[40.,200.,None])
def test_only_horizon_changes_every_array_metadata_ledger_and_resource_preserved(tmp_path,kind,duration):
 with patch.object(v.h.lineage,"SCRATCH",tmp_path):
  run,a,m,meta,log=(fresh_fixture if kind=="fresh" else inherited_fixture)(tmp_path)
 before=(run/"checkpoint_latest.npz").read_bytes()
 with sandbox(tmp_path),patch.object(w.h.lineage,"SCRATCH",tmp_path):
  cap=w.capture(run,tmp_path/"anchor",fresh_start_log=log)
  receipt=w.extend(tmp_path/"anchor",tmp_path/"new",duration)
 assert (run/"checkpoint_latest.npz").read_bytes()==before
 after=w.h.decode((tmp_path/"new/checkpoint_latest.npz").read_bytes());newmeta=w.h.metadata(after)
 for k,value in a.items():
  if k!="metadata_json":assert w.h.arrayhash(value)==w.h.arrayhash(after[k])
 for k,value in meta.items():
  if k!="parameter_identity":assert newmeta[k]==value
 assert {k:v for k,v in meta["parameter_identity"].items() if k!="duration_ps"}=={k:v for k,v in newmeta["parameter_identity"].items() if k!="duration_ps"}
 assert newmeta["parameter_identity"]["duration_ps"]==duration
 assert receipt["pins_before"]==receipt["pins_after"] and len(receipt["pins_after"])==13
 assert receipt["guard_ledger_sha256_before"]==receipt["guard_ledger_sha256_after"]
 assert receipt["resources_before"]==receipt["resources_after"]
 assert receipt["RHS_sources_identical"] and receipt["energy_guard_change"] is None
 assert receipt["conserved_metadata_before"]==receipt["conserved_metadata_after"]
 assert receipt["lineage_kind"]==kind
 assert receipt["physical_identity_verified"]["no_collision_RHS_evaluated"]
 assert receipt["physical_identity_verified"]["photon_source_reapplied"] is False
 assert receipt["numerical_segments"][-1]["duration_ps"]==duration
 if kind=="fresh":
  assert receipt["fresh_initialization_proof"]["no_prior_transfer_invented"]
  assert not (tmp_path/"new/provenance/parent_transfer_receipt.json").exists()
 else:assert (tmp_path/"new/provenance/parent_transfer_receipt.json").exists()
 assert "--resume" in receipt["resume_argv"]
 assert ("--duration-ps" in receipt["resume_argv"])==(duration is not None)
 newmanifest=json.loads((tmp_path/"new/manifest.json").read_text())
 assert newmanifest["run_complete"] is False and newmanifest["status"]=="EXPLICIT_WINDOW_EXTENSION_PREPARED"
 assert newmanifest["observation"]["development_duration_s"]==(None if duration is None else duration*1e-12)
@pytest.mark.parametrize("duration",[19.,20.,200.01,201.,float("nan"),float("inf"),True,"unchanged"])
def test_invalid_or_implicit_horizon_is_rejected(tmp_path,duration):
 run,a,m,meta,log=fresh_fixture(tmp_path)
 with pytest.raises(ValueError):w.horizon(m,meta,duration)
@pytest.mark.parametrize("fault",["running","status","candidate","crossing","confirmed","target","stopped","terminal","failure","past_crossing","prior_not20","ceiling_over200","ceiling_mismatch","clock"])
def test_pending_terminal_or_ineligible_control_cannot_reopen(tmp_path,fault):
 run,a,m,meta,log=fresh_fixture(tmp_path)
 if fault=="running":m["run_complete"]=False
 elif fault=="status":m["status"]="TRIGGER_MARGIN_COMPLETE"
 elif fault in ("candidate","crossing","confirmed","target","stopped"):
  key=dict(candidate="candidate_crossing_s",crossing="crossing_time_s",confirmed="confirmation_time_s",target="stop_target_s",stopped="actual_stop_s")[fault]
  meta["controller"]["state"][key]=19e-12
 elif fault=="terminal":meta["controller"]["state"]["status"]="TIME_CENSORED_NO_TRIGGER"
 elif fault=="failure":meta["controller"]["state"]["failure_reason"]="failed"
 elif fault=="past_crossing":meta["controller"]["state"]["rejected_candidate_count"]=1
 elif fault=="prior_not20":meta["parameter_identity"]["duration_ps"]=30.
 elif fault=="ceiling_over200":meta["parameter_identity"]["ceiling_ps"]=201.
 elif fault=="ceiling_mismatch":meta["controller"]["config"]["no_trigger_time_ceiling_s"]=300e-12
 elif fault=="clock":meta["time_ps"]=19.99
 with pytest.raises(ValueError):w.horizon(m,meta,40.)
@pytest.mark.parametrize("fault",["missing_log","wrong_start","wrong_source","wrong_output","duplicate_start","guard_not_zero","guard_origin","initial_gap","deposition","coords","history_crossing"])
def test_fresh_proof_must_be_explicit_and_complete(tmp_path,fault):
 run,a,m,meta,log=fresh_fixture(tmp_path)
 if fault=="missing_log":log=None
 elif fault in ("wrong_start","wrong_source","wrong_output","duplicate_start"):
  e=json.loads(log.read_text().splitlines()[0])
  if fault=="wrong_start":e["time_ps"]=.02
  if fault=="wrong_source":e["source"]["retained_energy_J"]=2.
  if fault=="wrong_output":e["output_dir"]=str(tmp_path/"other")
  log.write_text(json.dumps(e)+"\n"+(json.dumps(e)+"\n" if fault=="duplicate_start" else ""))
 elif fault=="guard_not_zero":meta["guard_active_since_time_ps"]=.02
 elif fault=="guard_origin":meta["energy_guard_ledger_origin"]="inherited"
 elif fault in ("initial_gap","deposition","coords"):
  z=w.h.decode((run/"source_initial_state.npz").read_bytes())
  if fault=="initial_gap":z["delta_bar"]*=1.1
  if fault=="deposition":z["node_deposition_J"]*=2.
  if fault=="coords":z["coordinates_m"]+=1.
  np.savez_compressed(run/"source_initial_state.npz",**z)
 elif fault=="history_crossing":a["trajectory__Vout_V"][0]=.001
 with sandbox(tmp_path),pytest.raises(ValueError):
  w.no_history_crossing(a,meta);w.fresh_evidence(run,a,m,meta,log)
@pytest.mark.parametrize("fault",["runtime","source_file","manifest","guard_missing"])
def test_capture_or_anchored_evidence_tampering_is_rejected_before_output(tmp_path,fault):
 run,a,m,meta,log=fresh_fixture(tmp_path)
 with sandbox(tmp_path):
  if fault=="guard_missing":
   meta.pop("local_energy_rejections");a["metadata_json"]=np.array(json.dumps(meta));save(run,a,m)
   with pytest.raises(ValueError):w.capture(run,tmp_path/"anchor",fresh_start_log=log)
   assert not (tmp_path/"anchor").exists();return
  w.capture(run,tmp_path/"anchor",fresh_start_log=log)
  if fault=="runtime":
   name=next(iter(meta["source_pins_sha256"]));(tmp_path/"anchor/runtime_sources"/name).write_text("tamper")
  elif fault=="source_file":(tmp_path/"anchor/source_initial_state.npz").write_bytes(b"tamper")
  else:(tmp_path/"anchor/manifest_original.json").write_text("{}")
  with pytest.raises(ValueError):w.extend(tmp_path/"anchor",tmp_path/"new",40.)
 assert not (tmp_path/"new").exists()
@pytest.mark.parametrize("key,value",[("threads",2),("actual_threads",1),("actual_affinity_cpus",[0,2,3,4])])
def test_resource_changes_or_ambiguity_are_rejected(tmp_path,key,value):
 run,a,m,meta,log=fresh_fixture(tmp_path)
 if key=="threads":meta["parameter_identity"]["threads"]=value
 else:m["resources"][key]=value
 with pytest.raises(ValueError):w.resource_start(m,meta)
def test_existing_output_is_preserved(tmp_path):
 run,a,m,meta,log=fresh_fixture(tmp_path)
 with sandbox(tmp_path):
  w.capture(run,tmp_path/"anchor",fresh_start_log=log)
  (tmp_path/"new").mkdir();(tmp_path/"new/keep").write_text("original")
  with pytest.raises(FileExistsError):w.extend(tmp_path/"anchor",tmp_path/"new",40.)
 assert (tmp_path/"new/keep").read_text()=="original"

def test_ambiguous_affinity_alias_is_rejected(tmp_path):
 run,a,m,meta,log=fresh_fixture(tmp_path)
 m["resources"]["affinity_cpus"]=[4,5,6,7]
 with pytest.raises(ValueError,match="Ambiguous"):w.resource_start(m,meta)
