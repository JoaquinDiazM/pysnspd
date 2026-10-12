"""Synthetic interruption/transfer checks; no physical evolution."""
import copy,importlib.util,json,types
from pathlib import Path
from unittest.mock import patch
import numpy as np
import pytest

HERE=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location("transfer_helper",HERE/"checkpoint_transfer_v3_v5.py")
h=importlib.util.module_from_spec(spec);spec.loader.exec_module(h)

def fixture(tmp_path,counters=False):
 run=tmp_path/"parent";run.mkdir();rt=run/"runtime_sources";rt.mkdir()
 repo=tmp_path/"new_repo";repo.mkdir()
 names=sorted(h.FRONTEND)+["pysnspd/experimental/core_%02d.py"%i for i in range(11)]
 pins={}
 for rel in names:
  old=rt/rel;old.parent.mkdir(parents=True,exist_ok=True);old.write_text("old "+rel)
  new=repo/rel;new.parent.mkdir(parents=True,exist_ok=True)
  new.write_text(("new " if rel in h.FRONTEND else "old ")+rel)
  pins[rel]=h.filehash(old)
 arc={"runtime_pins_sha256":pins,"package_sources_and_small_geometry":{r:{"sha256":v} for r,v in pins.items()}}
 h.writejson(rt/"archive_manifest.json",arc)
 (run/"dc").mkdir()
 np.savez(run/"dc/effective_dc.npz",delta_bar=np.array([1.]),p=np.array([[.1,.2]]),phi_V=np.array([0.]))
 h.writejson(run/"dc/effective_dc_receipt.json",{"fields_sha256":h.filehash(run/"dc/effective_dc.npz")})
 np.savez(run/"material.npz",x=np.array([1.]))
 h.writejson(run/"material.json",{"test":1})
 np.savez(run/"source_initial_state.npz",p=np.array([[.1,.2]]))
 identity={"dt_ps":.02,"step_tolerance":2e-9,"threads":4}
 t=.02e-12
 controller={"state":{"status":"ACTIVE","last_time_s":t,"last_excursion_V":0.,"sample_count":2},
  "config":{"baseline_V":0.,"polarity":"positive"}}
 meta={"schema":"fixture","time_ps":.02,"initial_total_energy_J":10.,"external_energy_integral_J":1.,
  "limited_event_energy_J":.5,"totalabs_event_energy_J":2.,"accepted_steps":1,"total_rejections":3,
  "next_step_ps":.02,"source":{"source_must_not_be_applied_again":True,"retained_energy_J":1.},
  "controller":controller,"source_pins_sha256":pins,"parameter_identity":identity}
 if counters:meta.update(stable_accepted_steps=7,consecutive_rejections=2)
 arrays={"delta_bar":np.array([1.+.1j]),"p":np.array([[.1,.2]]),"n":np.array([[.3,.4]]),
  "circuit_SI":np.array([1.,1.,0.]),"metadata_json":np.array(json.dumps(meta))}
 histories={"time_s":[0.,t],"Ib_A":[1.,1.],"Is_A":[1.,1.],"vc_V":[0.,0.],"Vout_V":[0.,0.],
  "energy_total_J":[10.,11.],"energy_defect_J":[0.,0.],
  "external_energy_integral_J":[0.,1.],"limited_event_energy_J":[0.,.5],"totalabs_event_energy_J":[0.,2.]}
 arrays.update({"trajectory__"+k:np.array(v) for k,v in histories.items()})
 arrays.update({"snapshots__time_s":np.array([0.,t]),"snapshots__delta_J":np.array([[1.+.1j],[1.+.1j]])})
 np.savez_compressed(run/"checkpoint_latest.npz",**arrays)
 manifest={"schema":"fixture","parameter_identity":identity,"source_pins_sha256":pins,
  "physicalparams":{"circuit":{"R_load_ohm":50.}},"geometry":{"mesh_sha256":"mesh"},
  "material":{"test":1},"numerical":{}}
 h.writejson(run/"manifest.json",manifest)
 return run,repo,arrays,manifest

class FakeParser:
 def parse_args(self,args):return types.SimpleNamespace(threads=4)
class FakeRunner:
 EffectiveState=staticmethod(lambda *a:a)
 TriggerWindow=types.SimpleNamespace(from_checkpoint=lambda c:c)
 parser=staticmethod(lambda:FakeParser())
 @staticmethod
 def build_system(a,output,resume=False):
  assert resume
  system=types.SimpleNamespace(evaluate=lambda state:types.SimpleNamespace(energy_J=11.))
  g=types.SimpleNamespace(receipt={"mesh_sha256":"mesh"},coordinates_m=np.zeros((1,2)),
   triangles=np.empty((0,3),int),graph=types.SimpleNamespace(edges=np.empty((0,2),int)))
  return system,None,g,types.SimpleNamespace(manifest={"test":1}),{"circuit":{"R_load_ohm":50.}}

def test_capture_once_preserves_packet(tmp_path):
 run,repo,arr,m=fixture(tmp_path);before=(run/"checkpoint_latest.npz").read_bytes()
 with patch.object(h,"SCRATCH",tmp_path):
  r=h.capture(run,tmp_path/"anchor")
 assert (tmp_path/"anchor/checkpoint_original.npz").read_bytes()==before
 assert r["conserved_metadata"]["time_ps"]==.02
 assert (run/"checkpoint_latest.npz").read_bytes()==before
 assert h.filehash(tmp_path/"anchor/runtime_sources/pipelines/05_effective_photon_run.py")==m["source_pins_sha256"]["pipelines/05_effective_photon_run.py"]

@pytest.mark.parametrize("fault",["clock","ledger","counter","controller","pauli","bose","reversal","nan"])
def test_rejects_incoherent_or_nonphysical_endpoint(tmp_path,fault):
 run,repo,a,m=fixture(tmp_path);meta=h.metadata(a)
 if fault=="clock":meta["time_ps"]+=.001
 if fault=="ledger":meta["external_energy_integral_J"]+=.1
 if fault=="counter":meta["accepted_steps"]+=1
 if fault=="controller":meta["controller"]["state"]["last_time_s"]+=1e-15
 if fault=="pauli":a["p"][0,0]=1.01
 if fault=="bose":a["n"][0,0]=-.01
 if fault=="reversal":a["circuit_SI"][1]=-.1
 if fault=="nan":meta["time_ps"]=float("nan")
 a["metadata_json"]=np.array(json.dumps(meta))
 with pytest.raises(ValueError):h.audit_arrays(a,m)

def test_runtime_transitives_are_required(tmp_path):
 run,repo,a,m=fixture(tmp_path)
 (run/"runtime_sources/pysnspd/experimental/core_00.py").write_text("corrupt")
 with pytest.raises(ValueError,match="runtime pin mismatch"):h.audit_runtime(run/"runtime_sources",m["source_pins_sha256"])

@pytest.mark.parametrize("counters",[False,True])
def test_explicit_transfer_arrays_source_ledger_clock_and_counters(tmp_path,counters):
 run,repo,a,m=fixture(tmp_path,counters=counters)
 with patch.object(h,"SCRATCH",tmp_path),patch.object(h,"load_runner",lambda p:FakeRunner):
  h.capture(run,tmp_path/"anchor")
  r=h.transfer(tmp_path/"anchor",tmp_path/"new",repo,.01,2e-7,2,12)
 assert r["state_arrays_sha256_before"]==r["state_arrays_sha256_after"]
 assert r["conserved_metadata_before"]==r["conserved_metadata_after"]
 assert r["RHS_sources_identical"] and set(r["frontend_changes"])==h.FRONTEND
 assert r["next_step_ps_after"]==.01
 assert r["numerical_segments"][1]=={"start_time_ps":.02,"end_time_ps":None,"origin":"v4","dt_ps":.01,"step_tolerance":2e-7,"threads":2}
 moved=h.decode((tmp_path/"new/checkpoint_latest.npz").read_bytes())
 mm=h.metadata(moved)
 assert mm["stable_accepted_steps"]==(7 if counters else 0)
 assert mm["consecutive_rejections"]==(2 if counters else 0)
 assert r["missing_old_adaptive_counters"]==([] if counters else ["stable_accepted_steps","consecutive_rejections"])
 assert (tmp_path/"new/provenance/parent_checkpoint.npz").read_bytes()==(run/"checkpoint_latest.npz").read_bytes()
 assert (tmp_path/"new/provenance/resume_command.txt").is_file()
 assert h.filehash(tmp_path/"new/source_initial_state.npz")==h.filehash(run/"source_initial_state.npz")

def test_transfer_refuses_rhs_change(tmp_path):
 run,repo,a,m=fixture(tmp_path)
 with patch.object(h,"SCRATCH",tmp_path):
  h.capture(run,tmp_path/"anchor")
  (repo/"pysnspd/experimental/core_00.py").write_text("new physics")
  with pytest.raises(ValueError,match="Unapproved RHS"):h.transfer(tmp_path/"anchor",tmp_path/"new",repo,.02,2e-7,2,12)
 assert not (tmp_path/"new").exists()

def test_refuses_existing_output_or_anchored_mutation(tmp_path):
 run,repo,a,m=fixture(tmp_path)
 with patch.object(h,"SCRATCH",tmp_path),patch.object(h,"load_runner",lambda p:FakeRunner):
  h.capture(run,tmp_path/"anchor")
  with pytest.raises(FileExistsError):h.capture(run,tmp_path/"anchor")
  h.transfer(tmp_path/"anchor",tmp_path/"new",repo,.02,2e-7,2,12)
  with pytest.raises(FileExistsError):h.transfer(tmp_path/"anchor",tmp_path/"new",repo,.02,2e-7,2,12)
  (tmp_path/"anchor/material.json").write_text("{}")
  with pytest.raises(ValueError,match="Anchored physical"):h.transfer(tmp_path/"anchor",tmp_path/"another",repo,.02,2e-7,2,12)

def completed_control_fixture(tmp_path):
 run,repo,a,m=fixture(tmp_path,counters=True)
 meta=h.metadata(a);identity=copy.deepcopy(meta["parameter_identity"])
 identity.update(duration_ps=20.,ceiling_ps=200.)
 meta["parameter_identity"]=identity;meta["time_ps"]=20.
 meta["controller"]["state"]["last_time_s"]=20e-12
 a["metadata_json"]=np.array(json.dumps(meta))
 a["trajectory__time_s"][-1]=20e-12;a["snapshots__time_s"][-1]=20e-12
 np.savez_compressed(run/"checkpoint_latest.npz",**a)
 m.update(parameter_identity=identity,run_complete=True,status="DEVELOPMENT_DURATION_REACHED",
  observation={"development_duration_s":20e-12,"status":"DEVELOPMENT_DURATION_REACHED"})
 h.writejson(run/"manifest.json",m)
 return run,repo,a,m

@pytest.mark.parametrize("new_duration",[40.,None])
def test_completed_active_control_extends_only_horizon(tmp_path,new_duration):
 run,repo,a,m=completed_control_fixture(tmp_path)
 with patch.object(h,"SCRATCH",tmp_path),patch.object(h,"load_runner",lambda p:FakeRunner):
  h.capture(run,tmp_path/"anchor")
  r=h.transfer(tmp_path/"anchor",tmp_path/"new",repo,.02,2e-9,4,0,duration_ps_override=new_duration)
 assert r["state_arrays_sha256_before"]==r["state_arrays_sha256_after"]
 assert r["conserved_metadata_before"]==r["conserved_metadata_after"]
 change=r["horizon_change"]
 assert change["runtime_key"]=="duration_ps" and change["before_ps"]==20.
 assert change["after_ps"]==new_duration and change["effective_time_ps"]==20.
 assert change["unchanged_no_trigger_ceiling_ps"]==200.
 assert change["controller_status_preserved"]=="ACTIVE" and change["no_event_reopened"]
 newmeta=h.metadata(h.decode((tmp_path/"new/checkpoint_latest.npz").read_bytes()))
 assert newmeta["time_ps"]==20. and newmeta["parameter_identity"]["duration_ps"]==new_duration
 assert newmeta["parameter_identity"]["ceiling_ps"]==200.
 assert newmeta["controller"]==h.metadata(a)["controller"]
 assert newmeta["source"]==h.metadata(a)["source"]
 for k in newmeta["parameter_identity"]:
  if k!="duration_ps":assert newmeta["parameter_identity"][k]==h.metadata(a)["parameter_identity"][k]
 assert newmeta["stable_accepted_steps"]==7 and newmeta["consecutive_rejections"]==2
 manifest=json.loads((tmp_path/"new/manifest.json").read_text())
 assert not manifest["run_complete"] and manifest["explicit_horizon_change"]==change
 assert manifest["observation"]["development_duration_s"]==(None if new_duration is None else new_duration*1e-12)
 assert r["numerical_segments"][0]["duration_ps"]==20.
 assert r["numerical_segments"][1]["duration_ps"]==new_duration
 argv=r["resume_argv"]
 assert ("--duration-ps" not in argv) if new_duration is None else argv[argv.index("--duration-ps")+1]=="40.0"

@pytest.mark.parametrize("requested",[19.,20.,float("nan"),float("inf"),True])
def test_completed_control_rejects_invalid_duration(tmp_path,requested):
 run,repo,a,m=completed_control_fixture(tmp_path)
 with pytest.raises(ValueError):h.extend_control_horizon(m,h.metadata(a),requested)

@pytest.mark.parametrize("fault",["not_complete","wrong_status","confirmed","terminal","missing_duration","wrong_clock","implicit","latent_terminal"])
def test_completed_control_extension_is_narrow_and_never_reopens_event(tmp_path,fault):
 run,repo,a,m=completed_control_fixture(tmp_path);meta=h.metadata(a);requested=40.
 if fault=="not_complete":m["run_complete"]=False
 if fault=="wrong_status":m["status"]="DC_HOLD_COMPLETED"
 if fault=="confirmed":meta["controller"]["state"]["status"]="TRIGGER_CONFIRMED"
 if fault=="terminal":meta["controller"]["state"]["status"]="TIME_CENSORED_NO_TRIGGER"
 if fault=="missing_duration":meta["parameter_identity"]["duration_ps"]=None
 if fault=="wrong_clock":meta["time_ps"]=19.99
 if fault=="implicit":requested="unchanged"
 if fault=="latent_terminal":meta["controller"]["state"]["actual_stop_s"]=20e-12
 with pytest.raises(ValueError):h.extend_control_horizon(m,meta,requested)

@pytest.mark.parametrize("terminal",["TRIGGER_MARGIN_COMPLETE","TIME_CENSORED_NO_TRIGGER","NUMERICAL_FAILURE","INVALID_INITIAL_EVENT"])
def test_no_terminal_event_ever_transfers_even_without_horizon_override(tmp_path,terminal):
 run,repo,a,m=fixture(tmp_path);meta=h.metadata(a)
 meta["controller"]["state"]["status"]=terminal
 with pytest.raises(ValueError):h.extend_control_horizon(m,meta,"unchanged")

def test_cli_horizon_argument_is_explicit():
 assert h.duration_argument("none") is None
 assert h.duration_argument("None") is None
 assert h.duration_argument("40") == 40.
 assert h.duration_argument("unchanged") == "unchanged"
 with pytest.raises(ValueError):h.duration_argument("reset")
