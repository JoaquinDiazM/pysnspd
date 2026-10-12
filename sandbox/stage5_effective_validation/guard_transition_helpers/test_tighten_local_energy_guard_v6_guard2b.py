"""Synthetic prospective guard transactions: no physics evolution or live capture."""
import copy,importlib.util,json
from pathlib import Path
from unittest.mock import patch
import numpy as np
import pytest
HERE=Path(__file__).parent
def load(name):
 s=importlib.util.spec_from_file_location(name,HERE/(name+".py"));m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m
t=load("tighten_local_energy_guard_v6_guard2b")
fixtures=load("test_extend_completed_control_v6")
def fixture(tmp_path,kind="fresh",time=21.375,oldtol=1e-5):
 run,a,parent,meta,log=(fixtures.fresh_fixture if kind=="fresh" else fixtures.inherited_fixture)(tmp_path)
 meta.update(time_ps=time,local_energy_rejections=17,stable_accepted_steps=7,consecutive_rejections=2,accepted_steps=2,total_rejections=19)
 meta["parameter_identity"]["local_energy_tolerance"]=oldtol
 meta["controller"]["state"].update(last_time_s=time*1e-12,sample_count=3)
 anchor=meta["guard_active_since_time_ps"]
 middle=time/2 if anchor==0 else anchor
 for key,value in list(a.items()):
  if key.startswith("trajectory__"):a[key]=np.array([value[0],.5*(value[0]+value[-1]),value[-1]])
 a["trajectory__time_s"]=np.array([0.,middle*1e-12,time*1e-12])
 delta=oldtol*(1.5 if anchor==0 else 2.)*.9
 a["trajectory__energy_total_J"]=np.array([10.,10.5-delta,11.])
 a["trajectory__external_energy_integral_J"]=np.array([0.,.5,1.])
 a["trajectory__energy_defect_J"]=np.array([0.,-delta,0.])
 a["snapshots__time_s"][-1]=time*1e-12
 dw=np.diff(a["trajectory__external_energy_integral_J"]);defect=np.diff(a["trajectory__energy_total_J"])-dw
 fractions=abs(defect)/(1.+np.cumsum(abs(dw)));post=a["trajectory__time_s"][:-1]>=anchor*1e-12
 meta.update(absolute_external_energy_integral_J=float(sum(abs(dw))),sum_abs_local_energy_defect_J=float(sum(abs(defect))),
  max_accepted_local_energy_fraction=float(max(fractions)),max_guarded_accepted_local_energy_fraction=float(max(fractions[post])))
 parent.update(run_complete=False,status="RUNNING",parameter_identity=copy.deepcopy(meta["parameter_identity"]))
 parent["numerical"].update(local_energy_tolerance=oldtol,local_energy_guard={"enabled":True})
 parent["local_energy_control"]={k:meta[k] for k in t.h.GUARD_KEYS}
 a["metadata_json"]=np.array(json.dumps(meta));fixtures.save(run,a,parent)
 return run,a,parent,meta,log
def sandbox(tmp_path):
 return patch.multiple(t.h,SCRATCH=tmp_path,load_runner=lambda p:fixtures.FakeRunner)
@pytest.mark.parametrize("kind",["fresh","inherited"])
@pytest.mark.parametrize("newtol",[2.5e-6,3.125e-7])
@pytest.mark.parametrize("time",[20.,21.375])
def test_every_array_counter_global_ledger_and_policy_except_declared_guard_window_is_exact(tmp_path,kind,newtol,time):
 with patch.object(fixtures.v.h.lineage,"SCRATCH",tmp_path):
  run,a,parent,old,log=fixture(tmp_path,kind,time)
 before=(run/"checkpoint_latest.npz").read_bytes()
 with sandbox(tmp_path),patch.object(t.h.lineage,"SCRATCH",tmp_path):
  t.capture(run,tmp_path/"anchor",fresh_start_log=log)
  receipt=t.tighten(tmp_path/"anchor",tmp_path/"new",newtol)
 assert (run/"checkpoint_latest.npz").read_bytes()==before
 after=t.h.decode((tmp_path/"new/checkpoint_latest.npz").read_bytes());meta=t.h.metadata(after)
 for k,value in a.items():
  if k!="metadata_json":assert t.h.arrayhash(value)==t.h.arrayhash(after[k])
 for k,value in old.items():
  if k not in t.WINDOW_KEYS|{"parameter_identity"}:assert meta[k]==value
 assert {k:v for k,v in old["parameter_identity"].items() if k!=t.IDENTITY_KEY}=={k:v for k,v in meta["parameter_identity"].items() if k!=t.IDENTITY_KEY}
 assert meta["parameter_identity"][t.IDENTITY_KEY]==newtol
 assert meta["guard_active_since_time_ps"]==time and meta["max_guarded_accepted_local_energy_fraction"]==0.
 assert receipt["global_energy_error_before"]==receipt["global_energy_error_after"]
 assert receipt["guard_ledger_before"]["max_guarded_accepted_local_energy_fraction"]==old["max_guarded_accepted_local_energy_fraction"]
 assert receipt["energy_guard_change"]["previous_guard_active_since_time_ps"]==old["guard_active_since_time_ps"]
 assert receipt["energy_guard_change"]["prospective_only"]
 assert receipt["energy_guard_change"]["prior_steps_certified_under_new_tolerance"] is False
 assert receipt["pins_before"]==receipt["pins_after"] and len(receipt["pins_after"])==13
 assert not receipt["frontend_changes"] and not receipt["implementation_changes"]
 assert receipt["resources_before"]==receipt["resources_after"]
 assert receipt["conserved_metadata_before"]==receipt["conserved_metadata_after"]
 assert receipt["next_step_ps_before"]==receipt["next_step_ps_after"]
 assert receipt["physical_identity_verified"]["no_collision_RHS_evaluated"]
 assert receipt["physical_identity_verified"]["photon_source_reapplied"] is False
 assert receipt["numerical_segments"][-2]["end_time_ps"]==time
 assert receipt["numerical_segments"][-2]["local_energy_tolerance"]==1e-5
 assert receipt["numerical_segments"][-1]["start_time_ps"]==time
 assert receipt["numerical_segments"][-1]["guard_active_since_time_ps"]==time
 assert receipt["numerical_segments"][-1]["local_energy_tolerance"]==newtol
 assert receipt["numerical_segments"][-1]["duration_ps"]==20.
 assert receipt["resume_argv"][receipt["resume_argv"].index("--local-energy-tolerance")+1]==str(newtol)
 with patch.object(t.h.lineage,"SCRATCH",tmp_path):
  verified,_,_,_=t.h.lineage.verify_receipt(tmp_path/"new/provenance/transfer_receipt.json",receipt["pins_after"])
  resolved,proofs=t.h.lineage.resolve_receipt(tmp_path/"new/provenance/transfer_receipt.json",receipt["pins_after"])
 assert verified["energy_guard_change"]["effective_time_ps"]==time
 assert resolved[-1]["start_time_ps"]==time and resolved[-1]["local_energy_tolerance"]==newtol
 assert resolved[-2]["local_energy_tolerance"]==1e-5
@pytest.mark.parametrize("newtol",[None,0.,-1.,1e-5,2e-5,float("inf"),float("nan"),True,"2.5e-6"])
def test_only_strictly_positive_strictly_smaller_guard_is_accepted(tmp_path,newtol):
 run,a,parent,meta,log=fixture(tmp_path)
 with pytest.raises(ValueError):t.guard_change(parent,meta,a,newtol)
@pytest.mark.parametrize("fault",["complete","terminal","candidate","crossing","confirmation","stop_target","actual_stop","failure","past_crossing","history_crossing","old_guard_missing","old_guard_future","guard_ledger_missing","old_guard_violation","pin_missing","resource_changed"])
def test_ineligible_state_cannot_be_reinterpreted(tmp_path,fault):
 run,a,parent,meta,log=fixture(tmp_path)
 if fault=="complete":parent["run_complete"]=True
 elif fault=="terminal":meta["controller"]["state"]["status"]="TIME_CENSORED_NO_TRIGGER"
 elif fault in ("candidate","crossing","confirmation","stop_target","actual_stop"):
  key={"candidate":"candidate_crossing_s","crossing":"crossing_time_s","confirmation":"confirmation_time_s","stop_target":"stop_target_s","actual_stop":"actual_stop_s"}[fault]
  meta["controller"]["state"][key]=20e-12
 elif fault=="failure":meta["controller"]["state"]["failure_reason"]="error"
 elif fault=="past_crossing":meta["controller"]["state"]["rejected_candidate_count"]=1
 elif fault=="history_crossing":a["trajectory__Vout_V"][0]=1e-3
 elif fault=="old_guard_missing":meta["parameter_identity"].pop("local_energy_tolerance")
 elif fault=="old_guard_future":meta["guard_active_since_time_ps"]=22.
 elif fault=="guard_ledger_missing":meta.pop("sum_abs_local_energy_defect_J")
 elif fault=="old_guard_violation":meta["max_guarded_accepted_local_energy_fraction"]=1.1e-5
 elif fault=="pin_missing":meta["source_pins_sha256"].pop(next(iter(meta["source_pins_sha256"])))
 elif fault=="resource_changed":parent["resources"]["actual_threads"]=1
 with pytest.raises((ValueError,KeyError)):t.guard_change(parent,meta,a,2.5e-6)
@pytest.mark.parametrize("fault",["runtime","source","checkpoint","manifest","lineage"])
def test_tampered_anchor_rejected_before_new_output(tmp_path,fault):
 run,a,parent,meta,log=fixture(tmp_path)
 with sandbox(tmp_path):
  t.capture(run,tmp_path/"anchor",fresh_start_log=log)
  anchor=tmp_path/"anchor"
  if fault=="runtime":path=anchor/"runtime_sources"/next(iter(meta["source_pins_sha256"]))
  else:path=anchor/{"source":"source_initial_state.npz","checkpoint":"checkpoint_original.npz","manifest":"manifest_original.json","lineage":"provenance/guard_parent_lineage.json"}[fault]
  path.write_bytes(b"tampered")
  with pytest.raises((ValueError,json.JSONDecodeError)):t.tighten(anchor,tmp_path/"new",2.5e-6)
 assert not (tmp_path/"new").exists()
def test_existing_destination_preserved(tmp_path):
 run,a,parent,meta,log=fixture(tmp_path)
 with sandbox(tmp_path):
  t.capture(run,tmp_path/"anchor",fresh_start_log=log)
  (tmp_path/"new").mkdir();(tmp_path/"new/keep").write_text("original")
  with pytest.raises(FileExistsError):t.tighten(tmp_path/"anchor",tmp_path/"new",2.5e-6)
 assert (tmp_path/"new/keep").read_text()=="original"
def test_pair_ratio_and_real_anchors_are_independent_prospective_policies(tmp_path):
 npath=tmp_path/"nominal";npath.mkdir();tpath=tmp_path/"temporal";tpath.mkdir()
 _,a,parent,meta,_=fixture(npath,time=21.375)
 _,ta,tp,tm,_=fixture(tpath,time=19.125,oldtol=1.25e-6)
 nominal_ledger,nominal=t.guard_change(parent,meta,a,2.5e-6)
 temporal_ledger,temporal=t.guard_change(tp,tm,ta,3.125e-7)
 assert nominal["local_energy_tolerance_after"]/temporal["local_energy_tolerance_after"]==8.
 assert nominal_ledger["guard_active_since_time_ps"]==21.375
 assert temporal_ledger["guard_active_since_time_ps"]==19.125
 assert nominal["effective_time_ps"]!=temporal["effective_time_ps"]
def test_fresh_proof_preserves_actual_duration_not_a_template_default(tmp_path):
 run,a,parent,meta,log=fixture(tmp_path)
 meta["parameter_identity"]["duration_ps"]=None;parent["parameter_identity"]["duration_ps"]=None
 a["metadata_json"]=np.array(json.dumps(meta));fixtures.save(run,a,parent)
 with sandbox(tmp_path):
  t.capture(run,tmp_path/"anchor",fresh_start_log=log)
  r=t.tighten(tmp_path/"anchor",tmp_path/"new",2.5e-6)
 assert all(s["duration_ps"] is None for s in r["numerical_segments"])
 assert "--duration-ps" not in r["resume_argv"]
def test_global_reference_and_work_origin_never_restart(tmp_path):
 run,a,parent,meta,log=fixture(tmp_path)
 before=t.global_error_at_anchor(meta,a)
 ledger,change=t.guard_change(parent,meta,a,2.5e-6)
 after=t.global_error_at_anchor(dict(meta,**ledger),a)
 assert before==after and before["global_origin_time_ps"]==0.
 assert change["global_error_origin_reset"] is False
 assert change["global_denominator_reset"] is False

@pytest.mark.parametrize("fault",["absolute_work","absolute_defect","maximum_all","maximum_guarded","hidden_history_violation","anchor_not_saved","manifest_tolerance","manifest_ledger"])
def test_falsified_old_guard_cannot_be_closed_even_with_plausible_saved_fields(tmp_path,fault):
 run,a,parent,meta,log=fixture(tmp_path)
 if fault=="absolute_work":meta["absolute_external_energy_integral_J"]+=.01
 elif fault=="absolute_defect":meta["sum_abs_local_energy_defect_J"]+=1e-3
 elif fault=="maximum_all":meta["max_accepted_local_energy_fraction"]*=.5
 elif fault=="maximum_guarded":meta["max_guarded_accepted_local_energy_fraction"]*=.5
 elif fault=="hidden_history_violation":
  a["trajectory__energy_total_J"][1]-=.1
  # Deliberately leave plausible reported maxima below tolerance.
 elif fault=="anchor_not_saved":meta["guard_active_since_time_ps"]=.1234567
 elif fault=="manifest_tolerance":parent["numerical"]["local_energy_tolerance"]=2e-5
 elif fault=="manifest_ledger":parent["local_energy_control"]["sum_abs_local_energy_defect_J"]+=.001
 if fault!="manifest_ledger":parent["local_energy_control"]={k:meta[k] for k in t.h.GUARD_KEYS}
 with pytest.raises(ValueError,match="Old guard history audit failed"):
  t.guard_change(parent,meta,a,2.5e-6)

@pytest.mark.parametrize("count",[-1,True])
def test_old_guard_rejection_count_nonnegative_integer_not_boolean(tmp_path,count):
 run,a,parent,meta,log=fixture(tmp_path)
 meta["local_energy_rejections"]=count
 parent["local_energy_control"]={k:meta[k] for k in t.h.GUARD_KEYS}
 with pytest.raises(ValueError,match="rejection_count_valid"):
  t.validate_old_guard_history(parent,meta,a)
 with pytest.raises(ValueError):t.guard_change(parent,meta,a,2.5e-6)

def test_new_segment_reanchors_guard_but_old_segment_is_never_rewritten():
 meta={"time_ps":22.435,"parameter_identity":{"dt_ps":.02,"step_tolerance":2e-7,"threads":4,"local_energy_tolerance":1e-5,"duration_ps":None}}
 old_segment={"start_time_ps":0.,"end_time_ps":None,"dt_ps":.02,"step_tolerance":2e-7,"threads":4,"local_energy_tolerance":1e-5,"duration_ps":None,"guard_active_since_time_ps":12.4,"origin":"historical"}
 proof={"numerical_segments":[old_segment]}
 result=t.segments_for_change(proof,meta,dict(meta["parameter_identity"],local_energy_tolerance=2.5e-6))
 assert result[-1]["guard_active_since_time_ps"]==22.435
 assert result[-2]["guard_active_since_time_ps"]==12.4
 assert old_segment["end_time_ps"] is None and proof["numerical_segments"]==[old_segment]
