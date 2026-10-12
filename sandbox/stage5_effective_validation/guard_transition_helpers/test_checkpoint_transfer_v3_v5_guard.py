"""Guard activation and numerical-segment provenance; synthetic only."""
import copy,importlib.util,json,types
from pathlib import Path
from unittest.mock import patch
import numpy as np
import pytest
HERE=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location("_guard_legacy_fixtures",HERE/"test_checkpoint_transfer_v3_v5_legacy.py")
fixtures=importlib.util.module_from_spec(spec);spec.loader.exec_module(fixtures)
h=fixtures.h

class GuardRunner(fixtures.FakeRunner):
 @staticmethod
 def build_system(a,output,resume=False):
  system,dc,g,material,params=fixtures.FakeRunner.build_system(a,output,resume)
  system.kinetic_residual_energy_scale_J=1.
  return system,dc,g,material,params

def guarded_repo(repo):
 p=repo/"pipelines/05_effective_photon_run.py"
 p.write_text(p.read_text()+"\n# supports --local-energy-tolerance\n")

def transfer_fixture(tmp_path,*,tolerance=1e-5):
 run,repo,arrays,manifest=fixtures.fixture(tmp_path,counters=True);guarded_repo(repo)
 with patch.object(h,"SCRATCH",tmp_path),patch.object(h,"load_runner",lambda p:GuardRunner):
  h.capture(run,tmp_path/"anchor")
  r=h.transfer(tmp_path/"anchor",tmp_path/"new",repo,.01,2e-7,2,12,local_energy_tolerance=tolerance)
 return run,repo,arrays,manifest,r

def test_explicit_guard_activation_preserves_every_physical_array_and_signed_ledger(tmp_path):
 run,repo,arrays,parent,r=transfer_fixture(tmp_path)
 assert r["state_arrays_sha256_before"]==r["state_arrays_sha256_after"]
 assert r["conserved_metadata_before"]==r["conserved_metadata_after"]
 change=r["energy_guard_change"]
 assert change["enabled_before"] is False and change["enabled_after"] is True
 assert change["local_energy_tolerance_before"] is None
 assert change["local_energy_tolerance_after"]==1e-5 and change["effective_time_ps"]==.02
 assert change["reference_energy_J"]==1. and not change["pre_guard_history_certified"]
 assert change["history_diagnostic_measured"]["absolute_external_energy_integral_J"]==1.
 assert change["history_diagnostic_measured"]["sum_abs_local_energy_defect_J"]==0.
 meta=h.metadata(h.decode((tmp_path/"new/checkpoint_latest.npz").read_bytes()))
 assert meta["parameter_identity"]["local_energy_tolerance"]==1e-5
 assert meta["guard_active_since_time_ps"]==.02
 assert meta["max_guarded_accepted_local_energy_fraction"]==0.
 assert meta["absolute_external_energy_integral_J"]==1.
 assert "pre_guard_history" in meta["energy_guard_ledger_origin"]
 assert r["numerical_segments"][0]["energy_guard_enabled"] is False
 assert r["numerical_segments"][0]["local_energy_tolerance"] is None
 assert r["numerical_segments"][1]["energy_guard_enabled"] is True
 assert r["numerical_segments"][1]["local_energy_tolerance"]==1e-5
 assert "--local-energy-tolerance" in r["resume_argv"]
 assert set(r["frontend_changes"])<=h.FRONTEND and r["RHS_sources_identical"]
 assert (tmp_path/"new/provenance/parent_checkpoint.npz").read_bytes()==(run/"checkpoint_latest.npz").read_bytes()

def test_guarded_runner_requires_explicit_tolerance_before_making_output(tmp_path):
 run,repo,arrays,m=fixtures.fixture(tmp_path);guarded_repo(repo)
 with patch.object(h,"SCRATCH",tmp_path):
  h.capture(run,tmp_path/"anchor")
  with pytest.raises(ValueError,match="requires explicit"):
   h.transfer(tmp_path/"anchor",tmp_path/"new",repo,.01,2e-7,2,12)
 assert not (tmp_path/"new").exists()

@pytest.mark.parametrize("tol",[0.,-1.,float("nan"),float("inf")])
def test_invalid_guard_tolerance_is_rejected(tmp_path,tol):
 run,repo,arrays,m=fixtures.fixture(tmp_path);guarded_repo(repo)
 with patch.object(h,"SCRATCH",tmp_path):
  h.capture(run,tmp_path/"anchor")
  with pytest.raises(ValueError,match="tolerance"):
   h.transfer(tmp_path/"anchor",tmp_path/"new",repo,.01,2e-7,2,12,local_energy_tolerance=tol)
 assert not (tmp_path/"new").exists()

def test_guard_reference_is_checked_against_reconstructed_system(tmp_path):
 class WrongReference(GuardRunner):
  @staticmethod
  def build_system(a,output,resume=False):
   system,dc,g,material,params=GuardRunner.build_system(a,output,resume)
   system.kinetic_residual_energy_scale_J=2.
   return system,dc,g,material,params
 run,repo,arrays,m=fixtures.fixture(tmp_path);guarded_repo(repo)
 with patch.object(h,"SCRATCH",tmp_path),patch.object(h,"load_runner",lambda p:WrongReference):
  h.capture(run,tmp_path/"anchor")
  with pytest.raises(ValueError,match="reference differs"):
   h.transfer(tmp_path/"anchor",tmp_path/"new",repo,.01,2e-7,2,12,local_energy_tolerance=1e-5)

def test_prior_numerical_segments_and_guard_ledger_remain_explicit(tmp_path):
 run,repo,arrays,m,first=transfer_fixture(tmp_path)
 with patch.object(h,"SCRATCH",tmp_path),patch.object(h,"load_runner",lambda p:GuardRunner):
  cr=h.capture(tmp_path/"new",tmp_path/"anchor_second")
  second=h.transfer(tmp_path/"anchor_second",tmp_path/"second",repo,.005,2e-7,1,12,local_energy_tolerance=1e-5)
 old=first["numerical_segments"]
 assert second["numerical_segments"][0]==old[0]
 assert len(second["numerical_segments"])==3
 assert second["numerical_segments"][1]["start_time_ps"]==old[1]["start_time_ps"]
 assert second["numerical_segments"][1]["end_time_ps"]==.02
 assert second["numerical_segments"][1]["local_energy_tolerance"]==1e-5
 assert second["numerical_segments"][2]["dt_ps"]==.005
 assert second["energy_guard_change"]["enabled_before"] is True
 assert second["energy_guard_change"]["ledger_after"]==first["energy_guard_change"]["ledger_after"]
 prior=tmp_path/"anchor_second/provenance/parent_transfer_receipt.json"
 assert h.filehash(prior)==h.filehash(tmp_path/"new/provenance/transfer_receipt.json")
 assert "provenance/parent_transfer_receipt.json" in cr["physical_file_sha256"]

def test_tolerance_change_resets_only_guarded_max_window_not_signed_work(tmp_path):
 run,repo,arrays,m,first=transfer_fixture(tmp_path)
 with patch.object(h,"SCRATCH",tmp_path),patch.object(h,"load_runner",lambda p:GuardRunner):
  h.capture(tmp_path/"new",tmp_path/"anchor_second")
  r=h.transfer(tmp_path/"anchor_second",tmp_path/"second",repo,.01,2e-7,2,12,local_energy_tolerance=5e-6)
 assert r["energy_guard_change"]["local_energy_tolerance_before"]==1e-5
 assert r["energy_guard_change"]["local_energy_tolerance_after"]==5e-6
 assert r["energy_guard_change"]["ledger_after"]["absolute_external_energy_integral_J"]==1.
 assert r["conserved_metadata_before"]==r["conserved_metadata_after"]
 assert "tolerance_changed" in r["energy_guard_change"]["ledger_after"]["energy_guard_ledger_origin"]

def test_old_source_and_new_guard_are_never_confused(tmp_path):
 run,repo,arrays,m=fixtures.fixture(tmp_path)
 with patch.object(h,"SCRATCH",tmp_path):
  h.capture(run,tmp_path/"anchor")
  with pytest.raises(ValueError,match="does not support"):
   h.transfer(tmp_path/"anchor",tmp_path/"new",repo,.01,2e-7,2,12,local_energy_tolerance=1e-5)
 assert not (tmp_path/"new").exists()

def test_diagnostic_actual_threads_and_helper_provenance_are_segmented(tmp_path):
 run,repo,arrays,m=fixtures.fixture(tmp_path);guarded_repo(repo)
 meta=h.metadata(arrays);meta["parameter_identity"]["max_iterations"]=40
 arrays["metadata_json"]=np.array(json.dumps(meta));np.savez_compressed(run/"checkpoint_latest.npz",**arrays)
 m["parameter_identity"]=meta["parameter_identity"];h.writejson(run/"manifest.json",m)
 diagnostic=dict(checkpoint_sha256=h.filehash(run/"checkpoint_latest.npz"),
  metadata_runtime_note="ancestral threads4; actual diagnostic threads1 CPU13",
  executed_A_numeric_settings=dict(start_time_ps=.01,end_time_ps=.02,max_step_ps=.02,
    step_tolerance=2e-9,max_iterations=40,threads=1,one_physical_CPU=13,BLAS_threads=1))
 h.writejson(run/"parent_preparation_receipt.json",diagnostic)
 (run/"diagnostic_provenance").mkdir();(run/"diagnostic_provenance/executed_helper.py").write_text("EXACT ORIGINAL")
 with patch.object(h,"SCRATCH",tmp_path),patch.object(h,"load_runner",lambda p:GuardRunner):
  cr=h.capture(run,tmp_path/"anchor")
  r=h.transfer(tmp_path/"anchor",tmp_path/"new",repo,.01,2e-7,2,12,local_energy_tolerance=1e-5)
 assert len(r["numerical_segments"])==3
 assert r["numerical_segments"][0]["end_time_ps"]==.01
 diag=r["numerical_segments"][1]
 assert diag["start_time_ps"]==.01 and diag["end_time_ps"]==.02
 assert diag["threads"]==1 and diag["actual_affinity_cpus"]==[13] and diag["step_tolerance"]==2e-9
 assert r["numerical_segments"][2]["threads"]==2 and r["numerical_segments"][2]["step_tolerance"]==2e-7
 assert h.filehash(tmp_path/"new/provenance/parent_preparation_receipt.json")==h.filehash(run/"parent_preparation_receipt.json")
 assert (tmp_path/"new/provenance/diagnostic_provenance/executed_helper.py").read_text()=="EXACT ORIGINAL"
 assert "parent_preparation_receipt.json" in cr["physical_file_sha256"]
