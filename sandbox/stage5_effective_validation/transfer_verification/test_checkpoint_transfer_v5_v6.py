"""Exact storage-only handoff: synthetic, no live jobs or physical evolution."""
import copy,importlib.util,json,shutil
from pathlib import Path
from unittest.mock import patch
import numpy as np
import pytest
HERE=Path(__file__).parent
def load(n):
 s=importlib.util.spec_from_file_location(n,HERE/(n+".py"));m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m
f=load("test_checkpoint_transfer_v3_v5_legacy")
g=load("test_checkpoint_transfer_v3_v5_guard")
h=load("checkpoint_transfer_v5_v6")

def fixture(tmp_path):
 run,repo,arrays,manifest=f.fixture(tmp_path,counters=True)
 # Insert the real two implementation paths into the 13-pin synthetic namespace.
 oldpins=manifest["source_pins_sha256"];newpins=copy.deepcopy(oldpins)
 for i,name in enumerate(sorted(h.IMPLEMENTATION_PATHS)):
  prior="pysnspd/experimental/core_%02d.py"%i
  for root in (run/"runtime_sources",repo):
   (root/name).parent.mkdir(parents=True,exist_ok=True);(root/prior).rename(root/name)
  newpins[name]=newpins.pop(prior)
 meta=f.h.metadata(arrays);meta["source_pins_sha256"]=newpins;manifest["source_pins_sha256"]=newpins
 arrays["metadata_json"]=np.array(json.dumps(meta))
 np.savez_compressed(run/"checkpoint_latest.npz",**arrays);f.h.writejson(run/"manifest.json",manifest)
 arc=f.h.readjson(run/"runtime_sources/archive_manifest.json") if hasattr(f.h,"readjson") else json.loads((run/"runtime_sources/archive_manifest.json").read_text())
 arc["runtime_pins_sha256"]=newpins
 arc["package_sources_and_small_geometry"]={k:{"sha256":v} for k,v in newpins.items()}
 f.h.writejson(run/"runtime_sources/archive_manifest.json",arc)
 g.guarded_repo(repo)
 with patch.object(f.h,"SCRATCH",tmp_path),patch.object(f.h,"load_runner",lambda p:g.GuardRunner):
  f.h.capture(run,tmp_path/"initial_anchor")
  r=f.h.transfer(tmp_path/"initial_anchor",tmp_path/"v5",repo,.02,2e-7,2,12,local_energy_tolerance=1e-5)
 current=tmp_path/"v5"
 # Emulate the resumed runner dropping only its descriptive receipt pointer.
 m=json.loads((current/"manifest.json").read_text());m.pop("explicit_transfer_receipt",None);f.h.writejson(current/"manifest.json",m)
 target=tmp_path/"v6_repo";shutil.copytree(current/"runtime_sources",target)
 pins_before=r["pins_after"];pins_after=copy.deepcopy(pins_before)
 for name in h.IMPLEMENTATION_PATHS:
  (target/name).write_text((target/name).read_text()+"\n# exact storage candidate\n")
  pins_after[name]=h.filehash(target/name)
 evidence=[];cases=[]
 for t in (2,4):
  oa=tmp_path/("oracle_n"+str(t));ca=tmp_path/("candidate_n"+str(t));oa.mkdir();ca.mkdir()
  for q in (oa,ca):
   np.savez_compressed(q/"DC_limited.npz",rhs=np.array([1.,2.]),diagnostics=np.array([3.]))
   f.h.writejson(q/"receipt.json",dict(status="COMPLETED"))
   for z in q.iterdir():evidence.append(dict(path=str(z),sha256=h.filehash(z)))
  cases.append(dict(case="DC_limited",thread_count=t,status="RESULT",passed=True,
   bitwise_equal=True,numerical_arrays_and_scalars_equal=True,oracle_path=str(oa/"DC_limited.npz"),
   candidate_path=str(ca/"DC_limited.npz"),oracle_npz_sha256=h.filehash(oa/"DC_limited.npz"),
   candidate_npz_sha256=h.filehash(ca/"DC_limited.npz"),
   oracle_receipt_path=str(oa/"receipt.json"),candidate_receipt_path=str(ca/"receipt.json")))
 certificate=dict(schema="stage5-node-storage-exact-equivalence-v1",
  sources_before={k:pins_before[k] for k in h.IMPLEMENTATION_PATHS},
  sources_after={k:pins_after[k] for k in h.IMPLEMENTATION_PATHS},
  full_RHS_and_diagnostics_bitwise_equal=True,numerical_core_or_physics_changed=False,
  evidence_files=evidence,cases=cases)
 cp=tmp_path/"certificate.json";f.h.writejson(cp,certificate)
 return current,target,cp,certificate

def prepare(tmp_path):
 current,target,cp,c=fixture(tmp_path)
 with patch.object(h,"SCRATCH",tmp_path):
  r=h.capture(current,tmp_path/"anchor")
 return current,target,cp,c,r

def execute(tmp_path,target,cp,**kwargs):
 args=dict(local_energy_tolerance=1e-5,certificate_path=cp,approved_equivalence_sha256=h.filehash(cp))
 args.update(kwargs)
 with patch.object(h,"SCRATCH",tmp_path),patch.object(h,"load_runner",lambda p:g.GuardRunner):
  return h.transfer(tmp_path/"anchor",tmp_path/"v6",target,.02,2e-7,2,12,**args)

def test_bitarrays_ledgers_controller_history_counters_guard_and_sources(tmp_path):
 current,target,cp,c,cr=prepare(tmp_path)
 original=(current/"checkpoint_latest.npz").read_bytes();receipt=execute(tmp_path,target,cp)
 assert (current/"checkpoint_latest.npz").read_bytes()==original
 assert receipt["state_arrays_sha256_before"]==receipt["state_arrays_sha256_after"]
 assert receipt["conserved_metadata_before"]==receipt["conserved_metadata_after"]
 assert receipt["guard_ledger_sha256_before"]==receipt["guard_ledger_sha256_after"]
 assert receipt["RHS_sources_identical"] is False and receipt["physical_RHS_bitwise_equivalence_certified"]
 assert set(receipt["implementation_changes"])==h.IMPLEMENTATION_PATHS and receipt["frontend_changes"]=={}
 before=h.decode(original);after=h.decode((tmp_path/"v6/checkpoint_latest.npz").read_bytes())
 for k in before:
  if k!="metadata_json":assert np.array_equal(before[k],after[k])
 mb=h.metadata(before);ma=h.metadata(after)
 assert ma["stable_accepted_steps"]==mb["stable_accepted_steps"]
 assert ma["consecutive_rejections"]==mb["consecutive_rejections"]
 assert all(ma[k]==mb[k] for k in h.GUARD_KEYS)
 assert receipt["numerical_segments"][0]["step_tolerance"]==2e-9
 assert receipt["numerical_segments"][-1]["origin"]=="v6_exact_storage_implementation"
 assert (tmp_path/"v6/provenance/parent_transfer_receipt.json").exists()
 assert (tmp_path/"v6/provenance/exact_implementation_equivalence_certificate.json").exists()
 assert receipt["physical_identity_verified"]["photon_source_reapplied"] is False

def test_capture_fallback_and_supplement_bind_original_receipt(tmp_path):
 current,target,cp,c,r=prepare(tmp_path)
 p=tmp_path/"anchor/provenance"
 assert (p/"parent_transfer_receipt.json").read_bytes()==(current/"provenance/transfer_receipt.json").read_bytes()
 s=json.loads((p/"numerical_lineage_supplement_v6.json").read_text())
 assert s["numerical_lineage_complete"]
 assert s["original_transfer_receipt_sha256"]==h.filehash(p/"parent_transfer_receipt.json")

@pytest.mark.parametrize("which",["approved_sha","missing_certificate","bad_evidence","bad_source_map","physics_flag","missing_n4"])
def test_equivalence_fail_closed_before_new_output(tmp_path,which):
 current,target,cp,c,_=prepare(tmp_path);args={}
 if which=="approved_sha":args["approved_equivalence_sha256"]="0"*64
 elif which=="missing_certificate":args["certificate_path"]=None
 elif which=="bad_evidence":Path(c["evidence_files"][0]["path"]).write_bytes(b"changed")
 elif which=="bad_source_map":c["sources_after"][next(iter(h.IMPLEMENTATION_PATHS))]="0"*64
 elif which=="physics_flag":c["numerical_core_or_physics_changed"]=True
 elif which=="missing_n4":c["cases"]=[x for x in c["cases"] if x["thread_count"]==2]
 if which in ("bad_source_map","physics_flag","missing_n4"):f.h.writejson(cp,c)
 with pytest.raises(ValueError):execute(tmp_path,target,cp,**args)
 assert not (tmp_path/"v6").exists()

@pytest.mark.parametrize("key,value",[("dt_ps",.01),("tolerance",2e-9),("threads",4),("local_energy_tolerance",1.25e-6),("duration_ps_override",None)])
def test_numerical_policy_may_not_change(tmp_path,key,value):
 current,target,cp,c,_=prepare(tmp_path)
 if key in ("local_energy_tolerance","duration_ps_override"):
  with pytest.raises(ValueError):execute(tmp_path,target,cp,**{key:value})
 else:
  vals=dict(dt_ps=.02,tolerance=2e-7,threads=2);vals[key]=value
  with patch.object(h,"SCRATCH",tmp_path),pytest.raises(ValueError):
   h.transfer(tmp_path/"anchor",tmp_path/"v6",target,vals["dt_ps"],vals["tolerance"],vals["threads"],12,
    local_energy_tolerance=1e-5,certificate_path=cp,approved_equivalence_sha256=h.filehash(cp))
 assert not (tmp_path/"v6").exists()

def test_frontend_or_other_rhs_change_is_rejected(tmp_path):
 current,target,cp,c,_=prepare(tmp_path)
 p=target/"pipelines/05_effective_photon_run.py";p.write_text(p.read_text()+"# forbidden\n")
 with pytest.raises(ValueError,match="two"):execute(tmp_path,target,cp)
 assert not (tmp_path/"v6").exists()

def test_mutated_supplement_is_rejected_by_capture_hash(tmp_path):
 current,target,cp,c,_=prepare(tmp_path)
 p=tmp_path/"anchor/provenance/numerical_lineage_supplement_v6.json";p.write_text("{}")
 with pytest.raises(ValueError,match="Anchored physical file"):execute(tmp_path,target,cp)
 assert not (tmp_path/"v6").exists()
