import copy,importlib.util,json,tempfile,unittest
from pathlib import Path
import numpy as np
spec=importlib.util.spec_from_file_location("lineage",Path(__file__).with_name("numerical_lineage_v6.py"))
m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
class LineageTest(unittest.TestCase):
 def setUp(self):
  self.t=tempfile.TemporaryDirectory(dir="/home/jdiaz/scratch/tmp/pysnspd")
  self.root=Path(self.t.name);self.count=0
 def tearDown(self):self.t.cleanup()
 def write(self,p,d):
  p.parent.mkdir(parents=True,exist_ok=True);p.write_text(json.dumps(d))
 def fixture(self,name="v4",parent=None,now=.73,before=None,after=None):
  before=before or dict(dt_ps=.02,step_tolerance=2e-9,threads=4)
  after=after or dict(dt_ps=.02,step_tolerance=2e-7,threads=2)
  run=self.root/name;run.mkdir();anchor=self.root/(name+"_anchor");anchor.mkdir()
  pr=parent or self.root/(name+"_original");pr.mkdir(exist_ok=True)
  rt=anchor/"runtime_sources";rt.mkdir()
  pins={}
  for i in range(13):
   p=rt/("pin"+str(i));p.write_text("source"+str(i));pins[p.name]=m.digest(p)
  self.write(rt/"archive_manifest.json",dict(runtime_pins_sha256=pins))
  cons=dict(time_ps=now,source=dict(source_must_not_be_applied_again=True),controller={"state":{"status":"ACTIVE"}})
  meta=dict(cons,parameter_identity=before,source_pins_sha256=pins)
  np.savez_compressed(anchor/"checkpoint_original.npz",metadata_json=np.array(json.dumps(meta)))
  pm=dict(parameter_identity=before,source_pins_sha256=pins)
  self.write(anchor/"manifest_original.json",pm)
  state={k:k+"_hash" for k in m.STATE}
  cap=dict(checkpoint_bytes_sha256=m.digest(anchor/"checkpoint_original.npz"),
   manifest_sha256=m.digest(anchor/"manifest_original.json"),source_pins_sha256=pins,
   parameter_identity=before,conserved_metadata=cons,state_arrays_sha256=state)
  self.write(anchor/"capture_receipt.json",cap)
  seg=[dict(start_time_ps=0.,end_time_ps=now,origin="original",**before),
   dict(start_time_ps=now,end_time_ps=None,origin=name,**after)]
  argv=[]
  for k,v in after.items():argv.extend(["--"+k.replace("_","-"),str(v)])
  r=dict(schema="pysnspd.explicit_checkpoint_transfer.v2",parent=dict(output_dir=str(pr),anchor=str(anchor),
   checkpoint_bytes_sha256=cap["checkpoint_bytes_sha256"],manifest_sha256=cap["manifest_sha256"],
   runtime_archive_manifest_sha256=m.digest(rt/"archive_manifest.json")),output_dir=str(run),
   pins_before=pins,pins_after=pins,state_arrays_sha256_before=state,state_arrays_sha256_after=state,
   conserved_metadata_before=cons,conserved_metadata_after=cons,numerical_segments=seg,resume_argv=argv)
  path=run/"provenance/transfer_receipt.json";self.write(path,r)
  return run,anchor,path,r
 def test_base_exact(self):
  run,_,_,_=self.fixture();s=m.supplement(run)
  self.assertTrue(s["numerical_lineage_complete"]);self.assertEqual(len(s["numerical_segments"]),2)
 def test_dropped_manifest_pointer_chain_recovers_strict_prefix(self):
  run,_,_,_=self.fixture()
  new,_,_,_=self.fixture("v5",parent=run,now=5.42,
    before=dict(dt_ps=.02,step_tolerance=2e-7,threads=2),
    after=dict(dt_ps=.02,step_tolerance=2e-7,threads=2,local_energy_tolerance=1e-5))
  s=m.supplement(new);self.assertTrue(s["numerical_lineage_complete"])
  self.assertEqual([x["start_time_ps"] for x in s["numerical_segments"]],[0.,.73,5.42])
  self.assertEqual(s["numerical_segments"][0]["step_tolerance"],2e-9)
 def test_receipts_remain_identical(self):
  run,anchor,path,r=self.fixture();before=m.digest(path);m.supplement(run)
  self.assertEqual(m.digest(path),before)
 def test_bad_checkpoint_hash(self):
  run,anchor,_,_=self.fixture();(anchor/"checkpoint_original.npz").write_bytes(b"bad")
  s=m.supplement(run);self.assertFalse(s["numerical_lineage_complete"]);self.assertIsNone(s["numerical_segments"])
 def test_bad_frozen_source_hash(self):
  run,anchor,_,_=self.fixture();(anchor/"runtime_sources/pin0").write_text("bad")
  self.assertFalse(m.supplement(run)["numerical_lineage_complete"])
 def test_state_change_rejected(self):
  run,_,p,r=self.fixture();r["state_arrays_sha256_after"]=dict(r["state_arrays_sha256_after"],p="bad");self.write(p,r)
  self.assertFalse(m.supplement(run)["numerical_lineage_complete"])
 def test_clock_gap_rejected(self):
  run,_,p,r=self.fixture();r["numerical_segments"][0]["end_time_ps"]=.5;self.write(p,r)
  self.assertFalse(m.supplement(run)["numerical_lineage_complete"])
 def test_resume_identity_mismatch_rejected(self):
  run,_,p,r=self.fixture();r["resume_argv"][1]=".1";self.write(p,r)
  self.assertFalse(m.supplement(run)["numerical_lineage_complete"])
 def test_missing_declared_receipt_rejected(self):
  parent=self.root/"parent";parent.mkdir()
  with self.assertRaisesRegex(ValueError,"missing"):m.prior_receipt(parent,dict(explicit_transfer_receipt="provenance/missing.json"))
 def test_two_copies_must_agree(self):
  run,anchor,p,r=self.fixture();saved=anchor/"provenance/parent_transfer_receipt.json";self.write(saved,{"bad":True})
  with self.assertRaisesRegex(ValueError,"disagree"):m.prior_receipt(run,{},anchor)
 def test_ancestor_settings_mismatch_rejected(self):
  run,_,_,_=self.fixture()
  new,_,_,_=self.fixture("v5",parent=run,now=5.42,
    before=dict(dt_ps=.02,step_tolerance=2e-7,threads=4),
    after=dict(dt_ps=.02,step_tolerance=2e-7,threads=2,local_energy_tolerance=1e-5))
  self.assertFalse(m.supplement(new)["numerical_lineage_complete"])
if __name__=="__main__":unittest.main()
