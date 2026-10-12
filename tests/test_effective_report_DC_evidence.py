"""Optional SHA-pinned current-code basal evidence; no PDF or evolution."""
import copy,importlib.util,json,tempfile,unittest
from pathlib import Path
from unittest.mock import patch
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location("effective_report_dc",ROOT/"sandbox/stage5_effective_validation/build_final_report.py")
report=importlib.util.module_from_spec(spec);spec.loader.exec_module(report)
from sandbox.stage5_effective_validation import audit_runs as audit

class DCEvidenceTests(unittest.TestCase):
 def setUp(self):
  scratch=Path("/home/jdiaz/scratch/stage5_external_audit_development_20261011");scratch.mkdir(exist_ok=True)
  self.tmp=tempfile.TemporaryDirectory(dir=scratch);self.addCleanup(self.tmp.cleanup)
  self.base=Path(self.tmp.name);self.root=self.base/"run";self.root.mkdir();(self.root/"dc").mkdir()
  self.patch=patch.object(report,"SCRATCH",self.base);self.patch.start();self.addCleanup(self.patch.stop)
  self.pins={str(i):str(i)*64 for i in range(13)}
  self.phys={"exterior_inductance":{"total_H":103.3055069878189e-9},"current_A":21.297206e-6}
  self.reference={"manifest":{"source_pins_sha256":self.pins,"physicalparams":copy.deepcopy(self.phys)}}
  self.run={"manifest":{"run_complete":True,"source_pins_sha256":copy.deepcopy(self.pins),"physicalparams":copy.deepcopy(self.phys)},
   "summary":{"passed":True,"declared_mode":"DC","status":"DC_HOLD_COMPLETED",
    "final_time_ps":10.,"checks":{"DC_baseline_within_nV":True},"energy":{"maximum_absolute_defect_J":2e-41}},
   "history":{"time_s":np.array([0.,5e-12,10e-12]),"Vout_V":np.zeros(3),"Vdev_V":np.ones(3)*1e-17,"Is_A":np.ones(3)*2e-5},
   "snap":{"delta_J":np.ones((2,2),complex)*1e-22,"Te_K":np.ones((2,2))*.9,"Tph_K":np.ones((2,2))*.9}}
  self.files=("manifest.json","trajectory.npz","snapshots.npz","checkpoint_latest.npz","dc/effective_dc.npz")
  for name in self.files:(self.root/name).write_bytes(name.encode())
  self.ap=self.base/"audit.json"
  self.audit={"schema":"pysnspd.effective_saved_data_audit.v2","all_requested_checks_passed":True,"errors":{},
   "auditor_source_sha256":report.sha(audit.__file__),
   "runs":{"dc":{"path":str(self.root),**copy.deepcopy(self.run["summary"])}}}
  self.receipt={"schema":"pysnspd.current_effective_DC_evidence.v1",
   "status":"DC_HOLD_COMPLETED_AND_INDEPENDENTLY_AUDITED","no_photon":True,"stage5_admission_inferred":False,
   "path":str(self.root),"global_reference_energy_J":2e-19,"source_pins_sha256":copy.deepcopy(self.pins),
   "input_sha256":{str(self.root/name):report.sha(self.root/name) for name in self.files},
   "accepted_steps":2,"max_Vout_V":0.,"max_Vdev_V":1e-17,"max_energy_defect_J":2e-41,
   "max_gap_drift_J":0.,"max_Te_drift_K":0.,"max_Tph_drift_K":0.,"max_Is_drift_A":0.}
  self.rp=self.base/"receipt.json";self.write()
 def write(self):
  self.ap.write_text(json.dumps(self.audit));self.receipt["audit"]={"path":str(self.ap),"sha256":report.sha(self.ap)}
  self.rp.write_text(json.dumps(self.receipt));self.entry={"path":str(self.rp),"sha256":report.sha(self.rp)}
 def load(self,loader=None):
  self.write();self.inputs={}
  return report.load_dc_validation(self.entry,self.reference,input_hashes=self.inputs,
   load_run=loader or (lambda *a,**kw:self.run))
 def test_current_completed_DC_has_same_pins_and_parameters(self):
  result,record=self.load()
  self.assertEqual(result["max_Vout_V"],0.);self.assertEqual(result["exterior_inductance_difference_fraction"],0.)
  self.assertEqual(result["duration_ps"],10.);self.assertFalse(result["stage5_admission_inferred"])
  self.assertEqual(record["sha256"],report.sha(self.rp));self.assertEqual(len(self.inputs),7)
 def test_required_receipt_SHA_and_each_saved_data_SHA(self):
  with self.assertRaises(ValueError):report.load_dc_validation({"path":str(self.rp)},self.reference,input_hashes={})
  for name in self.files:
   with self.subTest(name=name):
    value=self.receipt["input_sha256"].pop(str(self.root/name))
    with self.assertRaises(ValueError):self.load()
    self.receipt["input_sha256"][str(self.root/name)]=value
 def test_old_pins_do_not_become_current_through_a_label(self):
  self.run["manifest"]["source_pins_sha256"]["0"]="f"*64
  with self.assertRaisesRegex(ValueError,"full13"):self.load()
 def test_full_circuit_and_material_are_protected(self):
  self.run["manifest"]["physicalparams"]["exterior_inductance"]["total_H"]*=1.0001
  with self.assertRaisesRegex(ValueError,"full circuit"):self.load()
 def test_failed_or_different_auditor_is_rejected(self):
  self.audit["all_requested_checks_passed"]=False
  with self.assertRaises(ValueError):self.load()
  self.audit["all_requested_checks_passed"]=True;self.audit["auditor_source_sha256"]="0"*64
  with self.assertRaises(ValueError):self.load()
 def test_incomplete_or_photon_run_cannot_support_basal_evidence(self):
  for key,value in (("passed",False),("declared_mode","PHOTON"),("status","DEVELOPMENT_DURATION_REACHED")):
   old=self.run["summary"][key];self.run["summary"][key]=value
   with self.assertRaises(ValueError):self.load()
   self.run["summary"][key]=old
 def test_edited_data_after_loading_is_rejected(self):
  def loader(*a,**kw):
   (self.root/"trajectory.npz").write_bytes(b"changed");return self.run
  with self.assertRaisesRegex(ValueError,"checksum changed"):self.load(loader)
 def test_fabricated_zero_drift_is_rejected(self):
  self.run["snap"]["Te_K"][1,0]=1.
  with self.assertRaisesRegex(ValueError,"metrics differ"):self.load()
 def test_photon_or_global_admission_claim_is_rejected(self):
  self.receipt["stage5_admission_inferred"]=True
  with self.assertRaises(ValueError):self.load()
 def test_note_uses_nV_and_does_not_add_pages(self):
  result,_=self.load();text=report.dc_validation_latex(result)
  self.assertIn("nV",text);self.assertIn("sin fotón",text);self.assertNotIn("newpage",text)
  self.assertIn("nulas en los datos guardados",text);self.assertEqual(report.dc_validation_latex(None),"")
 def test_stage5_false_blocks_before_optional_DC_loading(self):
  p=self.base/"stage5false.json";p.write_text(json.dumps({"stage5":False,"dc_validation":self.entry}))
  target=self.base/"forbidden"
  with patch.object(report,"load_dc_validation",side_effect=AssertionError("must not load")):
   with self.assertRaisesRegex(ValueError,"stage5: true"):
    report.main(["--admission",str(p),"--build-root",str(target),"--compile","--cpu","4"])
  self.assertFalse(target.exists())

if __name__=="__main__":unittest.main()
