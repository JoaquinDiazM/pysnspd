"""Energy acceptance is a frontend transaction; no physical trajectory runs."""
import ast,importlib.util,json,math,os,sys,unittest
from pathlib import Path
from unittest.mock import patch
import numpy as np
ROOT=Path(os.environ.get("PYSNSPD_TEST_REPO",Path(__file__).resolve().parents[1]))
CANDIDATES=os.environ.get("PYSNSPD_FRONTEND_CANDIDATES")
SOURCE=Path(CANDIDATES)/"05_effective_photon_run_candidate.py" if CANDIDATES else ROOT/"pipelines/05_effective_photon_run.py"
tree=ast.parse(SOURCE.read_text())
names=("_local_energy_step","_energy_guard_from_history","_energy_guard_ledger")
functions=[node for node in tree.body if isinstance(node,ast.FunctionDef) and node.name in names]
namespace=dict(np=np,math=math)
exec(compile(ast.Module(body=functions,type_ignores=[]),str(SOURCE),"exec"),namespace)

class LocalEnergyUnitTests(unittest.TestCase):
 def test_norm_counts_absolute_work_and_uses_SI_without_a_floor(self):
  reference=1e-19;work=3e-20;dw=-2e-20;before=5e-19;after=before+dw+1e-25
  result=namespace["_local_energy_step"](before,after,dw,reference,work)
  self.assertAlmostEqual(result["local_energy_denominator_J"]/1.5e-19,1.)
  self.assertAlmostEqual(result["local_energy_defect_J"]/1e-25,1.,delta=1e-9)
  self.assertAlmostEqual(result["local_energy_fraction"],1e-25/1.5e-19,delta=1e-15)
 def test_bad_acceptance_inputs_are_rejected(self):
  for values in ((0.,1.,0.,0.,0.),(0.,1.,0.,1.,-1.),(0.,math.nan,0.,1.,0.)):
   with self.subTest(values=values),self.assertRaises((ValueError,ArithmeticError)):
    namespace["_local_energy_step"](*values)
 def test_legacy_sums_are_measured_from_every_accepted_increment(self):
  history=dict(energy_total_J=[10.,12.01,11.02],external_energy_integral_J=[0.,2.,1.])
  measured=namespace["_energy_guard_from_history"](history,1.)
  self.assertEqual(measured["absolute_external_energy_integral_J"],3.)
  self.assertAlmostEqual(measured["sum_abs_local_energy_defect_J"],.02)
  self.assertAlmostEqual(measured["max_accepted_local_energy_fraction"],.01/3.)
  legacy=namespace["_energy_guard_ledger"]({},history,1.,12.4)
  self.assertEqual(legacy["guard_active_since_time_ps"],12.4)
  self.assertEqual(legacy["max_guarded_accepted_local_energy_fraction"],0.)
  self.assertEqual(legacy["local_energy_rejections"],0)
  self.assertIn("pre_guard_history",legacy["energy_guard_ledger_origin"])
 def test_saved_ledger_is_preserved_without_recomputing_old_history(self):
  history=dict(energy_total_J=[10.],external_energy_integral_J=[0.])
  ledger=namespace["_energy_guard_ledger"](None,history,1.,0.)
  ledger.update(absolute_external_energy_integral_J=3.1,sum_abs_local_energy_defect_J=.01,
    local_energy_rejections=2,max_accepted_local_energy_fraction=.004,
    max_guarded_accepted_local_energy_fraction=.004)
  self.assertEqual(namespace["_energy_guard_ledger"](ledger,history,1.,12.),ledger)
  with self.assertRaisesRegex(ValueError,"reference changed"):
   namespace["_energy_guard_ledger"](ledger,history,1.1,12.)
  with self.assertRaisesRegex(ValueError,"Partial"):
   namespace["_energy_guard_ledger"]({"local_energy_rejections":2},history,1.,12.)

spec=importlib.util.spec_from_file_location("_energy_guard_transactions",ROOT/"tests/test_effective_frontend_transactions.py")
frontend=importlib.util.module_from_spec(spec);spec.loader.exec_module(frontend)
GUARD_KEYS=("absolute_external_energy_integral_J","sum_abs_local_energy_defect_J",
 "local_energy_rejections","max_accepted_local_energy_fraction",
 "max_guarded_accepted_local_energy_fraction","guard_active_since_time_ps",
 "local_energy_reference_energy_J","energy_guard_ledger_origin")

class EnergyGuardTransactionTests(frontend.FrontendTransactionTests):
 def assert_same_as_reference(self,out,reference):
  super().assert_same_as_reference(out,reference)
  left=json.loads(self.checkpoint(out)["metadata_json"].item())
  right=json.loads(self.checkpoint(reference)["metadata_json"].item())
  for key in GUARD_KEYS:self.assertEqual(left[key],right[key],key)
 def test_local_energy_rejects_before_state_clock_trigger_and_history_commit(self):
  out=self.work/"local-rejection";attempts=[];original=self.system.midpoint_step
  def biased(state,h_ps,**kw):
   attempts.append((h_ps,state.delta_bar.copy(),state.circuit_SI.copy()))
   candidate,receipt=original(state,h_ps,**kw)
   candidate.delta_bar+=10*h_ps**3
   if h_ps>.005*(1+1e-12):candidate.circuit_SI[1]=-1.
   return candidate,receipt
  with patch.object(self.system,"midpoint_step",side_effect=biased):
   result=self.execute(out,["--threshold-mV","10"])
  saved,meta=self.assert_coherent(out)
  self.assertEqual(attempts[0][0],.01);self.assertEqual(attempts[1][0],.005)
  np.testing.assert_array_equal(attempts[0][1],attempts[1][1])
  np.testing.assert_array_equal(attempts[0][2],attempts[1][2])
  self.assertEqual(float(saved["trajectory__time_s"][1]),.005e-12)
  self.assertEqual(meta["local_energy_rejections"],meta["total_rejections"])
  self.assertGreater(meta["local_energy_rejections"],0)
  self.assertLessEqual(meta["max_guarded_accepted_local_energy_fraction"],1e-5)
  self.assertTrue(np.all(saved["trajectory__Is_A"]>0.))
  self.assertFalse(result["observation"]["trigger_confirmed"])
  self.assertEqual(result["local_energy_control"]["local_energy_rejections"],meta["local_energy_rejections"])
  self.assertEqual(set(saved)&{"trajectory__"+key for key in self.r.HISTORY_KEYS},
    {"trajectory__"+key for key in self.r.HISTORY_KEYS})
  self.assertEqual(sum(key.startswith("trajectory__") for key in saved),len(self.r.HISTORY_KEYS))
 def test_endpoint_rhs_is_reused_once_per_trial(self):
  out=self.work/"endpoint-cached";count=0;original=self.system.evaluate
  def evaluate(state):
   nonlocal count
   count+=1
   return original(state)
  with patch.object(self.system,"evaluate",side_effect=evaluate):
   result=self.execute(out,["--threshold-mV","10"])
  self.assertEqual(count,2+result["execution"]["accepted_steps"])
 def test_partial_history_rollback_restores_all_energy_ledger_fields(self):
  out=self.work/"guard-rollback";original=self.r.append_history;calls=0
  def append(history,row):
   nonlocal calls
   calls+=1
   if calls==2:
    for key in self.r.HISTORY_KEYS[:4]:history[key].append(float(row[key]))
    raise KeyboardInterrupt("guard transaction rollback")
   original(history,row)
  with patch.object(self.r,"append_history",side_effect=append):
   with self.assertRaises(KeyboardInterrupt):self.execute(out)
  _,meta=self.assert_coherent(out)
  self.assertEqual(meta["absolute_external_energy_integral_J"],0.)
  self.assertEqual(meta["sum_abs_local_energy_defect_J"],0.)
  self.assertEqual(meta["max_accepted_local_energy_fraction"],0.)
  self.assertEqual(meta["local_energy_rejections"],0)
 def test_numerical_rejection_is_not_mislabeled_as_energy_rejection(self):
  self.max_step=.005
  out=self.work/"solver-rejection";self.execute(out,["--threshold-mV","10"])
  _,meta=self.assert_coherent(out)
  self.assertGreater(meta["total_rejections"],0);self.assertEqual(meta["local_energy_rejections"],0)
 def test_legacy_resume_reconstructs_measured_ledger_and_labels_activation(self):
  out=self.work/"legacy";self.signal_after_accepts=5
  extra=["--threshold-mV","10"]
  with self.assertRaises(KeyboardInterrupt):self.execute(out,extra,main=True)
  saved,oldmeta=self.assert_coherent(out);self.signal_after_accepts=None
  for key in GUARD_KEYS:oldmeta.pop(key)
  saved["metadata_json"]=np.array(json.dumps(oldmeta))
  self.r.atomic_npz(out/"checkpoint_latest.npz",saved)
  self.execute(out,extra+["--resume"])
  _,meta=self.assert_coherent(out)
  self.assertEqual(meta["guard_active_since_time_ps"],oldmeta["time_ps"])
  self.assertIn("pre_guard_history",meta["energy_guard_ledger_origin"])
  measured=self.r._energy_guard_from_history(
    {key:meta2 for key,meta2 in (
      ("energy_total_J",self.checkpoint(out)["trajectory__energy_total_J"]),
      ("external_energy_integral_J",self.checkpoint(out)["trajectory__external_energy_integral_J"]))},1.)
  self.assertAlmostEqual(meta["absolute_external_energy_integral_J"],measured["absolute_external_energy_integral_J"])
  self.assertAlmostEqual(meta["sum_abs_local_energy_defect_J"],measured["sum_abs_local_energy_defect_J"])
 def test_dc_reference_uses_system_volume_scale(self):
  out=self.work/"dc-reference"
  manifest=self.execute(out,["--dc-only"])
  self.assertEqual(manifest["numerical"]["local_energy_guard"]["reference_kind"],"sum N0 Eu^2 nodal volume")
  self.assertEqual(manifest["local_energy_control"]["local_energy_reference_energy_J"],1.)

if __name__=="__main__":unittest.main()
