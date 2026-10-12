"""Lightweight ETA/session regressions with a control-flow mock; no physical run."""
import ast,contextlib,importlib.util,io,json,os,unittest
from pathlib import Path
from unittest.mock import patch
ROOT=Path(os.environ.get("PYSNSPD_TEST_REPO",Path(__file__).resolve().parents[1]))
CANDIDATES=os.environ.get("PYSNSPD_FRONTEND_CANDIDATES")
SOURCE=(Path(CANDIDATES)/"05_effective_photon_run_candidate.py" if CANDIDATES
 else ROOT/"pipelines/05_effective_photon_run.py")
tree=ast.parse(SOURCE.read_text())
functions=[node for node in tree.body if isinstance(node,ast.FunctionDef)
 and node.name in ("_session_eta_s","_session_wall_info")]
namespace={}
exec(compile(ast.Module(body=functions,type_ignores=[]),str(SOURCE),"exec"),namespace)

class SessionEtaUnitTests(unittest.TestCase):
 def test_resume_uses_new_advance_not_absolute_photon_clock(self):
  elapsed=69.;start=12.42;now=12.54;remaining=200.-now
  eta=namespace["_session_eta_s"](elapsed,start,now,remaining)
  self.assertAlmostEqual(eta,107789.5,places=6)
  self.assertAlmostEqual(eta/(elapsed/now*remaining),now/(now-start),places=10)
 def test_no_new_advance_has_no_eta(self):
  self.assertIsNone(namespace["_session_eta_s"](69.,12.42,12.42,187.58))
 def test_fresh_session_retains_expected_estimate(self):
  self.assertEqual(namespace["_session_eta_s"](60.,0.,2.,8.),240.)
  self.assertEqual(namespace["_session_eta_s"](60.,0.,2.,0.),0.)
 def test_elapsed_alias_is_explicitly_one_session(self):
  data=namespace["_session_wall_info"](69.,12.42)
  self.assertEqual(data["elapsed_wall_s"],data["elapsed_session_wall_s"])
  self.assertEqual(data["session_start_time_ps"],12.42)
  self.assertIn("current_process_session",data["elapsed_wall_scope"])

spec=importlib.util.spec_from_file_location("_session_eta_transactions",
 ROOT/"tests/test_effective_frontend_transactions.py")
frontend=importlib.util.module_from_spec(spec);spec.loader.exec_module(frontend)

class ResumeEtaIntegrationTests(frontend.FrontendTransactionTests):
 def test_mock_resume_progress_uses_loaded_clock_and_preserves_alias(self):
  out=self.work/"eta-resume";extra=["--threshold-mV","10"]
  stopped={"ready":False};original=self.r.append_history
  def append(history,row):
   original(history,row)
   if row["time_s"]>=.05e-12:stopped["ready"]=True
  with self.resources(),patch.object(self.r,"append_history",side_effect=append):
   with self.assertRaises(KeyboardInterrupt):
    self.r.run(self.r.parser().parse_args(self.argv(out,extra)),stop_requested=lambda:stopped["ready"])
  _,meta=self.assert_coherent(out);start=meta["time_ps"]
  count={"n":0}
  def clock():
   count["n"]+=1
   return max(0.,(count["n"]-2)*69.)
  output=io.StringIO()
  with self.resources(),patch.object(self.r.time,"monotonic",side_effect=clock),contextlib.redirect_stdout(output):
   result=self.r.run(self.r.parser().parse_args(self.argv(out,extra+["--resume"])))
  records=[json.loads(line) for line in output.getvalue().splitlines() if line.startswith("{")]
  progress=[row for row in records if row.get("status")=="RUNNING"][0]
  self.assertEqual(progress["session_start_time_ps"],start)
  self.assertEqual(progress["elapsed_session_wall_s"],progress["elapsed_wall_s"])
  self.assertAlmostEqual(progress["eta_to_current_target_s"],
   progress["elapsed_session_wall_s"]/(progress["time_ps"]-start)*
   progress["remaining_to_current_observation_target_ps"])
  self.assertEqual(result["session_timing"]["start_time_ps"],start)
  self.assertEqual(result["execution"]["elapsed_session_wall_s"],result["execution"]["elapsed_wall_s"])
  _,final=self.assert_coherent(out)
  self.assertEqual(final["session_start_time_ps"],start)
  self.assertNotIn("accumulated_wall_s",final)

if __name__=="__main__":unittest.main()
