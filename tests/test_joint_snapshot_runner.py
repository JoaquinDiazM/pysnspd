"""Small end-to-end checkpoint/runner check; no cached detector matrix run."""
import argparse,json,tempfile,unittest
from pathlib import Path
from unittest.mock import patch
from test_joint_evolution import physical_small_system
from sandbox.stage5_photon import run_joint_snapshot as runner


class JointSnapshotRunnerTests(unittest.TestCase):
    def factory(self,*args,**kwargs):
        system=physical_small_system()
        system.provenance={"pinned_files_sha256":{},"test_fixture":True}
        system.initial_data={"coordinates_m":system.graph.coordinates_bar*1e-8}
        return system

    def args(self,path,**more):
        values=dict(initial=Path("/home/jdiaz/scratch/test_initial_fixture"),
            output=path,workers=1,pilot=False,execute=True,resume=False)
        values.update(more)
        return argparse.Namespace(**values)

    def test_whole_pipeline_solves_and_saves_without_claiming_stage5(self):
        with tempfile.TemporaryDirectory(prefix="stage51_runner_test_",
                dir="/home/jdiaz/scratch/tmp/pysnspd") as temp:
            output=Path(temp)/"small"
            with patch.object(runner,"build_initial_system",self.factory):
                runner.run(self.args(output))
            summary=json.loads((output/"summary.json").read_text())
            self.assertEqual(len(list((output/"columns").glob("*.npy"))),8)
            self.assertEqual(summary["physical_time_steps"],0)
            self.assertFalse(summary["stage51_complete"])
            self.assertFalse(summary["stage5_complete"])
            self.assertTrue((output/"joint_candidate_fields.npz").exists())
            with patch.object(runner,"build_initial_system",self.factory):
                with self.assertRaisesRegex(ValueError,"Completed result"):
                    runner.run(self.args(output,resume=True))
            with patch.object(runner,"build_initial_system",self.changed_factory):
                with self.assertRaisesRegex(ValueError,"identity changed: input_pins"):
                    runner.run(self.args(output,resume=True))

    def changed_factory(self,*args,**kwargs):
        system=self.factory()
        system.provenance["pinned_files_sha256"]={"reference_summary.json":"modified"}
        return system

    def test_pilot_never_starts_column_assembly(self):
        with tempfile.TemporaryDirectory(prefix="stage51_runner_test_",
                dir="/home/jdiaz/scratch/tmp/pysnspd") as temp:
            output=Path(temp)/"pilot"
            with patch.object(runner,"build_initial_system",self.factory):
                runner.run(self.args(output,pilot=True,execute=False))
            summary=json.loads((output/"pilot.json").read_text())
            self.assertEqual(summary["unknowns"],8)
            self.assertLess(summary["gauge_shift_scaled_row_difference"],1e-10)
            self.assertFalse((output/"columns").exists())
