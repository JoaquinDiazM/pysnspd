"""Receipt provenance only: tiny saved fixtures, no figures, PDF or evolution."""
import importlib.util
import inspect
import json
import os
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
SCRATCH = Path("/home/jdiaz/scratch")

def load(name, relative):
    spec = importlib.util.spec_from_file_location(name, ROOT/relative)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module

report = load("effective_report_dependency_receipt", "sandbox/stage5_effective_validation/build_final_report.py")
visual = load("effective_visual_dependency_receipt", "plot_pipelines/05_plot_effective_photon_run.py")

class VisualDependencyProvenanceTests(unittest.TestCase):
    def setUp(self):
        self.work = tempfile.TemporaryDirectory(dir=SCRATCH, prefix="visual_dependency_receipt_test_")
        self.addCleanup(self.work.cleanup)
        self.base = Path(self.work.name)
        self.helper = self.base/"complex_gap.py"
        self.helper.write_bytes(report.COMPLEX_GAP_SOURCE.read_bytes())

    def report_receipt(self):
        admission = self.base/"admission.json"
        admission.write_text(json.dumps({"stage5": True}))
        png = self.base/"existing_figure_fixture.png"
        png.write_bytes(b"Synthetic existing-figure fixture; not a rendered image.")
        thickness = self.base/"thickness.json"
        thickness.write_text(json.dumps({"outputs": {"png": {"path": str(png), "sha256": report.sha(png)}}}))
        failure = self.base/"failure.json"
        failure.write_text("{}")
        config = {"thickness_audit": str(thickness), "failure_audit": str(failure)}
        runs = {"775": {"summary": {}}}
        build = self.base/"report_build"
        replacements = {
            "admission_inputs": lambda _: ({"report": config}, runs, [], [], {}),
            "overview": lambda *a: "unused.png",
            "gap_maps": lambda *a: ("unused.png", [0.]),
            "energy_maps": lambda *a: "unused.png",
            "energy_balance": lambda *a: "unused.png",
            "convergence": lambda *a: "unused.png",
            "make_tex": lambda *a: "% Unit fixture: no figures or compiled document.\n",
            "expected_report_pages": lambda *a: 1,
        }
        with patch.object(report, "COMPLEX_GAP_SOURCE", self.helper), patch.multiple(report, **replacements):
            report.main(["--admission", str(admission), "--build-root", str(build),
                         "--cpu", str(min(os.sched_getaffinity(0)))])
        return admission, build, json.loads((build/"final_report_receipt.json").read_text())

    def test_recorded_helper_is_the_implementation_used_for_drawing(self):
        for module in (report, visual):
            self.assertEqual(module.COMPLEX_GAP_SOURCE,
                Path(inspect.getfile(module.draw_complex_gap)).resolve())

    def test_report_receipt_pins_shared_helper_without_rendering(self):
        _, _, receipt = self.report_receipt()
        self.assertEqual(receipt["input_sha256"][str(self.helper)], report.sha(self.helper))
        self.assertFalse(receipt["compiled"])
        self.assertEqual(receipt["simulations_launched"], 0)

    def test_publication_rejects_changed_shared_helper_before_copying(self):
        admission, build, receipt = self.report_receipt()
        # Exercise the existing publish verifier; this is deliberately no PDF.
        artifact = build/"not_a_pdf_fixture.bin"
        artifact.write_bytes(b"Synthetic artifact used only to reach the input checksum gate.")
        receipt.update(compiled=True, pdf_path=str(artifact), pdf_sha256=report.sha(artifact))
        receipt_path = build/"final_report_receipt.json"
        receipt_path.write_text(json.dumps(receipt))
        self.helper.write_bytes(self.helper.read_bytes()+b"\n# altered after receipt\n")
        target = self.base/"must_not_be_written.pdf"
        with self.assertRaisesRegex(ValueError, "Evidence changed after compilation"):
            report.publish_report(admission, build, target, self.base/"qa_never_read.json")
        self.assertFalse(target.exists())
        self.assertEqual(json.loads(receipt_path.read_text()), receipt)

    def test_visual_run_receipt_pins_shared_helper_without_rendering(self):
        run = self.base/"saved_run"; run.mkdir()
        manifest = {"source": {"retained_energy_J": 1e-19}, "observation": {}}
        (run/"manifest.json").write_text(json.dumps(manifest))
        np.savez(run/"trajectory.npz", time_s=np.array([0., 1e-12]))
        np.savez(run/"snapshots.npz", time_s=np.array([0., 1e-12]))
        history = {}
        maps = {"coordinates_nm": np.array([[0., 0.], [1., 0.], [0., 1.]]),
                "time_ps": np.array([0., 1.]), "selected_indices": np.array([0, 1])}
        from sandbox.stage4_core import parallel_runtime
        cpu = min(os.sched_getaffinity(0))
        budget = {"coordinator_cpu": cpu, "memory_limit_bytes": 2**30}
        with patch.object(visual, "COMPLEX_GAP_SOURCE", self.helper), \
             patch.object(visual, "prepare_view", return_value=(history, maps)), \
             patch.object(visual, "trace_figures") as trace, \
             patch.object(visual, "map_figure") as maps_draw, \
             patch.object(visual, "apply_thesis_style"), \
             patch.object(parallel_runtime, "linux_resources", return_value={}), \
             patch.object(parallel_runtime, "resource_budget", return_value=budget):
            receipt = visual.run(run, self.base/"visual_receipt", dpi=72)
        self.assertEqual(receipt["input_sha256"][str(self.helper)], visual.sha256(self.helper))
        self.assertEqual(len(receipt["input_sha256"]), 4)
        self.assertEqual(receipt["products"], [])
        self.assertFalse(receipt["model_admission_inferred"])
        trace.assert_called_once()
        self.assertEqual(maps_draw.call_count, 3)
        self.assertFalse(any((self.base/"visual_receipt").glob("*.png")))

if __name__ == "__main__":
    unittest.main()
