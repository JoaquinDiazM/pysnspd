"""Fast campaign safety/identity and censored-observable regressions; no runs."""
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
def module(name,path):
    spec = importlib.util.spec_from_file_location(name,ROOT/path)
    value = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(value)
    return value
campaign = module("sensitivity_campaign_test","pipelines/06_effective_sensitivity_campaign.py")
visual = module("sensitivity_visual_test","plot_pipelines/06_plot_effective_sensitivity_campaign.py")
runner = module("sensitivity_numerical_runner_test","pipelines/05_effective_photon_run.py")


class SensitivityCampaignTests(unittest.TestCase):
    def test_all_generated_commands_are_accepted_by_the_actual_runner(self):
        cases = campaign.case_catalog()
        self.assertEqual(len(cases),16)
        self.assertEqual(len({case["name"] for case in cases}),16)
        for case in cases:
            with self.subTest(case=case["name"]):
                command = campaign.runner_command(case,Path("/home/jdiaz/scratch/campaign_test"))
                args = runner.parser().parse_args(command[3:])
                runner.validate_args(args)
                self.assertEqual(runner.parameter_identity(args),case["parameters"])
                self.assertEqual((args.length_nm,args.local_energy_tolerance),(320.,2.5e-6))
                self.assertEqual((args.threshold_mV,args.confirmation_ps,args.extra_time_ps,args.ceiling_ps),
                    (.1,.5,5.,200.))
                self.assertFalse(args.dc_only)
                self.assertIsNone(args.duration_ps)
                self.assertLessEqual(len(case["changed_input"]),1)

    def test_groups_always_include_their_wavelength_matched_references(self):
        self.assertEqual(len(campaign.selected_cases("bias")),5)
        self.assertEqual(len(campaign.selected_cases("rates")),7)
        self.assertEqual(len(campaign.selected_cases("source")),7)
        for group in ("bias","rates","source","all"):
            cases = campaign.selected_cases(group)
            names = {case["name"] for case in cases}
            self.assertTrue(all(case["reference"] in names for case in cases))

    def terminal_fixture(self,directory):
        case = campaign.case_catalog()[0]
        pins = {"physical.py":"unchanged"}
        hashes = {}
        for name in ("trajectory.npz","snapshots.npz","checkpoint_latest.npz",
                "dc/effective_dc.npz","material.npz","material.json"):
            path = directory/name
            path.parent.mkdir(parents=True,exist_ok=True)
            path.write_bytes(name.encode())
            hashes[name] = campaign.sha256(path)
        manifest = dict(parameter_identity=case["parameters"],source_pins_sha256=pins,
            run_complete=True,status="TRIGGER_MARGIN_COMPLETE",output_sha256=hashes)
        (directory/"manifest.json").write_text(json.dumps(manifest))
        return case,pins,manifest

    def test_completed_skip_requires_exact_parameters_sources_and_terminal_bytes(self):
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            case,pins,manifest = self.terminal_fixture(directory)
            campaign.check_case_identity(case,directory,pins,require_complete=True)
            (directory/"trajectory.npz").write_bytes(b"corruption")
            with self.assertRaisesRegex(ValueError,"Terminal output changed"):
                campaign.check_case_identity(case,directory,pins,require_complete=True)
            manifest["parameter_identity"]["tau_ee_Tc_ps"] = 6.
            (directory/"manifest.json").write_text(json.dumps(manifest))
            with self.assertRaisesRegex(ValueError,"parameter identity"):
                campaign.check_case_identity(campaign.case_catalog()[0],directory,pins,require_complete=True)

    def test_incomplete_and_failed_runs_are_not_successful_terminals(self):
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            case,pins,manifest = self.terminal_fixture(directory)
            manifest["run_complete"] = False
            (directory/"manifest.json").write_text(json.dumps(manifest))
            campaign.check_case_identity(case,directory,pins,require_complete=False)
            with self.assertRaisesRegex(ValueError,"Successful finite-window"):
                campaign.check_case_identity(case,directory,pins,require_complete=True)

    def test_censoring_never_becomes_a_latency_or_slope(self):
        manifest = dict(status="TIME_CENSORED_NO_TRIGGER",observation=dict(
            status="TIME_CENSORED_NO_TRIGGER",trigger_confirmed=False,censored=True,
            crossing_time_s=None,latency_s=None,no_trigger_ceiling_s=200e-12))
        result = visual.event_metrics(manifest,np.array([0.,100.,200.]),np.array([0.,.0002,.0001]))
        self.assertTrue(result["censored"])
        self.assertEqual(result["censoring_ps"],200.)
        self.assertIsNone(result["relative_latency_ps"])
        self.assertIsNone(result["local_slope_mV_per_ps"])

    def test_confirmed_slope_uses_actual_saved_event_and_local_window(self):
        time = np.linspace(0,10,101)
        voltage = .02*time
        manifest = dict(status="TRIGGER_MARGIN_COMPLETE",observation=dict(
            status="TRIGGER_MARGIN_COMPLETE",trigger_confirmed=True,censored=False,
            crossing_time_s=5e-12,latency_s=5e-12))
        result = visual.event_metrics(manifest,time,voltage)
        self.assertAlmostEqual(result["relative_latency_ps"],5.)
        self.assertAlmostEqual(result["local_slope_mV_per_ps"],.02)
        self.assertLessEqual(result["slope_fit_window_ps"][1]-result["slope_fit_window_ps"][0],.5)

    def test_common_comparison_excludes_unshared_tail_and_censored_latency(self):
        reference = dict(time_ps=np.array([0.,5.,10.]),voltage_mV=np.array([0.,.1,9.]),
            metrics=dict(confirmed=True,relative_latency_ps=5.))
        candidate = dict(time_ps=np.array([0.,2.5,5.]),voltage_mV=np.array([0.,.05,.1]),
            metrics=dict(confirmed=False,relative_latency_ps=None))
        result = visual.common_metrics(reference,candidate)
        self.assertEqual(result["common_window_ps"],[0.,5.])
        self.assertAlmostEqual(result["maximum_absolute_Vout_difference_mV"],0.)
        self.assertAlmostEqual(result["reference_observed_common_max_mV"],.1)
        self.assertIsNone(result["delta_relative_latency_ps"])

    def test_execution_coordinator_writes_receipts_and_skips_exact_completed_case(self):
        """Exercise the coordinator with file fixtures, never a numerical child."""
        with tempfile.TemporaryDirectory() as temporary:
            scratch = Path(temporary)
            directory = scratch/"campaign"
            case = campaign.case_catalog()[0]
            calls = []
            def prepared_terminal(command,logfile,env,**kwargs):
                calls.append(command)
                output = Path(command[command.index("--output-dir")+1])
                output.mkdir()
                fixture_case,pins,manifest = self.terminal_fixture(output)
                manifest["source_pins_sha256"] = campaign.source_pins()
                manifest["observation"] = {"status":"TRIGGER_MARGIN_COMPLETE"}
                (output/"manifest.json").write_text(json.dumps(manifest))
            def audited_fixture(current,directory,env):
                path = directory/"audits"/"fixture_audit.json"
                path.write_text(json.dumps({"all_requested_checks_passed":True}))
                return dict(path=str(path),sha256=campaign.sha256(path),passed=True,
                    manifest_sha256=campaign.sha256(directory/"runs"/current["name"]/"manifest.json"))
            with patch.object(campaign,"SCRATCH",scratch),patch.object(campaign,"case_catalog",return_value=[case]), \
                    patch.object(campaign,"selected_cases",return_value=[case]), \
                    patch.object(campaign,"run_logged",side_effect=prepared_terminal), \
                    patch.object(campaign,"audit_case",side_effect=audited_fixture), \
                    patch.object(campaign,"child_environment",return_value={}):
                result = campaign.main(["--output-dir",str(directory),"--execute","--group","bias"])
                self.assertEqual(result["status"],"SELECTED_GROUP_COMPLETE_AUDITED")
                self.assertEqual(len(calls),1)
                result = campaign.main(["--output-dir",str(directory),"--execute","--resume","--group","bias"])
                self.assertEqual(len(calls),1)
                self.assertTrue(result["selected_group_complete"])
                self.assertTrue((directory/"aggregate_bias.json").is_file())
                self.assertTrue((directory/"campaign_execution.json").is_file())

    def test_memory_bound_uses_expanded_arrays_not_compressed_file_size(self):
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary)/"compressed.npz"
            data = {key:np.zeros(100000) for key in ("time_s","Vout_V","Vdev_V","Is_A")}
            np.savez_compressed(path,**data)
            expanded = visual.uncompressed_trace_bytes(path)
            self.assertGreaterEqual(expanded,sum(value.nbytes for value in data.values()))
            self.assertGreater(expanded,100*path.stat().st_size)

    def test_campaign_lock_refuses_a_second_writer_and_releases_on_exit(self):
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            with campaign.execution_lock(directory):
                with self.assertRaisesRegex(RuntimeError,"Another invocation"):
                    with campaign.execution_lock(directory):
                        self.fail("Concurrent writer admitted")
            with campaign.execution_lock(directory):
                pass

    def test_plan_only_does_not_invoke_a_numerical_process(self):
        with tempfile.TemporaryDirectory() as temporary:
            scratch = Path(temporary)
            with patch.object(campaign,"SCRATCH",scratch),patch.object(campaign,"run_logged") as process:
                plan = campaign.main(["--output-dir",str(scratch/"fresh"),"--group","bias","--plan-only"])
                self.assertEqual(len(plan["cases"]),16)
                process.assert_not_called()
                self.assertFalse((scratch/"fresh"/"runs").exists())

if __name__=="__main__":
    unittest.main()
