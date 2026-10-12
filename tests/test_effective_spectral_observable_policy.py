"""Synthetic regression of the analysis policy; no numerical evolution."""
import copy
import importlib.util
import json
import math
from pathlib import Path
import tempfile
import unittest
from unittest import mock

REPO = Path(__file__).resolve().parents[1]
BASE = REPO / "sandbox/stage5_effective_validation"
spec = importlib.util.spec_from_file_location("spectral_observable_policy", BASE / "spectral_observable_policy.py")
policy = importlib.util.module_from_spec(spec)
spec.loader.exec_module(policy)
# Immutable adopted decision and measured partial comparison copied as small test fixtures.
# Original audit SHA34052805b6af8a8ea024fd44266a82061a873b342a32b4b0dcd6635e4dd25b0c.
DECISION_PAYLOAD = '{\n  "schema": "pysnspd.prospective_spectral_observable_policy.v1",\n  "utc": "2026-10-11T17:21:31.679212+00:00",\n  "status": "ADOPTED_PRACTICAL_SPECTRAL_OBSERVABLE_POLICY_NOT_MODEL_ADMISSION",\n  "review_path": "/home/jdiaz/scratch/stage5_spectral_mesh_scope_review_20261011/proposal.json",\n  "review_sha256": "f251ac1aaa56e821b0ffb23116b40817a6b720e4dc64c5b7b4a64c2ae72bc02b",\n  "original_failures_preserved": {\n    "spectral_nodal_minimum_over_DC": 0.02276239242656607,\n    "spectral_limit": 0.02,\n    "mesh_deficit_over_peak": 0.06583534006348958,\n    "mesh_limit": 0.05\n  },\n  "spectral_limits": {\n    "nodal_minimum_over_DC": 0.03,\n    "central_field_L2_over_DC": 0.02,\n    "central_deficit_over_actual_reference_peak": 0.02\n  },\n  "unchanged_limits": {\n    "global_energy": 0.01,\n    "temporal": 0.02,\n    "mesh": 0.05,\n    "domain": 0.05,\n    "spectral_Vout": 0.02,\n    "spectral_observed_latency": 0.02\n  },\n  "scope": "775nm existing54/80state spectral comparison over the completed declared physical critical window; not universal microscopic core accuracy.",\n  "scientific_reason": "A vortex zero lies between mesh nodes. Both continuous P1 complex fields already vanish and share the same saved pair/onset. A submesh core displacement changes the nodal minimum and its minimizer. Retain a finite3percent nodal check and additionally require BOTH physical-field L2 and integrated-deficit2percent checks throughout the actual window.",\n  "full_physical_and_critical_coverage_required": true,\n  "completed_primary_runs_required": true,\n  "original_2percent_spectral_failure_must_remain_visible": true,\n  "no_automatic_further_margin_increase": true,\n  "next_action_if_any_new_gate_fails": "Review actual full critical-window evidence; a larger spectral run remains the alternative, not another automatic tolerance increase.",\n  "no_simulation_or_RHS_parameter_change": true,\n  "no_global_energy_correction": true,\n  "baseline_auditor_sha256": "3c08c2534efd2ba8d5803842a9e9cf347b6befc46e68ebddc15adb716156c27a",\n  "baseline_auditor_must_remain_unchanged": true,\n  "implementation": "Explicit SHA-bound analysis wrapper; retain complete baseline result and require original physical checks, real Vout and observed latency. No original failure becomes hidden.",\n  "stage5_admission": false,\n  "PDF_created": false\n}\n'
FROZEN_RESULT = json.loads('{"kind": "spectral", "reference": "nominal", "candidate": "spectral", "passed": false, "checks": {"reference_saved_window_passed": false, "candidate_saved_window_passed": false, "same_frozen_model_sources": true, "same_nonrefined_parameters": true, "actual_declared_refinement": true, "voltage_waveform": true, "electrical_latency_if_observed": true, "temporal_refinement_at_least_two": true, "gap_min_waveform": false}, "relative_limit": 0.02, "stage5_admission_inferred": false, "policy_version": "spatial_observables_v2", "legacy_nodal_policy_passed": false, "legacy_nodal_policy_checks": {"reference_saved_window_passed": false, "candidate_saved_window_passed": false, "same_frozen_model_sources": true, "same_nonrefined_parameters": true, "actual_declared_refinement": true, "voltage_waveform": true, "electrical_latency_if_observed": true, "temporal_refinement_at_least_two": true, "gap_min_waveform": false}, "robust_spatial": {"region_relative_nm": [-60, 60, -40, 40], "x_center": "saved source.center_m[0]", "y_coordinate": "physical", "complex_interpolation": "Saved linear complex triangulation first, magnitude afterwards; no vortex alignment", "time_interpolation": "Piecewise linear complex saved-field comparison only; winding and event are never computed on interpolated states", "quadrature": "Uniform midpoint spatial quadrature on fixed physical120x80nm ROI", "sample_spacing_actual_nm": [0.5, 0.5], "L2_field": {"maximum_relative_error": 0.01889458596621633, "DC_gap_scale_meV": 1.2184332629567942}, "deficit": {"maximum_absolute_difference_meV_nm2": 34.238701612783565, "reference_actual_peak_meV_nm2": 1842.107934696706, "relative_error_actual_peak": 0.018586696777038093, "relative_error_DCarea": 0.002927145479714926, "zero_peak_rule": "Zero difference/zero peak ->0; nonzero difference/zero peak ->failure, no arbitrary floor", "definition": "D_h,nodal=sum exact clipped linear-nodal weights * positive nodal gap deficit; spatial approximation to continuum, not integral of magnitude of complex interpolant"}}, "common_trajectory_interval_ps": [0.0, 12.964999999999817], "common_gap_snapshot_interval_ps": [0.0, 12.899999999999824], "interpolation": "Piecewise linear saved-data comparison only; no physical state or trigger is interpolated by this auditor", "changed_parameters": {"threads": {"reference": 4, "candidate": 2}, "xi_cutoff": {"reference": 32.0, "candidate": 64.0}, "duration_ps": {"reference": null, "candidate": 20.0}, "phonon_nodes": {"reference": 96, "candidate": 144}, "xi_order": {"reference": 6, "candidate": 8}, "transport_order": {"reference": 48, "candidate": 96}}, "forbidden_parameter_changes": {}, "step_refinement_ratio": null, "declared_numerical_control": null, "accepted_step_statistics_after_both_guard_anchors": {"available": true, "interval_ps": [12.399999999999821, 12.964999999999817], "reference": {"available": true, "count": 51, "minimum_ps": 1.8118839761882555e-13, "median_ps": 0.01, "maximum_ps": 0.02, "scope": "Steps ending in common saved trajectory interval; actual accepted h, diagnostic not nominal cap"}, "candidate": {"available": true, "count": 74, "minimum_ps": 0.0025, "median_ps": 0.005, "maximum_ps": 0.02, "scope": "Steps ending in common saved trajectory interval; actual accepted h, diagnostic not nominal cap"}, "scope": "Common interval after both guard anchors; no automatic2x ratio of every adaptive step is inferred"}, "accepted_step_statistics": {"reference": {"available": true, "count": 679, "minimum_ps": 1.8118839761882555e-13, "median_ps": 0.02, "maximum_ps": 0.02, "scope": "Steps ending in common saved trajectory interval; actual accepted h, diagnostic not nominal cap"}, "candidate": {"available": true, "count": 710, "minimum_ps": 0.0025, "median_ps": 0.02, "maximum_ps": 0.02, "scope": "Steps ending in common saved trajectory interval; actual accepted h, diagnostic not nominal cap"}}, "Vout": {"maximum_absolute_difference_mV": 7.979857892580901e-05, "reference_peak_mV": 0.0009448357677463414, "saved_threshold_mV": 0.1, "admission_scale_mV": 0.1, "relative_error": 0.0007979857892580901, "pulse_scaled_error": 0.08445761861465897}, "gap_min": {"maximum_absolute_difference_meV": 0.02773478031920723, "DC_reference_scale_meV": 1.2184475076019632, "relative_error": 0.02276239242656607, "diagnostic_only": false, "former_gate_passed": false}, "latency": {"available": false, "passed": null, "reason": "No pair of saved confirmed electrical events; no latency is assigned"}, "latency_comparison_complete": false}')


def example(gap=.02276239242656607, l2=.01889458596621633, deficit=.018586696777038093):
    checks = {name: True for name in policy.REQUIRED_CHECKS}
    checks["gap_min_waveform"] = gap <= .02
    return {
        "kind": "spectral", "relative_limit": .02,
        "stage5_admission_inferred": False, "policy_version": "spatial_observables_v2",
        "checks": checks, "passed": all(checks.values()),
        "legacy_nodal_policy_checks": copy.deepcopy(checks),
        "legacy_nodal_policy_passed": all(checks.values()),
        "gap_min": {"relative_error": gap, "former_gate_passed": gap <= .02, "diagnostic_only": False},
        "robust_spatial": {"L2_field": {"maximum_relative_error": l2},
                           "deficit": {"relative_error_actual_peak": deficit}},
        "Vout": {"relative_error": .000797985789, "pulse_scaled_error": .0844576},
        "latency": {"available": False, "passed": None},
        "latency_comparison_complete": False,
        "forbidden_parameter_changes": {},
        "sentinel_extra": {"arrays_untouched": [1, 2, 3]},
    }


class SpectralObservablePolicyTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.decision = Path(temporary.name) / "decision.json"
        self.decision.write_bytes(DECISION_PAYLOAD.encode("utf8"))
        baseline = mock.patch.object(policy, "BASELINE_AUDITOR_PATH", BASE / "audit_runs.py")
        baseline.start()
        self.addCleanup(baseline.stop)
    def apply(self, result):
        return policy.apply_spectral_observable_policy(result, decision_path=self.decision)

    def test_measured_nodal_failure_visible_after_prospective_pass(self):
        original = example()
        before = copy.deepcopy(original)
        output = self.apply(original)
        self.assertTrue(output["passed"])
        self.assertEqual(original, before)
        self.assertEqual(output["original_2percent_spectral_result"], before)
        self.assertFalse(output["original_2percent_spectral_result"]["checks"]["gap_min_waveform"])
        self.assertFalse(output["legacy_nodal_policy_passed"])
        self.assertFalse(output["gap_min"]["former_gate_passed"])
        self.assertEqual(output["relative_limit"], .02)
        self.assertEqual(output["gap_min"]["current_gate_relative_limit"], .03)
        self.assertFalse(output["stage5_admission_inferred"])
        output["sentinel_extra"]["arrays_untouched"].append(4)
        self.assertEqual(original, before)

    def test_all_unrelated_failures_are_preserved(self):
        for name in sorted(policy.REQUIRED_CHECKS - {"gap_min_waveform"}):
            with self.subTest(check=name):
                result = example()
                result["checks"][name] = False
                result["passed"] = False
                output = self.apply(result)
                self.assertFalse(output["passed"])
                self.assertFalse(output["checks"][name])
                for key, value in result["checks"].items():
                    if key != "gap_min_waveform":
                        self.assertIs(output["checks"][key], value)

    def test_extra_energy_coverage_guard_failures_preserved(self):
        for name in ["global_energy_independent", "complete_physical_critical_window",
                     "verified_bounded_prospective_guard_policy"]:
            with self.subTest(check=name):
                result = example()
                result["checks"][name] = False
                result["passed"] = False
                self.assertFalse(self.apply(result)["passed"])

    def test_each_new_gate_fails_at_its_declared_margin(self):
        for field, value in [("gap", .0300000001), ("l2", .0200000001), ("deficit", .0200000001)]:
            with self.subTest(field=field):
                self.assertFalse(self.apply(example(**{field: value}))["passed"])

    def test_exact_margins_are_inclusive(self):
        self.assertTrue(self.apply(example(gap=.03, l2=.02, deficit=.02))["passed"])

    def test_already_passing_original_is_not_labelled_failed(self):
        output = self.apply(example(gap=.019))
        self.assertTrue(output["original_2percent_spectral_result"]["passed"])
        self.assertTrue(output["gap_min"]["former_gate_passed"])

    def test_missing_spatial_evidence_fails_closed(self):
        for spatial in [None, {}, {"L2_field": {"maximum_relative_error": .01}}]:
            with self.subTest(spatial=spatial):
                result = example()
                result["robust_spatial"] = spatial
                output = self.apply(result)
                self.assertFalse(output["passed"])
                self.assertFalse(output["checks"]["spectral_common_physical_ROI_available"])

    def test_invalid_robust_numbers_fail_closed(self):
        for value in [math.nan, math.inf, -1, True, ".01", None]:
            for section, key in [("L2_field", "maximum_relative_error"), ("deficit", "relative_error_actual_peak")]:
                with self.subTest(value=value, section=section):
                    result = example()
                    result["robust_spatial"][section][key] = value
                    self.assertFalse(self.apply(result)["passed"])

    def test_wrong_kind_limit_schema_or_admission_rejected(self):
        for key, value in [("kind", "mesh"), ("relative_limit", .03),
                           ("policy_version", "unknown"), ("stage5_admission_inferred", True)]:
            with self.subTest(key=key):
                result = example()
                result[key] = value
                with self.assertRaises(ValueError):
                    self.apply(result)

    def test_missing_or_nonboolean_baseline_checks_rejected(self):
        for name, operation in [
            ("voltage_waveform", lambda checks: checks.pop("voltage_waveform")),
            ("voltage_waveform", lambda checks: checks.__setitem__("voltage_waveform", 1)),
            ("gap_min_waveform", lambda checks: checks.__setitem__("gap_min_waveform", True))]:
            with self.subTest(name=name, operation=operation):
                result = example()
                operation(result["checks"])
                with self.assertRaises(ValueError):
                    self.apply(result)

    def test_inconsistent_pass_or_legacy_rejected(self):
        for target in ["passed", "legacy_nodal_policy_passed", "former_gate_passed", "legacy_gap"]:
            with self.subTest(target=target):
                result = example()
                if target == "former_gate_passed":
                    result["gap_min"][target] = True
                elif target == "legacy_gap":
                    result["legacy_nodal_policy_checks"]["gap_min_waveform"] = True
                else:
                    result[target] = True
                with self.assertRaises(ValueError):
                    self.apply(result)

    def test_invalid_nodal_error_rejected(self):
        for value in [math.nan, math.inf, -1, True, None]:
            with self.subTest(value=value):
                result = example()
                result["gap_min"]["relative_error"] = value
                with self.assertRaises(ValueError):
                    self.apply(result)

    def test_double_application_rejected(self):
        with self.assertRaises(ValueError):
            self.apply(self.apply(example()))

    def test_decision_and_baseline_pins_enforced(self):
        with tempfile.TemporaryDirectory() as folder:
            changed = Path(folder) / "decision.json"
            changed.write_text(self.decision.read_text() + "\n")
            with self.assertRaises(ValueError):
                policy.apply_spectral_observable_policy(example(), decision_path=changed)
        with mock.patch.object(policy, "BASELINE_AUDITOR_SHA256", "0" * 64):
            with self.assertRaises(ValueError):
                self.apply(example())

    def test_frozen_partial_real_comparison_remains_failed(self):
        result = copy.deepcopy(FROZEN_RESULT)
        before = copy.deepcopy(result)
        output = self.apply(result)
        self.assertFalse(output["passed"])
        self.assertTrue(output["checks"]["gap_min_waveform"])
        self.assertTrue(output["checks"]["spectral_central_gap_magnitude_L2"])
        self.assertTrue(output["checks"]["spectral_central_deficit_over_reference_peak"])
        self.assertFalse(output["checks"]["reference_saved_window_passed"])
        self.assertFalse(output["checks"]["candidate_saved_window_passed"])
        self.assertEqual(output["original_2percent_spectral_result"], before)
        self.assertEqual(result, before)


if __name__ == "__main__":
    unittest.main(verbosity=2)
