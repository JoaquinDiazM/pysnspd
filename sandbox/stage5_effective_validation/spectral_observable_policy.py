"""SHA-bound prospective spectral observable policy; saved-data analysis only.

The caller must supply a comparison produced by the verified pinned comparator
(and its verified guard/window layers). Hashing that dependency here checks its
immutability, not the provenance of an arbitrary input dictionary. No state,
run status, physical source, numerical parameter or global admission is changed.
"""
from __future__ import annotations

import copy
import hashlib
import json
import math
from pathlib import Path

DECISION_SHA256 = "8c62658c36336a2e19cb0d37bef6b650cf5d44da734e586d3fe1d11b1ce6ed33"
BASELINE_AUDITOR_PATH = Path("/home/jdiaz/pysnspd/sandbox/stage5_effective_validation/audit_runs.py")
BASELINE_AUDITOR_SHA256 = "3c08c2534efd2ba8d5803842a9e9cf347b6befc46e68ebddc15adb716156c27a"
SCHEMA = "pysnspd.spectral_observable_policy_application.v1"
REQUIRED_CHECKS = {
    "reference_saved_window_passed", "candidate_saved_window_passed",
    "same_frozen_model_sources", "same_nonrefined_parameters",
    "actual_declared_refinement", "voltage_waveform",
    "electrical_latency_if_observed", "temporal_refinement_at_least_two",
    "gap_min_waveform",
}
NEW_CHECKS = {
    "spectral_common_physical_ROI_available",
    "spectral_central_gap_magnitude_L2",
    "spectral_central_deficit_over_reference_peak",
}
SPECTRAL_LIMITS = {
    "nodal_minimum_over_DC": .03,
    "central_field_L2_over_DC": .02,
    "central_deficit_over_actual_reference_peak": .02,
}
UNCHANGED_LIMITS = {
    "global_energy": .01, "temporal": .02, "mesh": .05,
    "domain": .05, "spectral_Vout": .02, "spectral_observed_latency": .02,
}


def _require(condition, message):
    if not condition:
        raise ValueError(message)


def _number(value, label):
    _require(not isinstance(value, bool) and isinstance(value, (int, float))
             and math.isfinite(value) and value >= 0,
             "A finite nonnegative number is required: " + label)
    return float(value)


def _sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def load_decision(decision_path):
    """Read exactly the adopted immutable decision and the unchanged comparator."""
    path = Path(decision_path).resolve()
    payload = path.read_bytes()
    _require(hashlib.sha256(payload).hexdigest() == DECISION_SHA256,
             "Prospective decision SHA mismatch")
    _require(_sha256(BASELINE_AUDITOR_PATH) == BASELINE_AUDITOR_SHA256,
             "Pinned baseline comparator changed")
    decision = json.loads(payload)
    expected = {
        "schema": "pysnspd.prospective_spectral_observable_policy.v1",
        "status": "ADOPTED_PRACTICAL_SPECTRAL_OBSERVABLE_POLICY_NOT_MODEL_ADMISSION",
        "spectral_limits": SPECTRAL_LIMITS,
        "unchanged_limits": UNCHANGED_LIMITS,
        "baseline_auditor_sha256": BASELINE_AUDITOR_SHA256,
        "baseline_auditor_must_remain_unchanged": True,
        "full_physical_and_critical_coverage_required": True,
        "completed_primary_runs_required": True,
        "original_2percent_spectral_failure_must_remain_visible": True,
        "no_automatic_further_margin_increase": True,
        "no_simulation_or_RHS_parameter_change": True,
        "no_global_energy_correction": True,
        "stage5_admission": False,
        "PDF_created": False,
    }
    _require(all(decision.get(key) == value for key, value in expected.items()),
             "Prospective decision does not match the fixed adopted contract")
    return path, decision


def _optional_spatial_metric(spatial, section, key):
    """Invalid/missing spatial evidence fails the new gates without inventing data."""
    try:
        return _number(spatial[section][key], section + "." + key)
    except (KeyError, TypeError, ValueError):
        return None


def apply_spectral_observable_policy(result, *, decision_path):
    """Return a deep copy with the explicit 3% nodal plus 2% field/deficit policy.

    All preexisting checks except gap_min_waveform retain their exact values.
    In particular incomplete windows, energy failures, source mismatches, signal
    or latency failures remain failures. This function never admits stage 5.
    """
    path, decision = load_decision(decision_path)
    _require(isinstance(result, dict), "A verified comparison dictionary is required")
    _require(result.get("kind") == "spectral", "Only a spectral comparison is eligible")
    _require(result.get("policy_version") == "spatial_observables_v2",
             "Unexpected baseline observable schema")
    _require(result.get("relative_limit") == .02
             and not isinstance(result.get("relative_limit"), bool),
             "The baseline spectral voltage/latency limit must remain 2 percent")
    _require(result.get("stage5_admission_inferred") is False,
             "Comparison may not infer model admission")
    _require("spectral_observable_policy" not in result
             and "original_2percent_spectral_result" not in result,
             "The prospective policy may be applied only once")
    checks = result.get("checks")
    _require(isinstance(checks, dict) and REQUIRED_CHECKS <= checks.keys(),
             "Complete baseline physical/completion checks are required")
    _require(not NEW_CHECKS.intersection(checks),
             "Prospective gates already exist in the input")
    _require(all(type(value) is bool for value in checks.values()),
             "Checks must be explicit Boolean decisions")
    _require(type(result.get("passed")) is bool
             and result["passed"] == all(checks.values()),
             "Input pass flag does not match its original checks")
    gap = result.get("gap_min")
    _require(isinstance(gap, dict), "Baseline nodal gap evidence is required")
    gap_error = _number(gap.get("relative_error"), "gap_min.relative_error")
    original_gap_passed = gap_error <= .02
    _require(checks["gap_min_waveform"] == original_gap_passed
             and gap.get("former_gate_passed") is original_gap_passed
             and gap.get("diagnostic_only") is False,
             "Original 2 percent nodal gate is inconsistent with its measured error")
    legacy_checks = result.get("legacy_nodal_policy_checks")
    _require(isinstance(legacy_checks, dict)
             and all(type(value) is bool for value in legacy_checks.values())
             and legacy_checks.get("gap_min_waveform") is original_gap_passed
             and result.get("legacy_nodal_policy_passed") is all(legacy_checks.values()),
             "The original legacy nodal policy must remain explicit and consistent")
    spatial = result.get("robust_spatial")
    l2_error = _optional_spatial_metric(spatial, "L2_field", "maximum_relative_error")
    deficit_error = _optional_spatial_metric(spatial, "deficit", "relative_error_actual_peak")
    spatial_available = l2_error is not None and deficit_error is not None

    original = copy.deepcopy(result)
    output = copy.deepcopy(result)
    output["original_2percent_spectral_result"] = original
    output["checks"]["gap_min_waveform"] = gap_error <= .03
    output["checks"]["spectral_common_physical_ROI_available"] = spatial_available
    output["checks"]["spectral_central_gap_magnitude_L2"] = (
        l2_error is not None and l2_error <= .02)
    output["checks"]["spectral_central_deficit_over_reference_peak"] = (
        deficit_error is not None and deficit_error <= .02)
    output["passed"] = all(output["checks"].values())
    output["policy_version"] = "prospective_spectral_observables_v1"
    output["gap_min"]["current_gate_relative_limit"] = .03
    output["gap_min"]["current_gate_passed"] = output["checks"]["gap_min_waveform"]
    output["spectral_observable_policy"] = {
        "schema": SCHEMA,
        "decision": {"path": str(path), "sha256": DECISION_SHA256},
        "baseline_dependency": {
            "path": str(BASELINE_AUDITOR_PATH),
            "sha256": BASELINE_AUDITOR_SHA256,
            "immutability_verified": True,
            "input_provenance": "Established by caller/dispatcher; not inferred from dictionary"},
        "spectral_limits": copy.deepcopy(SPECTRAL_LIMITS),
        "unchanged_limits": copy.deepcopy(UNCHANGED_LIMITS),
        "scope": decision["scope"],
        "original_2percent_nodal_gate_passed": original_gap_passed,
        "original_2percent_comparison_passed": original["passed"],
        "historical_failures_preserved": copy.deepcopy(decision["original_failures_preserved"]),
        "all_original_checks_except_nodal_gate_preserved": True,
        "legacy_nodal_policy_preserved": True,
        "full_physical_and_critical_coverage_required": True,
        "completed_primary_runs_required": True,
        "no_state_or_parameter_change": True,
        "stage5_admission_inferred": False,
    }
    return output
