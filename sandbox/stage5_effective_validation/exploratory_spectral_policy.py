"""Explicit SHA-bound exploratory spectral spatial budget; no model admission."""
from __future__ import annotations
import copy,hashlib,json,math
from pathlib import Path
DECISION_SHA="e683151e6e7b678d95eb6c133391b2ea84ca2ba81cd5a54b1a21cb9f16624e23"
SCHEMA="pysnspd.exploratory_spectral_spatial_policy_application.v1"
LIMITS=dict(nodal_minimum_over_DC=None,central_field_L2_over_DC=.05,
    central_deficit_over_actual_reference_peak=.05)
OLD_LIMITS=dict(nodal_minimum_over_DC=.03,central_field_L2_over_DC=.02,
    central_deficit_over_actual_reference_peak=.02)

def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def load_decision(path):
    path=Path(path).resolve();payload=path.read_bytes()
    if hashlib.sha256(payload).hexdigest()!=DECISION_SHA:
        raise ValueError("Adopted exploratory spectral decision SHA mismatch")
    d=json.loads(payload)
    if (d.get("schema")!="pysnspd.prospective_exploratory_spatial_budget.v1"
        or d.get("status")!="ADOPTED_EXPLORATORY_SPECTRAL_SPATIAL_POLICY_NOT_MODEL_ADMISSION"
        or d.get("approval_state")!="ROOT_REVIEWED_AND_ADOPTED"
        or d.get("policy_scope")!="spectral_775nm_current_54_to_80_bundle"
        or d.get("stage5_admission") is not False
        or d["scope_details"]["spectral_spatial_limits_proposed"]!=dict(
            field_L2_over_DC=.05,central_deficit_over_actual_reference_peak=.05)
        or d["scope_details"]["nodal_minimum"]["diagnostic_only"] is not True
        or d["window"]["all_prior_common_window_maxima_retained"] is not True
        or d["window"]["cannot_crop_away_prior_failure"] is not True
        or d["window"]["no_200ps_tail_for_denominator_dilution"] is not True):
        raise ValueError("Unexpected adopted exploratory spatial budget")
    return path,d

def apply_exploratory_spectral_budget(result,*,decision_path):
    path,decision=load_decision(decision_path)
    if (result.get("kind")!="spectral" or result.get("relative_limit")!=.02
        or result.get("stage5_admission_inferred") is not False
        or result.get("policy_version")!="prospective_spectral_observables_v1"
        or result.get("spectral_observable_policy",{}).get("spectral_limits")!=OLD_LIMITS):
        raise ValueError("Verified original3/2/2 spectral comparison required")
    if "exploratory_spectral_spatial_policy" in result:
        raise ValueError("Exploratory spatial budget may be applied only once")
    if "physical_gap_interpolation_correction" not in result:
        raise ValueError("Reviewed physical magnitude correction required")
    if "verified_spectral_application_scope" not in result:
        raise ValueError("Actual775nm54->80 manifest scope must be verified upstream")
    checks=result.get("checks",{})
    needed={"gap_min_waveform","spectral_common_physical_ROI_available",
        "spectral_central_gap_magnitude_L2","spectral_central_deficit_over_reference_peak",
        "voltage_waveform","electrical_latency_if_observed","same_frozen_model_sources",
        "same_nonrefined_parameters","reference_saved_window_passed","candidate_saved_window_passed"}
    if (not needed<=checks.keys() or any(type(v) is not bool for v in checks.values())
        or result.get("passed")!=all(checks.values())):
        raise ValueError("Complete consistent original spectral checks required")
    spatial=result["robust_spatial"];l2=spatial["L2_field"]["maximum_relative_error"]
    deficit=spatial["deficit"]["relative_error_actual_peak"]
    if any(isinstance(x,bool) or not isinstance(x,(float,int)) or not math.isfinite(x) or x<0 for x in (l2,deficit)):
        raise ValueError("Finite physical field and deficit errors required")
    if (checks["spectral_central_gap_magnitude_L2"]!=(l2<=.02)
        or checks["spectral_central_deficit_over_reference_peak"]!=(deficit<=.02)
        or checks["gap_min_waveform"]!=(result["gap_min"]["relative_error"]<=.03)):
        raise ValueError("Original3/2/2 decisions do not match measured errors")
    output=copy.deepcopy(result)
    output["original_3_2_2_spectral_result"]=copy.deepcopy(result)
    output["checks"].pop("gap_min_waveform")
    output["checks"]["spectral_central_gap_magnitude_L2"]=l2<=.05
    output["checks"]["spectral_central_deficit_over_reference_peak"]=deficit<=.05
    output["passed"]=all(output["checks"].values())
    output["gap_min"]["diagnostic_only"]=True
    output["gap_min"]["current_gate_relative_limit"]=None
    output["gap_min"]["current_gate_passed"]=None
    output["policy_version"]="exploratory_spectral_spatial_observables_v1"
    output["exploratory_spectral_spatial_policy"]=dict(schema=SCHEMA,
        decision=dict(path=str(path),sha256=DECISION_SHA),spectral_limits=copy.deepcopy(LIMITS),
        original_3_2_2_result_preserved=True,original_2percent_result_preserved=True,
        original_3_2_2_checks=copy.deepcopy(checks),
        historical_failures_preserved=copy.deepcopy(decision["historical_failures_preserved"]),
        all_original_checks_except_spectral_spatial_and_nodal_preserved=True,
        full_prior_window_maxima_retained=True,no_time_shift_or_vortex_alignment=True,
        no_200ps_tail_for_denominator_dilution=True,production_resolution="PENDING_ROOT_DECISION",
        no_RHS_or_source_or_model_parameter_change=True,stage5_admission_inferred=False)
    return output


VOLTAGE_DECISION_SHA256 = "bda09a6aaf095a175440dc77ac76faef8eaab98072fce2f5b5aa8984eb35efa6"
VOLTAGE_SCHEMA = "pysnspd.exploratory_spectral_voltage_policy_application.v1"

def apply_exploratory_voltage_budget(result, *, decision_path):
    """Post-data Vout5% only; preserve the complete prior2% result and all else."""
    path=Path(decision_path).resolve();payload=path.read_bytes()
    if hashlib.sha256(payload).hexdigest()!=VOLTAGE_DECISION_SHA256:
        raise ValueError("Adopted post-data voltage decision SHA mismatch")
    decision=json.loads(payload)
    expected=dict(schema="pysnspd.post_data_exploratory_spectral_voltage_decision.v1",
        status="ADOPTED_POST_DATA_VOLTAGE_POLICY_NOT_MODEL_ADMISSION",
        scope="775nm_current_54_to_80_bundle_only",policy_timing="AFTER_COMPLETED_DATA_REVIEW",
        voltage_relative_limit=.05,prior_voltage_relative_limit=.02,
        latency_relative_limit=.02,global_energy_relative_limit=.01,
        coverage_and_later_bounds_unchanged=True,prior_2percent_voltage_failure_retained=True,
        all_original_rows_and_maxima_retained=True,no_automatic_further_margin_increase=True,
        no_RHS_or_parameters_changed=True,stage5_admission=False)
    if any(decision.get(k)!=v for k,v in expected.items()):
        raise ValueError("Unexpected adopted exploratory voltage contract")
    dependencies=decision.get("preserved_dependencies",{})
    expected_dependencies=dict(baseline="3c08c2534efd2ba8d5803842a9e9cf347b6befc46e68ebddc15adb716156c27a",
        prior_observable_decision="8c62658c36336a2e19cb0d37bef6b650cf5d44da734e586d3fe1d11b1ce6ed33",
        prior_spatial_decision=DECISION_SHA)
    if set(dependencies)!=set(expected_dependencies):
        raise ValueError("Complete immutable prior voltage-policy dependencies required")
    for name,expected_sha in expected_dependencies.items():
        entry=dependencies[name]
        if entry.get("sha256")!=expected_sha or digest(entry["path"])!=expected_sha:
            raise ValueError("Preserved voltage-policy dependency changed: "+name)
    if "exploratory_spectral_voltage_policy" in result:
        raise ValueError("Exploratory voltage budget may be applied only once")
    if (result.get("kind")!="spectral" or result.get("relative_limit")!=.02
        or result.get("stage5_admission_inferred") is not False
        or result.get("policy_version")!="exploratory_spectral_spatial_observables_v1"):
        raise ValueError("Verified prior exploratory spatial comparison with Vout2% required")
    spatial_meta=result.get("exploratory_spectral_spatial_policy",{})
    entry=spatial_meta.get("decision",{})
    if entry.get("sha256")!=DECISION_SHA or spatial_meta.get("spectral_limits")!=LIMITS:
        raise ValueError("Prior adopted spatial policy must remain exact")
    load_decision(entry["path"])
    proof=result.get("verified_spectral_application_scope",{})
    cases=dict(reference=dict(wavelength_nm=775.,electron_states=54,xi_cutoff=32.,xi_order=6,phonon_nodes=96,transport_order=48),
        candidate=dict(wavelength_nm=775.,electron_states=80,xi_cutoff=64.,xi_order=8,phonon_nodes=144,transport_order=96))
    if (proof.get("cases")!=cases or proof.get("global_energy_relative_limit")!=.01
        or proof.get("full_original_comparison_verified_before_policy") is not True
        or proof.get("stage5_admission_inferred") is not False):
        raise ValueError("Actual verified775nm54->80 scope required")
    checks=result.get("checks",{});error=result.get("Vout",{}).get("relative_error")
    if (not isinstance(error,(float,int)) or isinstance(error,bool)
        or not math.isfinite(error) or error<0
        or not checks or any(type(v) is not bool for v in checks.values())
        or checks.get("voltage_waveform")!=(error<=.02)
        or result.get("passed")!=all(checks.values())):
        raise ValueError("Finite Vout error and consistent original2% checks required")
    output=copy.deepcopy(result)
    output["original_2percent_voltage_spectral_result"]=copy.deepcopy(result)
    output["checks"]["voltage_waveform"]=error<=.05
    output["passed"]=all(output["checks"].values())
    output["policy_version"]="exploratory_spectral_voltage_observables_v1"
    output["exploratory_spectral_voltage_policy"]=dict(schema=VOLTAGE_SCHEMA,
        decision=dict(path=str(path),sha256=VOLTAGE_DECISION_SHA256),
        policy_timing="AFTER_COMPLETED_DATA_REVIEW",voltage_relative_limit=.05,
        prior_voltage_relative_limit=.02,latency_relative_limit=.02,
        global_energy_relative_limit=.01,coverage_and_later_bounds_unchanged=True,
        only_active_check_changed="voltage_waveform",original_2percent_voltage_result_preserved=True,
        production_default_electron_states=54,sensitivity_control_electron_states=80,
        no_continuum_or_universal_accuracy_claim=True,no_RHS_or_parameter_change=True,
        stage5_admission_inferred=False)
    return output
