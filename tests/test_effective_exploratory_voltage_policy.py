"""Post-data Vout budget only: no evolution, spatial calculation, plots or PDF."""
import copy,hashlib,importlib.util,json,os,sys
from pathlib import Path
from unittest.mock import patch
import pytest

ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT/"sandbox/stage5_effective_validation"
def load(name,path):
    spec=importlib.util.spec_from_file_location(name,path)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module);return module
policy=load("post_data_voltage_policy",BASE/"exploratory_spectral_policy.py")
dispatch=load("sandbox.stage5_effective_validation.post_data_voltage_dispatch",BASE/"guard_protocol_dispatcher.py")
report=load("post_data_voltage_report",BASE/"build_final_report.py")
DECISION_PAYLOAD='{\n  "schema": "pysnspd.post_data_exploratory_spectral_voltage_decision.v1",\n  "status": "ADOPTED_POST_DATA_VOLTAGE_POLICY_NOT_MODEL_ADMISSION",\n  "created_UTC": "2026-10-12T05:57:27.126119+00:00",\n  "approval_state": "ROOT_EXPLICITLY_ADOPTED_CANDIDATE_ENCODING_FOR_REVIEW",\n  "scope": "775nm_current_54_to_80_bundle_only",\n  "policy_timing": "AFTER_COMPLETED_DATA_REVIEW",\n  "voltage_relative_limit": 0.05,\n  "prior_voltage_relative_limit": 0.02,\n  "latency_relative_limit": 0.02,\n  "global_energy_relative_limit": 0.01,\n  "spatial_limits": {\n    "central_field_L2_over_DC": 0.05,\n    "central_deficit_over_actual_reference_peak": 0.05,\n    "nodal_minimum_over_DC": null\n  },\n  "normalization": "maximum absolute saved Vout difference / max(actual reference pulse peak, saved trigger threshold)",\n  "coverage_and_later_bounds_unchanged": true,\n  "all_original_rows_and_maxima_retained": true,\n  "prior_2percent_voltage_failure_retained": true,\n  "no_automatic_further_margin_increase": true,\n  "no_RHS_or_parameters_changed": true,\n  "stage5_admission": false,\n  "production": {\n    "default_electron_states": 54,\n    "optional_spectral_sensitivity_control_states": 80,\n    "control_is_not_continuum_truth": true,\n    "scope_775nm_only": true\n  },\n  "closure_statement_only_after_four_axes_admitted": "Resolución54 validada en cuatro ejes a775nm dentro de estos márgenes;80 es un control de sensibilidad, no la verdad continua.",\n  "engineering_reason": "Revision explicitly made after observing a3.7792uV difference on a145.4982uV pulse (2.5974percent). Both saved triggers confirm, latency differs0.31353ps (about0.982percent), independent field/deficit and energy checks pass.5percent is our exploratory engineering signal budget, matching mesh/domain, not a published tolerance or an estimate of asymptotic error.",\n  "limitations": [\n    "No universal accuracy or continuum-error estimate.",\n    "No transfer of four-axis775nm validation to1550nm.",\n    "No quantitative guarantee near an untested detection boundary, other bias/material or source calibrators.",\n    "Critical-window coverage and all physical, source, support, current-sign and completion checks remain mandatory."\n  ],\n  "preserved_dependencies": {\n    "baseline": {\n      "path": "/home/jdiaz/pysnspd/sandbox/stage5_effective_validation/audit_runs.py",\n      "sha256": "3c08c2534efd2ba8d5803842a9e9cf347b6befc46e68ebddc15adb716156c27a"\n    },\n    "prior_observable_decision": {\n      "path": "/home/jdiaz/scratch/stage5_spectral_mesh_scope_review_20261011/root_prospective_spectral_policy_decision.json",\n      "sha256": "8c62658c36336a2e19cb0d37bef6b650cf5d44da734e586d3fe1d11b1ce6ed33"\n    },\n    "prior_spatial_decision": {\n      "path": "/home/jdiaz/scratch/stage5_exploratory_spatial_budget_candidate_20261011/root_prospective_exploratory_spatial_decision.json",\n      "sha256": "e683151e6e7b678d95eb6c133391b2ea84ca2ba81cd5a54b1a21cb9f16624e23"\n    }\n  },\n  "observed_closed_comparison": {\n    "window_path": "/home/jdiaz/scratch/stage5_closed_spectral_comparison_0_36p6_20261012/window.json",\n    "window_sha256": "143898ef77cd80fe77490aa4af720daccc63dfaa87939693212fcd17dc840fad",\n    "window_ps": [\n      0.0,\n      36.6\n    ],\n    "Vout_maximum_absolute_difference_uV": 3.779204888867732,\n    "Vout_reference_peak_uV": 145.49816411004997,\n    "Vout_relative_error": 0.025974244499808716,\n    "latency_difference_ps": 0.3135279760594457,\n    "latency_tolerance_ps": 0.6388263372329995,\n    "field_L2_over_DC": 0.03085254812735303,\n    "deficit_over_actual_reference_peak": 0.019712005758727242,\n    "original_voltage_2percent_passed": false\n  }\n}\n'

def decision(tmp_path):
    path=tmp_path/"voltage_decision.json";path.write_text(DECISION_PAYLOAD);return path

def example(error=.025974244499808716):
    d=json.loads(DECISION_PAYLOAD)
    checks={k:True for k in ["reference_saved_window_passed","candidate_saved_window_passed",
        "same_frozen_model_sources","same_nonrefined_parameters","actual_declared_refinement",
        "electrical_latency_if_observed","temporal_refinement_at_least_two",
        "spectral_common_physical_ROI_available","spectral_central_gap_magnitude_L2",
        "spectral_central_deficit_over_reference_peak"]}
    checks["voltage_waveform"]=error<=.02
    return dict(kind="spectral",relative_limit=.02,stage5_admission_inferred=False,
        policy_version="exploratory_spectral_spatial_observables_v1",checks=checks,
        passed=all(checks.values()),spectral_observable_policy={"historical":"unchanged"},
        exploratory_spectral_spatial_policy=dict(decision=d["preserved_dependencies"]["prior_spatial_decision"],spectral_limits=copy.deepcopy(policy.LIMITS)),
        verified_spectral_application_scope=dict(cases=copy.deepcopy(dispatch.SPECTRAL_CASES),
            global_energy_relative_limit=.01,full_original_comparison_verified_before_policy=True,
            stage5_admission_inferred=False),
        Vout=dict(relative_error=error,pulse_scaled_error=error,maximum_absolute_difference_mV=.003779204888867732),
        latency=dict(available=True,passed=True,absolute_difference_ps=.3135279760594457),
        robust_spatial=dict(L2_field=dict(maximum_relative_error=.03085254812735303),
            deficit=dict(relative_error_actual_peak=.019712005758727242)),
        gap_min=dict(relative_error=.12,diagnostic_only=True,former_gate_passed=False),
        energy=dict(maximum_relative_defect=.0066423,relative_limit=.01),
        sentinel=dict(rows=[1,2,3],maxima=[.1,.2]))

def scoped_example(error=.025974244499808716):
    row=example(error)
    for values in row["verified_spectral_application_scope"]["cases"].values():values["wavelength_nm"]=775.
    return row

def apply(row,tmp_path):return policy.apply_exploratory_voltage_budget(row,decision_path=decision(tmp_path))

def test_only_voltage_gate_changes_and_complete_old_failure_remains(tmp_path):
    row=scoped_example();before=copy.deepcopy(row);out=apply(row,tmp_path)
    assert out["passed"] is True and row==before
    assert {k for k in row["checks"] if out["checks"][k]!=row["checks"][k]}=={"voltage_waveform"}
    assert out["original_2percent_voltage_spectral_result"]==before
    assert out["original_2percent_voltage_spectral_result"]["passed"] is False
    assert out["relative_limit"]==.02 and out["stage5_admission_inferred"] is False
    for key in row:
        if key not in {"checks","passed","policy_version"}:assert out[key]==row[key]
    out["sentinel"]["rows"].append(4);assert row==before

@pytest.mark.parametrize("check",["electrical_latency_if_observed","global_energy_independent",
    "complete_physical_critical_window","same_frozen_model_sources","candidate_saved_window_passed"])
def test_unrelated_failures_remain_failures(tmp_path,check):
    row=scoped_example();row["checks"][check]=False;row["passed"]=False
    out=apply(row,tmp_path)
    assert out["checks"][check] is False and out["passed"] is False

@pytest.mark.parametrize("error,passed",[(.02,True),(.05,True),(.050000001,False)])
def test_finite_voltage_margin_is_inclusive(tmp_path,error,passed):
    assert apply(scoped_example(error),tmp_path)["passed"] is passed

@pytest.mark.parametrize("error",[float("nan"),float("inf"),-1.,True,".03",None])
def test_invalid_voltage_metric_rejected(tmp_path,error):
    row=scoped_example();row["Vout"]["relative_error"]=error
    with pytest.raises(ValueError):apply(row,tmp_path)

@pytest.mark.parametrize("fault",["kind","baseline_limit","scope_wavelength","scope_states",
    "energy_limit","spatial_policy","original_gate","inferred_admission","repeated"])
def test_invalid_scope_or_inconsistent_original_rejected(tmp_path,fault):
    row=scoped_example()
    if fault=="kind":row["kind"]="temporal"
    elif fault=="baseline_limit":row["relative_limit"]=.05
    elif fault=="scope_wavelength":row["verified_spectral_application_scope"]["cases"]["candidate"]["wavelength_nm"]=1550.
    elif fault=="scope_states":row["verified_spectral_application_scope"]["cases"]["candidate"]["electron_states"]=54
    elif fault=="energy_limit":row["verified_spectral_application_scope"]["global_energy_relative_limit"]=.02
    elif fault=="spatial_policy":row["exploratory_spectral_spatial_policy"]["decision"]["sha256"]="f"*64
    elif fault=="original_gate":row["checks"]["voltage_waveform"]=True;row["passed"]=True
    elif fault=="inferred_admission":row["stage5_admission_inferred"]=True
    elif fault=="repeated":row=apply(row,tmp_path)
    with pytest.raises(ValueError):apply(row,tmp_path)

def test_decision_tamper_rejected(tmp_path):
    path=decision(tmp_path);path.write_text(path.read_text()+" ")
    with pytest.raises(ValueError,match="SHA mismatch"):
        policy.apply_exploratory_voltage_budget(scoped_example(),decision_path=path)

def test_preserved_dependency_tamper_rejected_without_touching_original(tmp_path):
    payload=json.loads(DECISION_PAYLOAD)
    baseline=payload["preserved_dependencies"]["baseline"]
    copy_path=tmp_path/"baseline_copy.py";copy_path.write_bytes(Path(baseline["path"]).read_bytes())
    baseline["path"]=str(copy_path)
    path=tmp_path/"isolated_dependency_decision.json";path.write_text(json.dumps(payload))
    digest=hashlib.sha256(path.read_bytes()).hexdigest()
    # Authenticate this test-only path before altering its preserved dependency.
    with patch.object(policy,"VOLTAGE_DECISION_SHA256",digest):
        assert policy.apply_exploratory_voltage_budget(scoped_example(),decision_path=path)["passed"] is True
        copy_path.write_bytes(copy_path.read_bytes()+b"\n# test alteration\n")
        with pytest.raises(ValueError,match="dependency changed"):
            policy.apply_exploratory_voltage_budget(scoped_example(),decision_path=path)

def test_cli_forwards_optional_decision_to_finite_audit(tmp_path):
    recorded={}
    def finite(*a,**kw):recorded.update(kw);return dict(passed=False),{},{}
    with patch.object(dispatch,"audit_finite_window",side_effect=finite),patch.object(dispatch.audit,"write_table"):
        dispatch.main(["--kind","spectral","--reference",str(tmp_path/"r"),"--candidate",str(tmp_path/"c"),
            "--guard-protocol",str(tmp_path/"g"),"--window-declaration",str(tmp_path/"w"),
            "--spectral-policy-decision",str(tmp_path/"old"),"--exploratory-spectral-decision",str(tmp_path/"spatial"),
            "--exploratory-voltage-decision",str(tmp_path/"voltage"),"--output-dir",str(tmp_path/"out")])
    assert recorded["exploratory_voltage_decision_path"]==tmp_path/"voltage"
    assert not (tmp_path/"out").exists()

def test_voltage_decision_without_prior_policies_is_rejected_before_file_read():
    with pytest.raises(ValueError,match="preserve both"):
        dispatch.apply_verified_spectral_policy("spectral",{}, {},{},None,exploratory_voltage_decision_path="unused")

def test_report_reads_voltage5_separately_from_latency2_and_labels_post_data(tmp_path):
    out=apply(scoped_example(),tmp_path)
    margins=dict(policy.LIMITS,voltage_waveform_relative_limit=.05)
    with patch.dict(sys.modules,{"sandbox.stage5_effective_validation.guard_protocol_dispatcher":dispatch}),\
         patch.object(dispatch,"validate_spectral_policy_application",return_value=margins):
        displayed=report.comparison_admission_metrics(out)
    assert displayed["Vout_margin_ratio"]==pytest.approx(out["Vout"]["relative_error"]/.05)
    assert displayed["voltage_gate_relative_limit"]==.05 and displayed["latency_gate_relative_limit"]==.02
    text=report.convergence_caption_latex([out])
    for word in [r"señal 5\%",r"latencia 2\%","posterior a los datos","fallo previo","no verdad continua"]:assert word in text

def test_coverage_and_later_Vout_bound_remain2percent(tmp_path):
    from sandbox.stage5_effective_validation import window_evidence
    from tests.test_effective_window_evidence import window,review,primary,HASH,ENTRY
    out=apply(scoped_example(),tmp_path);w=window();w["original_comparison_result"]=out
    r=review();q=r["windows"][0]
    q.update(decision="QUANTIFIED_LATER_BOUND",critical_interval_end_ps=25.,
        later_observable_bounds={name:dict(maximum_relative_error=.01,declared_relative_limit=margin,
            expression="synthetic bound only",physical_argument="fixture only",evidence=[ENTRY])
            for name,margin in [("Vout",.02),("central_gap_L2",.05),("central_deficit",.05)]})
    margins=dict(policy.LIMITS,voltage_waveform_relative_limit=.05)
    physical=type("PhysicalFixture",(),dict(validate_physical_gap_correction=staticmethod(lambda row:{})))()
    with patch.dict(sys.modules,{"sandbox.stage5_effective_validation.guard_protocol_dispatcher":dispatch}),\
         patch.object(dispatch,"validate_spectral_policy_application",return_value=margins),patch.object(dispatch,"physical_gap_module",return_value=physical):
        assert window_evidence.validate_window_evidence(w,HASH,r,verify_evidence=lambda e:None,primary_summary=primary())["passed"] is True
        q["later_observable_bounds"]["Vout"]["declared_relative_limit"]=.05
        with pytest.raises(ValueError,match="later bound"):
            window_evidence.validate_window_evidence(w,HASH,r,verify_evidence=lambda e:None,primary_summary=primary())
        q["later_observable_bounds"]["Vout"]["declared_relative_limit"]=.02
        r["manual_scientific_review"]=False
        with pytest.raises(ValueError,match="manual"):
            window_evidence.validate_window_evidence(w,HASH,r,verify_evidence=lambda e:None,primary_summary=primary())

@pytest.mark.skipif(not os.environ.get("PYSNSPD_TEST_CLOSED_SPECTRAL_WINDOW"),reason="Optional immutable closed-JSON integration; no trajectory or RHS")
def test_real_closed_result_exact_validation_and_builder_input_pin(tmp_path):
    path=Path(os.environ["PYSNSPD_TEST_CLOSED_SPECTRAL_WINDOW"])
    assert hashlib.sha256(path.read_bytes()).hexdigest()=="143898ef77cd80fe77490aa4af720daccc63dfaa87939693212fcd17dc840fad"
    saved=json.loads(path.read_text());row=saved["original_comparison_result"];before=copy.deepcopy(row)
    output=apply(row,tmp_path)
    assert dispatch.validate_spectral_policy_application(output)["voltage_waveform_relative_limit"]==.05
    assert output["original_2percent_voltage_spectral_result"]==before and row==before
    assert saved["passed"] is False and output["passed"] is True
    altered=copy.deepcopy(output);altered["checks"]["electrical_latency_if_observed"]=False;altered["passed"]=False
    with pytest.raises(ValueError,match="transformation changed"):dispatch.validate_spectral_policy_application(altered)
    audit_path=tmp_path/"admitted_synthetic_audit.json"
    audit_path.write_text(json.dumps(dict(comparisons=[output]+[dict(kind=k,passed=True) for k in ["temporal","mesh","domain"]])))
    primary_runs={};fake={}
    for wave in ["775","1550"]:
        run=tmp_path/wave;run.mkdir();primary_runs[wave]=str(run)
        for filename in ["manifest.json","trajectory.npz","snapshots.npz","checkpoint_latest.npz"]:(run/filename).write_text("fixture only, not a trajectory")
        fake[str(run)]=dict(summary=dict(passed=True),manifest=dict(source=dict(wavelength_nm=float(wave)),run_complete=True))
    admission=tmp_path/"synthetic_admission.json";admission.write_text(json.dumps(dict(stage5=True,primary_runs=primary_runs,evidence_audits=[str(audit_path)])))
    from sandbox.stage5_effective_validation import audit_runs
    with patch.dict(sys.modules,{"sandbox.stage5_effective_validation.guard_protocol_dispatcher":dispatch}),\
         patch.object(audit_runs,"audit_run",side_effect=lambda wave,path:fake[str(path)]):
        _,_,_,_,hashes=report.admission_inputs(admission)
    p=tmp_path/"voltage_decision.json"
    assert hashes[str(p)]==policy.VOLTAGE_DECISION_SHA256
