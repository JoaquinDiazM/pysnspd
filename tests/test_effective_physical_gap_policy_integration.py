"""Saved-observable correction and explicit spectral-policy integration."""
import copy,importlib.util,json,sys
from pathlib import Path
from unittest.mock import patch
import numpy as np
import pytest

HERE=Path(__file__).resolve().parent
BASE=HERE if (HERE/"physical_gap_comparison.py").exists() else HERE.parent/"sandbox/stage5_effective_validation"
REPO=Path.cwd()
if not (REPO/"pysnspd").exists():REPO=Path(__file__).resolve().parents[1]
if str(REPO) not in sys.path:sys.path.insert(0,str(REPO))
def load(name):
    spec=importlib.util.spec_from_file_location("sandbox.stage5_effective_validation._physical_gap_test_"+name,BASE/(name+".py"))
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module);return module
physical=load("physical_gap_comparison")
dispatch=load("guard_protocol_dispatcher")
report=load("build_final_report")
windows=load("window_evidence")
from sandbox.stage5_effective_validation import audit_runs as audit
from tests.test_effective_window_evidence import window,review,primary,HASH
MEV_J=1.602176634e-22

def fixture(kind="mesh"):
    xy=np.array([[-60,-40],[60,-40],[60,40],[-60,40]])*1e-9
    tri=np.array([[0,1,2],[0,2,3]])
    a=dict(coordinates_m=xy,triangles=tri,time_s=np.array([0.,2e-12]),
        delta_J=MEV_J*np.array([[1,1,1,1],[-1,-1,-1,-1]],complex))
    b=dict(a,time_s=np.array([0.,1e-12,2e-12]),
        delta_J=MEV_J*np.ones((3,4),complex))
    runs=[dict(manifest=dict(source=dict(center_m=[0.,0.])),snap=x) for x in (a,b)]
    checks=dict(voltage_waveform=False,central_gap_magnitude_L2=False,central_deficit_over_reference_peak=True)
    row=dict(kind=kind,relative_limit=.05 if kind in {"mesh","domain"} else .02,
        stage5_admission_inferred=False,checks=checks,passed=False,
        gap_min=dict(relative_error=.01),common_gap_interval_ps=[0.,2.],
        robust_spatial=dict(L2_field=dict(maximum_relative_error=1.,DC_gap_scale_meV=1.),
            deficit=dict(relative_error_actual_peak=.01),time_interpolation="previous complex temporal chord"),
        event=dict(confirmed=False,latency_s=None),original_completion=False)
    gap=dict(time_ps=np.array([0.,1.,2.]),field_L2_relative=np.array([0.,1.,0.]),
        candidate_central_nodal_deficit_meV_nm2=np.zeros(3))
    return row,{"sentinel":np.array([3.])},gap,*runs

def corrected(kind="mesh"):
    return physical.correct_comparison(*fixture(kind),audit=audit)

def test_physical_correction_changes_only_L2_and_never_converts_other_failure():
    args=fixture();old=copy.deepcopy(args[0]);out,curves,gap=physical.correct_comparison(*args,audit=audit)
    assert out["original_complex_time_interpolation_result"]==old and args[0]==old
    assert out["checks"]["central_gap_magnitude_L2"] is True
    assert out["checks"]["voltage_waveform"] is False and out["passed"] is False
    assert out["robust_spatial"]["L2_field"]["maximum_relative_error"]==pytest.approx(0.,abs=1e-15)
    assert np.array_equal(gap["original_complex_time_field_L2_relative"],[0.,1.,0.])
    assert out["event"]==old["event"] and out["original_completion"] is False
    assert np.array_equal(curves["sentinel"],[3.])
    physical.validate_physical_gap_correction(out)

@pytest.mark.parametrize("fault",["voltage","deficit","event","original_max","corrected_max"])
def test_physical_correction_rejects_unrelated_or_inconsistent_changes(fault):
    out,_,_=corrected()
    if fault=="voltage":out["checks"]["voltage_waveform"]=True
    if fault=="deficit":out["robust_spatial"]["deficit"]["relative_error_actual_peak"]=0.
    if fault=="event":out["event"]["confirmed"]=True
    if fault=="original_max":out["original_complex_time_interpolation_result"]["robust_spatial"]["L2_field"]["maximum_relative_error"]=.9
    if fault=="corrected_max":out["robust_spatial"]["L2_field"]["maximum_relative_error"]=.1
    with pytest.raises(ValueError):physical.validate_physical_gap_correction(out)

def test_guard_original_verification_precedes_both_full_and_finite_correction():
    calls=[]
    runs=[dict(manifest={},history={},snap={}) for _ in range(2)]
    proof=dict(reference=dict(segments=[]),candidate=dict(segments=[]),guard_law=dict(nominal_anchor_time_ps=22.435))
    row=dict(declared_numerical_control=proof,checks=dict(reference_saved_window_passed=False,same_nonrefined_parameters=True),
        passed=False,legacy_nodal_policy_checks={},common_trajectory_interval_ps=[0.,20.],forbidden_parameter_changes={})
    def full(*a):calls.append("verify_original");return copy.deepcopy(row),{},{}
    def finite(*a):calls.append("finite_metrics");return copy.deepcopy(row),{},{}
    def correction(result,curves,gap,*runs):
        calls.append("correct_physical");result["correction_sentinel"]=True;return result,curves,gap
    def policy(kind,result,*a,**kw):
        calls.append("policy");assert result["correction_sentinel"];return result
    with patch.object(dispatch.guard,"compare_with_bounded_guard",side_effect=full),patch.object(dispatch.audit,"compare_runs",side_effect=finite),patch.object(dispatch.guard,"compare_guard_law",return_value={}),patch.object(dispatch,"correct_physical_gap",side_effect=correction),patch.object(dispatch,"apply_verified_spectral_policy",side_effect=policy):
        comparator,fullrow=dispatch.prepare_verified_comparator("temporal",*runs,"unused")
        out,_,_=comparator("temporal",dict(runs[0]),dict(runs[1]))
    assert calls==["verify_original","correct_physical","policy","finite_metrics","correct_physical","policy"]
    assert out["declared_numerical_control"]["full_original_comparison"]==fullrow
    assert out["checks"]["reference_saved_window_passed"] is False and out["passed"] is False

def test_exploratory_wrapper_keeps_original_failures_and_deficit_failure(monkeypatch):
    policy=load("exploratory_spectral_policy")
    decision=dict(historical_failures_preserved={"sentinel":"original immutable"})
    monkeypatch.setattr(policy,"load_decision",lambda p:(Path("decision.json"),decision))
    row=dict(kind="spectral",relative_limit=.02,stage5_admission_inferred=False,
        policy_version="prospective_spectral_observables_v1",
        spectral_observable_policy=dict(spectral_limits=copy.deepcopy(policy.OLD_LIMITS)),
        physical_gap_interpolation_correction={"sentinel":True},verified_spectral_application_scope={"sentinel":True},
        gap_min=dict(relative_error=.12),robust_spatial=dict(L2_field=dict(maximum_relative_error=.028),deficit=dict(relative_error_actual_peak=.0596)),
        checks={name:True for name in ["spectral_common_physical_ROI_available","voltage_waveform","electrical_latency_if_observed","same_frozen_model_sources","same_nonrefined_parameters"]})
    row["checks"].update(gap_min_waveform=False,spectral_central_gap_magnitude_L2=False,
        spectral_central_deficit_over_reference_peak=False,reference_saved_window_passed=False,candidate_saved_window_passed=False)
    row["passed"]=False
    out=policy.apply_exploratory_spectral_budget(row,decision_path="unused")
    assert out["original_3_2_2_spectral_result"]==row
    assert out["checks"]["spectral_central_gap_magnitude_L2"] is True
    assert out["checks"]["spectral_central_deficit_over_reference_peak"] is False
    assert "gap_min_waveform" not in out["checks"] and out["gap_min"]["diagnostic_only"] is True
    assert out["checks"]["reference_saved_window_passed"] is False and out["passed"] is False
    bad=copy.deepcopy(row);bad["kind"]="temporal"
    with pytest.raises(ValueError):policy.apply_exploratory_spectral_budget(bad,decision_path="unused")

def test_window_bound_keeps_2percent_electrical_5percent_spatial_and_manual_review(monkeypatch):
    w=window();r=review();row=w["original_comparison_result"]
    row["kind"]="spectral";row["robust_spatial"]={};row["spectral_observable_policy"]={}
    row["checks"].pop("common_physical_ROI_available");row["checks"].pop("central_gap_magnitude_L2");row["checks"].pop("central_deficit_over_reference_peak")
    row["checks"].update(spectral_common_physical_ROI_available=True,spectral_central_gap_magnitude_L2=True,spectral_central_deficit_over_reference_peak=True)
    stub=type("Stub",(),dict(physical_gap_module=staticmethod(lambda:type("Physical",(),dict(validate_physical_gap_correction=staticmethod(lambda r:{})))()),
        validate_spectral_policy_application=staticmethod(lambda r:dict(nodal_minimum_over_DC=None,central_field_L2_over_DC=.05,central_deficit_over_actual_reference_peak=.05))))()
    with patch.dict(sys.modules,{"sandbox.stage5_effective_validation.guard_protocol_dispatcher":stub}):
        q=r["windows"][0];q.update(decision="QUANTIFIED_LATER_BOUND",critical_interval_end_ps=25,
            later_observable_bounds={name:dict(maximum_relative_error=.01,declared_relative_limit=lim,expression="unit fixture",physical_argument="unit fixture",evidence=copy.deepcopy(q["evidence"])) for name,lim in [("Vout",.02),("central_gap_L2",.05),("central_deficit",.05)]})
        assert windows.validate_window_evidence(w,HASH,r,verify_evidence=lambda e:None,primary_summary=primary())["passed"]
        q["later_observable_bounds"]["Vout"]["declared_relative_limit"]=.05
        with pytest.raises(ValueError):windows.validate_window_evidence(w,HASH,r,verify_evidence=lambda e:None,primary_summary=primary())
        q["later_observable_bounds"]["Vout"]["declared_relative_limit"]=.02;r["manual_scientific_review"]=False
        with pytest.raises(ValueError):windows.validate_window_evidence(w,HASH,r,verify_evidence=lambda e:None,primary_summary=primary())

def test_report_stage5_false_still_blocks_before_output(tmp_path):
    path=tmp_path/"blocked.json";path.write_text('{"stage5":false}')
    output=tmp_path/"forbidden"
    with pytest.raises(ValueError,match="stage5: true"):report.main(["--admission",str(path),"--build-root",str(output),"--compile"])
    assert not output.exists()

def test_report_history_recomputes_prior_policies_on_current_final_window():
    old2=dict(checks=dict(gap_min_waveform=False,reference_saved_window_passed=False,candidate_saved_window_passed=False))
    old3=dict(checks=dict(gap_min_waveform=False,spectral_central_gap_magnitude_L2=False,spectral_central_deficit_over_reference_peak=True),
        original_2percent_spectral_result=old2,common_gap_snapshot_interval_ps=[0.,25.],
        gap_min=dict(relative_error=.08),robust_spatial=dict(L2_field=dict(maximum_relative_error=.04),deficit=dict(relative_error_actual_peak=.03)))
    row=dict(exploratory_spectral_spatial_policy=dict(decision=dict(path="fixture")),
        original_3_2_2_spectral_result=old3)
    decision=dict(current_provisional_metrics=dict(interval_gap_ps=[0.,15.7],
        nodal_minimum_over_DC=dict(value=.12),field_L2_over_DC=dict(value=.028),
        deficit_over_partial_actual_peak=dict(value=.0596)),
        historical_failures_preserved=dict(original_interpolation_metric=dict(temporal_original_L2=.381673535)))
    stub=type("Dispatch",(),dict(exploratory_policy_module=staticmethod(lambda:type("Policy",(),dict(load_decision=staticmethod(lambda p:(Path(p),decision))))())))()
    with patch.dict(sys.modules,{"sandbox.stage5_effective_validation.guard_protocol_dispatcher":stub}):
        text=report.spectral_policy_history_latex([row])
    assert r"0--25 & 8 & 4 & 3 & falla & falla" in text
    assert "antecedente provisional 0--15.7" in text
    assert "artefacto del análisis anterior, no como fallo físico" in text
