"""Spectral policy wiring: actual manifest scope, original proof and coverage."""
import copy,importlib.util,json,sys
from pathlib import Path
from unittest.mock import patch
import pytest
HERE=Path(__file__).resolve().parent
REPO=HERE.parent if (HERE.parent/"pysnspd").exists() else Path.cwd()
BASE=HERE if (HERE/"guard_protocol_dispatcher.py").exists() else REPO/"sandbox/stage5_effective_validation"
sys.path.insert(0,str(REPO))
def load(name,path):
    spec=importlib.util.spec_from_file_location(name,path)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    return module
dispatch=load("sandbox.stage5_effective_validation.guard_dispatch_candidate",BASE/"guard_protocol_dispatcher.py")
report=load("candidate_spectral_report",BASE/"build_final_report.py")
windows=load("candidate_spectral_window_evidence",BASE/"window_evidence.py")
# Existing scope/policy fixtures deliberately lack spatial arrays. Isolate the
# new physical layer here; its actual transformation has focused integration tests.
dispatch.correct_physical_gap=lambda result,curves,gap,*runs:(result,curves,gap)
dispatch.physical_gap_module=lambda:type("PhysicalFixture",(),dict(validate_physical_gap_correction=staticmethod(lambda r:{})))()
from tests.test_effective_spectral_observable_policy import example,DECISION_PAYLOAD
from tests.test_effective_window_evidence import window,review,primary,HASH,ENTRY
def inputs(tmp_path):
    decision=tmp_path/"decision.json";decision.write_bytes(DECISION_PAYLOAD.encode())
    runs={}
    for role,values in dispatch.SPECTRAL_CASES.items():
        location=tmp_path/role;location.mkdir()
        identity=dict(wavelength_nm=775.,**{k:v for k,v in values.items() if k!="electron_states"})
        manifest=dict(source=dict(wavelength_nm=775.),parameter_identity=identity,
            kinetics=dict(electron_states=values["electron_states"]),
            source_pins_sha256={str(i):"a"*64 for i in range(13)})
        (location/"manifest.json").write_text(json.dumps(manifest))
        runs[role]=dict(manifest=manifest,summary=dict(path=str(location),energy=dict(relative_limit=.01)),
            history={"sentinel":role},snap={"sentinel":role})
    return decision,runs
def apply(tmp_path,result=None):
    decision,runs=inputs(tmp_path);original=example() if result is None else result
    return dispatch.apply_verified_spectral_policy("spectral",original,runs["reference"],runs["candidate"],decision),runs
def test_actual_manifest_scope_and_original_failure_remain_bound(tmp_path):
    result,runs=apply(tmp_path);before=copy.deepcopy(result)
    assert dispatch.validate_spectral_policy_application(result)["nodal_minimum_over_DC"]==.03
    proof=result["verified_spectral_application_scope"]
    assert proof["cases"]["reference"]["electron_states"]==54
    assert proof["cases"]["candidate"]["electron_states"]==80
    assert result["original_2percent_spectral_result"]["passed"] is False
    assert result["passed"] is True and result["stage5_admission_inferred"] is False
    assert result==before
@pytest.mark.parametrize("fault",["wavelength_source","wavelength_identity","reference_states","candidate_states","cutoff","order","phonons","transport","pins","energy_limit","advanced_manifest"])
def test_policy_does_not_apply_outside_actual_adopted_scope(tmp_path,fault):
    decision,runs=inputs(tmp_path);r=runs["reference"];c=runs["candidate"]
    if fault=="wavelength_source":c["manifest"]["source"]["wavelength_nm"]=1550.
    if fault=="wavelength_identity":c["manifest"]["parameter_identity"]["wavelength_nm"]=1550.
    if fault=="reference_states":r["manifest"]["kinetics"]["electron_states"]=80
    if fault=="candidate_states":c["manifest"]["kinetics"]["electron_states"]=54
    if fault=="cutoff":c["manifest"]["parameter_identity"]["xi_cutoff"]=128.
    if fault=="order":c["manifest"]["parameter_identity"]["xi_order"]=10
    if fault=="phonons":c["manifest"]["parameter_identity"]["phonon_nodes"]=192
    if fault=="transport":c["manifest"]["parameter_identity"]["transport_order"]=128
    if fault=="pins":c["manifest"]["source_pins_sha256"]["0"]="b"*64
    if fault=="energy_limit":c["summary"]["energy"]["relative_limit"]=.02
    if fault=="advanced_manifest":
        (Path(c["summary"]["path"])/"manifest.json").write_text("{}")
    else:
        for run in (r,c):(Path(run["summary"]["path"])/"manifest.json").write_text(json.dumps(run["manifest"]))
    with pytest.raises(ValueError):
        dispatch.apply_verified_spectral_policy("spectral",example(),r,c,decision)
@pytest.mark.parametrize("kind",["mesh","domain","temporal"])
def test_no_other_axis_is_changed(tmp_path,kind):
    decision,runs=inputs(tmp_path)
    with pytest.raises(ValueError):dispatch.apply_verified_spectral_policy(kind,example(),runs["reference"],runs["candidate"],decision)
def test_no_explicit_decision_keeps_existing_behavior_without_file_reads():
    row=example();assert dispatch.apply_verified_spectral_policy("spectral",row,{}, {},None) is row
@pytest.mark.parametrize("fault",["gate","metric","scope","manifest","decision","old_failure"])
def test_display_cannot_accept_tampered_policy_evidence(tmp_path,fault):
    result,runs=apply(tmp_path)
    if fault=="gate":result["checks"]["gap_min_waveform"]=False
    if fault=="metric":result["robust_spatial"]["L2_field"]["maximum_relative_error"]=.0
    if fault=="scope":result["verified_spectral_application_scope"]["cases"]["candidate"]["electron_states"]=54
    if fault=="manifest":(Path(runs["candidate"]["summary"]["path"])/"manifest.json").write_text("{}")
    if fault=="decision":result["spectral_observable_policy"]["decision"]["sha256"]="f"*64
    if fault=="old_failure":result["original_2percent_spectral_result"]["gap_min"]["former_gate_passed"]=True
    with pytest.raises((ValueError,KeyError)):dispatch.validate_spectral_policy_application(result)
def test_finite_window_preserves_raw_transformation_before_manual_completion(tmp_path):
    result,_=apply(tmp_path)
    outer=copy.deepcopy(result);outer["comparison_scope"]="EXPLICIT_FINITE_COMPARISON_WINDOW"
    outer["original_full_run_comparison"]=result;outer["checks"]={"window_physically_covered":True};outer["passed"]=True
    assert dispatch.validate_spectral_policy_application(outer)["central_field_L2_over_DC"]==.02
def test_comparator_proves_originals_before_transforming_clipped_metrics(tmp_path):
    decision,runs=inputs(tmp_path);r=runs["reference"];c=runs["candidate"];row=example()
    control=dict(reference={"segments":[]},candidate={"segments":[]},guard_law={"nominal_anchor_time_ps":22.435})
    row["declared_numerical_control"]=control
    row.update(common_trajectory_interval_ps=[0.,20.],legacy_nodal_policy_checks=copy.deepcopy(row["checks"]))
    order=[]
    def full(*args):order.append("full_original_verified");return copy.deepcopy(row),{},{}
    def baseline(*args):order.append("clipped_metrics");return copy.deepcopy(row),{},{}
    original_apply=dispatch.apply_verified_spectral_policy
    def record(kind,result,ref,cand,path,**kwargs):
        order.append("apply");return original_apply(kind,result,ref,cand,path,**kwargs)
    with patch.object(dispatch.guard,"compare_with_bounded_guard",side_effect=full),patch.object(dispatch.audit,"compare_runs",side_effect=baseline),patch.object(dispatch.guard,"compare_guard_law",return_value={}),patch.object(dispatch,"apply_verified_spectral_policy",side_effect=record):
        comparator,proof=dispatch.prepare_verified_comparator("spectral",r,c,"unused",spectral_policy_decision_path=decision)
        out,_,_=comparator("spectral",dict(r),dict(c))
    assert order==["full_original_verified","apply","clipped_metrics","apply"]
    assert out["original_2percent_spectral_result"]["checks"]["gap_min_waveform"] is False
    assert out["declared_numerical_control"]["full_original_comparison"]==proof
    bad=dict(c);bad["snap"]=copy.deepcopy(c["snap"])
    with pytest.raises(ValueError):comparator("spectral",dict(r),bad)
def patched_dispatch():
    return patch.dict(sys.modules,{"sandbox.stage5_effective_validation.guard_protocol_dispatcher":dispatch})
def test_report_uses_all_four_current_margins_and_keeps_old_failure(tmp_path):
    result,_=apply(tmp_path)
    with patched_dispatch():metric=report.comparison_admission_metrics(result)
    assert metric["spectral_multiple_observables"] is True
    assert metric["primary_gap_margin_ratio"]==pytest.approx(.01889458596621633/.02)
    assert metric["deficit_margin_ratio"]==pytest.approx(.018586696777038093/.02)
    assert metric["nodal_margin_ratio"]==pytest.approx(.02276239242656607/.03)
    assert metric["nodal_former_gate_passed"] is False
    assert metric["Vout_margin_ratio"]==pytest.approx(.000797985789/.02)
def test_caption_is_explicit_about_mixed_spectral_limits_and_temporal_floor():
    text=report.convergence_caption_latex([{"spectral_observable_policy":{}}])
    for expected in [r"mínimo nodal 3\%",r"déficit/pico 2\%","fallo original",r"Sólo en el control temporal","no se suma a la latencia"]:assert expected in text
    assert r"Paso y espectro: 2\%" in report.convergence_caption_latex([])
def evidence_case(tmp_path):
    baseline=example();baseline["checks"]["reference_saved_window_passed"]=False
    baseline["checks"]["candidate_saved_window_passed"]=False;baseline["passed"]=False
    baseline["legacy_nodal_policy_checks"]=copy.deepcopy(baseline["checks"]);baseline["legacy_nodal_policy_passed"]=False
    result,_=apply(tmp_path,baseline);w=window();w["original_comparison_result"]=result
    return w,review()
def validate(w,r):
    with patched_dispatch():
        return windows.validate_window_evidence(w,HASH,r,verify_evidence=lambda e:None,primary_summary=primary())
def test_manual_covered_window_remains_required_and_preserves_original_failures(tmp_path):
    w,r=evidence_case(tmp_path);out=validate(w,r)
    assert out["passed"] is True and out["original_full_run_comparison"]["passed"] is False
    assert out["original_2percent_spectral_result"]["checks"]["gap_min_waveform"] is False
    assert out["latency_certified_in_this_window"] is False
    r["manual_scientific_review"]=False
    with pytest.raises(ValueError):validate(w,r)
def bounds():
    return {name:dict(maximum_relative_error=.01,declared_relative_limit=margin,
        expression="fixture quantified",physical_argument="fixture only",evidence=[ENTRY])
        for name,margin in [("Vout",.02),("gap_min",.03),("central_gap_L2",.02),("central_deficit",.02)]}
def test_quantified_later_bound_needs_every_current_spectral_observable(tmp_path):
    w,r=evidence_case(tmp_path);q=r["windows"][0]
    q.update(decision="QUANTIFIED_LATER_BOUND",critical_interval_end_ps=25,later_observable_bounds=bounds())
    assert validate(w,r)["passed"] is True
    for missing in ["central_gap_L2","central_deficit","gap_min","Vout"]:
        saved=q["later_observable_bounds"].pop(missing)
        with pytest.raises(ValueError):validate(w,r)
        q["later_observable_bounds"][missing]=saved
@pytest.mark.parametrize("observable",["central_gap_L2","central_deficit","gap_min","Vout"])
def test_later_bound_cannot_relax_any_current_margin(tmp_path,observable):
    w,r=evidence_case(tmp_path);q=r["windows"][0]
    q.update(decision="QUANTIFIED_LATER_BOUND",critical_interval_end_ps=25,later_observable_bounds=bounds())
    bound=q["later_observable_bounds"][observable];bound["declared_relative_limit"]*=2
    with pytest.raises(ValueError):validate(w,r)
def test_stage5_false_still_blocks_before_policy_or_output(tmp_path):
    path=tmp_path/"not_admitted.json";path.write_text('{"stage5":false}')
    output=tmp_path/"forbidden"
    with pytest.raises(ValueError,match="stage5: true"):report.main(["--admission",str(path),"--build-root",str(output),"--compile"])
    assert not output.exists()
