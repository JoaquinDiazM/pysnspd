import copy,importlib.util,sys
from pathlib import Path
import pytest
HERE=Path(__file__).parent
SOURCE=HERE/"guard_protocol_dispatcher.py"
if not SOURCE.exists():SOURCE=HERE.parent/"sandbox/stage5_effective_validation/guard_protocol_dispatcher.py"
ROOT=HERE.parent if (HERE.parent/"pysnspd").exists() else Path.cwd()
if ROOT.exists() and str(ROOT) not in sys.path:sys.path.insert(0,str(ROOT))
s=importlib.util.spec_from_file_location("sandbox.stage5_effective_validation.guard_protocol_dispatcher_test",SOURCE)
d=importlib.util.module_from_spec(s);s.loader.exec_module(d)
T=22.435
def ss(role):
    old=1.25e-6 if role=="temporal" else 1e-5;new=3.125e-7 if role=="temporal" else 2.5e-6
    anchor=T-.1 if role=="temporal" else T
    base=dict(step_tolerance=2e-7,dt_ps=.01 if role=="temporal" else .02,threads=2,energy_guard_enabled=True)
    return [dict(base,start_time_ps=0.,end_time_ps=12.4,origin="legacy",local_energy_tolerance=None),
        dict(base,start_time_ps=12.4,end_time_ps=anchor,origin="v6",local_energy_tolerance=old),
        dict(base,start_time_ps=anchor,end_time_ps=None,origin=d.guard.ORIGIN,
            local_energy_tolerance=new,guard_active_since_time_ps=anchor)]
def setup(monkeypatch):
    # This fixture isolates existing guard plumbing; real physical correction is tested separately.
    monkeypatch.setattr(d,"correct_physical_gap",lambda result,curves,gap,*runs:(result,curves,gap))
    ref=dict(manifest=dict(run_complete=False,status="RUNNING"),history={"original":"ref"},snap={"original":"ref"},summary=dict(final_time_ps=30.))
    cand=dict(manifest=dict(run_complete=False,status="RUNNING"),history={"original":"cand"},snap={"original":"cand"},summary=dict(final_time_ps=29.))
    proof=dict(guard_law=dict(nominal_anchor_time_ps=T),reference=dict(segments=ss("nominal")),
        candidate=dict(segments=ss("temporal")),original_completion_flags_preserved=True)
    full=dict(passed=False,declared_numerical_control=proof)
    result=dict(passed=False,common_trajectory_interval_ps=[0.,20.],
        forbidden_parameter_changes=dict(local_energy_tolerance=dict(reference=2.5e-6,candidate=3.125e-7)),
        checks=dict(same_nonrefined_parameters=False,reference_saved_window_passed=False,
            candidate_saved_window_passed=False,voltage_waveform=True,gap_min_waveform=True),
        legacy_nodal_policy_checks=dict(same_nonrefined_parameters=False,gap_min_waveform=True))
    calls=[]
    def verified(kind,a,b,protocol,audit):
        assert a is ref and b is cand
        assert a["summary"]["final_time_ps"]==30. and b["summary"]["final_time_ps"]==29.
        calls.append("full_original_verification")
        return copy.deepcopy(full),{},{}
    def baseline(kind,a,b):
        calls.append("finite_baseline_comparison");return copy.deepcopy(result),{},{}
    monkeypatch.setattr(d.guard,"compare_with_bounded_guard",verified)
    monkeypatch.setattr(d.audit,"compare_runs",baseline)
    return ref,cand,full,result,calls
def test_complete_originals_are_verified_before_any_window_crop(monkeypatch):
    ref,cand,full,base,calls=setup(monkeypatch)
    compare,got=d.prepare_verified_comparator("temporal",ref,cand,"protocol")
    rv=dict(ref);cv=dict(cand);rv["summary"]=dict(final_time_ps=20.);cv["summary"]=dict(final_time_ps=20.)
    result,_,_=compare("temporal",rv,cv)
    assert calls==["full_original_verification","finite_baseline_comparison"]
    assert result["checks"]["same_nonrefined_parameters"] is True
    assert result["passed"] is False
    assert result["checks"]["reference_saved_window_passed"] is False
    assert ref["summary"]["final_time_ps"]==30. and ref["manifest"]["run_complete"] is False
    assert result["declared_numerical_control"]["guard_law"]["guard_policy_intervals"][-1]["reference_guard"]==1e-5
    assert result["declared_numerical_control"]["full_original_comparison"]==full
    assert result["baseline_comparison_without_guard_protocol"]==base
@pytest.mark.parametrize("key",["manifest","history","snap"])
def test_cropped_table_may_not_replace_original_verified_table(monkeypatch,key):
    ref,cand,_,_,calls=setup(monkeypatch)
    compare,_=d.prepare_verified_comparator("temporal",ref,cand,"protocol")
    rv=dict(ref);rv[key]=copy.deepcopy(ref[key])
    with pytest.raises(ValueError,match="ORIGINAL"):compare("temporal",rv,dict(cand))
    assert calls==["full_original_verification"]
def test_no_other_forbidden_change_is_removed_in_finite_window(monkeypatch):
    ref,cand,_,base,_=setup(monkeypatch)
    base["forbidden_parameter_changes"]["sigma_nm"]=dict(reference=7.5,candidate=9.)
    compare,_=d.prepare_verified_comparator("temporal",ref,cand,"protocol")
    result,_,_=compare("temporal",dict(ref),dict(cand))
    assert result["checks"]["same_nonrefined_parameters"] is False
    assert set(result["forbidden_parameter_changes"])=={"sigma_nm"}
def test_guard_verification_failure_prevents_any_crop_comparison(monkeypatch):
    ref,cand,_,_,calls=setup(monkeypatch)
    def failed(*a):raise ValueError("Global energy failed")
    monkeypatch.setattr(d.guard,"compare_with_bounded_guard",failed)
    with pytest.raises(ValueError):d.prepare_verified_comparator("temporal",ref,cand,"protocol")
    assert calls==[]
def test_verified_axis_is_not_a_general_dispatch_exception(monkeypatch):
    ref,cand,_,_,_=setup(monkeypatch);compare,_=d.prepare_verified_comparator("temporal",ref,cand,"protocol")
    with pytest.raises(ValueError):compare("mesh",dict(ref),dict(cand))
