"""Synthetic bounded-guard contracts: no simulation, source or data mutation."""
import copy,importlib.util,json
from pathlib import Path
import numpy as np
import pytest
p=Path(__file__).with_name("bounded_guard_protocol.py")
if not p.exists():p=Path(__file__).resolve().parents[1]/"sandbox/stage5_effective_validation/bounded_guard_protocol.py"
s=importlib.util.spec_from_file_location("bounded_guard_protocol",p)
m=importlib.util.module_from_spec(s);s.loader.exec_module(m)
T=22.435
def segments(role,anchor=None):
    temporal=role=="temporal";old=m.OLD_TEMPORAL if temporal else m.OLD_NOMINAL
    new=m.NEW_TEMPORAL if temporal else m.NEW_NOMINAL
    anchor=T if anchor is None else anchor
    base=dict(dt_ps=.01 if temporal else .02,step_tolerance=2e-7,threads=2,
        energy_guard_enabled=True,local_energy_tolerance=old,guard_active_since_time_ps=12.4)
    return [dict(base,start_time_ps=0.,end_time_ps=12.4,origin="legacy",local_energy_tolerance=None,energy_guard_enabled=False),
        dict(base,start_time_ps=12.4,end_time_ps=anchor,origin="v6_exact_storage_implementation"),
        dict(base,start_time_ps=anchor,end_time_ps=None,origin=m.ORIGIN,local_energy_tolerance=new,guard_active_since_time_ps=anchor)]
@pytest.mark.parametrize("kind",["temporal","mesh","domain","spectral"])
@pytest.mark.parametrize("offset",[-.25,-.125,0.])
def test_exact_adopted_guard_pair_and_bounded_candidate_timing(kind,offset):
    r=m.compare_guard_law(segments("nominal"),segments(kind,T+offset),kind,T,[0.,30.])
    assert r["passed"] and not r["axis_isolation_claimed"]
    assert r["legacy_prefix_before_both_guards_ps"]==[0.,12.4]
    for row in r["guard_policy_intervals"]:
        assert row["candidate_guard"]<=row["reference_guard"]/(8. if kind=="temporal" else 1.)
@pytest.mark.parametrize("offset",[-.25001,.00001])
@pytest.mark.parametrize("kind",["temporal","mesh","domain","spectral"])
def test_no_candidate_anchor_outside_declared_band(kind,offset):
    with pytest.raises(ValueError,match="outside bounded"):
        m.compare_guard_law(segments("nominal"),segments(kind,T+offset),kind,T,[0.,30.])
@pytest.mark.parametrize("key,value",[
    ("guard_active_since_time_ps",12.4),("local_energy_tolerance",5e-6),
    ("energy_guard_enabled",False),("step_tolerance",2e-6),("dt_ps",.1)])
def test_new_segment_anchor_guard_nonlinear_tolerance_and_cap_strict(key,value):
    ss=segments("mesh");ss[-1][key]=value
    with pytest.raises(ValueError):m.guard_segment_contract(ss,"mesh",T)
def test_temporal_ratio8_is_not_general_permission_for_other_guard_pairs():
    ss=segments("temporal");ss[-1]["local_energy_tolerance"]=1e-7
    with pytest.raises(ValueError):m.compare_guard_law(segments("nominal"),ss,"temporal",T,[0.,30.])
def test_other_numerical_axes_never_use_temporal_guard_exception():
    with pytest.raises(ValueError):m.compare_guard_law(segments("nominal"),segments("temporal"),"mesh",T,[0.,30.])
def test_hidden_intermediate_guard_policy_rejected():
    ss=segments("mesh");ss[-2]["local_energy_tolerance"]=8e-6
    with pytest.raises(ValueError):m.guard_segment_contract(ss,"mesh",T)
def test_second_guard_transition_requires_new_review():
    ss=segments("mesh");ss[-2]["origin"]=m.ORIGIN
    with pytest.raises(ValueError):m.guard_segment_contract(ss,"mesh",T)
def fixture_arrays():
    old=dict(parameter_identity=dict(local_energy_tolerance=1e-5,dt_ps=.02),
        source_pins_sha256={"x":"x"},source={"retained_energy_J":1.},
        initial_total_energy_J=7.,guard_active_since_time_ps=12.4,
        local_energy_reference_energy_J=1.,accepted_steps=3,absolute_external_energy_integral_J=.5,
        sum_abs_local_energy_defect_J=.01,local_energy_rejections=4,
        max_accepted_local_energy_fraction=.01,max_guarded_accepted_local_energy_fraction=.009,
        energy_guard_ledger_origin="historical")
    ledger=dict(guard_active_since_time_ps=T,max_guarded_accepted_local_energy_fraction=0.,
        energy_guard_ledger_origin="historical;newprospective")
    new=copy.deepcopy(old);new["parameter_identity"]["local_energy_tolerance"]=2.5e-6;new.update(ledger)
    arrays=dict(metadata_json=np.array(json.dumps(old)),delta_bar=np.array([1.+2.j]),
        p=np.array([[.2]]),n=np.array([[.4]]),circuit_SI=np.array([1.,2.,3.]),
        trajectory__time_s=np.array([0.,T*1e-12]),
        trajectory__energy_total_J=np.array([7.,8.]),
        snapshots__time_s=np.array([0.,T*1e-12]),snapshots__delta_J=np.array([[1.+1.j],[.1+.2j]]))
    initial={k:v.copy() for k,v in arrays.items()};initial["metadata_json"]=np.array(json.dumps(new))
    live={k:v.copy() for k,v in initial.items()}
    for k in live:
        if k.startswith(("trajectory__","snapshots__")):
            live[k]=np.concatenate((live[k],live[k][-1:]),axis=0)
    lm=copy.deepcopy(new);lm["accepted_steps"]=4;live["metadata_json"]=np.array(json.dumps(lm))
    # Evolving state may change, immutable initial state and saved history cannot.
    live["delta_bar"]*=.8;live["p"]*=.9
    receipt=dict(energy_guard_change=dict(local_energy_tolerance_after=2.5e-6),guard_ledger_after=ledger)
    return arrays,initial,live,receipt
def test_advanced_checkpoint_compares_exact_saved_prefix_not_initial_full_state():
    a,b,c,r=fixture_arrays();m.validate_initial_and_live(a,b,c,r)
@pytest.mark.parametrize("fault",["state_transfer","initial_metadata","history_prefix","snapshot_prefix","live_source","live_global_origin","live_guard_anchor","namespace"])
def test_no_state_ledger_source_or_saved_history_rewrite(fault):
    a,b,c,r=fixture_arrays()
    if fault=="state_transfer":b["delta_bar"][0]+=1e-20+1j
    if fault=="initial_metadata":
        meta=json.loads(str(b["metadata_json"].item()));meta["absolute_external_energy_integral_J"]=0.;b["metadata_json"]=np.array(json.dumps(meta))
    if fault=="history_prefix":c["trajectory__energy_total_J"][0]+=.1
    if fault=="snapshot_prefix":c["snapshots__delta_J"][0,0]+=.1
    if fault.startswith("live_"):
        meta=json.loads(str(c["metadata_json"].item()))
        if fault=="live_source":meta["source"]={"retained_energy_J":2.}
        if fault=="live_global_origin":meta["initial_total_energy_J"]=0.
        if fault=="live_guard_anchor":meta["guard_active_since_time_ps"]=0.
        c["metadata_json"]=np.array(json.dumps(meta))
    if fault=="namespace":c["secret_new_state"]=np.array([1.])
    with pytest.raises(ValueError):m.validate_initial_and_live(a,b,c,r)
def test_root_declaration_and_actual_anchor_contract():
    d=dict(decision="ADOPT_PROSPECTIVE_GUARD_TIGHTENING_WITH_BOUNDED_TRANSITIONS",nominal_guard=2.5e-6,
        temporal_guard=3.125e-7,old_nominal_guard=1e-5,old_temporal_guard=1.25e-6,
        candidate_anchor_band_ps_relative_to_nominal=[-.25,0],candidate_anchors_must_not_be_after_nominal=True,
        global_energy_limit=.01,new_state_correction=False,global_ledger_reset=False,
        retrospective_new_guard_certification=False,admission=False,
        limits_unchanged=dict(temporal_spectral=.02,mesh_domain=.05,integrated_energy=.01))
    a=dict(schema="pysnspd.intentional_stop_final.v1",status="FINAL_COHERENT_INTERRUPTED",
        prospective_decision_sha256=m.DECLARATION_SHA,admission=False,t_star_ps=T,time_ps=T,
        controller_status="ACTIVE",control_band_ps=[T-.25,T])
    assert m.declaration_contract(d,a)==T
    for key,value in [("global_energy_limit",.02),("nominal_guard",1e-5),("global_ledger_reset",True)]:
        bad=copy.deepcopy(d);bad[key]=value
        with pytest.raises(ValueError):m.declaration_contract(bad,a)
    for key,value in [("t_star_ps",None),("controller_status","TRIGGERED"),("control_band_ps",[T-.5,T])]:
        bad=copy.deepcopy(a);bad[key]=value
        with pytest.raises(ValueError):m.declaration_contract(d,bad)
def test_sha_bound_evidence_mandatory(tmp_path,monkeypatch):
    monkeypatch.setattr(m,"SCRATCH",tmp_path)
    p=tmp_path/"a.json";p.write_text("{}")
    assert m.bound(dict(path=str(p),sha256=m.digest(p)))==p
    with pytest.raises(ValueError):m.bound(dict(path=str(p),sha256="0"*64))
    with pytest.raises(ValueError):m.bound(dict(path=str(p)))

def fake_comparison_fixture(tmp_path):
    from types import SimpleNamespace
    ap=dict(dt_ps=.02,local_energy_tolerance=2.5e-6,step_tolerance=2e-7)
    bp=dict(dt_ps=.01,local_energy_tolerance=3.125e-7,step_tolerance=2e-7)
    ref=dict(summary=dict(name="reference",path="reference"),manifest=dict(parameter_identity=ap,
        source_pins_sha256={str(i):str(i) for i in range(13)},run_complete=False,status="RUNNING"))
    cand=dict(summary=dict(name="candidate",path="candidate"),manifest=dict(parameter_identity=bp,
        source_pins_sha256=copy.deepcopy(ref["manifest"]["source_pins_sha256"]),run_complete=False,status="RUNNING"))
    base=dict(checks=dict(reference_saved_window_passed=False,candidate_saved_window_passed=False,
        same_frozen_model_sources=True,same_nonrefined_parameters=False,gap_min_waveform=True),
        forbidden_parameter_changes=dict(local_energy_tolerance=dict(reference=2.5e-6,candidate=3.125e-7)),
        common_trajectory_interval_ps=[0.,30.],passed=False,
        legacy_nodal_policy_checks=dict(same_nonrefined_parameters=False,gap_min_waveform=True),
        legacy_nodal_policy_passed=False)
    p=tmp_path/"auditor.py";p.write_text("# fixture")
    audit=SimpleNamespace(__file__=str(p),compare_runs=lambda *args:(copy.deepcopy(base),{},{}))
    protocol=tmp_path/"protocol.json";protocol.write_text("{}")
    cfg=dict(nominal={},controls=dict(temporal={}))
    return ref,cand,base,audit,protocol,cfg
def test_verified_exception_never_mutates_completion_identity_or_old_baseline(tmp_path,monkeypatch):
    ref,cand,base,audit,protocol,cfg=fake_comparison_fixture(tmp_path)
    before=copy.deepcopy((ref,cand))
    monkeypatch.setattr(m,"SCRATCH",tmp_path);monkeypatch.setattr(m,"AUDITOR_SHA",m.digest(audit.__file__))
    monkeypatch.setattr(m,"load_protocol",lambda p:(cfg,dict(comparison_label="bounded"),{},T,None,None))
    monkeypatch.setattr(m,"verify_run",lambda run,entry,role,ctx:dict(segments=segments(role)))
    result,_,_=m.compare_with_bounded_guard("temporal",ref,cand,protocol,audit)
    assert (ref,cand)==before and result["passed"] is False
    assert result["checks"]["same_nonrefined_parameters"] is True
    assert result["checks"]["reference_saved_window_passed"] is False
    assert result["baseline_comparison_without_guard_protocol"]==base
    assert result["declared_numerical_control"]["original_completion_flags_preserved"] is True
def test_guard_protocol_does_not_remove_other_forbidden_changes(tmp_path,monkeypatch):
    ref,cand,base,audit,protocol,cfg=fake_comparison_fixture(tmp_path)
    base["forbidden_parameter_changes"]["tau_ep_Tc_ps"]=dict(reference=24.7,candidate=5.)
    audit.compare_runs=lambda *args:(copy.deepcopy(base),{}, {})
    monkeypatch.setattr(m,"SCRATCH",tmp_path);monkeypatch.setattr(m,"AUDITOR_SHA",m.digest(audit.__file__))
    monkeypatch.setattr(m,"load_protocol",lambda p:(cfg,dict(comparison_label="bounded"),{},T,None,None))
    monkeypatch.setattr(m,"verify_run",lambda run,entry,role,ctx:dict(segments=segments(role)))
    result,_,_=m.compare_with_bounded_guard("temporal",ref,cand,protocol,audit)
    assert result["checks"]["same_nonrefined_parameters"] is False
    assert set(result["forbidden_parameter_changes"])=={"tau_ep_Tc_ps"}
def test_different13_source_model_still_rejected_with_guard_protocol(tmp_path,monkeypatch):
    ref,cand,base,audit,protocol,cfg=fake_comparison_fixture(tmp_path)
    cand["manifest"]["source_pins_sha256"]["0"]="modified"
    monkeypatch.setattr(m,"SCRATCH",tmp_path);monkeypatch.setattr(m,"AUDITOR_SHA",m.digest(audit.__file__))
    monkeypatch.setattr(m,"load_protocol",lambda p:(cfg,dict(comparison_label="bounded"),{},T,None,None))
    monkeypatch.setattr(m,"verify_run",lambda run,entry,role,ctx:dict(segments=segments(role)))
    with pytest.raises(ValueError,match="same complete13-source"):
        m.compare_with_bounded_guard("temporal",ref,cand,protocol,audit)

@pytest.mark.parametrize("limit,maximum,passed",[(.02,.005,True),(.01,.01001,True),(.01,float("nan"),True),(.01,.005,False)])
def test_global_energy_limit_cannot_be_relaxed_by_caller(limit,maximum,passed):
    run=dict(summary=dict(energy=dict(relative_limit=limit,maximum_relative_defect=maximum),
        checks=dict(energy_relative_defect=passed)))
    with pytest.raises(ValueError,match="global energy"):
        m.verify_run(run,{},"nominal",({},{},{},T,None,None))

def test_short_saved_window_does_not_invent_post_guard_coverage():
    r=m.compare_guard_law(segments("nominal"),segments("mesh",T-.1),"mesh",T,[0.,5.])
    assert r["guard_policy_intervals"]==[]
    assert r["legacy_prefix_before_both_guards_ps"]==[0.,5.]
