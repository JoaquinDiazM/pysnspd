"""Focused display/wiring tests for scratch candidates; never renders or evolves."""
import copy,hashlib,importlib.util,json,sys
from pathlib import Path
from unittest.mock import patch
import pytest
REPO=Path(__file__).resolve().parents[1]
BASE=REPO/"sandbox/stage5_effective_validation"
sys.path.insert(0,str(REPO))
def load(name,path):
    spec=importlib.util.spec_from_file_location(name,path)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    return module
report=load("candidate_guard_timeline_report",BASE/"build_final_report.py")
ps=load("candidate_guard_callback_sensitivity",BASE/"physical_sensitivity.py")
def segments(anchor,new,legacy_anchor,old):
    return [dict(start_time_ps=0.,end_time_ps=10.4,local_energy_tolerance=None),
        dict(start_time_ps=10.4,end_time_ps=legacy_anchor,local_energy_tolerance=None),
        dict(start_time_ps=legacy_anchor,end_time_ps=anchor,local_energy_tolerance=old),
        dict(start_time_ps=anchor,end_time_ps=None,local_energy_tolerance=new)]
def comparison():
    return dict(kind="temporal",reference="nominal",candidate="temporal",
        declared_numerical_control=dict(schema="pysnspd.bounded_guard_comparison_protocol.v1",
            reference=dict(segments=segments(22.435,2.5e-6,12.4,1e-5)),
            candidate=dict(segments=segments(22.25,3.125e-7,12.42,1.25e-6)),
            guard_law=dict(transition_band_ps=[22.185,22.435])),
        accepted_step_statistics_after_both_guard_anchors=dict(available=True,interval_ps=[22.435,30.],
            reference=dict(count=10,median_ps=.0025),candidate=dict(count=20,median_ps=.00125)))
def test_timeline_uses_verified_segments_and_keeps_legacy_prefix():
    row=comparison();old=copy.deepcopy(row)
    displayed=report.temporal_refinement_display_rows([row])[0]
    assert displayed["reference_guard"]==2.5e-6 and displayed["candidate_guard"]==3.125e-7
    assert displayed["reference_guard_segments"][0]==dict(start_time_ps=0.,end_time_ps=12.4,local_energy_tolerance=None)
    assert displayed["candidate_guard_segments"][0]["end_time_ps"]==12.42
    assert displayed["post_anchor_interval_ps"]==[22.435,30.]
    assert displayed["reference_accepted_steps"]["median_ps"]==.0025
    assert displayed["candidate_accepted_steps"]["median_ps"]==.00125
    assert row==old
def test_latex_captions_do_not_attribute_guard_retrospectively():
    text=report.guard_timeline_latex([comparison()])
    for expected in ["775 nm","sin guarda","0--12.4","12.4--22.435","desde 22.435",
        r"2.5\times10^{-6}",r"3.125\times10^{-7}","22.185--22.435","sin aislamiento"]:
        assert expected in text
def test_legacy_rows_still_use_old_declared_pair():
    row=comparison();row["declared_numerical_control"]={"reference_local_energy_tolerance":1e-5,"candidate_local_energy_tolerance":1.25e-6}
    out=report.temporal_refinement_display_rows([row])[0]
    assert out["reference_guard"]==1e-5 and out["candidate_guard"]==1.25e-6
    assert out["reference_guard_segments"] is None and report.guard_timeline_latex([row])==""
@pytest.mark.parametrize("fault",["gap","overlap","open_internal","negative","bool","nan","tolerance","missing_zero"])
def test_invalid_timelines_cannot_be_printed_as_verified(fault):
    data=segments(22.435,2.5e-6,12.4,1e-5)
    if fault=="gap":data[1]["start_time_ps"]=10.5
    if fault=="overlap":data[1]["start_time_ps"]=10.3
    if fault=="open_internal":data[1]["end_time_ps"]=None
    if fault=="negative":data[0]["start_time_ps"]=-1.
    if fault=="bool":data[0]["start_time_ps"]=False
    if fault=="nan":data[2]["end_time_ps"]=float("nan")
    if fault=="tolerance":data[-1]["local_energy_tolerance"]=-1.
    if fault=="missing_zero":data[0]["start_time_ps"]=1.
    with pytest.raises(ValueError):report.compact_guard_segments(data)
def test_pending_post_anchor_stats_are_not_fabricated():
    row=comparison();row["accepted_step_statistics_after_both_guard_anchors"]={"available":False}
    assert report.temporal_refinement_display_rows([row])[0]["post_anchor_interval_ps"] is None
def test_stage5_false_blocks_before_protocol_import_or_outputs(tmp_path):
    entry=tmp_path/"not_admitted.json";entry.write_text(json.dumps(dict(stage5=False,
        physical_sensitivity_runs={"bounded_guard_protocol":{"path":"/never/read","sha256":"f"*64}})))
    output=tmp_path/"forbidden"
    with pytest.raises(ValueError,match="stage5: true"):report.main(["--admission",str(entry),"--build-root",str(output),"--compile"])
    assert not output.exists()
def digest(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def configuration(tmp_path):
    protocol=tmp_path/"protocol.json";protocol.write_text("{}")
    cfg=dict(bounded_guard_protocol=dict(path=str(protocol),sha256=digest(protocol)),
        common_interval_ps=[0.,20.],cases={})
    for label in ["narrow","wide"]:
        log=tmp_path/(label+".log");log.write_text("STARTED0")
        cfg["cases"][label]=dict(original_STARTED_log=dict(path=str(log),sha256=digest(log)),
            original_output_dir="/home/jdiaz/scratch/"+label)
    return cfg
def test_wiring_calls_strict_full_original_verifier_before_result(tmp_path):
    from sandbox.stage5_effective_validation import physical_sensitivity as original_ps
    from sandbox.stage5_effective_validation import width_guard_window_candidate as width
    from sandbox.stage5_effective_validation import guard_protocol_dispatcher as dispatch
    cfg=configuration(tmp_path);ref={"id":"original reference"};cases={k:{"id":k} for k in cfg["cases"]};called=[]
    def loader(config,reference,**kwargs):
        assert reference is ref
        for label,run in cases.items():kwargs["parameter_validator"](reference,run,label)
        return {"numerical_admission_axis":False,"stage5_admission_inferred":False},{"nominal":ref,**cases},[]
    def verifier(r,c,label,path,entry,interval,audit):
        assert r is ref and c is cases[label] and interval==[0.,20.]
        called.append(label);return dict(passed=True,original_identities_preserved=True)
    with patch.object(original_ps,"load_sensitivity",side_effect=loader),patch.object(width,"verify_width_case",side_effect=verifier),patch.object(dispatch,"evidence_inputs",return_value={cfg["bounded_guard_protocol"]["path"]:cfg["bounded_guard_protocol"]["sha256"]}):
        hashes={};result,runs,records=report.load_physical_sensitivity(cfg,ref,input_hashes=hashes)
    assert called==["narrow","wide"] and set(result["declared_guard_window_policies"])==set(cases)
    assert result["stage5_admission_inferred"] is False and result["numerical_admission_axis"] is False
    assert runs["nominal"] is ref and records[0]["sha256"]==cfg["bounded_guard_protocol"]["sha256"]
    assert cfg["cases"]["wide"]["original_STARTED_log"]["path"] in hashes
def test_default_sensitivity_does_not_receive_guard_exception():
    from sandbox.stage5_effective_validation import physical_sensitivity as original_ps
    expected=({}, {}, [])
    with patch.object(original_ps,"load_sensitivity",return_value=expected) as loader:
        assert report.load_physical_sensitivity({}, {},input_hashes={})==expected
    assert "parameter_validator" not in loader.call_args.kwargs
@pytest.mark.parametrize("fault",["protocol_sha","missing_log_sha","failed_full_proof"])
def test_wiring_cannot_silently_ignore_guard_proof_faults(tmp_path,fault):
    from sandbox.stage5_effective_validation import physical_sensitivity as original_ps
    from sandbox.stage5_effective_validation import width_guard_window_candidate as width
    from sandbox.stage5_effective_validation import guard_protocol_dispatcher as dispatch
    cfg=configuration(tmp_path)
    if fault=="protocol_sha":cfg["bounded_guard_protocol"]["sha256"]="f"*64
    if fault=="missing_log_sha":cfg["cases"]["narrow"]["original_STARTED_log"].pop("sha256")
    def loader(config,reference,**kwargs):
        kwargs["parameter_validator"](reference,{},"narrow")
        raise AssertionError("Failed proof unexpectedly returned")
    with patch.object(original_ps,"load_sensitivity",side_effect=loader),patch.object(width,"verify_width_case",side_effect=ValueError("Full proof failed")),patch.object(dispatch,"evidence_inputs",return_value={}):
        with pytest.raises(ValueError):report.load_physical_sensitivity(cfg,{},input_hashes={})
def test_callback_is_explicit_and_default_parameter_guard_stays_strict():
    from tests.test_effective_physical_sensitivity import EffectivePhysicalSensitivityTests
    cfg,ref,cases=EffectivePhysicalSensitivityTests().inputs()
    ref["manifest"]["parameter_identity"]["local_energy_tolerance"]=2.5e-6
    with pytest.raises(ValueError,match="Only physical sigma"):ps.analyze_sensitivity(cfg,ref,cases)
    original=copy.deepcopy(ref);called=[]
    result=ps.analyze_sensitivity(cfg,ref,cases,parameter_validator=lambda r,c,label:called.append(label))
    assert called==["narrow","wide"] and result["stage5_admission_inferred"] is False
    assert ref["manifest"]==original["manifest"]
