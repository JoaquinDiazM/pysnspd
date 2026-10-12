"""Explicit dispatcher for the bounded guard protocol; baseline auditor unchanged.

The complete ORIGINAL saved runs are verified before any finite-window copy.
Only saved observables are clipped; physical arrays, ledgers, events and terminal
status remain original and visible. This tool never infers stage5 admission.
"""
from __future__ import annotations
import argparse,copy,json,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:sys.path.insert(0,str(ROOT))
from sandbox.stage5_effective_validation import audit_runs as audit
from sandbox.stage5_effective_validation import comparison_window as windows
try:
    from . import bounded_guard_protocol as guard
except ImportError:
    sys.path.insert(0,str(Path(__file__).resolve().parent))
    import bounded_guard_protocol as guard

SPECTRAL_POLICY_SOURCE_SHA="0ceb9ebcf48f3fc685ceb0637b75c66a91fabb677bdf4c7caa1f7970010425c1"
SPECTRAL_CASES={"reference":dict(electron_states=54,xi_cutoff=32.,xi_order=6,phonon_nodes=96,transport_order=48),
    "candidate":dict(electron_states=80,xi_cutoff=64.,xi_order=8,phonon_nodes=144,transport_order=96)}


PHYSICAL_GAP_SOURCE_SHA="8c8c3368b39ade2ed248587c38beff44e11fbfe25dbceb8cedd2bb6ef1e84888"
EXPLORATORY_POLICY_SOURCE_SHA="058bfd619bcd74157a647770e8ad13f2f4e7be3ccf9d272eb338d819ea19cf16"

def reviewed_analysis_module(name,digest):
    import importlib.util
    path=Path(__file__).with_name(name+".py")
    if guard.digest(path)!=digest:raise ValueError("Reviewed analysis source changed: "+name)
    spec=importlib.util.spec_from_file_location("_stage5_reviewed_"+name,path)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    return module

def physical_gap_module():
    return reviewed_analysis_module("physical_gap_comparison",PHYSICAL_GAP_SOURCE_SHA)

def exploratory_policy_module():
    return reviewed_analysis_module("exploratory_spectral_policy",EXPLORATORY_POLICY_SOURCE_SHA)

def correct_physical_gap(result,curves,gapcurves,reference,candidate):
    return physical_gap_module().correct_comparison(result,curves,gapcurves,
        reference,candidate,audit=audit)

def spectral_policy_module():
    from sandbox.stage5_effective_validation import spectral_observable_policy as policy
    if guard.digest(policy.__file__)!=SPECTRAL_POLICY_SOURCE_SHA:
        raise ValueError("Reviewed spectral policy source changed")
    return policy

def spectral_scope_from_manifests(manifests):
    """This adopted policy covers only the actual775nm54->80state pair."""
    pins=None;proof={}
    for role,expected in SPECTRAL_CASES.items():
        manifest=manifests[role];identity=manifest["parameter_identity"]
        wavelength=manifest["source"].get("wavelength_nm")
        if isinstance(wavelength,bool) or wavelength!=775. or identity.get("wavelength_nm")!=775.:
            raise ValueError("Spectral policy scope is775nm only")
        if manifest["kinetics"].get("electron_states")!=expected["electron_states"]:
            raise ValueError("Spectral policy scope requires54->80electron states")
        for key,value in expected.items():
            if key!="electron_states" and (isinstance(identity.get(key),bool) or identity.get(key)!=value):
                raise ValueError("Adopted spectral refinement differs: "+key)
        actual=manifest.get("source_pins_sha256",{})
        if len(actual)!=13 or (pins is not None and pins!=actual):
            raise ValueError("All13frozen physical sources must match")
        pins=actual;proof[role]=dict(wavelength_nm=775.,**expected)
    return proof

def apply_verified_spectral_policy(kind,result,reference,candidate,decision_path,*,exploratory_decision_path=None,exploratory_voltage_decision_path=None):
    if exploratory_voltage_decision_path is not None and (exploratory_decision_path is None or decision_path is None):
        raise ValueError("Voltage budget must preserve both explicit prior spectral policies")
    if exploratory_decision_path is not None and decision_path is None:
        raise ValueError("Exploratory policy must preserve the explicit original3/2/2 decision")
    if decision_path is None:return result
    if kind!="spectral":raise ValueError("A spectral decision cannot change another numerical axis")
    policy=spectral_policy_module();manifests={};entries={}
    for role,run in (("reference",reference),("candidate",candidate)):
        energy=run["summary"]["energy"]
        if isinstance(energy.get("relative_limit"),bool) or energy.get("relative_limit")!=.01:
            raise ValueError("Original full-run energy margin must remain1percent")
        source=guard.under(Path(run["summary"]["path"])/"manifest.json")
        entry=dict(path=str(source),sha256=guard.digest(source))
        saved=guard.read(guard.bound(entry))
        if saved!=run["manifest"]:raise ValueError("Original audited manifest advanced or changed")
        manifests[role]=saved;entries[role]=entry
    scope=spectral_scope_from_manifests(manifests)
    output=policy.apply_spectral_observable_policy(result,decision_path=decision_path)
    output["verified_spectral_application_scope"]=dict(schema="pysnspd.verified_spectral_scope.v1",
        original_manifests=entries,cases=scope,global_energy_relative_limit=.01,
        full_original_comparison_verified_before_policy=True,stage5_admission_inferred=False)
    if exploratory_decision_path is not None:
        output=exploratory_policy_module().apply_exploratory_spectral_budget(
            output,decision_path=exploratory_decision_path)
    if exploratory_voltage_decision_path is not None:
        output=exploratory_policy_module().apply_exploratory_voltage_budget(
            output,decision_path=exploratory_voltage_decision_path)
    return output

def validate_spectral_policy_application(result):
    """Recheck the exact transformation and its SHA-bound actual manifest scope."""
    if result.get("comparison_scope")=="EXPLICIT_FINITE_COMPARISON_WINDOW":
        result=result["original_full_run_comparison"]
    physical_gap_module().validate_physical_gap_correction(result)
    if "exploratory_spectral_voltage_policy" in result:
        original=result["original_2percent_voltage_spectral_result"]
        limits=validate_spectral_policy_application(original)
        meta=result["exploratory_spectral_voltage_policy"]
        decision=guard.bound(meta["decision"])
        expected=exploratory_policy_module().apply_exploratory_voltage_budget(original,decision_path=decision)
        if result!=expected:raise ValueError("Exploratory voltage transformation changed")
        return dict(limits,voltage_waveform_relative_limit=meta["voltage_relative_limit"])
    if "exploratory_spectral_spatial_policy" in result:
        original=result["original_3_2_2_spectral_result"]
        validate_spectral_policy_application(original)
        meta=result["exploratory_spectral_spatial_policy"]
        decision=guard.bound(meta["decision"])
        expected=exploratory_policy_module().apply_exploratory_spectral_budget(original,decision_path=decision)
        if result!=expected:raise ValueError("Exploratory spectral transformation changed")
        return meta["spectral_limits"]
    policy=spectral_policy_module()
    meta=result.get("spectral_observable_policy",{})
    entry=meta.get("decision",{})
    path=guard.bound(entry)
    expected=policy.apply_spectral_observable_policy(result["original_2percent_spectral_result"],decision_path=path)
    for key in ("checks","passed","policy_version","gap_min","relative_limit","robust_spatial",
        "legacy_nodal_policy_checks","legacy_nodal_policy_passed","spectral_observable_policy","original_2percent_spectral_result"):
        if result.get(key)!=expected.get(key):raise ValueError("Spectral policy result changed: "+key)
    proof=result.get("verified_spectral_application_scope",{})
    entries=proof.get("original_manifests",{})
    if set(entries)!={"reference","candidate"}:raise ValueError("Actual original spectral manifests required")
    scope=spectral_scope_from_manifests({role:guard.read(guard.bound(entry)) for role,entry in entries.items()})
    if proof!=dict(schema="pysnspd.verified_spectral_scope.v1",original_manifests=entries,cases=scope,
        global_energy_relative_limit=.01,full_original_comparison_verified_before_policy=True,stage5_admission_inferred=False):
        raise ValueError("Explicit original spectral scope proof changed")
    return meta["spectral_limits"]

RESOURCE_ADAPTER_SOURCE_SHA="283eb0067e881d586751b922fa9cc38f6849b658f5d91174d9d58aef98d3dc12"
CONTROL_RESOURCE_SOURCE_SHA="fc3aae98132d304743897c9151161d77b06b07c0b62b4395f4957b4937118a3c"
def compare_originals(kind,reference,candidate,protocol_path,*,resource_transition_receipt_path=None,control_resource_transition_receipt_path=None):
    if control_resource_transition_receipt_path is not None:
        if resource_transition_receipt_path is not None or kind not in {"domain","spectral"}:raise ValueError("Separate closed domain/spectral resource selector only; cannot combine with fine selector")
        adapter=reviewed_analysis_module("control_resource_transition",CONTROL_RESOURCE_SOURCE_SHA)
        return adapter.compare_with_resource_tail(kind,reference,candidate,protocol_path,audit,resource_receipt_path=control_resource_transition_receipt_path)
    if resource_transition_receipt_path is None:
        return guard.compare_with_bounded_guard(kind,reference,candidate,protocol_path,audit)
    if kind!="mesh":raise ValueError("Explicit resource-only continuation is limited to the fine mesh control")
    adapter=reviewed_analysis_module("fine_resource_guard_adapter",RESOURCE_ADAPTER_SOURCE_SHA)
    return adapter.compare_with_resource_tail(kind,reference,candidate,protocol_path,audit,
        resource_receipt_path=resource_transition_receipt_path)

def prepare_verified_comparator(kind,reference,candidate,protocol_path,*,spectral_policy_decision_path=None,exploratory_spectral_decision_path=None,exploratory_voltage_decision_path=None,resource_transition_receipt_path=None,control_resource_transition_receipt_path=None):
    """Verify originals once; then bind each analysis copy to those same arrays."""
    extra={} if control_resource_transition_receipt_path is None else dict(control_resource_transition_receipt_path=control_resource_transition_receipt_path)
    full,fullcurves,fullgap=compare_originals(kind,reference,candidate,protocol_path,resource_transition_receipt_path=resource_transition_receipt_path,**extra)
    full,fullcurves,fullgap=correct_physical_gap(full,fullcurves,fullgap,reference,candidate)
    full=apply_verified_spectral_policy(kind,full,reference,candidate,spectral_policy_decision_path,
        exploratory_decision_path=exploratory_spectral_decision_path,
        exploratory_voltage_decision_path=exploratory_voltage_decision_path)
    verified=copy.deepcopy(full["declared_numerical_control"])
    def compare_views(requested_kind,reference_view,candidate_view):
        if requested_kind!=kind:raise ValueError("Verified numerical axis changed")
        for original,view in ((reference,reference_view),(candidate,candidate_view)):
            for key in ("manifest","history","snap"):
                if view.get(key) is not original.get(key):
                    raise ValueError("Finite metrics must retain the verified ORIGINAL physical tables and identity")
        result,curves,gapcurves=audit.compare_runs(kind,reference_view,candidate_view)
        baseline=copy.deepcopy(result)
        policy=guard.compare_guard_law(verified["reference"]["segments"],verified["candidate"]["segments"],
            kind,verified["guard_law"]["nominal_anchor_time_ps"],result["common_trajectory_interval_ps"])
        forbidden=copy.deepcopy(result["forbidden_parameter_changes"])
        forbidden.pop("local_energy_tolerance",None)
        result["forbidden_parameter_changes"]=forbidden
        result["checks"]["same_nonrefined_parameters"]=not forbidden
        result["checks"]["verified_bounded_prospective_guard_policy"]=True
        if verified.get("resource_only_continuation_explicit"):
            result["checks"]["verified_explicit_resource_only_continuation"]=True
        result["passed"]=all(result["checks"].values())
        result["legacy_nodal_policy_checks"]["same_nonrefined_parameters"]=not forbidden
        result["legacy_nodal_policy_checks"]["verified_bounded_prospective_guard_policy"]=True
        if verified.get("resource_only_continuation_explicit"):
            result["legacy_nodal_policy_checks"]["verified_explicit_resource_only_continuation"]=True
        result["legacy_nodal_policy_passed"]=all(result["legacy_nodal_policy_checks"].values())
        result["baseline_comparison_without_guard_protocol"]=baseline
        control=copy.deepcopy(verified);control["guard_law"]=policy
        control["verification_scope"]="Complete ORIGINAL full saved runs verified before finite observable clipping"
        control["full_original_comparison"]=copy.deepcopy(full)
        result["declared_numerical_control"]=control
        result,curves,gapcurves=correct_physical_gap(result,curves,gapcurves,reference_view,candidate_view)
        result=apply_verified_spectral_policy(kind,result,reference,candidate,spectral_policy_decision_path,
            exploratory_decision_path=exploratory_spectral_decision_path,
        exploratory_voltage_decision_path=exploratory_voltage_decision_path)
        return result,curves,gapcurves
    return compare_views,full

def evidence_inputs(protocol_path):
    """Flatten SHA-bound configuration references for the window evidence record."""
    result={str(guard.under(protocol_path)):guard.digest(protocol_path)}
    def walk(value):
        if isinstance(value,dict):
            if set(value)=={"path","sha256"}:
                p=guard.bound(value);result[str(p)]=value["sha256"]
            else:
                for item in value.values():walk(item)
        elif isinstance(value,list):
            for item in value:walk(item)
    walk(guard.read(protocol_path))
    return result

def audit_finite_window(kind,reference_dir,candidate_dir,declaration_path,output_file,*,guard_protocol_path,spectral_policy_decision_path=None,exploratory_spectral_decision_path=None,exploratory_voltage_decision_path=None,resource_transition_receipt_path=None,control_resource_transition_receipt_path=None):
    output=guard.under(output_file)
    if output.exists():raise FileExistsError("Preserve existing window evidence")
    dp=guard.under(declaration_path);decl=guard.read(dp)
    files={str(dp):guard.digest(dp),**evidence_inputs(guard_protocol_path)}
    for path in (Path(reference_dir).resolve(),Path(candidate_dir).resolve()):
        for name in ("manifest.json","trajectory.npz","snapshots.npz","checkpoint_latest.npz"):
            p=path/name;files[str(p)]=audit.sha256(p)
    if spectral_policy_decision_path is not None:
        policy=spectral_policy_module();dp_policy,_=policy.load_decision(spectral_policy_decision_path)
        files[str(dp_policy)]=guard.digest(dp_policy);files[str(Path(policy.__file__).resolve())]=guard.digest(policy.__file__)
    physical=physical_gap_module();helper,_=physical.source_integrity(audit)
    for module in (physical,helper):files[str(Path(module.__file__).resolve())]=guard.digest(module.__file__)
    if exploratory_spectral_decision_path is not None:
        policy=exploratory_policy_module();ep,_=policy.load_decision(exploratory_spectral_decision_path)
        files[str(ep)]=guard.digest(ep);files[str(Path(policy.__file__).resolve())]=guard.digest(policy.__file__)
    if exploratory_voltage_decision_path is not None:
        policy=exploratory_policy_module()
        vp=guard.under(exploratory_voltage_decision_path)
        if guard.digest(vp)!=policy.VOLTAGE_DECISION_SHA256:raise ValueError("Voltage decision SHA mismatch")
        files.update(evidence_inputs(vp))
    if resource_transition_receipt_path is not None:
        files.update(evidence_inputs(resource_transition_receipt_path))
        for name in ("fine_resource_guard_adapter.py","fine_resource_transition.py"):
            p=Path(__file__).with_name(name);files[str(p)]=guard.digest(p)
    if control_resource_transition_receipt_path is not None:
        files.update(evidence_inputs(control_resource_transition_receipt_path))
        for name in ("control_resource_transition.py","num8_resource_equivalence.py","fine_resource_guard_adapter.py","fine_resource_transition.py"):
            p=Path(__file__).with_name(name);files[str(p)]=guard.digest(p)
    reference=audit.audit_run("reference",reference_dir)
    candidate=audit.audit_run("candidate",candidate_dir)
    compare,full=prepare_verified_comparator(kind,reference,candidate,guard_protocol_path,
        spectral_policy_decision_path=spectral_policy_decision_path,
        exploratory_spectral_decision_path=exploratory_spectral_decision_path,
        exploratory_voltage_decision_path=exploratory_voltage_decision_path,
        resource_transition_receipt_path=resource_transition_receipt_path,
        control_resource_transition_receipt_path=control_resource_transition_receipt_path)
    out,curves,gapcurves=windows.evaluate_window(kind,reference,candidate,
        start_ps=decl["start_ps"],end_ps=decl["end_ps"],compare_function=compare,declaration=decl)
    if any(audit.sha256(p)!=h for p,h in files.items()):
        raise ValueError("Saved inputs advanced during audit; use a coherent atomic capture")
    out.update(inputs_sha256=files,auditor_source_sha256=audit.sha256(audit.__file__),
        window_source_sha256=audit.sha256(windows.__file__),guard_dispatcher_source_sha256=guard.digest(__file__),
        bounded_guard_source_sha256=guard.digest(guard.__file__),
        original_full_runs_verified_before_window_crop=True,
        physical_gap_interpolation_method=physical.METHOD,
        interpolation="Saved scalar observables only in time; complex spatial P1 before magnitude; no electrical event or winding interpolated")
    output.parent.mkdir(parents=True,exist_ok=True)
    output.write_text(json.dumps(audit.json_ready(out),indent=2,allow_nan=False)+"\n")
    return out,curves,gapcurves

def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--kind",choices=sorted(guard.ROLES),required=True)
    p.add_argument("--reference",type=Path,required=True)
    p.add_argument("--candidate",type=Path,required=True)
    p.add_argument("--guard-protocol",type=Path,required=True)
    p.add_argument("--window-declaration",type=Path)
    p.add_argument("--spectral-policy-decision",type=Path,help="Explicit adopted775nm54->80state observable decision; no other axis")
    p.add_argument("--exploratory-spectral-decision",type=Path,
        help="Explicit adopted775nm54->80 spatial budget; original2percent and3/2/2 retained")
    p.add_argument("--exploratory-voltage-decision",type=Path,
        help="Explicit post-data775nm54->80 Vout5% budget; latency2% and coverage unchanged")
    p.add_argument("--resource-transition-receipt",type=Path,help="Explicit fine Num2-to-Num4 saved transition; mesh axis only")
    p.add_argument("--control-resource-transition-receipt",type=Path,help="Explicit closed domain/spectral2->4 or spectral4->8 saved continuation; never a fine alias")
    p.add_argument("--output-dir",type=Path,required=True)
    a=p.parse_args(argv);output=guard.under(a.output_dir)
    if output.exists():raise FileExistsError("Existing evidence preserved; choose a fresh output directory")
    if (a.spectral_policy_decision is not None or a.exploratory_spectral_decision is not None or a.exploratory_voltage_decision is not None) and a.kind!="spectral":
        raise ValueError("Spectral policy cannot change another axis")
    if a.resource_transition_receipt is not None and a.kind!="mesh":raise ValueError("Resource-only transition applies to the declared fine mesh control only")
    if a.control_resource_transition_receipt is not None and (a.resource_transition_receipt is not None or a.kind not in {"domain","spectral"}):raise ValueError("Closed domain/spectral resource selector cannot be combined or applied to another axis")
    if a.window_declaration:
        # audit_finite_window creates parent only after all verification.
        result,curves,gapcurves=audit_finite_window(a.kind,a.reference,a.candidate,
            a.window_declaration,output/"window.json",guard_protocol_path=a.guard_protocol,
            spectral_policy_decision_path=a.spectral_policy_decision,
            exploratory_spectral_decision_path=a.exploratory_spectral_decision,
            exploratory_voltage_decision_path=a.exploratory_voltage_decision,
            resource_transition_receipt_path=a.resource_transition_receipt,
            control_resource_transition_receipt_path=a.control_resource_transition_receipt)
    else:
        reference=audit.audit_run("reference",a.reference);candidate=audit.audit_run("candidate",a.candidate)
        result,curves,gapcurves=compare_originals(a.kind,reference,candidate,a.guard_protocol,resource_transition_receipt_path=a.resource_transition_receipt,control_resource_transition_receipt_path=a.control_resource_transition_receipt)
        result,curves,gapcurves=correct_physical_gap(result,curves,gapcurves,reference,candidate)
        result=apply_verified_spectral_policy(a.kind,result,reference,candidate,a.spectral_policy_decision,
            exploratory_decision_path=a.exploratory_spectral_decision,
            exploratory_voltage_decision_path=a.exploratory_voltage_decision)
        output.mkdir(parents=True)
        report=dict(schema="pysnspd.bounded_guard_saved_audit.v1",
            stage5_admission_inferred=False,comparisons=[result],
            runs=dict(reference=reference["summary"],candidate=candidate["summary"]),
            auditor_source_sha256=audit.sha256(audit.__file__),
            guard_dispatcher_source_sha256=guard.digest(__file__),bounded_guard_source_sha256=guard.digest(guard.__file__))
        (output/"audit.json").write_text(json.dumps(audit.json_ready(report),indent=2,allow_nan=False)+"\n")
    audit.write_table(output/"comparison_Vout.csv",curves)
    audit.write_table(output/"comparison_gap.csv",gapcurves)
    print(json.dumps(dict(output_dir=str(output),passed=result["passed"],stage5_admission_inferred=False)))
    return result
if __name__=="__main__":main()
