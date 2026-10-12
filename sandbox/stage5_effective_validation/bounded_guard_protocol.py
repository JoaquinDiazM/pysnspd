"""Candidate strict bounded prospective guard protocol; saved-data analysis only.

No model source, state, ledger, event or completion flag is changed. This bridge
requires explicit immutable adoption and anchor records, plus a verified transfer
for each compared run. Numerical observables retain the adopted auditor limits.
"""
from __future__ import annotations
import copy,hashlib,importlib.util,json,math
from pathlib import Path
import numpy as np
SCRATCH=Path("/home/jdiaz/scratch").resolve()
DECLARATION_SHA="4053e1a15560af4a079bae919de86081b9344e4c64492f7ef44755507afd74f5"
ANCHOR_SHA="fbe65cb96af882a82a1649f5dee6cb7dd80584a691abc2b16109df0fa556ec80"
AUDITOR_SHA="3c08c2534efd2ba8d5803842a9e9cf347b6befc46e68ebddc15adb716156c27a"
LINEAGE_SHA="a753432bf4ab76881b586bb7426bb7bcef4df8bd7941d048794146dcbcc19848"
GUARD2B_HELPER_SHA="567d90a84c9c89083b37d759e8a07ee04377c69d57810cc19f402850ddf3de41"
HELPER_APPROVAL_SHA="ca185b92256360d92275681e008fb2e36e518bc796d929cfc3de748aa22c2db5"
ROLES={"temporal","mesh","domain","spectral"}
OLD_NOMINAL=1e-5
OLD_TEMPORAL=1.25e-6
NEW_NOMINAL=2.5e-6
NEW_TEMPORAL=3.125e-7
SCHEMA="pysnspd.bounded_guard_comparison_protocol.v1"
ORIGIN="prospective_local_energy_guard_tightening"

def digest(path):
    h=hashlib.sha256()
    with Path(path).open("rb") as f:
        for block in iter(lambda:f.read(8*1024*1024),b""):h.update(block)
    return h.hexdigest()
def read(path):return json.loads(Path(path).read_text())
def under(path):
    p=Path(path).resolve()
    if p==SCRATCH or not p.is_relative_to(SCRATCH):raise ValueError("Evidence must remain under scratch")
    return p
def bound(entry,expected=None):
    if not isinstance(entry,dict) or set(entry)!={"path","sha256"}:raise ValueError("Explicit path/SHA evidence required")
    p=under(entry["path"]);actual=digest(p)
    if actual!=entry["sha256"] or (expected is not None and actual!=expected):raise ValueError("Evidence SHA mismatch")
    return p
def load_module(path,name,expected):
    if digest(path)!=expected:raise ValueError("Approved analysis dependency changed")
    spec=importlib.util.spec_from_file_location(name,path)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module);return module
def finite_number(x):
    return not isinstance(x,bool) and isinstance(x,(int,float)) and math.isfinite(x)
def close(a,b,label):
    if not finite_number(a) or not finite_number(b):raise ValueError("Nonfinite clock: "+label)
    allowance=64*max(math.ulp(float(a)),math.ulp(float(b)))
    if abs(a-b)>allowance:raise ValueError("Clock mismatch: "+label)
def declaration_contract(d,a):
    if d.get("decision")!="ADOPT_PROSPECTIVE_GUARD_TIGHTENING_WITH_BOUNDED_TRANSITIONS":
        raise ValueError("Prospective decision not adopted")
    values={"nominal_guard":NEW_NOMINAL,"temporal_guard":NEW_TEMPORAL,
        "old_nominal_guard":OLD_NOMINAL,"old_temporal_guard":OLD_TEMPORAL,
        "candidate_anchor_band_ps_relative_to_nominal":[-.25,0],
        "candidate_anchors_must_not_be_after_nominal":True,"global_energy_limit":.01,
        "new_state_correction":False,"global_ledger_reset":False,
        "retrospective_new_guard_certification":False,"admission":False}
    if any(d.get(k)!=v for k,v in values.items()):raise ValueError("Declared guard/physical invariants changed")
    if d.get("limits_unchanged")!={"temporal_spectral":.02,"mesh_domain":.05,"integrated_energy":.01}:
        raise ValueError("Observable tolerances may not change")
    if a.get("schema")!="pysnspd.intentional_stop_final.v1" or a.get("status")!="FINAL_COHERENT_INTERRUPTED":
        raise ValueError("Final accepted nominal anchor required")
    if a.get("prospective_decision_sha256")!=DECLARATION_SHA or a.get("admission") is not False:
        raise ValueError("Anchor does not bind prospective decision")
    t=a.get("t_star_ps")
    if not finite_number(t) or t<=.25 or a.get("time_ps")!=t or a.get("controller_status")!="ACTIVE":
        raise ValueError("Finite pre-trigger accepted anchor required")
    if a.get("control_band_ps")!=[t-.25,t]:raise ValueError("Anchor band differs from declared quarter ps")
    return t

def guard_segment_contract(segments,role,tstar):
    """Validate the declared last guard change, not generic looser tolerances."""
    if role not in ROLES|{"nominal"}:raise ValueError("Unknown physical control role")
    if not segments or segments[-1].get("origin")!=ORIGIN:raise ValueError("Explicit guard-only final segment required")
    final=segments[-1];anchor=final.get("start_time_ps")
    if not finite_number(anchor):raise ValueError("Finite guard activation required")
    close(final.get("guard_active_since_time_ps"),anchor,"prospective segment guard anchor")
    if role=="nominal":close(anchor,tstar,"nominal anchor")
    elif anchor<tstar-.25 or anchor>tstar:raise ValueError("Control activation outside bounded transition band")
    expected=NEW_TEMPORAL if role=="temporal" else NEW_NOMINAL
    old=OLD_TEMPORAL if role=="temporal" else OLD_NOMINAL
    if final.get("local_energy_tolerance")!=expected or final.get("energy_guard_enabled") is not True:
        raise ValueError("New guard differs from adopted role")
    if final.get("step_tolerance")!=2e-7 or final.get("dt_ps")!=(.01 if role=="temporal" else .02):
        raise ValueError("Nonlinear tolerance or step cap changed")
    if len(segments)<2 or segments[-2].get("local_energy_tolerance")!=old:
        raise ValueError("Previous guard differs from adopted legacy role")
    transitions=[s for s in segments if s.get("origin")==ORIGIN]
    if len(transitions)!=1:raise ValueError("Unreviewed extra guard transition")
    for s in segments[:-1]:
        g=s.get("local_energy_tolerance")
        if g is not None and g!=old:raise ValueError("Hidden intermediate guard policy")
    return anchor

def guard_policy_at(segments,time):
    rows=[s for s in segments if s["start_time_ps"]<=time and
        (s["end_time_ps"] is None or time<s["end_time_ps"])]
    if len(rows)!=1:raise ValueError("Lineage does not cover comparison clock")
    return rows[0].get("local_energy_tolerance")

def compare_guard_law(reference,candidate,kind,tstar,interval):
    """Check every constant-policy subinterval; retain earlier legacy scope."""
    if kind not in ROLES:raise ValueError("Unsupported numerical axis")
    ra=guard_segment_contract(reference,"nominal",tstar)
    ca=guard_segment_contract(candidate,kind,tstar)
    start,end=interval
    if not all(finite_number(v) for v in interval) or start<0 or end<=start:
        raise ValueError("Nonempty saved comparison interval required")
    first=lambda ss:min(s["start_time_ps"] for s in ss if s.get("local_energy_tolerance") is not None)
    guarded_start=max(start,first(reference),first(candidate))
    boundaries=[] if guarded_start>=end else sorted({guarded_start,end,*[s["start_time_ps"] for s in reference+candidate
        if guarded_start<s["start_time_ps"]<end]})
    rows=[];ratio=8. if kind=="temporal" else 1.
    for left,right in zip(boundaries[:-1],boundaries[1:]):
        if right<=left:continue
        mid=(left+right)/2;rg=guard_policy_at(reference,mid);cg=guard_policy_at(candidate,mid)
        if rg is None or cg is None:raise ValueError("Missing guard in common guarded interval")
        within=left>=tstar-.25 and right<=tstar
        if cg>rg/ratio:raise ValueError("Candidate guard became looser than declared reference ratio")
        if not within and cg!=rg/ratio:raise ValueError("Different guard outside explicit transition band")
        rows.append(dict(interval_ps=[left,right],reference_guard=rg,candidate_guard=cg,
            reference_to_candidate_ratio=rg/cg,bounded_transition=within))
    return dict(passed=True,nominal_anchor_time_ps=ra,candidate_anchor_time_ps=ca,
        transition_band_ps=[tstar-.25,tstar],guard_policy_intervals=rows,
        legacy_prefix_before_both_guards_ps=[start,min(end,guarded_start)],
        legacy_scope="Earlier unguarded or differently activated prefixes retain their own verified historical controls; no new guard is attributed retrospectively.",
        axis_isolation_claimed=False,stage5_admission_inferred=False)

def same_array(a,b):
    return a.shape==b.shape and a.dtype==b.dtype and np.array_equal(a,b,equal_nan=False)
def validate_initial_and_live(parent,initial,live,receipt):
    """Initial NPZ is immutable; advanced checkpoint only shares exact prefixes."""
    if set(parent)!=set(initial) or set(parent)!=set(live):raise ValueError("Checkpoint array namespace changed")
    old=json.loads(str(parent["metadata_json"].item()))
    new=json.loads(str(initial["metadata_json"].item()))
    expected=copy.deepcopy(old)
    expected["parameter_identity"]["local_energy_tolerance"]=receipt["energy_guard_change"]["local_energy_tolerance_after"]
    expected.update(receipt["guard_ledger_after"])
    if new!=expected:raise ValueError("Initial metadata changed beyond declared prospective guard window")
    for key in parent:
        if key=="metadata_json":continue
        if not same_array(parent[key],initial[key]):raise ValueError("Transfer altered physical/history array: "+key)
        if key.startswith(("trajectory__","snapshots__")):
            v=initial[key];current=live[key]
            if current.ndim!=v.ndim or current.shape[1:]!=v.shape[1:] or current.shape[0]<v.shape[0]:
                raise ValueError("Current saved prefix shape differs: "+key)
            if not same_array(v,current[:len(v)]):raise ValueError("Current saved prefix changed: "+key)
    lm=json.loads(str(live["metadata_json"].item()))
    for key in ("parameter_identity","source_pins_sha256","source","initial_total_energy_J"):
        if lm[key]!=new[key]:raise ValueError("Live identity/source/global energy origin changed: "+key)
    for key in ("guard_active_since_time_ps","local_energy_reference_energy_J"):
        if lm[key]!=new[key]:raise ValueError("Live prospective guard anchor/reference changed")
    return old,new,lm

def load_protocol(path):
    p=under(path);cfg=read(p)
    if cfg.get("schema")!=SCHEMA or cfg.get("admission") is not False:raise ValueError("Explicit non-admitting protocol required")
    d=read(bound(cfg["prospective_decision"],DECLARATION_SHA))
    a=read(bound(cfg["final_nominal_anchor"],ANCHOR_SHA))
    tstar=declaration_contract(d,a)
    helper_path=bound(cfg["approved_guard_helper"])
    if digest(helper_path) not in {d["helper_sha256"],GUARD2B_HELPER_SHA}:
        raise ValueError("Only the explicitly reviewed guard2 or guard2b helper may execute")
    # A successor is accepted only with its exact root-reviewed correction receipt.
    if digest(helper_path)!=d["helper_sha256"]:
        correction=read(bound(cfg["helper_correction_review"],HELPER_APPROVAL_SHA))
        if correction.get("status")!="APPROVED_DESCRIPTIVE_SEGMENT_ANCHOR_CORRECTION":
            raise ValueError("Guard helper successor requires explicit root approval")
        if correction.get("previous_helper_sha256")!=d["helper_sha256"] or correction.get("corrected_helper_sha256")!=digest(helper_path):
            raise ValueError("Guard helper correction binding differs")
        if correction.get("prospective_decision_sha256")!=DECLARATION_SHA or correction.get("admission") is not False:
            raise ValueError("Helper correction lacks prospective non-admitting binding")
        if correction.get("RHS_changed") is not False or correction.get("state_or_ledger_changed") is not False:
            raise ValueError("Helper correction may not alter physics or ledgers")
    helper=load_module(helper_path,"_approved_bounded_guard_helper",digest(helper_path))
    helper.helper_integrity()
    lineage_path=helper_path.with_name("numerical_lineage_v6.py")
    lineage=load_module(lineage_path,"_approved_bounded_guard_lineage",LINEAGE_SHA)
    return cfg,d,a,tstar,helper,lineage

def verify_run(run,entry,role,context):
    cfg,d,a,tstar,helper,lineage=context
    energy=run["summary"].get("energy",{})
    maximum=energy.get("maximum_relative_defect")
    if energy.get("relative_limit")!=.01 or not finite_number(maximum) or maximum>.01:
        raise ValueError("The original global energy bound is exactly1percent")
    if run["summary"]["checks"].get("energy_relative_defect") is not True:
        raise ValueError("Independent global energy audit failed")
    rp=bound(entry["transfer_receipt"]);r,old,parent_manifest,proof=lineage.verify_receipt(rp,run["manifest"]["source_pins_sha256"])
    if r.get("helper_sha256")!=cfg["approved_guard_helper"]["sha256"]:
        raise ValueError("Transfer helper is not the explicitly approved artifact")
    if r.get("action")!="EXPLICIT_LOCAL_ENERGY_GUARD_TIGHTENING_PREPARED_NO_EVOLUTION_NO_SIGNAL":
        raise ValueError("Only explicit prospective guard-only transfer accepted")
    if r.get("pins_before")!=r.get("pins_after") or len(r["pins_after"])!=13 or r.get("RHS_sources_identical") is not True:
        raise ValueError("All13 frozen sources must remain identical")
    if r.get("frontend_changes") or r.get("implementation_changes") or r.get("horizon_change") is not None:
        raise ValueError("Guard transition may not hide other changes")
    if r.get("global_energy_error_before")!=r.get("global_energy_error_after"):
        raise ValueError("Global energy origin/defect/denominator changed")
    if r.get("next_step_ps_before")!=r.get("next_step_ps_after") or r.get("resources_before")!=r.get("resources_after"):
        raise ValueError("Guard-only transition altered next step or resources")
    c=r["energy_guard_change"]
    for key in ("global_error_origin_reset","global_denominator_reset","heat_or_state_correction_applied","prior_steps_certified_under_new_tolerance"):
        if c.get(key) is not False:raise ValueError("Guard transition corrects or recertifies history")
    if c.get("prospective_only") is not True:raise ValueError("Prospective-only declaration required")
    ip=bound(entry["initial_checkpoint"],r["new_checkpoint_bytes_sha256"])
    mp=bound(entry["initial_manifest"]);initial_manifest=read(mp)
    initial=helper.h.decode(ip.read_bytes())
    parent=helper.h.decode((Path(r["parent"]["anchor"])/"checkpoint_original.npz").read_bytes())
    live_payload=(Path(run["summary"]["path"])/"checkpoint_latest.npz").read_bytes()
    live=helper.h.decode(live_payload)
    old,new,lm=validate_initial_and_live(parent,initial,live,r)
    if run["manifest"]["parameter_identity"]!=lm["parameter_identity"] or run["manifest"]["source_pins_sha256"]!=lm["source_pins_sha256"]:
        raise ValueError("Current manifest and checkpoint identity disagree")
    close(run["summary"]["final_time_ps"],lm["time_ps"],"current audited checkpoint")
    for prefix,table in (("trajectory__",run["history"]),("snapshots__",run["snap"])):
        for key,value in live.items():
            if key.startswith(prefix) and (key[len(prefix):] not in table or not same_array(value,table[key[len(prefix):]])):
                raise ValueError("Audit tables and current atomic checkpoint differ: "+key)
    current_guard=helper.validate_old_guard_history(run["manifest"],lm,live)
    if initial_manifest["parameter_identity"]!=new["parameter_identity"] or initial_manifest["source_pins_sha256"]!=new["source_pins_sha256"]:
        raise ValueError("Initial manifest and checkpoint disagree")
    if initial_manifest["physicalparams"]!=parent_manifest["physicalparams"] or run["manifest"]["physicalparams"]!=parent_manifest["physicalparams"]:
        raise ValueError("Guard transfer changed physical parameters")
    for label in ("source","energy_reference"):
        if initial_manifest.get(label)!=parent_manifest.get(label) or run["manifest"].get(label)!=parent_manifest.get(label):
            raise ValueError("Guard transition changed physical "+label)
    historical=helper.validate_old_guard_history(parent_manifest,old,parent)
    if historical!=c["previous_guard_history_validation"]:raise ValueError("Receipt's previous guard audit differs from independent reconstruction")
    if helper.global_error_at_anchor(old,parent)!=r["global_energy_error_before"]:
        raise ValueError("Original global energy ledger proof differs")
    segments,proofs=lineage.resolve_receipt(rp,new["source_pins_sha256"])
    anchor=guard_segment_contract(segments,role,tstar)
    close(c["effective_time_ps"],anchor,"receipt activation")
    close(new["guard_active_since_time_ps"],anchor,"checkpoint activation")
    if role=="nominal" and r["parent"]["checkpoint_bytes_sha256"]!=a["checkpoint_sha256"]:
        raise ValueError("Nominal transfer differs from declared final stopped checkpoint")
    if run["summary"]["checks"].get("independent_local_energy_guard") is not True:
        raise ValueError("Current independent local-energy audit failed")
    helper.h.audit_runtime(Path(run["summary"]["path"])/"runtime_sources",new["source_pins_sha256"])
    return dict(role=role,segments=segments,anchor_time_ps=anchor,receipt_sha256=digest(rp),
        initial_checkpoint_sha256=digest(ip),initial_manifest_sha256=digest(mp),
        current_checkpoint_sha256=hashlib.sha256(live_payload).hexdigest(),
        old_guard_independently_verified=historical,current_guard_independently_verified=current_guard,
        ancestry=proofs,exact_saved_prefix=True,
        global_origin_preserved=True,admission=False)

def compare_with_bounded_guard(kind,reference,candidate,protocol_path,audit):
    """Apply existing observable gates; allow only the SHA-verified guard policy."""
    if digest(audit.__file__)!=AUDITOR_SHA:raise ValueError("Baseline adopted auditor changed")
    context=load_protocol(protocol_path);cfg,d,a,tstar,helper,lineage=context
    if kind not in ROLES or kind not in cfg["controls"]:raise ValueError("Control not declared in protocol")
    refproof=verify_run(reference,cfg["nominal"],"nominal",context)
    candproof=verify_run(candidate,cfg["controls"][kind],kind,context)
    # We never replace run_complete/status or historical parameter identities.
    result,curves,gapcurves=audit.compare_runs(kind,reference,candidate)
    if reference["manifest"]["source_pins_sha256"]!=candidate["manifest"]["source_pins_sha256"]:
        raise ValueError("Numerical controls require the same complete13-source model")
    interval=result["common_trajectory_interval_ps"]
    policy=compare_guard_law(refproof["segments"],candproof["segments"],kind,tstar,interval)
    baseline_result=copy.deepcopy(result)
    forbidden=copy.deepcopy(result["forbidden_parameter_changes"])
    forbidden.pop("local_energy_tolerance",None)
    result["forbidden_parameter_changes"]=forbidden
    result["checks"]["same_nonrefined_parameters"]=not forbidden
    result["checks"]["verified_bounded_prospective_guard_policy"]=True
    result["passed"]=all(result["checks"].values())
    result["legacy_nodal_policy_checks"]["same_nonrefined_parameters"]=not forbidden
    result["legacy_nodal_policy_checks"]["verified_bounded_prospective_guard_policy"]=True
    result["legacy_nodal_policy_passed"]=all(result["legacy_nodal_policy_checks"].values())
    result["baseline_comparison_without_guard_protocol"]=baseline_result
    result["declared_numerical_control"]=dict(schema=SCHEMA,protocol_path=str(under(protocol_path)),
        protocol_sha256=digest(protocol_path),prospective_decision_sha256=DECLARATION_SHA,
        final_nominal_anchor_sha256=ANCHOR_SHA,comparison_label=d["comparison_label"],
        guard_law=policy,reference=refproof,candidate=candproof,
        observable_limits_unchanged=True,global_energy_origin_reset=False,
        original_completion_flags_preserved=True,stage5_admission_inferred=False)
    return result,curves,gapcurves
