"""Candidate window-specific guard validation for physical sigma sensitivity.

This is not a numerical control and never changes a simulation parameter. The
reference may already carry a later prospective guard; only its SHA-verified
historical policy in the explicit0--20ps sensitivity window is relevant.
"""
from __future__ import annotations
import math
from pathlib import Path
import numpy as np
try:
    from . import bounded_guard_protocol as guard
except ImportError:
    import bounded_guard_protocol as guard
SOURCE_KEYS=("wavelength_nm","photon_energy_J","retained_fraction","retained_energy_J","center_m","time_origin","initialization")
RATIOS={"narrow":1/math.sqrt(2),"wide":math.sqrt(2)}
SENSITIVITY_SHA="0732584ee5582b7ecb6a22dda66b7050ee32dc8d38278a2c81a88875c46609dc"

def window_guard_law(reference_segments,candidate_segments,interval):
    if interval!=[0.,20.]:raise ValueError("Only the explicitly declared original0--20ps physical sensitivity window is covered")
    boundaries=sorted({0.,20.,*[s["start_time_ps"] for s in reference_segments+candidate_segments
        if 0.<s["start_time_ps"]<20.]})
    rows=[]
    for left,right in zip(boundaries[:-1],boundaries[1:]):
        mid=(left+right)/2
        rg=guard.guard_policy_at(reference_segments,mid)
        cg=guard.guard_policy_at(candidate_segments,mid)
        if cg!=guard.OLD_NOMINAL:raise ValueError("Width cases must retain the original1e-5 guard throughout0--20ps")
        if rg not in (None,guard.OLD_NOMINAL):raise ValueError("Nominal policy changed inside the declared original sensitivity window")
        rows.append(dict(interval_ps=[left,right],nominal_guard=rg,width_case_guard=cg,
            earlier_legacy_nominal_guard_disabled=rg is None,candidate_never_looser=True))
    return dict(passed=True,interval_ps=interval,policy_intervals=rows,
        scope="Physical sigma sensitivity; common legacy nominal controls retained. Width guard starts at0 and is never looser; later guard2b is outside this window.",
        numerical_admission_axis=False,stage5_admission_inferred=False)

def validate_physical_parameters(reference,candidate,label,runtime_keys):
    am,bm=reference["manifest"],candidate["manifest"];ap,bp=am["parameter_identity"],bm["parameter_identity"]
    pins=am.get("source_pins_sha256",{})
    if len(pins)!=13 or pins!=bm.get("source_pins_sha256"):raise ValueError("All13 sources must match")
    if label not in RATIOS:raise ValueError("Expected only narrow or wide")
    changed={k for k in set(ap)|set(bp) if ap.get(k)!=bp.get(k)}
    forbidden=changed-{"sigma_nm","local_energy_tolerance"}-runtime_keys
    if forbidden:raise ValueError("Physical/spectral/solver parameters changed: "+str(sorted(forbidden)))
    if ap.get("local_energy_tolerance")!=guard.NEW_NOMINAL or bp.get("local_energy_tolerance")!=guard.OLD_NOMINAL:
        raise ValueError("This window exception covers only the explicit later nominal guard2b and original width guard")
    if not math.isclose(float(bp["sigma_nm"]),float(ap["sigma_nm"])*RATIOS[label],rel_tol=2e-9,abs_tol=0):
        raise ValueError("Physical sigma factor differs")
    if am["physicalparams"]!=bm["physicalparams"]:raise ValueError("Material/device/circuit changed")
    for key in SOURCE_KEYS:
        if am["source"].get(key)!=bm["source"].get(key):raise ValueError("Physical source changed: "+key)
    if bm["source"].get("sigma_nm")!=bp["sigma_nm"] or bm["source"].get("wavelength_nm")!=775.:
        raise ValueError("Source sigma/wavelength inconsistent")

def verify_width_case(reference,candidate,label,protocol_path,entry,interval,audit):
    """Hash-bound fresh STARTED proof plus exact current history/checkpoint."""
    if guard.digest(audit.__file__)!=guard.AUDITOR_SHA:raise ValueError("Adopted auditor changed")
    context=guard.load_protocol(protocol_path);cfg,d,a,tstar,helper,lineage=context
    refproof=guard.verify_run(reference,cfg["nominal"],"nominal",context)
    validate_physical_parameters(reference,candidate,label,audit.RUNTIME_KEYS)
    run=Path(candidate["summary"]["path"]).resolve()
    log=guard.bound(entry["original_STARTED_log"])
    manifest=candidate["manifest"];arrays=helper.h.decode((run/"checkpoint_latest.npz").read_bytes())
    meta=helper.h.audit_arrays(arrays,manifest)
    if manifest["parameter_identity"]!=meta["parameter_identity"]:raise ValueError("Width manifest/checkpoint identity differs")
    guard.close(candidate["summary"]["final_time_ps"],meta["time_ps"],"width current clock")
    for prefix,table in (("trajectory__",candidate["history"]),("snapshots__",candidate["snap"])):
        for key,value in arrays.items():
            if key.startswith(prefix) and (key[len(prefix):] not in table or not guard.same_array(value,table[key[len(prefix):]])):
                raise ValueError("Width tables/checkpoint incoherent")
    if meta["parameter_identity"].get("step_tolerance")!=2e-7 or meta["parameter_identity"].get("dt_ps")!=.02:
        raise ValueError("Width nonlinear tolerance or step cap changed")
    helper.h.audit_runtime(run/"runtime_sources",meta["source_pins_sha256"])
    energy=candidate["summary"]["energy"]
    if energy.get("relative_limit")!=.01 or not guard.finite_number(energy.get("maximum_relative_defect")) or energy["maximum_relative_defect"]>.01:
        raise ValueError("Original global energy1percent limit required")
    if candidate["summary"]["checks"].get("independent_local_energy_guard") is not True:
        raise ValueError("Width independent local guard failed")
    current_guard=helper.validate_old_guard_history(manifest,meta,arrays)
    proof,payload=helper.w.fresh_evidence(run,arrays,manifest,meta,log,expected_output_dir=entry["original_output_dir"])
    proof["numerical_segments"][0]["duration_ps"]=meta["parameter_identity"].get("duration_ps")
    law=window_guard_law(refproof["segments"],proof["numerical_segments"],interval)
    return dict(passed=True,nominal_transfer=refproof,fresh_width=proof,width_current_guard=current_guard,
        declared_window_policy=law,original_identities_preserved=True,
        original_run_complete=[reference["manifest"].get("run_complete"),manifest.get("run_complete")],
        stage5_admission_inferred=False)
