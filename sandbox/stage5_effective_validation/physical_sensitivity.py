"""Strict saved-data physical-width sensitivity; never a numerical admission axis."""
from __future__ import annotations
import math
from pathlib import Path
import numpy as np
from sandbox.stage5_effective_validation import audit_runs as audit

FILES={"manifest.json","trajectory.npz","snapshots.npz","checkpoint_latest.npz","source_initial_state.npz"}
CASE_RATIOS={"narrow":1/math.sqrt(2),"wide":math.sqrt(2)}
SOURCE_KEYS=("wavelength_nm","photon_energy_J","retained_fraction","retained_energy_J","center_m","time_origin","initialization")

def exact_saved_index(time_s,time_ps):
    times=np.asarray(time_s,float)*1e12
    tolerance=128*np.finfo(float).eps*max(1.,abs(float(time_ps)))
    found=np.flatnonzero(abs(times-time_ps)<=tolerance)
    if len(found)!=1:raise ValueError("Physical sensitivity maps require one actual saved snapshot at the declared time; no temporal interpolation")
    return int(found[0])

def validate_case_parameters(reference,candidate,label):
    am=reference["manifest"];bm=candidate["manifest"]
    ap=am["parameter_identity"];bp=bm["parameter_identity"]
    pins=am.get("source_pins_sha256",{})
    if len(pins)!=13 or pins!=bm.get("source_pins_sha256"):
        raise ValueError("Physical sensitivity requires identical full13-source pins")
    changed={key for key in set(ap)|set(bp) if ap.get(key)!=bp.get(key)}
    forbidden=changed-{"sigma_nm"}-audit.RUNTIME_KEYS
    if forbidden:raise ValueError("Only physical sigma and explicit runtime identifiers may differ: "+str(sorted(forbidden)))
    if label not in CASE_RATIOS:raise ValueError("Expected narrow and wide physical cases")
    expected=float(ap["sigma_nm"])*CASE_RATIOS[label]
    if not math.isclose(float(bp["sigma_nm"]),expected,rel_tol=2e-9,abs_tol=0):
        raise ValueError("Physical-width case differs from the declared sigma factor")
    if am.get("physicalparams")!=bm.get("physicalparams"):
        raise ValueError("Physical material, circuit and bias must remain identical")
    for key in SOURCE_KEYS:
        if am["source"].get(key)!=bm["source"].get(key):
            raise ValueError("Physical source contract changed: "+key)
    if float(bm["source"]["sigma_nm"])!=float(bp["sigma_nm"]):
        raise ValueError("Source and parameter sigma disagree")
    if float(bm["source"]["wavelength_nm"])!=775.:
        raise ValueError("This physical sensitivity contract is775nm only")

def analyze_sensitivity(config,reference,cases,*,parameter_validator=None):
    if set(cases)!=set(CASE_RATIOS):raise ValueError("Exactly narrow and wide physical cases required")
    interval=config.get("common_interval_ps")
    if not isinstance(interval,list) or len(interval)!=2 or interval[0]!=0 or not all(
        isinstance(v,(float,int)) and not isinstance(v,bool) and math.isfinite(v) for v in interval
    ) or interval[1]<=interval[0]:raise ValueError("Explicit physical common interval starting at0ps required")
    map_time=config.get("map_time_ps")
    if not isinstance(map_time,(int,float)) or isinstance(map_time,bool) or not math.isfinite(map_time) or not interval[0]<=map_time<=interval[1]:
        raise ValueError("Physical map time must lie inside the declared comparison interval")
    validate=validate_case_parameters if parameter_validator is None else parameter_validator
    all_runs={"nominal":reference,**cases};maps={};table=[]
    for label,run in all_runs.items():
        if label!="nominal":validate(reference,run,label)
        summary=run["summary"];manifest=run["manifest"]
        if summary.get("passed") is not True or manifest.get("run_complete") is not True:
            raise ValueError("Each physical sensitivity input must pass its own actual completed saved window")
        trajectory=np.asarray(run["traces"]["time_s"])*1e12
        if trajectory[0]>interval[0] or trajectory[-1]+128*np.finfo(float).eps*max(1.,interval[1])<interval[1]:
            raise ValueError("Physical sensitivity trajectory does not cover declared common interval")
        maps[label]=exact_saved_index(run["snap"]["time_s"],map_time)
        if run.get("spatial") is None:raise ValueError("Physical central deficit needs saved mesh and complex fields")
        st=np.asarray(run["spatial"]["time_s"])*1e12
        if st[0]>interval[0] or st[-1]+128*np.finfo(float).eps*max(1.,interval[1])<interval[1]:
            raise ValueError("Physical deficit snapshots do not cover declared common interval")
        event=summary["event"]
        latency=event.get("latency_s")*1e12 if event.get("confirmed") and event["confirmation_time_s"]*1e12<=interval[1] else None
        initial_Tmax=float(np.max(np.asarray(run["snap"]["Te_K"],float)[0]))
        declared_Tmax=float(manifest["source"]["maximum_initial_temperature_K"])
        if not math.isfinite(initial_Tmax) or initial_Tmax<=0 or not math.isclose(initial_Tmax,declared_Tmax,rel_tol=128*np.finfo(float).eps,abs_tol=0):
            raise ValueError("Initial equivalent temperature must match the actual saved post-cascade snapshot")
        table.append(dict(label=label,maximum_initial_temperature_K=initial_Tmax,sigma_nm=float(manifest["parameter_identity"]["sigma_nm"]),
            retained_energy_J=float(manifest["source"]["retained_energy_J"]),time_origin=manifest["source"]["time_origin"],
            status=summary["status"],run_complete=manifest["run_complete"],saved_final_time_ps=summary["final_time_ps"],
            common_interval_ps=list(interval),map_time_ps=float(np.asarray(run["snap"]["time_s"])[maps[label]]*1e12),
            latency_in_this_window_ps=latency,
            event_scope="Saved confirmation within common window only; no event invented or censored at a development horizon"))
    return dict(schema="pysnspd.physical_width_sensitivity.v1",comparison_scope="EXPLICIT_PHYSICAL_SENSITIVITY",
        numerical_admission_axis=False,experimental_calibration=False,stage5_admission_inferred=False,
        common_interval_ps=list(interval),map_time_ps=map_time,map_indices=maps,cases=table)


def load_sensitivity(config,reference,*,input_hashes,verify_entry,load_run=None,parameter_validator=None):
    """SHA-pin all actual saved arrays and the independent audit; no captures here."""
    if not isinstance(config,dict) or set(config.get("cases",{}))!=set(CASE_RATIOS):
        raise ValueError("Explicit physical_sensitivity_runs cases narrow/wide required")
    if config.get("reference_wavelength_nm")!=775:
        raise ValueError("Physical sensitivity reference must be the admitted775nm primary")
    load_run=audit.audit_run if load_run is None else load_run
    cases={};records=[]
    for label,entry in config["cases"].items():
        root=Path(entry["path"]).expanduser().resolve()
        if not root.is_relative_to(Path("/home/jdiaz/scratch")):
            raise ValueError("Physical sensitivity data must remain in scratch")
        pins=entry.get("sha256",{})
        if set(pins)!=FILES:raise ValueError("Explicit SHA pins for manifest/trajectory/snapshots/checkpoint/initial source required")
        for filename in FILES:
            path,digest=verify_entry({"path":str(root/filename),"sha256":pins[filename]})
            input_hashes[str(path)]=digest
        ae=entry.get("audit")
        if not isinstance(ae,dict) or len(ae.get("sha256",""))!=64:
            raise ValueError("Independent physical sensitivity audit must be SHA-pinned")
        path,digest=verify_entry(ae)
        import json,hashlib
        document=json.loads(path.read_text())
        if document.get("schema")!="pysnspd.effective_saved_data_audit.v2" or document.get("all_requested_checks_passed") is not True or document.get("errors"):
            raise ValueError("Completed physical sensitivity saved-data audit required")
        source_sha=hashlib.sha256(Path(audit.__file__).read_bytes()).hexdigest()
        if document.get("auditor_source_sha256")!=source_sha:
            raise ValueError("Physical sensitivity audit must use the adopted auditor")
        found=[r for r in document.get("runs",{}).values() if Path(r["path"]).resolve()==root]
        if len(found)!=1 or found[0].get("passed") is not True:
            raise ValueError("Audit must identify this exact physical sensitivity run")
        run=load_run(label,root)
        if found[0]["checks"]!=run["summary"]["checks"] or found[0]["final_time_ps"]!=run["summary"]["final_time_ps"]:
            raise ValueError("Saved physical audit differs from independent recomputation")
        for filename in FILES:
            verify_entry({"path":str(root/filename),"sha256":pins[filename]})
        cases[label]=run;input_hashes[str(path)]=digest
        records.append(dict(path=str(path),sha256=digest,scope="Physical sensitivity only; not a numerical convergence control"))
    result=analyze_sensitivity(config,reference,cases,parameter_validator=parameter_validator)
    for entry in config["cases"].values():
        for filename in FILES:verify_entry({"path":str(Path(entry["path"])/filename),"sha256":entry["sha256"][filename]})
    return result,{"nominal":reference,**cases},records
