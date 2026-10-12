"""Explicit saved-comparison correction; baseline comparator remains immutable."""
from __future__ import annotations
import copy,hashlib,importlib.util,math
from pathlib import Path
import numpy as np

BASELINE_SHA="3c08c2534efd2ba8d5803842a9e9cf347b6befc46e68ebddc15adb716156c27a"
HELPER_SHA="4702d5c01a0e61c8776bb1b35e9e32ec4086189a6db430d87273036a18b1e939"
SCHEMA="pysnspd.physical_gap_magnitude_time_correction.v1"
METHOD="Saved complex spatial P1 -> magnitude on each frame -> linear scalar interpolation in time"
MEV_J=1.602176634e-22

def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def helper_module():
    path=Path(__file__).with_name("physical_gap_interpolation.py")
    if sha(path)!=HELPER_SHA:
        raise ValueError("Reviewed physical magnitude helper changed")
    spec=importlib.util.spec_from_file_location("physical_gap_interpolation_bound",path)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    return module

def source_integrity(audit):
    if sha(audit.__file__)!=BASELINE_SHA:
        raise ValueError("Pinned baseline comparator changed")
    helper=helper_module()
    return helper,dict(
        baseline=dict(path=str(Path(audit.__file__).resolve()),sha256=BASELINE_SHA),
        helper=dict(path=str(Path(helper.__file__).resolve()),sha256=HELPER_SHA),
        wrapper=dict(path=str(Path(__file__).resolve()),sha256=sha(__file__)))

def correct_comparison(result,curves,gapcurves,reference,candidate,*,audit):
    """Correct ONLY robust_spatial.L2 and its curve, preserving full originals.

    The caller retains responsibility for authentic original-run/guard/window
    verification. No physical arrays, event, source, coverage, voltage, deficit,
    nodal minimum, numerical margin or completion decision is changed.
    """
    helper,sources=source_integrity(audit)
    if "physical_gap_interpolation_correction" in result:
        raise ValueError("Physical magnitude correction may be applied only once")
    if result.get("kind") not in {"temporal","mesh","domain","spectral"}:
        raise ValueError("Unknown comparison axis")
    if result.get("stage5_admission_inferred") is not False:
        raise ValueError("A numerical comparison cannot admit stage5")
    old=result.get("robust_spatial")
    if not isinstance(old,dict) or "L2_field" not in old:
        raise ValueError("Complete original spatial comparison required")
    times=np.asarray(gapcurves["time_ps"],float)*1e-12
    original=np.asarray(gapcurves["field_L2_relative"],float)
    # The baseline builds its union in SI, then exposes that clock in ps.
    # Distinct SI floats can share one ps float; recover every original row.
    if times.ndim==1 and np.all(np.isfinite(times)) and np.any(np.diff(times)<=0):
        grids=[np.asarray([row["time_s"] for row in run["hotbelt"]],float)
            for run in (reference,candidate)]
        if any(grid.ndim!=1 or len(grid)<2 or not np.all(np.isfinite(grid))
                or np.any(np.diff(grid)<=0) for grid in grids):
            raise ValueError("Original saved SI clocks must be finite and strictly increasing")
        authentic_times=audit.common_grid(*grids)
        if authentic_times.shape!=times.shape or not np.array_equal(
                authentic_times*1e12,np.asarray(gapcurves["time_ps"],float)):
            raise ValueError("Original SI comparison grid does not project exactly to the reported ps rows")
        times=authentic_times
    if (times.ndim!=1 or len(times)<2 or original.shape!=times.shape
        or not np.all(np.isfinite(times)) or np.any(np.diff(times)<=0)
        or not np.all(np.isfinite(original)) or np.any(original<0)):
        raise ValueError("Complete finite original L2 curve on its comparison grid required")
    if float(original.max())!=old["L2_field"]["maximum_relative_error"]:
        raise ValueError("Original curve maximum differs from original comparison")
    samplers=[]
    for run in (reference,candidate):
        center=run["manifest"]["source"]["center_m"][0]
        samplers.append(audit.complex_ROI_sampler(run["snap"],center))
    scale=float(old["L2_field"]["DC_gap_scale_meV"])*MEV_J
    def clock_time(snap,t):
        first=float(snap["time_s"][0]);last=float(snap["time_s"][-1])
        tolerance=32*np.finfo(float).eps*max(abs(first),abs(last),abs(float(t)),1e-12)
        if t<first and first-t<=tolerance:return first
        if t>last and t-last<=tolerance:return last
        return t
    # Match endpoint clock roundoff only; never extrapolate physical observations.
    values=[]
    for t in times:
        a=helper.spatial_magnitude_at(reference["snap"],samplers[0],clock_time(reference["snap"],t))
        b=helper.spatial_magnitude_at(candidate["snap"],samplers[1],clock_time(candidate["snap"],t))
        values.append(float(np.sqrt(np.mean((a-b)**2))/scale))
    new=np.asarray(values)
    output=copy.deepcopy(result);newcurves=copy.deepcopy(curves);newgap=copy.deepcopy(gapcurves)
    output["original_complex_time_interpolation_result"]=copy.deepcopy(result)
    output["robust_spatial"]["L2_field"]["maximum_relative_error"]=float(new.max())
    output["robust_spatial"]["time_interpolation"]=METHOD
    if result["kind"] in {"mesh","domain"}:
        if result.get("relative_limit")!=.05:
            raise ValueError("Spatial numerical margins must remain5percent")
        output["checks"]["central_gap_magnitude_L2"]=bool(new.max()<=.05)
        output["passed"]=all(output["checks"].values())
    newgap["original_complex_time_field_L2_relative"]=original.copy()
    newgap["field_L2_relative"]=new
    output["physical_gap_interpolation_correction"]=dict(
        schema=SCHEMA,method=METHOD,analysis_sources=sources,
        comparison_times_s=times.tolist(),original_L2_curve=original.tolist(),
        endpoint_clock_matching="32eps of saved physical clock; no physical extrapolation",
        corrected_L2_curve=new.tolist(),original_full_comparison_preserved=True,
        no_state_or_gate_or_tolerance_changed=True,no_event_or_winding_interpolated=True,
        changed_active_checks=["central_gap_magnitude_L2"] if result["kind"] in {"mesh","domain"} else [],
        stage5_admission_inferred=False)
    validate_physical_gap_correction(output)
    return output,newcurves,newgap

def underlying_corrected_comparison(result):
    """Normalize known reporting/policy layers without rewriting their originals."""
    for _ in range(6):
        if result.get("comparison_scope")=="EXPLICIT_FINITE_COMPARISON_WINDOW":
            result=result["original_full_run_comparison"]
        elif "original_3_2_2_spectral_result" in result:
            result=result["original_3_2_2_spectral_result"]
        elif "spectral_observable_policy" in result:
            result=result["original_2percent_spectral_result"]
        else:return result
    raise ValueError("Unexpected recursive comparison wrappers")

def validate_physical_gap_correction(result):
    """Validate the exact explicit transformation and its preserved original."""
    row=underlying_corrected_comparison(result)
    meta=row.get("physical_gap_interpolation_correction",{})
    if meta.get("schema")!=SCHEMA or meta.get("method")!=METHOD:
        raise ValueError("Physical magnitude interpolation metadata required")
    helper_module()
    sources=meta.get("analysis_sources",{})
    expected={"baseline":BASELINE_SHA,"helper":HELPER_SHA,"wrapper":sha(__file__)}
    if set(sources)!=set(expected):
        raise ValueError("Complete physical interpolation source bindings required")
    for name,digest in expected.items():
        entry=sources[name]
        if entry.get("sha256")!=digest or sha(entry["path"])!=digest:
            raise ValueError("Physical interpolation analysis dependency changed: "+name)
    original=row["original_complex_time_interpolation_result"]
    old=np.asarray(meta["original_L2_curve"],float)
    new=np.asarray(meta["corrected_L2_curve"],float)
    times=np.asarray(meta["comparison_times_s"],float)
    if (times.ndim!=1 or len(times)<2 or old.shape!=times.shape or new.shape!=times.shape
        or not np.all(np.isfinite(times)) or np.any(np.diff(times)<=0)
        or not np.all(np.isfinite(old)) or not np.all(np.isfinite(new))
        or np.any(old<0) or np.any(new<0)):
        raise ValueError("Invalid preserved original or corrected physical L2 curves")
    if float(old.max())!=original["robust_spatial"]["L2_field"]["maximum_relative_error"]:
        raise ValueError("Original full L2 maximum was altered")
    expected_spatial=copy.deepcopy(original["robust_spatial"])
    expected_spatial["L2_field"]["maximum_relative_error"]=float(new.max())
    expected_spatial["time_interpolation"]=METHOD
    if row["robust_spatial"]!=expected_spatial:
        raise ValueError("Correction may change only physical L2 and time-interpolation method")
    expected_checks=copy.deepcopy(original["checks"])
    changed=[]
    if row["kind"] in {"mesh","domain"}:
        if row["relative_limit"]!=.05:
            raise ValueError("Spatial numerical margin changed")
        expected_checks["central_gap_magnitude_L2"]=bool(new.max()<=.05)
        changed=["central_gap_magnitude_L2"]
    if row["checks"]!=expected_checks or row["passed"]!=all(expected_checks.values()):
        raise ValueError("Physical correction altered an unrelated check")
    for key,value in original.items():
        if key not in {"robust_spatial","checks","passed"} and row.get(key)!=value:
            raise ValueError("Physical correction changed original comparison: "+key)
    if (meta.get("changed_active_checks")!=changed
        or meta.get("original_full_comparison_preserved") is not True
        or meta.get("no_state_or_gate_or_tolerance_changed") is not True
        or meta.get("no_event_or_winding_interpolated") is not True
        or meta.get("stage5_admission_inferred") is not False):
        raise ValueError("Explicit physical correction invariants changed")
    return meta
