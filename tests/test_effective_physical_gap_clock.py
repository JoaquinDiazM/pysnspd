"""Exact SI clock recovery in saved-data analysis; no evolution or rendering."""
import copy
import importlib.util
from pathlib import Path
import numpy as np
import pytest
from sandbox.stage5_effective_validation import audit_runs as audit

ROOT=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location("effective_physical_gap_clock_candidate",
    ROOT/"sandbox/stage5_effective_validation/physical_gap_comparison.py")
physical=importlib.util.module_from_spec(spec);spec.loader.exec_module(physical)


def fixture():
    # These exact adjacent floats reproduce a real23.3ps unit-conversion collapse.
    a=2.3300000000000006e-11
    b=np.nextafter(a,np.inf)
    assert a<b and a*1e12==b*1e12
    first=np.array([0.,a,3e-11])
    second=np.array([0.,b,3e-11])
    grid=audit.common_grid(first,second)
    xy=np.array([[-60.,-40.],[60.,-40.],[60.,40.],[-60.,40.]])*1e-9
    triangle=np.array([[0,1,2],[0,2,3]])
    runs=[]
    for t,phase in ((first,1.),(second,1j)):
        snap=dict(coordinates_m=xy,triangles=triangle,time_s=t,
            delta_J=physical.MEV_J*np.full((len(t),4),phase,dtype=complex))
        runs.append(dict(manifest=dict(source=dict(center_m=[0.,0.])),snap=snap,
            hotbelt=[dict(time_s=float(q),gap_min_meV=1.) for q in t]))
    row=dict(kind="spectral",relative_limit=.02,stage5_admission_inferred=False,
        checks=dict(voltage_waveform=False,gap_min_waveform=False),passed=False,
        robust_spatial=dict(L2_field=dict(maximum_relative_error=1.,DC_gap_scale_meV=1.),
            deficit=dict(relative_error_actual_peak=.01),time_interpolation="original"),
        event=dict(confirmed=True,crossing_time_s=1e-11),original_completion=True)
    # Put the unique original maximum on the SECOND row of the collapsed pair.
    original=np.array([0.,.25,1.,0.])
    gap=dict(time_ps=grid*1e12,field_L2_relative=original,
        reference_gap_min_meV=np.ones(len(grid)))
    return (row,dict(sentinel=np.array([2.])),gap,*runs),grid


def test_exact_SI_recovery_keeps_all_rows_original_max_and_other_failures():
    args,raw=fixture();old=copy.deepcopy(args[0]);oldgap=copy.deepcopy(args[2])
    out,curves,gap=physical.correct_comparison(*args,audit=audit)
    meta=out["physical_gap_interpolation_correction"]
    np.testing.assert_array_equal(meta["comparison_times_s"],raw)
    assert np.all(np.diff(meta["comparison_times_s"])>0)
    assert len(gap["field_L2_relative"])==len(raw)==4
    np.testing.assert_array_equal(gap["time_ps"],oldgap["time_ps"])
    np.testing.assert_array_equal(meta["original_L2_curve"],oldgap["field_L2_relative"])
    np.testing.assert_array_equal(gap["original_complex_time_field_L2_relative"],oldgap["field_L2_relative"])
    assert max(meta["original_L2_curve"])==1.
    assert out["original_complex_time_interpolation_result"]==old and args[0]==old
    assert out["checks"]==old["checks"] and out["passed"] is False
    assert out["event"]==old["event"] and out["original_completion"] is True
    assert np.max(gap["field_L2_relative"])<1e-15
    np.testing.assert_array_equal(curves["sentinel"],[2.])
    physical.validate_physical_gap_correction(out)


@pytest.mark.parametrize("fault",["changed_ps_projection","shortened_ps_grid",
    "decreasing_saved_SI","duplicated_saved_SI","nonfinite_saved_SI",
    "nonfinite_reported_ps","negative_original_L2","nonfinite_original_L2"])
def test_clock_recovery_rejects_manipulated_or_invalid_originals(fault):
    args,_=fixture();row,curves,gap,reference,candidate=args
    if fault=="changed_ps_projection":gap["time_ps"][0]=1e-6
    elif fault=="shortened_ps_grid":
        for key in list(gap):gap[key]=gap[key][:-1]
    elif fault=="decreasing_saved_SI":reference["hotbelt"][1]["time_s"]=4e-11
    elif fault=="duplicated_saved_SI":reference["hotbelt"][1]["time_s"]=0.
    elif fault=="nonfinite_saved_SI":candidate["hotbelt"][1]["time_s"]=np.nan
    elif fault=="nonfinite_reported_ps":gap["time_ps"][0]=np.nan
    elif fault=="negative_original_L2":gap["field_L2_relative"][0]=-.1
    elif fault=="nonfinite_original_L2":gap["field_L2_relative"][0]=np.nan
    with pytest.raises(ValueError):physical.correct_comparison(*args,audit=audit)


def test_unambiguous_report_clock_preserves_previous_operation_exactly():
    args,_=fixture()
    for run in args[3:]:run.pop("hotbelt")
    gap=args[2]
    gap["time_ps"]=np.array([0.,10.,20.,30.])
    out,_,_=physical.correct_comparison(*args,audit=audit)
    np.testing.assert_array_equal(out["physical_gap_interpolation_correction"]["comparison_times_s"],
        gap["time_ps"]*1e-12)
