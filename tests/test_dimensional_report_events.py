"""Saved-state event selection; no solver evolution and no fitted event times."""
import importlib.util
import json
from pathlib import Path
import h5py
import numpy as np
import pytest
from scipy.spatial import Delaunay

ROOT=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location("event_report",ROOT/"plot_pipelines/08_dimensional_case_report.py")
R=importlib.util.module_from_spec(spec);spec.loader.exec_module(R)
L=R.legacy_module()


def archive(path,shift=0.,events=True):
    xy=np.array([[x,y] for x in (0,1,2) for y in (-1,0,1)],float)*1e-9
    tri=Delaunay(xy).simplices
    edges=np.array(sorted({tuple(sorted((int(a),int(b)))) for t in tri
        for a,b in ((t[0],t[1]),(t[1],t[2]),(t[2],t[0]))}))
    saved=np.array([0,1,2,2.1,2.2,3,5,6])+shift
    d=np.ones((len(saved)+1,len(xy)),complex)
    if events:
        z=xy[:,0]/1e-9+1j*xy[:,1]/1e-9
        d[2:5]=(z-(.73+.17j))/3
        d[5:8]=.1
    # Uncommitted data must never change the last available time/event selection.
    d[-1]=.01
    trace=np.zeros(7,dtype=[('time_ps','f8'),('Vout_V','f8')])
    trace['time_ps']=np.arange(7)+shift
    if events:trace['Vout_V']=[0,0,0,0,.05e-3,.15e-3,.2e-3]
    with h5py.File(path,'w') as f:
        q=f.create_group('2d')
        for key,value in dict(coordinates_m=xy,triangles=tri,edges=edges,
            dc_delta_bar=np.ones(len(xy)),snapshot_time_ps=np.r_[saved,99+shift],
            delta_bar=d,trace=trace).items():q[key]=value
    return dict(archive_counts={'2d':{'snapshot_time_ps':len(saved),'trace':len(trace)}},
        observation={'config':{'baseline_V':0.},'detected':events})


def plan(path,m,requested=None):
    with h5py.File(path,'r') as f:return R.event_snapshot_plan(f['2d'],m,requested,L)


def test_default_uses_physical_detectors_and_exactly_nine_semantic_slots(tmp_path):
    p=tmp_path/'run.h5';m=archive(p);slots,events,policy=plan(p,m)
    assert [v['label'] for v in slots]==['0 ps','1 ps','VS','VS + 0.1 ps','VS + 0.2 ps','HB','HB + 2 ps','Trigger','Final 2D']
    assert [v['requested_ps'] for v in slots]==pytest.approx([0,1,2,2.1,2.2,3,5,4.5,6])
    assert events['VS']['positive']==1 and events['VS']['negative']==0
    assert events['HB']['time_ps']==3
    assert events['Trigger']['time_ps']==pytest.approx(4.5)
    assert policy['mode']=='per-case events'
    assert slots[-1]['saved_ps']==6 # not uncommitted99
    assert policy['future_visualization_queue'] #4.5ps has no matching map


def test_default_is_recomputed_for_each_case(tmp_path):
    p=tmp_path/'a.h5';a=archive(p)
    q=tmp_path/'b.h5';b=archive(q,shift=10)
    sa,ea,_=plan(p,a);sb,eb,_=plan(q,b)
    assert ea['VS']['time_ps']==2 and eb['VS']['time_ps']==12
    assert sa[3]['requested_ps']!=sb[3]['requested_ps']
    assert sa[5]['requested_ps']!=sb[5]['requested_ps']


def test_missing_events_are_absent_and_use_labelled_final_state(tmp_path):
    p=tmp_path/'run.h5';m=archive(p,events=False)
    slots,events,_=plan(p,m)
    for name in ('VS','HB','Trigger'):assert events[name]['time_ps'] is None
    for j in (2,3,4,5,6,7):
        assert slots[j]['requested_ps'] is None
        assert slots[j]['fallback'] and slots[j]['saved_ps']==6
        assert 'not observed' in R.snapshot_panel_title(slots[j])
    json.dumps(dict(slots=slots,events=events),allow_nan=False)


def test_explicit_requests_preserve_order_pad_to_nine_and_ignore_extras(tmp_path):
    p=tmp_path/'run.h5';m=archive(p)
    slots,_,policy=plan(p,m,[5,.2])
    assert [v['requested_ps'] for v in slots]==[5]+[.2]*8
    assert len(slots)==9 and slots[-1]['repeated_saved_state']
    slots,_,policy=plan(p,m,list(range(9))+[float('nan'),-1])
    assert policy['ignored_extra_requests']==2 and len(slots)==9
    assert slots[-1]['outside_saved_window'] and slots[-1]['saved_ps']==6


def test_offset_past_end_is_flagged_and_earlier_nearest_tie_is_kept(tmp_path):
    p=tmp_path/'run.h5';m=archive(p)
    with h5py.File(p,'r+') as f:
        f['2d/delta_bar'][5]=1 #delayHB to5, soHB+2 exceeds6
    slots,events,_=plan(p,m)
    assert events['HB']['time_ps']==5
    assert slots[6]['requested_ps']==7 and slots[6]['saved_ps']==6
    assert slots[6]['outside_saved_window']
    explicit,_,_=plan(p,m,[.5])
    assert explicit[0]['saved_ps']==0 #tie chooses earlier


def test_one_dimensional_slots_stay_six_and_empty_invalid_inputs_fail():
    assert R.slots([1,2])==[1,2,2,2,2,2]
    assert R.slots([1,2],9)==[1]+[2]*8
    for values in ([],[-1],[float('nan')],[float('inf')]):
        with pytest.raises(ValueError):R.slots(values,9)


def test_threshold_is_positive_baseline_relative_and_not_abs_voltage():
    tr=np.array([(0,-.5e-3),(1,.05e-3),(2,.15e-3)],dtype=[('time_ps','f8'),('Vout_V','f8')])
    assert R.first_voltage_crossing(tr)==pytest.approx(1.5)
    assert R.first_voltage_crossing(tr,baseline_V=.025e-3)==pytest.approx(1.75)
    tr['Vout_V']=-1e-3
    assert R.first_voltage_crossing(tr) is None


def test_batch_preflight_refuses_live_cases_before_creating_output(tmp_path,monkeypatch):
    run=tmp_path/'live';run.mkdir();(run/'manifest.json').write_text('{"status":"RUNNING"}')
    output=tmp_path/'report'
    monkeypatch.setattr(R,'scratch',lambda p:Path(p))
    with pytest.raises(ValueError,match='No reports started'):
        R.main(['--run-dir',str(run),'--output-dir',str(output)])
    assert not output.exists()
