"""No native evolution: causal model, hook/guard contracts and strict modes."""
from pathlib import Path
from types import SimpleNamespace
import importlib.util,json,copy
import numpy as np
import pytest
ROOT=Path(__file__).resolve().parents[1]
from pysnspd.experimental.fast_latency import TailModel,HotbeltDetector,FastObserver,install,PROTOCOL,ESTIMATED
from pysnspd.experimental.trigger_window import TriggerWindow,TriggerWindowConfig

def load(name,path):
    s=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m
def model():
    artifact = ROOT / 'pysnspd/data/v126_fast_tail.json'
    if not artifact.is_file():
        pytest.skip('The lab calibration is required for this fast-mode test')
    return TailModel(artifact)
def example():return next(r for r in model().data['training_cases'] if r['case']=='I21p5_E1600_y00')
def rows():
    r=example();t=r['hb_ps'];v=r['Vout_HB_mV']/1e3;slope=r['slope_HB_mV_ps']/1e3;dt=.005
    def row(t,v):return dict(time_ps=t,Vout_V=v,Vdev_V=r['Vdev_HB_mV']/1e3,energy_defect_J=0.,absolute_external_J=0.,domain_fraction=0.)
    return row(t-dt,v-slope*dt),row(t,v),row(t+dt,v+slope*dt)
def config():return TriggerWindowConfig(0.,0.,'nativeDC',1e-4,'positive',.5e-12,.05,2e-12,100e-12)

def test_causal_one_step_prediction_and_horizon_rejection():
    m=model();p,h,n=rows();r=example();v=m.predict(21.5,1.6,0.,h,p,n)
    assert v['admitted'] and abs(v['estimated_latency_ps']-r['trigger_ps'])<1.
    assert v['derivative_bracket_ps']==[p['time_ps'],n['time_ps']]
    assert v['estimated_latency_ps']>n['time_ps']+1.
    assert not m.predict(15.5,.1,0.,h,p,n)['admitted']
    assert not m.predict(21.5,1.6,0.,h,p,n,ceiling_ps=20.)['admitted']
    assert not m.predict(21.5,1.6,0.,h,h,n)['admitted']
    assert len(m.data['training_cases'])==129
    assert all(abs(r['energy_eV']-.48)>1e-9 for r in m.data['training_cases'])

def test_hotbelt_is_exact_inherited_graph_criterion():
    legacy=load('_v126_inherited_detector',ROOT/'plot_pipelines/07_case_report.py')
    xy=np.array([[0.,0.],[1.,0.],[0.,1.],[1.,1.],[0.,2.],[1.,2.]])*1e-9
    edges=np.array([[0,2],[2,4],[1,3],[3,5],[0,1],[2,3],[4,5]])
    dc=np.ones(6,dtype=complex);det=HotbeltDetector(dc,edges,xy);rng=np.random.default_rng(7)
    fields=np.where(rng.random((100,6))>.5,.1,.9).astype(complex)
    assert [det(d) for d in fields]==legacy.hotbelts(fields,dc,edges,xy).tolist()

def test_guard_prevents_early_stop_and_checkpoint_rejects_mode_change():
    observer=FastObserver(model());observer.args=SimpleNamespace(current_uA=21.5,deposited_energy_eV=1.6,impact_y_nm=0.,snapshot_every_ps=.1,
        global_energy_tolerance=.01,guard_cut_energy_relative=.02,ceiling_ps=100.,confirmation_ps=.5)
    observer.reference=1.;observer.state=SimpleNamespace(delta_bar=np.zeros(2));observer.detector=lambda d:True
    p,h,n=rows();observer.previous=p;observer.process(h);assert observer.pending
    invalid=dict(n,domain_fraction=.03);observer.process(invalid)
    assert not observer.decision['admitted']
    saved=observer.checkpoint();clone=FastObserver(model());clone.restore(saved)
    assert clone.saw_hb and clone.decision==observer.decision
    with pytest.raises(ValueError):clone.restore(dict(saved,execution_mode='normal'))

def hooked(tmp_path):
    """Drive production hooks with a synthetic, non-evolving native interface."""
    class MemoryArchive:
        def __init__(self,*a,**k):self.file=SimpleNamespace(attrs={});self.file.attrs=Attrs();self.saved=[];self.trace=[]
        def append(self,name,row):self.trace.append(dict(row))
        def snapshot(self,name,s,state,rhs,t):self.saved.append(t)
    class Attrs(dict):
        def modify(self,k,v):self[k]=v
    def raw_write(path,value):Path(path).write_text(json.dumps(value))
    R=SimpleNamespace(ROOT=ROOT,B=SimpleNamespace(TriggerWindow=TriggerWindow),PhaseArchive=MemoryArchive,
        source_pins=lambda:{},monitor_metrics=lambda *a,**k:{},monitor_observe=lambda *a,**k:(None,None),
        OLD=SimpleNamespace(_ORIGINAL_WRITE_JSON=raw_write),run=lambda a:0)
    observer=install(R,model(),'pipelines/53_v126_fast_photon.py')
    observer.args=SimpleNamespace(current_uA=21.5,deposited_energy_eV=1.6,impact_y_nm=0.,snapshot_every_ps=.1,
        global_energy_tolerance=.01,guard_cut_energy_relative=.02,ceiling_ps=100.,confirmation_ps=.5)
    return R,observer

def test_controller_pause_resume_and_estimated_status_never_fakes_confirmation(tmp_path):
    R,o=hooked(tmp_path);c=R.B.TriggerWindow(config());c.observe_accepted(0.,0.)
    state=c.checkpoint();assert state['fast_mode_state']['execution_mode']=='fast'
    restored=R.B.TriggerWindow.from_checkpoint(state)
    with pytest.raises(ValueError):R.B.TriggerWindow.from_checkpoint(TriggerWindow(config()).checkpoint())
    p,h,n=rows();o.decision=model().predict(21.5,1.6,0.,h,p,n)
    restored.observe_accepted(n['time_ps']*1e-12,n['Vout_V'])
    assert restored.terminal and restored.status==ESTIMATED
    assert not restored.summary()['detected'] and restored.summary()['latency_s'] is None
    again=R.B.TriggerWindow.from_checkpoint(restored.checkpoint());assert again.status==ESTIMATED
    out=tmp_path/'manifest.json';R.write_json(out,dict(status=ESTIMATED,error=None,time_ps=n['time_ps'],run_complete=False))
    manifest=json.loads(out.read_text());assert manifest['run_complete'] and manifest['execution_mode']=='fast'
    R.write_json(out,dict(status='INTEGRATION_FAILURE',error='physical guard',time_ps=n['time_ps'],run_complete=False))
    assert not json.loads(out.read_text())['run_complete']

def test_sparse_history_saves_initial_hb_and_terminal_without_dense_maps(tmp_path):
    R,o=hooked(tmp_path);window=R.B.TriggerWindow(config());window.observe_accepted(0.,0.)
    system=SimpleNamespace(kinetic_residual_energy_scale_J=1.)
    geometry=SimpleNamespace(graph=SimpleNamespace(edges=np.array([[0,1]])),coordinates_m=np.array([[0.,0.],[0.,80e-9]]))
    dc=SimpleNamespace(delta_bar=np.ones(2));state=SimpleNamespace(delta_bar=np.ones(2));rhs=object()
    archive=R.PhaseArchive(tmp_path/'unused',{'2d':system},{'2d':geometry},{'2d':dc})
    zero=dict(time_ps=0.,Vout_V=0.,Vdev_V=0.,energy_defect_J=0.,absolute_external_J=0.,domain_fraction=0.)
    archive.append('2d',zero);archive.snapshot('2d',system,state,rhs,0.)
    archive.snapshot('2d',system,state,rhs,.1);assert archive.saved==[0.]
    p,h,n=rows();o.previous=p;state.delta_bar[:]=.1
    R.monitor_metrics(system,geometry,state,rhs);archive.append('2d',h)
    assert archive.saved[-1]==h['time_ps']
    archive.append('2d',n);window.observe_accepted(n['time_ps']*1e-12,n['Vout_V'])
    archive.snapshot('2d',system,state,rhs,n['time_ps'])
    assert len(archive.saved)==3 and archive.saved[-1]==n['time_ps']

def test_real_parser_all324_active_energy_position_inputs_without_evolution(tmp_path):
    model()
    C=load('_v126_cases_cli',ROOT/'pipelines/54_v126_fast_campaign.py');F=load('_v126_entry_cli',ROOT/'pipelines/53_v126_fast_photon.py')
    energies=(.1,.2,.3,.4,.5,.6,.7199082489669693,.84,.96,1.2,1.4398164979339385,1.6)
    a=SimpleNamespace(mode='spatial',source_runner=ROOT/'pipelines/13_dimensional_photon_run.py',native_duration_ps=None)
    for I in (15.5,18.5,21.5):
        for E in energies:
            for y in range(0,41,5):
                case=dict(current_uA=I,deposited_energy_eV=E,impact_y_nm=float(y),wavelength_nm=700 if E>1.5 else (775 if E>.8 or (I==15.5 and E==.6) else 1550))
                parsed=F.validate_arguments(C.command(case,Path('/not-executed'),[0,16,1],a)[3:])
                assert parsed['execution_mode']=='fast' and parsed['length_nm']==300. and parsed['latency_only']
                assert parsed['fast_protocol']==PROTOCOL and parsed['threads']==3
    reference=tmp_path/'normal';reference.mkdir()
    before = b'{"execution_mode": "normal"}\n'
    (reference / 'manifest.json').write_bytes(before)
    cmd=C.command(dict(current_uA=21.5,deposited_energy_eV=1.6,impact_y_nm=0.,wavelength_nm=700),reference,[0,16,1],a,True)
    with pytest.raises(ValueError,match='Normal-to-fast'):F.validate_arguments(cmd[3:])
    assert (reference/'manifest.json').read_bytes()==before

def test_user_verifier_compares_physical_preparation_up_to_global_phase():
    verifier=load('_v126_preparation_verifier',ROOT/'pipelines/58_v126_verify_result.py')
    prior={'delta_bar':np.array([1.+1j,1.-1j]),'p':np.array([.1,.2]),'temperature_K':np.array([.9,.9])}
    same={k:v.copy() for k,v in prior.items()};same['delta_bar']*=np.exp(1j*.4)
    verifier.compare_preparation(same,prior)
    wrong={k:v.copy() for k,v in same.items()};wrong['p'][0]+=.001
    with pytest.raises(AssertionError):verifier.compare_preparation(wrong,prior)
