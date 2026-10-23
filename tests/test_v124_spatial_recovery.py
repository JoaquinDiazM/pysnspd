"""Reflecting source, immutable latency reuse, and two concurrent resource leases."""
import importlib.util,json
from pathlib import Path
from types import SimpleNamespace
import h5py
import numpy as np
import pytest
from scipy.integrate import simpson
from pysnspd.experimental.reflected_photon import reflecting_shape
from pysnspd.experimental import latency_continuation as H

ROOT=Path(__file__).resolve().parents[1]
def load(name,rel):
    s=importlib.util.spec_from_file_location(name,ROOT/rel);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m
S=load('_v124_campaign_tests','pipelines/19_v124_campaign.py')

def transverse(y0,sigma=7.549834435,n=16001):
    y=np.linspace(-40,40,n)*1e-9
    xy=np.column_stack((np.zeros(n),y))
    shape,receipt=reflecting_shape(xy,80e-9,sigma*1e-9,y0*1e-9)
    return y,shape,receipt

@pytest.mark.parametrize('y0',[0,20,35,40,-40])
def test_reflection_conserves_continuum_energy_without_extra_photon(y0):
    y,p,r=transverse(y0)
    assert simpson(p,x=y)==pytest.approx(np.sqrt(2*np.pi)*7.549834435e-9,rel=1e-12)
    assert not r['virtual_images_are_additional_photons']

def test_edge_limit_double_half_profile_both_neumann_edges_and_mirror_symmetry():
    y,p,_=transverse(40)
    np.testing.assert_allclose(p,2*np.exp(-.5*((y-40e-9)/7.549834435e-9)**2),atol=1e-22)
    _,left,_=transverse(-40);np.testing.assert_allclose(p,left[::-1],atol=1e-13)
    _,off,_=transverse(35);_,mirror,_=transverse(-35)
    np.testing.assert_allclose(off,mirror[::-1],atol=1e-13)
    # The exact image derivative cancels at both sidewalls; a centered finite
    # difference on the reflected extension tests it without nodal truncation.
    eps=1e-13
    for boundary in [-40e-9,40e-9]:
        xy=np.array([[0,boundary-eps],[0,boundary+eps]])
        value,_=reflecting_shape(xy,80e-9,7.549834435e-9,35e-9)
        assert value[0]==pytest.approx(value[1],rel=1e-12,abs=1e-30)

def resources():return dict(logical_cpus=list(range(32)),physical_core_groups=[[i,i+16] for i in range(16)],
    available_memory_bytes=100*1024**3,total_memory_bytes=128*1024**3,cpu_quota_units=None)

def test_two_full_campaigns_allocate_disjoint_workers_and_two_reserved_cores():
    first=S.allocation(resources(),[],[])
    second=S.allocation(resources(),[],[dict(all_cpus=first['all_cpus'],memory_reservation_bytes=first['memory_reservation_bytes'])])
    assert len(first['slots'])==len(second['slots'])==4
    assert all(len(v)==3 for v in first['slots']+second['slots'])
    assert not set(first['all_cpus'])&set(second['all_cpus'])
    assert len(set(first['all_cpus'])|set(second['all_cpus']))==28
    assert len(second['reserved_whole_core_groups'])==2
    with pytest.raises(ValueError,match='FREE'):
        S.allocation(resources(),[],[dict(all_cpus=first['all_cpus']),dict(all_cpus=second['all_cpus'])])

def test_active_workers_and_memory_reservations_are_preserved():
    active=[dict(pid=123,affinity_cpus=list(range(8))+list(range(16,24)))]
    with pytest.raises(ValueError,match='FREE'):S.allocation(resources(),active,[])
    r=resources();r['available_memory_bytes']=25*1024**3
    with pytest.raises(ValueError,match='RAM'):S.allocation(r,[],[])

def fixture_source(tmp_path):
    source=tmp_path/'base';source.mkdir()
    pins={'native_kernel.py':'unchanged'}
    newpins={**pins,**{k:'new' for k in H.NEW_SOURCES}}
    oldid=dict(current_uA=21.5,retained_fraction=.9,threads=4,physical_cores=2,latency_only=True)
    identity={**oldid,'threads':3,'latency_only':False,'impact_y_nm':0.,'deposited_energy_eV':None,'source_boundary':'preserved'}
    np.savez(source/'dc1d.npz',delta_bar=np.ones(3))
    counts={}
    with h5py.File(source/'data.h5','w',libver='latest') as f:
        for name in ['2d','shadow1d','1d']:
            q=f.create_group(name);n=4 if name=='2d' else 0
            q['snapshot_time_ps']=np.array([0,1,2,3])[:n]
            for key in ['delta_bar','p','n']:q[key]=np.ones((n,3))
            trace=np.zeros(n,dtype=[(k,'f8') for k in ['time_ps','Ib_A','Is_A','vc_V']])
            trace['time_ps']=np.array([0,1,2,3])[:n];trace['Ib_A']=trace['Is_A']=21.5e-6
            q['trace']=trace;counts[name]={k:len(q[k]) for k in q}
    meta=dict(source_pins_sha256=pins,parameter_identity=oldid,phase='2d',time_ps=3.,terminal_diagnostic=False,
        shadow_active=False,archive_counts=counts,source={'retained_energy_J':1e-19,'source_must_not_be_applied_again':True})
    np.savez(source/'checkpoint_latest.npz',metadata_json=np.array(json.dumps(meta)),delta_bar=np.ones(3),
        p=np.ones(3),n=np.ones(3),circuit_SI=np.array([21.5e-6,21.5e-6,0.]))
    m=dict(model_version='1.2.3',status='TRIGGER_MARGIN_COMPLETE',run_complete=True,phase='2d',time_ps=3.,
        source_pins_sha256=pins,parameter_identity=oldid,archive_counts=counts,source=meta['source'],
        observation={'margin_complete':True,'crossing_time_s':1e-12},
        checkpoint_sha256=H.sha(source/'checkpoint_latest.npz'),dc1d_sha256=H.sha(source/'dc1d.npz'))
    (source/'manifest.json').write_text(json.dumps(m))
    return source,newpins,identity

def test_latency_fork_retains_exact_2d_source_and_replays_only_shadow(tmp_path,monkeypatch):
    source,pins,identity=fixture_source(tmp_path);dest=tmp_path/'continuation'
    monkeypatch.setattr(H,'SCRATCH',tmp_path)
    original={p.name:H.sha(p) for p in source.iterdir()}
    H.prepare_latency_continuation(source,dest,pins,identity)
    assert {p.name:H.sha(p) for p in source.iterdir()}==original
    assert H.sha(dest/'data.h5')==original['data.h5']
    m=json.loads((dest/'manifest.json').read_text())
    assert m['model_version']=='1.2.4' and not m['run_complete']
    assert not m['continuation']['photon_redeposited']
    assert not m['continuation']['original_2d_trajectory_reintegrated']
    with np.load(dest/'checkpoint_latest.npz') as z:
        meta=json.loads(str(z['metadata_json']))
        assert meta['replay_saved_overlap'] and not meta['shadow_active']
        assert meta['source']==json.loads((source/'manifest.json').read_text())['source']
        assert meta['time_ps']==3 and meta['parameter_identity']==identity

@pytest.mark.parametrize('damage',['source','parameter','complete','checkpoint','archive'])
def test_mismatched_or_corrupted_base_is_refused(tmp_path,damage):
    source,pins,identity=fixture_source(tmp_path)
    if damage=='source':pins['native_kernel.py']='changed'
    if damage=='parameter':identity['current_uA']=15.5
    if damage=='complete':
        p=source/'manifest.json';m=json.loads(p.read_text());m['run_complete']=False;p.write_text(json.dumps(m))
    if damage=='checkpoint':(source/'checkpoint_latest.npz').write_bytes(b'damaged')
    if damage=='archive':
        with h5py.File(source/'data.h5','r+') as f:f['2d/p'][-1]=2.
    with pytest.raises(ValueError):H.authenticate(source,pins,identity)

def test_commands_preserve_three_threads_and_recovery_budget_eight_hours():
    args=SimpleNamespace(source_runner=ROOT/'pipelines/13_dimensional_photon_run.py',mode='recovery',
        base_campaign=Path('/home/jdiaz/scratch/base'),max_1d_wall_hours=8.)
    cmd=S.command(S.catalog('recovery')[0],Path('/home/jdiaz/scratch/new'),[0,16,1],args)
    assert cmd[cmd.index('--threads')+1]=='3'
    assert cmd[cmd.index('--max-1d-wall-hours')+1]=='8.0'
    assert cmd[cmd.index('--final-time-ns')+1]=='15'
    assert '--recover-from-run' in cmd
    resumed=S.command(S.catalog('recovery')[0],Path('/home/jdiaz/scratch/new'),[0,16,1],args,True)
    assert '--resume' in resumed and '--recover-from-run' not in resumed
    args.mode='spatial';cmd=S.command(S.catalog('spatial')[-1],Path('/home/jdiaz/scratch/new'),[0,16,1],args)
    assert '--latency-only' in cmd and '--deposited-energy-eV' in cmd
    assert cmd[cmd.index('--extra-time-ps')+1]=='2'
