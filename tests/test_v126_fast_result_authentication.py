"""Reader roundtrip on explicitly synthetic files; no physical system evolves."""
from pathlib import Path
import json
import numpy as np,h5py,pytest
from pysnspd.experimental.fast_latency import TailModel,digest,PROTOCOL,ESTIMATED
from pysnspd.analysis.fast_inventory import validate_run
ROOT=Path(__file__).resolve().parents[1]

def fixture_run(tmp_path):
    artifact = ROOT / 'pysnspd/data/v126_fast_tail.json'
    if not artifact.is_file():
        pytest.skip('The lab calibration is required for this fast-mode test')
    model = TailModel(artifact)
    r = next(row for row in model.data['training_cases'] if row['case'] == 'I21p5_E1600_y00')
    p = dict(model.data['reference_parameter_identity'], current_uA=r['current_uA'],
             deposited_energy_eV=r['energy_eV'], impact_y_nm=r['y_nm'],
             execution_mode='fast', fast_protocol=PROTOCOL, fast_model_sha256=model.sha256)
    pins=dict(model.data['baseline_source_pins'])
    for rel in ('pipelines/18_spatial_photon_run.py','pipelines/35_spatial_energy_entry.py','pipelines/53_v126_fast_photon.py','pysnspd/experimental/fast_latency.py','pysnspd/data/v126_fast_tail.json'):pins[rel]=digest(ROOT/rel)
    t=r['hb_ps'];v=r['Vout_HB_mV']/1000;slope=r['slope_HB_mV_ps']/1000;dt=.005
    rows=[dict(time_ps=0.,Vout_V=0.,Vdev_V=0.),dict(time_ps=t-dt,Vout_V=v-slope*dt,Vdev_V=r['Vdev_HB_mV']/1000),
          dict(time_ps=t,Vout_V=v,Vdev_V=r['Vdev_HB_mV']/1000),dict(time_ps=t+dt,Vout_V=v+slope*dt,Vdev_V=r['Vdev_HB_mV']/1000)]
    decision=model.predict(21.5,1.6,0.,rows[2],rows[1],rows[3]);assert decision['admitted']
    trace=np.zeros(4,dtype=[(k,'f8') for k in rows[0]])
    for j,row in enumerate(rows):
        for k,value in row.items():trace[k][j]=value
    with h5py.File(tmp_path/'data.h5','w') as h:
        h.attrs['execution_mode']='fast';q=h.create_group('2d');q.create_dataset('trace',data=trace)
        q.create_dataset('snapshot_time_ps',data=[0.,t,t+dt]);q.create_dataset('dc_delta_bar',data=np.ones(2,complex))
        q.create_dataset('coordinates_m',data=[[0.,0.],[0.,80e-9]]);q.create_dataset('edges',data=[[0,1]])
        q.create_dataset('delta_bar',data=np.array([[1.,1.],[.1,.1],[.1,.1]],complex))
    meta=dict(parameter_identity=p,source_pins_sha256=pins,terminal_diagnostic=False,controller=dict(fast_mode_state=dict(execution_mode='fast')))
    np.savez(tmp_path/'checkpoint_latest.npz',metadata_json=np.array(json.dumps(meta)),delta_bar=np.ones(2,complex),p=np.ones((2,3)),n=np.ones((2,2)),circuit_SI=np.ones(3))
    np.savez(tmp_path/'dc1d.npz',delta_bar=np.ones(2))
    m=dict(execution_mode='fast',fast_protocol=PROTOCOL,model_version='1.2.6',run_complete=True,status=ESTIMATED,error=None,parameter_identity=p,
           fast_model_sha256=model.sha256,source_pins_sha256=pins,checkpoint_sha256=digest(tmp_path/'checkpoint_latest.npz'),dc1d_sha256=digest(tmp_path/'dc1d.npz'),
           max_global_energy_fraction=0.,domain_fraction=0.,time_ps=t+dt,archive_counts={'2d':{'trace':4,'snapshot_time_ps':3}},fast_result={'decision':decision})
    (tmp_path/'manifest.json').write_text(json.dumps(m));return m

def test_fast_estimate_independently_replays_native_features_and_hb_map(tmp_path):
    m=fixture_run(tmp_path);row=validate_run(tmp_path,ROOT)
    assert row['estimated'] and not row['censored'] and row['execution_mode']=='fast'
    assert row['tau_d_ps']==m['fast_result']['decision']['estimated_latency_ps']

def test_tampered_prediction_is_not_accepted(tmp_path):
    m=fixture_run(tmp_path);m['fast_result']['decision']['estimated_latency_ps']+=.1
    (tmp_path/'manifest.json').write_text(json.dumps(m))
    with pytest.raises(ValueError,match='reproduced'):validate_run(tmp_path,ROOT)

def test_renamed_normal_mode_is_not_accepted(tmp_path):
    m=fixture_run(tmp_path);m['parameter_identity']['execution_mode']='normal'
    (tmp_path/'manifest.json').write_text(json.dumps(m))
    with pytest.raises(ValueError,match='identity'):validate_run(tmp_path,ROOT)
