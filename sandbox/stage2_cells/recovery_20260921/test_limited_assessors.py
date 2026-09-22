"""Negative provenance/admission guards on tiny nonphysical file fixtures."""
from pathlib import Path
import copy
import json
import sys
import numpy as np
import pytest

ROOT=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(Path(__file__).resolve().parent))
import assess_limited_time as contract
import assess_limited_pair as pair


@pytest.fixture
def fixture_root(tmp_path,monkeypatch):
    criteria=json.loads((ROOT/contract.CRITERIA).read_text())
    files={'pysnspd/experimental/fixture.py':b'# Nonphysical guard fixture\n',
           'sandbox/stage2_cells/run_coupled.py':b'# Nonphysical runner fixture\n',
           contract.CATALOG:b'Nonphysical catalogue hash fixture',
           contract.CRITERIA:json.dumps(criteria).encode(),
           **{p:b'# Nonphysical integration source fixture\n' for p in contract.EXTRA}}
    for name,raw in files.items():
        path=tmp_path/name;path.parent.mkdir(parents=True,exist_ok=True);path.write_bytes(raw)
    monkeypatch.setattr(contract,'ROOT',tmp_path)
    monkeypatch.setattr(contract.original,'ROOT',tmp_path)
    return tmp_path


def make_run(root,name,steps=16,method=contract.LIMITED_METHOD):
    path=root/(name+'.json')
    x=np.array([.2,.8]);w=np.array([.3,.7]);omega=np.array([.1,.7]);cap=np.array([.2,.5])
    times=np.linspace(0.,.1,steps+1)
    initial=np.array([.7,.15,.05,.25,.1,0.,0.,0.,0.,0.])
    states=np.tile(initial,(steps+1,1))
    grids=dict(electron_count=x,electron_weights=w,phonon_energies=omega,phonon_capacities=cap)
    archive=path.with_suffix('.npz');start=path.with_name(path.stem+'_initial.npz')
    np.savez(archive,times=times,states=states,**grids);np.savez(start,state=initial,**grids)
    sources={'pysnspd/experimental/fixture.py','sandbox/stage2_cells/run_coupled.py'}
    if method==contract.LIMITED_METHOD:sources|=contract.EXTRA
    snapshot=dict(amplitudes=[.7],excitation_energy=[.2],electron_energy=[.2],phonon_energy=[.04],
        escape=[0.],input=[0.],electron_to_phonon=[0.],condensate_heat=[0.],transport=0.)
    record=dict(status='COMPLETED_NOT_YET_ADJUDICATED',parameters=dict(case='one',steps=steps,
        duration=.1,method=method,rtol=1e-9,output=str(path),heating=0.,infrared=.005),
        source_hashes={p:contract.sha(root/p) for p in sources},
        catalog_sha256=contract.sha(root/contract.CATALOG),criteria_sha256=contract.sha(root/contract.CRITERIA),
        initial_conditions_sha256=contract.sha(start),trajectory_sha256=contract.sha(archive),
        initial=copy.deepcopy(snapshot),final=copy.deepcopy(snapshot),
        snapshots=[copy.deepcopy(snapshot) for _ in times],times=times.tolist(),
        material={'scope':'NONPHYSICAL TEST FIXTURE'},python='fixture',numpy='fixture',scipy='fixture',
        energy_ledger_scaled_max=0.,instantaneous_residual_max=0.,
        minimum_electron=.05,maximum_electron=.15,minimum_phonon=.1)
    if method==contract.LIMITED_METHOD:
        record['time_integration']=dict(method='SSPRK3 with common donor-availability event-flux factors',
            population_clipping=False,posthoc_energy_repair=False,
            stats=dict(method='experimental conservative-event-limited SSPRK3',
                       status='COMPLETED_REQUIRES_CONVERGENCE_ASSESSMENT',completed_steps=steps))
    path.write_text(json.dumps(record))
    return path


def load(path):return contract.original.load_run(path)
def mutate_record(path,mutator):
    data=json.loads(path.read_text());mutator(data);path.write_text(json.dumps(data))


def test_limited_vs_frozen_rk4_contract_is_allowed(fixture_root):
    a=make_run(fixture_root,'limited');b=make_run(fixture_root,'rk',32,'rk4')
    contract.compare(load(a),load(b))


@pytest.mark.parametrize('method',['rk4','dop853'])
def test_only_declared_frozen_reference_methods_are_allowed(fixture_root,method):
    a=make_run(fixture_root,'limited');b=make_run(fixture_root,'reference',32,method)
    contract.compare(load(a),load(b))
    mutate_record(b,lambda r:r['parameters'].update(method='unrecorded_integrator'))
    with pytest.raises(ValueError,match='method'):contract.compare(load(a),load(b))


@pytest.mark.parametrize('physical_change',[
    lambda r:r['parameters'].update(heating=.1),
    lambda r:r['parameters'].update(infrared=.01),
    lambda r:r['material'].update(scope='OTHER PHYSICAL INPUT'),
    lambda r:r.update(numpy='other-version')])
def test_changed_physical_or_dependency_contract_is_rejected(fixture_root,physical_change):
    a=make_run(fixture_root,'a');b=make_run(fixture_root,'b',32)
    mutate_record(a,physical_change)
    with pytest.raises(ValueError):contract.compare(load(a),load(b))


@pytest.mark.parametrize('missing',['pysnspd/experimental/fixture.py',*sorted(contract.EXTRA)])
def test_even_matching_incomplete_source_tables_are_rejected(fixture_root,missing):
    paths=[make_run(fixture_root,'a'),make_run(fixture_root,'b',32)]
    for p in paths:mutate_record(p,lambda r:r['source_hashes'].pop(missing))
    with pytest.raises(ValueError,match='complete frozen source table'):contract.compare(*map(load,paths))


def test_matching_unknown_current_extra_source_is_rejected(fixture_root):
    extra='not_an_allowed_adapter.py';(fixture_root/extra).write_text('# fixture')
    paths=[make_run(fixture_root,'a'),make_run(fixture_root,'b',32)]
    for p in paths:mutate_record(p,lambda r:r['source_hashes'].update({extra:contract.sha(fixture_root/extra)}))
    with pytest.raises(ValueError,match='two declared'):contract.compare(*map(load,paths))


@pytest.mark.parametrize('which',sorted(contract.EXTRA))
def test_each_allowed_extra_must_have_current_hash(fixture_root,which):
    paths=[make_run(fixture_root,'a'),make_run(fixture_root,'b',32)]
    for p in paths:mutate_record(p,lambda r:r['source_hashes'].update({which:'0'*64}))
    with pytest.raises(ValueError,match='current source'):contract.compare(*map(load,paths))


@pytest.mark.parametrize('key',['catalog_sha256','criteria_sha256'])
def test_matching_stale_input_hashes_are_rejected(fixture_root,key):
    paths=[make_run(fixture_root,'a'),make_run(fixture_root,'b',32)]
    for p in paths:mutate_record(p,lambda r:r.update({key:'0'*64}))
    with pytest.raises(ValueError,match=key):contract.compare(*map(load,paths))


@pytest.mark.parametrize('change',[
    lambda r:r['time_integration'].update(population_clipping=True),
    lambda r:r['time_integration'].update(posthoc_energy_repair=True),
    lambda r:r['time_integration'].update(method='RK4 hidden as SSPRK3'),
    lambda r:r['time_integration']['stats'].update(method='Other SSP integrator'),
    lambda r:r['time_integration']['stats'].update(status='IN_PROGRESS_NOT_NUMERICALLY_ADMITTED'),
    lambda r:r['time_integration']['stats'].update(completed_steps=15)])
def test_false_integration_metadata_is_rejected(fixture_root,change):
    a=make_run(fixture_root,'a');b=make_run(fixture_root,'b',32)
    mutate_record(a,change)
    with pytest.raises(ValueError,match='integration source contract'):contract.compare(load(a),load(b))


def test_different_initial_populations_are_rejected(fixture_root):
    a=load(make_run(fixture_root,'a'));b=load(make_run(fixture_root,'b',32))
    a['states'][0,1]+=.001
    with pytest.raises(ValueError,match='initial populations'):contract.compare(a,b)


def test_changed_initial_archive_even_with_new_hash_is_rejected(fixture_root):
    a=make_run(fixture_root,'a');start=a.with_name(a.stem+'_initial.npz')
    with np.load(start) as data:arrays={k:data[k].copy() for k in data.files}
    arrays['state'][1]+=.001;np.savez(start,**arrays)
    mutate_record(a,lambda r:r.update(initial_conditions_sha256=contract.sha(start)))
    with pytest.raises(ValueError,match='declared initial state'):load(a)


def test_changed_trajectory_archive_hash_is_rejected(fixture_root):
    a=make_run(fixture_root,'a');archive=a.with_suffix('.npz')
    archive.write_bytes(archive.read_bytes()+b'changed')
    with pytest.raises(ValueError,match='trajectory NPZ'):load(a)


def invoke_pair(monkeypatch,a,b,output):
    monkeypatch.setattr(sys,'argv',['assess',str(a),str(b),'--output',str(output)])
    pair.main()


def test_pair_16_32_is_only_supplemental_not_three_level_admission(fixture_root,monkeypatch):
    a=make_run(fixture_root,'a',16);b=make_run(fixture_root,'b',32)
    output=fixture_root/'pair.json';invoke_pair(monkeypatch,a,b,output)
    record=json.loads(output.read_text())
    assert record['status']=='PASS_PAIR_DIAGNOSTIC'
    assert record['final_stage2_admission'] is False
    assert record['three_level_time_admission'] is False
    assert len(record['checkpoints'])==17
    assert record['tolerance']==2.5e-5


@pytest.mark.parametrize('steps',[16,24,64])
def test_pair_requires_exact_step_halving(fixture_root,monkeypatch,steps):
    a=make_run(fixture_root,'a',16);b=make_run(fixture_root,'b',steps)
    with pytest.raises(ValueError,match='twice'):invoke_pair(monkeypatch,a,b,fixture_root/'bad.json')


def test_pair_cannot_substitute_frozen_rk4_for_actual_limited_method(fixture_root,monkeypatch):
    a=make_run(fixture_root,'a',16);b=make_run(fixture_root,'b',32,'rk4')
    with pytest.raises(ValueError,match='same limited SSP'):invoke_pair(monkeypatch,a,b,fixture_root/'bad.json')


def test_three_resolution_time_wrapper_preserves_original_guards_and_restores_hook(fixture_root,monkeypatch):
    runs=[make_run(fixture_root,'run'+str(n),n) for n in [4,8,16]]
    reference=make_run(fixture_root,'reference',32,'dop853')
    output=fixture_root/'time.json'
    monkeypatch.setattr(sys,'argv',['assess',*map(str,runs),'--reference',str(reference),'--output',str(output)])
    original_hook=contract.original.require_same_problem
    contract.main()
    record=json.loads(output.read_text())
    assert record['status']=='PASS'
    assert record['integration_contract_comparison']['final_stage2_admission'] is False
    assert contract.original.require_same_problem is original_hook


def test_output_evidence_is_never_overwritten(fixture_root,monkeypatch):
    a=make_run(fixture_root,'a');b=make_run(fixture_root,'b',32)
    output=fixture_root/'exists.json';output.write_text('original')
    with pytest.raises(ValueError,match='overwrite'):invoke_pair(monkeypatch,a,b,output)
    assert output.read_text()=='original'
