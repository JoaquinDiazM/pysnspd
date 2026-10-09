"""End-to-end first-step pipeline test on actual small spatial components.

Inputs are synthesized by the component fixture and numeric provenance is
mocked accordingly. The expensive real-strip loader and finite-state
spectral reconstruction are not asserted by this test. Matrix columns,
parallel execution, QR, coupled recomputation, physical-field saves and
checkpoint identity checks use the actual pipeline implementation.
"""
from argparse import Namespace
from pathlib import Path
from tempfile import TemporaryDirectory
import json
import numpy as np
import pytest
from test_joint_rosenbrock import physical_fixture
from sandbox.stage5_photon import run_joint_first_step as runner


@pytest.fixture
def campaign(monkeypatch):
    snapshot,bundle=physical_fixture()
    snapshot.provenance={'pinned_files_sha256':{'synthetic_component_input':'fixture_v1'}}
    monkeypatch.setattr(runner,'build_initial_system',lambda *a,**k:snapshot)
    monkeypatch.setattr(runner,'initial_snapshot_check',lambda *a,**k:
        {'scope':'synthetic component pipeline check only','old_candidate_not_used':True})
    monkeypatch.setattr(runner,'build_thermal_bundle',lambda *a,**k:
        (bundle,{'synthetic_M16_input':'fixture_v1'}))
    parent=Path('/home/jdiaz/scratch/tmp/pysnspd')
    parent.mkdir(parents=True,exist_ok=True)
    with TemporaryDirectory(prefix='test_joint_first_step_',dir=parent) as directory:
        args=Namespace(initial=Path(directory)/'synthetic_initial',snapshot=Path(directory)/'synthetic_snapshot',
            output=Path(directory)/'campaign',dt_ps=[.001,.0005],workers=1,
            execute=True,resume=False,pilot=False,tau_escape_ps=15.)
        runner.run(args)
        yield args,snapshot


def test_actual_parallel_matrix_qr_save_and_two_dt_campaign(campaign):
    args,snapshot=campaign
    summary=json.loads((args.output/'summary.json').read_text())
    assert summary['status']=='COUPLED_FIRST_STEP_CAMPAIGN_COMPLETED_REVIEW_REQUIRED'
    assert summary['accepted_physical_time_steps']==0
    assert summary['stage51_complete'] is False
    assert summary['stage5_complete'] is False
    assert summary['time_admission'] is False
    assert summary['energy_admission'] is False
    assert summary['source_applied_again'] is False
    assert summary['full_reset'] is False
    assert summary['dt_ps']==args.dt_ps
    for index,step in enumerate(summary['steps']):
        folder=args.output/('step_%02d'%index)
        matrix=np.load(folder/'matrix.npy',allow_pickle=False)
        assert matrix.shape==(snapshot.size,snapshot.size)
        assert np.all(np.isfinite(matrix))
        assert len(list((folder/'columns').glob('*.npy')))==snapshot.size
        assert step['recomputed_scaled_rows_max']<1e-7
        assert step['accepted_physical_time_steps']==0
        assert step['general_material_rhs_available'] is False
        assert step['support']['exact_Pauli_support_ok']
        assert step['support']['no_clipping']
        # This legitimate q=0 component fixture has no positive bias. The
        # current-direction guard must prevent detector reconstruction.
        assert step['current_direction_ok'] is False
        assert step['finite_state_rebuild']['status']=='FINITE_STATE_NOT_REBUILT_OUTSIDE_DIAGNOSTIC_SUPPORT_LIMIT_OR_CURRENT_DIRECTION'
        with np.load(folder/'candidate_fields.npz',allow_pickle=False) as fields:
            assert fields['h_new'].shape==(4,snapshot.graph.n_nodes,2)
            assert fields['R0'].shape==(4,snapshot.graph.n_nodes,2,2)
            assert fields['phonon_occupation_new'].shape[1]==snapshot.graph.n_nodes
            assert np.all(fields['phonon_occupation_new']>=0)
            assert np.all(abs(fields['h_new'][...,0]+fields['h_new'][...,1])<=1+1e-14)
            np.testing.assert_allclose(fields['gap_bar_new'],snapshot.gap,atol=1e-12)
            assert float(fields['dt_s'])==pytest.approx(args.dt_ps[index]*1e-12)


def test_resume_rejects_changed_physical_dt(campaign):
    args,_=campaign
    args.resume=True;args.dt_ps=[.002,.001]
    with pytest.raises(ValueError,match='Campaign identity changed: dt_ps'):
        runner.run(args)


def test_resume_rejects_changed_producer_pins(campaign,monkeypatch):
    args,_=campaign
    args.resume=True
    original=runner.source_pins()
    monkeypatch.setattr(runner,'source_pins',lambda:dict(original,changed_test_producer='different_sha256'))
    with pytest.raises(ValueError,match='Campaign identity changed: source_pins'):
        runner.run(args)


def test_completed_campaign_is_preserved_instead_of_overwritten(campaign):
    args,_=campaign
    before=(args.output/'summary.json').read_bytes()
    args.resume=True
    with pytest.raises(ValueError,match='Campaign completed; preserve it'):
        runner.run(args)
    assert (args.output/'summary.json').read_bytes()==before


def test_resume_rejects_changed_phonon_escape_time(campaign):
    args,_=campaign
    args.resume=True;args.tau_escape_ps=30.
    with pytest.raises(ValueError,match='Campaign identity changed: tau_escape_ps'):
        runner.run(args)


def test_infinite_escape_control_preserves_json_and_saved_phonons(campaign):
    args,_=campaign
    args.output=args.output.parent/'closed_control';args.tau_escape_ps=np.inf
    runner.run(args)
    identity=json.loads((args.output/'identity.json').read_text())
    summary=json.loads((args.output/'summary.json').read_text())
    assert identity['tau_escape_ps'] is None
    assert summary['stage51_complete'] is False
    for index,step in enumerate(summary['steps']):
        assert step['tau_escape_s'] is None
        assert step['phonon_escape_power_W']==0.
        assert 'infinity' in step['phonon_escape_scope']
        with np.load(args.output/('step_%02d'%index)/'candidate_fields.npz',allow_pickle=False) as fields:
            assert np.isinf(float(fields['tau_escape_s']))
            np.testing.assert_array_equal(fields['phonon_occupation_new'],fields['phonon_occupation_initial'])
            np.testing.assert_array_equal(fields['phonon_rate_per_s'],0.)
