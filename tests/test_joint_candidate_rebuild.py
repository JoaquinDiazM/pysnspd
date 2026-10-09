"""Small finite endpoint roots and independent energy/support diagnostics."""
import json
from pathlib import Path
from tempfile import TemporaryDirectory
import numpy as np
import pytest
from scipy.constants import Boltzmann as kB
from test_joint_rosenbrock import physical_fixture
from pysnspd.experimental.bulk_current_reference import homogeneous_spectrum
from pysnspd.experimental.thermal_internal_energy import evaluate_internal_energy
from pysnspd.experimental import thermal_spatial_usadel as thermal
from sandbox.stage5_photon import joint_candidate_rebuild as rebuild


def test_fixed_energy_derivative_preserves_two_modes_without_temperature_fit():
    energy=np.array([.1,.31,.8,1.4,2.1])
    h=np.empty((len(energy),3,2))
    h[:,:,0]=energy[:,None]**2*np.array([.1,.2,.3])
    h[:,:,1]=energy[:,None]*np.array([.01,-.02,.03])
    result=rebuild.fixed_energy_distribution_derivative(energy*2e-22,2e-22,h)
    np.testing.assert_allclose(result[:,:,0],2*energy[:,None]*np.array([.1,.2,.3]),
        rtol=4e-14,atol=1e-15)
    np.testing.assert_allclose(result[:,:,1],np.broadcast_to([.01,-.02,.03],(len(energy),3)),
        rtol=4e-14,atol=1e-15)


def test_current_carrying_cut_temperature_tangent_matches_finite_difference():
    epsilon=np.array([.4,1.2,2.]);t=.3;q=.45
    gap=np.array([1.2*np.exp(.3j),1.25*np.exp(-.1j)])
    u,du=rebuild.thermal_cut_fields(epsilon,gap,q,t)
    eps=2e-6
    plus,_=rebuild.thermal_cut_fields(epsilon*(t+eps)/t,gap,q,t+eps)
    minus,_=rebuild.thermal_cut_fields(epsilon*(t-eps)/t,gap,q,t-eps)
    np.testing.assert_allclose((plus-minus)/(2*eps),du,rtol=3e-10,atol=2e-10)
    # q!=0 is the intended superconducting continuation, not the BCS law.
    assert np.max(abs(du+u/t))>.05


def test_support_rejection_never_runs_a_spectral_query(monkeypatch):
    snapshot,_=physical_fixture()
    candidate=dict(gap_bar_new=snapshot.gap,alpha_new=snapshot.alpha,
        h_new=np.zeros((len(snapshot.slices),snapshot.graph.n_nodes,2)),
        circuit_state_new_SI=np.array([1e-6,1e-6,0.]))
    candidate["h_new"][0,0,0]=1.01
    parent=Path('/home/jdiaz/scratch/tmp/pysnspd')
    with TemporaryDirectory(dir=parent,prefix="test_rebuild_support_") as directory:
        path=Path(directory)
        np.savez(path/'candidate.npz',**candidate)
        monkeypatch.setattr(rebuild,'solve_bulk_reference',lambda *a,**k:
            pytest.fail("No roots may run after a support rejection"))
        result=rebuild.rebuild_candidate(snapshot,path/'candidate.npz',path/'output')
        assert result['status']=='FINITE_STATE_REJECTED_BEFORE_SPECTRAL_REBUILD'
        assert result['screen']['reason']=='Pauli_support'
        assert result['finite_state_rebuilt'] is False
        assert result['energy_admission'] is False
        assert result['time_admission'] is False
        assert json.loads((path/'output/summary.json').read_text())==result
        assert candidate['h_new'][0,0,0]==1.01


def finite_fixture(directory):
    snapshot,bundle=physical_fixture()
    path=Path(directory);observer=path/'observer';(observer/'spectral_receipts').mkdir(parents=True)
    plan=dict(T_K=bundle.t*8.65,Tc_K=8.65,matsubara_count=16,width_m=80e-9,
        diffusion_m2_s=5e-5,sheet_resistance_ohm=608.)
    modes=[]
    for i,mode in enumerate(bundle.modes):
        np.savez(observer/'spectral_receipts'/('mode_%05d.npz'%i),
            u=mode.u,epsilon_bar=mode.epsilon)
        modes.append(thermal.solve_frequency(snapshot.graph,snapshot.gap,mode.epsilon,
            snapshot.alpha,fixed_nodes=snapshot.fixed,fixed_u=mode.u[snapshot.fixed],
            initial_u=mode.u,tol=1e-9))
    eps=np.array([mode.epsilon for mode in bundle.modes])
    _,cutT=rebuild.thermal_cut_fields(eps,snapshot.gap[snapshot.fixed],0.,bundle.t)
    internal=evaluate_internal_energy(snapshot.graph,snapshot.gap,bundle.t,modes,snapshot.alpha,
        energy_scale_J=bundle.scale,boundary_contract='fixture same q0 finite cut law',
        boundary_temperature_derivatives=cutT)
    snapshot.initial_data['inputs'].update(summary={'plan':plan},observed=dict(
        energy_scale_J=bundle.scale,internal_energy_J=internal.internal_energy_J))
    snapshot.initial_data['identity']={'plan':{'thermal_observer':str(observer)}}
    phononE=snapshot.Eu*np.array([.2,.6])
    bath=np.broadcast_to(1/np.expm1(phononE[:,None]/(kB*plan['T_K'])),
        (2,snapshot.graph.n_nodes)).copy()
    candidate=dict(gap_bar_new=snapshot.gap.copy(),alpha_new=snapshot.alpha.copy(),
        h_new=np.asarray([sl.h for sl in snapshot.slices]),
        # A positive small current perturbation selects an actual q>0 branch.
        # This component candidate is not claimed to be a solved detector step.
        circuit_state_new_SI=np.array([1e-7,1e-7,0.]),dt_s=np.asarray(1e-15),
        phonon_energy_J=phononE,phonon_weights_J=np.full(2,.2*snapshot.Eu),
        phonon_DOS_per_J_m3=np.array([1e65,2e65]),
        phonon_occupation_initial=bath*1.2,phonon_occupation_new=bath*1.1,
        phonon_bath=bath,tau_escape_s=np.asarray(15e-12))
    return snapshot,candidate


def test_small_real_roots_finite_Matsubara_U_and_neutral_energy_are_rebuilt():
    parent=Path('/home/jdiaz/scratch/tmp/pysnspd')
    with TemporaryDirectory(dir=parent,prefix="test_rebuild_actual_") as directory:
        path=Path(directory);snapshot,candidate=finite_fixture(path)
        before=candidate['gap_bar_new'].copy()
        np.savez(path/'candidate.npz',**candidate)
        result=rebuild.rebuild_candidate(snapshot,path/'candidate.npz',path/'output',workers=1)
        assert result['status']=='FINITE_STATE_ZERO_ORDER_REBUILT_DIAGNOSTIC_ONLY',result
        assert result['finite_state_rebuilt'] is True
        assert result['thermal_modes']==16 and result['real_energy_modes']==4
        assert result['bulk_continuation']['q_new_bar']>0
        assert result['bulk_continuation']['maximum_actual_cut_amplitude_mismatch_bar']>0
        assert result['charge_gauge_energy_completion_missing'] is True
        assert result['general_collision_provider_available'] is False
        assert result['new_first_gradient_fields_missing'] is True
        assert result['backwards_Euler_converged'] is False
        assert result['energy_admission'] is False and result['time_admission'] is False
        assert result['energy']['phonon_escape_out_new_W']>0
        assert result['energy']['dt_times_phonon_escape_out_new_J']==pytest.approx(
            1e-15*result['energy']['phonon_escape_out_new_W'])
        with np.load(result['fields_path']) as saved:
            np.testing.assert_array_equal(saved['gap_bar'],before)
            np.testing.assert_array_equal(saved['h'],candidate['h_new'])
            assert saved['R0'].shape==(4,snapshot.graph.n_nodes,2,2)
            np.testing.assert_allclose(saved['h_energy'],rebuild.fixed_energy_distribution_derivative(
                snapshot.rows.moments.energy,snapshot.Eu,candidate['h_new']),atol=0)
        assert rebuild.rebuild_candidate(snapshot,path/'candidate.npz',path/'output')==result


def test_failed_finite_root_is_a_returned_diagnostic(monkeypatch):
    parent=Path('/home/jdiaz/scratch/tmp/pysnspd')
    with TemporaryDirectory(dir=parent,prefix="test_rebuild_failure_") as directory:
        path=Path(directory);snapshot,candidate=finite_fixture(path)
        np.savez(path/'candidate.npz',**candidate)
        def failure(*args,**kwargs):
            raise RuntimeError('fixture stationary branch failure')
        monkeypatch.setattr(rebuild.ret,'solve',failure)
        result=rebuild.rebuild_candidate(snapshot,path/'candidate.npz',path/'output',workers=1)
        assert result['status']=='FINITE_STATE_SPECTRAL_REBUILD_REJECTED'
        assert result['failed_query']['kind']=='retarded'
        assert result['failed_query']['index']==0
        assert 'fixture stationary branch failure' in result['failed_query']['error']
        assert result['finite_state_rebuilt'] is False
        assert result['energy_admission'] is False and result['time_admission'] is False


def test_tiny_Pauli_excess_is_only_an_observer_not_an_admitted_state():
    snapshot,_=physical_fixture()
    h=np.zeros((len(snapshot.slices),snapshot.graph.n_nodes,2))
    h[0,0,0]=1.+1e-12
    candidate=dict(gap_bar_new=snapshot.gap,alpha_new=snapshot.alpha,h_new=h,
        circuit_state_new_SI=np.array([1e-6,1e-6,0.]))
    result=rebuild.candidate_state_screen(snapshot,candidate)
    assert result["ok_for_diagnostic_rebuild"] is True
    assert result["support_ok"] is False
    assert result["accepted_physical_time_steps"]==0
    assert result["no_clipping"] is True
    assert h[0,0,0]>1.
    rho=np.full(h.shape[:2],.4)
    observed=rebuild.support_excess_energy_J(snapshot.rows.moments.energy,
        snapshot.rows.moments.weights,snapshot.rows.moments.volume,
        snapshot.rows.moments.N0,rho,h)
    expected=2*snapshot.rows.moments.N0*snapshot.rows.moments.energy[0]*(
        snapshot.rows.moments.weights[0]*snapshot.rows.moments.volume[0]*.4*(h[0,0,0]-1.))
    assert observed==pytest.approx(expected,rel=1e-13)


def test_actual_two_worker_endpoint_pool_matches_serial_physical_fields():
    parent=Path('/home/jdiaz/scratch/tmp/pysnspd')
    with TemporaryDirectory(dir=parent,prefix="test_rebuild_parallel_") as directory:
        path=Path(directory);snapshot,candidate=finite_fixture(path)
        np.savez(path/'candidate.npz',**candidate)
        serial=rebuild.rebuild_candidate(snapshot,path/'candidate.npz',path/'serial',workers=1)
        parallel=rebuild.rebuild_candidate(snapshot,path/'candidate.npz',path/'parallel',workers=2)
        assert parallel["status"]=="FINITE_STATE_ZERO_ORDER_REBUILT_DIAGNOSTIC_ONLY",parallel
        assert parallel["resource_budget"]["workers"]==2
        assert parallel["resource_budget"]["threads_per_process"]==1
        assert parallel["resource_budget"]["cpu_fraction"]<=.9
        assert parallel["resource_budget"]["memory_fraction"]<=.9
        for key in ("thermal_internal_new_J","electronic_neutral_adiabatic_new_J",
                "total_neutral_adiabatic_new_J","support_excess_energy_new_DOS_J"):
            assert parallel["energy"][key]==pytest.approx(serial["energy"][key],rel=1e-13,abs=1e-30)
        with np.load(serial["fields_path"]) as first,np.load(parallel["fields_path"]) as second:
            for key in ("R0","K0","h_energy","zeroth_scaled_rows"):
                np.testing.assert_allclose(second[key],first[key],rtol=1e-13,atol=1e-13)
