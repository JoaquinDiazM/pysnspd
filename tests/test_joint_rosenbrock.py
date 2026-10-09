"""Physical sparse-operator tests for one coupled Rosenbrock-Euler candidate.

The small spatial fixture exercises actual Usadel/Matsubara/Nambu/CM
operators. It does not validate a photon trajectory or general material RHS.
"""
import numpy as np
import pytest
from scipy.constants import Boltzmann as kB, elementary_charge as e, hbar
from scipy.sparse import coo_matrix,diags
from scipy.sparse.linalg import spsolve
from pysnspd.experimental import adiabatic_matrix_usadel as ad
from pysnspd.experimental import causal_kinetic_regulator as reg
from pysnspd.experimental import retarded_spatial_usadel as ret
from pysnspd.experimental.bulk_current_reference import solve_bulk_reference,current_density_derivative
from pysnspd.experimental.thermal_spatial_usadel import rectangular_graph
from pysnspd.experimental.joint_evolution import BulkCutLaw,JointSpectralSlice,JointSnapshotSystem
from pysnspd.experimental.joint_rosenbrock import assemble_zeroth_kinetic_operator,JointRosenbrockSystem
from pysnspd.experimental.joint_state_tangent import ThermalBulkTangentBundle
from pysnspd.experimental.joint_nambu_rows import JointNambuRows
from pysnspd.experimental.nambu_moments import NambuMoments
from pysnspd.experimental.electrical_ports import ThesisCircuitParameters
from pysnspd.experimental.effective_material import EffectivePhononMaterial


def graph_laplacian(graph):
    i,j=graph.edges.T;c=graph.conductance
    return coo_matrix((np.r_[c,c,-c,-c],(np.r_[i,j,i,j],np.r_[i,j,j,i])),
                     shape=(graph.n_nodes,graph.n_nodes)).tocsr()


def test_zeroth_operator_matches_projected_full_nambu_residual():
    graph=rectangular_graph(4,3,2.4,1.3);n=graph.n_nodes
    rng=np.random.default_rng(251);gap=(1.2+.05*rng.normal(size=n))*np.exp(.1j*rng.normal(size=n))
    alpha=rng.normal(size=len(graph.edges))*.2;eta=.03;energy=1.7;kappa=3e-13
    a,b=ret.uniform_fields(gap,eta-1j*energy)
    f,ft,g=ret.fields(a,b);R=ad.Jet.of(ret.nambu_matrix(f,ft,g));B=ad.generator(gap,energy)
    operator=assemble_zeroth_kinetic_operator(graph,R,B,numerical_eta=eta,kappa_s=kappa,alpha=alpha)
    h=rng.uniform(-.3,.3,(n,2))
    result,_,_,_=reg.kinetic_residual(graph,R,ad.distribution(h[:,0],h[:,1]),B,
        numerical_eta=eta,kappa_s=kappa,alpha=alpha)
    np.testing.assert_allclose((operator@h.ravel()).reshape(n,2),ad.projections(result.value),
        rtol=2e-13,atol=2e-15)


def test_normal_limit_has_negative_laplacian_and_correct_temporal_mass():
    graph=rectangular_graph(5,3,3.,1.);n=graph.n_nodes;kappa=4e-13
    R=ad.Jet.of(np.broadcast_to(ad.TAU3,(n,2,2)));B=ad.generator(np.zeros(n),1.3)
    L=assemble_zeroth_kinetic_operator(graph,R,B,numerical_eta=.02,kappa_s=kappa)
    expected=-graph_laplacian(graph)
    np.testing.assert_allclose(L[0::2,0::2].toarray(),expected.toarray(),atol=1e-15)
    np.testing.assert_allclose(L[1::2,1::2].toarray(),expected.toarray(),atol=1e-15)
    np.testing.assert_array_equal(L[0::2,1::2].toarray(),0.)
    np.testing.assert_array_equal(L[1::2,0::2].toarray(),0.)
    M=reg.assemble_storage(graph,R,B,numerical_eta=.02,kappa_s=kappa)
    np.testing.assert_allclose(M.toarray(),-kappa*np.diag(np.repeat(graph.area_weights,2)),
        rtol=2e-15,atol=1e-28)


def test_m_plus_dt_l_solves_actual_normal_diffusion_with_no_clipping():
    graph=rectangular_graph(7,3,3.,1.);n=graph.n_nodes;kappa=4e-13;dt=5e-12
    R=ad.Jet.of(np.broadcast_to(ad.TAU3,(n,2,2)));B=ad.generator(np.zeros(n),1.3)
    L=assemble_zeroth_kinetic_operator(graph,R,B,numerical_eta=0.,kappa_s=kappa)
    M=reg.assemble_storage(graph,R,B,numerical_eta=0.,kappa_s=kappa)
    h=np.column_stack((.2+.35*np.cos(np.pi*graph.coordinates_bar[:,0]/3),np.zeros(n)))
    rate=spsolve(M+dt*L,-L@h.ravel()).reshape(n,2)
    endpoint=h+dt*rate
    expected=spsolve(kappa*diags(graph.area_weights)+dt*graph_laplacian(graph),
        kappa*graph.area_weights*h[:,0])
    np.testing.assert_allclose(endpoint[:,0],expected,rtol=2e-13,atol=2e-15)
    assert np.min(endpoint[:,0])>=np.min(h[:,0])-1e-14
    assert np.max(endpoint[:,0])<=np.max(h[:,0])+1e-14
    assert np.dot(graph.area_weights,endpoint[:,0])==pytest.approx(
        np.dot(graph.area_weights,h[:,0]),rel=2e-13)
    np.testing.assert_allclose(M@rate.ravel()+L@h.ravel()+dt*L@rate.ravel(),0.,atol=2e-14)


def physical_fixture(common_temperature_K=None):
    # Small two-cut superconducting strip with a real finite Matsubara
    # thermal Hessian and four real-energy kinetic modes, in SI scales.
    graph=rectangular_graph(3,2,2.,1.);n=graph.n_nodes
    fixed=np.array([0,1,4,5]);left=fixed[:2];right=fixed[2:]
    Tc=8.65;t=.25;Eu=kB*Tc;N0=9e46;v0=3e-25;I0=Eu/(2*e*608.)
    temperature=np.full(n,t*Tc) if common_temperature_K is None else np.asarray(common_temperature_K,float)
    if temperature.shape!=(n,) or np.any(temperature<=0):raise ValueError('One positive common FD/BE temperature per fixture node required')
    bulk=solve_bulk_reference(t,16,0.);gap=np.full(n,bulk.gap_bar,dtype=complex)
    law=BulkCutLaw(bulk.gap_bar,0.,0.,1/(I0*current_density_derivative(bulk)),
        'Small q=0 homogeneous superconducting branch; synthetic component fixture')
    energy=Eu*np.array([.3,.8,1.6,2.4]);weights=Eu*np.array([.4,.4,.4,.4])
    slices=[];Kth=[];Fth=[];alpha=np.zeros(len(graph.edges))
    for E in energy:
        a,b=ret.uniform_fields(gap,.03-1j*E/Eu)
        f,ft,g=ret.fields(a,b);R0=ret.nambu_matrix(f,ft,g)
        hl=np.tanh(E/(2*kB*temperature));h=np.column_stack((hl,np.zeros(n)))
        hE=np.column_stack(((1-hl**2)*Eu/(2*kB*temperature),np.zeros(n)))
        sl=JointSpectralSlice(graph,gap,alpha,R0,E,.03*Eu,Eu,h,hE,fixed,law,
            projected_source=np.zeros((n,2)),source_contract='Thermal FD/BE component fixture only')
        slices.append(sl)
        hthermal=np.full(n,np.tanh(E/(2*t*Eu)))
        _,flux,K=ad.kinetic_residual(graph,ad.Jet.of(R0),ad.distribution(hthermal,np.zeros(n)),
            ad.generator(gap,E/Eu),alpha=alpha,kappa=hbar/(2*Eu))
        Kth.append(K.value);Fth.append(flux.value)
    moments=NambuMoments(energy,weights,graph.area_weights*v0,N0,Eu,I0)
    rows=JointNambuRows(moments,graph.edges,left_nodes=left,right_nodes=right)
    system=JointSnapshotSystem(graph,gap,alpha,slices,rows,K_thermal_zero=Kth,
        thermal_gradient_J_per_gapbar=np.zeros(n,complex),thermal_edge_current_A=np.zeros(len(graph.edges)),
        thermal_edge_flux_zero=Fth,circuit=ThesisCircuitParameters(Lk_ext_H=1.3e-7,V_bias_V=0.),
        circuit_state=[0.,0.,0.],bulk_law=law,
        reference_contract='Synthetic finite-M16 equilibrium component fixture, same gap and links')
    material=EffectivePhononMaterial(Eu*np.array([0.,.3,1.]),np.array([0.,.001,.01]),
        np.array([0.,1e65,1e66]),N0,{'material_id':'component_fixture_only'})
    system.initial_data=dict(receipt={'source_must_not_be_applied_again':True},
        inputs={'summary':{'plan':{'T_K':t*Tc}}},effective_material=material,
        common_temperature_K=temperature)
    thermal=ThermalBulkTangentBundle(graph,gap,t,bulk.epsilon_bar,
        np.broadcast_to(bulk.u[:,None],(16,n)),alpha,fixed_nodes=fixed,bulk_q_bar=0.,
        energy_scale_J=N0*Eu**2*v0,current_unit_A=I0)
    return system,thermal


def test_zero_linearized_step_keeps_the_snapshot_operator_exactly():
    system,_=physical_fixture();sl=system.slices[1]
    dd=np.zeros(system.graph.n_nodes,complex);dd[system.free]=np.array([.1+.2j,-.04+.03j])/sl.kappa
    phi=np.linspace(-1e-6,2e-6,system.graph.n_nodes)
    before=sl.evaluate(dd,phi,0.)
    sl.configure_linearized_step(0.)
    after=sl.evaluate(dd,phi,0.)
    for j in range(4):
        if j<3:
            for field in ('value','energy','time','first'):
                np.testing.assert_array_equal(getattr(before[j],field),getattr(after[j],field))
        else:np.testing.assert_array_equal(before[j],after[j])


def test_real_spatial_rosenbrock_fixture_is_stationary_and_affine():
    system,thermal=physical_fixture();step=JointRosenbrockSystem(system,2e-15,thermal)
    zero=step.evaluate(np.zeros(step.size))
    np.testing.assert_allclose(zero.packed_rows,0.,atol=2e-12)
    assert zero.kinetic_residual_max<1e-12
    assert zero.full_nambu_residual_max_per_mass<1e-10
    rng=np.random.default_rng(351);a=rng.normal(size=step.size)*.0001;b=rng.normal(size=step.size)*.0001
    va=step.evaluate(a);vb=step.evaluate(b);both=step.evaluate(a+b)
    np.testing.assert_allclose(both.packed_rows,va.packed_rows+vb.packed_rows-zero.packed_rows,
        rtol=3e-10,atol=3e-12)
    assert np.max(abs(va.distribution_rate_per_s))>1e5
    assert np.max(abs(va.circuit_rate_SI_per_s))>0
    assert zero.scope.startswith('Single linearly implicit candidate')
    np.testing.assert_array_equal(zero.phonon_rate_per_s,0.)


def test_rosenbrock_steps_do_not_change_snapshot_or_each_others_dt():
    system,thermal=physical_fixture();a=JointRosenbrockSystem(system,1e-15,thermal)
    x=np.random.default_rng(514).normal(size=a.size)*.0001
    before=a.evaluate(x).packed_rows.copy()
    b=JointRosenbrockSystem(system,4e-15,thermal)
    assert not hasattr(system.slices[0],'linearized_step_s')
    assert a.slices[0].linearized_step_s==1e-15
    assert b.slices[0].linearized_step_s==4e-15
    np.testing.assert_array_equal(a.evaluate(x).packed_rows,before)


def test_phonon_escape_is_one_positive_backward_euler_row_with_energy_balance():
    bath=.25*8.65;temperatures=np.full(6,bath);temperatures[2]=4.;temperatures[3]=.8*bath
    system,thermal=physical_fixture(common_temperature_K=temperatures)
    dt=2e-12;tau=15e-12
    step=JointRosenbrockSystem(system,dt,thermal,tau_escape_s=tau)
    old=step.phonon_occupation_initial;new=step.phonon_occupation;nb=step.phonon_bath
    expected=nb+(old-nb)/(1+dt/tau)
    np.testing.assert_array_equal(new,expected)
    np.testing.assert_allclose((new-old)/dt,step.phonon_rate,rtol=3e-13,atol=2e-5)
    np.testing.assert_allclose(step.phonon_rate,-(new-nb)/tau,rtol=1e-15,atol=0.)
    assert np.all(new>=0)
    assert np.all(new>=np.minimum(old,nb))
    assert np.all(new<=np.maximum(old,nb))
    spectral_measure=(step.phonon_energy_J*step.phonon_weights_J*step.phonon_DOS_per_J_m3)[:,None]
    measure=spectral_measure*system.rows.moments.volume[None,:]
    energy_difference=float(np.sum(measure*(new-old)))
    assert energy_difference==pytest.approx(-dt*step.phonon_escape_power_W,rel=3e-13,abs=1e-36)
    # Transfer to the bath is signed: a node colder than the bath gains energy.
    assert np.all(step.phonon_rate[:,2]<0)
    assert np.all(step.phonon_rate[:,3]>0)
    result=step.evaluate(np.zeros(step.size))
    np.testing.assert_array_equal(result.phonon_rate_per_s,step.phonon_rate)
    assert 'linear phonon escape implicit' in result.material_rhs_contract
    assert 'not converged backward Euler' in result.scope


def test_phonons_at_bath_are_stationary_with_finite_escape_time():
    system,thermal=physical_fixture()
    step=JointRosenbrockSystem(system,1e-12,thermal,tau_escape_s=15e-12)
    np.testing.assert_array_equal(step.phonon_occupation,step.phonon_occupation_initial)
    np.testing.assert_array_equal(step.phonon_rate,0.)
    assert step.phonon_escape_power_W==0.


def test_infinite_phonon_escape_time_is_an_exact_closed_control():
    temperature=np.full(6,.25*8.65);temperature[2]=4.
    system,thermal=physical_fixture(common_temperature_K=temperature)
    step=JointRosenbrockSystem(system,2e-12,thermal,tau_escape_s=np.inf)
    np.testing.assert_array_equal(step.phonon_occupation,step.phonon_occupation_initial)
    np.testing.assert_array_equal(step.phonon_rate,0.)
    assert step.phonon_escape_power_W==0.
    assert np.isinf(step.tau_escape_s)


@pytest.mark.parametrize('tau',[0.,-1e-12,np.nan])
def test_invalid_phonon_escape_time_is_not_repaired(tau):
    system,thermal=physical_fixture()
    with pytest.raises(ValueError,match='Positive escape time'):
        JointRosenbrockSystem(system,2e-15,thermal,tau_escape_s=tau)
