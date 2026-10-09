"""DAE-discretization tests; these are not spatial detector admission tests."""
import numpy as np
import pytest
from scipy.constants import Boltzmann,hbar
from scipy.sparse import csr_matrix
from pysnspd.experimental.bulk_current_reference import (
    homogeneous_spectrum, solve_bulk_reference)
from pysnspd.experimental.joint_implicit_step import (
    ImplicitStepRejected, solve_backward_euler)


def options(n,*,mask=None,lower=None,upper=None):
    return dict(differential_mask=np.ones(n,bool) if mask is None else mask,
        variable_scales=np.ones(n),row_scales=np.ones(n),
        lower_bounds=np.full(n,-np.inf) if lower is None else lower,
        upper_bounds=np.full(n,np.inf) if upper is None else upper)


def test_cold_bulk_algebraic_gap_has_stiffness_without_radial_mobility():
    # The homogeneous finite-F256 dirty-Usadel self-consistency equation.
    # Its zero time-derivative coefficient does not remove state stiffness.
    # Deliberate off-manifold start is a nonlinear-solver contrast only;
    # the recorded projection is never admitted as photon time evolution.
    t=.9/8.65;q=.30;reference=solve_bulk_reference(t,256,q)
    ep=reference.epsilon_bar;coefficient=np.log(t)+2*np.pi*t*np.sum(1/ep)
    def residual(y,ydot):
        _,f,_=homogeneous_spectrum(y[0],ep,q*q)
        return np.array([coefficient*y[0]-2*np.pi*t*np.sum(f)])
    def jacobian(y,ydot):
        _,_,g=homogeneous_spectrum(y[0],ep,q*q)
        stiffness=coefficient-2*np.pi*t*np.sum(g**3/(ep+q*q*g**3))
        return np.array([[stiffness]]),np.zeros((1,1))
    kwargs=options(1,mask=np.zeros(1,bool),lower=[.5],upper=[2.])
    previous=np.array([.98*reference.gap_bar])
    result=solve_backward_euler(previous,1e-15,residual=residual,jacobian=jacobian,**kwargs)
    assert result.state[0]==pytest.approx(reference.gap_bar,rel=1e-11)
    assert result.initial_residual_max>1e-3
    assert result.maximum_increment_scaled>1e-2
    assert np.array_equal(result.derivative_per_s,[0.])
    assert not result.temporal_admission
    assert not result.energy_admission
    assert result.no_clipping_or_regularization
    np.testing.assert_array_equal(previous,[.98*reference.gap_bar])


def paired_normal_emission():
    # A physical normal-state golden-rule event: high <-> low + phonon.
    # u=1-2f enforces the same Pauli support as hL +/- hT.
    # Equal electronic cell counts and an explicit capacity ratio preserve
    # f_low+f_high and E_low*f_low+E_high*f_high+Omega*n/ratio.
    Eu=Boltzmann*8.65;weight=Eu*.13;alpha=.004
    gamma=2*np.pi*weight*alpha/hbar;ratio=.7
    vector=gamma*np.array([-2.,2.,ratio])
    def event(y):
        fl,fh=(1-y[:2])/2;n=y[2]
        flow=(n+1)*fh*(1-fl)-n*fl*(1-fh)
        gradient=np.array([.5*(n+fh),-.5*(n+1-fl),fh-fl])
        return flow,gradient
    def residual(y,ydot):return ydot-vector*event(y)[0]
    def jacobian(y,ydot):return -np.outer(vector,event(y)[1]),np.eye(3)
    return gamma,ratio,Eu,residual,jacobian


def test_joint_electron_phonon_backward_euler_stays_physical_and_conserves():
    gamma,ratio,Eu,residual,jacobian=paired_normal_emission()
    previous=np.array([.8,.2,0.]);dt=100/gamma
    # The same large explicit step would leave physical support.
    explicit=previous-dt*residual(previous,np.zeros(3))
    assert np.max(abs(explicit[:2]))>1
    kwargs=options(3,lower=[-1.,-1.,0.],upper=[1.,1.,np.inf])
    kwargs['row_scales']=np.full(3,gamma)
    result=solve_backward_euler(previous,dt,residual=residual,jacobian=jacobian,**kwargs)
    assert result.accepted_residual_max<1e-9
    assert np.all(abs(result.state[:2])<=1)
    assert result.state[2]>=0
    f0=(1-previous[:2])/2;f1=(1-result.state[:2])/2
    assert np.sum(f1)==pytest.approx(np.sum(f0),abs=2e-14)
    energies=Eu*np.array([.31,.73]);omega=Eu*.42
    energy0=np.dot(energies,f0)+omega*previous[2]/ratio
    energy1=np.dot(energies,f1)+omega*result.state[2]/ratio
    assert energy1==pytest.approx(energy0,rel=2e-13,abs=1e-40)
    assert not result.temporal_admission  # Accuracy still needs dt refinement.


def test_differential_and_algebraic_unknowns_share_one_discrete_solve():
    # Index-one DAE ydot+z=0, z=y^2. Algebraic state has no fictitious rate.
    def residual(y,ydot):return np.array([ydot[0]+y[1],y[1]-y[0]**2])
    def jacobian(y,ydot):
        return csr_matrix([[0.,1.],[-2*y[0],1.]]),csr_matrix([[1.,0.],[0.,0.]])
    result=solve_backward_euler([1.,1.],.2,residual=residual,jacobian=jacobian,
        **options(2,mask=np.array([True,False]),lower=[0.,0.],upper=[np.inf,np.inf]))
    exact=(-1+np.sqrt(1.8))/.4
    assert result.state[0]==pytest.approx(exact,rel=1e-10)
    assert result.state[1]==pytest.approx(exact**2,rel=1e-10)
    assert result.derivative_per_s[1]==0.
    assert np.max(abs(residual(result.state,result.derivative_per_s)))<1e-9


def test_no_solution_inside_support_is_rejected_without_clipping():
    previous=np.array([1.])
    def residual(y,ydot):return ydot-np.ones(1)
    def jacobian(y,ydot):return np.zeros((1,1)),np.eye(1)
    with pytest.raises(ImplicitStepRejected,match='exhausted support'):
        solve_backward_euler(previous,1.,residual=residual,jacobian=jacobian,
            **options(1,lower=[-1.],upper=[1.]))
    np.testing.assert_array_equal(previous,[1.])


def test_unclosed_discrete_dae_is_rejected_instead_of_regularized():
    with pytest.raises(ImplicitStepRejected,match='singular'):
        solve_backward_euler([0.],1e-12,residual=lambda y,d:np.array([1.]),
            jacobian=lambda y,d:(np.zeros((1,1)),np.zeros((1,1))),**options(1))


def test_support_violation_in_initial_guess_is_not_repaired():
    with pytest.raises(ValueError,match='Initial guess violates support'):
        solve_backward_euler([1.],1e-12,residual=lambda y,d:y-1.,
            jacobian=lambda y,d:(np.eye(1),np.zeros((1,1))),initial_guess=[1.01],
            **options(1,lower=[-1.],upper=[1.]))


def test_fraction_to_boundary_only_evaluates_physical_states():
    visited=[]
    def residual(y,ydot):
        assert 0<=y[0]<=1
        visited.append(y[0]);return np.array([np.exp(10*y[0])-100.])
    def jacobian(y,ydot):return np.array([[10*np.exp(10*y[0])]]),np.zeros((1,1))
    result=solve_backward_euler([0.],1e-12,residual=residual,jacobian=jacobian,
        **options(1,mask=np.zeros(1,bool),lower=[0.],upper=[1.]))
    assert result.state[0]==pytest.approx(np.log(100)/10,rel=1e-11)
    assert result.minimum_line_search_fraction<1
    assert all(0<=y<=1 for y in visited)


def test_full_row_residual_cannot_be_replaced_by_step_smallness():
    with pytest.raises(ImplicitStepRejected):
        solve_backward_euler([0.],1.,residual=lambda y,d:np.array([np.exp(y[0])+1.]),
            jacobian=lambda y,d:(np.array([[np.exp(y[0])]]),np.zeros((1,1))),
            max_iterations=5,**options(1))
