"""Independent common-energy tests for the complete effective spatial system."""
import numpy as np
import pytest
from scipy.constants import Boltzmann as kB, hbar
from pysnspd.experimental.thermal_spatial_usadel import rectangular_graph
from pysnspd.experimental.effective_bcs_functional import EffectiveBCSFunctional,bcs_count_quadrature
from pysnspd.experimental.effective_bcs_kinetics import EffectiveBCSKinetics
from pysnspd.experimental.effective_bcs_ports import EffectiveNeutralPorts
from pysnspd.experimental.electrical_ports import ThesisCircuitParameters
from pysnspd.experimental.kinetic_events import PhononGrid
from pysnspd.experimental.effective_coupled import EffectiveState,EffectiveCoupledSystem


def system(Tb=4.,coupling=.0005):
    Eu=kB*8.65;D=5e-5;ell=np.sqrt(hbar*D/(2*Eu));thickness=7e-9
    graph=rectangular_graph(5,3,20e-9/ell,10e-9/ell)
    xi,w=bcs_count_quadrature(4.,order=3)
    fun=EffectiveBCSFunctional(graph,energy_scale_J=Eu,N0_per_J_m3=9.153297821093998e46,
        ell0_m=ell,thickness_m=thickness,count_nodes=xi,count_weights=w,
        gradient_law='bcs_radial',bath_temperature_K=Tb,diffusion_m2_s=D)
    bulk=fun.bulk(temperature_K=Tb,current_A=1e-6,width_m=10e-9,phase_step_m=5e-9)
    x=graph.coordinates_bar[:,0]
    d=bulk.amplitude_bar*np.exp(1j*bulk.q_bar*(x-x[-1]))
    ph=PhononGrid(np.r_[1e-8,np.linspace(.02,10.,48)],np.full(49,.08),
        lambda q:coupling*q*q,(0.,10.),'synthetic coupled identity test')
    kin=EffectiveBCSKinetics.from_graph(graph,ell0_m=ell,thickness_m=thickness,
        xi_bar=xi,weights_bar=w,phonons=ph,Eu_J=Eu,N0_per_J_m3=fun.N0_per_J_m3,
        D_m2_s=D,bath_temperature_K=Tb,tau_escape_s=15e-12,transport_order=12)
    left=np.flatnonzero(x==min(x));right=np.flatnonzero(x==max(x))
    ports=EffectiveNeutralPorts(graph,608.,left,right)
    circuit=ThesisCircuitParameters(Lk_ext_H=100e-9,V_bias_V=.01)
    model=EffectiveCoupledSystem(fun,kin,ports,circuit,current_scale_A=1e-6)
    state=EffectiveState(d,kin.fermi_dirac(d,Tb),np.tile(kin.bath_phonons,(len(d),1)),
        np.array([1e-6,1e-6,0.]))
    return model,state


def test_full_dc_has_no_artificial_voltage_heating_or_gap_rate():
    m,s=system()
    r=m.evaluate(s)
    assert abs(r.electrical.Vdev_V)<1e-14
    assert np.max(abs(r.gap_coordinate_per_ps))<1e-11
    assert np.max(r.kwt_heat_W)<1e-30
    assert abs(r.instant_power_residual_W)<1e-22


def test_common_energy_rhs_and_independent_finite_direction():
    m,s=system()
    s.p[m.free[2]]=m.kinetics.fermi_dirac(s.delta_bar,8.)[m.free[2]]
    s.n[m.free[2]]=m.kinetics.bose(8.)
    r=m.evaluate(s)
    assert np.min(r.kwt_heat_W)>=0
    assert abs(r.film_power_residual_W)<1e-22
    assert abs(r.instant_power_residual_W)<1e-22
    # Central finite difference is independent of the analytic chain observer.
    h=1e-6
    plus=m.shifted(s,r.gap_coordinate_per_ps*h,r.dp_per_ps*h,r.dn_per_ps*h,
        r.circuit_scaled_per_ps*h)
    minus=m.shifted(s,-r.gap_coordinate_per_ps*h,-r.dp_per_ps*h,-r.dn_per_ps*h,
        -r.circuit_scaled_per_ps*h)
    ep=m.evaluate(plus).energy_J;em=m.evaluate(minus).energy_J
    fd=(ep-em)/(2*h*1e-12)
    assert fd==pytest.approx(r.external_power_W,rel=2e-6,abs=1e-19)


def test_rigid_superconducting_cuts_use_phase_coordinates():
    m,s=system()
    u=np.zeros(m.ngap);u[-2:]=[2.,-3.]
    shifted=m.shifted(s,u,np.zeros_like(s.p),np.zeros_like(s.n),np.zeros(3))
    for nodes in m.cuts:
        np.testing.assert_allclose(abs(shifted.delta_bar[nodes]),abs(s.delta_bar[nodes]),rtol=3e-16)
    T=m.coordinate_tangent(s.delta_bar)
    hd=1e-7
    numerical=(m.shifted(s,hd*u,np.zeros_like(s.p),np.zeros_like(s.n),np.zeros(3)).delta_bar-
        m.shifted(s,-hd*u,np.zeros_like(s.p),np.zeros_like(s.n),np.zeros(3)).delta_bar)/(2*hd)
    expected=(T@u).reshape(-1,2)
    np.testing.assert_allclose(numerical,expected[:,0]+1j*expected[:,1],atol=2e-9)


def test_midpoint_dc_and_populations_are_not_mutated():
    m,s=system()
    before=s.copy()
    out,receipt=m.midpoint_step(s,.01)
    np.testing.assert_allclose(out.delta_bar,s.delta_bar,rtol=1e-12,atol=1e-13)
    np.testing.assert_array_equal(s.p,before.p)
    np.testing.assert_array_equal(s.n,before.n)
    assert receipt['defect']<2e-9
