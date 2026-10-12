"""Physical residual regression for the full midpoint; no raw-p tolerance claim."""
import numpy as np
import pytest
from scipy.constants import Boltzmann as kB,hbar,Planck,speed_of_light
from pysnspd.experimental.thermal_spatial_usadel import rectangular_graph
from pysnspd.experimental.effective_bcs_functional import EffectiveBCSFunctional,bcs_count_quadrature
from pysnspd.experimental.effective_bcs_kinetics import EffectiveBCSKinetics
from pysnspd.experimental.effective_bcs_ports import EffectiveNeutralPorts
from pysnspd.experimental.electrical_ports import ThesisCircuitParameters
from pysnspd.experimental.kinetic_events import PhononGrid
from pysnspd.experimental.effective_coupled import EffectiveState,EffectiveCoupledSystem


def physical_system(length_nm=20.,nx=5,retained_energy_J=None):
    """Same film/spacing/current, with independently variable strip length."""
    Tb=4.;Eu=kB*8.65;D=5e-5;ell=np.sqrt(hbar*D/(2*Eu));thickness=7e-9
    graph=rectangular_graph(nx,3,length_nm*1e-9/ell,10e-9/ell)
    xi,w=bcs_count_quadrature(4.,order=3)
    f=EffectiveBCSFunctional(graph,energy_scale_J=Eu,N0_per_J_m3=9.153297821093998e46,
        ell0_m=ell,thickness_m=thickness,count_nodes=xi,count_weights=w,
        gradient_law='bcs_radial',bath_temperature_K=Tb,diffusion_m2_s=D)
    bulk=f.bulk(temperature_K=Tb,current_A=1e-6,width_m=10e-9,phase_step_m=5e-9)
    x=graph.coordinates_bar[:,0]
    d=bulk.amplitude_bar*np.exp(1j*bulk.q_bar*(x-x[-1]))
    ph=PhononGrid(np.r_[1e-8,np.linspace(.02,10.,48)],np.full(49,.08),
        lambda q:.0005*q*q,(0.,10.),'synthetic coupled physical-residual regression')
    k=EffectiveBCSKinetics.from_graph(graph,ell0_m=ell,thickness_m=thickness,
        xi_bar=xi,weights_bar=w,phonons=ph,Eu_J=Eu,N0_per_J_m3=f.N0_per_J_m3,
        D_m2_s=D,bath_temperature_K=Tb,tau_escape_s=15e-12,transport_order=12)
    left=np.flatnonzero(x==min(x));right=np.flatnonzero(x==max(x))
    ports=EffectiveNeutralPorts(graph,608.,left,right)
    circuit=ThesisCircuitParameters(Lk_ext_H=100e-9,V_bias_V=.01)
    m=EffectiveCoupledSystem(f,k,ports,circuit,current_scale_A=1e-6,
        kinetic_residual_energy_scale_J=retained_energy_J)
    state=EffectiveState(d,k.fermi_dirac(d,Tb),np.tile(k.bath_phonons,(len(d),1)),
        np.array([1e-6,1e-6,0.]))
    return m,state


@pytest.mark.parametrize('wavelength_nm',[775.,1550.])
def test_supplied_retained_energy_is_invariant_under_domain_length(wavelength_nm):
    # The supplied photon energy may EXCEED the small strip's condensation scale.
    retained=Planck*speed_of_light/(wavelength_nm*1e-9)*.667
    short,_=physical_system(retained_energy_J=retained)
    long,_=physical_system(length_nm=40.,nx=9,retained_energy_J=retained)
    assert np.sum(short.functional.node_energy_scale_J)<retained
    assert long.kinetic_residual_energy_scale_J==short.kinetic_residual_energy_scale_J==retained
    dc,_=physical_system()
    assert dc.kinetic_residual_energy_scale_J==pytest.approx(
        np.sum(dc.functional.node_energy_scale_J),rel=1e-15)


def test_accepted_midpoint_reports_independent_unsigned_energy_and_gap_guard():
    m,state=physical_system()
    cold=m.evaluate(state)
    node=m.free[len(m.free)//2]
    state.p[node]=m.kinetics.fermi_dirac(state.delta_bar,8.)[node]
    state.n[node]=m.kinetics.bose(8.)
    initial=m.evaluate(state)
    retained=initial.domain_energy_J-cold.domain_energy_J
    assert retained>0
    # Use the actual physical excitation energy, not the total circuit energy.
    m=EffectiveCoupledSystem(m.functional,m.kinetics,m.ports,m.circuit,current_scale_A=1e-6,
        kinetic_residual_energy_scale_J=retained)
    h=.005;tol=2e-9
    end,receipt=m.midpoint_step(state,h,tolerance=tol)
    du=np.zeros(m.ngap)
    dd=end.delta_bar-state.delta_bar
    du[:2*len(m.free)]=np.column_stack((dd[m.free].real,dd[m.free].imag)).ravel()
    for j,nodes in enumerate(m.cuts):
        du[2*len(m.free)+j]=np.angle(end.delta_bar[nodes[0]]/state.delta_bar[nodes[0]])
    dp=end.p-state.p;dn=end.n-state.n
    dc=(end.circuit_SI-state.circuit_SI)/m.circuit_scales
    mid=m.shifted(state,du,dp,dn,dc,.5)
    rhs=m.evaluate(mid,kinetic_budget=state,kinetic_dt_ps=h)
    gp=h*rhs.dp_per_ps-dp;gn=h*rhs.dn_per_ps-dn
    # Direct physical excitation capacities, without using returned residuals.
    f=m.functional
    native_energy=np.sqrt(f.xi_bar[None,:]**2+abs(mid.delta_bar[:,None])**2)
    nodal_scale=f.N0_per_J_m3*f.Eu_J**2*f.volumes_m3
    electron_J=float(np.sum(4*nodal_scale[:,None]*f.weights_bar*native_energy*abs(gp)))
    phonon_J=float(np.sum(nodal_scale[:,None]*m.kinetics.phonons.capacities
        *m.kinetics.phonons.energies*abs(gn)))
    assert receipt['electron_residual_L1_energy_J']==pytest.approx(electron_J,rel=2e-6,abs=1e-35)
    assert receipt['phonon_residual_L1_energy_J']==pytest.approx(phonon_J,rel=2e-6,abs=1e-35)
    assert electron_J<tol*retained and phonon_J<tol*retained
    assert electron_J>=abs(float(np.sum(4*nodal_scale[:,None]*f.weights_bar*native_energy*gp)))
    # A residual endpoint-population correction enters the midpoint at gp/2.
    derivative=f.direction(mid.delta_bar,mid.p,np.zeros_like(mid.delta_bar),p_direction=gp/2)
    force=derivative.gap_force_direction_J_per_gapbar
    force_cart=np.column_stack((force.real,force.imag)).ravel()
    induced=h*(m.mobility_blocks(mid.delta_bar,rhs.temperature_K)@force_cart)
    guard=float(np.max(abs(induced[m._free_cart]),initial=0.))
    assert receipt['induced_gap_defect']==pytest.approx(guard,rel=2e-6,abs=1e-15)
    assert guard<tol
    gu=h*rhs.gap_coordinate_per_ps-du
    gc=h*rhs.circuit_scaled_per_ps-dc
    score=max(float(np.max(abs(gu))),float(np.max(abs(gc))),electron_J/retained,
        phonon_J/retained,guard)
    assert receipt['defect']==pytest.approx(score,rel=2e-6,abs=1e-15)
    assert score<tol
    assert np.all((end.p>=0)&(end.p<=1)) and np.all(end.n>=0)
    # Conservation is independent of the nonlinear acceptance score.
    final=m.evaluate(end)
    balance=final.energy_J-initial.energy_J-receipt['external_energy_increment_J']
    assert abs(balance)<1e-7*retained
