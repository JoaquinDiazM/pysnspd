
'''Independent energy, units, support and admissibility checks of TE candidate.'''
from dataclasses import replace
import numpy as np
import pytest
from scipy.constants import Boltzmann as kB,hbar,elementary_charge as e
from pysnspd.experimental.thermal_spatial_usadel import ThermalGraph
from pysnspd.experimental.thermal_stiffness import ThermalStiffnessFunctional,entropy_compatible_carrier
from pysnspd.experimental.thermal_coupled import ThermalCoupledSystem,ThermalCapacityError
from pysnspd.experimental.effective_coupled import EffectiveState
from pysnspd.experimental.effective_bcs_kinetics import EffectiveBCSKinetics
from pysnspd.experimental.kinetic_events import PhononGrid
from pysnspd.experimental.effective_bcs_ports import EffectiveNeutralPorts
from pysnspd.experimental.electrical_ports import ThesisCircuitParameters

def small_system(conductance=.1):
    Tc=8.65;Eu=kB*Tc;D=5e-5;df=7e-9;Rs=608.;ell=np.sqrt(hbar*D/(2*Eu))
    N0=1/(2*e**2*D*Rs*df)
    graph=ThermalGraph(np.array([4.,8.,4.]),np.array([[0,1],[1,2]]),
        np.full(2,conductance),np.array([[0.,0.],[4.,0.],[8.,0.]]),np.array([0,2]))
    xi=np.linspace(.03,2.,16);w=np.full(16,.12)
    f=ThermalStiffnessFunctional(graph,energy_scale_J=Eu,N0_per_J_m3=N0,ell0_m=ell,
        thickness_m=df,count_nodes=xi,count_weights=w,bath_temperature_K=.9,diffusion_m2_s=D)
    ph=PhononGrid(np.r_[1e-9,np.linspace(.02,5.,64)],np.full(65,.04),
        lambda q:.002*q*q,(0.,5.),"synthetic conservative test quadrature")
    k=EffectiveBCSKinetics.from_graph(graph,ell0_m=ell,thickness_m=df,xi_bar=xi,
        weights_bar=w,phonons=ph,Eu_J=Eu,N0_per_J_m3=N0,D_m2_s=D,
        bath_temperature_K=.9,transport_order=16,accelerated=False)
    ports=EffectiveNeutralPorts(graph,sheet_resistance_ohm=Rs,left_nodes=np.array([0]),right_nodes=np.array([2]))
    circuit=ThesisCircuitParameters(Lk_ext_H=1e-7,V_bias_V=.011)
    m=ThermalCoupledSystem(f,k,ports,circuit,current_scale_A=1e-6)
    d=np.array([1.7,1.65*np.exp(.07j),1.7*np.exp(.14j)],complex)
    p=k.fermi_dirac(d,np.array([.9,7.,.9]))
    n=np.tile(k.bath_phonons,(3,1));n[1]=k.bose(5.)
    return m,EffectiveState(d,p,n,np.array([1.1e-6,1e-6,2e-5]))

def independent_energy(m,s):
    return (m.functional.evaluate(s.delta_bar,s.p).energy_J+
        np.sum(m.kinetics.phonon_energy_weights_J*s.n)+m.circuit.energy(s.circuit_SI))

def test_continuous_balance_matches_full_internal_direction_and_cut_reservoir():
    m,s=small_system();r=m.evaluate(s);f=m.functional;k=m.kinetics
    dd=r.delta_per_ps*1e12;dp=r.dp_per_ps*1e12;dn=r.dn_per_ps*1e12
    cmrate=r.circuit_scaled_per_ps*m.circuit_scales*1e12
    internal=f.direction(s.delta_bar,s.p,dd,dp).energy_direction_J
    ph=float(np.sum(k.phonon_energy_weights_J*dn))
    Ib,Is,vc=s.circuit_SI
    circuit=float(np.dot([m.circuit.L_bias_H*Ib,m.circuit.Lk_ext_H*Is,m.circuit.C_couple_F*vc],cmrate))
    total=internal+ph+circuit
    scale=abs(total)+abs(r.external_power_W)+1e-25
    assert abs(total-r.external_power_W)<3e-13*scale
    assert abs(r.instant_power_residual_W)<3e-13*scale
    assert abs(r.film_power_residual_W)<3e-13*scale
    np.testing.assert_array_equal(r.dp_per_ps[m.cut_nodes],0.)
    np.testing.assert_array_equal(r.dn_per_ps[m.cut_nodes],0.)
    np.testing.assert_array_equal(r.material_delta_per_ps[m.cut_nodes],0.)
    assert abs(r.reservoir_power_out_W)>1e-16
    assert r.kinetics.electron_population_power_W.shape==(3,)
    assert r.kinetics.heating_residual_W.shape==(3,)
    np.testing.assert_allclose(r.kinetics.electron_population_power_W,
        np.sum(f.thermal_data(s.delta_bar,s.p,r.temperature_K)["epsilon_J"]*r.kinetics.dp_dt_s,axis=1),
        rtol=3e-15,atol=1e-28)

def test_total_energy_finite_difference_uses_no_reported_power_observer():
    m,s=small_system();r=m.evaluate(s)
    dd=r.delta_per_ps*1e12;dp=r.dp_per_ps*1e12;dn=r.dn_per_ps*1e12
    dc=r.circuit_scaled_per_ps*m.circuit_scales*1e12
    derivatives=[]
    for h in [2e-17,1e-17]:
        plus=EffectiveState(s.delta_bar+h*dd,s.p+h*dp,s.n+h*dn,s.circuit_SI+h*dc)
        minus=EffectiveState(s.delta_bar-h*dd,s.p-h*dp,s.n-h*dn,s.circuit_SI-h*dc)
        derivatives.append((independent_energy(m,plus)-independent_energy(m,minus))/(2*h))
    scale=abs(r.external_power_W)+1e-25
    assert abs(derivatives[-1]-r.external_power_W)<2e-6*scale
    assert abs(derivatives[-1]-derivatives[0])<2e-6*scale

def test_entropy_compatible_carrier_has_dual_moments_for_nonthermal_populations_and_signed_power():
    m,s=small_system();f=m.functional;k=m.kinetics
    T=k.equivalent_temperature_K(s.delta_bar,s.p)
    en=f.energies_bar(s.delta_bar);eps=4*f.node_energy_scale_J[:,None]*f.weights_bar*en
    perturb=.02*s.p[1,0]
    s.p[1,0]+=perturb;s.p[1,1]-=eps[1,0]/eps[1,1]*perturb
    T=k.equivalent_temperature_K(s.delta_bar,s.p)
    h=entropy_compatible_carrier(s.p,eps,T,f.Eu_J,en)
    sigma=kB*eps/(f.Eu_J*en)*(np.log1p(-s.p)-np.log(s.p))
    np.testing.assert_allclose(np.sum(eps*h,axis=1),1.,rtol=3e-15,atol=0.)
    np.testing.assert_allclose(np.sum(sigma*h,axis=1)*T,1.,rtol=3e-15,atol=0.)
    assert np.all(h>0)
    for power in [2e-11,-2e-11]:
        dp=h*power
        np.testing.assert_allclose(np.sum(eps*dp,axis=1),power,rtol=3e-15,atol=0.)
        np.testing.assert_allclose(np.sum(sigma*dp,axis=1),power/T,rtol=3e-15,atol=0.)
        # A sufficiently small signed source step preserves Pauli support.
        bound=np.min(np.minimum(s.p,1-s.p)/np.maximum(abs(dp),1e-300))
        trial=s.p+.1*bound*dp
        assert np.all((trial>0)&(trial<1))

def test_negative_capacity_is_constitutive_failure_without_floor_or_mutation():
    m,s=small_system(conductance=100.)
    s.delta_bar[1]=1.1*np.exp(1.3j);s.p[1]=m.kinetics.fermi_dirac(s.delta_bar,9.)[1]
    before=s.copy()
    td=m.functional.thermal_data(s.delta_bar,s.p)
    assert np.any(td["capacity_effective_J_K"]<=0)
    with pytest.raises(ThermalCapacityError,match="no capacity floor"):
        m.evaluate(s)
    np.testing.assert_array_equal(s.delta_bar,before.delta_bar)
    np.testing.assert_array_equal(s.p,before.p)

def test_nonfinite_population_is_rejected_without_energy_repair():
    m,s=small_system();s.p[1,0]=np.nan
    with pytest.raises(ValueError):
        m.evaluate(s)

def test_small_coupled_midpoint_preserves_support_and_runs_full_residual():
    m,s=small_system();result,info=m.midpoint_step(s,1e-4,tolerance=2e-8,max_iterations=60)
    assert np.all((result.p>0)&(result.p<1))
    assert np.all(result.n>=0)
    assert np.all(np.isfinite(result.delta_bar))
    assert np.all(result.circuit_SI[:2]>0)
    assert info["defect"]<2e-8
    assert abs(info["midpoint_power_residual_W"])<1e-21


def test_total_entropy_rate_matches_native_production_and_bath_ledger():
    m,s=small_system();r=m.evaluate(s);f=m.functional;k=m.kinetics
    en=f.energies_bar(s.delta_bar);eps=4*f.node_energy_scale_J[:,None]*f.weights_bar*en
    sigma=kB*eps/(f.Eu_J*en)*(np.log1p(-s.p)-np.log(s.p))
    phcounts=k.phonon_energy_weights_J/(f.Eu_J*k.phonons.energies[None,:])
    sigph=kB*phcounts*np.log1p(1/s.n)
    native_dp=r.kinetics.collision_dp_dt_s+r.kinetics.transport_dp_dt_s+r.kinetics.bgk_dp_dt_s
    native_dn=r.kinetics.collision_dn_dt_s+r.kinetics.escape_dn_dt_s
    expected=(np.sum(sigma*native_dp)+np.sum(sigph*native_dn)+
        np.sum((r.kwt_heat_W+r.electrical.nodal_joule_W)/r.temperature_K)-
        r.reservoir_power_out_W/k.bath_temperature_K)
    def entropy(state):
        sqp=-kB*np.sum(4*f.node_energy_scale_J[:,None]*f.weights_bar/f.Eu_J*
            (state.p*np.log(state.p)+(1-state.p)*np.log1p(-state.p)))
        sph=kB*np.sum(phcounts*(np.log1p(state.n)+state.n*np.log1p(1/state.n)))
        gradient=float(np.sum(f.evaluate(state.delta_bar,state.p).gradient_entropy_J_per_K))
        return sqp+sph+gradient
    h=1e-17
    dd=r.delta_per_ps*1e12;dp=r.dp_per_ps*1e12;dn=r.dn_per_ps*1e12
    plus=EffectiveState(s.delta_bar+h*dd,s.p+h*dp,s.n+h*dn,s.circuit_SI.copy())
    minus=EffectiveState(s.delta_bar-h*dd,s.p-h*dp,s.n-h*dn,s.circuit_SI.copy())
    finite=(entropy(plus)-entropy(minus))/(2*h)
    assert finite==pytest.approx(expected,rel=2e-6,abs=1e-24)
    total_with_baths=finite+(r.escape_power_W+r.reservoir_power_out_W)/k.bath_temperature_K
    assert total_with_baths>=0.


def test_population_gap_response_matches_admissible_finite_direction_including_ports_and_mobility():
    m,s=small_system();r=m.evaluate(s)
    direction=np.zeros_like(s.p);direction[1]=.005*s.p[1]*(1-s.p[1])
    analytic,dphi,dtemp=m.population_gap_response(s,r,direction)
    step=1e-4
    plus=s.copy();minus=s.copy();plus.p+=step*direction;minus.p-=step*direction
    rp=m.evaluate(plus);rm=m.evaluate(minus)
    finite=(rp.delta_per_ps-rm.delta_per_ps)/(2*step)
    scale=max(float(np.max(abs(finite))),1e-30)
    assert np.max(abs(analytic-finite))<2e-6*scale
    phi_finite=(rp.electrical.phi_V-rm.electrical.phi_V)/(2*step)
    np.testing.assert_allclose(dphi,phi_finite,rtol=2e-6,atol=1e-15)
    temp_finite=(rp.temperature_K-rm.temperature_K)/(2*step)
    np.testing.assert_allclose(dtemp,temp_finite,rtol=2e-6,atol=1e-11)

def test_signed_population_residual_can_cross_empty_tail_without_forming_trial_state(monkeypatch):
    m,s=small_system()
    s.p[1]=m.kinetics.fermi_dirac(s.delta_bar,1.)[1]
    r=m.evaluate(s);direction=np.zeros_like(s.p)
    direction[1,-1]=-2*s.p[1,-1]
    assert s.p[1,-1]+direction[1,-1]<0.
    def forbidden(*args,**kwargs):
        raise AssertionError("A signed residual direction must not invert a trial FD state")
    monkeypatch.setattr(m.functional,"equivalent_temperature_K",forbidden)
    monkeypatch.setattr(m.kinetics,"equivalent_temperature_K",forbidden)
    response,dphi,dtemp=m.population_gap_response(s,r,direction)
    assert np.all(np.isfinite(response))
    assert np.all(np.isfinite(dphi))
    assert np.all(np.isfinite(dtemp))
    assert dtemp[1]<0

def test_population_response_respects_one_sided_mobility_floor_at_bath():
    m,s=small_system()
    s.p[1]=m.kinetics.fermi_dirac(s.delta_bar,m.kinetics.bath_temperature_K)[1]
    r=m.evaluate(s);direction=np.zeros_like(s.p)
    direction[1]=.01*s.p[1]*(1-s.p[1])
    analytic,dphi,dtemp=m.population_gap_response(s,r,direction)
    assert dtemp[1]>0
    step=1e-4;plus=s.copy();plus.p+=step*direction
    finite=(m.evaluate(plus).delta_per_ps-r.delta_per_ps)/step
    scale=max(float(np.max(abs(finite))),1e-30)
    assert np.max(abs(analytic-finite))<5e-4*scale
