"""Independent identities of the effective native BCS kinetic closure."""
from dataclasses import replace
import numpy as np
import pytest
from scipy.constants import Boltzmann as kB, hbar
from pysnspd.experimental.effective_bcs_kinetics import (
    EffectiveBCSKinetics,forward_euler_time_limit_s,acoustic_phonon_grid)
from pysnspd.experimental.kinetic_events import PhononGrid,ElectronPhononEvents
from pysnspd.experimental.effective_material import allmaras_acoustic_material


def model(**kw):
    xi=np.linspace(.03,2.,16);w=np.full(16,.12)
    ph=PhononGrid(np.r_[1e-9,np.linspace(.02,5.,64)],
        np.full(65,.04),lambda q:.002*q*q,(0.,5.),'synthetic numerical test')
    defaults=dict(xi_bar=xi,weights_bar=w,phonons=ph,Eu_J=1.2e-22,
        N0_per_J_m3=9e46,D_m2_s=5e-5,volumes_m3=np.array([1.,2.,3.])*1e-23,
        edges=np.array([[0,1],[1,2]]),face_areas_m2=np.array([1.,.7])*1e-16,
        edge_lengths_m=np.array([12.,17.])*1e-9,
        bath_temperature_K=.18*1.2e-22/kB,tau_ee_Tc_s=5e-12,
        tau_escape_s=15e-12,transport_order=48)
    defaults.update(kw)
    return EffectiveBCSKinetics(**defaults)


def initial(m,theta=(.35,.18,.27)):
    d=np.array([.45,.85,.65],complex)
    p=m.fermi_dirac(d,np.array(theta)*m.Eu_J/kB)
    n=np.tile(m.bath_phonons,(3,1))
    return d,p,n


def test_fd_be_stationarity_at_unequal_gaps_and_volumes():
    m=model();d,p,n=initial(m,(.18,.18,.18))
    r=m.evaluate(d,p,n)
    assert np.max(abs(r.dp_dt_s))<.03
    assert np.max(abs(r.dn_dt_s))<.03
    np.testing.assert_allclose(r.temperature_K,m.bath_temperature_K,rtol=2e-14)
    assert abs(r.balance_residual_W)<2e-24


def test_vectorized_collision_matches_independent_shared_event_implementation():
    m=model();d,p,n=initial(m)
    n[0]*=1.4
    r=m.evaluate(d,p,n)
    for k in range(3):
        e=r.electron_energies_J[k]/m.Eu_J
        independent=ElectronPhononEvents.from_arrays(e,m.weights_bar,
            abs(d[k])/e,m.phonons,rate_prefactor=8*np.pi*m.Eu_J/hbar,
            label='Independent scalar BCS event oracle')
        dp,dn=independent.rhs(p[k],n[k])
        np.testing.assert_allclose(r.collision_dp_dt_s[k],dp,rtol=3e-10,atol=.003)
        np.testing.assert_allclose(r.collision_dn_dt_s[k],dn,rtol=3e-10,atol=.003)


def test_native_energy_moments_conserve_every_channel_independently():
    m=model();d,p,n=initial(m)
    p[0]+=.03*np.exp(-((m.xi_bar-.5)/.2)**2)
    n[0]+=.02*np.exp(-((m.phonons.energies-1.2)/.3)**2)
    heat=np.array([2e-11,0.,1e-11])
    r=m.evaluate(d,p,n,heat)
    scale=np.sum(abs(r.collision_electron_power_W))+np.sum(abs(r.collision_phonon_power_W))
    assert np.max(abs(r.ep_energy_residual_W))<2e-14*scale
    ts=np.sum(abs(r.transport_face_power_W))
    assert abs(r.transport_energy_residual_W)<3e-14*ts
    assert abs(r.transport_number_residual_per_s)<.1
    assert np.max(abs(r.bgk_energy_residual_W))<1e-22
    assert np.max(abs(r.heating_residual_W))<1e-24
    assert abs(r.balance_residual_W)<1e-22
    independent=np.sum(r.electron_population_power_W+r.phonon_power_W+
        r.escape_power_W-heat)
    assert r.balance_residual_W==independent
    h=min(1e-18,.001*r.forward_euler_time_limit_s)
    finite=(m.energy_J(d,p+h*r.dp_dt_s,n+h*r.dn_dt_s)-m.energy_J(d,p,n))/h
    assert finite==pytest.approx(np.sum(r.electron_population_power_W+r.phonon_power_W),rel=1e-7,abs=1e-16)


def test_transport_direction_and_unequal_volume_scaling():
    m=model();d,p,n=initial(m)
    r=m.evaluate(d,p,n)
    assert r.transport_face_power_W[0]>0
    assert r.transport_face_power_W[1]<0
    doubled=model(volumes_m3=2*m.volumes_m3)
    other=doubled.evaluate(d,p,n)
    np.testing.assert_allclose(other.transport_dp_dt_s,r.transport_dp_dt_s/2,rtol=1e-14)
    np.testing.assert_array_equal(other.transport_face_power_W,r.transport_face_power_W)


def test_no_outward_vector_at_exhausted_native_support_and_no_mutation():
    m=model();d,p,n=initial(m)
    p[0,::7]=0.;p[1,::9]=1.;n[:,::3]=0.
    original_p=p.copy();original_n=n.copy()
    r=m.evaluate(d,p,n)
    assert np.all(r.dp_dt_s[p==0]>=0)
    assert np.all(r.dp_dt_s[p==1]<=0)
    assert np.all(r.dn_dt_s[n==0]>=0)
    assert r.forward_euler_time_limit_s>0
    np.testing.assert_array_equal(p,original_p);np.testing.assert_array_equal(n,original_n)


def test_tau_ee_calibrator_has_explicit_inverse_temperature_law():
    m=model(tau_escape_s=np.inf);d,p,n=initial(m)
    p[0]+=.01*np.exp(-((m.xi_bar-.6)/.2)**2)
    r=m.evaluate(d,p,n)
    fd=m.fermi_dirac(d,r.temperature_K)
    tau=m.tau_ee_Tc_s*m.Tc_K/np.maximum(r.temperature_K,m.bath_temperature_K)
    np.testing.assert_allclose(r.bgk_dp_dt_s,(fd-p)/tau[:,None],rtol=3e-14,atol=.002)
    assert m.metadata()['bgk_temperature_law']=='tau_ee(T)=tau_ee(Tc)*Tc/max(TE,Tb)'


def test_physical_heat_normalization_from_empty_electrons():
    m=model();d,p,n=initial(m);p[:]=0.
    source=np.array([3e-11,2e-11,0.])
    r=m.evaluate(d,p,n,source)
    assert np.all(r.heating_dp_dt_s>=0)
    np.testing.assert_allclose(r.heating_residual_W,0.,atol=1e-24)
    np.testing.assert_array_equal(r.temperature_K,np.zeros(3))


def test_support_bound_and_constitutive_rejection():
    assert forward_euler_time_limit_s(np.array([0.]),np.array([0.]),np.array([-1.]),np.array([0.]))==0.
    m=model();d,p,n=initial(m);p[0,0]=-1e-20
    with pytest.raises(ValueError):m.evaluate(d,p,n)


def test_acoustic_grid_uses_approved_material_without_fitted_multiplier():
    material=allmaras_acoustic_material()
    Eu=kB*8.65
    ph=acoustic_phonon_grid(material,Eu,nodes=96)
    assert len(ph.energies)==96
    assert ph.energies[-1]*Eu==pytest.approx(material.energy_J[-1])
    omega=np.array([.2,.8,2.])
    np.testing.assert_allclose(ph.coupling(omega),material.evaluate(omega*Eu)[0],rtol=1e-14)
    assert np.all(ph.capacities>0)


def test_shared_event_extent_step_preserves_support_and_energy_in_stiff_tails():
    m=model(tau_ee_Tc_s=np.inf,tau_escape_s=np.inf)
    d,p,n=initial(m)
    p[0,::7]=0.;p[1,::9]=1.;n[:,::3]=0.
    dt=1e-11
    r=m.evaluate_limited(d,p,n,budget_p=p,budget_n=n,dt_s=dt)
    pn=p+dt*r.dp_dt_s;nn=n+dt*r.dn_dt_s
    assert np.all((pn>=0)&(pn<=1))
    assert np.all(nn>=0)
    assert r.max_flux_reduction>0
    assert r.limited_event_energy_J>0
    assert 0<r.limited_event_fraction<=1
    assert r.limited_event_energy_J<=r.totalabs_event_energy_J
    assert abs(m.energy_J(d,pn,nn)-m.energy_J(d,p,n))<1e-32
    scale=np.sum(abs(r.collision_electron_power_W))+np.sum(abs(r.transport_face_power_W))
    assert abs(r.balance_residual_W)<1e-14*scale


def test_event_extent_correction_disappears_for_small_physical_step():
    m=model(tau_ee_Tc_s=np.inf,tau_escape_s=np.inf)
    d,p,n=initial(m)
    dt=1e-20
    raw=m.evaluate(d,p,n)
    limited=m.evaluate_limited(d,p,n,budget_p=p,budget_n=n,dt_s=dt)
    assert limited.max_flux_reduction==0.
    np.testing.assert_allclose(limited.dp_dt_s,raw.dp_dt_s,rtol=8e-16,atol=1e-7)
    np.testing.assert_allclose(limited.dn_dt_s,raw.dn_dt_s,rtol=8e-16,atol=1e-7)


def test_shared_event_limiter_does_not_rescale_supplied_heat():
    m=model();d,p,n=initial(m)
    source=np.array([2e-11,0.,1e-11])
    dt=1e-15
    r=m.evaluate_limited(d,p,n,budget_p=p,budget_n=n,dt_s=dt,heating_power_W=source)
    np.testing.assert_array_equal(r.heating_power_W,source)
    np.testing.assert_allclose(r.heating_residual_W,0.,atol=1e-24)
    assert np.all((p+dt*r.dp_dt_s>=0)&(p+dt*r.dp_dt_s<=1))
    assert np.all(n+dt*r.dn_dt_s>=0)
    assert abs(r.balance_residual_W)<1e-22


def test_reusing_actual_temperature_avoids_duplicate_inversion():
    m=model();d,p,n=initial(m)
    T=m.equivalent_temperature_K(d,p)
    a=m.evaluate(d,p,n,temperature_K=T)
    b=m.evaluate(d,p,n)
    np.testing.assert_allclose(a.dp_dt_s,b.dp_dt_s,rtol=1e-13,atol=.01)
    np.testing.assert_array_equal(a.collision_dp_dt_s,b.collision_dp_dt_s)


def test_midpoint_budget_allows_escaping_phonons_produced_in_current_step():
    m=model()
    d,p,n=initial(m)
    # Cold old tails, with newly emitted current-step high-energy phonons.
    oldn=n.copy();midn=n.copy()
    midn[:,-8:]+=1e-12
    dt=1e-15
    r=m.evaluate_limited(d,p,midn,budget_p=p,budget_n=oldn,dt_s=dt)
    assert np.all(np.isfinite(r.dn_dt_s))
    assert abs(r.balance_residual_W)<1e-22


def test_numba_matches_numpy_oracle_for_actual_acoustic_material():
    material=allmaras_acoustic_material()
    Eu=kB*8.65
    ph=acoustic_phonon_grid(material,Eu,nodes=48)
    fast=model(phonons=ph,Eu_J=Eu,accelerated=True)
    slow=model(phonons=ph,Eu_J=Eu,accelerated=False)
    d,p,n=initial(fast);dt=1e-15
    a=fast.evaluate_limited(d,p,n,budget_p=p,budget_n=n,dt_s=dt)
    b=slow.evaluate_limited(d,p,n,budget_p=p,budget_n=n,dt_s=dt)
    np.testing.assert_allclose(a.dp_dt_s,b.dp_dt_s,rtol=2e-13,atol=.003)
    np.testing.assert_allclose(a.dn_dt_s,b.dn_dt_s,rtol=2e-13,atol=.003)
    assert a.limited_event_energy_J==pytest.approx(b.limited_event_energy_J,rel=2e-13,abs=1e-300)


def test_equilibrium_identity_is_relative_and_reactivates_resolved_deviation():
    m=model()
    d,p,n=initial(m,(.18,.18,.18))
    r=m.evaluate(d,p,n)
    assert r.equilibrium_collision_nodes==3
    np.testing.assert_array_equal(r.collision_dp_dt_s,np.zeros_like(p))
    np.testing.assert_array_equal(r.collision_dn_dt_s,np.zeros_like(n))
    p[0,4]*=1+2e-12
    changed=m.evaluate(d,p,n)
    assert changed.equilibrium_collision_nodes==2
    assert np.max(abs(changed.collision_dp_dt_s[0]))>0


def test_equilibrium_identity_only_removes_roundoff_of_scalar_event_oracle():
    exact=model(recognize_equilibrium=True)
    raw=model(recognize_equilibrium=False)
    d,p,n=initial(exact,(.18,.18,.18))
    a=exact.evaluate(d,p,n);b=raw.evaluate(d,p,n)
    assert np.max(abs(a.collision_dp_dt_s-b.collision_dp_dt_s))<.003
    assert np.max(abs(a.collision_dn_dt_s-b.collision_dn_dt_s))<.003
    assert np.sum(abs(b.collision_electron_power_W))<1e-23


def test_compiled_parallel_node_rates_match_single_thread_and_keep_volume_energy():
    material=allmaras_acoustic_material()
    Eu=kB*8.65;ph=acoustic_phonon_grid(material,Eu,nodes=48)
    m=model(phonons=ph,Eu_J=Eu,accelerated=True)
    d,p,n=initial(m);n[0]*=1.3;dt=1e-15
    from pysnspd.experimental import effective_bcs_kinetics_numba as fast
    original=fast.get_num_threads()
    try:
        fast.set_num_threads(1)
        one=m.evaluate_limited(d,p,n,budget_p=p,budget_n=n,dt_s=dt)
        fast.set_num_threads(min(4,original))
        many=m.evaluate_limited(d,p,n,budget_p=p,budget_n=n,dt_s=dt)
    finally:fast.set_num_threads(original)
    np.testing.assert_array_equal(one.dp_dt_s,many.dp_dt_s)
    np.testing.assert_array_equal(one.dn_dt_s,many.dn_dt_s)
    assert one.limited_event_energy_J==pytest.approx(many.limited_event_energy_J,rel=1e-14,abs=1e-300)


@pytest.mark.parametrize("accelerated",[False,True])
def test_shared_event_fraction_is_continuous_at_availability_boundaries(accelerated):
    """One physical emission event crosses the formerly discontinuous budget."""
    m=model(accelerated=accelerated)
    d,p,n=initial(m)
    en=m.energies_bar(d)
    omega=np.where(m.recombination,en[:,m.event_i]+en[:,m.event_j],
                   en[:,m.event_j]-en[:,m.event_i])
    upper=np.searchsorted(m.phonons.energies,omega)
    lower=upper-1
    beta=(omega-m.phonons.energies[lower])/(m.phonons.energies[upper]-m.phonons.energies[lower])
    rate=np.zeros_like(omega)
    ev=int(np.flatnonzero(~m.recombination)[0])
    i,j=m.event_i[ev],m.event_j[ev]
    rate[0,ev]=1.
    events=dict(rate=rate,omega=omega,lower=lower,upper=upper,beta=beta)
    fractions=[]
    eps=1e-9
    ratios=[0.,.7,1.-eps,1.,1.+eps,1/.95-eps,1/.95,1/.95+eps,2.]
    for ratio in ratios:
        bp=np.full_like(p,10.);bh=np.full_like(p,10.);bn=np.full_like(n,10.)
        bp[0,j]=ratio
        dp,dn,_,_,_,metrics=m._limit_events(events,None,p.shape,n.shape,p,n,.1,bp,bh,bn)
        fraction=-dp[0,j]*4*m.weights_bar[j]
        fractions.append(fraction)
        # Every limited event keeps the SAME electron and phonon extent.
        assert dp[0,i]*4*m.weights_bar[i]==pytest.approx(fraction)
        moment=np.sum(4*m.weights_bar*en*dp,axis=1)+np.sum(m.phonons.capacities*m.phonons.energies*dn,axis=1)
        assert np.max(abs(moment))<2e-16
        assert metrics['max_flux_reduction']==pytest.approx(1-fraction)
    assert fractions[0]==0.
    assert fractions[1]==pytest.approx(.665)
    assert fractions[3]==pytest.approx(.95)
    assert fractions[-1]==1.
    assert abs(fractions[4]-fractions[2])<3*eps
    assert abs(fractions[7]-fractions[5])<3*eps
    assert np.all(np.diff(fractions)>=0)
