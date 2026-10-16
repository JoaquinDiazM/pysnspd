"""Regression of signed heat against Pauli limits and the second law."""
import numpy as np
import pytest
from scipy.constants import Boltzmann as kB
from pysnspd.experimental.thermal_stiffness import entropy_compatible_carrier
from test_thermal_stiffness_coupled_candidate import small_system

def data(m,s):
    f=m.functional
    T=m.kinetics.equivalent_temperature_K(s.delta_bar,s.p)
    en=f.energies_bar(s.delta_bar)
    eps=4*f.node_energy_scale_J[:,None]*f.weights_bar*en
    return T,en,eps

@pytest.mark.parametrize("sign",[-1.,0.,1.])
def test_FD_limit_is_exact_thermal_tangent_for_either_power_sign(sign):
    m,s=small_system();T,en,eps=data(m,s);f=m.functional
    power=np.full(len(s.p),sign*2e-11)
    h=entropy_compatible_carrier(s.p,eps,T,f.Eu_J,en,power_W=power)
    fd=f.fermi_dirac(s.delta_bar,T)
    tangent=f.Eu_J*en/(kB*T[:,None]**2)*fd*(1-fd)
    expected=tangent/np.sum(eps*tangent,axis=1)[:,None]
    np.testing.assert_allclose(h,expected,rtol=2e-13,atol=1e-100)

@pytest.mark.parametrize("sign",[-1.,1.])
def test_nonthermal_source_conserves_energy_and_produces_nonnegative_excess_entropy(sign):
    m,s=small_system();f=m.functional
    s.p[1]=np.linspace(.025,.29,s.p.shape[1])[::-1]
    T,en,eps=data(m,s);power=np.full(len(s.p),sign*2e-11)
    h=entropy_compatible_carrier(s.p,eps,T,f.Eu_J,en,power_W=power)
    sigma=kB*eps/(f.Eu_J*en)*(np.log1p(-s.p)-np.log(s.p))
    heat=h*power[:,None]
    np.testing.assert_allclose(np.sum(eps*heat,axis=1),power,rtol=3e-15)
    excess=np.sum(sigma*heat,axis=1)-power/T
    assert np.all(excess>=-5e-14*np.abs(power/T))
    assert excess[1]>1e-17
    assert np.all(h>=0)

def test_closed_Pauli_limits_receive_only_inward_heat_without_log_clipping():
    m,s=small_system();f=m.functional
    s.p[1,0]=1.;s.p[1,-1]=0.
    before=s.p.copy();T,en,eps=data(m,s)
    for sign in (-1.,1.):
        power=np.full(len(s.p),sign*2e-11)
        h=entropy_compatible_carrier(s.p,eps,T,f.Eu_J,en,power_W=power)
        heat=h*power[:,None]
        assert heat[1,0]<=0 and heat[1,-1]>=0
        if sign>0:assert heat[1,0]==0.
        else:assert heat[1,-1]==0.
        np.testing.assert_allclose(np.sum(eps*heat,axis=1),power,rtol=3e-15)
        assert np.all(np.isfinite(h))
    np.testing.assert_array_equal(s.p,before)

def test_cooling_loss_is_proportional_to_available_population_even_in_subnormal_tail():
    m,s=small_system();f=m.functional;rates=[]
    for value in (1e-30,1e-100,1e-240):
        s.p[1,-1]=value;T,en,eps=data(m,s)
        power=np.full(len(s.p),-2e-11)
        h=entropy_compatible_carrier(s.p,eps,T,f.Eu_J,en,power_W=power)
        heat=h*power[:,None]
        rates.append(-heat[1,-1]/value)
        assert np.isfinite(rates[-1]) and rates[-1]>0
    assert max(rates)/min(rates)<1.000001
    s.p[1,-1]=np.nextafter(0.,1.);T,en,eps=data(m,s)
    h=entropy_compatible_carrier(s.p,eps,T,f.Eu_J,en,power_W=np.full(len(s.p),-2e-11))
    assert np.all(np.isfinite(h))
    np.testing.assert_allclose(np.sum(eps*h,axis=1),1.,rtol=3e-15)

def test_actual_continuous_RHS_points_into_Pauli_domain_at_exact_endpoints():
    m,s=small_system()
    s.p[1,0]=1.;s.p[1,-1]=0.
    r=m.evaluate(s)
    assert r.dp_per_ps[1,0]<=0.
    assert r.dp_per_ps[1,-1]>=0.
    assert np.all(np.isfinite(r.dp_per_ps))
    scale=abs(r.external_power_W)+1e-20
    assert abs(r.instant_power_residual_W)<1e-12*scale

@pytest.mark.parametrize("invalid",[-1e-12,1.00000001,np.nan])
def test_invalid_occupations_are_rejected_and_never_repaired(invalid):
    m,s=small_system();T,en,eps=data(m,s)
    s.p[1,0]=invalid
    with pytest.raises(ValueError):
        entropy_compatible_carrier(s.p,eps,T,m.functional.Eu_J,en,power_W=np.ones(len(s.p)))
    assert np.isnan(s.p[1,0]) if np.isnan(invalid) else s.p[1,0]==invalid
