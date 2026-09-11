"""Independent checks for the opt-in uniform vacuum catalogue."""
import json

import numpy as np
import pytest
from scipy.integrate import quad
from scipy.optimize import brentq

from pysnspd.experimental.energy_catalog import (
    BCS_GAP_RATIO, UniformVacuumCatalog, build_vacuum_catalog,
    retarded_spectrum, vacuum_state, spectral_count, energy_at_count,
    build_occupation_catalog, OccupationEnergyCatalog,
)


def _catalog(n=33):
    return build_vacuum_catalog(np.linspace(.04, 1.6, n), np.linspace(0, 1.2, n),
                               Tc_K=8.65, N0_per_J_m3=1e47, D_m2_s=1.581e-4)


def _derivative(fun, x, h=1e-4):
    return (fun(x-2*h)-8*fun(x-h)+8*fun(x+h)-fun(x+2*h))/(12*h)


@pytest.mark.parametrize("delta,gamma", [(.4,0),(.8,.2),(.4,.4),(.4,.8),(.02,1.2)])
def test_vacuum_derivatives_and_independent_matsubara_integral(delta, gamma):
    value, force, conjugate = vacuum_state(delta, gamma)
    # Independent angle parameterization: omega=delta*cot(theta)-gamma*cos(theta).
    # It avoids both the closed-form expression and retarded polynomial roots.
    theta_max = np.pi/2 if gamma <= delta else np.arcsin(delta/gamma)
    reference = quad(lambda th: delta-gamma*np.sin(th)**3, 0, theta_max,
                     epsabs=1e-12, epsrel=1e-12)[0]
    assert conjugate == pytest.approx(reference, rel=1e-10, abs=1e-12)
    assert force == pytest.approx(_derivative(lambda d: vacuum_state(d,gamma)[0],delta),
                                  rel=1e-7, abs=1e-8)
    if gamma > .001:
        assert conjugate == pytest.approx(_derivative(lambda g: vacuum_state(delta,g)[0],gamma),
                                          rel=1e-7, abs=1e-8)
    if gamma == 0:
        assert value == pytest.approx(delta**2*(np.log(delta)-.5))


def test_normal_limit_and_gapless_small_amplitude():
    assert vacuum_state(0, 0) == (0,0,0)
    assert vacuum_state(0, 1) == (0,0,0)
    delta, gamma = 1e-10, .8
    value, force, conjugate = vacuum_state(delta,gamma)
    assert value/delta**2 == pytest.approx(np.log(2*gamma), rel=1e-12)
    assert force/(2*delta) == pytest.approx(np.log(2*gamma), rel=1e-12)
    assert conjugate/(delta**2/gamma) == pytest.approx(1, rel=1e-12)


def test_spline_uses_own_potential_and_preserves_maxwell():
    table = _catalog()
    for delta,gamma in ((.17,.1),(.55,.45),(1.13,.32),(.32,1.03)):
        _, force, conjugate = table.evaluate(delta,gamma)
        assert force == pytest.approx(_derivative(lambda d:table.evaluate(d,gamma)[0],delta),abs=1e-10)
        assert conjugate == pytest.approx(_derivative(lambda g:table.evaluate(delta,g)[0],gamma),abs=1e-10)
        mixed_d = _derivative(lambda d:table.evaluate(d,gamma)[2],delta)
        mixed_g = _derivative(lambda g:table.evaluate(delta,g)[1],gamma)
        assert mixed_d == pytest.approx(mixed_g, abs=1e-9)


def test_off_node_refinement_reduces_reference_errors():
    points = np.random.default_rng(8031).uniform([.1,.01],[1.5,1.15],size=(60,2))
    errors=[]
    for n in (17,33,65):
        table=_catalog(n)
        errors.append(np.max([np.abs(np.array(table.evaluate(*point))-vacuum_state(*point))
                              for point in points],axis=0))
    assert np.all(errors[-1] < errors[0])
    assert np.max(errors[-1]) < .003


def test_spectral_bcs_normal_and_original_equation():
    energy=np.unique(np.r_[0,np.linspace(.001,3,81),1])
    eta=1e-4
    for delta,gamma in ((0,0),(0,.5),(1,0),(1,.2),(1,1.3),(.1,1.2)):
        c,s=retarded_spectrum(energy,delta=delta,gamma=gamma,eta=eta)
        assert np.min(c.real)>=0
        lhs=delta*c
        rhs=(gamma*c-1j*(energy+1j*eta))*s
        assert np.max(np.abs(lhs-rhs)/np.maximum(1,np.abs(lhs)))<1e-8
        assert np.max(np.abs(c*c+s*s-1)/np.maximum(1,np.abs(c)**2+np.abs(s)**2))<1e-10
        if delta==0:
            assert np.array_equal(c,np.ones_like(c))
            assert np.array_equal(s,np.zeros_like(s))
        if gamma==0 and delta>0:
            z=energy+1j*eta
            assert np.allclose(c,z/np.sqrt(z*z-delta*delta),rtol=1e-12)


def test_eta_current_integral_has_independent_exact_reference():
    # The finite regulator has a known bias even with exact quadrature.
    eta=.003
    def integrand(energy):
        denominator=(energy*energy-eta*eta-1)**2+4*energy*energy*eta*eta
        return 2*energy*eta/denominator
    numerical=quad(integrand,0,1,epsabs=1e-11)[0]+quad(integrand,1,np.inf,epsabs=1e-11)[0]
    assert numerical==pytest.approx(np.pi/2-np.arctan(eta),abs=1e-10)
    assert abs(numerical-np.pi/2)>1e-3


@pytest.mark.parametrize("delta,gamma",[(np.nan,0),(-1,0),(1,-1),(1,np.inf)])
def test_invalid_vacuum_coordinates_rejected(delta,gamma):
    with pytest.raises(ValueError):
        vacuum_state(delta,gamma)


def test_support_rejects_extrapolation_and_normal_is_separate():
    table=_catalog()
    for point in ((.01,.1),(1.61,.1),(.1,1.21),(.1,-.01)):
        with pytest.raises(ValueError):
            table.evaluate(*point)
    assert table.evaluate(0, .2)==(0,0,0)
    with pytest.raises(ValueError):
        retarded_spectrum([0,1],delta=1,gamma=.1,eta=0)
    with pytest.raises(ValueError):
        retarded_spectrum([-1,1],delta=1,gamma=.1,eta=.001)


def test_npz_roundtrip_without_pickle_and_si_parity(tmp_path):
    table=_catalog()
    path=table.save(tmp_path/'catalog_without_extension')
    assert path.is_file()
    with np.load(path,allow_pickle=False) as arrays:
        assert json.loads(str(arrays['metadata_json']))['production_connected'] is False
    loaded=UniformVacuumCatalog.load(path)
    assert loaded.evaluate(.61,.18)==pytest.approx(table.evaluate(.61,.18))
    assert BCS_GAP_RATIO==pytest.approx(1.7638769888620456)
    positive=loaded.evaluate_si(.7*loaded.delta0_J,1e7)
    negative=loaded.evaluate_si(.7*loaded.delta0_J,-1e7)
    assert positive['Uvac_J_m3']==negative['Uvac_J_m3']
    assert positive['js_A_m2']==-negative['js_A_m2']


def test_count_primitive_inversion_and_spectral_derivatives():
    count=np.array([.002,.05,.2,1,3])
    for delta,gamma in ((.8,0),(.8,.2),(.2,.8)):
        energy=energy_at_count(count,delta=delta,gamma=gamma,eta=.001)
        assert np.allclose(spectral_count(energy,delta=delta,gamma=gamma,eta=.001),count,rtol=1e-8,atol=1e-10)
        if gamma==0:
            assert np.allclose(energy,count*np.sqrt(1+delta**2/(count**2+1e-6)))
        c,s=retarded_spectrum(energy,delta=delta,gamma=gamma,eta=.001)
        derivative=_derivative(lambda d:energy_at_count(count,delta=d,gamma=gamma,eta=.001),delta,h=1e-5)
        assert np.allclose(derivative,s.imag/c.real,rtol=1e-5,atol=1e-7)
        if gamma>0:
            derivative_g=_derivative(lambda g:energy_at_count(count,delta=delta,gamma=g,eta=.001),gamma,h=1e-5)
            assert np.allclose(derivative_g,-s.real*s.imag/c.real,rtol=1e-5,atol=1e-7)
    difficult_x=.01097964968413434
    energy=energy_at_count(difficult_x,delta=.72,gamma=1e-6,eta=.001)
    assert float(spectral_count(energy,delta=.72,gamma=1e-6,eta=.001))==pytest.approx(difficult_x,abs=1e-10)


def test_nonthermal_fixed_p_energy_derivatives_and_roundtrip(tmp_path):
    vacuum=build_vacuum_catalog(np.linspace(.15,1.4,9),np.linspace(0,1.2,9),
                               Tc_K=8.65,N0_per_J_m3=1e47,D_m2_s=1.581e-4)
    table=build_occupation_catalog(vacuum,count_order=24,count_max=5,eta=.003)
    occupation=.12*np.exp(-((table.count_nodes-1.1)/.35)**2)
    delta,gamma=.67,.43
    value,force,conjugate=table.evaluate(delta,gamma,occupation)
    assert force==pytest.approx(_derivative(lambda d:table.evaluate(d,gamma,occupation)[0],delta),abs=1e-9)
    assert conjugate==pytest.approx(_derivative(lambda g:table.evaluate(delta,g,occupation)[0],gamma),abs=1e-9)
    mixed_d=_derivative(lambda d:table.evaluate(d,gamma,occupation)[2],delta)
    mixed_g=_derivative(lambda g:table.evaluate(delta,g,occupation)[1],gamma)
    assert mixed_d==pytest.approx(mixed_g,abs=1e-8)
    path=table.save(tmp_path/'nonthermal')
    loaded=OccupationEnergyCatalog.load(path)
    assert np.array_equal(loaded.gamma_energy_derivatives,table.gamma_energy_derivatives)
    assert loaded.evaluate(delta,gamma,occupation)==pytest.approx((value,force,conjugate))
    normal=loaded.evaluate(0,.4,occupation)
    assert normal[0]==pytest.approx(4*np.dot(loaded.count_nodes*occupation,loaded.count_weights))
    assert normal[1:]==(0,0)
    # At a field node, the Hermite slope must reproduce its own physical
    # low-energy current kernel, including the Gamma=0 endpoint.
    delta_node=float(loaded.vacuum.delta_axis[3])
    p_low=.2*np.exp(-loaded.count_nodes/.25)
    c,s=retarded_spectrum(loaded.excitation_energies[3,0],delta=delta_node,gamma=0,eta=loaded.eta)
    expected=loaded.vacuum.evaluate(delta_node,0)[2]-4*np.dot(s.real*s.imag/c.real*p_low,loaded.count_weights)
    assert loaded.evaluate(delta_node,0,p_low)[2]==pytest.approx(expected,abs=1e-10)
    for invalid in (occupation*0+1.01,occupation*0-.1,np.ones(3),occupation*0+np.nan,occupation*0+np.inf):
        with pytest.raises(ValueError):
            loaded.evaluate(delta,gamma,invalid)
