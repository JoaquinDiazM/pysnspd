"""Static spatial contracts, with independent BCS and phase-curvature oracles."""
from pathlib import Path

import numpy as np
import pytest

from pysnspd.experimental.energy_catalog import (
    OccupationEnergyCatalog, BCS_GAP_RATIO, HBAR_J_S, K_B_J_K, E_CHARGE_C)
from pysnspd.experimental.refined_cells import refined_count_catalog
from pysnspd.experimental.spatial_functional import (
    PeriodicSpatialFunctional, SpatialAdmissibilityError)


@pytest.fixture(scope="module")
def catalog():
    path = (Path(__file__).resolve().parents[1]/
            "docs/implementation/stage1_r2/catalogs/occupation_catalog.npz")
    return OccupationEnergyCatalog.load(path)


def model(catalog, cells=5, length_bar=8.):
    tc = catalog.vacuum.metadata["Tc_K"]
    ell = np.sqrt(HBAR_J_S*catalog.vacuum.D_m2_s/(2*K_B_J_K*tc))
    return PeriodicSpatialFunctional(catalog, length_bar*ell, 120e-9*7e-9, cells)


def populations(catalog, cells):
    x = catalog.count_nodes
    return np.array([(.04+.002*i)*np.exp(-x/.35)
                     +.035*np.exp(-((x-.95)/.24)**2) for i in range(cells)])


def complex_gradient(result):
    return result.gradient_cartesian_bar[:, 0]+1j*result.gradient_cartesian_bar[:, 1]


def fd5(function, step=2e-5):
    return (function(-2*step)-8*function(-step)+8*function(step)-function(2*step))/(12*step)


def test_scaling_and_uniform_BCS_vacuum_are_independent(catalog):
    system = model(catalog)
    a = .83
    p = np.zeros((system.cells, len(catalog.count_nodes)))
    result = system.evaluate(np.full(system.cells, a), p)
    expected = system.length_bar*a*a*(np.log(a)-.5)
    assert result.energy_bar == pytest.approx(expected, rel=2e-14)
    np.testing.assert_allclose(complex_gradient(result), system.h_bar*2*a*np.log(a),
                               rtol=2e-14, atol=2e-14)
    np.testing.assert_allclose(result.wirtinger_force_bar, a*np.log(a), atol=2e-14)
    np.testing.assert_array_equal(result.link_current_bar, np.zeros(system.cells))
    assert system.kappa == np.pi/4
    assert system.gap_ratio == pytest.approx(BCS_GAP_RATIO)
    expected_scale = (catalog.vacuum.N0_per_J_m3*catalog.vacuum.delta0_J**2
                      *system.cross_section_m2*system.ell0_m)
    assert result.energy_J == pytest.approx(expected_scale*expected)
    with pytest.raises(ValueError, match="BCS"):
        PeriodicSpatialFunctional(catalog, 360e-9, 120e-9*7e-9, 8,
                                  Tc_K=2*system.Tc_K)


def test_helical_energy_and_current_match_closed_form_discrete_oracle(catalog):
    system = model(catalog, cells=8, length_bar=12.)
    amplitude, q = .87, .13
    h, r = system.h_bar, system.gap_ratio
    alpha = q*h
    t = alpha/2
    a = amplitude*np.cos(t)
    k = 2*amplitude*np.sin(t)/h
    ap, kp = -amplitude*np.sin(t)/2, amplitude*np.cos(t)/h
    s = a*a+.1**2
    qd = a*k/s
    qp = ((ap*k+a*kp)*s-a*k*2*a*ap)/s**2
    gamma = qd*qd/r
    assert gamma < a  # Closed gapped formula below, independent of the engine.
    uvac = a*a*(np.log(a)-.5)+np.pi*a*gamma/2-gamma*gamma/3
    ua = 2*a*np.log(a)+np.pi*gamma/2
    ug = np.pi*a/2-2*gamma/3
    density = uvac+np.pi/4*(k*k-a*a*qd*qd)
    j = h*(ua*ap+ug*2*qd*qp/r+np.pi/4*(2*k*kp-2*a*ap*qd*qd-2*a*a*qd*qp))
    p = np.zeros((system.cells, len(catalog.count_nodes)))
    z = amplitude*np.exp(1j*q*h*np.arange(system.cells))
    result = system.evaluate(z, p, twist=q*system.length_bar)
    assert result.energy_bar == pytest.approx(system.length_bar*density, abs=2e-13)
    np.testing.assert_allclose(result.link_current_bar, j, rtol=2e-13, atol=2e-13)
    np.testing.assert_allclose(result.gamma_links_bar, gamma, rtol=2e-14)
    np.testing.assert_allclose(result.current_A,
                               (2*E_CHARGE_C/HBAR_J_S)*system.energy_scale_J*j)
    equivalent = system.evaluate(np.full(system.cells, amplitude), p,
                                 link_phases=np.full(system.cells, alpha))
    assert result.energy_bar == pytest.approx(equivalent.energy_bar, abs=2e-13)
    np.testing.assert_allclose(result.current_A, equivalent.current_A, atol=1e-19)


@pytest.mark.parametrize("use_complementary", [False, True])
def test_cartesian_force_and_link_current_differentiate_energy_at_fixed_p(catalog, use_complementary):
    actual = refined_count_catalog(catalog) if use_complementary else catalog
    system = model(actual, cells=4, length_bar=7.)
    z = np.array([.83+.02j, .89+.09j, .86+.18j, .82+.22j])
    p = populations(actual, 4)
    original_p = p.copy()
    links = np.array([.01, -.02, .015, .31])
    result = system.evaluate(z, p, link_phases=links, require_stability=False)
    dz = np.array([.2+.3j, -.1+.2j, .1-.3j, -.15+.1j])
    dl = np.array([.2, -.3, .15, .1])
    derivative = fd5(lambda t: system.evaluate(
        z+t*dz, p, link_phases=links+t*dl, require_stability=False).energy_bar)
    analytic = np.real(np.vdot(complex_gradient(result), dz))+np.dot(result.link_current_bar, dl)
    assert abs(derivative-analytic)/max(1., abs(analytic)) < 2e-7
    # A second, radial direction verifies the user-facing amplitude derivative.
    da = np.array([.2, -.1, .15, -.13])
    radial = fd5(lambda t: system.evaluate(
        z+t*da*z/abs(z), p, link_phases=links, require_stability=False).energy_bar)
    assert abs(radial-np.dot(result.amplitude_gradient_bar, da)) < 2e-7
    np.testing.assert_array_equal(p, original_p)


def test_gauge_covariance_and_discrete_current_divergence(catalog):
    system = model(catalog)
    z = np.array([.81+.03j, .88+.05j, .87+.12j, .84+.19j, .8+.2j])
    links = np.array([.01, -.015, .02, -.005, .25])
    chi = np.array([.3, -.2, .18, .7, -.45])
    p = populations(catalog, system.cells)
    a = system.evaluate(z, p, link_phases=links, require_stability=False)
    b = system.evaluate(z*np.exp(1j*chi), p,
                        link_phases=links+chi-np.roll(chi, -1), require_stability=False)
    assert a.energy_bar == pytest.approx(b.energy_bar, abs=2e-13)
    np.testing.assert_allclose(complex_gradient(b), complex_gradient(a)*np.exp(1j*chi), atol=3e-12)
    np.testing.assert_allclose(a.link_current_bar, b.link_current_bar, atol=3e-12)
    np.testing.assert_allclose(np.imag(np.conj(z)*complex_gradient(a)),
                               np.roll(a.link_current_bar, 1)-a.link_current_bar, atol=3e-14)


@pytest.mark.parametrize("q,stable", [(1., False), (.65, True)])
def test_D36_against_independent_documented_phase_curvature(catalog, q, stable):
    system = model(catalog)
    amplitude = .6
    p = np.zeros(len(catalog.count_nodes))
    symbol = system.principal_symbol(amplitude, 1j*amplitude*q, p)
    # C.8.1 is in kBTc units: convert its phase curvature to the Cartesian
    # gradient Hessian in Delta0 units, rather than comparing unlike factors.
    r = BCS_GAP_RATIO
    d = amplitude*r
    m = amplitude**2/(amplitude**2+.1**2)
    qd = m*q
    phase_curvature_kBTc = m*m*(np.pi*d-4*qd*qd)+np.pi/2*d*d*(1-m*m)
    phase_eigenvalue = phase_curvature_kBTc/(r*r*amplitude*amplitude)
    np.testing.assert_allclose(symbol.matrix, np.diag([np.pi/2, phase_eigenvalue]), atol=3e-10)
    assert symbol.stable is stable
    assert symbol.uncertainty < 1e-7
    if not stable:
        assert phase_curvature_kBTc == pytest.approx(-.343431, abs=6e-7)


def test_negative_or_unresolved_D36_is_not_an_evolution_permission(catalog):
    system = model(catalog, cells=8, length_bar=8.)
    # Realize the exact local negative-control amplitude and q on every link.
    alpha = 2*np.arctan(system.h_bar/2)
    z = np.full(system.cells, .6/np.cos(alpha/2))
    p = np.zeros((system.cells, len(catalog.count_nodes)))
    with pytest.raises(SpatialAdmissibilityError) as caught:
        system.evaluate(z, p, link_phases=np.full(system.cells, alpha))
    assert caught.value.symbol.eigenvalues[0] < 0
    diagnosis = system.evaluate(z, p, link_phases=np.full(system.cells, alpha),
                                require_stability=False)
    assert diagnosis.principal_symbols is None
    # Analytic zero of the gapped phase curvature tests the uncertain-sign gate.
    a, r, k = .6, system.gap_ratio, system.kappa
    rho, s = a*a, a*a+.01
    qd2 = r*r/4*(np.pi*a/r-2*k*rho+2*k*s*s/rho)
    q = np.sqrt(qd2)*s/rho
    uncertain = system.principal_symbol(a, 1j*a*q, p[0])
    assert abs(uncertain.eigenvalues[0]) < 1e-8
    assert not uncertain.stable


def test_principal_symbol_is_local_gradient_hessian_and_rotates(catalog):
    system = model(catalog)
    p = populations(catalog, 1)[0]
    z, d = .83+.11j, .03+.16j
    symbol = system.principal_symbol(z, d, p)
    measured = np.column_stack([
        [(fd5(lambda t: system.local_density(z, d+t*direction, p).gradient_derivative_bar).real),
         (fd5(lambda t: system.local_density(z, d+t*direction, p).gradient_derivative_bar).imag)]
        for direction in (1., 1j)])
    np.testing.assert_allclose(symbol.matrix, measured, atol=2e-7)
    angle = .7
    rotation = np.array([[np.cos(angle), -np.sin(angle)], [np.sin(angle), np.cos(angle)]])
    other = system.principal_symbol(z*np.exp(1j*angle), d*np.exp(1j*angle), p)
    np.testing.assert_allclose(other.matrix, rotation@symbol.matrix@rotation.T, atol=2e-9)


def test_exact_normal_point_and_no_low_amplitude_extrapolation(catalog):
    system = model(catalog)
    p = populations(catalog, system.cells)
    normal = system.evaluate(np.zeros(system.cells), p)
    expected = 4*system.h_bar*np.sum(p*(catalog.count_weights*catalog.count_nodes)[None, :])
    assert normal.energy_bar == pytest.approx(expected, rel=2e-14)
    np.testing.assert_array_equal(normal.gradient_cartesian_bar, np.zeros((system.cells, 2)))
    assert normal.amplitude_gradient_bar is None
    for symbol in normal.principal_symbols:
        np.testing.assert_array_equal(symbol.matrix, (np.pi/2)*np.eye(2))
    with pytest.raises(ValueError, match="amplitude"):
        system.evaluate(np.full(system.cells, .04), p, require_stability=False)
    with pytest.raises(ValueError, match="Gamma"):
        system.local_density(.8, 10j, p[0])


@pytest.mark.parametrize("bad", [-1e-16, 1.+1e-15, np.nan, np.inf])
def test_occupation_guards_and_shape_are_not_clipped(catalog, bad):
    system = model(catalog)
    p = np.zeros((system.cells, len(catalog.count_nodes)))
    p[2, 0] = bad
    with pytest.raises(ValueError, match="p_links"):
        system.evaluate(np.full(system.cells, .8), p, require_stability=False)


def test_progress_callback_and_vacuum_shortcut_never_call_spectral_kernel(catalog, monkeypatch):
    system = model(catalog)
    def forbidden(*args, **kwargs):
        raise AssertionError("an empty population does not require spectral kernels")
    monkeypatch.setattr(type(catalog), "energy_kernel", forbidden)
    events = []
    result = system.evaluate(np.full(system.cells, .8),
                             np.zeros((system.cells, len(catalog.count_nodes))),
                             on_cell=lambda i, n: events.append((i, n)))
    assert events == [(i, system.cells) for i in range(system.cells)]
    assert result.principal_symbols is not None


def test_amplitude_gradient_has_the_expected_laplacian_limit(catalog):
    errors = []
    for cells in (16, 32, 64):
        system = model(catalog, cells=cells, length_bar=12.)
        x = system.h_bar*np.arange(cells)
        wavenumber = 2*np.pi/system.length_bar
        a = .9+.03*np.cos(wavenumber*x)
        result = system.evaluate(a, np.zeros((cells, len(catalog.count_nodes))),
                                 require_stability=False)
        mid = (a+np.roll(a, -1))/2
        u_prime_mid = 2*mid*np.log(mid)
        gradient_only = (result.amplitude_gradient_bar/system.h_bar
                         -(u_prime_mid+np.roll(u_prime_mid, 1))/2)
        reference = 2*(np.pi/4)*.03*wavenumber**2*np.cos(wavenumber*x)
        errors.append(np.max(abs(gradient_only-reference)))
    assert errors[0]/errors[1] > 3.9
    assert errors[1]/errors[2] > 3.9
