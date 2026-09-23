"""Independent dispersion and BCS references for the nodal spatial energy."""
from pathlib import Path

import numpy as np
import pytest

from pysnspd.experimental.energy_catalog import OccupationEnergyCatalog, HBAR_J_S, K_B_J_K
from pysnspd.experimental.refined_cells import refined_count_catalog
from pysnspd.experimental.spatial_nodal import PeriodicNodalSpatialFunctional


@pytest.fixture(scope="module")
def catalog():
    return OccupationEnergyCatalog.load(Path(__file__).resolve().parents[1]/
        "docs/implementation/stage1_r2/catalogs/occupation_catalog.npz")


def system(catalog, cells=8, length_bar=12.):
    ell = np.sqrt(HBAR_J_S*catalog.vacuum.D_m2_s/
                  (2*K_B_J_K*catalog.vacuum.metadata["Tc_K"]))
    return PeriodicNodalSpatialFunctional(catalog, ell*length_bar, 120e-9*7e-9, cells)


def gradient(result):
    return result.gradient_cartesian_bar[:, 0]+1j*result.gradient_cartesian_bar[:, 1]


def fd5(function, h=2e-5):
    return (function(-2*h)-8*function(-h)+8*function(h)-function(2*h))/(12*h)


@pytest.mark.parametrize("phase_step", [.02, .4, 1.1, np.pi])
def test_dispersion_matches_independent_Fourier_symbols(catalog, phase_step):
    model = system(catalog)
    z = np.exp(1j*phase_step*np.arange(model.cells))
    twist = phase_step*model.cells
    h = model.h_bar
    wave = (8*np.sin(phase_step)-np.sin(2*phase_step))/(6*h)
    t = np.sin(phase_step/2)**2
    stiff = 4*t*(1+t/3)/h**2
    np.testing.assert_allclose(model.covariant_derivative(z, twist=twist), 1j*wave*z, atol=2e-15)
    np.testing.assert_allclose(model.positive_laplacian(z, twist=twist), stiff*z, atol=2e-15)
    assert stiff-wave*wave == pytest.approx(16*t**3*(2+t)/(9*h*h), abs=3e-15)
    if phase_step == np.pi:
        assert stiff == pytest.approx(16/(3*h*h))  # No checkerboard stiffness loss.


def test_unitary_gauge_shift_preserves_the_global_gradient_inequality(catalog):
    model = system(catalog, cells=9)
    rng = np.random.default_rng(7231)
    z = rng.normal(size=9)+1j*rng.normal(size=9)
    links = rng.normal(size=9)
    d = model.covariant_derivative(z, twist=.7, link_phases=links)
    lap = model.positive_laplacian(z, twist=.7, link_phases=links)
    assert np.real(np.vdot(z, lap))-np.vdot(d, d).real > 0
    rho = abs(z)**2
    q = np.imag(np.conj(z)*d)/(rho+.01)
    assert np.real(np.vdot(z, lap))-np.sum(rho*q*q) > 0
    # Independent check of the quadratic via its two positive/negative path terms.
    phases = links.copy(); phases[-1] += .7
    shift = np.exp(1j*phases)*np.roll(z, -1)
    shift2 = np.exp(1j*(phases+np.roll(phases, -1)))*np.roll(z, -2)
    expected = (4/3*np.sum(abs(shift-z)**2)-np.sum(abs(shift2-z)**2)/12)/model.h_bar**2
    assert np.real(np.vdot(z, lap)) == pytest.approx(expected, rel=2e-15)


def test_helical_BCS_energy_force_and_current_have_no_midpoint_amplitude_term(catalog):
    model = system(catalog)
    a, q0 = .9, .13
    h, r, kappa = model.h_bar, model.gap_ratio, model.kappa
    alpha = h*q0
    wave = (8*np.sin(alpha)-np.sin(2*alpha))/(6*h)
    wave_prime = (8*np.cos(alpha)-2*np.cos(2*alpha))/(6*h)
    stiff = (30-32*np.cos(alpha)+2*np.cos(2*alpha))/(12*h*h)
    stiff_prime = (32*np.sin(alpha)-4*np.sin(2*alpha))/(12*h*h)
    m = a*a/(a*a+.01)
    ma = 2*a*.01/(a*a+.01)**2
    q, qa, qp = m*wave, ma*wave, m*wave_prime
    gamma = q*q/r
    assert gamma < a
    u = a*a*(np.log(a)-.5)+np.pi*a*gamma/2-gamma*gamma/3
    ua = 2*a*np.log(a)+np.pi*gamma/2
    ug = np.pi*a/2-2*gamma/3
    expected_energy = model.length_bar*(u+kappa*a*a*(stiff-q*q))
    expected_force_density = ua+ug*2*q*qa/r+kappa*(2*a*(stiff-q*q)-2*a*a*q*qa)
    expected_current = h*(ug*2*q*qp/r+kappa*a*a*(stiff_prime-2*q*qp))
    z = a*np.exp(1j*alpha*np.arange(model.cells))
    p = np.zeros((model.cells, len(catalog.count_nodes)))
    result = model.evaluate(z, p, twist=q0*model.length_bar)
    assert result.energy_bar == pytest.approx(expected_energy, abs=2e-13)
    np.testing.assert_allclose(result.amplitude_nodes_bar, a, atol=2e-16)
    np.testing.assert_allclose(result.amplitude_gradient_bar/h, expected_force_density, atol=4e-14)
    np.testing.assert_allclose(result.link_current_bar, expected_current, atol=3e-14)
    np.testing.assert_allclose(result.gamma_nodes_bar, gamma, atol=2e-16)
    assert result.gradient_remainder_bar >= 0


@pytest.mark.parametrize("complementary", [False, True])
def test_nonthermal_force_and_all_original_face_currents_differentiate_one_energy(catalog, complementary):
    cat = refined_count_catalog(catalog) if complementary else catalog
    model = system(cat, cells=5, length_bar=9.)
    z = np.array([.86+.02j, .91+.06j, .88+.16j, .82+.24j, .79+.3j])
    links = np.array([.02, -.015, .008, -.01, .34])
    x = cat.count_nodes
    p = np.array([(.06+.003*i)*np.exp(-x/.3)+.02*np.exp(-((x-.9)/.2)**2) for i in range(5)])
    saved = p.copy()
    result = model.evaluate(z, p, link_phases=links, require_stability=False)
    dz = np.array([.2+.1j, -.1+.2j, -.15-.1j, .2-.15j, .1+.15j])
    dl = np.array([.3, -.2, .1, .2, -.15])
    measured = fd5(lambda t: model.evaluate(z+t*dz, p, link_phases=links+t*dl,
                                          require_stability=False).energy_bar)
    expected = np.real(np.vdot(gradient(result), dz))+np.dot(result.link_current_bar, dl)
    assert abs(measured-expected) < 3e-7
    # Gauge/phase identity exercises deposition of every two-hop path on both faces.
    np.testing.assert_allclose(np.imag(np.conj(z)*gradient(result)),
                               np.roll(result.link_current_bar, 1)-result.link_current_bar,
                               atol=2e-14)
    np.testing.assert_array_equal(p, saved)


def test_local_gauge_covariance_twist_equivalence_and_progress(catalog):
    model = system(catalog)
    a, alpha = .88, .2
    z = a*np.exp(1j*alpha*np.arange(model.cells))
    p = np.broadcast_to(.02*np.exp(-catalog.count_nodes/.4), (model.cells, len(catalog.count_nodes))).copy()
    events = []
    first = model.evaluate(z, p, twist=alpha*model.cells, require_stability=False,
                           on_cell=lambda i, n: events.append((i, n)))
    chi = np.array([.2, -.4, .1, .6, -.15, -.3, .25, .18])
    other = model.evaluate(z*np.exp(1j*chi), p, twist=alpha*model.cells,
                           link_phases=chi-np.roll(chi, -1), require_stability=False)
    assert first.energy_bar == pytest.approx(other.energy_bar, abs=5e-14)
    np.testing.assert_allclose(gradient(other), gradient(first)*np.exp(1j*chi), atol=5e-14)
    np.testing.assert_allclose(first.link_current_bar, other.link_current_bar, atol=2e-14)
    uniform = model.evaluate(np.full(model.cells, a), p, link_phases=np.full(model.cells, alpha),
                             require_stability=False)
    assert uniform.energy_bar == pytest.approx(first.energy_bar, abs=5e-14)
    assert events == [(i, model.cells) for i in range(model.cells)]
    with pytest.raises(ValueError, match="alias"):
        model.evaluate(z, p, on_node=lambda *_: None, on_cell=lambda *_: None)


def test_normal_and_zero_node_are_cartesian_but_unsupported_positive_core_is_rejected(catalog):
    model = system(catalog)
    p = np.broadcast_to(.05*np.exp(-catalog.count_nodes/.3), (model.cells, len(catalog.count_nodes))).copy()
    normal = model.evaluate(np.zeros(model.cells), p)
    expected = 4*model.length_bar*np.dot(catalog.count_weights*catalog.count_nodes, p[0])
    assert normal.energy_bar == pytest.approx(expected, rel=2e-15)
    np.testing.assert_array_equal(normal.gradient_cartesian_bar, np.zeros((model.cells, 2)))
    assert normal.amplitude_gradient_bar is None
    z = np.full(model.cells, .8); z[2] = 0
    isolated_zero = model.evaluate(z, p, require_stability=False)
    assert np.all(np.isfinite(isolated_zero.gradient_cartesian_bar))
    with pytest.raises(ValueError, match="amplitude"):
        model.evaluate(np.full(model.cells, .04), p, require_stability=False)


@pytest.mark.parametrize("bad", [-1e-16, 1+1e-15, np.nan, np.inf])
def test_population_guard_has_no_clipping(catalog, bad):
    model = system(catalog)
    p = np.zeros((model.cells, len(catalog.count_nodes))); p[0, 0] = bad
    with pytest.raises(ValueError, match="p_nodes"):
        model.evaluate(np.full(model.cells, .8), p, require_stability=False)


def test_vacuum_shortcut_and_fourth_order_amplitude_gradient(catalog, monkeypatch):
    def forbidden(*args, **kwargs):
        raise AssertionError("zero population must not request a spectral kernel")
    monkeypatch.setattr(type(catalog), "energy_kernel", forbidden)
    errors = []
    for n in (16, 32, 64):
        model = system(catalog, cells=n)
        x = model.h_bar*np.arange(n)
        k = 2*np.pi/model.length_bar
        a = .9+.03*np.cos(k*x)
        result = model.evaluate(a, np.zeros((n, len(catalog.count_nodes))), require_stability=False)
        # Pure real fields have q=0. The nodal local potential contributes
        # exactly 2a ln(a), with no midpoint averaging or phase contamination.
        measured = result.amplitude_gradient_bar/model.h_bar-2*a*np.log(a)
        expected = 2*model.kappa*.03*k*k*np.cos(k*x)
        errors.append(np.max(abs(measured-expected)))
    assert errors[0]/errors[1] > 15.5
    assert errors[1]/errors[2] > 15.5
