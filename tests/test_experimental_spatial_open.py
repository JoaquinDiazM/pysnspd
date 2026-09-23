"""Fast algebraic contracts for the open SEM energy, using analytic data.

The small catalogue below is synthetic, differentiable and internally
consistent. These tests isolate assembly, chain rules, boundary signs and
thermal Legendre structure; they do not admit a material spectrum or transient.
The unittest form also runs without an installed pytest/conftest environment.
"""
from types import SimpleNamespace
from pathlib import Path
import importlib.util
import json
import tempfile
import unittest

import numpy as np

from pysnspd.experimental.energy_catalog import BCS_GAP_RATIO, HBAR_J_S, K_B_J_K, E_CHARGE_C
from pysnspd.experimental.spatial_open import OpenSpatialFunctional, lobatto_operators


class AnalyticCatalogue:
    """Four independent positive excitation energies with exact derivatives."""

    def __init__(self):
        self.count_nodes = np.array([.07, .31, .8, 1.7])
        self.count_weights = np.array([.12, .23, .31, .34])
        tc = 8.65
        self.vacuum = SimpleNamespace(metadata={'Tc_K': tc}, D_m2_s=1.58e-4,
            N0_per_J_m3=2.1e47, delta0_J=BCS_GAP_RATIO*K_B_J_K*tc,
            gamma_axis=np.array([0., 2.]), delta_axis=np.array([0., 2.]), evaluate=self.vacuum_values)

    @staticmethod
    def vacuum_values(a, gamma):
        if not 0 <= a <= 2 or not 0 <= gamma <= 2:
            raise ValueError('synthetic catalogue support')
        return (.5*(a*a-1)**2+.3*gamma*a*a+.07*gamma*gamma,
                2*a*(a*a-1)+.6*gamma*a, .3*a*a+.14*gamma)

    def energy_kernel(self, a, gamma):
        self.vacuum_values(a, gamma)
        root = np.sqrt(self.count_nodes**2+a*a)
        alpha = .13+.03*self.count_nodes
        return (root+alpha*gamma*(1+.2*a*a),
                a/root+.4*alpha*gamma*a,
                alpha*(1+.2*a*a))


def complex_gradient(result):
    return result.gradient_cartesian_bar[:, 0]+1j*result.gradient_cartesian_bar[:, 1]


def fd5(function, h=3e-5):
    return (function(-2*h)-8*function(-h)+8*function(h)-function(2*h))/(12*h)


class TestOpenSpatialOperators(unittest.TestCase):
    def test_lobatto_positive_mass_polynomial_exactness_and_sbp(self):
        for elements in (1, 2, 4):
            with self.subTest(elements=elements):
                x, mass, derivative, stiffness = lobatto_operators(elements, 2., 4)
                self.assertTrue(np.all(np.diff(x) > 0))
                self.assertTrue(np.all(mass > 0))
                self.assertAlmostEqual(mass.sum(), 2., places=14)
                for power in range(5):
                    expected = np.zeros_like(x) if power == 0 else power*x**(power-1)
                    np.testing.assert_allclose(derivative@x**power, expected, rtol=2e-13, atol=2e-13)
                for power in range(8):
                    self.assertAlmostEqual(np.dot(mass, x**power), 2.**(power+1)/(power+1), delta=2e-12)
                boundary = np.zeros_like(derivative)
                boundary[0, 0], boundary[-1, -1] = -1., 1.
                np.testing.assert_allclose(mass[:, None]*derivative+derivative.T*mass[None, :], boundary, atol=3e-14)
                np.testing.assert_allclose(stiffness, stiffness.T, atol=3e-14)

    def test_stiffness_and_compatible_remainder_are_psd_without_checkerboard(self):
        for elements in (1, 2, 4):
            with self.subTest(elements=elements):
                x, mass, derivative, stiffness = lobatto_operators(elements, 7., 4)
                scale = max(1., np.linalg.norm(stiffness, 2))
                eig = np.linalg.eigvalsh(stiffness)
                self.assertGreaterEqual(eig[0], -1e-13*scale)
                self.assertGreater(eig[1], 1e-3)
                np.testing.assert_allclose(stiffness@np.ones(len(x)), 0., atol=1e-13*scale)
                compatible = stiffness-derivative.T@(mass[:, None]*derivative)
                self.assertGreaterEqual(np.linalg.eigvalsh(compatible)[0], -1e-13*scale)
                checkerboard = (-1.)**np.arange(len(x))
                self.assertGreater(checkerboard@stiffness@checkerboard, 1.)
                if elements == 1:
                    np.testing.assert_allclose(compatible, 0., atol=1e-13*scale)
                else:
                    # The assembly remainder is nonzero at shared endpoints;
                    # averaging the two derivatives must not discard its cost.
                    self.assertGreater(np.linalg.norm(compatible, 2), .1)

    def test_bad_operator_geometry_is_rejected(self):
        for args in ((0, 1.), (1, 0.), (1, np.inf), (True, 1.), (1, 1., 1), (1, 1., 2.5)):
            with self.subTest(args=args), self.assertRaises(ValueError):
                lobatto_operators(*args)


class TestOpenSpatialVariations(unittest.TestCase):
    def setUp(self):
        self.catalog = AnalyticCatalogue()
        ell = np.sqrt(HBAR_J_S*self.catalog.vacuum.D_m2_s/(2*K_B_J_K*8.65))
        self.model = OpenSpatialFunctional(self.catalog, 12.*ell, 120e-9*7e-9, 2)
        x = self.model.x_bar
        self.z = (.88+.03*np.cos(2*np.pi*x/12))*np.exp(1j*(.035*x+.02*np.sin(2*np.pi*x/12)))
        self.links = .012*np.sin(np.arange(self.model.cells-1)+.4)
        self.p = np.array([(.055+.002*i)*np.exp(-self.catalog.count_nodes/.7)
                           for i in range(self.model.cells)])

    def evaluate(self, z=None, p=None, links=None):
        return self.model.evaluate(self.z if z is None else z,
            self.p if p is None else p, link_phases=self.links if links is None else links,
            require_stability=False)

    def test_every_cartesian_force_component_differentiates_the_fixed_population_energy(self):
        original_p, original_z = self.p.copy(), self.z.copy()
        result = self.evaluate()
        for node in range(self.model.cells):
            for component, direction in enumerate((1., 1j)):
                perturbation = np.zeros(self.model.cells, complex)
                perturbation[node] = direction
                measured = fd5(lambda t: self.evaluate(z=self.z+t*perturbation).energy_bar)
                self.assertAlmostEqual(measured, result.gradient_cartesian_bar[node, component], delta=2e-9)
        np.testing.assert_array_equal(self.p, original_p)
        np.testing.assert_array_equal(self.z, original_z)
        self.assertIsNone(result.free_energy_bar)

    def test_each_link_current_is_the_link_phase_derivative_including_endpoint_links(self):
        result = self.evaluate()
        for link in range(self.model.cells-1):
            perturbation = np.zeros(self.model.cells-1)
            perturbation[link] = 1.
            measured = fd5(lambda t: self.evaluate(links=self.links+t*perturbation).energy_bar)
            self.assertAlmostEqual(measured, result.link_current_bar[link], delta=2e-9)
        np.testing.assert_allclose(result.current_A,
            2*E_CHARGE_C/HBAR_J_S*self.model.energy_scale_J*result.link_current_bar, rtol=1e-15)
        phase_gradient = np.imag(np.conj(self.z)*complex_gradient(result))
        exterior_and_internal = np.r_[0., result.link_current_bar, 0.]
        np.testing.assert_allclose(phase_gradient, -np.diff(exterior_and_internal), atol=4e-14)
        self.assertAlmostEqual(result.global_phase_residual_bar, 0., delta=4e-14)

    def test_nonuniform_local_gauge_transformation_preserves_energy_currents_and_force(self):
        chi = .35*np.cos(np.arange(self.model.cells)*1.7)
        first = self.evaluate()
        transformed = self.evaluate(z=self.z*np.exp(1j*chi), links=self.links+chi[:-1]-chi[1:])
        self.assertAlmostEqual(first.energy_bar, transformed.energy_bar, delta=4e-14)
        np.testing.assert_allclose(complex_gradient(transformed), complex_gradient(first)*np.exp(1j*chi), atol=7e-14)
        np.testing.assert_allclose(first.link_current_bar, transformed.link_current_bar, atol=5e-14)
        np.testing.assert_allclose(self.model.covariant_derivative(self.z*np.exp(1j*chi), link_phases=self.links+chi[:-1]-chi[1:]),
            self.model.covariant_derivative(self.z, link_phases=self.links)*np.exp(1j*chi), atol=3e-15)

    def test_gauged_remainder_nonnegative_and_shared_nodes_are_not_double_counted(self):
        z, derivative, stiffness = self.model._operators(self.z, self.links)
        mass = self.model.mass_bar
        remainder = stiffness-derivative.conj().T@(mass[:, None]*derivative)
        self.assertGreaterEqual(np.linalg.eigvalsh(remainder)[0], -3e-13)
        result = self.evaluate()
        self.assertGreaterEqual(result.gradient_remainder_bar, 0.)
        self.assertAlmostEqual(mass.sum(), self.model.length_bar, places=13)
        self.assertEqual(self.model.cells, 2*4+1)

    def test_real_linear_field_has_only_outward_endpoint_gradient_flux(self):
        # Integration by parts for a''=0: K*a=(-a',0,...,0,+a').
        slope = .008
        amplitude = .8+slope*self.model.x_bar
        zero_p = np.zeros_like(self.p)
        result = self.model.evaluate(amplitude, zero_p, require_stability=False)
        local_force = 2*amplitude*(amplitude**2-1)
        residual = complex_gradient(result)-self.model.mass_bar*local_force
        expected = np.zeros(self.model.cells)
        expected[0], expected[-1] = -2*self.model.kappa*slope, 2*self.model.kappa*slope
        np.testing.assert_allclose(residual, expected, atol=7e-15)
        np.testing.assert_array_equal(result.link_current_bar, np.zeros(self.model.cells-1))
        self.assertAlmostEqual(result.gradient_remainder_bar, self.model.kappa*slope*slope*self.model.length_bar, delta=4e-15)

    def test_thermal_force_and_current_differentiate_free_energy_not_internal_energy(self):
        theta = .3
        result = self.model.evaluate_thermal(self.z, theta, link_phases=self.links, require_stability=False)
        dz = .2*np.cos(np.arange(self.model.cells))+.15j*np.sin(np.arange(self.model.cells)+.2)
        dl = .2*np.cos(np.arange(self.model.cells-1)+.3)
        def varied(t):
            return self.model.evaluate_thermal(self.z+t*dz, theta,
                link_phases=self.links+t*dl, require_stability=False)
        measured = fd5(lambda t: varied(t).free_energy_bar)
        expected = np.real(np.vdot(complex_gradient(result), dz))+np.dot(result.link_current_bar, dl)
        self.assertAlmostEqual(measured, expected, delta=2e-9)
        internal_derivative = fd5(lambda t: varied(t).energy_bar)
        self.assertGreater(abs(internal_derivative-expected), 1e-4)
        # Envelope identity: freeze FD occupations at the unvaried state.
        frozen = self.evaluate(p=result.p_nodes)
        np.testing.assert_allclose(frozen.gradient_cartesian_bar, result.gradient_cartesian_bar, atol=3e-14)
        np.testing.assert_allclose(frozen.link_current_bar, result.link_current_bar, atol=3e-14)
        self.assertAlmostEqual(frozen.energy_bar, result.energy_bar, delta=3e-14)
        entropy = -4*np.sum(self.model.mass_bar[:, None]*self.catalog.count_weights[None, :]*(
            result.p_nodes*np.log(result.p_nodes)+(1-result.p_nodes)*np.log1p(-result.p_nodes)))
        self.assertAlmostEqual(result.free_energy_bar, result.energy_bar-theta*entropy, delta=3e-14)

    def test_zero_temperature_and_vacuum_do_not_request_excitation_kernels(self):
        def forbidden(*args):
            raise AssertionError('vacuum/zero temperature must not request excitation kernels')
        self.catalog.energy_kernel = forbidden
        vacuum = self.evaluate(p=np.zeros_like(self.p))
        thermal = self.model.evaluate_thermal(self.z, 0., link_phases=self.links, require_stability=False)
        self.assertAlmostEqual(vacuum.energy_bar, thermal.free_energy_bar, delta=2e-15)
        np.testing.assert_array_equal(vacuum.gradient_cartesian_bar, thermal.gradient_cartesian_bar)
        np.testing.assert_array_equal(thermal.p_nodes, np.zeros_like(self.p))

    def test_stability_and_progress_cover_each_global_node_once(self):
        events = []
        result = self.model.evaluate(self.z, self.p, link_phases=self.links,
            on_node=lambda i, n: events.append((i, n)))
        self.assertEqual(events, [(i, self.model.cells) for i in range(self.model.cells)])
        self.assertEqual(len(result.principal_symbols), self.model.cells)
        self.assertTrue(all(symbol.stable for symbol in result.principal_symbols))

    def test_registered_open_solver_on_synthetic_catalogue_keeps_flux_and_gauge_scope(self):
        root = Path(__file__).resolve().parents[1]
        path = root/'sandbox/stage3_spatial/ports_20260923/run_open_batch.py'
        spec = importlib.util.spec_from_file_location('open_batch_synthetic_test', path)
        runner = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(runner)
        reg = json.loads((root/'docs/implementation/stage3/ports_20260923/open_registration.json').read_text(encoding='utf-8'))
        case = dict(reg['pilot'], id='synthetic_open_L360_E4')
        progress = SimpleNamespace(checkpoint=lambda *args, **kwargs: None)
        with tempfile.TemporaryDirectory(prefix='pysnspd-open-synthetic-') as folder:
            record = runner.solve_case(self.catalog, reg, case, progress, Path(folder))
            self.assertEqual(record['status'], 'PASS_STATIC_OPEN_CASE')
            self.assertFalse(record['solver']['active_bounds'])
            self.assertEqual(record['phase_rad'][-1], 0.)
            self.assertGreater(abs(record['phase_rad'][0]), .5)  # Left phase was not clamped to zero.
            np.testing.assert_allclose(record['amplitude_bar'][[0, -1]], record['reference_amplitude_bar'], atol=1e-15)
            self.assertFalse(record['stage3_closed'])
            self.assertFalse(record['production'])
            self.assertTrue((Path(folder)/'iterations.json').is_file())
            self.assertTrue((Path(folder)/'reservoir.json').is_file())

    def test_invalid_fields_populations_links_and_temperature_are_rejected(self):
        for value in (-1e-14, 1+1e-14, np.nan, np.inf, 1j):
            p = self.p.astype(complex) if isinstance(value, complex) else self.p.copy()
            p[0, 0] = value
            with self.subTest(population=value), self.assertRaises(ValueError):
                self.evaluate(p=p)
        for links in (np.zeros(self.model.cells), np.full(self.model.cells-1, np.inf), 1j*np.ones(self.model.cells-1)):
            with self.subTest(links=links), self.assertRaises(ValueError):
                self.evaluate(links=links)
        for theta in (-1., np.inf, np.nan):
            with self.subTest(theta=theta), self.assertRaises(ValueError):
                self.model.evaluate_thermal(self.z, theta)
        bad = self.z.copy(); bad[0] = np.nan
        with self.assertRaises(ValueError):
            self.evaluate(z=bad)


if __name__ == '__main__':
    unittest.main()
