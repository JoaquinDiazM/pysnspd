"""Independent finite-action temperature and boundary-work checks."""
import unittest

import numpy as np
from numpy.polynomial.legendre import leggauss

from pysnspd.experimental import thermal_spatial_usadel as thermal
from pysnspd.experimental.thermal_internal_energy import evaluate_internal_energy
from pysnspd.experimental.adiabatic_energy_observer import AdiabaticElectronicEnergy


EMPTY = np.empty(0, dtype=int)


class ThermalInternalEnergyTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.graph = thermal.rectangular_graph(4, 3, 1.2, .7)
        coordinates = cls.graph.coordinates_bar
        cls.gap = (.85+.06*np.cos(3*coordinates[:, 0]))*np.exp(
            1j*(.3*coordinates[:, 0]+.07*np.sin(4*coordinates[:, 1])))
        cls.alpha = .04*np.sin(np.arange(len(cls.graph.edges)))
        cls.t = .37
        cls.modes = 7

    def solve(self, gap, t, *, law='free', frozen_boundaries=None):
        solutions, derivatives = [], []
        fixed = EMPTY if law == 'free' else self.graph.boundary_nodes
        for n in range(self.modes):
            epsilon = 2*np.pi*t*(n+.5)
            if law == 'free':
                boundary, boundary_rate = np.empty(0, complex), np.empty(0, complex)
            elif law == 'bcs':
                boundary = gap[fixed]/epsilon
                boundary_rate = -boundary/t
            elif law == 'frozen':
                boundary = frozen_boundaries[n]
                boundary_rate = np.zeros(len(fixed), complex)
            else:
                raise AssertionError('unknown test boundary law')
            solution = thermal.solve_frequency(self.graph, gap, epsilon, self.alpha,
                fixed_nodes=fixed, fixed_u=boundary, tol=1e-11)
            solutions.append(solution)
            derivatives.append(boundary_rate)
        result = evaluate_internal_energy(self.graph, gap, t, solutions, self.alpha,
            energy_scale_J=2.3e-21, boundary_contract=law,
            boundary_temperature_derivatives=derivatives)
        return result, tuple(solutions), tuple(derivatives)

    def test_uniform_free_graph_matches_analytic_finite_sum_and_temperature_tangent(self):
        graph = thermal.rectangular_graph(3, 2, 1.1, .8)
        gap = np.full(graph.n_nodes, .91*np.exp(.3j))
        solutions = tuple(thermal.solve_frequency(graph, gap, 2*np.pi*self.t*(n+.5),
            fixed_nodes=EMPTY, fixed_u=np.empty(0, complex), tol=1e-11) for n in range(self.modes))
        result = evaluate_internal_energy(graph, gap, self.t, solutions,
            energy_scale_J=1.9e-21, boundary_contract='natural zero flux')
        eps = 2*np.pi*self.t*(np.arange(self.modes)+.5)
        root = np.hypot(eps, abs(gap[0]))
        area = graph.area_weights.sum()
        a2 = abs(gap[0])**2
        spectral = 2*eps-2*root
        f_expected = area*(a2*np.log(self.t)+2*np.pi*self.t*np.sum(spectral+a2/eps))
        df_expected = area*(a2/self.t+2*np.pi*np.sum(
            spectral+2*eps*(1-eps/root)))
        self.assertAlmostEqual(result.free_energy_bar, f_expected, delta=5e-14)
        self.assertAlmostEqual(result.temperature_derivative_bar, df_expected, delta=1e-12)
        self.assertAlmostEqual(result.internal_energy_bar, f_expected-self.t*df_expected, delta=5e-13)
        for solution, tangent in zip(solutions, result.temperature_tangents):
            np.testing.assert_allclose(tangent, -solution.u/self.t, rtol=4e-12, atol=1e-13)
        self.assertEqual(result.boundary_reaction_temperature_derivative_bar, 0.)
        self.assertAlmostEqual(result.internal_energy_J,
            result.internal_energy_bar*result.energy_scale_J, delta=1e-34)

    def test_bcs_boundary_reactions_match_independent_temperature_difference(self):
        result, solutions, _ = self.solve(self.gap, self.t, law='bcs')
        dt = 2e-5
        before, before_solutions, _ = self.solve(self.gap, self.t-dt, law='bcs')
        after, after_solutions, _ = self.solve(self.gap, self.t+dt, law='bcs')
        measured = (after.free_energy_bar-before.free_energy_bar)/(2*dt)
        self.assertAlmostEqual(result.temperature_derivative_bar, measured, delta=3e-8)
        self.assertGreater(abs(result.boundary_reaction_temperature_derivative_bar), 1e-5)
        # Dropping the boundary reaction produces an observable derivative defect.
        self.assertGreater(abs(result.explicit_temperature_derivative_bar-measured), 1e-5)
        for tangent, low, high in zip(result.temperature_tangents, before_solutions, after_solutions):
            np.testing.assert_allclose(tangent, (high.u-low.u)/(2*dt), rtol=3e-8, atol=5e-9)
        self.assertLess(result.maximum_free_tangent_residual, 2e-14)
        self.assertLess(abs(result.free_stationarity_temperature_derivative_bar), 2e-10)

    def test_externally_frozen_boundaries_are_distinct_from_bcs_temperature_law(self):
        physical, solutions, _ = self.solve(self.gap, self.t, law='bcs')
        frozen = tuple(s.u[s.fixed_nodes].copy() for s in solutions)
        result, _, _ = self.solve(self.gap, self.t, law='frozen', frozen_boundaries=frozen)
        dt = 2e-5
        before, _, _ = self.solve(self.gap, self.t-dt, law='frozen', frozen_boundaries=frozen)
        after, _, _ = self.solve(self.gap, self.t+dt, law='frozen', frozen_boundaries=frozen)
        measured = (after.free_energy_bar-before.free_energy_bar)/(2*dt)
        self.assertAlmostEqual(result.temperature_derivative_bar, measured, delta=3e-8)
        self.assertEqual(result.boundary_reaction_temperature_derivative_bar, 0.)
        self.assertAlmostEqual(result.free_energy_bar, physical.free_energy_bar, delta=1e-14)
        self.assertGreater(abs(result.internal_energy_bar-physical.internal_energy_bar), 1e-5)

    def test_boundary_contract_and_derivatives_are_required_and_validated(self):
        _, solutions, derivatives = self.solve(self.gap, self.t, law='bcs')
        with self.assertRaisesRegex(ValueError, 'explicit temperature derivatives'):
            evaluate_internal_energy(self.graph, self.gap, self.t, solutions, self.alpha,
                energy_scale_J=1e-21, boundary_contract='no implicit law')
        with self.assertRaisesRegex(ValueError, 'boundary temperature contract'):
            evaluate_internal_energy(self.graph, self.gap, self.t, solutions, self.alpha,
                energy_scale_J=1e-21, boundary_contract='',
                boundary_temperature_derivatives=derivatives)
        with self.assertRaisesRegex(ValueError, 'one boundary temperature derivative vector'):
            evaluate_internal_energy(self.graph, self.gap, self.t, solutions, self.alpha,
                energy_scale_J=1e-21, boundary_contract='bcs',
                boundary_temperature_derivatives=derivatives[:-1])
        with self.assertRaisesRegex(ValueError, 'per prescribed spectral boundary'):
            evaluate_internal_energy(self.graph, self.gap, self.t, solutions, self.alpha,
                energy_scale_J=1e-21, boundary_contract='bcs',
                boundary_temperature_derivatives=tuple(np.zeros(1) for _ in solutions))

    def test_energy_scale_rejects_bool_complex_nonfinite_and_nonscalar_values(self):
        _, solutions, derivatives = self.solve(self.gap, self.t, law='bcs')
        for scale in (True, np.bool_(True), 1+0j, 1+2j, float('nan'),
                      float('inf'), 0., -1., '1', np.array(1.), 10**400):
            with self.subTest(scale_type=type(scale).__name__):
                with self.assertRaisesRegex(ValueError, 'finite positive real'):
                    evaluate_internal_energy(self.graph, self.gap, self.t, solutions,
                        self.alpha, energy_scale_J=scale, boundary_contract='bcs',
                        boundary_temperature_derivatives=derivatives)
        for scale in (np.float64(1e-21), np.float32(1e-21), 1):
            result = evaluate_internal_energy(self.graph, self.gap, self.t, solutions,
                self.alpha, energy_scale_J=scale, boundary_contract='bcs',
                boundary_temperature_derivatives=derivatives)
            self.assertEqual(result.energy_scale_J, float(scale))

    def test_normal_state_reference_and_global_gauge_rotation(self):
        zero = np.zeros(self.graph.n_nodes, complex)
        reference, _, _ = self.solve(zero, self.t)
        self.assertEqual(reference.free_energy_bar, 0.)
        self.assertEqual(reference.internal_energy_bar, 0.)
        first, _, _ = self.solve(self.gap, self.t, law='bcs')
        rotated, _, _ = self.solve(self.gap*np.exp(.47j), self.t, law='bcs')
        self.assertAlmostEqual(first.free_energy_bar, rotated.free_energy_bar, delta=2e-13)
        self.assertAlmostEqual(first.internal_energy_bar, rotated.internal_energy_bar, delta=2e-13)

    def test_off_shell_gap_direction_connects_to_nonthermal_energy_observer(self):
        # Synthetic dimensionless quadrature and SI scale, not admitted NbN.
        graph = thermal.rectangular_graph(2, 2, 1., 1.)
        x, w = leggauss(96)
        energy = 1.3+1.4*(x+1)/2
        weights = 1.4*w/2
        h0 = np.tanh(energy/(2*self.t))
        hL = h0[:, None]-.03*np.exp(-((energy[:, None]-2.)/.3)**2)*np.ones((1, graph.n_nodes))
        observer = AdiabaticElectronicEnergy(energy, weights, graph.area_weights, .7, h0)

        def snapshot(amplitude):
            gap = np.full(graph.n_nodes, amplitude*np.exp(.2j))
            solutions = tuple(thermal.solve_frequency(graph, gap, 2*np.pi*self.t*(n+.5),
                fixed_nodes=EMPTY, fixed_u=np.empty(0, complex), tol=1e-11) for n in range(self.modes))
            reference = evaluate_internal_energy(graph, gap, self.t, solutions,
                energy_scale_J=.7, boundary_contract='natural zero flux')
            rho = (energy/np.sqrt(energy**2-amplitude**2))[:, None]*np.ones((1, graph.n_nodes))
            return reference, rho, observer.energy(reference.internal_energy_J, rho, hL)

        amplitude, step = .9, 1e-5
        reference, rho, actual = snapshot(amplitude)
        low, _, low_energy = snapshot(amplitude-step)
        high, _, high_energy = snapshot(amplitude+step)
        thermal_direction = (high.internal_energy_J-low.internal_energy_J)/(2*step)
        rho_direction = (energy*amplitude/(energy**2-amplitude**2)**1.5)[:, None]*np.ones((1, graph.n_nodes))
        rate = observer.rate(thermal_direction, rho, hL, rho_direction, np.zeros_like(hL))
        measured = (high_energy.total_J-low_energy.total_J)/(2*step)
        self.assertAlmostEqual(rate.total_rate_W, measured, delta=2e-10)
        self.assertGreater(abs(rate.moving_dos_rate_W), 1e-3)
        omitted = observer.rate(thermal_direction, rho, hL, np.zeros_like(rho), np.zeros_like(hL))
        self.assertGreater(abs(omitted.total_rate_W-measured), 1e-3)
        # Temperature is a fixed bath reference in this gap direction.
        self.assertEqual(reference.temperature_ratio, self.t)
        self.assertEqual(actual.thermal_internal_J, reference.internal_energy_J)


if __name__ == '__main__':
    unittest.main()
