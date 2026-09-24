"""Analytic limits and graph consistency of intrinsic DC strip reservoirs."""
import unittest

import numpy as np
from scipy.constants import Boltzmann, elementary_charge, hbar

from pysnspd.experimental import bulk_current_reference as bulk
from pysnspd.experimental import thermal_spatial_usadel as thermal


class BulkCurrentReferenceTests(unittest.TestCase):
    t = .9/8.65
    count = 256
    material = dict(width_m=80e-9, Tc_K=8.65, diffusion_m2_s=.5e-4,
                    sheet_resistance_ohm=608.)

    def test_zero_current_bcs_and_finite_gap_equation(self):
        ref = bulk.solve_bulk_reference(self.t, self.count, 0.)
        exact_f = ref.gap_bar/np.hypot(ref.epsilon_bar, ref.gap_bar)
        np.testing.assert_allclose(ref.f, exact_f, atol=1e-15)
        np.testing.assert_allclose(ref.u, ref.gap_bar/ref.epsilon_bar, atol=0.)
        self.assertLess(abs(ref.gap_equation_residual), 1e-12)
        self.assertAlmostEqual(ref.current_density_bar, 0.)
        self.assertAlmostEqual(ref.gap_bar, 1.7639, delta=.001)
        # Full Matsubara BCS Lk/square = hbar*Rsheet/[pi*Delta*tanh(Delta/2T)].
        # Correct analytically for the intentionally finite current sum, rather
        # than mistaking its known omitted tail for a units discrepancy.
        actual = bulk.differential_inductance_per_length_H_m(ref, **self.material)
        exact = hbar*self.material['sheet_resistance_ohm']/(np.pi*ref.gap_bar
                *Boltzmann*self.material['Tc_K']*self.material['width_m']
                *np.tanh(ref.gap_bar/(2*self.t)))
        infinite_sum = ref.gap_bar/(4*self.t)*np.tanh(ref.gap_bar/(2*self.t))
        self.assertAlmostEqual(actual/exact, infinite_sum/np.sum(ref.f**2), delta=2e-14)

    def test_spectral_depairing_gap_suppression_and_reversal(self):
        zero = bulk.solve_bulk_reference(self.t, self.count, 0.)
        ref = bulk.solve_bulk_reference(self.t, self.count, .3)
        opposite = bulk.solve_bulk_reference(self.t, self.count, -.3)
        residual = ref.epsilon_bar*ref.u+ref.q_bar**2*ref.f-ref.gap_bar
        np.testing.assert_allclose(residual, 0., atol=4e-14)
        np.testing.assert_allclose(ref.f**2+ref.g**2, 1., atol=4e-16)
        self.assertLess(ref.gap_bar, zero.gap_bar)
        self.assertLess(abs(ref.gap_equation_residual), 1e-12)
        self.assertGreater(np.max(abs(ref.u-ref.gap_bar/ref.epsilon_bar)), .1)
        self.assertAlmostEqual(opposite.gap_bar, ref.gap_bar)
        self.assertAlmostEqual(opposite.current_density_bar, -ref.current_density_bar)
        with self.assertRaises(ValueError):
            ref.u[0] = 0.

    def test_current_target_selects_ascending_branch_and_rejects_overbias(self):
        peak = bulk.find_depairing_current(self.t, self.count)
        for target in (0., 15.5e-6, 21.5e-6, -15.5e-6):
            ref = bulk.solve_bulk_current(self.t, self.count, target, **self.material)
            self.assertAlmostEqual(ref.total_current_A(**self.material), target, delta=1e-14)
            self.assertLess(abs(ref.q_bar), peak.q_bar)
            self.assertGreater(bulk.current_density_derivative(ref), 0.)
        maximum = peak.total_current_A(**self.material)
        with self.assertRaisesRegex(ValueError, 'exceeds homogeneous'):
            bulk.solve_bulk_current(self.t, self.count, 1.001*maximum, **self.material)
        normal = bulk.solve_bulk_reference(self.t, self.count,
                                           1.001*bulk.normal_endpoint_q(self.t, self.count))
        self.assertEqual(normal.gap_bar, 0.)
        self.assertEqual(normal.current_density_bar, 0.)

    def test_current_derivative_and_differential_inductance(self):
        q, step = .33, 2e-5
        ref = bulk.solve_bulk_reference(self.t, self.count, q)
        plus = bulk.solve_bulk_reference(self.t, self.count, q+step)
        minus = bulk.solve_bulk_reference(self.t, self.count, q-step)
        derivative = (plus.current_density_bar-minus.current_density_bar)/(2*step)
        self.assertAlmostEqual(bulk.current_density_derivative(ref), derivative, delta=2e-8)
        scales = bulk.physical_scales(**{key: value for key, value in self.material.items()
                                        if key != 'width_m'})
        physical_derivative = (plus.total_current_A(**self.material)
                               -minus.total_current_A(**self.material))/(2*step)
        expected = hbar/(2*elementary_charge)/(scales['ell0_m']*physical_derivative)
        actual = bulk.differential_inductance_per_length_H_m(ref, **self.material)
        self.assertAlmostEqual(actual/expected, 1., delta=2e-8)

    def test_matched_reservoir_and_co_moving_gauge(self):
        graph = thermal.rectangular_graph(9, 3, 2., 1.)
        ref = bulk.solve_bulk_reference(self.t, self.count, .2)
        fields = bulk.intrinsic_boundary_fields(graph, ref)
        ids = np.array([0, 1, 2, 24, 25, 26])
        fixed = bulk.intrinsic_boundary_fields(graph, ref, nodes=ids)
        np.testing.assert_allclose(fixed['u'], fields['u'][:, ids])
        np.testing.assert_allclose(abs(fields['delta_bar']), ref.gap_bar)
        x = graph.coordinates_bar[:, 0]
        tail, head = graph.edges.T
        alpha = -ref.q_bar*(x[head]-x[tail])
        original = thermal.spectral_energy_gradient(graph, fields['delta_bar'],
                        ref.epsilon_bar[0], fields['u'][0])
        moved = thermal.spectral_energy_gradient(graph, np.full(graph.n_nodes, ref.gap_bar),
                        ref.epsilon_bar[0], np.full(graph.n_nodes, ref.u[0]), alpha)
        self.assertAlmostEqual(original.energy, moved.energy, delta=2e-14)
        np.testing.assert_allclose(moved.link_derivative_alpha,
                                   original.link_derivative_alpha, atol=3e-15)

    def test_cartesian_discrete_symbol_and_continuum_current_refinement(self):
        # Uniform f*exp(i q x) is an EXACT graph spectral state when its scalar
        # depairing uses the graph's sine symbol, with matching end reservoirs.
        # This independently checks q² and the factor in the integrated current.
        gap, ep, q, width = 1.4, .37, .35, 2.
        errors = []
        for nx in (9, 17):
            length = 4.
            graph = thermal.rectangular_graph(nx, 4, length, width)
            h = length/(nx-1)
            a = 2*(1-np.cos(q*h))/h**2
            u, f, _ = bulk.homogeneous_spectrum(gap, np.array([ep]), a)
            x = graph.coordinates_bar[:, 0]
            phase = np.exp(1j*q*x)
            value = thermal.spectral_energy_gradient(graph, gap*phase, ep, u[0]*phase)
            free = (x > x.min()) & (x < x.max())
            np.testing.assert_allclose(value.residual[free], 0., atol=4e-15)
            # One Matsubara prefactor; links crossing an interior transverse cut.
            tail, head = graph.edges.T
            cut = (x[tail] <= length/2) & (x[head] > length/2)
            graph_current = -2*np.pi*self.t*np.sum(value.link_derivative_alpha[cut])
            exact_graph_current = 4*np.pi*self.t*width*np.sin(q*h)/h*f[0]**2
            self.assertAlmostEqual(graph_current, exact_graph_current, delta=2e-15)
            continuum_f = bulk.homogeneous_spectrum(gap, np.array([ep]), q*q)[1][0]
            continuum_current = 4*np.pi*self.t*width*q*continuum_f**2
            errors.append(abs(graph_current-continuum_current))
        self.assertGreater(errors[0]/errors[1], 3.9)

    def test_invalid_material_and_cutoff_rejected(self):
        for arguments in ((0., 256, 0.), (1., 256, 0.), (.01, 1, 0.),
                          (.3, 3.5, 0.), (.3, 256, np.nan)):
            with self.assertRaises(ValueError):
                bulk.solve_bulk_reference(*arguments)
        with self.assertRaises(ValueError):
            bulk.solve_bulk_current(self.t, self.count, 1e-6,
                                   **dict(self.material, width_m=0.))


if __name__ == '__main__':
    unittest.main()
