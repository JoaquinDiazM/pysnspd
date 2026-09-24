import unittest
import numpy as np

from sandbox.stage4_core.harmonic_uniform_reference import (
    TAU3, uniform_reference, thermal_hessian_symbols, shifted_branches,
    modal_response, graph_laplacian, end_contacts, dirichlet_modes,
    harmonic_lift, project_field, korzh_scales, uniform_inductance_reference)
from pysnspd.experimental.thermal_spatial_usadel import rectangular_graph
from pysnspd.experimental.energy_catalog import HBAR_J_S, E_CHARGE_C


class HarmonicUniformReferenceTest(unittest.TestCase):
    def test_finite_sum_gap_scale_and_static_inductance_match_the_same_reference(self):
        reference = uniform_reference()
        self.assertAlmostEqual(reference.gap_bar, 1.76392592757718, delta=1e-12)
        self.assertLess(abs(reference.gap_equation_residual), 1e-13)
        radial, phase = thermal_hessian_symbols(np.array([0., .01, .3]), reference)
        self.assertEqual(phase[0], 0.)
        self.assertTrue(np.all(radial > 0))
        self.assertTrue(np.all(np.diff(phase) > 0))
        scales = korzh_scales()
        self.assertAlmostEqual(scales['graph_current_unit_A']/
            (scales['conductivity_S_m']*scales['thickness_m']),
            scales['E0_J']/(2*E_CHARGE_C), delta=1e-18)
        ell = scales['ell0_m']
        graph = rectangular_graph(9, 5, 160e-9/ell, 80e-9/ell)
        branch = uniform_inductance_reference(graph, reference)
        self.assertAlmostEqual(branch['geometric_conductance'], .5, delta=1e-13)
        expected = HBAR_J_S*scales['R_sheet_ohm']*2/(scales['E0_J']*reference.phase_stiffness_bar)
        self.assertAlmostEqual(branch['resolved_differential_H'], expected, delta=1e-22)
        self.assertAlmostEqual(branch['external_H']+branch['resolved_differential_H'], 10e-9, delta=1e-22)
        self.assertLess(branch['free_current_residual_bar'], 1e-14)
        # The infinite-sum London/BCS result is a convergence reference, not
        # silently substituted for this finite-Matsubara implementation.
        infinite_stiffness = np.pi*reference.gap_bar*np.tanh(reference.gap_bar/(2*reference.temperature_ratio))
        self.assertGreater(infinite_stiffness, reference.phase_stiffness_bar)
        self.assertLess((infinite_stiffness-reference.phase_stiffness_bar)/infinite_stiffness, .01)

    def test_analytic_branches_and_mixed_normal_limit_retain_finite_frequency(self):
        plus, minus = shifted_branches(.7, .6, 1.1+.3j, .25, kinds=('R', 'A'))
        for branch in (plus, minus):
            np.testing.assert_allclose(branch.matrix@branch.matrix, np.eye(2), atol=2e-15)
            np.testing.assert_allclose(branch.generator, branch.root*branch.matrix, atol=2e-15)
        for force in (np.array([[0., 1.], [1., 0.]]),
                      np.array([[0., 1j], [-1j, 0.]]), 1j*np.eye(2)):
            response = modal_response(plus, minus, .2, force)
            self.assertLess(response.normalization_residual_max, 5e-15)
            self.assertLess(response.equation_residual_max, 5e-15)
        normal_plus, normal_minus = shifted_branches(.7, .6, 0., .25, kinds=('R', 'A'))
        normal = modal_response(normal_plus, normal_minus, .2, 1j*np.eye(2))
        self.assertAlmostEqual(normal.denominator, .9-.6j)
        np.testing.assert_allclose(normal.matrix, 2j*np.eye(2)/(.9-.6j), atol=1e-15)
        np.testing.assert_allclose(normal_plus.matrix, TAU3, atol=1e-15)
        np.testing.assert_allclose(normal_minus.matrix, -TAU3, atol=1e-15)

    def test_graph_modes_preserve_natural_sides_and_report_omitted_fields(self):
        graph = rectangular_graph(9, 5, 8., 4.)
        left, right = end_contacts(graph)
        modes = dirichlet_modes(graph, 8)
        self.assertEqual(len(modes.fixed_nodes), 10)
        self.assertLess(len(modes.fixed_nodes), len(graph.boundary_nodes))
        self.assertLess(np.max(modes.relative_residuals), 1e-12)
        self.assertLess(modes.mass_orthogonality_error, 1e-13)
        np.testing.assert_allclose(modes.vectors[modes.fixed_nodes], 0., atol=0.)
        # Lowest mode is uniform across width, but the basis includes modes
        # with transverse variation and remains a two-dimensional discretization.
        rows = modes.vectors[:, 0].reshape(9, 5)
        self.assertLess(np.max(np.ptp(rows, axis=1)), 1e-14)
        self.assertGreater(np.max(np.ptp(modes.vectors[:, 2].reshape(9, 5), axis=1)), .1)
        fixed = np.r_[left, right]
        values = np.r_[np.full(len(left), -.3j), np.full(len(right), 1.+.4j)]
        lift = harmonic_lift(graph, fixed, values)
        expected = -.3j+(1.+.7j)*graph.coordinates_bar[:, 0]/8.
        np.testing.assert_allclose(lift, expected, atol=1e-14)
        lap = graph_laplacian(graph)
        self.assertLess(np.max(abs((lap@lift)[modes.free_nodes])), 1e-14)
        field = (1.+2j)*modes.vectors[:, 0]-.3*modes.vectors[:, 6]
        projection = project_field(graph, modes, field)
        self.assertLess(projection['omitted_mass_norm'], 1e-13)
        np.testing.assert_allclose(projection['reconstructed'], field, atol=1e-13)
        short = dirichlet_modes(graph, 1)
        self.assertAlmostEqual(project_field(graph, short, field)['omitted_mass_norm'], .3, delta=1e-13)
        with self.assertRaises(ValueError):
            project_field(graph, modes, field+lift)


if __name__ == '__main__':
    unittest.main()
