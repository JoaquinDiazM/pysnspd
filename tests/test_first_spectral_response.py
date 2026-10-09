"""Analytic moving-gap/gauge solutions and equivalence to the uncached solve."""
import unittest
import numpy as np
from dataclasses import replace
from pysnspd.experimental import adiabatic_matrix_usadel as ad
from pysnspd.experimental.thermal_spatial_usadel import rectangular_graph
from pysnspd.experimental.first_spectral_response import FirstSpectralResponse, ChartFirstSpectralResponse
from pysnspd.experimental import retarded_spatial_usadel as ret


def uniform_jet(nodes, d, energy, eta, d_time=0.):
    z = eta-1j*energy
    root = np.sqrt(z*z+d*d)
    tau1 = np.array([[0., 1.], [1., 0.]])
    R = (z*ad.TAU3+d*tau1)/root
    RE = -1j*ad.TAU3/root+1j*z*R/(root*root)
    Rt = d_time*(tau1/root-d*R/(root*root))
    return ad.Jet.of(np.broadcast_to(R, (nodes, 2, 2)), RE, Rt)


class FirstSpectralResponseTests(unittest.TestCase):
    def setUp(self):
        self.graph = rectangular_graph(4, 3, 2., 1.5)
        self.n = self.graph.n_nodes
        self.fixed = self.graph.boundary_nodes
        self.kappa = .4

    def test_uniform_gap_motion_has_zero_first_correction_and_reuses_factor(self):
        R = uniform_jet(self.n, 1.2, 2., 0.)
        B = ad.generator(np.full(self.n, 1.2), 2.)
        response = FirstSpectralResponse(self.graph, R, B,
            fixed_nodes=self.fixed, kappa_s=self.kappa)
        factor = response.factor
        self.assertIsNotNone(factor)
        for number, rate in enumerate((.01, -.07, .4), 1):
            moving = uniform_jet(self.n, 1.2, 2., 0., rate)
            Bt = ad.generator(np.full(self.n, 1.2), 2., d_time=rate)
            result, metrics = response.apply(moving, Bt,
                fixed_first=np.zeros((len(self.fixed), 2, 2)), tolerance=1e-13)
            np.testing.assert_allclose(result.first, 0., atol=2e-14)
            self.assertLess(metrics['spectral_first_max'], 2e-14)
            self.assertLess(metrics['normalization_first_all_max'], 2e-14)
            self.assertEqual(metrics['apply_calls'], number)
            self.assertEqual(metrics['fallback_calls'], 0)
            self.assertFalse(metrics['admission'])
            self.assertIs(response.factor, factor)

    def test_temporal_gauge_preserves_required_identity_correction(self):
        R = uniform_jet(self.n, 1.2, .8, .03)
        B = ad.generator(np.full(self.n, 1.2), .8, eta=.03)
        U = ad.links(np.full(self.n, -.31), np.full(self.n, -.17))
        target, Bg = ad.transport(U, R, self.kappa), ad.transport(U, B, self.kappa)
        input_R = replace(target, first=np.zeros_like(target.first))
        response = FirstSpectralResponse(self.graph, input_R, Bg,
            fixed_nodes=self.fixed, kappa_s=self.kappa)
        result, metrics = response.apply(input_R, Bg,
            fixed_first=target.first[self.fixed], tolerance=1e-13)
        np.testing.assert_allclose(result.first, target.first, atol=2e-11)
        self.assertGreater(metrics['identity_first_max'], .01)
        self.assertLess(metrics['spectral_first_max'], 1e-11)
        self.assertLess(metrics['normalization_first_all_max'], 1e-11)
        original, old = ad.solve_first_spectral(self.graph, input_R, Bg,
            fixed_nodes=self.fixed, fixed_first=target.first[self.fixed],
            kappa=self.kappa, tolerance=1e-13)
        np.testing.assert_allclose(result.first, original.first, atol=3e-11)

    def test_nonuniform_temporal_gauge_and_link_rates_match_analytic_solution(self):
        rng = np.random.default_rng(10926)
        R = uniform_jet(self.n, 1.2, 2., 0.)
        B = ad.generator(np.full(self.n, 1.2), 2.)
        chi, chi_time = rng.normal(size=(2, self.n))*.1
        U = ad.links(-chi, -chi_time)
        target, Bg = ad.transport(U, R, self.kappa), ad.transport(U, B, self.kappa)
        tail, head = self.graph.edges.T
        alpha = chi[head]-chi[tail]
        alpha_time = chi_time[head]-chi_time[tail]
        input_R = replace(target, first=np.zeros_like(target.first))
        response = FirstSpectralResponse(self.graph, input_R, Bg, alpha,
            self.fixed, kappa_s=self.kappa)
        result, metrics = response.apply(input_R, Bg, alpha_time_per_s=alpha_time,
            fixed_first=target.first[self.fixed], tolerance=1e-13)
        np.testing.assert_allclose(result.first, target.first, atol=1e-11)
        self.assertLess(metrics['spectral_first_max'], 1e-11)
        self.assertLess(metrics['normalization_first_all_max'], 1e-11)

    def test_explicit_original_lsmr_fallback_has_independent_residual(self):
        R = uniform_jet(self.n, 1.2, .8, .03)
        B = ad.generator(np.full(self.n, 1.2), .8, eta=.03)
        U = ad.links(np.zeros(self.n), np.full(self.n, -.11))
        target, Bg = ad.transport(U, R, self.kappa), ad.transport(U, B, self.kappa)
        input_R = replace(target, first=np.zeros_like(target.first))
        response = FirstSpectralResponse(self.graph, input_R, Bg,
            fixed_nodes=self.fixed, kappa_s=self.kappa, factor_normal_equations=False)
        result, metrics = response.apply(input_R, Bg,
            fixed_first=target.first[self.fixed], tolerance=1e-13)
        self.assertEqual(metrics['solver'], 'lsmr_original_rectangular_system')
        self.assertIsNotNone(metrics['lsmr_stop'])
        self.assertEqual(metrics['fallback_calls'], 1)
        np.testing.assert_allclose(result.first, target.first, atol=2e-11)

    def test_inaccurate_factor_answer_falls_back_to_original_system(self):
        R = uniform_jet(self.n, 1.2, .8, .03)
        B = ad.generator(np.full(self.n, 1.2), .8, eta=.03)
        U = ad.links(np.zeros(self.n), np.full(self.n, -.11))
        target, Bg = ad.transport(U, R, self.kappa), ad.transport(U, B, self.kappa)
        input_R = replace(target, first=np.zeros_like(target.first))
        response = FirstSpectralResponse(self.graph, input_R, Bg,
            fixed_nodes=self.fixed, kappa_s=self.kappa)
        class BadFactor:
            def solve(self, rhs):
                return np.zeros_like(rhs)
        response.factor = BadFactor()
        result, metrics = response.apply(input_R, Bg,
            fixed_first=target.first[self.fixed], tolerance=1e-13)
        self.assertEqual(metrics['solver'], 'lsmr_original_rectangular_system')
        np.testing.assert_allclose(result.first, target.first, atol=2e-11)

    def test_all_prescribed_reports_supplied_boundary_normalization_defect(self):
        R = uniform_jet(self.n, 1.2, 2., 0.)
        B = ad.generator(np.full(self.n, 1.2), 2.)
        response = FirstSpectralResponse(self.graph, R, B,
            fixed_nodes=np.arange(self.n), kappa_s=self.kappa)
        wrong = np.broadcast_to(ad.TAU3*.1, R.value.shape).copy()
        result, metrics = response.apply(R, B, fixed_first=wrong)
        np.testing.assert_array_equal(result.first, wrong)
        self.assertEqual(metrics['normalization_first_max'], 0.)
        self.assertGreater(metrics['normalization_first_all_max'], .1)
        self.assertFalse(metrics['admission'])

    def test_scalar_identity_generator_shift_preserves_cached_commutator(self):
        R = uniform_jet(self.n, 1.2, 2., 0., .07)
        B = ad.generator(np.full(self.n, 1.2), 2., d_time=.07)
        response = FirstSpectralResponse(self.graph, R, B,
            fixed_nodes=self.fixed, kappa_s=self.kappa)
        shifted = replace(B, value=B.value+.17j*ad.I2)
        result, metrics = response.apply(R, shifted,
            fixed_first=np.zeros((len(self.fixed), 2, 2)))
        self.assertLess(metrics['spectral_first_max'], 1e-12)

    def test_cached_operator_matches_full_matrix_direction_on_nonuniform_fields(self):
        rng = np.random.default_rng(5109)
        R = uniform_jet(self.n, 1.2, .8, .03)
        B = ad.generator(np.full(self.n, 1.2), .8, eta=.03)
        U = ad.links(rng.normal(size=self.n)*.3)
        R = ad.transport(U, R, self.kappa)
        B = ad.transport(U, B, self.kappa)
        alpha = rng.normal(size=len(self.graph.edges))*.1
        response = FirstSpectralResponse(self.graph, R, B, alpha, kappa_s=self.kappa)
        direction = rng.normal(size=R.value.shape)+1j*rng.normal(size=R.value.shape)
        base, norm0 = ad.spectral_residual(self.graph, R, B, alpha=alpha, kappa=self.kappa)
        moved, norm1 = ad.spectral_residual(self.graph, replace(R, first=direction), B,
            alpha=alpha, kappa=self.kappa)
        difference = np.concatenate((((moved.first-base.first)
            /self.graph.area_weights[:, None, None]).ravel(), (norm1.first-norm0.first).ravel()))
        np.testing.assert_allclose(response.full_matrix@direction.ravel(), difference,
            rtol=3e-14, atol=3e-14)
        self.assertGreater(response.estimated_storage_bytes, response.R0_value.nbytes)

    def test_stale_state_and_missing_cut_time_data_are_rejected(self):
        R = uniform_jet(self.n, 1.2, 2., 0.)
        B = ad.generator(np.full(self.n, 1.2), 2.)
        response = FirstSpectralResponse(self.graph, R, B,
            fixed_nodes=self.fixed, kappa_s=self.kappa)
        with self.assertRaisesRegex(ValueError, 'cuts'):
            response.apply(R, B)
        with self.assertRaisesRegex(ValueError, 'identical'):
            response.apply(replace(R, value=R.value+.001*ad.I2), B,
                fixed_first=np.zeros((len(self.fixed), 2, 2)))
        with self.assertRaisesRegex(ValueError, 'commutator'):
            response.apply(R, replace(B, value=B.value+.001*ad.TAU3),
                fixed_first=np.zeros((len(self.fixed), 2, 2)))


class ChartSpectralResponseTests(unittest.TestCase):
    def setUp(self):
        self.graph = rectangular_graph(4, 3, 2., 1.5)
        self.n, self.kappa = self.graph.n_nodes, .4
        self.fixed = self.graph.boundary_nodes

    def test_temporal_gauge_identity_matches_full_solver(self):
        R = uniform_jet(self.n, 1.2, .8, .03)
        B = ad.generator(np.full(self.n, 1.2), .8, eta=.03)
        U = ad.links(np.full(self.n, -.31), np.full(self.n, -.17))
        target, Bg = ad.transport(U, R, self.kappa), ad.transport(U, B, self.kappa)
        raw = replace(target, first=np.zeros_like(target.first))
        chart = ChartFirstSpectralResponse(self.graph, raw, Bg,
            fixed_nodes=self.fixed, kappa_s=self.kappa)
        result, metrics = chart.apply(raw, Bg, fixed_first=target.first[self.fixed], tolerance=1e-13)
        np.testing.assert_allclose(result.first, target.first, atol=3e-13)
        self.assertGreater(metrics['identity_first_max'], .01)
        self.assertLess(metrics['spectral_first_max'], 1e-12)
        self.assertLess(metrics['normalization_first_all_max'], 1e-12)
        self.assertFalse(hasattr(chart, 'full_matrix'))
        self.assertEqual(metrics['fallback_calls'], 0)

    def test_nonuniform_physical_gap_links_and_full_generator_source(self):
        rng = np.random.default_rng(1095109)
        gap = 1.2*np.exp(1j*rng.normal(size=self.n)*.1)
        gap[[4, 7]] *= .82
        alpha = rng.normal(size=len(self.graph.edges))*.06
        z = .03-.8j
        fixed_a, fixed_b = ret.uniform_fields(gap[self.fixed], z)
        root = ret.solve(self.graph, gap, z, alpha, fixed_nodes=self.fixed,
            fixed_a=fixed_a, fixed_b=fixed_b, tolerance=1e-13)
        tangent = ret.RetardedTangent(self.graph, gap, root, alpha, fixed_nodes=self.fixed,
            spectral_tolerance=1e-12)
        cut_root = np.sqrt(z*z+abs(gap[self.fixed])**2)
        RE = tangent.apply(z_direction=-1j, fixed_a_direction=1j*root.a[self.fixed]/cut_root,
            fixed_b_direction=1j*root.b[self.fixed]/cut_root).matrix
        rate = rng.normal(size=self.n)*.07+1j*rng.normal(size=self.n)*.04
        rate[self.fixed] = 0.
        Rt = tangent.apply(gap_direction=rate).matrix
        R = ad.Jet.of(tangent.matrix, RE, Rt)
        B = ad.generator(gap, .8, eta=.03, d_time=rate)
        source = rng.normal(size=B.value.shape)*.03+1j*rng.normal(size=B.value.shape)*.02
        B = replace(B, first=source)
        chart = ChartFirstSpectralResponse(self.graph, R, B, alpha, self.fixed,
            kappa_s=self.kappa, retarded_tangent=tangent)
        full = FirstSpectralResponse(self.graph, R, B, alpha, self.fixed, kappa_s=self.kappa)
        known = np.zeros((len(self.fixed), 2, 2), complex)
        result, metrics = chart.apply(R, B, fixed_first=known, tolerance=1e-13)
        expected, previous = full.apply(R, B, fixed_first=known, tolerance=1e-13)
        np.testing.assert_allclose(result.first, expected.first, atol=1e-10, rtol=1e-10)
        self.assertLess(metrics['spectral_first_max'], 1e-11)
        self.assertLess(metrics['normalization_first_all_max'], 1e-12)
        self.assertIs(chart.factor, tangent.factor)
        self.assertTrue(metrics['reuses_retarded_tangent'])
        self.assertLess(chart.estimated_storage_bytes, full.estimated_storage_bytes)

    def test_spatial_temporal_gauge_and_link_velocity_are_retained(self):
        rng = np.random.default_rng(315109)
        R = uniform_jet(self.n, 1.2, 2., .03)
        B = ad.generator(np.full(self.n, 1.2), 2., eta=.03)
        chi, chi_t = rng.normal(size=(2, self.n))*.1
        U = ad.links(-chi, -chi_t)
        target, Bg = ad.transport(U, R, self.kappa), ad.transport(U, B, self.kappa)
        tail, head = self.graph.edges.T
        alpha, alpha_t = chi[head]-chi[tail], chi_t[head]-chi_t[tail]
        raw = replace(target, first=np.zeros_like(target.first))
        chart = ChartFirstSpectralResponse(self.graph, raw, Bg, alpha, self.fixed, kappa_s=self.kappa)
        result, metrics = chart.apply(raw, Bg, alpha_time_per_s=alpha_t,
            fixed_first=target.first[self.fixed], tolerance=1e-13)
        np.testing.assert_allclose(result.first, target.first, atol=1e-12)
        self.assertLess(metrics['spectral_first_max'], 1e-11)
        self.assertLess(metrics['normalization_first_all_max'], 1e-12)

    def test_bad_cut_normalization_is_preserved_and_reported(self):
        R = uniform_jet(self.n, 1.2, 2., .03)
        B = ad.generator(np.full(self.n, 1.2), 2., eta=.03)
        chart = ChartFirstSpectralResponse(self.graph, R, B, fixed_nodes=self.fixed,
            kappa_s=self.kappa)
        wrong = np.broadcast_to(ad.TAU3*.1, (len(self.fixed), 2, 2)).copy()
        result, metrics = chart.apply(R, B, fixed_first=wrong)
        np.testing.assert_array_equal(result.first[self.fixed], wrong)
        self.assertGreater(metrics['normalization_first_all_max'], .1)
        self.assertFalse(metrics['full_equations_within_requested_scaled_tolerance'])
        self.assertFalse(metrics['admission'])


if __name__ == '__main__':
    unittest.main()
