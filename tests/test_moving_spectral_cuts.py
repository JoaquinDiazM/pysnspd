"""Exact moving superconducting cuts: FD, SI rates, gauge and R1 controls."""
import unittest
import numpy as np
from scipy.constants import Boltzmann, hbar
from scipy.optimize import root as scipy_root

from pysnspd.experimental import thermal_spatial_usadel as thermal
from pysnspd.experimental import retarded_spatial_usadel as retarded
from pysnspd.experimental import adiabatic_matrix_usadel as ad
from pysnspd.experimental import bulk_current_reference as bulk


def bulk_theta(z, d, q, initial=None):
    seed = np.arctan(d/z) if initial is None else initial
    def equations(raw):
        value = raw[0]+1j*raw[1]
        residual = z*np.sin(value)-d*np.cos(value)+q*q*np.sin(value)*np.cos(value)
        return [residual.real, residual.imag]
    answer = scipy_root(equations, [seed.real, seed.imag], tol=1e-12)
    theta = answer.x[0]+1j*answer.x[1]
    if np.max(np.abs(equations(answer.x))) > 1e-11:
        raise RuntimeError('Test bulk root did not converge')
    return theta


class MovingSpectralCutTests(unittest.TestCase):
    def setUp(self):
        self.graph = thermal.rectangular_graph(5, 3, 2., 1.5)
        xy = self.graph.coordinates_bar
        self.fixed = np.flatnonzero((xy[:, 0] == xy[:, 0].min()) |
                                    (xy[:, 0] == xy[:, 0].max()))
        self.d = (.8+.03*np.cos(xy[:, 1]))*np.exp(.12j*xy[:, 0])
        self.z = .22-.7j
        self.alpha = .015*np.sin(np.arange(len(self.graph.edges)))
        aa, bb = retarded.uniform_fields(self.d[self.fixed], self.z)
        self.solution = retarded.solve(self.graph, self.d, self.z, self.alpha,
            fixed_nodes=self.fixed, fixed_a=aa, fixed_b=bb, tolerance=2e-13)
        self.tangent = retarded.RetardedTangent(self.graph, self.d, self.solution,
            self.alpha, fixed_nodes=self.fixed)
        self.rng = np.random.default_rng(51009)

    def test_joint_gap_links_and_independent_moving_cuts_match_finite_difference(self):
        n = self.graph.n_nodes
        dd = .07*(self.rng.normal(size=n)+1j*self.rng.normal(size=n))
        dz = .04-.09j
        dalpha = .03*self.rng.normal(size=len(self.graph.edges))
        da = .025*(self.rng.normal(size=len(self.fixed))+1j*self.rng.normal(size=len(self.fixed)))
        db = .025*(self.rng.normal(size=len(self.fixed))+1j*self.rng.normal(size=len(self.fixed)))
        direction = self.tangent.apply(dd, z_direction=dz, alpha_direction=dalpha,
            fixed_a_direction=da, fixed_b_direction=db)
        h = 2e-6
        results, residuals = [], []
        for sign in (1, -1):
            delta = sign*h
            result = retarded.solve(self.graph, self.d+delta*dd, self.z+delta*dz,
                self.alpha+delta*dalpha, initial_a=self.solution.a, initial_b=self.solution.b,
                fixed_nodes=self.fixed, fixed_a=self.solution.a[self.fixed]+delta*da,
                fixed_b=self.solution.b[self.fixed]+delta*db, tolerance=2e-13)
            results.append(result)
            residuals.append(retarded.evaluate(self.graph, self.d+delta*dd,
                self.z+delta*dz, result.a, result.b, self.alpha+delta*dalpha).residual)
        for field, expected in (('a', direction.da), ('b', direction.db),
                                ('f', direction.df), ('f_tilde', direction.df_tilde),
                                ('g', direction.dg)):
            np.testing.assert_allclose((getattr(results[0], field)-getattr(results[1], field))/(2*h),
                                       expected, rtol=3e-7, atol=3e-9)
        np.testing.assert_allclose((residuals[0]-residuals[1])/(2*h),
                                   direction.residual_derivative, rtol=3e-7, atol=3e-9)
        self.assertLess(direction.free_equation_residual_relative, 2e-14)
        self.assertEqual(direction.fixed_law_residual_max, 0.)
        self.assertGreater(np.max(abs(direction.boundary_reaction)), .01)

    def test_prescribed_full_R_derivative_matches_chart_coordinates(self):
        da = .01*(1+2j)*np.ones(len(self.fixed))
        db = .02*(3-1j)*np.ones(len(self.fixed))
        direct = self.tangent.apply(.03j*self.d, fixed_a_direction=da, fixed_b_direction=db)
        cut_matrix = retarded.field_direction(self.solution.a[self.fixed],
            self.solution.b[self.fixed], da, db)[3]
        matrix = self.tangent.apply(.03j*self.d, fixed_matrix_direction=cut_matrix)
        np.testing.assert_allclose(matrix.matrix, direct.matrix, atol=2e-16)
        np.testing.assert_allclose(matrix.boundary_reaction, direct.boundary_reaction, atol=2e-16)
        self.assertLess(matrix.normalization_derivative_max, 1e-15)

    def test_energy_derivative_of_moving_BCS_anchor_matches_finite_difference(self):
        fixed = self.fixed
        d = self.d[fixed]
        w = np.sqrt(self.z*self.z+abs(d)**2)
        R = retarded.nambu_matrix(self.solution.f[fixed],
            self.solution.f_tilde[fixed], self.solution.g[fixed])
        RE = -1j*ad.TAU3/w[:, None, None]+1j*self.z*R/(w*w)[:, None, None]
        direction = self.tangent.apply(z_direction=-1j, fixed_matrix_direction=RE)
        h = 2e-6
        matrices = []
        for sign in (1, -1):
            z = self.z-sign*1j*h
            aa, bb = retarded.uniform_fields(d, z)
            solved = retarded.solve(self.graph, self.d, z, self.alpha,
                initial_a=self.solution.a, initial_b=self.solution.b,
                fixed_nodes=fixed, fixed_a=aa, fixed_b=bb, tolerance=2e-13)
            matrices.append(retarded.nambu_matrix(solved.f, solved.f_tilde, solved.g))
        np.testing.assert_allclose((matrices[0]-matrices[1])/(2*h),
            direction.matrix, rtol=2e-7, atol=3e-9)

    def test_global_phase_rate_per_second_is_covariant(self):
        omega_per_s = 2.7e11
        rate = self.tangent.apply(1j*omega_per_s*self.d,
            fixed_a_direction=1j*omega_per_s*self.solution.a[self.fixed],
            fixed_b_direction=-1j*omega_per_s*self.solution.b[self.fixed])
        expected = .5j*omega_per_s*(ad.TAU3@self.tangent.matrix-self.tangent.matrix@ad.TAU3)
        np.testing.assert_allclose(rate.matrix, expected, rtol=4e-14, atol=.01)
        self.assertLess(rate.free_equation_residual_relative, 2e-14)
        self.assertLess(rate.normalization_derivative_max/omega_per_s, 3e-14)

    def test_static_local_gauge_transports_tangents_with_cut_reactions(self):
        chi = .2*self.rng.normal(size=self.graph.n_nodes)
        phase = np.exp(1j*chi)
        tail, head = self.graph.edges.T
        gauged_alpha = self.alpha+chi[head]-chi[tail]
        gauged = retarded.solve(self.graph, self.d*phase, self.z, gauged_alpha,
            initial_a=self.solution.a*phase, initial_b=self.solution.b/phase,
            fixed_nodes=self.fixed, fixed_a=self.solution.a[self.fixed]*phase[self.fixed],
            fixed_b=self.solution.b[self.fixed]/phase[self.fixed], tolerance=2e-13)
        other = retarded.RetardedTangent(self.graph, self.d*phase, gauged,
            gauged_alpha, fixed_nodes=self.fixed)
        dd = .03+.07j*self.d
        da = .02j*self.solution.a[self.fixed]
        db = -.02j*self.solution.b[self.fixed]
        base = self.tangent.apply(dd, fixed_a_direction=da, fixed_b_direction=db)
        rotated = other.apply(dd*phase, fixed_a_direction=da*phase[self.fixed],
            fixed_b_direction=db/phase[self.fixed])
        expected = base.matrix.copy()
        expected[:, 0, 1] *= phase
        expected[:, 1, 0] /= phase
        np.testing.assert_allclose(rotated.matrix, expected, atol=2e-15)
        np.testing.assert_allclose(rotated.boundary_reaction[:, 0],
            base.boundary_reaction[:, 0]*phase[self.fixed], atol=2e-15)
        np.testing.assert_allclose(rotated.boundary_reaction[:, 1],
            base.boundary_reaction[:, 1]/phase[self.fixed], atol=2e-15)

    def test_current_carrying_cut_tangent_matches_energy_gap_q_phase_difference(self):
        d, q, phase, z = .8, .31, np.array([-.2, .4]), .09-.6j
        theta = bulk_theta(z, d, q)
        dd, dq, dp, dz = .07, -.025, np.array([.03, -.04]), .02-.11j
        da, db = retarded.current_carrying_cut_tangent(theta, z, d, q, phase,
            gap_direction=dd, q_direction=dq, phase_direction=dp, z_direction=dz)
        h = 2e-6
        cuts = []
        for sign in (1, -1):
            epsilon = sign*h
            th = bulk_theta(z+epsilon*dz, d+epsilon*dd, q+epsilon*dq, theta)
            chart = np.tan(th/2)
            cuts.append((chart*np.exp(1j*(phase+epsilon*dp)),
                         chart*np.exp(-1j*(phase+epsilon*dp))))
        np.testing.assert_allclose((cuts[0][0]-cuts[1][0])/(2*h), da, rtol=1e-8, atol=1e-10)
        np.testing.assert_allclose((cuts[0][1]-cuts[1][1])/(2*h), db, rtol=1e-8, atol=1e-10)
        # The same linear map gives physical s^-1 when rate inputs are s^-1.
        rate_da, rate_db = retarded.current_carrying_cut_tangent(theta, z, d, q, phase,
            gap_direction=dd*1e12, q_direction=dq*1e12,
            phase_direction=dp*1e12, z_direction=dz*1e12)
        np.testing.assert_allclose(rate_da, da*1e12, rtol=2e-15)
        np.testing.assert_allclose(rate_db, db*1e12, rtol=2e-15)

    def test_R1_boundary_follows_temporal_gauge_in_physical_seconds(self):
        graph = self.graph
        n = graph.n_nodes
        energy, gap, eta = .6, .8, .09
        z = eta-1j*energy
        w = np.sqrt(z*z+gap*gap)
        R0 = (z*ad.TAU3+gap*np.array([[0., 1.], [1., 0.]]))/w
        RE = -1j*ad.TAU3/w+1j*z*R0/(w*w)
        R = ad.Jet.of(np.broadcast_to(R0, (n, 2, 2)), RE)
        B = ad.generator(np.full(n, gap), energy, eta=eta)
        kappa_s = hbar/(2*Boltzmann*10.)
        U = ad.links(np.full(n, -.31), np.full(n, -2e11))
        target, BG = ad.transport(U, R, kappa_s), ad.transport(U, B, kappa_s)
        supplied = ad.Jet.of(target.value, target.energy, target.time)
        result, metrics = ad.solve_first_spectral(graph, supplied, BG,
            fixed_nodes=self.fixed, fixed_first=target.first[self.fixed],
            kappa=kappa_s, tolerance=1e-13)
        np.testing.assert_allclose(result.first, target.first, atol=3e-12)
        self.assertLess(metrics['normalization_first_boundary_max'], 1e-13)
        self.assertLess(metrics['normalization_first_max'], 1e-11)
        self.assertTrue(metrics['fixed_first_was_explicit'])
        with self.assertRaisesRegex(ValueError, 'moving-cut star normalization'):
            ad.solve_first_spectral(graph, supplied, BG,
                fixed_nodes=self.fixed, kappa=kappa_s, tolerance=1e-13)

    def test_bulk_R1_recovers_temporal_gauge_with_current(self):
        d, q, z = .8, .31, .09-.6j
        theta = bulk_theta(z, d, q)
        a = np.asarray([np.tan(theta/2)])
        b = a.copy()
        f, ft, g = retarded.fields(a, b)
        R0 = retarded.nambu_matrix(f, ft, g)
        ae, be = retarded.current_carrying_cut_tangent(theta, z, d, q, 0., z_direction=-1j)
        RE = retarded.field_direction(a, b, ae, be)[3]
        original = ad.Jet.of(R0, RE)
        B = ad.generator(np.asarray([d], complex), .6, eta=.09)
        kappa_s = hbar/(2*Boltzmann*10.)
        U = ad.links(np.asarray([-.31]), np.asarray([-2e11]))
        target, BG = ad.transport(U, original, kappa_s), ad.transport(U, B, kappa_s)
        supplied = ad.Jet.of(target.value, target.energy, target.time)
        result, metrics = ad.bulk_first_spectral(supplied, BG, q, kappa_s=kappa_s, tolerance=1e-12)
        np.testing.assert_allclose(result.first, target.first, atol=2e-14)
        self.assertTrue(metrics['consistent'])
        self.assertEqual(metrics['rank_min'], 4)
        self.assertLess(metrics['spectral_first_max'], 1e-13)
        self.assertLess(metrics['normalization_first_max'], 1e-13)

    def test_bulk_R1_moving_gap_q_and_phase_keeps_all_temporal_terms(self):
        d, q, z, phase = .8, .31, .09-.6j, np.asarray([-.2, .4])
        theta = bulk_theta(z, d, q)
        chart = np.tan(theta/2)
        a, b = chart*np.exp(1j*phase), chart*np.exp(-1j*phase)
        f, ft, g = retarded.fields(a, b)
        R0 = retarded.nambu_matrix(f, ft, g)
        ae, be = retarded.current_carrying_cut_tangent(theta, z, d, q, phase, z_direction=-1j)
        RE = retarded.field_direction(a, b, ae, be)[3]
        dd, dq, dp = .04e12, .02e12, np.array([.03e12, -.04e12])
        at, bt = retarded.current_carrying_cut_tangent(theta, z, d, q, phase,
            gap_direction=dd, q_direction=dq, phase_direction=dp)
        Rt = retarded.field_direction(a, b, at, bt)[3]
        B = ad.generator(d*np.exp(1j*phase), .6, eta=.09,
            d_time=np.exp(1j*phase)*(dd+1j*d*dp))
        R = ad.Jet.of(R0, RE, Rt)
        kappa_s = hbar/(2*Boltzmann*10.)
        result, metrics = ad.bulk_first_spectral(R, B, q, q_time_per_s=dq,
            kappa_s=kappa_s, tolerance=1e-12)
        self.assertTrue(metrics['consistent'])
        self.assertLess(metrics['spectral_first_max'], 1e-13)
        self.assertLess(metrics['normalization_first_max'], 1e-13)
        # If Qdot is omitted while Rt contains qdot, no R1 can repair the
        # inconsistent bulk cut. Its residual must stay visible.
        broken, bad = ad.bulk_first_spectral(R, B, q, q_time_per_s=0.,
            kappa_s=kappa_s, require_consistent=False)
        self.assertFalse(bad['consistent'])
        self.assertGreater(bad['spectral_first_max'], 1e-4)
        with self.assertRaisesRegex(FloatingPointError, 'Incompatible moving bulk'):
            ad.bulk_first_spectral(R, B, q, q_time_per_s=0., kappa_s=kappa_s)
        U = ad.links(np.full(2, -.37))
        rotated, rotated_B = ad.transport(U, R, kappa_s), ad.transport(U, B, kappa_s)
        same, other = ad.bulk_first_spectral(rotated, rotated_B, q,
            q_time_per_s=dq, kappa_s=kappa_s)
        np.testing.assert_allclose(same.first,
            ad.transport(U, result, kappa_s).first, atol=1e-14)

    def test_q_zero_bulk_R1_matches_uniform_graph_and_single_matrix_batch(self):
        n = self.graph.n_nodes
        d, energy, eta = .8, .6, .09
        z, w = eta-1j*energy, np.sqrt((eta-1j*energy)**2+d*d)
        R0 = (z*ad.TAU3+d*np.array([[0., 1.], [1., 0.]]))/w
        RE = -1j*ad.TAU3/w+1j*z*R0/(w*w)
        original = ad.Jet.of(R0, RE)
        B0 = ad.generator(np.asarray(d, complex), energy, eta=eta)
        kappa_s = hbar/(2*Boltzmann*10.)
        U0 = ad.links(np.asarray(-.21), np.asarray(-2e11))
        R = ad.transport(U0, original, kappa_s)
        B = ad.transport(U0, B0, kappa_s)
        one, metrics = ad.bulk_first_spectral(ad.Jet.of(R.value, R.energy, R.time),
            B, 0., kappa_s=kappa_s)
        self.assertEqual(one.first.shape, (2, 2))
        repeated_R = ad.Jet.of(np.broadcast_to(R.value, (n, 2, 2)),
            R.energy, R.time)
        repeated_B = ad.Jet.of(np.broadcast_to(B.value, (n, 2, 2)),
            B.energy, B.time, B.first)
        graph, graph_metrics = ad.solve_first_spectral(self.graph, repeated_R,
            repeated_B, fixed_nodes=self.fixed,
            fixed_first=np.broadcast_to(one.first, (len(self.fixed), 2, 2)),
            kappa=kappa_s, tolerance=1e-13)
        np.testing.assert_allclose(graph.first, np.broadcast_to(one.first, graph.first.shape), atol=2e-12)
        self.assertTrue(metrics['consistent'])

    def test_bulk_current_tangent_matches_same_finite_sum_DC_branch(self):
        t, count, q = .12, 256, .27
        reference = bulk.solve_bulk_reference(t, count, q)
        params = dict(width_m=80e-9, Tc_K=10., diffusion_m2_s=5e-5,
                      sheet_resistance_ohm=450.)
        tangent = bulk.bulk_current_tangent(reference, **params)
        h = 2e-6
        lower, upper = [bulk.solve_bulk_reference(t, count, q+sign*h)
                        for sign in (-1, 1)]
        derivative_I = (upper.total_current_A(**params)-lower.total_current_A(**params))/(2*h)
        derivative_gap = (upper.gap_bar-lower.gap_bar)/(2*h)
        self.assertAlmostEqual(tangent.current_A_derivative_per_q_bar/derivative_I, 1., places=8)
        self.assertAlmostEqual(tangent.gap_bar_derivative_per_q_bar/derivative_gap, 1., places=8)
        self.assertAlmostEqual(tangent.q_bar_derivative_per_A*derivative_I, 1., places=8)
        self.assertAlmostEqual(tangent.gap_bar_derivative_per_A*derivative_I/derivative_gap, 1., places=8)
        self.assertEqual(tangent.gap_J_derivative_per_A,
                         Boltzmann*params['Tc_K']*tangent.gap_bar_derivative_per_A)
        self.assertAlmostEqual(tangent.current_A_derivative_per_q_bar,
            bulk.current_density_derivative(reference)*bulk.physical_scales(
                Tc_K=params['Tc_K'],diffusion_m2_s=params['diffusion_m2_s'],
                sheet_resistance_ohm=params['sheet_resistance_ohm'])['current_unit_A']*
            params['width_m']/bulk.physical_scales(Tc_K=params['Tc_K'],
                diffusion_m2_s=params['diffusion_m2_s'],
                sheet_resistance_ohm=params['sheet_resistance_ohm'])['ell0_m'], places=16)

    def test_fixed_cut_default_and_malformed_boundary_inputs(self):
        self.assertFalse(self.tangent.jacobian.data.flags.writeable)
        with self.assertRaises(ValueError):
            self.tangent.jacobian.data[0] = 0.
        old = self.tangent.apply(.03j*self.d)
        np.testing.assert_array_equal(old.da[self.fixed], np.zeros(len(self.fixed)))
        np.testing.assert_array_equal(old.db[self.fixed], np.zeros(len(self.fixed)))
        with self.assertRaises(ValueError):
            self.tangent.apply(fixed_a_direction=np.zeros(len(self.fixed)))
        with self.assertRaises(ValueError):
            self.tangent.apply(fixed_matrix_direction=np.broadcast_to(ad.I2, (len(self.fixed), 2, 2)))
        with self.assertRaises(ValueError):
            self.tangent.apply(alpha_direction=1j)
        with self.assertRaises(ValueError):
            retarded.current_carrying_cut_tangent(.1, self.z, .8, .2, 0.)
        with self.assertRaises(ValueError):
            retarded.current_carrying_cut_tangent(.1, self.z, .8, .2, 0.,
                phase_direction=1j)


if __name__ == '__main__':
    unittest.main()
