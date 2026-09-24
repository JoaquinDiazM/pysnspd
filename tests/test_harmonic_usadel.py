"""Analytic finite-frequency Nambu controls; no trajectory or long sweep."""
import unittest
import numpy as np
from scipy.linalg import eigh
from scipy.sparse import coo_matrix, diags
from scipy.sparse.linalg import spsolve

from pysnspd.experimental.thermal_spatial_usadel import rectangular_graph
from pysnspd.experimental.harmonic_usadel import (
    HarmonicUsadelFactor, I2, TAU3, response_residual,
)


def generator(gap, z):
    return np.array([[z, gap], [np.conj(gap), -z]], complex)


def uniform(gap, z):
    B = generator(gap, z)
    root = np.sqrt(z*z+abs(gap)**2)
    if root.real < 0:
        root = -root
    return B/root, B, root


def rotate(matrices, phases):
    result = np.broadcast_to(matrices, (len(phases), 2, 2)).copy()
    result[:, 0, 1] *= np.exp(1j*phases)
    result[:, 1, 0] *= np.exp(-1j*phases)
    return result


class HarmonicUsadelTests(unittest.TestCase):
    def setUp(self):
        self.graph = rectangular_graph(5, 4, 5., 3.)
        self.fixed = self.graph.boundary_nodes
        self.free = np.ones(self.graph.n_nodes, bool)
        self.free[self.fixed] = False

    def test_zero_frequency_matches_full_complex_uniform_gap_derivative_and_reuses_LU(self):
        R, B, s = uniform(1.3+.2j, .3-.7j)
        delta = np.array([[0., .2+.15j], [.1-.05j, 0.]], complex)
        exact = delta/s-B*np.trace(B@delta)/(2*s**3)
        factor = HarmonicUsadelFactor(self.graph, R, R, B, B, fixed_nodes=self.fixed)
        saved_factor = factor.factor
        result = factor.solve(delta, exact)
        np.testing.assert_allclose(result.delta_R, np.broadcast_to(exact, result.delta_R.shape), atol=2e-13)
        np.testing.assert_allclose(factor.solve(2*delta, 2*exact).delta_R, 2*result.delta_R, atol=2e-13)
        self.assertIs(saved_factor, factor.factor)
        self.assertLess(result.metrics['full_spectral_scaled_max'], 1e-12)
        self.assertLess(result.metrics['normalization_absolute_max'], 1e-12)
        self.assertGreater(np.max(abs(result.delta_R.imag)), .01)
        with self.assertRaises(ValueError):
            factor.solve(delta)

    def test_normal_gap_response_at_finite_frequency_matches_complex_diffusion_resolvent(self):
        graph = self.graph
        center, omega, eta = .4, .7, .3
        z = eta-1j*center
        Bp = (eta-1j*(center+omega/2))*TAU3
        Bm = (eta-1j*(center-omega/2))*TAU3
        source = np.zeros((graph.n_nodes, 2, 2), complex)
        x, y = graph.coordinates_bar.T
        shape = np.cos(np.pi*x/5)*np.cos(np.pi*y/3)
        source[:, 0, 1] = (.3+.2j)*shape
        source[:, 1, 0] = (-.1+.25j)*shape
        factor = HarmonicUsadelFactor(graph, TAU3, TAU3, Bp, Bm, fixed_nodes=self.fixed)
        result = factor.solve(source, np.zeros((len(self.fixed), 2, 2), complex))
        tail, head = graph.edges.T; c = graph.conductance
        L = coo_matrix((np.r_[c, c, -c, -c],
            (np.r_[tail, head, tail, head], np.r_[tail, head, head, tail])),
            shape=(graph.n_nodes, graph.n_nodes)).tocsr()
        resolvent = (L+z*diags(graph.area_weights))[self.free][:, self.free]
        expected = np.zeros_like(result.delta_R)
        for a, b in ((0, 1), (1, 0)):
            expected[self.free, a, b] = spsolve(resolvent,
                graph.area_weights[self.free]*source[self.free, a, b])
        np.testing.assert_allclose(result.delta_R, expected, atol=2e-13)
        self.assertGreater(np.max(abs(result.delta_R.imag)), .01)
        self.assertLess(np.max(abs(result.normalization_residual)), 1e-12)
        self.assertLess(np.max(abs(result.spectral_residual[self.free])), 1e-12)

    def test_finite_frequency_pure_gauge_retains_identity_and_solves_advanced_independently(self):
        graph = self.graph
        center, omega, eta, gap = .6, .8, .2, 1.1+.3j
        theta = .2+.03j
        Rp0, Bp0, _ = uniform(gap, eta-1j*(center+omega/2))
        Rm0, Bm0, _ = uniform(gap, eta-1j*(center-omega/2))
        phases = .07*graph.coordinates_bar[:, 0]-.04*graph.coordinates_bar[:, 1]
        alpha = phases[graph.edges[:, 1]]-phases[graph.edges[:, 0]]
        Rp, Rm, Bp, Bm = (rotate(value, phases) for value in (Rp0, Rm0, Bp0, Bm0))
        D = rotate(generator(gap, 0), phases)
        delta_D = .5j*theta*(TAU3@D-D@TAU3)
        delta_v = .5j*omega*theta
        delta_B = delta_D+1j*delta_v*I2
        exact = .5j*theta*(TAU3@Rm-Rp@TAU3)
        factor = HarmonicUsadelFactor(graph, Rp, Rm, Bp, Bm, alpha=alpha, fixed_nodes=self.fixed)
        retarded = factor.solve(delta_B, exact[self.fixed])
        np.testing.assert_allclose(retarded.delta_R, exact, atol=5e-13)
        self.assertGreater(retarded.metrics['maximum_identity_response'], .01)
        Ap, Am = (-TAU3@np.swapaxes(np.conj(value), -1, -2)@TAU3 for value in (Rp, Rm))
        Bap = rotate(generator(gap, -eta-1j*(center+omega/2)), phases)
        Bam = rotate(generator(gap, -eta-1j*(center-omega/2)), phases)
        exact_advanced = .5j*theta*(TAU3@Am-Ap@TAU3)
        advanced_factor = HarmonicUsadelFactor(graph, Ap, Am, Bap, Bam, alpha=alpha, fixed_nodes=self.fixed)
        advanced = advanced_factor.solve(delta_B, exact_advanced[self.fixed])
        np.testing.assert_allclose(advanced.delta_R, exact_advanced, atol=5e-13)
        naive_same_frequency = -TAU3@np.swapaxes(np.conj(retarded.delta_R), -1, -2)@TAU3
        self.assertGreater(np.max(abs(advanced.delta_R-naive_same_frequency)), .01)
        for result in (retarded, advanced):
            self.assertLess(np.max(abs(result.spectral_residual[self.free])), 1e-12)
            self.assertLess(np.max(abs(result.normalization_residual)), 1e-12)
            self.assertTrue(result.metrics['full_four_entry_residual_checked'])

    def test_uniform_modal_closed_response_matches_full_spatial_factor_for_three_forcings(self):
        from sandbox.stage4_core.harmonic_uniform_reference import shifted_branches, modal_response
        graph = self.graph
        tail, head = graph.edges.T; c = graph.conductance
        L = coo_matrix((np.r_[c, c, -c, -c],
            (np.r_[tail, head, tail, head], np.r_[tail, head, head, tail])),
            shape=(graph.n_nodes, graph.n_nodes)).tocsr()
        eigenvalues, eigenvectors = eigh(L[self.free][:, self.free].toarray(),
                                        np.diag(graph.area_weights[self.free]))
        lam = eigenvalues[1]
        mode = np.zeros(graph.n_nodes)
        mode[self.free] = eigenvectors[:, 1]
        plus, minus = shifted_branches(.7, .6, 1.1+.3j, .25)
        Rp, Bp = plus.matrix, plus.generator
        Rm, Bm = minus.matrix, minus.generator
        factor = HarmonicUsadelFactor(graph, Rp, Rm, Bp, Bm, fixed_nodes=self.fixed)
        for forcing in (np.array([[0, 1], [1, 0]], complex),
                        np.array([[0, 1j], [-1j, 0]], complex), 1j*I2):
            # Independent exact modal solution follows by multiplying the
            # mixed-normalized matrix equation by R_plus; this is not used
            # by the general sparse spatial factor under test.
            modal = modal_response(plus, minus, lam, forcing)
            amplitude = modal.matrix
            self.assertLess(modal.equation_residual_max, 1e-12)
            self.assertLess(modal.normalization_residual_max, 1e-12)
            exact = mode[:, None, None]*amplitude
            solved = factor.solve(mode[:, None, None]*forcing,
                                  np.zeros((len(self.fixed), 2, 2), complex))
            np.testing.assert_allclose(solved.delta_R, exact, atol=5e-13)
            self.assertLess(solved.metrics['full_spectral_scaled_max'], 1e-12)
            self.assertLess(solved.metrics['normalization_absolute_max'], 1e-12)


if __name__ == '__main__':
    unittest.main()
