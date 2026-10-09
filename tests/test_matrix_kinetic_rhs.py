"""Analytic storage, diffusion and spectral-work tests; no detector fit."""
import unittest
import numpy as np

from pysnspd.experimental import adiabatic_matrix_usadel as ad
from pysnspd.experimental import thermal_spatial_usadel as thermal
from pysnspd.experimental.matrix_kinetic_rhs import assemble_storage, solve_distribution_rhs


def uniform_jet(nodes, d, energy, eta, d_time=0.):
    z = eta-1j*energy
    root = np.sqrt(z*z+d*d)
    X = np.array([[0., 1.], [1., 0.]])
    R = (z*ad.TAU3+d*X)/root
    RE = -1j*ad.TAU3/root+1j*z*R/(root*root)
    Rt = d_time*(X/root-d*R/(root*root))
    return ad.Jet.of(np.broadcast_to(R, (nodes, 2, 2)), RE, Rt)


class MatrixKineticRHSTest(unittest.TestCase):
    def setUp(self):
        self.graph = thermal.rectangular_graph(4, 3, 2., 1.5)
        self.n = self.graph.n_nodes
        self.rng = np.random.default_rng(100526)

    def solve(self, R, h, hE, B, **kwargs):
        return solve_distribution_rhs(self.graph, R, h, hE, B,
            kappa_s=.4, projected_source=np.zeros_like(h), source_contract='collisionless test', **kwargs)

    def test_sparse_storage_matches_full_matrix_direction_at_nonzero_links(self):
        R = uniform_jet(self.n, 1.2, 2., .03, .07)
        B = ad.generator(np.full(self.n, 1.2), 2., d_time=.07)
        alpha = self.rng.normal(size=len(self.graph.edges))*.1
        M = assemble_storage(self.graph, R, B, alpha=alpha, kappa_s=.4)
        h = self.rng.normal(size=(self.n, 2))*.1
        hE = self.rng.normal(size=(self.n, 2))*.1
        direction = self.rng.normal(size=(self.n, 2))
        def residual(v):
            H = ad.distribution(h[:, 0], h[:, 1], energy_L=hE[:, 0], energy_T=hE[:, 1],
                                time_L=v[:, 0], time_T=v[:, 1])
            r, _, _ = ad.kinetic_residual(self.graph, R, H, B, alpha=alpha, kappa=.4)
            return ad.projections(r.value+r.first)
        np.testing.assert_allclose((M@direction.ravel()).reshape(self.n, 2),
            residual(direction)-residual(np.zeros_like(h)), atol=4e-15)
        self.assertLess(M.nnz, (2*self.n)**2)

    def test_normal_two_mode_rhs_is_mass_weighted_graph_diffusion(self):
        R = ad.Jet.of(np.broadcast_to(ad.TAU3, (self.n, 2, 2)))
        B = ad.generator(np.zeros(self.n), .7)
        h = self.rng.normal(size=(self.n, 2))*.1
        result = self.solve(R, h, np.zeros_like(h), B)
        expected = result.frozen_residual/(.4*self.graph.area_weights[:, None])
        np.testing.assert_allclose(result.distribution_rhs_per_s, expected, atol=1e-14)
        self.assertLess(result.free_projected_max, 2e-15)
        self.assertAlmostEqual(np.dot(self.graph.area_weights, result.distribution_rhs_per_s[:, 0]), 0., delta=2e-15)
        self.assertAlmostEqual(np.dot(self.graph.area_weights, result.distribution_rhs_per_s[:, 1]), 0., delta=2e-15)
        self.assertGreater(np.linalg.norm(result.distribution_rhs_per_s[:, 1]), .1)

    def test_uniform_moving_gap_recovers_level_drift(self):
        d, E, dt = 1.2, 2., .07
        R = uniform_jet(self.n, d, E, 0., dt)
        B = ad.generator(np.full(self.n, d), E, d_time=dt)
        h = np.column_stack((np.full(self.n, .4), np.zeros(self.n)))
        hE = np.column_stack((np.full(self.n, .13), np.zeros(self.n)))
        result = self.solve(R, h, hE, B)
        expected = -d*dt/E*hE[:, 0]
        np.testing.assert_allclose(result.distribution_rhs_per_s[:, 0], expected, atol=3e-15)
        np.testing.assert_allclose(result.distribution_rhs_per_s[:, 1], 0., atol=3e-15)
        self.assertLess(result.free_projected_max, 5e-16)
        self.assertGreater(np.linalg.norm(result.known_temporal_residual), .002)

    def test_prescribed_distribution_reactions_are_retained(self):
        R = ad.Jet.of(np.broadcast_to(ad.TAU3, (self.n, 2, 2)))
        B = ad.generator(np.zeros(self.n), .7)
        h = self.rng.normal(size=(self.n, 2))*.1
        fixed = self.graph.boundary_nodes
        with self.assertRaisesRegex(ValueError, 'Explicit time law'):
            self.solve(R, h, np.zeros_like(h), B, fixed_nodes=fixed)
        result = self.solve(R, h, np.zeros_like(h), B, fixed_nodes=fixed,
                            fixed_rhs_per_s=np.zeros((len(fixed), 2)))
        self.assertLess(result.free_projected_max, 1e-15)
        self.assertGreater(np.linalg.norm(result.boundary_reaction), .1)

    def test_physical_source_sign_and_seconds(self):
        R = ad.Jet.of(np.broadcast_to(ad.TAU3, (self.n, 2, 2)))
        B = ad.generator(np.zeros(self.n), .7)
        h = np.column_stack((np.full(self.n, .4), np.zeros(self.n)))
        physical_direction = np.column_stack((np.full(self.n, -2e9), np.full(self.n, 3e8)))
        kappa = 1e-13
        result = solve_distribution_rhs(self.graph, R, h, np.zeros_like(h), B,
            kappa_s=kappa, projected_source=-kappa*self.graph.area_weights[:, None]*physical_direction,
            source_contract='normal collision residual -kappa*m*hdot')
        np.testing.assert_allclose(result.distribution_rhs_per_s, physical_direction, rtol=1e-14)

    def test_no_pauli_clip_or_eta_bath(self):
        R = uniform_jet(self.n, 1.2, 2., .03)
        B = ad.generator(np.full(self.n, 1.2), 2.)
        h = np.column_stack((np.full(self.n, 1.01), np.zeros(self.n)))
        with self.assertRaisesRegex(ValueError, 'Pauli support'):
            self.solve(R, h, np.zeros_like(h), B)
        h[:, 0] = .4
        with self.assertRaisesRegex(ValueError, 'eta is forbidden'):
            self.solve(R, h, np.zeros_like(h), ad.generator(np.full(self.n, 1.2), 2., eta=.03))

    def test_singular_storage_is_rejected_without_dos_floor(self):
        R = uniform_jet(self.n, 1.2, .7, 0.)
        B = ad.generator(np.full(self.n, 1.2), .7)
        h = np.column_stack((np.full(self.n, .4), np.zeros(self.n)))
        with self.assertRaisesRegex(FloatingPointError, 'Singular|Nonfinite'):
            self.solve(R, h, np.zeros_like(h), B)


if __name__ == '__main__':
    unittest.main()
