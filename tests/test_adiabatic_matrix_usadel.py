"""Analytic normal, moving-gap and temporal-gauge checks (no detector claim)."""
import unittest
import numpy as np
from pysnspd.experimental import adiabatic_matrix_usadel as ad
from pysnspd.experimental import thermal_spatial_usadel as thermal
from pysnspd.experimental import frozen_kinetic_usadel as frozen


def uniform_jet(nodes, d, energy, eta, d_time=0.):
    z = eta-1j*energy
    root = np.sqrt(z*z+d*d)
    R = (z*ad.TAU3+d*np.array([[0., 1.], [1., 0.]]))/root
    RE = -1j*ad.TAU3/root+1j*z*R/(root*root)
    Rt = d_time*(np.array([[0., 1.], [1., 0.]])/root-d*R/(root*root))
    return ad.Jet.of(np.broadcast_to(R, (nodes, 2, 2)), RE, Rt)


class AdiabaticMatrixTests(unittest.TestCase):
    def setUp(self):
        self.graph = thermal.rectangular_graph(4, 3, 2., 1.5)
        self.n = self.graph.n_nodes
        self.rng = np.random.default_rng(92426)

    def test_normal_storage_and_frozen_diffusion(self):
        h = self.rng.normal(size=(self.n, 2))
        ht = self.rng.normal(size=(self.n, 2))
        R = ad.Jet.of(np.broadcast_to(ad.TAU3, (self.n, 2, 2)))
        H = ad.distribution(h[:, 0], h[:, 1], time_L=ht[:, 0], time_T=ht[:, 1])
        B = ad.generator(np.zeros(self.n), .7)
        residual, flux, _ = ad.kinetic_residual(self.graph, R, H, B, kappa=.4)
        op = frozen.assemble(self.graph, np.zeros(self.n), np.ones(self.n), np.zeros(self.n), np.zeros(self.n))
        np.testing.assert_allclose(ad.projections(residual.value), op.evaluate(h).residual, atol=1e-14)
        np.testing.assert_allclose(ad.projections(residual.first), -.4*self.graph.area_weights[:, None]*ht, atol=1e-14)
        np.testing.assert_allclose(ad.projections(flux.value), op.evaluate(h).edge_flux, atol=1e-14)

    def test_uniform_moving_gap_gives_published_spectral_work(self):
        d, energy, d_time, eta, kappa = 1.2, 2., .07, 0., .4
        R = uniform_jet(self.n, d, energy, eta, d_time)
        B = ad.generator(np.full(self.n, d), energy, d_time=d_time)
        h, hE, ht = .4, .13, -.08
        H = ad.distribution(np.full(self.n, h), np.zeros(self.n), energy_L=hE, time_L=ht)
        residual, _, _ = ad.kinetic_residual(self.graph, R, H, B, kappa=kappa)
        rho, R2 = R.value[:, 0, 0].real, R.value[:, 0, 1].imag
        expected = -self.graph.area_weights*kappa*(rho*ht+R2*d_time*hE)
        np.testing.assert_allclose(ad.projections(residual.first)[:, 0], expected, atol=1e-14)
        spectral, norm = ad.spectral_residual(self.graph, R, B, kappa=kappa)
        self.assertLess(np.max(abs(spectral.first)), 1e-14)
        self.assertLess(np.max(abs(norm.first)), 1e-14)

    def test_temporal_gauge_requires_identity_correction(self):
        d, energy, eta, kappa = 1.2, .8, .03, .4
        R = uniform_jet(self.n, d, energy, eta)
        B = ad.generator(np.full(self.n, d), energy, eta=eta)
        chi, chi_time = .31, .17
        U = ad.links(np.full(self.n, -chi), np.full(self.n, -chi_time))
        rotated = ad.transport(U, R, kappa)
        transformed_B = ad.transport(U, B, kappa)
        np.testing.assert_allclose(transformed_B.first, np.broadcast_to(-1j*kappa*chi_time*ad.I2, B.value.shape), atol=1e-14)
        self.assertGreater(np.max(abs(np.trace(rotated.first, axis1=-2, axis2=-1))), .01)
        residual, norm = ad.spectral_residual(self.graph, rotated, transformed_B, kappa=kappa)
        self.assertLess(np.max(abs(residual.first)), 1e-14)
        self.assertLess(np.max(abs(norm.first)), 1e-14)
        without = ad.Jet.of(rotated.value, rotated.energy, rotated.time)
        _, broken = ad.spectral_residual(self.graph, without, transformed_B, kappa=kappa)
        self.assertGreater(np.max(abs(broken.first)), .01)

    def test_sparse_first_correction_recovers_temporal_gauge(self):
        kappa = .4
        R = uniform_jet(self.n, 1.2, .8, .03)
        B = ad.generator(np.full(self.n, 1.2), .8, eta=.03)
        U = ad.links(np.full(self.n, -.31), np.full(self.n, -.17))
        target, gauged_B = ad.transport(U, R, kappa), ad.transport(U, B, kappa)
        input_R = ad.Jet.of(target.value, target.energy, target.time)
        result, metrics = ad.solve_first_spectral(self.graph, input_R, gauged_B,
            fixed_nodes=self.graph.boundary_nodes,
            fixed_first=target.first[self.graph.boundary_nodes], kappa=kappa, tolerance=1e-13)
        np.testing.assert_allclose(result.first, target.first, atol=2e-11)
        self.assertLess(metrics['normalization_first_max'], 1e-10)
        self.assertLess(metrics['spectral_first_max'], 1e-10)

    def test_spatially_varying_temporal_gauge_transports_link_derivatives(self):
        kappa = .4
        R = uniform_jet(self.n, 1.2, 2., 0.)
        B = ad.generator(np.full(self.n, 1.2), 2.)
        h = ad.distribution(np.full(self.n, .4), np.zeros(self.n), energy_L=.12)
        chi = self.rng.normal(size=self.n)*.1
        chi_time = self.rng.normal(size=self.n)*.08
        U = ad.links(-chi, -chi_time)
        Rg, Bg, hg = (ad.transport(U, field, kappa) for field in (R, B, h))
        tail, head = self.graph.edges.T
        alpha, alpha_time = chi[head]-chi[tail], chi_time[head]-chi_time[tail]
        spectral, norm = ad.spectral_residual(self.graph, Rg, Bg, alpha=alpha,
            alpha_time=alpha_time, kappa=kappa)
        kinetic, _, _ = ad.kinetic_residual(self.graph, Rg, hg, Bg, alpha=alpha,
            alpha_time=alpha_time, kappa=kappa)
        self.assertLess(np.max(abs(spectral.first)), 1e-13)
        self.assertLess(np.max(abs(norm.first)), 1e-13)
        self.assertLess(np.max(abs(kinetic.first)), 1e-13)

    def test_star_associativity_to_the_retained_order(self):
        def random_jet():
            arrays = [self.rng.normal(size=(2, 2))+1j*self.rng.normal(size=(2, 2)) for _ in range(4)]
            return ad.Jet.of(*arrays)
        a, b, c = random_jet(), random_jet(), random_jet()
        left = ad.star(ad.star(a, b, .3), c, .3)
        right = ad.star(a, ad.star(b, c, .3), .3)
        for field in ('value', 'energy', 'time', 'first'):
            np.testing.assert_allclose(getattr(left, field), getattr(right, field), atol=2e-14)


if __name__ == '__main__':
    unittest.main()
