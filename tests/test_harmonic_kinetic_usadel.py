"""Exact normal, gauge and graph-mode controls of finite-frequency kinetics."""
import unittest
import numpy as np
from scipy.linalg import eigh
from pysnspd.experimental import harmonic_kinetic_usadel as hk
from pysnspd.experimental import thermal_spatial_usadel as thermal
from pysnspd.experimental import frozen_kinetic_usadel as frozen
from pysnspd.experimental.harmonic_usadel import HarmonicUsadelFactor


def uniform(d, energy, eta):
    z = eta-1j*energy
    R = np.array([[z, d], [np.conj(d), -z]], complex)/np.sqrt(z*z+abs(d)**2)
    A = -hk.TAU3@R.conj().T@hk.TAU3
    B = np.array([[-1j*energy, d], [np.conj(d), 1j*energy]], complex)
    return R, A, B


class HarmonicKineticTests(unittest.TestCase):
    def setUp(self):
        self.graph = thermal.rectangular_graph(4, 3, 2., 1.5)
        self.n = self.graph.n_nodes
        self.rng = np.random.default_rng(4926)
        self.zero = np.zeros((2, 2), complex)

    def test_normal_harmonic_diffusion_and_storage(self):
        omega = .7
        Rp, Ap, Bp = uniform(0., 1+omega/2, .03)
        Rm, Am, Bm = uniform(0., 1-omega/2, .03)
        op = hk.assemble(self.graph, Rp, Rm, Ap, Am, Bp, Bm,
            self.zero, self.zero, self.zero, .8, .7)
        h = self.rng.normal(size=(self.n, 2))+1j*self.rng.normal(size=(self.n, 2))
        obs = op.evaluate(h)
        tail, head = self.graph.edges.T
        flux = self.graph.conductance[:, None]*(h[head]-h[tail])
        expected = .5j*omega*self.graph.area_weights[:, None]*h
        np.add.at(expected, tail, flux); np.add.at(expected, head, -flux)
        np.testing.assert_allclose(obs.residual, expected, atol=1e-14)
        np.testing.assert_allclose(obs.edge_flux, flux, atol=1e-14)
        np.testing.assert_allclose(op.matrix@h.ravel()+op.source.ravel(), obs.residual.ravel(), atol=1e-14)

    def test_uniform_temporal_pure_gauge_and_contact_solve(self):
        omega, theta, d, energy, temperature = .3, .4+.2j, 1.2, .8, .25
        Rp, Ap, Bp = uniform(d, energy+omega/2, .03)
        Rm, Am, Bm = uniform(d, energy-omega/2, .03)
        D = np.array([[0., d], [d, 0.]], complex)
        voltage = .5j*omega*theta
        dB = .5j*theta*(hk.TAU3@D-D@hk.TAU3)+1j*voltage*hk.I2
        dR = .5j*theta*(hk.TAU3@Rm-Rp@hk.TAU3)
        dA = .5j*theta*(hk.TAU3@Am-Ap@hk.TAU3)
        hp, hm = np.tanh((energy+omega/2)/(2*temperature)), np.tanh((energy-omega/2)/(2*temperature))
        expected = np.tile([0., -.5j*theta*(hp-hm)], (self.n, 1))
        op = hk.assemble(self.graph, Rp, Rm, Ap, Am, Bp, Bm, dR, dA, dB, hp, hm)
        observed = op.evaluate(expected)
        self.assertLess(np.max(abs(observed.residual)), 1e-14)
        self.assertLess(np.max(abs(observed.edge_flux)), 1e-14)
        solved, _, _ = op.solve(self.graph.boundary_nodes, expected)
        np.testing.assert_allclose(solved, expected, atol=2e-13)
        # Exercise the actual spectral-to-kinetic connection with contact data.
        retarded = HarmonicUsadelFactor(self.graph, Rp, Rm,
            Bp+.03*hk.TAU3, Bm+.03*hk.TAU3,
            fixed_nodes=self.graph.boundary_nodes).solve(dB, fixed_response=dR).delta_R
        advanced = HarmonicUsadelFactor(self.graph, Ap, Am,
            Bp-.03*hk.TAU3, Bm-.03*hk.TAU3,
            fixed_nodes=self.graph.boundary_nodes).solve(dB, fixed_response=dA).delta_R
        joined = hk.assemble(self.graph, Rp, Rm, Ap, Am, Bp, Bm,
            retarded, advanced, dB, hp, hm)
        joined_h, joined_result, _ = joined.solve(self.graph.boundary_nodes, expected)
        np.testing.assert_allclose(joined_h, expected, atol=2e-12)
        self.assertLess(np.max(abs(joined_result.residual)), 2e-12)

    def test_zero_frequency_matches_frozen_complex_extension(self):
        R, A, B = uniform(1.2, .8, .03)
        op = hk.assemble(self.graph, R, R, A, A, B, B, self.zero, self.zero, self.zero, .8, .8)
        static = frozen.assemble(self.graph, np.full(self.n, 1.2),
            np.full(self.n, R[0, 0]), np.full(self.n, R[0, 1]), np.full(self.n, R[1, 0]))
        np.testing.assert_allclose(op.matrix.toarray(), static.matrix.toarray(), atol=1e-14)

    def test_modal_restriction_equals_full_graph(self):
        omega = .4
        Rp, Ap, Bp = uniform(1.2, .8+omega/2, .03)
        Rm, Am, Bm = uniform(1.2, .8-omega/2, .03)
        dR = self.rng.normal(size=(2, 2))+1j*self.rng.normal(size=(2, 2))
        dA = self.rng.normal(size=(2, 2))+1j*self.rng.normal(size=(2, 2))
        dB = hk.gap_vertex(.01+.02j, -.03+.01j, .007j)
        tail, head = self.graph.edges.T
        L = np.zeros((self.n, self.n))
        for i, j, c in zip(tail, head, self.graph.conductance):
            L[i, i] += c; L[j, j] += c; L[i, j] -= c; L[j, i] -= c
        values, vectors = eigh(L, np.diag(self.graph.area_weights))
        lam, mode = values[3], vectors[:, 3]
        hp, hm = .8, .7
        full = hk.assemble(self.graph, Rp, Rm, Ap, Am, Bp, Bm,
            mode[:, None, None]*dR, mode[:, None, None]*dA,
            mode[:, None, None]*dB, hp, hm)
        reduced = hk.modal_coefficients(lam, Rp, Rm, Ap, Am, Bp, Bm, dR, dA, dB, hp, hm)
        h = np.array([.06+.02j, -.04+.03j])
        observed = full.evaluate(mode[:, None]*h)
        expected = self.graph.area_weights[:, None]*mode[:, None]*(reduced['matrix']@h+reduced['source'])
        np.testing.assert_allclose(observed.residual, expected, atol=3e-14)

    def test_gap_cartesian_keeps_both_complex_fourier_amplitudes(self):
        K = self.rng.normal(size=(self.n, 2, 2))+1j*self.rng.normal(size=(self.n, 2, 2))
        result = hk.gap_cartesian_moments(K, self.graph.area_weights)
        np.testing.assert_allclose(result[:, 0]+1j*result[:, 1], 1j*self.graph.area_weights*K[:, 0, 1])
        np.testing.assert_allclose(result[:, 0]-1j*result[:, 1], 1j*self.graph.area_weights*K[:, 1, 0])
        self.assertGreater(np.max(abs(result.imag)), .01)

    def test_dynamic_charge_ward_at_each_energy(self):
        d, energy, omega, hp, hm = 1.2, .8, .4, .8, .7
        Rp, Ap, Bp = uniform(d, energy+omega/2, .03)
        Rm, Am, Bm = uniform(d, energy-omega/2, .03)
        x, y, voltage = (self.rng.normal(size=self.n)+1j*self.rng.normal(size=self.n) for _ in range(3))
        dR, dA = (self.rng.normal(size=(self.n, 2, 2))+1j*self.rng.normal(size=(self.n, 2, 2)) for _ in range(2))
        dB = hk.gap_vertex(x, y, voltage)
        op = hk.assemble(self.graph, Rp, Rm, Ap, Am, Bp, Bm, dR, dA, dB, hp, hm)
        h = self.rng.normal(size=(self.n, 2))+1j*self.rng.normal(size=(self.n, 2))
        result = op.evaluate(h)
        tail, head = self.graph.edges.T
        divergence = np.zeros(self.n, complex)
        np.add.at(divergence, tail, result.edge_flux[:, 1])
        np.add.at(divergence, head, -result.edge_flux[:, 1])
        Kp, Km = hp*(Rp-Ap), hm*(Rm-Am)
        average = np.broadcast_to((Kp+Km)/2, (self.n, 2, 2))
        base_moments = hk.gap_cartesian_moments(average, self.graph.area_weights)
        torque = (d*result.gap_cartesian[:, 1]+x*base_moments[:, 1]-y*base_moments[:, 0])/2
        storage = 1j*omega*self.graph.area_weights*np.trace(result.delta_K, axis1=-2, axis2=-1)/8
        voltage_term = -1j*self.graph.area_weights*voltage*np.trace(hk.TAU3@(Km-Kp))/8
        np.testing.assert_allclose(result.residual[:, 1], divergence+torque+storage+voltage_term, atol=3e-14)


if __name__ == '__main__':
    unittest.main()
