"""Physical sign/factor controls for the vectorized final coupling experiment."""
import unittest
import numpy as np
from sandbox.stage4_core import coupled_response as campaign
from pysnspd.experimental import harmonic_kinetic_usadel as kinetic


class CoupledResponseTests(unittest.TestCase):
    def setUp(self):
        self.rng = np.random.default_rng(42924)

    def test_vectorized_kinetic_response_matches_validated_modal_operator(self):
        d, temperature, omega, eta, lam = 1.2, .3, .27, .035, .16
        energy = np.array([-2.4, -1.1, -.21, .17, 1.04, 2.3])
        Rp, Ap, bp, Bp = campaign.backgrounds(d, energy+omega/2, eta)
        Rm, Am, bm, Bm = campaign.backgrounds(d, energy-omega/2, eta)
        hp, hm = np.tanh((energy+omega/2)/(2*temperature)), np.tanh((energy-omega/2)/(2*temperature))
        dR, dA = campaign.spectral_response(Rp, Rm, Ap, Am, bp, bm, lam)
        K, h, error, coefficient = campaign.kinetic_modal(Rp, Rm, Ap, Am, Bp, Bm, dR, dA, hp, hm, lam)
        self.assertLess(error, 1e-13)
        for index in range(len(energy)):
            for forcing in range(3):
                reduced = kinetic.modal_coefficients(lam, Rp[index], Rm[index], Ap[index], Am[index],
                    Bp[index], Bm[index], dR[index, forcing], dA[index, forcing],
                    campaign.FORCING[forcing], hp[index], hm[index])
                expected_h = np.linalg.solve(reduced['matrix'], -reduced['source'])
                expected_K = reduced['K_source']+np.einsum('c,cij->ij', expected_h, reduced['K_basis'])
                np.testing.assert_allclose(coefficient[index], reduced['matrix'], atol=1e-13)
                np.testing.assert_allclose(h[index, :, forcing], expected_h, atol=1e-12)
                np.testing.assert_allclose(K[index, forcing], expected_K, atol=1e-12)

    def test_normal_film_port_has_positive_normal_conductance(self):
        for omega in (.07, .7, 3.):
            task = dict(d=0., t=.2, omega=omega, eta=.02, order=32, cutoff=12.,
                matsubara_count=256, graph_conductance=.47, sheet_resistance_ohm=608.)
            result = campaign.film_port_coefficients(task)
            self.assertAlmostEqual(result['admittance'].real, .47/608., places=12)
            self.assertLess(abs(result['admittance'].imag), 1e-14)
            self.assertLess(abs(result['neutral']), 2e-12)
            self.assertAlmostEqual(result['current'].real, -2., places=11)
            self.assertEqual(result['stiffness'], 0.)

    def test_normal_modes_dynamic_ward_closes_with_neutrality(self):
        omega, lam, temperature = .3, .2, .25
        energy, weights = campaign.energy_rule(0., omega, .02, 32, 12.)
        Rp, Ap, bp, Bp = campaign.backgrounds(0., energy+omega/2, .02)
        Rm, Am, bm, Bm = campaign.backgrounds(0., energy-omega/2, .02)
        hp, hm = np.tanh((energy+omega/2)/(2*temperature)), np.tanh((energy-omega/2)/(2*temperature))
        dR, dA = campaign.spectral_response(Rp, Rm, Ap, Am, bp, bm, lam)
        K, h, _, _ = campaign.kinetic_modal(Rp, Rm, Ap, Am, Bp, Bm, dR, dA, hp, hm, lam)
        # Voltage forcing only: gap terms vanish exactly in the normal state.
        charge = np.dot(weights/2, campaign.projection(K[:, 2])[:, 0])
        divergence = 0j
        integrated_residual = 0j
        for index in range(len(energy)):
            mode = kinetic.modal_coefficients(lam, Rp[index], Rm[index], Ap[index], Am[index],
                Bp[index], Bm[index], dR[index, 2], dA[index, 2], campaign.FORCING[2], hp[index], hm[index])
            divergence += weights[index]*(mode['divergence_matrix']@h[index, :, 2]+mode['divergence_source'])[1]
            integrated_residual += weights[index]*(mode['matrix']@h[index, :, 2]+mode['source'])[1]
        np.testing.assert_allclose(integrated_residual, divergence+1j*omega*(charge+1), atol=2e-12)
        self.assertLess(abs(integrated_residual), 2e-12)

    def test_uniform_temporal_gauge_cancels_voltage_charge(self):
        d, t, omega, eta = 1.2, .25, .4, .015
        energy, weights = campaign.energy_rule(d, omega, eta, 32, 30.)
        Rp, Ap, bp, Bp = campaign.backgrounds(d, energy+omega/2, eta)
        Rm, Am, bm, Bm = campaign.backgrounds(d, energy-omega/2, eta)
        hp, hm = np.tanh((energy+omega/2)/(2*t)), np.tanh((energy-omega/2)/(2*t))
        dR, dA = campaign.spectral_response(Rp, Rm, Ap, Am, bp, bm, 0.)
        # Unit voltage has a Josephson phase delta_y=-2i*d/omega.
        forcing = np.array([0., -2j*d/omega, 1.])
        expected_R = (campaign.Z@Rm-Rp@campaign.Z)/omega
        expected_A = (campaign.Z@Am-Ap@campaign.Z)/omega
        np.testing.assert_allclose(np.einsum('efij,f->eij', dR, forcing), expected_R, atol=2e-12)
        np.testing.assert_allclose(np.einsum('efij,f->eij', dA, forcing), expected_A, atol=2e-12)
        Kp, Km = hp[:, None, None]*(Rp-Ap), hm[:, None, None]*(Rm-Am)
        expected_K = (campaign.Z@Km-Kp@campaign.Z)/omega
        # Finite-cutoff tail remains explicit; it falls as Delta^2/Emax^2.
        neutral = 1+np.dot(weights/2, campaign.projection(expected_K)[:, 0])
        self.assertLess(abs(neutral), .002)

    def test_static_matsubara_stiffness_matches_known_uniform_limits(self):
        d, temperature = 1.76392592757718, .9/8.65
        ep = 2*np.pi*temperature*(np.arange(256)+.5)
        hx, hy = campaign.thermal_stiffness(d, temperature, ep, 0.)
        self.assertGreater(hx, 0.)
        self.assertEqual(hy, 0.)
        g = ep/np.hypot(ep, d)
        h = 1e-6
        actual_slope = campaign.thermal_stiffness(d, temperature, ep, h)[1]/h
        expected_slope = 4*np.pi*temperature*np.sum(g*g/(ep*ep))
        self.assertAlmostEqual(actual_slope/expected_slope, 1., places=5)


if __name__ == '__main__':
    unittest.main()
