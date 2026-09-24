"""Independent analytic checks of the radial reference; no full campaign."""
import importlib.util
from pathlib import Path
import sys
import unittest
import numpy as np
from scipy.special import iv, ivp, kv, kvp

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location('stage4_radial_reference', ROOT/'sandbox/stage4_core/radial_usadel_reference.py')
radial = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = radial
SPEC.loader.exec_module(radial)


class TestRadialReference(unittest.TestCase):
    def test_uniform_no_winding_has_exact_angle(self):
        r = np.linspace(.001, 4., 81)
        eps, d = 1.7, .8
        result = radial.solve_mode(eps, d, 4., kind='uniform', winding=0,
                                   grid=r, tolerance=1e-7)
        np.testing.assert_allclose(result['theta'], np.arctan(d/eps), atol=2e-12, rtol=0.)
        np.testing.assert_allclose(result['derivative'], 0., atol=2e-12, rtol=0.)

    def test_small_amplitude_vortex_matches_independent_bessel_solution(self):
        # L1 theta-eps theta=-alpha*r has particular solution alpha*r/eps.
        # The two homogeneous solutions are I1(sqrt(eps)r) and K1(...).
        # Impose the SAME finite-rmin Robin and finite-R Dirichlet conditions.
        rmin, outer, eps, alpha = .001, 4., 2., 1e-5
        r = np.r_[np.geomspace(rmin, .2, 60, endpoint=False), np.linspace(.2, outer, 121)]
        result = radial.solve_mode(eps, alpha, outer, rmin=rmin, kind='linear',
                                   grid=r, tolerance=1e-8)
        k = np.sqrt(eps)
        matrix = np.array([[k*ivp(1,k*rmin)-iv(1,k*rmin)/rmin,
                            k*kvp(1,k*rmin)-kv(1,k*rmin)/rmin],
                           [iv(1,k*outer), kv(1,k*outer)]])
        coefficients = np.linalg.solve(matrix, [0., result['outer_angle']-alpha*outer/eps])
        exact = alpha*r/eps+coefficients[0]*iv(1,k*r)+coefficients[1]*kv(1,k*r)
        np.testing.assert_allclose(result['theta'], exact, atol=2e-10, rtol=2e-4)
        self.assertLess(max(abs(np.array(result['boundary_residual']))), 1e-12)

    def test_bracketed_thermal_root_in_strong_depairing_and_zero_current(self):
        for d, eps, gamma in ((.01,.3,100.), (1.4,.3,20.), (1.7,30.,.01), (1.7,.3,0.)):
            theta = radial.uniform_angle(d, eps, gamma)
            residual = eps*np.sin(theta)+gamma*np.sin(theta)*np.cos(theta)-d*np.cos(theta)
            self.assertGreaterEqual(theta, 0.)
            self.assertLess(theta, np.pi/2)
            self.assertAlmostEqual(residual, 0., delta=2e-11)
        self.assertEqual(radial.uniform_angle(0., 1., 10.), 0.)

    def test_uniform_equilibrium_force_and_regular_complex_laplacian(self):
        t = .9/8.65
        gap, residual = radial.thermal_gap(t)
        force, moment, root_residual = radial.uniform_moments([gap], [0.], 600, t)
        self.assertLess(abs(residual), 1e-11)
        self.assertLess(abs(force[0]), 2e-8)
        self.assertEqual(moment[0], 0.)
        r = np.array([1e-6, 1e-5, 1e-4])
        lap = radial.profile(r, gap)[3]
        np.testing.assert_allclose(lap/r, -8*gap/3, rtol=2e-8)


if __name__ == '__main__':
    unittest.main()
