import unittest
import numpy as np
from numpy.polynomial.legendre import leggauss

from pysnspd.experimental.adiabatic_energy_observer import (
    AdiabaticElectronicEnergy, GlobalPower, integrated_balance,
    instantaneous_balance_residual_W, phonon_energy_J)
from pysnspd.experimental.retarded_spatial_usadel import uniform_fields, fields


def quadrature(a, b, n=640):
    x, w = leggauss(n)
    return a+(b-a)*(x+1)/2, (b-a)*w/2


def packet(x):
    """Smooth nonthermal occupation with compact support away from the gap."""
    r = (x-2.)/.7
    inside = abs(r) < 1
    p, derivative = np.zeros_like(x), np.zeros_like(x)
    p[inside] = .1*np.exp(-1/(1-r[inside]**2))
    derivative[inside] = p[inside]*(-2*r[inside]/(.7*(1-r[inside]**2)**2))
    return p, derivative


class AdiabaticEnergyObserverTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.energy, cls.weights = quadrature(1.25, 4.)
        cls.x, cls.xweights = quadrature(1.3, 2.7)

    def setup_bcs(self, gap):
        # Numerical SI example; N0 and volume fix an arbitrary scale, not NbN.
        energy, weights = self.energy, self.weights
        count = np.sqrt(energy**2-gap**2)
        population, p_prime = packet(count)
        observer = AdiabaticElectronicEnergy(energy, weights, [1.], .7, np.ones_like(energy))
        rho = (energy/count)[:, None]
        hL = (1-2*population)[:, None]
        vacuum = .7*gap**2*(np.log(gap/1.76)-.5)
        vacuum_prime = 1.4*gap*np.log(gap/1.76)
        return observer, rho, hL, population, p_prime, count, vacuum, vacuum_prime

    def test_bcs_fixed_count_energy_and_work_match_fixed_energy_observer(self):
        gap = 1.1
        observer, rho, hL, _, p_prime, count, vacuum, vacuum_prime = self.setup_bcs(gap)
        population_x, _ = packet(self.x)
        count_energy = vacuum+2.8*np.sum(self.xweights*np.hypot(self.x, gap)*population_x)
        actual = observer.energy(vacuum, rho, hL)
        self.assertAlmostEqual(actual.total_J, count_energy, delta=2e-10)
        rho_dot = (self.energy*gap/count**3)[:, None]
        h_dot = (2*p_prime*gap/count)[:, None]
        rate = observer.rate(vacuum_prime, rho, hL, rho_dot, h_dot)
        count_rate = vacuum_prime+2.8*np.sum(self.xweights*gap/np.hypot(self.x, gap)*population_x)
        self.assertAlmostEqual(rate.total_rate_W, count_rate, delta=2e-9)
        difference = 2e-5
        side_energies = []
        for side in (gap-difference, gap+difference):
            o, r, h, _, _, _, u, _ = self.setup_bcs(side)
            side_energies.append(o.energy(u, r, h).total_J)
        finite_difference = (side_energies[1]-side_energies[0])/(2*difference)
        self.assertAlmostEqual(rate.total_rate_W, finite_difference, delta=2e-8)
        omitted = observer.rate(vacuum_prime, rho, hL, np.zeros_like(rho_dot), h_dot)
        self.assertGreater(abs(omitted.total_rate_W-count_rate), .02)
        self.assertAlmostEqual(rate.total_rate_W-omitted.total_rate_W, rate.moving_dos_rate_W)

    def test_gauge_phase_changes_bcs_propagators_but_not_observed_energy(self):
        energy, weights = quadrature(.2, 4., n=64)
        observer = AdiabaticElectronicEnergy(energy, weights, [2., 3.], .7, np.tanh(energy/.2))
        hL = np.tanh(energy[:, None]/.2)-.01*np.exp(-((energy[:, None]-2.)/.5)**2)*np.ones((1, 2))
        values = []
        spectral = []
        for phase in (np.array([0., 0.]), np.array([.3, -.7])):
            d = 1.1*np.exp(1j*phase)
            rho = []
            anomalous = []
            for e in energy:
                a, b = uniform_fields(d, .03-1j*e)
                f, _, g = fields(a, b)
                rho.append(g.real); anomalous.append(f)
            values.append(observer.energy(-.5, np.asarray(rho), hL).total_J)
            spectral.append(np.asarray(anomalous))
        self.assertGreater(np.max(abs(spectral[1]-spectral[0])), .1)
        self.assertAlmostEqual(values[0], values[1], delta=2e-14)

    def test_global_observer_detects_missing_external_power_and_has_no_energy_repair(self):
        times = np.array([0., .1, .3, .7])
        exact_energy = 7.+3*times+2*times**2
        power = 3.+4*times
        _, residual = integrated_balance(times, exact_energy, power)
        np.testing.assert_allclose(residual, 0., atol=2e-15)
        work_bad, residual_bad = integrated_balance(times, exact_energy, power-.2)
        np.testing.assert_allclose(residual_bad, .2*times, atol=2e-15)
        np.testing.assert_allclose(exact_energy, 7.+3*times+2*times**2)
        self.assertLess(work_bad[-1], exact_energy[-1]-exact_energy[0])
        external = GlobalPower(10., 2., 1., 3., 2., 1.)
        self.assertEqual(instantaneous_balance_residual_W(1., 1., 1., external), 0.)
        self.assertAlmostEqual(phonon_energy_J([1., 2.], [.2, .3], [2., 3.], [[1., 2.], [3., 4.]], [2., 1.]), 19.6)
        with self.assertRaises(ValueError):
            GlobalPower(1., -1., 0.).net_into_W()


if __name__ == '__main__':
    unittest.main()
