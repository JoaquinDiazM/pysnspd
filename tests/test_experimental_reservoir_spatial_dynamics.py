"""Algebraic fixed-radius load tests; no spectral queries or trajectories."""
from types import SimpleNamespace
import unittest
import numpy as np

from pysnspd.experimental.energy_catalog import E_CHARGE_C, HBAR_J_S
from pysnspd.experimental.reservoir_spatial_dynamics import (
    fixed_radius_current_load, fixed_radius_reservoir_kwt,
)


class Mobility:
    scales = SimpleNamespace(t_ref_s=1e-12)

    @staticmethod
    def coefficients(amplitude, temperature):
        return dict(Abar=.4+.1*temperature, tau0bar=.8, R=np.sqrt(1+.7*amplitude**2))


class TestReservoirLoads(unittest.TestCase):
    def setUp(self):
        self.z = np.array([.9*np.exp(.4j), .8*np.exp(.2j), .9*np.exp(-.1j)])
        self.g = np.array([[.15, -.3], [.2, .01], [-.14, .28]])
        self.mass = np.array([.4, .6, 1.2, .7])
        self.mapping = np.array([0, 0, 1, 2])
        self.temperature = np.array([.12, .2, .15, .12])
        self.energy = 2e-20
        self.current = 8e-6

    def response(self, g=None, phi=None):
        return fixed_radius_reservoir_kwt(self.z, self.g if g is None else g,
            self.mass, Mobility(), self.temperature, np.zeros(3) if phi is None else phi,
            quadrature_to_node=self.mapping, terminal_nodes=[0, 2],
            terminal_injection_A=[self.current, -self.current], fixed_radii_bar=[.9, .9],
            energy_scale_J=self.energy)

    def test_phase_sign_and_current_scale_are_prescribed_not_measured(self):
        result = self.response()
        scale = 2*E_CHARGE_C/HBAR_J_S*self.energy
        np.testing.assert_allclose(result.load.phase_load_bar, [-self.current/scale, self.current/scale])
        load = result.load.load_cartesian_bar[:, 0]+1j*result.load.load_cartesian_bar[:, 1]
        np.testing.assert_allclose(np.imag(np.conj(self.z[[0, 2]])*load[[0, 2]]), result.load.phase_load_bar, atol=3e-17)
        changed = self.response(g=3*self.g)
        np.testing.assert_array_equal(result.load.phase_load_bar, changed.load.phase_load_bar)
        np.testing.assert_allclose(changed.load.radial_reaction_bar, 3*result.load.radial_reaction_bar)
        self.assertGreater(np.max(abs(result.kwt.material_velocity_bar[[0, 2]])), .01)

    def test_zero_radius_rate_is_enforced_inside_tensor_solve_and_work_closes(self):
        result = self.response(phi=[1e-4, 4e-5, 0.])
        np.testing.assert_allclose(result.terminal_radius_rates_bar, 0., atol=2e-16)
        radius_velocity = np.real(np.conj(self.z[[0, 2]])/abs(self.z[[0, 2]])*result.kwt.field_velocity_bar[[0, 2]])
        np.testing.assert_allclose(radius_velocity, 0., atol=2e-16)
        self.assertLess(abs(result.radial_reaction_work_rate_bar), 2e-16)
        self.assertLess(abs(result.decomposition_residual_bar), 2e-16)
        self.assertLess(abs(result.kwt.identity_residual_bar), 2e-15)
        self.assertGreaterEqual(np.min(result.kwt.heat_density_bar), 0.)
        np.testing.assert_allclose(result.kwt.node_mass_bar, [1., 1.2, .7])

    def test_independent_polar_tensor_solution_includes_unequal_quadrature_temperatures(self):
        result = self.response()
        load = result.load.load_cartesian_bar[:, 0]+1j*result.load.load_cartesian_bar[:, 1]
        g = self.g[:, 0]+1j*self.g[:, 1]
        expected = np.empty(3, complex)
        for node in range(3):
            matrix = np.zeros((2, 2))
            unit = np.array([self.z[node].real, self.z[node].imag])/abs(self.z[node])
            radial = np.outer(unit, unit)
            for q in np.flatnonzero(self.mapping == node):
                c = Mobility.coefficients(abs(self.z[node]), self.temperature[q])
                matrix += self.mass[q]*2*c["Abar"]*c["tau0bar"]*(c["R"]*radial+(np.eye(2)-radial)/c["R"])
            force = g[node]-load[node]
            velocity = np.linalg.solve(matrix, -np.array([force.real, force.imag]))
            expected[node] = velocity[0]+1j*velocity[1]
        np.testing.assert_allclose(result.kwt.material_velocity_bar, expected, atol=3e-16)

    def test_phase_covariance_does_not_fix_two_dirichlet_phases(self):
        original = self.response()
        angle = np.array([.31, -.12, .4])
        gradient = (self.g[:, 0]+1j*self.g[:, 1])*np.exp(1j*angle)
        transformed = fixed_radius_reservoir_kwt(self.z*np.exp(1j*angle),
            np.column_stack((gradient.real, gradient.imag)), self.mass, Mobility(), self.temperature,
            np.zeros(3), quadrature_to_node=self.mapping, terminal_nodes=[0, 2],
            terminal_injection_A=[self.current, -self.current], fixed_radii_bar=[.9, .9], energy_scale_J=self.energy)
        np.testing.assert_allclose(transformed.kwt.material_velocity_bar,
                                   original.kwt.material_velocity_bar*np.exp(1j*angle), atol=3e-16)
        self.assertAlmostEqual(transformed.reservoir_work_rate_bar, original.reservoir_work_rate_bar, places=14)

    def test_incompatible_radius_or_injections_rejected_without_projection(self):
        args = dict(delta_bar=self.z, gradient_cartesian_bar=self.g, terminal_nodes=[0, 2],
                    terminal_injection_A=[self.current, -self.current], fixed_radii_bar=[.9, .9], energy_scale_J=self.energy)
        for changes in ({"terminal_injection_A": [self.current, -.99*self.current]},
                        {"fixed_radii_bar": [.89, .9]}, {"fixed_radii_bar": [0., .9]},
                        {"terminal_nodes": [0, 0]}, {"terminal_nodes": [0., 2.]},
                        {"energy_scale_J": 0.}):
            with self.subTest(changes=changes), self.assertRaises(ValueError):
                fixed_radius_current_load(**(args|changes))


if __name__ == "__main__":
    unittest.main()
