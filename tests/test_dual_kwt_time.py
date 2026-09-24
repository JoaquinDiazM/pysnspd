"""Focused checks of the dual KWT campaign's physical normalization and CFL."""
import unittest
import numpy as np
from scipy.linalg import eigvalsh

from sandbox.stage4_core.dual_kwt_time import (
    uniform_reference, explicit_step_bound, compare_observations,
    response, material_options)
from pysnspd.experimental import thermal_spatial_usadel as thermal
from pysnspd.experimental.thermal_weak_response import ThermalKWTNormal, SpectralTangent
from pysnspd.experimental.thermal_stable_newton import solve_frequency


class DualKWTTimeTests(unittest.TestCase):
    def setUp(self):
        self.plan = dict(T_K=.9, Tc_K=8.65, tau_ee_Tc_ps=.5,
            tau_ep_Tc_ps=2.47, euler_safety=.5,
            refinement_absolute_norm=1e-8, refinement_relative_limit=.02)
        self.graph = thermal.rectangular_graph(5, 4, 5., 3.)

    def test_uniform_equilibrium_is_stationary_on_the_actual_graph(self):
        graph = self.graph
        d, G, frequencies = uniform_reference(graph, .9/8.65, 32)
        self.assertLess(np.max(abs(G)), 1e-12)
        for e in frequencies:
            value = thermal.spectral_energy_gradient(graph, d, e, d/e)
            self.assertLess(np.max(abs(value.residual)), 1e-12)
        model = ThermalKWTNormal(graph, d, **material_options(self.plan))
        result = model.response(G, np.zeros(len(graph.edges)))
        self.assertLess(np.max(abs(result['velocity'])), 1e-12)

    def test_uniform_phase_spectrum_and_normal_potential_fit_analytic_bound(self):
        graph = self.graph
        d, G, frequencies = uniform_reference(graph, .9/8.65, 8)
        model = ThermalKWTNormal(graph, d, **material_options(self.plan))
        free = np.flatnonzero(model.free)
        tangents = []
        for e in frequencies:
            solution = solve_frequency(graph, d, e, fixed_nodes=graph.boundary_nodes,
                fixed_u=d[graph.boundary_nodes]/e, initial_u=d/e, tol=1e-10)
            tangents.append(SpectralTangent(graph, d, solution))
        t = self.plan['T_K']/self.plan['Tc_K']
        jacobian = np.zeros((len(free), len(free)))
        for j, node in enumerate(free):
            direction = np.zeros(graph.n_nodes, complex)
            direction[node] = 1j
            force = 2*graph.area_weights*direction*np.log(t)
            current = np.zeros(len(graph.edges))
            for e, tangent in zip(frequencies, tangents):
                item = tangent.apply(direction)
                force += 4*np.pi*t*graph.area_weights*(direction/e-item.df)
                current += 2*np.pi*t*item.current_derivative
            response_j = model.rhs_tangent(direction, force, G,
                np.zeros(len(graph.edges)), current)['velocity_direction']
            jacobian[:, j] = -response_j[free].imag/model.tD_ps
        rootmass = np.sqrt(graph.area_weights[free])
        symmetric = rootmass[:, None]*jacobian/rootmass[None, :]
        self.assertLess(np.max(abs(symmetric-symmetric.T)), 1e-10)
        eigenvalues = eigvalsh(symmetric)
        bound = explicit_step_bound(graph, d, frequencies, self.plan)
        self.assertGreater(eigenvalues[0], 0)
        self.assertLessEqual(eigenvalues[-1], bound['maximum_rate_bound_per_ps']*(1+1e-10))
        self.assertLess(bound['primary_step_ps']*eigenvalues[-1], 1.)

    def test_zero_cross_signal_uses_shared_initial_scale_and_power_uses_ps(self):
        graph = self.graph
        d, G, _ = uniform_reference(graph, .9/8.65, 8)
        force = np.zeros(graph.n_nodes, complex)
        force[6] = .01+.02j
        current = np.zeros(len(graph.edges))
        model, result = response(graph, d, [force, current, .5], G, self.plan)
        expected = (result['kwt_loss']+result['normal_loss'])/model.tD_ps
        self.assertEqual(result['loss_per_ps'], expected)
        fields = {}
        for method in ('primary', 'refined'):
            for probe in ('amplitude', 'angular_phase'):
                fields[method+'_'+probe] = dict(gap_difference=np.ones(graph.n_nodes)*.001,
                    force_density=np.ones(graph.n_nodes)*.001,
                    phase_torque_density=np.ones(graph.n_nodes)*(.001 if probe == 'angular_phase' else 0),
                    current=np.ones(len(graph.edges))*.001)
        later = {name: {key: value.copy() for key, value in row.items()} for name, row in fields.items()}
        later['primary_amplitude']['phase_torque_density'] += 1e-6
        result = compare_observations(graph, [(0., fields), (1., later)], self.plan)
        self.assertTrue(result['all_admitted'])


if __name__ == '__main__':
    unittest.main()
