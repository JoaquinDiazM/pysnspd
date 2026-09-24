"""Physical checks for the DC preparation ports, phase convention and hold."""
import contextlib
import io
import json
from pathlib import Path
import tempfile
import unittest

import numpy as np
from scipy.constants import Boltzmann, elementary_charge, hbar
from scipy.integrate import solve_ivp
from scipy.sparse import csr_matrix

from pysnspd.experimental import thermal_spatial_usadel as thermal
from pysnspd.experimental.bulk_current_reference import solve_bulk_reference
from pysnspd.experimental.dc_current_ports import DCCurrentPorts, advance_bulk_contact_phases
from pysnspd.experimental.stable_kwt_euler import incremental_local_update
from pysnspd.solver.core import TDGLSolver
from sandbox.stage5_prephoton.dc_hold import SerialSpectralEvaluator, cut_currents, execute_hold


def end_fixed_graph(nx=5, ny=4, length=4., width=3.):
    raw = thermal.rectangular_graph(nx, ny, length, width)
    x = raw.coordinates_bar[:, 0]
    ends = np.flatnonzero((x == x.min()) | (x == x.max()))
    return thermal.ThermalGraph(raw.area_weights, raw.edges, raw.conductance,
        raw.coordinates_bar, ends)


def plan(bias):
    return dict(schema='pysnspd.stage5.prephoton_dc_hold.v1', T_K=.9, Tc_K=8.65,
        tau_ee_Tc_ps=6., tau_ep_Tc_ps=24.7, dt_ps=.002, duration_ps=.01,
        spectral_tolerance=1e-10, delta0_over_kBTc=1.764, sheet_resistance_ohm=608.,
        gap_drift_relative_limit=.005, gap_ripple_relative_limit=.01,
        current_cut_relative_limit=.02, current_absolute_floor_A=1e-9,
        maximum_relative_displacement=.05, maximum_newton_iterations=40,
        observation_cadence=2, bias_current_A=bias,
        circuit=dict(R_load_ohm=50., R_bias_ohm=1e4, L_bias_H=1e-6,
            Lk_ext_H=10e-9, C_couple_F=100e-12))


class PrephotonDCHoldTests(unittest.TestCase):
    def test_quotient_ports_recover_sheet_resistance_and_no_spurious_side_injection(self):
        graph = end_fixed_graph()
        ports = DCCurrentPorts(graph, Tc_K=8.65, sheet_resistance_ohm=608.)
        imposed = 15e-6
        result = ports.solve(np.zeros(len(graph.edges)), imposed)
        self.assertAlmostEqual(result.Vdev_V, 608.*4/3*imposed, places=13)
        self.assertLess(np.ptp(result.phi_V[ports.left]), 1e-15)
        self.assertLess(np.ptp(result.phi_V[ports.right]), 1e-15)
        _, currents = cut_currents(graph, result.normal_current_bar*ports.current_scale_A)
        np.testing.assert_allclose(currents, imposed, rtol=2e-14, atol=1e-19)
        self.assertLess(np.max(abs(result.current_residual_A)), 1e-18)
        self.assertAlmostEqual(result.normal_joule_W, imposed*result.Vdev_V, places=15)
        np.testing.assert_allclose(result.potential_v*Boltzmann*8.65/(2*elementary_charge), result.phi_V)

    def test_contact_phase_difference_has_the_passive_josephson_voltage_sign(self):
        graph = end_fixed_graph()
        gap = np.ones(graph.n_nodes, complex)
        phi = np.zeros(graph.n_nodes)
        x = graph.coordinates_bar[:, 0]
        phi[x == x.min()] = 1e-6
        v = 2*elementary_charge*phi/(Boltzmann*8.65)
        tD_ps = hbar/(2*Boltzmann*8.65)*1e12
        dt_ps = .1
        advanced = advance_bulk_contact_phases(gap, v, graph.boundary_nodes,
            dt_ps=dt_ps, tD_ps=tD_ps)
        delta_phase = np.angle(advanced[x == x.max()][0]/advanced[x == x.min()][0])
        expected = 2*elementary_charge*1e-6/hbar*dt_ps*1e-12
        self.assertAlmostEqual(delta_phase, expected, places=15)
        np.testing.assert_allclose(abs(advanced), 1., atol=1e-15, rtol=0)

    def test_homogeneous_biased_discrete_state_is_preserved_without_force_subtraction(self):
        graph = end_fixed_graph()
        q = .2; hx = 1.
        # On an orthogonal mesh the exact uniform spectral depairing parameter
        # is 2(1-cos(q*hx))/hx^2, not the continuum q^2 approximation.
        reference = solve_bulk_reference(.9/8.65, 32,
            np.sqrt(2*(1-np.cos(q*hx)))/hx)
        phase = np.exp(1j*q*graph.coordinates_bar[:, 0])
        gap = reference.gap_bar*phase
        u = reference.u[:, None]*phase[None, :]
        test_plan = plan(0.)
        evaluator = SerialSpectralEvaluator(graph, reference.epsilon_bar, u, test_plan)
        value = evaluator(gap)
        I0 = Boltzmann*8.65/(2*elementary_charge*608.)
        _, current = cut_currents(graph, value['current']*I0)
        test_plan['bias_current_A'] = float(current.mean())
        with tempfile.TemporaryDirectory() as directory:
            with contextlib.redirect_stdout(io.StringIO()):
                result = execute_hold(graph, gap, reference.epsilon_bar, test_plan, evaluator,
                    output=Path(directory))
            persisted = json.loads((Path(directory)/'summary.json').read_text())
            self.assertEqual(persisted['checks'], result['checks'])
            with np.load(Path(directory)/'fields.npz') as fields:
                np.testing.assert_allclose(fields['time_ps'][-1], test_plan['duration_ps'])
                self.assertEqual(fields['gap'].shape[-1], graph.n_nodes)
        self.assertTrue(result['dc_hold_passed'])
        self.assertFalse(result['source_force_subtracted'])
        self.assertLess(result['maxima']['gap_drift_relative'], 1e-12)
        self.assertLess(result['maxima']['gap_ripple_relative'], 1e-12)
        self.assertLess(result['maxima']['dark_readout_absolute_V'], 1e-15)
        self.assertLess(result['maxima']['complete_power_residual_W'], 1e-18)
        self.assertFalse(result['nonlinear_heat_closed'])

    def test_moving_spectral_contacts_keep_the_bulk_magnitude_and_phase(self):
        graph = end_fixed_graph()
        reference = solve_bulk_reference(.9/8.65, 16, .1)
        phase = np.exp(.1j*graph.coordinates_bar[:, 0])
        gap = reference.gap_bar*phase
        u = reference.u[:, None]*phase[None, :]
        evaluator = SerialSpectralEvaluator(graph, reference.epsilon_bar, u, plan(0.))
        # A constant phase rotation is a gauge transformation, not a physical
        # disturbance; both spectra and their boundary data must follow it.
        rotated = gap*np.exp(.013j)
        evaluator(rotated)
        np.testing.assert_allclose(evaluator.u[:, graph.boundary_nodes],
            u[:, graph.boundary_nodes]*np.exp(.013j), rtol=1e-13, atol=1e-13)


class StableInheritedEulerTests(unittest.TestCase):
    def test_equivalent_to_original_quadratic_at_moderate_gamma(self):
        rng = np.random.default_rng(23)
        psi = (.7+.3*rng.random(12))*np.exp(1j*rng.random(12))
        force = .02*(rng.normal(size=12)+1j*rng.normal(size=12))
        mu = .1*rng.normal(size=12)
        for gamma in (0., .2, 2., 10.):
            old = TDGLSolver.solve_for_psi_squared(psi=psi, abs_sq_psi=abs(psi)**2,
                forcing_dimensionless=force, mu=mu, gamma=gamma, u=1., dt=.001,
                epsilon=np.zeros(12), psi_laplacian=csr_matrix((12, 12)))
            new = incremental_local_update(psi, force, gamma=gamma,u=1.,dt=.001,mu=mu)
            self.assertIsNotNone(old)
            np.testing.assert_allclose(new[0], old[0], rtol=2e-11, atol=2e-12)
            np.testing.assert_allclose(new[1], old[1], rtol=2e-11, atol=2e-12)

    def test_stationary_state_is_preserved_at_large_gamma_without_freezing_force(self):
        psi = np.array([.8+.3j, 1., -.2+.7j])
        zero = np.zeros(3)
        for gamma in (200., 1000., 1e5):
            new, square = incremental_local_update(psi, zero, gamma=gamma,u=1.,dt=.001,mu=zero)
            np.testing.assert_array_equal(new, psi)
            np.testing.assert_array_equal(square, abs(psi)**2)
        advanced, _ = incremental_local_update(psi, 1e-6*psi, gamma=200.,u=1.,dt=.001,mu=zero)
        self.assertGreater(np.max(abs(advanced-psi)), 1e-12)

    def test_update_satisfies_the_same_quadratic_root_equation(self):
        psi = np.array([.8+.2j, .1-.6j])
        force = np.array([-.03+.04j, .015+.02j])
        mu = np.array([.1, -.2]); dt = .0004
        for gamma in (0., 5., 300.):
            new, rnew = incremental_local_update(psi,force,gamma=gamma,u=1.,dt=dt,mu=mu)
            r = abs(psi)**2
            a = dt*np.sqrt(1+gamma**2*r)*force
            U = np.exp(-1j*mu*dt)
            residual = new+.5*U*gamma**2*psi*(rnew-r)-U*(psi+a)
            self.assertLess(np.max(abs(residual)), 3e-12)
            np.testing.assert_allclose(abs(new)**2, rnew, rtol=1e-12, atol=1e-12)

    def test_radial_observable_converges_at_first_order_under_step_refinement(self):
        gamma = 3.; horizon = .2; kappa = 2.
        oracle = solve_ivp(lambda t,y: -kappa*y/np.sqrt(1+gamma**2*y*y),
            (0.,horizon), [.9], rtol=1e-12, atol=1e-14).y[0,-1]
        errors = []
        for count in (20,40,80):
            psi = np.array([.9+0j]); dt=horizon/count
            for _ in range(count):
                psi,_ = incremental_local_update(psi,-kappa*psi,gamma=gamma,
                    u=1.,dt=dt,mu=np.zeros(1))
            errors.append(abs(abs(psi[0])-oracle))
        self.assertTrue(1.8 < errors[0]/errors[1] < 2.2, errors)
        self.assertTrue(1.8 < errors[1]/errors[2] < 2.2, errors)


if __name__ == '__main__':
    unittest.main()
