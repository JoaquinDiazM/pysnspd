"""Manufactured mixed-dimensional algebra; no material or transient admission."""
from types import SimpleNamespace
import unittest

import numpy as np

from pysnspd.experimental.energy_catalog import BCS_GAP_RATIO, HBAR_J_S, K_B_J_K
from pysnspd.experimental.mixed_spatial import MixedSpatialFunctional
from pysnspd.experimental.spatial_functional import SpatialAdmissibilityError
from pysnspd.experimental.spatial_open import OpenSpatialFunctional
from pysnspd.experimental.electrical_ports import solve_potential


class SyntheticCatalogue:
    def __init__(self):
        self.count_nodes = np.array([.1, .5, 1.2])
        self.count_weights = np.array([.2, .3, .5])
        tc = 8.65
        self.vacuum = SimpleNamespace(metadata={"Tc_K": tc}, D_m2_s=1.58e-4,
            N0_per_J_m3=2.1e47, delta0_J=BCS_GAP_RATIO*K_B_J_K*tc,
            gamma_axis=np.array([0., 2.]), delta_axis=np.array([0., 2.]), evaluate=self.vacuum_values)

    @staticmethod
    def vacuum_values(a, gamma):
        if not 0 <= a <= 2 or not 0 <= gamma <= 2:
            raise ValueError("synthetic support")
        return (.5*(a*a-1)**2+.3*gamma*a*a+.07*gamma**2,
                2*a*(a*a-1)+.6*gamma*a, .3*a*a+.14*gamma)

    def energy_kernel(self, a, gamma):
        self.vacuum_values(a, gamma)
        root = np.sqrt(self.count_nodes**2+a*a)
        alpha = .13+.03*self.count_nodes
        return (root+alpha*gamma*(1+.2*a*a), a/root+.4*alpha*gamma*a, alpha*(1+.2*a*a))


def gradient(result):
    return result.gradient_cartesian_bar[:, 0]+1j*result.gradient_cartesian_bar[:, 1]


def fd5(function, h=2e-5):
    return (function(-2*h)-8*function(-h)+8*function(h)-function(2*h))/(12*h)


class TestMixedSpatial(unittest.TestCase):
    def setUp(self):
        self.catalog = SyntheticCatalogue()
        self.ell = np.sqrt(HBAR_J_S*self.catalog.vacuum.D_m2_s/(2*K_B_J_K*8.65))
        self.model = MixedSpatialFunctional(self.catalog, 8*self.ell, 4*self.ell, 7e-9,
            4*self.ell, 3*self.ell, elements_x=1, elements_y=1, left_elements=1, right_elements=1)
        m = self.model
        x, y = m.dof_coordinates_bar.T
        interior = (x > 0)&(x < 8)
        envelope = np.where(interior, x*(8-x)/16, 0.)
        self.z = (.9+.012*np.cos(x/3)+.008*envelope*y)*np.exp(1j*(.025*x+.012*envelope*y))
        self.links = .006*np.sin(np.arange(len(m.graph_edges))*.7)
        self.p = np.array([(.045+.0002*i)*np.exp(-self.catalog.count_nodes)
                           for i in range(m.quadrature_size)])

    def evaluate(self, z=None, links=None, p=None, stability=False):
        return self.model.evaluate(self.z if z is None else z,
            self.p if p is None else p, link_phases=self.links if links is None else links,
            require_stability=stability)

    def test_positive_measures_and_one_hot_trace_constraint(self):
        m = self.model
        self.assertTrue(np.all(m.mass_bar > 0))
        self.assertTrue(np.all(m.quadrature_mass_bar > 0))
        self.assertAlmostEqual(m.mass_bar.sum(), 15., places=12)
        self.assertAlmostEqual(m.quadrature_mass_bar.sum(), 15., places=12)
        self.assertAlmostEqual(m.node_volumes_m3.sum()/(m.cross_section_m2*self.ell), 15., places=12)
        np.testing.assert_array_equal(np.diff(m.prolongation.indptr), np.ones(m.quadrature_size, dtype=int))
        qfield = m.prolongation@self.z
        np.testing.assert_array_equal(qfield[m.rectangle_quadrature[0]],
                                      np.full(m.rectangle_shape[1], self.z[m.interface_dofs["left"]]))
        np.testing.assert_array_equal(qfield[m.rectangle_quadrature[-1]],
                                      np.full(m.rectangle_shape[1], self.z[m.interface_dofs["right"]]))
        # Quadrature endpoint weights are adjacent subvolume contributions.
        np.testing.assert_allclose(m.mass_bar, m.prolongation.T@m.quadrature_mass_bar, atol=0.)
        self.assertFalse(np.any(m.graph_edges[:, 0] == m.graph_edges[:, 1]))

    def test_tensor_polynomial_derivatives_and_positive_compatible_stiffness(self):
        m = self.model
        x, y = m.dof_coordinates_bar.T
        # Interior rectangle polynomial vanishes on the two restricted traces.
        field = np.where((x >= 0)&(x <= 8), .9+.001*x*(8-x)*y, .9)
        fields = m.sample_fields(field)
        qx, qy = m.quadrature_coordinates_bar[m.sector_slices["rectangle"]].T
        measured = fields.derivative_quadrature_bar[m.sector_slices["rectangle"]]
        np.testing.assert_allclose(measured[:, 0], .001*(8-2*qx)*qy, atol=5e-15)
        np.testing.assert_allclose(measured[:, 1], .001*qx*(8-qx), atol=5e-15)
        _, zq, derivatives, stiffnesses = m._operators(self.z, self.links)
        for d, k in zip(derivatives, stiffnesses):
            compatible = k-d.conj().T@(m.quadrature_mass_bar[:, None]*d)
            self.assertGreaterEqual(np.linalg.eigvalsh(compatible)[0], -1e-12)
        self.assertGreaterEqual(self.evaluate().gradient_remainder_bar, 0.)
        assembled = (m.prolongation.T@(stiffnesses[0]+stiffnesses[1]))@m.prolongation
        self.assertGreaterEqual(np.linalg.eigvalsh(assembled)[0], -1e-12)

    def test_all_interface_forces_and_selected_bulk_forces_differentiate_energy(self):
        m = self.model
        result = self.evaluate()
        original_p = self.p.copy()
        for node in (0, m.interface_dofs["left"], m.interface_dofs["right"], m.cells//2, m.cells-1):
            for direction in (1., 1j):
                variation = np.zeros(m.cells, complex)
                variation[node] = direction
                numerical = fd5(lambda t: self.evaluate(z=self.z+t*variation).energy_bar)
                analytic = np.real(np.vdot(gradient(result), variation))
                self.assertAlmostEqual(numerical, analytic, delta=2e-9)
        raw = result.quadrature_gradient_cartesian_bar
        np.testing.assert_allclose(result.gradient_cartesian_bar, m.prolongation.T@raw, atol=1e-15)
        r = result.interface_reactions_bar
        self.assertAlmostEqual(r["left_lead"]+r["left_rectangle"], gradient(result)[m.interface_dofs["left"]])
        self.assertAlmostEqual(r["right_lead"]+r["right_rectangle"], gradient(result)[m.interface_dofs["right"]])
        np.testing.assert_array_equal(self.p, original_p)

    def test_x_y_lead_currents_are_energy_derivatives_and_noether_includes_restriction(self):
        m = self.model
        result = self.evaluate()
        xlead = m._lines[0][2][0]
        xrectangle = m._lines[2][2][1]
        yrectangle = next(ids[1] for _, axis, ids in m._lines if axis == 1 and len(ids))
        for edge in (xlead, xrectangle, yrectangle, len(m.graph_edges)-1):
            variation = np.zeros(len(m.graph_edges))
            variation[edge] = 1
            numerical = fd5(lambda t: self.evaluate(links=self.links+t*variation).energy_bar)
            self.assertAlmostEqual(numerical, result.link_current_bar[edge], delta=2e-9)
        np.testing.assert_allclose(result.noether_residual_bar, 0., atol=2e-14)
        np.testing.assert_allclose(result.line_phase_residuals_bar, 0., atol=2e-14)
        self.assertAlmostEqual(result.global_phase_residual_bar, 0., delta=2e-14)

    def test_compatible_local_gauge_changes_preserve_current_and_force(self):
        m = self.model
        chi = .17*np.cos(np.arange(m.cells)*.9)
        newlinks = self.links+chi[m.graph_edges[:, 0]]-chi[m.graph_edges[:, 1]]
        first = self.evaluate()
        transformed = self.evaluate(z=self.z*np.exp(1j*chi), links=newlinks)
        self.assertAlmostEqual(first.energy_bar, transformed.energy_bar, delta=8e-14)
        np.testing.assert_allclose(gradient(transformed), gradient(first)*np.exp(1j*chi), atol=8e-14)
        np.testing.assert_allclose(transformed.link_current_bar, first.link_current_bar, atol=7e-14)
        np.testing.assert_allclose(transformed.sampled_fields.gamma_quadrature_bar,
                                   first.sampled_fields.gamma_quadrature_bar, atol=3e-15)

    def test_linear_real_field_interface_tractions_cancel_with_correct_section_weights(self):
        m = self.model
        x = m.dof_coordinates_bar[:, 0]
        field = .85+.005*x
        result = m.evaluate(field, np.zeros_like(self.p), require_stability=False)
        qfield = m.prolongation@field
        local = m.prolongation.T@(m.quadrature_mass_bar*2*qfield*(qfield*qfield-1))
        residual = gradient(result)-local
        expected = np.zeros(m.cells)
        expected[0], expected[-1] = -2*m.kappa*.005, 2*m.kappa*.005
        np.testing.assert_allclose(residual, expected, atol=2e-14)
        for side in ("left", "right"):
            self.assertAlmostEqual(residual[m.interface_dofs[side]], 0., delta=2e-14)
        self.assertAlmostEqual(result.gradient_remainder_bar, m.kappa*.005**2*15., delta=1e-14)

    def test_electrical_graph_has_uniform_wire_solution_and_same_port_work(self):
        m = self.model
        injection = np.zeros(m.cells)
        injection[0], injection[-1] = 30e-6, -30e-6
        sigma = 2.5e5
        solution = solve_potential(m.graph_edges, m.conductance(sigma), np.zeros(len(m.graph_edges)),
                                    injection, left_node=0, right_node=m.cells-1)
        expected = 30e-6*(3-m.dof_coordinates_bar[:, 0]+8)*m.ell0_m/(sigma*m.cross_section_m2)
        np.testing.assert_allclose(solution.phi_V, expected, atol=5e-17)
        self.assertLess(np.max(abs(solution.current_residual_A)), 1e-18)
        # This is floating point algebra, scaled by its nonzero physical power.
        self.assertLess(abs(solution.port_power_W-solution.normal_joule_W),
                        64*np.finfo(float).eps*solution.port_power_W)
        for invalid in (0., -1., np.nan, 1j, np.ones(2)):
            with self.assertRaises(ValueError):
                m.conductance(invalid)

    def test_thermal_envelope_is_fixed_population_force_and_free_energy_derivative(self):
        m = self.model
        result = m.evaluate_thermal(self.z, .25, link_phases=self.links, require_stability=False)
        fixed = self.evaluate(p=result.p_quadrature)
        np.testing.assert_allclose(fixed.gradient_cartesian_bar, result.gradient_cartesian_bar, atol=2e-14)
        np.testing.assert_allclose(fixed.link_current_bar, result.link_current_bar, atol=2e-14)
        dz = .1*np.cos(np.arange(m.cells))+.1j*np.sin(np.arange(m.cells))
        dl = .05*np.cos(np.arange(len(m.graph_edges)))
        def varied(t):
            return m.evaluate_thermal(self.z+t*dz, .25, link_phases=self.links+t*dl, require_stability=False)
        expected = np.real(np.vdot(gradient(result), dz))+np.dot(result.link_current_bar, dl)
        self.assertAlmostEqual(fd5(lambda t: varied(t).free_energy_bar), expected, delta=2e-9)
        self.assertGreater(abs(fd5(lambda t: varied(t).energy_bar)-expected), 1e-5)

    def test_local_2D_principal_symbol_matches_independent_gradient_hessian(self):
        m = self.model
        z = .88+.04j
        d = np.array([.015+.025j, -.01+.018j])
        p = np.array([.04, .03, .01])
        symbol = m.principal_symbol(z, d, p)
        self.assertEqual(symbol.matrix.shape, (4, 4))
        self.assertTrue(symbol.stable)
        x = np.array([d[0].real, d[0].imag, d[1].real, d[1].imag])
        def energy(v):
            derivatives = v[::2]+1j*v[1::2]
            q = np.imag(np.conj(z)*derivatives)/(abs(z)**2+.01)
            gamma = np.dot(q, q)/m.gap_ratio
            u = self.catalog.vacuum.evaluate(abs(z), gamma)[0]
            u += 4*np.dot(self.catalog.count_weights*p, self.catalog.energy_kernel(abs(z), gamma)[0])
            return u+m.kappa*(np.vdot(derivatives, derivatives).real-abs(z)**2*np.dot(q, q))
        h = 1e-4
        eye = h*np.eye(4)
        numerical = np.array([[(energy(x+ei+ej)-energy(x+ei-ej)-energy(x-ei+ej)+energy(x-ei-ej))/(4*h*h)
                                for ej in eye] for ei in eye])
        np.testing.assert_allclose(symbol.matrix, numerical, atol=3e-7, rtol=3e-7)
        open_model = OpenSpatialFunctional(self.catalog, 8*self.ell, m.cross_section_m2, 1)
        one = m.principal_symbol(z, d[:1], p, spatial_dimensions=1)
        reference = open_model.principal_symbol(z, d[0], p)
        np.testing.assert_allclose(one.matrix, reference.matrix, atol=1e-15)

    def test_symbols_cover_each_sector_and_explicit_negative_control_rejects(self):
        m = self.model
        events = []
        result = m.evaluate(self.z, self.p, link_phases=self.links,
                            on_quadrature=lambda i, total: events.append((i, total)))
        self.assertEqual(events, [(i, m.quadrature_size) for i in range(m.quadrature_size)])
        self.assertEqual([s.spatial_dimensions for s in result.principal_symbols], m.spatial_dimensions.tolist())
        self.catalog.vacuum.evaluate = lambda a, g: (-10*g, 0., -10.)
        negative = m.principal_symbol(.9, np.zeros(2), np.zeros(3))
        self.assertEqual(negative.matrix.shape, (4, 4))
        self.assertFalse(negative.stable)
        self.assertLess(negative.eigenvalues[0], 0.)
        with self.assertRaises(SpatialAdmissibilityError):
            m.evaluate(np.full(m.cells, .9+0j), np.zeros_like(self.p))

    def test_invalid_states_and_implicit_population_remapping_are_rejected(self):
        m = self.model
        with self.assertRaises(ValueError):
            m.evaluate(self.z, np.zeros((m.cells, len(self.catalog.count_nodes))))
        for bad in (-1e-12, 1+1e-12, np.nan, 1j):
            p = self.p.astype(complex) if isinstance(bad, complex) else self.p.copy()
            p[0, 0] = bad
            with self.assertRaises(ValueError):
                self.evaluate(p=p)
        for links in (np.zeros(len(m.graph_edges)+1), np.full(len(m.graph_edges), np.inf), 1j*self.links):
            with self.assertRaises(ValueError):
                self.evaluate(links=links)
        for temperature in (-1., np.nan, np.inf):
            with self.assertRaises(ValueError):
                m.evaluate_thermal(self.z, temperature)


if __name__ == "__main__":
    unittest.main()
