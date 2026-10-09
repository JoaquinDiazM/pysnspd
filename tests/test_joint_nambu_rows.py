"""Independent SI integrals, current orientation, directional and gauge controls."""
import unittest
import numpy as np
from scipy.constants import elementary_charge
from scipy.sparse import csr_matrix
from numpy.polynomial.legendre import leggauss
from pysnspd.experimental import adiabatic_matrix_usadel as ad
from pysnspd.experimental.nambu_moments import NambuMoments
from pysnspd.experimental.joint_nambu_rows import JointNambuRows, checked_gauge_reduction


class JointNambuRowsTests(unittest.TestCase):
    def setUp(self):
        points, weights = leggauss(128)
        self.Eu, self.N0, self.I0 = 2e-22, 8e46, 4e-6
        self.E, self.w = 2.5*(points+1), 2.5*weights
        self.volume = np.array([2e-24, 3e-24, 4e-24])
        self.moments = NambuMoments(self.E*self.Eu, self.w*self.Eu,
            self.volume, self.N0, self.Eu, self.I0)
        self.rows = JointNambuRows(self.moments, [[0, 1], [1, 2]], left_nodes=[0], right_nodes=[2])
        self.zeros = np.zeros(self.moments.shape, complex)
        self.flux = np.zeros((len(self.E), 2, 2, 2), complex)
        h0 = np.tanh(self.E/.4)
        self.reference = np.broadcast_to(2*h0[:, None, None, None]*ad.TAU3,
            self.moments.shape).copy()
        self.arguments = dict(K_zero=self.reference, K_first=self.zeros,
            K_thermal_zero=self.reference, phi_V=np.zeros(3),
            thermal_gradient_J_per_gapbar=np.zeros(3, complex),
            thermal_edge_current_A=np.full(2, 21e-6), edge_flux_zero=self.flux,
            edge_flux_first=self.flux, thermal_edge_flux_zero=self.flux,
            imposed_current_A=21e-6,
            reference_contract='same instantaneous gap, links and gauge; independent Matsubara baseline')

    def test_independent_pair_charge_and_completed_current_integrals(self):
        zero = self.reference.copy()
        zero += .04*ad.I2
        zero[:, :, 0, 1] = .2+.3j
        first = self.zeros.copy()
        first[:, :, 0, 1] = -.07+.02j
        flux0, flux1 = self.flux.copy(), self.flux.copy()
        a, b = np.array([.01, .03]), np.array([.004, -.002])
        flux0 += a[None, :, None, None]*ad.TAU3
        flux1 += b[None, :, None, None]*ad.TAU3
        phi = np.full(3, -.02*5*self.Eu/elementary_charge)
        thermal = np.array([1e-22+2e-22j, 0, -1e-22])
        result = self.rows.evaluate(**dict(self.arguments, K_zero=zero, K_first=first,
            edge_flux_zero=flux0, edge_flux_first=flux1, phi_V=phi,
            thermal_gradient_J_per_gapbar=thermal))
        expected = thermal+1j*self.N0*self.Eu*self.volume*5*self.Eu*(.13+.32j)
        np.testing.assert_allclose(result.gap_force_J_per_gapbar, expected, rtol=3e-15)
        np.testing.assert_allclose(result.neutrality_V, 0, atol=2e-18)
        currents = 21e-6+5*self.I0*(a+b)
        np.testing.assert_allclose(result.edge_current_A, currents, rtol=3e-15)
        np.testing.assert_allclose(result.port_current_residual_A, currents-21e-6, atol=1e-20)
        np.testing.assert_allclose(result.interior_current_divergence_A, [currents[1]-currents[0]], atol=1e-20)
        np.testing.assert_allclose(result.thermal_real_axis_discrepancy_A, [-21e-6]*2)
        self.assertEqual(result.port_voltage_V, 0.)

    def test_current_reference_and_orientation_do_not_add_an_ohmic_face(self):
        # Equilibrium subtraction restores the independent baseline exactly.
        thermal = self.flux.copy()
        thermal += .3*ad.TAU3
        result = self.rows.evaluate(**dict(self.arguments, edge_flux_zero=thermal,
            thermal_edge_flux_zero=thermal, phi_V=np.array([.002, .001, 0.])))
        np.testing.assert_array_equal(result.edge_current_A, [21e-6, 21e-6])
        np.testing.assert_array_equal(result.port_current_residual_A, [0., 0.])
        self.assertEqual(result.port_voltage_V, .002)
        reverse = JointNambuRows(self.moments, [[1, 0], [2, 1]], left_nodes=[0], right_nodes=[2])
        flipped = reverse.evaluate(**dict(self.arguments,
            thermal_edge_current_A=np.full(2, -21e-6)))
        np.testing.assert_array_equal(flipped.port_current_residual_A, [0., 0.])

    def test_direction_matches_independent_affine_difference(self):
        rng = np.random.default_rng(42017)
        k1 = rng.normal(size=self.zeros.shape)+1j*rng.normal(size=self.zeros.shape)
        flux1 = rng.normal(size=self.flux.shape)+1j*rng.normal(size=self.flux.shape)
        phi = rng.normal(size=3)*1e-3
        baseline = self.rows.evaluate(**self.arguments)
        direction = self.rows.direction(K_first_direction=k1,
            edge_flux_first_direction=flux1, phi_direction_V=phi,
            imposed_current_direction_A=2e-6,
            reference_contract=self.arguments['reference_contract'])
        step = .17
        moved = self.rows.evaluate(**dict(self.arguments, K_first=step*k1,
            edge_flux_first=step*flux1, phi_V=step*phi,
            imposed_current_A=21e-6+step*2e-6))
        for name in ('neutrality_V', 'gap_force_J_per_gapbar', 'edge_current_A',
                     'port_current_residual_A', 'interior_current_divergence_A'):
            np.testing.assert_allclose((getattr(moved, name)-getattr(baseline, name))/step,
                getattr(direction, name), rtol=2e-13, atol=1e-30)

    def test_physical_density_scales_remove_volume_but_not_rows(self):
        first = self.zeros.copy()
        first[:, :, 0, 1] = .1+.2j
        direction = self.rows.direction(K_first_direction=first,
            edge_flux_first_direction=self.flux, phi_direction_V=np.zeros(3),
            reference_contract=self.arguments['reference_contract'])
        scaled = direction.pack(neutrality_nodes=[0, 1, 2], force_nodes=[0, 1, 2],
            scales=self.rows.scales, include_port_currents=True)
        expected_force = np.tile([-1., .5], 3)
        np.testing.assert_allclose(scaled[3:9], expected_force, rtol=3e-15)
        self.assertEqual(len(scaled), 11)
        expanded = direction.pack(neutrality_nodes=[], force_nodes=[],
            scales=self.rows.scales, include_port_currents=False,
            include_interior_current=True)
        self.assertEqual(len(expanded), 1)

    def test_temporal_gauge_neutrality_uses_distribution_correction(self):
        kappa, chi_t = .07, np.array([.13, -.04, .02])
        h0 = np.tanh(self.E/.4)
        derivative = (1-h0*h0)/.4
        jet = ad.Jet.of(self.reference, energy=np.broadcast_to(
            2*derivative[:, None, None, None]*ad.TAU3, self.moments.shape))
        U = ad.links(np.zeros(3), -chi_t)
        moved = ad.transport(U, jet, kappa)
        phi = -kappa*chi_t*self.Eu/elementary_charge
        value = self.rows.evaluate(**dict(self.arguments, K_zero=moved.value,
            K_first=moved.first, phi_V=phi))
        np.testing.assert_allclose(value.neutrality_V, 0., atol=1e-12*self.Eu/elementary_charge)
        wrong = self.rows.evaluate(**dict(self.arguments, phi_V=phi))
        self.assertGreater(np.linalg.norm(wrong.neutrality_V), 1e-7)

    def test_equipotential_residuals_do_not_hide_cut_node_disagreement(self):
        moments = NambuMoments(self.E*self.Eu, self.w*self.Eu,
            np.full(4, 2e-24), self.N0, self.Eu, self.I0)
        rows = JointNambuRows(moments, [[0, 2], [1, 3], [0, 1], [2, 3]],
            left_nodes=[0, 1], right_nodes=[2, 3])
        zero = np.zeros(moments.shape, complex)
        flux = np.zeros((len(self.E), 4, 2, 2), complex)
        value = rows.evaluate(K_zero=zero, K_first=zero, K_thermal_zero=zero,
            phi_V=np.array([.001, .002, 0., .0004]),
            thermal_gradient_J_per_gapbar=np.zeros(4, complex),
            thermal_edge_current_A=np.array([2e-6, 3e-6, 1e-6, 1e-6]),
            edge_flux_zero=flux, edge_flux_first=flux, thermal_edge_flux_zero=flux,
            imposed_current_A=5e-6, reference_contract='same spectrum and gauge')
        np.testing.assert_allclose(value.port_current_residual_A, 0., atol=1e-20)
        self.assertAlmostEqual(value.port_voltage_V, .0013)
        np.testing.assert_allclose(value.left_equipotential_residual_V, [-.0005, .0005])
        np.testing.assert_allclose(value.right_equipotential_residual_V, [-.0002, .0002])

    def test_bad_topology_and_empty_reference_rejected(self):
        with self.assertRaises(ValueError):
            JointNambuRows(self.moments, [[0, 0]], left_nodes=[0], right_nodes=[2])
        with self.assertRaises(ValueError):
            JointNambuRows(self.moments, [[0, 1]], left_nodes=[0], right_nodes=[0])
        with self.assertRaises(ValueError):
            self.rows.evaluate(**dict(self.arguments, reference_contract=''))


class CheckedGaugeTests(unittest.TestCase):
    def setUp(self):
        rng = np.random.default_rng(271828)
        self.columns = np.array([1e12, 1e-3, 1e12, 1e-3])
        self.rows = np.array([1e-3, 1e-21, 1e-6, 1e-3, 1e-21])
        self.q = np.array([2., -1., 3., 4.])
        A = rng.normal(size=(5, 4))
        A[:, 3] = -(A[:, :3]@self.q[:3])/self.q[3]
        self.A = A*self.rows[:, None]/self.columns[None, :]
        self.g = self.q*self.columns

    def test_one_gauge_removed_without_changing_any_physical_row(self):
        for matrix in (self.A, csr_matrix(self.A)):
            reduction = checked_gauge_reduction(matrix, self.g, reference_index=3,
                row_scales=self.rows, column_scales=self.columns)
            self.assertEqual(reduction.matrix.shape, (5, 3))
            value = self.columns*np.array([.4, .7, -.2, .9])
            representative = reduction.representative(value)
            self.assertEqual(representative[3], 0.)
            np.testing.assert_allclose((self.A@representative)/self.rows,
                (self.A@value)/self.rows, rtol=1e-14, atol=2e-15)
            np.testing.assert_array_equal(reduction.expand(representative[[0, 1, 2]]), representative)

    def test_tiny_si_gap_row_cannot_be_ignored_by_a_voltage_norm(self):
        wrong = self.A.copy()
        wrong[1, 3] += self.rows[1]/self.columns[3]
        with self.assertRaisesRegex(FloatingPointError, 'Ward'):
            checked_gauge_reduction(wrong, self.g, reference_index=3,
                row_scales=self.rows, column_scales=self.columns)

    def test_frozen_potential_alone_is_not_a_gauge_vector(self):
        wrong = np.zeros(4)
        wrong[3] = 1.
        with self.assertRaisesRegex(FloatingPointError, 'Ward'):
            checked_gauge_reduction(self.A, wrong, reference_index=3,
                row_scales=self.rows, column_scales=self.columns)

    def test_extra_physical_nullspace_is_preserved(self):
        # One gauge may be removed, but a zero physical rate row is not repaired.
        reduction = checked_gauge_reduction(np.zeros((2, 3)), np.array([1., 1., 0.]),
            reference_index=0, row_scales=np.ones(2), column_scales=np.ones(3))
        self.assertEqual(np.linalg.matrix_rank(reduction.matrix), 0)


if __name__ == '__main__':
    unittest.main()
