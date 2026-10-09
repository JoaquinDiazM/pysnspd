"""Coordinate changes must preserve physical rows and expose inconsistency."""
import unittest
import numpy as np
from pysnspd.experimental.joint_linear_solve import (
    equilibrate_matrix, solve_equilibrated_affine_candidate)


class JointLinearSolveTests(unittest.TestCase):
    def test_scaled_unique_system_retains_the_small_physical_mode(self):
        A = np.diag([1e9, 1e-6, 2.])
        expected = np.array([.2, -.3, .4])
        base = -A@expected
        before = A.copy()
        x, receipt = solve_equilibrated_affine_candidate(A, base,
            relative_rank_tolerance=1e-7, residual_tolerance=1e-7)
        np.testing.assert_allclose(x, expected, rtol=2e-14, atol=2e-14)
        np.testing.assert_array_equal(A, before)
        self.assertEqual(receipt["rank"], 3)
        self.assertTrue(receipt["algebraically_solved_full_rank"])
        self.assertFalse(receipt["temporal_admission"])
        self.assertFalse(receipt["energy_admission"])
        self.assertTrue(receipt["no_regularization"])

    def test_original_rows_are_not_admitted_by_a_small_balanced_residual(self):
        # A complete gauge mode leaves a genuine inconsistent dependent row.
        A = np.array([[1e6, -1e6], [1e6, -1e6]])
        base = np.array([-1e6, -1e6+2.])
        x, receipt = solve_equilibrated_affine_candidate(A, base,
            gauge_direction=np.ones(2), gauge_reference_index=1,
            residual_tolerance=1e-7)
        self.assertEqual(receipt["rank"], 1)
        self.assertTrue(receipt["full_column_rank"])
        self.assertFalse(receipt["original_rows_ok"])
        self.assertFalse(receipt["algebraically_solved_full_rank"])
        self.assertGreater(receipt["original_scaled_residual_max"], .9)
        self.assertLess(receipt["balanced_residual_max"], .01)
        self.assertEqual(receipt["physical_rows"], 2)
        np.testing.assert_allclose(A@x+base, [-1., 1.], atol=1e-8)

    def test_checked_gauge_removes_one_column_and_preserves_all_rows(self):
        A = np.array([[2., -2., 3.], [1., -1., -2.], [3., -3., 1.]])
        representative = np.array([.4, 0., -.2])
        x, receipt = solve_equilibrated_affine_candidate(A, -A@representative,
            gauge_direction=np.array([1., 1., 0.]), gauge_reference_index=1)
        np.testing.assert_allclose(x, representative, atol=3e-15)
        self.assertEqual(receipt["physical_rows"], 3)
        self.assertEqual(receipt["physical_columns"], 2)
        self.assertTrue(receipt["no_physical_rows_removed"])
        self.assertTrue(receipt["gauge_reduced"])
        wrong = A.copy()
        wrong[2, 1] += .1
        with self.assertRaisesRegex(FloatingPointError, "Ward"):
            solve_equilibrated_affine_candidate(wrong, -A@representative,
                gauge_direction=np.array([1., 1., 0.]), gauge_reference_index=1)

    def test_an_additional_null_mode_stays_deficient(self):
        A = np.array([[1., 2., -3.], [2., 4., -6.], [0., 0., 0.]])
        x, receipt = solve_equilibrated_affine_candidate(A, np.array([1., 2., 0.]),
            gauge_direction=np.ones(3), gauge_reference_index=2)
        self.assertEqual(receipt["rank"], 1)
        self.assertEqual(receipt["physical_columns"], 2)
        self.assertFalse(receipt["algebraically_solved_full_rank"])
        self.assertFalse(receipt["full_column_rank"])
        np.testing.assert_allclose(A@x+[1., 2., 0.], 0., atol=1e-14)

    def test_equilibration_identity_keeps_zero_rows_and_columns(self):
        A = np.array([[0., 0., 0.], [1e8, 1e-8, 0.], [2e8, -3e-8, 0.]])
        B, row, column = equilibrate_matrix(A)
        np.testing.assert_allclose(B*row[:, None]*column[None, :], A, rtol=5e-15)
        self.assertEqual(row[0], 1.)
        self.assertEqual(column[2], 1.)
        np.testing.assert_array_equal(B[0], 0.)
        np.testing.assert_array_equal(B[:, 2], 0.)

    def test_invalid_input_rejected_without_silent_real_projection(self):
        for matrix, base in [(np.eye(2)*1j, np.zeros(2)),
                (np.eye(2), np.zeros(3)), (np.eye(2), [0., np.nan])]:
            with self.assertRaises(ValueError):
                solve_equilibrated_affine_candidate(matrix, base)
        with self.assertRaises(ValueError):
            equilibrate_matrix(np.eye(2), iterations=True)
        with self.assertRaises(ValueError):
            solve_equilibrated_affine_candidate(np.eye(2), np.zeros(2),
                gauge_reference_index=0)


if __name__ == "__main__":
    unittest.main()
