"""Test the real thesis update against the admitted continuous thermal law."""
import unittest
from sandbox.stage4_core.kwt_bridge_check import bridge_diagnostic


class InheritedKWTBridgeTests(unittest.TestCase):
    def test_real_local_solver_has_correct_continuous_limit_and_fixed_contacts(self):
        result = bridge_diagnostic()
        self.assertEqual(result['status'], 'PASSED', result)
        self.assertEqual(result['direct_inherited_solver_calls'], 3)
        self.assertLessEqual(result['maximum_spectral_residual'], 1e-7)


if __name__ == '__main__':
    unittest.main()
