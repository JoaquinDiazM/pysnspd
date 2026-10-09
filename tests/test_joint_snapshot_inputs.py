"""Pinned cache reuse, changed-producer provenance and current-root checks."""
import hashlib
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch
import numpy as np
from pysnspd.experimental import adiabatic_matrix_usadel as ad
from pysnspd.experimental.thermal_spatial_usadel import rectangular_graph
from sandbox.stage5_photon import joint_snapshot_inputs as inputs


def homogeneous_matrix(nodes, gap, energy, eta):
    z = eta-1j*energy
    root = np.sqrt(z*z+gap*gap)
    value = (z*ad.TAU3+gap*np.array([[0., 1.], [1., 0.]]))/root
    return np.broadcast_to(value, (nodes, 2, 2)).copy()


class JointInputPinTests(unittest.TestCase):
    def test_data_pin_remains_exact_when_a_producer_changes(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            data = root/'numeric.bin'
            producer = root/'producer.py'
            data.write_bytes(b'unchanged numerical fields')
            producer.write_bytes(b'new current producer')
            digest = hashlib.sha256(data.read_bytes()).hexdigest()
            ledger = {}
            inputs._verify_pin(data, digest, ledger)
            with patch.object(inputs, 'ROOT', root):
                history = inputs._producer_provenance({'producer.py': 'old producer pin'}, 'fixture')
            self.assertFalse(history[0]['unchanged'])
            self.assertEqual(history[0]['saved_sha256'], 'old producer pin')
            self.assertEqual(history[0]['current_sha256'], hashlib.sha256(producer.read_bytes()).hexdigest())
            self.assertEqual(ledger[str(data)], digest)
            data.write_bytes(b'changed numerical fields')
            with self.assertRaisesRegex(ValueError, 'Pinned'):
                inputs._verify_pin(data, digest, {})

    def test_missing_producer_is_visible_provenance(self):
        with tempfile.TemporaryDirectory() as directory, patch.object(inputs, 'ROOT', Path(directory)):
            history = inputs._producer_provenance({'missing.py': 'saved'}, 'fixture')
        self.assertFalse(history[0]['unchanged'])
        self.assertIsNone(history[0]['current_sha256'])


class CurrentCachedRootTests(unittest.TestCase):
    def setUp(self):
        self.graph = rectangular_graph(4, 3, 2., 1.5)
        self.gap = np.full(self.graph.n_nodes, 1.2, complex)
        self.alpha = np.zeros(len(self.graph.edges))
        self.Eu = 2e-22
        self.energy, self.eta = .8, .03
        self.R = homogeneous_matrix(self.graph.n_nodes, 1.2, self.energy, self.eta)

    def evaluate(self, value):
        return inputs._validate_cached_R0(self.graph, self.gap, self.alpha, value,
            self.energy*self.Eu, self.eta*self.Eu, self.Eu)

    def test_unchanged_causal_spectrum_passes_current_equations(self):
        result = self.evaluate(self.R)
        self.assertLess(result['matrix_per_mass_max'], 1e-13)
        self.assertLess(result['chart_scaled_max'], 1e-13)
        self.assertLess(result['normalization_max'], 1e-13)
        self.assertGreater(result['minimum_DOS'], 0.)

    def test_wrong_current_physical_gap_rejects_stale_root(self):
        with self.assertRaisesRegex(ValueError, 'spectral'):
            inputs._validate_cached_R0(self.graph, self.gap*.97, self.alpha, self.R,
                self.energy*self.Eu, self.eta*self.Eu, self.Eu)

    def test_noncausal_root_is_rejected_without_dos_clipping(self):
        with self.assertRaisesRegex(ValueError, 'DOS'):
            self.evaluate(-self.R)

    def test_discarded_full_matrix_component_is_detected(self):
        wrong = self.R.copy()
        wrong[:, 1, 1] += .1
        with self.assertRaisesRegex(ValueError, 'stereographic'):
            self.evaluate(wrong)

    def test_static_spatial_gauge_keeps_same_physical_root(self):
        rng = np.random.default_rng(5109)
        chi = rng.normal(size=self.graph.n_nodes)*.3
        U = ad.links(-chi).value
        moved = U@self.R@U.conj().transpose(0, 2, 1)
        tail, head = self.graph.edges.T
        alpha = chi[head]-chi[tail]
        result = inputs._validate_cached_R0(self.graph, self.gap*np.exp(1j*chi),
            alpha, moved, self.energy*self.Eu, self.eta*self.Eu, self.Eu)
        self.assertLess(result['matrix_per_mass_max'], 1e-13)


class ActualInitialCacheTests(unittest.TestCase):
    @unittest.skipUnless(Path('/home/jdiaz/pysnspd/docs/implementation/stage5/photon_readiness_20261008/fd_be_775_variance_1p0_v2_uv24/receipt.json').exists(),
        'Optional remote numeric fixture is not shipped with the source repository')
    def test_all_112_actual_roots_and_pins_are_reused_without_recomputation(self):
        data = inputs.load_joint_initial('/home/jdiaz/pysnspd/docs/implementation/stage5/photon_readiness_20261008/fd_be_775_variance_1p0_v2_uv24')
        provenance = data['provenance']
        self.assertEqual(provenance['full_energy_count'], 112)
        self.assertEqual(len(provenance['current_R0_validation']), 112)
        self.assertGreater(len(provenance['pinned_files_sha256']), 112)
        self.assertFalse(provenance['source_reapplied'])
        self.assertEqual(provenance['spectral_roots_repeated'], 0)
        self.assertEqual(provenance['physical_time_steps'], 0)
        self.assertEqual(data['roots_R0'].shape, (112, 1712, 2, 2))
        self.assertEqual(data['inputs']['observed']['count'], 256)


if __name__ == '__main__':
    unittest.main()
