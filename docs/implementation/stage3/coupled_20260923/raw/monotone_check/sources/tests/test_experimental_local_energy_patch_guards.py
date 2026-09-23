"""Cheap independent guards and smooth oracles for the local energy patch.

Analytic spectra below are synthetic test data, not material admission.
No causal catalogue is evaluated and no trajectory is integrated.
"""
from pathlib import Path
from types import SimpleNamespace
import tempfile
import unittest

import numpy as np

from pysnspd.experimental.local_energy_patch import LocalEnergyPatch


class SmoothSource:
    def __init__(self):
        self.count_nodes = np.array([.03, .12, .5, 1.4, 3.])
        self.count_weights = np.array([.04, .1, .3, .6, .8])
        self.eta = 1e-5
        self.samples = []
        self.vacuum = SimpleNamespace(delta0_J=1., N0_per_J_m3=2., D_m2_s=3.,
            evaluate=lambda a, g: (.5*a*a+g*g, a, 2*g))

    def energy_kernel(self, a, g):
        if not .985 <= a <= 1.005 or not .003 <= g <= .009:
            raise ValueError('outside synthetic box')
        self.samples.append((a, g))
        root = np.sqrt(self.count_nodes**2+a*a)
        factor = (.03+.002*self.count_nodes)*np.exp(g/.05)
        return root+factor*a**3, a/root+3*factor*a*a, factor*a**3/.05


def make_patch(degree=4):
    source = SmoothSource()
    return source, LocalEnergyPatch.build(source, (.985, 1.005), (.003, .009), degree)


class TestLocalEnergyPatchAdditionalGuards(unittest.TestCase):
    def test_off_grid_smooth_values_and_both_scaled_derivatives(self):
        # Independent analytic values at points that are not construction nodes.
        source, patch = make_patch(6)
        for a, g in ((.99137, .00611), (.98723, .00821), (1.00217, .00357)):
            with self.subTest(a=a, gamma=g):
                np.testing.assert_allclose(patch.energy_kernel(a, g), source.energy_kernel(a, g),
                                           rtol=3e-9, atol=2e-10)

    def test_higher_degree_improves_smooth_derivative_error(self):
        source, coarse = make_patch(4)
        _, fine = make_patch(6)
        a, g = .99013, .00433
        exact = np.asarray(source.energy_kernel(a, g))
        coarse_error = np.max(abs(np.asarray(coarse.energy_kernel(a, g))-exact))
        fine_error = np.max(abs(np.asarray(fine.energy_kernel(a, g))-exact))
        self.assertGreater(coarse_error, 1e-10)
        self.assertLess(fine_error, coarse_error/20)

    def test_mixed_derivatives_obey_maxwell_identity_at_fixed_population(self):
        _, patch = make_patch(6)
        p = np.array([.4, .3, .12, .07, .008])
        a, g, h = .99411, .00623, 2e-6
        dgamma_ua = (patch.evaluate(a, g+h, p)[1]-patch.evaluate(a, g-h, p)[1])/(2*h)
        da_ug = (patch.evaluate(a+h, g, p)[2]-patch.evaluate(a-h, g, p)[2])/(2*h)
        self.assertAlmostEqual(dgamma_ua, da_ug, delta=5e-8)

    def test_sampling_stays_inside_exact_box_and_progress_is_complete(self):
        source = SmoothSource()
        progress = []
        patch = LocalEnergyPatch.build(source, (.985, 1.005), (.003, .009), 4,
            on_sample=lambda i, n: progress.append((i, n)))
        self.assertEqual(progress, [(i, 25) for i in range(1, 26)])
        self.assertEqual(len(source.samples), 25)
        np.testing.assert_array_equal(patch.count_nodes, source.count_nodes)
        np.testing.assert_array_equal(patch.count_weights, source.count_weights)
        self.assertEqual({a for a, _ in source.samples if a in (.985, 1.005)}, {.985, 1.005})
        self.assertEqual({g for _, g in source.samples if g in (.003, .009)}, {.003, .009})

    def test_even_one_ulp_outside_box_is_rejected_without_clipping(self):
        _, patch = make_patch()
        points = [(np.nextafter(.985, -np.inf), .006), (np.nextafter(1.005, np.inf), .006),
                  (.995, np.nextafter(.003, -np.inf)), (.995, np.nextafter(.009, np.inf))]
        for a, g in points:
            with self.subTest(a=a, gamma=g), self.assertRaises(ValueError):
                patch.energy_kernel(a, g)
            with self.subTest(zero_occupation=(a, g)), self.assertRaises(ValueError):
                patch.evaluate(a, g, np.zeros(5))

    def test_saved_file_is_not_overwritten_with_missing_npz_suffix(self):
        _, patch = make_patch()
        with tempfile.TemporaryDirectory(prefix='patch-extension-guard-') as folder:
            path = Path(folder)/'patch'
            try:
                patch.save(path)
            except ValueError:
                # Rejecting a non-.npz filename is also a safe documented API.
                self.assertFalse(path.exists())
                self.assertFalse(path.with_suffix('.npz').exists())
                return
            actual = path if path.exists() else path.with_suffix('.npz')
            before = actual.read_bytes()
            with self.assertRaises(FileExistsError):
                patch.save(path)
            self.assertEqual(actual.read_bytes(), before)

    def test_complex_coefficients_are_rejected_instead_of_discarding_imaginary_values(self):
        source, patch = make_patch()
        coefficients = patch.coefficients.astype(complex)
        coefficients[0, 0, 0] += 1e-3j
        with self.assertRaises(ValueError):
            LocalEnergyPatch(source, patch.amplitude_bounds, patch.gamma_bounds, coefficients)

    def test_complex_bounds_are_rejected_instead_of_discarding_imaginary_values(self):
        source, patch = make_patch()
        with self.assertRaises(ValueError):
            LocalEnergyPatch(source, np.array([.985+.01j, 1.005]), patch.gamma_bounds, patch.coefficients)

    def test_patch_coefficients_and_derivatives_are_immutable(self):
        _, patch = make_patch()
        for array in (patch.coefficients, patch._da, patch._dg):
            with self.assertRaises(ValueError):
                array.flat[0] = 1.

    def test_positive_but_nonmonotone_energy_is_rejected(self):
        source, patch = make_patch()
        coefficients = np.zeros_like(patch.coefficients)
        coefficients[0, 0] = [1., 2., 1.5, 3., 4.]
        broken = LocalEnergyPatch(source, patch.amplitude_bounds, patch.gamma_bounds, coefficients)
        with self.assertRaises(FloatingPointError):
            broken.energy_kernel(.995, .006)

    def test_load_rejects_changed_count_nodes_weights_or_physical_scales(self):
        _, patch = make_patch()
        with tempfile.TemporaryDirectory(prefix='patch-load-guard-') as folder:
            path = Path(folder)/'patch.npz'
            patch.save(path)
            for kind in ('count_nodes', 'count_weights', 'delta0_J', 'N0_per_J_m3', 'D_m2_s', 'eta'):
                other = SmoothSource()
                if kind in ('count_nodes', 'count_weights'):
                    setattr(other, kind, getattr(other, kind)*1.01)
                elif kind == 'eta':
                    other.eta *= 2
                else:
                    setattr(other.vacuum, kind, getattr(other.vacuum, kind)*1.01)
                with self.subTest(changed=kind), self.assertRaises(ValueError):
                    LocalEnergyPatch.load(path, other)


if __name__ == '__main__':
    unittest.main()
