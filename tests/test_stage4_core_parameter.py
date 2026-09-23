"""Changed-core identities on analytic data, without physical admission.

These tests expose stale regularizer constants in assembled and thermal paths.
The catalogue is synthetic; it does not extend the physical near-zero support.
"""
from types import SimpleNamespace
import unittest

import numpy as np

from pysnspd.experimental.energy_catalog import BCS_GAP_RATIO, HBAR_J_S, K_B_J_K
from pysnspd.experimental.spatial_functional import PeriodicSpatialFunctional
from pysnspd.experimental.spatial_nodal import PeriodicNodalSpatialFunctional
from pysnspd.experimental.spatial_open import OpenSpatialFunctional
from pysnspd.experimental.mixed_spatial import MixedSpatialFunctional


class AnalyticCatalogue:
    def __init__(self):
        self.count_nodes = np.array([.1, .5, 1.2])
        self.count_weights = np.array([.2, .3, .5])
        self.vacuum = SimpleNamespace(metadata={"Tc_K": 8.65}, D_m2_s=5e-5,
            N0_per_J_m3=2.1e47, delta0_J=BCS_GAP_RATIO*K_B_J_K*8.65,
            gamma_axis=np.array([0., 2.]), evaluate=self.vacuum_values)

    @staticmethod
    def vacuum_values(a, gamma):
        return (.5*(a*a-1)**2+.3*gamma*a*a+.07*gamma**2,
                2*a*(a*a-1)+.6*gamma*a, .3*a*a+.14*gamma)

    def energy_kernel(self, a, gamma):
        root = np.sqrt(self.count_nodes**2+a*a)
        alpha = .13+.03*self.count_nodes
        return (root+alpha*gamma*(1+.2*a*a), a/root+.4*alpha*gamma*a,
                alpha*(1+.2*a*a))


def gradient(result):
    return result.gradient_cartesian_bar[:, 0]+1j*result.gradient_cartesian_bar[:, 1]


def fd5(function, h=3e-5):
    return (function(-2*h)-8*function(-h)+8*function(h)-function(2*h))/(12*h)


class TestStage4CoreParameter(unittest.TestCase):
    schemes = ("midpoint", "nodal", "open", "mixed")

    def setUp(self):
        self.catalog = AnalyticCatalogue()
        self.ell = np.sqrt(HBAR_J_S*self.catalog.vacuum.D_m2_s/(2*K_B_J_K*8.65))

    def model(self, scheme, **kwargs):
        if scheme == "mixed":
            return MixedSpatialFunctional(self.catalog, 8*self.ell, 4*self.ell, 7e-9,
                3*self.ell, 3*self.ell, elements_x=1, elements_y=1, left_elements=1,
                right_elements=1, degree=2, **kwargs)
        if scheme == "open":
            return OpenSpatialFunctional(self.catalog, 12*self.ell, 80e-9*7e-9,
                                         2, degree=2, **kwargs)
        cls = PeriodicSpatialFunctional if scheme == "midpoint" else PeriodicNodalSpatialFunctional
        return cls(self.catalog, 12*self.ell, 80e-9*7e-9, 5, **kwargs)

    def state(self, scheme, model):
        if scheme == "mixed":
            x, y = model.dof_coordinates_bar.T
            envelope = np.where((x > 0)&(x < 8), x*(8-x)/16, 0.)
            z = (.22+.012*np.cos(x/3)+.008*envelope*y)*np.exp(1j*(.06*x+.012*envelope*y))
            count, edges = model.quadrature_size, len(model.graph_edges)
        else:
            x = model.x_bar if scheme == "open" else np.arange(model.cells)*model.h_bar
            z = (.22+.02*np.cos(2*np.pi*x/12))*np.exp(1j*(.06*x+.02*np.sin(x)))
            count, edges = model.cells, model.cells-(scheme == "open")
        p = np.array([(.04+.001*i)*np.exp(-self.catalog.count_nodes) for i in range(count)])
        links = .007*np.sin(np.arange(edges)+.2)
        return z, p, links

    def test_default_and_explicit_point_one_are_identical(self):
        for scheme in self.schemes:
            with self.subTest(scheme=scheme):
                default, explicit = self.model(scheme), self.model(scheme, delta_regularizer_bar=.1)
                self.assertEqual(default.delta_regularizer_bar, .1)
                z, p, links = self.state(scheme, default)
                first = default.evaluate(z, p, link_phases=links, require_stability=False)
                second = explicit.evaluate(z, p, link_phases=links, require_stability=False)
                self.assertEqual(first.energy_bar, second.energy_bar)
                np.testing.assert_array_equal(first.gradient_cartesian_bar, second.gradient_cartesian_bar)
                np.testing.assert_array_equal(first.link_current_bar, second.link_current_bar)

    def test_nonfinite_nonpositive_and_nonscalar_regularizers_rejected(self):
        for scheme in self.schemes:
            for delta in (0., -.1, np.nan, np.inf, None, True, .1+0j, [.1], ".1"):
                with self.subTest(scheme=scheme, delta=delta), self.assertRaisesRegex(ValueError, "delta_regularizer_bar"):
                    self.model(scheme, delta_regularizer_bar=delta)

    def test_changed_regularizer_force_and_current_are_energy_derivatives(self):
        for scheme in self.schemes:
            results = []
            for delta in (.05, .2):
                with self.subTest(scheme=scheme, delta=delta):
                    model = self.model(scheme, delta_regularizer_bar=delta)
                    z, p, links = self.state(scheme, model)
                    original = p.copy()
                    result = model.evaluate(z, p, link_phases=links, require_stability=False)
                    dz = .13*np.cos(np.arange(model.cells))+.1j*np.sin(np.arange(model.cells)+.2)
                    dl = .12*np.cos(np.arange(len(links))+.4)
                    evaluate = lambda t, zz, ll: model.evaluate(z+t*zz, p,
                        link_phases=links+t*ll, require_stability=False).energy_bar
                    self.assertAlmostEqual(fd5(lambda t: evaluate(t, dz, 0)),
                        np.real(np.vdot(gradient(result), dz)), delta=3e-8)
                    self.assertAlmostEqual(fd5(lambda t: evaluate(t, 0, dl)),
                        np.dot(result.link_current_bar, dl), delta=3e-8)
                    np.testing.assert_array_equal(p, original)
                    results.append(result)
            self.assertGreater(abs(results[0].energy_bar-results[1].energy_bar), 1e-6)
            self.assertGreater(np.linalg.norm(results[0].link_current_bar-results[1].link_current_bar), 1e-5)

    def test_changed_thermal_paths_preserve_free_energy_and_envelope_identities(self):
        for scheme in ("open", "mixed"):
            for delta in (.05, .2):
                with self.subTest(scheme=scheme, delta=delta):
                    model = self.model(scheme, delta_regularizer_bar=delta)
                    z, _, links = self.state(scheme, model)
                    result = model.evaluate_thermal(z, .3, link_phases=links, require_stability=False)
                    dz = .13*np.cos(np.arange(model.cells))+.1j*np.sin(np.arange(model.cells)+.2)
                    dl = .12*np.cos(np.arange(len(links))+.4)
                    measured = fd5(lambda t: model.evaluate_thermal(z+t*dz, .3,
                        link_phases=links+t*dl, require_stability=False).free_energy_bar)
                    expected = np.real(np.vdot(gradient(result), dz))+np.dot(result.link_current_bar, dl)
                    self.assertAlmostEqual(measured, expected, delta=3e-8)
                    p = result.p_quadrature if scheme == "mixed" else result.p_nodes
                    frozen = model.evaluate(z, p, link_phases=links, require_stability=False)
                    self.assertAlmostEqual(frozen.energy_bar, result.energy_bar, delta=1e-13)
                    np.testing.assert_allclose(frozen.gradient_cartesian_bar, result.gradient_cartesian_bar, atol=1e-13)
                    np.testing.assert_allclose(frozen.link_current_bar, result.link_current_bar, atol=1e-13)

    def test_D36_is_the_changed_core_local_gradient_hessian(self):
        # Differentiate the independently assembled scalar density with fixed z,p.
        z, p = .18+.11j, np.array([.04, .02, .01])
        for dimensions, scheme in ((1, "midpoint"), (2, "mixed")):
            for delta in (.05, .2):
                with self.subTest(dimensions=dimensions, delta=delta):
                    model = self.model(scheme, delta_regularizer_bar=delta)
                    d = np.array([.017+.029j, -.012+.018j])[:dimensions]
                    symbol = (model.principal_symbol(z, d[0], p) if dimensions == 1
                              else model.principal_symbol(z, d, p))
                    def density(vector):
                        derivative = vector[::2]+1j*vector[1::2]
                        q = np.imag(np.conj(z)*derivative)/(abs(z)**2+delta**2)
                        gamma = np.dot(q, q)/BCS_GAP_RATIO
                        energy = self.catalog.vacuum_values(abs(z), gamma)[0]
                        energy += 4*np.dot(self.catalog.count_weights*p,
                                            self.catalog.energy_kernel(abs(z), gamma)[0])
                        return energy+np.pi/4*(np.vdot(derivative, derivative).real-abs(z)**2*np.dot(q, q))
                    point = np.column_stack((d.real, d.imag)).ravel()
                    h, size = 2e-5, 2*dimensions
                    measured = np.empty((size, size))
                    for i in range(size):
                        ei = np.eye(size)[i]*h
                        for j in range(size):
                            ej = np.eye(size)[j]*h
                            measured[i, j] = (density(point+ei+ej)-density(point+ei-ej)
                                -density(point-ei+ej)+density(point-ei-ej))/(4*h*h)
                    np.testing.assert_allclose(symbol.matrix, measured, atol=1e-6, rtol=1e-6)


if __name__ == "__main__":
    unittest.main()
