"""Light diagnostics of stage4 follow-up metrics, with no causal spectra."""
import importlib.util
from pathlib import Path
import sys
import tempfile
import unittest

import numpy as np

from test_stage4_core_parameter import AnalyticCatalogue
from pysnspd.experimental.rectangular_spatial import RectangularSpatialFunctional
from pysnspd.experimental.cell_closures import CellScales, KWTMobility
from pysnspd.experimental.spatial_dynamics import kwt_spatial_response


ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location('stage4_review_controls', ROOT/'sandbox/stage4_core/controls.py')
controls = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = controls
SPEC.loader.exec_module(controls)


class TestStage4ReviewControls(unittest.TestCase):
    def model(self, nx=2, ny=1, degree=2, delta=.1):
        return RectangularSpatialFunctional(AnalyticCatalogue(), 160e-9, 80e-9, 7e-9,
            elements_x=nx, elements_y=ny, degree=degree, delta_regularizer_bar=delta)

    def test_exact_profile_derivatives_match_independent_coordinate_variations(self):
        coordinates = np.array([[0., -2.], [4., .5], [10., 0.], [16., 2.]])
        for profile in ('smooth', 'suppressed'):
            _, exact = controls.spatial_profile(coordinates, 20., profile)
            for axis in (0, 1):
                step = np.zeros_like(coordinates)
                step[:, axis] = 1e-5
                plus = controls.spatial_profile(coordinates+step, 20., profile)[0]
                minus = controls.spatial_profile(coordinates-step, 20., profile)[0]
                np.testing.assert_allclose((plus-minus)/2e-5, exact[:, axis], atol=2e-10, rtol=2e-8)

    def test_oriented_section_sum_uses_all_transverse_weights_and_reversed_edges(self):
        model = self.model()
        x = model.dof_coordinates_bar[:, 0]
        currents = np.zeros(len(model.graph_edges))
        longitudinal = x[model.graph_edges[:, 1]] > x[model.graph_edges[:, 0]]
        density = 2.3e7
        currents[longitudinal] = density*model.face_areas_m2[longitudinal]
        _, flux = controls.section_currents(model, currents)
        np.testing.assert_allclose(flux, density*model.cross_section_m2, rtol=1e-14)
        # The same physical flux with every edge orientation reversed.
        model.graph_edges = model.graph_edges[:, ::-1]
        _, reversed_flux = controls.section_currents(model, -currents)
        np.testing.assert_allclose(reversed_flux, flux, rtol=1e-14)

    def test_volume_norm_is_independent_of_nodal_measure_for_uniform_force(self):
        for nx, ny in ((2, 1), (4, 2)):
            model = self.model(nx, ny)
            gradient = model.mass_bar[:, None]*np.array([3., 4.])
            metrics, density = controls.volume_force_metrics(model, gradient, controls.boundary_mask(model))
            np.testing.assert_allclose(density, np.broadcast_to([3., 4.], density.shape))
            for scope in metrics.values():
                self.assertAlmostEqual(scope['rms_volume_bar'], 5., places=13)

    def test_fixed_conjugate_boundary_load_removes_boundary_response_with_zero_work(self):
        model = self.model()
        z = controls.spatial_field(model, 'suppressed')
        mask = controls.boundary_mask(model)
        self.assertEqual(int(mask.sum()), 2*sum(model.rectangle_shape)-4)
        p = np.zeros((model.cells, len(model.catalog.count_nodes)))
        result = model.evaluate(z, p, require_stability=False)
        load = np.zeros_like(result.gradient_cartesian_bar)
        load[mask] = result.gradient_cartesian_bar[mask]
        scales = CellScales(model.catalog.vacuum.delta0_J, model.catalog.vacuum.N0_per_J_m3, 8.65, .9, 1.)
        mob = KWTMobility(scales)
        args = (z, result.gradient_cartesian_bar, model.mass_bar, mob,
                np.zeros(model.cells), np.zeros(model.cells))
        free = kwt_spatial_response(*args)
        fixed = kwt_spatial_response(*args, boundary_load_bar=load)
        np.testing.assert_array_equal(fixed.material_velocity_bar[mask], np.zeros(mask.sum()))
        np.testing.assert_allclose(fixed.material_velocity_bar[~mask], free.material_velocity_bar[~mask], atol=0.)
        self.assertEqual(fixed.boundary_work_rate_bar, 0.)
        self.assertEqual(float(np.sum(fixed.heat_density_bar[mask])), 0.)
        self.assertAlmostEqual(fixed.field_energy_rate_bar+fixed.condensate_heat_rate_bar, 0., delta=1e-11)

    def test_analytic_center_and_discrete_center_are_explicitly_different(self):
        errors = []
        for nx, ny in ((4, 2), (8, 4)):
            model = self.model(nx, ny, degree=4, delta=.05)
            fields = model.sample_fields(controls.spatial_field(model, 'suppressed'))
            report, _, _ = controls.derivative_diagnostics(model, 'suppressed', fields, controls.boundary_mask(model))
            center = report['center']
            self.assertAlmostEqual(center['gamma_analytic'], .0004252166, delta=1e-10)
            errors.append(center['derivative_error_norm'])
        self.assertLess(errors[1], errors[0]/2)

    def test_spatial_output_keeps_history_and_explicit_fixed_control(self):
        plan = dict(material=dict(Tc_K=8.65, Tb_K=.9, width_m=80e-9, thickness_m=7e-9),
            geometry=dict(length_m=160e-9), mobility_pairs_ps=dict(memory=[.5, 2.47]))
        case = dict(profile='suppressed', delta_reg=.1, elements_x=2, elements_y=1,
            degree=2, population='vacuum', finite_differences=False,
            boundary_response='fixed_prescribed_all_sides')
        with tempfile.TemporaryDirectory() as path:
            result = controls.spatial_case(plan, case, AnalyticCatalogue(), lambda *args: None, Path(path))
            fixed = result['mobility_fixed_boundary']['memory']
            self.assertEqual(fixed['maximum_boundary_velocity_per_ps'], 0.)
            self.assertEqual(fixed['boundary_work_bar'], 0.)
            self.assertGreater(result['mobility']['memory']['boundary_heat_bar'], 0.)
            self.assertEqual(result['response_boundary_policy'], 'fixed_prescribed_all_sides')
            with np.load(Path(path)/'fields.npz') as fields:
                self.assertIn('section_current_A', fields)
                self.assertIn('derivative_analytic', fields)
                self.assertEqual(fields['boundary_mask'].dtype, np.dtype(bool))
        old_case = {key: value for key, value in case.items() if key != 'boundary_response'}
        with tempfile.TemporaryDirectory() as path:
            historical = controls.spatial_case(plan, old_case, AnalyticCatalogue(), lambda *args: None, Path(path))
            self.assertEqual(historical['response_boundary_policy'], 'unconstrained_instantaneous')
            self.assertEqual(historical['mobility_fixed_boundary'], {})
            self.assertEqual(historical['energy_bar'], result['energy_bar'])
            self.assertEqual(historical['mobility'], result['mobility'])


if __name__ == '__main__':
    unittest.main()
