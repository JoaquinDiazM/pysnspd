"""Tiny full-graph end-to-end controls of the biased energy worker."""
import tempfile
import json
import multiprocessing as mp
from pathlib import Path
import time
import unittest
import numpy as np
from scipy.linalg import eigh
from sandbox.stage4_core.biased_harmonic_query import energy_query
from sandbox.stage4_core.biased_coupled_response import thermal_query, fail_fast_pool
from sandbox.stage4_core import coupled_response as uniform_campaign
from pysnspd.experimental import thermal_spatial_usadel as thermal
from pysnspd.experimental import retarded_spatial_usadel as spectral


def make_reference(path, phase_bias=.2, gap=1.2):
    base = thermal.rectangular_graph(5, 3, 2., 1.)
    x = base.coordinates_bar[:, 0]
    fixed = np.flatnonzero((x == 0.) | (x == 2.))
    graph = thermal.ThermalGraph(base.area_weights, base.edges, base.conductance,
        base.coordinates_bar, fixed)
    d = gap*np.exp(1j*phase_bias*(x/2-.5))
    epsilon = .5
    fa, fb = spectral.uniform_fields(d[fixed], epsilon)
    result = spectral.solve(graph, d, epsilon, fixed_nodes=fixed, fixed_a=fa,
        fixed_b=fb, tolerance=1e-12)
    L = np.zeros((graph.n_nodes, graph.n_nodes))
    for (i, j), c in zip(graph.edges, graph.conductance):
        L[i, i] += c; L[j, j] += c; L[i, j] -= c; L[j, i] -= c
    free = np.ones(graph.n_nodes, bool); free[fixed] = False
    values, vectors = eigh(L[free][:, free], np.diag(graph.area_weights[free]))
    modes = np.zeros((graph.n_nodes, 2)); modes[free] = vectors[:, :2]
    np.savez_compressed(path, delta_bar=d, epsilon_bar=np.array([epsilon]),
        f=result.f[None], g=result.g[None], u=(result.f/result.g)[None],
        area_weights=graph.area_weights, edges=graph.edges,
        conductance=graph.conductance, coordinates_bar=graph.coordinates_bar,
        boundary_nodes=fixed, alpha=np.zeros(len(graph.edges)), modes=modes,
        lift=.5-x/2)
    return graph


class BiasedHarmonicQueryTests(unittest.TestCase):
    def test_biased_worker_connects_spectrum_kinetics_and_all_forcings(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory)/'state.npz'
            graph = make_reference(path)
            result = energy_query(dict(reference_npz=str(path), energy=.7,
                omega=.2, eta=.03, t=.3, spectral_tolerance=1e-10,
                continuation_steps=6))
            thermal_force, thermal_current, thermal_residual = thermal_query(dict(
                reference_npz=str(path), modes_npz=str(path), n=0, t=.3))
            self.assertEqual(thermal_force.shape, (graph.n_nodes, 2, 5))
            self.assertEqual(thermal_current.shape, (len(graph.edges), 5))
            self.assertLess(thermal_residual, 1e-10)
        self.assertEqual(result['anomaly_correction'].shape, (graph.n_nodes, 2, 7))
        self.assertEqual(result['current_correction'].shape, (len(graph.edges), 7))
        self.assertEqual(result['charge_density'].shape, (graph.n_nodes, 7))
        self.assertTrue(np.isfinite(result['anomaly_correction']).all())
        self.assertTrue(np.isfinite(result['distributions']).all())
        np.testing.assert_allclose(result['anomaly_density_correction']*graph.area_weights[:, None, None],
            result['anomaly_correction'], atol=1e-14)
        np.testing.assert_array_equal(result['charge'], result['charge_density'])
        np.testing.assert_array_equal(result['anomaly_static'][:, :, 4:6], 0.)
        np.testing.assert_array_equal(result['current_static'][:, 4:6], 0.)
        self.assertGreater(np.max(abs(result['distributions'][:, 1, :2])), 1e-7)
        for row in result['metadata']['response_metrics']:
            self.assertLess(row['kinetic_scaled_residual'], 1e-8)
            self.assertLess(row['full_spectral_scaled_max'], 1e-7)
        for history in result['metadata']['continuation'].values():
            self.assertEqual(len(history), 6)
            self.assertGreaterEqual(min(row['minimum_DOS'] for row in history), -1e-10)

    def test_failed_query_cancels_queued_campaign_instead_of_waiting_for_sleep(self):
        with tempfile.TemporaryDirectory() as directory:
            manifest = {'status': 'RUNNING', 'active_query': {'phase': 'test', 'index': 1}}
            started = time.monotonic()
            with self.assertRaises(ValueError):
                with fail_fast_pool(directory, manifest, max_workers=2,
                        mp_context=mp.get_context('spawn')) as pool:
                    pool.submit(time.sleep, 20.)
                    failed = pool.submit(int, 'deliberate-worker-failure')
                    for _ in range(6):
                        pool.submit(time.sleep, 20.)
                    failed.result(timeout=10.)
            elapsed = time.monotonic()-started
            receipt = json.loads((Path(directory)/'failure.json').read_text())
            self.assertEqual(receipt['status'], 'FAILED')
            self.assertEqual(receipt['exception'], 'ValueError')
            self.assertTrue(receipt['terminated_owned_worker_pids'])
            self.assertLess(elapsed, 10.)
            self.assertEqual(manifest['status'], 'FAILED')

    def test_uniform_worker_port_matches_exact_gauge_and_thermal_current(self):
        gap, omega, temperature = 1.76392592757718, .2, .3
        eta = .01*gap
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory)/'uniform.npz'
            graph = make_reference(path, phase_bias=0., gap=gap)
            lift = .5-graph.coordinates_bar[:, 0]/2
            tail, head = graph.edges.T
            gradient = graph.conductance*(lift[head]-lift[tail])
            _, thermal_current, _ = thermal_query(dict(reference_npz=str(path),
                modes_npz=str(path), n=0, t=temperature))
            f0 = gap/np.hypot(.5, gap)
            np.testing.assert_allclose(thermal_current[:, -1],
                4*np.pi*temperature*f0*f0*gradient, atol=1e-12)
            for energy in (.5, gap-omega/2+.1*eta, gap+omega/2+.1*eta):
                result = energy_query(dict(reference_npz=str(path), energy=energy,
                    omega=omega, eta=eta, t=temperature, spectral_tolerance=1e-10,
                    continuation_steps=6))
                Rp, Ap, _, _ = uniform_campaign.backgrounds(gap, np.array([energy+omega/2]), eta)
                Rm, Am, _, _ = uniform_campaign.backgrounds(gap, np.array([energy-omega/2]), eta)
                Rp, Ap, Rm, Am = Rp[0], Ap[0], Rm[0], Am[0]
                hp, hm = np.tanh((energy+omega/2)/(2*temperature)), np.tanh((energy-omega/2)/(2*temperature))
                Z = uniform_campaign.Z
                Kp, Km = hp*(Rp-Ap), hm*(Rm-Am)
                dr, da = (Z@Rm-Rp@Z)/omega, (Z@Am-Ap@Z)/omega
                dk = (Z@Km-Kp@Z)/omega
                flux = .5*(Rp@dk-dr@Km+Kp@da-dk@Am)
                expected_current = 2*gradient*uniform_campaign.projection(flux)[1]
                expected_charge = lift*uniform_campaign.projection(dk)[0]
                expected_hT = -lift*(hp-hm)/omega
                np.testing.assert_allclose(result['current_dynamic'][:, -1], expected_current, atol=2e-8, rtol=2e-8)
                np.testing.assert_allclose(result['charge'][:, -1], expected_charge, atol=2e-8, rtol=2e-8)
                np.testing.assert_allclose(result['distributions'][:, 1, -1], expected_hT, atol=2e-8, rtol=2e-8)

    def test_signed_high_energy_moments_follow_cubic_tail_without_noise_growth(self):
        gap = 1.76392592757718
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory)/'tail_reference.npz'
            graph = make_reference(path, phase_bias=.4, gap=gap)
            with np.load(path) as archive:
                modes = archive['modes'].copy()
            def signed_pair(energy, tolerance):
                values = []
                for sign in (-1., 1.):
                    result = energy_query(dict(reference_npz=str(path), energy=sign*energy,
                        omega=6., eta=.01*gap, t=.9/8.65,
                        spectral_tolerance=tolerance, continuation_steps=12))
                    values.append(np.vstack((modes.T@result['anomaly_correction'][:, 0],
                        modes.T@result['anomaly_correction'][:, 1],
                        modes.T@(graph.area_weights[:, None]*result['charge']))))
                return values[0]+values[1]
            first, second = signed_pair(128., 1e-8), signed_pair(256., 1e-8)
            relative = np.linalg.norm(first*128**3-second*256**3)/np.linalg.norm(second*256**3)
            self.assertLess(relative, .01)
            high, tight = signed_pair(12000., 1e-8), signed_pair(12000., 1e-11)
            self.assertTrue(np.isfinite(high).all())
            self.assertLess(np.linalg.norm(high-tight)/np.linalg.norm(tight), .02)
            # The semi-infinite change E=C/x multiplies the integrand by E^2/C.
            self.assertLess(np.max(abs(high-tight))*12000**2/32, 1e-6)


if __name__ == '__main__':
    unittest.main()
