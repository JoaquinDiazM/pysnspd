import unittest
from types import SimpleNamespace
import numpy as np
from pysnspd.mesh.pytdgl_like import PyTDGLLikeMeshParameters, generate_rectangular_pytdgl_fvm_mesh_from_parameters
from pysnspd.mesh.operators import FVOperators, divergence_from_edge_scalar
from pysnspd.gtdgl.tdgl_operators import build_laplacian
from pysnspd.experimental.thermal_mesh_adapter import graph_from_mesh, graph_from_fv_operators


class ThermalMeshAdapterTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        params = PyTDGLLikeMeshParameters(length_m=160e-9, width_m=80e-9,
            target_spacing_m=20e-9, max_edge_length_m=20e-9, seed=12345,
            boundary_points=21)
        cls.mesh = generate_rectangular_pytdgl_fvm_mesh_from_parameters(params)
        cls.nodes = cls.mesh.sites
        tolerance = 64*np.finfo(float).eps*params.length_m
        cls.fixed = np.flatnonzero((abs(cls.nodes[:, 0]) <= tolerance)
            | (abs(cls.nodes[:, 0]-params.length_m) <= tolerance))
        cls.ell0 = 2.5e-9

    def test_exact_geometry_and_coordinate_unit_invariance(self):
        graph = graph_from_mesh(self.mesh, self.ell0, fixed_nodes=self.fixed, origin_m=[1e-9, 0])
        np.testing.assert_array_equal(graph.edges, self.mesh.edge_mesh.edges)
        np.testing.assert_array_equal(graph.boundary_nodes, self.fixed)
        self.assertLess(len(graph.boundary_nodes), len(self.mesh.boundary_indices))
        scale = 3.7
        e = self.mesh.edge_mesh
        scaled = SimpleNamespace(sites=self.mesh.sites*scale, areas=self.mesh.areas*scale**2,
            edge_mesh=SimpleNamespace(edges=e.edges, edge_lengths=e.edge_lengths*scale,
                                     dual_edge_lengths=e.dual_edge_lengths*scale))
        other = graph_from_mesh(scaled, self.ell0*scale, fixed_nodes=self.fixed, origin_m=[scale*1e-9, 0])
        np.testing.assert_allclose(other.area_weights, graph.area_weights, rtol=5e-15)
        np.testing.assert_allclose(other.conductance, graph.conductance, rtol=5e-15)
        np.testing.assert_allclose(other.coordinates_bar, graph.coordinates_bar, atol=1e-16)

    def test_graph_dirichlet_energy_matches_production_laplacian_and_fv_divergence(self):
        graph = graph_from_mesh(self.mesh, self.ell0, fixed_nodes=self.fixed)
        e = self.mesh.edge_mesh; tail, head = graph.edges.T
        v = np.sin(np.arange(graph.n_nodes)*.4)+1j*np.cos(np.arange(graph.n_nodes)*.9)
        flux = graph.conductance*(v[head]-v[tail])
        numerator = np.zeros(graph.n_nodes, complex)
        np.add.at(numerator, tail, flux); np.add.at(numerator, head, -flux)
        lap, _ = build_laplacian(self.mesh)
        np.testing.assert_allclose(numerator/graph.area_weights, self.ell0**2*(lap@v), atol=2e-14)
        energy = np.sum(graph.conductance*abs(v[head]-v[tail])**2)
        self.assertAlmostEqual(energy, -np.vdot(v, numerator).real, delta=1e-13)
        vec = self.mesh.sites[head]-self.mesh.sites[tail]
        ops = FVOperators(e.edges, tail, head, vec, e.edge_lengths,
                          vec/e.edge_lengths[:, None], e.dual_edge_lengths, self.mesh.areas)
        other = graph_from_fv_operators(ops, self.mesh.sites, self.ell0, fixed_nodes=self.fixed)
        np.testing.assert_array_equal(other.area_weights, graph.area_weights)
        np.testing.assert_array_equal(other.conductance, graph.conductance)
        direct = divergence_from_edge_scalar((v.real[head]-v.real[tail])/e.edge_lengths, ops)
        np.testing.assert_allclose(numerator.real/graph.area_weights, self.ell0**2*direct, atol=2e-14)

    def test_no_implicit_boundary_or_geometry_repair(self):
        with self.assertRaises(TypeError):
            graph_from_mesh(self.mesh, self.ell0)
        with self.assertRaises(ValueError):
            graph_from_mesh(self.mesh, 0, fixed_nodes=self.fixed)
        with self.assertRaises(ValueError):
            graph_from_mesh(SimpleNamespace(areas=None, edge_mesh=None), self.ell0, fixed_nodes=self.fixed)


if __name__ == '__main__':
    unittest.main()
