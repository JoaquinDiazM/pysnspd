from pathlib import Path
import unittest
import numpy as np
from pysnspd.experimental import thermal_spatial_usadel as thermal
from pysnspd.experimental import thermal_stable_newton as stable
from pysnspd.experimental.thermal_snapshot import spectral_difference


class StableThermalNewtonTests(unittest.TestCase):
    def test_same_stationary_solution_and_contacts_as_original(self):
        graph=thermal.rectangular_graph(7,6,3.,2.);xy=graph.coordinates_bar
        d=(1+.1*np.sin(xy[:,0]))*np.exp(.2j*xy[:,1]);eps=.9;fixed=graph.boundary_nodes
        options=dict(fixed_nodes=fixed,fixed_u=d[fixed]/eps,tol=1e-7)
        old=thermal.solve_frequency(graph,d,eps,**options);new=stable.solve_frequency(graph,d,eps,**options)
        np.testing.assert_allclose(new.f,old.f,atol=2e-7,rtol=0)
        self.assertLessEqual(new.residual,1e-7)
        np.testing.assert_array_equal(new.u[fixed],d[fixed]/eps)

    def test_real_roundoff_regression_keeps_original_tolerance(self):
        root=Path(__file__).resolve().parents[1]
        with np.load(root/'docs/implementation/stage4/moment_review_20260924/thermal_weak/raw/full_node_operator_checks.npz') as a:
            graph=thermal.ThermalGraph(a['area_weights'],a['edges'],a['conductance'],a['coordinates_bar'],a['boundary_nodes'])
        with np.load(root/'docs/implementation/stage4/time_review_20260924/nonlinear_snapshots/roundoff_fixture.npz') as a:
            d,u,eps=a['d'],a['u'],float(a['epsilon'])
        scale=graph.area_weights*np.maximum(1.,abs(d));free=np.ones(len(d),bool);free[graph.boundary_nodes]=False
        before=thermal.spectral_energy_gradient(graph,d,eps,u)
        self.assertGreater(np.max(abs(before.residual[free])/scale[free]),1e-7)
        result=stable.solve_frequency(graph,d,eps,initial_u=u,fixed_nodes=graph.boundary_nodes,fixed_u=u[graph.boundary_nodes],tol=1e-7)
        self.assertLessEqual(result.residual,1e-7)
        self.assertLessEqual(result.iterations,2)
        self.assertLess(spectral_difference(graph,d,u,d,result.u,eps)['renormalized_energy_difference'],0.)


if __name__=='__main__':unittest.main()
