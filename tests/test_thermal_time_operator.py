import unittest
import numpy as np
from pysnspd.experimental import thermal_spatial_usadel as thermal
from pysnspd.experimental.thermal_weak_response import SpectralTangent
from pysnspd.experimental.thermal_time_operator import SavedSpectralMode,FullNodeCoordinates


class ThermalTimeOperatorTests(unittest.TestCase):
    def test_saved_factor_matches_live_full_node_tangent(self):
        graph=thermal.rectangular_graph(6,5,3.,2.);xy=graph.coordinates_bar
        d=(1+.05*xy[:,0])*np.exp(.1j*xy[:,1]);t=.35;epsilon=np.pi*t
        s=thermal.solve_frequency(graph,d,epsilon,fixed_nodes=graph.boundary_nodes,fixed_u=d[graph.boundary_nodes]/epsilon,tol=1e-11)
        live=SpectralTangent(graph,d,s);J=live.jacobian
        arrays=dict(u=s.u,f=s.f,g=s.g,free_real_components=live.components,jacobian_data=J.data,
            jacobian_indices=J.indices,jacobian_indptr=J.indptr,jacobian_shape=J.shape)
        rng=np.random.default_rng(23);v=.1*(rng.normal(size=len(d))+1j*rng.normal(size=len(d)));v[graph.boundary_nodes]=0
        actual=SavedSpectralMode(graph,arrays,epsilon,t).apply(v);expected=live.apply(v)
        np.testing.assert_allclose(actual[0],4*np.pi*t*graph.area_weights*(v/epsilon-expected.df),atol=1e-14)
        np.testing.assert_allclose(actual[1],2*np.pi*t*expected.current_derivative,atol=1e-14)
        self.assertLess(actual[2],1e-13)
        coord=FullNodeCoordinates(graph,1.76)
        np.testing.assert_allclose(coord.decode(coord.encode(v)),v,atol=1e-15)
        self.assertAlmostEqual(np.linalg.norm(coord.encode(v))**2,np.sum(graph.area_weights*abs(v/1.76)**2))


if __name__=='__main__':unittest.main()
