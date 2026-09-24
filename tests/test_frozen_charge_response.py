"""Controls for causal resolvent and two independent phasor quadratures."""
import unittest
import numpy as np
from pysnspd.experimental import thermal_spatial_usadel as thermal
from pysnspd.experimental import retarded_spatial_usadel as retarded
from pysnspd.experimental import frozen_kinetic_usadel as kinetic
from pysnspd.experimental import frozen_charge_response as response


class ChargeFrequencyTests(unittest.TestCase):
    def setUp(self):
        self.graph=thermal.rectangular_graph(5,5,2.,2.)
        xy=self.graph.coordinates_bar
        d=1.2*np.exp(.4j*xy[:,0])
        s=retarded.solve(self.graph,d,.05-.8j,tolerance=1e-11)
        self.op=kinetic.assemble(self.graph,d,s.g,s.f,s.f_tilde)
        self.hL=np.maximum(0.,1-np.sum((xy-[1.,0.])**2,axis=1))**3

    def test_zero_frequency_is_static_elimination(self):
        actual=response.solve(self.op,self.hL,self.graph.boundary_nodes,0.)
        expected=self.op.charge_response(self.hL,self.graph.boundary_nodes)
        np.testing.assert_allclose(actual.directions,expected,atol=1e-14)

    def test_resolvent_storage_and_contacts(self):
        actual=response.solve(self.op,self.hL,self.graph.boundary_nodes,.1)
        obs=response.observe(self.op,actual.directions)
        dynamic=obs['charge_residual']+1j*.1*actual.storage*actual.directions[:,1]
        self.assertLess(np.max(abs(dynamic[actual.free_nodes])),1e-13)
        np.testing.assert_array_equal(actual.directions[self.graph.boundary_nodes],0.)

    def test_phasor_linearity_retains_four_force_coordinates(self):
        state=response.solve(self.op,self.hL,self.graph.boundary_nodes,.1).directions
        first=response.observe(self.op,state)
        second=response.observe(self.op,1j*state)
        for key in first:
            np.testing.assert_allclose(second[key],1j*first[key],atol=1e-14)
        self.assertEqual(first['force_cartesian'].shape,(self.graph.n_nodes,2))

    def test_dynamic_ward_identity_in_each_time_quadrature(self):
        solution=response.solve(self.op,self.hL,self.graph.boundary_nodes,.1)
        obs=response.observe(self.op,solution.directions)
        divergence=np.zeros(self.graph.n_nodes,complex)
        np.add.at(divergence,self.graph.edges[:,0],obs['current'])
        np.add.at(divergence,self.graph.edges[:,1],-obs['current'])
        torque=self.op.d.real*obs['force_cartesian'][:,1]-self.op.d.imag*obs['force_cartesian'][:,0]
        np.testing.assert_allclose(divergence+torque,2*obs['charge_residual'],atol=1e-14)

    def test_normal_limit_physical_diffusion_mass(self):
        n=self.graph.n_nodes
        op=kinetic.assemble(self.graph,np.zeros(n),np.ones(n),np.zeros(n),np.zeros(n))
        result=response.solve(op,self.hL,self.graph.boundary_nodes,.1)
        np.testing.assert_array_equal(result.storage,self.graph.area_weights)
        np.testing.assert_array_equal(result.directions[:,1],np.zeros(n))
        # With ell0^2=hbar D/(2 kBTc), storage time tD=ell0^2/D.
        from scipy.constants import hbar,k
        Tc,D=8.65,5e-5
        self.assertAlmostEqual((hbar*D/(2*k*Tc))/D,hbar/(2*k*Tc),places=25)


if __name__=='__main__':unittest.main()
