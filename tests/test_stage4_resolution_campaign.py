"""Cheap controls for registration and physically consistent mesh comparison."""
import json
from pathlib import Path
import sys
import unittest
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'sandbox/stage4_core'))
import resolution_campaign as campaign
from resolution_campaign_analysis import transfer_flux, rectangular_interpolate
from pysnspd.experimental.thermal_spatial_usadel import rectangular_graph


class ResolutionCampaignTests(unittest.TestCase):
    def test_nested_registration_reuses_each_spectrum_once(self):
        plan=json.loads(campaign.DEFAULT_PLAN.read_text())
        jobs=campaign.spectral_jobs(plan)
        self.assertEqual(len(jobs),181)
        self.assertEqual(len({job['id'] for job in jobs}),181)
        self.assertEqual(len(campaign.spectral_jobs(plan,True)),1)
        self.assertTrue(set(plan['energy_grids']['base31'])<=set(plan['energy_grids']['fine50']))

    def test_thermal_susceptibility_quadrature_is_reported_without_renormalizing(self):
        plan=json.loads(campaign.DEFAULT_PLAN.read_text())
        t=plan['spectral']['T_K']/plan['spectral']['Tc_K']
        estimates=[]
        for key in ('base31','fine50'):
            energy=plan['spectral']['gap_reference_kBTc']*np.array(plan['energy_grids'][key])
            q=np.exp(-energy/t);chi=2*q/(t*(1+q)**2)
            weights=np.zeros(len(energy));weights[:-1]+=np.diff(energy)/2;weights[1:]+=np.diff(energy)/2
            estimates.append(float(weights@chi))
        self.assertGreater(estimates[0],estimates[1])
        self.assertAlmostEqual(estimates[0],1.016162,places=5)
        self.assertAlmostEqual(estimates[1],1.004243,places=5)

    def test_mesh_comparison_preserves_linear_physical_flux_density(self):
        source,target=rectangular_graph(9,9,2.,2.),rectangular_graph(5,5,2.,2.)
        def fields(graph):
            return {key:getattr(graph,key) for key in ('coordinates_bar','edges','conductance')}
        def exact(graph):
            xy,edges=graph.coordinates_bar,graph.edges
            vector=xy[edges[:,1]]-xy[edges[:,0]]
            midpoint=(xy[edges[:,0]]+xy[edges[:,1]])/2
            component=np.where(abs(vector[:,0])>0,2+midpoint[:,0],3-midpoint[:,1])
            return component*graph.conductance*np.linalg.norm(vector,axis=1)
        np.testing.assert_allclose(transfer_flux(fields(source),fields(target),exact(source)),exact(target),atol=1e-14)

    def test_node_density_interpolation_preserves_linear_complex_field(self):
        source,target=rectangular_graph(9,9,2.,2.),rectangular_graph(5,5,2.,2.)
        value=lambda xy:2+xy[:,0]+1j*(3-xy[:,1])
        np.testing.assert_allclose(rectangular_interpolate(source.coordinates_bar,value(source.coordinates_bar),target.coordinates_bar),value(target.coordinates_bar),atol=1e-14)


if __name__=='__main__':
    unittest.main()
