"""Fresh effective DC preparation; no historical state relabeling."""
import tempfile
from pathlib import Path
import unittest
import numpy as np
from scipy.constants import Boltzmann as kB,hbar
from pysnspd.experimental.thermal_spatial_usadel import rectangular_graph
from pysnspd.experimental.effective_bcs_functional import EffectiveBCSFunctional,bcs_count_quadrature
from pysnspd.experimental.effective_bcs_preparation import (
    load_dual_geometry,get_effective_geometry,prepare_effective_dc,save_effective_dc)


class TestEffectiveBCSPreparation(unittest.TestCase):
    def model(self,temperature=.9):
        eu=kB*8.65;ell=np.sqrt(hbar*5e-5/(2*eu))
        graph=rectangular_graph(7,5,30e-9/ell,80e-9/ell)
        xi,w=bcs_count_quadrature(32.,6)
        model=EffectiveBCSFunctional(graph,energy_scale_J=eu,
            N0_per_J_m3=9.153297821e46,ell0_m=ell,thickness_m=7e-9,
            count_nodes=xi,count_weights=w,gradient_law="bcs_radial",
            bath_temperature_K=temperature,diffusion_m2_s=5e-5)
        x=graph.coordinates_bar[:,0]
        return model,np.flatnonzero(x==x.min()),np.flatnonzero(x==x.max())

    def test_new_dc_gap_current_and_nv_baseline(self):
        model,left,right=self.model()
        dc=prepare_effective_dc(model,left_nodes=left,right_nodes=right,
            temperature_K=.9,current_A=21.297206e-6,width_m=80e-9)
        self.assertFalse(dc.receipt["old_DC_state_reused"])
        self.assertLess(dc.receipt["scaled_gap_force_max"],2e-10)
        self.assertLess(dc.receipt["current_relative_defect"],2e-10)
        self.assertLess(np.max(abs(dc.phi_V)),1e-9)
        np.testing.assert_array_equal(dc.p,model.fermi_dirac(dc.delta_bar,.9))
        np.testing.assert_allclose(dc.delta_bar[right],dc.bulk.amplitude_bar,rtol=0,atol=0)

    def test_thermal_hessian_response_at_populated_gap(self):
        model,left,right=self.model(6.)
        dc=prepare_effective_dc(model,left_nodes=left,right_nodes=right,
            temperature_K=6.,current_A=5e-6,width_m=80e-9)
        self.assertGreater(dc.p.max(),.05)
        self.assertLess(dc.receipt["scaled_gap_force_max"],2e-10)
        self.assertLess(dc.receipt["current_relative_defect"],2e-10)
        self.assertLess(dc.receipt["potential_max_abs_V"],1e-9)

    def test_inherited_geometry_arrays_and_hash_are_preserved(self):
        model,left,right=self.model()
        graph=model.graph;xy=graph.coordinates_bar*model.ell0_m
        t,h=graph.edges.T
        lengths=np.linalg.norm(xy[h]-xy[t],axis=1)
        with tempfile.TemporaryDirectory() as directory:
            path=Path(directory)/"fixture.npz"
            np.savez(path,coordinates_m=xy,coordinates_bar=graph.coordinates_bar,
                triangles=np.array([[0,1,5]],int),edges=graph.edges,
                area_weights=graph.area_weights,conductance=graph.conductance,
                edge_length_m=lengths,dual_face_length_m=lengths*graph.conductance,
                fixed_nodes=np.r_[left,right],left_nodes=left,right_nodes=right,
                ell0_m=np.array(model.ell0_m))
            loaded=load_dual_geometry(path)
            np.testing.assert_array_equal(loaded.graph.edges,graph.edges)
            np.testing.assert_array_equal(loaded.graph.area_weights,graph.area_weights)
            np.testing.assert_allclose(loaded.face_areas_m2(7e-9),
                lengths*graph.conductance*7e-9,rtol=0,atol=0)
            self.assertLess(abs(loaded.receipt["area_relative_error"]),1e-14)
            with self.assertRaisesRegex(ValueError,"Pinned"):
                load_dual_geometry(path,expected_sha256="0"*64)

    def test_new_numeric_outputs_are_guarded_to_scratch(self):
        with self.assertRaisesRegex(ValueError,"scratch"):
            get_effective_geometry(489e-9,5e-9,output_root="/home/jdiaz/not_scratch")
        model,left,right=self.model()
        dc=prepare_effective_dc(model,left_nodes=left,right_nodes=right,
            temperature_K=.9,current_A=0.,width_m=80e-9)
        with self.assertRaisesRegex(ValueError,"scratch"):
            save_effective_dc("/home/jdiaz/not_scratch",dc,None)


if __name__=="__main__":
    unittest.main()
