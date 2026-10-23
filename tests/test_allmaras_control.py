"""Isolated historical inputs and concurrent allocation; no latency admission."""
import importlib.util
from pathlib import Path
from types import SimpleNamespace
import unittest
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
def load(name,path):
    spec=importlib.util.spec_from_file_location(name,path)
    m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m
P=load('_historical_producer_test',ROOT/'pipelines/16_allmaras_parameter_control.py')
S=load('_historical_scheduler_test',ROOT/'pipelines/17_allmaras_reference_cases.py')

class HistoricalControlTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.runner=P.configure(ROOT/'pipelines/13_dimensional_photon_run.py',1)

    def resources(self):
        return dict(logical_cpus=list(range(32)),physical_core_groups=[[i,i+16] for i in range(16)],
            available_memory_bytes=70*1024**3,total_memory_bytes=128*1024**3,cpu_quota_units=None)

    def records(self,last608=False):
        old=Path('/home/jdiaz/scratch/photon_runs/identified_original')
        rows=[]
        for index,cpus in enumerate([[0,16,1,17],[4,20,5],[7,23,8],[21,6,22]]):
            rows.append(dict(pid=100+index,affinity_cpus=cpus,
                words=['python','pipelines/13_dimensional_photon_run.py','--output-dir',str(old/'runs'/str(index))]))
        rows.append(dict(pid=105,affinity_cpus=[13,29],
            words=['python','pipelines/15_calibration_basecases.py','--campaign-dir',str(old)]))
        if last608:rows.append(dict(pid=106,affinity_cpus=[27,12,28],
            words=['python','pipelines/13_dimensional_photon_run.py','--output-dir',str(old/'runs'/'last608')]))
        return old,rows

    def test_combined_resource_masks_preserve_current_workers_and_two_whole_cores(self):
        old,rows=self.records();r=S.inventory(old,self.resources(),rows)
        cpus=[v for values in r['slots'].values() for v in values]
        busy={v for row in rows for v in row['affinity_cpus']}
        self.assertEqual(len(cpus),12);self.assertEqual(len(set(cpus)),12)
        self.assertFalse(set(cpus)&busy)
        self.assertTrue(r['coordinator_shared_with_existing'])
        self.assertLessEqual(len(r['combined_admitted_logical_cpus']),28)
        self.assertGreaterEqual(len(r['reserved_whole_core_groups']),2)

    def test_early_or_unidentified_campaign_is_refused_without_repinning(self):
        old,rows=self.records(True)
        with self.assertRaises(ValueError):S.inventory(old,self.resources(),rows)
        old,rows=self.records();rows[0]['words'][-1]='/home/jdiaz/scratch/unidentified/run'
        with self.assertRaisesRegex(ValueError,'Unidentified'):S.inventory(old,self.resources(),rows)

    def test_hotbelt_has_exact_volume_uniform_interior_and_no_cut_deposition(self):
        from scipy.spatial import Delaunay
        from pysnspd.experimental.dimensional_geometry import clipped_voronoi
        xy=np.array([[x,y] for x in np.arange(0,301,10) for y in np.arange(-40,41,20)],float)*1e-9
        triangles=Delaunay(xy).simplices
        edges=np.array(sorted({tuple(sorted((int(a),int(b)))) for tri in triangles
            for a,b in [(tri[0],tri[1]),(tri[1],tri[2]),(tri[2],tri[0])]}))
        areas,_=clipped_voronoi(xy,edges,300e-9,80e-9);vol=areas*7e-9
        g=SimpleNamespace(coordinates_m=xy,width_m=80e-9,graph=SimpleNamespace(edges=edges),
            left_nodes=np.flatnonzero(xy[:,0]==0),right_nodes=np.flatnonzero(xy[:,0]==xy[:,0].max()))
        weights=P.hotbelt_weights(g,vol,7e-9)
        self.assertAlmostEqual(weights.sum(),1.)
        self.assertEqual(weights[g.left_nodes].sum()+weights[g.right_nodes].sum(),0.)
        interior=(xy[:,0]>135e-9)&(xy[:,0]<165e-9)
        np.testing.assert_allclose(weights[interior]/vol[interior],1/(40e-9*80e-9*7e-9),rtol=1e-10)
        self.assertTrue(np.all(weights[(xy[:,0]<125e-9)|(xy[:,0]>175e-9)]==0))

    def test_fixed_D_is_not_silently_recalibrated_to_measured_depairing(self):
        from scipy.constants import elementary_charge as e
        c=self.runner.OLD.calibration_result(bulk_factory=self.runner.ResummedThermalFunctional,
            sheet_resistance_ohm=608.,target_depairing_current_A=38.8*80/120*1e-6,
            width_m=80e-9,requested_current_A=21.5e-6,bath_temperature_K=2.,
            critical_temperature_K=8.65,count_cutoff_bar=32.,count_order=6)
        self.assertEqual(c['diffusion_m2_s'],5e-5)
        self.assertFalse(c['measured_depairing_fit_enforced'])
        self.assertFalse(c['calibration']['diffusion_is_fitted_effective_parameter'])
        self.assertNotAlmostEqual(c['achieved_depairing_current_A']*1e6,38.8*80/120,places=3)
        self.assertAlmostEqual(c['N0_per_J_m3']/(1/(2*e**2*5e-5*608*7e-9)),1.)

    def test_gamma_override_changes_DOS_and_capacity_but_not_electronic_alpha2F(self):
        from pysnspd.experimental.effective_material import allmaras_acoustic_material
        original=allmaras_acoustic_material(Tc_K=8.65,tau0_s=1.868435901588345e-9,
            diffusivity_m2_s=5e-5,sheet_resistance_ohm=608.,thickness_m=7e-9,nodes=129)
        changed=self.runner.B.allmaras_acoustic_material(Tc_K=8.65,tau0_s=1.868435901588345e-9,
            diffusivity_m2_s=5e-5,sheet_resistance_ohm=608.,thickness_m=7e-9,nodes=129)
        scale=original.manifest['gamma_A20']/60.
        np.testing.assert_array_equal(original.alpha2F,changed.alpha2F)
        np.testing.assert_allclose(changed.g_ph_per_J_m3,original.g_ph_per_J_m3*scale)
        self.assertEqual(changed.N0_per_J_m3,original.N0_per_J_m3)
        self.assertAlmostEqual(changed.heat_capacity_J_per_m3_K(8.65)/original.heat_capacity_J_per_m3_K(8.65),scale)
        self.assertFalse(changed.manifest['acoustic_speed_alone_determines_DOS'])
        self.assertEqual(changed.manifest['gamma_A20'],60.)

if __name__=='__main__':unittest.main()
