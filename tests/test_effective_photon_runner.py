"""Lightweight physical source and full-state checkpoint controls for pipeline05."""
import importlib.util
import os
import json
from pathlib import Path
from types import SimpleNamespace
import tempfile
import unittest
import numpy as np
from scipy.constants import Boltzmann as kB,hbar
from pysnspd.experimental.thermal_spatial_usadel import rectangular_graph
from pysnspd.experimental.effective_bcs_functional import EffectiveBCSFunctional,bcs_count_quadrature
from pysnspd.experimental.effective_bcs_kinetics import EffectiveBCSKinetics,acoustic_phonon_grid
from pysnspd.experimental.effective_material import allmaras_acoustic_material
from pysnspd.experimental.effective_coupled import EffectiveState
ROOT=Path(os.environ.get("PYSNSPD_TEST_REPO",Path(__file__).resolve().parents[1]))
spec=importlib.util.spec_from_file_location("effective_photon_runner",ROOT/"pipelines/05_effective_photon_run.py")
runner=importlib.util.module_from_spec(spec)
candidate=os.environ.get("PYSNSPD_FRONTEND_CANDIDATES")
if candidate:exec(compile((Path(candidate)/"05_effective_photon_run_candidate.py").read_text(),str(ROOT/"pipelines/05_effective_photon_run.py"),"exec"),runner.__dict__)
else:spec.loader.exec_module(runner)

class EffectivePhotonRunnerTests(unittest.TestCase):
 def fixtures(self):
  Eu=kB*8.65;ell=np.sqrt(hbar*5e-5/(2*Eu));graph=rectangular_graph(4,3,4.,2.)
  xi,w=bcs_count_quadrature(32.,6)
  f=EffectiveBCSFunctional(graph,energy_scale_J=Eu,N0_per_J_m3=9.153297821e46,
   ell0_m=ell,thickness_m=7e-9,count_nodes=xi,count_weights=w,
   gradient_law="bcs_radial",bath_temperature_K=.9,diffusion_m2_s=5e-5)
  xy=graph.coordinates_bar*ell;xy[:,0]-=xy[:,0].min()
  left=np.flatnonzero(xy[:,0]==xy[:,0].min());right=np.flatnonzero(xy[:,0]==xy[:,0].max())
  g=SimpleNamespace(graph=graph,coordinates_m=xy,left_nodes=left,right_nodes=right,
   triangles=np.empty((0,3),int))
  ph=acoustic_phonon_grid(allmaras_acoustic_material(nodes=257),Eu,nodes=24)
  k=EffectiveBCSKinetics.from_graph(graph,ell0_m=ell,thickness_m=7e-9,
   xi_bar=xi,weights_bar=w,phonons=ph,Eu_J=Eu,N0_per_J_m3=f.N0_per_J_m3,
   D_m2_s=5e-5,bath_temperature_K=.9,transport_order=4,accelerated=False)
  d=np.full(graph.n_nodes,1.65+0j);p=f.fermi_dirac(d,.9)
  return f,k,g,d,p
 def test_defaults_are_physical_allmaras_calibrators_and_short_window(self):
  a=runner.parser().parse_args(["--output-dir","/home/jdiaz/scratch/new_run"])
  runner.validate_args(a)
  self.assertEqual((a.tau_ee_Tc_ps,a.tau_ep_Tc_ps),(5.,24.7))
  self.assertEqual((a.dt_ps,a.threshold_mV,a.confirmation_ps,a.extra_time_ps),(.02,.1,.5,5.))
  self.assertEqual((a.step_tolerance,a.local_energy_tolerance),(2e-7,1e-5))
  self.assertIn("--cpu",runner.parser().format_help())
  a.retained_fraction=1.01
  with self.assertRaises(ValueError):runner.validate_args(a)
 def test_new_source_energy_center_and_cut_reservoirs(self):
  f,k,g,d,p0=self.fixtures()
  p,n,T,target,r=runner.prepare_photon(f,k,d,p0,g,wavelength_nm=775.,
   sigma_nm=7.549834435,retained_fraction=.667,bath_temperature_K=.9)
  oldn=np.broadcast_to(k.bath_phonons,n.shape)
  actual=k.energy_J(d,p,n)-k.energy_J(d,p0,oldn)
  self.assertAlmostEqual(actual/r["retained_energy_J"],1.,delta=2e-14)
  self.assertAlmostEqual(r["center_m"][0],np.ptp(g.coordinates_m[:,0])/2,delta=1e-23)
  cuts=np.r_[g.left_nodes,g.right_nodes]
  np.testing.assert_array_equal(target[cuts],0.)
  np.testing.assert_array_equal(T[cuts],.9)
  np.testing.assert_array_equal(p[cuts],p0[cuts])
  np.testing.assert_allclose(n[cuts],oldn[cuts],rtol=3e-15)
  self.assertTrue(np.all((p>=0)&(p<=1)));self.assertTrue(np.all(n>=0))
 def test_phonon_temperature_is_energy_equivalent(self):
  f,k,g,d,p0=self.fixtures();T=np.linspace(.9,30.,g.graph.n_nodes)
  n=runner.bose_matrix(k.phonons.energies,f.Eu_J,T)
  np.testing.assert_allclose(runner.phonon_equivalent_temperature_K(k,n),T,rtol=2e-14)
 def test_atomic_checkpoint_preserves_full_state_and_histories(self):
  f,k,g,d,p=self.fixtures()
  n=np.broadcast_to(k.bath_phonons,(len(d),k.phonon_size)).copy()
  state=EffectiveState(d,p,n,np.array([21e-6,21e-6,0.]))
  base=Path("/home/jdiaz/scratch/stage5_runner_development_20261011");base.mkdir(exist_ok=True)
  with tempfile.TemporaryDirectory(dir=base) as tmp:
   output=Path(tmp);history={key:[0.] for key in runner.HISTORY_KEYS}
   snapshots={key:[np.zeros(len(d))] for key in runner.SNAPSHOT_KEYS};snapshots["time_s"]=[0.]
   runner.save_checkpoint(output,state,{"source_applied_once":True},history,snapshots,g)
   with np.load(output/"checkpoint_latest.npz",allow_pickle=False) as saved:
    for key in ("delta_bar","p","n","circuit_SI"):
     np.testing.assert_array_equal(saved[key],getattr(state,key))
    self.assertTrue(json.loads(saved["metadata_json"].item())["source_applied_once"])
    self.assertIn("trajectory__energy_defect_J",saved.files)
   self.assertFalse(any(output.glob("*.part")))
 def test_full_physical_assembly_DC_photon_RHS_and_one_step(self):
  base=Path("/home/jdiaz/scratch/stage5_runner_development_20261011");base.mkdir(exist_ok=True)
  with tempfile.TemporaryDirectory(dir=base) as tmp:
   output=Path(tmp)
   a=runner.parser().parse_args(["--output-dir",str(output),"--threads","1"])
   system,dc,g,material,params=runner.build_system(a,output)
   f,k=system.functional,system.kinetics
   self.assertLess(np.max(abs(dc.phi_V)),1e-9)
   p,n,T,target,source=runner.prepare_photon(f,k,dc.delta_bar,dc.p,g,
    wavelength_nm=775.,sigma_nm=a.sigma_nm,retained_fraction=a.retained_fraction,bath_temperature_K=.9)
   state=EffectiveState(dc.delta_bar.copy(),p,n,np.array([a.current_uA*1e-6,a.current_uA*1e-6,0.]))
   before=system.evaluate(state)
   after,receipt=system.midpoint_step(state,.005,max_iterations=40)
   observed=system.evaluate(after)
   defect=observed.energy_J-before.energy_J-receipt["external_energy_increment_J"]
   self.assertLess(abs(defect),1e-7*source["retained_energy_J"])
   self.assertLess(abs(observed.instant_power_residual_W),1e-16)
   self.assertTrue(np.all((after.p>=0)&(after.p<=1)))
   self.assertTrue(np.all(after.n>=0));self.assertTrue(np.all(after.circuit_SI[:2]>0))
   cuts=np.r_[g.left_nodes,g.right_nodes]
   np.testing.assert_array_equal(after.p[cuts],state.p[cuts])
   np.testing.assert_array_equal(after.n[cuts],state.n[cuts])
   np.testing.assert_allclose(abs(after.delta_bar[cuts]),abs(state.delta_bar[cuts]),rtol=3e-15)
   resumed,loaded_dc,loaded_g,loaded_material,loaded_params=runner.build_system(a,output,resume=True)
   np.testing.assert_array_equal(loaded_dc.delta_bar,dc.delta_bar)
   self.assertEqual(loaded_g.receipt["mesh_sha256"],g.receipt["mesh_sha256"])
   self.assertEqual(loaded_params,params)
 def test_development_output_never_claims_admission_or_latency(self):
  a=runner.parser().parse_args(["--output-dir","/home/jdiaz/scratch/new_run","--duration-ps",".05"])
  observation=runner.observation_manifest(None,a,"DEVELOPMENT_DURATION_REACHED")
  self.assertFalse(observation["trigger_confirmed"])
  self.assertIsNone(observation["trigger_time_s"])
  self.assertFalse(observation["full_reset_requested"])

if __name__=="__main__":unittest.main()
