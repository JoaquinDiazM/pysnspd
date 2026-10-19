"""Short physical/conservation regressions; no photon-latency admission."""
import os
os.environ.setdefault("NUMBA_NUM_THREADS","1")
import tempfile,unittest
from pathlib import Path
import numpy as np
from scipy.constants import Boltzmann
from pysnspd.experimental.dimensional_transition import *
from pysnspd.experimental.dimensional_archive import PhaseArchive
from pysnspd.experimental.effective_bcs_functional import EffectiveBCSFunctional,bcs_count_quadrature
from pysnspd.experimental.effective_bcs_preparation import prepare_effective_dc
from pysnspd.experimental.effective_bcs_kinetics import EffectiveBCSKinetics,acoustic_phonon_grid
from pysnspd.experimental.effective_bcs_ports import EffectiveNeutralPorts
from pysnspd.experimental.effective_material import allmaras_acoustic_material
from pysnspd.experimental.electrical_ports import ThesisCircuitParameters
from pysnspd.experimental.resummed_thermal import ResummedThermalFunctional
from pysnspd.experimental.thermal_spatial_usadel import rectangular_graph

class DimensionalTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        Eu=Boltzmann*8.65;D=1.4103757147386595e-4
        from scipy.constants import hbar,elementary_charge
        ell=np.sqrt(hbar*D/(2*Eu));df=7e-9;W=80e-9;L=80e-9
        N0=1/(2*elementary_charge**2*D*425*df)
        graph=rectangular_graph(9,3,L/ell,W/ell)
        xy=graph.coordinates_bar*ell;ed=graph.edges
        dx=np.linalg.norm(xy[ed[:,1]]-xy[ed[:,0]],axis=1)
        left=np.where(xy[:,0]==xy[:,0].min())[0];right=np.where(xy[:,0]==xy[:,0].max())[0]
        graph=ThermalGraph(graph.area_weights,ed,graph.conductance,graph.coordinates_bar,np.r_[left,right])
        # Structured triangle fixture is only used for the winding diagnostic.
        triangles=[]
        for j in range(2):
            for i in range(8):
                v=j*9+i;triangles.extend([[v,v+1,v+9],[v+1,v+10,v+9]])
        g2=EffectiveBCSGeometry(graph,xy,np.array(triangles),dx,graph.conductance*dx,
            left,right,ell,L,W,"fixture",{})
        g1=chain_geometry(L,W,10e-9,ell)
        xi,w=bcs_count_quadrature(12,3)
        material=allmaras_acoustic_material(Tc_K=8.65,tau0_s=1.8e-9,diffusivity_m2_s=D,
            sheet_resistance_ohm=425,thickness_m=df,nodes=257)
        phonons=acoustic_phonon_grid(material,Eu,nodes=24)
        systems=[];states=[]
        for g in (g2,g1):
            kw=dict(graph=g.graph,energy_scale_J=Eu,N0_per_J_m3=N0,ell0_m=ell,
                thickness_m=df,gradient_law="bcs_radial",count_nodes=xi,count_weights=w,bath_temperature_K=.9,diffusion_m2_s=D)
            f=ResummedThermalFunctional(**kw);bath=EffectiveBCSFunctional(**kw)
            dc=prepare_effective_dc(bath,left_nodes=g.left_nodes,right_nodes=g.right_nodes,
                temperature_K=.9,current_A=15.5e-6,width_m=W,sheet_resistance_ohm=425)
            k=EffectiveBCSKinetics(xi_bar=xi,weights_bar=w,phonons=phonons,Eu_J=Eu,
                N0_per_J_m3=N0,D_m2_s=D,volumes_m3=f.volumes_m3,edges=g.graph.edges,
                face_areas_m2=g.face_areas_m2(df),edge_lengths_m=g.edge_lengths_m,bath_temperature_K=.9,
                tau_ee_Tc_s=4e-12,tau_escape_s=20e-12,transport_order=8,numba_threads=1)
            ports=EffectiveNeutralPorts(g.graph,sheet_resistance_ohm=425,left_nodes=g.left_nodes,right_nodes=g.right_nodes)
            s=observed(ThermalCoupledSystem(f,k,ports,ThesisCircuitParameters(
                Lk_ext_H=100e-9,V_bias_V=.155),tau_ee_Tc_ps=4,tau_ep_Tc_ps=24.7,
                current_scale_A=15.5e-6,kinetic_residual_energy_scale_J=1e-19))
            state=EffectiveState(dc.delta_bar,dc.p,np.broadcast_to(k.bath_phonons,(g.graph.n_nodes,k.phonon_size)).copy(),
                                 np.array([15.5e-6,15.5e-6,0.]))
            systems.append(s);states.append(state)
        cls.s2,cls.s1=systems;cls.dc2,cls.dc1=states;cls.g2,cls.g1=g2,g1
        cls.projector=DimensionalProjection(cls.s2,g2,cls.dc2,cls.s1,g1,cls.dc1)
        cls.tol=TransitionTolerances()

    def test_clipped_voronoi_boundary_volumes_and_faces(self):
        from pysnspd.experimental.dimensional_geometry import clipped_voronoi
        xy=np.array([[0,-1],[2,-1],[0,1],[2,1]])*1e-9
        edges=np.array([[0,1],[0,2],[1,3],[2,3],[0,3]])
        area,face=clipped_voronoi(xy,edges,2e-9,2e-9)
        np.testing.assert_allclose(area,1e-18,rtol=1e-12)
        np.testing.assert_allclose(face[:4],1e-9,rtol=1e-12)
        self.assertLess(face[-1],1e-20)

    def test_volume_and_constant_marginals(self):
        R=self.projector.R
        np.testing.assert_allclose(R@np.ones(self.g2.graph.n_nodes),1,atol=1e-12)
        np.testing.assert_allclose(self.s1.functional.volumes_m3@R,self.s2.functional.volumes_m3,rtol=1e-12)
        self.assertTrue(np.all(R>=0))
        with self.assertRaises(ValueError):volume_restriction(np.arange(2),[1,1],[1,2])

    def test_exact_dc_and_circuit(self):
        dc=self.dc2.copy();dc.circuit_SI=np.array([15e-6,14e-6,.001])
        result,report=self.projector.project(dc,1e-19,self.tol)
        np.testing.assert_array_equal(result.delta_bar,self.dc1.delta_bar)
        np.testing.assert_array_equal(result.circuit_SI,dc.circuit_SI)
        self.assertLess(report["film_energy_relative"],1e-10)

    def test_dc_relative_phase_is_not_reset(self):
        st=self.dc2.copy();st.delta_bar*=np.exp(.23j)
        result,report=self.projector.project(st,1e-19,self.tol)
        np.testing.assert_allclose(result.delta_bar,self.dc1.delta_bar*np.exp(.23j),rtol=1e-9,atol=1e-9)
        self.assertLess(report["film_energy_relative"],2e-9)

    def test_uniform_transverse_nonequilibrium_retained(self):
        st=self.dc2.copy();f=self.s2.functional
        # Native non-FD spectral perturbation, uniform across y at each x.
        x=self.g2.coordinates_m[:,0];profile=np.sin(np.pi*x/self.g2.length_m)**2
        st.p+=profile[:,None]*np.exp(-f.xi_bar[None,:])*1e-5
        st.n+=profile[:,None]*np.exp(-self.s2.kinetics.phonons.energies[None,:])*1e-5
        result,report=self.projector.project(st,1e-19,self.tol)
        self.assertLess(report["film_energy_relative"],2e-9)
        self.assertLess(report["gap_relative"],1e-6)
        self.assertGreater(np.max(abs(result.p-self.dc1.p)),1e-7)
        np.testing.assert_array_equal(result.circuit_SI,st.circuit_SI)

    def test_same_continuous_zero_mode_rhs(self):
        a=self.dc2.copy();b=self.dc1.copy()
        x2=self.g2.coordinates_m[:,0];x1=self.g1.coordinates_m[:,0]
        for st,ss,x,g in ((a,self.s2,x2,self.g2),(b,self.s1,x1,self.g1)):
            profile=np.sin(np.pi*x/g.length_m)**2
            st.p+=profile[:,None]*np.exp(-ss.functional.xi_bar[None,:])*1e-6
            st.n+=profile[:,None]*np.exp(-ss.kinetics.phonons.energies[None,:])*1e-6
        ra=self.s2.evaluate(a);rb=self.s1.evaluate(b)
        np.testing.assert_allclose(self.projector.R@ra.dp_per_ps,rb.dp_per_ps,rtol=2e-5,atol=1e-10)
        np.testing.assert_allclose(self.projector.R@ra.dn_per_ps,rb.dn_per_ps,rtol=2e-5,atol=1e-10)
        self.assertLess(abs(ra.domain_energy_J-rb.domain_energy_J),1e-27)
        self.assertLess(abs(ra.electrical.Vdev_V-rb.electrical.Vdev_V),1e-9)

    def test_midpoint_instrumentation_and_temporal_guard(self):
        st=self.dc1.copy();s=self.s1
        st.p[s.free,0]+=1e-7
        rhs=s.evaluate(st)
        kw=dict(tolerance=2e-7,max_iterations=60,energy_tolerance=2.5e-6,
                reference=1e-19,absolute_work=0.)
        new,rr,rec=doubled_step(s,st,rhs,1e-3,time_rtol=2e-4,time_atol=2e-6,**kw)
        self.assertLessEqual(rec["temporal_error"],1)
        self.assertGreaterEqual(rec["absolute_cut_energy_increment_J"],0)
        self.assertGreaterEqual(rec["absolute_device_work_increment_J"],0)
        self.assertTrue(np.all((new.p>=0)&(new.p<=1)))
        self.assertTrue(np.all(new.n>=0))

    def test_phase_archive_commit_truncation(self):
        with tempfile.TemporaryDirectory(dir="/home/jdiaz/scratch/tmp/pysnspd") as d:
            path=Path(d)/"test.h5"
            systems={"2d":self.s2,"shadow1d":self.s1,"1d":self.s1}
            geometries={"2d":self.g2,"shadow1d":self.g1,"1d":self.g1}
            dcs={"2d":self.dc2,"shadow1d":self.dc1,"1d":self.dc1}
            ar=PhaseArchive(path,systems,geometries,dcs)
            ar.append("2d",dict(time_ps=0));counts=ar.counts()
            ar.append("2d",dict(time_ps=1));ar.close()
            ar=PhaseArchive(path,systems,geometries,dcs,resume=True,counts=counts)
            self.assertEqual(len(ar.file["2d/trace"]),1);ar.close()

    def test_runner_handoff_and_exact_resume(self):
        import argparse,importlib.util,json
        from unittest.mock import patch
        from types import SimpleNamespace
        spec=importlib.util.spec_from_file_location("_transition_runner_test",
            "/home/jdiaz/pysnspd/pipelines/13_dimensional_photon_run.py")
        runner=importlib.util.module_from_spec(spec);spec.loader.exec_module(runner)
        original=runner.B.TriggerWindow
        from sandbox.stage4_core.parallel_runtime import linux_resources
        prior_affinity=os.sched_getaffinity(0)
        os.sched_setaffinity(0,os.sched_getaffinity(1))
        inventory=linux_resources()
        os.sched_setaffinity(0,prior_affinity)
        class ArtificialVoltageProbe(original):
            def observe_accepted(self,time_s,V_out_V):
                # Wiring test only: artificial probe, actual DC equations unchanged.
                return super().observe_accepted(time_s,self.config.threshold_V*time_s/1e-15)
        dc=SimpleNamespace(delta_bar=self.dc2.delta_bar,p=self.dc2.p)
        with tempfile.TemporaryDirectory(dir="/home/jdiaz/scratch/tmp/pysnspd") as directory:
            out=Path(directory)/"run"
            words=["--output-dir",str(out),"--threads","1","--physical-cores","1","--cpu-start","12",
                "--dc-only","--duration-ps",".006","--dt-ps",".001",
                "--confirmation-ps",".0005","--extra-time-ps",".002",
                "--snapshot-every-ps",".001","--shadow-check-every-ps",".0005",
                "--xi-cutoff","12","--xi-order","3","--phonon-nodes","24","--transport-order","8",
                "--current-uA","15.5","--spacing-1d-nm","10"]
            build=lambda *args,**kw:(self.s2,dc,self.g2,SimpleNamespace(manifest={}),{})
            # 1D DC cache/resume follows the real implementation.
            with patch.object(runner.OLD,"build_system",side_effect=build),patch.object(runner.B,"TriggerWindow",ArtificialVoltageProbe),patch("sandbox.stage4_core.parallel_runtime.linux_resources",return_value=inventory):
                code=runner.run(runner.parser().parse_args(words))
                self.assertEqual(code,0)
                m=json.loads((out/"manifest.json").read_text())
                self.assertEqual(m["phase"],"1d");self.assertGreaterEqual(len(m["shadow_checks"]),2)
                words[words.index(".006")]=".009";words.append("--resume")
                code=runner.run(runner.parser().parse_args(words))
                self.assertEqual(code,0)
                m=json.loads((out/"manifest.json").read_text())
                self.assertAlmostEqual(m["time_ps"],.009)
                self.assertLess(abs(m["energy_defect_J"]),1e-25)
                # Authentication must reject a modified discrete reference.
                with (out/"dc1d.npz").open("ab") as f:f.write(b"changed")
                with self.assertRaisesRegex(ValueError,"1D DC reference changed"):
                    runner.run(runner.parser().parse_args(words))

    def test_campaign_shared_budget_and_busy_refusal(self):
        import importlib.util
        spec=importlib.util.spec_from_file_location("_dimensional_campaign_test",
            "/home/jdiaz/pysnspd/pipelines/14_dimensional_basecases.py")
        module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
        inventory=dict(logical_cpus=list(range(32)),
            physical_core_groups=[[i,i+16] for i in range(16)],
            available_memory_bytes=64*1024**3,total_memory_bytes=128*1024**3,cpu_quota_units=None)
        plan=module.allocation(inventory,[])
        self.assertEqual(plan["worker_threads"],24)
        self.assertGreaterEqual(len(plan["reserved"]),2)
        used=[c for slot in plan["slots"] for g in slot["core_groups"] for c in g]
        self.assertEqual(len(used),len(set(used)))
        with self.assertRaises(ValueError):module.allocation(inventory,[dict(affinity_cpus=list(range(12)))])
        with self.assertRaises(ValueError):module.allocation(dict(inventory,available_memory_bytes=8*1024**3),[])

    def test_vortex_pair_does_not_cancel(self):
        z=np.array([1,np.exp(2j*np.pi/3),np.exp(4j*np.pi/3),
                    1,np.exp(-2j*np.pi/3),np.exp(-4j*np.pi/3)])
        self.assertEqual(resolved_winding(z,np.array([[0,1,2],[3,4,5]]),np.ones(6)),2)
        self.assertEqual(resolved_winding(np.zeros(6,complex),np.array([[0,1,2]]),np.ones(6)),0)

if __name__=="__main__":unittest.main()
