"""Control-flow regressions of effective-run interruption and resume.
No physical photon trajectory is simulated; geometry checks use small mesh builds.
"""
import contextlib,importlib.util,io,json,os,signal,sys,tempfile,types,unittest
from pathlib import Path
from types import SimpleNamespace as NS
from unittest.mock import patch
import numpy as np
ROOT=Path(os.environ.get("PYSNSPD_TEST_REPO",Path(__file__).resolve().parents[1]))
CANDIDATES=os.environ.get("PYSNSPD_FRONTEND_CANDIDATES")
BASE=Path(os.environ.get("PYSNSPD_SCRATCH_ROOT","/home/jdiaz/scratch"))/"tmp/pysnspd/frontend_transaction_tests"
BASE.mkdir(parents=True,exist_ok=True)
def source_path(filename,original):
    return Path(CANDIDATES)/filename if CANDIDATES else original
sys.path.insert(0,str(ROOT))
from pysnspd.experimental.effective_coupled import EffectiveState
import sandbox.stage4_core.parallel_runtime as runtime

def candidate_module(name,filename,original,package=None):
    module=types.ModuleType(name);module.__file__=str(original)
    if package is not None:module.__package__=package
    sys.modules[name]=module
    exec(compile(source_path(filename,original).read_text(),str(original),"exec"),module.__dict__)
    return module

class FrontendTransactionTests(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory(dir=BASE)
        self.work=Path(self.temp.name)
        self.runner=candidate_module("candidate_runner_"+self.work.name,
            "05_effective_photon_run_candidate.py",ROOT/"pipelines/05_effective_photon_run.py")
        self.r=self.runner
        self.hard_interrupt_line=next(i for i,s in enumerate(
            source_path("05_effective_photon_run_candidate.py",ROOT/"pipelines/05_effective_photon_run.py").read_text().splitlines(),1)
            if "external+=float(receipt" in s)
        self.max_step=None;self.signal_after_accepts=None;self.signalled=False
        graph=NS(n_nodes=2,edges=np.empty((0,2),int))
        f=NS(graph=graph,volumes_m3=np.ones(2),Eu_J=1.,
            N0_per_J_m3=1.,weights_bar=np.ones(1),node_energy_scale_J=np.ones(2))
        f.evaluate=lambda d,p:NS(node_energy_J=np.zeros(2))
        self.fake_k=NS(bath_phonons=np.ones(1),phonon_size=1,electron_size=1,
            phonon_energy_weights_J=np.ones((2,1)),accelerated=False,numba_threads=1)
        self.fake_k.metadata=lambda:dict(accelerated=False)
        self.g=NS(graph=graph,receipt=dict(control_flow_mock=True),
            coordinates_m=np.zeros((2,2)),triangles=np.empty((0,3),int))
        self.dc=NS(delta_bar=np.ones(2,complex),p=np.full((2,1),.1))
        owner=self
        class FakeSystem:
            functional=f;kinetics=owner.fake_k;circuit=NS(R_load_ohm=50.)
            kinetic_residual_energy_scale_J=1.
            def evaluate(self,state):
                return NS(energy_J=float(state.delta_bar.real.sum()),
                    domain_energy_J=0.,escape_power_W=0.,reservoir_power_out_W=0.,
                    instant_power_residual_W=0.,electrical=NS(Vdev_V=0.))
            def midpoint_step(self,state,h_ps,**kw):
                if owner.max_step is not None and h_ps>owner.max_step*(1+1e-12):
                    raise ArithmeticError("Mock convergence-only rejection")
                out=state.copy()
                out.delta_bar+=h_ps
                out.circuit_SI[1]-=h_ps*50e-6
                energy=float(out.delta_bar.real.sum()-state.delta_bar.real.sum())
                count=int(round((out.delta_bar.real.mean()-1.)/.01))
                if owner.signal_after_accepts is not None and count>=owner.signal_after_accepts and not owner.signalled:
                    owner.signalled=True
                    os.kill(os.getpid(),signal.SIGTERM)
                return out,dict(external_energy_increment_J=energy,
                    limited_event_energy_J=0.,totalabs_event_energy_J=0.)
        self.system=FakeSystem()
        def build(a,out,resume=False):
            self.fake_k.numba_threads=a.threads
            if not resume:
                (out/"dc").mkdir(exist_ok=True)
                for name in ("dc/effective_dc.npz","material.npz","material.json"):
                    (out/name).write_text("CONTROL_FLOW_MOCK")
            return self.system,self.dc,self.g,NS(manifest={}),{}
        self.r.build_system=build
        self.r.prepare_photon=lambda *a,**kw:(np.full((2,1),.1),np.ones((2,1)),
            np.ones(2),np.zeros(2),dict(retained_energy_J=1.))
        def snapshot(snapshots,system,state,rhs,t_ps,reference):
            for key in self.r.SNAPSHOT_KEYS:
                snapshots[key].append(t_ps*1e-12 if key=="time_s" else np.zeros(2))
        self.r.append_snapshot=snapshot
        self.r.source_pins=lambda:dict(control_flow_mock="stable-source")
    def tearDown(self):
        sys.settrace(None)
        self.temp.cleanup()
    def argv(self,out,extra=()):
        return ["--output-dir",str(out),"--threads","1","--duration-ps",".12",
            "--dt-ps",".01","--confirmation-ps",".005","--extra-time-ps",".03"]+list(extra)
    @contextlib.contextmanager
    def resources(self):
        with patch.object(runtime,"linux_resources",return_value=dict(
                physical_core_groups=[[13],[14],[15],[16]])),\
             patch.object(runtime,"resource_budget",return_value=dict(
                selected_affinity_cpus=[13,14,15,16])),\
             patch.object(self.r.os,"sched_getaffinity",return_value={13}),\
             patch.object(self.r.os,"sched_setaffinity"),contextlib.redirect_stdout(io.StringIO()):
            yield
    def execute(self,out,extra=(),main=False):
        argv=self.argv(out,extra)
        with self.resources():
            if main:return self.r.main(argv)
            return self.r.run(self.r.parser().parse_args(argv))
    def checkpoint(self,out):
        with np.load(out/"checkpoint_latest.npz",allow_pickle=False)as z:
            return {k:z[k].copy() for k in z.files}
    def assert_coherent(self,out):
        saved=self.checkpoint(out);meta=json.loads(saved["metadata_json"].item())
        t=meta["time_ps"]*1e-12
        for key in self.r.HISTORY_KEYS:
            self.assertEqual(len(saved["trajectory__"+key]),len(saved["trajectory__time_s"]))
        self.assertEqual(float(saved["trajectory__time_s"][-1]),t)
        self.assertEqual(float(saved["snapshots__time_s"][-1]),t)
        self.assertEqual(meta["controller"]["state"]["last_time_s"],t)
        self.assertEqual(meta["external_energy_integral_J"],
            float(saved["trajectory__external_energy_integral_J"][-1]))
        self.assertEqual(meta["accepted_steps"],len(saved["trajectory__time_s"])-1)
        return saved,meta
    def assert_same_as_reference(self,out,reference):
        left=self.checkpoint(out);right=self.checkpoint(reference)
        for key in left:
            if key!="metadata_json" and not key.startswith("snapshots__"):
                np.testing.assert_array_equal(left[key],right[key],err_msg=key)
        # An interruption adds a diagnostic snapshot at its accepted boundary.
        # All scheduled snapshots and all physical trajectory entries remain exact.
        lt=left["snapshots__time_s"];rt=right["snapshots__time_s"]
        for i,t in enumerate(rt):
            matched=np.flatnonzero(lt==t)
            self.assertEqual(len(matched),1)
            for key in right:
                if key.startswith("snapshots__"):
                    np.testing.assert_array_equal(left[key][matched[0]],right[key][i],err_msg=key)
        lm=json.loads(left["metadata_json"].item());rm=json.loads(right["metadata_json"].item())
        for key in ("time_ps","external_energy_integral_J","accepted_steps","total_rejections",
                "stable_accepted_steps","consecutive_rejections","next_step_ps","controller"):
            self.assertEqual(lm[key],rm[key],key)
    def test_hard_interrupt_between_state_and_ledger_rolls_back_then_resumes_exactly(self):
        out=self.work/"interrupt";reference=self.work/"reference"
        def trace(frame,event,arg):
            if event=="line" and frame.f_code.co_filename==self.r.__file__ and frame.f_lineno==self.hard_interrupt_line:
                sys.settrace(None);raise KeyboardInterrupt("Mock immediate interrupt")
            return trace
        try:
            sys.settrace(trace)
            with self.assertRaises(KeyboardInterrupt):self.execute(out)
        finally:sys.settrace(None)
        _,meta=self.assert_coherent(out)
        self.assertEqual(meta["accepted_steps"],0)
        self.execute(out,["--resume"]);self.execute(reference)
        self.assert_same_as_reference(out,reference)
    def test_partial_history_append_is_rolled_back_for_every_field(self):
        out=self.work/"partial"
        original=self.r.append_history;count=0
        def append(history,row):
            nonlocal count
            count+=1
            if count==2:
                for key in self.r.HISTORY_KEYS[:3]:history[key].append(float(row[key]))
                raise KeyboardInterrupt("Mock interruption inside history append")
            original(history,row)
        with patch.object(self.r,"append_history",side_effect=append):
            with self.assertRaises(KeyboardInterrupt):self.execute(out)
        _,meta=self.assert_coherent(out)
        self.assertEqual(meta["accepted_steps"],0)
        self.execute(out,["--resume"])
    def test_sigterm_is_deferred_until_state_ledger_and_trigger_commit(self):
        out=self.work/"signal";reference=self.work/"reference"
        self.signal_after_accepts=5
        with self.assertRaises(KeyboardInterrupt):self.execute(out,main=True)
        _,meta=self.assert_coherent(out)
        self.assertGreaterEqual(meta["accepted_steps"],5)
        self.assertIsNotNone(meta["controller"]["state"]["candidate_crossing_s"])
        self.signal_after_accepts=None
        self.execute(out,["--resume"]);self.execute(reference)
        self.assert_same_as_reference(out,reference)
    def test_resume_preserves_adaptive_regrowth_counters(self):
        out=self.work/"adaptive";reference=self.work/"reference"
        self.max_step=.01;self.signal_after_accepts=6
        extra=["--dt-ps",".02","--threshold-mV","10"]
        with self.assertRaises(KeyboardInterrupt):self.execute(out,extra,main=True)
        _,meta=self.assert_coherent(out)
        self.assertEqual(meta["stable_accepted_steps"],6)
        self.assertEqual(meta["total_rejections"],1)
        self.signal_after_accepts=None
        self.execute(out,extra+["--resume"]);self.execute(reference,extra)
        self.assert_same_as_reference(out,reference)
    def test_numpy_fallback_records_actual1_and_reserves_requested4(self):
        out=self.work/"fallback"
        manifest=self.execute(out,["--threads","4"])
        self.assertEqual(manifest["resources"]["actual_threads"],1)
        self.assertEqual(manifest["resources"]["requested_threads"],4)
        self.assertEqual(manifest["resources"]["reserved_physical_core_count"],4)
        self.assertIn("Serial NumPy",manifest["resources"]["threading"])

class GeometryIdentityTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.prep=candidate_module("pysnspd.experimental._review_preparation",
            "effective_bcs_preparation_candidate.py",
            ROOT/"pysnspd/experimental/effective_bcs_preparation.py","pysnspd.experimental")
    def test_80nm_references_remain_bitwise_identical(self):
        for length,spacing in ((160.,5.),(160.,3.5),(320.,5.)):
            with self.subTest(length=length,spacing=spacing):
                current=__import__("pysnspd.experimental.effective_bcs_preparation",fromlist=["get_effective_geometry"])
                old=current.load_dual_geometry(current.HISTORICAL_MESHES[(length,spacing)])
                new=self.prep.get_effective_geometry(length*1e-9,spacing*1e-9)
                self.assertEqual(new.receipt["mesh_sha256"],old.receipt["mesh_sha256"])
                for attr in ("coordinates_m","triangles","left_nodes","right_nodes","edge_lengths_m"):
                    np.testing.assert_array_equal(getattr(new,attr),getattr(old,attr))
                self.assertEqual(new.ell0_m,old.ell0_m)
    def test_85nm_does_not_reuse_80nm_historical_mesh(self):
        with tempfile.TemporaryDirectory(dir=BASE)as tmp:
            g=self.prep.get_effective_geometry(160e-9,5e-9,width_m=85e-9,output_root=tmp)
            np.testing.assert_allclose(g.width_m,85e-9,rtol=1e-13,atol=0.)
            self.assertNotEqual(g.source_path,str(self.prep.BASE_MESH))
    def test_historical_scale_rejects_changed_D_or_Tc(self):
        for kw in (dict(diffusion_m2_s=6e-5),dict(Tc_K=9.),
                dict(diffusion_m2_s=5e-5*(1+1e-6)),dict(Tc_K=8.65*(1+1e-6))):
            with self.subTest(**kw),self.assertRaisesRegex(ValueError,"different physical length scale"):
                self.prep.get_effective_geometry(160e-9,5e-9,**kw)
    def test_cached_scale_rejects_changed_D_or_Tc(self):
        with tempfile.TemporaryDirectory(dir=BASE)as tmp:
            self.prep.get_effective_geometry(120e-9,10e-9,output_root=tmp)
            for kw in (dict(diffusion_m2_s=6e-5),dict(Tc_K=9.),
                dict(diffusion_m2_s=5e-5*(1+1e-6)),dict(Tc_K=8.65*(1+1e-6))):
                with self.subTest(**kw),self.assertRaisesRegex(ValueError,"different physical length scale"):
                    self.prep.get_effective_geometry(120e-9,10e-9,output_root=tmp,**kw)

if __name__=="__main__":unittest.main()
