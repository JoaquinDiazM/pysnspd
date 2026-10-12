"""Small physical-unit and observation-label controls for visual pipeline 05."""
import importlib.util
from pathlib import Path
import unittest
from unittest.mock import patch
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location("effective_photon_visual",ROOT/"plot_pipelines/05_plot_effective_photon_run.py")
visual=importlib.util.module_from_spec(spec)
spec.loader.exec_module(visual)


class EffectiveVisualFormatTests(unittest.TestCase):
    def fixtures(self):
        t=np.array([0.,1e-12,3e-12])
        history={"time_s":t}
        for key in ("Vout_V","Vdev_V","Ib_A","Is_A","vc_V","energy_total_J",
                "energy_defect_J","escape_power_W","domain_energy_J","instant_power_residual_W",
                "external_energy_integral_J"):
            history[key]=np.zeros(3)
        history["Vout_V"]=np.array([2e-6,1.2e-4,4e-5])
        history["Ib_A"]=np.array([21e-6,20e-6,19e-6])
        history["Is_A"]=np.array([21e-6,18e-6,17e-6])
        history["energy_defect_J"]=np.array([0.,2e-22,-1e-22])
        manifest={"source":{"retained_energy_J":2e-19},"observation":{"baseline_Vout_V":2e-6}}
        maps={"coordinates_m":np.array([[0.,0.],[1e-9,0.],[0.,1e-9],[1e-9,1e-9]]),
            "time_s":t,"delta_J":np.full((3,4),visual.MEV_J,complex),
            "Te_K":np.full((3,4),.9),"Tph_K":np.full((3,4),.9),
            "u_excess_J_m3":np.zeros((3,4))}
        return manifest,history,maps

    def test_physical_units_and_retained_energy_normalization(self):
        manifest,history,maps=self.fixtures()
        h,m=visual.prepare_view(manifest,history,maps)
        np.testing.assert_allclose(h["time_ps"],[0.,1.,3.])
        np.testing.assert_allclose(h["Vout_relative_mV"],[0.,.118,.038])
        np.testing.assert_allclose(h["Ib_uA"],[21.,20.,19.])
        np.testing.assert_allclose(h["energy_defect_percent_retained"],[0.,.1,-.05])
        np.testing.assert_allclose(m["gap_meV"],1.)
        np.testing.assert_allclose(m["coordinates_nm"],maps["coordinates_m"]*1e9)
        np.testing.assert_array_equal(history["Vout_V"],[2e-6,1.2e-4,4e-5])

    def test_global_energy_normalization_uses_absolute_signed_work_increments(self):
        manifest,history,maps=self.fixtures()
        retained=manifest["source"]["retained_energy_J"]
        history["external_energy_integral_J"]=np.array([0.,retained,0.])
        h,_=visual.prepare_view(manifest,history,maps)
        np.testing.assert_allclose(h["absolute_external_work_J"],[0.,retained,2*retained])
        np.testing.assert_allclose(h["energy_defect_percent_global"],[0.,.05,.05/3])
        # The separate signed/Eret curve is retained, not silently relabeled.
        np.testing.assert_allclose(h["energy_defect_percent_retained"],[0.,.1,-.05])
        history.pop("external_energy_integral_J")
        with self.assertRaises(KeyError):
            visual.prepare_view(manifest,history,maps)

    def test_absolute_external_work_keeps_nonzero_first_saved_integral(self):
        manifest,history,maps=self.fixtures()
        retained=manifest["source"]["retained_energy_J"]
        history["external_energy_integral_J"]=retained*np.array([-.5,.5,-.5])
        history["energy_defect_J"]=np.array([1e-22,2e-22,-1e-22])
        h,_=visual.prepare_view(manifest,history,maps)
        np.testing.assert_allclose(h["absolute_external_work_J"],retained*np.array([.5,1.5,2.5]))
        np.testing.assert_allclose(h["energy_defect_percent_global"],[.05/1.5,.1/2.5,.05/3.5])

    def test_snapshots_keep_initial_final_and_nearest_unique_saved_times(self):
        t=np.array([0.,.2e-12,.7e-12,2e-12,5e-12])
        np.testing.assert_array_equal(visual.nearest_snapshots(t,[.3,.6,1.8]),[0,1,2,3,4])
        np.testing.assert_array_equal(visual.nearest_snapshots(t,[.3,.31]),[0,1,4])
        with self.assertRaisesRegex(ValueError,"outside"):
            visual.nearest_snapshots(t,[6.])

    def test_no_confirmed_trigger_is_censored_without_invented_latency(self):
        m={"observation":{"status":"CENSORED_TIME_CEILING","trigger_confirmed":False,"trigger_time_s":1e-12}}
        text=visual.observation_text(m,200e-12)
        self.assertIn("censura",text)
        self.assertIn("200 ps",text)
        self.assertNotIn("Gatillo confirmado",text)
        failed={"observation":{"status":"INTEGRATION_FAILURE","trigger_confirmed":False}}
        self.assertIn("incompleta",visual.observation_text(failed,1e-12))

    def test_confirmed_trigger_uses_saved_time_only(self):
        m={"observation":{"trigger_confirmed":True,"trigger_time_s":2.75e-12}}
        self.assertEqual(visual.observation_text(m,10e-12),"Gatillo confirmado: 2.75 ps")
        self.assertEqual(visual.time_label_ps(.5e-15),"0.0005 ps")

    def test_threshold_line_reads_runner_manifest_and_legacy_fallback(self):
        runner_spec=importlib.util.spec_from_file_location(
            "effective_photon_runner_for_visual",ROOT/"pipelines/05_effective_photon_run.py")
        runner=importlib.util.module_from_spec(runner_spec)
        runner_spec.loader.exec_module(runner)
        args=runner.parser().parse_args(["--output-dir","/home/jdiaz/scratch/visual_threshold_test"])
        manifest,trajectory,snapshots=self.fixtures()
        manifest["observation"]=runner.observation_manifest(None,args,"RUNNING")
        self.assertEqual(manifest["observation"]["threshold_voltage_V"],1e-4)
        for legacy in (False,True):
            with self.subTest(legacy=legacy):
                observation=manifest["observation"].copy()
                if legacy:
                    observation["trigger_voltage_V"]=observation.pop("threshold_voltage_V")
                current=dict(manifest,observation=observation)
                history,_=visual.prepare_view(current,trajectory,snapshots)
                figures=[]
                def capture(fig,*args,**kwargs):
                    figures.append(fig)
                try:
                    with patch.object(visual,"save_figure",side_effect=capture):
                        visual.trace_figures(current,history,Path("."),[],72)
                    lines=[line for line in figures[0].axes[0].lines
                        if line.get_label().startswith("Umbral:")]
                    self.assertEqual(len(lines),1)
                    np.testing.assert_array_equal(lines[0].get_ydata(),[.1,.1])
                finally:
                    for fig in figures:
                        visual.plt.close(fig)

    def test_missing_energy_reference_or_invalid_data_are_not_repaired(self):
        manifest,history,maps=self.fixtures()
        with self.assertRaisesRegex(ValueError,"retained_energy_J"):
            visual.prepare_view({},history,maps)
        history["Vout_V"][1]=np.nan
        with self.assertRaisesRegex(ValueError,"Vout_V"):
            visual.prepare_view(manifest,history,maps)
        manifest,history,maps=self.fixtures()
        maps["Te_K"][1,1]=-1.
        with self.assertRaisesRegex(ValueError,"negative"):
            visual.prepare_view(manifest,history,maps)


if __name__=="__main__":
    unittest.main()
