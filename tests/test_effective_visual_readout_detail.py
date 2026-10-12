"""Readout visibility without rescaling saved values or inventing a trigger."""
import importlib.util
import os
from pathlib import Path
import unittest
from unittest.mock import patch
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
path=Path(os.environ.get("PYSNSPD_VISUAL_SOURCE", ROOT/"plot_pipelines/05_plot_effective_photon_run.py"))
spec=importlib.util.spec_from_file_location("effective_visual_readout_detail",path)
visual=importlib.util.module_from_spec(spec); spec.loader.exec_module(visual)

class ReadoutDetailTests(unittest.TestCase):
    def make(self, pulse, device, times=(0.,100.,200.), obs=None):
        fig,ax=visual.plt.subplots()
        self.addCleanup(lambda:visual.plt.close(fig))
        history={"time_ps":np.array(times),"Vout_relative_mV":np.array(pulse),
                 "Vdev_mV":np.array(device)}
        manifest={"observation":dict(threshold_voltage_V=1e-4, **(obs or {}))}
        return ax,manifest,history

    def test_331_nanovolt_pulse_is_point331_microvolt_with_saved_window(self):
        ax,m,h=self.make([0.,.000331,.0002],[0.,.037,-.001])
        before={k:v.copy() for k,v in h.items()}
        detail=visual.add_readout_detail(ax,m,h)
        np.testing.assert_array_equal(detail.lines[0].get_ydata(),[0.,.331,.2])
        self.assertIn(r"$\mu$V",detail.get_ylabel())
        self.assertEqual(detail.get_xlim(),(0.,200.))
        self.assertEqual(len(detail.lines),1)  # No out-of-scale trigger or fabricated crossing.
        for k in h:np.testing.assert_array_equal(h[k],before[k])

    def test_millivolt_pulse_and_real_confirmed_times_are_kept(self):
        ax,m,h=self.make([0.,.1,.149],[0.,19.,22.],times=(0.,31.9,36.94),
            obs={"trigger_confirmed":True,"first_crossing_time_s":31.9e-12,
                 "trigger_time_s":32.4e-12})
        detail=visual.add_readout_detail(ax,m,h)
        np.testing.assert_array_equal(detail.lines[0].get_ydata(),[0.,.1,.149])
        self.assertTrue(detail.get_ylabel().endswith("(mV)"))
        self.assertAlmostEqual(detail.lines[1].get_xdata()[0],31.9)
        self.assertAlmostEqual(detail.lines[2].get_xdata()[0],32.4)

    def test_readable_signal_does_not_add_unnecessary_inset(self):
        ax,m,h=self.make([0.,.1,.15],[0.,.08,.1])
        self.assertIsNone(visual.add_readout_detail(ax,m,h))

    def test_trigger_alone_can_flatten_a_small_readout(self):
        ax,m,h=self.make([0.,.001,.002],[0.,.001,.002])
        self.assertIsNotNone(visual.add_readout_detail(ax,m,h))

if __name__=="__main__":unittest.main()
