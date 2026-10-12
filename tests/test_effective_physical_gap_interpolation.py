import importlib.util
from pathlib import Path
import unittest
import numpy as np

HERE=Path(__file__).resolve().parent
SOURCE=HERE/"physical_gap_interpolation.py"
if not SOURCE.exists():SOURCE=HERE.parent/"sandbox/stage5_effective_validation/physical_gap_interpolation.py"
spec = importlib.util.spec_from_file_location("physical_gap_interpolation",SOURCE)
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)


class PhysicalMagnitudeTimeInterpolation(unittest.TestCase):
    def setUp(self):
        self.sampler = dict(nodes=np.array([[0, 1, 2], [0, 1, 2]]),
            barycentric=np.array([[.2, .3, .5], [.6, .1, .3]]))
        self.snap = dict(time_s=np.array([0., 1., 2.]),
            delta_J=np.array([[1+1j, 2-.3j, -.2+.7j],
                              [.6+.3j, 1.7+1.2j, -.4+.2j],
                              [.9-.1j, 1.4-.2j, .7+.5j]]))

    def test_arbitrary_global_phase_per_frame(self):
        rotated = dict(self.snap,
            delta_J=self.snap["delta_J"] * np.exp(1j*np.array([.7, 2.6, -1.9]))[:, None])
        for t in (0., .17, .5, 1., 1.6, 2.):
            np.testing.assert_allclose(m.spatial_magnitude_at(self.snap,self.sampler,t),
                m.spatial_magnitude_at(rotated,self.sampler,t),rtol=2e-15,atol=2e-15)
        # The replaced complex-chord operation is demonstrably phase dependent.
        def old(s):
            frame=.5*s["delta_J"][0]+.5*s["delta_J"][1]
            return abs(np.sum(frame[self.sampler["nodes"]]*self.sampler["barycentric"],axis=1))
        self.assertGreater(float(max(abs(old(self.snap)-old(rotated)))),.1)

    def test_exact_saved_frames(self):
        for j,t in enumerate(self.snap["time_s"]):
            expected=abs(np.sum(self.snap["delta_J"][j][self.sampler["nodes"]]
                        * self.sampler["barycentric"],axis=1))
            np.testing.assert_array_equal(m.spatial_magnitude_at(self.snap,self.sampler,t),expected)

    def test_known_scalar_magnitude_interpolation(self):
        snap=dict(time_s=np.array([2.,4.]),
            delta_J=np.array([[2,2,2],[6j,6j,6j]],complex))
        np.testing.assert_allclose(m.spatial_magnitude_at(snap,self.sampler,2.5),[3.,3.])
        np.testing.assert_allclose(m.spatial_magnitude_at(snap,self.sampler,3.),[4.,4.])

    def test_complex_spatial_P1_precedes_magnitude(self):
        phase=np.exp(2j*np.pi*np.arange(3)/3)
        snap=dict(time_s=np.array([0.,1.]),delta_J=np.array([phase,2*phase]))
        sampler=dict(nodes=np.array([[0,1,2]]),barycentric=np.full((1,3),1/3))
        self.assertLess(float(m.spatial_magnitude_at(snap,sampler,.5)[0]),1e-15)
        self.assertAlmostEqual(float(np.sum(abs(phase)*sampler["barycentric"][0])),1.)

    def test_no_extrapolation_or_nonfinite_data(self):
        for t in (-.1,2.1,float("nan")):
            with self.subTest(t=t),self.assertRaises(ValueError):
                m.spatial_magnitude_at(self.snap,self.sampler,t)
        bad=dict(self.snap,delta_J=self.snap["delta_J"].copy());bad["delta_J"][0,0]=np.nan
        with self.assertRaises(ValueError):m.spatial_magnitude_at(bad,self.sampler,.5)

    def test_field_norm_invariance(self):
        rotated=dict(self.snap,delta_J=self.snap["delta_J"]*np.exp(1j*np.array([2.,-.9,.4]))[:,None])
        z=m.field_L2_saved_magnitude(self.snap,rotated,self.sampler,self.sampler,[0.,.3,1.7,2.],1.)
        self.assertLess(float(z.max()),1e-15)

if __name__=="__main__":unittest.main()
