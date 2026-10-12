"""Numerical saved-field visualization checks; no figures or simulation."""
import unittest
import numpy as np
from matplotlib.tri import Triangulation
from pysnspd.plotting.complex_gap import complex_gap_at, complex_gap_raster, draw_complex_gap
from sandbox.stage5_effective_validation.physical_gap_interpolation import spatial_magnitude_at


class ComplexGapRasterTests(unittest.TestCase):
    def setUp(self):
        self.tri=Triangulation([0.,2.,0.],[0.,0.,2.],[[0,1,2]])
        self.delta=np.exp(2j*np.pi*np.arange(3)/3)

    def test_interior_zero_with_nonzero_nodal_gap(self):
        value=complex_gap_at(self.tri,self.delta,2/3,2/3)
        self.assertLess(float(value),5e-16)
        self.assertGreater(float(abs(self.delta).min()),.999)
        self.assertGreater(float(np.mean(abs(self.delta))),.999)

    def test_numeric_agreement_with_adopted_complex_P1_sampler(self):
        bary=np.array([[1/3,1/3,1/3],[.1,.7,.2],[.6,.1,.3]])
        snap=dict(time_s=np.array([0.]),delta_J=self.delta[None,:])
        sampler=dict(nodes=np.tile([0,1,2],(len(bary),1)),barycentric=bary)
        expected=spatial_magnitude_at(snap,sampler,0.)
        actual=complex_gap_at(self.tri,self.delta,2*bary[:,1],2*bary[:,2])
        np.testing.assert_allclose(actual,expected,rtol=2e-15,atol=5e-16)

    def test_constant_complex_field_is_constant_inside_and_masked_outside(self):
        delta=np.full(3,3+4j)
        actual=complex_gap_at(self.tri,delta,[0.,.2,2/3,2.],[0.,.2,2/3,2.])
        np.testing.assert_allclose(actual[:3],[5.,5.,5.],rtol=0,atol=1e-15)
        self.assertTrue(np.ma.getmaskarray(actual)[3])

    def test_affine_complex_field_reconstructed_before_magnitude(self):
        def affine(x,y):
            return (1+.3j)+(2-.7j)*x+(-.5+.2j)*y
        delta=affine(self.tri.x,self.tri.y)
        x=np.array([.1,.3,.9]);y=np.array([.1,.6,.4])
        actual=complex_gap_at(self.tri,delta,x,y)
        np.testing.assert_allclose(actual,abs(affine(x,y)),rtol=2e-15,atol=2e-15)

    def test_phase_rotation_and_physical_energy_units(self):
        x=np.array([.2,2/3,.6]);y=np.array([.2,2/3,.3])
        a=complex_gap_at(self.tri,self.delta,x,y)
        b=complex_gap_at(self.tri,self.delta*np.exp(1.7j)*1.602176634e-22,x,y)/1.602176634e-22
        np.testing.assert_allclose(a,b,rtol=2e-15,atol=5e-16)

    def test_mask_and_raster_preserve_domain_and_interior_zero(self):
        x,y,gap=complex_gap_raster(self.tri,self.delta,samples=7)
        self.assertEqual(gap.shape,(7,7))
        self.assertTrue(np.ma.getmaskarray(gap)[-1,-1])
        self.assertLess(float(gap[2,2]),5e-16)
        self.assertEqual((x[0],x[-1],y[0],y[-1]),(0.,2.,0.,2.))

    def test_renderer_passes_complex_sampled_values_without_a_figure(self):
        class Ax:
            def pcolormesh(self,x,y,values,**kwargs):
                self.received=(x,y,values,kwargs)
                return "sentinel"
        ax=Ax();unit=1.602176634e-22
        result=draw_complex_gap(ax,self.tri,self.delta*unit,unit_J=unit,
            vmin=0.,vmax=1.,samples=7)
        self.assertEqual(result,"sentinel")
        self.assertLess(float(ax.received[2][2,2]),5e-16)
        self.assertEqual(ax.received[3]["shading"],"nearest")

    def test_invalid_fields_bounds_and_units_rejected(self):
        for delta in ([1.,2.],[1.,np.nan,1.]):
            with self.subTest(delta=delta),self.assertRaises(ValueError):
                complex_gap_at(self.tri,delta,0.,0.)
        for kwargs in ({"samples":2},{"samples":True},{"bounds":(0,0,0,1)}):
            with self.subTest(kwargs=kwargs),self.assertRaises(ValueError):
                complex_gap_raster(self.tri,self.delta,**kwargs)
        with self.assertRaises(ValueError):
            draw_complex_gap(None,self.tri,self.delta,unit_J=0,vmin=0,vmax=1)

if __name__=="__main__":unittest.main()
