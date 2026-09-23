"""Analytic interpolation oracle; physical catalogue assessed separately."""
from types import SimpleNamespace
import numpy as np
import pytest
from pysnspd.experimental.local_energy_patch import LocalEnergyPatch


class PolynomialSource:
    count_nodes=np.array([.1,.2,.5,1.,2.])
    count_weights=np.ones(5)/5
    eta=1e-5
    vacuum=SimpleNamespace(delta0_J=1.,N0_per_J_m3=2.,D_m2_s=3.,
                           evaluate=lambda a,g:(a*a+g*g,2*a,2*g))
    def energy_kernel(self,a,g):
        x=self.count_nodes
        return x+a*a+.2*a*g+.3*g*g, np.full_like(x,2*a+.2*g),np.full_like(x,.2*a+.6*g)


def test_patch_matches_polynomial_energy_and_both_derivatives(tmp_path):
    source=PolynomialSource();seen=[]
    patch=LocalEnergyPatch.build(source,(.9,1.1),(.001,.02),4,on_sample=lambda i,n:seen.append((i,n)))
    assert seen[-1]==(25,25)
    for a,g in [(1.013,.0117),(.9,.001),(1.1,.02)]:
        np.testing.assert_allclose(patch.energy_kernel(a,g),source.energy_kernel(a,g),rtol=2e-10,atol=1e-10)
    path=tmp_path/'patch.npz';patch.save(path);loaded=LocalEnergyPatch.load(path,source)
    np.testing.assert_array_equal(loaded.coefficients,patch.coefficients)
    with pytest.raises(FileExistsError):patch.save(path)
    source.eta=2e-5
    with pytest.raises(ValueError):LocalEnergyPatch.load(path,source)


def test_one_energy_produces_force_at_fixed_population():
    patch=LocalEnergyPatch.build(PolynomialSource(),(.9,1.1),(.001,.02),4)
    p=np.array([.2,.17,.12,.04,.01]);a,g=1.013,.0076;h=1e-6
    u,ua,ug=patch.evaluate(a,g,p)
    measured=[(patch.evaluate(a+h,g,p)[0]-patch.evaluate(a-h,g,p)[0])/(2*h),
              (patch.evaluate(a,g+h,p)[0]-patch.evaluate(a,g-h,p)[0])/(2*h)]
    np.testing.assert_allclose(measured,[ua,ug],rtol=1e-8,atol=1e-9)


@pytest.mark.parametrize('a,g',[(.899,.01),(1.101,.01),(1.,0.),(1.,.02001),(np.nan,.01),(1.,np.inf)])
def test_no_extrapolation(a,g):
    patch=LocalEnergyPatch.build(PolynomialSource(),(.9,1.1),(.001,.02),2)
    with pytest.raises(ValueError):patch.energy_kernel(a,g)


def test_population_and_spectral_guards():
    patch=LocalEnergyPatch.build(PolynomialSource(),(.9,1.1),(.001,.02),2)
    for p in [np.full(5,-1e-10),np.full(5,1+1e-10),np.full(5,1j),np.zeros(4)]:
        with pytest.raises(ValueError):patch.evaluate(1.,.01,p)
    bad=LocalEnergyPatch(PolynomialSource(),(.9,1.1),(.001,.02),-patch.coefficients)
    with pytest.raises(FloatingPointError):bad.energy_kernel(1.,.01)
