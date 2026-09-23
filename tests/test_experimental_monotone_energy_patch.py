import numpy as np
import pytest
from test_experimental_local_energy_patch import PolynomialSource
from pysnspd.experimental.local_energy_patch import LocalEnergyPatch
from pysnspd.experimental.monotone_energy_patch import MonotoneEnergyPatch


def test_exact_same_interpolant_derivatives_and_order(tmp_path):
    source=PolynomialSource();patch=MonotoneEnergyPatch.build(source,(.9,1.1),(.001,.02),6)
    a,g=1.013,.0076;h=1e-6
    e,ea,eg=patch.energy_kernel(a,g)
    assert np.all(e>0) and np.all(np.diff(e)>0)
    np.testing.assert_allclose([e,ea,eg],source.energy_kernel(a,g),rtol=2e-7,atol=2e-7)
    for analytic,plus,minus in [(ea,(a+h,g),(a-h,g)),(eg,(a,g+h),(a,g-h))]:
        np.testing.assert_allclose(analytic,(patch.energy_kernel(*plus)[0]-patch.energy_kernel(*minus)[0])/(2*h),rtol=1e-7,atol=1e-8)
    path=tmp_path/'no_suffix';patch.save(path)
    with pytest.raises(FileExistsError):patch.save(path)
    loaded=MonotoneEnergyPatch.load(path,source)
    np.testing.assert_array_equal(loaded.energy_kernel(a,g),patch.energy_kernel(a,g))
    with pytest.raises(ValueError):patch.energy_kernel(.899,g)


def test_representations_cannot_be_confused(tmp_path):
    source=PolynomialSource();path=tmp_path/'raw.npz'
    LocalEnergyPatch.build(source,(.9,1.1),(.001,.02),2).save(path)
    with pytest.raises(ValueError):MonotoneEnergyPatch.load(path,source)


def test_underresolved_increment_rejected_not_rounded_up():
    source=PolynomialSource();coeff=np.zeros((3,3,5));coeff[0,0]=[0.,-1000.,0.,0.,0.]
    patch=MonotoneEnergyPatch(source,(.9,1.1),(.001,.02),coeff)
    with pytest.raises(FloatingPointError):patch.energy_kernel(1.,.01)
