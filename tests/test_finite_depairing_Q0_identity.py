"""The zero-depairing reference is an exact algebraic identity."""
import numpy as np
from pysnspd.experimental.finite_depairing_sums import finite_depairing_partials

def test_zero_depairing_rho_temperature_partials_are_exactly_zero_even_when_cold():
    rho=np.array([.1,1.,4.,4.,4.])
    theta=np.array([.01,.05,.01,.05,1.])
    values=np.array(finite_depairing_partials(rho,theta,np.zeros_like(rho)))
    np.testing.assert_array_equal(values[[0,1,3,4,7,9]],np.zeros((6,len(rho))))
    # Q derivatives must remain present: the reference has finite rigidity.
    expected=np.pi*np.sqrt(rho)/2*np.tanh(np.sqrt(rho)/(2*theta))
    np.testing.assert_allclose(values[2],expected,rtol=2e-12,atol=2e-12)
    assert np.all(values[5]>0.)
    assert np.all(values[6]<0.)
