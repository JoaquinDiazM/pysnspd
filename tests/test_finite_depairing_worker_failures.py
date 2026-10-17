"""Parallel summation must never return uninitialized or partial results."""
import numpy as np
import pytest
from pysnspd.experimental.finite_depairing_sums import finite_depairing_partials


def test_cold_mixed_batch_matches_scalar_values_after_root_roundoff_stagnation():
    rho=np.array([2.1,2.1,.7428571428571429,1.8285714285714285,3.7285714285714286])
    theta=np.array([.3,.05,.05,.2,.05])/8.65
    Q=np.array([0.,1.4084745762711863,.9050847457627118,1.156779661016949,1.9118644067796606])
    batch=np.array(finite_depairing_partials(rho,theta,Q))
    scalar=np.array([finite_depairing_partials(r,t,q)for r,t,q in zip(rho,theta,Q)]).T
    assert np.all(np.isfinite(batch))
    np.testing.assert_allclose(batch,scalar,rtol=0.,atol=1e-10)
    assert np.all(batch[2]>0.)
    expected=np.pi*np.sqrt(rho[0])/2
    assert abs(batch[2,0]-expected)<1e-11


def test_unresolved_worker_aborts_entire_batch_instead_of_returning_zero_blocks():
    # The second state requires more than the declared 100000 explicit modes.
    # This must be a numerical-resolution failure, never an accepted zero RHS.
    with pytest.raises(ArithmeticError,match="Matsubara|resolution"):
        finite_depairing_partials(np.array([2.1,2.1]),
            np.array([.3/8.65,1e-8]),np.array([0.,1.]))


def test_invalid_member_is_rejected_before_parallel_dispatch():
    with pytest.raises(ValueError):
        finite_depairing_partials(np.array([2.1,2.1]),
            np.array([.3/8.65,0.]),np.array([0.,1.]))


@pytest.mark.parametrize("kwargs",[dict(min_terms=64.5),dict(tail_scale=np.nan)])
def test_nonfinite_or_fractional_work_controls_are_rejected(kwargs):
    with pytest.raises(ValueError):
        finite_depairing_partials(2.1,.3/8.65,0.,**kwargs)
