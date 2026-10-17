"""Independent summation fixtures, weak-depairing limit and domain guards."""
import json
from pathlib import Path
import numpy as np
import pytest
from pysnspd.experimental.finite_depairing_sums import finite_depairing_partials

_REFERENCE_PATH=Path(__file__).parent/"fixtures"/"finite_depairing_reference.json"
_REFERENCE=json.loads(_REFERENCE_PATH.read_text())
_FIELDS=("phi","rho","Q","t","rhorho","rhoQ","QQ","rhot","Qt","tt")
assert tuple(_REFERENCE["fields"])==_FIELDS


@pytest.mark.parametrize("case",_REFERENCE["cases"],
    ids=lambda c:f"rho{c['rho']:g}_theta{c['theta']:g}_Q{c['Q']:g}")
def test_ten_partials_match_independent_scipy_matsubara_reference(case):
    actual=np.array(finite_depairing_partials(case["rho"],case["theta"],case["Q"]),float)
    np.testing.assert_allclose(actual,case["values"],rtol=0.,atol=1e-10)


@pytest.mark.parametrize("rho,theta",[(0.,.1),(.04,.1),(.3,.3),(1.,1.),(4.,2.)])
def test_zero_depairing_recovers_dirty_limit_stiffness_and_temperature_derivative(rho,theta):
    values=finite_depairing_partials(rho,theta,0.)
    a=np.sqrt(rho)
    if a==0:
        current_coefficient=0.
        rho_current_coefficient=np.pi/(4*theta)
        temperature_current_coefficient=0.
    else:
        x=a/(2*theta);th=np.tanh(x);sech2=1/np.cosh(x)**2
        current_coefficient=np.pi*a/2*th
        rho_current_coefficient=np.pi/(4*a)*(th+x*sech2)
        temperature_current_coefficient=-np.pi*rho/(4*theta**2)*sech2
    np.testing.assert_allclose([values[2],values[5],values[8]],
        [current_coefficient,rho_current_coefficient,temperature_current_coefficient],
        rtol=2e-11,atol=1e-11)
    # The reference correction vanishes identically at Q=0.
    np.testing.assert_allclose(np.array(values)[[0,1,3,4,7,9]],0.,rtol=0.,atol=1e-12)


@pytest.mark.parametrize("rho,theta",[(.04,.1),(.3,.3),(1.,1.),(4.,2.)])
def test_weak_depairing_recovers_stiffness_temperature_curvature(rho,theta):
    a=np.sqrt(rho);x=a/(2*theta);th=np.tanh(x);sech2=1/np.cosh(x)**2
    kappa=np.pi/(2*a)*th
    kappa_t=-np.pi/(4*theta**2)*sech2
    kappa_tt=np.pi/(2*theta**3)*sech2*(1-x*th)
    expected=rho*np.array([kappa,kappa_t,kappa_tt])
    Q=1e-6
    values=np.array(finite_depairing_partials(rho,theta,Q))
    actual=values[[0,3,9]]/Q
    # The finite-Q correction has a physical O(Q) departure from this limit.
    np.testing.assert_allclose(actual,expected,rtol=8e-5,atol=2e-6)


@pytest.mark.parametrize("rho,theta,Q",[(0.,.05,1.),(.1,.104,.01),(1.,.104,1.),(4.,5.,1e6)])
def test_explicit_sum_refinement_preserves_all_ten_partials(rho,theta,Q):
    coarse=np.array(finite_depairing_partials(rho,theta,Q,min_terms=64))
    fine=np.array(finite_depairing_partials(rho,theta,Q,min_terms=128))
    np.testing.assert_allclose(coarse,fine,rtol=0.,atol=1e-10)


@pytest.mark.parametrize("rho,theta,Q",[
    (-1.,1.,0.),(1.,0.,0.),(1.,-1.,0.),(1.,1.,-1.),
    (np.nan,1.,0.),(1.,np.inf,0.),(1.,1.,np.inf)])
def test_invalid_physical_domain_is_rejected(rho,theta,Q):
    with pytest.raises(ValueError):
        finite_depairing_partials(rho,theta,Q)


@pytest.mark.parametrize("kwargs",[dict(min_terms=32),dict(tail_scale=8.)])
def test_declared_tail_resolution_is_not_silently_lowered(kwargs):
    with pytest.raises(ValueError):
        finite_depairing_partials(1.,1.,1.,**kwargs)
