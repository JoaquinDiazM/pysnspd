"""Physical clock progress and guarded thermal midpoint convergence."""
import numpy as np
from test_effective_midpoint_physical_norm import physical_system
from test_thermal_stiffness_coupled_candidate import small_system


def test_small_nonzero_rhs_is_not_accepted_as_an_unchanged_predictor():
    m,s=physical_system()
    i=m.free[len(m.free)//2]
    s.p[i]=m.kinetics.fermi_dirac(s.delta_bar,8.)[i]
    rhs=m.evaluate(s)
    h=1e-10
    result,receipt=m.midpoint_step(s,h,tolerance=2e-7,max_iterations=40)
    assert receipt["iterations"]>=1
    assert np.max(abs(result.delta_bar-s.delta_bar))>0
    assert np.max(abs(result.p-s.p))>0
    np.testing.assert_allclose((result.delta_bar-s.delta_bar)/h,
        rhs.delta_per_ps,rtol=2e-4,atol=3e-6)
    np.testing.assert_allclose((result.p-s.p)/h,
        rhs.dp_per_ps,rtol=2e-4,atol=3e-6)


def test_thermal_solve_tolerance_scales_with_physical_step():
    m,s=small_system()
    _,receipt=m.midpoint_step(s,1e-4,tolerance=2e-7,max_iterations=60)
    assert receipt["iterations"]>=1
    assert receipt["nonlinear_actual_tolerance"]==1e-9
    assert receipt["defect"]<receipt["nonlinear_actual_tolerance"]
    assert receipt["thermal_induced_gap_defect"]<receipt["nonlinear_actual_tolerance"]


def test_internal_capacity_line_search_retains_positive_domain():
    m,s=small_system()
    zeros=(np.zeros(m.ngap),np.zeros_like(s.p),np.zeros_like(s.n),np.zeros(3))
    cu=np.full(m.ngap,10.);cp=np.zeros_like(s.p);cn=np.zeros_like(s.n);cc=np.zeros(3)
    fraction=m.coupled_iteration_damping(s,*zeros,cu,cp,cn,cc)
    assert 0<fraction<=1
    for scale in (.5,1.):
        trial=m.shifted(s,fraction*cu,cp,cn,cc,scale)
        ce=m.functional.thermal_data(trial.delta_bar,trial.p)["capacity_effective_J_K"]
        assert np.all(ce>0)

def test_full_thermal_gap_guard_participates_in_iteration(monkeypatch):
    m,s=small_system()
    original=m.thermal_population_gap_defect
    calls=[]
    def guard(state,rhs,residual,h):
        calls.append(1)
        if len(calls)<=2:return 1.
        return original(state,rhs,residual,h)
    monkeypatch.setattr(m,"thermal_population_gap_defect",guard)
    _,receipt=m.midpoint_step(s,1e-4,tolerance=2e-7,max_iterations=60)
    assert len(calls)>=3
    assert receipt["iterations"]>=2
