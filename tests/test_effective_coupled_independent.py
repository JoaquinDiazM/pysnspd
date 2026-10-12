"""Independent residual audit of the coupled phase/potential/CM Newton block."""
import numpy as np
import pytest
from test_effective_coupled import system


def test_joint_midpoint_block_matches_independent_finite_residual():
    """Differentiate the defining residual, freezing only KWT mobility."""
    m,s=system()
    s.delta_bar[m.free[2]]*=.95*np.exp(.04j)
    s.circuit_SI*=np.array([1.01,.98,1.])
    s.circuit_SI[2]=2e-5
    r=m.evaluate(s)
    step=.025
    factor,nphi=m.midpoint_preconditioner(s,r,step)
    mob=m.mobility_blocks(s.delta_bar,r.temperature_K)
    ng=m.ngap
    ports=m.ports
    P=ports.prolongation[:,ports.unknown]
    phi0=np.r_[r.electrical.phi_V[m.free],
        r.electrical.phi_V[ports.left[0]],r.electrical.phi_V[ports.right[0]]][ports.unknown]
    B=ports.reduced_incidence[ports.unknown]
    lap=ports.normal_laplacian[ports.unknown][:,ports.unknown]
    injection=ports.group_injection_per_A[ports.unknown]

    def residual(du,dphi,dc):
        mid=m.shifted(s,du,np.zeros_like(s.p),np.zeros_like(s.n),dc,.5)
        phi=P@(phi0+dphi)
        ev=m.functional.evaluate(mid.delta_bar,mid.p)
        force=np.column_stack((ev.gap_force_J_per_gapbar.real,
                               ev.gap_force_J_per_gapbar.imag)).ravel()
        material=(-mob@force).reshape(-1,2)
        dz=material[:,0]+1j*material[:,1]-1j*m.kappa_per_ps_V*phi*mid.delta_bar
        rates=np.r_[np.column_stack((dz[m.free].real,dz[m.free].imag)).ravel(),
            -m.kappa_per_ps_V*phi[ports.left[0]],
            -m.kappa_per_ps_V*phi[ports.right[0]]]
        gap=du-step*rates
        potential=B@ev.edge_current_A+lap@(phi0+dphi)-injection*mid.circuit_SI[1]
        vdev=phi[ports.left[0]]-phi[ports.right[0]]
        circuit=dc-step*1e-12*m.circuit.rhs(mid.circuit_SI,vdev)/m.circuit_scales
        return np.r_[gap,potential,circuit]

    rng=np.random.default_rng(20261011)
    du=rng.normal(size=ng)*.1
    dphi=rng.normal(size=nphi)*1e-4
    dc=rng.normal(size=3)*.1
    eps=1e-5
    finite=(residual(eps*du,eps*dphi,eps*dc)-
            residual(-eps*du,-eps*dphi,-eps*dc))/(2*eps)
    # Solving the tested assembled derivative must recover every unknown,
    # including phi_mid (not the endpoint potential) and both cut phases.
    recovered=factor.solve(finite)
    expected=np.r_[du,dphi,dc]
    np.testing.assert_allclose(recovered,expected,rtol=2e-6,atol=2e-9)


def test_hot_cut_reservoir_heat_sign_and_zero_material_rates():
    m,s=system()
    s.p[m.cut_nodes]=m.kinetics.fermi_dirac(s.delta_bar,8.)[m.cut_nodes]
    s.n[m.cut_nodes]=m.kinetics.bose(8.)
    r=m.evaluate(s)
    # A prescribed hot reservoir injects energy into the colder interior.
    assert r.reservoir_power_out_W<0.
    np.testing.assert_array_equal(r.material_delta_per_ps[m.cut_nodes],0.)
    np.testing.assert_array_equal(r.dp_per_ps[m.cut_nodes],0.)
    np.testing.assert_array_equal(r.dn_per_ps[m.cut_nodes],0.)
    assert abs(r.film_power_residual_W)<1e-22
    assert abs(r.instant_power_residual_W)<1e-22
