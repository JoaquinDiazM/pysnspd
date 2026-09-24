"""The inherited local KWT Euler quadratic, evaluated in amplitude increments.

This is the same first-order update as TDGLSolver.solve_for_psi_squared.  It
does not change the continuous law or add a solver. With r=|psi|^2 and
a=dt*sqrt(1+gamma^2*r)*forcing/u, write

  psi_next = U*(psi+a-gamma^2*psi*delta_r/2),  U=exp(-i*mu*dt).
  A*delta_r^2-B*delta_r+C=0,
  A=gamma^4*r/4,
  B=1+gamma^2*(r+Re(conj(psi)*a)),
  C=2*Re(conj(psi)*a)+|a|^2.

The root is algebraically identical to the inherited smaller-amplitude root.
Computing delta_r directly avoids subtracting two large w and z*r terms at
equilibrium. No force, baseline drift or physical residual is subtracted.
The historical production implementation and previous campaign bridge remain
unchanged; only the experimental DC hold selects this arithmetic form.
"""
import numpy as np

from .heredado_kwt_bridge import InheritedKWTStep
from .thermal_weak_response import ThermalKWTNormal


def incremental_local_update(psi, forcing, *, gamma, u, dt, mu):
    psi = np.asarray(psi, complex)
    forcing = np.asarray(forcing, complex)
    mu = np.asarray(mu, float)
    if psi.ndim != 1 or forcing.shape != psi.shape or mu.shape != psi.shape:
        raise ValueError('One-dimensional psi, forcing and potential arrays of equal shape required')
    if any(np.any(~np.isfinite(a)) for a in (psi, forcing, mu)):
        raise ValueError('Finite local fields required')
    if any(not np.isfinite(v) or v <= 0 for v in (u, dt)) or not np.isfinite(gamma) or gamma < 0:
        raise ValueError('Positive u,dt and nonnegative finite gamma required')
    r = abs(psi)**2
    g2 = gamma**2
    a = dt/u*np.sqrt(1+g2*r)*forcing
    c = np.real(np.conj(psi)*a)
    A = .25*g2*g2*r
    B = 1+g2*(r+c)
    C = 2*c+abs(a)**2
    discriminant = B*B-4*A*C
    roundoff = 64*np.finfo(float).eps*np.maximum.reduce((B*B, abs(4*A*C), np.ones_like(B)))
    if np.any(discriminant < -roundoff):
        raise RuntimeError('Inherited KWT quadratic has no real admitted amplitude root')
    root = np.sqrt(np.maximum(discriminant, 0.))
    delta = np.zeros_like(r)
    stable_denominator = B+root
    normal = (B >= 0) & (stable_denominator != 0)
    delta[normal] = 2*C[normal]/stable_denominator[normal]
    # With B<0 the non-rationalized numerator has no cancellation; it also
    # preserves the inherited root for large, potentially rejected updates.
    negative = B < 0
    delta[negative] = (B[negative]-root[negative])/(2*A[negative])
    updated = np.exp(-1j*mu*dt)*(psi+a-.5*g2*psi*delta)
    amplitude_squared = r+delta
    if (np.any(~np.isfinite(updated)) or np.any(~np.isfinite(amplitude_squared))
            or np.any(amplitude_squared < -64*np.finfo(float).eps*np.maximum(1., r))):
        raise RuntimeError('Nonfinite or negative inherited KWT amplitude')
    return updated, amplitude_squared


def stable_inherited_kwt_step(model, gradient, current, *, dt_ps, delta0_over_kBTc):
    """Apply the original bridge's physical normalization to the same quadratic."""
    if not isinstance(model, ThermalKWTNormal):
        raise TypeError('A ThermalKWTNormal state is required')
    if any(not np.isfinite(v) or v <= 0 for v in (dt_ps, delta0_over_kBTc)):
        raise ValueError('Positive finite time step and Delta0/(kBTc) required')
    force = np.asarray(gradient, complex)
    if force.shape != model.d.shape or np.any(~np.isfinite(force)):
        raise ValueError('One finite integrated complex gradient per node required')
    potential = model.potential(current)
    ratio = np.pi/4
    tau0_ps = ratio*model.tD_ps
    dt = float(dt_ps/tau0_ps)
    gamma = float(np.sqrt(model.kappa)*delta0_over_kBTc)
    free = model.free
    gap = model.d.copy()
    if np.any(free):
        updated, _ = incremental_local_update(model.d[free]/delta0_over_kBTc,
            -ratio*force[free]/(model.denominator[free]*delta0_over_kBTc),
            gamma=gamma, u=1., dt=dt, mu=ratio*potential[free]/2)
        gap[free] = delta0_over_kBTc*updated
    return InheritedKWTStep(gap, potential, dt, float(tau0_ps), gamma)
