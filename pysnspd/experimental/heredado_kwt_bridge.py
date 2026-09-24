"""Use the thesis/pyTDGL local Euler update for the admitted thermal KWT law.

The graph gap d is Delta/(kB*Tc), its gradient is area-integrated, and dt_ps
is a physical time in picoseconds. The fixed material normalization
delta0_over_kBTc is Delta0/(kB*Tc), not the local or equilibrium thermal gap.
Only the local algebraic update is reused. Spectral current and the normal
potential remain those supplied by ThermalKWTNormal; no GL current is added.
"""
from dataclasses import dataclass

import numpy as np
from scipy.sparse import csr_matrix

from pysnspd.solver.core import TDGLSolver
from .thermal_weak_response import ThermalKWTNormal


@dataclass(frozen=True)
class InheritedKWTStep:
    gap: np.ndarray
    potential_v: np.ndarray
    backend_dt: float
    backend_tau0_ps: float
    backend_gamma: float


def inherited_kwt_step(model: ThermalKWTNormal, gradient, current, *,
                       dt_ps: float, delta0_over_kBTc: float) -> InheritedKWTStep:
    """Advance one first-order local KWT step with the same continuous law.

    The model supplies the state, graph, fixed contacts and physical KWT
    coefficients. The caller supplies the integrated gap gradient and spectral
    edge current at that same state. The potential is solved by the model's
    existing normal-current block. No temporal error estimator is implied.

    With r=tau0/tD=pi/4 and d*=Delta0/(kB*Tc), the inherited call uses
      psi=d/d*, gamma**2=kappa*d*², mu=r*v/2, dt=dt_ps/(r*tD_ps),
      forcing=-r*gradient/(model.denominator*d*), u=1.
    model.denominator contains the node area times the KWT temperature factor.
    Thus forcing is before mobility inversion and excludes the gauge velocity.
    Fixed contacts are not passed to the local solve and remain bitwise fixed.
    """
    if not isinstance(model, ThermalKWTNormal):
        raise TypeError('A ThermalKWTNormal state is required')
    if (not np.isfinite(dt_ps) or dt_ps <= 0 or
            not np.isfinite(delta0_over_kBTc) or delta0_over_kBTc <= 0):
        raise ValueError('Positive finite dt_ps and Delta0/(kB*Tc) required')
    force = np.asarray(gradient, dtype=complex)
    if force.shape != model.d.shape or np.any(~np.isfinite(force)):
        raise ValueError('One finite integrated complex gradient per node required')
    potential = model.potential(current)
    ratio = np.pi / 4
    tau0_ps = ratio * model.tD_ps
    dt = float(dt_ps / tau0_ps)
    gamma = float(np.sqrt(model.kappa) * delta0_over_kBTc)
    free = model.free
    psi = model.d[free] / delta0_over_kBTc
    n = len(psi)
    gap = model.d.copy()
    if n:
        result = TDGLSolver.solve_for_psi_squared(
            psi=psi,
            abs_sq_psi=np.abs(psi)**2,
            mu=ratio * potential[free] / 2,
            epsilon=np.zeros(n),
            gamma=gamma,
            u=1.,
            dt=dt,
            psi_laplacian=csr_matrix((n, n)),
            forcing_dimensionless=-ratio * force[free] /
                (model.denominator[free] * delta0_over_kBTc),
        )
        if result is None:
            raise RuntimeError('Inherited local KWT update rejected this time step')
        gap[free] = delta0_over_kBTc * result[0]
    return InheritedKWTStep(gap, potential, dt, float(tau0_ps), gamma)
