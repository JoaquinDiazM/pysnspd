"""Leading adiabatic charge susceptibility at fixed gap and fixed spectrum.

The phasor convention is exp(-i omega t), nu=omega*hbar/(2*kB*Tc).
This diagnostic does not introduce a potential, capacitance, collision rate,
moving-gap energy law or exact finite-frequency Usadel spectral response.
"""
from __future__ import annotations
from dataclasses import dataclass
import warnings
import numpy as np
from scipy.sparse import diags
from scipy.sparse.linalg import MatrixRankWarning, spsolve
from .thermal_spatial_usadel import _fixed_indices


@dataclass(frozen=True)
class ChargeResponse:
    directions: np.ndarray
    storage: np.ndarray
    free_nodes: np.ndarray
    equation_residual: np.ndarray


def solve(operator, hL, fixed_nodes, nu):
    """Solve (-L_TT-i nu M) hT=L_TL hL with homogeneous contacts.

    M=m Re(g) is the DOS storage of the fixed-spectrum adiabatic kinetic
    closure. Numerical eta appears only through the supplied spectrum.
    """
    n=operator.graph.n_nodes
    hL=np.asarray(hL,complex)
    if hL.shape != (n,) or np.any(~np.isfinite(hL)) or not np.isfinite(nu) or nu<0:
        raise ValueError('Finite phasor per node and nonnegative nu required')
    fixed=_fixed_indices(fixed_nodes,n)
    if not len(fixed):
        raise ValueError('Explicit valid homogeneous contacts required')
    if np.any(hL[fixed]!=0):
        raise ValueError('Energy-mode direction must vanish at contacts')
    free=np.ones(n,bool);free[fixed]=False
    storage=operator.graph.area_weights*operator.R[:,0,0].real
    if np.any(storage<0) or np.any(~np.isfinite(storage)):
        raise ValueError('Nonnegative finite DOS storage required; no clipping applied')
    charge=np.arange(1,2*n,2)
    matrix=-operator.matrix[charge][:,charge]-1j*nu*diags(storage)
    rhs=operator.matrix[charge][:,charge-1]@hL
    hT=np.zeros(n,complex)
    with warnings.catch_warnings():
        warnings.simplefilter('error',MatrixRankWarning)
        hT[free]=spsolve(matrix[free][:,free],rhs[free])
    if np.any(~np.isfinite(hT)):
        raise RuntimeError('Nonfinite response; no leakage or regularizer added')
    return ChargeResponse(np.column_stack((hL,hT)),storage,free,matrix@hT-rhs)


def observe(operator, directions, *, eta=0.):
    """Apply real kinetic maps to two quadratures, keeping Nambu axes distinct.

    force_cartesian[n,0/1] are independent complex phasors of the real and
    imaginary gap-force coordinates. They must never be collapsed into one
    complex number: its imaginary unit would mix Nambu and time quadratures.
    """
    directions=np.asarray(directions,complex)
    first=operator.evaluate(directions.real,eta=eta)
    second=operator.evaluate(directions.imag,eta=eta)
    force=np.column_stack((first.gap_force_increment.real+1j*second.gap_force_increment.real,
                           first.gap_force_increment.imag+1j*second.gap_force_increment.imag))
    return dict(force_cartesian=force,
        current=first.charge_current_increment+1j*second.charge_current_increment,
        energy_flux=first.edge_flux[:,0]+1j*second.edge_flux[:,0],
        charge_residual=first.residual[:,1]+1j*second.residual[:,1],
        energy_residual=first.residual[:,0]+1j*second.residual[:,0],
        artificial_eta_leakage=first.artificial_eta_leakage+1j*second.artificial_eta_leakage)
