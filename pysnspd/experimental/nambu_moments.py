"""Finite SI moments of the same full Nambu Keldysh field.

The caller supplies the causal spectrum, both distributions, first-gradient
correction and fixed-energy quadrature.  This module performs no evolution,
collision inference, eta regularization, calibration or residual cancellation.
It makes the charge, pair and current moments available to a joint closure.
The energy moment of K1 is reported separately from the existing thermal
internal-energy observer: adopting it requires the same approximation order,
normal background and charge/gauge energy contract, not merely summing ledgers.
"""
from __future__ import annotations
from dataclasses import dataclass
import numpy as np


def _real(value, name):
    raw=np.asarray(value)
    if np.iscomplexobj(raw) or raw.dtype.kind=='b' or np.any(~np.isfinite(raw)):
        raise ValueError(name+' must be finite real')
    return np.asarray(raw,float)


def _positive_scalar(value,name):
    raw=_real(value,name)
    if raw.ndim or raw<=0:raise ValueError(name+' must be a positive scalar')
    return float(raw)


def _matrix(value,shape,name):
    raw=np.asarray(value,complex)
    if raw.shape!=shape or np.any(~np.isfinite(raw)):
        raise ValueError(name+' must have finite (energy,node-or-edge,2,2) entries')
    return raw


@dataclass(frozen=True)
class NambuMomentResult:
    neutrality_residual_V: np.ndarray
    charge_moment_zero_J: np.ndarray
    charge_moment_first_J: np.ndarray
    gap_force_J_per_gapbar: np.ndarray
    gap_force_increment_zero_J_per_gapbar: np.ndarray
    gap_force_increment_first_J_per_gapbar: np.ndarray
    edge_current_A: np.ndarray
    edge_current_zero_A: np.ndarray
    edge_current_first_A: np.ndarray
    energy_increment_zero_J: float
    energy_increment_first_J: float
    reference_contract: str
    scope: str='finite moments only; no closed detector or material admission'


@dataclass(frozen=True)
class StableCurrentCompletion:
    edge_current_A: np.ndarray
    thermal_current_A: np.ndarray
    frozen_increment_A: np.ndarray
    first_increment_A: np.ndarray
    direct_spectral_current_A: np.ndarray
    direct_thermal_current_A: np.ndarray
    thermal_quadrature_discrepancy_A: np.ndarray
    reference_contract: str


class NambuMoments:
    """Positive-energy moments in the repository's electron-hole convention.

    N0 is the single-spin normal DOS [1/(J m^3)], Eu=kB Tc [J], volume [m^3],
    I0=Eu/(2e R_sheet) [A].  K has the normalization R H-H A and retained
    first-gradient correction. All energy measures are physical joules.
    The reference is the *frozen* thermal K0 of the same instantaneous gap,
    spectrum and gauge; it has no K1. Its off-shell free-energy gradient must
    be supplied independently. K1 therefore enters the charge and pair moment
    once. No Poisson ohmic current is added to the returned full current.
    """
    def __init__(self,energy_J,weights_J,volume_m3,N0_per_J_m3,energy_unit_J,current_unit_A):
        self.energy=_real(energy_J,'energy_J').copy()
        self.weights=_real(weights_J,'weights_J').copy()
        self.volume=_real(volume_m3,'volume_m3').copy()
        self.N0=_positive_scalar(N0_per_J_m3,'N0_per_J_m3')
        self.Eu=_positive_scalar(energy_unit_J,'energy_unit_J')
        self.I0=_positive_scalar(current_unit_A,'current_unit_A')
        if (self.energy.ndim!=1 or not self.energy.size or self.weights.shape!=self.energy.shape
            or np.any(self.energy<0) or np.any(np.diff(self.energy)<=0) or np.any(self.weights<=0)
            or self.volume.ndim!=1 or not self.volume.size or np.any(self.volume<=0)):
            raise ValueError('Ordered nonnegative energies and positive weights/volumes required')
        self.shape=(len(self.energy),len(self.volume),2,2)
        for array in (self.energy,self.weights,self.volume):array.setflags(write=False)

    def evaluate(self,*,K_zero,K_first,K_thermal_zero,phi_V,
                 thermal_gradient_J_per_gapbar,edge_flux_zero,edge_flux_first,reference_contract):
        """Return the retained moments, without solving neutrality or gap.

        Neutrality: e phi + integral Re Tr K/4 dE =0.
        Pair increment: i N0 Eu V integral [K-Kthermal]01 dE.
        Current: 2 I0 integral Pi_T(J_K) dE/Eu.
        Neutral energy increment: -N0 V/2 integral E Re Tr(tau3 deltaK) dE.
        A supplied K1 cannot silently be discarded in the first three laws.
        The energy increment's first term is a diagnostic, not a new validated
        total energy or a source from which missing heat is inferred.
        """
        if not isinstance(reference_contract,str) or not reference_contract.strip():
            raise ValueError('Explicit same-gap frozen thermal reference contract required')
        zero=_matrix(K_zero,self.shape,'K_zero')
        first=_matrix(K_first,self.shape,'K_first')
        reference=_matrix(K_thermal_zero,self.shape,'K_thermal_zero')
        phi=_real(phi_V,'phi_V')
        thermal=np.asarray(thermal_gradient_J_per_gapbar,complex)
        if phi.shape!=self.volume.shape or thermal.shape!=self.volume.shape or np.any(~np.isfinite(thermal)):
            raise ValueError('One finite voltage and complex thermal gradient per node required')
        flux0=np.asarray(edge_flux_zero,complex);flux1=np.asarray(edge_flux_first,complex)
        if (flux0.ndim!=4 or flux0.shape[0]!=len(self.energy) or flux0.shape[-2:]!=(2,2)
                or flux1.shape!=flux0.shape or np.any(~np.isfinite(flux0)) or np.any(~np.isfinite(flux1))):
            raise ValueError('Matching finite full Nambu edge fluxes per energy required')
        q0=np.einsum('e,en->n',self.weights,np.trace(zero,axis1=-2,axis2=-1).real/4)
        q1=np.einsum('e,en->n',self.weights,np.trace(first,axis1=-2,axis2=-1).real/4)
        # e is positive; potential follows the repository's B=Delta-iE tau3+i(e phi/Eu)I.
        from scipy.constants import elementary_charge
        neutral=phi+(q0+q1)/elementary_charge
        increment=zero-reference
        pair0=1j*self.N0*self.Eu*self.volume*np.einsum('e,en->n',self.weights,increment[:,:,0,1])
        pair1=1j*self.N0*self.Eu*self.volume*np.einsum('e,en->n',self.weights,first[:,:,0,1])
        def current(flux):
            traceT=(flux[:,:,0,0]-flux[:,:,1,1]).real/4
            return 2*self.I0/self.Eu*np.einsum('e,ek->k',self.weights,traceT)
        i0,i1=current(flux0),current(flux1)
        def energy(matrix):
            traceT=(matrix[:,:,0,0]-matrix[:,:,1,1]).real
            return -self.N0/2*float(np.einsum('e,n,en->',self.weights*self.energy,self.volume,traceT))
        return NambuMomentResult(neutral,q0,q1,thermal+pair0+pair1,pair0,pair1,
            i0+i1,i0,i1,energy(increment),energy(first),reference_contract)


    def complete_current(self,*,thermal_edge_current_A,edge_flux_zero,
                         thermal_edge_flux_zero,edge_flux_first,reference_contract):
        """Preserve an independently computed thermal current and its tails.

        I=Ithermal+integral(J0-Jthermal) + integral J1. The subtraction cancels
        a common ultraviolet equilibrium tail. It does not certify agreement
        of a direct real-axis current with the Matsubara oracle: both direct
        finite-quadrature currents and their discrepancy are returned. In DC,
        reproducing the baseline by this subtraction is an identity, never a
        new acceptance test. The reference must have the same instantaneous
        gap/phase/links and bath temperature, not merely the original DC gap.
        """
        if not isinstance(reference_contract,str) or not reference_contract.strip():
            raise ValueError('Explicit same-gap thermal current reference contract required')
        base=_real(thermal_edge_current_A,'thermal_edge_current_A')
        flux0=np.asarray(edge_flux_zero,complex);reference=np.asarray(thermal_edge_flux_zero,complex)
        flux1=np.asarray(edge_flux_first,complex)
        if (base.ndim!=1 or flux0.shape!=(len(self.energy),len(base),2,2)
                or reference.shape!=flux0.shape or flux1.shape!=flux0.shape
                or np.any(~np.isfinite(flux0)) or np.any(~np.isfinite(reference))
                or np.any(~np.isfinite(flux1))):
            raise ValueError('Matching finite full fluxes and one thermal current per edge required')
        def integrate(flux):
            projection=(flux[:,:,0,0]-flux[:,:,1,1]).real/4
            return 2*self.I0/self.Eu*np.einsum('e,ek->k',self.weights,projection)
        direct0=integrate(flux0);directref=integrate(reference)
        delta=integrate(flux0-reference);first=integrate(flux1)
        return StableCurrentCompletion(base+delta+first,base.copy(),delta,first,
            direct0+first,directref,directref-base,reference_contract)
