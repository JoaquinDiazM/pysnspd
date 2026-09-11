"""Opt-in uniform electronic catalogue for the model documented in v0.4.

This module is deliberately not imported by PRE, the kinetic solver, or gTDGL.
The amplitude is an independent coordinate. The vacuum potential and its two
derivatives share ONE interpolant; the complex retarded factors are evaluated
directly, not interpolated independently or silently repaired.

Dimensionless coordinates are a=|Delta|/Delta0, g=Gamma/Delta0 and E/Delta0.
The symbol ``a`` is confined to implementation details. Potentials are in
N0*Delta0**2 and N0 is the one-spin normal electronic DOS per joule and volume.
"""
from __future__ import annotations

from dataclasses import dataclass, replace
import hashlib
import json
from pathlib import Path
from typing import Any

import numpy as np
from scipy.interpolate import RectBivariateSpline, CubicSpline, CubicHermiteSpline
from scipy.optimize import brentq

HBAR_J_S = 1.054571817e-34
E_CHARGE_C = 1.602176634e-19
K_B_J_K = 1.380649e-23
BCS_GAP_RATIO = float(np.pi * np.exp(-np.euler_gamma))
SCHEMA = "pysnspd.experimental.uniform_vacuum_catalog.v1"


def _nonnegative(value: float, name: str) -> float:
    value = float(value)
    if not np.isfinite(value) or value < 0:
        raise ValueError(f"{name} must be finite and nonnegative")
    return value


def vacuum_state(delta: float, gamma: float) -> tuple[float, float, float]:
    """Return Uvac, dU/d|Delta| and dU/dGamma in normalized units.

    The zero-temperature Matsubara change of variable
    omega=u-Gamma*u/sqrt(u*u+Delta*Delta) gives
    dGamma Uvac=integral_0^infinity sin(Theta)**2 d(omega).
    Integrating in Gamma from the BCS vacuum gives the piecewise expression
    below, including the gapless Gamma>|Delta| branch. No gap equation is
    imposed. At exactly zero amplitude the normal reference is U=0.
    """
    delta, gamma = _nonnegative(delta, "delta"), _nonnegative(gamma, "gamma")
    if delta == 0:
        return 0.0, 0.0, 0.0
    if gamma <= delta:
        value = delta**2*(np.log(delta)-0.5) + np.pi*delta*gamma/2 - gamma**2/3
        force = 2*delta*np.log(delta) + np.pi*gamma/2
        conjugate = np.pi*delta/2 - 2*gamma/3
    else:
        ratio = delta/gamma
        t = np.sqrt((1-ratio)*(1+ratio))
        # Using log(Gamma), rather than cancelling log(delta)+log(Gamma/delta),
        # also controls the small-amplitude limit on the gapless branch.
        asin_over_ratio = np.arcsin(ratio)/ratio if ratio else 1.0
        k = (1+t+t*t)/(3*(1+t))
        h_without_log_ratio = asin_over_ratio + np.log1p(t) - k
        value = delta**2*(np.log(gamma)-0.5+h_without_log_ratio)
        # (1-t) = ratio**2/(1+t), avoiding subtractive cancellation.
        g_prime_over_ratio = (asin_over_ratio - ratio**2*(t+2)/(3*(1+t)**2))
        conjugate = delta*ratio*g_prime_over_ratio
        force = delta*(2*(np.log(gamma)+h_without_log_ratio)-g_prime_over_ratio)
    return float(value), float(force), float(conjugate)


def retarded_spectrum(energy: np.ndarray | float, *, delta: float,
                      gamma: float, eta: float) -> tuple[np.ndarray, np.ndarray]:
    """Causal c=N1+iR1, s=N2+iR2, sharing the robust scalar/batch solver.

    eta is a positive numerical regulator, not physical broadening. All roots
    are checked in the original unsquared spectral equation; invalid states
    raise. The earlier independent algorithm is preserved in release history.
    """
    return retarded_spectrum_batch(energy,delta=delta,gamma=gamma,eta=eta)


def _axis(values: np.ndarray, name: str, *, positive: bool = False) -> np.ndarray:
    axis = np.array(values, dtype=float, copy=True)
    if (axis.ndim != 1 or len(axis) < 4 or np.any(~np.isfinite(axis)) or
            np.any(np.diff(axis) <= 0) or np.any(axis <= 0 if positive else axis < 0)):
        raise ValueError(f"{name} needs at least four strictly increasing valid nodes")
    axis.setflags(write=False)
    return axis


@dataclass(frozen=True)
class UniformVacuumCatalog:
    """Vacuum spline with explicit positive-amplitude support and SI scales.

    The exact normal point is available separately; 0<|Delta|<delta_axis[0]
    is rejected. This pilot is not a licence to extrapolate into a vortex core.
    Spectral quadrature and E(x,p) interpolation remain separate admission
    checks before coupling this object to a nonequilibrium kinetic solver.
    """
    delta_axis: np.ndarray
    gamma_axis: np.ndarray
    potential: np.ndarray
    delta0_J: float
    N0_per_J_m3: float
    D_m2_s: float
    metadata: dict[str, Any]

    def __post_init__(self) -> None:
        deltas = _axis(self.delta_axis, "delta_axis", positive=True)
        gammas = _axis(self.gamma_axis, "gamma_axis")
        potential = np.array(self.potential, dtype=float, copy=True)
        if potential.shape != (len(deltas), len(gammas)) or np.any(~np.isfinite(potential)):
            raise ValueError("potential must be finite and match the two axes")
        for name in ("delta0_J", "N0_per_J_m3", "D_m2_s"):
            if not np.isfinite(getattr(self, name)) or getattr(self, name) <= 0:
                raise ValueError(f"{name} must be finite and positive")
        potential.setflags(write=False)
        object.__setattr__(self, "delta_axis", deltas)
        object.__setattr__(self, "gamma_axis", gammas)
        object.__setattr__(self, "potential", potential)
        object.__setattr__(self, "_spline", RectBivariateSpline(deltas, gammas, potential, s=0))

    def evaluate(self, delta: float, gamma: float) -> tuple[float, float, float]:
        """Return potential and its actual spline derivatives, without clipping."""
        delta, gamma = _nonnegative(delta, "delta"), _nonnegative(gamma, "gamma")
        if not self.gamma_axis[0] <= gamma <= self.gamma_axis[-1]:
            raise ValueError("Gamma query is outside catalogue support")
        if delta == 0:
            return 0.0, 0.0, 0.0
        if not self.delta_axis[0] <= delta <= self.delta_axis[-1]:
            raise ValueError("amplitude query is outside catalogue support")
        if self.metadata.get("vacuum_evaluation") == "closed_form":
            return vacuum_state(delta, gamma)
        return tuple(float(self._spline.ev(delta, gamma, dx=dx, dy=dy))
                     for dx, dy in ((0, 0), (1, 0), (0, 1)))

    def evaluate_si(self, delta_J: float, q_m_inv: float) -> dict[str, float]:
        """Map the uniform energy, force, and conjugate current to SI units."""
        if not np.isfinite(q_m_inv):
            raise ValueError("q_m_inv must be finite")
        gamma_J = HBAR_J_S*self.D_m2_s*q_m_inv**2/2
        potential, force, conjugate = self.evaluate(delta_J/self.delta0_J, gamma_J/self.delta0_J)
        pi = self.N0_per_J_m3*self.delta0_J*HBAR_J_S*self.D_m2_s*q_m_inv*conjugate
        return {"Uvac_J_m3": self.N0_per_J_m3*self.delta0_J**2*potential,
                "amplitude_force_per_m3": self.N0_per_J_m3*self.delta0_J*force,
                "Pi_J_m2": pi, "js_A_m2": 2*E_CHARGE_C/HBAR_J_S*pi,
                "Gamma_J": gamma_J}

    def save(self, path: str | Path) -> Path:
        """Persist arrays and JSON text only; loading never enables pickle."""
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        metadata = dict(self.metadata, schema=SCHEMA)
        with path.open("wb") as stream:
            np.savez_compressed(stream, delta_axis=self.delta_axis, gamma_axis=self.gamma_axis,
                                potential=self.potential, delta0_J=self.delta0_J,
                                N0_per_J_m3=self.N0_per_J_m3, D_m2_s=self.D_m2_s,
                                metadata_json=np.array(json.dumps(metadata, sort_keys=True, allow_nan=False)))
        return path

    @classmethod
    def load(cls, path: str | Path) -> "UniformVacuumCatalog":
        with np.load(path, allow_pickle=False) as values:
            metadata = json.loads(str(values["metadata_json"]))
            if metadata.get("schema") != SCHEMA:
                raise ValueError("unsupported experimental catalogue schema")
            return cls(values["delta_axis"], values["gamma_axis"], values["potential"],
                       float(values["delta0_J"]), float(values["N0_per_J_m3"]),
                       float(values["D_m2_s"]), metadata)


def build_vacuum_catalog(delta_axis: np.ndarray, gamma_axis: np.ndarray, *,
                         Tc_K: float, N0_per_J_m3: float, D_m2_s: float,
                         reference_hashes: dict[str, str] | None = None,
                         analytic: bool = False) -> UniformVacuumCatalog:
    """Build a small static table; callers determine and verify resolution."""
    deltas, gammas = _axis(delta_axis, "delta_axis", positive=True), _axis(gamma_axis, "gamma_axis")
    if not np.isfinite(Tc_K) or Tc_K <= 0:
        raise ValueError("Tc_K must be finite and positive")
    potential = np.array([[vacuum_state(d, g)[0] for g in gammas] for d in deltas])
    source_hash = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    metadata = {
        "schema": SCHEMA, "experimental": True, "production_connected": False,
        "model_documentation": "0.4", "reference_hashes": reference_hashes or {},
        "builder_sha256": source_hash, "Tc_K": float(Tc_K), "gap_ratio": BCS_GAP_RATIO,
        "units": {"axes": "Delta0", "potential": "N0 Delta0^2", "N0": "one spin / (J m^3)"},
        "amplitude_self_consistency": "not imposed",
        "vacuum_evaluation": "closed_form" if analytic else "bicubic",
        "interpolant": ("closed-form vacuum energy and its derivatives" if analytic else
                        "one bicubic potential; derivatives obtained from that spline"),
        "normal_endpoint": "exact separate state; no interpolation from zero to first positive node",
        "scope": "uniform vacuum only; no spatial PDE or nonthermal E(x,p) interpolation",
        "material_status": "reference electronic scales; not admission of absolute NbN phonon data",
        "spectral_policy": "direct complex causal solve; no independent c/s interpolation",
    }
    return UniformVacuumCatalog(deltas, gammas, potential, BCS_GAP_RATIO*K_B_J_K*Tc_K,
                               N0_per_J_m3, D_m2_s, metadata)


def spectral_count(energy: np.ndarray | float, *, delta: float,
                   gamma: float, eta: float) -> np.ndarray:
    """Return x(E)=integral_0^E Re(c) dE without numerical DOS quadrature.

    With z=u-i*Gamma*u/w and w=sqrt(u^2-Delta^2), a primitive of c dz is
    w-i*Gamma*Delta^2/(2*w^2). Its real part at E=0 is zero on the causal
    branch. This is a state-count coordinate, not the excitation energy.
    """
    c, s = retarded_spectrum(energy, delta=delta, gamma=gamma, eta=eta)
    if delta == 0:
        return np.asarray(energy, float)
    w = 1j*delta/s
    return np.real(w-1j*gamma*delta**2/(2*w*w))


def energy_at_count(count: np.ndarray | float, *, delta: float,
                    gamma: float, eta: float) -> np.ndarray:
    """Invert the causal count coordinate, preserving the finite eta support."""
    delta, gamma = _nonnegative(delta, "delta"), _nonnegative(gamma, "gamma")
    eta = _nonnegative(eta, "eta")
    counts = np.asarray(count, float)
    if eta == 0 or np.any(~np.isfinite(counts)) or np.any(counts < 0):
        raise ValueError("count must be finite/nonnegative and eta positive")
    if delta == 0:
        return counts.copy()
    if gamma == 0:
        return counts*np.sqrt(1+delta**2/(counts**2+eta**2))
    result = np.empty(counts.shape)
    for index in np.ndindex(counts.shape):
        x = counts[index]
        if x == 0:
            result[index] = 0
            continue
        upper = float(np.hypot(x, delta)+2*gamma+1)
        def residual(energy: float) -> float:
            return float(spectral_count(energy, delta=delta, gamma=gamma, eta=eta))-x
        result[index] = brentq(residual, 0, upper, xtol=2e-12, rtol=2e-13)
    return result


@dataclass(frozen=True)
class OccupationEnergyCatalog:
    """Finite-support fixed-p energy and its own interpolated derivatives.

    A single set of state-count quadrature nodes is shared by all amplitudes
    and depairing energies. p is held fixed at those nodes when differentiating.
    No Fermi-Dirac closure is imposed. The caller must check count cutoff,
    quadrature, eta and field-grid convergence before using this pilot.
    """
    vacuum: UniformVacuumCatalog
    count_nodes: np.ndarray
    count_weights: np.ndarray
    excitation_energies: np.ndarray
    eta: float
    gamma_energy_derivatives: np.ndarray
    gamma_log_offsets: np.ndarray | None = None
    gamma_coordinate_scale: float | None = None
    delta_energy_derivatives: np.ndarray | None = None
    gamma_energy_offsets: np.ndarray | None = None
    gamma_ratio_axis: np.ndarray | None = None

    def __post_init__(self) -> None:
        count = np.array(self.count_nodes, float, copy=True)
        weights = np.array(self.count_weights, float, copy=True)
        energies = np.array(self.excitation_energies, float, copy=True)
        gamma_derivatives = np.array(self.gamma_energy_derivatives, float, copy=True)
        axis=self.vacuum.gamma_axis
        if self.gamma_ratio_axis is not None:
            axis=_axis(self.gamma_ratio_axis,'gamma_ratio_axis')
            axis.setflags(write=False)
            object.__setattr__(self,'gamma_ratio_axis',axis)
            if self.gamma_energy_offsets is None or axis[0]!=0 or axis[-1]<self.vacuum.gamma_axis[-1]/self.vacuum.delta_axis[0]:
                raise ValueError('ratio representation requires direct energy and full physical Gamma support')
        object.__setattr__(self,'_parameter_axis',axis)
        shape = (len(self.vacuum.delta_axis),len(axis),len(count))
        if (count.ndim != 1 or len(count)<2 or weights.shape!=count.shape or
                np.any(~np.isfinite(count)) or np.any(np.diff(count)<=0) or np.any(count<=0) or
                np.any(~np.isfinite(weights)) or np.any(weights<=0) or
                energies.shape!=shape or np.any(~np.isfinite(energies)) or np.any(energies<=0) or
                gamma_derivatives.shape!=shape or np.any(~np.isfinite(gamma_derivatives)) or
                np.any(np.diff(energies,axis=-1)<=0) or not np.isfinite(self.eta) or self.eta<=0):
            raise ValueError("invalid fixed-count energy catalogue")
        for values in (count,weights,energies,gamma_derivatives):
            values.setflags(write=False)
        object.__setattr__(self,"count_nodes",count)
        object.__setattr__(self,"count_weights",weights)
        object.__setattr__(self,"excitation_energies",energies)
        object.__setattr__(self,"gamma_energy_derivatives",gamma_derivatives)
        if self.gamma_energy_offsets is not None:
            offsets=np.array(self.gamma_energy_offsets,float,copy=True)
            derivative=np.array(self.delta_energy_derivatives,float,copy=True)
            if (axis[0]!=0 or offsets.shape!=shape or derivative.shape!=shape or
                    np.any(~np.isfinite(offsets)) or np.any(~np.isfinite(derivative)) or np.any(offsets[:,0,:]!=0)):
                raise ValueError("direct energy offsets require Gamma=0 anchor and both nodal derivatives")
            offsets.setflags(write=False); derivative.setflags(write=False)
            object.__setattr__(self,'gamma_energy_offsets',offsets)
            object.__setattr__(self,'delta_energy_derivatives',derivative)
            a=self.vacuum.delta_axis[:,None,None]
            x=count[None,None,:]
            base_a=x*a/(x*x+self.eta*self.eta)/np.sqrt(1+a*a/(x*x+self.eta*self.eta))
            parameter_derivative=derivative
            slope_data=gamma_derivatives
            if self.gamma_ratio_axis is not None:
                parameter_derivative=derivative+axis[None,:,None]*gamma_derivatives
                slope_data=gamma_derivatives*a
            value_coefficients=CubicHermiteSpline(self.vacuum.delta_axis,offsets,parameter_derivative-base_a,axis=0).c
            if self.gamma_coordinate_scale is not None:
                if not np.isfinite(self.gamma_coordinate_scale) or self.gamma_coordinate_scale<=0:
                    raise ValueError("Gamma coordinate scale must be finite and positive")
                slope_data=slope_data*(axis[None,:,None]+self.gamma_coordinate_scale)
            slope_coefficients=CubicSpline(self.vacuum.delta_axis,slope_data,axis=0).c
            object.__setattr__(self,'_log_coefficients',value_coefficients)
            object.__setattr__(self,'_gamma_log_coefficients',slope_coefficients)
            object.__setattr__(self,'_base_log_coefficients',None)
            return
        increments=np.diff(energies,axis=-1,prepend=np.zeros((*energies.shape[:-1],1)))
        gamma_increments=np.diff(gamma_derivatives,axis=-1,prepend=np.zeros((*energies.shape[:-1],1)))
        if self.gamma_log_offsets is not None:
            offsets=np.array(self.gamma_log_offsets,float,copy=True)
            if (self.vacuum.gamma_axis[0]!=0 or offsets.shape!=shape or
                    np.any(~np.isfinite(offsets)) or np.any(offsets[:,0,:]!=0)):
                raise ValueError("compensated log offsets require Gamma=0 reference and finite matching data")
            offsets.setflags(write=False)
            object.__setattr__(self,"gamma_log_offsets",offsets)
            base=np.log(increments[:,0,:])
            increments=np.exp(base[:,None,:]+offsets)
            value_data=offsets
            base_coefficients=CubicSpline(self.vacuum.delta_axis,base,axis=0).c
        else:
            value_data=np.log(increments)
            base_coefficients=None
        # Exact nodal slopes resolve the narrow Gamma~eta^(3/2)/sqrt(Delta)
        # boundary layer. Both values and slopes refer to the SAME log increments.
        # Store polynomial coefficients so a query reads just its local cell.
        value_coefficients=CubicSpline(self.vacuum.delta_axis,value_data,axis=0).c
        slope_data=gamma_increments/increments
        if self.gamma_coordinate_scale is not None:
            if not np.isfinite(self.gamma_coordinate_scale) or self.gamma_coordinate_scale<=0:
                raise ValueError("Gamma coordinate scale must be finite and positive")
            slope_data=slope_data*(self.vacuum.gamma_axis[None,:,None]+self.gamma_coordinate_scale)
        slope_coefficients=CubicSpline(self.vacuum.delta_axis,slope_data,axis=0).c
        object.__setattr__(self,"_log_coefficients",value_coefficients)
        object.__setattr__(self,"_gamma_log_coefficients",slope_coefficients)
        object.__setattr__(self,"_base_log_coefficients",base_coefficients)

    def energy_kernel(self, delta: float, gamma: float) -> tuple[np.ndarray,np.ndarray,np.ndarray]:
        """Return E and physical field derivatives of the persisted interpolant."""
        self.vacuum.evaluate(delta,gamma)  # common support guard
        if delta==0:
            return self.count_nodes.copy(),np.zeros_like(self.count_nodes),np.zeros_like(self.count_nodes)
        i=min(max(int(np.searchsorted(self.vacuum.delta_axis,delta,side='right')-1),0),len(self.vacuum.delta_axis)-2)
        parameter=gamma/delta if self.gamma_ratio_axis is not None else gamma
        axis=self._parameter_axis
        j=min(max(int(np.searchsorted(axis,parameter,side='right')-1),0),len(axis)-2)
        offset=delta-self.vacuum.delta_axis[i]
        def delta_polynomial(coefficients):
            c=coefficients[:,i,j:j+2,:]
            value=((c[0]*offset+c[1])*offset+c[2])*offset+c[3]
            derivative=(3*c[0]*offset+2*c[1])*offset+c[2]
            return value,derivative
        values,values_d=delta_polynomial(self._log_coefficients)
        slopes,slopes_d=delta_polynomial(self._gamma_log_coefficients)
        if self.gamma_coordinate_scale is None:
            h=axis[j+1]-axis[j]
            t=(parameter-axis[j])/h
            chain=1.0
        else:
            scale=self.gamma_coordinate_scale
            # log1p of local differences retains precision in the first cell.
            h=np.log1p((axis[j+1]-axis[j])/(axis[j]+scale))
            t=np.log1p((parameter-axis[j])/(axis[j]+scale))/h
            chain=1/(parameter+scale)
        def hermite(v,m):
            difference=v[1]-v[0]
            return v[0]+t*(h*m[0]+t*(3*difference-h*(2*m[0]+m[1])+
                                    t*(-2*difference+h*(m[0]+m[1]))))
        log_g=(6*t*(1-t)*(values[1]-values[0])/h+(1-4*t+3*t*t)*slopes[0]+
               (-2*t+3*t*t)*slopes[1])*chain
        value_log=hermite(values,slopes)
        derivative_log=hermite(values_d,slopes_d)
        if self._base_log_coefficients is not None:
            c=self._base_log_coefficients[:,i,:]
            value_log+=((c[0]*offset+c[1])*offset+c[2])*offset+c[3]
            derivative_log+=(3*c[0]*offset+2*c[1])*offset+c[2]
        if self.gamma_energy_offsets is not None:
            if self.gamma_ratio_axis is not None:
                derivative_log-=parameter*log_g/delta
                log_g=log_g/delta
            x=self.count_nodes
            factor=np.sqrt(1+delta*delta/(x*x+self.eta*self.eta))
            arrays=(x*factor+value_log,x*delta/(x*x+self.eta*self.eta)/factor+derivative_log,log_g)
            if np.any(~np.isfinite(arrays)) or np.any(arrays[0]<=0) or np.any(np.diff(arrays[0])<=0):
                raise FloatingPointError(f"direct interpolated energy lost positive monotone support: Delta={delta}, Gamma={gamma}, minE={np.min(arrays[0]):.6g}, min increment={np.min(np.diff(arrays[0])):.6g}")
            return arrays
        logs=(value_log,derivative_log,log_g)
        increments=np.exp(logs[0])
        arrays=(np.cumsum(increments),np.cumsum(increments*logs[1]),np.cumsum(increments*logs[2]))
        if np.any(~np.isfinite(arrays)) or np.any(arrays[0]<=0) or np.any(np.diff(arrays[0])<=0):
            raise FloatingPointError("interpolated excitation energies lost positive monotone support")
        return arrays

    def evaluate(self, delta: float, gamma: float, occupation: np.ndarray) -> tuple[float,float,float]:
        p=np.asarray(occupation,float)
        if p.shape!=self.count_nodes.shape or np.any(~np.isfinite(p)) or np.any((p<0)|(p>1)):
            raise ValueError("occupation must match count nodes and lie in [0,1]")
        vacuum=np.asarray(self.vacuum.evaluate(delta,gamma))
        kernels=self.energy_kernel(delta,gamma)
        return tuple(float(base+4*np.dot(kernel*p,self.count_weights))
                     for base,kernel in zip(vacuum,kernels))

    def evaluate_si(self,delta_J: float,q_m_inv: float,occupation: np.ndarray) -> dict[str,float]:
        """SI energy, amplitude force and uniform current at fixed p(x)."""
        if not np.isfinite(q_m_inv):
            raise ValueError("q_m_inv must be finite")
        scale=self.vacuum.delta0_J
        N0,D=self.vacuum.N0_per_J_m3,self.vacuum.D_m2_s
        gamma_J=HBAR_J_S*D*q_m_inv**2/2
        energy,force,conjugate=self.evaluate(delta_J/scale,gamma_J/scale,occupation)
        pi=N0*scale*HBAR_J_S*D*q_m_inv*conjugate
        return {"energy_J_m3":N0*scale**2*energy,"amplitude_force_per_m3":N0*scale*force,
                "Pi_J_m2":pi,"js_A_m2":2*E_CHARGE_C/HBAR_J_S*pi,"Gamma_J":gamma_J}

    def save(self,path: str | Path) -> Path:
        path=Path(path)
        path.parent.mkdir(parents=True,exist_ok=True)
        schema=("pysnspd.experimental.occupation_catalog.v5" if self.gamma_ratio_axis is not None else "pysnspd.experimental.occupation_catalog.v4" if self.gamma_energy_offsets is not None else "pysnspd.experimental.occupation_catalog.v3" if self.gamma_log_offsets is not None or self.gamma_coordinate_scale is not None else
                "pysnspd.experimental.occupation_catalog.v2")
        metadata=dict(self.vacuum.metadata, schema=schema,
                      scope="uniform fixed-p finite-count quadrature; no spatial PDE",
                      eta_over_Delta0=self.eta,count_max=float(np.sum(self.count_weights)),
                      count_quadrature="explicit nodes and weights (builder uses Gauss-Legendre)",count_order=len(self.count_nodes),
                      array_units={"excitation_energies":"Delta0","count_nodes":"Delta0",
                                   "count_weights":"Delta0","potential":"N0 Delta0^2",
                                   "gamma_energy_derivatives":"d(E/Delta0)/d(Gamma/Delta0) at fixed x"},
                      occupation_policy="p(x) held fixed; no thermal reduction",
                      gamma_slope_reference="-N2 R2/N1 from the same causal finite-eta spectrum at each node",
                      gamma_coordinate_scale=self.gamma_coordinate_scale,
                      interpolation_policy=("one Uvac energy ("+self.vacuum.metadata.get('vacuum_evaluation','bicubic')+"); log positive E increments: cubic in Delta, Hermite in "+('Gamma' if self.gamma_coordinate_scale is None else 'log1p(Gamma/scale)')+" with exact nodal spectral slopes; all forces differentiated"))
        if self.gamma_log_offsets is not None:
            metadata['gamma_value_representation']='separate Gamma=0 base plus compensated log-increment offsets'
        extra={} if self.gamma_log_offsets is None else {'gamma_log_offsets':self.gamma_log_offsets}
        if self.gamma_energy_offsets is not None:
            extra.update(gamma_energy_offsets=self.gamma_energy_offsets,delta_energy_derivatives=self.delta_energy_derivatives)
            metadata['array_units'].update(gamma_energy_offsets='Delta0; E(Gamma)-E(0) at fixed count and amplitude',
                                           delta_energy_derivatives='d(E/Delta0)/d(absDelta/Delta0) at fixed physical Gamma and count')
            metadata['energy_representation']='direct E_BCS analytic plus compensated Gamma energy correction'
            if self.gamma_ratio_axis is not None:
                extra['gamma_ratio_axis']=self.gamma_ratio_axis
                metadata['interpolation_field_coordinate']='r=Gamma/absDelta; physical support from vacuum axes'
                metadata['ratio_derivative_chain']='E_Delta|r=E_Delta|Gamma+r E_Gamma; E_r=Delta E_Gamma; physical forces restore chain'
                metadata['array_units']['gamma_ratio_axis']='dimensionless Gamma/absDelta'
            metadata['interpolation_policy']='common energy; Hermite Delta with exact nodal amplitude slope; Hermite Gamma coordinate with causal nodal slope; positivity and monotonicity checked without repair'
        with path.open("wb") as stream:
            np.savez_compressed(stream,delta_axis=self.vacuum.delta_axis,gamma_axis=self.vacuum.gamma_axis,
                                potential=self.vacuum.potential,delta0_J=self.vacuum.delta0_J,
                                N0_per_J_m3=self.vacuum.N0_per_J_m3,D_m2_s=self.vacuum.D_m2_s,
                                count_nodes=self.count_nodes,count_weights=self.count_weights,
                                excitation_energies=self.excitation_energies,
                                gamma_energy_derivatives=self.gamma_energy_derivatives,
                                metadata_json=np.array(json.dumps(metadata,sort_keys=True,allow_nan=False)),**extra)
        return path

    @classmethod
    def load(cls,path: str | Path) -> "OccupationEnergyCatalog":
        with np.load(path,allow_pickle=False) as values:
            metadata=json.loads(str(values["metadata_json"]))
            if metadata.get("schema") not in ("pysnspd.experimental.occupation_catalog.v2",
                                               "pysnspd.experimental.occupation_catalog.v3", "pysnspd.experimental.occupation_catalog.v4", "pysnspd.experimental.occupation_catalog.v5"):
                raise ValueError("unsupported occupation catalogue schema")
            vacuum=UniformVacuumCatalog(values["delta_axis"],values["gamma_axis"],values["potential"],
                                        float(values["delta0_J"]),float(values["N0_per_J_m3"]),
                                        float(values["D_m2_s"]),metadata)
            return cls(vacuum,values["count_nodes"],values["count_weights"],
                       values["excitation_energies"],float(metadata["eta_over_Delta0"]),
                       values["gamma_energy_derivatives"],
                       values["gamma_log_offsets"] if 'gamma_log_offsets' in values else None,
                       metadata.get('gamma_coordinate_scale'),
                       values['delta_energy_derivatives'] if 'delta_energy_derivatives' in values else None,
                       values['gamma_energy_offsets'] if 'gamma_energy_offsets' in values else None,
                       values['gamma_ratio_axis'] if 'gamma_ratio_axis' in values else None)


def build_occupation_catalog(vacuum: UniformVacuumCatalog, *, count_order: int=48,
                             count_max: float=6, eta: float=1e-3,
                             count_nodes: np.ndarray | None=None,
                             count_weights: np.ndarray | None=None,
                             inverse_method: str="scalar", compensated_gamma: bool=False,
                             gamma_coordinate_scale: float | None=None,
                             gamma_offset_threshold: float=1e-10,
                             gamma_offset_order: int=8, direct_energy: bool=False,
                             gamma_ratio_axis: np.ndarray | None=None) -> OccupationEnergyCatalog:
    """Build a small Gauss-Legendre fixed-count pilot, with no production hooks."""
    if not isinstance(count_order,int) or count_order<2 or not np.isfinite(count_max) or count_max<=0:
        raise ValueError("positive count support and count_order>=2 required")
    if (count_nodes is None)!=(count_weights is None):
        raise ValueError("provide both count_nodes and count_weights")
    if inverse_method not in ("scalar","batch"):
        raise ValueError("inverse_method must be scalar or batch")
    if count_nodes is None:
        nodes,weights=np.polynomial.legendre.leggauss(count_order)
        count=(nodes+1)*count_max/2
        weights=weights*count_max/2
    else:
        count,weights=np.asarray(count_nodes,float),np.asarray(count_weights,float)
        if (count.ndim!=1 or len(count)<2 or weights.shape!=count.shape or
                np.any(~np.isfinite(count)) or np.any(count<=0) or np.any(np.diff(count)<=0) or
                np.any(~np.isfinite(weights)) or np.any(weights<=0)):
            raise ValueError("custom state-count nodes and weights must be positive, finite and ordered")
    if gamma_ratio_axis is not None and (not direct_energy or not compensated_gamma):
        raise ValueError('ratio coordinate requires direct compensated energy')
    energies,gamma_derivatives,delta_derivatives,offsets,offset_diagnostics=[],[],[],[],[]
    for delta in vacuum.delta_axis:
        energy_row,derivative_row,delta_row=[],[],[]
        field_gammas=vacuum.gamma_axis if gamma_ratio_axis is None else delta*np.asarray(gamma_ratio_axis)
        for gamma in field_gammas:
            inverse=energy_at_count_batch if inverse_method=="batch" else energy_at_count
            spectral=retarded_spectrum_batch if inverse_method=="batch" else retarded_spectrum
            energy=inverse(count,delta=delta,gamma=gamma,eta=eta)
            c,s=spectral(energy,delta=delta,gamma=gamma,eta=eta)
            energy_row.append(energy)
            derivative_row.append(-s.real*s.imag/c.real)
            delta_row.append(s.imag/c.real)
        energies.append(energy_row)
        gamma_derivatives.append(derivative_row)
        delta_derivatives.append(delta_row)
        if compensated_gamma:
            correction,diagnostics=integrated_gamma_offsets(count,field_gammas,delta=delta,eta=eta,direct_energy=direct_energy,
                energies=np.asarray(energy_row),threshold=gamma_offset_threshold,order=gamma_offset_order)
            offsets.append(correction)
            offset_diagnostics.append(diagnostics)
    if compensated_gamma:
        vacuum.metadata['gamma_offset_diagnostics']=offset_diagnostics
    return OccupationEnergyCatalog(vacuum,count,weights,np.asarray(energies),eta,np.asarray(gamma_derivatives),
                                   np.asarray(offsets) if compensated_gamma and not direct_energy else None,gamma_coordinate_scale,
                                   np.asarray(delta_derivatives) if direct_energy else None,
                                   np.asarray(offsets) if direct_energy else None,gamma_ratio_axis)


def segmented_count_quadrature(breaks: np.ndarray, orders: int | np.ndarray) -> tuple[np.ndarray,np.ndarray]:
    """Positive Gauss rules on fixed global count intervals, not field-dependent nodes."""
    breaks=np.asarray(breaks,float)
    if breaks.ndim!=1 or len(breaks)<2 or breaks[0]!=0 or np.any(~np.isfinite(breaks)) or np.any(np.diff(breaks)<=0):
        raise ValueError("count breaks must start at zero and strictly increase")
    raw_orders=np.asarray(orders)
    if np.any(~np.isfinite(raw_orders)) or np.any(raw_orders!=np.floor(raw_orders)):
        raise ValueError("count quadrature orders must be integers")
    orders=np.broadcast_to(np.asarray(raw_orders,int),(len(breaks)-1,))
    if np.any(orders<2):
        raise ValueError("each count interval requires at least two nodes")
    nodes,weights=[],[]
    for lo,hi,order in zip(breaks[:-1],breaks[1:],orders):
        x,w=np.polynomial.legendre.leggauss(int(order))
        nodes.append(lo+(x+1)*(hi-lo)/2)
        weights.append(w*(hi-lo)/2)
    return np.concatenate(nodes),np.concatenate(weights)


def retarded_spectrum_batch(energy: np.ndarray | float, *, delta: float,
                            gamma: float, eta: float) -> tuple[np.ndarray,np.ndarray]:
    """Causal batched quartic solve with residual-checked original-equation refinement.

    This opt-in implementation reduces the cost of constructing a refined
    catalogue. Scalar and array entry points share this validated algorithm.
    No spectral roots or DOS values are clipped into admissibility.
    """
    delta,gamma=_nonnegative(delta,"delta"),_nonnegative(gamma,"gamma")
    eta=_nonnegative(eta,"eta")
    energies=np.asarray(energy,float)
    if eta==0 or np.any(~np.isfinite(energies)) or np.any(energies<0):
        raise ValueError("energy must be finite/nonnegative and eta positive")
    shape=energies.shape
    z=energies.reshape(-1).astype(complex)+1j*eta
    if delta==0:
        return np.ones(shape,complex),np.zeros(shape,complex)
    if gamma==0:
        w=np.sqrt(z*z-delta*delta)
        w=np.where(w.imag<0,-w,w)
        return (z/w).reshape(shape),(1j*delta/w).reshape(shape)
    zero=energies.reshape(-1)==0
    if np.any(zero):
        alpha=brentq(lambda value:value-eta-gamma*value/np.hypot(value,delta),
                     eta,eta+gamma,xtol=1e-30,rtol=1e-15)
        radius=np.hypot(alpha,delta)
        c=np.full(len(z),alpha/radius,complex)
        s=np.full(len(z),delta/radius,complex)
        if np.any(~zero):
            c[~zero],s[~zero]=retarded_spectrum_batch(energies.reshape(-1)[~zero],delta=delta,gamma=gamma,eta=eta)
        return c.reshape(shape),s.reshape(shape)
    scale=np.maximum(np.maximum(abs(z),delta),gamma)
    matrix=np.zeros((len(z),4,4),complex)
    matrix[:,0,:]=-np.stack([-2*z/scale,(z*z-delta*delta+gamma*gamma)/scale**2,
                             2*z*delta*delta/scale**3,-z*z*delta*delta/scale**4],axis=1)
    matrix[:,1,0]=1;matrix[:,2,1]=1;matrix[:,3,2]=1
    roots=np.linalg.eigvals(matrix)*scale[:,None]
    w_bcs=np.sqrt(z*z-delta*delta)
    w_bcs=np.where(w_bcs.imag<0,-w_bcs,w_bcs)
    # At tiny Gamma the quartic has almost repeated roots. A causal BCS
    # perturbation seed retains information lost by a polynomial eigensolve.
    # It is screened and refined in the same original equation as every root.
    bcs_seed=z+1j*gamma*z/w_bcs
    roots=np.column_stack([roots,bcs_seed])
    w=np.sqrt(roots*roots-delta*delta)
    w=np.where(w.imag<0,-w,w)
    residual=roots-z[:,None]-1j*gamma*roots/w
    valid=(roots.imag>0)&(roots.real>=-1e-10*scale[:,None])
    scores=np.where(valid,abs(residual),np.inf)
    if np.any(np.all(~valid,axis=1)):
        raise FloatingPointError("no causal root in batch")
    # Refine every causal seed before selecting. Near the spectral edge a seed
    # with a smaller initial residual can lie in the wrong Newton basin.
    u=roots.copy()
    for _ in range(24):
        w=np.sqrt(u*u-delta*delta)
        w=np.where(w.imag<0,-w,w)
        residual=u-z[:,None]-1j*gamma*u/w
        active=valid&(abs(residual)>2e-14*np.maximum(1,scale[:,None]))
        if not np.any(active):
            break
        step=np.zeros_like(u)
        denominator=np.ones_like(u)
        denominator[active]=1+1j*gamma*delta*delta/w[active]**3
        active &= np.isfinite(denominator)&(abs(denominator)>0)
        step[active]=residual[active]/denominator[active]
        accepted=np.zeros(u.shape,bool)
        for damping in (1,.5,.25,.125,.0625,.03125):
            candidate=u-damping*step
            wc=np.sqrt(candidate*candidate-delta*delta)
            wc=np.where(wc.imag<0,-wc,wc)
            rc=candidate-z[:,None]-1j*gamma*candidate/wc
            use=active&~accepted&(candidate.imag>0)&(abs(rc)<abs(residual))
            u=np.where(use,candidate,u)
            accepted|=use
        if not np.any(accepted):
            break
    w=np.sqrt(u*u-delta*delta)
    w=np.where(w.imag<0,-w,w)
    u_candidates_c,u_candidates_s=u/w,1j*delta/w
    u_candidates_valid=valid.copy()
    # Near-gap cubic seeds retain the small w scale when the quartic in u
    # becomes nearly degenerate. They are only seeds: every accepted value
    # is refined and checked against the original, unapproximated equation.
    cubic_p=2*delta*(delta-z)
    cubic_q=-2j*gamma*delta*delta
    w_scale=np.maximum(np.sqrt(abs(cubic_p)),abs(cubic_q)**(1/3))
    cubic=np.zeros((len(z),3,3),complex)
    cubic[:,0,1]=-cubic_p/w_scale**2
    cubic[:,0,2]=-cubic_q/w_scale**3
    cubic[:,1,0]=1; cubic[:,2,1]=1
    w_seeds=np.linalg.eigvals(cubic)*w_scale[:,None]
    w=np.column_stack([w,w_seeds])
    valid=np.column_stack([valid,(w_seeds.real>=0)&(w_seeds.imag>0)])
    # Refine in w near the edge: computing sqrt(u*u-Delta**2) loses the
    # small gap-edge quantity once u and Delta almost coincide. The stable
    # residual below factors u-z, and keeps w as the independent unknown.
    def w_residual(w_value):
        uv=np.sqrt(w_value*w_value+delta*delta)
        return (w_value*w_value+(delta-z[:,None])*(delta+z[:,None]))/(uv+z[:,None])-1j*gamma*uv/w_value
    for _ in range(32):
        u=np.sqrt(w*w+delta*delta)
        residual=w_residual(w)
        active=valid&(abs(u)>0)&(abs(w)>0)&np.isfinite(residual)&(abs(residual)>2e-16*np.maximum(1,scale[:,None]))
        if not np.any(active): break
        step=np.zeros_like(w)
        derivative=np.ones_like(w)
        derivative[active]=w[active]/u[active]+1j*gamma*delta*delta/(u[active]*w[active]**2)
        active &= np.isfinite(derivative)&(abs(derivative)>0)
        step[active]=residual[active]/derivative[active]
        accepted=np.zeros(w.shape,bool)
        for damping in (1,.5,.25,.125,.0625,.03125):
            candidate=w-damping*step
            rc=w_residual(candidate)
            use=active&~accepted&(candidate.real>=0)&(candidate.imag>0)&(abs(rc)<abs(residual))
            w=np.where(use,candidate,w)
            accepted|=use
        if not np.any(accepted): break
    # The u coordinate is accurate around E=0, and w around E=Delta. Select
    # using the physical equation rather than discarding either representation.
    wc=np.sqrt(w*w+delta*delta)/w
    ws=1j*delta/w
    all_c=np.column_stack([u_candidates_c,wc])
    all_s=np.column_stack([u_candidates_s,ws])
    all_valid=np.column_stack([u_candidates_valid,valid])
    lh=delta*all_c; rh=(gamma*all_c-1j*z[:,None])*all_s
    scores=abs(lh-rh)/np.maximum(1,np.maximum(abs(lh),abs(rh)))
    all_valid &= np.isfinite(all_c+all_s)&(all_c.real>=0)&(all_s.real>=0)
    scores=np.where(all_valid,scores,np.inf)
    selected=np.argmin(scores,axis=1)
    c=all_c[np.arange(len(z)),selected]
    s=all_s[np.arange(len(z)),selected]
    lhs,rhs=delta*c,(gamma*c-1j*z)*s
    equation=abs(lhs-rhs)/np.maximum(1,np.maximum(abs(lhs),abs(rhs)))
    normalization=abs(c*c+s*s-1)/np.maximum(1,abs(c)**2+abs(s)**2)
    if (np.any(~np.isfinite(c+s)) or np.any(c.real<0) or
            np.any(equation>1e-8) or np.any(normalization>1e-10)):
        worst=int(np.argmax(equation))
        raise FloatingPointError(f"inadmissible batch spectrum: delta={delta:g}, gamma={gamma:g}, "
                                 f"eta={eta:g}, E={energies.reshape(-1)[worst]:.17g}, residual={np.max(equation):.3g}")
    return c.reshape(shape),s.reshape(shape)


def energy_at_count_batch(count: np.ndarray | float, *, delta: float,
                          gamma: float, eta: float) -> np.ndarray:
    """Bracketed vector Newton inverse of the causal state count.

    At fixed count, depairing lowers E monotonically from its BCS value toward
    the normal value x. Thus [x,E_BCS(x)] is a tight physical bracket. Each
    Newton step uses dx/dE=N1 and is replaced by bisection if it leaves it.
    """
    delta,gamma=_nonnegative(delta,"delta"),_nonnegative(gamma,"gamma")
    eta=_nonnegative(eta,"eta")
    counts=np.asarray(count,float)
    if eta==0 or np.any(~np.isfinite(counts)) or np.any(counts<0):
        raise ValueError("count must be finite/nonnegative and eta positive")
    if delta==0 or gamma==0:
        return energy_at_count(counts,delta=delta,gamma=gamma,eta=eta)
    shape=counts.shape
    x=counts.reshape(-1)
    lo=x.copy()
    hi=energy_at_count(x,delta=delta,gamma=0,eta=eta)
    energy=hi.copy()
    active=x>0
    for _ in range(80):
        if not np.any(active):
            return energy.reshape(shape)
        indices=np.flatnonzero(active)
        e=energy[indices]
        c,s=retarded_spectrum_batch(e,delta=delta,gamma=gamma,eta=eta)
        w=1j*delta/s
        actual=np.real(w-1j*gamma*delta*delta/(2*w*w))
        residual=actual-x[indices]
        hi[indices]=np.where(residual>=0,e,hi[indices])
        lo[indices]=np.where(residual<0,e,lo[indices])
        step=residual/c.real
        tolerance=4e-14*np.maximum(1,abs(e))
        done=(abs(step)<=tolerance)|((hi[indices]-lo[indices])<=tolerance)
        active[indices[done]]=False
        candidate=e-step
        safe=(candidate>lo[indices])&(candidate<hi[indices])&np.isfinite(candidate)
        # A small but representable final correction must not be discarded.
        finish=done&safe
        energy[indices[finish]]=candidate[finish]
        candidate=np.where(safe,candidate,(lo[indices]+hi[indices])/2)
        energy[indices[~done]]=candidate[~done]
    raise FloatingPointError("state-count batch inversion did not converge within 80 iterations")


def integrated_gamma_offsets(count_nodes, gamma_axis, *, delta, eta,
                             energies=None, gamma_derivatives=None,
                             threshold=1e-10, order=8, direct_energy=False):
    """Return offsets[nGamma,nX], diagnostics with a Gamma=0 anchor.

    Up to threshold integrate dGamma log(E_i-E_(i-1)). Above threshold use
    direct log ratios. The endpoint discrepancy at the threshold is reported,
    never removed by a fitted constant or hidden renormalization. Optional
    nodal arrays avoid redundant nodal inversions; quadrature nodes still need
    their own spectrum. E_(minus one)=0, dGamma E_(minus one)=0.
    """
    import time
    started = time.perf_counter()
    x, gammas = np.asarray(count_nodes, float), np.asarray(gamma_axis, float)
    if (x.ndim != 1 or np.any(~np.isfinite(x)) or np.any(x <= 0)
            or np.any(np.diff(x) <= 0) or gammas.ndim != 1 or len(gammas) < 2
            or gammas[0] != 0 or np.any(~np.isfinite(gammas))
            or np.any(np.diff(gammas) <= 0) or threshold <= 0 or order < 2):
        raise ValueError("ordered positive counts, Gamma=0 anchor and positive quadrature required")
    if energies is None:
        energies = np.array([energy_at_count_batch(x, delta=delta, gamma=g, eta=eta)
                             for g in gammas])
    energies = np.asarray(energies, float)
    if energies.shape != (len(gammas), len(x)) or np.any(~np.isfinite(energies)):
        raise ValueError("invalid nodal energy shape/values")
    inc = np.diff(energies, axis=-1, prepend=0)
    if np.any(inc <= 0):
        raise FloatingPointError("nonpositive nodal energy increment")
    if gamma_derivatives is not None:
        deriv = np.asarray(gamma_derivatives, float)
        if deriv.shape != energies.shape or np.any(~np.isfinite(deriv)):
            raise ValueError("invalid supplied nodal Gamma derivative")
    # log1p does not recover lost endpoint differences; it merely avoids adding
    # an additional log subtraction when the direct diagnostic is resolvable.
    direct = energies-energies[0] if direct_energy else np.log1p((inc - inc[0]) / inc[0])
    offsets = direct.copy()
    offsets[0] = 0
    nodes, weights = np.polynomial.legendre.leggauss(order)
    accum = np.zeros_like(x)
    compensation = np.zeros_like(x)
    calls = 0
    min_inc = float(np.min(inc))
    last = 0.0
    threshold = min(float(threshold), float(gammas[-1]))

    def integrate(lo, hi):
        nonlocal calls, min_inc
        values = []
        for g in lo + (nodes + 1) * (hi - lo) / 2:
            energy = energy_at_count_batch(x, delta=delta, gamma=float(g), eta=eta)
            c, s = retarded_spectrum_batch(energy, delta=delta, gamma=float(g), eta=eta)
            increments = np.diff(energy, prepend=0)
            if np.any(increments <= 0):
                raise FloatingPointError("nonpositive quadrature-node energy increment")
            min_inc = min(min_inc, float(np.min(increments)))
            slope = -s.real*s.imag/c.real if direct_energy else np.diff(-s.real * s.imag / c.real, prepend=0) / increments
            if np.any(~np.isfinite(slope)):
                raise FloatingPointError("nonfinite Gamma log slope")
            values.append(slope)
            calls += 1
        return (hi - lo) / 2 * np.einsum("j,jk->k", weights, np.array(values))

    for i, endpoint in enumerate(gammas[1:], 1):
        endpoint = min(float(endpoint), threshold)
        if endpoint <= last:
            break
        addition = integrate(last, endpoint) - compensation
        updated = accum + addition
        compensation = (updated - accum) - addition
        accum = updated
        last = endpoint
        if gammas[i] <= threshold:
            offsets[i] = accum
        if endpoint == threshold:
            break

    threshold_energy = energy_at_count_batch(x, delta=delta, gamma=threshold, eta=eta)
    threshold_inc = np.diff(threshold_energy, prepend=0)
    if np.any(threshold_inc <= 0):
        raise FloatingPointError("nonpositive transition energy increment")
    direct_at_threshold = threshold_energy-energies[0] if direct_energy else np.log1p((threshold_inc - inc[0]) / inc[0])
    transition_discrepancy = accum - direct_at_threshold
    reconstructed = energies[0]+offsets if direct_energy else np.cumsum(inc[0] * np.exp(offsets), axis=-1)
    diagnostics = {
        "delta": float(delta), "eta": float(eta), "order": int(order),
        "offset_units": "Delta0" if direct_energy else "log dimensionless energy increment",
        "threshold": threshold, "shape": list(offsets.shape),
        "quadrature_spectrum_calls": calls,
        "minimum_energy_increment": min_inc,
        "anchor_max_abs": float(np.max(abs(offsets[0]))),
        "max_abs_transition_log_discrepancy": float(np.max(abs(transition_discrepancy))),
        "transition_worst_count": float(x[np.argmax(abs(transition_discrepancy))]),
        "max_abs_reconstructed_nodal_energy_discrepancy": float(np.max(abs(reconstructed-energies))),
        "max_abs_direct_log_discrepancy_below_threshold": float(np.max(abs(offsets-direct)[gammas <= threshold])),
        "runtime_seconds": time.perf_counter()-started,
        "interpretation": "Nodal potential/quadrature diagnostic; not final integrated-current error.",
    }
    return offsets, diagnostics


def refine_occupation_catalog(seed: OccupationEnergyCatalog, ratio_axis: np.ndarray,
                              *, reference_hash: str | None=None) -> OccupationEnergyCatalog:
    """Insert ratio nodes while retaining every existing nodal value exactly.

    New nodes solve the same causal spectrum. Below the compensation threshold
    their energy offset is integrated from the nearest existing left anchor;
    above it a direct energy difference is well resolved. No value or derivative
    is adjusted to make an acceptance check pass.
    """
    if seed.gamma_ratio_axis is None or seed.gamma_energy_offsets is None:
        raise ValueError('refinement requires a direct-energy ratio catalogue')
    axis=_axis(ratio_axis,'ratio_axis')
    old=seed.gamma_ratio_axis
    if axis[0]!=old[0] or axis[-1]!=old[-1] or not np.all(np.isin(old,axis)):
        raise ValueError('refinement must retain every old node and the same support')
    shape=(len(seed.vacuum.delta_axis),len(axis),len(seed.count_nodes))
    energy,da,dg,offset=[np.empty(shape) for _ in range(4)]
    old_positions=np.searchsorted(axis,old)
    for new,original in zip((energy,da,dg,offset),(seed.excitation_energies,
                seed.delta_energy_derivatives,seed.gamma_energy_derivatives,seed.gamma_energy_offsets)):
        new[:,old_positions,:]=original
    new_positions=np.flatnonzero(~np.isin(axis,old))
    nodes,weights=np.polynomial.legendre.leggauss(8)
    threshold=1e-10
    for i,delta in enumerate(seed.vacuum.delta_axis):
        for j in new_positions:
            gamma=delta*axis[j]
            e=energy_at_count_batch(seed.count_nodes,delta=delta,gamma=gamma,eta=seed.eta)
            c,s=retarded_spectrum_batch(e,delta=delta,gamma=gamma,eta=seed.eta)
            energy[i,j]=e; da[i,j]=s.imag/c.real; dg[i,j]=-s.real*s.imag/c.real
            if gamma>threshold:
                offset[i,j]=e-seed.excitation_energies[i,0]
            else:
                left=int(np.searchsorted(old,axis[j])-1)
                start=delta*old[left]
                integral=np.zeros_like(seed.count_nodes)
                for point,weight in zip(start+(nodes+1)*(gamma-start)/2,weights):
                    eq=energy_at_count_batch(seed.count_nodes,delta=delta,gamma=point,eta=seed.eta)
                    cq,sq=retarded_spectrum_batch(eq,delta=delta,gamma=point,eta=seed.eta)
                    integral+=weight*(-sq.real*sq.imag/cq.real)
                offset[i,j]=seed.gamma_energy_offsets[i,left]+integral*(gamma-start)/2
    metadata=dict(seed.vacuum.metadata)
    metadata.update(builder_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                    refinement_seed_sha256=reference_hash,
                    refinement_new_ratio_nodes=axis[new_positions].tolist(),
                    refinement_policy='old nodes retained exactly; new nodes from the same causal equations',
                    refinement_offset_policy='below physical Gamma=1e-10 integrate E_Gamma from existing left anchor with Gauss8; otherwise use direct E-E_BCS')
    vacuum=replace(seed.vacuum,metadata=metadata)
    return OccupationEnergyCatalog(vacuum,seed.count_nodes,seed.count_weights,energy,seed.eta,dg,
                                   None,seed.gamma_coordinate_scale,da,offset,axis)

