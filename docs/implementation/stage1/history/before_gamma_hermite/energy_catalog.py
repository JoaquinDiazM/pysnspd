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

from dataclasses import dataclass
import hashlib
import json
from pathlib import Path
from typing import Any

import numpy as np
from scipy.interpolate import RectBivariateSpline
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
    """Return c=N1+iR1, s=N2+iR2 on the causal positive-energy branch.

    eta is strictly positive and numerical, not a physical broadening claim.
    Squared-equation roots are screened against the unsquared equation.
    Invalid states raise instead of applying clipping or nan_to_num.
    """
    delta, gamma = _nonnegative(delta, "delta"), _nonnegative(gamma, "gamma")
    eta = _nonnegative(eta, "eta")
    if eta == 0:
        raise ValueError("eta must be strictly positive; converge its limit separately")
    energies = np.asarray(energy, dtype=float)
    if np.any(~np.isfinite(energies)) or np.any(energies < 0):
        raise ValueError("energy must be finite and nonnegative")
    if delta == 0:
        return np.ones(energies.shape, complex), np.zeros(energies.shape, complex)
    c_out, s_out = np.empty(energies.shape, complex), np.empty(energies.shape, complex)
    for index in np.ndindex(energies.shape):
        z = complex(energies[index], eta)
        if gamma == 0:
            u = z
        else:
            roots = np.roots([1, -2*z, z*z-delta*delta+gamma*gamma,
                              2*z*delta*delta, -z*z*delta*delta])
            candidates = []
            for u_candidate in roots:
                w = np.sqrt(u_candidate*u_candidate-delta*delta)
                if w.imag < 0:
                    w = -w
                if u_candidate.imag > 0 and u_candidate.real >= -1e-10*max(1, delta, gamma):
                    score = abs(u_candidate-z-1j*gamma*u_candidate/w)
                    candidates.append((score, u_candidate))
            if not candidates:
                raise FloatingPointError("no causal retarded root")
            u = min(candidates, key=lambda item: item[0])[1]
        w = np.sqrt(u*u-delta*delta)
        if w.imag < 0:
            w = -w
        c, s = u/w, 1j*delta/w
        lhs, rhs = delta*c, (gamma*c-1j*z)*s
        residual = abs(lhs-rhs)/max(1, abs(lhs), abs(rhs))
        normalization = abs(c*c+s*s-1)/max(1, abs(c)**2+abs(s)**2)
        if (not np.isfinite(c+s) or c.real < -1e-10 or
                residual > 1e-8 or normalization > 1e-10):
            raise FloatingPointError(f"inadmissible retarded solution: residual={residual:g}")
        c_out[index], s_out[index] = c, s
    return c_out, s_out


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
                         reference_hashes: dict[str, str] | None = None) -> UniformVacuumCatalog:
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
        "interpolant": "one bicubic potential; derivatives obtained from that spline",
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

    def __post_init__(self) -> None:
        count = np.array(self.count_nodes, float, copy=True)
        weights = np.array(self.count_weights, float, copy=True)
        energies = np.array(self.excitation_energies, float, copy=True)
        shape = (len(self.vacuum.delta_axis),len(self.vacuum.gamma_axis),len(count))
        if (count.ndim != 1 or len(count)<2 or weights.shape!=count.shape or
                np.any(~np.isfinite(count)) or np.any(np.diff(count)<=0) or np.any(count<=0) or
                np.any(~np.isfinite(weights)) or np.any(weights<=0) or
                energies.shape!=shape or np.any(~np.isfinite(energies)) or np.any(energies<=0) or
                np.any(np.diff(energies,axis=-1)<=0) or not np.isfinite(self.eta) or self.eta<=0):
            raise ValueError("invalid fixed-count energy catalogue")
        for values in (count,weights,energies):
            values.setflags(write=False)
        object.__setattr__(self,"count_nodes",count)
        object.__setattr__(self,"count_weights",weights)
        object.__setattr__(self,"excitation_energies",energies)
        increments=np.diff(energies,axis=-1,prepend=np.zeros((*energies.shape[:-1],1)))
        object.__setattr__(self,"_energy_splines",[
            RectBivariateSpline(self.vacuum.delta_axis,self.vacuum.gamma_axis,np.log(increments[:,:,k]),s=0)
            for k in range(len(count))])

    def energy_kernel(self, delta: float, gamma: float) -> tuple[np.ndarray,np.ndarray,np.ndarray]:
        """Return E and derivatives from positive, cumulatively summed increments."""
        self.vacuum.evaluate(delta,gamma)  # common support guard
        if delta==0:
            return self.count_nodes.copy(),np.zeros_like(self.count_nodes),np.zeros_like(self.count_nodes)
        logs=tuple(np.array([float(s.ev(delta,gamma,dx=dx,dy=dy)) for s in self._energy_splines])
                   for dx,dy in ((0,0),(1,0),(0,1)))
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
        metadata=dict(self.vacuum.metadata, schema="pysnspd.experimental.occupation_catalog.v1",
                      scope="uniform fixed-p finite-count quadrature; no spatial PDE",
                      eta_over_Delta0=self.eta,count_max=float(np.sum(self.count_weights)),
                      count_quadrature="explicit nodes and weights (builder uses Gauss-Legendre)",count_order=len(self.count_nodes),
                      array_units={"excitation_energies":"Delta0","count_nodes":"Delta0",
                                   "count_weights":"Delta0","potential":"N0 Delta0^2"},
                      occupation_policy="p(x) held fixed; no thermal reduction",
                      interpolation_policy="one Uvac spline; cubic splines of log positive E increments in x; reconstructed E and all forces differentiated")
        with path.open("wb") as stream:
            np.savez_compressed(stream,delta_axis=self.vacuum.delta_axis,gamma_axis=self.vacuum.gamma_axis,
                                potential=self.vacuum.potential,delta0_J=self.vacuum.delta0_J,
                                N0_per_J_m3=self.vacuum.N0_per_J_m3,D_m2_s=self.vacuum.D_m2_s,
                                count_nodes=self.count_nodes,count_weights=self.count_weights,
                                excitation_energies=self.excitation_energies,
                                metadata_json=np.array(json.dumps(metadata,sort_keys=True,allow_nan=False)))
        return path

    @classmethod
    def load(cls,path: str | Path) -> "OccupationEnergyCatalog":
        with np.load(path,allow_pickle=False) as values:
            metadata=json.loads(str(values["metadata_json"]))
            if metadata.get("schema")!="pysnspd.experimental.occupation_catalog.v1":
                raise ValueError("unsupported occupation catalogue schema")
            vacuum=UniformVacuumCatalog(values["delta_axis"],values["gamma_axis"],values["potential"],
                                        float(values["delta0_J"]),float(values["N0_per_J_m3"]),
                                        float(values["D_m2_s"]),metadata)
            return cls(vacuum,values["count_nodes"],values["count_weights"],
                       values["excitation_energies"],float(metadata["eta_over_Delta0"]))


def build_occupation_catalog(vacuum: UniformVacuumCatalog, *, count_order: int=48,
                             count_max: float=6, eta: float=1e-3) -> OccupationEnergyCatalog:
    """Build a small Gauss-Legendre fixed-count pilot, with no production hooks."""
    if not isinstance(count_order,int) or count_order<2 or not np.isfinite(count_max) or count_max<=0:
        raise ValueError("positive count support and count_order>=2 required")
    nodes,weights=np.polynomial.legendre.leggauss(count_order)
    count=(nodes+1)*count_max/2
    weights=weights*count_max/2
    energies=np.array([[energy_at_count(count,delta=d,gamma=g,eta=eta) for g in vacuum.gamma_axis]
                       for d in vacuum.delta_axis])
    return OccupationEnergyCatalog(vacuum,count,weights,energies,eta)
