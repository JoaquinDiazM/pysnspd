"""Explicit scales, inherited KWT mobility, and experimental phonon inputs.

The common units are E/Delta0, U/(N0 Delta0**2), and t/t_ref. KWT D.11--13
is an effective dynamical closure, not derived from the static catalogue.
Its tau_psi is deliberately unrelated to the separately supplied BGK tau_kin.
No production loader or solver imports this module.
"""
from __future__ import annotations
from dataclasses import dataclass
import hashlib
import json
from pathlib import Path
import numpy as np
from scipy.special import eval_legendre, roots_jacobi

from .energy_catalog import HBAR_J_S, K_B_J_K


def _number(value, name, *, zero=False):
    value = float(value)
    if not np.isfinite(value) or value < 0 or (value == 0 and not zero):
        raise ValueError(f"{name} must be finite and {'nonnegative' if zero else 'positive'}")
    return value


def _finite_array(value, name):
    value = np.asarray(value, float)
    if np.any(~np.isfinite(value)):
        raise ValueError(f"{name} must be finite")
    return value


@dataclass(frozen=True)
class CellScales:
    delta0_J: float
    N0_per_J_m3: float
    Tc_K: float
    Tb_K: float
    t_ref_ps: float

    def __post_init__(self):
        for name in ("delta0_J", "N0_per_J_m3", "Tc_K", "Tb_K", "t_ref_ps"):
            object.__setattr__(self, name, _number(getattr(self, name), name))

    @property
    def t_ref_s(self):
        return self.t_ref_ps * 1e-12

    @property
    def energy_density_scale_J_m3(self):
        return self.N0_per_J_m3 * self.delta0_J**2

    @property
    def power_density_scale_W_m3(self):
        return self.energy_density_scale_J_m3 / self.t_ref_s

    @property
    def bath_temperature_bar(self):
        return K_B_J_K * self.Tb_K / self.delta0_J

    @property
    def quantum_rate_scale(self):
        """Delta0*t_ref/hbar, without a collision-specific pi or degeneracy."""
        return self.delta0_J * self.t_ref_s / HBAR_J_S

    def temperature_to_bar(self, temperature_K):
        return K_B_J_K * _number(temperature_K, "temperature_K", zero=True) / self.delta0_J

    def temperature_bar(self, temperature_K):
        """Driver-facing spelling of temperature_to_bar."""
        return self.temperature_to_bar(temperature_K)

    def temperature_to_K(self, temperature_bar):
        return _number(temperature_bar, "temperature_bar", zero=True) * self.delta0_J / K_B_J_K

    def metadata(self):
        return {name: getattr(self, name) for name in self.__dataclass_fields__} | {
            "energy_density_scale_J_m3": self.energy_density_scale_J_m3,
            "temperature_convention": "temperature_bar = k_B T / Delta0",
            "force_convention": "X = d[U/(N0 Delta0^2)] / d[Delta/Delta0]",
            "time_convention": "dimensionless time = physical time / t_ref",
        }


@dataclass(frozen=True)
class KWTResponse:
    velocity: float | np.ndarray
    heat: float
    mobility: float | np.ndarray
    R: float
    taupsi_ps: float
    Tmob_K: float


@dataclass(frozen=True)
class KWTMobility:
    scales: CellScales

    def coefficients(self, amplitude, temperature_bar):
        amplitude = _number(amplitude, "amplitude", zero=True)
        temperature_K = self.scales.temperature_to_K(temperature_bar)
        tmob = max(temperature_K, self.scales.Tb_K)
        ratio = tmob / self.scales.Tc_K
        try:
            taupsi_ps = 1.0 / (ratio/.50 + ratio**3/2.47)
            tau0_ps = np.pi*HBAR_J_S/(8*K_B_J_K*self.scales.Tc_K)*1e12
            abar = np.sqrt((1+ratio)/2)
            c = 4*(self.scales.delta0_J*taupsi_ps*1e-12/HBAR_J_S)**2
            stretch = np.sqrt(1+c*amplitude**2)
            tau0bar = tau0_ps/self.scales.t_ref_ps
        except (OverflowError, ZeroDivisionError) as exc:
            raise ValueError("KWT coefficients are outside finite representable scales") from exc
        values = dict(Abar=float(abar), tau0_ps=float(tau0_ps), tau0bar=float(tau0bar),
                      c=float(c), R=float(stretch), taupsi_ps=float(taupsi_ps), Tmob_K=float(tmob))
        if not all(np.isfinite(v) and v > 0 for k,v in values.items() if k != "c") or not np.isfinite(c):
            raise ValueError("KWT coefficients are not finite and positive")
        return values

    def amplitude_response(self, amplitude, force, temperature_bar):
        """Uniform radial D.12: a_dot=-X/(2*Abar*tau0bar*R).

        The factor two follows from the complex variational force F=X/2.
        Returns power per dimensionless time; Q=-X*a_dot, deposited once.
        """
        force = float(_finite_array(force, "force"))
        coefficients = self.coefficients(amplitude, temperature_bar)
        mobility = 1/(2*coefficients["Abar"]*coefficients["tau0bar"]*coefficients["R"])
        velocity = -mobility*force
        heat = mobility*force**2
        if not np.isfinite(velocity) or not np.isfinite(heat):
            raise ValueError("KWT response is not finite")
        return KWTResponse(velocity,float(heat),float(mobility),coefficients["R"],
                           coefficients["taupsi_ps"],coefficients["Tmob_K"])

    def tensor_response(self, z, force, temperature_bar):
        """Cartesian D.12 for real energy gradient force=(dUbar/dz1,dUbar/dz2).

        Eigenvalues are positive radial 1/(2*Abar*tau0bar*R) and tangential
        R/(2*Abar*tau0bar). Resolving these directions avoids subtracting two
        nearly equal matrices in the inverse of I+c*z*z.T. At z=0 use I.
        """
        z,force = _finite_array(z,"z"),_finite_array(force,"force")
        if z.shape != (2,) or force.shape != (2,):
            raise ValueError("z and force must each contain two real components")
        amplitude = float(np.hypot(*z))
        coefficients = self.coefficients(amplitude,temperature_bar)
        denominator = 2*coefficients["Abar"]*coefficients["tau0bar"]
        radial = 1/(denominator*coefficients["R"])
        tangential = coefficients["R"]/denominator
        if amplitude == 0:
            matrix=radial*np.eye(2); velocity=-radial*force
            heat=radial*float(np.dot(force,force))
        else:
            direction=z/amplitude
            transverse=np.array([-direction[1],direction[0]])
            xr,xt=float(np.dot(direction,force)),float(np.dot(transverse,force))
            matrix=radial*np.outer(direction,direction)+tangential*np.outer(transverse,transverse)
            velocity=-radial*xr*direction-tangential*xt*transverse
            heat=radial*xr*xr+tangential*xt*xt
        if np.any(~np.isfinite(velocity)) or not np.isfinite(heat):
            raise ValueError("KWT tensor response is not finite")
        matrix.setflags(write=False); velocity.setflags(write=False)
        return KWTResponse(velocity,float(heat),matrix,coefficients["R"],
                           coefficients["taupsi_ps"],coefficients["Tmob_K"])


def bose_occupation(energy_bar, temperature_bar):
    """Bose occupation at strictly positive energy, including exact T=0."""
    energies=_finite_array(energy_bar,"energy_bar")
    temperature=_number(temperature_bar,"temperature_bar",zero=True)
    if np.any(energies<=0):
        raise ValueError("Bose occupation requires positive energies; represent the infrared cutoff explicitly")
    if temperature==0:
        return np.zeros_like(energies)
    x=energies/temperature
    with np.errstate(over="ignore"):
        result=np.exp(-x)/(-np.expm1(-x))
    if np.any(~np.isfinite(result)):
        raise ValueError("Bose occupation is outside finite representable scales")
    return result


def _lobatto(nodes, low, high):
    if not isinstance(nodes,(int,np.integer)) or isinstance(nodes,bool) or nodes<3:
        raise ValueError("at least three Gauss-Lobatto nodes are required")
    low,high=_number(low,"infrared_cutoff_bar"),_number(high,"cutoff_energy_bar")
    if low>=high:
        raise ValueError("infrared cutoff must lie strictly below the upper cutoff")
    middle,_=roots_jacobi(nodes-2,1,1)
    rule=np.r_[-1.,middle,1.]
    weights=2/(nodes*(nodes-1)*eval_legendre(nodes-1,rule)**2)
    energies=low+(rule+1)*(high-low)/2
    energies[0],energies[-1]=low,high
    return energies,weights*(high-low)/2


@dataclass(frozen=True)
class PhononQuadrature:
    energies: np.ndarray
    weights: np.ndarray
    dos: np.ndarray
    metadata: dict

    def __post_init__(self):
        e,w,g=(_finite_array(v,n).copy() for v,n in [(self.energies,"energies"),(self.weights,"weights"),(self.dos,"dos")])
        if e.ndim!=1 or len(e)<2 or w.shape!=e.shape or g.shape!=e.shape or np.any(e<=0) or np.any(np.diff(e)<=0) or np.any(w<=0) or np.any(g<0):
            raise ValueError("positive ordered phonon energies/weights and nonnegative matching DOS required")
        for name,array in [("energies",e),("weights",w),("dos",g)]:
            array.setflags(write=False); object.__setattr__(self,name,array)
        try:
            with np.errstate(over="raise",invalid="raise"):
                cap=w*g
        except FloatingPointError as exc:
            raise ValueError("phonon capacities are not finite after applying quadrature weights") from exc
        cap.setflags(write=False); object.__setattr__(self,"capacities",cap)

    def energy(self, occupation):
        n=_finite_array(occupation,"phonon occupation")
        if n.shape!=self.energies.shape or np.any(n<0):
            raise ValueError("phonon occupation must match its grid and be nonnegative")
        with np.errstate(over="ignore",invalid="ignore"):
            value=float(np.dot(self.capacities,self.energies*n))
        if not np.isfinite(value): raise ValueError("phonon energy is not finite")
        return value

    def thermal_energy(self, temperature_bar):
        return self.energy(bose_occupation(self.energies,temperature_bar))

    def thermal_capacity(self, temperature_bar):
        """dUbar/d(k_B T/Delta0), not SI J/(m^3 K)."""
        temperature=_number(temperature_bar,"temperature_bar",zero=True)
        if temperature==0: return 0.0
        x=self.energies/temperature
        with np.errstate(over="ignore",invalid="ignore"):
            kernel=x*x*np.exp(-x)/(-np.expm1(-x))**2
        # x=inf is the limiting zero integrand, not a population repair.
        kernel=np.where(np.isinf(x),0.,kernel)
        if np.any(~np.isfinite(kernel)): raise ValueError("thermal capacity is not finite")
        return float(np.dot(self.capacities,kernel))


@dataclass(frozen=True)
class DebyePhonons:
    scales: CellScales
    cutoff_energy_bar: float
    atom_density_m3: float
    lambda_eph: float
    synthetic_label: str
    infrared_cutoff_bar: float=0.0

    def __post_init__(self):
        for name in ("cutoff_energy_bar","atom_density_m3","lambda_eph"):
            object.__setattr__(self,name,_number(getattr(self,name),name))
        low=_number(self.infrared_cutoff_bar,"infrared_cutoff_bar",zero=True)
        if low>=self.cutoff_energy_bar: raise ValueError("infrared cutoff must be below the Debye cutoff")
        object.__setattr__(self,"infrared_cutoff_bar",low)
        if not isinstance(self.synthetic_label,str) or not self.synthetic_label.strip():
            raise ValueError("a nonempty synthetic scenario label is required")

    @property
    def coupling_support(self):
        return (self.infrared_cutoff_bar,self.cutoff_energy_bar)

    def _support(self, energy_bar):
        energy=_finite_array(energy_bar,"energy_bar")
        if np.any(energy<0): raise ValueError("phonon energies must be nonnegative")
        return energy,(energy>=self.infrared_cutoff_bar)&(energy<=self.cutoff_energy_bar)

    def alpha2F(self, energy_bar):
        energy,support=self._support(energy_bar)
        return np.where(support,self.lambda_eph*(energy/self.cutoff_energy_bar)**2,0.)

    def dos_bar(self, energy_bar):
        """gSI(Delta0*Omega_bar)/N0; includes the essential 1/Delta0."""
        energy,support=self._support(energy_bar)
        prefactor=9*self.atom_density_m3/(self.scales.N0_per_J_m3*self.scales.delta0_J)
        return np.where(support,prefactor*energy**2/self.cutoff_energy_bar**3,0.)

    def dos_SI(self, energy_J):
        return self.scales.N0_per_J_m3*self.dos_bar(np.asarray(energy_J)/self.scales.delta0_J)

    def metadata(self):
        ratio=self.infrared_cutoff_bar/self.cutoff_energy_bar
        return {"status":"SYNTHETIC_ANALYTIC_NOT_NBN","label":self.synthetic_label,
            "scales":self.scales.metadata(),"cutoff_energy_bar":self.cutoff_energy_bar,
            "infrared_cutoff_bar":self.infrared_cutoff_bar,"atom_density_m3":self.atom_density_m3,
            "lambda_eph_full_Debye":self.lambda_eph,
            "lambda_eph_retained":self.lambda_eph*(1-ratio**2),
            "modes_per_m3_retained":3*self.atom_density_m3*(1-ratio**3),
            "infrared_omitted_fractions":{"lambda":ratio**2,"modes":ratio**3,"first_DOS_energy_moment":ratio**4},
            "support_policy":"compact declared interval, endpoints included; zero outside is a model cutoff, not extrapolation",
            "normalization":"gSI=9*n_atom*Omega^2/OmegaD^3; gbar=gSI/N0; measure dOmega_bar",
        }

    def quadrature(self,nodes=65):
        """Endpoint-inclusive Lobatto on the explicitly truncated interval.

        For event grids choose a positive infrared cutoff before construction.
        With no cutoff the analytic functions remain available, but this rule
        refuses a Bose-singular node at zero. Nothing silently drops that node.
        """
        energies,weights=_lobatto(nodes,self.infrared_cutoff_bar,self.cutoff_energy_bar)
        return PhononQuadrature(energies,weights,self.dos_bar(energies),
                               self.metadata()|{"quadrature":"Gauss-Lobatto","nodes":int(nodes)})


class ConditionalPhononShape:
    """A hash-bound derived shape under two explicitly synthetic scale factors.

    source_axis_to_energy_bar maps the original abscissa to Omega/Delta0.
    dos_ordinate_to_dos_bar is the complete independent ordinate scale gSI/N0.
    Neither factor is inferred from a target mode count or called a certified
    NbN conversion. This class intentionally has no absolute-SI admission API.
    """
    def __init__(self,csv_path,manifest_path,*,source_axis_to_energy_bar,
                 dos_ordinate_to_dos_bar,synthetic_label,infrared_cutoff_bar=0.0):
        factor=_number(source_axis_to_energy_bar,"source_axis_to_energy_bar")
        density_factor=_number(dos_ordinate_to_dos_bar,"dos_ordinate_to_dos_bar")
        if not isinstance(synthetic_label,str) or not synthetic_label.strip():
            raise ValueError("a declared synthetic mapping label is required")
        raw=Path(csv_path).read_bytes(); manifest_raw=Path(manifest_path).read_bytes()
        manifest=json.loads(manifest_raw)
        if manifest.get("schema")!="pysnspd.derived-phonon-shape.v1" or manifest.get("derived_csv_sha256")!=hashlib.sha256(raw).hexdigest():
            raise ValueError("derived material schema or CSV hash mismatch")
        data=np.loadtxt(csv_path,delimiter=",",skiprows=1,ndmin=2)
        if data.shape[1]!=4 or len(data)<2 or np.any(~np.isfinite(data)) or np.any(np.diff(data[:,0])<=0) or np.any(data[:,:3]<0):
            raise ValueError("invalid derived shape table")
        if np.any(data[data[:,2]==0,1]!=0) or (data[0,0]==0 and data[0,1]!=0):
            raise ValueError("derived shape must obey common support and origin integrability")
        try:
            with np.errstate(over="raise",invalid="raise"):
                self.energies=data[:,0]*factor; self.dos=data[:,2]*density_factor
        except FloatingPointError as exc:
            raise ValueError("conditional shape mapping produces nonfinite coordinates or DOS") from exc
        self.alpha=data[:,1].copy()
        if np.any(np.diff(self.energies)<=0) or np.any((self.alpha>0)&(self.dos<=0)):
            raise ValueError("conditional shape mapping loses ordered coordinates or common support")
        self.infrared_cutoff_bar=_number(infrared_cutoff_bar,"infrared_cutoff_bar",zero=True)
        if self.infrared_cutoff_bar>=self.energies[-1]: raise ValueError("infrared cutoff outside source support")
        for v in (self.energies,self.alpha,self.dos): v.setflags(write=False)
        self.coupling_support=(max(self.infrared_cutoff_bar,float(self.energies[0])),float(self.energies[-1]))
        self._metadata={"status":"CONDITIONAL_SYNTHETIC_MAPPING_NOT_ABSOLUTE_NBN",
            "synthetic_label":synthetic_label,"source_axis_to_energy_bar":factor,
            "dos_ordinate_to_dos_bar":density_factor,"infrared_cutoff_bar":self.infrared_cutoff_bar,
            "csv_sha256":hashlib.sha256(raw).hexdigest(),"manifest_sha256":hashlib.sha256(manifest_raw).hexdigest(),
            "source_sha256":manifest.get("source_sha256"),"source_url":manifest.get("source_url"),
            "original_operations":manifest.get("operations"),"SI_admitted":False,
            "scope":"Finite source shape plus explicit infrared truncation only; no inferred DOS units, atom/cell basis, physical density or absolute rate calibration."}

    def _evaluate(self,energy_bar,values):
        q=_finite_array(energy_bar,"energy_bar")
        if np.any(q<0): raise ValueError("phonon energies must be nonnegative")
        result=np.interp(q,self.energies,values,left=0,right=0)
        return np.where(q>=self.infrared_cutoff_bar,result,0.)

    def alpha2F(self,energy_bar): return self._evaluate(energy_bar,self.alpha)
    def dos_bar(self,energy_bar): return self._evaluate(energy_bar,self.dos)
    def metadata(self): return dict(self._metadata)
