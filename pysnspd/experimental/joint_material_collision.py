"""Full Nambu Born electron-phonon building block for joint stage 5.1.

Primary: Vargunin/Silaev, PRB96,214507(2017), Eqs7,8,15--17,
https://arxiv.org/abs/1709.01444. The dissipative acoustic signed kernel is
written with supplied dimensionless Eliashberg alpha2F. Its pi/4 prefactor
is fixed by the normal golden-rule limit (2*pi/hbar)*alpha2F*dOmega.
No principal-value dispersive renormalization is supplied. This is a weak
Born material hypothesis, not an admitted quantitative NbN photon model.

Actual R,H at E+/-Omega must be supplied. No extrapolation, eta bath, Pauli
clipping or temperature replacement is performed. An independent general
phonon bubble is supplied below. Pairing its physical shifted-energy
quadrature with this self-energy must be validated before material evolution.
"""
from __future__ import annotations
from dataclasses import dataclass
import numpy as np
from scipy.constants import hbar
from scipy.special import zeta
from . import adiabatic_matrix_usadel as ad

_FIELDS=("value","energy","time","first")
SOURCE="https://arxiv.org/abs/1709.01444"


def _positive(value,name):
    raw=np.asarray(value)
    if (raw.ndim or np.iscomplexobj(raw) or raw.dtype.kind=="b"
            or not np.isfinite(raw) or raw<=0):
        raise ValueError(name+" must be finite positive real")
    return float(raw)


def _array(values,name,nonnegative=False):
    raw=np.asarray(values)
    if np.iscomplexobj(raw) or raw.dtype.kind=="b":
        raise ValueError(name+" must be real")
    result=np.array(raw,float,copy=True)
    if np.any(~np.isfinite(result)) or (nonnegative and np.any(result<0)):
        raise ValueError(name+" must be finite"+(" and nonnegative" if nonnegative else ""))
    return result


def _multiply(jet,value,time):
    value,time=np.asarray(value)[...,None,None],np.asarray(time)[...,None,None]
    return ad.Jet(jet.value*value,jet.energy*value,
                  jet.time*value+jet.value*time,jet.first*value)


def _integrate(jet,weights):
    return ad.Jet(*(np.tensordot(weights,getattr(jet,name),axes=(0,0)) for name in _FIELDS))


def self_energy_advanced(retarded_J):
    """Sigma_A=+tau3 Sigma_R^dagger tau3, unlike the g_A minus sign."""
    adj=retarded_J.dagger()
    return ad.Jet(*(ad.TAU3@getattr(adj,name)@ad.TAU3 for name in _FIELDS))


@dataclass(frozen=True)
class InelasticCalibrators:
    """Allmaras2020 nominal FREE reference scales; no range imposed."""
    tau_ee_Tc_s: float=5e-12
    tau_ep_Tc_s: float=24.7e-12

    def __post_init__(self):
        for name in ("tau_ee_Tc_s","tau_ep_Tc_s"):
            object.__setattr__(self,name,_positive(getattr(self,name),name))

    @property
    def tau0_s(self):
        return float(720*zeta(5)/np.pi**2*self.tau_ep_Tc_s)

    def tau_sc_s(self,electron_temperature_K,critical_temperature_K):
        ratio=_positive(electron_temperature_K,"electron_temperature_K")/_positive(
            critical_temperature_K,"critical_temperature_K")
        return 1/(ratio/self.tau_ee_Tc_s+ratio**3/self.tau_ep_Tc_s)


@dataclass(frozen=True)
class RAKSelfEnergy:
    """Self energies in J; dE is w.r.t E/Eu, dt per physical second."""
    retarded_J: ad.Jet
    advanced_J: ad.Jet
    keldysh_J: ad.Jet
    energy_unit_J: float
    approximation: str
    source: str=SOURCE
    phonon_feedback_complete: bool=False

    def __post_init__(self):
        object.__setattr__(self,"energy_unit_J",_positive(self.energy_unit_J,"energy_unit_J"))
        jets=(self.retarded_J,self.advanced_J,self.keldysh_J)
        if any(not isinstance(jet,ad.Jet) for jet in jets):
            raise ValueError("full retarded, advanced and Keldysh Jets required")
        if any(jet.value.shape!=jets[0].value.shape for jet in jets):
            raise ValueError("R/A/K shapes must agree")
        if not isinstance(self.approximation,str) or not self.approximation.strip():
            raise ValueError("explicit approximation required")

    @property
    def retarded_generator(self):
        """+i Sigma_R/Eu added to B_R=D-iE tau3."""
        return self.retarded_J.scale(1j/self.energy_unit_J)

    @property
    def advanced_generator(self):
        return self.advanced_J.scale(1j/self.energy_unit_J)

    @property
    def keldysh_generator(self):
        return self.keldysh_J.scale(1j/self.energy_unit_J)

    def symmetry_errors_J(self):
        expected=self_energy_advanced(self.retarded_J)
        adj=self.keldysh_J.dagger()
        return {
            "advanced_symmetry_max_J":max(float(np.max(abs(getattr(self.advanced_J,f)
                -getattr(expected,f)),initial=0)) for f in _FIELDS),
            "keldysh_symmetry_max_J":max(float(np.max(abs(getattr(self.keldysh_J,f)
                +ad.TAU3@getattr(adj,f)@ad.TAU3),initial=0)) for f in _FIELDS)}


class ElectronPhononBornQuadrature:
    """Full local dissipative R/A/K kernel on supplied shifted full Jets.

    Inputs have shape (phonon,*target,2,2), n shape (phonon,) or
    (phonon,*target). The caller identifies actual physical E+/-Omega
    queries, with no silently repeated endpoint or artificial FD tails.
    """
    def __init__(self,phonon_energy_J,quadrature_weights_J,alpha2F,*,
                 energy_unit_J,shift_contract):
        omega=_array(phonon_energy_J,"phonon_energy_J")
        weights=_array(quadrature_weights_J,"quadrature_weights_J")
        alpha=_array(alpha2F,"alpha2F",True)
        if (omega.ndim!=1 or not len(omega) or np.any(omega<=0)
                or np.any(np.diff(omega)<=0) or weights.shape!=omega.shape
                or alpha.shape!=omega.shape or np.any(weights<=0)):
            raise ValueError("ordered positive energies and matching positive quadrature required")
        if not isinstance(shift_contract,str) or not shift_contract.strip():
            raise ValueError("explicit physical energy-shift contract required")
        self.energy_unit_J=_positive(energy_unit_J,"energy_unit_J")
        self.phonon_energy_J,self.quadrature_weights_J,self.alpha2F=omega,weights,alpha
        self.shift_contract=shift_contract
        for value in (omega,weights,alpha):value.setflags(write=False)

    @staticmethod
    def _field(values,shape,name,nonnegative):
        field=_array(values,name,nonnegative)
        if field.shape==(shape[0],):
            field=field.reshape((shape[0],)+(1,)*(len(shape)-1))
        try:return np.broadcast_to(field,shape)
        except ValueError as error:
            raise ValueError(name+" must broadcast to (phonon,*target)") from error

    def evaluate(self,retarded_plus,distribution_plus,retarded_minus,
                 distribution_minus,phonon_occupation,*,
                 phonon_time_per_s=0.,kappa_s=None):
        jets=(retarded_plus,distribution_plus,retarded_minus,distribution_minus)
        if any(not isinstance(jet,ad.Jet) for jet in jets):
            raise ValueError("shifted full Nambu Jets required")
        shape=retarded_plus.value.shape
        if shape[0]!=len(self.phonon_energy_J) or any(j.value.shape!=shape for j in jets):
            raise ValueError("matching (phonon,*target,2,2) Jets required")
        expected=hbar/(2*self.energy_unit_J)
        kappa=expected if kappa_s is None else _positive(kappa_s,"kappa_s")
        if not np.isclose(kappa,expected,rtol=4e-15,atol=0):
            raise ValueError("kappa_s must be hbar/(2*energy_unit_J)")
        n=self._field(phonon_occupation,shape[:-2],"phonon_occupation",True)
        ndot=self._field(phonon_time_per_s,shape[:-2],"phonon_time_per_s",False)
        b,bdot=2*n+1,2*ndot
        Rp,Rm=retarded_plus,retarded_minus
        Ap,Am=ad.advanced(Rp),ad.advanced(Rm)
        Kp,Km=ad.keldysh(Rp,distribution_plus,kappa),ad.keldysh(Rm,distribution_minus,kappa)
        r=_multiply(Rp+Rm,2*b,2*bdot)-(Kp-Km)
        a=_multiply(Ap+Am,2*b,2*bdot)+(Kp-Km)
        k=(_multiply(Kp+Km,b,bdot)-((Rp-Ap)-(Rm-Am))).scale(2)
        weights=self.quadrature_weights_J*self.alpha2F
        sigma=tuple(_integrate(j,weights).scale(-1j*np.pi/4) for j in (r,a,k))
        return RAKSelfEnergy(*sigma,self.energy_unit_J,
            "Local dissipative Born e-ph; supplied full shifted spectra and alpha2F; "
            "normal golden-rule normalization; dispersive part omitted; independent phonon bubble pairing required")


@dataclass(frozen=True)
class RAKCollisionSource:
    source: ad.Jet
    collision_rate_per_s: ad.Jet
    self_energy: RAKSelfEnergy
    source_contract: str

    @property
    def projected_source(self):
        return ad.projections(self.source.value+self.source.first)


def collision_source(retarded,distribution,self_energy,area_weights,*,kappa_s=None):
    """C=m*i/(2Eu)*(Sigma_R K+Sigma_K A-R Sigma_K-K Sigma_A).

    All products are truncated star products. Equivalently C=-kappa*m*I_K
    with I_K=i/hbar*(R Sigma_K+K Sigma_A-Sigma_R K-Sigma_K A).
    Normal homogeneous Pi(I_K)=hdot_collision. Retain R/A separately.
    """
    if not isinstance(self_energy,RAKSelfEnergy):
        raise ValueError("explicit physical full R/A/K self-energy required")
    expected=hbar/(2*self_energy.energy_unit_J)
    kappa=expected if kappa_s is None else _positive(kappa_s,"kappa_s")
    if not np.isclose(kappa,expected,rtol=4e-15,atol=0):
        raise ValueError("kappa_s must be hbar/(2*energy_unit_J)")
    if retarded.value.shape!=self_energy.retarded_J.value.shape:
        raise ValueError("spectrum and integrated self-energy shapes must agree")
    area=_array(area_weights,"area_weights")
    if area.shape!=retarded.value.shape[:-2] or np.any(area<=0):
        raise ValueError("one positive graph mass per target required")
    K,A=ad.keldysh(retarded,distribution,kappa),ad.advanced(retarded)
    sr,sa,sk=self_energy.retarded_J,self_energy.advanced_J,self_energy.keldysh_J
    rate=(ad.star(retarded,sk,kappa)+ad.star(K,sa,kappa)
        -ad.star(sr,K,kappa)-ad.star(sk,A,kappa)).scale(1j/hbar)
    return RAKCollisionSource(rate.scale(-kappa*area),rate,self_energy,
        "Full Nambu C=-hbar/(2Eu)*m*I_K; SAME Sigma enters spectral R/A; no eta bath")


@dataclass(frozen=True)
class WeakBornFirstOrder:
    """Born counts as first slow order; input spectra/occupations zeroth order."""
    spectral_rhs_first: np.ndarray
    kinetic_rhs_first: np.ndarray
    projected_kinetic_rhs_first: np.ndarray
    max_self_energy_over_Eu: float
    approximation: str


def weak_born_first_order(retarded,distribution,self_energy,area_weights):
    """Return RHS additions for the first spectral and kinetic equations.

    Spectral: L_R R1 = existing_rhs + m/2*[i Sigma_R/Eu,R0].
    Kinetic: existing_residual = C with zeroth-order full R/A/K formula.
    Star gradients of Sigma (already first order) are second order here.
    The self-energy is NOT placed in R0 as an arbitrary eta regulator.
    """
    if not isinstance(self_energy,RAKSelfEnergy):
        raise ValueError("explicit physical R/A/K self-energy required")
    zeros=lambda j:ad.Jet.of(j.value)
    sigma=RAKSelfEnergy(*(zeros(j) for j in (self_energy.retarded_J,
        self_energy.advanced_J,self_energy.keldysh_J)),
        self_energy.energy_unit_J,self_energy.approximation)
    R,H=zeros(retarded),zeros(distribution)
    result=collision_source(R,H,sigma,area_weights)
    spectral=ad.commutator(sigma.retarded_generator,R).value*np.asarray(area_weights)[...,None,None]/2
    size=max(float(np.max(abs(j.value),initial=0)) for j in
        (sigma.retarded_J,sigma.advanced_J,sigma.keldysh_J))/sigma.energy_unit_J
    return WeakBornFirstOrder(spectral,result.source.value,
        ad.projections(result.source.value),size,
        "Explicit weak Born order: Sigma[R0,K0,n0] is O(epsilon), "
        "Moyal Sigma derivatives and Sigma*R1 are omitted O(epsilon^2); "
        "instantaneous supplied phonon snapshot; no total energy/evolution admission")


class GeneralPhononFeedbackRequired(NotImplementedError):
    pass


def require_general_phonon_feedback(self_energy):
    if not isinstance(self_energy,RAKSelfEnergy) or not self_energy.phonon_feedback_complete:
        raise GeneralPhononFeedbackRequired(
            "General superconducting phonon bubble and paired conserving quadrature "
            "required; neutral dp/dn or residual heat cannot supply it.")


def require_general_bgk_self_energy(*,charge_imbalance,moving_spectrum):
    if np.any(np.asarray(charge_imbalance)!=0) or bool(moving_spectrum):
        raise GeneralPhononFeedbackRequired(
            "Neutral frozen-spectrum BGK is not a general conserving Nambu R/A/K self-energy.")


@dataclass(frozen=True)
class PhononPolarizationRHS:
    """Independent phonon collision evaluated from the electron bubble.

    dn/dt has shape (phonon,*node), per s. F_integral has units J.
    The power density is in W/m^3 after the supplied phonon quadrature.
    There is no total-energy residual or adjustment anywhere in this formula.
    """
    dn_dt_per_s: np.ndarray
    bubble_integral_J: np.ndarray
    phonon_power_density_W_per_m3: np.ndarray
    source: str="https://arxiv.org/abs/2011.04471"
    approximation: str=(
        "Zeroth-gradient superconducting polarization bubble paired with weak "
        "local Born alpha2F; signed electronic energies; no dispersive correction")


def phonon_polarization_rhs(retarded_high,distribution_high,
        retarded_low,distribution_low,*,electron_weights_J,phonon_energy_J,
        phonon_weights_J,alpha2F,g_ph_per_J_m3,N0_per_J_m3,
        phonon_occupation,energy_unit_J):
    """Full Nambu bubble, at E_high-E_low=Omega, not a neutral dp/dn adapter.

    Primary universal matrix factor: Nikolic/Basko/Belzig PRB102,214514(2020),
    Eqs33,36, https://arxiv.org/abs/2011.04471. For Q=gR-gA and b=2n+1,
      F=ReTr{b*(Qhigh*Klow-Khigh*Qlow)-Khigh*Klow+Qhigh*Qlow}.
    Its coefficient is derived from the SAME Born self-energy above:
      ndot=pi*N0*alpha2F/(8*hbar*gph)*integral_signed_E(F dE).
    The local vertex is our acoustic/Eliashberg hypothesis; the paper's dirty
    angular Y(q*l) prefactor is not silently substituted for alpha2F.

    Arrays have shape (phonon,electron,*node,2,2). Signed E integration and
    ACTUAL high/low shifted spectral queries must be provided. Finite-cut Ward,
    positivity, and quadrature accuracy remain diagnostics of the caller.
    """
    jets=(retarded_high,distribution_high,retarded_low,distribution_low)
    if any(not isinstance(j,ad.Jet) for j in jets):
        raise ValueError("full shifted Nambu Jets required for phonon bubble")
    shape=jets[0].value.shape
    if len(shape)<4 or any(j.value.shape!=shape for j in jets):
        raise ValueError("matching (phonon,electron,*node,2,2) Jets required")
    omega=_array(phonon_energy_J,"phonon_energy_J")
    wp=_array(phonon_weights_J,"phonon_weights_J")
    we=_array(electron_weights_J,"electron_weights_J")
    alpha=_array(alpha2F,"alpha2F",True)
    dos=_array(g_ph_per_J_m3,"g_ph_per_J_m3")
    if (omega.shape!=(shape[0],) or wp.shape!=omega.shape or alpha.shape!=omega.shape
            or dos.shape!=omega.shape or we.shape!=(shape[1],)
            or np.any(omega<=0) or np.any(wp<=0) or np.any(we<=0) or np.any(dos<=0)):
        raise ValueError("positive physical paired electron/phonon quadratures and phonon DOS required")
    _positive(energy_unit_J,"energy_unit_J")
    N0=_positive(N0_per_J_m3,"N0_per_J_m3")
    n=ElectronPhononBornQuadrature._field(phonon_occupation,
        (shape[0],)+shape[2:-2],"phonon_occupation",True)
    b=(2*n+1)[:,None,...,None,None]
    Rh,Rl=retarded_high.value,retarded_low.value
    Ah,Al=ad.advanced(retarded_high).value,ad.advanced(retarded_low).value
    Hh,Hl=distribution_high.value,distribution_low.value
    Kh,Kl=Rh@Hh-Hh@Ah,Rl@Hl-Hl@Al
    Qh,Ql=Rh-Ah,Rl-Al
    F=np.trace(b*(Qh@Kl-Kh@Ql)-Kh@Kl+Qh@Ql,axis1=-2,axis2=-1).real
    integral=np.tensordot(we,F,axes=(0,1))
    # tensordot puts phonon first among the surviving axes of F.
    prefactor=np.pi*N0*alpha/(8*hbar*dos)
    ndot=integral*prefactor.reshape((shape[0],)+(1,)*(integral.ndim-1))
    power_weights=wp*dos*omega
    power=np.tensordot(power_weights,ndot,axes=(0,0))
    if any(np.any(~np.isfinite(v)) for v in (ndot,integral,power)):
        raise FloatingPointError("nonfinite phonon bubble; no rate/energy repair")
    return PhononPolarizationRHS(ndot,integral,power)


@dataclass(frozen=True)
class BornMaterialSnapshot:
    """One physically paired weak-Born material provider for a signed E grid.

    Complete means all local R/A/K and phonon equations are evaluated from
    one alpha2F, n, and quadrature. It does NOT mean the finite energy cut has
    passed Ward/convergence checks or a temporal NbN detector is admitted.
    """
    retarded_generator_first: np.ndarray
    full_nambu_source: np.ndarray
    projected_source: np.ndarray
    phonon_rate_per_s: np.ndarray
    self_energy: RAKSelfEnergy
    electron_collision_power_density_W_per_m3: np.ndarray
    phonon_collision_power_density_W_per_m3: np.ndarray
    collision_energy_defect_W_per_m3: np.ndarray
    signed_charge_collision_moment_J_per_s: np.ndarray
    max_self_energy_over_Eu: float
    shift_error_max_J: float
    source_contract: str
    general_material_provider_complete: bool=True
    energy_admission: bool=False
    temporal_admission: bool=False


def _zeroth(jet):
    if not isinstance(jet,ad.Jet):
        raise ValueError("full Nambu Jet required")
    return ad.Jet.of(jet.value)


def _validate_distribution_value(jet,name):
    value=jet.value
    off=value.copy()
    off[...,0,0]=0;off[...,1,1]=0
    if np.any(off!=0) or np.any(value.imag!=0):
        raise ValueError(name+" must be real diagonal hL I+hT tau3")
    if np.any(abs(value.real)>1+32*np.finfo(float).eps):
        raise ValueError(name+" leaves Pauli support; no clipping applied")


def evaluate_born_material_snapshot(*,electron_energy_J,electron_weights_J,
        phonon_energy_J,phonon_weights_J,alpha2F,g_ph_per_J_m3,
        N0_per_J_m3,energy_unit_J,retarded,distribution,
        retarded_plus,distribution_plus,retarded_minus,distribution_minus,
        shifted_energy_plus_J,shifted_energy_minus_J,phonon_occupation,
        area_weights,shift_contract):
    """Evaluate R/A/K AND phonon collision from one explicit physical state.

    E and wE are a signed symmetric quadrature. R,H have shape
    (electron,node,2,2). Shifted Jets have (phonon,electron,node,2,2) shape.
    The supplied energy-query tags MUST equal E+/-Omega at machine precision;
    no extrapolation/interpolation is constructed. Both signs and the same
    wE are reused for the independent phonon bubble and collision moments.
    n has (phonon,node) shape or a spatially uniform (phonon,) shape.

    Output sources are in the ORIGINAL signed E order; a positive-energy
    joint solver may select E>0 afterwards, but must not halve the bubble.
    Sigma is O(Born), so R0/H0 are used throughout and no first-order Sigma
    derivatives are mixed into that order. The general e-e BGK closure is not
    silently enabled. Rates/balance errors are returned without repair.
    """
    from dataclasses import replace
    E=_array(electron_energy_J,"electron_energy_J")
    we=_array(electron_weights_J,"electron_weights_J")
    Eu=_positive(energy_unit_J,"energy_unit_J")
    omega=_array(phonon_energy_J,"phonon_energy_J")
    wp=_array(phonon_weights_J,"phonon_weights_J")
    alpha=_array(alpha2F,"alpha2F",True)
    dos=_array(g_ph_per_J_m3,"g_ph_per_J_m3")
    N0=_positive(N0_per_J_m3,"N0_per_J_m3")
    if (E.ndim!=1 or len(E)<2 or np.any(np.diff(E)<=0) or E[0]>=0 or E[-1]<=0
            or we.shape!=E.shape or np.any(we<=0)):
        raise ValueError("ordered signed energies and matching positive electron weights required")
    energy_scale=max(Eu,float(np.max(abs(E))))
    if (not np.allclose(E,-E[::-1],rtol=32*np.finfo(float).eps,atol=32*np.finfo(float).eps*Eu)
            or not np.allclose(we,we[::-1],rtol=32*np.finfo(float).eps,atol=0)):
        raise ValueError("signed electron quadrature must have symmetric E and weights")
    if (omega.ndim!=1 or not len(omega) or np.any(omega<=0)
            or np.any(np.diff(omega)<=0) or wp.shape!=omega.shape or alpha.shape!=omega.shape
            or dos.shape!=omega.shape or np.any(wp<=0) or np.any(dos<=0)):
        raise ValueError("matching positive phonon quadrature and DOS required")
    if not isinstance(shift_contract,str) or not shift_contract.strip():
        raise ValueError("explicit shifted spectral query provenance required")
    plus=_array(shifted_energy_plus_J,"shifted_energy_plus_J")
    minus=_array(shifted_energy_minus_J,"shifted_energy_minus_J")
    expected_plus=E[None,:]+omega[:,None]
    expected_minus=E[None,:]-omega[:,None]
    if plus.shape!=expected_plus.shape or minus.shape!=expected_minus.shape:
        raise ValueError("one physical energy tag per phonon/electron shifted query required")
    shift_error=max(float(np.max(abs(plus-expected_plus))),float(np.max(abs(minus-expected_minus))))
    tolerance=32*np.finfo(float).eps*max(energy_scale,float(np.max(omega)))
    if shift_error>tolerance:
        raise ValueError("shifted queries do not equal E +/- Omega; no endpoint/interpolation repair")
    base=(_zeroth(retarded),_zeroth(distribution))
    shifted=tuple(_zeroth(j) for j in (
        retarded_plus,distribution_plus,retarded_minus,distribution_minus))
    R,H=base;Rp,Hp,Rm,Hm=shifted
    if (R.value.ndim!=4 or R.value.shape[0]!=len(E) or H.value.shape!=R.value.shape
            or R.value.shape[1]<1):
        raise ValueError("base spectra/distributions need (electron,node,2,2) shape")
    nn=R.value.shape[1]
    shape=(len(omega),len(E),nn,2,2)
    if any(j.value.shape!=shape for j in shifted):
        raise ValueError("shifted spectra/distributions need (phonon,electron,node,2,2) shape")
    for field,name in ((H,"distribution"),(Hp,"distribution_plus"),(Hm,"distribution_minus")):
        _validate_distribution_value(field,name)
    mass=_array(area_weights,"area_weights")
    if mass.shape!=(nn,) or np.any(mass<=0):
        raise ValueError("one positive graph area/mass per spatial node required")
    n=ElectronPhononBornQuadrature._field(
        phonon_occupation,(len(omega),nn),"phonon_occupation",True)
    operator=ElectronPhononBornQuadrature(omega,wp,alpha,
        energy_unit_J=Eu,shift_contract=shift_contract)
    sigma=operator.evaluate(Rp,Hp,Rm,Hm,n[:,None,:])
    mass_grid=np.broadcast_to(mass,R.value.shape[:-2])
    weak=weak_born_first_order(R,H,sigma,mass_grid)
    Rlow=ad.Jet.of(np.broadcast_to(R.value,shape))
    Hlow=ad.Jet.of(np.broadcast_to(H.value,shape))
    bubble=phonon_polarization_rhs(Rp,Hp,Rlow,Hlow,
        electron_weights_J=we,phonon_energy_J=omega,phonon_weights_J=wp,
        alpha2F=alpha,g_ph_per_J_m3=dos,N0_per_J_m3=N0,
        phonon_occupation=n,energy_unit_J=Eu)
    sigma=replace(sigma,phonon_feedback_complete=True)
    # Independent electron moment of I_K; never infer phonon rate from it.
    rate=collision_source(R,H,sigma,mass_grid).collision_rate_per_s.value
    projected_rate=ad.projections(rate)
    electron_power=-N0*np.tensordot(we*E,projected_rate[...,0],axes=(0,0))
    phonon_power=np.asarray(bubble.phonon_power_density_W_per_m3)
    charge=np.tensordot(we,projected_rate[...,1],axes=(0,0))
    return BornMaterialSnapshot(sigma.retarded_generator.value,
        weak.kinetic_rhs_first,weak.projected_kinetic_rhs_first,
        bubble.dn_dt_per_s,sigma,electron_power,phonon_power,
        electron_power+phonon_power,charge,weak.max_self_energy_over_Eu,
        shift_error,
        "Paired weak-Born local material: full Nambu R/A/K and independently "
        "evaluated phonon polarization, identical alpha2F/n/signed quadrature; "
        "sources at first order, no eta bath/BGK/residual heat; "
        "finite-cut energy/charge defects remain diagnostics; "+shift_contract)
