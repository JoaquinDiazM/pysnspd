"""User-authorized effective NbN phonons, anchored to acoustic heat capacity.

This is an EXPLORATORY material contract, not recovered DFPT metadata or a
calibration of the Korzh film. An arbitrary phonon-shape ordinate obtains its
volumetric scale from Cph(Tref)=B*Tref^3, with a primary-source sound speed.
alpha2F retains the declared dimensionless Eliashberg convention; no fitted
rate multiplier is introduced. The chosen reference temperature and shapes
remain explicit model choices, and their sensitivity is measured separately.
"""
from __future__ import annotations
from dataclasses import dataclass
import hashlib
import json
from pathlib import Path
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.constants import Boltzmann as kB, elementary_charge as e, hbar, Planck


VODOLAZOV_SOURCE='https://arxiv.org/pdf/1611.06060'
ALLMARAS_SOURCE='https://thesis.caltech.edu/13748/08/Allmaras_Thesis_Final.pdf'
QNN_SOURCE='https://github.com/qnngroup/proj-KE-solver/tree/5b6bd747f80016da5ccd51db73c110a8ecc6abf6'


def _positive(value,name):
    raw=np.asarray(value)
    if raw.ndim or np.iscomplexobj(raw) or raw.dtype.kind=='b' or not np.isfinite(raw) or raw<=0:
        raise ValueError(name+' must be a finite positive real scalar')
    return float(raw)


def _axis(value,name,allow_zero=False):
    raw=np.asarray(value)
    if np.iscomplexobj(raw) or raw.dtype.kind=='b':raise ValueError(name+' must be real')
    out=np.asarray(raw,float)
    if (out.ndim!=1 or len(out)<2 or np.any(~np.isfinite(out))
            or np.any(out<0 if allow_zero else out<=0) or np.any(np.diff(out)<=0)):
        raise ValueError(name+' must be ordered, finite and nonnegative' if allow_zero else name+' must be ordered, finite and positive')
    return out.copy()


def acoustic_coefficient_J_per_m3_K4(sound_speed_m_per_s):
    speed=_positive(sound_speed_m_per_s,'sound_speed_m_per_s')
    return 2*np.pi**2*kB**4/(5*hbar**3*speed**3)


def _bose_cv_kernel(x):
    result=np.zeros_like(x)
    small=x<1e-4
    result[small]=1-x[small]**2/12+x[small]**4/240
    regular=(~small)&(x<700)
    # x^2 exp(x)/(exp(x)-1)^2 with no large exponential or zero division.
    result[regular]=(x[regular]/(2*np.sinh(x[regular]/2)))**2
    return result


@dataclass(frozen=True)
class EffectivePhononMaterial:
    energy_J: np.ndarray
    alpha2F: np.ndarray
    g_ph_per_J_m3: np.ndarray
    N0_per_J_m3: float
    manifest: dict

    def __post_init__(self):
        energy=_axis(self.energy_J,'energy_J',True)
        alpha=np.asarray(self.alpha2F,float).copy();dos=np.asarray(self.g_ph_per_J_m3,float).copy()
        if (alpha.shape!=energy.shape or dos.shape!=energy.shape or np.any(~np.isfinite(alpha))
                or np.any(~np.isfinite(dos)) or np.any(alpha<0) or np.any(dos<0)
                or np.any((dos==0)&(alpha>0))):
            raise ValueError('Finite nonnegative compatible alpha2F and volumetric DOS required')
        _positive(self.N0_per_J_m3,'N0_per_J_m3')
        for array in (energy,alpha,dos):array.setflags(write=False)
        object.__setattr__(self,'energy_J',energy);object.__setattr__(self,'alpha2F',alpha)
        object.__setattr__(self,'g_ph_per_J_m3',dos)

    def _quadrature(self):
        # Integrate the declared piecewise-linear DOS against smooth thermal
        # kernels; using the native nodes' trapezoid is biased at low T.
        x,w=leggauss(8);width=np.diff(self.energy_J)
        energies=((self.energy_J[:-1,None]+self.energy_J[1:,None])/2+width[:,None]*x/2).ravel()
        weights=(width[:,None]*w/2).ravel()
        fraction=(x+1)/2
        density=(self.g_ph_per_J_m3[:-1,None]*(1-fraction)+self.g_ph_per_J_m3[1:,None]*fraction).ravel()
        if self.manifest.get('analytic_acoustic_Debye'):
            density=self.manifest['g_ph_coefficient_per_J3_m3']*energies**2
        return energies,weights,density

    def heat_capacity_J_per_m3_K(self,temperature_K):
        temperature=_positive(temperature_K,'temperature_K')
        energy,weights,density=self._quadrature()
        return kB*float(np.dot(weights,density*_bose_cv_kernel(energy/(kB*temperature))))

    def energy_density_J_per_m3(self,temperature_K):
        raw=np.asarray(temperature_K)
        if raw.ndim or np.iscomplexobj(raw) or raw.dtype.kind=='b' or not np.isfinite(raw) or raw<0:
            raise ValueError('temperature_K must be finite real nonnegative')
        if raw==0:return 0.
        energy,weights,density=self._quadrature();x=energy/(kB*float(raw));factor=np.zeros_like(x)
        small=x<1e-4;factor[small]=1-x[small]/2+x[small]**2/12-x[small]**4/720
        regular=(~small)&(x<700);factor[regular]=x[regular]/np.expm1(x[regular])
        return kB*float(raw)*float(np.dot(weights,density*factor))

    def evaluate(self,energy_J):
        raw=np.asarray(energy_J)
        if np.iscomplexobj(raw) or raw.dtype.kind=='b':raise ValueError('Queries must be real energies')
        energy=np.asarray(raw,float)
        if np.any(~np.isfinite(energy)) or np.any(energy<self.energy_J[0]) or np.any(energy>self.energy_J[-1]):
            raise ValueError('Query outside explicit material energy support; no extrapolation')
        alpha=np.interp(energy,self.energy_J,self.alpha2F)
        dos=np.interp(energy,self.energy_J,self.g_ph_per_J_m3)
        if self.manifest.get('analytic_acoustic_Debye'):
            alpha=self.manifest['alpha2F_coefficient_per_J2']*energy**2
            dos=self.manifest['g_ph_coefficient_per_J3_m3']*energy**2
        coupling=np.divide(alpha,dos,out=np.zeros_like(alpha),where=dos>0)
        return alpha,dos,coupling

    def phonon_grid(self,energy_J,weights_J,delta0_J):
        """Convert an explicit physical phonon quadrature to the native events.

        Capacities are gph/N0 * weights_J/Delta0_J. No event-rate fit, cutoff
        extrapolation, empty-mode filling or phonon preparation is performed.
        The event builder still rejects any coupled event lacking brackets.
        """
        from .kinetic_events import PhononGrid
        energy=_axis(energy_J,'phonon quadrature energy_J');weights=np.asarray(weights_J,float)
        delta=_positive(delta0_J,'delta0_J')
        if weights.shape!=energy.shape or np.any(~np.isfinite(weights)) or np.any(weights<=0):
            raise ValueError('Positive physical quadrature weights required')
        _,dos,_=self.evaluate(energy)
        if np.any(dos<=0):raise ValueError('Phonon quadrature includes an empty physical mode; supply its active support explicitly')
        capacities=dos/self.N0_per_J_m3*weights/delta
        return PhononGrid(energy/delta,capacities,
            lambda omega:self.evaluate(np.asarray(omega)*delta)[0],
            (self.energy_J[0]/delta,self.energy_J[-1]/delta),
            self.manifest['material_id']+'; acoustic anchored effective full shape; not film calibration')

    @staticmethod
    def event_rate_prefactor(delta0_J,reference_time_s):
        # D.15 event convention; dividing the resulting RHS by tref gives s^-1.
        return 8*np.pi*_positive(delta0_J,'delta0_J')*_positive(reference_time_s,'reference_time_s')/hbar

    @classmethod
    def load(cls,npz_path):
        with np.load(npz_path,allow_pickle=False) as z:
            return cls(z['energy_J'],z['alpha2F'],z['g_ph_per_J_m3'],float(z['N0_per_J_m3']),json.loads(str(z['manifest_json'])))

    def save(self,npz_path,manifest_path):
        npz_path=Path(npz_path);manifest_path=Path(manifest_path)
        with npz_path.open('xb') as stream:
            np.savez_compressed(stream,energy_J=self.energy_J,alpha2F=self.alpha2F,
                g_ph_per_J_m3=self.g_ph_per_J_m3,N0_per_J_m3=self.N0_per_J_m3,
                manifest_json=np.array(json.dumps(self.manifest,sort_keys=True,allow_nan=False)))
        with manifest_path.open('x',encoding='utf8') as stream:
            json.dump(self.manifest,stream,indent=2,allow_nan=False);stream.write('\n')


def acoustic_anchored_material(energy_J,alpha2F,phonon_shape,*,reference_temperature_K,
                               sound_speed_m_per_s,N0_per_J_m3,shape_contract,alpha2F_contract,
                               source_sha256=None,material_id='NbN_effective_acoustic_20261008'):
    energy=_axis(energy_J,'energy_J',True);alpha=np.asarray(alpha2F,float);shape=np.asarray(phonon_shape,float)
    for contract in (shape_contract,alpha2F_contract):
        if not isinstance(contract,str) or not contract.strip():raise ValueError('Explicit shape and alpha2F convention contracts required')
    temperature=_positive(reference_temperature_K,'reference_temperature_K');N0=_positive(N0_per_J_m3,'N0_per_J_m3')
    speed=_positive(sound_speed_m_per_s,'sound_speed_m_per_s');coefficient=acoustic_coefficient_J_per_m3_K4(speed)
    # The unknown original ordinate unit disappears after the heat-capacity
    # normalization; it is not declared to be states/THz or states/meV.
    temporary=EffectivePhononMaterial(energy,alpha,shape,N0,{})
    shape_capacity=temporary.heat_capacity_J_per_m3_K(temperature)
    if not np.isfinite(shape_capacity) or shape_capacity<=0:raise ValueError('Shape carries no finite thermal heat capacity')
    scale=coefficient*temperature**3/shape_capacity
    manifest=dict(schema='pysnspd.effective_phonon_material.v1',status='EFFECTIVE_EXPLORATORY_MATERIAL',
        material_id=material_id,reference_temperature_K=temperature,sound_speed_m_per_s=speed,
        acoustic_B_J_per_m3_K4=coefficient,shape_capacity_at_reference_in_arbitrary_units=shape_capacity,
        shape_to_gph_multiplier=scale,shape_contract=shape_contract,alpha2F_contract=alpha2F_contract,
        source_sha256=source_sha256,N0_per_J_m3=N0,
        units=dict(energy='J',g_ph='1/(J m^3)',alpha2F='dimensionless Eliashberg',alpha2='J m^3',N0='single spin / (J m^3)'),
        sound_speed_source=ALLMARAS_SOURCE,sound_speed_locator='Allmaras2020 material parameter discussion, printed80-81; heat-capacity equations printed13-16',
        user_authorization='2026-10-08: material effective for exploratory tests',
        calibration='C_ph(Tref)=2*pi^2*kB^4/(5*hbar^3*vavg^3)*Tref^3',
        mode_count_normalized=False,DFPT_metadata_recovered=False,Korzh_film_rates_calibrated=False,
        alpha2F_rate_multiplier=1.,rate_prefactor='8*pi*Delta0*tref/hbar for shared conservative S/R events; no arbitrary rate fitted',
        phonon_feedback_convention='g_ph fixes occupation-energy conversion and N0*alpha2F/g_ph feedback together',
        limitation='Acoustic heat-capacity anchor at one temperature, extended full-band shape is exploratory; no measured confidence interval or optical-DOS calibration')
    material=EffectivePhononMaterial(energy,alpha,shape*scale,N0,manifest)
    manifest['heat_capacity_at_reference_J_per_m3_K']=material.heat_capacity_J_per_m3_K(temperature)
    manifest['target_capacity_at_reference_J_per_m3_K']=coefficient*temperature**3
    manifest['total_modes_per_m3']=float(np.trapezoid(material.g_ph_per_J_m3,energy))
    return material


def from_nbn_shape_csv(path,*,Tc_K=8.65,diffusivity_m2_s=5e-5,sheet_resistance_ohm=608.,
                       thickness_m=7e-9,sound_speed_m_per_s=4912.,expected_sha256=None):
    path=Path(path);raw=path.read_bytes();digest=hashlib.sha256(raw).hexdigest()
    if expected_sha256 is not None and digest!=expected_sha256:raise ValueError('Phonon shape source hash changed')
    data=np.loadtxt(path,delimiter=',',skiprows=1)
    if data.ndim!=2 or data.shape[1]!=4:raise ValueError('Expected existing stage1 derived four-column shape CSV')
    D=_positive(diffusivity_m2_s,'diffusivity_m2_s');R=_positive(sheet_resistance_ohm,'sheet_resistance_ohm');d=_positive(thickness_m,'thickness_m')
    N0=1/(2*e**2*D*R*d)
    material=acoustic_anchored_material(data[:,0]*Planck*1e12,data[:,1],data[:,2],
        reference_temperature_K=Tc_K,sound_speed_m_per_s=sound_speed_m_per_s,N0_per_J_m3=N0,
        shape_contract='Existing stable-first common-support stage1 shape; ordinate is arbitrary, axis is explicitly adopted as ordinary frequency THz under the intended upstream reader. No atom/cell count inferred.',
        alpha2F_contract='Column2 adopted as dimensionless Eliashberg alpha2F of QNN25: lambda=2 integral alpha2F dOmega/Omega; ordinary-frequency change Omega=h nu leaves this ordinate unchanged. Conditional effective model input, not recovered file metadata.',
        source_sha256=digest)
    material.manifest.update(source_path=str(path),shape_primary_source=QNN_SOURCE,
        diffusivity_m2_s=D,sheet_resistance_ohm=R,film_thickness_m=d,
        Einstein_relation='N0=1/(2 e^2 D R_sheet d)',Tc_K=float(Tc_K))
    # lambda quadrature is a diagnostic of the adopted dimensionless A; it is
    # not used to tune the acoustic DOS scale or infer optical lifetimes.
    mask=material.energy_J>0
    material.manifest['lambda_ep']=float(2*np.trapezoid(material.alpha2F[mask]/material.energy_J[mask],material.energy_J[mask]))
    return material



def allmaras_acoustic_material(*,Tc_K=8.65,tau0_s=1.87e-9,Debye_temperature_K=664.,
        sound_speed_m_per_s=4912.,diffusivity_m2_s=5e-5,sheet_resistance_ohm=608.,
        thickness_m=7e-9,nodes=4097):
    """Explicit approved acoustic A20 completion over the full Debye band.

    Omega is physical ENERGY in J. Allmaras' dimensionless variable is
    omega=Omega/(kB*Tc). Matching its normal-state collision coefficient gives
    (2*pi/hbar)*alpha2F(Omega)*dOmega = 2*omega^2*domega/tau0.
    The factor two comes from expanding the FD/Bose brackets in Vodolazov
    (2017) Eq3 / Allmaras Eq2.6. Vodolazov Eq4 retains pair+2*scattering
    phonon feedback, unlike the missing factor in the printed A20 Eq2.7.
    gph is the three-branch Debye DOS per energy and volume. Both absolute
    normal-state rates and occupation feedback follow these same functions.
    This is an effective exploratory A20 acoustic model, not a fit to the
    experimental 608-ohm film or a reconstructed raw optical phonon table.
    """
    Tc=_positive(Tc_K,'Tc_K');tau=_positive(tau0_s,'tau0_s');TD=_positive(Debye_temperature_K,'Debye_temperature_K')
    speed=_positive(sound_speed_m_per_s,'sound_speed_m_per_s');D=_positive(diffusivity_m2_s,'diffusivity_m2_s')
    resistance=_positive(sheet_resistance_ohm,'sheet_resistance_ohm');thick=_positive(thickness_m,'thickness_m')
    if isinstance(nodes,bool) or not isinstance(nodes,(int,np.integer)) or nodes<3:raise ValueError('At least three energy nodes required')
    Eu=kB*Tc;cutoff=kB*TD;N0=1/(2*e**2*D*resistance*thick)
    ga=3/(2*np.pi**2*hbar**3*speed**3);aa=hbar/(np.pi*tau*Eu**3)
    energy=np.linspace(0,cutoff,int(nodes));B=acoustic_coefficient_J_per_m3_K4(speed)
    E0=4*N0*Eu**2;gamma=(4*np.pi**4/15)*E0/(B*Tc**4)
    manifest=dict(schema='pysnspd.effective_phonon_material.v1',status='EFFECTIVE_EXPLORATORY_MATERIAL',
        material_id='NbN_A20_acoustic_V17_conservative_effective_v2_20261008',analytic_acoustic_Debye=True,
        Tc_K=Tc,tau0_s=tau,Debye_temperature_K=TD,Debye_cutoff_J=cutoff,sound_speed_m_per_s=speed,
        diffusivity_m2_s=D,sheet_resistance_ohm=resistance,film_thickness_m=thick,N0_per_J_m3=N0,
        Einstein_relation='N0=1/(2 e^2 D R_sheet d), single-spin DOS',
        acoustic_B_J_per_m3_K4=B,g_ph_coefficient_per_J3_m3=ga,alpha2F_coefficient_per_J2=aa,
        g_ph_equation='3 Omega_J^2/(2*pi^2*hbar^3*v_avg^3), 0<Omega<kB*T_D',
        alpha2F_equation='hbar Omega_J^2/(pi*tau0*(kB*Tc)^3), dimensionless Eliashberg ordinate',
        normal_collision_identity='(2*pi/hbar)*alpha2F(Omega_J)*dOmega_J=(2*omega^2/tau0)*domega; omega=Omega_J/(kB*Tc)',
        event_prefactor='8*pi*Delta0*tref/hbar, native shared event convention with electron capacity 4; per-second RHS divides by tref',
        phonon_feedback_convention='g_ph and the SAME alpha2F determine N0*alpha2F/g_ph, rather than a separately fitted phonon rate',
        collision_primary_source=VODOLAZOV_SOURCE,
        collision_source_locator='Vodolazov2017 equations3-4: electronic brackets expand to2 activities; phonon pair+2scattering',
        A20_printed_collision_discrepancy='Allmaras Eq2.7 printed13 omits factor2 on phonon scattering relative to Vodolazov Eq4; adopt the conservative primary Eq3-4 pair+2scattering, preserving A20 measured tau0 and sound-speed anchor',
        previous_v1_rate_scale_superseded='Earlier 20261008 material export had alpha2F half the tau0-normalized value; gph,N0,Cph,Uph unchanged',
        gamma_A20=gamma,gamma_independent=(8*np.pi**2/3)*N0*hbar**3*speed**3/Eu**2,
        reference_temperature_K=Tc,source=ALLMARAS_SOURCE,
        source_locator='Allmaras 2020 thesis printed13-16, equations2.6-2.8; printed15 material discussion: tau0=1.87ns, vavg=4912m/s, DebyeT=664K',
        source_pdf_sha256='2554e9a2b047f93f7c7f50a2897887afaf522eb530e057307c5a31ef4e1e0e96',
        source_energy_frequency_convention='Omega is energy in J, not ordinary frequency or angular frequency; source collision equations use Omega/(kB Tc)',
        units=dict(energy='J',g_ph='1/(J m^3)',alpha2F='dimensionless Eliashberg',alpha2='J m^3',N0='single spin / (J m^3)'),
        user_authorization='2026-10-08: effective exploratory material and full A20 acoustic band authorized',
        raw_QNN_shape_used=False,raw_QNN_low_energy_silently_completed=False,DFPT_metadata_recovered=False,
        Korzh_film_rates_calibrated=False,alpha2F_rate_multiplier=1.,mode_count_normalized=False,
        limitation='Effective three-dimensional acoustic phonons through Debye cutoff, clean-limit normal-state electron-phonon anchor; thin-film and optical structure uncertainties remain explicit')
    material=EffectivePhononMaterial(energy,aa*energy**2,ga*energy**2,N0,manifest)
    manifest['heat_capacity_at_reference_J_per_m3_K']=material.heat_capacity_J_per_m3_K(Tc)
    manifest['target_capacity_at_reference_J_per_m3_K']=B*Tc**3
    manifest['total_modes_per_m3']=ga*cutoff**3/3
    manifest['lambda_ep']=aa*cutoff**2
    return material
