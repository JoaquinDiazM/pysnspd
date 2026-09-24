"""Full three-state thesis circuit driven by a declared linear device port.

This algebra never replaces the film by an invented resistor or certifies a
device admittance. Phasors use an explicit sign convention and peak, not RMS,
amplitudes. The historical default here is exp(+i omega t); finite-frequency
Usadel currently uses exp(-i omega t).
"""
from dataclasses import dataclass
import numpy as np

from .electrical_ports import ThesisCircuitParameters
from .energy_catalog import HBAR_J_S, E_CHARGE_C


def _finite_scalar(value, name, *, complex_allowed=False):
    raw = np.asarray(value)
    if raw.ndim or not np.isfinite(raw) or (not complex_allowed and np.iscomplexobj(raw)):
        raise ValueError(name + ' must be a finite scalar')
    return complex(raw) if complex_allowed else float(raw)


@dataclass(frozen=True)
class InductancePartition:
    total_reference_H: float
    resolved_differential_H: float
    external_H: float


def partition_inductance(total_reference_H, phase_slope_rad_per_A):
    """Fix Lext from a supplied stable-branch d(thetaR-thetaL)/dI.

    The caller must establish branch stability and use the same voltage port
    planes. This function checks positivity, not the physical provenance.
    """
    total = _finite_scalar(total_reference_H, 'total_reference_H')
    slope = _finite_scalar(phase_slope_rad_per_A, 'phase_slope_rad_per_A')
    resolved = HBAR_J_S/(2*E_CHARGE_C)*slope
    external = total-resolved
    if total <= 0 or resolved <= 0 or external <= 0:
        raise ValueError('Positive resolved and exterior inductance required; no clipping')
    return InductancePartition(total, resolved, external)


@dataclass(frozen=True)
class HarmonicCircuitResponse:
    omega_rad_per_s: float
    bias_voltage_peak_V: complex
    device_impedance_ohm: complex
    state_peak: np.ndarray
    device_voltage_peak_V: complex
    readout_voltage_peak_V: complex
    rf_current_peak_A: complex
    equation_residual: np.ndarray
    source_complex_power_VA: complex
    device_complex_power_VA: complex
    bias_loss_W: float
    load_loss_W: float
    external_reactive_power_var: float
    complex_power_residual_VA: complex
    average_power_residual_W: float
    time_convention: str


def solve_harmonic(parameters, omega_rad_per_s, device_impedance_ohm, *, bias_voltage_peak_V=1.,
                   time_convention='exp(+iwt)'):
    """Solve CM.4 without shortening any physical circuit time constant.

    device_impedance is supplied by an independent film calculation, in the
    passive orientation delta Vdev = Zdev delta Is. A negative real part is
    retained and reported, not clipped: a biased device may be active.
    The DC operating point and its powers are not part of these AC increments.
    Zdev and the forcing must use the declared time_convention. No implicit
    conjugation converts a device result from the opposite convention.
    """
    if not isinstance(parameters, ThesisCircuitParameters):
        raise TypeError('Full ThesisCircuitParameters required')
    omega = _finite_scalar(omega_rad_per_s, 'omega_rad_per_s')
    impedance = _finite_scalar(device_impedance_ohm, 'device_impedance_ohm', complex_allowed=True)
    drive = _finite_scalar(bias_voltage_peak_V, 'bias_voltage_peak_V', complex_allowed=True)
    if omega <= 0:
        raise ValueError('Positive angular frequency required for an AC average')
    if time_convention not in ('exp(+iwt)', 'exp(-iwt)'):
        raise ValueError('Explicit exp(+iwt) or exp(-iwt) time convention required')
    sign = 1. if time_convention == 'exp(+iwt)' else -1.
    p = parameters
    matrix = np.array([
        [sign*1j*omega*p.L_bias_H+p.R_bias_ohm+p.R_load_ohm, -p.R_load_ohm, 1.],
        [-p.R_load_ohm, sign*1j*omega*p.Lk_ext_H+p.R_load_ohm+impedance, -1.],
        [-1., 1., sign*1j*omega*p.C_couple_F]], complex)
    forcing = np.array([drive, 0., 0.], complex)
    state = np.linalg.solve(matrix, forcing)
    if np.any(~np.isfinite(state)):
        raise FloatingPointError('Nonfinite circuit response; no fallback impedance used')
    ib, detector, vc = state
    rf = ib-detector
    vdev = impedance*detector
    vout = p.R_load_ohm*rf
    source = .5*drive*np.conj(ib)
    device = .5*vdev*np.conj(detector)
    bias_loss = .5*p.R_bias_ohm*abs(ib)**2
    load_loss = .5*p.R_load_ohm*abs(rf)**2
    reactive = sign*.5*omega*(p.L_bias_H*abs(ib)**2+p.Lk_ext_H*abs(detector)**2-p.C_couple_F*abs(vc)**2)
    residual = source-device-bias_loss-load_loss-1j*reactive
    return HarmonicCircuitResponse(omega, drive, impedance, state, vdev, vout, rf,
        matrix@state-forcing, complex(source), complex(device), float(bias_loss),
        float(load_loss), float(reactive), complex(residual), float(residual.real), time_convention)
