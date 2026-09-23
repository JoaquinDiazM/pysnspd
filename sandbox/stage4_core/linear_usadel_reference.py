"""Independent thermal linear-response reference, C.29--C.31; no trajectory.

Uses only NumPy/SciPy, not the experimental catalogue or spatial operators.
The Matsubara sums include analytic large-frequency tails. Fixed T is essential:
these Hessians are NOT the frozen-population Hessian of a nonequilibrium state.
"""
from __future__ import annotations

import argparse
import json
import math
import time
from pathlib import Path

import numpy as np
from scipy.constants import hbar, Boltzmann
from scipy.optimize import brentq
from scipy.special import zeta


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_OUTPUT = ROOT / 'docs/implementation/stage4/review_20260923/linear_usadel_reference.json'


def frequency_sum_tail(power, count, t):
    """Sum epsilon_n**(-power), n=count..infinity, epsilon=2*pi*t*(n+1/2)."""
    return float(zeta(power, count + .5) / (2 * np.pi * t) ** power)


def gap_equation(d, count, t, *, include_tail=True):
    eps = 2 * np.pi * t * (np.arange(count) + .5)
    energy = np.hypot(eps, d)
    # Algebraically stable 1/eps - 1/E, without subtracting two near-equal terms.
    total = float(np.sum(d*d / (eps * energy * (eps + energy))))
    if include_tail:
        total += (d*d/2 * frequency_sum_tail(3, count, t)
                  - 3*d**4/8 * frequency_sum_tail(5, count, t)
                  + 5*d**6/16 * frequency_sum_tail(7, count, t))
    return float(math.log(t) + 2*np.pi*t*total)


def amplitude_stiffness(d, kbar, count, t):
    eps = 2*np.pi*t*(np.arange(count)+.5)
    energy = np.hypot(eps, d)
    a, s = d*d, kbar*kbar
    h0_raw = float(4*np.pi*t*np.sum(a/energy**3))
    h0_tail = 4*np.pi*t*(
        a*frequency_sum_tail(3, count, t)
        - 1.5*a*a*frequency_sum_tail(5, count, t)
        + 15*a**3/8*frequency_sum_tail(7, count, t))
    increment_raw = float(4*np.pi*t*np.sum((eps/energy)**2*s/(energy*(energy+s))))
    # s*eps**2/[E**3*(E+s)] expanded through eps**(-6).
    coefficients = {2: s, 3: -s*s, 4: s**3-2*s*a,
                    5: -s**4+2.5*s*s*a,
                    6: s**5-3*s**3*a+3*s*a*a}
    increment_tail = 4*np.pi*t*sum(
        coefficient*frequency_sum_tail(power, count, t)
        for power, coefficient in coefficients.items())
    hu = h0_raw + h0_tail + increment_raw + increment_tail
    hk = h0_raw + h0_tail + np.pi/2*s
    return dict(k_ell0=float(kbar), H_U0_raw=h0_raw, H_U0_tail=float(h0_tail),
                H_U0=float(h0_raw+h0_tail), increment_U_raw=increment_raw,
                increment_U_tail=float(increment_tail), H_U=float(hu),
                H_K0=float(hk), stiffness_excess_fraction=float(hk/hu-1),
                conditional_time_K0_over_U=float(hu/hk))


def calculate():
    start = time.monotonic()
    tc, temperature, diffusion = 8.65, .9, 5e-5
    t = temperature/tc
    ell0 = float(np.sqrt(hbar*diffusion/(2*Boltzmann*tc)))
    counts = (2500, 5000)
    grids = []
    for count in counts:
        d = float(brentq(lambda d: gap_equation(d, count, t), .01, 4.))
        rows = [amplitude_stiffness(d, k, count, t) for k in (0., .15, .3, .5, 1.)]
        for row in rows:
            k = row['k_ell0']
            row['spatial_wavelength_nm'] = None if k == 0 else float(2*np.pi*ell0/k*1e9)
        grids.append(dict(matsubara_count=count, gap_kBTc=d,
                          gap_equation_residual=gap_equation(d, count, t), rows=rows))
    low, high = grids
    differences = [abs(a['H_U']-b['H_U'])/abs(b['H_U'])
                   for a, b in zip(low['rows'], high['rows'])]
    return dict(schema='pysnspd.stage4.thermal_linear_usadel_reference.v1',
        status='COMPUTED_STATIC_REFERENCE_NOT_CORE_OR_TRAJECTORY_ADMISSION',
        source='docs/modelo_v0_4/C_condensado_corriente_y_senal_v0_4.md C.29-C.31',
        independent_of_experimental_catalogue_and_spatial_assembly=True,
        physical_time_steps=0, new_photon_runs=0,
        ensemble='Fixed imposed temperature, thermal BCS stationary gap, zero background current; not fixed occupation counts',
        parameters=dict(Tc_K=tc, T_K=temperature, D_m2_s=diffusion,
                        ell0_nm=ell0*1e9, T_over_Tc=t,
                        gap_units='kB*Tc', zero_temperature_weak_coupling_ratio=float(np.pi*np.exp(-np.euler_gamma))),
        definitions=dict(H='Derivative of force normalized by N0*kB*Tc with respect to d=abs(Delta)/(kB*Tc)',
            HU='Spatial-linearized Usadel response after thermal spectral relaxation',
            HK0='Same homogeneous stiffness with the local K0 gradient completion',
            time_ratio='t_K0/t_U=H_U/H_K0 only if identical positive mobility is imposed',
            wavelength='2*pi/k, not 1/k; must additionally satisfy diffusive k*mean_free_path << 1'),
        tail_method=dict(gap='Inverse odd frequency powers 3,5,7 with Hurwitz zeta',
            H0='Inverse odd frequency powers 3,5,7 with Hurwitz zeta',
            HU_increment='Inverse frequency powers 2..6 with Hurwitz zeta',
            interpretation='A cutoff comparison checks numerical stability, not physical uncertainty'),
        cutoff_runs=grids,
        cutoff_comparison=dict(maximum_corrected_HU_relative_difference=max(differences),
                               gap_absolute_difference=abs(low['gap_kBTc']-high['gap_kBTc'])),
        results=high['rows'],
        restrictions=['Linear amplitude perturbations only; no core, vortex or phase-slip event',
                      'No nonequilibrium frozen-population susceptibility equivalence',
                      'No mobility calibration or absolute detector latency',
                      'Unknown elastic mean free path limits experimental interpretation of short wavelengths',
                      'No validation of finite delta regularization: q_delta is zero here'],
        runtime_seconds=time.monotonic()-start)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError('Choose a new output path; preserve previous reference')
    result = calculate()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, ensure_ascii=False, allow_nan=False)+'\n', encoding='utf8')
    print(json.dumps(dict(output=str(args.output), runtime_seconds=result['runtime_seconds'],
                         cutoff_comparison=result['cutoff_comparison'], results=result['results']), indent=2))


if __name__ == '__main__':
    main()
