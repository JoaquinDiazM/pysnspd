"""Independent causal FD thermodynamics and ultraviolet-converged Matsubara sums.

No catalogue implementation is imported. The real-angle quadrature and real
Matsubara stationary angle are two different spectral representations. The
forces are partial derivatives at fixed occupations; only the value of the
energy is converted to free energy by subtracting occupation entropy.
"""
from pathlib import Path
import hashlib
import json
import math
import time

import numpy as np
from scipy.integrate import quad_vec
from scipy.optimize import brentq

from causal_reference import CRITERIA, parametric_state, t_at_count, vacuum_reference


def thermal_causal(a, g, temperature, *, cutoff=28., tolerance=2e-12):
    def values(energy, ea, eg, jacobian):
        occupation = 1/(math.exp(energy/temperature)+1)
        return 4*jacobian*np.array([-temperature*np.logaddexp(0., -energy/temperature),
                                   ea*occupation, eg*occupation])
    if g == 0:
        def integrand(x):
            energy = math.hypot(x, a)
            return values(energy, a/energy, 0., 1.)
        result, error = quad_vec(integrand, 0., cutoff, epsabs=tolerance, epsrel=tolerance)
        result += vacuum_reference(a, g)
        result[2] = math.pi*a/2*math.tanh(a/(2*temperature))
        return result, error
    breaks = [0., .001, .005, .015, .035, .075, .1, .25, .5, 1., 2., 4., 8., 12., cutoff]
    ts = sorted(set(t_at_count(x, a, g) for x in breaks))
    def integrand(t):
        x, energy, ea, eg, jacobian = parametric_state(t, a, g)
        return values(energy, ea, eg, jacobian)
    result = vacuum_reference(a, g)
    error = 0.
    for lower, upper in zip(ts[:-1], ts[1:]):
        value, estimate = quad_vec(integrand, lower, upper, epsabs=tolerance, epsrel=tolerance)
        result += value
        error += estimate
    return result, error


def thermal_matsubara(a, g, temperature, counts=(256, 512, 1024, 2048, 4096)):
    frequency = math.pi*temperature*(2*np.arange(max(counts))+1)
    if g == 0:
        sine = a/np.hypot(frequency, a)
        cosine = frequency/np.hypot(frequency, a)
    else:
        theta = np.array([brentq(lambda angle: w*math.sin(angle)-a*math.cos(angle)
                                +g*math.sin(angle)*math.cos(angle), 0., math.pi/2,
                                xtol=5e-15, rtol=5e-15) for w in frequency])
        sine, cosine = np.sin(theta), np.cos(theta)
    bcs_ratio = math.pi/math.exp(np.euler_gamma)
    base = np.array([-math.pi**2*temperature**2/3+a*a*math.log(temperature*bcs_ratio),
                     2*a*math.log(temperature*bcs_ratio), 0.])
    summands = np.array([a*a/frequency+2*frequency*sine*sine/(1+cosine)-2*a*sine+g*sine*sine,
                        2*(a/frequency-sine), sine*sine]).T
    cumulative = np.cumsum(summands, axis=0)
    finite = np.array([base+2*math.pi*temperature*cumulative[n-1] for n in counts])
    # Smooth ultraviolet tail has an inverse-count expansion. Two independent
    # four-cutoff extrapolations test the omitted next coefficient.
    inverse = 1/np.asarray(counts, dtype=float)
    first = np.polynomial.polynomial.polyfit(inverse[:-1], finite[:-1], 3)[0]
    final = np.polynomial.polynomial.polyfit(inverse[1:], finite[1:], 3)[0]
    return final, np.abs(final-first), finite


def gamma_zero_finite_eta(a, temperature, eta):
    """Analytic finite-eta inverse and kernels, integrated at fixed count."""
    breaks = np.unique(np.r_[0., eta*np.array([.01, .1, 1, 10, 100]),
                             np.sqrt(a*eta)*np.array([.01, .1, 1, 10, 100]), .1, 1, 10, 28])
    def integrand(x):
        denominator = x*x+eta*eta
        energy = x*np.sqrt(1+a*a/denominator)
        p = np.exp(-energy/temperature)/(1+np.exp(-energy/temperature))
        ea = x*a/np.sqrt(denominator*(denominator+a*a))
        eg = -a*a*eta*x/(denominator**2+a*a*eta*eta)
        return 4*np.array([-temperature*np.logaddexp(0, -energy/temperature), ea*p, eg*p])
    result, error = vacuum_reference(a, 0), 0.
    for lower, upper in zip(breaks[:-1], breaks[1:]):
        value, estimate = quad_vec(integrand, lower, upper, epsabs=1e-14, epsrel=1e-12)
        result += value
        error += estimate
    return result, error


def main():
    started = time.perf_counter()
    criteria = json.loads(CRITERIA.read_text(encoding='utf-8'))
    rows, eta_controls = [], []
    for a, g in criteria['thermal_equilibrium_checks']['field_points']:
        for temperature in criteria['thermal_equilibrium_checks']['kBT_over_Delta0']:
            causal, estimate = thermal_causal(a, g, temperature)
            coarse, coarse_estimate = thermal_causal(a, g, temperature, cutoff=20., tolerance=1e-10)
            matsubara, tail_error, finite = thermal_matsubara(a, g, temperature)
            delta = np.abs(causal-matsubara)
            rows.append({'amplitude': a, 'gamma': g, 'temperature': temperature,
                         'causal_reference': causal.tolist(), 'causal_quadrature_error_estimate': estimate,
                         'causal_cutoff_and_quadrature_change': np.abs(causal-coarse).tolist(),
                         'matsubara_reference': matsubara.tolist(),
                         'matsubara_tail_extrapolation_change': tail_error.tolist(),
                         'matsubara_finite_cutoffs': finite.tolist(),
                         'independent_representation_absolute_difference': delta.tolist()})
            if g == 0:
                for eta in (1e-7, 1e-8, 3e-9):
                    finite_eta, finite_eta_error = gamma_zero_finite_eta(a, temperature, eta)
                    eta_controls.append(dict(amplitude=a, gamma=g, temperature=temperature, eta=eta,
                                             finite_eta_result=finite_eta.tolist(), ideal_reference=causal.tolist(),
                                             quadrature_error_estimate=finite_eta_error,
                                             relative_regulator_bias=(abs(finite_eta-causal)/abs(causal)).tolist()))
    result = {'schema': 'pysnspd.stage1_r2.independent_thermal_reference.v1',
              'criteria_sha256': hashlib.sha256(CRITERIA.read_bytes()).hexdigest(),
              'reference_code_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              'causal_code_sha256': hashlib.sha256(Path(__file__).with_name('causal_reference.py').read_bytes()).hexdigest(),
              'observable_order': ['free_energy', 'amplitude_force', 'Gamma_response'],
              'causal_count_cutoffs': [20., 28.], 'matsubara_counts': [256, 512, 1024, 2048, 4096],
              'method': 'Independent real causal spectral-angle quadrature versus real Matsubara angle sums; two cubic inverse-count ultraviolet extrapolations; exact Gamma=0 boundary response.',
              'rows': rows, 'Gamma_zero_FD_regulator_controls': eta_controls,
              'runtime_seconds': time.perf_counter()-started}
    Path(__file__).with_suffix('.json').write_text(json.dumps(result, indent=2)+'\n', encoding='utf-8')
    print(json.dumps({'cases': len(rows), 'runtime_seconds': result['runtime_seconds'],
                      'max_representation_difference': np.max([row['independent_representation_absolute_difference'] for row in rows]),
                      'max_matsubara_tail_change': np.max([row['matsubara_tail_extrapolation_change'] for row in rows])}, indent=2))


if __name__ == '__main__':
    main()
