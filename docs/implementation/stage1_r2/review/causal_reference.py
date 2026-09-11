"""Independent eta=0 uniform references using a real spectral-angle coordinate.

The retarded angle has sin(Re Theta)=s and cosh(Im Theta)=y>=1. The real
Usadel equation gives y^3+(s^2-1)y-(Delta/Gamma)s=0. This parametrization
contains no quartic complex roots and no finite numerical broadening. It
integrates all occupation profiles over the exact count coordinate with an
analytic Jacobian. Gamma=0 uses the separate exact BCS endpoint formulas.
"""
from pathlib import Path
import hashlib
import json
import math
import sys
import time

import numpy as np
from scipy.integrate import quad, quad_vec
from scipy.optimize import brentq

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from gamma_zero_reference import population, reference as gamma_zero

ROOT = HERE.parents[3]
CRITERIA = ROOT/'docs/implementation/stage1_r2/acceptance_criteria.json'


def vacuum_reference(amplitude, gamma):
    if amplitude == 0:
        return np.zeros(3)
    a, g = amplitude, gamma
    if g <= a:
        return np.array([a*a*(math.log(a)-.5)+math.pi*a*g/2-g*g/3,
                         2*a*math.log(a)+math.pi*g/2, math.pi*a/2-2*g/3])
    def conjugate(kappa):
        ratio = a/kappa
        cosine = math.sqrt((1-ratio)*(1+ratio))
        integral_sine_cubed = ratio**4*(cosine+2)/(3*(1+cosine)**2)
        return a*math.asin(ratio)-kappa*integral_sine_cubed
    def derivative(kappa):
        ratio = a/kappa
        return math.asin(ratio)+ratio*math.sqrt((1-ratio)*(1+ratio))
    value_at_a = a*a*(math.log(a)-.5)+math.pi*a*a/2-a*a/3
    value = value_at_a+quad(conjugate, a, g, epsabs=2e-13, epsrel=2e-13)[0]
    force = 2*a*math.log(a)+math.pi*a/2+quad(derivative, a, g, epsabs=2e-13, epsrel=2e-13)[0]
    return np.array([value, force, conjugate(g)])


def parametric_state(t, amplitude, gamma):
    """Return count x, E, E_Delta, E_Gamma and dx/dt at t>0."""
    a, g = amplitude, gamma
    b = a/g
    cap = min(1., b)
    s = cap*math.exp(-t)
    rhs = s*((b-cap)+cap*(-math.expm1(-t)))
    r = min(rhs/(s*s+2), rhs**(1/3))
    for _ in range(35):
        slope = 3*r*(r+2)+s*s+2
        change = (r*(r*(r+3)+s*s+2)-rhs)/slope
        r -= change
        if abs(change) <= 4e-15*max(r, 1e-100):
            break
    else:
        raise ArithmeticError('real cubic failed to converge')
    y = 1+r
    h = r*(2+r)
    k = (1-cap*cap)+cap*cap*(-math.expm1(-2*t))
    v = math.sqrt(h)
    count = g*y*v*k**1.5/s
    energy = g*v*(h+k)/s
    rp = s*(2*s*y-b)/(3*y*y+s*s-1)
    jacobian = count*(1+rp*(1/y+y/h)+3*s*s/k)
    if not all(math.isfinite(item) for item in (count, energy, jacobian)) or jacobian <= 0:
        raise ArithmeticError('invalid real spectral coordinate')
    return count, energy, v/y, -s*v, jacobian


def t_at_count(count, amplitude, gamma):
    if count == 0:
        return 0.
    # The endpoint is known analytically; do not evaluate its singular Jacobian.
    return brentq(lambda t: -count if t == 0 else parametric_state(t, amplitude, gamma)[0]-count,
                  0., 100., xtol=5e-15, rtol=5e-15)


def integrate_profiles(amplitude, gamma, profiles, *, cutoff=20., tolerance=1e-11):
    if gamma == 0 or amplitude == 0:
        rows = [gamma_zero(amplitude, profile, precision=tolerance/10) for profile in profiles]
        return np.array([row[0] for row in rows]), max(max(row[1]) for row in rows)
    count_breaks = [0., .001, .005, .015, .035, .05, .075, .1, .25, .5, 1., 2., 3.2, 5., 8., 12., cutoff]
    ts = sorted(set(t_at_count(x, amplitude, gamma) for x in count_breaks if x <= cutoff))
    tmax = ts[-1]
    ts = sorted(set(ts+[t for t in (.001, .01, .1, .5, 1., 2., 4., 8., 16.) if t < tmax]))
    def integrand(t):
        x, energy, force, response, jacobian = parametric_state(t, amplitude, gamma)
        weights = np.array([population(profile, x) for profile in profiles])
        return 4*jacobian*weights[:, None]*np.array([energy, force, response])[None, :]
    integral = np.zeros((len(profiles), 3))
    estimate = 0.
    for lower, upper in zip(ts[:-1], ts[1:]):
        values, error = quad_vec(integrand, lower, upper, epsabs=tolerance, epsrel=tolerance, limit=150)
        integral += values
        estimate += error
    return integral+vacuum_reference(amplitude, gamma)[None, :], estimate


def main():
    start = time.perf_counter()
    criteria = json.loads(CRITERIA.read_text(encoding='utf-8'))
    profiles = [profile['id'] for profile in criteria['mandatory_population_profiles']]
    results = []
    for amplitude, gamma in criteria['mandatory_field_points']:
        coarse, estimate = integrate_profiles(amplitude, gamma, profiles, cutoff=12., tolerance=1e-10)
        fine, fine_estimate = integrate_profiles(amplitude, gamma, profiles, cutoff=20., tolerance=2e-12)
        for i, profile in enumerate(profiles):
            delta = np.abs(fine[i]-coarse[i])
            assert np.max(delta) < 1e-8, (amplitude, gamma, profile, delta)
            results.append({'amplitude': amplitude, 'gamma': gamma, 'profile': profile,
                            'reference': fine[i].tolist(), 'reference_error_estimate': fine_estimate,
                            'cutoff_and_quadrature_refinement_absolute_difference': delta.tolist()})
    result = {'schema': 'pysnspd.stage1_r2.independent_causal_reference.v1',
              'criteria_sha256': hashlib.sha256(CRITERIA.read_bytes()).hexdigest(),
              'reference_code_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              'method': 'eta=0 real spectral-angle parametrization; adaptive integration in exact state count; no experimental implementation imported',
              'observable_order': ['energy', 'amplitude_force', 'Gamma_response'],
              'rows': results, 'runtime_seconds': time.perf_counter()-start,
              'cutoff_comparison': [12, 20], 'quadrature_tolerances': [1e-10, 2e-12]}
    Path(__file__).with_suffix('.json').write_text(json.dumps(result, indent=2)+'\n', encoding='utf-8')
    print(json.dumps({'cases': len(results), 'runtime_seconds': result['runtime_seconds'],
                      'max_refinement_difference': max(max(row['cutoff_and_quadrature_refinement_absolute_difference']) for row in results)}, indent=2))


if __name__ == '__main__':
    main()
