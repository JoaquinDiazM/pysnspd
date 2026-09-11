"""Independent analytic BCS/count references for the frozen r2 acceptance set.

No experimental catalogue or spectral solver is imported. Adaptive integrals
resolve both the specified occupation bands and the finite-regulator layer.
The Gamma=0 causal current slope is an analytic endpoint result, not a
pointwise eta=0 spectral quadrature.
"""
from pathlib import Path
import hashlib
import json
import math
import time

from scipy.integrate import quad

ROOT = Path(__file__).resolve().parents[4]
CRITERIA = ROOT/'docs/implementation/stage1_r2/acceptance_criteria.json'


def population(name, x):
    if name == 'vacuum':
        return 0.0
    if name == 'low_energy':
        return .2*math.exp(-x/.25)
    if name == 'finite_energy_band':
        return .16*math.exp(-((x-.95)/.3)**2)
    if name == 'near_edge_band':
        return .2*math.exp(-((x-.035)/.01)**2)
    if name == 'high_energy_band':
        return .12*math.exp(-((x-3.2)/.55)**2)
    raise ValueError(name)


def integrate(function, amplitude, eta=0, precision=1e-12):
    points = {0., .001, .005, .015, .035, .05, .075, .1, .25, .5, 1., 2., 3.2, 5., 8.}
    if eta > 0 and amplitude > 0:
        layer = math.sqrt(amplitude*eta)
        points.update((eta, layer, 3*layer, 10*layer))
    bounds = sorted(points)+[math.inf]
    results = [quad(function, low, high, epsabs=precision, epsrel=precision, limit=120)
               for low, high in zip(bounds[:-1], bounds[1:])]
    return sum(value for value, _ in results), sum(error for _, error in results)


def reference(amplitude, profile, eta=0, precision=1e-12):
    p = lambda x: population(profile, x)
    vacuum = amplitude**2*(math.log(amplitude)-.5) if amplitude else 0.
    vacuum_force = 2*amplitude*math.log(amplitude) if amplitude else 0.
    if amplitude == 0:
        excitation, error = integrate(lambda x: x*p(x), amplitude, precision=precision)
        return [4*excitation, 0., 0.], [4*error, 0., 0.]
    if eta == 0:
        excitation, ee = integrate(lambda x: math.hypot(x, amplitude)*p(x), amplitude, precision=precision)
        force, fe = integrate(lambda x: amplitude/math.hypot(x, amplitude)*p(x), amplitude, precision=precision)
        gamma = math.pi*amplitude/2*(1-2*p(0))
        return [vacuum+4*excitation, vacuum_force+4*force, gamma], [4*ee, 4*fe, 0.]
    def energy(x):
        return x*math.sqrt(1+amplitude**2/(x*x+eta*eta))
    def amplitude_derivative(x):
        return x*amplitude/(math.sqrt(x*x+eta*eta)*math.sqrt(x*x+eta*eta+amplitude**2))
    def gamma_derivative(x):
        return -amplitude**2*eta*x/((x*x+eta*eta)**2+amplitude**2*eta*eta)
    excitation, ee = integrate(lambda x: energy(x)*p(x), amplitude, eta, precision)
    force, fe = integrate(lambda x: amplitude_derivative(x)*p(x), amplitude, eta, precision)
    gamma, ge = integrate(lambda x: gamma_derivative(x)*p(x), amplitude, eta, precision)
    return [vacuum+4*excitation, vacuum_force+4*force, math.pi*amplitude/2+4*gamma], [4*ee, 4*fe, 4*ge]


def main():
    start = time.perf_counter()
    criteria = json.loads(CRITERIA.read_text(encoding='utf-8'))
    profiles = [profile['id'] for profile in criteria['mandatory_population_profiles']]
    amplitudes = sorted({point[0] for point in criteria['mandatory_field_points'] if point[1] == 0} | {0.})
    rows = []
    for amplitude in amplitudes:
        for profile in profiles:
            ideal, ideal_error = reference(amplitude, profile)
            refined, _ = reference(amplitude, profile, precision=2e-13)
            refinement = [abs(a-b) for a, b in zip(ideal, refined)]
            assert max(ideal_error+refinement) < 1e-9
            row = {'amplitude': amplitude, 'gamma': 0, 'profile': profile,
                   'causal_reference': ideal, 'quadrature_error_estimate': ideal_error,
                   'reference_refinement_absolute_difference': refinement,
                   'finite_eta_references': []}
            for eta in (3e-6, 3e-7, 3e-8):
                values, errors = reference(amplitude, profile, eta)
                row['finite_eta_references'].append({
                    'eta': eta, 'values': values, 'quadrature_error_estimate': errors,
                    'scaled_bias_to_causal': [abs(v-r)/max(1, abs(r)) for v, r in zip(values, ideal)],
                    'relative_bias_to_causal': [abs(v/r-1) if r else None for v, r in zip(values, ideal)]})
            rows.append(row)
    result = {'schema': 'pysnspd.stage1_r2.independent_gamma_zero.v1',
              'criteria_sha256': hashlib.sha256(CRITERIA.read_bytes()).hexdigest(),
              'reference_code_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              'units': ['N0 Delta0^2', 'N0 Delta0', 'N0 Delta0'],
              'observable_order': ['energy', 'amplitude_force', 'Gamma_response'],
              'purpose': 'Independent references only; does not evaluate or admit the r2 catalogue.',
              'rows': rows, 'runtime_seconds': time.perf_counter()-start}
    Path(__file__).with_suffix('.json').write_text(json.dumps(result, indent=2)+'\n', encoding='utf-8')
    print(json.dumps({'reference_rows': len(rows), 'runtime_seconds': result['runtime_seconds'],
                      'criteria_sha256': result['criteria_sha256']}, indent=2))


if __name__ == '__main__':
    main()
