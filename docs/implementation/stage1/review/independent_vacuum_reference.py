"""Independent zero-temperature benchmarks for the stage-1 electronic catalogue.

This review calculation does not import the experimental implementation. Energy
units are Delta0 and density units N0*Delta0**2. Its reference is a direct
imaginary-frequency integral, evaluated after subtracting the zero-current BCS
integrand; the closed gapped formula supplies a second, analytic result.
"""

from pathlib import Path
import json
import time

import numpy as np
from scipy.optimize import brentq
from numpy.polynomial.legendre import leggauss


def angle(omega, amplitude, gamma):
    if amplitude == 0:
        return 0.0
    upper = np.pi / 2 if gamma <= amplitude else np.arcsin(amplitude / gamma)
    return brentq(
        lambda theta: omega * np.sin(theta) - amplitude * np.cos(theta)
        + gamma * np.sin(theta) * np.cos(theta),
        0.0, upper, xtol=5e-15, rtol=1e-14,
    )


def vacuum_by_frequency_integral(amplitude, gamma, order):
    if amplitude == 0:
        return np.zeros(3)
    nodes, weights = leggauss(order)
    t = (nodes + 1) / 2
    omega = t / (1 - t)
    weights = weights / (2 * (1 - t)**2)
    theta = np.array([angle(w, amplitude, gamma) for w in omega])
    sin = np.sin(theta)
    cos = np.cos(theta)
    theta0 = np.arctan(amplitude / omega)
    sin0 = np.sin(theta0)
    cos0 = np.cos(theta0)
    integrand = 2 * omega * sin**2 / (1 + cos) - 2 * amplitude * sin + gamma * sin**2
    integrand0 = 2 * omega * sin0**2 / (1 + cos0) - 2 * amplitude * sin0
    energy = amplitude**2 * (np.log(amplitude) - 0.5) + weights @ (integrand - integrand0)
    force = 2 * amplitude * np.log(amplitude) + 2 * weights @ (sin0 - sin)
    gamma_derivative = weights @ sin**2
    return np.array([energy, force, gamma_derivative])


def analytic_gapped(amplitude, gamma):
    assert 0 <= gamma <= amplitude and amplitude > 0
    return np.array([
        amplitude**2 * (np.log(amplitude) - 0.5)
        + np.pi * amplitude * gamma / 2 - gamma**2 / 3,
        2 * amplitude * np.log(amplitude) + np.pi * gamma / 2,
        np.pi * amplitude / 2 - 2 * gamma / 3,
    ])


def main():
    start = time.perf_counter()
    rows = []
    for amplitude in (0.15, 0.6, 1.0, 1.4):
        for fraction in (0.0, 0.2, 0.8, 1.0, 1.5, 4.0):
            gamma = amplitude * fraction
            coarse = vacuum_by_frequency_integral(amplitude, gamma, 128)
            fine = vacuum_by_frequency_integral(amplitude, gamma, 256)
            row = {
                'amplitude_over_Delta0': amplitude,
                'Gamma_over_Delta0': gamma,
                'energy_over_N0_Delta0_squared': float(fine[0]),
                'amplitude_force_over_N0_Delta0': float(fine[1]),
                'Gamma_derivative_over_N0_Delta0': float(fine[2]),
                'refinement_error_scaled': float(np.max(abs(fine-coarse) / np.maximum(1, abs(fine)))),
            }
            if fraction <= 1:
                exact = analytic_gapped(amplitude, gamma)
                row['analytic_error_scaled'] = float(np.max(abs(fine-exact) / np.maximum(1, abs(exact))))
            rows.append(row)
    max_refinement = max(row['refinement_error_scaled'] for row in rows)
    max_analytic = max(row.get('analytic_error_scaled', 0) for row in rows)
    # Gamma=amplitude is a nonanalytic threshold; convergence is slower there.
    assert max_refinement < 2e-5, max_refinement
    assert max_analytic < 2e-5, max_analytic
    result = {
        'scope': 'Independent imaginary-frequency integration, no production solver or experimental-module import',
        'normalization': 'Delta0=1, N0=1; Gamma=hbar*D*q^2/2',
        'state_count': len(rows), 'orders': [128, 256],
        'max_refinement_error_scaled': max_refinement,
        'max_gapped_analytic_error_scaled': max_analytic,
        'runtime_seconds': time.perf_counter()-start,
        'rows': rows,
    }
    path = Path(__file__).with_suffix('.json')
    path.write_text(json.dumps(result, indent=2)+'\n', encoding='utf-8')
    print(json.dumps({k: v for k, v in result.items() if k != 'rows'}, indent=2))


if __name__ == '__main__':
    main()
