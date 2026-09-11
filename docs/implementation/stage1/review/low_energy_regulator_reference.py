"""Resolve the regulator bias independently of field interpolation/quadrature.

At Gamma=0, the finite-eta BCS inverse count is known analytically. Its Gamma
derivative is -Delta^2*eta*x / ((x^2+eta^2)^2+Delta^2*eta^2). For continuous p,
the eta->0 conjugate tends pi*Delta/2*(1-2*p(0)). This script compares that limit
with adaptive integrals after resolving the x~sqrt(Delta*eta) boundary layer.
No experimental implementation is imported.
"""
from pathlib import Path
import json
import time

import numpy as np
from scipy.integrate import quad


def main():
    started = time.perf_counter()
    amplitude = .72
    population_zero = .2
    population_scale = .25
    causal_limit = np.pi*amplitude/2*(1-2*population_zero)
    rows = []
    for eta in (1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7):
        scale = np.sqrt(amplitude*eta)
        def integrand(coordinate):
            count = scale*coordinate
            kernel = -amplitude**2*eta*count/((count**2+eta**2)**2+amplitude**2*eta**2)
            population = population_zero*np.exp(-count/population_scale)
            return kernel*population*scale
        integral, error = quad(integrand, 0, np.inf, epsabs=1e-12, epsrel=1e-11, limit=150)
        conjugate = np.pi*amplitude/2+4*integral
        rows.append({'eta_over_Delta0': eta, 'Gamma_conjugate_normalized': conjugate,
                     'causal_limit': causal_limit, 'relative_regulator_bias': conjugate/causal_limit-1,
                     'quad_estimated_absolute_error_after_factor4': 4*error})
    result = {
        'scope': 'Independent BCS Gamma=0 reference; no field interpolation, no fixed-count quadrature error',
        'units': 'Delta0=1, N0=1', 'amplitude': amplitude,
        'occupation': 'p(x)=0.2 exp(-x/0.25), held fixed',
        'causal_limit': causal_limit, 'rows': rows,
        'decision': 'eta=1e-3 has an 8.61 percent regulator bias for this occupation; not converged for its small-current response',
        'runtime_seconds': time.perf_counter()-started,
    }
    assert rows[1]['relative_regulator_bias'] > .08
    assert all(rows[i+1]['relative_regulator_bias'] < rows[i]['relative_regulator_bias'] for i in range(len(rows)-1))
    Path(__file__).with_suffix('.json').write_text(json.dumps(result, indent=2)+'\n', encoding='utf-8')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
