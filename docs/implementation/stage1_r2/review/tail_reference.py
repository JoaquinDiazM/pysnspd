"""Adaptive causal tail integrals to infinity for the actual acceptance profiles."""
from pathlib import Path
import hashlib
import json
import math
import time

import numpy as np
from scipy.integrate import quad_vec

from causal_reference import CRITERIA, parametric_state, t_at_count
from gamma_zero_reference import population


def main():
    started = time.perf_counter()
    criteria = json.loads(CRITERIA.read_text(encoding='utf-8'))
    fields = {tuple(field):[profile['id'] for profile in criteria['mandatory_population_profiles']]
              for field in criteria['mandatory_field_points']}
    for field in criteria['thermal_equilibrium_checks']['field_points']:
        fields.setdefault(tuple(field), []).extend('FD_'+str(t) for t in criteria['thermal_equilibrium_checks']['kBT_over_Delta0'])
    rows = []
    for (a, g), profiles in fields.items():
        for cutoff in (12., 20.):
            lower = cutoff if g == 0 else t_at_count(cutoff, a, g)
            def integrand(coordinate):
                if g == 0:
                    x = coordinate
                    energy = math.hypot(x, a)
                    ea, eg, jacobian = a/energy, 0., 1.
                else:
                    x, energy, ea, eg, jacobian = parametric_state(coordinate, a, g)
                values = []
                for name in profiles:
                    if name.startswith('FD_'):
                        temperature = float(name[3:])
                        exponential = np.exp(-energy/temperature)
                        p = exponential/(1+exponential)
                        values.append([-4*temperature*jacobian*np.logaddexp(0., -energy/temperature),
                                       4*jacobian*p*ea, 4*jacobian*p*eg])
                    else:
                        p = population(name, x)
                        values.append(4*jacobian*p*np.array([energy, ea, eg]))
                return np.asarray(values)
            tail, estimate = quad_vec(integrand, lower, np.inf, epsabs=1e-14, epsrel=1e-10, limit=100)
            for index, profile in enumerate(profiles):
                rows.append(dict(amplitude=a, gamma=g, profile=profile, cutoff=cutoff,
                                 causal_tail_to_infinity=tail[index].tolist(), error_estimate=estimate))
    result = dict(schema='pysnspd.stage1_r2.independent_infinite_tail.v1',
                  criteria_sha256=hashlib.sha256(CRITERIA.read_bytes()).hexdigest(),
                  code_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                  method='Adaptive integration to infinity using the independent real causal-angle parametrization; Gamma=0 uses the analytic BCS count energy. No catalogue implementation imported.',
                  observable_order=['energy_or_FD_free_energy', 'amplitude_force', 'Gamma_response'],
                  scope='Ideal causal tail, used together with the separately measured finite-eta cutoff differences and regulator convergence.',
                  rows=rows, runtime_seconds=time.perf_counter()-started)
    Path(__file__).with_suffix('.json').write_text(json.dumps(result, indent=2)+'\n', encoding='utf-8')
    print(json.dumps(dict(rows=len(rows), max_absolute_tail=np.max(np.abs([r['causal_tail_to_infinity'] for r in rows])),
                         runtime_seconds=result['runtime_seconds']), indent=2))


if __name__ == '__main__':
    main()
