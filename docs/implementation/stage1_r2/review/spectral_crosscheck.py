"""Independent API audit of batched causal roots and state-count inversion.

The unsquared Usadel equation and normalization are recomputed directly. The
scalar API now delegates to the batch implementation: their agreement is an
API check, never independent physical evidence. The real causal-angle and
Matsubara reference scripts provide that independence. Every failure is kept.
"""
from pathlib import Path
import hashlib
import json
import sys
import time

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
sys.path.insert(0, str(ROOT))
from pysnspd.experimental.energy_catalog import retarded_spectrum, retarded_spectrum_batch, energy_at_count_batch


def main():
    started = time.perf_counter()
    module_path = ROOT/'pysnspd/experimental/energy_catalog.py'
    source_at_start = hashlib.sha256(module_path.read_bytes()).hexdigest()
    criteria_path = ROOT/'docs/implementation/stage1_r2/acceptance_criteria.json'
    criteria = json.loads(criteria_path.read_text(encoding='utf-8'))
    rows, failures, scalar_comparator_failures = [], [], []
    for eta in (1e-7, 1e-8, 3e-9):
        for a, g in criteria['mandatory_field_points']:
            identifier = dict(amplitude=a, gamma=g, eta=eta)
            gap = a*max(0., 1-(g/a)**(2/3))**1.5
            energy = np.unique(np.r_[0, np.geomspace(1e-12, 20, 45),
                                     np.maximum(0, a+eta*np.array([-100, -10, -1, 0, 1, 10, 100])),
                                     gap*np.r_[1., 1-np.geomspace(1e-8, .1, 8), 1+np.geomspace(1e-8, .1, 8)]])
            try:
                c, s = retarded_spectrum_batch(energy, delta=a, gamma=g, eta=eta)
                z = energy+1j*eta
                lhs, rhs = a*c, (g*c-1j*z)*s
                equation = abs(lhs-rhs)/np.maximum(1, np.maximum(abs(lhs), abs(rhs)))
                normalization = abs(c*c+s*s-1)/np.maximum(1, abs(c)**2+abs(s)**2)
                scalar_difference = None
                try:
                    cs, ss = retarded_spectrum(energy, delta=a, gamma=g, eta=eta)
                    scalar_difference = np.maximum(abs(c-cs), abs(s-ss))/np.maximum(1, np.maximum(abs(cs), abs(ss)))
                except Exception as exc:
                    rejected_energies = []
                    for e in energy:
                        try:
                            retarded_spectrum(e, delta=a, gamma=g, eta=eta)
                        except Exception as scalar_exc:
                            rejected_energies.append(dict(energy=float(e), exception=f'{type(scalar_exc).__name__}: {scalar_exc}'))
                    scalar_comparator_failures.append(dict(identifier, rejected_energies=rejected_energies,
                                                           exception=f'{type(exc).__name__}: {exc}'))
                counts = np.r_[0, np.geomspace(1e-11, 12, 35)]
                inverted = energy_at_count_batch(counts, delta=a, gamma=g, eta=eta)
                ci, si = retarded_spectrum_batch(inverted, delta=a, gamma=g, eta=eta)
                w = 1j*a/si
                reconstructed = np.real(w-1j*g*a*a/(2*w*w))
                inverse_error = abs(reconstructed-counts)/np.maximum(1, counts)
                row = dict(identifier, energy_points=len(energy), count_points=len(counts),
                           minimum_DOS=float(min(c.real)), max_unsquared_scaled_residual=float(max(equation)),
                           max_normalization_scaled_residual=float(max(normalization)),
                           max_scalar_batch_scaled_difference=float(max(scalar_difference)) if scalar_difference is not None else None,
                           max_count_inverse_scaled_error=float(max(inverse_error)))
                rows.append(row)
                if (row['minimum_DOS'] < 0 or row['max_unsquared_scaled_residual'] > 1e-8 or
                        row['max_normalization_scaled_residual'] > 1e-8 or row['max_count_inverse_scaled_error'] > 1e-8):
                    failures.append(dict(row, reason='mandatory causal or count gate'))
            except Exception as exc:
                failures.append(dict(identifier, exception=f'{type(exc).__name__}: {exc}'))
    source_at_end = hashlib.sha256(module_path.read_bytes()).hexdigest()
    result = dict(schema='pysnspd.stage1_r2.spectral_crosscheck.v1',
                  status=('INCONCLUSIVE_SOURCE_CHANGED' if source_at_start != source_at_end else 'PASS' if not failures else 'FAIL'),
                  criteria_sha256=hashlib.sha256(criteria_path.read_bytes()).hexdigest(),
                  implementation_sha256=source_at_start, implementation_sha256_at_end=source_at_end,
                  reviewer_code_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                  eta_values=[1e-7, 1e-8, 3e-9], rows=rows, failures=failures,
                  scalar_comparator_failures=scalar_comparator_failures,
                  comparison_is_independent=False,
                  comparator_scope='The scalar public API delegates to the batch implementation. Equality checks API behavior only; independent physics comes from the real causal-angle and Matsubara reference scripts. Historical scalar rejections are retained in pre_scalar_unification_spectral_crosscheck.json.',
                  runtime_seconds=time.perf_counter()-started)
    Path(__file__).with_suffix('.json').write_text(json.dumps(result, indent=2)+'\n', encoding='utf-8')
    print(json.dumps(dict(status=result['status'], parameter_cases=len(rows), failures=failures,
                         scalar_comparator_failures=scalar_comparator_failures,
                         maxima={key:max(row[key] for row in rows if row[key] is not None) for key in
                                 ('max_unsquared_scaled_residual', 'max_normalization_scaled_residual',
                                  'max_scalar_batch_scaled_difference', 'max_count_inverse_scaled_error')},
                         runtime_seconds=result['runtime_seconds']), indent=2))


if __name__ == '__main__':
    main()
