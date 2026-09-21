"""Independent BCS check of the fields actually visited by zero-current cells."""
from pathlib import Path
import argparse
import hashlib
import json
import sys
import time
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
import numpy as np
from pysnspd.experimental.energy_catalog import OccupationEnergyCatalog
from pysnspd.experimental.refined_cells import ComplementaryCountCatalog


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('trajectories', nargs='+', type=Path)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    started = time.perf_counter()
    path = ROOT/'docs/implementation/stage1_r2/catalogs/occupation_catalog.npz'
    catalog = OccupationEnergyCatalog.load(path)
    base = catalog
    eta = base.eta
    rows = []
    for item in args.trajectories:
        data = json.loads(item.read_text())
        if sha(item.with_suffix('.npz')) != data['trajectory_sha256']:
            raise ValueError('trajectory archive hash mismatch')
        if any(sha(ROOT/key) != value for key,value in data['source_hashes'].items()):
            raise ValueError('current code differs from trajectory source hashes')
        if data['parameters']['scenario'] == 'equilibrium':
            raise ValueError('finite-Gamma equilibrium requires its own causal reference')
        with np.load(item.with_suffix('.npz')) as stored:
            states, times = stored['states'], stored['times']
            phonon_size = len(stored['phonon_energies'])
            x, weights = stored['electron_count'], stored['electron_weights']
        catalog = (ComplementaryCountCatalog(base, x, weights)
                   if data['parameters'].get('electron_refinement', 0) else base)
        cells = len(data['initial']['amplitudes'])
        block = 1+len(x)+phonon_size
        samples = np.unique(np.r_[np.linspace(0, len(states)-1, 11, dtype=int),
                                  [np.argmin(states[:, i*block]) for i in range(cells)],
                                  [np.argmax(states[:, i*block]) for i in range(cells)]])
        for k in samples:
            for i in range(cells):
                amplitude = states[k, i*block]
                p = states[k, i*block+1:i*block+1+len(x)]
                # Solve Re sqrt((E+i eta)^2-a^2)=x algebraically; no spline,
                # stored energy-kernel coefficients or production force used.
                e = x*np.sqrt((x*x+eta*eta+amplitude*amplitude)/(x*x+eta*eta))
                da = x*amplitude/np.sqrt((x*x+eta*eta)*(x*x+eta*eta+amplitude*amplitude))
                reference_energy = amplitude**2*(np.log(amplitude)-.5)+4*np.dot(weights*e, p)
                reference_force = 2*amplitude*np.log(amplitude)+4*np.dot(weights*da, p)
                value, force, _ = catalog.evaluate(amplitude, 0., p)
                h = 1e-5
                finite_difference = (catalog.evaluate(amplitude+h, 0., p)[0]
                                     -catalog.evaluate(amplitude-h, 0., p)[0])/(2*h)
                ideal_e = np.sqrt(x*x+amplitude*amplitude)
                regulator_shift = 4*np.dot(weights*(e-ideal_e), p)
                rows.append(dict(trajectory=item.as_posix(), trajectory_sha256=sha(item),
                                 time=float(times[k]), cell=i, amplitude=float(amplitude),
                                 energy_relative_error=float(abs(value-reference_energy)/max(abs(reference_energy),1e-12)),
                                 force_relative_error=float(abs(force-reference_force)/max(abs(reference_force),1e-12)),
                                 derivative_identity_scaled=float(abs(force-finite_difference)/max(1., abs(force))),
                                 regulator_energy_shift=float(regulator_shift)))
    maxima = {key: max(r[key] for r in rows) for key in
              ('energy_relative_error','force_relative_error','derivative_identity_scaled')}
    passed = (maxima['energy_relative_error'] <= 1e-3 and maxima['force_relative_error'] <= 1e-3
              and maxima['derivative_identity_scaled'] <= 1e-7)
    result = dict(status='PASS' if passed else 'FAIL', scope='Actual zero-current trajectory points only; independent finite-eta BCS algebra on the same count quadrature.',
                  catalog_sha256=sha(path), runner_sha256=sha(__file__),
                  criteria_sha256=sha(ROOT/'docs/implementation/stage2/acceptance_criteria.json'),
                  amendment_sha256=sha(ROOT/'docs/implementation/stage2/review/complementary_grid_amendment.json'),
                  complementary_source_sha256=sha(ROOT/'pysnspd/experimental/refined_cells.py'),
                  samples=len(rows), maxima=maxima, cases=rows,
                  runtime_seconds=time.perf_counter()-started)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({k: result[k] for k in ('status','samples','maxima','runtime_seconds')},indent=2))
    if not passed: raise SystemExit(1)


if __name__ == '__main__':
    main()
