"""Independent omitted absorption bounds at actual moving-cell checkpoints."""
from pathlib import Path
import argparse
import hashlib
import json
import sys
import time
ROOT = Path(__file__).resolve().parents[2]
sys.path[:0] = [str(ROOT), str(ROOT/'sandbox/stage2_cells'),
                str(ROOT/'docs/implementation/stage2/review')]
import numpy as np
from run_coupled import setup, jsonable
from pysnspd.experimental.cell_validation import ElectronicCell
from upper_support_reference import upper_absorption_bound


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('trajectories', nargs='+', type=Path)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    started = time.perf_counter()
    rows = []
    for path in args.trajectories:
        data = json.loads(path.read_text())
        if sha(path.with_suffix('.npz')) != data['trajectory_sha256']:
            raise ValueError('trajectory hash mismatch')
        if any(sha(ROOT/key) != value for key,value in data['source_hashes'].items()):
            raise ValueError('trajectory sources differ from current validation code')
        par = data['parameters']
        system, _, debye = setup(par['case'], par['phonon_nodes'], par['infrared'],
            par['face_order'], float(par['escape']), par['heating'], par['scenario'],
            par['reaction_order'], par['reaction_method'], par.get('electron_refinement', 0),
            par.get('reaction_layout','global'), par.get('reaction_outer_order',2),
            par.get('reaction_max_panel',.5), par.get('reaction_max_energy_panel',.125))
        with np.load(path.with_suffix('.npz')) as saved:
            states, times = saved['states'], saved['times']
            for key, expected in (('electron_count',system.catalog.count_nodes),
                ('electron_weights',system.catalog.count_weights),
                ('phonon_energies',system.phonons.energies),
                ('phonon_capacities',system.phonons.capacities)):
                if not np.array_equal(saved[key],expected):
                    raise ValueError('trajectory grid mismatch: '+key)
        extrema = []
        for cell in range(system.cell_count):
            for quantity in ('temperatures','phonon_energy','amplitudes'):
                values = [s[quantity][cell] for s in data['snapshots']]
                extrema += [int(np.argmin(values)),int(np.argmax(values))]
            offset = cell*system.block_size+1
            tail = states[:, offset:offset+system.electron_size][:,system.catalog.count_nodes > 8.]
            tail_count = tail @ system.catalog.count_weights[system.catalog.count_nodes > 8.]
            extrema.append(int(np.argmax(tail_count)))
        indices = np.unique(np.r_[np.linspace(0, len(states)-1, 11, dtype=int),extrema])
        for k in indices:
            amplitudes, populations, phonons = system.unpack(states[k])
            for i, amplitude in enumerate(amplitudes):
                if system.gammas[i] != 0:
                    raise ValueError('this independent bound supports zero-Gamma trajectories')
                cell = ElectronicCell(system.catalog, amplitude, 0.)
                event = system.reaction_events(cell)
                logf, logb = event.log_activities(populations[i], phonons[i])
                gross = event.coefficients*(np.exp(logf)+np.exp(logb))
                gross_number, gross_power = np.sum(gross), np.dot(event.omega, gross)
                bounds = []
                for order in (16, 32):
                    bound = upper_absorption_bound(float(amplitude), 0., cell.energies[-1],
                        lambda energy: np.interp(energy, cell.energies, populations[i]),
                        lambda omega: np.interp(omega, system.phonons.energies, phonons[i]),
                        order=order, omega_min=par['infrared'], omega_max=4.,
                        coupling=debye.alpha2F)
                    bounds.append(bound)
                fractions = {}
                for key, denominator in (
                    ('absorption_number_rate_upper_bound', gross_number),
                    ('absorption_phonon_power_upper_bound', gross_power)):
                    fractions[key] = bounds[-1][key]*system.rate_prefactor/max(denominator, 1e-300)
                rows.append(dict(trajectory=path.as_posix(), time=float(times[k]), cell=i,
                    amplitude=float(amplitude), source_sha256=sha(path),
                    gross_event_number=float(gross_number), gross_event_power=float(gross_power),
                    reference_orders=[16, 32], bounds=bounds, relative_bounds=fractions))
    worst = max(value for row in rows for value in row['relative_bounds'].values())
    result = dict(status='PASS' if worst <= 1e-3 else 'FAIL', rows=rows,
        worst_relative_upper_bound=worst,
        scope='Absorption to unresolved higher electron energies, external holes bounded by1. '
              'Linear p and n bound the geometric positive activities from above. '
              'No assumption of a populated external electron reservoir is introduced.',
        runner_sha256=sha(__file__), runtime_seconds=time.perf_counter()-started)
    result['amendment_sha256'] = sha(ROOT/'docs/implementation/stage2/review/complementary_grid_amendment.json')
    result['complementary_source_sha256'] = sha(ROOT/'pysnspd/experimental/refined_cells.py')
    args.output.write_text(json.dumps(jsonable(result), indent=2, allow_nan=False)+'\n')
    print(json.dumps({k: result[k] for k in ('status','worst_relative_upper_bound','runtime_seconds')}))


if __name__ == '__main__':
    main()
