"""At most 17 fixed-state RHS samples; no time integration or kernel edits."""
from pathlib import Path
import argparse
import hashlib
import json
import sys
import time

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT/'sandbox/stage2_cells'))
import numpy as np
from run_coupled import setup
from pysnspd.experimental.cell_validation import ElectronicCell


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--initial', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    start = time.perf_counter()
    system, expected, _ = setup('one', phonon_nodes=1025, infrared=.005,
                                 escape=np.inf, heating=0.)
    archived = np.load(args.initial)
    y = archived['state']
    if not np.array_equal(y, expected):
        raise ValueError('The archived initial vector must exactly match setup.')
    for key, expected_grid in [('electron_count', system.catalog.count_nodes),
                               ('electron_weights', system.catalog.count_weights),
                               ('phonon_energies', system.phonons.energies),
                               ('phonon_capacities', system.phonons.capacities)]:
        if not np.array_equal(archived[key], expected_grid):
            raise ValueError('Grid mismatch: '+key)
    ne = system.electron_size
    np_ = system.phonon_size
    es, ps = slice(1, 1+ne), slice(1+ne, 1+ne+np_)
    capacities = np.r_[4*system.catalog.count_weights, system.phonons.capacities]
    populations = slice(1, 1+ne+np_)
    def norm(v):
        return float(np.dot(capacities, abs(v[populations])))
    rhs_times = []
    def evaluate(state):
        before = time.perf_counter()
        result = system.rhs(0., state)
        rhs_times.append(time.perf_counter()-before)
        return result
    f0 = evaluate(y)
    fscale = norm(f0)
    hlist = [1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8]
    central = []
    for h in hlist:
        plus, minus = y.copy(), y.copy()
        plus[0] += h; minus[0] -= h
        fp, fm = evaluate(plus), evaluate(minus)
        slope = (fp-fm)/(2*h)
        second = fp+fm-2*f0
        central.append((h, slope, second, fp, fm))
    amplitude = []
    for i, (h, slope, second, fp, fm) in enumerate(central):
        row = dict(h=h, slope_weighted_L1=norm(slope),
                   second_difference_weighted_L1_relative=norm(second)/fscale,
                   second_derivative_weighted_L1=norm(second)/(h*h),
                   max_second_difference_raw=float(np.max(abs(second[populations]))),
                   argmax_second_difference_population=int(np.argmax(abs(second[populations]))),
                   plus_relative_RHS_change=norm(fp-f0)/fscale)
        if i:
            row['slope_relative_change_from_previous'] = norm(slope-central[i-1][1])/max(norm(slope), norm(central[i-1][1]))
        amplitude.append(row)
    frozen_field_direction = f0.copy()
    frozen_field_direction[0] = 0.
    frozen_field_direction[system.population_size:] = 0.
    direction_samples = []
    for step in [1e-5, 1e-6]:
        evaluated = evaluate(y+step*frozen_field_direction)
        jv = (evaluated-f0)/step
        direction_samples.append(dict(probe_scale=step,
            weighted_Jv_over_v_per_time=norm(jv)/norm(frozen_field_direction),
            max_raw_Jv=float(np.max(abs(jv[populations])))))

    # This exact frozen-event phonon diagonal uses the same finite operator,
    # not a time integration or an assumption about the full spectrum.
    build_start = time.perf_counter()
    cell = ElectronicCell(system.catalog, y[0], 0.)
    events = system.reaction_events(cell)
    lf, lr = events.log_activities(y[es], y[ps])
    forward = events.coefficients*np.exp(lf)
    reverse = events.coefficients*np.exp(lr)
    diagonal = np.zeros(np_)
    for index, beta in [(events.phonon_lower, events.beta_lower),
                        (events.phonon_upper, events.beta_upper)]:
        diagonal += np.bincount(index, weights=beta*beta*(
            forward/(1+y[ps][index])-reverse/y[ps][index]), minlength=np_)
    diagonal /= system.phonons.capacities
    diagonal_seconds = time.perf_counter()-build_start
    selected = int(np.argmin(diagonal))
    perturb = max(1e-12, y[ps][selected]*1e-4)
    plus, minus = y.copy(), y.copy()
    plus[1+ne+selected] += perturb; minus[1+ne+selected] -= perturb
    column = (evaluate(plus)-evaluate(minus))/(2*perturb)
    source_paths = [Path(__file__), ROOT/'sandbox/stage2_cells/run_coupled.py',
                    *sorted((ROOT/'pysnspd/experimental').glob('*.py'))]
    result = dict(schema='pysnspd.stage2.rhs_smoothness.v1',
        scope='Static fixed-state samples only; no trajectory, no eigenvalue claim.',
        source_hashes={p.relative_to(ROOT).as_posix():sha(p) for p in source_paths},
        initial_sha256=sha(args.initial), initial_path=str(args.initial),
        criteria_sha256=sha(ROOT/'docs/implementation/stage2/acceptance_criteria.json'),
        event_count=int(len(events.omega)), rhs_calls=system.rhs_calls,
        rhs_seconds=rhs_times, total_seconds=time.perf_counter()-start,
        base=dict(amplitude=float(y[0]), velocity=float(f0[0]),
            population_weighted_L1_RHS=fscale,
            maximum_raw_electron_derivative=float(np.max(abs(f0[es]))),
            maximum_raw_phonon_derivative=float(np.max(abs(f0[ps]))),
            min_electron_capacity=float(capacities[:ne].min()),
            min_phonon_capacity=float(capacities[ne:].min())),
        amplitude_finite_differences=amplitude,
        frozen_amplitude_flow_direction=direction_samples,
        phonon_diagonal=dict(minimum=float(diagonal.min()), maximum=float(diagonal.max()),
            worst_node=selected, energy=float(system.phonons.energies[selected]),
            capacity=float(system.phonons.capacities[selected]),
            numerical_diagonal=float(column[1+ne+selected]),
            weighted_column_L1_over_capacity=norm(column)/system.phonons.capacities[selected],
            evaluation_seconds=diagonal_seconds,
            limitation='Frozen-field phonon diagonal/one column and one flow direction do not bound all eigenvalues.'))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2)+'\n', encoding='utf-8')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
