"""Small independent reviewer checks for stage-1 catalogue queries and persistence.

These checks exercise public APIs and independent limiting cases. They do not
run the detector, admit the material data, or establish convergence of the
candidate dynamical model. Run from the repository root with NumPy/SciPy.
"""
from pathlib import Path
import hashlib
import json
import sys
import tempfile
import time

import numpy as np

ROOT = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(ROOT))
from pysnspd.experimental.energy_catalog import (
    E_CHARGE_C, HBAR_J_S, OccupationEnergyCatalog, build_occupation_catalog,
    build_vacuum_catalog, retarded_spectrum, vacuum_state,
)


def main():
    start = time.perf_counter()
    reference = json.loads(Path(__file__).with_name('independent_vacuum_reference.json').read_text())
    errors = []
    for row in reference['rows']:
        actual = np.array(vacuum_state(row['amplitude_over_Delta0'], row['Gamma_over_Delta0']))
        expected = np.array([row['energy_over_N0_Delta0_squared'],
                             row['amplitude_force_over_N0_Delta0'],
                             row['Gamma_derivative_over_N0_Delta0']])
        errors.append(float(np.max(abs(actual-expected)/np.maximum(1, abs(expected)))))
    assert max(errors) < 2e-7
    spectral_states = 0
    min_dos = 1e99
    max_normalization = 0
    for amplitude in (0.0, 1e-8, .01, .2, 1., 2.):
        for gamma in (0., 1e-5, .01, .3, 1., 3., 100.):
            for eta in (1e-3, 1e-6, 1e-9):
                energy = np.unique(np.r_[0, np.geomspace(1e-7, 100, 80),
                                         amplitude*(1+np.linspace(-.1, .1, 41))])
                c, s = retarded_spectrum(energy, delta=amplitude, gamma=gamma, eta=eta)
                min_dos = min(min_dos, float(min(c.real)))
                max_normalization = max(max_normalization, float(np.max(
                    abs(c*c+s*s-1)/np.maximum(1, abs(c)**2+abs(s)**2))))
                spectral_states += 1
    assert min_dos >= 0 and max_normalization < 1e-10

    vacuum = build_vacuum_catalog(np.linspace(.15, 1.4, 9), np.linspace(0, 1.2, 9),
                                  Tc_K=8.65, N0_per_J_m3=1e47, D_m2_s=1.581e-4,
                                  reference_hashes={'independent_fixture': 'query-contract-v1'})
    table = build_occupation_catalog(vacuum, count_order=24, count_max=5, eta=.003)
    population = .12*np.exp(-((table.count_nodes-1.1)/.35)**2)
    directory = Path(tempfile.mkdtemp(prefix='review-query-contract-', dir=ROOT/'tmp'))
    path = table.save(directory/'catalog_without_extension')
    assert path.is_file()
    loaded = OccupationEnergyCatalog.load(path)
    for field in ('count_nodes', 'count_weights', 'excitation_energies', 'gamma_energy_derivatives'):
        assert np.array_equal(getattr(table, field), getattr(loaded, field))
    for field in ('delta0_J', 'N0_per_J_m3', 'D_m2_s'):
        assert getattr(table.vacuum, field) == getattr(loaded.vacuum, field)
    assert loaded.eta == table.eta
    assert loaded.vacuum.metadata['reference_hashes'] == table.vacuum.metadata['reference_hashes']
    rejected = []
    cases = [('NaN occupation', .6, .2, np.full(24, np.nan)),
             ('infinite occupation', .6, .2, np.full(24, np.inf)),
             ('negative occupation', .6, .2, np.full(24, -1e-12)),
             ('occupation above one', .6, .2, np.full(24, 1+1e-12)),
             ('low positive amplitude', .01, .2, population),
             ('Gamma above support', .6, 1.3, population),
             ('negative Gamma', .6, -.1, population)]
    for label, amplitude, gamma, occupation in cases:
        try:
            loaded.evaluate(amplitude, gamma, occupation)
        except ValueError:
            rejected.append(label)
    assert len(rejected) == len(cases)
    normal = loaded.evaluate(0, .2, population)
    normal_reference = 4*np.dot(table.count_nodes*population, table.count_weights)
    assert normal == (normal_reference, 0, 0)
    def centered(fun, coordinate, step=1e-5):
        return (fun(coordinate-2*step)-8*fun(coordinate-step)
                +8*fun(coordinate+step)-fun(coordinate+2*step))/(12*step)
    derivative_errors = []
    maxwell_errors = []
    for amplitude, gamma in ((.23, .04), (.67, .43), (1.2, .85)):
        value, force, conjugate = loaded.evaluate(amplitude, gamma, population)
        derivative_errors.extend([
            abs(force-centered(lambda d: loaded.evaluate(d, gamma, population)[0], amplitude)),
            abs(conjugate-centered(lambda g: loaded.evaluate(amplitude, g, population)[0], gamma)),
        ])
        maxwell_errors.append(abs(
            centered(lambda d: loaded.evaluate(d, gamma, population)[2], amplitude)
            - centered(lambda g: loaded.evaluate(amplitude, g, population)[1], gamma)))
        assert loaded.evaluate(amplitude, gamma, population) == table.evaluate(amplitude, gamma, population)
    endpoint_errors = []
    step = 1e-6
    for amplitude in (.15, .67, 1.4):
        values = [loaded.evaluate(amplitude, offset*step, population)[0] for offset in range(4)]
        derivative = (-11*values[0]+18*values[1]-9*values[2]+2*values[3])/(6*step)
        endpoint_errors.append(abs(derivative-loaded.evaluate(amplitude, 0, population)[2]))
    assert max(derivative_errors+maxwell_errors+endpoint_errors) < 1e-6
    minimum = 1e99
    spacing = 1e99
    points = np.random.default_rng(89341).uniform([.15, 0], [1.4, 1.2], size=(200, 2))
    for amplitude, gamma in points:
        energies = loaded.energy_kernel(amplitude, gamma)[0]
        minimum = min(minimum, float(min(energies)))
        spacing = min(spacing, float(min(np.diff(energies))))
    assert minimum > 0 and spacing > 0
    delta = .63*vacuum.delta0_J
    q = np.sqrt(2*.73*vacuum.delta0_J/(HBAR_J_S*vacuum.D_m2_s))
    hq, hd = q*1e-5, delta*1e-5
    state = vacuum.evaluate_si(delta, q)
    derivative_q = (vacuum.evaluate_si(delta, q+hq)['Uvac_J_m3']
                    - vacuum.evaluate_si(delta, q-hq)['Uvac_J_m3'])/(2*hq)
    derivative_delta = (vacuum.evaluate_si(delta+hd, q)['Uvac_J_m3']
                        - vacuum.evaluate_si(delta-hd, q)['Uvac_J_m3'])/(2*hd)
    current_error = abs(2*E_CHARGE_C/HBAR_J_S*derivative_q/state['js_A_m2']-1)
    force_error = abs(derivative_delta/state['amplitude_force_per_m3']-1)
    assert max(current_error, force_error) < 1e-7
    nonthermal = loaded.evaluate_si(delta, q, population)
    nonthermal_dq = (loaded.evaluate_si(delta, q+hq, population)['energy_J_m3']
                    - loaded.evaluate_si(delta, q-hq, population)['energy_J_m3'])/(2*hq)
    nonthermal_current_error = abs(2*E_CHARGE_C/HBAR_J_S*nonthermal_dq/nonthermal['js_A_m2']-1)
    assert nonthermal_current_error < 1e-7
    assert loaded.evaluate_si(delta, 0, population)['js_A_m2'] == 0
    assert loaded.evaluate_si(delta, -q, population)['js_A_m2'] == -nonthermal['js_A_m2']
    result = {
        'status': 'passed', 'scope': 'independent numerical/API review; synthetic query populations',
        'implementation_sha256': hashlib.sha256((ROOT/'pysnspd/experimental/energy_catalog.py').read_bytes()).hexdigest(),
        'vacuum_reference_states': len(errors), 'max_vacuum_reference_scaled_error': max(errors),
        'spectral_parameter_combinations': spectral_states, 'minimum_DOS': min_dos,
        'maximum_scaled_spectral_normalization_residual': max_normalization,
        'roundtrip': 'arrays exact; eta, SI scales, provenance retained; extensionless filename works',
        'invalid_cases_rejected': rejected, 'normal_energy_factor4_error': abs(normal[0]-normal_reference),
        'nonthermal_derivative_identity_max_absolute_error': max(derivative_errors),
        'nonthermal_Maxwell_max_absolute_error': max(maxwell_errors),
        'Gamma_zero_one_sided_derivative_max_absolute_error': max(endpoint_errors),
        'normal_forces': list(normal[1:]), 'random_query_points': len(points),
        'minimum_excitation_energy': minimum, 'minimum_adjacent_energy_difference': spacing,
        'vacuum_SI_current_relative_derivative_error': current_error,
        'vacuum_SI_force_relative_derivative_error': force_error,
        'nonthermal_SI_current_relative_derivative_error': nonthermal_current_error,
        'runtime_seconds': time.perf_counter()-start,
        'limits': 'Does not replace grid/cutoff/eta convergence or validate phonon input and detector transients.',
    }
    Path(__file__).with_suffix('.json').write_text(json.dumps(result, indent=2)+'\n', encoding='utf-8')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
