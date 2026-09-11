"""Apply frozen physical/API gates to a saved catalogue, without fitting it.

Usage: python assess_candidate.py CATALOG.npz --output assessment.json
Reference thermodynamics come from independent real causal-angle integration.
The complex-step check differentiates a local equivalent Hermite polynomial;
it tests implementation arithmetic, not spectral/field convergence. Very small
Gamma intervals are never tested by subtracting nearly equal total energies.
"""
from pathlib import Path
import argparse
import hashlib
import json
import math
import sys
import tempfile
import time
from unittest.mock import patch

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(HERE))
from gamma_zero_reference import population, reference as normal_reference
from causal_reference import integrate_profiles
from pysnspd.experimental.energy_catalog import OccupationEnergyCatalog, HBAR_J_S, E_CHARGE_C

CRITERIA = ROOT/'docs/implementation/stage1_r2/acceptance_criteria.json'


def local_complex_energy(table, amplitude, gamma, *, variable, step=1e-30):
    """Differentiate the fixed real cell, preserving its rounded coefficients."""
    a_axis = table.vacuum.delta_axis
    g_axis = getattr(table, '_parameter_axis', table.vacuum.gamma_axis)
    ratio_coordinates = getattr(table, 'gamma_ratio_axis', None) is not None
    real_parameter = gamma/amplitude if ratio_coordinates else gamma
    i = np.clip(np.searchsorted(a_axis, amplitude, side='right')-1, 0, len(a_axis)-2)
    j = np.clip(np.searchsorted(g_axis, real_parameter, side='right')-1, 0, len(g_axis)-2)
    a = amplitude+(1j*step if variable == 'amplitude' else 0)
    g = gamma+(1j*step if variable == 'Gamma' else 0)
    if ratio_coordinates:
        g = g/a  # vary amplitude at fixed physical Gamma, not fixed ratio
    offset = a-a_axis[i]
    def polynomial(coefficients):
        c = coefficients[:, i, j:j+2, :]
        return ((c[0]*offset+c[1])*offset+c[2])*offset+c[3]
    values = polynomial(table._log_coefficients)
    slopes = polynomial(table._gamma_log_coefficients)
    coordinate_scale = getattr(table, 'gamma_coordinate_scale', None)
    if coordinate_scale is not None and coordinate_scale > 0:
        gamma_coordinate = np.log1p(g/coordinate_scale)
        axis_coordinate = np.log1p(g_axis/coordinate_scale)
    else:
        gamma_coordinate, axis_coordinate = g, g_axis
    h = axis_coordinate[j+1]-axis_coordinate[j]
    t = (gamma_coordinate-axis_coordinate[j])/h
    difference = values[1]-values[0]
    # Algebraically equivalent local Hermite basis. Group the common endpoints
    # before dividing by a tiny Gamma interval, so a derivative of an O(1)
    # constant cannot contaminate the result by O(eps/h).
    log_increment = (values[0]+t*h*slopes[0]
                     +t*t*(3*difference-h*(2*slopes[0]+slopes[1]))
                     +t*t*t*(-2*difference+h*(slopes[0]+slopes[1])))
    if getattr(table, 'gamma_energy_offsets', None) is not None:
        # The v4 coefficients represent E-E_BCS directly. The exact finite-eta
        # BCS base is analytic in positive amplitude, including complex step.
        base_energy = table.count_nodes*np.sqrt(1+a*a/(table.count_nodes**2+table.eta**2))
        return np.imag(base_energy+log_increment)/step
    if getattr(table, '_base_log_coefficients', None) is not None:
        base = table._base_log_coefficients[:, i, :]
        log_increment += ((base[0]*offset+base[1])*offset+base[2])*offset+base[3]
    return np.imag(np.cumsum(np.exp(log_increment)))/step


def physical_row(values, reference, identifier, uncertainty, gates):
    values, reference = np.asarray(values), np.asarray(reference)
    absolute = abs(values-reference)
    scaled = absolute/np.maximum(1, abs(reference))
    relative = [None if abs(ref) <= uncertainty else float(err/abs(ref))
                for ref, err in zip(reference, absolute)]
    failures = []
    limits = [gates['energy_max_scaled_error'], gates['amplitude_force_max_scaled_error'],
              gates['gamma_response_max_scaled_error']]
    names = ['energy', 'amplitude_force', 'Gamma_response']
    for index, name in enumerate(names):
        if scaled[index] > limits[index]:
            failures.append(name+'_scaled')
        if index and relative[index] is not None and relative[index] > gates['nonzero_force_and_gamma_max_relative_error']:
            failures.append(name+'_relative')
        if index and relative[index] is None and identifier['kind'] != 'normal_point':
            failures.append(name+'_reference_zero_inconclusive')
        applicable = limits[index]*max(1, abs(reference[index]))
        if index and relative[index] is not None:
            applicable = min(applicable, gates['nonzero_force_and_gamma_max_relative_error']*abs(reference[index]))
        if uncertainty > applicable*gates['reference_error_budget_fraction_of_applicable_tolerance']:
            failures.append(name+'_reference_inconclusive')
    return dict(identifier, candidate=values.tolist(), reference=reference.tolist(),
                reference_error_bound=uncertainty, absolute_error=absolute.tolist(),
                scaled_error=scaled.tolist(), relative_error=relative, failures=failures)


def summarize(rows):
    result = {}
    for index, name in enumerate(('energy', 'amplitude_force', 'Gamma_response')):
        worst = max(rows, key=lambda row: row['scaled_error'][index])
        relative_rows = [row for row in rows if row['relative_error'][index] is not None]
        worst_relative = max(relative_rows, key=lambda row: row['relative_error'][index]) if relative_rows else None
        identifiers = ('kind', 'amplitude', 'gamma', 'profile', 'temperature')
        result[name] = {'max_scaled_error': worst['scaled_error'][index],
                        'worst_scaled_case': {key: worst[key] for key in identifiers if key in worst},
                        'max_relative_error': worst_relative['relative_error'][index] if worst_relative else None,
                        'worst_relative_case': {key: worst_relative[key] for key in identifiers if key in worst_relative} if worst_relative else None}
    return result


def assess(path):
    start = time.perf_counter()
    source_at_start = hashlib.sha256((ROOT/'pysnspd/experimental/energy_catalog.py').read_bytes()).hexdigest()
    criteria = json.loads(CRITERIA.read_text(encoding='utf-8'))
    gates = criteria['numerical_gates']
    references = json.loads((HERE/'causal_reference.json').read_text(encoding='utf-8'))
    thermal = json.loads((HERE/'thermal_reference.json').read_text(encoding='utf-8'))
    table = OccupationEnergyCatalog.load(path)
    rows, exceptions, derivative_rows = [], [], []
    for ref in references['rows']:
        a, g, profile = ref['amplitude'], ref['gamma'], ref['profile']
        identifier = dict(kind='fixed_occupation', amplitude=a, gamma=g, profile=profile)
        try:
            p = np.array([population(profile, x) for x in table.count_nodes])
            value = table.evaluate(a, g, p)
            uncertainty = max(ref['reference_error_estimate'], max(ref['cutoff_and_quadrature_refinement_absolute_difference']), 5e-13)
            rows.append(physical_row(value, ref['reference'], identifier, uncertainty, gates))
        except Exception as exc:
            exceptions.append(dict(identifier, exception=f'{type(exc).__name__}: {exc}'))
    # Deterministic extra points stay unchanged when field grids are refined.
    # They match the builder's supplemental samples but use an independent
    # ideal causal reference, not the builder's finite-eta spectral evaluator.
    extra_rng = np.random.default_rng(7431)
    names = [profile['id'] for profile in criteria['mandatory_population_profiles']]
    extra_references = []
    for index in range(12):
        a = float(extra_rng.uniform(.085, 1.49))
        g = float(10**extra_rng.uniform(-12, -1.5) if index < 6 else extra_rng.uniform(.031, 1.19))
        expected, estimate = integrate_profiles(a, g, names, cutoff=20., tolerance=2e-12)
        for profile_index, name in enumerate(names):
            identifier = dict(kind='offnode_fixed_occupation', amplitude=a, gamma=g, profile=name)
            extra_references.append(dict(identifier, reference=expected[profile_index].tolist(), reference_error_estimate=estimate))
            try:
                p = np.array([population(name, x) for x in table.count_nodes])
                value = table.evaluate(a, g, p)
                rows.append(physical_row(value, expected[profile_index], identifier, max(estimate, 5e-13), gates))
            except Exception as exc:
                exceptions.append(dict(identifier, exception=f'{type(exc).__name__}: {exc}'))
    for ref in thermal['rows']:
        a, g, temperature = ref['amplitude'], ref['gamma'], ref['temperature']
        identifier = dict(kind='thermal_equilibrium', amplitude=a, gamma=g, temperature=temperature)
        try:
            energy = table.energy_kernel(a, g)[0]
            p = np.exp(-energy/temperature)/(1+np.exp(-energy/temperature))
            result = np.array(table.evaluate(a, g, p))
            # -T*S = 4*T*sum[p log p+(1-p) log(1-p)] w, with 0 log 0=0.
            entropy_term = np.zeros_like(p)
            interior = (p > 0)&(p < 1)
            entropy_term[interior] = p[interior]*np.log(p[interior])+(1-p[interior])*np.log1p(-p[interior])
            result[0] += 4*temperature*np.dot(entropy_term, table.count_weights)
            uncertainty = max(ref['causal_quadrature_error_estimate'], max(ref['matsubara_tail_extrapolation_change']), 5e-13)
            rows.append(physical_row(result, ref['causal_reference'], identifier, uncertainty, gates))
        except Exception as exc:
            exceptions.append(dict(identifier, exception=f'{type(exc).__name__}: {exc}'))
    for g in (0., .4, 1.2):
        for profile in criteria['mandatory_population_profiles']:
            name = profile['id']
            p = np.array([population(name, x) for x in table.count_nodes])
            reference, uncertainty = normal_reference(0., name)
            result = table.evaluate(0., g, p)
            rows.append(physical_row(result, reference, dict(kind='normal_point', amplitude=0., gamma=g, profile=name),
                                     max(max(uncertainty), 5e-13), gates))
    for a, g in criteria['mandatory_field_points']:
        try:
            kernels = table.energy_kernel(a, g)
            derivatives = [local_complex_energy(table, a, g, variable=variable) for variable in ('amplitude', 'Gamma')]
            for profile in criteria['mandatory_population_profiles']:
                p = np.array([population(profile['id'], x) for x in table.count_nodes])
                actual = np.array([4*np.dot(kernel*p, table.count_weights) for kernel in kernels[1:]])
                expected = np.array([4*np.dot(kernel*p, table.count_weights) for kernel in derivatives])
                err = abs(actual-expected)/np.maximum(1, abs(expected))
                derivative_rows.append(dict(amplitude=a, gamma=g, profile=profile['id'],
                                            candidate_excitation_derivative=actual.tolist(),
                                            complex_step_excitation_derivative=expected.tolist(), scaled_error=err.tolist()))
        except Exception as exc:
            exceptions.append(dict(kind='complex_step_identity', amplitude=a, gamma=g, exception=f'{type(exc).__name__}: {exc}'))
    invalid = []
    p = np.full_like(table.count_nodes, .2)
    bad_cases = [('p_nan', .72, 0, p*np.nan), ('p_inf', .72, 0, p*np.inf),
                 ('p_negative', .72, 0, p*0-1e-15), ('p_above_one', .72, 0, p*0+1+1e-15),
                 ('p_wrong_shape', .72, 0, p[:-1]), ('amplitude_below_support', .04, 0, p),
                 ('amplitude_above_support', 1.6, 0, p), ('negative_gamma', .72, -1e-20, p),
                 ('gamma_above_support', .72, 1.3, p), ('nonfinite_amplitude', math.nan, 0, p),
                 ('nonfinite_gamma', .72, math.inf, p)]
    for name, a, g, occupation in bad_cases:
        try:
            table.evaluate(a, g, occupation)
            invalid.append(dict(case=name, rejected=False))
        except ValueError:
            invalid.append(dict(case=name, rejected=True))
    pauli_endpoints_accepted = True
    for occupation in (np.zeros_like(p), np.ones_like(p), (np.arange(len(p))%2).astype(float)):
        try:
            pauli_endpoints_accepted &= bool(np.all(np.isfinite(table.evaluate(.72, 0, occupation))))
        except Exception:
            pauli_endpoints_accepted = False
    support_rows = []
    rng = np.random.default_rng(831317)
    fields = [(a, g) for a, g in criteria['mandatory_field_points']]
    fields += list(zip(rng.uniform(.08, 1.5, 120), np.r_[rng.uniform(0, 1.2, 60), 10**rng.uniform(-13, -2, 60)]))
    regression_path = HERE/'support_regression_points.json'
    regression_points = set()
    if regression_path.is_file():
        historical = json.loads(regression_path.read_text(encoding='utf-8'))
        regression_points = {(point['amplitude'], point['gamma']) for point in historical['points']}
        fields += sorted(regression_points)
    for a, g in fields:
        try:
            energy = table.energy_kernel(a, g)[0]
            support_rows.append(dict(amplitude=float(a), gamma=float(g), minimum_energy=float(min(energy)),
                                     minimum_adjacent_spacing=float(min(np.diff(energy))), finite=bool(np.all(np.isfinite(energy))),
                                     historical_regression=(a, g) in regression_points))
        except Exception as exc:
            exceptions.append(dict(kind='support', amplitude=float(a), gamma=float(g), exception=f'{type(exc).__name__}: {exc}'))
    dense_support = dict(enabled=table.gamma_energy_offsets is not None, queries=0, failures=[],
                         minimum_energy=None, minimum_adjacent_spacing=None,
                         sampling='Three deterministic interior diagonal pairs at fractions 0.2113248654, 0.5, 0.7886751346 in every represented field cell. Use the actual interpolation coordinate; ratio cells are intersected with physical Gamma support before sampling.')
    if dense_support['enabled']:
        point_digest = hashlib.sha256()
        parameter_axis = getattr(table, '_parameter_axis', table.vacuum.gamma_axis)
        ratio_coordinates = getattr(table, 'gamma_ratio_axis', None) is not None
        for a0, a1 in zip(table.vacuum.delta_axis[:-1], table.vacuum.delta_axis[1:]):
            for g0, g1 in zip(parameter_axis[:-1], parameter_axis[1:]):
                for fraction in (.21132486540518713, .5, .7886751345948129):
                    a = float(a0+fraction*(a1-a0))
                    upper = min(g1, table.vacuum.gamma_axis[-1]/a) if ratio_coordinates else g1
                    if upper <= g0:
                        continue
                    coordinate_scale = table.gamma_coordinate_scale
                    parameter = float(g0+fraction*(upper-g0) if coordinate_scale is None else
                                      g0+(g0+coordinate_scale)*np.expm1(fraction*np.log1p((upper-g0)/(g0+coordinate_scale))))
                    g = a*parameter if ratio_coordinates else parameter
                    point_digest.update(np.asarray([a, g], dtype='<f8').tobytes())
                    dense_support['queries'] += 1
                    try:
                        energy = table.energy_kernel(a, g)[0]
                        smallest, spacing = float(min(energy)), float(min(np.diff(energy)))
                        if dense_support['minimum_energy'] is None or smallest < dense_support['minimum_energy']:
                            dense_support['minimum_energy'] = smallest
                            dense_support['minimum_energy_case'] = dict(amplitude=a, gamma=g)
                        if dense_support['minimum_adjacent_spacing'] is None or spacing < dense_support['minimum_adjacent_spacing']:
                            dense_support['minimum_adjacent_spacing'] = spacing
                            dense_support['minimum_spacing_case'] = dict(amplitude=a, gamma=g)
                    except Exception as exc:
                        dense_support['failures'].append(dict(amplitude=a, gamma=g, exception=f'{type(exc).__name__}: {exc}'))
        dense_support['sample_coordinates_sha256'] = point_digest.hexdigest()
    temp = Path(tempfile.mkdtemp(prefix='stage1-r2-independent-', dir=ROOT/'tmp'))
    saved = table.save(temp/'catalog_without_extension')
    with patch('pysnspd.experimental.energy_catalog.retarded_spectrum', side_effect=AssertionError('reload resolved spectra')), \
            patch('pysnspd.experimental.energy_catalog.retarded_spectrum_batch', side_effect=AssertionError('reload resolved spectra')):
        reloaded = OccupationEnergyCatalog.load(saved)
        roundtrip = all(np.array_equal(getattr(table, key), getattr(reloaded, key))
                        for key in ('count_nodes', 'count_weights', 'excitation_energies', 'gamma_energy_derivatives'))
        for key in ('gamma_log_offsets', 'gamma_coordinate_scale', 'gamma_energy_offsets', 'delta_energy_derivatives', 'gamma_ratio_axis'):
            if hasattr(table, key):
                roundtrip &= np.array_equal(getattr(table, key), getattr(reloaded, key))
        def query_signature(catalog, amplitude, gamma):
            try:
                return ('value', catalog.evaluate(amplitude, gamma, p))
            except Exception as exc:
                return ('rejection', type(exc).__name__, str(exc))
        roundtrip &= all(query_signature(table, a, g) == query_signature(reloaded, a, g) for a, g in fields[:30])
    metadata_retained = all(reloaded.vacuum.metadata[key] == value for key, value in table.vacuum.metadata.items())
    scales_retained = all(getattr(table.vacuum, key) == getattr(reloaded.vacuum, key)
                          for key in ('delta0_J', 'N0_per_J_m3', 'D_m2_s')) and table.eta == reloaded.eta
    incompatible = temp/'incompatible_schema.npz'
    np.savez_compressed(incompatible, metadata_json=np.array(json.dumps({'schema':'pysnspd.experimental.occupation_catalog.v1'})))
    incompatible_schema_rejected = False
    try:
        OccupationEnergyCatalog.load(incompatible)
    except ValueError:
        incompatible_schema_rejected = True
    def derivative(function, coordinate, h=2e-6):
        return (function(coordinate-2*h)-8*function(coordinate-h)+8*function(coordinate+h)-function(coordinate+2*h))/(12*h)
    maxwell_rows = []
    for a, g in ((.23, .045), (.67, .431), (1.23, .851)):
        p = np.array([population('low_energy', x) for x in table.count_nodes])
        try:
            ag = derivative(lambda coordinate: table.evaluate(coordinate, g, p)[2], a)
            ga = derivative(lambda coordinate: table.evaluate(a, coordinate, p)[1], g)
            maxwell_rows.append(dict(amplitude=a, gamma=g, derivative_amplitude_Gamma_response=float(ag),
                                     derivative_Gamma_amplitude_force=float(ga),
                                     scaled_error=float(abs(ag-ga)/max(1, abs(ag), abs(ga)))))
        except Exception as exc:
            exceptions.append(dict(kind='Maxwell', amplitude=a, gamma=g, exception=f'{type(exc).__name__}: {exc}'))
    a, g = .67, .431
    scale, N0, D = table.vacuum.delta0_J, table.vacuum.N0_per_J_m3, table.vacuum.D_m2_s
    q = math.sqrt(2*g*scale/(HBAR_J_S*D))
    si = table.evaluate_si(a*scale, q, p)
    normalized = table.evaluate(a, g, p)
    expected_si = dict(energy_J_m3=N0*scale*scale*normalized[0], amplitude_force_per_m3=N0*scale*normalized[1],
                       Pi_J_m2=N0*scale*HBAR_J_S*D*q*normalized[2],
                       js_A_m2=2*E_CHARGE_C*N0*scale*D*q*normalized[2], Gamma_J=g*scale)
    conversion_error = max(abs(si[key]/value-1) for key, value in expected_si.items() if value != 0)
    d_energy_q = derivative(lambda coordinate: table.evaluate_si(a*scale, coordinate, p)['energy_J_m3'], q, h=q*1e-5)
    current_derivative_error = abs(2*E_CHARGE_C/HBAR_J_S*d_energy_q/si['js_A_m2']-1)
    current_parity = table.evaluate_si(a*scale, -q, p)['js_A_m2'] == -si['js_A_m2']
    current_zero = table.evaluate_si(a*scale, 0, p)['js_A_m2'] == 0
    si_result = dict(max_conversion_relative_error=float(conversion_error),
                     current_energy_derivative_relative_error=float(current_derivative_error),
                     odd_signed_q_current=bool(current_parity), zero_q_current=bool(current_zero))
    reference_hashes = {name: hashlib.sha256((HERE/name).read_bytes()).hexdigest()
                        for name in ('causal_reference.json', 'thermal_reference.json', 'gamma_zero_reference.py')}
    failures = [row for row in rows if row['failures']]
    worst_derivative = max(derivative_rows, key=lambda row: max(row['scaled_error'])) if derivative_rows else None
    identity_failed = worst_derivative is None or max(worst_derivative['scaled_error']) > gates['derivative_identity_max_scaled_error']
    passed = not (failures or exceptions or identity_failed or not roundtrip or not metadata_retained or
                  not scales_retained or not all(row['rejected'] for row in invalid) or
                  bool(dense_support['failures']) or
                  not incompatible_schema_rejected or not pauli_endpoints_accepted or
                  max((row['scaled_error'] for row in maxwell_rows), default=math.inf) > gates['Maxwell_identity_max_scaled_error'] or
                  max(conversion_error, current_derivative_error) > gates['derivative_identity_max_scaled_error'] or
                  not current_parity or not current_zero)
    source_at_end = hashlib.sha256((ROOT/'pysnspd/experimental/energy_catalog.py').read_bytes()).hexdigest()
    result = dict(schema='pysnspd.stage1_r2.independent_candidate_assessment.v1',
                  status=('INCONCLUSIVE_SOURCE_CHANGED' if source_at_start != source_at_end else 'PASS_TESTED_PHYSICAL_AND_API_GATES' if passed else 'FAIL'),
                  scope='Mandatory causal thermodynamics and API checks only; full acceptance additionally requires separated interpolation/quadrature/cutoff/eta and causal solver evidence.',
                  criteria_sha256=hashlib.sha256(CRITERIA.read_bytes()).hexdigest(),
                  catalog_path=str(path), catalog_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                  implementation_sha256=source_at_start, implementation_sha256_at_end=source_at_end,
                  reviewer_code_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                  reference_hashes=reference_hashes, eta=table.eta,
                  counts=dict(physical_cases=len(rows), derivative_cases=len(derivative_rows), support_cases=len(support_rows)),
                  historical_support_regression_points=len(regression_points),
                  support_regression_points_sha256=hashlib.sha256(regression_path.read_bytes()).hexdigest() if regression_path.is_file() else None,
                  summary=summarize(rows), failed_physical_cases=failures, exceptions=exceptions,
                  physical_rows=rows, supplemental_offnode_causal_references=extra_references,
                  derivative_rows=derivative_rows,
                  derivative_method='Complex-step 1e-30 of an algebraically equivalent local Hermite polynomial with rounded stored coefficients. Tests differentiation of that representation, not its physical accuracy; no subtraction of total energies over tiny Gamma intervals.',
                  worst_derivative_identity=worst_derivative, derivative_identity_passed=not identity_failed,
                  support_rows=support_rows, invalid_input_checks=invalid,
                  dense_cell_support=dense_support,
                  pauli_endpoints_accepted=bool(pauli_endpoints_accepted), incompatible_schema_rejected=bool(incompatible_schema_rejected),
                  Maxwell_rows=maxwell_rows, SI_checks=si_result,
                  roundtrip_arrays_and_queries_exact=bool(roundtrip), roundtrip_metadata_retained=bool(metadata_retained),
                  roundtrip_scales_eta_retained=bool(scales_retained), extensionless_filename_exists=saved.is_file(),
                  runtime_seconds=time.perf_counter()-start)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('catalog', type=Path)
    parser.add_argument('--output', type=Path, required=True)
    arguments = parser.parse_args()
    result = assess(arguments.catalog)
    arguments.output.parent.mkdir(parents=True, exist_ok=True)
    arguments.output.write_text(json.dumps(result, indent=2)+'\n', encoding='utf-8')
    print(json.dumps({key: result[key] for key in ('status', 'eta', 'counts', 'summary', 'worst_derivative_identity', 'runtime_seconds')}, indent=2))


if __name__ == '__main__':
    main()
