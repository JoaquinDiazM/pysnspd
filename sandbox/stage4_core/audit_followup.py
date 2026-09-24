"""Audit two completed stage4A follow-up fields without new spectral queries.

All original user results remain untouched. Uses only JSON/NPZ data, geometric
quadrature, and the analytical derivatives of the prescribed input field.
"""
from pathlib import Path
import argparse
import hashlib
import json
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
REVIEW = ROOT/'docs/implementation/stage4/review_20260923'
DELIVERY = ROOT/'docs/implementation/stage4/followup_20260923'


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read(path):
    return json.loads(Path(path).read_text(encoding='utf8'))


def source_location(path, digest):
    for kind, root in [('current', ROOT),
                       ('followup_archive', DELIVERY/'previous_delivery_exact'),
                       ('review_archive', REVIEW/'previous_delivery_exact')]:
        candidate = root/path
        if candidate.is_file() and sha(candidate) == digest:
            return kind
    return 'unverified'


def metrics(directory, identity):
    result = read(directory/'result.json')
    with np.load(directory/'fields.npz') as archive:
        fields = {key: archive[key].copy() for key in archive.files}
    material = identity['material']
    mass = fields['quadrature_volume_m3']/(material['width_m']*material['thickness_m']*material['ell0_m'])
    shape = fields['shape']
    ij = np.indices(shape)
    boundary = ((ij[0] == 0)|(ij[1] == 0)|(ij[0] == shape[0]-1)|(ij[1] == shape[1]-1)).ravel()
    interior = ~boundary
    force = fields['cartesian_force']/mass[:, None]
    xx, yy = (fields['coordinates_m']/material['ell0_m']).T
    xx = xx-160e-9/material['ell0_m']/2
    phase_envelope = np.exp(-(xx/10)**2)
    phase_dx = .07-.0005*xx*np.sin(yy/5)*phase_envelope
    phase_dy = .005*np.cos(yy/5)*phase_envelope
    amplitude = abs(fields['delta'])
    scales = next(iter(result['mobility'].values()))['metadata']['scales']
    gap_ratio = scales['delta0_J']/(1.380649e-23*material['Tc_K'])
    exact_gamma = (amplitude**2/(amplitude**2+result['delta_reg']**2))**2*(phase_dx**2+phase_dy**2)/gap_ratio
    gamma_error = fields['gamma']-exact_gamma
    center = int(np.argmin(xx*xx+yy*yy))
    region_force = float(np.sqrt(np.sum(mass[interior]*np.sum(force[interior]**2, axis=1))))
    heat = {name: float(np.dot(mass[interior], fields['heat_'+name][interior]))
            for name in result['mobility']}
    smallest = fields['symbol_eigenvalues'][:, 0]
    uncertainty = fields['symbol_uncertainty']
    output = dict(nodes=int(len(mass)), elements=result['elements'], delta_reg=result['delta_reg'],
        energy_bar=result['energy_bar'], interior_force_L2=region_force,
        interior_heat=heat, symbol_minimum=float(smallest.min()),
        positive_symbols=int(np.count_nonzero(smallest > uncertainty)),
        negative_symbols=int(np.count_nonzero(smallest < -uncertainty)),
        unresolved_symbols=int(np.count_nonzero(abs(smallest) <= uncertainty)),
        gamma_relative_weighted_L2_error=float(np.sqrt(np.dot(mass, gamma_error**2)/np.dot(mass, exact_gamma**2))),
        gamma_maximum_absolute_error=float(max(abs(gamma_error))),
        center_gamma_discrete=float(fields['gamma'][center]),
        center_gamma_analytic=float(exact_gamma[center]),
        center_gamma_absolute_error=float(gamma_error[center]),
        center_gamma_ratio=float(fields['gamma'][center]/exact_gamma[center]),
        runtime_seconds=result['runtime_seconds'], result_sha256=sha(directory/'result.json'),
        fields_sha256=sha(directory/'fields.npz'))
    if 'derivative_diagnostics' in result:
        output['derivative_diagnostics'] = result['derivative_diagnostics']
    if 'section_currents' in result:
        currents = np.array(result['section_currents']['current_A'])
        output['section_current_summary_A'] = dict(minimum=float(currents.min()), maximum=float(currents.max()),
            leftmost=float(currents[0]), rightmost=float(currents[-1]),
            meaning='Nonstationary prescribed field section response, not an imposed detector bias or circuit solution')
    if result.get('mobility_fixed_boundary'):
        verification = {}
        for name, fixed in result['mobility_fixed_boundary'].items():
            free_velocity = fields['velocity_unconstrained_'+name]
            fixed_velocity = fields['velocity_fixed_boundary_'+name]
            fixed_heat = fields['heat_fixed_boundary_'+name]
            verification[name] = dict(boundary_velocity_max=float(max(abs(fixed_velocity[boundary]))),
                interior_velocity_max_change=float(max(abs(fixed_velocity[interior]-free_velocity[interior]))),
                boundary_heat=float(np.dot(mass[boundary], fixed_heat[boundary])),
                total_heat=float(np.dot(mass, fixed_heat)), recorded_total_heat=fixed['total_heat_bar'],
                heat_record_absolute_difference=float(abs(np.dot(mass, fixed_heat)-fixed['total_heat_bar'])),
                field_work=fixed['field_work_bar'], constraint_work=fixed['boundary_work_bar'],
                recorded_power_identity_residual=fixed['identity_residual_bar'])
        output['fixed_boundary_checks'] = verification
    return output


def relative_comparison(coarse, fine):
    return dict(energy_relative_change=(fine['energy_bar']-coarse['energy_bar'])/abs(coarse['energy_bar']),
        interior_force_L2_relative_change=fine['interior_force_L2']/coarse['interior_force_L2']-1,
        interior_heat_relative_change={name: fine['interior_heat'][name]/coarse['interior_heat'][name]-1
                                       for name in coarse['interior_heat']},
        gamma_weighted_error_coarse=coarse['gamma_relative_weighted_L2_error'],
        gamma_weighted_error_fine=fine['gamma_relative_weighted_L2_error'])


def calculate(raw):
    identity = read(raw/'identity.json')
    summary = read(raw/'summary.json')
    plan = REVIEW/'next_campaign_plan.json'
    plan_content = read(plan)
    source_records = {name: source_location(name, digest) for name, digest in identity['sources'].items()}
    checked = []
    for case in summary['cases']:
        for record in case['files']:
            path = raw/record['path']
            checked.append(dict(path=record['path'], sha256_match=sha(path)==record['sha256'],
                                size_match=path.stat().st_size==record['bytes']))
    new = {case['id']:metrics(raw/case['id'], identity) for case in summary['cases']}
    old_raw = REVIEW/'raw/stage4A_spatial_20260923'
    old_identity = read(old_raw/'identity.json')
    baseline_names = ('2d_suppressed_d0.05_4x2', '2d_suppressed_d0.1_4x2', '2d_suppressed_d0.1_8x4')
    baseline = {name: metrics(old_raw/name, old_identity) for name in baseline_names}
    d05fine = new['2d_suppressed_d0.05_8x4_fixed_boundary']
    d10finest = new['2d_suppressed_d0.1_16x8_fixed_boundary']
    d10fine = baseline['2d_suppressed_d0.1_8x4']
    return dict(schema='pysnspd.stage4A.followup_audit.v1', date='2026-09-23',
        source_commit='e8f23fab7a0431cc0917bf0ab3bf73684623e6ee',
        scope='Prescribed static non-photon fields; no trajectory or physical core admission',
        identity=dict(plan_sha256=sha(plan), plan_hash_match=sha(plan)==identity['plan_sha256'],
            sources=source_records, all_sources_verified=all(v!='unverified' for v in source_records.values()),
            cases_expected=len(plan_content['cases']), cases_completed=len(summary['cases']),
            task_ids_match=[case['id'] for case in summary['cases']]==identity['task_ids'],
            output_checks=checked, all_output_checks_pass=all(c['sha256_match'] and c['size_match'] for c in checked),
            runtime_seconds=summary['runtime_seconds'], summary_sha256=sha(raw/'summary.json'),
            identity_sha256=sha(raw/'identity.json')),
        new_cases=new, baseline_cases=baseline,
        comparisons=dict(delta_005_mesh_153_to_561=relative_comparison(baseline['2d_suppressed_d0.05_4x2'], d05fine),
            delta_010_mesh_153_to_561=relative_comparison(baseline['2d_suppressed_d0.1_4x2'], d10fine),
            delta_010_mesh_561_to_2145=relative_comparison(d10fine, d10finest),
            regularizer_010_to_005_same_561_nodes=relative_comparison(d10fine, d05fine)),
        diagnosis=dict(coarse_negative_center_removed=True,
            sampled_D36_positive_new_nodes=sum(case['positive_symbols'] for case in new.values()),
            all_new_D36_positive=all(case['negative_symbols']==case['unresolved_symbols']==0 for case in new.values()),
            bulk_numerical_resolution_sufficient_for_next_diagnostic=True,
            pointwise_center_gamma_fully_converged=False,
            no_new_mesh_only_batch_recommended=True,
            stage4_complete=False, physical_core_admitted=False, photon_admitted=False,
            production_promotion=False, time_trajectories=0, new_spectral_queries_in_audit=0),
        interpretation=[
            'The delta=.05 negative coarse center disappears at 561 nodes; its exact analytic reference is positive. This confirms the earlier spatial discretization diagnosis for this field, not universal physical stability.',
            'For delta=.10 the next mesh changes energy by about .009%, weighted interior force by .57%, and constrained Korzh heat by -2.55%; this is sufficient to move to a physical reference and controlled weak dynamics, without claiming a strict asymptotic error bound.',
            'The center Gamma is still 23.7% above an unusually small analytic value; the absolute difference and integrated error must also be reported. Do not require pointwise vanishing relative error before any further scientific work.',
            'At equal 561-node resolution changing delta=.10 to .05 changes Korzh heat by +21.2%. This is candidate-closure sensitivity in a prescribed synthetic population; it is not a measured uncertainty or a reason to calibrate delta from numerical positivity.',
            'D.4.4 still requires adequate spatial core reference and mobility/heat comparison; no vortex barrier, phase event, kinetic time identification or material detector validation is supplied by these snapshots.'
        ],
        next_action='Retain delta=.10 as the declared candidate; proceed to a same-ensemble physical spatial reference and an explicitly controlled weak non-photon transient/energy ledger. Do not schedule another blanket mesh refinement merely to tighten pointwise Gamma.',
        formula_conventions=dict(force='F_i=G_i/m_i; ||F||_L2=sqrt(sum_i m_i |F_i|^2)',
            gamma_error='sqrt(sum_i m_i (Gamma_i-Gamma_exact_i)^2 / sum_i m_i Gamma_exact_i^2)',
            heat='Interior weighted algebraic KWT response, equivalent to the explicitly clamped static boundary result',
            comparison='(new-old)/abs(old) for energy; new/old-1 for positive force and heat'),
        audit_script_sha256=sha(Path(__file__)))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--raw-root', type=Path, default=DELIVERY/'raw/stage4A_core_followup_20260923')
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError('Choose a new output; previous audits remain preserved')
    value = calculate(args.raw_root)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(value, indent=2, ensure_ascii=False, allow_nan=False)+'\n', encoding='utf8')
    print(json.dumps(dict(identity=value['identity'], comparisons=value['comparisons'], diagnosis=value['diagnosis']), indent=2))


if __name__ == '__main__':
    main()
