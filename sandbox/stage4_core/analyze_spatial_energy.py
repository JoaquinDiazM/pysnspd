"""Read-only postprocessing of the frozen thermal spatial campaign.

No spectral solves, source changes, new mesh calculation or time evolution.
The original runner's radial mask excludes the origin; this analysis includes
the origin in Cartesian force norms and records its finite force separately.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT/'docs/implementation/stage4/spatial_energy_20260924'
DEFAULT_RAW = DATA/'raw/stage4_spatial_energy_reference_20260924'


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def complex_pair(value):
    return [float(np.real(value)), float(np.imag(value))]


def l2(values, weights, mask):
    return float(np.sqrt(np.sum(weights[mask]*np.abs(values[mask])**2)))


def nodal_current(arrays):
    """Cartesian-only diagnostic: mean adjacent face current densities.

    Edge current is integrated over its dual face, of length w_ij*l_ij.
    This reconstruction is for comparisons, not an additional solver field.
    """
    xy, edges = arrays['coordinates_bar'], arrays['edges']
    tail, head = edges.T
    delta = xy[head]-xy[tail]
    length = np.linalg.norm(delta, axis=1)
    direction = delta/length[:, None]
    if not np.all(np.isclose(np.max(np.abs(direction), axis=1), 1.)):
        raise ValueError('This current reconstruction requires Cartesian edges')
    face = arrays['conductance']*length
    active = face > 0
    density = np.divide(arrays['current_bar'], face,
                        out=np.zeros(len(face)), where=active)
    numerator, denominator = np.zeros_like(xy), np.zeros_like(xy)
    for end in (tail, head):
        np.add.at(numerator, end, density[:, None]*direction)
        np.add.at(denominator, end, active[:, None]*np.abs(direction))
    answer = np.divide(numerator, denominator, out=np.zeros_like(xy), where=denominator > 0)
    return answer[:, 0]+1j*answer[:, 1]


def process(raw):
    summary = json.loads((raw/'summary.json').read_text(encoding='utf8'))
    identity = json.loads((raw/'identity.json').read_text(encoding='utf8'))
    plan_path = DATA/'campaign_plan.json'
    plan = json.loads(plan_path.read_text(encoding='utf8'))
    if sha(plan_path) != identity['plan_sha256']:
        raise ValueError('Plan differs from the frozen campaign')
    source_records = {}
    for name, expected in identity['sources'].items():
        actual = sha(ROOT/name)
        if actual != expected:
            raise ValueError('Campaign source no longer matches identity: '+name)
        source_records[name] = actual
    expected_jobs = {(case['id'], n) for case in plan['cases'] for n in range(case['matsubara_count'])}
    observed_jobs = [(record['case_id'], record['n']) for record in summary['records']]
    if set(observed_jobs) != expected_jobs or len(observed_jobs) != len(expected_jobs):
        raise ValueError('Missing, duplicate or unexpected spectral record')
    if not all(record['metadata']['converged'] for record in summary['records']):
        raise ValueError('Unconverged spectral record')
    fields, hashes, cases = {}, {}, {}
    for cid, case in summary['cases'].items():
        cases[cid] = {}
        for count, previous in case['counts'].items():
            path = raw/previous['fields_file']
            actual = sha(path)
            if actual != previous['fields_sha256']:
                raise ValueError('Map hash mismatch: '+str(path))
            hashes[path.name] = actual
            with np.load(path) as stored:
                arrays = {key: stored[key].copy() for key in stored.files}
            if any(np.any(~np.isfinite(value)) for value in arrays.values()):
                raise ValueError('Nonfinite map: '+path.name)
            fields[cid, count] = arrays
            xy, mass, force = arrays['coordinates_bar'], arrays['area_weights'], arrays['force']
            radius = np.linalg.norm(xy, axis=1)
            core = radius <= 4.
            center = np.flatnonzero(radius == 0.)
            if len(center) != 1:
                raise ValueError('Expected exactly one origin node')
            center = int(center[0])
            phase_mask = core & (abs(arrays['d']) > 0) & (abs(arrays['f_lowest']) > 0)
            phase = np.angle(arrays['f_lowest'][phase_mask]*np.conj(arrays['d'][phase_mask]))
            peak = int(np.flatnonzero(phase_mask)[np.argmax(abs(phase))])
            direction = np.divide(xy[:, 0]+1j*xy[:, 1], radius,
                                  out=np.zeros(len(radius), complex), where=radius > 0)
            current = nodal_current(arrays)
            force_norm = l2(force, mass, core)
            record = dict(
                energy_finite_sum=previous['energy'],
                force_core_L2_including_origin=force_norm,
                raw_force_core_L2_excluding_origin=previous['core_force_L2'],
                current_core_L2_reconstructed=l2(current, mass, core),
                center_gap=complex_pair(arrays['d'][center]),
                center_force=complex_pair(force[center]),
                center_force_magnitude=float(abs(force[center])),
                center_lowest_spectral_f=complex_pair(arrays['f_lowest'][center]),
                center_lowest_spectral_g=float(arrays['g_lowest'][center]),
                center_dual_area=float(mass[center]),
                phase_lowest_max_abs_rad=float(np.max(abs(phase))),
                phase_lowest_max_abs_deg=float(np.rad2deg(np.max(abs(phase)))),
                phase_lowest_area_weighted_rms_rad=float(np.sqrt(np.sum(mass[phase_mask]*phase**2)/np.sum(mass[phase_mask]))),
                phase_peak_coordinates_ell0=xy[peak].tolist(),
                phase_peak_gap_magnitude=float(abs(arrays['d'][peak])),
                phase_peak_f_magnitude=float(abs(arrays['f_lowest'][peak])),
                tangential_force_L2_excluding_origin=l2(np.imag(np.conj(direction)*force), mass, core & (radius > 0)),
                interior_condensate_torque_plus_current_divergence_max=previous['maximum_interior_gauge_identity_residual'],
            )
            if case['case']['profile'] == 'radial_vortex':
                reference = arrays['reference_radial_force']
                record['radial_force_relative_L2_including_origin'] = l2(force-reference, mass, core)/l2(reference, mass, core)
            cases[cid][count] = record
    cutoffs = {}
    for cid in ('radial_33', 'radial_65', 'radial_129', 'asymmetric_65'):
        low, high = fields[cid, '128'], fields[cid, '256']
        core = np.linalg.norm(high['coordinates_bar'], axis=1) <= 4.
        mass = high['area_weights']
        jlow, jhigh = nodal_current(low), nodal_current(high)
        elo = cases[cid]['128']['energy_finite_sum']
        ehi = cases[cid]['256']['energy_finite_sum']
        cutoffs[cid] = dict(
            force_relative_L2=l2(high['force']-low['force'], mass, core)/l2(high['force'], mass, core),
            current_relative_L2_reconstructed=l2(jhigh-jlow, mass, core)/l2(jhigh, mass, core),
            energy_change=ehi-elo, energy_relative_change=abs(ehi-elo)/abs(ehi),
            interpretation='Observed raw finite-sum change N128 to N256, not a rigorous remainder bound')
    mesh = {}
    for small, large in ((33, 65), (65, 129)):
        low, high = fields[f'radial_{small}', '256'], fields[f'radial_{large}', '256']
        factor = (large-1)//(small-1)
        sampled_force = high['force'].reshape(large, large)[::factor, ::factor].ravel()
        sampled_current = nodal_current(high).reshape(large, large)[::factor, ::factor].ravel()
        core = np.linalg.norm(low['coordinates_bar'], axis=1) <= 4.
        mass = low['area_weights']
        elo, ehi = cases[f'radial_{small}']['256']['energy_finite_sum'], cases[f'radial_{large}']['256']['energy_finite_sum']
        mesh[f'{small}_to_{large}'] = dict(
            force_relative_L2_on_coarse_nodes=l2(low['force']-sampled_force, mass, core)/l2(sampled_force, mass, core),
            current_relative_L2_on_coarse_nodes=l2(nodal_current(low)-sampled_current, mass, core)/l2(sampled_current, mass, core),
            radial_error_reduction_factor=cases[f'radial_{small}']['256']['radial_force_relative_L2_including_origin']/cases[f'radial_{large}']['256']['radial_force_relative_L2_including_origin'],
            energy_change=ehi-elo, energy_relative_change=abs(ehi-elo)/abs(ehi))
    directional = {}
    for check in plan['directional_checks']:
        count = str(check['count'])
        center = summary['cases'][check['base']]['counts'][count]
        high = summary['cases'][check['plus']]['counts'][count]
        low = summary['cases'][check['minus']]['counts'][count]
        measured = (high['energy']-low['energy'])/(2*check['step'])
        predicted = center[check['kind']+'_directional_work']
        directional[check['kind']] = dict(
            energy_directional_derivative=measured,force_directional_work=predicted,
            absolute_difference=abs(measured-predicted),relative_difference=abs(measured-predicted)/abs(predicted),
            perturbation_step=check['step'],spectral_boundary_trace_fixed=True,
            same_finite_sum=True)
    return dict(
        schema='pysnspd.stage4.thermal_spatial_postprocessing.v1',
        provenance=dict(raw_directory=raw.relative_to(ROOT).as_posix(),
            summary_sha256=sha(raw/'summary.json'),identity_sha256=sha(raw/'identity.json'),
            plan_sha256=sha(plan_path),analysis_script_sha256=sha(__file__),
            sources_matching_identity=source_records,map_sha256=hashes,
            locally_verified_maps=len(hashes),expected_spectral_records=len(expected_jobs),
            spectral_records_unique_complete=True,all_records_converged=True,
            mode_files='2048 individual mode files retained on Geminga; record hashes present in summary, mode contents not reverified locally'),
        execution=dict(runtime_seconds=summary['runtime_seconds'],
            maximum_spectral_residual=max(record['metadata']['residual'] for record in summary['records']),
            maximum_newton_iterations=max(record['metadata']['iterations'] for record in summary['records']),
            maximum_worker_seconds=max(record['worker_seconds'] for record in summary['records']),
            process_budget=identity['budget'],physical_time_steps=0),
        cases=cases,raw_cutoff_changes=cutoffs,mesh_changes=mesh,directional_derivatives=directional,
        interpretation=dict(
            admitted='Consistent finite-sum thermal Cartesian spatial energy oracle for prescribed fields and fixed spectral contacts',
            spectral_phase='The anomalous Matsubara propagator phase is solved independently of the condensate phase; it is not a new independent dynamical quantum state',
            current_reconstruction='Cartesian arithmetic average of adjacent face current densities I/(w*edge_length); used only for weighted diagnostic comparisons',
            raw_origin_policy='Frozen runner excluded origin from core norm; this postprocessing includes it and preserves both values',
            next_action='Adopt this experimental auxiliary-spectral energy as the thermal spatial reference; derive and verify any tail from the same action, then connect the dissipative/kinetic contract. Do not fit another q_delta or refine a rejected local closure.',
            limitations=['Finite Matsubara sum still has a measurable cutoff dependence',
                'Square 12 ell0 local core test, not the full Korzh detector geometry',
                'Prescribed gap profile, not a self-consistent stationary vortex or activation barrier',
                'Thermal oracle does not establish a closure for nonthermal occupations',
                'No mobility, time integration, photon, normal current or circuit evolution has been validated'],
            stage4_complete=False,production_solver_changed=False))


def markdown(result):
    lines = ['# Energía térmica espacial: análisis independiente', '',
        'El campo espectral espacial resuelve el desacuerdo de la fuerza del núcleo del cierre local anterior. La prueba admite este oráculo térmico experimental para campos prescritos; la etapa 4 continúa abierta.', '',
        f"Se completaron {result['provenance']['expected_spectral_records']} problemas espectrales en {result['execution']['runtime_seconds']:.3f} s. Los {result['provenance']['locally_verified_maps']} mapas finales coinciden con sus SHA-256 y las fuentes con la identidad de ejecución. Los archivos individuales por modo permanecen en Geminga; su contenido no se volvió a verificar localmente.", '',
        '## Comparación con la referencia radial independiente', '',
        '| Nodos | Error L2 de fuerza frente al BVP, N=256 | Cambio de fuerza N128→256 |',
        '|---:|---:|---:|']
    for nodes in (33, 65, 129):
        cid = f'radial_{nodes}'
        lines.append(f"| {nodes}² | {100*result['cases'][cid]['256']['radial_force_relative_L2_including_origin']:.5f}% | {100*result['raw_cutoff_changes'][cid]['force_relative_L2']:.5f}% |")
    lines += ['', 'El error de malla disminuye aproximadamente cuatro veces al dividir el paso por dos. Las comparaciones usan el mismo corte en la malla y en el BVP, de modo que esta convergencia no elimina la dependencia de corte. La norma cartesiana incluye el origen y utiliza áreas duales; se conserva por separado la norma histórica que lo excluía.', '',
        '## Energía, fuerza y fase espectral', '',
        '| Variación compacta | Derivada de energía | Trabajo de fuerza | Diferencia relativa |',
        '|---|---:|---:|---:|']
    for key, record in result['directional_derivatives'].items():
        lines.append(f"| {key} | {record['energy_directional_derivative']:.9g} | {record['force_directional_work']:.9g} | {record['relative_difference']:.3g} |")
    asymmetric = result['cases']['asymmetric_65']['256']
    lines += ['', 'Estas variaciones mantienen la traza espectral de contacto fija y vuelven a resolver el interior. No se añadió una cola solo a la fuerza ni un término K0.', '',
        f"En el perfil asimétrico, la fase del propagador anómalo del modo más bajo difiere hasta {asymmetric['phase_lowest_max_abs_deg']:.4f}° de la fase local del condensado. El máximo ocurre en x,y={asymmetric['phase_peak_coordinates_ell0']} ℓ0, con |Δ|/(kBTc)={asymmetric['phase_peak_gap_magnitude']:.6g} y |f₀|={asymmetric['phase_peak_f_magnitude']:.6g}. La fase espectral es una variable auxiliar estacionaria; no se está introduciendo otra variable dinámica cuántica.", '',
        f"El centro tiene Δ=0 y fuerza cartesiana finita {asymmetric['center_force']} en unidades adimensionales. Allí f₀={asymmetric['center_lowest_spectral_f']}; no existe una fase local del gap con la cual comparar. Este punto se conserva en la nueva norma completa, {asymmetric['force_core_L2_including_origin']:.8g}, frente a {asymmetric['raw_force_core_L2_excluding_origin']:.8g} en el registro anterior que excluía el origen.", '',
        '## Decisión y trabajo pendiente', '',
        'Usar la energía espacial con campos espectrales auxiliares como referencia térmica experimental de la continuación. Su acción finita produce fuerza y corriente compatibles y converge hacia el BVP independiente. No se justifica ajustar otro regularizador del momento ni invertir más malla en el cierre local descartado.', '',
        'El siguiente trabajo es obtener cualquier aceleración de la cola desde esta misma energía y conectar explícitamente el contrato disipativo/cinético. El cambio observado con el corte orienta ese esfuerzo. Estos resultados no validan la movilidad temporal, un transiente, el transporte no térmico, un fotón, el circuito, un vórtice autoconsistente ni el dispositivo Korzh completo.', '',
        'El archivo analysis.json contiene las cifras completas, la reconstrucción diagnóstica de corriente, la comparación de mallas, los datos del origen y la procedencia verificable. El postproceso no ejecutó nuevos problemas espectrales.', '']
    return '\n'.join(lines)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--raw', type=Path, default=DEFAULT_RAW)
    parser.add_argument('--output-directory', type=Path, default=DATA)
    args = parser.parse_args()
    result = process(args.raw)
    args.output_directory.mkdir(parents=True, exist_ok=True)
    with (args.output_directory/'analysis.json').open('x', encoding='utf8') as stream:
        json.dump(result, stream, indent=2, ensure_ascii=False, allow_nan=False)
        stream.write('\n')
    with (args.output_directory/'analysis.md').open('x', encoding='utf8') as stream:
        stream.write(markdown(result))
    print(json.dumps(dict(status='POSTPROCESSED_NO_SOLVES',maps=result['provenance']['locally_verified_maps'],
        output_directory=str(args.output_directory)),ensure_ascii=False))


if __name__ == '__main__':
    main()
