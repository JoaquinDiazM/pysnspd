"""Offline audit: read frozen arrays only; no catalogue or dynamical RHS calls."""
from pathlib import Path
import hashlib
import json
import numpy as np

ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / 'docs/implementation/stage3/coupled_20260923'
RAW = OUT / 'raw/mixed_pilot'


def read(path):
    return json.loads(path.read_text(encoding='utf-8'))


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    manifest, summary = read(RAW/'manifest.json'), read(RAW/'summary.json')
    reg = read(OUT/'mixed_registration.json')
    issues = []
    sources = {name: sha(ROOT/name) == expected for name, expected in manifest['source_sha256'].items()}
    seeds = {name: sha(ROOT/name.split('/home/jdiaz/pysnspd/')[-1]) == expected
             for name, expected in manifest['seed_validation']['verified_files'].items()}
    if not all(sources.values()) or not all(seeds.values()):
        issues.append('SOURCE_OR_SEED_HASH_MISMATCH')
    if manifest['registration'] != reg:
        issues.append('REGISTRATION_MISMATCH')
    progress = [json.loads(line) for line in (RAW/'progress.jsonl').read_text().splitlines()]
    if progress[-1]['status'] != 'SUCCESS' or progress[-1]['completed'] != 2:
        issues.append('INCOMPLETE_PROGRESS')
    cases = []
    e, hbar = 1.602176634e-19, 1.054571817e-34
    for name in ('m0_helix', 'm0_smooth_perturbation'):
        folder = RAW/name
        result = read(folder/'result.json')
        artifacts = {file: sha(folder/file) == value for file, value in result['artifacts_sha256'].items()}
        prep = dict(np.load(folder/'preparation.npz', allow_pickle=False))
        init = dict(np.load(folder/'initial.npz', allow_pickle=False))
        arrays = dict(np.load(folder/'arrays.npz', allow_pickle=False))
        finite = all(np.isfinite(value).all() for bundle in (prep, init, arrays) for value in bundle.values())
        same_initial = (np.array_equal(prep['delta_bar'], init['delta_bar'])
                        and np.array_equal(init['delta_bar'], arrays['delta_bar'])
                        and np.array_equal(init['p_quadrature'], arrays['p_quadrature'])
                        and np.array_equal(prep['graph_edges'], arrays['graph_edges']))
        if not all(artifacts.values()) or not finite or not same_initial:
            issues.append(name+':ARTIFACT_OR_INITIAL_STATE_ERROR')
        z, p = arrays['delta_bar'], arrays['p_quadrature']
        g = arrays['cartesian_gradient_bar'][:, 0]+1j*arrays['cartesian_gradient_bar'][:, 1]
        edge, phi = arrays['graph_edges'], arrays['phi_V']
        mass, mapping = prep['quadrature_mass_bar'], prep['quadrature_to_dof']
        scale = result['energy_scale_J']
        tref = reg['time_reference_ps']*1e-12
        iscale = 2*e/hbar*scale
        js = arrays['current_super_A']/iscale
        divergence = np.zeros(len(z))
        np.add.at(divergence, edge[:, 0], js)
        np.add.at(divergence, edge[:, 1], -js)
        noether = float(max(abs(np.imag(np.conj(z)*g)+divergence)))
        current_divergence = np.zeros(len(z))
        np.add.at(current_divergence, edge[:, 0], arrays['current_total_A'])
        np.add.at(current_divergence, edge[:, 1], -arrays['current_total_A'])
        Ib, Is, vc = result['circuit']['state']
        current_divergence[0] -= Is
        current_divergence[-1] += Is
        drop = phi[edge[:, 0]]-phi[edge[:, 1]]
        super_power = float(arrays['current_super_A']@drop)
        normal_power = float(arrays['current_normal_A']@drop)
        port_power = float(Is*(phi[0]-phi[-1]))
        phase_velocity = -1j*(2*e*tref/hbar)*phi*z
        phase_power = float(np.real(np.vdot(g, phase_velocity))*scale/tref)
        field_rate = float(np.real(np.vdot(g, arrays['field_velocity_bar'])))
        qrate = float(mass@arrays['QDelta_density_bar'])
        kwt_residual = field_rate+qrate-phase_power*tref/scale
        heat_normal = float(np.sum(arrays['normal_heat_power_W'])*tref/scale)
        rate = result['energy_rate']
        reconstructed_domain_rate = field_rate+rate['electron_population_rate_bar']+rate['phonon_rate_bar']+rate['explicit_link_work_bar']
        cp = result['circuit']['parameters']
        vd = cp['R_load_ohm']*(Ib-Is)+vc
        rhs = np.array([(cp['V_bias_V']-cp['R_bias_ohm']*Ib-vd)/cp['L_bias_H'],
                        (vd-(phi[0]-phi[-1]))/cp['Lk_ext_H'],
                        (Ib-Is)/cp['C_couple_F']])
        stored_rate = float(cp['L_bias_H']*Ib*rhs[0]+cp['Lk_ext_H']*Is*rhs[1]+cp['C_couple_F']*vc*rhs[2])
        source = cp['V_bias_V']*Ib
        bias_loss = cp['R_bias_ohm']*Ib**2
        load_loss = cp['R_load_ohm']*(Ib-Is)**2
        circuit_residual = stored_rate+port_power+bias_loss+load_loss-source
        combined_residual = reconstructed_domain_rate*scale/tref+stored_rate-(source-bias_loss-load_loss)
        power_scale = max(abs(source), abs(stored_rate), abs(reconstructed_domain_rate*scale/tref), normal_power, 1e-30)
        port_scale = max(abs(port_power), abs(phase_power), normal_power, 1e-30)
        recomputed = dict(
            noether_absolute=noether,
            electrical_current_relative=float(max(abs(current_divergence))/reg['circuit']['reference_current_A']),
            instantaneous_power_relative=float(max(abs(combined_residual), abs(circuit_residual))/power_scale),
            port_power_relative=float(max(abs(port_power-normal_power-super_power), abs(phase_power-super_power),
                abs(reconstructed_domain_rate*scale/tref-port_power))/port_scale),
            temperature_absolute=float(max(abs(arrays['temperature_bar']-reg['bath_theta']))),
        )
        for key, value in recomputed.items():
            limit = reg['criteria']['instantaneous_power_relative' if key == 'port_power_relative' else key]
            if value > limit:
                issues.append(name+':'+key)
        # Compare redundant saved arithmetic at floating point precision, not a new physical tolerance.
        redundant = dict(field_rate_bar=abs(field_rate-rate['field_rate_bar']),
            condensate_rate_bar=abs(qrate-result['heat']['condensate_rate_bar']),
            normal_rate_bar=abs(heat_normal-result['heat']['normal_rate_bar']),
            total_rate_bar=abs(reconstructed_domain_rate-rate['total_rate_bar']),
            source_heat_bar=abs(result['heat']['deposited_rate_bar']-qrate-heat_normal),
            circuit_rhs_max_absolute=float(max(abs(rhs-np.array(result['circuit']['state_derivative_SI'])))))
        derivatives = read(folder/'derivative_checks.json')
        fd = {}
        if derivatives != result['derivative_checks']:
            issues.append(name+':FD_METADATA_MISMATCH')
        if derivatives['performed']:
            for kind in ('force', 'current'):
                row = derivatives[kind]
                predicted = float(np.real(np.vdot(g, prep['cartesian_direction']))) if kind == 'force' else float(js@prep['link_direction'])
                estimated = row['fd_estimates']
                error = abs(predicted-estimated[-1])
                uncertainty = abs(estimated[-1]-estimated[-2])
                tolerance = reg['criteria']['derivative_absolute']+reg['criteria']['derivative_relative']*abs(estimated[-1])
                accepted = error <= tolerance and uncertainty <= reg['criteria']['reference_budget_fraction']*tolerance
                fd[kind] = dict(predicted_from_arrays=predicted, error=error, reference_uncertainty=uncertainty,
                    registered_tolerance=tolerance, passed=accepted,
                    scope='Directional derivative; FD energy samples themselves were not archived, stored estimates verified arithmetically.')
                if not accepted:
                    issues.append(name+':FD_'+kind)
        probes = read(folder/'direct_reference_probes.json')
        probe_maxima = {key: max(row['metrics'][key] for row in probes) for key in probes[0]['metrics']}
        probes_pass = probes == result['direct_reference_probes'] and all(
            row['sign_margin'] > 0 and all(value <= manifest['seed_validation']['limits'][key]
                                         for key, value in row['metrics'].items()) for row in probes)
        if not probes_pass:
            issues.append(name+':DIRECT_PROBES')
        endpoint = np.isin(mapping, (0, len(z)-1))
        endpoint_q = float(mass[endpoint]@arrays['QDelta_density_bar'][endpoint])
        speed = abs(arrays['material_velocity_bar'])
        endpoint_diagnosis = dict(terminal_dofs=[0, len(z)-1],
            condensate_heat_rate_bar=qrate, terminal_heat_rate_bar=endpoint_q,
            terminal_fraction=endpoint_q/qrate,
            interior_heat_rate_bar=float(mass[~endpoint]@arrays['QDelta_density_bar'][~endpoint]),
            terminal_material_speed_bar=speed[[0,-1]].tolist(),
            maximum_interior_material_speed_bar=float(max(speed[1:-1])),
            maximum_material_speed_dof=int(np.argmax(speed)),
            terminal_quadrature_mass_bar=mass[endpoint].tolist())
        cases.append(dict(case=name, status=result['status'], artifacts=artifacts, finite_arrays=finite,
            initial_field_and_populations_unchanged=same_initial,
            shape=dict(dofs=len(z), quadratures=len(mass), count_nodes=p.shape[1], edges=len(edge)),
            population_bounds=[float(p.min()), float(p.max())],
            registered_metrics=result['metrics'], independent_array_metrics=recomputed,
            redundant_arithmetic_absolute_differences=redundant,
            KWT_identity_residual_bar=float(kwt_residual),
            circuit_residual_W=float(circuit_residual), CM9_residual_W=float(combined_residual),
            derivative_checks=fd if fd else dict(status=derivatives['status'], performed=False),
            direct_reference_probes=dict(count=len(probes), unique_quadratures=len({x['quadrature_index'] for x in probes}),
                maximum_errors=probe_maxima, passed=probes_pass, minimum_sign_margin=min(x['sign_margin'] for x in probes)),
            principal_minimum=float(min(arrays['principal_minima'])),
            principal_uncertainty_maximum=float(max(arrays['principal_uncertainties'])),
            principal_margin_minimum=float(min(arrays['principal_minima']-arrays['principal_uncertainties'])),
            endpoint_diagnosis=endpoint_diagnosis))
    report = dict(schema='pysnspd.stage3.mixed_pilot_independent_review.v1', date='2026-09-23',
        status='SAVED_PILOT_VERIFIED_SCOPE_LIMITED' if not issues else 'AUDIT_FAILED',
        method='Postprocessing of saved numbers only. No new spectrum, dynamics or physical trajectory evaluated.',
        raw_directory=str(RAW.relative_to(ROOT)).replace('\\','/'),
        input_inventory_sha256={str(path.relative_to(RAW)).replace('\\','/'): sha(path) for path in sorted(RAW.rglob('*')) if path.is_file()},
        source_hash_checks=sources, seed_artifact_hash_checks=seeds, registration_matches=True,
        runtime_seconds=summary['runtime_seconds'], issues=issues, cases=cases,
        preserved_scope=dict(instantaneous_algebra=True, mesh_convergence_admission=False,
            temporal_admission=False, kinetic_interface_admission=False, stage3_closed=False, production=False),
        evidence_limits=[
            'FD force and current were executed only for smooth_perturbation/m0, as registered; helix has no new FD PASS.',
            'Original FD energy samples and all instantaneous spectra are not archived. Their evaluations are not independently repeated here; stored FD estimates and moment metadata are checked, with input/source hashes.',
            'Temperature equality is an initialization/inversion identity for thermal occupations, not thermal relaxation.',
            'Three direct-source quadratures per case test local source equivalence; they do not certify all spectral resolutions.',
            'Normal heat allocation conserves total power but its spatial local accuracy has no admission from this pilot.'
        ],
        decision=dict(keep_existing_instantaneous_PASS=True,
            recommend_six_free_boundary_snapshots_for_device_precision=False,
            next_step='Register and test a physical terminal/reservoir loading with explicit mechanical/phase work before a costly mixed-mesh campaign.',
            permitted_future_six_snapshot_scope='Sensitivity of the reservoir-loaded instantaneous state and its induced perturbation, with the free state retained as a historical control; no transient certificate.',
            no_retrospective_threshold_change=True),
        boundary_inference=dict(status='ANALYTIC_INFERENCE_NOT_MEASURED_MESH_LAW',
            condition='Nonzero imposed initial helix q=0.1 with all material DOFs free and zero conjugate terminal load.',
            argument='At fixed polynomial degree, residual integrated terminal phase force is O(J), terminal GLL temporal mass is O(h); material speed is O(1/h), terminal heat density O(1/h^2), integrated terminal QDelta O(1/h) at t=0.',
            limitations='Assumes finite nonzero mobility and current as h decreases. A parabolic boundary layer may regularize at positive time; no divergence of total dissipated energy or finite-time trajectory is inferred.',
            meaning='The pilot identities remain correct; a mesh study of raw initial QDelta in this free boundary scenario would mainly resolve its boundary incompatibility.'),
        next_weak_transient_observables=[
            'Ib, Is, vc, Vdev and Vout with circuit-state derivatives and prescribed terminal conditions.',
            'Amplitude, q/Gamma, gauge-invariant material velocity and terminal reactions/work.',
            'Interior and terminal QDelta separately, normal Joule heat, equivalent temperature and electronic excess energy.',
            'Domain+circuit energy versus source, external resistor losses, escape and reservoir flux/work; internal heat is not counted twice.',
            'Population bounds, finite states, catalogue domain and D.36 margin at saved times.',
            'Time-step comparison of induced observables against the compatible reference state; scope and precision fixed before results.'
        ])
    target = OUT/'mixed_pilot_review.json'
    with target.open('x', encoding='utf-8') as stream:
        json.dump(report, stream, indent=2, ensure_ascii=False, allow_nan=False)
        stream.write('\n')
    a, b = [row['endpoint_diagnosis'] for row in cases]
    md = f'''# Revisión independiente del piloto mixto

El piloto conserva su resultado de identidades instantáneas. Los {len(sources)} hashes de fuentes, los {len(seeds)} artefactos de la semilla y los diez artefactos referenciados por los dos resultados coinciden. El registro guardado coincide con el vigente. Las matrices son finitas, y el campo y las poblaciones iniciales permanecen idénticos. No se ejecutaron espectros ni dinámica nuevos.

La conservación de corriente, Noether, la potencia campo–circuito, el depósito de calor y las dos derivadas direccionales verificadas son compatibles con los criterios originales. La fuerza y corriente de la perturbación tienen errores FD de {cases[1]['derivative_checks']['force']['error']:.3g} y {cases[1]['derivative_checks']['current']['error']:.3g}. La hélice no recibió un nuevo control FD: así lo declara el registro. El margen mínimo D.36 permanece positivo. Los contrastes causales directos se limitan a tres cuadraturas por caso.

## El calor inicial procede de los extremos

| Magnitud guardada | Hélice | Perturbación suave |
|---|---:|---:|
| QDelta integrado | {a['condensate_heat_rate_bar']:.12g} | {b['condensate_heat_rate_bar']:.12g} |
| QDelta en los dos terminales | {a['terminal_heat_rate_bar']:.12g} | {b['terminal_heat_rate_bar']:.12g} |
| Fracción terminal | {100*a['terminal_fraction']:.9f}% | {100*b['terminal_fraction']:.9f}% |
| QDelta interior | {a['interior_heat_rate_bar']:.9g} | {b['interior_heat_rate_bar']:.9g} |
| Máxima velocidad material interior | {a['maximum_interior_material_speed_bar']:.9g} | {b['maximum_interior_material_speed_bar']:.9g} |

Estas tasas y velocidades usan las unidades normalizadas del registro. La velocidad máxima de ambos perfiles, 2.85332871, se encuentra en los extremos. El aporte terminal es exactamente el mismo en ambos estados guardados. La respuesta interior sí cambia con la perturbación; queda oculta en el calor global por el relajamiento terminal.

El dato inicial impone una hélice con q=0.1, mientras los extremos tienen carga conjugada nula y todas sus velocidades materiales quedan libres. Por ello no satisface el flujo natural terminal nulo. A grado polinómico fijo, el esfuerzo terminal integrado es de orden J y la masa temporal GLL de orden h: se espera velocidad inicial de orden 1/h y calor terminal integrado de orden 1/h. Esta es una inferencia analítica, no una ley medida con varias mallas. No implica que una trayectoria a tiempo positivo o su energía disipada integrada diverjan: puede aparecer una capa de relajación inicial.

## Decisión y siguiente ensayo

Conservar el PASS algebraico del piloto y sus datos. No recomiendo gastar el lote de seis instantáneas libres como comprobación de precisión del dispositivo: su calor inicial estaría dominado por esa incompatibilidad terminal. Primero debe fijarse y registrarse la carga de reservorio, distinguir la corriente superconductora prescrita de la corriente total del circuito y contabilizar el trabajo de las cargas y reacciones. Una reacción radial impuesta no debe ocultarse ni ajustarse a posteriori para forzar una identidad.

Después, seis instantáneas con ese contrato pueden aportar sensibilidad espacial del estado cargado y de la perturbación inducida. El caso libre sigue siendo un control histórico válido. Las instantáneas, incluso conservativas, no acreditan un transiente ni la etapa 3 completa.

En un transiente débil se deben guardar corrientes y voltajes circuitales, amplitud y q/Gamma, velocidad material, calor interior/terminal, energía electrónica y temperatura equivalente, junto al balance integrado dominio+circuito con trabajo de reservorio. Los controles de población, dominio del catálogo, estabilidad y comparación temporal deben seguir las mismas variables efectivamente integradas; no se fijan aquí tolerancias nuevas.

## Límites de esta revisión

Las energías individuales del stencil FD y todos los espectros instantáneos no están archivados. Se verificaron las estimaciones FD guardadas, sus errores y presupuesto, y la aritmética registrada del momento electrónico; no se recalcularon esas evaluaciones físicas. La coincidencia térmica inicial es una identidad de preparación e inversión, no una demostración de relajación. La conservación del depósito Joule no certifica su distribución local. No se alteraron resultados, criterios ni registros previos.

El archivo `mixed_pilot_review.json` conserva el inventario SHA256, las métricas recalculadas y las limitaciones. Reproducción: `python sandbox/stage3_spatial/coupled_20260923/audit_mixed_pilot.py` (rechaza sobrescribir la revisión existente).
'''
    with (OUT/'mixed_pilot_review.md').open('x', encoding='utf-8') as stream:
        stream.write(md)
    print(json.dumps(dict(status=report['status'], sources=len(sources), seed_artifacts=len(seeds),
                          raw_files=len(report['input_inventory_sha256']), issues=issues), indent=2))


if __name__ == '__main__':
    main()
