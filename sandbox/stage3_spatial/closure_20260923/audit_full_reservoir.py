"""Full six-case offline audit. No spectral query, dynamical RHS or time integration."""
from pathlib import Path
import hashlib
import json
import numpy as np

ROOT = Path(__file__).resolve().parents[3]
BASE = ROOT/'docs/implementation/stage3/closure_20260923'
COUPLED = ROOT/'docs/implementation/stage3/coupled_20260923'
RAW = BASE/'raw/stage3_reservoir_full_20260923'
DATA = RAW/'reservoir'
INPUT = RAW/'mixed_control'


def read(path):
    return json.loads(path.read_text(encoding='utf-8'))


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    summary, manifest = read(DATA/'summary.json'), read(DATA/'manifest.json')
    registration = read(COUPLED/'reservoir_postprocess_registration.json')
    old_registration = read(COUPLED/'reservoir_postprocess_registration_before_preparation_recovery.json')
    sources = {p: sha(ROOT/p) == digest for p, digest in manifest['source_sha256'].items()}
    inputs = {p: sha(INPUT/p) == digest for p, digest in manifest['input_artifacts_sha256'].items()}
    issues = []
    if not all(sources.values()) or not all(inputs.values()):
        issues.append('HASH_MISMATCH')
    if sha(DATA/'manifest.json') != summary['manifest_sha256'] or manifest['registration'] != registration:
        issues.append('MANIFEST_OR_REGISTRATION_MISMATCH')
    if registration['criteria'] != old_registration['criteria']:
        issues.append('PHYSICAL_CRITERIA_CHANGED')
    cases = []
    e, hbar = 1.602176634e-19, 1.054571817e-34
    for item in summary['cases']:
        name = item['case_id']
        record = read(DATA/name/'result.json')
        artifacts = {p: sha(DATA/name/p) == digest for p, digest in record['artifacts_sha256'].items()}
        if record != item or not all(artifacts.values()) or sha(INPUT/name/'result.json') != record['input_record_sha256']:
            issues.append(name+':RECORD_OR_ARTIFACT_MISMATCH')
        old = read(INPUT/name/'result.json')
        a = dict(np.load(DATA/name/'arrays.npz', allow_pickle=False))
        before = dict(np.load(INPUT/name/'arrays.npz', allow_pickle=False))
        initial = dict(np.load(INPUT/name/'initial.npz', allow_pickle=False))
        prep = dict(np.load(INPUT/name/'preparation.npz', allow_pickle=False))
        z, p, E = a['delta_bar'], a['p_quadrature'], a['recovered_saved_energies_bar']
        mass, mapping = a['quadrature_mass_bar'], a['quadrature_to_dof']
        weights = initial['count_weights']
        unchanged = all(np.array_equal(a[k], before[k]) for k in ('delta_bar', 'p_quadrature'))
        mass_exact = np.array_equal(mass, prep['quadrature_mass_bar'])
        mapping_exact = np.array_equal(mapping, prep['quadrature_to_dof'])
        finite = all(np.isfinite(v).all() for v in a.values())
        if not (unchanged and mass_exact and mapping_exact and finite):
            issues.append(name+':ARCHIVED_STATE_OR_MASS_CHANGED')
        recovered = .12*(np.log1p(-p)-np.log(p))
        E_error = float(max(abs(E-recovered).flat))
        E_relative_roundoff = float(np.max(abs(E-recovered)/E))
        ordering = bool(np.all(E > 0) and np.all(np.diff(E, axis=1) > 0))
        old_moment = float(mass@(4*np.sum(weights*E*before['electron_heating_rhs'], axis=1)))
        new_moment = float(mass@(4*np.sum(weights*E*a['electron_heating_rhs'], axis=1)))
        scale, tref = record['energy_scale_J'], record['time_scale_s']
        iscale = 2*e/hbar*scale
        g = before['cartesian_gradient_bar'][:, 0]+1j*before['cartesian_gradient_bar'][:, 1]
        load = a['boundary_load_cartesian_bar'][:, 0]+1j*a['boundary_load_cartesian_bar'][:, 1]
        radial = a['radial_reaction_cartesian_bar'][:, 0]+1j*a['radial_reaction_cartesian_bar'][:, 1]
        material, velocity = a['material_velocity_bar'], a['field_velocity_bar']
        terminal = np.array(record['terminal']['nodes'])
        reference = registration['reservoir']['reference_current_A']
        expected_phase_load = np.array([-reference, reference])/iscale
        phase_load_error = float(max(abs(np.imag(np.conj(z[terminal])*load[terminal])-expected_phase_load)))
        radius_rates = np.real(np.conj(z[terminal])/abs(z[terminal])*material[terminal])
        load_work = float(np.real(np.vdot(load, material)))
        radial_work = float(np.real(np.vdot(radial, material)))
        field_work = float(np.real(np.vdot(g, velocity)))
        phase_velocity = -1j*(2*e*tref/hbar)*a['phi_V']*z
        gauge_work = float(np.real(np.vdot(g, phase_velocity)))
        q = float(mass@a['QDelta_density_bar'])
        qn = float(np.sum(a['normal_heat_power_W'])*tref/scale)
        field_heat_residual = field_work+q-load_work-gauge_work
        heat_residual = new_moment-q-qn
        edge = before['graph_edges']
        drop = a['phi_V'][edge[:, 0]]-a['phi_V'][edge[:, 1]]
        current_res = np.zeros(len(z))
        np.add.at(current_res, edge[:, 0], a['current_total_A'])
        np.add.at(current_res, edge[:, 1], -a['current_total_A'])
        current_res[terminal] -= [reference, -reference]
        cp, state = record['circuit']['parameters'], record['circuit']['state']
        Ib, Is, vc = state
        expected_state = [reference, reference, 0.]
        if state != expected_state:
            issues.append(name+':CIRCUIT_STATE_NOT_REGISTERED')
        Vdev = float(a['phi_V'][terminal[0]]-a['phi_V'][terminal[1]])
        vd = cp['R_load_ohm']*(Ib-Is)+vc
        rhs = np.array([(cp['V_bias_V']-cp['R_bias_ohm']*Ib-vd)/cp['L_bias_H'],
                        (vd-Vdev)/cp['Lk_ext_H'], (Ib-Is)/cp['C_couple_F']])
        stored_rate = float(cp['L_bias_H']*Ib*rhs[0]+cp['Lk_ext_H']*Is*rhs[1]+cp['C_couple_F']*vc*rhs[2])
        external = cp['V_bias_V']*Ib-cp['R_bias_ohm']*Ib**2-cp['R_load_ohm']*(Ib-Is)**2+load_work*scale/tref
        domain = (field_work+new_moment)*scale/tref
        cm9 = domain+stored_rate-external
        port = Is*Vdev
        power_norm = max(abs(domain), abs(port), abs(load_work*scale/tref), record['potential']['normal_joule_W'], 1e-30)
        recomputed = dict(original_heat_moment_relative=abs(old_moment-old['heat']['deposited_rate_bar'])/abs(old['heat']['deposited_rate_bar']),
            terminal_radius_rate_absolute=float(max(abs(radius_rates))),
            radial_reaction_work_relative=abs(radial_work)/max(abs(load_work-radial_work), q, 1e-30),
            instantaneous_power_relative=abs(cm9)/power_norm,
            heat_moment_relative=abs(heat_residual)/max(abs(new_moment), q, qn, 1e-30),
            current_continuity_relative=float(max(abs(current_res))/reference),
            KWT_identity_relative=abs(field_heat_residual)/max(abs(field_work), q, abs(load_work), abs(gauge_work), 1e-30))
        registered_pass = all(np.isfinite(value) and value <= registration['criteria'][k] for k, value in record['metrics'].items())
        independent_pass = all(np.isfinite(value) and value <= registration['criteria'][k] for k, value in recomputed.items())
        if not registered_pass or not independent_pass or not ordering or E_relative_roundoff > 64*np.finfo(float).eps:
            issues.append(name+':ARITHMETIC_CHECK_FAILED')
        cases.append(dict(case=name, artifact_hash_checks=artifacts, unchanged_field_and_populations=unchanged,
            exact_archived_mass=mass_exact, exact_mapping=mapping_exact, finite_arrays=finite,
            decoded_energies_positive_and_ordered=ordering, decoded_energies_formula_max_error=E_error, decoded_energies_formula_relative_roundoff=E_relative_roundoff,
            registered_metrics=record['metrics'], independently_recomputed_metrics=recomputed,
            phase_load_formula_max_error=phase_load_error, circuit_rhs_max_error=float(max(abs(rhs-record['circuit']['state_derivative_SI']))),
            CM9_recomputed_residual_W=cm9, KWT_recomputed_residual_bar=field_heat_residual,
            geometry_mass_relative_roundoff=record['recovery']['geometry_mass_relative_roundoff'],
            condensate_heat_free_bar=old['heat']['condensate_rate_bar'], condensate_heat_loaded_bar=q,
            condensate_heat_terminal_bar=record['reservoir']['terminal_condensate_rate_bar'],
            maximum_material_speed_bar=float(max(abs(material))),
            reservoir_work_power_W=load_work*scale/tref, radial_reaction_work_rate_bar=radial_work,
            Vdev_V=Vdev, adjacent_normal_current_over_Iref=record['terminal']['adjacent_internal_normal_current_over_Iref'],
            scope='Snapshot with fixed terminal radius and Iref load; not a new spectrum or trajectory.'))
    # Audit the generation phase independently of its postprocess receipt.
    mixed_summary = read(INPUT/'summary.json')
    mixed_manifest = read(INPUT/'manifest.json')
    batch_summary, batch_manifest = read(RAW/'summary.json'), read(RAW/'manifest.json')
    mixed_reg = mixed_manifest['registration']
    expected_cases = {f'{mesh}_{profile}' for mesh in ('m0', 'm1', 'm2')
                      for profile in ('helix', 'smooth_perturbation')}
    by_name = {row['case']: row for row in cases}
    mixed_records = {row['mesh']['id']+'_'+row['profile']['id']: row for row in mixed_summary['cases']}
    if set(by_name) != expected_cases or len(cases) != 6 or set(mixed_records) != expected_cases:
        issues.append('SIX_CASE_COVERAGE_MISMATCH')
    source_groups = {'postprocess': sources}
    for group, values in (('mixed', mixed_manifest['source_sha256']), ('batch', batch_manifest['source_sha256'])):
        source_groups[group] = {path: sha(ROOT/path) == digest for path, digest in values.items()}
        if not all(source_groups[group].values()):
            issues.append(group+':SOURCE_HASH_MISMATCH')
    seed_checks = {path: sha(ROOT/path.split('/home/jdiaz/pysnspd/')[-1]) == digest
                   for path, digest in mixed_manifest['seed_validation']['verified_files'].items()}
    if not all(seed_checks.values()):
        issues.append('SEED_ARTIFACT_HASH_MISMATCH')
    if mixed_manifest['source_sha256'] != mixed_summary['source_sha256']:
        issues.append('MIXED_SOURCE_MANIFEST_DISAGREEMENT')
    if mixed_reg != read(COUPLED/'mixed_registration.json'):
        issues.append('MIXED_REGISTRATION_CHANGED')
    if batch_summary['status'] != 'RESERVOIR_INSTANTANEOUS_BATCH_PASS_SCOPE_LIMITED':
        issues.append('BATCH_NOT_COMPLETE')
    direct_maxima = {}
    fd_checks = {}
    for name, m in mixed_records.items():
        folder = INPUT/name
        if read(folder/'result.json') != m:
            issues.append(name+':MIXED_SUMMARY_MISMATCH')
        if m['status'] != 'PASS_MIXED_INSTANTANEOUS_IDENTITIES':
            issues.append(name+':MIXED_CASE_NOT_PASS')
        for key, value in m['metrics'].items():
            threshold_key = 'instantaneous_power_relative' if key == 'port_power_relative' else key
            if not np.isfinite(value) or value > mixed_reg['criteria'][threshold_key]:
                issues.append(name+':MIXED_'+key)
        for artifact, digest in m['artifacts_sha256'].items():
            if sha(folder/artifact) != digest:
                issues.append(name+':MIXED_ARTIFACT_'+artifact)
        a = dict(np.load(folder/'arrays.npz', allow_pickle=False))
        prep = dict(np.load(folder/'preparation.npz', allow_pickle=False))
        g = a['cartesian_gradient_bar'][:, 0]+1j*a['cartesian_gradient_bar'][:, 1]
        js = a['current_super_A']/(2*e/hbar*m['energy_scale_J'])
        divergence = np.zeros(len(a['delta_bar']))
        np.add.at(divergence, a['graph_edges'][:, 0], js)
        np.add.at(divergence, a['graph_edges'][:, 1], -js)
        noether = float(max(abs(np.imag(np.conj(a['delta_bar'])*g)+divergence)))
        if noether > mixed_reg['criteria']['noether_absolute']:
            issues.append(name+':NOETHER_ARRAY_CHECK')
        principal_margin = float(min(a['principal_minima']-a['principal_uncertainties']))
        if principal_margin <= 0:
            issues.append(name+':PRINCIPAL_SIGN_NOT_ADMITTED')
        derivative = m['derivative_checks']
        if read(folder/'derivative_checks.json') != derivative:
            issues.append(name+':DERIVATIVE_FILE_DISAGREEMENT')
        expected_fd = name == 'm0_smooth_perturbation'
        if derivative['performed'] != expected_fd:
            issues.append(name+':FD_SCOPE_MISMATCH')
        fd_checks[name] = {'performed': derivative['performed'], 'status': derivative['status']}
        if derivative['performed']:
            for kind in ('force', 'current'):
                d = derivative[kind]
                predicted = float(np.real(np.vdot(g, prep['cartesian_direction']))) if kind == 'force' else float(js@prep['link_direction'])
                error = abs(predicted-d['fd_estimates'][-1])
                uncertainty = abs(d['fd_estimates'][-1]-d['fd_estimates'][-2])
                tolerance = mixed_reg['criteria']['derivative_absolute']+mixed_reg['criteria']['derivative_relative']*abs(d['fd_estimates'][-1])
                ok = error <= tolerance and uncertainty <= mixed_reg['criteria']['reference_budget_fraction']*tolerance
                fd_checks[name][kind] = dict(error=error, reference_uncertainty=uncertainty, tolerance=tolerance, passed=ok)
                if not ok:
                    issues.append(name+':FD_'+kind)
        probes = read(folder/'direct_reference_probes.json')
        if probes != m['direct_reference_probes'] or len(probes) != 3:
            issues.append(name+':DIRECT_PROBE_COVERAGE')
        for p in probes:
            for key, value in p['metrics'].items():
                direct_maxima[key] = max(direct_maxima.get(key, 0.), value)
                if value > mixed_manifest['seed_validation']['limits'][key] or p['sign_margin'] <= 0:
                    issues.append(name+':DIRECT_SOURCE_PROBE')
        c = by_name[name]
        c.update(mesh=m['mesh']['id'], profile=m['profile']['id'], dofs=m['dofs'], quadratures=m['quadrature_points'],
                 energy_bar=m['internal_energy_bar'], free_energy_bar=m['free_energy_bar'],
                 gradient_remainder_bar=m['gradient_remainder_bar'], principal_minimum=float(min(a['principal_minima'])),
                 principal_uncertainty_maximum=float(max(a['principal_uncertainties'])), principal_margin_minimum=principal_margin,
                 noether_recomputed_absolute=noether,
                 amplitude_range=[float(min(a['amplitude_quadrature_bar'])), float(max(a['amplitude_quadrature_bar']))],
                 gamma_range=[float(min(a['gamma_quadrature_bar'])), float(max(a['gamma_quadrature_bar']))],
                 adjacent_normal_current_max_abs_over_Iref=max(abs(x) for x in c['adjacent_normal_current_over_Iref']))

    def sensitivity(values):
        values = list(map(float, values))
        d01, d12 = abs(values[1]-values[0]), abs(values[2]-values[1])
        return dict(values=values, absolute_coarse_to_medium=d01, absolute_medium_to_fine=d12,
                    medium_to_fine_relative_to_fine=None if values[2] == 0 else d12/abs(values[2]),
                    successive_difference_decreases=d12 < d01,
                    observed_order=None if d01 == 0 or d12 == 0 else float(np.log2(d01/d12)),
                    scope='Descriptive only; no registered mesh-accuracy threshold or convergence certificate.')

    observables = ('energy_bar', 'condensate_heat_loaded_bar', 'condensate_heat_terminal_bar',
                   'maximum_material_speed_bar', 'Vdev_V', 'adjacent_normal_current_max_abs_over_Iref',
                   'principal_minimum', 'principal_uncertainty_maximum', 'reservoir_work_power_W')
    comparison = {profile: {key: sensitivity([by_name[f'm{i}_{profile}'][key] for i in range(3)])
                            for key in observables} for profile in ('helix', 'smooth_perturbation')}
    excess = {key: sensitivity([by_name[f'm{i}_smooth_perturbation'][key]-by_name[f'm{i}_helix'][key]
                               for i in range(3)])
              for key in ('energy_bar', 'condensate_heat_loaded_bar', 'Vdev_V')}
    supported = not issues
    report = dict(schema='pysnspd.stage3.full_reservoir_independent_review.v1', date='2026-09-23',
        status='VERIFIED_SIX_SNAPSHOT_CAMPAIGN_SCOPE_LIMITED' if supported else 'AUDIT_FAILED', issues=issues,
        method='Offline arithmetic on saved arrays only; no catalogue query, PDE RHS, time integration or new material calculation.',
        source_hash_checks=source_groups, seed_artifact_hash_checks=seed_checks, input_hash_checks=inputs,
        raw_inventory_sha256={p.relative_to(RAW).as_posix(): sha(p) for p in sorted(RAW.rglob('*')) if p.is_file()},
        runtime_seconds=dict(total=batch_summary['runtime_seconds'], material=mixed_summary['runtime_seconds'], postprocess=summary['runtime_seconds']),
        cases=cases, finite_difference_scope=fd_checks, direct_source_probe_maximum_errors=direct_maxima,
        comparison=comparison, induced_perturbation_minus_helix=excess,
        mesh_scope='Three meshes of the same geometry/profile. The registration requested sensitivity only; no quantitative accuracy certificate is added after results.',
        decision=dict(development_closure_supported=supported,
            authority='The parent reports explicit user choice: cerrar desarrollo e iniciar3.5. This review supports only that bounded development closure.',
            closed_development_scope='Common spatial energy, forces/current, mixed2D+1D uniform-trace interface, potential and instantaneous KWT/heat/three-state-circuit/reservoir work in the registered weak states.',
            stage3_complete_physical_dynamic_admission=False, temporal_admission=False, full_D27_admission=False,
            kinetic_interface_admission=False, production=False,
            proceed_to_stage3_5_research=True,
            stage3_5_scope='Research physical/numerical operating ranges, margin definitions and required future evidence; no simulation sweep implied.'),
        pending_before_claiming_complete_transients=[
            {'stage': '3B', 'requirement': 'Dynamic reservoir amplitude a_b(Is), derivative da_b/dIs and consistent reservoir work/energy exchange.'},
            {'stage': '3B', 'requirement': 'Kinetic f(E) matching and conservative energy flux through the2D–1D interface; the present map admits the uniform field trace only.'},
            {'stage': '3B', 'requirement': 'External normal-boundary trace and continuation-length influence in the intended dynamic scenario; an adjacent internal face is not that certificate.'},
            {'stage': '3D/3E', 'requirement': 'A weak coupled time trajectory with spatial field, occupations and the memory circuit, assessed at common physical times with integrated CM9.'},
            {'stage': '3E', 'requirement': 'Prospectively registered localized deposition and applicable time/space/spectral accuracy for its observables.'}
        ],
        interpretation_limits=[
            'No six-case FD claim: new finite-difference force/current controls exist only for m0/smooth_perturbation.',
            'Three direct-source samples per case and D.36 at saved quadratures are sampled evidence, not a global material/range certificate.',
            'FD energy decoding is valid only for the saved thermal preparation. Fermi–Dirac inversion is not new microscopic validation.',
            'CM9 is instantaneous; no energy history, current pulse, latency or detector threshold was simulated.',
            'Small pairwise changes describe mesh sensitivity; near-zero baseline observables need absolute reporting, not relative claims.',
            'The raw stage3_closed=false values are preserved; a separate user-authorized development decision does not rewrite them.',
            'The initial entry_contract.json is historical preparation. Any current status should identify the later records rather than reinterpret its unstarted flags.'
        ])
    BASE.mkdir(parents=True, exist_ok=True)
    with (BASE/'full_reservoir_review.json').open('x', encoding='utf-8') as stream:
        json.dump(report, stream, ensure_ascii=False, indent=2, allow_nan=False)
        stream.write('\n')
    lines = ['# Auditoría independiente del lote mixto completo', '',
        f"Dictamen: **{report['status']}**. El lote terminó en {batch_summary['runtime_seconds']:.3f} s; contiene seis estados instantáneos, dos perfiles sobre tres mallas. La revisión usa sólo archivos guardados y no ejecuta física nueva.", '',
        '## Integridad e identidades', '',
        f"Se cotejaron las fuentes del corredor ({len(source_groups['batch'])}), del cálculo espacial ({len(source_groups['mixed'])}) y del postproceso ({len(sources)}), cinco artefactos de la semilla y {len(inputs)} entradas del postproceso. Los casos coinciden con sus resúmenes; el inventario contiene {len(report['raw_inventory_sha256'])} archivos. Incidencias: {len(issues)}.", '',
        'Campo, poblaciones y masas permanecen idénticos entre la instantánea original y su lectura con reservorio. Se recompusieron momentos energéticos, carga de fase, trabajo de borde, continuidad, reacción radial, ecuaciones circuitales y CM9. Los criterios originales se mantienen. El control de diferencias finitas se ejecutó sólo en m0 perturbada, como estaba registrado; no se atribuye a los otros cinco casos.', '',
        '## Comparación de las mallas', '',
        '| Perfil/malla | Energía [bar] | QDelta [bar] | Velocidad máxima [bar] | Vdev [V] | max abs(Jn)/Iref adyacente |',
        '|---|---:|---:|---:|---:|---:|']
    for c in cases:
        lines.append(f"| {c['case']} | {c['energy_bar']:.12g} | {c['condensate_heat_loaded_bar']:.9g} | {c['maximum_material_speed_bar']:.9g} | {c['Vdev_V']:.9g} | {c['adjacent_normal_current_max_abs_over_Iref']:.9g} |")
    lines += ['',
        'Las etiquetas [bar] indican magnitudes normalizadas del modelo; no representan la unidad de presión. La hélice cargada reduce fuertemente su respuesta residual al refinar. La perturbación retiene una respuesta interior finita. Para evitar que el error de la hélice o una gran energía de fondo oculten esa respuesta, se compara además perturbación menos hélice:', '',
        '| Exceso inducido | m0 | m1 | m2 | Cambio m1→m2 / valor m2 |', '|---|---:|---:|---:|---:|']
    for key, s in excess.items():
        lines.append(f"| {key} | {s['values'][0]:.10g} | {s['values'][1]:.10g} | {s['values'][2]:.10g} | {100*s['medium_to_fine_relative_to_fine']:.6g}% |")
    worst_cm9 = max(abs(c['CM9_recomputed_residual_W']) for c in cases)
    margin = min(c['principal_margin_minimum'] for c in cases)
    lines += ['',
        f"El menor margen D.36 observado es {margin:.9g}; el residuo máximo de CM9 recompuesta es {worst_cm9:.4g} W. El trabajo de reservorio queda incluido con su signo. El JSON conserva los valores de trabajo, calor terminal, incertidumbre de D.36 y las diferencias por malla.", '',
        'Estos cambios son sensibilidad observada. El registro no fijó un umbral de precisión multimalla y aquí no se inventa uno. Una disminución del defecto de corriente en la cara interior próxima al extremo tampoco acredita por sí sola la traza normal externa de D.27. Las energías de fondo se informan, pero no se usan para ocultar el tamaño relativo del exceso inducido.', '',
        '## Cierre de desarrollo y pendientes', '',
        ('La evidencia respalda el cierre de desarrollo acotado que el usuario autorizó y el inicio de 3.5 como investigación de rangos y márgenes. ' if supported else 'La auditoría tiene incidencias que deben resolverse antes de un cierre apoyado en estos datos. ')+
        'Ese cierre conserva los resultados y sus límites. No equivale a completar una validación temporal, toda D.27 o la interfaz cinética; los campos stage3_closed=false de los datos originales permanecen intactos.', '',
        'Antes de atribuir resultados a transientes completos siguen pendientes el radio de reservorio dependiente de Is y su derivada, el intercambio cinético y energético en la unión 2D–1D, la condición normal externa apropiada y una trayectoria débil acoplada con balance integrado y comparación temporal. El depósito localizado necesita su propio ensayo. Estas obligaciones se registran para las etapas dependientes; no se ejecutan barridos ni se bloquea la investigación 3.5 por no disponer aún de ellos.', '',
        'El contrato de entrada histórico describía la preparación inicial y contiene flags de etapa no iniciada. No debe citarse como estado actual ignorando la evidencia posterior. La decisión actual debe referenciar esta auditoría y el cierre de desarrollo separado.', '',
        'Reproducción: `python sandbox/stage3_spatial/closure_20260923/audit_full_reservoir.py`. Sólo aritmética e integridad sobre archivos guardados; rechaza sobrescribir esta revisión.', '']
    with (BASE/'full_reservoir_review.md').open('x', encoding='utf-8') as stream:
        stream.write('\n'.join(lines))
    print(json.dumps(dict(status=report['status'], issues=issues, total_cases=len(cases),
                          raw_files=len(report['raw_inventory_sha256']), maximum_CM9_residual_W=worst_cm9,
                          minimum_D36_margin=margin, induced_comparisons=excess), indent=2))


if __name__ == '__main__':
    main()
