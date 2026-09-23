"""Independent archived-array arithmetic for the reservoir snapshot postprocess."""
from pathlib import Path
import hashlib
import json
import numpy as np

ROOT = Path(__file__).resolve().parents[3]
BASE = ROOT/'docs/implementation/stage3/coupled_20260923'
DATA = BASE/'reservoir_postprocess_v2'
INPUT = BASE/'raw/mixed_pilot'


def read(path):
    return json.loads(path.read_text(encoding='utf-8'))


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    summary, manifest = read(DATA/'summary.json'), read(DATA/'manifest.json')
    registration = read(BASE/'reservoir_postprocess_registration.json')
    old_registration = read(BASE/'reservoir_postprocess_registration_before_preparation_recovery.json')
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
        if not registered_pass or not independent_pass or not ordering or E_error != 0:
            issues.append(name+':ARITHMETIC_CHECK_FAILED')
        cases.append(dict(case=name, artifact_hash_checks=artifacts, unchanged_field_and_populations=unchanged,
            exact_archived_mass=mass_exact, exact_mapping=mapping_exact, finite_arrays=finite,
            decoded_energies_positive_and_ordered=ordering, decoded_energies_formula_max_error=E_error,
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
    review = dict(schema='pysnspd.stage3.reservoir_postprocess_review.v1', date='2026-09-23',
        status='VERIFIED_RESERVOIR_SNAPSHOT_POSTPROCESS' if not issues else 'AUDIT_FAILED', issues=issues,
        source_hash_checks=sources, input_artifact_hash_checks=inputs,
        output_inventory_sha256={p.relative_to(DATA).as_posix(): sha(p) for p in sorted(DATA.rglob('*')) if p.is_file()},
        criteria_unchanged_since_failed_preparation=registration['criteria'] == old_registration['criteria'],
        first_attempt_preserved=dict(failure_sha256=sha(BASE/'reservoir_postprocess/failure.json'),
            reason='Bitwise cross-platform GLL mass comparison failed before KWT; recovery uses exact archived mass, with geometry consistency bounded by64eps.'),
        registration_scope='Explicitly retrospective to exploratory algebra; formal results follow the revised registration. No claim of blind preregistration.',
        cases=cases,
        decision=dict(recommend_six_loaded_snapshots=True,
            purpose='Measure mesh sensitivity of the compatible-current/fixed-radius snapshot, induced interior response and residual adjacent normal current.',
            no_new_admission_threshold=True, stage3_closed=False, temporal_admission=False,
            full_D27_admission=False, normal_boundary_trace_admission=False, production=False),
        limitations=[
            'FD means Fermi–Dirac decoding here, not a new finite-difference or spectral validation.',
            'Recovered energies apply only to the unchanged archived thermal populations/fields; cannot serve as a spectrum at changed fields or arbitrary nonequilibrium populations.',
            'Boundary load is prescribed from Iref independently of the discrete gradient. Radius reactions are solved inside KWT and their work is explicit.',
            'Bias state changed from original Is=.99Iref to Ib=Is=Iref; the voltage/Joule change cannot be attributed solely to a numerical improvement.',
            'External normal injection zero follows the imposed contract; adjacent internal normal current is a measured discretization diagnostic, not proof of D.27 boundary-trace convergence.',
            'Fixed radius excludes da_b/dIs, reservoir kinetics, time integration and a detector pulse.'
        ])
    with (BASE/'reservoir_postprocess_review.json').open('x', encoding='utf-8') as stream:
        json.dump(review, stream, ensure_ascii=False, indent=2, allow_nan=False)
        stream.write('\n')
    c0, c1 = cases
    md = f'''# Revisión del postproceso con reservorios

El postproceso queda verificado como ensayo instantáneo: {len(sources)} fuentes, {len(inputs)} artefactos de entrada y los cuatro artefactos de salida de los casos coinciden con sus hashes. Los campos, poblaciones y masas usadas permanecen idénticos a los archivados. No se consultó un espectro nuevo ni se integró un paso temporal.

La recuperación E=theta[log(1-p)-log(p)] es la inversión de la preparación de Fermi–Dirac registrada. Se exige 0<p<0.5, energías positivas y ordenadas, y reproducción del momento de calor original. La revisión recompone además el depósito nuevo y CM9 directamente de los arrays guardados. Esta recuperación no permite usar esas energías para otros campos o poblaciones fuera del estado archivado.

La carga de fase procede de Iref, no del gradiente discreto. La reacción radial impone radio fijo dentro del solve KWT y mantiene explícito su trabajo, numéricamente nulo. El potencial se resuelve con Iref y el circuito usa Ib=Is=Iref, vc=0. El trabajo material de reservorio entra una vez en CM9 y no se suma al calor electrónico.

| Magnitud | Hélice | Perturbación suave |
|---|---:|---:|
| QDelta con extremos libres | {c0['condensate_heat_free_bar']:.9g} | {c1['condensate_heat_free_bar']:.9g} |
| QDelta con carga de reservorio | {c0['condensate_heat_loaded_bar']:.9g} | {c1['condensate_heat_loaded_bar']:.9g} |
| QDelta terminal con carga | {c0['condensate_heat_terminal_bar']:.9g} | {c1['condensate_heat_terminal_bar']:.9g} |
| Velocidad material máxima | {c0['maximum_material_speed_bar']:.9g} | {c1['maximum_material_speed_bar']:.9g} |
| Vdev [microvoltios] | {1e6*c0['Vdev_V']:.9g} | {1e6*c1['Vdev_V']:.9g} |

QDelta y velocidad usan las unidades normalizadas del registro. El trabajo de reservorio es -6.82937323 pW en ambos estados. CM9 recompuesta deja residuos de orden1e-27 W. El radio terminal cambia menos de1.4e-18 por unidad de tiempo normalizada. El defecto normal en la primera cara interior es aproximadamente -0.04193495% de Iref: se conserva como diagnóstico, no se fuerza a cero.

El descenso de QDelta revela la importancia de imponer una condición terminal compatible. El circuito también cambió de Is=0.99Iref a Is=Iref; por ello la comparación de voltaje o calor Joule no es una medida aislada de mejora numérica.

Recomiendo las seis instantáneas cargadas para estudiar sensibilidad espacial de QDelta, de la respuesta interior inducida y del defecto normal adyacente. El control libre permanece archivado. No se obtiene de ello un transiente, radio dinámico a_b(Is), intercambio cinético de reservorio ni admisión completa de D.27. La corriente normal externa nula es parte del contrato impuesto; la cara interior no certifica por sí sola esa traza continua.

El primer intento, fallido antes de KWT por comparar masas Linux/Windows bit a bit, permanece íntegro. La diferencia geométrica fue2.58e-15 relativa; el nuevo postproceso valida geometría dentro de64eps y usa exactamente las masas archivadas. Los criterios físicos no cambiaron. El registro declara la exploración algebraica anterior y no se presenta como un prerregistro ciego.

Reproducción sólo por aritmética guardada: `python sandbox/stage3_spatial/coupled_20260923/audit_saved_reservoir.py`. Rechaza sobrescribir esta revisión.
'''
    with (BASE/'reservoir_postprocess_review.md').open('x', encoding='utf-8') as stream:
        stream.write(md)
    print(json.dumps(dict(status=review['status'], issues=issues, source_count=len(sources), input_count=len(inputs)), indent=2))


if __name__ == '__main__':
    main()
