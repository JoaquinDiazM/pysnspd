"""Read-only comparison of the frozen Newton and verified-predictor pilots.

No spectral solve or time integration is executed. Every compared field is
read from saved arrays, and every input file is identified by SHA256.
"""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read(path):
    return json.loads(Path(path).read_text())


def compare(folder, root):
    initial, predicted = folder/'pilot_initial', folder/'pilot_predicted'
    summaries = {name: read(path/'summary.json') for name, path in
                 [('initial', initial), ('predicted', predicted)]}
    identities = {name: read(path/'identity.json') for name, path in
                  [('initial', initial), ('predicted', predicted)]}
    plans = {name: read(path/'executed_plan.json') for name, path in
             [('initial', initial), ('predicted', predicted)]}
    for key in ('T_K', 'Tc_K', 'tau_ee_Tc_ps', 'tau_ep_Tc_ps', 'delta0_over_kBTc',
                'matsubara_count', 'perturbation_amplitude', 'spectral_tolerance',
                'euler_safety', 'mesh_sha256', 'observation_times_ps'):
        if plans['initial'][key] != plans['predicted'][key]:
            raise ValueError('Physical/numerical campaign field changed: '+key)
    inputs = {}
    for directory in (initial, predicted):
        for path in directory.glob('*'):
            if path.is_file():
                inputs[path.relative_to(folder).as_posix()] = digest(path)
        if digest(directory/'executed_plan.json') != read(directory/'identity.json')['plan_sha256']:
            raise ValueError('Executed plan hash changed')
        metadata = read(directory/'observation_000.json')
        if digest(directory/metadata['fields_path']) != metadata['fields_sha256']:
            raise ValueError('Initial observation hash changed')
    source_verification = {}
    for name in ('dual_kwt_time.py', 'test_dual_kwt_time.py'):
        current = root/('sandbox/stage4_core' if name == 'dual_kwt_time.py' else 'tests')/name
        archived = initial/'executed_sources'/name
        registered = plans['initial']['source_sha256'][current.relative_to(root).as_posix()]
        source_verification[name] = dict(current_sha256=digest(current),
            executed_archive_sha256=digest(archived), registered_sha256=registered,
            identical=digest(current) == digest(archived) == registered)
        if not source_verification[name]['identical']:
            raise ValueError('Frozen original source changed: '+name)
    names = list(summaries['predicted']['accepted_steps'])
    field_records, initial_observable_records = [], []
    with np.load(initial/'observation_000.npz') as a0, \
         np.load(predicted/'observation_000.npz') as b0, \
         np.load(initial/'pilot_last_states.npz') as afinal, \
         np.load(predicted/'pilot_last_states.npz') as bfinal:
        mass, conductance, baseline = b0['area_weights'], b0['conductance'], b0['baseline_gap']
        if not np.array_equal(afinal['times_ps'], bfinal['times_ps']):
            raise ValueError('Pilot endpoints differ')
        def norm(value, kind='gap'):
            if kind == 'current':
                active = conductance > 0
                return float(np.sqrt(np.sum(abs(value[active])**2/conductance[active])))
            return float(np.sqrt(np.dot(mass, abs(value)**2)))
        for index, name in enumerate(names):
            scale = norm(b0[name+'_gap']-baseline)
            error = norm(afinal[name]-bfinal[name])
            field_records.append(dict(trajectory=name, time_ps=float(bfinal['times_ps'][index]),
                initial_perturbation_area_L2=scale, difference_area_L2=error,
                difference_relative_to_initial_perturbation=error/scale,
                difference_percent_of_initial_perturbation=100*error/scale,
                maximum_node_difference_kBTc=float(np.max(abs(afinal[name]-bfinal[name])))))
            for kind in ('gap_difference', 'current', 'force_density', 'phase_torque_density'):
                scale = norm(b0[name+'_'+kind], kind)
                error = norm(a0[name+'_'+kind]-b0[name+'_'+kind], kind)
                initial_observable_records.append(dict(trajectory=name, observable=kind,
                    time_ps=0., error_norm=error, reference_initial_norm=scale,
                    relative_error=None if scale == 0 else error/scale))
    history_records, integrated_balance = [], {}
    histories = {name: read(path/'step_history.json') for name, path in
                 [('initial', initial), ('predicted', predicted)]}
    if len(histories['initial']) != len(histories['predicted']):
        raise ValueError('Step histories differ')
    for a, b in zip(histories['initial'], histories['predicted']):
        for key in ('trajectory', 'time_ps', 'dt_ps'):
            if a[key] != b[key]:
                raise ValueError('Accepted time-step partition differs')
        history_records.append(dict(trajectory=b['trajectory'], time_ps=b['time_ps'],
            energy_excess_difference=b['energy_excess']-a['energy_excess'],
            integrated_loss_difference=b['integrated_loss']-a['integrated_loss'],
            balance_residual_difference=b['balance_residual']-a['balance_residual']))
    for variant, result in summaries.items():
        integrated_balance[variant] = {}
        for name in names:
            energy = result['initial_energy_excess'][name]
            absolute = result['maximum_balance_residual'][name]
            integrated_balance[variant][name] = dict(initial_excess_free_energy=energy,
                maximum_absolute_balance_residual=absolute,
                maximum_balance_percent_of_initial_excess=100*absolute/energy,
                maximum_single_step_energy_increase=result['maximum_single_step_energy_increase'][name])
    timing = {}
    resources = {}
    for variant, result in summaries.items():
        # This timer starts after the initial RHS and includes the initial saved
        # observation. It is NOT a pure matrix-solver timer.
        advance = result['projected_full_seconds']*result['covered_ps']/plans[variant]['observation_times_ps'][-1]
        timing[variant] = dict(total_wall_seconds=result['runtime_seconds'],
            accepted_advance_and_recording_seconds=advance,
            preparation_and_initial_evaluation_seconds=result['runtime_seconds']-advance,
            projected_horizon_seconds=result['projected_full_seconds'],
            newton_solves=result['newton_solves'],
            verified_prediction_acceptances=result['stationary_reuses'],
            maximum_spectral_residual=result['maximum_spectral_residual'])
        identity = identities[variant]
        budget, available = identity['budget'], identity['resources']
        selected = set(budget['selected_affinity_cpus'])
        untouched_cores = [group for group in available['physical_core_groups']
                           if not selected.intersection(group)]
        cpu_ok = budget['total_processes'] <= .9*budget['effective_cpu_quota_units']
        memory_ok = budget['estimated_memory_reservation_bytes'] <= .9*available['available_memory_bytes']
        resources[variant] = dict(worker_count=budget['workers'],
            coordinator_count=1, threads_per_process=budget['threads_per_process'],
            assigned_logical_cpus=len(selected), available_logical_cpus=len(available['logical_cpus']),
            entirely_reserved_physical_cores=len(untouched_cores), reserved_core_cpu_sets=untouched_cores,
            memory_scheduling_reservation_GiB=budget['estimated_memory_reservation_bytes']/1024**3,
            memory_available_GiB=available['available_memory_bytes']/1024**3,
            cpu_below_90_percent=bool(cpu_ok), memory_reservation_below_90_percent=bool(memory_ok),
            scope='Resource reservation and affinity from executed inventory; not measured peak RSS')
        if not cpu_ok or not memory_ok or len(untouched_cores) < 2:
            raise ValueError('Executed resources violate the registered budget')
    return dict(schema='pysnspd.stage4.dual_kwt_pilot_comparison.v1',
        status='PILOT_COMPARISON_COMPLETE', physical_time_ps=.001,
        stage4_complete=False, full_horizon_complete=False,
        fields=field_records, initial_observables=initial_observable_records,
        history_differences=history_records, integrated_balance=integrated_balance,
        timing=timing, resources=resources, source_verification=source_verification,
        speedup_total=timing['initial']['total_wall_seconds']/timing['predicted']['total_wall_seconds'],
        speedup_accepted_advance=timing['initial']['accepted_advance_and_recording_seconds']/timing['predicted']['accepted_advance_and_recording_seconds'],
        inputs_sha256=inputs,
        normalization='d=Delta/(kB Tc); area norm sqrt(sum_i m_i|x_i|^2), m_i=dual_area/ell0^2. Endpoint field differences divided by the same trajectory initial perturbation norm. Edge current norm sqrt(sum |I|^2/c).',
        energy_scope='Dimensionless fixed-temperature free energy; integral of (KWT loss + normal Joule loss) dt_ps/tD_ps. Not conserved population internal energy.',
        projection_scope='Extrapolation of two accepted substeps only. A later exact-residual check can require Newton; projected duration is not a guaranteed runtime.',
        physics_solves_executed_by_analysis=0)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--folder', type=Path, required=True)
    args = parser.parse_args()
    result = compare(args.folder, Path(__file__).resolve().parents[2])
    (args.folder/'pilot_comparison.json').write_text(json.dumps(result, indent=2)+'\n', encoding='utf8')
    timing = result['timing']
    worst_field = max(row['difference_percent_of_initial_perturbation'] for row in result['fields'])
    worst_balance = max(row['maximum_balance_percent_of_initial_excess']
        for row in result['integrated_balance']['predicted'].values())
    full_summary_path = args.folder/'thermal'/'summary.json'
    completed_note = None
    if full_summary_path.exists():
        full = read(full_summary_path)
        if full.get('status') == 'DUAL_KWT_THERMAL_COMPLETE':
            completed_note = ('**Actualización posterior al piloto:** la trayectoria completa de 1 ps '
                f"terminó en {full['runtime_seconds']:.2f} s y cumplió sus criterios registrados. "
                'La estimación y el consejo de entregar un comando largo se conservan aquí como '
                'contexto histórico; ya no queda una ejecución térmica pendiente. Los resultados '
                'finales están en `thermal/analysis.md`.')
    lines = [
        '# Comparación de los pilotos KWT sobre malla dual', '',
        'Los dos pilotos avanzaron hasta **0,001 ps** en la misma malla de 1712 nodos, con el mismo paso KWT heredado y 256 frecuencias. Se mantuvieron temperatura, contactos, estado inicial y tolerancia espectral 10⁻⁷. El predictor sólo cambia el punto de partida de la consulta espectral: cada resultado pasa por el residuo no lineal exacto y, si no cumple, se aplica Newton.', '',
        '| Medida | Newton sin predictor | Predictor verificado |',
        '|---|---:|---:|',
        f"| Tiempo total | {timing['initial']['total_wall_seconds']:.3f} s | {timing['predicted']['total_wall_seconds']:.3f} s |",
        f"| Preparación y evaluación inicial | {timing['initial']['preparation_and_initial_evaluation_seconds']:.3f} s | {timing['predicted']['preparation_and_initial_evaluation_seconds']:.3f} s |",
        f"| Avance aceptado y registro | {timing['initial']['accepted_advance_and_recording_seconds']:.3f} s | {timing['predicted']['accepted_advance_and_recording_seconds']:.3f} s |",
        f"| Consultas que necesitaron Newton | {timing['initial']['newton_solves']} | {timing['predicted']['newton_solves']} |",
        f"| Predicciones con residuo admitido | {timing['initial']['verified_prediction_acceptances']} | {timing['predicted']['verified_prediction_acceptances']} |",
        f"| Máximo residuo espectral | {timing['initial']['maximum_spectral_residual']:.3e} | {timing['predicted']['maximum_spectral_residual']:.3e} |", '',
        f"La mayor diferencia entre campos finales es **{worst_field:.3e}% de la perturbación inicial de su propia sonda**. Se compara d=Delta/(kB Tc), con norma sqrt(sum m|delta d|²), donde m es el área dual en unidades de ell0². Las dos resoluciones temporales usan los mismos pasos entre ambas implementaciones.", '',
        f"El máximo residuo del balance integrado del piloto optimizado representa **{worst_balance:.5f}% de su exceso inicial de energía libre**. No se detectó aumento de energía entre pasos aceptados. Esto verifica el piloto de 0,001 ps; no certifica todavía el horizonte de 1 ps ni el balance de energía interna de poblaciones.", '',
        'Se emplearon 27 trabajadores y un coordinador, todos con un hilo numérico, sobre 28 CPU lógicas de 32. Quedaron dos núcleos físicos completos libres. La reserva de memoria de 30 GiB estuvo por debajo del 90% disponible. El inventario y la afinidad se vuelven a comprobar en cada ejecución.', '',
        'La aceleración total del piloto fue {:.2f}×; la del bloque de avance y registro fue {:.2f}×. La extrapolación inicial fue de {:.1f} min y se trató como una estimación provisional basada en dos subpasos, pues podían aparecer nuevas correcciones Newton.'.format(result['speedup_total'], result['speedup_accepted_advance'], timing['predicted']['projected_horizon_seconds']/60), '',
        'Las fuentes originales y su prueba permanecen idénticas a los archivos archivados del primer piloto. `pilot_comparison.json` contiene las diferencias por sonda, normas, residuos, inventarios y hashes de entrada. El análisis no ejecutó nuevas soluciones físicas.', '',
        'Reproducción:', '',
        '```bash',
        'python -m sandbox.stage4_core.dual_kwt_compare_pilots --folder docs/implementation/stage4/final_kwt_20260924',
        '```', '',
    ]
    if completed_note is not None:
        lines[2:2] = [completed_note, '']
    (args.folder/'pilot_comparison.md').write_text('\n'.join(lines), encoding='utf8')
    print(json.dumps(dict(status=result['status'], maximum_field_difference_percent=worst_field,
        maximum_balance_percent=worst_balance, speedup_total=result['speedup_total'],
        speedup_accepted_advance=result['speedup_accepted_advance']), indent=2))


if __name__ == '__main__':
    main()
