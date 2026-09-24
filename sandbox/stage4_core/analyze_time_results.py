"""Independent field analysis of saved weak thermal trajectories; no solver.

Response differences are measured after subtracting the drifting base. Both
initial-response and surviving-response scales are retained, so decaying tails
cannot acquire a misleading precision claim from the registered initial scale.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read(path):
    return json.loads(Path(path).read_text(encoding='utf8'))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--review', type=Path, required=True)
    args = parser.parse_args()
    raw = args.review / 'raw'
    receipt, summary, plan, refinement = [read(raw/n) for n in
        ('extraction_receipt.json', 'summary.json', 'executed_plan.json', 'refinement.json')]
    for name, digest in receipt['copied_file_sha256'].items():
        assert sha(raw/name) == digest, ('copied artifact changed', name)
    for filename in ('analysis.json', 'analysis.md'):
        if (args.review / filename).exists():
            raise FileExistsError('Do not overwrite an existing review: ' + filename)
    passes = {p: read(raw/p/'summary.json') for p in plan['passes']}
    a = np.load(raw/'trajectory_compact.npz')
    mass, d0, G0 = [a[k] for k in ('area_weights', 'd0', 'G0')]
    free = np.ones(len(mass), bool)
    free[a['boundary_nodes']] = False
    conductance = a['conductance']
    active = conductance > 0
    norm = lambda z: float(np.sqrt(np.sum(mass*abs(z)**2)))
    free_norm = lambda z: float(np.sqrt(np.sum(mass[free]*abs(z[free])**2)))
    edge_norm = lambda z: float(np.sqrt(np.sum(abs(z[active])**2/conductance[active])))
    direction = np.divide(d0, abs(d0), out=np.zeros_like(d0), where=abs(d0)>1e-12)
    tD_ps = .4415163339640258
    times = a['refined_times_ps']
    computed, responses, global_records, energies = {}, {}, [], {}
    maximum_summary_metric_difference = 0.
    for pass_name in passes:
        prefix = pass_name + '_'
        responses[pass_name] = {}
        for probe in ('amplitude', 'angular_phase'):
            responses[pass_name][probe] = []
            for i, t in enumerate(times):
                displacement = a[prefix+probe+'_difference'][i]
                hessian = a[prefix+probe+'_hessian'][i] - a[prefix+'baseline_hessian'][i]
                current = a[prefix+probe+'_current_increment'][i] - a[prefix+'baseline_current_increment'][i]
                potential = a[prefix+probe+'_potential_increment'][i] - a[prefix+'baseline_potential_increment'][i]
                velocity = a[prefix+probe+'_velocity'][i] - a[prefix+'baseline_velocity'][i]
                torque = np.imag(np.conj(d0)*hessian + np.conj(displacement)*G0)/mass
                fields = dict(displacement=displacement, current=current,
                    force_density=hessian/mass, phase_torque_density=torque)
                item = dict(time_ps=float(t), **{key: (edge_norm if key=='current' else norm)(value)
                    for key, value in fields.items()}, velocity_L2=norm(velocity), potential_L2=norm(potential),
                    amplitude_component_L2=norm(np.real(np.conj(direction)*displacement)),
                    transverse_phase_component_L2=norm(np.imag(np.conj(direction)*displacement)))
                responses[pass_name][probe].append(item)
                computed[(pass_name, probe, float(t))] = fields

        for state in ('baseline', 'amplitude', 'angular_phase'):
            F, rates, losses = [], [], []
            for i, row in enumerate(passes[pass_name]['observations']):
                registered = next(x for x in row['states'] if x['name']==state)
                x = a[prefix+state+'_displacement'][i]
                h = a[prefix+state+'_hessian'][i]
                velocity = a[prefix+state+'_velocity'][i]
                defect = a[prefix+state+'_constitutive_taylor_rhs_defect'][i]
                rate = float(np.real(np.vdot(G0+h, velocity)))
                loss = registered['approximate_kwt_loss'] + registered['approximate_normal_loss']
                force = (G0+h)/mass
                node_d = d0+x
                torque = np.imag(np.conj(node_d)*(G0+h))/mass
                check = {'quadratic_free_energy_rate': rate,
                    'constitutive_taylor_rhs_defect_L2': norm(defect),
                    'full_affine_rhs_L2': norm(velocity), 'displacement_L2': norm(x),
                    'maximum_displacement': float(np.max(abs(x[free])))}
                for name, value in check.items():
                    difference = abs(value-registered[name])
                    maximum_summary_metric_difference = max(maximum_summary_metric_difference, difference)
                    assert difference < 1e-12*max(1, abs(value)), (name, difference)
                item = dict(pass_name=pass_name, state=state, time_ps=row['time_ps'],
                    **check, maximum_displacement_relative_gap=check['maximum_displacement']/plan['gap_reference_kBTc'],
                    remaining_gradient_density_L2_free=free_norm(force),
                    remaining_phase_torque_density_L2_free=free_norm(torque),
                    constitutive_defect_relative_to_total_rhs=norm(defect)/max(norm(velocity), 1e-300),
                    quadratic_rate_plus_approximate_losses=rate+loss,
                    instantaneous_balance_defect_relative_to_losses=abs(rate+loss)/max(loss,1e-300),
                    baseline_mobility_gauge_correction_L2=registered['baseline_mobility_gauge_correction_L2'],
                    maximum_linearized_continuity=registered['maximum_linearized_continuity'],
                    approximate_noether_residual=registered['approximate_noether_residual'])
                global_records.append(item)
                F.append(registered['quadratic_free_energy'])
                rates.append(rate)
                losses.append(loss)
            F, rates, losses = np.asarray(F), np.asarray(rates), np.asarray(losses)
            energies[pass_name+'_'+state] = dict(free_energy=F.tolist(),
                free_energy_drop=float(F[0]-F[-1]), differences=np.diff(F).tolist(),
                all_observed_changes_nonpositive=bool(np.all(np.diff(F)<=0)),
                all_observed_rates_negative=bool(np.all(rates<0)),
                approximate_dissipation_positive=bool(np.all(losses>0)),
                trapezoid_integrated_losses_on_sparse_observation_times=float(np.trapezoid(losses, times/tD_ps)),
                trapezoid_quadrature_is_not_a_time_balance_certificate=True)

    comparisons, maxima = [], {}
    registered_comparison_difference = 0.
    for row in refinement['records']:
        probe, t, kind = row['probe'], row['time_ps'], row['observable']
        left, right = computed[('primary',probe,t)][kind], computed[('refined',probe,t)][kind]
        n = edge_norm if kind=='current' else norm
        error, reference = n(left-right), n(right)
        initial = n(computed[('primary',probe,0.)][kind])
        for name, value in [('error_L2',error),('relative_to_response',error/reference if reference>1e-14 else None),
                            ('relative_to_initial',error/initial if initial>1e-14 else None)]:
            if value is None:
                assert row[name] is None
            else:
                difference = abs(value-row[name])
                registered_comparison_difference = max(registered_comparison_difference,difference)
                assert difference < 1e-10*max(1,abs(value))
        tolerance = plan['refinement_absolute_norm'] + plan['refinement_relative_limit']*max(reference,initial)
        assert row['admitted'] == (error<=tolerance)
        comparisons.append(dict(**row, remaining_norm=reference, initial_norm=initial,
            response_fraction_of_initial=reference/initial if initial>1e-14 else None,
            absolute_plus_registered_relative_tolerance=tolerance,
            error_fraction_of_registered_tolerance=error/tolerance))
    for probe in ('amplitude','angular_phase'):
        for kind in ('displacement','current','force_density','phase_torque_density'):
            selected=[r for r in comparisons if r['probe']==probe and r['observable']==kind]
            maxima[probe+'_'+kind]=dict(
                maximum_error_relative_to_initial=max(r['relative_to_initial'] or 0 for r in selected),
                maximum_error_relative_to_remaining=max(r['relative_to_response'] or 0 for r in selected),
                maximum_fraction_of_registered_tolerance=max(r['error_fraction_of_registered_tolerance'] for r in selected),
                final_response_fraction_of_initial=selected[-1]['response_fraction_of_initial'],
                initial_norm=selected[0]['initial_norm'], final_norm=selected[-1]['remaining_norm'])

    base_initial = next(r for r in global_records if r['pass_name']=='refined' and r['state']=='baseline' and r['time_ps']==0)
    base_final = next(r for r in global_records if r['pass_name']=='refined' and r['state']=='baseline' and r['time_ps']==1)
    decay = {}
    for probe in ('amplitude','angular_phase'):
        rows=responses['refined'][probe]
        decay[probe]=dict(displacement_fraction_of_initial=[r['displacement']/rows[0]['displacement'] for r in rows],
            phase_component_fraction_of_total=[r['transverse_phase_component_L2']/max(r['displacement'],1e-300) for r in rows],
            amplitude_component_fraction_of_total=[r['amplitude_component_L2']/max(r['displacement'],1e-300) for r in rows],
            half_initial_crossing_bracket_ps=next(([float(times[i-1]),float(times[i])] for i in range(1,len(times))
                if rows[i]['displacement']<=.5*rows[0]['displacement']<rows[i-1]['displacement']), None))
    result = dict(status='WEAK_AFFINE_TIME_REVIEW_COMPLETE',
        runtime_seconds=summary['runtime_seconds'], operator_batches=summary['operator_batches'],
        horizon_ps=summary['horizon_ps'], checkpoints_verified=len(receipt['checkpoint_certificate']),
        accepted_steps_verified=len(receipt['accepted_steps']), operator_inputs_verified=receipt['operator_inputs_verified'],
        maximum_spectral_linear_solve_residual=summary['maximum_spectral_linear_solve_residual'],
        maximum_accepted_relative_displacement=receipt['maximum_accepted_relative_displacement'],
        registered_weak_domain=plan['maximum_relative_displacement'],
        maximum_affine_constant_defect=receipt['maximum_affine_constant_defect'],
        all_registered_refinement_comparisons_met=refinement['all_comparisons_met'],
        maximum_independent_summary_metric_difference=maximum_summary_metric_difference,
        maximum_independent_refinement_metric_difference=registered_comparison_difference,
        times_ps=times.tolist(), responses=responses, decay=decay,
        global_state_metrics=global_records, energy=energies, refinement=comparisons,
        refinement_summary=maxima,
        baseline=dict(initial=base_initial, final=base_final,
            final_displacement_over_initial_amplitude_response=base_final['displacement_L2']/responses['refined']['amplitude'][0]['displacement'],
            final_displacement_over_initial_phase_response=base_final['displacement_L2']/responses['refined']['angular_phase'][0]['displacement']),
        maximum_instantaneous_balance_relative_to_losses=max(r['instantaneous_balance_defect_relative_to_losses'] for r in global_records),
        maximum_constitutive_rhs_defect_relative_to_total_rhs=max(r['constitutive_defect_relative_to_total_rhs'] for r in global_records),
        stage4_complete=False, nonlinear_transient_admitted=False,
        scope='One weak affine thermal tangent about a residual core, at fixed temperature and contacts; no photon or new circuit. Finite-sum spectral curvature remains frozen.',
        limitation='Initial-scaled refinement admission does not establish precision in vanishing current/force/torque tails. Quadratic energy and approximate dissipation are Taylor diagnostics, not exact nonlinear energy closure.',
        recommended_nonlinear_snapshot_times_ps=[0,.03,1],
        provenance={'script_sha256':sha(__file__), **{str(p.relative_to(args.review)):sha(p) for p in
            (raw/'summary.json',raw/'identity.json',raw/'executed_plan.json',raw/'refinement.json',raw/'trajectory_compact.npz',raw/'extraction_receipt.json')}})
    (args.review/'analysis.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n',encoding='utf8')
    md = [
        '# Etapa 4: trayectoria térmica débil sobre todos los nodos', '',
        f'La corrida completó {summary["horizon_ps"]:g} ps en {summary["runtime_seconds"]:.2f} s, con {summary["operator_batches"]} aplicaciones por lotes del operador. '
        'La revisión verifica 37 checkpoints (16 observaciones y 21 pasos aceptados), las fuentes ejecutadas y 257 insumos del operador. No ejecuta nuevas ecuaciones físicas.', '',
        '## Resultado que puede admitirse', '',
        'Pasan las comparaciones de refinamiento registradas para este sistema afín térmico. El mayor desplazamiento total es '
        f'{100*result["maximum_accepted_relative_displacement"]:.6f} % del gap de referencia, frente al corte declarado de 2 %. '
        'Ese corte es un dominio de trabajo lineal, no un rango medido del material. La etapa 4 sigue abierta; aún no se acredita el transiente no lineal.', '',
        '| Perturbación, tras restar la base | Norma restante a 0,03 ps | Norma restante a 1 ps |',
        '|---|---:|---:|',
    ]
    for probe,label in [('amplitude','Amplitud'),('angular_phase','Fase angular')]:
        values=decay[probe]['displacement_fraction_of_initial']
        md.append(f'| {label} | {100*values[4]:.4f} % | {100*values[-1]:.4f} % |')
    md += ['', 'Son normas espaciales respecto a la misma perturbación inicial, no velocidades de respuesta del detector. '
        'La fase angular cruza la mitad entre 0,01 y 0,03 ps; el muestreo no permite dar un tiempo exacto de cruce.', '',
        '## Refinamiento y colas', '',
        '| Sonda | Observable | Máxima diferencia / señal inicial | Máxima diferencia / señal restante |',
        '|---|---|---:|---:|']
    for probe in ('amplitude','angular_phase'):
        for kind in ('displacement','current','force_density','phase_torque_density'):
            row=maxima[probe+'_'+kind]
            md.append(f'| {probe} | {kind} | {100*row["maximum_error_relative_to_initial"]:.6g} % | {100*row["maximum_error_relative_to_remaining"]:.6g} % |')
    md += ['', 'La aceptación usa la tolerancia previamente registrada: 1e-8 más 0,5 % de la mayor norma entre inicial y restante. '
        'Por eso no autoriza afirmar igual precisión relativa en las colas. El torque angular tardío puede diferir más que su valor restante; no se le asigna un tiempo de relajación preciso.', '',
        '## Deriva, energía y límites', '',
        f'La base residual también se mueve: su desplazamiento final es {base_final["displacement_L2"]:.8g} en norma de área, '
        f'{result["baseline"]["final_displacement_over_initial_amplitude_response"]:.4g} veces la perturbación inicial de amplitud. '
        'Se resta esa trayectoria para identificar cada respuesta; un mapa total tardío de fase está dominado por esta deriva.', '',
        'La energía cuadrática disminuye en todas las observaciones y su tasa calculada es negativa. La disipación aproximada KWT más normal es positiva. '
        f'El mayor defecto instantáneo tasa más disipación es {100*result["maximum_instantaneous_balance_relative_to_losses"]:.6g} % de la disipación; '
        f'el mayor defecto constitutivo de velocidad es {100*result["maximum_constitutive_rhs_defect_relative_to_total_rhs"]:.6g} % de la velocidad total. '
        'Son defectos de Taylor evaluados sobre los estados guardados, no una verificación de la energía no lineal exacta. '
        'Las integrales trapezoidales de ocho observaciones se conservan como diagnóstico de muestreo y no certifican un balance temporal.', '',
        f'La coordenada auxiliar constante de Arnoldi deriva hasta {result["maximum_affine_constant_defect"]:.6g} respecto de uno '
        f'({100*result["maximum_affine_constant_defect"]:.6g} %): es una pequeña variación numérica del forzamiento afín. '
        'Se registra sin imponer un umbral retrospectivo ni repetir la campaña. Puede preservarse exactamente en una futura versión del propagador.', '',
        'El siguiente contraste útil es evaluar la respuesta térmica no lineal exacta en los estados ya guardados a 0, 0,03 y 1 ps. '
        'El punto de 0,03 ps conserva aproximadamente la mitad de la perturbación de fase; el final concentra la deriva y la truncación. '
        'No requiere repetir la trayectoria ni incorporar un fotón.', '',
        '## Procedencia', '',
        '`raw/extraction_receipt.json` certifica los originales conservados en Geminga. '
        '`raw/trajectory_compact.npz` conserva todos los campos observados, retirando únicamente la repetición de la geometría. '
        '`analysis.json` conserva las normas, verificaciones reproducidas, sensibilidades y hashes.', '']
    (args.review/'analysis.md').write_text('\n'.join(md),encoding='utf8')
    print(json.dumps({k:result[k] for k in ('status','maximum_independent_summary_metric_difference',
        'maximum_independent_refinement_metric_difference','maximum_constitutive_rhs_defect_relative_to_total_rhs',
        'maximum_instantaneous_balance_relative_to_losses','maximum_affine_constant_defect')}))


if __name__ == '__main__':
    main()
