"""Reproduce the failed ETD2 comparison and distinguish its physical scales.

This is postprocessing, not a relaxation of the registered acceptance criterion.
All 96 original comparisons and the failed result are retained verbatim.
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
    raw = args.review/'raw'
    receipt, summary, plan, refinement = [read(raw/name) for name in
        ('extraction_receipt.json','summary.json','executed_plan.json','refinement.json')]
    for name, digest in receipt['copied_file_sha256'].items():
        assert sha(raw/name) == digest, ('copied artifact changed', name)
    for name in ('analysis.json','analysis.md'):
        if (args.review/name).exists():
            raise FileExistsError('Fresh analysis required: '+name)
    passes = {p: read(raw/p/'summary.json') for p in plan['passes']}
    a = np.load(raw/'trajectory_compact.npz')
    mass,d0,G0,c = [a[key] for key in ('area_weights','d0','G0','conductance')]
    active = c>0
    norm = lambda z: float(np.sqrt(np.sum(mass*abs(z)**2)))
    edge_norm = lambda z: float(np.sqrt(np.sum(abs(z[active])**2/c[active])))
    times = a['refined_times_ps']
    assert np.array_equal(times,a['primary_times_ps'])
    fields = {}
    base_torque = np.imag(np.conj(d0)*G0)/mass
    global_metrics,responses = [],{}
    for label in passes:
        responses[label] = {}
        for probe in ('baseline','amplitude','angular_phase'):
            observations=[]
            for index,t in enumerate(times):
                key=lambda state,kind:a[label+'_'+state+'_'+kind][index]
                x,g,j = [key(probe,k) for k in ('displacement','gap_gradient_increment','current_increment')]
                xb,gb,jb = [key('baseline',k) for k in ('displacement','gap_gradient_increment','current_increment')]
                total_torque = np.imag(np.conj(d0+x)*(G0+g))/mass
                if probe=='baseline':
                    values = dict(displacement=x,current=j,force_density=g/mass,
                        phase_torque_density=total_torque-base_torque)
                else:
                    values = dict(displacement=x-xb,current=j-jb,force_density=(g-gb)/mass,
                        phase_torque_density=np.imag(np.conj(d0+x)*(G0+g)-np.conj(d0+xb)*(G0+gb))/mass)
                fields[(label,probe,float(t))] = values
                states=passes[label]['observations'][index]['states']
                state=next(row for row in states if row['name']==probe)
                rate=float(np.real(np.vdot(G0+g,key(probe,'velocity'))))
                assert abs(rate-state['free_energy_rate'])<1e-14
                losses=state['kwt_loss']+state['normal_loss']
                observations.append(dict(time_ps=float(t),**{kind:(edge_norm if kind=='current' else norm)(value)
                    for kind,value in values.items()}))
                global_metrics.append(dict(pass_name=label,probe=probe,time_ps=float(t),
                    total_force_density_L2=norm((G0+g)/mass),total_phase_torque_density_L2=norm(total_torque),
                    maximum_relative_displacement=float(np.max(abs(x))/plan['gap_reference_kBTc']),
                    free_energy_change=state['free_energy_change'],free_energy_rate=rate,
                    kwt_loss=state['kwt_loss'],normal_loss=state['normal_loss'],power_residual=state['power_residual'],
                    power_residual_relative_to_loss=abs(state['power_residual'])/losses,
                    continuity_max=state['continuity_max'],noether_max=state['noether_max']))
            responses[label][probe]=observations
    comparisons=[];maximum_reproduction_difference=0.
    for registered in refinement['records']:
        probe,t,kind=registered['probe'],registered['time_ps'],registered['observable']
        n=edge_norm if kind=='current' else norm
        left,right=[fields[(p,probe,t)][kind] for p in ('primary','refined')]
        error,remaining=n(left-right),n(right)
        common_scale=max(n(fields[('refined',p,0.)][kind]) for p in ('amplitude','angular_phase'))
        scale=common_scale if probe=='baseline' else n(fields[('refined',probe,0.)][kind])
        for name,value in (('error_L2',error),('relative_to_initial',error/scale if scale>1e-14 else None),
                           ('relative_to_remaining',error/remaining if remaining>1e-14 else None)):
            if value is None:
                assert registered[name] is None
            else:
                difference=abs(value-registered[name]);maximum_reproduction_difference=max(maximum_reproduction_difference,difference)
                assert difference<1e-10*max(1,abs(value)),(probe,t,kind,name,value,registered[name])
        tolerance=plan['refinement_absolute_norm']+plan['refinement_relative_limit']*max(scale,remaining)
        assert registered['admitted']==bool(error<=tolerance)
        maximum_observed=max(n(fields[('refined',probe,float(time))][kind]) for time in times)
        comparisons.append(dict(**registered,initial_norm=scale,remaining_norm=remaining,
            primary_remaining_norm=n(left),initial_same_observable_all_probe_scale=common_scale,
            error_relative_to_common_probe_scale=error/common_scale if common_scale>0 else None,
            maximum_observed_signal_norm=maximum_observed,error_relative_to_maximum_observed_signal=error/maximum_observed if maximum_observed>0 else None,
            remaining_fraction_of_initial=remaining/scale if scale>0 else None,
            registered_absolute_plus_relative_tolerance=tolerance,error_fraction_of_registered_tolerance=error/tolerance))
    failed=[row for row in comparisons if not row['admitted']]
    final_time=float(times[-1]);maxima={}
    for probe in ('baseline','amplitude','angular_phase'):
        for kind in ('displacement','current','force_density','phase_torque_density'):
            selected=[row for row in comparisons if row['probe']==probe and row['observable']==kind]
            maxima[probe+'_'+kind]=dict(maximum_error_relative_to_initial=max(row['relative_to_initial'] or 0 for row in selected),
                maximum_error_relative_to_remaining=max(row['relative_to_remaining'] or 0 for row in selected),
                maximum_error_relative_to_common_probe_scale=max(row['error_relative_to_common_probe_scale'] or 0 for row in selected),
                final_response_fraction_of_initial=selected[-1]['remaining_fraction_of_initial'],
                initial_norm=selected[0]['initial_norm'],final_norm=selected[-1]['remaining_norm'])
    history={label:dict(attempts=len(data['history']),accepted=sum(row['accepted'] for row in data['history']),
        rejected=sum(not row['accepted'] for row in data['history']),runtime_seconds=data['runtime_seconds'],
        maximum_embedded_correction=max(row['embedded_correction'] or 0 for row in data['history']),
        maximum_first_krylov_dimension=max(row['first_dimension'] for row in data['history']),
        maximum_second_krylov_dimension=max(row['second_dimension'] or 0 for row in data['history']))
        for label,data in passes.items()}
    energies={}
    for label in passes:
        for probe in ('baseline','amplitude','angular_phase'):
            rows=[row for row in global_metrics if row['pass_name']==label and row['probe']==probe]
            values=np.array([row['free_energy_change'] for row in rows])
            energies[label+'_'+probe]=dict(free_energy_change=values.tolist(),
                differences=np.diff(values).tolist(),all_observed_changes_nonpositive=bool(np.all(np.diff(values)<=0)),
                all_observed_rates_negative=all(row['free_energy_rate']<0 for row in rows),
                maximum_power_residual_relative_to_loss=max(row['power_residual_relative_to_loss'] for row in rows))
    spatial_failure=[]
    radius=np.linalg.norm(a['coordinates_bar'],axis=1)
    for row in failed:
        error=fields[('primary',row['probe'],row['time_ps'])][row['observable']]-fields[('refined',row['probe'],row['time_ps'])][row['observable']]
        if row['observable']=='current':
            continue
        density=mass*abs(error)**2
        spatial_failure.append(dict(probe=row['probe'],time_ps=row['time_ps'],observable=row['observable'],
            cumulative_squared_error_fraction_by_radius_ell0={str(r):float(np.sum(density[radius<=r])/np.sum(density)) for r in (1,2,4,8,12)},
            maximum_error_density=float(np.max(abs(error))),maximum_error_coordinate_ell0=a['coordinates_bar'][np.argmax(abs(error))].tolist()))
    identical_times=[float(time) for time in times if all(np.array_equal(fields[('primary',probe,float(time))][kind],fields[('refined',probe,float(time))][kind])
        for probe in ('baseline','amplitude','angular_phase') for kind in ('displacement','current','force_density','phase_torque_density'))]
    distinct_comparisons=[row for row in comparisons if row['time_ps'] not in identical_times]
    plot_metadata=dict(
        physical_case='Fixed bath T=0.9 K, Tc=8.65 K; 65x65 graph vortex-core control on square [-6,6] ell0 in x and y, spacing 0.1875 ell0 and 256 fixed boundary nodes; 256 Matsubara terms; no photon, no strip experiment and no circuit trajectory.',
        baseline='d_b(t) is the unperturbed evolving core. It is not stationary and is subtracted at the same time.',
        amplitude_probe='Initial d_A-d_b = 0.001 d0 b(r), b=max(0,1-r^2/R^2)^3; R=4 ell0 from admitted operator plan.',
        angular_phase_probe='Initial d_P-d_b = 0.001 i d0 (x/R) b(r). Cartesian perturbation, not a global phase shift.',
        node_norm='||z||_M=sqrt(sum_i m_i |z_i|^2), including fixed nodes, m_i dimensionless cell area.',
        current_norm='||j||=sqrt(sum_edges_positive_c j_e^2/c_e), c_e graph conductance in dimensionless action.',
        displacement='d=Delta/(k_B Tc). Curve: ||d_probe(t)-d_b(t)||_M / ||d_probe(0)-d_b(0)||_M; unitless, not a point amplitude.',
        force_density='q_i=G_i/m_i, where G is discrete thermal free-energy gradient. Curve uses probe-minus-base force; dimensionless variational density, not newtons.',
        phase_torque_density='tau_i=Im(conj(d_i) G_i)/m_i. Curve uses tau_probe-tau_base; dimensionless force for condensate phase, not mechanical torque.',
        current='Derivative of the same dimensionless spectral action with respect to edge gauge link. Probe-minus-base graph current; no ampere conversion is asserted.',
        time='Real physical time in ps; tau=t/[hbar/(2 k_B Tc)] is internal solver time only.',
        refinement='Primary minus refined trajectory in the same observable and norm. Original gate: error <= 1e-8 + 0.005 max(initial norm, remaining refined norm).',
        contextual_scale='Maximum initial amplitude/phase signal for the SAME observable. A diagnostic contextualization, not a replacement gate.',
        free_energy='Change in dimensionless finite-Matsubara thermal free-energy action from d0 at fixed bath. It is not conserved internal energy.',
        plots_must_distinguish='Response magnitude, refinement difference, and remaining-signal denominator. State which field, probe, norm, normalization, units and baseline appear in each caption.')
    result=dict(status='NONLINEAR_TIME_FAILURE_REPRODUCED_NO_SOLVES',original_status=summary['status'],
        failure_reason=read(raw/'failure.json')['reason'],runtime_seconds=summary['runtime_seconds'],
        horizon_ps=summary['horizon_ps'],spectral_roots=summary['spectral_roots'],nonlinear_batches=summary['nonlinear_batches'],linear_batches=summary['linear_batches'],
        maximum_spectral_residual=summary['maximum_spectral_residual'],maximum_linear_residual=summary['maximum_linear_residual'],
        checkpoints_verified=len(receipt['checkpoint_certificate']),accepted_steps_verified=len(receipt['accepted_steps']),
        original_gate=dict(relative_limit=plan['refinement_relative_limit'],absolute_norm=plan['refinement_absolute_norm'],changed=False,
            passed_comparisons=len(comparisons)-len(failed),failed_comparisons=len(failed),all_passed=False),
        temporal_refinement_scope=dict(bitwise_identical_observation_times_ps=identical_times,
            comparisons_of_identical_trajectories=12*len(identical_times),
            comparisons_of_distinct_trajectories=len(distinct_comparisons),
            passed_distinct_comparisons=sum(row['admitted'] for row in distinct_comparisons),
            interpretation='Both passes share their first 14 accepted steps through 0.16 ps. The maximum-step bounds 0.1 ps and 0.05 ps then separate the grids. Only observations 0.3 ps and 1 ps compare distinct trajectories; this is not a systematic convergence-order study.'),
        implementation_observation=dict(primary_final_step_ps=passes['primary']['history'][-1]['step_ps'],
            recommendation='Future drivers should identify the observation endpoint within floating-point rounding and avoid a redundant last 1.11e-16 ps step; original results are retained.'),
        maximum_registered_metric_reproduction_difference=maximum_reproduction_difference,
        comparisons=comparisons,failed_comparisons=failed,maxima=maxima,responses=responses,global_metrics=global_metrics,
        energies=energies,history=history,spatial_failure=spatial_failure,plot_metadata=plot_metadata,
        original_failure_preserved=True,new_physics_solves=0,analysis_source_sha256=sha(__file__),
        raw_extraction_receipt_sha256=sha(raw/'extraction_receipt.json'))
    with (args.review/'analysis.json').open('x',encoding='utf8') as stream:
        json.dump(result,stream,indent=2,allow_nan=False);stream.write('\n')
    f=failed[0]
    lines=['# Diagnostico de la corrida termica no lineal','',
        f'Las dos trayectorias llegaron a {final_time:g} ps en {summary["runtime_seconds"]:.2f} s. La interrupcion fue el control posterior de refinamiento: {len(comparisons)-len(failed)}/{len(comparisons)} comparaciones cumplen y {len(failed)} no. No fallo Newton ni la integracion.',
        'Hasta 0.1 ps ambas salidas son identicas porque utilizaron los mismos pasos. De las 24 comparaciones en tiempos con trayectorias diferentes (0.3 y 1 ps), 23 cumplen. No son 96 pruebas independientes ni un estudio de orden de convergencia.','',
        '## Comparacion que no cumple','',
        f'- Campo: torque variacional de fase inducido por la perturbacion radial de amplitud, restando el nucleo base al mismo tiempo; t={f["time_ps"]:g} ps.',
        f'- Diferencia primaria-refinada: {f["error_L2"]:.9g}; tolerancia registrada: {f["registered_absolute_plus_relative_tolerance"]:.9g}.',
        f'- Escala inicial de ese torque: {f["initial_norm"]:.9g}; diferencia: {100*f["relative_to_initial"]:.6g}% de esa escala.',
        f'- Torque restante refinado: {f["remaining_norm"]:.9g}; diferencia: {100*f["relative_to_remaining"]:.6g}% de la senal restante.',
        f'- Diferencia contextual respecto al torque inicial principal (sonda angular): {100*f["error_relative_to_common_probe_scale"]:.6g}%. Esta escala adicional no cambia el resultado registrado.',
        '', '## Alcance fisico de la evidencia','',
        'El torque de fase generado por una sonda de amplitud es una componente cruzada pequena. Su precision relativa tardia no esta acreditada; no debe presentarse como un resultado cuantitativo resuelto. La amplitud, la corriente y la fuerza principales si concuerdan entre los dos pasos de tiempo. Esto permite continuar el desarrollo pertinente al sistema final sin repetir esta bateria solo para resolver una cola pequena. El estado historico sigue siendo TEMPORAL_REFINEMENT_NOT_MET.',
        '', '## Magnitudes y figuras','',*['- '+key+': '+value for key,value in plot_metadata.items()],
        '',f'Verificados {len(receipt["checkpoint_certificate"])} archivos NPZ y {len(receipt["accepted_steps"])} pasos aceptados. Reproduccion independiente desde campos: diferencia maxima {maximum_reproduction_difference:.3g}. No se ejecutaron nuevas ecuaciones fisicas.','']
    (args.review/'analysis.md').write_text('\n'.join(lines),encoding='utf8')
    print(json.dumps(dict(status=result['status'],failed=result['failed_comparisons'],maxima=maxima,history=history,
        maximum_power_residual_relative_to_loss=max(row['power_residual_relative_to_loss'] for row in global_metrics)),indent=2))


if __name__=='__main__':
    main()
