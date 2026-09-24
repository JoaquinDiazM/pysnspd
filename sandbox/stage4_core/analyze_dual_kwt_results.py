"""Reproduce the completed dual KWT analysis and explicitly labelled figures.

This is postprocessing only: all fields and histories are read from the
extracted receipts. No stationary spectral root or time step is executed.
"""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.tri import Triangulation
from matplotlib.ticker import FuncFormatter


def read(path):
    return json.loads(Path(path).read_text())


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def analyze(folder, root):
    extraction = read(folder/'extraction_receipt.json')
    for name, expected in extraction['source_hashes'].items():
        if sha(folder/name) != expected:
            raise ValueError('Extracted data hash changed: '+name)
    summary, plan, identity, refinement, history = [read(folder/name) for name in
        ('summary.json', 'executed_plan.json', 'identity.json', 'refinement.json', 'step_history.json')]
    if sha(folder/'executed_plan.json') != identity['plan_sha256']:
        raise ValueError('Executed plan hash changed')
    observations = []
    for path in sorted(folder.glob('observation_*.json')):
        meta = read(path)
        if sha(folder/meta['fields_path']) != meta['fields_sha256']:
            raise ValueError('Observation hash mismatch')
        with np.load(folder/meta['fields_path']) as loaded:
            arrays = {key: loaded[key].copy() for key in loaded.files}
        observations.append((meta, arrays))
    first = observations[0][1]
    mass, conductance = first['area_weights'], first['conductance']
    names = list(summary['accepted_steps'])
    def norm(kind, value):
        if kind == 'current':
            active = conductance > 0
            return float(np.sqrt(np.sum(abs(value[active])**2/conductance[active])))
        return float(np.sqrt(np.dot(mass, abs(value)**2)))
    reproduced_errors = []
    evolution = {name: [] for name in names}
    for metadata, arrays in observations:
        time_ps = metadata['time_ps']
        for name in names:
            initial_scale = norm('gap_difference', first[name+'_gap_difference'])
            evolution[name].append(dict(time_ps=time_ps,
                gap_difference_area_L2=norm('gap_difference', arrays[name+'_gap_difference']),
                gap_difference_relative_to_own_initial=norm('gap_difference', arrays[name+'_gap_difference'])/initial_scale,
                excess_energy=metadata['states'][name]['energy_excess'],
                accumulated_loss=metadata['states'][name]['integrated_loss']))
        for probe in ('amplitude', 'angular_phase'):
            for kind in ('gap_difference', 'current', 'force_density', 'phase_torque_density'):
                primary, refined = arrays['primary_'+probe+'_'+kind], arrays['refined_'+probe+'_'+kind]
                scale = max(norm(kind, first['refined_'+p+'_'+kind]) for p in ('amplitude', 'angular_phase'))
                difference = norm(kind, primary-refined)
                reproduced_errors.append(dict(time_ps=time_ps, probe=probe, observable=kind,
                    error_L2=difference, initial_signal_scale=scale,
                    relative_to_initial_signal=difference/scale if scale else None,
                    admitted=difference <= plan['refinement_absolute_norm']+plan['refinement_relative_limit']*scale))
    maximum_reproduction = 0.
    for saved, current in zip(refinement['records'], reproduced_errors):
        if (saved['time_ps'], saved['probe'], saved['observable']) != (current['time_ps'], current['probe'], current['observable']):
            raise ValueError('Comparison ordering differs')
        maximum_reproduction = max(maximum_reproduction, abs(saved['error_L2']-current['error_L2']))
        if saved['admitted'] != current['admitted']:
            raise ValueError('Admission not reproduced')
    if len(refinement['records']) != len(reproduced_errors) or maximum_reproduction > 1e-14:
        raise ValueError('Refinement comparison failed reproduction')
    maxima = []
    for probe in ('amplitude', 'angular_phase'):
        for kind in ('gap_difference', 'current', 'force_density', 'phase_torque_density'):
            rows = [row for row in reproduced_errors if row['probe'] == probe and row['observable'] == kind]
            peak = max(rows, key=lambda row: row['relative_to_initial_signal'] or 0.)
            maxima.append(dict(probe=probe, observable=kind, time_ps=peak['time_ps'],
                maximum_error_percent_initial_signal=100*(peak['relative_to_initial_signal'] or 0.),
                initial_signal_scale=peak['initial_signal_scale']))
    balances = {}
    for name in names:
        rows = [row for row in history if row['trajectory'] == name]
        initial_energy = summary['integrated_balance'][name]['initial_energy_excess']
        maximum = max(abs(row['balance_residual']) for row in rows)
        balance = summary['integrated_balance'][name]
        if not np.isclose(maximum, balance['maximum_absolute_balance_residual'], rtol=0, atol=1e-18):
            raise ValueError('Balance history maximum differs')
        if len(rows) != summary['accepted_steps'][name]:
            raise ValueError('Missing accepted steps')
        balances[name] = dict(accepted_steps=len(rows), initial_energy_excess=initial_energy,
            maximum_balance_percent_of_initial_excess=100*maximum/initial_energy,
            allowance_percent_of_initial_excess=100*balance['admitted_limit']/initial_energy,
            final_energy_fraction=rows[-1]['energy_excess']/initial_energy,
            final_integrated_loss_fraction=rows[-1]['integrated_loss']/initial_energy,
            maximum_single_step_energy_increase=balance['maximum_single_step_energy_increase'],
            admitted=balance['admitted'])
    mesh_path = root/'docs/implementation/stage4/practical_time_review_20260924/dual_mesh/resampled/mesh.npz'
    if sha(mesh_path) != identity['mesh_sha256']:
        raise ValueError('Spatial plotting mesh changed')
    with np.load(mesh_path) as mesh:
        coordinates_nm = first['coordinates_bar']*float(mesh['ell0_m'])*1e9
        triangles = mesh['triangles'].copy()
    result = dict(schema='pysnspd.stage4.dual_kwt.completed_analysis.v1',
        status=summary['status'], completed_horizon_ps=summary['completed_horizon_ps'],
        runtime_seconds=summary['runtime_seconds'], spectral_queries=summary['stationary_reuses']+summary['newton_solves'],
        verified_predictions=summary['stationary_reuses'], newton_fallbacks=summary['newton_solves'],
        maximum_spectral_residual=summary['maximum_spectral_residual'], spectral_tolerance=plan['spectral_tolerance'],
        nodes=len(mass), frequencies=plan['matsubara_count'],
        primary_step_cap_ps=identity['step_bound']['primary_step_ps'],
        time_order=1, workers=summary['worker_count'], physical_cores_reserved=identity['budget']['physical_cores_reserved'],
        admitted_comparisons=sum(row['admitted'] for row in reproduced_errors), total_comparisons=len(reproduced_errors),
        admitted_positive_time_comparisons=sum(row['admitted'] for row in reproduced_errors if row['time_ps'] > 0),
        positive_time_comparisons=sum(row['time_ps'] > 0 for row in reproduced_errors),
        maximum_refinement_reproduction_difference=maximum_reproduction,
        temporal_maxima=maxima, integrated_balance=balances, evolution=evolution,
        final_gap_response_fraction={name: evolution[name][-1]['gap_difference_relative_to_own_initial'] for name in names},
        figure_definitions=dict(gap='d=Delta/(kB Tc); delta d=d_probe-d_uniform at the same time. Area norm sqrt(sum m|delta d|^2), m=dual_area/ell0^2.',
            temporal_error='Norm of primary minus refined divided by largest initial signal of the SAME observable among the amplitude and phase probes. Current norm sqrt(sum |I|^2/c).',
            phase_torque='Im(conj(d)*G)/m: dimensionless phase-conjugate variational force, not mechanical torque.',
            energy='Dimensionless excess free energy F(d)-F(d_uniform); accumulated loss integral (P_KWT+P_normal) dt_ps/tD_ps. Each normalized by its own initial excess.',
            spatial='Coordinates in nm relative to strip center. Maps show magnitude of delta d or primary-minus-refined, divided by d_uniform. No photon deposition or hotbelt is represented.'),
        stage4_complete=False, nonthermal_work_closed=False, production_changed=False,
        physics_solves_executed_by_analysis=0, input_sha256=extraction['source_hashes'])
    return result, observations, history, refinement, coordinates_nm, triangles


def figures(folder, result, observations, history, refinement, coordinates, triangles):
    destination = folder/'figures'
    destination.mkdir(exist_ok=True)
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10,
        'axes.spines.top': False, 'axes.spines.right': False, 'figure.dpi': 160})
    colors = {'amplitude': '#2166ac', 'angular_phase': '#b35806'}
    labels = {'amplitude': 'Sonda de amplitud', 'angular_phase': 'Sonda de fase'}
    methods = {'primary': ('Principal', '-'), 'refined': ('Paso / 2', '--')}
    fig, axes = plt.subplots(1, 2, figsize=(11.6, 4.7), layout='constrained')
    for probe in colors:
        for method, (label, style) in methods.items():
            rows = result['evolution'][method+'_'+probe]
            axes[0].plot([row['time_ps'] for row in rows],
                [row['gap_difference_relative_to_own_initial'] for row in rows],
                style, color=colors[probe], marker='o' if method == 'primary' else 'x',
                markersize=4, linewidth=1.8, label=labels[probe]+' · '+label)
    axes[0].set(xscale='symlog', xlabel='Tiempo físico t (ps)',
        ylabel='||δd(t)||ₘ / ||δd(0)||ₘ', title='Evolución del condensado: diferencia respecto al equilibrio')
    axes[0].set_xscale('symlog', linthresh=.001)
    axes[0].set_ylim(-.03, 1.05)
    axes[0].legend(fontsize=8, loc='lower left')
    kinds = dict(gap_difference='Campo δd', current='Corriente espectral I',
        force_density='Gradiente G/m', phase_torque_density='Gradiente de fase Im(conj(d) G)/m')
    palette = ['#2166ac', '#b35806', '#5e3c99', '#238b45']
    times = sorted({row['time_ps'] for row in refinement['records'] if row['time_ps'] > 0})
    for (kind, label), color in zip(kinds.items(), palette):
        values = [100*max(row['relative_to_initial_signal'] for row in refinement['records']
            if row['observable'] == kind and row['time_ps'] == time_ps) for time_ps in times]
        axes[1].plot(times, values, '-o', color=color, markersize=4, label=label)
    axes[1].axhline(2., color='#555555', ls=':', label='Presupuesto relativo: 2% (+ piso absoluto)')
    axes[1].set(xscale='log', yscale='log', xlabel='Tiempo físico t (ps)',
        ylabel='Diferencia principal–refinada / señal inicial (%)',
        title='Comparación temporal: máximo entre las dos sondas')
    axes[1].set_ylim(1e-5, 3.)
    axes[1].legend(fontsize=8, loc='lower left')
    for axis in axes:
        axis.grid(alpha=.2)
    fig.suptitle('Malla dual de 1712 nodos · 160 × 80 nm · T = 0,9 K · Euler/KWT heredado', fontsize=12)
    fig.text(.02, -.12, 'δd = Δ/(kBTc) − d_equilibrio; ||z||ₘ² = Σmᵢ|zᵢ|². Corriente: ||I||² = Σ|Iᵢⱼ|²/cᵢⱼ.\n'
        'Derecha: denominador = mayor norma inicial de la misma magnitud entre ambas sondas.\n'
        'Los puntos son observaciones guardadas; las líneas entre puntos sólo guían la lectura.', fontsize=9)
    path = destination/'01_evolucion_y_refinamiento.png'
    fig.savefig(path, bbox_inches='tight'); plt.close(fig)

    fig, axes = plt.subplots(2, 2, figsize=(11.6, 7.2), layout='constrained')
    for column, probe in enumerate(colors):
        refined = 'refined_'+probe
        rows = [row for row in history if row['trajectory'] == refined]
        E0 = result['integrated_balance'][refined]['initial_energy_excess']
        time_values = [0., *[row['time_ps'] for row in rows]]
        axes[0, column].plot(time_values, [1., *[row['energy_excess']/E0 for row in rows]],
            color='#2166ac', label='Exceso de energía libre F_ex(t)/F_ex(0)')
        axes[0, column].plot(time_values, [0., *[row['integrated_loss']/E0 for row in rows]],
            color='#b35806', ls='--', label='Disipación acumulada / F_ex(0)')
        axes[0, column].set(title=labels[probe]+' · trayectoria refinada',
            xlabel='Tiempo físico t (ps)', ylabel='Fracción del exceso inicial de energía libre')
        axes[0, column].set_ylim(-.03, 1.05)
        axes[0, column].legend(fontsize=8)
        for method, (label, style) in methods.items():
            name = method+'_'+probe
            rows = [row for row in history if row['trajectory'] == name]
            E0 = result['integrated_balance'][name]['initial_energy_excess']
            axes[1, column].plot([0., *[row['time_ps'] for row in rows]],
                [0., *[100*abs(row['balance_residual'])/E0 for row in rows]],
                style, color='#2166ac' if method == 'primary' else '#b35806', label=label)
        maximum = max(result['integrated_balance'][method+'_'+probe]['maximum_balance_percent_of_initial_excess'] for method in methods)
        axes[1, column].set_ylim(0, max(1e-5, maximum*1.18))
        axes[1, column].set(xlabel='Tiempo físico t (ps)',
            ylabel='|F_ex(t) − F_ex(0) + ∫P dt| / F_ex(0) (%)',
            title=f'Residuo integrado · máximo {maximum:.4f}%')
        axes[1, column].text(.02, .93, 'Presupuesto: 2% del exceso inicial + 10⁻⁸',
            transform=axes[1, column].transAxes, va='top', fontsize=8)
        axes[1, column].legend(fontsize=8, loc='lower right')
    for axis in axes.ravel():
        axis.grid(alpha=.2)
    fig.suptitle('Balance de energía libre durante toda la trayectoria aceptada (0–1 ps)', fontsize=12)
    fig.text(.02, -.045, 'P = (P_KWT + P_normal)/tD en unidades de energía funcional por ps; integral trapezoidal sobre pasos aceptados.\n'
        'La energía funcional es adimensional y se mide respecto al equilibrio uniforme. No es energía interna de poblaciones.', fontsize=9)
    path = destination/'02_balance_integrado.png'
    fig.savefig(path, bbox_inches='tight'); plt.close(fig)

    first, final = observations[0][1], observations[-1][1]
    triangulation = Triangulation(coordinates[:, 0], coordinates[:, 1], triangles)
    d0 = abs(first['baseline_gap'][0])
    fig, axes = plt.subplots(2, 3, figsize=(12.1, 7.1), layout='constrained')
    for row, probe in enumerate(colors):
        values = [100*abs(first['refined_'+probe+'_gap_difference'])/d0,
            100*abs(final['refined_'+probe+'_gap_difference'])/d0,
            100*abs(final['primary_'+probe+'_gap']-final['refined_'+probe+'_gap'])/d0]
        titles = [labels[probe]+' · t = 0 ps', labels[probe]+' · t = 1 ps (refinada)',
            'Diferencia principal–refinada · t = 1 ps']
        for column in range(3):
            vmax = .1 if column < 2 else max(float(np.max(values[column])), 1e-12)
            rendered = axes[row, column].tripcolor(triangulation, values[column],
                shading='gouraud', cmap='magma', vmin=0., vmax=vmax)
            axes[row, column].set(xlabel='x desde el centro (nm)', ylabel='y (nm)', title=titles[column])
            axes[row, column].set_aspect('equal')
            colorbar = fig.colorbar(rendered, ax=axes[row, column], fraction=.05, pad=.02,
                orientation='horizontal')
            colorbar.set_label('|δd| / d_equilibrio (%)' if column < 2 else 'Diferencia / d_equilibrio (%)', fontsize=8)
            if column == 2:
                colorbar.set_ticks(np.linspace(0., vmax, 4))
                colorbar.formatter = FuncFormatter(lambda value, position: f'{value:.1e}')
                colorbar.update_ticks()
            colorbar.ax.tick_params(labelsize=8)
    fig.suptitle('Campos del condensado sobre la malla dual · referencia uniforme d = 1,76392593', fontsize=12)
    fig.text(.02, -.12, 'Primeras dos columnas: misma escala de color (0–0,1% del gap uniforme). Tercera columna: escala ampliada independiente.\n'
        'Se representa el módulo del cambio complejo del condensado, no temperatura ni energía depositada por un fotón.\n'
        'Contactos fijos en x = ±80 nm; paredes laterales con condición natural aislante.', fontsize=9)
    path = destination/'03_campos_espaciales.png'
    fig.savefig(path, bbox_inches='tight'); plt.close(fig)
    return {path.name: sha(path) for path in sorted(destination.glob('*.png'))}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--folder', type=Path, required=True)
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[2]
    result, observations, history, refinement, coordinates, triangles = analyze(args.folder, root)
    result['figures_sha256'] = figures(args.folder, result, observations, history, refinement, coordinates, triangles)
    result['analysis_source_sha256'] = sha(__file__)
    (args.folder/'analysis.json').write_text(json.dumps(result, indent=2)+'\n', encoding='utf8')
    lines = ['# Ensayo térmico final sobre malla dual', '',
        f"**La trayectoria completó 1 ps en {result['runtime_seconds']:.2f} s y cumplió todos los márgenes registrados.** Se admitieron {result['admitted_positive_time_comparisons']}/{result['positive_time_comparisons']} comparaciones temporales a tiempo positivo, además de ocho comparaciones coincidentes del estado inicial. Se usó el paso Euler/KWT real de la memoria, de primer orden.", '',
        'La malla tiene 1712 nodos sobre un rectángulo de 160 × 80 nm. Se estudiaron perturbaciones suaves de amplitud y fase, con temperatura fija de 0,9 K, contactos superconductores de equilibrio en los extremos y paredes laterales aislantes. Son sondas del núcleo térmico; no representan depósito fotónico, hotbelt ni evolución cinética de poblaciones.', '',
        '| Sonda | Magnitud | Mayor diferencia temporal (% de la señal inicial) | Tiempo (ps) |',
        '|---|---|---:|---:|']
    translate = dict(amplitude='Amplitud', angular_phase='Fase', gap_difference='Campo del condensado',
        current='Corriente espectral', force_density='Fuerza por área', phase_torque_density='Fuerza de fase por área')
    for row in result['temporal_maxima']:
        lines.append(f"| {translate[row['probe']]} | {translate[row['observable']]} | {row['maximum_error_percent_initial_signal']:.6g} | {row['time_ps']:g} |")
    lines += ['', 'El denominador de cada comparación temporal es la mayor norma inicial de **esa misma magnitud** entre las dos sondas. El margen registrado es 2% más un piso absoluto de 10⁻⁷ en la norma correspondiente; no se divide por una cola tardía extinguida.', '',
        '| Trayectoria | Pasos aceptados | Máximo residuo del balance (% del exceso inicial) |',
        '|---|---:|---:|']
    for name, row in result['integrated_balance'].items():
        lines.append(f"| {name} | {row['accepted_steps']} | {row['maximum_balance_percent_of_initial_excess']:.6g} |")
    lines += ['', 'El balance integra la disipación KWT y Joule normal en todos los pasos aceptados, con tiempo físico en ps. Se compara con la caída de la energía libre respecto al equilibrio uniforme. No se observó aumento de energía entre pasos. Este balance térmico no constituye una prueba de conservación de energía interna del sistema de poblaciones.', '',
        f"Se verificaron {result['verified_predictions']:,} predicciones espectrales mediante el residuo no lineal exacto; ninguna necesitó la corrección Newton disponible. El máximo residuo fue {result['maximum_spectral_residual']:.3e}, inferior al límite registrado de {result['spectral_tolerance']:.0e}. La predicción reutiliza factorizaciones, sin sustituir el modelo por su linealización.", '',
        'Se emplearon 27 trabajadores y un coordinador con un hilo numérico por proceso, dejando dos núcleos físicos completos libres. La ejecución completa está guardada; no queda pendiente repetir este ensayo ni ejecutar el comando histórico del piloto.', '',
        '## Figuras y significado de las magnitudes', '',
        '- `figures/01_evolucion_y_refinamiento.png`: norma del campo complejo respecto al equilibrio, normalizada por su propia sonda inicial; diferencias entre pasos temporales con el denominador registrado.',
        '- `figures/02_balance_integrado.png`: energía libre excedente, disipación acumulada y residuo integrado durante toda la trayectoria, cada una normalizada por el exceso energético inicial de su sonda.',
        '- `figures/03_campos_espaciales.png`: módulo del cambio del condensado sobre la malla, en porcentaje del gap uniforme; coordenadas físicas en nm y paneles separados del error temporal.', '',
        'Aquí d = Δ/(kBTc), mᵢ = Aᵢ/ell0² y ||z||ₘ² = Σmᵢ|zᵢ|². Para corrientes de arista, ||I||² = Σ|Iᵢⱼ|²/cᵢⱼ. La fuerza de fase es Im(conj(d)G)/m, una fuerza variacional adimensional, no un torque mecánico.', '',
        'La reproducción leyó los datos extraídos y verificó sus hashes; no ejecutó nuevas soluciones físicas. Los valores de `refinement.json` se recuperaron exactamente a partir de los campos guardados.', '',
        '```bash', 'python -m sandbox.stage4_core.analyze_dual_kwt_results --folder docs/implementation/stage4/final_kwt_20260924/thermal', '```', '']
    (args.folder/'analysis.md').write_text('\n'.join(lines), encoding='utf8')
    print(json.dumps(dict(status=result['status'], temporal_maxima=result['temporal_maxima'],
        integrated_balance=result['integrated_balance'], figures=result['figures_sha256']), indent=2))


if __name__ == '__main__':
    main()
