"""Explicitly labelled plots of measured stage4 weak harmonic responses.

Reads the numerical records without assigning acceptance or interpolating a
missing query. Multiple regulators/orders are distinguished and a difference
panel exposes overlapping admittance curves. Figures are scientific PNG/PDF
artifacts; the accompanying captions describe each observable and normalization.
"""
from __future__ import annotations
import argparse
from collections import defaultdict
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator, FuncFormatter
from scipy.constants import hbar, k as KB, e as EC


def decode(value):
    if not isinstance(value, dict) or set(value) != {'real', 'imag'}:
        raise ValueError('Expected explicit real/imag encoded complex values')
    result = np.asarray(value['real'])+1j*np.asarray(value['imag'])
    if np.any(~np.isfinite(result)):
        raise ValueError('Nonfinite record cannot be silently omitted from a plot')
    return result


def configuration(row):
    return row['eta'], row['quadrature_order'], row['energy_cutoff'], row.get('tail_order', 0)


def label_config(config, gap):
    eta, order, cutoff, tail = config
    tail_label = f'E*={cutoff:g} kBTc; cola={tail}' if tail else f'corte={cutoff:g} kBTc'
    return f'η/|Δ|={eta/gap:g}; orden={order}; {tail_label}'


def grouped(records):
    groups = defaultdict(list)
    for row in records:
        groups[configuration(row)].append(row)
    for key in groups:
        groups[key].sort(key=lambda r:r['omega'])
        frequencies = [r['omega'] for r in groups[key]]
        if len(frequencies) != len(set(frequencies)):
            raise ValueError('Duplicate physical query in a plotted configuration')
    return dict(sorted(groups.items()))


def appearance(index):
    colors = plt.get_cmap('tab10').colors
    markers = ('o', 's', '^', 'D', 'v', 'P', 'X', '<', '>')
    lines = ('-', '--', '-.', ':')
    return dict(color=colors[index % len(colors)], marker=markers[index % len(markers)],
        linestyle=lines[(index//len(colors)) % len(lines)], linewidth=1.45,
        markersize=4.5, markerfacecolor='none', alpha=.85)


def finish(fig, axes, output, stem, title, handles, labels):
    for axis in np.ravel(axes):
        axis.set_xscale('log')
        axis.set_xlabel('Frecuencia f = ω/(2π) [GHz]')
        axis.grid(True, which='both', alpha=.23)
    fig.suptitle(title, fontsize=12, y=.995)
    fig.legend(handles, labels, loc='lower center', bbox_to_anchor=(.5, 0.),
               ncol=min(2, len(labels)), fontsize=7.5, frameon=False)
    # Reserve enough vertical room for complete configuration descriptions.
    legend_rows = int(np.ceil(len(labels)/2))
    fig.tight_layout(rect=(0., .055+.020*legend_rows, 1., .97))
    files = []
    for extension in ('png', 'pdf'):
        target = output/(stem+'.'+extension)
        fig.savefig(target, dpi=190, facecolor='white', bbox_inches='tight')
        files.append(str(target))
    plt.close(fig)
    return files


def plot_results(results_path, output, manifest_path=None):
    source = Path(results_path)
    payload = json.loads(source.read_text(encoding='utf8'))
    manifest_path = Path(manifest_path) if manifest_path else source.with_name('manifest.json')
    manifest = json.loads(manifest_path.read_text(encoding='utf8'))
    plan = manifest['plan']
    temperature = plan['Tc_K']
    gap = manifest['d']
    E0 = KB*temperature
    frequency_GHz = lambda rows: np.array([r['omega'] for r in rows])*E0/hbar/(2*np.pi)*1e-9
    records = payload['records']
    if not records:
        raise ValueError('No completed records available to plot')
    output = Path(output); output.mkdir(parents=True, exist_ok=True)
    captions, files = [], []
    plt.rcParams.update({'font.family':'DejaVu Sans', 'axes.titlesize':10,
                         'axes.labelsize':9, 'font.size':9})

    port_rows = [r for r in records if r['kind'] == 'port']
    if port_rows:
        groups = grouped(port_rows)
        reference_key = min(groups, key=lambda x:(x[0], -x[1], -x[2], -x[3]))
        reference = {r['omega']:complex(decode(r['admittance_S'])) for r in groups[reference_key]}
        fig, axes = plt.subplots(2, 2, figsize=(11.6, 7.4))
        handles, labels = [], []
        drives = set()
        differences = []
        for index, (config, rows) in enumerate(groups.items()):
            x = frequency_GHz(rows)
            Y = np.array([complex(decode(r['admittance_S'])) for r in rows])
            style = appearance(index)
            line, = axes[0, 0].plot(x, Y.real*1e6, **style)
            handles.append(line); labels.append(label_config(config, gap))
            axes[0, 1].plot(x, Y.imag*1e3, **style)
            common = [(row['omega'], value) for row, value in zip(rows, Y) if row['omega'] in reference]
            if common:
                fx = np.array([o for o, _ in common])*E0/hbar/(2*np.pi)*1e-9
                relative = [100*abs(v-reference[o])/abs(reference[o]) if abs(reference[o]) else np.nan
                            for o, v in common]
                differences.extend(value for value in relative if np.isfinite(value))
                axes[1, 0].plot(fx, relative, **style)
            if all('bias_voltage_peak_V' in row for row in rows):
                for row in rows:
                    value = row['bias_voltage_peak_V']
                    drives.add(complex(decode(value)) if isinstance(value, dict) else complex(value))
                voltage = np.array([abs(complex(decode(r['readout_peak_V']))) for r in rows])
                axes[1, 1].plot(x, voltage*1e6, **style)
        axes[0, 0].set(ylabel='Re Ydev [µS]', title='Conductancia: componente disipativa, con su signo')
        axes[0, 0].set_yscale('symlog', linthresh=.01)
        axes[0, 1].set(ylabel='Im Ydev [mS]', title='Susceptancia: componente reactiva, con su signo')
        axes[0, 1].set_yscale('symlog', linthresh=.1)
        axes[1, 0].set(ylabel='100 |Y − Yref| / |Yref| [%]', title='Diferencias que una superposición puede ocultar')
        axes[1, 0].set_yscale('symlog', linthresh=.001)
        if differences and max(differences)<.001:
            axes[1, 0].set_yscale('linear')
            axes[1, 0].yaxis.set_major_locator(MaxNLocator(4))
            axes[1, 0].yaxis.set_major_formatter(FuncFormatter(lambda value, _: f'{value:.2g}'))
        axes[1, 0].set_ylim(bottom=0.)
        axes[1, 1].set(ylabel='Amplitud máxima |δVout| [µV]', title='Lectura después del condensador de la memoria')
        if not drives:
            axes[1, 1].text(.5, .5, 'Falta amplitud de fuente en los registros.\nNo se atribuye una escala de lectura.',
                            ha='center', va='center', transform=axes[1, 1].transAxes)
        else:
            axes[1, 1].set_yscale('log')
            if len(drives)==1:
                axes[1, 1].set_title('Lectura después del condensador de la memoria\n'
                    f'|δVbias|={abs(next(iter(drives)))*1e6:.2g} µV; amplitudes máximas')
        files += finish(fig, axes, output, 'port_admittance_and_readout',
            'Cinta uniforme sin fotón: respuesta lineal del puerto y circuito completo', handles, labels)
        captions.append('port_admittance_and_readout: Ydev=δIs/δVdev, con puerto pasivo. '
            'Arriba se grafican sus partes real (µS) e imaginaria (mS), sin valor absoluto. '
            'Abajo izquierda: norma compleja de la diferencia relativa respecto de '+label_config(reference_key, gap)+
            '; se comparan sólo frecuencias comunes, sin interpolar. Abajo derecha: amplitud máxima de δVout '
            'después del condensador, en µV. Fuente de polarización armónica registrada: '+str(sorted(drives, key=abs))+'. '
            'Son incrementos armónicos sobre equilibrio, no un pulso fotónico ni una latencia. '
            'Todas las curvas muestran marcadores huecos y su configuración completa; coincidencia visual '
            'se comprueba en el panel de diferencias. Convención cinética exp(−iωt).')

        fig, axes = plt.subplots(1, 2, figsize=(11.6, 4.3))
        handles, labels = [], []
        for index, (config, rows) in enumerate(groups.items()):
            x = frequency_GHz(rows); style = appearance(index)
            neutral = np.array([abs(complex(decode(r['neutral_gauge_residual']))) for r in rows])
            line, = axes[0].plot(x, neutral, **style)
            handles.append(line); labels.append(label_config(config, gap))
            axes[1].plot(x, [abs(r['circuit_power_residual_W']) for r in rows], **style)
            # The optional correction must have been computed upstream from
            # the physical endpoint integral; the plotting code never fits it.
            if all('neutral_finite_cutoff_target' in r for r in rows):
                error = [abs(complex(decode(r['neutral_gauge_residual']))-
                             complex(decode(r['neutral_finite_cutoff_target']))) for r in rows]
                axes[0].plot(x, error, color=style['color'], linestyle=':', linewidth=1.1)
        all_tails = all(r.get('tail_order', 0)>0 for r in port_rows)
        axes[0].set(ylabel='|1 + ∫ δK₀ dE / 2| [adimensional]',
                    title='Residuo de calibre: integral con colas resueltas' if all_tails else
                          'Residuo de calibre: incluye corte de energía finito')
        axes[0].set_yscale('symlog', linthresh=1e-10)
        axes[1].set(ylabel='|⟨Pfuente − Pdev − PRb − PRL⟩| [W]',
                    title='Balance medio del circuito de tres estados')
        axes[1].set_yscale('symlog', linthresh=1e-35)
        files += finish(fig, axes, output, 'port_balance_diagnostics',
            'Diagnósticos independientes: un residuo no es una decisión automática de aceptación', handles, labels)
        captions.append('port_balance_diagnostics: izquierda, módulo del residuo de neutralidad para '
            'la dirección de calibre de amplitud potencial v=1; '+
            ('las dos colas infinitas se integran con x=E*/|E|. ' if all_tails else 'incluye el corte espectral finito. ')+
            'Si el registro incluye neutral_finite_cutoff_target, líneas punteadas muestran el error '
            'contra esa integral de extremo, sin ajustar constantes. Derecha: módulo del residuo de '
            'potencia media del circuito en W; usa amplitudes máximas, por ello ⟨VI⟩=Re(VI*)/2. '
            'La pequeña identidad circuital no constituye un balance energético completo de la película.')

    for mode in sorted({r['mode'] for r in records if r['kind'] == 'mode'}):
        rows_mode = [r for r in records if r['kind'] == 'mode' and r['mode'] == mode]
        groups = grouped(rows_mode)
        fig, axes = plt.subplots(2, 2, figsize=(11.6, 7.4))
        handles, labels = [], []
        for index, (config, rows) in enumerate(groups.items()):
            x = frequency_GHz(rows); style = appearance(index)
            matrices = [decode(r['kernel']) for r in rows]
            reduced = [m[0, 0]-m[0, 1:]@np.linalg.solve(m[1:, 1:], m[1:, 0]) for m in matrices]
            line, = axes[0, 0].plot(x, np.real(reduced), **style)
            handles.append(line); labels.append(label_config(config, gap))
            axes[0, 1].plot(x, [r['resolved_over_kwt'] for r in rows], **style)
            axes[1, 0].plot(x, [abs(m[2, 1]) for m in matrices], **style)
            response = [decode(r['candidate_response']) for r in rows]
            axes[1, 1].plot(x, [abs(r[2, 1])*E0/EC*1e6 for r in response], **style)
        axes[0, 0].set(ylabel='Re Keff,radial [adimensional]',
                      title='Fuerza radial tras resolver fase y neutralidad')
        axes[0, 1].set(ylabel='γresuelta / γKWT [adimensional]',
                      title='Disipación ya resuelta frente a movilidad heredada')
        axes[0, 1].axhline(1., color='.4', linewidth=.8, linestyle=':')
        axes[0, 1].axhline(0., color='.4', linewidth=.6)
        axes[0, 1].set_yscale('symlog', linthresh=.1)
        axes[1, 0].set(ylabel='|∂Q / ∂dy| [adimensional]',
                      title='Momento de carga inducido por la dirección de fase')
        axes[1, 1].set(ylabel='|δφ| por fuerza hy en fase [µV]',
                      title='Potencial del cierre candidato, con KWT añadido')
        eigenvalue = rows_mode[0]['lambda_value']
        files += finish(fig, axes, output, f'mode_{mode:03d}_coupling',
            f'Modo espacial {mode}; λ={eigenvalue:.6g}: respuesta débil, no transiente no lineal', handles, labels)
        captions.append(f'mode_{mode:03d}_coupling: modo propio completo de la malla, Lψ=λMψ, '
            f'λ={eigenvalue:.12g}, con contactos homogéneos. Arriba izquierda: parte real del '
            'complemento de Schur del núcleo radial al eliminar las dos variables de fase/potencial. '
            'Arriba derecha: γresuelta=−2 Im(Keff)/Ω, dividida por γKWT heredada; el signo se conserva '
            'y la línea 1 señala igualdad de escalas, no aceptación. Abajo izquierda: módulo del '
            'coeficiente que convierte la perturbación cartesiana dy en el momento de carga Q, '
            'tercera fila/segunda columna del núcleo. Abajo derecha: potencial φ=(kBTc/e)v '
            'obtenido del cierre candidato cuando la fuerza de fase tiene amplitud hy; se expresa '
            'en µV por esa fuerza de prueba, no como voltaje de un fotón. El candidato añade KWT '
            'a la respuesta resuelta; este gráfico mide la posible superposición y no la valida. '
            'No hay normalización por una cola temporal extinguida.')

    summary = dict(source=str(source), manifest=str(manifest_path), results_status=payload['status'],
        completed_records=len(records), failed_records=len(payload.get('failures', [])),
        outputs=files, captions=captions, acceptance_decision='NONE_ASSIGNED_BY_PLOTTER')
    (output/'plot_manifest.json').write_text(json.dumps(summary, indent=2, ensure_ascii=False)+'\n', encoding='utf8')
    text = '# Qué muestra cada figura\n\n'+f"Estado de los datos: {payload['status']}; {len(records)} consultas guardadas, "+\
           f"{len(payload.get('failures', []))} fallidas. Ninguna figura asigna aceptación física o numérica.\n\n"+\
           '\n\n'.join(captions)+'\n'
    (output/'CAPTIONS.md').write_text(text, encoding='utf8')
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--results', required=True, type=Path)
    parser.add_argument('--output', required=True, type=Path)
    parser.add_argument('--manifest', type=Path)
    args = parser.parse_args()
    summary = plot_results(args.results, args.output, args.manifest)
    print(json.dumps({'figures':len(summary['outputs'])//2,
                      'failed_queries':summary['failed_records'], 'output':str(args.output)}))


if __name__ == '__main__':
    main()
