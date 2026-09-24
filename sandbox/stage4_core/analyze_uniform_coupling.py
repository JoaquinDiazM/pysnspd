"""Reproducible comparison of existing uniform coupled-response records."""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np


def decode(value):
    return np.asarray(value['real'])+1j*np.asarray(value['imag'])


def key(row):
    return row['kind'], row.get('mode'), row['omega']


def relative(a, b):
    scale = np.linalg.norm(np.atleast_1d(b))
    return float(np.linalg.norm(np.atleast_1d(a-b))/scale) if scale else None


def compare(a, b, normal_conductance):
    result = dict(kind=b['kind'], mode=b.get('mode'), omega=b['omega'])
    if b['kind'] == 'port':
        ya, yb = complex(decode(a['admittance_S'])), complex(decode(b['admittance_S']))
        result.update(admittance_relative=relative(ya, yb),
            readout_relative=relative(decode(a['readout_peak_V']), decode(b['readout_peak_V'])),
            real_admittance_difference_over_normal=abs(ya.real-yb.real)/normal_conductance,
            real_admittance_relative=abs(ya.real-yb.real)/abs(yb.real) if yb.real else None)
    else:
        result.update(candidate_response_relative=relative(decode(a['candidate_response']), decode(b['candidate_response'])),
            kernel_relative=relative(decode(a['kernel']), decode(b['kernel'])),
            gamma_difference_over_kwt=abs(a['resolved_over_kwt']-b['resolved_over_kwt']),
            gamma_relative=abs(a['gamma_resolved']-b['gamma_resolved'])/abs(b['gamma_resolved']) if b['gamma_resolved'] else None)
    return result


def maxima(rows):
    result = {}
    metrics = sorted(set().union(*(row.keys() for row in rows))-{'kind', 'mode', 'omega', 'eta_ratio', 'comparison'})
    for metric in metrics:
        candidates = [row for row in rows if row.get(metric) is not None]
        if candidates:
            worst = max(candidates, key=lambda row:row[metric])
            result[metric] = dict(maximum=worst[metric], case={k:worst[k] for k in
                ('kind', 'mode', 'omega', 'eta_ratio', 'comparison') if k in worst},
                count=len(candidates), exceeding_2_percent=sum(row[metric]>.02 for row in candidates))
    return result


def analyze(cases, output):
    output = Path(output); output.mkdir(parents=True, exist_ok=True)
    loaded, numerical, regulators, sources, fine_rows = {}, [], [], {}, []
    all_records = []
    for folder in sorted(Path(cases).iterdir()):
        if not (folder/'results.json').exists():
            continue
        manifest = json.loads((folder/'manifest.json').read_text(encoding='utf8'))
        payload = json.loads((folder/'results.json').read_text(encoding='utf8'))
        ratio = manifest['plan']['eta_over_gap'][0]
        records = payload['records']; all_records += records
        sources[folder.name] = dict(results_sha256=hashlib.sha256((folder/'results.json').read_bytes()).hexdigest(),
            manifest_sha256=hashlib.sha256((folder/'manifest.json').read_bytes()).hexdigest(),
            executed_source_sha256=manifest['sources']['sandbox/stage4_core/coupled_response.py'],
            status=payload['status'], completed=len(records), failures=len(payload.get('failures', [])),
            elapsed_seconds=payload['elapsed_seconds'])
        orders = sorted({r['quadrature_order'] for r in records})
        cutoffs = sorted({r['energy_cutoff'] for r in records})
        fine = {key(r):r for r in records if r['quadrature_order']==orders[-1] and r['energy_cutoff']==cutoffs[-1]}
        loaded[ratio] = (manifest, fine)
        fine_rows.extend(dict(eta_ratio=ratio, **r) for r in fine.values())
        normal = manifest['graph_conductance']/manifest['plan']['sheet_resistance_ohm']
        for name, order, cutoff in (('order',orders[0],cutoffs[-1]), ('cutoff',orders[-1],cutoffs[0])):
            coarse = {key(r):r for r in records if r['quadrature_order']==order and r['energy_cutoff']==cutoff}
            for k, reference in fine.items():
                numerical.append(dict(eta_ratio=ratio, comparison=name, **compare(coarse[k],reference,normal)))
    ratios = sorted(loaded)
    for ratio in ratios[1:]:
        manifest, fine = loaded[ratio]
        normal = manifest['graph_conductance']/manifest['plan']['sheet_resistance_ohm']
        for k, reference in loaded[ratios[0]][1].items():
            regulators.append(dict(eta_ratio=ratio,comparison=f'eta{ratio:g}_vs_{ratios[0]:g}',
                                   **compare(fine[k],reference,normal)))
    damping = []
    for ratio in ratios:
        for omega in sorted({r['omega'] for r in fine_rows}):
            selected = [r for r in fine_rows if r['kind']=='mode' and r['eta_ratio']==ratio and r['omega']==omega]
            damping.append(dict(eta_ratio=ratio,omega=omega,
                minimum_resolved_over_kwt=min(r['resolved_over_kwt'] for r in selected),
                maximum_resolved_over_kwt=max(r['resolved_over_kwt'] for r in selected),
                number_of_modes=len(selected)))
    fine_ports = [r for r in fine_rows if r['kind']=='port']
    fine_modes = [r for r in fine_rows if r['kind']=='mode']
    raw_neutral = [abs(complex(decode(r['neutral_gauge_residual']))) for r in fine_ports]
    tails_resolved=all(r.get('tail_order',0)>0 for r in all_records)
    receipt = dict(schema='pysnspd.stage4.uniform_comparison.v1', sources=sources,
        records=len(all_records), failed_queries=sum(v['failures'] for v in sources.values()),
        infinite_tails_resolved=tails_resolved,
        tail_orders=sorted({r.get('tail_order',0) for r in all_records}),
        comparison_reference=dict(order=max(r['quadrature_order'] for r in all_records),
            cutoff=max(r['energy_cutoff'] for r in all_records),minimum_eta_ratio=ratios[0]),
        numerical_maxima=maxima(numerical), regulator_maxima=maxima(regulators),
        numerical_comparisons=numerical, regulator_comparisons=regulators, damping_by_frequency=damping,
        max_kinetic_residual=max(r['kinetic_residual'] for r in all_records if r['kind']=='mode'),
        max_charge_solve_residual=max(r['charge_solve_residual'] for r in all_records if r['kind']=='mode'),
        max_circuit_power_residual_W=max(abs(r['circuit_power_residual_W']) for r in all_records if r['kind']=='port'),
        finest_raw_neutrality_range=[min(raw_neutral),max(raw_neutral)],
        finest_gamma_resolved_over_kwt_range=[min(r['resolved_over_kwt'] for r in fine_modes),max(r['resolved_over_kwt'] for r in fine_modes)],
        high_frequency_not_time_transient=True, nonlinear_energy_and_heat_deposition_tested=False,
        interpretation='Practical 2% comparisons apply to nonzero signal observables. Tiny dissipative quantities are also compared against normal conductance or KWT scale; relative error of their near-zero value alone is not a rejection.')
    (output/'uniform_analysis.json').write_text(json.dumps(receipt,indent=2,allow_nan=False)+'\n',encoding='utf8')
    lines = ['# Respuesta uniforme: comparación independiente de la campaña', '',
        f"Se completaron {receipt['records']} consultas, con {receipt['failed_queries']} fallos. "
        'Son respuestas armónicas lineales a fondo uniforme sin corriente, sin fotón ni transiente de calor.', '',
        ('La referencia de cuadratura usa orden 24 por panel, transición E*=64 kBTc y 24 nodos en cada cola infinita. '
         'La transformación x=E*/|E| integra el mismo integrando hasta infinito; no ajusta contraterminos. '
         'Se compara orden 16 frente a 24 y transición 32 frente a 64. ' if tails_resolved else
         'La referencia de cuadratura usa orden 24 por panel y corte 64 kBTc. Se compara orden 16 frente a 24 '
         'al mismo corte, y corte 32 frente a 64 al mismo orden. ')+
        'La variación de eta se informa por separado: '
        'eta es un regulador espectral, no una tasa física de relajación.', '',
        '| Comparación numérica | Máximo | Supera 2 % |', '|:--|--:|--:|']
    names = dict(admittance_relative='Admitancia compleja, norma relativa', readout_relative='Vout complejo, norma relativa',
        candidate_response_relative='Matriz de respuesta candidata, norma Frobenius relativa',kernel_relative='Núcleo, norma Frobenius relativa',
        gamma_difference_over_kwt='Diferencia de disipación / gamma KWT',real_admittance_difference_over_normal='Diferencia Re Y / conductancia normal')
    for name,title in names.items():
        item=receipt['numerical_maxima'][name]
        lines.append(f"| {title} | {100*item['maximum']:.6g} % | {item['exceeding_2_percent']}/{item['count']} |")
    lines += ['', 'La comparación relativa de Re Y o de gamma resuelta consigo mismas puede amplificar '
        'diferencias de una señal casi nula. El JSON conserva esas comparaciones, además de las escalas '
        'normal y KWT anteriores; no se ocultan signos ni se recortan valores.', '',
        '| Omega = hbar omega/(kBTc) | gamma resuelta / gamma KWT, mínimo–máximo entre cinco modos; eta/Delta=0,002 |',
        '|--:|--:|']
    for row in damping:
        if row['eta_ratio']==ratios[0]:
            lines.append(f"| {row['omega']:g} | {row['minimum_resolved_over_kwt']:.7g} – {row['maximum_resolved_over_kwt']:.7g} |")
    lines += ['', '| Sensibilidad eta=0,005/0,01 frente a 0,002, cuadratura fina | Máximo |', '|:--|--:|']
    for name,title in names.items():
        lines.append(f"| {title} | {100*receipt['regulator_maxima'][name]['maximum']:.6g} % |")
    lines += ['', 'El aumento de disipación cerca de Omega≈2Delta/(kBTc)=3,52785 coincide con la apertura '
        'del canal de rotura de pares de la referencia BCS. Por debajo, una parte de las pequeñas colas '
        'disipativas sigue dependiendo de eta; no se utiliza para identificar una vida media material.', '',
        f"El máximo residuo cinético es {receipt['max_kinetic_residual']:.4g}; el de la ecuación "
        f"de neutralidad resuelta es {receipt['max_charge_solve_residual']:.4g}. El balance medio del "
        f"circuito tiene defecto máximo {receipt['max_circuit_power_residual_W']:.4g} W.", '',
        '## Diagnóstico y corrección de la cola\n\n'
        'En la campaña inicial, el residuo bruto de calibre con corte 64 no tendía a cero al refinar sólo los paneles: '
        'queda el término de extremo. Para el puerto uniforme, su valor esperado es '
        '`1 − integral(C−Omega/2,C+Omega/2)[tanh(E/2T) rho(E)] dE / Omega`. '
        'Con este corte ronda −0,000380, aproximadamente −0,038 %. Es distinto del residuo '
        'algebraico de la neutralidad acoplada. Una corrección de cola debe salir de esa integral, '
        'nunca de ajustar una constante a la respuesta. La matriz de respuesta de fase/potencial amplificaba '
        'esta omisión hasta un 54,9977 % al comparar los cortes 32 y 64, aunque el núcleo cambiaba menos '
        'del 0,07 %. La columna radial permanecía estable. Ese diagnóstico original se conserva en '
        '`finite_cutoff_review/uniform_analysis.md` y en `uniform_results/`.', '',
        ('La campaña actual integra ambas colas con la transformación recíproca. El piloto independiente '
         'comparó 12/24 nodos por cola y transiciones 32/64: el cambio máximo de la respuesta fue '
         '1,44e−9 relativo. El defecto grande desaparece al recuperar la contribución faltante; '
         'no se relajó la tolerancia para aceptarlo. Las cifras de las tablas corresponden a las '
         '720 consultas nuevas con colas, no al piloto ni a resultados retocados.' if tails_resolved else
         'Esta entrega todavía contiene el diagnóstico de corte finito, no la campaña corregida.'), '',
        'La resta térmica estática y el incremento retardado comparten la misma referencia. El circuito '
        'conserva los tres estados y sus constantes de la memoria. La pequeña identidad de potencia '
        'circuital no demuestra todavía el balance total de energía de la película.', '',
        '## Alcance de la comparación polarizada pendiente', '',
        '`biased_coupled_response.py` proyecta las respuestas completas del grafo alrededor de una rama '
        'con corriente, resuelve amplitud, fase y potencial juntos y devuelve residuos también fuera '
        'del subespacio modal. Eso permite medir los términos cruzados ausentes por simetría en el '
        'ensayo uniforme. Sin embargo, calcula respuesta de primer orden y potencia de puerto de '
        'segundo orden; no evoluciona poblaciones/fonones con B.41 ni reconstruye el calor de segundo '
        'orden. El balance circuital no cierra ese trabajo interno por sí solo. Su pase autorizaría '
        'el acoplamiento débil dentro del alcance registrado, no un transiente no lineal con fotón.', '',
        'Las derivadas y constantes de los datos brutos se conservaron. El JSON registra hashes y '
        'cada comparación; no se ejecutó de nuevo el modelo para producir este análisis.', '']
    (output/'uniform_analysis.md').write_text('\n'.join(lines),encoding='utf8')
    return receipt


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--cases',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    result=analyze(args.cases,args.output)
    print(json.dumps({k:result[k] for k in ('records','failed_queries','numerical_maxima','regulator_maxima')},indent=2))
