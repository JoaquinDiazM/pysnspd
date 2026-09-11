"""Prepare the reviewed stage-1 verdict and result figures from measured data."""
from __future__ import annotations
import hashlib
import json
from pathlib import Path
import re

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import NullLocator
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT/'docs/implementation/stage1'


def read(path):
    return json.loads((DATA/path).read_text(encoding='utf-8'))


def main():
    material = read('material_results.json')
    electronic = read('electronic/electronic_diagnostics.json')
    benchmark = read('query_benchmark.json')
    boundary = read('electronic/boundary_comparison.json')
    regulator = read('review/low_energy_regulator_reference.json')
    review = read('review/verify_query_contract.json')
    tests = (DATA/'pytest_geminga.txt').read_text(encoding='utf-8')
    match = re.search(r'(\d+) passed in ([\d.]+)s', tests)
    if not match:
        raise RuntimeError('A successful final regression log is required.')
    passed, seconds = int(match[1]), float(match[2])
    case = next(r for r in boundary['rows'] if r['delta']==.72 and r['gamma']==0)
    bias = next(r['relative_regulator_bias'] for r in regulator['rows'] if r['eta_over_Delta0']==.001)
    assert material['nbn_audit']['status']=='REJECTED'
    assert bias > .08 and case['new_relative_error'] < 1e-5
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,
                         'axes.spines.top':False,'axes.spines.right':False})
    fig, axes = plt.subplots(1,2,figsize=(10.2,3.5),layout='constrained')
    old = np.loadtxt(DATA/'history/before_gamma_hermite/low_energy_grid_refinement.csv',delimiter=',',skiprows=1)
    new = np.loadtxt(DATA/'electronic/low_energy_grid_refinement.csv',delimiter=',',skiprows=1)
    axes[0].loglog(old[:,0],old[:,3],'o--',color='#ad5c45',label='Interpolación anterior')
    axes[0].loglog(new[:,0],new[:,3],'o-',color='#007c83',label='Hermite: pendientes espectrales')
    axes[0].set(xlabel='Nodos por eje de campo',ylabel='Error absoluto normalizado',
                title='Respuesta al superflujo: 5 campos')
    axes[0].set_xticks([9,17,33],['9','17','33'])
    axes[0].xaxis.set_minor_locator(NullLocator())
    axes[0].legend(fontsize=8)
    eta = np.array([r['eta_over_Delta0'] for r in regulator['rows']])
    biases = 100*np.array([r['relative_regulator_bias'] for r in regulator['rows']])
    axes[1].loglog(eta,biases,'o-',color='#16334b',label='Integral adaptativa / límite causal')
    axes[1].scatter([.001],[100*bias],s=65,color='#ad5c45',zorder=5,label='Piloto actual: 8.61%')
    axes[1].set(xlabel='Regulador η / Δ₀',ylabel='Sesgo relativo en la pendiente j/q (%)',
                title='Error que permanece sin interpolación')
    axes[1].legend(fontsize=8)
    (DATA/'figures').mkdir(exist_ok=True)
    for suffix in ('png','pdf'):
        fig.savefig(DATA/f'figures/catalogue_decision.{suffix}',dpi=180)
    plt.close(fig)

    admission = {
        'schema':'pysnspd.stage1.catalog-admission.v1',
        'release_tag':'v1.0.0',
        'release_commit':'5ea0cd6a08d2cae414c82944f8230469ca5aa7d6',
        'implementation_status':'implemented_opt_in_experimental',
        'stage1_acceptance_status':'NOT_CLOSED',
        'promotion_to_stage2':'NOT_ADMITTED',
        'vacuum_status':'sampled_uniform_references_verified',
        'occupation_catalog_status':'diagnostic_pilot_not_admitted_for_coupled_dynamics',
        'material_status':'REJECTED',
        'regulator_case_relative_bias':bias,
        'regulator_reference':'review/low_energy_regulator_reference.json',
        'interpolation_case_relative_error_after_fix':case['new_relative_error'],
        'interpolation_reference':'electronic/boundary_comparison.json',
        'pytest':{'passed':passed,'new_tests':passed-175,'seconds':seconds,'host':'geminga'},
        'catalog_build_seconds':electronic['runtime_seconds'],
        'files_sha256':{p.name:hashlib.sha256(p.read_bytes()).hexdigest()
                        for p in sorted((DATA/'catalogs').glob('*.npz'))},
        'next_required_work':[
            'Converge the causal regulator with quadrature concentrated at the low-energy edge and a Gamma mesh resolving its narrow scale.',
            'Repeat fixed-p force/current and independent thermal references using those resolved tables before coupling cells.',
            'Obtain single-valued, positive phonon data with verified ordinate units, complete basis and matching density.',
        ],
        'physical_interpretation':'The regulator discrepancy is a numerical approximation error; this test does not falsify the underlying uniform functional.',
        'production_model_updated':False,'full_transient_executed':False,
        'long_calculation_running_or_required_for_this_report':False,
    }
    (DATA/'catalog_admission.json').write_text(json.dumps(admission,indent=2,ensure_ascii=False)+'\n',encoding='utf-8')

    def b(text,kind='body'): return {'type':kind,'text':text}
    def table(rows,widths): return {'type':'table','rows':rows,'widths':widths}
    def picture(path,caption): return {'type':'image','path':path,'caption':caption}
    query = benchmark['catalogs']['occupation']['queries']['median_ms_per_query']
    vac_query = benchmark['catalogs']['vacuum']['queries']['median_ms_per_query']
    docs_url = 'https://github.com/JoaquinDiazM/pysnspd/tree/v1.0.0'
    source_url = material['nbn_audit']['metadata']['source_url']
    pages = [
        {'title':'Etapa 1: datos y catálogo', 'blocks':[
            b('<b>Resultado: implementación experimental terminada; admisión de la etapa todavía pendiente.</b> Recomiendo corregir la resolución espectral antes de pasar al acoplamiento de una y dos celdas.'),
            b(f'La <link href="{docs_url}" color="#007c83">versión 1.0.0</link> quedó publicada primero, en el commit 5ea0cd6, con el solver de la memoria y la documentación 0.4 / E-r02. El código de esta etapa se incorpora después del tag y se activa sólo por importación explícita.'),
            table([
                ['Verificación','Resultado','Decisión'],
                ['Regresión del paquete',f'{passed} pruebas aprobadas; {passed-175} nuevas. {seconds:.2f} s en Geminga.','Implementación verificable'],
                ['Vacío y espectro uniforme','24 referencias independientes; 10 836 estados espectrales. Residuo espectral 2.63 × 10⁻¹⁴.','Referencias muestreadas correctas'],
                ['Ocupaciones no térmicas','Energía, fuerzas y corriente proceden de la misma tabla. Conserva orden de energías y rechaza consultas inválidas.','Piloto numérico'],
                ['Respuesta a corriente pequeña',f'Sesgo del regulador: {100*bias:.2f}% en el caso crítico, aun con cuadratura resuelta.','No promover a etapa 2'],
                ['Datos fonónicos NbN','Unidades/base sin certificar; DOS negativa y frecuencias repetidas.','Material rechazado'],
            ],[.26,.48,.26]),
            b('Costo medido en Geminga','subheading'),
            b(f'Construcción y diagnósticos completos: <b>{electronic["runtime_seconds"]:.2f} s</b>, un hilo. Consulta de energía y sus dos derivadas: <b>{query:.3f} ms</b> con población y <b>{vac_query:.3f} ms</b> en el vacío; mediana de tres lotes de 1000 consultas. El archivo del catálogo con población ocupa 1.03 MB.'),
            b('El benchmark mide consultas locales. El transporte, las colisiones, la malla, el circuito y los pasos temporales no están incluidos en ese costo. No se ejecutó un transitorio completo ni quedó un cálculo largo en marcha.','small'),
            b('Apreciación para continuar','subheading'),
            b('El próximo trabajo debe concentrarse en la cuadratura cerca del borde espectral y en la convergencia conjunta de sus mallas y del regulador. Este resultado no exige reformular ahora el funcional uniforme; sí exige resolver su aproximación numérica antes del acoplamiento. La entrada fonónica necesita una fuente material admisible.'),
        ]},
        {'title':'Datos: rechazo con causas identificadas', 'blocks':[
            b('Se recuperó la procedencia exacta del contenido numérico de NbN: al retirar la primera línea añadida localmente, sus bytes coinciden con el archivo público. Los seis archivos usados en la comparación de procedencia pasaron sus hashes. Esa coincidencia no certifica las unidades ni la base de normalización.'),
            picture('docs/implementation/stage1/figures/material_admission.png',
                    'Datos originales, sin recorte ni renormalización. La referencia de tres modos presupone la base por átomo; esa base sigue sin certificar para este archivo.'),
            table([
                ['Hallazgo en NbN','Resultado observado'],
                ['DOS negativa','361 muestras; peso negativo / positivo = 1.601511%.'],
                ['Eje no unívoco','Cuatro frecuencias repetidas con ordenadas distintas.'],
                ['Normalización bajo cabecera local','Integral 0.708420; no coincide con tres modos por átomo. Las unidades de la cabecera no provienen del archivo público.'],
                ['Soporte de interacción','Hay peso de interacción donde la DOS no es positiva y un valor no nulo de α²F en energía cero.'],
            ],[.35,.65]),
            b('Admisión comprobada con una referencia sintética','subheading'),
            b('Un espectro Debye generado analíticamente se admite como dato sintético. Sus cinco representaciones —THz, meV, eV, J y rad/s— conservan la energía hasta precisión numérica. La conversión por átomo / por celda mantiene la densidad volumétrica. Con 1001 nodos, el error energético frente a integración adaptativa es <b>9.36 × 10⁻⁹</b>.'),
            b('Decisión: mantener NbN en cuarentena. Se requiere una tabla positiva y unívoca, unidades de ambas columnas espectrales y base/densidad documentadas. Multiplicar por una constante, promediar duplicados o recortar valores no resolvería esa trazabilidad.'),
            b(f'Fuente y licencia: <link href="{source_url}" color="#007c83">datos públicos de Simon et al., revisión 5b6bd747</link>. Los hashes, rechazos y conversiones están en material_results.json; la licencia de los datos derivados se conserva en MATERIAL_SOURCE_NOTICE.md.','small'),
        ]},
        {'title':'Catálogo: corrección y límite pendiente', 'blocks':[
            b(f'En el caso |Δ|/Δ₀ = 0.72 y Γ = 0, la interpolación Hermite redujo el error de la respuesta al superflujo de <b>{100*case["previous_relative_error"]:.2f}% a {100*case["new_relative_error"]:.7f}%</b>, frente al espectro directo con los mismos 64 nodos y el mismo regulador. Aquí se compara la pendiente j/q al aproximarse a corriente nula; en q = 0 la corriente vale cero.'),
            picture('docs/implementation/stage1/figures/catalogue_decision.png',
                    'Izquierda: máximo error de interpolación en cinco campos con población de baja energía. Derecha: sesgo del regulador calculado por una referencia independiente, sin interpolación y con cuadratura adaptativa.'),
            table([
                ['Comprobación independiente','Resultado'],
                ['Vacío: referencia sobre frecuencias imaginarias','Error máximo escalado 4.15 × 10⁻⁹ en 24 estados.'],
                ['Derivadas de la tabla con población','Errores de energía/fuerza de orden 10⁻¹¹; reciprocidad de derivadas cruzadas de orden 10⁻⁹.'],
                ['Baja energía: sesgo del regulador','η/Δ₀ = 10⁻³: 8.61%; 10⁻⁴: 2.99%; 10⁻⁵: 0.986%.'],
                ['Cuadratura separada, η/Δ₀ = 10⁻³','La respuesta pasa de 0.736723 a 0.737013 al refinar de 64 a 128 nodos. Límite causal: 0.678584.'],
            ],[.42,.58]),
            b('Alcance que queda documentado','subheading'),
            b('El piloto con población cubre 0.08 ≤ |Δ|/Δ₀ ≤ 1.5 y 0 ≤ Γ/Δ₀ ≤ 1.2, con 33 × 33 campos y 64 nodos de estados. El punto normal se trata por separado; amplitudes positivas menores de 0.08 se rechazan. Dentro de los cinco campos de baja energía, persisten errores relativos de interpolación de hasta 0.35% fuera de los nodos.'),
            b('<b>No cerrar todavía la etapa 1.</b> Reducir el regulador sobre la misma malla de 64 nodos puede empeorar la integral de corriente. El siguiente ensayo debe concentrar nodos cerca del borde y resolver la estrecha escala de Γ, repetir estos contrastes y después decidir la promoción. Los ensayos sintéticos de contabilidad que se adelanten no equivalen a admitir el catálogo.'),
            b('Reproducción: sandbox/stage1_catalog y docs/implementation/stage1/review. Comandos y tiempos en /home/jdiaz/GEMINGA_COMMANDS.md. Los archivos de esta entrega se verifican por SHA-256.','small'),
        ]},
    ]
    content = {'title':'pySNSPD · Informe de resultados de la etapa 1',
               'subtitle':'11 de septiembre de 2026 · Modelo documentado 0.4 · Implementación posterior a v1.0.0',
               'pages':pages}
    (DATA/'report_content.json').write_text(json.dumps(content,indent=2,ensure_ascii=False)+'\n',encoding='utf-8')
    print(json.dumps({'pytest':passed,'promotion':'NOT_ADMITTED','regulator_bias_percent':100*bias}))


if __name__=='__main__':
    main()
