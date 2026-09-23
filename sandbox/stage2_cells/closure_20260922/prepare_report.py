"""Final development-stage report from saved trajectories; zero RHS calls.

Run with numpy, matplotlib and reportlab. The original result files and criteria
are read-only. Regeneration is lightweight and does not rerun the acceptance batch.
"""
from pathlib import Path
import hashlib
import json
import shutil
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter

ROOT = Path(__file__).resolve().parents[3]
DATA = ROOT/'docs/implementation/stage2/closure_20260922'
PREVIOUS = ROOT/'docs/implementation/stage2/practical_review_20260922'
RAW = PREVIOUS/'raw'
FIG = DATA/'figures'
FIG.mkdir(parents=True, exist_ok=True)
read = lambda p: json.loads(p.read_text(encoding='utf-8'))
sha = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()
audit = read(PREVIOUS/'saved_results_audit.json')
for name, digest in audit['raw_file_sha256'].items():
    assert sha(RAW/name) == digest, name
assert audit['trajectory_count'] == 13 and audit['completed_task_count'] == 21
runs = {case: {n: read(RAW/f'{case}_guarded_{n}.json')
               for n in (160, 320, 640, 1280)} for case in ('one', 'two')}
assessments = {case: read(RAW/f'{case}_guarded_time_assessment.json') for case in runs}
plt.rcParams.update({'font.family':'DejaVu Sans', 'font.size':10,
    'axes.titlesize':11, 'axes.spines.top':False, 'axes.spines.right':False,
    'savefig.dpi':200})
BLUE, TEAL, ORANGE = '#3274b2', '#007f73', '#cc6e13'
ref, coarse = runs['two'][1280], runs['two'][320]
t, tc = np.asarray(ref['times']), np.asarray(coarse['times'])
series = lambda run, name: np.array([s[name] for s in run['snapshots']])
assert np.allclose(t[::4], tc, atol=1e-14, rtol=0)
fig, ax = plt.subplots(2, 2, figsize=(10, 6.0), layout='constrained')
for i, color in enumerate((ORANGE, BLUE)):
    ax[0,0].plot(t, series(ref,'amplitudes')[:,i], color=color, lw=2,
        label=f'Celda {i+1}: 1280 pasos')
    ax[0,0].plot(tc[::20], series(coarse,'amplitudes')[::20,i], ls='none',
        marker='o' if i==0 else 's', ms=4, mfc='white', mec=color,
        label=f'Celda {i+1}: 320 pasos')
    ax[1,0].plot(tc, 1e8*(series(coarse,'amplitudes')[:,i]-series(ref,'amplitudes')[::4,i]),
        color=color, lw=1.6, label=f'Celda {i+1}')
    ax[0,1].plot(t, series(ref,'excitation_energy')[:,i], color=color, lw=2,
        label=f'Electrones, celda {i+1}')
    ax[0,1].plot(t, series(ref,'phonon_energy')[:,i], color=color, lw=1.6, ls='--',
        label=f'Fonones, celda {i+1}')
ax[0,0].set(title='Respuesta del condensado', ylabel=r'$|\Delta|/\Delta_0$')
ax[1,0].set(title='Lo que oculta la superposición', ylabel=r'Amplitud: (320 - 1280) $\times 10^8$')
ax[1,0].axhline(0, color='#888888', lw=.7)
ax[0,1].set(title='La excitación se redistribuye', ylabel=r'Energía / $(N_0\Delta_0^2)$')
for name, label, color in [('transport','Entre celdas (neto)',TEAL),
                          ('input','Fuente externa',ORANGE),('escape','Escape al baño',BLUE)]:
    values=series(ref,name)
    if values.ndim>1: values=values.sum(axis=1)
    ax[1,1].plot(t, values, color=color, lw=2, label=label)
ax[1,1].set(title='Transferencias acumuladas', ylabel=r'Energía / $(N_0\Delta_0^2)$')
for a in ax.flat:
    a.set_xlabel(r'Tiempo $t/t_{\mathrm{ref}}$'); a.grid(alpha=.2); a.legend(fontsize=8)
fig.savefig(FIG/'dynamics_guarded.png'); plt.close(fig)

steps=np.array([160,320,640,1280])
fig,ax=plt.subplots(1,2,figsize=(10,3.8),layout='constrained')
for case,label,color,marker in [('one','Una celda',BLUE,'o'),('two','Dos celdas',TEAL,'s')]:
    assessment=assessments[case]
    assert assessment['status']=='PASS'
    errors=[row['max_error'] for row in assessment['cases']]
    ax[0].plot(steps[:3],100*np.array(errors),marker+'-',color=color,lw=2,label=label)
    ax[1].plot(steps,[runs[case][int(n)]['energy_ledger_scaled_max'] for n in steps],
        marker+'-',color=color,lw=2,label=label)
ax[0].axhline(.01,ls='--',color='#aa4444',label='Criterio registrado: 0,01 %')
ax[1].axhline(1e-7,ls='--',color='#aa4444',label='Criterio registrado: 1e-7')
ax[0].set(title='Diferencia frente a 1280 pasos',ylabel='Máximo relativo registrado (%)',ylim=(.00018,.02))
ax[1].set(title='Balance energético con fuente y baño',ylabel='Máximo defecto escalado',ylim=(1e-11,3e-7))
for i,a in enumerate(ax):
    a.set(xscale='log',yscale='log',xlabel='Número de pasos')
    a.set_xticks(steps[:3] if i==0 else steps,steps[:3] if i==0 else steps)
    a.xaxis.set_minor_formatter(NullFormatter());a.grid(alpha=.2,which='both');a.legend(fontsize=8)
fig.savefig(FIG/'time_and_energy.png');plt.close(fig)

copied={
 'precision_context.png':PREVIOUS/'figures/precision_context.png',
 'circuit_comparison.png':PREVIOUS/'figures/circuit_comparison.png',
 'guarded_rounding_and_escape.png':ROOT/'docs/implementation/stage2/closure_prep_20260922/figures/guarded_rounding_and_escape.png',
}
for name, source in copied.items(): shutil.copy2(source,FIG/name)

def text(value,kind='body'):return dict(type=kind,text=value)
def image(name,caption):return dict(type='image',path=(FIG/name).relative_to(ROOT).as_posix(),caption=caption)
def table(rows,widths):return dict(type='table',rows=rows,widths=widths)
pages=[dict(title='Etapa 2 cerrada: celdas acopladas',blocks=[
 text('<b>Decisión: cerrar esta etapa de desarrollo y avanzar al ensayo espacial con bordes y circuito.</b> La evidencia obtenida permite conservar los operadores cinéticos y la estrategia numérica, sin reformular ahora el modelo físico. El cierre tiene el alcance de los ensayos sintéticos de una y dos celdas.'),
 table([['Resultado de cierre','Evidencia disponible'],
 ['Último lote ejecutado','21 tareas completas; 13 trayectorias válidas; 34 archivos fuente verificados'],
 ['Ocupaciones físicas','Las 13 trayectorias mantienen electrones entre 0 y 1 y fonones no negativos'],
 ['Balance energético','Máximo defecto escalado: 5,44e-8; criterio registrado: 1e-7'],
 ['Consistencia instantánea','Residuo máximo registrado: 1,64e-14 en los tiempos muestreados'],
 ['Refinamiento temporal propio','Una y dos celdas pasan en la malla candidata de 630/1025'],
 ['Ensayo fonónico más fino','Trayectorias válidas; comparación auxiliar conserva su FAIL'],
 ['Paso siguiente','Desarrollo espacial secuencial habilitado; implementación todavía no iniciada']],[.34,.66]),
 text('Qué queda disponible','subheading'),
 text('Eventos electrón-fonón conservativos; transporte entre espectros distintos a energía fija; evolución del condensado y sus intercambios de energía; redistribución BGK, calentamiento y escape fonónico. El catálogo electrónico de etapa 1 se conserva como entrada. La malla complementaria seleccionada se usa explícitamente, sin reemplazar silenciosamente ese catálogo.'),
 text('Alcance del cierre','subheading'),
 text('Se acepta el paquete de desarrollo para construir la etapa 3. <b>No se declara completo el certificado dinámico de todas las mallas</b>, ni la calibración material de NbN, ni una predicción de pulsos del dispositivo. El solver de producción y el tag v1.0.0 permanecen preservados.'),
 text('La decisión operativa cambia por el cierre solicitado; los resultados, los criterios originales y los fallos históricos quedan intactos. No se exige otra corrida larga para entregar esta etapa. Los datos ausentes se investigarán sólo si afectan una decisión concreta del ensayo siguiente.','small'),
 ]),dict(title='1. La dinámica acoplada funciona en el ensayo',blocks=[
 text('Las dos celdas parten con amplitudes y poblaciones distintas. Durante el intervalo normalizado de 0 a 2 evolucionan el condensado, las excitaciones y el intercambio entre celdas. La prueba incluye fuente externa y escape al baño; no es un estado estacionario trivial.'),
 image('dynamics_guarded.png','Datos del integrador SSP protegido: líneas de 1280 pasos y marcadores de 320. El panel inferior izquierdo separa sus diferencias, invisibles a la escala de la curva principal. Energías y tiempos son normalizados; no equivalen todavía a un pulso de NbN en unidades experimentales.'),
 text('La energía electrónica de excitación se muestra separada de la fonónica. El transporte acumulado redistribuye energía interna. La contabilidad total incluye además el condensado y su trabajo. Estas trayectorias mantienen Gamma fijo; fase, potencial y trabajo espectral con Gamma(q) variable se comprobarán en etapa 3.'),
 text('Se adopta 320 pasos como referencia económica de desarrollo para este ensayo de duración 2; 640 pasos permite contrastar su sensibilidad temporal. Esta elección no fija el paso temporal de una futura malla espacial.','small'),
 ]),dict(title='2. Precisión medida y límite que permanece',blocks=[
 image('time_and_energy.png','Refinamiento propio del método protegido en 630 estados electrónicos y 1025 nodos fonónicos. La referencia de 1280 pasos es numérica, no una solución exacta. Se usan los puntos comunes registrados y la norma original de cada familia de observables.'),
 text('Con 320 pasos, la diferencia máxima registrada frente a 1280 es <b>0,001570 % en una celda y 0,001126 % en dos</b>. Con 640 baja a 0,000374 % y 0,000282 %, respectivamente. El defecto del balance disminuye también al refinar; conservar energía y aproximar bien la trayectoria son controles distintos.'),
 image('precision_context.png','Izquierda: diferencias del par temporal 320/640 en la malla de 2049 fonones. Derecha: sensibilidad entre 1025 y 2049 fonones a 320 pasos; incluye error temporal residual. Dos mallas no certifican convergencia.'),
 text('El par de 2049 fonones mide 0,004565 % frente al presupuesto auxiliar de 0,0025 %. Ese FAIL se conserva. La diferencia está principalmente en energías normalizadas superiores a 1, no en una cola infrarroja prescindible. Al cambiar 1025 por 2049 nodos, la amplitud cambia 0,000420 %, la energía fonónica 0,0263 % y el intercambio electrón-fonón 0,122 %.','small'),
 ]),dict(title='3. La mejora necesaria fue numérica',blocks=[
 text('El refinamiento electrónico de RK4 había producido una población no física en una etapa interna. Se adoptó SSPRK3 con factores comunes para cada evento conservativo. Después, una prueba aislada reveló un segundo problema: aritmética de números extremadamente pequeños, próxima al mínimo representable por float64.'),
 image('guarded_rounding_and_escape.png','Izquierda: reproducción de la etapa aislada que antes daba -5 unidades del mínimo positivo float64; el método protegido produce +7. Derecha: el ensayo de escape conserva la reducción cercana a ocho al duplicar los pasos, compatible con orden tres en ese caso.'),
 text('La corrección usa aritmética ampliada para inventarios y flujos sólo cerca del subdesbordamiento. No recorta poblaciones ni repara la energía después de integrar. En las 13 trayectorias finales no se activó esa protección ni el limitador: la rama ordinaria produjo los resultados. Los controles aislados comprueban la rama protegida cuando sí se necesita.'),
 table([['Costo registrado en Geminga','Una celda','Dos celdas'],
 ['320 pasos','3,28 min','6,45 min'],['640 pasos','6,54 min','12,89 min'],
 ['1280 pasos','13,02 min','25,72 min']],[.46,.27,.27]),
 text('El último lote terminó a los 94,5 minutos, con 21 de sus 35 tareas completas. La siguiente optimización debe reducir el costo de evaluar eventos y reutilizar datos compatibles antes de multiplicar celdas. El tiempo de una o dos celdas no permite prometer el costo del detector completo.'),
 text('El prototipo protegido requiere un tipo longdouble con mayor alcance que float64, disponible en Geminga. Los entornos sin esa capacidad se rechazan explícitamente. Esta dependencia debe resolverse o mantenerse declarada al promover el integrador.','small'),
 ]),dict(title='4. El circuito previsto es el de la memoria',blocks=[
 text('Se fija la red usada en los resultados de la memoria, sección 4.4.1, ecuación 4.16. Sus estados son la corriente de polarización, la corriente total del ramal del detector y el voltaje del capacitor. La producción ya contiene estas tres ecuaciones; el nuevo acoplamiento espacial las reutilizará con un puerto de energía consistente.'),
 image('circuit_comparison.png','Diagnóstico de la red con resistencia prescrita: 0/150/0 ohm durante 0-1/1-4/4-16 ns; corriente inicial 20 microamperios e inductancia serie exterior de prueba 10 nH. La red completa añade Rb = 10 kohm, Lb = 1 microhenrio, Cc = 100 pF y Vbias = 0,2 V. Es un diagnóstico circuital, no un transitorio del SNSPD.'),
 text('En este ejemplo, la red de la memoria alcanza 0,738 mV frente a 0,750 mV con la fuente ideal y presenta una cola negativa de -0,380 mV. La diferencia justifica mantener la red completa para comparar formas de pulso. La identidad instantánea de potencia cierra con residuo escalado de 3,19e-16.'),
 table([['Decisión de acoplamiento','Consecuencia para etapa 3'],
 ['Corriente y voltaje con signo pasivo','El mismo producto I·V entra en las cuentas del dispositivo y del circuito con signos opuestos'],
 ['Puerto sobre todo el dominio resuelto','El voltaje central de 100 nm queda como diagnóstico separado'],
 ['Inductancia exterior identificada una vez','No sumar otra vez la energía inductiva ya resuelta en el dominio'],
 ['Tres estados continuos al depositar energía','Inicializar el régimen estacionario antes de la perturbación']],[.42,.58]),
 text('Referencia histórica: Rb = 10 kohm, Lb = 1 microhenrio, RL = 50 ohm, Cc = 100 pF, Vbias = 0,300 V y Lk total = 10 nH. Los 10 nH son totales; no se asignan automáticamente al exterior del nuevo dominio. Ecuaciones y signos completos: adenda CM.1-CM.11.','small'),
 ]),dict(title='5. Etapa 3 preparada: espacio, bordes y circuito',blocks=[
 text('<b>Objetivo siguiente:</b> resolver un ensayo espacial dentro del dominio admitido del modelo, incluyendo terminales, intercambio con reservorios y la red circuital de la memoria. La etapa está preparada; aún no hay un resultado espacial nuevo.'),
 table([['Secuencia de implementación','Salida que permite continuar'],
 ['3A. Energía y campos espaciales','Energía, fuerza y corriente coherentes; prueba uniforme y perturbaciones pequeñas; dominio de estabilidad declarado'],
 ['3B. Bordes y reservorios','Balance de carga y energía con flujos de entrada y salida; caso abierto distinguido del control periódico'],
 ['3C. Potencial y corriente','Continuidad de carga, signos de terminales y voltaje del dominio con referencia de potencial definida'],
 ['3D. Circuito de la memoria','Régimen estacionario acoplado, tres estados circuitales, puerto pasivo y reparto de inductancia'],
 ['3E. Ensayo dinámico admitido','Respuesta débil y depósito localizado; balance global y precisión de los observables registrados']],[.35,.65]),
 text('La progresión exige comprobaciones propias de cada acoplamiento, no repetir automáticamente el lote de etapa 2. Se conserva el control negativo D.36: si la respuesta espacial pierde estabilidad en el estado probado, debe detectarse y excluirse ese dominio. No se oculta esa limitación refinando una malla.'),
 text('Antes de una corrida espacial costosa se fijan geometría, excitación, ventana temporal y precisión útil de amplitud, energías, corriente y voltaje de salida. El error numérico debe permitir decidir sobre la respuesta del dispositivo. Un umbral auxiliar previo no se convierte en requisito universal.'),
 text('Entrega y reproducción','subheading'),
 text('El cierre se registra en stage2/closure_20260922/closure_decision.json; la entrada y secuencia en stage3/entry_contract.json. Los datos originales, incluidos el FAIL y las tareas ausentes, permanecen en stage2/practical_review_20260922/raw/. saved_results_audit.json verifica los archivos y conserva las comparaciones. verify_delivery.py comprueba la entrega sin integrar trayectorias.'),
 text('La libreta /home/jdiaz/GEMINGA_COMMANDS.md queda sin corridas largas pendientes de esta etapa y sin instrucciones de screen. Los cálculos futuros de más de cinco minutos se prepararán allí para ejecución del usuario. El informe se regeneró con datos guardados, sin nuevas evaluaciones físicas.','small'),
 ])]
content=dict(title='Informe final de la etapa 2: celdas acopladas',
    subtitle='Cierre de desarrollo y preparación de etapa 3 | 22 de septiembre de 2026',pages=pages)
(DATA/'report_content.json').write_text(json.dumps(content,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
builder=(ROOT/'sandbox/stage2_cells/practical_review_20260922/build_report.py').read_text(encoding='utf-8')
builder=builder.replace('docs/implementation/stage2/practical_review_20260922','docs/implementation/stage2/closure_20260922')
builder=builder.replace('Informe_revision_criterios_y_circuito_20260922','Informe_cierre_etapa_2')
(Path(__file__).parent/'build_report.py').write_text(builder,encoding='utf-8')
record=dict(status='PREPARED_FROM_SAVED_RESULTS',new_rhs_evaluations=0,new_trajectories=0,
    source_audit_sha256=sha(PREVIOUS/'saved_results_audit.json'),
    inputs={str(p.relative_to(ROOT)).replace('\\','/'):sha(p) for p in
        [*(RAW/f'{case}_guarded_{n}.json' for case in runs for n in (160,320,640,1280)),
         *(RAW/f'{case}_guarded_time_assessment.json' for case in runs),*copied.values()]},
    figures={p.name:sha(p) for p in sorted(FIG.glob('*.png'))})
(DATA/'report_provenance.json').write_text(json.dumps(record,indent=2)+'\n',encoding='utf-8')
print(json.dumps(dict(pages=len(pages),figures=len(record['figures']),new_rhs_evaluations=0)))
