"""Build a results-only report from audited temporal and closure diagnostics."""
from pathlib import Path
import json
ROOT=Path(__file__).resolve().parents[3]
DATA=ROOT/'docs/implementation/stage2/closure_prep_20260922'
load=lambda p:json.loads(p.read_text(encoding='utf-8'))
audit=load(DATA/'two_cell_audit.json');state=load(DATA/'closure_status.json')
assert audit['checks_passed']==audit['checks_total']==98
runs={n:load(DATA/f'raw_two/two_ssp_{n}.json') for n in (160,320,640,1280)}
def text(value,kind='body'):return dict(type=kind,text=value)
def image(name,caption):return dict(type='image',path=(DATA/'figures'/name).relative_to(ROOT).as_posix(),caption=caption)
def table(rows,widths):return dict(type='table',rows=rows,widths=widths)
pages=[dict(title='Dos celdas: tiempo validado',blocks=[
text('<b>El lote completó sus seis tareas y pasa 98 de 98 controles de auditoría.</b> Las resoluciones de 160, 320 y 640 pasos cumplen su propio criterio temporal frente a 1280. Los balances incluyen escape al baño y calentamiento externo una sola vez.'),
table([['Pasos','Error máximo','Defecto energético','Tiempo (min)']]+[[str(n),f"{audit['rows'][i]['worst']['error']:.3e}",f"{runs[n]['energy_ledger_scaled_max']:.3e}",f"{runs[n]['runtime_seconds']/60:.2f}"] for i,n in enumerate((160,320,640))]+[['1280','Referencia',f"{runs[1280]['energy_ledger_scaled_max']:.3e}",f"{runs[1280]['runtime_seconds']/60:.2f}"]],[.13,.28,.32,.27]),
image('two_time_convergence.png','Las curvas representan errores medidos, no resultados experimentales. Cada familia se compara con su propia referencia de 1280 pasos. La energía usa el criterio registrado de 1e-7.'),
text('<b>Se eligen 320 pasos como base para comparar mallas.</b> La pareja 160/320 tiene un error de 4,78e-5 y excede el presupuesto suplementario de 2,5e-5. La pareja 320/640 queda en 1,24e-5 y sí pasa. Esto conserva la tolerancia general de 1e-4 y reserva precisión para distinguir el error de malla.'),
text('Lote observado: 48,54 min, un hilo, 630 estados electrónicos y 1025 nodos fonónicos por celda. Intervalo t/t_ref = 0 a 2; escape 15 y fuente externa 0,01 en las unidades declaradas. El catálogo, las leyes físicas y el integrador que produjeron estos datos se preservan por hash.','small')
]),dict(title='Intercambio entre celdas resuelto',blocks=[
text('Las dos amplitudes evolucionan desde estados distintos y el transporte redistribuye energía. El cambio de amplitud máximo es <b>0,199</b>, por encima del mínimo de 0,02. El intercambio electrónico-fonónico y el transporte también superan sus umbrales de actividad.'),
image('two_cell_dynamics.png','Líneas: referencia SSP1280; círculos superiores: SSP320. El panel de diferencias muestra la separación que oculta la superposición. El transporte es interno al sistema; la fuente y el escape intervienen en su balance externo.'),
text('En SSP320, las transferencias equivalen al <b>4,78% y 11,08%</b> de la excitación inicial. Las ocupaciones son físicas, sin recortes ni reparación energética. Los controles de campos y soporte pasan en 22 muestras: identidad derivativa 5,39e-11 y cota de absorción fuera del soporte 6,66e-10.'),
text('<b>El limitador no actuó en estas cuatro corridas.</b> Por ello, la validación temporal de la malla candidata no acredita todavía la trayectoria fina de 2520 estados, donde podría activarse. Los resultados corresponden a entradas Debye sintéticas; no predicen aún tasas absolutas de NbN ni un pulso de detector.','small')
]),dict(title='Un fallo de redondeo, una corrección verificable',blocks=[
text('<b>El ensayo de transporte aislado rechazó correctamente una ocupación negativa.</b> Ocurrió al multiplicar y dividir números próximos al mínimo positivo de float64: el inventario se redondeó antes de calcular cuánto podía salir de un estado electrónico. El margen de seguridad dejó de protegerlo.'),
image('guarded_rounding_and_escape.png','Izquierda: valores medidos en el electrón que provocó el rechazo. La corrección evalúa de nuevo esa misma etapa. Derecha: error del escape frente a la solución exponencial exacta; el mapa protegido reproduce bit a bit las tres trayectorias anteriores.'),
text('La versión protegida conserva el mapa anterior en pasos ordinarios y usa <b>aritmética extendida en los inventarios, límites y actualizaciones</b> cuando se aproxima el subdesbordamiento. Mantiene un flujo común por evento. No cambia los kernels físicos, recorta poblaciones ni repara energía después del paso.'),
table([['Verificación de la corrección','Resultado medido'],['Trayectoria ordinaria de 5 pasos, 630/1025','Idéntica bit a bit; 15 etapas delegadas'],['Etapa antes rechazada','Pasa de -5q a +7q; sin clipping'],['Transporte aislado completo, 20 pasos','Energía: 2,02e-15; conteo: 1,25e-15'],['Activación de aritmética extendida','27 de 60 etapas de transporte'],['Escape, 10/20/40 pasos','Reducciones de error 8,09 y 8,05']],[.57,.43]),
text('El redondeo final a float64 puede representar como cero un valor positivo menor que su rango; se registra expresamente y no se usa como reparación. La protección requiere que longdouble tenga rango y precisión superiores a float64: se verificó en Geminga y se rechazó explícitamente en este Windows.','small'),
text('<b>Estos ensayos acreditan la corrección localizada, no su convergencia global.</b> Los PASS temporales de las páginas anteriores pertenecen al mapa histórico; no se transfieren automáticamente al nuevo.','small')
]),dict(title='Comprobaciones adicionales y dictamen',blocks=[
text(state['summary']),
table([['Comprobación','Resultado','Alcance']]+state['report_gate_rows'],[.27,.23,.50]),
text('Decisión para continuar','subheading'),
text(state['next_action']),
text(state['remaining_runtime_note']),
text('Estado de entrega','subheading'),
text('La etapa 2 permanece sin admisión global mientras falte una prueba obligatoria. El push de cierre sigue condicionado a esa admisión; la etapa 3 está preparada y no implementada. Los datos, diagnósticos y comandos se sincronizan con Geminga.','small'),
text('Trazabilidad: closure_prep_20260922/two_cell_audit.json, closure_status.json y sus evidencias; criterio fijo stage2/acceptance_criteria.json; estado global stage2/stage2_admission.json. Los archivos fallidos y los informes anteriores se conservan.','small')
])]
content=dict(title='Etapa 2 - Validación de dos celdas y pruebas de cierre',subtitle='Resultados y decisión de continuación | 22 de septiembre de 2026',pages=pages)
(DATA/'report_content.json').write_text(json.dumps(content,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
builder=(ROOT/'sandbox/stage2_cells/review_20260922/build_report.py').read_text(encoding='utf-8').replace("DATA = ROOT / 'docs/implementation/stage2/review_20260922'","DATA = ROOT / 'docs/implementation/stage2/closure_prep_20260922'").replace('Informe_revision_etapa_2_20260922','Informe_validacion_dos_celdas_etapa_2_20260922')
(Path(__file__).parent/'build_report.py').write_text(builder,encoding='utf-8')
print('Report content prepared from measured closure status.')
