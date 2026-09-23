"""Prepare figures and the stage-3 start report from the saved pilot only."""
from pathlib import Path
import hashlib
import json
import re
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[2]
DATA=ROOT/'docs/implementation/stage3/iteration_20260923'
FIG=DATA/'figures';FIG.mkdir(exist_ok=True)
read=lambda p:json.loads(p.read_text(encoding='utf-8'))
record=read(DATA/'pilot/weak_phase_thermal_n8.json')
negative=read(DATA/'pilot/negative_control.json')
test_match=re.search(r'(\d+) passed in ([\d.]+)s', (DATA/'checks/pytest.log').read_text())
test_count,test_seconds=test_match.groups()
assert record['status']=='PASS_STATIC_IDENTITIES_AND_SAMPLED_SYMBOL'
a=record['observable_arrays'];x=np.array(a['x_over_ell0'])*record['ell0_m']*1e9
mid=x+.5*record['h_bar']*record['ell0_m']*1e9
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,
    'axes.spines.right':False,'savefig.dpi':200,'axes.titlesize':11})
fig,ax=plt.subplots(1,2,figsize=(10,3.7),layout='constrained')
phase=np.array(a['phase'])-.1*np.array(a['x_over_ell0'])
ax[0].plot(x,phase,'o-',lw=1.8,color='#3274b2',label='Modulación de fase nodal')
ax[0].set(title='Perfil débil prescrito, ocho celdas',xlabel='Posición longitudinal (nm)',ylabel='Fase adicional (rad)')
ax[1].plot(mid,np.array(a['current_A'])*1e6,'s-',lw=1.8,color='#007f73',label='Corriente conjugada de la energía')
ax[1].set(title='Respuesta estática calculada',xlabel='Punto medio del enlace (nm)',ylabel='Corriente (µA)')
for b in ax:b.grid(alpha=.2);b.legend(fontsize=8,loc='best')
fig.savefig(FIG/'pilot_response.png');plt.close(fig)
fig,ax=plt.subplots(1,2,figsize=(10,3.7),layout='constrained')
for i,(key,label) in enumerate((('force','Fuerza'),('current','Corriente'))):
    v=record['derivative_checks'][key]
    ax[0].bar(i-.16,v['absolute_error'],width=.3,color='#007f73',label='Diferencia medida' if i==0 else None)
    ax[0].bar(i+.16,v['tolerance'],width=.3,color='#cc6e13',label='Presupuesto registrado' if i==0 else None)
ax[0].set(yscale='log',ylim=(1e-14,1e-3),xticks=[0,1],xticklabels=['Fuerza','Corriente'],
    title='Derivadas de la misma energía',ylabel='Diferencia absoluta normalizada')
ax[0].legend(fontsize=8);ax[0].grid(axis='y',alpha=.2)
values=[min(s['eigenvalues'][0] for s in record['principal_symbols']),negative['measured_eigenvalues'][0]]
ax[1].bar(['Piloto suave','Control inestable'],values,color=['#007f73','#b54646'],width=.52)
ax[1].axhline(0,color='#555555',lw=.9)
for i,value in enumerate(values):ax[1].text(i,value+(.07 if value>0 else -.07),f'{value:.4f}',ha='center',va='bottom' if value>0 else 'top')
ax[1].set(ylim=(-.65,1.95),title='El control distingue el signo esperado',ylabel='Menor valor propio local D.36')
fig.savefig(FIG/'pilot_checks.png');plt.close(fig)

def text(s,kind='body'):return dict(type=kind,text=s)
def image(name,caption):return dict(type='image',path=(FIG/name).relative_to(ROOT).as_posix(),caption=caption)
def table(rows,widths):return dict(type='table',rows=rows,widths=widths)
content=dict(title='Etapa 3: implementación inicial y lote preparado',
 subtitle='Piloto espacial estático | 23 de septiembre de 2026',pages=[
 dict(title='Etapa 3 iniciada: el funcional espacial',blocks=[
 text('<b>Se implementó el primer bloque 3A:</b> una energía espacial discreta que produce tanto la fuerza cartesiana como la corriente de los enlaces. El flujo regularizado modifica el espectro local; sus contribuciones se incluyen en las dos derivadas. La implementación permanece experimental.'),
 table([['Comprobación realizada en Geminga','Resultado'],
 ['Pruebas del funcional, progreso y reuso',f'{test_count} pruebas aprobadas en {test_seconds.replace(".",",")} s'],
 ['Piloto registrado de ocho celdas',f'Aprobado en {record["runtime_seconds"]:.1f} s'],
 ['Fuerza frente a variación de energía','Diferencia absoluta 7,08e-12; presupuesto 1,80e-4'],
 ['Corriente frente a variación de energía','Diferencia absoluta 3,42e-13; presupuesto 9,83e-5'],
 ['Fase global y referencia helicoidal','Invariancias y signos comprobados']],[.49,.51]),
 image('pilot_response.png','Tira de 360 nm, sección de 120 nm × 7 nm. Campo prescrito de amplitud 0,9 Delta0, flujo base q ell0 = 0,1 y modulación de fase de 0,1 rad. Población preparada térmicamente y luego congelada. Las corrientes usan las escalas electrónicas de referencia del catálogo.'),
 text('El perfil es una prueba estática prescrita, no una solución estacionaria del detector. Por eso la corriente puede variar entre enlaces: todavía no se resolvió el potencial que impondrá continuidad de corriente total. Los bordes físicos y el circuito se incorporarán después de este primer bloque.'),
 text('Las fuerzas y la corriente se verificaron con diferencias independientes de la energía, sin termalizar de nuevo las poblaciones al perturbar los campos. La prueba distingue la corriente conjugada de la discretización de una fórmula continua simplemente muestreada.','small')]),
 dict(title='Validación pendiente y ejecución con progreso',blocks=[
 image('pilot_checks.png','Izquierda: las diferencias observadas quedan por debajo de presupuestos fijados antes de ejecutar el piloto. Derecha: el estado suave tiene signo positivo en la malla espectral candidata y el control negativo conserva su signo inestable. Estos valores no certifican todos los estados del modelo.'),
 text('El menor valor propio del piloto es 1,5708; la mayor incertidumbre por diferencias finitas del símbolo es 1,04e-6. El control de amplitud 0,6 Delta0 y q ell0 = 1 da -0,30662 en estas unidades y se rechaza. Queda pendiente el contraste espectral 630/1260 del lote completo; la incertidumbre indicada aquí no lo sustituye.'),
 table([['Lote preparado para el usuario','Alcance'],
 ['Seis casos en 8, 16 y 32 celdas','Vacío, población térmica y no térmica; gradientes débiles de amplitud o fase'],
 ['Piloto de ocho celdas','Se reutiliza tras comprobar fuentes, condiciones y hashes'],
 ['Precisión espacial','Objetivo del 1 % en respuestas registradas; incertidumbre y orden observado informados'],
 ['Costo estimado','Aproximadamente 15-20 min; margen orientativo 12-25 min, un proceso CPU y menos de 1 GB de RAM estimado'],
 ['Salida de terminal','Barra del caso y del lote, avance completado, tiempo transcurrido y ETA aproximada'],
 ['Archivos persistentes','Resultados, estados iniciales, recibos de integridad y progress.jsonl']],[.37,.63]),
 text('<b>El lote largo queda preparado y no fue lanzado por el agente.</b> La orden exacta está al inicio de /home/jdiaz/GEMINGA_COMMANDS.md y no usa instrucciones de screen. La estimación se ajusta con el trabajo terminado y nunca se interpreta una evaluación en curso como avance aceptado.'),
 text('Al terminar, basta comunicar que concluyó; los resultados se leerán desde Geminga. Su revisión decidirá el alcance admitido para continuar con bordes y reservorios. La etapa 3 completa permanece abierta: todavía no hay evolución espacial con circuito ni un pulso de detección.','small')])])
(DATA/'report_content.json').write_text(json.dumps(content,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
builder=(ROOT/'sandbox/stage2_cells/closure_20260922/build_report.py').read_text(encoding='utf-8')
builder=builder.replace('parents[3]','parents[2]').replace('docs/implementation/stage2/closure_20260922','docs/implementation/stage3/iteration_20260923')
builder=builder.replace('Informe_cierre_etapa_2','Informe_inicio_etapa_3_20260923')
builder=builder.replace('pySNSPD · Etapa 2 - Celdas acopladas · 22 septiembre 2026','pySNSPD · Etapa 3 - Inicio espacial · 23 septiembre 2026')
(Path(__file__).parent/'build_report.py').write_text(builder,encoding='utf-8')
print('Two-page start report prepared from saved data; no physical calculations.')
