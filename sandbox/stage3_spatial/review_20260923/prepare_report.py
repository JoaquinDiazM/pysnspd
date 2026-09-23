"""Result-focused review from preserved files; no physical computations."""
from pathlib import Path
import json
import re
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[3]
DATA=ROOT/'docs/implementation/stage3/review_20260923'
NEW=ROOT/'docs/implementation/stage3/nodal_20260923'
FIG=DATA/'figures';FIG.mkdir(exist_ok=True)
read=lambda p:json.loads(p.read_text(encoding='utf-8'))
old=read(DATA/'raw/static_campaign/summary.json')
pilot=read(NEW/'pilot/weak_phase_thermal_n8.json')
tests=re.search(r'(\d+) passed in ([\d.]+)s',(NEW/'checks/pytest.log').read_text()).groups()
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,
                    'axes.spines.right':False,'savefig.dpi':190})
fig,ax=plt.subplots(1,2,figsize=(10,3.7),layout='constrained')
rows=[v for v in old['grid_comparisons'] if 'metrics' in v]
labels=['Amplitud térmica','Fase térmica','Fase no térmica','Amplitud vacío']
for j,(key,label,color) in enumerate([
    ('excess_energy_density_bar','Energía','#3478ac'),
    ('induced_force_modes','Fuerza','#db7922'),
    ('induced_current_modes','Corriente','#00897a')]):
    values=[100*v['metrics'][key]['relative_estimate']
            if v['metrics'][key]['relative_estimate'] is not None else np.nan for v in rows]
    ax[0].bar(np.arange(4)+(j-1)*.23,values,width=.21,label=label,color=color)
ax[0].axhline(1,color='#333333',ls='--',label='Objetivo 1 %')
ax[0].set(yscale='log',xticks=range(4),xticklabels=labels,ylabel='Estimación espacial (%)',
          title='Campaña terminada: precisión pendiente',ylim=(.75,100))
ax[0].tick_params(axis='x',rotation=20);ax[0].legend(fontsize=8);ax[0].grid(alpha=.2)
meshes=[8,16,32];errors=[]
for n in meshes:
    result=read(DATA/f'raw/static_campaign/uniform_thermal_n{n}.json')
    ref=result['uniform_reference']
    errors.append(100*abs(ref['link_current_bar']/ref['continuum_current_bar']-1))
ref=pilot['uniform_reference']
new_error=100*abs(ref['link_current_bar']/ref['continuum_current_bar']-1)
ax[1].plot(meshes,errors,'o-',lw=2,label='Punto medio: lote terminado',color='#b44743')
ax[1].scatter([8],[new_error],s=85,marker='D',color='#007f73',label='Nodal: piloto nuevo N=8',zorder=3)
ax[1].axhline(1,color='#333333',ls='--',label='Referencia visual 1 %')
ax[1].set(yscale='log',xticks=meshes,xlabel='Celdas de la misma tira',
          ylabel='Diferencia con corriente continua (%)',title='Sesgo uniforme de corriente')
ax[1].legend(fontsize=8);ax[1].grid(alpha=.2)
fig.savefig(FIG/'spatial_review.png');plt.close(fig)

def text(value,kind='body'):return dict(type=kind,text=value)
def table(rows,widths):return dict(type='table',rows=rows,widths=widths)
content=dict(title='Etapa 3: revisión del lote y siguiente ejecución',
 subtitle='Resultados del ensayo espacial | 23 de septiembre de 2026',pages=[
 dict(title='El lote terminó; falta precisión espacial',blocks=[
 text(f'Los <b>18 casos terminaron en {old["runtime_seconds"]/60:.1f} minutos</b>. Pasaron las derivadas de energía, los controles de fase, los signos locales y los contrastes espectrales muestreados. El dictamen final fue precisión espacial insuficiente: no hubo un fallo de integración temporal ni una ejecución interrumpida.'),
 table([['Respuesta en la malla de 32 celdas','Estimación de error'],
        ['Amplitud térmica: energía / fuerza / corriente','1,86 % / 2,00 % / 69,78 %'],
        ['Fase térmica: energía / fuerza / corriente','5,89 % / 8,66 % / 6,09 %'],
        ['Fase no térmica: energía / fuerza / corriente','5,87 % / 8,76 % / 6,07 %'],
        ['Control de gradiente aislado','0,321 %: aprobado']],[.64,.36]),
 dict(type='image',path=(FIG/'spatial_review.png').relative_to(ROOT).as_posix(),
      caption='Izquierda: estimaciones condicionales al refinar 8/16/32 celdas; la corriente nula del vacío no lleva porcentaje. Derecha: corriente uniforme frente al límite continuo, con población congelada. El punto verde es el nuevo piloto; aún no es una curva de convergencia nodal.'),
 text('La causa principal es geométrica: promediar dos números complejos con distinta fase acorta su amplitud, aunque ambos nodos tengan el mismo módulo. En el estado prescrito esa reducción introduce corriente numérica adicional. Una expansión analítica explica el 98,23 % del sesgo uniforme de vacío a 32 celdas.'),
 text('Los 11 fallos espaciales se conservan. Los 10 contrastes electrónicos 630/1260 mostraron cambios de matriz de hasta 4,20e-9; no aportan evidencia de que aumentar globalmente el espectro sea el siguiente paso útil.','small')]),
 dict(title='La siguiente corrida cambia la discretización',blocks=[
 text('Se implementó una energía <b>nodal de cuarto orden</b>: consulta la amplitud en cada nodo y usa operadores covariantes compatibles para gradiente y corriente. Fuerza y corriente siguen siendo derivadas de la misma energía. El modelo continuo, el catálogo, la regularización y el objetivo del 1 % permanecen iguales.'),
 table([['Verificación nueva','Resultado'],
        ['Suite del esquema nodal y controles históricos',f'{tests[0]} pruebas aprobadas en {tests[1]} s'],
        ['Piloto nodal de ocho celdas',f'Aprobado en {pilot["runtime_seconds"]:.1f} s'],
        ['Corriente uniforme del piloto frente al continuo',f'Diferencia {new_error:.3f} %'],
        ['Convergencia espacial del esquema nodal','Pendiente del lote del usuario']],[.64,.36]),
 text('La energía de gradiente también controla oscilaciones alternantes entre nodos. Usar solamente una derivada central de orden alto dejaría ese modo sin penalización; el operador positivo complementario evita ese defecto.'),
 text('El nuevo lote mantiene los seis perfiles y las mallas de 8/16/32 celdas. Las diferencias independientes de energía se ejecutan para cada caso en 16 celdas, además del piloto de ocho; no se repiten las 24 consultas derivativas en cada malla. Los controles de fase, referencias uniformes, signos y respuestas espaciales se mantienen en las tres mallas.'),
 text('<b>Próximo paso:</b> ejecutar el bloque activo de /home/jdiaz/GEMINGA_COMMANDS.md. La terminal muestra barras, tiempo transcurrido y ETA; progress.jsonl guarda el avance. Se reutiliza únicamente el piloto nodal con fuentes y estados verificados.'),
 text('Este lote decide el cierre del bloque estático 3A. Después quedan bordes y reservorios, potencial eléctrico, las tres ecuaciones del circuito de la memoria y la dinámica débil con depósito sintético. La etapa 3 completa sigue abierta; todavía no se ha calculado una señal de detección.','small')])])
(DATA/'report_content.json').write_text(json.dumps(content,ensure_ascii=False,indent=2)+'\n',encoding='utf-8',newline='\n')
builder=(ROOT/'sandbox/stage3_spatial/build_report.py').read_text(encoding='utf-8')
builder=builder.replace('parents[2]','parents[3]').replace('stage3/iteration_20260923','stage3/review_20260923')
builder=builder.replace('Informe_inicio_etapa_3_20260923','Informe_revision_espacial_etapa_3_20260923')
builder=builder.replace('Inicio espacial','Revisión espacial')
(Path(__file__).parent/'build_report.py').write_text(builder,encoding='utf-8',newline='\n')
print(json.dumps(dict(new_pilot_uniform_current_error_percent=new_error,tests=int(tests[0]))))
