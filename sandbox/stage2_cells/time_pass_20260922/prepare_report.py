"""Postprocess the completed one-cell batch; no RHS or new dynamics."""
from pathlib import Path
import hashlib
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter

ROOT=Path(__file__).resolve().parents[3]
DATA=ROOT/'docs/implementation/stage2/time_pass_20260922'
FIG=DATA/'figures'
FIG.mkdir(exist_ok=True)
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
load=lambda p:json.loads(p.read_text(encoding='utf-8'))
audit=load(DATA/'one_cell_audit.json')
assert audit['checks_passed']==audit['checks_total'] and not audit['stage2_closure']
runs={n:load(DATA/f'raw/one_ssp_{n}.json') for n in (160,320,640,1280)}
assessment=load(DATA/'raw/one_ssp_time_assessment.json')
errors=np.array([r['max_error'] for r in assessment['cases']])
steps=np.array([160,320,640])
ledgers=np.array([runs[n]['energy_ledger_scaled_max'] for n in (*steps,1280)])
colors=['#cc6e13','#3274b2','#007f73']
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.titlesize':11,
                     'axes.spines.top':False,'axes.spines.right':False,'savefig.dpi':190})
fig,ax=plt.subplots(1,2,figsize=(10,3.8),layout='constrained')
ax[0].plot(steps,errors,'o-',color=colors[1],lw=2,label='Error frente a SSP1280')
ax[0].axhline(1e-4,color='#b23737',ls='--',label='Límite registrado: 1e-4')
for n,e in zip(steps,errors):ax[0].annotate(f'{e:.2e}',(n,e),xytext=(0,-18),textcoords='offset points',ha='center',fontsize=9)
ax[0].set(yscale='log',xscale='log',ylim=(1.5e-6,2.5e-4),title='Precisión de observables y poblaciones',ylabel='Máximo error relativo',xlabel='Pasos en el intervalo [0, 2]')
ax[0].set_xticks(steps,steps);ax[0].legend(fontsize=8,loc='lower left')
ax[1].plot([160,320,640,1280],ledgers,'s-',color=colors[2],lw=2,label='Defecto del balance energético')
ax[1].axhline(1e-7,color='#b23737',ls='--',label='Límite registrado: 1e-7')
ax[1].plot([160,320,640,1280],ledgers[0]*(160/np.array([160,320,640,1280]))**3,':',color='#777777',label='Pendiente cúbica, guía visual')
ax[1].set(yscale='log',xscale='log',ylim=(3e-11,4e-7),title='Energía: reducción cercana a factor 8',ylabel='Máximo defecto escalado',xlabel='Pasos en el intervalo [0, 2]')
ax[1].set_xticks([160,320,640,1280],[160,320,640,1280]);ax[1].legend(fontsize=8,loc='lower left')
for a in ax:
    a.xaxis.set_minor_formatter(NullFormatter())
    a.grid(True,alpha=.2,which='both')
fig.savefig(FIG/'temporal_convergence.png');plt.close(fig)

ref=runs[1280];t=np.array(ref['times']);sn=ref['snapshots']
series=lambda name:np.array([s[name] for s in sn]).reshape(len(t),-1)[:,0]
fig,ax=plt.subplots(2,2,figsize=(10,6.0),layout='constrained')
ax[0,0].plot(t,series('amplitudes'),color=colors[2],lw=2,label='Referencia SSP1280')
coarse=runs[160];tc=np.array(coarse['times']);sc=coarse['snapshots']
ax[0,0].plot(tc[::10],np.array([s['amplitudes'][0] for s in sc])[::10],ls='none',marker='o',ms=4,mfc='white',mec=colors[0],label='SSP160 (marcas espaciadas)')
ax[0,0].set(title='Recuperación del condensado',ylabel=r'$|\Delta|/\Delta_0$');ax[0,0].legend(fontsize=8)
for name,label,col in [('electron_energy','Sector electrónico',colors[1]),('phonon_energy','Fonones',colors[0]),('total','Total',colors[2])]:
    ax[0,1].plot(t,series(name)-series(name)[0],label=label,color=col,lw=2)
ax[0,1].set(title='Redistribución de energía',ylabel=r'Variación / $(N_0\Delta_0^2)$');ax[0,1].legend(fontsize=8)
for n,col in zip(steps,colors):
    current=runs[n];values=np.array([s['amplitudes'][0] for s in current['snapshots']]);truth=series('amplitudes')[::1280//n]
    ax[1,0].plot(current['times'],1e8*(values-truth),color=col,label=f'{n} pasos',lw=1.4)
ax[1,0].axhline(0,color='#777777',lw=.6)
ax[1,0].set(title='Diferencias que la superposición oculta',ylabel=r'$(|\Delta|_N-|\Delta|_{1280})/\Delta_0$  ($10^{-8}$)');ax[1,0].legend(fontsize=8)
for name,label,col in [('electron_to_phonon','Electrones hacia fonones',colors[1]),('condensate_heat','Calor del condensado',colors[0])]:
    ax[1,1].plot(t,series(name),color=col,lw=2,label=label)
ax[1,1].set(title='Transferencias internas acumuladas',ylabel=r'Energía / $(N_0\Delta_0^2)$');ax[1,1].legend(fontsize=8)
for a in ax.flat:a.set_xlabel(r'Tiempo $t/t_{\mathrm{ref}}$');a.grid(True,alpha=.2)
fig.savefig(FIG/'one_cell_results.png');plt.close(fig)

def body(s):return {'type':'body','text':s}
def small(s):return {'type':'small','text':s}
def heading(s):return {'type':'subheading','text':s}
def table(rows,widths):return {'type':'table','rows':rows,'widths':widths}
def image(name,caption):return {'type':'image','path':(FIG/name).relative_to(ROOT).as_posix(),'caption':caption}
pages=[{'title':'Una celda: tiempo validado','blocks':[
body('<b>El lote termina correctamente y supera 85 de 85 controles de auditoría.</b> Las tres resoluciones medidas cumplen por separado el criterio de precisión temporal y de balance de energía. La etapa 2 completa sigue pendiente: este lote sólo ensaya una celda en la malla candidata.'),
table([['Pasos','Error máximo frente a 1280','Defecto energético','Tiempo (s)']]+[[str(n),f'{errors[i]:.3e}',f'{ledgers[i]:.3e}',f"{runs[n]['runtime_seconds']:.1f}"] for i,n in enumerate(steps)]+[['1280','Referencia',f'{ledgers[-1]:.3e}',f"{ref['runtime_seconds']:.1f}"]],[.13,.34,.28,.25]),
image('temporal_convergence.png','Cada resolución se compara en los mismos 161 tiempos físicos. Error máximo: observables y norma L1 de poblaciones ponderada por sus capacidades. La referencia de 1280 pasos no se certifica a sí misma.'),
body('Al duplicar los pasos, el error de la solución disminuye por factores <b>3,45 y 4,20</b>; el defecto energético disminuye aproximadamente por 8. El resultado rescata SSP160 para este escenario. El fallo anterior de SSP40 se conserva como evidencia y no se ha relajado ninguna tolerancia.'),
small('Lote: stage2_ssp_time_20260922; duración total 24,37 min en Geminga. Malla: 630 estados electrónicos y 1025 nodos fonónicos; intervalo t/t_ref = 0 a 2; Debye sintético; sin escape ni calentamiento externo.')
]},{'title':'Qué muestran las trayectorias','blocks':[
body(f"La amplitud pasa de <b>0,600000 a {ref['final']['amplitudes'][0]:.6f}</b>. La energía se redistribuye entre el sector electrónico y los fonones con un defecto global pequeño. El panel de diferencias hace visible la separación entre curvas que parecen idénticas al superponerlas."),
image('one_cell_results.png','Arriba se ve la evolución; abajo, la diferencia con la referencia y las transferencias internas. La transferencia electrónica-fonónica negativa indica absorción neta de energía por los electrones. El calor del condensado se contabiliza dentro del sistema una sola vez.'),
body('<b>El limitador de eventos estuvo inactivo en las cuatro corridas.</b> Las ocupaciones electrónicas permanecieron entre 0 y 1 y las fonónicas fueron no negativas, sin recortes ni reparación energética. Este resultado no comprueba el régimen fino con limitación activa.'),
small('Dos controles posteriores sobre SSP160 también pasan: energía y fuerzas en 11 tiempos físicos (defecto de identidad derivativa 5,07e-11), y soporte superior (cota relativa máxima 9,28e-10). No sustituyen la comprobación de las futuras trayectorias de dos celdas.'),
small('Las escalas de tiempo y energía son las declaradas por el ensayo sintético. No se infieren latencias absolutas de NbN, pulsos de salida ni comportamiento de un detector espacial completo. La forma de la DOS y los núcleos físicos no cambiaron en esta revisión.')
]},{'title':'Dictamen y trabajo restante','blocks':[
table([['Bloque','Resultado actual','Acción necesaria'],
['Una celda / tiempo SSP','Aprobado para 160, 320 y 640 pasos','Conservar estas corridas; no repetir.'],
['Dos celdas / tiempo SSP','Pendiente','160/320/640 frente a 1280, más control pareado para las futuras mallas.'],
['Mallas de energía','Pendiente','Tres mallas electrónicas y tres fonónicas; controlar el error temporal de cada malla.'],
['Equilibrio y poblaciones límite','Pendiente','Comprobar estacionariedad del sistema completo y fronteras en canales aislados.'],
['Campos y soporte','Una celda: aprobado. Dos celdas: pendiente','Al menos 10 tiempos distintos por celda; incluir referencia causal para cada campo empleado.'],
['Etapa 3 espacial','Preparada, no iniciada','Entrada condicionada a la admisión global de etapa 2.']],[.27,.29,.44]),
heading('Decisión'),
body('La evidencia nueva respalda continuar con el método numérico candidato. No aporta una razón para reformular las ecuaciones físicas. Tampoco permite declarar cerrada la etapa 2: una prueba temporal aprobada no sustituye la convergencia de malla ni los casos de dos celdas y frontera.'),
body('El siguiente comando manual ensaya dos celdas; se estima entre 45 y 65 minutos con un hilo. No se inicia una cadena costosa de mallas hasta conocer el paso temporal admisible en ese sistema. El plan completo de pendientes identifica además las comprobaciones que necesitan ampliar sus evaluadores; no se presentan como ensayos ya preparados o aprobados.'),
heading('Entrega y trazabilidad'),
small('Informe editable y figuras: docs/implementation/stage2/time_pass_20260922/. Auditoría: one_cell_audit.json. Resultados originales: raw/. Criterio vigente: stage2/acceptance_criteria.json. Estado global: stage2/stage2_admission.json. Preparación espacial: stage3/README.md y entry_contract.json. El push de cierre solicitado queda pendiente del cumplimiento global; los archivos de trabajo se sincronizan con Geminga.'),
small('Se preservan el catálogo R2, los archivos físicos que originaron las corridas, los fracasos anteriores y v1.0.0. Esta revisión usa postprocesamiento; no ejecuta nuevas trayectorias.')
]}]
(DATA/'report_content.json').write_text(json.dumps({'title':'Etapa 2 - Resultado del lote temporal','subtitle':'Conclusión del lote de una celda | 22 de septiembre de 2026 | Admisión global pendiente','pages':pages},ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
builder=(ROOT/'sandbox/stage2_cells/review_20260922/build_report.py').read_text(encoding='utf-8')
builder=builder.replace("DATA = ROOT / 'docs/implementation/stage2/review_20260922'", "DATA = ROOT / 'docs/implementation/stage2/time_pass_20260922'")
builder=builder.replace('Informe_revision_etapa_2_20260922','Informe_resultados_temporales_etapa_2_20260922')
(Path(__file__).parent/'build_report.py').write_text(builder,encoding='utf-8')
metrics={'schema':'pysnspd.stage2.one-cell-report-metrics.v1','new_RHS_evaluations':0,'stage2_closure':False,
         'audit_sha256':sha(DATA/'one_cell_audit.json'),'batch_minutes':audit['batch_elapsed_seconds']/60,
         'time_errors':errors.tolist(),'energy_defects':ledgers.tolist(),'limiter_inactive':True,
         'energy_defect_reductions':(ledgers[:-1]/ledgers[1:]).tolist(),
         'delta_initial':ref['initial']['amplitudes'][0],'delta_final':ref['final']['amplitudes'][0]}
(DATA/'report_metrics.json').write_text(json.dumps(metrics,indent=2)+'\n',encoding='utf-8')
print(json.dumps(metrics,indent=2))
