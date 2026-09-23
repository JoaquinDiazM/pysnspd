"""Results-only figures and report content from saved outputs; no physical solve."""
from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[3]
DATA=ROOT/'docs/implementation/stage3/ports_20260923'
read=lambda p:json.loads(p.read_text(encoding='utf-8'))
figures=DATA/'figures';figures.mkdir(exist_ok=True)
summary=read(DATA/'raw/nodal_campaign/summary.json')
pilot=read(DATA/'open_pilot/open_L360_E4/result.json')
electrical=read(DATA/'electrical/electrical_diagnostics.json')
plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})

fig,ax=plt.subplots(figsize=(9,4.3),layout='constrained')
names=['Amplitud térmica','Fase térmica','Fase no térmica']
for shift,key,label,color,marker in [(-.17,'excess_energy_density_bar','Energía de exceso','#b75a28','s'),
    (0,'induced_force_modes','Fuerza inducida','#146b83','o'),
    (.17,'induced_current_modes','Corriente inducida','#71449a','D')]:
    vals=[100*c['metrics'][key]['relative_estimate'] for c in summary['grid_comparisons'][2:5]]
    ax.scatter(vals,np.arange(3)+shift,label=label,color=color,marker=marker,s=60)
    for x,y in zip(vals,np.arange(3)+shift):ax.annotate(f'{x:.4g}%',(x,y),xytext=(5,5),textcoords='offset points',fontsize=8,color=color)
ax.axvline(1,color='#a32525',ls='--',label='Objetivo registrado: 1 %')
ax.set(xscale='log',xlim=(1e-4,3),yticks=np.arange(3),yticklabels=names,
       xlabel='Estimación relativa espacial (%) - escala logarítmica',title='3A: las nueve respuestas térmicas/no térmicas cumplen el objetivo')
ax.invert_yaxis();ax.grid(axis='x',alpha=.17);ax.legend(loc='lower right',fontsize=8)
fig.savefig(figures/'nodal_precision.png',dpi=200);plt.close(fig)

x=np.asarray(pilot['x_m'])*1e9;i=np.asarray(pilot['current_A']);iref=pilot['reference_current_A']
a=pilot['reference_amplitude_bar'];q=np.asarray(pilot['q_delta_bar'])/(a*a/(a*a+.01))
fig,axes=plt.subplots(1,2,figsize=(9,3.3),layout='constrained')
axes[0].plot((x[:-1]+x[1:])/2,1e6*(i/iref-1),'-o',color='#146b83',ms=4)
axes[0].set(xlabel='Posición longitudinal (nm)',ylabel='Desviación de corriente (ppm)',title=f'Referencia: {iref*1e6:.4f} µA')
axes[1].plot(x,100*(q/pilot['reference_q_bare_bar']-1),'-s',color='#b75a28',ms=4)
axes[1].set(xlabel='Posición longitudinal (nm)',ylabel='Desviación del gradiente de fase (%)',title='Extremos abiertos incluidos')
for ax in axes:ax.axhline(0,color='.5',lw=.8);ax.grid(alpha=.17)
fig.savefig(figures/'open_pilot.png',dpi=200);plt.close(fig)

def body(t):return {'type':'body','text':t}
def small(t):return {'type':'small','text':t}
def picture(name,caption):return {'type':'image','path':(DATA/name).relative_to(ROOT).as_posix(),'caption':caption}
def table(rows,widths):return {'type':'table','rows':rows,'widths':widths}
m=pilot['metrics'];cm=electrical['circuit']['metrics']
content={'title':'Etapa 3: de la malla a los bordes','subtitle':'Informe de avance con resultados · 23 septiembre 2026 · Modelo 0.4 experimental','pages':[
 {'title':'3A completado; etapa 3 en desarrollo','blocks':[
 body('<b>La nueva discretización supera el problema de precisión del ensayo anterior.</b> Los 18 casos finalizaron en 7,85 minutos. La auditoría verificó ocho fuentes, 76 archivos y 336 nodos. Hay diez estimaciones relativas aceptadas, ninguna fallida y una diferencia demasiado pequeña para dar un certificado relativo.'),
 picture('figures/nodal_precision.png','Estimadores conservadores a partir de 8/16/32 celdas. Se representan errores, no curvas superpuestas. Son estimaciones condicionales de estos perfiles; no cotas universales ni una prueba de orden cuatro para todas las corrientes.'),
 body('La peor estimación es <b>0,613 % en corriente</b>, por debajo del objetivo de 1 %. Su orden observado es 0,608: permite continuar el desarrollo en este alcance, pero debe revisarse si cambia la geometría o el estado. La energía de exceso del caso vacío conserva el dictamen <b>sin certificado relativo</b>: diferencia 1,835e-9 frente al piso registrado de 1e-8.'),
 table([['Comprobación independiente','Resultado'],['Gradiente espacial frente a referencia analítica','0,001646 % en la malla fina'],['12 contrastes espectrales: 630 / 1260 nodos','Signo positivo conservado; cambio máximo de matriz 4,05e-9'],['Estado de control físicamente inestable','Rechazado como se esperaba'],['Nuevos módulos y regresiones en Geminga','168 pruebas y 23 subpruebas pasaron; 9,87 s']], [.57,.43]),
 small('La primera campaña permanece archivada con sus fallos. Este cierre de desarrollo de 3A no declara terminado el detector espacial ni promueve el modelo a producción.') ]},
 {'title':'Bordes abiertos: primer piloto físico','blocks':[
 body('Se implementó un tramo abierto de 360 nm con 17 nodos. Los reservorios fijan la amplitud mediante el equilibrio de la misma energía; la corriente se impone por su trabajo en el borde. El piloto térmico pasó en <b>20,9 segundos</b> con corriente de referencia <b>8,634 µA</b>.'),
 picture('figures/open_pilot.png','Diferencias respecto del equilibrio uniforme. La corriente se obtiene de la energía discreta. Los puntos terminales están incluidos en el control del gradiente; no se identifican con un cierre periódico.'),
 table([['Observable','Piloto','Objetivo'],['Error máximo de corriente',f'{100*m["current_relative"]:.6f} %','1 %'],['Error del gradiente en extremos',f'{100*m["endpoint_q_relative"]:.5f} %','1 %'],['Residuo estacionario normalizado',f'{m["stationarity_absolute"]:.3g}','0,0001'],['Amplitud / escala del catálogo',f'{a:.8f}','Rama estable'],['Inductancia resuelta / exterior fija',f'{pilot["L_res_H"]*1e9:.4f} / {pilot["L_ext_H"]*1e9:.4f} nH','Total de referencia: 10 nH']], [.43,.34,.23]),
 body('La hélice inicial ya satisface la tolerancia: <b>cero iteraciones de relajación</b>. Este resultado comprueba consistencia del equilibrio y sus bordes, no demuestra relajación de un pulso. El intercambio electrónico con un baño fijo pasó pruebas de equilibrio, soporte y balances; las ocupaciones se comparan a la misma energía física.'),
 small('La partición de inductancia corresponde a esta longitud y a esta rama. Se fija antes de la dinámica. El sesgo de este piloto débil es distinto de los 30 µA históricos. Falta contrastar longitudes y resoluciones; tampoco se ha resuelto todavía la interfaz 2D-1D.') ]},
 {'title':'Circuito de la memoria y siguiente ejecución','blocks':[
 body('Están implementados el potencial con continuidad de corriente y los tres estados circuitales de la memoria: corriente de polarización, corriente del dispositivo y tensión del condensador. El control independiente impone un salto de resistencia de 0 a 1000 ohmios y compara la integración con una solución matricial exacta.'),
 picture('electrical/electrical_diagnostics.png','Ensayo del circuito aislado con resistencia prescrita. Las curvas no son una predicción de detección ni proceden de un transiente del material. Los 7 nH exteriores son una entrada de este control algebraico, independiente de la partición identificada en el piloto abierto.'),
 body(f'La diferencia máxima en la salida es <b>{cm["maximum_Vout_error_V"]:.2e} V</b>. El balance instantáneo de potencia y la invariancia frente al origen del potencial pasan. Las pruebas incluyen que el trabajo del puerto tiene el mismo signo y magnitud en dispositivo y circuito.'),
 body('<b>Siguiente paso preparado:</b> seis equilibrios abiertos con longitudes 360, 720 y 1080 nm, dos resoluciones por longitud. Estimación: 6-10 minutos, un núcleo, menos de 1 GB de RAM. El comando está en <b>/home/jdiaz/GEMINGA_COMMANDS.md</b>, con barra de avance y ETA; debe ejecutarlo el usuario según la política acordada.'),
 small('Tras revisar ese lote siguen el empalme 2D-1D y el ensayo débil acoplado en el tiempo, con sus balances y observables. No se requieren repetir los lotes estáticos ya aceptados sin que haya cambios. Datos, criterios previos, fuentes y pruebas: docs/implementation/stage3/ports_20260923/.') ]}
]}
(DATA/'report_content.json').write_text(json.dumps(content,indent=2,ensure_ascii=False)+'\n',encoding='utf-8')
print('Prepared saved-data figures and three-page content.')
