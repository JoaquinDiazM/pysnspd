"""Illustrated review of the completed focused core controls; no new physics."""
from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import NullLocator
from matplotlib.font_manager import findfont
from reportlab.pdfgen import canvas
from reportlab.lib.colors import HexColor
from reportlab.lib.styles import ParagraphStyle
from reportlab.platypus import Paragraph
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont

ROOT=Path(__file__).resolve().parents[2]
DATA=ROOT/'docs/implementation/stage4/followup_20260923'
OLD=ROOT/'docs/implementation/stage4/review_20260923'
FIG=DATA/'figures'
PDF=ROOT/'output/pdf/implementation/Informe_avance_etapa_4_nucleo_y_paralelismo.pdf'
COLORS=['#2267a9','#cf6b2f','#278378']
plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,'savefig.dpi':180})


def read(path):
    return json.loads(path.read_text(encoding='utf8'))


def figures():
    FIG.mkdir(parents=True,exist_ok=True)
    raw=DATA/'raw/stage4A_core_followup_20260923'
    fine=read(raw/'2d_suppressed_d0.1_16x8_fixed_boundary/result.json')
    reduced=read(raw/'2d_suppressed_d0.05_8x4_fixed_boundary/result.json')
    old=read(OLD/'boundary_response_review.json')['cases']
    fig,axes=plt.subplots(1,2,figsize=(10,3.2),layout='constrained')
    for key,label,color in zip(['memory_effective','Allmaras_reference','Korzh_reference'],['Heredado','Allmaras','Korzh'],COLORS):
        values=[old[f'2d_suppressed_d0.1_{m}']['mobility'][key]['constrained_total_heat'] for m in ['4x2','8x4']]
        values.append(fine['mobility_fixed_boundary'][key]['total_heat_bar'])
        axes[0].plot([153,561,2145],values,'o-',color=color,label=label)
        for x,v in zip([153,561,2145],values):
            if key=='Korzh_reference':axes[0].annotate(f'{v:.4f}',(x,v),xytext=(0,9),textcoords='offset points',ha='center',fontsize=8,color=color)
        ref=values[1]
        d05=reduced['mobility_fixed_boundary'][key]['total_heat_bar']
        i=['memory_effective','Allmaras_reference','Korzh_reference'].index(key)
        axes[1].bar(i-.18,100*(values[2]/ref-1),width=.34,color=COLORS[0],label='Refinar malla: δ = 0,1' if i==0 else None)
        axes[1].bar(i+.18,100*(d05/ref-1),width=.34,color=COLORS[1],label='Cambiar δ: misma malla' if i==0 else None)
    axes[0].set(xscale='log',xlabel='Nodos 2D',ylabel='Disipación interior normalizada',title='La respuesta integrada se estabiliza')
    axes[0].set_xticks([153,561,2145],['153','561','2145']);axes[0].xaxis.set_minor_locator(NullLocator())
    axes[0].set_ylim(.32,1.25);axes[0].legend(fontsize=8);axes[0].grid(alpha=.2)
    axes[1].set_xticks([0,1,2],['Heredado','Allmaras','Korzh'])
    axes[1].set(ylabel='Cambio relativo (%)',title='Ahora domina la sensibilidad del cierre',ylim=(-7,27))
    axes[1].axhline(0,color='gray',lw=.8);axes[1].legend(fontsize=8);axes[1].grid(axis='y',alpha=.2)
    fig.savefig(FIG/'heat_and_closure.png');plt.close(fig)

    fig,axes=plt.subplots(1,2,figsize=(10,3.4),layout='constrained')
    for path,color,marker,label in [(OLD/'raw/stage4A_spatial_20260923/2d_suppressed_d0.1_8x4/fields.npz',COLORS[0],'s','561 nodos'),
                                    (raw/'2d_suppressed_d0.1_16x8_fixed_boundary/fields.npz',COLORS[1],'o','2145 nodos')]:
        with np.load(path) as a:
            xy=a['coordinates_m']*1e9
            ids=np.flatnonzero(abs(xy[:,1])<1e-9);ids=ids[np.argsort(xy[ids,0])]
            axes[0].semilogy(xy[ids,0]-80,a['gamma'][ids],color=color,marker=marker,ms=3,label=label)
            if 'gamma_analytic' in a:
                axes[0].semilogy(xy[ids,0]-80,a['gamma_analytic'][ids],color='#222222',ls='--',label='Perfil analítico')
    axes[0].set(xlim=(-22,22),xlabel='x respecto del centro (nm)',ylabel='Γ/Δ₀, corte y = 0',title='Núcleo mejor resuelto')
    axes[0].legend(fontsize=8);axes[0].grid(alpha=.2)
    dx=[-.004550200013554329,fine['derivative_diagnostics']['center']['derivative_discrete'][0][1],.0028]
    axes[1].bar(['561 nodos','2145 nodos','Exacta'],dx,color=[COLORS[0],COLORS[1],'#666666'])
    for i,v in enumerate(dx):axes[1].annotate(f'{v:.5f}',(i,v),xytext=(0,6 if v>0 else -14),textcoords='offset points',ha='center',fontsize=9)
    axes[1].set(ylabel='Im(∂Δ̄/∂X) central',title='También importa el signo de la derivada',ylim=(-.006,.0045))
    axes[1].axhline(0,color='gray',lw=.8);axes[1].grid(axis='y',alpha=.2)
    fig.savefig(FIG/'core_gradient.png');plt.close(fig)

    # Pedagogical schematic, expressly not a fitted or simulated detector pulse.
    t=np.linspace(0,12,600);v=(1-np.exp(-np.maximum(t-1,0)/.6))*np.exp(-np.maximum(t-1,0)/4)
    threshold=.35;cross=t[np.flatnonzero(v>=threshold)[0]];stop=cross+1.6
    fig,ax=plt.subplots(figsize=(10,2.65),layout='constrained')
    ax.plot(t,v,color='#aab3bf',lw=2,ls='--',label='Continuación del mismo sistema')
    keep=t<=stop;ax.plot(t[keep],v[keep],color=COLORS[0],lw=2.8,label='Ventana inicial propuesta')
    ax.axhline(threshold,color='#555555',ls=':',label='Umbral por definir')
    ax.axvline(cross,color=COLORS[2],ls='--');ax.axvline(stop,color=COLORS[1],ls='--')
    ax.axvspan(cross,stop,color=COLORS[1],alpha=.13)
    ax.text(cross+.06,.08,'Cruce',color=COLORS[2]);ax.text(stop+.1,.14,'Cruce + margen',color=COLORS[1])
    ax.set(xlabel='Tiempo ilustrativo, sin escala física',ylabel='Vout ilustrativo',xticks=[],yticks=[],title='Se recorta el horizonte; no se cambia el circuito ni la evolución previa')
    ax.legend(fontsize=8,loc='upper right')
    fig.savefig(FIG/'observation_horizon.png');plt.close(fig)
    with np.load(DATA/'raw/stage4_radial_reference_20260923/R8_N256.npz') as data:
        rr=data['radius'];mask=rr<=4
        fig,axes=plt.subplots(1,2,figsize=(10,3.3),layout='constrained')
        for ax in axes:
            ax.plot(rr[mask],data['force_raw'][mask],color='#151e2e',lw=2.5,label='Usadel espacial, N = 256')
            ax.plot(rr[mask],data['force_with_leading_tail_estimate'][mask],color='#151e2e',ls='--',lw=1,label='Con estimación de cola')
            for fraction,color in zip([.05,.1,.2],COLORS):
                ax.plot(rr[mask],data[f'candidate_force_delta_{fraction:g}'][mask],color=color,label=f'Candidato δ/Δ₀ = {fraction:g}')
            ax.axhline(0,color='gray',lw=.6);ax.grid(alpha=.15)
            ax.set(xlabel='r / ℓ₀',ylabel='Fuerza de amplitud normalizada')
        axes[0].set(xlim=(0,.6),title='El pico espurio queda dentro del núcleo')
        axes[1].set(xlim=(.6,4),ylim=(-.1,1.9),title='Comparación fuera del centro')
        axes[1].legend(fontsize=7)
        fig.savefig(FIG/'radial_force.png');plt.close(fig)
    return fine,reduced


def main():
    fine,reduced=figures();PDF.parent.mkdir(parents=True,exist_ok=True)
    pdfmetrics.registerFont(TTFont('Report',findfont('DejaVu Sans')))
    pdfmetrics.registerFont(TTFont('ReportBold',findfont('DejaVu Sans:weight=bold')))
    pdfmetrics.registerFontFamily('Report',normal='Report',bold='ReportBold',italic='Report',boldItalic='ReportBold')
    c=canvas.Canvas(str(PDF),pagesize=(595.276,841.89));c.setTitle('Etapa 4: núcleo, paralelismo y horizonte de observación')
    style=ParagraphStyle('body',fontName='Report',fontSize=10.2,leading=15,textColor=HexColor('#253343'))
    small=ParagraphStyle('small',parent=style,fontSize=8.6,leading=12)
    def text(value,y,st=style):
        p=Paragraph(value,st);_,h=p.wrap(499,730);p.drawOn(c,48,y-h);return y-h-13
    def start(n,title,subtitle):
        c.setFillColor(HexColor(COLORS[0]));c.rect(0,832,596,10,fill=1,stroke=0)
        c.setFont('Report',9);c.drawString(48,799,'pySNSPD  |  ETAPA 4 SIN FOTÓN  |  23 SEPTIEMBRE 2026')
        c.setFillColor(HexColor('#162c46'));c.setFont('ReportBold',21);c.drawString(48,758,title)
        text(subtitle,738,small)
        c.setFillColor(HexColor('#718096'));c.setFont('Report',8)
        c.drawString(48,29,'Avance experimental. Producción y v1.0.0 sin cambios.');c.drawRightString(547,29,f'{n} / 5')
    def picture(name,y,h):
        c.drawImage(str(FIG/name),42,y-h,width=511,height=h,preserveAspectRatio=True,anchor='c',mask='auto');return y-h-16
    start(1,'Del error de malla al cierre físico','Dos casos completados en 46,09 min. Resultados preservados con sus fuentes y mapas originales.')
    y=text('<b>Los dos campos nuevos son espacialmente positivos:</b> 561 nodos con δ = 0,05 y 2145 con δ = 0,1. Desapareció el único signo negativo de la malla gruesa. La evidencia permite avanzar; no hace falta otra campaña dedicada sólo a refinar el mismo perfil.',695)
    y=picture('heat_and_closure.png',y,213)
    y=text('<b>De 561 a 2145 nodos, con δ = 0,1:</b> la energía cambia +0,00889 %, la fuerza interior +0,572 % y la disipación Korzh -2,551 %. Los otros pares de movilidad cambian -1,38 % y -1,96 %. Son sensibilidades medidas, no intervalos de confianza del dispositivo.',y)
    y=text('<b>A igual malla, δ = 0,05 frente a 0,1:</b> el calor cambia +21,24 % para Korzh, aunque la energía sólo cambia +0,00332 %. La movilidad pesa de forma diferente las componentes locales de fuerza. Una energía total cercana no garantiza una dinámica cercana.',y)
    y=text('<b>Decisión:</b> el contraste físico independiente ya descarta admitir el cierre actual en el núcleo (página 5). δ = 0,1 se conserva como evidencia del candidato ensayado. No seleccionar δ para mejorar un signo o un tiempo de respuesta. Las poblaciones de estos dos casos son sintéticas y fijas; todavía no representan un fotón.',y)
    text('Fuentes: independent_audit.json y raw/stage4A_core_followup_20260923. Los cuatro archivos de resultados/mapas y las once fuentes coinciden con sus identidades registradas. No se integró el tiempo.',y,small);c.showPage()

    start(2,'Qué resolvió la malla fina','El diagnóstico puntual y los observables integrados se interpretan juntos.')
    y=picture('core_gradient.png',695,235)
    y=text('El error L² relativo de Γ respecto de la derivada exacta baja de <b>4,255 % a 0,181 %</b>. El error L² de la derivada compleja es 0,349 %. La malla intermedia incluso tenía el signo incorrecto de la derivada central; la fina recupera el signo y se acerca a su magnitud.',y)
    y=text('Γ central aún supera al exacto en 23,7 %, pero el error absoluto es 1,26 × 10<super>-5</super>. Se conserva esa limitación puntual. La estabilidad y los cambios de fuerza/calor ya permiten separar el error espacial dominante anterior de la dependencia del cierre.',y)
    y=text('<b>Bordes controlados:</b> su velocidad, calor y trabajo de restricción son cero, con balance de potencia a precisión de redondeo. La corriente se informa sumada por sección. Estos campos prescritos no son una solución estacionaria de carga ni el circuito acoplado del detector.',y)
    y=text('<b>Referencia ya calculada:</b> un núcleo con amplitud prescrita y fase que gira alrededor de su centro. Se resolvió cómo se ajusta el espectro de Usadel en el espacio, y se comparó la fuerza con el cierre local a la misma temperatura. Es un problema auxiliar circular; no sustituye la geometría de la cinta de 80 nm.',y)
    y=text('Una comprobación adicional cerca de T<sub>c</sub> recupera el límite GL: la diferencia de energía de un núcleo real prescrito baja de 2,75 % a 0,53 % al pasar de T/T<sub>c</sub> = 0,95 a 0,99. Es un control separado; no cambia la temperatura Korzh de 0,9 K ni valida el regularizador de fase.',y)
    text('La comparación térmica y los controles anteriores a población congelada pertenecen a regímenes estadísticos diferentes. No se superponen sus fuerzas como si la diferencia fuese un error físico. La referencia no calcula aún una barrera, un vórtice autoconsistente o su cruce.',y,small);c.showPage()

    start(3,'Paralelismo con recursos reservados','16 núcleos físicos, 32 hilos y aproximadamente 123 GiB RAM; un único nodo NUMA.')
    y=text('La ejecución anterior era sucesiva: 8,63 min y 37,47 min. La nueva planificación limita el conjunto de procesos a <b>28 CPU lógicas</b>, dejando libres dos núcleos físicos completos. Se consulta afinidad, topología y memoria disponibles al iniciar; un límite menor reduce la concurrencia.',695)
    y=text('<b>Un presupuesto compartido:</b> las consultas espectrales independientes se distribuyen entre trabajadores y entre casos. No se multiplica un grupo de procesos por otro grupo de hilos BLAS sin contarlos. Se restringen los hilos de las bibliotecas y se muestra progreso/ETA.',y)
    y=text('Para los controles cartesianos se preparan las mismas consultas exactas que pide el operador, sin redondear estados ni cambiar las derivadas de energía. Para la referencia radial, las frecuencias de Matsubara se resuelven independientemente en la misma cola. La agregación posterior conserva un orden definido.',y)
    y=text('<b>El 90 % es un techo de recursos, no una garantía de uso ni de aceleración.</b> También hay tareas de preparación y ensamblaje que no llenan todos los trabajadores. La memoria disponible limita el grupo, y la reserva se recalcula en cada ejecución. No se promete una aceleración de 28 veces.',y)
    receipt=DATA/'pilot/parallel_equivalence_receipt.json'
    if receipt.exists():
        p=read(receipt)
        y=text('<b>Piloto medido:</b> 19,56 s serial frente a 8,72 s con cuatro trabajadores: 2,24 veces más rápido, con resultados JSON y 33 matrices por caso exactamente iguales. La referencia radial completa usó 27 trabajadores y un coordinador: 518 tareas en 4,84 s. Es otra carga de trabajo; no se mezclan sus tiempos como una aceleración del mismo cálculo.',y)
    y=text('<b>Estado de la etapa:</b> la referencia física requiere revisar el cierre del núcleo antes del transiente. Siguen pendientes las capacidades dinámicas que permitan interpretar una trayectoria, la movilidad material y los límites de las tasas fonónicas. Este avance no promueve el modelo a producción.',y)
    text('Inventario: geminga_resources.json. Política persistente: AGENTS.md. Las corridas previstas de más de cinco minutos siguen siendo manuales, con comando en la libreta y en el chat. No se lanzan trabajos de fondo desde el agente.',y,small);c.showPage()

    start(4,'Observar el inicio del mismo sistema','Prioridad futura: latencia relativa de 775 frente a 1550 nm en la cinta de 80 nm de Korzh.')
    y=picture('observation_horizon.png',695,182)
    y=text('<b>Hasta el gatillo en Vout más un margen.</b> La figura es pedagógica, sin escala ni parámetros físicos: la curva corta es exactamente el principio de la larga. Se acorta sólo el horizonte de integración y observación; no se eliminan estados lentos ni se altera su influencia sobre el comienzo del pulso.',y)
    y=text('Se conservan ecuaciones, material, depósito, malla, bordes y circuito completo de la memoria, incluidas inductancias y capacitancias. No se acortan sus constantes de tiempo ni se cambia la precisión para adaptar el dispositivo a una ventana breve. El estado final se guarda para continuar después si corresponde.',y)
    y=text('El cruce se calcula a partir de tiempos aceptados, con umbral y polaridad comunes para ambos colores. Una confirmación evita contar oscilaciones, pero el tiempo de latencia sigue siendo el del cruce interpolado. Umbral, confirmación y margen se registrarán cuando se defina ese ensayo; no se inventan aquí.',y)
    y=text('<b>Si no hay cruce:</b> se detiene en un techo temporal finito y se registra que la ventana no observó el disparo. No se asigna una latencia ficticia ni se confunde un fallo numérico con no detección. El retardo óptico previo a la transferencia continúa explícitamente abierto.',y)
    text('El modo heredado latency espera un máximo confirmado; no implementa todavía gatillo + margen. Esta entrega documenta el comportamiento futuro, sin cambiar producción. Véanse horizon_policy.md/json. La recuperación total de nanosegundos se pospone hasta estudiar la respuesta inicial pertinente.',y,small)
    c.showPage()
    start(5,'Un límite físico ya identificado','Referencia térmica a 0,9 K; ℓ₀ = 4,70 nm. Mismo perfil prescrito para ambos cálculos.')
    y=picture('radial_force.png',695,217)
    y=text('<b>La fuerza del cierre actual no supera el contraste del núcleo.</b> La diferencia L² con Usadel espacial es 233 %, 158 % y 75 % para δ/Δ₀ = 0,05, 0,10 y 0,20. Con δ = 0,10 aparece una fuerza de -11,88 donde Usadel da +0,203: las direcciones de relajación son opuestas. La discrepancia sigue siendo 27,2 % fuera de r = ℓ₀ (comparación con cola estimada); no se limita al pico subnanométrico.',y)
    y=text('Duplicar las frecuencias de la referencia cambia su fuerza integrada 1,87 %; con la estimación explícita de cola, ese cambio baja a 0,080 %. Ampliar el radio exterior cambia sólo 0,000106 % en el núcleo comparado. Esos controles son mucho menores que la discrepancia del candidato; no se presentan como cotas rigurosas de error.',y)
    y=text('<b>Decisión:</b> no ajustar δ para disimular el pico ni repetir mallas costosas. El siguiente desarrollo debe corregir la respuesta espacial del cierre cerca de amplitud pequeña, contrastando cada propuesta con estos mismos campos espectrales. Quitar ese aporte y conservar sólo BCS más gradiente deja una discrepancia de 83,7 % frente a esa referencia con cola estimada. Se necesita revisar el cierre de energía, no reescalar la movilidad.',y)
    y=text('Este resultado corresponde a un vórtice circular auxiliar con amplitud prescrita. No resuelve una barrera, un vórtice autoconsistente, la cinta completa ni la movilidad. La etapa 4 permanece abierta; no se justifica todavía una corrida pesada con fotón.',y)
    text('Datos: raw/stage4_radial_reference_20260923. Análisis: radial_reference_analysis.md. Los términos de cola se conservan separados y los problemas espectrales comparten el presupuesto de recursos.',y,small)
    c.save();print(json.dumps(dict(pdf=str(PDF),pages=5)))


if __name__=='__main__':
    main()
