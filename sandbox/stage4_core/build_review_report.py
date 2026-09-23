"""Build the stage4A results report from frozen raw fields and audit records."""
from pathlib import Path
import json
import sys

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.font_manager import findfont
from reportlab.pdfgen import canvas
from reportlab.lib.colors import HexColor
from reportlab.lib.styles import ParagraphStyle
from reportlab.platypus import Paragraph
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT/'docs/implementation/stage4/review_20260923'
FIG = DATA/'figures'
OUTPUT = ROOT/'output/pdf/implementation/Informe_avance_etapa_4A.pdf'
BLUE, ORANGE, GREEN = '#2267a9', '#cf6b2f', '#278378'
plt.rcParams.update({'font.size': 10, 'axes.spines.top': False,
                     'axes.spines.right': False, 'savefig.dpi': 170})


def read(path):
    return json.loads(path.read_text(encoding='utf8'))


def plots():
    audit = read(DATA/'independent_audit.json')
    boundary = read(DATA/'boundary_response_review.json')['cases']
    local = [read(p) for p in sorted((DATA/'raw/stage4A_local_20260923').glob('*/result.json'))]
    fig, axes = plt.subplots(1, 2, figsize=(10, 3.25), sharey=True, layout='constrained')
    controls = ['weak', 'current', 'small', 'high_gradient', 'normal', 'near_zero']
    controls = [x for x in controls if any(r['case']['control'] == x for r in local)]
    # Retain any additional registered names instead of silently dropping them.
    controls += sorted({r['case']['control'] for r in local}-set(controls))
    labels = {'weak':'Débil', 'current':'Corriente', 'small':'Amplitud\nbaja',
              'high_gradient':'Gradiente\nfuerte', 'normal':'Normal', 'near_zero':'Casi cero',
              'near_normal':'Casi\nnormal'}
    for ax, pop, title in zip(axes, ['vacuum','synthetic_fixed_count'], ['Vacío electrónico', 'Población sintética fija']):
        for shift, delta, color in zip([-.23,0,.23],[.05,.1,.2],[BLUE,ORANGE,GREEN]):
            values = [next(r['symbol_eigenvalues'][0] for r in local if r['case']['control']==x
                and r['population']==pop and r['delta_reg']==delta and r['case']['resolution']=='reference') for x in controls]
            ax.bar(np.arange(len(controls))+shift,values,width=.22,color=color,label=f'δ = {delta:g}')
        ax.axhline(0,color='#3c4653',lw=.8)
        ax.set_xticks(np.arange(len(controls)),[labels.get(x,x) for x in controls], fontsize=8)
        ax.set_title(title,fontsize=11)
        ax.grid(axis='y',alpha=.2)
    axes[0].set_ylabel('Menor autovalor de D.36 (adimensional)')
    axes[1].legend(fontsize=8,loc='lower right',ncol=3)
    fig.savefig(FIG/'local_stability.png');plt.close(fig)

    fig, axes = plt.subplots(1,3,figsize=(10,3.2),layout='constrained')
    for profile, marker, color, label in [('smooth','o',BLUE,'Suave'),('suppressed','s',ORANGE,'Centro suprimido')]:
        records=[boundary[f'2d_{profile}_d0.1_{mesh}']['mobility']['Korzh_reference'] for mesh in ['4x2','8x4']]
        for ax,key in zip(axes,['original_total_heat','original_boundary_heat','constrained_total_heat']):
            values=[r[key] for r in records]
            ax.plot([153,561],values,color=color,marker=marker,lw=2,label=label,
                    ls='--' if profile=='suppressed' else '-',
                    markerfacecolor='white' if profile=='suppressed' else color)
            for x,y in zip([153,561],values):
                ax.annotate(f'{y:.4g}',(x,y),xytext=(0,9 if profile=='suppressed' else -15),
                            textcoords='offset points',ha='center',fontsize=8,color=color)
            ax.set_xticks([153,561]);ax.set_xlabel('Nodos 2D');ax.grid(alpha=.2)
            ax.margins(x=.2,y=.3)
    axes[0].set_title('Descenso original: todos libres',fontsize=10)
    axes[1].set_title('Contribución de los bordes',fontsize=10)
    axes[2].set_title('Respuesta con bordes fijos',fontsize=10);axes[2].set_yscale('log')
    axes[0].set_ylabel('Potencia normalizada de condensado')
    axes[1].legend(fontsize=8)
    fig.savefig(FIG/'boundary_heat.png');plt.close(fig)

    fig, axes = plt.subplots(1,2,figsize=(10,3.7),layout='constrained')
    for mesh,color,marker,label in [('4x2',BLUE,'o','153 nodos'),('8x4',ORANGE,'s','561 nodos')]:
        path=DATA/'raw/stage4A_spatial_20260923'/f'2d_suppressed_d0.1_{mesh}'/'fields.npz'
        with np.load(path) as a:
            xy=a['coordinates_m']*1e9
            take=np.abs(xy[:,1])<1e-9
            indices=np.flatnonzero(take);indices=indices[np.argsort(xy[indices,0])]
            axes[0].semilogy(xy[indices,0]-80,a['gamma'][indices],marker=marker,ms=4,color=color,label=label)
    # Analytic field on a dense line, not a new spectral or physical solution.
    identity=read(DATA/'raw/stage4A_spatial_20260923/identity.json')
    ell=identity['material']['ell0_m'];xx=np.linspace(-80,80,1001)*1e-9/ell
    amp=.95-.91*np.exp(-(xx/2.5)**2)
    phase_dx=np.full_like(xx,.07);phase_dy=.005*np.exp(-(xx/10)**2)
    gap_ratio=np.pi/np.exp(np.euler_gamma)
    gamma=(amp**2/(amp**2+.1**2))**2*(phase_dx**2+phase_dy**2)/gap_ratio
    axes[0].semilogy(xx*ell*1e9,gamma,color='#222222',ls='--',lw=1.4,label='Derivada analítica del perfil')
    axes[0].set(xlim=(-30,30),xlabel='x respecto del centro (nm)',ylabel='Γ/Δ₀, corte y = 0',title='El gradiente del núcleo cambia con la malla')
    axes[0].legend(fontsize=8);axes[0].grid(alpha=.2)
    anchor=audit['analytic_center_anchor'] if 'analytic_center_anchor' in audit else next(v for v in audit.values() if isinstance(v,dict) and 'analytic_center' in v)
    exact=next(r['gamma'] for r in anchor['results'] if r['delta_reg']==.1)
    values=[audit['spatial'][f'2d_suppressed_d0.1_{mesh}']['center_gamma_discrete']/exact for mesh in ['4x2','8x4']]+[1]
    axes[1].bar(['153 nodos','561 nodos','Perfil exacto'],values,color=[BLUE,ORANGE,'#777777'])
    axes[1].set(yscale='log',ylabel='Γ central / Γ del perfil exacto',title='La derivada exacta permite medir el error')
    for i,v in enumerate(values):axes[1].text(i,v*1.15,f'{v:.3g} ×',ha='center')
    axes[1].set_ylim(.6,220);axes[1].grid(axis='y',alpha=.2)
    fig.savefig(FIG/'core_resolution.png');plt.close(fig)
    ref=read(DATA/'linear_usadel_reference.json')['results'][1:]
    fig,axes=plt.subplots(1,2,figsize=(10,3.25),layout='constrained')
    k=[r['k_ell0'] for r in ref]
    axes[0].plot(k,[r['H_U'] for r in ref],color=BLUE,marker='o',label='Usadel espacial linealizado')
    axes[0].plot(k,[r['H_K0'] for r in ref],color=ORANGE,marker='s',ls='--',label='Gradiente local K₀')
    axes[0].set(xlabel='Número de onda kℓ₀',ylabel='Rigidez térmica normalizada',title='Mismo estado térmico; distinta respuesta espacial')
    axes[0].legend(fontsize=8);axes[0].grid(alpha=.2)
    excess=[100*r['stiffness_excess_fraction'] for r in ref]
    axes[1].bar(np.arange(len(ref)),excess,color=ORANGE)
    axes[1].set_xticks(np.arange(len(ref)),[f"{r['spatial_wavelength_nm']:.1f}" for r in ref])
    axes[1].set(xlabel='Longitud de onda 2π/k (nm)',ylabel='Exceso de rigidez local (%)',title='Error físico que la malla no elimina',ylim=(0,37))
    for i,v in enumerate(excess):axes[1].text(i,v+.8,f'{v:.2f} %',ha='center',fontsize=9)
    axes[1].grid(axis='y',alpha=.2)
    fig.savefig(FIG/'linear_usadel_reference.png');plt.close(fig)
    return audit,boundary


def build():
    FIG.mkdir(parents=True,exist_ok=True);OUTPUT.parent.mkdir(parents=True,exist_ok=True)
    audit,boundary=plots()
    pdfmetrics.registerFont(TTFont('Report',findfont('DejaVu Sans')))
    pdfmetrics.registerFont(TTFont('ReportBold',findfont('DejaVu Sans:weight=bold')))
    pdfmetrics.registerFontFamily('Report',normal='Report',bold='ReportBold',italic='Report',boldItalic='ReportBold')
    c=canvas.Canvas(str(OUTPUT),pagesize=(595.276,841.89))
    c.setTitle('Etapa 4A: resultados y continuación enfocada')
    style=ParagraphStyle('body',fontName='Report',fontSize=10.2,leading=15,textColor=HexColor('#253343'))
    small=ParagraphStyle('small',parent=style,fontSize=8.5,leading=12)
    def text(t,y,st=style):
        p=Paragraph(t,st);_,h=p.wrap(499,720);p.drawOn(c,48,y-h);return y-h-13
    def start(number,title,subtitle):
        c.setFillColor(HexColor('#2267a9'));c.rect(0,832,596,10,fill=1,stroke=0)
        c.setFont('Report',9);c.drawString(48,799,'pySNSPD  |  ETAPA 4 SIN FOTÓN  |  23 SEPTIEMBRE 2026')
        c.setFillColor(HexColor('#162c46'));c.setFont('ReportBold',22);c.drawString(48,758,title)
        text(subtitle,738,small)
        c.setFillColor(HexColor('#718096'));c.setFont('Report',8)
        c.drawString(48,29,'Controles estáticos e instantáneos. No predicción de un detector.')
        c.drawRightString(547,29,f'{number} / 5')
    def figure(name,y,height):
        c.drawImage(str(FIG/name),42,y-height,width=511,height=height,preserveAspectRatio=True,anchor='c',mask='auto')
        return y-height-16
    start(1,'46 controles completados','Informe de avance. La etapa 4 continúa abierta; los resultados permiten reducir y orientar la próxima corrida.')
    y=text('<b>40 controles locales en 27,2 s y seis estados 2D en 31,3 min.</b> Los archivos y las fuentes coinciden con el plan ejecutado. No hubo transientes ni fotones. La población electrónica espacial es sintética, fijada por número de estados, y el rectángulo de 160 × 80 nm es un ensayo prescrito.',695)
    y=text('<b>Lo aceptado:</b> energía, fuerza y corriente conjugada son coherentes en los controles realizados. Las variaciones espectrales cambian la fuerza relativa como máximo 3,2 × 10<super>-8</super>; no se repite esa batería. Las identidades de potencia KWT cierran a precisión numérica.',y)
    y=figure('local_stability.png',y,198)
    y=text('<b>36 resultados locales positivos y cuatro negativos.</b> La figura muestra los 36 casos base; las otras cuatro anclas repiten la resolución espectral. D.36 mide la rigidez de variaciones espaciales cortas: un autovalor negativo limita el dominio de evolución. Los cuatro negativos corresponden a gradiente fuerte con δ = 0,05 o 0,1, en vacío y con población sintética; no son trayectorias fallidas.',y)
    y=text('<b>Decisión:</b> conservar esas exclusiones locales. No aumentar δ por conveniencia ni declarar que todo núcleo físico está validado. El problema espacial observado requiere separar bordes y resolución antes de interpretar disipación.',y)
    text('Fuentes de resultados: raw/stage4A_local_20260923, raw/stage4A_spatial_20260923 e independent_audit.json. δ es el regularizador adimensional de la fase en zonas donde la amplitud del condensado se aproxima a cero.',y,small)
    c.showPage()

    start(2,'La disipación aparente del borde','La reacción que sostiene un campo impuesto no debe interpretarse como movimiento libre del condensado.')
    y=text('El cálculo original dejó moverse todos los nodos, incluidos los bordes de un perfil con fase impuesta. Al refinar, la fuerza de borde se divide por una masa nodal menor y genera una velocidad artificialmente grande. Esto explica que el calor total casi se duplique.',695)
    y=figure('boundary_heat.png',y,208)
    y=text('<b>Perfil suave, movilidad Korzh:</b> el calor total original pasa de 4,234 a 8,462. El interior pasa de 0,006724 a 0,006825: <b>+1,50 %</b>. En la malla fina, el 99,92 % del calor original estaba en los nodos exteriores. La norma volumétrica de fuerza interior cambia sólo 0,051 %.',y)
    y=text('<b>Corrección aplicada al diagnóstico:</b> se añade la carga conjugada que mantiene los cuatro lados fijos. Su velocidad y trabajo son cero; el interior conserva su respuesta y el balance energía-calor sigue cerrando. Se reprocesaron los seis estados sin consultar nuevos espectros. Esta condición sirve para comparar el ensayo prescrito; no representa los contactos ni los bordes aislantes del nanohilo real.',y)
    y=text('<b>Cómo leer los resultados:</b> las potencias están normalizadas con la misma escala en ambas mallas. El tercer panel usa escala logarítmica para mostrar también el perfil suave. Una menor potencia instantánea no selecciona por sí sola la movilidad física correcta.',y)
    text('Se retiran como indicadores de convergencia la norma euclídea de fuerzas integradas y la corriente máxima de una arista: ambas cambian al repartirse el volumen o la sección entre más nodos. La próxima salida incorpora fuerza volumétrica y corriente sumada en cortes transversales.',y,small)
    c.showPage()

    start(3,'El centro aún está poco resuelto','La energía total puede parecer convergente mientras las derivadas en una región pequeña siguen cambiando.')
    y=text('Para el perfil con amplitud central 0,04, la energía cambia un 0,20 % al pasar de 153 a 561 nodos, pero la disipación interior cambia <b>-28,9 %</b> y la fuerza interior un -7,98 %. La región estrecha del centro necesita más resolución.',695)
    y=figure('core_resolution.png',y,235)
    y=text('<b>No confundir un error de malla con una falla constitutiva.</b> Con δ = 0,05, la malla gruesa da un único autovalor negativo, -0,785, en el centro. Al usar la derivada analítica del mismo perfil, el menor autovalor es +1,571 para los tres δ. El Γ central grueso es unas 110 veces el exacto; con δ = 0,1, la malla fina reduce esa razón a 2,65, pero aún no llega a 1.',y)
    y=text('<b>Sólo dos estados nuevos:</b> δ = 0,05 con 561 nodos, para comprobar el signo con una mejor derivada; y δ = 0,1 con 2145 nodos, para comprobar la resolución del núcleo. Se mantienen población, perfil, material y resolución espectral; los cuatro bordes quedan fijos en el diagnóstico de movilidad.',y)
    text('Coste previsto: 40-65 min, un proceso y un hilo; reservar 2 GB RAM. La estimación usa los 518 s medidos para 561 nodos y el crecimiento del número de consultas. No se añade una batería nueva de tolerancias: se comparan signo, error de derivada y magnitud del cambio en la respuesta interior.',y,small)
    c.showPage()

    start(4,'Qué solver se está usando','La implementación de producción y el banco experimental tienen discretizaciones espaciales diferentes.')
    y=text('<b>Producción conserva la malla dual Delaunay-Voronoi.</b> Los operadores adaptados de pyTDGL usan volúmenes finitos y flujos por aristas, con fases de enlace para la derivada covariante. No se ha sustituido ese núcleo en producción.',695)
    y=text('<b>El condensado usa Euler adaptativo, de primer orden.</b> pyTDGL lo denomina Euler implícito: resuelve localmente la dependencia de la amplitud nueva, mientras evalúa el término espacial con el estado anterior. Adaptar el paso no convierte el método en uno de segundo orden. El circuito sí contiene un paso RK2 de punto medio; eso no hace de segundo orden a todo el acoplamiento.',y)
    y=text('<b>Los ensayos 4A usan elementos nodales GLL en un rectángulo 2D.</b> GLL significa Gauss-Lobatto-Legendre: los nodos y pesos permiten integrar y derivar el mismo funcional discreto. Aquí se evalúan campos estáticos y respuestas instantáneas. No se avanzó el tiempo ni se eligió el futuro integrador del nuevo modelo.',y)
    y=text('<b>Por qué importa para la siguiente etapa.</b> Una formulación física y un método de malla son decisiones distintas. Conservar energía, corriente y condiciones de borde debe verificarse al pasar al solver que se vaya a usar en el dispositivo. Estos controles GLL son un banco experimental; no autorizan por sí solos su promoción a producción.',y)
    y=text('<b>Coste observado:</b> prácticamente todo el tiempo de los casos finos se consume en consultas espectrales directas. Repetir ese coste en cada evaluación de un transiente sería poco eficiente. Antes de trayectorias extensas corresponde preparar reutilización o interpolación de energía y derivadas en el dominio realmente visitado, sin desacoplar fuerza y corriente de la energía común.',y)
    text('Referencia externa: documentación oficial de pyTDGL, Theoretical Background, secciones Finite volume method e Implicit Euler method; py-tdgl.readthedocs.io/en/latest/background.html. Procedencia y rutas exactas del código: solver_scope.md. El circuito previsto sigue siendo el de la memoria.',y,small)
    c.showPage()

    start(5,'Qué falta para cerrar la etapa','La comparación física y la precisión numérica responden preguntas diferentes.')
    reference=DATA/'linear_usadel_reference.json'
    if reference.exists():
        # The independently computed figure is produced by its own reference script.
        external=DATA/'figures/linear_usadel_reference.png'
        if external.exists():
            y=figure('linear_usadel_reference.png',695,245)
        else:y=695
        y=text('Se recalculó una referencia espacial independiente de Usadel linealizado, a temperatura impuesta de 0,9 K y T<sub>c</sub> = 8,65 K, sin corriente. Compara la rigidez frente a una pequeña variación sinusoidal de amplitud con la aproximación local de gradiente del modelo. Es un control de física espacial, no un transiente a ocupación fija ni una validación de vórtices.',y)
        y=text('<b>Exceso de rigidez local: 3,03 % a 98,4 nm y 30,88 % a 29,5 nm.</b> Esta diferencia no se elimina refinando la malla. Con idéntica movilidad, el segundo caso daría un tiempo de relajación 23,6 % menor; no es un error demostrado en la latencia del detector. La escala material usada es ℓ₀ = 4,6985 nm.',y)
    else:y=695
    y=text('<b>Después de los dos estados nuevos:</b> si el núcleo y su respuesta interior quedan suficientemente resueltos, continuar con un transiente débil sin fotón y condiciones de borde compatibles. Si el signo negativo persiste con derivadas resueltas, revisar el dominio constitutivo antes de evolucionarlo. No se continúa refinando a ciegas.',y)
    y=text('<b>Alcance aún abierto:</b> referencia física de núcleo, evolución temporal y las dependencias dinámicas pendientes de la etapa 3. El ensayo actual no mide latencia, hotbelt, jitter ni un pulso de lectura. La transferencia fotónica, su ancho y las tasas NbN absolutas siguen sin admitirse.',y)
    text('Reproducción: plan next_campaign_plan.json y libreta GEMINGA_COMMANDS.md. El corredor muestra barras y ETA, guarda fuentes y mapas, y rechaza sobrescribir una corrida. Los resultados originales y las fuentes de su versión se preservan íntegros.',y,small)
    c.save()
    print(json.dumps(dict(pdf=str(OUTPUT),pages=5,figures=str(FIG))))


if __name__ == '__main__':
    build()
