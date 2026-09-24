"""Three-page result report; figures identify fields, norms, units and cases.

This report reads certified saved fields only. It does not change the original
temporal gate or represent the diagnostic core as a photon detector transient.
"""
from pathlib import Path
import hashlib
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.font_manager import findfont
from reportlab.pdfgen import canvas
from reportlab.lib.colors import HexColor
from reportlab.lib.styles import ParagraphStyle
from reportlab.platypus import Paragraph, Table, TableStyle
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont

ROOT=Path(__file__).resolve().parents[2]
DATA=ROOT/'docs/implementation/stage4/practical_time_review_20260924'
FIG=DATA/'figures'
PDF=ROOT/'output/pdf/implementation/Informe_etapa_4_revision_practica.pdf'
BLUE,ORANGE,GREEN,RED,GREY='#1b6595','#c66f23','#22836e','#ac4149','#617281'
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9.2,
    'axes.spines.top':False,'axes.spines.right':False,'savefig.dpi':220})


def read(path):
    return json.loads(path.read_text(encoding='utf8'))


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def time_axis(ax):
    ax.set_xscale('symlog',linthresh=.001)
    ax.set_xticks([0,.001,.01,.1,1],['0','0,001','0,01','0,1','1'])
    ax.set_xlim(0,1.05)
    ax.set_xlabel('Tiempo físico t (ps)')
    ax.grid(alpha=.2)


def make_figures(data):
    FIG.mkdir(exist_ok=True,parents=True)
    times=[row['time_ps'] for row in data['responses']['refined']['amplitude']]
    fig,axes=plt.subplots(1,2,figsize=(7.5,3.2),layout='constrained')
    for probe,label,color in [('amplitude','Sonda de amplitud',BLUE),('angular_phase','Sonda angular',ORANGE)]:
        for run,style,marker in [('primary','--',None),('refined','-','o')]:
            rows=data['responses'][run][probe]
            y=100*np.array([row['displacement'] for row in rows])/rows[0]['displacement']
            axes[0].plot(times,y,style,color=color,marker=marker,markersize=3.5,
                markerfacecolor='white',label=label if run=='refined' else None)
        rows=[row for row in data['comparisons'] if row['probe']==probe and row['observable']=='displacement' and row['time_ps']>=.3]
        axes[1].plot([row['time_ps'] for row in rows],
            [1e6*row['relative_to_initial'] for row in rows],'o-',color=color,markersize=4,label=label)
    axes[0].set(title='Campo del condensado tras restar la base',ylabel=r'$100\,\|\delta d(t)\|_M/\|\delta d(0)\|_M$ (%)',ylim=(0,105))
    axes[0].legend(fontsize=8,loc='lower left');time_axis(axes[0])
    axes[1].set(title='Diferencia primaria - refinada',xlabel='Tiempo físico t (ps)',
        ylabel=r'$10^6\|\delta d_P-\delta d_R\|_M/\|\delta d(0)\|_M$ (ppm)',xlim=(.25,1.05),ylim=(0,3.6))
    axes[1].set_xticks([.3,1],['0,3','1']);axes[1].grid(alpha=.2)
    fig.suptitle('Núcleo térmico 65 × 65; cuadrado [−6, 6] ℓ₀; T = 0,9 K; 0–1 ps',fontsize=10)
    fig.savefig(FIG/'condensate_response_and_difference.png');plt.close(fig)

    fig,axes=plt.subplots(1,2,figsize=(7.5,3.15),layout='constrained')
    for run,label,style,color in [('primary','Primaria','--',ORANGE),('refined','Refinada','-',BLUE)]:
        rows=data['responses'][run]['amplitude']
        axes[0].plot(times,[row['phase_torque_density'] for row in rows],style,marker='o',
            color=color,markersize=3.5,markerfacecolor='white',label=label)
    axes[0].set_yscale('log');axes[0].set_ylim(1e-8,1e-5)
    axes[0].set(title='Torque de fase de la sonda de amplitud',ylabel=r'$\|\tau_A(t)-\tau_b(t)\|_M$ (adimensional)')
    axes[0].legend(fontsize=8);time_axis(axes[0])
    failed=data['failed_comparisons'][0]
    values=[failed['error_L2'],failed['registered_absolute_plus_relative_tolerance'],failed['remaining_norm']]
    bars=axes[1].bar(['Diferencia\nP − R','Límite\nregistrado','Señal\nrefinada'],np.array(values)*1e8,color=[RED,GREY,BLUE],width=.63)
    for bar,value in zip(bars,values):
        axes[1].text(bar.get_x()+bar.get_width()/2,bar.get_height()+.16,f'{value*1e8:.2f}'.replace('.',','),ha='center',fontsize=9)
    axes[1].set(title='Comparación que no cumple a 1 ps',ylabel='Norma del torque × 10⁸ (adimensional)',ylim=(0,9.8))
    axes[1].grid(axis='y',alpha=.2)
    fig.suptitle('Núcleo 65 × 65; cuadrado [−6, 6] ℓ₀; T = 0,9 K; 0–1 ps',fontsize=10)
    fig.savefig(FIG/'failed_cross_torque_scale.png');plt.close(fig)


def report(data):
    font=findfont('DejaVu Sans');bold=findfont('DejaVu Sans:weight=bold')
    pdfmetrics.registerFont(TTFont('DV',font));pdfmetrics.registerFont(TTFont('DVB',bold))
    pdfmetrics.registerFontFamily('DV',normal='DV',bold='DVB',italic='DV',boldItalic='DVB')
    W,H=595.276,841.89;L=39;WIDTH=W-2*L
    c=canvas.Canvas(str(PDF),pagesize=(W,H),pageCompression=1)
    c.setTitle('Etapa 4: revisión práctica del transiente térmico no lineal')
    c.setAuthor('pySNSPD - revisión de implementación')
    body=ParagraphStyle('body',fontName='DV',fontSize=9.5,leading=13,textColor=HexColor('#233647'))
    small=ParagraphStyle('small',parent=body,fontSize=8.2,leading=11)
    def p(text,y,style=body,x=L,width=WIDTH):
        q=Paragraph(text,style);_,height=q.wrap(width,1000);q.drawOn(c,x,y-height);return y-height
    def heading(number,title,subtitle):
        c.setFillColor(HexColor('#15576b'));c.rect(0,H-29,W,29,fill=1,stroke=0)
        c.setFillColor(HexColor('#ffffff'));c.setFont('DVB',9)
        c.drawString(L,H-19,'pySNSPD | ETAPA 4 | REVISIÓN PRÁCTICA')
        c.setFillColor(HexColor('#183449'));c.setFont('DVB',19);c.drawString(L,H-66,title)
        y=p(subtitle,H-79,small)
        c.setFillColor(HexColor('#677985'));c.setFont('DV',7.6)
        c.drawString(L,23,'24 septiembre 2026 | Evidencia conservada; producción intacta')
        c.drawRightString(W-L,23,f'{number} / 3')
        return y-14
    def image(name,y,height):
        c.drawImage(str(FIG/name),L,y-height,width=WIDTH,height=height,mask='auto');return y-height
    def table(rows,y,widths):
        converted=[[Paragraph(str(cell),small) for cell in row] for row in rows]
        t=Table(converted,colWidths=widths,hAlign='LEFT')
        t.setStyle(TableStyle([('BACKGROUND',(0,0),(-1,0),HexColor('#e8f0f3')),('VALIGN',(0,0),(-1,-1),'TOP'),
            ('LINEBELOW',(0,0),(-1,0),.6,HexColor('#a6bac4')),('BOTTOMPADDING',(0,0),(-1,-1),6),
            ('TOPPADDING',(0,0),(-1,-1),6),('LEFTPADDING',(0,0),(-1,-1),6),
            ('ROWBACKGROUNDS',(0,1),(-1,-1),[HexColor('#ffffff'),HexColor('#f6f8f9')])]))
        _,height=t.wrap(WIDTH,1000);t.drawOn(c,L,y-height);return y-height

    y=heading(1,'La integración llegó hasta 1 ps',
        'El proceso se detuvo al evaluar un criterio posterior. Las dos trayectorias y todos sus estados guardados están disponibles.')
    y=p('<b>Resultado:</b> 18,16 minutos de cálculo; 24 pasos primarios y 31 refinados, sin rechazos. '
        'El control final dio 95/96 comparaciones favorables. Sin embargo, 72 comparaban salidas idénticas porque los pasos coincidieron hasta 0,16 ps: '
        '<b>son 23/24 comparaciones favorables entre trayectorias realmente distintas</b>, observadas a 0,3 y 1 ps.',y)-10
    y=image('condensate_response_and_difference.png',y,221)-8
    y=p('<b>Figura 1.</b> d = Δ/(k<sub>B</sub>T<sub>c</sub>) es el condensado adimensional; δd = d<sub>sonda</sub> − d<sub>base</sub>. '
        'La norma espacial es ||z||<sub>M</sub> = √Σ m<sub>i</sub>|z<sub>i</sub>|², con áreas nodales adimensionales m<sub>i</sub>. '
        'Izquierda: respuesta restante respecto a su propia norma inicial. Líneas discontinuas: primaria; continuas y círculos: refinada. '
        'Se superponen porque concuerdan. Derecha: P y R denotan primaria y refinada; su diferencia se muestra en partes por millón de la respuesta inicial. '
        'No es una latencia ni una comparación con un experimento.',y,small)-11
    y=table([['Magnitud a 1 ps','Sonda de amplitud','Sonda angular'],
        ['Respuesta espacial ||δd|| / inicial','73,5503 %','1,73882 %'],
        ['Diferencia P − R / inicial: condensado','0,000315 %','0,000279 %'],
        ['Diferencia P − R / inicial: corriente','0,000363 %','0,000227 %'],
        ['Diferencia P − R / inicial: fuerza','0,000827 %','0,001251 %']],y,[WIDTH*.49,WIDTH*.255,WIDTH*.255])-12
    y=p('<b>Caso físico.</b> Núcleo con vórtice a T = 0,9 K, T<sub>c</sub> = 8,65 K; 256 frecuencias de Matsubara. '
        'Cuadrado x,y ∈ [−6,6]ℓ<sub>0</sub>, 65 × 65 nodos, contactos fijos y sin fotón. '
        'Sondas iniciales: δd<sub>A</sub> = 10<super>−3</super>d<sub>0</sub>b(r); '
        'δd<sub>F</sub> = 10<super>−3</super>i d<sub>0</sub>(x/R)b(r), donde b(r) = max(0,1−r²/R²)³ y R = 4ℓ<sub>0</sub>. '
        'ℓ<sub>0</sub> = √(ℏD/2k<sub>B</sub>T<sub>c</sub>) = 4,70 nm es la unidad de longitud. '
        'La base también evoluciona y se resta al mismo instante.',y,small)-8
    p('La acción térmica disminuye en las tres trayectorias. El defecto máximo del balance instantáneo entre su tasa y la disipación KWT + normal es '
        '1,74 × 10<super>−9</super> relativo a la disipación. Es energía libre a baño fijo; no acredita conservación de la energía interna del detector.',y,small)
    c.showPage()

    y=heading(2,'Un límite real en una cola pequeña',
        'El criterio original permanece incumplido. Se conserva ese resultado y se delimita qué magnitud no está resuelta.')
    y=p('La única comparación que no cumple es el <b>torque de fase inducido por la sonda de amplitud</b>, a 1 ps. '
        'Es una respuesta cruzada pequeña: cambiar la amplitud genera aquí mucho menos torque que perturbar directamente la fase. '
        'El error relativo de esa cola no sirve para rechazar toda la trayectoria, pero tampoco puede ocultarse.',y)-10
    y=image('failed_cross_torque_scale.png',y,218)-9
    y=p('<b>Figura 2.</b> G<sub>i</sub> es el gradiente de la acción térmica discreta respecto al condensado; '
        'τ<sub>i</sub> = Im(d<sub>i</sub>*G<sub>i</sub>)/m<sub>i</sub> es su fuerza para la fase, no un torque mecánico (* denota conjugación compleja). '
        'Izquierda: ||τ<sub>A</sub> − τ<sub>base</sub>||<sub>M</sub>, adimensional, en ambas resoluciones temporales. '
        'Derecha: norma de su diferencia P − R, límite registrado y norma restante refinada a 1 ps. '
        'El límite es 10<super>−8</super> + 0,005 max(norma inicial, norma restante). No se ha cambiado.',y,small)-12
    y=table([['Escala de la misma diferencia a 1 ps','Valor','Qué informa'],
        ['Norma absoluta P − R','8,305 × 10<super>−8</super>','Supera 2,204 veces el límite.'],
        ['Respecto a su torque inicial propio','1,500 %','No cumple la precisión registrada.'],
        ['Respecto al torque restante refinado','471,3 %','La cola tardía no está cuantificada.'],
        ['Respecto al torque inicial de la sonda angular','0,004566 %','Contextualiza su pequeño tamaño; no reemplaza el criterio.']],y,[WIDTH*.45,WIDTH*.19,WIDTH*.36])-12
    y=p('<b>Decisión práctica.</b> Se cierra el desarrollo de este ensayo térmico con una limitación explícita en la cola del torque cruzado. '
        'No se repite el lote para reducir esa componente. Las respuestas principales justifican avanzar al dominio y al método espacial finales; '
        'no certifican todavía la latencia, el hotbelt ni el detector.',y)-8
    p('El 79,55 % de la diferencia cuadrática de este torque está fuera del radio inicial de la sonda, r &gt; 4ℓ<sub>0</sub>. '
        'Es una localización observada, no prueba de que la geometría sea la causa ni de que otra malla eliminará la discrepancia. '
        'El estado histórico sigue siendo TEMPORAL_REFINEMENT_NOT_MET.',y,small)
    c.showPage()

    y=heading(3,'Siguiente avance: el método espacial final',
        'Menos controles sintéticos nuevos; más integración con los operadores y la geometría que usará el dispositivo.')
    y=p('La implementación continuará reutilizando la malla dual Delaunay-Voronoi y las estructuras de pyTDGL ya presentes en el repositorio. '
        'Su paso temporal es <b>Euler con resolución algebraica local de la amplitud, de primer orden</b>. '
        'La adaptación conserva la acción Usadel y la movilidad KWT; la publicación del método no valida automáticamente esta extensión.',y)-11
    y=p('<b>Avance concreto, sin otra corrida larga.</b> La interfaz con el paso KWT de la memoria pasó un control suave de 63 nodos y 16 frecuencias '
        'en 0,172 s. A Δt = 2,5 × 10<super>−5</super> ps, la velocidad estimada por un paso difiere <b>0,004697 %</b> del RHS de referencia, '
        'en norma de área. El error se divide aproximadamente por dos al dividir el paso por dos: es consistencia de primer orden. '
        'Este resultado verifica la conversión local de fuerza y tiempo; no una trayectoria ni el circuito.',y,small)-12
    y=table([['Trabajo siguiente','Qué se reutiliza o conserva','Evidencia necesaria'],
        ['Adaptar la acción al grafo dual','Áreas Voronoi, enlaces Delaunay, longitudes duales y contactos del Mesh existente.','Concordancia de geometría y operadores; balance discreto en un caso suave.'],
        ['Integración temporal más directa','Actualización KWT algebraica; misma fuerza microscópica, misma movilidad y mismas escalas físicas.','Un caso representativo y una comparación temporal útil; no una cadena de colas sintéticas.'],
        ['Completar el acoplamiento físico','Transporte cinético, devolución del calor de condensado y circuito completo de la memoria.','Balance del sistema acoplado y respuesta con bordes admitidos.'],
        ['Preparar la medición final','Cinta de 80 nm y ventana hasta gatillo V_out + margen.','Transferencia del fotón y definición del gatillo pendientes; sin acortar constantes del circuito.']],y,[WIDTH*.24,WIDTH*.38,WIDTH*.38])-17
    y=p('<b>Qué queda acreditado por esta corrida.</b> Evolución térmica no lineal estable en todo el campo durante 1 ps; '
        'Newton espectral funcionando; conservación de la identidad instantánea de disipación; buena coincidencia temporal de los observables principales. '
        'Los estados fueron revisados desde 71 archivos NPZ, con hashes de las fuentes y resultados. No se ejecutaron nuevas simulaciones para este diagnóstico.',y)-13
    y=p('<b>Qué sigue abierto.</b> La etapa 4 completa, el calor fuera de equilibrio, el transporte general y el circuito acoplado al transiente real. '
        'La malla dual puede mejorar la representación espacial, pero no sustituye la física faltante. '
        'Este control no reproduce el experimento de Korzh ni incluye un fotón.',y)-13
    y=p('<b>Regla para los próximos informes.</b> Cada figura identificará el campo u observable, la sonda y el dominio, la temperatura, las unidades, '
        'la resta de la base, la norma y su normalización. Una diferencia entre pasos temporales se distinguirá de una incertidumbre física y de una discrepancia experimental.',y)-15
    y=p('<b>Trazabilidad.</b> Corrida: stage4_nonlinear_time_20260924. Datos y extracción certificados: '
        'docs/implementation/stage4/practical_time_review_20260924/raw. Reconstrucción independiente: analysis.json. '
        'La tolerancia espectral, las escalas KWT y el criterio histórico no se han cambiado. Producción y v1.0.0 permanecen intactos.',y,small)-11
    y=p('Fuentes del método: <link href="https://py-tdgl.readthedocs.io/en/latest/background.html" color="#15576b">documentación oficial de pyTDGL</link>; '
        '<link href="https://github.com/loganbvh/py-tdgl/blob/main/tdgl/solver/solver.py" color="#15576b">implementación oficial del solver</link>; '
        '<link href="https://arxiv.org/abs/2302.03812" color="#15576b">artículo, arXiv:2302.03812</link>. '
        'Control nuevo: kwt_bridge/receipt.json.',y,small)-8
    p('Observación de implementación: el pase primario hizo un último paso de 1,11 × 10<super>−16</super> ps por redondeo del extremo temporal. '
        'El próximo controlador debe reconocer ese extremo sin un cálculo redundante. Los datos originales se conservan.',y,small)
    c.save()


def main():
    data=read(DATA/'analysis.json')
    make_figures(data);report(data)
    metadata=dict(pdf_sha256=sha(PDF),analysis_sha256=sha(DATA/'analysis.json'),
        figures={p.name:sha(p) for p in FIG.glob('*.png')},plot_definitions=data['plot_metadata'],
        geometry_verified_from_fields=dict(node_count=4225,x_bounds_ell0=[-6,6],y_bounds_ell0=[-6,6],probe_radius_ell0=4),
        original_gate_changed=False,new_simulations_for_etd_diagnosis=0,
        separate_lightweight_kwt_control_receipt_sha256=sha(DATA/'kwt_bridge/receipt.json'))
    (DATA/'report_context.json').write_text(json.dumps(metadata,indent=2,ensure_ascii=False)+'\n',encoding='utf8')
    print(json.dumps(dict(pdf=str(PDF),sha256=metadata['pdf_sha256'],figures=len(metadata['figures']))))


if __name__=='__main__':
    main()
