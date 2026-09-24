"""Result-focused report of resolved charge moments and the next dynamic test."""
from pathlib import Path
import json
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

ROOT=Path(__file__).resolve().parents[2]
DATA=ROOT/'docs/implementation/stage4/moment_review_20260924'
RAW=DATA/'raw'; FIG=DATA/'figures'
OLD=ROOT/'docs/implementation/stage4/self_consistent_review_20260924'
PDF=ROOT/'output/pdf/implementation/Informe_etapa_4_respuesta_de_carga.pdf'
COLORS=['#1766a4','#d17b2f','#288c72','#8e60ac']
plt.rcParams.update({'font.size':9,'axes.spines.top':False,'axes.spines.right':False,'savefig.dpi':190})

def read(p):return json.loads(p.read_text(encoding='utf8'))

def figures():
    FIG.mkdir(parents=True,exist_ok=True)
    analysis=read(RAW/'moment_comparisons.json')
    old=read(OLD/'potential_projection/summary.json')
    previous=next(row for row in old['records'] if row['case_id']=='radial_65_N256' and row['eta_relative']==.02)
    old_difference=previous['probes']['angular']['projection_errors']['charge_current']['relative_difference']
    current=[100*old_difference]
    for grid in ['base31','fine50']:
        m=analysis['moments']['radial_65_N256_eta0p0200_'+grid]['angular']['charge_current']
        current.append(100*m['projection_difference_D']/m['full_norm'])
    nodes=[12,31,50]
    chi=[100*(previous['chi_quadrature']-1)]
    chi += [100*(analysis['chi_quadrature']['radial_65_N256_eta0p0200_'+g]['numerical']-1) for g in ['base31','fine50']]
    fig,axes=plt.subplots(1,2,figsize=(10,3),layout='constrained')
    axes[0].semilogy(nodes,chi,'o-',color=COLORS[0]);axes[0].set(ylim=(.2,35),
        xlabel='Número de energías',ylabel='Exceso de la integral conocida (%)',title='La cuadratura mejora al resolver su escala')
    for x,y in zip(nodes,chi):axes[0].annotate(f'{y:.2f} %',(x,y),xytext=(0,8),textcoords='offset points',ha='center')
    axes[1].plot(nodes,current,'o-',color=COLORS[1]);axes[1].set(ylim=(0,24),
        xlabel='Número de energías',ylabel='Diferencia de corriente (%)',title='El potencial único conserva una discrepancia')
    for x,y in zip(nodes,current):axes[1].annotate(f'{y:.2f} %',(x,y),xytext=(0,8),textcoords='offset points',ha='center')
    for ax in axes:ax.set_xticks(nodes);ax.grid(alpha=.2)
    fig.savefig(FIG/'quadrature_and_current.png');plt.close(fig)
    budget=analysis['diagnostic_variation_budget']['angular']
    names=['charge_current','phase_torque_density','gap_force_density','energy_weighted_energy_flux']
    labels=['Corriente','Fuerza de fase','Fuerza compleja','Flujo de energía']
    differences=np.array([100*budget[k]['projection_difference_D']/budget[k]['full_norm'] for k in names])
    variation=np.array([100*budget[k]['observed_variation_U_D']/budget[k]['full_norm'] for k in names])
    fig,axes=plt.subplots(1,2,figsize=(10,3),layout='constrained')
    x=np.arange(len(names))
    axes[0].bar(x-.17,differences,width=.34,color=COLORS[1],label='Diferencia de representación')
    axes[0].bar(x+.17,variation,width=.34,color=COLORS[0],label='Variación numérica observada')
    axes[0].set(yscale='log',ylim=(.1,160),xticks=x,xticklabels=labels,
        ylabel='Porcentaje de cada respuesta de referencia',title='La corriente y la fase requieren atención')
    axes[0].tick_params(axis='x',labelsize=7);axes[0].legend(fontsize=6.5);axes[0].grid(axis='y',alpha=.2)
    for xi,di in zip(x,differences):axes[0].text(xi-.17,di*1.1,f'{di:.2f}',ha='center',fontsize=7)
    with np.load(RAW/'projection/fields/radial_65_N256_eta0p0100_fine50.npz') as fields:
        xy=fields['coordinates_bar'];coordinates=np.unique(xy[:,0]);ell_nm=4.698490897269843
        for state,color,label in [('full',COLORS[0],'Respuesta espectral'),('projected',COLORS[1],'Un potencial')]:
            values=fields['angular_'+state+'_phase_torque_density']
            # A y-directed section through x=0 crosses the angular probe torque lobes.
            select=(abs(xy[:,0])<1e-12)&(abs(xy[:,1])<=4)
            axes[1].plot(xy[select,1]*ell_nm,values[select],'-',color=color,label=label)
        axes[1].set(xlabel='y en el corte x = 0 (nm)',ylabel='Fuerza de fase (unidades del ensayo)',
            title='La misma perturbación produce otra respuesta')
        axes[1].legend(fontsize=7);axes[1].grid(alpha=.2)
    fig.savefig(FIG/'resolved_charge_error.png');plt.close(fig)
    rows=read(DATA/'charge_frequency/comparisons.json')
    fig,axes=plt.subplots(1,2,figsize=(10,3),layout='constrained')
    for ax,key,title in zip(axes,['current','phase_torque_density'],['Corriente','Fuerza de fase']):
        subset=sorted([row for row in rows if row['case_id']=='radial_65_N256'
            and row['eta_relative']==.01 and row['probe']=='angular'
            and row['moment']==key and row['nu']>0],key=lambda x:x['nu'])
        nu=[row['nu'] for row in subset]
        ax.axvspan(.0032,1.4,color='#edf0f3',label='Exploración del cierre aproximado')
        ax.loglog(nu,[100*row['potential_relative_error'] for row in subset],
            's--',color=COLORS[1],label='Comprimir a un potencial')
        ax.loglog(nu,[100*row['dynamic_change']/row['dynamic_full_norm'] for row in subset],
            'o-',color=COLORS[0],label='Usar respuesta espectral instantánea')
        ax.set(xlim=(.00007,1.4),ylim=(.001,150),xlabel='Frecuencia reducida ν = ω tD',
            ylabel='Diferencia relativa (%)',title=title)
        ax.grid(alpha=.2);ax.legend(fontsize=6,loc='lower right')
    fig.savefig(FIG/'frequency_response.png');plt.close(fig)

def main():
    figures();PDF.parent.mkdir(parents=True,exist_ok=True)
    pdfmetrics.registerFont(TTFont('MomentReport',findfont('DejaVu Sans')))
    pdfmetrics.registerFont(TTFont('MomentReportBold',findfont('DejaVu Sans:weight=bold')))
    pdfmetrics.registerFontFamily('MomentReport',normal='MomentReport',bold='MomentReportBold',italic='MomentReport',boldItalic='MomentReportBold')
    c=canvas.Canvas(str(PDF),pagesize=(595.276,841.89));c.setTitle('Etapa 4: respuesta de carga y siguiente ensayo')
    style=ParagraphStyle('body',fontName='MomentReport',fontSize=10.2,leading=15,textColor=HexColor('#253343'))
    small=ParagraphStyle('small',parent=style,fontSize=8.5,leading=12)
    def text(value,y,st=style):
        p=Paragraph(value,st);_,h=p.wrap(499,730);p.drawOn(c,48,y-h);return y-h-13
    def start(n,title,subtitle):
        c.setFillColor(HexColor('#1766a4'));c.rect(0,832,596,10,fill=1,stroke=0)
        c.setFont('MomentReport',9);c.drawString(48,799,'pySNSPD  |  ETAPA 4 SIN FOTÓN  |  24-09-2026 UTC')
        c.setFillColor(HexColor('#162c46'));c.setFont('MomentReportBold',20);c.drawString(48,758,title)
        text(subtitle,738,small)
        c.setFillColor(HexColor('#718096'));c.setFont('MomentReport',8)
        c.drawString(48,29,'Resultados numéricos. Producción y v1.0.0 sin cambios.');c.drawRightString(547,29,f'{n} / 3')
    def picture(name,y,h=200):
        c.drawImage(str(FIG/name),42,y-h,width=511,height=h,preserveAspectRatio=True,anchor='c',mask='auto');return y-h-16
    start(1,'La diferencia de corriente persiste','La nueva malla energética separa el error de integración de la respuesta que falta representar.')
    y=text('<b>La campaña terminó correctamente en 3 min 18 s.</b> Completó 181 consultas espectrales, 362 respuestas a dos sondas cinéticas y siete comparaciones integradas. Reutilizó los núcleos térmicos y ejecutó los trabajos independientes en paralelo.',696)
    y=picture('quadrature_and_current.png',y)
    y=text('El peso térmico describe cuánto contribuye cada energía a un pequeño desplazamiento de la distribución. Su integral exacta es conocida: el exceso numérico bajó de <b>18,02 % a 1,62 % y luego a 0,424 %</b>. La concentración de puntos a baja energía redujo el defecto identificado.',y)
    y=text('Sin embargo, la diferencia de corriente entre un único potencial y la respuesta espectral completa permanece cerca del <b>20 %</b>. La figura conserva la misma malla espacial 65² y el mismo contorno η = 0,02 del gap de referencia: el cambio de cuadratura no oculta un cambio de caso físico.',y)
    y=text('<b>Decisión:</b> no repetir el núcleo ni seguir refinando la DOS para intentar eliminar esta discrepancia. La respuesta de carga requiere una representación más rica o una reducción con su dominio temporal justificado.',y)
    text('El ensayo prescribe direcciones infinitesimales de población con gap, fase y contactos espectrales fijos. No representa un fotón, una trayectoria del detector ni un rechazo de toda la dinámica de fase y potencial de la memoria.',y,small);c.showPage()
    start(2,'La fase revela la limitación','Referencia: núcleo radial 65², cincuenta energías y η = 0,01 del gap. Sonda que rompe la simetría radial.')
    y=picture('resolved_charge_error.png',696)
    y=text('La diferencia es <b>20,55 % en corriente</b>, frente a una variación numérica observada de 2,87 %. En la fuerza de fase es <b>77,53 %</b>, frente a 3,64 %. Esas variaciones suman los cambios de cuadratura, contorno y malla; son referencias observadas, no cotas rigurosas de error.',y)
    y=text('La fuerza compleja total difiere sólo 3,23 %, porque domina su componente de amplitud. Esta última cambia aproximadamente 0,0015 %, por debajo de la variación numérica observada. Una comparación centrada sólo en amplitud o calentamiento habría ocultado la discrepancia de fase, relevante para la señal eléctrica.',y)
    y=text('El flujo de energía integrado difiere 0,438 %. El potencial único sí recupera parte de la corrección de corriente: reduce su error respecto de omitir la respuesta de carga desde aproximadamente 50 % hasta 20,55 %. Esa mejoría no lo convierte todavía en un cierre admitido para Vout.',y)
    text('Cada porcentaje se refiere a la norma espacial de su propia respuesta de prueba. El 77,53 % de fuerza de fase no es una predicción de error en latencia. El control radial casi simétrico tiene un efecto de carga demasiado pequeño para decidir por sí solo la reducción.',y,small);c.showPage()
    start(3,'Conservar el detalle energético','La eliminación algebraica conserva el detalle espectral y evita añadir un estado temporal por defecto.')
    y=picture('frequency_response.png',696,195)
    body=read(DATA/'report_content.json')
    for paragraph in body['next_paragraphs']:y=text(paragraph,y)
    text(body['footer'],y,small)
    c.save();print(json.dumps(dict(pdf=str(PDF),pages=3)))

if __name__=='__main__':main()
