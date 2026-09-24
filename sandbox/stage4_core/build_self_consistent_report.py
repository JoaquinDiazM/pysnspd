"""Illustrated results of the completed thermal core and its spectral follow-up."""
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
DATA=ROOT/'docs/implementation/stage4/self_consistent_review_20260924'
RAW=DATA/'raw'
FIG=DATA/'figures'
PDF=ROOT/'output/pdf/implementation/Informe_etapa_4_nucleo_autoconsistente.pdf'
ELL_NM=4.698490897269843
CASES=['radial_65_N128','radial_65_N256','radial_129_N256','asymmetric_65_N256']
COLORS=['#d17b2f','#8e60ac','#1766a4','#288c72']
LABELS=['65², 128 frecuencias','65², 256 frecuencias','129², 256 frecuencias','Asimétrico: 65², 256']
plt.rcParams.update({'font.size':9,'axes.spines.top':False,'axes.spines.right':False,'savefig.dpi':190})

def read(path):return json.loads(path.read_text(encoding='utf8'))

def figures():
    FIG.mkdir(parents=True,exist_ok=True)
    fig,axes=plt.subplots(1,2,figsize=(10,3.1),layout='constrained')
    for case,color,label in zip(CASES,COLORS,LABELS):
        history=[read(p) for p in sorted((RAW/'histories').glob(f'{case}_sweep*.json'))]
        n=np.asarray([v['completed_sweeps'] for v in history])
        residual=np.asarray([v['metrics']['core_mass_rms_relative'] for v in history])*100
        energy=np.asarray([v['metrics']['energy'] for v in history])
        axes[0].semilogy(n,residual,'o-',ms=3,color=color,label=label)
        axes[1].plot(n,energy-energy[0],'o-',ms=3,color=color,label=label)
    axes[0].axhline(.1,color='#283b4f',ls='--',label='Objetivo: 0,1 %')
    axes[0].set(xlabel='Barrido de optimización (sin tiempo físico)',ylabel='Residuo RMS del núcleo (%)',title='Los cuatro casos alcanzan el objetivo')
    axes[1].set(xlabel='Barrido de optimización (sin tiempo físico)',ylabel='Energía - inicial de cada caso (adimensional)',title='La energía disminuye durante el ajuste')
    for ax in axes:ax.grid(alpha=.2);ax.legend(fontsize=7)
    fig.savefig(FIG/'relaxation_history.png');plt.close(fig)
    identity=read(RAW/'identity.json');gap=identity['gap_reference']
    fig,axes=plt.subplots(1,2,figsize=(10,3.1),layout='constrained')
    radius=np.linspace(0,4,300)
    axes[0].plot(radius*ELL_NM,np.tanh(radius),'--',color='#767676',label='Perfil impuesto al iniciar')
    fine=None
    for case,color,label in zip(CASES[:3],COLORS[:3],LABELS[:3]):
        with np.load(RAW/f'{case}.npz') as a:
            xy=a['coordinates_bar'];ids=np.flatnonzero((abs(xy[:,1])<1e-12)&(xy[:,0]>=0)&(xy[:,0]<=4))
            x=xy[ids,0];v=abs(a['d'][ids])/gap
            axes[0].plot(x*ELL_NM,v,'o-',ms=2.5,color=color,label=label)
            if case=='radial_129_N256':fine=(x,v)
    with np.load(RAW/'radial_65_N256.npz') as a:
        xy=a['coordinates_bar'];ids=np.flatnonzero((abs(xy[:,1])<1e-12)&(xy[:,0]>=0)&(xy[:,0]<=4))
        coarse256=(xy[ids,0],abs(a['d'][ids])/gap)
    for case,color,label in zip(CASES[:2],COLORS[:2],['Corte: 128 - 256, malla 65²','Malla: 65² - 129², corte 256']):
        with np.load(RAW/f'{case}.npz') as a:
            xy=a['coordinates_bar'];ids=np.flatnonzero((abs(xy[:,1])<1e-12)&(xy[:,0]>=0)&(xy[:,0]<=4))
            x=xy[ids,0];v=abs(a['d'][ids])/gap
            reference=coarse256 if case=='radial_65_N128' else fine
            axes[1].plot(x*ELL_NM,100*(v-np.interp(x,*reference)),'o-',ms=3,color=color,label=label)
    axes[0].set(xlabel='Distancia al centro (nm)',ylabel='Amplitud / referencia homogénea',title='El núcleo se determina por la energía')
    axes[1].axhline(0,ls='--',color=COLORS[2],label='Diferencia cero')
    axes[1].set(xlabel='Distancia al centro (nm)',ylabel='Diferencia (% de amplitud homogénea)',title='Se muestra la diferencia entre curvas')
    for ax in axes:ax.grid(alpha=.2);ax.legend(fontsize=7)
    fig.savefig(FIG/'core_profiles.png');plt.close(fig)
    spectra=DATA/'retarded_oracle'
    with np.load(spectra/'radial_129_N256_DOS_maps.npz') as a:
        xy=a['coordinates_bar'];mask=np.isclose(a['eta_relative'],.02)
        energies=a['energy_relative'][mask];dos=a['DOS'][mask]
        fig,axes=plt.subplots(1,2,figsize=(10,3.1),layout='constrained')
        for radius,color in zip([0.,1.,2.,4.],['#283b4f','#1766a4','#d17b2f','#288c72']):
            idx=int(np.argmin(np.sum((xy-[radius,0])**2,axis=1)))
            actual=float(np.linalg.norm(xy[idx]))*ELL_NM
            axes[0].plot(energies,dos[:,idx],'o-',ms=3,color=color,label=f'r = {actual:.1f} nm')
        axes[0].set(xlabel='Energía / gap de referencia',ylabel='DOS electrónica / DOS normal',
            title='El núcleo cambia los estados disponibles',xlim=(0,2))
        axes[0].grid(alpha=.2);axes[0].legend(fontsize=7)
        field=dos[np.flatnonzero(np.isclose(energies,.5))[0]].reshape(129,129).T
        extent=np.r_[xy[:,0].min(),xy[:,0].max(),xy[:,1].min(),xy[:,1].max()]*ELL_NM
        image=axes[1].imshow(field,origin='lower',extent=extent,cmap='viridis',vmin=0,vmax=1)
        axes[1].set(xlabel='x (nm)',ylabel='y (nm)',title='E = 0,5 del gap de referencia')
        fig.colorbar(image,ax=axes[1],label='DOS electrónica / DOS normal',shrink=.86)
        fig.savefig(FIG/'electronic_DOS.png');plt.close(fig)
    projection=read(DATA/'potential_projection/summary.json')
    fig,axes=plt.subplots(1,2,figsize=(10,3.1),layout='constrained')
    names=['charge_current','gap_force_density','energy_weighted_energy_flux']
    labels=['Corriente','Fuerza','Flujo de energía']
    for k,(name,label,color) in enumerate(zip(names,labels,COLORS)):
        values=100*np.array([row['probes']['angular']['projection_errors'][name]['relative_difference']
            for row in projection['records']])
        axes[0].bar(k,values.mean(),color=color,width=.6)
        axes[0].plot([k,k],[values.min(),values.max()],color='#162c46',lw=2)
        axes[0].text(k,values.max()*1.13,f'{values.min():.2f}–{values.max():.2f} %',
            ha='center',fontsize=8)
    axes[0].set(xticks=range(3),xticklabels=labels,yscale='log',ylim=(.2,40),
        ylabel='Diferencia relativa (%)',
        title='Un potencial: respuesta angular integrada')
    axes[0].grid(axis='y',alpha=.2)
    temperature=projection['records'][0]['temperature_ratio']/gap
    energy=np.linspace(0,.5,500)
    def susceptibility(x):
        e=np.exp(-np.asarray(x)/temperature)
        return 2*e/(temperature*(1+e)**2)
    coarse=np.array([0.,.25,.5])
    axes[1].plot(energy,susceptibility(energy),color='#1766a4',label='Respuesta térmica exacta')
    axes[1].plot(coarse,susceptibility(coarse),'o--',color='#d17b2f',label='Puntos de la malla actual')
    axes[1].fill_between(coarse,susceptibility(coarse),alpha=.12,color='#d17b2f')
    axes[1].set(xlabel='Energía / gap de referencia',ylabel='Peso de respuesta electroquímica',
        title='Faltan puntos donde pesa la respuesta',xlim=(0,.5))
    axes[1].text(.18,6.2,'Integral en 0–5: 1,180\nValor exacto: 1',fontsize=9)
    axes[1].legend(fontsize=7);axes[1].grid(alpha=.2)
    fig.savefig(FIG/'potential_projection.png');plt.close(fig)

def main():
    figures();summary=read(RAW/'summary.json');PDF.parent.mkdir(parents=True,exist_ok=True)
    pdfmetrics.registerFont(TTFont('CoreReport',findfont('DejaVu Sans')))
    pdfmetrics.registerFont(TTFont('CoreReportBold',findfont('DejaVu Sans:weight=bold')))
    pdfmetrics.registerFontFamily('CoreReport',normal='CoreReport',bold='CoreReportBold',italic='CoreReport',boldItalic='CoreReportBold')
    c=canvas.Canvas(str(PDF),pagesize=(595.276,841.89));c.setTitle('Etapa 4: núcleo autoconsistente y siguiente prueba')
    style=ParagraphStyle('body',fontName='CoreReport',fontSize=10.2,leading=15,textColor=HexColor('#253343'))
    small=ParagraphStyle('small',parent=style,fontSize=8.5,leading=12)
    def text(value,y,st=style):
        p=Paragraph(value,st);_,h=p.wrap(499,730);p.drawOn(c,48,y-h);return y-h-13
    def start(n,title,subtitle):
        c.setFillColor(HexColor('#1766a4'));c.rect(0,832,596,10,fill=1,stroke=0)
        c.setFont('CoreReport',9);c.drawString(48,799,'pySNSPD  |  ETAPA 4 SIN FOTÓN  |  24-09-2026 UTC')
        c.setFillColor(HexColor('#162c46'));c.setFont('CoreReportBold',20);c.drawString(48,758,title)
        text(subtitle,738,small)
        c.setFillColor(HexColor('#718096'));c.setFont('CoreReport',8)
        c.drawString(48,29,'Resultados numéricos. Producción y v1.0.0 sin cambios.')
        c.drawRightString(547,29,f'{n} / 4')
    def picture(name,y,h):
        c.drawImage(str(FIG/name),42,y-h,width=511,height=h,preserveAspectRatio=True,anchor='c',mask='auto');return y-h-16
    start(1,'El núcleo ya es autoconsistente','Cuatro casos térmicos, con dos mallas y dos cortes de frecuencias. Sin fotón ni evolución temporal.')
    y=text('<b>La continuación terminó bien:</b> los cuatro casos cumplieron el objetivo en 15 a 17 barridos totales. La corrida completó 12.672 problemas espectrales en <b>8,08 minutos</b>; aprovechó los checkpoints del piloto anterior.',696)
    y=picture('relaxation_history.png',y,202)
    y=text('El residuo del núcleo quedó entre <b>0,084 % y 0,098 %</b>, frente al objetivo de 0,1 %. El mayor residuo nodal global fue 0,195 %. La norma del núcleo y el máximo global se informan por separado: un promedio pequeño no garantiza que todos los nodos tengan ese mismo residuo.',y)
    y=text('En cada barrido, el espectro responde al condensado y después se ajusta el condensado a ese espectro. Ambos usan la misma energía. El descenso observado concuerda con la identidad analítica de la actualización; las iteraciones no se interpretan como picosegundos de respuesta del dispositivo.',y)
    y=text('<b>Decisión:</b> aceptar esta referencia térmica autoconsistente y continuar hacia su descripción espectral a energías reales. No se requiere repetir la relajación con tolerancias más estrictas para dar ese paso.',y)
    text('Los resultados usan la misma energía de corte finito, los bordes espectrales prescritos y el condensado emparejado con sus propios espectros. Los datos completos permanecen en scratch; la entrega publica mapas compactos, historial y comprobaciones de integridad.',y,small);c.showPage()
    start(2,'Perfiles y sensibilidad observable','El mismo núcleo se compara en unidades físicas; ℓ₀ = 4,70 nm y T = 0,9 K.')
    y=picture('core_profiles.png',696,202)
    body=read(DATA/'report_content.json')
    for paragraph in body['profile_paragraphs']:y=text(paragraph,y)
    text('Comparaciones entre mallas: se usa densidad de corriente por ancho dual, no corriente integrada de un enlace. Comparaciones entre cortes: cada número de frecuencias tiene su propio condensado relajado. No se atribuye a estas diferencias una incertidumbre experimental del material.',y,small);c.showPage()
    start(3,'El espectro espacial ya está calculado','DOS electrónica local: estados electrónicos disponibles a cada energía. Distinta de la DOS fonónica del catálogo.')
    y=text('<b>Pasaron las 72 consultas:</b> tres núcleos, doce energías y dos valores del desplazamiento numérico η. La campaña terminó en <b>88,48 s</b>, con 27 trabajadores y un coordinador. Los 1.728 pasos internos de continuación mantuvieron la rama causal.',696)
    y=picture('electronic_DOS.png',y,202)
    y=text('En el centro del vórtice la DOS recupera el valor normal; al alejarse, cambian tanto la supresión de estados de baja energía como el máximo cercano al gap. La figura usa el núcleo radial de malla 129² y η = 0,02 del gap de referencia. Las líneas conectan las energías calculadas; no representan una resolución continua del pico.',y)
    y=text('<b>Correspondencia comprobada:</b> al volver al eje térmico, ambas representaciones recuperan la misma acción, con diferencia máxima 2,85 × 10<super>-13</super>. El cambio máximo de DOS entre mallas es 0,455 % en la norma del núcleo. Se mantiene el mismo entorno radial que suministra los bordes.',y)
    y=text('<b>Resolución espectral todavía explícita:</b> reducir η de 0,04 a 0,02 cambia la DOS hasta un 9,92 % en esa norma. η es un desplazamiento numérico del contorno, no una tasa material medida. Este resultado admite estudiar el siguiente operador cinético; todavía no acredita el límite η → 0 ni integrales precisas con sólo doce energías.',y)
    text('Cada espectro conserva el condensado autoconsistente calculado anteriormente. Sus archivos completos y hashes se mantienen en Geminga; la entrega incluye los mapas de DOS de las 72 consultas y campos complejos representativos.',y,small);c.showPage()
    start(4,'Qué simplificación permite el espectro','Pruebas estáticas de respuesta pequeña, con el condensado fijo. No son un transiente ni una preparación fotónica.')
    y=picture('potential_projection.png',696,202)
    for paragraph in body['next_paragraphs']:y=text(paragraph,y)
    text(body['footer'],y,small)
    c.save();print(json.dumps(dict(pdf=str(PDF),pages=4)))

if __name__=='__main__':main()
