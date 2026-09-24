"""Illustrated results report of the saved weak thermal trajectory.

No physical solve. Plots preserve initial normalization and show separate
refinement differences rather than hiding almost-coincident trajectories.
"""
from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.font_manager import findfont
from matplotlib.lines import Line2D
from reportlab.pdfgen import canvas
from reportlab.lib.colors import HexColor
from reportlab.lib.styles import ParagraphStyle
from reportlab.platypus import Paragraph, Table, TableStyle
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont

ROOT=Path(__file__).resolve().parents[2]
DATA=ROOT/'docs/implementation/stage4/time_review_20260924'
RAW=DATA/'raw'
FIG=DATA/'figures'
PDF=ROOT/'output/pdf/implementation/Informe_etapa_4_evolucion_termica.pdf'
BLUE,ORANGE,GREEN,PURPLE,GREY='#1766a4','#d17b2f','#288c72','#8e60ac','#76818c'
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9.2,
    'axes.spines.top':False,'axes.spines.right':False,'savefig.dpi':210})


def read(path):
    return json.loads(path.read_text(encoding='utf8'))


def time_axis(ax, zero=True):
    ax.set_xscale('symlog' if zero else 'log', **({'linthresh':.001} if zero else {}))
    ax.set_xticks([0,.001,.01,.1,1] if zero else [.001,.01,.1,1])
    ax.set_xticklabels(['0','0,001','0,01','0,1','1'] if zero else ['0,001','0,01','0,1','1'])
    ax.set_xlabel('Tiempo (ps)')
    ax.set_xlim(0 if zero else .0008,1.05)
    ax.grid(alpha=.2)


def figures(analysis, arrays):
    FIG.mkdir(parents=True,exist_ok=True)
    times=np.asarray(analysis['times_ps'])
    fig,axes=plt.subplots(1,2,figsize=(8,3.15),layout='constrained')
    for probe,color,label in [('amplitude',BLUE,'Amplitud'),('angular_phase',ORANGE,'Fase angular')]:
        for name,style in [('primary','--'),('refined','-')]:
            rows=analysis['responses'][name][probe]
            values=np.array([r['displacement'] for r in rows])/rows[0]['displacement']*100
            axes[0].plot(times,values,style,color=color,lw=1.5,
                marker='o' if name=='refined' else None,markersize=4,
                markerfacecolor='white',label=label if name=='refined' else None)
        initial=analysis['responses']['refined'][probe][0]['displacement']
        base=[r['displacement_L2']/initial*100 for r in analysis['global_state_metrics']
            if r['pass_name']=='refined' and r['state']=='baseline']
        axes[1].plot(times,base,'o-',color=color,markersize=3,label='Inicial de '+label.lower())
    axes[0].set(ylim=(0,105),ylabel='Norma / su valor inicial (%)',title='Respuesta tras restar la base')
    axes[0].legend(fontsize=8,loc='lower left')
    axes[1].set(ylabel='Deriva / perturbación inicial (%)',title='La base residual también se mueve')
    axes[1].legend(fontsize=7.5,loc='upper left')
    for ax in axes:time_axis(ax)
    fig.savefig(FIG/'response_and_baseline.png');plt.close(fig)

    d0=arrays['d0'];xy=arrays['coordinates_bar']
    unit=np.divide(d0,abs(d0),out=np.zeros_like(d0),where=abs(d0)>1e-12)
    transverse=np.imag(np.conj(unit)[None,:]*arrays['refined_angular_phase_difference'])
    gap=read(RAW/'executed_plan.json')['gap_reference_kBTc']
    transverse=100*transverse/gap
    transverse[:,abs(d0)<=1e-12]=np.nan
    xs,ys=np.unique(xy[:,0]),np.unique(xy[:,1])
    xi,yi=np.searchsorted(xs,xy[:,0]),np.searchsorted(ys,xy[:,1])
    indices=[0,3,7]
    limit=np.nanmax(abs(transverse[indices]))
    fig,axes=plt.subplots(1,3,figsize=(8,2.65),layout='constrained',sharex=True,sharey=True)
    cmap=plt.get_cmap('RdBu_r').copy();cmap.set_bad('#aab1b8')
    for ax,i in zip(axes,indices):
        grid=np.full((len(ys),len(xs)),np.nan);grid[yi,xi]=transverse[i]
        im=ax.pcolormesh(xs,ys,grid,cmap=cmap,vmin=-limit,vmax=limit,shading='nearest',rasterized=True)
        ax.set(xlim=(-4,4),ylim=(-4,4),aspect='equal',xlabel='x / ℓ₀',title=f'{times[i]:g} ps'.replace('.',','))
        ax.set_xticks([-4,0,4]);ax.set_yticks([-4,0,4])
    axes[0].set_ylabel('y / ℓ₀')
    bar=fig.colorbar(im,ax=axes,shrink=.82,pad=.025)
    bar.set_label('Componente transversal / Δref (%)',fontsize=8)
    fig.savefig(FIG/'transverse_phase_maps.png');plt.close(fig)

    labels={'displacement':'Desplazamiento','current':'Corriente',
        'force_density':'Fuerza','phase_torque_density':'Torque de fase'}
    colors=dict(zip(labels,[BLUE,GREEN,PURPLE,ORANGE]))
    fig,axes=plt.subplots(1,2,figsize=(8,2.85),layout='constrained')
    for kind,label in labels.items():
        rows=[r for r in analysis['refinement'] if r['probe']=='angular_phase' and r['observable']==kind and r['time_ps']>0]
        t=[r['time_ps'] for r in rows]
        for ax,denom in zip(axes,['relative_to_initial','relative_to_response']):
            values=np.array([r[denom] for r in rows])*100
            ax.loglog(t,np.maximum(values,1e-12),'o-',color=colors[kind],markersize=3,label=label)
    axes[0].axhline(.5,color=GREY,ls=':',lw=1)
    axes[0].text(.0011,.57,'Escala relativa registrada: 0,5 %',fontsize=7,color=GREY)
    axes[0].set(ylim=(1e-8,1.4),ylabel='Diferencia / señal inicial (%)',title='La refinación cumple su contrato')
    axes[1].set(ylim=(1e-8,400),ylabel='Diferencia / señal restante (%)',title='Las colas tienen otra precisión')
    axes[1].legend(fontsize=6.8,loc='upper left')
    for ax in axes:time_axis(ax,False)
    fig.savefig(FIG/'refinement_two_scales.png');plt.close(fig)

    fig,axes=plt.subplots(2,2,figsize=(8,5.2),layout='constrained')
    for state,color,label in [('baseline',GREY,'Base residual'),('amplitude',BLUE,'Base + amplitud'),('angular_phase',ORANGE,'Base + fase')]:
        E=np.array(analysis['energy']['refined_'+state]['free_energy'])
        axes[0,0].plot(times,1000*(E[0]-E),'o-',color=color,markersize=3,
            markerfacecolor='white' if state=='angular_phase' else color,label=label)
        rows=[r for r in analysis['global_state_metrics'] if r['pass_name']=='refined' and r['state']==state]
        axes[1,0].plot(times,[100*r['constitutive_defect_relative_to_total_rhs'] for r in rows],'o-',color=color,markersize=3)
        axes[1,1].plot(times,[100*r['instantaneous_balance_defect_relative_to_losses'] for r in rows],'o-',color=color,markersize=3)
        if state=='baseline':
            rate=np.array([-r['quadratic_free_energy_rate'] for r in rows])
            residual=np.array([r['quadratic_rate_plus_approximate_losses'] for r in rows])
            axes[0,1].plot(times,1000*rate,'-',color=GREEN,lw=2,label='Energía liberada')
            axes[0,1].plot(times,1000*(rate+residual),'o',color=PURPLE,markersize=4,markerfacecolor='none',label='Disipación KWT + normal')
    axes[0,0].set(title='La energía cuadrática disminuye',ylabel='Energía liberada ×10³')
    axes[0,0].legend(fontsize=6.7)
    axes[0,1].set(title='Tasas en la base residual',ylabel='Tasa ×10³ (u. del ensayo)')
    axes[0,1].legend(fontsize=6.7)
    axes[1,0].set(title='Truncación de la ley de velocidad',ylabel='Defecto / velocidad total (%)')
    axes[1,1].set(title='Desajuste instantáneo del balance',ylabel='Desajuste / disipación (%)')
    for ax in axes.ravel():time_axis(ax)
    fig.savefig(FIG/'energy_and_taylor_limits.png');plt.close(fig)

    receipt=read(RAW/'extraction_receipt.json')
    refined=[r for r in receipt['accepted_steps'] if r['pass_name']=='refined']
    fig,axes=plt.subplots(1,2,figsize=(8,2.5),layout='constrained')
    axes[0].bar(['Resta directa','Diferencia estable'],[14.21,-3.153],color=[ORANGE,GREEN],width=.55)
    axes[0].axhline(0,color=GREY,lw=.8)
    axes[0].set(ylabel='Cambio de energía ×10¹⁵',title='Newton: recuperar el signo correcto')
    axes[0].text(0,15.2,'+1,421 × 10⁻¹⁴',ha='center',fontsize=7)
    axes[0].text(1,-5.6,'−3,153 × 10⁻¹⁵',ha='center',fontsize=7)
    axes[0].set_ylim(-8,20)
    axes[1].plot([r['time_ps'] for r in refined],[1e6*r['affine_constant_defect'] for r in refined],'o-',color=ORANGE,markersize=3,label='Propagación guardada')
    axes[1].axhline(0,color=GREY,ls='--',label='Constante exacta')
    axes[1].set(ylabel='Desviación de la constante (ppm)',title='Forzamiento afín: evitar su deriva')
    time_axis(axes[1]);axes[1].legend(fontsize=7)
    fig.savefig(FIG/'numerical_repairs.png');plt.close(fig)

    exact=read(DATA/'nonlinear_snapshots/summary.json')
    observables=['rhs','force_density','current','phase_torque_density']
    labels=['Velocidad','Fuerza','Corriente','Torque de fase']
    x=np.arange(len(observables))
    fig,axes=plt.subplots(1,2,figsize=(8,3.1),layout='constrained')
    for probe,color,shift,label in [('amplitude',BLUE,-.17,'Sonda de amplitud'),('angular_phase',ORANGE,.17,'Sonda angular')]:
        rows={r['observable']:r for r in exact['comparisons'] if r['probe']==probe and r['time_ps']==1.}
        for ax,key in zip(axes,['relative_to_initial','relative_to_response']):
            values=[100*rows[name][key] for name in observables]
            ax.bar(x+shift,values,width=.32,color=color,label=label)
    axes[0].axhline(1.,color=GREY,ls=':',lw=1)
    axes[0].text(-.45,1.17,'Margen registrado: 1 %',fontsize=7,color=GREY)
    axes[0].set(ylim=(.001,2.2),title='Diferencia / respuesta inicial',ylabel='Diferencia a 1 ps (%)')
    axes[1].set(ylim=(.01,250),title='Diferencia / respuesta restante',ylabel='Diferencia a 1 ps (%)')
    for ax in axes:
        ax.set_yscale('log');ax.set_xticks(x,labels);ax.tick_params(axis='x',labelsize=7.4)
        ax.grid(axis='y',alpha=.2)
    axes[1].legend(fontsize=7,loc='upper left')
    fig.savefig(FIG/'nonlinear_snapshot_differences.png');plt.close(fig)


def main():
    analysis=read(DATA/'analysis.json')
    context=read(DATA/'report_context.json')
    with np.load(RAW/'trajectory_compact.npz') as arrays:figures(analysis,arrays)
    pdfmetrics.registerFont(TTFont('TimeReport',findfont('DejaVu Sans')))
    pdfmetrics.registerFont(TTFont('TimeReportBold',findfont('DejaVu Sans:weight=bold')))
    pdfmetrics.registerFontFamily('TimeReport',normal='TimeReport',bold='TimeReportBold',italic='TimeReport',boldItalic='TimeReportBold')
    PDF.parent.mkdir(parents=True,exist_ok=True)
    c=canvas.Canvas(str(PDF),pagesize=(595.276,841.89))
    c.setTitle('Etapa 4: evolución térmica y límites numéricos')
    c.setAuthor('pySNSPD')
    style=ParagraphStyle('body',fontName='TimeReport',fontSize=10,leading=14.6,textColor=HexColor('#253343'))
    small=ParagraphStyle('small',parent=style,fontSize=8.3,leading=11.6)
    def text(value,y,st=style):
        p=Paragraph(value,st);_,h=p.wrap(499,750)
        if y-h<51:raise ValueError(f'Page overflow: bottom={y-h}: '+value[:100])
        p.drawOn(c,48,y-h);return y-h-12
    def start(n,title,subtitle):
        c.setFillColor(HexColor(BLUE));c.rect(0,832,596,10,fill=1,stroke=0)
        c.setFont('TimeReport',9);c.drawString(48,799,'pySNSPD  |  ETAPA 4 SIN FOTÓN  |  24-09-2026 UTC')
        c.setFillColor(HexColor('#162c46'));c.setFont('TimeReportBold',20);c.drawString(48,758,title)
        text(subtitle,738,small)
        c.setFillColor(HexColor('#718096'));c.setFont('TimeReport',8)
        c.drawString(48,29,'Resultados numéricos. Producción, circuito y v1.0.0 sin cambios.')
        c.drawRightString(547,29,f'{n} / 4')
    def picture(name,y,height):
        c.drawImage(str(FIG/name),40,y-height,width=515,height=height,preserveAspectRatio=True,anchor='c',mask='auto')
        return y-height-14

    start(1,'Un primer avance temporal completo','Dos perturbaciones pequeñas sobre el mismo núcleo; evolución térmica hasta 1 ps.')
    y=text('<b>La corrida terminó en 71,10 s.</b> Evolucionó todos los nodos libres y comparó dos precisiones temporales, con 27 trabajadores y un coordinador. Reutilizó 256 factorizaciones espectrales: 246 aplicaciones por lotes del operador.',696)
    y=picture('response_and_baseline.png',y,210)
    y=text('<b>La amplitud conserva el 73,59 % de su perturbación inicial a 1 ps; la sonda angular, el 1,751 %.</b> En ambas curvas se resta primero la evolución de la base. Cada norma se divide por su propio valor inicial: no se compara una perturbación grande con otra pequeña.',y)
    y=text('La base no era un equilibrio exacto: también se relaja. Al final, su desplazamiento equivale a <b>5,20 veces</b> la sonda inicial de amplitud y <b>17,25 veces</b> la de fase. Esta deriva es pequeña frente al condensado completo, pero podría ocultar la respuesta buscada si sólo se mirara el mapa total.',y)
    y=text('La perturbación máxima de cualquier estado fue <b>0,2217 % del gap de referencia</b>, dentro del dominio lineal declarado del 2 %. Se verificaron los 37 checkpoints y todos los insumos del operador.',y)
    text('Lectura de las líneas: discontinuas = corrida principal; círculos y línea continua = refinada. Casi coinciden; sus diferencias se muestran separadamente en la página siguiente. Este resultado admite el ensayo térmico afín de 1 ps. La etapa 4 sigue abierta y aún no representa un fotón ni la señal del dispositivo.',y,small)
    c.showPage()

    start(2,'Ver la fase sin esconder la cola','Mapas de la componente transversal de la perturbación angular, después de restar la base.')
    y=picture('transverse_phase_maps.png',696,178)
    y=text('El color mide cuánto se separa la perturbación en la dirección de fase local: es perpendicular al condensado inicial en su plano complejo. Los tres mapas usan <b>la misma escala</b>. El gris central indica que la dirección de fase no se define donde el condensado se anula; no se divide por ese cero.',y,small)
    y=picture('refinement_two_scales.png',y,191)
    y=text('<b>El refinamiento cumple el criterio registrado</b> con diferencias máximas de 0,0383 % en desplazamiento y 0,1155 % en torque angular, medidas contra las señales iniciales. La parte relativa del criterio era 0,5 %, más una tolerancia absoluta.',y)
    y=text('La escala cambia cuando la señal desaparece. A 1 ps queda sólo el <b>0,0163 % del torque inicial</b>; la diferencia entre corridas puede alcanzar el 160 % de ese pequeño resto. No se atribuye precisión fina a esa cola ni se extrae de ella un tiempo exacto de relajación.',y)
    text('Ventana mostrada: |x|, |y| ≤ 4 ℓ₀; ℓ₀ es la longitud de referencia de la malla. La sonda angular pasa por aproximadamente la mitad de su norma entre 0,01 y 0,03 ps. Son respuestas térmicas de prueba, no latencias de detección.',y,small)
    c.showPage()

    start(3,'La disipación tiene el signo esperado','Energía y velocidad evaluadas en los estados guardados; se muestra también el error de aproximación.')
    y=picture('energy_and_taylor_limits.png',696,318)
    y=text('<b>La energía cuadrática disminuye y la disipación KWT más normal es positiva.</b> La base y el estado angular casi coinciden en el panel de energía: su cercanía es un resultado de este ensayo, no una curva ausente.',y)
    y=text('El desajuste instantáneo entre energía liberada y disipación llega a <b>0,07692 %</b> de la disipación. Al volver a evaluar la ley de velocidad sobre el estado desplazado, la diferencia respecto de su aproximación lineal llega a <b>0,2316 %</b> de la velocidad total.',y)
    y=text('Estos valores describen la aproximación de Taylor usada para avanzar. <b>No son todavía el balance no lineal exacto</b>: el condensado se mueve, pero la curvatura espectral permanece fijada en la base. El contraste de la página siguiente examina directamente esa limitación.',y)
    text('Las tasas están expresadas por unidad de tiempo reducida t/tD; tD = 0,4415 ps. Las líneas de tasas casi coincidentes se distinguen por símbolos y por el panel de desajuste. Se usaron los mismos parámetros y contactos; no se ajustó ninguna constante física para facilitar el avance.',y,small)
    c.showPage()

    start(4,'El contraste no lineal también pasa','Las ecuaciones espectrales se vuelven a resolver en estados ya guardados; no se repite la trayectoria.')
    y=picture('nonlinear_snapshot_differences.png',696,197)
    for paragraph in context['final_paragraphs']:y=text(paragraph,y)
    y=text('La dificultad de Newton era numérica: al restar energías próximas a −62,296, un descenso real de −3,153 × 10⁻¹⁵ aparecía como +1,421 × 10⁻¹⁴. La evaluación estable recupera el signo de <b>la misma acción</b>, sin cambiar la condición de descenso ni la tolerancia. El fallo original queda documentado.',y,small)
    y=text('Además, la coordenada auxiliar constante del propagador guardado derivaba un 0,00221 %. La nueva aplicación del avance afín evita esa coordenada. Los resultados originales y sus pruebas se conservan; esta mejora no altera los datos del informe.',y,small)
    text('Datos y hashes: docs/implementation/stage4/time_review_20260924. Los gráficos usan snapshots guardados; las líneas entre ellos guían la lectura. La etapa 4 continúa abierta por la evolución no lineal acoplada y el trabajo espectral no térmico. Etapa 5 y producción no se activan.',y,small)
    c.save()
    print(json.dumps(dict(pdf=str(PDF),pages=4,context_status=context['status'])))


if __name__=='__main__':main()
