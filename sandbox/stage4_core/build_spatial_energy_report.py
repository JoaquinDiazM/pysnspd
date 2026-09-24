"""Results report for the spatial thermal-energy reference, with static plots."""
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
DATA=ROOT/'docs/implementation/stage4/spatial_energy_20260924'
RAW=DATA/'raw/stage4_spatial_energy_reference_20260924'
FIG=DATA/'figures'
PDF=ROOT/'output/pdf/implementation/Informe_etapa_4_energia_espacial.pdf'
COLORS=['#2166ac','#dd792e','#22846c']
plt.rcParams.update({'font.size':9,'axes.spines.top':False,'axes.spines.right':False,'savefig.dpi':190})


def figures():
    FIG.mkdir(parents=True,exist_ok=True)
    fig,ax=plt.subplots(1,2,figsize=(10,3.1),layout='constrained')
    old=ROOT/'docs/implementation/stage4/followup_20260923/raw/stage4_radial_reference_20260923/R8_N256.npz'
    with np.load(old) as a:
        ids=a['radius']<=4
        ax[0].plot(a['radius'][ids],a['candidate_force_delta_0.1'][ids],color='#b66a45',label='Cierre anterior, δ/Δ₀ = 0,1')
        ax[0].plot(a['radius'][ids],a['force_raw'][ids],color='#1d2d41',lw=2,label='Usadel radial independiente')
    for n,color in zip([33,65,129],COLORS):
        with np.load(RAW/f'radial_{n}_N256.npz') as a:
            xy=a['coordinates_bar'];ids=np.flatnonzero((abs(xy[:,1])<1e-10)&(xy[:,0]>0)&(xy[:,0]<=4))
            ref=np.real(a['reference_radial_force'][ids]);val=a['amplitude_force'][ids]
            if n==129:ax[0].plot(xy[ids,0],val,'o',ms=2.8,mfc='none',color=COLORS[2],label='Nueva energía 2D: 129² nodos')
            ax[1].plot(xy[ids,0],100*(val-ref)/max(abs(ref)),'o-',ms=2.5,color=color,label=f'{n}² nodos')
    ax[0].set(xlabel='Radio / ℓ₀',ylabel='Fuerza de amplitud normalizada',title='Se elimina el pico de signo opuesto',ylim=(-13,2),xlim=(0,4))
    ax[1].axhline(0,color='#1d2d41',lw=1,ls='--',label='Usadel radial: diferencia cero')
    ax[1].set(xlabel='Radio / ℓ₀',ylabel='Diferencia / máximo de referencia (%)',title='El error restante disminuye con la malla')
    for axis in ax:axis.grid(alpha=.2);axis.legend(fontsize=7)
    fig.savefig(FIG/'force_recovery.png');plt.close(fig)
    with np.load(RAW/'asymmetric_65_N256.npz') as a:
        xy=a['coordinates_bar']*4.698490897269843
        d=a['d'];f=a['f_lowest'];phase=np.angle(f*np.conj(d))*180/np.pi
        phase[np.abs(d)<1e-8]=np.nan
        fig,axes=plt.subplots(1,3,figsize=(10,3.1),layout='constrained')
        extent=[xy[:,0].min(),xy[:,0].max(),xy[:,1].min(),xy[:,1].max()]
        for axis,field,title,cmap,label in zip(axes,[abs(d),abs(f),phase],
            ['Condensado impuesto','Correlación de pares','Diferencia de fase espectral'],
            ['viridis','viridis','RdBu_r'],['|Δ| / kBTc','|f₀|','Grados']):
            image=axis.imshow(field.reshape(65,65).T,origin='lower',extent=extent,cmap=cmap,
                              **({'vmin':-18,'vmax':18} if cmap=='RdBu_r' else {}))
            axis.set(title=title,xlabel='x (nm)',ylabel='y (nm)')
            fig.colorbar(image,ax=axis,label=label,shrink=.83)
        fig.savefig(FIG/'spectral_environment.png');plt.close(fig)


def main():
    figures();PDF.parent.mkdir(parents=True,exist_ok=True)
    pdfmetrics.registerFont(TTFont('Report',findfont('DejaVu Sans')))
    pdfmetrics.registerFont(TTFont('ReportBold',findfont('DejaVu Sans:weight=bold')))
    pdfmetrics.registerFontFamily('Report',normal='Report',bold='ReportBold',italic='Report',boldItalic='ReportBold')
    c=canvas.Canvas(str(PDF),pagesize=(595.276,841.89));c.setTitle('Etapa 4: recuperación de la energía espacial')
    style=ParagraphStyle('body',fontName='Report',fontSize=10.2,leading=15,textColor=HexColor('#253343'))
    small=ParagraphStyle('small',parent=style,fontSize=8.5,leading=12)
    def text(value,y,st=style):
        p=Paragraph(value,st);_,h=p.wrap(499,730);p.drawOn(c,48,y-h);return y-h-13
    def start(n,title,subtitle):
        c.setFillColor(HexColor(COLORS[0]));c.rect(0,832,596,10,fill=1,stroke=0)
        c.setFont('Report',9);c.drawString(48,799,'pySNSPD  |  ETAPA 4 SIN FOTÓN  |  ACTUALIZACIÓN 24-09-2026 UTC')
        c.setFillColor(HexColor('#162c46'));c.setFont('ReportBold',20);c.drawString(48,758,title)
        text(subtitle,738,small)
        c.setFillColor(HexColor('#718096'));c.setFont('Report',8)
        c.drawString(48,29,'Referencia térmica experimental. Producción y v1.0.0 sin cambios.');c.drawRightString(547,29,f'{n} / 3')
    def picture(name,y,h):
        c.drawImage(str(FIG/name),42,y-h,width=511,height=h,preserveAspectRatio=True,anchor='c',mask='auto');return y-h-16
    start(1,'Una respuesta espacial que sí coincide','T = 0,9 K. Perfil de núcleo prescrito; ℓ₀ = 4,70 nm. Cálculo térmico, sin evolución temporal.')
    y=text('<b>El contraste térmico que fallaba queda recuperado.</b> El cierre anterior difería de Usadel en 158 % y daba una fuerza de signo opuesto cerca del núcleo. La nueva referencia 2D llega a <b>0,18 %</b> de diferencia sin ajustar δ ni añadir otra fuerza de gradiente.',696)
    y=picture('force_recovery.png',y,202)
    y=text('<b>Al refinar:</b> 33 × 33 nodos → 2,90 %; 65 × 65 → 0,72 %; 129 × 129 → 0,18 %. El error cae aproximadamente cuatro veces al dividir por dos el paso espacial. Se compara el mismo perfil y la misma suma de 256 frecuencias con la referencia radial independiente.',y)
    y=text('La solución conserva el espectro como campos espaciales adicionales y lo ajusta antes de evaluar energía, fuerza y corriente. En el régimen térmico declarado es una representación del funcional de Usadel, no un ajuste a esta curva. La discretización introduce el error de malla que muestra el gráfico derecho. Por separado, pasar de 128 a 256 frecuencias cambia la fuerza 1,87 %: el 0,18 % no es una cota del error total.',y)
    y=text('<b>Coste medido:</b> 2.048 problemas espectrales en 60,35 s, con 27 trabajadores y un coordinador. Se usaron como máximo 28 hilos de Geminga, reservando dos núcleos completos. Los ocho casos compartieron el grupo de procesos.',y)
    text('La figura izquierda superpone deliberadamente las soluciones coincidentes; el panel derecho hace visible su diferencia. Las fronteras espectrales de este cuadrado auxiliar proceden de la referencia radial. No es todavía un ensayo de la cinta completa de 80 nm.',y,small);c.showPage()
    start(2,'El espectro responde al entorno','Un segundo perfil rompe la simetría radial y permite comprobar la respuesta compleja en dos dimensiones.')
    y=picture('spectral_environment.png',696,202)
    y=text('Las correlaciones de pares describen cómo responde el estado electrónico al condensado y a sus vecinos. Su fase puede diferir de la fase local de Δ: en el caso asimétrico la diferencia del modo inferior llega a <b>18,0°</b>. Obligar a ambas fases a coincidir eliminaría una parte de esa respuesta espacial.',y)
    y=text('La representación usa coordenadas que permanecen regulares cuando Δ = 0. El centro forma parte de la malla; no hay que dividir allí por la amplitud del condensado ni insertar un radio artificial. En el centro asimétrico, Δ = 0 pero la correlación del modo inferior vale |f₀| = 0,0705: la influencia de los vecinos sigue presente y la fuerza es finita.',y)
    y=text('<b>Una sola energía comprobada:</b> al perturbar la amplitud, la variación independiente de energía y el trabajo de la fuerza coinciden con error relativo 9,34 × 10<super>-8</super>; al perturbar la fase, 2,26 × 10<super>-7</super>. La corriente también se deriva de esa energía, y las pruebas separan el trabajo de los bordes de la identidad interior.',y)
    y=text('<b>Alcance:</b> los espectros son estacionarios para un condensado todavía impuesto. Que su fuerza esté bien calculada no significa que ese condensado ya haya alcanzado su propio equilibrio. Tampoco demuestra la estabilidad de un transiente fuera del equilibrio.',y)
    text('Los 2.048 archivos espectrales completos permanecen en Geminga con hashes verificados. La entrega conserva identidades, comparaciones, mapas y pruebas. Las sumas finitas y los cortes 128/256 se mantienen explícitos: no se añade una cola sólo a la fuerza.',y,small);c.showPage()
    start(3,'Siguiente dato: el núcleo autoconsistente','El ensayo siguiente determina la forma del núcleo a partir de la misma energía; no asigna tiempo físico a las iteraciones.')
    y=text('<b>Dos pasos alternados:</b> primero se ajustan los campos espectrales para el condensado actual; después se minimiza exactamente el bloque de energía que depende del condensado, manteniendo esos espectros. Ambos pasos conservan el funcional térmico y sus condiciones de borde. La actualización del condensado tiene un descenso energético demostrado algebraicamente.',696)
    y=text('Se compararán dos mallas y dos cortes de frecuencias; cada corte se relajará por separado. También se incluye el núcleo inicialmente asimétrico. Se registrarán perfil, tamaño del núcleo, corriente, energía y residuo de autoconsistencia. El objetivo inicial del residuo es 0,1 % en norma sobre el núcleo, no una exigencia de precisión arbitraria sobre cada punto.',y)
    y=text('<b>La campaña comparte los recursos:</b> las frecuencias y los casos activos utilizan el mismo grupo. Cada caso puede terminar y liberar capacidad. Se mostrarán progreso y una estimación del trabajo máximo restante; esa estimación no promete cuándo se alcanzará la autoconsistencia. Los checkpoints guardarán condensado y espectros del mismo estado.',y)
    pilot=DATA/'self_consistent_pilot_receipt.json'
    if pilot.exists():
        receipt=json.loads(pilot.read_text(encoding='utf8'))
        y=text(receipt['report_sentence'],y)
    y=text('<b>La etapa 4 permanece abierta.</b> Este cálculo puede admitir una referencia térmica autoconsistente del núcleo. Antes de usarla en el detector todavía debe justificarse su acoplamiento con las distribuciones fuera del equilibrio. Reemplazarlas por una temperatura equivalente cambiaría el modelo y no se hace aquí.',y)
    y=text('Se revisaron los métodos de la memoria: se conserva la posibilidad de reutilizar la malla dual, los enlaces de calibre, las matrices dispersas y el circuito completo. Este bloque nuevo es experimental; no reemplaza el solver temporal de producción. La ventana futura hasta Vout más margen sigue siendo únicamente un cambio del horizonte.',y)
    text('Fundamento: especialización térmica del funcional de Virtanen, Vargunin y Silaev, arXiv:1909.00992, y revisión de los anexos de la memoria. Detalles y límites: memory_and_sources.md, variational_contract.md y self_consistent_plan.json. El comando vigente se copia en el chat y en GEMINGA_COMMANDS.md.',y,small)
    c.save();print(json.dumps(dict(pdf=str(PDF),pages=3)))


if __name__=='__main__':main()
