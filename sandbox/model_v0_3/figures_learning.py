"""Original pedagogical figures and small algebra checks for B/E revision 0.3.
No material parameters, stochastic hotspot rates or DFPT runs are inferred.
"""
from pathlib import Path
import json, platform
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Rectangle, Circle

ROOT=Path(__file__).resolve().parents[2]
FIG=ROOT/'docs/modelo_v0_3/figuras'; DATA=ROOT/'docs/modelo_v0_3/verificaciones'
FIG.mkdir(parents=True,exist_ok=True);DATA.mkdir(parents=True,exist_ok=True)
BLUE='#236482';TEAL='#32866e';ORANGE='#bf7038';PURPLE='#745f98';INK='#22333f'
plt.rcParams.update({'font.size':10,'axes.titlesize':11,'axes.labelsize':10,'font.family':'DejaVu Sans',
 'axes.spines.top':False,'axes.spines.right':False,'savefig.dpi':200,'figure.dpi':130,'legend.frameon':False})

def save(fig,name):
    for ext in ('png','pdf'):fig.savefig(FIG/f'{name}.{ext}',bbox_inches='tight',facecolor='white')
    plt.close(fig)

def box(ax,x,y,w,h,label,color=BLUE,size=10):
    ax.add_patch(FancyBboxPatch((x,y),w,h,boxstyle='round,pad=0.02,rounding_size=0.035',
                    linewidth=1.3,edgecolor=color,facecolor=color+'14'))
    ax.text(x+w/2,y+h/2,label,ha='center',va='center',color=INK,fontsize=size)

def arrow(ax,xy0,xy1,color=INK):
    ax.annotate('',xy=xy1,xytext=xy0,arrowprops={'arrowstyle':'->','lw':1.5,'color':color})

def main():
    fig,axs=plt.subplots(1,3,figsize=(10.8,3.3),layout='constrained')
    x=np.linspace(0,2*np.pi,300)
    axs[0].plot(x,np.sin(x),color=BLUE,lw=2.5);axs[0].axhline(0,color='#aab4bc',lw=.7)
    axs[0].set(xlabel='Posición',ylabel='Amplitud del modo',title='Un modo espacial (analogía)')
    for ax,nmax,title,col in [(axs[1],4,'Ocupación bosónica',TEAL),(axs[2],2,'Ocupación fermiónica',ORANGE)]:
        for n in range(nmax):
            ax.plot([.15,.82],[n,n],lw=2,color=col)
            ax.text(.85,n,f'n = {n}',va='center')
            for j in range(n):ax.plot(.24+.16*j,n,'o',color=col,markersize=7)
        ax.set(xlim=(0,1.3),ylim=(-.4,4.35),title=title)
        ax.axis('off');ax.text(.6,3.85,'Un estado completamente\nespecificado',ha='center',va='center',fontsize=9)
    save(fig,'E_01_campos_ocupacion')

    fig,ax=plt.subplots(figsize=(8.8,3.4),layout='constrained');ax.axis('off');ax.set(xlim=(0,10),ylim=(0,4))
    ax.text(5.6,3.72,'Índice de espín de la base',ha='center',weight='bold')
    ax.text(4,3.25,r'$\uparrow$',ha='center',fontsize=16);ax.text(7.25,3.25,r'$\downarrow$',ha='center',fontsize=16)
    ax.text(.1,2.4,'Nambu: destrucción',fontsize=10,va='center')
    ax.text(.1,1.15,'Nambu: conjugado\ncon inversión temporal',fontsize=10,va='center')
    labels=[r'$c_{\mathbf{k}\uparrow}$',r'$c_{\mathbf{k}\downarrow}$',r'$c^\dagger_{-\mathbf{k}\downarrow}$',r'$-c^\dagger_{-\mathbf{k}\uparrow}$']
    for i,l in enumerate(labels):box(ax,2.5+(i%2)*3.25,1.9-(i//2)*1.25,3,1,l,BLUE if i<2 else TEAL,17)
    ax.text(5.6,.12,r'$\tau$: entre filas      $\sigma$: entre columnas',ha='center',fontsize=11)
    save(fig,'E_02_nambu_spin')

    fig,axs=plt.subplots(1,2,figsize=(10.2,3.65),layout='constrained')
    x=np.linspace(0,np.pi,500);base=np.sin(x);pert=np.sin(x)*np.exp(-((x-1.5)/.25)**2)
    axs[0].plot(x,base,lw=2,label=r'$y(x)$',color=BLUE)
    axs[0].plot(x,base+.3*pert,lw=2,label=r'$y(x)+\epsilon\eta(x)$',color=ORANGE)
    axs[0].fill_between(x,base,base+.3*pert,color=ORANGE,alpha=.2)
    axs[0].set(xlabel=r'$x$',ylabel='Perfil',title='Cambiar localmente una función');axs[0].legend(fontsize=9)
    yy,zz=np.meshgrid(np.linspace(-1.25,1.25,150),np.linspace(-3,3,170));f=.5*yy**2+.5*(zz-2*yy)**2
    axs[1].contourf(yy,zz,f,levels=[0,.08,.2,.4,.8,1.5,3,6,12],cmap='Blues')
    y=np.linspace(-1.2,1.2,150);axs[1].plot(y,2*y,color=ORANGE,lw=2.5,label=r'$z^*(y)=2y$')
    axs[1].set(xlabel=r'$y$',ylabel=r'$z$',title='Eliminar una variable estacionaria');axs[1].legend(loc='upper left',fontsize=9)
    save(fig,'E_03_variacion_envolvente')

    fig,axs=plt.subplots(1,2,figsize=(10.6,3.6),layout='constrained')
    ax=axs[0];ax.set(xlim=(-.6,7.6),ylim=(-.5,2.5));ax.axis('off');ax.set_title('El electrón y la vacante se desplazan al revés')
    for row,hole in [(1.8,3),(.35,4)]:
        for i in range(8):
            ax.add_patch(Circle((i,row),.18,facecolor='white' if i==hole else BLUE,edgecolor=ORANGE if i==hole else BLUE,lw=2))
    arrow(ax,(4,1.37),(3,1.37),BLUE);ax.text(3.5,1.08,'salto del electrón',ha='center',fontsize=9,color=BLUE)
    arrow(ax,(3,.8),(4,.8),ORANGE);ax.text(5,.82,'hueco',fontsize=9,color=ORANGE)
    ax.text(.1,-.25,'Puntos llenos: estados ocupados; círculo: vacante',fontsize=8.8)
    ax=axs[1];ax.set(xlim=(0,4),ylim=(-3,4.1));ax.set_title('Energía respecto de una referencia llena')
    for y,col in [(3,ORANGE),(0,'#89939b'),(-2,BLUE)]:ax.hlines(y,.4,3.5,color=col,lw=1.6)
    ax.text(3.6,0,r'$\mu$',va='center');ax.text(.5,3.3,'Añadir electrón: +3 meV, carga -e',fontsize=9)
    ax.text(.5,-2.6,'Retirar electrón: +2 meV, carga +e',fontsize=9)
    ax.plot(2,3,'o',color=ORANGE,ms=9);ax.plot(2,-2,'x',color=BLUE,ms=10,mew=2)
    ax.set_ylabel(r'$\varepsilon-\mu$ (meV)');ax.set_xticks([])
    save(fig,'E_04_electron_hueco')

    fig,axs=plt.subplots(1,3,figsize=(11,3.5),layout='constrained');ax=axs[0];ax.axis('off');ax.set(xlim=(-.3,3.3),ylim=(-1,1.5));ax.set_title('Dos coordenadas, dos modos')
    for center in [1,2]:ax.add_patch(Circle((center,1),.11,color=BLUE))
    for lo,hi in [(0,.85),(1.15,1.85),(2.15,3)]:
        xx=np.linspace(lo,hi,15);yy=1+.04*np.where(np.arange(15)%2,1,-1);yy[[0,-1]]=1;ax.plot(xx,yy,color=INK)
    ax.vlines([0,3],.8,1.2,color=INK,lw=3)
    ax.plot([0,1,2,3],[0,.25,.25,0],'-o',color=TEAL);ax.text(.2,.4,r'Modo 1: $\omega/\sqrt{K/m}=1$',fontsize=9)
    ax.plot([0,1,2,3],[-.65,-.4,-.9,-.65],'-o',color=ORANGE);ax.text(.2,-.23,r'Modo 2: $\omega/\sqrt{K/m}=\sqrt{3}$',fontsize=9)
    omegas=np.array([1,np.sqrt(3)])
    for ax,weights,title,col in [(axs[1],[1,1],'DOS: área por modo',BLUE),(axs[2],[.1,.2],'Interacción: pesos distintos',ORANGE)]:
        ax.vlines(omegas,0,weights,color=col,lw=3);ax.plot(omegas,weights,'o',color=col)
        ax.set(xlim=(.4,2.1),ylim=(0,max(weights)*1.3),xlabel=r'$\Omega/(\hbar\sqrt{K/m})$',title=title)
    axs[1].set_ylabel('Área de la delta');axs[2].set_ylabel(r'$W_i$ (unidad de energía)')
    save(fig,'E_05_modos_dos')

    fig,ax=plt.subplots(figsize=(10.8,4.2),layout='constrained');ax.axis('off');ax.set(xlim=(0,12),ylim=(0,4.4))
    box(ax,.15,2.65,2.4,1.2,'Estructura y\nequilibrio DFT',BLUE)
    box(ax,3.2,2.65,2.6,1.2,'Respuesta electrónica\na desplazamientos\nDFPT',TEAL)
    box(ax,6.55,2.8,2.3,.9,'Curvaturas\nmodos y frecuencias',BLUE)
    box(ax,6.55,1.4,2.3,.9,'Perturbación potencial\nelementos de interacción',ORANGE,9)
    box(ax,9.5,2.8,2.2,.9,r'DOS $F(\Omega)$',BLUE)
    box(ax,9.5,1.4,2.2,.9,r'Pesos $\alpha^2F(\Omega)$',ORANGE)
    box(ax,3.2,.05,5.65,.95,'B: espectro superconductor + ocupaciones\n+ entradas materiales → tasas y balances',PURPLE,10)
    arrow(ax,(2.57,3.25),(3.18,3.25));arrow(ax,(5.82,3.25),(6.53,3.25));arrow(ax,(5.4,2.63),(6.6,2.32))
    arrow(ax,(8.87,3.25),(9.48,3.25));arrow(ax,(8.87,1.85),(9.48,1.85));arrow(ax,(10.55,1.35),(8.88,.55))
    ax.text(.25,.65,'No se calcula todavía\nla absorción ni la gTDGL.',fontsize=9,color=INK)
    save(fig,'E_06_dfpt_a_cinetica')

    fig,axs=plt.subplots(1,2,figsize=(10.7,3.75),layout='constrained');ax=axs[0];ax.axis('off');ax.set(xlim=(0,6),ylim=(0,4))
    box(ax,.15,2.65,1.65,.95,'Fotón\nabsorbido',BLUE)
    box(ax,2.3,2.65,2.5,.95,'Cascada y reparto\nfluctuante',ORANGE)
    box(ax,.6,.55,2.15,1.0,'Energía retenida\nelectrones + fonones',TEAL,9.5)
    box(ax,3.3,.55,2.15,1.,'Energía que escapa\nal sustrato',PURPLE,9.5)
    arrow(ax,(1.85,3.12),(2.25,3.12));arrow(ax,(3.25,2.6),(1.65,1.6));arrow(ax,(3.85,2.6),(4.35,1.6))
    ax.set_title('Absorción y retención son etapas distintas')
    ax=axs[1];x=np.linspace(-4,4,700)
    smooth=np.exp(-x*x/2)/np.sqrt(2*np.pi)
    lobes=(np.exp(-(x-1.3)**2/(2*.45**2))+np.exp(-(x+1.3)**2/(2*.45**2)))/(2*np.sqrt(2*np.pi)*.45)
    ax.plot(x,smooth,color=BLUE,lw=2,label='Cierre espacial liso')
    ax.plot(x,lobes,color=ORANGE,lw=2,label='Evento hipotético de dos lóbulos')
    ax.set(xlabel='Posición ilustrativa',ylabel='Densidad normalizada',title='Igual energía, distinta forma espacial')
    ax.set_ylim(-.13,.61)
    ax.legend(fontsize=8.5,loc='upper left');ax.text(.5,.04,'No fija probabilidades ni distancias del detector',transform=ax.transAxes,ha='center',fontsize=8.5)
    save(fig,'B_05_burbuja_y_retencion')

    stiffness=np.array([[2.,-1.],[-1.,2.]])
    eig=np.linalg.eigvalsh(stiffness)
    H=np.array([[3.,4.],[4.,-3.]])
    ev,vec=np.linalg.eigh(H);u2=vec[:,1]**2
    assert np.allclose(eig,[1,3]);assert np.allclose(ev,[-5,5]);assert np.allclose(u2,[.8,.2])
    checks={'scope':'Sólo ejemplos resueltos; no contiene soluciones de las actividades evaluadas',
      'host':platform.node(),'python':platform.python_version(),'E2_BdG_eigenvalues':ev.tolist(),'E2_positive_weights':u2.tolist(),
      'E3_reduced_curvature':5-(-2)**2,'E5_stiffness_eigenvalues':eig.tolist(),'E5_lambda_example':float(2*(.1+.2/np.sqrt(3))),
      'B5_scope':'Esquema construido, no simulación estocástica ni tasa de eventos'}
    (DATA/'E_ejemplos_verificados.json').write_text(json.dumps(checks,indent=2,ensure_ascii=False),encoding='utf-8')
    print(json.dumps(checks,indent=2,ensure_ascii=False))

if __name__=='__main__':main()
