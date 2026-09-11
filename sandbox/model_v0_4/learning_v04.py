"""New v0.4 learning figures and independent checks of worked examples only."""
from pathlib import Path
import json,platform
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch,Circle
ROOT=Path(__file__).resolve().parents[2];FIG=ROOT/'docs/modelo_v0_4/figuras';DATA=ROOT/'docs/modelo_v0_4/verificaciones'
FIG.mkdir(parents=True,exist_ok=True);DATA.mkdir(parents=True,exist_ok=True)
BLUE='#236482';GREEN='#32866e';ORANGE='#bf7038';INK='#22333f'
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False,'savefig.dpi':200})
def box(ax,x,y,w,h,label,c=BLUE,size=11):
    ax.add_patch(FancyBboxPatch((x,y),w,h,boxstyle='round,pad=.015',facecolor=c+'12',edgecolor=c,lw=1.4))
    ax.text(x+w/2,y+h/2,label,ha='center',va='center',fontsize=size,color=INK)
def arr(ax,start,end):ax.annotate('',xy=end,xytext=start,arrowprops={'arrowstyle':'->','lw':1.5,'color':INK})
def save(fig,name):
    for ext in ['png','pdf']:fig.savefig(FIG/f'{name}.{ext}',bbox_inches='tight',facecolor='white')
    plt.close(fig)

fig,ax=plt.subplots(figsize=(10,3.2));ax.axis('off');ax.set(xlim=(0,12),ylim=(0,4))
for x,title,detail,c in [(0.1,'Campos elementales','Electrón, fotón\nEspecie y estado cuántico',BLUE),(4.1,'Excitaciones del sólido','Fonón, cuasipartícula\nModo y ocupación',GREEN),(8.1,'Variables del dispositivo','Corriente y tensión\nEstado del circuito',ORANGE)]:
    box(ax,x,.8,3.7,2.3,title+'\n\n'+detail,c,10.5)
ax.text(6,.2,'Cambiar de escala cambia las variables que bastan para responder una pregunta.',ha='center',fontsize=10)
save(fig,'E01_escalas')

fig,axs=plt.subplots(1,2,figsize=(9.5,3.2),layout='constrained')
ax=axs[0];ax.axis('off');ax.set(xlim=(-.5,5.5),ylim=(-.4,2.6));ax.set_title('Un electrón llena la vacante vecina')
for y,hole in [(2,2),(.3,3)]:
    for i in range(6):ax.add_patch(Circle((i,y),.15,facecolor='white' if i==hole else BLUE,edgecolor=ORANGE if i==hole else BLUE,lw=2))
arr(ax,(3,1.55),(2,1.55));ax.text(2.5,1.15,'electrón',ha='center',fontsize=9)
arr(ax,(2,.85),(3,.85));ax.text(4,.85,'vacante',va='center',fontsize=9)
ax=axs[1];ax.set(xlim=(0,3),ylim=(-2.5,3.7),ylabel=r'$\varepsilon-\mu$ (meV)',xticks=[],title='La referencia decide el signo')
ax.axhline(0,color='grey',lw=1);ax.hlines([3,-2],.3,2.7,colors=[ORANGE,BLUE],lw=2)
ax.text(.4,3.25,'Añadir: cuesta 3 meV',fontsize=9);ax.text(.4,-1.6,'Retirar: cuesta 2 meV',fontsize=9)
ax.text(2.5,.15,r'$\mu$',fontsize=11);save(fig,'E04_referencia_hueco')

fig,ax=plt.subplots(figsize=(9.2,3.2));ax.axis('off');ax.set(xlim=(0,10),ylim=(0,4))
box(ax,.4,1.15,2.3,1.6,'Vacío\n'+r'$|0\rangle=(1,0)^T$',BLUE)
box(ax,7.2,1.15,2.3,1.6,'Ocupado\n'+r'$|1\rangle=(0,1)^T$',GREEN)
arr(ax,(2.9,2.45),(7,2.45));ax.text(5,2.85,r'$c^\dagger$: añadir',ha='center')
arr(ax,(7,1.4),(2.9,1.4));ax.text(5,.95,r'$c$: retirar',ha='center')
ax.text(5,.22,r'$c|0\rangle=0$ y $c^\dagger|1\rangle=0$: vector cero, no otro estado físico.',ha='center',fontsize=10)
save(fig,'E06_ocupacion_operadores')

fig,ax=plt.subplots(figsize=(9.8,3.1));ax.axis('off');ax.set(xlim=(0,12),ylim=(0,4))
box(ax,.2,1,3,2,r'$\Psi=\binom{c_1}{c_2^\dagger}$'+'\nColumna de operadores',BLUE)
box(ax,4.35,1,3.3,2,'H = [[ξ, Δ], [Δ*, −ξ]]\n\nMezcla por emparejamiento',ORANGE,10)
box(ax,8.75,1,3,2,r'$\binom{\xi c_1+\Delta c_2^\dagger}{\Delta^*c_1-\xi c_2^\dagger}$',GREEN,12)
arr(ax,(3.25,2),(4.3,2));arr(ax,(7.7,2),(8.7,2));ax.text(6,.3,'Las dos entradas son operaciones relacionadas; no son dos electrones adicionales.',ha='center',fontsize=10)
save(fig,'E07_columna_nambu')

fig,axs=plt.subplots(1,2,figsize=(9.3,3.2),layout='constrained')
x=np.linspace(0,1,301);y=np.ones_like(x);eta=np.sin(np.pi*x)
axs[0].plot(x,y,color=BLUE,label=r'$y(x)=1$');axs[0].plot(x,y+.2*eta,color=ORANGE,label=r'$y+0.2\eta$')
axs[0].set(xlabel='Posición x',ylabel='Perfil y',title='Elegir una dirección de cambio');axs[0].legend(fontsize=9)
eps=np.linspace(-.3,.3,201);energy=.5+2*eps/np.pi+(np.pi**2+1)*eps**2/4
axs[1].plot(eps,energy,color=BLUE,label=r'$J(\epsilon)=\mathcal{F}[y+\epsilon\eta]$')
axs[1].plot(eps,.5+2*eps/np.pi,'--',color=ORANGE,label='Valor inicial + término lineal')
axs[1].set(xlabel=r'Tamaño $\epsilon$',ylabel='Energía del ejemplo',title='Derivar una curva ordinaria');axs[1].legend(fontsize=8.5)
save(fig,'E03_curva_variacion')

fig,axs=plt.subplots(1,2,figsize=(9.4,3.4),layout='constrained');ax=axs[0];ax.axis('off');ax.set(xlim=(-.2,3.2),ylim=(-1,1.8),title='Dos masas: dos formas independientes')
ax.plot([0,1,2,3],[1.3]*4,color=INK,lw=1);ax.plot([1,2],[1.3]*2,'o',color=BLUE,ms=12);ax.vlines([0,3],1.1,1.5,lw=3,color=INK)
ax.plot([0,1,2,3],[.2,.5,.5,.2],'-o',color=GREEN);ax.text(.1,.75,r'En fase: $\omega_1=\sqrt{K/m}$',fontsize=10)
ax.plot([0,1,2,3],[-.6,-.3,-.9,-.6],'-o',color=ORANGE);ax.text(.1,-.1,r'En oposición: $\omega_2=\sqrt{3K/m}$',fontsize=10)
ax=axs[1];ax.vlines([1,np.sqrt(3)],0,1,colors=[GREEN,ORANGE],lw=3);ax.plot([1,np.sqrt(3)],[1,1],'o',color=BLUE)
ax.set(xlim=(.5,2),ylim=(0,1.3),xlabel=r'$\omega/\sqrt{K/m}$',ylabel='Área de cada línea: 1 modo',title='DOS: conservar el conteo')
save(fig,'E05_modos_conteo')

fig,axs=plt.subplots(1,2,figsize=(9.5,3.4),layout='constrained');u=np.linspace(-.4,.4,200)
axs[0].plot(u,u*u,color=BLUE,label='Energía de equilibrio')
axs[0].plot([-.2,0,.2],[.04,0,.04],'o',color=ORANGE)
axs[0].set(xlabel='Desplazamiento atómico u',ylabel='Energía relativa',title='Curvatura: cuánto cuesta desplazar');axs[0].legend(fontsize=9)
ax=axs[1];ax.vlines([1,2],0,[.1,.3],colors=[GREEN,ORANGE],lw=3);ax.set(xlim=(.5,2.5),ylim=(0,.4),xlabel=r'Energía fonónica $\Omega$ (unidad del ejemplo)',ylabel=r'Peso $W$ (misma unidad de energía)',title='Interacción: cuánto pesa cada modo')
ax.text(1,.12,r'$W_1=0.1$',ha='center',fontsize=9);ax.text(2,.32,r'$W_2=0.3$',ha='center',fontsize=9)
save(fig,'E08_curvatura_acoplamiento')

# Matrix checks do not duplicate the evaluated exercises.
c=np.array([[0.,1.],[0.,0.]]);cd=c.T;I=np.eye(2)
assert np.allclose(c@cd+cd@c,I)
assert np.allclose(cd@c,np.diag([0,1]))
H=np.array([[2.,1.],[1.,-2.]])
ev=np.linalg.eigvalsh(H);assert np.allclose(ev,[-np.sqrt(5),np.sqrt(5)])
# Independent central differentiation of the scalar curve behind the first variation.
analytic=2/np.pi;errors=[]
for h in [1e-2,1e-3,1e-4]:
    J=lambda z: .5+2*z/np.pi+(np.pi**2+1)*z*z/4
    errors.append(abs((J(h)-J(-h))/(2*h)-analytic))
assert max(errors)<1e-10
out={'revision':'0.4','host':platform.node(),'python':platform.python_version(),
 'scope':'Sólo ejemplos resueltos del cuaderno; no respuestas de actividades evaluadas',
 'anticommutator_residual':float(np.max(abs(c@cd+cd@c-I))),
 'nambu_example_eigenvalues':ev.tolist(),'gateaux_central_errors':errors,
 'two_mass_stiffness_eigenvalues':np.linalg.eigvalsh([[2,-1],[-1,2]]).tolist(),
 'lambda_two_mode_example':2*(.1/1+.3/2)}
(DATA/'E_ejemplos_v04.json').write_text(json.dumps(out,indent=2,ensure_ascii=False),encoding='utf-8')
print(json.dumps(out,indent=2))
