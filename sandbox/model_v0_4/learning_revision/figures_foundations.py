"""Pedagogical systems and fields; explicit numerical checks for variations."""
from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Circle

ROOT=Path(__file__).resolve().parents[3]
FIG=ROOT/'docs/modelo_v0_4/figuras'
INK='#243644'; BLUE='#236482'; GREEN='#32866e'; ORANGE='#bf7038'
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':11,'text.color':INK,
                     'axes.spines.top':False,'axes.spines.right':False,'savefig.dpi':220})
fig, axs=plt.subplots(2,2,figsize=(10,9.1),layout='constrained')
for ax in axs.flat:
    ax.set(xlim=(0,10),ylim=(0,10));ax.axis('off')

def describe(ax, title, lines):
    ax.text(.25,9.8,title,fontsize=13,weight='bold',va='top',color=BLUE)
    for y,line in zip([5.4,4.35,3.3,2.25,1.2],lines):
        ax.text(.25,y,line,fontsize=10.5,va='top',linespacing=1.4)

ax=axs[0,0]
describe(ax,'a. Masa y resorte',[
    'Dominio: tiempo t; una coordenada q(t).',
    'Campo espacial: no se necesita.',
    'Estado: posición q y velocidad dq/dt.',
    'Modo: oscilación de la única masa.',
    'Excitación clásica: desplazar y soltar.\nRol: almacenar energía mecánica.'])
ax.plot([.5,.5],[6.6,8.5],lw=3,color=INK)
x=np.linspace(1,5.5,70);ax.plot(x,7.5+.3*np.sin(12*np.pi*(x-1)/4.5),color=BLUE)
ax.add_patch(Rectangle((5.5,6.8),1.3,1.3,facecolor=BLUE+'25',edgecolor=BLUE,lw=2))
ax.annotate('',xy=(8.8,6.4),xytext=(5.4,6.4),arrowprops={'arrowstyle':'->','color':INK})
ax.text(8.9,6.4,'q',va='center')

ax=axs[0,1]
describe(ax,'b. Cavidad electromagnética',[
    'Dominio: posición r y tiempo t.',
    'Campos: E(r,t) y B(r,t), vectoriales.',
    'Configuración: perfiles de E y B.',
    'Modo: onda estacionaria permitida.',
    'Excitación cuántica: un fotón del modo.\nRol: aportar energía electromagnética.'])
x=np.linspace(1,9,200)
ax.add_patch(Rectangle((1,6.4),8,2.2,fill=False,edgecolor=INK,lw=1.4))
ax.plot(x,7.5+np.sin(np.pi*(x-1)/8),color=ORANGE,lw=2)
ax.text(5,6.65,'perfil de una componente',ha='center',fontsize=9)

ax=axs[1,0]
describe(ax,'c. Sólido elástico',[
    'Dominio: posición r y tiempo t.',
    'Campo: desplazamiento vectorial u(r,t).',
    'Estado clásico: u y su velocidad.',
    'Modo: vibración colectiva con\nun patrón espacial elegido.',
    'Excitación cuántica: un fonón del modo.\nRol: transportar y almacenar energía.'])
for x in np.linspace(1,9,9):
    y=7.5+.6*np.sin(2*np.pi*x/8)
    ax.plot([x,x],[7.5,y],color=ORANGE,lw=1)
    ax.add_patch(Circle((x,y),.16,color=GREEN))
ax.axhline(7.5,xmin=.09,xmax=.92,linestyle=':',color=INK,lw=.8)

ax=axs[1,1]
describe(ax,'d. Circuito LC concentrado',[
    'Dominio: tiempo t; variables I(t), V(t).',
    'Campo espacial: eliminado al reducir\nel dispositivo a elementos ideales.',
    'Estado: corriente I y tensión V.',
    'Modo: oscilación LC, con ω = 1/√(LC).',
    'Excitación clásica: cargar el capacitor.\nRol: intercambiar energía eléctrica.'])
ax.plot([1,1,4.3],[6.5,8.2,8.2],color=INK,lw=1.5)
ax.plot([5.5,9,9,6.6],[8.2,8.2,6.5,6.5],color=INK,lw=1.5)
ax.plot([1,3.3],[6.5,6.5],color=INK,lw=1.5)
ax.plot([4.5,4.5],[7.7,8.7],color=BLUE,lw=2)
ax.plot([5.2,5.2],[7.7,8.7],color=BLUE,lw=2)
ax.plot([4.3,4.5],[8.2,8.2],color=INK)
ax.plot([5.2,5.5],[8.2,8.2],color=INK)
ax.text(4.8,9,'C',ha='center')
ax.add_patch(Rectangle((3.3,6.2),3.3,.6,facecolor='white',edgecolor='none'))
x=np.linspace(3.3,6.6,120);ax.plot(x,6.5+.25*np.sin(8*np.pi*(x-3.3)/3.3),color=GREEN,lw=1.5)
ax.text(5,5.85,'L',ha='center',fontsize=10)
for ext in ['png','pdf']:
    fig.savefig(FIG/f'E01_sistemas_campos.{ext}',bbox_inches='tight',facecolor='white')
plt.close(fig)

# Numerical integration checks the functional directly, not its expanded polynomial.
x=np.linspace(0,1,40001);eta=np.sin(np.pi*x);eta_prime=np.pi*np.cos(np.pi*x)
def energy(eps):
    return np.trapezoid(((eps*eta_prime)**2+(1+eps*eta)**2)/2,x)
slopes=[(energy(h)-energy(-h))/(2*h) for h in [1e-2,1e-3,1e-4]]
assert np.allclose(slopes,2/np.pi,atol=1e-8,rtol=0)
assert np.allclose((2+3*x+.2*eta)[[0,-1]],[2,5])
q=np.array([1.,2.]);v=np.array([3.,-1.]);h=1e-4
direction=(np.sum((q+h*v)**2)-np.sum((q-h*v)**2))/(4*h)
assert abs(direction-1)<1e-10
out={'functional_slopes_by_quadrature':slopes,'expected':2/np.pi,
     'fixed_endpoint_values':[2.,5.],'finite_dimensional_directional_derivative':direction,
     'scope':'Comprobaciones independientes de ejemplos; actividades sin resolver.'}
(ROOT/'docs/modelo_v0_4/verificaciones/E_fundamentos_r02.json').write_text(json.dumps(out,indent=2),encoding='utf-8')
print(json.dumps(out,indent=2))
