"""Draw the three-state thesis circuit with a passive device port."""
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle,Circle
ROOT=Path(__file__).resolve().parents[3]
OUT=ROOT/'docs/implementation/stage2/practical_review_20260922/figures/circuit_topology.png'
fig,ax=plt.subplots(figsize=(10,4.4),layout='constrained')
ax.set(xlim=(-.1,8.1),ylim=(-.6,4.2));ax.axis('off')
def wire(xs,ys):ax.plot(xs,ys,color='#263c4a',lw=1.8)
def resistor(x,y,label,vertical=False):
    ax.add_patch(Rectangle((x-.17,y-.45) if vertical else (x-.45,y-.17),.34 if vertical else .9,.9 if vertical else .34,fill=False,lw=1.6))
    ax.text(x+.35,y,label,va='center',fontsize=12) if vertical else ax.text(x,y+.3,label,ha='center',fontsize=12)
def coil(x,y,label,vertical=False):
    u=np.linspace(-.45,.45,150);v=.12*np.sin(np.linspace(0,6*np.pi,150))
    wire(x+v,y+u) if vertical else wire(x+u,y+v)
    ax.text(x+.35,y,label,va='center',fontsize=12) if vertical else ax.text(x,y+.35,label,ha='center',fontsize=12)
def arrow(a,b,label,offset=(0,.16)):
    ax.annotate('',xy=b,xytext=a,arrowprops=dict(arrowstyle='->',color='#007e80',lw=1.8))
    ax.text((a[0]+b[0])/2+offset[0],(a[1]+b[1])/2+offset[1],label,color='#007e80',fontsize=12,ha='center')
wire([.5,.5],[0,1.08]);ax.add_patch(Circle((.5,1.5),.42,fill=False,lw=1.8));wire([.5,.5,1.05],[1.92,3,3])
ax.text(.5,1.64,'+',ha='center',va='center');ax.text(.5,1.35,'−',ha='center',va='center')
ax.text(1.08,1.5,r'$V_{bias}$',ha='left',va='center',fontsize=12)
resistor(1.5,3,r'$R_b$');wire([1.95,2.35],[3,3]);coil(2.8,3,r'$L_b$');wire([3.25,5.35],[3,3])
ax.plot([4.2],[3],'o',color='#263c4a');ax.text(4.2,3.32,r'$V_d$',ha='center',fontsize=12)
wire([5.35,5.35],[2.64,3.36]);wire([5.53,5.53],[2.64,3.36]);wire([5.53,7],[3,3]);wire([7,7],[3,1.95])
ax.text(5.44,3.57,r'$C_c$',ha='center',fontsize=12);ax.text(5.44,2.38,r'$+\;v_c\;-$',ha='center',fontsize=12)
resistor(7,1.5,r'$R_L$',True);wire([7,7,.5],[1.05,0,0])
wire([4.2,4.2],[3,2.65]);coil(4.2,2.2,r'$L_{k,ext}$',True);wire([4.2,4.2],[1.75,1.4])
ax.add_patch(Rectangle((3.47,.57),1.46,.83,facecolor='#e2f2f1',edgecolor='#007e80',lw=1.6))
ax.text(4.2,1.12,'Dominio SNSPD',ha='center',va='center',fontsize=10);ax.text(4.2,.81,r'$V_{dev}$',ha='center',va='center',fontsize=12)
wire([4.2,4.2],[.57,0]);ax.plot([4.2],[0],'o',color='#263c4a')
arrow((3.4,3.85),(4.0,3.85),r'$I_b$');arrow((3.16,2.62),(3.16,1.88),r'$I_s$',(-.28,0))
arrow((5.96,3.85),(6.64,3.85),r'$I_{rf}=I_b-I_s$')
ax.text(6.27,1.48,r'$V_{out}$',ha='center',fontsize=12);ax.text(6.27,2.5,'+',ha='center');ax.text(6.27,.45,'−',ha='center')
ax.text(2.4,-.42,r'$V_d=v_c+V_{out}$',ha='center',fontsize=13)
ax.text(5.7,-.42,r'$V_{out}=R_L(I_b-I_s)$',ha='center',fontsize=13)
OUT.parent.mkdir(parents=True,exist_ok=True);fig.savefig(OUT,dpi=190);plt.close(fig)
print(OUT)
