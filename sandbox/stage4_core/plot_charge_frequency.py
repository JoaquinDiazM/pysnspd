"""Scientific result plot of the frozen-spectrum charge response diagnostic."""
from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[2]
DATA=ROOT/'docs/implementation/stage4/moment_review_20260924/charge_frequency'
rows=json.loads((DATA/'comparisons.json').read_text())
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False})
fig,axes=plt.subplots(1,2,figsize=(10,3.8),sharex=True)
for ax,key,title in zip(axes,('current','phase_torque_density'),('Corriente integrada','Torque de fase')):
    selected=sorted([r for r in rows if r['case_id']=='radial_65_N256' and r['eta_relative']==.01 and r['probe']=='angular' and r['moment']==key and r['nu']>0],key=lambda r:r['nu'])
    x=np.array([r['nu'] for r in selected])
    low=[];high=[]
    for nu in x:
        vals=[100*r['relative_dynamic_change'] for r in rows if r['probe']=='angular' and r['moment']==key and r['nu']==nu]
        low.append(min(vals));high.append(max(vals))
    ax.axvspan(.01,1,color='#eee8db',alpha=.75,zorder=0)
    ax.fill_between(x,low,high,color='#226f9c',alpha=.22)
    ax.loglog(x,[100*r['relative_dynamic_change'] for r in selected],'o-',color='#226f9c',label='Cambio dinámico frente a hT instantáneo')
    ax.loglog(x,[100*r['potential_relative_error'] for r in selected],'s--',color='#b44f27',label='Error al reducir hT(E) a un potencial')
    ax.set_title(title,fontweight='bold');ax.set_xlabel(r'$\nu=\omega\hbar/(2k_BT_c)$');ax.set_ylabel('Diferencia relativa (%)')
    ax.grid(True,which='major',alpha=.23);ax.set_xlim(.000075,1.3)
axes[0].legend(frameon=False,fontsize=8,loc='lower right')
fig.text(.5,.015,'Anclas lentas a la izquierda. Zona beige: exploración del cierre aproximado; respuesta AC no acreditada.',ha='center',fontsize=9,color='#554b3c')
fig.suptitle('La compresión espectral domina el error en las anclas lentas',fontsize=13,fontweight='bold')
fig.tight_layout(rect=(0,.06,1,.93))
fig.savefig(DATA/'frequency_response.png',dpi=180);fig.savefig(DATA/'frequency_response.pdf')
print(str(DATA/'frequency_response.png'))
