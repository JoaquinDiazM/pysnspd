"""Plot the measured two-cell trajectories, without evaluating dynamics."""
from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter
ROOT=Path(__file__).resolve().parents[3]
DATA=ROOT/'docs/implementation/stage2/closure_prep_20260922'
FIG=DATA/'figures';FIG.mkdir(exist_ok=True)
read=lambda p:json.loads(p.read_text(encoding='utf-8'))
runs={n:read(DATA/f'raw_two/two_ssp_{n}.json') for n in (160,320,640,1280)}
audit=read(DATA/'two_cell_audit.json');one=read(ROOT/'docs/implementation/stage2/time_pass_20260922/one_cell_audit.json')
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.titlesize':11,'axes.spines.top':False,'axes.spines.right':False,'savefig.dpi':190})
colors=['#cc6e13','#3274b2','#007f73'];steps=np.array([160,320,640])
fig,ax=plt.subplots(1,2,figsize=(10,3.8),layout='constrained')
for data,label,color,marker in [(one,'Una celda',colors[1],'o'),(audit,'Dos celdas',colors[2],'s')]:
    errors=[r['worst']['error'] for r in data['rows']]
    ax[0].plot(steps,errors,marker+'-',color=color,label=label,lw=2)
ax[0].axhline(1e-4,ls='--',color='#b23737',label='Criterio temporal: 1e-4')
ax[0].axhline(2.5e-5,ls=':',color='#777777',label='Presupuesto de mallas: 2.5e-5')
ax[0].set(yscale='log',xscale='log',ylim=(1.7e-6,2e-4),title='Precisión frente a 1280 pasos',ylabel='Máximo error relativo',xlabel='Número de pasos')
ax[0].set_xticks(steps,steps);ax[0].legend(fontsize=8,loc='lower left')
ledger=[runs[n]['energy_ledger_scaled_max'] for n in (*steps,1280)]
ax[1].plot([160,320,640,1280],ledger,'s-',lw=2,color=colors[2],label='Dos celdas')
ax[1].axhline(1e-7,ls='--',color='#b23737',label='Criterio energético: 1e-7')
ax[1].set(yscale='log',xscale='log',ylim=(1e-11,3e-7),title='Balance de energía con baño y fuente',ylabel='Máximo defecto escalado',xlabel='Número de pasos')
ax[1].set_xticks([160,320,640,1280],[160,320,640,1280]);ax[1].legend(fontsize=8,loc='lower left')
for a in ax:a.xaxis.set_minor_formatter(NullFormatter());a.grid(True,alpha=.2,which='both')
fig.savefig(FIG/'two_time_convergence.png');plt.close(fig)

ref=runs[1280];t=np.array(ref['times']);sn=ref['snapshots']
series=lambda name:np.array([s[name] for s in sn])
fig,ax=plt.subplots(2,2,figsize=(10,6.1),layout='constrained')
for i,col in enumerate(colors[:2]):
    ax[0,0].plot(t,series('amplitudes')[:,i],color=col,lw=2,label=f'Celda {i+1}, referencia')
    coarse=runs[320];tc=np.array(coarse['times']);sc=np.array([s['amplitudes'][i] for s in coarse['snapshots']])
    ax[0,0].plot(tc[::20],sc[::20],ls='none',marker='o',ms=4,mfc='white',mec=col)
    ax[1,0].plot(tc,1e8*(sc-series('amplitudes')[::4,i]),color=col,label=f'Celda {i+1}')
ax[0,0].set(title='Condensados inicialmente distintos',ylabel=r'$|\Delta|/\Delta_0$');ax[0,0].legend(fontsize=8)
ax[1,0].set(title='Diferencia SSP320 menos SSP1280',ylabel=r'Diferencia de amplitud ($10^{-8}$)');ax[1,0].legend(fontsize=8)
ax[1,0].axhline(0,color='#777777',lw=.6)
for name,label,col in [('transport','Transporte acumulado',colors[2]),('input','Fuente externa acumulada',colors[0]),('escape','Escape acumulado',colors[1])]:
    data=series(name)
    if data.ndim>1:data=data.sum(axis=1)
    ax[0,1].plot(t,data,color=col,lw=2,label=label)
ax[0,1].set(title='Canales activos y medidos',ylabel=r'Energía / $(N_0\Delta_0^2)$');ax[0,1].legend(fontsize=8)
activity=next(row for row in audit['trajectories'] if row['steps']==320)['activity']
fractions=[100*activity['electron_to_phonon_fraction'],100*activity['face_energy_fraction']]
ax[1,1].barh(['Electrones-fonones','Entre celdas'],fractions,color=[colors[1],colors[2]],height=.48)
ax[1,1].axvline(1,ls='--',color='#b23737',label='Actividad mínima: 1%')
for i,v in enumerate(fractions):ax[1,1].text(v+.15,i,f'{v:.2f}%',va='center',fontsize=10)
ax[1,1].set(title='La prueba resuelve intercambio real',xlabel='Transferencia / excitación inicial (%)',xlim=(0,13))
ax[1,1].legend(fontsize=8,loc='lower right')
for a in [ax[0,0],ax[0,1],ax[1,0]]:a.set_xlabel(r'Tiempo $t/t_{\mathrm{ref}}$');a.grid(True,alpha=.2)
fig.savefig(FIG/'two_cell_dynamics.png');plt.close(fig)
diagnosis=read(DATA/'isolated_transport_diagnostic/isolated_transport_diagnosis.json')
guard=read(DATA/'isolated_guarded/isolated_results.json')
escape=read(DATA/'isolated_transport_diagnostic/isolated_escape_order.json')
witness=diagnosis['rejected_electrons'][0]
quantum=np.nextafter(0.,1.)
values=[witness['p_before']/quantum,witness['p_rejected']/quantum,
        guard['evidence']['previous_rejected_stage']['witness_guarded_value']/quantum]
fig,ax=plt.subplots(1,2,figsize=(10,3.9),layout='constrained')
ax[0].bar(['Entrada','SSP anterior','SSP protegido'],values,color=['#657888','#b23737',colors[2]],width=.58)
ax[0].axhline(0,color='#333333',lw=.8)
for i,v in enumerate(values):ax[0].text(i,v+(2 if v>=0 else -3),f'{v:.0f} q',ha='center',va='bottom' if v>=0 else 'top')
ax[0].set(ylim=(-17,84),title='La misma etapa que produjo el rechazo',ylabel='Ocupación / q',xlabel='q = 4.94e-324 (mínimo float64 positivo)')
ns=[10,20,40];err=escape['errors_capacity_L1']
ax[1].plot(ns,err,'o-',color=colors[1],lw=2,label='Error frente a solución exacta')
ax[1].plot(ns,err[0]*(np.array(ns)/10.)**-3,'--',color=colors[0],label='Pendiente de orden 3')
ax[1].set(xscale='log',yscale='log',title='Escape: reducción esperada de 8 veces',xlabel='Número de pasos',ylabel='Error L1 ponderado')
ax[1].set_xticks(ns,ns);ax[1].xaxis.set_minor_formatter(NullFormatter());ax[1].legend(fontsize=8);ax[1].grid(True,alpha=.2)
fig.savefig(FIG/'guarded_rounding_and_escape.png');plt.close(fig)
print('Three measured-result figures prepared.')
