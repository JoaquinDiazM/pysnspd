"""Plot already measured differences; no RHS or new trajectories."""
from pathlib import Path
import csv,json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[3];DATA=ROOT/'docs/implementation/stage2/practical_review_20260922'
pair=json.loads((DATA/'raw/two_guarded_ph2049_pair_320_640.json').read_text())
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False,'savefig.dpi':185})
fig,ax=plt.subplots(1,2,figsize=(10,3.6),layout='constrained')
times=[row['time'] for row in pair['checkpoints']]
for name,label,col in [('phonon_population','Distribución fonónica','#007e80'),('electron_to_phonon','Transferencia e-fonón','#3274b2'),('amplitudes','Amplitud','#c7791c')]:
    vals=[100*row['errors'][name]['assessment_error'] for row in pair['checkpoints']]
    ax[0].plot(times,vals,color=col,label=label,lw=1.8)
ax[0].axhline(.0025,ls='--',color='#b74343',label='Presupuesto auxiliar: 0.0025%')
ax[0].set_yscale('symlog',linthresh=1e-7)
ax[0].set_ylim(0,.012)
ax[0].set(xlabel='Tiempo normalizado',ylabel='Diferencia temporal (%)',title='320 frente a 640 pasos, 2049 fonones')
ax[0].legend(fontsize=7.4,loc='lower right');ax[0].grid(alpha=.2)
rows=list(csv.DictReader((DATA/'comparison_metrics.csv').open()))
comparisons=sorted({r['comparison'] for r in rows if 'mesh' in r['comparison']})
assert len(comparisons)==1,comparisons
chosen={r['family']:float(r['relative_percent']) for r in rows if r['comparison']==comparisons[0] and r['cell']=='all' and r['relative_percent']!=''}
names=['amplitudes','phonon_energy','electron_to_phonon'];labels=['Amplitud','Energía fonónica','Transferencia e-fonón']
vals=[chosen[n] for n in names]
ax[1].barh(labels,vals,color=['#c7791c','#3274b2','#007e80'],height=.55)
for i,v in enumerate(vals):ax[1].text(v+.002,i,f'{v:.5f}%',va='center',fontsize=9)
ax[1].set(xlim=(0,.16),xlabel='Sensibilidad entre dos mallas (%)',title='1025 frente a 2049 fonones, 320 pasos')
fig.savefig(DATA/'figures/precision_context.png');plt.close(fig)
print('Saved-data figure written; no physical evaluation.')
