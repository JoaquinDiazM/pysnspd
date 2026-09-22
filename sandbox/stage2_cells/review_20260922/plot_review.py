"""Result plots from archived trajectories only; no RHS or integration calls."""
from pathlib import Path
import json,hashlib
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[3]
DATA=ROOT/'docs/implementation/stage2/review_20260922'
OUT=DATA/'figures'
BLUE='#2369a1';ORANGE='#dc7737';TEAL='#07817c';RED='#b33840';GRAY='#596775'
def read(path):return json.loads(Path(path).read_text(encoding='utf-8'))
def snapshots(record,key):return np.asarray([r[key] for r in record['snapshots']])
def save(fig,name):
    fig.savefig(OUT/(name+'.png'),dpi=220,bbox_inches='tight',facecolor='white')
    plt.close(fig)
def main():
    OUT.mkdir(parents=True,exist_ok=True)
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.titlesize':12,
        'axes.labelsize':10,'axes.spines.top':False,'axes.spines.right':False,
        'axes.grid':True,'grid.alpha':.18,'legend.frameon':False,'savefig.pad_inches':.14})
    one=read(DATA/'raw/rk4/one_time_320.json');two=read(DATA/'raw/rk4/two_time_320.json')
    ssp=read(DATA/'raw/ssp/one_ssp_40.json')
    # Distinct symbols/line styles plus the numerical difference expose agreement.
    fig,axes=plt.subplots(1,2,figsize=(10.0,4.25),layout='constrained')
    for case,color,marker in [('one',BLUE,'o'),('two',ORANGE,'s')]:
        result=read(DATA/f'raw/rk4/{case}_time_assessment.json')
        steps=[r['steps'] for r in result['cases']];errors=[r['max_error'] for r in result['cases']]
        axes[0].loglog(steps,errors,color=color,marker=marker,lw=2,label='Una celda' if case=='one' else 'Dos celdas')
        axes[0].annotate(f'{errors[-1]:.2e}',(steps[-1],errors[-1]),xytext=(8,3 if case=='one' else -3),
            textcoords='offset points',color=color,fontsize=9)
    axes[0].axhline(1e-4,color=RED,ls='--',label='Límite: 1e-4')
    axes[0].set(xlabel='Pasos en el intervalo 0 a 2',ylabel='Error relativo máximo',title='A. Precisión temporal RK4')
    axes[0].set_xticks([40,80,160],['40','80','160']);axes[0].set_xticks([],minor=True)
    axes[0].set_xlim(36,218);axes[0].set_ylim(3.8e-5,6.8e-4)
    axes[0].legend(loc='upper right',fontsize=9)
    counts=np.arange(4);labels=['RK4\n1 celda\n160','RK4\n2 celdas\n160','SSP\n1 celda\n40','SSP\nintervalo\nfino']
    probe=read(ROOT/'docs/implementation/stage2/recovery_20260921/limited_continuous_first_interval.json')
    records=[read(DATA/'raw/rk4/one_time_160.json'),read(DATA/'raw/rk4/two_time_160.json'),ssp,probe]
    ratios=[r['energy_ledger_scaled_max']/1e-7 for r in records]
    axes[1].bar(counts,ratios,color=[BLUE,ORANGE,RED,TEAL],width=.62)
    axes[1].axhline(1,color=RED,ls='--');axes[1].set_yscale('log')
    axes[1].set_xticks(counts,labels);axes[1].set(ylabel='Error energético / límite',title='B. Balance acumulado')
    axes[1].set_ylim(1e-4,300)
    for x,y in zip(counts,ratios):axes[1].text(x,y*1.3,f'{y:.3g}×',ha='center',fontsize=9)
    axes[1].text(.02,.95,'Menor que 1: pasa',transform=axes[1].transAxes,va='top',color=GRAY,fontsize=9)
    save(fig,'acceptance_overview')
    # Compare actual saved states. The reference is temporal, not a material result.
    t=np.asarray(two['times']);a=snapshots(two,'amplitudes');ph=snapshots(two,'phonon_energy')
    fig,axes=plt.subplots(2,2,figsize=(10.0,6.1),layout='constrained')
    for i,(color,style) in enumerate([(ORANGE,'-'),(BLUE,'--')]):
        axes[0,0].plot(t,a[:,i],color=color,ls=style,lw=2,label=f'Celda {i+1}')
        axes[0,1].plot(t,ph[:,i],color=color,ls=style,lw=2,label=f'Celda {i+1}')
    axes[0,0].set(title='A. Condensado',ylabel=r'$|\Delta|/\Delta_0$');axes[0,0].legend()
    axes[0,1].set(title='B. Energía fonónica',ylabel=r'$U_{ph}/(N_0\Delta_0^2)$');axes[0,1].legend()
    for field,label,color,style in [('transport','Transporte interno',TEAL,'-'),
            ('electron_to_phonon','Electrón-fonón neto',BLUE,'--'),('input','Fuente externa',ORANGE,':')]:
        v=snapshots(two,field)
        if v.ndim>1:v=v.sum(axis=1)
        axes[1,0].plot(t,v,label=label,color=color,ls=style,lw=2)
    axes[1,0].set(title='C. Transferencias acumuladas',ylabel=r'$Q/(N_0\Delta_0^2)$');axes[1,0].legend(fontsize=8)
    with np.load(DATA/'raw/rk4/two_time_320.npz') as archive:
        states=archive['states'];x=archive['electron_count'];ne=len(x);block=1+ne+len(archive['phonon_capacities'])
        for i,(color,style) in enumerate([(ORANGE,'-'),(BLUE,'--')]):
            start=i*block+1;delta=states[-1,start:start+ne]-states[0,start:start+ne]
            axes[1,1].plot(x,delta,color=color,ls=style,lw=1.8,label=f'Celda {i+1}')
    axes[1,1].axhline(0,color=GRAY,lw=.6);axes[1,1].set(xlim=(0,3),
        title='D. Cambio de ocupaciones',xlabel=r'Coordenada de conteo $x/\Delta_0$',ylabel=r'$p(x,2)-p(x,0)$')
    axes[1,1].legend(fontsize=8)
    for ax in axes.flat:
        if ax is not axes[1,1]:ax.set_xlabel(r'Tiempo $t/t_{ref}$')
    save(fig,'two_cell_results')
    fig,axes=plt.subplots(1,2,figsize=(10.0,4.2),layout='constrained')
    for r,label,color,style in [(ssp,'SSP, 40 pasos',RED,'-'),
            (read(DATA/'raw/rk4/one_time_40.json'),'RK4, 40 pasos',BLUE,'--')]:
        balance=snapshots(r,'conserved');residual=balance-balance[0]
        axes[0].plot(r['times'],residual/1e-7,color=color,ls=style,lw=2,label=label)
    axes[0].axhspan(-1,1,color=TEAL,alpha=.15,label='Banda admitida')
    axes[0].set(xlabel=r'Tiempo $t/t_{ref}$',ylabel='Residuo energético / límite',title='A. Mismo ensayo completo, distinto integrador')
    axes[0].legend(fontsize=8)
    fine_steps=np.array([10,20,40]);values=[]
    for count in fine_steps:
        r=read(ROOT/f'docs/implementation/stage2/recovery_20260921/one_ssp_short_{count}.json')
        values.append(r['energy_ledger_scaled_max'])
    axes[1].loglog(fine_steps,values,'o-',color=TEAL,lw=2,label='SSP: medido en 0 a 0,1')
    axes[1].loglog(fine_steps,values[0]*(10/fine_steps)**3,'--',color=GRAY,label='Guía de tercer orden')
    axes[1].set_xticks(fine_steps,[str(n) for n in fine_steps])
    axes[1].set(xlabel='Pasos en el ensayo corto',ylabel='Error energético escalado',title='B. El paso sí importa')
    axes[1].legend(fontsize=8)
    save(fig,'ssp_failure_diagnosis')
    receipt={'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'scope':'Postprocessing of stored trajectories only; no RHS evaluation, interpolation of populations or new dynamics.',
        'plots':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(OUT.glob('*.png'))},
        'ledger_over_tolerance':dict(zip(['rk4_one_160','rk4_two_160','ssp_one_40','ssp_fine_short'],ratios))}
    (DATA/'plot_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(receipt,indent=2))
if __name__=='__main__':main()
