"""Plot measured accuracy and cost without implying complete stage admission."""
from pathlib import Path
import hashlib,json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter
ROOT=Path(__file__).resolve().parents[3]
DATA=ROOT/'docs/implementation/stage2/resume_20260921'

def main():
    result=json.loads((DATA/'short_time_assessment.json').read_text())
    steps=np.asarray([r['steps'] for r in result['cases']])
    runs=[json.loads((DATA/f'short_time/one_rk4_{n}.json').read_text()) for n in steps]
    ref=json.loads((DATA/'manual_reference/one_reference_probe.json').read_text())
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
    fig,axes=plt.subplots(1,2,figsize=(10.4,4.8))
    for key,label,color,marker in (
        ('electronic_population','Distribución electrónica','#2166ac','o'),
        ('phonon_population','Distribución fonónica','#d66024','s'),
        ('electron_to_phonon','Transferencia e → fonones','#00897b','^')):
        values=[r['errors'][key]['assessment_error']*100 for r in result['cases']]
        axes[0].loglog(steps,values,color=color,marker=marker,label=label,lw=1.6)
    axes[0].axhline(1e-2,color='#a61c35',ls=':',label='Límite temporal: 0,01 %')
    axes[0].set(xticks=steps,xticklabels=[str(n) for n in steps],xlabel='Pasos RK4 en el intervalo 0–0,1',
        ylabel='Error relativo máximo (%)',title='Precisión frente a DOP853')
    axes[0].grid(alpha=.2,which='both')
    axes[0].xaxis.set_minor_formatter(NullFormatter())
    fig.legend(*axes[0].get_legend_handles_labels(),fontsize=8,loc='lower center',
               bbox_to_anchor=(.5,.065),ncol=4,frameon=False)
    elapsed=[r['runtime_seconds'] for r in runs]+[ref['runtime_seconds']]
    names=[f'RK4 · {n} pasos' for n in steps]+['DOP853 · referencia']
    bars=axes[1].barh(names,elapsed,color=['#8bb4d5','#518bb8','#2166ac','#69757e'],height=.6)
    axes[1].set_xscale('log');axes[1].invert_yaxis()
    for bar,value in zip(bars,elapsed):
        axes[1].text(value*1.12,bar.get_y()+bar.get_height()/2,f'{value:.1f} s',va='center',fontsize=9)
    axes[1].set(xlim=(1,3500),xlabel='Tiempo de cálculo (s; escala logarítmica)',
        title='Coste medido en Geminga · un hilo')
    axes[1].grid(axis='x',alpha=.2)
    fig.suptitle('Una celda · 630 estados electrónicos / 1025 fonónicos',fontsize=13,fontweight='bold',y=.98)
    fig.text(.06,.025,'Comparación en t = 0; 0,05; 0,1. Datos sintéticos. Este ensayo corto no cierra la etapa 2.',fontsize=9,color='#555d65')
    fig.subplots_adjust(left=.08,right=.96,top=.83,bottom=.27,wspace=.68)
    folder=DATA/'figures';folder.mkdir(exist_ok=True)
    for extension in ('png','svg'):
        fig.savefig(folder/f'short_time_accuracy_cost.{extension}',dpi=190)
    svg=folder/'short_time_accuracy_cost.svg'
    svg.write_text('\n'.join(line.rstrip() for line in svg.read_text(encoding='utf-8').splitlines())+'\n',encoding='utf-8')
    inputs=['short_time_assessment.json','manual_reference/one_reference_probe.json']+[f'short_time/one_rk4_{n}.json' for n in steps]
    record={'inputs':{p:hashlib.sha256((DATA/p).read_bytes()).hexdigest() for p in inputs},
            'runner_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            'outputs':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(folder.glob('short_time_accuracy_cost.*'))}}
    (folder/'short_time_accuracy_cost.json').write_text(json.dumps(record,indent=2)+'\n')
    print(json.dumps(dict(status=result['status'],fastest_seconds=elapsed[0],
        finest_seconds=elapsed[2],reference_seconds=elapsed[3],reference_to_finest_speedup=elapsed[3]/elapsed[2]),indent=2))

if __name__=='__main__':main()
