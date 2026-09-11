"""Result plots from the independent assessment; no physics is recalculated."""
from pathlib import Path
import argparse,hashlib,json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[2]
DATA=ROOT/'docs/implementation/stage1_r2'

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('assessment',type=Path)
    args=parser.parse_args()
    assessment=json.loads(args.assessment.read_text(encoding='utf-8'))
    criteria=json.loads((DATA/'acceptance_criteria.json').read_text(encoding='utf-8'))
    previous=json.loads((ROOT/'docs/implementation/stage1/electronic/boundary_comparison.json').read_text(encoding='utf-8'))
    before=next(r for r in previous['rows'] if r['delta']==.72 and r['gamma']==0)
    candidates=[r for r in assessment['physical_rows'] if r.get('amplitude')==.72
                and r.get('gamma')==0 and r.get('profile')=='low_energy']
    if len(candidates)!=1: raise ValueError('One unique critical case is required.')
    current=candidates[0]
    reference=criteria['critical_analytic_reference']['Gamma_conjugate_value']
    errors=100*np.abs(np.array([before['new_gamma_conjugate'],current['candidate'][2]])-reference)/abs(reference)
    gates=criteria['numerical_gates']
    plt.rcParams.update({'font.size':11,'axes.spines.top':False,'axes.spines.right':False})
    blue,teal,orange='#2463a6','#008b86','#d05f17'
    fig,ax=plt.subplots(figsize=(8.6,3.3),layout='constrained')
    ax.plot([0,1],errors,'--',color='#a7b4bf',lw=1.2)
    ax.scatter([0,1],errors,color=[blue,teal],s=90,zorder=4)
    for i,value in enumerate(errors):
        ax.annotate(f'{value:.5f}%',(i,value),xytext=(15,8),textcoords='offset points',fontsize=12,weight='bold')
    ax.axhline(.1,color=orange,lw=1.5,ls=':',label='Tolerancia fijada: 0,1%')
    ax.set_yscale('log')
    ax.set_ylim(.02, 30)
    ax.set_xticks([0,1],['R1: catálogo entregado','R2: catálogo evaluado'])
    ax.set_xlim(-.35,1.6)
    ax.set_ylabel('Error total de la pendiente j/q (%)')
    ax.set_title('Caso crítico: comparación con el mismo límite causal')
    ax.legend(loc='upper right',fontsize=11)
    ax.grid(axis='y',which='major',alpha=.16)
    figures=DATA/'figures';figures.mkdir(exist_ok=True)
    for ext in ('png','pdf'):fig.savefig(figures/f'critical_before_after.{ext}',dpi=180)
    plt.close(fig)

    fig,ax=plt.subplots(figsize=(9.5,4.25),layout='constrained')
    groups=[('scaled_error',0,gates['energy_max_scaled_error'],'Energía\nescalada'),
            ('scaled_error',1,gates['amplitude_force_max_scaled_error'],'Fuerza\nescalada'),
            ('scaled_error',2,gates['gamma_response_max_scaled_error'],'Respuesta Γ\nescalada'),
            ('relative_error',1,gates['nonzero_force_and_gamma_max_relative_error'],'Fuerza\nrelativa'),
            ('relative_error',2,gates['nonzero_force_and_gamma_max_relative_error'],'Respuesta Γ\nrelativa')]
    rng=np.random.default_rng(302)
    summaries=[]
    for k,(key,index,limit,label) in enumerate(groups):
        values=np.array([r[key][index]/limit for r in assessment['physical_rows'] if r[key][index] is not None])
        # Exact zeros are explicitly identified and drawn at a stated floor,
        # never passed off as a measured nonzero error.
        visible=np.maximum(values,1e-10)
        ax.scatter(k+rng.uniform(-.18,.18,len(values)),visible,s=12,alpha=.45,color=blue)
        maximum=float(values.max())
        ax.scatter([k],[max(maximum,1e-10)],marker='D',s=48,color=orange,zorder=4)
        summaries.append({'group':label.replace('\n',' '),'samples':len(values),'maximum_fraction_of_tolerance':maximum,'exact_zeros':int(np.sum(values==0))})
    ax.axhline(1,color='#a13636',ls='--',lw=1.5,label='Límite de aceptación')
    ax.set_yscale('log');ax.set_ylim(5e-11,max(3,10*max(r['maximum_fraction_of_tolerance'] for r in summaries)))
    ax.set_xticks(range(len(groups)),[r[3] for r in groups])
    ax.set_ylabel('Error / tolerancia aplicable')
    ax.set_title('Cada punto es un caso; el rombo marca el peor resultado')
    ax.legend(loc='upper left',fontsize=11)
    ax.grid(axis='y',which='major',alpha=.16)
    for ext in ('png','pdf'):fig.savefig(figures/f'acceptance_cases.{ext}',dpi=180)
    plt.close(fig)
    (DATA/'plotted_results.json').write_text(json.dumps({
        'assessment_sha256':hashlib.sha256(args.assessment.read_bytes()).hexdigest(),
        'assessment_status':assessment['status'],
        'before_after_critical_total_error_percent':errors.tolist(),
        'critical_reference':reference,'current_critical_case':current,
        'groups':summaries,'plotting_floor_fraction_of_tolerance':1e-10,
        'zero_display_policy':'Exact zero errors share the stated plotting floor; this does not change acceptance.',
        'scaled_definition':'abs(candidate-reference)/max(1,abs(reference))',
        'historical_comparison':'R1 is the delivered Hermite catalogue at commit d17d7c3; its total error differs from its independently resolved 8.61% regulator-only bias.'
    },indent=2)+'\n',encoding='utf-8')

if __name__=='__main__':main()
