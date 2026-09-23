"""Summarize digitized historical fractions and phonon percentiles only."""
from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[3]
OUT=ROOT/'docs/implementation/stage3_5/assessment_r2_20260923'

def main():
    v=json.loads((OUT/'cascade/cascade_partition.json').read_text(encoding='utf8'))
    rows=v['finite_radius_and_phonon_percentiles'];x=np.arange(len(rows))
    fig,axes=plt.subplots(1,2,figsize=(8,3.8))
    e=np.array([r['electron_global_fraction_bound'] for r in rows])*100
    mid=e.mean(axis=1)
    axes[0].bar(x,mid,color='#b75b26',label='Electrones')
    axes[0].bar(x,100-mid,bottom=mid,color='#207c91',label='Fonones')
    axes[0].errorbar(x,mid,yerr=np.array([mid-e[:,0],e[:,1]-mid]),fmt='none',ecolor='black',capsize=5)
    for i,r in enumerate(rows):
        axes[0].text(i,35,'~92 %',ha='center',color='white',fontsize=12)
        axes[0].text(i,15,f'{e[i,0]:.1f}–{e[i,1]:.1f} % e',ha='center',fontsize=9)
    axes[0].set_title('La energía no es sólo fonónica');axes[0].set_ylabel('Porcentaje de energía del fotón')
    axes[0].set_ylim(0,112);axes[0].legend(loc='upper center',ncol=2,fontsize=8)
    ratio=np.array([r['phonon_R90_over_R50_envelope'] for r in rows]);m=ratio.mean(axis=1)
    axes[1].errorbar(x,m,yerr=np.array([m-ratio[:,0],ratio[:,1]-m]),fmt='o',color='#207c91',capsize=6,label='Perfil fonónico A20')
    axes[1].axhline(v['gaussian2D_R90_over_R50'],color='#b75b26',linestyle='--',label='Cualquier gaussiana 2D')
    axes[1].set_ylim(1.6,2.7);axes[1].set_ylabel('Radio del 90 % / radio del 50 %')
    axes[1].set_title('El perfil fonónico también tiene cola');axes[1].legend(fontsize=8,loc='upper left')
    for ax in axes:
        ax.set_xticks(x,[f"{r['time_ps']:.4f}" for r in rows]);ax.set_xlabel('Tiempo del cálculo A20 (ps)')
        ax.spines[['top','right']].set_visible(False);ax.grid(axis='y',alpha=.15)
    fig.text(.5,.012,'Datos leídos de Allmaras, figuras 2.6/2.8. Barras: lectura y cola no mostrada; no intervalos físicos de Korzh.',ha='center',fontsize=8)
    fig.tight_layout(rect=[0,.05,1,1]);(OUT/'figures').mkdir(exist_ok=True)
    fig.savefig(OUT/'figures/cascade_partition_shape.png',dpi=220);plt.close(fig)
    print('Cascade summary figure created; no new trajectories.')

if __name__=='__main__':main()
