"""Pedagogical equal-energy profiles; not a detector simulation."""
from pathlib import Path
import json,hashlib
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[3]
OUT=ROOT/'docs/implementation/stage3_5/research_20260923'

def main():
    x=np.linspace(0,2.5,1001)
    profiles=[('Disco uniforme',np.where(x<=1,1.,0.),np.minimum(x*x,1.),'#707880','-'),
      ('Gauss: mismo segundo momento',2*np.exp(-2*x*x),1-np.exp(-2*x*x),'#207c91','-'),
      ('Gauss: mismo máximo',np.exp(-x*x),1-np.exp(-x*x),'#b75b26','--')]
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
    fig,axes=plt.subplots(1,2,figsize=(7.3,3.8))
    for label,density,cumulative,color,style in profiles:
        axes[0].plot(x,density,color=color,ls=style,label=label)
        axes[1].plot(x,cumulative,color=color,ls=style)
    axes[0].set(xlabel='Distancia al centro / radio del disco',ylabel='Densidad / densidad del disco',title='La forma cambia la concentración')
    axes[1].set(xlabel='Radio de integración / radio del disco',ylabel='Fracción de energía encerrada',title='La energía total es la misma',ylim=(0,1.05))
    for ax in axes:ax.grid(alpha=.2);ax.title.set_fontsize(11)
    fig.legend(*axes[0].get_legend_handles_labels(),loc='lower center',bbox_to_anchor=(.5,.015),ncol=1,fontsize=9)
    fig.tight_layout(rect=[0,.2,1,1]);fig.savefig(OUT/'figures/03_disk_gaussian.png',dpi=220);plt.close(fig)
    data=dict(status='ANALYTICAL_PROFILE_IDENTITY',new_physical_transients=0,
       disk=dict(energy_normalization='E/(pi R^2)',mean_squared_radius='R^2/2'),
       gaussian=dict(energy_normalization='E/(2 pi s^2)',mean_squared_radius='2s^2',
         enclosed_fraction='1-exp(-r^2/(2s^2))',R90_over_s=float(np.sqrt(2*np.log(10)))),
       same_energy_and_second_moment=dict(s_over_R=.5,peak_over_disk=2),
       same_energy_and_peak=dict(s_over_R=float(1/np.sqrt(2)),second_moment_over_disk=2),
       meaning='Conversion conventions, not empirical calibration or proof of identical forces.',
       script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    (OUT/'profile_comparison.json').write_text(json.dumps(data,indent=2)+'\n',encoding='utf8')
    print('Equal-energy profile comparison saved.')

if __name__=='__main__':main()
