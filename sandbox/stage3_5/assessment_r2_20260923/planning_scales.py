"""Algebraic mobility/domain planning; no new spectra or physical trajectories."""
from pathlib import Path
from statistics import NormalDist
import csv,hashlib,json,math
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[3]
OUT=ROOT/'docs/implementation/stage3_5/assessment_r2_20260923'
KB=1.380649e-23;HBAR=1.054571817e-34;TC=8.65;TB=.9
DELTA=1.7638769888620456*KB*TC
PAIRS={'memory_effective':(.5,2.47),'Allmaras_reference':(5.,24.7),'Korzh_reference':(6.,24.7)}

def mobility(T,amplitude,pair):
    t=max(T,TB)/TC
    tee,tep=pair
    tau=1/(t/tee+t**3/tep)
    R=math.sqrt(1+(2*DELTA*tau*1e-12/HBAR*amplitude)**2)
    tau0=math.pi*HBAR/(8*KB*TC)*1e12
    A=math.sqrt((1+t)/2)
    return dict(taupsi_ps=tau,R=R,radial_per_ps=1/(2*A*tau0*R),
                tangential_per_ps=R/(2*A*tau0))

def main():
    OUT.mkdir(parents=True,exist_ok=True);figdir=OUT/'figures';figdir.mkdir(exist_ok=True)
    rows=[]
    for name,pair in PAIRS.items():
        for T in [TB,TC/2,TC]:
            for amplitude in [0.,.1,.5,1.]:
                rows.append(dict(scenario=name,T_K=T,amplitude_bar=amplitude,**mobility(T,amplitude,pair)))
    with (OUT/'mobility_comparison.csv').open('w',newline='',encoding='utf8') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    temperature=np.linspace(TB,TC,250)
    fig,axes=plt.subplots(1,2,figsize=(8,3.6),sharex=True)
    for key,label,color in [('Allmaras_reference','Allmaras / heredado','#007c83'),('Korzh_reference','Korzh / heredado','#b75b26')]:
        for i,component in enumerate(['radial_per_ps','tangential_per_ps']):
            ratio=[mobility(T,1.,PAIRS[key])[component]/mobility(T,1.,PAIRS['memory_effective'])[component] for T in temperature]
            axes[i].plot(temperature,ratio,color=color,label=label)
    axes[0].set_title('Respuesta de la amplitud');axes[1].set_title('Respuesta tangencial (fase)')
    for ax in axes:
        ax.set_xlabel('Temperatura usada por KWT (K)');ax.set_ylabel('Movilidad / movilidad heredada')
        ax.legend(fontsize=8);ax.grid(alpha=.2);ax.spines[['top','right']].set_visible(False)
    fig.text(.5,.012,'Comparación algebraica a |Δ|/Δ₀ = 1 y fuerza cartesiana fija; no predice latencias.',ha='center',fontsize=9)
    fig.tight_layout(rect=[0,.05,1,1]);fig.savefig(figdir/'mobility_components.png',dpi=220);plt.close(fig)
    # Geometric planning is parametrized by source extent, not a new width choice.
    z=NormalDist().inv_cdf(.995)
    domain=[]
    for time in [20,50,100,200]:
        for length in [320,480,640,960]:
            room=(length/(2*z))**2-2*50*time
            domain.append(dict(window_ps=time,L_nm=length,L_over_width=length/80,
                largest_initial_sigma_for_linear_one_percent_filter_nm=math.sqrt(room) if room>=0 else None,
                source_width_selected=False))
    result={'schema':'pysnspd.stage3_5.algebraic_planning.r2','new_physical_transients':0,
      'code_parameters_changed':False,'mobility_pairs_ps':PAIRS,'mobility_rows':rows,
      'mobility_units':'inverse picoseconds per dimensionless real Cartesian energy gradient',
      'interpretation':'Historical mobility pairs act differently on amplitude and phase; no global time-rescaling equivalence.',
      'domain_prefilter':domain,'domain_scope':'Auxiliary normal free diffusion D50nm2/ps; instantaneous outside mass is not accumulated escape or coupled-observable error.',
      'width_nm':None,'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    (OUT/'planning_scales.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf8')
    print(json.dumps({'mobility_rows':len(rows),'domain_cases':len(domain),'new_transients':0}))

if __name__=='__main__':main()
