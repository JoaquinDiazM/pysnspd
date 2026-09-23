"""Analytical planning scales only: no Usadel queries or transient simulation."""
from pathlib import Path
from statistics import NormalDist
import hashlib,json,math
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[3]
OUT=ROOT/'docs/implementation/stage3_5/research_20260923'

def main():
    OUT.mkdir(parents=True,exist_ok=True)
    figures=OUT/'figures';figures.mkdir(exist_ok=True)
    width_nm=80.;sigma_nm=20.;D_nm2_ps=50.
    rows=[]
    for D in [50.,158.1]:
        for tau in [10.,20.,50.,100.,200.]:
            spread=math.sqrt(sigma_nm**2+2*D*tau)
            for epsilon in [.01,.001]:
                z=NormalDist().inv_cdf(1-epsilon/2)
                length=2*z*spread
                rows.append(dict(D_nm2_ps=D,tau_ps=tau,sigma_initial_nm=sigma_nm,
                    exterior_mass_fraction=epsilon,length_between_planes_nm=length,
                    length_over_width=length/width_nm))
    values=[]
    for D in [50.,158.1]:
        for ratio in [1.5,3.,4.,6.]:
            L=ratio*width_nm
            z=NormalDist().inv_cdf(.995)
            limit=((L/(2*z))**2-sigma_nm**2)/(2*D)
            values.append(dict(D_nm2_ps=D,L_over_W=ratio,
                maximum_tau_for_one_percent_mass_ps=max(0.,limit),
                initially_within_budget=limit>=0))
    hbar=1.054571817e-34;kB=1.380649e-23;Tc=8.65
    scales=[]
    for D in [5e-5,1.581e-4]:
        ell=math.sqrt(hbar*D/(2*kB*Tc))*1e9
        xi=math.sqrt(2)*ell
        scales.append(dict(D_m2_s=D,ell0_nm=ell,xi_c_nm=xi,
            normal_transverse_slowest_mode_ps=(80e-9)**2/(math.pi**2*D)*1e12,
            normal_width_diffusion_w2_over_4D_ps=(80e-9)**2/(4*D)*1e12,
            tau0_ps=math.pi*hbar/(8*kB*Tc)*1e12,
            hbar_over_BCS_gap_ps=hbar/(1.7638769888620456*kB*Tc)*1e12))
    cost=[]
    for L in [320,480,640,960]:
        for h in ([2,10,20] if L==640 else [10,20]):
            nx=math.ceil(L/h);ny=math.ceil(width_nm/h)
            nodes=(4*nx+1)*(4*ny+1)
            cost.append(dict(L_nm=L,element_nm=h,degree=4,nodes_2D=nodes,
                population_arrays_MiB=nodes*(630+1025)*8/2**20,
                meaning='One electron+phonon array pair only; no stages, spectra, Jacobians or runtime estimate.'))
    plt.rcParams.update({'font.size':11,'axes.spines.top':False,'axes.spines.right':False})
    time=np.linspace(0,200,401)
    fig,ax=plt.subplots(figsize=(7.3,4.2))
    colors=['#207c91','#b75b26']
    for D,color,label in [(50,colors[0],'D = 0,5 cm²/s (modelo K20)'),(158.1,colors[1],'D = 1,581 cm²/s (heredado)')]:
        L=2*NormalDist().inv_cdf(.995)*np.sqrt(sigma_nm**2+2*D*time)
        ax.plot(time,L/width_nm,color=color,label=label)
    ax.axhline(6,color='#777777',ls='--',label='6 anchos: ejemplo inicial')
    ax.set(xlabel='Tiempo desde la transferencia (ps)',ylabel='Distancia entre planos / ancho',
        title='La longitud necesaria depende de la ventana temporal')
    ax.legend(fontsize=9,loc='upper left');ax.grid(alpha=.2)
    fig.text(.5,.015,'Gaussiana libre: masa exterior 1 %, ancho 80 nm, σ inicial 20 nm supuesto.',ha='center',fontsize=9)
    fig.tight_layout(rect=[0,.05,1,1]);fig.savefig(figures/'01_length_window.png',dpi=220);plt.close(fig)
    fig,axes=plt.subplots(1,2,figsize=(7.3,3.5))
    x=np.linspace(-600,600,1201)
    for tau,color in [(20,'#207c91'),(100,'#b75b26')]:
        spread=math.sqrt(sigma_nm**2+2*D_nm2_ps*tau)
        y=np.exp(-.5*(x/spread)**2)/(math.sqrt(2*math.pi)*spread)
        axes[0].plot(x,y,label=f'{tau} ps',color=color)
    for sign in [-1,1]:axes[0].axvline(sign*240,color='#666666',ls='--')
    axes[0].set(xlabel='Distancia desde el impacto (nm)',ylabel='Densidad longitudinal (1/nm)',title='Mismo dominio: 6 anchos')
    axes[0].legend(fontsize=9)
    for eps,color in [(.01,'#207c91'),(.001,'#b75b26')]:
        vals=[r['length_over_width'] for r in rows if r['D_nm2_ps']==50 and r['exterior_mass_fraction']==eps]
        axes[1].plot([10,20,50,100,200],vals,'o-',label=f'Masa exterior {100*eps:g} %',color=color)
    axes[1].set(xlabel='Tiempo desde transferencia (ps)',ylabel='Distancia entre planos / ancho',title='Elegir el margen tiene coste')
    axes[1].legend(fontsize=8)
    for ax in axes:ax.grid(alpha=.2);ax.title.set_fontsize(11)
    fig.text(.5,.01,'Problema auxiliar lineal; no mide el error del voltaje ni el escape acumulado.',ha='center',fontsize=9)
    fig.tight_layout(rect=[0,.055,1,1]);fig.savefig(figures/'02_boundary_pedagogy.png',dpi=220);plt.close(fig)
    record=dict(status='ANALYTICAL_PLANNING_ONLY',physical_simulations=0,source='K20 supplementary table1 D=.5cm2/s; repository inventory inheritedD=1.581cm2/s; analytical Gaussian heat kernel',
        user_selected_width_nm=80,initial_sigma_nm=20,initial_sigma_status='ILLUSTRATIVE_ASSUMPTION_NOT_MEASUREMENT',
        epsilon_status='Illustrations, not adopted coupled-observable tolerance',
        distance_definition='Symmetric planes at +/-L/2 from deposition center; tau=t-t0',
        restrictions=['Exterior instantaneous mass is not accumulated first passage or lost heat.',
            'D is the normal diffusion parameter, not a certified upper bound on coupled spectral transport.',
            'No reaction, condensate, electrostatic or phonon feedback is included.',
            'L is the span to the external planes; does not by itself locate a valid 2D-to-1D interface.'],
        length_table=rows,fixed_length_time_limits=values,physical_scales=scales,memory_lower_bounds=cost,
        script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    (OUT/'geometry_scales.json').write_text(json.dumps(record,indent=2)+'\n',encoding='utf8')
    print(json.dumps(dict(scales=scales,time_limits=values,example_50ps=[r for r in rows if r['tau_ps']==50 and r['exterior_mass_fraction']==.01])))

if __name__=='__main__':main()
