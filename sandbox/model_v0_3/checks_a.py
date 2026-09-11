"""Lightweight uniform Usadel checks for document A, independent of production.

Energy unit k_B T_c, Q=q sqrt(hbar D/(2 k_B T_c)), Gamma/(k_B T_c)=Q^2.
Only equilibrium algebraic problems are solved. No detector transient is run.
"""
from pathlib import Path
import json
import time
import platform
import numpy as np
from scipy.optimize import brentq, minimize_scalar
from scipy.special import zeta
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / 'docs/modelo_v0_3'
FIG = OUT / 'figuras'
DATA = OUT / 'verificaciones'
for folder in (FIG, DATA):
    folder.mkdir(parents=True, exist_ok=True)
plt.rcParams.update({'font.size': 11, 'axes.spines.top': False, 'axes.spines.right': False,
                     'figure.dpi': 150, 'savefig.dpi': 200, 'font.family': 'DejaVu Sans'})
COLORS = ['#146a8a', '#bf6334', '#59894b']
KB, HBAR, ECH = 1.380649e-23, 1.054571817e-34, 1.602176634e-19
TC, TEMP, D, SIGMA, AREA = 8.65, .9, 1.581e-4, 4.2e5, 120e-9*7e-9
T = TEMP/TC
ENERGY = KB*TC
QL = np.sqrt(HBAR*D/(2*ENERGY))
IPREF = 2*np.pi*SIGMA*KB*TEMP/ECH/QL*AREA*1e6
XI2 = np.pi*HBAR*D/(4*np.sqrt(2)*KB*TC*np.sqrt(1+T))
BCS_AL = 1.764*np.tanh(1.74*np.sqrt(TC/TEMP-1))

def spectrum(delta, gamma, n=800):
    eps = np.pi*T*(2*np.arange(n)+1)
    theta = np.arctan(delta/(eps+gamma))
    for _ in range(25):
        s, c = np.sin(theta), np.cos(theta)
        residual = delta*c-(eps+gamma*c)*s
        step = residual/(delta*s+eps*c+gamma*(c*c-s*s))
        theta = np.clip(theta+step, 0, np.pi/2)
        if np.max(np.abs(step)) < 3e-15:
            break
    return eps, np.sin(theta), np.cos(theta)

def tails(n):
    return [zeta(p,n+.5)/(2*np.pi*T)**p for p in (2,3,4)]

def gap_residual(delta, gamma, n=800):
    eps,s,c=spectrum(delta,gamma,n)
    t2,t3,t4=tails(n)
    tail=delta*gamma*t2+(delta**3/2-delta*gamma**2)*t3+(delta*gamma**3-2*delta**3*gamma)*t4
    return delta*np.log(T)+2*np.pi*T*(np.sum(delta/eps-s)+tail)

def gap(Q,n=800):
    if gap_residual(1e-7,Q*Q,n)>0:
        return 0.
    return brentq(lambda delta:gap_residual(delta,Q*Q,n),1e-7,2.2,xtol=3e-13)

def sum_s2(delta,gamma,n=800):
    eps,s,c=spectrum(delta,gamma,n)
    t2,t3,t4=tails(n)
    return np.sum(s*s)+delta**2*t2-2*delta**2*gamma*t3+(3*delta**2*gamma**2-delta**4)*t4

def free(delta,Q,n=800):
    g=Q*Q
    eps,s,c=spectrum(delta,g,n)
    t2,t3,t4=tails(n)
    tail=delta**2*g*t2+(delta**4/4-delta**2*g*g)*t3+(delta**2*g**3-delta**4*g)*t4
    terms=delta**2/eps+2*eps*s*s/(1+c)-2*delta*s+g*s*s
    return delta**2*np.log(T)+2*np.pi*T*(np.sum(terms)+tail)

def delta_al(Q):
    return BCS_AL*np.sqrt(max(0.,1-XI2*(Q/QL)**2/(1-T)))

def current(Q,closure='usadel',n=800):
    d=gap(Q,n) if closure=='usadel' else delta_al(Q)
    if closure=='vodolazov':
        return np.pi*SIGMA/(2*ECH)*ENERGY*d*np.tanh(d/(2*T))*Q/QL*AREA*1e6
    return IPREF*Q*sum_s2(d,Q*Q,n)

def maximum(closure,n=800):
    r=minimize_scalar(lambda Q:-current(Q,closure,n),bounds=(.03,.86),method='bounded',options={'xatol':2e-10})
    d=gap(r.x,n) if closure=='usadel' else delta_al(r.x)
    eg=d*max(0,1-(r.x*r.x/d)**(2/3))**1.5 if d>0 else 0
    return {'I_uA':float(-r.fun),'Q':float(r.x),'delta_meV':float(d*ENERGY/ECH*1000),
            'Eg_meV':float(eg*ENERGY/ECH*1000),'Eg_over_delta':float(eg/d)}

def save(fig,name):
    fig.savefig(FIG/f'{name}.png',bbox_inches='tight')
    fig.savefig(FIG/f'{name}.pdf',bbox_inches='tight')
    plt.close(fig)

def main():
    start=time.perf_counter()
    maxima={k:maximum(k) for k in ('usadel','allmaras_usadel','vodolazov')}
    conv={str(n):maximum('usadel',n) for n in (400,800,1600,3200)}
    points=[]
    for delta,Q in ((.8,.2),(1.1,.4),(1.5,.5)):
        h=1e-4
        fd=(free(delta,Q-2*h)-8*free(delta,Q-h)+8*free(delta,Q+h)-free(delta,Q+2*h))/(12*h)
        expected=4*np.pi*T*Q*sum_s2(delta,Q*Q)
        fd_d=(free(delta-2*h,Q)-8*free(delta-h,Q)+8*free(delta+h,Q)-free(delta+2*h,Q))/(12*h)
        expected_d=2*gap_residual(delta,Q*Q)
        points.append({'delta_kBTc':delta,'Q':Q,'current_derivative_relative_error':float(abs(fd/expected-1)),
                       'amplitude_derivative_relative_error':float(abs(fd_d/expected_d-1))})
    assert max(p['current_derivative_relative_error'] for p in points)<1e-7
    assert max(p['amplitude_derivative_relative_error'] for p in points)<1e-7
    assert max(v['I_uA'] for v in conv.values())-min(v['I_uA'] for v in conv.values())<1e-6
    assert 0<maxima['usadel']['Eg_over_delta']<1
    q=np.linspace(0,.94,160)
    values=np.array([[Q,gap(Q),delta_al(Q),current(Q),current(Q,'allmaras_usadel'),current(Q,'vodolazov')] for Q in q])
    np.savetxt(DATA/'A_ramas_uniformes.csv',values,delimiter=',',header='Q,delta_usadel_kBTc,delta_allmaras_kBTc,I_usadel_uA,I_allmaras_usadel_uA,I_vodolazov_uA',comments='')
    fig,axes=plt.subplots(1,2,figsize=(10,3.8),layout='constrained')
    for i,(label,key) in enumerate((('Usadel autoconsistente','usadel'),('Amplitud Allmaras + Usadel','allmaras_usadel'),('Vodolazov aproximado','vodolazov'))):
        axes[0].plot(q,values[:,3+i],label=label,color=COLORS[i],lw=2)
        axes[0].plot(maxima[key]['Q'],maxima[key]['I_uA'],'o',color=COLORS[i])
    axes[0].set(xlabel=r'$Q=q\sqrt{\hbar D/(2k_BT_c)}$',ylabel=r'$I$ ($\mu$A)',title='Mismo material, distintos cierres',ylim=(0,46))
    axes[0].legend(fontsize=8,loc='upper left')
    axes[1].plot(q,values[:,1],color=COLORS[0],label='Usadel',lw=2)
    axes[1].plot(q,values[:,2],color=COLORS[1],label='Allmaras',lw=2)
    axes[1].set(xlabel=r'$Q$',ylabel=r'$|\Delta|/(k_BT_c)$',title='La amplitud también cambia')
    axes[1].legend(fontsize=9)
    save(fig,'C_04_corriente_uniforme')

    ratio=np.linspace(0,1.2,400)
    eg=np.maximum(0,1-ratio**(2/3))**1.5
    fig,ax=plt.subplots(figsize=(7.8,3.6),layout='constrained')
    ax.plot(ratio,eg,lw=2.5,color=COLORS[0],label=r'$E_g/|\Delta|$')
    ax.axhline(1,color=COLORS[1],ls='--',label=r'Amplitud $|\Delta|$ como referencia')
    ax.fill_between(ratio,eg,1,color=COLORS[0],alpha=.10)
    ax.annotate('Hay estados por debajo de la amplitud',xy=(.4,.31),xytext=(.43,.64),arrowprops={'arrowstyle':'->','color':'#45535d'},fontsize=10)
    ax.set(xlabel=r'$\Gamma/|\Delta|$',ylabel='Energía normalizada',ylim=(-.05,1.13),title='Un gap espectral menor no implica amplitud nula')
    ax.legend(loc='upper right',fontsize=9)
    save(fig,'A_03_gap_espectral')

    # Two actual sections of the reduced uniform functional; no artificial potential.
    fig,axes=plt.subplots(1,2,figsize=(10,3.8),layout='constrained')
    delta_grid=np.linspace(.15,2.15,200)
    for Q,col in ((0,COLORS[0]),(.55,COLORS[1])):
        f=np.array([free(d,Q) for d in delta_grid])
        axes[0].plot(delta_grid,f,lw=2,color=col,label=f'Q = {Q:.2f}')
        deq=gap(Q)
        axes[0].plot(deq,free(deq,Q),'o',color=col)
    axes[0].set(xlabel=r'$|\Delta|/(k_BT_c)$',ylabel=r'$\delta f/[N_0(k_BT_c)^2]$',title='Pendiente en amplitud: fuerza')
    axes[0].legend(fontsize=9)
    q2=np.linspace(0,.65,140); d=1.4
    f=np.array([free(d,Q)-free(d,0) for Q in q2])
    axes[1].plot(q2,f,color=COLORS[0],lw=2)
    Qp=.43; y=free(d,Qp)-free(d,0); slope=4*np.pi*T*Qp*sum_s2(d,Qp*Qp)
    xx=np.linspace(.28,.57,20)
    axes[1].plot(xx,y+slope*(xx-Qp),ls='--',color=COLORS[1],label='Tangente: proporcional a corriente')
    axes[1].plot(Qp,y,'o',color=COLORS[1])
    axes[1].set(xlabel=r'$Q$ (amplitud fija)',ylabel='Incremento de energía libre normalizado',title='Pendiente en flujo: corriente')
    axes[1].legend(fontsize=8)
    save(fig,'A_01_energia_comun')
    report={'scope':'Equilibrio uniforme de acoplamiento débil; sin transitorio ni calibración experimental nueva',
            'parameters':{'T_K':TEMP,'Tc_K':TC,'D_m2_s':D,'sigma_S_m':SIGMA,'w_nm':120,'d_nm':7,'n_matsubara':800},
            'maxima':maxima,'convergence':conv,'derivative_checks':points,'runtime_seconds':time.perf_counter()-start,
            'host':platform.node(),'python':platform.python_version()}
    (DATA/'A_resumen.json').write_text(json.dumps(report,indent=2,ensure_ascii=False),encoding='utf-8')
    print(json.dumps(report,indent=2,ensure_ascii=False))

if __name__=='__main__':main()
