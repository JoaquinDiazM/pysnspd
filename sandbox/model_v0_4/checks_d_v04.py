"""Independent circuit diagnostic for revision 0.4; synthetic imposed resistance.

No detector transient and no production imports. The current is integrated
segment by segment with RK4 and compared with DOP853. Energy integrals use
the trapezoidal rule independently of the RK stages. An exact eta=0 solution
checks each segment without importing a reference trajectory.
"""
from pathlib import Path
import json
import platform
import time
import numpy as np
from scipy.integrate import solve_ivp
import scipy
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT/'docs/modelo_v0_4'
FIG, DATA = OUT/'figuras', OUT/'verificaciones'
for folder in (FIG, DATA):
    folder.mkdir(parents=True, exist_ok=True)
IB, RL, L0, ETA = 20e-6, 50., 100e-9, 1.
TAU0, EUNIT, PUNIT = L0/RL, L0*IB**2, RL*IB**2
BREAKS = np.array([0., 1., 4., 16.])*1e-9
RESISTANCES = np.array([0., 150., 0.])
STEPS_PS = np.array([100., 50., 25., 12.5, 6.25])
BLUE, ORANGE, GREEN, PURPLE = '#236482', '#bf7038', '#32866e', '#8663a0'
plt.rcParams.update({'font.family':'DejaVu Sans', 'font.size':10,
    'axes.spines.top':False, 'savefig.dpi':210})


def rhs(x, r, eta):
    return (1-(1+r)*x)/(1+eta*x*x)


def energy(x, eta):
    """Integral I dPhi/dI dI, divided by L0 Ib^2."""
    return .5*x*x+.25*eta*x**4


def integrate(step_ps, eta):
    rows, segment_rows = [], []
    x, x_exact, x_reference = 1., 1., 1.
    work = np.zeros(3)
    for left, right, resistance in zip(BREAKS[:-1], BREAKS[1:], RESISTANCES):
        r = resistance/RL
        intervals = int(round((right-left)/(step_ps*1e-12)))
        assert intervals >= 4
        times = np.linspace(left, right, intervals+1)
        tau = (times-left)/TAU0
        dtau = tau[1]-tau[0]
        current = np.empty(intervals+1)
        current[0] = x
        for i in range(intervals):
            k1 = rhs(current[i], r, eta)
            k2 = rhs(current[i]+dtau*k1/2, r, eta)
            k3 = rhs(current[i]+dtau*k2/2, r, eta)
            k4 = rhs(current[i]+dtau*k3, r, eta)
            current[i+1] = current[i]+dtau*(k1+2*k2+2*k3+k4)/6
        reference = solve_ivp(lambda _, y: rhs(y, r, eta), (tau[0],tau[-1]),
            [x_reference], method='DOP853', t_eval=tau, rtol=2e-13, atol=2e-15)
        assert reference.success
        current_reference = reference.y[0]
        exact = 1/(1+r)+(x_exact-1/(1+r))*np.exp(-(1+r)*tau)
        x_exact = exact[-1]
        x, x_reference = current[-1], current_reference[-1]
        assert np.all(current > 0) and np.all(current <= 1+1e-12)
        # Normalized powers: Ib Vport, load Joule heat, patch Joule heat.
        powers = np.column_stack((1-current, (1-current)**2, r*current**2))
        integrated = np.vstack((work, work+np.cumsum(
            .5*(powers[:-1]+powers[1:])*dtau, axis=0)))
        work = integrated[-1].copy()
        u = energy(current, eta)
        residual = integrated[:,0]-integrated[:,1]-integrated[:,2]-(u-energy(1.,eta))
        wrong_u = .5*(1+eta*current**2)*current**2
        wrong_residual = integrated[:,0]-integrated[:,1]-integrated[:,2]-(wrong_u-.5*(1+eta))
        # Independent five-point derivatives; do not cross resistance switches.
        dcurrent = (current[:-4]-8*current[1:-3]+8*current[3:-1]-current[4:])/(12*dtau)
        du = (u[:-4]-8*u[1:-3]+8*u[3:-1]-u[4:])/(12*dtau)
        center = current[2:-2]
        kvl = 1-center-r*center-(1+eta*center**2)*dcurrent
        power_residual = powers[2:-2,0]-powers[2:-2,1]-powers[2:-2,2]-du
        segment_rows.append([np.max(abs(kvl)), np.max(abs(power_residual))])
        rows.append(np.column_stack((times,current,np.full_like(times,resistance),
            current_reference,exact,powers,integrated,u,residual,wrong_residual)))
    trajectory = np.concatenate(rows)
    return trajectory, np.max(segment_rows,axis=0)


def save_figure(trajectory, convergence):
    t, x, resistance = trajectory[:,0]*1e9, trajectory[:,1], trajectory[:,2]
    fig, axes = plt.subplots(2,2,figsize=(10.,6.8),layout='constrained')
    ax=axes[0,0]
    ax.plot(t, x*IB*1e6, color=BLUE, label='Corriente')
    ax.set(xlabel='Tiempo (ns)',ylabel='Corriente (µA)',title='Resistencia impuesta: tres tramos')
    ax.set_ylim(0,22)
    ax2=ax.twinx()
    ax2.plot(t,resistance,'--',color=ORANGE,label='Resistencia prescrita')
    ax2.set(ylabel='Resistencia (Ω)',ylim=(0,185))
    ax2.spines['top'].set_visible(False)
    lines1,labels1=ax.get_legend_handles_labels()
    lines2,labels2=ax2.get_legend_handles_labels()
    ax.legend(lines1+lines2,labels1+labels2,fontsize=8,loc='upper right')
    ax=axes[0,1]
    ps,pl,pp=trajectory[:,5:8].T*PUNIT*1e9
    pu=(x+ETA*x**3)*np.array([rhs(xx,rr/RL,ETA) for xx,rr in zip(x,resistance)])*PUNIT*1e9
    for value,color,label in [(ps,BLUE,'Fuente'),(pl+pp,ORANGE,'Calor de ambas resistencias'),(pu,GREEN,'Cambio de energía inductiva')]:
        ax.plot(t,value,color=color,label=label)
    ax.axhline(0,color='grey',lw=.7)
    ax.set(xlabel='Tiempo (ns)',ylabel='Potencia (nW)',title='La energía almacenada también participa')
    ax.legend(fontsize=7.5,loc='lower right')
    ax=axes[1,0]
    ax.plot(t,trajectory[:,12]*EUNIT*1e18,color=BLUE,label='Con energía correcta')
    ax.plot(t,trajectory[:,13]*EUNIT*1e18,color=ORANGE,label=r'Usando erróneamente $L_{\rm diff}I^2/2$')
    ax.axhline(0,color='grey',lw=.7)
    ax.set(xlabel='Tiempo (ns)',ylabel='Defecto del balance integrado (aJ)',title=r'Energía: integrar $I\,d\Phi_{\rm ext}$')
    ax.legend(fontsize=8,loc='center right')
    ax=axes[1,1]
    for column,color,label in [(1,BLUE,'Corriente frente a DOP853'),(2,GREEN,'Ecuación de circuito'),(4,ORANGE,'Energía integrada, trapecios')]:
        ax.loglog(convergence[:,0],np.maximum(convergence[:,column],1e-16),'o-',color=color,label=label)
    ax.set(xlabel='Paso temporal (ps)',ylabel='Error adimensional máximo',title='Refinar separa discretización y fórmula')
    ax.legend(fontsize=7.7)
    for ax in axes.flat:
        if ax is not axes[0,0]:
            ax.spines['right'].set_visible(False)
    for extension in ('png','pdf'):
        fig.savefig(FIG/f'D_01_circuito_y_balance.{extension}',bbox_inches='tight',facecolor='white')
    plt.close(fig)


def main():
    started=time.perf_counter()
    convergence=[]
    for step in STEPS_PS:
        trajectory,residuals=integrate(step,ETA)
        linear,linear_residuals=integrate(step,0.)
        convergence.append([step,np.max(abs(trajectory[:,1]-trajectory[:,3])),
            *residuals,np.max(abs(trajectory[:,12])),np.max(abs(linear[:,1]-linear[:,4])),
            np.max(abs(linear[:,12]))])
    convergence=np.array(convergence)
    orders=np.log2(convergence[:-1,1:]/convergence[1:,1:])
    assert convergence[-1,1] < 2e-9, ('DOP853 current',convergence[-1])
    assert convergence[-1,2] < 1e-8, ('KVL derivative',convergence[-1])
    assert convergence[-1,3] < 1e-8, ('power derivative',convergence[-1])
    assert convergence[-1,4] < 2e-5, ('integrated energy',convergence[-1])
    assert convergence[-1,5] < 2e-10, ('eta zero exact',convergence[-1])
    assert np.all(orders[:3,[0,1,2,4]] > 3.4), ('fourth order',orders)
    assert np.all(orders[:,[3,5]] > 1.85), ('quadrature order',orders)
    # Correct expression suppresses energy defect; the mistaken energy does not.
    wrong=np.max(abs(trajectory[:,13]))
    assert wrong > .15 and wrong > 1e4*np.max(abs(trajectory[:,12]))
    vport=RL*IB*(1-trajectory[:,1])
    ld=L0*(1+ETA*trajectory[:,1]**2)
    output=np.column_stack((trajectory[:,0],trajectory[:,1]*IB,trajectory[:,2],vport,ld,
        trajectory[:,5:8]*PUNIT,trajectory[:,8:11]*EUNIT,trajectory[:,11:14]*EUNIT,
        trajectory[:,3]*IB))
    np.savetxt(DATA/'D04_trayectoria_circuito.csv',output,delimiter=',',comments='',
        header='time_s,current_A,Rpatch_ohm,Vport_V,Ldiff_H,Psource_W,Pload_W,Ppatch_W,Wsource_J,Qload_J,Qpatch_J,Uext_J,balance_defect_J,wrong_energy_balance_defect_J,current_DOP853_A')
    np.savetxt(DATA/'D04_convergencia_circuito.csv',convergence,delimiter=',',comments='',
        header='step_ps,current_DOP853_error_over_Ib,KVL_residual_over_RLIb,power_residual_over_RLIb2,energy_defect_over_L0Ib2,eta0_exact_current_error_over_Ib,eta0_energy_defect_over_L0Ib2')
    save_figure(trajectory,convergence)
    report={'revision':'0.4','diagnostic':'D04','host':platform.node(),
        'python':platform.python_version(),'numpy':np.__version__,'scipy':scipy.__version__,
        'scope':'Circuito con resistencia sintética prescrita; no es un pulso SNSPD predicho, ni una solución de A–C acoplados.',
        'parameters':{'Ib_A':IB,'RL_ohm':RL,'L0_H':L0,'eta':ETA,'initial_current_A':IB,
            'breaks_s':BREAKS.tolist(),'Rpatch_segments_ohm':RESISTANCES.tolist()},
        'equation':'Vport=RL(Ib-I)=Rpatch I+L0[1+eta(I/Ib)^2] dI/dt',
        'energy':'Uext=L0 I^2/2+L0 eta I^4/(4 Ib^2), reference U(0)=0',
        'balance':'Ib Vport=RL(Ib-I)^2+Rpatch I^2+dUext/dt',
        'scales':{'time_s':TAU0,'energy_J':EUNIT,'power_W':PUNIT,'voltage_V':RL*IB},
        'method':'RK4 en cada tramo constante; referencia DOP853; derivadas independientes de 5 puntos sin cruzar saltos; trabajo por trapecios; eta=0 exacto por tramos.',
        'steps_ps':STEPS_PS.tolist(),
        'convergence_columns':['step_ps','current_DOP853','KVL','instantaneous_power','integrated_energy','eta0_exact_current','eta0_integrated_energy'],
        'convergence':convergence.tolist(),'observed_orders':orders.tolist(),
        'finest_max_errors':dict(zip(['current_over_Ib','KVL_over_RLIb','power_over_RLIb2','energy_over_L0Ib2','eta0_current_over_Ib','eta0_energy_over_L0Ib2'],convergence[-1,1:].tolist())),
        'finest_max_energy_defect_J':float(convergence[-1,4]*EUNIT),
        'wrong_energy_max_defect_J':float(wrong*EUNIT),
        'final_energy_ledger_J':{'source_work':float(trajectory[-1,8]*EUNIT),
            'load_heat':float(trajectory[-1,9]*EUNIT),'patch_heat':float(trajectory[-1,10]*EUNIT),
            'change_Uext':float((trajectory[-1,11]-energy(1,ETA))*EUNIT)},
        'current_range_A':[float(np.min(trajectory[:,1])*IB),float(np.max(trajectory[:,1])*IB)],
        'max_port_voltage_V':float(np.max(vport)),'all_checks_passed':True,
        'runtime_seconds':time.perf_counter()-started}
    (DATA/'D04_circuito_y_balance.json').write_text(json.dumps(report,indent=2,ensure_ascii=False)+'\n',encoding='utf-8')
    print(json.dumps(report,indent=2,ensure_ascii=False))


if __name__=='__main__':
    main()
