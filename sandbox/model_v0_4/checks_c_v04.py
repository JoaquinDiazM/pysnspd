"""Independent C v0.4 diagnostics: no production solver or detector transient.

Units: energy k_B Tc, length ell0=sqrt(hbar D/(2 k_B Tc)),
free-energy density N0(k_B Tc)^2; dimensionless time t/tau0.
Thermal uniform Usadel is evaluated by a digamma subtraction valid at large q.
"""
from pathlib import Path
import json, csv, time, platform
import numpy as np
from scipy.special import digamma, polygamma, zeta, expit
from scipy.optimize import brentq, minimize_scalar
from scipy.integrate import solve_ivp, trapezoid
from numpy.polynomial.legendre import leggauss
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'docs/modelo_v0_4'
FIG=OUT/'figuras'; DATA=OUT/'verificaciones'
FIG.mkdir(parents=True,exist_ok=True); DATA.mkdir(parents=True,exist_ok=True)
TC=8.65; TB=.9; TEMP=TB/TC
KB=1.380649e-23; HBAR=1.054571817e-34; D=1.581e-4
ELL0=np.sqrt(HBAR*D/(2*KB*TC)); TAU0=np.pi*HBAR/(8*KB*TC)
K0=np.pi/4
DELTA0=np.pi*np.exp(-np.euler_gamma)
BLUE,RED,GREEN,GOLD='#235789','#c4503c','#278473','#b27a19'
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,
 'axes.spines.top':False,'axes.spines.right':False,'legend.frameon':False,
 'axes.titleweight':'bold','savefig.dpi':200,'grid.alpha':.18})

def save(fig,name):
    fig.savefig(FIG/(name+'.png'),bbox_inches='tight',facecolor='white')
    fig.savefig(FIG/(name+'.pdf'),bbox_inches='tight',facecolor='white')
    plt.close(fig)

def csvout(name,columns,rows):
    with (DATA/name).open('w',encoding='utf-8',newline='') as stream:
        writer=csv.writer(stream); writer.writerow(columns); writer.writerows(rows)

def uniform(a,q,n=600,temp=TEMP):
    """Thermal free energy, amplitude derivative, q derivative; scalars.

    The entire linearized large-Gamma tail is analytic, rather than expanding
    Gamma/epsilon where the core test makes that expansion invalid.
    """
    a=float(a); q=float(q); gamma=q*q
    eps=2*np.pi*temp*(np.arange(n)+.5); den=eps+gamma
    theta=np.arctan(a/den)
    for _ in range(30):
        s=np.sin(theta); c=np.cos(theta)
        step=(a*c-den*s+gamma*(1-c)*s)/(a*s+eps*c+gamma*(c*c-s*s))
        theta=np.clip(theta+step,0,np.pi/2)
        if np.max(np.abs(step))<2e-15:break
    s=np.sin(theta); c=np.cos(theta)
    lin=np.log(temp)+digamma(.5+gamma/(2*np.pi*temp))-digamma(.5)
    tail=lambda p:zeta(p,n+.5+gamma/(2*np.pi*temp))/(2*np.pi*temp)**p
    t3=tail(3); t4=tail(4)
    correction=a*a/den+2*eps*s*s/(1+c)-2*a*s+gamma*s*s
    free=a*a*lin+2*np.pi*temp*(np.sum(correction)+a**4/4*(t3-gamma*t4))
    force=2*a*lin+4*np.pi*temp*(np.sum(a/den-s)+a**3/2*(t3-gamma*t4))
    pi=4*np.pi*temp*q*(np.sum(s*s)+a*a*tail(2)-a**4*(t4-gamma*tail(5)))
    return np.array([free,force,pi])

def array_uniform(a,q,n=600):
    values=np.array([uniform(ai,qi,n) for ai,qi in zip(np.ravel(a),np.ravel(q))])
    return tuple(values[:,j].reshape(np.shape(a)) for j in range(3))

def check_variation():
    """Periodic Cartesian field containing exact zeros, with independent variations."""
    n=128; length=2*np.pi; x=length*np.arange(n)/n; dx=length/n
    modes=2*np.pi*np.fft.fftfreq(n,d=dx)
    deriv=lambda v:np.fft.ifft(1j*modes*np.fft.fft(v))
    psi=DELTA0*(.7*np.cos(x)+.6j*np.sin(2*x))
    psi[[32,96]]=0
    avec=.07*np.cos(x); delta=.1*DELTA0
    def eval_field(z,av,with_force=False):
        cov=lambda v:deriv(v)-1j*av*v
        dz=cov(z); rho=abs(z)**2; a=np.sqrt(rho); denom=rho+delta**2
        pair=np.imag(np.conj(z)*dz); q=pair/denom
        f,X,pi=array_uniform(a,q)
        energy=float(dx*np.sum(f+K0*(abs(dz)**2-rho*q*q)))
        b=(pi-2*K0*rho*q)/denom
        current=2*K0*pair+rho*b
        if not with_force:return energy
        local=np.zeros(n,dtype=complex)
        nz=a>0; local[nz]=X[nz]*z[nz]/a[nz]
        h=local-2*z*(b*q+K0*q*q)-1j*b*dz-cov(2*K0*dz+1j*b*z)
        return energy,h,current
    energy,h,current=eval_field(psi,avec,True)
    direction=.27*np.cos(3*x)+.41j*np.sin(4*x)+.1
    adir=.19*np.cos(2*x)+.11
    step=2e-5
    five=lambda func:(func(-2*step)-8*func(-step)+8*func(step)-func(2*step))/(12*step)
    fd=five(lambda u:eval_field(psi+u*direction,avec))
    pred=float(dx*np.sum(np.real(np.conj(h)*direction)))
    fdA=five(lambda u:eval_field(psi,avec+u*adir))
    predA=-float(dx*np.sum(current*adir))
    chi=.15*np.sin(x)
    gauge=eval_field(psi*np.exp(1j*chi),avec+.15*np.cos(x))
    # Cartesian KWT positive mobility and its analytic inverse.
    tau_sc=1/(TEMP/(.5e-12)+TEMP**3/(2.47e-12))
    beta=4*(KB*TC*tau_sc/HBAR)**2; rho=abs(psi)**2
    rr=np.sqrt(1+beta*rho); abar=np.sqrt((1+TEMP)/2)
    velocity=-rr/(2*abar)*(h-beta*psi*np.real(np.conj(psi)*h)/(1+beta*rho))
    diss=2*abar/rr*(abs(velocity)**2+beta*np.real(np.conj(psi)*velocity)**2)
    power=float(dx*np.sum(np.real(np.conj(h)*velocity)))
    qtot=float(dx*np.sum(diss))
    result={'force_directional_relative_error':abs(fd-pred)/max(1,abs(pred)),
     'current_variation_relative_error':abs(fdA-predA)/max(1,abs(predA)),
     'gauge_energy_relative_error':abs(gauge-energy)/max(1,abs(energy)),
     'kwt_power_relative_error':abs(power+qtot)/max(1,qtot),
     'zero_node_count':2,'zero_node_force_finite':bool(np.isfinite(h[[32,96]]).all()),
     'kwt_min_eigenvalue':float(np.min(2*abar/rr)),
     'discretization':'128 periodic collocation nodes; independent five-point directional differences'}
    assert max(result[k] for k in ('force_directional_relative_error','current_variation_relative_error','gauge_energy_relative_error','kwt_power_relative_error'))<2e-6,result
    return result

def core_diagnostics():
    aeq=brentq(lambda a:uniform(a,0)[1],.5,2.1,xtol=1e-12)
    radius=np.geomspace(1e-5,8,1400)
    a=aeq*np.tanh(radius); ap=aeq/np.cosh(radius)**2
    app=-2*aeq*np.tanh(radius)/np.cosh(radius)**2
    lapcomplex=app+ap/radius-a/radius**2
    small=radius<1e-3
    lapcomplex[small]=aeq*(-8*radius[small]/3+16*radius[small]**3/5)
    rho=a*a; pair=rho/radius; g=ap*ap+rho/radius**2
    fbulk=uniform(aeq,0)[0]
    rows=[]; profiles={}; energies={}; norms={}
    for fraction in (.05,.1,.2):
        delta=fraction*DELTA0; denom=rho+delta**2; q=pair/denom; m=rho/denom
        f,X,pi=array_uniform(a,q)
        dq_da=2*a*delta**2/(denom**2*radius)
        force=X-2*K0*lapcomplex-2*K0*a*q*q+(pi-2*K0*rho*q)*dq_da
        density=f-fbulk+K0*(g-rho*q*q)
        current=m*pi+2*K0*(1-m*m)*pair
        energy=float(trapezoid(2*np.pi*radius*density,radius))
        energies[str(fraction)]=energy
        norms[str(fraction)]={'max_abs_force':float(np.max(abs(force))),
          'force_L2':float(np.sqrt(trapezoid(2*np.pi*radius*force**2,radius))),
          'radius_max_force_ell0':float(radius[np.argmax(abs(force))])}
        profiles[fraction]=(force,density,current,q)
        rows.extend(zip(np.full(len(radius),fraction),radius,a,q,force,density,current))
    csvout('C04_vortex_profiles.csv',['delta_over_Delta0','r_over_ell0','a_over_kBTc','q_delta_ell0','force_dimensionless','energy_excess_density','current_moment_dimensionless'],rows)
    # Exact original local polar closure on the same smooth vortex profile.
    f_orig,X_orig,_=array_uniform(a,1/radius)
    original_force=X_orig-2*K0*(app+ap/radius)
    original_energy=float(trapezoid(2*np.pi*radius*(f_orig-fbulk+K0*ap*ap),radius))
    # Forced regulator-to-zero diagnostic is restricted to the analytic core.
    tiny=np.geomspace(1e-7,.2,1600); aa=tiny
    limit=[]
    for delta in (.04,.02,.01,.005,.0025,.00125):
        ss=aa*aa+delta*delta; qq=aa*aa/(ss*tiny)
        ll=np.log(TEMP)+digamma(.5+qq*qq/(2*np.pi*TEMP))-digamma(.5)
        xx=2*aa*ll
        pp=aa*aa*polygamma(1,.5+qq*qq/(2*np.pi*TEMP))*qq/(np.pi*TEMP)
        qq_a=2*aa*delta*delta/(ss*ss*tiny)
        hh=xx-2*K0*aa*qq*qq+(pp-2*K0*aa*aa*qq)*qq_a
        density=aa*aa*ll+K0*(2-aa*aa*qq*qq)
        limit.append((delta,float(np.max(abs(hh))),float(trapezoid(2*np.pi*tiny*density,tiny))))
    csvout('C04_core_delta_limit.csv',['delta_over_kBTc','max_abs_force','energy_on_r_le_0p2'],limit)
    fig,ax=plt.subplots(1,3,figsize=(11.8,3.7),layout='constrained')
    for frac,col in zip((.05,.1,.2),(BLUE,RED,GREEN)):
        ax[0].plot(radius,profiles[frac][0],color=col,label=rf'$\delta/\Delta_0={frac:g}$')
        ax[1].plot(radius,profiles[frac][2],color=col)
    ax[0].plot(radius,original_force,'--',color=GOLD,label='Polar original')
    ax[0].set(xscale='log',xlim=(.002,2),ylim=(-65,15),xlabel=r'$r/\ell_0$',ylabel='Fuerza radial normalizada',title='Un núcleo revela la sensibilidad')
    ax[0].legend(fontsize=8)
    ax[1].set(xscale='log',xlim=(.002,8),xlabel=r'$r/\ell_0$',ylabel='Momento de corriente',title='La corriente también se deriva')
    lim=np.array(limit)
    ax[2].loglog(lim[:,0],lim[:,1],'o-',color=RED,label='Máxima fuerza')
    ax[2].loglog(lim[:,0],lim[0,1]*lim[0,0]/lim[:,0],'--',color=GOLD,label=r'Guía $1/\delta$')
    ax[2].set(xlabel=r'$\delta/(k_BT_c)$',ylabel='Máxima fuerza normalizada',title=r'No se toma $\delta\to0$')
    ax[2].legend(fontsize=8)
    save(fig,'C_05_auditoria_nucleo')
    # Independent sum truncation check, including strong depairing.
    conv=[]
    for av,qv in ((.1,20),(1.3,.4),(aeq,0)):
        coarse=uniform(av,qv,600); fine=uniform(av,qv,1200)
        conv.append({'a':av,'q':qv,'max_scaled_difference':float(np.max(abs(coarse-fine)/np.maximum(1,abs(fine))))})
    return {'profile':'Prescribed Delta=aeq*tanh(r/ell0)*exp(i theta), disk r<=8ell0; not a stationary vortex',
      'aeq_over_kBTc':aeq,'ell0_nm':ELL0*1e9,'regularized_energy_dimensionless':energies,
      'original_energy_dimensionless':original_energy,'force_norms':norms,
      'energy_spread_over_reference':(max(energies.values())-min(energies.values()))/abs(energies['0.1']),
      'max_force_spread_over_reference':(max(v['max_abs_force'] for v in norms.values())-min(v['max_abs_force'] for v in norms.values()))/norms['0.1']['max_abs_force'],
      'delta_to_zero_core_test':limit,'matsubara_convergence':conv,
      'decision':'delta=0.1Delta0 is the fixed reference of an experimental core closure, not a calibrated material parameter; no production promotion'}

def finite_modes():
    """Independent q=0 spatial Usadel linear response at equilibrium gap."""
    a=brentq(lambda av:uniform(av,0)[1],.5,2.1)
    wave=np.linspace(0,3,101)
    def stiffness(n):
        eps=2*np.pi*TEMP*(np.arange(n)+.5); en=np.hypot(eps,a); cc=eps/en
        h0=(uniform(a+1e-4,0,n)[1]-uniform(a-1e-4,0,n)[1])/2e-4
        response=np.array([h0+4*np.pi*TEMP*np.sum(cc*cc*kv*kv/(en*(en+kv*kv))) for kv in wave])
        # Leading omitted tail for k^2/epsilon^2.
        response+=4*np.pi*TEMP*wave*wave*zeta(2,n+.5)/(2*np.pi*TEMP)**2
        return h0,response
    h0,full=stiffness(2500); _,refined=stiffness(5000)
    local=h0+2*K0*wave*wave
    tau_sc=1/(TEMP/(.5e-12)+TEMP**3/(2.47e-12))
    rr=np.sqrt(1+4*(a*KB*TC*tau_sc/HBAR)**2); gamma=2*np.sqrt((1+TEMP)/2)*rr
    tau_full=gamma/full*TAU0*1e12; tau_local=gamma/local*TAU0*1e12
    csvout('C04_finite_wavevector.csv',['k_ell0','H_spatial_Usadel','H_local_K0','tau_spatial_with_same_KWT_ps','tau_local_with_same_KWT_ps'],zip(wave,full,local,tau_full,tau_local))
    fig,ax=plt.subplots(1,2,figsize=(9.2,3.7),layout='constrained')
    ax[0].plot(wave,full,color=BLUE,label='Usadel espacial linealizado')
    ax[0].plot(wave,local,color=RED,label='Gradiente local K0')
    ax[0].set(xlabel=r'$k\ell_0$',ylabel='Rigidez de amplitud normalizada',title='La rigidez depende de la longitud')
    ax[0].legend(fontsize=8)
    ax[1].plot(wave,tau_full,color=BLUE);ax[1].plot(wave,tau_local,color=RED)
    ax[1].set(xlabel=r'$k\ell_0$',ylabel='Tiempo lineal condicional (ps)',title='Misma movilidad: distinta relajación')
    save(fig,'C_06_modos_espaciales')
    probes=[]
    for kval in (.3,1.,2.):
        i=np.argmin(abs(wave-kval)); probes.append({'k_ell0':float(wave[i]),'relative_stiffness_excess':float(local[i]/full[i]-1),'tau_full_ps':float(tau_full[i]),'tau_local_ps':float(tau_local[i])})
    return {'scope':'Isothermal linear amplitude perturbations around uniform q=0 BCS equilibrium; identical assumed KWT mobility, no kinetic feedback',
      'H0':h0,'KWT_R':rr,'tau0_ps':TAU0*1e12,'probes':probes,
      'sum_refinement_max_relative_difference':float(np.max(abs(full-refined)/refined))}

def phase_principal_diagnostic():
    """Exact vacuum constitutive derivative; detects loss of phase parabolicity.

    In the gapped vacuum regime Gamma=q_delta^2<a, integration of the
    zero-temperature Matsubara equation gives Pi=pi*a*q_delta-4*q_delta^3/3.
    Derivatives here hold amplitude and p=0 fixed, not the equilibrium gap.
    """
    def vacuum_pi_quadrature(a,q,n=192):
        points,weights=leggauss(n); y=(points+1)/2
        eps=a*y/(1-y); ww=weights*a/(2*(1-y)**2); gamma=q*q
        theta=np.arctan2(a,eps)
        for _ in range(40):
            ss=np.sin(theta); cc=np.cos(theta)
            step=(a*cc-eps*ss-gamma*ss*cc)/(a*ss+eps*cc+gamma*(cc*cc-ss*ss))
            theta=np.clip(theta+step,0,np.pi/2)
            if np.max(abs(step))<3e-15:break
        return float(2*q*np.dot(ww,np.sin(theta)**2))
    rows=[]
    for ratio,q in ((.6,.65),(.6,1.),(1.,1.)):
        a=ratio*DELTA0; delta=.1*DELTA0
        m=a*a/(a*a+delta*delta); qeff=m*q
        assert qeff*qeff<a
        phase=m*m*(np.pi*a-4*qeff*qeff)+2*K0*a*a*(1-m*m)
        current=lambda qq:m*vacuum_pi_quadrature(a,m*qq)+2*K0*a*a*qq*(1-m*m)
        step=1e-4
        numerical=(current(q-2*step)-8*current(q-step)+8*current(q+step)-current(q+2*step))/(12*step)
        rows.append({'amplitude_over_Delta0':ratio,'q_ell0':q,'q_delta_ell0':qeff,
          'phase_principal_coefficient':phase,'independent_Matsubara_quadrature_derivative_absolute_error':abs(phase-numerical),
          'quadrature_refinement_absolute_difference':abs(vacuum_pi_quadrature(a,qeff,192)-vacuum_pi_quadrature(a,qeff,384))})
    assert rows[0]['phase_principal_coefficient']>0
    assert rows[1]['phase_principal_coefficient']<0
    wave=np.linspace(.4,1.03,200); a=.6*DELTA0
    fig,ax=plt.subplots(figsize=(7.4,3.7),layout='constrained')
    data=[]
    for fraction,col in zip((.05,.1,.2),(BLUE,RED,GREEN)):
        delta=fraction*DELTA0; m=a*a/(a*a+delta*delta)
        assert np.max((m*wave)**2)<a
        principal=m*m*(np.pi*a-4*(m*wave)**2)+2*K0*a*a*(1-m*m)
        ax.plot(wave,principal,color=col,label=rf'$\delta/\Delta_0={fraction:g}$')
        data.extend(zip(np.full(len(wave),fraction),wave,principal))
    ax.axhline(0,color='#444444',linewidth=1)
    ax.axhspan(-1,0,color='#c4503c',alpha=.09)
    ax.plot([1],[rows[1]['phase_principal_coefficient']],'o',color=RED)
    ax.annotate(r'$H_\theta=-0.3434$',(1,rows[1]['phase_principal_coefficient']),
      xytext=(.65,-.61),arrowprops={'arrowstyle':'->','color':RED},color=RED)
    ax.set(xlabel=r'$q\ell_0$',ylabel=r'Rigidez de fase $H_\theta$ a amplitud fija',
      title=r'Una movilidad positiva no corrige una rigidez negativa',ylim=(-.85,2.8))
    ax.legend(fontsize=9)
    save(fig,'C_09_admisibilidad_espacial')
    csvout('C04_phase_principal.csv',['delta_over_Delta0','q_ell0','H_phase_vacuum_at_fixed_a_0p6Delta0'],data)
    return {'scope':'Vacuum p=0, fixed amplitude. Exact gapped zero-temperature constitutive derivative, not derivative along a self-consistent branch.',
      'probes':rows,
      'decision':'The local finite-core PDE is not globally parabolic. A positive principal symbol is an admission condition; loss requires stopping the candidate, not interpreting ultraviolet growth as detector physics.'}

def uniform_regularized():
    """Uniform branches of the selected finite-core functional, not old j alone."""
    sigma=4.2e5; area=120e-9*7e-9; charge=1.602176634e-19
    ipref=sigma*(KB*TC)/(2*charge*ELL0)*area*1e6
    def values(a,q,fraction):
        if fraction==0:return uniform(a,q)
        delta=fraction*DELTA0; rho=a*a; den=rho+delta*delta
        m=rho/den; m_a=2*a*delta*delta/(den*den)
        ff,xx,pp=uniform(a,m*q)
        free=ff+K0*rho*q*q*(1-m*m)
        force=xx+pp*q*m_a+2*K0*a*q*q*(1-m*m)-2*K0*rho*m*q*q*m_a
        current=m*pp+2*K0*rho*q*(1-m*m)
        return np.array([free,force,current])
    multiple=0
    def state(q,fraction):
        nonlocal multiple
        trials=np.r_[0.,np.geomspace(1e-7,.15,12),np.linspace(.17,2.1,30)]
        forces=[values(aa,q,fraction)[1] for aa in trials[1:]]
        roots=[]
        for j in range(len(forces)-1):
            if forces[j]*forces[j+1]<0:
                roots.append(brentq(lambda aa:values(aa,q,fraction)[1],trials[j+1],trials[j+2],xtol=1e-11))
        multiple+=int(len(roots)>1)
        candidates=[(0.,0.)]+[(values(aa,q,fraction)[0],aa) for aa in roots]
        a=min(candidates)[1]
        return a,values(a,q,fraction)[2]*ipref
    qs=np.linspace(0,1.8,121)
    allrows=[];report={};curves={}
    for fraction in (0.,.05,.1,.2):
        vals=np.array([state(q,fraction) for q in qs]);curves[fraction]=vals
        peak=minimize_scalar(lambda q:-state(q,fraction)[1],bounds=(.3,.85),method='bounded',options={'xatol':2e-8})
        pa,pi=state(peak.x,fraction)
        report[str(fraction)]={'maximum_current_uA':float(pi),'q_peak_ell0':float(peak.x),'gap_at_peak_over_Delta0':float(pa/DELTA0),'gap_at_q1_over_Delta0':float(state(1.,fraction)[0]/DELTA0)}
        allrows.extend((fraction,q,*vv) for q,vv in zip(qs,vals))
    for key,item in report.items():item['current_max_change_percent']=100*(item['maximum_current_uA']/report['0.0']['maximum_current_uA']-1)
    csvout('C04_uniform_regularized.csv',['delta_over_Delta0','q_ell0','gap_over_kBTc','current_uA'],allrows)
    fig,ax=plt.subplots(1,2,figsize=(9.4,3.7),layout='constrained')
    for fraction,col in zip((0.,.05,.1,.2),(GOLD,BLUE,RED,GREEN)):
        label='Usadel uniforme' if fraction==0 else rf'$\delta/\Delta_0={fraction:g}$'
        ax[0].plot(qs,curves[fraction][:,1],color=col,label=label)
        ax[1].plot(qs,curves[fraction][:,0]/DELTA0,color=col)
    ax[0].set(xlabel=r'$q\ell_0$',ylabel=r'$I$ ($\mu$A)',title='Rama del mismo funcional regularizado')
    ax[0].legend(fontsize=8)
    ax[1].set(xlabel=r'$q\ell_0$',ylabel=r'$|\Delta|/\Delta_0$',title='Una corriente próxima no valida el núcleo')
    ax[1].axvline(1.,color='#777777',lw=.7,ls=':')
    save(fig,'C_08_ramas_regularizadas')
    return {'reference_material':{'sigma_S_m':sigma,'w_nm':120,'d_nm':7},'branches':report,
     'multiple_stationary_root_search_events':multiple,
     'selection':'lowest free energy among positive stationary roots and normal state at fixed q; 42-point bracketing grid; maximum on main branch',
     'scope':'Uniform equilibrium at Tb, using current and force from the same regularized functional. Not device switching.'}

def isolated_cell(order=96):
    """Nonthermal BCS cell, selected effective Q_Delta source, finite KWT."""
    nodes,weights=leggauss(order); x=(nodes+1)*10; w=weights*10
    ainit=.6*DELTA0
    pinit=expit(-np.hypot(x,ainit)/TEMP)
    initial=np.r_[ainit,pinit]
    def vacuum(a):return a*a*(np.log(a/DELTA0)-.5)
    def quantities(y):
        a=y[0];p=y[1:];en=np.hypot(x,a)
        uq=4*np.dot(w,en*p)
        temperature=brentq(lambda tt:4*np.dot(w,en*expit(-en/tt))-uq,1e-6,5,xtol=1e-12)
        force=2*a*np.log(a/DELTA0)+4*np.dot(w,a/en*p)
        mobility_temperature=max(temperature,TEMP)
        tau_sc=1/(mobility_temperature/(.5e-12)+mobility_temperature**3/(2.47e-12))
        rr=np.sqrt(1+4*(a*KB*TC*tau_sc/HBAR)**2)
        gamma=2*np.sqrt((1+mobility_temperature)/2)*rr
        adot=-force/gamma; heat=force*force/gamma
        seed=expit(-en/max(temperature,TEMP)); source=en*seed*(1-p)
        pdot=heat*source/(4*np.dot(w,en*source))
        return adot,pdot,temperature,vacuum(a)+uq,force,heat,uq
    def rhs(t,y):
        adot,pdot,*_=quantities(y);return np.r_[adot,pdot]
    times=np.linspace(0,250,351)
    sol=solve_ivp(rhs,(times[0],times[-1]),initial,t_eval=times,method='DOP853',rtol=3e-9,atol=2e-12)
    assert sol.success,sol.message
    outputs=np.array([[y[0],*quantities(y)[2:]] for y in sol.y.T])
    return times,sol.y,outputs

def cell_diagnostics():
    tt,state,vals=isolated_cell(96)
    _,ref_state,ref_vals=isolated_cell(192)
    energy=vals[:,2];drift=float(np.max(abs(energy-energy[0])))
    csvout('C04_isolated_BCS_cell.csv',['t_over_tau0','t_ps','a_over_kBTc','TE_over_Tc','energy_total','force','QDelta','u_qp'],((t,t*TAU0*1e12,*v) for t,v in zip(tt,vals)))
    fig,ax=plt.subplots(1,3,figsize=(11.6,3.5),layout='constrained')
    ax[0].plot(tt*TAU0*1e12,vals[:,0]/DELTA0,color=BLUE,label=r'$|\Delta|/\Delta_0$')
    ax[0].plot(tt*TAU0*1e12,vals[:,1],color=RED,label=r'$T_E/T_c$')
    ax[0].set(xlabel='Tiempo condicional (ps)',ylabel='Magnitud normalizada',title='Recuperación en una celda BCS')
    ax[0].legend(fontsize=8)
    ax[1].plot(tt*TAU0*1e12,energy-vals[:,5],color=BLUE,label='Fondo emparejado')
    ax[1].plot(tt*TAU0*1e12,vals[:,5],color=RED,label='Excitaciones')
    ax[1].plot(tt*TAU0*1e12,energy,color=GREEN,label='Suma')
    ax[1].set(xlabel='Tiempo condicional (ps)',ylabel='Energía normalizada',title='La disipación transfiere energía')
    ax[1].legend(fontsize=8)
    ax[2].plot(tt*TAU0*1e12,energy-energy[0],color=GREEN)
    ax[2].set(xlabel='Tiempo condicional (ps)',ylabel='Cambio de la energía total',title='Residuo de integración')
    ax[2].ticklabel_format(axis='y',style='sci',scilimits=(0,0))
    save(fig,'C_07_celda_BCS_cinetica')
    assert drift<2e-7,drift
    return {'scope':'Isolated q=0 BCS cell: evolving spectral p, no e-ph/ee/escape/Joule; effective QDelta source from B and specified KWT times. A subsystem diagnostic, not a detector transient.',
     'initial_a_over_Delta0':.6,'initial_TE_K':TB,'duration_ps':float(tt[-1]*TAU0*1e12),
     'final_a_over_Delta0':float(vals[-1,0]/DELTA0),'final_TE_K':float(vals[-1,1]*TC),
     'energy_drift_absolute_normalized':drift,'refinement_max_gap_difference':float(np.max(abs(vals[:,0]-ref_vals[:,0]))),
     'refinement_max_TE_difference':float(np.max(abs(vals[:,1]-ref_vals[:,1]))),
     'population_range':[float(np.min(state[1:])),float(np.max(state[1:]))],
     'solver':'DOP853 rtol3e-9 atol2e-12, 96 and192 Gauss nodes on x/(kBTc) in[0,20]'}

def main():
    start=time.perf_counter()
    report={'version':'0.4','host':platform.node(),'python':platform.python_version(),
      'parameters':{'Tc_K':TC,'Tb_K':TB,'D_m2_s':D,'delta_reference_over_Delta0':.1,
      'tau_ee_KWT_Tc_ps':.5,'tau_ep_KWT_Tc_ps':2.47},
      'variation':check_variation(),'core':core_diagnostics(),'finite_modes':finite_modes(),'cell':cell_diagnostics(),'uniform_regularized':uniform_regularized(),
      'phase_principal':phase_principal_diagnostic()}
    report['runtime_seconds']=time.perf_counter()-start
    (DATA/'C04_checks.json').write_text(json.dumps(report,indent=2,ensure_ascii=False)+'\n',encoding='utf-8')
    print(json.dumps(report,indent=2,ensure_ascii=False))

if __name__=='__main__':main()
