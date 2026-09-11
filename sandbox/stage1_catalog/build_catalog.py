"""Build and independently diagnose a small uniform catalogue (<5 min expected).

Run with OMP_NUM_THREADS=OPENBLAS_NUM_THREADS=MKL_NUM_THREADS=1.
No production solver, mesh, transient, or material-collision catalogue is used.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import platform
import sys
import time

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT))
import numpy as np
import scipy
from scipy.integrate import quad
from scipy.optimize import brentq
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from pysnspd.experimental.energy_catalog import (
    BCS_GAP_RATIO, E_CHARGE_C, build_vacuum_catalog, build_occupation_catalog,
    retarded_spectrum, vacuum_state, energy_at_count, spectral_count,
)


def derivative(fun,x,h=1e-4):
    return (fun(x-2*h)-8*fun(x-h)+8*fun(x+h)-fun(x+2*h))/(12*h)


def matsubara_reference(T,delta,gamma,n):
    """Independent bracketing in spectral angle, not the production formula."""
    eps=np.pi*T*(2*np.arange(n)+1)
    theta=np.array([brentq(lambda th:e*np.sin(th)-delta*np.cos(th)+gamma*np.sin(th)*np.cos(th),
                          0,np.pi/2,xtol=1e-15) for e in eps])
    s,c=np.sin(theta),np.cos(theta)
    # 1-c=2sin(theta/2)^2 avoids cancellation of the normal reference.
    terms=delta*delta/eps+4*eps*np.sin(theta/2)**2-2*delta*s+gamma*s*s
    free=-np.pi**2*T*T/3+delta*delta*np.log(T*BCS_GAP_RATIO)+2*np.pi*T*np.sum(terms)
    force=2*delta*np.log(T*BCS_GAP_RATIO)+4*np.pi*T*np.sum(delta/eps-s)
    conjugate=2*np.pi*T*np.sum(s*s)
    return np.array([free,force,conjugate])


def matsubara_converged(T,delta,gamma):
    cutoffs=np.array([256,512,1024,2048])
    values=np.array([matsubara_reference(T,delta,gamma,int(n)) for n in cutoffs])
    extrap=np.polynomial.polynomial.polyfit(cutoffs[0]/cutoffs,values,3)[0]
    return extrap,values


def direct_fixed_p(delta,gamma,nodes,weights,p,eta):
    energy=energy_at_count(nodes,delta=delta,gamma=gamma,eta=eta)
    c,s=retarded_spectrum(energy,delta=delta,gamma=gamma,eta=eta)
    kernels=np.array([energy,s.imag/c.real,-s.real*s.imag/c.real])
    return np.asarray(vacuum_state(delta,gamma))+4*np.sum(kernels*(weights*p)[None,:],axis=1)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,default=ROOT/'docs/implementation/stage1/electronic')
    parser.add_argument('--catalogs',type=Path,default=ROOT/'docs/implementation/stage1/catalogs')
    args=parser.parse_args()
    out=args.output
    figures=out/'figures'
    figures.mkdir(parents=True,exist_ok=True)
    started=time.perf_counter()
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,
                         'savefig.dpi':180,'font.family':'DejaVu Sans'})
    sources={p.relative_to(ROOT).as_posix():hashlib.sha256(p.read_bytes()).hexdigest()
             for p in sorted((ROOT/'docs/modelo_v0_4').glob('[ABD]_*.md'))}
    D,sigma=1.581e-4,4.2e5
    N0=sigma/(2*E_CHARGE_C**2*D)
    kwargs=dict(Tc_K=8.65,N0_per_J_m3=N0,D_m2_s=D,reference_hashes=sources)
    points=np.random.default_rng(41003).uniform([.09,.008],[1.48,1.16],size=(80,2))
    grid_rows=[]
    for n in (17,33,65):
        table=build_vacuum_catalog(np.linspace(.04,1.6,n),np.linspace(0,1.2,n),**kwargs)
        errors=np.array([np.abs(np.array(table.evaluate(*p))-vacuum_state(*p)) for p in points])
        grid_rows.append([n,*np.max(errors,axis=0)])
    table.save(args.catalogs/'vacuum_catalog.npz')
    grid_rows=np.asarray(grid_rows)
    np.savetxt(out/'vacuum_grid_refinement.csv',grid_rows,delimiter=',',
               header='field_nodes,max_energy_error,max_amplitude_error,max_gamma_error',comments='')
    identity=[]
    for delta,gamma in points:
        _,force,conjugate=table.evaluate(delta,gamma)
        identity.append([
            abs(derivative(lambda d:table.evaluate(d,gamma)[0],delta)-force),
            abs(derivative(lambda g:table.evaluate(delta,g)[0],gamma)-conjugate),
            abs(derivative(lambda d:table.evaluate(d,gamma)[2],delta)-
                derivative(lambda g:table.evaluate(delta,g)[1],gamma))])

    # A genuinely non-FD occupation: narrow band in the state-count coordinate.
    # Reference derivatives use complex spectral factors, not spline derivatives.
    nonthermal_points=np.array([[.17,.012],[.37,.14],[.72,.3],[1.17,.72],[.43,.8],[1.2,.002]])
    nonthermal_rows=[]
    low_energy_rows=[]
    low_energy_points=np.array([[.17,0],[.72,0],[1.2,.0001],[.72,.002],[.43,.02]])
    nonthermal_identity=[]
    for n in (9,17,33):
        grid_started=time.perf_counter()
        n_low=n//3
        gamma_axis=np.r_[0,np.geomspace(1e-4,.03,n_low),np.linspace(.03,1.2,n-n_low)[1:]]
        vacuum=build_vacuum_catalog(np.linspace(.08,1.5,n),gamma_axis,**kwargs)
        occupation_table=build_occupation_catalog(vacuum,count_order=64,count_max=6,eta=.001)
        p=.16*np.exp(-((occupation_table.count_nodes-.95)/.3)**2)
        errors=[]
        for delta,gamma in nonthermal_points:
            value=np.array(occupation_table.evaluate(delta,gamma,p))
            reference=direct_fixed_p(delta,gamma,occupation_table.count_nodes,occupation_table.count_weights,p,.001)
            errors.append(np.abs(value-reference))
        nonthermal_rows.append([n,*np.max(errors,axis=0)])
        p_low=.2*np.exp(-occupation_table.count_nodes/.25)
        low_errors=[]
        for delta,gamma in low_energy_points:
            value=np.array(occupation_table.evaluate(delta,gamma,p_low))
            reference=direct_fixed_p(delta,gamma,occupation_table.count_nodes,occupation_table.count_weights,p_low,.001)
            low_errors.append(np.abs(value-reference))
        low_energy_rows.append([n,*np.max(low_errors,axis=0)])
        print(f'Occupation catalogue {n}x{n}x64 completed in {time.perf_counter()-grid_started:.2f}s',flush=True)
    occupation_table.save(args.catalogs/'occupation_catalog.npz')
    # Keep p constant during variation, including when its origin is nonthermal.
    for delta,gamma in nonthermal_points:
        _,force,conjugate=occupation_table.evaluate(delta,gamma,p)
        h=min(1e-5,gamma/5)
        nonthermal_identity.append([
            abs(derivative(lambda d:occupation_table.evaluate(d,gamma,p)[0],delta,h)-force),
            abs(derivative(lambda g:occupation_table.evaluate(delta,g,p)[0],gamma,h)-conjugate),
            abs(derivative(lambda d:occupation_table.evaluate(d,gamma,p)[2],delta,h)-
                derivative(lambda g:occupation_table.evaluate(delta,g,p)[1],gamma,h))])
    nonthermal_rows=np.asarray(nonthermal_rows)
    np.savetxt(out/'nonthermal_grid_refinement.csv',nonthermal_rows,delimiter=',',
               header='field_nodes,max_energy_error,max_amplitude_error,max_gamma_error',comments='')
    np.savetxt(out/'low_energy_grid_refinement.csv',low_energy_rows,delimiter=',',
               header='field_nodes,max_energy_error,max_amplitude_error,max_gamma_error',comments='')
    # Independent count-quadrature/cutoff checks use no interpolated field values.
    quadrature=[]
    low_quadrature=[]
    for support in (4,6,8):
        for order in (32,64,128):
            nodes,weights=np.polynomial.legendre.leggauss(order)
            nodes=(nodes+1)*support/2
            weights=weights*support/2
            p_band=.16*np.exp(-((nodes-.95)/.3)**2)
            result=direct_fixed_p(.72,.3,nodes,weights,p_band,.001)
            quadrature.append([support,order,*result])
            p_low=.2*np.exp(-nodes/.25)
            low_result=direct_fixed_p(.72,0,nodes,weights,p_low,.001)
            low_quadrature.append([support,order,*low_result])
    np.savetxt(out/'nonthermal_quadrature.csv',quadrature,delimiter=',',
               header='count_cutoff,count_order,energy,amplitude_force,gamma_conjugate',comments='')
    np.savetxt(out/'low_energy_quadrature.csv',low_quadrature,delimiter=',',
               header='count_cutoff,count_order,energy,amplitude_force,gamma_conjugate',comments='')

    # Auxiliary-temperature identity, compared against independently solved Matsubara sums.
    thermal_rows=[]
    thermal_force_rows=[]
    matsubara_rows=[]
    for delta,gamma in ((1,0),(.7,.2),(.35,.65)):
        for T in (.1,.25,.5):
            reference,raw=matsubara_converged(T,delta,gamma)
            for n,row in zip((256,512,1024,2048),raw):
                matsubara_rows.append([delta,gamma,T,n,*np.abs(row-reference)])
            for eta in (.001,.0001):
                # Integrate energy directly with the analytic DOS and adaptive quadrature.
                # Split around the ideal spectral edge; eta only supplies the regulator.
                edge=delta*max(0,1-(gamma/delta)**(2/3))**1.5
                def thermal_integrand(E):
                    c,_=retarded_spectrum(E,delta=delta,gamma=gamma,eta=eta)
                    return float(c.real*np.logaddexp(0,-E/T))
                breaks=sorted(set([0,edge,delta,12*T+2*delta]))+[np.inf]
                entropy_correction=4*T*sum(quad(thermal_integrand,a,b,epsabs=2e-10,
                                              epsrel=2e-9,limit=150)[0]
                                             for a,b in zip(breaks[:-1],breaks[1:]))
                thermal_rows.append([delta,gamma,T,eta,reference[0]+entropy_correction,
                                     vacuum_state(delta,gamma)[0],
                                     abs(reference[0]+entropy_correction-vacuum_state(delta,gamma)[0])])
                if T==.25:
                    for order in (64,128,256):
                        x,w=np.polynomial.legendre.leggauss(order)
                        x=(x+1)*4
                        w=w*4
                        energies_fd=energy_at_count(x,delta=delta,gamma=gamma,eta=eta)
                        p_fd=1/(np.exp(energies_fd/T)+1)
                        direct=direct_fixed_p(delta,gamma,x,w,p_fd,eta)
                        free=vacuum_state(delta,gamma)[0]-4*T*np.dot(np.logaddexp(0,-energies_fd/T),w)
                        thermal_force_rows.append([delta,gamma,T,eta,order,free-reference[0],
                                                   direct[1]-reference[1],direct[2]-reference[2]])
    np.savetxt(out/'vacuum_auxiliary_temperature.csv',thermal_rows,delimiter=',',
               header='delta,gamma,T_over_Delta0,eta,reconstructed_vacuum,exact_vacuum,error',comments='')
    np.savetxt(out/'matsubara_cutoffs.csv',matsubara_rows,delimiter=',',
               header='delta,gamma,T_over_Delta0,cutoff,free_error,amplitude_error,gamma_error',comments='')
    np.savetxt(out/'thermal_fixed_p_reference.csv',thermal_force_rows,delimiter=',',
               header='delta,gamma,T_over_Delta0,eta,count_order,free_difference,amplitude_difference,gamma_difference',comments='')

    # Preserve all complex factors for a compact spectral pilot at independent fields.
    energies=np.unique(np.r_[np.linspace(0,6,301),.1,.4,.8,1.2,1.6])
    deltas=np.array([0,.1,.4,.8,1.2,1.6])
    gammas=np.array([0,.05,.2,.5,.8,1.2])
    c_grid,s_grid=[],[]
    residuals=[]
    for delta in deltas:
        c_row,s_row=[],[]
        for gamma in gammas:
            c,s=retarded_spectrum(energies,delta=delta,gamma=gamma,eta=.001)
            c_row.append(c);s_row.append(s)
            lhs=delta*c;rhs=(gamma*c-1j*(energies+1j*.001))*s
            residuals.append([np.max(np.abs(lhs-rhs)/np.maximum(1,np.abs(lhs))),
                              np.max(np.abs(c*c+s*s-1)),np.min(c.real)])
        c_grid.append(c_row);s_grid.append(s_row)
    np.savez_compressed(args.catalogs/'complex_spectral_pilot.npz',energy=energies,delta=deltas,gamma=gammas,
                        c=np.asarray(c_grid),s=np.asarray(s_grid),eta=np.array(.001),
                        metadata_json=np.array(json.dumps({'energy_unit':'Delta0','causal':True,
                            'reference_hashes':sources,'builder_sha256':table.metadata['builder_sha256'],
                            'scope':'direct spectra at independent fields; no c/s interpolation'})))

    # eta convergence: the q=0 kernel integral has an exact finite-eta value.
    eta_rows=[]
    for eta in (.01,.003,.001,.0003):
        exact=np.pi/2-np.arctan(eta)
        for order in (201,801,3201):
            e=np.linspace(0,4,order)
            _,s=retarded_spectrum(e,delta=1,gamma=0,eta=eta)
            integral=np.trapezoid(2*s.real*s.imag,e)
            # Tail integral beyond E=4 is retained in the adaptive reference.
            def integrand(E):
                return 2*E*eta/((E*E-eta*eta-1)**2+4*E*E*eta*eta)
            adaptive=sum(quad(integrand,a,b,epsabs=1e-11,epsrel=1e-10,limit=200)[0]
                         for a,b in ((0,1),(1,4),(4,np.inf)))
            eta_rows.append([eta,order,integral,adaptive,exact,abs(integral-exact),
                             abs(adaptive-exact),abs(exact-np.pi/2)])
    np.savetxt(out/'spectral_eta_quadrature.csv',eta_rows,delimiter=',',
               header='eta,uniform_nodes,uniform_integral,adaptive_integral,finite_eta_exact,uniform_error,adaptive_error,eta_bias',comments='')

    fig,axes=plt.subplots(1,2,figsize=(10,3.8),layout='constrained')
    for i,label in enumerate(('Energía','Fuerza de amplitud','Conjugada a Γ')):
        axes[0].loglog(grid_rows[:,0],grid_rows[:,i+1],'o-',label=label)
        axes[1].loglog(nonthermal_rows[:,0],nonthermal_rows[:,i+1],'o-',label=label)
    axes[0].set(title='Vacío: 80 puntos fuera de los nodos',xlabel='Nodos por eje',ylabel='Error absoluto normalizado')
    axes[1].set(title='Ocupación no térmica: 6 campos',xlabel='Nodos por eje',ylabel='Error absoluto normalizado')
    for ax in axes:ax.legend(fontsize=8)
    fig.savefig(figures/'catalogue_refinement.png');fig.savefig(figures/'catalogue_refinement.pdf');plt.close(fig)
    fig,axes=plt.subplots(1,2,figsize=(10,3.8),layout='constrained')
    for gamma in (0,.2,.5,.8):
        c,_=retarded_spectrum(energies,delta=1,gamma=gamma,eta=.001)
        axes[0].plot(energies,c.real,label=f'Γ/Δ₀ = {gamma:g}')
    axes[0].set(xlim=(0,2.5),ylim=(0,4),xlabel='E/Δ₀',ylabel='DOS normalizada',title='Amplitud fija |Δ| = Δ₀')
    axes[0].legend(fontsize=8)
    eta_array=np.array(eta_rows)
    for order in (201,801,3201):
        rows=eta_array[eta_array[:,1]==order]
        axes[1].loglog(rows[:,0],rows[:,5],'o-',label=f'{order} nodos uniformes')
    unique=eta_array[eta_array[:,1]==201]
    axes[1].loglog(unique[:,0],unique[:,7],'k--',label='Sesgo de η con integral exacta')
    axes[1].set(xlabel='η/Δ₀',ylabel='Error de integral de corriente',title='Reducir η exige resolver el pico')
    axes[1].legend(fontsize=8)
    fig.savefig(figures/'spectral_eta.png');fig.savefig(figures/'spectral_eta.pdf');plt.close(fig)
    fig,ax=plt.subplots(figsize=(7.5,3.5),layout='constrained')
    thermal=np.array(thermal_rows)
    for eta in (.001,.0001):
        rows=thermal[thermal[:,3]==eta]
        ax.semilogy(np.arange(len(rows)),rows[:,6],'o-',label=f'η/Δ₀={eta:g}')
    ax.set(xlabel='Caso auxiliar: 3 espectros × 3 temperaturas',ylabel='Error absoluto en Uvac / (N₀Δ₀²)',
           title='Contraste independiente con Matsubara')
    ax.legend()
    fig.savefig(figures/'vacuum_reference.png');fig.savefig(figures/'vacuum_reference.pdf');plt.close(fig)
    report={
        'schema':'pysnspd.stage1.electronic_diagnostics.v1','host':platform.node(),
        'python':platform.python_version(),'numpy':np.__version__,'scipy':scipy.__version__,
        'runtime_seconds':time.perf_counter()-started,'production_connected':False,
        'scope':'uniform electronic pilot: vacuum and arbitrary fixed-p energy; no kinetics/PDE/transient',
        'field_normalization':'Delta0; U in N0 Delta0^2, amplitude/gamma derivatives in N0 Delta0',
        'occupation_field_support':{'delta':[.08,1.5],'gamma':[0,1.2],'normal':'separate exact delta=0 state',
                         'count_max':6,'count_order':64,'eta':.001},
        'catalogue_axes':{
            'vacuum':{'delta_nodes':table.delta_axis.tolist(),'gamma_nodes':table.gamma_axis.tolist()},
            'occupation':{'delta_nodes':occupation_table.vacuum.delta_axis.tolist(),
                          'gamma_nodes':occupation_table.vacuum.gamma_axis.tolist(),
                          'count_nodes':occupation_table.count_nodes.tolist(),
                          'count_weights':occupation_table.count_weights.tolist()},
        },
        'material':{'Tc_K':8.65,'D_m2_s':D,'sigma_n_S_m':sigma,'N0_per_J_m3':N0,
                    'Delta0_J':table.delta0_J,'BCS_ratio_exact':BCS_GAP_RATIO},
        'vacuum_grid_refinement':grid_rows.tolist(),
        'vacuum_identity_max_abs_errors':np.max(identity,axis=0).tolist(),
        'nonthermal_grid_refinement':nonthermal_rows.tolist(),
        'low_energy_grid_refinement':np.asarray(low_energy_rows).tolist(),
        'low_energy_occupation':'0.2 exp(-x/(0.25 Delta0)); includes Gamma=0 and 1e-4 Delta0',
        'low_energy_points':low_energy_points.tolist(),
        'thermal_fixed_p_reference':{'scope':'FD occupation evaluated at each state, held fixed for force; compared with independent Matsubara',
                                    'csv':'thermal_fixed_p_reference.csv','energy_cutoff_in_count':8},
        'nonthermal_identity_max_abs_errors':np.max(nonthermal_identity,axis=0).tolist(),
        'nonthermal_occupation':'0.16 exp(-((x/Delta0-0.95)/0.3)^2), held fixed under field variations',
        'nonthermal_reference':'same p and quadrature, direct causal E(x), R2/N1 and -N2R2/N1',
        'nonthermal_points':nonthermal_points.tolist(),
        'spectral_state_count':int(len(deltas)*len(gammas)*len(energies)),
        'spectral_max_equation_residual':float(np.max(np.array(residuals)[:,0])),
        'spectral_max_absolute_normalization_residual':float(np.max(np.array(residuals)[:,1])),
        'spectral_min_DOS':float(np.min(np.array(residuals)[:,2])),
        'vacuum_auxiliary_T_max_errors':{str(eta):float(np.max(thermal[thermal[:,3]==eta,6])) for eta in (.001,.0001)},
        'eta_adaptive_max_error':float(np.max(eta_array[:,6])),
        'reference_hashes':sources,
        'limitations':['finite count support and numerical eta; refine together before material prediction',
                       'field interpolation is tested at specified points, not globally certified',
                       'phonon admission, collision operators and spatial stability are independent'],
    }
    (out/'electronic_diagnostics.json').write_text(json.dumps(report,indent=2,ensure_ascii=False,allow_nan=False)+'\n',encoding='utf-8')
    print(json.dumps(report,indent=2,ensure_ascii=False))


if __name__=='__main__':
    main()
