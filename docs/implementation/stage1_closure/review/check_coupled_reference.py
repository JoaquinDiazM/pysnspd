"""Independent Gamma=0 moving-amplitude reference on the native count rule.

Reference energies, amplitude kernels and vacuum force are analytic finite-eta
BCS expressions, never catalogue queries. DOP853 is separate from the runner's
RK4. Synthetic mobility/time do not certify the material KWT law.
"""
from pathlib import Path
import hashlib
import json
import sys
import time

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
from scipy.special import expit

ROOT=Path(__file__).resolve().parents[4]
sys.path.insert(0,str(ROOT))
from pysnspd.experimental.energy_catalog import OccupationEnergyCatalog
from pysnspd.experimental.cell_validation import ElectronicCell


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    started=time.perf_counter()
    path=ROOT/'docs/implementation/stage1_r2/catalogs/occupation_catalog.npz'
    cat=OccupationEnergyCatalog.load(path)
    x,w,eta=cat.count_nodes,cat.count_weights,cat.eta
    mobility,tau,bath,duration=.2,.7,.08,1.5
    initial=np.r_[.60,.035*np.exp(-((x-.5)/.2)**2)]

    def analytic_fields(a):
        d=x*x+eta*eta
        e=x*np.sqrt(1+a*a/d)
        ea=x*a/np.sqrt(d*(d+a*a))
        return e,ea

    def analytic_energy(state):
        a,p=state[0],state[1:]
        e,_=analytic_fields(a)
        return a*a*(np.log(a)-.5)+4*np.dot(w*e,p)

    def reference_rhs(time,state):
        a,p=state[0],state[1:]
        e,ea=analytic_fields(a)
        target=4*np.dot(w*e,p)
        temperature=brentq(lambda t: -target if t==0 else
                          4*np.dot(w*e,expit(-e/t))-target,0,10,xtol=1e-14,rtol=1e-14)
        fd=expit(-e/temperature)
        force=2*a*np.log(a)+4*np.dot(w*ea,p)
        adot=-mobility*force
        power=mobility*force*force
        # Direct FD formula is well conditioned for this declared trajectory.
        shape=e*expit(-e/max(temperature,bath))*(1-p)
        heat=power*shape/(4*np.dot(w*e,shape))
        return np.r_[adot,(fd-p)/tau+heat]

    def catalogue_rhs(time,state):
        a,p=state[0],state[1:]
        cell=ElectronicCell(cat,float(a),0.)
        force=cell.moments(p)[1]
        power=mobility*force*force
        return np.r_[-mobility*force,cell.bgk(p,tau)[0]+cell.heating(p,power,bath)]

    times=np.linspace(0,duration,81)
    reference=solve_ivp(reference_rhs,(0,duration),initial,method='DOP853',t_eval=times,
                        rtol=2e-11,atol=2e-13,max_step=.08)
    measured=solve_ivp(catalogue_rhs,(0,duration),initial,method='DOP853',t_eval=times,
                       rtol=2e-11,atol=2e-13,max_step=.08)
    if not reference.success or not measured.success:
        raise RuntimeError('independent coupled trajectory did not finish')
    e0=analytic_energy(initial)
    residual=max(abs(analytic_energy(state)-e0) for state in measured.y.T)
    population_error=float(np.max(abs(reference.y[1:]-measured.y[1:])))
    amplitude_error=float(np.max(abs(reference.y[0]-measured.y[0])))
    pauli=(float(np.min(measured.y[1:])),float(np.max(measured.y[1:])))
    passed=residual<=1e-7 and population_error<=1e-7 and amplitude_error<=1e-7 and pauli[0]>=0 and pauli[1]<=1
    output=dict(schema='pysnspd.stage1_closure.independent_coupled_reference.v1',status='PASS' if passed else 'FAIL',
                scope='One synthetic Gamma=0 cell; analytic BCS reference at the same native quadrature and eta, not cutoff convergence or NbN dynamics.',
                settings=dict(mobility=mobility,tau=tau,bath_temperature=bath,duration=duration,
                              solver='DOP853',relative_tolerance=2e-11,absolute_tolerance=2e-13,max_step=.08),
                criteria_sha256=sha(ROOT/'docs/implementation/stage1_closure/acceptance_criteria.json'),
                catalog_sha256=sha(path),cell_source_sha256=sha(ROOT/'pysnspd/experimental/cell_validation.py'),
                reviewer_source_sha256=sha(Path(__file__)),
                summary=dict(energy_absolute_residual=residual,population_absolute_error=population_error,
                             amplitude_absolute_error=amplitude_error,pauli_minimum=pauli[0],pauli_maximum=pauli[1],
                             initial_amplitude=float(initial[0]),final_amplitude=float(measured.y[0,-1])),
                time=times.tolist(),amplitude=measured.y[0].tolist(),
                analytic_energy=[analytic_energy(state) for state in measured.y.T],
                runtime_seconds=time.perf_counter()-started)
    Path(__file__).with_name('coupled_reference.json').write_text(json.dumps(output,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({key:output[key] for key in ('status','summary','runtime_seconds')},indent=2))


if __name__=='__main__':
    main()
