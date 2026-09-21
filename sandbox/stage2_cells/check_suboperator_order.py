"""Check RK4 order against exact exponentials for two smooth suboperators.

BGK uses ElectronicCell.bgk, including its energy-temperature inversion at
every RK stage. Escape isolates the same linear term and energy ledger as
CoupledCellSystem.rhs. No coupled transient, population repair or modified
experimental module is used.
"""
from __future__ import annotations
import argparse
import hashlib
import inspect
import json
from pathlib import Path
import platform
import sys
import time

import numpy as np
import scipy

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT))
from pysnspd.experimental import cell_validation,cell_closures,energy_catalog
from pysnspd.experimental.cell_validation import ElectronicCell,rk4_trajectory
from pysnspd.experimental.cell_closures import CellScales,DebyePhonons,bose_occupation


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def write(path,value):
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(value,indent=2,ensure_ascii=False,allow_nan=False)+'\n',encoding='utf-8')


def summarize(rows,key):
    ratios=[rows[i][key]/rows[i+1][key] for i in range(len(rows)-1)]
    return {'error_reduction_ratios':ratios,'observed_orders':[float(np.log2(v)) for v in ratios],
            'passes_RK4_reduction_with_25_percent_margin':all(12<=v<=20 for v in ratios)}


def bgk_case(catalog):
    cell=ElectronicCell(catalog,.72,.2)
    initial=.4*cell.fermi_dirac(.15)+.6*cell.fermi_dirac(.5)
    temperature=cell.equivalent_temperature(initial)
    equilibrium=cell.fermi_dirac(temperature)
    duration,tau=1.4,.7
    exact=equilibrium+(initial-equilibrium)*np.exp(-duration/tau)
    capacity=4*cell.weights
    scale=float(np.dot(capacity,np.abs(initial-equilibrium)))
    initial_energy=cell.excitation_energy(initial)
    rows=[]
    for steps in [10,20,40]:
        bounds=[float(initial.min()),float(initial.max())]; temperatures=[]
        def rhs(t,state):
            bounds[0]=min(bounds[0],float(state.min()));bounds[1]=max(bounds[1],float(state.max()))
            value,theta=cell.bgk(state,tau)
            temperatures.append(theta)
            return value
        _,states=rk4_trajectory(rhs,initial,duration,steps)
        absolute=float(np.dot(capacity,np.abs(states[-1]-exact)))
        energy_drift=max(abs(cell.excitation_energy(p)-initial_energy) for p in states)
        rows.append({'steps':steps,'time_step':duration/steps,'weighted_L1_error':absolute,
            'relative_to_initial_relaxing_population':absolute/scale,
            'maximum_energy_drift_scaled':energy_drift/max(1.,abs(initial_energy)),
            'maximum_equivalent_temperature_drift':max(abs(t-temperature) for t in temperatures),
            'minimum_stage_population':bounds[0],'maximum_stage_population':bounds[1],
            'RHS_calls':len(temperatures)})
    floor=100*np.finfo(float).eps*max(1.,float(np.dot(capacity,np.abs(exact))))
    order=summarize(rows,'weighted_L1_error')
    return {'field':[.72,.2],'tau':tau,'duration':duration,'energy_matched_temperature':temperature,
        'exact_solution':'p(t)=f_FD(T_equal_energy)+(p(0)-f_FD(T_equal_energy))*exp(-t/tau_kin)',
        'initial_relaxing_population_norm':scale,'roundoff_floor':floor,'rows':rows,**order,
        'above_roundoff_floor':all(r['weighted_L1_error']>floor for r in rows),
        'valid_populations_at_every_stage':all(0<=r['minimum_stage_population']<=r['maximum_stage_population']<=1 for r in rows),
        'energy_conservation_pass':all(r['maximum_energy_drift_scaled']<=1e-11 for r in rows)}


def escape_case(catalog):
    vacuum=catalog.vacuum
    scales=CellScales(vacuum.delta0_J,vacuum.N0_per_J_m3,8.65,.12*vacuum.delta0_J/energy_catalog.K_B_J_K,1.)
    model=DebyePhonons(scales,4.,10*vacuum.N0_per_J_m3*vacuum.delta0_J,.1,'RK4 order: synthetic Debye',.005)
    grid=model.quadrature(65)
    bath=bose_occupation(grid.energies,.12)
    initial=bose_occupation(grid.energies,.35)
    duration,tau=30.,15.
    exact=bath+(initial-bath)*np.exp(-duration/tau)
    exact_escaped=grid.energy(initial-exact)
    initial_energy=grid.energy(initial)
    norm=float(np.dot(grid.capacities,np.abs(initial-bath)))
    rows=[]
    for steps in [10,20,40]:
        minimum=[float(initial.min())];calls=[0]
        def rhs(t,state):
            n=state[:-1]
            if np.any(n<0) or np.any(~np.isfinite(state)):
                raise ValueError('invalid stage phonon occupation; no clipping')
            minimum[0]=min(minimum[0],float(n.min()));calls[0]+=1
            escape=(n-bath)/tau
            # This sign convention matches the recorded escaping-energy ledger.
            return np.r_[-escape,np.dot(grid.capacities*grid.energies,escape)]
        _,states=rk4_trajectory(rhs,np.r_[initial,0.],duration,steps)
        absolute=float(np.dot(grid.capacities,np.abs(states[-1,:-1]-exact)))
        ledger_error=abs(float(states[-1,-1])-exact_escaped)
        energy_drift=max(abs(grid.energy(state[:-1])+state[-1]-initial_energy) for state in states)
        rows.append({'steps':steps,'time_step':duration/steps,'weighted_L1_error':absolute,
            'relative_to_initial_relaxing_population':absolute/norm,
            'escaped_energy_absolute_error':ledger_error,'escaped_energy_relative_error':ledger_error/abs(exact_escaped),
            'maximum_energy_ledger_residual_scaled':energy_drift/max(1.,abs(initial_energy)),
            'minimum_stage_population':minimum[0],'RHS_calls':calls[0]})
    floor=100*np.finfo(float).eps*max(1.,float(np.dot(grid.capacities,np.abs(exact))))
    order=summarize(rows,'weighted_L1_error')
    ledger_order=summarize(rows,'escaped_energy_absolute_error')
    return {'tau':tau,'duration':duration,'initial_temperature_bar':.35,'bath_temperature_bar':.12,
        'exact_solution':'n(t)=n_bath+(n(0)-n_bath)*exp(-t/tau_escape); U_escape=U_ph(0)-U_ph(t)',
        'implementation_scope':'The isolated escape term and escaping-energy ledger from CoupledCellSystem.rhs; complete coupled dynamics are not invoked.',
        'material':model.metadata(),'roundoff_floor':floor,'rows':rows,**order,
        'escaped_energy_order':ledger_order,'above_roundoff_floor':all(r['weighted_L1_error']>floor for r in rows),
        'valid_populations_at_every_stage':all(r['minimum_stage_population']>=0 for r in rows),
        'energy_ledger_pass':all(r['maximum_energy_ledger_residual_scaled']<=1e-11 for r in rows)}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,default=ROOT/'docs/implementation/stage2/time_suboperator_order.json')
    args=parser.parse_args();started=time.perf_counter()
    catalogue=ROOT/'docs/implementation/stage1_r2/catalogs/occupation_catalog.npz'
    sources=[Path(cell_validation.__file__),Path(cell_closures.__file__),Path(energy_catalog.__file__),
             ROOT/'pysnspd/experimental/coupled_cells.py',ROOT/'docs/implementation/stage2/acceptance_criteria.json']
    registration={'schema':'pysnspd.stage2.smooth-time-scenarios.v1','steps':[10,20,40],
        'bgk':{'field':[.72,.2],'initial':'0.4*FD(.15)+0.6*FD(.5)','tau':.7,'duration':1.4},
        'escape':{'initial':'Bose(.35)','bath':.12,'tau':15.,'duration':30.,'phonon_nodes':65,
                  'synthetic_Debye':{'cutoff':4.,'IR':.005,'lambda':.1,'atom_density':'10*N0*Delta0'}},
        'RK4_expected_error_ratio':16.,'allowed_ratio':[12.,20.],
        'roundoff_policy':'100 machine epsilons times max(1, exact weighted-population norm)',
        'sources':{p.relative_to(ROOT).as_posix():sha(p) for p in sources},
        'RK4_function_sha256':hashlib.sha256(inspect.getsource(rk4_trajectory).encode()).hexdigest(),
        'catalog_sha256':sha(catalogue),'runner_sha256':sha(__file__)}
    scenario_path=args.output.with_name(args.output.stem+'_scenarios.json')
    write(scenario_path,registration)
    catalog=energy_catalog.OccupationEnergyCatalog.load(catalogue)
    bgk,escape=bgk_case(catalog),escape_case(catalog)
    gates={'BGK_order':bgk['passes_RK4_reduction_with_25_percent_margin'],
        'escape_order':escape['passes_RK4_reduction_with_25_percent_margin'],
        'escape_ledger_order':escape['escaped_energy_order']['passes_RK4_reduction_with_25_percent_margin'],
        'BGK_above_roundoff':bgk['above_roundoff_floor'],'escape_above_roundoff':escape['above_roundoff_floor'],
        'BGK_valid_stages':bgk['valid_populations_at_every_stage'],'escape_valid_stages':escape['valid_populations_at_every_stage'],
        'BGK_energy':bgk['energy_conservation_pass'],'escape_energy':escape['energy_ledger_pass']}
    result={'schema':'pysnspd.stage2.smooth-time-order.v1','status':'PASS' if all(gates.values()) else 'FAIL',
        'gates':gates,'BGK':bgk,'escape':escape,'scenario_sha256':sha(scenario_path),
        'sources':registration['sources'],'RK4_function_sha256':registration['RK4_function_sha256'],
        'runner_sha256':sha(__file__),'catalog_sha256':sha(catalogue),
        'runtime_seconds':time.perf_counter()-started,'python':platform.python_version(),
        'numpy':np.__version__,'scipy':scipy.__version__,'platform':platform.platform()}
    write(args.output,result)
    print(json.dumps({'status':result['status'],'gates':gates,'BGK_orders':bgk['observed_orders'],
        'escape_orders':escape['observed_orders'],'escape_ledger_orders':escape['escaped_energy_order']['observed_orders'],
        'runtime_seconds':result['runtime_seconds']},indent=2))
    if not all(gates.values()):raise SystemExit(1)


if __name__=='__main__':main()
