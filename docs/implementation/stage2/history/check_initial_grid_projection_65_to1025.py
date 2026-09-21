"""Cheap initial-state projection diagnostic; no time integration or admission."""
from pathlib import Path
import hashlib
import json
import time
import numpy as np
from assess_grid_refinement import projected_mass
from run_coupled import setup

ROOT=Path(__file__).resolve().parents[2]


def main():
    started=time.perf_counter();levels=[65,129,257,513,1025]
    output=ROOT/'docs/implementation/stage2/review/initial_phonon_projection.json'
    profile='run_coupled.setup: one driven cell, Bose(.12) plus fixed0.15 injected phonon energy'
    source_paths=[Path(__file__),ROOT/'sandbox/stage2_cells/run_coupled.py',
        ROOT/'sandbox/stage2_cells/assess_grid_refinement.py',ROOT/'pysnspd/experimental/cell_closures.py']
    registration={'phonon_nodes':levels,'infrared':.005,'electron_refinement':1,
        'profile':profile,'projection':'65 positive linear hats on Omega[0,4]',
        'scope':'Initial-state diagnostic only, no time integration or admission',
        'sources':{p.relative_to(ROOT).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in source_paths}}
    output.parent.mkdir(parents=True,exist_ok=True)
    output.with_name(output.stem+'_scenarios.json').write_text(json.dumps(registration,indent=2)+'\n',encoding='utf-8')
    masses=[];rows=[]
    for nodes in levels:
        system,initial,_=setup('one',phonon_nodes=nodes,infrared=.005,electron_refinement=1)
        _,_,n=system.unpack(initial)
        mass=projected_mass(system.phonons.energies,n,system.phonons.capacities,np.linspace(0,4,65))
        masses.append(mass)
        rows.append({'phonon_nodes':nodes,'total_projected_phonon_number':float(np.sum(mass)),
            'total_projected_phonon_energy':float(np.sum(mass*np.linspace(0,4,65)))})
    ref=masses[-1]
    for row,mass in zip(rows,masses):row['relative_L1_projection_error_vs_1025']=float(np.sum(abs(mass-ref))/np.sum(abs(ref)))
    result={**registration,'rows':rows,'runtime_seconds':time.perf_counter()-started}
    output.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({'rows':rows,'runtime_seconds':result['runtime_seconds']},indent=2))


if __name__=='__main__':main()
