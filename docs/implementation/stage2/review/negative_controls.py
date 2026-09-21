"""Intentional diagnostic-only errors prove that the frozen gates detect them.

No incorrect variant is inserted into an operator or accepted trajectory.
"""
from pathlib import Path
import json
import sys
import time
import numpy as np
from continuous_reactions import ROOT,sha
sys.path.insert(0,str(ROOT))
from pysnspd.experimental.energy_catalog import OccupationEnergyCatalog,K_B_J_K
from pysnspd.experimental.cell_validation import ElectronicCell
from pysnspd.experimental.cell_closures import CellScales,KWTMobility


def main():
    start=time.perf_counter();here=Path(__file__).resolve().parent
    source=ROOT/'docs/implementation/stage1_r2/catalogs/occupation_catalog.npz'
    cat=OccupationEnergyCatalog.load(source)
    scales=CellScales(cat.vacuum.delta0_J,cat.vacuum.N0_per_J_m3,8.65,
                      .12*cat.vacuum.delta0_J/K_B_J_K,1.)
    mobility=KWTMobility(scales);rows=[]
    for a,g in ((.6,0.),(.55,0.),(.95,0.),(.72,.2)):
        cell=ElectronicCell(cat,a,g)
        p=.2*np.exp(-((cat.count_nodes-.8)/.3)**2)+cell.fermi_dirac(.12)
        force=cell.moments(p)[1]
        move=mobility.amplitude_response(a,force,cell.equivalent_temperature(p))
        dp=cell.heating(p,move.heat,.12)
        deposited=float(4*np.dot(cell.weights*cell.energies,dp))
        work=force*move.velocity
        scale=max(1.,abs(work)+abs(deposited))
        correct=abs(work+deposited)/scale
        doubled=abs(work+2*deposited)/scale
        missing_work=abs(deposited)/scale
        rows.append(dict(amplitude=a,gamma=g,force=force,velocity=move.velocity,
                         condensate_heat=move.heat,actual_deposited_heat=float(deposited),
                         correct_energy_residual=correct,
                         incorrect_double_condensate_heat_residual=doubled,
                         incorrect_omitted_field_work_residual=missing_work,
                         correct_pass=correct<=1e-10,
                         double_heat_rejected=doubled>100*1e-10,
                         omitted_work_rejected=missing_work>100*1e-10))
    # An independently computed nonzero continuous rate checks the absolute
    # collision factor. Energy conservation alone cannot detect a common factor.
    ref=json.loads((here/'continuous_reactions.json').read_text())
    hot=next(r for r in ref['rows'] if r['amplitude']==0 and r['profile']=='hot_electrons')
    expected=hot['channels']['scattering']['reference'][0]
    wrong=.5*expected
    factor=dict(reference_scattering_power=expected,intentionally_wrong_power=wrong,
                wrong_relative_error=abs(wrong-expected)/abs(expected),
                gate=1e-3,rejected=abs(wrong-expected)/abs(expected)>1e-3,
                caveat='Diagnostic factor perturbation of an independent reference; this is not a new production operator')
    gates=dict(correct_balance=all(r['correct_pass'] for r in rows),
               double_heat_detected=all(r['double_heat_rejected'] for r in rows),
               omitted_work_detected=all(r['omitted_work_rejected'] for r in rows),
               common_collision_factor_detected=factor['rejected'])
    data=dict(schema='pysnspd.stage2.independent_negative_controls.v1',
        status='PASS' if all(gates.values()) else 'FAIL',gates=gates,rows=rows,
        factor_control=factor,
        criteria_sha256=sha(ROOT/'docs/implementation/stage2/acceptance_criteria.json'),
        catalog_sha256=sha(source),source_sha256=sha(Path(__file__)),
        module_sha256={name:sha(ROOT/'pysnspd/experimental'/name) for name in ('energy_catalog.py','cell_validation.py','cell_closures.py')},
        runtime_seconds=time.perf_counter()-start)
    (here/'negative_controls.json').write_text(json.dumps(data,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(dict(status=data['status'],gates=gates,runtime_seconds=data['runtime_seconds'],
                          minimum_bad_residual=min(r['incorrect_double_condensate_heat_residual'] for r in rows))))
    if not all(gates.values()):raise SystemExit(1)


if __name__=='__main__':main()
