"""Apply the frozen closure gates to exact saved evidence; never fit tolerances."""
from pathlib import Path
import hashlib
import json
import math

ROOT=Path(__file__).resolve().parents[4]
OUT=ROOT/'docs/implementation/stage1_closure'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load(path):
    return json.loads(path.read_text(encoding='utf-8'))


def main():
    criteria=load(OUT/'acceptance_criteria.json')
    data=load(OUT/'cells_results.json')
    gates=[]
    def gate(name,value,limit=None,*,passed=None,detail=None):
        okay=bool(passed) if passed is not None else math.isfinite(value) and value<=limit
        gates.append(dict(name=name,passed=okay,value=value,limit=limit,detail=detail))
    def pauli(name,rows):
        minimum=min(row['minimum_population'] for row in rows)
        maximum=max(row['maximum_population'] for row in rows)
        gate(name,[minimum,maximum],passed=minimum>=-1e-13 and maximum<=1+1e-13,
             detail='No population clipping; tolerance only distinguishes arithmetic roundoff.')
    def reduction(name,rows,error_key,minimum):
        ordered=sorted(rows,key=lambda row:row['steps'])
        for coarse,fine in zip(ordered[:-1],ordered[1:]):
            c,f=coarse[error_key],fine[error_key]
            floor=100*2.220446049250313e-16
            limited=c<=floor or f<=floor
            ratio=c/f if f else None
            gate(f'{name}_{coarse["steps"]}_to_{fine["steps"]}',ratio,
                 passed=limited or (ratio is not None and ratio>=minimum),
                 detail=dict(minimum_reduction=minimum,coarse_error=c,fine_error=f,roundoff_limited=limited))

    chash=sha(OUT/'acceptance_criteria.json')
    gate('fresh_cells_criteria_hash',data['criteria_sha256'],passed=data['criteria_sha256']==chash)
    for path_key,hash_key,data_key in [('catalog_path','catalog_sha256','catalog_sha256'),
                                     ('query_source_path','query_source_sha256','query_code_sha256')]:
        expected=criteria['immutable_inputs'][hash_key]
        got=sha(ROOT/criteria['immutable_inputs'][path_key])
        gate(hash_key,got,passed=got==expected==data[data_key])
    cell_hash=sha(ROOT/'pysnspd/experimental/cell_validation.py')
    gate('measured_cell_code_hash',data['cell_code_sha256'],passed=data['cell_code_sha256']==cell_hash)
    gate('measured_runner_hash',data['runner_sha256'],passed=data['runner_sha256']==sha(ROOT/'sandbox/stage1_closure/run_cells.py'))
    gate('production_disconnected',data['production_connected'],passed=data['production_connected'] is False)
    recheck=load(OUT/'catalog_recheck.json')
    gate('fresh_R2_recheck',recheck['status'],passed=recheck['status']=='PASS_TESTED_PHYSICAL_AND_API_GATES')

    b=data['bgk']; bc=criteria['electronic']['bgk']
    finest_steps=max(row['steps'] for row in b['resolutions'])
    bfinest=[row for row in b['resolutions'] if row['steps']==finest_steps]
    for key,limit in [('energy_scaled_error',bc['energy_scaled_max']),
                      ('instantaneous_energy_moment_max',bc['instantaneous_energy_moment_scaled_max']),
                      ('final_population_error',bc['exact_solution_population_abs_max']),
                      ('entropy_decrease_max',bc['entropy_decrease_abs_max'])]:
        gate('bgk_'+key,max(row[key] for row in bfinest),limit,detail='Finest time grid; coarse grids retained for independent refinement ratios.')
    pauli('bgk_Pauli',b['resolutions'])
    for profile in sorted(set(row['profile'] for row in b['resolutions'])):
        reduction('bgk_'+profile,[row for row in b['resolutions'] if row['profile']==profile],
                  'final_population_error',bc['time_refinement']['minimum_error_reduction_per_halving'])
    gate('bgk_vacuum_stationarity',b['stationarity']['vacuum_rhs_max'],1e-12)
    gate('bgk_FD_stationarity',b['stationarity']['thermal_rhs_max'],1e-12)

    m=data['moving_spectrum']; mc=criteria['electronic']['moving_spectrum_and_heat']
    gate('moving_energy_ledger',max(m['resolutions'],key=lambda row:row['steps'])['energy_ledger_scaled_max'],mc['energy_ledger_scaled_max'])
    gate('moving_heat_energy_moment',m['source_energy_moment_max'],mc['source_energy_moment_scaled_max'])
    gate('moving_thermal_direction',m['thermal_direction_abs_error'],mc['thermal_direction_abs_max'])
    negative=min(row['negative_control_final_abs'] for row in m['resolutions'])
    gate('moving_negative_control',negative,passed=negative>=mc['negative_control_minimum_multiple_of_gate']*mc['energy_ledger_scaled_max'])
    reduction('moving_time',m['resolutions'],'energy_ledger_scaled_max',mc['time_refinement_minimum_reduction'])
    pauli('moving_Pauli',m['resolutions'])

    s=data['self_consistent_cell']; sc=criteria['electronic']['self_consistent_single_cell']
    gate('self_consistent_energy',max(s['resolutions'],key=lambda row:row['steps'])['energy_ledger_scaled_max'],sc['energy_ledger_scaled_max'])
    gate('self_consistent_instantaneous_balance',s['instantaneous_balance_max'],sc['instantaneous_balance_scaled_max'])
    gate('self_consistent_positive_dissipation',s['minimum_dissipation'],passed=s['minimum_dissipation']>=-sc['negative_dissipation_tolerance'])
    reduction('self_consistent_time',s['resolutions'],'energy_ledger_scaled_max',sc['time_refinement_minimum_reduction'])
    pauli('self_consistent_Pauli',s['resolutions'])

    t=data['two_cell_transport']; tc=criteria['electronic']['two_cell_transport']
    rows=sorted(t['resolutions'],key=lambda row:row['nodes']); finest=rows[-1]
    gate('transport_reconstructed_energy',max(row['reconstructed_energy_scaled_max'] for row in rows),tc['energy_scaled_max'])
    gate('transport_same_physical_distribution',max(row['common_distribution_rhs_max'] for row in rows),tc['common_distribution_flux_abs_max'])
    gate('transport_native_R2_energy_bias',finest['native_R2_reconstruction_relative_max'],tc['native_R2_reconstruction_energy_relative_max'])
    for key in ('FD_native_remap_relative_max','quasiparticle_count_relative_drift'):
        gate('transport_'+key,finest[key],tc['remap_and_count_relative_max'])
        errors=[row[key] for row in rows]
        gate('transport_refinement_'+key,errors,passed=all(f<c or c<100*2.220446049250313e-16 for c,f in zip(errors[:-1],errors[1:])))
    gate('transport_hot_to_cold',finest['energy_transferred_left_to_right'],passed=finest['energy_transferred_left_to_right']>0)
    pauli('transport_Pauli',rows)
    time_rows=t['time_resolutions']
    gate('transport_time_error',max(time_rows,key=lambda row:row['steps'])['final_state_scaled_error'],tc['time_final_state_scaled_max'])
    reduction('transport_time',time_rows,'final_state_scaled_error',tc['time_refinement_minimum_reduction'])
    pauli('transport_time_Pauli',time_rows)

    evidence={}
    for name in ('operator_checks','transport_checks','coupled_reference'):
        path=OUT/'review'/f'{name}.json'
        result=load(path)
        gate('independent_'+name,result['status'],passed=result['status']=='PASS')
        gate('independent_'+name+'_code_hash',result['cell_source_sha256'],passed=result['cell_source_sha256']==cell_hash)
        gate('independent_'+name+'_criteria_hash',result['criteria_sha256'],passed=result['criteria_sha256']==chash)
        evidence[name]=dict(path=path.relative_to(ROOT).as_posix(),sha256=sha(path))
    for name,relative in [('criteria','acceptance_criteria.json'),('cells','cells_results.json'),('catalog_recheck','catalog_recheck.json')]:
        path=OUT/relative
        evidence[name]=dict(path=path.relative_to(ROOT).as_posix(),sha256=sha(path))
    passed=all(row['passed'] for row in gates)
    result=dict(schema='pysnspd.stage1_closure.electronic_assessment.v1',
                status='PASS_REDUCED_ELECTRONIC_USE_TESTS' if passed else 'FAIL_REDUCED_ELECTRONIC_USE_TESTS',
                scope='Actual R2 native energy for BGK, heat and amplitude; reconstructed energy only for the transport prototype. No simultaneous electron-phonon or full spatial dynamics.',
                criteria_sha256=chash,cell_source_sha256=cell_hash,reviewer_source_sha256=sha(Path(__file__)),
                gates=gates,gate_count=len(gates),failed_gates=[row['name'] for row in gates if not row['passed']],
                evidence=evidence,
                next_permitted_work='Complete self-consistent electronic/phonon one-and-two-cell coupling, with a shared event energy ledger, before spatial and circuit integration.')
    (OUT/'review'/'electronic_assessment.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    print(json.dumps({key:result[key] for key in ('status','gate_count','failed_gates')},indent=2))


if __name__=='__main__':
    main()
