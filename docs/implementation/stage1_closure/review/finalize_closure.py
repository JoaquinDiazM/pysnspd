"""Bind the independent electronic verdict and restricted material verdict."""
from pathlib import Path
from datetime import datetime,timezone
import hashlib
import json

ROOT=Path(__file__).resolve().parents[4]
OUT=ROOT/'docs/implementation/stage1_closure'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read(path):
    return json.loads(path.read_text(encoding='utf-8'))


def main():
    criteria=read(OUT/'acceptance_criteria.json')
    electronic=read(OUT/'review/electronic_assessment.json')
    cells=read(OUT/'cells_results.json')
    material=read(OUT/'material_results.json')
    replay=read(OUT/'review/material_replay.json')
    assert electronic['criteria_sha256']==sha(OUT/'acceptance_criteria.json')
    assert electronic['evidence']['cells']['sha256']==sha(OUT/'cells_results.json')
    assert electronic['cell_source_sha256']==sha(ROOT/'pysnspd/experimental/cell_validation.py')
    assert replay['material_results_sha256']==sha(OUT/'material_results.json')
    assert replay['derived_csv_sha256']==sha(OUT/'material_nbn_shape_v1.csv')
    assert replay['module_sha256']==material['module_sha256']==sha(ROOT/'pysnspd/experimental/material_preprocessing.py')
    assert material['runner_sha256']==sha(ROOT/'sandbox/stage1_closure/review_material.py')
    changes=[]
    def add(group,name,effect):
        changes.append(dict(group=group,name=name,absolute_change=effect['absolute_change'],
                            relative_change=effect['relative_change'],small_effect=effect['small_effect']))
    for row in material['thermal']:
        for name,effect in row['derived_vs_signed_reference'].items():
            add('thermal',f'T={row["T_K_assuming_THz"]} K conditional; {name}',effect)
    for row in material['nonthermal']:
        add('nonthermal_DOS',row['profile'],row['derived_vs_signed_reference'])
        if 'coupling_moment_derived_vs_reference' in row:
            add('nonthermal_coupling',row['profile'],row['coupling_moment_derived_vs_reference'])
    for name,effect in material['alpha_moment_changes']['common_support_v1'].items():
        add('alpha_moments',name,effect)
    for row in material['normal_exchange']:
        add('exchange_proxy',f'Te={row["Te_K"]} K conditional',row['derived_vs_reference'])
    epass=electronic['status']=='PASS_REDUCED_ELECTRONIC_USE_TESTS'
    shape_pass=replay['status']=='PASS'
    all_small=all(row['small_effect'] for row in changes)
    evidence={}
    paths={
        'criteria':'docs/implementation/stage1_closure/acceptance_criteria.json',
        'electronic_assessment':'docs/implementation/stage1_closure/review/electronic_assessment.json',
        'cell_results':'docs/implementation/stage1_closure/cells_results.json',
        'cell_source':'pysnspd/experimental/cell_validation.py',
        'material_results':'docs/implementation/stage1_closure/material_results.json',
        'material_replay':'docs/implementation/stage1_closure/review/material_replay.json',
        'derived_material_shape':'docs/implementation/stage1_closure/material_nbn_shape_v1.csv',
        'derived_material_manifest':'docs/implementation/stage1_closure/material_nbn_shape_v1.manifest.json',
        'material_source':'pysnspd/experimental/material_preprocessing.py',
        'R2_catalogue':criteria['immutable_inputs']['catalog_path'],
        'R2_query_source':criteria['immutable_inputs']['query_source_path']}
    for name,path in paths.items():
        evidence[name]=dict(path=path,sha256=sha(ROOT/path))
    fine=max(cells['two_cell_transport']['resolutions'],key=lambda row:row['nodes'])
    result=dict(schema='pysnspd.stage1_closure.admission.v1',created_utc=datetime.now(timezone.utc).isoformat(),
                status='CLOSED_ELECTRONIC_STAGE1_WITH_REDUCED_USE_TESTS; MATERIAL_RESTRICTED' if epass and shape_pass else 'NOT_CLOSED',
                electronic_status=electronic['status'],electronic_gate_count=electronic['gate_count'],
                electronic_stage1_closed=epass,
                material_status='CONDITIONAL_NUMERICAL_SHAPE_ONLY' if shape_pass else 'DERIVED_SHAPE_REJECTED',
                global_material_small_impact_pass=all_small,
                absolute_NbN_admitted=False,global_physical_phonon_kinetics_admitted=False,
                complete_D4_item2_dynamics=False,production_promotion=False,
                electronic_scope='Native R2 energy for BGK, heating and one self-consistent amplitude. Two-cell fixed-field transport conserves explicitly reconstructed energy; native energy bias, remapping and QP-count defect are distinct measured errors.',
                sampled_transport_summary=dict(native_energy_relative_bias=fine['native_R2_reconstruction_relative_max'],
                                               FD_remapping_relative_error=fine['FD_native_remap_relative_max'],
                                               QP_count_relative_drift=fine['quasiparticle_count_relative_drift'],
                                               reconstructed_energy_scaled_error=fine['reconstructed_energy_scaled_max']),
                material_scope='The authorized derived shape is numerically traceable and has a finite common-support ratio. It is not an equivalent replacement for arbitrary raw populations: global moment and tail tests fail the frozen 0.1% small-impact criterion. Temperature/meV statements remain conditional on the THz-axis assumption; DOS normalization and absolute rates are not certified.',
                material_small_impact_threshold=criteria['material']['small_impact_relative_max'],
                material_observable_changes=changes,
                material_failed_small_impact_observables=[dict(group=row['group'],name=row['name']) for row in changes if not row['small_effect']],
                material_restricted_examples='Low-temperature DOS-weighted energy/capacity and the explicitly inspected low-frequency region may be used for conditional shape diagnostics. No inference extends this to arbitrary nonthermal profiles or photon cascades.',
                remaining_work=['Shared superconducting electron-phonon events with exact reaction energy and Pauli/phonon positivity',
                                'Simultaneous two-cell moving-condensate, collision and fixed-energy transport coupling',
                                'Material mobility and kinetic-time identification, phonon-unit/basis provenance',
                                'Spatial principal-symbol/core checks, interfaces, reservoirs and circuit before full transients'],
                next_permitted_implementation='Complete one-and-two-cell self-consistent electron-phonon coupling using an explicitly synthetic admissible phonon input; retain material-shape sensitivity as a separate diagnostic.',
                evidence=evidence,reviewer_source_sha256=sha(Path(__file__)),
                final_regression_policy='The delivery manifest binds the separate full-suite log; this scientific verdict does not turn that software regression into microscopic or material validation.')
    (OUT/'closure_admission.json').write_text(json.dumps(result,indent=2,ensure_ascii=False)+'\n',encoding='utf-8')
    print(json.dumps({key:result[key] for key in ('status','electronic_gate_count','material_status','global_material_small_impact_pass')},indent=2))


if __name__=='__main__':
    main()
