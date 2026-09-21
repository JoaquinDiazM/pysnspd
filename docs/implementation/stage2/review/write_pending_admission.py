"""Assemble an evidence-bound checkpoint verdict; no numerical jobs are run."""
from pathlib import Path
import hashlib
import json

ROOT=Path(__file__).resolve().parents[4]
STAGE=ROOT/'docs/implementation/stage2'


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read(relative):
    return json.loads((STAGE/relative).read_text(encoding='utf-8'))


def artifact(relative):
    path=STAGE/relative
    return dict(path=path.relative_to(ROOT).as_posix(),sha256=digest(path))


def main():
    criteria=read('acceptance_criteria.json')
    criteria_hash=digest(STAGE/'acceptance_criteria.json')
    if criteria_hash!='48a56a76b64b58b81b176325a6a34cc2e6c690bfde3b4eb739a8cd2723c702cf':
        raise ValueError('The frozen contract changed.')
    immutable=[]
    for prefix in ('R2_catalogue','R2_query_source','stage1_closure_certificate','stage1_criteria'):
        path=ROOT/criteria['immutable_inputs'][prefix+'_path']
        expected=criteria['immutable_inputs'][prefix+'_sha256']
        actual=digest(path)
        immutable.append(dict(path=path.relative_to(ROOT).as_posix(),sha256=actual,matches_frozen_input=actual==expected))
    if not all(row['matches_frozen_input'] for row in immutable):
        raise ValueError('An immutable input changed.')
    event=read('event_summary.json')
    energy=read('review/complementary_energy_reference.json')
    moments=read('review/complementary_moment_quadrature.json')
    transport=read('transport_results.json')
    closure=read('closure_results.json')
    regression=read('regression_result.json')
    if regression['exit_code']!=0 or regression['hashes_before']!=regression['hashes_after']:
        raise ValueError('Regression suite failed or its sources changed during execution.')
    if any(digest(ROOT/path)!=expected for path,expected in regression['hashes_after'].items()):
        raise ValueError('Current regression sources differ from the completed suite.')
    algebra=read('projected_phonon_1025.json')['algebra']
    source_paths=['pysnspd/experimental/'+name+'.py' for name in
        ('cell_closures','cell_transport','cell_validation','coupled_cells','energy_catalog','kinetic_events','refined_cells')]
    source_hashes={name:digest(ROOT/name) for name in source_paths}
    if source_hashes['pysnspd/experimental/refined_cells.py']!=event['factory_sha256']:
        raise ValueError('Complementary factory differs from event evidence.')
    if source_hashes['pysnspd/experimental/kinetic_events.py']!=event['kernel_sha256']:
        raise ValueError('Event kernel differs from event evidence.')
    artifacts=[
        'acceptance_criteria.json','review/complementary_grid_amendment.json',
        'review/complementary_energy_reference.json','review/complementary_moment_quadrature.json',
        'review/continuous_reactions.json','review/continuous_projection_check.json',
        'review/infrared_reaction_reference.json','review/upper_support_reference.json',
        'review/negative_controls.json','review/weak_population_reference_final_all_1.json',
        'review/weak_population_reference_final_all_2.json','review/weak_population_reference_final_all_4.json',
        'review/weak_population_reference_finite_eta_bcs_1.json',
        'closure_results.json','transport_results.json','event_summary.json',
        'complementary_event_final.json','projected_phonon_convergence.json','projected_phonon_1025.json',
        'selected_event_measurements.json','weak_reference_final_comparison_1.json',
        'weak_reference_final_comparison_2.json','weak_reference_final_comparison_4.json',
        'time_suboperator_order.json','review/initial_phonon_projection_fine.json',
        'review/time_assessment_selfcheck.json','review/grid_assessment_selftest.json',
        'reference_handoff.json','trajectories/one_reference_probe.json','regression_result.json']
    finite_comparison=STAGE/'weak_reference_finite_eta_comparison.json'
    if finite_comparison.exists(): artifacts.append(finite_comparison.name)
    completed=[]
    for name in ('one_40','one_80','two_40','vacuum'):
        path=STAGE/f'trajectories/{name}.json'
        if not path.exists(): continue
        run=json.loads(path.read_text(encoding='utf-8'))
        npz=path.with_suffix('.npz')
        if digest(npz)!=run['trajectory_sha256']: raise ValueError('Trajectory hash mismatch: '+name)
        artifacts.append(f'trajectories/{name}.json')
        completed.append(dict(name=name,status='HISTORICAL_COMPLETED_NOT_ADMITTED',
            phonon_nodes=run['parameters']['phonon_nodes'],electron_states=run['occupation_mesh']['electron_states'],
            method=run['parameters']['method'],steps=run['parameters']['steps'],
            duration=run['parameters']['duration'],ledger_scaled_max=run['energy_ledger_scaled_max'],
            trajectory_sha256=run['trajectory_sha256'],
            reason='513-node phonon representation does not meet the selected static population gate; no complete time-reference or three-grid dynamic assessment.'))
    gates=[
        dict(id='immutable_parent',status='PASS',detail='Exact hashes of R2 source/catalog and stage1 contract/certificate preserved.'),
        dict(id='complementary_energy_and_forces',status='PASS_SAMPLED_STATIC_DOMAIN',
            independent_maxima=energy['maxima'],moment_errors_by_refinement=moments['maxima_by_refinement'],
            detail='630/1260/2520 fixed count states; high-precision implicit spectral reference and independently refined moment quadrature. This is not uniform validation of arbitrary distributions.'),
        dict(id='KWT_and_synthetic_Debye',status=closure['status'],maxima=closure['maxima'],
            detail='Inherited mobility algebra and synthetic SI conversions; no microscopic mobility calibration or NbN admission.'),
        dict(id='event_equilibrium_faces_entropy',status='PASS_DECLARED_STATIC_CASES',
            maximum_relative_equilibrium_event_imbalance=max(e['maximum_relative_event_imbalance'] for row in algebra for e in row['equilibria']),
            maximum_weighted_equilibrium_RHS=max(e['weighted_rhs_L1'] for row in algebra for e in row['equilibria']),
            all_faces_inward=all(row['faces_inward'] for row in algebra),
            minimum_isolated_entropy_production=min(row['entropy_production'] for row in algebra),
            detail='Four static fields and three temperatures on the selected 1025-node phonon grid. Boundary trajectories are a separate gate.'),
        dict(id='transport_continuum',status=transport['selected_precision_status'],
            mesh_summary=transport['mesh_summary'],reference_budget=transport['reference_max_relative_RHS_refinement'],
            detail='Independent continuous-energy response and power on normal/gapped/gapless thermal and nonthermal profiles; not full coupled trajectories.'),
        dict(id='reaction_continuum_IR_001',status='PARTIAL_WITHIN_ERROR_BUDGET_CONVERGENCE_UNRESOLVED',
            maximum_full_RHS_relative_errors=[row['relative_L1'] for row in event['electronic_weak_maxima'] if row['representation']=='full_operator'],
            detail='All three-grid measured weak errors are below 1e-3, but a small nondecreasing BCS floor must be separated from eta and quadrature. Magnitude alone does not satisfy the monotonic convergence gate.'),
        dict(id='selected_phonon_interpolation',status='PASS_DECLARED_STATIC_CASES',phonon_nodes=1025,
            maxima=event['phonon_interpolation_maxima'][-1]['maxima'],
            detail='Same electronic representation, analytic populations evaluated at exact transition energy as the interpolation reference; not a full continuous-reaction or time check.'),
        dict(id='selected_collision_quadrature',status='PASS_INTERNAL_REFINEMENT',
            maximum_electronic_RHS_relative_L1=event['selected_inner_quadrature_maximum_relative_L1'],
            detail='630 electrons, 1025 phonons, IR=.005; inner orders2/3/4. Self-convergence does not replace an independent continuous reference.'),
        dict(id='selected_IR_0005_independent_reactions',status='PENDING',
            detail='Exact same-cut independent continuous powers and electronic/phonon weak responses have not been compared. IR omission bounds are separate evidence.'),
        dict(id='infrared_and_upper_support_static',status='PASS_DECLARED_ANALYTIC_PROFILES_ONLY',
            detail='Independent positive-cut omission and upper absorption bounds exist for nine analytic profiles. Actual selected coupled trajectories, external populations and arbitrary narrow distributions remain outside this result.'),
        dict(id='negative_controls',status=read('review/negative_controls.json')['status'],
            detail='Tests detect double condensate heat, omitted spectral work and a common incorrect collision factor.'),
        dict(id='smooth_RK4_suboperators',status=read('time_suboperator_order.json')['status'],
            detail='BGK and escape exact exponential problems show fourth order; this cannot certify the full nonsmooth event system.'),
        dict(id='simultaneous_trajectories_and_time_convergence',status='PENDING_COMPUTE',
            detail='No completed independent temporal reference and no accepted three-time-level comparisons for the selected candidate. A bounded reference probe is INCOMPLETE.'),
        dict(id='dynamic_electron_phonon_resolution',status='PENDING_COMPUTE',
            detail='Three-grid full-trajectory electron and phonon population convergence remains unmeasured. Initial-state projection and static kernel refinement are insufficient.'),
        dict(id='selected_boundary_equilibrium_and_activity',status='PENDING_COMPUTE',
            detail='Repeat complete equilibrium, empty/sparse boundary, escape/heating, resolved-amplitude and channel-transfer checks on the selected 1025-phonon candidate with accepted time accuracy.'),
        dict(id='actual_trajectory_fields_and_support',status='PENDING_COMPUTE',
            detail='At least ten actual points plus endpoints/extrema, energy-force identities, and upper-support bounds must refer to the final selected and time-converged trajectories.'),
        dict(id='final_regression_suite',status='PASS',exit_code=regression['exit_code'],
            runtime_seconds=regression['runtime_seconds'],source_hashes_unchanged=True,
            detail='509 tests passed on Geminga, pytest46.76s, wrapper47.09s. All recorded source hashes match before/after and the current workspace. This cannot close missing physical convergence.')]
    data=dict(schema='pysnspd.stage2.admission_checkpoint.v1',status='PENDING_COMPUTE',
        stage2_status='NOT_CLOSED',numerical_admission=False,production_promotion=False,
        recommendation='Finish the user-run temporal reference and remaining registered validations before considering closure or the spatial stage.',
        date='2026-09-21',criteria=artifact('acceptance_criteria.json'),
        amendment=artifact('review/complementary_grid_amendment.json'),
        immutable_inputs=immutable,audited_source_sha256=source_hashes,
        candidate=dict(electron_states=630,phonon_nodes=1025,infrared=.005,phonon_cutoff=4.,
            electron_mesh='Complementary causal fixed-count quadrature; R2 vacuum and parent catalog unchanged',
            eta=1e-8,reaction_layout='resolved_panels',reaction_inner_order=2,reaction_outer_order=2,
            face_order=2,scope='Synthetic one/two cells; proposed configuration, not yet admitted'),
        gates=gates,historical_completed_trajectories=completed,
        computation_handoff=read('reference_handoff.json'),
        material=dict(status='SYNTHETIC_ONLY_NBN_ABSOLUTE_NOT_ADMITTED',
            detail='Debye inputs are declared synthetic. The traceable NbN shape is conditional; unresolved absolute units/basis/density and material preprocessing defects are not inferred away.'),
        model_scope=dict(status='EXPERIMENTAL_CELL_ALGORITHM_ONLY',
            excludes=['absolute material rates or latency','spatial condensate amplitude/phase gradients','principal-symbol/core stability','reservoir boundary problem','Poisson/current/circuit coupling','full SNSPD transient','production replacement']),
        historical_failures_preserved=True,tolerances_relaxed=False,
        evidence=[artifact(name) for name in artifacts],
        generator_sha256=digest(Path(__file__)))
    (STAGE/'stage2_admission.json').write_text(json.dumps(data,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    print(json.dumps(dict(status=data['status'],stage2_status=data['stage2_status'],gates=len(gates),evidence=len(artifacts))))


if __name__=='__main__':main()
