"""Runner contracts using synthetic files and mocked work, never physical RHS."""
from pathlib import Path
import copy
import hashlib
import importlib.util
import json
import sys
from types import SimpleNamespace

import numpy as np
import pytest


ROOT = Path(__file__).resolve().parents[1]
RUNNER = ROOT / 'sandbox/stage3_spatial/run_static_batch.py'


@pytest.fixture(scope='module')
def runner():
    # The script intentionally uses a sibling ``from progress import Progress``.
    directory = str(RUNNER.parent)
    sys.path.insert(0, directory)
    try:
        spec = importlib.util.spec_from_file_location('stage3_runner_semantic_tests', RUNNER)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
    finally:
        sys.path.remove(directory)
    return module


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


CASE = dict(id='weak_phase_thermal', amplitude_base=.9,
            amplitude_modulation=0., phase_modulation=.1,
            q_bare_ell0=.1, occupation='thermal')
HASHES = {'synthetic_test_source.py': '1'*64}
PASS = 'PASS_STATIC_IDENTITIES_AND_SAMPLED_SYMBOL'


def artifact(directory, *, task='weak_phase_thermal_n8', case=None, cells=8,
             hashes=None, status=PASS):
    """Small synthetic NPZ envelopes test provenance, not physical populations."""
    directory.mkdir(parents=True, exist_ok=True)
    path = directory/(task+'.json')
    data = path.with_suffix('.npz')
    initial = directory/(task+'_initial.npz')
    np.savez(data, synthetic_fixture=np.array([1.,2.]))
    np.savez(initial, synthetic_fixture=np.array([3.,4.]))
    record = dict(schema='pysnspd.stage3.static-case.v1',
                  case=copy.deepcopy(CASE if case is None else case), cells=cells,
                  source_hashes=dict(HASHES if hashes is None else hashes), status=status,
                  arrays_sha256=sha(data), initial_sha256=sha(initial),
                  runtime_seconds=0., scope='NONPHYSICAL SEMANTIC TEST FIXTURE')
    path.write_text(json.dumps(record), encoding='utf-8')
    receipt = dict(json_sha256=sha(path), arrays_sha256=sha(data),
                   initial_sha256=sha(initial), source_hashes=record['source_hashes'])
    path.with_suffix('.receipt.json').write_text(json.dumps(receipt), encoding='utf-8')
    return path, record


def test_compatible_reuse_accepts_exact_case_grid_sources_and_artifacts(runner, tmp_path):
    _, expected = artifact(tmp_path)
    assert runner.verify_reuse(tmp_path, 'weak_phase_thermal_n8', HASHES, CASE, 8) == expected


def test_absent_case_is_not_claimed_as_reused(runner, tmp_path):
    assert runner.verify_reuse(tmp_path, 'weak_phase_thermal_n8', HASHES, CASE, 8) is None


@pytest.mark.parametrize('mismatch', ['case_id', 'same_id_different_parameters', 'cells'])
def test_renamed_or_different_case_is_rejected_even_with_valid_receipt(runner, tmp_path, mismatch):
    wrong_case = copy.deepcopy(CASE)
    wrong_cells = 8
    if mismatch == 'case_id':
        wrong_case['id'] = 'weak_phase_nonthermal'
    elif mismatch == 'same_id_different_parameters':
        wrong_case['q_bare_ell0'] = .11
    else:
        wrong_cells = 16
    # All hashes are valid: filename alone must not establish case identity.
    artifact(tmp_path, case=wrong_case, cells=wrong_cells)
    with pytest.raises(ValueError):
        runner.verify_reuse(tmp_path, 'weak_phase_thermal_n8', HASHES, CASE, 8)


@pytest.mark.parametrize('suffix', ['.npz', '_initial.npz', '.json'])
def test_tampered_arrays_initial_or_result_are_rejected(runner, tmp_path, suffix):
    path, _ = artifact(tmp_path)
    target = tmp_path/(path.stem+suffix)
    if suffix == '.json':
        record = json.loads(target.read_text())
        record['runtime_seconds'] = 10.
        target.write_text(json.dumps(record), encoding='utf-8')
    else:
        # Appending preserves an otherwise readable archive but changes its hash.
        target.write_bytes(target.read_bytes()+b'tampered')
    with pytest.raises(ValueError):
        runner.verify_reuse(tmp_path, path.stem, HASHES, CASE, 8)


def test_source_mismatch_cannot_reuse_old_result(runner, tmp_path):
    artifact(tmp_path)
    with pytest.raises(ValueError):
        runner.verify_reuse(tmp_path, 'weak_phase_thermal_n8',
                            {'synthetic_test_source.py': '2'*64}, CASE, 8)


@pytest.mark.parametrize('status', ['FAIL_STATIC_CHECK', 'RUNNING'])
def test_failed_or_incomplete_case_is_not_reusable(runner, tmp_path, status):
    artifact(tmp_path, status=status)
    with pytest.raises(ValueError):
        runner.verify_reuse(tmp_path, 'weak_phase_thermal_n8', HASHES, CASE, 8)


def test_missing_receipt_is_not_a_completed_reuse(runner, tmp_path):
    path, _ = artifact(tmp_path)
    path.with_suffix('.receipt.json').unlink()
    with pytest.raises((ValueError, OSError)):
        runner.verify_reuse(tmp_path, path.stem, HASHES, CASE, 8)


def registration(case_count=1):
    cases = [dict(id='synthetic_response_a')]
    if case_count == 2:
        cases.append(dict(id='synthetic_response_b'))
    return dict(cases=cases, cell_counts=[8,16,32], thresholds=dict(
        spatial_response_absolute_floor=1e-8,
        nonlinear_mesh_response_target_relative=.01,
        gradient_finest_relative_error=.01))


def synthetic_response(case, cells):
    value = 1+1/cells**2
    vector = dict(real=[value,0.,0.,0.,0.], imag=[0.]*5)
    return dict(case=copy.deepcopy(case), cells=cells, status=PASS,
                excess_energy_density_bar=value, induced_force_modes=vector,
                induced_current_modes=copy.deepcopy(vector), analytic_gradient=None)


def complete_responses(reg):
    return [synthetic_response(case, cells)
            for case in reg['cases'] for cells in reg['cell_counts']]


def test_complete_three_grid_fixture_is_assessed(runner):
    reg = registration()
    rows = runner.assess_grids(complete_responses(reg), reg, require_complete=True)
    assert len(rows) == 1
    assert rows[0]['status'] == 'PASS_REGISTERED_SPATIAL_RESPONSE'


def set_current_response(rows,values):
    for row,value in zip(rows,values):
        row['induced_current_modes']=dict(real=[value,0.,0.,0.,0.],imag=[0.]*5)


@pytest.mark.parametrize('values', [
    [5.4e-9,5.1e-9,5.025e-9],  # Response itself is below the registered floor.
    [1.+1.5625e-8,1.+3.90625e-9,1.+9.765625e-10],  # Only differences are unresolved.
    [1.,1.,1.],  # Exact equality of a nonzero result is no convergence certificate.
])
def test_unresolved_response_has_no_relative_pass(runner,values):
    reg=registration();data=complete_responses(reg)
    set_current_response(data,values)
    row=runner.assess_grids(data,reg,require_complete=True)[0]
    metric=row['metrics']['induced_current_modes']
    assert metric['status']=='ROUND_OFF_LIMITED'
    assert metric['passed'] is None
    assert metric['relative_certificate'] is False
    assert metric['relative_estimate'] is None
    assert row['status']=='SPATIAL_RESPONSE_SCOPE_LIMITED'
    assert 'induced_current_modes' in row['unresolved_observables']
    assert 'induced_current_modes' not in row['relative_certificate_observables']


def test_nondecreasing_differences_are_inconclusive_not_a_relative_pass(runner):
    reg=registration();data=complete_responses(reg)
    set_current_response(data,[1.,1.1,1.3])
    row=runner.assess_grids(data,reg,require_complete=True)[0]
    metric=row['metrics']['induced_current_modes']
    assert metric['status']=='INCONCLUSIVE_SPATIAL_ESTIMATE'
    assert metric['passed'] is None
    assert metric['richardson_conservative_estimate'] is None
    assert row['status']=='SPATIAL_RESPONSE_SCOPE_LIMITED'


def test_resolved_difference_exceeding_registered_budget_retains_failure(runner):
    reg=registration();data=complete_responses(reg)
    set_current_response(data,[2.,1.5,1.2])
    row=runner.assess_grids(data,reg,require_complete=True)[0]
    metric=row['metrics']['induced_current_modes']
    assert metric['status']=='FAIL_PRECISION_TARGET'
    assert metric['passed'] is False
    assert row['status']=='PRECISION_TARGET_NOT_MET'


@pytest.mark.parametrize('values,expected', [
    ([0.,0.,0.],'PASS_ANALYTIC_ZERO_ABSOLUTE'),
    ([2e-8,0.,0.],'FAIL_ANALYTIC_ZERO_ABSOLUTE'),
    ([2e-8,2e-8,2e-8],'FAIL_ANALYTIC_ZERO_ABSOLUTE'),
])
def test_zero_flow_real_field_current_has_only_an_absolute_control(runner,values,expected):
    reg=registration()
    reg['cases'][0].update(q_bare_ell0=0.,phase_modulation=0.)
    data=complete_responses(reg);set_current_response(data,values)
    row=runner.assess_grids(data,reg,require_complete=True)[0]
    metric=row['metrics']['induced_current_modes']
    assert metric['status']==expected
    assert metric['analytic_zero'] is True
    assert metric['relative_certificate'] is False
    assert metric['relative_estimate'] is None
    assert metric['absolute_control_passed']==(expected=='PASS_ANALYTIC_ZERO_ABSOLUTE')
    assert row['status']==('PASS_REGISTERED_SPATIAL_RESPONSE' if metric['absolute_control_passed']
                           else 'PRECISION_TARGET_NOT_MET')


def test_spectral_record_exposes_changed_phase_stiffness_with_unchanged_minimum(runner):
    # Pure synthetic arrays: no catalogue query, spectrum, or spatial evaluation.
    coarse=dict(eigenvalues=[1.5,2.],matrix=[[1.5,0.],[0.,2.]],uncertainty=1e-6)
    fine=SimpleNamespace(eigenvalues=np.array([1.5,2.3]),
                         matrix=np.diag([1.5,2.3]),uncertainty=2e-6)
    record=runner.spectral_probe_record(4,coarse,fine,[1.,2.,3.],[1.01,1.98,3.03])
    assert record['sampled_spectral_shift']==0.
    assert record['coarse_eigenvalues']==[1.5,2.]
    assert record['fine_eigenvalues']==[1.5,2.3]
    assert record['eigenvalue_absolute_differences']==pytest.approx([0.,.3])
    assert record['matrix_operator_norm_difference']==pytest.approx(.3)
    assert record['fine_matrix'][1][1]==2.3
    assert record['electronic_moments']['names']==['u','u_a','u_Gamma']
    assert record['electronic_moments']['absolute_differences']==pytest.approx([.01,.02,.03])
    assert record['coarse_uncertainty']==1e-6
    assert record['fine_uncertainty']==2e-6
    assert record['remaining_positive_margin']==pytest.approx(1.5-3e-6)
    json.dumps(record,allow_nan=False)


@pytest.mark.parametrize('missing', ['all', 'one_grid', 'one_case'])
def test_full_campaign_rejects_incomplete_coverage(runner, missing):
    reg = registration(case_count=2)
    rows = complete_responses(reg)
    if missing == 'all':
        rows = []
    elif missing == 'one_grid':
        rows.pop()
    else:
        rows = [row for row in rows if row['case']['id'] == reg['cases'][0]['id']]
    with pytest.raises(ValueError):
        runner.assess_grids(rows, reg, require_complete=True)


def test_full_campaign_rejects_duplicate_case_grid(runner):
    reg = registration()
    rows = complete_responses(reg)
    rows.append(copy.deepcopy(rows[0]))
    with pytest.raises(ValueError):
        runner.assess_grids(rows, reg, require_complete=True)


def test_partial_pilot_does_not_claim_three_grid_precision(runner):
    reg = registration()
    rows = runner.assess_grids([synthetic_response(reg['cases'][0],8)],
                               reg, require_complete=False)
    assert not any(row.get('status') == 'PASS_REGISTERED_SPATIAL_RESPONSE' for row in rows)


def test_pilot_main_preserves_limited_scope_without_running_physics(runner, tmp_path, monkeypatch):
    # Reuse the declaration as text only. Every computational entry is replaced.
    source = ROOT/'docs/implementation/stage3/iteration_20260923/registration.json'
    declared = json.loads(source.read_text(encoding='utf-8'))
    declared_path = tmp_path/'registration.json'
    declared_path.write_text(json.dumps(declared), encoding='utf-8')
    out = tmp_path/'pilot_output'
    calls = []

    def forbidden_physics(*args, **kwargs):
        raise AssertionError('semantic runner test must not construct a physical model')

    def fake_case(catalog, case, cells, registration, output, source_hashes, log_callback):
        calls.append((case['id'], cells))
        _, record = artifact(output.parent, task=output.stem, case=case,
                             cells=cells, hashes=source_hashes)
        return record

    monkeypatch.setattr(runner, 'PeriodicSpatialFunctional', forbidden_physics)
    monkeypatch.setattr(runner, 'OccupationEnergyCatalog', SimpleNamespace(load=lambda _: object()))
    monkeypatch.setattr(runner, 'refined_count_catalog', lambda catalog: catalog)
    monkeypatch.setattr(runner, 'frozen_sources', lambda _: dict(HASHES))
    monkeypatch.setattr(runner, 'negative_control', lambda _: dict(status='PASS_EXPECTED_REJECTION'))
    monkeypatch.setattr(runner, 'run_case', fake_case)
    monkeypatch.setattr(sys, 'argv', [str(RUNNER), '--pilot', '--registration', str(declared_path),
                                    '--output-root', str(out)])
    runner.main()
    summary = json.loads((out/'summary.json').read_text(encoding='utf-8'))
    assert calls == [('weak_phase_thermal',8)]
    assert summary['status'] == 'PILOT_PASS_FULL_CAMPAIGN_PENDING'
    assert summary['completed_cases'] == 1
    assert summary['stage3_closed'] is False
    assert summary['production_promotion'] is False
    assert summary['new_trajectories'] == 0
    assert summary['circuit_and_physical_boundaries_implemented'] is False
    assert not any(r.get('status') == 'PASS_REGISTERED_SPATIAL_RESPONSE'
                   for r in summary['grid_comparisons'])


def test_completed_scope_limited_campaign_is_not_promoted_to_resolved_precision(runner,tmp_path,monkeypatch):
    source=ROOT/'docs/implementation/stage3/iteration_20260923/registration.json'
    declared=json.loads(source.read_text(encoding='utf-8'))
    declared_path=tmp_path/'registration.json'
    declared_path.write_text(json.dumps(declared),encoding='utf-8')
    out=tmp_path/'campaign_output'

    def forbidden_physics(*args,**kwargs):
        raise AssertionError('semantic test must not construct a physical model')

    def fake_case(catalog,case,cells,registration,output,source_hashes,log_callback):
        return artifact(output.parent,task=output.stem,case=case,cells=cells,
                        hashes=source_hashes)[1]

    def fake_assess(results,registration,require_complete=False):
        assert require_complete and len(results)==18
        return [dict(case='synthetic_scope_only',status='SPATIAL_RESPONSE_SCOPE_LIMITED',
                     unresolved_observables=['induced_current_modes'])]

    monkeypatch.setattr(runner,'PeriodicSpatialFunctional',forbidden_physics)
    monkeypatch.setattr(runner,'OccupationEnergyCatalog',SimpleNamespace(load=lambda _:object()))
    monkeypatch.setattr(runner,'refined_count_catalog',lambda catalog:catalog)
    monkeypatch.setattr(runner,'frozen_sources',lambda _:dict(HASHES))
    monkeypatch.setattr(runner,'negative_control',lambda _:dict(status='PASS_EXPECTED_REJECTION'))
    monkeypatch.setattr(runner,'run_case',fake_case)
    monkeypatch.setattr(runner,'assess_grids',fake_assess)
    monkeypatch.setattr(sys,'argv',[str(RUNNER),'--execute','--registration',str(declared_path),
                                  '--output-root',str(out)])
    runner.main()
    summary=json.loads((out/'summary.json').read_text(encoding='utf-8'))
    assert summary['status']=='STATIC_CAMPAIGN_COMPLETE_SCOPE_LIMITED'
    assert summary['completed_cases']==18
    assert summary['spatial_precision_fully_resolved'] is False
    assert summary['stage3_closed'] is False
    assert summary['production_promotion'] is False
