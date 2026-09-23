"""Nodal runner/oracle semantics with synthetic polynomials and no catalogue RHS."""
from pathlib import Path
import copy
import importlib.util
import json
import sys
from types import SimpleNamespace

import numpy as np
import pytest


ROOT=Path(__file__).resolve().parents[1]
RUNNER=ROOT/'sandbox/stage3_spatial/run_nodal_batch.py'


@pytest.fixture(scope='module')
def runner():
    directory=str(RUNNER.parent)
    sys.path.insert(0,directory)
    try:
        spec=importlib.util.spec_from_file_location('nodal_runner_semantic_tests',RUNNER)
        module=importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
    finally:
        sys.path.remove(directory)
    return module


class PolynomialCatalogue:
    """Synthetic smooth u(a,Gamma,p), deliberately unrelated to a material."""
    def __init__(self):
        self.calls=[]
        self.vacuum=SimpleNamespace(evaluate=lambda a,g:self.evaluate(a,g,np.zeros(2)))

    def evaluate(self,a,g,p):
        self.calls.append((float(a),float(g),np.asarray(p).copy()))
        occupation=float(np.sum(p))
        u=a**4+.3*g*g+.7*a*g+occupation*(.2*a+.1*g)
        ua=4*a**3+.7*g+.2*occupation
        ug=.6*g+.7*a+.1*occupation
        return u,ua,ug


def polynomial_energy(a,q,alpha,p,*,cells,h,ratio,kappa):
    """Independent covariant stencil sum, not the oracle's trigonometric formula."""
    x=np.arange(cells)*h
    z=a*np.exp(1j*q*x)
    links=np.full(cells,np.exp(1j*alpha),complex)
    links[-1]*=np.exp(1j*q*cells*h)
    shift=lambda v:links*np.roll(v,-1)
    adjoint=lambda v:np.roll(np.conj(links)*v,1)
    derivative=(8*(shift(z)-adjoint(z))-(shift(shift(z))-adjoint(adjoint(z))))/(12*h)
    laplacian=(30*z-16*(shift(z)+adjoint(z))+shift(shift(z))+adjoint(adjoint(z)))/(12*h*h)
    amplitude=abs(z)
    qdelta=np.imag(np.conj(z)*derivative)/(amplitude**2+.01)
    gamma=qdelta*qdelta/ratio
    occupation=float(np.sum(p))
    u=amplitude**4+.3*gamma**2+.7*amplitude*gamma+occupation*(.2*amplitude+.1*gamma)
    return float(h*np.sum(u-kappa*amplitude**2*qdelta**2)+kappa*h*np.real(np.vdot(z,laplacian)))


def fd5(function,h=1e-4):
    return (function(-2*h)-8*function(-h)+8*function(h)-function(2*h))/(12*h)


@pytest.mark.parametrize('population', [np.zeros(2),np.array([.2,.1])])
def test_uniform_oracle_matches_synthetic_stencil_energy_and_derivatives(runner,population):
    catalog=PolynomialCatalogue()
    cells=16;h=.4;a=.9;q=.21
    model=SimpleNamespace(h_bar=h,length_bar=cells*h,gap_ratio=1.764,kappa=np.pi/4)
    case=dict(amplitude_base=a,q_bare_ell0=q)
    oracle=runner.uniform_oracle(model,catalog,case,population)
    energy=lambda a_,alpha:polynomial_energy(a_,q,alpha,population,cells=cells,h=h,
                                           ratio=model.gap_ratio,kappa=model.kappa)
    assert oracle['energy_bar']==pytest.approx(energy(a,0),rel=2e-12,abs=1e-12)
    assert oracle['force_radial_density_bar']==pytest.approx(
        fd5(lambda da:energy(a+da,0))/model.length_bar,rel=1e-8,abs=1e-9)
    # All link phases vary together: dU/dalpha=N*one_link_current.
    assert oracle['link_current_bar']==pytest.approx(
        fd5(lambda alpha:energy(a,alpha))/cells,rel=1e-8,abs=1e-9)
    assert all(entry[0]==a for entry in catalog.calls)
    assert all(np.array_equal(entry[2],population) for entry in catalog.calls)
    assert oracle['negative_laplacian_symbol_K']!=oracle['first_derivative_symbol_Q']**2


def test_uniform_zero_phase_current_is_exact_zero(runner):
    model=SimpleNamespace(h_bar=.5,length_bar=4.,gap_ratio=1.764,kappa=np.pi/4)
    result=runner.uniform_oracle(model,PolynomialCatalogue(),
        dict(amplitude_base=.9,q_bare_ell0=0.),np.zeros(2))
    assert result['link_current_bar']==0
    assert result['first_derivative_symbol_Q']==0
    assert result['negative_laplacian_symbol_K']==0


def fd_registration():
    return dict(fd_cells=[16],declared_spatial_order=4,richardson_denominator=15)


@pytest.mark.parametrize('cells,pilot,performed',[(8,False,False),(16,False,True),(32,False,False),(8,True,True)])
def test_fd_scope_is_explicit_and_does_not_invent_passes(runner,cells,pilot,performed):
    reg=fd_registration()
    checks={name:dict(status='PASS',passed=True) for name in ('force','current')} if performed else {}
    meta=runner.derivative_metadata(cells,reg,checks,pilot=pilot)
    assert meta['performed'] is performed
    if performed:
        assert meta['status']=='PASS' and set(checks)<=set(meta)
    else:
        assert meta['status']=='NOT_REPEATED_UNIT_AND_N16_REFERENCES'
        assert 'force' not in meta and 'current' not in meta
    record=dict(status='PASS_STATIC_IDENTITIES_AND_SAMPLED_SYMBOL' if performed else
                'PASS_STATIC_STRUCTURE_WITHOUT_REPEATED_FD',derivative_checks=meta)
    assert runner.accepted_case(record)


def test_skipped_metadata_cannot_be_presented_as_fd_pass(runner):
    meta=runner.derivative_metadata(32,fd_registration(),{})
    assert not runner.accepted_case(dict(status='PASS_STATIC_IDENTITIES_AND_SAMPLED_SYMBOL',
                                        derivative_checks=meta))
    with pytest.raises(ValueError):
        runner.derivative_metadata(16,fd_registration(),{})
    with pytest.raises(ValueError):
        runner.derivative_metadata(8,fd_registration(),{'force':dict(status='PASS',passed=True)})


def test_cost_weights_reflect_selective_fd_and_verified_reuse(runner):
    case=dict(occupation='thermal');reg=fd_registration()
    assert runner.cost_weight(case,8,reg)==8*8
    assert runner.cost_weight(case,16,reg)==16*32
    assert runner.cost_weight(case,32,reg)==32*8
    assert runner.cost_weight(case,8,reg,pilot=True)==8*32
    assert runner.cost_weight(case,8,reg,reused={'derivative_checks':{'performed':True}})==8*32
    assert runner.cost_weight(dict(occupation='vacuum'),16,reg)==pytest.approx(16*32*.005)


def test_fourth_order_richardson_denominator_and_coverage(runner):
    case=dict(id='synthetic_nonuniform',q_bare_ell0=.1,phase_modulation=.1)
    reg=dict(cases=[case],cell_counts=[8,16,32],richardson_denominator=15,
        richardson_denominator_by_observable={'induced_current_modes':3},
        thresholds=dict(spatial_response_absolute_floor=1e-12,nonlinear_mesh_response_target_relative=.01,
                        gradient_finest_relative_error=.01))
    values=[]
    for n in reg['cell_counts']:
        value=1+256/n**4
        vector=dict(real=[value,0.,0.,0.,0.],imag=[0.]*5)
        values.append(dict(case=case,cells=n,excess_energy_density_bar=value,
            induced_force_modes=vector,induced_current_modes=copy.deepcopy(vector),analytic_gradient=None))
    rows=runner.assess_grids(values,reg,require_complete=True)
    metric=rows[0]['metrics']['excess_energy_density_bar']
    assert metric['observed_order']==pytest.approx(4.)
    assert metric['richardson_conservative_estimate']==pytest.approx(metric['medium_fine_difference']/15)
    current=rows[0]['metrics']['induced_current_modes']
    assert metric['richardson_denominator']==15
    assert current['richardson_denominator']==3
    assert current['richardson_conservative_estimate']==pytest.approx(current['medium_fine_difference']/3)
    for invalid in (values[:-1],values+[copy.deepcopy(values[0])]):
        with pytest.raises(ValueError):runner.assess_grids(invalid,reg,require_complete=True)
    assert runner.assess_grids(values[:1],reg,require_complete=False)==[]


def test_help_does_not_load_a_catalogue_or_run_a_case(runner,monkeypatch,capsys):
    def forbidden(*args,**kwargs):raise AssertionError('help must not execute physics')
    monkeypatch.setattr(runner,'OccupationEnergyCatalog',SimpleNamespace(load=forbidden))
    monkeypatch.setattr(runner,'run_case',forbidden)
    monkeypatch.setattr(sys,'argv',[str(RUNNER),'--help'])
    with pytest.raises(SystemExit) as exit_info:runner.main()
    assert exit_info.value.code==0
    text=capsys.readouterr().out
    assert '--pilot' in text and '--execute' in text and '--reuse' in text


@pytest.mark.parametrize('pilot',[False,True])
def test_mocked_campaign_records_performed_and_skipped_scopes(runner,tmp_path,monkeypatch,pilot):
    declaration=ROOT/'docs/implementation/stage3/nodal_20260923/registration.json'
    reg=json.loads(declaration.read_text(encoding='utf-8'))
    out=tmp_path/'synthetic_output'
    hashes={'synthetic': '1'*64}
    calls=[]
    def forbidden(*args,**kwargs):raise AssertionError('mocked campaign must not construct a physical model')
    def fake_case(catalog,case,cells,registration,output,source_hashes,callback,*,pilot=False):
        performed=runner.fd_required(cells,registration,pilot=pilot)
        checks={key:dict(status='PASS',passed=True) for key in ('force','current')} if performed else {}
        metadata=runner.derivative_metadata(cells,registration,checks,pilot=pilot)
        record=dict(schema='pysnspd.stage3.nodal-static-case.v1',case=case,cells=cells,
            derivative_checks=metadata,runtime_seconds=.001,source_hashes=source_hashes,
            status='PASS_STATIC_IDENTITIES_AND_SAMPLED_SYMBOL' if performed else
                   'PASS_STATIC_STRUCTURE_WITHOUT_REPEATED_FD')
        np.savez(output.with_suffix('.npz'),synthetic=np.array([1.]))
        initial=output.with_name(output.stem+'_initial.npz')
        np.savez(initial,synthetic=np.array([2.]))
        record.update(arrays_sha256=runner.digest(output.with_suffix('.npz')),initial_sha256=runner.digest(initial))
        runner.write_json(output,record)
        calls.append((case['id'],cells,performed))
        return record
    def fake_assess(results,registration,require_complete=False):
        assert require_complete is (not pilot)
        assert len(results)==(1 if pilot else 18)
        return [] if pilot else [dict(status='SPATIAL_RESPONSE_SCOPE_LIMITED')]
    monkeypatch.setattr(runner,'PeriodicNodalSpatialFunctional',forbidden)
    monkeypatch.setattr(runner,'OccupationEnergyCatalog',SimpleNamespace(load=lambda _:object()))
    monkeypatch.setattr(runner,'refined_count_catalog',lambda value:value)
    monkeypatch.setattr(runner,'negative_control',lambda _:dict(status='PASS_EXPECTED_REJECTION'))
    monkeypatch.setattr(runner,'frozen_sources',lambda _:hashes)
    monkeypatch.setattr(runner,'run_case',fake_case)
    monkeypatch.setattr(runner,'assess_grids',fake_assess)
    monkeypatch.setattr(sys,'argv',[str(RUNNER),'--pilot' if pilot else '--execute',
        '--registration',str(declaration),'--output-root',str(out)])
    runner.main()
    result=json.loads((out/'summary.json').read_text(encoding='utf-8'))
    assert result['stage3_closed'] is False and result['production_promotion'] is False
    assert result['spatial_precision_fully_resolved'] is False
    assert len(result['cases_with_directional_fd'])==(1 if pilot else 6)
    assert len(result['cases_without_repeated_fd'])==(0 if pilot else 12)
    assert result['status']==('PILOT_PASS_FULL_CAMPAIGN_PENDING' if pilot else 'STATIC_CAMPAIGN_COMPLETE_SCOPE_LIMITED')
    assert all(n==16 or pilot for _,n,performed in calls if performed)
