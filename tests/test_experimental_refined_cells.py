"""Complementary occupation-grid contracts; the archived catalogue is immutable."""
from pathlib import Path
import hashlib
import numpy as np
import pytest
from scipy.integrate import quad
from pysnspd.experimental.energy_catalog import OccupationEnergyCatalog
from pysnspd.experimental.refined_cells import ComplementaryCountCatalog,refined_count_catalog


@pytest.fixture(scope='module')
def base():
    path=Path(__file__).resolve().parents[1]/'docs/implementation/stage1_r2/catalogs/occupation_catalog.npz'
    assert hashlib.sha256(path.read_bytes()).hexdigest()=='86dd59209220dbe907b9da3608818d5cd462d9df4f97ef39ee6ef681b68ce970'
    return OccupationEnergyCatalog.load(path)


def test_new_grid_changes_degrees_of_freedom_without_mutating_parent(base):
    original=base.count_nodes.copy()
    coarse=refined_count_catalog(base,1);fine=refined_count_catalog(base,2)
    assert len(fine.count_nodes)>len(coarse.count_nodes)>len(original)
    np.testing.assert_array_equal(base.count_nodes,original)
    assert coarse.vacuum is base.vacuum and coarse.eta==base.eta
    for mesh in (coarse,fine):
        assert not mesh.count_nodes.flags.writeable and not mesh.count_weights.flags.writeable
        assert np.sum(mesh.count_weights)==pytest.approx(12.,rel=1e-13)
        assert np.dot(mesh.count_weights,mesh.count_nodes)==pytest.approx(72.,rel=1e-13)


@pytest.mark.parametrize('x,w',[
    ([0.,1.],[1.,1.]),([1.,1.],[1.,1.]),([2.,1.],[1.,1.]),
    ([1.,np.nan],[1.,1.]),([1.,2.],[1.,np.inf]),([1.,2.],[0.,1.]),
    ([1.,2.],[1.]),([[1.,2.]],[[1.,1.]])])
def test_invalid_count_coordinates_rejected(base,x,w):
    with pytest.raises(ValueError):ComplementaryCountCatalog(base,np.array(x),np.array(w))


@pytest.mark.parametrize('bad',[-1e-16,1+1e-15,np.nan,np.inf])
def test_invalid_pauli_occupations_rejected(base,bad):
    mesh=refined_count_catalog(base)
    p=np.zeros(len(mesh.count_nodes));p[len(p)//2]=bad
    with pytest.raises(ValueError):mesh.evaluate(.72,.2,p)


@pytest.mark.parametrize('amplitude,gamma',[(.6,0.),(.72,.2),(.72,.72),(.35,.65)])
def test_force_derivatives_belong_to_the_same_potential(base,amplitude,gamma):
    # Keep populations fixed while changing fields. Differentiating a freshly
    # thermalized p at each field would check a different thermodynamic function.
    mesh=refined_count_catalog(base)
    p=.13*np.exp(-mesh.count_nodes/.3)+.08*np.exp(-((mesh.count_nodes-1.2)/.3)**2)
    energy,fa,fg=mesh.evaluate(amplitude,gamma,p)
    kernel=mesh.energy_kernel(amplitude,gamma)
    expected=np.asarray(mesh.vacuum.evaluate(amplitude,gamma))+4*np.asarray([np.dot(mesh.count_weights*v,p) for v in kernel])
    np.testing.assert_allclose([energy,fa,fg],expected,rtol=1e-13,atol=1e-13)
    h=2e-5
    def derivative(function,point):
        return (function(point-2*h)-8*function(point-h)+8*function(point+h)-function(point+2*h))/(12*h)
    da=derivative(lambda a:mesh.evaluate(a,gamma,p)[0],amplitude)
    assert abs(da-fa)/max(1.,abs(fa))<1e-7
    # The independent60-digit audit covers Gamma0 and tinyGamma without a
    # finite difference crossing the domain or rounding away its increments.
    if gamma>0:
        dg=derivative(lambda g:mesh.evaluate(amplitude,g,p)[0],gamma)
        assert abs(dg-fg)/max(1.,abs(fg))<1e-7


def test_normal_limit_and_low_amplitude_domain_are_explicit(base):
    mesh=refined_count_catalog(base)
    e,da,dg=mesh.energy_kernel(0.,.2)
    np.testing.assert_array_equal(e,mesh.count_nodes)
    np.testing.assert_array_equal(da,np.zeros_like(e))
    np.testing.assert_array_equal(dg,np.zeros_like(e))
    with pytest.raises(ValueError):mesh.energy_kernel(.04,.2)


def test_Gamma_moment_resolves_its_regulator_layer_and_refines(base):
    amplitude=.72;eta=base.eta
    # Independent analytic derivative of the finite-eta BCS count inversion.
    # Its narrow sqrt(amplitude*eta) layer is invisible to an energy-only test.
    def integrand(x):
        derivative=-amplitude**2*eta*x/((x*x+eta*eta)**2+amplitude**2*eta*eta)
        return 4*.2*np.exp(-x/.25)*derivative
    layer=np.sqrt(amplitude*eta)
    breaks=[0.,eta,layer,3*layer,10*layer,.01,.1,1.,12.]
    expected=np.pi*amplitude/2+sum(quad(integrand,lo,hi,epsabs=1e-12,epsrel=1e-12)[0]
                                  for lo,hi in zip(breaks[:-1],breaks[1:]))
    errors=[]
    for level in (1,2,4):
        mesh=refined_count_catalog(base,level)
        measured=mesh.evaluate(amplitude,0.,.2*np.exp(-mesh.count_nodes/.25))[2]
        errors.append(abs(measured-expected)/abs(expected))
    assert max(errors)<1e-3
    assert errors[2]<errors[1]<errors[0]
    unresolved=refined_count_catalog(base,edge_order=2)
    wrong=unresolved.evaluate(amplitude,0.,.2*np.exp(-unresolved.count_nodes/.25))[2]
    assert abs(wrong-expected)/abs(expected)>1e-3


@pytest.mark.parametrize('options',[{'edge_order':1},{'edge_order':True},{'order':2.0},
                                    {'order':0},{'refinement':True}])
def test_quadrature_orders_and_levels_are_explicit_integers(base,options):
    with pytest.raises(ValueError):refined_count_catalog(base,**options)
