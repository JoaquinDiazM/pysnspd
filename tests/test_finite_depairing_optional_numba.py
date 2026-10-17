"""Optional Numba: pure-Python oracle parity without production edits."""
import builtins
import importlib.util
import json
import os
from pathlib import Path

import numpy as np
import pytest

ROOT=Path(os.environ.get('PYSNSPD_TEST_REPO',Path(__file__).resolve().parents[1]))
TARGET=Path(os.environ.get('PYSNSPD_FINITE_DEPAIRING_TEST_SOURCE',
                          ROOT/'pysnspd/experimental/finite_depairing_sums.py'))
REFERENCE=json.loads((ROOT/'tests/fixtures/finite_depairing_reference.json').read_text())


def without_numba(monkeypatch):
    original=builtins.__import__
    def blocked(name,*args,**kwargs):
        if name=='numba' or name.startswith('numba.'):
            raise ModuleNotFoundError('Numba intentionally unavailable in optional-dependency test')
        return original(name,*args,**kwargs)
    spec=importlib.util.spec_from_file_location('_finite_depairing_without_numba',TARGET)
    module=importlib.util.module_from_spec(spec)
    with monkeypatch.context() as context:
        context.setattr(builtins,'__import__',blocked)
        spec.loader.exec_module(module)
    assert module.NUMBA_AVAILABLE is False
    assert module.prange is range
    return module


def fixture_inputs():
    cases=REFERENCE['cases']
    return tuple(np.array([case[key] for case in cases]) for key in ('rho','theta','Q'))


def test_import_without_numba_matches_all_ten_independent_reference_partials(monkeypatch):
    module=without_numba(monkeypatch)
    actual=np.array(module.finite_depairing_partials(*fixture_inputs()))
    expected=np.array([case['values'] for case in REFERENCE['cases']]).T
    assert actual.shape[0]==10
    np.testing.assert_allclose(actual,expected,rtol=0.,atol=1e-10)


def test_without_numba_matches_compiled_oracle(monkeypatch):
    pytest.importorskip('numba')
    from pysnspd.experimental.finite_depairing_sums import finite_depairing_partials
    module=without_numba(monkeypatch)
    inputs=fixture_inputs()
    actual=np.array(module.finite_depairing_partials(*inputs))
    compiled=np.array(finite_depairing_partials(*inputs))
    np.testing.assert_allclose(actual,compiled,rtol=2e-12,atol=1e-12)


def test_fallback_preserves_zero_Q_identity_and_declared_resolution_failure(monkeypatch):
    module=without_numba(monkeypatch)
    actual=np.array(module.finite_depairing_partials(np.array([0.,.1,1.,4.]),.1,0.))
    np.testing.assert_array_equal(actual[[0,1,3,4,7,9]],0.)
    assert np.all(actual[5]>0)
    with pytest.raises(ArithmeticError,match='resolution exceeded'):
        module.finite_depairing_partials(1.,1e-10,1.)
    with pytest.raises(ValueError):
        module.finite_depairing_partials(1.,1.,-1.)