"""NPZ attachment tests for PRE Usadel supercurrent tables."""
from __future__ import annotations

from types import SimpleNamespace

import numpy as np

from pysnspd.gtdgl.usadel_current import (
    UsadelCatalogWithSupercurrentTable,
    attach_usadel_supercurrent_table_from_npz,
    load_usadel_supercurrent_table_arrays_npz,
)


def test_load_usadel_supercurrent_table_arrays_npz(tmp_path):
    path = tmp_path / "catalog.npz"
    np.savez(
        path,
        js_A_m2=np.arange(6.0).reshape(2, 3),
        q_axis_m_inv=np.array([0.0, 1.0, 2.0]),
        delta_axis_J=np.array([1.0, 2.0]),
        object_metadata=np.array(["ignored"], dtype=object),
    )
    arrays = load_usadel_supercurrent_table_arrays_npz(path)
    assert set(arrays) >= {"js_A_m2", "q_axis_m_inv", "delta_axis_J"}
    assert arrays["js_A_m2"].shape == (2, 3)


def test_attach_usadel_supercurrent_table_from_npz_wraps_base_catalog(tmp_path):
    path = tmp_path / "catalog.npz"
    np.savez(
        path,
        j_s_A_m2=np.arange(4.0),
        q_values_m_inv=np.arange(4.0),
    )
    base = SimpleNamespace(files=["rho"], rho=np.ones(3))
    wrapped = attach_usadel_supercurrent_table_from_npz(base, path)
    assert isinstance(wrapped, UsadelCatalogWithSupercurrentTable)
    assert "js_A_m2" in wrapped.files
    assert np.array_equal(wrapped.rho, base.rho)
    assert np.array_equal(wrapped["js_A_m2"], np.arange(4.0))


def test_attach_without_table_returns_original_catalog(tmp_path):
    path = tmp_path / "catalog.npz"
    np.savez(path, rho=np.ones(3))
    base = SimpleNamespace(files=["rho"], rho=np.ones(3))
    wrapped = attach_usadel_supercurrent_table_from_npz(base, path)
    assert wrapped is base
