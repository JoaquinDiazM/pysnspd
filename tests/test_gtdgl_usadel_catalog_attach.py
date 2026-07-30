"""NPZ attachment tests for PRE Usadel stiffness resources."""
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
        js_A_m2=np.zeros((1, 2, 3)),
        js_stiffness_A_per_m_J2=np.ones((1, 2, 3)),
        q_axis_m_inv=np.array([0.0, 1.0, 2.0]),
        delta_axis_J=np.array([0.0, 2.0]),
        delta2_axis_J2=np.array([0.0, 4.0]),
        Te_axis_K=np.array([1.0]),
        object_metadata=np.array(["ignored"], dtype=object),
    )
    arrays = load_usadel_supercurrent_table_arrays_npz(path)
    assert set(arrays) >= {
        "js_A_m2",
        "js_stiffness_A_per_m_J2",
        "q_axis_m_inv",
        "delta2_axis_J2",
    }
    assert arrays["js_stiffness_A_per_m_J2"].shape == (1, 2, 3)


def test_attach_usadel_supercurrent_table_from_npz_wraps_base_catalog(tmp_path):
    path = tmp_path / "catalog.npz"
    np.savez(
        path,
        js_A_m2=np.zeros((1, 2, 2)),
        js_stiffness_A_per_m_J2=np.ones((1, 2, 2)),
        Te_axis_K=np.array([1.0]),
        delta_axis_J=np.array([0.0, 1.0]),
        delta2_axis_J2=np.array([0.0, 1.0]),
        q_axis_m_inv=np.array([0.0, 1.0]),
    )
    base = SimpleNamespace(files=["rho"], rho=np.ones(3))
    wrapped = attach_usadel_supercurrent_table_from_npz(base, path)
    assert isinstance(wrapped, UsadelCatalogWithSupercurrentTable)
    assert "js_stiffness_A_per_m_J2" in wrapped.files
    assert np.array_equal(wrapped.rho, base.rho)
    assert np.array_equal(wrapped["js_stiffness_A_per_m_J2"], np.ones((1, 2, 2)))


def test_attach_without_table_returns_original_catalog(tmp_path):
    path = tmp_path / "catalog.npz"
    np.savez(path, rho=np.ones(3))
    base = SimpleNamespace(files=["rho"], rho=np.ones(3))
    wrapped = attach_usadel_supercurrent_table_from_npz(base, path)
    assert wrapped is base
