from __future__ import annotations

from types import SimpleNamespace

import numpy as np

from pysnspd.gtdgl.pytdgl_like.usadel_current import (
    attach_usadel_supercurrent_table_from_npz,
    load_usadel_supercurrent_table_arrays_npz,
)


def test_attach_usadel_table_arrays_preserves_base_catalog(tmp_path):
    path = tmp_path / "catalog.npz"
    np.savez(
        path,
        q_axis_m_inv=np.array([0.0, 1.0]),
        delta_axis_J=np.array([0.0, 2.0]),
        js_A_m2=np.array([[0.0, 1.0], [0.0, 2.0]]),
    )
    base = SimpleNamespace(delta_axis_J=np.array([123.0]), old_value="ok")
    wrapped = attach_usadel_supercurrent_table_from_npz(base, path)
    assert wrapped.old_value == "ok"
    assert np.array_equal(wrapped.q_axis_m_inv, np.array([0.0, 1.0]))
    assert np.array_equal(wrapped["js_A_m2"], np.array([[0.0, 1.0], [0.0, 2.0]]))
    assert "js_A_m2" in wrapped.files


def test_numeric_usadel_table_loader_ignores_unrelated_object_arrays(tmp_path):
    path = tmp_path / "catalog_with_object.npz"
    np.savez(
        path,
        q_axis_m_inv=np.array([0.0, 1.0]),
        delta_axis_J=np.array([0.0, 2.0]),
        js_A_m2=np.array([[0.0, 1.0], [0.0, 2.0]]),
        metadata=np.array([{"object": True}], dtype=object),
    )
    arrays = load_usadel_supercurrent_table_arrays_npz(path)
    assert set(["q_axis_m_inv", "delta_axis_J", "js_A_m2"]).issubset(arrays)
    assert "metadata" not in arrays
