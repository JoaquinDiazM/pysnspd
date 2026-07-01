"""Smoke tests for strict 3D Usadel supercurrent tables."""
from __future__ import annotations

import numpy as np

from pysnspd.usadel.supercurrent_table import build_matsubara_supercurrent_table_3d
from pysnspd.gtdgl.usadel_current import (
    interpolate_strict_usadel_current_table,
    validate_strict_usadel_supercurrent_table_npz,
)


def test_build_matsubara_supercurrent_table_3d_shape_and_axes():
    table = build_matsubara_supercurrent_table_3d(
        Te_axis_K=np.array([0.9]),
        delta_axis_J=np.array([0.0, 2.0e-22]),
        q_axis_m_inv=np.array([0.0, 1.0e7, 2.0e7]),
        D_m2_s=1.58e-4,
        sigma_n_S_m=4.2e5,
        n_matsubara=8,
        workers=1,
    )

    assert table.js_T_delta_q_A_m2.shape == (1, 2, 3)
    assert np.allclose(table.js_T_delta_q_A_m2[:, 0, :], 0.0)
    assert np.allclose(table.js_T_delta_q_A_m2[:, :, 0], 0.0)
    assert table.js_T_delta_q_A_m2[0, 1, 2] > table.js_T_delta_q_A_m2[0, 1, 1] > 0.0


def test_interpolator_requires_full_3d_layout_and_is_odd_in_q():
    Te = np.array([1.0, 2.0])
    Delta = np.array([0.0, 1.0])
    q = np.array([0.0, 10.0])
    js = np.zeros((2, 2, 2), dtype=float)
    js[:, 1, 1] = np.array([100.0, 200.0])

    out = interpolate_strict_usadel_current_table(
        table=js,
        Te_axis_K=Te,
        delta_axis_J=Delta,
        q_axis_m_inv=q,
        q_edge_m_inv=np.array([5.0, -5.0]),
        delta_edge_J=np.array([0.5, 0.5]),
        Te_edge_K=np.array([1.5, 1.5]),
    )

    assert np.allclose(out, [37.5, -37.5])


def test_validate_strict_table_npz_rejects_legacy_1d(tmp_path):
    path = tmp_path / "legacy.npz"
    np.savez(path, js_A_m2=np.zeros(4), q_axis_m_inv=np.linspace(0.0, 1.0, 4))

    try:
        validate_strict_usadel_supercurrent_table_npz(path)
    except RuntimeError as exc:
        assert "not strict 3D" in str(exc) or "layout" in str(exc)
    else:  # pragma: no cover
        raise AssertionError("legacy 1D table was accepted")


def test_validate_strict_table_npz_accepts_T_delta_q(tmp_path):
    path = tmp_path / "strict.npz"
    np.savez(
        path,
        js_A_m2=np.zeros((1, 2, 3)),
        Te_axis_K=np.array([0.9]),
        delta_axis_J=np.array([0.0, 1.0]),
        q_axis_m_inv=np.array([0.0, 1.0, 2.0]),
    )

    summary = validate_strict_usadel_supercurrent_table_npz(path)
    assert summary["valid"] is True
    assert summary["layout"] == "Te,delta,q"
    assert summary["shape"] == [1, 2, 3]
