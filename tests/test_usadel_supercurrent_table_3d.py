"""Tests for the strict 3D Usadel stiffness/current PRE resources."""
from __future__ import annotations

import numpy as np
import pytest

from pysnspd.usadel.calibration import matsubara_energy_axis_J
from pysnspd.usadel.parameters import E_CHARGE_C, HBAR_J_S, K_B_J_K
from pysnspd.usadel.supercurrent_table import build_matsubara_supercurrent_table_3d
from pysnspd.gtdgl.usadel_current import (
    interpolate_strict_usadel_stiffness_table,
    validate_strict_usadel_supercurrent_table_npz,
)


def test_build_matsubara_table_contains_finite_stiffness_and_current():
    table = build_matsubara_supercurrent_table_3d(
        Te_axis_K=np.array([0.9]),
        delta_axis_J=np.array([0.0, 2.0e-22]),
        q_axis_m_inv=np.array([0.0, 1.0e7, 2.0e7]),
        D_m2_s=1.58e-4,
        sigma_n_S_m=4.2e5,
        n_matsubara=8,
        workers=1,
    )

    assert table.shape == (1, 2, 3)
    assert np.array_equal(table.delta2_axis_J2, table.delta_axis_J**2)
    assert np.all(np.isfinite(table.js_stiffness_T_delta2_q_A_per_m_J2))
    assert np.all(table.js_stiffness_T_delta2_q_A_per_m_J2 > 0.0)
    assert np.allclose(table.js_T_delta_q_A_m2[:, 0, :], 0.0)
    assert np.allclose(table.js_T_delta_q_A_m2[:, :, 0], 0.0)


def test_exact_zero_amplitude_stiffness_matches_matsubara_limit():
    temperature = 1.7
    diffusion = 1.58e-4
    conductivity = 4.2e5
    q_axis = np.array([0.0, 1.5e7])
    n_matsubara = 24
    table = build_matsubara_supercurrent_table_3d(
        Te_axis_K=np.array([temperature]),
        delta_axis_J=np.array([0.0, 1.0e-22]),
        q_axis_m_inv=q_axis,
        D_m2_s=diffusion,
        sigma_n_S_m=conductivity,
        n_matsubara=n_matsubara,
        workers=1,
    )

    eps = matsubara_energy_axis_J(T_K=temperature, n_matsubara=n_matsubara)
    gamma = 0.5 * HBAR_J_S * diffusion * q_axis**2
    prefactor = 2.0 * np.pi * K_B_J_K * temperature * conductivity / E_CHARGE_C
    expected = prefactor * np.sum(1.0 / (eps[:, None] + gamma[None, :]) ** 2, axis=0)

    assert np.allclose(
        table.js_stiffness_T_delta2_q_A_per_m_J2[0, 0],
        expected,
        rtol=2.0e-14,
        atol=0.0,
    )


def test_stiffness_interpolator_is_even_in_q_and_linear_in_delta_squared():
    Te = np.array([1.0, 2.0])
    delta2 = np.array([0.0, 4.0])
    q = np.array([0.0, 10.0])
    stiffness = np.empty((2, 2, 2), dtype=float)
    stiffness[0] = np.array([[2.0, 4.0], [6.0, 8.0]])
    stiffness[1] = 2.0 * stiffness[0]

    out = interpolate_strict_usadel_stiffness_table(
        table=stiffness,
        Te_axis_K=Te,
        delta2_axis_J2=delta2,
        q_axis_m_inv=q,
        q_edge_m_inv=np.array([5.0, -5.0]),
        delta2_edge_J2=np.array([2.0, 2.0]),
        Te_edge_K=np.array([1.5, 1.5]),
    )

    assert np.allclose(out, [7.5, 7.5])


def test_reconstructed_current_is_continuous_at_first_positive_delta2_node():
    first = 1.0
    points = np.array([np.nextafter(first, 0.0), first, np.nextafter(first, np.inf)])
    stiffness = interpolate_strict_usadel_stiffness_table(
        table=np.array([[[2.0, 2.0], [3.0, 3.0], [5.0, 5.0]]]),
        Te_axis_K=np.array([1.0]),
        delta2_axis_J2=np.array([0.0, first, 4.0]),
        q_axis_m_inv=np.array([0.0, 2.0]),
        q_edge_m_inv=np.ones(3),
        delta2_edge_J2=points,
        Te_edge_K=np.ones(3),
    )
    current = stiffness * points

    assert np.allclose(current, np.full(3, 3.0), rtol=4.0 * np.finfo(float).eps, atol=0.0)


def test_validate_strict_table_rejects_current_only_npz(tmp_path):
    path = tmp_path / "current_only.npz"
    np.savez(
        path,
        js_A_m2=np.zeros((1, 2, 3)),
        Te_axis_K=np.array([0.9]),
        delta_axis_J=np.array([0.0, 1.0]),
        q_axis_m_inv=np.array([0.0, 1.0, 2.0]),
    )

    with pytest.raises(RuntimeError, match="stiffness|v2"):
        validate_strict_usadel_supercurrent_table_npz(path)


def test_validate_strict_table_accepts_v2_contract(tmp_path):
    path = tmp_path / "strict_v2.npz"
    np.savez(
        path,
        js_A_m2=np.zeros((1, 2, 3)),
        js_stiffness_A_per_m_J2=np.ones((1, 2, 3)),
        Te_axis_K=np.array([0.9]),
        delta_axis_J=np.array([0.0, 1.0]),
        delta2_axis_J2=np.array([0.0, 1.0]),
        q_axis_m_inv=np.array([0.0, 1.0, 2.0]),
        js_table_backend=np.array("matsubara_usadel_stiffness_table_3d_v2"),
    )

    summary = validate_strict_usadel_supercurrent_table_npz(path)
    assert summary["valid"] is True
    assert summary["layout"] == "Te,delta2,q"
    assert summary["shape"] == [1, 2, 3]
