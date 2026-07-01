"""Parallel invariance smoke tests for PRE Usadel catalogue pieces."""
from __future__ import annotations

import numpy as np

from pysnspd.usadel.solver import compute_dos_grid
from pysnspd.usadel.supercurrent_table import build_matsubara_supercurrent_table_3d


def test_compute_dos_grid_parallel_matches_serial():
    energy = np.linspace(0.0, 5.0e-22, 64)
    deltas = np.array([0.0, 1.0e-22, 2.0e-22])
    gammas = np.array([0.0, 1.0e-23, 2.0e-23])

    rho_s, anom_s = compute_dos_grid(
        energy,
        deltas,
        gammas,
        eta_J=1.0e-25,
        workers=1,
        backend="serial",
    )
    rho_p, anom_p = compute_dos_grid(
        energy,
        deltas,
        gammas,
        eta_J=1.0e-25,
        workers=2,
        backend="process",
    )

    assert np.allclose(rho_p, rho_s, rtol=0.0, atol=0.0)
    assert np.allclose(anom_p, anom_s, rtol=0.0, atol=0.0)


def test_supercurrent_table_parallel_matches_serial():
    kwargs = dict(
        Te_axis_K=np.array([0.9, 2.0]),
        delta_axis_J=np.array([0.0, 1.0e-22, 2.0e-22]),
        q_axis_m_inv=np.array([0.0, 1.0e7, 2.0e7]),
        D_m2_s=1.58e-4,
        sigma_n_S_m=4.2e5,
        n_matsubara=8,
    )

    serial = build_matsubara_supercurrent_table_3d(
        **kwargs,
        workers=1,
        backend="serial",
    )
    parallel = build_matsubara_supercurrent_table_3d(
        **kwargs,
        workers=2,
        backend="process",
    )

    assert parallel.metadata["parallel_tasks"] == 6
    assert parallel.metadata["workers"] == 2
    assert parallel.metadata["parallel_backend"] == "process"
    assert np.allclose(parallel.js_T_delta_q_A_m2, serial.js_T_delta_q_A_m2)
