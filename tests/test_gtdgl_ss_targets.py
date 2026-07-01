"""Tests for the first stationary-state target diagnostics."""
from __future__ import annotations

import numpy as np

from pysnspd.gtdgl.material import GTDGLMaterial
from pysnspd.gtdgl.ss_targets import (
    apply_terminal_proximity_seed,
    contact_recovery_diagnostics,
    stationarity_diagnostics,
)


def _material() -> GTDGLMaterial:
    return GTDGLMaterial(
        Tc_K=8.65,
        D_m2_s=1.58e-4,
        sigma_n_S_m=4.2e5,
        delta0_J=2.10667708314e-22,
        thickness_m=7e-9,
        width_m=120e-9,
        tau_ee_Tc_s=5.0e-12,
        tau_ep_Tc_s=24.7e-12,
        tau_scale=0.10,
    )


def test_terminal_proximity_seed_reaches_bulk_after_few_xi():
    material = _material()
    x = np.linspace(0.0, 240e-9, 401)
    nodes = np.column_stack([x, np.zeros_like(x)])
    psi = np.ones(x.size, dtype=np.complex128)
    Te = np.full(x.size, 0.9)

    out, seed_diag = apply_terminal_proximity_seed(
        psi,
        nodes_m=nodes,
        material=material,
        Te_K=Te,
        healing_target_xi=2.5,
        target_bulk_fraction=0.95,
    )
    rec = contact_recovery_diagnostics(
        psi_dimensionless=out,
        nodes_m=nodes,
        material=material,
        Te_K=Te,
        threshold_fraction=0.95,
        min_allowed_xi=1.5,
        max_allowed_xi=4.0,
        bin_width_m=0.5e-9,
    )

    assert seed_diag.enabled
    assert np.isclose(np.abs(out[0]), 0.0)
    assert np.max(np.abs(out)) <= 1.0 + 1e-12
    assert rec.passes
    assert 1.5 <= rec.left_recovery_length_xi <= 4.0
    assert 1.5 <= rec.right_recovery_length_xi <= 4.0


def test_stationarity_diagnostics_detects_steady_snapshots():
    material = _material()
    psi0 = np.ones((3, 5), dtype=np.complex128) * material.delta0_J
    psi0[-1] *= 1.0 + 1.0e-7
    phi = np.zeros((3, 5), dtype=float)
    phi[-1] = 1.0e-12
    history = {
        "psi_snapshot_real_J": np.real(psi0),
        "psi_snapshot_imag_J": np.imag(psi0),
        "phi_snapshot_V": phi,
        "eta_R": np.array([1.0e-3, 1.0e-6, 1.0e-6]),
    }

    diag = stationarity_diagnostics(
        history=history,
        material=material,
        delta_rel_tol=1.0e-4,
        phi_rel_tol=2.0,
        eta_tol=1.0e-5,
        eta_window=2,
    )

    assert diag.passes
    assert diag.delta_rel_change < 1.0e-4
    assert diag.eta_R_final < 1.0e-5
