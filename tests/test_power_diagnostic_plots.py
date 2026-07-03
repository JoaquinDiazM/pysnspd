"""Smoke tests for PRE power-table diagnostic plotting helpers."""

from __future__ import annotations

from pathlib import Path

import numpy as np

from pysnspd.plotting.power_diagnostics import write_power_table_diagnostic_plots


def test_write_power_table_diagnostic_plots_smoke(tmp_path: Path) -> None:
    Te = np.linspace(0.9, 8.0, 5)
    Tph = np.linspace(0.9, 8.0, 5)
    delta = np.linspace(0.0, 2.0e-22, 4)
    q = np.linspace(0.0, 8.0e7, 3)
    omega = np.linspace(0.0, 6.0e-22, 12)

    shape = (Te.size, Tph.size, delta.size, q.size)
    Te4 = Te[:, None, None, None]
    Tph4 = Tph[None, :, None, None]
    delta4 = delta[None, None, :, None]
    q4 = q[None, None, None, :]
    scale = 1.0 + delta4 / max(delta[-1], 1.0e-30) + q4 / max(q[-1], 1.0)
    P_S = 1.0e10 * (Te4 - Tph4) * scale
    P_R = 2.0e10 * (Te4 - Tph4) * (delta4 / max(delta[-1], 1.0e-30))
    P_total = P_S + P_R

    u_e = np.empty((Te.size, delta.size, q.size), dtype=float)
    C_e = np.empty_like(u_e)
    for i, T in enumerate(Te):
        u_e[i] = 100.0 * T**2 + 1.0e22 * delta[:, None] - 1.0e-7 * q[None, :]
        C_e[i] = 200.0 * T + 10.0

    path = tmp_path / "power_table_catalog.npz"
    np.savez_compressed(
        path,
        Te_values_K=Te,
        Tph_values_K=Tph,
        delta_values_J=delta,
        q_values_m_inv=q,
        P_S_W_m3=P_S,
        P_R_W_m3=P_R,
        P_total_W_m3=P_total,
        u_e_J_m3=u_e,
        C_e_J_m3_K=C_e,
        u_ph_weighted_J=Tph**4,
        C_ph_weighted_J_K=4.0 * Tph**3,
        omega_values_J=omega,
        alpha2F=np.linspace(0.0, 0.1, omega.size),
        phdos_states_per_THz=np.linspace(1.0, 2.0, omega.size),
        metadata=np.array({"Tc_K": 8.65}, dtype=object),
    )

    outputs = write_power_table_diagnostic_plots(
        power_table_npz=path,
        output_dir=tmp_path,
        dpi=80,
    )

    assert set(outputs) == {
        "power_channels_Te_Tph_maps_png",
        "power_total_Delta_q_maps_png",
        "power_total_Te_curves_png",
        "energy_heat_capacity_curves_png",
        "power_equal_temperature_residual_png",
    }
    for value in outputs.values():
        out = Path(value)
        assert out.exists()
        assert out.stat().st_size > 0
