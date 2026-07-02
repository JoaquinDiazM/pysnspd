"""Smoke tests for PRE diagnostic plotting helpers."""

from __future__ import annotations

from types import SimpleNamespace

import numpy as np

from pysnspd.mesh.delaunay import MeshData
from pysnspd.mesh.edges import build_edge_data
from pysnspd.plotting.pre_diagnostics import write_pre_diagnostic_plots


def test_write_pre_diagnostic_plots_smoke(tmp_path):
    nodes = np.array(
        [
            [0.0, -0.5],
            [1.0, -0.5],
            [1.0, 0.5],
            [0.0, 0.5],
            [0.5, 0.0],
        ],
        dtype=float,
    )
    triangles = np.array(
        [
            [0, 1, 4],
            [1, 2, 4],
            [2, 3, 4],
            [3, 0, 4],
        ],
        dtype=np.int64,
    )

    mesh = MeshData(
        nodes=nodes,
        triangles=triangles,
        length_m=1.0,
        width_m=1.0,
        target_spacing_m=0.5,
        seed=12345,
        triangulation_method="test",
        boundary_guard_layers=0,
    )
    edge_data = build_edge_data(nodes, triangles, length_m=1.0, width_m=1.0)

    q_values = np.linspace(0.0, 1.2e8, 7)
    delta_values = np.linspace(0.0, 2.2e-22, 6)
    energy_values = np.linspace(0.0, 6.0e-22, 32)
    calibration_q = np.linspace(0.0, 1.0e8, 8)
    calibration_delta = 2.0e-22 * np.clip(1.0 - (calibration_q / 1.0e8) ** 2, 0.0, None)
    js_table_temperature_values_K = np.array([0.9, 4.0, 8.0], dtype=float)
    js_table_delta_eq_values_J = np.vstack([
        calibration_delta,
        0.65 * calibration_delta,
        0.15 * calibration_delta,
    ])

    rho = np.empty((delta_values.size, q_values.size, energy_values.size), dtype=float)
    anomalous = np.empty_like(rho)
    for i_delta, delta in enumerate(delta_values):
        for i_q, q in enumerate(q_values):
            gap = max(delta * (1.0 - 0.15 * q / max(q_values[-1], 1.0)), 0.0)
            e_scaled = energy_values / max(gap, 1.0e-24)
            rho[i_delta, i_q, :] = 0.02 + np.sqrt(np.maximum(e_scaled - 1.0, 0.0))
            anomalous[i_delta, i_q, :] = gap / (energy_values + gap + 1.0e-24)

    usadel_catalog = SimpleNamespace(
        calibration_q_values_m_inv=calibration_q,
        calibration_current_values_A=np.array(
            [0.0, 5.0e-6, 1.0e-5, 2.0e-5, 3.0e-5, 3.8e-5, 3.6e-5, 3.2e-5],
            dtype=float,
        ),
        calibration_delta_eq_values_J=calibration_delta,
        js_table_temperature_values_K=js_table_temperature_values_K,
        js_table_q_values_m_inv=calibration_q,
        js_table_delta_eq_values_J=js_table_delta_eq_values_J,
        q_values_m_inv=q_values,
        delta_values_J=delta_values,
        energy_values_J=energy_values,
        rho_delta_gamma_E=rho,
        anomalous_delta_gamma_E=anomalous,
        metadata={
            "Ic_target_A": 3.88e-5,
            "T_bias_K": 0.9,
            "Tc_K": 8.65,
            "D_m2_s": 1.58e-4,
            "sigma_n_S_m": 4.2e5,
            "delta0_meV": 1.315,
        },
    )

    outputs = write_pre_diagnostic_plots(
        mesh=mesh,
        edge_data=edge_data,
        usadel_catalog=usadel_catalog,
        output_dir=tmp_path,
        dpi=80,
    )

    assert set(outputs) == {
        "mesh_triangle_area_hist_png",
        "mesh_edge_length_hist_png",
        "usadel_supercurrent_curve_png",
        "usadel_equilibrium_dos_map_png",
        "usadel_zero_energy_dos_map_png",
        "usadel_equilibrium_anomalous_map_png",
        "usadel_equilibrium_gap_Tq_map_png",
    }
    for value in outputs.values():
        path = tmp_path / value.split("/")[-1]
        assert path.exists()
        assert path.stat().st_size > 0
