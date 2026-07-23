from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from types import SimpleNamespace

import numpy as np

from pysnspd.thermal.evolution import (
    PowerTableRuntimeInterpolator,
    ThermalRuntimeConfig,
    ThermalRuntimeController,
    build_central_thermal_mask,
    thermal_stationarity_diagnostics,
)


@dataclass
class DummyOps:
    edge_i: np.ndarray
    edge_j: np.ndarray
    edge_length_m: np.ndarray
    dual_face_length_m: np.ndarray
    node_area_m2: np.ndarray


def _power_table(path: Path) -> Path:
    Te = np.array([0.9, 2.0, 4.0])
    Tph = np.array([0.9, 2.0, 4.0])
    delta = np.array([0.0, 2.0e-22])
    q = np.array([0.0, 1.0e8])
    shape4 = (Te.size, Tph.size, delta.size, q.size)
    shape3 = (Te.size, delta.size, q.size)
    P_total = np.zeros(shape4)
    for i, te in enumerate(Te):
        for j, tph in enumerate(Tph):
            P_total[i, j, :, :] = 1.0e10 * (te - tph)
    np.savez(
        path,
        Te_values_K=Te,
        Tph_values_K=Tph,
        delta_values_J=delta,
        q_values_m_inv=q,
        P_S_W_m3=P_total,
        P_R_W_m3=np.zeros(shape4),
        P_total_W_m3=P_total,
        u_e_J_m3=np.zeros(shape3),
        C_e_J_m3_K=np.ones(shape3) * 1.0e3,
        kappa_s_W_m_K=np.ones((Te.size, delta.size)) * 1.0e-4,
        u_ph_J_m3=np.zeros(Tph.size),
        C_ph_J_m3_K=np.ones(Tph.size) * 2.0e3,
        P_esc_W_m3=np.zeros(Tph.size),
    )
    return path


def test_central_thermal_mask_default_width():
    nodes = np.column_stack([np.linspace(0.0, 360e-9, 7), np.zeros(7)])
    mask = build_central_thermal_mask(nodes, window_m=100e-9)
    assert mask.tolist() == [False, False, False, True, False, False, False]


def test_power_table_runtime_interpolator_constant_shapes(tmp_path: Path):
    interp = PowerTableRuntimeInterpolator(_power_table(tmp_path / "power_table_catalog.npz"))
    out = interp.evaluate(
        Te_K=np.array([1.45, 3.0]),
        Tph_K=np.array([0.9, 2.0]),
        delta_J=np.array([1.0e-22, 1.5e-22]),
        q_abs_m_inv=np.array([2.0e7, 7.0e7]),
    )
    assert out.P_total_W_m3.shape == (2,)
    assert out.C_e_J_m3_K.shape == (2,)
    assert np.all(out.C_e_J_m3_K > 0.0)
    assert np.all(out.C_ph_J_m3_K > 0.0)


def test_thermal_runtime_controller_heats_only_active_nodes(tmp_path: Path):
    nodes = np.column_stack([np.linspace(0.0, 200e-9, 5), np.zeros(5)])
    ops = DummyOps(
        edge_i=np.array([0, 1, 2, 3]),
        edge_j=np.array([1, 2, 3, 4]),
        edge_length_m=np.ones(4) * 50e-9,
        dual_face_length_m=np.ones(4) * 50e-9,
        node_area_m2=np.ones(5) * (50e-9 * 100e-9),
    )
    material = SimpleNamespace(delta0_J=2.0e-22, sigma_n_S_m=4.0e5)
    Te = np.ones(5) * 0.9
    Tph = np.ones(5) * 0.9
    controller = ThermalRuntimeController(
        nodes_m=nodes,
        ops=ops,
        material=material,
        Te_K=Te,
        Tph_K=Tph,
        power_table_npz=_power_table(tmp_path / "power_table_catalog.npz"),
        config=ThermalRuntimeConfig(enabled=True, window_m=100e-9, start_time_s=0.0, bath_K=0.9),
    )
    diag = controller.step(
        time_s=1.0e-12,
        dt_s=1.0e-15,
        psi_dimensionless=np.ones(5, dtype=np.complex128),
        native_normal_current=np.ones(4) * 1.0e3,
        current_scale_A_m2=1.0,
    )
    assert diag["thermal_active"] == 1.0
    assert np.allclose(Te[~controller.mask], 0.9)
    assert float(np.max(Te[controller.mask])) >= 0.9


def test_thermal_stationarity_diagnostics_tail_rate():
    hist = {
        "t_s": np.linspace(0.0, 10e-12, 20),
        "thermal_max_rate_K_per_ps": np.r_[np.ones(10) * 1.0, np.ones(10) * 1.0e-3],
    }
    diag = thermal_stationarity_diagnostics(
        hist,
        enabled=True,
        start_time_s=2e-12,
        requested_total_time_s=10e-12,
        rate_tol_K_per_ps=1.0e-2,
    )
    assert diag["passes"] is True
