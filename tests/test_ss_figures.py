from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np

from pysnspd.analysis.ss_run import SSRunData, build_ss_plot_dataset
from pysnspd.plotting.ss_figures import make_ss_run_figures


@dataclass
class DummyMesh:
    nodes: np.ndarray
    triangles: np.ndarray
    length_m: float = 2.0e-7
    width_m: float = 1.0e-7
    target_spacing_m: float = 5.0e-8

    @property
    def n_nodes(self) -> int:
        return int(self.nodes.shape[0])


def _dummy_mesh() -> DummyMesh:
    x = np.array([0.0, 50e-9, 100e-9, 150e-9, 200e-9])
    y0 = np.zeros_like(x)
    y1 = np.ones_like(x) * 100e-9
    nodes = np.column_stack([np.r_[x, x], np.r_[y0, y1]])
    triangles = np.array(
        [
            [0, 1, 5],
            [0, 5, 4],
            [1, 2, 6],
            [1, 6, 5],
            [2, 3, 7],
            [2, 7, 6],
            [3, 4, 8],
            [3, 8, 7],
            [4, 9, 8],
        ],
        dtype=np.int64,
    )
    return DummyMesh(nodes=nodes, triangles=triangles)


def test_make_ss_run_figures_no_x_profiles(tmp_path: Path):
    mesh = _dummy_mesh()
    n = mesh.nodes.shape[0]
    dataset = {
        "run_name": "dummy",
        "x_nm": mesh.nodes[:, 0] * 1.0e9,
        "y_nm": mesh.nodes[:, 1] * 1.0e9,
        "triangles": mesh.triangles,
        "delta_over_delta0": np.linspace(0.7, 1.0, n),
        "phi_mV": np.linspace(-0.5, 0.5, n),
        "javg_A_m2": 2.0e10,
        "jtot_mag_A_m2": np.ones(n) * 2.0e10,
        "js_us_mag_A_m2": np.ones(n) * 1.5e10,
        "jn_mag_A_m2": np.ones(n) * 0.5e10,
        "pairbreaking_ratio": np.linspace(0.0, 0.2, n),
        "normal_terminal_node_mask": np.zeros(n, dtype=bool),
        "bulk_node_mask": np.ones(n, dtype=bool),
        "t_ps": np.array([0.0, 0.5, 1.0]),
        "eta_R": np.array([1.0e-2, 1.0e-3, 1.0e-4]),
        "pairbreaking_max_history": np.array([0.3, 0.25, 0.2]),
        "tdgl_probe_voltage_t_ps": np.array([0.0, 0.5, 1.0]),
        "tdgl_probe_voltage_mV": np.array([0.0, 0.8, 1.0]),
        "normal_current_fraction": np.array([0.0, 0.1, 0.2]),
        "dt_accepted_fs": np.array([0.3, 0.3, 0.4]),
        "dt_next_fs": np.array([0.3, 0.4, 0.5]),
        "dt_attempt_fs": np.array([0.3, 0.4, 0.5]),
        "adaptive_target_dt_fs": np.array([0.3, 0.4, 0.5]),
        "adaptive_retries": np.array([0, 1, 0]),
        "adaptive_rejected_attempts": np.array([0, 1, 1]),
        "adaptive_window_mean_d_abs_sq": np.array([1.0e-2, 1.0e-3, 1.0e-4]),
    }

    saved = make_ss_run_figures(
        mesh=mesh,
        dataset=dataset,
        output_dir=tmp_path,
        dpi=80,
    )

    assert "profiles" not in saved
    assert not (tmp_path / "ss_x_profiles.png").exists()
    assert set(saved) == {"overview", "relaxation", "adaptive", "masks"}
    for path in saved.values():
        assert path.exists()
        assert path.stat().st_size > 0


def test_build_ss_plot_dataset_uses_center_probe_voltage(tmp_path: Path):
    mesh = _dummy_mesh()
    n = mesh.nodes.shape[0]
    raw_ss = tmp_path / "raw" / "dummy" / "ss"
    raw_ss.mkdir(parents=True)

    t_s = np.array([0.0, 1.0e-12, 2.0e-12])
    x_m = mesh.nodes[:, 0]
    phi_snap = np.vstack([(k + 1.0) * 1.0e4 * x_m for k in range(t_s.size)])
    np.savez(raw_ss / "stationary_snapshots.npz", snapshot_t_s=t_s, phi_snapshot_V=phi_snap)

    state = {
        "psi_real_J": np.ones(n) * 2.0e-22,
        "psi_imag_J": np.zeros(n),
        "phi_V": phi_snap[-1],
        "node_jtot_x_A_m2": np.ones(n) * 2.0e10,
        "node_jtot_y_A_m2": np.zeros(n),
        "node_js_us_x_A_m2": np.ones(n) * 1.5e10,
        "node_js_us_y_A_m2": np.zeros(n),
        "node_jn_x_A_m2": np.ones(n) * 0.5e10,
        "node_jn_y_A_m2": np.zeros(n),
        "node_pairbreaking_ratio": np.zeros(n),
        "node_div_jtot_A_m3": np.zeros(n),
    }
    history = {
        "t_s": t_s,
        "dt_s": np.ones_like(t_s) * 1.0e-12,
        "eta_R": np.ones_like(t_s) * 1.0e-3,
        "terminal_voltage_V": np.ones_like(t_s) * 99.0,
        "normal_current_max_A_m2": np.ones_like(t_s),
        "total_current_max_A_m2": np.ones_like(t_s) * 2.0,
        "pairbreaking_max": np.zeros_like(t_s),
        "javg_A_m2": np.array([2.0e10]),
        "normal_terminal_node_mask": np.zeros(n, dtype=bool),
    }
    run = SSRunData(
        run_name="dummy",
        pre_run_name="pre_dummy",
        raw_ss=raw_ss,
        figures_dir=tmp_path / "figures",
        mesh=mesh,
        edge_data=None,
        state=state,
        history=history,
        summary={"solver": {"delta0_meV": 1.3, "javg_A_m2": 2.0e10}},
    )

    dataset = build_ss_plot_dataset(run)

    # Mesh center is x=100 nm, so probes are x=50 nm and x=150 nm.
    # With phi=(k+1)*1e4*x, right-left is (k+1)*1e-3 V = (k+1) mV.
    np.testing.assert_allclose(dataset["tdgl_probe_voltage_mV"], np.array([1.0, 2.0, 3.0]))
    np.testing.assert_allclose(dataset["tdgl_probe_voltage_t_ps"], np.array([0.0, 1.0, 2.0]))
    assert dataset["tdgl_probe_left_node_count"] == 2
    assert dataset["tdgl_probe_right_node_count"] == 2
    assert "profiles" not in dataset
    assert "x_profile_nm" not in dataset
