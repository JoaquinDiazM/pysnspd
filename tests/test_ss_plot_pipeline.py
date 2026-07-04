from __future__ import annotations

from types import SimpleNamespace

import numpy as np

from pysnspd.analysis.ss_run import build_ss_plot_dataset
from pysnspd.plotting.ss_figures import make_ss_run_figures


class _Run(SimpleNamespace):
    pass


def _small_mesh():
    nodes = np.array(
        [
            [0.0, -0.5e-8],
            [1.0e-8, -0.5e-8],
            [2.0e-8, -0.5e-8],
            [0.0, 0.5e-8],
            [1.0e-8, 0.5e-8],
            [2.0e-8, 0.5e-8],
        ],
        dtype=float,
    )
    triangles = np.array([[0, 1, 3], [1, 4, 3], [1, 2, 4], [2, 5, 4]], dtype=np.int64)
    return SimpleNamespace(
        nodes=nodes,
        triangles=triangles,
        length_m=2.0e-8,
        width_m=1.0e-8,
        target_spacing_m=1.0e-8,
    )


def test_build_ss_plot_dataset_and_figures(tmp_path):
    mesh = _small_mesh()
    n = mesh.nodes.shape[0]
    psi = np.ones(n, dtype=float) * 2.0e-22
    state = {
        "psi_real_J": psi,
        "psi_imag_J": np.zeros(n),
        "phi_V": np.linspace(-1.0e-3, 1.0e-3, n),
        "node_jtot_x_A_m2": np.ones(n) * 2.0,
        "node_jtot_y_A_m2": np.zeros(n),
        "node_js_us_x_A_m2": np.ones(n) * 1.5,
        "node_js_us_y_A_m2": np.zeros(n),
        "node_jn_x_A_m2": np.ones(n) * 0.5,
        "node_jn_y_A_m2": np.zeros(n),
        "node_pairbreaking_ratio": np.linspace(0.0, 1.0, n),
        "node_div_jtot_A_m3": np.zeros(n),
    }
    history = {
        "t_s": np.array([0.0, 1.0e-12, 2.0e-12]),
        "dt_s": np.array([1.0e-15, 1.5e-15, 2.0e-15]),
        "eta_R": np.array([1.0e-2, 1.0e-3, 1.0e-4]),
        "terminal_voltage_V": np.array([0.0, 1.0e-3, 2.0e-3]),
        "normal_current_max_A_m2": np.array([1.0, 1.0, 1.0]),
        "total_current_max_A_m2": np.array([2.0, 2.0, 2.0]),
        "pairbreaking_max": np.array([1.0, 1.0, 1.0]),
        "delta0_meV": np.array([1.3]),
        "javg_A_m2": np.array([2.0]),
        "normal_terminal_node_mask": np.array([True, False, True, True, False, True]),
    }
    summary = {
        "solver": {
            "first_magic_ready": True,
            "target_current_A": 1.0e-6,
            "stationarity": {"passes": True, "bulk_exclusion_length_m": 2.0e-9},
            "contact_recovery": {"passes": True},
            "continuity": {"passes": True},
        }
    }
    run = _Run(
        run_name="unit_plot",
        pre_run_name="pre_unit",
        raw_ss=tmp_path / "raw" / "unit_plot" / "ss",
        figures_dir=tmp_path / "plots" / "unit_plot" / "figures",
        mesh=mesh,
        edge_data=SimpleNamespace(),
        state=state,
        history=history,
        summary=summary,
    )
    dataset = build_ss_plot_dataset(run)
    assert dataset["delta_over_delta0"].shape == (n,)
    assert "profiles" not in dataset

    saved = make_ss_run_figures(mesh=mesh, dataset=dataset, output_dir=run.figures_dir, dpi=80)
    assert set(saved) == {"overview", "relaxation", "adaptive", "masks"}
    for path in saved.values():
        assert path.exists()
