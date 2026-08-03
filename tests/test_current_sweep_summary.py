from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace

import numpy as np

from pysnspd.analysis import ss_run as ss_run_module
from pysnspd.analysis.ss_run import build_ss_plot_dataset
from pysnspd.plotting.current_sweep import (
    _summary_diagnostics,
    collect_current_sweep_iv_points,
)
from pysnspd.plotting.current_sweep_iv import _build_iv_point
from pysnspd.plotting.current_sweep_summary import (
    build_current_sweep_regime_summary,
    plot_current_sweep_regime_summary,
)


def _point(current: float, **extra):
    point = {
        "run_name": f"sweep_I{current:g}uA",
        "current_uA": current,
        "complete": True,
        "strict_stationarity_passes": False,
        "dynamic_stationarity_passes": False,
        "photon_ready": False,
        "approximately_ohmic": False,
        "normal_voltage_ratio": 0.2,
        "terminal_normal_voltage_ratio": 0.2,
        "mean_delta_over_delta0": 0.9,
        "normal_like_fraction_final": 0.0,
        "accepted_steps": 1000,
        "rejected_steps": 250,
        "rejected_over_accepted": 0.25,
        "final_time_ps": 200.0,
    }
    point.update(extra)
    return point


def test_regime_summary_keeps_sampled_ranges_and_incomplete_cases(tmp_path: Path) -> None:
    points = [
        {"run_name": "synthetic_origin", "current_uA": 0.0},
        _point(
            20.0,
            dynamic_stationarity_passes=True,
            photon_ready=True,
            approximately_ohmic=True,
            normal_voltage_ratio=0.96,
            terminal_normal_voltage_ratio=1.03,
        ),
        _point(30.0),
    ]
    skipped = [
        {
            "run_name": "sweep_I40uA",
            "current_uA": 40.0,
            "complete": False,
            "reason": "missing ss_summary.yaml",
        }
    ]

    summary = build_current_sweep_regime_summary(
        points,
        skipped,
        ohmic_relative_tolerance=0.10,
    )
    assert summary["counts"]["discovered"] == 3
    assert summary["counts"]["complete"] == 2
    assert summary["sampled_ranges"]["photon_ready"]["currents_uA"] == [20.0]
    assert summary["sampled_ranges"]["approximately_ohmic"]["currents_uA"] == [20.0]
    assert summary["incomplete_currents_uA"] == [40.0]

    output = plot_current_sweep_regime_summary(
        points,
        skipped,
        tmp_path / "summary.pdf",
        ohmic_relative_tolerance=0.10,
    )
    assert output.exists()
    assert output.stat().st_size > 1000


def test_incomplete_record_is_rejected_without_loading_npz() -> None:
    record = {
        "run_name": "sweep_I20to100uA_dI_plus13uA_I33uA",
        "stages": {
            "ss": {
                "exists": True,
                "npz_files": [{"relative_path": "ss/stationary_seed.npz"}],
                "summary_files": [],
            }
        },
    }
    points, skipped, meta = collect_current_sweep_iv_points(
        config_path="not-read.yaml",
        project_config={},
        records=[record],
        include_origin=False,
    )
    assert points == []
    assert meta["n_runs_loaded"] == 0
    assert skipped[0]["current_uA"] == 33.0
    assert "ss_summary.yaml" in skipped[0]["reason"]
    assert "stationary_state.npz" in skipped[0]["reason"]


def test_legacy_final_photon_gate_is_reconstructed_from_summary() -> None:
    summary = {
        "solver": {
            "requested_time_ps": 200.0,
            "final_time_ps": 200.0,
            "requested_time_reached": True,
            "stationarity": {"passes": False},
            "dynamic_stationarity": {
                "passes": True,
                "normal_like_fraction_final": 0.25,
            },
            "contact_recovery": {"passes": True},
            "continuity": {"passes": True},
            "thermal_stationarity": {"passes": True},
            "circuit_stationarity": {"passes": True},
            "allmaras_phase_continuation": {"final_converged": True},
            "accepted_steps": 100,
            "rejected_steps": 25,
        }
    }
    diagnostics = _summary_diagnostics(summary)
    assert diagnostics["complete"] is True
    assert diagnostics["dynamic_stationarity_passes"] is True
    assert diagnostics["photon_ready"] is True
    assert diagnostics["photon_ready_source"] == "reconstructed_final_gate"
    assert diagnostics["rejected_over_accepted"] == 0.25


def test_endpoint_dataset_does_not_touch_snapshot_archive(monkeypatch, tmp_path: Path) -> None:
    def fail_if_called(*args, **kwargs):
        raise AssertionError("snapshot archive must not be opened")

    monkeypatch.setattr(
        ss_run_module,
        "_center_probe_voltage_from_snapshots",
        fail_if_called,
    )
    run = SimpleNamespace(
        run_name="sweep_I20uA",
        pre_run_name="pre",
        raw_ss=tmp_path,
        mesh=SimpleNamespace(
            nodes=np.asarray([[0.0, 0.0], [1.0e-7, 0.0], [2.0e-7, 0.0]]),
            triangles=np.asarray([[0, 1, 2]]),
            width_m=1.0e-7,
        ),
        state={
            "psi_real_J": np.full(3, 1.0e-22),
            "psi_imag_J": np.zeros(3),
            "phi_V": np.asarray([0.0, 1.0e-3, 2.0e-3]),
        },
        history={},
        summary={"solver": {"delta0_meV": 1.0, "target_current_A": 20.0e-6}},
    )
    dataset = build_ss_plot_dataset(run, load_snapshots=False)
    assert np.allclose(dataset["phi_mV"], [0.0, 1.0, 2.0])


def test_terminal_voltage_prefers_solver_result_over_zero_seed_value() -> None:
    nodes = np.asarray([[0.0, 0.0], [1.0e-7, 0.0], [2.0e-7, 0.0]])
    run = SimpleNamespace(
        summary={
            "seed": {"terminal_voltage_V": 0.0},
            "solver": {"target_current_A": 20.0e-6, "terminal_voltage_V": 3.0e-3},
        },
        mesh=SimpleNamespace(nodes=nodes, width_m=1.0e-7, length_m=2.0e-7),
        pre_run_name="pre",
        raw_ss=Path("/tmp/run/ss"),
    )
    dataset = {
        "x_nm": nodes[:, 0] * 1.0e9,
        "phi_mV": np.asarray([0.0, 0.1, 0.2]),
        "summary_scalars": {"target_current_A": 20.0e-6},
    }
    point, _, _, _ = _build_iv_point(
        run_name="sweep_I20uA",
        run=run,
        dataset=dataset,
        project_config={
            "material": {"thickness_m": 5.0e-9, "sigma_n_S_m": 1.0e7}
        },
        voltage_probe_offset_nm=50.0,
        voltage_probe_half_window_nm=1.0,
    )
    assert point["terminal_voltage_mV"] == 3.0
    assert point["source"].endswith("ss_summary_terminal_voltage_V")
