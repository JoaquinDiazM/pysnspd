from __future__ import annotations

from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

from pysnspd.io.manager import create_run_layout


def _minimal_config(root: Path) -> dict:
    return {
        "project": {
            "name": "storage_parallel_test",
            "big_data_root": str(root),
            "default_run_name": "default_run",
        },
        "parallel": {"enabled": True, "workers": 4, "backend": "process"},
        "material": {
            "name": "NbN",
            "Tc_K": 8.65,
            "sigma_n_S_m": 4.2e5,
            "lambda_L_m": 5.4e-7,
            "thickness_m": 7.0e-9,
            "width_m": 120.0e-9,
            "tau_ee_Tc_ps": 0.5,
            "tau_ep_Tc_ps": 2.47,
        },
        "calibration": {"Ic_target_A": 38.8e-6},
        "bias": {"T_bias_K": 0.9, "I_bias_A": 20.0e-6},
        "mesh": {"type": "delaunay", "target_spacing_m": 4.1e-9, "seed": 222222},
        "catalogs": {
            "dos": {"n_delta": 3, "n_q": 3, "n_energy": 32, "n_matsubara": 16},
            "phase_space": {"n_Te": 3, "n_Tph": 3, "n_delta": 3, "n_q": 3, "n_omega": 16},
        },
        "ss_run": {"max_steps": 10, "dt_s": 1.0e-15, "convergence_tol": 1.0e-7},
        "photon_run": {
            "photon_wavelength_m": 1064.0e-9,
            "max_steps": 10,
            "dt_s": 1.0e-15,
            "bubble_radius_m": 10.0e-9,
        },
        "circuit": {"R_load_ohm": 50.0, "L_bias_H": 1.0e-6, "C_rf_F": 1.0e-12},
    }


def _create_layout_worker(root_text: str, index: int) -> str:
    cfg = _minimal_config(Path(root_text))
    layout = create_run_layout(cfg, f"parallel_run_{index:02d}")
    return layout["raw_ss"]


def test_create_run_layout_is_process_safe(tmp_path: Path) -> None:
    root = tmp_path / "big_data"
    n = 8
    with ProcessPoolExecutor(max_workers=4) as pool:
        paths = list(pool.map(_create_layout_worker, [str(root)] * n, range(n)))

    assert len(paths) == n
    for path in paths:
        assert Path(path).is_dir()

    # The write probe uses unique temporary names and should leave no fixed
    # shared marker behind after many concurrent layout creations.
    assert not (root / ".pysnspd_write_test").exists()
