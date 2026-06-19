from pathlib import Path

from pysnspd.config import load_config, validate_config
from pysnspd.io.manager import (
    initialize_project_storage,
    create_run_layout,
    write_manifest,
    read_manifest,
    resolve_stage_path,
    resolve_plot_path,
)


def test_config_and_storage_manager(tmp_path):
    config_path = tmp_path / "project.yaml"
    big_data_root = tmp_path / "big_data"

    config_path.write_text(
        f"""
project:
  name: test_project
  big_data_root: {big_data_root}
  default_run_name: NbN_test_run

parallel:
  enabled: true
  workers: 2
  backend: process

material:
  name: NbN
  Tc_K: 8.65
  sigma_n_S_m: 4.2e5
  lambda_L_m: 5.4e-7
  thickness_m: 7.0e-9
  width_m: 120.0e-9

calibration:
  Ic_target_A: 38.8e-6
  n_gamma_sweep: 40
  gamma_max_fraction: 0.80
  D_warn_min_m2_s: 5.0e-5
  D_warn_max_m2_s: 5.0e-4

bias:
  T_bias_K: 0.9
  I_bias_A: 35.0e-6

mesh:
  type: delaunay
  target_spacing_m: 4.0e-9
  seed: 12345

catalogs:
  dos:
    n_delta: 21
    n_q: 26
    n_energy: 2000
    n_matsubara: 500
  phase_space:
    n_Te: 40
    n_Tph: 40
    n_delta: 21
    n_q: 26
    n_omega: 1200

ss_run:
  max_steps: 10000
  dt_s: 1.0e-15
  convergence_tol: 1.0e-7

photon_run:
  photon_wavelength_m: 1064.0e-9
  max_steps: 20000
  dt_s: 1.0e-15
  bubble_radius_m: 10.0e-9

circuit:
  R_load_ohm: 50.0
  L_bias_H: 1.0e-6
  C_rf_F: 1.0e-12
""",
        encoding="utf-8",
    )

    cfg = load_config(config_path)
    cfg = validate_config(cfg, require_big_data_root_exists=False)

    base = initialize_project_storage(cfg)
    assert Path(base["raw"]).is_dir()
    assert Path(base["plots"]).is_dir()
    assert Path(base["logs"]).is_dir()
    assert Path(base["catalogs"]).is_dir()
    assert Path(base["tmp"]).is_dir()

    layout = create_run_layout(cfg)
    assert Path(layout["raw_run"]).is_dir()
    assert Path(layout["raw_pre"]).is_dir()
    assert Path(layout["raw_ss"]).is_dir()
    assert Path(layout["raw_photon"]).is_dir()
    assert Path(layout["plots_run"]).is_dir()
    assert Path(layout["plots_figures"]).is_dir()
    assert Path(layout["plots_mesh"]).is_dir()
    assert Path(layout["plots_diagnostics"]).is_dir()
    assert Path(layout["plots_comparisons"]).is_dir()
    assert Path(layout["logs_run"]).is_dir()

    manifest_path = write_manifest(
        cfg,
        stage="project",
        extra={"test": "pytest OE1"},
    )
    assert manifest_path.exists()

    manifest = read_manifest(cfg, stage="project")
    assert manifest["run"]["name"] == "NbN_test_run"
    assert manifest["run"]["stage"] == "project"
    assert manifest["extra"]["test"] == "pytest OE1"

    assert resolve_stage_path(cfg, stage="pre").is_dir()
    assert resolve_stage_path(cfg, stage="ss").is_dir()
    assert resolve_stage_path(cfg, stage="photon").is_dir()
    assert resolve_plot_path(cfg).is_dir()
