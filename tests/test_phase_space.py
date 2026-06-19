from pathlib import Path

import numpy as np

from pysnspd.config import validate_config
from pysnspd.usadel.catalog import build_usadel_catalog_from_config
from pysnspd.kinetic.phase_space import (
    build_phase_space_catalog_from_usadel_catalog,
    fermi_positive_energy,
    load_phase_space_catalog_npz,
    pair_recombination_thermal_factor,
    phase_space_summary,
    save_phase_space_catalog_npz,
)


def _minimal_config(tmp_path: Path) -> dict:
    return {
        "project": {
            "name": "test_project",
            "big_data_root": str(tmp_path / "big_data"),
            "default_run_name": "phase_space_test",
        },
        "parallel": {
            "enabled": True,
            "workers": 2,
            "backend": "process",
        },
        "material": {
            "name": "NbN",
            "Tc_K": 8.65,
            "sigma_n_S_m": "4.2e5",
            "lambda_L_m": 5.4e-7,
            "thickness_m": 7.0e-9,
            "width_m": 120.0e-9,
        },
        "calibration": {
            "Ic_target_A": 38.8e-6,
            "n_gamma_sweep": 20,
            "gamma_max_fraction": 0.80,
            "D_warn_min_m2_s": 5.0e-5,
            "D_warn_max_m2_s": 5.0e-4,
        },
        "bias": {
            "T_bias_K": 0.9,
            "I_bias_A": 35.0e-6,
        },
        "mesh": {
            "type": "delaunay",
            "length_m": 240.0e-9,
            "target_spacing_m": 20.0e-9,
            "seed": 12345,
        },
        "catalogs": {
            "dos": {
                "n_delta": 5,
                "n_q": 5,
                "n_energy": 120,
                "n_matsubara": 30,
            },
            "phase_space": {
                "n_Te": 4,
                "n_Tph": 4,
                "n_delta": 4,
                "n_q": 4,
                "n_omega": 40,
            },
        },
        "ss_run": {
            "max_steps": 10,
            "dt_s": 1.0e-15,
            "convergence_tol": 1.0e-7,
        },
        "photon_run": {
            "photon_wavelength_m": 1064.0e-9,
            "max_steps": 10,
            "dt_s": 1.0e-15,
            "bubble_radius_m": 10.0e-9,
        },
        "circuit": {
            "R_load_ohm": 50.0,
            "L_bias_H": 1.0e-6,
            "C_rf_F": 1.0e-12,
        },
    }


def test_fermi_and_pair_factor_are_finite():
    E = np.linspace(0.0, 1.0e-21, 100)
    f = fermi_positive_energy(E, 5.0)

    assert f.shape == E.shape
    assert np.all(np.isfinite(f))
    assert np.all(f >= 0.0)
    assert np.all(f <= 0.5)

    factor = pair_recombination_thermal_factor(E + 1.0e-24, E + 2.0e-24, 5.0)

    assert factor.shape == E.shape
    assert np.all(np.isfinite(factor))
    assert np.all(factor >= 0.0)


def test_build_save_load_phase_space_catalog(tmp_path):
    cfg = validate_config(
        _minimal_config(tmp_path),
        require_big_data_root_exists=False,
    )

    usadel_catalog = build_usadel_catalog_from_config(
        cfg,
        gamma_max_fraction=0.80,
        energy_max_factor=6.0,
    )

    phase_catalog = build_phase_space_catalog_from_usadel_catalog(
        usadel_catalog,
        cfg,
        n_Te=3,
        n_delta=3,
        n_q=3,
        n_omega=30,
        Te_min_K=0.9,
        Te_max_K=20.0,
    )

    assert phase_catalog.J_S_TdqO_J.shape == (3, 3, 3, 30)
    assert phase_catalog.J_R_TdqO_J.shape == (3, 3, 3, 30)
    assert np.all(np.isfinite(phase_catalog.J_S_TdqO_J))
    assert np.all(np.isfinite(phase_catalog.J_R_TdqO_J))
    assert np.min(phase_catalog.J_S_TdqO_J) >= 0.0
    assert np.min(phase_catalog.J_R_TdqO_J) >= 0.0

    summary = phase_space_summary(phase_catalog)
    assert summary["shape"] == [3, 3, 3, 30]
    assert summary["J_S_is_finite"] is True
    assert summary["J_R_is_finite"] is True

    path = save_phase_space_catalog_npz(
        phase_catalog,
        tmp_path / "phase_space_catalog.npz",
    )
    loaded = load_phase_space_catalog_npz(path)

    assert np.allclose(phase_catalog.Te_values_K, loaded.Te_values_K)
    assert np.allclose(phase_catalog.omega_values_J, loaded.omega_values_J)
    assert np.allclose(phase_catalog.J_S_TdqO_J, loaded.J_S_TdqO_J)
    assert np.allclose(phase_catalog.J_R_TdqO_J, loaded.J_R_TdqO_J)
    assert loaded.metadata["backend"] == "phase_space_from_usadel_dos_oe4_v1_2"


def test_recombination_threshold_and_normal_limit_are_explicit():
    from pysnspd.kinetic.phase_space import recombination_phase_space_spectrum

    meV_J = 1.602176634e-22
    E = np.linspace(0.0, 10.0, 500) * meV_J
    delta = 1.0 * meV_J
    rho = np.ones_like(E)

    omega = np.array([0.0, 1.0, 1.999, 2.0, 2.5, 4.0]) * delta
    JR = recombination_phase_space_spectrum(
        E,
        rho,
        omega,
        Te_K=6.0,
        delta_J=delta,
    )

    assert np.all(np.isfinite(JR))
    assert np.all(JR >= 0.0)
    assert np.allclose(JR[:4], 0.0)
    assert np.any(JR[4:] > 0.0)

    JR_normal = recombination_phase_space_spectrum(
        E,
        rho,
        omega,
        Te_K=6.0,
        delta_J=0.0,
    )
    assert np.allclose(JR_normal, 0.0)


def test_scattering_zero_omega_is_zero_and_positive_omega_is_finite():
    from pysnspd.kinetic.phase_space import scattering_phase_space_spectrum

    meV_J = 1.602176634e-22
    E = np.linspace(0.0, 10.0, 500) * meV_J
    delta = 1.0 * meV_J
    rho = np.ones_like(E)

    omega = np.array([0.0, 0.5, 2.0]) * meV_J
    JS = scattering_phase_space_spectrum(
        E,
        rho,
        omega,
        Te_K=6.0,
        delta_J=delta,
    )

    assert np.all(np.isfinite(JS))
    assert np.all(JS >= 0.0)
    assert np.isclose(JS[0], 0.0, atol=1.0e-35)


def test_phase_space_uses_trapezoid_not_deprecated_trapz():
    source = Path("pysnspd/kinetic/phase_space.py").read_text(encoding="utf-8")
    assert "np.trapz" not in source
    assert "np.trapezoid" in source


def test_scattering_has_finite_dos_support_cutoff():
    from pysnspd.kinetic.phase_space import scattering_phase_space_spectrum

    meV_J = 1.602176634e-22
    E = np.linspace(0.0, 10.0, 501) * meV_J
    rho = np.ones_like(E)
    delta = 2.0 * meV_J

    # For scattering E' = E + Omega. With E_max = 10 meV and Delta = 2 meV,
    # the numerical integration support collapses when Omega > 8 meV.
    # This is a finite-catalogue limitation relative to the Appendix-A integral,
    # whose upper energy limit is formally infinite.
    omega = np.array([0.0, 4.0, 7.9, 8.1]) * meV_J
    JS = scattering_phase_space_spectrum(
        E,
        rho,
        omega,
        Te_K=6.0,
        delta_J=delta,
    )

    assert np.all(np.isfinite(JS))
    assert np.all(JS >= 0.0)
    assert JS[-1] == 0.0
    assert np.any(JS[1:-1] > 0.0)


def test_phase_space_energy_window_summary_reports_truncation():
    from types import SimpleNamespace

    from pysnspd.kinetic.phase_space import phase_space_energy_window_summary

    meV_J = 1.602176634e-22
    catalog = SimpleNamespace(
        omega_values_J=np.linspace(0.0, 10.0, 11) * meV_J,
        delta_values_J=np.array([0.0, 2.0]) * meV_J,
        metadata={
            "energy_max_J": 10.0 * meV_J,
            "js_hard_cutoff_by_delta_J": [10.0 * meV_J, 8.0 * meV_J],
            "jr_threshold_by_delta_J": [0.0, 4.0 * meV_J],
        },
    )

    summary = phase_space_energy_window_summary(catalog)

    assert summary["energy_max_meV"] == 10.0
    assert summary["J_S_hard_cutoff_min_meV"] == 8.0
    assert summary["J_R_threshold_max_meV"] == 4.0
    assert summary["scattering_window_is_truncated"] is True


def test_phase_space_source_documents_pdf_appendix_caveats():
    source = Path("pysnspd/kinetic/phase_space.py").read_text(encoding="utf-8")

    assert "PDF Appendix A mismatch" in source
    assert "finite Usadel DOS grid" in source
    assert "N1(E)N1(E')" in source
    assert "f(E)->f_FD(E,Te)" in source

