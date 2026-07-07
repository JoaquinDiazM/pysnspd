from pathlib import Path

import numpy as np

from pysnspd.plotting.usadel_gap import (
    interpolate_gap_curves,
    load_usadel_gap_catalog,
    plot_gap_eq_vs_temperature,
)


def _write_official_like_usadel_npz(path: Path) -> None:
    Tc = 8.65
    D = 1.58e-4
    qcrit = 8.0e7
    q = np.linspace(0.0, 1.1 * qcrit, 9)
    gamma = np.linspace(0.0, 1.0e-22, q.size)
    delta = np.linspace(0.0, 2.0e-22, 6)
    energy = np.linspace(0.0, 8.0e-22, 12)
    metadata = {
        "Tc_K": Tc,
        "T_bias_K": 0.9,
        "D_m2_s": D,
        "n_matsubara_configured": 80,
        "calibration": {
            "q_critical_m_inv": qcrit,
        },
    }
    np.savez_compressed(
        path,
        energy_values_J=energy,
        delta_values_J=delta,
        gamma_values_J=gamma,
        q_values_m_inv=q,
        rho_delta_gamma_E=np.ones((delta.size, gamma.size, energy.size)),
        anomalous_delta_gamma_E=np.zeros((delta.size, gamma.size, energy.size)),
        eta_J=np.array(1.0e-25),
        calibration_gamma_values_J=gamma,
        calibration_q_values_m_inv=q,
        calibration_delta_eq_values_J=np.maximum(delta[-1] * (1.0 - q / q[-1]), 0.0),
        calibration_current_values_A=np.linspace(0.0, 4.0e-5, q.size),
        calibration_current_density_values_A_m2=np.linspace(0.0, 4.0e10, q.size),
        calibration_sum_s2_values=np.ones(q.size),
        metadata=np.array(metadata, dtype=object),
    )


def test_reconstruct_gap_from_official_like_usadel_catalog(tmp_path: Path) -> None:
    npz = tmp_path / "usadel_dos_catalog.npz"
    _write_official_like_usadel_npz(npz)

    catalog = load_usadel_gap_catalog(npz, n_curves=6, n_temperature=12, n_matsubara=60)
    q_targets, T_out, curves = interpolate_gap_curves(catalog)

    assert catalog.source_key == "computed_from_matsubara_self_consistency"
    assert q_targets.shape == (6,)
    assert T_out.shape == (12,)
    assert curves.shape == (6, 12)
    assert np.all(np.isfinite(curves))
    assert curves[0, 0] > curves[0, -1]
    assert np.isclose(q_targets[-1], catalog.q_critical_m_inv)


def test_plot_gap_eq_vs_temperature_writes_pdf(tmp_path: Path) -> None:
    npz = tmp_path / "usadel_dos_catalog.npz"
    _write_official_like_usadel_npz(npz)

    catalog = load_usadel_gap_catalog(npz, n_curves=6, n_temperature=10, n_matsubara=50)
    output = plot_gap_eq_vs_temperature(catalog, tmp_path / "gap.pdf")

    assert output.exists()
    assert output.suffix == ".pdf"
    assert output.stat().st_size > 0
