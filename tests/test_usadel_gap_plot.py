from pathlib import Path

import numpy as np

from pysnspd.plotting.usadel_gap import (
    interpolate_gap_curves,
    load_usadel_gap_catalog,
    plot_gap_eq_vs_temperature,
)


def test_load_and_interpolate_synthetic_usadel_gap_catalog(tmp_path: Path) -> None:
    T = np.linspace(0.9, 8.65, 9)
    q = np.linspace(0.0, 9.0e7, 10)
    Tc = 8.65
    qcrit = q[-1]
    gap_meV = 1.30 * np.sqrt(np.clip(1.0 - T[:, None] / Tc, 0.0, None)) * np.sqrt(
        np.clip(1.0 - (q[None, :] / qcrit) ** 2, 0.0, None)
    )
    npz = tmp_path / "usadel_catalog.npz"
    np.savez(npz, T_values_K=T, q_values_m_inv=q, delta_eq_meV=gap_meV)

    catalog = load_usadel_gap_catalog(npz)
    q_targets, T_out, curves = interpolate_gap_curves(catalog, n_curves=6)

    assert catalog.source_key == "delta_eq_meV"
    assert q_targets.shape == (6,)
    assert T_out.shape == T.shape
    assert curves.shape == (6, T.size)
    assert np.all(np.isfinite(curves[:, 0]))


def test_plot_gap_eq_vs_temperature_writes_pdf(tmp_path: Path) -> None:
    T = np.linspace(0.9, 8.65, 9)
    q = np.linspace(0.0, 9.0e7, 10)
    gap_J = (1.0e-3 / 1.0e3 * 1.602176634e-19) * np.ones((T.size, q.size))
    npz = tmp_path / "usadel_catalog.npz"
    np.savez(npz, temperature_values_K=T, q_values_m_inv=q, delta_eq_J=gap_J)

    catalog = load_usadel_gap_catalog(npz)
    output = plot_gap_eq_vs_temperature(catalog, tmp_path / "gap.pdf", n_curves=6)

    assert output.exists()
    assert output.suffix == ".pdf"
    assert output.stat().st_size > 0
