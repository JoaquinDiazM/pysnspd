from __future__ import annotations

import numpy as np

from pysnspd.analysis.usadel_low_amplitude_diagnostic import (
    DiagnosticCurrentCatalog,
    build_candidate_stiffness_table,
    build_notch_diagnostic,
    compare_constitutive_curves,
)
from pysnspd.usadel.supercurrent_table import build_matsubara_supercurrent_table_3d


def _small_catalog() -> DiagnosticCurrentCatalog:
    Te = np.array([1.5, 4.0])
    delta0 = 2.0e-22
    delta = delta0 * np.array([0.0, 0.08, 0.25, 0.60, 1.0])
    q = np.array([0.0, 2.0e6, 6.0e6])
    D = 1.58e-4
    sigma = 4.2e5
    n_matsubara = 96
    table = build_matsubara_supercurrent_table_3d(
        Te_axis_K=Te,
        delta_axis_J=delta,
        q_axis_m_inv=q,
        D_m2_s=D,
        sigma_n_S_m=sigma,
        n_matsubara=n_matsubara,
        workers=1,
        backend="serial",
    )
    stiffness = build_candidate_stiffness_table(
        Te_axis_K=Te,
        delta_axis_J=delta,
        q_axis_m_inv=q,
        js_A_m2=table.js_T_delta_q_A_m2,
        D_m2_s=D,
        sigma_n_S_m=sigma,
        n_matsubara=n_matsubara,
    )
    return DiagnosticCurrentCatalog(
        Te_axis_K=Te,
        delta_axis_J=delta,
        q_axis_m_inv=q,
        js_A_m2=table.js_T_delta_q_A_m2,
        stiffness_A_m_inv_J_inv2=stiffness,
        D_m2_s=D,
        sigma_n_S_m=sigma,
        Tc_K=8.65,
        delta0_J=delta0,
        n_matsubara=n_matsubara,
        metadata={},
    )


def test_stiffness_reconstruction_restores_quadratic_low_amplitude_limit() -> None:
    catalog = _small_catalog()
    first = catalog.first_positive_delta_J
    delta = np.geomspace(1.0e-4 * first, 0.25 * first, 120)
    curves = compare_constitutive_curves(
        catalog,
        Te_K=1.5,
        q_m_inv=2.0e6,
        delta_J=delta,
    )

    assert 0.98 < curves.metrics["current_low_amplitude_exponent"] < 1.02
    assert 1.98 < curves.metrics["stiffness_low_amplitude_exponent"] < 2.02
    assert 1.98 < curves.metrics["direct_low_amplitude_exponent"] < 2.02
    assert (
        curves.metrics["stiffness_max_relative_error_below_first_node"]
        < curves.metrics["current_max_relative_error_below_first_node"]
    )


def test_candidate_table_has_finite_positive_zero_delta_and_zero_q_anchors() -> None:
    catalog = _small_catalog()
    stiffness = catalog.stiffness_A_m_inv_J_inv2

    assert np.all(np.isfinite(stiffness))
    assert np.all(stiffness > 0.0)
    assert np.all(stiffness[:, 0, :] > 0.0)
    assert np.all(stiffness[:, :, 0] > 0.0)


def test_current_interpolation_amplifies_smooth_notch_source() -> None:
    catalog = _small_catalog()
    notch = build_notch_diagnostic(
        catalog,
        Te_K=1.5,
        q_m_inv=2.0e6,
        n_points=401,
    )

    assert notch.metrics["current_source_peak_over_direct_peak"] > 2.0
    assert notch.metrics["stiffness_source_peak_over_direct_peak"] < 1.5
    assert (
        notch.metrics["stiffness_source_rms_relative_to_direct"]
        < notch.metrics["current_source_rms_relative_to_direct"]
    )
