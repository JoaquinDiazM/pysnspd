from __future__ import annotations

from dataclasses import fields
from pathlib import Path
from typing import Any

import numpy as np
import pytest
from scipy.constants import Boltzmann

from pysnspd.plotting.power_diagnostics import (
    AllmarasReferenceParameters,
    PowerTablePlotCatalog,
    phonon_dos_normalization_diagnostic,
    resolve_allmaras_reference_parameters,
)


N0_J_M3 = 2.0e47
TC_K = 8.65
REPO_ROOT = Path(__file__).resolve().parents[1]


def _catalog(*, energy_scale: float = 1.0, capacity_scale: float = 1.0) -> PowerTablePlotCatalog:
    Te = np.linspace(0.9, 34.6, 121)
    Tph = np.linspace(0.9, 34.6, 41)
    delta = np.asarray([0.0, 1.0e-22])
    q = np.asarray([0.0, 1.0e7])
    shape3 = (Te.size, delta.size, q.size)
    shape4 = (Te.size, Tph.size, delta.size, q.size)
    omega = np.asarray([1.0, 2.0, 3.0]) * 1.0e-22

    analytic_u = np.pi**2 * N0_J_M3 * Boltzmann**2 * Te**2 / 3.0
    analytic_C = 2.0 * np.pi**2 * N0_J_M3 * Boltzmann**2 * Te / 3.0
    u_e = np.broadcast_to(analytic_u[:, None, None], shape3).copy()
    C_e = np.broadcast_to(analytic_C[:, None, None], shape3).copy()
    u_e[:, 0, 0] *= energy_scale
    C_e[:, 0, 0] *= capacity_scale

    return PowerTablePlotCatalog(
        Te_values_K=Te,
        Tph_values_K=Tph,
        delta_values_J=delta,
        q_values_m_inv=q,
        P_S_W_m3=np.zeros(shape4),
        P_R_W_m3=np.zeros(shape4),
        P_total_W_m3=np.zeros(shape4),
        u_e_J_m3=u_e,
        C_e_J_m3_K=C_e,
        kappa_s_W_m_K=np.ones((Te.size, delta.size)),
        u_ph_J_m3=np.ones(Tph.size),
        C_ph_J_m3_K=np.ones(Tph.size),
        u_ph_weighted_J=np.ones(Tph.size),
        C_ph_weighted_J_K=np.ones(Tph.size),
        omega_values_J=omega,
        alpha2F=np.asarray([0.1, 0.2, 0.15]),
        phdos_states_per_THz=np.ones(omega.size),
        metadata={
            "material": {"N0_J_m3": N0_J_M3},
            "delta0_J": 1.764 * Boltzmann * TC_K,
        },
    )


def _reference(**overrides: float) -> AllmarasReferenceParameters:
    kwargs = {
        "Tc_K": TC_K,
        "normal_dos_spin_convention": "single_spin",
        "gamma_heat_capacity": 60.0,
        "tau_ep_Tc_s": 24.7e-12,
        "tau_esc_s": 20.0e-12,
    }
    kwargs.update(overrides)
    return resolve_allmaras_reference_parameters(_catalog(), **kwargs)


def test_allmaras_electron_energy_derivative_is_heat_capacity() -> None:
    reference = _reference()
    temperature = np.linspace(0.9, 34.6, 5001)
    numerical = np.gradient(reference.u_e_normal_J_m3(temperature), temperature)
    analytic = reference.C_e_normal_J_m3_K(temperature)
    assert np.allclose(numerical[1:-1], analytic[1:-1], rtol=2.0e-11)


def test_allmaras_phonon_energy_derivative_is_heat_capacity() -> None:
    reference = _reference()
    temperature = np.linspace(0.9, 34.6, 5001)
    numerical = np.gradient(reference.u_ph_J_m3(temperature), temperature)
    analytic = reference.C_ph_J_m3_K(temperature)
    assert np.allclose(numerical[1:-1], analytic[1:-1], rtol=1.0e-4)


def test_allmaras_gamma_is_reconstructed_from_capacities_at_Tc() -> None:
    reference = _reference()
    recovered = (
        8.0
        * np.pi**2
        / 5.0
        * reference.C_e_normal_J_m3_K(reference.Tc_K)
        / reference.C_ph_J_m3_K(reference.Tc_K)
    )
    assert np.isclose(recovered, reference.gamma_heat_capacity, rtol=2.0e-15)


def test_allmaras_power_has_equilibrium_zero_and_declared_sign() -> None:
    reference = _reference()
    assert reference.P_eph_W_m3(4.0, 4.0) == pytest.approx(0.0)
    assert reference.P_eph_W_m3(5.0, 4.0) > 0.0
    assert reference.P_eph_W_m3(4.0, 5.0) < 0.0


def test_allmaras_fourth_power_equation_matches_phonon_energy_balance() -> None:
    reference = _reference()
    Te_K = np.asarray([3.0, 8.0, 15.0])
    Tph_K = np.asarray([2.0, 5.0, 10.0])
    Tb_K = 0.9
    A_A = reference.phonon_storage_coefficient_J_m3_K4
    P_escape = A_A * (Tph_K**4 - Tb_K**4) / reference.tau_esc_s
    rate_from_energy = (reference.P_eph_W_m3(Te_K, Tph_K) - P_escape) / A_A
    rate_from_temperature = reference.phonon_temperature_fourth_power_rate_K4_s(
        Te_K, Tph_K, Tb_K
    )
    assert np.allclose(rate_from_energy, rate_from_temperature, rtol=2.0e-15)


def test_allmaras_reference_validates_catalogue_normal_electronic_limit() -> None:
    reference = _reference()
    validation = reference.manifest_dict()["normal_limit_validation"]
    assert validation["passed"] is True
    assert validation["max_relative_error_u_e"] < 1.0e-14
    assert validation["max_relative_error_C_e"] < 1.0e-14


def test_allmaras_reference_rejects_catalogue_normalization_mismatch() -> None:
    with pytest.raises(ValueError, match="normal-like production storage"):
        resolve_allmaras_reference_parameters(
            _catalog(energy_scale=2.0, capacity_scale=2.0),
            Tc_K=TC_K,
            normal_dos_spin_convention="single_spin",
        )


def test_allmaras_reference_rejects_an_undeclared_spin_conversion() -> None:
    with pytest.raises(ValueError, match="single-spin"):
        resolve_allmaras_reference_parameters(
            _catalog(),
            Tc_K=TC_K,
            normal_dos_spin_convention="two_spin",
        )


def test_allmaras_reference_cross_checks_Tc_against_persisted_delta0() -> None:
    with pytest.raises(ValueError, match="configured Tc_K does not match"):
        resolve_allmaras_reference_parameters(
            _catalog(),
            Tc_K=TC_K + 0.1,
            normal_dos_spin_convention="single_spin",
        )


def test_allmaras_override_sources_do_not_claim_evaluated_values_are_published() -> None:
    reference = _reference(
        gamma_heat_capacity=61.0,
        tau_ep_Tc_s=25.0e-12,
        tau_esc_s=21.0e-12,
    )
    sources = reference.manifest_dict()["sources"]
    assert sources["gamma_heat_capacity"].startswith("evaluated CLI override")
    assert "published gamma=60" in sources["gamma_heat_capacity"]
    assert sources["tau_ep_Tc_s"].startswith("evaluated CLI override")
    assert "published tau_ep(Tc)=24.7 ps" in sources["tau_ep_Tc_s"]
    assert sources["tau_esc_s"].startswith("evaluated CLI override")
    assert "published tau_esc=20 ps" in sources["tau_esc_s"]


def test_allmaras_reference_has_no_debye_or_material_density_dependency() -> None:
    field_names = {item.name for item in fields(AllmarasReferenceParameters)}
    forbidden = {
        "ion_density_m3",
        "omega_D_J",
        "debye_cutoff_meV",
        "lambda_ep",
        "lambda_provenance",
        "tau0_s",
        "runtime_catalog_spectrum",
    }
    assert field_names.isdisjoint(forbidden)
    manifest_keys = set(_reference().manifest_dict())
    assert manifest_keys.isdisjoint(forbidden)


def test_allmaras_manifest_is_complete_and_finite() -> None:
    manifest = _reference().manifest_dict()
    assert manifest["model"] == "Allmaras 2019 normal-state two-temperature reference"
    assert manifest["material_normalization"] == "production N0 and Tc"
    assert manifest["gamma_heat_capacity"] == pytest.approx(60.0)
    assert manifest["tau_ep_Tc_s"] == pytest.approx(24.7e-12)
    assert manifest["tau_esc_s"] == pytest.approx(20.0e-12)
    assert manifest["published_gamma"] == pytest.approx(60.0)
    assert manifest["published_tau_ep_Tc_ps"] == pytest.approx(24.7)
    assert manifest["published_tau_esc_ps"] == pytest.approx(20.0)
    assert manifest["evaluated_gamma"] == pytest.approx(60.0)
    assert manifest["evaluated_tau_ep_Tc_ps"] == pytest.approx(24.7)
    assert manifest["evaluated_tau_esc_ps"] == pytest.approx(20.0)
    assert manifest["normal_dos_spin_convention"] == "single_spin"
    assert manifest["doi"] == "10.1103/PhysRevApplied.11.034062"
    _assert_serializable_and_finite(manifest)


def test_phonon_dos_normalization_is_a_separate_nonrescaling_audit() -> None:
    diagnostic = phonon_dos_normalization_diagnostic(_catalog())
    assert diagnostic["integral_F_dnu_states_per_material_entity"] > 0.0
    assert diagnostic["material_entity"].startswith("ion")
    assert "not rescaled" in diagnostic["policy"]
    forbidden = {"gamma_heat_capacity", "tau_ep_Tc_s", "lambda_ep", "omega_D_J"}
    assert set(diagnostic).isdisjoint(forbidden)


def test_E1_exposes_only_the_allmaras_reference_cli() -> None:
    source = (REPO_ROOT / "plot_pipelines" / "E1_plot_prerun.py").read_text(
        encoding="utf-8"
    )
    for option in (
        "--allmaras-reference",
        "--no-allmaras-reference",
        "--allmaras-gamma",
        "--allmaras-tau-ep-Tc-ps",
    ):
        assert option in source
    for obsolete in (
        "--debye-reference",
        "--no-debye-reference",
        "--debye-cutoff-meV",
    ):
        assert obsolete not in source


def test_eliashberg_axis_treats_alpha2F_as_dimensionless() -> None:
    source = (
        REPO_ROOT / "pysnspd" / "plotting" / "eliashberg_spectrum.py"
    ).read_text(encoding="utf-8")
    assert 'set_ylabel(r"$\\alpha^2F(\\Omega)$"' in source
    assert '$\\alpha^2F(\\Omega)$ [meV]' not in source


def _assert_serializable_and_finite(value: Any) -> None:
    if isinstance(value, dict):
        for key, nested in value.items():
            assert isinstance(key, str)
            _assert_serializable_and_finite(nested)
        return
    if isinstance(value, (list, tuple)):
        for nested in value:
            _assert_serializable_and_finite(nested)
        return
    if isinstance(value, (bool, str)):
        return
    if isinstance(value, (int, float, np.integer, np.floating)):
        assert np.isfinite(value)
        return
    raise AssertionError(f"non-serializable manifest value: {value!r}")
