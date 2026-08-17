from __future__ import annotations

import numpy as np
import pytest

from pysnspd.plotting.power_diagnostics import (
    PowerTablePlotCatalog,
    resolve_debye_reference_parameters,
)


def _catalog() -> PowerTablePlotCatalog:
    Te = np.asarray([1.0, 2.0])
    Tph = np.asarray([1.0, 2.0])
    delta = np.asarray([0.0, 1.0])
    q = np.asarray([0.0, 1.0])
    shape3 = (Te.size, delta.size, q.size)
    shape4 = (Te.size, Tph.size, delta.size, q.size)
    omega = np.asarray([1.0, 2.0, 3.0]) * 1.0e-22
    return PowerTablePlotCatalog(
        Te_values_K=Te,
        Tph_values_K=Tph,
        delta_values_J=delta,
        q_values_m_inv=q,
        P_S_W_m3=np.zeros(shape4),
        P_R_W_m3=np.zeros(shape4),
        P_total_W_m3=np.zeros(shape4),
        u_e_J_m3=np.ones(shape3),
        C_e_J_m3_K=np.ones(shape3),
        kappa_s_W_m_K=np.ones((Te.size, delta.size)),
        u_ph_J_m3=np.ones(Tph.size),
        C_ph_J_m3_K=np.ones(Tph.size),
        u_ph_weighted_J=np.ones(Tph.size),
        C_ph_weighted_J_K=np.ones(Tph.size),
        omega_values_J=omega,
        alpha2F=np.asarray([0.1, 0.2, 0.15]),
        phdos_states_per_THz=np.ones(omega.size),
        metadata={"material": {"N0_J_m3": 2.0e47}},
    )


def test_debye_reference_recovers_and_persists_all_constants() -> None:
    canonical_lambda = 1.2156653376542017
    reference = resolve_debye_reference_parameters(
        _catalog(),
        ion_density_m3=48.0e27,
        Tc_K=8.65,
        omega_D_J=30.0 * 1.602176634e-22,
        normal_dos_spin_convention="single_spin",
        lambda_ep=canonical_lambda,
        lambda_provenance=_complete_spectrum_provenance(),
    )

    manifest = reference.manifest_dict()
    assert manifest["enabled"] is True
    assert manifest["normal_dos_spin_convention"] == "single_spin"
    assert np.isclose(manifest["omega_D_meV"], 30.0)
    assert reference.lambda_ep == canonical_lambda
    assert manifest["lambda_provenance"]["n_points"] == 11999
    assert manifest["lambda_provenance"]["frequency_max_THz"] == 20.0
    assert manifest["lambda_provenance"]["sha256"] == "a" * 64
    assert manifest["runtime_catalog_spectrum"]["role"].startswith("truncated")
    assert not np.isclose(
        manifest["runtime_catalog_spectrum"]["lambda_if_reintegrated_diagnostic_only"],
        canonical_lambda,
    )
    assert reference.tau0_s > 0.0


def test_debye_reference_rejects_an_undeclared_spin_conversion() -> None:
    with pytest.raises(ValueError, match="single-spin"):
        resolve_debye_reference_parameters(
            _catalog(),
            ion_density_m3=48.0e27,
            Tc_K=8.65,
            omega_D_J=30.0 * 1.602176634e-22,
            normal_dos_spin_convention="two_spin",
            lambda_ep=1.2156653376542017,
            lambda_provenance=_complete_spectrum_provenance(),
        )


def test_debye_reference_requires_auditable_complete_spectrum_provenance() -> None:
    with pytest.raises(ValueError, match="missing required fields"):
        resolve_debye_reference_parameters(
            _catalog(),
            ion_density_m3=48.0e27,
            Tc_K=8.65,
            omega_D_J=30.0 * 1.602176634e-22,
            normal_dos_spin_convention="single_spin",
            lambda_ep=1.2156653376542017,
            lambda_provenance={"source_path": "nbn-a2f-ph.dat"},
        )


def _complete_spectrum_provenance() -> dict[str, object]:
    return {
        "definition": "2*integral(alpha2F(nu)/nu dnu) over the complete source spectrum",
        "loader": "pysnspd.kinetic.eliashberg.load_simon_eliashberg_dat",
        "source_path": "/catalogs/simon_2025/nbn-a2f-ph.dat",
        "sha256": "a" * 64,
        "n_points": 11999,
        "frequency_min_THz": 0.0,
        "frequency_max_THz": 20.0,
    }
