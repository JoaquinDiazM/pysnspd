"""No-run fixtures for candidate source energy, normalized geometry and identity."""
import importlib.util
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest
from scipy.constants import Boltzmann as kB, elementary_charge as e, hbar

from pysnspd.experimental.effective_bcs_functional import EffectiveBCSFunctional, bcs_count_quadrature
from pysnspd.experimental.effective_bcs_kinetics import EffectiveBCSKinetics, acoustic_phonon_grid
from pysnspd.experimental.effective_material import allmaras_acoustic_material
from pysnspd.experimental.thermal_spatial_usadel import rectangular_graph
from pysnspd.experimental.thermal_stiffness import ThermalStiffnessFunctional

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location("candidate_thermal_runner_test", ROOT/"pipelines/09_thermal_stiffness_photon_run.py")
runner = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(runner)


@pytest.fixture(scope="module")
def physical_fixture():
    D, Rs, Tc, thickness = 1.4104e-4, 425., 8.65, 7e-9
    Eu = kB*Tc
    ell = np.sqrt(hbar*D/(2*Eu))
    N0 = 1/(2*e**2*D*Rs*thickness)
    graph = rectangular_graph(5, 3, 160e-9/ell, 80e-9/ell)
    xi, weights = bcs_count_quadrature()
    kwargs = dict(graph=graph, energy_scale_J=Eu, N0_per_J_m3=N0, ell0_m=ell,
        thickness_m=thickness, count_nodes=xi, count_weights=weights,
        gradient_law="bcs_radial", bath_temperature_K=.9, diffusion_m2_s=D)
    functional = ThermalStiffnessFunctional(**kwargs)
    bath = EffectiveBCSFunctional(**kwargs)
    xy = graph.coordinates_bar*ell
    left, right = np.flatnonzero(xy[:, 0] == xy[:, 0].min()), np.flatnonzero(xy[:, 0] == xy[:, 0].max())
    edge_lengths = np.linalg.norm(xy[graph.edges[:, 1]]-xy[graph.edges[:, 0]], axis=1)
    geometry = SimpleNamespace(graph=graph, coordinates_m=xy, left_nodes=left, right_nodes=right,
        triangles=np.empty((0, 3), dtype=int), edge_lengths_m=edge_lengths,
        dual_face_lengths_m=graph.conductance*edge_lengths, ell0_m=ell,
        source_path="fixture physical mesh", receipt={"mesh_sha256": "fixture"})
    material = allmaras_acoustic_material(Tc_K=Tc, tau0_s=1e-9, diffusivity_m2_s=D,
        sheet_resistance_ohm=Rs, thickness_m=thickness, nodes=129)
    phonons = acoustic_phonon_grid(material, Eu, nodes=24)
    kinetics = EffectiveBCSKinetics(xi_bar=xi, weights_bar=weights, phonons=phonons,
        Eu_J=Eu, N0_per_J_m3=N0, D_m2_s=D, volumes_m3=functional.volumes_m3,
        edges=graph.edges, face_areas_m2=geometry.dual_face_lengths_m*thickness,
        edge_lengths_m=edge_lengths, bath_temperature_K=.9, tau_ee_Tc_s=5e-12,
        tau_escape_s=20e-12, transport_order=4, numba_threads=1)
    gap = 1.65*np.exp(.03j*graph.coordinates_bar[:, 0])
    p0 = functional.fermi_dirac(gap, .9)
    return functional, bath, kinetics, geometry, gap, p0


def test_bath_dc_force_and_current_preserved_with_new_internal_energy(physical_fixture):
    f, bath, _, _, gap, p0 = physical_fixture
    result = f.evaluate(gap, p0, temperature_K=.9)
    old = bath.evaluate(gap, p0)
    np.testing.assert_allclose(result.gap_force_J_per_gapbar, old.gap_force_J_per_gapbar, rtol=2e-14, atol=1e-35)
    np.testing.assert_allclose(result.edge_current_A, old.edge_current_A, rtol=2e-14, atol=1e-25)
    data = f.thermal_data(gap, p0, np.full(len(gap), .9))
    assert np.all(data["capacity_effective_J_K"] > 0)


def test_photon_source_deposits_new_internal_energy_plus_phonons(physical_fixture):
    f, _, k, g, gap, p0 = physical_fixture
    old_gap = gap.copy()
    p, n, temperature, target, receipt = runner.prepare_photon(f, k, gap, p0, g,
        wavelength_nm=775., sigma_nm=7.549834435, retained_fraction=.8, bath_temperature_K=.9)
    n0 = runner.BASE.bose_matrix(k.phonons.energies, f.Eu_J, np.full(len(gap), .9))
    initial = f.evaluate(gap, p0, temperature_K=.9).node_energy_J+np.sum(k.phonon_energy_weights_J*n0, axis=1)
    final = f.evaluate(gap, p, temperature_K=temperature).node_energy_J+np.sum(k.phonon_energy_weights_J*n, axis=1)
    np.testing.assert_allclose(final-initial, target, rtol=1e-9, atol=1e-29)
    assert float(np.sum(final-initial)) == pytest.approx(receipt["retained_energy_J"], rel=1e-10)
    np.testing.assert_array_equal(gap, old_gap)
    np.testing.assert_array_equal(p, f.fermi_dirac(gap, temperature))
    cuts = np.r_[g.left_nodes, g.right_nodes]
    np.testing.assert_array_equal(temperature[cuts], .9)
    assert receipt["cut_deposition_J"] == 0
    assert receipt["minimum_preparation_effective_capacity_ratio"] > 0
    assert receipt["source_must_not_be_applied_again"]
    assert not receipt["candidate_physical_admission"]


def test_nonpositive_capacity_source_fails_without_floor(physical_fixture, monkeypatch):
    f, _, k, g, gap, p0 = physical_fixture
    original = f.thermal_data
    def bad(*args, **kwargs):
        value = original(*args, **kwargs)
        value["capacity_effective_J_K"] = -abs(value["capacity_effective_J_K"])
        return value
    monkeypatch.setattr(f, "thermal_data", bad)
    with pytest.raises(ArithmeticError, match="positive"):
        runner.prepare_photon(f, k, gap, p0, g, wavelength_nm=775., sigma_nm=7.55,
            retained_fraction=.8, bath_temperature_K=.9)


def test_calibrated_geometry_preserves_physical_mesh_and_exact_resume(physical_fixture, monkeypatch, tmp_path):
    f, _, _, old, _, _ = physical_fixture
    monkeypatch.setattr(runner.BASE, "get_effective_geometry", lambda *args, **kwargs: old)
    args = SimpleNamespace(length_nm=160., spacing_nm=40.)
    ell = old.ell0_m*1.3
    result = runner.calibrated_geometry(args, tmp_path, ell)
    np.testing.assert_array_equal(result.coordinates_m, old.coordinates_m)
    np.testing.assert_array_equal(result.graph.edges, old.graph.edges)
    np.testing.assert_array_equal(result.graph.conductance, old.graph.conductance)
    np.testing.assert_allclose(result.graph.area_weights*ell**2,
                               old.graph.area_weights*old.ell0_m**2, rtol=3e-16)
    np.testing.assert_allclose(result.graph.coordinates_bar*ell,
                               old.graph.coordinates_bar*old.ell0_m, rtol=3e-16)
    assert result.receipt["physical_triangulation_unchanged"]
    again = runner.calibrated_geometry(args, tmp_path, ell)
    assert again.receipt["mesh_sha256"] == result.receipt["mesh_sha256"]
    with pytest.raises(ValueError, match="resume"):
        runner.calibrated_geometry(args, tmp_path, ell*1.01)


def test_cli_and_resume_identity_capture_every_new_physical_parameter():
    base = runner.parser().parse_args(["--output-dir", "/home/jdiaz/scratch/candidate_fixture"])
    assert base.sheet_resistance_ohm == 425.
    assert base.target_depairing_uA == pytest.approx(38.8*80/120)
    assert base.tau_escape_ps == 20.
    assert base.sigma_nm == 7.549834435
    assert base.stiffness_temperature == "electron"
    original = runner.BASE.parameter_identity(base)
    for name, value in (("sheet_resistance_ohm", 608.), ("target_depairing_uA", 25.),
                        ("stiffness_temperature", "bath")):
        changed = SimpleNamespace(**vars(base))
        setattr(changed, name, value)
        assert runner.BASE.parameter_identity(changed) != original
    base.current_uA = base.target_depairing_uA
    with pytest.raises(ValueError, match="below"):
        runner.validate_args(base)


def test_manifest_remains_report_compatible_and_marks_candidate(tmp_path, monkeypatch):
    monkeypatch.setattr(runner, "_MODE", "electron")
    manifest = dict(closure={"microscopic_exact_Usadel": False}, admission={"stage5": False})
    target = tmp_path/"manifest.json"
    runner.write_json(target, manifest)
    loaded = runner.BASE.json.loads(target.read_text())
    assert loaded["schema"] == "pysnspd.effective_photon_run.v1"
    assert loaded["model_version"] == "1.1.0-candidate"
    assert loaded["admission"]["candidate"]
    assert not loaded["admission"]["physical_admission"]
    assert not loaded["closure"]["microscopic_exact_Usadel"]
    assert "U_BCS" in loaded["energy_reference"]["u_excess"]
    assert loaded["numerical"]["positive_effective_capacity_required"]