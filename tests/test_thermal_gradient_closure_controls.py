"""Declared gradient controls retain calibrated material and case identity."""
import importlib.util
from pathlib import Path
from types import SimpleNamespace
import numpy as np
import pytest
from scipy.constants import Boltzmann as kB, elementary_charge as e, hbar
from pysnspd.experimental.thermal_spatial_usadel import rectangular_graph

ROOT=Path(__file__).resolve().parents[1]

def _load(name,relative):
    spec=importlib.util.spec_from_file_location(name,ROOT/relative)
    module=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module

runner=_load("thermal_gradient_control_runner","pipelines/09_thermal_stiffness_photon_run.py")
report=_load("thermal_gradient_control_report","plot_pipelines/07_case_report.py")

def _physical_geometry():
    ell=np.sqrt(hbar*5e-5/(2*kB*8.65))
    graph=rectangular_graph(9,17,40e-9/ell,80e-9/ell)
    xy=graph.coordinates_bar*ell
    lengths=np.linalg.norm(xy[graph.edges[:,1]]-xy[graph.edges[:,0]],axis=1)
    left=np.flatnonzero(xy[:,0]==xy[:,0].min())
    right=np.flatnonzero(xy[:,0]==xy[:,0].max())
    return SimpleNamespace(graph=graph,coordinates_m=xy,left_nodes=left,right_nodes=right,
        triangles=np.empty((0,3),dtype=int),edge_lengths_m=lengths,
        dual_face_lengths_m=graph.conductance*lengths,ell0_m=ell,
        source_path="small physical control fixture",receipt={"mesh_sha256":"fixture"})

@pytest.mark.parametrize("gradient",["quadratic","thermal-depairing"])
@pytest.mark.parametrize("mode",["bath","electron"])
def test_small_fresh_DC_build_uses_calibrated_material_for_each_control(tmp_path,monkeypatch,gradient,mode):
    args=runner.parser().parse_args(["--output-dir",str(tmp_path),"--gradient-closure",gradient,
        "--stiffness-temperature",mode,"--dc-only","--length-nm","40","--spacing-nm","5",
        "--xi-cutoff","8","--xi-order","4","--phonon-nodes","24","--transport-order","4",
        "--threads","1"])
    monkeypatch.setattr(runner.BASE,"get_effective_geometry",lambda *a,**k:_physical_geometry())
    system,dc,geometry,material,params=runner.build_system(args,tmp_path)
    f=system.functional
    calibration=params["diffusion_calibration"]
    assert params["gradient_closure"]==gradient
    assert params["diffusion_m2_s"]==calibration["diffusion_m2_s"]
    assert f.ell0_m==calibration["ell0_m"]
    assert f.N0_per_J_m3==calibration["N0_per_J_m3"]
    assert params["N0_single_spin_per_J_m3"]==pytest.approx(
        1/(2*e**2*params["diffusion_m2_s"]*425*7e-9),rel=2e-15)
    expected=("TE_DC_matched_finite_depairing_thermal_gradient_v1" if gradient=="thermal-depairing" else
        "TE_free_gradient_internal_entropy_candidate_v1") if mode=="electron" else None
    if expected is not None:
        assert f.closure_id==expected
        assert calibration["calibration"]["thermal_closure_bulk_verified"]
    assert dc.receipt["current_relative_defect"]<2e-10
    assert dc.receipt["potential_max_abs_V"]<1e-9
    assert (tmp_path/"dc/effective_dc.npz").is_file()
    assert np.all(np.isfinite(dc.delta_bar))

@pytest.mark.parametrize("family",["closure","physicalparams","parameter_identity"])
@pytest.mark.parametrize("gradient",["quadratic","thermal-depairing"])
def test_report_respects_explicit_gradient_identity_in_each_saved_family(family,gradient):
    manifest={"model_version":"1.1.3-candidate","closure":{"stiffness_temperature":"electron"},
        "physicalparams":{},"parameter_identity":{}}
    manifest[family]["gradient_closure"]=gradient
    contract=report.model_contract(manifest)
    assert contract["thermal_depairing"]==(gradient=="thermal-depairing")
    assert ("resummed_thermal.py" in " ".join(report.reconstruction_sources(manifest)))==(gradient=="thermal-depairing")
    assert ("DC-matched" in contract["model_label"])==(gradient=="thermal-depairing")

def test_report_rejects_conflicting_gradient_identity():
    manifest={"model_version":"1.1.3-candidate",
        "closure":{"stiffness_temperature":"electron","gradient_closure":"quadratic"},
        "physicalparams":{"gradient_closure":"thermal-depairing"}}
    with pytest.raises(ValueError,match="identity fields disagree"):
        report.model_contract(manifest)

def test_report_rejects_unknown_gradient_identity():
    manifest={"model_version":"1.1.3-candidate",
        "closure":{"stiffness_temperature":"electron","gradient_closure":"unspecified"}}
    with pytest.raises(ValueError,match="Unsupported gradient_closure"):
        report.model_contract(manifest)

def test_early_113_identity_and_historical_bath_controls_remain_supported():
    manifest={"model_version":"1.1.3-candidate","closure":{"stiffness_temperature":"electron"}}
    assert report.thermal_depairing_candidate(manifest)
    for version in report.THERMAL_CANDIDATE_VERSIONS:
        manifest={"model_version":version,"closure":{"stiffness_temperature":"bath","gradient_closure":"thermal-depairing"}}
        assert not report.thermal_depairing_candidate(manifest)
        assert report.model_contract(manifest)["bath_control"]
        assert report.reconstruction_sources(manifest)==report.RECON_SOURCES


@pytest.mark.parametrize("gradient",["quadratic","thermal-depairing"])
@pytest.mark.parametrize("mode",["bath","electron"])
def test_report_reconstructs_declared_control_functional(gradient,mode):
    from pysnspd.experimental.effective_bcs_functional import EffectiveBCSFunctional
    from pysnspd.experimental.thermal_stiffness import ThermalStiffnessFunctional
    from pysnspd.experimental.resummed_thermal import ResummedThermalFunctional
    geometry=_physical_geometry()
    eu=kB*8.65;diffusion=5e-5;ell=geometry.ell0_m
    mesh={"area_weights":geometry.graph.area_weights,"edges":geometry.graph.edges,
        "conductance":geometry.graph.conductance,"coordinates_bar":geometry.graph.coordinates_bar,
        "fixed_nodes":np.r_[geometry.left_nodes,geometry.right_nodes]}
    manifest={"model_version":"1.1.3-candidate",
        "closure":{"stiffness_temperature":mode,"gradient_closure":gradient},
        "parameter_identity":{"xi_cutoff":8.,"xi_order":4,"gradient_closure":gradient},
        "physicalparams":{"energy_scale_J":eu,"N0_single_spin_per_J_m3":1/(2*e**2*diffusion*425*7e-9),
            "film_thickness_m":7e-9,"bath_temperature_K":.9,"diffusion_m2_s":diffusion,
            "sheet_resistance_ohm":425.,"gradient_closure":gradient},
        "geometry":{"ell0_m":ell}}
    functional=report.reconstruct_functional({"m":manifest,"mesh":mesh})
    expected=(ResummedThermalFunctional if gradient=="thermal-depairing" else
        ThermalStiffnessFunctional) if mode=="electron" else EffectiveBCSFunctional
    assert type(functional) is expected
    assert functional.ell0_m==ell
    assert functional.N0_per_J_m3==manifest["physicalparams"]["N0_single_spin_per_J_m3"]


def test_normal_runner_iteration_budget_matches_validated_anchors():
    args=runner.parser().parse_args(["--output-dir","/home/jdiaz/scratch/unused_iteration_test"])
    assert args.max_iterations==60
    assert args.step_tolerance==2e-7
    assert args.local_energy_tolerance==2.5e-6
    assert args.ceiling_ps==100.
