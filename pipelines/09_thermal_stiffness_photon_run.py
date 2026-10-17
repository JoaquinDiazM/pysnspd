"""Candidate electron-temperature rigidity runner; calibrated bath control available.

Reuses the frozen runner's transaction, implicit-midpoint and trigger machinery
through a private imported module. The shared runner retains its transaction and explicit resolution-stop hook.
All numerical products stay in a new scratch run directory. Completion is not
physical admission; positive effective capacity is a required candidate guard.
"""
from __future__ import annotations

from dataclasses import asdict
import importlib.util
import math
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
_spec = importlib.util.spec_from_file_location("_thermal_candidate_frozen_runner", ROOT/"pipelines/05_effective_photon_run.py")
BASE = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(BASE)
import numpy as np
from scipy.constants import Boltzmann as kB, elementary_charge as e, hbar, Planck, speed_of_light
from scipy.special import zeta
from pysnspd.experimental.depairing_calibration import calibration_result
from pysnspd.experimental.effective_bcs_preparation import load_dual_geometry

_ORIGINAL_PARSER = BASE.parser
_ORIGINAL_SOURCE_PINS = BASE.source_pins
_ORIGINAL_WRITE_JSON = BASE.write_json
_ORIGINAL_PREPARE_PHOTON = BASE.prepare_photon
_ORIGINAL_VALIDATE_ARGS = BASE.validate_args
_MODE = "electron"
_LAST_SYSTEM = None
SCHEMA = "pysnspd.effective_photon_run.v1"


def parser():
    result = _ORIGINAL_PARSER()
    result.description = __doc__
    result.add_argument("--gradient-closure", choices=("thermal-depairing","quadratic"), default="thermal-depairing",
        help="DC-matched finite-depairing thermal gradient; quadratic is a historical control")
    result.add_argument("--sheet-resistance-ohm", type=float, default=425.)
    result.add_argument("--target-depairing-uA", type=float, default=38.8*80/120)
    result.add_argument("--stiffness-temperature", choices=("electron", "bath"), default="electron")
    result.add_argument("--minimum-adaptive-step-ps", type=float, default=1e-6,
        help="Stop and preserve checkpoint below this adaptive resolution; never floor a timestep")
    result.set_defaults(tau_escape_ps=20., current_uA=21.5, local_energy_tolerance=2.5e-6, ceiling_ps=100., max_iterations=60)
    return result


def validate_args(args):
    _ORIGINAL_VALIDATE_ARGS(args)
    for name in ("sheet_resistance_ohm", "target_depairing_uA"):
        value = float(getattr(args, name))
        if not math.isfinite(value) or value <= 0:
            raise ValueError(name+" must be finite and positive")
    if not math.isfinite(args.minimum_adaptive_step_ps) or not 0<args.minimum_adaptive_step_ps<=args.dt_ps:
        raise ValueError("minimum_adaptive_step_ps must be positive and no larger than dt_ps")
    if args.current_uA >= args.target_depairing_uA:
        raise ValueError("current_uA must remain below target_depairing_uA")
    if args.bath_K >= 8.65:
        raise ValueError("Candidate bath temperature must be below Tc=8.65 K")
    return args


def source_pins():
    result = _ORIGINAL_SOURCE_PINS()
    paths = [Path(__file__), ROOT/"pysnspd/experimental/depairing_calibration.py"]
    if _MODE == "electron":
        paths += [ROOT/"pysnspd/experimental/thermal_stiffness.py",
                  ROOT/"pysnspd/experimental/thermal_coupled.py",
                  ROOT/"pysnspd/experimental/general_thermal_gradient.py",
                  ROOT/"pysnspd/experimental/finite_depairing_sums.py",
                  ROOT/"pysnspd/experimental/resummed_thermal.py"]
    if _MODE=="bath":
        paths += [ROOT/"pysnspd/experimental/general_thermal_gradient.py",
                  ROOT/"pysnspd/experimental/finite_depairing_sums.py",
                  ROOT/"pysnspd/experimental/resummed_thermal.py",
                  ROOT/"pysnspd/experimental/thermal_stiffness.py"]
    for path in paths:
        if not path.is_file():
            raise FileNotFoundError("Candidate implementation is incomplete: "+str(path))
        result[str(path.relative_to(ROOT))] = BASE.sha256(path)
    return result


def write_json(path, value):
    if Path(path).name == "manifest.json" and isinstance(value, dict):
        value["schema"] = SCHEMA
        value["model_version"] = "1.1.3-candidate"
        if value.get("status")=="INTEGRATION_FAILURE" and "Declared adaptive-step resolution reached" in value.get("error",""):
            value["status"]="NUMERICAL_RESOLUTION_LIMIT"
            value["resolution_limit_is_detection"]=False
        value["candidate_schema"] = "pysnspd.thermal_stiffness_photon_candidate.v1"
        value.setdefault("admission", {}).update(candidate=True, stage5=False,
            physical_admission=False, policy="Candidate execution is not numerical or physical admission")
        closure = value.setdefault("closure", {})
        closure.update(name="effective neutral thermal-rigidity BCS/KWT/native-kinetic/CM candidate",
            stiffness_temperature=_MODE,
            gradient_rigidity=("G(Delta,Te); Ugrad=G-sum(Te*dG/dTe); positive Ceff guard" if _MODE == "electron"
                              else "Calibrated-D bath-temperature rigidity control; fixed Tb"),
            calibrated_experimental_model=False,
            diffusion_calibration="Homogeneous depairing-target fit; not independent experimental validation",
            signed_heat_carrier=("f*(1-p) heating / p*(1-f) cooling; energy moment and entropy inequality; no clipping" if _MODE == "electron" else "Inherited fixed-bath dissipative heat carrier; no Te reversible-gradient source"))
        value.setdefault("energy_reference", {}).update(
            u_excess=("Local internal U_BCS+G-sum(Te*G_Te)+phonons minus fresh DC per physical nodal volume"
                      if _MODE == "electron" else "Local BCS+fixed-bath gradient+phonons minus fresh DC per nodal volume"))
        selected=value.get("parameter_identity",{}).get("gradient_closure","thermal-depairing")
        closure.update(gradient_closure=selected,
            finite_depairing_thermal_completion=selected=="thermal-depairing",
            matching_reference="Frozen bath/DC reference; no free latency parameter",
            free_gradient=("G=sum[kappa(rho,Tb)*L+u*(phi(rho,Te/Tc,Qbar)-phi(rho,Tb/Tc,Qbar))], Qbar=L/(u*rho)"
                           if selected=="thermal-depairing" and _MODE=="electron" else
                           "Frozen validated bath gradient G=sum kappa(rho,Tb)*L" if _MODE=="bath" else
                           "Historical quadratic G=sum L*kappa(rho,Te)"),
            native_kinetic_spectrum="Local q=0 BCS; not a microscopic finite-q Usadel kinetic evolution",
            effective_Q_includes_amplitude_gradients=selected=="thermal-depairing",
            persistent_elastic_matching_background=(selected=="thermal-depairing" and _MODE=="electron"),
            matching_background_limitation=("H=kappa(rho,Tb)*L-u*phi(rho,Tb/Tc,Q) is independent of evolved Te and persists at fixed rho,L as Te grows; this is a strong mesoscopic hypothesis, not exact microscopic high-temperature stiffness"
                if selected=="thermal-depairing" and _MODE=="electron" else None))
        if _MODE == "electron":
            value.setdefault("numerical", {}).update(positive_effective_capacity_required=True,
                capacity_floor_used=False, signed_reversible_gradient_source=True,
                thermal_induced_gap_guard="Full thermal response checked inside every nonlinear convergence test and at returned endpoint",
                nonlinear_tolerance_policy="nominal tolerance * min(1, h_ps/0.02); at least one coupled correction",
                internal_capacity_line_search=True,
                full_coupled_anderson=True,
                elliptic_gradient_required=selected=="thermal-depairing",
                stiffness_floor_used=False,
                adaptive_resolution_stop_ps=value.get("parameter_identity",{}).get("minimum_adaptive_step_ps",1e-6))
            if _LAST_SYSTEM is not None and hasattr(_LAST_SYSTEM, "last_thermal_diagnostics"):
                value["thermal_candidate_diagnostics"] = dict(_LAST_SYSTEM.last_thermal_diagnostics)
                if hasattr(_LAST_SYSTEM, "minimum_evaluated_capacity_ratio"):
                    value["thermal_candidate_diagnostics"].update(
                        min_successful_rhs_evaluations=float(_LAST_SYSTEM.minimum_evaluated_capacity_ratio),
                        capacity_minimum_scope="All successful RHS evaluations, including rejected positive-capacity trials; not accepted-only")
    return _ORIGINAL_WRITE_JSON(path, value)


def calibrated_geometry(args, output, ell0_m):
    """Keep physical triangulation/areas/faces; rebuild dimensionless coordinates."""
    old = BASE.get_effective_geometry(args.length_nm*1e-9, args.spacing_nm*1e-9,
                                     output_root=output/"meshes/reference_geometry")
    physical_area = old.graph.area_weights*old.ell0_m**2
    normalized_coordinates = old.graph.coordinates_bar*(old.ell0_m/ell0_m)
    directory = output/"meshes/calibrated_scale"
    path = directory/"mesh.npz"
    fixed = old.graph.boundary_nodes
    if fixed is None:
        fixed = np.unique(np.r_[old.left_nodes, old.right_nodes])
    arrays = dict(coordinates_m=old.coordinates_m, coordinates_bar=normalized_coordinates,
        triangles=old.triangles, edges=old.graph.edges,
        area_weights=physical_area/ell0_m**2, conductance=old.graph.conductance,
        dual_area_m2=physical_area, edge_length_m=old.edge_lengths_m,
        dual_face_length_m=old.dual_face_lengths_m, fixed_nodes=fixed,
        left_nodes=old.left_nodes, right_nodes=old.right_nodes, ell0_m=np.array(ell0_m))
    if not path.exists():
        directory.mkdir(parents=True, exist_ok=True)
        BASE.atomic_npz(path, arrays)
    else:
        with np.load(path, allow_pickle=False) as saved:
            if set(saved.files) != set(arrays) or any(not np.array_equal(saved[name], array)
                                                     for name, array in arrays.items()):
                raise ValueError("Saved calibrated geometry differs; exact resume refused")
    geometry = load_dual_geometry(path)
    geometry.receipt.update(reference_mesh_sha256=old.receipt["mesh_sha256"],
        reference_mesh_path=old.source_path, reference_ell0_m=old.ell0_m,
        normalization_rebuilt=True, physical_triangulation_unchanged=True,
        physical_dual_areas_unchanged=True, physical_edge_lengths_unchanged=True,
        calibrated_ell0_m=ell0_m, new_physical_mesh_generated=False)
    return geometry


def build_system(args, output, *, resume=False):
    global _LAST_SYSTEM
    Tc, thickness = 8.65, 7e-9
    from pysnspd.experimental.resummed_thermal import ResummedThermalFunctional
    finite_depairing=getattr(args,"gradient_closure","thermal-depairing")=="thermal-depairing"
    calibration = calibration_result(bulk_factory=ResummedThermalFunctional if finite_depairing else None,
        sheet_resistance_ohm=args.sheet_resistance_ohm,
        target_depairing_current_A=(None if math.isclose(args.target_depairing_uA, 38.8*80/120, rel_tol=8*np.finfo(float).eps) else args.target_depairing_uA*1e-6),
        width_m=80e-9, requested_current_A=args.current_uA*1e-6,
        bath_temperature_K=args.bath_K, critical_temperature_K=Tc,
        count_cutoff_bar=args.xi_cutoff, count_order=args.xi_order)
    if finite_depairing:
        calibration["sources"]["closure"]="pysnspd.experimental.resummed_thermal.ResummedThermalFunctional; finite native BCS plus DC-matched equilibrium finite-depairing thermal correction"
        calibration["sources"]["finite_depairing"]="https://arxiv.org/html/1207.6421"
        calibration["assumptions"]["future_closure"]="Actual finite-depairing bulk callback rebuilt and re-evaluated; no old-closure D reused"
    D, N0, ell = (calibration[name] for name in ("diffusion_m2_s", "N0_per_J_m3", "ell0_m"))
    Eu, Rs = kB*Tc, args.sheet_resistance_ohm
    geometry = calibrated_geometry(args, output, ell)
    xi, weights = BASE.bcs_count_quadrature(args.xi_cutoff, args.xi_order)
    kwargs = dict(graph=geometry.graph, energy_scale_J=Eu, N0_per_J_m3=N0,
        ell0_m=ell, thickness_m=thickness, count_nodes=xi, count_weights=weights,
        gradient_law="bcs_radial", bath_temperature_K=args.bath_K, diffusion_m2_s=D)
    bath_functional = BASE.EffectiveBCSFunctional(**kwargs)
    if args.stiffness_temperature == "electron":
        from pysnspd.experimental.thermal_stiffness import ThermalStiffnessFunctional
        from pysnspd.experimental.thermal_coupled import ThermalCoupledSystem
        functional = (ResummedThermalFunctional(**kwargs) if finite_depairing else ThermalStiffnessFunctional(**kwargs))
        System = ThermalCoupledSystem
        # Check the actual new closure at the calibrated parameters; inheritance
        # or a shared name is not sufficient evidence of the same DC functional.
        actual_bulk = functional.bulk(temperature_K=args.bath_K, current_A=0., width_m=geometry.width_m)
        residual = actual_bulk.depairing_current_A/calibration["target_depairing_current_A"]-1
        if not math.isfinite(residual) or abs(residual) > 2e-8:
            raise ArithmeticError("Thermal closure depairing differs from calibration; recalibrate that closure")
        calibration["calibration"].update(requires_verification_with_new_closure=False,
            thermal_closure_bulk_verified=True, thermal_closure_relative_residual=residual)
    else:
        functional, System = bath_functional, BASE.EffectiveCoupledSystem
    dcpath = output/"dc/effective_dc.npz"
    if resume:
        receipt = BASE.json.loads((output/"dc/effective_dc_receipt.json").read_text())
        if BASE.sha256(dcpath) != receipt["fields_sha256"] or receipt["geometry"]["mesh_sha256"] != geometry.receipt["mesh_sha256"]:
            raise ValueError("Saved DC or calibrated geometry changed")
        with np.load(dcpath, allow_pickle=False) as data:
            dc = BASE.EffectiveBCSDC(data["delta_bar"].copy(), data["p"].copy(), data["phi_V"].copy(),
                BASE.EffectiveBCSBulk(**receipt["bulk"]), float(receipt["phase_left_rad"]), receipt)
    else:
        dc = BASE.prepare_effective_dc(bath_functional, left_nodes=geometry.left_nodes,
            right_nodes=geometry.right_nodes, temperature_K=args.bath_K,
            current_A=args.current_uA*1e-6, width_m=geometry.width_m, sheet_resistance_ohm=Rs)
        BASE.save_effective_dc(output/"dc", dc, geometry)
    if args.stiffness_temperature == "electron":
        temperatures = np.full(geometry.graph.n_nodes, args.bath_K)
        fresh = functional.evaluate(dc.delta_bar, dc.p, temperature_K=temperatures)
        old = bath_functional.evaluate(dc.delta_bar, dc.p)
        if not np.allclose(fresh.gap_force_J_per_gapbar, old.gap_force_J_per_gapbar,
                           rtol=2e-8, atol=1e-10*np.max(functional.node_energy_scale_J)):
            raise ArithmeticError("Thermal closure does not retain the prepared bath DC force")
        thermal = functional.thermal_data(dc.delta_bar, dc.p, temperatures)
        ce = thermal["capacity_effective_J_K"]
        if np.any(~np.isfinite(ce)) or np.any(ce <= 0):
            raise ArithmeticError("Prepared candidate DC has nonpositive effective capacity; no floor allowed")
    tau0 = 720*float(zeta(5))/np.pi**2*args.tau_ep_Tc_ps*1e-12
    material = BASE.allmaras_acoustic_material(Tc_K=Tc, tau0_s=tau0, diffusivity_m2_s=D,
        sheet_resistance_ohm=Rs, thickness_m=thickness, nodes=1025)
    if not resume:
        material.save(output/"material.npz", output/"material.json")
    phonons = BASE.acoustic_phonon_grid(material, Eu, nodes=args.phonon_nodes)
    kinetics = BASE.EffectiveBCSKinetics(xi_bar=xi, weights_bar=weights, phonons=phonons, Eu_J=Eu,
        N0_per_J_m3=N0, D_m2_s=D, volumes_m3=functional.volumes_m3, edges=geometry.graph.edges,
        face_areas_m2=geometry.face_areas_m2(thickness), edge_lengths_m=geometry.edge_lengths_m,
        bath_temperature_K=args.bath_K, tau_ee_Tc_s=args.tau_ee_Tc_ps*1e-12,
        tau_escape_s=args.tau_escape_ps*1e-12, transport_order=args.transport_order,
        numba_threads=args.threads)
    ports = BASE.EffectiveNeutralPorts(geometry.graph, sheet_resistance_ohm=Rs,
        left_nodes=geometry.left_nodes, right_nodes=geometry.right_nodes)
    exterior = BASE.freeze_exterior_inductance(hbar/(2*e)*dc.bulk.q_derivative_m_inv_per_A, geometry.length_m)
    circuit = BASE.ThesisCircuitParameters(Lk_ext_H=exterior.total_H, V_bias_V=1e4*args.current_uA*1e-6)
    system = System(functional, kinetics, ports, circuit, tau_ee_Tc_ps=args.tau_ee_Tc_ps,
        tau_ep_Tc_ps=args.tau_ep_Tc_ps, current_scale_A=args.current_uA*1e-6,
        kinetic_residual_energy_scale_J=None if args.dc_only else
            Planck*speed_of_light/(args.wavelength_nm*1e-9)*args.retained_fraction)
    _LAST_SYSTEM = system
    params = dict(Tc_K=Tc, bath_temperature_K=args.bath_K, diffusion_m2_s=D,
        sheet_resistance_ohm=Rs, film_thickness_m=thickness, N0_single_spin_per_J_m3=N0,
        ell0_m=ell, energy_scale_J=Eu, current_A=args.current_uA*1e-6,
        tau_ee_Tc_s=args.tau_ee_Tc_ps*1e-12, tau_ep_Tc_s=args.tau_ep_Tc_ps*1e-12,
        tau0_s=tau0, tau_escape_s=args.tau_escape_ps*1e-12,
        stiffness_temperature=args.stiffness_temperature, diffusion_calibration=calibration,
        gradient_closure="thermal-depairing" if finite_depairing else "quadratic",
        gradient_hypothesis=("Qbar=L/(u*rho) dimensionless; Qeff_J=Eu*Qbar; theta=Te/Tc; G=kappa(rho,Tb)*L+u*(phi(rho,theta,Qbar)-phi(rho,Tb/Tc,Qbar)); native q=0 kinetics retained" if finite_depairing else "Historical quadratic small-gradient free energy"),
        persistent_elastic_matching_background=(finite_depairing and args.stiffness_temperature=="electron"),
        matching_background_limitation=("Own temperature-independent H hypothesis can dominate hot strong-gradient states; no quantitative latency admission"
            if finite_depairing and args.stiffness_temperature=="electron" else None),
        calibrators="Free effective hypotheses; depairing-target fit is not independent material measurement",
        exterior_inductance=asdict(exterior), circuit=asdict(circuit))
    physics_identity = dict(schema="pysnspd.thermal_candidate_physics_identity.v1",
        model_version="1.1.3-candidate", source_pins_sha256=source_pins(),
        parameter_identity=BASE.parameter_identity(args), calibration=calibration,
        calibrated_mesh_sha256=geometry.receipt["mesh_sha256"],
        candidate_physical_admission=False)
    for path, content in ((output/"calibration.json", calibration),
                          (output/"candidate_physics_identity.json", physics_identity)):
        if path.exists():
            if BASE.json.loads(path.read_text()) != content:
                raise ValueError("Saved candidate calibration/physics identity changed")
        else:
            _ORIGINAL_WRITE_JSON(path, content)
    return system, dc, geometry, material, params


def prepare_photon(functional, kinetics, gap, p0, geometry, *, wavelength_nm,
                   sigma_nm, retained_fraction, bath_temperature_K):
    if not hasattr(functional, "thermal_data"):
        result = _ORIGINAL_PREPARE_PHOTON(functional, kinetics, gap, p0, geometry,
            wavelength_nm=wavelength_nm, sigma_nm=sigma_nm,
            retained_fraction=retained_fraction, bath_temperature_K=bath_temperature_K)
        result[-1].update(preparation_width_independent_of_calibrated_D=True,
            energy_functional="BCS+fixed-bath gradient+phonons; calibrated-D control")
        return result
    count = functional.graph.n_nodes
    bath = np.full(count, bath_temperature_K)
    n0 = BASE.bose_matrix(kinetics.phonons.energies, functional.Eu_J, bath)
    baseline = functional.evaluate(gap, p0, temperature_K=bath).node_energy_J
    baseline = baseline+np.sum(kinetics.phonon_energy_weights_J*n0, axis=1)
    center = np.array([.5*(geometry.coordinates_m[:, 0].min()+geometry.coordinates_m[:, 0].max()), 0.])
    shape = np.exp(-np.sum((geometry.coordinates_m-center)**2, axis=1)/(2*(sigma_nm*1e-9)**2))
    cuts = np.r_[geometry.left_nodes, geometry.right_nodes]
    shape[cuts] = 0.
    weighted = shape*functional.volumes_m3
    if not np.isfinite(weighted.sum()) or weighted.sum() <= 0:
        raise ValueError("Photon profile lacks finite interior support")
    photon = Planck*speed_of_light/(wavelength_nm*1e-9)
    retained = photon*retained_fraction
    target = retained*weighted/weighted.sum()
    active = target > 0
    minimum_ratio = math.inf

    def values(temperature):
        nonlocal minimum_ratio
        p = functional.fermi_dirac(gap, temperature)
        n = BASE.bose_matrix(kinetics.phonons.energies, functional.Eu_J, temperature)
        thermal = functional.thermal_data(gap, p, temperature)
        ce, cfd = thermal["capacity_effective_J_K"], thermal["capacity_FD_J_K"]
        cph = np.sum(kinetics.phonon_energy_weights_J*n*(n+1)*
            (kinetics.phonons.energies[None, :]*functional.Eu_J)/(kB*temperature[:, None]**2), axis=1)
        if np.any(~np.isfinite(ce+cph)) or np.any(ce <= 0) or np.any(cfd <= 0) or np.any(ce+cph <= 0):
            raise ArithmeticError("Photon preparation left positive Ceff/Ctotal domain; no capacity floor")
        minimum_ratio = min(minimum_ratio, float(np.min(ce/cfd)))
        energy = functional.evaluate(gap, p, temperature_K=temperature).node_energy_J
        energy = energy+np.sum(kinetics.phonon_energy_weights_J*n, axis=1)
        return energy-baseline, p, n

    values(bath)
    lo = bath.copy()
    hi = np.where(active, np.maximum(2*bath, 2.), bath)
    for _ in range(32):
        added, _, _ = values(hi)
        need = (added < target)&active
        if not np.any(need):
            break
        hi[need] *= 2
    else:
        raise ArithmeticError("No positive-capacity photon internal-energy bracket")
    for _ in range(60):
        middle = (lo+hi)/2
        added, _, _ = values(middle)
        below = added < target
        lo = np.where(below&active, middle, lo)
        hi = np.where((~below)&active, middle, hi)
    temperature = np.where(active, (lo+hi)/2, bath)
    added, p, n = values(temperature)
    defect = added-target
    if abs(float(defect.sum())) > 1e-10*retained or float(np.max(abs(defect))) > 1e-10*retained:
        raise ArithmeticError("New thermal photon preparation failed retained internal-energy balance")
    receipt = dict(wavelength_nm=wavelength_nm, photon_energy_J=photon,
        retained_fraction=retained_fraction, retained_energy_J=retained,
        sigma_nm=sigma_nm, center_m=center.tolist(),
        initialization="Common FD/BE post-cascade transfer; fresh calibrated DC gap unchanged",
        energy_functional="U_BCS+G-sum(Te*dG/dTe)+phonons; fixed DC gap during transfer",
        time_origin="Effective post-cascade transfer, no extra cascade delay",
        cascade_reference="Independent exploratory width; unchanged by calibrated diffusion",
        preparation_width_independent_of_calibrated_D=True,
        energy_deposited_J=float(added.sum()), energy_defect_J=float(defect.sum()),
        maximum_nodal_energy_defect_J=float(np.max(abs(defect))),
        relative_total_energy_defect=float(abs(defect.sum())/retained),
        maximum_initial_temperature_K=float(temperature.max()),
        minimum_preparation_effective_capacity_ratio=minimum_ratio,
        positive_capacity_checked_throughout_preparation=True, capacity_floor_used=False,
        cut_deposition_J=float(target[cuts].sum()), source_must_not_be_applied_again=True,
        old_prepared_photon_state_reused=False, candidate_physical_admission=False)
    return p, n, temperature, target, receipt


BASE.build_system = build_system
BASE.prepare_photon = prepare_photon
BASE.source_pins = source_pins
BASE.write_json = write_json
BASE.validate_args = validate_args
BASE.parser = parser
BASE.SCHEMA = SCHEMA


def main(argv=None):
    global _MODE
    args = parser().parse_args(argv)
    _MODE = args.stiffness_temperature
    requested = {"signal": None}
    def terminate(signum, frame):
        requested["signal"] = int(signum)
    previous = BASE.signal.signal(BASE.signal.SIGTERM, terminate)
    try:
        return BASE.run(args, stop_requested=lambda: requested["signal"] is not None)
    finally:
        BASE.signal.signal(BASE.signal.SIGTERM, previous)


if __name__ == "__main__":
    main()