"""Registered mixed-field and instantaneous power diagnostic; no time evolution.

Default describes the frozen experiment without evaluating a spectrum. --pilot
runs the two profiles on the smallest mesh; --execute runs all six snapshots.
Every output directory is new. Failure is preserved and stops the sequence.
Energy patches are disabled: rejected interpolants cannot be selected silently.
"""
from pathlib import Path
import argparse
from dataclasses import asdict
from datetime import datetime, timezone
import hashlib
import json
import platform
import sys
import time
import traceback

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from progress import Progress
from pysnspd.experimental.energy_catalog import OccupationEnergyCatalog, K_B_J_K
from pysnspd.experimental.refined_cells import refined_count_catalog
from pysnspd.experimental.local_energy_patch import LocalEnergyPatch
from pysnspd.experimental.seeded_count_catalog import SeededCountCatalog
from pysnspd.experimental.mixed_spatial import MixedSpatialFunctional
from pysnspd.experimental.cell_closures import CellScales, KWTMobility
from pysnspd.experimental.cell_validation import ElectronicCell
from pysnspd.experimental.electrical_ports import solve_potential, ThesisCircuitParameters
from pysnspd.experimental.spatial_dynamics import (
    equivalent_temperatures, kwt_spatial_response, gauge_power_check,
    normal_heat_distribution, deposit_spatial_heat, spatial_energy_rate, cm9_power_balance,
)


def plain(value):
    if isinstance(value, np.ndarray):
        return plain(value.tolist())
    if isinstance(value, np.generic):
        return plain(value.item())
    if isinstance(value, complex):
        return dict(real=value.real, imag=value.imag)
    if isinstance(value, dict):
        return {k: plain(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [plain(v) for v in value]
    return value


def save(path, value):
    Path(path).write_text(json.dumps(plain(value), indent=2, ensure_ascii=False,
                                    allow_nan=False)+"\n", encoding="utf-8")


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def sources(registration):
    modules = ("energy_catalog", "refined_cells", "spatial_functional", "spatial_open",
               "mixed_spatial", "cell_closures", "cell_validation", "electrical_ports", "spatial_dynamics",
               "local_energy_patch", "seeded_count_catalog")
    paths = [ROOT/("pysnspd/experimental/"+name+".py") for name in modules]
    paths += [Path(__file__).resolve(), ROOT/"sandbox/stage3_spatial/progress.py", registration,
              ROOT/"docs/modelo_v0_4/D_sintesis_plan_y_verificaciones_v0_4.md",
              ROOT/"docs/modelo_v0_4/actualizaciones/circuito_memoria_20260922.md"]
    reg = json.loads(registration.read_text(encoding="utf-8"))
    paths += [ROOT/reg["catalogue"], ROOT/reg["reference_reservoir"]["path"]]
    return {str(p.resolve().relative_to(ROOT)).replace("\\", "/"): sha(p) for p in paths}


def validate_registration(reg):
    if sha(ROOT/reg["catalogue"]) != reg["catalogue_sha256"]:
        raise ValueError("registered catalogue changed")
    reference = reg["reference_reservoir"]
    if sha(ROOT/reference["path"]) != reference["sha256"]:
        raise ValueError("registered reservoir reference changed")
    reservoir = json.loads((ROOT/reference["path"]).read_text(encoding="utf-8"))
    partition = reservoir["inductance_partition"]
    expected = {"exterior_fixed_H": reg["circuit"]["Lk_ext_H"],
                "reference_current_A": reg["circuit"]["reference_current_A"],
                "reference_amplitude_bar": reg["field"]["amplitude_base"],
                "reference_q_bare_bar": reg["field"]["q_bare_bar"], "bath_theta": reg["bath_theta"]}
    if any(partition[k] != v for k, v in expected.items()):
        raise ValueError("registered circuit/field reference differs from the archived partition")


def verify_seed_directory(folder, mixed_reg):
    """Read-only acceptance/provenance checks before enabling exact Newton seeding."""
    folder = (ROOT/folder).resolve()
    manifest_path, summary_path = folder/"manifest.json", folder/"summary.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    reg = manifest["registration"]
    if (summary.get("status") != "PASS_SEEDED_CAUSAL_CATALOG"
            or summary.get("representation") != "predictor_then_exact_causal_Newton"
            or reg.get("schema") != "pysnspd.seeded-causal-catalog.checks.v3"
            or reg["catalogue"] != mixed_reg["catalogue"] or reg["refinement"] != mixed_reg["refinement"]):
        raise ValueError("seeded catalogue has no compatible exact-causal PASS")
    expected_sources = {"pysnspd/experimental/"+name+".py" for name in
                        ("local_energy_patch", "seeded_count_catalog", "refined_cells", "energy_catalog", "spatial_functional")}
    registration_path = "docs/implementation/stage3/coupled_20260923/seeded_registration.json"
    expected_sources.update(("sandbox/stage3_spatial/coupled_20260923/check_seeded_catalog.py",
                             registration_path, reg["catalogue"], reg["predictor"]))
    if set(manifest["sources"]) != expected_sources:
        raise ValueError("seeded manifest source coverage is incomplete or unexpected")
    for path, expected in manifest["sources"].items():
        if sha(ROOT/path) != expected:
            raise ValueError("seeded source mismatch: "+path)
    if json.loads((ROOT/registration_path).read_text(encoding="utf-8")) != reg:
        raise ValueError("seeded registration contents mismatch")
    predictor_manifest = ROOT/reg["predictor_manifest"]
    previous = json.loads(predictor_manifest.read_text(encoding="utf-8"))
    for path, expected in previous["sources"].items():
        if sha(ROOT/path) != expected:
            raise ValueError("predictor provenance changed: "+path)
    degree = reg["selected_degree"]
    if summary["selected_degree"] != degree or degree not in reg["degrees"]:
        raise ValueError("unregistered seeded candidate")
    artifact = folder/f"patch_degree{degree}.npz"
    if sha(artifact) != summary["patch_sha256"] or sha(artifact) != reg["predictor_sha256"]:
        raise ValueError("seed predictor artifact changed")
    probes_path = folder/"probes.json"
    probes = json.loads(probes_path.read_text(encoding="utf-8"))
    if probes != summary["probes"] or len(probes) != len(reg["fractional_probes"]):
        raise ValueError("incomplete or mismatched seeded probe coverage")
    for row, fractions in zip(probes, reg["fractional_probes"]):
        expected_a = reg["amplitude_bounds"][0]+fractions[0]*np.diff(reg["amplitude_bounds"])[0]
        expected_g = reg["gamma_bounds"][0]+fractions[1]*np.diff(reg["gamma_bounds"])[0]
        candidate = row["degrees"][str(degree)]
        if row["amplitude"] != expected_a or row["gamma"] != expected_g:
            raise ValueError("seeded probe coordinates mismatch")
        if (candidate["passed"] is not True or candidate["sign_margin"] <= 0
                or not np.isfinite(candidate["sign_margin"])
                or set(candidate["metrics"]) != set(reg["limits"])
                or any(not np.isfinite(candidate["metrics"][key]) or candidate["metrics"][key] < 0
                       or candidate["metrics"][key] > limit for key, limit in reg["limits"].items())):
            raise ValueError("seeded probe verdict is not supported by its metrics")
    return dict(artifact=str(artifact), limits=reg["limits"], selected_degree=degree,
        representation=summary["representation"],
        verified_files={str(path): sha(path) for path in (manifest_path, summary_path, probes_path, artifact, predictor_manifest)})


def direct_reference_probes(direct, model, baseline, z, populations, limits, progress):
    """Three labelled local comparisons, never a whole-field bias certificate."""
    fields = baseline.sampled_fields
    selected = [("minimum_amplitude", int(np.argmin(fields.amplitude_quadrature_bar))),
                ("maximum_Gamma", int(np.argmax(fields.gamma_quadrature_bar))),
                ("maximum_transverse_q", int(np.argmax(abs(fields.q_delta_quadrature_bar[:, 1]))))]
    old_catalog = model.catalog
    records = []
    for reason, index in selected:
        progress.checkpoint("contraste causal directo: "+reason, force=True)
        a, gamma = fields.amplitude_quadrature_bar[index], fields.gamma_quadrature_bar[index]
        p = populations[index]
        exact = np.asarray(direct.energy_kernel(a, gamma))
        corrected = np.asarray(old_catalog.energy_kernel(a, gamma))
        difference = corrected-exact
        w = direct.count_weights
        kernel = np.sum(abs(difference)*w, axis=1)/np.sum(abs(exact)*w, axis=1)
        moment = abs(difference@(w*p))/np.sum(abs(exact)*w*p, axis=1)
        dimensions = int(model.spatial_dimensions[index])
        # Only the local principal-symbol catalogue changes, at identical z,d,p.
        try:
            model.catalog = direct
            reference = model.principal_symbol(fields.delta_quadrature_bar[index],
                fields.derivative_quadrature_bar[index, :dimensions], p, spatial_dimensions=dimensions)
        finally:
            model.catalog = old_catalog
        measured = baseline.principal_symbols[index]
        difference_norm = np.linalg.norm(measured.matrix-reference.matrix, 2)
        margin = min(measured.eigenvalues[0], reference.eigenvalues[0])-measured.uncertainty-reference.uncertainty-difference_norm
        metrics = dict(pointwise_relative_energy=float(np.max(abs(difference[0]/exact[0]))),
            weighted_relative_kernel_L1=float(max(kernel)),
            moment_error_over_absolute_weighted_reference=float(max(moment)),
            relative_principal_matrix=float(difference_norm/np.linalg.norm(reference.matrix, 2)))
        passed = bool(margin > 0 and all(np.isfinite(metrics[k]) and metrics[k] <= v for k, v in limits.items()))
        records.append(dict(reason=reason, quadrature_index=index, amplitude=float(a), gamma=float(gamma),
                            metrics=metrics, sign_margin=float(margin), passed=passed))
    return records


class ExactKernelCache:
    """Per-snapshot exact-key memoization only; no rounding/interpolation/retries."""
    def __init__(self, source):
        self.source = source
        self.cache = {}
        self.evaluations = 0
        self.hits = 0

    def __getattr__(self, name):
        return getattr(self.source, name)

    def energy_kernel(self, amplitude, gamma):
        key = (float(amplitude), float(gamma))
        if key not in self.cache:
            self.cache[key] = tuple(np.array(v, copy=True) for v in self.source.energy_kernel(*key))
            for v in self.cache[key]:
                v.setflags(write=False)
            self.evaluations += 1
        else:
            self.hits += 1
        return self.cache[key]

    def evaluate(self, amplitude, gamma, p):
        return tuple(float(v+4*np.dot(self.count_weights*p, k)) for v, k in
                     zip(self.vacuum.evaluate(amplitude, gamma), self.energy_kernel(amplitude, gamma)))


def heat_allocation(model):
    """Conservative declared geometry rule; no local-accuracy claim."""
    allocation = np.zeros((model.quadrature_size, len(model.graph_edges)))
    for edge, endpoints in enumerate(model.graph_edges):
        for node in endpoints:
            rows = np.flatnonzero(model.quadrature_to_dof == node)
            allocation[rows, edge] += .5*model.quadrature_mass_bar[rows]/model.mass_bar[node]
    return allocation


def solve_case(direct, reg, mesh, profile, out, callback, *, seeded=None, seed_metadata=None):
    started = time.monotonic()
    catalog = ExactKernelCache(direct if seeded is None else seeded)
    geometry = reg["geometry"]
    options = {k: mesh[k] for k in ("elements_x", "elements_y", "left_elements", "right_elements")}
    model = MixedSpatialFunctional(catalog, **geometry, **options, degree=reg["degree"])
    x, y = (model.dof_coordinates_bar*model.ell0_m).T
    inside = (x >= 0)&(x <= geometry["rectangle_length_m"])
    envelope = np.zeros(model.cells)
    envelope[inside] = np.sin(np.pi*x[inside]/geometry["rectangle_length_m"])**2*np.cos(2*np.pi*y[inside]/geometry["width_m"])
    phase = reg["field"]["q_bare_bar"]*x/model.ell0_m+profile["phase_perturbation"]*envelope
    z = (reg["field"]["amplitude_base"]+profile["amplitude_perturbation"]*envelope)*np.exp(1j*phase)
    links = np.zeros(len(model.graph_edges))
    # Deterministic directions are saved before any material evaluation.
    direction = (.2+np.cos(np.arange(model.cells)*.73)+.7j*np.sin(np.arange(model.cells)*.43))*np.exp(1j*phase)
    direction /= max(abs(direction))
    link_direction = .3+.7*np.cos(np.arange(len(links))*.37)
    link_direction /= max(abs(link_direction))
    np.savez_compressed(out/"preparation.npz", delta_bar=z, link_phases=links,
        cartesian_direction=direction, link_direction=link_direction,
        quadrature_to_dof=model.quadrature_to_dof, quadrature_mass_bar=model.quadrature_mass_bar,
        graph_edges=model.graph_edges, dof_coordinates_m=model.dof_coordinates_bar*model.ell0_m,
        quadrature_coordinates_m=model.quadrature_coordinates_bar*model.ell0_m)
    perform_fd = mesh["id"] == reg["fd_case"]["mesh"] and profile["id"] == reg["fd_case"]["profile"]
    progress = Progress(4+(8*len(reg["fd_steps"]) if perform_fd else 0)+(1 if seeded is not None else 0), mesh["id"]+" / "+profile["id"],
                        min_interval=2., callback=callback)
    try:
        progress.start_task("preparación térmica y D.36 en todas las cuadraturas")
        baseline = model.evaluate_thermal(z, reg["bath_theta"], link_phases=links,
            require_stability=True, on_quadrature=lambda i, n: progress.checkpoint(f"cuadratura {i+1}/{n}"))
        p = baseline.p_quadrature
        np.savez_compressed(out/"initial.npz", delta_bar=z, link_phases=links, p_quadrature=p,
                            count_nodes=catalog.count_nodes, count_weights=catalog.count_weights)
        progress.advance()
        fields = baseline.sampled_fields
        direct_probes = []
        if seeded is not None:
            progress.start_task("tres contrastes causales directos a campos y poblaciones idénticos")
            direct_probes = direct_reference_probes(direct, model, baseline, z, p, seed_metadata["limits"], progress)
            save(out/"direct_reference_probes.json", direct_probes)
            if not all(row["passed"] for row in direct_probes):
                raise RuntimeError("exact-corrected catalogue failed direct local comparison; no fallback")
            progress.advance()
        progress.start_task("espectros instantáneos y temperaturas equivalentes")
        cells = []
        for i, (amplitude, gamma) in enumerate(zip(fields.amplitude_quadrature_bar, fields.gamma_quadrature_bar)):
            cells.append(ElectronicCell(catalog, float(amplitude), float(gamma)))
            progress.checkpoint(f"temperatura/espectro {i+1}/{model.quadrature_size}")
        temperature = equivalent_temperatures(cells, p)
        scales = CellScales(catalog.vacuum.delta0_J, catalog.vacuum.N0_per_J_m3,
            model.Tc_K, reg["bath_theta"]*catalog.vacuum.delta0_J/K_B_J_K, reg["time_reference_ps"])
        progress.advance()
        circuit_reg = reg["circuit"]
        reference_current = circuit_reg["reference_current_A"]
        circuit = ThesisCircuitParameters(**{k: circuit_reg[k] for k in
            ("Lk_ext_H", "R_bias_ohm", "L_bias_H", "R_load_ohm", "C_couple_F")},
            V_bias_V=circuit_reg["R_bias_ohm"]*reference_current)
        circuit_state = np.array([circuit_reg["bias_current_fraction"]*reference_current,
            circuit_reg["detector_current_fraction"]*reference_current, circuit_reg["initial_capacitor_voltage_V"]])
        progress.start_task("potencial, KWT y depósito conservativo del calor")
        injection = np.zeros(model.cells)
        left, right = model.terminal_dofs
        injection[left], injection[right] = circuit_state[1], -circuit_state[1]
        conductance = model.conductance(reg["normal_conductivity_S_m"])
        potential = solve_potential(model.graph_edges, conductance, baseline.current_A, injection,
                                    left_node=left, right_node=right)
        kwt = kwt_spatial_response(z, baseline.gradient_cartesian_bar, model.quadrature_mass_bar,
            KWTMobility(scales), temperature, potential.phi_V, quadrature_to_node=model.quadrature_to_dof)
        normal = normal_heat_distribution(conductance, potential.delta_phi_V, heat_allocation(model))
        heating = deposit_spatial_heat(cells, p, kwt, normal.quadrature_power_W,
                                      energy_scale_J=model.energy_scale_J, scales=scales)
        progress.advance()
        progress.start_task("interfaz, Noether y balance instantáneo CM.9")
        gauge = gauge_power_check(z, baseline.gradient_cartesian_bar, model.graph_edges,
            baseline.link_current_bar, potential.phi_V, energy_scale_J=model.energy_scale_J,
            time_scale_s=scales.t_ref_s)
        rate = spatial_energy_rate(baseline.gradient_cartesian_bar, kwt.field_velocity_bar,
                                  cells, heating.electron_rhs, model.quadrature_mass_bar)
        circuit_balance = circuit.power_balance(circuit_state, potential.Vdev_V)
        balance = cm9_power_balance(rate.total_rate_bar, energy_scale_J=model.energy_scale_J,
            time_scale_s=scales.t_ref_s, circuit_balance=circuit_balance)
        progress.advance()
        complex_gradient = baseline.gradient_cartesian_bar[:, 0]+1j*baseline.gradient_cartesian_bar[:, 1]
        derivatives = {}
        checks = (("force", float(np.real(np.vdot(complex_gradient, direction)))),
                  ("current", float(np.dot(baseline.link_current_bar, link_direction)))) if perform_fd else ()
        for kind, predicted in checks:
            estimates = []
            for h in reg["fd_steps"]:
                energy = []
                for multiplier in (-2, -1, 1, 2):
                    progress.start_task(f"FD {kind}: h={h:g}, desplazamiento {multiplier:+d}")
                    value = model.evaluate(z+multiplier*h*direction if kind == "force" else z, p,
                        link_phases=multiplier*h*link_direction if kind == "current" else links,
                        require_stability=False, on_quadrature=lambda i, n: progress.checkpoint(f"cuadratura {i+1}/{n}"))
                    energy.append(value.energy_bar)
                    progress.advance()
                estimates.append((energy[0]-8*energy[1]+8*energy[2]-energy[3])/(12*h))
            tolerance = reg["criteria"]["derivative_absolute"]+reg["criteria"]["derivative_relative"]*abs(estimates[-1])
            error = abs(predicted-estimates[-1])
            uncertainty = abs(estimates[-1]-estimates[-2])
            status = ("INCONCLUSIVE_REFERENCE" if uncertainty > reg["criteria"]["reference_budget_fraction"]*tolerance
                      else "PASS" if error <= tolerance else "FAIL_IDENTITY")
            derivatives[kind] = dict(status=status, predicted=predicted, fd_estimates=estimates,
                                     error=error, reference_uncertainty=uncertainty, tolerance=tolerance)
            save(out/"derivative_checks.json", derivatives)
        derivative_summary = dict(performed=perform_fd,
            status=("PASS" if all(v["status"] == "PASS" for v in derivatives.values()) else
                    "INCONCLUSIVE_REFERENCE" if any(v["status"] == "INCONCLUSIVE_REFERENCE" for v in derivatives.values()) else
                    "FAIL_IDENTITY") if perform_fd else "NOT_REPEATED_IDENTITIES_UNIT_AND_PILOT",
            registered_case=reg["fd_case"], **derivatives)
        save(out/"derivative_checks.json", derivative_summary)
        reactions = baseline.interface_reactions_bar
        interface_force_error = max(abs(reactions["left_lead"]+reactions["left_rectangle"]-complex_gradient[model.interface_dofs["left"]]),
                                    abs(reactions["right_lead"]+reactions["right_rectangle"]-complex_gradient[model.interface_dofs["right"]]))
        trace_error = max(np.max(abs(fields.delta_quadrature_bar[model.rectangle_quadrature[0]]-z[model.interface_dofs["left"]])),
                          np.max(abs(fields.delta_quadrature_bar[model.rectangle_quadrature[-1]]-z[model.interface_dofs["right"]])))
        # Normalize identities by their actual participating rates/powers.
        power_scale = max(abs(circuit_balance.source_power_W), abs(circuit_balance.stored_energy_rate_W),
                          abs(balance.domain_energy_rate_W), potential.normal_joule_W, 1e-30)
        port_scale = max(abs(potential.port_power_W), abs(gauge.field_phase_power_W), potential.normal_joule_W, 1e-30)
        heat_scale = max(abs(heating.deposited_rate_bar), abs(heating.condensate_rate_bar), abs(heating.normal_rate_bar), 1e-30)
        metrics = dict(noether_absolute=float(np.max(abs(gauge.noether_residual_bar))),
            electrical_current_relative=float(np.max(abs(potential.current_residual_A))/reference_current),
            instantaneous_power_relative=float(max(abs(balance.residual_W), abs(circuit_balance.residual_W))/power_scale),
            port_power_relative=float(max(abs(balance.device_port_residual_W), abs(gauge.power_residual_W),
                abs(potential.power_residual_W), abs(potential.decomposition_residual_W))/port_scale),
            heat_moment_relative=float(abs(heating.moment_residual_bar)/heat_scale),
            temperature_absolute=float(np.max(abs(temperature-reg["bath_theta"]))),
            interface_trace_absolute=float(trace_error), interface_force_absolute=float(interface_force_error))
        passed = all(metrics[k] <= reg["criteria"][k] for k in reg["criteria"] if k in metrics)
        passed = passed and metrics["port_power_relative"] <= reg["criteria"]["instantaneous_power_relative"]
        passed = passed and all(v["status"] == "PASS" for v in derivatives.values())
        np.savez_compressed(out/"arrays.npz", delta_bar=z, p_quadrature=p,
            cartesian_gradient_bar=baseline.gradient_cartesian_bar, graph_edges=model.graph_edges,
            phi_V=potential.phi_V, current_super_A=baseline.current_A, current_normal_A=potential.I_normal_A,
            current_total_A=potential.I_total_A, material_velocity_bar=kwt.material_velocity_bar,
            field_velocity_bar=kwt.field_velocity_bar, QDelta_density_bar=kwt.heat_density_bar,
            normal_heat_power_W=normal.quadrature_power_W, electron_heating_rhs=np.array(heating.electron_rhs),
            temperature_bar=temperature, amplitude_quadrature_bar=fields.amplitude_quadrature_bar,
            gamma_quadrature_bar=fields.gamma_quadrature_bar, q_delta_quadrature_bar=fields.q_delta_quadrature_bar,
            principal_minima=np.array([s.eigenvalues[0] for s in baseline.principal_symbols]),
            principal_uncertainties=np.array([s.uncertainty for s in baseline.principal_symbols]))
        record = dict(status="PASS_MIXED_INSTANTANEOUS_IDENTITIES" if passed else "MIXED_INSTANTANEOUS_NOT_ACCEPTED",
            mesh=mesh, profile=profile, scope=reg["scope"], runtime_seconds=time.monotonic()-started,
            dofs=model.cells, quadrature_points=model.quadrature_size, edges=len(model.graph_edges),
            catalogue_mode="DIRECT_CAUSAL_EXACT_KEY_CACHE" if seeded is None else "EXACT_CORRECTED_CAUSAL_EXACT_KEY_CACHE",
            causal_evaluations=catalog.evaluations, cache_hits=catalog.hits, direct_reference_probes=direct_probes,
            metrics=metrics, derivative_checks=derivative_summary, minimum_principal_eigenvalue=min(s.eigenvalues[0] for s in baseline.principal_symbols),
            maximum_principal_uncertainty=max(s.uncertainty for s in baseline.principal_symbols),
            internal_energy_bar=baseline.energy_bar, free_energy_bar=baseline.free_energy_bar,
            gradient_remainder_bar=baseline.gradient_remainder_bar, energy_scale_J=model.energy_scale_J,
            maximum_material_speed_bar=float(max(abs(kwt.material_velocity_bar))),
            maximum_field_speed_bar=float(max(abs(kwt.field_velocity_bar))),
            potential=dict(Vdev_V=potential.Vdev_V, port_power_W=potential.port_power_W,
                normal_joule_W=potential.normal_joule_W, superconducting_power_W=potential.super_power_W),
            interface_reactions_bar=reactions, energy_rate=asdict(rate), CM9=asdict(balance),
            circuit=dict(parameters=asdict(circuit), state=circuit_state, state_derivative_SI=circuit.rhs(circuit_state, potential.Vdev_V),
                         outputs=asdict(circuit.outputs(circuit_state, potential.Vdev_V)), CM8=asdict(circuit_balance)),
            heat=dict(condensate_rate_bar=heating.condensate_rate_bar, normal_rate_bar=heating.normal_rate_bar,
                      deposited_rate_bar=heating.deposited_rate_bar, moment_residual_bar=heating.moment_residual_bar),
            artifacts_sha256={name: sha(out/name) for name in ("preparation.npz", "initial.npz", "arrays.npz", "derivative_checks.json")},
            stage3_closed=False, temporal_admission=False, kinetic_interface_admission=False, production=False)
        if seeded is not None:
            record["artifacts_sha256"]["direct_reference_probes.json"] = sha(out/"direct_reference_probes.json")
        save(out/"result.json", record)
        progress.finish(success=bool(passed))
        return record
    except BaseException:
        progress.finish(success=False)
        raise


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--registration", default="docs/implementation/stage3/coupled_20260923/mixed_registration.json")
    parser.add_argument("--output-root")
    parser.add_argument("--patch-dir", help="Disabled: no interpolated energy candidate was admitted")
    parser.add_argument("--seed-dir", help="Directory with compatible PASS_SEEDED_CAUSAL_CATALOG evidence; exact Newton correction only")
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--pilot", action="store_true")
    mode.add_argument("--execute", action="store_true")
    args = parser.parse_args()
    if args.patch_dir:
        raise ValueError("Energy patches were rejected; --patch-dir is disabled. No fallback is allowed.")
    registration = (ROOT/args.registration).resolve()
    reg = json.loads(registration.read_text(encoding="utf-8"))
    validate_registration(reg)
    seed_metadata = verify_seed_directory(args.seed_dir, reg) if args.seed_dir else None
    meshes = [m for m in reg["meshes"] if not args.pilot or m["id"] == reg["pilot_mesh"]]
    tasks = [(m, p) for m in meshes for p in reg["profiles"]]
    hashes = sources(registration)
    if not (args.pilot or args.execute):
        print(json.dumps(dict(status="READY_NOT_EXECUTED", cases=[m["id"]+"_"+p["id"] for m, p in tasks],
            scope=reg["scope"], source_sha256=hashes, criteria=reg["criteria"],
            long_run_policy=reg["execution"], seed_validation=seed_metadata), indent=2, ensure_ascii=False))
        return
    if not args.output_root:
        parser.error("--output-root is required for execution")
    output = (ROOT/args.output_root).resolve()
    if output.exists():
        raise FileExistsError("Output exists; choose a new directory. No overwrite, fallback or retry.")
    output.mkdir(parents=True)
    save(output/"manifest.json", dict(registration=reg, source_sha256=hashes, arguments=vars(args),
        started_at_utc=datetime.now(timezone.utc).isoformat(), python=platform.python_version(), numpy=np.__version__,
        seed_validation=seed_metadata))
    started = time.monotonic()
    records = []
    with (output/"progress.jsonl").open("w", encoding="utf-8", buffering=1) as log:
        def event(row):
            log.write(json.dumps(row, allow_nan=False)+"\n")
            log.flush()
        def weight(mesh):
            return ((mesh["elements_x"]*reg["degree"]+1)*(mesh["elements_y"]*reg["degree"]+1)
                    +(mesh["left_elements"]+mesh["right_elements"])*reg["degree"]+2)
        progress = Progress(len(tasks), "Mixto: identidades instantáneas", callback=event, min_interval=2.,
                            total_weight=sum(weight(m) for m, p in tasks))
        try:
            direct = refined_count_catalog(OccupationEnergyCatalog.load(ROOT/reg["catalogue"]), refinement=reg["refinement"])
            if len(direct.count_nodes) != reg["count_nodes"]:
                raise ValueError("unexpected direct count quadrature")
            seeded = (None if seed_metadata is None else
                      SeededCountCatalog(direct, LocalEnergyPatch.load(seed_metadata["artifact"], direct)))
            for mesh, profile in tasks:
                name = mesh["id"]+"_"+profile["id"]
                progress.start_task(name, weight=weight(mesh))
                target = output/name
                target.mkdir()
                def case_event(row):
                    event(row)
                    progress.checkpoint("lote en ejecución; detalle por cuadratura arriba")
                record = solve_case(direct, reg, mesh, profile, target, case_event,
                                    seeded=seeded, seed_metadata=seed_metadata)
                records.append(record)
                if record["status"] != "PASS_MIXED_INSTANTANEOUS_IDENTITIES":
                    raise RuntimeError("Mixed-case identity gate not met; result preserved and sequence stopped")
                progress.advance()
            summary = dict(status="PILOT_INSTANTANEOUS_PASS_SCOPE_LIMITED" if args.pilot else "MIXED_INSTANTANEOUS_BATCH_PASS_SCOPE_LIMITED",
                runtime_seconds=time.monotonic()-started, cases=records, source_sha256=hashes,
                sensitivity_only=True, mesh_convergence_admission=False, temporal_admission=False,
                stage3_closed=False, production=False, exclusions=reg["excluded"])
            save(output/"summary.json", summary)
            progress.finish(success=True)
        except BaseException as error:
            save(output/"failure.json", dict(exception=type(error).__name__, reason=str(error),
                runtime_seconds=time.monotonic()-started, completed_cases=len(records), traceback=traceback.format_exc()))
            progress.finish(success=False)
            raise


if __name__ == "__main__":
    main()
