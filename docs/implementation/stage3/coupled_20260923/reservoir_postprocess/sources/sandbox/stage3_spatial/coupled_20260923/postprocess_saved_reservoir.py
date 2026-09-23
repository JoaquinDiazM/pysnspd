"""Re-evaluate instantaneous reservoir loads using saved fields and FD spectra.

No new spectral query and no time step. The original free-boundary run remains
unchanged. The potential solve, mobility, heat deposition and ledger are cheap
algebra on archived arrays. E is recovered ONLY from the documented saved FD
preparation and verified against its original heat moment before further use.
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
from types import SimpleNamespace

import numpy as np
from scipy.special import expit

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from progress import Progress
from pysnspd.experimental.energy_catalog import K_B_J_K, E_CHARGE_C, HBAR_J_S
from pysnspd.experimental.cell_closures import CellScales, KWTMobility
from pysnspd.experimental.cell_validation import ElectronicCell
from pysnspd.experimental.mixed_spatial import MixedSpatialFunctional
from pysnspd.experimental.electrical_ports import solve_potential, ThesisCircuitParameters
from pysnspd.experimental.reservoir_spatial_dynamics import fixed_radius_reservoir_kwt
from pysnspd.experimental.spatial_dynamics import (
    normal_heat_distribution, deposit_spatial_heat, spatial_energy_rate, cm9_power_balance, gauge_power_check,
)


def plain(v):
    if isinstance(v, np.ndarray):
        return plain(v.tolist())
    if isinstance(v, np.generic):
        return plain(v.item())
    if isinstance(v, dict):
        return {k: plain(x) for k, x in v.items()}
    if isinstance(v, (list, tuple)):
        return [plain(x) for x in v]
    return v


def save(path, value):
    Path(path).write_text(json.dumps(plain(value), ensure_ascii=False, indent=2, allow_nan=False)+"\n", encoding="utf-8")


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read(path):
    return json.loads(Path(path).read_text(encoding="utf-8"))


def forbidden_spectral_query(*args, **kwargs):
    raise RuntimeError("This postprocessor may not evaluate or reconstruct a new causal spectrum")


class RecoveredThermalCell:
    """Saved-state energy values, exposing unchanged D.18 arithmetic only."""
    validate = ElectronicCell.validate
    excitation_energy = ElectronicCell.excitation_energy
    equivalent_temperature = ElectronicCell.equivalent_temperature
    heating = ElectronicCell.heating

    def __init__(self, catalog, amplitude, gamma, energies, weights):
        self.catalog, self.amplitude, self.gamma = catalog, float(amplitude), float(gamma)
        self.energies, self.weights = np.array(energies, copy=True), np.array(weights, copy=True)
        self.energies.setflags(write=False)
        self.weights.setflags(write=False)


def verify_input(input_root, registration):
    manifest, summary = read(input_root/"manifest.json"), read(input_root/"summary.json")
    if (summary.get("status") not in ("PILOT_INSTANTANEOUS_PASS_SCOPE_LIMITED", "MIXED_INSTANTANEOUS_BATCH_PASS_SCOPE_LIMITED")
            or summary.get("stage3_closed") is not False or summary.get("production") is not False):
        raise ValueError("input batch is not a completed scoped mixed diagnostic")
    sources = manifest["source_sha256"]
    if summary["source_sha256"] != sources:
        raise ValueError("input source manifests disagree")
    for path, digest in sources.items():
        if sha(ROOT/path) != digest:
            raise ValueError("archived source changed: "+path)
    runner_path = "sandbox/stage3_spatial/coupled_20260923/run_mixed_check.py"
    reg_path = "docs/implementation/stage3/coupled_20260923/mixed_registration.json"
    if (sources.get(runner_path) != registration["input_mixed_runner_sha256"]
            or sources.get(reg_path) != registration["input_mixed_registration_sha256"]
            or read(ROOT/reg_path) != manifest["registration"]):
        raise ValueError("unregistered mixed generator or preparation")
    oldreg = manifest["registration"]
    expected = {(m["id"], p["id"]) for m in oldreg["meshes"] for p in oldreg["profiles"]
                if summary["status"] != "PILOT_INSTANTANEOUS_PASS_SCOPE_LIMITED" or m["id"] == oldreg["pilot_mesh"]}
    actual = [(r["mesh"]["id"], r["profile"]["id"]) for r in summary["cases"]]
    if set(actual) != expected or len(actual) != len(expected):
        raise ValueError("input case coverage incomplete or duplicated")
    records, inventory = [], {"manifest.json": sha(input_root/"manifest.json"), "summary.json": sha(input_root/"summary.json")}
    for item in summary["cases"]:
        name = item["mesh"]["id"]+"_"+item["profile"]["id"]
        path = input_root/name
        record = read(path/"result.json")
        if record != item or record["status"] != "PASS_MIXED_INSTANTANEOUS_IDENTITIES":
            raise ValueError("input case has no matching PASS record: "+name)
        for artifact, expected_sha in record["artifacts_sha256"].items():
            if sha(path/artifact) != expected_sha:
                raise ValueError("input artifact changed: "+name+"/"+artifact)
            inventory[name+"/"+artifact] = expected_sha
        inventory[name+"/result.json"] = sha(path/"result.json")
        records.append(record)
    return manifest, records, inventory


def geometry_catalog(path, nodes, weights):
    """Read only unit metadata; all spectral APIs deliberately fail."""
    with np.load(path, allow_pickle=False) as raw:
        metadata = json.loads(str(raw["metadata_json"]))
        vacuum = SimpleNamespace(delta0_J=float(raw["delta0_J"]), N0_per_J_m3=float(raw["N0_per_J_m3"]),
            D_m2_s=float(raw["D_m2_s"]), metadata=metadata, evaluate=forbidden_spectral_query)
    return SimpleNamespace(vacuum=vacuum, count_nodes=nodes, count_weights=weights, energy_kernel=forbidden_spectral_query)


def allocation(model):
    value = np.zeros((model.quadrature_size, len(model.graph_edges)))
    for edge, ends in enumerate(model.graph_edges):
        for node in ends:
            rows = np.flatnonzero(model.quadrature_to_dof == node)
            value[rows, edge] += .5*model.quadrature_mass_bar[rows]/model.mass_bar[node]
    return value


def postprocess_case(input_root, original, oldreg, reg, output):
    started = time.monotonic()
    name = original["mesh"]["id"]+"_"+original["profile"]["id"]
    src = input_root/name
    with np.load(src/"arrays.npz", allow_pickle=False) as handle:
        arrays = {k: handle[k] for k in handle.files}
    with np.load(src/"initial.npz", allow_pickle=False) as handle:
        initial = {k: handle[k] for k in handle.files}
    with np.load(src/"preparation.npz", allow_pickle=False) as handle:
        preparation = {k: handle[k] for k in handle.files}
    z, p = arrays["delta_bar"], arrays["p_quadrature"]
    if not np.array_equal(z, initial["delta_bar"]) or not np.array_equal(p, initial["p_quadrature"]):
        raise ValueError("saved field/population disagree with the initial snapshot")
    theta = oldreg["bath_theta"]
    if np.any(~np.isfinite(p)) or np.any((p <= 0)|(p >= .5)):
        raise ValueError("FD energies cannot be recovered as finite positive values from this snapshot")
    energies = theta*(np.log1p(-p)-np.log(p))
    if np.any(~np.isfinite(energies)) or np.any(energies <= 0) or np.any(np.diff(energies, axis=1) <= 0):
        raise ValueError("recovered saved energies are not positive and strictly ordered")
    fd_error = float(np.max(abs(expit(-energies/theta)-p)/p))
    mass, mapping = preparation["quadrature_mass_bar"], preparation["quadrature_to_dof"]
    weights = initial["count_weights"]
    old_heat_moments = 4*np.sum(weights*energies*arrays["electron_heating_rhs"], axis=1)
    recovered_old_heat = float(np.dot(mass, old_heat_moments))
    original_heat = original["heat"]["deposited_rate_bar"]
    moment_error = abs(recovered_old_heat-original_heat)/max(abs(original_heat), 1e-30)
    temperature_error = float(np.max(abs(arrays["temperature_bar"]-theta)))
    preconditions = dict(original_heat_moment_relative=moment_error, FD_recovery_relative=fd_error,
                         temperature_absolute=temperature_error)
    save(output/"recovery_preconditions.json", preconditions)
    if any(preconditions[k] > reg["criteria"][k] for k in preconditions):
        raise ValueError("saved FD decoding did not reproduce the archived state/energy moment")
    catalog = geometry_catalog(ROOT/oldreg["catalogue"], initial["count_nodes"], weights)
    mesh = original["mesh"]
    model = MixedSpatialFunctional(catalog, **oldreg["geometry"], degree=oldreg["degree"],
        **{k: mesh[k] for k in ("elements_x", "elements_y", "left_elements", "right_elements")})
    for label, stored, measured in (("identification", mapping, model.quadrature_to_dof),
                                     ("mass", mass, model.quadrature_mass_bar),
                                     ("edges", preparation["graph_edges"], model.graph_edges)):
        if not np.array_equal(stored, measured):
            raise ValueError("reassembled geometry changed: "+label)
    if model.energy_scale_J != original["energy_scale_J"]:
        raise ValueError("energy units changed")
    reference_current = reg["reservoir"]["reference_current_A"]
    if reference_current != oldreg["circuit"]["reference_current_A"]:
        raise ValueError("reservoir current is not the independent original reference")
    terminal = np.array(model.terminal_dofs)
    injection = np.zeros(model.cells)
    injection[terminal] = [reference_current, -reference_current]
    conductance = model.conductance(oldreg["normal_conductivity_S_m"])
    potential = solve_potential(model.graph_edges, conductance, arrays["current_super_A"], injection,
                                left_node=int(terminal[0]), right_node=int(terminal[1]))
    scales = CellScales(catalog.vacuum.delta0_J, catalog.vacuum.N0_per_J_m3, model.Tc_K,
                       theta*catalog.vacuum.delta0_J/K_B_J_K, oldreg["time_reference_ps"])
    response = fixed_radius_reservoir_kwt(z, arrays["cartesian_gradient_bar"], mass,
        KWTMobility(scales), arrays["temperature_bar"], potential.phi_V,
        quadrature_to_node=mapping, terminal_nodes=terminal, terminal_injection_A=injection[terminal],
        fixed_radii_bar=np.full(2, reg["reservoir"]["fixed_terminal_amplitude_bar"]), energy_scale_J=model.energy_scale_J)
    cells = [RecoveredThermalCell(catalog, amplitude, gamma, values, weights) for amplitude, gamma, values in
             zip(arrays["amplitude_quadrature_bar"], arrays["gamma_quadrature_bar"], energies)]
    normal = normal_heat_distribution(conductance, potential.delta_phi_V, allocation(model))
    heating = deposit_spatial_heat(cells, p, response.kwt, normal.quadrature_power_W,
                                  energy_scale_J=model.energy_scale_J, scales=scales)
    rate = spatial_energy_rate(arrays["cartesian_gradient_bar"], response.kwt.field_velocity_bar,
                              cells, heating.electron_rhs, mass)
    params = dict(original["circuit"]["parameters"])
    params["V_bias_V"] = params["R_bias_ohm"]*reference_current
    circuit = ThesisCircuitParameters(**params)
    state = np.array([reference_current, reference_current, 0.])
    circuit_balance = circuit.power_balance(state, potential.Vdev_V)
    power_scale_factor = model.energy_scale_J/scales.t_ref_s
    reservoir_power = response.reservoir_work_rate_bar*power_scale_factor
    balance = cm9_power_balance(rate.total_rate_bar, energy_scale_J=model.energy_scale_J, time_scale_s=scales.t_ref_s,
                                circuit_balance=circuit_balance, reservoir_power_into_W=reservoir_power)
    gauge = gauge_power_check(z, arrays["cartesian_gradient_bar"], model.graph_edges,
        arrays["current_super_A"]/(2*E_CHARGE_C/HBAR_J_S*model.energy_scale_J),
        potential.phi_V, energy_scale_J=model.energy_scale_J, time_scale_s=scales.t_ref_s)
    power_norm = max(abs(balance.domain_energy_rate_W), abs(circuit_balance.device_power_W),
                     abs(reservoir_power), potential.normal_joule_W, 1e-30)
    work_norm = max(abs(response.phase_reservoir_work_rate_bar), response.kwt.condensate_heat_rate_bar, 1e-30)
    heat_norm = max(abs(heating.deposited_rate_bar), abs(heating.condensate_rate_bar), abs(heating.normal_rate_bar), 1e-30)
    kwt_norm = max(abs(response.kwt.field_energy_rate_bar), response.kwt.condensate_heat_rate_bar,
                   abs(response.kwt.boundary_work_rate_bar), abs(response.kwt.electrical_phase_work_rate_bar), 1e-30)
    metrics = dict(**preconditions, terminal_radius_rate_absolute=float(max(abs(response.terminal_radius_rates_bar))),
        radial_reaction_work_relative=abs(response.radial_reaction_work_rate_bar)/work_norm,
        instantaneous_power_relative=max(abs(balance.residual_W), abs(balance.device_port_residual_W),
            abs(circuit_balance.residual_W), abs(gauge.power_residual_W), abs(potential.power_residual_W))/power_norm,
        heat_moment_relative=abs(heating.moment_residual_bar)/heat_norm,
        current_continuity_relative=float(max(abs(potential.current_residual_A))/reference_current),
        KWT_identity_relative=abs(response.kwt.identity_residual_bar)/kwt_norm)
    passed = all(np.isfinite(metrics[k]) and metrics[k] <= limit for k, limit in reg["criteria"].items())
    adjacent = []
    for node in terminal:
        edge = np.flatnonzero(np.any(model.graph_edges == node, axis=1))
        if len(edge) != 1:
            raise ValueError("the selected remote lead terminal must have one adjacent internal graph edge")
        adjacent.append(int(edge[0]))
    complex_gradient = arrays["cartesian_gradient_bar"][:, 0]+1j*arrays["cartesian_gradient_bar"][:, 1]
    free_node_heat = np.bincount(mapping, weights=mass*arrays["QDelta_density_bar"], minlength=model.cells)
    node_heat = np.bincount(mapping, weights=mass*response.kwt.heat_density_bar, minlength=model.cells)
    np.savez_compressed(output/"arrays.npz", delta_bar=z, recovered_saved_energies_bar=energies,
        p_quadrature=p, phi_V=potential.phi_V, current_normal_A=potential.I_normal_A,
        current_total_A=potential.I_total_A, material_velocity_bar=response.kwt.material_velocity_bar,
        field_velocity_bar=response.kwt.field_velocity_bar, QDelta_density_bar=response.kwt.heat_density_bar,
        electron_heating_rhs=np.array(heating.electron_rhs), boundary_load_cartesian_bar=response.load.load_cartesian_bar,
        radial_reaction_cartesian_bar=response.load.radial_load_cartesian_bar,
        normal_heat_power_W=normal.quadrature_power_W, quadrature_mass_bar=mass,
        quadrature_to_dof=mapping, dof_coordinates_m=preparation["dof_coordinates_m"])
    result = dict(case_id=name, status="PASS_SAVED_RESERVOIR_INSTANTANEOUS_IDENTITIES" if passed else "SAVED_RESERVOIR_NOT_ACCEPTED",
        runtime_seconds=time.monotonic()-started, reference_current_A=reference_current, mesh=mesh,
        energy_scale_J=model.energy_scale_J, time_scale_s=scales.t_ref_s, metrics=metrics,
        free=dict(condensate_rate_bar=original["heat"]["condensate_rate_bar"],
                  maximum_material_speed_bar=original["maximum_material_speed_bar"],
                  terminal_condensate_rate_bar=float(sum(free_node_heat[terminal]))),
        reservoir=dict(condensate_rate_bar=heating.condensate_rate_bar,
            maximum_material_speed_bar=float(max(abs(response.kwt.material_velocity_bar))),
            terminal_condensate_rate_bar=float(sum(node_heat[terminal])),
            work_rate_bar=response.reservoir_work_rate_bar, work_power_W=reservoir_power,
            phase_work_rate_bar=response.phase_reservoir_work_rate_bar,
            radial_reaction_work_rate_bar=response.radial_reaction_work_rate_bar,
            decomposition_residual_bar=response.decomposition_residual_bar),
        terminal=dict(nodes=terminal, injection_A=injection[terminal], phase_load_bar=response.load.phase_load_bar,
            radial_reaction_bar=response.load.radial_reaction_bar, fixed_radius_bar=response.load.fixed_radii_bar,
            radius_rate_bar=response.terminal_radius_rates_bar,
            phase_balance_residual_bar=np.imag(np.conj(z[terminal])*complex_gradient[terminal])-response.load.phase_load_bar,
            adjacent_internal_edges=adjacent, adjacent_internal_normal_current_A=potential.I_normal_A[adjacent],
            adjacent_internal_normal_current_over_Iref=potential.I_normal_A[adjacent]/reference_current,
            normal_current_scope="Interior graph face next to remote terminal; not an admitted continuous external-boundary normal flux.",
            external_normal_injection_algebraic_A=injection[terminal]+response.load.phase_load_bar*response.load.current_scale_A),
        potential=dict(Vdev_V=potential.Vdev_V, port_power_W=potential.port_power_W,
            normal_joule_W=potential.normal_joule_W, super_power_W=potential.super_power_W),
        circuit=dict(parameters=params, state=state, state_derivative_SI=circuit.rhs(state, potential.Vdev_V),
                     CM8=asdict(circuit_balance)), CM9=asdict(balance), energy_rate=asdict(rate),
        heat=dict(deposited_rate_bar=heating.deposited_rate_bar, condensate_rate_bar=heating.condensate_rate_bar,
                  normal_rate_bar=heating.normal_rate_bar, residual_bar=heating.moment_residual_bar),
        recovery=dict(archived_deposited_rate_bar=original_heat, reconstructed_original_deposited_rate_bar=recovered_old_heat,
                      spectrum_scope="Decoded ONLY from saved thermal populations and their declared temperature; no new causal spectrum."),
        input_record_sha256=sha(src/"result.json"), artifacts_sha256={n: sha(output/n) for n in ("arrays.npz", "recovery_preconditions.json")},
        new_spectral_queries=0, time_steps=0, stage3_closed=False, temporal_admission=False,
        full_D27_admission=False, normal_boundary_trace_admission=False, production=False)
    save(output/"result.json", result)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-root", default="docs/implementation/stage3/coupled_20260923/raw/mixed_pilot")
    parser.add_argument("--output-root", required=True)
    parser.add_argument("--registration", default="docs/implementation/stage3/coupled_20260923/reservoir_postprocess_registration.json")
    args = parser.parse_args()
    registration = (ROOT/args.registration).resolve()
    reg = read(registration)
    input_root, output = (ROOT/args.input_root).resolve(), (ROOT/args.output_root).resolve()
    manifest, records, inventory = verify_input(input_root, reg)
    if output.exists():
        raise FileExistsError("Output directory exists; no overwrite or retry")
    output.mkdir(parents=True)
    names = ("reservoir_spatial_dynamics", "spatial_dynamics", "cell_closures", "cell_validation", "electrical_ports", "mixed_spatial", "spatial_open", "spatial_functional", "energy_catalog")
    paths = [ROOT/("pysnspd/experimental/"+n+".py") for n in names]
    paths += [Path(__file__).resolve(), registration, ROOT/"sandbox/stage3_spatial/progress.py"]
    save(output/"manifest.json", dict(registration=reg, source_sha256={str(p.relative_to(ROOT)).replace("\\", "/"): sha(p) for p in paths},
        input_root=str(input_root), input_artifacts_sha256=inventory,
        started_at_utc=datetime.now(timezone.utc).isoformat(), python=platform.python_version(), numpy=np.__version__,
        new_spectral_queries=0, time_steps=0))
    started = time.monotonic()
    results = []
    progress = Progress(len(records), "Reservorios: postproceso instantáneo", min_interval=1.)
    try:
        for record in records:
            name = record["mesh"]["id"]+"_"+record["profile"]["id"]
            progress.start_task(name)
            target = output/name
            target.mkdir()
            result = postprocess_case(input_root, record, manifest["registration"], reg, target)
            results.append(result)
            if not result["status"].startswith("PASS_"):
                raise RuntimeError("reservoir postprocess criteria failed; evidence preserved")
            progress.advance()
        summary = dict(status="PASS_SAVED_RESERVOIR_SNAPSHOTS_SCOPE_LIMITED", runtime_seconds=time.monotonic()-started,
            cases=results, manifest_sha256=sha(output/"manifest.json"), new_spectral_queries=0, time_steps=0,
            sensitivity_only=True, stage3_closed=False, temporal_admission=False, full_D27_admission=False,
            normal_boundary_trace_admission=False, production=False, exclusions=reg["excluded"])
        save(output/"summary.json", summary)
        progress.finish(success=True)
        print(json.dumps(dict(status=summary["status"], runtime_seconds=summary["runtime_seconds"], output_root=str(output))), flush=True)
    except BaseException as error:
        save(output/"failure.json", dict(exception=type(error).__name__, reason=str(error), traceback=traceback.format_exc(),
                                        completed_cases=len(results), runtime_seconds=time.monotonic()-started))
        progress.finish(success=False)
        raise


if __name__ == "__main__":
    main()
