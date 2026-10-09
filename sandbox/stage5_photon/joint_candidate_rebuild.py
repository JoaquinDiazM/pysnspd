"""Finite ROS1 endpoint reconstruction, without time/material admission.

Called only by the explicit user-run campaign: independent retarded and
Matsubara queries can use a bounded fork pool. The evolved Wilson gap, links,
and full hL/hT fields are never replaced by a fitted temperature/profile.
Expected physical/root failures return diagnostics, so dt refinement proceeds.
"""
from __future__ import annotations
from concurrent.futures import ProcessPoolExecutor, as_completed
import hashlib
import json
import multiprocessing as mp
import os
from pathlib import Path
import time
import numpy as np
from scipy.constants import Boltzmann as kB
from pysnspd.experimental import adiabatic_matrix_usadel as ad
from pysnspd.experimental import retarded_spatial_usadel as ret
from pysnspd.experimental import thermal_spatial_usadel as thermal
from pysnspd.experimental import causal_kinetic_regulator as reg
from pysnspd.experimental.bulk_current_reference import (
    solve_bulk_reference, solve_bulk_current, homogeneous_spectrum)
from pysnspd.experimental.thermal_internal_energy import evaluate_internal_energy
from pysnspd.experimental.adiabatic_energy_observer import (
    AdiabaticElectronicEnergy, phonon_energy_J)
from sandbox.stage4_core.parallel_runtime import linux_resources, resource_budget
from sandbox.stage5_rhs.run_matrix_rhs_audit import outer_angle

_CONTEXT = None


def fixed_energy_distribution_derivative(energy_J, energy_unit_J, h):
    """Three-point nonuniform finite difference dh/d(E/Eu); no FD/Te fit."""
    energy = np.asarray(energy_J, float)/float(energy_unit_J)
    values = np.asarray(h, float)
    if (energy.ndim != 1 or len(energy) < 2 or np.any(np.diff(energy) <= 0)
            or values.shape[0] != len(energy) or np.any(~np.isfinite(values))):
        raise ValueError("Ordered fixed-energy nodes and finite full distributions required")
    return np.gradient(values, energy, axis=0, edge_order=2 if len(energy) > 2 else 1)


def thermal_cut_fields(epsilon_bar, gap_boundary, q_bar, temperature_ratio):
    """Actual cut gaps with current-continuation q; exact fixed-d,q T tangent."""
    d = np.asarray(gap_boundary, complex)
    epsilon = np.asarray(epsilon_bar, float)
    if np.any(abs(d) <= 0):
        raise ValueError("The superconducting continuation has a zero cut gap")
    u = np.empty((len(epsilon), len(d)), complex)
    du = np.empty_like(u)
    for j, gap in enumerate(d):
        radial, _, g = homogeneous_spectrum(abs(gap), epsilon, q_bar*q_bar)
        u[:, j] = gap/abs(gap)*radial
        du[:, j] = -(epsilon/temperature_ratio)*u[:, j]/(epsilon+q_bar*q_bar*g**3)
    return u, du


def candidate_state_screen(snapshot, candidate, *, diagnostic_excess_limit=1e-7):
    """Preserve exact support flag; tiny violations allow an observer ONLY."""
    n = snapshot.graph.n_nodes
    shape = (len(snapshot.slices), n, 2)
    common = dict(ok=False, ok_for_diagnostic_rebuild=False, support_ok=False,
        maximum_excess=None, accepted_physical_time_steps=0,
        diagnostic_excess_limit=float(diagnostic_excess_limit), no_clipping=True)
    required = {"gap_bar_new": (n,), "alpha_new": (len(snapshot.graph.edges),),
        "h_new": shape, "circuit_state_new_SI": (3,)}
    for key, expected in required.items():
        if key not in candidate or candidate[key].shape != expected:
            return dict(common, reason="candidate_shape", field=key)
        if np.any(~np.isfinite(candidate[key])):
            return dict(common, reason="candidate_nonfinite", field=key)
        if key != "gap_bar_new" and np.iscomplexobj(candidate[key]):
            return dict(common, reason="candidate_complex_real_field", field=key)
    h = candidate["h_new"]
    excess = float(np.max(np.maximum(abs(h[..., 0]+h[..., 1]),
        abs(h[..., 0]-h[..., 1]))-1., initial=0.))
    common.update(maximum_excess=excess,
        support_ok=bool(excess <= 32*np.finfo(float).eps))
    current = candidate["circuit_state_new_SI"]
    if excess > diagnostic_excess_limit:
        return dict(common, reason="Pauli_support")
    if current[0] <= 0 or current[1] <= 0:
        return dict(common, reason="global_current_direction")
    if np.any(abs(candidate["gap_bar_new"][snapshot.fixed]) <= 0):
        return dict(common, reason="zero_superconducting_cut_gap")
    return dict(common, ok=True, ok_for_diagnostic_rebuild=True,
        decision="diagnostic reconstruction only; exact support is independent and no step is accepted")


def support_excess_energy_J(energy_J, weights_J, volume_m3, N0_per_J_m3, rho, h):
    """Energy-weighted support excess on the NEW DOS; diagnostic, no bound."""
    values = np.asarray(h, float)
    excess_plus = np.maximum(abs(values[..., 0]+values[..., 1])-1., 0.)
    excess_minus = np.maximum(abs(values[..., 0]-values[..., 1])-1., 0.)
    return float(N0_per_J_m3*np.einsum("e,n,en,en->",
        np.asarray(energy_J)*np.asarray(weights_J), np.asarray(volume_m3),
        np.asarray(rho), excess_plus+excess_minus))


def _write_json(path, value):
    with Path(path).open("x") as stream:
        json.dump(value, stream, indent=2, allow_nan=False)
        stream.write("\n")


def _save(path, **fields):
    path = Path(path)
    temporary = path.with_suffix(path.suffix+".partial")
    with temporary.open("wb") as stream:
        np.savez_compressed(stream, **fields)
    temporary.replace(path)


def _sha(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024*1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _worker_init(cpus, counter):
    with counter.get_lock():
        index = counter.value
        counter.value += 1
    os.sched_setaffinity(0, {cpus[index % len(cpus)]})


def _query(task):
    kind, index = task
    c = _CONTEXT
    try:
        if kind == "retarded":
            sl = c["snapshot"].slices[index]
            cache = c["output"]/"retarded"/("mode_%05d.npz" % index)
            if cache.exists():
                with np.load(cache, allow_pickle=False) as z:
                    R = z["R0"].copy()
                    receipt = json.loads(str(z["receipt"]))
            else:
                fixed = c["snapshot"].fixed
                zeta = sl.eta-1j*sl.energy
                theta = np.array([outer_angle(zeta, abs(d), c["bulk"].q_bar, seed)
                    for d, seed in zip(c["gap"][fixed], sl.theta)])
                phase = c["gap"][fixed]/abs(c["gap"][fixed])
                initial = sl.Rbase.value
                root = ret.solve(c["graph"], c["gap"], zeta, c["alpha"],
                    initial_a=initial[:, 0, 1]/(1+initial[:, 0, 0]),
                    initial_b=initial[:, 1, 0]/(1+initial[:, 0, 0]),
                    fixed_nodes=fixed, fixed_a=np.tan(theta/2)*phase,
                    fixed_b=np.tan(theta/2)*phase.conj(), tolerance=1e-9,
                    max_iterations=100, causal_tolerance=1e-8)
                R = ret.nambu_matrix(root.f, root.f_tilde, root.g)
                receipt = dict(index=index, energy_J=sl.energy*sl.Eu,
                    chart_scaled_residual=root.residual,
                    normalization_residual=root.normalization_residual,
                    minimum_DOS=root.minimum_DOS, iterations=root.iterations)
                _save(cache, R0=R, receipt=np.asarray(json.dumps(receipt)))
            H = ad.distribution(c["h"][index, :, 0], c["h"][index, :, 1])
            Hth = ad.distribution(np.full(c["graph"].n_nodes, c["h0"][index]),
                np.zeros(c["graph"].n_nodes))
            B = ad.generator(c["gap"], sl.energy)
            residual, flux, K, _ = reg.kinetic_residual(c["graph"], ad.Jet.of(R), H, B,
                alpha=c["alpha"], numerical_eta=sl.eta, kappa_s=sl.kappa)
            _, fluxth, Kth = ad.kinetic_residual(c["graph"], ad.Jet.of(R), Hth, B,
                alpha=c["alpha"], kappa=0.)
            return kind, index, dict(R0=R, K0=K.value, flux=flux.value,
                Kth=Kth.value, fluxth=fluxth.value, residual=residual.value,
                receipt=receipt)
        cache = c["output"]/"matsubara"/("mode_%05d.npz" % index)
        if cache.exists():
            with np.load(cache, allow_pickle=False) as z:
                solution = thermal.SpectralSolution(float(z["epsilon"]), z["u"].copy(),
                    z["f"].copy(), z["g"].copy(), float(z["energy"]), float(z["residual"]),
                    float(z["gradient_residual"]), int(z["iterations"]),
                    tuple(z["line_search_steps"].tolist()), z["fixed_nodes"].copy(),
                    str(z["source_signature"]))
        else:
            with np.load(c["observer"]/"spectral_receipts"/("mode_%05d.npz" % index),
                    allow_pickle=False) as z:
                initial_u = z["u"].copy()
                if float(z["epsilon_bar"]) != c["epsilon"][index]:
                    raise ValueError("The cached Matsubara temperature/frequency changed")
            solution = thermal.solve_frequency(c["graph"], c["gap"], c["epsilon"][index],
                c["alpha"], fixed_nodes=c["snapshot"].fixed,
                fixed_u=c["thermal_cut"][index], initial_u=initial_u,
                tol=1e-9, max_iterations=100)
            _save(cache, epsilon=solution.epsilon, u=solution.u, f=solution.f, g=solution.g,
                energy=solution.energy, residual=solution.residual,
                gradient_residual=solution.gradient_residual, iterations=solution.iterations,
                line_search_steps=np.asarray(solution.line_search_steps, int),
                fixed_nodes=solution.fixed_nodes, source_signature=np.asarray(solution.source_signature))
        return kind, index, dict(solution=solution)
    except Exception as error:
        return kind, index, dict(error=repr(error), error_type=type(error).__name__)


def _parallel_queries(tasks, budget):
    """Stop scheduling after a physical query failure; preserve completed roots."""
    if budget["workers"] == 1:
        for task in tasks:
            value = _query(task)
            yield value
            if "error" in value[2]:
                return
        return
    counter = mp.get_context("fork").Value("i", 0)
    usable = set(budget["selected_affinity_cpus"])
    groups = linux_resources()["physical_core_groups"]
    cpus = [group[0] for group in groups if group[0] in usable]
    cpus += [cpu for group in groups for cpu in group[1:] if cpu in usable]
    cpus = [cpu for cpu in cpus if cpu != budget["coordinator_cpu"]][:budget["workers"]]
    original = os.sched_getaffinity(0)
    os.sched_setaffinity(0, {budget["coordinator_cpu"]})
    try:
        with ProcessPoolExecutor(budget["workers"], mp_context=mp.get_context("fork"),
                initializer=_worker_init, initargs=(cpus, counter)) as pool:
            futures = [pool.submit(_query, task) for task in tasks]
            for future in as_completed(futures):
                value = future.result()
                if "error" in value[2]:
                    # Cancel before yielding: the caller can immediately
                    # return on this diagnostic and close our generator.
                    for other in futures:
                        other.cancel()
                    yield value
                    return
                yield value
    finally:
        os.sched_setaffinity(0, original)


def rebuild_candidate(snapshot, candidate_npz_path, output_dir, *, workers=1):
    """Re-solve finite R0 and F_N/U_N; always return explicit admission limits.

    The bulk branch is anchored to its ORIGINAL ideal-bulk current plus the
    change in circuit Is. Graph port discretization is not silently rebased.
    The actual evolved cut gaps remain untouched; their nonlinear mismatch
    with that bulk branch is a measured result.
    """
    global _CONTEXT
    output = Path(output_dir).resolve()
    if not output.is_relative_to(Path("/home/jdiaz/scratch")):
        raise ValueError("Finite-state numerical data must remain in scratch")
    with np.load(candidate_npz_path, allow_pickle=False) as saved:
        candidate = {key: saved[key].copy() for key in saved.files}
    identity = dict(candidate_sha256=_sha(candidate_npz_path),
        graph_signature=snapshot.graph.signature,
        energy_J=snapshot.rows.moments.energy.tolist(),
        producer_sha256={str(Path(module.__file__).resolve()): _sha(module.__file__)
            for module in (ret, thermal, ad, reg)},
        rebuild_producer_sha256=_sha(__file__))
    if output.exists():
        if json.loads((output/"identity.json").read_text()) != identity:
            raise ValueError("A preserved finite-state rebuild has a different candidate")
        if (output/"summary.json").exists():
            return json.loads((output/"summary.json").read_text())
    else:
        output.mkdir(parents=True)
        _write_json(output/"identity.json", identity)
    result = dict(finite_state_rebuilt=False, source_applied_again=False,
        energy_admission=False, time_admission=False, stage51_complete=False,
        charge_gauge_energy_completion_missing=True,
        new_first_gradient_fields_missing=True, general_collision_provider_available=False,
        finite_state_collision_source_computed=False, backwards_Euler_converged=False)
    started = time.monotonic()
    try:
        screen = candidate_state_screen(snapshot, candidate)
        result["screen"] = screen
        if not screen["ok"]:
            result["status"] = "FINITE_STATE_REJECTED_BEFORE_SPECTRAL_REBUILD"
            _write_json(output/"summary.json", result)
            return result
        data = snapshot.initial_data
        w = data["inputs"]
        plan = w["summary"]["plan"]
        t = plan["T_K"]/plan["Tc_K"]
        count = int(plan["matsubara_count"])
        material = dict(width_m=plan["width_m"], Tc_K=plan["Tc_K"],
            diffusion_m2_s=plan["diffusion_m2_s"],
            sheet_resistance_ohm=plan["sheet_resistance_ohm"])
        oldbulk = solve_bulk_reference(t, count, snapshot.bulk.q_bar)
        oldbulkI = oldbulk.total_current_A(**material)
        deltaI = float(candidate["circuit_state_new_SI"][1]-snapshot.y[1])
        target = oldbulkI+deltaI
        if target <= 0:
            raise ValueError("The anchored exterior bulk current would reverse")
        bulk = solve_bulk_current(t, count, target, **material)
        gap = np.asarray(candidate["gap_bar_new"], complex)
        alpha = np.asarray(candidate["alpha_new"], float)
        h = np.asarray(candidate["h_new"], float)
        hE = fixed_energy_distribution_derivative(snapshot.rows.moments.energy, snapshot.Eu, h)
        fixed = snapshot.fixed
        epsilon = 2*np.pi*t*(np.arange(count)+.5)
        cut, cutT = thermal_cut_fields(epsilon, gap[fixed], bulk.q_bar, t)
        result["bulk_continuation"] = dict(q_new_bar=bulk.q_bar, gap_new_bar=bulk.gap_bar,
            initial_bulk_current_A=oldbulkI, initial_graph_port_current_A=float(snapshot.y[1]),
            initial_bulk_graph_offset_A=oldbulkI-float(snapshot.y[1]),
            imposed_bulk_current_A=target, circuit_current_change_A=deltaI,
            maximum_actual_cut_amplitude_mismatch_bar=float(np.max(abs(abs(gap[fixed])-bulk.gap_bar))),
            maximum_cut_mismatch_from_linear_law_bar=float(np.max(abs(abs(gap[fixed])-
                (snapshot.bulk.gap_bar+snapshot.bulk.gap_bar_derivative_per_A*deltaI)))),
            evolved_gap_not_replaced=True, branch="same ascending finite Matsubara branch anchored to admitted DC")
        budget = resource_budget(linux_resources(), .9, .9, workers)
        rss = 0
        for line in Path("/proc/self/status").read_text().splitlines():
            if line.startswith("VmRSS:"):
                rss = int(line.split()[1])*1024
        reservation = (budget["workers"]+1)*rss+1024**3
        if reservation > budget["memory_limit_bytes"]:
            raise MemoryError("Snapshot/fork endpoint reservation exceeds90% available RAM; reduce workers")
        result["resource_budget"] = budget
        result["scratch_memory_reservation_bytes"] = reservation
        for directory in ("retarded", "matsubara"):
            (output/directory).mkdir(exist_ok=True)
        observer = Path(data["identity"]["plan"]["thermal_observer"])
        h0 = np.tanh(snapshot.rows.moments.energy/(2*kB*plan["T_K"]))
        _CONTEXT = dict(snapshot=snapshot, output=output, graph=snapshot.graph,
            gap=gap, alpha=alpha, h=h, h0=h0, bulk=bulk, epsilon=epsilon,
            thermal_cut=cut, observer=observer)
        tasks = [("retarded", i) for i in range(len(snapshot.slices))]
        tasks += [("matsubara", i) for i in range(count)]
        roots = [None]*len(snapshot.slices)
        modes = [None]*count
        completed = 0
        for kind, index, value in _parallel_queries(tasks, budget):
            completed += 1
            if "error" in value:
                result.update(status="FINITE_STATE_SPECTRAL_REBUILD_REJECTED", failed_query=dict(
                    kind=kind, index=index, **value), completed_queries=completed-1)
                _write_json(output/"summary.json", result)
                return result
            (roots if kind == "retarded" else modes)[index] = value
            if completed == 1 or completed % 16 == 0 or completed == len(tasks):
                print(json.dumps(dict(status="FINITE_STATE_REBUILD_PROGRESS", completed=completed,
                    total=len(tasks), workers=budget["workers"], elapsed_s=time.monotonic()-started)), flush=True)
        solutions = tuple(item["solution"] for item in modes)
        reference = thermal.evaluate_thermal(snapshot.graph, gap, t, solutions, alpha)
        internal = evaluate_internal_energy(snapshot.graph, gap, t, solutions, alpha,
            energy_scale_J=w["observed"]["energy_scale_J"],
            boundary_contract="actual evolved gap and anchored current-carrying q held fixed in thermodynamic T derivative",
            boundary_temperature_derivatives=cutT)
        R0 = np.asarray([item["R0"] for item in roots])
        K0 = np.asarray([item["K0"] for item in roots])
        flux = np.asarray([item["flux"] for item in roots])
        Kth = np.asarray([item["Kth"] for item in roots])
        fluxth = np.asarray([item["fluxth"] for item in roots])
        Gth = reference.gap_gradient*w["observed"]["energy_scale_J"]
        Ith = reference.current_bar*snapshot.rows.moments.I0
        rows = snapshot.rows.evaluate(K_zero=K0, K_first=np.zeros_like(K0), K_thermal_zero=Kth,
            phi_V=np.zeros(snapshot.graph.n_nodes), thermal_gradient_J_per_gapbar=Gth,
            thermal_edge_current_A=Ith, edge_flux_zero=flux, edge_flux_first=np.zeros_like(flux),
            thermal_edge_flux_zero=fluxth, imposed_current_A=candidate["circuit_state_new_SI"][1],
            reference_contract="actual finite evolved gap/links, full evolved h, same bath and finite Matsubara thermal reference; local Wilson coordinate; K0 only")
        packed = rows.pack(neutrality_nodes=snapshot.free, force_nodes=snapshot.free,
            scales=snapshot.rows.scales)
        observer_energy = AdiabaticElectronicEnergy(snapshot.rows.moments.energy,
            snapshot.rows.moments.weights, snapshot.rows.moments.volume,
            snapshot.rows.moments.N0, h0)
        electronic = observer_energy.energy(internal.internal_energy_J, R0[:, :, 0, 0].real, h[:, :, 0])
        initial_h = np.asarray([sl.h for sl in snapshot.slices])
        initial_R0 = np.asarray([sl.Rbase.value for sl in snapshot.slices])
        initial_electronic = observer_energy.energy(w["observed"]["internal_energy_J"],
            initial_R0[:, :, 0, 0].real, initial_h[:, :, 0])
        phonon = phonon_energy_J(candidate["phonon_energy_J"], candidate["phonon_weights_J"],
            candidate["phonon_DOS_per_J_m3"], candidate["phonon_occupation_new"],
            snapshot.rows.moments.volume)
        ninitial = candidate.get("phonon_occupation_initial", data.get("phonon_occupation"))
        initial_phonon = None if ninitial is None else phonon_energy_J(
            candidate["phonon_energy_J"], candidate["phonon_weights_J"],
            candidate["phonon_DOS_per_J_m3"], ninitial, snapshot.rows.moments.volume)
        cm = snapshot.circuit.energy(candidate["circuit_state_new_SI"])
        initial_cm = snapshot.circuit.energy(snapshot.y)
        energy = dict(thermal_internal_new_J=internal.internal_energy_J,
            electronic_distribution_increment_new_J=electronic.distribution_increment_J,
            electronic_neutral_adiabatic_new_J=electronic.total_J,
            electronic_neutral_adiabatic_initial_J=initial_electronic.total_J,
            electronic_neutral_adiabatic_change_J=electronic.total_J-initial_electronic.total_J,
            phonon_new_J=phonon, phonon_initial_J=initial_phonon,
            circuit_new_J=cm, circuit_initial_J=initial_cm,
            total_neutral_adiabatic_new_J=electronic.total_J+phonon+cm,
            total_neutral_adiabatic_initial_J=None if initial_phonon is None else initial_electronic.total_J+initial_phonon+initial_cm,
            charge_gauge_completion_missing=True, K1_energy_not_included=True,
            normal_background="same fixed-bath normal-state background omitted; finite action reference held consistent",
            scope="independent neutral adiabatic energy diagnostic, not a conserved full Nambu energy")
        excess_energy = support_excess_energy_J(snapshot.rows.moments.energy,
            snapshot.rows.moments.weights, snapshot.rows.moments.volume,
            snapshot.rows.moments.N0, R0[:, :, 0, 0].real, h)
        retained = data.get("receipt", {}).get("retained_energy_J")
        energy.update(support_excess_energy_new_DOS_J=excess_energy,
            support_excess_fraction_retained_energy=None if retained is None or retained <= 0 else excess_energy/retained,
            support_excess_scope="energy-weighted diagnostic using actual new DOS; not a priori bound or acceptance")
        if "phonon_bath" in candidate and "tau_escape_s" in candidate:
            bath = phonon_energy_J(candidate["phonon_energy_J"], candidate["phonon_weights_J"],
                candidate["phonon_DOS_per_J_m3"], candidate["phonon_bath"], snapshot.rows.moments.volume)
            escape = (phonon-bath)/float(candidate["tau_escape_s"])
            energy.update(phonon_escape_out_new_W=escape,
                dt_times_phonon_escape_out_new_J=float(candidate["dt_s"])*escape)
        _save(output/"finite_state_fields.npz", gap_bar=gap, alpha=alpha, h=h, h_energy=hE,
            R0=R0, K0=K0, K_thermal_zero=Kth, edge_flux_zero=flux,
            thermal_edge_flux_zero=fluxth, thermal_gradient_J_per_gapbar=Gth,
            thermal_edge_current_A=Ith, zeroth_scaled_rows=packed,
            neutrality_covariant_zero_V=rows.neutrality_V,
            gap_force_zero_J_per_gapbar=rows.gap_force_J_per_gapbar,
            edge_current_zero_A=rows.edge_current_A,
            kinetic_full_residual_zero=np.asarray([item["residual"] for item in roots]))
        _write_json(output/"retarded_receipts.json", dict(records=[item["receipt"] for item in roots]))
        result.update(status="FINITE_STATE_ZERO_ORDER_REBUILT_DIAGNOSTIC_ONLY",
            finite_state_rebuilt=True, completed_queries=completed,
            energy=energy, zeroth_scaled_rows_max=float(np.max(abs(packed), initial=0.)),
            zeroth_port_current_residual_A=rows.port_current_residual_A.tolist(),
            zeroth_interior_current_divergence_max_A=float(np.max(abs(rows.interior_current_divergence_A), initial=0.)),
            thermal_modes=count, real_energy_modes=len(roots),
            thermal_maximum_spectral_residual=reference.maximum_spectral_residual,
            thermal_internal_energy_maximum_tangent_residual=internal.maximum_free_tangent_residual,
            distribution_energy_derivative="nonuniform three-point derivative of actual h_new at fixed E/Eu; no temperature fit",
            distribution_energy_derivative_convergence_admitted=False,
            fields_path=str(output/"finite_state_fields.npz"),
            independent_new_first_order_and_collision_rows_not_computed=True)
    except Exception as error:
        result.update(status="FINITE_STATE_REBUILD_REJECTED_OR_INCOMPLETE",
            error=repr(error), error_type=type(error).__name__)
    finally:
        _CONTEXT = None
    result["elapsed_s"] = time.monotonic()-started
    _write_json(output/"summary.json", result)
    return result
