"""Six simultaneous v1.1.3 latency controls; no change to the physical producer.

Two low-current escape20 references and matched low/high-current tests of tau_ee.
All cases retain eta=.8, the depairing-calibrated material and the full circuit.
One Numba thread per distinct physical core; exactly two cores per case.
"""
from __future__ import annotations

import argparse
import copy
import importlib.util
import json
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SHARED = ROOT / "pipelines/10_thermal_stiffness_basecases.py"


def load_driver():
    spec = importlib.util.spec_from_file_location("_latency_controls_shared", SHARED)
    driver = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(driver)
    driver.KNOWN_RUNNERS.update({"11_bias_escape_controls.py", Path(__file__).name})
    return driver


def allocation(driver, stage, *, anchor_threads=None, max_worker_cores=12,
               resources=None, active=None):
    """Admit all six slots together, or fail before starting any photon run."""
    if stage != "anchors" or anchor_threads not in (None, 2) or max_worker_cores != 12:
        raise ValueError("This campaign requires six simultaneous two-core slots")
    resources = driver.linux_resources() if resources is None else resources
    budget = driver.resource_budget(resources, .9, .9, 1)
    groups = resources["physical_core_groups"]
    allowed = set(budget["selected_affinity_cpus"])
    active = driver.active_known_processes() if active is None else active
    occupied = {cpu for row in active for cpu in row["affinity_cpus"]}
    free = {i for i, group in enumerate(groups)
            if set(group) <= allowed and not set(group) & occupied}
    slots = []
    for _ in range(6):
        start = next((i for i in sorted(free) if i+1 in free), None)
        if start is None:
            raise ValueError("Six simultaneous cases need twelve free physical cores; preserve other jobs and retry later")
        slots.append(dict(cpu_start=start, threads=2,
                          core_groups=[groups[start], groups[start+1]]))
        free -= {start, start+1}
    if not free:
        raise ValueError("An additional complete coordinator core is required")
    coordinator = groups[max(free)]
    assigned = set(coordinator)
    for slot in slots:
        for group in slot["core_groups"]:
            assigned.update(group)
    reserved = [group for group in groups if not set(group) & (assigned | occupied)]
    if len(reserved) < 2:
        raise ValueError("The shared allocation must reserve two whole physical cores")
    reservation = (3 + 4*6)*1024**3
    if reservation > budget["memory_limit_bytes"]:
        raise ValueError("Available memory cannot admit all six cases")
    return dict(inventory=resources, budget=budget, active_known_processes=active,
                coordinator_core_group=coordinator, slots=slots,
                global_reserved_core_groups=reserved, worker_physical_cores=12,
                memory_reservation_bytes=reservation,
                memory_reservation_scope="Conservative scheduling estimate, not an OS memory limit",
                allocation_scope="Six simultaneous cases; two distinct physical cores each; no SMT performance assumption")


def configure(driver, args):
    reference = driver.scratch_path(args.reference_campaign)
    definition = json.loads((reference/"campaign_plan.json").read_text())
    execution = json.loads((reference/"campaign_execution.json").read_text())
    driver.require_anchors(reference, definition, execution)
    manifests = {}
    for wavelength in (775, 1550):
        path = reference/"runs"/f"anchor425_eta0800_{wavelength}"/"manifest.json"
        m = json.loads(path.read_text())
        if (not m["run_complete"] or m["model_version"] != "1.1.3-candidate"
                or m["parameter_identity"]["stiffness_temperature"] != "electron"
                or m["parameter_identity"]["gradient_closure"] != "thermal-depairing"
                or m["source_pins_sha256"] != driver.pins("electron")):
            raise ValueError("Authenticated current v1.1.3 reference pair required")
        manifests[wavelength] = (path, m)
    base = dict(manifests[775][1]["parameter_identity"])
    for key, value in manifests[1550][1]["parameter_identity"].items():
        if key != "wavelength_nm" and value != base[key]:
            raise ValueError("Reference pair differs in a non-wavelength input: "+key)
    if (base["tau_ee_Tc_ps"] != 5. or base["tau_escape_ps"] != 20.
            or base["retained_fraction"] != .8 or base["current_uA"] != 21.5):
        raise ValueError("Expected the ee5/escape20/eta0800/21.5uA anchors")
    base.pop("threads")
    if not math.isfinite(args.tau_ee_Tc_ps) or not 0 < args.tau_ee_Tc_ps < 5:
        raise ValueError("The finite tau_ee test must be smaller than the 5 ps reference")
    if not math.isfinite(args.ceiling_ps) or args.ceiling_ps <= 0:
        raise ValueError("Positive finite observation ceiling required")
    cases = []
    for current, tau in ((15.5, 5.), (15.5, args.tau_ee_Tc_ps),
                         (21.5, args.tau_ee_Tc_ps)):
        label = ("bias"+format(current, "g").replace(".", "p")+"_escape20_ee"
                 +format(tau, "g").replace(".", "p")+"_eta0800")
        for wavelength in (775., 1550.):
            params = dict(base, current_uA=current, tau_ee_Tc_ps=tau,
                          tau_escape_ps=20., wavelength_nm=wavelength,
                          ceiling_ps=args.ceiling_ps)
            driver.RUNNER.validate_args(argparse.Namespace(
                **params, threads=2, cpu_start=0, output_dir=driver.SCRATCH/"unused_latency_validation",
                resume=False))
            cases.append(dict(name=f"{label}_{int(wavelength)}", pair=label,
                              stage="anchors", parameters=params))
    old_plan = driver.plan
    driver.catalog = lambda *, ceiling_ps=100.: copy.deepcopy(cases)

    def plan(campaign, *, ceiling_ps=100.):
        result = old_plan(campaign, ceiling_ps=ceiling_ps)
        result.update(schema="pysnspd.low_bias_latency_controls.v1",
                      orchestrator_sha256=driver.sha(__file__),
                      shared_campaign_driver_sha256=driver.sha(SHARED),
                      stage_order=["anchors"],
                      reference_manifests={str(w): dict(path=str(p), sha256=driver.sha(p))
                                           for w, (p, m) in manifests.items()},
                      commands={c["name"]: driver.command(c, campaign, threads=2, cpu_start=0)
                                for c in cases})
        result["scope"].pop("current_uA", None)
        result["scope"].update(
            bias_currents_uA=[15.5, 21.5], tau_escape_ps=20.,
            tau_ee_test_Tc_ps=args.tau_ee_Tc_ps,
            retained_fraction=.8, reference_campaign=str(reference),
            only_new_physical_control="tau_ee(Tc), shared by BGK and KWT mobility",
            no_per_bias_material_retuning=True,
            low_current_reference="ee5/escape20 for both wavelengths",
            high_current_reference="existing authenticated ee5/escape20 pair",
            hypothesis="A shorter condensate relaxation can reduce energy loss before instability, especially near the low-current detection threshold",
            experimental_comparison="Fixed-energy deterministic crossings are not IRF peak latencies; no jitter or fit claimed")
        result["scheduling"] = dict(
            simultaneous_cases=6, threads_per_case=2, worker_physical_cores=12,
            one_numba_thread_per_distinct_physical_core=True,
            affinity_assigned_at_runtime=True, coordinator_separate_core=True,
            minimum_whole_cores_reserved=2, no_silent_concurrency_fallback=True,
            no_background_service=True)
        result["gate"] = dict(reference_pair_authenticated=True,
                              each_case_original_terminal_audit_required=True,
                              no_quantitative_admission_inferred=True)
        return result

    driver.plan = plan
    driver.allocation = lambda stage, **kwargs: allocation(driver, stage, **kwargs)
    return cases


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference-campaign", type=Path, required=True)
    parser.add_argument("--campaign-dir", type=Path, required=True)
    parser.add_argument("--tau-ee-Tc-ps", type=float, default=4.)
    parser.add_argument("--ceiling-ps", type=float, default=100.)
    parser.add_argument("--threads-per-case", type=int, choices=(2,), default=2)
    parser.add_argument("--max-worker-cores", type=int, choices=(12,), default=12)
    parser.add_argument("--prepare-only", action="store_true")
    parser.add_argument("--resume", action="store_true")
    args = parser.parse_args(argv)
    driver = load_driver()
    cases = configure(driver, args)
    campaign = driver.scratch_path(args.campaign_dir)
    if args.prepare_only:
        definition = driver.plan(campaign, ceiling_ps=args.ceiling_ps)
        resources = driver.allocation("anchors", max_worker_cores=args.max_worker_cores)
        campaign.mkdir(parents=True, exist_ok=True)
        path = campaign/"campaign_plan.json"
        if path.exists() and json.loads(path.read_text()) != definition:
            raise ValueError("Existing plan differs; choose a fresh campaign directory")
        if not path.exists():
            driver.write_json(path, definition)
        driver.write_json(campaign/"preparation_resources.json",
                          dict(scope="Preparation only; runtime must recheck capacity",
                               resources=resources))
        print(json.dumps(dict(status="PREPARED", cases=len(cases), simultaneous_cases=6,
                              threads_per_case=2, worker_physical_cores=12,
                              reserved_whole_cores=resources["global_reserved_core_groups"],
                              case_names=[c["name"] for c in cases],
                              no_simulation_executed=True)), flush=True)
        return 0
    words = ["--campaign-dir", str(campaign), "--stage", "anchors",
             "--max-worker-cores", "12", "--ceiling-ps", str(args.ceiling_ps)]
    if args.resume:
        words.append("--resume")
    return driver.main(words)


if __name__ == "__main__":
    raise SystemExit(main())
