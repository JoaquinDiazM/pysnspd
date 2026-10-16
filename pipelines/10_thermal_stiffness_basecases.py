"""Staged 12-case thermal-rigidity candidate, with a shared foreground budget.

No --stage means prepare only. Anchors must pass finite-window source, capacity,
identity and energy-record audits before controls can start. These checks never
claim physical calibration or temporal/spatial convergence. No full reset.
"""
from __future__ import annotations

from concurrent.futures import ThreadPoolExecutor, as_completed
import argparse
import fcntl
import importlib.util
import json
import math
import os
from pathlib import Path
import queue
import signal
import subprocess
import sys
import threading

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from sandbox.stage4_core.parallel_runtime import linux_resources, resource_budget
SPEC = importlib.util.spec_from_file_location("_thermal_campaign_runner", ROOT/"pipelines/09_thermal_stiffness_photon_run.py")
RUNNER = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(RUNNER)
PYTHON = Path("/home/jdiaz/.conda/envs/snspd/bin/python")
SCRATCH = Path("/home/jdiaz/scratch")
TERMINAL = {"TRIGGER_MARGIN_COMPLETE", "TIME_CENSORED_NO_TRIGGER"}
KNOWN_RUNNERS = {"05_effective_photon_run.py", "08_effective_basecase_sweep.py",
                 "09_thermal_stiffness_photon_run.py", "10_thermal_stiffness_basecases.py"}


def sha(path):
    return RUNNER.BASE.sha256(path)


def write_json(path, value):
    return RUNNER._ORIGINAL_WRITE_JSON(path, value)


def scratch_path(value):
    result = Path(value).expanduser().resolve()
    if not result.is_relative_to(SCRATCH.resolve()):
        raise ValueError("All campaign products must stay under /home/jdiaz/scratch")
    return result


def catalog(*, ceiling_ps=100.):
    defaults = vars(RUNNER.parser().parse_args(["--output-dir", str(SCRATCH/"unused_plan_placeholder")]))
    common = {k: v for k, v in defaults.items() if k not in ("output_dir", "resume", "cpu_start", "threads")}
    common.update(ceiling_ps=float(ceiling_ps), length_nm=320., spacing_nm=5., current_uA=21.5, tau_escape_ps=20.,
                  retained_fraction=.8, sheet_resistance_ohm=425., stiffness_temperature="electron")
    cases = []
    def pair(name, stage, overrides):
        for wavelength in (775., 1550.):
            parameters = dict(common, **dict(overrides, wavelength_nm=wavelength))
            cases.append(dict(name=f"{name}_{int(wavelength)}", pair=name, stage=stage,
                              parameters=RUNNER.BASE.safe_json(parameters)))
    pair("anchor425_eta0800", "anchors", {})
    pair("bath425_eta0800", "controls", dict(stiffness_temperature="bath"))
    pair("electron608_eta0800", "controls", dict(sheet_resistance_ohm=608.))
    pair("electron425_eta0667", "controls", dict(retained_fraction=.667))
    pair("electron425_eta0900", "controls", dict(retained_fraction=.9))
    pair("electron425_eta0800_eehalf", "controls", dict(tau_ee_Tc_ps=2.5))
    return cases


def pins(mode):
    old = RUNNER._MODE
    try:
        RUNNER._MODE = mode
        return RUNNER.source_pins()
    finally:
        RUNNER._MODE = old


def command(case, campaign, *, threads, cpu_start, resume=False):
    result = [str(PYTHON), "-u", str(ROOT/"pipelines/09_thermal_stiffness_photon_run.py"),
              "--output-dir", str(campaign/"runs"/case["name"])]
    flags = {"bath_K": "--bath-K", "tau_ee_Tc_ps": "--tau-ee-Tc-ps", "tau_ep_Tc_ps": "--tau-ep-Tc-ps",
             "target_depairing_uA": "--target-depairing-uA", "current_uA": "--current-uA"}
    for key, value in case["parameters"].items():
        flag = flags.get(key, "--"+key.replace("_", "-"))
        if isinstance(value, bool):
            if value:
                result.append(flag)
        elif value is not None:
            result.extend((flag, str(value)))
    result.extend(("--threads", str(threads), "--cpu-start", str(cpu_start)))
    if resume:
        result.append("--resume")
    return result


def plan(campaign, *, ceiling_ps=100.):
    cases = catalog(ceiling_ps=ceiling_ps)
    return dict(schema="pysnspd.thermal_candidate_basecases.v1", model_version="1.1.2-candidate",
        campaign_dir=str(campaign), orchestrator_sha256=sha(__file__),
        source_pins_sha256={mode: pins(mode) for mode in ("electron", "bath")},
        cases=cases, stage_order=["anchors", "controls"],
        commands={case["name"]: command(case, campaign,
            threads=6 if case["stage"] == "anchors" else 1, cpu_start=0) for case in cases},
        scheduling=dict(maximum_worker_physical_cores=12, minimum_whole_cores_reserved=2,
            anchors="Two independent cases, 1-6 distinct physical cores each; prefer 6 when free global budget permits, preserve exact resume thread identity",
            controls="Ten independent cases, one physical core each; bounded queue",
            no_background_service=True, no_minimum_utilization_requirement=True),
        gate=dict(anchors_both_complete_audited=True, global_energy_record_fraction_max=.01,
            source_energy_relative_defect_max=1e-10, positive_effective_capacity_required=True,
            accepted_local_energy_guard_required=True, original_outputs_authenticated=True,
            censored_anchor="Permitted only at declared finite ceiling; latency remains unassigned"),
        scope=dict(length_nm=320., width_nm=80., spacing_nm=5., current_uA=21.5,
            tau_escape_ps=20., fixed_sigma_nm=7.549834435, target_depairing_uA=38.8*80/120,
            retention_eta0800="Exploratory assumed retained fraction; not fitted to a four-ps latency",
            time_origin="effective post-cascade transfer", no_full_reset=True, no_statistical_jitter=True,
            physical_admission=False, new_convergence_not_inherited=True,
            performance_scope="No campaign ETA inferred from the first femtoseconds; short stiff-branch timing is not full-window throughput",
            thermal_control="Bath mode is a calibrated-D control, not the electron-temperature candidate"))


def known_runner(words):
    """Count the Python owner, not its timeout/taskset/shell wrappers."""
    decoded = [os.fsdecode(word) for word in words if word]
    if not decoded or not Path(decoded[0]).name.startswith("python"):
        return None
    return next((Path(word).name for word in decoded[1:]
                 if Path(word).name in KNOWN_RUNNERS), None)


def active_known_processes(exclude=()):
    excluded = set(exclude)|{os.getpid()}
    records = []
    for directory in Path("/proc").iterdir():
        if not directory.name.isdigit() or int(directory.name) in excluded:
            continue
        try:
            if directory.stat().st_uid != os.getuid():
                continue
            words = (directory/"cmdline").read_bytes().split(b"\0")
            matched = known_runner(words)
            if matched:
                records.append(dict(pid=int(directory.name), runner=matched,
                                    affinity_cpus=sorted(os.sched_getaffinity(int(directory.name)))))
        except (FileNotFoundError, ProcessLookupError, PermissionError):
            continue
    return records


def allocation(stage, *, anchor_threads=None, max_worker_cores=12, resources=None, active=None):
    resources = linux_resources() if resources is None else resources
    budget = resource_budget(resources, .9, .9, 1)
    groups = resources["physical_core_groups"]
    allowed = set(budget["selected_affinity_cpus"])
    admitted = [(index, group) for index, group in enumerate(groups) if set(group) <= allowed]
    if len(groups)-len(admitted) < 2:
        raise ValueError("Runtime topology cannot reserve two complete physical cores within the 90% budget")
    active = active_known_processes() if active is None else active
    occupied = {cpu for item in active for cpu in item["affinity_cpus"]}
    occupied_physical = sum(bool(set(group)&occupied) for group in groups)
    available_worker_count = min(max_worker_cores, len(admitted)-occupied_physical-1)
    if available_worker_count < 1:
        raise ValueError("Known active jobs already consume the shared physical-core ceiling")
    free = [(index, group) for index, group in admitted if not set(group)&occupied]
    if len(free) < 2:
        raise ValueError("Active known runs leave no complete coordinator/worker budget")
    coordinator = free[-1]
    worker_indices = {index for index, group in free[:-1]}
    slots = []
    if stage == "anchors":
        options = (anchor_threads,) if anchor_threads is not None else (6, 5, 4, 3, 2, 1)
        for threads in options:
            if not 1 <= threads <= 6 or 2*threads > available_worker_count:
                continue
            for possible_coordinator in reversed(free):
                remaining = {index for index, group in free if index != possible_coordinator[0]}
                trial = []
                for _ in range(2):
                    start = next((i for i in sorted(remaining)
                                  if set(range(i, i+threads)) <= remaining), None)
                    if start is None:
                        break
                    trial.append(dict(cpu_start=start, threads=threads,
                                      core_groups=[groups[i] for i in range(start, start+threads)]))
                    remaining -= set(range(start, start+threads))
                if len(trial) == 2:
                    coordinator = possible_coordinator
                    slots = trial
                    break
            if slots:
                break
        if not slots:
            raise ValueError("Anchors need two independent worker slots plus coordinator; preserve active runs and retry after capacity is free")
    else:
        slots = [dict(cpu_start=i, threads=1, core_groups=[groups[i]])
                 for i in sorted(worker_indices)[:min(10, available_worker_count)]]
    assigned_cpus = set(coordinator[1])
    for slot in slots:
        for group in slot["core_groups"]:
            assigned_cpus.update(group)
    global_reserved = [group for group in groups if not set(group)&(assigned_cpus|occupied)]
    if len(global_reserved) < 2:
        raise ValueError("Shared allocation failed to leave two whole physical cores free")
    reservation = (3+4*len(slots))*1024**3
    if reservation > budget["memory_limit_bytes"]:
        if stage == "anchors":
            raise ValueError("Available memory cannot admit two anchor processes")
        fit = max(0, int((budget["memory_limit_bytes"]-3*1024**3)//(4*1024**3)))
        slots = slots[:fit]
        assigned_cpus = set(coordinator[1])
    for slot in slots:
        for group in slot["core_groups"]:
            assigned_cpus.update(group)
    global_reserved = [group for group in groups if not set(group)&(assigned_cpus|occupied)]
    if len(global_reserved) < 2:
        raise ValueError("Shared allocation failed to leave two whole physical cores free")
    reservation = (3+4*len(slots))*1024**3
    if not slots:
        raise ValueError("Available memory cannot fit a candidate worker")
    return dict(inventory=resources, budget=budget, active_known_processes=active,
        coordinator_core_group=coordinator[1], slots=slots, global_reserved_core_groups=global_reserved,
        worker_physical_cores=sum(slot["threads"] for slot in slots),
        memory_reservation_bytes=reservation, memory_reservation_scope="Scheduling estimate, not an OS memory limit",
        allocation_scope="Known 05/08/09/10 affinity counted in one budget; unrelated CPU activity is not a guaranteed utilization measure")


def audit_case(case, campaign, definition, *, expected_threads):
    run = campaign/"runs"/case["name"]
    manifest = json.loads((run/"manifest.json").read_text())
    identity = dict(case["parameters"], threads=expected_threads)
    mode = case["parameters"]["stiffness_temperature"]
    expected_pins = definition["source_pins_sha256"][mode]
    if manifest.get("parameter_identity") != identity or manifest.get("source_pins_sha256") != expected_pins:
        raise ValueError("Terminal candidate identity/source mismatch: "+case["name"])
    if not manifest.get("run_complete") or manifest.get("status") not in TERMINAL:
        raise ValueError("Successful declared-window terminal required: "+case["name"])
    if manifest.get("model_version") != "1.1.2-candidate" or manifest.get("closure", {}).get("stiffness_temperature") != mode:
        raise ValueError("Candidate closure metadata mismatch")
    hashes = manifest.get("output_sha256", {})
    required = {"trajectory.npz", "snapshots.npz", "checkpoint_latest.npz", "dc/effective_dc.npz", "material.npz", "material.json"}
    if not required <= set(hashes):
        raise ValueError("Terminal numerical output hashes missing")
    for name, expected in hashes.items():
        path = (run/name).resolve()
        if not path.is_relative_to(run.resolve()) or sha(path) != expected:
            raise ValueError("Terminal output hash mismatch: "+name)
    calibration = json.loads((run/"calibration.json").read_text())
    physics = json.loads((run/"candidate_physics_identity.json").read_text())
    if calibration != manifest["physicalparams"]["diffusion_calibration"] or physics["calibration"] != calibration:
        raise ValueError("Saved calibration identity differs")
    if physics["parameter_identity"] != identity or physics["source_pins_sha256"] != expected_pins:
        raise ValueError("Saved candidate physics identity differs")
    observation = manifest["observation"]
    for key, expected in (("threshold_voltage_V", identity["threshold_mV"]*1e-3),
                          ("confirmation_time_s_required", identity["confirmation_ps"]*1e-12),
                          ("extra_time_s", identity["extra_time_ps"]*1e-12),
                          ("no_trigger_ceiling_s", identity["ceiling_ps"]*1e-12)):
        if not math.isclose(float(observation[key]), expected, rel_tol=8*sys.float_info.epsilon, abs_tol=1e-27):
            raise ValueError("Readout threshold/window identity mismatch: "+key)
    if observation.get("full_reset_requested") is not False:
        raise ValueError("Unexpected full-reset request")
    if manifest["status"] == "TIME_CENSORED_NO_TRIGGER":
        if observation.get("trigger_confirmed") or observation.get("first_crossing_time_s") is not None:
            raise ValueError("Censored run must not have an assigned detection latency")
    elif observation.get("trigger_confirmed") is not True:
        raise ValueError("Triggered terminal lacks a confirmed crossing")
    source = manifest["source"]
    retained = float(source["retained_energy_J"])
    if retained <= 0 or abs(float(source["relative_total_energy_defect"])) > 1e-10:
        raise ValueError("Initial retained-energy receipt failed")
    if abs(float(source["maximum_nodal_energy_defect_J"]))/retained > 1e-10:
        raise ValueError("Initial nodal-energy receipt failed")
    capacity = None
    if mode == "electron":
        capacity = float(manifest["thermal_candidate_diagnostics"]["min_successful_rhs_evaluations"])
        if not math.isfinite(capacity) or capacity <= 0 or not source.get("positive_capacity_checked_throughout_preparation"):
            raise ValueError("Positive effective-capacity receipts missing or failed")
        if source["minimum_preparation_effective_capacity_ratio"] <= 0:
            raise ValueError("Initial positive effective-capacity guard failed")
    guard = manifest["local_energy_control"]
    maximum_local = float(guard["max_guarded_accepted_local_energy_fraction"])
    if not math.isfinite(maximum_local) or maximum_local > manifest["numerical"]["local_energy_tolerance"]*(1+1e-8):
        raise ValueError("Accepted local-energy guard failed")
    np = RUNNER.np
    with np.load(run/"trajectory.npz", allow_pickle=False) as data:
        t, energy, external = (data[name] for name in ("time_s", "energy_total_J", "external_energy_integral_J"))
        if len(t) < 2 or any(not np.all(np.isfinite(value)) for value in (t, energy, external)) or np.any(np.diff(t) <= 0):
            raise ValueError("Finite increasing trajectory required")
        if manifest["status"] == "TIME_CENSORED_NO_TRIGGER":
            target_end = identity["ceiling_ps"]*1e-12
            if not math.isclose(float(t[-1]), target_end, rel_tol=1e-9, abs_tol=1e-18):
                raise ValueError("Censored anchor did not reach its declared finite ceiling")
        else:
            crossing = float(observation["first_crossing_time_s"])
            confirmed = float(observation["trigger_time_s"])
            if not all(math.isfinite(value) for value in (crossing, confirmed)) or confirmed < crossing+identity["confirmation_ps"]*1e-12-1e-18 or t[-1] < crossing+identity["extra_time_ps"]*1e-12-1e-18:
                raise ValueError("Triggered terminal lacks confirmation or full post-crossing margin")
        defect = energy-energy[0]-external
        denominator = retained+np.r_[0., np.cumsum(abs(np.diff(external)))]
        fraction = float(np.max(abs(defect)/denominator))
        if fraction > .01:
            raise ValueError("Global stored energy/work record exceeds 1% envelope")
    result = dict(schema="pysnspd.thermal_candidate_terminal_audit.v1", case=case["name"],
        all_requested_checks_passed=True, physical_admission=False,
        manifest_sha256=sha(run/"manifest.json"), calibration_sha256=sha(run/"calibration.json"),
        physics_identity_sha256=sha(run/"candidate_physics_identity.json"),
        source_initial_state_sha256=sha(run/"source_initial_state.npz"),
        maximum_global_energy_record_fraction=fraction, maximum_accepted_local_energy_fraction=maximum_local,
        minimum_successful_rhs_capacity_ratio=capacity,
        scope="Authenticated terminal, source, capacity and stored energy/work records; no new convergence or independent microscopic validation")
    target = campaign/"audits"/(case["name"]+".json")
    target.parent.mkdir(exist_ok=True)
    if target.exists():
        if json.loads(target.read_text()) != result:
            raise ValueError("Preserved terminal audit differs; original receipt retained")
    else:
        write_json(target, result)
    return dict(path=str(target), sha256=sha(target), manifest_sha256=result["manifest_sha256"])


def require_anchors(campaign, definition, receipt):
    for case in definition["cases"]:
        if case["stage"] != "anchors":
            continue
        item = receipt.get("cases", {}).get(case["name"], {})
        if item.get("status") != "COMPLETE_AUDITED":
            raise ValueError("Controls require both completed audited anchors; missing "+case["name"])
        audit = item["audit"]
        audit_path = scratch_path(audit["path"])
        if sha(audit_path) != audit["sha256"] or sha(campaign/"runs"/case["name"]/"manifest.json") != audit["manifest_sha256"]:
            raise ValueError("Recorded anchor audit changed")
        audit_case(case, campaign, definition, expected_threads=item["threads"])


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--campaign-dir", type=Path, required=True)
    parser.add_argument("--stage", choices=("anchors", "controls"))
    parser.add_argument("--prepare-only", action="store_true")
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--ceiling-ps", type=float, default=100.)
    parser.add_argument("--max-worker-cores", type=int, default=12)
    args = parser.parse_args(argv)
    if not 1 <= args.max_worker_cores <= 12:
        raise ValueError("max-worker-cores must lie in 1..12")
    campaign = scratch_path(args.campaign_dir)
    if not math.isfinite(args.ceiling_ps) or args.ceiling_ps<=0:
        raise ValueError("Positive finite ceiling-ps required")
    definition = plan(campaign,ceiling_ps=args.ceiling_ps)
    campaign.mkdir(parents=True, exist_ok=True)
    with (campaign/".execution.lock").open("a+") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX|fcntl.LOCK_NB)
        planpath = campaign/"campaign_plan.json"
        if planpath.exists():
            if json.loads(planpath.read_text()) != definition:
                raise ValueError("Existing plan/source identity differs; preserve directory and use a fresh one")
        else:
            write_json(planpath, definition)
        if args.prepare_only or args.stage is None:
            print(json.dumps(dict(status="PREPARED", cases=12, stage_order=["anchors", "controls"],
                campaign_dir=str(campaign), resources=linux_resources(), active=active_known_processes())), flush=True)
            return 0
        path = campaign/"campaign_execution.json"
        receipt = json.loads(path.read_text()) if path.exists() else dict(cases={})
        if receipt.get("plan_sha256", sha(planpath)) != sha(planpath):
            raise ValueError("Execution receipt belongs to a different authenticated plan")
        selected = [case for case in definition["cases"] if case["stage"] == args.stage]
        if args.stage == "controls":
            require_anchors(campaign, definition, receipt)
        if any(case["name"] in receipt["cases"] for case in selected) and not args.resume:
            raise ValueError("Existing stage progress requires --resume")
        existing_threads = {receipt["cases"][case["name"]]["threads"] for case in selected
                            if case["name"] in receipt["cases"]}
        if args.stage == "anchors" and len(existing_threads) > 1:
            raise ValueError("Anchor thread identities disagree")
        resources = allocation(args.stage, anchor_threads=next(iter(existing_threads), None),
                               max_worker_cores=args.max_worker_cores)
        own_cores = set(resources["coordinator_core_group"])
        for slot in resources["slots"]:
            for group in slot["core_groups"]:
                own_cores.update(group)
        full_affinity = set(os.sched_getaffinity(0))
        slots = queue.Queue()
        for slot in resources["slots"]:
            slots.put(slot)
        mutex = threading.RLock()
        stop = threading.Event()
        children = {}
        receipt.update(plan_sha256=sha(planpath), stage=args.stage, status="RUNNING", resources=resources)
        write_json(path, receipt)
        env = dict(os.environ, PYTHONDONTWRITEBYTECODE="1",
            PYTHONPYCACHEPREFIX=str(SCRATCH/"tmp/pysnspd/pycache"),
            NUMBA_CACHE_DIR=str(SCRATCH/"tmp/pysnspd/numba_cache"),
            MPLCONFIGDIR=str(SCRATCH/"tmp/pysnspd/matplotlib"), TMPDIR=str(SCRATCH/"tmp/pysnspd"))
        for key in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS", "BLIS_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
            env[key] = "1"
        (campaign/"logs").mkdir(exist_ok=True)
        def persist(name, value):
            with mutex:
                receipt["cases"][name] = value
                write_json(path, receipt)
        def terminate(signum, frame):
            stop.set()
            with mutex:
                active_children = list(children.values())
            for child in active_children:
                try:
                    child.send_signal(signal.SIGTERM)
                except ProcessLookupError:
                    pass
        previous = {signum: signal.signal(signum, terminate) for signum in (signal.SIGINT, signal.SIGTERM)}
        os.sched_setaffinity(0, set(resources["coordinator_core_group"]))
        def work(case):
            name = case["name"]
            item = receipt["cases"].get(name, {})
            if item.get("status") == "COMPLETE_AUDITED":
                audit_case(case, campaign, definition, expected_threads=item["threads"])
                return True
            slot = slots.get()
            try:
                if stop.is_set():
                    return False
                with mutex:
                    known = active_known_processes(exclude=children)
                    if any(set(row["affinity_cpus"])&own_cores for row in known):
                        raise RuntimeError("Another known campaign entered the reserved CPU set; refusing oversubscription")
                    if pins(case["parameters"]["stiffness_temperature"]) != definition["source_pins_sha256"][case["parameters"]["stiffness_temperature"]]:
                        raise ValueError("Candidate source changed during campaign; refusing new work")
                    resumed = (campaign/"runs"/name).exists()
                    if resumed and not args.resume:
                        raise ValueError("Existing run requires --resume")
                    if item and item["threads"] != slot["threads"]:
                        raise ValueError("Resume numerical thread identity differs")
                    cmd = command(case, campaign, threads=slot["threads"], cpu_start=slot["cpu_start"], resume=resumed)
                    persist(name, dict(status="RUNNING", threads=slot["threads"], command=cmd))
                    launched = ["taskset", "-c", ",".join(map(str, sorted(full_affinity))), *cmd]
                    child = subprocess.Popen(launched, cwd=ROOT, env=env,
                        stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, bufsize=1)
                    children[child.pid] = child
                with (campaign/"logs"/(name+".log")).open("a") as stream:
                    stream.write("\nCOMMAND "+json.dumps(cmd)+"\n");stream.flush()
                    for line in child.stdout:
                        stream.write(line);stream.flush()
                        print(json.dumps(dict(case=name, progress=line.rstrip())), flush=True)
                result = child.wait()
                with mutex:
                    children.pop(child.pid, None)
                if result != 0 or stop.is_set():
                    raise RuntimeError("Runner interrupted or failed; exit="+str(result))
                audit = audit_case(case, campaign, definition, expected_threads=slot["threads"])
                persist(name, dict(status="COMPLETE_AUDITED", threads=slot["threads"], audit=audit))
                return True
            except Exception as exc:
                persist(name, dict(status="INTERRUPTED" if stop.is_set() else "FAILED",
                                  threads=slot["threads"], error=type(exc).__name__+": "+str(exc)))
                print(json.dumps(dict(case=name, status="FAILED", error=str(exc))), flush=True)
                return False
            finally:
                slots.put(slot)
        global_path = SCRATCH/"thermal_candidate_runtime.lock"
        try:
            with global_path.open("a+") as global_lock:
                fcntl.flock(global_lock, fcntl.LOCK_EX|fcntl.LOCK_NB)
                with ThreadPoolExecutor(max_workers=len(resources["slots"])) as executor:
                    futures = [executor.submit(work, case) for case in selected]
                    results = [future.result() for future in as_completed(futures)]
            receipt.update(status="INTERRUPTED" if stop.is_set() else
                           "COMPLETE_AUDITED" if all(results) else "COMPLETE_WITH_FAILURES")
            write_json(path, receipt)
            print(json.dumps(dict(status=receipt["status"], stage=args.stage,
                                  completed=sum(results), total=len(selected))), flush=True)
            return 0 if all(results) and not stop.is_set() else 1
        finally:
            os.sched_setaffinity(0, full_affinity)
            for signum, handler in previous.items():
                signal.signal(signum, handler)


if __name__ == "__main__":
    raise SystemExit(main())