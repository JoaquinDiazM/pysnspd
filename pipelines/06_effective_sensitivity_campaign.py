"""Prepare or execute a finite, one-input-at-a-time effective SNSPD campaign.

Execution is sequential and foreground. This orchestrator does not modify the
numerical model, its accuracy profile, or the trigger-plus-margin controller.
"""
from __future__ import annotations

import argparse
from contextlib import contextmanager
import fcntl
import datetime as dt
import hashlib
import json
import os
from pathlib import Path
import shlex
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
SCRATCH = Path("/home/jdiaz/scratch")
PYTHON = Path("/home/jdiaz/.conda/envs/snspd/bin/python")
RUNNER = ROOT / "pipelines/05_effective_photon_run.py"
AUDITOR = ROOT / "sandbox/stage5_effective_validation/audit_runs.py"
IDEP_REFERENCE_UA = 30.36744305465
SIGMA_REFERENCE_NM = 7.549834435
COMMON = dict(length_nm=320., spacing_nm=5., xi_cutoff=32., xi_order=6,
    phonon_nodes=96, transport_order=48, wavelength_nm=775.,
    sigma_nm=SIGMA_REFERENCE_NM, retained_fraction=.667, dt_ps=.02,
    ceiling_ps=200., threshold_mV=.1, confirmation_ps=.5, extra_time_ps=5.,
    tau_ee_Tc_ps=5., tau_ep_Tc_ps=24.7, tau_escape_ps=15.,
    current_uA=21.29720614489798, bath_K=.9, dc_only=False, duration_ps=None,
    snapshot_every_ps=.1, checkpoint_every_ps=1., checkpoint_wall_s=60.,
    step_tolerance=2e-7, local_energy_tolerance=2.5e-6,
    max_iterations=40, max_rejections=50, threads=4)
PHYSICAL_NAMES = ("effective_bcs_functional.py", "effective_bcs_preparation.py",
    "effective_bcs_kinetics.py", "effective_bcs_ports.py", "effective_coupled.py",
    "effective_material.py", "electrical_ports.py", "thermal_spatial_usadel.py",
    "kinetic_events.py", "trigger_window.py", "effective_bcs_kinetics_numba.py")
TERMINAL_STATUSES = {"TRIGGER_MARGIN_COMPLETE", "TIME_CENSORED_NO_TRIGGER"}


def sha256(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024**2), b""):
            h.update(block)
    return h.hexdigest()


def write_json(path, value):
    path = Path(path)
    temporary = path.with_name(path.name + ".part")
    temporary.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")
    os.replace(temporary, path)


def scratch_path(value):
    path = Path(value).expanduser().resolve()
    if not path.is_relative_to(SCRATCH.resolve()):
        raise ValueError("Campaign data, logs, and plots must remain in scratch")
    return path


def source_pins():
    paths = [ROOT/"pysnspd/experimental"/name for name in PHYSICAL_NAMES]
    paths += [RUNNER, ROOT/"sandbox/stage4_core/parallel_runtime.py"]
    return {str(path.relative_to(ROOT)): sha256(path) for path in paths}


def case_catalog():
    """Fixed exploratory probes; these are not measured parameter intervals."""
    specs = [
        ("nominal775", "bias", "775 nm nominal", {}),
        ("bias065_775", "bias", "775 nm; I=0.65 I_dep", {"current_uA":.65*IDEP_REFERENCE_UA}),
        ("bias075_775", "bias", "775 nm; I=0.75 I_dep", {"current_uA":.75*IDEP_REFERENCE_UA}),
        ("nominal1550", "bias", "1550 nm nominal", {"wavelength_nm":1550.}),
        ("bias075_1550", "bias", "1550 nm; I=0.75 I_dep",
            {"wavelength_nm":1550., "current_uA":.75*IDEP_REFERENCE_UA}),
        ("ee_half775", "rates", "775 nm; tau_ee=2.5 ps", {"tau_ee_Tc_ps":2.5}),
        ("ee_double775", "rates", "775 nm; tau_ee=10 ps", {"tau_ee_Tc_ps":10.}),
        ("ep_half775", "rates", "775 nm; tau_ep=12.35 ps", {"tau_ep_Tc_ps":12.35}),
        ("ep_double775", "rates", "775 nm; tau_ep=49.4 ps", {"tau_ep_Tc_ps":49.4}),
        ("escape_half775", "rates", "775 nm; tau_esc=7.5 ps", {"tau_escape_ps":7.5}),
        ("escape_double775", "rates", "775 nm; tau_esc=30 ps", {"tau_escape_ps":30.}),
        ("retention050_775", "source", "775 nm; eta=0.5", {"retained_fraction":.5}),
        ("retention080_775", "source", "775 nm; eta=0.8", {"retained_fraction":.8}),
        ("sigma_narrow775", "source", "775 nm; sigma/sqrt(2)",
            {"sigma_nm":SIGMA_REFERENCE_NM/2**.5}),
        ("sigma_wide775", "source", "775 nm; sigma*sqrt(2)",
            {"sigma_nm":SIGMA_REFERENCE_NM*2**.5}),
        ("retention080_1550", "source", "1550 nm; eta=0.8",
            {"wavelength_nm":1550., "retained_fraction":.8}),
    ]
    result = []
    for name, group, label, overrides in specs:
        parameters = dict(COMMON, **overrides)
        reference = "nominal1550" if parameters["wavelength_nm"] == 1550. else "nominal775"
        changes = {key:value for key,value in parameters.items()
            if value != dict(COMMON, wavelength_nm=parameters["wavelength_nm"])[key]}
        result.append(dict(name=name, group=group, label=label, parameters=parameters,
            reference=reference, changed_input=changes))
    return result


def selected_cases(group):
    catalog = case_catalog()
    if group == "all":
        return catalog
    names = {case["name"] for case in catalog if case["group"] == group}
    for case in catalog:
        if case["name"] in names:
            names.add(case["reference"])
    return [case for case in catalog if case["name"] in names]


def runner_command(case, campaign_dir, *, resume=False):
    args = [str(PYTHON), "-u", str(RUNNER), "--output-dir",
        str(Path(campaign_dir)/"runs"/case["name"])]
    for key, value in case["parameters"].items():
        if value is None or value is False:
            continue
        option = "--" + key.replace("_", "-")
        args.append(option)
        if value is not True:
            args.append(str(value))
    args += ["--cpu-start", "0"]
    if resume:
        args.append("--resume")
    return args


def build_plan(campaign_dir):
    cases = case_catalog()
    return dict(schema="pysnspd.effective_sensitivity_campaign.v1",
        campaign_dir=str(campaign_dir), source_pins_sha256=source_pins(),
        orchestrator_sha256=sha256(__file__), auditor_sha256=sha256(AUDITOR),
        cases=cases, commands={case["name"]:runner_command(case,campaign_dir) for case in cases},
        resources=dict(simultaneous_runs=1, numba_threads=4,
            blas_threads=1, global_ceiling_fraction=.9,
            actual_cpu_inventory_rechecked_by_runner_for_every_case=True),
        scope=dict(model="neutral effective BCS/KWT", fresh_dc_per_case=True,
            exploratory_parameter_points_not_experimental_intervals=True,
            time_origin="effective post-cascade transfer, not optical absorption",
            observation="positive baseline-subtracted Vout; threshold0.1mV; confirmation0.5ps; stop crossing+5ps; no-trigger ceiling200ps",
            no_full_reset=True, no_statistical_jitter=True,
            reference_Idep_uA=IDEP_REFERENCE_UA,
            reference_Idep_scope="effective bulk prediction for the fixed7nm material; not measured critical current",
            tau_ee_effects="BGK thermalization and KWT scattering mobility vary together",
            tau_ep_effects="acoustic alpha2F collision normalization and KWT scattering mobility vary together",
            source_effects="common FD/BE initialization changes spatial width AND electron/phonon partition when sigma or eta changes",
            spectral_and_time_convergence_not_inherited_for_new_physical_cases=True))


def check_case_identity(case, run_dir, pins, *, require_complete):
    """Reject stale parameters, changed code, incomplete or corrupted terminals."""
    run_dir = Path(run_dir)
    manifest = json.loads((run_dir/"manifest.json").read_text())
    if manifest.get("parameter_identity") != case["parameters"]:
        raise ValueError("Case parameter identity differs: " + case["name"])
    if manifest.get("source_pins_sha256") != pins:
        raise ValueError("Case physical sources differ: " + case["name"])
    if require_complete:
        if manifest.get("run_complete") is not True or manifest.get("status") not in TERMINAL_STATUSES:
            raise ValueError("Successful finite-window terminal required: " + case["name"])
        hashes = manifest.get("output_sha256", {})
        required = {"trajectory.npz","snapshots.npz","checkpoint_latest.npz",
            "dc/effective_dc.npz","material.npz","material.json"}
        if not required <= set(hashes):
            raise ValueError("Terminal output hashes missing")
        for name, expected in hashes.items():
            path = (run_dir/name).resolve()
            if not path.is_relative_to(run_dir.resolve()) or sha256(path) != expected:
                raise ValueError("Terminal output changed: " + name)
    return manifest


def child_environment():
    env = os.environ.copy()
    for key in ("OPENBLAS_NUM_THREADS","OMP_NUM_THREADS","MKL_NUM_THREADS","NUMEXPR_NUM_THREADS"):
        env[key] = "1"
    base = SCRATCH/"tmp/pysnspd"
    for key, leaf in (("PYTHONPYCACHEPREFIX","pycache"),("NUMBA_CACHE_DIR","numba_cache"),
            ("MPLCONFIGDIR","matplotlib")):
        path = base/leaf
        path.mkdir(parents=True,exist_ok=True)
        env[key] = str(path)
    env["TMPDIR"] = str(base)
    env["NUMBA_NUM_THREADS"] = "4"
    return env


def run_logged(command, logfile, env, *, cpu=None):
    """Stream visible foreground progress and preserve every session's log."""
    with Path(logfile).open("a") as stream:
        stream.write("\nCOMMAND "+shlex.join(command)+"\n")
        stream.flush()
        action = None if cpu is None else lambda:os.sched_setaffinity(0,{cpu})
        with subprocess.Popen(command, cwd=ROOT, env=env, stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT, text=True, bufsize=1, preexec_fn=action) as child:
            try:
                for line in child.stdout:
                    print(line, end="", flush=True)
                    stream.write(line)
                    stream.flush()
                code = child.wait()
            except BaseException:
                child.terminate()
                try:
                    child.wait(timeout=75)
                except subprocess.TimeoutExpired:
                    # Preserve the process and its latest checkpoint; never SIGKILL a state write.
                    print("Child still shutting down; check foreground PID "+str(child.pid),flush=True)
                raise
    if code:
        raise subprocess.CalledProcessError(code, command)


def audit_case(case, campaign_dir, env):
    sys.path.insert(0,str(ROOT))
    from sandbox.stage4_core.parallel_runtime import linux_resources, resource_budget
    budget = resource_budget(linux_resources(),.9,.9,1)
    stamp = dt.datetime.now(dt.timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
    audit_dir = campaign_dir/"audits"/(case["name"]+"_"+stamp)
    command = [str(PYTHON),str(AUDITOR),"--run",
        case["name"]+"="+str(campaign_dir/"runs"/case["name"]),
        "--output-dir",str(audit_dir)]
    run_logged(command,campaign_dir/"logs"/(case["name"]+".audit.log"),env,cpu=budget["coordinator_cpu"])
    report = json.loads((audit_dir/"audit.json").read_text())
    if report.get("all_requested_checks_passed") is not True:
        raise ValueError("Saved-state audit failed for "+case["name"]+"; inspect "+str(audit_dir))
    return dict(path=str(audit_dir/"audit.json"), sha256=sha256(audit_dir/"audit.json"),
        passed=True, manifest_sha256=sha256(campaign_dir/"runs"/case["name"]/"manifest.json"))


def aggregate_evidence(plan, receipt, campaign_dir, selected, group):
    """Authenticate finite terminals and individual audits; no physics admission."""
    authenticated = {}
    catalog = {case["name"]:case for case in plan["cases"]}
    for name,item in receipt["cases"].items():
        if name not in catalog or item.get("status") != "COMPLETE_AUDITED":
            raise ValueError("Unexpected execution entry")
        run_dir = campaign_dir/"runs"/name
        check_case_identity(catalog[name],run_dir,plan["source_pins_sha256"],require_complete=True)
        if sha256(run_dir/"manifest.json") != item["manifest_sha256"]:
            raise ValueError("Recorded terminal manifest changed: "+name)
        path = scratch_path(item["audit"]["path"])
        if sha256(path) != item["audit"]["sha256"]:
            raise ValueError("Recorded audit changed: "+name)
        audit = json.loads(path.read_text())
        if audit.get("all_requested_checks_passed") is not True:
            raise ValueError("Recorded audit failed: "+name)
        authenticated[name] = dict(manifest_sha256=item["manifest_sha256"],
            audit_sha256=item["audit"]["sha256"],status="COMPLETE_AUDITED")
    selected_names = [case["name"] for case in selected]
    if not set(selected_names) <= set(authenticated):
        raise ValueError("Selected group is incomplete")
    return dict(schema="pysnspd.effective_sensitivity_aggregate.v1",group=group,
        plan_sha256=sha256(campaign_dir/"campaign_plan.json"),
        selected_names=selected_names,selected_group_complete_audited=True,
        all16_complete_audited=set(authenticated)==set(catalog),
        authenticated_cases=authenticated,
        scope="Individual finite-window conservation/support/sign/trigger audits; no new parameter convergence or experimental admission")


@contextmanager
def execution_lock(campaign_dir):
    """One campaign writer/foreground child at a time, including --resume."""
    path = Path(campaign_dir)/".execution.lock"
    with path.open("a+") as stream:
        try:
            fcntl.flock(stream.fileno(),fcntl.LOCK_EX|fcntl.LOCK_NB)
        except BlockingIOError:
            raise RuntimeError("Another invocation owns this campaign; do not execute groups concurrently") from None
        try:
            stream.seek(0)
            stream.truncate()
            stream.write(json.dumps(dict(pid=os.getpid(),host=os.uname().nodename))+"\n")
            stream.flush()
            yield
        finally:
            fcntl.flock(stream.fileno(),fcntl.LOCK_UN)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir",type=Path,required=True)
    parser.add_argument("--group",choices=("bias","rates","source","all"),default="bias")
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--plan-only",action="store_true")
    mode.add_argument("--execute",action="store_true")
    parser.add_argument("--resume",action="store_true",
        help="Require the existing exact campaign plan; resume incomplete cases and verify/audit completed ones")
    args = parser.parse_args(argv)
    campaign_dir = scratch_path(args.output_dir)
    plan = build_plan(campaign_dir)
    plan_path = campaign_dir/"campaign_plan.json"
    if campaign_dir.exists():
        if not args.resume:
            raise FileExistsError("Preserve existing campaign; add --resume with identical plan")
        if not plan_path.is_file() or json.loads(plan_path.read_text()) != plan:
            raise ValueError("Campaign source/parameter plan changed; choose a fresh directory")
    else:
        if args.resume:
            raise ValueError("Existing campaign is required for --resume")
        campaign_dir.mkdir(parents=True)
        write_json(plan_path,plan)
    selected = selected_cases(args.group)
    print(json.dumps(dict(status="PLANNED",group=args.group,
        case_names=[case["name"] for case in selected],maximum_distinct_cases=16,
        output_dir=str(campaign_dir),long_runs_started=False)),flush=True)
    if args.plan_only:
        return plan
    with execution_lock(campaign_dir):
        return execute_campaign(plan,campaign_dir,selected,args)


def execute_campaign(plan,campaign_dir,selected,args):
    plan_path = campaign_dir/"campaign_plan.json"
    for name in ("runs","logs","audits"):
        (campaign_dir/name).mkdir(exist_ok=True)
    receipt_path = campaign_dir/"campaign_execution.json"
    receipt = json.loads(receipt_path.read_text()) if receipt_path.exists() else dict(
        schema="pysnspd.effective_sensitivity_execution.v1",
        plan_sha256=sha256(plan_path), cases={},model_admission_inferred=False)
    if receipt.get("plan_sha256") != sha256(plan_path):
        raise ValueError("Execution receipt belongs to a different plan")
    env = child_environment()
    try:
        for index,case in enumerate(selected,1):
            run_dir = campaign_dir/"runs"/case["name"]
            complete = False
            if run_dir.exists():
                manifest = check_case_identity(case,run_dir,plan["source_pins_sha256"],require_complete=False)
                complete = manifest.get("run_complete") is True
                if complete:
                    check_case_identity(case,run_dir,plan["source_pins_sha256"],require_complete=True)
                elif not args.resume:
                    raise ValueError("Incomplete case requires --resume")
            print(json.dumps(dict(status="CASE_AUDIT" if complete else "CASE_RUN",
                group=args.group,index=index,selected_count=len(selected),case=case["name"])),flush=True)
            if not complete:
                run_logged(runner_command(case,campaign_dir,resume=run_dir.exists()),
                    campaign_dir/"logs"/(case["name"]+".console.log"),env)
            manifest = check_case_identity(case,run_dir,plan["source_pins_sha256"],require_complete=True)
            audit = audit_case(case,campaign_dir,env)
            receipt["cases"][case["name"]] = dict(status="COMPLETE_AUDITED",
                run_dir=str(run_dir), manifest_sha256=sha256(run_dir/"manifest.json"),
                observation=manifest["observation"], audit=audit)
            receipt.update(last_group=args.group,status="GROUP_IN_PROGRESS")
            write_json(receipt_path,receipt)
        aggregate = aggregate_evidence(plan,receipt,campaign_dir,selected,args.group)
        aggregate_path = campaign_dir/("aggregate_"+args.group+".json")
        write_json(aggregate_path,aggregate)
        receipt.update(status="SELECTED_GROUP_COMPLETE_AUDITED",last_group=args.group,
            all16_complete=aggregate["all16_complete_audited"],
            aggregate=dict(path=str(aggregate_path),sha256=sha256(aggregate_path)),
            selected_group_complete=True)
        write_json(receipt_path,receipt)
        print(json.dumps(dict(status=receipt["status"],group=args.group,
            all16_complete=receipt["all16_complete"],receipt=str(receipt_path))),flush=True)
        return receipt
    except BaseException as error:
        receipt.update(status="STOPPED_WITH_EVIDENCE_PRESERVED",error=repr(error),
            last_group=args.group,selected_group_complete=False)
        write_json(receipt_path,receipt)
        raise


if __name__ == "__main__":
    main()
