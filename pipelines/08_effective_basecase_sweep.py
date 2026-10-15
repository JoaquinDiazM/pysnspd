"""Paired 775/1550 nm base-case sweep, with bounded foreground concurrency.

All cases use the unchanged effective runner and accepted-step trigger.
No statistical jitter, readout filtering, PDF generation or full reset.
"""
from __future__ import annotations
import argparse, csv, datetime as dt, fcntl, hashlib, importlib.util, json
import math, os, queue, signal, subprocess, sys, threading, time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
SPEC = importlib.util.spec_from_file_location("basecase_existing_campaign", ROOT/"pipelines/06_effective_sensitivity_campaign.py")
BASE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(BASE)
from sandbox.stage4_core.parallel_runtime import linux_resources, resource_budget

ETAS = (.667, .7, .8, .9, 1.)
COMMON = dict(BASE.COMMON, current_uA=21.5, tau_escape_ps=20., threads=1)
TERMINAL = BASE.TERMINAL_STATUSES


def catalog():
    result = []
    def pair(name, group, overrides):
        for wavelength in (775., 1550.):
            result.append(dict(name=f"{name}_{int(wavelength)}", pair=name, group=group,
                parameters=dict(COMMON, **dict(overrides, wavelength_nm=wavelength))))
    for eta in ETAS:
        pair("eta"+str(int(round(eta*1000))).zfill(4), "retention", dict(retained_fraction=eta))
    for name, value in (("ee_half", 2.5), ("ee_double", 10.)):
        pair(name, "mobility_and_BGK", dict(retained_fraction=.8, tau_ee_Tc_ps=value))
    for name, factor in (("sigma_narrow", 2**-.5), ("sigma_wide", 2**.5)):
        pair(name, "preparation", dict(retained_fraction=.8, sigma_nm=BASE.SIGMA_REFERENCE_NM*factor))
    pair("length640", "domain", dict(retained_fraction=.8, length_nm=640.))
    return result


def command(case, campaign, slot, resume=False):
    result = BASE.runner_command(case, campaign, resume=resume)
    result[result.index("--cpu-start")+1] = str(slot)
    return result


def inventory(workers):
    resources = linux_resources()
    budget = resource_budget(resources, .9, .9, 1)
    allowed = set(budget["selected_affinity_cpus"])
    cores = [g for g in resources["physical_core_groups"] if set(g) <= allowed]
    # One independent core per run, one complete admitted core for coordination.
    memory_slots = (budget["memory_limit_bytes"]-3*1024**3)//(4*1024**3)
    maximum = min(len(cores)-1, int(memory_slots))
    minimum = math.ceil(.5*len(resources["physical_core_groups"]))
    if not minimum <= workers <= maximum:
        raise ValueError(f"Require {minimum}..{maximum} workers for >=50% physical-core allocation and <=90% CPU/RAM budget")
    return dict(inventory=resources, budget=budget, simultaneous_runs=workers,
        physical_cores_per_run=1, worker_core_groups=cores[:workers],
        coordinator_core_group=cores[-1], physical_worker_fraction=workers/len(resources["physical_core_groups"]),
        reserved_outside_budget=[g for g in resources["physical_core_groups"] if not set(g)<=allowed],
        memory_reservation_bytes=(3+4*workers)*1024**3,
        reservation_scope="Scheduling estimate, not an OS memory limit or a CPU utilization guarantee",
        logical_monitor_note="One hardware thread per physical core; SMT siblings are not counted as additional physical capacity")


def plan(campaign):
    cases = catalog()
    return dict(schema="pysnspd.paired_basecase_sweep.v1", campaign_dir=str(campaign),
        source_pins_sha256=BASE.source_pins(), orchestrator_sha256=BASE.sha256(__file__),
        helper_sha256=BASE.sha256(BASE.__file__), auditor_sha256=BASE.sha256(BASE.AUDITOR),
        cases=cases, commands={c["name"]:command(c,campaign,0) for c in cases},
        scope=dict(current_uA=21.5, tau_escape_ps=20., sheet_resistance_ohm=608.,
            eta_definition="Fraction retained in electron+phonon excitations at post-cascade transfer; conditional on absorption",
            eta_one="Upper-bound initial retention; later 20 ps phonon escape remains active",
            eta_points=list(ETAS), paired_wavelengths_nm=[775.,1550.],
            controls_reference="eta0800, only the named input changes",
            width_control="Common FD/BE initialization changes spatial profile AND initial electron/phonon partition",
            tau_ee_control="BGK and KWT mobility change together",
            domain_control="Same 5 um detector; exterior inductance subtracts the newly resolved length once",
            observation="Positive Vout0.1mV; confirmation0.5ps; stop at crossing+5ps or censor at200ps",
            no_full_reset=True, no_statistical_jitter=True, no_PDF=True,
            no_physical_admission_from_execution_or_balance=True))


def assert_no_existing_runner():
    uid=os.getuid()
    for directory in Path("/proc").iterdir():
        if not directory.name.isdigit() or int(directory.name)==os.getpid():
            continue
        try:
            if directory.stat().st_uid != uid:
                continue
            words=(directory/"cmdline").read_bytes().split(b"\0")
        except (FileNotFoundError, PermissionError, ProcessLookupError):
            continue
        if any(Path(os.fsdecode(word)).name=="05_effective_photon_run.py" for word in words if word):
            raise RuntimeError("An effective photon run is already active; global resource budget must be shared explicitly")


def write_pairs(campaign, definition, receipt):
    pairs={}
    for case in definition["cases"]:
        item=receipt["cases"].get(case["name"],{})
        if item.get("status") != "COMPLETE_AUDITED":
            continue
        path=campaign/"runs"/case["name"]
        manifest=BASE.check_case_identity(case,path,definition["source_pins_sha256"],require_complete=True)
        obs=manifest["observation"]
        crossing=obs.get("first_crossing_time_s")
        if not obs.get("trigger_confirmed"):
            crossing=None
        pairs.setdefault(case["pair"],{})[int(case["parameters"]["wavelength_nm"])]=(case,manifest,crossing)
    target=campaign/"latency_pairs.csv"
    temporary=target.with_suffix(".csv.part")
    with temporary.open("w",newline="") as stream:
        writer=csv.writer(stream)
        writer.writerow(["pair","eta","tau_ee_Tc_ps","sigma_nm","length_nm","status775","status1550",
            "crossing775_ps","crossing1550_ps","difference1550_minus775_ps"])
        for name,data in pairs.items():
            if set(data)!={775,1550}:
                continue
            c,m,t=data[775];_,m2,t2=data[1550];p=c["parameters"]
            writer.writerow([name,p["retained_fraction"],p["tau_ee_Tc_ps"],p["sigma_nm"],p["length_nm"],
                m["status"],m2["status"],"" if t is None else t*1e12,
                "" if t2 is None else t2*1e12,"" if t is None or t2 is None else (t2-t)*1e12])
    os.replace(temporary,target)


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--campaign-dir",type=Path,required=True)
    parser.add_argument("--workers",type=int,default=13)
    parser.add_argument("--prepare-only",action="store_true")
    parser.add_argument("--resume",action="store_true")
    args=parser.parse_args(argv)
    campaign=BASE.scratch_path(args.campaign_dir)
    resources=inventory(args.workers)
    definition=plan(campaign)
    campaign.mkdir(parents=True,exist_ok=True)
    with (campaign/".execution.lock").open("a") as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        planpath=campaign/"campaign_plan.json"
        if planpath.exists():
            if json.loads(planpath.read_text()) != definition:
                raise ValueError("Saved plan/source identity differs; preserve this directory and choose a new one")
        else:
            if (campaign/"runs").exists():
                raise ValueError("Run data exists without an authenticated plan")
            BASE.write_json(planpath,definition)
        print(json.dumps(dict(status="PREPARED" if args.prepare_only else "RESOURCE_PLAN",
            cases=len(definition["cases"]),campaign_dir=str(campaign),resources=resources)),flush=True)
        if args.prepare_only:
            return 0
        assert_no_existing_runner()
        for leaf in ("logs","audits","runs"):
            (campaign/leaf).mkdir(exist_ok=True)
        path=campaign/"campaign_execution.json"
        receipt=json.loads(path.read_text()) if path.exists() else dict(cases={})
        if receipt.get("plan_sha256",BASE.sha256(planpath)) != BASE.sha256(planpath):
            raise ValueError("Execution receipt belongs to a different plan")
        if path.exists() and not args.resume:
            raise ValueError("Existing execution requires --resume")
        receipt.update(plan_sha256=BASE.sha256(planpath),resources=resources,status="RUNNING")
        BASE.write_json(path,receipt)
        env=BASE.child_environment()
        env["NUMBA_NUM_THREADS"]="1"
        for key in ("BLIS_NUM_THREADS","VECLIB_MAXIMUM_THREADS","NUMEXPR_MAX_THREADS"):
            env[key]="1"
        slots=queue.Queue()
        for slot in range(args.workers):
            slots.put(slot)
        original_affinity=set(os.sched_getaffinity(0))
        mutex=threading.RLock()
        stop=threading.Event()
        children={}
        stamp=dt.datetime.now(dt.timezone.utc).strftime("%Y%m%dT%H%M%SZ")
        def emit(value):
            with mutex:
                print(json.dumps(value),flush=True)
        def update(name,value):
            with mutex:
                receipt["cases"][name]=value
                BASE.write_json(path,receipt)
        def run_logged(cmd, logfile, name):
            if stop.is_set():
                raise InterruptedError("Campaign interruption requested")
            # Restore the full inherited CPU set before the runner rechecks topology.
            # The coordinator itself remains on its separate admitted core.
            cmd=["taskset","-c",",".join(map(str,sorted(original_affinity))),*cmd]
            with logfile.open("a") as stream:
                stream.write("\nCOMMAND "+__import__("shlex").join(cmd)+"\n");stream.flush()
                with subprocess.Popen(cmd,cwd=ROOT,env=env,stdout=subprocess.PIPE,
                        stderr=subprocess.STDOUT,text=True,bufsize=1) as child:
                    with mutex:
                        children[child.pid]=child
                        stopping=stop.is_set()
                    if stopping:child.terminate()
                    try:
                        for line in child.stdout:
                            stream.write(line);stream.flush()
                            emit(dict(case=name,output=line.rstrip()))
                        code=child.wait()
                    finally:
                        with mutex:children.pop(child.pid,None)
                    if code:raise subprocess.CalledProcessError(code,cmd)
        def work(case):
            slot=slots.get()
            name=case["name"];run_dir=campaign/"runs"/name
            try:
                if stop.is_set():raise InterruptedError("Campaign interruption requested")
                complete=False
                if run_dir.exists():
                    m=BASE.check_case_identity(case,run_dir,definition["source_pins_sha256"],require_complete=False)
                    complete=m.get("run_complete") is True and m.get("status") in TERMINAL
                    if complete:
                        BASE.check_case_identity(case,run_dir,definition["source_pins_sha256"],require_complete=True)
                    elif not args.resume:
                        raise ValueError("Incomplete case requires --resume")
                update(name,dict(status="AUDITING" if complete else "RUNNING",slot=slot))
                if not complete:
                    run_logged(command(case,campaign,slot,resume=run_dir.exists()),
                        campaign/"logs"/f"{name}_{stamp}.console.log",name)
                manifest=BASE.check_case_identity(case,run_dir,definition["source_pins_sha256"],require_complete=True)
                audit_dir=campaign/"audits"/f"{name}_{stamp}"
                audit_cmd=["taskset","-c",str(resources["worker_core_groups"][slot][0]),
                    str(BASE.PYTHON),str(BASE.AUDITOR),"--run",name+"="+str(run_dir),
                    "--output-dir",str(audit_dir)]
                run_logged(audit_cmd,campaign/"logs"/f"{name}_{stamp}.audit.log",name)
                audit=json.loads((audit_dir/"audit.json").read_text())
                if audit.get("all_requested_checks_passed") is not True:
                    raise ValueError("Saved-data audit failed: "+str(audit_dir))
                value=dict(status="COMPLETE_AUDITED",run_dir=str(run_dir),
                    manifest_sha256=BASE.sha256(run_dir/"manifest.json"),
                    observation=manifest["observation"],audit_path=str(audit_dir/"audit.json"),
                    audit_sha256=BASE.sha256(audit_dir/"audit.json"))
                update(name,value)
                emit(dict(status="CASE_COMPLETE_AUDITED",case=name))
                return True
            except Exception as error:
                update(name,dict(status="INTERRUPTED" if stop.is_set() else "FAILED",error=repr(error)))
                emit(dict(status="CASE_FAILED",case=name,error=repr(error)))
                return False
            finally:
                slots.put(slot)
        previous={}
        def interrupt(signum,frame):
            stop.set()
            with mutex:
                running=list(children.values())
            for child in running:
                if child.poll() is None:
                    try:child.terminate()
                    except ProcessLookupError:pass
        for signum in (signal.SIGINT,signal.SIGTERM):
            previous[signum]=signal.signal(signum,interrupt)
        try:
            coordinator=resources.get("coordinator_core_group",[min(original_affinity)])[0]
            os.sched_setaffinity(0,{coordinator})
            with ThreadPoolExecutor(max_workers=args.workers) as executor:
                futures=[executor.submit(work,c) for c in definition["cases"]]
                outcomes=[future.result() for future in as_completed(futures)]
            write_pairs(campaign,definition,receipt)
            receipt.update(status="INTERRUPTED" if stop.is_set() else
                ("COMPLETE_AUDITED" if all(outcomes) else "COMPLETE_WITH_FAILURES"),
                completed_cases=sum(item.get("status")=="COMPLETE_AUDITED" for item in receipt["cases"].values()))
            BASE.write_json(path,receipt)
            emit(dict(status=receipt["status"],completed_cases=receipt["completed_cases"],
                total_cases=len(definition["cases"]),pairs_csv=str(campaign/"latency_pairs.csv")))
            return 0 if all(outcomes) and not stop.is_set() else 1
        finally:
            os.sched_setaffinity(0,original_affinity)
            for signum,handler in previous.items():signal.signal(signum,handler)


if __name__=="__main__":
    raise SystemExit(main())
