"""Six simultaneous v1.2 cases, four SMT threads on two physical cores each.

Preparation never launches a simulation. Actual execution requires twelve
free worker cores plus a coordinator, two completely reserved cores, and the
same runtime lock as v1.1 campaigns. A requested CPU percentage is not a
performance guarantee.
"""
from __future__ import annotations
import argparse,fcntl,importlib.util,json,os,signal,subprocess,sys,time
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor,as_completed
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
for key in ("OPENBLAS_NUM_THREADS","OMP_NUM_THREADS","MKL_NUM_THREADS","NUMEXPR_NUM_THREADS"):
    os.environ[key]="1"
os.environ.setdefault("NUMBA_NUM_THREADS","4")
spec=importlib.util.spec_from_file_location("_v12_single",ROOT/"pipelines/13_dimensional_photon_run.py")
RUN=importlib.util.module_from_spec(spec);spec.loader.exec_module(RUN)
spec=importlib.util.spec_from_file_location("_v12_budget",ROOT/"pipelines/10_thermal_stiffness_basecases.py")
DRIVER=importlib.util.module_from_spec(spec);spec.loader.exec_module(DRIVER)
DRIVER.KNOWN_RUNNERS.update({"11_bias_escape_controls.py","12_low_bias_latency_controls.py",
                            "13_dimensional_photon_run.py","14_dimensional_basecases.py"})
PYTHON="/home/jdiaz/.conda/envs/snspd/bin/python"


def cases():
    return [dict(name="bias"+str(i).replace(".","p")+"_ee"+str(int(t))+"_"+str(w),
                 current_uA=i,tau_ee_Tc_ps=t,wavelength_nm=w)
            for i,t in ((15.5,5.),(15.5,4.),(21.5,4.)) for w in (775,1550)]


def allocation(resources=None,active=None):
    res=DRIVER.linux_resources() if resources is None else resources
    budget=DRIVER.resource_budget(res,.9,.9,1)
    active=DRIVER.active_known_processes() if active is None else active
    occupied={c for row in active for c in row["affinity_cpus"]}
    groups=res["physical_core_groups"];allowed=set(budget["selected_affinity_cpus"])
    free={i for i,g in enumerate(groups) if set(g)<=allowed and not set(g)&occupied}
    slots=[]
    for _ in range(6):
        first=next((i for i in sorted(free) if i+1 in free),None)
        if first is None:raise ValueError("Six cases need twelve FREE worker physical cores; active jobs are preserved")
        slots.append(dict(cpu_start=first,core_groups=groups[first:first+2],threads=4))
        free-={first,first+1}
    if not free:raise ValueError("A separate complete coordinator core is required")
    coordinator=groups[max(free)]
    assigned=set(coordinator)|{c for slot in slots for g in slot["core_groups"] for c in g}
    reserved=[g for g in groups if not set(g)&(assigned|occupied)]
    if len(reserved)<2:raise ValueError("Global budget must reserve two complete physical cores")
    if 27*1024**3>budget["memory_limit_bytes"]:raise ValueError("Insufficient RAM for six conservative worker reservations")
    return dict(inventory=res,slots=slots,coordinator=coordinator,reserved=reserved,
                worker_threads=24,worker_physical_cores=12,active_known_processes=active)


def command(case,out,slot,a,resume=False):
    cmd=[PYTHON,str(ROOT/"pipelines/13_dimensional_photon_run.py"),"--output-dir",str(out),
         "--threads","4","--physical-cores","2","--cpu-start",str(slot["cpu_start"]),
         "--current-uA",str(case["current_uA"]),"--tau-ee-Tc-ps",str(case["tau_ee_Tc_ps"]),
         "--wavelength-nm",str(case["wavelength_nm"]),"--length-nm","300",
         "--tau-escape-ps","20","--retained-fraction","0.8","--sheet-resistance-ohm","425",
         "--confirmation-ps","0.5","--extra-time-ps","2","--ceiling-ps",str(a.ceiling_ps),
         "--final-time-ns",str(a.final_time_ns)]
    if resume:cmd.append("--resume")
    return cmd


def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--campaign-dir",type=Path,required=True)
    p.add_argument("--prepare-only",action="store_true");p.add_argument("--resume",action="store_true")
    p.add_argument("--threads-per-case",type=int,choices=(4,),default=4)
    p.add_argument("--physical-cores-per-case",type=int,choices=(2,),default=2)
    p.add_argument("--ceiling-ps",type=float,default=100.);p.add_argument("--final-time-ns",type=float,default=15.)
    a=p.parse_args(argv);out=RUN.B.scratch_path(a.campaign_dir)
    if (out/"campaign_execution.json").exists() and not (a.resume or a.prepare_only):raise FileExistsError("Executed campaign exists: choose --resume")
    out.mkdir(parents=True,exist_ok=True)
    pins=RUN.source_pins();pins["pipelines/14_dimensional_basecases.py"]=RUN.B.sha256(__file__)
    pins["pipelines/10_thermal_stiffness_basecases.py"]=RUN.B.sha256(ROOT/"pipelines/10_thermal_stiffness_basecases.py")
    definition=dict(schema="pysnspd.dimensional_campaign.v1",model_version="1.2.0",cases=cases(),
        source_pins_sha256=pins,ceiling_ps=a.ceiling_ps,final_time_ns=a.final_time_ns,
        scheduling=dict(simultaneous_cases=6,threads_per_case=4,physical_cores_per_case=2,
            worker_threads=24,worker_physical_cores=12,SMT=True,minimum_reserved_physical_cores=2),
        scope="300nm same strip; guarded overlap; no amplifier; native kinetics/full circuit; no quantitative admission")
    plan=out/"campaign_plan.json"
    if plan.exists() and json.loads(plan.read_text())!=definition:raise ValueError("Existing campaign sources/parameters differ")
    RUN.write_json(plan,definition)
    if a.prepare_only:
        try:resources=allocation();status="PREPARED"
        except ValueError as exc:resources=dict(reason=str(exc));status="PREPARED_WAITING_CAPACITY"
        RUN.write_json(out/"preparation_resources.json",resources)
        print(json.dumps(dict(status=status,cases=6,no_simulations_started=True,resources=resources)),flush=True)
        return 0
    with (out/".execution.lock").open("a+") as local_lock, (DRIVER.SCRATCH/"thermal_candidate_runtime.lock").open("a+") as global_lock:
        try:
            fcntl.flock(local_lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
            fcntl.flock(global_lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        except BlockingIOError:raise RuntimeError("Another campaign is running; preserve it and retry later")
        resources=allocation();original_affinity=set(os.sched_getaffinity(0))
        os.sched_setaffinity(0,set(resources["coordinator"]))
        RUN.write_json(out/"execution_resources.json",resources)
        stop={"requested":False};processes={}
        def request(sig,frame):
            stop["requested"]=True
            for proc in list(processes.values()):
                if proc.poll() is None:proc.send_signal(signal.SIGTERM)
        signal.signal(signal.SIGTERM,request);signal.signal(signal.SIGINT,request)
        def worker(case,slot):
            dest=out/"runs"/case["name"];dest.parent.mkdir(exist_ok=True)
            if (dest/"manifest.json").exists():
                prior=json.loads((dest/"manifest.json").read_text())
                if prior.get("run_complete"):return dict(case=case["name"],status=prior["status"],reused_complete=True)
                if prior.get("status","").startswith("POSSIBLE_"):
                    return dict(case=case["name"],status=prior["status"],diagnostic_preserved=True)
            cmd=command(case,dest,slot,a,resume=a.resume and dest.exists())
            if stop["requested"]:return dict(case=case["name"],status="NOT_STARTED")
            env=dict(os.environ,NUMBA_NUM_THREADS="4",PYTHONPYCACHEPREFIX="/home/jdiaz/scratch/tmp/pysnspd/pycache",
                     MPLCONFIGDIR="/home/jdiaz/scratch/tmp/pysnspd/matplotlib",
                     NUMBA_CACHE_DIR="/home/jdiaz/scratch/tmp/pysnspd/numba_cache")
            with (out/(case["name"]+".log")).open("a") as log:
                # Child first inventories the full admitted topology, then pins its own SMT group.
                proc=subprocess.Popen(["taskset","-c",",".join(map(str,sorted(original_affinity)))]+cmd,
                    cwd=ROOT,env=env,stdout=log,stderr=subprocess.STDOUT)
                processes[case["name"]]=proc
                code=proc.wait()
            path=dest/"manifest.json";m=json.loads(path.read_text()) if path.exists() else {}
            return dict(case=case["name"],exit_code=code,status=m.get("status","FAILED_BEFORE_MANIFEST"),
                        time_ps=m.get("time_ps"),phase=m.get("phase"),run_complete=m.get("run_complete",False))
        results=[]
        with ThreadPoolExecutor(max_workers=6) as pool:
            jobs=[pool.submit(worker,c,s) for c,s in zip(cases(),resources["slots"])]
            for job in as_completed(jobs):
                result=job.result();results.append(result)
                RUN.write_json(out/"campaign_execution.json",dict(results=results,total=6,complete=len(results)))
                print(json.dumps(result),flush=True)
        good=all(x.get("run_complete") for x in results)
        print(json.dumps(dict(status="COMPLETE" if good else "COMPLETE_WITH_DIAGNOSTICS_OR_FAILURES",results=results)),flush=True)
        return 0 if good else 2


if __name__=="__main__":raise SystemExit(main())
