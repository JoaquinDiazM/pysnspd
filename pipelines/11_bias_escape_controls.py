"""Four bias/escape controls of authenticated photon anchors; no solver changes.

Reuses pipeline10's locking, CPU/RAM admission, queued workers, checkpoints and
per-case terminal audits. Only current, escape time, wavelength and the declared
observation ceiling may differ from the authenticated reference pair.
"""
from __future__ import annotations
import argparse
import copy
import importlib.util
import json
import math
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
SHARED=ROOT/"pipelines/10_thermal_stiffness_basecases.py"


def load_driver():
    spec=importlib.util.spec_from_file_location("_bias_escape_shared_driver",SHARED)
    driver=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(driver)
    return driver


def configure(driver,args):
    reference=driver.scratch_path(args.reference_campaign)
    definition=json.loads((reference/"campaign_plan.json").read_text())
    execution=json.loads((reference/"campaign_execution.json").read_text())
    driver.require_anchors(reference,definition,execution)
    manifests={}
    for wavelength in (775,1550):
        path=reference/"runs"/f"anchor425_eta0800_{wavelength}"/"manifest.json"
        m=json.loads(path.read_text())
        if (not m["run_complete"] or m["model_version"]!="1.1.3-candidate"
                or m["parameter_identity"]["stiffness_temperature"]!="electron"
                or m["parameter_identity"]["gradient_closure"]!="thermal-depairing"):
            raise ValueError("Reference must be the completed current thermal-depairing pair")
        if m["source_pins_sha256"]!=driver.pins("electron"):
            raise ValueError("Reference physics sources differ from current producer")
        manifests[wavelength]=(path,m)
    base=dict(manifests[775][1]["parameter_identity"])
    if base["threads"]!=args.threads_per_case:
        raise ValueError("Keep the reference thread count per case")
    for key,value in manifests[1550][1]["parameter_identity"].items():
        if key!="wavelength_nm" and value!=base[key]:
            raise ValueError("Reference pair has non-wavelength parameter differences: "+key)
    base.pop("threads")
    if len(set(args.bias_uA))!=2:
        raise ValueError("Two distinct bias currents are required")
    for current in args.bias_uA:
        if not math.isfinite(current) or not 0<current<base["target_depairing_uA"]:
            raise ValueError("Bias must be positive and below the unchanged depairing target")
    if not math.isfinite(args.tau_escape_ps) or args.tau_escape_ps<=0:
        raise ValueError("Positive finite escape time required")
    if not math.isfinite(args.ceiling_ps) or args.ceiling_ps<=0:
        raise ValueError("Positive finite observation ceiling required")
    cases=[]
    for current in args.bias_uA:
        label="bias"+format(current,"g").replace(".","p")+"_escape"+format(args.tau_escape_ps,"g").replace(".","p")
        for wavelength in (775,1550):
            params=dict(base,current_uA=current,wavelength_nm=float(wavelength),
                        tau_escape_ps=args.tau_escape_ps,ceiling_ps=args.ceiling_ps)
            driver.RUNNER.validate_args(argparse.Namespace(**params,threads=args.threads_per_case,
                cpu_start=0,output_dir=driver.SCRATCH/"unused_bias_escape_validation",resume=False))
            cases.append(dict(name=f"{label}_{wavelength}",pair=label,stage="anchors",parameters=params))
    old_plan=driver.plan
    old_allocation=driver.allocation
    driver.catalog=lambda *,ceiling_ps=100.:copy.deepcopy(cases)

    def plan(campaign,*,ceiling_ps=100.):
        result=old_plan(campaign,ceiling_ps=ceiling_ps)
        result.update(schema="pysnspd.bias_escape_controls.v1",orchestrator_sha256=driver.sha(__file__),
            shared_campaign_driver_sha256=driver.sha(SHARED),stage_order=["anchors"],
            reference_manifests={str(w):dict(path=str(path),sha256=driver.sha(path))
                                 for w,(path,m)in manifests.items()})
        result["scope"].pop("current_uA",None)
        result["scope"].update(bias_currents_uA=args.bias_uA,tau_escape_ps=args.tau_escape_ps,
            reference_campaign=str(reference),same_reference_material_and_calibrators=True,
            no_per_bias_diffusion_retuning=True)
        result["scheduling"].update(anchors="Four cases queued in two slots; exactly six threads per case",
            commands_are_argument_templates_affinity_assigned_at_runtime=True)
        result["scheduling"].pop("controls",None)
        result["gate"]=dict(each_completed_case_must_pass_original_audit=True,required_cases=4,
            no_quantitative_admission_inferred=True)
        return result

    def allocation(stage,*,anchor_threads=None,**kwargs):
        if stage!="anchors" or anchor_threads not in (None,args.threads_per_case):
            raise ValueError("Unexpected stage or changed resume thread identity")
        return old_allocation(stage,anchor_threads=args.threads_per_case,**kwargs)

    driver.plan=plan
    driver.allocation=allocation
    return cases


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference-campaign",type=Path,required=True)
    parser.add_argument("--campaign-dir",type=Path,required=True)
    parser.add_argument("--bias-uA",type=float,nargs=2,default=[15.5,21.5])
    parser.add_argument("--tau-escape-ps",type=float,default=15.)
    parser.add_argument("--threads-per-case",type=int,choices=(6,),default=6)
    parser.add_argument("--max-worker-cores",type=int,default=12,choices=(12,))
    parser.add_argument("--ceiling-ps",type=float,default=100.)
    parser.add_argument("--prepare-only",action="store_true")
    parser.add_argument("--resume",action="store_true")
    args=parser.parse_args(argv)
    driver=load_driver()
    cases=configure(driver,args)
    campaign=driver.scratch_path(args.campaign_dir)
    if args.prepare_only:
        prepared=driver.plan(campaign,ceiling_ps=args.ceiling_ps)
        resources=driver.allocation("anchors",max_worker_cores=args.max_worker_cores)
        campaign.mkdir(parents=True,exist_ok=True)
        path=campaign/"campaign_plan.json"
        if path.exists() and json.loads(path.read_text())!=prepared:
            raise ValueError("Existing plan differs; preserve it and choose a fresh directory")
        if not path.exists():driver.write_json(path,prepared)
        driver.write_json(campaign/"preparation_resources.json",
            dict(scope="Runtime availability at preparation only; not an active reservation",resources=resources))
        print(json.dumps(dict(status="PREPARED",cases=len(cases),case_names=[c["name"]for c in cases],
            simultaneous_cases=len(resources["slots"]),threads_per_case=args.threads_per_case,
            campaign_dir=str(campaign),no_simulation_executed=True)),flush=True)
        return 0
    words=["--campaign-dir",str(campaign),"--stage","anchors",
           "--max-worker-cores",str(args.max_worker_cores),"--ceiling-ps",str(args.ceiling_ps)]
    if args.resume:words.append("--resume")
    return driver.main(words)


if __name__=="__main__":
    raise SystemExit(main())
