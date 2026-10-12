"""Compare saved sensitivity runs in physical units, without evolving a model.

Electrical latency is taken only from confirmed accepted-step events. Maxima
are labelled as finite-window observed maxima; no full pulse peak is inferred.
"""
from __future__ import annotations
import argparse
import csv
import importlib.util
import json
import os
from pathlib import Path
import sys
import zipfile

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from sandbox.stage4_core.parallel_runtime import limit_thread_environment
limit_thread_environment()
os.environ.setdefault("MPLCONFIGDIR","/home/jdiaz/scratch/tmp/pysnspd/matplotlib")
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pysnspd.plotting.style import apply_thesis_style
spec = importlib.util.spec_from_file_location("effective_sensitivity_campaign",
    ROOT/"pipelines/06_effective_sensitivity_campaign.py")
campaign = importlib.util.module_from_spec(spec)
spec.loader.exec_module(campaign)


def uncompressed_trace_bytes(path):
    """NPZ file bytes are not a bound on array memory; inspect NPY ZIP sizes."""
    names = {"time_s.npy","Vout_V.npy","Vdev_V.npy","Is_A.npy"}
    with zipfile.ZipFile(path) as archive:
        sizes = {item.filename:item.file_size for item in archive.infolist()}
    if not names <= set(sizes):
        raise ValueError("Required saved trace arrays are missing")
    return sum(sizes[name] for name in names)


def event_metrics(manifest, time_ps, voltage_mV):
    observation = manifest["observation"]
    confirmed = observation.get("trigger_confirmed") is True
    cross = observation.get("crossing_time_s")
    latency = observation.get("latency_s")
    result = dict(status=observation.get("status"),confirmed=confirmed,
        censored=observation.get("censored") is True,crossing_ps=None,
        relative_latency_ps=None,local_slope_mV_per_ps=None,
        slope_fit_window_ps=None,slope_fit_sample_count=0)
    if not confirmed:
        if manifest.get("status") not in campaign.TERMINAL_STATUSES:
            raise ValueError("Incomplete states cannot be presented as censored trials")
        result["censoring_ps"] = float(observation["no_trigger_ceiling_s"])*1e12
        return result
    if cross is None or latency is None or not np.isfinite([cross,latency]).all():
        raise ValueError("Confirmed event requires finite saved crossing and latency")
    cross_ps = float(cross)*1e12
    if not time_ps[0] <= cross_ps <= time_ps[-1]:
        raise ValueError("Saved crossing is outside saved trajectory")
    mask = (time_ps >= cross_ps-.25)&(time_ps <= cross_ps+.25)
    if np.count_nonzero(mask) >= 3:
        times = time_ps[mask]
        slope = np.polyfit(times-cross_ps,voltage_mV[mask],1)[0]
        result.update(local_slope_mV_per_ps=float(slope),
            slope_fit_window_ps=[float(times[0]),float(times[-1])],
            slope_fit_sample_count=int(mask.sum()))
    result.update(crossing_ps=cross_ps,relative_latency_ps=float(latency)*1e12)
    return result


def load_case(case, campaign_dir, plan, receipt):
    name = case["name"]
    item = receipt["cases"].get(name)
    if not item or item.get("status") != "COMPLETE_AUDITED":
        raise ValueError("Complete audited case required: "+name)
    run_dir = campaign_dir/"runs"/name
    manifest = campaign.check_case_identity(case,run_dir,plan["source_pins_sha256"],require_complete=True)
    if campaign.sha256(run_dir/"manifest.json") != item["manifest_sha256"]:
        raise ValueError("Campaign receipt no longer authenticates terminal "+name)
    audit_path = campaign.scratch_path(item["audit"]["path"])
    if campaign.sha256(audit_path) != item["audit"]["sha256"]:
        raise ValueError("Audit changed for "+name)
    audit = json.loads(audit_path.read_text())
    if not audit.get("all_requested_checks_passed"):
        raise ValueError("Failed audit cannot be hidden in comparison")
    with np.load(run_dir/"trajectory.npz",allow_pickle=False) as data:
        time_ps = np.asarray(data["time_s"],float)*1e12
        voltage = np.asarray(data["Vout_V"],float)*1e3
        device = np.asarray(data["Vdev_V"],float)*1e3
        current = np.asarray(data["Is_A"],float)*1e6
    if time_ps.ndim != 1 or not len(time_ps) or np.any(np.diff(time_ps)<=0):
        raise ValueError("Increasing accepted trajectory required")
    if any(v.shape != time_ps.shape or not np.isfinite(v).all()
            for v in (time_ps,voltage,device,current)):
        raise ValueError("Finite coherent physical traces required")
    baseline = float(manifest["observation"]["baseline_Vout_V"])*1e3
    voltage -= baseline
    metric = dict(name=name,label=case["label"],group=case["group"],
        changed_input=case["changed_input"],reference=case["reference"],
        saved_start_ps=float(time_ps[0]),saved_end_ps=float(time_ps[-1]),
        observed_max_Vout_mV=float(voltage.max()),
        maximum_scope="Observed saved finite interval; not a full pulse peak",
        minimum_Is_uA=float(current.min()),
        **event_metrics(manifest,time_ps,voltage))
    return dict(case=case,manifest=manifest,time_ps=time_ps,voltage_mV=voltage,
        Vdev_mV=device,Is_uA=current,metrics=metric,
        pins=dict(manifest_sha256=campaign.sha256(run_dir/"manifest.json"),
            trajectory_sha256=campaign.sha256(run_dir/"trajectory.npz"),
            audit_sha256=item["audit"]["sha256"]))


def common_metrics(reference,candidate):
    start = max(reference["time_ps"][0],candidate["time_ps"][0])
    end = min(reference["time_ps"][-1],candidate["time_ps"][-1])
    if end <= start:
        raise ValueError("Nonzero common physical window required")
    times = np.unique(np.r_[start,end,
        reference["time_ps"][(reference["time_ps"]>start)&(reference["time_ps"]<end)],
        candidate["time_ps"][(candidate["time_ps"]>start)&(candidate["time_ps"]<end)]])
    a = np.interp(times,reference["time_ps"],reference["voltage_mV"])
    b = np.interp(times,candidate["time_ps"],candidate["voltage_mV"])
    detected = reference["metrics"]["confirmed"] and candidate["metrics"]["confirmed"]
    delta_latency = candidate["metrics"]["relative_latency_ps"]-reference["metrics"]["relative_latency_ps"] if detected else None
    return dict(common_window_ps=[float(start),float(end)],
        reference_observed_common_max_mV=float(a.max()),
        candidate_observed_common_max_mV=float(b.max()),
        maximum_absolute_Vout_difference_mV=float(np.max(abs(b-a))),
        delta_relative_latency_ps=delta_latency,
        latency_comparison="Saved confirmed events only; null if either trajectory is censored",
        curve_interpolation="Linear accepted-trace interpolation restricted to this common interval; does not create trigger events")


def families(runs):
    result = []
    choices = [
        ("Sesgo, 775 nm","current_uA",775.,"nominal775"),
        ("Sesgo, 1550 nm","current_uA",1550.,"nominal1550"),
        ("Calibrador tau_ee, 775 nm","tau_ee_Tc_ps",775.,"nominal775"),
        ("Calibrador tau_ep, 775 nm","tau_ep_Tc_ps",775.,"nominal775"),
        ("Escape, 775 nm","tau_escape_ps",775.,"nominal775"),
        ("Retencion, 775 nm","retained_fraction",775.,"nominal775"),
        ("Anchura postcascada, 775 nm","sigma_nm",775.,"nominal775"),
        ("Retencion, 1550 nm","retained_fraction",1550.,"nominal1550"),
    ]
    for title,key,wavelength,reference in choices:
        selected = [run for run in runs.values()
            if key in run["case"]["changed_input"] and run["case"]["parameters"]["wavelength_nm"]==wavelength]
        if selected:
            result.append((title,[runs[reference]]+selected))
    if not result:
        for name in ("nominal775","nominal1550"):
            if name in runs:
                result.append((runs[name]["case"]["label"],[runs[name]]))
    return result


def save_figure(fig,output,name,products,caption):
    paths = {}
    for suffix in ("png","pdf"):
        path = output/(name+"."+suffix)
        fig.savefig(path,dpi=170,facecolor="white",bbox_inches="tight")
        paths[suffix] = dict(path=str(path),sha256=campaign.sha256(path))
    products.append(dict(name=name,caption=caption,files=paths))
    plt.close(fig)


def curve_figures(runs,output,products):
    panels = families(runs)
    rows = (len(panels)+1)//2
    fig,axes = plt.subplots(rows,2,figsize=(12,3.65*rows),squeeze=False,layout="constrained")
    for ax,(title,members) in zip(axes.ravel(),panels):
        start = max(run["time_ps"][0] for run in members)
        end = min(run["time_ps"][-1] for run in members)
        for run in members:
            times = run["time_ps"]
            mask = (times>=start)&(times<=end)
            t = np.r_[start,times[mask],end]
            v = np.interp(t,times,run["voltage_mV"])
            line, = ax.plot(t,v,label=run["case"]["label"])
            cross = run["metrics"]["crossing_ps"]
            if cross is not None and start<=cross<=end:
                ax.axvline(cross,color=line.get_color(),ls=":",alpha=.7,lw=.8)
        ax.axhline(.1,color=".3",ls="--",lw=.8,label="Umbral 0.1 mV")
        ax.set(title=title+"\nVentana comun "+f"{start:.3g}–{end:.4g} ps",
            xlabel="Tiempo postcascada (ps)",ylabel="Vout menos basal (mV)")
        ax.legend(fontsize=7,loc="upper left")
        if max(run["metrics"]["observed_max_Vout_mV"] for run in members)<.02:
            inset = ax.inset_axes([.46,.22,.49,.35])
            for run in members:
                mask = (run["time_ps"]>=start)&(run["time_ps"]<=end)
                inset.plot(run["time_ps"][mask],run["voltage_mV"][mask]*1e3)
            inset.set(xlabel="ps",ylabel="µV")
            inset.tick_params(labelsize=7)
    for ax in axes.ravel()[len(panels):]:
        ax.set_visible(False)
    fig.suptitle("Sensibilidad electrica del SNSPD efectivo\nCada panel termina en el minimo de sus ventanas guardadas",fontsize=13)
    save_figure(fig,output,"01_Vout_common_windows",products,
        "Vout=50ohm*(Ib-Is), basal DC sustraido. Curvas solo en su ventana comun; lineas verticales solo cruces confirmados reales dentro de ella. El recuadro usa µV cuando las curvas son pequenas. No picos completos, ancho total, jitter ni reset inferidos.")


def event_figures(runs,output,products):
    values = list(runs.values())
    fig,axes = plt.subplots(2,1,figsize=(12,7),sharex=True,layout="constrained")
    x = np.arange(len(values))
    for i,run in enumerate(values):
        metric = run["metrics"]
        if metric["confirmed"]:
            axes[0].plot(i,metric["relative_latency_ps"],"o",color="C0")
            if metric["local_slope_mV_per_ps"] is not None:
                axes[1].plot(i,metric["local_slope_mV_per_ps"],"o",color="C1")
        else:
            for ax in axes:
                ax.text(i,.02,"cens.\n"+f'{metric["censoring_ps"]:g} ps',
                    transform=ax.get_xaxis_transform(),ha="center",va="bottom",fontsize=6)
    axes[0].set(ylabel="Latencia relativa (ps)",
        title="Eventos confirmados; los censurados no reciben una latencia ficticia")
    axes[1].set(ylabel="Pendiente local (mV/ps)",
        title="Ajuste lineal de muestras aceptadas en cruce ±0.25 ps; sin filtro del setup")
    axes[1].set_xticks(x,[run["case"]["name"] for run in values],rotation=60,ha="right",fontsize=7)
    axes[1].set_xlim(-.6,len(values)-.4)
    save_figure(fig,output,"02_confirmed_latency_and_local_slope",products,
        "Tiempo relativo a transferencia postcascada. Sólo controlador aceptado asigna cruce; censura a200ps se etiqueta sin punto numerico de latencia. Pendiente observada local con ancho de ajuste±0.25ps, no respuesta filtrada por ancho de banda experimental ni estimacion de jitter.")


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--campaign-dir",required=True,type=Path)
    parser.add_argument("--output-dir",required=True,type=Path)
    parser.add_argument("--group",choices=("bias","rates","source","all"),default="bias")
    parser.add_argument("--individual",action="store_true",
        help="Also call visual05 for each selected trajectory's maps (saved data only)")
    args = parser.parse_args(argv)
    directory,output = campaign.scratch_path(args.campaign_dir),campaign.scratch_path(args.output_dir)
    if output.exists():
        raise FileExistsError("Preserve plots; choose a fresh --output-dir")
    from sandbox.stage4_core.parallel_runtime import linux_resources,resource_budget
    budget = resource_budget(linux_resources(),.9,.9,1)
    plan = json.loads((directory/"campaign_plan.json").read_text())
    receipt = json.loads((directory/"campaign_execution.json").read_text())
    if receipt["plan_sha256"] != campaign.sha256(directory/"campaign_plan.json"):
        raise ValueError("Campaign plan authentication failed")
    catalog = {case["name"]:case for case in plan["cases"]}
    wanted = [case["name"] for case in campaign.selected_cases(args.group)]
    if any(name not in catalog for name in wanted):
        raise ValueError("Selected catalog absent from the recorded plan")
    numeric_size = sum(uncompressed_trace_bytes(directory/"runs"/name/"trajectory.npz") for name in wanted)
    if 6*numeric_size+512*1024**2 > budget["memory_limit_bytes"]:
        raise MemoryError("Conservative comparison input budget exceeds90% available memory")
    old_affinity = os.sched_getaffinity(0)
    os.sched_setaffinity(0,{budget["coordinator_cpu"]})
    try:
        runs = {name:load_case(catalog[name],directory,plan,receipt) for name in wanted}
        comparisons = {}
        for name,run in runs.items():
            if name != run["case"]["reference"]:
                comparisons[name] = common_metrics(runs[run["case"]["reference"]],run)
        output.mkdir(parents=True)
        apply_thesis_style()
        products = []
        curve_figures(runs,output,products)
        event_figures(runs,output,products)
        metrics = {name:run["metrics"] for name,run in runs.items()}
        columns = ["name","saved_end_ps","observed_max_Vout_mV","confirmed","censored",
            "relative_latency_ps","local_slope_mV_per_ps","minimum_Is_uA"]
        with (output/"metrics.csv").open("w",newline="") as stream:
            writer = csv.DictWriter(stream,fieldnames=columns,extrasaction="ignore")
            writer.writeheader()
            writer.writerows(metrics.values())
        result = dict(schema="pysnspd.effective_sensitivity_plots.v1",
            status="SELECTED_AUDITED_SAVED_DATA_COMPARED",group=args.group,
            all_selected_individual_audits_passed=True,
            plan_sha256=campaign.sha256(directory/"campaign_plan.json"),
            execution_receipt_sha256=campaign.sha256(directory/"campaign_execution.json"),
            analysis_source_sha256=campaign.sha256(__file__),
            run_inputs={name:run["pins"] for name,run in runs.items()},
            metrics=metrics,common_window_comparisons=comparisons,products=products,
            resources=dict(processes=1,selected_cpu=budget["coordinator_cpu"],
                memory_limit_bytes=budget["memory_limit_bytes"]),
            scope="Saved physical sensitivity, not automatic convergence, experimental calibration, absolute optical latency, jitter, peak pulse, or total reset")
        campaign.write_json(output/"comparison_manifest.json",result)
    finally:
        os.sched_setaffinity(0,old_affinity)
    if args.individual:
        # Visual05 needs the broad inventory first, then pins its renderer itself.
        # Pre-pinning this child to one CPU would invalidate its own resource gate.
        for name in wanted:
            command = [str(campaign.PYTHON),str(ROOT/"plot_pipelines/05_plot_effective_photon_run.py"),
                "--run-dir",str(directory/"runs"/name),
                "--output-dir",str(output/"individual"/name)]
            env = campaign.child_environment()
            env["NUMBA_NUM_THREADS"] = "1"
            campaign.run_logged(command,output/"individual_plots.log",env)
    print(json.dumps(dict(status=result["status"],cases=len(runs),output_dir=str(output))))
    return result


if __name__ == "__main__":
    main()
