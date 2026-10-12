"""Physical scientific figures from one effective neutral photon trajectory.

Consumes saved data only. Maps retain physical units and common color limits
across time. This visual pipeline does not infer model admission, fabricate a
trigger, interpolate time snapshots or launch a simulation.
"""
from __future__ import annotations

import argparse
import hashlib
import inspect
import json
import os
from pathlib import Path
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from sandbox.stage4_core.parallel_runtime import limit_thread_environment
limit_thread_environment()
import numpy as np
os.environ.setdefault("MPLCONFIGDIR", "/home/jdiaz/scratch/tmp/pysnspd/matplotlib")
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.tri import Triangulation
from pysnspd.plotting.style import apply_thesis_style, THESIS_DPI
from pysnspd.plotting.complex_gap import draw_complex_gap
COMPLEX_GAP_SOURCE = Path(inspect.getfile(draw_complex_gap)).resolve()

SCRATCH = Path("/home/jdiaz/scratch")
MEV_J = 1.602176634e-22
EV_J = 1.602176634e-19
SCHEMA = "pysnspd.effective_photon_plots.v1"


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024**2), b""):
            digest.update(block)
    return digest.hexdigest()


def in_scratch(path):
    result = Path(path).expanduser().resolve()
    if not result.is_relative_to(SCRATCH.resolve()):
        raise ValueError("Numeric products and renders must remain under /home/jdiaz/scratch")
    return result


def time_label_ps(time_s):
    value = float(time_s)*1e12
    if not np.isfinite(value):
        raise ValueError("A finite physical time is required")
    return f"{value:.5g} ps"


def nearest_snapshots(time_s, requested_ps=None):
    """Nearest saved states, always including initial and final physical state."""
    t = np.asarray(time_s, float)
    if t.ndim != 1 or not len(t) or np.any(~np.isfinite(t)) or np.any(np.diff(t) <= 0):
        raise ValueError("Distinct increasing finite snapshot times required")
    if requested_ps is None:
        if len(t) <= 12:
            return np.arange(len(t), dtype=int)
        request = np.array([0., .1, .25, .5, 1., 2., 5., 10., 20., 50., 100., 200.])
        request = request[(request >= t[0]*1e12) & (request <= t[-1]*1e12)]
    else:
        request = np.asarray(requested_ps, float)
        if request.ndim != 1 or not len(request) or np.any(~np.isfinite(request)):
            raise ValueError("Finite requested snapshot times in ps required")
        tolerance = 64*np.finfo(float).eps*max(1., abs(t[0]*1e12), abs(t[-1]*1e12))
        if np.any(request < t[0]*1e12-tolerance) or np.any(request > t[-1]*1e12+tolerance):
            raise ValueError("Requested snapshot lies outside the saved physical window")
    selected = [0, len(t)-1]
    selected += [int(np.argmin(abs(t*1e12-r))) for r in request]
    return np.asarray(sorted(set(selected)), dtype=int)


def retained_energy_J(manifest):
    source = manifest.get("source", {})
    for obj in (source, manifest):
        for key in ("retained_energy_J", "retained_photon_energy_J"):
            if key in obj:
                value = float(obj[key])
                if np.isfinite(value) and value > 0:
                    return value
                raise ValueError("Retained photon energy must be finite and positive")
    raise ValueError("manifest.source.retained_energy_J is required for physical energy normalization")


def observation_text(manifest, final_time_s):
    obs = manifest.get("observation", {})
    status = str(obs.get("status", manifest.get("status", ""))).upper()
    if obs.get("trigger_confirmed") is True:
        value = obs.get("trigger_time_s", obs.get("confirmed_trigger_time_s"))
        if value is None or not np.isfinite(float(value)):
            return "Gatillo confirmado; tiempo no guardado"
        return "Gatillo confirmado: "+time_label_ps(value)
    if "FAIL" in status or "ERROR" in status or "INCOMPLETE" in status:
        return "Integracion incompleta: sin latencia admitida"
    if "CENSOR" in status or "CEILING" in status:
        return "Sin gatillo confirmado; censura a "+time_label_ps(final_time_s)
    return "Sin gatillo confirmado; ventana guardada "+time_label_ps(final_time_s)


def prepare_view(manifest, trajectory, snapshots, requested_ps=None):
    """Validate saved physical arrays and convert units only, without smoothing."""
    history = {}
    t = np.asarray(trajectory["time_s"], float)
    if t.ndim != 1 or not len(t) or np.any(~np.isfinite(t)) or np.any(np.diff(t) <= 0):
        raise ValueError("Finite increasing trajectory time_s required")
    history["time_ps"] = t*1e12
    for key in ("Vout_V", "Vdev_V", "Ib_A", "Is_A", "vc_V", "energy_total_J",
            "energy_defect_J", "escape_power_W", "domain_energy_J", "instant_power_residual_W",
            "external_energy_integral_J"):
        value = np.asarray(trajectory[key], float)
        if value.shape != t.shape or np.any(~np.isfinite(value)):
            raise ValueError("Trajectory field must match saved finite times: "+key)
        history[key] = value
    baseline = float(manifest.get("observation", {}).get("baseline_Vout_V", 0.))
    if not np.isfinite(baseline):
        raise ValueError("Finite declared Vout baseline required")
    history["Vout_relative_mV"] = (history["Vout_V"]-baseline)*1e3
    history["Vdev_mV"] = history["Vdev_V"]*1e3
    history["Ib_uA"], history["Is_uA"] = history["Ib_A"]*1e6, history["Is_A"]*1e6
    retained = retained_energy_J(manifest)
    history["energy_defect_percent_retained"] = history["energy_defect_J"]/retained*100
    absolute_work = np.r_[abs(history["external_energy_integral_J"][0]),
        np.abs(np.diff(history["external_energy_integral_J"]))].cumsum()
    history["absolute_external_work_J"] = absolute_work
    history["energy_defect_percent_global"] = np.abs(history["energy_defect_J"])/(retained+absolute_work)*100
    xy = np.asarray(snapshots["coordinates_m"], float)
    ts = np.asarray(snapshots["time_s"], float)
    if xy.ndim != 2 or xy.shape[1] != 2 or len(xy) < 3 or np.any(~np.isfinite(xy)):
        raise ValueError("Finite planar coordinates_m required")
    indices = nearest_snapshots(ts, requested_ps)
    if ts[0] < t[0]-1e-25 or ts[-1] > t[-1]+1e-25:
        raise ValueError("Snapshots must lie within the saved trajectory")
    maps = {"coordinates_nm": xy*1e9, "time_ps": ts*1e12, "selected_indices": indices}
    for key in ("delta_J", "Te_K", "Tph_K", "u_excess_J_m3", "electron_energy_J_m3", "phonon_energy_J_m3"):
        if key not in snapshots:
            if key in ("delta_J", "Te_K", "Tph_K"):
                raise ValueError("Required physical snapshot field absent: "+key)
            continue
        value = np.asarray(snapshots[key])
        if value.shape != (len(ts), len(xy)) or np.any(~np.isfinite(value)):
            raise ValueError("Snapshot field must match (time,node): "+key)
        if key != "delta_J" and np.iscomplexobj(value):
            raise ValueError("Only gap snapshots may be complex: "+key)
        maps[key] = value
    maps["gap_meV"] = abs(maps["delta_J"])/MEV_J
    if np.any(maps["Te_K"] < 0) or np.any(maps["Tph_K"] < 0):
        raise ValueError("Physical temperature cannot be negative")
    if "triangles" in snapshots:
        triangles = np.asarray(snapshots["triangles"])
        if triangles.ndim != 2 or triangles.shape[1] != 3 or triangles.dtype.kind not in "iu" or np.any(triangles < 0) or np.any(triangles >= len(xy)):
            raise ValueError("Valid saved triangle indices required")
        maps["triangles"] = triangles.astype(int)
    return history, maps


def figure_header(fig, title, manifest, final_time_s):
    source = manifest.get("source", {})
    wavelength = source.get("wavelength_nm")
    prefix = "" if wavelength is None else f"{float(wavelength):g} nm | "
    if manifest.get("admission", {}).get("synthetic_plot_smoke_only") is True:
        prefix = "FORMATO SINTETICO | "+prefix
    fig.suptitle(prefix+title+"\n"+observation_text(manifest, final_time_s), fontsize=11)


def save_figure(fig, name, output, products, caption, quantities, *, dpi, limits=None):
    files = {}
    for suffix in ("png", "pdf"):
        path = output/(name+"."+suffix)
        fig.savefig(path, dpi=dpi, facecolor="white", bbox_inches="tight")
        files[suffix] = {"path": str(path), "sha256": sha256(path)}
    products.append({"name": name, "files": files, "caption": caption,
        "quantities": quantities, "common_color_limits": limits})
    plt.close(fig)


def add_trigger_lines(axes, manifest):
    obs = manifest.get("observation", {})
    for ax in np.asarray(axes, dtype=object).ravel():
        if obs.get("trigger_confirmed") is True:
            crossing = obs.get("first_crossing_time_s")
            confirmed = obs.get("trigger_time_s", obs.get("confirmed_trigger_time_s"))
            if crossing is not None and np.isfinite(float(crossing)):
                ax.axvline(float(crossing)*1e12, color="tab:green", ls=":", lw=.9, label="Primer cruce")
            if confirmed is not None and np.isfinite(float(confirmed)):
                ax.axvline(float(confirmed)*1e12, color="tab:green", ls="--", lw=.9, label="Confirmacion")



def add_readout_detail(ax, manifest, history):
    """Keep the physical readout visible when the device voltage or trigger dominates."""
    pulse = np.asarray(history["Vout_relative_mV"])
    peak = float(np.max(np.abs(pulse)))
    device_peak = float(np.max(np.abs(history["Vdev_mV"])))
    obs = manifest.get("observation", {})
    threshold = obs.get("threshold_voltage_V", obs.get("trigger_voltage_V"))
    threshold_mV = None if threshold is None else abs(float(threshold))*1e3
    small_pulse = threshold_mV is not None and np.isfinite(threshold_mV) and peak < .2*threshold_mV
    if not (device_peak > 5*max(peak, np.finfo(float).tiny) or small_pulse):
        return None
    scale, unit = (1e3, r"$\mu$V") if peak < .01 else (1., "mV")
    detail = ax.inset_axes([.08, .53, .40, .35])
    detail.plot(history["time_ps"], pulse*scale, color="blue", lw=.9)
    detail.set_xlim(float(history["time_ps"][0]), float(history["time_ps"][-1]))
    detail.set_title("Lectura menos basal", fontsize=8)
    detail.set_ylabel(f"Vout - Vout,DC ({unit})", fontsize=7)
    detail.set_xlabel("Tiempo (ps)", fontsize=7)
    detail.tick_params(labelsize=7)
    detail.grid(alpha=.18)
    add_trigger_lines([detail], manifest)
    return detail

def trace_figures(manifest, history, output, products, dpi):
    t = history["time_ps"]
    fig, axes = plt.subplots(2, 1, figsize=(7.3, 6.), sharex=True, layout="constrained")
    axes[0].plot(t, history["Vout_relative_mV"], label=r"$V_{\rm out}-V_{\rm out,DC}$")
    axes[0].plot(t, history["Vdev_mV"], label=r"$V_{\rm dev}$", alpha=.75)
    observation = manifest.get("observation", {})
    threshold = observation.get("threshold_voltage_V", observation.get("trigger_voltage_V"))
    if threshold is not None and np.isfinite(float(threshold)):
        axes[0].axhline(float(threshold)*1e3, color=".4", ls="--", lw=.8,
            label=f"Umbral: {float(threshold)*1e3:g} mV")
    axes[0].set_ylabel("Voltaje (mV)")
    axes[0].text(.02, .97, f"Vout inicial absoluto = {history['Vout_V'][0]*1e9:.4g} nV",
        transform=axes[0].transAxes, va="top", fontsize=8)
    axes[1].plot(t, history["Ib_uA"], label=r"$I_b$")
    axes[1].plot(t, history["Is_uA"], label=r"$I_s$")
    axes[1].set_ylabel(r"Corriente ($\mu$A)")
    axes[1].set_xlabel("Tiempo desde transferencia efectiva del foton (ps)")
    add_trigger_lines(axes, manifest)
    detail = add_readout_detail(axes[0], manifest, history)
    for ax in axes:
        ax.grid(alpha=.18); ax.legend(fontsize=8)
    if detail is not None:
        axes[0].legend(fontsize=8, loc="lower left")
    figure_header(fig, "Pulso electrico y redistribucion de corriente", manifest, t[-1]*1e-12)
    save_figure(fig, "01_pulse_and_currents", output, products,
        "El pulso de lectura es R_load*(Ib-Is), menos la referencia DC declarada. Vdev es la diferencia de potencial entre cortes completos. Las corrientes son estados del circuito CM completo. Las lineas de deteccion proceden del manifiesto; no se infiere latencia a partir del maximo del pulso. El tiempo comienza en la transferencia efectiva tras la cascada inicial. Si el eje principal aplana la lectura, el recuadro muestra Vout menos basal con escala propia (mV o microV), sin alterar valores, tiempos o gatillo.",
        ["Vout-relative (mV)", "Vdev (mV)", "Ib (microA)", "Is (microA)", "time (ps)"], dpi=dpi)

    fig, axes = plt.subplots(2, 2, figsize=(9., 6.4), sharex=True, layout="constrained")
    axes[0,0].plot(t, (history["domain_energy_J"]-history["domain_energy_J"][0])/EV_J)
    axes[0,0].set_ylabel(r"$U_{\rm film}(t)-U_{\rm film}(0)$ (eV)")
    axes[0,0].set_title("Energia del dominio resuelto")
    axes[0,1].plot(t, (history["energy_total_J"]-history["energy_total_J"][0])/EV_J)
    axes[0,1].set_ylabel(r"$E_{\rm total}(t)-E_{\rm total}(0)$ (eV)")
    axes[0,1].set_title("Energia total registrada")
    axes[1,0].plot(t, history["energy_defect_percent_global"], label="|defecto| / (Eret + Σ|ΔWext|)")
    axes[1,0].plot(t, history["energy_defect_percent_retained"], ls="--", label="defecto firmado / Eret")
    axes[1,0].axhline(1., color=".4", lw=.7, ls=":", label="Referencia global 1%")
    axes[1,0].set_ylabel("Defecto energético (%)")
    axes[1,0].set_title(f"Eret = {retained_energy_J(manifest)/EV_J:.5g} eV; Wext firmado guardado")
    axes[1,0].legend(fontsize=7)
    axes[1,1].plot(t, history["escape_power_W"]*1e9, label="Escape fononico")
    axes[1,1].plot(t, history["instant_power_residual_W"]*1e9, label="Residuo de potencia")
    axes[1,1].set_ylabel("Potencia (nW)"); axes[1,1].legend(fontsize=8)
    for ax in axes.ravel():
        ax.grid(alpha=.18); ax.set_xlabel("Tiempo (ps)")
    add_trigger_lines(axes, manifest)
    figure_header(fig, "Energia y balance independiente guardado", manifest, t[-1]*1e-12)
    save_figure(fig, "02_energy_balance", output, products,
        "Las energias se restan respecto del primer estado guardado. Se muestran el valor absoluto del defecto independiente dividido por Eret + Σ|ΔWext|, como en el auditor global, y su valor firmado dividido por Eret como referencia adicional. ΔWext son diferencias del trabajo externo acumulado guardado, no trapecios de potencias ni el valor absoluto del trabajo neto. Las curvas no conceden admision. La potencia de escape describe transferencia hacia el sustrato. Un grafico por si solo no acredita conservacion o convergencia.",
        ["domain energy increment (eV)", "saved total energy increment (eV)",
         "absolute independent defect / (Eret + cumulative absolute saved external work increments) (%)",
          "signed independent defect / retained photon energy (%)",
         "phonon escape and instantaneous power residual (nW)"], dpi=dpi)


def map_figure(manifest, maps, key, label, name, output, products, dpi, caption,
        *, cmap="viridis", signed=False, center_width_nm=None):
    values = np.asarray(maps[key], float)
    indices = maps["selected_indices"]
    low, high = float(values.min()), float(values.max())
    if signed and low < 0:
        high = max(abs(low), abs(high)); low = -high; cmap = "coolwarm"
    elif key == "gap_meV" or key.endswith("J_m3"):
        low = min(0., low)
    if high <= low:
        high = low+max(abs(low)*1e-12, 1e-30)
    columns = min(3, len(indices)); rows = (len(indices)+columns-1)//columns
    fig, axes = plt.subplots(rows, columns, figsize=(3.55*columns, 2.6*rows+.65),
        squeeze=False, layout="constrained")
    xy = maps["coordinates_nm"]
    tri = Triangulation(xy[:,0], xy[:,1], maps.get("triangles"))
    active = []
    for ax, index in zip(axes.ravel(), indices):
        if key == "gap_meV":
            bounds = None
            if center_width_nm is not None:
                center = .5*(xy[:,0].min()+xy[:,0].max())
                bounds = (center-center_width_nm/2, center+center_width_nm/2, xy[:,1].min(), xy[:,1].max())
            image = draw_complex_gap(ax, tri, maps["delta_J"][index], unit_J=MEV_J,
                vmin=low, vmax=high, cmap=cmap, bounds=bounds)
        else:
            image = ax.tripcolor(tri, values[index], shading="gouraud", cmap=cmap,
                vmin=low, vmax=high, rasterized=True)
        ax.set_aspect("equal"); ax.set_xlabel("x (nm)"); ax.set_ylabel("y (nm)")
        ax.set_title(time_label_ps(maps["time_ps"][index]*1e-12), fontsize=10)
        if center_width_nm is not None:
            center = .5*(xy[:,0].min()+xy[:,0].max())
            ax.set_xlim(center-center_width_nm/2, center+center_width_nm/2)
        active.append(ax)
    for ax in axes.ravel()[len(indices):]:
        ax.set_visible(False)
    fig.colorbar(image, ax=active, label=label, shrink=.85)
    figure_header(fig, label+" | estados guardados", manifest, maps["time_ps"][-1]*1e-12)
    extra = " Limites absolutos comunes sobre TODOS los estados guardados; no se normaliza cada mapa por su maximo. Se usan los instantes realmente guardados, sin interpolacion temporal."
    if key == "gap_meV":
        extra += " Gap: |P1 Delta complejo| muestreado en raster de hasta 601 puntos por eje; celdas sin suavizado. La resolucion visual no resuelve un nucleo fisico ni sustituye el diagnostico de hotbelt."
    if center_width_nm is not None:
        extra += f" Ampliacion central de {center_width_nm:g} nm; no define una ventana termica ni modifica la simulacion."
    save_figure(fig, name, output, products, caption+extra,
        [label, "coordinates (nm)", "physical time (ps)"], dpi=dpi, limits=[low, high])


def run(run_dir, output_dir, snapshots_ps=None, *, dpi=THESIS_DPI):
    start = time.monotonic()
    run_dir, output = in_scratch(run_dir), in_scratch(output_dir)
    if output.exists():
        raise FileExistsError("Preserve prior plots; use a fresh --output-dir")
    if not 72 <= int(dpi) <= 600:
        raise ValueError("dpi must be between 72 and 600")
    from sandbox.stage4_core.parallel_runtime import linux_resources, resource_budget
    budget = resource_budget(linux_resources(), .9, .9, 1)
    old_affinity = os.sched_getaffinity(0)
    os.sched_setaffinity(0, {budget["coordinator_cpu"]})
    try:
        inputs = [run_dir/"manifest.json", run_dir/"trajectory.npz", run_dir/"snapshots.npz",
            COMPLEX_GAP_SOURCE]
        manifest = json.loads(inputs[0].read_text())
        with np.load(inputs[1], allow_pickle=False) as z:
            trajectory = {key:z[key].copy() for key in z.files}
        with np.load(inputs[2], allow_pickle=False) as z:
            wanted = {"coordinates_m", "time_s", "delta_J", "Te_K", "Tph_K",
                "u_excess_J_m3", "electron_energy_J_m3", "phonon_energy_J_m3", "triangles"}
            snapshots = {key:z[key].copy() for key in z.files if key in wanted}
        numeric_bytes = sum(v.nbytes for v in trajectory.values())+sum(v.nbytes for v in snapshots.values())
        if 2*numeric_bytes+512*1024**2 > budget["memory_limit_bytes"]:
            raise MemoryError("Input arrays plus plotting headroom exceed 90% available RAM")
        history, maps = prepare_view(manifest, trajectory, snapshots, snapshots_ps)
        output.mkdir(parents=True, exist_ok=False)
        apply_thesis_style()
        products = []
        trace_figures(manifest, history, output, products, int(dpi))
        map_figure(manifest, maps, "gap_meV", r"$|\Delta|$ (meV)", "03_gap_montage",
            output, products, int(dpi),
            "El gap mide la integridad del condensado. Su depresion y eventual conexion transversal permiten estudiar formacion de una hotbelt y paso al estado resistivo.")
        if np.ptp(maps["coordinates_nm"][:,0]) > 160:
            map_figure(manifest, maps, "gap_meV", r"$|\Delta|$ (meV)", "04_gap_central_montage",
                output, products, int(dpi),
                "La ampliacion muestra el condensado alrededor del punto central de absorcion y toda la anchura de la cinta.",
                center_width_nm=120.)
        map_figure(manifest, maps, "Te_K", r"$T_e^{eq}$ (K)", "05_electron_temperature_montage",
            output, products, int(dpi),
            "Temperatura electronica equivalente asociada a la energia de poblaciones. Es un diagnostico fisico; esta etiqueta no afirma que la distribucion no termica sea Fermi-Dirac.", cmap="inferno")
        map_figure(manifest, maps, "Tph_K", r"$T_{ph}^{eq}$ (K)", "06_phonon_temperature_montage",
            output, products, int(dpi),
            "Temperatura fononica equivalente. Permite localizar la energia transferida al subsistema fononico antes de su escape al sustrato.", cmap="inferno")
        for key, name, label, caption in (
            ("u_excess_J_m3", "07_excess_energy_density_montage", "Exceso de energia (J/m3)",
             "Exceso local de energia guardado respecto de la referencia declarada por la simulacion. Muestra expansion energetica espacial sin imponer una ventana de dinamica termica."),
            ("electron_energy_J_m3", "08_electron_energy_density_montage", "Energia electronica (J/m3)",
             "Densidad electronica registrada, con la referencia y convencion del manifiesto de la corrida. Incluye solo las contribuciones definidas en ese campo."),
            ("phonon_energy_J_m3", "09_phonon_energy_density_montage", "Energia fononica (J/m3)",
             "Densidad de energia fononica registrada. Su expansion y reduccion permiten separar redistribucion interna y escape hacia el sustrato.")):
            if key in maps:
                map_figure(manifest, maps, key, label, name, output, products, int(dpi),
                    caption, cmap="magma", signed=True)
        receipt = {"schema": SCHEMA, "status":"SAVED_PHYSICAL_FIGURES_NO_NEW_MODEL_ADMISSION",
            "run_dir":str(run_dir), "output_dir":str(output), "pipeline":str(Path(__file__).resolve()),
            "source_sha256":sha256(__file__),
            "input_sha256":{str(path):sha256(path) for path in inputs},
            "requested_snapshots_ps":snapshots_ps,
            "selected_saved_snapshots_ps":maps["time_ps"][maps["selected_indices"]].tolist(),
            "observation":manifest.get("observation", {}), "source_admission":manifest.get("admission"),
            "model_admission_inferred":False, "retained_energy_J":retained_energy_J(manifest),
            "snapshot_interpolation":"None in time; gap: magnitude after complex spatial P1 sampled on nearest-cell raster; scalar maps: nodal Gouraud rendering",
            "gap_raster_max_samples_per_axis":601,
            "energy_defect_normalizations":["absolute defect / (Eret + cumulative absolute saved external work increments)", "signed defect / Eret"],
            "style":"pysnspd.plotting.style.apply_thesis_style", "dpi":int(dpi),
            "resources":{"processes":1, "selected_cpu":budget["coordinator_cpu"],
                "numeric_array_bytes":numeric_bytes, "memory_limit_bytes":budget["memory_limit_bytes"]},
            "products":products, "elapsed_s":time.monotonic()-start}
        (output/"plot_manifest.json").write_text(json.dumps(receipt, indent=2)+"\n")
        print(json.dumps({"status":receipt["status"], "products":len(products),
            "output_dir":str(output), "elapsed_s":receipt["elapsed_s"]}))
        return receipt
    finally:
        os.sched_setaffinity(0, old_affinity)


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--snapshots-ps", nargs="+", type=float, default=None,
        help="Nearest saved times; initial and final states are always included")
    parser.add_argument("--dpi", type=int, default=THESIS_DPI)
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    run(args.run_dir, args.output_dir, args.snapshots_ps, dpi=args.dpi)
