"""Physical, explicitly normalized figures for the stage-4 closure report.

This is a read-only reduction of saved fields: no extra physical solve, no
fitted time constant and no reconstructed photon transient. Harmonic curves
are periodic real parts of peak phasors with the stored exp(-i omega t)
convention. The voltage specified in the campaign is the *source* voltage.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import numpy as np


KB = 1.380649e-23
HBAR = 1.054571817e-34
EC = 1.602176634e-19
TC = 8.65
D = 0.5e-4
RS = 608.0
L0_NM = np.sqrt(HBAR * D / (2 * KB * TC)) * 1e9
V0 = KB * TC / EC
I0 = V0 / (2 * RS)
COLORS = ["#2363a2", "#d17424", "#347e63", "#9564ad"]
INPUTS: dict[str, str] = {}


def register(path):
    path = Path(path)
    INPUTS[path.as_posix()] = hashlib.sha256(path.read_bytes()).hexdigest()
    return path


def read_json(path):
    return json.loads(register(path).read_text(encoding="utf-8"))


def read_npz(path):
    with np.load(register(path)) as data:
        return {key: data[key] for key in data.files}


def phasor(value):
    return np.asarray(value["real"]) + 1j * np.asarray(value["imag"])


def tri_for(fields):
    xy = fields["coordinates_bar"] * L0_NM
    return mtri.Triangulation(xy[:, 0], xy[:, 1])


def line(fields, values):
    tri = tri_for(fields)
    x = np.linspace(tri.x.min(), tri.x.max(), 401)
    y = np.zeros_like(x)
    if np.iscomplexobj(values):
        out = (mtri.LinearTriInterpolator(tri, values.real)(x, y)
               + 1j * mtri.LinearTriInterpolator(tri, values.imag)(x, y))
    else:
        out = mtri.LinearTriInterpolator(tri, values)(x, y)
    return x, np.asarray(out)


def style(ax, title, xlabel=None, ylabel=None):
    ax.set_title(title, loc="left", fontweight="bold", pad=9)
    if xlabel:
        ax.set_xlabel(xlabel)
    if ylabel:
        ax.set_ylabel(ylabel)
    ax.grid(alpha=0.19)
    ax.spines[["top", "right"]].set_visible(False)


def map_field(ax, fields, values, title, label, limits=None, cmap="viridis"):
    image = ax.tripcolor(tri_for(fields), values, shading="gouraud", cmap=cmap,
                         vmin=None if limits is None else limits[0],
                         vmax=None if limits is None else limits[1], rasterized=True)
    style(ax, title, "x [nm]", "y [nm]")
    ax.grid(False)
    ax.set_aspect("equal")
    ax.axhline(0, color="white", linestyle=":", linewidth=0.8, alpha=0.7)
    plt.colorbar(image, ax=ax, pad=0.025, shrink=0.78, label=label)
    return image


def finish(fig, path, title, footer, bottom=0.14):
    fig.suptitle(title, x=0.06, ha="left", y=0.985, fontsize=16, fontweight="bold")
    fig.text(0.06, 0.025, footer, va="bottom", fontsize=9, color="#38475a")
    fig.subplots_adjust(top=0.89, bottom=bottom, left=0.095, right=0.95,
                        hspace=0.52, wspace=0.42)
    fig.savefig(path, dpi=180, facecolor="white")
    plt.close(fig)


def load_case(data_root, name):
    root = data_root / "responses" / "cases" / name
    fields = read_npz(root / "coupled_fields.npz")
    results = read_json(root / "results.json")
    state = phasor(results["circuit_state_peak"])
    voltage = complex(phasor(results["impedance_ohm"])) * state[1]
    assert results["circuit_time_convention"] == "exp(-iwt)"
    return dict(name=name, fields=fields, results=results, state=state, vdev=voltage)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-root", type=Path,
                        default=Path("docs/implementation/stage4/closure_20260924/data"))
    parser.add_argument("--thermal-root", type=Path,
                        default=Path("docs/implementation/stage4/final_kwt_20260924/thermal"))
    parser.add_argument("--output-root", type=Path,
                        default=Path("docs/implementation/stage4/closure_20260924/figures"))
    args = parser.parse_args()
    out = args.output_root
    out.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10.5,
                         "axes.titlesize": 11, "axes.labelsize": 11,
                         "legend.fontsize": 9.5, "lines.linewidth": 2.0,
                         "axes.formatter.useoffset": False})
    refs = {}
    for name in ("phase2", "phase8"):
        root = args.data_root / "references" / "cases" / name
        refs[name] = read_npz(root / "reference.npz")
    cases = {name: load_case(args.data_root, name) for name in
             ("phase2_w020", "phase8_w005", "phase8_w020", "phase8_w080", "phase8_refined")}
    gap0 = float(np.max(abs(refs["phase2"]["delta_bar"])))
    metrics = {"scales": {"length_nm": L0_NM, "voltage_V": V0,
                           "current_A": I0, "gap0_over_kBTc": gap0},
               "case_metrics": {}, "thermal_metrics": {}}
    captions = []

    fig, axes = plt.subplots(2, 2, figsize=(11.6, 8.4))
    for col, (name, phase) in enumerate((("phase2", 2), ("phase8", 8))):
        fields = refs[name]
        current = cases[f"{name}_w020"]["results"]["reference_dc_current_A"] * 1e6
        map_field(axes[0, col], fields, abs(fields["delta_bar"]) / gap0,
                  f"{'ab'[col]}) Fase total {phase} rad; $I_{{DC}}$ = {current:.3f} µA",
                  r"$|\Delta|/\Delta_0$ [1]", (0.974, 1.0))
        x, amplitude = line(fields, abs(fields["delta_bar"]) / gap0)
        axes[1, 0].plot(x, amplitude, color=COLORS[col], label=f"{phase} rad")
        x, delta = line(fields, fields["delta_bar"])
        q = np.gradient(np.unwrap(np.angle(delta)), x)
        axes[1, 1].plot(x, q, color=COLORS[col], label=f"{phase} rad")
        metrics["case_metrics"][name] = dict(
            current_uA=current, minimum_gap_over_gap0=float(min(abs(fields["delta_bar"])) / gap0),
            centerline_mean_q_per_nm=float(np.mean(q)),
            phase_advance_rad=float(np.unwrap(np.angle(delta))[-1]-np.unwrap(np.angle(delta))[0]))
    style(axes[1, 0], "c) Supresión del condensado en el eje", "x [nm]",
          r"$|\Delta(x,0)|/\Delta_0$ [1]")
    style(axes[1, 1], "d) Gradiente longitudinal de fase", "x [nm]", r"$q_x(x,0)$ [nm$^{-1}$]")
    for ax in axes[1]:
        ax.legend(frameon=False)
    finish(fig, out / "01_condensado_polarizado.png", "La corriente deprime el condensado en el interior",
           "Cinta 160 × 80 nm; T = 0,9 K; contactos con amplitud fija y bordes laterales naturales; A = 0.\n"
           "Δ₀ es el gap uniforme de la misma suma de Matsubara. Cortes y qₓ = ∂ₓarg Δ en y = 0; interpolación triangular.")
    captions.append(dict(file="01_condensado_polarizado.png",
        title="Condensado estacionario bajo corriente",
        quantities="Mapas y corte de |Delta|/Delta0, sin sustracción; q_x=d(arg Delta)/dx para A=0, en nm^-1. Corriente DC de contacto en microamperios.",
        case="Referencias phase2 y phase8; 160 x 80 nm, 0.9 K, sin fotón.",
        interpretation="Cuadruplicar la fase casi cuadruplica I, pero aumenta más la depresión del gap: el transporte ya modifica la rigidez del condensado. No es una curva de corriente crítica."))

    observations = [read_npz(args.thermal_root / f"observation_{i:03d}.npz") for i in range(6)]
    time = np.array([float(z["time_ps"]) for z in observations])
    amp, phase = [], []
    for z in observations:
        d0 = z["baseline_gap"]
        amp.append((abs(z["refined_amplitude_gap"])-abs(d0))/abs(d0))
        phase.append(np.angle(z["refined_angular_phase_gap"] / d0))
    amp = np.asarray(amp)
    phase = np.asarray(phase)
    fig, axes = plt.subplots(2, 2, figsize=(11.6, 8.5))
    for color, idx in zip(COLORS, (0, 2, 3, 5)):
        x, v = line(observations[idx], amp[idx] * 1e6)
        axes[0, 0].plot(x, v, color=color, label=f"{time[idx]:g} ps")
        x, v = line(observations[idx], phase[idx] * 1e3)
        axes[0, 1].plot(x, v, color=color, label=f"{time[idx]:g} ps")
    style(axes[0, 0], "a) Perturbación inicial de amplitud", "x [nm]", r"$(|\Delta|-\Delta_0)/\Delta_0$ [ppm]")
    style(axes[0, 1], "b) Perturbación inicial de fase", "x [nm]", r"$\arg(\Delta/\Delta_0)$ [mrad]")
    for ax in axes[0]:
        ax.legend(frameon=False, ncol=2)
    for color, values, label in ((COLORS[0], amp, "Amplitud"), (COLORS[1], phase, "Fase")):
        mass = observations[0]["area_weights"]
        norm = np.sqrt(np.sum(values**2 * mass[None, :], axis=1) / mass.sum())
        axes[1, 0].semilogy(time, norm / norm[0], "o-", color=color, label=label)
        metrics["thermal_metrics"][label] = dict(final_norm_over_initial=float(norm[-1]/norm[0]))
    history = read_json(args.thermal_root / "step_history.json")
    for color, trajectory, label in ((COLORS[0], "refined_amplitude", "Amplitud"),
                                      (COLORS[1], "refined_angular_phase", "Fase")):
        rows = [v for v in history if v["trajectory"] == trajectory]
        initial = read_json(args.thermal_root / "observation_000.json")["states"][trajectory]["energy_excess"]
        t = np.r_[0, [v["time_ps"] for v in rows]]
        e = np.r_[initial, [v["energy_excess"] for v in rows]]
        axes[1, 1].semilogy(t, e / initial, color=color, label=label)
        metrics["thermal_metrics"][label]["final_free_energy_over_initial"] = float(e[-1]/initial)
    style(axes[1, 0], "c) Dos escalas de relajación", "Tiempo aceptado [ps]", r"Norma RMS / su valor inicial [1]")
    style(axes[1, 1], "d) Pérdida de energía libre", "Tiempo aceptado [ps]", r"$(F-F_0)/(F(0)-F_0)$ [1]")
    for ax in axes[1]:
        ax.legend(frameon=False)
    finish(fig, out / "02_relajacion_real.png", "La fase se relaja antes que la amplitud",
           "Dos ensayos independientes a T fija = 0,9 K, sin corriente DC ni fotón; Euler KWT, paso refinado.\n"
           "RMS ponderada por áreas de Voronoi, respecto al equilibrio uniforme; cada norma y exceso de F usa su propio valor inicial.")
    captions.append(dict(file="02_relajacion_real.png", title="Relajación temporal real de perturbaciones suaves",
        quantities="Cortes y=0: exceso relativo de amplitud en ppm y fase relativa en mrad. RMS ponderada por áreas y exceso de energía libre, ambos divididos por su propio valor inicial.",
        case="Dos excitaciones distintas sobre equilibrio uniforme a 0.9 K, sin corriente ni fotón. Pasos aceptados de Euler KWT; curvas de energía del historial completo; normas de seis observaciones.",
        interpretation="La redistribución de fase casi termina en 1 ps, mientras la amplitud conserva una fracción importante. F es energía libre a temperatura fija, no energía interna de poblaciones acopladas."))

    fig, axes = plt.subplots(2, 2, figsize=(11.6, 8.4))
    for col, name in enumerate(("phase2_w020", "phase8_w020")):
        c = cases[name]
        values = abs(c["fields"]["radial_response"] * c["vdev"]) * 1e12
        map_field(axes[0, col], c["fields"], values,
                  f"{'ab'[col]}) {2 if col == 0 else 8} rad; Ω = 0,20", r"$|\delta|\Delta||_{ac}$ [peV]", (0, 0.40))
    for color, name in zip(COLORS, ("phase8_w005", "phase8_w020", "phase8_w080")):
        c = cases[name]
        omega = c["results"]["omega"]
        f = c["fields"]
        x, radial = line(f, f["radial_response"] * c["vdev"])
        axes[1, 0].semilogy(x[1:-1], abs(radial[1:-1])*1e12, color=color, label=f"Ω = {omega:.2f}")
        x, potential = line(f, f["potential_response"] * c["vdev"])
        axes[1, 1].plot(x, abs(potential)*1e12, color=color, label=f"Ω = {omega:.2f}")
        metrics["case_metrics"].setdefault(name, {}).update(dict(
            device_voltage_peak_pV=float(abs(c["vdev"])*1e12),
            maximum_gap_amplitude_peV=float(max(abs(f["radial_response"]*c["vdev"]))*1e12),
            maximum_potential_amplitude_pV=float(max(abs(f["potential_response"]*c["vdev"]))*1e12),
            maximum_phase_amplitude_urad=float(max(abs(f["angular_response"]*c["vdev"]/V0/abs(f["gap_reference"]))))*1e6))
    style(axes[1, 0], "c) Respuesta de amplitud; sesgo de 8 rad", "x [nm]", r"$|\delta|\Delta|(x,0)|_{ac}$ [peV]")
    style(axes[1, 1], "d) Potencial alterno; sesgo de 8 rad", "x [nm]", r"$|\delta\phi(x,0)|_{ac}$ [pV]")
    for ax in axes[1]:
        ax.legend(frameon=False)
    finish(fig, out / "03_respuesta_espacial.png", "Sesgo y frecuencia cambian la respuesta del condensado",
           "Amplitudes pico de perturbaciones lineales respecto al estado DC; fuente alterna = 1 µV, con el circuito completo.\n"
           "Ω = ħω/(kBTc). 1 peV = 10⁻¹² eV; 1 pV = 10⁻¹² V. |δφ| muestra amplitud, no signo ni campo eléctrico.")
    captions.append(dict(file="03_respuesta_espacial.png", title="Respuesta espacial armónica del condensado y potencial",
        quantities="Módulo del fasor de la variación radial de Delta en peV y del potencial escalar en pV; referencias DC sustraídas. Amplitudes pico de la misma fuente alterna de 1 microvoltio, transferida mediante el circuito completo.",
        case="Dos mapas a Omega=0.2 sobre phase2/phase8; perfiles y=0 sobre phase8 a Omega=0.05, 0.2, 0.8. Respuesta candidata en 6 modos; sin fotón.",
        interpretation="La corriente facilita la perturbación de la amplitud; la respuesta de la fuente cae con frecuencia. Esa caída mezcla la física intrínseca con el filtrado circuital y no se interpreta como una relajación propia del material."))

    c = cases["phase8_w020"]
    r = c["results"]
    period = 2*np.pi/r["angular_frequency_per_ps"]
    t = np.linspace(0, 2*period, 600)
    rotation = np.exp(-1j*r["angular_frequency_per_ps"]*t)
    fig, axes = plt.subplots(2, 2, figsize=(11.6, 8.4))
    axes[0, 0].plot(t, np.real(r["bias_voltage_peak_V"]*rotation)*1e6, color="#7c8794")
    style(axes[0, 0], "a) Fuente: modulación sobre el sesgo DC", "Tiempo periódico [ps]", r"$\delta V_b$ [µV]")
    for color, name, ls in ((COLORS[0], "phase8_w020", "-"), (COLORS[1], "phase8_refined", "--")):
        ac = cases[name]
        wave = np.real(phasor(ac["results"]["readout_peak_V"])*rotation)*1e12
        axes[0, 1].plot(t, wave, color=color, linestyle=ls, label="Base" if ls == "-" else "Refinado")
    style(axes[0, 1], "b) Lectura del circuito completo", "Tiempo periódico [ps]", r"$\delta V_{out}$ [pV]")
    axes[0, 1].legend(frameon=False)
    axes[1, 0].plot(t, np.real(c["vdev"]*rotation)*1e12, color=COLORS[2])
    style(axes[1, 0], "c) Caída en el tramo espacial resuelto", "Tiempo periódico [ps]", r"$\delta V_{dev}$ [pV]")
    names = ("phase8_w005", "phase8_w020", "phase8_w080")
    frequency = np.array([cases[n]["results"]["angular_frequency_per_ps"]/(2*np.pi) for n in names])*1e3
    readout = np.array([abs(phasor(cases[n]["results"]["readout_peak_V"])) for n in names])*1e12
    axes[1, 1].loglog(frequency, readout, "o-", color=COLORS[0], label=r"$|\delta V_{out}|$")
    axes[1, 1].loglog(frequency, [abs(cases[n]["vdev"])*1e12 for n in names], "s-", color=COLORS[2], label=r"$|\delta V_{dev}|$")
    style(axes[1, 1], "d) La misma fuente a tres frecuencias", "Frecuencia f [GHz]", "Amplitud pico [pV]")
    axes[1, 1].legend(frameon=False)
    finish(fig, out / "04_circuito_completo.png", "El circuito filtra fuertemente la excitación de alta frecuencia",
           "Estado periódico, no encendido ni pulso de fotón. Paneles a-c: Ω = 0,20; fase DC = 8 rad; I_DC = 13,002 µA.\n"
           "Se grafica Re[fasor × exp(−iωt)] tras restar el sesgo DC. Vout = 50 Ω (Ib − Is); Vdev = Zdev Is.")
    captions.append(dict(file="04_circuito_completo.png", title="Transferencia del circuito completo de la memoria",
        quantities="Variaciones temporales periódicas de fuente (microvoltios), lectura Vout y caída Vdev (picovoltios), con el sesgo DC restado; exp(-i omega t), fasores pico. Amplitudes frente a frecuencia física en GHz.",
        case="Referencia phase8, fuente pico de 1 microvoltio. Rb=10 kiloohmios, Lb=1 microhenrio, Rload=50 ohmios, C=100 picofaradios, inductancia total=10 nanohenrios; Lexterna descuenta una sola vez la porción espacial resuelta.",
        interpretation="La concordancia de Vout base/refinado está favorecida por el fuerte filtrado externo. Esta salida no certifica por sí sola la dissipación interna y no representa una señal de detección."))

    fig, (cut_ax, ax) = plt.subplots(1, 2, figsize=(11.6, 6.9))
    for color, name in zip(COLORS, ("phase8_w005", "phase8_w020", "phase8_w080", "phase8_refined")):
        c = cases[name]
        fields = c["fields"]
        xn = fields["coordinates_bar"][:, 0] * L0_NM
        tail, head = fields["edges"].T
        current = fields["current_response"]
        unique_x = np.unique(xn)
        cuts = (unique_x[:-1] + unique_x[1:]) / 2
        icut = np.array([np.sum(((xn[tail] <= x).astype(int)-(xn[head] <= x).astype(int))*current) for x in cuts])
        contact = abs(xn-xn.min()) <= 1e-10*max(L0_NM, np.ptp(xn))
        ileft = np.sum((contact[tail].astype(int)-contact[head].astype(int))*current)
        relative = abs(icut-ileft)/abs(ileft)*100
        label = f"Ω={c['results']['omega']:.2f}" + (" refinado" if "refined" in name else "")
        cut_ax.plot(cuts, relative, color=color, label=label, linestyle="--" if "refined" in name else "-")
        metrics["case_metrics"].setdefault(name, {})["maximum_crosscut_current_difference_from_left_percent"] = float(max(relative))
    style(cut_ax, "a) Corriente frente al contacto izquierdo", "Posición del corte x [nm]",
          r"$100|I(x)-I_{izq}|/|I_{izq}|$ [%]")
    cut_ax.legend(frameon=False)
    names = ("phase2_w020", "phase8_w005", "phase8_w020", "phase8_refined", "phase8_w080")
    labels = ("2 rad\nΩ=0,20", "8 rad\nΩ=0,05", "8 rad\nΩ=0,20", "8 rad\nΩ=0,20\nrefinado", "8 rad\nΩ=0,80")
    pport = [cases[n]["results"]["film_port_power_per_voltage_squared_W"]*1e6 for n in names]
    pkwt = [cases[n]["results"]["candidate_radial_heat_per_voltage_squared_W"]*1e6 for n in names]
    index = np.arange(len(names))
    ax.bar(index-.19, pport, .36, color=COLORS[0], label=r"Puerto: $\langle P\rangle/|V_{dev}|^2 = \mathrm{Re}(Y)/2$")
    ax.bar(index+.19, pkwt, .36, color=COLORS[1], label="Término radial KWT candidato")
    ax.set_xticks(index, labels)
    style(ax, "b) Potencia media: puerto y condensado", ylabel=r"$\langle P\rangle/|V_{dev,pico}|^2$ [µW/V²]")
    ax.legend(frameon=False, loc="upper right", fontsize=8)
    ax.tick_params(axis="x", labelsize=8.5)
    ax.set_ylim(0, max(pkwt)*1.36)
    for i, (a, b) in enumerate(zip(pport, pkwt)):
        metrics["case_metrics"].setdefault(names[i], {}).update(dict(
            radial_KWT_over_port_power=float(b/a), port_power_per_voltage_squared_uW_per_V2=a,
            radial_heat_per_voltage_squared_uW_per_V2=b))
    finish(fig, out / "05_potencia_interna.png", "La unión disipativa aún no satisface un cierre energético",
           "Izquierda: corriente compleja integrada en cortes transversales, referida al corte junto al contacto izquierdo.\n"
           "Derecha: potencia por voltaje pico AL DISPOSITIVO al cuadrado; no es la potencia real de la fuente de 1 µV.\n"
           "La resta puerto − KWT no es un flujo de calor calculado. El desequilibrio exige resolver el acoplamiento antes de usarlo con fotones.", bottom=0.22)
    captions.append(dict(file="05_potencia_interna.png", title="Desajuste físico del cierre disipativo candidato",
        quantities="Izquierda: desviación de corriente compleja integrada sobre cada corte transversal frente al contacto izquierdo, módulo relativo en porcentaje. Cortes exhaustivos en los puntos medios entre todas las coordenadas x nodales distintas; suma orientada exacta de aristas que cruzan el corte, sin interpolación. Derecha: potencia AC media de puerto Re(Y)/2 y calor radial KWT candidato, ambos por el cuadrado de la amplitud pico de Vdev; unidades microvatios/voltio^2. No se usa normalización por la fuente externa.",
        case="Cinco respuestas armónicas guardadas; todos los casos del candidato y su refinamiento.",
        interpretation="El calor radial candidato supera la potencia de puerto en todos los casos. La diferencia no es un calor de reservorio evaluado y estos datos no establecen balance energético total; no puede ocultarse mediante una tolerancia de Vout."))
    register(Path(__file__))
    receipt = dict(schema="pysnspd.stage4.physical_figures.v1", source_inputs_sha256=INPUTS,
                   physical_parameters=dict(T_K=0.9, Tc_K=TC, D_cm2_per_s=0.5,
                        sheet_resistance_ohm=RS, tau_ee_at_Tc_ps=0.5,
                        tau_ep_at_Tc_ps=2.47, mobility_origin="Inherited phenomenological KWT parameters; not a fully calibrated Korzh-2020 device or photon preparation."),
                   spatial_interpolation="Piecewise-linear interpolation over a display triangulation of saved mesh vertices. No new physical solve.",
                   harmonic_conversion="Vdev_peak=Zdev*Is_peak; delta_Delta_eV=radial_response*Vdev_peak; delta_phi_V=potential_response*Vdev_peak; delta_theta=angular_response*Vdev_peak/(V0*abs(gap_reference)).",
                   metrics=metrics, captions=captions,
                   files=[dict(file=p.name, sha256=hashlib.sha256(p.read_bytes()).hexdigest()) for p in sorted(out.glob("*.png"))])
    (out / "figure_manifest.json").write_text(json.dumps(receipt, indent=2, ensure_ascii=False)+"\n", encoding="utf-8")
    print(json.dumps({"figures": len(captions), "output_root": str(out), "metrics": metrics}, ensure_ascii=False))


if __name__ == "__main__":
    main()
