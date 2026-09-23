"""Cheap graph-port and prescribed-input thesis-circuit diagnostics.

Example (Lext is an explicitly selected test value, not a measured partition):
  python sandbox/stage3_spatial/ports_20260923/diagnose_electrical.py \
    --output-dir tmp/electrical_ports_example --lk-ext-h 7e-9

No superconducting fields, catalogue queries or detector transients are run.
The output directory must be new. Outputs: JSON, NPZ, and one PNG/PDF figure.
"""
from __future__ import annotations

import argparse
from dataclasses import asdict
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import platform
import sys
import time

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))

import numpy as np
import scipy
from scipy.integrate import solve_ivp
from scipy.linalg import expm
import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt

from pysnspd.experimental.electrical_ports import (
    ThesisCircuitParameters, solve_potential,
)


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def run(output_dir, lk_ext_h):
    started = time.perf_counter()
    circuit = ThesisCircuitParameters(Lk_ext_H=lk_ext_h)
    out = Path(output_dir).resolve()
    if out.exists():
        raise FileExistsError("Choose a fresh output directory; existing evidence is not overwritten")
    out.mkdir(parents=True)
    sources = [
        Path(__file__).resolve(),
        ROOT / "pysnspd/experimental/electrical_ports.py",
        ROOT / "docs/modelo_v0_4/actualizaciones/circuito_memoria_20260922.md",
    ]
    hashes = {str(p.relative_to(ROOT)).replace("\\", "/"): sha256(p) for p in sources}
    # Numerical checks are fixed before either diagnostic is evaluated. They do
    # not replace stage 3 spatial, boundary or coupled temporal admission gates.
    tolerances = dict(graph_relative=1e-11, circuit_state_scaled=1e-8,
                      circuit_instantaneous_power_relative=1e-12)
    (out / "manifest.json").write_text(json.dumps(dict(
        schema="pysnspd.stage3.electrical-prescribed-diagnostic.v1",
        started_at_utc=datetime.now(timezone.utc).isoformat(), source_sha256=hashes,
        parameters=asdict(circuit), tolerances=tolerances,
        scope="Graph algebra and a prescribed resistor driving the external circuit; no detector evolution.",
        Lk_ext_scope="Explicit diagnostic input only; resolved/exterior inductance partition is not identified.",
        environment=dict(python=platform.python_version(), numpy=np.__version__,
                         scipy=scipy.__version__, matplotlib=matplotlib.__version__),
    ), indent=2, allow_nan=False) + "\n", encoding="utf-8")

    # Uniform open wire with a known, fixed superconducting current. Conductivity
    # is an illustrative input, not material data inferred from a catalogue.
    faces = 12
    length = 360e-9
    area = 120e-9 * 7e-9
    sigma = 2.5e5
    imposed_current = 30e-6
    supercurrent = 18e-6
    edges = np.column_stack((np.arange(faces), np.arange(1, faces + 1)))
    conductance = np.full(faces, sigma * area / (length / faces))
    superconducting = np.full(faces, supercurrent)
    injection = np.zeros(faces + 1)
    injection[0], injection[-1] = imposed_current, -imposed_current
    port = solve_potential(edges, conductance, superconducting, injection,
                           left_node=0, right_node=faces)
    gauge = solve_potential(edges, conductance, superconducting, injection,
        left_node=0, right_node=faces, reference_node=faces, reference_potential_V=.125)
    x = np.linspace(0., length, faces + 1)
    exact_phi = (imposed_current - supercurrent) * (length - x) / (sigma * area)
    voltage_scale = abs(exact_phi[0])
    power_scale = max(abs(port.port_power_W), abs(port.normal_joule_W), abs(port.super_power_W))
    graph_metrics = dict(
        potential_error_relative=float(np.max(abs(port.phi_V - exact_phi)) / voltage_scale),
        current_residual_relative=float(np.max(abs(port.current_residual_A)) / imposed_current),
        total_current_error_relative=float(np.max(abs(port.I_total_A - imposed_current)) / imposed_current),
        gauge_potential_error_relative=float(np.max(abs(gauge.phi_V - .125 - port.phi_V)) / voltage_scale),
        gauge_current_error_relative=float(np.max(abs(gauge.I_total_A - port.I_total_A)) / imposed_current),
        gauge_voltage_error_relative=float(abs(gauge.Vdev_V - port.Vdev_V) / voltage_scale),
        gauge_power_error_relative=float(abs(gauge.port_power_W - port.port_power_W) / power_scale),
        port_face_power_residual_relative=float(abs(port.power_residual_W) / power_scale),
        power_decomposition_residual_relative=float(abs(port.decomposition_residual_W) / power_scale),
    )
    graph_passed = all(v <= tolerances["graph_relative"] for v in graph_metrics.values())

    # A resistor is switched from 0 to 1000 ohm at t=0. This is a prescribed
    # circuit input, not a prediction for a hotspot or detector resistance.
    resistance = 1000.
    initial = circuit.stationary_for_prescribed_voltage(0.)
    stationary = circuit.stationary_for_prescribed_resistance(resistance)
    Rb, Lb = circuit.R_bias_ohm, circuit.L_bias_H
    Rl, C, Le = circuit.R_load_ohm, circuit.C_couple_F, circuit.Lk_ext_H
    matrix = np.array([[-(Rb + Rl) / Lb, Rl / Lb, -1 / Lb],
                       [Rl / Le, -(Rl + resistance) / Le, 1 / Le],
                       [1 / C, -1 / C, 0.]])
    scale = np.array([initial[0], initial[0], circuit.V_bias_V])
    matrix_scaled = matrix * scale[None, :] / scale[:, None]
    times = np.unique(np.r_[0., np.geomspace(1e-13, 1e-9, 80), np.linspace(0., 300e-9, 501)])
    exact = np.array([stationary + scale * (expm(matrix_scaled * t) @ ((initial - stationary) / scale))
                      for t in times])
    measured = solve_ivp(lambda t, state: circuit.rhs(state * scale, resistance * state[1] * scale[1]) / scale,
        (times[0], times[-1]), initial / scale, t_eval=times, method="Radau",
        jac=matrix_scaled, rtol=2e-10, atol=1e-12)
    if not measured.success or measured.y.shape != (3, len(times)):
        raise RuntimeError("Prescribed-input circuit solve failed: " + measured.message)
    states = measured.y.T * scale
    exact_vout = Rl * (exact[:, 0] - exact[:, 1])
    measured_vout = Rl * (states[:, 0] - states[:, 1])
    power = [circuit.power_balance(y, resistance * y[1]) for y in states]
    circuit_power_scale = max(max(abs(p.source_power_W), abs(p.device_power_W),
                                  abs(p.stored_energy_rate_W)) for p in power)
    circuit_metrics = dict(
        maximum_scaled_state_error=float(np.max(abs(states - exact) / scale)),
        maximum_state_error_SI=np.max(abs(states - exact), axis=0).tolist(),
        state_units=["A", "A", "V"],
        scaling_SI=scale.tolist(),
        maximum_Vout_error_V=float(np.max(abs(measured_vout - exact_vout))),
        instantaneous_power_residual_W=float(max(abs(p.residual_W) for p in power)),
        instantaneous_power_residual_relative=float(max(abs(p.residual_W) for p in power) / circuit_power_scale),
    )
    circuit_passed = (circuit_metrics["maximum_scaled_state_error"] <= tolerances["circuit_state_scaled"]
                     and circuit_metrics["instantaneous_power_residual_relative"] <= tolerances["circuit_instantaneous_power_relative"])
    arrays = out / "electrical_arrays.npz"
    np.savez_compressed(arrays, node_position_m=x, phi_V=port.phi_V, exact_phi_V=exact_phi,
        gauge_phi_V=gauge.phi_V, I_total_A=port.I_total_A, I_normal_A=port.I_normal_A,
        I_super_A=superconducting, time_s=times, circuit_exact=exact, circuit_numerical=states,
        Vout_exact_V=exact_vout, Vout_numerical_V=measured_vout,
        CM8_residual_W=np.array([p.residual_W for p in power]))

    plt.rcParams.update({"font.size": 10, "axes.titlesize": 12, "axes.labelsize": 11,
                         "pdf.fonttype": 42, "ps.fonttype": 42})
    fig, axes = plt.subplots(2, 2, figsize=(11.5, 8.))
    fig.subplots_adjust(left=.09, right=.98, top=.86, bottom=.14, hspace=.48, wspace=.32)
    fig.suptitle("Puerto eléctrico y circuito de lectura", fontsize=17, y=.97)
    fig.text(.5, .915, "Pruebas analíticas y entrada resistiva prescrita · sin evolución del detector",
             ha="center", fontsize=11, color="#444444")
    blue, orange, green = "#0072B2", "#D55E00", "#009E73"
    ax = axes[0, 0]
    ax.plot(x * 1e9, exact_phi * 1e3, color="#444444", lw=2, label="Referencia analítica")
    ax.plot(x * 1e9, port.phi_V * 1e3, "o", mfc="white", mec=blue, mew=1.5, ms=5,
            label="Solución del grafo")
    ax.set(title="a  Caída pasiva en el puerto abierto", xlabel="Posición (nm)", ylabel="Potencial (mV)")
    ax.legend(frameon=False)
    ax.text(.98, .52, f"$V_{{dev}}$ = {port.Vdev_V * 1e3:.3f} mV\n"
            f"Error relativo: {graph_metrics['potential_error_relative']:.1e}",
            transform=ax.transAxes, ha="right", va="center", fontsize=9)
    ax = axes[0, 1]
    vals = [supercurrent * 1e6, np.mean(port.I_normal_A) * 1e6, np.mean(port.I_total_A) * 1e6]
    bars = ax.bar(["Superconductora\nprescrita", "Normal\ncalculada", "Total\ninyectada"], vals,
                  color=[green, orange, blue], width=.6)
    for bar, value in zip(bars, vals):
        ax.text(bar.get_x() + bar.get_width() / 2, value + .8, f"{value:g}", ha="center")
    ax.set(title="b  Una sola corriente total por cara", ylabel="Corriente (µA)", ylim=(0, 36))
    ax.text(.98, .95, "$I_{tot}=I_{sc}+I_n$", ha="right", va="top", transform=ax.transAxes)
    ax = axes[1, 0]
    ax.plot(times * 1e9, exact_vout * 1e3, color="#333333", lw=2, label="Exponencial matricial")
    markers = np.unique(np.r_[np.searchsorted(times, np.linspace(0, times[-1], 25)), np.argmax(exact_vout)])
    ax.plot(times[markers] * 1e9, measured_vout[markers] * 1e3, "o", ms=4, mfc="white", mec=orange,
            label="Integración Radau")
    ax.set(title="c  Pulso del circuito ante $R_{dev}$ prescrita", xlabel="Tiempo (ns)", ylabel="$V_{out}$ (mV)")
    ax.legend(frameon=False)
    zoom = ax.inset_axes([.56, .42, .38, .28])
    early = times <= .2e-9
    zoom.plot(times[early] * 1e12, exact_vout[early] * 1e3, color="#333333", lw=1.)
    zoom.set(title="Inicio continuo", xlabel="Tiempo (ps)", xlim=(0., 200.))
    zoom.tick_params(labelsize=7)
    zoom.xaxis.label.set_size(7)
    zoom.title.set_size(8)
    zoom.spines[["top", "right"]].set_visible(False)
    ax = axes[1, 1]
    ax.plot(times * 1e9, (measured_vout - exact_vout) * 1e15, color=blue, lw=1.2)
    ax.axhline(0., color="#999999", lw=.7)
    ax.set(title="d  Diferencia entre los dos métodos", xlabel="Tiempo (ns)",
           ylabel="$V_{out}^{num}-V_{out}^{exacta}$ (fV)")
    ax.text(.98, .95, f"Máximo: {circuit_metrics['maximum_Vout_error_V'] * 1e15:.2g} fV\n"
            f"Residuo CM.8: {circuit_metrics['instantaneous_power_residual_W']:.2g} W",
            transform=ax.transAxes, ha="right", va="top", fontsize=9)
    for ax in axes.ravel():
        ax.spines[["top", "right"]].set_visible(False)
        ax.grid(axis="y", alpha=.18)
        ax.set_axisbelow(True)
    fig.text(.09, .065, f"Circuito histórico: 0.300 V; 10 kΩ; 1 µH; 50 Ω; 100 pF.  "
             f"$L_{{k,ext}}$ = {lk_ext_h * 1e9:g} nH: entrada explícita del ensayo.", fontsize=9)
    fig.text(.09, .035, "Rdev: salto de 0 a 1000 Ω impuesto. No identifica una resistencia del detector ni la partición de inductancia.", fontsize=9)
    png, pdf = out / "electrical_diagnostics.png", out / "electrical_diagnostics.pdf"
    fig.savefig(png, dpi=180)
    fig.savefig(pdf)
    plt.close(fig)
    result = dict(
        schema="pysnspd.stage3.electrical-prescribed-diagnostic-result.v1",
        status="PASS_PRESCRIBED_ELECTRICAL_CHECKS" if graph_passed and circuit_passed else "FAIL_PRESCRIBED_ELECTRICAL_CHECKS",
        scope="Electrical graph algebra and circuit with prescribed resistance; no field, occupation or detector transient.",
        detector_simulated=False, coupled_admission=False, production_promotion=False,
        source_sha256=hashes, manifest_sha256=sha256(out / "manifest.json"), tolerances=tolerances,
        graph=dict(passed=graph_passed, nodes=faces + 1, faces=faces,
            geometry=dict(length_m=length, cross_section_m2=area),
            conductivity_S_m=sigma, conductivity_scope="Illustrative explicitly selected input, not a material characterization.",
            superconducting_current_A=supercurrent, injected_current_A=imposed_current,
            voltage_V=port.Vdev_V, reference_shift_V=.125, metrics=graph_metrics,
            power_W=dict(port=port.port_power_W, normal_joule=port.normal_joule_W,
                         superconducting=port.super_power_W, face=port.face_power_W)),
        circuit=dict(passed=circuit_passed, parameters=asdict(circuit),
            inductance_scope="Explicit exterior test value; no identification of the resolved/exterior partition.",
            prescribed_resistance_ohm=resistance, resistance_switch_time_s=0.,
            initial_state=initial.tolist(), final_stationary_oracle=stationary.tolist(),
            state_order=["Ib_A", "Is_A", "vc_V"], duration_s=float(times[-1]), samples=len(times),
            linear_matrix_SI=matrix.tolist(), exact_reference="y(t)=y* + exp(A*t)*(y0-y*)",
            numerical_method="scipy.solve_ivp Radau on scaled circuit states with analytic Jacobian",
            solver=dict(nfev=measured.nfev, njev=measured.njev, nlu=measured.nlu, rtol=2e-10, atol_scaled=1e-12),
            metrics=circuit_metrics,
            sampled_Vout_peak_V=float(max(exact_vout)),
            sampled_Vout_peak_time_s=float(times[np.argmax(exact_vout)]),
            latency_or_detector_pulse_claim=False,
            power_scope="Instantaneous CM.8 identity at saved states; not an integrated coupled energy ledger."),
        artifacts_sha256={p.name: sha256(p) for p in (arrays, png, pdf)},
        runtime_seconds=time.perf_counter() - started,
    )
    (out / "electrical_diagnostics.json").write_text(json.dumps(result, indent=2, ensure_ascii=False, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps(dict(status=result["status"], runtime_seconds=result["runtime_seconds"],
                         graph_voltage_V=port.Vdev_V, maximum_Vout_error_V=circuit_metrics["maximum_Vout_error_V"],
                         output_dir=str(out)), ensure_ascii=False), flush=True)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--lk-ext-h", type=float, required=True,
                        help="Explicit external inductance in H; not the historical total inductance")
    args = parser.parse_args()
    result = run(args.output_dir, args.lk_ext_h)
    if not result["status"].startswith("PASS_"):
        raise SystemExit(2)


if __name__ == "__main__":
    main()
