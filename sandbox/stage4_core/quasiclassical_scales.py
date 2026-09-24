"""Algebraic scale estimates, not a material fit or a physical transient."""
from pathlib import Path
import hashlib
import json
import math

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "docs/implementation/stage4/final_kwt_20260924/scales"


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    # Exact SI constants; hbar derived from h, energies reported in meV.
    h = 6.62607015e-34
    kb = 1.380649e-23
    e = 1.602176634e-19
    hbar = h / (2 * math.pi)
    tc, tb, diffusion = 8.65, 0.9, 5e-5
    delta = 1.764 * kb * tc
    tau_gap = hbar / delta * 1e12
    # Other NbN films (Lomakin Table I), NOT measured Korzh parameters.
    elastic_fs = [0.15, 0.7]
    ef_ev = [5.0, 6.9]
    vf_other = 5.3e5
    conditional_tau = 3 * diffusion / vf_other**2 * 1e15
    conditional_ell = 3 * diffusion / vf_other * 1e9
    records = {
        "schema": "pysnspd.stage4.algebraic_scale_assessment.v1",
        "physical_transients": 0,
        "production_parameters_changed": False,
        "scope": "Scale hierarchy only; no prediction of latency, jitter or charge-imbalance relaxation.",
        "sources": {
            "Korzh2020_SI_table1": "https://eprints.lancs.ac.uk/id/eprint/140252/3/Binder1.pdf",
            "Lomakin2023_tableI_appendixC": "https://arxiv.org/pdf/2207.05012",
            "Belzig1999_sections2_2_2_4": "https://arxiv.org/abs/cond-mat/9812297",
        },
        "Korzh_reference": {
            "Tc_K": tc, "Tb_K": tb, "D_m2_s": diffusion,
            "tau_ee_at_Tc_ps": 6.0, "tau_ep_at_Tc_ps": 24.7,
            "status": "Model-reference inputs; fitted times are not branch-imbalance time at Tb.",
        },
        "BCS_scale": {
            "gap_ratio_assumed": 1.764, "Delta0_meV": delta/e*1e3,
            "hbar_over_Delta0_ps": tau_gap,
            "hbar_over_kBTb_ps": hbar/(kb*tb)*1e12,
            "sqrt_hbar_D_over_Delta0_nm": math.sqrt(hbar*diffusion/delta)*1e9,
            "meaning": "Reference weak-coupling gap scale, not a measured gap or a relaxation time.",
        },
        "comparison_NbN_films": {
            "elastic_time_range_fs": elastic_fs,
            "Fermi_energy_range_eV": ef_ev,
            "hbar_over_EF_range_fs": [hbar/(ef_ev[1]*e)*1e15, hbar/(ef_ev[0]*e)*1e15],
            "kF_ell_range": [1.6, 6.3],
            "vF_fit_m_s": vf_other,
            "conditional_tau_el_fs_using_Korzh_D": conditional_tau,
            "conditional_ell_nm_using_Korzh_D": conditional_ell,
            "scope": "Different 2.5-nm NbN films; illustrative compatibility, not sample-specific calibration.",
        },
        "adiabatic_indicators": [
            {"variation_time_ps": t, "gap_fraction": f, "hbar_over_gap_time": tau_gap/f/t}
            for f in [1., .5, .1] for t in [1., 3., 10., 20.]
        ],
        "weak_amplitude_example": {
            "relative_gap_change": 0.01, "duration_ps": 1.0,
            "hbar_abs_dDelta_dt_over_Delta_squared": 0.01*tau_gap,
            "meaning": "A weak 1% change over 1 ps is not an order-one change over 1 ps.",
        },
        "diffusion_scales": [
            {"length_nm": length, "L_squared_over_D_ps": (length*1e-9)**2/diffusion*1e12}
            for length in [5., 10., 40., 80.]
        ],
        "definitions": {
            "variation_time": "Hypothetical relative-amplitude time |Delta|/|dDelta/dt|; not jitter, oscillation period, latency, timestep or final horizon. For an order-one change only, comparable to event duration.",
            "indicator": "hbar/(gap*t) is an order-of-magnitude spectral indicator, NOT an error bound. Broadenings, gradients, energy edges and moving phase also matter.",
            "diffusion": "L^2/D convention; neither mean first-passage time nor lowest-mode decay L^2/(pi^2 D).",
            "tau_Q": "Undetermined for the selected wire; never identified with tau_el, hbar/Delta, tau_ee or tau_ep.",
        },
        "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    }
    (OUT / "scales.json").write_text(json.dumps(records, indent=2, ensure_ascii=False)+"\n", encoding="utf8")
    plt.rcParams.update({"font.size": 10, "axes.titlesize": 11})
    fig, axes = plt.subplots(1, 2, figsize=(13, 5.5), gridspec_kw={"width_ratios": [1.15, 1]})
    ax = axes[0]
    spans = [records["comparison_NbN_films"]["hbar_over_EF_range_fs"], elastic_fs]
    for row, span in zip([4, 3], spans):
        ax.plot(np.array(span)*1e-3, [row, row], lw=7, color="#4c78a8", solid_capstyle="round")
    ax.scatter([tau_gap, 6., 24.7], [2, 1, 0], color="#d17c29", s=55, zorder=3)
    ax.set_yticks(range(5), ["e–fonón a Tc (referencia)", "e–e a Tc (referencia)",
                              "ħ/Δ₀ (cálculo BCS)", "Colisión elástica (otras películas)", "ħ/E_F (otras películas)"])
    ax.set_xscale("log"); ax.set_xlim(5e-5, 100)
    ax.set_xlabel("Escala temporal [ps]; eje logarítmico")
    ax.set_title("A. Escalas diferentes, procesos diferentes")
    ax.grid(axis="x", alpha=.25)
    ax.text(.02, -.22, "τ_Q, relajación electrón–hueco: no determinada.\nLos rangos azules no son mediciones de la cinta de Korzh.",
            transform=ax.transAxes, fontsize=9)
    ax = axes[1]
    times = np.geomspace(.5, 50, 200)
    for fraction, color in [(1., "#238b45"), (.5, "#df8f27"), (.1, "#a23e48")]:
        ax.loglog(times, tau_gap/fraction/times, label=f"|Δ|/Δ₀ = {fraction:g}", color=color)
    ax.axhline(1., color="gray", ls="--", lw=1)
    ax.set_xlabel("Tiempo |Δ| / |∂tΔ| supuesto [ps]")
    ax.set_ylabel("Indicador ħ / (|Δ| · tiempo), adimensional")
    ax.set_title("B. Suprimir el gap debilita la separación")
    ax.grid(which="major", alpha=.25); ax.legend(loc="upper right")
    ax.text(.02, -.22, "Un cambio del 1 % en 1 ps tiene tiempo relativo ≈100 ps.\nEl indicador no es una cota del error del modelo.",
            transform=ax.transAxes, fontsize=9)
    fig.suptitle("Separación de escalas: estimaciones algebraicas, sin transiente", fontsize=13)
    fig.subplots_adjust(left=.23, right=.98, bottom=.25, top=.86, wspace=.35)
    fig.savefig(OUT/"scale_hierarchy.png", dpi=190)
    plt.close(fig)
    print(json.dumps({"status": "ALGEBRAIC_ONLY", "Delta0_meV": delta/e*1e3, "hbar_over_Delta0_ps": tau_gap,
                      "conditional_tau_el_fs": conditional_tau, "physical_transients": 0}))


if __name__ == "__main__":
    main()
