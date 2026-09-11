"""Light, deterministic checks for document B revision 0.2.

No material data, solver imports or production transient are used.  Energies
are dimensionless (gap = k_B = 1); stored measures specify their normalization.
Run from any directory; output is resolved relative to this script.
"""
from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "docs" / "modelo_v0_2"
FIG = OUT / "figuras"
DATA = OUT / "verificaciones"
FIG.mkdir(parents=True, exist_ok=True)
DATA.mkdir(parents=True, exist_ok=True)

BLUE, RED, TEAL, GRAY = "#2463A6", "#CF6745", "#218A83", "#657181"
plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 11, "axes.titlesize": 12,
    "axes.labelsize": 11, "axes.spines.top": False, "axes.spines.right": False,
    "axes.grid": True, "grid.alpha": .16, "figure.facecolor": "white",
    "savefig.facecolor": "white", "pdf.fonttype": 42,
})


def save(fig, name):
    fig.savefig(FIG / (name + ".png"), dpi=185, bbox_inches="tight")
    fig.savefig(FIG / (name + ".pdf"), bbox_inches="tight")
    plt.close(fig)


def csv_write(name, fields, rows):
    with (DATA / name).open("w", encoding="utf-8", newline="") as out:
        writer = csv.writer(out)
        writer.writerow(fields)
        writer.writerows(rows)


def fd(energy, temperature):
    return 1. / (1. + np.exp(np.clip(np.asarray(energy) / temperature, 0, 700)))


def bose(energy, temperature):
    if temperature == 0:
        return np.zeros_like(np.asarray(energy), dtype=float)
    return 1. / np.expm1(np.clip(np.asarray(energy) / temperature, 1e-300, 700))


def gauss_interval(n, lower, upper):
    nodes, weights = np.polynomial.legendre.leggauss(n)
    return lower + (upper - lower) * (nodes + 1) / 2, weights * (upper - lower) / 2


def reaction_checks():
    rng = np.random.default_rng(290802)
    e = rng.uniform(.04, 8, 3000)
    omega = rng.uniform(.04, 8, 3000)
    f_low, f_high = rng.uniform(.0, .49, (2, len(e)))
    nph = rng.uniform(.0, .9, len(e))
    bs = f_high * (1 - f_low) * (nph + 1) - f_low * (1 - f_high) * nph
    # An arbitrary positive quadrature/rate weight tests the stoichiometry.
    rate_s = rng.uniform(.1, 2, len(e)) * bs
    e_pair = omega * rng.uniform(.02, .98, len(e))
    fa, fb = rng.uniform(.0, .49, (2, len(e)))
    br = fa * fb * (nph + 1) - (1 - fa) * (1 - fb) * nph
    rate_r = .5 * rng.uniform(.1, 2, len(e)) * br
    de_s = (e - (e + omega)) * rate_s
    de_r = -(e_pair + (omega - e_pair)) * rate_r
    dp_s, dp_r = omega * rate_s, omega * rate_r
    scale = np.sum(abs(de_s)) + np.sum(abs(de_r))
    residual = float(abs(np.sum(de_s + dp_s) + np.sum(de_r + dp_r)) / scale)
    max_event = float(max(np.max(abs(de_s + dp_s)), np.max(abs(de_r + dp_r))))
    # Detailed balance and thermal reduction are tested independently of rates.
    temp = .7
    f0, f1 = fd(e, temp), fd(e + omega, temp)
    n_eq = bose(omega, temp)
    a, b = fd(e_pair, temp), fd(omega - e_pair, temp)
    terms_s = (f1 * (1 - f0) * (n_eq + 1), f0 * (1 - f1) * n_eq)
    terms_r = (a * b * (n_eq + 1), (1 - a) * (1 - b) * n_eq)
    eq_s = float(np.max(abs(terms_s[0] - terms_s[1])))
    eq_r = float(np.max(abs(terms_r[0] - terms_r[1])))
    direct_s = f1 * (1 - f0) * (nph + 1) - f0 * (1 - f1) * nph
    direct_r = a * b * (nph + 1) - (1 - a) * (1 - b) * nph
    reduced_s = (f0 - f1) * (n_eq - nph)
    reduced_r = (1 - a - b) * (n_eq - nph)
    reduction_error = float(max(np.max(abs(direct_s - reduced_s)), np.max(abs(direct_r - reduced_r))))
    fig, axes = plt.subplots(2, 1, figsize=(10, 5.7))
    for ax in axes:
        ax.set_axis_off()
        ax.set_xlim(0, 10)
        ax.set_ylim(0, 1.6)
    cases = [
        (axes[0], "Dispersión: una cuasipartícula cambia de nivel", "QP: 3", "QP: 1", "Fonón: 2", r"$\Delta u_e=-2,\quad \Delta u_{ph}=+2$"),
        (axes[1], "Recombinación: dos cuasipartículas desaparecen", "QP: 1.4 + QP: 1.6", "Condensado*", "Fonón: 3", r"$\Delta u_e=-3,\quad \Delta u_{ph}=+3$"),
    ]
    for ax, title, left, right, ph, balance in cases:
        ax.text(0, 1.45, title, fontweight="bold", color="#23344D")
        ax.text(.35, .88, left, color=BLUE, bbox=dict(boxstyle="round,pad=.4", fc="#EAF1FA", ec="none"))
        ax.annotate("", xy=(5.35, .91), xytext=(3.15, .91), arrowprops=dict(arrowstyle="->", color=GRAY, lw=2))
        ax.text(5.5, .91, right, color=BLUE, va="center")
        ax.text(8, .91, "+ " + ph, color=RED, va="center")
        ax.text(.35, .16, balance + r"$\quad\Rightarrow\quad\Delta(u_e+u_{ph})=0$", fontsize=13)
    fig.text(.13, .025, "FIGURA PEDAGÓGICA · Energías de ejemplo. *Se mantiene fijo el espectro durante la colisión.\nLa reacción inversa invierte todas las flechas y todos los signos.", color=GRAY, fontsize=9)
    fig.subplots_adjust(hspace=.24, bottom=.16)
    save(fig, "B_01_reacciones_balance")
    csv_write("B_reacciones.csv", ["reaction", "electronic_energy_change_arbitrary", "phonon_energy_change_arbitrary"],
              [("scattering", float(np.sum(de_s)), float(np.sum(dp_s))), ("recombination", float(np.sum(de_r)), float(np.sum(dp_r)))])
    assert residual < 1e-14
    assert reduction_error < 1e-13
    assert max(eq_s, eq_r) < 1e-13
    return {"sample_count_each_reaction": 3000, "rng_seed": 290802, "collision_balance_relative_residual": residual,
            "max_event_absolute_residual": max_event, "equilibrium_BS_max_absolute": eq_s,
            "equilibrium_BR_max_absolute": eq_r, "thermal_reduction_max_absolute_error": reduction_error,
            "scope": "arbitrary bounded populations and positive weights; algebraic reaction test, no material kernel"}


def energy_checks():
    x, w = gauss_interval(1024, 0, 60)
    x2, w2 = gauss_interval(512, 0, 60)
    energy = np.sqrt(1 + x*x)
    energy2 = np.sqrt(1 + x2*x2)

    def u(temp, delta=1):
        en = np.sqrt(delta*delta + x*x)
        return float(np.sum(w * en * fd(en, temp)))

    def cv(temp):
        p = fd(energy, temp)
        return float(np.sum(w * energy**2 * p * (1-p)) / temp**2)

    def inverse(target):
        lo, hi = 1e-4, 2.
        for _ in range(80):
            mid = (lo+hi)/2
            if u(mid) < target:
                lo = mid
            else:
                hi = mid
        return (lo+hi)/2

    temps = np.linspace(.055, 1.5, 220)
    energies = np.array([u(t) for t in temps])
    caps = np.array([cv(t) for t in temps])
    coarse = np.array([np.sum(w2 * energy2 * fd(energy2, t)) for t in temps])
    convergence = float(np.max(abs(coarse - energies) / energies))
    recovered = np.array([inverse(v) for v in energies])
    inversion_error = float(np.max(abs(recovered - temps)))
    dt = 1e-5
    central_caps = np.array([(u(t+dt)-u(t-dt))/(2*dt) for t in temps])
    derivative_error = float(np.max(abs(central_caps-caps)/caps))
    target = .025
    te = inverse(target)
    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    ax = axes[0]
    ax.plot(temps, energies, color=BLUE, lw=2.3, label="BCS, amplitud fija")
    ax.plot(temps, np.pi**2/12*temps**2, color=GRAY, ls="--", label="Sin gap (referencia)")
    ax.scatter([te], [target], color=RED, zorder=4)
    ax.plot([0, te, te], [target, target, 0], ls=":", color=RED, lw=1.5)
    ax.annotate(r"Una energía $\to$ un $T_E$" + f"\n$T_E={te:.3f}$", xy=(te, target), xytext=(.53, .34),
                fontsize=10, arrowprops=dict(arrowstyle="->", color=RED), color=RED)
    ax.set(xlabel=r"Temperatura $k_BT/|\Delta|$", ylabel=r"$u_{qp}/(4N_0|\Delta|^2)$", xlim=(0, 1.5), title="El mapa depende del espectro")
    ax.legend(frameon=False, fontsize=9)
    ax = axes[1]
    ax.semilogy(temps, 1/caps, color=TEAL, lw=2.3)
    ax.set(xlabel=r"Temperatura $k_BT/|\Delta|$", ylabel=r"Sensibilidad absoluta $dT/du_{qp}$ (normalizada)", title="El extremo frío exige precisión")
    ax.text(.46, .82, "Pequeños errores de energía\npueden mover la temperatura\nsi la capacidad es pequeña.", transform=ax.transAxes, fontsize=10, color=GRAY)
    fig.text(.1, -.025, "RESULTADO CUANTITATIVO SIMPLIFICADO · BCS homogéneo, q = 0, amplitud fijada; sin datos de película.", fontsize=9, color=GRAY)
    fig.tight_layout()
    save(fig, "B_02_energia_temperatura")
    csv_write("B_energia_temperatura.csv", ["kBT_over_gap", "u_over_4N0gap2", "du_dT_normalized", "inverse_absolute_sensitivity", "recovered_temperature"],
              zip(temps, energies, caps, 1/caps, recovered))
    # Equal-energy packet comparison. All occupations lie between 0 and 1/2.
    p_low = np.exp(-.5*((energy-1.2)/.045)**2)
    p_high = np.exp(-.5*((energy-6)/.2)**2)
    p_low *= target / np.sum(w * energy * p_low)
    p_high *= target / np.sum(w * energy * p_high)
    p_thermal = fd(energy, te)
    dist = [p_low, p_high, p_thermal]
    names = ["Paquete cerca del gap", "Paquete de alta energía", "Fermi–Dirac equivalente"]
    colors = [BLUE, RED, TEAL]
    packet_u = [float(np.sum(w * energy * p)) for p in dist]
    force = [float(np.sum(w / energy * p)) for p in dist]
    moment_3 = [float(np.sum(w * energy**3 * p)) for p in dist]
    ratio = force[0]/force[1]
    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    e_plot = np.linspace(1, 8, 1800)
    norms = [float(p_low.max()/np.exp(-.5*((energy-1.2)/.045)**2).max()),
             float(p_high.max()/np.exp(-.5*((energy-6)/.2)**2).max())]
    y_plot = [norms[0]*np.exp(-.5*((e_plot-1.2)/.045)**2), norms[1]*np.exp(-.5*((e_plot-6)/.2)**2), fd(e_plot, te)]
    for name, color, pp in zip(names, colors, y_plot):
        axes[0].semilogy(e_plot, np.maximum(pp, 1e-12), color=color, lw=2, label=name)
    axes[0].set(xlabel=r"Energía $E/|\Delta|$", ylabel="Ocupación f(E)", ylim=(1e-5, .3), title="Tres repartos de la misma energía")
    axes[0].legend(frameon=False, fontsize=8.5)
    pos = np.arange(3)
    axes[1].bar(pos-.18, np.array(packet_u)/target, .34, color="#C8D4E3", label="Energía / objetivo")
    axes[1].bar(pos+.18, np.array(force)/force[2], .34, color=colors, label="Fuerza QP / fuerza térmica")
    axes[1].set_xticks(pos, ["Cerca del gap", "Alta energía", "Térmica"], fontsize=9)
    axes[1].set(ylabel="Momento relativo", ylim=(0,1.3), title=f"La fuerza de los paquetes difiere ×{ratio:.2f}")
    axes[1].legend(frameon=False, fontsize=8.5)
    fig.text(.1, -.025, r"RESULTADO CUANTITATIVO SIMPLIFICADO · $u_{qp}/(4N_0|\Delta|^2)=0.025$; fuerza: contribución QP, no fuerza total.", fontsize=9, color=GRAY)
    fig.tight_layout()
    save(fig, "B_03_misma_energia")
    csv_write("B_paquetes.csv", ["population", "u_over_4N0gap2", "Xqp_over_4N0gap", "third_energy_moment", "T_E", "max_occupation"],
              [(name, uu, ff, mm, te, float(pp.max())) for name, uu, ff, mm, pp in zip(names, packet_u, force, moment_3, dist)])
    csv_write("B_distribuciones.csv", ["E_over_gap", "f_near_gap", "f_high_energy", "f_FD"], zip(e_plot, *y_plot))
    # Instantaneous energy-transfer check at a deliberately nonequilibrium gap.
    # u_total/(4N0 Delta0^2) = d^2/4*(log(d)-1/2) + integral E p dx.
    d, t = .4, .08
    en = np.sqrt(d*d+x*x)
    pop = fd(en, t)
    xforce = .5*d*np.log(d) + float(np.sum(w*d/en*pop))
    d_dot = -xforce  # positive mobility = 1 in declared arbitrary time units
    q_delta = xforce*xforce
    ce = float(np.sum(w*en**2*pop*(1-pop))/t**2)
    p_temp = en/t**2 * pop*(1-pop)
    # At thermal fixed T, u_d = vacuum_d + integral[(d/E)p + E dp/dd].
    u_d = .5*d*np.log(d) + float(np.sum(w*(d/en*pop - d/t*pop*(1-pop))))
    x_T = float(np.sum(w*d/en*p_temp))
    t_dot = (q_delta + t*x_T*d_dot) / ce
    chain_residual = u_d*d_dot + ce*t_dot
    assert np.all(np.diff(energies) > 0)
    assert inversion_error < 1e-10
    assert convergence < 1e-8
    assert abs(chain_residual) < 1e-12
    assert max(float(pp.max()) for pp in dist) < .5
    return {"quadrature": {"variable": "x/Delta; E=sqrt(x^2+Delta^2)", "x_max": 60, "nodes_fine": 1024, "nodes_coarse": 512},
            "temperature_range_kBT_over_gap": [float(temps[0]), float(temps[-1])],
            "quadrature_refinement_max_relative_change": convergence,
            "inversion_max_absolute_error": inversion_error,
            "central_difference_capacity_max_relative_error": derivative_error,
            "min_capacity_in_test": float(caps.min()), "target_energy_normalized": target,
            "equivalent_temperature": te, "force_packet_ratio": ratio,
            "packet_energies": packet_u, "packet_force_moments": force,
            "isolated_instantaneous_check": {"delta_over_Delta0": d, "kBT_over_Delta0": t,
                "amplitude_force_normalized": xforce, "delta_dot_arbitrary_time": d_dot,
                "Q_delta_normalized": q_delta, "T_dot_arbitrary_time": t_dot,
                "energy_chain_rule_residual": chain_residual,
                "scope": "local derivative identity, not integrated dynamics or material prediction"}}


def normal_limit_checks():
    temp = .7
    omega_over_t = np.linspace(.04, 9, 160)
    omega = temp*omega_over_t
    e, w = gauss_interval(512, 0, 45)
    j_s = np.array([np.sum(w*(fd(e,temp)-fd(e+om,temp))) for om in omega])
    j_r = []
    for om in omega:
        ep, wp = gauss_interval(48, 0, om)
        j_r.append(float(np.sum(wp*(1-fd(ep,temp)-fd(om-ep,temp)))))
    j_r = np.array(j_r)
    identity_error = float(np.max(abs(2*j_s+j_r-omega)/omega))
    ph, wph = gauss_interval(512, 0, 90)
    ratios = np.linspace(0, 1.4, 100)
    powers = np.array([np.sum(wph*ph**4*(bose(ph,temp)-bose(ph,temp*rr))) for rr in ratios])
    scaled_power = powers/powers[0]
    law_error = float(np.max(abs(scaled_power-(1-ratios**5))))
    fig, axes = plt.subplots(1, 2, figsize=(10,4))
    axes[0].plot(omega_over_t, 2*j_s/omega, color=BLUE, label=r"$2J_S/\Omega$")
    axes[0].plot(omega_over_t, j_r/omega, color=RED, label=r"$J_R/\Omega$")
    axes[0].plot(omega_over_t, (2*j_s+j_r)/omega, color=TEAL, ls="--", label="Suma = 1")
    axes[0].set(xlabel=r"$\Omega/(k_BT_e)$", ylabel="Fracción del núcleo normal", ylim=(-.03,1.12), title="Dos canales forman el límite normal")
    axes[0].legend(frameon=False, fontsize=9)
    axes[1].plot(ratios, 1-ratios**5, color=GRAY, lw=3, label=r"Referencia $1-(T_{ph}/T_e)^5$")
    axes[1].plot(ratios[::5], scaled_power[::5], "o", color=BLUE, ms=4, label="Integral numérica")
    axes[1].axhline(0,color=GRAY,lw=.7)
    axes[1].axvline(1,color=GRAY,lw=.7,ls=":")
    axes[1].text(.32,.33,"P > 0: electrones ceden energía",transform=axes[1].transAxes,fontsize=9,color=TEAL)
    axes[1].text(.40,.08,"P < 0: electrones reciben energía",transform=axes[1].transAxes,fontsize=9,color=RED)
    axes[1].set(xlabel=r"$T_{ph}/T_e$", ylabel=r"$P(T_e,T_{ph})/P(T_e,0)$", title="El signo surge del balance detallado")
    axes[1].legend(frameon=False,fontsize=8.5)
    fig.text(.1,-.025,r"RESULTADO CUANTITATIVO SIMPLIFICADO · Estado normal, $\alpha^2F\propto\Omega^2$, espectro acústico ideal sin corte.",fontsize=9,color=GRAY)
    fig.tight_layout()
    save(fig,"B_04_limite_normal")
    csv_write("B_limite_normal.csv", ["Omega_over_kBT", "two_JS_over_Omega", "JR_over_Omega", "sum_over_Omega"], zip(omega_over_t,2*j_s/omega,j_r/omega,(2*j_s+j_r)/omega))
    csv_write("B_potencia_normal.csv", ["Tph_over_Te", "P_over_P_Te_0", "fifth_power_reference"], zip(ratios,scaled_power,1-ratios**5))
    assert identity_error < 1e-10
    assert law_error < 1e-9
    return {"normal_2JS_plus_JR_relative_error": identity_error, "normal_fifth_power_max_absolute_error": law_error,
            "scope": "normal-state constant DOS and ideal alpha2F proportional to Omega^2; not DFT spectrum"}


def main():
    result = {"document": "B", "revision": "0.2", "date": "2026-09-08", "normalization": "gap=k_B=1; u_qp unit 4N0 gap^2; X_QP unit 4N0 gap",
              "numpy_version": np.__version__, "matplotlib_version": matplotlib.__version__}
    result["reactions"] = reaction_checks()
    result["energy_and_populations"] = energy_checks()
    result["normal_limit"] = normal_limit_checks()
    (DATA/"B_verificaciones.json").write_text(json.dumps(result,indent=2,ensure_ascii=False),encoding="utf-8")
    print(json.dumps(result,indent=2,ensure_ascii=False))


if __name__ == "__main__":
    main()
