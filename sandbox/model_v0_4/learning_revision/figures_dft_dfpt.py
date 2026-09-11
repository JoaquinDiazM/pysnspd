"""Synthetic teaching figures for E09/E10, and checks of solved examples only.

No DFT/DFPT solver is run. The electronic DOS uses two invented cosine bands;
the BCS DOS is analytic; the phonon/interaction spectra use finite box profiles.
Output is restricted to the two figures and a reproducibility JSON beside them.
"""
from pathlib import Path
import json
import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch

ROOT = Path(__file__).resolve().parents[3]
FIG = ROOT / "docs/modelo_v0_4/figuras"
FIG.mkdir(parents=True, exist_ok=True)
BLUE = "#236482"
GREEN = "#32866e"
ORANGE = "#bf7038"
INK = "#22333f"
GRAY = "#6e7780"
plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 10.5,
    "axes.spines.top": False, "axes.spines.right": False,
    "axes.labelcolor": INK, "text.color": INK,
    "axes.titlesize": 11.5, "axes.titleweight": "bold",
    "savefig.dpi": 220,
})


def box(ax, x, y, w, h, label, color=BLUE, size=10):
    ax.add_patch(FancyBboxPatch(
        (x, y), w, h, boxstyle="round,pad=.015", linewidth=1.25,
        facecolor=color + "12", edgecolor=color,
    ))
    ax.text(x + w / 2, y + h / 2, label, ha="center", va="center",
            fontsize=size)


def arrow(ax, start, end, color=INK, connectionstyle="arc3,rad=0"):
    ax.annotate("", xy=end, xytext=start, arrowprops={
        "arrowstyle": "->", "lw": 1.3, "color": color,
        "connectionstyle": connectionstyle,
    })


def save(fig, name):
    for ext in ("png", "pdf"):
        fig.savefig(FIG / f"{name}.{ext}", facecolor="white", bbox_inches="tight")
    plt.close(fig)


# Two spin-degenerate invented bands, not a Kohn-Sham calculation.
xi = np.linspace(-2.3, 2.3, 1801)
k = np.linspace(-np.pi, np.pi, 1800, endpoint=False)
bands = [-0.65 - 1.10 * np.cos(k), 0.75 - 0.95 * np.cos(k)]
sigma = 0.055  # eV, solely for finite plotting width.
dos = np.zeros_like(xi)
for band in bands:
    dos += 2 * np.exp(-0.5 * ((xi[:, None] - band[None, :]) / sigma) ** 2).mean(axis=1)
dos /= sigma * np.sqrt(2 * np.pi)
dos_area = np.trapezoid(dos, xi)
assert np.isclose(dos_area, 4.0, atol=1e-10)

fig = plt.figure(figsize=(10.5, 6.3), layout="constrained")
gs = fig.add_gridspec(2, 2, height_ratios=[1, 2.4])
ax = fig.add_subplot(gs[0, :])
ax.set(xlim=(0, 12), ylim=(-0.1, 2.25))
ax.axis("off")
ax.set_title("a  Autoconsistencia: cerrar el círculo entre densidad y modos", loc="left")
box(ax, .12, .75, 2.5, .88, "Densidad propuesta\n" + r"$n_e(\mathbf{r})$", BLUE)
box(ax, 3.15, .75, 2.5, .88, "Potencial efectivo\n" + r"$v_{\mathrm{eff}}[n_e]$", BLUE)
box(ax, 6.2, .75, 2.5, .88, "Orbitales y energías\n" + r"$\phi_i,\ \varepsilon_i$", GREEN)
box(ax, 9.25, .75, 2.6, .88, "Densidad reconstruida\n" + r"$\sum_i f_i|\phi_i|^2$", GREEN)
for x in (2.7, 5.75, 8.8):
    arrow(ax, (x, 1.19), (x + .34, 1.19))
ax.plot([10.55, 10.55, 1.37], [.66, .29, .29], color=ORANGE, lw=1.3)
arrow(ax, (1.37, .29), (1.37, .66), ORANGE)
ax.text(5.95, -.02, "Repetir hasta que ambas densidades coincidan", ha="center", color=ORANGE)

ax = fig.add_subplot(gs[1, 0])
ax.plot(xi, dos, color=BLUE, lw=2)
ax.fill_between(xi, 0, dos, where=xi <= 0, color=BLUE, alpha=.2)
ax.axvline(0, color=GRAY, lw=1, linestyle="--")
ax.set(xlim=(-2.2, 2.2), ylim=(0, 4.9), xlabel=r"Energía normal $\xi=\varepsilon-\mu$ (eV)",
       ylabel="DOS (modos / eV / celda)", title="b  DOS normal: agrupar modos por energía")
ax.text(.03, .95, "Dos bandas de juguete\nÁrea total: 4 modos/celda", transform=ax.transAxes,
        va="top", fontsize=9.5)
ax.text(.035, .52, "Sombreado:\nocupados a T = 0", transform=ax.transAxes, fontsize=9,
        bbox={"facecolor": "white", "edgecolor": "none", "alpha": .85})

ax = fig.add_subplot(gs[1, 1])
delta = 1.0  # meV, synthetic choice.
energy = delta + np.geomspace(.0001, 4, 1801)
rho = energy / np.sqrt(energy ** 2 - delta ** 2)
ax.plot(energy, rho, color=GREEN, lw=2)
ax.plot([0, delta], [0, 0], color=GREEN, lw=2)
ax.axhline(1, color=GRAY, linestyle="--", lw=1)
ax.axvspan(0, delta, color=GREEN, alpha=.09)
ax.axvline(delta, color=GREEN, linestyle=":", lw=1)
ax.scatter([2], [2 / np.sqrt(3)], color=ORANGE, zorder=5, s=27)
ax.annotate(r"$\rho(2\Delta)=1.155$", (2, 2 / np.sqrt(3)), (2.5, 2.1), fontsize=9.5,
            arrowprops={"arrowstyle": "->", "color": ORANGE})
ax.text(2.8, 1.06, "Referencia normal", fontsize=8.8, color=GRAY)
ax.text(.48, 2.2, "Sin\nexcitaciones", ha="center", fontsize=9, color=GREEN)
ax.text(1.65, 4.45, "El borde ideal\ndiverge", fontsize=9)
arrow(ax, (1.62, 4.35), (1.04, 5.0), GRAY)
ax.text(.97, -.1, r"$\Delta$", ha="center", va="top", transform=ax.get_xaxis_transform(),
        fontsize=11, color=GREEN)
ax.set(xlim=(0, 5), ylim=(-.05, 5.3), xlabel="Energía de excitación E (meV)",
       ylabel=r"DOS superconductora normalizada $\rho(E)$", title="c  Añadir emparejamiento: modelo BCS")
ax.text(.98, .96, r"$\Delta=1$ meV", ha="right", va="top", transform=ax.transAxes, fontsize=9.5)
fig.suptitle("Modelos ilustrativos: no son resultados ab initio ni mediciones", fontsize=11, color=GRAY)
save(fig, "E09_dft_dos")


def log_window_integral(cutoff, lower, upper):
    """Integral from zero to cutoff of 1/Omega over the positive window."""
    return np.log(np.clip(cutoff, lower, upper) / lower)


def cumulative(cutoff, weights):
    return 2 * (weights[0] / 4 * log_window_integral(cutoff, 8, 12)
                + weights[1] / 4 * log_window_integral(cutoff, 28, 32))


edges = np.array([0, 8, 12, 28, 32, 40.])
phonon = np.array([0, .25, 0, .25, 0.])
a2f_a = np.array([0, .5, 0, .5, 0.])
a2f_b = np.array([0, .75, 0, .25, 0.])
omega = np.linspace(0, 40, 2001)
lambda_a = float(cumulative(40, [2, 2]))
lambda_b = float(cumulative(40, [3, 1]))
assert np.isclose(np.sum(phonon * np.diff(edges)), 2)
assert np.isclose(np.sum(a2f_a * np.diff(edges)), 4)
assert np.isclose(np.sum(a2f_b * np.diff(edges)), 4)

fig, axs = plt.subplots(2, 2, figsize=(10.5, 7.0), layout="constrained")
ax = axs[0, 0]
ax.axis("off")
ax.set(xlim=(0, 10), ylim=(0, 6))
ax.set_title("a  Respuesta electrónica a un desplazamiento", loc="left")
box(ax, 1.1, 4.05, 7.8, 1.07, "Mover átomos → cambia la densidad electrónica", BLUE, 9.8)
arrow(ax, (3.1, 3.98), (3.1, 3.33))
arrow(ax, (7, 3.98), (7, 3.33))
box(ax, .2, 2.05, 5.15, 1.2, "Fuerzas y rigideces\n→ modos y frecuencias", BLUE)
box(ax, 5.75, 2.05, 4.05, 1.2, "Potencial perturbado\n→ interacción g", ORANGE)
arrow(ax, (2.77, 1.98), (2.77, 1.33))
arrow(ax, (7.77, 1.98), (7.77, 1.33))
box(ax, .2, .25, 5.15, 1.0, r"$F(\Omega)$: conteo de modos", BLUE)
box(ax, 5.75, .25, 4.05, 1.0, r"$\alpha^2F(\Omega)$: interacción", ORANGE)

ax = axs[0, 1]
ax.stairs(phonon, edges, color=INK, lw=2, label="Misma DOS en A y B")
ax.stairs(phonon, edges, color=GRAY, fill=True, alpha=.12)
ax.set(xlim=(0, 40), ylim=(0, .36), xlabel=r"Energía fonónica $\Omega$ (meV)",
       ylabel="F (modos / meV)", title="b  Dos grupos de modos disponibles")
ax.text(10, .28, "Área 1", ha="center", fontsize=10)
ax.text(30, .28, "Área 1", ha="center", fontsize=10)
ax.legend(loc="lower center", fontsize=9, frameon=False)

ax = axs[1, 0]
ax.stairs(a2f_a, edges, color=BLUE, lw=2.2, label="A: pesos 2 y 2 meV")
ax.stairs(a2f_b, edges, color=ORANGE, lw=2, linestyle="--", label="B: pesos 3 y 1 meV")
ax.set(xlim=(0, 40), ylim=(0, 1.01), xlabel=r"Energía fonónica $\Omega$ (meV)",
       ylabel=r"$\alpha^2F(\Omega)$ (adimensional)", title="c  Misma área total; distintos pesos")
ax.legend(loc="upper right", fontsize=9, frameon=False)
ax.text(.51, .32, "Área total:\n4 meV en ambos", transform=ax.transAxes, ha="center", fontsize=9.2)

ax = axs[1, 1]
ax.plot(omega, cumulative(omega, [2, 2]), color=BLUE, lw=2.2, label=f"A: λ = {lambda_a:.4f}")
ax.plot(omega, cumulative(omega, [3, 1]), color=ORANGE, lw=2.0, linestyle="--", label=f"B: λ = {lambda_b:.4f}")
ax.set(xlim=(0, 40), ylim=(0, .84), xlabel=r"Límite de integración $\Omega_c$ (meV)",
       ylabel=r"$\lambda(<\Omega_c)$ (adimensional)", title="d  Acumular la interacción con peso 1/Ω")
ax.legend(loc="lower right", fontsize=9.5, frameon=False)
ax.text(.035, .93, r"$\lambda(<\Omega_c)=2\int_0^{\Omega_c}\alpha^2F(\Omega)\,d\Omega/\Omega$",
        transform=ax.transAxes, va="top", fontsize=10.2)
fig.suptitle("Espectros sintéticos de dos intervalos; las curvas no representan ocupaciones", fontsize=11, color=GRAY)
save(fig, "E10_dfpt_a2f")

# Independent numerical quadrature inside each smooth positive window.
numerical_lambda = []
for weights in ([2, 2], [3, 1]):
    val = 0
    for weight, lo, hi in zip(weights, (8, 28), (12, 32)):
        q = np.linspace(lo, hi, 10001)
        val += 2 * np.trapezoid(weight / 4 / q, q)
    numerical_lambda.append(float(val))
assert np.allclose(numerical_lambda, [lambda_a, lambda_b], rtol=1e-9)

def emission_integral(weights, energy=16, gap=1):
    total = 0.
    for weight, lo, hi in zip(weights, (8, 28), (12, 32)):
        upper = min(hi, energy - gap)
        if upper <= lo:
            continue
        q = np.linspace(lo, upper, 10001)
        final = energy - q
        spectral = final / np.sqrt(final ** 2 - gap ** 2)
        coherence = 1 - gap ** 2 / (energy * final)
        total += np.trapezoid(weight / 4 * spectral * coherence, q)
    return float(total)

emission_ratio = emission_integral([3, 1]) / emission_integral([2, 2])
assert np.isclose(emission_ratio, 1.5, rtol=1e-12)
p = np.array([.5, .75, 1])
energies = 2 * p ** 2 - 3 * p + 2
assert np.allclose(energies, [1, .875, 1])

checks = {
    "scope": "Solved teaching examples; no DFT/DFPT computation or measured material data",
    "dft_discrete_example": {"p": p.tolist(), "energy_eV": energies.tolist(), "minimum_p": .75},
    "synthetic_normal_dos_integral_modes_per_cell": float(dos_area),
    "BCS_rho_at_2Delta": float(2 / np.sqrt(3)),
    "phonon_example_integral_modes": float(np.sum(phonon * np.diff(edges))),
    "alpha2F_integral_meV_A_B": [float(np.sum(a2f_a * np.diff(edges))), float(np.sum(a2f_b * np.diff(edges)))],
    "lambda_exact_A_B": [lambda_a, lambda_b],
    "lambda_numeric_A_B": numerical_lambda,
    "lambda_after_first_interval_A_B": [float(cumulative(12, [2, 2])), float(cumulative(12, [3, 1]))],
    "emission_loss_integral_ratio_B_A_E16_Delta1": emission_ratio,
    "source_equation_map": {"Simon_version": "arXiv:2501.13791v3", "emission": "(1a), second term; dilute zero-phonon limit of (6b)", "lambda": "definition following (2)"},
}
check_path = FIG / "E10_dfpt_a2f.checks.json"
check_path.write_text(json.dumps(checks, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
print(json.dumps(checks, indent=2, ensure_ascii=False))
