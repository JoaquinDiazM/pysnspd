"""Four microscopic circuit insets and p-n versus N-S, original teaching schematics.

Coordinates are diagram coordinates, not atomistic simulation outputs.
Run from repository root with the repository's Matplotlib PYTHONPATH.
"""
from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Circle, FancyBboxPatch, Ellipse

ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / "docs/modelo_v0_4/figuras"
BLUE, ORANGE, GREEN, INK, GRAY = "#236482", "#b36428", "#287663", "#203744", "#aab8be"
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 12.5,
                     "text.color": INK, "axes.labelcolor": INK,
                     "savefig.dpi": 250, "pdf.fonttype": 42})

def arrow(ax, start, end, color=INK, lw=1.7, style="->"):
    ax.annotate("", xy=end, xytext=start,
                arrowprops={"arrowstyle": style, "color": color, "lw": lw})

def text(ax, x, y, s, **kw):
    ax.text(x, y, s, ha=kw.pop("ha", "center"), va=kw.pop("va", "center"), **kw)

def electron(ax, x, y, hole=False, size=52):
    ax.scatter([x], [y], s=size, color="white" if hole else BLUE,
               edgecolors=ORANGE if hole else BLUE, linewidths=1.6, zorder=5)

def save(fig, name):
    OUT.mkdir(parents=True, exist_ok=True)
    for ext in ("png", "pdf"):
        fig.savefig(OUT / f"{name}.{ext}", bbox_inches="tight", facecolor="white")
    plt.close(fig)

def circuit(ax, kind):
    ax.plot([2, 7, 7], [12, 12, 12], color=INK, lw=1.7)
    ax.plot([19, 24], [12, 12], color=INK, lw=1.7)
    ax.scatter([2, 24], [12, 12], s=20, facecolor="white", edgecolor=INK, zorder=8)
    if kind == "sc":
        x = np.linspace(7, 19, 300)
        ax.plot(x, 12 + 1.25*np.sin(2*np.pi*(x-7)/3), color=GREEN, lw=1.8)
        label = r"$L_k$"
    else:
        ax.add_patch(Rectangle((7, 10.1), 12, 3.8, fc="white", ec=INK, lw=1.7))
        label = {"metal": r"$R$", "intrinsic": "Si", "doped": "Si:n"}[kind]
    text(ax, 13, 7.8, label)
    arrow(ax, (6, 17), (20, 17))
    text(ax, 13, 19.4, r"$I$ convencional", fontsize=12)
    ax.plot([19, 31], [14, 20], "--", color=GRAY, lw=1)
    ax.plot([19, 31], [10, 6], "--", color=GRAY, lw=1)

def lattice(ax, sc=False):
    ax.add_patch(FancyBboxPatch((31, 6), 34, 14, boxstyle="round,pad=0.15",
                  fc="#edf5f2" if sc else "#f4f7f8", ec=GRAY, lw=1.2))
    for x in [35, 43, 51, 59]:
        for y in [9, 13, 17]:
            ax.scatter([x], [y], s=24, c=GRAY, zorder=2)

def band(ax, kind):
    # Panels intentionally diagrammatic; vertical axis is energy.
    arrow(ax, (73, 6), (73, 21))
    text(ax, 73, 23.2, "Energía", fontsize=12)
    if kind == "metal":
        ax.add_patch(Rectangle((79, 6.8), 14, 12, fc="white", ec=BLUE, lw=1.4))
        ax.add_patch(Rectangle((79, 6.8), 14, 6, fc=BLUE+"35", ec="none"))
        ax.plot([77, 96], [12.8]*2, "--", color=INK, lw=1.1)
        text(ax, 98, 12.8, r"$\mu$", ha="left")
        text(ax, 86, 2.2, "Banda parcial", fontsize=12)
    elif kind in ("intrinsic", "doped"):
        ax.add_patch(Rectangle((79, 15), 14, 4, fc="#e8f1f5", ec=BLUE, lw=1.4))
        ax.add_patch(Rectangle((79, 6.5), 14, 4, fc=BLUE+"45", ec=BLUE, lw=1.4))
        text(ax, 97, 16.8, "c", ha="left")
        text(ax, 97, 8.5, "v", ha="left")
        arrow(ax, (95.5, 10.5), (95.5, 15), style="<->", lw=1.2)
        text(ax, 89, 12.8, r"$E_g$")
        if kind == "intrinsic":
            electron(ax, 86, 17)
            electron(ax, 86, 8.5, True)
        else:
            for x in [82, 87, 91]:
                electron(ax, x, 17, size=35)
        text(ax, 86, 2.2, "c: conducción\nv: valencia", fontsize=12)
    else:
        ax.add_patch(Rectangle((79, 15), 14, 4.4, fc=GREEN+"30", ec=GREEN, lw=1.4))
        ax.plot([77, 94], [7]*2, "--", color=INK, lw=1)
        text(ax, 97, 7, "0", ha="left")
        text(ax, 95, 15, r"$\Delta$", ha="left")
        text(ax, 86, 11, "Sin QP", fontsize=12)
        text(ax, 86, 2.2, "Brecha de QP", fontsize=12)

fig, axes = plt.subplots(4, 1, figsize=(7.6, 9.2))
fig.subplots_adjust(left=.015, right=.97, top=.985, bottom=.012, hspace=.12)
kinds = ["metal", "intrinsic", "doped", "sc"]
titles = ["a  Metal resistivo", "b  Semiconductor intrínseco",
          "c  Semiconductor dopado tipo n", "d  Superconductor convencional"]
for i, (ax, kind, title) in enumerate(zip(axes, kinds, titles)):
    ax.set(xlim=(0, 104), ylim=(-1, 27)); ax.axis("off")
    text(ax, 1, 25.2, title, ha="left", weight="bold", fontsize=14)
    circuit(ax, kind); lattice(ax, kind == "sc"); band(ax, kind)
    if kind == "metal":
        ax.plot([59, 55, 56, 47, 48, 39], [16, 17, 14, 15, 10, 12], color=BLUE, lw=1.5)
        arrow(ax, (42, 11), (38, 12), color=BLUE)
        electron(ax, 59, 16)
        text(ax, 48, 2.2, "Deriva y dispersión\nEnergía a la red", fontsize=12)
    elif kind == "intrinsic":
        electron(ax, 53, 15)
        arrow(ax, (50, 15), (38, 15), color=BLUE)
        electron(ax, 39, 10, hole=True)
        arrow(ax, (42, 10), (59, 10), color=ORANGE)
        text(ax, 48, 2.2, "Electrón y hueco\nIgual concentración", fontsize=12)
    elif kind == "doped":
        text(ax, 43, 13, r"$D^+$", color=ORANGE,
             bbox={"boxstyle":"round,pad=.06", "fc":"white", "ec":"none"})
        electron(ax, 56, 16); electron(ax, 51, 10)
        arrow(ax, (54, 16), (45, 16), color=BLUE)
        arrow(ax, (49, 10), (36, 10), color=BLUE)
        text(ax, 48, 2.2, "Donante fijo $D^+$\nElectrones mayoritarios", fontsize=12)
    else:
        for x, y in [(38, 15), (48, 11), (57, 16)]:
            ax.add_patch(Ellipse((x, y), 7, 3.1, fc=GREEN+"15", ec=GREEN, lw=1.3))
            electron(ax, x-1.4, y, size=20); electron(ax, x+1.4, y, size=20)
        ax.scatter([59], [8.3], marker="*", s=95, color=ORANGE)
        text(ax, 48, 2.2, "Condensado y QP\nQP: estrella naranja", fontsize=12)
save(fig, "E04_circuitos_micro")

# Two separate spatial-interface and spectral descriptions.
fig, axes = plt.subplots(2, 1, figsize=(7.6, 7.3))
fig.subplots_adjust(left=.025, right=.98, top=.98, bottom=.04, hspace=.27)
for ax in axes:
    ax.set(xlim=(0, 100), ylim=(0, 35)); ax.axis("off")

ax = axes[0]
text(ax, 0, 33.6, "a  Unión p–n: semiconductor en equilibrio", ha="left", weight="bold", fontsize=14)
ax.add_patch(Rectangle((1, 11), 52, 13, fc="#f6f7f8", ec=GRAY))
ax.add_patch(Rectangle((20, 11), 16, 13, fc="#e7ddd2", ec="none"))
text(ax, 10, 26.5, "p"); text(ax, 45, 26.5, "n")
for x in [6, 13]:
    for y in [15, 20]: electron(ax, x, y, hole=True)
for x in [41, 48]:
    for y in [15, 20]: electron(ax, x, y)
for y in [15, 20]:
    text(ax, 24, y, "−", color=ORANGE, fontsize=18)
    text(ax, 32, y, "+", color=ORANGE, fontsize=16)
arrow(ax, (36, 28), (21, 28), color=ORANGE)
text(ax, 28, 30.6, r"Campo $\mathcal{E}$", fontsize=12)
arrow(ax, (45, 8.1), (12, 8.1), color=BLUE)
text(ax, 28, 5.6, "Difusión de electrones", fontsize=12)
arrow(ax, (12, 2.6), (45, 2.6), color=ORANGE)
text(ax, 28, -.1, "Difusión de huecos", fontsize=12)
arrow(ax, (61, 9), (61, 29))
text(ax, 61, 31, r"$\varepsilon$", fontsize=13)
arrow(ax, (61, 9), (98, 9))
text(ax, 99, 7, "x")
x = np.linspace(65, 94, 200)
ev = 15 - 2.8*np.tanh((x-79)/4)
ax.plot(x, ev, color=BLUE, lw=2)
ax.plot(x, ev+7, color=BLUE, lw=2)
text(ax, 96, ev[-1], r"$E_v$", ha="left")
text(ax, 96, ev[-1]+7, r"$E_c$", ha="left")
ax.plot([64, 94], [18.5, 18.5], "--", color=INK, lw=1)
text(ax, 97, 17.2, r"$\mu$", ha="left")
arrow(ax, (91, ev[-1]), (91, ev[-1]+7), style="<->", lw=1.2)
text(ax, 86, ev[-1]+3.5, r"$E_g$")
text(ax, 79, 4.5, "Bordes de banda\nfrente a posición", fontsize=12)

ax = axes[1]
text(ax, 0, 33.6, "b  Interfaz N–S: reflexión de Andreev", ha="left", weight="bold", fontsize=14)
ax.add_patch(Rectangle((1, 10), 28, 15, fc="#eaf1f5", ec=GRAY))
ax.add_patch(Rectangle((29, 10), 24, 15, fc="#eaf3ef", ec=GRAY))
text(ax, 14, 28, "N: metal normal", fontsize=12)
text(ax, 42, 28, "S", fontsize=13)
electron(ax, 7, 21)
arrow(ax, (9, 21), (28, 21), color=BLUE)
text(ax, 17, 23.7, r"$e^-$", color=BLUE)
electron(ax, 7, 14, hole=True)
arrow(ax, (28, 14), (9, 14), color=ORANGE)
text(ax, 17, 11.9, r"$h^+$", color=ORANGE)
arrow(ax, (29, 18), (47, 18), color=GREEN)
text(ax, 41, 22, r"Par: $-2e$", fontsize=12, color=GREEN)
text(ax, 27, 6.7, r"$0<E<\Delta$; contacto transparente", fontsize=12)
text(ax, 27, 2.6, "El condensado recibe un par", fontsize=12)
arrow(ax, (61, 9), (61, 29))
text(ax, 61, 31, "E", fontsize=13)
arrow(ax, (61, 9), (98, 9))
text(ax, 99, 7, "x")
ax.add_patch(Rectangle((65, 12), 13, 15, fc=BLUE+"20", ec="none"))
ax.add_patch(Rectangle((78, 23), 16, 4, fc=GREEN+"30", ec="none"))
ax.plot([65, 94], [12, 12], "--", color=INK, lw=1)
ax.plot([78, 94], [23, 23], color=GREEN, lw=1.7)
text(ax, 96, 23, r"$\Delta$", ha="left")
text(ax, 96, 12, "0", ha="left")
text(ax, 71, 7, "N", fontsize=12)
text(ax, 87, 7, "S", fontsize=12)
text(ax, 86, 18, "Sin QP\npropagantes", fontsize=12)
text(ax, 79, 2.6, r"Energía de excitación", fontsize=12)
save(fig, "E04_interfaces")

# Check the worked examples, not the unanswered activities.
n = np.ones(12, dtype=int); n[[3, 9]] = 0
h = 1-n
v = np.array([-3, -1, 1, 3])
occupied = np.array([1, 1, 0, 1])
charge_velocity_electrons = -int(np.dot(v, occupied))
charge_velocity_hole = int(np.dot(v, 1-occupied))
assert n.sum() == 10 and h.sum() == 2
assert charge_velocity_electrons == charge_velocity_hole == 1
assert 3 + (-2)*(-1) == 5
print(json.dumps({"occupied": int(n.sum()), "holes": int(h.sum()),
                  "hole_indices": (np.flatnonzero(h)+1).tolist(),
                  "electron_qv_in_e_v0": charge_velocity_electrons,
                  "hole_qv_in_e_v0": charge_velocity_hole,
                  "normal_excitation_energy_meV": 5,
                  "figures": ["E04_circuitos_micro", "E04_interfaces"]}, indent=2))
