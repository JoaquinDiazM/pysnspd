"""Scientific figures from saved JSON only; never imports a physical solver.

Example: python sandbox/stage3_spatial/coupled_20260923/plot_results.py \
    --figures open seeded
The mixed figure requires the completed pilot's summary and manifest. Missing
data are an error, never replaced by a prediction or a placeholder curve.
"""
from pathlib import Path
import argparse
import hashlib
import json

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[3]
BASE = ROOT/"docs/implementation/stage3/coupled_20260923"
BLUE, ORANGE = "#236587", "#B65B2C"


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read(path, sources):
    path = Path(path)
    data = json.loads(path.read_text(encoding="utf-8"))
    sources[str(path.resolve().relative_to(ROOT)).replace("\\", "/")] = sha(path)
    return data


def finite(values, name, positive=False):
    result = np.asarray(values, float)
    if np.any(~np.isfinite(result)) or (positive and np.any(result <= 0)):
        raise ValueError(f"Invalid saved {name}; figures do not repair results")
    return result


def style():
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 11,
        "axes.labelsize": 11, "axes.titlesize": 11, "xtick.labelsize": 11,
        "ytick.labelsize": 11, "legend.fontsize": 11, "axes.spines.top": False,
        "axes.spines.right": False, "pdf.fonttype": 42, "ps.fonttype": 42,
        "savefig.facecolor": "white", "axes.axisbelow": True})


def panels():
    fig, axes = plt.subplots(1, 2, figsize=(7., 3.2))
    fig.subplots_adjust(left=.105, right=.98, bottom=.21, top=.86, wspace=.42)
    for ax in axes:
        ax.grid(axis="y", color="#D8DDE1", linewidth=.65)
    return fig, axes


def open_figure(folder):
    sources = {}
    manifest = read(folder/"manifest.json", sources)
    summary = read(folder/"summary.json", sources)
    cases = summary["cases"]
    if len(cases) != 6 or {c["case"]["id"] for c in cases} != {
            c["id"] for c in manifest["registration"]["cases"]}:
        raise ValueError("Expected all six registered open cases")
    rows = []
    for case in cases:
        record = read(folder/case["case"]["id"]/"result.json", sources)
        if record != case:
            raise ValueError("Open summary and case record disagree")
        length = record["case"]["length_m"]*1e9
        width = length/record["case"]["elements"]
        rows.append(dict(case=record["case"]["id"], length_nm=length,
            element_length_nm=width,
            current_relative=record["metrics"]["current_relative"],
            endpoint_q_relative=record["metrics"]["endpoint_q_relative"]))
    fig, axes = panels()
    for width, color, marker in ((90., ORANGE, "s"), (45., BLUE, "o")):
        group = sorted((r for r in rows if np.isclose(r["element_length_nm"], width)),
                       key=lambda r: r["length_nm"])
        if len(group) != 3:
            raise ValueError("Expected three lengths at each 90/45 nm element size")
        x = finite([r["length_nm"] for r in group], "length", True)
        for ax, key in zip(axes, ("current_relative", "endpoint_q_relative")):
            y = finite([r[key] for r in group], key, True)*100
            ax.plot(x, y, marker=marker, color=color, label=f"{width:g} nm", lw=1.5, ms=5)
    criteria = manifest["registration"]["criteria"]
    for ax, key, title in zip(axes, ("current_relative", "endpoint_q_relative"),
                            ("Corriente respecto al reservorio", r"Gradiente $q$ en los terminales")):
        limit = float(criteria[key])*100
        ax.axhline(limit, color="#77828A", ls="--", lw=1.)
        ax.set_yscale("log")
        ax.set_xticks([360, 720, 1080])
        ax.set_xlabel("Longitud del dominio (nm)")
        ax.set_ylabel("Error relativo (%)")
        ax.set_title(title)
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", bbox_to_anchor=(.51, 1.02),
               ncol=2, frameon=False, title=None)
    fig.text(.5, .025, "Longitud de elemento GLL; línea discontinua: criterio estático del 1 %.",
             ha="center", fontsize=11)
    return fig, dict(sources=sources, rows=rows,
        scope="Six saved static open cases. Element length is not uniform nodal spacing; no temporal accuracy claim.")


def seeded_figure(folder):
    sources = {}
    result = read(folder/"summary.json", sources)
    probes = read(folder/"probes.json", sources)
    if probes != result["probes"] or not probes:
        raise ValueError("Seeded summary and probe records disagree or are empty")
    if result["representation"] != "predictor_then_exact_causal_Newton":
        raise ValueError("Only the corrected causal kernel benchmark belongs in this figure")
    candidate = str(result["selected_degree"])
    times = finite([result["direct_median_kernel_seconds"],
                    result["patch_median_kernel_seconds"]], "kernel times", True)*1e3
    energy = finite([r["degrees"][candidate]["metrics"]["pointwise_relative_energy"]
                     for r in probes], "energy differences")
    derivative = finite([max(r["degrees"][candidate]["kernel_errors"][1:])
                         for r in probes], "derivative differences")
    if np.any(energy < 0) or np.any(derivative < 0):
        raise ValueError("Reported difference norms cannot be negative")
    fig, axes = panels()
    axes[0].bar([0, 1], times, width=.58, color=[ORANGE, BLUE])
    for i, value in enumerate(times):
        axes[0].text(i, value*1.035, f"{value:.1f}", ha="center", va="bottom")
    axes[0].set_ylim(0, max(times)*1.25)
    axes[0].set_xticks([0, 1], ["Causal directo", "Con semilla"])
    axes[0].set_ylabel("Tiempo mediano por kernel (ms)")
    axes[0].set_title(f"Kernel causal: aceleración {times[0]/times[1]:.2f}×")
    points = np.arange(1, len(probes)+1)
    zero_count = 0
    for values, offset, color, marker, label in ((energy, -.09, ORANGE, "s", "Energía"),
            (derivative, .09, BLUE, "o", "Derivadas (L1)")):
        positive = values > 0
        axes[1].scatter(points[positive]+offset, values[positive], color=color,
                        marker=marker, s=27, label=label)
        zero_count += int(np.count_nonzero(~positive))
    axes[1].set_yscale("log")
    axes[1].set_xlim(.5, len(probes)+.5)
    axes[1].set_xticks(points[::2])
    axes[1].set_xlabel("Punto de control fuera de nodos")
    axes[1].set_ylabel("Diferencia relativa")
    axes[1].set_title("Verificación fuera de nodos")
    handles, labels = axes[1].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", bbox_to_anchor=(.755, 1.025),
               frameon=False, ncol=2, handletextpad=.2, columnspacing=.5)
    if zero_count:
        axes[1].text(.98, .03, f"{zero_count} ceros exactos fuera del eje log",
                     transform=axes[1].transAxes, ha="right", fontsize=11)
    fig.text(.5, .025, "Velocidad del kernel; no mide el RHS ni un transiente completo.",
             ha="center", fontsize=11)
    return fig, dict(sources=sources, source_ms=float(times[0]), seeded_ms=float(times[1]),
        speedup=float(times[0]/times[1]), energy_differences=energy.tolist(),
        derivative_weighted_L1_differences=derivative.tolist(),
        scope="Same causal equations with different Newton initialization. Predictor loading is not predictor construction time.")


def mixed_figure(folder):
    sources = {}
    result = read(folder/"summary.json", sources)
    manifest = read(folder/"manifest.json", sources)
    reg = manifest["registration"]
    rows = result["cases"]
    expected = ("helix", "smooth_perturbation")
    if len(rows) != 2 or {r["profile"]["id"] for r in rows} != set(expected):
        raise ValueError("Mixed figure requires exactly the two completed pilot profiles")
    rows = sorted(rows, key=lambda r: expected.index(r["profile"]["id"]))
    if any(r["mesh"]["id"] != reg["pilot_mesh"] for r in rows):
        raise ValueError("Mixed figure cannot silently mix mesh resolutions")
    for row in rows:
        saved = read(folder/(row["mesh"]["id"]+"_"+row["profile"]["id"])/"result.json", sources)
        if saved != row:
            raise ValueError("Mixed summary and individual result disagree")
    tref_s = float(reg["time_reference_ps"])*1e-12
    if not np.isfinite(tref_s) or tref_s <= 0:
        raise ValueError("Invalid recorded reference time")
    scale = finite([r["energy_scale_J"]/tref_s for r in rows], "power scale", True)
    # bar powers -> W -> nW. No spectrum is loaded or re-evaluated.
    qdelta = finite([r["heat"]["condensate_rate_bar"] for r in rows], "QDelta")*scale*1e9
    joule = finite([r["potential"]["normal_joule_W"] for r in rows], "Joule")*1e9
    material = finite([r["maximum_material_speed_bar"] for r in rows], "material velocity")
    field = finite([r["maximum_field_speed_bar"] for r in rows], "field velocity")
    if np.any(np.r_[qdelta, joule, material, field] < 0):
        raise ValueError("Heat and speed norms must be nonnegative")
    fig, axes = panels()
    fig.subplots_adjust(bottom=.31, top=.89, wspace=.58)
    xpos = np.arange(2)
    for index, offset, color, name in ((0, -.18, BLUE, "Hélice"),
                                      (1, .18, ORANGE, "Perturbación")):
        axes[0].bar(xpos+offset, [qdelta[index], joule[index]], width=.34,
                    color=color, label=name)
    if np.any(np.r_[qdelta, joule] == 0):
        raise ValueError("This logarithmic heat panel requires positive saved powers")
    axes[0].set_yscale("log")
    axes[0].set_ylim(min(np.r_[qdelta, joule])*.45, max(np.r_[qdelta, joule])*7)
    axes[0].set_xticks(xpos, [r"$Q_\Delta$", "Joule"])
    axes[0].set_ylabel("Potencia depositada (nW)")
    axes[0].set_title("Extremos libres: calor")
    axes[0].legend(frameon=False, loc="upper right", handlelength=1., handletextpad=.4)
    base = np.array([qdelta[0], joule[0], material[0], field[0]])
    perturbed = np.array([qdelta[1], joule[1], material[1], field[1]])
    if np.any(base <= 0):
        raise ValueError("A relative change requires positive helix reference observables")
    change = 100*(perturbed-base)/base
    axes[1].barh(np.arange(4), change, height=.5,
                 color=[BLUE, ORANGE, BLUE, ORANGE])
    axes[1].scatter(change, np.arange(4), s=21, color="#253945", zorder=3)
    axes[1].axvline(0., color="#77828A", lw=.8)
    axes[1].set_xscale("symlog", linthresh=1e-4)
    axes[1].set_yticks(np.arange(4), [r"$Q_\Delta$", "Joule", r"$v_{\mathrm{mat}}$", r"$v_{\mathrm{campo}}$"])
    axes[1].set_ylim(3.65, -.45)
    axes[1].set_xlim(-.01, 200)
    axes[1].set_xticks([-.001, 0., .01, 1., 100.], ["−0.001", "0", "0.01", "1", "100"])
    axes[1].set_xlabel("Variación (%) · symlog")
    axes[1].set_title("Cambio frente a hélice")
    axes[1].grid(False, axis="y")
    axes[1].grid(True, axis="x", color="#D8DDE1", linewidth=.65)
    for y, value in enumerate(change):
        label = "0 %" if value == 0 else f"{value:+.3g} %"
        axes[1].text(.97, y+.23, label, transform=axes[1].get_yaxis_transform(),
                     fontsize=11, ha="right", va="top")
    residual = finite([r["CM9"]["residual_W"] for r in rows], "CM9 residual")
    fig.text(.5, .115,
        f"Hélice: velocidad máxima material {material[0]:.3f}; campo {field[0]:.3f} "
        f"(τ=t/{reg['time_reference_ps']:g} ps).", ha="center", fontsize=11)
    fig.text(.5, .035,
        f"CM.9: residuos {residual[0]:.1e} / {residual[1]:.1e} W. Sin evolución temporal.",
        ha="center", fontsize=11)
    return fig, dict(sources=sources, profiles=list(expected), QDelta_nW=qdelta.tolist(),
        Joule_nW=joule.tolist(), material_speed_bar=material.tolist(), field_speed_bar=field.tolist(),
        CM9_residual_W=residual.tolist(), perturbation_relative_changes_percent=change.tolist(),
        scope=reg["scope"],
        warning="Instantaneous free boundary response. Neither profile is a DC solution or an evolved trajectory.")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--figures", nargs="+", choices=("open", "seeded", "mixed"),
                        default=["open", "seeded", "mixed"])
    parser.add_argument("--open-root", type=Path, default=BASE/"raw/stage3b_open_full_20260923")
    parser.add_argument("--seed-root", type=Path, default=BASE/"seeded_check")
    parser.add_argument("--mixed-root", type=Path, default=BASE/"raw/mixed_pilot")
    parser.add_argument("--output-dir", type=Path, default=BASE/"figures")
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()
    style()
    sources = {"open": args.open_root, "seeded": args.seed_root, "mixed": args.mixed_root}
    builders = {"open": open_figure, "seeded": seeded_figure, "mixed": mixed_figure}
    filenames = {"open": "01_open_errors", "seeded": "02_causal_kernel", "mixed": "03_mixed_instantaneous"}
    chosen = list(dict.fromkeys(args.figures))
    # Validate required inputs and target paths before writing any figure.
    for name in chosen:
        if not (sources[name]/"summary.json").is_file():
            raise FileNotFoundError(f"Saved {name} summary missing: {sources[name]}/summary.json")
        for suffix in ("png", "pdf"):
            target = args.output_dir/(filenames[name]+"."+suffix)
            if target.exists() and not args.overwrite:
                raise FileExistsError(f"Figure exists; use --overwrite deliberately: {target}")
    args.output_dir.mkdir(parents=True, exist_ok=True)
    manifest_path = args.output_dir/"figures_manifest.json"
    manifest = (json.loads(manifest_path.read_text(encoding="utf-8")) if manifest_path.exists()
                else {"method": "Saved data only; no physical calculations", "figures": {}})
    for name in chosen:
        fig, info = builders[name](sources[name])
        artifacts = {}
        for suffix in ("png", "pdf"):
            target = args.output_dir/(filenames[name]+"."+suffix)
            fig.savefig(target, dpi=300)  # Fixed 7 x 3.2 inch scientific canvas.
            artifacts[target.name] = sha(target)
        plt.close(fig)
        manifest["figures"][name] = dict(info, artifacts=artifacts,
            script_sha256=sha(__file__), canvas_inches=[7., 3.2], dpi=300)
    manifest_path.write_text(json.dumps(manifest, indent=2, ensure_ascii=False, allow_nan=False)+"\n",
                             encoding="utf-8")
    print(json.dumps({"figures": chosen, "output": str(args.output_dir), "physics_evaluated": False}))


if __name__ == "__main__":
    main()
