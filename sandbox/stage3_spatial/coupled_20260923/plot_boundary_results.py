"""Fourth figure: stored reservoir postprocess versus the unchanged free pilot.

No physical module is imported. Requires completed summary JSON from the
postprocess; this figure does not grant admission to external normal flux D.27.
"""
from pathlib import Path
import argparse
import json

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from plot_results import ROOT, BASE, BLUE, ORANGE, read, sha, finite, style, panels


def build(postprocess_root, pilot_root):
    sources = {}
    result = read(postprocess_root/"summary.json", sources)
    original = read(pilot_root/"summary.json", sources)
    original_manifest = read(pilot_root/"manifest.json", sources)
    time_s = float(original_manifest["registration"]["time_reference_ps"])*1e-12
    if not np.isfinite(time_s) or time_s <= 0:
        raise ValueError("Invalid recorded time normalization")
    originals = {r["mesh"]["id"]+"_"+r["profile"]["id"]: r for r in original["cases"]}
    by_id = {row["case_id"]: row for row in result["cases"]}
    ids = [key for profile in ("helix", "smooth_perturbation")
           for key, row in originals.items() if row["profile"]["id"] == profile]
    if len(ids) != 2 or set(by_id) != set(ids):
        raise ValueError("Expected both unchanged pilot profiles and matching postprocess cases")
    rows = [by_id[key] for key in ids]
    energy_scale = finite([originals[key]["energy_scale_J"] for key in ids], "energy scale", True)
    free = finite([row["free"]["condensate_rate_bar"] for row in rows], "free heat", True)
    reservoir = finite([row["reservoir"]["condensate_rate_bar"] for row in rows], "reservoir heat", True)
    expected_free = finite([originals[key]["heat"]["condensate_rate_bar"] for key in ids], "pilot heat", True)
    if not np.allclose(free, expected_free, rtol=2e-12, atol=0.):
        raise ValueError("Postprocess free reference no longer matches the archived pilot")
    free_nW, reservoir_nW = free*energy_scale/time_s*1e9, reservoir*energy_scale/time_s*1e9
    reference = finite([row["reference_current_A"] for row in rows], "reference current", True)
    signed_current = finite([row["terminal"]["adjacent_internal_normal_current_A"]
                            for row in rows], "adjacent interior normal currents")
    if signed_current.shape != (2, 2):
        raise ValueError("One signed left/right adjacent interior current per profile required")
    normal_percent = 100*np.abs(signed_current)/reference[:, None]
    residual = finite([row["CM9"]["residual_W"] for row in rows], "CM9 residual")

    fig, axes = panels()
    fig.subplots_adjust(bottom=.29, top=.82, wspace=.46)
    positions = np.arange(2)
    for values, shift, color, name in ((free_nW, -.18, ORANGE, "Libres"),
            (reservoir_nW, .18, BLUE, "Con reservorio")):
        axes[0].bar(positions+shift, values, width=.34, color=color, label=name)
    axes[0].set_yscale("log")
    axes[0].set_ylim(min(reservoir_nW)*.4, max(free_nW)*10)
    axes[0].set_ylabel(r"Disipación integrada $Q_\Delta$ (nW)")
    axes[0].set_title("Efecto de las condiciones de borde")
    for terminal, shift, color, label in ((0, -.18, BLUE, "Izquierda"),
                                         (1, .18, ORANGE, "Derecha")):
        values = normal_percent[:, terminal]
        axes[1].bar(positions+shift, values, width=.34, color=color, label=label)
        # A zero remains zero rather than becoming an arbitrary logarithmic floor.
        axes[1].scatter(positions+shift, values, color=color, marker="o", s=20)
    axes[1].set_ylim(0, max(float(np.max(normal_percent))*1.3, 1e-5))
    axes[1].set_ylabel(r"$|I_{\rm n}|/I_{\rm ref}$ (%)")
    axes[1].set_title("Aristas interiores junto al borde")
    for ax, center in zip(axes, (.28, .785)):
        ax.set_xticks(positions, ["Hélice", "Perturbación"])
        handles, labels = ax.get_legend_handles_labels()
        fig.legend(handles, labels, loc="upper center", bbox_to_anchor=(center, 1.01),
                   ncol=2, frameon=False, handlelength=.9, handletextpad=.3, columnspacing=.7)
    fig.text(.5, .12, "Corriente interior adyacente; flujo normal externo D.27 aún no admitido.",
             ha="center", fontsize=11)
    fig.text(.5, .04, f"Residuos CM.9: {residual[0]:.1e} / {residual[1]:.1e} W. Postproceso instantáneo.",
             ha="center", fontsize=11)
    return fig, dict(sources=sources, cases=ids, free_QDelta_nW=free_nW.tolist(),
        reservoir_QDelta_nW=reservoir_nW.tolist(), reservoir_over_free=(reservoir/free).tolist(),
        signed_adjacent_internal_normal_current_A=signed_current.tolist(),
        absolute_adjacent_internal_normal_current_percent=normal_percent.tolist(),
        CM9_residual_W=residual.tolist(),
        scope="Saved instantaneous boundary postprocess. Adjacent interior current is not the external normal flux D.27. No temporal or boundary admission inferred.")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--postprocess-root", type=Path, default=BASE/"reservoir_postprocess_v2")
    parser.add_argument("--pilot-root", type=Path, default=BASE/"raw/mixed_pilot")
    parser.add_argument("--output-dir", type=Path, default=BASE/"figures")
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()
    stem = "04_reservoir_boundary"
    targets = [args.output_dir/(stem+suffix) for suffix in (".png", ".pdf", "_manifest.json")]
    if any(path.exists() for path in targets) and not args.overwrite:
        raise FileExistsError("Boundary figure exists; use --overwrite deliberately")
    style()
    fig, info = build(args.postprocess_root, args.pilot_root)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    artifacts = {}
    for path in targets[:2]:
        fig.savefig(path, dpi=300)
        artifacts[path.name] = sha(path)
    plt.close(fig)
    info.update(artifacts=artifacts, script_sha256=sha(__file__),
        shared_plotting_script_sha256=sha(Path(__file__).with_name("plot_results.py")),
        canvas_inches=[7., 3.2], dpi=300, physics_evaluated=False)
    targets[2].write_text(json.dumps(info, indent=2, ensure_ascii=False, allow_nan=False)+"\n", encoding="utf-8")
    print(json.dumps({"figure": stem, "output": str(args.output_dir), "physics_evaluated": False}))


if __name__ == "__main__":
    main()
