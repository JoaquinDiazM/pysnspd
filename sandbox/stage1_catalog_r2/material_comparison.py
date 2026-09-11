"""Independent R2 comparison of published phonon inputs and the actual legacy loader.

This is a phonon-DOS comparison, not an electronic Usadel-DOS comparison. It
calls the unchanged production reader with its PRE defaults, reports its
transformations, and calls strict experimental admission on the original bytes.
No material is repaired or used to run a detector transient.
"""

from __future__ import annotations

import argparse
import hashlib
import io
import json
import math
from pathlib import Path
import platform
import sys
import time
import urllib.request

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
import pysnspd.kinetic.eliashberg as legacy_module
import pysnspd.experimental.material_admission as admission_module
from pysnspd.experimental.material_admission import (
    MaterialSpec, MaterialAdmissionError, audit_material_table, admit_material_table,
)

PIN = "5b6bd747f80016da5ccd51db73c110a8ecc6abf6"
PUBLIC_SHA = "9aeea0948033d771deedae40da0cb4dc59fef80ac6ea41ac4d3a67b180efc610"
LOCAL_SHA = "e94f11273e6c42b7c8fad978b78d119166e241bec1614ee77e4176906fa8f5df"
PUBLIC_REPO = "https://github.com/qnngroup/proj-KE-solver"


def digest(raw):
    return hashlib.sha256(raw).hexdigest()


def fetch(url):
    request = urllib.request.Request(url, headers={"User-Agent": "pysnspd-source-audit-r2"})
    return urllib.request.urlopen(request, timeout=12).read()


def area_independent(x, y):
    """Scalar compensated sum, independently checked against NumPy trapezoid."""
    return math.fsum((float(b)-float(a))*(float(v)+float(u))/2
                     for a,b,u,v in zip(x[:-1],x[1:],y[:-1],y[1:]))


def draw_pair(axis, x, raw, legacy, *, every=280):
    # A wide pale line plus narrow dark dashes makes agreement visible.
    axis.plot(x, legacy, color="#d78532", lw=2.7, alpha=.9, label="Legado (PRE)")
    axis.plot(x, raw, color="#1d3449", lw=1.4, ls=(0,(5,3)), label="Simón 2025 (archivo)")
    axis.plot(x[::every], legacy[::every], ls="none", marker="o", ms=4,
              mfc="white", mec="#b96917", mew=1, label="_nolegend_")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--nbn-path", type=Path, required=True)
    parser.add_argument("--output-root", type=Path, default=ROOT / "docs/implementation/stage1_r2")
    parser.add_argument("--offline", action="store_true", help="Use pinned SHA evidence without a fresh upstream query.")
    args = parser.parse_args()
    start = time.perf_counter()
    out = args.output_root
    figures = out / "figures"
    figures.mkdir(parents=True, exist_ok=True)
    watched = {
        "source": args.nbn_path,
        "legacy_loader": Path(legacy_module.__file__),
        "experimental_admission": Path(admission_module.__file__),
        "legacy_power_consumer": ROOT / "pysnspd/kinetic/power_table.py",
        "legacy_pre_pipeline": ROOT / "pipelines/01_prerun_template.py",
    }
    initial_hashes = {name:digest(path.read_bytes()) for name,path in watched.items()}
    raw = args.nbn_path.read_bytes()
    first, separator, payload = raw.partition(b"\n")
    table = np.loadtxt(io.BytesIO(raw), comments="#")
    checks = {}

    def check(name, condition):
        checks[name] = bool(condition)
        if not condition:
            raise AssertionError(name)

    check("exact_local_source_hash", digest(raw) == LOCAL_SHA)
    check("published_payload_hash_after_removing_added_header", digest(payload) == PUBLIC_SHA)
    check("local_header_added_48_bytes", len(first+separator) == 48)
    online = {"checked": False}
    if not args.offline:
        published = fetch(f"https://raw.githubusercontent.com/qnngroup/proj-KE-solver/{PIN}/nbn-a2f-ph.dat")
        check("fresh_public_bytes_equal_local_numeric_payload", published == payload)
        public_table = np.loadtxt(io.BytesIO(published))
        check("fresh_public_numeric_table_identical", np.array_equal(public_table, table))
        head = json.loads(fetch("https://api.github.com/repos/qnngroup/proj-KE-solver/commits/main"))
        tree = json.loads(fetch(f"https://api.github.com/repos/qnngroup/proj-KE-solver/git/trees/{head['sha']}"))
        nbn_names = [item["path"] for item in tree["tree"] if "nbn" in item["path"].lower()]
        online = {"checked":True, "main_sha":head["sha"], "main_commit_date":head["commit"]["committer"]["date"],
                  "nbn_paths":nbn_names, "missing_solver_input_present": "nbn-a2f-ph_2.dat" in nbn_names,
                  "public_numeric_file_has_header":published.splitlines()[0].lstrip().startswith(b"#")}

    # Invoke precisely the defaults used by pipelines/01_prerun_template.py.
    legacy = legacy_module.load_simon_eliashberg_dat(args.nbn_path)
    unclipped = legacy_module.load_simon_eliashberg_dat(args.nbn_path, clip_negative=False)
    order = np.argsort(table[:,0])
    axis, alpha_raw, dos_raw = table[order].T
    dos_legacy, alpha_legacy = legacy.phdos_states_per_THz, legacy.alpha2F
    check("legacy_axis_is_numpy_argsort_of_raw", np.array_equal(legacy.frequency_THz, axis))
    check("unclipped_loader_preserves_sorted_source_dos", np.array_equal(unclipped.phdos_states_per_THz,dos_raw))
    check("legacy_dos_is_exact_zero_clipping", np.array_equal(dos_legacy,np.maximum(dos_raw,0)))
    check("legacy_alpha_is_exact_zero_clipping", np.array_equal(alpha_legacy,np.maximum(alpha_raw,0)))
    check("no_phdos_renormalization", np.array_equal(dos_legacy[dos_raw>=0],dos_raw[dos_raw>=0]))
    check("alpha2F_identical_in_this_file", np.array_equal(alpha_raw,alpha_legacy))
    duplicate_count = int(np.sum(np.diff(axis)==0))
    check("duplicate_intervals_preserved", duplicate_count == int(np.sum(np.diff(legacy.frequency_THz)==0)) == 4)
    delta_dos = dos_legacy-dos_raw
    raw_area = area_independent(table[:,0],table[:,2])
    sorted_raw_area = area_independent(axis,dos_raw)
    clipped_area = area_independent(axis,dos_legacy)
    removed_negative_area = area_independent(axis,np.minimum(dos_raw,0))
    for name,values in (("raw",dos_raw),("legacy",dos_legacy)):
        check(f"independent_{name}_quadrature", abs(area_independent(axis,values)-np.trapezoid(values,axis)) < 5e-15)
    check("clipping_area_change_equals_removed_negative_area", abs(clipped_area-sorted_raw_area+removed_negative_area)<5e-15)
    check("loader_metadata_matches_measured_integral", abs(legacy.metadata["integral_phdos_dTHz"]-clipped_area)<5e-15)

    spec = MaterialSpec(
        material_id="NbN-Simon-R2-quarantine", axis_unit="THz", phdos_unit="states/THz",
        normalization_basis="atom", atoms_per_basis=1,basis_density_m3=48e27,
        source_url=f"{PUBLIC_REPO}/blob/{PIN}/nbn-a2f-ph.dat",source_revision=PIN,raw_sha256=LOCAL_SHA,
        unit_evidence="Local added header is absent from the public source; not certified.",
        normalization_evidence="Atom/cell basis and complete phonon mode count are not certified.",
        provenance_evidence="Independent R2 byte comparison against pinned public numeric payload.",
        provenance_verified=True,
    )
    audit = audit_material_table(args.nbn_path,spec)
    try:
        admit_material_table(args.nbn_path,spec)
    except MaterialAdmissionError:
        experimental_rejected=True
    else:
        experimental_rejected=False
    check("experimental_original_input_rejected", experimental_rejected)
    check("experimental_repairs_empty", audit.to_dict()["repairs_applied"]==[])
    check("experimental_audit_repeats_raw_integral", abs(audit.metrics_assuming_declared_units["signed_modes_per_basis"]-raw_area)<5e-15)

    # Check the actual interpolation wrapper on a fresh dense axis, without
    # computing the expensive PRE kernels or claiming this is a saved PRE grid.
    probe_axis = np.linspace(float(axis[0]),float(axis[-1]),4001)
    probe_omega = legacy_module.thz_to_j(probe_axis)
    queried_dos = legacy.phdos_on_omega_J(probe_omega)
    expected_query = np.interp(probe_omega,legacy.omega_J,dos_legacy,left=0.,right=0.)
    check("legacy_energy_argument_keeps_dos_ordinate_per_THz",np.array_equal(queried_dos,expected_query))
    check("omega_axis_uses_h_times_nu",np.array_equal(legacy.omega_J,legacy_module.thz_to_j(axis)))
    check("outside_source_support_returns_zero",np.array_equal(legacy.phdos_on_omega_J(np.array([-1e-22,2*legacy.omega_J[-1]])),np.zeros(2)))

    negative = dos_raw<0
    tail_start = max(17.0,float(axis[negative].min())-.2)
    plt.rcParams.update({"font.size":11,"axes.labelsize":11,"axes.titlesize":12,
                         "xtick.labelsize":10.5,"ytick.labelsize":10.5,
                         "axes.spines.top":False,"axes.spines.right":False})
    fig,ax=plt.subplots(2,2,figsize=(9.5,8.5))
    fig.subplots_adjust(left=.10,right=.98,bottom=.15,top=.875,hspace=.40,wspace=.30)
    fig.suptitle("DOS fonónica: archivo de Simón (2025) y legado",fontsize=14,y=.975)
    fig.text(.5,.925,"El legado recorta la cola negativa; las demás ordenadas coinciden.",ha="center",fontsize=11)
    draw_pair(ax[0,0],axis,dos_raw,dos_legacy)
    ax[0,0].set(title="a) Espectro completo",xlabel="Frecuencia declarada [THz]",ylabel="DOS tabulada")
    ax[0,0].legend(fontsize=10.5,loc="upper right")
    ax[0,0].text(.50,.61,"Coincidencia\npara DOS ≥ 0",transform=ax[0,0].transAxes,fontsize=11,
                 bbox=dict(facecolor="white",edgecolor="none",alpha=.9))
    draw_pair(ax[0,1],axis,dos_raw,dos_legacy,every=15)
    ax[0,1].set(xlim=(tail_start,20),ylim=(-.0072,.0012),title="b) Zoom de la cola negativa",xlabel="Frecuencia declarada [THz]",ylabel="DOS tabulada")
    ax[0,1].axhline(0,color="#888888",lw=.6,zorder=0)
    ax[0,1].text(.04,.25,"Archivo: cola negativa\nLegado: cero",transform=ax[0,1].transAxes,fontsize=11,
                 bbox=dict(facecolor="white",edgecolor="none",alpha=.9))
    ax[0,1].legend(fontsize=10.5,loc="lower left",bbox_to_anchor=(.0,-.015))
    ax[1,0].plot(axis,1000*delta_dos,color="#a13636",lw=1.8)
    ax[1,0].fill_between(axis,1000*delta_dos,0,color="#c66d62",alpha=.18)
    ax[1,0].set(title="c) Residual: legado − archivo",xlabel="Frecuencia declarada [THz]",ylabel="Diferencia × 1000")
    ax[1,0].text(.03,.76,f"361 muestras recortadas\nÁrea: {raw_area:.5f} → {clipped_area:.5f}\nSin renormalizar",transform=ax[1,0].transAxes,fontsize=11,
                 bbox=dict(facecolor="white",edgecolor="none",alpha=.9))
    draw_pair(ax[1,1],axis,alpha_raw,alpha_legacy)
    ax[1,1].set(title=r"d) Interacción $\alpha^2F$: sin cambio",xlabel="Frecuencia declarada [THz]",ylabel=r"$\alpha^2F$ tabulada")
    ax[1,1].text(.55,.85,"Residual = 0\nCoinciden",transform=ax[1,1].transAxes,fontsize=11,
                 bbox=dict(facecolor="white",edgecolor="none",alpha=.9))
    fig.text(.5,.045,"Fonones: no representa la DOS electrónica de Usadel.\nEje según cabecera local; unidad de ordenada y normalización pendientes.",
             ha="center",fontsize=11)
    for suffix in ("png","pdf"):
        fig.savefig(figures/f"material_legacy_comparison.{suffix}",dpi=200)
    plt.close(fig)

    cumulative_raw=np.r_[0,np.cumsum(np.diff(axis)*(dos_raw[1:]+dos_raw[:-1])/2)]
    cumulative_legacy=np.r_[0,np.cumsum(np.diff(axis)*(dos_legacy[1:]+dos_legacy[:-1])/2)]
    np.savetxt(out/"material_curves.csv",np.c_[axis,dos_raw,dos_legacy,delta_dos,alpha_raw,alpha_legacy,cumulative_raw,cumulative_legacy],
               delimiter=",",header="axis_as_supplied,phdos_raw_sorted,phdos_legacy,phdos_legacy_minus_raw,alpha2F_raw_sorted,alpha2F_legacy,cumulative_raw,cumulative_legacy",comments="")
    np.savetxt(out/"material_interpolation_probe.csv",np.c_[probe_axis,probe_omega,queried_dos],delimiter=",",
               header="frequency_THz_declared,omega_J_legacy,phdos_returned_in_legacy_states_per_THz",comments="")
    prior_path=ROOT/"docs/implementation/stage1/material_results.json"
    prior=json.loads(prior_path.read_text(encoding="utf-8"))
    prior_metrics=prior["nbn_audit"]["metrics_assuming_declared_units"]
    repeat={
        "stage1_result_sha256":digest(prior_path.read_bytes()),
        "delta_raw_integral":raw_area-prior_metrics["signed_modes_per_basis"],
        "delta_negative_fraction":(-removed_negative_area/area_independent(table[:,0],np.maximum(table[:,2],0)))-prior_metrics["negative_to_positive_fraction"],
        "duplicate_pairs_before":prior_metrics["zero_width_intervals"],"duplicate_pairs_r2":duplicate_count,
    }
    check("independent_r2_reproduces_stage1_raw_integral",abs(repeat["delta_raw_integral"])<5e-15)
    check("independent_r2_reproduces_stage1_negative_fraction",abs(repeat["delta_negative_fraction"])<5e-15)
    final_hashes={name:digest(path.read_bytes()) for name,path in watched.items()}
    check("raw_input_and_production_files_unchanged",initial_hashes==final_hashes)
    result={
        "schema":"pysnspd.stage1_r2.material-comparison.v1","runtime_seconds":time.perf_counter()-start,
        "python":platform.python_version(),"numpy":np.__version__,"platform":platform.platform(),
        "hashes":initial_hashes,"script_sha256":digest(Path(__file__).read_bytes()),
        "fresh_upstream":online,"primary_repository":PUBLIC_REPO,"pinned_revision":PIN,
        "public_file_sha256":PUBLIC_SHA,"added_header_bytes":len(first+separator),
        "physical_quantity":"phonon DOS and alpha2F input; NOT the electronic Usadel DOS",
        "legacy_transforms":{"uses_actual_default_loader":True,"sorts_axis":True,"merges_duplicates":False,
                             "rows_moved_by_argsort":int(np.sum(order!=np.arange(len(order)))),
                             "plot_alignment":"Source columns reordered with the same argsort permutation for pointwise comparison; original-order integral reported separately.",
                             "phdos_negative_samples_clipped":int(legacy.metadata["n_phdos_negative_clipped"]),
                             "alpha2F_negative_samples_clipped":int(legacy.metadata["n_alpha2F_negative_clipped"]),
                             "renormalizes_phdos":False,"duplicate_axis_pairs_retained":duplicate_count,
                             "energy_argument_change":"Omega=h*nu; ordinate remains in the declared states/THz convention",
                             "thermal_consumer":"power_table.py integrates over nu_THz and multiplies by the configured ion density; no automatic mode-count repair",
                             "PRE_call":"pipelines/01_prerun_template.py:365 calls load_simon_eliashberg_dat without a clip_negative override"},
        "metrics":{"raw_integral":raw_area,"sorted_raw_integral":sorted_raw_area,
                   "sorting_area_change":sorted_raw_area-raw_area,"legacy_integral":clipped_area,"negative_integral":removed_negative_area,
                   "negative_fraction_of_positive":-removed_negative_area/clipped_area,
                   "legacy_area_increase_relative_to_raw":clipped_area/raw_area-1,
                   "max_phdos_change":float(delta_dos.max()),"max_alpha2F_change":float(np.max(np.abs(alpha_legacy-alpha_raw))),
                   "max_change_where_raw_dos_nonnegative":float(np.max(np.abs(delta_dos[dos_raw>=0]))),
                   "lambda_legacy":float(legacy.metadata["lambda_ep"]),
                   "normalization_interpretation":"integrals under the local added header, not certified modes per atom"},
        "experimental":{"status":audit.status,"curve_produced":False,"repairs_applied":[],"audit":audit.to_dict(),
                        "decision":"The raw source is quarantined; no new phonon spectrum is manufactured."},
        "comparison_to_stage1":repeat,"checks":checks,"passed_checks":sum(checks.values()),
        "figure":"figures/material_legacy_comparison.png",
        "limitations":["No raw DFPT inputs or corrected nbn-a2f-ph_2.dat located in the published repository.",
                       "A near-integer integral after a unit hypothesis does not certify units or basis.",
                       "The interpolation probe is a fresh diagnostic grid, not a saved production PRE catalog.",
                       "This comparison establishes data transformations, not physically validated phonon rates."],
    }
    (out/"material_comparison.json").write_text(json.dumps(result,indent=2,ensure_ascii=False,allow_nan=False)+"\n",encoding="utf-8")
    notice=ROOT/"docs/implementation/stage1/MATERIAL_SOURCE_NOTICE.md"
    license_text="MIT License"+notice.read_text(encoding="utf-8").split("MIT License",1)[1]
    (out/"material_source_notice.md").write_text(
        "The R2 material CSV and figures derive from the published NbN numeric table. Column names identify original values and our diagnostic transformations. "
        "The raw third-party file remains outside this repository; experimental admission applies no repair.\n\n"
        f"Source: {PUBLIC_REPO}/blob/{PIN}/nbn-a2f-ph.dat\n\nLicense: {PUBLIC_REPO}/blob/{PIN}/LICENSE\n\n"+license_text,encoding="utf-8")
    (out/"material_findings.md").write_text(
        f"La comparación es de **DOS fonónica**, no de la DOS electrónica de Usadel. Se invocó el loader real del legado con los argumentos usados por PRE.\n\n"
        f"El legado conserva la parte positiva y recorta {int(negative.sum())} valores negativos. La integral pasa de {raw_area:.9f} a {clipped_area:.9f}; no hay renormalización. "
        f"Conserva cuatro pares de abscisas repetidas; el ordenamiento invierte un par y modifica la integral bruta en {sorted_raw_area-raw_area:.2e}, antes del recorte. "
        f"α²F coincide punto a punto después del mismo ordenamiento: por eso sus dos trazos se superponen. "
        f"La nueva figura muestra las coincidencias mediante estilos distintos y separa el residual y la cola para visualizar el cambio.\n\n"
        f"**R2 mantiene la entrada rechazada.** No existe una nueva curva fonónica experimental que superponer. La procedencia numérica se verifica contra "
        f"[el archivo público fijado por revisión]({PUBLIC_REPO}/blob/{PIN}/nbn-a2f-ph.dat), que carece de la cabecera añadida localmente. "
        f"La unidad y la base de normalización de la ordenada siguen sin certificar.\n\n"
        f"El [solver de la fuente]({PUBLIC_REPO}/blob/{PIN}/solver.m#L93-L113) referencia `nbn-a2f-ph_2.dat`, ausente del árbol publicado consultado. "
        f"Su conversión de unidades no demuestra que este archivo sin sufijo tenga las mismas convenciones. No se adopta una reinterpretación por conveniencia.\n\n"
        f"La verificación independiente R2 pasa {len(checks)} controles y reproduce la integral de la primera etapa con diferencia {repeat['delta_raw_integral']:.2e}. "
        f"Tiempo de ejecución: {result['runtime_seconds']:.3f} s.\n",encoding="utf-8")
    print(json.dumps({"runtime_seconds":result["runtime_seconds"],"passed_checks":len(checks),"raw_area":raw_area,
                      "legacy_area":clipped_area,"alpha2F_change":result["metrics"]["max_alpha2F_change"],"experimental":audit.status}))


if __name__=="__main__":
    main()
