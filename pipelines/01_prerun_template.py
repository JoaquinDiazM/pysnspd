"""PRE-run template for pySNSPD.

Official PRE-run responsibilities
---------------------------------
OE2:
    mesh, edges, boundary tags and mesh plots.

OE3:
    Usadel/material calibration and quasiparticle DOS catalogue.

OE4:
    phase-space catalogues J_S and J_R.

OE5:
    thermal Usadel grid Delta_eq(Te,q), projected electron-phonon powers,
    Eliashberg/Debye comparisons and spectral-support diagnostics.

PRE-run is the natural place for parallel work: catalogue points are
independent and reusable. Later SS/PHOTON time evolutions should load these
objects rather than recomputing them.
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Any

import numpy as np
import yaml

from pysnspd.config import load_config, validate_config
from pysnspd.io.manager import create_run_layout, write_manifest

from pysnspd.mesh.delaunay import (
    generate_rectangular_delaunay_mesh,
    mesh_summary,
    save_mesh_npz,
)
from pysnspd.mesh.edges import (
    assert_edge_data_consistent,
    build_edge_data,
    edge_summary,
    save_edges_npz,
)

from pysnspd.plotting.figures import (
    plot_boundary_tags,
    plot_mesh_geometry,
    plot_phase_space_slices,
    plot_usadel_calibration_sweep,
    plot_usadel_dos_slices,
)
from pysnspd.plotting.kinetic import (
    plot_eliashberg_spectrum,
    plot_low_energy_recombination_scattering_band,
    plot_power_curve_thermal_usadel,
    plot_power_ratios_thermal_usadel,
    plot_power_scan_thermal_usadel,
    plot_spectral_support,
    plot_thermal_usadel_gap_grid,
)

from pysnspd.usadel.catalog import (
    build_usadel_catalog_from_config,
    catalog_summary,
    save_usadel_catalog_npz,
)

from pysnspd.kinetic.eliashberg import (
    load_simon_eliashberg_dat,
    spectrum_summary,
)
from pysnspd.kinetic.phase_space import (
    build_phase_space_catalog_from_usadel_catalog,
    phase_space_summary,
    save_phase_space_catalog_npz,
)
from pysnspd.kinetic.powers import (
    TAU0_OVER_TAU_EP_TC,
    compute_power_curve_thermal_usadel_state,
    compute_power_scan_thermal_usadel,
    compute_projected_powers,
    cumulative_spectral_support,
    electronic_density_of_states_from_sigma_D,
    select_thermal_usadel_q_state,
    tau0_from_tau_ep_Tc,
    tau_ep_Tc_from_tau0,
    thermal_usadel_delta_at_state,
)
from pysnspd.kinetic.thermal_usadel import (
    build_thermal_usadel_grid_parallel,
    save_thermal_usadel_grid_npz,
    thermal_usadel_grid_summary,
)


MEV_J = 1.602176634e-22


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Generate official PRE-run catalogues and diagnostics."
    )

    parser.add_argument("--config", required=True, help="Path to YAML project config.")
    parser.add_argument(
        "--run-name",
        default=None,
        help="Optional run name. If omitted, project.default_run_name is used.",
    )

    parser.add_argument(
        "--workers",
        type=int,
        default=16,
        help="Number of process workers for PRE-run catalogue tasks.",
    )

    parser.add_argument("--jitter-fraction", type=float, default=0.10)
    parser.add_argument("--boundary-guard-layers", type=int, default=1)

    parser.add_argument("--eta-fraction", type=float, default=1.0e-3)
    parser.add_argument("--gamma-max-fraction", type=float, default=0.80)
    parser.add_argument("--energy-max-factor", type=float, default=30.0)

    parser.add_argument("--phase-n-Te", type=int, default=6)
    parser.add_argument("--phase-n-delta", type=int, default=6)
    parser.add_argument("--phase-n-q", type=int, default=6)
    parser.add_argument("--phase-n-omega", type=int, default=480)
    parser.add_argument("--phase-omega-max-meV", type=float, default=35.0)
    parser.add_argument("--phase-Te-min-K", type=float, default=None)
    parser.add_argument("--phase-Te-max-K", type=float, default=None)

    parser.add_argument(
        "--skip-oe5",
        action="store_true",
        help="Skip OE5 thermal Usadel and projected-power diagnostics.",
    )
    parser.add_argument(
        "--eliashberg-dat",
        default=None,
        help=(
            "Path to Simon/MIT nbn-a2f-ph.dat. If omitted, the PRE-run tries "
            "big_data_root/catalogs/simon_2025/nbn-a2f-ph.dat."
        ),
    )
    parser.add_argument(
        "--oe5-Tph-K",
        type=float,
        default=None,
        help=(
            "OE5 phonon temperature in K. If omitted, use bias.T_bias_K "
            "from the config, falling back to 0.9 K."
        ),
    )
    parser.add_argument("--oe5-Te-min-K", type=float, default=0.9)
    parser.add_argument("--oe5-Te-max-K", type=float, default=34.6)
    parser.add_argument("--oe5-n-Te", type=int, default=180)
    parser.add_argument(
        "--oe5-q-scan-Te-K",
        type=float,
        nargs="*",
        default=None,
        help="Temperatures used for the thermal Usadel q-scan plot.",
    )
    parser.add_argument("--oe5-n-q-thermal", type=int, default=140)
    parser.add_argument("--oe5-n-matsubara-thermal", type=int, default=500)
    parser.add_argument(
        "--oe5-state-current-fraction",
        type=float,
        default=None,
        help=(
            "Select representative q by reference I(q,T_bias)/Ic(T_bias). "
            "If omitted, use I_bias/Ic from the Usadel catalogue."
        ),
    )
    parser.add_argument("--oe5-omega-max-meV", type=float, default=35.0)
    parser.add_argument("--oe5-tau-ep-Tc-ps", type=float, default=24.7)
    parser.add_argument("--oe5-tau0-ps", type=float, default=None)
    parser.add_argument(
        "--oe5-support-min-delta-fraction",
        type=float,
        default=0.05,
        help=(
            "Choose the spectral-support representative Te as the highest "
            "temperature where Delta_eq(Te,q) remains above this fraction of "
            "its maximum along the selected q-state."
        ),
    )
    parser.add_argument(
        "--oe5-low-energy-max-meV",
        type=float,
        default=None,
        help="Optional x-limit for low-energy S/R partial-power plot.",
    )

    return parser.parse_args()


def main() -> int:
    args = parse_args()

    cfg = validate_config(load_config(args.config))
    layout = create_run_layout(cfg, args.run_name)

    run_name = layout["run_name"]
    raw_pre = Path(layout["raw_pre"])
    plots_mesh = Path(layout["plots_mesh"])
    plots_diagnostics = Path(layout["plots_diagnostics"])
    plots_comparisons = Path(layout["plots_comparisons"])

    workers = max(1, int(args.workers))

    # OE2: mesh
    mesh = generate_rectangular_delaunay_mesh(
        cfg,
        jitter_fraction=args.jitter_fraction,
        boundary_guard_layers=args.boundary_guard_layers,
    )
    edge_data = build_edge_data(
        mesh.nodes,
        mesh.triangles,
        length_m=mesh.length_m,
        width_m=mesh.width_m,
    )
    assert_edge_data_consistent(edge_data)

    mesh_npz = save_mesh_npz(mesh, raw_pre / "mesh.npz")
    edges_npz = save_edges_npz(edge_data, raw_pre / "edges.npz")

    mesh_edge_summary = {
        "run_name": run_name,
        "mesh": mesh_summary(mesh),
        "edges": edge_summary(edge_data),
    }

    mesh_summary_path = raw_pre / "mesh_summary.yaml"
    with mesh_summary_path.open("w", encoding="utf-8") as f:
        yaml.safe_dump(mesh_edge_summary, f, sort_keys=False, allow_unicode=True)

    mesh_plot = plot_mesh_geometry(mesh, edge_data, plots_mesh / "mesh_nodes_edges.png")
    tags_plot = plot_boundary_tags(mesh, edge_data, plots_mesh / "mesh_boundary_tags.png")

    # OE3: Usadel/DOS
    usadel_catalog = build_usadel_catalog_from_config(
        cfg,
        eta_fraction=args.eta_fraction,
        gamma_max_fraction=args.gamma_max_fraction,
        energy_max_factor=args.energy_max_factor,
    )
    usadel_npz = save_usadel_catalog_npz(
        usadel_catalog,
        raw_pre / "usadel_dos_catalog.npz",
    )

    usadel_summary = catalog_summary(usadel_catalog)
    usadel_summary_path = raw_pre / "usadel_dos_summary.yaml"
    with usadel_summary_path.open("w", encoding="utf-8") as f:
        yaml.safe_dump(
            {
                "run_name": run_name,
                "usadel": usadel_summary,
                "metadata": usadel_catalog.metadata,
            },
            f,
            sort_keys=False,
            allow_unicode=True,
        )

    usadel_plot = plot_usadel_dos_slices(
        usadel_catalog,
        plots_diagnostics / "usadel_dos_slices.png",
    )
    calibration_plot = plot_usadel_calibration_sweep(
        usadel_catalog,
        plots_diagnostics / "usadel_calibration_sweep.png",
    )

    # OE4: phase-space
    phase_catalog = build_phase_space_catalog_from_usadel_catalog(
        usadel_catalog,
        cfg,
        n_Te=args.phase_n_Te,
        n_delta=args.phase_n_delta,
        n_q=args.phase_n_q,
        n_omega=args.phase_n_omega,
        Te_min_K=args.phase_Te_min_K,
        Te_max_K=args.phase_Te_max_K,
        omega_max_meV=args.phase_omega_max_meV,
    )

    phase_npz = save_phase_space_catalog_npz(
        phase_catalog,
        raw_pre / "phase_space_catalog.npz",
    )
    phase_summary_data = phase_space_summary(phase_catalog)
    phase_summary_path = raw_pre / "phase_space_summary.yaml"
    with phase_summary_path.open("w", encoding="utf-8") as f:
        yaml.safe_dump(
            {
                "run_name": run_name,
                "phase_space": phase_summary_data,
                "metadata": phase_catalog.metadata,
            },
            f,
            sort_keys=False,
            allow_unicode=True,
        )

    phase_plot = plot_phase_space_slices(
        phase_catalog,
        plots_diagnostics / "phase_space_slices.png",
    )

    oe5_outputs: dict[str, str] = {}
    oe5_summary: dict[str, Any] | None = None

    if not args.skip_oe5:
        oe5_summary, oe5_outputs = _run_oe5_prerun_block(
            cfg=cfg,
            args=args,
            raw_pre=raw_pre,
            plots_diagnostics=plots_diagnostics,
            plots_comparisons=plots_comparisons,
            usadel_catalog=usadel_catalog,
            phase_catalog=phase_catalog,
            workers=workers,
        )

    outputs = {
        "mesh_npz": str(mesh_npz),
        "edges_npz": str(edges_npz),
        "mesh_summary": str(mesh_summary_path),
        "mesh_plot": str(mesh_plot),
        "boundary_tags_plot": str(tags_plot),
        "usadel_npz": str(usadel_npz),
        "usadel_summary": str(usadel_summary_path),
        "usadel_plot": str(usadel_plot),
        "calibration_plot": str(calibration_plot),
        "phase_space_npz": str(phase_npz),
        "phase_space_summary": str(phase_summary_path),
        "phase_space_plot": str(phase_plot),
    }
    outputs.update(oe5_outputs)

    manifest_path = write_manifest(
        cfg,
        run_name,
        stage="pre",
        extra={
            "pipeline": "01_prerun_template.py",
            "purpose": (
                "Official PRE-run: OE2 mesh, OE3 Usadel/DOS, OE4 phase-space, "
                "OE5 thermal Usadel projected powers."
            ),
            "workers": workers,
            "outputs": outputs,
            "mesh_edge_summary": mesh_edge_summary,
            "usadel_summary": usadel_summary,
            "phase_space_summary": phase_summary_data,
            "oe5_summary": oe5_summary,
        },
    )

    print("PRE-run generation")
    print(f"run_name              : {run_name}")
    print(f"workers               : {workers}")
    print(f"raw_pre               : {raw_pre}")
    print(f"plots_mesh            : {plots_mesh}")
    print(f"plots_diagnostics     : {plots_diagnostics}")
    print(f"plots_comparisons     : {plots_comparisons}")
    print()

    print("Mesh summary")
    _print_dict(mesh_edge_summary["mesh"])
    print()

    print("Edge summary")
    _print_dict(mesh_edge_summary["edges"])
    print()

    print("Usadel/DOS summary")
    _print_dict(usadel_summary)
    print()

    calibration_warnings = usadel_summary.get("calibration_warnings", [])
    if calibration_warnings:
        print("Calibration warnings")
        for warning in calibration_warnings:
            print(f"  WARNING: {warning}")
    else:
        print("Calibration warnings: none")
    print()

    print("Phase-space summary")
    _print_dict(phase_summary_data)
    print()

    if oe5_summary is not None:
        print("OE5 projected-power summary")
        _print_dict(oe5_summary.get("representative_power", {}))
        print()

    print("Outputs")
    for key, value in outputs.items():
        print(f"  {key}: {value}")
    print(f"  pre_manifest: {manifest_path}")
    print("Status: OK")

    return 0


def _run_oe5_prerun_block(
    *,
    cfg: dict[str, Any],
    args: argparse.Namespace,
    raw_pre: Path,
    plots_diagnostics: Path,
    plots_comparisons: Path,
    usadel_catalog,
    phase_catalog,
    workers: int,
) -> tuple[dict[str, Any], dict[str, str]]:
    eliashberg_path = _resolve_eliashberg_path(cfg, args.eliashberg_dat)
    spectrum = load_simon_eliashberg_dat(eliashberg_path)

    Tph_K = _resolve_oe5_Tph_K(cfg, args)

    D_m2_s = float(usadel_catalog.metadata["D_m2_s"])
    sigma_n = float(usadel_catalog.metadata["sigma_n_S_m"])
    Tc_K = float(usadel_catalog.metadata["Tc_K"])
    T_bias_K = float(usadel_catalog.metadata["T_bias_K"])
    I_bias_A = float(usadel_catalog.metadata.get("I_bias_A", 0.0))

    calibration = usadel_catalog.metadata.get("calibration", {})
    Ic_A = float(calibration.get("Ic_model_A", calibration.get("Ic_target_A", np.nan)))

    N0 = electronic_density_of_states_from_sigma_D(sigma_n, D_m2_s)

    if args.oe5_tau0_ps is not None:
        tau0_s = float(args.oe5_tau0_ps) * 1.0e-12
        tau_ep_Tc_s = tau_ep_Tc_from_tau0(tau0_s)
        tau_source = "direct_tau0"
    else:
        tau_ep_Tc_s = float(args.oe5_tau_ep_Tc_ps) * 1.0e-12
        tau0_s = tau0_from_tau_ep_Tc(tau_ep_Tc_s)
        tau_source = "converted_from_tau_ep_Tc"

    Te_values = np.linspace(
        float(args.oe5_Te_min_K),
        float(args.oe5_Te_max_K),
        int(args.oe5_n_Te),
    )

    if args.oe5_q_scan_Te_K:
        q_scan_Te = np.asarray(args.oe5_q_scan_Te_K, dtype=float)
    else:
        q_scan_Te = np.asarray(
            [
                Tph_K,
                0.8 * Tc_K,
                Tc_K,
                min(2.0 * Tc_K, float(args.oe5_Te_max_K)),
                float(args.oe5_Te_max_K),
            ],
            dtype=float,
        )

    q_scan_Te = np.unique(
        np.clip(q_scan_Te, float(args.oe5_Te_min_K), float(args.oe5_Te_max_K))
    )
    thermal_Te_axis = np.unique(np.concatenate([Te_values, q_scan_Te]))

    thermal_grid = build_thermal_usadel_grid_parallel(
        usadel_catalog,
        thermal_Te_axis,
        n_q=int(args.oe5_n_q_thermal),
        n_matsubara=int(args.oe5_n_matsubara_thermal),
        stable_lowT_branch_only=True,
        workers=workers,
    )
    thermal_grid_npz = save_thermal_usadel_grid_npz(
        thermal_grid,
        raw_pre / "thermal_usadel_grid.npz",
    )
    thermal_grid_summary = thermal_usadel_grid_summary(thermal_grid)

    if args.oe5_state_current_fraction is None:
        if np.isfinite(Ic_A) and Ic_A > 0.0:
            target_fraction = float(np.clip(I_bias_A / Ic_A, 0.0, 1.0))
        else:
            target_fraction = 0.90
    else:
        target_fraction = float(args.oe5_state_current_fraction)

    state = select_thermal_usadel_q_state(thermal_grid, target_fraction)

    power_curve = compute_power_curve_thermal_usadel_state(
        Te_values,
        Tph_K=Tph_K,
        state=state,
        thermal_grid=thermal_grid,
        phase_space_catalog=phase_catalog,
        spectrum=spectrum,
        N0_J_m3=N0,
        tau0_s=tau0_s,
        Tc_K=Tc_K,
        omega_max_meV=args.oe5_omega_max_meV,
    )

    q_scan = compute_power_scan_thermal_usadel(
        q_scan_Te,
        Tph_K=Tph_K,
        thermal_grid=thermal_grid,
        phase_space_catalog=phase_catalog,
        spectrum=spectrum,
        N0_J_m3=N0,
        omega_max_meV=args.oe5_omega_max_meV,
    )

    support_Te, support_delta = _select_support_state_below_gap_cutoff(
        thermal_grid,
        state,
        Te_values,
        min_delta_fraction=float(args.oe5_support_min_delta_fraction),
    )

    support_result = compute_projected_powers(
        support_Te,
        Tph_K,
        support_delta,
        float(state["q_m_inv"]),
        phase_catalog,
        spectrum,
        N0_J_m3=N0,
        omega_max_meV=args.oe5_omega_max_meV,
    )
    support = cumulative_spectral_support(support_result)

    power_npz = raw_pre / "oe5_power_catalog.npz"
    np.savez_compressed(
        power_npz,
        Te_values_K=power_curve["Te_values_K"],
        state_delta_values_J=power_curve["delta_values_J"],
        state_q_values_m_inv=power_curve["q_values_m_inv"],
        state_P_S_W_m3=power_curve["P_S_W_m3"],
        state_P_R_W_m3=power_curve["P_R_W_m3"],
        state_P_total_W_m3=power_curve["P_total_W_m3"],
        state_P_normal_Eliashberg_W_m3=power_curve["P_normal_Eliashberg_W_m3"],
        state_P_Debye_Vodolazov_W_m3=power_curve["P_Debye_Vodolazov_W_m3"],
        q_scan_Te_values_K=q_scan["Te_values_K"],
        q_scan_q_values_m_inv=q_scan["q_values_m_inv"],
        q_scan_reference_current_fraction=q_scan["reference_current_fraction"],
        q_scan_delta_values_J=q_scan["delta_values_J"],
        q_scan_P_S_W_m3=q_scan["P_S_W_m3"],
        q_scan_P_R_W_m3=q_scan["P_R_W_m3"],
        q_scan_P_total_W_m3=q_scan["P_total_W_m3"],
        support_omega_meV=support["omega_meV"],
        cumulative_alpha_omega=support["cumulative_alpha_omega"],
        cumulative_scattering=support["cumulative_scattering"],
        cumulative_recombination=support["cumulative_recombination"],
        support_integrand_S_J2=support_result.integrand_S_J2,
        support_integrand_R_J2=support_result.integrand_R_J2,
    )

    spectrum_plot = plot_eliashberg_spectrum(
        spectrum,
        plots_diagnostics / "eliashberg_spectrum.png",
    )
    thermal_gap_plot = plot_thermal_usadel_gap_grid(
        thermal_grid,
        plots_diagnostics / "thermal_usadel_gap_grid.png",
        target_fraction=state["reference_current_fraction"],
    )
    power_plot = plot_power_curve_thermal_usadel(
        power_curve,
        plots_comparisons / "electron_phonon_power_vs_Te_thermal_usadel.png",
        tau_label=r"$\tau_0$ from $\tau_{ep}(T_c)$",
        title_suffix=(
            rf"thermal Usadel state: "
            rf"$I/I_c^{{bias}}={state['reference_current_fraction']:.3f}$"
        ),
    )
    ratio_plot = plot_power_ratios_thermal_usadel(
        power_curve,
        plots_comparisons / "electron_phonon_power_ratios_vs_Te.png",
    )
    scan_plot = plot_power_scan_thermal_usadel(
        q_scan,
        plots_comparisons / "electron_phonon_power_vs_thermal_usadel_current.png",
    )
    support_plot = plot_spectral_support(
        support,
        plots_comparisons / "spectral_support_thermal_usadel_state.png",
    )
    low_energy_plot = plot_low_energy_recombination_scattering_band(
        support_result,
        plots_comparisons / "low_energy_recombination_scattering_band.png",
        omega_max_meV=args.oe5_low_energy_max_meV,
    )

    summary = {
        "backend": "oe5_prerun_projected_power_v5_integrated",
        "inputs": {
            "eliashberg_dat": str(eliashberg_path),
        },
        "material": {
            "sigma_n_S_m": sigma_n,
            "D_m2_s": D_m2_s,
            "N0_J_m3": N0,
            "Tc_K": Tc_K,
            "T_bias_K": T_bias_K,
            "I_bias_A": I_bias_A,
            "Ic_A": Ic_A,
            "tau_source": tau_source,
            "tau_ep_Tc_ps": float(tau_ep_Tc_s / 1.0e-12),
            "tau0_ps": float(tau0_s / 1.0e-12),
            "tau0_over_tau_ep_Tc": float(TAU0_OVER_TAU_EP_TC),
        },
        "thermal_usadel_grid": thermal_grid_summary,
        "selected_thermal_usadel_state": state,
        "support_state": {
            "Te_K": float(support_Te),
            "Tph_K": float(Tph_K),
            "delta_meV": float(support_delta / MEV_J),
            "q_m_inv": float(state["q_m_inv"]),
            "P_S_W_m3": support_result.P_S_W_m3,
            "P_R_W_m3": support_result.P_R_W_m3,
            "P_total_W_m3": support_result.P_total_W_m3,
            "policy": (
                "Chosen as the highest Te on the OE5 axis with Delta_eq still "
                "above the configured min_delta_fraction. This keeps the "
                "spectral-support and low-energy-band plots below the local "
                "thermal superconducting cutoff."
            ),
        },
        "power_axis": {
            "Tph_K": float(Tph_K),
            "Te_min_K": float(Te_values[0]),
            "Te_max_K": float(Te_values[-1]),
            "n_Te": int(Te_values.size),
            "q_scan_Te_values_K": [float(v) for v in q_scan_Te],
            "omega_max_meV": float(args.oe5_omega_max_meV),
        },
        "eliashberg": spectrum_summary(spectrum),
        "curve_extrema": {
            "state_P_S_max_W_m3": float(np.max(power_curve["P_S_W_m3"])),
            "state_P_R_max_W_m3": float(np.max(power_curve["P_R_W_m3"])),
            "state_P_total_max_W_m3": float(np.max(power_curve["P_total_W_m3"])),
            "state_P_normal_Eliashberg_max_W_m3": float(
                np.max(power_curve["P_normal_Eliashberg_W_m3"])
            ),
            "state_P_Debye_max_W_m3": float(
                np.max(power_curve["P_Debye_Vodolazov_W_m3"])
            ),
            "q_scan_P_S_max_W_m3": float(np.max(q_scan["P_S_W_m3"])),
            "q_scan_P_R_max_W_m3": float(np.max(q_scan["P_R_W_m3"])),
            "q_scan_P_total_max_W_m3": float(np.max(q_scan["P_total_W_m3"])),
        },
        "interpretation": {
            "pre_run_policy": (
                "OE5 projected-power diagnostics are now generated during "
                "01_prerun_template.py. The old 02_oe5_power_debug.py pipeline "
                "is obsolete."
            ),
            "parallel_policy": (
                "PRE-run is the correct stage for parallel catalogue work. "
                "SS/PHOTON time evolution should load these catalogues rather "
                "than recomputing them."
            ),
            "recombination_policy": (
                "P_R is a superconducting recombination/pair-breaking channel. "
                "It is set to zero when Delta_eq(Te,q)=0 and must not be treated "
                "as an additional normal-state power."
            ),
        },
        "outputs": {
            "thermal_usadel_grid_npz": str(thermal_grid_npz),
            "power_npz": str(power_npz),
            "power_summary": str(raw_pre / "oe5_power_summary.yaml"),
            "eliashberg_plot": str(spectrum_plot),
            "thermal_gap_plot": str(thermal_gap_plot),
            "power_plot": str(power_plot),
            "ratio_plot": str(ratio_plot),
            "scan_plot": str(scan_plot),
            "support_plot": str(support_plot),
            "low_energy_plot": str(low_energy_plot),
        },
    }

    summary_path = raw_pre / "oe5_power_summary.yaml"
    with summary_path.open("w", encoding="utf-8") as f:
        yaml.safe_dump(summary, f, sort_keys=False, allow_unicode=True)

    outputs = {
        "thermal_usadel_grid_npz": str(thermal_grid_npz),
        "oe5_power_npz": str(power_npz),
        "oe5_power_summary": str(summary_path),
        "eliashberg_plot": str(spectrum_plot),
        "thermal_gap_plot": str(thermal_gap_plot),
        "oe5_power_plot": str(power_plot),
        "oe5_ratio_plot": str(ratio_plot),
        "oe5_scan_plot": str(scan_plot),
        "oe5_support_plot": str(support_plot),
        "oe5_low_energy_plot": str(low_energy_plot),
    }

    return summary, outputs


def _resolve_eliashberg_path(cfg: dict[str, Any], user_path: str | None) -> Path:
    if user_path is not None:
        path = Path(user_path).expanduser().resolve()
    else:
        path = (
            Path(cfg["project"]["big_data_root"]).expanduser().resolve()
            / "catalogs"
            / "simon_2025"
            / "nbn-a2f-ph.dat"
        )

    if not path.exists():
        raise FileNotFoundError(
            "OE5 requires the Simon/MIT Eliashberg data file. Expected: "
            f"{path}. Provide --eliashberg-dat or copy the file under "
            "big_data_root/catalogs/simon_2025/."
        )

    return path


def _select_support_state_below_gap_cutoff(
    thermal_grid,
    state: dict[str, float],
    Te_values: np.ndarray,
    *,
    min_delta_fraction: float,
) -> tuple[float, float]:
    q = float(state["q_m_inv"])

    deltas = np.array(
        [
            thermal_usadel_delta_at_state(thermal_grid, Te_K=float(Te), q_m_inv=q)
            for Te in Te_values
        ],
        dtype=float,
    )

    delta_max = float(np.max(deltas))
    if delta_max <= 0.0:
        idx = int(np.argmin(np.abs(Te_values - 0.8 * float(thermal_grid.metadata["Tc_K"]))))
        return float(Te_values[idx]), float(deltas[idx])

    threshold = max(0.0, float(min_delta_fraction)) * delta_max
    valid = np.where(deltas > threshold)[0]

    if valid.size == 0:
        idx = int(np.argmax(deltas))
    else:
        idx = int(valid[-1])

    return float(Te_values[idx]), float(deltas[idx])


def _print_dict(data: dict[str, Any]) -> None:
    for key, value in data.items():
        print(f"  {key}: {value}")


if __name__ == "__main__":
    raise SystemExit(main())