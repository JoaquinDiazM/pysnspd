"""OE5 debug pipeline: thermal-Usadel projected electron-phonon powers.

This final OE5 audit consumes:

1. An OE4 phase-space catalogue:
       phase_space_catalog.npz

2. The parent OE3 Usadel catalogue:
       usadel_dos_catalog.npz

3. A Simon et al. 2025 Eliashberg/PhDOS data file:
       nbn-a2f-ph.dat

It constructs a thermal Usadel equilibrium grid,

    Delta_eq(Te,q),

using the same Matsubara self-consistency solver as OE3. The projected powers
are then evaluated as

    P_ep(Te; Delta_eq(Te,q), q).

This removes the earlier fixed-gap and BCS-like diagnostic plots.
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Any

import numpy as np
import yaml

from pysnspd.config import load_config, validate_config
from pysnspd.io.manager import create_run_layout, write_manifest
from pysnspd.kinetic.eliashberg import (
    load_simon_eliashberg_dat,
    spectrum_summary,
)
from pysnspd.kinetic.phase_space import load_phase_space_catalog_npz
from pysnspd.kinetic.powers import (
    TAU0_OVER_TAU_EP_TC,
    build_thermal_usadel_grid,
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
from pysnspd.plotting.kinetic import (
    plot_eliashberg_spectrum,
    plot_power_curve_thermal_usadel,
    plot_power_ratios_thermal_usadel,
    plot_power_scan_thermal_usadel,
    plot_spectral_support,
    plot_thermal_usadel_gap_grid,
)
from pysnspd.usadel.catalog import load_usadel_catalog_npz


MEV_J = 1.602176634e-22


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="OE5 projected powers with thermal Usadel self-consistency."
    )

    parser.add_argument("--config", required=True, help="Path to YAML config.")
    parser.add_argument("--run-name", required=True, help="OE5 debug run name.")

    parser.add_argument(
        "--phase-space-npz",
        required=True,
        help="Path to OE4 phase_space_catalog.npz.",
    )
    parser.add_argument(
        "--usadel-npz",
        required=True,
        help="Path to parent OE3 usadel_dos_catalog.npz.",
    )
    parser.add_argument(
        "--eliashberg-dat",
        required=True,
        help="Path to Simon/MIT nbn-a2f-ph.dat.",
    )

    parser.add_argument("--Tph-K", type=float, default=0.9)
    parser.add_argument("--Te-min-K", type=float, default=0.9)
    parser.add_argument("--Te-max-K", type=float, default=34.6)
    parser.add_argument("--n-Te", type=int, default=160)

    parser.add_argument(
        "--q-scan-Te-K",
        type=float,
        nargs="*",
        default=None,
        help=(
            "Electron temperatures used for the thermal Usadel q-scan. "
            "If omitted, uses Tph, 0.8Tc, Tc and 2Tc clipped to Te range."
        ),
    )
    parser.add_argument(
        "--n-q-thermal",
        type=int,
        default=120,
        help="Number of q points in the thermal Usadel grid.",
    )
    parser.add_argument(
        "--n-matsubara-thermal",
        type=int,
        default=500,
        help="Number of Matsubara frequencies for Delta_eq(Te,q).",
    )
    parser.add_argument(
        "--state-current-fraction",
        type=float,
        default=None,
        help=(
            "Select representative q by reference I(q,T_bias)/Ic(T_bias). "
            "If omitted, use configured I_bias/Ic from the OE3 metadata."
        ),
    )
    parser.add_argument(
        "--omega-max-meV",
        type=float,
        default=None,
        help="Optional maximum Omega used for projected power integrals.",
    )

    tau_group = parser.add_mutually_exclusive_group()
    tau_group.add_argument(
        "--tau-ep-Tc-ps",
        type=float,
        default=24.7,
        help=(
            "Linear electron-phonon relaxation time at Tc in ps. The pipeline "
            "converts it to Vodolazov tau0."
        ),
    )
    tau_group.add_argument(
        "--tau0-ps",
        type=float,
        default=None,
        help="Direct Vodolazov/Allmaras tau0 in ps.",
    )

    return parser.parse_args()


def main() -> int:
    args = parse_args()

    cfg = validate_config(load_config(args.config))
    layout = create_run_layout(cfg, args.run_name)

    raw_ss = Path(layout["raw_ss"])
    plots_diag = Path(layout["plots_diagnostics"])
    plots_comp = Path(layout["plots_comparisons"])

    phase_path = Path(args.phase_space_npz).expanduser().resolve()
    usadel_path = Path(args.usadel_npz).expanduser().resolve()
    eliashberg_path = Path(args.eliashberg_dat).expanduser().resolve()

    phase_catalog = load_phase_space_catalog_npz(phase_path)
    usadel_catalog = load_usadel_catalog_npz(usadel_path)
    spectrum = load_simon_eliashberg_dat(eliashberg_path)

    D_m2_s = float(usadel_catalog.metadata["D_m2_s"])
    sigma_n = float(usadel_catalog.metadata["sigma_n_S_m"])
    Tc_K = float(usadel_catalog.metadata["Tc_K"])
    T_bias_K = float(usadel_catalog.metadata["T_bias_K"])
    I_bias_A = float(usadel_catalog.metadata.get("I_bias_A", 0.0))

    calibration = usadel_catalog.metadata.get("calibration", {})
    Ic_A = float(calibration.get("Ic_model_A", calibration.get("Ic_target_A", np.nan)))

    N0 = electronic_density_of_states_from_sigma_D(sigma_n, D_m2_s)

    if args.tau0_ps is not None:
        tau0_s = float(args.tau0_ps) * 1.0e-12
        tau_ep_Tc_s = tau_ep_Tc_from_tau0(tau0_s)
        tau_source = "direct_tau0"
    else:
        tau_ep_Tc_s = float(args.tau_ep_Tc_ps) * 1.0e-12
        tau0_s = tau0_from_tau_ep_Tc(tau_ep_Tc_s)
        tau_source = "converted_from_tau_ep_Tc"

    Te_values = np.linspace(float(args.Te_min_K), float(args.Te_max_K), int(args.n_Te))

    if args.q_scan_Te_K:
        q_scan_Te = np.asarray(args.q_scan_Te_K, dtype=float)
    else:
        q_scan_Te = np.asarray(
            [
                float(args.Tph_K),
                0.8 * Tc_K,
                Tc_K,
                min(2.0 * Tc_K, float(args.Te_max_K)),
            ],
            dtype=float,
        )
    q_scan_Te = np.unique(
        np.clip(q_scan_Te, float(args.Te_min_K), float(args.Te_max_K))
    )

    thermal_Te_axis = np.unique(np.concatenate([Te_values, q_scan_Te]))

    thermal_grid = build_thermal_usadel_grid(
        usadel_catalog,
        thermal_Te_axis,
        n_q=int(args.n_q_thermal),
        n_matsubara=int(args.n_matsubara_thermal),
        stable_lowT_branch_only=True,
    )

    if args.state_current_fraction is None:
        if np.isfinite(Ic_A) and Ic_A > 0.0:
            target_fraction = float(np.clip(I_bias_A / Ic_A, 0.0, 1.0))
        else:
            target_fraction = 0.90
    else:
        target_fraction = float(args.state_current_fraction)

    state = select_thermal_usadel_q_state(thermal_grid, target_fraction)

    power_curve = compute_power_curve_thermal_usadel_state(
        Te_values,
        Tph_K=float(args.Tph_K),
        state=state,
        thermal_grid=thermal_grid,
        phase_space_catalog=phase_catalog,
        spectrum=spectrum,
        N0_J_m3=N0,
        tau0_s=tau0_s,
        Tc_K=Tc_K,
        omega_max_meV=args.omega_max_meV,
    )

    q_scan = compute_power_scan_thermal_usadel(
        q_scan_Te,
        Tph_K=float(args.Tph_K),
        thermal_grid=thermal_grid,
        phase_space_catalog=phase_catalog,
        spectrum=spectrum,
        N0_J_m3=N0,
        omega_max_meV=args.omega_max_meV,
    )

    representative_Te = float(q_scan_Te[-1])
    representative_delta = thermal_usadel_delta_at_state(
        thermal_grid,
        Te_K=representative_Te,
        q_m_inv=float(state["q_m_inv"]),
    )
    representative = compute_projected_powers(
        representative_Te,
        float(args.Tph_K),
        representative_delta,
        float(state["q_m_inv"]),
        phase_catalog,
        spectrum,
        N0_J_m3=N0,
        omega_max_meV=args.omega_max_meV,
    )
    support = cumulative_spectral_support(representative)

    spectrum_plot = plot_eliashberg_spectrum(
        spectrum,
        plots_diag / "eliashberg_spectrum.png",
    )
    thermal_gap_plot = plot_thermal_usadel_gap_grid(
        thermal_grid,
        plots_diag / "thermal_usadel_gap_grid.png",
        target_fraction=state["reference_current_fraction"],
    )
    power_plot = plot_power_curve_thermal_usadel(
        power_curve,
        plots_comp / "electron_phonon_power_vs_Te_thermal_usadel.png",
        tau_label=r"$\tau_0$ from $\tau_{ep}(T_c)$",
        title_suffix=(
            rf"thermal Usadel state: "
            rf"$I/I_c^{{bias}}={state['reference_current_fraction']:.3f}$"
        ),
    )
    ratio_plot = plot_power_ratios_thermal_usadel(
        power_curve,
        plots_comp / "electron_phonon_power_ratios_vs_Te.png",
    )
    scan_plot = plot_power_scan_thermal_usadel(
        q_scan,
        plots_comp / "electron_phonon_power_vs_thermal_usadel_current.png",
    )
    support_plot = plot_spectral_support(
        support,
        plots_comp / "spectral_support_thermal_usadel_state.png",
    )

    power_npz = raw_ss / "oe5_power_catalog.npz"
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
        thermal_grid_Te_values_K=thermal_grid.Te_values_K,
        thermal_grid_q_values_m_inv=thermal_grid.q_values_m_inv,
        thermal_grid_gamma_values_J=thermal_grid.gamma_values_J,
        thermal_grid_delta_eq_Tq_J=thermal_grid.delta_eq_Tq_J,
        thermal_grid_current_Tq_A=thermal_grid.current_Tq_A,
        thermal_grid_reference_current_fraction=thermal_grid.reference_current_fraction,
        q_scan_Te_values_K=q_scan["Te_values_K"],
        q_scan_q_values_m_inv=q_scan["q_values_m_inv"],
        q_scan_reference_current_fraction=q_scan["reference_current_fraction"],
        q_scan_delta_values_J=q_scan["delta_values_J"],
        q_scan_P_S_W_m3=q_scan["P_S_W_m3"],
        q_scan_P_R_W_m3=q_scan["P_R_W_m3"],
        q_scan_P_total_W_m3=q_scan["P_total_W_m3"],
        representative_omega_J=representative.omega_J,
        representative_alpha2F=representative.alpha2F,
        representative_bose_difference=representative.bose_difference,
        representative_J_S_J=representative.J_S_J,
        representative_J_R_J=representative.J_R_J,
        representative_integrand_S_J2=representative.integrand_S_J2,
        representative_integrand_R_J2=representative.integrand_R_J2,
        support_omega_meV=support["omega_meV"],
        cumulative_alpha_omega=support["cumulative_alpha_omega"],
        cumulative_scattering=support["cumulative_scattering"],
        cumulative_recombination=support["cumulative_recombination"],
    )

    summary = {
        "backend": "oe5_projected_power_debug_v4_thermal_usadel",
        "inputs": {
            "phase_space_npz": str(phase_path),
            "usadel_npz": str(usadel_path),
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
        "thermal_usadel_grid": {
            **thermal_grid.metadata,
            "Te_min_K": float(thermal_grid.Te_values_K[0]),
            "Te_max_K": float(thermal_grid.Te_values_K[-1]),
            "delta_max_meV": float(np.max(thermal_grid.delta_eq_Tq_J) / MEV_J),
            "delta_min_meV": float(np.min(thermal_grid.delta_eq_Tq_J) / MEV_J),
        },
        "selected_thermal_usadel_state": state,
        "power_axis": {
            "Tph_K": float(args.Tph_K),
            "Te_min_K": float(Te_values[0]),
            "Te_max_K": float(Te_values[-1]),
            "n_Te": int(Te_values.size),
            "q_scan_Te_values_K": [float(v) for v in q_scan_Te],
            "omega_max_meV_requested": args.omega_max_meV,
        },
        "eliashberg": spectrum_summary(spectrum),
        "representative_power": {
            "Te_K": representative_Te,
            "Tph_K": float(args.Tph_K),
            "delta_meV": float(representative_delta / MEV_J),
            "q_m_inv": float(state["q_m_inv"]),
            "P_S_W_m3": representative.P_S_W_m3,
            "P_R_W_m3": representative.P_R_W_m3,
            "P_total_W_m3": representative.P_total_W_m3,
            "omega_max_meV_used": representative.metadata["omega_max_meV_used"],
        },
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
            "self_consistency_policy": (
                "OE5-v4 uses Delta_eq(Te,q) from the Matsubara Usadel gap "
                "equation. Fixed-gap and BCS-like diagnostic plots were removed."
            ),
            "normal_reference_policy": (
                "The normal Eliashberg and Debye curves are references for the "
                "Delta -> 0 limit. They are not replacements for the "
                "superconducting projected powers."
            ),
            "future_coupled_policy": (
                "In PHOTON-runs, Delta(r,t) and q(r,t) from gTDGL will replace "
                "this local thermal-equilibrium Usadel grid."
            ),
        },
        "outputs": {
            "power_npz": str(power_npz),
            "power_summary": str(raw_ss / "oe5_power_summary.yaml"),
            "eliashberg_plot": str(spectrum_plot),
            "thermal_gap_plot": str(thermal_gap_plot),
            "power_plot": str(power_plot),
            "ratio_plot": str(ratio_plot),
            "scan_plot": str(scan_plot),
            "support_plot": str(support_plot),
        },
    }

    summary_path = raw_ss / "oe5_power_summary.yaml"
    with summary_path.open("w", encoding="utf-8") as f:
        yaml.safe_dump(summary, f, sort_keys=False, allow_unicode=True)

    manifest = write_manifest(
        cfg,
        args.run_name,
        stage="ss",
        extra={
            "pipeline": "02_oe5_power_debug.py",
            "purpose": (
                "OE5 local projected electron-phonon powers evaluated with "
                "thermal Usadel self-consistency Delta_eq(Te,q)."
            ),
            "summary": summary,
        },
    )

    print("OE5 projected-power debug")
    print(f"run_name                  : {args.run_name}")
    print(f"phase_space_npz           : {phase_path}")
    print(f"usadel_npz                : {usadel_path}")
    print(f"eliashberg_dat            : {eliashberg_path}")
    print(f"raw_ss                    : {raw_ss}")
    print(f"plots_diagnostics         : {plots_diag}")
    print(f"plots_comparisons         : {plots_comp}")
    print()
    print("Material")
    print(f"  D_m2_s                  : {D_m2_s}")
    print(f"  sigma_n_S_m             : {sigma_n}")
    print(f"  N0_J_m3                 : {N0}")
    print(f"  Tc_K                    : {Tc_K}")
    print(f"  T_bias_K                : {T_bias_K}")
    print(f"  I_bias_A                : {I_bias_A}")
    print(f"  Ic_A                    : {Ic_A}")
    print(f"  tau_source              : {tau_source}")
    print(f"  tau_ep_Tc_ps            : {tau_ep_Tc_s / 1.0e-12}")
    print(f"  tau0_ps                 : {tau0_s / 1.0e-12}")
    print()
    print("Selected thermal Usadel q-state")
    for key, value in state.items():
        print(f"  {key}: {value}")
    print()
    print("Representative projected power")
    print(f"  Te_K                    : {representative_Te}")
    print(f"  Tph_K                   : {args.Tph_K}")
    print(f"  delta_meV               : {representative_delta / MEV_J}")
    print(f"  q_m_inv                 : {state['q_m_inv']}")
    print(f"  P_S_W_m3                : {representative.P_S_W_m3}")
    print(f"  P_R_W_m3                : {representative.P_R_W_m3}")
    print(f"  P_total_W_m3            : {representative.P_total_W_m3}")
    print()
    print("Outputs")
    print(f"  power_npz               : {power_npz}")
    print(f"  power_summary           : {summary_path}")
    print(f"  eliashberg_plot         : {spectrum_plot}")
    print(f"  thermal_gap_plot        : {thermal_gap_plot}")
    print(f"  power_plot              : {power_plot}")
    print(f"  ratio_plot              : {ratio_plot}")
    print(f"  scan_plot               : {scan_plot}")
    print(f"  support_plot            : {support_plot}")
    print(f"  manifest                : {manifest}")
    print("Status: OK")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())