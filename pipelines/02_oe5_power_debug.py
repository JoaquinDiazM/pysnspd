"""OE5 debug pipeline: projected electron-phonon powers.

This final OE5 debug pipeline consumes:

1. An OE4 phase-space catalogue:
       phase_space_catalog.npz

2. The parent OE3 Usadel catalogue:
       usadel_dos_catalog.npz

3. A Simon et al. 2025 Eliashberg/PhDOS data file:
       nbn-a2f-ph.dat

Unlike the intermediate OE5 attempts, this version does not plot a fixed,
non-self-consistent gap trajectory. Instead, it extracts the stable branch of
the OE3 Usadel calibration sweep,

    Delta_eq(q,T_bias),

and evaluates the projected powers along that self-consistent branch.
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
    build_usadel_self_consistent_trajectory,
    compute_power_curve_at_usadel_state,
    compute_projected_powers,
    compute_usadel_q_power_scan,
    cumulative_spectral_support,
    electronic_density_of_states_from_sigma_D,
    select_usadel_state_by_current_fraction,
    tau0_from_tau_ep_Tc,
    tau_ep_Tc_from_tau0,
)
from pysnspd.plotting.kinetic import (
    plot_eliashberg_spectrum,
    plot_power_curve,
    plot_power_vs_usadel_current,
    plot_spectral_support,
    plot_usadel_self_consistent_trajectory,
)
from pysnspd.usadel.catalog import load_usadel_catalog_npz


MEV_J = 1.602176634e-22


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="OE5 projected power debug using Usadel self-consistent branch."
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
    parser.add_argument(
        "--Te-max-K",
        type=float,
        default=None,
        help=(
            "Maximum Te for power-vs-Te diagnostic. If omitted, the pipeline "
            "uses 0.95 Tc to avoid interpreting a low-temperature Usadel "
            "self-consistent gap above Tc."
        ),
    )
    parser.add_argument("--n-Te", type=int, default=100)

    parser.add_argument(
        "--q-scan-Te-K",
        type=float,
        nargs="*",
        default=None,
        help=(
            "Electron temperatures used for the self-consistent Usadel q-scan. "
            "If omitted, uses Tph, 0.5Tc and 0.8Tc."
        ),
    )
    parser.add_argument(
        "--n-q-trajectory",
        type=int,
        default=120,
        help="Number of q points in the interpolated stable Usadel branch.",
    )
    parser.add_argument(
        "--state-current-fraction",
        type=float,
        default=None,
        help=(
            "Select representative Usadel state by I/Ic. If omitted, use the "
            "configured bias current divided by the calibrated Ic."
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
    N0 = electronic_density_of_states_from_sigma_D(sigma_n, D_m2_s)

    if args.tau0_ps is not None:
        tau0_s = float(args.tau0_ps) * 1.0e-12
        tau_ep_Tc_s = tau_ep_Tc_from_tau0(tau0_s)
        tau_source = "direct_tau0"
    else:
        tau_ep_Tc_s = float(args.tau_ep_Tc_ps) * 1.0e-12
        tau0_s = tau0_from_tau_ep_Tc(tau_ep_Tc_s)
        tau_source = "converted_from_tau_ep_Tc"

    trajectory = build_usadel_self_consistent_trajectory(
        usadel_catalog,
        n_q=int(args.n_q_trajectory),
        stable_branch_only=True,
    )

    Ic_A = float(trajectory.metadata["Ic_A"])
    if args.state_current_fraction is None:
        I_bias_A = float(usadel_catalog.metadata.get("I_bias_A", 0.0))
        target_fraction = I_bias_A / Ic_A if Ic_A > 0.0 else 0.80
        target_fraction = float(np.clip(target_fraction, 0.0, 1.0))
    else:
        target_fraction = float(args.state_current_fraction)

    state = select_usadel_state_by_current_fraction(trajectory, target_fraction)

    Te_max = float(args.Te_max_K) if args.Te_max_K is not None else 0.95 * Tc_K
    Te_values = np.linspace(float(args.Te_min_K), Te_max, int(args.n_Te))

    if args.q_scan_Te_K:
        q_scan_Te = np.asarray(args.q_scan_Te_K, dtype=float)
    else:
        q_scan_Te = np.asarray(
            [
                float(args.Tph_K),
                0.5 * Tc_K,
                0.8 * Tc_K,
            ],
            dtype=float,
        )
    q_scan_Te = np.unique(np.clip(q_scan_Te, float(args.Tph_K), Te_max))

    power_curve = compute_power_curve_at_usadel_state(
        Te_values,
        Tph_K=float(args.Tph_K),
        state=state,
        phase_space_catalog=phase_catalog,
        spectrum=spectrum,
        N0_J_m3=N0,
        tau0_s=tau0_s,
        Tc_K=Tc_K,
        omega_max_meV=args.omega_max_meV,
    )

    q_scan = compute_usadel_q_power_scan(
        q_scan_Te,
        Tph_K=float(args.Tph_K),
        trajectory=trajectory,
        phase_space_catalog=phase_catalog,
        spectrum=spectrum,
        N0_J_m3=N0,
        omega_max_meV=args.omega_max_meV,
    )

    representative = compute_projected_powers(
        float(q_scan_Te[-1]),
        float(args.Tph_K),
        float(state["delta_J"]),
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
    trajectory_plot = plot_usadel_self_consistent_trajectory(
        trajectory,
        plots_diag / "usadel_self_consistent_trajectory.png",
    )
    power_plot = plot_power_curve(
        power_curve,
        plots_comp / "electron_phonon_power_vs_Te_usadel_state.png",
        tau_label=r"$\tau_0$ from $\tau_{ep}(T_c)$",
        title_suffix=(
            rf"Usadel state: $I/I_c={state['current_fraction']:.3f}$, "
            rf"$\Delta={state['delta_meV']:.3f}$ meV"
        ),
    )
    q_power_plot = plot_power_vs_usadel_current(
        q_scan,
        plots_comp / "electron_phonon_power_vs_usadel_current.png",
    )
    support_plot = plot_spectral_support(
        support,
        plots_comp / "spectral_support_usadel_state.png",
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
        state_P_Debye_Vodolazov_W_m3=power_curve["P_Debye_Vodolazov_W_m3"],
        q_scan_Te_values_K=q_scan["Te_values_K"],
        q_scan_q_values_m_inv=q_scan["q_values_m_inv"],
        q_scan_gamma_values_J=q_scan["gamma_values_J"],
        q_scan_delta_values_J=q_scan["delta_values_J"],
        q_scan_current_values_A=q_scan["current_values_A"],
        q_scan_current_fraction=q_scan["current_fraction"],
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
        "backend": "oe5_projected_power_debug_v3_usadel_self_consistent",
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
            "tau_source": tau_source,
            "tau_ep_Tc_ps": float(tau_ep_Tc_s / 1.0e-12),
            "tau0_ps": float(tau0_s / 1.0e-12),
            "tau0_over_tau_ep_Tc": float(TAU0_OVER_TAU_EP_TC),
        },
        "usadel_trajectory": trajectory.metadata,
        "selected_usadel_state": state,
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
            "Te_K": float(q_scan_Te[-1]),
            "Tph_K": float(args.Tph_K),
            "delta_meV": float(state["delta_meV"]),
            "q_m_inv": float(state["q_m_inv"]),
            "P_S_W_m3": representative.P_S_W_m3,
            "P_R_W_m3": representative.P_R_W_m3,
            "P_total_W_m3": representative.P_total_W_m3,
            "omega_max_meV_used": representative.metadata["omega_max_meV_used"],
        },
        "curve_extrema": {
            "state_P_S_min_W_m3": float(np.min(power_curve["P_S_W_m3"])),
            "state_P_S_max_W_m3": float(np.max(power_curve["P_S_W_m3"])),
            "state_P_R_min_W_m3": float(np.min(power_curve["P_R_W_m3"])),
            "state_P_R_max_W_m3": float(np.max(power_curve["P_R_W_m3"])),
            "state_P_total_min_W_m3": float(np.min(power_curve["P_total_W_m3"])),
            "state_P_total_max_W_m3": float(np.max(power_curve["P_total_W_m3"])),
            "state_P_Debye_max_W_m3": float(
                np.max(power_curve["P_Debye_Vodolazov_W_m3"])
            ),
            "q_scan_P_S_max_W_m3": float(np.max(q_scan["P_S_W_m3"])),
            "q_scan_P_R_max_W_m3": float(np.max(q_scan["P_R_W_m3"])),
            "q_scan_P_total_max_W_m3": float(np.max(q_scan["P_total_W_m3"])),
        },
        "interpretation": {
            "self_consistency_policy": (
                "OE5-v3 uses Delta_eq(q,T_bias) from the OE3 Matsubara Usadel "
                "calibration sweep. No artificial fixed-gap or BCS-like diagnostic "
                "trajectory is plotted."
            ),
            "future_coupled_policy": (
                "In the coupled gTDGL/thermal model, Delta(r,t) and q(r,t) will "
                "replace this static Usadel calibration trajectory."
            ),
        },
        "outputs": {
            "power_npz": str(power_npz),
            "power_summary": str(raw_ss / "oe5_power_summary.yaml"),
            "eliashberg_plot": str(spectrum_plot),
            "trajectory_plot": str(trajectory_plot),
            "power_plot": str(power_plot),
            "q_power_plot": str(q_power_plot),
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
                "OE5 local projected electron-phonon powers evaluated on the "
                "OE3 Usadel self-consistent stable branch."
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
    print(f"  tau_source              : {tau_source}")
    print(f"  tau_ep_Tc_ps            : {tau_ep_Tc_s / 1.0e-12}")
    print(f"  tau0_ps                 : {tau0_s / 1.0e-12}")
    print()
    print("Selected Usadel state")
    for key, value in state.items():
        print(f"  {key}: {value}")
    print()
    print("Representative projected power")
    print(f"  Te_K                    : {q_scan_Te[-1]}")
    print(f"  Tph_K                   : {args.Tph_K}")
    print(f"  delta_meV               : {state['delta_meV']}")
    print(f"  q_m_inv                 : {state['q_m_inv']}")
    print(f"  P_S_W_m3                : {representative.P_S_W_m3}")
    print(f"  P_R_W_m3                : {representative.P_R_W_m3}")
    print(f"  P_total_W_m3            : {representative.P_total_W_m3}")
    print()
    print("Outputs")
    print(f"  power_npz               : {power_npz}")
    print(f"  power_summary           : {summary_path}")
    print(f"  eliashberg_plot         : {spectrum_plot}")
    print(f"  trajectory_plot         : {trajectory_plot}")
    print(f"  power_plot              : {power_plot}")
    print(f"  q_power_plot            : {q_power_plot}")
    print(f"  support_plot            : {support_plot}")
    print(f"  manifest                : {manifest}")
    print("Status: OK")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())