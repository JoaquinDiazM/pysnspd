"""OE5 debug pipeline: projected electron-phonon powers.

This pipeline consumes:

1. An OE4 phase-space catalogue:
       phase_space_catalog.npz

2. A Simon et al. 2025 Eliashberg/PhDOS data file:
       nbn-a2f-ph.dat

and produces local 0D projected powers

    P_ep^S(Te,Tph,Delta,q),
    P_ep^R(Te,Tph,Delta,q),

plus spectral-support diagnostics.

This is not yet the thermal evolution solver. It is the OE5 validation layer
that checks signs, magnitudes, spectral support and comparison against the
Vodolazov/Allmaras Debye T^5 scale.

Important:
    The Debye/Vodolazov tau0 is not tau_ep(Tc). By default this pipeline takes
    tau_ep(Tc) and converts it internally into tau0.
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
    compute_power_curve,
    compute_projected_powers,
    cumulative_spectral_support,
    diagnostic_bcs_gap_factor,
    electronic_density_of_states_from_sigma_D,
    tau0_from_tau_ep_Tc,
    tau_ep_Tc_from_tau0,
)
from pysnspd.plotting.kinetic import (
    plot_eliashberg_spectrum,
    plot_gap_policy_delta_curves,
    plot_gap_policy_power_curves,
    plot_power_curve,
    plot_spectral_support,
)


MEV_J = 1.602176634e-22


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="OE5 projected electron-phonon power debug run."
    )

    parser.add_argument("--config", required=True, help="Path to YAML config.")
    parser.add_argument("--run-name", required=True, help="OE5 debug run name.")

    parser.add_argument(
        "--phase-space-npz",
        required=True,
        help="Path to OE4 phase_space_catalog.npz.",
    )
    parser.add_argument(
        "--eliashberg-dat",
        required=True,
        help="Path to Simon/MIT nbn-a2f-ph.dat.",
    )
    parser.add_argument(
        "--usadel-summary-yaml",
        default=None,
        help=(
            "Optional path to usadel_dos_summary.yaml. If omitted, the pipeline "
            "looks next to the phase-space NPZ."
        ),
    )

    parser.add_argument("--Tph-K", type=float, default=0.9)
    parser.add_argument("--Te-min-K", type=float, default=0.9)
    parser.add_argument("--Te-max-K", type=float, default=34.6)
    parser.add_argument("--n-Te", type=int, default=100)

    parser.add_argument(
        "--delta-index",
        type=int,
        default=-1,
        help="Delta index in the OE4 phase-space catalogue. Default: largest Delta.",
    )
    parser.add_argument(
        "--q-index",
        type=int,
        default=-1,
        help="q index in the OE4 phase-space catalogue. Default: middle q.",
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
            "converts it to Vodolazov tau0 using "
            "tau0 = [720 zeta(5)/pi^2] tau_ep(Tc)."
        ),
    )
    tau_group.add_argument(
        "--tau0-ps",
        type=float,
        default=None,
        help=(
            "Direct Vodolazov/Allmaras tau0 in ps. Use only if tau0 itself is "
            "known, not tau_ep(Tc)."
        ),
    )

    parser.add_argument(
        "--representative-Te-K",
        type=float,
        default=None,
        help=(
            "Representative Te for spectral-support plot. Default: min(0.8 Tc, "
            "Te_max), which avoids using Te >> Tc with a fixed full gap."
        ),
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
    eliashberg_path = Path(args.eliashberg_dat).expanduser().resolve()

    phase_catalog = load_phase_space_catalog_npz(phase_path)
    spectrum = load_simon_eliashberg_dat(eliashberg_path)

    D_m2_s = _resolve_D_m2_s(
        cfg,
        phase_path=phase_path,
        usadel_summary_yaml=args.usadel_summary_yaml,
    )
    sigma_n = float(cfg["material"]["sigma_n_S_m"])
    Tc_K = float(cfg["material"]["Tc_K"])
    N0 = electronic_density_of_states_from_sigma_D(sigma_n, D_m2_s)

    if args.tau0_ps is not None:
        tau0_s = float(args.tau0_ps) * 1.0e-12
        tau_ep_Tc_s = tau_ep_Tc_from_tau0(tau0_s)
        tau_source = "direct_tau0"
    else:
        tau_ep_Tc_s = float(args.tau_ep_Tc_ps) * 1.0e-12
        tau0_s = tau0_from_tau_ep_Tc(tau_ep_Tc_s)
        tau_source = "converted_from_tau_ep_Tc"

    q_index = args.q_index
    if q_index < 0:
        q_index = phase_catalog.q_values_m_inv.size // 2

    delta_index = args.delta_index
    if delta_index < 0:
        delta_index = phase_catalog.delta_values_J.size - 1

    delta0_J = float(phase_catalog.delta_values_J[delta_index])
    q_m_inv = float(phase_catalog.q_values_m_inv[q_index])

    Te_values = np.linspace(float(args.Te_min_K), float(args.Te_max_K), int(args.n_Te))

    fixed_curve = compute_power_curve(
        Te_values,
        Tph_K=float(args.Tph_K),
        delta_J=delta0_J,
        q_m_inv=q_m_inv,
        phase_space_catalog=phase_catalog,
        spectrum=spectrum,
        N0_J_m3=N0,
        tau0_s=tau0_s,
        Tc_K=Tc_K,
        omega_max_meV=args.omega_max_meV,
    )

    gap_factor = diagnostic_bcs_gap_factor(Te_values, Tc_K)
    bcs_delta_curve = delta0_J * gap_factor
    bcs_q_curve = q_m_inv * (gap_factor > 0.0).astype(float)

    bcs_curve = compute_power_curve(
        Te_values,
        Tph_K=float(args.Tph_K),
        delta_J=delta0_J,
        q_m_inv=q_m_inv,
        phase_space_catalog=phase_catalog,
        spectrum=spectrum,
        N0_J_m3=N0,
        tau0_s=tau0_s,
        Tc_K=Tc_K,
        omega_max_meV=args.omega_max_meV,
        delta_values_J=bcs_delta_curve,
        q_values_m_inv=bcs_q_curve,
    )

    normal_delta_curve = np.zeros_like(Te_values)
    normal_q_curve = np.zeros_like(Te_values)

    normal_curve = compute_power_curve(
        Te_values,
        Tph_K=float(args.Tph_K),
        delta_J=0.0,
        q_m_inv=0.0,
        phase_space_catalog=phase_catalog,
        spectrum=spectrum,
        N0_J_m3=N0,
        tau0_s=tau0_s,
        Tc_K=Tc_K,
        omega_max_meV=args.omega_max_meV,
        delta_values_J=normal_delta_curve,
        q_values_m_inv=normal_q_curve,
    )

    if args.representative_Te_K is None:
        representative_Te = min(0.8 * Tc_K, float(args.Te_max_K))
        representative_Te = max(representative_Te, float(args.Tph_K))
    else:
        representative_Te = float(args.representative_Te_K)

    representative_delta = float(
        delta0_J * diagnostic_bcs_gap_factor(np.array([representative_Te]), Tc_K)[0]
    )
    representative_q = q_m_inv if representative_delta > 0.0 else 0.0

    representative = compute_projected_powers(
        representative_Te,
        float(args.Tph_K),
        representative_delta,
        representative_q,
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
    power_plot = plot_power_curve(
        fixed_curve,
        plots_comp / "electron_phonon_power_vs_Te_fixed_gap.png",
        tau_label=r"$\tau_0$ from $\tau_{ep}(T_c)$",
        title_suffix="fixed gap diagnostic",
    )
    gap_policy_plot = plot_gap_policy_power_curves(
        {
            "fixed gap": fixed_curve,
            "BCS-like gap": bcs_curve,
            "normal": normal_curve,
        },
        plots_comp / "electron_phonon_power_gap_policy.png",
        Tc_K=Tc_K,
    )
    delta_policy_plot = plot_gap_policy_delta_curves(
        {
            "fixed gap": fixed_curve,
            "BCS-like gap": bcs_curve,
            "normal": normal_curve,
        },
        plots_diag / "gap_policy_delta_vs_Te.png",
        Tc_K=Tc_K,
    )
    support_plot = plot_spectral_support(
        support,
        plots_comp / "spectral_support_cumulative.png",
    )

    power_npz = raw_ss / "oe5_power_catalog.npz"
    np.savez_compressed(
        power_npz,
        Te_values_K=fixed_curve["Te_values_K"],
        fixed_delta_values_J=fixed_curve["delta_values_J"],
        fixed_q_values_m_inv=fixed_curve["q_values_m_inv"],
        fixed_P_S_W_m3=fixed_curve["P_S_W_m3"],
        fixed_P_R_W_m3=fixed_curve["P_R_W_m3"],
        fixed_P_total_W_m3=fixed_curve["P_total_W_m3"],
        fixed_P_Debye_Vodolazov_W_m3=fixed_curve["P_Debye_Vodolazov_W_m3"],
        bcs_delta_values_J=bcs_curve["delta_values_J"],
        bcs_q_values_m_inv=bcs_curve["q_values_m_inv"],
        bcs_P_S_W_m3=bcs_curve["P_S_W_m3"],
        bcs_P_R_W_m3=bcs_curve["P_R_W_m3"],
        bcs_P_total_W_m3=bcs_curve["P_total_W_m3"],
        normal_delta_values_J=normal_curve["delta_values_J"],
        normal_q_values_m_inv=normal_curve["q_values_m_inv"],
        normal_P_S_W_m3=normal_curve["P_S_W_m3"],
        normal_P_R_W_m3=normal_curve["P_R_W_m3"],
        normal_P_total_W_m3=normal_curve["P_total_W_m3"],
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
        "backend": "oe5_projected_power_debug_v2",
        "inputs": {
            "phase_space_npz": str(phase_path),
            "eliashberg_dat": str(eliashberg_path),
        },
        "state": {
            "Tph_K": float(args.Tph_K),
            "Te_min_K": float(Te_values[0]),
            "Te_max_K": float(Te_values[-1]),
            "n_Te": int(Te_values.size),
            "representative_Te_K": float(representative_Te),
            "representative_delta_meV": float(representative_delta / MEV_J),
            "representative_q_m_inv": float(representative_q),
            "delta_index": int(delta_index),
            "delta_fixed_meV": float(delta0_J / MEV_J),
            "q_index": int(q_index),
            "q_fixed_m_inv": q_m_inv,
            "omega_max_meV_requested": args.omega_max_meV,
        },
        "material": {
            "sigma_n_S_m": sigma_n,
            "D_m2_s": float(D_m2_s),
            "N0_J_m3": float(N0),
            "Tc_K": Tc_K,
            "tau_source": tau_source,
            "tau_ep_Tc_ps": float(tau_ep_Tc_s / 1.0e-12),
            "tau0_ps": float(tau0_s / 1.0e-12),
            "tau0_over_tau_ep_Tc": float(TAU0_OVER_TAU_EP_TC),
            "tau_comment": (
                "Vodolazov tau0 is not tau_ep(Tc). The default comparison uses "
                "tau0 = [720 zeta(5)/pi^2] tau_ep(Tc)."
            ),
        },
        "eliashberg": spectrum_summary(spectrum),
        "representative_power": {
            "P_S_W_m3": representative.P_S_W_m3,
            "P_R_W_m3": representative.P_R_W_m3,
            "P_total_W_m3": representative.P_total_W_m3,
            "omega_max_meV_used": representative.metadata["omega_max_meV_used"],
        },
        "curve_extrema": {
            "fixed_P_S_min_W_m3": float(np.min(fixed_curve["P_S_W_m3"])),
            "fixed_P_S_max_W_m3": float(np.max(fixed_curve["P_S_W_m3"])),
            "fixed_P_R_min_W_m3": float(np.min(fixed_curve["P_R_W_m3"])),
            "fixed_P_R_max_W_m3": float(np.max(fixed_curve["P_R_W_m3"])),
            "fixed_P_total_min_W_m3": float(np.min(fixed_curve["P_total_W_m3"])),
            "fixed_P_total_max_W_m3": float(np.max(fixed_curve["P_total_W_m3"])),
            "bcs_P_R_max_W_m3": float(np.max(bcs_curve["P_R_W_m3"])),
            "bcs_P_total_max_W_m3": float(np.max(bcs_curve["P_total_W_m3"])),
            "normal_P_R_max_W_m3": float(np.max(normal_curve["P_R_W_m3"])),
            "normal_P_total_max_W_m3": float(np.max(normal_curve["P_total_W_m3"])),
            "P_Debye_min_W_m3": float(np.min(fixed_curve["P_Debye_Vodolazov_W_m3"])),
            "P_Debye_max_W_m3": float(np.max(fixed_curve["P_Debye_Vodolazov_W_m3"])),
        },
        "diagnostic_interpretation": {
            "fixed_gap_warning": (
                "The fixed-gap curve is useful for debugging the OE4/OE5 tables, "
                "but it is not a self-consistent thermal trajectory for Te>Tc."
            ),
            "recombination_policy": (
                "P_R is retained as a Simon superconducting recombination/pair-breaking "
                "channel. It must be interpreted together with the superconducting "
                "energy u_e^S(Te,Delta,q); it should not be compared directly to the "
                "normal Debye T^5 scattering-only limit."
            ),
        },
        "outputs": {
            "power_npz": str(power_npz),
            "power_summary": str(raw_ss / "oe5_power_summary.yaml"),
            "eliashberg_plot": str(spectrum_plot),
            "fixed_gap_power_plot": str(power_plot),
            "gap_policy_power_plot": str(gap_policy_plot),
            "gap_policy_delta_plot": str(delta_policy_plot),
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
                "OE5 local projected electron-phonon power debug from OE4 "
                "phase-space and Simon/MIT Eliashberg material data."
            ),
            "summary": summary,
        },
    )

    print("OE5 projected-power debug")
    print(f"run_name                  : {args.run_name}")
    print(f"phase_space_npz           : {phase_path}")
    print(f"eliashberg_dat            : {eliashberg_path}")
    print(f"raw_ss                    : {raw_ss}")
    print(f"plots_diagnostics         : {plots_diag}")
    print(f"plots_comparisons         : {plots_comp}")
    print()
    print("Material")
    print(f"  D_m2_s                  : {D_m2_s}")
    print(f"  sigma_n_S_m             : {sigma_n}")
    print(f"  N0_J_m3                 : {N0}")
    print(f"  tau_source              : {tau_source}")
    print(f"  tau_ep_Tc_ps            : {tau_ep_Tc_s / 1.0e-12}")
    print(f"  tau0_ps                 : {tau0_s / 1.0e-12}")
    print(f"  tau0/tau_ep_Tc          : {TAU0_OVER_TAU_EP_TC}")
    print()
    print("Eliashberg")
    for key in [
        "header",
        "frequency_max_THz",
        "omega_max_meV",
        "lambda_ep",
        "alpha2F_peak_meV",
        "phdos_peak_meV",
        "n_phdos_negative_clipped",
    ]:
        print(f"  {key}: {spectrum.metadata.get(key)}")
    print()
    print("Representative projected power")
    print(f"  Te_K                    : {representative_Te}")
    print(f"  Tph_K                   : {args.Tph_K}")
    print(f"  delta_meV               : {representative_delta / MEV_J}")
    print(f"  q_m_inv                 : {representative_q}")
    print(f"  P_S_W_m3                : {representative.P_S_W_m3}")
    print(f"  P_R_W_m3                : {representative.P_R_W_m3}")
    print(f"  P_total_W_m3            : {representative.P_total_W_m3}")
    print()
    print("Outputs")
    print(f"  power_npz               : {power_npz}")
    print(f"  power_summary           : {summary_path}")
    print(f"  eliashberg_plot         : {spectrum_plot}")
    print(f"  fixed_gap_power_plot    : {power_plot}")
    print(f"  gap_policy_power_plot   : {gap_policy_plot}")
    print(f"  gap_policy_delta_plot   : {delta_policy_plot}")
    print(f"  support_plot            : {support_plot}")
    print(f"  manifest                : {manifest}")
    print("Status: OK")

    return 0


def _resolve_D_m2_s(
    cfg: dict[str, Any],
    *,
    phase_path: Path,
    usadel_summary_yaml: str | None,
) -> float:
    candidates: list[Path] = []

    if usadel_summary_yaml is not None:
        candidates.append(Path(usadel_summary_yaml).expanduser().resolve())

    candidates.append(phase_path.parent / "usadel_dos_summary.yaml")

    for candidate in candidates:
        if candidate.exists():
            with candidate.open("r", encoding="utf-8") as f:
                data = yaml.safe_load(f)
            if isinstance(data, dict):
                if "D_m2_s" in data:
                    return float(data["D_m2_s"])

                usadel = data.get("usadel", {})
                if isinstance(usadel, dict) and "D_m2_s" in usadel:
                    return float(usadel["D_m2_s"])

    material = cfg.get("material", {})
    if "D_m2_s" in material:
        return float(material["D_m2_s"])

    raise RuntimeError(
        "Could not resolve D_m2_s. Provide --usadel-summary-yaml or keep "
        "usadel_dos_summary.yaml next to phase_space_catalog.npz."
    )


if __name__ == "__main__":
    raise SystemExit(main())