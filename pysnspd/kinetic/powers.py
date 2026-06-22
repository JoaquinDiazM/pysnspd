"""Projected microscopic powers for the thermal model.

OE5 converts the expensive OE4 phase-space catalogues,

    J_S(Omega; Te, Delta, q)
    J_R(Omega; Te, Delta, q),

into projected electron-phonon powers by multiplying by the material
Eliashberg spectral function and Bose imbalance.

The implemented formulas follow Appendix A of the pySNSPD thesis draft:

    P_ep^S = (8 pi N(0)/hbar) int dOmega alpha^2F(Omega) Omega
             [n_e(Omega,Te)-n_ph(Omega,Tph)] J_S(Omega),

    P_ep^R = (4 pi N(0)/hbar) int dOmega alpha^2F(Omega) Omega
             [n_e(Omega,Te)-n_ph(Omega,Tph)] J_R(Omega).

The microscopic starting point is the kinetic formulation of Simon et al.,
Physical Review B 112, 174512 (2025). The Vodolazov/Allmaras T^5 form is
used only as a normal-state Debye consistency check, not as the microscopic
starting point.

Important convention
--------------------
The Vodolazov/Allmaras parameter tau0 is not the linear electron-phonon
relaxation time at Tc. Linearizing the Debye T^5 power around T=Tc gives

    tau0 = [720 zeta(5) / pi^2] tau_ep(Tc).

Therefore a material value such as tau_ep(Tc)=24.7 ps must not be inserted
directly as tau0 in the T^5 comparison.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np

from pysnspd.kinetic.eliashberg import EliashbergSpectrum, j_to_mev


HBAR_J_S = 1.054571817e-34
KB_J_K = 1.380649e-23
E_CHARGE_C = 1.602176634e-19
ZETA_5 = 1.03692775514337
MEV_J = 1.602176634e-22

TAU0_OVER_TAU_EP_TC = 720.0 * ZETA_5 / (np.pi**2)


@dataclass(frozen=True)
class ProjectedPowerResult:
    """Projected local electron-phonon powers for one thermodynamic state."""

    Te_K: float
    Tph_K: float
    delta_J: float
    q_m_inv: float
    N0_J_m3: float
    P_S_W_m3: float
    P_R_W_m3: float
    P_total_W_m3: float
    omega_J: np.ndarray
    alpha2F: np.ndarray
    bose_difference: np.ndarray
    J_S_J: np.ndarray
    J_R_J: np.ndarray
    integrand_S_J2: np.ndarray
    integrand_R_J2: np.ndarray
    metadata: dict[str, Any]


@dataclass(frozen=True)
class UsadelSelfConsistentTrajectory:
    """Stable-branch Usadel trajectory Delta_eq(q,T_bias).

    This object is extracted from the OE3 Usadel calibration sweep. It is
    self-consistent in the sense used by the current pySNSPD OE3 block:
    Delta_eq is obtained from the Matsubara Usadel gap equation for each
    depairing energy Gamma_q at the configured bias temperature T_bias.

    It is not a substitute for the future time-dependent gTDGL field
    Delta(r,t). It is the correct OE5 diagnostic trajectory for avoiding the
    artificial fixed-gap plots used during earlier debugging.
    """

    q_values_m_inv: np.ndarray
    gamma_values_J: np.ndarray
    delta_eq_values_J: np.ndarray
    current_values_A: np.ndarray
    current_density_values_A_m2: np.ndarray
    current_fraction: np.ndarray
    metadata: dict[str, Any]


def tau0_from_tau_ep_Tc(tau_ep_Tc_s: float) -> float:
    """Convert the linear relaxation time at Tc into Vodolazov tau0."""
    if tau_ep_Tc_s <= 0.0:
        raise ValueError("tau_ep_Tc_s must be positive.")
    return float(TAU0_OVER_TAU_EP_TC * tau_ep_Tc_s)


def tau_ep_Tc_from_tau0(tau0_s: float) -> float:
    """Inverse of :func:`tau0_from_tau_ep_Tc`."""
    if tau0_s <= 0.0:
        raise ValueError("tau0_s must be positive.")
    return float(tau0_s / TAU0_OVER_TAU_EP_TC)


def electronic_density_of_states_from_sigma_D(
    sigma_n_S_m: float,
    D_m2_s: float,
) -> float:
    """Return single-spin N(0) from the dirty-limit Einstein relation.

    Simon et al. use

        rho_N = 1 / [2 e^2 D N(0)],

    so, with sigma_n = 1/rho_N,

        N(0) = sigma_n / [2 e^2 D].

    Units: J^-1 m^-3.
    """
    if sigma_n_S_m <= 0.0:
        raise ValueError("sigma_n_S_m must be positive.")
    if D_m2_s <= 0.0:
        raise ValueError("D_m2_s must be positive.")
    return float(sigma_n_S_m / (2.0 * E_CHARGE_C**2 * D_m2_s))


def bose_positive_energy(omega_J: np.ndarray, T_K: float) -> np.ndarray:
    """Bose-Einstein occupation for positive energy Omega."""
    if T_K <= 0.0:
        raise ValueError("Temperature must be positive.")

    omega = np.asarray(omega_J, dtype=float)
    out = np.zeros_like(omega, dtype=float)

    positive = omega > 0.0
    x = omega[positive] / (KB_J_K * float(T_K))

    vals = np.zeros_like(x)
    small = x < 1.0e-6
    vals[small] = 1.0 / x[small] - 0.5 + x[small] / 12.0
    vals[~small] = 1.0 / np.expm1(np.clip(x[~small], 0.0, 700.0))

    out[positive] = vals
    return out


def bose_difference(omega_J: np.ndarray, Te_K: float, Tph_K: float) -> np.ndarray:
    """Return n_e(Omega,Te)-n_ph(Omega,Tph)."""
    return bose_positive_energy(omega_J, Te_K) - bose_positive_energy(omega_J, Tph_K)


def build_usadel_self_consistent_trajectory(
    usadel_catalog,
    *,
    n_q: int = 120,
    stable_branch_only: bool = True,
) -> UsadelSelfConsistentTrajectory:
    """Build an interpolated stable-branch Usadel trajectory.

    The OE3 Usadel block stores a calibration sweep containing
    ``calibration_q_values_m_inv``, ``calibration_delta_eq_values_J`` and
    ``calibration_current_values_A``. This function extracts the stable branch
    up to the maximum supercurrent and interpolates it onto a smooth q-axis.
    """
    q_raw = np.asarray(usadel_catalog.calibration_q_values_m_inv, dtype=float)
    gamma_raw = np.asarray(usadel_catalog.calibration_gamma_values_J, dtype=float)
    delta_raw = np.asarray(usadel_catalog.calibration_delta_eq_values_J, dtype=float)
    current_raw = np.asarray(usadel_catalog.calibration_current_values_A, dtype=float)
    current_density_raw = np.asarray(
        usadel_catalog.calibration_current_density_values_A_m2,
        dtype=float,
    )

    finite = (
        np.isfinite(q_raw)
        & np.isfinite(gamma_raw)
        & np.isfinite(delta_raw)
        & np.isfinite(current_raw)
        & np.isfinite(current_density_raw)
    )
    finite &= q_raw >= 0.0
    finite &= delta_raw >= 0.0
    finite &= current_raw >= 0.0

    if np.sum(finite) < 3:
        raise ValueError("Usadel calibration sweep does not contain enough valid points.")

    q_raw = q_raw[finite]
    gamma_raw = gamma_raw[finite]
    delta_raw = delta_raw[finite]
    current_raw = current_raw[finite]
    current_density_raw = current_density_raw[finite]

    order = np.argsort(q_raw)
    q_raw = q_raw[order]
    gamma_raw = gamma_raw[order]
    delta_raw = delta_raw[order]
    current_raw = current_raw[order]
    current_density_raw = current_density_raw[order]

    idx_ic = int(np.argmax(current_raw))
    Ic_A = float(current_raw[idx_ic])
    q_critical = float(q_raw[idx_ic])
    delta_critical = float(delta_raw[idx_ic])

    if stable_branch_only:
        q_raw = q_raw[: idx_ic + 1]
        gamma_raw = gamma_raw[: idx_ic + 1]
        delta_raw = delta_raw[: idx_ic + 1]
        current_raw = current_raw[: idx_ic + 1]
        current_density_raw = current_density_raw[: idx_ic + 1]

    if Ic_A <= 0.0:
        raise ValueError("Usadel calibration critical current is not positive.")

    if n_q <= 1:
        q_values = q_raw
    else:
        q_values = np.linspace(float(q_raw[0]), float(q_raw[-1]), int(n_q))

    gamma_values = np.interp(q_values, q_raw, gamma_raw)
    delta_values = np.interp(q_values, q_raw, delta_raw)
    current_values = np.interp(q_values, q_raw, current_raw)
    current_density_values = np.interp(q_values, q_raw, current_density_raw)
    current_fraction = current_values / Ic_A

    metadata = {
        "backend": "usadel_self_consistent_trajectory_from_oe3_calibration",
        "stable_branch_only": bool(stable_branch_only),
        "n_q": int(q_values.size),
        "Ic_A": Ic_A,
        "q_critical_m_inv": q_critical,
        "delta_critical_J": delta_critical,
        "delta_critical_meV": float(delta_critical / MEV_J),
        "T_bias_K": float(usadel_catalog.metadata.get("T_bias_K", np.nan)),
        "source": (
            "Extracted from OE3 Matsubara Usadel calibration sweep. "
            "Delta_eq(q,T_bias) is used as the self-consistent OE5 diagnostic "
            "trajectory."
        ),
    }

    return UsadelSelfConsistentTrajectory(
        q_values_m_inv=q_values,
        gamma_values_J=gamma_values,
        delta_eq_values_J=delta_values,
        current_values_A=current_values,
        current_density_values_A_m2=current_density_values,
        current_fraction=current_fraction,
        metadata=metadata,
    )


def select_usadel_state_by_current_fraction(
    trajectory: UsadelSelfConsistentTrajectory,
    target_fraction: float,
) -> dict[str, float]:
    """Select the self-consistent trajectory point closest to I/Ic."""
    if target_fraction < 0.0:
        raise ValueError("target_fraction must be non-negative.")

    frac = np.asarray(trajectory.current_fraction, dtype=float)
    idx = int(np.argmin(np.abs(frac - float(target_fraction))))

    return {
        "index": idx,
        "current_fraction": float(trajectory.current_fraction[idx]),
        "current_A": float(trajectory.current_values_A[idx]),
        "current_density_A_m2": float(trajectory.current_density_values_A_m2[idx]),
        "q_m_inv": float(trajectory.q_values_m_inv[idx]),
        "gamma_J": float(trajectory.gamma_values_J[idx]),
        "gamma_meV": float(trajectory.gamma_values_J[idx] / MEV_J),
        "delta_J": float(trajectory.delta_eq_values_J[idx]),
        "delta_meV": float(trajectory.delta_eq_values_J[idx] / MEV_J),
    }


def phase_space_spectra_at_state(
    phase_space_catalog,
    *,
    Te_K: float,
    delta_J: float,
    q_m_inv: float,
) -> tuple[np.ndarray, np.ndarray]:
    """Trilinearly interpolate OE4 J_S and J_R at a local state."""
    t_axis = np.asarray(phase_space_catalog.Te_values_K, dtype=float)
    d_axis = np.asarray(phase_space_catalog.delta_values_J, dtype=float)
    q_axis = np.asarray(phase_space_catalog.q_values_m_inv, dtype=float)

    JS = np.asarray(phase_space_catalog.J_S_TdqO_J, dtype=float)
    JR = np.asarray(phase_space_catalog.J_R_TdqO_J, dtype=float)

    t0, t1, wt = _bracket_weight(t_axis, Te_K)
    d0, d1, wd = _bracket_weight(d_axis, delta_J)
    q0, q1, wq = _bracket_weight(q_axis, q_m_inv)

    out_S = _trilinear(JS, t0, t1, wt, d0, d1, wd, q0, q1, wq)
    out_R = _trilinear(JR, t0, t1, wt, d0, d1, wd, q0, q1, wq)
    return out_S, out_R


def compute_projected_powers(
    Te_K: float,
    Tph_K: float,
    delta_J: float,
    q_m_inv: float,
    phase_space_catalog,
    spectrum: EliashbergSpectrum,
    *,
    N0_J_m3: float,
    omega_max_meV: float | None = None,
) -> ProjectedPowerResult:
    """Compute OE5 projected powers for one local state.

    Positive power means energy leaves the electronic system and enters the
    phonon system.
    """
    if N0_J_m3 <= 0.0:
        raise ValueError("N0_J_m3 must be positive.")

    omega = np.asarray(phase_space_catalog.omega_values_J, dtype=float)
    if omega_max_meV is not None:
        omega_max_J = float(omega_max_meV) * MEV_J
        mask = omega <= omega_max_J
    else:
        mask = np.ones_like(omega, dtype=bool)

    if np.sum(mask) < 2:
        raise ValueError("Need at least two Omega points in the selected range.")

    omega_sel = omega[mask]
    alpha = spectrum.alpha2F_on_omega_J(omega_sel)
    d_bose = bose_difference(omega_sel, Te_K, Tph_K)

    JS_all, JR_all = phase_space_spectra_at_state(
        phase_space_catalog,
        Te_K=Te_K,
        delta_J=delta_J,
        q_m_inv=q_m_inv,
    )
    JS = np.asarray(JS_all[mask], dtype=float)
    JR = np.asarray(JR_all[mask], dtype=float)

    integrand_S = alpha * omega_sel * d_bose * JS
    integrand_R = alpha * omega_sel * d_bose * JR

    P_S = (8.0 * np.pi * N0_J_m3 / HBAR_J_S) * float(
        np.trapezoid(integrand_S, omega_sel)
    )
    P_R = (4.0 * np.pi * N0_J_m3 / HBAR_J_S) * float(
        np.trapezoid(integrand_R, omega_sel)
    )

    metadata = {
        "backend": "projected_powers_oe5_v3_usadel_self_consistent",
        "sign_convention": (
            "Positive P_S/P_R means energy leaves electrons and enters phonons."
        ),
        "source": "pySNSPD Appendix A based on Simon et al. 2025 kinetic equations.",
        "omega_max_meV_used": float(j_to_mev(omega_sel[-1])),
        "alpha2F_source": spectrum.metadata.get("source", ""),
        "alpha2F_path": spectrum.metadata.get("path", ""),
        "alpha2F_policy": spectrum.metadata.get("alpha2F_policy", ""),
    }

    return ProjectedPowerResult(
        Te_K=float(Te_K),
        Tph_K=float(Tph_K),
        delta_J=float(delta_J),
        q_m_inv=float(q_m_inv),
        N0_J_m3=float(N0_J_m3),
        P_S_W_m3=float(P_S),
        P_R_W_m3=float(P_R),
        P_total_W_m3=float(P_S + P_R),
        omega_J=omega_sel,
        alpha2F=alpha,
        bose_difference=d_bose,
        J_S_J=JS,
        J_R_J=JR,
        integrand_S_J2=integrand_S,
        integrand_R_J2=integrand_R,
        metadata=metadata,
    )


def compute_scattering_power(
    Te: float,
    Tph: float,
    delta: float,
    q: float,
    phase_space_catalog,
    alpha2F: EliashbergSpectrum,
    *,
    N0_J_m3: float,
    omega_max_meV: float | None = None,
) -> float:
    """Compute ``P_ep^S`` from the projected Simon scattering channel."""
    return compute_projected_powers(
        Te,
        Tph,
        delta,
        q,
        phase_space_catalog,
        alpha2F,
        N0_J_m3=N0_J_m3,
        omega_max_meV=omega_max_meV,
    ).P_S_W_m3


def compute_recombination_power(
    Te: float,
    Tph: float,
    delta: float,
    q: float,
    phase_space_catalog,
    alpha2F: EliashbergSpectrum,
    *,
    N0_J_m3: float,
    omega_max_meV: float | None = None,
) -> float:
    """Compute ``P_ep^R`` from the projected recombination/pair-breaking channel."""
    return compute_projected_powers(
        Te,
        Tph,
        delta,
        q,
        phase_space_catalog,
        alpha2F,
        N0_J_m3=N0_J_m3,
        omega_max_meV=omega_max_meV,
    ).P_R_W_m3


def compute_power_curve_at_usadel_state(
    Te_values_K: np.ndarray,
    *,
    Tph_K: float,
    state: dict[str, float],
    phase_space_catalog,
    spectrum: EliashbergSpectrum,
    N0_J_m3: float,
    tau0_s: float,
    Tc_K: float,
    omega_max_meV: float | None = None,
) -> dict[str, np.ndarray]:
    """Compute projected powers versus Te at one self-consistent Usadel state."""
    Te_values = np.asarray(Te_values_K, dtype=float)

    P_S = np.zeros_like(Te_values)
    P_R = np.zeros_like(Te_values)
    P_total = np.zeros_like(Te_values)
    P_D = np.zeros_like(Te_values)

    for i, Te in enumerate(Te_values):
        result = compute_projected_powers(
            float(Te),
            float(Tph_K),
            float(state["delta_J"]),
            float(state["q_m_inv"]),
            phase_space_catalog,
            spectrum,
            N0_J_m3=float(N0_J_m3),
            omega_max_meV=omega_max_meV,
        )
        P_S[i] = result.P_S_W_m3
        P_R[i] = result.P_R_W_m3
        P_total[i] = result.P_total_W_m3
        P_D[i] = compute_vodolazov_debye_power_density(
            float(Te),
            float(Tph_K),
            N0_J_m3=float(N0_J_m3),
            tau0_s=float(tau0_s),
            Tc_K=float(Tc_K),
        )

    return {
        "Te_values_K": Te_values,
        "delta_values_J": np.full_like(Te_values, float(state["delta_J"])),
        "q_values_m_inv": np.full_like(Te_values, float(state["q_m_inv"])),
        "P_S_W_m3": P_S,
        "P_R_W_m3": P_R,
        "P_total_W_m3": P_total,
        "P_Debye_Vodolazov_W_m3": P_D,
    }


def compute_usadel_q_power_scan(
    Te_values_K: np.ndarray,
    *,
    Tph_K: float,
    trajectory: UsadelSelfConsistentTrajectory,
    phase_space_catalog,
    spectrum: EliashbergSpectrum,
    N0_J_m3: float,
    omega_max_meV: float | None = None,
) -> dict[str, np.ndarray]:
    """Compute projected powers along the self-consistent Usadel q trajectory."""
    Te_values = np.asarray(Te_values_K, dtype=float)
    q_values = np.asarray(trajectory.q_values_m_inv, dtype=float)
    delta_values = np.asarray(trajectory.delta_eq_values_J, dtype=float)

    shape = (Te_values.size, q_values.size)
    P_S = np.zeros(shape, dtype=float)
    P_R = np.zeros(shape, dtype=float)
    P_total = np.zeros(shape, dtype=float)

    for iT, Te in enumerate(Te_values):
        for iq, (q, delta) in enumerate(zip(q_values, delta_values, strict=True)):
            result = compute_projected_powers(
                float(Te),
                float(Tph_K),
                float(delta),
                float(q),
                phase_space_catalog,
                spectrum,
                N0_J_m3=float(N0_J_m3),
                omega_max_meV=omega_max_meV,
            )
            P_S[iT, iq] = result.P_S_W_m3
            P_R[iT, iq] = result.P_R_W_m3
            P_total[iT, iq] = result.P_total_W_m3

    return {
        "Te_values_K": Te_values,
        "q_values_m_inv": q_values,
        "gamma_values_J": np.asarray(trajectory.gamma_values_J, dtype=float),
        "delta_values_J": delta_values,
        "current_values_A": np.asarray(trajectory.current_values_A, dtype=float),
        "current_fraction": np.asarray(trajectory.current_fraction, dtype=float),
        "P_S_W_m3": P_S,
        "P_R_W_m3": P_R,
        "P_total_W_m3": P_total,
    }


def compute_vodolazov_debye_power_density(
    Te_K: float,
    Tph_K: float,
    *,
    N0_J_m3: float,
    tau0_s: float,
    Tc_K: float,
) -> float:
    """Normal-state Debye/Vodolazov electron-phonon power density."""
    if N0_J_m3 <= 0.0:
        raise ValueError("N0_J_m3 must be positive.")
    if tau0_s <= 0.0:
        raise ValueError("tau0_s must be positive.")
    if Tc_K <= 0.0:
        raise ValueError("Tc_K must be positive.")

    prefactor = 96.0 * ZETA_5 * N0_J_m3 * KB_J_K**2 / (tau0_s * Tc_K**3)
    return float(prefactor * (Te_K**5 - Tph_K**5))


def compute_total_electron_phonon_power(
    Te: float,
    Tph: float,
    delta: float,
    q: float,
    catalogs,
) -> float:
    """Compatibility wrapper for future thermal solver."""
    result = compute_projected_powers(
        Te,
        Tph,
        delta,
        q,
        catalogs["phase_space_catalog"],
        catalogs["eliashberg_spectrum"],
        N0_J_m3=float(catalogs["N0_J_m3"]),
    )
    return result.P_total_W_m3


def compute_escape_power(Tph, Tbath, phonon_dos, tau_esc):
    """Placeholder for OE6 phonon escape power."""
    raise NotImplementedError("Phonon escape belongs to the next thermal-balance OE.")


def cumulative_spectral_support(result: ProjectedPowerResult) -> dict[str, np.ndarray]:
    """Return cumulative support curves for the OE5 spectral integrands."""
    omega = result.omega_J
    alpha_weight = np.abs(result.alpha2F * omega)
    S_weight = np.abs(result.integrand_S_J2)
    R_weight = np.abs(result.integrand_R_J2)

    return {
        "omega_J": omega,
        "omega_meV": np.asarray(j_to_mev(omega), dtype=float),
        "cumulative_alpha_omega": _cumulative_fraction(omega, alpha_weight),
        "cumulative_scattering": _cumulative_fraction(omega, S_weight),
        "cumulative_recombination": _cumulative_fraction(omega, R_weight),
    }


def _bracket_weight(axis: np.ndarray, value: float) -> tuple[int, int, float]:
    axis = np.asarray(axis, dtype=float)
    if axis.ndim != 1 or axis.size == 0:
        raise ValueError("Interpolation axis must be one-dimensional and non-empty.")

    if value <= axis[0]:
        return 0, 0, 0.0
    if value >= axis[-1]:
        last = axis.size - 1
        return last, last, 0.0

    hi = int(np.searchsorted(axis, value, side="right"))
    lo = hi - 1
    denom = axis[hi] - axis[lo]
    if denom <= 0.0:
        return lo, hi, 0.0
    w = float((value - axis[lo]) / denom)
    return lo, hi, w


def _trilinear(
    arr: np.ndarray,
    t0: int,
    t1: int,
    wt: float,
    d0: int,
    d1: int,
    wd: float,
    q0: int,
    q1: int,
    wq: float,
) -> np.ndarray:
    c000 = arr[t0, d0, q0, :]
    c001 = arr[t0, d0, q1, :]
    c010 = arr[t0, d1, q0, :]
    c011 = arr[t0, d1, q1, :]
    c100 = arr[t1, d0, q0, :]
    c101 = arr[t1, d0, q1, :]
    c110 = arr[t1, d1, q0, :]
    c111 = arr[t1, d1, q1, :]

    c00 = c000 * (1.0 - wq) + c001 * wq
    c01 = c010 * (1.0 - wq) + c011 * wq
    c10 = c100 * (1.0 - wq) + c101 * wq
    c11 = c110 * (1.0 - wq) + c111 * wq

    c0 = c00 * (1.0 - wd) + c01 * wd
    c1 = c10 * (1.0 - wd) + c11 * wd

    return c0 * (1.0 - wt) + c1 * wt


def _cumulative_fraction(x: np.ndarray, weight: np.ndarray) -> np.ndarray:
    x = np.asarray(x, dtype=float)
    w = np.maximum(np.asarray(weight, dtype=float), 0.0)
    increments = 0.5 * (w[1:] + w[:-1]) * np.diff(x)
    cumulative = np.concatenate([[0.0], np.cumsum(increments)])
    total = float(cumulative[-1])
    if total <= 0.0:
        return np.zeros_like(cumulative)
    return cumulative / total