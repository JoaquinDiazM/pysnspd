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


def tau0_from_tau_ep_Tc(tau_ep_Tc_s: float) -> float:
    """Convert the linear relaxation time at Tc into Vodolazov tau0.

    The Debye/Vodolazov power is

        P_D = 96 zeta(5) N(0) k_B^2/(tau0 Tc^3) * (Te^5 - Tph^5).

    Linearizing this expression around T=Tc and using

        C_e^N(T) = (2 pi^2/3) N(0) k_B^2 T

    gives

        tau_ep(Tc) = [pi^2/(720 zeta(5))] tau0.

    Hence tau0 is approximately 75.6 times larger than tau_ep(Tc).
    """
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
    """Bose-Einstein occupation for positive energy Omega.

    At Omega=0 the mathematical occupation diverges. For OE5 power integrals
    the point Omega=0 has zero measure and alpha^2F is zero or negligible.
    We set the exact zero point to zero to avoid inf-inf cancellation in
    [n_e - n_ph].
    """
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


def diagnostic_bcs_gap_factor(Te_K: np.ndarray, Tc_K: float) -> np.ndarray:
    """Return a BCS-like equilibrium gap factor Delta(T)/Delta(0).

    This is a diagnostic closure only. It is not the final coupled gTDGL
    prescription. It is used in OE5 to demonstrate that the large fixed-gap
    recombination power at Te >> Tc is not a physically self-consistent state.
    """
    Te = np.asarray(Te_K, dtype=float)
    if Tc_K <= 0.0:
        raise ValueError("Tc_K must be positive.")

    out = np.zeros_like(Te, dtype=float)
    mask = (Te > 0.0) & (Te < Tc_K)
    x = np.sqrt(np.maximum(Tc_K / Te[mask] - 1.0, 0.0))
    out[mask] = np.tanh(1.74 * x)
    return np.clip(out, 0.0, 1.0)


def phase_space_spectra_at_state(
    phase_space_catalog,
    *,
    Te_K: float,
    delta_J: float,
    q_m_inv: float,
) -> tuple[np.ndarray, np.ndarray]:
    """Trilinearly interpolate OE4 J_S and J_R at a local state.

    The interpolation is over the catalogue axes Te, Delta and q. The returned
    spectra live on the catalogue Omega axis.
    """
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
    result = compute_projected_powers(
        Te,
        Tph,
        delta,
        q,
        phase_space_catalog,
        alpha2F,
        N0_J_m3=N0_J_m3,
        omega_max_meV=omega_max_meV,
    )
    return result.P_S_W_m3


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
    result = compute_projected_powers(
        Te,
        Tph,
        delta,
        q,
        phase_space_catalog,
        alpha2F,
        N0_J_m3=N0_J_m3,
        omega_max_meV=omega_max_meV,
    )
    return result.P_R_W_m3


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
    phonon system. In the electron-temperature equation, these terms enter
    with a minus sign.
    """
    if N0_J_m3 <= 0.0:
        raise ValueError("N0_J_m3 must be positive.")

    omega = np.asarray(phase_space_catalog.omega_values_J, dtype=float)
    if omega_max_meV is not None:
        omega_max_J = float(omega_max_meV) * 1.602176634e-22
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
        "backend": "projected_powers_oe5_v2",
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


def compute_total_electron_phonon_power(
    Te: float,
    Tph: float,
    delta: float,
    q: float,
    catalogs,
) -> float:
    """Compute the net projected electron-phonon power used by thermal solver.

    This compatibility wrapper expects ``catalogs`` to contain:

        phase_space_catalog
        eliashberg_spectrum
        N0_J_m3
    """
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
    """Placeholder for OE6 phonon escape power.

    Escape is not part of OE5 projected electron-phonon collision powers. It
    will be implemented when the phonon energy balance is activated.
    """
    raise NotImplementedError("Phonon escape belongs to the next thermal-balance OE.")


def compute_power_curve(
    Te_values_K: np.ndarray,
    *,
    Tph_K: float,
    delta_J: float,
    q_m_inv: float,
    phase_space_catalog,
    spectrum: EliashbergSpectrum,
    N0_J_m3: float,
    tau0_s: float,
    Tc_K: float,
    omega_max_meV: float | None = None,
    delta_values_J: np.ndarray | None = None,
    q_values_m_inv: np.ndarray | None = None,
) -> dict[str, np.ndarray]:
    """Compute projected powers as a function of Te for diagnostics.

    ``delta_values_J`` and ``q_values_m_inv`` may be supplied to evaluate a
    non-fixed trajectory through the catalogue, for example a BCS-like
    diagnostic Delta(Te) curve. This is only a 0D diagnostic; the final coupled
    model will use Delta(r,t) and q(r,t) from gTDGL.
    """
    Te_values = np.asarray(Te_values_K, dtype=float)
    n = Te_values.size

    if delta_values_J is None:
        delta_curve = np.full(n, float(delta_J), dtype=float)
    else:
        delta_curve = np.asarray(delta_values_J, dtype=float)
        if delta_curve.shape != Te_values.shape:
            raise ValueError("delta_values_J must have the same shape as Te_values_K.")

    if q_values_m_inv is None:
        q_curve = np.full(n, float(q_m_inv), dtype=float)
    else:
        q_curve = np.asarray(q_values_m_inv, dtype=float)
        if q_curve.shape != Te_values.shape:
            raise ValueError("q_values_m_inv must have the same shape as Te_values_K.")

    P_S = np.zeros_like(Te_values)
    P_R = np.zeros_like(Te_values)
    P_total = np.zeros_like(Te_values)
    P_D = np.zeros_like(Te_values)

    for i, Te in enumerate(Te_values):
        result = compute_projected_powers(
            float(Te),
            float(Tph_K),
            float(delta_curve[i]),
            float(q_curve[i]),
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
        "delta_values_J": delta_curve,
        "q_values_m_inv": q_curve,
        "P_S_W_m3": P_S,
        "P_R_W_m3": P_R,
        "P_total_W_m3": P_total,
        "P_Debye_Vodolazov_W_m3": P_D,
    }


def compute_vodolazov_debye_power_density(
    Te_K: float,
    Tph_K: float,
    *,
    N0_J_m3: float,
    tau0_s: float,
    Tc_K: float,
) -> float:
    """Normal-state Debye/Vodolazov electron-phonon power density.

    P_D = 96 zeta(5) N(0) k_B^2 / [tau0 Tc^3] * (Te^5 - Tph^5).

    This is used as a magnitude and limiting-form check for the scattering
    channel. It is not used to overwrite the Simon/Eliashberg projected power.
    """
    if N0_J_m3 <= 0.0:
        raise ValueError("N0_J_m3 must be positive.")
    if tau0_s <= 0.0:
        raise ValueError("tau0_s must be positive.")
    if Tc_K <= 0.0:
        raise ValueError("Tc_K must be positive.")

    prefactor = 96.0 * ZETA_5 * N0_J_m3 * KB_J_K**2 / (tau0_s * Tc_K**3)
    return float(prefactor * (Te_K**5 - Tph_K**5))


def cumulative_spectral_support(result: ProjectedPowerResult) -> dict[str, np.ndarray]:
    """Return cumulative support curves for the OE5 spectral integrands.

    These curves are not powers; they diagnose which Omega intervals dominate
    the projected power integrals.
    """
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