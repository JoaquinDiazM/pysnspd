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

OE5-v4 policy
-------------
OE5-v4 removes the earlier fixed-gap and BCS-like diagnostic closures. It
builds a thermal self-consistent Usadel grid,

    Delta_eq(Te, q),

using the existing OE3 Matsubara Usadel self-consistency solver. Projected
powers are then evaluated as

    P_ep(Te; Delta_eq(Te,q), q).

This still is not the final gTDGL-coupled dynamics, because future PHOTON-runs
will use Delta(r,t) and q(r,t). It is, however, the correct local equilibrium
audit for OE5.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np

from pysnspd.kinetic.eliashberg import EliashbergSpectrum, j_to_mev
from pysnspd.usadel.calibration import (
    matsubara_energy_axis_J,
    solve_gap_for_gamma_J,
    solve_matsubara_s_values,
)


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
class ThermalUsadelGrid:
    """Thermal Usadel equilibrium grid Delta_eq(Te,q).

    The q-axis is inherited from the OE3 low-temperature current calibration
    branch, usually up to the low-temperature critical point. For each
    temperature Te and depairing energy Gamma_q, Delta_eq is recomputed from
    the Matsubara self-consistency equation.

    The field ``reference_current_fraction`` is I(q,T_bias)/Ic(T_bias), used
    only as a stable horizontal coordinate for diagnostics.
    """

    Te_values_K: np.ndarray
    q_values_m_inv: np.ndarray
    gamma_values_J: np.ndarray
    delta_eq_Tq_J: np.ndarray
    current_Tq_A: np.ndarray
    current_density_Tq_A_m2: np.ndarray
    current_fraction_Tq: np.ndarray
    reference_current_fraction: np.ndarray
    metadata: dict[str, Any]


def tau0_from_tau_ep_Tc(tau_ep_Tc_s: float) -> float:
    """Convert linear tau_ep(Tc) into the Vodolazov/Allmaras tau0."""
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


def build_thermal_usadel_grid(
    usadel_catalog,
    Te_values_K: np.ndarray,
    *,
    n_q: int = 160,
    n_matsubara: int = 500,
    stable_lowT_branch_only: bool = True,
) -> ThermalUsadelGrid:
    """Build Delta_eq(Te,q) using the OE3 Matsubara Usadel solver.

    Parameters
    ----------
    usadel_catalog:
        OE3 catalogue loaded from ``usadel_dos_catalog.npz``.
    Te_values_K:
        Electron-temperature axis where Delta_eq is recomputed.
    n_q:
        Number of q points in the interpolated low-temperature branch.
    n_matsubara:
        Number of Matsubara frequencies used for the gap solve.
    stable_lowT_branch_only:
        If True, use q in [0,q_c(T_bias)] from the OE3 calibration sweep.

    Notes
    -----
    This function is the OE5-v4 replacement for all fixed-gap diagnostic plots.
    It evaluates the same self-consistency equation as OE3, but at each Te.
    """
    Te_values = np.asarray(Te_values_K, dtype=float)
    if Te_values.ndim != 1 or Te_values.size < 2:
        raise ValueError("Te_values_K must be a one-dimensional array with >=2 values.")
    if np.any(Te_values <= 0.0):
        raise ValueError("All Te values must be positive.")
    if np.any(np.diff(Te_values) < 0.0):
        raise ValueError("Te_values_K must be sorted.")

    D_m2_s = float(usadel_catalog.metadata["D_m2_s"])
    sigma_n = float(usadel_catalog.metadata["sigma_n_S_m"])
    width_m = float(usadel_catalog.metadata["width_m"])
    thickness_m = float(usadel_catalog.metadata["thickness_m"])
    area_m2 = width_m * thickness_m
    Tc_K = float(usadel_catalog.metadata["Tc_K"])
    T_bias_K = float(usadel_catalog.metadata["T_bias_K"])

    q_cal = np.asarray(usadel_catalog.calibration_q_values_m_inv, dtype=float)
    current_cal = np.asarray(usadel_catalog.calibration_current_values_A, dtype=float)

    valid = np.isfinite(q_cal) & np.isfinite(current_cal)
    valid &= q_cal >= 0.0
    valid &= current_cal >= 0.0
    if np.sum(valid) < 3:
        raise ValueError("Invalid OE3 calibration branch.")

    q_cal = q_cal[valid]
    current_cal = current_cal[valid]
    order = np.argsort(q_cal)
    q_cal = q_cal[order]
    current_cal = current_cal[order]

    idx_ic = int(np.argmax(current_cal))
    Ic_bias_A = float(current_cal[idx_ic])
    q_c_bias = float(q_cal[idx_ic])

    if stable_lowT_branch_only:
        q_min = 0.0
        q_max = q_c_bias
    else:
        q_min = float(np.min(q_cal))
        q_max = float(np.max(q_cal))

    if n_q <= 1:
        q_values = np.asarray([q_max], dtype=float)
    else:
        q_values = np.linspace(q_min, q_max, int(n_q))

    gamma_values = 0.5 * HBAR_J_S * D_m2_s * q_values * q_values

    current_reference = np.interp(q_values, q_cal, current_cal)
    if Ic_bias_A > 0.0:
        reference_fraction = current_reference / Ic_bias_A
    else:
        reference_fraction = np.zeros_like(q_values)

    delta_grid = np.zeros((Te_values.size, q_values.size), dtype=float)
    current_grid = np.zeros_like(delta_grid)
    current_density_grid = np.zeros_like(delta_grid)
    current_fraction_grid = np.zeros_like(delta_grid)

    for iT, Te in enumerate(Te_values):
        eps_n = matsubara_energy_axis_J(T_K=float(Te), n_matsubara=int(n_matsubara))

        for iq, (q, gamma) in enumerate(zip(q_values, gamma_values, strict=True)):
            delta = solve_gap_for_gamma_J(
                gamma_J=float(gamma),
                T_K=float(Te),
                Tc_K=Tc_K,
                eps_n_J=eps_n,
            )
            delta_grid[iT, iq] = delta

            if delta > 0.0 and q > 0.0:
                s = solve_matsubara_s_values(
                    delta_J=float(delta),
                    gamma_J=float(gamma),
                    eps_n_J=eps_n,
                )
                sum_s2 = float(np.sum(s * s))
                current = (
                    area_m2
                    * (2.0 * np.pi * KB_J_K * float(Te) / E_CHARGE_C)
                    * sigma_n
                    * float(q)
                    * sum_s2
                )
            else:
                current = 0.0

            current_grid[iT, iq] = current
            current_density_grid[iT, iq] = current / area_m2

        row_max = float(np.max(current_grid[iT, :]))
        if row_max > 0.0:
            current_fraction_grid[iT, :] = current_grid[iT, :] / row_max

    metadata = {
        "backend": "thermal_usadel_grid_oe5_v4",
        "description": (
            "Delta_eq(Te,q) recomputed with the OE3 Matsubara Usadel "
            "self-consistency solver. This replaces fixed-gap OE5 diagnostics."
        ),
        "n_Te": int(Te_values.size),
        "n_q": int(q_values.size),
        "n_matsubara": int(n_matsubara),
        "Tc_K": Tc_K,
        "T_bias_K": T_bias_K,
        "D_m2_s": D_m2_s,
        "sigma_n_S_m": sigma_n,
        "width_m": width_m,
        "thickness_m": thickness_m,
        "Ic_bias_A": Ic_bias_A,
        "q_c_bias_m_inv": q_c_bias,
        "stable_lowT_branch_only": bool(stable_lowT_branch_only),
    }

    return ThermalUsadelGrid(
        Te_values_K=Te_values,
        q_values_m_inv=q_values,
        gamma_values_J=gamma_values,
        delta_eq_Tq_J=delta_grid,
        current_Tq_A=current_grid,
        current_density_Tq_A_m2=current_density_grid,
        current_fraction_Tq=current_fraction_grid,
        reference_current_fraction=reference_fraction,
        metadata=metadata,
    )


def select_thermal_usadel_q_state(
    grid: ThermalUsadelGrid,
    target_reference_current_fraction: float,
) -> dict[str, float]:
    """Select a fixed q state by low-temperature reference I/Ic."""
    if target_reference_current_fraction < 0.0:
        raise ValueError("target_reference_current_fraction must be non-negative.")

    frac = np.asarray(grid.reference_current_fraction, dtype=float)
    idx = int(np.argmin(np.abs(frac - float(target_reference_current_fraction))))

    return {
        "q_index": idx,
        "reference_current_fraction": float(frac[idx]),
        "q_m_inv": float(grid.q_values_m_inv[idx]),
        "gamma_J": float(grid.gamma_values_J[idx]),
        "gamma_meV": float(grid.gamma_values_J[idx] / MEV_J),
    }


def thermal_usadel_delta_at_state(
    grid: ThermalUsadelGrid,
    *,
    Te_K: float,
    q_m_inv: float,
) -> float:
    """Bilinearly interpolate Delta_eq(Te,q) from a ThermalUsadelGrid."""
    return float(
        _interp2(
            grid.Te_values_K,
            grid.q_values_m_inv,
            grid.delta_eq_Tq_J,
            float(Te_K),
            float(q_m_inv),
        )
    )


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
    """Compute OE5 projected powers for one local state."""
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
        "backend": "projected_powers_oe5_v4_thermal_usadel",
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


def compute_power_curve_thermal_usadel_state(
    Te_values_K: np.ndarray,
    *,
    Tph_K: float,
    state: dict[str, float],
    thermal_grid: ThermalUsadelGrid,
    phase_space_catalog,
    spectrum: EliashbergSpectrum,
    N0_J_m3: float,
    tau0_s: float,
    Tc_K: float,
    omega_max_meV: float | None = None,
) -> dict[str, np.ndarray]:
    """Compute projected powers versus Te at fixed q with Delta_eq(Te,q)."""
    Te_values = np.asarray(Te_values_K, dtype=float)

    delta_values = np.zeros_like(Te_values)
    q_values = np.full_like(Te_values, float(state["q_m_inv"]))
    P_S = np.zeros_like(Te_values)
    P_R = np.zeros_like(Te_values)
    P_total = np.zeros_like(Te_values)
    P_D = np.zeros_like(Te_values)
    P_N = np.zeros_like(Te_values)

    for i, Te in enumerate(Te_values):
        delta = thermal_usadel_delta_at_state(
            thermal_grid,
            Te_K=float(Te),
            q_m_inv=float(state["q_m_inv"]),
        )
        delta_values[i] = delta

        result = compute_projected_powers(
            float(Te),
            float(Tph_K),
            float(delta),
            float(state["q_m_inv"]),
            phase_space_catalog,
            spectrum,
            N0_J_m3=float(N0_J_m3),
            omega_max_meV=omega_max_meV,
        )
        normal = compute_projected_powers(
            float(Te),
            float(Tph_K),
            0.0,
            0.0,
            phase_space_catalog,
            spectrum,
            N0_J_m3=float(N0_J_m3),
            omega_max_meV=omega_max_meV,
        )

        P_S[i] = result.P_S_W_m3
        P_R[i] = result.P_R_W_m3
        P_total[i] = result.P_total_W_m3
        P_N[i] = normal.P_S_W_m3
        P_D[i] = compute_vodolazov_debye_power_density(
            float(Te),
            float(Tph_K),
            N0_J_m3=float(N0_J_m3),
            tau0_s=float(tau0_s),
            Tc_K=float(Tc_K),
        )

    return {
        "Te_values_K": Te_values,
        "delta_values_J": delta_values,
        "q_values_m_inv": q_values,
        "P_S_W_m3": P_S,
        "P_R_W_m3": P_R,
        "P_total_W_m3": P_total,
        "P_normal_Eliashberg_W_m3": P_N,
        "P_Debye_Vodolazov_W_m3": P_D,
    }


def compute_power_scan_thermal_usadel(
    Te_values_K: np.ndarray,
    *,
    Tph_K: float,
    thermal_grid: ThermalUsadelGrid,
    phase_space_catalog,
    spectrum: EliashbergSpectrum,
    N0_J_m3: float,
    omega_max_meV: float | None = None,
) -> dict[str, np.ndarray]:
    """Compute projected powers on the full Delta_eq(Te,q) grid."""
    Te_values = np.asarray(Te_values_K, dtype=float)
    q_values = np.asarray(thermal_grid.q_values_m_inv, dtype=float)

    shape = (Te_values.size, q_values.size)
    delta = np.zeros(shape, dtype=float)
    P_S = np.zeros(shape, dtype=float)
    P_R = np.zeros(shape, dtype=float)
    P_total = np.zeros(shape, dtype=float)

    for iT, Te in enumerate(Te_values):
        for iq, q in enumerate(q_values):
            d = thermal_usadel_delta_at_state(
                thermal_grid,
                Te_K=float(Te),
                q_m_inv=float(q),
            )
            delta[iT, iq] = d

            result = compute_projected_powers(
                float(Te),
                float(Tph_K),
                float(d),
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
        "gamma_values_J": np.asarray(thermal_grid.gamma_values_J, dtype=float),
        "reference_current_fraction": np.asarray(
            thermal_grid.reference_current_fraction,
            dtype=float,
        ),
        "delta_values_J": delta,
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


def _interp2(
    x_axis: np.ndarray,
    y_axis: np.ndarray,
    values_xy: np.ndarray,
    x: float,
    y: float,
) -> float:
    x0, x1, wx = _bracket_weight(x_axis, x)
    y0, y1, wy = _bracket_weight(y_axis, y)

    v00 = values_xy[x0, y0]
    v01 = values_xy[x0, y1]
    v10 = values_xy[x1, y0]
    v11 = values_xy[x1, y1]

    v0 = v00 * (1.0 - wy) + v01 * wy
    v1 = v10 * (1.0 - wy) + v11 * wy
    return float(v0 * (1.0 - wx) + v1 * wx)


def _cumulative_fraction(x: np.ndarray, weight: np.ndarray) -> np.ndarray:
    x = np.asarray(x, dtype=float)
    w = np.maximum(np.asarray(weight, dtype=float), 0.0)
    increments = 0.5 * (w[1:] + w[:-1]) * np.diff(x)
    cumulative = np.concatenate([[0.0], np.cumsum(increments)])
    total = float(cumulative[-1])
    if total <= 0.0:
        return np.zeros_like(cumulative)
    return cumulative / total