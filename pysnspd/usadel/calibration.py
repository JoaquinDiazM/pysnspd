"""
Critical-current calibration for the pySNSPD Usadel block.

This module estimates the electronic diffusion coefficient D from a user-given
critical current Ic_target_A, the nanowire dimensions, sigma_n and Tc.

The calculation uses the dirty-limit Matsubara Usadel equation with depairing.
It does not use D as an input. Instead, the current sweep is written in terms
of the depairing energy Gamma, and the scaling Ic(D) ~ 1/sqrt(D) is used to
obtain D.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Mapping
import math

import numpy as np

from pysnspd.config import validate_config
from pysnspd.usadel.parameters import (
    E_CHARGE_C,
    HBAR_J_S,
    K_B_J_K,
    bcs_gap_zero_J,
)


@dataclass(frozen=True)
class DiffusionCalibrationResult:
    """Result of the critical-current calibration."""
    D_m2_s: float
    Ic_target_A: float
    jc_target_A_m2: float
    Ic_model_A: float
    jc_model_A_m2: float
    q_critical_m_inv: float
    gamma_critical_J: float
    delta_critical_J: float
    gamma_values_J: np.ndarray
    q_values_m_inv: np.ndarray
    delta_eq_values_J: np.ndarray
    current_values_A: np.ndarray
    current_density_values_A_m2: np.ndarray
    sum_s2_values: np.ndarray
    warnings: list[str]

    @property
    def gamma_critical_meV(self) -> float:
        """Critical depairing energy in meV."""
        return float(self.gamma_critical_J / 1.602176634e-22)

    @property
    def delta_critical_meV(self) -> float:
        """Critical gap in meV."""
        return float(self.delta_critical_J / 1.602176634e-22)


def calibrate_diffusion_from_config(
    config: Mapping[str, Any],
) -> DiffusionCalibrationResult:
    """
    Calibrate D from ``calibration.Ic_target_A`` in the project config.
    """
    cfg = validate_config(config, require_big_data_root_exists=False)

    Tc_K = float(cfg["material"]["Tc_K"])
    T_K = float(cfg["bias"]["T_bias_K"])
    sigma_n = float(cfg["material"]["sigma_n_S_m"])
    width_m = float(cfg["material"]["width_m"])
    thickness_m = float(cfg["material"]["thickness_m"])
    area_m2 = width_m * thickness_m

    Ic_target_A = float(cfg["calibration"]["Ic_target_A"])
    jc_target_A_m2 = Ic_target_A / area_m2

    n_matsubara = int(cfg["catalogs"]["dos"]["n_matsubara"])
    n_gamma = int(cfg["calibration"]["n_gamma_sweep"])
    gamma_max_fraction = float(cfg["calibration"]["gamma_max_fraction"])

    D_warn_min = float(cfg["calibration"]["D_warn_min_m2_s"])
    D_warn_max = float(cfg["calibration"]["D_warn_max_m2_s"])

    delta0_J = bcs_gap_zero_J(Tc_K)

    gamma_values_J = np.linspace(
        0.0,
        gamma_max_fraction * delta0_J,
        n_gamma,
    )

    eps_n_J = matsubara_energy_axis_J(T_K=T_K, n_matsubara=n_matsubara)

    delta_eq = np.zeros_like(gamma_values_J)
    sum_s2 = np.zeros_like(gamma_values_J)

    for i, gamma in enumerate(gamma_values_J):
        delta_i = solve_gap_for_gamma_J(
            gamma_J=float(gamma),
            T_K=T_K,
            Tc_K=Tc_K,
            eps_n_J=eps_n_J,
        )
        delta_eq[i] = delta_i

        if delta_i > 0.0:
            s = solve_matsubara_s_values(
                delta_J=delta_i,
                gamma_J=float(gamma),
                eps_n_J=eps_n_J,
            )
            sum_s2[i] = float(np.sum(s * s))

    current_prefactor = current_prefactor_A_sqrt_D(
        gamma_values_J=gamma_values_J,
        sum_s2_values=sum_s2,
        T_K=T_K,
        sigma_n_S_m=sigma_n,
        cross_section_m2=area_m2,
    )

    idx_critical = int(np.argmax(current_prefactor))
    max_prefactor = float(current_prefactor[idx_critical])

    if max_prefactor <= 0.0:
        raise RuntimeError(
            "Could not calibrate D: Usadel Matsubara sweep produced zero critical current."
        )

    D_m2_s = (max_prefactor / Ic_target_A) ** 2

    q_values_m_inv = np.sqrt(
        np.maximum(0.0, 2.0 * gamma_values_J / (HBAR_J_S * D_m2_s))
    )

    current_values_A = current_prefactor / math.sqrt(D_m2_s)
    current_density_values_A_m2 = current_values_A / area_m2

    warnings = []

    if D_m2_s < D_warn_min:
        warnings.append(
            "Calibrated D is below the configured NbN warning range: "
            f"D={D_m2_s:.4e} m^2/s < {D_warn_min:.4e} m^2/s."
        )

    if D_m2_s > D_warn_max:
        warnings.append(
            "Calibrated D is above the configured NbN warning range: "
            f"D={D_m2_s:.4e} m^2/s > {D_warn_max:.4e} m^2/s."
        )

    if idx_critical == len(gamma_values_J) - 1:
        warnings.append(
            "Critical current maximum occurred at the largest Gamma in the sweep. "
            "Increase calibration.gamma_max_fraction."
        )

    Ic_model_A = float(current_values_A[idx_critical])
    jc_model_A_m2 = float(current_density_values_A_m2[idx_critical])

    return DiffusionCalibrationResult(
        D_m2_s=float(D_m2_s),
        Ic_target_A=float(Ic_target_A),
        jc_target_A_m2=float(jc_target_A_m2),
        Ic_model_A=Ic_model_A,
        jc_model_A_m2=jc_model_A_m2,
        q_critical_m_inv=float(q_values_m_inv[idx_critical]),
        gamma_critical_J=float(gamma_values_J[idx_critical]),
        delta_critical_J=float(delta_eq[idx_critical]),
        gamma_values_J=gamma_values_J,
        q_values_m_inv=q_values_m_inv,
        delta_eq_values_J=delta_eq,
        current_values_A=current_values_A,
        current_density_values_A_m2=current_density_values_A_m2,
        sum_s2_values=sum_s2,
        warnings=warnings,
    )


def matsubara_energy_axis_J(
    *,
    T_K: float,
    n_matsubara: int,
) -> np.ndarray:
    """
    Return positive fermionic Matsubara energies.

    epsilon_n = hbar omega_n = pi k_B T (2n+1).
    """
    if T_K <= 0.0:
        raise ValueError("T_K must be positive.")

    if n_matsubara <= 0:
        raise ValueError("n_matsubara must be positive.")

    n = np.arange(int(n_matsubara), dtype=float)
    return math.pi * K_B_J_K * T_K * (2.0 * n + 1.0)


def solve_matsubara_s_values(
    *,
    delta_J: float,
    gamma_J: float,
    eps_n_J: np.ndarray,
    n_iter: int = 80,
) -> np.ndarray:
    """
    Solve the Matsubara Usadel equation for s_n = sin(theta_n).

    We solve

        Delta sqrt(1-s^2) = (eps_n + Gamma sqrt(1-s^2)) s

    by vectorized bisection on 0 <= s <= 1.
    """
    eps = np.asarray(eps_n_J, dtype=float)

    if delta_J <= 0.0:
        return np.zeros_like(eps)

    if gamma_J < 0.0:
        raise ValueError("gamma_J must be nonnegative.")

    lo = np.zeros_like(eps)
    hi = np.ones_like(eps) * (1.0 - 1.0e-14)

    for _ in range(n_iter):
        mid = 0.5 * (lo + hi)
        c = np.sqrt(np.maximum(0.0, 1.0 - mid * mid))
        f = delta_J * c - (eps + gamma_J * c) * mid

        lo = np.where(f > 0.0, mid, lo)
        hi = np.where(f > 0.0, hi, mid)

    return 0.5 * (lo + hi)


def self_consistency_residual_J(
    *,
    delta_J: float,
    gamma_J: float,
    T_K: float,
    Tc_K: float,
    eps_n_J: np.ndarray,
) -> float:
    """
    BCS self-consistency residual for fixed Gamma.
    """
    if delta_J <= 0.0:
        return 0.0

    s = solve_matsubara_s_values(
        delta_J=delta_J,
        gamma_J=gamma_J,
        eps_n_J=eps_n_J,
    )

    residual = (
        delta_J * math.log(T_K / Tc_K)
        +
        2.0
        * math.pi
        * K_B_J_K
        * T_K
        * float(np.sum(delta_J / eps_n_J - s))
    )

    return float(residual)


def solve_gap_for_gamma_J(
    *,
    gamma_J: float,
    T_K: float,
    Tc_K: float,
    eps_n_J: np.ndarray,
    n_scan: int = 80,
    n_bisect: int = 80,
) -> float:
    """
    Solve the nonzero self-consistent gap for fixed depairing Gamma.

    The trivial Delta=0 solution is ignored. If no nonzero superconducting
    solution exists, return 0.
    """
    if T_K >= Tc_K:
        return 0.0

    delta0 = bcs_gap_zero_J(Tc_K)
    delta_min = delta0 * 1.0e-8
    delta_max = delta0 * 1.05

    grid = np.linspace(delta_min, delta_max, n_scan)
    values = np.array(
        [
            self_consistency_residual_J(
                delta_J=float(delta),
                gamma_J=gamma_J,
                T_K=T_K,
                Tc_K=Tc_K,
                eps_n_J=eps_n_J,
            )
            for delta in grid
        ],
        dtype=float,
    )

    bracket = None
    for i in range(len(grid) - 1):
        if values[i] <= 0.0 and values[i + 1] >= 0.0:
            bracket = (float(grid[i]), float(grid[i + 1]))
            break

    if bracket is None:
        return 0.0

    lo, hi = bracket

    for _ in range(n_bisect):
        mid = 0.5 * (lo + hi)
        f_mid = self_consistency_residual_J(
            delta_J=mid,
            gamma_J=gamma_J,
            T_K=T_K,
            Tc_K=Tc_K,
            eps_n_J=eps_n_J,
        )

        if f_mid <= 0.0:
            lo = mid
        else:
            hi = mid

    return float(0.5 * (lo + hi))


def current_prefactor_A_sqrt_D(
    *,
    gamma_values_J: np.ndarray,
    sum_s2_values: np.ndarray,
    T_K: float,
    sigma_n_S_m: float,
    cross_section_m2: float,
) -> np.ndarray:
    """
    Return the current prefactor A(D=1)*sqrt(D).

    The supercurrent density is

        j_s = (2 pi k_B T / e) sigma_n k sum_n s_n^2,

    and

        Gamma = hbar D k^2 / 2.

    Hence

        I(Gamma,D) = prefactor(Gamma) / sqrt(D).
    """
    gamma = np.asarray(gamma_values_J, dtype=float)
    sum_s2 = np.asarray(sum_s2_values, dtype=float)

    k_prefactor = np.sqrt(np.maximum(0.0, 2.0 * gamma / HBAR_J_S))

    return (
        cross_section_m2
        * (2.0 * math.pi * K_B_J_K * T_K / E_CHARGE_C)
        * sigma_n_S_m
        * k_prefactor
        * sum_s2
    )


def calibration_summary(result: DiffusionCalibrationResult) -> dict[str, Any]:
    """
    Convert a calibration result to a compact manifest-friendly dictionary.
    """
    return {
        "D_m2_s": float(result.D_m2_s),
        "Ic_target_A": float(result.Ic_target_A),
        "jc_target_A_m2": float(result.jc_target_A_m2),
        "Ic_model_A": float(result.Ic_model_A),
        "jc_model_A_m2": float(result.jc_model_A_m2),
        "q_critical_m_inv": float(result.q_critical_m_inv),
        "gamma_critical_J": float(result.gamma_critical_J),
        "gamma_critical_meV": float(result.gamma_critical_meV),
        "delta_critical_J": float(result.delta_critical_J),
        "delta_critical_meV": float(result.delta_critical_meV),
        "warnings": list(result.warnings),
    }