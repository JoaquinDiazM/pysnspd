"""
Material and spectral parameters for the pySNSPD Usadel block.

D is not a user input. It is calibrated from the user-provided critical current.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Mapping
import math

import numpy as np

from pysnspd.config import validate_config


E_CHARGE_C = 1.602176634e-19
K_B_J_K = 1.380649e-23
HBAR_J_S = 1.054571817e-34


@dataclass(frozen=True)
class MaterialParameters:
    """Material parameters required by the Usadel catalogue."""
    name: str
    Tc_K: float
    T_bias_K: float
    D_m2_s: float
    sigma_n_S_m: float
    thickness_m: float
    width_m: float
    I_bias_A: float
    Ic_target_A: float
    jc_target_A_m2: float

    @property
    def delta0_J(self) -> float:
        """Weak-coupling zero-temperature BCS gap."""
        return bcs_gap_zero_J(self.Tc_K)

    @property
    def delta_bias_J(self) -> float:
        """Weak-coupling BCS gap at the configured bias temperature."""
        return bcs_gap_J(self.T_bias_K, self.Tc_K)

    @property
    def cross_section_m2(self) -> float:
        """Nanowire cross section."""
        return self.width_m * self.thickness_m

    @property
    def bias_current_density_A_m2(self) -> float:
        """Bias current density estimated from the configured geometry."""
        return self.I_bias_A / self.cross_section_m2


def material_parameters_from_config(config: Mapping[str, Any]) -> MaterialParameters:
    """
    Build :class:`MaterialParameters` from a validated config dictionary.

    The diffusion coefficient is calibrated from ``calibration.Ic_target_A``.
    """
    cfg = validate_config(config, require_big_data_root_exists=False)

    from pysnspd.usadel.calibration import calibrate_diffusion_from_config

    calibration = calibrate_diffusion_from_config(cfg)

    return MaterialParameters(
        name=str(cfg["material"]["name"]),
        Tc_K=float(cfg["material"]["Tc_K"]),
        T_bias_K=float(cfg["bias"]["T_bias_K"]),
        D_m2_s=float(calibration.D_m2_s),
        sigma_n_S_m=float(cfg["material"]["sigma_n_S_m"]),
        thickness_m=float(cfg["material"]["thickness_m"]),
        width_m=float(cfg["material"]["width_m"]),
        I_bias_A=float(cfg["bias"]["I_bias_A"]),
        Ic_target_A=float(calibration.Ic_target_A),
        jc_target_A_m2=float(calibration.jc_target_A_m2),
    )


def bcs_gap_zero_J(Tc_K: float) -> float:
    """
    Return the weak-coupling BCS zero-temperature gap.

    Uses ``Delta_0 = 1.764 k_B T_c``.
    """
    if Tc_K <= 0.0:
        raise ValueError("Tc_K must be positive.")

    return 1.764 * K_B_J_K * Tc_K


def bcs_gap_J(T_K: float, Tc_K: float) -> float:
    """
    Approximate weak-coupling BCS gap at temperature ``T_K``.
    """
    if Tc_K <= 0.0:
        raise ValueError("Tc_K must be positive.")

    if T_K < 0.0:
        raise ValueError("T_K must be nonnegative.")

    if T_K >= Tc_K:
        return 0.0

    if T_K == 0.0:
        return bcs_gap_zero_J(Tc_K)

    delta0 = bcs_gap_zero_J(Tc_K)
    argument = 1.74 * math.sqrt(Tc_K / T_K - 1.0)
    return float(delta0 * math.tanh(argument))


def depairing_energy_grid_J(
    *,
    delta_ref_J: float,
    n_q: int,
    gamma_max_fraction: float = 0.35,
) -> np.ndarray:
    """
    Build a depairing-energy grid.
    """
    if delta_ref_J <= 0.0:
        raise ValueError("delta_ref_J must be positive.")

    if n_q <= 0:
        raise ValueError("n_q must be positive.")

    if gamma_max_fraction < 0.0:
        raise ValueError("gamma_max_fraction must be nonnegative.")

    gamma_max = gamma_max_fraction * delta_ref_J
    return np.linspace(0.0, gamma_max, int(n_q))


def q_axis_from_depairing_energy_m_inv(
    gamma_values_J: np.ndarray,
    *,
    D_m2_s: float,
) -> np.ndarray:
    """
    Convert depairing energies to phase-gradient values.

        Gamma_q = hbar D q^2 / 2.
    """
    if D_m2_s <= 0.0:
        raise ValueError("D_m2_s must be positive.")

    gamma = np.asarray(gamma_values_J, dtype=float)
    gamma = np.maximum(gamma, 0.0)

    return np.sqrt(2.0 * gamma / (HBAR_J_S * D_m2_s))


def energy_axis_J(
    *,
    delta_ref_J: float,
    n_energy: int,
    energy_max_factor: float = 6.0,
) -> np.ndarray:
    """
    Build the quasiparticle energy axis for the DOS catalogue.
    """
    if delta_ref_J <= 0.0:
        raise ValueError("delta_ref_J must be positive.")

    if n_energy < 2:
        raise ValueError("n_energy must be at least 2.")

    if energy_max_factor <= 0.0:
        raise ValueError("energy_max_factor must be positive.")

    return np.linspace(0.0, energy_max_factor * delta_ref_J, int(n_energy))