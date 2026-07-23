"""Material closures for the pySNSPD gTDGL sector.

This module contains only local material functions used by the mesoscopic
gTDGL/Poisson solver. It deliberately does not evolve the thermal model or the
external circuit.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Mapping
import math

import numpy as np

E_CHARGE_C = 1.602176634e-19
K_B_J_K = 1.380649e-23
HBAR_J_S = 1.054571817e-34


@dataclass(frozen=True)
class GTDGLMaterial:
    """Local material parameters used by the gTDGL stationary relaxation."""

    Tc_K: float
    D_m2_s: float
    sigma_n_S_m: float
    delta0_J: float
    thickness_m: float
    width_m: float
    tau_ee_Tc_s: float
    tau_ep_Tc_s: float
    D_base_m2_s: float | None = None
    D_effective_factor: float = 1.0

    @property
    def tau0_GL_s(self) -> float:
        """Return the GL time scale pi*hbar/(8*kB*Tc)."""
        return math.pi * HBAR_J_S / (8.0 * K_B_J_K * self.Tc_K)

    def tau_ee_s(self, Te_K: np.ndarray | float) -> np.ndarray:
        """Electron-electron relaxation time from the YAML value at Tc."""
        Te = _positive_temperature(Te_K)
        return self.tau_ee_Tc_s * (self.Tc_K / Te)

    def tau_ep_s(self, Te_K: np.ndarray | float) -> np.ndarray:
        """Electron-phonon relaxation time from the YAML value at Tc."""
        Te = _positive_temperature(Te_K)
        return self.tau_ep_Tc_s * (self.Tc_K / Te) ** 3

    def tau_sc_s(self, Te_K: np.ndarray | float) -> np.ndarray:
        """Combined superconducting relaxation time."""
        tee = self.tau_ee_s(Te_K)
        tep = self.tau_ep_s(Te_K)
        return 1.0 / (1.0 / tee + 1.0 / tep)

    def rho_kwt(self, Te_K: np.ndarray | float, R_J: np.ndarray | float) -> np.ndarray:
        """Kramer-Watts-Tobin relaxation factor."""
        R = np.maximum(np.asarray(R_J, dtype=float), 0.0)
        tau = self.tau_sc_s(Te_K)
        return np.sqrt(1.0 + 4.0 * R**2 * tau**2 / HBAR_J_S**2)


    def xi_mod_squared_m2(self, Te_K: np.ndarray | float) -> np.ndarray:
        """Modified coherence-length scale appearing in Appendix B."""
        Te = _positive_temperature(Te_K)
        return (
            math.pi
            * HBAR_J_S
            * self.D_m2_s
            / (4.0 * math.sqrt(2.0) * K_B_J_K * self.Tc_K * np.sqrt(1.0 + Te / self.Tc_K))
        )

    def delta_mod_squared_J2(self, Te_K: np.ndarray | float) -> np.ndarray:
        """GL-compatible gap scale Delta_mod^2(Te)."""
        Te = _positive_temperature(Te_K)
        out = np.empty_like(Te, dtype=float)
        below = Te < self.Tc_K
        one_minus_t = np.maximum(1.0 - Te / self.Tc_K, 1.0e-12)

        out[below] = bcs_gap_J_array(Te[below], self.Tc_K, self.delta0_J) ** 2 / one_minus_t[below]
        out[~below] = self.delta0_J**2
        return np.maximum(out, (1.0e-12 * self.delta0_J) ** 2)


def build_gtdgl_material(
    config: Mapping[str, Any],
    usadel_catalog: Any,
    *,
    diffusion_factor: float = 1.0,
) -> GTDGLMaterial:
    """Build gTDGL material parameters from the project config and OE3 catalogue."""

    material_cfg = dict(config.get("material", {}))
    metadata = dict(getattr(usadel_catalog, "metadata", {}) or {})

    Tc_K = float(metadata.get("Tc_K", material_cfg.get("Tc_K")))
    D_base_m2_s = float(metadata.get("D_m2_s", material_cfg.get("D_m2_s", np.nan)))
    diffusion_factor = float(diffusion_factor)
    if not math.isfinite(diffusion_factor) or diffusion_factor <= 0.0:
        raise ValueError(f"diffusion_factor must be finite and positive, got {diffusion_factor!r}.")
    D_m2_s = D_base_m2_s * diffusion_factor
    sigma_n = float(metadata.get("sigma_n_S_m", material_cfg.get("sigma_n_S_m")))
    delta0_J = float(metadata.get("delta0_J", 1.764 * K_B_J_K * Tc_K))
    thickness_m = float(metadata.get("thickness_m", material_cfg.get("thickness_m")))
    width_m = float(metadata.get("width_m", material_cfg.get("width_m")))

    tau_ee_Tc_s = _read_time_seconds(
        material_cfg,
        names_s=("tau_ee_Tc_s", "tau_ee_s", "tau_ee_at_Tc_s"),
        names_ps=("tau_ee_Tc_ps", "tau_ee_ps", "tau_ee_at_Tc_ps"),
        default_s=5.0e-12,
    )
    tau_ep_Tc_s = _read_time_seconds(
        material_cfg,
        names_s=("tau_ep_Tc_s", "tau_ep_s", "tau_ep_at_Tc_s"),
        names_ps=("tau_ep_Tc_ps", "tau_ep_ps", "tau_ep_at_Tc_ps"),
        default_s=24.7e-12,
    )

    values = {
        "Tc_K": Tc_K,
        "D_m2_s": D_m2_s,
        "sigma_n_S_m": sigma_n,
        "delta0_J": delta0_J,
        "thickness_m": thickness_m,
        "width_m": width_m,
        "tau_ee_Tc_s": tau_ee_Tc_s,
        "tau_ep_Tc_s": tau_ep_Tc_s,
    }
    for key, value in values.items():
        if not np.isfinite(value) or value <= 0.0:
            raise ValueError(f"{key} must be finite and positive, got {value!r}.")

    return GTDGLMaterial(
        Tc_K=Tc_K,
        D_m2_s=D_m2_s,
        sigma_n_S_m=sigma_n,
        delta0_J=delta0_J,
        thickness_m=thickness_m,
        width_m=width_m,
        tau_ee_Tc_s=tau_ee_Tc_s,
        tau_ep_Tc_s=tau_ep_Tc_s,
        D_base_m2_s=D_base_m2_s,
        D_effective_factor=diffusion_factor,
    )


def bcs_gap_J_array(Te_K: np.ndarray, Tc_K: float, delta0_J: float) -> np.ndarray:
    """Weak-coupling BCS gap interpolation for Te < Tc."""
    Te = np.asarray(Te_K, dtype=float)
    out = np.zeros_like(Te, dtype=float)
    below = (Te > 0.0) & (Te < Tc_K)
    arg = np.zeros_like(Te, dtype=float)
    arg[below] = 1.74 * np.sqrt(Tc_K / Te[below] - 1.0)
    out[below] = delta0_J * np.tanh(arg[below])
    out[Te == 0.0] = delta0_J
    return out


def _positive_temperature(Te_K: np.ndarray | float) -> np.ndarray:
    Te = np.asarray(Te_K, dtype=float)
    return np.maximum(Te, 1.0e-12)


def _read_time_seconds(
    mapping: Mapping[str, Any],
    *,
    names_s: tuple[str, ...],
    names_ps: tuple[str, ...],
    default_s: float,
) -> float:
    """Read one time parameter from seconds or picoseconds aliases.

    The YAML validator normally catches invalid values first, but this helper is
    intentionally defensive because :func:`build_gtdgl_material` is also used in
    tests with small dictionaries.
    """
    found: list[tuple[str, float]] = []

    for name in names_s:
        if name in mapping:
            found.append((name, float(mapping[name])))

    for name in names_ps:
        if name in mapping:
            found.append((name, float(mapping[name]) * 1.0e-12))

    if len(found) > 1:
        names = ", ".join(name for name, _ in found)
        raise ValueError(
            f"Multiple aliases were supplied for one relaxation time ({names}). "
            "Use exactly one seconds or picoseconds key."
        )

    value = float(found[0][1]) if found else float(default_s)
    if not math.isfinite(value) or value <= 0.0:
        raise ValueError(f"Relaxation time must be finite and positive, got {value!r}.")
    return value
