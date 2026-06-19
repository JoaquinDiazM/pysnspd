"""
Usadel/DOS catalogue construction for pySNSPD.

OE3 v2:
- Build a quasistatic spectral catalogue over |Delta| and depairing energy.
- The DOS backend is now the real-axis uniform dirty-limit Usadel quartic.
- The catalogue interface remains compatible with OE3 v1.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any, Mapping

import numpy as np

from pysnspd.config import validate_config
from pysnspd.usadel.parameters import (
    material_parameters_from_config,
    depairing_energy_grid_J,
    q_axis_from_depairing_energy_m_inv,
    energy_axis_J,
)
from pysnspd.usadel.solver import compute_dos_grid, dos_diagnostics


@dataclass(frozen=True)
class UsadelCatalog:
    """
    Container for the OE3 quasistatic spectral catalogue.
    """
    energy_values_J: np.ndarray
    delta_values_J: np.ndarray
    gamma_values_J: np.ndarray
    q_values_m_inv: np.ndarray
    rho_delta_gamma_E: np.ndarray
    anomalous_delta_gamma_E: np.ndarray
    eta_J: float
    metadata: dict[str, Any]

    @property
    def shape(self) -> tuple[int, int, int]:
        """Return catalogue shape ``(n_delta, n_gamma, n_energy)``."""
        return tuple(int(v) for v in self.rho_delta_gamma_E.shape)


def build_usadel_catalog_from_config(
    config: Mapping[str, Any],
    *,
    eta_fraction: float = 1.0e-3,
    gamma_max_fraction: float = 0.35,
    energy_max_factor: float = 6.0,
) -> UsadelCatalog:
    """
    Build the OE3 v2 DOS catalogue from the project config.

    This function does not yet solve the Matsubara self-consistency problem for
    ``Delta(q,T)`` or ``j_s(q,T)``. Instead, it builds the real-axis spectral
    catalogue over independent ``|Delta|`` and depairing ``Gamma_q`` axes.
    """
    cfg = validate_config(config, require_big_data_root_exists=False)
    mat = material_parameters_from_config(cfg)

    dos_cfg = cfg["catalogs"]["dos"]
    n_delta = int(dos_cfg["n_delta"])
    n_q = int(dos_cfg["n_q"])
    n_energy = int(dos_cfg["n_energy"])
    n_matsubara = int(dos_cfg["n_matsubara"])

    delta_ref_J = mat.delta0_J
    delta_bias_J = mat.delta_bias_J

    delta_values_J = np.linspace(0.0, delta_ref_J, n_delta)

    gamma_values_J = depairing_energy_grid_J(
        delta_ref_J=delta_ref_J,
        n_q=n_q,
        gamma_max_fraction=gamma_max_fraction,
    )

    q_values_m_inv = q_axis_from_depairing_energy_m_inv(
        gamma_values_J,
        D_m2_s=mat.D_m2_s,
    )

    energy_values = energy_axis_J(
        delta_ref_J=delta_ref_J,
        n_energy=n_energy,
        energy_max_factor=energy_max_factor,
    )

    eta_J = eta_fraction * delta_ref_J

    rho, anomalous = compute_dos_grid(
        energy_values,
        delta_values_J,
        gamma_values_J,
        eta_J=eta_J,
    )

    metadata = {
        "backend": "uniform_dirty_usadel_quartic_oe3_v2",
        "description": (
            "OE3 v2 spectral catalogue. The DOS is computed from the real-axis "
            "uniform dirty-limit Usadel quartic with depairing Gamma_q. "
            "The catalogue is built over independent |Delta| and Gamma_q axes. "
            "The Matsubara self-consistent Delta(q,T), j_s(q,T), and Ic sweep "
            "are not implemented in this backend yet."
        ),
        "spectral_equation": "(Gamma*c - i*z)^2*(1-c^2) - Delta^2*c^2 = 0",
        "rho_definition": "rho(E)=Re[c(E)]",
        "gamma_definition": "Gamma_q = hbar*D*q^2/2",
        "material": mat.name,
        "Tc_K": mat.Tc_K,
        "T_bias_K": mat.T_bias_K,
        "D_m2_s": mat.D_m2_s,
        "sigma_n_S_m": mat.sigma_n_S_m,
        "width_m": mat.width_m,
        "thickness_m": mat.thickness_m,
        "I_bias_A": mat.I_bias_A,
        "bias_current_density_A_m2": mat.bias_current_density_A_m2,
        "delta0_J": mat.delta0_J,
        "delta0_meV": J_to_meV(mat.delta0_J),
        "delta_bias_J": delta_bias_J,
        "delta_bias_meV": J_to_meV(delta_bias_J),
        "eta_J": eta_J,
        "eta_fraction": eta_fraction,
        "gamma_max_fraction": gamma_max_fraction,
        "energy_max_factor": energy_max_factor,
        "n_delta": n_delta,
        "n_q": n_q,
        "n_energy": n_energy,
        "n_matsubara_configured": n_matsubara,
    }

    return UsadelCatalog(
        energy_values_J=energy_values,
        delta_values_J=delta_values_J,
        gamma_values_J=gamma_values_J,
        q_values_m_inv=q_values_m_inv,
        rho_delta_gamma_E=rho,
        anomalous_delta_gamma_E=anomalous,
        eta_J=float(eta_J),
        metadata=metadata,
    )


def catalog_summary(catalog: UsadelCatalog) -> dict[str, Any]:
    """
    Build a compact summary for manifests and console output.
    """
    diag = dos_diagnostics(
        catalog.rho_delta_gamma_E,
        catalog.energy_values_J,
    )

    summary = {
        "backend": str(catalog.metadata.get("backend", "unknown")),
        "shape": list(catalog.shape),
        "n_delta": int(catalog.delta_values_J.size),
        "n_q": int(catalog.gamma_values_J.size),
        "n_energy": int(catalog.energy_values_J.size),
        "delta_min_J": float(np.min(catalog.delta_values_J)),
        "delta_max_J": float(np.max(catalog.delta_values_J)),
        "delta_max_meV": J_to_meV(float(np.max(catalog.delta_values_J))),
        "gamma_min_J": float(np.min(catalog.gamma_values_J)),
        "gamma_max_J": float(np.max(catalog.gamma_values_J)),
        "gamma_max_meV": J_to_meV(float(np.max(catalog.gamma_values_J))),
        "q_min_m_inv": float(np.min(catalog.q_values_m_inv)),
        "q_max_m_inv": float(np.max(catalog.q_values_m_inv)),
        "energy_min_J": float(np.min(catalog.energy_values_J)),
        "energy_max_J": float(np.max(catalog.energy_values_J)),
        "energy_max_meV": J_to_meV(float(np.max(catalog.energy_values_J))),
        "eta_J": float(catalog.eta_J),
        "eta_meV": J_to_meV(float(catalog.eta_J)),
    }
    summary.update(diag)

    return summary


def save_usadel_catalog_npz(catalog: UsadelCatalog, path: str | Path) -> Path:
    """
    Save a Usadel catalogue to compressed ``.npz``.
    """
    output = Path(path)
    output.parent.mkdir(parents=True, exist_ok=True)

    np.savez_compressed(
        output,
        energy_values_J=catalog.energy_values_J,
        delta_values_J=catalog.delta_values_J,
        gamma_values_J=catalog.gamma_values_J,
        q_values_m_inv=catalog.q_values_m_inv,
        rho_delta_gamma_E=catalog.rho_delta_gamma_E,
        anomalous_delta_gamma_E=catalog.anomalous_delta_gamma_E,
        eta_J=np.array(catalog.eta_J),
        metadata=np.array(catalog.metadata, dtype=object),
    )

    return output


def load_usadel_catalog_npz(path: str | Path) -> UsadelCatalog:
    """
    Load a catalogue saved by :func:`save_usadel_catalog_npz`.
    """
    source = Path(path)

    with np.load(source, allow_pickle=True) as data:
        metadata = data["metadata"].item()
        return UsadelCatalog(
            energy_values_J=np.asarray(data["energy_values_J"], dtype=float),
            delta_values_J=np.asarray(data["delta_values_J"], dtype=float),
            gamma_values_J=np.asarray(data["gamma_values_J"], dtype=float),
            q_values_m_inv=np.asarray(data["q_values_m_inv"], dtype=float),
            rho_delta_gamma_E=np.asarray(data["rho_delta_gamma_E"], dtype=float),
            anomalous_delta_gamma_E=np.asarray(data["anomalous_delta_gamma_E"], dtype=float),
            eta_J=float(data["eta_J"]),
            metadata=dict(metadata),
        )


def J_to_meV(value_J: float) -> float:
    """
    Convert Joules to meV.
    """
    return float(value_J / 1.602176634e-22)


def meV_axis(values_J: np.ndarray) -> np.ndarray:
    """
    Convert a Joule energy axis to meV.
    """
    return np.asarray(values_J, dtype=float) / 1.602176634e-22