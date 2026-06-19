"""
Phase-space catalogues for the pySNSPD kinetic block.

OE4 first implementation:
- Build tabulated phase-space integrals J_S and J_R from the Usadel DOS catalogue.
- Keep this independent from alpha^2 F(Omega), phonon DOS and T_ph.
- Store only the expensive electronic/superconducting part of the collision integrals.

The implemented objects are

    J_S(Omega; Te, Delta, q)

for quasiparticle-phonon scattering, and

    J_R(Omega; Te, Delta, q)

for recombination / pair-breaking phase space.

The final powers will be built later as

    P_ep^S ~ int dOmega alpha^2F(Omega) Omega [n_e - n_ph] J_S,
    P_ep^R ~ int dOmega alpha^2F(Omega) Omega [n_e - n_ph] J_R.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any, Mapping

import numpy as np

from pysnspd.config import validate_config
from pysnspd.usadel.catalog import UsadelCatalog, J_to_meV
from pysnspd.usadel.parameters import K_B_J_K


@dataclass(frozen=True)
class PhaseSpaceCatalog:
    """
    Container for phase-space integrals.

    Arrays use the shape

        (n_Te, n_delta, n_q, n_omega).

    Attributes
    ----------
    Te_values_K:
        Electron-temperature axis.
    omega_values_J:
        Phonon-energy axis.
    delta_values_J:
        Selected gap axis inherited from the Usadel catalogue.
    gamma_values_J:
        Selected depairing-energy axis inherited from the Usadel catalogue.
    q_values_m_inv:
        Selected q axis inherited from the Usadel catalogue.
    J_S_TdqO_J:
        Scattering phase-space integral.
    J_R_TdqO_J:
        Recombination / pair-breaking phase-space integral.
    delta_indices:
        Indices selected from the parent Usadel catalogue.
    q_indices:
        Indices selected from the parent Usadel catalogue.
    metadata:
        Metadata dictionary.
    """
    Te_values_K: np.ndarray
    omega_values_J: np.ndarray
    delta_values_J: np.ndarray
    gamma_values_J: np.ndarray
    q_values_m_inv: np.ndarray
    J_S_TdqO_J: np.ndarray
    J_R_TdqO_J: np.ndarray
    delta_indices: np.ndarray
    q_indices: np.ndarray
    metadata: dict[str, Any]

    @property
    def shape(self) -> tuple[int, int, int, int]:
        """Return catalogue shape ``(n_Te, n_delta, n_q, n_omega)``."""
        return tuple(int(v) for v in self.J_S_TdqO_J.shape)


def build_phase_space_catalog_from_usadel_catalog(
    usadel_catalog: UsadelCatalog,
    config: Mapping[str, Any],
    *,
    n_Te: int | None = None,
    n_delta: int | None = None,
    n_q: int | None = None,
    n_omega: int | None = None,
    Te_min_K: float | None = None,
    Te_max_K: float | None = None,
) -> PhaseSpaceCatalog:
    """
    Build a phase-space catalogue from a Usadel DOS catalogue.

    For the first OE4 attempt, it is recommended to use moderate values such as
    ``n_Te=6``, ``n_delta=6``, ``n_q=6`` and ``n_omega=160``. The full grid in
    the YAML can be much heavier and should be reserved for later production
    runs.
    """
    cfg = validate_config(config, require_big_data_root_exists=False)

    phase_cfg = cfg["catalogs"]["phase_space"]

    if n_Te is None:
        n_Te = int(phase_cfg["n_Te"])
    if n_delta is None:
        n_delta = int(phase_cfg["n_delta"])
    if n_q is None:
        n_q = int(phase_cfg["n_q"])
    if n_omega is None:
        n_omega = int(phase_cfg["n_omega"])

    if Te_min_K is None:
        Te_min_K = float(phase_cfg.get("Te_min_K", cfg["bias"]["T_bias_K"]))
    if Te_max_K is None:
        Te_max_K = float(phase_cfg.get("Te_max_K", max(4.0 * cfg["material"]["Tc_K"], 30.0)))

    if n_Te <= 0 or n_delta <= 0 or n_q <= 0 or n_omega <= 1:
        raise ValueError("Phase-space grid sizes must be positive and n_omega > 1.")

    if Te_min_K <= 0.0 or Te_max_K <= Te_min_K:
        raise ValueError("Require 0 < Te_min_K < Te_max_K.")

    delta_indices = _select_axis_indices(usadel_catalog.delta_values_J.size, n_delta)
    q_indices = _select_axis_indices(usadel_catalog.q_values_m_inv.size, n_q)

    Te_values_K = np.linspace(float(Te_min_K), float(Te_max_K), int(n_Te))
    omega_values_J = np.linspace(
        0.0,
        float(np.max(usadel_catalog.energy_values_J)),
        int(n_omega),
    )

    delta_values_J = usadel_catalog.delta_values_J[delta_indices]
    gamma_values_J = usadel_catalog.gamma_values_J[q_indices]
    q_values_m_inv = usadel_catalog.q_values_m_inv[q_indices]

    shape = (
        Te_values_K.size,
        delta_values_J.size,
        q_values_m_inv.size,
        omega_values_J.size,
    )

    J_S = np.empty(shape, dtype=float)
    J_R = np.empty(shape, dtype=float)

    energy = usadel_catalog.energy_values_J

    for iT, Te_K in enumerate(Te_values_K):
        for id_local, id_parent in enumerate(delta_indices):
            delta_J = float(usadel_catalog.delta_values_J[id_parent])

            for iq_local, iq_parent in enumerate(q_indices):
                rho_E = usadel_catalog.rho_delta_gamma_E[id_parent, iq_parent, :]

                J_S[iT, id_local, iq_local, :] = scattering_phase_space_spectrum(
                    energy,
                    rho_E,
                    omega_values_J,
                    Te_K=float(Te_K),
                    delta_J=delta_J,
                )

                J_R[iT, id_local, iq_local, :] = recombination_phase_space_spectrum(
                    energy,
                    rho_E,
                    omega_values_J,
                    Te_K=float(Te_K),
                    delta_J=delta_J,
                )

    metadata = {
        "backend": "phase_space_from_usadel_dos_oe4_v1",
        "description": (
            "OE4 first phase-space catalogue. It tabulates J_S and J_R from the "
            "Usadel DOS catalogue. It does not yet include alpha^2F(Omega), "
            "phonon DOS, T_ph, escape, or power integrals."
        ),
        "units": {
            "Te_values_K": "K",
            "omega_values_J": "J",
            "delta_values_J": "J",
            "gamma_values_J": "J",
            "q_values_m_inv": "m^-1",
            "J_S_TdqO_J": "J",
            "J_R_TdqO_J": "J",
        },
        "parent_usadel_backend": str(usadel_catalog.metadata.get("backend", "unknown")),
        "parent_usadel_shape": list(usadel_catalog.shape),
        "delta_indices": delta_indices.tolist(),
        "q_indices": q_indices.tolist(),
        "grid_is_downsampled": bool(
            delta_indices.size < usadel_catalog.delta_values_J.size
            or q_indices.size < usadel_catalog.q_values_m_inv.size
        ),
        "Te_min_K": float(Te_values_K[0]),
        "Te_max_K": float(Te_values_K[-1]),
        "omega_max_J": float(omega_values_J[-1]),
        "omega_max_meV": J_to_meV(float(omega_values_J[-1])),
    }

    return PhaseSpaceCatalog(
        Te_values_K=Te_values_K,
        omega_values_J=omega_values_J,
        delta_values_J=delta_values_J,
        gamma_values_J=gamma_values_J,
        q_values_m_inv=q_values_m_inv,
        J_S_TdqO_J=J_S,
        J_R_TdqO_J=J_R,
        delta_indices=delta_indices,
        q_indices=q_indices,
        metadata=metadata,
    )


def scattering_phase_space_spectrum(
    energy_values_J: np.ndarray,
    rho_E: np.ndarray,
    omega_values_J: np.ndarray,
    *,
    Te_K: float,
    delta_J: float,
) -> np.ndarray:
    """
    Compute J_S(Omega) for one (Te, Delta, q) slice.

    Implemented expression:

        J_S(Omega) =
            int dE rho(E) rho(E+Omega)
            [1 - Delta^2/(E(E+Omega))]
            [f(E,Te) - f(E+Omega,Te)].

    The lower integration limit follows the appendix expression. For
    ``Delta > 0`` we use ``E >= Delta`` and ``E+Omega >= Delta``. For
    ``Delta = 0`` we use the positive-energy normal-state limit.
    """
    E = np.asarray(energy_values_J, dtype=float)
    rho = np.asarray(rho_E, dtype=float)
    omega = np.asarray(omega_values_J, dtype=float)

    if E.ndim != 1 or rho.shape != E.shape:
        raise ValueError("energy_values_J and rho_E must be one-dimensional arrays with the same shape.")

    f_E = fermi_positive_energy(E, Te_K)

    out = np.zeros_like(omega, dtype=float)
    E_max = float(E[-1])
    E_floor = max(float(E[1]) * 0.5, 1.0e-300)

    lower = max(delta_J, E_floor) if delta_J > 0.0 else E_floor

    for i, Om in enumerate(omega):
        Ep = E + Om
        mask = (E >= lower) & (Ep <= E_max)

        if delta_J > 0.0:
            mask &= Ep >= delta_J

        if np.count_nonzero(mask) < 2:
            out[i] = 0.0
            continue

        E_m = E[mask]
        Ep_m = Ep[mask]

        rho_E_m = rho[mask]
        rho_Ep = np.interp(Ep_m, E, rho, left=0.0, right=0.0)

        denom = np.maximum(E_m * Ep_m, E_floor * E_floor)
        coherence = 1.0 - delta_J * delta_J / denom

        f_Ep = fermi_positive_energy(Ep_m, Te_K)

        integrand = rho_E_m * rho_Ep * coherence * (f_E[mask] - f_Ep)
        integrand = np.nan_to_num(integrand, nan=0.0, posinf=0.0, neginf=0.0)

        out[i] = float(np.trapz(integrand, E_m))

    return np.maximum(out, 0.0)


def recombination_phase_space_spectrum(
    energy_values_J: np.ndarray,
    rho_E: np.ndarray,
    omega_values_J: np.ndarray,
    *,
    Te_K: float,
    delta_J: float,
) -> np.ndarray:
    """
    Compute J_R(Omega) for one (Te, Delta, q) slice.

    Implemented expression:

        J_R(Omega) =
            int_Delta^{Omega-Delta} dE rho(E) rho(Omega-E)
            [1 + Delta^2/(E(Omega-E))]
            f(E) f(Omega-E) [exp(Omega/kBT)-1].

    For ``Delta = 0`` the superconducting recombination channel is set to zero.
    """
    E = np.asarray(energy_values_J, dtype=float)
    rho = np.asarray(rho_E, dtype=float)
    omega = np.asarray(omega_values_J, dtype=float)

    if E.ndim != 1 or rho.shape != E.shape:
        raise ValueError("energy_values_J and rho_E must be one-dimensional arrays with the same shape.")

    out = np.zeros_like(omega, dtype=float)

    if delta_J <= 0.0:
        return out

    E_max = float(E[-1])
    E_floor = max(float(E[1]) * 0.5, 1.0e-300)

    for i, Om in enumerate(omega):
        if Om <= 2.0 * delta_J:
            out[i] = 0.0
            continue

        Ep = Om - E

        mask = (
            (E >= delta_J)
            & (E <= Om - delta_J)
            & (Ep >= delta_J)
            & (Ep <= E_max)
        )

        if np.count_nonzero(mask) < 2:
            out[i] = 0.0
            continue

        E_m = E[mask]
        Ep_m = Ep[mask]

        rho_E_m = rho[mask]
        rho_Ep = np.interp(Ep_m, E, rho, left=0.0, right=0.0)

        denom = np.maximum(E_m * Ep_m, E_floor * E_floor)
        coherence = 1.0 + delta_J * delta_J / denom

        thermal_factor = pair_recombination_thermal_factor(E_m, Ep_m, Te_K)

        integrand = rho_E_m * rho_Ep * coherence * thermal_factor
        integrand = np.nan_to_num(integrand, nan=0.0, posinf=0.0, neginf=0.0)

        out[i] = float(np.trapz(integrand, E_m))

    return np.maximum(out, 0.0)


def fermi_positive_energy(energy_J: np.ndarray, T_K: float) -> np.ndarray:
    """
    Fermi function for positive quasiparticle energies.
    """
    if T_K <= 0.0:
        raise ValueError("T_K must be positive.")

    x = np.asarray(energy_J, dtype=float) / (K_B_J_K * T_K)
    x = np.clip(x, 0.0, 700.0)

    return 1.0 / (np.exp(x) + 1.0)


def pair_recombination_thermal_factor(
    E_J: np.ndarray,
    Ep_J: np.ndarray,
    T_K: float,
) -> np.ndarray:
    """
    Stable version of

        f(E) f(E') [exp((E+E')/kBT) - 1].

    Direct evaluation can overflow at low temperature. With
    ``a=E/kBT`` and ``b=E'/kBT``,

        (exp(a+b)-1)/[(exp(a)+1)(exp(b)+1)]
        =
        [1-exp(-(a+b))]/[(1+exp(-a))(1+exp(-b))].
    """
    if T_K <= 0.0:
        raise ValueError("T_K must be positive.")

    a = np.asarray(E_J, dtype=float) / (K_B_J_K * T_K)
    b = np.asarray(Ep_J, dtype=float) / (K_B_J_K * T_K)

    exp_minus_a = np.exp(-np.minimum(a, 700.0))
    exp_minus_b = np.exp(-np.minimum(b, 700.0))
    exp_minus_sum = np.exp(-np.minimum(a + b, 700.0))

    numerator = 1.0 - exp_minus_sum
    denominator = (1.0 + exp_minus_a) * (1.0 + exp_minus_b)

    return numerator / denominator


def phase_space_summary(catalog: PhaseSpaceCatalog) -> dict[str, Any]:
    """
    Return a compact summary dictionary.
    """
    JS = catalog.J_S_TdqO_J
    JR = catalog.J_R_TdqO_J

    return {
        "backend": str(catalog.metadata.get("backend", "unknown")),
        "shape": list(catalog.shape),
        "n_Te": int(catalog.Te_values_K.size),
        "n_delta": int(catalog.delta_values_J.size),
        "n_q": int(catalog.q_values_m_inv.size),
        "n_omega": int(catalog.omega_values_J.size),
        "Te_min_K": float(np.min(catalog.Te_values_K)),
        "Te_max_K": float(np.max(catalog.Te_values_K)),
        "omega_max_J": float(np.max(catalog.omega_values_J)),
        "omega_max_meV": J_to_meV(float(np.max(catalog.omega_values_J))),
        "delta_min_meV": J_to_meV(float(np.min(catalog.delta_values_J))),
        "delta_max_meV": J_to_meV(float(np.max(catalog.delta_values_J))),
        "gamma_min_meV": J_to_meV(float(np.min(catalog.gamma_values_J))),
        "gamma_max_meV": J_to_meV(float(np.max(catalog.gamma_values_J))),
        "q_min_m_inv": float(np.min(catalog.q_values_m_inv)),
        "q_max_m_inv": float(np.max(catalog.q_values_m_inv)),
        "J_S_min_J": float(np.min(JS)),
        "J_S_max_J": float(np.max(JS)),
        "J_S_is_finite": bool(np.all(np.isfinite(JS))),
        "J_R_min_J": float(np.min(JR)),
        "J_R_max_J": float(np.max(JR)),
        "J_R_is_finite": bool(np.all(np.isfinite(JR))),
        "grid_is_downsampled": bool(catalog.metadata.get("grid_is_downsampled", False)),
    }


def save_phase_space_catalog_npz(
    catalog: PhaseSpaceCatalog,
    path: str | Path,
) -> Path:
    """
    Save phase-space catalogue to compressed ``.npz``.
    """
    output = Path(path)
    output.parent.mkdir(parents=True, exist_ok=True)

    np.savez_compressed(
        output,
        Te_values_K=catalog.Te_values_K,
        omega_values_J=catalog.omega_values_J,
        delta_values_J=catalog.delta_values_J,
        gamma_values_J=catalog.gamma_values_J,
        q_values_m_inv=catalog.q_values_m_inv,
        J_S_TdqO_J=catalog.J_S_TdqO_J,
        J_R_TdqO_J=catalog.J_R_TdqO_J,
        delta_indices=catalog.delta_indices,
        q_indices=catalog.q_indices,
        metadata=np.array(catalog.metadata, dtype=object),
    )

    return output


def load_phase_space_catalog_npz(path: str | Path) -> PhaseSpaceCatalog:
    """
    Load a phase-space catalogue saved by :func:`save_phase_space_catalog_npz`.
    """
    source = Path(path)

    with np.load(source, allow_pickle=True) as data:
        metadata = data["metadata"].item()

        return PhaseSpaceCatalog(
            Te_values_K=np.asarray(data["Te_values_K"], dtype=float),
            omega_values_J=np.asarray(data["omega_values_J"], dtype=float),
            delta_values_J=np.asarray(data["delta_values_J"], dtype=float),
            gamma_values_J=np.asarray(data["gamma_values_J"], dtype=float),
            q_values_m_inv=np.asarray(data["q_values_m_inv"], dtype=float),
            J_S_TdqO_J=np.asarray(data["J_S_TdqO_J"], dtype=float),
            J_R_TdqO_J=np.asarray(data["J_R_TdqO_J"], dtype=float),
            delta_indices=np.asarray(data["delta_indices"], dtype=np.int64),
            q_indices=np.asarray(data["q_indices"], dtype=np.int64),
            metadata=dict(metadata),
        )


def _select_axis_indices(n_total: int, n_requested: int) -> np.ndarray:
    """
    Select approximately evenly spaced indices from an axis.
    """
    if n_total <= 0:
        raise ValueError("n_total must be positive.")

    if n_requested <= 0:
        raise ValueError("n_requested must be positive.")

    if n_requested >= n_total:
        return np.arange(n_total, dtype=np.int64)

    return np.unique(
        np.round(np.linspace(0, n_total - 1, n_requested)).astype(np.int64)
    )