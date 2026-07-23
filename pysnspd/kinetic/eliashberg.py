"""Eliashberg spectral data and phonon-density utilities.

The external ``*.dat`` files used here are material data from

    Simon et al.,
    "Ab initio modeling of nonequilibrium dynamics in superconducting
    detectors and qubits", Physical Review B 112, 174512 (2025).

For NbN, the MIT/Simon data file used in the OE5 debug runs has header

    #E (THz)    a^2F    PhDOS (st/THz)

The file is intentionally not stored in the git repository. It should live
under ``project.big_data_root/catalogs/simon_2025/``. This keeps the source
code light and makes the provenance of external material data explicit.

Conventions
-----------
The first column is a frequency axis in THz. The thermal projection in the
pySNSPD appendix writes Omega as a phonon energy. Therefore this module
converts the THz axis to an energy axis using

    Omega_J = h * f_THz * 1e12.

The column ``a^2F`` has no explicit unit in the file header. We treat it as
the tabulated Eliashberg spectral weight evaluated on the THz axis. This
preserves the dimensionless coupling check

    lambda = 2 int df a^2F(f) / f,

and is consistent with evaluating the same spectral function on the energy
axis in the projected-power integrals.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import numpy as np


PLANCK_J_S = 6.62607015e-34
THZ_HZ = 1.0e12
EV_J = 1.602176634e-19
MEV_J = 1.602176634e-22


def thz_to_j(f_THz: np.ndarray | float) -> np.ndarray | float:
    """Convert ordinary frequency in THz to energy in joule."""
    return PLANCK_J_S * THZ_HZ * np.asarray(f_THz)


def j_to_mev(E_J: np.ndarray | float) -> np.ndarray | float:
    """Convert energy in joule to meV."""
    return np.asarray(E_J) / MEV_J


@dataclass(frozen=True)
class EliashbergSpectrum:
    """Tabulated Eliashberg function and phonon DOS.

    Parameters
    ----------
    frequency_THz:
        Frequency axis from the source file.
    omega_J:
        Energy axis obtained from ``omega = h f``.
    alpha2F:
        Tabulated Eliashberg spectral function. The source header does not
        assign an explicit unit to this column.
    phdos_states_per_THz:
        Phonon density of states in states/THz.
    metadata:
        Source, units, normalization checks and clipping information.
    """

    frequency_THz: np.ndarray
    omega_J: np.ndarray
    alpha2F: np.ndarray
    phdos_states_per_THz: np.ndarray
    metadata: dict[str, Any] = field(default_factory=dict)

    @property
    def omega_meV(self) -> np.ndarray:
        return np.asarray(j_to_mev(self.omega_J), dtype=float)


    def alpha2F_on_omega_J(self, omega_J: np.ndarray) -> np.ndarray:
        """Interpolate alpha^2F onto an energy axis in joule."""
        return np.interp(
            np.asarray(omega_J, dtype=float),
            self.omega_J,
            self.alpha2F,
            left=0.0,
            right=0.0,
        )

    def phdos_on_omega_J(self, omega_J: np.ndarray) -> np.ndarray:
        """Interpolate PhDOS onto an energy axis in joule."""
        return np.interp(
            np.asarray(omega_J, dtype=float),
            self.omega_J,
            self.phdos_states_per_THz,
            left=0.0,
            right=0.0,
        )


def load_simon_eliashberg_dat(
    path: str | Path,
    *,
    clip_negative: bool = True,
) -> EliashbergSpectrum:
    """Load a Simon et al. 2025 ``a^2F``/PhDOS material file.

    The loader expects the three-column format used by the MIT data files:

        #E (THz)    a^2F    PhDOS (st/THz)

    Negative values can appear as small interpolation artifacts in the phonon
    DOS tail. By default, negative ``a^2F`` and PhDOS values are clipped to
    zero and the number of clipped samples is recorded in metadata.
    """
    data_path = Path(path).expanduser().resolve()
    if not data_path.exists():
        raise FileNotFoundError(f"Eliashberg data file does not exist: {data_path}")

    header = ""
    with data_path.open("r", encoding="utf-8", errors="replace") as f:
        first = f.readline().strip()
        if first.startswith("#"):
            header = first

    data = np.loadtxt(data_path, comments="#")
    if data.ndim != 2 or data.shape[1] < 3:
        raise ValueError(
            "Expected at least three columns: frequency_THz, alpha2F, PhDOS."
        )

    frequency_THz = np.asarray(data[:, 0], dtype=float)
    alpha2F = np.asarray(data[:, 1], dtype=float)
    phdos = np.asarray(data[:, 2], dtype=float)

    order = np.argsort(frequency_THz)
    frequency_THz = frequency_THz[order]
    alpha2F = alpha2F[order]
    phdos = phdos[order]

    if np.any(~np.isfinite(frequency_THz)):
        raise ValueError("Frequency axis contains non-finite values.")
    if np.any(np.diff(frequency_THz) < 0.0):
        raise ValueError("Frequency axis could not be sorted monotonically.")
    if float(frequency_THz[-1]) <= 0.0:
        raise ValueError("Frequency axis must extend above zero.")

    n_alpha_negative = int(np.sum(alpha2F < 0.0))
    n_phdos_negative = int(np.sum(phdos < 0.0))

    if clip_negative:
        alpha2F = np.maximum(alpha2F, 0.0)
        phdos = np.maximum(phdos, 0.0)

    positive = frequency_THz > 0.0
    lambda_ep = 2.0 * np.trapezoid(
        alpha2F[positive] / frequency_THz[positive],
        frequency_THz[positive],
    )

    omega_J = np.asarray(thz_to_j(frequency_THz), dtype=float)
    omega_meV = np.asarray(j_to_mev(omega_J), dtype=float)

    int_alpha = float(np.trapezoid(alpha2F, omega_meV))
    int_omega_alpha = float(np.trapezoid(omega_meV * alpha2F, omega_meV))
    int_phdos = float(np.trapezoid(phdos, frequency_THz))

    metadata = {
        "source": "Simon et al. 2025, Physical Review B 112, 174512",
        "source_note": (
            "External material data. Store under big_data/catalogs/simon_2025; "
            "do not commit raw data to the pySNSPD repository."
        ),
        "path": str(data_path),
        "header": header,
        "frequency_unit": "THz",
        "omega_energy_unit": "J",
        "omega_plot_unit": "meV",
        "alpha2F_unit_from_header": "not explicitly specified",
        "alpha2F_policy": (
            "Treated as the tabulated Eliashberg spectral weight on the THz "
            "axis and evaluated on the energy axis after Omega=h f conversion."
        ),
        "phdos_unit_from_header": "states/THz",
        "clip_negative": bool(clip_negative),
        "n_alpha2F_negative_clipped": n_alpha_negative if clip_negative else 0,
        "n_phdos_negative_clipped": n_phdos_negative if clip_negative else 0,
        "n_points": int(frequency_THz.size),
        "frequency_min_THz": float(frequency_THz[0]),
        "frequency_max_THz": float(frequency_THz[-1]),
        "omega_max_meV": float(omega_meV[-1]),
        "alpha2F_max": float(np.max(alpha2F)),
        "alpha2F_peak_meV": float(omega_meV[int(np.argmax(alpha2F))]),
        "phdos_max_states_per_THz": float(np.max(phdos)),
        "phdos_peak_meV": float(omega_meV[int(np.argmax(phdos))]),
        "lambda_ep": float(lambda_ep),
        "integral_alpha2F_dmeV": int_alpha,
        "integral_omega_alpha2F_dmeV2": int_omega_alpha,
        "integral_phdos_dTHz": int_phdos,
    }

    return EliashbergSpectrum(
        frequency_THz=frequency_THz,
        omega_J=omega_J,
        alpha2F=alpha2F,
        phdos_states_per_THz=phdos,
        metadata=metadata,
    )
