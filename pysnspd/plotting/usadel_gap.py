"""Extra plotting utilities for Usadel equilibrium-gap diagnostics.

This module is intentionally light on pySNSPD internals: it reads the
pre-run Usadel ``.npz`` catalog directly, extracts ``Delta_eq(T, q)``, and
builds a publication-ready PDF plot.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Mapping, Sequence

import numpy as np


_E_CHARGE_C = 1.602176634e-19


@dataclass(frozen=True)
class UsadelGapCatalog:
    """Minimal representation needed for E1 pre-run gap plots.

    Attributes
    ----------
    temperature_K:
        One-dimensional electron-temperature grid in K.
    q_m_inv:
        One-dimensional superfluid momentum grid in m^-1.
    gap_meV:
        Equilibrium gap array with shape ``(n_temperature, n_q)`` in meV.
    source_key:
        Name of the field used from the ``.npz`` catalog.
    """

    temperature_K: np.ndarray
    q_m_inv: np.ndarray
    gap_meV: np.ndarray
    source_key: str


def _npz_mapping(npz: np.lib.npyio.NpzFile | Mapping[str, np.ndarray]) -> Mapping[str, np.ndarray]:
    return npz  # type: ignore[return-value]


def _first_existing_array(
    data: Mapping[str, np.ndarray],
    candidates: Sequence[str],
) -> tuple[str, np.ndarray] | None:
    for key in candidates:
        if key in data:
            return key, np.asarray(data[key])
    return None


def _one_dimensional(name: str, array: np.ndarray) -> np.ndarray:
    out = np.asarray(array, dtype=float).squeeze()
    if out.ndim != 1:
        raise ValueError(f"Expected {name} to be one-dimensional, got shape {array.shape}.")
    return out


def _find_temperature_grid(data: Mapping[str, np.ndarray]) -> tuple[str, np.ndarray]:
    explicit = _first_existing_array(
        data,
        (
            "T_values_K",
            "temperature_values_K",
            "Te_values_K",
            "te_values_K",
            "T_grid_K",
            "temperature_grid_K",
            "temperatures_K",
        ),
    )
    if explicit is not None:
        key, arr = explicit
        return key, _one_dimensional(key, arr)

    for key in data:
        key_l = key.lower()
        if ("temp" in key_l or key_l in {"t", "te"}) and key_l.endswith("_k"):
            arr = np.asarray(data[key])
            if arr.squeeze().ndim == 1:
                return key, _one_dimensional(key, arr)

    raise KeyError(
        "Could not find a temperature grid in the Usadel catalog. Expected a key like "
        "T_values_K, temperature_values_K or Te_values_K."
    )


def _find_q_grid(data: Mapping[str, np.ndarray]) -> tuple[str, np.ndarray]:
    explicit = _first_existing_array(
        data,
        (
            "q_values_m_inv",
            "q_grid_m_inv",
            "q_values_inv_m",
            "q_m_inv",
            "q_values",
        ),
    )
    if explicit is not None:
        key, arr = explicit
        return key, _one_dimensional(key, arr)

    for key in data:
        key_l = key.lower()
        if key_l.startswith("q") and ("m_inv" in key_l or "inv_m" in key_l or key_l == "q"):
            arr = np.asarray(data[key])
            if arr.squeeze().ndim == 1:
                return key, _one_dimensional(key, arr)

    raise KeyError(
        "Could not find a q grid in the Usadel catalog. Expected a key like "
        "q_values_m_inv or q_grid_m_inv."
    )


def _candidate_gap_keys(data: Mapping[str, np.ndarray]) -> list[str]:
    preferred = [
        "delta_eq_J",
        "Delta_eq_J",
        "gap_eq_J",
        "Delta_equilibrium_J",
        "delta_equilibrium_J",
        "delta_eq_meV",
        "Delta_eq_meV",
        "gap_eq_meV",
        "delta_grid_J",
        "Delta_grid_J",
        "gap_grid_J",
    ]
    out: list[str] = [key for key in preferred if key in data]

    banned = ("rho", "dos", "energy", "omega", "anomalous", "gamma", "phase")
    for key in data:
        if key in out:
            continue
        key_l = key.lower()
        if any(word in key_l for word in banned):
            continue
        if ("delta" in key_l or "gap" in key_l) and ("eq" in key_l or "equil" in key_l or "grid" in key_l):
            out.append(key)
    return out


def _orient_gap_array(
    key: str,
    array: np.ndarray,
    n_temperature: int,
    n_q: int,
) -> np.ndarray:
    gap = np.asarray(array, dtype=float).squeeze()

    if gap.ndim != 2:
        raise ValueError(
            f"Candidate gap field {key!r} has shape {array.shape}; expected a 2-D "
            "array compatible with (n_temperature, n_q)."
        )

    if gap.shape == (n_temperature, n_q):
        return gap
    if gap.shape == (n_q, n_temperature):
        return gap.T

    raise ValueError(
        f"Candidate gap field {key!r} has shape {gap.shape}; expected either "
        f"({n_temperature}, {n_q}) or ({n_q}, {n_temperature})."
    )


def _gap_to_mev(key: str, gap: np.ndarray) -> np.ndarray:
    key_l = key.lower()
    out = np.asarray(gap, dtype=float)

    if "mev" in key_l:
        return out
    if key_l.endswith("_j") or "_j_" in key_l:
        return out / _E_CHARGE_C * 1.0e3

    finite = np.abs(out[np.isfinite(out)])
    if finite.size == 0:
        return out

    max_abs = float(np.nanmax(finite))
    # NbN gaps stored in joules are around 1e-22 J.  This conservative
    # fallback keeps old catalogs usable even if the key did not include _J.
    if max_abs < 1.0e-18:
        return out / _E_CHARGE_C * 1.0e3
    return out


def load_usadel_gap_catalog(npz_path: str | Path) -> UsadelGapCatalog:
    """Load ``Delta_eq(T, q)`` from a pySNSPD pre-run Usadel catalog.

    The loader accepts several historical field names so that the extra plot
    pipeline is not tied to one exact catalog version.
    """

    path = Path(npz_path)
    if not path.exists():
        raise FileNotFoundError(f"Usadel catalog not found: {path}")

    with np.load(path, allow_pickle=False) as npz:
        data = _npz_mapping(npz)
        _, temperature_K = _find_temperature_grid(data)
        _, q_m_inv = _find_q_grid(data)

        errors: list[str] = []
        for key in _candidate_gap_keys(data):
            try:
                gap = _orient_gap_array(key, np.asarray(data[key]), temperature_K.size, q_m_inv.size)
            except ValueError as exc:
                errors.append(str(exc))
                continue
            gap_meV = _gap_to_mev(key, gap)
            return UsadelGapCatalog(
                temperature_K=np.asarray(temperature_K, dtype=float),
                q_m_inv=np.asarray(q_m_inv, dtype=float),
                gap_meV=np.asarray(gap_meV, dtype=float),
                source_key=key,
            )

    detail = "\n".join(f"  - {err}" for err in errors)
    raise KeyError(
        "Could not find a usable Delta_eq(T, q) field in the Usadel catalog. "
        "Expected a 2-D array with a name like delta_eq_J, Delta_eq_J or gap_eq_J."
        + (f"\nRejected candidates:\n{detail}" if detail else "")
    )


def estimate_q_critical(
    q_m_inv: np.ndarray,
    gap_meV: np.ndarray,
    *,
    gap_fraction_threshold: float = 1.0e-3,
) -> float:
    """Estimate q_c as the largest q with a non-negligible low-T gap."""

    q = np.asarray(q_m_inv, dtype=float)
    gap = np.asarray(gap_meV, dtype=float)
    if gap.ndim != 2:
        raise ValueError(f"gap_meV must be 2-D, got shape {gap.shape}.")
    if q.ndim != 1 or q.size != gap.shape[1]:
        raise ValueError("q_m_inv must be 1-D and match the second gap dimension.")

    reference = gap[0, :]
    finite = np.isfinite(reference)
    if not np.any(finite):
        return float(np.nanmax(q))

    max_gap = float(np.nanmax(np.abs(reference[finite])))
    if max_gap <= 0.0:
        return float(np.nanmax(q))

    active = finite & (np.abs(reference) >= gap_fraction_threshold * max_gap)
    if not np.any(active):
        return float(np.nanmax(q))
    return float(np.nanmax(q[active]))


def interpolate_gap_curves(
    catalog: UsadelGapCatalog,
    *,
    n_curves: int = 6,
    q_critical_m_inv: float | None = None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return target q values and interpolated ``Delta_eq(T)`` curves.

    Returns
    -------
    q_targets_m_inv:
        Shape ``(n_curves,)``.
    temperature_K:
        Sorted temperature grid.
    curves_meV:
        Shape ``(n_curves, n_temperature)``.
    """

    if n_curves < 2:
        raise ValueError("n_curves must be at least 2.")

    temperature = np.asarray(catalog.temperature_K, dtype=float)
    q = np.asarray(catalog.q_m_inv, dtype=float)
    gap = np.asarray(catalog.gap_meV, dtype=float)

    t_order = np.argsort(temperature)
    q_order = np.argsort(q)
    temperature = temperature[t_order]
    q = q[q_order]
    gap = gap[np.ix_(t_order, q_order)]

    qcrit = float(q_critical_m_inv) if q_critical_m_inv is not None else estimate_q_critical(q, gap)
    qcrit = min(qcrit, float(np.nanmax(q)))
    if qcrit <= float(np.nanmin(q)):
        raise ValueError(f"Invalid q critical value: {qcrit:.6e} m^-1.")

    q_targets = np.linspace(float(np.nanmin(q)), qcrit, int(n_curves))
    curves = np.empty((q_targets.size, temperature.size), dtype=float)
    for i_t in range(temperature.size):
        row = gap[i_t, :]
        valid = np.isfinite(row)
        if np.count_nonzero(valid) < 2:
            curves[:, i_t] = np.nan
            continue
        curves[:, i_t] = np.interp(q_targets, q[valid], row[valid], left=np.nan, right=np.nan)

    return q_targets, temperature, curves


def plot_gap_eq_vs_temperature(
    catalog: UsadelGapCatalog,
    output_path: str | Path,
    *,
    n_curves: int = 6,
    q_critical_m_inv: float | None = None,
    dpi: int = 480,
    title: str | None = None,
) -> Path:
    """Plot ``Delta_eq`` versus temperature for multiple q values.

    The figure is saved as PDF when ``output_path`` has suffix ``.pdf``.  A
    non-PDF suffix is allowed, but the E1 pipeline uses PDF by default.
    """

    import matplotlib

    matplotlib.use("Agg", force=True)
    import matplotlib.pyplot as plt

    output = Path(output_path)
    output.parent.mkdir(parents=True, exist_ok=True)

    q_targets, temperature, curves = interpolate_gap_curves(
        catalog,
        n_curves=n_curves,
        q_critical_m_inv=q_critical_m_inv,
    )
    qcrit = float(q_targets[-1])

    fig, ax = plt.subplots(figsize=(6.4, 4.2), constrained_layout=True)
    for q_value, curve in zip(q_targets, curves):
        ratio = q_value / qcrit if qcrit != 0.0 else np.nan
        ax.plot(
            temperature,
            curve,
            linewidth=1.7,
            label=rf"$q/q_c={ratio:.2f}$",
        )

    ax.set_xlabel(r"Electron temperature $T_e$ (K)")
    ax.set_ylabel(r"Equilibrium gap $\Delta_{\mathrm{eq}}$ (meV)")
    ax.set_xlim(float(np.nanmin(temperature)), float(np.nanmax(temperature)))
    ax.set_ylim(bottom=0.0)
    ax.grid(True, alpha=0.25)
    ax.legend(frameon=False, fontsize=8, ncol=2)
    if title:
        ax.set_title(title)

    fig.savefig(output, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    return output


__all__ = [
    "UsadelGapCatalog",
    "estimate_q_critical",
    "interpolate_gap_curves",
    "load_usadel_gap_catalog",
    "plot_gap_eq_vs_temperature",
]
