"""Explicit, traceable shape preprocessing; this module cannot admit SI materials.

The input ordinates and axis retain their original, uncertified units. Nothing
here chooses an atom/cell basis, multiplies the DOS to enforce a mode count, or
changes the production loader. Every modification is recorded in the manifest.
"""
from __future__ import annotations

from dataclasses import dataclass
import hashlib
import io
from pathlib import Path
import numpy as np


POLICY = "stable-first-common-support-origin-v1"


@dataclass(frozen=True)
class DerivedPhononShape:
    axis: np.ndarray
    alpha2F: np.ndarray
    phdos: np.ndarray
    source_rows: np.ndarray
    manifest: dict

    def evaluate(self, coordinates):
        """Linear interpolants and their ratio, with zero coupling at zero DOS.

        Interpolate alpha2F and F on the same support, then divide. Because both
        vanish at a support boundary, the ratio has a finite one-sided limit.
        This differs intentionally from interpolating their nodal ratio first.
        No extrapolation is allowed; even shape-only queries must stay in range.
        """
        q = np.asarray(coordinates, dtype=float)
        if np.any(~np.isfinite(q)) or np.any(q < self.axis[0]) or np.any(q > self.axis[-1]):
            raise ValueError("query outside the finite source support")
        alpha = np.interp(q, self.axis, self.alpha2F)
        dos = np.interp(q, self.axis, self.phdos)
        try:
            with np.errstate(over="raise", invalid="raise", divide="raise"):
                ratio = np.divide(alpha, dos, out=np.zeros_like(alpha), where=dos > 0)
        except FloatingPointError as exc:
            raise ValueError("the interpolated coupling ratio is not representable as a finite float") from exc
        return alpha, dos, ratio


def derive_phonon_shape(path, *, source_url: str, source_revision: str,
                        dos_floor: float = 0.0) -> DerivedPhononShape:
    """Stable-sort, keep first duplicates, and explicitly regularize support.

    ``dos_floor`` is in the *original ordinate units*. It is never inferred from
    a desired normalization. The default removes only nonpositive DOS. Setting
    alpha(0)=0 is an explicit integrability condition for integral alpha/axis;
    it is not evidence that the raw values were correct at the origin.
    """
    if not isinstance(source_url, str) or not source_url.startswith("https://"):
        raise ValueError("an explicit HTTPS primary-source URL is required")
    if not isinstance(source_revision, str) or not source_revision.strip():
        raise ValueError("an explicit source revision is required")
    if not np.isfinite(dos_floor) or dos_floor < 0:
        raise ValueError("dos_floor must be finite and nonnegative")
    raw = Path(path).read_bytes()
    data = np.loadtxt(io.BytesIO(raw), comments="#", ndmin=2)
    if data.shape[1] != 3 or data.shape[0] < 2 or not np.all(np.isfinite(data)):
        raise ValueError("expected at least two finite rows: axis, alpha2F, phdos")
    if np.any(data[:, 0] < 0) or np.any(data[:, 1] < 0):
        raise ValueError("negative frequencies or negative alpha2F are not covered by this policy")
    order = np.argsort(data[:, 0], kind="stable")
    ordered = data[order]
    first = np.r_[True, np.diff(ordered[:, 0]) != 0]
    kept = order[first]
    dropped = order[~first]
    unique = data[kept].copy()
    if len(unique) < 2:
        raise ValueError("at least two distinct frequencies are required")
    mask = unique[:, 2] <= dos_floor
    changed_dos = mask & (unique[:, 2] != 0)
    changed_alpha = mask & (unique[:, 1] != 0)
    origin = (unique[:, 0] == 0) & (unique[:, 1] != 0)
    unique[mask, 2] = 0.0
    unique[mask | origin, 1] = 0.0
    arrays = [unique[:, i].copy() for i in range(3)]
    kept = kept.copy()
    for array in [*arrays, kept]:
        array.flags.writeable = False
    manifest = {
        "schema": "pysnspd.derived-phonon-shape.v1", "policy": POLICY,
        "status": "DERIVED_SHAPE_ONLY_NOT_SI_ADMISSION",
        "source_sha256": hashlib.sha256(raw).hexdigest(),
        "source_url": source_url, "source_revision": source_revision,
        "input_rows": int(len(data)), "output_rows": int(len(unique)),
        "operations": [
            {"name": "stable_sort_ascending_axis", "moved_rows": int(np.sum(order != np.arange(len(order))))},
            {"name": "keep_first_duplicate", "dropped_source_rows_zero_based": dropped.tolist()},
            {"name": "zero_dos_at_or_below_explicit_floor", "floor_original_ordinate_units": float(dos_floor),
             "changed_source_rows_zero_based": kept[changed_dos].tolist()},
            {"name": "zero_alpha2F_where_DOS_is_zero", "changed_source_rows_zero_based": kept[changed_alpha].tolist()},
            {"name": "set_alpha2F_at_zero_frequency_to_zero", "changed_source_rows_zero_based": kept[origin].tolist()},
        ],
        "interpolation": "piecewise linear alpha2F and F; quotient only where interpolated F>0; zero at F=0",
        "renormalization_factor": 1.0,
        "certified_axis_unit": None, "certified_DOS_unit": None,
        "certified_normalization_basis": None, "certified_density_m3": None,
        "limitation": "Preprocessing produces a numerical shape, not a validated physical material or a transient prediction.",
    }
    return DerivedPhononShape(*arrays, kept, manifest)
