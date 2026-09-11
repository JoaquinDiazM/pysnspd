"""Explicit, non-repairing admission of experimental phonon material tables.

This module is independent of the production Eliashberg loader. A successful
audit verifies the supplied contract, not the truth of a bibliographic claim.
The caller must establish the units, normalization basis and provenance before
marking their evidence as verified. Synthetic data remain labelled synthetic.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
import hashlib
import io
from pathlib import Path
from typing import Any

import numpy as np


PLANCK_J_S = 6.62607015e-34
EV_J = 1.602176634e-19
KB_J_K = 1.380649e-23
ENERGY_PER_AXIS_UNIT = {
    "THz": PLANCK_J_S * 1e12,
    "rad/s": PLANCK_J_S / (2 * np.pi),
    "meV": EV_J * 1e-3,
    "eV": EV_J,
    "J": 1.0,
}
ENERGY_PER_DOS_UNIT = {
    f"states/{unit}": value
    for unit, value in ENERGY_PER_AXIS_UNIT.items()
    if unit != "rad/s"
}
ENERGY_PER_DOS_UNIT["states/(rad/s)"] = ENERGY_PER_AXIS_UNIT["rad/s"]


@dataclass(frozen=True)
class MaterialSpec:
    """Declared contract for columns (axis, alpha2F, phonon DOS).

    ``basis_density_m3`` counts the same objects as ``normalization_basis``:
    atoms for ``atom``, crystallographic cells for ``cell``. A cell may contain
    several atoms. A complete DOS must integrate to three modes per atom.
    ``alpha2F`` uses the dimensionless Eliashberg convention in B.10/B.17;
    unlike a DOS, its ordinate acquires no Jacobian when its argument changes.
    Evidence strings identify sources/sections, not merely a plausible sum.
    """

    material_id: str
    axis_unit: str
    phdos_unit: str
    normalization_basis: str
    atoms_per_basis: int
    basis_density_m3: float
    source_url: str
    source_revision: str
    raw_sha256: str
    unit_evidence: str
    normalization_evidence: str
    provenance_evidence: str
    units_verified: bool = False
    normalization_verified: bool = False
    provenance_verified: bool = False
    synthetic: bool = False
    alpha2F_convention: str = "dimensionless_eliashberg"


@dataclass(frozen=True)
class AdmissionPolicy:
    """A declared quadrature tolerance, never an automatic renormalization."""

    relative_mode_tolerance: float = 5e-3

    def __post_init__(self) -> None:
        if not np.isfinite(self.relative_mode_tolerance) or not (
            0 <= self.relative_mode_tolerance < 1
        ):
            raise ValueError("relative_mode_tolerance must be finite in [0, 1).")


@dataclass(frozen=True)
class MaterialAudit:
    status: str
    source_sha256: str
    metadata: dict[str, Any]
    policy: dict[str, Any]
    issues: tuple[dict[str, str], ...]
    metrics_assuming_declared_units: dict[str, Any]
    header: str | None

    @property
    def admitted(self) -> bool:
        return self.status.startswith("ADMITTED_")

    def to_dict(self) -> dict[str, Any]:
        result = asdict(self)
        result["repairs_applied"] = []
        result["scope"] = "input contract only; no transient or absolute-rate validation"
        return result


class MaterialAdmissionError(ValueError):
    def __init__(self, audit: MaterialAudit):
        self.audit = audit
        super().__init__("Material rejected: " + ", ".join(i["code"] for i in audit.issues))


@dataclass(frozen=True)
class AdmittedMaterial:
    """SI spectral data available only after successful strict admission."""

    omega_J: np.ndarray
    alpha2F: np.ndarray
    phdos_per_basis_J: np.ndarray
    g_ph_per_m3_J: np.ndarray
    audit: MaterialAudit

    def thermal_energy_density(self, temperature_K: float) -> float:
        """Return phonon energy above its zero-point reference, in J/m^3.

        The finite table is integrated linearly by the trapezoid rule. At zero
        energy, Omega*n_B has its continuous limit k_B*T; there is no 0/0.
        """
        if not np.isfinite(temperature_K) or temperature_K < 0:
            raise ValueError("Temperature must be finite and nonnegative.")
        if temperature_K == 0:
            return 0.0
        reduced = self.omega_J / (KB_J_K * temperature_K)
        factor = np.zeros_like(reduced)
        small = reduced < 1e-5
        factor[small] = 1 - reduced[small] / 2 + reduced[small] ** 2 / 12
        regular = (~small) & (reduced < 700)
        factor[regular] = reduced[regular] / np.expm1(reduced[regular])
        return float(np.trapezoid(
            self.g_ph_per_m3_J * KB_J_K * temperature_K * factor, self.omega_J
        ))


def _inspect(
    raw: bytes, spec: MaterialSpec, policy: AdmissionPolicy
) -> tuple[MaterialAudit, np.ndarray | None]:
    issues: list[dict[str, str]] = []

    def reject(code: str, detail: str) -> None:
        issues.append({"code": code, "detail": detail})

    def finish() -> MaterialAudit:
        return MaterialAudit(
            status="REJECTED" if issues else ("ADMITTED_SYNTHETIC" if spec.synthetic else "ADMITTED_MATERIAL"),
            source_sha256=actual_hash, metadata=asdict(spec), policy=asdict(policy),
            issues=tuple(issues), metrics_assuming_declared_units=metrics, header=header,
        )

    actual_hash = hashlib.sha256(raw).hexdigest()
    metrics: dict[str, Any] = {}
    first = raw.splitlines()[0].decode("utf-8", errors="replace").strip() if raw else ""
    header = first if first.startswith("#") else None
    if actual_hash != spec.raw_sha256:
        reject("SHA256_MISMATCH", "Exact file bytes do not match the declared source hash.")
    for field in (
        "material_id", "source_url", "source_revision", "unit_evidence",
        "normalization_evidence", "provenance_evidence",
    ):
        value = getattr(spec, field)
        if not isinstance(value, str) or not value.strip():
            reject("MISSING_EVIDENCE", f"A nonempty {field} is required.")
    for field in ("units_verified", "normalization_verified", "provenance_verified"):
        if getattr(spec, field) is not True:
            reject("UNVERIFIED_EVIDENCE", f"{field} has not been established.")
    if spec.alpha2F_convention != "dimensionless_eliashberg":
        reject("ALPHA2F_CONVENTION", "Unsupported alpha2F convention; no implicit conversion.")
    axis_scale = ENERGY_PER_AXIS_UNIT.get(spec.axis_unit)
    dos_scale = ENERGY_PER_DOS_UNIT.get(spec.phdos_unit)
    if axis_scale is None or dos_scale is None:
        reject("UNKNOWN_UNITS", "Axis and DOS units must both be explicitly supported.")
    valid_basis = (
        spec.normalization_basis in ("atom", "cell")
        and type(spec.atoms_per_basis) is int and spec.atoms_per_basis >= 1
        and (spec.normalization_basis != "atom" or spec.atoms_per_basis == 1)
    )
    if not valid_basis:
        reject("INVALID_BASIS", "Use one atom, or a cell with its integer atom count.")
    valid_density = bool(np.isfinite(spec.basis_density_m3) and spec.basis_density_m3 > 0)
    if not valid_density:
        reject("INVALID_DENSITY", "The density of the declared basis must be positive in m^-3.")

    data = None
    try:
        data = np.loadtxt(io.BytesIO(raw), comments="#", ndmin=2)
    except (ValueError, UnicodeError) as exc:
        reject("TABLE_PARSE_ERROR", str(exc))
    if data is not None:
        if data.ndim != 2 or data.shape[1] != 3 or data.shape[0] < 3:
            reject("TABLE_SHAPE", "Require at least three rows and exactly three columns.")
            data = None
        elif not np.all(np.isfinite(data)):
            reject("NONFINITE_TABLE", "All three columns must be finite.")
            data = None
    if data is not None:
        axis, alpha, dos = data.T
        metrics.update(rows=len(data), raw_axis_min=float(axis.min()), raw_axis_max=float(axis.max()),
                       negative_alpha2F_samples=int(np.sum(alpha < 0)),
                       negative_phdos_samples=int(np.sum(dos < 0)),
                       zero_width_intervals=int(np.sum(np.diff(axis) == 0)),
                       decreasing_intervals=int(np.sum(np.diff(axis) < 0)))
        monotone = bool(np.all(np.diff(axis) > 0) and axis[0] >= 0 and axis[-1] > 0)
        if not monotone:
            reject("INVALID_AXIS", "Require a nonnegative, strictly increasing axis; no sorting or merging.")
        if axis[0] != 0:
            reject("INCOMPLETE_LOW_ENERGY_SUPPORT", "The full material table must include zero energy.")
        if np.any(alpha < 0):
            reject("NEGATIVE_ALPHA2F", "Negative interaction weights are not repaired.")
        if np.any(dos < 0):
            reject("NEGATIVE_DOS", "Negative phonon modes are not repaired.")
        if np.any((dos <= 0) & (alpha > 0)):
            reject("INCOMPATIBLE_SUPPORT", "Positive alpha2F requires positive phonon DOS at the same energy.")
        if axis[0] == 0 and alpha[0] != 0:
            reject("NONZERO_ALPHA2F_AT_ZERO", "A finite lambda needs a zero alpha2F intercept.")
        # An ordered table with duplicates still has a diagnostic path integral:
        # zero-width intervals carry no area. It remains inadmissible as a
        # single-valued interpolant, and no duplicate is silently merged.
        integrable_path = bool(np.all(np.diff(axis) >= 0) and axis[0] >= 0 and axis[-1] > 0)
        if integrable_path and axis_scale is not None and dos_scale is not None:
            with np.errstate(over="ignore", invalid="ignore"):
                omega = axis * axis_scale
                density = dos / dos_scale
            if not np.all(np.isfinite(omega)) or not np.all(np.isfinite(density)):
                reject("SI_OVERFLOW", "Axis or DOS conversion to SI is not representable finitely.")
                return finish(), None
            signed = float(np.trapezoid(density, omega))
            # The sign split measures rejected mass; neither array is repaired.
            positive = float(np.trapezoid(np.maximum(density, 0), omega))
            negative = float(np.trapezoid(np.minimum(density, 0), omega))
            mask = omega > 0
            ratio = alpha[mask] / omega[mask]
            lambda_ep = float(2 * np.trapezoid(ratio, omega[mask]))
            if omega[0] == 0 and alpha[0] == 0:
                # Piecewise-linear alpha2F is proportional to E in the first bin.
                lambda_ep += float(2 * alpha[1])
            metrics.update(
                omega_max_J=float(omega[-1]), signed_modes_per_basis=signed,
                quadrature_axis_interpretation="ordered samples; zero-width intervals retained for diagnosis",
                positive_modes_per_basis=positive, negative_modes_per_basis=negative,
                negative_to_positive_fraction=(-negative / positive if positive > 0 else None),
                lambda_ep=lambda_ep,
                lambda_interpretation=("finite-grid quadrature with linear first-bin limit"
                                       if alpha[0] == 0 else "diagnostic excluding singular zero endpoint"),
            )
            if valid_basis:
                expected = 3 * spec.atoms_per_basis
                error = abs(signed / expected - 1)
                metrics.update(expected_modes_per_basis=expected, relative_mode_error=error)
                if error > policy.relative_mode_tolerance:
                    reject("MODE_COUNT", "Integrated DOS disagrees with the declared complete basis.")
            if valid_density:
                with np.errstate(over="ignore", invalid="ignore"):
                    finite_volume_dos = np.all(np.isfinite(density * spec.basis_density_m3))
                if not finite_volume_dos:
                    reject("SI_OVERFLOW", "SI mode density cannot be represented finitely.")
    return finish(), data


def audit_material_table(
    path: str | Path, spec: MaterialSpec, *, policy: AdmissionPolicy = AdmissionPolicy()
) -> MaterialAudit:
    """Inspect and report a table even when its physics/metadata are rejected."""
    return _inspect(Path(path).read_bytes(), spec, policy)[0]


def admit_material_table(
    path: str | Path, spec: MaterialSpec, *, policy: AdmissionPolicy = AdmissionPolicy()
) -> AdmittedMaterial:
    """Return SI arrays only for an admitted input; otherwise raise with its audit.

    Bytes are read once, so the validated hash and returned arrays necessarily
    describe the same read. No clipping, sorting or normalization is performed.
    """
    audit, data = _inspect(Path(path).read_bytes(), spec, policy)
    if not audit.admitted:
        raise MaterialAdmissionError(audit)
    assert data is not None
    omega = data[:, 0] * ENERGY_PER_AXIS_UNIT[spec.axis_unit]
    alpha = data[:, 1].copy()
    dos = data[:, 2] / ENERGY_PER_DOS_UNIT[spec.phdos_unit]
    volume_dos = dos * spec.basis_density_m3
    for array in (omega, alpha, dos, volume_dos):
        array.flags.writeable = False
    return AdmittedMaterial(omega, alpha, dos, volume_dos, audit)
