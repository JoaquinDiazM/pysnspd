"""Scalar diagnostics and seed/material helpers for OE7 relaxation."""
from __future__ import annotations

import numpy as np

from pysnspd.gtdgl.material import GTDGLMaterial
from pysnspd.gtdgl.state import CurrentFields


def current_residual(
    currents: CurrentFields,
    mesh,
    material: GTDGLMaterial | None = None,
    target_current_A: float | None = None,
) -> float:
    """Dimensionless RMS residual of div(j_tot).

    New notebook-order runs pass ``material`` and ``target_current_A`` so the
    scale is the imposed average current density divided by the mesh spacing.
    Older smoke tests called this with only ``(currents, mesh)``; in that case
    we fall back to the reconstructed total-current RMS scale.
    """
    div = np.asarray(currents.node_div_jtot_A_m3, dtype=float)
    h = float(getattr(mesh, "target_spacing_m", getattr(mesh, "xi_mesh_m", 1.0e-9)))
    if material is not None and target_current_A is not None:
        jscale = abs(target_current_density_A_m2(material, float(target_current_A)))
    else:
        jscale = float(
            np.sqrt(
                np.nanmean(
                    currents.node_jtot_x_A_m2**2 + currents.node_jtot_y_A_m2**2
                )
            )
        )
    scale = max(jscale / max(h, 1.0e-300), 1.0)
    return float(np.sqrt(np.nanmean(div * div)) / scale)


def max_current_residual(
    currents: CurrentFields,
    mesh,
    material: GTDGLMaterial | None = None,
    target_current_A: float | None = None,
) -> float:
    div = np.asarray(currents.node_div_jtot_A_m3, dtype=float)
    h = float(getattr(mesh, "target_spacing_m", getattr(mesh, "xi_mesh_m", 1.0e-9)))
    if material is not None and target_current_A is not None:
        jscale = abs(target_current_density_A_m2(material, float(target_current_A)))
    else:
        jscale = float(
            np.sqrt(
                np.nanmean(
                    currents.node_jtot_x_A_m2**2 + currents.node_jtot_y_A_m2**2
                )
            )
        )
    scale = max(jscale / max(h, 1.0e-300), 1.0)
    return float(np.nanmax(np.abs(div)) / scale)


def current_density_maxima_A_m2(currents: CurrentFields) -> tuple[float, float]:
    """Return max |j_n| and max |j_tot| from edge fields."""
    jn_max = float(np.nanmax(np.abs(currents.edge_jn_A_m2))) if currents.edge_jn_A_m2.size else 0.0
    jt_max = float(np.nanmax(np.abs(currents.edge_jtot_A_m2))) if currents.edge_jtot_A_m2.size else 0.0
    return jn_max, jt_max


def seed_target_current_A(seed) -> float:
    """Extract imposed transport current from an OE6 seed-like object."""
    for name in ("I_bias_A", "target_current_A", "current_A"):
        if hasattr(seed, name):
            value = float(getattr(seed, name))
            if np.isfinite(value):
                return value
    metadata = getattr(seed, "metadata", None)
    if isinstance(metadata, dict):
        for name in ("I_bias_A", "target_current_A", "current_A"):
            if name in metadata:
                value = float(metadata[name])
                if np.isfinite(value):
                    return value
    return 0.0


def target_current_density_A_m2(material: GTDGLMaterial, target_current_A: float) -> float:
    return float(target_current_A) / max(material.width_m * material.thickness_m, 1.0e-300)
