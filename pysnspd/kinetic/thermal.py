"""Future two-temperature thermal evolution layer.

The active repository currently provides kinetic catalogues and projected
electron--phonon power diagnostics, but not yet the final spatially coupled
two-temperature time evolution.

This module intentionally exposes no public runtime helpers until the thermal
OE is implemented with a physically audited energy functional, heat capacities,
diffusion, escape and coupling to the gTDGL fields.

Previous placeholder functions returned ``0`` and were removed because a
zero-valued thermal update can silently hide missing physics.
"""

from __future__ import annotations

__all__: list[str] = []