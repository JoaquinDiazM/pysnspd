"""Future photon-energy-deposition layer.

The active repository does not yet implement the spatial phonon-bubble profile
or the thermal injection rule. The official PHOTON-run is still a placeholder
interface and records photon inputs directly from command-line parameters.

This module intentionally exposes no public runtime helpers until the photon
deposition OE is implemented and validated against the coupled no-photon
stationary state.

Previous placeholder functions returned ``0`` and were removed because a
zero-valued photon-deposition model can silently hide missing physics.
"""

from __future__ import annotations

__all__: list[str] = []