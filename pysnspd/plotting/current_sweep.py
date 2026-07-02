"""Current-sweep plotting entry points.

This module is intentionally a placeholder for the first Z1 audit step.  The
current task is to verify that a Z-series pipeline can discover and index all
raw run files, summaries, manifests and NPZ keys before any physics-level sweep
figure is implemented.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Mapping, Sequence


def plot_current_sweep_placeholder(
    records: Sequence[Mapping[str, Any]],
    output_dir: str | Path,
    *,
    dpi: int = 480,
) -> int:
    """Placeholder for future current-sweep figures.

    Parameters are accepted so the Z1 pipeline already has the final call shape,
    but no figure is produced yet.
    """
    _ = records
    _ = Path(output_dir)
    _ = int(dpi)
    return 0


__all__ = ["plot_current_sweep_placeholder"]
