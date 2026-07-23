"""Apply the reviewed package-move import map across source and tests."""

from __future__ import annotations

import argparse
from pathlib import Path


REPLACEMENTS = {
    "pysnspd.gtdgl.finite_volume": "pysnspd.mesh.finite_volume",
    "pysnspd.gtdgl.geometry": "pysnspd.mesh.geometry",
    "pysnspd.gtdgl.tdgl_compat": "pysnspd.mesh.tdgl_compat",
    "pysnspd.gtdgl.device": "pysnspd.mesh.device",
    "pysnspd.gtdgl.operators": "pysnspd.mesh.operators",
    "pysnspd.gtdgl.adapter": "pysnspd.solver.stationary",
    "pysnspd.gtdgl.solver": "pysnspd.solver.core",
    "pysnspd.gtdgl.ss_targets": "pysnspd.solver.targets",
    "pysnspd.gtdgl.diagnostics": "pysnspd.solver.diagnostics",
    "pysnspd.gtdgl.options": "pysnspd.solver.options",
    "pysnspd.gtdgl.seed": "pysnspd.solver.seed",
    "pysnspd.gtdgl.state_io": "pysnspd.solver.state_io",
    "pysnspd.gtdgl.photon_transient": "pysnspd.solver.transient",
    "pysnspd.gtdgl.thermal": "pysnspd.thermal.evolution",
    "pysnspd.gtdgl.photon": "pysnspd.excitation.photon",
    "pysnspd.gtdgl.snapshot_diagnostics": "pysnspd.analysis.snapshots",
    "pysnspd.plotting.ss_power_figures import _snapshot_diffusion_power_density": (
        "pysnspd.plotting.ss_power_helpers import _snapshot_diffusion_power_density"
    ),
    "pysnspd.kinetic.powers": "pysnspd.kinetic.power_table",
    "from .currents import": "from pysnspd.gtdgl.currents import",
    "from .usadel_current import": "from pysnspd.gtdgl.usadel_current import",
    "from .allmaras import": "from pysnspd.gtdgl.allmaras import",
    "from .device import": "from pysnspd.mesh.device import",
    "from .tdgl_operators import": "from pysnspd.gtdgl.tdgl_operators import",
    "from .solver import TDGLSolver": "from .core import TDGLSolver",
    "from .thermal import": "from pysnspd.thermal.evolution import",
    "from .ss_targets import": "from .targets import",
    "from pysnspd.gtdgl import solve_stationary_pytdgl_like": (
        "from pysnspd.solver.stationary import solve_stationary_pytdgl_like"
    ),
}


def rewrite(root: Path) -> None:
    paths: list[Path] = []
    for directory in ("pysnspd", "pipelines", "plot_pipelines", "tests"):
        paths.extend((root / directory).rglob("*.py"))
    for path in sorted(paths):
        original = path.read_text(encoding="utf-8")
        updated = original
        for old, new in REPLACEMENTS.items():
            updated = updated.replace(old, new)
        if updated != original:
            path.write_text(updated, encoding="utf-8")
            print(path.relative_to(root).as_posix())


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path.cwd())
    args = parser.parse_args()
    rewrite(args.root.resolve())


if __name__ == "__main__":
    main()
