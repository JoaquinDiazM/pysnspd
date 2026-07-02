"""Deprecated stationary-initial-guess module.

The active stationary SS seed is implemented in

    pysnspd.gtdgl.seed.build_stationary_seed

and is called directly by ``pipelines/02_ss_run_template.py``.

This module used to contain placeholder helpers that returned ``0``.
Those placeholders were intentionally removed during repository cleanup
because silent zero-valued physics stubs are dangerous in a simulation code.
"""

from __future__ import annotations

__all__: list[str] = []