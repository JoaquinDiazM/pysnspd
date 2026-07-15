"""Shared thesis plotting style and canonical figure dimensions."""

from __future__ import annotations

from math import sqrt
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib import font_manager

THESIS_STYLE_FILE = Path(__file__).with_name("thesis_sty.mplstyle")
GOLDEN_RATIO = (1.0 + sqrt(5.0)) / 2.0
THESIS_WIDTH_IN = 6.5
THESIS_BASE_HEIGHT_IN = THESIS_WIDTH_IN / GOLDEN_RATIO / 2.0
THESIS_SINGLE_FIGSIZE = (THESIS_WIDTH_IN, THESIS_BASE_HEIGHT_IN)
THESIS_DOUBLE_FIGSIZE = (THESIS_WIDTH_IN, 2.0 * THESIS_BASE_HEIGHT_IN)
THESIS_DPI = 200


def apply_thesis_style() -> None:
    """Load the packaged thesis style without requiring a LaTeX installation."""
    plt.style.use(str(THESIS_STYLE_FILE))
    plt.rcParams["text.usetex"] = False
    try:
        font_manager.findfont("Computer Modern Roman", fallback_to_default=False)
    except ValueError:
        plt.rcParams["font.serif"] = ["STIXGeneral", "DejaVu Serif"]
        plt.rcParams["mathtext.fontset"] = "stix"


__all__ = [
    "GOLDEN_RATIO",
    "THESIS_BASE_HEIGHT_IN",
    "THESIS_DOUBLE_FIGSIZE",
    "THESIS_DPI",
    "THESIS_SINGLE_FIGSIZE",
    "THESIS_STYLE_FILE",
    "THESIS_WIDTH_IN",
    "apply_thesis_style",
]
