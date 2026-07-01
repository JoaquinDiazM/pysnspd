"""Analysis helpers for pySNSPD raw simulation outputs.

This package contains data-loading and numerical reduction utilities only.
Figure styling and matplotlib code live in :mod:`pysnspd.plotting`.
"""

from .ss_run import SSRunData, build_ss_plot_dataset, load_ss_run, summarize_ss_npz_contents

__all__ = [
    "SSRunData",
    "build_ss_plot_dataset",
    "load_ss_run",
    "summarize_ss_npz_contents",
]
