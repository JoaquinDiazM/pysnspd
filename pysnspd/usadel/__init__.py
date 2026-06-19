"""
Usadel/material spectral tools for pySNSPD.
"""

from pysnspd.usadel.parameters import (
    E_CHARGE_C,
    HBAR_J_S,
    K_B_J_K,
    MaterialParameters,
    bcs_gap_J,
    bcs_gap_zero_J,
    depairing_energy_grid_J,
    energy_axis_J,
    material_parameters_from_config,
    q_axis_from_depairing_energy_m_inv,
)
from pysnspd.usadel.solver import (
    anomalous_proxy,
    compute_dos_grid,
    dos_diagnostics,
    dynes_bcs_dos,
)
from pysnspd.usadel.catalog import (
    UsadelCatalog,
    J_to_meV,
    build_usadel_catalog_from_config,
    catalog_summary,
    load_usadel_catalog_npz,
    meV_axis,
    save_usadel_catalog_npz,
)

__all__ = [
    "E_CHARGE_C",
    "HBAR_J_S",
    "K_B_J_K",
    "MaterialParameters",
    "bcs_gap_J",
    "bcs_gap_zero_J",
    "depairing_energy_grid_J",
    "energy_axis_J",
    "material_parameters_from_config",
    "q_axis_from_depairing_energy_m_inv",
    "anomalous_proxy",
    "compute_dos_grid",
    "dos_diagnostics",
    "dynes_bcs_dos",
    "UsadelCatalog",
    "J_to_meV",
    "build_usadel_catalog_from_config",
    "catalog_summary",
    "load_usadel_catalog_npz",
    "meV_axis",
    "save_usadel_catalog_npz",
]