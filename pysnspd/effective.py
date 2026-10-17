"""Shared baseline-compatible API of the effective conserving SNSPD evolution.

The evolved state contains the complex gap, native quasiparticle occupations,
phonon occupations and the complete three-state bias/readout circuit. This
mesoscopic neutral closure uses one shared BCS/radial-gradient energy. It is an
explicit effective hypothesis, not an exact spatial Usadel or charge-imbalance
solver. Electron and phonon energy-equivalent temperatures are diagnostics;
the evolving distributions are not replaced by thermal populations.

Allmaras-2020 tau_ee(Tc) and tau_ep(Tc) are free material calibrators. A
deterministic post-cascade transient does not by itself predict stochastic
jitter or the optical cascade delay.

Current temperature-dependent workflows (v1.1.3):
    pipelines/09_thermal_stiffness_photon_run.py
    pipelines/10_thermal_stiffness_basecases.py
    plot_pipelines/07_case_report.py

The classes exported here retain the fixed-bath baseline API for compatibility.
The current thermal functional is ResummedThermalFunctional and its coupled
system is ThermalCoupledSystem in pysnspd.experimental; the current workflows
select them explicitly. Baseline pipeline05 is a historical/control workflow.

The implementation modules retain their experimental location for historical
imports; this module is the public entry point for the effective model.
"""
from .experimental.effective_bcs_functional import (
    EffectiveBCSFunctional, EffectiveBCSBulk, bcs_count_quadrature,
)
from .experimental.effective_bcs_preparation import (
    EffectiveBCSGeometry, EffectiveBCSDC, get_effective_geometry,
    prepare_effective_dc, save_effective_dc,
)
from .experimental.effective_bcs_kinetics import (
    EffectiveBCSKinetics, EffectiveBCSKineticResult, acoustic_phonon_grid,
)
from .experimental.effective_bcs_ports import (
    EffectiveNeutralPorts, freeze_exterior_inductance,
)
from .experimental.effective_coupled import (
    EffectiveState, EffectiveRHS, EffectiveCoupledSystem,
)
from .experimental.effective_material import (
    EffectivePhononMaterial, allmaras_acoustic_material,
)
from .experimental.electrical_ports import ThesisCircuitParameters
from .experimental.trigger_window import TriggerWindow, TriggerWindowConfig

__all__ = [
    "EffectiveBCSFunctional", "EffectiveBCSBulk", "bcs_count_quadrature",
    "EffectiveBCSGeometry", "EffectiveBCSDC", "get_effective_geometry",
    "prepare_effective_dc", "save_effective_dc", "EffectiveBCSKinetics",
    "EffectiveBCSKineticResult", "acoustic_phonon_grid", "EffectiveNeutralPorts",
    "freeze_exterior_inductance", "EffectiveState", "EffectiveRHS",
    "EffectiveCoupledSystem", "EffectivePhononMaterial",
    "allmaras_acoustic_material", "ThesisCircuitParameters",
    "TriggerWindow", "TriggerWindowConfig",
]
