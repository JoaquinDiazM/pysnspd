# Pending audit: native low-temperature recovery

Recorded 2026-10-29; deferred by the user. No thermal equation, tolerance,
capacity guard, boundary condition or population was modified for v1.2.5.

The low-bias Rs608 reference continuations remain superconducting and their
electrical tails are closely reproduced by a lumped inductive three-state
circuit, while maximum effective electron temperature remains above the 0.9 K
bath. This electrical agreement is not thermal equilibration.

| Case (reference only) | Native endpoint, ns | Maximum effective Te, K |
|---|---:|---:|
| 15.5 uA, Edep 1.4398164979339385 eV (775 nm) |16.066646941749794|2.277898209488907|
| 15.5 uA, Edep 0.7199082489669693 eV (1550 nm)|15.54047788384909|2.2978768800280758|

Source: `/home/jdiaz/scratch/photon_runs/v124_base_recovery_20261023_01/runs/`;
resumable session coordinator: `v124_lowbias_recovery50ns_20261028_01`.
The exact native manifests, checkpoints, gap/p/n fields, work/energy balances,
monitor values and authenticated source pins remain unchanged. Every v1.2.5 0D
continuation archives copies of the final native manifest/checkpoint and records
`thermal_audit_pending=true`, `thermal_evolved=false`, and
`thermal_reset_certified=false`.

Future audit questions include the effective-temperature diagnostic in sparse
gapped populations, its relation to current unloading and stiffness, prescribed
instantaneous-current boundary work, low-temperature rates, energy/population
relaxation and sensitivity to the native spectral/domain discretization. These
are questions, not identified causes. No audit or remediation is part of this
release. Never replace the stored Te curve by an invented decay to the bath or
read a completed electrical horizon as full physical reset.
