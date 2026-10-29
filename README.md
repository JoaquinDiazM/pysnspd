# pySNSPD 1.2.5

Conserving effective neutral BCS/KWT simulations of photon detection in SNSPDs,
with native electron/phonon populations, temperature-dependent stiffness,
depairing-calibrated diffusion and the thesis bias/readout circuit.

The publication reference is the80nm Rs608ohm/square case family at0.9K,
eta0.9,tau_ee(Tc)5ps and phonon escape20ps. The completed central775/1550nm
pairs give deterministic relative delays of3.4521ps at21.5uA and13.6161ps
at15.5uA. These are finite-window intrinsic readout results, not a full
experimental IRF, exact experimental fit or certified reset.

## Run and inspect

- [New electrical 0D phase](pipelines/25_v125_electrical_phase.py) and
  [low-bias recovery pair](pipelines/26_v125_recovery_pair.py): authenticate a
  native1D checkpoint, require independent electrical overlap and append a
  distinct electrical-only phase with the resolved DC film inductance restored.
  Historical1.2.3/1.2.4 native data stay unchanged; bath recovery is not implied.
- [Three-phase history report](plot_pipelines/10_v125_phase_history.py): read
  native2D/1D and the separate0D trace, without inventing thermal/spatial fields.
- [1.2.5 scope and checkpoint protocol](docs/releases/v1.2.5.md);
  [deferred thermal audit](docs/audits/v125_thermal_recovery_pending.md).

- [Source-frozen native single-case entry](pipelines/18_spatial_photon_run.py): reflecting
  transverse Gaussian preparation with explicit retained energy/impact position,
  or authenticated continuation of completed1.2.3 latency cases.
- [Source-frozen native campaign entry](pipelines/19_v124_campaign.py): four spatial latency
  probes or four baseline1D continuations;3threads/case; independent screens
  coordinate atomic CPU leases and print readable progress tables.
- [Individual case reports](plot_pipelines/08_dimensional_case_report.py):
  nine per-case event-based2D snapshots and six longitudinal slots.
- [Exact commands](docs/GEMINGA_COMMANDS_PRESENT.md);
  [command history](docs/GEMINGA_COMMANDS.md).
- [Release physics, assumptions and continuation protocol](docs/releases/v1.2.4.md).
- [Mathematical documents](Matematica/) describe the underlying native equations;
  the reflecting PREPARATION formula and its limitations are in the release note.

All heavy data and caches belong under /home/jdiaz/scratch. New latency cases
stop at confirmed positive0.1mV crossing+2ps, with a100ps no-trigger ceiling.
The separate recovery campaign reuses all2D data and validated1D overlap,
retains a1000nm1D domain/15ns physical target and an8h1Dsession wall budget.
Physical/domain guards may stop it earlier. A wall-budget stop is checkpointed.

The new electrical phase defaults to a100ns global horizon, using the exact
current1D checkpoint rather than rewinding to12ns. Its1ns admission overlap and
SC/electrical guards are distinct from full thermal recovery. Native p/n/gap
fields and balances remain archived;0D evolves only the three circuit states.
The historical native15/50ns horizons above do not limit the new0D observation.

## Scope and reproducibility

The native evolution kernel is unchanged from1.2.3. The frozen pipeline13 is
retained under its historical identity for source-authenticated reuse; use18/19
for the unchanged native stages, then25/26 for the new electrical phase.
Historical native manifests deliberately retain their authenticated producer
version. No population/capacity floors, tolerance relaxation, new
latency fit, statistical emulator or promised full-reset runtime are introduced.
Historical Allmaras/Korzh input controls remain explicitly isolated in16/17.

See [changelog](CHANGELOG.md), [earlier constitutive revision](docs/releases/v1.1.3.md),
and the source-hashed manifests stored alongside every run.
