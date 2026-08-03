# pySNSPD publication status

Last updated: 2026-08-03

Publication window: 2026-07-23 to 2026-10-23

Current phase: Week 2 closure - stationarity, localized artifacts, and E1/E2/E3
figure consolidation

## Executive status

pySNSPD is a functional multiscale research prototype coupling dirty-limit
Usadel material data, finite-volume gTDGL/Poisson dynamics, two-temperature
evolution, photon deposition, and circuit observables in physical units.

The stationarity and localized low-amplitude artifact tasks have reached a
satisfactory technical closure. The corrected stiffness catalogue has been
used by the production PRE and long SS runs, and the former temporary isolated
diagnostic pipeline has been removed. The remaining Week 2 acceptance step is
visual inspection of the I-V curve after the current sweep running in screen
`code3` completes.

The inspected 30 uA, 200 ps stationary result is photon-ready under the new
strict-fixed-point OR weak-dynamic-attractor policy. Contact recovery, current
continuity, phase continuation, thermal fields, and the circuit remain hard
validity gates.

## Stationarity closure

- Dynamic stationarity now uses a physical 5 ps tail rather than a fixed count
  of snapshots. The reference run supplies 51 frames over 5.002 ps.
- Its late profile drift is 0.0082% against a 1% limit and its voltage span is
  1.697% against a 2% limit.
- Thermal readiness is evaluated from stored Te/Tph fields rather than the
  instantaneous explicit RHS. The reference RMS thermal drift is 0.1018%
  against a 0.3% limit; p99 nodal drift is 1.376 mK against 3 mK.
- Early termination additionally requires five consecutive successful
  evaluations separated by 0.5 ps. The old eta-residual stop and obsolete
  command/config aliases have been removed.
- Reanalysis of the legacy 30 uA run with the current policy first reaches and
  latches `photon_ready` at 104.452 ps. Its stored legacy value remains `null`
  and is reported separately rather than rewritten.

## Low-amplitude constitutive result

For fixed temperature and superfluid momentum, dirty-limit Matsubara theory
requires `j_s = O(|Delta|^2 q)` as the condensate vanishes. Interpolating the
current itself between zero and the first positive amplitude node instead
forces an unphysical linear law.

The isolated diagnostic established:

- exponent `1` for the former current interpolation and exponent `2` for both
  direct Matsubara evaluation and the stiffness formulation;
- less than `0.1%` stiffness error below the first positive amplitude node
  across the audited temperature/momentum cases;
- reduction of the synthetic notch phase-source amplification from about
  `117x` to `1.0018x` the direct result, with relative RMS error reduced from
  about `20.8` to `1.8e-3`.

Together with the long-run field inspection, this closes the low-amplitude
artifact mechanism at the current scope. Mesh/time convergence remains a
separate publication requirement.

## Implemented regularization

- PRE now stores the finite Matsubara stiffness
  `kappa(Te, |Delta|^2, |q|)` and its `|Delta|^2` axis, while retaining the
  simultaneously computed current-density table as a diagnostic resource.
- The exact analytic Matsubara limit is used at `|Delta|=0`; nonfinite or
  nonpositive stiffness data are rejected rather than repaired silently.
- SS and photon runs require the new stiffness contract. Current-only PRE
  catalogues are intentionally incompatible and must be regenerated.
- Temporal solvers interpolate stiffness in `(Te, |Delta|^2, |q|)` and form
  current from the regular gauge-invariant edge pair flow
  `Im(conj(Delta_i) U_ij Delta_j) / ell_ij`.
- An edge incident on an exact zero of the order parameter carries exactly zero
  pair flow. The gauge-invariant phase gradient is used only as the stiffness
  table coordinate.
- The Usadel–GL current difference is formed on each edge before applying one
  finite-volume divergence, avoiding cancellation between two large nodal
  divergences.
- Harmonic continuation of the phase quotient now begins at
  `64*sqrt(machine epsilon)` by default; the audited sensitivity range is
  `16–256`, replacing the former physical-amplitude cutoff.
- Nonfinite Allmaras forcing raises an error. The adaptive integrator rejects
  that attempt, reduces `dt`, and retries instead of inserting zeros into the
  field.
- The Laplacian, finite-volume geometry, and local quadratic `|psi|^2` update
  were deliberately left unchanged.

E1 plotting now reports the Dynes spectral broadening used by each DOS curve,
marks q_c on the Usadel supercurrent branch, and stores the mesh presentation
and quality diagnostics under the run's `mesh/` folder. The temporary
low-amplitude plotting path and its dedicated analysis code/tests were removed.

E2 plotting now separates solver procedure from physical stationarity, removes
dense snapshot markers from histories, and uses paired electron/phonon axes
where their scales differ. The final longitudinal profile includes total,
superfluid, normal-current fits, order parameter, potential, Te, and Tph.
Snapshot atlases use the seven standard times from 0 to 200 ps and include
current-direction arrows, four power-density channels (including finite-volume
thermal diffusion), and electronic/phononic energy densities and heat
capacities. New runs persist diffusion power directly; old runs reconstruct it
exactly from stored fields and operators without rerunning the simulation.

E3 photon plotting now derives every field-atlas color scale from all persisted
snapshots, while rendering only the requested seven times. The exact full-run
extrema are fingerprinted in the figure manifest and reused on unchanged data,
so later plotting passes do not repeat the multi-gigabyte scan. Photon power
and energy/heat-capacity atlases use the same reconstructed constitutive and
finite-volume quantities as E2. A conditional four-panel recovery diagnostic
is emitted only for detected but right-censored runs. For the 30 uA reference,
it reports `t_lat = 5.872 ps`, `t_rec > 1450 ps`, final normalized recovery
residuals of `0.0097x` (bias current), `2.1x` (strip current), `12x` (readout
current), `3.0x` (terminal voltage), `4.7x` (capacitor voltage), and `0.2x`
(TDGL voltage), and overdamped circuit decay times of `4791.24`, `210.70`, and
`99.06 ps`. The post-photon window covers only `30.3%` of the slowest mode, so
the absence of `t_rec` is diagnosed as insufficient horizon rather than silently
treated as failed physics.

Z2 current-sweep analysis now performs a shallow inventory by default and
loads only `ss_summary.yaml`, `stationary_state.npz`, and the shared PRE mesh
for completed endpoints. It no longer opens multi-gigabyte relaxation or
snapshot archives merely to construct a final I-V point. On the `_01` sweep,
the complete pipeline dropped from more than three minutes without finishing
to `0.92 s` wall time and `84.9 MB` peak RSS. The new four-panel regime summary
separates missing data from physical gate failures, reports exact sampled
currents for strict SS, dynamic SS, photon readiness, and an electrical ohmic
approximation, and stores the same classification in YAML.

The completed `code3` sweep is not a completed physical sweep: only the 20 uA
base case has a summary and final state. The 25 uA history is not a valid NPZ;
the 30 uA history reaches 200 ps but has no summary; the remaining ten cases
retain only their seeds. All twelve are classified as unavailable rather than
failed stationarity. Consequently there is no defensible I-V curve or current
range yet. At 20 uA, the central and terminal voltages are only `0.00470` and
`0.1916` of their respective normal-state references, mean
`|Delta|/Delta0 = 0.9193`, and strict SS, dynamic SS, photon readiness, and the
10% two-voltage ohmic criterion all fail. The source of the previously zero
terminal ratio was also corrected: Z2 now prioritizes the final solver voltage
(`3.91089 mV`) over the zero-voltage analytic seed stored earlier in the YAML.

The `code3` failure has now been correlated with Geminga's host telemetry. The
machine did not reboot, but from 06:20 to 06:32 UTC its RAM usage rose from
76.4% to 98.1%, swap reached 100%, load reached 51.5, and committed memory
reached 103.3% while several cases serialized multi-gigabyte NPZ archives. The
kernel counters contain OOM kills, although the unprivileged account cannot
timestamp the exact victim from the journal. This makes an OOM-killed worker,
followed by `ProcessPoolExecutor` invalidating every outstanding future, the
high-confidence cause of `BrokenProcessPool`. The simultaneous Grafana
`DatasourceNoData` alert is a symptom of host memory/I/O starvation, not a
credible cause: an SSH disconnect cannot terminate work running in `screen`.
The earlier 16 ps sweeps succeeded because their completed case directories are
only about 241--247 MB, versus about 26 GB for the completed 200 ps base case.

## Validation completed

- Complete suite on Geminga after the Z2 update: `154 passed`; the focused
  photon/E2 plotting set passes all `5` tests.
- The SS early-stop callback now forwards the public phase-gradient tolerance
  names correctly. Its focused adapter/stop/target regression set passes all
  `18` tests without launching a production run.
- E1 was regenerated from the production PRE and its six affected PDFs were
  rendered and visually inspected without clipping, overlap, or illegible text.
- E2 was regenerated from
  `ss_phasecg_I30uA_200ps_circuitthermal_stiffness_01`; all seven replacement
  PDFs were rendered and visually inspected. The three obsolete combined
  outputs are removed only after their replacements save successfully.
- E3 was regenerated from
  `photon_phasecg_I30uA_0p8eV_sigma10nm_t50ps_1500ps_stiffness_01`; the scalar
  history, globally scaled field atlas, power atlas, energy/heat-capacity atlas,
  and conditional censored-recovery diagnostic were rendered and visually
  inspected. E2 snapshot arrows were also regenerated and checked with one
  shared current scale.
- Z2 focused tests pass all `10` cases. Its central I-V, terminal I-V, and
  four-panel coverage/regime PDF were regenerated and visually inspected from
  the partial `_01` sweep without opening its large histories.
- Focused constitutive suite: exact zero-amplitude Matsubara limit, quadratic
  amplitude power, q parity, stiffness interpolation, discrete plane wave,
  exact-zero edge current, gauge invariance, difference-before-divergence,
  finite phase-drive continuation over the `16–256` tolerance range, and
  nonfinite-step rejection.
- Static reachability audit: 85 library modules, 851 definitions, all 851
  definitions reachable, and zero resolution failures. Only empty package
  `__init__.py` modules remain unreachable.
- A new small PRE built stiffness, current, phase-space, and power resources in
  parallel.
- A 3 fs SS smoke completed with circuit and thermal coupling, the v2
  stiffness closure, three accepted steps, zero rejected steps, exact requested
  final time, and complete persistence.
- A 6 fs photon smoke inherited that SS state, deposited a 0.8 eV photon at
  3 fs, evolved circuit and thermal state, reached the requested horizon, and
  wrote all state, history, snapshot, timing, summary, and manifest outputs.
- E3 generated its complete diagnostic figure set from those archives.

These runs validate interfaces, discrete identities, orchestration, and
serialization only. Their femtosecond horizons are not physical stationarity,
detection, recovery, or notch-removal evidence.

## Immediate work and acceptance gates

1. Rerun the isothermal/no-circuit sweep under the new 150 ps name recorded in
   `GEMINGA_COMMANDS.md`: 15 currents (adding 40 and 44 uA around the transition),
   301 snapshots, 14 child workers plus the base process, and one thread per
   case. This leaves one of Geminga's 16 physical cores nominally free while
   reducing retained field history to 15% of the failed 2000-snapshot run. The
   13-case `_01` run ended with `BrokenProcessPool` and is not reusable as a
   sweep.
2. Generate and visually validate the I-V curve only after Z2 reports adequate
   completed-current coverage; the present `_01` PDFs are coverage diagnostics,
   not an accepted I-V result.
3. Continue subsequent development directly on `main`; the former feature
   branch history is already preserved upstream.
4. Start the planned mesh/time/thermal convergence campaign using the new mesh
   edge-length and triangle-quality baselines.

The wider publication gates remain unchanged: mesh/time/thermal convergence,
time-resolved current continuity, accumulated energy closure, material
provenance, threshold sensitivity, and long-run latency/recovery acceptance
are still open.
