# pySNSPD publication status

Last updated: 2026-09-23

Publication window: 2026-07-23 to 2026-10-23

Current phase: experimental stage4A has started without photons. Explicit KWT/core
parameters and a complete2D static geometry are implemented. The changed paths pass
177 tests and79 subtests on Geminga; a bounded local pilot completed in145.2s.
Two manual campaigns (40local controls,6spatial states) are prepared with progress
and ETA. No time trajectories, physical core admission or production promotion.
See [current stage4](docs/implementation/stage4/README.md) and the command notebook.

## Historical stage 3.5 closure and stage4 preparation

The user selected “Cerrar investigación 3.5 y preparar etapa 4 sin fotón”. The
[current delivery](docs/implementation/stage3_5/CURRENT.md) covers 127 variables in
15 families with coupled restrictions, observable relevance and next evidence.
The historical Allmaras extraction leaves about 7–9% electronic energy at
0.187 ps; even its phonon-only profile is not an exact single Gaussian. No
physical width is selected for Korzh. NbN retains its conditional shape, while
absolute volumetric mode normalization remains open. The
[final report](docs/implementation/stage3_5/assessment_r2_20260923/Informe_cierre_investigacion_etapa_3_5.md)
and [stage 4 contract](docs/implementation/stage4/entry_contract.json) define the
result and next tasks. No new physical trajectories or solver-parameter changes
were made, and no long job is pending. Dynamic stage 3 dependencies remain
required for the capabilities actually used; static controls can proceed.

## Release 1.0.0 and next-stage boundary

The implementation baseline is `f3c26b95ff4e4a93504371e78b46ad3a20e06273`, including
the latest persisted-data diagnostics already on GitHub. Release changes cover
version/package metadata, documentation, reproducible checks and published artifacts.
They do not replace the thesis kinetic, condensate or circuit solvers with model 0.4.

The model 0.4 candidate is not admitted for production: the core force is sensitive
to its effective scale, the condensate principal symbol is negative in demonstrated
states, and the supplied phonon DOS is not normalized with verified units. These
findings remain explicit in A-D. The pedagogical notebook is E-r02, with eight active
classes and 24 unanswered activities, independently versioned from the physics.

No long calculation is currently pending. Future computations expected to exceed
five minutes must be recorded in `/home/jdiaz/GEMINGA_COMMANDS.md` and supplied
in the chat as exact copyable commands, with their purpose, outputs and estimated
resources. The user launches them; completed batches are not automatically repeated.
The notebook preserves historical commands separately and lists lightweight checks.

## Experimental stage 2 — retained development closure

The user-authorized [development closure](docs/implementation/stage2/closure_20260922/closure_decision.json)
accepts the implemented conservative events, fixed-energy transport, KWT, BGK,
heating and escape as the basis for subsequent development. The
[final report](docs/implementation/stage2/closure_20260922/Informe_cierre_etapa_2.md)
([PDF](output/pdf/implementation/Informe_cierre_etapa_2.pdf)) consolidates the
results and remaining numerical limits.

The latest guarded SSP batch completed 21 tasks, including 13 trajectories.
All retain physical stored populations and valid balances. The largest scaled
energy defect is 5.43766e-8; the largest recorded instantaneous residual is
1.63498e-14. ONE/TWO temporal comparisons at 630 electronic and 1025 phonon
states pass against separately refined references. The 2049-phonon pair differs
by 0.004565%, exceeding its auxiliary 0.0025% budget. Its original FAIL remains;
the [strict numerical admission record](docs/implementation/stage2/stage2_admission.json)
does not become a complete dynamic mesh certificate.

## Historical stage 3 closure — preserved evidence

The [user-authorized closure](docs/implementation/stage3/closure_20260923/README.md)
accepts the static spatial work and instantaneous material/potential/KWT/heating,
three-state circuit and prescribed-reservoir balances within their recorded scope.
The completed batch contains six mixed-domain snapshots and took 355.651 s.
For the perturbed profile, medium-to-fine differences are 0.043961% in total
condensate heating and 0.96757% in maximum material speed. Those are descriptive
mesh sensitivities; no retrospective accuracy threshold was introduced.
The focused regression passes 243 tests and 46 subtests in 10.29 s.

The [final report](docs/implementation/stage3/closure_20260923/Informe_cierre_desarrollo_etapa_3_y_apertura_3_5.md)
([PDF](output/pdf/implementation/Informe_cierre_desarrollo_etapa_3_y_apertura_3_5.pdf))
and [current sequence](docs/implementation/SECUENCIA_VIGENTE.md) distinguish
development closure from completion of the original stage 3 contract. The latter,
temporal admission, full D.27, kinetic-interface admission and production promotion
remain false. A current-dependent reservoir, conservative equal-energy interface
transport and a weak coupled trajectory with integrated balances remain required
before dynamic studies using those capabilities. The current stage 4 entry permits static non-photon diagnostics first. No detector pulse, latency or hotbelt formation is claimed.

The original [stage 3.5 opening](docs/implementation/stage3_5/README.md) began documentary research into
the model's physical domain of confidence: 127 parameter entries in 15 families,
with separate Korzh and Allmaras experimental/model sources. No ranges or sweep
were adopted at that opening. L2D/W = 1.5–6 remains only an example, not a universal validity interval; transverse 1D
sufficiency must be tested during the future transient. Physical scope is distinct
from a statistical confidence interval. No long calculation is pending. That historical delivery checker is
`sandbox/stage3_spatial/closure_20260923/verify_delivery.py`;
the completed six-case batch must not be repeated. `v1.0.0` remains unchanged.

The earlier [RK4/SSP review](docs/implementation/stage2/review_20260922/README.md)
and [scope review](docs/implementation/stage2/practical_review_20260922/validation_scope_review.md)
are historical checkpoints. Their failures, incomplete attempts and then-pending
manual plans are preserved rather than presented as the current work queue.

## Experimental stage 1 closure — retained result

The [final report](docs/implementation/stage1_closure/Informe_cierre_etapa_1.md)
and [admission record](docs/implementation/stage1_closure/closure_admission.json)
close electronic stage 1 with reduced catalogue-use tests: BGK, heating,
prescribed spectral work, a self-consistent synthetic condensate cell and
fixed-field two-cell transport. All 53 electronic gates pass. A fresh independent
R2 reassessment also passes; its 87-file historical delivery remains unchanged.
The full repository regression passes 317 tests in 34.69 s on Geminga.

Transport conserves its explicitly reconstructed energy to 2.22e-16. At 513
shared-energy nodes its native-R2 energy bias is 0.02119%, thermal remapping error
0.06879% and quasiparticle count drift 0.00233%; these distinct defects are not
silently projected away. Native R2 energy is used directly in BGK, heating and
the synthetic condensate test. Time is not calibrated to a physical NbN rate.

The NbN preprocessing decision is updated: an explicit, traceable common-support
derived shape is retained for restricted experimental work. It is not admitted
as absolute SI material or as an equivalent global physical phonon kernel.
Thermal U/C changes are small at the tested points through 80 K, but the coupling
changes in an internal phonon gap and for nonthermal populations are significant.
The primary-source audit does not settle this file's absolute DOS normalization.

The [next-stage contract](docs/implementation/stage1_closure/NEXT_STAGE.md)
starts with shared electron-phonon reaction events in cells, then simultaneous
condensate/transport coupling. Full D.4 item 2, spatial admission and circuit
validation remain outstanding. No production activation, long calculation or
change to v1.0.0 was made. All useful commands are recorded in Geminga's notebook.

## Experimental stage 1, R2 — retained historical result

The [R2 report](docs/implementation/stage1_r2/Informe_etapa_1_r2.md) and
[source-bound admission certificate](docs/implementation/stage1_r2/catalog_admission.json)
close the sampled uniform electronic catalogue for synthetic-cell validation.
The final file passes 234 physical cases, 150 derivative checks, 476 support checks
(including all 326 retained failure points), and 6,855 dense queries. Exact cubic
ordering minima are positive in 9,260 cells at 65 sampled amplitudes. These are
sampled numerical guarantees, not a uniform bound for arbitrary populations.

R2 represents energy in Gamma/absDelta, retains compensated differences near
zero pair breaking and inserts seven diagnosed ratio nodes. The critical slope
error falls from 8.56756% total in R1 to 0.03196%; the worst tested relative current
response error is 0.05134%, below the fixed 0.1% tolerance. Interpolation,
quadrature, cutoff and causal-regulator errors have separate evidence. No change
to the uniform physical functional was required by this numerical diagnosis.

The full suite passes 278 tests in 34.09 s. Geminga construction took 143.89 s
plus 4.61 s and 1.09 s for local refinements; occupied queries cost 0.0537 ms.
These timings exclude coupled dynamics. No calculation exceeded the compute
handoff limit. The next permitted step is synthetic one- or two-cell conservation,
relaxation and time-refinement checks from D.4.2, not a full production transient.

The material input remains **REJECTED**. The improved plots distinguish Simon's
phonon DOS from the electronic Usadel DOS and expose what the legacy loader clips.
Public provenance of the numerical body is verified; units and normalization
remain unverified, and negative/duplicate phonon data are unresolved. Spatial
admissibility and core sensitivity remain later-stage physical checks.

The historical first-iteration section below is retained for traceability. Its
manifest must be checked against source commit d17d7c3, since R2 changes the shared
experimental module. The release tag v1.0.0 remains unchanged at 5ea0cd6.

## Experimental stage 1 after v1.0.0 — historical first iteration

Explicit material admission and uniform vacuum/fixed-occupation catalogues are
implemented in `pysnspd.experimental`, independently of production. The final
Geminga regression has 232 passing tests (57 new), in 22.41 s; catalogue construction
and diagnostics took 27.37 s. Query timing is 0.048 ms for the occupied catalogue.

Stage 1 is **not closed for promotion to stage 2**. The input NbN data have 361
negative DOS samples and four conflicting duplicate frequency pairs; source units
and normalization remain unverified. Hermite interpolation fixes the identified
small-current slope error, but a separate, resolved reference demonstrates an
8.61% numerical-regulator bias for a low-energy population at eta/Delta0=0.001.
The next required work is joint edge-quadrature, Gamma-grid and regulator
convergence, plus admissible phonon data. This result does not require changing
the underlying uniform functional on its own.

See `docs/implementation/stage1/Informe_etapa_1.md` and `catalog_admission.json`.
The published tag v1.0.0 remains at 5ea0cd6; no new production transient was run.
Remote-only responses to the old v0.3 learning exercises were preserved separately
under `docs/learning_history`, without transferring answers to new questions.

The sections below preserve the earlier operational and scientific checkpoints.

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

## D3 energy-projection checkpoint

D3 now diagnoses the electronic spectral-storage term directly from the saved
full photon trajectory; it does not rerun or alter the solver. In the central
100 nm, the full spectral term is transiently large and almost reversible, its
finite-step split is dominated by the gap contribution rather than the
spectral-q contribution, and the omitted-energy residual follows the full
spectral term closely. These observations justify further diagnosis, but not
activation of independent `P_delta` and `P_q` sources.

The present kinetic catalogue reaches only `q*xi = 1.221`. D3 records a maximum
catalogue-clipped central fraction of `28.8%` for `41.2 ps`; no clipping is
seen in `Te`, `Tph`, or `Delta`. The snapshot-wise maximum reaches an isolated
`q*xi = 9.832`, while 99% of those temporal maxima remain below `2.792`.
`GEMINGA_COMMANDS.md` therefore queues the user-run-only PRE
`pre_qwide_01`, with `gamma_max_fraction=5.0`, nominal coverage to
`q*xi ~= 3.05`, 64 DOS/phase q nodes, and 351 stiffness q nodes. This expands
the range by 2.5 while retaining the nominal q-point density. It deliberately
does not size the expensive phase/power catalogue around the single extreme;
the follow-up D3 run must report the residual clipping before any physics is
changed.

The framework update is archived at this checkpoint. When resumed, the first
candidate is a conservative electronic-energy update in `u_e`, compared
against the existing temperature update without changing the rest of the
splitting. Independent `P_delta`/`P_q` sources remain blocked until the wider
catalogue, temporal refinement, closed-cycle energy balance, and
gTDGL/circuit double-counting checks are complete.

Before that work, the commission changes currently pending in
`/home/jdiaz/memoria` must be incorporated. Changing the scientific model while
the submitted manuscript is being corrected would obscure which equations and
results the supervisors are reviewing. The thesis working tree is therefore
left untouched by this checkpoint.

The undergraduate-thesis simulation reference is commit `523c53d` (the last
pySNSPD commit before the 20 July delivery). The later publication checkpoint
`23ea557` changes plotting and adds `status.md`, but does not change the
simulation implementation relative to `523c53d`. The commit that archives D3
is a later diagnostic/implementation descendant and is intentionally not the
final thesis-reference commit: pending visualization changes and the commission
edits must be incorporated before that reference is frozen.

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

The E3 center/edge comparison now uses four default snapshots at `50, 51, 52,
53 ps`, adds `|q| xi` to the matched physical fields, and emits matched power
and energy/heat-capacity atlases. Every color scale is shared between positions
and uses exact extrema from all persisted snapshots in both runs; unchanged
single-run E3 manifests supply fingerprint-validated limits so only the four
selected maps are reconstructed. A conditional four-panel numerical comparison
shows recovery margins, final residual/tolerance ratios, circuit modes, and the
available window. On the reference pair, `xi = 8.591 nm`, the central and edge
latencies are `5.872` and `6.743 ps`, and both recoveries remain censored. The
documented lateral run name was corrected from the nonexistent `_eg` suffix to
the completed `_edge` run.

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

- Complete suite on Geminga at the D3 checkpoint: `157 passed`. The focused
  energy-projection diagnostic set passes both tests, the PRE CLI exposes the
  explicit `--dos-n-q` override, and `compileall` succeeds without launching a
  solver or production catalogue.
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
- The E3 position comparison was regenerated from the central and lateral
  30 uA references. Its five one-page PDFs (fields, circuit, power,
  energy/heat capacity, and four-panel censored recovery) were rendered and
  visually inspected; the focused comparison/photon diagnostic set passes all
  `4` tests.
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

1. The user runs `pre_qwide_01` from `GEMINGA_COMMANDS.md`; no agent launches
   this production PRE. Verify `Status: OK`, the realized q axes, file sizes,
   and absence of host-memory pressure before using it downstream.
2. Incorporate and review the commission changes in `/home/jdiaz/memoria`,
   including the pending visualization edits, before changing framework
   science or freezing the final undergraduate-thesis reference commit.
3. Resume D3 only after the wider PRE exists. Recompute catalogue coverage,
   the `P_delta`/`P_q` split, path sensitivity, accumulated energy closure, and
   the possible overlap with gTDGL/circuit energy. Decide between a conservative
   `u_e` update and explicit sources only from those checks.
4. Rerun the isothermal/no-circuit sweep under the new 150 ps name recorded in
   `GEMINGA_COMMANDS.md`: 15 currents (adding 40 and 44 uA around the transition),
   301 snapshots, 14 child workers plus the base process, and one thread per
   case. This leaves one of Geminga's 16 physical cores nominally free while
   reducing retained field history to 15% of the failed 2000-snapshot run. The
   13-case `_01` run ended with `BrokenProcessPool` and is not reusable as a
   sweep.
5. Generate and visually validate the I-V curve only after Z2 reports adequate
   completed-current coverage; the present `_01` PDFs are coverage diagnostics,
   not an accepted I-V result.
6. Continue subsequent development directly on `main`; the former feature
   branch history is already preserved upstream.
7. Start the planned mesh/time/thermal convergence campaign using the new mesh
   edge-length and triangle-quality baselines.

The wider publication gates remain unchanged: mesh/time/thermal convergence,
time-resolved current continuity, accumulated energy closure, material
provenance, threshold sensitivity, and long-run latency/recovery acceptance
are still open.
