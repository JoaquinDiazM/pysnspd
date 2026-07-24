# pySNSPD publication status

Last updated: 2026-07-24

Publication window: 2026-07-23 to 2026-10-23

Current phase: Week 1 - detection and recovery criteria

Baseline branch: `main`

Baseline parent: `523c53d`

Baseline commit: `23ea557657d890b1c902f5962a669d3fb845fd93`

## Executive status

pySNSPD is a functional research prototype that couples material-specific
Usadel catalogues, electron-phonon power tables, finite-volume gTDGL/Poisson
dynamics, two-temperature evolution, photon deposition, and circuit-level
observables in physical units. The current scientific contribution is the
auditable multiscale chain; the publication objective is to turn that chain
into a numerically defensible and experimentally comparable prediction of
detection threshold and latency.

The production tree on Geminga was frozen before the publication campaign in
commit `23ea557`. The first post-baseline change reorganized the library
without introducing a new physics model or running an expensive PRE, SS, or
photon case. It removed pipeline-unreachable implementation, made package
ownership explicit, and preserved the complete production-reachable call
graph. Plotting was then unified under the thesis style in commit `618103c`.
The baseline-freeze task is complete; the active Week 1 task is to freeze
operational definitions of detection latency and recovery before implementing
them.

## Week 1 architecture cleanup record

The deletion and placement policy is now versioned in
`docs/ARCHITECTURE_POLICY.md` and enforced by the audit utilities under
`tools/`. The reviewed reachability graph starts from all scripts in
`pipelines/` and `plot_pipelines/`, follows direct and transitive references,
and conservatively accounts for dynamically typed method/property access.

The cleanup:

- removed 205 unreachable definitions and seven fully dead modules rather than
  moving them to a legacy folder;
- consolidated all finite-volume meshes, geometry, compatibility objects, and
  discrete operators under `pysnspd.mesh`;
- separated orchestration into `pysnspd.solver`, thermal evolution into
  `pysnspd.thermal`, and photon deposition into `pysnspd.excitation`;
- retained gTDGL constitutive physics under `pysnspd.gtdgl`;
- split the oversized plotting and solver modules by responsibility;
- replaced package-level compatibility re-exports with imports from defining
  modules;
- merged the two surviving helpers from `kinetic/powers.py` into
  `kinetic/power_table.py`.

The final static audit covers 80 library modules and 765 definitions. All 765
production definitions are reachable, Jedi reports zero resolution failures,
and the only unreachable modules are intentionally empty package
`__init__.py` files. No production library module exceeds 800 lines; the
largest is 784 lines.

Test maintenance removed five obsolete test cases: three that protected the
flat `gtdgl` facade and two that protected unused mesh/device convenience
behavior. Assertions for two deleted operator helpers and the deleted
power-table reload API were also removed from otherwise surviving tests.
Imports in the remaining tests were updated to their owning packages.

Validation on Geminga:

- `compileall`: passed for library, pipelines, plot pipelines, tests, and
  tools;
- pipeline smoke tests: all 11 production and plotting entry points accept
  `--help`;
- pytest: `103 passed in 13.26s`;
- no PRE, SS, photon, sweep, or publication plot was regenerated.

## Plot-style unification record

All figure-producing modules under `pysnspd/plotting/` now load
`thesis_sty.mplstyle` through `apply_thesis_style()`. Public plotting functions
and all seven scripts in `plot_pipelines/` use `THESIS_DPI` as their default,
and figure widths come from the shared thesis constants rather than independent
numeric widths. The deprecated style TODOs in pipelines 01, 02, and 03 were
removed.

The migration covered the remaining legacy producers for photon, stationary
summary, memory, snapshot-power, snapshot-grid, and adaptive-step figures.
Explicit 14-point overrides in the legacy photon and SS-memory plots were
removed so labels, ticks, and legends inherit the common typography. Compact
multi-panel labels remain explicit only where the canonical page width
requires them.

`tests/test_plotting_style_policy.py` now rejects figure producers that do not
apply the shared style, numeric DPI defaults, independent numeric figure
widths, divergent pipeline DPI defaults, and reintroduced migration TODOs. A
new pipeline-03 plotting smoke test covers its five-column scalar snapshot
figure.

Validation on Geminga:

- style policy: `6 passed`;
- plotting smoke tests: all seven entry points accept `--help`;
- complete pytest suite: `110 passed in 13.31s`;
- synthetic renders for SS overview, SS relaxation/adaptive diagnostics,
  snapshot power maps, runtime metrics, and photon scalar maps were visually
  checked for clipping, overlap, and legibility;
- no production PRE, SS, photon, sweep, or publication dataset was modified.

## Baseline acceptance update

On 2026-07-24 the production smoke and test results were accepted by the
project owner. The additional manual validation passed except for the long
16 ps, 80 ps, and 800 ps normal-pipeline runs. Those runs are intentionally
deferred until the latency/recovery implementation, when their outputs will
exercise the new criteria. This is a planned validation deferral rather than a
baseline-freeze blocker.

## Detection and recovery definition study

Status: criteria proposed from literature and an existing 800 ps production
run; no implementation has been started.

The production case
`photon_phasecg_I30uA_0p8eV_sigma10nm_t50ps_800ps_01` is sufficient for
post-processing. Its lightweight history stores time, circuit observables,
thermal and condensate summaries, and a one-row `photon_applied` marker; its
snapshot archive stores the spatial mesoscopic and thermal fields. Therefore
latency and recovery can be recomputed under revised criteria without rerunning
the coupled simulation.

### Proposed latency contract

The canonical operational latency is

`t_lat = t_cross(V_out, V_threshold) - t_gamma`,

where `t_gamma` is the recorded photon-deposition time and `t_cross` is the
first leading-edge crossing of a declared physical output-voltage threshold.
The extractor should:

- estimate the pre-photon baseline robustly and apply an explicit pulse
  polarity rather than silently using an absolute value;
- linearly interpolate between stored samples;
- require the expected edge direction plus a short confirmation
  window/hysteresis to reject chatter;
- record the threshold value, reference plane, baseline window, polarity,
  interpolation method, and confirmation rule beside every result;
- return a detected flag and a null/right-censored latency when the threshold
  is not crossed before the stored trajectory ends.

An absolute leading-edge threshold is required to classify detection because
it represents the comparator and can reject a physically negligible pulse.
Constant-fraction times should also be reported for waveform comparison and
future jitter studies, but not used alone as the detection classifier: a
fraction of a pulse's own future maximum is noncausal for early termination
and would classify arbitrarily small pulses as detections. Constant-fraction
timing is nevertheless valuable because it reduces amplitude-dependent
time-walk.

The same event extractor naturally extends to stochastic physics. A
deterministic run yields one latency; an ensemble will yield a detection
probability and the conditional distribution of `t_lat`, including quantiles
and jitter, while non-detections remain censored observations rather than
being assigned the simulation end time.

In the inspected run, `t_gamma = 50 ps` and the baseline-subtracted
`V_out` peak is approximately `0.942 mV` at `115.1 ps`. The resulting latency
depends materially on the declared threshold: about `2.26`, `4.09`, `5.66`,
and `19.49 ps` for `10`, `50`, `100`, and `500 microV`, respectively; the
10% and 50% constant-fraction values are about `5.49 ps` and `18.24 ps`.
Consequently, a latency number without its threshold definition is not an
auditable result.

### Proposed recovery contract

The literature's operational reset time generally tracks current recovery or,
more directly, recovery of detection efficiency. pySNSPD should publish that
comparable quantity and keep it distinct from a stricter full-state recovery:

1. `t_rec,DE90`: first time the inferred or directly sampled detection
   efficiency reaches 90% of its pre-event value. This is the primary
   experiment-comparable reset metric once a defensible `DE(I)` map exists.
2. `t_rec,electrical`: return of detector/readout observables to their declared
   circuit tolerance, including `I_b`, `I_s`, `I_rf`, `V_out`, `V_tdgl`, the
   internal capacitor state, and circuit-equation residuals.
3. `t_rec,state`: return of circuit, thermal, and mesoscopic state families to
   the accepted steady-state neighbourhood. This is the strict metric for
   early termination and numerical claims.

For each state variable or field diagnostic `x_i`, define a dimensionless
acceptance residual from a declared absolute-plus-relative tolerance around
the pre-photon stationary reference. Tolerances must be larger than resolved
baseline noise and discretization error and must be accompanied by a
sensitivity table; normalization only by the photon-induced peak is unsuitable
as the final rule because peaks can be singular or arbitrarily small. Spatial
fields should use an RMS/quantile norm plus a hard maximum guard, and raw phase
or scalar potential should be replaced by gauge-invariant or centered
quantities.

Let the family residual be the maximum of its normalized component residuals.
Recovery is the earliest post-event time at which every required family is
inside tolerance for a declared hold window and its rates/residual equations
also indicate stationarity. If the stored run ends before this can be
confirmed, the result is `not_recovered` with
`t_rec > t_end - t_gamma`; the last timestamp must not be reported as a
recovery time. A non-detecting transient may still have a relaxation time, but
should not silently be labeled detector reset.

Online early termination should use the same contract on lightweight scalar
history and solver residuals, with sparse spatial checks. Plotting pipelines
03 and E3 should independently recompute the more complete post-processed
metrics from saved history/snapshots, report the criterion and censored state,
and annotate the threshold crossing and recovery acceptance interval.

The 800 ps example illustrates why the tiers are necessary. `V_out` finishes
only about `6.5 microV` from its pre-photon baseline and the plotted thermal
maxima are within 1% of their photon-induced excursions about `210 ps` after
absorption. However, the unplotted capacitor state `v_c` finishes roughly
`41.2 microV` from its initial value, near its largest excursion, and the mean
gap remains about 6.2% of its peak excursion away from baseline. The existing
E3 figure therefore supports near-recovery of the displayed output variables,
but the strict all-state result is presently right-censored:
`t_rec,state > 750 ps` under even the provisional peak-normalized audit.

### Literature basis

- [Allmaras et al., *Intrinsic Timing Jitter and Latency in Superconducting
  Nanowire Single-photon Detectors* (2019)](https://doi.org/10.1103/PhysRevApplied.11.034062)
  connects intrinsic jitter to fluctuations of microscopic detection latency.
- [Korzh et al., *Demonstration of sub-3 ps temporal resolution with a
  superconducting nanowire single-photon detector*
  (2020)](https://doi.org/10.1038/s41566-020-0589-x) experimentally probes
  detection latency and its material dependence.
- [Schuck et al., *Waveguide integrated low noise NbTiN nanowire
  single-photon detectors with milli-Hz dark count rate*
  (2013)](https://doi.org/10.1038/srep01893) uses a trigger near half pulse
  height at the maximum rising-edge slope.
- [Gras et al., *Fast single-photon detectors and real-time key distillation
  enable high secret-key-rate quantum key distribution systems*
  (2023)](https://doi.org/10.1038/s41566-023-01168-2) shows why
  constant-fraction discrimination reduces amplitude-dependent timing error.
- [Kerman et al., *Kinetic-inductance-limited reset time of superconducting
  nanowire photon counters*
  (2006)](https://doi.org/10.1063/1.2183810) relates current recovery to kinetic
  inductance and reports reset through 90% recovery of detection efficiency.
- [Annunziata et al., *Reset dynamics and latching in niobium superconducting
  nanowire single-photon detectors*
  (2010)](https://doi.org/10.1063/1.3498809) separates hotspot cooling from the
  slower inductive current reset and identifies the latching condition.
- [Burenkov et al., *Investigations of afterpulsing and detection efficiency
  recovery in superconducting nanowire single-photon detectors*
  (2013)](https://doi.org/10.1063/1.4807833) demonstrates that readout dynamics
  can make recovery non-monotonic.
- [Autebert et al., *Direct measurement of the recovery time of SNSPDs and its
  application for quantum communication*
  (2019)](https://doi.org/10.1364/QIM.2019.S1D.3) treats recovery of detection
  efficiency as the operational quantity.
- [Wang et al., *Timing Jitter Induced by Stochastic Baseline Fluctuations in
  High-Count-Rate SNSPDs*
  (accepted 2026)](https://doi.org/10.1103/8yf7-blyh) shows that finite-memory
  readout baselines directly perturb threshold-extracted times in the
  stochastic/high-rate regime.

## Frozen change inventory

The production diff at the freeze contains nine plotting modules, with 83
insertions and 68 deletions:

- `current_sweep.py`, `eliashberg_spectrum.py`, `mesh.py`,
  `pre_diagnostics.py`, `usadel_dos_curves.py`, and `usadel_gap.py` standardize
  scientific units using bracketed labels;
- `photon_comparison.py` and `ss_phasecg_figures.py` replace repeated panel
  labels with shared figure axes, adjust layout margins, and clarify
  dimensionless/count/rate labels;
- `power_diagnostics.py` standardizes units, simplifies legends and state
  labels, and introduces explicit plotting ranges and a symmetric-log power
  scale;
- current-conversion and transient figures receive small color, legend, and
  annotation refinements.

The canonical plotting entry point remains
`pysnspd/plotting/thesis_sty.mplstyle`, loaded through
`pysnspd.plotting.style.apply_thesis_style()`. Migration is intentionally
incomplete at this baseline.

## Freeze acceptance record

| Criterion | Status | Evidence or limitation |
| --- | --- | --- |
| Production tree captured in Git | Done | This baseline commit contains all nine pre-existing production modifications and this report. |
| Baseline published on `origin/main` | Done in this delivery | This commit is the publication baseline; verify it with `git status --branch` and `git log -1`. |
| Scope documented | Done | The change inventory above separates existing plot work from future refactors. |
| Defect register initialized | Done | Scientific, numerical, reproducibility, and maintenance gaps are listed below. |
| Simulation results regenerated | Not run | No PRE, SS, PHOTON, sweep, or plotting pipeline was executed for this documentation-only freeze. |
| Numerical validation completed | Open | Convergence and conservation targets belong to Weeks 1-2. |

## Defect and risk register

Priority meanings: P0 blocks defensible publication results; P1 can invalidate
or materially weaken a central claim; P2 is reproducibility or maintenance
debt; P3 is cleanup that can wait until the scientific path is stable.

| ID | Priority | Status | Defect, risk, or limitation | Closure criterion |
| --- | --- | --- | --- | --- |
| NUM-001 | P0 | Definition proposed | Detection latency, experiment-comparable reset, and strict all-state recovery are not yet implemented as separate, auditable metrics. | Accept and implement the documented threshold, censoring, hold-window, tolerance, and sensitivity contracts; classifications remain stable under accepted numerical refinements. |
| NUM-002 | P0 | Open | Convergence is not yet demonstrated jointly for mesh, time step, thermal subcycling, thermal-domain size, and terminal length. | Refined runs change detection current by less than 2-3% and latency by less than 5%. |
| CONS-001 | P0 | Open | Accumulated energy closure is not yet demonstrated; omitted or reduced terms such as `P_Delta` and `P_q` must be included or bounded. | Energy imbalance remains below 1-2%, with every omitted term quantified. |
| CONS-002 | P1 | Open | Current-continuity and energy errors are not yet reported through the most violent part of the transient. | Time-resolved residuals are plotted and remain below declared tolerances. |
| THERM-001 | P1 | Open | A true thermal steady state before photon injection has not yet been demonstrated for the publication cases. | Pre-injection thermal rates and residuals satisfy a documented stationary threshold. |
| GTDGL-001 | P1 | Open | Center and edge cases have used non-identical phase-solver tolerances, leaving a numerical confounder in spatial comparisons. | All spatial cases use one justified tolerance policy and retain their classification after refinement. |
| PHYS-001 | P1 | Open | Small notch-like structures appear before or far from the photon impact in the thesis diagnostics. | The structures disappear under correction/refinement or are explained and shown to be physical. |
| MAT-001 | P1 | Open | `D_eff = 1.581 cm^2 s^-1` is calibrated to the critical-current scale and then reused in transport predictions, creating potential circularity. | Parameter provenance is explicit and `D_eff` is independently constrained or propagated as an uncertainty interval. |
| MAT-002 | P1 | Open | Cross-consistency among `D`, `sigma_n`, sheet resistance, `N(0)`, `T_c`, thickness, and `I_c` is not yet presented as one audit. | A single material table reports source, uncertainty, inference path, and consistency checks for every parameter. |
| MAT-003 | P1 | Open | The phonon escape time is not yet constrained by data or a defensible material/interface range. | `tau_esc` is tied to evidence or treated in a sensitivity analysis. |
| PERF-001 | P1 | Open | A representative trajectory costs roughly 40-45 hours, making the full threshold map infeasible without acceleration. | Early termination, continuation, bisection, and coarse-to-fine search provide an effective 5-10x campaign speedup. |
| PLOT-001 | P2 | Closed | Every plotting producer uses the shared thesis style API; all plotting pipelines share one DPI default and no deprecated style TODO remains. | Closed by the plot-style unification record above, including static enforcement and representative visual review. |
| PLOT-002 | P2 | Open | `power_diagnostics.py` contains hard-coded axis ranges, commented alternatives, and simplified labels whose general validity is not documented. | Dataset-driven defaults are restored or the publication-specific choices become named, documented options. |
| DOC-001 | P2 | Open | README examples and implementation have drifted, including the documented `Z1_current_sweep_analysis.py` versus tracked `Z2_current_sweep_analysis.py`. | All documented commands are checked against tracked entry points. |
| DOC-002 | P2 | Open | README, thesis, appendix, and implementation still require synchronization. | Every publication equation and workflow statement maps to current code and a stable reference. |
| PKG-001 | P2 | Open | `pyproject.toml` has an incomplete dependency declaration and is not a complete installation recipe. | A clean environment can install the package and run the lightweight validation path from documented commands. |
| REPRO-001 | P2 | Open | Machine-specific Geminga configurations and data paths remain mixed with public examples. | Public templates are portable; machine-local overrides are clearly isolated and ignored where appropriate. |
| REPRO-002 | P2 | Open | No small public demo dataset or lightweight end-to-end smoke case exists. | A versioned, documented example exercises the official workflow without production-scale data. |
| ARCH-001 | P3 | Closed | Production reachability is classified, dead code is deleted, package ownership is explicit, and all production library modules are below 800 lines. | Closed by the Week 1 architecture cleanup; policy and validation evidence are recorded above. |
| PHYS-002 | P3 | Scoped out | Full spatial nonthermal kinetics and stochastic observables are not implemented. | Reserve for a second article after the deterministic threshold/latency result is established. |

## Publication acceptance targets

The primary claim to support is:

> A material-specific, low-temperature model based on Usadel-consistent gTDGL
> can predict detection threshold, latency, and absorption-position dependence
> while closing the circuit loop, without a phenomenological hotspot
> resistance or independent fitting at each wavelength.

The campaign should target:

- detection-current variation below 2-3% under numerical refinement;
- latency variation below 5%;
- accumulated energy error below 1-2%;
- median experimental detection-current error below 5-10%, using one global
  parameter set rather than per-curve fits;
- correct trends with bias current, photon energy, and absorption position;
- an explicit ablation showing what changes when microscopic spectra,
  Usadel-consistent gTDGL, thermal closure, and circuit feedback are enabled.

## Thirteen-week roadmap

| Week | Dates | Main objective | Exit artifact |
| --- | --- | --- | --- |
| 1 | Jul 23-29 | Freeze baseline; initialize defects; define automatic detection/recovery criteria. | Published baseline and living `status.md`. |
| 2 | Jul 30-Aug 5 | Build convergence matrix and material/provenance audit. | Numerical protocol and parameter table. |
| 3 | Aug 6-12 | Add early termination, threshold bisection, continuation, and coarse search. | Timed threshold-search workflow. |
| 4 | Aug 13-19 | Sweep photon energy and absorption position on the coarse grid. | Preliminary `I_det(E_gamma, y_gamma)` map. |
| 5 | Aug 20-26 | Refine boundary points and freeze the defensible threshold map. | Converged central result. |
| 6 | Aug 27-Sep 2 | Extract internal latency, first phase slip, dissipative-band formation, and circuit pulse time. | Latency dataset and definitions. |
| 7 | Sep 3-9 | Run controlled physical ablations with fixed geometry and parameters. | Ablation figure and causal interpretation. |
| 8 | Sep 10-16 | Assemble the Allmaras/Korzh comparison and uncertainty model. | Experimental benchmark table. |
| 9 | Sep 17-23 | Complete comparison and uncertainty analysis; keep work robust to reduced meeting time around Chilean holidays. | Comparison figure and decision memo. |
| 10 | Sep 24-30 | Write the paper around the quantitative result, not the software architecture. | Full main-text draft and six-figure plan. |
| 11 | Oct 1-7 | Build the reproducible supplement and move detailed audits out of the main narrative. | Supplement draft and reproduction guide. |
| 12 | Oct 8-14 | Internal scientific review; remove unsupported claims. | Reviewer-style issue list resolved. |
| 13 | Oct 15-21 | Final revision, figure consistency, references, and submission package. | Submission candidate. |
| Close | Oct 22-23 | Protected buffer for final checks and submission. | Submitted manuscript. |

## Scope guardrails

Do not delay the first paper by adding full stochastic physics, intrinsic jitter,
dark counts, a full spatial evolution of `f(E, r, t)` and `n(Omega, r, t)`,
many materials, or an unsupported ab initio claim. Establishing the
deterministic, material-specific threshold and latency result first creates the
validated baseline needed for those additions to form a coherent second paper.

## Update protocol

Update this file whenever a defect changes state, a numerical criterion is
accepted, a roadmap gate is passed, or a publication claim changes. Every
update should identify the evidence (commit, figure, table, run manifest, or
analysis note) and should preserve closed items rather than deleting their
history.
