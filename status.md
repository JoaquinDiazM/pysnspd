# pySNSPD publication status

Last updated: 2026-07-23

Publication window: 2026-07-23 to 2026-10-23

Current phase: Week 1 - baseline freeze

Baseline branch: `main`

Baseline parent: `523c53d`

Baseline commit: the commit that first introduces this file

## Executive status

pySNSPD is a functional research prototype that couples material-specific
Usadel catalogues, electron-phonon power tables, finite-volume gTDGL/Poisson
dynamics, two-temperature evolution, photon deposition, and circuit-level
observables in physical units. The current scientific contribution is the
auditable multiscale chain; the publication objective is to turn that chain
into a numerically defensible and experimentally comparable prediction of
detection threshold and latency.

This baseline freezes the production tree on Geminga before the publication
campaign begins. It includes the current plotting work but does not include a
new physics run, a refactor, a directory reorganization, or a complete
migration of every plot to the thesis style. Those are tracked below and must
be introduced as explicit post-baseline changes.

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
| NUM-001 | P0 | Open | Detection/recovery classification is not yet frozen as one automatic, tolerance-independent criterion. | A documented classifier gives the same outcome under the accepted numerical refinements. |
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
| PLOT-001 | P2 | Active | Thesis-style migration is partial; several plotting pipelines explicitly retain deprecated local styling. | Every publication figure uses the shared style API and passes a visual consistency review. |
| PLOT-002 | P2 | Open | `power_diagnostics.py` contains hard-coded axis ranges, commented alternatives, and simplified labels whose general validity is not documented. | Dataset-driven defaults are restored or the publication-specific choices become named, documented options. |
| DOC-001 | P2 | Open | README examples and implementation have drifted, including the documented `Z1_current_sweep_analysis.py` versus tracked `Z2_current_sweep_analysis.py`. | All documented commands are checked against tracked entry points. |
| DOC-002 | P2 | Open | README, thesis, appendix, and implementation still require synchronization. | Every publication equation and workflow statement maps to current code and a stable reference. |
| PKG-001 | P2 | Open | `pyproject.toml` has an incomplete dependency declaration and is not a complete installation recipe. | A clean environment can install the package and run the lightweight validation path from documented commands. |
| REPRO-001 | P2 | Open | Machine-specific Geminga configurations and data paths remain mixed with public examples. | Public templates are portable; machine-local overrides are clearly isolated and ignored where appropriate. |
| REPRO-002 | P2 | Open | No small public demo dataset or lightweight end-to-end smoke case exists. | A versioned, documented example exercises the official workflow without production-scale data. |
| ARCH-001 | P3 | Open | Official, diagnostic, exploratory, and legacy code is not yet classified; several modules exceed 1,000 lines. | Files are inventoried first, then reorganized without changing baseline results. |
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
