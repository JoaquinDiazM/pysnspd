# Repository working policy

- Work directly on `main` for routine `pysnspd` tasks.
- Do not create a feature branch or pull request unless the user explicitly asks for one.

## Geminga computation policy

- Use the non-administrator `jdiaz` account. No `sudo` or administrator access.
- Parallelize independent cases and expensive independent spectral queries.
  Recheck CPU topology, effective affinity/quotas and available memory at runtime.
  Allocate at most 90% of the available computational resources; reserve complete
  physical cores where possible. The 2026-09-23 inventory is 16 physical cores,
  32 logical CPUs and about 123 GiB RAM: at most 28 logical CPUs, leaving two
  complete cores free. This is a ceiling, not a requirement to fill idle slots.
- Count nested workers and numerical-library threads within one shared budget.
  Pin the task to the admitted CPU set, limit BLAS/OpenMP threads per worker,
  and reduce concurrency when memory or a quota is tighter. Report the actual
  worker count and progress/ETA. Never assume that independent cases need to run
  successively or that an untested parallel speedup makes a long job lightweight.
- Do not launch computations known or reasonably expected to take more than five minutes.
- Write those commands, their purpose, expected outputs and estimated resources to
  `/home/jdiaz/GEMINGA_COMMANDS.md`. Preserve the existing entries.
- Also include the exact long-run command explicitly in the chat, in a copyable
  code block. Do not refer the user only to the command notebook.
- Leave useful commands the user may want to try in that same file, including
  reproducible lightweight diagnostics.
- For a required long calculation, stop the dependent work and wait for the user
  to run it and return the terminal output with `reinicia` or an equivalent request.
  Do not use polling loops, scheduled wakeups or background jobs to bypass this handoff.
- Use a bounded timeout for lightweight checks with uncertain runtime. A timed-out
  check is incomplete; record the unrestricted user-run command instead of restarting
  it repeatedly or splitting a long calculation to evade the five-minute limit.
- By explicit user instruction on2026-10-15, thesis-v1.0.0 preserves the former thesis/model0.4 tag and v1.0.0 freezes the effective neutral BCS/KWT baseline. Never move the preserved thesis alias. Version1.1.0 is the subsequent temperature-dependent-stiffness/calibrated-diffusion revision.

## Practical numerical acceptance

- Prefer revising an unnecessarily strict tolerance to postponing a small
  discrepancy. Tie the revised margin to the observable, its signal scale and
  the purpose of the calculation; record both the original result and the new
  acceptance decision. A tiny residual signal need not block development.
- Do not use tolerance changes to conceal a large physical error, nonfinite
  states, instability, missing coupling terms or material conservation defects.
  Explain such failures and repair the cause before promoting that capability.
- Reuse published and inherited algorithms before adding a new solver. Prefer
  the thesis/pyTDGL KWT local Euler update for the next temporal comparison.
  Its quadratic amplitude solve does not make it second-order in time; establish
  accuracy by step refinement of the relevant observables.
- Reports must label each plotted quantity, units, reference subtraction, norm,
  normalization and physical case explicitly.

## Initial photon-transient observation window

- The first target is relative detection latency for the 80 nm Korzh (2020)
  device. Evolve through the agreed `V_out` trigger plus a declared extra window;
  full nanosecond-scale recovery is not the initial target.
- This changes only the observation/integration horizon. Preserve the same
  physical equations, material and device parameters, deposition, mesh, boundary
  conditions, full thesis circuit and numerical accuracy policy when comparing
  short and extended windows. Do not shorten circuit time constants or discard
  slow physical states to make the window shorter.
- Record a first crossing from accepted time steps with an identified voltage
  observable, threshold and crossing confirmation; record the extra margin.
  Do not invent these physical measurement choices. An event that never triggers
  before its finite time ceiling is censored, not assigned a fictitious latency.
- The existing production `latency` early-stop mode waits for a peak; it is not
  yet the requested trigger-plus-margin semantics. Production changes require
  the later implementation step, not a silent change in these static controls.

## Interior-strip DC preparation and device scope

- The next campaign prepares a photon-free DC equilibrium in an interior section
  of the 80 nm strip. Its longitudinal cuts are continuations of the same
  current-carrying superconductor, not metallic contacts or zero-current BCS
  reservoirs. Prepare the gap, spectral fields, distributions, potential and
  circuit consistently at the same operating current.
- For the ideal straight strip, require a spatially uniform gap magnitude and
  longitudinal current up to measured discretization error. The phase has a
  finite gradient q; compare lengths at equal current, not equal total phase.
  Preserve the Delaunay-Voronoi mesh and inherited first-order KWT Euler update.
- Use a constant DC source and a future single localized photon. Do not add new
  AC campaigns or unrelated artificial perturbations as mandatory intermediate
  stages. Existing AC results remain historical coupling diagnostics.
- Exclude longitudinal geometric jitter. Hold the longitudinal absorption
  coordinate fixed; future position dependence concerns the strip width.
- Preserve the thesis three-state circuit. For the Korzh-motivated scenario,
  freeze the exterior inductance once per DC reference as 96 nH plus the
  differential bulk inductance of the unmodeled part of the 5 micrometer active
  strip. The 96 nH describes only the added exterior inductor; never subtract
  the resolved segment from it or double count resolved condensate dynamics.
  Do not introduce a dynamic lumped Lk(t) or recompute the partition during a pulse.
- A photon-free equilibrium with zero dissipated power does not establish the
  missing finite-energy heat balance, absolute NbN rates or photon preparation.
  Report that scope explicitly before admitting a photon transient.

## Heavy-data storage policy (user, 2026-10-08)

- Store all new heavy numerical data, temporary builds, caches, rendering images
  and campaign logs under `/home/jdiaz/scratch`. Source code, small manifests,
  scientific decisions and final report PDFs may remain in the repository.
- `/home/jdiaz/pysnspd/tmp` redirects to `/home/jdiaz/scratch/tmp/pysnspd`;
  `/home/jdiaz/tmp` redirects to `/home/jdiaz/scratch/tmp/jdiaz`. Preserve these
  compatibility symlinks; do not recreate real temporary directories in home.
- Preserve historical results. Relocations must verify contents and retain
  compatible references; never use `git reset` or `git clean` as storage cleanup.

## Scoped user authorization, 2026-10-11

The user explicitly authorizes this task to execute runs longer than five minutes, iterate until the complete improved-model update, prepare a final PDF report and push the completed work. This supersedes the earlier handoff rule for this task only. Keep the global 90% resource limit, jdiaz-only access and scratch storage. Heavy jobs remain foreground, observable and checkpointed; no automation is needed. A conserving effective neutral BCS/KWT closure is being developed after the full first-Moyal candidate failed; it requires its own DC and photon validation and must not be called an exact microscopic Usadel evolution.


## Active model standard (user, 2026-10-15)

- Use the latest implemented model version for new normal simulations and case reports; the current standard is v1.1.3, with DC-matched finite-depairing thermal correction and depairing-calibrated diffusion.
- Use older versions only for explicitly identified historical analyses or scientific controls. Keep version/closure identity visible; the current v1.1.3 quantitative photon-latency admission remains pending.

## Signed heat correction (2026-10-15, v1.1.1)

- Electron-temperature runs use directional FD/Pauli heat activities, energy normalization and nonnegative excess entropy. Do not restore the former exact dual-moment heat carrier: it can cool empty states.
- Preserve the positive-capacity guard, actual circuit, calibrated diffusion and source hashes. Never resume a v1.1.0 trajectory under the changed v1.1.1 equations; initialize fresh cases in a new scratch campaign.

## Numerical clock and observation policy (user task, 2026-10-16)

- New thermal campaigns use v1.1.3 and a 100 ps no-trigger ceiling; retain crossing confirmation and 5 ps extra time. This is observation censoring, not full reset.
- Never accept an untouched nonlinear predictor merely because h times the RHS fits an absolute tolerance. Require at least one coupled correction and scale thermal nonlinear tolerance by min(1,h/0.02 ps).
- Test the complete thermal population-to-gap response within convergence, not solely after returning from the solver. Internal iterates require Ceff>0 without a floor.
- The declared minimum adaptive step is 1e-6 ps by default. Stop with checkpoint below that resolution; never replace the required step with the limit or count it as detection. Retain numerical-resolution failures as failed/incomplete model evidence.
- Preserve 1.1.1 frozen trajectories; restart new cases from fresh DC/source under the changed numerical acceptance policy.


## Active closure and performance scope (2026-10-16, v1.1.3)

- Use the DC-matched finite-depairing electron-temperature correction for new normal photon cases; quadratic and fixed-bath rigidity are explicit scientific controls.
- State the temperature-independent matching elasticity and nonzero fixed-field high-Te limit. Native kinetic DOS remains q=0 BCS; Qeff includes amplitude gradients and does not imply microscopic Usadel evolution.
- Keep physical residual,tolerances,support,positive capacity/principal stiffness and independent work/energy ledgers when optimizing. Conservation alone does not establish latency validation.
- Changed physics requires fresh cases with exact sources. Independent diagnostic segments identify their scope and never redeposit the photon.
- Set NUMBA_CACHE_DIR and all heavy outputs/caches under scratch.
- Approximately ten wall-clock hours per100 ps after DC is a performance objective. Report measured horizon,setup/JIT inclusion,affinity and hardware; a linear projection is not a guarantee.


## Dimensional evolution standard (user task, 2026-10-19, v1.2.0)

- New normal cases use pipelines/13_dimensional_photon_run.py and the v1.2 phase-aware archive. The v1.1.3 physical kernel remains source-frozen for its running historical campaign.
- Default interior length is 300 nm. Confirm the first positive 0.1 mV crossing for 0.5 ps; official 2D ends at crossing + 2 ps, including confirmation.
- A native-state 1D shadow overlaps the final 2D interval. Require projection/reconstruction, independent voltage/current/energy and resolved-winding checks. Initialize official 1D by conservatively reprojecting the last 2D, including all native occupations and the actual circuit; never redeposit or force FD/BE.
- The same 300 nm strip and frozen exterior inductance remain after reduction. No tails or new hotbelt constitutive law are implied. Full-1D evolution to 15 ns is explicitly authorized, subject to dimensional/domain/physical guards; it is not a certified full reset.
- Integrate absolute left/right cut energies separately and monitor their fraction of retained photon energy plus absolute device work; the default 2% bound is an engineering truncation indicator, not a rigorous exterior-state error bound.
- Use temporal step doubling in 1D with independent local energy acceptance. Save actual 1D snapshot times without forcing tiny time steps for plotting. Preserve diagnostic states rather than loosening a failed guard.
- Six-case campaigns admit 24 Numba threads on twelve whole physical cores (four SMT threads / two cores per case), plus a separate coordinator and at least two fully reserved cores. Inventory quotas/RAM and known active runs before admission; the legacy global campaign lock prevents overlapping campaigns. Do not claim a measured 70% use or an untested speedup.


## Extended 1D scope (user task, 2026-10-19, v1.2.1)

- The default 1D domain is 1000 nm, centered on the 300 nm 2D lens (350 nm per added side); spacing remains 5 nm /201nodes. Keep the full native central state. Added regions use an explicit cold, instantaneous-current stable-branch hypothesis with matched DC-relative phase; initial exterior history is unknown and must pass overlap/domain/energy checks.
- Freeze both exterior-inductance representations from the same original DC. Subtract only L'_DC times added resolved length from the lumped circuit at representation change; preserve all added96nH. Transfer the corresponding excess inductive energy to the resolved film; never hide it in arbitrary reference work.
- Seed the 1D DC with its exact finite-link bulk dispersion, not a continuum phase guess across a long chain; keep force/current/nV tolerances unchanged.
- Check free25nm outer collars, accepted current/voltage ceilings, native support/capacity, independent energy and global branch polarity. Moderate negative readout voltage is allowed for capacitive undershoot.
- Preserve state and explicit stop reason for domain escape, electrical guard or possible latching. Latching is a possible physical outcome, not proof of a numerical error. Early recovery requires1ns of jointly stable gap/native populations/electron and phonon moments/full circuit after2ns; otherwise the15ns target remains subject to all guards.
- Never resume a v1.2.0 checkpoint under v1.2.1. Use a fresh campaign directory and record all geometry, initialization, partition and monitor inputs.


## Active v1.2.2 standard (user task, 2026-10-21)

- New normal cases use v1.2.2. Preserve the20-source native material/kinetic/circuit kernel and historical data. The1D outer condition is an explicitly declared isothermal stable bulk continuation at instantaneous Is, including prescribed boundary work.
- Native1D uses mixed polar/Cartesian midpoint coordinates, a zero-mean temporal gauge and physical-current/voltage/amplitude/native-energy temporal checks, with unchanged nonlinear/capacity/local/global energy gates. Exact gradient reuse must never round states.
- A0.5K collar rise is a visible warning; thermal truncation additionally requires significant native collar energy (default0.1% retained photon). Gap and integrated cut-energy guards remain independent.
- Late cold recovery is an additional physical approximation, never a preset-time substitution. Require full-domain admission and independent10ps circuit overlap. Preserve the actual circuit, derive cold inductance/energy/bath work from the native functional, and bound omitted normal Joule heat. Label reconstructed fields and do not claim late quasiparticle/retrigger validation. --no-cold-recovery is the full-native control.
- Pair commands use a3h1Dsession budget. WALL_BUDGET_CHECKPOINTED is resumable, not physical failure or completed15ns. Report actual warm/cold validation horizons and wall times; a fast cold control is not proof of a full-photon runtime.
- The old1.2.1 rigid-cut shadow is invalid for the new exterior condition. Recover authenticated2D/DC data into a new directory and replay only1D overlap; never redeposit a photon.
