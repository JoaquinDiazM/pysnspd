# Changelog

## 1.1.2 - 2026-10-16

- Prevent an accepted zero nonlinear predictor from advancing a nonstationary physical clock.
- Scale the thermal nonlinear error budget with adaptive step, include the full Te-induced gap response within the convergence loop, and backtrack internal corrections to retain positive capacity.
- Preserve all physical equations, calibrated diffusion, Te stiffness, signed Pauli heat and the independent energy guard; no physical capacity floor or tolerance relaxation.
- Set new thermal campaigns to100 ps without trigger and add an explicit checkpointed numerical-resolution stop below the declared adaptive limit.
- Preserve stalled1.1.1 checkpoints and record that the critical branch remains subject to fresh trajectory validation.


## 1.1.1 - 2026-10-15

- Fix a continuous-model Pauli violation in signed thermal heat: use FD birth/death activities selected by physical power sign, exact energy normalization and a nonnegative excess-entropy inequality.
- Preserve the Te-dependent energy ledger, depairing-calibrated diffusion, real circuit and positive-capacity guard; add no fitted timescale, clipping or heat from a residual.
- Add endpoint, subnormal-tail, signed entropy and full-RHS regressions; validate fresh 775/1550 nm cases through 0.2 ps on the nominal mesh.
- Replace five mathematical PDFs in place, align report equations/version and plotting DPI, and prepare a fresh anchor/control campaign. Old checkpoints remain historical; full photon admission remains pending.

## 1.1.0 - preparada, 2026-10-15

- Add an opt-in electron-temperature stiffness candidate with gradient free energy, internal energy/entropy, derived signed reversible heating and a positive-effective-capacity guard.
- Calibrate diffusion to the width-scaled Korzh depairing estimate, updating DOS and geometry normalization consistently; retain Rs425 and Rs608 scenarios.
- Replace the four mathematical PDFs and add a cross-referenced assumptions register.
- Prepare paired 775/1550 nm anchors and bounded mechanistic controls; extend per-case report reconstruction to the actual thermal closure.
- Preserve the published v1.0.0 solver and active-run source hashes. Candidate completion does not confer photon or experimental admission.

## 1.0.0 - 2026-10-15

- Freeze the conserving neutral effective BCS/KWT/native-kinetic/three-state-circuit implementation and its bounded photon pipelines.
- Include paired base-case sweep08, independent per-case reports07 and current Geminga command provenance.
- Preserve exploratory admission scope, numerical limitations, failed controls, post-cascade source hypotheses and scratch-only heavy data.
- Preserve the original thesis commit as thesis-v1.0.0 and explicitly reassign v1.0.0 at the user request.
- Reserve the next1.1.0 for temperature-dependent stiffness and depairing-calibrated diffusion.


## thesis-v1.0.0 - 2026-09-11

- Freeze the thesis-based solver and the latest previously published diagnostics.
- Publish model documentation 0.2, 0.3 and 0.4, including editable sources, PDFs,
  figures, independent lightweight calculations, data tables and verification records.
- Include the latest notebook E-r02 and its class/history registry; no unanswered
  activity is recorded as approved.
- Document the complete continuous candidate, boundaries, circuit, numerical
  diagnostics and the physical/mathematical restrictions on model 0.4.
- Set package version 1.0.0, declare runtime requirements and ship the plotting style.
- Record the five-minute computation policy and the Geminga command handoff.

Model 0.4 is documentation in this release. Its first experimental catalogue and
data-admission implementation belongs to subsequent commits. Historical transient
results are preserved, not rerun as part of the release.
