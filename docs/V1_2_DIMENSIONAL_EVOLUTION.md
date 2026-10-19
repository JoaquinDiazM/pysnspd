# pysnspd v1.2.0: dimensional evolution

The effective physical closure from v1.1.3 is unchanged. New code restricts its own functional, collision distributions, electroquasineutral solve and complete three-state circuit to the transverse constant mode. It does not import the Allmaras hotbelt constitutive equations.

- 300nm 2D strip; 0.1mV positive threshold,0.5ps confirmation,2ps total margin from crossing.
- Independent1D shadow during the final2ps, with finite accepted-state crossing lag recorded. All overlap checks and final restriction must pass.
- Conservative volume remap of the DC-relative complex gap and native electron/phonon populations. Copy circuit exactly. A bounded scalar activity shift reconciles excess energy between discrete DC references; no FD/BE reset or repeated photon.
- After admission, reproject the final official2D state and switch off2D entirely. The1D chain has61nodes over the same300nm. The exterior inductance is unchanged.15ns is an absolute final-time target, not a promise of complete reset.
- Integrate absolute left/right cut energy independently. Its2% activity-normalized margin is an explicit truncation indicator, not a theorem about the unobserved exterior. Significant cut flow stops either phase with POSSIBLE_DOMAIN_ERROR.
- Projection/independent-overlap errors or resolved vortex winding retain the native2D checkpoint and terminate POSSIBLE_DIMENSIONAL_ERROR. No failed criterion is silently waived.
- Coupled implicit midpoint is retained.1D adds full-step/two-half-step temporal error control, support/capacity/stiffness guards and the independent energy ledger with the inherited1% global envelope; no physical state is discarded for speed.
- Append-only phase-specific HDF5 stores native states and actual times. Atomic accepted-state checkpoints authenticate sources, parameters and the discrete DC.1D visualization cadence does not force a short integration step.
- A deterministic convex Voronoi reconstruction repairs the malformed300nm boundary cell while retaining Delaunay sites and triangles. Areas and common faces are recomputed geometrically, never globally renormalized.

Starting margins are engineering assumptions and CLI inputs: gap RMS3%, current-pattern residual3%, unsigned population energy5%, raw projection energy0.5%, reconciliation1%, shadow voltage2uV,current0.5%,film energy1%,cut activity2%. Reference is retained photon energy plus absolute device work except for temporal/local solver residuals, which keep the inherited photon-based normalization.

## Verification scope

Short regressions cover the exact uniform transverse operator, positive conservative remap, nonthermal occupation and relative-phase inheritance, winding of both signs, independent energy acceptance, archive commit/truncation, synthetic-probe control-flow handoff and exact resume/hash rejection, and shared-resource refusal. A synthetic voltage probe verifies wiring only, not photon detection.

Read-only projection of the saved775nm v1.1.3 crossing+5ps state (320nm) gave gap RMS0.0034852, population-energy error0.0046328 and raw projection-energy error0.00041949, with no resolved winding. A0.1ps65-node hot-state evolution on one physical core used80accepted steps,4rejections,29.25s and an independent energy defect9.34e-24J. This is not the300nm crossing+2ps admission test or a15ns benchmark. Production validation belongs to the new six-case campaign.

## Sources and assumptions

The derivation and tolerances are in the five stable Matematica PDFs, with editable sources in sandbox/model_v1_1_math_20261015. The former unimplemented spatial-tail proposal is replaced by the requested temporal whole-domain switch. No thermal tails are implemented in this release.

Primary context: [Allmaras et al.,1D hotbelt geometry](https://arxiv.org/pdf/1805.00130), [Allmaras2020 thesis](https://thesis.caltech.edu/13748/), [Vodolazov2017](https://arxiv.org/pdf/1611.06060), [Annunziata et al.,reset/latching](https://arxiv.org/abs/1008.0895). Their different closures, initial states and circuit assumptions are not silently substituted.

## Execution

Use docs/GEMINGA_COMMANDS_PRESENT.md. Six simultaneous cases use24Numba SMT threads on12physical worker cores plus a separate coordinator. Current v1.1.3 workers and their20pinned sources are preserved. Actual throughput and long-time dimensional/domain admissibility remain to be measured.
