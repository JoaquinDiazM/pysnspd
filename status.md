# pySNSPD publication status

Last updated: 2026-07-30

Publication window: 2026-07-23 to 2026-10-23

Current phase: Week 2 — stationarity, convergence, and localized-artifact
regularization

## Executive status

pySNSPD is a functional multiscale research prototype coupling dirty-limit
Usadel material data, finite-volume gTDGL/Poisson dynamics, two-temperature
evolution, photon deposition, and circuit observables in physical units.

Week 2 is active. The stationarity diagnostics and the investigation of
localized low-amplitude artifacts are both advanced. The constitutive
regularization identified by the isolated diagnostic has now been implemented
through PRE, stationary, and photon execution paths. Engineering validation is
complete at smoke-test scale; long simulations have not yet been run with the
new catalogues, so neither the weekly scientific task nor the disappearance of
the 2D notches is considered resolved.

The SS stage also evolves the readout circuit together with the thermal state,
providing photon runs with a more faithful initial circuit condition. The
inspected 30 uA stationary result changed very little visually relative to the
previous thermal-only initialization, but a quantitative equivalence and
convergence claim is still pending.

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

This validates the numerical mechanism and motivates the correction, but does
not prove that it is the sole cause of every notch in the coupled 2D maps.

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

No plotting source or plotting pipeline was modified. Existing E2 and
single-run E3 diagnostics successfully consumed the new smoke-run archives, so
no plotting migration task is currently required.

## Validation completed

- Complete test suite on Geminga: `143 passed`.
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
- The unchanged E2 and E3 pipelines generated their complete diagnostic figure
  sets from those archives.

These runs validate interfaces, discrete identities, orchestration, and
serialization only. Their femtosecond horizons are not physical stationarity,
detection, recovery, or notch-removal evidence.

## Immediate work and acceptance gates

1. Run the prepared refined, parallel PRE command on Geminga.
2. Run the prepared classical 30 uA, 200 ps SS case from that PRE.
3. Run the prepared classical 0.8 eV photon case from the new stationary state.
4. Compare stationary and photon 2D fields against the previous baseline under
   identical numerical controls.
5. Close the artifact task only if the notches disappear or are quantitatively
   explained; close stationarity only after current, thermal, circuit, phase,
   and dynamic-attractor criteria pass.

The wider publication gates remain unchanged: mesh/time/thermal convergence,
time-resolved current continuity, accumulated energy closure, material
provenance, threshold sensitivity, and long-run latency/recovery acceptance
are still open.
