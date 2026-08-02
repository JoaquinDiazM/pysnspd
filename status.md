# pySNSPD publication status

Last updated: 2026-08-02

Publication window: 2026-07-23 to 2026-10-23

Current phase: Week 2 closure - stationarity, localized artifacts, and E1
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

## Validation completed

- Complete suite on Geminga after the E1/mesh update: `147 passed`.
- E1 was regenerated from the production PRE and its six affected PDFs were
  rendered and visually inspected without clipping, overlap, or illegible text.
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

1. Let the isothermal/no-circuit SS current sweep in screen `code3` finish.
2. Generate and visually validate its I-V curve and investigate any outlier
   before accepting the sweep.
3. Merge `codex/ss-photon-ready-stationarity` into `main` after that review.
4. Start the planned mesh/time/thermal convergence campaign using the new mesh
   edge-length and triangle-quality baselines.

The wider publication gates remain unchanged: mesh/time/thermal convergence,
time-resolved current continuity, accumulated energy closure, material
provenance, threshold sensitivity, and long-run latency/recovery acceptance
are still open.
