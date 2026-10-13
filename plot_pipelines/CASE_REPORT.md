# Individual case report

Pipeline: /home/jdiaz/pysnspd/plot_pipelines/07_case_report.py

One self-contained PDF per effective photon run (pipeline05 schema).
Multiple input directories only batch independent reports; no inter-run comparison.
The report and all CLI/report text are in English. Current layout: **19 landscape A4 pages**.

## Sections

1. Case parameters: geometry, thickness versus mesh spacing, material, DC bias,
   deposition/retained energy, free Allmaras calibration parameters, circuit,
   discriminator and numerical settings.
2. Laboratory-related observables: baseline-subtracted Vout, readout current
   Ib-Is, confirmed crossing/latency, confirmation and final observation time.
   Circuit states, Vdev, coupling voltage and readout power are shown separately.
3. Internal fields: condensate magnitude, relative phase, equivalent electronic
   and phonon temperatures, excess energy, reconstructed super/normal/total
   current density, algebraic potential and **signed phase-winding vorticity**.
4. Independent energy balance.
5. Precision, provenance and limits.
6. Model summary: boxed coupled evolution equations, shared BCS energy,
   algebraic neutrality, circuit diagram, circuit equations and energy balance.
   Explains why the closure is called effective neutral BCS/KWT.

Vout has no filter, amplifier, noise or measurement-bandwidth model.
Only confirmed crossings define relative postcascade latency; censored cases
have no assigned latency. The observed maximum is not a complete pulse peak.
No full reset, jitter, statistical efficiency or absolute optical latency is inferred.

## Exactly six snapshot slots

--snapshots-ps accepts requested times in ps.

- Keep the first six entries if more are supplied.
- Repeat the last entry until there are six if fewer are supplied.
- Preserve input order and repetitions. Do not add initial/final states.
- Select the nearest actually saved state for each request; ties choose the earlier.
- A request outside the saved interval selects its nearest endpoint, without
  inventing data or extending the simulation.
- Display requested and saved time in every panel; retain both in the JSON receipt.

Examples: --snapshots-ps 0 1 gives [0,1,1,1,1,1] as target times.
--snapshots-ps 0 1 5 10 20 35 100 ignores 100.
Without the flag, six automatic targets cover the stored evolution and available
detection/hotbelt events. Repeated nearest states are retained even in short runs.

Every snapshot page remains a fixed 2-by-3 mosaic. Main axis labels are 14 pt,
ticks 12 pt; map axes 13 pt, map titles/ticks 12 pt, colorbar labels 14 pt.
Requested time is not substituted for saved time.

The receipt records cadence min/median/max, requested/saved offsets, repeated states
and ignored extra inputs. If too few states or too coarse a cadence prevent the
requested distinct times, future_visualization_queue records a follow-up to
save states with a stable, sufficiently fine physical-time cadence. No acquisition
or simulation source is changed by this visualization task. The two nominal
examples contain enough data for six distinct automatic states.

## Run several saved cases

Run on geminga. Choose a NEW output directory for each invocation.

~~~bash
cd /home/jdiaz/pysnspd
CASE_RUNS=/home/jdiaz/scratch/photon_runs/sensitivity_20261012_01/runs
env PYTHONPYCACHEPREFIX=/home/jdiaz/scratch/tmp/pysnspd/pycache \
  MPLCONFIGDIR=/home/jdiaz/scratch/tmp/pysnspd/matplotlib \
  NUMBA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
  MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 \
  /home/jdiaz/.conda/envs/snspd/bin/python \
  plot_pipelines/07_case_report.py \
  --run-dir "$CASE_RUNS/nominal775" "$CASE_RUNS/bias065_775" \
    "$CASE_RUNS/bias075_775" "$CASE_RUNS/nominal1550" \
  --snapshots-ps 0 1 5 10 20 35 \
  --output-dir /home/jdiaz/scratch/case_reports/english_20261013_01 \
  --workers 2
~~~

Omit --snapshots-ps for six automatic targets. One input directory also works.
--allow-incomplete captures one atomic checkpoint, without stopping or writing
to the running solver. Every page is provisional until a successful terminal
exists. A checkpoint without crossing is not terminal censoring at the ceiling.
The trajectory end and last saved snapshot are recorded separately.

Two independent workers use one physical core each, one-thread numerical
libraries and a shared resource budget capped at 90 %. No new evolution,
kinetic RHS evaluation or DC preparation is launched.

## Vorticity definition and limits

For each physical counterclockwise triangle C:

$$
n_C=\frac{1}{2\pi}\sum_{(i,j)\in C}
\operatorname{Arg}\!\left[\Delta_j\Delta_i^*\right].
$$

Positive/negative winding labels vortex/antivortex chirality. The core marker
is the zero of the same complex P1 interpolant, found by barycentric inversion.
This is phase topology, rather than curl(J): a current curl can arise from
inhomogeneous condensate density and is not itself a quantized vortex count.
The present states have A=0; this diagnostic does not reconstruct magnetic flux.

A contour is indeterminate (gray) if its minimum P1 edge amplitude is no larger
than 1e-6 of the minimum local DC vertex amplitude, or if its core cannot be
resolved consistently. This diagnostic resolution is not a physical model
parameter. White means resolved zero winding, not proof that no event occurred
between saved states. Short events, unresolved structures, edge cores and phase
in normal regions cannot be recovered as reliable vortex trajectories.

Primary reference: [Guo et al., Extracting, Tracking and Visualizing Vortices
in 3D Complex-Valued Superconductor Simulation Data, Eq. (1)](https://www.mcs.anl.gov/papers/P5324-0315.pdf).
The report caption links that source. Our positive sign explicitly denotes
counterclockwise physical phase winding; it is not a magnetic-flux sign claim.

## Quantity-specific palettes

| Quantity | Palette |
|---|---|
| Gap magnitude | viridis |
| Relative phase | twilight |
| Equivalent electronic temperature | inferno |
| Equivalent phonon temperature | YlOrBr |
| Signed excess energy density | PuOr_r |
| Signed super/normal current components | coolwarm |
| Total current magnitude | cividis |
| Algebraic potential | GnBu |
| Integer winding | discrete teal / white / magenta, gray indeterminate |

Color ranges are fixed across the six slots of a quantity. Energy and current
maps use different palettes and distinct signed physical units.

## Reconstruction and validation

Original terminal inputs, DC, mesh and reconstruction sources are authenticated.
Gap magnitude is evaluated after complex P1 spatial interpolation. Phase angles
are not interpolated. No temporal interpolation creates snapshots or events.

Supercurrent face fluxes use the exact gauge derivative of the pinned gradient,
independent of unsaved occupations. The same algebraic neutrality recovers
potential and normal fluxes at accepted states matching each snapshot exactly.
Vdev and current conservation are checked against saved data.

Current vectors by triangle are weighted least-squares display approximations
to conservative dual-face projections divided by physical face area.
Arrows show direction with a fixed 6 nm display length, omitting magnitudes
below 0.01 % of the selected maximum. Color shows signed x components or total
magnitude. Local signs do not establish reversal of the global bias.
Complete distributions, collision maps and KWT-heating maps are not invented.

Te/Tph are energy-equivalent observers. Saved quasiparticle energy excludes
condensation/gradient; the independent total energy balance includes them.

The 2026-10-13 refresh passes **31 quick tests**, including analytic signed
vortices/antivortices, global-phase/orientation invariance, zero/near-zero
contours, ordered six-slot selection, padding, truncation and source integrity.
The actual two-time CLI was exercised: six panels, 19 pages, repeated final
state. Two final examples take about six seconds per case; this is not a
guarantee for larger archives. Final PDFs and visual QA remain in
/home/jdiaz/scratch/case_report_refresh_20261013/release.
The 13 physical model sources remain unchanged.
