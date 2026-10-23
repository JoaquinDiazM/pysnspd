## 2026-10-23 | v1.2.4: four reflecting spatial probes and four seeded baseline1D continuations

Current publication-reference material: Rs608ohm/square, calibratedD0.6891369283
cm2/s, bath0.9K,Tc8.65K,width80nm,thickness7nm,sigma7.549834435nm,
tau_ee(Tc)5ps,tau_ep(Tc)24.7ps,escape20ps; all native quadratures, physical
tolerances and actual circuit retained. Evolution RHS unchanged. The new
entry19/producer18 identifies1.2.4; source-runner13 is the source-frozen
native kernel marker, NOT an instruction to choose the old release.

A. Four NEW LATENCY-ONLY cases (all3threads):
edge35_I15p5_1550:15.5uA,1550nm,retainedE0.71990825eV,y0=35nm.
edge35_I21p5_1550:21.5uA,1550nm,retainedE0.71990825eV,y0=35nm.
center_I18p5_1550:18.5uA,1550nm,retainedE0.71990825eV,y0=0nm.
center_I15p5_E096_775:15.5uA,775nm,retainedE0.96eV,y0=0nm.
First three retain eta0.9; explicit0.96eV corresponds to eta~0.60008.
Reflecting Neumann images at BOTH transverse edges; sum shapes BEFORE a single
native internal-energy inversion. Exact discrete total-energy normalization;
positive-capacity guards and unchanged initialDC gap. This is one photon and
an explicit postcascade sidewall hypothesis, not a measured reflection law.
At the boundary the coincident profile doubles the interior half-Gaussian,
without doubling total energy; substrate escape remains20ps.
Source: https://doi.org/10.1007/s00020-023-02749-0
Stop: confirmed+0.1mV crossing+2ps, confirmation0.5ps, no-trigger ceiling100ps.
No1D/reset evolution. From a separate shell open screen -S v124_photons,
then paste the following command:

```bash
cd /home/jdiaz/pysnspd
env PYTHONPYCACHEPREFIX=/home/jdiaz/scratch/tmp/pysnspd/pycache \
  MPLCONFIGDIR=/home/jdiaz/scratch/tmp/pysnspd/matplotlib \
  NUMBA_CACHE_DIR=/home/jdiaz/scratch/tmp/pysnspd/numba_cache \
  NUMBA_NUM_THREADS=3 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
  MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 \
  /home/jdiaz/.conda/envs/snspd/bin/python -u \
  pipelines/19_v124_campaign.py \
  --mode spatial \
  --campaign-dir /home/jdiaz/scratch/photon_runs/v124_spatial_probes_20261023_01 \
  --source-runner /home/jdiaz/pysnspd/pipelines/13_dimensional_photon_run.py \
  --execute --refresh-seconds 10
```

B. Four BASE Rs608 cases, full1D continuation, all3threads.
The immutable original2D data/source/DC/circuit were authenticated and copied
independently into the NEW recovery campaign. Only the native1D shadow was
replayed, starting at the last saved2D state <=trigger (<0.1ps earlier)
through trigger+2ps;22 physical overlap checks per case and all unchanged
handoff gates passed. Each advanced0.02ps into official1D during bounded
checks. The command RESUMES those already validated checkpoints; it does NOT
reintegrate2D, redeposit a photon, reset the bath or redo the overlap.
Current seeded1D times:
rs608_bias15p5_775 35.0690532541ps;
rs608_bias15p5_1550 48.6851464076ps;
rs608_bias21p5_775 22.5133684858ps;
rs608_bias21p5_1550 25.9654489148ps.
A root validated_seed_receipt.json records exact provenance/checksums.

Keep1000nm1D and15ns physical target;8hours NATIVE1D wall time per case/session,
excluding the short setup. The existing domain/current/voltage/capacity/
energy/recovery/cold-admission guards may stop earlier. A wall-budget stop
is resumable, not a physical failure or completed reset. From another shell
open screen -S v124_recovery, then paste:

```bash
cd /home/jdiaz/pysnspd
env PYTHONPYCACHEPREFIX=/home/jdiaz/scratch/tmp/pysnspd/pycache \
  MPLCONFIGDIR=/home/jdiaz/scratch/tmp/pysnspd/matplotlib \
  NUMBA_CACHE_DIR=/home/jdiaz/scratch/tmp/pysnspd/numba_cache \
  NUMBA_NUM_THREADS=3 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
  MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 \
  /home/jdiaz/.conda/envs/snspd/bin/python -u \
  pipelines/19_v124_campaign.py \
  --mode recovery \
  --campaign-dir /home/jdiaz/scratch/photon_runs/v124_base_recovery_20261023_01 \
  --base-campaign /home/jdiaz/scratch/photon_runs/calibration_eta090_ee5_rs425_rs608_20261021_01 \
  --source-runner /home/jdiaz/pysnspd/pipelines/13_dimensional_photon_run.py \
  --max-1d-wall-hours 8 \
  --execute --resume --refresh-seconds 10
```

Both screens can run together: each campaign has4x3Numba worker threads,
BLAS/OpenMP1,6physical worker cores with declared SMT sharing within pairs,
plus one separate coordinator core. Combined28/32logical CPUs admitted,
two COMPLETE cores reserved; conservative46GiB combined RAM reservation.
Runtime inventory/quotas/RAM/active workers and atomic live-PID/start-time CPU
leases are checked. Full startup masks of leased workers/coordinators are
accounted through their valid reservation. No existing process is repinned
or interrupted; insufficient capacity means refusal, not oversubscription.
Tested with two real coordinators and8synthetic workers starting concurrently;
assigned CPU sets were disjoint and the two-core reserve was preserved.
Allocation is not measured CPU utilization or a promised speedup.

Each screen shows phase,time,Vout,dt,accepted/rejected counts,rolling speed,
ETA and remaining1D wall budget, with CPU/RAM in wide terminals. Replay clock
is distinct from the saved final2D time. Complete logs are in campaign/logs/.
Ctrl-C/SIGTERM requests own checkpoints only. Add --resume to the first
command for its later checkpoint pauses; the second ALREADY includes --resume.
Own diagnostics are preserved and not blindly resumed.

Output: campaign plans/resource receipts/logs; each run's HDF5 archive,
manifest,source receipt/DC files and checkpoint. Nine event-based2D and
six longitudinal report slots remain available through pipeline08.
No statistical emulator or analytic latency fit is implemented.
Validation:91tests passed in178.51s; all4real reflected preparations +0.0025ps
passed; all4real overlap/handoffs+0.02ps1D passed in100-142s each without
relaxed tolerances. Exact-state wall checkpoint/resume passed;17-page1.2.4
report rendered and visually checked. All28native source hashes preserved.
Long production simulations have NOT been started by this preparation.
Detailed audit: /home/jdiaz/scratch/v124_validation_20261022/
Release notes: docs/releases/v1.2.4.md.

## 2026-10-22 | Rs425 individual reports: nine event-based 2D snapshots

Saved data only, current dimensional report pipeline08; no evolution, newDC,
model/version change, physical tolerance change or interruption of active jobs.
One English PDF and JSON receipt per case; nine2D panels in a3x3 mosaic with
a shared color range per quantity. The1D curve schedule remains SIX slots.

Without --snapshots-2d-ps, select independently WITHIN EACH CASE, in this
semantic order (do not sort or share event times across cases):
0ps,1ps,VS,VS+0.1ps,VS+0.2ps,HB,HB+2ps,Trigger,Final2D.
VS is the first SAVED resolved nonzero condensate-phase-winding P1 core;
uses the existing report detector, contour floor1e-6 relative to localDC.
It is not curl(J), and no sub-snapshot spawn time or vortex trajectory is inferred.
HB is the first SAVED connected transverse path of nodes with
abs(Delta)<=0.2*abs(DeltaDC), using the existing report definition.
Previous sample time and unresolved contour samples are recorded in the receipt;
short-lived vortices between stored frames may be missed.
Trigger is the first positive baseline-relative0.1mV crossing of accepted2D
Vout samples, linearly timestamped; controller-confirmed timestamp/detected
status remain separately recorded. Fields are NEVER interpolated in time.
Final2D is the last COMMITTED2D snapshot, not a later1D/reset state.

Missing VS/HB/Trigger: keep the requested semantic slots, repeat final2D and
label event not observed; requested_ps remains null. No fictitious onset.
Offsets outside2D window: use nearest endpoint, retain original target and
flag truncation. Requests/saved times are shown on every panel. Repeated
frames remain visible. Exact nearest ties choose earlier saved sample.
Manual --snapshots-2d-ps (alias --snapshots-ps): first NINE inputs; ignore
extras, repeat last input if fewer than nine; finite nonnegative used inputs
required. --snapshots-1d-ps retains first SIX/pad-last behavior.
Receipt schema v4 stores each event definition/time/provenance, nine slots,
fallbacks/offsets/truncation/repeats, cadence and future visualization queue.

Run the following command once ALL FOUR inputs are terminal. Batch preflight
checks every manifest before creating any report, so an active case causes a
clear refusal without partially generating the batch. At preparation readback,
775nm pair was complete, whereas both1550nm Rs425 processes were still active.
They were preserved. Rerun the same command after completion; don't use an
already populated report destination from an earlier successful batch.
Two report workers, one numerical-library thread each; no full-reset run.
All generated outputs/caches/build files stay in scratch.

```bash
cd /home/jdiaz/pysnspd
env PYTHONPYCACHEPREFIX=/home/jdiaz/scratch/tmp/pysnspd/pycache \
  MPLCONFIGDIR=/home/jdiaz/scratch/tmp/pysnspd/matplotlib \
  NUMBA_CACHE_DIR=/home/jdiaz/scratch/tmp/pysnspd/numba_cache \
  NUMBA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
  MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 \
  /home/jdiaz/.conda/envs/snspd/bin/python \
  plot_pipelines/08_dimensional_case_report.py \
  --run-dir \
    /home/jdiaz/scratch/photon_runs/calibration_eta090_ee5_rs425_rs608_20261021_01/runs/rs425_bias15p5_775 \
    /home/jdiaz/scratch/photon_runs/calibration_eta090_ee5_rs425_rs608_20261021_01/runs/rs425_bias15p5_1550 \
    /home/jdiaz/scratch/photon_runs/calibration_eta090_ee5_rs425_rs608_20261021_01/runs/rs425_bias21p5_775 \
    /home/jdiaz/scratch/photon_runs/calibration_eta090_ee5_rs425_rs608_20261021_01/runs/rs425_bias21p5_1550 \
  --output-dir /home/jdiaz/scratch/case_reports/rs425_events9_20261022_01 \
  --workers 2
```

Validation:39tests passed in0.68s, including actual phase-winding/hotbelt
fixtures, per-case defaults, missing events, uncommitted data exclusion,
explicit slots, endpoint/tie rules, voltage polarity and live-batch preflight.
One real completed Rs42521.5uA/775nm report generated13pages in7.14s;
nine distinct2D slots. VS6.8ps,HB7.6ps,Trigger22.0220092991ps,
Final2D24.0220092991ps. Gap/vorticity3x3 pages rendered and visually checked.
Rs42515.5uA/775nm readback independently selectedVS14.9ps,HB15.1ps.
The29 campaign source hashes, completed manifest/checkpoint hashes remain
unchanged. No production physics or new simulation was executed.
Audit: /home/jdiaz/scratch/report_nine_20261022/validation_receipt.json

## 2026-10-22 | Four Allmaras/Korzh historical-parameter controls alongside the active Rs425 campaign

Purpose: answer whether the current native2D effective neutral BCS/KWT model has
been tried with the historical fitted hotbelt inputs. Four deterministic cases:
15.5/21.5uA x775/1550nm; width80nm, thickness7nm, Tc8.65K.
Rs608ohm/square; fixed D0.5cm^2/s (=5e-5m^2/s), gamma60, eta_ret0.667,
tau_ee(Tc)6ps, tau_ep(Tc)24.7ps, tau_escape20ps.
Initial uniform transverse hotbelt40nm long, common FD/BE energy preparation,
native current-carrying DC gap. Bath2K is the historical simulated bath
(thesispp95-98), not the0.9K experiment/current calibration-campaign bath.
Primary parameter source: Korzh2020 Supplementary Note3/Table1, PDF33:
https://eprints.lancs.ac.uk/id/eprint/140252/3/Binder1.pdf#page=33
Thesis: https://thesis.caltech.edu/13748/08/Allmaras_Thesis_Final.pdf#page=107

Explicit scope: current v1.2.3 equations, kinetics, circuit and tolerances;
NOT an implementation of the old Allmaras1D solver or reproduction of its IRF.
The published Fano92meV and nonuniformity40meV widths are recorded but not
sampled in these four mean-energy cases. Gamma60 independently rescales the
phonon DOS/capacity consistently within our native kernel; nominal sound
speed alone no longer determines its normalization. D is fixed, NOT retuned
to the experimental inferred depairing current. Native predicted bulk Idep
is about29.4282uA at2K;21.5/15.5uA are about0.731/0.527 of that prediction.
These are declared historical controls, not a new physical parameter admission.

Observation: positive Vout0.1mV crossing;0.5ps confirmation; extra2ps;
100ps no-trigger censoring. Latency-only, no1D evolution/full reset. Mesh300nm
long,80nm wide, target spacing5nm. Retain native circuit/detection definitions.

When the last Rs608 process exits, open a NEW screen from a separate shell:
screen -S allmaras2020
Then execute the command below. It does not control, interrupt or repin the
existing screen/campaign. Resource checks run again on launch. With four Rs425
workers present:12new Numba worker threads (4x3),13existing worker threads,
two lightweight coordinators sharing the admitted coordination core;27/32
logical CPUs admitted (84.375%), at least two complete physical cores reserved.
New numerical-worker masks are disjoint from existing ones. Some SMT siblings
are shared only between NEW cases. BLAS/OpenMP1; conservative23GiB additional
RAM reservation, checked against90%currently available RAM.
This allocation is not a promise of measured84%CPU utilization or speedup.
If Rs608/unknown jobs still occupy the needed cores, the command refuses
without starting workers or affecting those jobs. No automatic waiting/polling.
A separate resource lease is used; the running campaign's exclusive lock stays
untouched. New coordinator/worker argv explicitly expose the native producer
to the unchanged old resource scanner, so a future overlapping eight-case
launch is refused while these historical controls run.

```bash
cd /home/jdiaz/pysnspd
env PYTHONPYCACHEPREFIX=/home/jdiaz/scratch/tmp/pysnspd/pycache \
  MPLCONFIGDIR=/home/jdiaz/scratch/tmp/pysnspd/matplotlib \
  NUMBA_CACHE_DIR=/home/jdiaz/scratch/tmp/pysnspd/numba_cache \
  NUMBA_NUM_THREADS=3 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
  MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 \
  /home/jdiaz/.conda/envs/snspd/bin/python -u \
  pipelines/17_allmaras_reference_cases.py \
  --campaign-dir /home/jdiaz/scratch/photon_runs/allmaras2020_parameter_control_20261022_01 \
  --coexist-with-campaign /home/jdiaz/scratch/photon_runs/calibration_eta090_ee5_rs425_rs608_20261021_01 \
  --source-runner /home/jdiaz/pysnspd/pipelines/13_dimensional_photon_run.py \
  --execute --refresh-seconds 10
```

Terminal: four-row English progress board every10s; simulated time,Vout,
accepted/rejected steps, rolling ps/h, CPU,RAM and ETA where width permits.
Before a trigger ETA refers to100ps ceiling, NOT a detection forecast.
Expected wall time: multi-hour and parameter-dependent; no completed-run
runtime forecast for this historical set yet. Full raw logs stay in logs/.
Outputs: campaign_plan.json, execution_resources.json, campaign_execution.json;
runs/allmaras2020_bias{15p5,21p5}_{775,1550}/ with data.h5, manifest.json,
checkpoint_latest.npz, source_initial_state.npz, DC/material/source receipts.
Some1D setup metadata remains for archive compatibility; no1D evolution runs.
Run in a screen to survive SSH disconnect. Ctrl-C/SIGTERM affects ONLY these
new producers and requests own checkpoints. Resume the SAME command with
--resume appended; native/adapter source hashes and exact inputs must match.
Completed own cases are reused. Diagnostic failures remain preserved and are
not blindly resumed. Never use a pre-existing calibration campaign as output.

Preparation validation:5tests passed; real native3260-node DC+photon+0.001ps
probes at21.5uA/775nm and15.5uA/1550nm passed in44.15/35.72s. Source-energy
relative defects7.94e-17/1.05e-16; positive capacity guards retained; deposition
does not touch cuts. The29 live-campaign source hashes remain unchanged.
Probe runtime excludes prediction of latency-run wall time; no production
historical case was started. Terminal board and producer/plan source identity
were verified. Audit receipt:
 /home/jdiaz/scratch/allmaras_reference_preparation_20261022/validation_receipt.json


# Current commands: pysnspd v1.2.3


## 2026-10-21 | Eight-case calibration: eta0.9, tau_ee5, Rs425 vs Rs608

Purpose: one common material per four-case set, each crossing 15.5/21.5 uA with
775/1550 nm, all eta_ret0.9, tau_ee(Tc)5ps, tau_ep(Tc)24.7ps, escape20ps,
80nm width, 7nm thickness, 300nm2D, sigma7.549834435nm, Tb0.9K.
Use the latest v1.2.3 native physical kernel. Optional latency-only execution
stops at confirmed positive0.1mV crossing +2ps (confirmation0.5ps); no trigger
at100ps is censored. No shadow1D or official1D time evolution/full reset.

Control rationale: Rs425 is conditional340/RRR0.8 inference; Rs608 is the
alternative used in K20 Supplementary Table1, not a measured resistance of our
particular nanowire. Both recalibrate D against the same width-scaled depairing
target25.8666667uA. D425=1.4103757 and D608=0.6891369 cm2/s.
N0, ell0, physical mesh and frozen exterior inductance are derived consistently.
Gamma is not independently fitted; for the acoustic model gamma~11.907/17.034.
Lower D improves localization, while larger N0 changes heat capacity and
initial energy partition; the net latency effect is a test, not a promised result.
Primary sources:
https://eprints.lancs.ac.uk/id/eprint/140252/3/Binder1.pdf#page=33
https://www.nature.com/articles/s41598-018-20035-7 (other NbN films, not this sample).

Resources:16physical/32logical CPUs, no CPU quota in current inventory.
26Numba threads: six cases x3 plus two low-current1550nm cases x4.
Workers81.25% logical capacity; coordinator on a separate complete core gives
28/32 admitted logical CPUs (87.5%). Reserve two whole physical cores.
Three pairs of cases share SMT siblings on one core, with disjoint logical
CPU masks; this is explicitly not whole-core exclusivity per case.
BLAS/OpenMP1 per process, eight simultaneous workers, conservative43GiB RAM
reservation (<90% available). Runtime rechecks topology/quotas/RAM/active
simulation processes and refuses overlapping campaigns or insufficient capacity.
Allocation percentages are not measured CPU utilization or guaranteed speedup.

```bash
cd /home/jdiaz/pysnspd
env PYTHONPYCACHEPREFIX=/home/jdiaz/scratch/tmp/pysnspd/pycache \
  MPLCONFIGDIR=/home/jdiaz/scratch/tmp/pysnspd/matplotlib \
  NUMBA_CACHE_DIR=/home/jdiaz/scratch/tmp/pysnspd/numba_cache \
  OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 \
  /home/jdiaz/.conda/envs/snspd/bin/python -u \
  pipelines/15_calibration_basecases.py \
  --campaign-dir /home/jdiaz/scratch/photon_runs/calibration_eta090_ee5_rs425_rs608_20261021_01 \
  --execute --refresh-seconds 10
```

The coordinator redraws an English table every10s in a TTY (screen supported);
narrow terminals use a compact table. It shows DC/source preparation, simulated
time, Vout, accepted step, accepted/rejected counts, rolling ps/h, ETA, actual
process CPU% (100%=one logical CPU), and RSS when width permits.
Before crossing, ETA is to the100ps search ceiling, not predicted detection;
after crossing, ETA is to crossing+2ps. --plain prints periodic tables without
ANSI control. Complete raw producer JSON/tracebacks remain in campaign/logs/.
At exit it prints the four1550-minus775 relative latencies when available.

Output: campaign_plan.json, execution_resources.json, campaign_execution.json,
logs/<case>.log and runs/<case>/{data.h5,manifest.json,checkpoint_latest.npz,
setup2d/,dc1d.npz,source_initial_state.npz}. 1D setup metadata is retained for
archive/report compatibility, but no1D evolution is requested.
Ctrl-C/SIGTERM requests checkpoints from all active producers; for a fully
checkpointed interruption, rerun the same command with --resume added. Failed
diagnostic states are preserved and not automatically resumed.
Exact source identities are checked; no old data or checkpoints are overwritten.
New --latency-only and --cpu-list options change observation/scheduling only;
the24 audited physical evolution source hashes and all physical tolerances
are unchanged. Default dimensional execution still follows the existing1D path.
Preparation performed here does not launch any production case.



## 2026-10-21 | Saved-data reports: two snapshot schedules and tabulated design

Both v1.2.3 cases stopped intentionally at the three-hour 1D session budget, with resumable checkpoints: 775 nm at 338.068921 ps; 1550 nm at 304.909834 ps. This command only renders saved data; it neither resumes a simulation nor computes a new DC state. Two workers, one numerical-library thread per worker. Original manifest/checkpoint hashes remain unchanged.

```bash
cd /home/jdiaz/pysnspd
env PYTHONPYCACHEPREFIX=/home/jdiaz/scratch/tmp/pysnspd/pycache \
  MPLCONFIGDIR=/home/jdiaz/scratch/tmp/pysnspd/matplotlib \
  NUMBA_CACHE_DIR=/home/jdiaz/scratch/tmp/pysnspd/numba_cache \
  NUMBA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
  MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 \
  /home/jdiaz/.conda/envs/snspd/bin/python \
  plot_pipelines/08_dimensional_case_report.py \
  --run-dir \
    /home/jdiaz/scratch/photon_runs/dimensional_v1_2_3_recovered_20261021_01/runs/bias21p5_ee4_775 \
    /home/jdiaz/scratch/photon_runs/dimensional_v1_2_3_recovered_20261021_01/runs/bias21p5_ee4_1550 \
  --snapshots-2d-ps 1 11 20 23 25 27 \
  --snapshots-1d-ps 1 20 50 100 200 340 \
  --output-dir /home/jdiaz/scratch/case_reports/dimensional_v1_2_3_design_20261021_01 \
  --workers 2
```

Each report has 17 English pages. --snapshots-2d-ps selects six nearest saved colormaps (--snapshots-ps remains an alias). --snapshots-1d-ps selects six longitudinal profiles: pre-handoff slots are volume-conserving scalar averages of saved 2D fields over the actual 300 nm lens, with no exterior filling; later slots use native 1000 nm 1D data. Both flags retain the first six inputs or repeat the last to reach six; no temporal interpolation. Requested 340 ps uses the last available state. Repeated colormap frames near handoff are intentional and explicitly labelled.

No model evolution, tolerances, version or release changed. The cold constitutive/circuit reduction remains inactive in both cases; a new hot 0D closure has not been implemented. Internal render QA and source-integrity receipt: /home/jdiaz/scratch/report_design_20261021_01/report_plot_validation.json.


## 2026-10-21 | v1.2.3 exact-state continuation after monitor audit

Purpose: keep every valid2D/1D accepted state and continue the stopped v1.2.2
775/1550nm pair from62.2983964514/54.8984934464ps. Same21.5uA, escape20ps,
eta0.8, Rs425ohm, tau_ee(Tc)4ps,300nm2D/1000nm1D,4threads/2exclusive physical
cores per case. No photon redeposition, freshDC, shadow replay or bath reset.
The new source directory must be the v1.2.2 campaign below.

Diagnosis: unloading increases the exterior gap; Te rise outside the hotbelt
can include transport and KWT/normal dissipation. The old collar criterion
used photon energy alone in a strongly driven circuit. New active native
collar normalization uses retained photon plus cumulative absolute device
work (same activity ledger as cuts), retaining photon-only values visibly.
Gap suppression uses instantaneous-current cold bulk rather than initialDC.
Margins0.1%,90%,2ps and all physical/numerical/cut/electrical guards remain.
This is a declared engineering diagnostic revision, not an RHS repair or
proof of finite-domain convergence. See docs/V1_2_DIMENSIONAL_EVOLUTION.md.

Reauthentication requires exact published v1.2.2 source pins, unchanged RHS,
identical inputs, checkpoint/DC/archive hashes and NEW monitor readmission
before stepping. Original folders/failure records remain unchanged.
Validated:27tests, both saved photon states +3.2ps, and .04/.02ps refinement.
No full15ns runtime/temperature/domain admission is inferred.

Resources:8worker SMT threads/4whole physical cores for this pair, plus
separate coordinator; at least2whole cores reserved; conservative11GiB RAM.
All heavy data/caches/logs in scratch.3h1Dsession cap saves a resumable
WALL_BUDGET_CHECKPOINTED. Add --resume to continue the NEW campaign later;
add --prepare-only to inspect its plan without simulations.

```bash
cd /home/jdiaz/pysnspd
env PYTHONPYCACHEPREFIX=/home/jdiaz/scratch/tmp/pysnspd/pycache \
  MPLCONFIGDIR=/home/jdiaz/scratch/tmp/pysnspd/matplotlib \
  NUMBA_CACHE_DIR=/home/jdiaz/scratch/tmp/pysnspd/numba_cache \
  NUMBA_NUM_THREADS=4 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
  MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 \
  /home/jdiaz/.conda/envs/snspd/bin/python -u \
  pipelines/14_dimensional_basecases.py \
  --campaign-dir /home/jdiaz/scratch/photon_runs/dimensional_v1_2_3_recovered_20261021_01 \
  --recover-from-campaign /home/jdiaz/scratch/photon_runs/dimensional_v1_2_2_recovered_20261021_01 \
  --cases bias21p5_ee4_775 bias21p5_ee4_1550 \
  --length-1d-nm 1000 \
  --threads-per-case 4 --physical-cores-per-case 2 \
  --ceiling-ps 100 --final-time-ns 15 \
  --max-1d-wall-hours 3
```

# Current report command: stopped v1.2.2 pair

## 2026-10-21 | Report-only diagnostic template for both stopped v1.2.2 cases

Saved data only; no simulation, no kinetic RHS, no DC solve and no model/version update.
Use phase-aware pipeline08 for 2D + official1D. One independent English PDF per
case,16pages with observables, phase-relative energy/power, numerical resolution,
six2D gap/Te/Tph/energy/total-current/winding mosaics,1D space-time fields and
longitudinal profiles, separate free left/right collar temperatures/gaps/native
energy with the actual stop thresholds, last native spectra, independent overlap
and circuit/model scope. Remove normal and supercurrent snapshot pages from
both templates07/08; retain total current. Zero-area Voronoi faces are excluded
only from display-vector fitting, not altered in the simulation mesh.

Default2D slots follow the saved detection/handoff times; default1D profiles
span handoff through the terminal state. --snapshots-ps sets six2D slots;
--snapshots-1d-ps independently sets six1D profile slots. Each selects nearest
stored states; extra requests are ignored and the last request is repeated when
fewer than six are supplied. No temporal interpolation.
Collar curves use stored snapshots; the report does not recreate accepted-step
persistence. A JSON receipt stores exact chosen frames and reconstruction checks.
No amplifier/filter is added. Outputs do not overwrite prior reports.

Resources:2workers,1Numba and1BLAS thread each; no numerical evolution.
Internal validation: both actual cases rendered successfully,16pages each;
last reconstructed collar metrics agree with the simulator; visual QA performed.
The original simulation sources/checkpoints/archives remain unchanged.

```bash
cd /home/jdiaz/pysnspd
env PYTHONPYCACHEPREFIX=/home/jdiaz/scratch/tmp/pysnspd/pycache \
  MPLCONFIGDIR=/home/jdiaz/scratch/tmp/pysnspd/matplotlib \
  NUMBA_CACHE_DIR=/home/jdiaz/scratch/tmp/pysnspd/numba_cache \
  NUMBA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
  MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 \
  /home/jdiaz/.conda/envs/snspd/bin/python -u \
  plot_pipelines/08_dimensional_case_report.py \
  --run-dir \
    /home/jdiaz/scratch/photon_runs/dimensional_v1_2_2_recovered_20261021_01/runs/bias21p5_ee4_775 \
    /home/jdiaz/scratch/photon_runs/dimensional_v1_2_2_recovered_20261021_01/runs/bias21p5_ee4_1550 \
  --output-dir /home/jdiaz/scratch/case_reports/dimensional_v1_2_2_diagnostics_20261021_01 \
  --workers 2
```

# Current commands: pysnspd v1.2.2


## 2026-10-21 | v1.2.2: cold cuts, native polar solver and gated late recovery

Purpose: reuse valid2D775/1550nm at21.5uA, escape20ps, eta0.8, Rs425ohm,
tau_ee(Tc)4ps. Keep300nm2D/1000nm1D. Replay only1D overlap with current-following
cold cuts; the old fixed-cut1D heating is not a valid new exterior history.
Preserve every original directory and never redeposit the photon.

Hot1D remains native. Polar coordinates, zero-mean temporal gauge and exact
gradient reuse accelerate it. A late cold approximation additionally requires
full-domain gap/native-energy/current/voltage bounds AND10ps independent overlap.
It is never forced at100ps. The actual3-state circuit, cold differential
inductance, bath/energy accounting and omitted-Joule budget remain explicit.
Read docs/V1_2_DIMENSIONAL_EVOLUTION.md. --no-cold-recovery selects a full-native
control in a fresh campaign directory.

Resources:4SMT threads/2exclusive physical cores per case;8threads/4worker cores
plus coordinator for this pair; BLAS/OpenMP1;>=2whole cores reserved,11GiB
conservative RAM reservation. Heavy files/logs/caches stay in scratch.
The3h1Dsession budget saves WALL_BUDGET_CHECKPOINTED rather than claiming
physical completion. Same command + --resume continues the NEW campaign.

Validated scopes: actual checkpoints beyond the old2ps stop; gauge invariance,
boundary work, exact cache, cold admission/restart/labels and guards; a
nontrivial15ns cold-circuit control in about4s after preparation. Full15ns
photon recovery and<3h end-to-end remain pending. COLD_RECOVERY_ADMITTED
identifies entry into the fast regime.

Simulation (user-run):
```bash
cd /home/jdiaz/pysnspd
env PYTHONPYCACHEPREFIX=/home/jdiaz/scratch/tmp/pysnspd/pycache \
  MPLCONFIGDIR=/home/jdiaz/scratch/tmp/pysnspd/matplotlib \
  NUMBA_CACHE_DIR=/home/jdiaz/scratch/tmp/pysnspd/numba_cache \
  NUMBA_NUM_THREADS=4 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
  MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 \
  /home/jdiaz/.conda/envs/snspd/bin/python -u \
  pipelines/14_dimensional_basecases.py \
  --campaign-dir /home/jdiaz/scratch/photon_runs/dimensional_v1_2_2_recovered_20261021_01 \
  --recover-from-campaign /home/jdiaz/scratch/photon_runs/dimensional_v1_2_1_pair_20261019_01 \
  --cases bias21p5_ee4_775 bias21p5_ee4_1550 \
  --length-1d-nm 1000 \
  --threads-per-case 4 --physical-cores-per-case 2 \
  --ceiling-ps 100 --final-time-ns 15 \
  --max-1d-wall-hours 3
```

Add --prepare-only for the exact plan/capacity check without simulation.
The100ps no-trigger ceiling does not truncate admitted1D evolution.

Individual phase-aware visualization:
```bash
cd /home/jdiaz/pysnspd
env PYTHONPYCACHEPREFIX=/home/jdiaz/scratch/tmp/pysnspd/pycache \
  MPLCONFIGDIR=/home/jdiaz/scratch/tmp/pysnspd/matplotlib \
  OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 \
  /home/jdiaz/.conda/envs/snspd/bin/python \
  plot_pipelines/08_dimensional_case_report.py \
  --run-dir /home/jdiaz/scratch/photon_runs/dimensional_v1_2_2_recovered_20261021_01/runs/* \
  --snapshots-ps 1 11 20 25 30 40 \
  --output-dir /home/jdiaz/scratch/case_reports/dimensional_v1_2_2_recovered_20261021_01 \
  --workers 2
```

# Current commands: pysnspd v1.2.1 handoff correction


## 2026-10-20 | v1.2.1 handoff repair and authenticated checkpoint reuse

Purpose: continue the 775/1550nm pair at21.5uA from its valid pre-handoff2D checkpoints (23.0809897591ps /26.9709634588ps). Preserve original run directories, DC states, trajectories, energy ledgers and native shadow histories. Never redeposit the photon or re-run the valid2D interval.

Repair: reconstruct the nodal complex gap at physical sites rather than treating FV partial-cell averages as site values; match each disconnected exterior phase at its own interface, selecting interface nodes by index to avoid floating-coordinate exclusion. Keep native p,n volume restriction and small bounded energy reconciliation. Energy and projection gates remain0.5% raw/1% reconciliation/3% complex projected gap. Independent overlap uses3% gap-amplitude and3% actual supercurrent norms, with full complex bank-phase discrepancy explicitly retained as a diagnostic. The actual phase-dependent current, native populations, energy, circuit and vortex checks remain active. No material parameter, circuit or kinetic kernel changed.

The continuation authenticates the known original runner/projection hashes, every unchanged physical-source hash, parameters, checkpoint andDC. It copies and hashes the archive; the originals remain untouched. Stored overlap is re-audited using exact2D snapshots and BOTH bracketing native1D snapshots. Work/energy and circuit errors use stored simultaneous comparisons, never different-time energies. No missing state history is invented. Rejection remains a diagnostic outcome.

Resources: same4Numba SMT threads on2exclusive physical cores per case (8threads/4worker cores plus coordinator); two full cores reserved, conservative11GiB RAM, all heavy files/caches under scratch. Global lock/budget remain. Live progress now reaches the terminal as well as per-case logs. This command starts no old1.1.3 cases. Full15ns validation and throughput remain pending.

Outputs: fresh scratch campaign, full reused2D archive and new1D evolution; continuation provenance, pre_update_manifest.json and pre_update_checkpoint.npz, native HDF5/checkpoints, projection/overlap/energy/monitor records. Use the new directory below. Add --prepare-only to check the plan without starting simulations. For an interrupted NEW continuation campaign, add --resume to this same command; the original-source flag is retained in its plan and no second copy occurs.

```bash
cd /home/jdiaz/pysnspd
env PYTHONPYCACHEPREFIX=/home/jdiaz/scratch/tmp/pysnspd/pycache \
  MPLCONFIGDIR=/home/jdiaz/scratch/tmp/pysnspd/matplotlib \
  NUMBA_CACHE_DIR=/home/jdiaz/scratch/tmp/pysnspd/numba_cache \
  NUMBA_NUM_THREADS=4 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
  MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 \
  /home/jdiaz/.conda/envs/snspd/bin/python -u \
  pipelines/14_dimensional_basecases.py \
  --campaign-dir /home/jdiaz/scratch/photon_runs/dimensional_v1_2_1_recovered_20261020_01 \
  --recover-from-campaign /home/jdiaz/scratch/photon_runs/dimensional_v1_2_1_pair_20261019_01 \
  --cases bias21p5_ee4_775 bias21p5_ee4_1550 \
  --length-1d-nm 1000 \
  --threads-per-case 4 --physical-cores-per-case 2 \
  --ceiling-ps 100 --final-time-ns 15
```

# Current commands: pysnspd v1.2.1

## 2026-10-19 | v1.2.1: 300nm 2D to 1000nm full1D

Purpose: test the longer post-detection domain, native tail history, correct film/circuit inductance transfer and explicit physical stop monitors. Six simultaneous cases: currents15.5uA with tau_ee(Tc)5/4ps, and21.5uA with4ps, each775/1550nm. Shared escape20ps,eta.8,Rs425ohm,width80nm. Same physical BCS/KWT/native-kinetic closure.

The2D lens stays300nm; the full1D default is1000nm/201nodes,350nm added per side. First crossing starts an independent shadow with a cold instantaneous-current exterior hypothesis. At crossing+2ps (0.5ps confirmation included), retain official2D core data and evolved native shadow tails. Do not reinitialize tails cold at handoff. Compare gap modulo a global phase because potential references use different right cuts.

Freeze both inductance partitions from the original DC, preserving96nH. Transfer the removed lumped excess energy into the resolved film, including full circuit energy in the reference shift. No double count. A1um domain is a starting hypothesis, not a guarantee through15ns.

Stops: 2% activity-normalized cut exchange; persistent outer25nm collar heating/suppression (0.5K or gap<90%DC for2ps); global current/voltage limits; support/capacity/stiffness/energy guards. Possible latching means>=20nm resistive segment and1ns electrical/thermal plateau after5ns. Early recovery requires all declared gap/native/electron+phonon/full-circuit criteria for1ns after2ns. Moderate negative Vout is allowed as AC undershoot; global branch reversal is guarded. Check manifest monitor_config/monitor_reason; coordinator prints explicit reasons.

Resources:24Numba SMT threads on12physical worker cores plus coordinator, at least2whole cores reserved; BLAS/OpenMP1; conservative27GiB RAM. Runtime lock preserves existing jobs. Wait until the current v1.1.3 campaign finishes; no automatic wait/service. Wall-time to15ns is not benchmarked, and CPU70% is not guaranteed.

Outputs: scratch campaign logs/resources, phase-aware HDF5 including native p,n, state checkpoints, source/DC/material/mesh identity, inductance_partition and monitor history/reasons. Never resume1.2.0 under1.2.1. Use a fresh directory. --length-1d-nm 2000 is an optional future domain control, with a different campaign path; do not silently change a running case.

Execution (user-run, foreground):
```bash
cd /home/jdiaz/pysnspd
env PYTHONPYCACHEPREFIX=/home/jdiaz/scratch/tmp/pysnspd/pycache \
  MPLCONFIGDIR=/home/jdiaz/scratch/tmp/pysnspd/matplotlib \
  NUMBA_CACHE_DIR=/home/jdiaz/scratch/tmp/pysnspd/numba_cache \
  NUMBA_NUM_THREADS=4 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
  MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 \
  /home/jdiaz/.conda/envs/snspd/bin/python \
  pipelines/14_dimensional_basecases.py \
  --campaign-dir /home/jdiaz/scratch/photon_runs/dimensional_v1_2_1_20261019_01 \
  --length-1d-nm 1000 \
  --threads-per-case 4 --physical-cores-per-case 2 \
  --ceiling-ps 100 --final-time-ns 15
```

Append --prepare-only for plan/capacity checking without simulation. Append --resume only for an interrupted identical1.2.1 case/campaign.

Individual phase-aware visualization (intrinsic Vout without amplifier/filter):
```bash
cd /home/jdiaz/pysnspd
env PYTHONPYCACHEPREFIX=/home/jdiaz/scratch/tmp/pysnspd/pycache \
  MPLCONFIGDIR=/home/jdiaz/scratch/tmp/pysnspd/matplotlib \
  OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 \
  /home/jdiaz/.conda/envs/snspd/bin/python \
  plot_pipelines/08_dimensional_case_report.py \
  --run-dir /home/jdiaz/scratch/photon_runs/dimensional_v1_2_1_20261019_01/runs/* \
  --snapshots-ps 1 11 20 30 50 70 \
  --output-dir /home/jdiaz/scratch/case_reports/dimensional_v1_2_1_20261019_01 \
  --workers 2
```


## 2026-10-19 | Two priority v1.2.1 cases, exclusive cores

Selected 775/1550nm at21.5uA, tau_ee(Tc)4ps, escape20ps, eta0.8, Rs425ohm; 300nm2D/1000nm1D, 0.5ps confirmation included in2ps extra, no-trigger ceiling100ps, guarded1D target15ns. Same4Numba SMT threads /2physical cores per case, eight threads/four worker cores total plus a separate coordinator; reserve two complete cores and estimate11GiB RAM. Outputs: per-case logs, manifests, phase-aware HDF5/checkpoints and resources under the scratch campaign.

The user chose exclusive cores and waiting over shared SMT. The existing global campaign lock still refuses concurrent campaigns, even if individual workers finish first. Run after the old coordinator exits. This command does not wait automatically or interrupt current jobs. --prepare-only checks planning/capacity without a simulation. Added --cases selection changes orchestration only; physical source hashes and v1.2.1 equations are unchanged. Short allocation checks pass for two and six cases, partial occupancy and busy refusal; no heavy run launched. Wall-time unmeasured.

```bash
cd /home/jdiaz/pysnspd
env PYTHONPYCACHEPREFIX=/home/jdiaz/scratch/tmp/pysnspd/pycache \
  MPLCONFIGDIR=/home/jdiaz/scratch/tmp/pysnspd/matplotlib \
  NUMBA_CACHE_DIR=/home/jdiaz/scratch/tmp/pysnspd/numba_cache \
  NUMBA_NUM_THREADS=4 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
  MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 \
  /home/jdiaz/.conda/envs/snspd/bin/python \
  pipelines/14_dimensional_basecases.py \
  --campaign-dir /home/jdiaz/scratch/photon_runs/dimensional_v1_2_1_pair_20261019_01 \
  --cases bias21p5_ee4_775 bias21p5_ee4_1550 \
  --length-1d-nm 1000 \
  --threads-per-case 4 --physical-cores-per-case 2 \
  --ceiling-ps 100 --final-time-ns 15
```
