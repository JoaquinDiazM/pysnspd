# Current commands: pysnspd v1.2.3

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
