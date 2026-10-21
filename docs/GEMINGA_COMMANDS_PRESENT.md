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
