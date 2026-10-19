# Current commands: pysnspd v1.2.0

The previous v1.1.3 commands/results remain historical. Do not interrupt their active campaign.

## 2026-10-19 | v1.2.0 guarded full 2D-to-1D campaign

Purpose: test dimensional admission and post-detection recovery of the SAME effective neutral BCS/KWT/native-kinetic/circuit model. Six cases: (15.5 uA, tau_ee(Tc)=5 ps), (15.5 uA, 4 ps), (21.5 uA, 4 ps), each at 775/1550 nm. Shared Rs=425 ohm, eta=.8, escape=20 ps, width=80 nm, length=300 nm, D rechecked from the inherited depairing calibration.

Official 2D through positive 0.1 mV crossing +2 ps (0.5 ps confirmation included). Independent native 1D shadow over the final interval. Only admitted cases switch to 61-node full 1D, with final target 15 ns from photon. Same finite strip and frozen exterior inductance; no tails. Domain-cut fraction >2%, failed transverse/overlap guards, physical support/capacity/stiffness or current reversal preserve a diagnostic checkpoint and stop. Final-time completion is not physical or reset certification.

Resources: six concurrent processes, four Numba SMT threads on two physical cores each (24 worker threads / 12 cores), one separate coordinator, at least two whole physical cores reserved. BLAS/OpenMP=1. Conservative RAM reservation27 GiB. Runtime lock and active-affinity inventory prevent overlap with the still-running v1.1.3 campaign. Wait for that campaign to finish. Around70% CPU use is an aim, not a measured speedup.

Outputs: campaign plan/resources/execution/logs in scratch; each run has phase-aware data.h5, native accepted-state checkpoint_latest.npz, manifest, DC/source/mesh/material identity. Exact --resume only, same sources/parameters/DC/mesh. A stopped diagnostic is preserved and cannot be resumed as if validated.

Estimate: long; 2D detection may still take hours. Short hot-state 1D control only covered0.1 ps on one core; it cannot establish total15ns runtime. A 150 s DC hold benchmark was intentionally bounded and incomplete. Do not extrapolate it as photon performance.

Execution (foreground, user-run; not launched by preparation):
```bash
cd /home/jdiaz/pysnspd
env PYTHONPYCACHEPREFIX=/home/jdiaz/scratch/tmp/pysnspd/pycache \
  MPLCONFIGDIR=/home/jdiaz/scratch/tmp/pysnspd/matplotlib \
  NUMBA_CACHE_DIR=/home/jdiaz/scratch/tmp/pysnspd/numba_cache \
  NUMBA_NUM_THREADS=4 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
  MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 \
  /home/jdiaz/.conda/envs/snspd/bin/python \
  pipelines/14_dimensional_basecases.py \
  --campaign-dir /home/jdiaz/scratch/photon_runs/dimensional_v1_2_20261019_01 \
  --threads-per-case 4 --physical-cores-per-case 2 \
  --ceiling-ps 100 --final-time-ns 15
```

To prepare without computing, append --prepare-only. To resume an interrupted campaign, append --resume. A prepared plan may be executed without --resume if no execution exists.

Phase-aware individual-case visualization (intrinsic Vout without amplifier/filter; requires committed manifests, no comparison of distinct cases):
```bash
cd /home/jdiaz/pysnspd
env PYTHONPYCACHEPREFIX=/home/jdiaz/scratch/tmp/pysnspd/pycache \
  MPLCONFIGDIR=/home/jdiaz/scratch/tmp/pysnspd/matplotlib \
  OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 \
  /home/jdiaz/.conda/envs/snspd/bin/python \
  plot_pipelines/08_dimensional_case_report.py \
  --run-dir /home/jdiaz/scratch/photon_runs/dimensional_v1_2_20261019_01/runs/* \
  --snapshots-ps 1 11 20 30 50 70 \
  --output-dir /home/jdiaz/scratch/case_reports/dimensional_v1_2_20261019_01 \
  --workers 2
```
