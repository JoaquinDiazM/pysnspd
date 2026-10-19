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
