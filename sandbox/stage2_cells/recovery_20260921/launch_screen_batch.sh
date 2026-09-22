#!/usr/bin/env bash
# Foreground work inside the explicitly authorized screen. No retry or polling.
set -u
set -o pipefail
cd /home/jdiaz/pysnspd || exit 1
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
runroot=/home/jdiaz/pysnspd/tmp/stage2_recovery_20260921_ssp
console="${runroot}.screen.log"
if [[ -e "$runroot" || -e "$console" ]]; then
    printf 'Existing recovery evidence preserved. No relaunch or overwrite performed.\n'
    exec bash -i
fi
/home/jdiaz/.conda/envs/snspd/bin/python -u \
    sandbox/stage2_cells/recovery_20260921/run_limited_acceptance_batch.py \
    --plan docs/implementation/stage2/recovery_20260921/limited_acceptance_plan.json \
    --output-root "$runroot" --execute 2>&1 | tee "$console"
code=${PIPESTATUS[0]}
printf '\nRECOVERY_BATCH_EXIT_CODE=%s\nResults: %s\n' "$code" "$runroot"
exec bash -i
