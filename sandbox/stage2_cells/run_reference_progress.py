"""User-run reference with flushed progress; the physical RHS is unchanged.

This observer is intended for the computation handoff after the bounded probe
remained incomplete. It does not resume, subdivide or relaunch that calculation.
Every line reports an RHS evaluation time, which may move backwards when the
adaptive integrator rejects a trial step. No asynchronous/background job is used.
"""
from pathlib import Path
import hashlib
import json
import sys
import time

ROOT = Path(__file__).resolve().parents[2]
sys.path[:0] = [str(ROOT), str(Path(__file__).resolve().parent)]
import run_coupled
from pysnspd.experimental.coupled_cells import CoupledCellSystem


def main():
    if '--output' not in sys.argv:
        raise SystemExit('An explicit --output path is required for this user-run diagnostic.')
    output = Path(sys.argv[sys.argv.index('--output')+1])
    output.parent.mkdir(parents=True, exist_ok=True)
    progress = output.with_suffix('.progress.jsonl')
    if progress.exists() or output.exists():
        raise SystemExit('Choose a new output path; existing diagnostic evidence is preserved.')
    started = time.perf_counter()
    original = CoupledCellSystem.rhs
    observer_hash = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    with progress.open('w', encoding='utf-8', buffering=1) as log:
        def emit(row):
            line = json.dumps(row, allow_nan=False)
            print(line, flush=True)
            log.write(line+'\n')
            log.flush()

        emit(dict(event='START', observer_sha256=observer_hash, arguments=sys.argv[1:],
                  meaning='t_rhs is a trial evaluation time, not a guaranteed accepted step'))

        def observed_rhs(system, time_value, state):
            result = original(system, time_value, state)
            if system.rhs_calls == 1 or system.rhs_calls % 100 == 0:
                amplitudes, p, n = system.unpack(state)
                emit(dict(event='RHS', evaluations=system.rhs_calls,
                    elapsed_seconds=time.perf_counter()-started, t_rhs=float(time_value),
                    amplitudes=amplitudes.tolist(), minimum_electron=float(p.min()),
                    maximum_electron=float(p.max()), minimum_phonon=float(n.min())))
            return result

        CoupledCellSystem.rhs = observed_rhs
        try:
            run_coupled.main()
        finally:
            CoupledCellSystem.rhs = original
            if output.exists():
                record = json.loads(output.read_text(encoding='utf-8'))
                record['progress_observer'] = dict(source_sha256=observer_hash,
                    path=progress.as_posix(), physical_rhs_unchanged=True)
                output.write_text(json.dumps(record, indent=2, allow_nan=False)+'\n', encoding='utf-8')
            emit(dict(event='STOP', elapsed_seconds=time.perf_counter()-started))


if __name__ == '__main__':
    main()
