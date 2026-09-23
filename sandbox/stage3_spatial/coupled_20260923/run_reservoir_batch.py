"""Six registered spatial snapshots plus a saved-data reservoir-work diagnostic.

Run without --execute to describe the work. This runner never starts a background
job, retries a failure or overwrites an output. The detailed progress/ETA is
printed by the underlying mixed runner and retained in its progress.jsonl.
"""
from pathlib import Path
import argparse
from datetime import datetime, timezone
import hashlib
import json
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[3]
HERE = Path(__file__).resolve().parent
DATA = ROOT/'docs/implementation/stage3/coupled_20260923'


def save(path, value):
    path.write_text(json.dumps(value, indent=2, ensure_ascii=False, allow_nan=False)+'\n', encoding='utf-8')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-root')
    parser.add_argument('--execute', action='store_true')
    args = parser.parse_args()
    mixed = HERE/'run_mixed_check.py'
    post = HERE/'postprocess_saved_reservoir.py'
    seed = DATA/'seeded_check'
    files = [Path(__file__).resolve(), mixed, post, DATA/'mixed_registration.json',
             DATA/'reservoir_postprocess_registration.json', ROOT/'pysnspd/experimental/reservoir_spatial_dynamics.py']
    source_hashes = {p.relative_to(ROOT).as_posix(): hashlib.sha256(p.read_bytes()).hexdigest() for p in files}
    for path in files:
        if path.suffix == '.json':
            json.loads(path.read_text(encoding='utf-8'))
    if not args.execute:
        check = subprocess.run([sys.executable, str(mixed), '--seed-dir', str(seed)],
                               cwd=ROOT, text=True, capture_output=True, check=True)
        ready = json.loads(check.stdout)
        print(json.dumps(dict(status='READY_NOT_EXECUTED', cases=ready['cases'],
            scope='Six instantaneous meshes/profiles followed by prescribed reservoir load and work accounting; no time integration.',
            purpose='Measure reservoir-compatible instantaneous response and terminal normal-current defect; do not converge free-end heating.',
            output_layout=['mixed_control/', 'reservoir/', 'manifest.json', 'summary.json or failure.json'],
            expected_minutes=[8, 15], source_sha256=source_hashes,
            user_execution_required=True), indent=2))
        return
    if not args.output_root:
        parser.error('--output-root is required with --execute')
    out = (ROOT/args.output_root).resolve()
    if out.exists():
        raise FileExistsError('Choose a new output directory; no overwrite')
    out.mkdir(parents=True)
    started = time.monotonic()
    save(out/'manifest.json', dict(started_utc=datetime.now(timezone.utc).isoformat(),
        source_sha256=source_hashes,
        stage3_closed=False, temporal_admission=False, production=False))
    commands = [
        [sys.executable, '-u', str(mixed), '--seed-dir', str(seed), '--output-root', str(out/'mixed_control'), '--execute'],
        [sys.executable, '-u', str(post), '--input-root', str(out/'mixed_control'),
         '--output-root', str(out/'reservoir'), '--registration', str(DATA/'reservoir_postprocess_registration.json')],
    ]
    try:
        for index, command in enumerate(commands, 1):
            print(json.dumps(dict(event='PHASE_START', phase=index, phases=2,
                name='Material snapshots' if index == 1 else 'Saved-data reservoir work',
                elapsed_seconds=time.monotonic()-started)), flush=True)
            subprocess.run(command, cwd=ROOT, check=True)
        result = json.loads((out/'reservoir/summary.json').read_text(encoding='utf-8'))
        if not result['status'].startswith('PASS'):
            raise RuntimeError('Reservoir diagnostic did not pass; saved outputs retained')
        save(out/'summary.json', dict(status='RESERVOIR_INSTANTANEOUS_BATCH_PASS_SCOPE_LIMITED',
            runtime_seconds=time.monotonic()-started, reservoir_status=result['status'],
            stage3_closed=False, temporal_admission=False, production=False))
        print(json.dumps(dict(event='SUCCESS', elapsed_seconds=time.monotonic()-started)), flush=True)
    except BaseException as error:
        save(out/'failure.json', dict(exception=type(error).__name__, reason=str(error),
            runtime_seconds=time.monotonic()-started, policy='No retry, fallback, overwrite or continuation'))
        raise


if __name__ == '__main__':
    main()
