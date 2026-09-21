"""Reproduce the batch's seven offline guards without any physical calculation.

The actual test implementation is run_acceptance_batch.self_test(). This entry
point records both that implementation's hash and the launcher hash, and never
overwrites an existing test result.
"""
from pathlib import Path
import argparse
import hashlib
import json
import sys

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(Path(__file__).resolve().parent))
import run_acceptance_batch


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    output = run_acceptance_batch.workspace_path(args.output)
    if output.exists():
        raise SystemExit('Choose a new output path; existing guard evidence is preserved.')
    result = run_acceptance_batch.self_test()
    result['guard_script_path'] = Path(__file__).relative_to(ROOT).as_posix()
    result['guard_script_sha256'] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    result['test_implementation'] = 'run_acceptance_batch.self_test'
    run_acceptance_batch.write_new(output, result)
    print(json.dumps(result, indent=2))
    if result['status'] != 'PASS':
        raise SystemExit(1)


if __name__ == '__main__':
    main()
