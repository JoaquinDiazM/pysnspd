"""Preserve one atomic full checkpoint and its exact runtime source files.

This performs one read/copy operation; it never starts, waits for or polls a
simulation. All preserved numeric states and code receipts remain in scratch.
"""
from pathlib import Path
import argparse
import hashlib
import json
import os
import shutil
import numpy as np

SCRATCH = Path("/home/jdiaz/scratch").resolve()
REPO = Path(__file__).resolve().parents[2]

def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def preserve(run_dir):
    run = Path(run_dir).resolve()
    if not run.is_relative_to(SCRATCH):
        raise ValueError("Only scratch trajectories may be archived")
    manifest = json.loads((run/"manifest.json").read_text())
    archive = run/"runtime_sources"
    for relative, sha in manifest["source_pins_sha256"].items():
        rel = Path(relative)
        if rel.is_absolute() or ".." in rel.parts:
            raise ValueError("Invalid relative source path")
        source = (REPO/rel).resolve()
        if not source.is_relative_to(REPO.resolve()):
            raise ValueError("Source escaped repository")
        target = archive/rel
        if target.exists():
            if digest(target) != sha:
                raise ValueError("Preserved source changed")
            continue
        if digest(source) != sha:
            raise ValueError("Current source differs from the run; preserve earlier source manually")
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source, target)
        if digest(target) != sha:
            raise ValueError("Source copy verification failed")
    with np.load(run/"checkpoint_latest.npz", allow_pickle=False) as saved:
        meta = json.loads(str(saved["metadata_json"].item()))
        destination = run/"archived_full_states"/("state_%013.6fps.npz"%meta["time_ps"])
        if not destination.exists():
            destination.parent.mkdir(exist_ok=True)
            temporary = destination.with_suffix(".part")
            arrays = {key: saved[key].copy() for key in
                      ("delta_bar", "p", "n", "circuit_SI", "metadata_json")}
            with temporary.open("wb") as stream:
                np.savez_compressed(stream, **arrays)
            os.replace(temporary, destination)
    return dict(run_dir=str(run), time_ps=meta["time_ps"],
                full_state_path=str(destination), state_sha256=digest(destination),
                runtime_sources_preserved=len(manifest["source_pins_sha256"]))

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-dir", required=True, action="append", type=Path)
    args = parser.parse_args()
    print(json.dumps([preserve(p) for p in args.run_dir], indent=2))

if __name__ == "__main__":
    main()
