"""Restore explicitly archived report/data paths without overwriting user work.

The manifest is intentionally fixed. A byte-identical retained copy is preferred;
otherwise the committed Git object is read locally. No network request is made.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path, PurePosixPath
import re
import subprocess


ROOT = Path(__file__).resolve().parents[2]
MANIFEST = ROOT / "docs/implementation/archive_20260924/manifest.json"


def digest(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def inside(root: Path, relative: str) -> Path:
    """Validate lexical and resolved containment, including existing symlinks."""
    part = PurePosixPath(relative)
    if (not relative or part.is_absolute() or ".." in part.parts
            or "\\" in relative or ":" in relative):
        raise ValueError(f"Unsafe archived path: {relative!r}")
    resolved = (root / relative).resolve()
    if not resolved.is_relative_to(root.resolve()):
        raise ValueError(f"Archived path escapes its destination: {relative!r}")
    return resolved


def source_bytes(entry: dict) -> tuple[bytes, str]:
    expected = entry["sha256"]
    if not re.fullmatch(r"[0-9a-f]{64}", expected):
        raise ValueError(f"Invalid SHA256 for {entry['path']}")
    canonical = entry.get("canonical_copy")
    if canonical:
        candidate = inside(ROOT, canonical)
        if candidate.is_file():
            payload = candidate.read_bytes()
            if digest(payload) == expected:
                return payload, canonical
    commit = entry["source_commit"]
    if not re.fullmatch(r"[0-9a-f]{40}", commit):
        raise ValueError(f"Invalid source commit for {entry['path']}")
    result = subprocess.run(
        ["git", "show", f"{commit}:{entry['path']}"], cwd=ROOT,
        capture_output=True, check=False,
    )
    if result.returncode:
        raise RuntimeError(
            f"Cannot read archived Git object for {entry['path']}: "
            + result.stderr.decode("utf-8", errors="replace").strip()
        )
    if digest(result.stdout) != expected:
        raise ValueError(f"Git object SHA256 mismatch for {entry['path']}")
    return result.stdout, f"git:{commit}"


def restore(entry: dict, destination: Path) -> dict:
    target = inside(destination, entry["path"])
    if target.exists():
        if target.is_file() and digest(target.read_bytes()) == entry["sha256"]:
            return {"path": entry["path"], "status": "ALREADY_IDENTICAL"}
        raise FileExistsError(f"Refusing to overwrite different content: {target}")
    payload, source = source_bytes(entry)
    if len(payload) != entry["bytes"]:
        raise ValueError(f"Size mismatch for {entry['path']}")
    target.parent.mkdir(parents=True, exist_ok=True)
    # Recheck after making parents; exclusive creation also protects a concurrent
    # writer from accidental replacement between existence check and creation.
    target = inside(destination, entry["path"])
    with target.open("xb") as stream:
        stream.write(payload)
    return {"path": entry["path"], "status": "RESTORED", "source": source,
            "bytes": len(payload), "sha256": entry["sha256"]}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    selection = parser.add_mutually_exclusive_group(required=True)
    selection.add_argument("--list", action="store_true", help="List archived paths")
    selection.add_argument("--all", action="store_true", help="Restore all listed paths")
    selection.add_argument("--path", action="append", help="Restore one listed path; repeatable")
    parser.add_argument(
        "--destination-root", default=".",
        help="Destination within this repository (default: repository root)",
    )
    args = parser.parse_args()
    manifest = json.loads(MANIFEST.read_text(encoding="utf-8"))
    records = {entry["path"]: entry for entry in manifest["artifacts"]}
    if len(records) != len(manifest["artifacts"]):
        parser.error("Archive manifest contains duplicate paths")
    for relative in records:
        inside(ROOT, relative)
    if args.list:
        for entry in records.values():
            print(f"{entry['bytes']:>10}  {entry['category']:<34} {entry['path']}")
        return
    destination = (ROOT / args.destination_root).resolve()
    if not destination.is_relative_to(ROOT):
        parser.error("Destination must remain within the repository")
    paths = list(records) if args.all else args.path
    unknown = sorted(set(paths) - records.keys())
    if unknown:
        parser.error("Paths are not listed in the archive: " + ", ".join(unknown))
    for relative in paths:
        print(json.dumps(restore(records[relative], destination), ensure_ascii=False))


if __name__ == "__main__":
    main()
