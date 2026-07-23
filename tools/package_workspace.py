"""Create a portable source ZIP without VCS, cache, or local report files."""

from __future__ import annotations

import argparse
import stat
import zipfile
from pathlib import Path


EXCLUDED_PARTS = {
    ".git",
    ".pytest_cache",
    "__pycache__",
    "tmp",
}


def package(root: Path, output: Path) -> None:
    files = [
        path
        for path in root.rglob("*")
        if path.is_file()
        and not EXCLUDED_PARTS.intersection(path.relative_to(root).parts)
        and path.resolve() != output.resolve()
    ]
    output.parent.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(output, "w", compression=zipfile.ZIP_DEFLATED) as archive:
        for path in sorted(files):
            relative = path.relative_to(root).as_posix()
            info = zipfile.ZipInfo(relative)
            info.compress_type = zipfile.ZIP_DEFLATED
            info.external_attr = (stat.S_IFREG | 0o644) << 16
            archive.writestr(info, path.read_bytes())
    print(f"{len(files)} files -> {output}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path.cwd())
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    package(args.root.resolve(), args.output.resolve())


if __name__ == "__main__":
    main()
