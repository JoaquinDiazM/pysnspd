"""Apply a reviewed reachability report to the pySNSPD library.

The report must be produced by ``tools/reachability_audit.py`` from the same
Git revision. Package ``__init__.py`` files are retained even when they are not
pipeline roots; their exports are cleaned separately after the move plan.
"""

from __future__ import annotations

import argparse
import ast
import json
from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class Range:
    start: int
    end: int
    qualname: str


class RangeCollector(ast.NodeVisitor):
    def __init__(self) -> None:
        self.stack: list[str] = []
        self.ranges: dict[str, Range] = {}

    def visit_ClassDef(self, node: ast.ClassDef) -> None:
        self._visit_definition(node)

    def visit_FunctionDef(self, node: ast.FunctionDef) -> None:
        self._visit_definition(node)

    def visit_AsyncFunctionDef(self, node: ast.AsyncFunctionDef) -> None:
        self._visit_definition(node)

    def _visit_definition(
        self,
        node: ast.ClassDef | ast.FunctionDef | ast.AsyncFunctionDef,
    ) -> None:
        qualname = ".".join([*self.stack, node.name])
        decorator_lines = [decorator.lineno for decorator in node.decorator_list]
        start = min([node.lineno, *decorator_lines])
        self.ranges[qualname] = Range(
            start=start,
            end=node.end_lineno or node.lineno,
            qualname=qualname,
        )
        self.stack.append(node.name)
        self.generic_visit(node)
        self.stack.pop()


def _collapse_blank_lines(lines: list[str]) -> list[str]:
    output: list[str] = []
    blank_count = 0
    for line in lines:
        if line.strip():
            blank_count = 0
            output.append(line)
            continue
        blank_count += 1
        if blank_count <= 2:
            output.append(line)
    while output and not output[-1].strip():
        output.pop()
    return [*output, ""]


def _outermost_ranges(ranges: list[Range]) -> list[Range]:
    selected: list[Range] = []
    for candidate in sorted(ranges, key=lambda item: (item.start, -item.end)):
        if any(
            parent.start <= candidate.start and candidate.end <= parent.end
            for parent in selected
        ):
            continue
        selected.append(candidate)
    return selected


def _prune_definitions(path: Path, qualnames: set[str]) -> None:
    source = path.read_text(encoding="utf-8")
    tree = ast.parse(source, filename=str(path))
    collector = RangeCollector()
    collector.visit(tree)
    missing = sorted(qualnames.difference(collector.ranges))
    if missing:
        raise RuntimeError(f"{path}: report definitions not found: {missing}")
    ranges = _outermost_ranges(
        [collector.ranges[qualname] for qualname in qualnames]
    )
    removed_lines: set[int] = set()
    for item in ranges:
        removed_lines.update(range(item.start, item.end + 1))
    lines = source.splitlines()
    retained = [
        line
        for line_number, line in enumerate(lines, start=1)
        if line_number not in removed_lines
    ]
    path.write_text(
        "\n".join(_collapse_blank_lines(retained)),
        encoding="utf-8",
    )


def prune(root: Path, report_path: Path) -> None:
    report = json.loads(report_path.read_text(encoding="utf-8"))
    if Path(report["root"]).name != root.name:
        raise RuntimeError(
            f"Report root {report['root']!r} does not match {str(root)!r}."
        )
    module_paths = {
        item["path"]
        for item in report["unreachable_modules"]
        if Path(item["path"]).name != "__init__.py"
    }
    definitions_by_path: dict[str, set[str]] = {}
    for item in report["unreachable_definitions"]:
        definitions_by_path.setdefault(item["path"], set()).add(item["qualname"])

    for relative_path in sorted(module_paths):
        target = (root / relative_path).resolve()
        target.relative_to(root.resolve())
        if not relative_path.startswith("pysnspd/"):
            raise RuntimeError(f"Refusing non-library deletion: {relative_path}")
        target.unlink()
        print(f"deleted module: {relative_path}")

    for relative_path, qualnames in sorted(definitions_by_path.items()):
        if relative_path in module_paths:
            continue
        target = (root / relative_path).resolve()
        target.relative_to(root.resolve())
        _prune_definitions(target, qualnames)
        print(f"pruned {len(qualnames):3d} definitions: {relative_path}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path.cwd())
    parser.add_argument("--report", type=Path, required=True)
    args = parser.parse_args()
    prune(args.root.resolve(), args.report.resolve())


if __name__ == "__main__":
    main()
