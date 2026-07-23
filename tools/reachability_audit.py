"""Audit production reachability from pySNSPD pipeline entry points.

The tool combines Python's AST with Jedi's project-wide reference resolver.
It reports definitions and modules that cannot be reached from module-level
execution in ``pipelines`` or ``plot_pipelines``. Tests are deliberately not
roots: the deletion policy is defined by production entry points.
"""

from __future__ import annotations

import argparse
import ast
import json
from collections import defaultdict, deque
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Iterable

import jedi


SOURCE_DIRS = ("pysnspd", "pipelines", "plot_pipelines")
ROOT_DIRS = ("pipelines", "plot_pipelines")


@dataclass(frozen=True)
class Definition:
    key: str
    module: str
    path: str
    qualname: str
    kind: str
    line: int
    end_line: int
    name_column: int
    parent_class: str | None
    is_special_method: bool


@dataclass(frozen=True)
class Module:
    key: str
    name: str
    path: str
    is_pipeline: bool


class DefinitionCollector(ast.NodeVisitor):
    def __init__(self, module: str, path: str) -> None:
        self.module = module
        self.path = path
        self.stack: list[tuple[str, str]] = []
        self.definitions: list[Definition] = []

    def visit_ClassDef(self, node: ast.ClassDef) -> None:
        self._record(node, "class")
        self.stack.append((node.name, "class"))
        self.generic_visit(node)
        self.stack.pop()

    def visit_FunctionDef(self, node: ast.FunctionDef) -> None:
        self._visit_function(node)

    def visit_AsyncFunctionDef(self, node: ast.AsyncFunctionDef) -> None:
        self._visit_function(node)

    def _visit_function(
        self, node: ast.FunctionDef | ast.AsyncFunctionDef
    ) -> None:
        self._record(node, "function")
        self.stack.append((node.name, "function"))
        self.generic_visit(node)
        self.stack.pop()

    def _record(
        self,
        node: ast.ClassDef | ast.FunctionDef | ast.AsyncFunctionDef,
        kind: str,
    ) -> None:
        qualname = ".".join([name for name, _ in self.stack] + [node.name])
        parent_class = next(
            (name for name, parent_kind in reversed(self.stack) if parent_kind == "class"),
            None,
        )
        self.definitions.append(
            Definition(
                key=f"def:{self.module}:{qualname}",
                module=self.module,
                path=self.path,
                qualname=qualname,
                kind=kind,
                line=node.lineno,
                end_line=node.end_lineno or node.lineno,
                name_column=_definition_name_column(node),
                parent_class=parent_class,
                is_special_method=(
                    kind == "function"
                    and parent_class is not None
                    and node.name.startswith("__")
                    and node.name.endswith("__")
                ),
            )
        )


def _definition_name_column(
    node: ast.ClassDef | ast.FunctionDef | ast.AsyncFunctionDef,
) -> int:
    if isinstance(node, ast.ClassDef):
        return node.col_offset + len("class ")
    if isinstance(node, ast.AsyncFunctionDef):
        return node.col_offset + len("async def ")
    return node.col_offset + len("def ")


def _module_name(root: Path, path: Path) -> str:
    relative = path.relative_to(root).with_suffix("")
    parts = list(relative.parts)
    if parts[-1] == "__init__":
        parts.pop()
    return ".".join(parts)


def _python_files(root: Path) -> list[Path]:
    files: list[Path] = []
    for directory in SOURCE_DIRS:
        files.extend(sorted((root / directory).rglob("*.py")))
    return files


def _resolve_import_module(
    current_module: str,
    imported_module: str | None,
    level: int,
    current_is_package: bool,
) -> str | None:
    if level == 0:
        return imported_module
    package_parts = current_module.split(".")
    if not current_is_package:
        package_parts.pop()
    remove = level - 1
    if remove > len(package_parts):
        return None
    base = package_parts[: len(package_parts) - remove]
    if imported_module:
        base.extend(imported_module.split("."))
    return ".".join(base)


def _containing_node(
    module: Module,
    definitions: list[Definition],
    line: int,
) -> str:
    candidates = [
        definition
        for definition in definitions
        if definition.line <= line <= definition.end_line
    ]
    if not candidates:
        return module.key
    return min(
        candidates,
        key=lambda definition: definition.end_line - definition.line,
    ).key


def _iter_import_edges(
    root: Path,
    module: Module,
    tree: ast.AST,
    definitions: list[Definition],
    module_names: set[str],
) -> Iterable[tuple[str, str]]:
    is_package = Path(module.path).name == "__init__.py"
    for node in ast.walk(tree):
        if not isinstance(node, (ast.Import, ast.ImportFrom)):
            continue
        caller = _containing_node(module, definitions, node.lineno)
        targets: set[str] = set()
        if isinstance(node, ast.Import):
            for alias in node.names:
                parts = alias.name.split(".")
                for index in range(len(parts), 0, -1):
                    candidate = ".".join(parts[:index])
                    if candidate in module_names:
                        targets.add(candidate)
                        break
        else:
            base = _resolve_import_module(
                module.name,
                node.module,
                node.level,
                is_package,
            )
            if base:
                if base in module_names:
                    targets.add(base)
                for alias in node.names:
                    candidate = f"{base}.{alias.name}"
                    if candidate in module_names:
                        targets.add(candidate)
        for target in targets:
            yield caller, f"module:{target}"


def _collect_project(root: Path) -> tuple[
    dict[str, Module],
    dict[str, Definition],
    dict[str, ast.AST],
    dict[str, list[Definition]],
]:
    modules: dict[str, Module] = {}
    definitions: dict[str, Definition] = {}
    trees: dict[str, ast.AST] = {}
    definitions_by_module: dict[str, list[Definition]] = {}
    for path in _python_files(root):
        module_name = _module_name(root, path)
        relative_path = path.relative_to(root).as_posix()
        module = Module(
            key=f"module:{module_name}",
            name=module_name,
            path=relative_path,
            is_pipeline=relative_path.split("/", 1)[0] in ROOT_DIRS,
        )
        source = path.read_text(encoding="utf-8")
        tree = ast.parse(source, filename=str(path))
        collector = DefinitionCollector(module_name, relative_path)
        collector.visit(tree)
        modules[module_name] = module
        trees[module_name] = tree
        definitions_by_module[module_name] = collector.definitions
        definitions.update(
            (definition.key, definition) for definition in collector.definitions
        )
    return modules, definitions, trees, definitions_by_module


def _reference_edges(
    root: Path,
    modules: dict[str, Module],
    definitions: dict[str, Definition],
    definitions_by_module: dict[str, list[Definition]],
) -> tuple[set[tuple[str, str]], list[dict[str, str | int]]]:
    project = jedi.Project(path=str(root))
    scripts: dict[str, jedi.Script] = {}
    module_by_path = {
        str((root / module.path).resolve()).lower(): module
        for module in modules.values()
    }
    edges: set[tuple[str, str]] = set()
    failures: list[dict[str, str | int]] = []
    for definition in definitions.values():
        script = scripts.get(definition.path)
        if script is None:
            script = jedi.Script(
                path=str(root / definition.path),
                project=project,
            )
            scripts[definition.path] = script
        try:
            references = script.get_references(
                line=definition.line,
                column=definition.name_column,
                include_builtins=False,
                scope="project",
            )
        except Exception as exc:  # pragma: no cover - diagnostic safeguard
            failures.append(
                {
                    "definition": definition.key,
                    "line": definition.line,
                    "error": repr(exc),
                }
            )
            continue
        for reference in references:
            if reference.is_definition() or reference.module_path is None:
                continue
            path_key = str(Path(reference.module_path).resolve()).lower()
            caller_module = module_by_path.get(path_key)
            if caller_module is None:
                continue
            caller = _containing_node(
                caller_module,
                definitions_by_module[caller_module.name],
                reference.line,
            )
            if caller != definition.key:
                edges.add((caller, definition.key))
    return edges, failures


def _class_special_method_edges(
    definitions: dict[str, Definition],
) -> set[tuple[str, str]]:
    classes = {
        (definition.module, definition.qualname): definition
        for definition in definitions.values()
        if definition.kind == "class"
    }
    edges: set[tuple[str, str]] = set()
    for definition in definitions.values():
        if not definition.is_special_method:
            continue
        class_qualname = definition.qualname.rsplit(".", 1)[0]
        parent = classes.get((definition.module, class_qualname))
        if parent is not None:
            edges.add((parent.key, definition.key))
    return edges


def _attribute_method_edges(
    modules: dict[str, Module],
    trees: dict[str, ast.AST],
    definitions: dict[str, Definition],
    definitions_by_module: dict[str, list[Definition]],
) -> set[tuple[str, str]]:
    """Resolve indirect instance calls that Jedi cannot always type-infer.

    A production expression such as ``diagnostics.as_dict()`` may not resolve
    to a concrete dataclass method when the object crosses a callback or
    factory boundary.  Attribute-name fallback edges are deliberately
    conservative: they only become live when the containing production node is
    itself reachable.
    """

    methods_by_name: dict[str, list[Definition]] = defaultdict(list)
    for definition in definitions.values():
        if definition.parent_class is not None:
            methods_by_name[definition.qualname.rsplit(".", 1)[-1]].append(
                definition
            )

    edges: set[tuple[str, str]] = set()
    for module_name, tree in trees.items():
        module = modules[module_name]
        module_definitions = definitions_by_module[module_name]
        for node in ast.walk(tree):
            if not isinstance(node, ast.Attribute) or not isinstance(
                node.ctx, ast.Load
            ):
                continue
            caller = _containing_node(
                module,
                module_definitions,
                node.lineno,
            )
            for target in methods_by_name.get(node.attr, ()):
                if caller != target.key:
                    edges.add((caller, target.key))
    return edges


def _reachable_nodes(
    roots: set[str],
    edges: set[tuple[str, str]],
) -> set[str]:
    adjacency: dict[str, set[str]] = defaultdict(set)
    for source, target in edges:
        adjacency[source].add(target)
    reachable = set(roots)
    queue = deque(sorted(roots))
    while queue:
        source = queue.popleft()
        for target in adjacency.get(source, ()):
            if target not in reachable:
                reachable.add(target)
                queue.append(target)
    return reachable


def audit(root: Path) -> dict[str, object]:
    modules, definitions, trees, definitions_by_module = _collect_project(root)
    module_names = set(modules)
    edges: set[tuple[str, str]] = set()
    for module_name, tree in trees.items():
        edges.update(
            _iter_import_edges(
                root,
                modules[module_name],
                tree,
                definitions_by_module[module_name],
                module_names,
            )
        )
    reference_edges, failures = _reference_edges(
        root,
        modules,
        definitions,
        definitions_by_module,
    )
    edges.update(reference_edges)
    edges.update(_class_special_method_edges(definitions))
    edges.update(
        _attribute_method_edges(
            modules,
            trees,
            definitions,
            definitions_by_module,
        )
    )
    roots = {
        module.key for module in modules.values() if module.is_pipeline
    }
    reachable = _reachable_nodes(roots, edges)
    library_modules = [
        module for module in modules.values() if module.path.startswith("pysnspd/")
    ]
    library_definitions = [
        definition
        for definition in definitions.values()
        if definition.path.startswith("pysnspd/")
    ]
    unreachable_modules = sorted(
        (
            module
            for module in library_modules
            if module.key not in reachable
        ),
        key=lambda module: module.path,
    )
    unreachable_definitions = sorted(
        (
            definition
            for definition in library_definitions
            if definition.key not in reachable
        ),
        key=lambda definition: (definition.path, definition.line),
    )
    reachable_definitions = sorted(
        (
            definition
            for definition in library_definitions
            if definition.key in reachable
        ),
        key=lambda definition: (definition.path, definition.line),
    )
    module_definition_counts: list[dict[str, object]] = []
    for module in sorted(library_modules, key=lambda item: item.path):
        module_definitions = definitions_by_module[module.name]
        module_definition_counts.append(
            {
                "path": module.path,
                "reachable_module": module.key in reachable,
                "definitions": len(module_definitions),
                "reachable_definitions": sum(
                    definition.key in reachable
                    for definition in module_definitions
                ),
                "unreachable_definitions": sum(
                    definition.key not in reachable
                    for definition in module_definitions
                ),
            }
        )
    return {
        "root": str(root),
        "roots": sorted(roots),
        "counts": {
            "modules": len(library_modules),
            "definitions": len(library_definitions),
            "reachable_modules": sum(
                module.key in reachable for module in library_modules
            ),
            "reachable_definitions": sum(
                definition.key in reachable for definition in library_definitions
            ),
            "unreachable_modules": len(unreachable_modules),
            "unreachable_definitions": len(unreachable_definitions),
            "edges": len(edges),
            "resolution_failures": len(failures),
        },
        "unreachable_modules": [asdict(module) for module in unreachable_modules],
        "unreachable_definitions": [
            asdict(definition) for definition in unreachable_definitions
        ],
        "reachable_definitions": [
            asdict(definition) for definition in reachable_definitions
        ],
        "module_definition_counts": module_definition_counts,
        "resolution_failures": failures,
        "edges": sorted([list(edge) for edge in edges]),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--root",
        type=Path,
        default=Path.cwd(),
        help="Repository root (default: current directory).",
    )
    parser.add_argument(
        "--output",
        type=Path,
        help="Optional JSON output path; stdout is used otherwise.",
    )
    args = parser.parse_args()
    report = audit(args.root.resolve())
    rendered = json.dumps(report, indent=2, sort_keys=True)
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered + "\n", encoding="utf-8")
    else:
        print(rendered)


if __name__ == "__main__":
    main()
