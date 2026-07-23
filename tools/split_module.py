"""Split a Python module before a named top-level definition.

The command preserves the original preamble in both files and creates explicit
cross-module imports for names referenced across the split. It is intended for
reviewed structural refactors, not arbitrary source generation.
"""

from __future__ import annotations

import argparse
import ast
from pathlib import Path


def _top_level_definition(tree: ast.Module, name: str) -> ast.AST:
    for node in tree.body:
        if isinstance(
            node,
            (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef),
        ) and node.name == name:
            return node
    raise RuntimeError(f"Top-level definition {name!r} was not found.")


def _first_definition_line(tree: ast.Module) -> int:
    lines = []
    for node in tree.body:
        if not isinstance(
            node,
            (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef),
        ):
            continue
        decorator_lines = [decorator.lineno for decorator in node.decorator_list]
        lines.append(min([node.lineno, *decorator_lines]))
    if not lines:
        raise RuntimeError("Module has no top-level definitions.")
    return min(lines)


def _defined_names(tree: ast.Module) -> set[str]:
    names: set[str] = set()
    for node in tree.body:
        if isinstance(
            node,
            (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef),
        ):
            names.add(node.name)
        elif isinstance(node, (ast.Assign, ast.AnnAssign)):
            targets = node.targets if isinstance(node, ast.Assign) else [node.target]
            for target in targets:
                if isinstance(target, ast.Name):
                    names.add(target.id)
    return names


def _loaded_names(tree: ast.AST) -> set[str]:
    return {
        node.id
        for node in ast.walk(tree)
        if isinstance(node, ast.Name) and isinstance(node.ctx, ast.Load)
    }


def _module_name(root: Path, path: Path) -> str:
    relative = path.resolve().relative_to(root.resolve()).with_suffix("")
    parts = list(relative.parts)
    if parts[-1] == "__init__":
        parts.pop()
    return ".".join(parts)


def _render_import(module: str, names: set[str]) -> list[str]:
    if not names:
        return []
    if len(names) == 1:
        return [f"from {module} import {next(iter(names))}"]
    return [
        f"from {module} import (",
        *[f"    {name}," for name in sorted(names)],
        ")",
    ]


def split_module(
    root: Path,
    source: Path,
    target: Path,
    before: str,
) -> None:
    source = source.resolve()
    target = target.resolve()
    source.relative_to(root.resolve())
    target.relative_to(root.resolve())
    if target.exists():
        raise RuntimeError(f"Target already exists: {target}")

    text = source.read_text(encoding="utf-8")
    lines = text.splitlines()
    tree = ast.parse(text, filename=str(source))
    split_node = _top_level_definition(tree, before)
    preamble_end = _first_definition_line(tree) - 1
    preamble = lines[:preamble_end]
    head_body = lines[preamble_end : split_node.lineno - 1]
    tail_body = lines[split_node.lineno - 1 :]

    head_text = "\n".join([*preamble, *head_body, ""]) + "\n"
    target_text = "\n".join([*preamble, *tail_body, ""]) + "\n"
    head_tree = ast.parse(head_text, filename=str(source))
    target_tree = ast.parse(target_text, filename=str(target))

    head_defined = _defined_names(head_tree)
    target_defined = _defined_names(target_tree)
    head_needs = _loaded_names(head_tree).intersection(target_defined)
    target_needs = _loaded_names(target_tree).intersection(head_defined)

    source_module = _module_name(root, source)
    target_module = _module_name(root, target)
    head_imports = _render_import(target_module, head_needs)
    target_imports = _render_import(source_module, target_needs)

    source.write_text(
        "\n".join([head_text.rstrip(), "", *head_imports, ""]).rstrip() + "\n",
        encoding="utf-8",
    )
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(
        "\n".join(
            [
                "\n".join(preamble).rstrip(),
                "",
                *target_imports,
                "",
                "\n".join(tail_body).rstrip(),
                "",
            ]
        ).lstrip(),
        encoding="utf-8",
    )
    print(f"split {source.relative_to(root)} -> {target.relative_to(root)}")
    if head_needs:
        print(f"  source imports: {', '.join(sorted(head_needs))}")
    if target_needs:
        print(f"  target imports: {', '.join(sorted(target_needs))}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path.cwd())
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--target", type=Path, required=True)
    parser.add_argument("--before", required=True)
    args = parser.parse_args()
    root = args.root.resolve()
    split_module(
        root,
        root / args.source,
        root / args.target,
        args.before,
    )


if __name__ == "__main__":
    main()
