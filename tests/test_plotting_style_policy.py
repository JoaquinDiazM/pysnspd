"""Static and runtime checks for the shared thesis plotting policy."""

from __future__ import annotations

import ast
from pathlib import Path

import matplotlib.pyplot as plt

from pysnspd.plotting.style import THESIS_DPI, apply_thesis_style


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
PLOTTING_ROOT = REPOSITORY_ROOT / "pysnspd" / "plotting"
PLOT_PIPELINES_ROOT = REPOSITORY_ROOT / "plot_pipelines"


def _python_sources(directory: Path) -> list[Path]:
    return sorted(directory.glob("*.py"))


def _creates_figure(tree: ast.AST) -> bool:
    for node in ast.walk(tree):
        if not isinstance(node, ast.Call) or not isinstance(node.func, ast.Attribute):
            continue
        if (
            isinstance(node.func.value, ast.Name)
            and node.func.value.id == "plt"
            and node.func.attr in {"figure", "subplots"}
        ):
            return True
    return False


def _style_imports(tree: ast.Module) -> set[str]:
    imported: set[str] = set()
    for node in tree.body:
        if isinstance(node, ast.ImportFrom) and node.module == "pysnspd.plotting.style":
            imported.update(alias.name for alias in node.names)
    return imported


def _applies_style_at_module_scope(tree: ast.Module) -> bool:
    for node in tree.body:
        if not isinstance(node, ast.Expr) or not isinstance(node.value, ast.Call):
            continue
        if isinstance(node.value.func, ast.Name) and node.value.func.id == "apply_thesis_style":
            return True
    return False


def _dpi_defaults(function: ast.FunctionDef | ast.AsyncFunctionDef):
    positional = list(function.args.posonlyargs) + list(function.args.args)
    positional_defaults = [None] * (len(positional) - len(function.args.defaults))
    positional_defaults.extend(function.args.defaults)
    yield from zip(positional, positional_defaults)
    yield from zip(function.args.kwonlyargs, function.args.kw_defaults)


def test_every_figure_producer_applies_shared_style():
    failures: list[str] = []
    for path in _python_sources(PLOTTING_ROOT):
        tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
        if not _creates_figure(tree):
            continue
        imports = _style_imports(tree)
        if "apply_thesis_style" not in imports:
            failures.append(f"{path.name}: missing apply_thesis_style import")
        if not _applies_style_at_module_scope(tree):
            failures.append(f"{path.name}: missing module-scope apply_thesis_style()")
    assert not failures, "\n".join(failures)


def test_plotting_dpi_defaults_use_shared_constant():
    failures: list[str] = []
    for path in _python_sources(PLOTTING_ROOT):
        tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
        for node in ast.walk(tree):
            if not isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
                continue
            for argument, default in _dpi_defaults(node):
                if argument.arg != "dpi":
                    continue
                if default is None:
                    continue
                if not isinstance(default, ast.Name) or default.id != "THESIS_DPI":
                    failures.append(f"{path.name}:{node.lineno} {node.name}(dpi=...)")
    assert not failures, "\n".join(failures)


def test_plot_pipeline_dpi_defaults_use_shared_constant():
    failures: list[str] = []
    for path in _python_sources(PLOT_PIPELINES_ROOT):
        tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
        for node in ast.walk(tree):
            if not isinstance(node, ast.Call) or not node.args:
                continue
            if not isinstance(node.args[0], ast.Constant) or node.args[0].value != "--dpi":
                continue
            default = next(
                (keyword.value for keyword in node.keywords if keyword.arg == "default"),
                None,
            )
            if not isinstance(default, ast.Name) or default.id != "THESIS_DPI":
                failures.append(f"{path.name}:{node.lineno} --dpi")
    assert not failures, "\n".join(failures)


def test_figures_do_not_embed_an_independent_numeric_width():
    failures: list[str] = []
    for path in _python_sources(PLOTTING_ROOT):
        tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
        for node in ast.walk(tree):
            if not isinstance(node, ast.Call):
                continue
            for keyword in node.keywords:
                if keyword.arg != "figsize" or not isinstance(keyword.value, ast.Tuple):
                    continue
                width = keyword.value.elts[0] if keyword.value.elts else None
                if isinstance(width, ast.Constant) and isinstance(width.value, (int, float)):
                    failures.append(f"{path.name}:{node.lineno} numeric figsize width")
    assert not failures, "\n".join(failures)


def test_deprecated_style_todos_are_gone():
    failures = [
        path.name
        for path in _python_sources(PLOT_PIPELINES_ROOT)
        if "TODO(plot-style)" in path.read_text(encoding="utf-8")
    ]
    assert not failures


def test_apply_thesis_style_sets_runtime_defaults():
    apply_thesis_style()
    assert plt.rcParams["figure.dpi"] == THESIS_DPI
    assert plt.rcParams["text.usetex"] is False
    assert plt.rcParams["font.family"] == ["serif"]
    assert plt.rcParams["axes.grid"] is False
