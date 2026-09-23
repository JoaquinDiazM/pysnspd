"""Read-only stage-4 delivery check with one preserved predecessor traversal.

Historical claims apply to historical bytes. The stage-3.5 decision that stage
4 had not started is not a restriction on the newly authorized implementation.
No physical calculations, tolerance gates or manifest generation occur here.
"""
from contextlib import redirect_stdout
from io import StringIO
from pathlib import Path
from types import SimpleNamespace
import argparse
import inspect
import json
import runpy
import sys


ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "docs/implementation/stage4/start_20260923"
PREVIOUS = DATA / "previous_delivery_exact"
MANIFEST = DATA / "delivery_manifest.json"
HISTORICAL = ROOT / "sandbox/stage3_5/assessment_r2_20260923/verify_delivery.py"
BASE_VERIFIER = ROOT / "sandbox/stage3_spatial/closure_20260923/verify_delivery.py"
ARCHIVED_NAMES = (
    "pysnspd/experimental/cell_closures.py",
    "pysnspd/experimental/spatial_functional.py",
    "pysnspd/experimental/mixed_spatial.py",
    "pysnspd/experimental/spatial_open.py",
    "docs/GEMINGA_COMMANDS.md",
    "docs/implementation/stage4/README.md",
)
EXPECTED_FLAGS = dict(stage4_started=True, stage4_complete=False,
    photon_transfer_admitted=False, production_promotion=False)


def historical_delivery():
    """Reuse frozen validators, changing only their historical file lookup.

    Existing predecessor overrides take precedence over this later archive:
    different generations have different historical command notebooks. Records
    already rooted in an older snapshot retain their original interpretation.
    The proxy is local to each loaded verifier, never a global runpy patch.
    """
    archived = {name: PREVIOUS / name for name in ARCHIVED_NAMES}
    for path in archived.values():
        if not path.is_file():
            raise RuntimeError("Missing preserved predecessor file: " + str(path))

    def load(path):
        namespace = runpy.run_path(str(path))
        globals_ = namespace["main"].__globals__
        original = namespace.get("verify")
        if original is not None:
            signature = inspect.signature(original)

            def verify(records, *args, **kwargs):
                bound = signature.bind(records, *args, **kwargs)
                bound.apply_defaults()
                lookup_root = Path(bound.arguments.get("root", ROOT)).resolve()
                if lookup_root == ROOT:
                    older = bound.arguments.get("overrides") or {}
                    bound.arguments["overrides"] = archived | older
                return original(*bound.args, **bound.kwargs)

            globals_["verify"] = namespace["verify"] = verify
        proxy = SimpleNamespace(run_path=load)
        globals_["runpy"] = namespace["runpy"] = proxy
        return namespace

    previous = load(HISTORICAL)
    captured = StringIO()
    arguments = sys.argv
    try:
        sys.argv = [str(HISTORICAL)]
        with redirect_stdout(captured):
            previous["main"]()
    finally:
        sys.argv = arguments
    return json.loads(captured.getvalue())


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.parse_args()
    previous = historical_delivery()
    # Fresh, unmodified helpers check current bytes without historical overrides.
    base = runpy.run_path(str(BASE_VERIFIER))
    manifest = base["read"](MANIFEST)
    base["require_flags"](manifest, EXPECTED_FLAGS, "Stage4 implementation start")
    records = base["entries"](manifest)
    names = {record["path"] for record in records}
    if MANIFEST.relative_to(ROOT).as_posix() in names:
        raise RuntimeError("The current manifest must exclude itself")
    required = {
        Path(__file__).resolve().relative_to(ROOT).as_posix(),
        "docs/GEMINGA_COMMANDS.md",
        "docs/implementation/stage4/README.md",
    }
    required.update((PREVIOUS / name).relative_to(ROOT).as_posix() for name in ARCHIVED_NAMES)
    if not required.issubset(names):
        raise RuntimeError("Delivery omits preserved/current files: " + ", ".join(sorted(required - names)))
    for name in names:
        if name.startswith("pysnspd/") and not name.startswith("pysnspd/experimental/"):
            raise RuntimeError("Stage4 delivery may change only experimental solver modules: " + name)
    count = base["verify"](records, require_size=True)
    references = sum(base["verify"](manifest.get(key, []), require_size=True)
                     for key in ("references", "retained_references"))
    external = Path("/home/jdiaz/GEMINGA_COMMANDS.md")
    if sys.platform.startswith("linux") and external.exists():
        if base["sha"](external) != base["sha"](ROOT / "docs/GEMINGA_COMMANDS.md"):
            raise RuntimeError("External Geminga command notebook differs from the current repository")
    print(json.dumps(dict(verified=True, files=count, references=references,
        manifest_sha256=base["sha"](MANIFEST), historical_delivery=previous,
        **EXPECTED_FLAGS), indent=2))


if __name__ == "__main__":
    main()
