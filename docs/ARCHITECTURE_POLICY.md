# pySNSPD architecture and deletion policy

This policy applies to production library code under `pysnspd/`. The executable
roots that decide whether code is production-reachable are every Python entry
point under `pipelines/` and `plot_pipelines/`.

## Deletion rule

A function, class, method, or complete module may be deleted only when it is
not reachable from any production entry point.

Reachability includes:

- direct imports and calls;
- transitive calls through library functions;
- constructors, methods, properties, callbacks, and protocol methods used by
  a reachable object;
- imports executed inside functions;
- module-level execution caused by a reachable import.

Tests are evidence, but they are not production roots. A test whose only
purpose is to protect unreachable legacy behavior is removed with that
behavior. Deleted implementation is not copied to a `deprecated`, `legacy`, or
sandbox package; Git history is the archive.

The reproducible audit command is:

```bash
python tools/reachability_audit.py --output tmp/reports/reachability.json
```

The auditor combines AST import/call structure, Jedi project references, and a
conservative attribute-name fallback for dynamically typed method/property
access. Its report must still be reviewed before applying
`tools/prune_unreachable.py`.

## Placement rule

Modules are grouped by scientific responsibility:

- `usadel`: microscopic superconducting material catalogues;
- `kinetic`: electron-phonon spectra, phase space, and power catalogues;
- `mesh`: geometry, finite-volume data structures, and discrete operators;
- `gtdgl`: mesoscopic gTDGL physics and constitutive laws;
- `thermal`: electron/phonon temperature evolution;
- `excitation`: external perturbations such as photon deposition;
- `circuit`: readout and lumped electrical evolution;
- `solver`: orchestration, runtime options, stationary/transient evolution,
  targets, callbacks, history, and state persistence;
- `analysis`: derived run diagnostics;
- `plotting`: figure construction only;
- `io`: run/catalogue storage and retrieval.

Package `__init__.py` files remain intentionally small and do not recreate a
flat compatibility API. Production code imports the defining module directly,
which makes ownership and reachability explicit.

## File-size rule

Production library modules must not exceed 800 physical lines. A module below
roughly 100 lines is acceptable only when its responsibility is naturally
small, such as configuration, models, state IO, handlers, style, or package
initialization. Large cohesive implementations are split by responsibility,
not by arbitrary line ranges.

## Change acceptance

An architecture change is accepted only after:

1. `compileall` succeeds for the library, pipelines, plot pipelines, tests, and
   maintenance tools;
2. every pipeline entry point accepts `--help`;
3. the complete surviving pytest suite passes on Geminga;
4. the reachability audit reports no unresolved references and no unreachable
   production definitions;
5. no production library module exceeds 800 lines;
6. `status.md` records removed behavior, obsolete tests, and validation
   evidence.
