# Implementation plan

This document is intentionally short. The detailed architectural explanation lives in `README.md`.

Recommended first implementation targets:

1. `pysnspd.config`
2. `pysnspd.io.manager`
3. `pysnspd.mesh.delaunay`
4. `pysnspd.mesh.edges`
5. `pysnspd.usadel.catalog`
6. `pysnspd.kinetic.phase_space`
7. `pipelines/01_prerun_template.py`

Do not implement PHOTON-run before the SS-run diagnostics are reliable.
