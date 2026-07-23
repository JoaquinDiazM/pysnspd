# pySNSPD

**Multiscale research code for superconducting nanowire single-photon detector (SNSPD) simulations.**

This README documents the production research workflow of `pySNSPD` used to
obtain the main thermal-photon result of the thesis work. The publication
baseline, open defects, and current validation record live in `status.md`.

```text
Publication baseline : 23ea557657d890b1c902f5962a669d3fb845fd93
Historical tag       : dirty-functional-thermal-photon-v1
Repository           : github.com/JoaquinDiazM/pysnspd
Status               : research prototype, production pipeline validated
Date                 : 2026-07-23
```

---

## 0. Honest status

`pySNSPD` is not yet a polished detector simulator. It is a thesis/research codebase that grew quickly while integrating several physical layers that are usually studied separately:

1. dirty-limit Usadel superconductivity;
2. quasiparticle density-of-states catalogues;
3. electron-phonon phase-space integrals;
4. projected power/energy catalogues;
5. finite-volume mesoscopic gTDGL/Poisson dynamics;
6. two-temperature electron/phonon thermal dynamics;
7. photon perturbation;
8. external circuit readout;
9. plotting and current-sweep diagnostics.

The code at this SHA is **dirty but functional**. It reproduces the central result needed for the thesis: a photon-induced transient in a NbN nanowire where the order parameter, thermal fields, electric potential, current redistribution and output voltage evolve consistently over picosecond scales.

The repository still contains historical naming, partially redundant diagnostics, exploratory plotting code and configuration files created during development. The goal of this README is not to pretend that the project is clean, but to explain what is currently reliable and how the present result was obtained.

---

## 1. Main result of this frozen version

The current working pipeline produces a coupled photon response for a short NbN nanowire device. The representative result includes:

- transient suppression and recovery of the superconducting order parameter `Delta_s`;
- RF current response `I_RF`;
- central TDGL voltage `V_TDGL^center`;
- external output signal `V_out`;
- spatial maps of `|Delta|/Delta0`, electrostatic potential `phi`, gauge-invariant momentum magnitude `|q|`, electron temperature `T_e` and phonon temperature `T_ph`;
- time-resolved diagnostics at selected times after photon injection.

---

## 2. What this repository is trying to do

The purpose of `pySNSPD` is to build a reproducible multiscale framework for SNSPD modelling, especially for NbN devices. The code is designed around physical SI units and explicit intermediate catalogues. It is not a purely dimensionless TDGL toy model.

The central modelling idea is:

```text
microscopic material information
  -> mesoscopic superconducting dynamics
  -> thermal electron/phonon response
  -> electrical detector output
```

The thesis contribution is not a single new closed-form equation. The contribution is the integration of multiple layers into one auditable numerical workflow:

- Usadel-like superconducting material calibration;
- Allmaras/Vodolazov-inspired gTDGL dynamics;
- Simon/MIT-style electron-phonon spectral data and kernels;
- finite-volume current continuity and Poisson projection;
- two-temperature local thermal dynamics;
- photon perturbation and circuit-level observables.

---

## 3. Physical model, in implementation language

### 3.1 Material, geometry and calibration philosophy

The current working material is NbN, using physical parameters close to the thesis configuration:

```text
Tc                ~ 8.65 K
D_eff             ~ 1.58e-4 m^2/s   # calibrated effective diffusion constant for Ic = 38.8 uA
sigma_n           ~ 4.2e5 S/m
lambda_L          ~ 5.4e-7 m
thickness         ~ 7 nm
width             ~ 120 nm
```

The diffusion constant used in the functional runs should be read as an **effective calibrated dirty-limit diffusion constant**, not merely as a raw value copied from literature. During the development of the main result, `D_eff` was tuned so that the mesoscopic superconducting model reproduces the expected experimental cutoff/depairing current scale for the nominal NbN strip.

This calibration matters because `D_eff` enters the superconducting coherence scale, the Usadel pair-breaking scale, the current tables and the final critical-current scale of the gTDGL/Poisson state. In practice, the repo treats `D_eff` as part of the material-calibration layer connecting microscopic material assumptions to the experimentally relevant bias-current window.

The code assumes a thin superconducting strip and works with a 2D mesh representing the nanowire plane. All numerical quantities are stored and passed in SI units unless explicitly labelled otherwise. The calibration does not change this convention: the repo still evolves real physical quantities, but with a calibrated material parameter set.

### 3.2 Usadel, DOS and current catalogues

The PRE-run stage builds reusable superconducting catalogues. These include equilibrium gap information, quasiparticle density-of-states information and current-related tables over physical axes such as temperature, order-parameter amplitude and gauge-invariant momentum.

Conceptually, this layer provides quantities such as

```text
Delta_eq(T, q; D_eff)
rho(E; Delta, q, D_eff)
j_s(T, Delta, q; D_eff)
```

which are later queried by the stationary and dynamic solvers.

The current catalogue is not only a diagnostic output. It is part of the bridge between the calibrated material model and the mesoscopic solver: the stationary and photon stages inherit the current scale implied by the calibrated `D_eff`, `Tc`, `sigma_n`, geometry and gap model.

### 3.3 Electron-phonon phase space and power/energy catalogues

Earlier versions focused mainly on phase-space integrals. The current functional version goes further: it builds and uses projected **power/energy catalogues** suitable for thermal dynamics.

This is one of the main implementation changes relative to the original plan. Instead of recomputing expensive microscopic integrals inside the dynamic solver, expensive spectral quantities are precomputed and compressed into interpolation-ready catalogues.

The practical benefit is that the photon/thermal run can evolve local fields without repeatedly solving the full microscopic electron-phonon problem at every node and time step.

### 3.4 Mesh and finite-volume layer

The spatial domain is represented by a physical 2D mesh. The mesh stores nodes, triangles, edges and boundary tags. The solver uses finite-volume-like operators to evaluate gradients, currents, divergence diagnostics and boundary conditions.

This layer is important because the photon response is not treated as a purely 0D box. The plots show real spatial structure across the nanowire width and along the device.

### 3.5 Stationary gTDGL/Poisson state

Before injecting a photon, the code constructs a stationary superconducting state at the selected bias current. This state provides:

```text
Delta(r)
theta(r)
phi(r)
q(r)
j_s(r)
j_n(r)
```

The stationary solver uses a gTDGL/Allmaras-like backend coupled to a Poisson/current-continuity projection. Its material scale is inherited from the calibrated PRE-run catalogues. Therefore, the stationary state is not an arbitrary mathematical initial condition: it is the calibrated superconducting operating point from which the photon perturbation starts.

### 3.6 Photon, thermal dynamics and circuit response

The photon run evolves the coupled response after local energy deposition. The active dynamic variables include, schematically:

```text
Delta(r,t)
phi(r,t)
q(r,t)
T_e(r,t)
T_ph(r,t)
I_RF(t)
V_TDGL(t)
V_out(t)
```

The resulting diagnostic plots connect the local mesoscopic dynamics to circuit-level observables. In the frozen functional version, this final layer should be interpreted as the result of the full calibrated chain:

```text
calibrated NbN material set
  -> Usadel/DOS/current catalogues
  -> power/energy thermal catalogues
  -> stationary gTDGL/Poisson operating point
  -> photon/thermal/circuit transient
```

---

## 4. Official workflow at this SHA

The dirty functional workflow has three main pipeline stages plus plotting scripts.

### 4.1 PRE-run

Script:

```bash
python pipelines/01_prerun_template.py --help
```

Role:

- build/load mesh;
- build Usadel/DOS catalogue;
- build phase-space catalogue;
- build thermal power/energy interpolation objects;
- store raw reusable outputs under the external `big_data_root`.

### 4.2 SS-run

Script:

```bash
python pipelines/02_ss_run_template.py --help
```

Role:

- load PRE-run outputs;
- create the superconducting stationary seed;
- relax the stationary gTDGL/Poisson state;
- write stationary fields and diagnostics.

### 4.3 PHOTON-run

Script:

```bash
python pipelines/03_photon_run_template.py --help
```

Role:

- load PRE-run catalogues;
- load the SS stationary state;
- inject photon-like localized energy;
- evolve coupled thermal/gTDGL/circuit variables;
- write dynamic histories and raw fields.

### 4.4 Plot pipelines

Scripts:

```bash
python plot_pipelines/01_plot_prerun.py --help
python plot_pipelines/02_plot_ss_run.py --help
python plot_pipelines/03_plot_photon_run.py --help
python plot_pipelines/Z2_current_sweep_analysis.py --help
```

Role:

- regenerate plots without rerunning the expensive physics;
- inspect PRE, SS and PHOTON outputs;
- analyze current sweeps and overcritical branches.

---

## 5. Representative commands

These commands are written for the development machine used during the thesis work. Paths and run names may need adjustment on another system.

### 5.1 Checkout the frozen functional version

```bash
cd ~/pysnspd
git fetch --all --tags
git checkout dirty-functional-thermal-photon-v1
git rev-parse HEAD
```

Expected SHA:

```text
9aeeab333dd01952c94bc28286d028ccef1d7445
```

### 5.2 Basic validation

```bash
cd ~/pysnspd
python -m compileall -q pysnspd pipelines plot_pipelines tests
python -m pytest -q
```

At the frozen functional state, the local development run had the full test suite passing.

### 5.3 Representative photon plotting command

This command regenerates the main photon plots from an already completed PRE/PHOTON run:

```bash
cd ~/pysnspd
python plot_pipelines/03_plot_photon_run.py \
  --config configs/geminga_local_v3.yaml \
  --run-name photon_1p6eV_center_sigma12nm_I33uA_500ps_01 \
  --pre-run-name pre_oe6_v3_ultra_L360nm_mesh4p0nm_smooth50_js81T101D121Q_phase200T31D41Q2400W_power200Tph_01 \
  --scalar-times-ps 25 50 75 100 125 500 \
  --center-width-nm 100.0 \
  --dpi 640
```

This is a plotting command, not a full recomputation of the photon dynamics.

---

## 6. Data policy

Large outputs are intentionally kept outside git. The repository should not store raw `.npz` simulation outputs, large catalogues, logs or generated high-resolution plots.

The development configuration uses an external data root similar to:

```text
/home/jdiaz/scratch/big_data
```

Expected layout:

```text
big_data_root/
  raw/
    <run_name>/
      pre/
      ss/
      photon/
  plots/
    <run_name>/
  logs/
  catalogs/
    simon_2025/
      nbn-a2f-ph.dat
  tmp/
```

The same `run_name` should be used to link raw outputs, diagnostic plots and logs.

---

## 7. External data

The electron-phonon spectral input is external material data and is not committed to this repository. The local development path is:

```text
/home/jdiaz/scratch/big_data/catalogs/simon_2025/nbn-a2f-ph.dat
```

Expected text-file header:

```text
#E (THz) a^2F PhDOS (st/THz)
```

The code converts the frequency axis to an energy axis internally. Keep the raw external data outside the git repository.

---

## 8. Repository structure after the Week 1 cleanup

The top-level structure is approximately:

```text
configs/
docs/
pipelines/
  00_configure_project.py
  01_prerun_template.py
  02_ss_run_template.py
  03_photon_run_template.py
plot_pipelines/
  01_plot_prerun.py
  02_plot_ss_run.py
  03_plot_photon_run.py
  Z2_current_sweep_analysis.py
pysnspd/
  config.py
  analysis/
  circuit/
  excitation/
  io/
  kinetic/
  mesh/
  gtdgl/
  plotting/
  solver/
  thermal/
  usadel/
tests/
tools/
pyproject.toml
README.md
LICENSE
```

Finite-volume infrastructure now lives entirely under `mesh/`; stationary and
transient orchestration lives under `solver/`; thermal evolution and photon
deposition have their own packages. The governing placement and deletion rules
are documented in `docs/ARCHITECTURE_POLICY.md`.

---

## 9. Development objectives and status

The internal development was organized as objective blocks. The exact numbering evolved during the semester, but the current conceptual status is:

```text
[done] OE0  Project skeleton, configuration and external-data policy
[done] OE1  Material parameters, run manager and reproducible IO
[done] OE2  Physical 2D mesh, edge table and diagnostics
[done] OE3  Usadel/DOS catalogues and superconducting current tables
[done] OE4  Electron-phonon phase-space catalogues
[done] OE5  Thermal power/energy catalogues and Debye/Eliashberg diagnostics
[done] OE6  Stationary seed, plotting infrastructure and SS diagnostics
[done] OE7  Stationary gTDGL/Poisson solver and current-sweep analysis
[done] OE8  Photon, thermal dynamics and circuit-level output diagnostics
[*]    OE9  Calibration, parameter sweeps, validation and physical interpretation
```

The star intentionally remains on the last OE. Calibration and parameter sweeps are open-ended: there will always be better material parameters, better meshes, better bias points, better photon models and better validation datasets. The purpose of the star is to mark the active frontier, not to imply that the earlier layers are useless or incomplete.

---

## 10. What is reliable right now

At this SHA, the following parts are considered reliable enough for thesis-level use, with the normal caution required for a research prototype:

- the PRE/SS/PHOTON workflow concept;
- the use of physical SI units across the repository;
- the external-data policy using `big_data_root`;
- the mesh and finite-volume diagnostic philosophy;
- the Usadel/DOS catalogue approach;
- the power/energy catalogue approach for thermal dynamics;
- the stationary gTDGL/Poisson state as a photon initial condition;
- the photon plotting pipeline for the main dynamic figures;
- the current-sweep plotting approach for comparing bias regimes.

The plots produced from the current pipeline are suitable as the basis for thesis figures, after final captioning, discussion and consistency checks.

---

## 11. What is still dirty

This version still needs a full audit. Known issues include:

- some names still refer to earlier development objectives;
- some plotting modules contain thesis-specific choices;
- some configuration files are machine-specific;
- the dependency list in `pyproject.toml` is not yet a complete user-facing install recipe;
- some commands assume the development environment `snspd` and the `geminga` data layout;
- some diagnostic files are exploratory and should later be separated from the official workflow;
- the README, appendix and thesis text still need to be synchronized with the actual implementation;
- future users need clearer minimal examples and smaller test datasets.

This is expected. The current tag exists precisely to preserve the functional state before the repository is cleaned.

---

## 12. What should not be claimed yet

This repository should not yet be advertised as:

- a production-ready SNSPD simulator;
- a complete replacement for experimental calibration;
- a universal model for all SNSPD materials and geometries;
- a fully validated commercial detector-design tool;
- a polished package for non-expert installation.

The honest claim is stronger and more defensible:

> `pySNSPD` is a functional research framework that demonstrates a coupled micro/meso/thermal/circuit SNSPD simulation workflow in physical units, suitable for thesis-level analysis and future refinement.

---

## 13. Suggested next repository tasks

After this dirty functional tag, the next steps should be done in this order:

1. update the thesis appendix so that every equation corresponds to actual code;
2. write the main thesis narrative around the achieved coupled result;
3. audit every file and classify it as official, diagnostic, legacy or removable;
4. rename/reorganize source modules and pipelines;
5. create a clean public README for new users;
6. add a small reproducible demo dataset or lightweight smoke example;
7. improve installation instructions and dependency tracking;
8. document citation expectations and thesis references.

Do not refactor the full repository before the functional result has been fully documented.

---

## 14. Thesis positioning

The scientific value of this work is in the integrated modelling workflow. The project combines ideas from superconductivity, nonequilibrium kinetics, numerical PDEs, thermal modelling and circuit readout.

For the thesis, the main message should be:

```text
This work builds local computational capacity for modelling SNSPDs from material physics to detector-level signals, using an auditable Python framework in physical units.
```

This is relevant because SNSPDs are important for quantum communication, astronomical instrumentation, low-light sensing and future photonic technologies. A local framework like this helps Chilean research groups study detector physics without relying only on black-box tools or external codebases.

---

## 15. Minimal citation note

The final thesis/manuscript should cite the relevant literature for:

- dirty-limit Usadel superconductivity;
- gTDGL/KWT and Allmaras/Vodolazov SNSPD modelling;
- electron-phonon kinetic kernels and the Simon/MIT NbN spectral data;
- SNSPD detector physics and circuit readout;
- numerical finite-volume/mesh methods where appropriate.

This README intentionally avoids pretending that the code is independently derived from first principles without external literature. The implementation is a synthesis of published physical models and thesis-specific numerical engineering.

---

## 16. License

See `LICENSE` in the repository.

---

## 17. Maintainer

Joaquin Andres Diaz Monge  
Universidad de Chile, Departamento de Ingenieria Electrica  
GitHub: `JoaquinDiazM`

