# pySNSPD

`pySNSPD` is a modular research code for building a multiscale simulation framework for NbN superconducting nanowire single-photon detectors (SNSPDs).

The repository contains the source code developed for an Electrical Engineering thesis. The goal is not to provide a monolithic detector simulator, but to build a reproducible pipeline where each physical layer can be audited before being coupled to the next one.

The active model is organized around:

1. dirty-limit Usadel material calibration;
2. quasiparticle density-of-states catalogues;
3. electron–phonon kinetic phase-space catalogues;
4. projected electron–phonon power diagnostics;
5. stationary mesoscopic gTDGL/Poisson relaxation;
6. future thermal, photon and circuit coupling.

Raw simulation outputs are intentionally kept outside the repository.

---

## 1. Physical model

### 1.1 Material and Usadel layer

The material block calibrates the dirty-limit diffusion coefficient (D) from a target critical current (I_c), the wire geometry, (\sigma_n), (T_c) and the bias temperature.

The Usadel block provides local-equilibrium superconducting information:

$$
\Delta_{\rm eq}(q,T), \qquad
s_n(q,T), \qquad
j_s(q,T), \qquad
\rho(E;|\Delta|,q).
$$

The dirty-limit Matsubara current relation used for calibration is

$$
j_s(q,T)
========

\frac{2\pi k_B T}{|e|\hbar}
\sigma_n q
\sum_{n=0}^{\infty} s_n^2(q,T).
$$

The density-of-states catalogue is built from the retarded dirty-limit Usadel equation

$$
(\Gamma_q c-i z)^2(1-c^2)-|\Delta|^2c^2=0,
\qquad
\rho(E)=\operatorname{Re}c(E),
$$

with

$$
\Gamma_q=\frac{\hbar Dq^2}{2}.
$$

### 1.2 Kinetic electron–phonon layer

The kinetic layer follows the Simon/MIT microscopic structure and projects the expensive electronic integrals into reusable phase-space catalogues:

$$
\mathcal J_S(\Omega;T_e,|\Delta|,q),
$$

$$
\mathcal J_R(\Omega;T_e,|\Delta|,q).
$$

These catalogues are not full thermal powers by themselves. They are superconducting electronic phase-space factors that can later be combined with external Eliashberg material data, phonon occupation factors and local thermal fields.

### 1.3 Projected power diagnostics

The projected electron–phonon powers are constructed from the phase-space catalogues and the external Eliashberg spectrum:

$$
P_{ep}^{S}
==========

\frac{8\pi N(0)}{\hbar}
\int d\Omega,
\alpha^2F(\Omega)\Omega
\left[
n_e(\Omega,T_e)-n_{ph}(\Omega,T_{ph})
\right]
\mathcal J_S(\Omega;T_e,|\Delta|,q),
$$

$$
P_{ep}^{R}
==========

\frac{4\pi N(0)}{\hbar}
\int d\Omega,
\alpha^2F(\Omega)\Omega
\left[
n_e(\Omega,T_e)-n_{ph}(\Omega,T_{ph})
\right]
\mathcal J_R(\Omega;T_e,|\Delta|,q).
$$

The sign convention is:

$$
P_{ep}>0
$$

meaning energy leaves the electronic system and enters the phonon system.

The recombination channel (P_{ep}^{R}) is a superconducting channel. It is not an additional normal-state power. In the thermal Usadel audit it is set to zero when

$$
\Delta_{\rm eq}(T_e,q)=0.
$$

The normal-state comparison is made against

$$
P_{\Delta=0}^{\rm Eliashberg}
$$

and against the Vodolazov/Allmaras Debye reference

$$
P_D
===

\frac{96\zeta(5)N(0)k_B^2}{\tau_0T_c^3}
\left(T_e^5-T_{ph}^5\right).
$$

The parameter (\tau_0) is not the same as the linear electron–phonon relaxation time at (T_c). The conversion used in the code is

$$
\tau_0
======

\frac{720\zeta(5)}{\pi^2}
\tau_{ep}(T_c).
$$

For (\tau_{ep}(T_c)=24.7,{\rm ps}),

$$
\tau_0\simeq1.868,{\rm ns}.
$$

### 1.4 Mesoscopic gTDGL/Poisson layer

The stationary mesoscopic layer evolves the complex order parameter

$$
\Psi(\mathbf r,t)
=================

R(\mathbf r,t)e^{i\theta(\mathbf r,t)},
$$

the electrostatic potential, the supercurrent, the normal current and the current-continuity diagnostics.

The current stationary backend is the promoted pyTDGL-like finite-volume core, adapted to pySNSPD while keeping SI physical units at the repository level. The active SS-run loads PRE-run objects and relaxes a stationary state using:

* the mesh and edge table;
* the strict 3D Matsubara Usadel supercurrent table;
* Allmaras/gTDGL material factors;
* metallic left/right terminal conditions;
* finite-volume current-continuity diagnostics.

The practical SS result now used as baseline is:

* a short physical device, roughly twice the wire width, is sufficient to obtain a useful superconducting bulk;
* first phase-slip-line formation appears around (38)–(39,\mu{\rm A});
* multiple PSL states appear at higher current, with about three PSLs near (50,\mu{\rm A});
* the strongly overcritical branch recovers an Ohmic-like behavior around (60)–(70,\mu{\rm A}).

This validated stationary map is the reference before adding thermal and circuit coupling.

### 1.5 Future thermal, photon and circuit coupling

The future coupled model should evolve

$$
T_e(\mathbf r,t), \qquad
T_{ph}(\mathbf r,t), \qquad
\Psi(\mathbf r,t), \qquad
\phi(\mathbf r,t),
$$

and then couple the nanowire voltage to an external circuit to obtain

$$
I_{\rm SNSPD}(t), \qquad
V_{\rm TDGL}(t), \qquad
V_{\rm out}(t).
$$

The next major physical target is not photon injection yet. The next target is to verify that the fully coupled no-photon stationary state reproduces the already validated SS gTDGL branch.

---

## 2. Data policy

Raw data, catalogues, logs and plots must not be stored inside the git repository.

The config must define an external data root:

```yaml
project:
  big_data_root: /home/jdiaz/scratch/big_data
```

The expected layout is:

```text
big_data_root/
  raw/
    <run_name>/
      pre/
      ss/
      photon/
  plots/
    <run_name>/
      mesh/
      diagnostics/
      comparisons/
      figures/
  logs/
  catalogs/
    simon_2025/
      nbn-a2f-ph.dat
      README.txt
  tmp/
```

The same `run_name` links raw outputs, plots and logs.

---

## 3. External material data

The NbN Eliashberg file used by the kinetic diagnostics is external material data and should live under:

```text
/home/jdiaz/scratch/big_data/catalogs/simon_2025/nbn-a2f-ph.dat
```

The source file header is expected to be:

```text
#E (THz) a^2F PhDOS (st/THz)
```

The code converts the THz axis to an energy axis using

$$
\Omega=hf.
$$

The file is not committed to the repository.

---

## 4. Official workflow

The official workflow has three stages.

### 4.1 PRE-run

The PRE-run generates reusable expensive objects:

* mesh;
* edge table;
* boundary tags;
* Usadel/DOS catalogue;
* strict 3D Matsubara supercurrent table;
* superconducting phase-space catalogue;
* PRE diagnostic plots;
* manifests.

The PRE-run is the natural place for parallelization because catalogue points are independent.

### 4.2 SS-run

The SS-run constructs and relaxes the stationary detector state before photon absorption.

It loads PRE-run outputs and initializes:

$$
T_e=T_{ph}=T_{\rm bath},
$$

$$
|\Delta|\approx \Delta_{\rm eq}(T_{\rm bias},q_{\rm bias}),
$$

$$
\theta(\mathbf r)\approx qx.
$$

It then relaxes the gTDGL/Poisson system and writes stationary fields, histories, diagnostics and plots.

### 4.3 PHOTON-run

The PHOTON-run is still a future physical target.

It should eventually load the PRE-run catalogues and the SS-run state, inject a photon-induced perturbation, and evolve the coupled thermal, gTDGL and circuit variables.

The current repository should not pretend that this layer is implemented before the coupled no-photon stationary state is validated.

---

## 5. Current repository structure

```text
configs/
pipelines/
  00_configure_project.py
  01_prerun_template.py
  02_ss_run_template.py
  03_photon_run_template.py

plot_pipelines/
  01_plot_ss_run.py

pysnspd/
  config.py
  io/
    manager.py
  mesh/
    delaunay.py
    edges.py
    pytdgl_like.py
    quality.py
  usadel/
    calibration.py
    catalog.py
    parameters.py
    solver.py
    supercurrent_table.py
  kinetic/
    eliashberg.py
    phase_space.py
    powers.py
    thermal_usadel.py
  gtdgl/
    adapter.py
    allmaras.py
    currents.py
    device.py
    diagnostics.py
    finite_volume/
    geometry.py
    material.py
    operators.py
    options.py
    seed.py
    solver.py
    ss_targets.py
    state.py
    state_io.py
    tdgl_compat.py
    tdgl_operators.py
    usadel_current.py
  plotting/
    figures.py
    kinetic.py
    pre_diagnostics.py
    ss_figures.py
    ss_run.py

tests/
```

Legacy orchestration, circuit placeholders, validation placeholders, utility placeholders and obsolete plotting modules were removed in Audit 01.

---

## 6. PRE-run outputs

A useful PRE-run writes:

```text
raw/<run_name>/pre/
  mesh.npz
  edges.npz
  mesh_summary.yaml
  usadel_dos_catalog.npz
  usadel_dos_summary.yaml
  phase_space_catalog.npz
  phase_space_summary.yaml
  manifest.yaml
  plots_diagnostics/
```

The Usadel NPZ also stores the strict 3D Matsubara current table:

```text
js_A_m2[Te, delta, q]
Te_axis_K
delta_axis_J
q_axis_m_inv
```

This table is required by the active SS `usadel-poisson` current law.

---

## 7. Representative commands

### PRE-run

```bash
cd ~/pysnspd

RUN_NAME=NbN_120nm_35uA_1064nm_pre

python pipelines/01_prerun_template.py \
  --config configs/geminga_local.yaml \
  --run-name "$RUN_NAME" \
  --workers 16 \
  --eta-fraction 1.0e-3 \
  --gamma-max-fraction 0.80 \
  --energy-max-factor 30.0 \
  --phase-omega-max-meV 35.0 \
  --phase-n-Te 12 \
  --phase-n-delta 12 \
  --phase-n-q 12 \
  --phase-n-omega 480 \
  --js-table-n-Te 3
```

### SS-run

```bash
cd ~/pysnspd

PRE_RUN=NbN_120nm_35uA_1064nm_pre
SS_RUN=NbN_120nm_35uA_1064nm_ss_35uA

python pipelines/02_ss_run_template.py \
  --config configs/geminga_local.yaml \
  --run-name "$SS_RUN" \
  --pre-run-name "$PRE_RUN" \
  --ss-target-current-uA 35 \
  --ss-time-ps 20 \
  --ss-dt-fs 0.30 \
  --ss-snapshots 8 \
  --ss-progress \
  --dpi 480
```

### SS current sweep

```bash
cd ~/pysnspd

PRE_RUN=NbN_120nm_35uA_1064nm_pre
SS_RUN=ss_sweep_Dgtdgl100_tau_v2_long_base34uA_60ps_01

python pipelines/02_ss_run_template.py \
  --config configs/geminga_local.yaml \
  --run-name "$SS_RUN" \
  --pre-run-name "$PRE_RUN" \
  --ss-target-current-uA 34 \
  --extra-currents-uA +3 +6 +9 +12 +16 +26 +36 \
  --ss-sweep-workers 4 \
  --ss-time-ps 60 \
  --ss-dt-fs 0.30 \
  --ss-snapshots 8 \
  --ss-progress \
  --dpi 480
```

### SS plot pipeline

```bash
cd ~/pysnspd

python plot_pipelines/01_plot_ss_run.py \
  --config configs/geminga_local.yaml \
  --run-name "$SS_RUN" \
  --pre-run-name "$PRE_RUN" \
  --dpi 480
```

---

## 8. Development roadmap

The roadmap is organized by physical milestones. Early objectives were compacted because they became small infrastructure layers rather than full standalone physical stages.

### OE1 — Repository, configuration and data layout

Status: done.

Implemented:

* external `big_data_root`;
* run-specific `raw/`, `plots/`, `logs/`, `catalogs/`, `tmp/` policy;
* config validation;
* manifest writing;
* process-safe storage checks.

### OE2 — PRE-run foundation

Status: done.

Implemented:

* pyTDGL-style rectangular finite-volume mesh;
* boundary and edge table;
* dirty-limit Usadel calibration;
* DOS catalogue;
* strict 3D Matsubara supercurrent table;
* superconducting phase-space catalogue;
* PRE diagnostic plots.

This objective compactly replaces the older split between mesh-only, Usadel-only and phase-space-only objectives.

### OE3 — Projected electron–phonon power audit

Status: done.

Implemented:

* Simon/MIT Eliashberg loader;
* (\alpha^2F(\Omega)) and phonon-DOS support;
* (N(0)=\sigma_n/(2e^2D));
* projected scattering and recombination powers;
* Debye/Vodolazov comparison with correct (\tau_0);
* thermal Usadel grid (\Delta_{\rm eq}(T_e,q));
* normal-state Eliashberg reference;
* spectral-support diagnostics.

This objective is now considered a diagnostic/catalogue layer, not the final thermal solver.

### OE4 — Stationary gTDGL/Poisson backend

Status: done.

Implemented:

* promoted pyTDGL-like finite-volume backend;
* SI material adapter;
* metallic terminal treatment;
* Allmaras/gTDGL current-mismatch correction;
* strict `usadel-poisson` supercurrent law;
* adaptive stationary relaxation;
* stationary state, history and diagnostic outputs;
* SS plot pipeline.

### OE5 — Stationary current sweep and PSL map

Status: done.

Validated:

* stable superconducting bulk in a short wire of length about twice the width;
* first PSL formation around (38)–(39,\mu{\rm A});
* multiple PSL states at higher current, including roughly three PSLs near (50,\mu{\rm A});
* Ohmic-like overcritical behavior recovered around (60)–(70,\mu{\rm A}).

This is the current reference stationary map.

### OE6 ★ — Coupled no-photon stationary state

Status: next.

Goal:

Validate that the complete system, including thermal variables and circuit variables, has a no-photon stationary state consistent with the already validated SS gTDGL branch.

Expected checks:

* (T_e\approx T_{ph}\approx T_{\rm bath}) in subcritical superconducting states;
* no artificial Joule heating in a zero-voltage superconducting state;
* same qualitative branch structure as the OE5 stationary map;
* consistent (I_{\rm SNSPD}), (V_{\rm TDGL}), (|\Delta|), (T_e), (T_{ph}) and circuit state.

### OE7 — Photon energy deposition

Status: planned.

Goal:

Implement the photon-induced perturbation as a controlled local energy deposition model.

This objective should not be started before the coupled no-photon stationary state is validated.

### OE8 — Decoupled thermal PHOTON-run

Status: planned.

Goal:

Evolve (T_e) and (T_{ph}) with projected powers, diffusion and escape, without yet evolving the full gTDGL/circuit feedback.

### OE9 — Coupled thermal + gTDGL PHOTON-run

Status: planned.

Goal:

Couple (T_e,T_{ph}) to (\Psi,\phi,V) and test whether a local perturbation becomes dissipative.

### OE10 — Circuit coupling and readout pulse

Status: planned.

Goal:

Couple the nanowire voltage to an external readout model and generate

$$
I_{\rm SNSPD}(t), \qquad
V_{\rm TDGL}(t), \qquad
V_{\rm out}(t).
$$

### OE11 — Comparison and publication plots

Status: planned.

Goal:

Generate reproducible plots comparing bias current, photon energy, geometry, thermal assumptions and microscopic catalogues.

---

## 9. Design rules

1. No heavy raw data inside the repository.
2. Every heavy object must be generated in PRE-run and registered in a manifest.
3. SS-run and PHOTON-run should load catalogues, not rebuild them.
4. PRE-run can and should be parallelized.
5. Time evolution is not the first target for parallelization.
6. The same `run_name` links raw data, plots and logs.
7. Physical shortcuts must be explicit in metadata and plots.
8. No placeholder function should silently return `0`.
9. The README star should move only after tests, smoke tests and representative outputs are verified.
