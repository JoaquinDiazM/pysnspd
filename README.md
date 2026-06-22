# pySNSPD

`pySNSPD` is a modular research code for building a multiscale simulation framework for NbN superconducting nanowire single-photon detectors (SNSPDs).

The repository contains the source code developed for an Electrical Engineering thesis. The goal is not to implement a monolithic detector simulator, but to build a reproducible pipeline where each physical layer can be tested before being coupled to the next one.

The model is organized around:

1. Dirty-limit Usadel material calibration.
2. Quasiparticle density-of-states catalogues.
3. Electron–phonon kinetic phase-space catalogues.
4. Projected electron–phonon powers.
5. Mesoscopic gTDGL dynamics.
6. External readout circuit.
7. PRE-run, SS-run and PHOTON-run workflows.

Raw simulation outputs are intentionally kept outside the repository.

---

## 1. Physical idea

The framework separates the SNSPD problem into reusable physical layers.

### 1.1 Material and Usadel layer

The material block calibrates the dirty-limit diffusion coefficient $D$ from a target critical current $I_c$, the wire dimensions, $\sigma_n$, $T_c$ and the bias temperature.

The Usadel block provides:

$$
\Delta_{\rm eq}(q,T),
\qquad
s_n(q,T),
\qquad
I_s(q,T),
\qquad
\rho(E;|\Delta|,q).
$$

The current relation used for calibration is the dirty-limit Matsubara expression

$$
j_s(q,T)
=
\frac{2\pi k_B T}{e}
\sigma_n q
\sum_n s_n^2(q,T).
$$

The DOS catalogue is then constructed from the retarded dirty-limit Usadel quartic

$$
(\Gamma_q c-i z)^2(1-c^2)-|\Delta|^2c^2=0,
\qquad
\rho(E)=\operatorname{Re}c(E),
$$

with

$$
\Gamma_q=\frac{\hbar Dq^2}{2}.
$$

### 1.2 Kinetic QP–phonon layer

The kinetic layer follows the Simon/MIT microscopic structure and projects the expensive electronic integrals into reusable phase-space catalogues:

$$
\mathcal J_S(\Omega;T_e,|\Delta|,q),
$$

$$
\mathcal J_R(\Omega;T_e,|\Delta|,q).
$$

These objects are not yet full powers. They do not include $\alpha^2F(\Omega)$, phonon DOS, $T_{ph}$, escape, diffusion or thermal evolution. They are expensive superconducting electronic phase-space integrals that can be reused later.

### 1.3 Projected power layer

The projected electron–phonon powers are built from the phase-space catalogues and the external Eliashberg material data:

$$
P_{ep}^{S}
=
\frac{8\pi N(0)}{\hbar}
\int d\Omega\,
\alpha^2F(\Omega)\Omega
\left[n_e(\Omega,T_e)-n_{ph}(\Omega,T_{ph})\right]
\mathcal J_S(\Omega;T_e,|\Delta|,q),
$$

$$
P_{ep}^{R}
=
\frac{4\pi N(0)}{\hbar}
\int d\Omega\,
\alpha^2F(\Omega)\Omega
\left[n_e(\Omega,T_e)-n_{ph}(\Omega,T_{ph})\right]
\mathcal J_R(\Omega;T_e,|\Delta|,q).
$$

The sign convention is:

$$
P_{ep}>0
$$

means energy leaves the electronic system and enters the phonon system.

The recombination channel $P_{ep}^{R}$ is a superconducting channel. It is not an additional normal-state power. In the OE5 implementation it is set to zero when

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
=
\frac{96\zeta(5)N(0)k_B^2}{\tau_0T_c^3}
\left(T_e^5-T_{ph}^5\right).
$$

The parameter $\tau_0$ is not the same as the linear electron–phonon relaxation time at $T_c$. The conversion used in the code is

$$
\tau_0
=
\frac{720\zeta(5)}{\pi^2}
\tau_{ep}(T_c).
$$

For $\tau_{ep}(T_c)=24.7\,{\rm ps}$,

$$
\tau_0\simeq1.868\,{\rm ns}.
$$

### 1.4 Thermal Usadel self-consistency

OE5 does not use an artificial fixed gap. The PRE-run now constructs a thermal Usadel grid

$$
\Delta_{\rm eq}(T_e,q)
$$

by recomputing the Matsubara Usadel self-consistency equation over a temperature and $q$ grid.

Projected powers are then evaluated as

$$
P_{ep}(T_e;\Delta_{\rm eq}(T_e,q),q).
$$

This is still not the final gTDGL-coupled PHOTON-run. In the final coupled model, the local fields

$$
\Delta(\mathbf r,t),
\qquad
q(\mathbf r,t)
$$

will come from the gTDGL sector. However, the thermal Usadel grid is the correct local-equilibrium audit for projected powers.

### 1.5 Mesoscopic gTDGL layer

The gTDGL layer will evolve the order parameter

$$
\Psi(\mathbf r,t)=R(\mathbf r,t)e^{i\phi(\mathbf r,t)},
$$

the electrostatic potential, supercurrent, normal current and current continuity. This layer determines whether a thermal perturbation turns into a dissipative detection event.

### 1.6 Circuit layer

The external circuit will couple the internal SNSPD voltage to a load and readout chain, producing

$$
I_{\rm SNSPD}(t),
\qquad
V_{\rm TDGL}(t),
\qquad
V_{\rm out}(t).
$$

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
    <run_name>/
  catalogs/
    simon_2025/
      nbn-a2f-ph.dat
      README.txt
  tmp/
```

The same `run_name` must be used for raw outputs and plots.

---

## 3. External material data

The NbN Eliashberg file used in OE5 is external material data and should live under:

```text
/home/jdiaz/scratch/big_data/catalogs/simon_2025/nbn-a2f-ph.dat
```

The source file header is:

```text
#E (THz)    a^2F    PhDOS (st/THz)
```

The code converts the THz axis to an energy axis using

$$
\Omega=hf.
$$

The file is not committed to the repository. The source must be documented in the external data folder and in code docstrings as Simon et al. 2025 material data.

---

## 4. Official workflow

The official workflow has three stages.

### 4.1 PRE-run

The PRE-run generates all expensive reusable objects:

* mesh;
* edges;
* boundary tags;
* Usadel/DOS catalogue;
* phase-space catalogue;
* thermal Usadel grid $\Delta_{\rm eq}(T_e,q)$;
* projected power diagnostics;
* manifests;
* plots.

The PRE-run is the natural place for parallelization, because catalogue points are independent. In contrast, SS-run and PHOTON-run are time evolutions where parallelizing explicit Euler over a spatial mesh would add technical complexity that is not currently useful.

### 4.2 SS-run

The SS-run will construct the stationary detector state before photon absorption. It should load PRE-run outputs and initialize:

$$
T_e=T_{ph}=T_{\rm bath},
$$

$$
|\Delta|\approx \Delta_{\rm eq}(T_{\rm bias},q_{\rm bias}),
$$

$$
\phi(\mathbf r)\approx qx,
$$

with current and potential fields compatible with the bias condition.

### 4.3 PHOTON-run

The PHOTON-run will load the PRE-run catalogues and the SS-run state, inject a phonon bubble or local energy source, and evolve:

$$
T_e(\mathbf r,t),
\quad
T_{ph}(\mathbf r,t),
\quad
\Psi(\mathbf r,t),
\quad
\phi_{\rm electric}(\mathbf r,t),
\quad
I_{\rm SNSPD}(t),
\quad
V_{\rm out}(t).
$$

---

## 5. Current repository structure

```text
configs/
pipelines/
  00_configure_project.py
  01_prerun_template.py
  02_oe5_power_debug.py        # obsolete after OE5 integration into PRE-run
pysnspd/
  config.py
  io/
    manager.py
  mesh/
    delaunay.py
    edges.py
  usadel/
    calibration.py
    catalog.py
    parameters.py
    solver.py
  kinetic/
    eliashberg.py
    phase_space.py
    powers.py
    thermal_usadel.py
  plotting/
    figures.py
    kinetic.py
tests/
```

---

## 6. PRE-run outputs

A complete PRE-run should write:

```text
raw/<run_name>/pre/
  mesh.npz
  edges.npz
  mesh_summary.yaml
  usadel_dos_catalog.npz
  usadel_dos_summary.yaml
  phase_space_catalog.npz
  phase_space_summary.yaml
  thermal_usadel_grid.npz
  oe5_power_catalog.npz
  oe5_power_summary.yaml
  manifest.yaml
```

and plots:

```text
plots/<run_name>/mesh/
  mesh_nodes_edges.png
  mesh_boundary_tags.png

plots/<run_name>/diagnostics/
  usadel_dos_slices.png
  usadel_calibration_sweep.png
  phase_space_slices.png
  eliashberg_spectrum.png
  thermal_usadel_gap_grid.png

plots/<run_name>/comparisons/
  electron_phonon_power_vs_Te_thermal_usadel.png
  electron_phonon_power_ratios_vs_Te.png
  electron_phonon_power_vs_thermal_usadel_current.png
  spectral_support_thermal_usadel_state.png
  low_energy_recombination_scattering_band.png
```

---

## 7. Basic command

A representative PRE-run command is:

```bash
cd ~/pysnspd

RUN_NAME=NbN_120nm_35uA_1064nm_PRE_OE5_final

python pipelines/01_prerun_template.py \
  --config configs/geminga_local.yaml \
  --run-name $RUN_NAME \
  --workers 16 \
  --jitter-fraction 0.05 \
  --boundary-guard-layers 1 \
  --eta-fraction 1.0e-3 \
  --gamma-max-fraction 0.80 \
  --energy-max-factor 30.0 \
  --phase-omega-max-meV 35.0 \
  --phase-n-Te 12 \
  --phase-n-delta 12 \
  --phase-n-q 12 \
  --phase-n-omega 480 \
  --eliashberg-dat /home/jdiaz/scratch/big_data/catalogs/simon_2025/nbn-a2f-ph.dat \
  --oe5-Te-min-K 0.9 \
  --oe5-Te-max-K 34.6 \
  --oe5-n-Te 180 \
  --oe5-q-scan-Te-K 0.9 6.92 8.65 17.3 34.6 \
  --oe5-n-q-thermal 140 \
  --oe5-n-matsubara-thermal 500 \
  --oe5-omega-max-meV 35.0 \
  --oe5-tau-ep-Tc-ps 24.7 \
  --oe5-support-min-delta-fraction 0.05
```

The flag

```bash
--workers 16
```

is intended for PRE-run catalogue work. It should not be interpreted as the strategy for future time evolution.

---

## 8. Development roadmap

The roadmap is organized as objectives. The star marks the current objective.

### OE1 — Configuration and file management

Status: done.

Implemented:

* required `big_data_root`;
* creation of `raw/`, `plots/`, `logs/`, `catalogs/`, `tmp/`;
* run-specific folders;
* manifests.

### OE2 — Mesh, edges, boundaries and plots

Status: done.

Implemented:

* protected rectangular Delaunay-like mesh;
* boundary guard layers without jitter;
* interior jitter;
* all nodes used;
* edge list;
* tags: `left`, `right`, `top`, `bottom`, `interior`;
* mesh plots.

### OE3 — Usadel/material block

Status: done.

Implemented:

* $D$ calibrated from target $I_c$;
* Matsubara Usadel sweep;
* $\Delta_{\rm eq}(q,T_{\rm bias})$;
* $I_s(q,T_{\rm bias})$;
* DOS catalogue $\rho(E;|\Delta|,q)$;
* calibration plots.

Reference result:

$$
D_{\rm fit}=1.5813\times10^{-4}\ {\rm m^2/s},
$$

$$
I_c=38.8\,\mu{\rm A}.
$$

### OE4 — Phase-space catalogues

Status: done.

Implemented:

$$
\mathcal J_S(\Omega;T_e,|\Delta|,q),
$$

$$
\mathcal J_R(\Omega;T_e,|\Delta|,q).
$$

The $\Omega$ axis is decoupled from the DOS energy axis, with finite-window diagnostics. For the current NbN tests, $\Omega_{\max}=35\,{\rm meV}$ is sufficient for the thermal projected powers.

### OE5 — Projected electron–phonon powers

Status: closing.

Implemented:

* Simon/MIT Eliashberg loader;
* $\alpha^2F(\Omega)$ from THz data converted to energy;
* $N(0)=\sigma_n/(2e^2D)$;
* projected powers $P_{ep}^{S}$ and $P_{ep}^{R}$;
* Debye/Vodolazov comparison with correct $\tau_0$;
* thermal Usadel grid $\Delta_{\rm eq}(T_e,q)$;
* $P_R=0$ when $\Delta_{\rm eq}=0$;
* normal Eliashberg reference;
* spectral-support diagnostics;
* low-energy gap-band comparison of scattering and recombination.

OE5 is considered complete once the PRE-run produces all OE5 plots and tests pass.

### OE6 — SS-run analytic seed and stationary initialization

Status: current.

Goal:

Build a physically consistent stationary initial condition for the nanowire before photon absorption.

Expected outputs:

* stationary initial $R(x,y)$;
* phase ramp $\phi(x,y)$;
* current density fields;
* potential field close to zero;
* compatibility with $I_{\rm bias}$;
* validation plots for current conservation and boundary behavior.

### OE7 ★ — Stationary gTDGL relaxation

Goal:

Relax the analytic SS seed with the gTDGL solver until a numerically clean stationary state is obtained.

Expected outputs:

* converged $|\Psi|$;
* converged phase;
* current continuity diagnostics;
* absence of spurious voltage;
* stationary metadata.

### OE8 — Electrical consistency validation

Goal:

Validate supercurrent, normal current, current crowding, $\nabla\cdot\mathbf j$, voltage and boundary conditions.

### OE9 — Phonon bubble

Goal:

Implement the photon-induced phonon bubble as a controlled energy deposition profile.

### OE10 — Decoupled thermal PHOTON-run

Goal:

Evolve $T_e$ and $T_{ph}$ with projected powers, diffusion and escape, without yet evolving gTDGL.

### OE11 — Coupled thermal + gTDGL PHOTON-run

Goal:

Couple $T_e,T_{ph}$ to $\Psi,\phi,V$ and test whether a local perturbation becomes dissipative.

### OE12 — Circuit coupling

Goal:

Couple the nanowire voltage to an external readout model and generate $I_{\rm SNSPD}(t)$, $V_{\rm TDGL}(t)$, $V_{\rm out}(t)$.

### OE13 — Comparison and publication plots

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
8. The README star should move only after tests, pipeline and outputs are verified.