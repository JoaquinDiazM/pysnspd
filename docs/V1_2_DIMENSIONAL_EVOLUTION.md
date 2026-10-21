# pysnspd v1.2.3: driven lead diagnostics and exact-state continuation

The v1.2.2 evolution RHS, geometry, native kinetic/material/circuit kernels and
all nonlinear/local/global energy tolerances remain unchanged. This revision
corrects a domain-monitor interpretation, not an identified evolution-solver bug.
The archived2D and1D histories can be retained when the authenticated current
state independently passes the new monitor. A new folder preserves the source.

## Continuous physics and measured diagnosis

A hotbelt is not an insulated material boundary. Electronic transport and
normal-current/KWT dissipation can change native occupations outside it.
Furthermore unloading the current reduces depairing: the cold stable bulk gap
increases. An exterior gap above its original DC value is not a global
suppression. Negative excess internal-energy density relative to the initial
biased DC can include release of its current-carrying energy.

For the two v1.2.2 terminal states the free collar gap/initialDC was1.0711
(1550nm) and1.0934(775nm), while Te reached2.1743/2.3671K and Tph stayed
near the bath. Static native-RHS decomposition finds positive condensate
relaxation power and transport, rather than a uniform injected photon source.
At x-x0=475nm(1550nm), native power densities in W/m3 were approximately
KWT2.40e11, transport5.36e10, normal Joule1.73e10 and reversible1.43e10.
The independent instantaneous total-power residual was below7e-21W.
This is evidence inside our effective closure, not microscopic validation of
low-temperature KWT/constant-normal-conductivity heating.

Primary context: Vodolazov2017, equations30-32 and the discussion of diffusion,
Joule heating and normal-domain growth:
https://arxiv.org/pdf/1611.06060 .
Neither that reference nor a small energy residual certifies our closure's
late quantitative temperature, decay time or finite-domain independence.

## Explicit diagnostic correction

For retained energy Eret, use the already-existing driven activity ledger
$$
A(t)=E_{\rm ret}+\int_0^t|I_s V_{\rm dev}|\,dt,\qquad
E_{\rm col}^{\rm abs}=\sum_{i\in{\rm col}}
\left[\sum_k\epsilon_{ik}|p_{ik}-f_{\rm FD}(E_{ik},T_b)|
+\sum_l c_{il}\Omega_l|n_{il}-n_{B,l}|\right].
$$
Here f_FD denotes the electronic Fermi occupation at bath temperature;
epsilon is the actual node spectral energy weight, and c Omega is the native
phonon energy weight. This finite-volume definition is the same unsigned
native collar metric as1.2.2. The new active normalization is
$$
\varepsilon_{\rm col}=E_{\rm col}^{\rm abs}/A(t),\qquad
r_{\Delta,\rm col}
=\min_{i\in{\rm col}}\frac{|\Delta_i|}
{s_i a_b(I_s,T_b)} .
$$
The fixed matching s_i reproduces the discrete initialDC and a_b is the same
finite-link stable bulk branch already used at the exterior cuts.

Keep the0.5K warning,0.1% energy margin,90% gap margin and2ps persistence.
A stop requires current-relative gap suppression below90%, OR a thermal
warning AND collar energy above0.1% of the driven activity. The old
photon-only ratio and initial-DC gap ratio remain visible diagnostics.
At the old stops, epsilon_col was6.6883e-6/1.1027e-5 and rDelta was
0.96484/0.97385(1550/775nm). The photon-only ratios were0.11637/0.11179%.

This normalization is an explicit engineering choice, changed because the
system is driven far beyond the original photon energy. It is not a rigorous
exterior error bound. The independent2% integrated absolute-cut/activity
guard, current-relative gap guard, support/capacity/voltage/polarity/energy
guards remain active. Real front escape can still stop a subsequent run;
domain-length convergence through15ns is NOT established by this patch.

## Data integrity and scope

A source-frozen v1.2.2 thermal-only diagnostic stop may be forked byte-for-byte
into a new scratch directory. Authenticate the published revisions, unchanged
RHS/material/circuit hashes, exact inputs, checkpoint/DC and copied HDF5.
Retain every prior state and energy ledger, preserve source provenance and
apply the new monitor BEFORE advancing. Do not replay overlap, reset a bath,
redeposit a photon, replace a spectrum, modify Is or remove recorded failures.
Only the old monitor persistence timer is reset for the changed criterion.
Actual equations,1um domain and15ns conditional target remain unchanged.

The late cold approximation and its independent10ps qualification are unchanged.
A3h1Dsession budget remains a resumable runtime checkpoint, not completed reset.
No full15ns validation or<3h end-to-end claim is made.
The report labels both collar normalizations and the instantaneous-current
bulk reference; old reports continue to show their historical active basis.

---

# pysnspd v1.2.2: cold continuation cuts and gated late recovery

The hot-stage equations remain the neutral effective thermal BCS/KWT functional,
native electron/phonon populations and actual three-state circuit. Keep300nm2D,
1000nm/201nodes1D, the material/photon inputs, positive0.1mV trigger,0.5ps
confirmation and crossing+2ps handoff. The final target is15ns;100ps is ONLY
the no-trigger censoring ceiling. No absorbed photon is deposited twice.

## Diagnosis and exterior assumption

The v1.2.1 handoff was admitted. Its two continuations stopped after about2ps
when their outer collars reached1.62K, despite a gap above the original DC gap.
The235nm suppressed segment did not establish escape from the1000nm domain.
The added tails followed the instantaneous lower current but terminal gap and
populations remained fixed to the original higher-current DC, creating a
terminal relaxation/heating layer.

The new explicit exterior hypothesis is an isothermal stable superconducting
continuation at the instantaneous detector current. For $d=\Delta/(k_BT_c)$,
bulk amplitude $a_b(I_s,T_b)$ and a discrete DC matching factor $s_\partial$,
$$
d_\partial=s_\partial a_b(I_s,T_b)e^{i\theta_\partial},\qquad
p_\partial=f_{\rm FD}(|d_\partial|,T_b),\qquad n_\partial=n_B(T_b).
$$
The prescribed radial rate is $s_\partial a'_b(I_s)\dot I_s$; the FD rate uses
the chain rule. Its work is explicitly included, not hidden in an energy offset:
$$
P_\partial^{\rm prescribed}
=2\operatorname{Re}\sum_\partial F_d^*\dot d_\partial^{\rm radial}
+\sum_{\partial,k}\epsilon_k\dot p_{\partial,k},\qquad
P_{\rm res,out}=P_{\rm raw,cuts}-P_\partial^{\rm prescribed}.
$$
The full thermal chain rule, signed fluxes, support, capacity and independent
energy audits remain. The monotone finite-link bulk table is computed from
our functional, not fitted to pulse timing, and never extrapolated.

A0.5K collar excursion remains a visible warning. Its thermal truncation stop
also requires unsigned native collar excess energy above0.1% retained photon
energy. Gap suppression below90%DC and cumulative2% absolute-cut exchange
remain independent stops. The margins are declared engineering choices,
not experimental material constants.

## Native solver

Solve midpoint equations in amplitude/phase above $|d|=10^{-3}$; smaller
condensates retain Cartesian coordinates. This is not a gap/population floor.
A zero-mean temporal gauge removes only an unobservable common phase:
$$
\widetilde\phi_i=\phi_i-\frac{\sum_jV_j\phi_j}{\sum_jV_j},\qquad
\dot a_i=\operatorname{Re}(e^{-i\theta_i}\dot d_i^{\rm material}),\qquad
\dot\theta_i=\operatorname{Im}(\dot d_i^{\rm material}/d_i)
-\frac{2e}{\hbar}\widetilde\phi_i .
$$
Newton and population-response tests contain the same rank-one gauge
response. Electrostatic output stays right-cut referenced. Use the actual
converged polar midpoint, not the Cartesian endpoint average.

Step doubling controls amplitude, physical superconducting face currents,
device voltage, unsigned native population energy and every circuit state.
The complex-gap error remains diagnostic. Keep rtol/atol, full nonlinear
thermal-response acceptance, local energy2.5e-6, global energy1%, support and
positive capacity. Repeated gradients use exact full d,T,alpha byte keys,
with eight cache entries and no rounded or approximate state substitution.

## Additional late cold approximation

Native1D remains official while the film is hot or significantly nonequilibrium.
Only after100ps can the following full-domain bounds qualify a cold candidate:

- gap amplitude within1% of the instantaneous-current cold branch;
- unsigned native energy and film projection exchange each<=1% retained photon;
- superconducting face-current discrepancy<=0.5% DC bias;
- full device-voltage agreement<=1uV;
- omitted normal-current fraction<=1e-3.

An independent actual-circuit cold shadow must additionally agree with native1D
for10ps, within1uV Vout and0.5% branch currents. A failed qualification leaves
native1D running. A hot/latched wire cannot qualify; no preset-time replacement
is allowed.

After admission, the approximation neglects the small native residual and
reconstructs the isothermal BCS film at the bath. The signed projection bath
exchange is recorded. The actual circuit at activation is preserved. Late
quasiparticle-tail spectroscopy and quantitative retrigger claims are outside
this approximation.

With resolved length $L_1$, bulk phase gradient $q_b(I_s)$ and native cold-film
energy $U_b(I_s)$,
$$
L_f(I_s)=\frac{\hbar L_1}{2e}\frac{dq_b}{dI_s},\qquad
\dot I_s=\frac{v_c+R_L(I_b-I_s)}{L_{\rm ext}+L_f(I_s)}.
$$
The other two CM equations and physical circuit parameters stay unchanged.
The frozen exterior retains all added96nH. The resolved film is replaced
by its cold constitutive branch; it is not also evolved or counted twice.
For $L_U=(1/I_s)dU_b/dI_s$,
$$
P_{b,\rm out}=(L_f-L_U)I_s\dot I_s,\qquad
\frac{d}{dt}(U_b+E_{\rm CM})
=V_bI_b-R_bI_b^2-R_L(I_b-I_s)^2-P_{b,\rm out}.
$$
This power is independently integrated and compared against actual
film+circuit energy, not inferred from its observed residual.
DOP853 uses scaled rtol2e-9/atol2e-11, accepted steps<=50ps and internal
steps<=10ps. The omitted normal Joule estimate
$\int V_f^2/R_n\,dt$ is independently accumulated and bounded by1% retained
photon energy.

This conditional separation of film cooling and electrical recovery has
primary physical context in [Kerman et al.](https://arxiv.org/pdf/1302.2852)
and [Annunziata et al.](https://arxiv.org/pdf/1008.0895).
Neither supplies our admission margins or validates this particular NbN closure.
HDF5/manifest data label the activation time and regime_id=1 for reconstructed
late fields; native data have regime_id=0. Reports label the transition.
Use --no-cold-recovery for a fresh full-native control.

## Data reuse and measured performance

Recover the original valid v1.2.1 2D pair, not its later contaminated fixed-cut
1D outer layer. Preserve the20 native material/kinetic/circuit source hashes,
DC, checkpoint and complete original archive. In a new scratch directory replay
only1D overlap from the stored2D snapshot immediately before crossing (at most
0.1ps earlier), record every comparison and never rerun the accepted2D history.

Each case retains4SMT threads on2exclusive physical cores, BLAS/OpenMP1; a pair
uses8worker threads/4cores plus coordinator and reserves>=2whole cores globally.
The optional --max-1d-wall-hours budget saves a resumable
WALL_BUDGET_CHECKPOINTED, not physical completion. A later --resume starts a
new wall session without resetting physical time.

Actual saved-photon controls advanced3.2ps beyond handoff and passed the old2ps
stop. A nontrivial15ns cold-circuit recovery control took about4s after system
preparation, with relative energy defect<2e-10. This is not a full15ns photon
validation or a guarantee of3h end-to-end: admission can dominate. Persistent
hot states stay native and can reach wall/domain/latching guards. Exact scopes
and commands are in the machine validation record and command notebook.
