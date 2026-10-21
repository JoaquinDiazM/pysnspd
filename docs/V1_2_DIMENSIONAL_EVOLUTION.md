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
