# pysnspd v1.2.1: longer full1D and physical stop monitors

The default full1D domain is1000nm (201nodes), centered on the300nm2D lens with350nm added on each side. The new --length-1d-nm input can extend up to the5um active wire. This supersedes the same-domain300nm scope of1.2.0; its historical source/validation remains in that tag.

## Evolution and initial state

Keep the effective neutral thermal BCS/KWT/native-kinetic closure and actual three-state circuit. An independent1D shadow begins at the first accepted threshold crossing. Only initially unobserved tails use a cold instantaneous-current stable-branch hypothesis, with phase matched to the core and the prepared DC profile. The core retains its nonequilibrium native populations. This is a declared adiabatic exterior hypothesis, not recovery of an unobserved history.

At admitted crossing+2ps, reproject the last official2D core and inherit the already evolved native shadow tails; do not reset them to a cold FD/BE state. Conservative partial-volume weights account for cells straddling both interfaces. Bound the small native activity reconciliation by the same energy margins. Compare the final projected complex gap modulo one global U(1) phase. For independent long/short-domain overlap, retain that complex norm diagnostically and gate amplitude plus the actual phase-dependent supercurrent, populations, energy and circuit. Bank-phase evolution across a near-normal belt and different exterior cuts must not be mislabeled as a transverse-amplitude discrepancy; phase-gradient errors remain guarded through the closure's supercurrent and energy.

Seed1D DC with exact finite-link bulk dispersion to avoid accumulated phase error across a long chain. Bath/DC force/current/nV tolerances remain unchanged.

## Inductance and energy

Freeze Lext2=96nH+Lprime_DC*(5um-L2) and Lext1=96nH+Lprime_DC*(5um-L1) from the same original continuous DC. Never reduce the added96nH. Copy Ib,Is,vc exactly. Resolved-film excess energy receives deltaEfilm2+0.5*(Lext2-Lext1)*(Is^2-IDC^2), keeping total excess film+circuit energy continuous. Record only the constant total-DC reference shift; no duplicated inductance or invented ongoing heat correction.

## Stop monitors

- Retain2% absolute cut-energy/activity and1% independent global-energy guards.
- Watch free25nm outer collars:0.5K electron/phonon moment rise or gap below90%DC for2ps gives POSSIBLE_DOMAIN_ERROR.
- Accepted current limits derive from1.1Idep; voltage limits from Rnormal, the actual bias and load. They are explicit engineering ceilings, not measured hardware ratings.
- Global branch reversal is guarded; moderate bipolar Vout is allowed for capacitive undershoot.
- POSSIBLE_LATCHING requires a>=20nm resistive segment and1ns electrical/thermal plateau after5ns. It is a possible physical state, not proof of a numerical defect.
- RECOVERY_CONFIRMED can stop before15ns only after1ns of jointly close full-circuit, gap, native population and electron/phonon-moment diagnostics, after2ns.
- Store accepted native state and explicit monitor reason at termination. Default15ns remains an upper target, not a validated reset or runtime promise.

Primary context: [Kerman et al.2013,Fig.1](https://arxiv.org/pdf/1302.2852), [Annunziata et al.2010](https://arxiv.org/abs/1008.0895). The own restricted functional is not replaced by the Allmaras hotbelt equations.

## Validation and execution

Short tests include conservative core volume and total-energy transfer, cold-tail initialization, finite-link DC, inherited nonthermal shadow-tail populations, gauge-invariant overlap, support/capacity/temporal energy acceptance, exact restart, physical stop distinctions and shared-resource refusal. Read-only saved-hot-state projection to201nodes and one short coupled1D advance are controls, not a15ns photon validation.

The five stable Matematica PDFs and their sources describe the implementation, new P50-P56 monitors and S33-S35 assumptions. Use docs/GEMINGA_COMMANDS_PRESENT.md. Six cases use24SMT threads on12physical worker cores, preserving current jobs and at least2whole cores. Changes of physical parameters or geometry require fresh cases. The known1.2.1 join failure has the authenticated same-physics continuation described below.


## 2026-10-20 handoff correction (same package version1.2.1)

The failed775/1550nm pair completed its official2D detection windows. Its raw final projection errors were0.51986%/0.50980%, but a single exterior phase alignment also created artificial join gradients. Each disconnected exterior component must instead be attached to its own interface phase. These constant component rotations preserve each tail's internal phase differences and native p,n; this is an explicit initial-state joining operation, not a space-dependent gauge transformation of an already connected wire. Choose interface nodes by grid index: floating coordinate inequalities can exclude a nominally coincident boundary node.

The conservative occupation matrix R still satisfies V1 R=V2. A nodal condensate must instead be reconstructed at its physical sites from the remap's physical longitudinal centroids and actual cut values. Cell averages at a partial interface are not nodal values. The small existing activity correction then conserves total film+circuit excess energy within its original envelope.

For the actual long/short independent overlap, let d2 be the2D gap and Ld1 the lifted1D gap. In addition to the recorded globally aligned complex norm, admission uses

    eps_amp^2 = sum_i V_i (|d2_i|-|Ld1_i|)^2 / sum_i V_i |dDC_i|^2
    eps_Is^2 = sum_e (Rs/g_e) (Is2_e-IsLift1_e)^2 / [IDC^2 Rs L2/W].

Is_e is the superconducting FACE CURRENT in amperes, from the same functional variation; g_e is the dimensionless positive face conductance factor. Both thresholds are3%. Zero-conductance faces have zero weight. The final instantaneous projection still requires the3% complex norm as well. Resolved vortices, native population discrepancy, current pattern,2uV readout,0.5% circuit current,1% shadow energy,0.5% raw energy and1% reconciliation gates remain. This changes the interpretation of the independent phase comparison; it does not claim the old complex20% discrepancy was small.

The fixed saved-state projections give raw energy0.23780%/0.23699%, complex gap0.37063%/0.36940%, and supercurrent1.07363%/1.06748% (775/1550). Conservative saved-shadow time brackets pass with amplitude<=1.804%, supercurrent<=2.931%, and stored simultaneous film-energy<=0.560%. No tolerance was widened. These figures are implementation/finite-overlap controls, not15ns physical validation.

Recovery is a new-directory, source-hash-authenticated fork of the specific original1.2.1 pre-handoff failure. Every material/kinetic/circuit/geometry source stays identical. Preserve the full original archive/DC/checkpoint and re-audit saved native overlap; never redeposit the photon. Exact2D snapshots are compared against both neighboring saved shadow states for spatial norms; recorded simultaneous scalars audit energy/work/circuit, avoiding a false defect from real Joule work between different times.

The runner also preserves its pre-existing proposed physical step after a successful step shortened solely to land on a snapshot target. A plotting remainder must not seed a long sequence of microscopic future steps. Rejected steps still halve and retain all physical/energy tolerances. Campaign progress is now mirrored to the terminal.
