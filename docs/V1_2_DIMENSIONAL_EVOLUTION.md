# pysnspd v1.2.1: longer full1D and physical stop monitors

The default full1D domain is1000nm (201nodes), centered on the300nm2D lens with350nm added on each side. The new --length-1d-nm input can extend up to the5um active wire. This supersedes the same-domain300nm scope of1.2.0; its historical source/validation remains in that tag.

## Evolution and initial state

Keep the effective neutral thermal BCS/KWT/native-kinetic closure and actual three-state circuit. An independent1D shadow begins at the first accepted threshold crossing. Only initially unobserved tails use a cold instantaneous-current stable-branch hypothesis, with phase matched to the core and the prepared DC profile. The core retains its nonequilibrium native populations. This is a declared adiabatic exterior hypothesis, not recovery of an unobserved history.

At admitted crossing+2ps, reproject the last official2D core and inherit the already evolved native shadow tails; do not reset them to a cold FD/BE state. Conservative partial-volume weights account for cells straddling both interfaces. Bound the small native activity reconciliation by the same energy margins. Compare gap fields modulo a global U(1) phase, since domains fix potential at different right cuts; preserve all relative gradients.

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

The five stable Matematica PDFs and their sources describe the implementation, new P50-P56 monitors and S33-S35 assumptions. Use docs/GEMINGA_COMMANDS_PRESENT.md. Six cases use24SMT threads on12physical worker cores, preserving current jobs and at least2whole cores. Fresh1.2.1 cases are required.
