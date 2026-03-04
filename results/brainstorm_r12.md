# 10 Most Important Actionable Items for the sBootstrap SUSY Lagrangian Program

**Date**: 2026-03-04
**Scope**: Framework as described: SU(3) SQCD N_f=N_c=3 / NMSSM / Koide / bloom / CKM.
**Context**: Complete Lagrangian assembled (Rounds 10–11). North Star = SUSY Lagrangian
with Yukawa couplings from Higgs mechanism.

---

## Item 1: Baryon Stabilization — Fix the Global Minimum Problem

**The problem precisely**: The scalar potential V = V_F + V_D + V_soft has a global minimum
at M^i_j = 0 (all mesons vanish), X = 0, BB~ = -Lambda^6, H_u = -m_d/y_c, H_d = -m_s/y_b,
with V_min = m_u^2 = 4.67 MeV^2. The Seiberg seesaw vacuum (M_i = C/m_i, B = B~ = 0) has
V_seesaw ~ 10^15 MeV^2, which is 10^14 times higher. All mass predictions (m_c, m_b, m_h)
and the Koide structure live at the seesaw vacuum, not the true ground state (offdiag_vacuum.md).

**Why it matters**: Every prediction in the paper is computed at the seesaw vacuum. If this
vacuum is not even metastable, or if the tunneling lifetime to M=0 is shorter than the age
of the universe, the entire predictive framework is operating at the wrong vacuum and the
predictions are meaningless. A referee will identify this immediately.

**Specific calculation**: Add a baryon mass term W_B = m_B B B~ to the superpotential. This
lifts the baryonic flat direction, gives B and B~ masses of order m_B, and removes the M=0
escape route (the constraint det M - BB~ = Lambda^6 can no longer be satisfied with M=0 if
B is massive). Redo the vacuum analysis with this term and verify that the seesaw becomes
the unique minimum with V > 0. Show that m_B > f_pi (any value above the soft scale works).
If the baryon mass term is too ad hoc, compute instead the Coleman-Weinberg potential for
B along the baryonic flat direction and show that quantum corrections lift it.

Alternatively, compute the bounce action for vacuum decay from seesaw to M=0. Thin-wall
approximation: S_bounce = 27 pi^2 sigma^4 / (2 epsilon^3) with Delta V ~ 10^15 MeV^4 and
field displacement Delta phi ~ C/m_u ~ 4 x 10^5 MeV. Even a rough estimate that gives
S_bounce >> 1 would establish metastability.

**Difficulty**: 2/5 (analytic estimate is straightforward; exact calculation requires numerics)
**Impact**: 5/5 (foundation of the entire predictive framework)

---

## Item 2: FCNC — State the Correct Physical Identification and Compute the New Contribution

**The problem precisely**: fcnc_constraints.md shows that if off-diagonal mesons M^d_s, M^u_d,
M^u_s are new scalar particles at mass m_S = f_pi = 92 MeV, the K-K̄ mixing Wilson coefficient
C_1 = g^2/(2m_S^2) exceeds the UTfit bound 9 x 10^-13 GeV^-2 by a factor of 66 (at g=1),
rising to a factor of 3460 at NDA coupling. This appears to exclude the model outright.

**Why it matters**: If this analysis is correct, the framework is ruled out by decades-old
experimental data. It is the single most dangerous objection to publishability. The paper
currently has no answer (assembly_round11.md notes the FCNC is "reframed" but the
resolution is not explicit).

**The resolution and the specific calculation**: The off-diagonal mesons M^i_j (i≠j) ARE
the SM pseudoscalar mesons (pions, kaons), not new particles. The K^0 meson is M^d_s bar
(or equivalently, the kaon is the ds component of the meson matrix in the 3-flavor sector
u,d,s). Their K-K̄ contribution is already in the SM calculation. There is no additional
tree-level exchange from new particles.

The new FCNC that does exist must be computed: loop diagrams from the NMSSM coupling
lambda X (H_u H_d) and the Yukawa terms y_c M^d_d H_u, y_b M^s_s H_d. The effective
four-fermion operator (d̄s)(d̄s) from these loops has coefficient:

    C_1^{loop} ~ lambda^2 y_c^2 / (16 pi^2 m_X^2 m_{mesino}^2)

With lambda = 0.72, y_c = 1.03 x 10^-2, m_X = lambda|⟨X⟩|v^2 ~ 10^5 MeV (Higgsino mass),
m_mesino ~ |X_vac| M_u M_d ~ 10^-4 MeV (extremely small), this coefficient is either
negligible or the mesino is so light that it raises a different question (see Item 5).

Write down the meson-SM-particle dictionary explicitly, compute the genuine new FCNC
contribution from NMSSM loops, and show it satisfies the UTfit bound.

**Difficulty**: 3/5 (requires identifying the dictionary clearly; loop computation is standard)
**Impact**: 5/5 (without this, the paper appears excluded by existing data)

---

## Item 3: CKM Derivation — Analytic Perturbative Treatment at the Seesaw Vacuum

**The problem precisely**: The Oakes/GST relation sin(theta_C) = sqrt(m_d/m_s) gives
theta_C = 12.92 degrees (PDG: 13.04 degrees, deviation -0.9%). This is numerically accurate
but currently treated as a coincidence or an empirical observation, not a derived result.
The full numerical minimization (assembly_round10.md, 10D) finds theta_12 = 2.6 degrees,
inconsistent with the Oakes result. The discrepancy arises because the numerical code finds
the M=0 global minimum or a saddle point, not the physical seesaw vacuum.

**Why it matters**: The Cabibbo angle is the only CKM element the model has any chance of
predicting from first principles. If the derivation succeeds, it upgrades the Oakes connection
from "observed consistency" to "structural prediction." If it fails, the CKM connection must
be demoted from prediction to coincidence.

**Specific calculation**: With the baryonic flat direction stabilized (Item 1), expand
M^i_j = (C/m_i) delta^i_j + epsilon^i_j at the seesaw vacuum, treating epsilon as small.
The off-diagonal potential at second order in epsilon reads:

    V^{(2)} = sum_{i≠j} (f_pi^2 + W_{ij}^2 / something) |epsilon^i_j|^2 + cross terms

The tachyonic off-diagonal mass squared at the EWSB vacuum (F_X = -lambda v^2/2) is
(from assembly_round10.md): m^2_{ds} = f_pi^2 - 2 lambda^2 v^2 M_u / Lambda^2 < 0.
The condensate epsilon^d_s satisfies a gap equation. At leading order in the condensate,
the mixing angle between mass eigenstates is:

    theta_{12} ~ |epsilon^d_s| / M_s

Show that |epsilon^d_s| / M_s = sqrt(m_d/m_s) follows from the gap equation, recovering
the Oakes relation. The key ratio is determined by the tachyonic mass hierarchy:
m^2_{ds} : m^2_{us} : m^2_{ud} = M_u : M_d : M_s = 1/m_u : 1/m_d : 1/m_s (seesaw inverted).

This is a concrete 3-day calculation (3x3 matrix perturbation theory plus one gap equation).

**Difficulty**: 3/5 (second-order perturbation theory in a 3x3 matrix; algebra is tractable)
**Impact**: 5/5 (a first-principles derivation of the Cabibbo angle would be a major result)

---

## Item 4: mu-Problem — Explicit Statement and Literature Survey of Solutions

**The problem precisely**: After the canonical rescaling X → S = X Lambda^4,
the effective mu-term is mu_eff = lambda ⟨S⟩ = 0.72 x (-9.80 MeV) = -7.06 MeV
(mu_term_analysis.md). This is three orders of magnitude below the ~100 GeV needed for
viable electroweak symmetry breaking. The Higgs mass prediction m_h = lambda v/sqrt(2) = 125.4
GeV is independent of mu_eff (it comes from the NMSSM quartic), so the Higgs mass is
fine. But the physical Higgs spectrum (m_A, m_H^+, the mixing matrix) depends on mu_eff,
and EWSB requires B_mu ~ mu_eff^2 ~ (100 GeV)^2 not (7 MeV)^2.

**Why it matters**: Without a viable mu-source, the electroweak spectrum is wrong. The W
and Z masses are generated by v = 246 GeV (assumed), but the SUSY partner spectrum
(Higgsinos at mu_eff = 7 MeV) would produce unobserved light charged particles.

**Specific analysis**: Survey the three viable resolutions:

(A) Giudice-Masiero: K_GM = kappa (X H_u H_d + h.c.) / M_Pl generates mu_GM = kappa F_X / M_Pl.
With F_X = -lambda v^2/2 = -2.18 x 10^10 MeV^2 from the EWSB vacuum (assembly_round10.md),
mu_GM = kappa x 2.18 x 10^10 / (1.22 x 10^22) = 1.79 kappa x 10^-12 MeV. For mu_GM ~ 100 GeV
= 10^5 MeV, need kappa ~ 5.6 x 10^16. Clearly does not work; the F_X in this model is
suppressed by the meson sector, not by the electroweak scale.

(B) Different identification: X is NOT the NMSSM singlet. X generates the Koide constraint
via the Seiberg term. A SEPARATE singlet S with ⟨S⟩ ~ 100 GeV provides the mu-term.
This adds one parameter (the S VEV) but decouples the mu-problem from the Koide structure.
Cost: the elegant X = S identification is lost. The Higgs mass prediction remains (it comes
from the same coupling structure). Benefit: mu is a free parameter, solvable by any standard
NMSSM mechanism.

(C) EWSB is generated by the tachyonic condensation of off-diagonal mesons (Item 3), not
by the standard Higgs mechanism. In this scenario v is not 246 GeV; it is determined by
the tachyonic mass scale ~ sqrt(lambda^2 v^2 M_u / Lambda^2). This is a genuinely
nonstandard Higgs mechanism that would change the entire EWSB structure.

Write these three options down explicitly in the paper with their numerical consequences.
Choose one and state it as the framework's answer to the mu-problem.

**Difficulty**: 2/5 (this is largely a literature survey and numerical estimate; no new physics required)
**Impact**: 4/5 (referees will ask; honest acknowledgment is required)

---

## Item 5: Mesino Mass Scale — Phenomenological Bounds on the Lightest Charged Fermion

**The problem precisely**: The off-diagonal mesino masses at the seesaw vacuum are
extremely small: m_mesino(d-s) = |X_vac| M_u M_d = 4.94 x 10^-4 MeV = 0.494 keV
(fcnc_constraints.md). These are neutral or charged fermions (the charged mesinos M^u_d
and M^u_s have electric charges +1 and +1). A charged fermion at 0.5 keV would have
been detected in precision measurements of the hydrogen energy levels, beam-dump experiments,
or atomic physics since the 1960s.

**Why it matters**: If charged fermions at 0.5 keV exist in the model and are not identified
with known particles, the model is ruled out by atomic physics alone. This is a show-stopper
that needs to be addressed before any other calculation.

**Specific analysis**: The mesino masses at the seesaw vacuum scale as |X_vac| x M_i x M_j.
With X_vac ~ C/Lambda^6 ~ 10^-9 MeV^-4 and M_i ~ C/m_i ~ 10^5 MeV:

    m_mesino(u-d) = |X_vac| M_s ~ (10^-9)(10^4) ~ 10^-5 MeV  (electric charge +1)
    m_mesino(u-s) = |X_vac| M_d ~ (10^-9)(10^5) ~ 10^-4 MeV  (electric charge +1)
    m_mesino(d-s) = |X_vac| M_u ~ (10^-9)(10^5) ~ 10^-4 MeV  (electric charge 0)

These are sub-eV fermions with SM charges. The question is: are they the SM fermions?
Specifically, are the charged mesinos M^u_d identified with the muon (as suggested by
lepton_susy.md, which finds the pi-mu pair works), or are they new particles?

Compute the experimental bounds on new charged fermions below 1 MeV:
- Measurement of the electron magnetic moment (g-2) at 0.3 ppb precision: a new charged
  fermion at m ~ 10^-4 MeV contributes Delta a_e ~ alpha/(pi) (m_e/m_mesino)^2 ~ 10^7 alpha/pi.
  This is excluded by 10^7 sigma.
- The only resolution: the charged mesinos ARE the SM leptons (electron, muon), not new particles.
  The u-d mesino at 0.01 eV is identified with the electron neutrino; the off-diagonal M^u_d
  fermion is identified with the muon.

Verify or falsify the identification of mesinos with SM leptons by checking that the
mesino quantum numbers (electric charge, SU(2)_L representation, baryon number) match the
SM leptons.

**Difficulty**: 2/5 (quantum number comparison is straightforward; bound computation is literature)
**Impact**: 5/5 (a charged fermion at sub-eV mass is excluded unless it is the electron)

---

## Item 6: Lepton Sector — Specify the Sp(2) Model or Acknowledge the Gap

**The problem precisely**: Lepton masses are absent from the Lagrangian. The Kahler approach
is definitively ruled out (rank invariance; best Q = 0.499, not 2/3; complete_lagrangian.md
Section 9). The Sp(2) SQCD proposal from assembly_round11.md (Agent 11D) is the most
promising path: Sp(2) with N_f = 3 antisymmetric flavors produces 3 Pfaffian eigenvalues,
the O'Raifeartaigh vacuum at gv/m = sqrt(3), Q = 2/3 at the seed, and bloom rotation
gives the hierarchy. But the Sp(2) vacuum has not been verified, the mediation scale
M_* ~ 43.5 GeV is not derived, and the model is not written down.

**Why it matters**: The paper claims to explain fermion masses. Without a mechanism for
lepton masses, this claim is false. At minimum, the paper must state clearly that lepton
masses are not explained by the quark/Higgs sector of the model.

**Specific calculation**: For the Sp(2) route, write down the superpotential:
W_Sp = sum_{a<b} m_ab Q^a Q^b + X_L (Pf M_L - Lambda_L^3)
where M_L^{ab} = Q^a Q^b (antisymmetric, 6 entries for 4 flavors, or 3 entries for Sp(2)
with N_f = 3 fundamentals). Identify the Pfaffian constraint, find the ADS vacuum,
check gv/m = sqrt(3) condition, and compute the fermion mass eigenvalues. Compare to
(m_e, m_mu, m_tau) = (0.511, 105.7, 1777) MeV. The mediation scale (connecting the
Sp(2) sector to the lepton-Higgs coupling) must produce the overall lepton mass scale
v_0(lepton) ~ 1 MeV^(1/2), which requires a specific relation between Lambda_L and
the electroweak scale.

If this calculation fails or is too involved for the current paper, demote to: "The
model explains the Koide structure for leptons by postulating a parallel confining sector.
The overall lepton mass scale is an input (m_tau = 1776.86 MeV, with m_mu and m_e predicted
by Koide)." One extra input, zero claims about derivation.

**Difficulty**: 4/5 (full Sp(2) model); 1/5 (honest demotion)
**Impact**: 4/5 (a paper about SUSY fermion masses needs an answer to the lepton question)

---

## Item 7: Bloom Mechanism — Derive or Justify the Nonperturbative Origin

**The problem precisely**: The Q-preserving bloom is a pure delta-rotation in the Koide
parametrization (bloom_mechanism.md). The seed at delta = 3pi/4 transitions to delta ~ 3pi/4
+ 22 degrees = 157.2 degrees during the quark seed-to-full transition, with v_0 doubling
(v_0(full)/v_0(seed) = 2.0005). The delta-rotation is NOT driven by the perturbative
Coleman-Weinberg potential (m_CW ~ 0.073 MeV vs m_b = 4180 MeV, a discrepancy of 5
orders of magnitude). The bloom mechanism must be nonperturbative.

The explicit R-symmetry breaking analysis (bloom_mechanism.md, Result 3) shows that any
nonzero R-breaking lifts the goldstino mass as m_0 ~ epsilon^2/m but pushes Q away from
2/3 continuously. The seed is an unstable fixed point. No nonzero epsilon preserves the
Koide condition.

**Why it matters**: The v_0-doubling prediction (m_b = 4177 MeV, 0.1 sigma from PDG) is
among the most precise predictions in the model. Its physical basis is the bloom from seed
to full triple. If the bloom mechanism has no derivation — just the observation that the ratio
is 2.0005 — then the v_0-doubling is a numerical coincidence masquerading as a prediction.

**Specific calculation**: The bion Kähler potential K_bion = (zeta^2/Lambda^2)|sum_k s_k sqrt(M_k)|^2
provides a potential on the Koide manifold. Compute the minimum of K_bion restricted to the
Koide constraint (Q = 2/3, fixed v_0). The minimum is at the value of delta that minimizes
the bion potential. Show that this minimum corresponds to delta ~ 157.2 degrees (the full
quark triple) rather than delta ~ 135.2 degrees (the seed).

The bion potential on the constraint surface is:
V_bion(delta) = (zeta^2/Lambda^2) |s_s sqrt(M_s(delta)) + s_c sqrt(M_c(delta)) + s_b sqrt(M_b(delta))|^2

where M_k(delta) = v_0^2 (1 + sqrt(2) cos(2pi k/3 + delta))^2 and the signs s_k flip at
M_k = 0. Map the minima of this potential and check that the physical bloom corresponds to
a local minimum (or a transition between the seed minimum and the full-triple minimum driven
by the bion sign flip).

**Difficulty**: 3/5 (numerical minimization of a one-variable function; analytic approximation feasible)
**Impact**: 4/5 (transforms v_0-doubling from numerology to a derived consequence of bion dynamics)

---

## Item 8: Renormalization Scheme — On-Shell vs MS-bar for sin^2(theta_W)

**The problem precisely**: The de Vries angle R = (sqrt(19)-3)(sqrt(19)-sqrt(3))/16 = 0.22310
matches the on-shell value sin^2(theta_W)(on-shell) = 0.22306 to 0.13 sigma. But the MS-bar
value sin^2(theta_W)(MS-bar, M_Z) = 0.23122 is 36 sigma away. The prediction is scheme-dependent
at the level of 36 sigma difference. A referee will immediately ask: why is the on-shell
scheme the natural comparison for this Casimir formula?

**Why it matters**: If the scheme choice is arbitrary, the sin^2(theta_W) prediction is
not a prediction. If the on-shell scheme is physically motivated by the model, this needs
to be explained. This is the single most vulnerable number in the paper from a referee's
perspective.

**Specific analysis**: The de Vries angle arises from a Casimir energy equation for the
SO(32) adjoint 496, decomposed under SU(5) x SU(3). The eigenvalue condition is derived
from a spin-statistics sum that depends on on-shell degrees of freedom (physical polarizations).
The on-shell scheme counts physical polarizations; the MS-bar scheme includes unphysical
degrees of freedom. If the Casimir formula uses physical propagating states (not gauge-fixed),
then the on-shell scheme is the natural comparison.

Write this argument explicitly:
1. State that the Casimir formula involves only transverse (physical) polarizations.
2. Define sin^2(theta_W)(on-shell) = 1 - M_W^2/M_Z^2 = 0.22290 (from M_W = 80.377 GeV,
   M_Z = 91.188 GeV).
3. Note that the difference between on-shell and MS-bar at M_Z is a precisely computed
   electroweak radiative correction: Delta sin^2 theta_W = 0.23122 - 0.22290 = 0.00832,
   which depends on the top mass, Higgs mass, and alpha_s.
4. Ask: does the model predict any of these radiative corrections, or does it predict
   the tree-level on-shell value before electroweak corrections?

If the answer is "tree-level on-shell value," state it explicitly and note that radiative
corrections will shift the prediction by ~4%, within the current uncertainty.

**Difficulty**: 2/5 (this is a physics argument, not a new calculation)
**Impact**: 4/5 (the sin^2(theta_W) prediction is among the model's strongest; must be defended)

---

## Item 9: Top Quark as External — Derive the Decoupling Scale

**The problem precisely**: The top quark is treated as an elementary chiral superfield in
W_Yukawa = y_t H_u^0 t t̄, external to the SQCD sector (complete_lagrangian.md, Gap 2).
This is the central structural asymmetry: charm and bottom are composites (meson eigenstates),
but top is elementary. No derivation of this asymmetry is provided. The MEMORY.md note is
"top = elephant (non-composite)" but no threshold scale is computed.

**Why it matters**: The claim that charm and bottom are seesaw composites (M^d_d and M^s_s)
requires that the SQCD sector confines at a scale Lambda where charm and bottom can be
bound states but top cannot. The confining scale is Lambda ~ 300 MeV. But the charm quark
has m_c ~ 1270 MeV >> Lambda, so charm is not lighter than the confinement scale. If charm
is a composite at Lambda = 300 MeV, why not top (m_t ~ 172 GeV >> Lambda even more so)?

**Specific analysis**: The condition for a quark to be "composited" by the confining sector
is not simply m_q < Lambda; in Seiberg duality, all N_f quark flavors appear in the meson
matrix regardless of their mass. The distinction between composite and elementary should be:

(A) The top quark's mass m_t >> Lambda means the electric description (quark degrees of
    freedom) is appropriate for the top, while the magnetic description (meson degrees of
    freedom) is appropriate for (u,d,s) with m_q < Lambda.

(B) Specifically, for the N_f = N_c = 3 SQCD to confine and produce the meson matrix M^i_j,
    the quark masses m_i must satisfy m_i < Lambda for the confining phase to be valid.
    The charm quark (m_c = 1270 MeV >> Lambda = 300 MeV) violates this condition.

Compute the mass threshold at which a quark decouples from the confined description:
m_decouple = Lambda x (some function of N_f, N_c). Check whether m_c is above or below
this threshold. If charm is above the threshold, the entire identification of M^d_d with
the charm-generating condensate is in the wrong phase of the theory.

The resolution may be that the confinement scale Lambda is not 300 MeV but some other
value determined self-consistently from the masses m_i. Compute this self-consistency
condition: Lambda must be such that the lightest three quark flavors (u,d,s) are in the
confined phase (m_u, m_d, m_s < Lambda) while charm decouples (m_c > Lambda). This gives
300 MeV < Lambda < 1270 MeV, consistent with Lambda = 300 MeV = f_pi sqrt(4pi).

**Difficulty**: 3/5 (requires careful analysis of the Seiberg phase structure)
**Impact**: 3/5 (clarifies the physical basis of the composite vs elementary distinction)

---

## Item 10: Full 3x3 CKM Matrix — Beyond the Cabibbo Angle

**The problem precisely**: The current CKM analysis produces only theta_12 (the Cabibbo
angle) from the Oakes/GST relation sin(theta_C) = sqrt(m_d/m_s). The full CKM matrix
has four parameters: three angles theta_12, theta_13, theta_23 and one CP-violating phase
delta. The model has no prediction for theta_13, theta_23, or delta. Furthermore, the
parameter count shows that after the Koide constraints fix (m_c, m_b, m_t) and Oakes
fixes m_d/m_s, the remaining free parameters are m_u and m_d independently. The ratio
m_u/m_d (or equivalently m_u) connects to the other CKM angles through the Fritzsch
mass-matrix texture.

**Why it matters**: CKM is one of the deepest structures in the SM. A model that claims
to explain fermion masses should at least have an approach to CKM mixing beyond the
Cabibbo angle. The North Star (complete SUSY Lagrangian with Yukawa couplings) requires
a complete treatment of flavor mixing.

**Specific calculation**: The Fritzsch texture (1977) gives:

    V_us = sqrt(m_d/m_s) - sqrt(m_u/m_c) e^{i phi}
    V_cb = sqrt(m_s/m_b)
    V_ub = sqrt(m_u/m_t) e^{i delta}

Compute these with PDG values m_u = 2.16, m_d = 4.67, m_s = 93.4, m_c = 1270, m_b = 4180,
m_t = 172760 MeV:

    sqrt(m_s/m_b) = sqrt(93.4/4180) = 0.1495 (PDG: |V_cb| = 0.0413)  — fails by factor 3.6
    sqrt(m_u/m_t) = sqrt(2.16/172760) = 3.54 x 10^-3 (PDG: |V_ub| = 3.82 x 10^-3)  — matches 7%

The V_cb prediction fails. Check whether the extended Koide constraints (v_0-doubling
gives sqrt(m_b) = 3 sqrt(m_s) + sqrt(m_c)) affect the Fritzsch formula. With m_b from
v_0-doubling (4177 MeV) vs PDG (4180 MeV), the change in V_cb is negligible (0.07%).

Identify which of the four CKM parameters are predicted by the model and which remain
free. Compare to the PDG CKM matrix with current precision.

**Difficulty**: 2/5 (numerical evaluation of known formulas with the model's mass predictions)
**Impact**: 3/5 (gives a complete CKM picture; the failure of V_cb in Fritzsch must be explained)

---

## Summary Table

| # | Item | Difficulty | Impact | Core issue |
|---|------|------------|--------|------------|
| 1 | Baryon stabilization / vacuum lifetime | 2/5 | 5/5 | Seesaw is not the ground state |
| 2 | FCNC identification (mesons = SM mesons) + new FCNC | 3/5 | 5/5 | Model appears excluded by K-K̄ mixing |
| 3 | CKM: analytic derivation of Cabibbo from tachyonic condensation | 3/5 | 5/5 | Best quantitative prediction needs derivation |
| 4 | mu-problem: options and honest statement | 2/5 | 4/5 | EWSB spectrum is wrong without mu ~ 100 GeV |
| 5 | Mesino mass scale: charged fermions at sub-keV | 2/5 | 5/5 | Sub-eV charged fermions excluded unless = SM leptons |
| 6 | Lepton sector: Sp(2) or honest demotion | 1–4/5 | 4/5 | Lepton masses missing from Lagrangian |
| 7 | Bloom mechanism: bion potential derivation of delta | 3/5 | 4/5 | v_0-doubling needs physical basis |
| 8 | Renormalization scheme for sin^2(theta_W) | 2/5 | 4/5 | 36-sigma MS-bar discrepancy needs explanation |
| 9 | Top decoupling: derive the threshold scale | 3/5 | 3/5 | Composite vs elementary distinction is assumed |
| 10 | Full 3x3 CKM beyond Cabibbo | 2/5 | 3/5 | Complete CKM matrix needed for North Star |

---

## Priority Ordering for the North Star (Complete SUSY Lagrangian with Yukawa)

**Critical path** (must do before any other work):
1 → 5 → 2 (establish which vacuum is physical; identify SM particles correctly)

**Core predictions** (strongest results, most important to defend):
3 → 7 → 8 (Cabibbo derivation; bloom physics; sin^2 theta_W scheme)

**Completeness** (needed for a complete Lagrangian):
6 → 4 → 9 → 10 (lepton sector; mu-problem; top decoupling; full CKM)

The key observation: Items 1, 2, and 5 are logically prior to everything else.
Without knowing which vacuum is physical and which particles are which, the Lagrangian
is not yet written for the physical world. Items 3, 7, and 8 elevate the model's
strongest numerical coincidences to derived results. Items 4, 6, 9, 10 complete the
picture but can be stated as known open problems without invalidating the rest.

---

## What a Physically Motivated Referee Will Kill the Paper On

In decreasing order of lethality:

1. "Your charged mesinos at 5 x 10^-4 MeV are excluded by atomic physics. The electron
   g-2 is measured to 0.3 ppb; a new charged fermion at 0.5 keV contributes at the 10^7
   ppb level. What are these particles?" (Item 5)

2. "Your global minimum is M = 0 at V = (2.16 MeV)^2. Every prediction you make is at
   a metastable vacuum 10^14 times higher. Compute the tunneling rate or add a term that
   makes the seesaw the minimum." (Item 1)

3. "Off-diagonal scalars at 92 MeV violate K-K̄ mixing by a factor of 66 at g = 1.
   Your table in fcnc_constraints.md shows this clearly. Is the model excluded?" (Item 2)

4. "The mu-term from your NMSSM coupling is 7 MeV. The W boson is at 80 GeV. How does
   electroweak symmetry breaking work in your model?" (Item 4)

5. "sin^2(theta_W) in the MS-bar scheme is 0.23122, 36 sigma from your prediction.
   Why do you compare to the on-shell value?" (Item 8)

---

*Generated: 2026-03-04*
*Based on: complete_lagrangian.md, offdiag_vacuum.md, fcnc_constraints.md,
mu_term_analysis.md, lepton_alternatives.md, lepton_susy.md, assembly_round10.md,
assembly_round11.md, bloom_mechanism.md, v0_doubling.md, ckm_koide.md,
cabibbo_oakes.md, pauli_15_10bar.md, iss_cw_koide.md, round2_status.md*
