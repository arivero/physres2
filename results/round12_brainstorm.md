# Round 12 Brainstorm: Open Problems and Priorities

**Date**: 2026-03-04
**Scope**: What remains before a final paper; referee questions; single most important calculation.

---

## Model Status (what this brainstorm assumes as established)

The model is N=1 SUSY, confined SU(3) SQCD, N_f = N_c = 3 (u,d,s). Low-energy fields:
mesons M^i_j, baryons B, B-tilde, Lagrange multiplier X (= NMSSM singlet after rescaling),
Higgs H_u, H_d. The complete Lagrangian is assembled in complete_lagrangian.md.

**Verified predictions**: m_c = 1301 MeV (1.6 sigma), m_b = 4177 MeV (0.1 sigma),
m_h = 125.4 GeV (0.6 sigma), sin^2 theta_W = 0.22310 (0.4 sigma), V_us = 0.224 (0.9 sigma).

**Vacuum structure** (offdiag_vacuum.md, assembly_round10.md): The true global minimum is
M = 0 with BB-tilde = -Lambda^6. F_X = -lambda v^2 / 2 dominates once EWSB is included,
generating 8 tachyonic modes whose condensation reduces V by 99.9997 pct and produces
CKM mixing. The Seiberg seesaw is a metastable false vacuum.

**Known failures**: lepton masses not generated (Q = 0.48 from mesinos, Kahler exhausted),
mu_eff = -7 MeV (3 orders below needed), FCNC from off-diagonal scalars at m_S = f_pi
is excluded by K-Kbar by factor 8 (g=1) to 59 (NDA).

---

## 10 Specific, Actionable Items

---

### Item 1: Resolve the vacuum question — which vacuum is physical?

**The problem**: The global minimum (offdiag_vacuum.md) is M = 0, V = m_u^2. The Seiberg
seesaw vacuum has V ~ 10^15 MeV^2. The mass predictions (m_c, m_b, m_h) come from the
seesaw vacuum. The CKM mixing comes from tachyonic condensation starting at the seesaw
vacuum. The global minimum has no meson condensate and no CKM.

**What needs to be done**: Determine whether the seesaw vacuum is the cosmologically
selected one, and compute the tunneling lifetime from the seesaw to the M = 0 minimum.
If the lifetime exceeds the age of the universe, the model is consistent. If not, the
vacuum structure is wrong.

The key input: the bubble nucleation rate Gamma ~ exp(-S_bounce) where S_bounce is the
bounce action for the tunneling path in field space from M_i = C/m_i to M_i = 0. A rough
estimate uses the thin-wall approximation: S_bounce = 27 pi^2 sigma^4 / (2 epsilon^3) where
sigma ~ Delta V^(1/4) * Delta phi (surface tension) and epsilon ~ Delta V (energy gap).
With Delta V ~ 10^15 MeV^4 and Delta phi ~ C/m_u ~ 4 x 10^5 MeV, a crude estimate of
S_bounce can be obtained analytically in an afternoon.

**Without this calculation, the paper's vacuum structure is unjustified.**

- Difficulty: medium (analytic estimate possible; full numerical bounce hard)
- Impact: HIGH (touches the logical foundation of every prediction)
- Type: computational (bounce action estimate is analytic; full calculation is numerical)

---

### Item 2: FCNC resolution — identify which physical states are the off-diagonal mesons

**The problem**: fcnc_constraints.md shows K-Kbar mixing is violated by factor 8-59 if
off-diagonal scalars M^d_s sit at m_S = f_pi = 92 MeV. The paper currently has no answer.

**The resolution that is already implicit in the model but not stated explicitly**: The
off-diagonal mesons M^i_j (i != j) ARE the SM kaons, pions, and B mesons. They are not
new particles. K^0 = M^d_s bar, pi^+ = M^u_d, B^0 = M^d_b (extended to 5 flavors), etc.
If this identification is correct, there is no new FCNC source: the standard K-Kbar mixing
is already accounted for by SM QCD. The SUSY contribution to FCNC comes only from
loops involving the NMSSM coupling lambda X H_u H_d and the diagonal Yukawa terms y_c, y_b.

**What needs to be done**: Write down explicitly which SM meson is identified with which
element of M^i_j. Compute the NMSSM-induced FCNC (the loop diagrams involving X and H)
and show they are suppressed by lambda^2 / (16 pi^2) relative to the SM, safely below
bounds. This is an explicit 2-3 page computation that removes the most alarming-looking
objection.

**Key formula to check**: The off-diagonal meson SUSY coupling generates a four-fermion
operator (ds)(ds) with coefficient C_1 ~ lambda^4 X^2 / (m_mesino^2 m_H^2) from a
box diagram with X, H propagators. With lambda = 0.72, X_vac ~ 10^-9 MeV^-4,
m_mesino ~ 10^-4 MeV, this coefficient needs to be computed and compared to the UTfit
bound 9 x 10^-13 GeV^-2.

- Difficulty: medium (requires identifying the meson-SM particle dictionary + 1-loop calc)
- Impact: HIGH (currently looks like model is excluded; this removes the objection)
- Type: computational + conceptual (the identification is the key step)

---

### Item 3: mu-term — accept as open problem or identify a mechanism

**The problem**: mu_eff = lambda * vev(S) = -7 MeV (mu_term_analysis.md). This is 10^4
times too small for viable EWSB. The Higgs mass is correctly predicted by the NMSSM lambda
term independently (m_h = lambda v / sqrt(2) = 125.4 GeV), but EWSB at v = 246 GeV also
requires the bilinear B_mu term and mu ~ 100 GeV.

**Options** (each with different paper implications):

Option A: Giudice-Masiero. Add a Kahler coupling K_GM = kappa (X H_u H_d + h.c.) / M_Pl.
This generates mu_GM = kappa F_X / M_Pl. With F_X = -lambda v^2 / 2 ~ -(0.72)(246 GeV)^2 / 2
~ -2.2 x 10^10 MeV, mu_GM = kappa * 2.2 x 10^10 / 1.2 x 10^22 = 1.8 kappa x 10^-12 MeV.
Still tiny. Does not work.

Option B: Mu from a different NMSSM singlet. Introduce a second singlet S' uncharged under
the Seiberg constraint with W' = lambda' S' H_u H_d + kappa' S'^3. Standard NMSSM with
no connection to X. This restores the mu-problem to standard NMSSM territory (solvable but
adds a free parameter and loses the X=S identification).

Option C: Accept the mu-problem as open and note that mu_eff = -7 MeV could be relevant
if EWSB is not the conventional mechanism. Specifically, if the electroweak symmetry is
broken by a condensate at the QCD scale (as in technicolor-like scenarios), v could be
much smaller than 246 GeV in this sector. But the Higgs mass prediction uses v = 246 GeV
explicitly.

**Recommended action for the paper**: Be explicit that mu_eff = -7 MeV is physical (not
a dimensional mismatch) and that the model requires a separate mu-source. Cite this as
a known open problem, not a solved one. Do NOT claim the NMSSM coupling solves the
mu-problem in its standard form.

- Difficulty: hard (no clean solution known; best option is honest acknowledgment)
- Impact: HIGH (referees will ask; currently the paper appears to claim NMSSM solves mu)
- Type: requires new physics ideas (no purely computational solution)

---

### Item 4: CKM derivation from tachyonic condensation — complete the analytic treatment

**The problem**: assembly_round10.md shows F_X = -lambda v^2 / 2 = -2.18 x 10^10 MeV
dominates, making all 8 off-diagonal sectors tachyonic. The numerical minimization gives
theta_12 = 2.6 degrees (vs PDG 13.04 degrees), theta_23 = 14.1 degrees, theta_13 = 6.1
degrees. The Oakes relation tan(theta_C) = sqrt(m_d / m_s) = 0.2236 is separately accurate
to 0.3 pct. But the two results are inconsistent: the numerical minimization gives the
wrong Cabibbo angle, while the analytic Oakes relation gives the right one.

**What needs to be done**: Identify analytically which vacuum the Oakes relation
corresponds to, and why the numerical minimization finds a different result. The most
likely explanation: the numerical minimization converges to a saddle point or to the
M = 0 global minimum, not to the physical vacuum. The Oakes relation arises from
the mass-matrix diagonalization at the seesaw vacuum, before tachyonic condensation.
The perturbative off-diagonal VEVs at second order in (f_pi / M_i) give
delta M^i_j ~ (f_pi^2 / m_q) which, when rotated to diagonalize the mass matrix,
yield the Oakes angle.

**The calculation**: At the Seiberg seesaw vacuum with baryon sector stabilized (B = B-tilde = 0
enforced, removing the escape route to M = 0), expand M^i_j = (C/m_i) delta^i_j +
epsilon^i_j and minimize V_soft + V_F perturbatively in epsilon. Show that the off-diagonal
VEVs epsilon^d_s, epsilon^u_s, epsilon^u_d satisfy epsilon^d_s / (C/m_s) ~ sqrt(m_d/m_s),
leading to the Cabibbo angle. This is a 3x3 linear algebra problem and is completely doable.

- Difficulty: medium (second-order perturbation theory in 3x3 matrix)
- Impact: HIGH (either confirms the Oakes prediction as structural, or shows it is coincidental)
- Type: computational (analytic perturbation theory)

---

### Item 5: Lepton sector — commit to Sp(2) or demote leptons to input

**The problem**: lepton_alternatives.md identifies Sp(2) SQCD with N_f = 3 as the most
promising lepton sector. Round 11 Agent D confirms this is feasible. But the Sp(2) model
has not been verified to produce the O'Raifeartaigh vacuum, and the mediation scale
M_* ~ 43.5 GeV (required to get the right lepton mass scale) is uncomfortably low and
not derived from first principles.

**Option A (ambitious)**: Specify the Sp(2) model completely. Write down the superpotential,
identify the ADS vacuum, check gv/m = sqrt(3) condition, compute the mediation coupling.
Show m_e, m_mu, m_tau are predicted to within the known Koide precision. This is a
significant additional model and would substantially extend the paper.

**Option B (honest)**: Acknowledge that the lepton sector is not explained by the
meson-baryon sector alone. State that a separate confining sector (Sp(2) or otherwise)
is needed, and that the Koide relation for leptons is a PARALLEL structure to the quark
sector, not a derived consequence. Treat m_tau (or m_e) as an input and predict the others
from the Koide condition. This costs one free parameter but does not require a new model.

**Recommended path**: Option B for the current paper, with a note that Option A is in
development. The Sp(2) sector needs its own paper; forcing it into the current one would
make the paper too long and dilute focus.

- Difficulty: hard (Option A requires a new model); easy (Option B is a paragraph)
- Impact: medium (lepton masses are not unique to this model; the quark+Higgs sector is the advance)
- Type: requires new physics ideas (Option A); editorial decision (Option B)

---

### Item 6: Verify the soft-breaking identification V_soft = f_pi^2 Tr(M†M)

**The problem**: The soft SUSY-breaking term V_soft = f_pi^2 Tr(M†M) is identified with the
QCD soft scale, but no mediation mechanism is specified. The scale f_pi = 92 MeV is far
below the electroweak scale; this is unusual (soft masses are typically at the messenger scale
or higher). Furthermore, offdiag_vacuum.md shows this term DRIVES the vacuum to M = 0,
undermining the seesaw structure.

**What needs to be done**: Check whether V_soft = f_pi^2 Tr(M†M) can be derived from a
gauge-mediation or anomaly-mediation scenario in the confined theory. Specifically:
in anomaly mediation, soft masses for composites arise from superconformal anomalies:
m_soft^2 ~ (alpha_s / pi)^2 m_{3/2}^2. For m_soft ~ f_pi, this requires the gravitino
mass m_{3/2} ~ f_pi * pi / alpha_s ~ 3 GeV. This is a specific prediction for the
gravitino mass that can be compared to cosmological constraints.

Alternatively, if the soft masses arise from direct gauge mediation via a messenger sector,
the messenger scale F_mess / M_mess ~ f_pi gives a relation between the messenger mass and
F-term that can be constrained.

**Key question for the paper**: Is V_soft = f_pi^2 Tr(M†M) a prediction of a specific
mediation mechanism, or an assumption? The paper currently treats it as an assumption.
If a mediation mechanism can be identified that fixes f_pi^2, this becomes a prediction.

- Difficulty: medium (literature check + parametric estimates)
- Impact: medium (changes parameter count and narrative coherence)
- Type: requires new physics ideas (mediation mechanism)

---

### Item 7: Clarify the baryon sector — what enforces B = B-tilde = 0 at the physical vacuum?

**The problem**: The seesaw vacuum has B = B-tilde = 0 and M_i = C/m_i. The global minimum
has BB-tilde = -Lambda^6 and M = 0 (offdiag_vacuum.md). The mass predictions require the
seesaw vacuum. But nothing in the Lagrangian forces B = B-tilde = 0; the baryonic flat
direction is a genuine modulus.

**Physical interpretation**: In real QCD with N_f = 3 and N_c = 3, there is no baryon
condensate in the vacuum (baryon number is conserved, and the ground state has B = 0 by
Lorentz invariance). The SUSY case is different: BB-tilde is a Lorentz scalar, so
a nonzero BB-tilde is not forbidden by Lorentz invariance. However, it does break baryon
number U(1)_B, which is a global symmetry of the Seiberg superpotential.

**Key question**: Is U(1)_B spontaneously broken in the sBootstrap model? If so, there
should be a Goldstone boson (a baryon-number axion). If not, what stabilizes B = B-tilde = 0?

A possible stabilization: add a baryon mass term W_B = m_B B B-tilde to the superpotential.
This breaks U(1)_B explicitly and gives a potential m_B^2 |B|^2 + m_B^2 |B-tilde|^2,
lifting the baryonic flat direction. The constraint det M = Lambda^6 is then enforced by
F_X = 0 directly, restoring the seesaw as the only minimum.

**Recommended action**: Add W_B = m_B B B-tilde (or equivalently, integrate out B, B-tilde
as massive at some scale m_B > Lambda) and recheck whether the seesaw vacuum is stabilized.
This removes the global minimum problem and the BB-tilde escape route.

- Difficulty: easy (add one term; check vacuum stability analytically)
- Impact: HIGH (fixes the global minimum problem; makes seesaw the true vacuum)
- Type: computational (add term, redo vacuum analysis)

---

### Item 8: Notation and convention consistency in the paper

**The problem**: consistency_review.md (reviewed sbootstrap_v4d.tex) found:
(a) Q is overloaded: Koide ratio, electric charge, and SUSY generator all called Q;
(b) The Q = 2/3 vs Q = 3/2 convention switch footnote contains an arithmetic error
    ("Q_here = 3 Q_prev" is wrong; it should be Q_here = 1/Q_prev = 3/2 when Q_prev = 2/3);
(c) v_0^2 = 29,580 MeV (table) vs 29,600 MeV (text) -- 0.07 pct rounding discrepancy.

**What needs to be done**: Before submission:
- Replace Q_em for electric charge throughout (or use e_Q or Q_{em}).
- Replace script-Q for the SUSY generator.
- Fix the footnote: Q_here = 1 / Q_prev, not 3 Q_prev.
- Synchronize the v_0^2 rounding to one value (29,580 MeV, the more precise one).

This is mechanical editing but will draw referee attention if not fixed.

- Difficulty: easy (editorial; 1-2 hours)
- Impact: medium (professional quality; wrong formulas attract attention disproportionate to significance)
- Type: computational (find-and-replace + arithmetic check)

---

### Item 9: Statistical significance — quark Koide triples in the meson pool

**The problem**: open_problems_brainstorm.md notes that when meson masses are included in
the pool, (c,b,t) ranks 113/4060 and (-s,c,b) ranks 324/4060. Only (e,mu,tau) survives
at rank 1. This undermines the claim that the three Koide triples are simultaneously
significant.

**Two responses the paper can take**:

Response A (structural): The significance of (c,b,t) and (-s,c,b) is NOT from their
raw Q proximity to 2/3 but from their structural role in the sBootstrap: they are
PREDICTED to have Q = 2/3 by the model (the O'Raifeartaigh seed + bloom), and the
fact that they approximately satisfy it (within 1.2 pct) is a consistency check, not
a discovery. The model predicts these three triples specifically; the look-elsewhere
problem does not apply to predicted targets.

Response B (quantitative): Compute the probability that a random assignment of 6 quark
masses (drawn from the same distribution as PDG values) produces two overlapping triples
both within 1.5 pct of Q = 2/3 while sharing two common members (b, c). This combinatorial
question has a sharper answer than the brute-force rank comparison, because the structural
constraint (overlapping triples, shared members) dramatically reduces the look-elsewhere
factor.

**Recommended action**: Use Response A in the paper (one paragraph clarifying the
prediction vs discovery distinction) and include Response B as a quantitative footnote.
This reframes the statistics correctly without overclaiming.

- Difficulty: easy to medium (Response A is editorial; Response B requires a short calculation)
- Impact: medium (pre-empts a common referee objection)
- Type: computational (short Monte Carlo for Response B)

---

### Item 10: Write the Yukawa Lagrangian for quarks explicitly, showing quark-meson identification

**The problem**: The Yukawa couplings y_c H_u M^d_d and y_b H_d M^s_s are written in terms
of mesons M^d_d and M^s_s. But the identification of M^d_d with the charm-quark Yukawa
(rather than, say, M^u_u or M^s_s) is ASSUMED, not derived from the quark quantum numbers.
Specifically: why does the down-type-strange meson M^d_d couple to the charm quark through
H_u, and the strange-strange meson M^s_s couple to bottom through H_d?

**The derivation**: In the ISS seesaw, M^i_i = Lambda^2 / m_i at the vacuum. The meson
M^d_d has VEV proportional to 1/m_d (large), and M^s_s has VEV proportional to 1/m_s
(smaller). The NMSSM coupling X H_u H_d gives Higgs an effective mass from F_X. The Yukawa
y_c M^d_d H_u generates a charm mass y_c * (Lambda^2 / m_d) = y_c * C / (m_d m_s m_u)^{1/3} / m_d
at the vacuum. For this to equal m_c, y_c = m_c m_d / (Lambda^2).

This derivation shows that the flavor assignment (which meson couples to which Higgs) is
determined by a mass hierarchy argument: the meson with the larger VEV (M^d_d > M^s_s,
since m_d < m_s) couples to the heavier quark (charm > bottom). This is a structural
argument, not a coincidence, but it has not been written down explicitly in the paper.

**What needs to be done**: Add a half-page derivation showing the flavor assignment follows
from the seesaw VEV hierarchy and the requirement that mass eigenvalues match the observed
spectrum. This closes a gap that every careful reader will notice.

- Difficulty: easy (the derivation is already implicit in the numbers; needs to be made explicit)
- Impact: medium-high (removes an apparent assumption, showing it is a consequence of the structure)
- Type: computational (algebraic; 1-2 pages)

---

## The Single Most Important Calculation Not Yet Done

**Item 7 (baryon stabilization) is the most important, but it is also the quickest**.
The single most important CALCULATION (not editorial fix) that has not been done is:

### Perturbative CKM from the Seiberg seesaw with B = B-tilde = 0 enforced

(A combination of Items 4 and 7.)

**The setup**: Stabilize the baryonic flat direction by adding W_B = m_B B B-tilde with
m_B >> f_pi. Integrate out B, B-tilde at the scale m_B. The effective theory below m_B
has only mesons M^i_j, X, H_u, H_d. At this level, the constraint det M = Lambda^6 is
enforced by F_X = 0 exactly.

Now include V_soft = f_pi^2 Tr(M†M) as a perturbation. At second order in f_pi^2 / M_i^2,
the off-diagonal meson VEVs epsilon^i_j appear. Their values satisfy a linear system:

    (m_ij^2) epsilon^i_j = - (∂^2 V / ∂M^i_j ∂M^k_l) epsilon^k_l - f_pi^2 (∂^2 K / ∂M^i_j ∂M^k_l) epsilon^k_l

At leading order in the perturbation, epsilon^d_s / (C/m_s) ~ f_pi^2 / (X_vac M_u M_d).
With X_vac ~ -1.2 x 10^-9 MeV^-4 and M_u, M_d ~ 10^5 MeV:

    epsilon^d_s ~ f_pi^2 / (|X_vac| M_u M_d) ~ 8464 / (1.2e-9 * 4e5 * 2e5) ~ tiny

But with F_X = -lambda v^2 / 2 (the NMSSM contribution from EWSB), the effective coupling
between off-diagonal mesons and the vacuum is much larger. The tachyonic mass-squared for
the ds sector is m^2_ds = f_pi^2 - lambda^2 v^2 (M_u / Lambda^2). The condensate epsilon^d_s
from the tachyonic instability satisfies a gap equation whose solution gives the Cabibbo angle.

**Why this calculation is decisive**: If the condensate epsilon^d_s / (C/m_s) = sqrt(m_d/m_s)
follows analytically from this gap equation, then the Oakes relation is a DERIVED result of
the model, not a coincidence. This would be the strongest analytic result in the paper:
a first-principles derivation of the Cabibbo angle from the quark mass hierarchy.

If the gap equation gives a different answer, the CKM connection is broken and must be demoted
from "prediction" to "observation that the model is consistent with."

**Estimated effort**: 3-5 days for a careful analytic treatment. The gap equation is a 3x3
(or 9x9) system of quadratic equations in the off-diagonal VEVs, linearizable at small epsilon.
The key algebraic step is the connection between the tachyonic mass matrix eigenvalues and
the quark mass ratios.

**This is the calculation that determines whether CKM is truly explained by the model.**

---

## Summary Table

| Item | Topic | Difficulty | Impact | Type |
|------|-------|------------|--------|------|
| 1 | Vacuum tunneling lifetime | medium | HIGH | computational |
| 2 | FCNC resolution (meson = SM particle) | medium | HIGH | computational + conceptual |
| 3 | mu-term: accept as open problem | hard | HIGH | new physics ideas |
| 4 | CKM analytic from tachyonic condensation | medium | HIGH | computational |
| 5 | Lepton sector: Sp(2) or demote | hard / easy | medium | new physics ideas |
| 6 | V_soft mediation mechanism | medium | medium | new physics ideas |
| 7 | Baryon stabilization (W_B term) | easy | HIGH | computational |
| 8 | Notation and convention fixes | easy | medium | editorial |
| 9 | Statistical significance framing | easy-medium | medium | computational |
| 10 | Yukawa flavor assignment derivation | easy | medium-high | computational |

**Priority ordering for a final paper**: 7 → 2 → 4 → 1 → 10 → 8 → 9 → 3 → 5 → 6

Items 7, 2, 4 together constitute the "CKM problem": stabilize the vacuum, identify the
mesons with SM particles, derive the Cabibbo angle. These three together would be the
primary new result of the paper beyond what was in arXiv:2407.05397.

Items 3 and 5 (mu-problem, leptons) require new ideas and should be stated as open problems
with a clear description of what is missing and why, rather than hand-waved away.

---

## What a Referee Will Ask

A careful referee in EPJC or JHEP will ask, in roughly this order:

1. "The global minimum is M = 0, not the seesaw. Why is the seesaw vacuum physical?"
   (Item 1 + Item 7 -- the most fundamental objection.)

2. "Off-diagonal scalars at 92 MeV violate K-Kbar mixing by a factor of 8-59.
    How does the model survive this constraint?"
   (Item 2 -- currently has no answer in the paper.)

3. "The mu-term from the NMSSM coupling gives mu_eff = -7 MeV, three orders too small.
    Where does the physical mu come from?"
   (Item 3 -- currently acknowledged but unsolved.)

4. "The CKM mixing angles from your minimization (theta_12 = 2.6 degrees) do not match
    PDG. How do you derive the Cabibbo angle from first principles?"
   (Item 4 -- the Oakes relation is known but not derived.)

5. "You have no mechanism for lepton masses. Does the theory explain them or not?"
   (Item 5 -- needs a clear answer: either the Sp(2) sector or honest demotion to input.)

6. "The notation uses Q for three different quantities. The Q=2/3 vs Q=3/2 footnote
    has an arithmetic error. Please fix."
   (Item 8 -- trivial but will be noted.)

7. "The statistical significance of the quark Koide triples is only 113/4060 and 324/4060
    in a fair sampling. How are these 'predicted' rather than selected after the fact?"
   (Item 9 -- needs the prediction-vs-discovery reframing.)

---

*Generated: 2026-03-04*
*Based on: complete_lagrangian.md, offdiag_vacuum.md, fcnc_constraints.md,*
*mu_term_analysis.md, lepton_alternatives.md, assembly_round10.md, assembly_round11.md,*
*open_problems_brainstorm.md, consistency_review.md, ckm_koide.md, bloom_mechanism.md,*
*v0_doubling.md*
