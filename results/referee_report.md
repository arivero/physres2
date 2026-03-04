# Referee Report

**Manuscript:** "A Tale of Two Tuples: Seed Triples, Signed Charges, and SUSY Breaking in the sBootstrap"
**Author:** A. Rivero
**Submitted to:** Physical Review D
**Referee:** Anonymous
**Date:** March 2026

---

## Summary of the Paper

The author extends the "sBootstrap" program, which proposes that Standard Model pseudoscalar mesons and charged leptons fill common N=1 supermultiplets with an SU(5) flavor symmetry. The paper introduces "signed mass charges" z = z_0 + z_i with m = z^2, reinterprets the Koide formula Q = 2/3 as an energy-balance condition, identifies "seed triples" with one massless member as pre-breaking configurations, discusses a meson Koide triple (-pi, D_s, B), revisits the Diophantine generation-counting argument, addresses the Pauli exclusion conflict for the symmetric 15 diquark representation, embeds the representations in SO(32), derives an electroweak mixing angle from a Casimir eigenvalue equation, and sketches how SUSY breaking proceeds through O'Raifeartaigh and Volkov-Akulov mechanisms. The paper claims to reduce the SM quark mass spectrum to two or three free parameters.

---

## Major Criticisms

**1. The central claim -- that mesons and leptons are superpartners -- has no dynamical derivation and faces a fundamental obstruction.**

The paper asserts that the pion-muon mass relation m_pi^2 - m_mu^2 = f_pi^2 (Eq. 1, labeled eq:pimu) is a "Volkov-Akulov mass formula for a pseudo-Goldstino." This is presented as the foundational equation of the program. However:

(a) The Volkov-Akulov construction describes the nonlinear realization of broken spacetime SUSY, where the goldstino is a neutral Majorana fermion. The muon is a charged Dirac fermion. The paper acknowledges this mismatch (Sec. 8.2, paragraph on quantum numbers) but dismisses it by asserting the sBootstrap SUSY is "emergent hadronic" rather than spacetime. This is circular: the Miyazawa superalgebra is invoked to justify the identification, but no Miyazawa superalgebra Lagrangian is constructed. The claim that "the quantum numbers are inherited from the broken hadronic generator" is not demonstrated; it is stated.

(b) The numerical relation m_pi^2 - m_mu^2 = f_pi^2 holds to 2.2%. A skeptical referee notes that m_pi = 139.57 MeV, m_mu = 105.66 MeV, and f_pi = 92.2 MeV are all of order 100 MeV. The dimensionless coincidence being tested is of order unity, and a 2.2% match among quantities of comparable magnitude is not extraordinary. No statistical analysis of the significance of this particular relation is provided, unlike the careful Monte Carlo treatment given to the lepton Koide formula (Sec. 6.1).

**2. The signed-charge formalism is a rewriting, not a derivation.**

The paper presents z = z_0 + z_i with m = z^2 as if it were a physical framework. But this is just the standard parametrization of Koide triples dressed up with the word "charge." The derivation in Sec. 2.2 is correct but trivially follows from algebra: given any three positive quantities whose Koide ratio Q = 2/3, one can define z_k = sqrt(m_k) and decompose into mean plus deviation. Calling z a "signed mass charge" does not promote this decomposition to physics. The energy-balance interpretation (perturbation energies equal vacuum energy) is suggestive language, not a physical mechanism. The paper never shows what field z corresponds to in any Lagrangian, what symmetry protects the balance, or what dynamics enforces it.

The extension to negative z (Sec. 5) is the genuinely new ingredient, allowing the meson triple (-pi, D_s, B). But the freedom to choose signs introduces ambiguity: for each triple of three masses, there are 2^3 = 8 sign assignments. The paper does not systematically address how many of the 84 x 8 = 672 possible signed triples satisfy Q ~ 2/3, which would be necessary to assess the look-elsewhere effect for the meson triple claim.

**3. The "predictions" are predominantly post-dictions, and the claimed parameter reduction is not demonstrated.**

The paper claims that after applying Koide constraints and vacuum-charge doubling, "only two to three free parameters remain: m_u, m_d, and m_s." This would be a spectacular reduction from the SM's 6 quark masses to 3 inputs. But the logic is not self-consistent:

(a) The seed Koide on (0, s, c) with exact ratio gives m_c = 1301 MeV. The PDG value is 1270 MeV, a 2.4% deviation. The v_0-doubling with the *physical* m_c (not the predicted one) gives m_b = 4177 MeV (0.07% from PDG). The paper openly acknowledges (Sec. 4.1, "A tension exists...") that these two constraints "slightly overconstrain the system." This is a serious problem: if the seed Koide is exact, the v_0-doubling fails at 1.8 sigma. If the v_0-doubling is exact, the seed ratio deviates. These are not small refinements -- they are incompatible constraints at the stated level of precision. The paper does not resolve this; it presents both as successes while acknowledging they contradict each other.

(b) The charm mass prediction is consistently off by 7-8% across multiple routes (Chain A: 1357 MeV, Chain B: 1360 MeV, vs PDG 1270 MeV). A 7% deviation is not a "prediction" in any meaningful sense -- it is a rough order-of-magnitude match. For comparison, the Gell-Mann-Okubo mass formula achieves sub-percent accuracy for baryon octets. Claiming a parameter reduction while being 7% off on one of the "predicted" parameters undermines the central narrative.

(c) The v_0-doubling "prediction" sqrt(m_b) = 3*sqrt(m_s) + sqrt(m_c) is presented as if it follows from a theoretical principle ("the bloom exactly doubles the vacuum charge"). But no mechanism for exact doubling is derived. The ratio v_0(full)/v_0(seed) = 2.0005 is observed and then elevated to a principle. This is pattern-fitting, not prediction. The paper should state clearly: "We observe that v_0 approximately doubles and explore the consequences," rather than "the bloom doubles v_0, predicting..."

**4. The look-elsewhere problem is inadequately addressed for the key claims.**

(a) The meson Koide triple (-pi, D_s, B) at Q = 0.6674 (0.10% from 2/3): No look-elsewhere correction is performed. The pseudoscalar meson spectrum contains pi, K, eta, eta', D, D_s, B, B_s, B_c, eta_c, eta_b -- at least 11 states. With signed charges, the number of possible triples is C(11,3) x 8 = 1320. Finding one at 0.1% is not evaluated against this background. This omission is particularly striking given that the paper performs an explicit look-elsewhere analysis for the lepton Koide (Sec. 6.1). The asymmetric treatment suggests the author is aware that the meson triple might not survive the correction.

(b) The 2/9 phase coincidence (Sec. 6.2): The paper quotes 3.1 sigma after look-elsewhere correction over 128 "nice" fractions with q <= 20. But why q <= 20? The cutoff is arbitrary. With q <= 100, the number of fractions increases dramatically, and the significance drops. The choice of q <= 20 is not justified by any physical argument.

(c) The v_0-doubling: The ratio 2.0005 is presented without any error analysis. What is the uncertainty on this ratio, given uncertainties in m_s, m_c, m_b? Given that m_s has ~5% uncertainty at the PDG level, the uncertainty on v_0(seed) is not negligible.

**5. The RG issue is fatal unless properly addressed.**

The paper states (Sec. 8.7) that one-loop SUSY RG drives Q -> 1/3, not 2/3, and that therefore the Koide condition must be a "boundary condition on the superpotential." This is presented in one paragraph and then set aside. But this is a catastrophic problem for the program. If the dynamics of the theory actively drives masses *away* from the Koide surface, then:

(a) Any Koide-satisfying boundary condition at a UV scale will be destroyed by running to the IR.
(b) The claim that QCD mass ratios are RG invariants (Sec. 3.2 footnote) is true only at leading order, and the higher-order corrections are of the same order as the claimed Koide precision (~1%).
(c) The paper's resolution -- that Koide is a boundary condition at the SUSY-breaking scale, and mass ratios don't run -- requires specifying the scale. The quark masses used (m_s at 2 GeV, m_c at m_c, m_b at m_b) are evaluated at *different* scales. The paper acknowledges this in a footnote (Sec. 3.2) and says the meson seed "is free of this scheme issue." But the quark seed and the quark Koide triples are not free of this issue, and these are the central quantitative claims.

**6. The Pauli theorem is correctly stated but its consequences are not honestly assessed.**

The theorem in Sec. 7.3 correctly proves that no J=0 color-3bar diquark can have symmetric flavor. The paper then says this "is not a weakness of the sBootstrap but a prediction: the diquark sector requires new physics beyond the qq picture." This is a remarkable rhetorical move. The sBootstrap claims that diquarks are superpartners of antiquarks. The Pauli theorem proves these diquarks cannot exist as qq composites. The paper's response is to declare this a feature rather than a bug, by invoking either the Miyazawa supercharge (which has not been constructed for this system) or SO(32) string theory (which has not been compactified to produce the required spectrum). Both escape routes are promissory notes, not resolutions.

A more honest assessment: the Pauli theorem proves that the sBootstrap diquark sector, as currently formulated, is inconsistent with QCD. The proposed resolutions (string theory, algebraic supercharges) would need to be developed into concrete models before this inconsistency can be considered resolved.

**7. The SO(32) embedding and the Casimir equation are disconnected from the rest of the paper.**

Section 7 ("Two Algebraic Structures") presents two independent algebraic observations. The paper explicitly states (Sec. 7.5) that "No algebraic homomorphism linking the embedding chain to the eigenvalue equation has been found." These are therefore two separate observations appended to a paper about Koide triples and SUSY breaking. They share a common set of particles but are not integrated into a unified argument. A referee must ask: do they strengthen the paper's case, or do they dilute it?

The Casimir eigenvalue equation x^2 + C_2 x - C_2 = 0 is presented with three "structural constraints" that fix the polynomial. But these constraints (f_0 = 0, g_0 = 0, g_1 = -f_1) are not derived from any physical principle. They are reverse-engineered to produce the desired result. The paper calls them "structural" rather than "phenomenological," which obscures their ad hoc character. Furthermore, the equation uses the on-shell weak mixing angle, not the MSbar value. The paper correctly notes this (Sec. 7.4, "Renormalization scheme dependence") but does not resolve it. The on-shell and MSbar values differ by 0.008, much larger than experimental errors. Choosing the scheme that matches is not a prediction -- it is a fit.

**8. The superpotential sketch (Sec. 8.8) is not a superpotential.**

The paper's stated north star is to construct a SUSY Lagrangian for the full theory. Section 8.8 writes down the standard Seiberg SQCD effective superpotential with a mass matrix, which is well-known. The new content is the claim that F-term breaking at one massless flavor reproduces the seed triple, and that the O'Raifeartaigh condition gv/m = sqrt(3) fixes the Koide ratio. However:

(a) The bloom from seed to full triple is stated to be "non-holomorphic" and therefore cannot arise from the superpotential (Sec. 8.8, "Fixed by the vacuum-charge doubling"). This means the actual dynamics responsible for the key mass predictions is not captured by the superpotential. The Lagrangian is therefore incomplete precisely where it matters most.

(b) The Seiberg SQCD is for N_f = N_c = 3, but the paper works with N_f = 5 flavors. For N_f = 5, N_c = 3, Seiberg duality gives a free magnetic phase, not confinement. The paper mentions this (Sec. 8.5) but does not resolve the tension -- it gestures at "scale-dependent N_f" without a concrete implementation.

(c) The Coleman-Weinberg mass for the pseudo-modulus is calculated at m_CW ~ 13 MeV, while the physical m_b = 4180 MeV. The paper calls this "five orders of magnitude" (it is actually a factor ~320, not 10^5, so this appears to be an error). Regardless, the bloom requires a mass enhancement far beyond any known perturbative effect, and the claim that it must be "nonperturbative" is a placeholder, not an explanation.

---

## Minor Criticisms

**1.** The two conventions for the Koide ratio (Q = 2/3 in Sec. 2 and Q = 3/2 in Sec. 6) are confusing. The footnote in Sec. 6 attempts to clarify but actually makes things worse by introducing Q_here and Q_prev. The author should pick one convention and use it throughout.

**2.** The vacuum energy table in Sec. 3.5 lists v_0^2 values with units of MeV. But v_0 has dimensions MeV^{1/2}, so v_0^2 has dimensions MeV, and calling it "vacuum energy" is misleading. The actual energy density or energy scale associated with this quantity is not defined.

**3.** The paper states the lepton Koide formula has significance "2.8 sigma" after look-elsewhere correction. This is below the standard 3-sigma evidence threshold. The paper should not call this "suggestive but not conclusive" (which implies the result is interesting) but rather "not significant" (which is the standard statistical assessment at 2.8 sigma for a look-elsewhere-corrected result).

**4.** The "Nature Abhors Chirality" section (Sec. 10) is speculative editorial commentary that does not advance the technical content. It belongs in a discussion or conclusion, not as a standalone section.

**5.** The reference list is sparse. Many claims lack proper citations:
- The "dual Koide" observation Q(1/m_d, 1/m_s, 1/m_b) = 0.665 (Sec. 8.5) -- is this new or previously known?
- The CKM-Koide connection (Sec. 8.7, "CKM mixing and the Koide triples") -- has this been observed before?
- The v_0-doubling result -- is this first reported here?
- H. de Vries is credited for the Casimir eigenvalue observation (Sec. 7.4, "Origin") but no publication is cited. A Physics Forums post is not a citable reference.
Without clear attribution, it is difficult to assess what is new in this paper versus what was in the prior EPJC publication.

**6.** In Sec. 3.2, the footnote on quark mass conventions is important but buried. The fact that m_s is evaluated at 2 GeV while m_c is at m_c means the ratio sqrt(m_s)/sqrt(m_c) = 0.271 is not physically meaningful as stated. The footnote says "at a common scale the ratio is ~0.29, further from 2 - sqrt(3) = 0.268." This worsens the quark seed match from 1.2% to ~8%, which significantly undermines the claim. This should be in the main text, not a footnote.

**7.** Table 2 (Sec. 4.2, "Precision summary"): The (-s, c, b) triple has Q = 0.675, deviating from 2/3 by 1.24%. This is quite poor compared to the other triples. Yet this is the "full triple" that the entire seed-to-bloom narrative hinges upon. If the bloom is supposed to preserve Q = 2/3, a 1.24% deviation needs explanation.

**8.** The "cyclic echo" result (Sec. 6.4, Step 4) is presented as a "striking algebraic consistency check," but the paper also states it is "an algebraic necessity." If it is necessary, it is not striking -- it is a tautology. The presentation conflates the two.

**9.** Section 8.1 states that for N_f = N_c = 3 the quantum-modified constraint is det M - B Bbar = Lambda^{2N_c}. The exponent should be 2N_c = 6, which is correctly used later, but the notation Lambda^{2N_c} in Eq. (eq:seiberg) is potentially confusing when N_c was not yet defined in that equation's context.

**10.** The factor-of-20 discrepancy between <F> ~ |qqbar|/f_pi ~ 1.7 x 10^5 MeV^2 and f_pi^2 ~ 8500 MeV^2 (Sec. 8.2) is significant and inadequately addressed. The paper says the sBootstrap identification f_VA = f_pi is a "low-energy effective statement," but does not explain why the microscopic F-term is 20 times larger than the effective one.

**11.** The claim that the Koide condition is "immune to D-terms" (Sec. 8.8) because D-terms shift scalar masses but leave fermion masses untouched applies only if the Koide triples are conditions on fermion masses. But the meson Koide triple (-pi, D_s, B) involves meson masses, which are scalar masses and *are* affected by D-terms. The paper's own framework requires Koide to hold for both fermion and scalar triples, so this immunity argument is inconsistent.

**12.** Equation numbering is inconsistent. Some key equations are numbered, others (including several predictions and derived relations) are unnumbered equation* environments. This makes cross-referencing difficult.

---

## Questions for the Author

1. Can you perform and report a look-elsewhere-corrected significance for the meson Koide triple (-pi, D_s, B), analogous to the Monte Carlo analysis performed for the lepton Koide formula? Specifically: scan all pseudoscalar meson triples with all sign assignments, and compute the probability that the best one matches Q = 2/3 at or below the observed 0.10%.

2. You state that the v_0-doubling gives m_b = 4177 MeV with PDG masses as input, but gives m_b = 4233 MeV (1.8 sigma from PDG) when the exact seed value m_c = 1301 MeV is used as input. Which is the actual prediction of the theory? If the theory predicts the exact seed (m_c = 1301 MeV), then the predicted m_b is 4233 MeV, which is 1.8 sigma off. If the theory uses the physical m_c (1270 MeV), then what fixes m_c?

3. What is the renormalization scale at which the Koide condition is supposed to hold for quarks? If it is scheme-independent (as claimed via mass-ratio invariance), then why does running the quark masses to a common scale worsen the quark seed match from 1.2% to ~8%?

4. The O'Raifeartaigh model with gv/m = sqrt(3) gives a seed with eigenvalues (0, (2-sqrt(3))m, (2+sqrt(3))m). But the bloom from this seed to the full triple is claimed to be "nonperturbative" and to break Q = 2/3 under any perturbative R-breaking deformation. How, then, does the bloom proceed while preserving Q? If the bloom is truly nonperturbative, what theoretical control do you have over it?

5. The ISS mechanism requires N_c < N_f < 3N_c/2, i.e., 3 < N_f < 4.5 for N_c = 3. You state that with charm included, N_f = 4 falls in this window. But the ISS vacuum is metastable, and the lifetime depends on details of the potential. Have you estimated whether the ISS metastable vacuum has a cosmologically acceptable lifetime for the parameter values relevant to QCD?

6. The Casimir eigenvalue equation x^2 + C_2 x - C_2 = 0 has the "codimension analysis" claiming that the equation is pinned by three structural constraints. But these constraints are stated axiomatically. Can you derive any of them from a known principle? In particular, why should g_1 = -f_1 (the "coefficient antisymmetry")? This is the constraint that produces the specific numerical value of R.

7. You assert that "the sBootstrap diquark in the 15 cannot be a qq composite" and that this is a "prediction." What experimental signature would confirm or falsify this prediction? Exotic +4/3-charged states are mentioned -- at what mass scale? Are they within reach of the LHC?

8. The "dual Koide" result Q(1/m_d, 1/m_s, 1/m_b) = 0.665 is intriguing. Have you checked whether this is stable under variation of m_d (which has a ~20% PDG uncertainty)?

---

## Assessment of Novelty

Compared to the author's prior work (arXiv:2407.05397, EPJC 84, 1058 (2024)), the genuinely new elements appear to be:

- The signed-charge formalism and the meson Koide triple (-pi, D_s, B) -- **new and interesting**
- The seed-triple concept with the universal ratio 2 - sqrt(3) -- **new framing of known observations**
- The O'Raifeartaigh connection with gv/m = sqrt(3) -- **new and the strongest technical result in the paper**
- The v_0-doubling relation sqrt(m_b) = 3*sqrt(m_s) + sqrt(m_c) -- **new observation, but post-dictive**
- The Pauli theorem for all-L diquarks -- **new and clean result**
- The bloom-instability observation (R-breaking pushes Q away from 2/3) -- **new negative result, important**
- The dual Koide and CKM-Koide connection -- **new observations**
- The superpotential sketch in Sec. 8.8 -- **new but incomplete**

The Diophantine counting, SO(32) decomposition, Casimir equation, and basic Koide phenomenology appear to have been in the prior paper (or in the earlier viXra versions). The Volkov-Akulov connection was implicit before and is made explicit here.

---

## Overall Assessment

The paper contains a mix of genuinely interesting algebraic observations, creative physical interpretation, and insufficiently rigorous analysis. The strongest elements are: (1) the O'Raifeartaigh seed construction with its single-constraint derivation; (2) the Pauli theorem closing all field-theoretic diquark escape routes; and (3) the honest reporting of negative results (scaling fails at second step, chain stalls, RG drives Q -> 1/3, seed and doubling overconstrain). The weakest elements are: (1) the absence of any dynamical mechanism for the bloom; (2) the inconsistent treatment of statistical significance (careful for leptons, absent for mesons); (3) the persistent 7% charm mass discrepancy presented alongside claims of parameter reduction; and (4) the confusion between observation and prediction throughout.

The paper would benefit enormously from: (a) a clear separation of what is derived vs. observed vs. assumed; (b) a systematic look-elsewhere analysis for all claimed Koide triples, not just the lepton one; (c) an honest reckoning with the quark mass scale ambiguity; and (d) either dropping the superpotential section or developing it to the point where it actually produces the claimed spectrum.

---

## Recommendation

**Major revision required.** The paper is not suitable for publication in PRD in its present form, but it contains enough novel content to warrant revision rather than outright rejection. The following must be addressed:

1. Perform and report look-elsewhere-corrected significances for *all* claimed Koide triples, not just the lepton one. The meson triple in particular requires this analysis.

2. Resolve the overconstrained system: either the seed Koide determines m_c = 1301 MeV (and the v_0-doubling prediction for m_b is 4233 MeV, 1.8 sigma off), or the v_0-doubling uses the physical m_c (and the seed ratio is not exact). The paper cannot have both. State clearly which constraint is fundamental and which is approximate.

3. Address the quark mass scale issue in the main text, not a footnote. If the quark seed ratio worsens to ~8% at a common scale, this must be acknowledged as a significant limitation, and the claim that "QCD mass ratios are RG invariants" must be qualified.

4. Either develop the superpotential (Sec. 8.8) to the point where it actually reproduces the bloom -- specifying the Kahler potential or nonperturbative mechanism that generates the v_0-doubling -- or clearly label it as a sketch that does not yet produce the claimed spectrum.

5. Separate the two algebraic structures (SO(32), Casimir equation) into a clearly delineated section that does not imply connection to the Koide/SUSY breaking parts of the paper, since no connection has been found.

6. Provide a clear, tabulated distinction between: (a) mathematical identities that hold by construction; (b) empirical observations that happen to hold approximately; (c) genuine predictions that could be falsified. The current text blurs these categories throughout.

If these revisions are made honestly and thoroughly, the paper could make a contribution to PRD as a phenomenological observation paper. It should not, however, claim to have constructed a theory or a Lagrangian until the dynamical mechanisms (the bloom, the Kahler origin of v_0-doubling, the Miyazawa supercharge) are actually built.
