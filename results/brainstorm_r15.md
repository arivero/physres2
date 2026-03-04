# Round 15 Brainstorm: N=1 SQCD with Electroweak Coupling

**Date:** 2026-03-04

---

## (a) CKM mixing: dead or relocated?

The tree-level tachyonic CKM mechanism is dead at the pre-EWSB Seiberg vacuum. At M_i = C/m_i with B = B-tilde = 0, we have F_X = 0 identically. The off-diagonal meson mass matrix is diagonal with eigenvalues m^2 = f_pi^2 + 2 X_0^2 M_k^2 > 0 for all sectors (the X_0^2 M_k^2 correction is at most 5 x 10^{-7} MeV^2, negligible against f_pi^2 = 8464 MeV^2). There are no tachyons and no CKM mixing at tree level.

Three alternative mechanisms:

### Alternative 1: One-loop Coleman-Weinberg corrections

The CW potential generates off-diagonal meson couplings at one loop through diagrams where the virtual X propagator connects pairs of off-diagonal mesons. The relevant quantity is the second derivative of V_CW with respect to the off-diagonal entries epsilon^a_b, evaluated at epsilon = 0.

The CW mass contribution to the off-diagonal sector is

    delta m^2_CW ~ (1/16 pi^2) * |W_{X,M^a_b,M^b_a}|^2 * log(M_heavy^2/Lambda^2)
                 ~ (1/16 pi^2) * M_k^2 * log(M_heavy^2/Lambda^2)

where M_heavy ~ X_0 * M_u * M_d ~ 7.7 x 10^{10} MeV is the heaviest fermion mass in the central block. Taking k = u (the ds sector):

    delta m^2_CW ~ (1/16 pi^2) * (4.08 x 10^5)^2 * log((7.7e10)^2/(300)^2)
                 ~ (1/16 pi^2) * 1.67 x 10^{11} * 38.7
                 ~ 4.1 x 10^{10} MeV^2

This is formally large, but there is a subtlety: the CW potential at the seesaw vacuum involves a supertrace over the full spectrum. The supertrace STr[M^2] = 18 f_pi^2 comes entirely from the soft breaking, and STr[M^4] is what drives the CW corrections. The precise sign and magnitude depend on the full block structure.

**Assessment:** This calculation is doable but has not been done cleanly for the off-diagonal sector. If delta m^2_CW is negative and larger than f_pi^2 = 8464 MeV^2 for the ds sector but not the ud sector, it would produce tachyonic condensation with a mass hierarchy tracking the seesaw. The hierarchy would automatically be M_u^2 : M_d^2 : M_s^2 because the coupling W_{X,M^a_b,M^b_a} = M_k is flavor-dependent. The condensate ratio would still give tan(theta_C) = sqrt(m_d/m_s) by the same algebraic argument as in the B-term derivation.

**Difficulty:** Straightforward but numerically demanding (need the full 16 x 16 fermion and 32 x 32 scalar mass matrices). The real question is the sign.

### Alternative 2: Higher-dimensional Kahler operators

Non-canonical Kahler terms of the form

    delta K = c_1 / Lambda_UV^2 * |M^d_s|^2 |M^s_d|^2 + ...

generate quartic scalar couplings in the potential. At leading order, the Kahler correction to the off-diagonal meson mass-squared is

    delta m^2 ~ c_1 * M_k^2 / Lambda_UV^2

where Lambda_UV is the UV scale where the Kahler correction is generated. For c_1 ~ O(1) and Lambda_UV ~ Lambda = 300 MeV, this gives delta m^2 ~ M_k^2 / Lambda^2, which for k = u is (408471/300)^2 ~ 1.85 x 10^6 MeV^2. This is large and positive, stabilizing the off-diagonal directions further.

Negative Kahler corrections (c_1 < 0) could destabilize, but there is no generic reason to expect c_1 < 0. In a strongly coupled confining theory, the sign of c_1 is unknown.

**Assessment:** This is the least predictive option. The Kahler potential is not protected by holomorphy; any operator consistent with the symmetries can appear with unknown coefficient. Without a UV completion that computes c_1, this mechanism is unfalsifiable.

**Difficulty:** Requires new ideas (or lattice SQCD data).

### Alternative 3: Quark mass texture -- theta_C as an input relation

This is the most honest option. The framework after Koide constraints has 3 free parameters: m_s, m_u, m_d (or equivalently m_s, theta_C, m_u/m_d via the Oakes relation). The Koide seed + v_0-doubling + (c,b,t) overlap determine m_c and m_b from m_s. The top mass m_t is fixed by the top Yukawa at tan(beta) = 1. That leaves {m_u, m_d, m_s} free.

The Oakes relation tan(theta_C) = sqrt(m_d/m_s) has been known since 1969. It holds to 3.3% (Weinberg-Oakes) or 0.9% (GST variant sin(theta_C) = sqrt(m_d/m_s)). This is an empirical fact of the Standard Model, not a prediction of this framework.

**Parameter counting:** The model has 3 free quark mass parameters. The quark sector has 6 masses + 3 CKM angles + 1 CP phase = 10 observables. The Koide constraints fix 3 ratios (m_c/m_s, m_b/m_s via v_0-doubling, m_t from overlap or top Yukawa). That gives 10 - 3 = 7 remaining observables for 3 parameters, so there are 4 predictions. These predictions are: (i) m_c = 1301 MeV, (ii) m_b = 4177 MeV, (iii) the dual Koide Q(1/m_d, 1/m_s, 1/m_b) = 2/3, and (iv) any CKM relation derived from {m_u, m_d, m_s} alone.

If the Oakes relation is taken as empirical, then it is a consistency check: given m_d and m_s, the Cabibbo angle is determined by a well-known formula. The sBootstrap does not add anything beyond the observation that the seesaw M_j = C/m_j algebraically ensures that any off-diagonal condensate ratio reproduces sqrt(m_d/m_s). But without a mechanism that actually generates the condensate, the seesaw provides the scaffolding for Oakes, not its derivation.

**Is this a prediction, a consistency check, or neither?** It is a consistency check. The model has 3 free parameters and claims to describe a sector with 10 observables. With 3 Koide-type constraints, 4 predictions remain. If one of those predictions is the Oakes relation -- which was already known for 57 years -- it cannot honestly be called a new prediction. It can be called a structural consistency: the seesaw naturally produces the mass-ratio scaling that underlies Oakes. This is meaningful but not predictive.

The correct way to present this: "The Seiberg seesaw, which is required for the Koide mass predictions, automatically encodes the Oakes relation in the off-diagonal meson sector. This does not add new predictive power, but it shows that the CKM structure is compatible with -- and natural in -- the seesaw framework."

---

## (b) The Oakes coincidence

### Explicit look-elsewhere calculation

**Simple dimensionless ratios from {m_u, m_d, m_s}:**

Taking ratios, square roots, and cube roots of the three light quark masses gives these simple dimensionless quantities:

    m_d/m_s = 0.0500
    m_u/m_d = 0.4625
    m_u/m_s = 0.02313
    sqrt(m_d/m_s) = 0.2236
    sqrt(m_u/m_d) = 0.6801
    sqrt(m_u/m_s) = 0.1521
    (m_d/m_s)^{1/3} = 0.3684
    (m_u/m_d)^{1/3} = 0.7734
    (m_u/m_s)^{1/3} = 0.2853

That gives 9 simple ratios. Adding mixed products like m_u m_d / m_s^2 or (m_u/m_s)^{1/2} etc. quickly doubles the list. A generous estimate is ~15 independent simple ratios.

**CKM-related observables:**

    |V_us| = 0.2253
    |V_cb| = 0.0408
    |V_ub| = 0.00382
    sin^2 theta_W = 0.2312 (MS-bar) or 0.2231 (on-shell)
    J = 3.18 x 10^{-5}

That gives ~5 targets.

**Number of trials:** 15 x 5 = 75 combinations to check.

**Probability of accidental match to 1%:** Each dimensionless ratio is a random number on [0, 1] (after suitable normalization). A "match to 1%" means the fractional difference |r_1 - r_2|/r_2 < 0.01. If both r_1 and r_2 are in the range [0.05, 1], the probability of a match to 1% is approximately 2 x 0.01 = 0.02 (the target occupies a window of relative width 2% centered on r_2). For ratios near zero, the probability is higher because the absolute window is smaller but the relative criterion is easier to satisfy.

For the Oakes case: sqrt(m_d/m_s) = 0.2236, |V_us| = 0.2253. The fractional deviation is |0.2236 - 0.2253|/0.2253 = 0.0075, which passes a 1% criterion. However, if we use the less precise Weinberg-Oakes form (theta_C rather than sin(theta_C)), the deviation is 3.3%, which fails.

**P(at least one match in N trials) = 1 - (1 - p)^N:**

With N = 75 trials and p = 0.02 per trial:

    P(at least one) = 1 - (1 - 0.02)^75 = 1 - 0.98^75 = 1 - 0.218 = 0.782

With N = 75 trials and p = 0.01 per trial (stricter criterion):

    P(at least one) = 1 - 0.99^75 = 1 - 0.471 = 0.529

So the probability of finding at least one match to 1% accuracy among ~75 mass-ratio/CKM combinations is between 53% and 78%.

**The Oakes relation, taken in isolation, is not statistically significant.** Any scan of O(100) simple dimensionless ratios against O(5) CKM observables will produce at least one 1%-level match about half the time.

**What makes it non-trivial is not the numerics but the physics.** The Oakes relation has a physical derivation: in the limit of two light quarks mixing with a heavier one through a mass matrix of the form [[m_d, delta], [delta, m_s]], the mixing angle satisfies tan(theta) = delta/(m_s - m_d). If delta = sqrt(m_d * m_s) (geometric mean), then tan(theta) = sqrt(m_d/m_s). This is the Fritzsch texture zero, and it follows from the assumption that the 12 element of the down-type mass matrix is sqrt(m_d * m_s). The Seiberg seesaw provides exactly this structure, since the off-diagonal condensate at the geometric mean of the diagonal VEVs gives delta_quark = C/sqrt(m_d m_s) -> sqrt(m_d * m_s) after rescaling.

**Conclusion:** The numerical coincidence alone has P ~ 50-80% of being accidental. The physics (seesaw structure -> Fritzsch texture -> Oakes) is what gives it meaning. If the seesaw is real, the Oakes relation is a consequence. If the seesaw is not real, the Oakes relation is a moderately common numerological coincidence.

---

## (c) What does X != S mean for the Higgs sector?

If the Seiberg Lagrange multiplier X is NOT the NMSSM singlet S, the following questions arise.

### 1. What provides the Higgs quartic coupling?

Three options:

**(a) MSSM radiative correction from stops.** At tan(beta) = 1, the tree-level MSSM Higgs mass is zero (the D-term quartic vanishes along the D-flat direction |H_u| = |H_d|). The one-loop stop correction is

    delta m_h^2 = (3 y_t^4 v^2 / 32 pi^2) sin^4(beta) * ln(m_stop^2 / m_t^2)

At sin(beta) = 1/sqrt(2) (tan(beta) = 1), the sin^4(beta) = 1/4 suppression means m_stop ~ 2 x 10^8 GeV is required for m_h = 125 GeV. This is unphysical.

**(b) D-terms from an extended gauge group.** If SU(2)_L x U(1)_Y is embedded in a larger group (e.g., SU(2)_L x SU(2)_R x U(1)_{B-L}), additional D-term quartics arise. At tan(beta) = 1, these can provide the needed quartic if the extended gauge symmetry breaking scale is O(TeV). However, this introduces new gauge bosons and significantly enlarges the model.

**(c) A separate NMSSM singlet S.** A field S unrelated to X, with superpotential W = lambda_S S H_u . H_d + (kappa/3) S^3 and lambda_S = 0.72, gives m_h = lambda_S v/sqrt(2) = 125.4 GeV at tree level. The field S is a garden-variety NMSSM singlet. Its VEV provides the mu-term: mu = lambda_S <S>. Since S is unrelated to the SQCD sector, <S> is a free parameter set by the NMSSM soft terms.

**Assessment:** Option (c) is the only viable one at tan(beta) = 1 without exotic gauge structure. The MSSM with stops cannot do it; D-terms from extended gauge symmetry are possible but over-engineered.

### 2. What provides the mu-term?

If S is a separate singlet, the mu-term is mu = lambda_S <S>, with <S> determined by the NMSSM soft potential for S. This is the standard NMSSM solution to the mu-problem: mu ~ m_soft ~ O(100 GeV to 1 TeV). The numerical value of mu is a free parameter.

Alternative: Giudice-Masiero mechanism. A Kahler term K = c (H_u . H_d)^dagger / M_Pl would generate mu ~ F / M_Pl after SUSY breaking. With F ~ m_soft^2 ~ (1 TeV)^2, mu ~ 10^{-12} GeV, far too small. This fails.

Another alternative: Radiative generation. If S has a soft mass m_S^2 > 0 at the UV scale, RG running with the top Yukawa can drive m_S^2 < 0, triggering <S> != 0. This is the standard radiative EWSB mechanism applied to the singlet. It works but does not predict <S>.

### 3. How many additional parameters does separating X from S introduce?

The X = S identification removed one field from the spectrum (S is identified with X, so no new field). Separating them re-introduces S as a new chiral superfield with its own Kahler metric, its own superpotential couplings (lambda_S, kappa), and its own soft terms (m_S^2, A_lambda, A_kappa).

New parameters: lambda_S (1), kappa (1), m_S^2 (1), A_lambda (1), A_kappa (1) = 5 new parameters.

Old parameters removed: lambda_0 (the X-H coupling) is no longer present, saving 1 parameter but introducing 5.

Net cost: 4 additional parameters.

Of these, lambda_S = sqrt(2) m_h / v = 0.72 is determined by the Higgs mass. kappa is constrained by stability and the singlino mass. m_S^2, A_lambda, A_kappa are soft parameters with O(TeV) natural scale. So effectively 3 free soft parameters are added.

### 4. Any residual connection between SQCD dynamics and the Higgs sector?

Yes, through the Yukawa couplings W_Yukawa = y_c H_u M^d_d + y_b H_d M^s_s. These couple the Higgs sector to the meson sector regardless of whether X = S. The ISS flavor-universality identity y_j <M_j> = 2C/v holds, ensuring that the Higgs F-terms are flavor-blind. The meson-Higgs couplings still transmit SUSY breaking from the meson sector to the Higgs sector, generating soft masses of order F/v ~ C/v^2 ~ 10^{-2} MeV (far too small to be the dominant SUSY-breaking source).

The Seiberg seesaw structure (M_j = C/m_j) still determines the pattern of meson VEVs, and through the Yukawa couplings, it still determines the effective quark masses at tree level: m_q^{eff} = y_q <M_q>. The quark mass predictions (m_c, m_b from Koide) are unaffected by the X/S separation.

What is lost: the elegant connection between the SQCD quantum constraint (enforced by X) and the Higgs quartic (provided by S). The Higgs mass prediction m_h = lambda v/sqrt(2) = 125.4 GeV survives because it depends on lambda_S, which is the same numerical value regardless of the identity of S. But the prediction is less powerful because lambda_S is now a parameter of the separate NMSSM sector rather than a consequence of the SQCD dynamics.

---

## (d) The mesino problem revisited

### X as a Lagrange multiplier: kinetic terms and the X-ino

For N_f = N_c, the Seiberg effective theory has X as a Lagrange multiplier with NO kinetic term in the classical Kahler potential. The classical K = |X|^2 term is absent because X is not a composite of the microscopic quarks; it is introduced to enforce the quantum-modified constraint det M - B B-tilde = Lambda^{2 N_c}.

However, at one loop, the meson loops generate a kinetic term for X through the wavefunction renormalization:

    K_eff = Z_X(M, M^dagger) |X|^2 + ...

where Z_X ~ (1/16 pi^2) sum_k M_k^2 / Lambda_UV^2. At the seesaw vacuum, Z_X ~ (1/16 pi^2) * M_u^2 / Lambda^2 ~ (1/16 pi^2) * (408471/300)^2 ~ 1.2 x 10^4. This is a large one-loop correction, reflecting the strong coupling of the confined theory.

The question of whether the X-ino (psi_X) is a propagating degree of freedom depends entirely on whether X has a kinetic term:

- **If Z_X = 0 (strictly a Lagrange multiplier):** The equation of motion for X is algebraic: det M - B B-tilde = Lambda^6. There is no X propagator, no psi_X fermion in the spectrum, and the fermion count drops by 1 (from 16 Weyl to 15). The constraint removes 1 complex scalar and 1 Weyl fermion simultaneously, as befits a supersymmetric constraint.

- **If Z_X != 0 (dynamically generated kinetic term):** X is a propagating field with mass m_X^2 ~ W_{XX} / Z_X. The psi_X fermion exists and mixes with the diagonal meson fermions and Higgsinos through the central 6x6 block. Its mass is ~ 7.73 x 10^{10} MeV (the heaviest state in the fermion spectrum, from the X-M_s^s mixing).

In the Seiberg framework, the convention is to INCLUDE X as a dynamical field with Z_X = 1, even though it is "technically" a Lagrange multiplier. This is because the quantum constraint is enforced by the dynamics (not by hand), and the X field's equation of motion is the constraint equation. The X-ino exists as a physical fermion.

### Fermion spectrum counting

With X dynamical (Z_X = 1):
- 16 Weyl fermions total (9 meson + 1 X + 2 baryon + 4 Higgs)
- Form 8 Dirac fermions (paired by the superpotential mass matrix)

Without X (Z_X = 0):
- 15 Weyl fermions (9 meson + 2 baryon + 4 Higgs)
- The constraint det M = Lambda^6 + B B-tilde removes one real equation on the scalar manifold and its fermionic partner

### The mesino problem with and without X

The off-diagonal meson fermions (mesinos) acquire their masses from the superpotential coupling

    W_{M^a_b, M^b_a} = X_0 * cofactor(M, a, b)|_{diag} = X_0 * M_k

where k is the third flavor. The mesino mass is |X_0| * M_k, which is proportional to X_0 = -C/Lambda^6 ~ -1.21 x 10^{-9} MeV^{-4}.

This structure exists regardless of whether X is dynamical or not. Even if X has no kinetic term, it still appears in the superpotential, and its VEV X_0 is determined by the F-term equations. The mesino masses |X_0| M_k are set by the algebraic value of X_0, not by the propagation of X.

**Distinction between "mesino" and "X-ino":**

- **Mesino** (psi_{M^a_b}): fermion in the M^a_b chiral supermultiplet. Off-diagonal mesinos have masses |X_0| M_k ~ 10^{-5} to 10^{-4} MeV. These are physical propagating fermions at the seesaw vacuum, regardless of the status of X.

- **X-ino** (psi_X): fermion in the X supermultiplet. If X is dynamical, psi_X has mass ~ 7.73 x 10^{10} MeV (from mixing with psi_{M^s_s}). If X is a Lagrange multiplier, psi_X does not exist.

The mesino catastrophe is NOT resolved by making X a Lagrange multiplier. The off-diagonal mesino masses come from W_{M^a_b, M^b_a}, which involves X_0 as a numerical coefficient, not as a propagating field. Whether X propagates is irrelevant to the mesino mass scale.

### What does resolve the mesino problem?

The off-diagonal mesino masses are

    m_mesino(a,b) = |X_0| * M_k = |C/Lambda^6| * C/m_k = C^2 / (Lambda^6 m_k)

These are holomorphic masses (from the superpotential second derivatives). The physical masses are

    m_phys = Z_M^{-1/2} * m_hol

where Z_M is the Kahler metric for the off-diagonal meson field. In a strongly coupled confined theory, Z_M is not calculable from first principles. For m_phys to be phenomenologically viable:

For the charged mesinos (Q = +/- 1) to avoid g-2 exclusion, we need m_phys > m_e = 0.511 MeV (absolute minimum). More conservatively, LEP constraints require m_phys > M_Z/2 ~ 45 GeV for new charged fermions.

The required Kahler suppression:

    Z_M > (m_hol / m_phys)^2

For the ud mesino: m_hol = 1.14 x 10^{-5} MeV, m_phys > 45000 MeV:
    Z_M > (1.14e-5 / 45000)^2 = 6.4 x 10^{-20}

Wait -- this is backwards. Z_M^{-1/2} amplifies the mass if Z_M < 1 (strong wavefunction renormalization). We need Z_M^{-1/2} * m_hol > m_phys, i.e., Z_M < (m_hol/m_phys)^2. For the ud mesino, Z_M < (1.14e-5 / 45000)^2 = 6.4 x 10^{-20}, i.e., the canonically normalized field has a Kahler metric of order 10^{-20}.

Alternatively: the physical mass is m_phys = m_hol / sqrt(Z_M). To get m_phys = 45 GeV from m_hol = 1.14 x 10^{-5} MeV:

    sqrt(Z_M) = m_hol / m_phys = 1.14e-5 / 45000 = 2.5 x 10^{-10}
    Z_M = 6.4 x 10^{-20}

This is an extremely small Kahler metric, corresponding to an enormous anomalous dimension for the composite operator Q^u Q-bar_d. In a confining theory, Z_M ~ (Lambda/M_{UV})^{2 gamma} where gamma is the anomalous dimension. For Lambda = 300 MeV and M_{UV} = M_Pl:

    Z_M ~ (300 / 2.4e21)^{2 gamma} = (1.25e-19)^{2 gamma}

To get Z_M = 6.4 x 10^{-20}: 2 gamma * log(1.25e-19) = log(6.4e-20), so 2 gamma * (-43.5) = -44.2, giving gamma = 0.508. An anomalous dimension of 1/2 is natural at a conformal fixed point.

This is actually not unreasonable. SQCD near the conformal window (N_f ~ 3 N_c / 2) has anomalous dimensions gamma ~ 1/2. For N_f = N_c = 3, we are below the conformal window, but the theory is strongly coupled. The anomalous dimension at the confining scale could plausibly be O(1/2).

**The mesino problem is not solved, but it may not be as catastrophic as it appears.** If the anomalous dimension of the off-diagonal meson operators is gamma ~ 1/2, the physical mesino masses could be O(TeV) rather than O(eV). This requires a computation in strongly coupled SQCD that is beyond current analytic control.

---

## (e) Five most important calculations to do next

### Calculation 1: Off-diagonal meson mass matrix at the combined SQCD + EWSB vacuum

**What:** Start from the complete superpotential W = W_Seiberg + W_Yukawa + W_NMSSM (with S as a separate singlet if needed). Set H_u = H_d = v/sqrt(2). At the EWSB vacuum, the F-term F_X picks up a correction from the NMSSM coupling: F_X = det M - Lambda^6 + lambda_S <S> v^2/2. With the Seiberg vacuum det M = Lambda^6, this gives F_X = lambda_S <S> v^2/2. If X = S, this is F_X = lambda X_0 v^2/2 = -lambda C v^2 / (2 Lambda^6). Compute the full 6x6 off-diagonal meson mass matrix including the soft term f_pi^2 and the F_X-driven B-term. Determine whether any eigenvalue is negative.

The critical comparison is:

    m^2_tach = -|F_X|^2 * M_k^2 vs. m^2_soft = f_pi^2

For the ds sector (k = u): |F_X|^2 * M_u^2 = (lambda C v^2 / 2 Lambda^6)^2 * M_u^2. With C = 882297, v = 246220, Lambda = 300, lambda = 0.72:

    F_X = 0.72 * 882297 * 246220^2 / (2 * 300^6) = 0.72 * 882297 * 6.06e10 / (2 * 7.29e14)
         = 0.72 * 5.35e16 / 1.458e15 = 0.72 * 36.7 = 26.4 MeV

Wait, let me redo this carefully. The dimensions must track properly.

X_0 = -C/Lambda^6 = -882297 / (300^6) = -882297 / 7.29e14 = -1.210e-9 MeV^{-4}

F_X at the EWSB vacuum with NMSSM coupling: the constraint equation is modified to det M - B B-tilde - Lambda^6 + lambda X H_u . H_d = 0. At the seesaw with det M = Lambda^6 and B = 0:

    F_X = lambda * (H_u . H_d) = lambda * v^2/2

But F_X = dW/dX and W has [X] such that X (det M) has mass dimension. Since det M ~ [mass]^3 (mesons have mass dimension 1 in the confined theory -- actually, [M] = mass^2 since M = Q Q-bar/Lambda has dimension 2 in mass units if Q has dimension 3/2 in the UV)... Let me not redo the dimensional analysis here. The key numerical comparison was done in ckm_analytic.md: F_X = -lambda v^2/2 = -2.18e10 MeV^5.

The tachyonic mass-squared for the ds sector is |F_X|^2 * M_u^2 / Lambda^{12} (from the B-term coupling W_{X, M^d_s, M^s_d} which involves the cofactor). From ckm_analytic.md: m^2_tach(ds) = -7.95e31 MeV^2, while f_pi^2 = 8464 MeV^2.

The B-term completely overwhelms the soft mass. If this computation is correct, there are tachyons.

But ckm_perturbative.md says: at the exact seesaw vacuum with H = 0, F_X = 0 and all off-diagonal masses are positive (= f_pi^2). The discrepancy is that the two computations use different vacua. The ckm_analytic.md computation uses the EWSB vacuum (H != 0) and takes F_X from the NMSSM coupling. The ckm_perturbative.md computation uses the pre-EWSB vacuum (H = 0) where F_X = 0 by the seesaw.

**The reconciliation:** The physical vacuum is the EWSB one (H = v/sqrt(2)). At this vacuum, F_X != 0 if the NMSSM coupling is present. The computation in ckm_perturbative.md is correct for the pre-EWSB vacuum but not the physical one. The computation in ckm_analytic.md is correct for the EWSB vacuum but makes approximations about the meson condensation dynamics.

What is needed is a single computation that:
1. Sets H = v/sqrt(2) from the start
2. Computes the meson vacuum (diagonal + off-diagonal) self-consistently including both the soft term and the NMSSM-induced F_X
3. Diagonalizes the resulting 6x6 mass matrix for the off-diagonal mesons
4. If tachyons exist, computes the condensate and extracts theta_12

**Expected result:** Tachyons in the ds and us sectors, with condensate ratio epsilon_us/epsilon_ud = sqrt(m_d/m_s), giving theta_C ~ 12.6 deg (Weinberg-Oakes).

**Cost if negative:** If no tachyons (because the F_X contribution is smaller than f_pi^2 after self-consistent treatment), then the Cabibbo angle is NOT derived from the seesaw. The model would have theta_C as a free parameter, equivalent to m_d/m_s.

**Difficulty:** Straightforward. This is a classical minimization problem with 18 real variables (9 complex meson entries). The key is to do it at the EWSB vacuum, not the pre-EWSB one.

### Calculation 2: Off-diagonal mesino physical masses from anomalous dimensions

**What:** Compute the anomalous dimension of the off-diagonal meson operator M^a_b = Q^a Q-bar_b / Lambda in the confining regime. Use the NSVZ exact beta function for N_f = N_c = 3 to extract gamma_M at the confining scale. Then compute the physical mesino mass m_phys = m_hol * Lambda^{gamma_M} / M_UV^{gamma_M}, where m_hol = |X_0| M_k are the holomorphic masses.

The NSVZ formula relates the anomalous dimension to the beta function:

    beta(g) = -g^3/(16 pi^2) * [3 N_c - N_f(1 - gamma)] / [1 - N_c g^2/(8 pi^2)]

At the confining scale (g -> infinity), the anomalous dimension approaches the unitarity bound gamma -> 1 (for free magnetic quarks in the Seiberg dual). For N_f = N_c, there is no magnetic dual, so the anomalous dimension at confinement is less constrained. Lattice studies of N_f = N_c = 3 SQCD would be the definitive input.

**Expected result:** gamma_M ~ 0.3 to 0.7 for off-diagonal mesons. This would give Z_M ~ 10^{-15} to 10^{-25}, lifting the mesino masses from eV to keV-MeV range. Whether this reaches the LEP-safe range (> 45 GeV) depends on the precise value of gamma.

**Cost if negative:** If gamma_M is too small (< 0.3), the physical mesino masses remain below 1 MeV and the model is catastrophically excluded. This would be a fatal problem requiring either identification of mesinos with SM particles or removal from the spectrum by additional dynamics.

**Difficulty:** Requires new ideas. The anomalous dimension in a confining theory away from the conformal window is not analytically computable. Lattice SQCD or holographic models would be needed.

### Calculation 3: Lepton sector mediation mechanism

**What:** Construct an explicit coupling between the Sp(2) confining sector (which produces the lepton Koide seed) and the SM lepton Yukawa sector. The Sp(2) theory has mesons Phi^a_b with VEVs that satisfy Q = 2/3. These mesons must couple to L . H_d . e_R through some mediator field. Identify the mediator, compute the lepton Yukawa couplings, and derive the overall lepton mass scale M_0 = 313.84 MeV^{1/2} from the Sp(2) parameters (Lambda_L, mu_L).

**Expected result:** A superpotential of the form W_med = sum_a y_a^L L_a . H_d . e_R / M_* with y_a^L = Phi_a / M_*, giving m_l_a = <Phi_a> v / (sqrt(2) M_*). If M_* ~ 10^5 GeV and Lambda_L ~ O(GeV), this could work.

**Cost if negative:** If no natural mediator exists, the lepton sector must be treated as phenomenological (3 free lepton masses). This does not kill the quark sector predictions but weakens the claim of a unified SUSY framework.

**Difficulty:** Requires new ideas. The mediator must be a heavy field that couples to both the Sp(2) sector and the SM leptons. No such field has been identified in the current framework.

### Calculation 4: Full one-loop CW correction to the off-diagonal meson mass matrix

**What:** Compute the one-loop Coleman-Weinberg effective potential for the off-diagonal meson entries at the seesaw vacuum, including all 16 chiral superfields and the soft breaking. The CW potential is

    V_CW = (1/64 pi^2) STr[M^4(log(M^2/Lambda^2) - 3/2)]

evaluated as a function of the off-diagonal entries epsilon^a_b. The mass-squared correction is d^2 V_CW / d epsilon^a_b d epsilon^{b*}_a at epsilon = 0.

This is distinct from Calculation 1 (which uses the B-term mechanism at the EWSB vacuum). Here we stay at the pre-EWSB seesaw vacuum and ask whether quantum corrections alone can destabilize the off-diagonal directions.

**Expected result:** The CW correction is positive (stabilizing) because the spectrum at the seesaw vacuum is non-tachyonic and the supertrace STr[M^2] = 18 f_pi^2 > 0. Positive STr[M^2] generally implies positive CW masses (by the supertrace sum rule). If so, the CW mechanism does not generate CKM mixing, reinforcing the conclusion that the EWSB-induced F_X is the necessary ingredient.

**Cost if negative:** If the CW correction is negative and large enough to overcome f_pi^2, it would provide a CKM mechanism without requiring the NMSSM coupling. This would be a significant simplification but would also decouple the CKM from the Higgs sector.

**Difficulty:** Straightforward. Requires diagonalizing the 32x32 scalar mass matrix and the 16x16 fermion mass matrix as functions of the off-diagonal entries, then computing the supertrace-weighted sum. Numerically demanding but conceptually standard.

### Calculation 5: Complete Yukawa Lagrangian with identified fermion spectrum

**What:** Assemble the full Yukawa coupling matrix connecting the meson fermion spectrum to the Higgs sector. At the seesaw vacuum, the fermion mass eigenstates are mixtures of diagonal mesinos and Higgsinos (from the central 6x6 block) plus off-diagonal mesinos (from the 2x2 blocks). Identify which mass eigenstates correspond to which SM quarks and leptons (if possible). Write the Yukawa matrix Y^{u,d}_{ij} in the mass eigenstate basis. Compute the CKM matrix as V_CKM = U_L^{u dagger} U_L^d, where U_L^{u,d} are the left rotations that diagonalize the up-type and down-type Yukawa matrices.

This is the ultimate deliverable: a Yukawa Lagrangian with numerical entries that reproduces the SM quark masses and CKM matrix (at least at leading order) from the seesaw framework.

**Expected result:** The Yukawa matrix has a hierarchical structure set by the seesaw inversion. The up-type Yukawa matrix is approximately diagonal (charm and top masses from y_c and y_t). The down-type Yukawa matrix has off-diagonal entries from the meson condensate, giving a CKM matrix with theta_12 ~ arctan(sqrt(m_d/m_s)) ~ 12.6 deg, theta_23 ~ 0 (no mechanism for V_cb at leading order), theta_13 ~ 0.

**Cost if negative:** If the fermion identification fails (quantum numbers don't match, masses are wrong by orders of magnitude), the seesaw-to-Yukawa program is incomplete. The most likely failure mode is the mesino masses being too light (the same problem as in (d)).

**Difficulty:** Straightforward in principle, but requires resolving the fermion identification problem (which SM fermions correspond to which meson fermion eigenstates). The diagonal meson-Higgsino mixtures at 10^{-4} to 10^{-2} MeV are 4-5 orders of magnitude below the physical quark masses, which means additional dynamics (Kahler corrections, or the bloom mechanism) must be invoked to reach the physical spectrum. This makes the "straightforward" calculation contingent on solving the non-perturbative Kahler problem.

### Priority ranking

1. **Calculation 1** (off-diagonal mass matrix at EWSB vacuum) -- highest priority because it resolves the tachyon vs perturbative vacuum conflict, which has been the central open question since Round 12. It is also the most straightforward.

2. **Calculation 5** (Yukawa Lagrangian) -- the North Star goal of the project. Even a partial result (tree-level Yukawa matrix with identified fermion content) would be the most important deliverable for the paper.

3. **Calculation 2** (mesino anomalous dimensions) -- critical for phenomenological viability. If the mesino masses are wrong by 10 orders of magnitude, nothing else matters. But the calculation may require non-perturbative input that is currently unavailable.

4. **Calculation 4** (CW correction to off-diagonal masses) -- important as a cross-check on Calculation 1. If the CW correction is the dominant effect (rather than the B-term), the CKM mechanism is different and more robust (it does not depend on the NMSSM coupling).

5. **Calculation 3** (lepton mediation) -- important for completeness but not for the quark sector predictions. Can be deferred.

---

## Summary

The framework has three sharp results and three open problems:

**Sharp results:**
- The Koide seed + v_0-doubling predicts m_b = 4177 MeV (0.1 sigma from PDG)
- The Seiberg seesaw algebraically encodes the Oakes relation tan(theta_C) = sqrt(m_d/m_s)
- The NMSSM coupling lambda = 0.72 gives m_h = 125.4 GeV at tree level for tan(beta) = 1

**Open problems:**
- The CKM mechanism (tachyonic condensation vs perturbative stability) is unresolved because the two vacuum computations use different starting points (pre-EWSB vs post-EWSB). Calculation 1 is the priority.
- The mesino masses at the seesaw vacuum are 10 orders of magnitude below experimental bounds for charged fermions. The resolution may come from large anomalous dimensions (gamma ~ 1/2), but this is not calculable with current techniques.
- The X = S identification fails for the mu-term (mu_eff = -7 MeV vs ~100 GeV needed). Separating X from S costs 4 additional parameters.

The Oakes relation, taken in isolation, has a look-elsewhere probability of ~50-80% of being accidental (among ~75 simple mass-ratio/CKM comparisons). Its significance comes from the physical mechanism (seesaw -> Fritzsch texture) rather than the numerical coincidence.

---

*Generated: 2026-03-04*
*Based on: complete_lagrangian_opus.md, ckm_analytic.md, ckm_perturbative.md, brainstorm_r13.md, mesino_phenomenology.md, offdiag_vacuum.md, nmssm_spectrum.md, higgs_potential.md, full_spectrum.md, bloom_bion_potential.md, iss_cw_koide.md, cabibbo_oakes.md, v0_doubling.md, bloom_mechanism.md, ckm_koide.md, baryon_stabilization.md*
