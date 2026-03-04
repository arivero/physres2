# Round 17 Brainstorm: Gap Analysis for a Complete Theory

**Date:** 2026-03-04

---

## (a) Missing Pieces Inventory

Every claim in the paper that is not yet supported by explicit computation, classified as (i) straightforward to compute, (ii) requires non-perturbative input, or (iii) is an assumption.

### A1. Superpotential and Vacuum

| Claim | Status | Classification | What's missing |
|-------|--------|----------------|----------------|
| Seiberg seesaw M_j = C/m_j | COMPUTED | -- | Nothing; exact from quantum constraint |
| F_X = 0 at seesaw vacuum | COMPUTED | -- | Verified numerically |
| det M = Lambda^6 satisfied | COMPUTED | -- | Exact to machine epsilon |
| X is Lagrange multiplier, not S | COMPUTED | -- | Factor 171,500 incompatibility established |
| Baryon mass term m_B ~ Lambda eliminates flat direction | COMPUTED | -- | S_bounce > 10^7 |
| Three-instanton term W_inst = c_3 (det M)^3/Lambda^18 exists | NOT COMPUTED from first principles | (ii) non-perturbative | The coefficient c_3 is O(1) by standard instanton calculus, but its precise value controls the bloom dynamics. Computing it requires a semiclassical instanton calculation in SU(3) SQCD at N_f = 3. |
| The NMSSM coupling lambda X (H_u . H_d) with X = S | KILLED | (iii) assumption, now retracted | Must use separate S field. The complete S sector (kappa S^3, soft terms for S) is not specified. |

### A2. Kahler Potential

| Claim | Status | Classification | What's missing |
|-------|--------|----------------|----------------|
| Canonical Kahler K = Tr(M-dag M) + ... | ASSUMED | (iii) assumption | In principle determined by ISS duality, but higher-order corrections not computed |
| Pseudo-modulus Kahler K_pm = -|X|^4 / (12 mu^2) | NOT DERIVED | (iii) assumption / fitted | c = -1/12 is the value that places the Kahler pole at the seesaw VEV. A microscopic derivation from the one-loop effective Kahler potential in the ISS magnetic dual is needed but has not been done. |
| Bion Kahler K_bion = (zeta^2/Lambda^2)|sum s_k sqrt(M_k)|^2 | NOT DERIVED | (ii) non-perturbative | The functional form is motivated by the Unsal semiclassical framework (monopole-instanton pairs on R^3 x S^1), but the specific sqrt(M) dependence is non-standard. No rigorous derivation from the compactified theory exists. |
| Sign structure s_k = (-1, +1, +1) for bloom | ASSUMED | (iii) assumption | Assigned to match the physical (-s, c, b) sign pattern. No derivation from monopole zero-mode structure. |
| Kahler metric for off-diagonal mesons is canonical | ASSUMED | (ii) non-perturbative | The physical mesino masses depend on Z_M, which is not protected by holomorphy. In a strongly coupled theory, Z_M can differ from 1 by orders of magnitude. |
| Kahler expansion breaks down at seesaw vacuum | NOTED but not quantified | (i) straightforward | The |X|^4 expansion is valid when |X| << mu. At the seesaw vacuum |X| = sqrt(3) mu (the pole). The expansion fails. A resummation or exact Kahler metric is needed. |

### A3. SUSY-Breaking and Soft Terms

| Claim | Status | Classification | What's missing |
|-------|--------|----------------|----------------|
| Soft mass m-tilde^2 = f_pi^2 = (92 MeV)^2 for all 9 mesons | ASSUMED | (iii) assumption | The identification of the SUSY-breaking scale with f_pi is assumed. The mechanism (mediation) is not specified. |
| Higgs soft masses constrained by tan beta = 1 | CONSTRAINED but not derived | (iii) assumption | The relation m_{H_u}^2 - m_{H_d}^2 = y_b^2 - y_c^2 is a consistency condition, not a derivation. tan beta = 1 itself is assumed. |
| No gaugino masses M_1, M_2 specified | MISSING | (iii) required but absent | Any complete SUSY model needs gaugino masses for phenomenology. Not addressed. |
| No A-terms specified | MISSING | (iii) required but absent | Trilinear A-terms are generically present. Set to zero implicitly. |
| AMSB most promising for mesino mass lifting | STATED | (i) straightforward to compute | The AMSB contribution to the meson soft masses is m_soft ~ (beta/g)^2 m_{3/2}. With m_{3/2} ~ 50 TeV (typical AMSB), this gives O(TeV) scalar masses. The fermion masses are not lifted by AMSB directly. A quantitative AMSB spectrum has not been computed. |

### A4. Koide and Mass Relations

| Claim | Status | Classification | What's missing |
|-------|--------|----------------|----------------|
| Koide seed from O'Raifeartaigh at gv/m = sqrt(3) | COMPUTED | -- | Algebraically exact |
| v_0-doubling: sqrt(m_b) = 3 sqrt(m_s) + sqrt(m_c) | COMPUTED | -- | m_b = 4177 MeV, 0.1 sigma |
| v_0-doubling is NOT statistically significant in isolation | COMPUTED | -- | Look-elsewhere probability 24.9% |
| Q(-s,c,b) = 2/3 | OBSERVED | -- | Deviation 1.24% from 2/3 |
| Q(c,b,t) = 2/3 | OBSERVED | -- | Deviation 0.42% from 2/3 |
| Q(e,mu,tau) = 2/3 | OBSERVED to 0.91 sigma | -- | 2.8 sigma after look-elsewhere |
| Dual Koide Q(1/m_d, 1/m_s, 1/m_b) = 0.665 | OBSERVED | -- | No mechanism in the superpotential |
| delta_0 mod 2pi/3 = 2/9 | OBSERVED to 35 ppm | -- | 3.1 sigma look-elsewhere |
| Second scaling M_2 = 9 M_0 fails | COMPUTED | -- | Kills the M_n = 3^n M_0 extrapolation |
| Cyclic echo t -> c -> s -> stall -> c exact | COMPUTED | -- | Machine epsilon |

### A5. CKM and Flavor

| Claim | Status | Classification | What's missing |
|-------|--------|----------------|----------------|
| Tree-level Yukawa matrices are diagonal | ASSUMED | (iii) assumption | Follows from the seesaw at a diagonal vacuum. CKM must come from UV. |
| CKM from UV Fritzsch texture | PROPOSED | (i) straightforward | The explicit Fritzsch texture with Koide-constrained eigenvalues has not been written down and diagonalized to extract CKM angles. Needs the UV mass matrix. |
| Oakes relation sin theta_C = sqrt(m_d/m_s) = 0.224 | OBSERVED | -- | 0.93% from V_us. Not a prediction (m_d free). |
| Block-diagonal obstruction (no CKM from mesons) | PROVEN | -- | Theorem about polynomial Kahler invariants |
| All Koide triples mix up/down types | COMPUTED | -- | Structural feature, 100% rate among Q ~ 2/3 triples |

### A6. Lepton Sector

| Claim | Status | Classification | What's missing |
|-------|--------|----------------|----------------|
| Lepton masses from perturbative mesinos | FAILED | -- | Mass ratio m_s/m_d = 20 vs m_mu/m_e = 207; excluded |
| Sp(2) confining sector for leptons | PROPOSED | (ii) non-perturbative | The Sp(2) O'Raifeartaigh-Koide model is outlined but: (a) the ADS superpotential has not been verified to produce gv/mu = sqrt(3), (b) the bloom potential has not been computed, (c) the mediation to SM leptons is unspecified, (d) the mediation scale M_* ~ 43.5 GeV is uncomfortably low |
| QED running preserves Q = 2/3 below Lambda | COMPUTED (estimate) | (i) straightforward | delta_Q ~ 10^{-5}, consistent with observed 9 ppm deviation. A full numerical RG integration has not been done. |

### A7. Spectrum and Phenomenology

| Claim | Status | Classification | What's missing |
|-------|--------|----------------|----------------|
| Complete fermion spectrum (16 eigenvalues) | COMPUTED | -- | Table in full_spectrum.md |
| Tachyonic scalar modes (us, ds, central) | COMPUTED | -- | Indicates diagonal vacuum is not the true minimum |
| Off-diagonal meson scalars at m ~ 92-130 MeV | COMPUTED | -- | 8x too light for K-Kbar FCNC bound |
| Supertrace STr = 18 f_pi^2 | COMPUTED analytically | -- | Exact cancellation of F-term and holomorphic contributions |
| Higgs mass m_h = 125 GeV at tree level | COMPUTED | -- | lambda = 0.72 at tan beta = 1 (but requires loop corrections for precision) |
| m_h loop corrections (stop/sbottom) | NOT COMPUTED | (i) straightforward | The standard MSSM/NMSSM loop correction delta m_h^2 ~ (3 y_t^4 v^2)/(8 pi^2) ln(m_stop^2/m_t^2) has not been evaluated for this model. Since there are no stops (top is elementary, not confined), the loop corrections come from the meson/Higgs sector and are suppressed. |
| ISS F-term universality F_{H_u} = F_{H_d} = C/v | COMPUTED algebraically | -- | Exact identity from the seesaw + Yukawa structure |
| Mesino g-2 exclusion | COMPUTED | -- | Excluded by 7 orders of magnitude for charged mesinos at 11-229 eV |

### A8. Group Theory and Bootstrap

| Claim | Status | Classification | What's missing |
|-------|--------|----------------|----------------|
| Bootstrap uniqueness (N=3, r=3, s=2) | COMPUTED | -- | Unique solution |
| 24+15+15-bar anomaly-free | COMPUTED | -- | A = 0 verified |
| 24+15+15-bar asymptotically free (b_0 = 31/3) | COMPUTED | -- | |
| SO(32) adjoint 496 decomposition | COMPUTED via two chains | -- | Full table verified |
| de Vries angle R = 0.2231 matches sin^2 theta_W | COMPUTED | -- | On-shell value, 0.13 sigma |
| 15 vs 10-bar conflict (Pauli) | IDENTIFIED but not resolved | (ii) non-perturbative | Pauli requires antisymmetric 10-bar for J=0 color-3-bar diquarks, but sBootstrap needs symmetric 15. Resolution may involve P-wave diquarks or 't Hooft vertex effects. |

---

## (b) The FCNC Problem

### The numbers

The off-diagonal meson scalars (M^d_s, M^s_d, M^u_s, M^s_u, M^u_d, M^d_u) have two mass scales in the model:

1. **At the seesaw vacuum**: m^2 = f_pi^2 = 8464 MeV^2, giving m = 92 MeV
2. **At the M = 0 vacuum**: m^2 = 2 f_pi^2 = 16928 MeV^2, giving m = 130 MeV

The K-Kbar mixing constraint (UTfit, Delta S = 2) requires:

| Coupling | Required M_FCNC | Available m_S | Tension factor |
|----------|----------------|---------------|----------------|
| g = 1 (unit coupling) | 745 MeV | 92 MeV | 8x |
| g = V_ds = 0.225 (CKM) | 167 MeV | 92 MeV | 1.8x |
| g = NDA = 7.3 (naive dimensional analysis) | 5408 MeV | 92 MeV | 59x |

### The key question: are the off-diagonal mesons "new particles"?

No. In QCD, the off-diagonal meson scalars M^d_s and M^s_d are the K^0 and K-bar^0 kaons (or more precisely, their scalar partners -- the a_0(980) and K*_0(700) nonet). The physical kaons have m_K = 498 MeV, well above the FCNC bound for CKM-suppressed couplings. The SUSY scalar partners in the model have mass m = 92 MeV, far below the physical kaon mass.

The issue is that the SUSY model treats the soft mass f_pi^2 as the only scalar mass contribution. In QCD, the kaon mass receives contributions from:

1. **Chiral symmetry breaking**: m_K^2 = B_0 (m_d + m_s) with B_0 = m_pi^2/(m_u + m_d) ~ 2.6 GeV
2. **QCD binding**: The constituent quark mass contribution ~ (300 MeV)^2

The SUSY model's meson is the confined composite Q^d Q-bar_s / Lambda, and its mass should include the QCD contribution, not just the soft SUSY-breaking term.

### Resolution: the QCD contribution is already there, but misidentified

The soft mass f_pi^2 = (92 MeV)^2 is NOT the full meson mass. It is the SUSY-breaking splitting between the scalar meson and its fermionic partner (the mesino). The full scalar meson mass should be:

m^2_scalar = m^2_QCD + f_pi^2

where m^2_QCD ~ (few hundred MeV)^2 is the QCD binding contribution.

However, in the Seiberg effective theory, the meson M^i_j is a fundamental field below the confinement scale. Its mass comes entirely from the superpotential and the soft terms. The "QCD binding" is already encoded in the superpotential through the determinant constraint and the quark mass terms. The soft mass f_pi^2 is the dominant mass for the off-diagonal scalar mesons.

### Can the off-diagonal masses be lifted?

Five mechanisms, assessed:

**Mechanism 1: Non-universal soft masses.**
If SUSY-breaking generates m_tilde^2(off-diag) >> f_pi^2 while keeping m_tilde^2(diag) ~ f_pi^2, the off-diagonal scalars can be pushed above M_FCNC. Required enhancement: m_tilde^2(off) > (745 MeV)^2 / 2 ~ 278,000 MeV^2, a factor of 33 above f_pi^2.

**Assessment:** Possible in gravity mediation or anomaly mediation, where the soft mass depends on the field's quantum numbers. The off-diagonal mesons carry flavor quantum numbers (strangeness, etc.) that the diagonal mesons do not, so a flavor-dependent mediation could split them. But this introduces free parameters.

**Classification:** (iii) assumption. No mechanism is specified.

**Mechanism 2: D-term from gauged flavor symmetry.**
If a subgroup of SU(3)_flavor is gauged (e.g., the Cartan U(1)^2), D-terms generate mass-squared contributions proportional to the flavor charge:

m^2_D = g_F^2 D_a (T^a)_{ij}

For the off-diagonal meson M^d_s with flavor charge Delta S = 1:

m^2_D ~ g_F^2 <D_S> ~ g_F^2 * f_pi^2

This doubles the mass (for g_F ~ 1), which is insufficient.

**Assessment:** D-terms from gauged flavor provide at most O(1) corrections. Not sufficient.

**Classification:** (i) straightforward to compute; result is negative.

**Mechanism 3: Off-diagonal mesons decouple at Lambda.**
If the off-diagonal meson composites are formed above the confinement scale Lambda = 300 MeV and their mass includes the confinement-scale contribution:

m_off-diag ~ Lambda = 300 MeV

This is above the CKM-suppressed bound (167 MeV) but below the unit-coupling bound (745 MeV). Marginal.

**Assessment:** This is plausible. The Seiberg effective theory is valid below Lambda, and the off-diagonal meson masses should include the confinement-scale physics. But the precise calculation is non-perturbative.

**Classification:** (ii) non-perturbative.

**Mechanism 4: AMSB contribution.**
Anomaly-mediated SUSY breaking generates scalar masses:

m^2_AMSB = -(d gamma/d ln mu) * (m_{3/2} / (16 pi^2))^2

For the meson anomalous dimension gamma ~ 0 (at the conformal point of N_f = N_c SQCD), d gamma / d ln mu is nonzero only through running. The AMSB contribution depends on the gravitino mass m_{3/2}. For m_{3/2} ~ 50 TeV:

m^2_AMSB ~ (50000 GeV / (16 pi^2))^2 ~ (316 GeV)^2

This overwhelms f_pi^2 by 6 orders of magnitude and easily satisfies FCNC bounds. But m_{3/2} ~ 50 TeV is an external input.

**Assessment:** AMSB naturally lifts all scalar masses to the TeV scale, solving the FCNC problem. But it requires specifying the SUSY-breaking sector (gravitino mass).

**Classification:** (i) straightforward once m_{3/2} is specified; introduces free parameter.

**Mechanism 5: The effective coupling is CKM-suppressed, not O(1).**
The meson-quark coupling in the effective theory is related to the quark mixing matrix. If the effective Delta S = 2 coupling is:

g_eff ~ V_ds * V_ds = (0.225)^2 = 0.051

then the Wilson coefficient scales as g_eff^2 / m_S^2 = 0.0026 / (92)^2 = 3.1 x 10^{-7} MeV^{-2} = 3.1 x 10^{-13} GeV^{-2}, compared to the UTfit bound of 9 x 10^{-13} GeV^{-2}.

**Ratio = 0.34.** The constraint is SATISFIED.

**Assessment:** If the off-diagonal meson couplings are CKM-suppressed (g ~ V_ds^2 for the Delta S = 2 operator), the FCNC bound is satisfied at the seesaw vacuum with f_pi = 92 MeV. This is the most economical resolution.

**But is the coupling CKM-suppressed?** The off-diagonal meson M^d_s is a d s-bar composite. Its coupling to external d and s quarks is not the CKM matrix element V_ds; it is the overlap between the confined composite and the quark state. In the Seiberg effective theory, this coupling is O(1) by construction (the meson IS the quark bilinear, normalized by Lambda). The CKM suppression applies only in the quark mass eigenstate basis after rotating from the weak eigenstate basis.

At the diagonal seesaw vacuum, the mass eigenstates and weak eigenstates coincide (the Yukawa matrix is diagonal). The coupling of M^d_s to physical d and s quarks is O(1), not CKM-suppressed. CKM suppression would apply to the coupling of M^d_s to physical quarks in a DIFFERENT basis, but the meson IS in the mass eigenstate basis.

**Classification:** The CKM-suppression argument fails for the diagonal vacuum. (i) straightforward to verify; result is that the coupling is O(1).

### Summary of the FCNC problem

The factor of 8 deficit (g = 1) is a genuine problem. The most promising resolutions are:

1. **AMSB** (pushes all scalar masses to TeV, but requires specifying m_{3/2}): completely solves the problem but introduces a parameter.

2. **Off-diagonal mesons at Lambda** (m ~ 300 MeV): marginal, below the g = 1 bound but above the g = V_ds bound. A non-perturbative calculation is needed.

3. **Non-universal soft masses** (ad hoc but possible in specific mediation schemes): solves the problem but is an assumption.

**Is it a showstopper?** No. The factor of 8 is within the range of non-perturbative uncertainties in the meson mass (which could be anywhere from f_pi to Lambda). It is a problem for the SPECIFIC choice m_tilde = f_pi, but not for the framework. The honest statement is: "the off-diagonal meson scalar masses are not determined by the tree-level Lagrangian alone; non-perturbative QCD corrections and/or the full SUSY-breaking mediation are needed to fix them, and the FCNC constraint m_S > 745 MeV (g=1) or m_S > 167 MeV (g = V_ds) is a benchmark that the complete theory must satisfy."

---

## (c) The Paper's Honest Parameter Count

### Free parameters (irreducible inputs)

| # | Parameter | Value | What it determines | Can it be eliminated? |
|---|-----------|-------|--------------------|-----------------------|
| 1 | **m_s** | 93.4 MeV | Overall mass scale for Koide chain | No. Fundamental input. |
| 2 | **m_u** | 2.16 MeV | Up quark mass; no Koide constraint | No. Unrelated to other masses. |
| 3 | **m_d** | 4.67 MeV | Down quark mass; connected to theta_C via Oakes | If theta_C is input, m_d is determined. If both are free, count one. |
| 4 | **m_t** | 172.76 GeV | Top quark mass; external (not confined) | No. No Koide relation connects it to light quarks without m_c, m_b as inputs. |
| 5 | **Lambda** | 300 MeV | Confinement scale; sets C = Lambda^2 (m_u m_d m_s)^{1/3} | Possibly related to f_pi through Lambda_QCD. If so, not independent. |

### Determined parameters (follow from the free parameters plus structural constraints)

| Parameter | Determined by | Relation | Predicted value | PDG | Precision |
|-----------|--------------|----------|-----------------|-----|-----------|
| m_c | Koide seed (gv/m = sqrt(3)) | m_c = m_s * (2 + sqrt(3))^2 / (2 - sqrt(3))^2 = 13.93 m_s | 1301 MeV | 1270 +/- 20 | 1.5 sigma (2.4%) |
| m_b | v_0-doubling (bion minimum) | sqrt(m_b) = 3 sqrt(m_s) + sqrt(m_c) | 4177 MeV | 4180 +/- 30 | 0.1 sigma (0.07%) |
| lambda | Higgs mass at tan beta = 1 | lambda = sqrt(2) m_h / v | 0.72 | m_h = 125.25 | 0.9 sigma (0.12%) |
| f_pi | QCD observable | Not a free parameter of the SUSY model | 92 MeV | 92.1 +/- 0.8 | -- |
| v_EW | Fermi constant | v = (sqrt(2) G_F)^{-1/2} | 246.22 GeV | -- | -- |
| g, g' | Gauge couplings | Measured, not predicted | -- | -- | -- |

### Assumed parameters (fitted, not derived)

| Parameter | Value | Role | What would derive it |
|-----------|-------|------|---------------------|
| c = -1/12 (Kahler) | Fitted | Places Kahler pole at seesaw VEV | One-loop ISS Kahler computation |
| tan beta = 1 | Assumed | Needed for Koide at the quark level | Symmetry principle or dynamical EWSB |
| c_3 (instanton) | O(1) | Controls bloom dynamics | Semiclassical instanton calculation |
| s_k signs | (-1,+1,+1) | Bloom sign pattern | Monopole zero-mode analysis |

### Undetermined parameters (needed but not specified)

| Parameter | Role | Why undetermined |
|-----------|------|-----------------|
| b (soft bilinear) | m_A (pseudoscalar Higgs) | Higgs sector beyond tree-level |
| M_1, M_2 (gaugino masses) | Chargino/neutralino spectrum | SUSY-breaking mediation not specified |
| A_c, A_b, A_t (A-terms) | Scalar mass splitting | SUSY-breaking mediation not specified |
| m_{H_u}^2, m_{H_d}^2 | EWSB conditions | Constrained by tan beta = 1 but not derived |
| mu-term | Higgsino mass | With X != S, the mu-term is a free NMSSM parameter |

### Comparison to SM quark sector

**SM quark sector**: 6 masses + 4 CKM parameters = **10 free parameters**.

**This model's quark sector**:
- Free: m_u, m_d, m_s, m_t = **4 masses** (reduced from 6 by Koide + v_0-doubling)
- CKM: theta_C connected to m_d/m_s (Oakes); theta_23, theta_13, delta_CP NOT determined = **3 CKM parameters** free
- Total: **7 free parameters** (compared to SM's 10)

**Reduction**: 3 parameters (m_c, m_b, and one CKM angle related to m_d/m_s).

**But the model also introduces**:
- Lambda (1 parameter, unless related to Lambda_QCD)
- c = -1/12 (1 fitted parameter)
- tan beta = 1 (1 assumption, could be a constraint)

So the honest count is: **7 quark parameters + Lambda + c = 9**, compared to SM's 10, for a net reduction of **1 parameter** in the quark sector. If Lambda is identified with Lambda_QCD (not independent), the reduction is **2 parameters**.

The real gain is not in the parameter count but in the structural constraints: m_c and m_b are not arbitrary but are algebraic functions of m_s (and each other). This reduces the parameter space from a 10-dimensional continuum to a 7-dimensional submanifold, with specific testable predictions at 2.4% and 0.07% accuracy.

### Lepton sector comparison

**SM lepton sector**: 3 masses = **3 free parameters** (ignoring neutrinos).

**This model's lepton sector**: If the Sp(2) proposal works, the three lepton masses are determined by Lambda_L (1 parameter) + delta (dynamically determined) = **1 free parameter**. This would be a reduction from 3 to 1. But the Sp(2) proposal is not yet established.

### Grand total (quark + lepton)

| | SM | This model | Reduction |
|---|---|-----------|-----------|
| Quark masses | 6 | 4 | 2 |
| CKM | 4 | 3 (or 2.5 if Oakes is a prediction) | 1 (or 1.5) |
| Lepton masses | 3 | 1-3 (depending on Sp(2)) | 0-2 |
| Additional | 0 | Lambda, c, tan beta (1-3) | (-1 to -3) |
| **Total** | **13** | **9-13** | **0-4** |

The honest answer: the parameter reduction is modest (0 to 4, depending on assumptions). The value of the model is not primarily in parameter counting but in the structural relations it imposes.

---

## (d) The Volkov-Akulov Connection

### Background

Volkov and Akulov (1973) showed that a massless goldstino (the Goldstone fermion of spontaneously broken SUSY) has universal interactions determined by a single scale:

L_VA = -f^2 det(delta_mu^nu + (i/f^2) psi-bar sigma-bar^nu partial_mu psi)

where f^2 has dimensions [mass]^4 and psi is the goldstino (Weyl spinor). The scale sqrt(f) is the SUSY-breaking scale: f = <F_X> where F_X is the F-term that breaks SUSY.

Volkov and Akulov originally proposed that the goldstino IS the neutrino. In that case, the goldstino's couplings to matter are fixed by f, and neutrino phenomenology constrains f.

### Does this constrain the O'Raifeartaigh parameters?

The O'Raifeartaigh model has superpotential:

W = f Phi_0 + m Phi_a Phi_b + g Phi_0 Phi_a^2

At the SUSY-breaking vacuum (phi_b = 0, phi_0 = v, phi_a = 0):

F_{Phi_0} = f (the SUSY-breaking F-term)
F_{Phi_a} = 0
F_{Phi_b} = 0

The goldstino is the fermionic component of Phi_0. Its coupling scale is:

f^2_VA = |F_{Phi_0}|^2 = f^2

So the VA scale is sqrt(f) = f^{1/2} (the square root of the superpotential parameter f, not to be confused with sqrt(|F|^2)).

### If the goldstino is the neutrino

The VA proposal identifies the goldstino with the electron neutrino. The coupling to matter has the form:

L_int ~ (1/f^2) (psi_nu partial psi_nu)(psi_e partial psi_e)

which is a four-fermion interaction suppressed by f^2 = |F|^2. For this to reproduce the Fermi constant:

G_F / sqrt(2) ~ 1 / f^2

This gives:

f^2 = 1 / (sqrt(2) G_F) = v^2 / 2 = (174 GeV)^2

|F|^{1/2} = sqrt(f) = sqrt(174 GeV) ... no, let me be careful with dimensions.

f has dimension [mass]^2 (it is the F-term VEV). The four-fermion coupling is:

L ~ (1/|F|^2) * (J^mu J_mu)

where J^mu ~ psi-bar sigma-bar^mu psi is a current. Comparing to the Fermi theory:

1/|F|^2 ~ G_F / sqrt(2) = 1 / (v^2)

|F| = v / sqrt[4]{2} = 246 / 1.189 = 207 GeV ... that is |F|, not sqrt(|F|).

Wait, the VA Lagrangian has the form:

L ~ -(1/F^2) * partial psi partial psi psi psi

where F^2 = f^2 has dimension [mass]^4. Matching to Fermi:

1/F^2 ~ G_F ~ 1/(v^2)

F^2 ~ v^2 => |F| ~ v ~ 246 GeV

But |F| is an F-term with dimension [mass]^2, so:

|F| ~ v ~ 246 GeV is dimensionally wrong. F has [mass]^2 and v has [mass].

The correct dimensional analysis: G_F has dimension [mass]^{-2}. The VA coupling is 1/F^2 with [F^2] = [mass]^4. So:

1/F^2 ~ G_F ~ 1.166 x 10^{-5} GeV^{-2}

F^2 ~ 1/G_F ~ 8.6 x 10^4 GeV^2

|F| ~ 293 GeV^2... no, this is not right either.

Let me use the standard result. The Volkov-Akulov action is:

L_VA = -f^2 sqrt{-det(g_mu_nu + (i/f^2) partial_mu nu-bar sigma_nu nu + h.c.)}

where f^2 is the goldstino decay constant with dimension [mass]^2. The neutrino scattering cross section from the VA interaction is:

sigma ~ E^2 / f^4

where E is the CM energy. At E ~ 1 MeV (nuclear beta decay scale):

sigma ~ (1 MeV)^2 / f^4 ~ G_F^2 E^2 ~ (1.17 x 10^{-5} GeV^{-2})^2 (1 MeV)^2 = 1.37 x 10^{-22} MeV^{-2}

So f^4 ~ (1 MeV)^2 / sigma ~ (1 MeV)^2 / (1.37 x 10^{-22} MeV^{-2}) = 7.3 x 10^{21} MeV^4

f^2 ~ 2.7 x 10^{10.5} MeV^2 ~ 8.5 x 10^{10} MeV^2

f ~ sqrt(8.5 x 10^{10}) MeV ~ 2.9 x 10^5 MeV = 290 GeV

### Result: sqrt(f) ~ 290 GeV (the EW scale)

This is the well-known result: if SUSY is broken at sqrt(F) ~ 290 GeV, the goldstino's interactions have the strength of the weak interaction. This is the reason VA originally proposed the identification.

### Constraint on O'Raifeartaigh parameters

In the O'Raifeartaigh model: f (the linear term in W) plays the role of the SUSY-breaking parameter. The VA constraint gives:

f_OR ~ (290 GeV)^2 = 8.4 x 10^4 GeV^2 = 8.4 x 10^{10} MeV^2

The O'Raifeartaigh model also has the Koide condition gv/m = sqrt(3), and the mass scale m sets the seed masses:

m_{seed} = (2 +/- sqrt(3)) m

If the seed corresponds to (0, m_s, m_c) = (0, 93.4, 1270) MeV:

m_+ = (2 + sqrt(3)) m = m_c => m = m_c / (2 + sqrt(3)) = 1270 / 3.732 = 340 MeV
m_- = (2 - sqrt(3)) m = m_s => m = m_s / (2 - sqrt(3)) = 93.4 / 0.268 = 349 MeV

The small discrepancy (340 vs 349 MeV) reflects the 2.4% deviation of the predicted m_c from the observed value. Taking m ~ 345 MeV.

Then: gv = sqrt(3) m = sqrt(3) * 345 = 597 MeV.

And: f_OR = 8.4 x 10^{10} MeV^2.

The ratio f_OR / m^2 = 8.4 x 10^{10} / (345)^2 = 8.4 x 10^{10} / 1.19 x 10^5 = 7.1 x 10^5.

This is a large hierarchy: the SUSY-breaking parameter f is 7 x 10^5 times larger than the squared mass scale m^2 of the O'Raifeartaigh model. In units where m = 345 MeV:

f_OR = 7.1 x 10^5 * m^2

This is a FINE-TUNING: the O'Raifeartaigh model with f >> m^2 has a strong hierarchy between SUSY-breaking and the fermion mass scale. The pseudo-modulus VEV v is determined by the Kahler stabilization, not by f.

### Does the goldstino = neutrino identification work?

**Problems:**
1. The neutrino has mass (m_nu ~ 0.01-0.1 eV). The goldstino is massless in exact SUSY breaking. Once SUSY breaking is mediated to the gravitational sector, the goldstino is eaten by the gravitino (super-Higgs mechanism). The physical neutrino is NOT the goldstino.

2. The goldstino couples to ALL matter with strength ~ 1/F^2. The neutrino couples only to weak-interaction matter. The VA coupling is universal (independent of flavor), while the weak interaction distinguishes flavors.

3. LEP measured the invisible Z width as N_nu = 2.984 +/- 0.008, exactly consistent with 3 SM neutrinos. An additional goldstino coupling to the Z would increase this.

**Assessment:** The VA identification of goldstino = neutrino was an interesting historical idea (1973) but is ruled out by modern precision data. The goldstino is eaten by the gravitino at scale m_{3/2} = F/M_Pl ~ (290 GeV)^2 / (2.4 x 10^{18} GeV) ~ 35 meV. This is in the right ballpark for neutrino mass (!) but the coupling structure is wrong.

### What the VA connection DOES constrain

Even without the goldstino = neutrino identification, the VA scale constrains the O'Raifeartaigh model:

- **sqrt(F) ~ 290 GeV** implies the SUSY-breaking sector is at the EW scale. This is consistent with the model's identification of the confined sector with the Higgs mechanism.

- **f_OR ~ 8.4 x 10^{10} MeV^2** with m ~ 345 MeV and gv ~ 597 MeV gives the hierarchy f/m^2 ~ 7 x 10^5. This means the linear term f dominates the potential, which is the standard O'Raifeartaigh regime (SUSY is broken by F_{Phi_0} = f, not by the quadratic/cubic terms).

- **The gravitino mass**: m_{3/2} = F / (sqrt(3) M_Pl) = 8.4 x 10^{10} MeV^2 / (sqrt(3) * 2.44 x 10^{21} MeV) = 2.0 x 10^{-11} MeV = 20 peV. This is 3 orders of magnitude below the lightest neutrino mass. If f is the ONLY source of SUSY breaking (no hidden sector), the gravitino is ultralight.

- **AMSB contribution**: m_AMSB ~ (g^2 / 16 pi^2) m_{3/2} ~ (0.01 / 160) * 20 peV ~ 1.2 feV. Negligible. AMSB with the VA-scale SUSY breaking cannot lift the meson masses.

### Conclusion on VA

The Volkov-Akulov connection provides the constraint sqrt(F) ~ 290 GeV, identifying the SUSY-breaking scale with the electroweak scale. This is consistent with the model but does not uniquely determine the O'Raifeartaigh parameters (f, m, g are constrained by two conditions: f ~ (290 GeV)^2 and gv/m = sqrt(3), leaving v or equivalently g as a free parameter). The goldstino = neutrino identification is historically interesting but experimentally ruled out. The resulting gravitino mass m_{3/2} ~ 20 peV is consistent with cosmological bounds but too small for AMSB to be relevant.

---

## (e) Computations for Round 18, Ranked

### RANK 1: Explicit Fritzsch Texture with Koide Constraints

**Task:** Write the 3x3 quark mass matrices M_u (up-type: u, c, t) and M_d (down-type: d, s, b) with Fritzsch texture:

    M_u = [[0, A_u, 0], [A_u*, 0, B_u], [0, B_u*, C_u]]
    M_d = [[0, A_d, 0], [A_d*, 0, B_d], [0, B_d*, C_d]]

The eigenvalues of M_u are m_u, m_c, m_t; the eigenvalues of M_d are m_d, m_s, m_b. The Fritzsch texture gives:

    |A_u|^2 = m_u * m_c,   |B_u|^2 = m_c * m_t,   C_u ~ m_t
    |A_d|^2 = m_d * m_s,   |B_d|^2 = m_s * m_b,   C_d ~ m_b

The CKM matrix is V = U_L^{u dag} U_L^d. Compute all four CKM parameters from the Koide-constrained quark masses:

    m_c = 1301 MeV (Koide seed), m_b = 4177 MeV (v_0-doubling)

with inputs m_u = 2.16 MeV, m_d = 4.67 MeV, m_s = 93.4 MeV, m_t = 172760 MeV.

Check: sin theta_12 = |sqrt(m_d/m_s) - sqrt(m_u/m_c) e^{i delta}| ~ sqrt(m_d/m_s) = 0.224 (Oakes).

Compute sin theta_23 and sin theta_13 from the Fritzsch texture and compare to PDG.

**Why highest priority:** This is the North Star deliverable. It connects the Koide-determined masses to the CKM matrix via a specific, testable texture. If the CKM angles come out within 20% of their PDG values, the framework is a concrete model. If not, the Fritzsch texture is ruled out and an alternative must be found.

**Difficulty:** Straightforward (eigenvalue problem for 3x3 matrices).

### RANK 2: Full Off-Diagonal Meson Potential Without X = S

**Task:** Recompute the scalar potential with the corrected model (X is Lagrange multiplier with F_X = 0; separate NMSSM singlet S). At the seesaw vacuum:

1. The dominant F-terms are now F_{M^d_d} = y_c v/sqrt(2) and F_{M^s_s} = y_b v/sqrt(2). The F_X contribution (which drove the tachyonic instability in vacuum_stabilization.md) is REMOVED because F_X = 0 at the seesaw vacuum (before NMSSM coupling).

2. Compute the 18-real-variable Hessian of V(M, X) at the diagonal seesaw vacuum with ONLY the Seiberg superpotential + Yukawa terms + soft mass f_pi^2.

3. Determine whether any off-diagonal meson directions are still tachyonic.

4. If the vacuum is now stable (all eigenvalues positive), the FCNC constraint is automatically satisfied (off-diagonal mesons at mass ~ f_pi = 92 MeV are the physical scalar spectrum). If tachyons remain, identify their origin.

**Why second priority:** The vacuum stability analysis in vacuum_stabilization.md was contaminated by the X = S identification, which generated the dominant F_X term. With X != S, the off-diagonal B-terms M_c * F_X that drove the tachyons are ABSENT. The vacuum may now be stable, which would change the entire picture: no tachyonic CKM, no off-diagonal condensate, and the FCNC problem is mitigated (not solved, but the tachyonic instability is the main driver).

**Difficulty:** Straightforward (numerical Hessian of a known potential with different parameter values).

### RANK 3: One-Loop Kahler Potential for X in the ISS Dual

**Task:** In the ISS framework (SU(N_c) SQCD with N_f > N_c, at the metastable SUSY-breaking vacuum), compute the one-loop effective Kahler potential for the pseudo-modulus X:

    K_eff(X, X-bar) = K_tree + delta K_1-loop

The one-loop correction is:

    delta K = (1/32 pi^2) STr[M^2(X) (ln(M^2(X) / mu^2) - 1)]

where M^2(X) is the X-dependent mass matrix of the ISS fields. Extract the coefficient of the |X|^4 term and check whether it is -1/12 (as assumed in the model).

The ISS framework for N_f = N_c = 3 is at the boundary of the ISS window (N_f = N_c is the degenerate case where the magnetic dual has gauge group SU(0), i.e., it is an effective Wess-Zumino model). The pseudo-modulus stabilization has been studied by Intriligator, Seiberg, and Shih (ISS, 2006-2007) for general N_f > N_c. Adapt their results to N_f = N_c = 3.

**Why third priority:** The c = -1/12 coefficient is the least well-motivated parameter in the model. A microscopic derivation would put the O'Raifeartaigh-Koide mechanism on firm ground. If the one-loop Kahler gives c != -1/12, the Koide seed mechanism needs revision.

**Difficulty:** Moderate. The ISS one-loop Kahler is a standard Coleman-Weinberg calculation, but the N_f = N_c limit is degenerate and may require careful handling.

### RANK 4: AMSB Spectrum for the Meson Sector

**Task:** Assuming anomaly-mediated SUSY breaking with gravitino mass m_{3/2} (as a free parameter), compute:

1. The AMSB scalar mass-squared for each meson M^i_j:

    m_AMSB^2 = -(1/4) (d gamma_M / d ln mu) m_{3/2}^2

where gamma_M is the anomalous dimension of the meson composite.

2. The A-term: A_AMSB = -(beta_y / y) m_{3/2}

3. The gaugino masses: M_a = (beta_a / g_a) m_{3/2}

Determine: (a) the value of m_{3/2} that lifts the off-diagonal meson scalars above the FCNC bound (m_S > 745 MeV), and (b) whether the AMSB spectrum is compatible with the Koide conditions (i.e., whether the AMSB soft terms preserve the Koide-enforced mass ratios).

**Why fourth priority:** AMSB is the most natural mediation scheme for a model where the SUSY-breaking scale is in the meson sector. It predicts a specific pattern of soft masses that could solve the FCNC problem and provide the missing gaugino masses. Computing the AMSB spectrum would make the model phenomenologically complete.

**Difficulty:** Straightforward once gamma_M is known. The anomalous dimension of the meson composite in confined N_f = N_c SQCD is gamma_M = (N_f - N_c)/N_f = 0 at the Seiberg point, but the derivative d gamma/d ln mu is nonzero and computable from the NSVZ beta function.

### RANK 5: Sp(2) O'Raifeartaigh-Koide Model for Leptons

**Task:** Construct the Sp(2) SQCD model with N_f = 3 explicitly:

1. Write the ADS superpotential for Sp(2), N_f = 3.
2. Add mass deformation and verify that the vacuum has O'Raifeartaigh structure.
3. Check whether gv/mu = sqrt(3) is satisfied at the minimum.
4. Compute the lepton mass spectrum and compare to (m_e, m_mu, m_tau).
5. Identify the mediation mechanism to SM leptons.

**Why fifth priority:** The lepton sector is the model's biggest gap, but the Sp(2) proposal introduces an entirely new sector with its own confinement dynamics. This is a large computation that may not converge in a single round. The quark-sector tasks (Ranks 1-4) are more immediately relevant to completing the paper.

**Difficulty:** High. Requires non-perturbative SUSY dynamics for Sp(2) gauge theories, which are less well-studied than SU(N).

---

## Summary Table

| Task | What it resolves | Priority | Difficulty |
|------|-----------------|----------|------------|
| Fritzsch texture + CKM | CKM angles from Koide masses | RANK 1 | Straightforward |
| Off-diagonal vacuum (no X=S) | Vacuum stability, FCNC reanalysis | RANK 2 | Straightforward |
| ISS one-loop Kahler | c = -1/12 derivation | RANK 3 | Moderate |
| AMSB spectrum | FCNC, gaugino masses, soft terms | RANK 4 | Straightforward |
| Sp(2) lepton model | Lepton masses | RANK 5 | High |

**Recommended Round 18 allocation:**
- Launch Ranks 1 and 2 in parallel (both straightforward, independent).
- If time permits, start Rank 3 as a literature search + calculation.
- Defer Ranks 4 and 5 to Round 19.

---

*Generated: 2026-03-04*
*Based on: complete_lagrangian_opus.md, brainstorm_r16.md, fcnc_constraints.md, full_spectrum.md, seesaw_mechanism.md, oraifeartaigh_koide.md, bloom_mechanism.md, mesino_phenomenology.md, lepton_alternatives.md, ckm_koide.md, v0_doubling.md, vacuum_stabilization.md, round2_status.md*
