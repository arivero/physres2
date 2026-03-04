# Coleman-Weinberg Corrections to Off-Diagonal Meson Masses

## Setup

N=1 SU(3) SQCD with N_f = N_c = 3 in the Seiberg confined phase. The low-energy effective theory contains mesons M^i_j (3x3 complex matrix), baryons B, B-tilde, and a Lagrange multiplier X enforcing the quantum-modified moduli constraint.

Superpotential:

    W = m_i M^i_i + X(det M - B B-tilde - Lambda^6)

Soft SUSY-breaking potential:

    V_soft = f_pi^2 Tr(M^dag M)

The seesaw vacuum (B = B-tilde = 0, M diagonal):

    <M^i_i> = C / m_i,    C = Lambda^2 (m_u m_d m_s)^{1/3}

    X_0 = -C / Lambda^6

| Parameter | Value |
|-----------|-------|
| m_u | 2.16 MeV |
| m_d | 4.67 MeV |
| m_s | 93.4 MeV |
| Lambda | 300 MeV |
| f_pi | 92 MeV |
| C | 882297 MeV^2 |
| M_u = C/m_u | 408471 MeV |
| M_d = C/m_d | 188929 MeV |
| M_s = C/m_s | 9446 MeV |
| X_0 | -1.210 x 10^{-9} MeV^{-4} |


## Part (a): Field-Dependent Mass Matrices

The six off-diagonal meson fields are phi_{ab} = M^a_b for a != b. The scalar potential V = sum |F_I|^2 + f_pi^2 Tr(M^dag M) gives a 6x6 mass-squared matrix at the seesaw vacuum:

    m^2_{(ab),(cd)} = 2 sum_{i,j} [dF_{M^i_j}/dM^a_b] [dF_{M^i_j}/dM^c_d]
                    + 2 [dF_X/dM^a_b][dF_X/dM^c_d]
                    + 2 f_pi^2 delta_{(ab),(cd)}

At the diagonal seesaw vacuum, all F-terms vanish. The derivatives of F-terms with respect to off-diagonal fields are:

    dF_{M^i_j}/dM^a_b = X_0 * d(cof(i,j))/dM^a_b

The nonzero cofactor derivatives at diagonal M are:

    d cof(a,b) / d M^b_a = -M_c    (c is the spectator: c = 3 - a - b)

All other cofactor derivatives with respect to off-diagonal fields vanish at the diagonal vacuum. Crucially:

    dF_X/dM^a_b = cof(a,b)|_diag = 0    for a != b

The F_X constraint does not contribute to the off-diagonal mass matrix because the off-diagonal cofactors vanish at a diagonal matrix.

The resulting 6x6 mass matrix is **exactly block-diagonal**, decomposing into three 2x2 blocks for conjugate pairs (M^a_b, M^b_a):

    m^2_{(ab),(ab)} = 2 f_pi^2 + 2 X_0^2 M_c^2
    m^2_{(ab),(ba)} = 0

The off-diagonal entry within each 2x2 block vanishes because M^a_b and M^b_a couple to *different* F-terms: dF_{M^b_a}/dM^a_b is nonzero, but dF_{M^b_a}/dM^b_a is zero (the variable M^b_a does not appear in its own cofactor minor). This prevents the (M^a_b, M^b_a) mass matrix from having off-diagonal entries.

| Pair | Spectator M_c | m^2 [MeV^2] | F-term/soft ratio |
|------|---------------|-------------|-------------------|
| u-d | M_s = 9446 | 16928 + 2.6e-10 | 1.5e-14 |
| u-s | M_d = 188929 | 16928 + 1.0e-7 | 6.2e-12 |
| d-s | M_u = 408471 | 16928 + 4.9e-7 | 2.9e-11 |

The F-term contribution X_0^2 M_c^2 is negligible (ratio 10^{-11} to 10^{-14}) because |X_0| = 1.21 x 10^{-9} MeV^{-4} is extremely small. All six eigenvalues are degenerate at 2 f_pi^2 = 16928 MeV^2.

**Cross-pair couplings** (between the u-d, u-s, and d-s sectors) are **exactly zero** at quadratic order. This follows from the structure of det M: the cubic terms M^u_d M^d_s M^s_u + M^u_s M^s_d M^d_u are trilinear in off-diagonal fields and do not contribute to the quadratic mass matrix.


## Part (b): One-Loop Effective Potential

### Structure of the CW correction

In exact SUSY, V_CW = 0 by the supertrace cancellation. The soft breaking V_soft = f_pi^2 Tr(M^dag M) shifts boson masses by f_pi^2 relative to fermion masses. For each off-diagonal meson field:

    m^2_boson = X_0^2 M_c^2 + f_pi^2
    m^2_fermion = X_0^2 M_c^2

The SUSY fermion mass for the off-diagonal mesons is m_F = |X_0 M_c|, which is extremely small:

| Pair | m_F [MeV] | m_F^2 [MeV^2] |
|------|-----------|---------------|
| u-d | 1.14e-5 | 1.31e-10 |
| u-s | 2.29e-4 | 5.23e-8 |
| d-s | 4.94e-4 | 2.44e-7 |

These tiny fermion masses arise because the only coupling of off-diagonal mesons to other fields at the seesaw vacuum goes through X_0, which has dimensions MeV^{-4} and magnitude 10^{-9}.

### Leading contribution

The CW correction to the off-diagonal meson mass-squared is proportional to:

    delta m^2_CW ~ (g_eff^2 / 16 pi^2) * f_pi^2 * ln(Lambda^2/m^2)

where g_eff = |X_0 M_c| is the effective coupling strength. The key factors:

1. **Proportional to**: g_eff^2 * f_pi^2, where g_eff = |X_0 M_c| ~ 10^{-5} to 10^{-4}
2. **Sign**: POSITIVE (stabilizing). The soft breaking increases the boson mass relative to the fermion mass, and the CW potential pushes scalars back toward the SUSY vacuum.
3. **Numerical magnitude**:

| Pair | g_eff | delta m^2_CW [MeV^2] | Ratio to tree |
|------|-------|---------------------|---------------|
| u-d | 1.14e-5 | 1.7e-8 | 1.0e-12 |
| u-s | 2.29e-4 | 6.6e-6 | 3.9e-10 |
| d-s | 4.94e-4 | 3.1e-5 | 1.8e-9 |

The CW correction is suppressed by a factor of 10^{-9} to 10^{-12} relative to the tree-level mass 2 f_pi^2 = 16928 MeV^2. This double suppression comes from:
- Loop factor 1/(16 pi^2) ~ 6.3e-3
- Tiny coupling g_eff^2 ~ 10^{-8} to 10^{-10}

### Cross-pair mixing at one loop

The off-diagonal mesons M^a_b do NOT couple to X at the seesaw vacuum (because cof(a,b)|_diag = 0 for a != b). Therefore, there is no one-loop diagram that connects different off-diagonal pairs through X loops. The cubic det M vertex (M^u_d M^d_s M^s_u) involves three off-diagonal fields and contributes only at two loops.

**The one-loop CW preserves the block-diagonal structure. Cross-pair mixing is zero at one-loop order.**


## Part (c): Yukawa Contributions

The NMSSM coupling W += lambda_0 X H_u . H_d adds a Higgs sector coupled to X. The key properties:

1. **X is a flavor singlet.** Its coupling to the Higgs carries no flavor indices.

2. **Off-diagonal mesons decouple from X.** At the diagonal seesaw vacuum, cof(a,b) = 0 for a != b, so there is no M^a_b--X vertex. The one-loop diagram M^a_b -> X -> (H_u, H_d loop) -> X -> M^c_d has a zero vertex at each end.

3. **delta m^2(M^a_b) from NMSSM Higgs = 0 (exact at one loop).**

Even with explicit flavor-dependent Yukawa couplings (e.g., y_c H_u M^d_d + y_b H_d M^s_s), the Yukawa terms enter only the *diagonal* meson F-terms. They modify the seesaw condition but preserve the diagonal structure of the vacuum. The off-diagonal meson F-terms F_{M^a_b} = X cof(a,b) are unaffected by diagonal Yukawa couplings.

At two loops, Yukawa-Higgs corrections could contribute through the chain M^a_b -> (det M) -> M^a_a, M^b_b -> (Yukawa) -> H loop. The estimate:

    delta m^2 ~ (y^2/(16 pi^2))^2 * v^2

For y_b ~ 0.034: delta m^2 ~ 1.6 MeV^2 (ratio to tree: 9.5e-5).

This is still too small to overcome the tree-level mass, and it is a diagonal correction (not a cross-pair mixing).


## Part (d): Can Radiative Corrections Generate CKM Mixing?

### Tachyonic instability

For delta m^2 = -(alpha/4pi) M^2 to overcome 2 f_pi^2 = 16928 MeV^2:

| Scale M | Required alpha |
|---------|---------------|
| Lambda = 300 MeV | 2.36 |
| 1 GeV | 0.21 |
| v_ew = 174 GeV | 7e-6 |

At the QCD scale, alpha > 1 is not perturbatively achievable. At the EW scale, alpha ~ 10^{-5} would suffice, but the Higgs coupling to off-diagonal mesons is exactly zero (flavor-blindness of the NMSSM vertex + vanishing cofactor).

**A tachyonic instability in the off-diagonal meson sector is not physically plausible.**

### Non-tachyonic mixing

The mixing angle between different flavor pairs is theta ~ m^2_{cross} / Delta m^2_{diag}. Two problems:

1. **m^2_{cross} = 0 at one loop.** The block-diagonal structure is preserved by CW.

2. **Delta m^2_{diag} ~ 0.** All six off-diagonal mesons are degenerate at 2 f_pi^2 with splittings ~ 10^{-7} MeV^2 from the F-term. This is a degenerate perturbation theory problem.

In degenerate perturbation theory, the mixing angle is determined by the DIRECTION of the perturbation that lifts the degeneracy, not by the ratio of cross-coupling to splitting. The first perturbation that distinguishes between pairs controls the result:

- The F-term X_0^2 M_c^2 splits the pairs by ~ 10^{-7} MeV^2 (hierarchy: d-s > u-s > u-d, following M_u > M_d > M_s).
- Any Kahler correction that is larger than 10^{-7} MeV^2 will dominate.

But regardless of which perturbation dominates, it only lifts the degeneracy WITHIN each conjugate pair block -- it does not mix BETWEEN pairs.

**Radiative corrections cannot generate CKM mixing in the Seiberg effective theory at the seesaw vacuum.**


## Part (e): Kahler Potential Corrections

### c_2 term: Tr(M^dag M M^dag M)/Lambda^6

At quadratic order in off-diagonal fields around the seesaw vacuum:

    delta_2 Tr(MMMM) = sum_{a!=b} |phi_{ab}|^2 * 2(M_a^2 + M_b^2)
                      + sum_{a!=b} phi_{ab} phi_{ba} * 2 M_a M_b

This is **also block-diagonal**, coupling only conjugate pairs (phi_{ab}, phi_{ba}). It does not mix different off-diagonal sectors.

The c_2 term shifts the diagonal mass entries in a flavor-dependent way:

| Pair | 2(M_a^2 + M_b^2) [MeV^2] | Relative shift |
|------|--------------------------|----------------|
| u-d | 4.05e11 | 5.56e-4 * c_2 |
| u-s | 3.34e11 | 4.58e-4 * c_2 |
| d-s | 7.16e10 | 9.82e-5 * c_2 |

This creates a mass splitting between pairs, but no mixing.

### c_3 term: |det M|^2/Lambda^10

The Hessian of |det M|^2 in the off-diagonal sector is also block-diagonal, with nonzero entries only for conjugate pairs (phi_{ab}, phi_{ba}). The coupling is:

    d^2|det M|^2 / d(M^a_b) d(M^b_a) = -2 * det(M) * M_c

This provides a negative (potentially tachyonic) contribution within each pair block, with magnitude proportional to M_c. For c_3 < 0, the (d,s) pair (spectator M_u, the largest) would become tachyonic first.

The c_3 term does NOT generate cross-pair mixing.

### Fundamental obstruction

At quadratic order in off-diagonal fields around a diagonal vacuum, ALL single-trace polynomial Kahler invariants of the form Tr(M^dag M ... M^dag M) preserve the block-diagonal structure. This is a consequence of index contraction: each off-diagonal entry M^a_b in a trace contracts with entries sharing index a or b. At a diagonal vacuum, the only contraction producing a nonzero result pairs M^a_b with M^b_a (the conjugate), not with M^a_c or M^c_b for c != a,b.

This was verified numerically for:
- Tr(M^dag M M^dag M)
- |det M|^2
- Tr(M^T M^T M M)  [non-Hermitian trace variant]

All give block-diagonal Hessians in the off-diagonal sector.

### What would work

To generate cross-pair mixing at quadratic order, one needs:

1. **Multi-trace Kahler terms**: e.g., Tr(M^dag M) Tr(M^dag M^dag M M). These are higher order in 1/Lambda and typically subdominant.

2. **Non-perturbative Kahler corrections**: instanton-generated terms in the Kahler potential can have non-polynomial dependence on M that breaks the trace structure.

3. **Going beyond the confined description**: including magnetic quarks (as in the ISS framework) introduces additional fields whose loops could generate flavor-changing effects.

### Required magnitude for Cabibbo mixing

For a tachyonic instability in the u-s sector from the c_2 term:

    |c_2| > Lambda^6 / (M_u^2 + M_s^2) ~ 4400

This is highly unnatural for a perturbative Kahler correction.

For non-tachyonic mixing from a hypothetical cross-pair coupling c_cross:

    theta_12 ~ c_cross * f(M_diag) / [Delta m^2 / f_pi^2]

The energy cost of turning on a Cabibbo-scale off-diagonal VEV is:

    |M^u_s| ~ sin(13 deg) * M_u ~ 92000 MeV
    delta V ~ f_pi^2 * |M^u_s|^2 ~ 7.1 x 10^{13} MeV^2

This is 14 orders of magnitude above the SUSY-breaking scale m_u^2 ~ 4.7 MeV^2, confirming that the diagonal seesaw vacuum is deeply stable against off-diagonal perturbations.


## Summary Table

| Quantity | Value | Assessment |
|----------|-------|------------|
| Tree-level m^2(offdiag) | 16928 MeV^2 | All positive, degenerate |
| F-term splitting | 10^{-7} MeV^2 | Negligible (10^{-11} of tree) |
| CW correction (one-loop) | 10^{-5} MeV^2 | Positive, negligible |
| Cross-pair coupling (one-loop) | 0 | Exact, block-diagonal preserved |
| NMSSM Higgs correction | 0 | Exact, flavor-blind |
| Yukawa two-loop correction | ~1 MeV^2 | Small, still diagonal |
| c_2 Kahler: tachyonic threshold | c_2 > 4400 | Unnatural |
| c_2 Kahler: cross-pair mixing | 0 | Fundamental obstruction |
| c_3 Kahler: cross-pair mixing | 0 | Fundamental obstruction |


## Conclusion

One-loop Coleman-Weinberg corrections cannot generate off-diagonal meson VEVs or CKM mixing in the Seiberg effective theory for N_f = N_c = 3 SQCD at the seesaw vacuum. Three independent mechanisms prevent this:

1. **Magnitude**: The CW correction is suppressed by (X_0 M_c)^2 / (16 pi^2) ~ 10^{-12} relative to the tree-level soft mass, which is positive definite.

2. **Block-diagonal structure**: The mass matrix for off-diagonal mesons decomposes into three independent 2x2 blocks for conjugate pairs (M^a_b, M^b_a). No single-trace polynomial Kahler invariant or one-loop CW diagram generates cross-pair couplings at quadratic order around the diagonal vacuum.

3. **Flavor blindness**: The NMSSM coupling lambda_0 X H_u H_d is flavor-blind, and the off-diagonal mesons decouple from X at the seesaw vacuum (vanishing cofactors), so Higgs loops contribute nothing.

CKM mixing in this framework requires a mechanism that breaks the fundamental block-diagonal obstruction: either non-perturbative Kahler corrections (instanton-generated), going beyond the confined-meson effective theory, or a qualitatively different origin for flavor mixing.

---

*Numerical verification: `results/cw_offdiag_meson.py`*
