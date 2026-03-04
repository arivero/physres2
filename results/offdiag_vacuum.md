# Off-Diagonal Meson Vacuum Analysis

## Setup

Superpotential for N_f = N_c = 3 SQCD with two Higgs doublets:

```
W = sum_i m_i M^i_i + X(det M - B B~ - Lambda^6) + y_c H_u^0 M^d_d + y_b H_d^0 M^s_s
```

with soft SUSY-breaking term V_soft = m~^2 Tr(M†M), where m~^2 = f_pi^2 = (92 MeV)^2 = 8464 MeV^2.

**Inputs**: m_u = 2.16, m_d = 4.67, m_s = 93.4 MeV, m_c = 1270, m_b = 4180 MeV, Lambda = 300 MeV, v = 246220 MeV, y_c = 2m_c/v = 1.032e-2, y_b = 2m_b/v = 3.395e-2.

The scalar potential is V = sum |F_I|^2 + m~^2 Tr(M†M) with 14 real degrees of freedom (9 meson entries + X + B + B~ + H_u^0 + H_d^0).


## Reference: Diagonal (Seiberg) Vacuum

The Seiberg seesaw gives M^i_j = delta^i_j C/m_i with C = Lambda^2 (m_u m_d m_s)^{1/3} = 882297 MeV^2:

| Field | Value |
|-------|-------|
| M_u   | 408471 MeV |
| M_d   | 188929 MeV |
| M_s   | 9446 MeV |
| X     | -1.210e-9 MeV^{-4} |
| det M / Lambda^6 | 1.000000000000 |

At this vacuum with Hu = Hd = 0, the scalar potential is V_0 = 1.715e15 MeV^2, dominated by:
- F-term |F_{Hu}|^2 + |F_{Hd}|^2 from the Yukawa couplings (F_{Hu} = y_c M_d, F_{Hd} = y_b M_s)
- The soft term m~^2(M_u^2 + M_d^2 + M_s^2) ~ 1.4e15 MeV^2, overwhelmingly from M_u^2

The pure Seiberg vacuum has a *huge* soft-term cost because M_u = C/m_u ~ 4e5 MeV is enormous.


## Diagonal Minimum (with soft terms)

Optimizing V over diagonal M, X, with H_u, H_d integrated out analytically:

```
H_u = -(m_d + X M_u M_s)/y_c
H_d = -(m_s + X M_u M_d)/y_b
```

The remaining V_diag depends on (M_u, M_d, M_s, X). Minimizing:

| Field | Seesaw | Soft-optimized | Shift |
|-------|--------|---------------|-------|
| M_u   | 408471 | 344532 | -15.6% |
| M_d   | 188929 | 195492 | +3.5% |
| M_s   | 9446   | 10824  | +14.6% |
| X     | -1.21e-9 | 1.53e-5 | flipped sign |

V_diag = 1.329e15 MeV^2 (8% below the pure seesaw V_0, but still large).

The soft term pulls M_u down (it dominates Tr M†M) while the determinant constraint det M = Lambda^6 and the F-term for M^u_u maintain a balance. The X VEV flips sign, indicating the soft terms substantially perturb the Seiberg vacuum.


## Full 14D Minimum (off-diagonal VEVs allowed)

Searched with 46 initial conditions (30 small perturbation + 15 large perturbation + 1 diagonal), using Powell + Nelder-Mead two-stage optimization with adaptive refinement. The best minimum was polished with three rounds of high-precision iteration.

**Result: V_best = 2.169e13 MeV^2** (factor ~60 below the diagonal minimum).

### Meson matrix at the minimum

```
          u              d              s
u    -5816.11       3108.36       7660.15
d    18742.29      30560.77      -6503.04
s       21.08      10506.64      10010.27
```

The diagonal VEVs are completely rearranged from the Seiberg seesaw:
- M^u_u = -5816 (vs seesaw 408471, a factor -70x change)
- M^d_d = 30561 (vs 188929, down by ~84%)
- M^s_s = 10010 (vs 9446, up by ~6%)

Off-diagonal entries are of the same order as (or larger than) the diagonals.

### Other field VEVs

| Field | Value | Comment |
|-------|-------|---------|
| X | -6.20e-5 | 5 orders larger than seesaw |
| B | 6.51e5 | nonzero |
| B~ | -1.12e9 | nonzero, large |
| H_u^0 | 2.69e8 | very large |
| H_d^0 | 6.41e6 | large |
| det M | -1.26e12 | NOT equal to Lambda^6 |

The baryons B, B~ acquire nonzero VEVs, which absorb part of the F_X constraint (det M - BB~ - Lambda^6 = 0 is approximately satisfied via a large BB~ ~ -7.3e14 that compensates the sign-flipped det M).

### Potential decomposition

| Contribution | Value (MeV^2) | Fraction |
|-------------|---------------|----------|
| sum F_M^2   | 7.80e12 | 36% |
| F_X^2       | 1.95e6 | negligible |
| F_B^2 + F_Bt^2 | 4.84e9 | 0.02% |
| F_Hu^2 + F_Hd^2 | 2.15e5 | negligible |
| m~^2 Tr(MM) | 1.39e13 | 64% |
| **Total V** | **2.17e13** | 100% |

The minimum is a balance between F-term cost (~36%) and soft term cost (~64%). Compared to the diagonal vacuum where the soft term alone was ~1.4e15, this minimum dramatically reduces Tr(M†M) by making all meson entries O(10^4) rather than having M_u ~ 4e5.


## CKM-Like Mixing Angles

Singular value decomposition M = U diag(sigma) V^T:

Singular values: 37212, 15842, 2132 MeV (vs seesaw: 408471, 188929, 9446).

The CKM-like matrix V_CKM = U^T V:

```
[ 0.799   -0.533    0.278]
[-0.093    0.348    0.933]
[-0.594   -0.771    0.228]
```

Extracted mixing angles (standard parametrization):

| Angle | This vacuum | PDG value |
|-------|------------|-----------|
| theta_12 | 33.7 deg | 13.04 deg (Cabibbo) |
| theta_23 | 76.3 deg | 2.38 deg |
| theta_13 | 16.2 deg | 0.20 deg |

The mixing angles are very large -- the meson matrix is far from diagonal. This is order-one mixing, not the small-angle CKM pattern.


## Stability Analysis

**Scaled Hessian**: All 14 eigenvalues are non-negative relative to the largest eigenvalue (1.84e33), so the point appears stable in the scaled coordinates.

**Physical Hessian**: After converting to physical field units, the most negative eigenvalue is -1.33e27, which is dominated by the extreme conditioning of the B, B~, H fields. The 8 negative eigenvalues of the physical Hessian indicate this is actually a **saddle point** in the physical field space, not a true minimum. The apparent stability in scaled coordinates is an artifact of the per-field scaling absorbing the tachyonic directions into flat-looking directions.


## Interpretation

1. **The soft term dominates the vacuum structure.** With m~^2 = f_pi^2, the soft term cost m~^2 Tr(M†M) overwhelms the seesaw F-term structure. The optimizer finds configurations where all meson entries are as small as possible (O(10^4 MeV)) while still approximately satisfying the determinant constraint through the baryon VEVs.

2. **The Seiberg seesaw is destabilized.** The diagonal seesaw vacuum M_i = C/m_i has M_u ~ 4e5 MeV, which costs (4e5)^2 * 8464 ~ 1.4e15 MeV^2 in soft-term energy. This is the dominant contribution at the seesaw point. Any configuration with smaller |M| entries is energetically preferred, even at the cost of nonzero F-terms.

3. **Off-diagonal VEVs are preferred** because they allow det M to be satisfied (approximately, via BB~ absorption) with smaller total Tr(M†M). The minimum has |M^i_j| ~ 10^4 for all entries, whereas the seesaw concentrates the VEVs hierarchically.

4. **The CKM mixing angles are order-one**, not small. This is because the vacuum is completely rearranged from the seesaw, not a perturbation around it. The SVD mixing matrix shows maximal mixing, not small Cabibbo-like angles.

5. **Saddle point nature.** The physical Hessian has tachyonic directions, indicating the true global minimum may be at even lower V, possibly at the origin M = 0 (with V = Lambda^12 from the F_X constraint, but that gives V ~ 5.3e29 which is much larger). The landscape is complex.


## Implications for the sBootstrap

The key physical question is whether the **combined** soft-term + F-term potential has a minimum that:
(a) preserves the approximate Seiberg seesaw structure, and
(b) generates small off-diagonal perturbations that can be identified with CKM mixing.

The answer from this analysis is that m~^2 = f_pi^2 is too large relative to the Seiberg VEV scale. The soft term wants to collapse M toward zero, and the resulting vacuum has no resemblance to the seesaw. For the seesaw vacuum to be (approximately) preserved, one needs m~^2 << m_i^2 (the quark masses squared), i.e., the soft breaking must be treated as a perturbation on the F-flat directions rather than as a competing O(1) effect.

If instead one works in the regime m~^2 << m_i X cof(M,i) ~ m_i^2, the off-diagonal tachyonic modes (identified in the known diagonal-vacuum analysis) would lead to small off-diagonal VEVs proportional to m~/m_i, which could produce parametrically small CKM-like angles. That perturbative regime would require a different analysis from the brute-force global minimization done here.

## Code

Script: `results/offdiag_vacuum.py`
