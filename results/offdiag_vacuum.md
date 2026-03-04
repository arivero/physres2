# Off-Diagonal Meson Vacuum Analysis

## Setup

Superpotential for N_f = N_c = 3 SQCD with two Higgs doublets:

```
W = sum_i m_i M^i_i + X(det M - B B~ - Lambda^6) + y_c H_u^0 M^d_d + y_b H_d^0 M^s_s
```

with soft SUSY-breaking term V_soft = m~^2 Tr(M†M), where m~^2 = f_pi^2 = (92 MeV)^2 = 8464 MeV^2.

**Inputs**: m_u = 2.16, m_d = 4.67, m_s = 93.4 MeV; m_c = 1270, m_b = 4180 MeV; Lambda = 300 MeV; v = 246220 MeV; y_c = 2m_c/v = 1.032e-2, y_b = 2m_b/v = 3.395e-2.

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

At this vacuum with Hu = Hd = 0:
- V_0 = 1.715e15 MeV^2
- Dominated by the soft term m~^2(M_u^2 + M_d^2 + M_s^2) ~ 1.71e15 (99.7% of V)
- M_u = C/m_u ~ 4.1e5 MeV is enormous, contributing ~1.41e15 alone


## Key Result: The M = 0 Global Minimum

The true global minimum of the scalar potential is analytically solvable.

**Vacuum**: M = 0 (all meson VEVs vanish), X = 0, BB~ = -Lambda^6, H_u = -m_d/y_c, H_d = -m_s/y_b.

**Potential**: V_min = m_u^2 = 4.666 MeV^2.

### Derivation

At M = 0:
- All 3x3 cofactors vanish, so X drops out of F_{M^i_j}
- F_{M^i_j} = m_i delta_{ij} + y_c H_u delta_{i,d} delta_{j,d} + y_b H_d delta_{i,s} delta_{j,s}
- Setting H_u = -m_d/y_c and H_d = -m_s/y_b zeros F_{M^d_d} and F_{M^s_s}
- F_{M^u_u} = m_u remains uncancellable (no Yukawa coupling for up quark)
- det M = 0, but B B~ = -Lambda^6 satisfies F_X = det M - BB~ - Lambda^6 = 0
- X = 0 ensures F_B = -X B~ = 0, F_{B~} = -X B = 0
- F_{H_u} = y_c M^d_d = 0, F_{H_d} = y_b M^s_s = 0
- Soft term m~^2 Tr(M†M) = 0

The only nonzero F-term is F_{M^u_u} = m_u = 2.16 MeV, giving V = m_u^2 = 4.666 MeV^2.

### Stability

All Hessian eigenvalues at M = 0 are non-negative:

| Direction | d^2V/dphi^2 | Status |
|-----------|-------------|--------|
| M^i_j (all 9) | 2 m~^2 = 16928 | Stable |
| X | 2(B^2 + B~^2) = 2.92e15 | Stable |
| H_u | 2 y_c^2 = 2.13e-4 | Stable |
| H_d | 2 y_b^2 = 2.31e-3 | Stable |
| B, B~ | eigenvalues 0, 4 Lambda^6 | Flat + stable |

The single zero eigenvalue is the flat direction along the baryonic moduli space B B~ = -Lambda^6.

### Comparison of all vacua

| Vacuum | V (MeV^2) | Ratio to M=0 |
|--------|-----------|--------------|
| Seiberg seesaw (B=0, H=0) | 1.715e15 | 3.68e14 |
| Diagonal (soft-optimized) | 1.329e15 | 2.85e14 |
| 14D numerical (off-diag) | 2.169e13 | 4.65e12 |
| **M = 0 (analytical)** | **4.666** | **1** |


## Physical Interpretation

1. **The soft term dominates.** With m~^2 = f_pi^2 = 8464 MeV^2, the soft-breaking term m~^2 Tr(M†M) overwhelms the F-term structure. At the Seiberg seesaw vacuum, the soft term costs ~1.7e15 MeV^2, while the F-terms cost nothing (SUSY is unbroken there, before adding soft terms). The potential landscape is a steep bowl centered at M = 0, with F-term barriers preventing exact SUSY restoration.

2. **The baryonic flat direction resolves the constraint.** The determinant constraint det M = Lambda^6 (from F_X = 0) seems to require nonzero meson VEVs. But the full superpotential has F_X = det M - BB~ - Lambda^6, and BB~ = -Lambda^6 with det M = 0 satisfies the constraint. The baryons B, B~ act as a "pressure valve" that decouples the meson VEVs from the confinement scale.

3. **SUSY is broken by a single F-term: F_{M^u_u} = m_u.** The up quark mass is the irreducible source of SUSY breaking. The d and s masses are cancelled by the Higgs VEVs (H_u = -m_d/y_c, H_d = -m_s/y_b), but the up quark has no Yukawa coupling in this superpotential to cancel its mass term.

4. **No CKM mixing at the global minimum.** Since all meson VEVs vanish, M is the zero matrix. There are no off-diagonal entries and no CKM-like mixing angles. The SVD of the zero matrix is trivially the identity.

5. **The Seiberg seesaw is a false vacuum.** The seesaw vacuum M_i = C/m_i sits at V ~ 10^{15} MeV^2, some 14 orders of magnitude above the global minimum. The enormous hierarchy is entirely due to the soft term cost of the large M_u VEV.


## Implications for the sBootstrap

### The problem with m~^2 = f_pi^2

The identification m~^2 = f_pi^2 makes the soft term a strong competitor to the F-terms. Since f_pi ~ Lambda/3, the soft scale is comparable to the confinement scale, and the meson VEVs (which scale as Lambda^2/m_quark >> Lambda) are heavily penalized. This drives the vacuum to M = 0, which is physically unacceptable (no quark condensate, no chiral symmetry breaking).

### What would be needed for a nontrivial meson vacuum

For the Seiberg seesaw to survive as the true vacuum, one needs:

1. **Much smaller soft breaking**: m~^2 << m_quark^2 for the lightest quark. Since m_u ~ 2 MeV, this would require m~ << 2 MeV. With f_pi = 92 MeV, the soft term is 40x too large.

2. **Remove the baryonic flat direction**: If BB~ is somehow stabilized at zero (e.g., by a baryon number symmetry that prevents the BB~ = -Lambda^6 solution), then the constraint becomes det M = Lambda^6, and M cannot be zero. This would restore the seesaw vacuum. However, in the standard Seiberg-Intriligator-Shenker framework, BB~ is a modulus.

3. **Treat the soft term perturbatively**: If m~^2 << Lambda^2, the soft term is a small perturbation on the Seiberg vacuum. The meson mass matrix at the seesaw vacuum has off-diagonal tachyonic directions (ds, us sectors) as was known. The perturbative off-diagonal VEVs would scale as:
   delta M^i_j ~ m~^2 M^i_i / (X * cofactor scale) ~ m~^2 / m_quark

   This could give small, calculable CKM-like angles.


## Numerical exploration (intermediate minima)

For completeness, the numerical exploration (46 initial conditions) found local minima in the V ~ 10^{13} range when starting from perturbed seesaw vacua. At these intermediate minima:

- Off-diagonal VEVs are order-one relative to diagonal VEVs
- CKM-like mixing angles are 16--76 degrees (order-one mixing, not small-angle CKM)
- The diagonal VEVs are completely rearranged from the seesaw (M_u negative, M_d reduced by 84%)
- The physical Hessian shows tachyonic directions (these are saddle points, not true minima)

These intermediate stationary points are en route from the seesaw vacuum to the global M = 0 minimum.


## Conclusion

The scalar potential V = sum|F|^2 + f_pi^2 Tr(M†M) has its global minimum at M = 0 with V = m_u^2 = 4.67 MeV^2. This vacuum has no meson condensate and no CKM mixing angles. The result is a consequence of the baryonic flat direction BB~ = -Lambda^6 decoupling the determinant constraint from the meson sector, combined with the soft term driving M to zero.

For CKM mixing to emerge from off-diagonal meson VEVs, the model requires either:
(a) eliminating the baryonic escape route (BB~ = 0 enforced), or
(b) a soft-breaking scale small enough to treat as a perturbation on the Seiberg vacuum.

## Code

Script: `results/offdiag_vacuum.py`
