# Combined Vacuum Analysis: SQCD + NMSSM X-Higgs Coupling

**Computation date:** 2026-03-04

## Setup

N=1 SU(3) SQCD with N_f = N_c = 3 and superpotential:

W = sum_i m_i M^i_i + X(det M - BB~ - Lambda^6) + lambda X H_u . H_d + m_B BB~

Soft SUSY breaking: V_soft = f_pi^2 Tr(M^dag M)

Parameters: m_u = 2.16 MeV, m_d = 4.67 MeV, m_s = 93.4 MeV, Lambda = 300 MeV, f_pi = 92 MeV, lambda = 0.72, v = 246220 MeV, m_B = 300 MeV.

## Problem 1: EWSB-Shifted Seesaw Vacuum

### Dimensional Analysis

The NMSSM coupling lambda X H_u H_d requires careful dimensional analysis:

- [M^i_j] = MeV^2 (composite meson)
- [X] = MeV^{-3} (Lagrange multiplier, so that [X det M] = MeV^3 = [W])
- [lambda X H_u H_d] must equal MeV^3

For [H_u] = [H_d] = MeV (canonical scalars), [lambda] must be MeV^4. Writing lambda_phys = lambda_hat Lambda^4 with lambda_hat = 0.72 dimensionless:

| Quantity | Value | Units |
|----------|-------|-------|
| Lambda^6 | 7.290e+14 | MeV^6 |
| lambda_hat | 0.72 | -- |
| lambda_phys = lambda_hat Lambda^4 | 5.832e+09 | MeV^4 |
| lambda_phys v^2/2 | 1.768e+20 | MeV^6 |
| Delta = lambda_phys v^2/(2 Lambda^6) | 2.425e+05 | -- |

**Delta >> 1.** The EWSB shift completely overwhelms Lambda^6. This means:
- If the constraint is det M = Lambda^6 - lambda_phys v^2/2, the shifted det M is deeply negative
- For the constraint to make physical sense with positive meson VEVs, lambda_hat must be extremely small (lambda_hat < 3.0e-06)
- This is a strong constraint on the NMSSM coupling: the X-Higgs coupling must be suppressed by at least (Lambda/v)^2 ~ 1.5e-06

### Unshifted Seiberg Seesaw

| Quantity | Value |
|----------|-------|
| C = Lambda^2 (m_u m_d m_s)^{1/3} | 882297 MeV^2 |
| M_u = C/m_u | 408471 MeV^2 |
| M_d = C/m_d | 188929 MeV^2 |
| M_s = C/m_s | 9446 MeV^2 |
| X_0 = -C/Lambda^6 | -1.210e-09 MeV^{-4} |
| det M / Lambda^6 | 1.000000000000 |

### Shifted Seesaw (benchmark Delta = 0.01)

| Quantity | Shifted | Original | Fractional shift |
|----------|---------|----------|-----------------|
| det M | 7.217e+14 | 7.290e+14 | -1.00% |
| C' | 879347 | 882297 | -0.33% |
| M_u' | 407105 | 408471 | -0.33% |
| M_d' | 188297 | 188929 | -0.33% |
| M_s' | 9415 | 9446 | -0.33% |

For small Delta, C' = C(1-Delta)^{1/3} and all meson VEVs shift uniformly by factor (1-Delta)^{1/3}.

## Problem 2-4: Off-Diagonal Meson Mass-Squared Matrix

### Two Interpretations of X

**(A) X = Lagrange multiplier** (Seiberg electric, N_f = N_c = 3):
- X has NO kinetic term
- Its equation of motion IS the constraint det M - BB~ = Lambda^6
- The constraint is modified by EWSB: det M = Lambda^6 - lambda_phys v^2/2
- F_X = 0 at the shifted vacuum (constraint satisfied by construction)

**(B) X = dynamical field** (hypothetical, ISS-like):
- X has a canonical Kahler potential
- The vacuum is the original seesaw with det M = Lambda^6
- F_X = lambda_phys v^2/2 is nonzero (SUSY broken by X)

### Canonical Variables

To properly compare terms with different field dimensions, define:
- Phi_{ij} = M^i_j / Lambda (dim MeV, canonical)
- X_c = X Lambda^4 (dim MeV, canonical)

| Canonical field | Value | Units |
|----------------|-------|-------|
| Phi_u | 1361.6 | MeV |
| Phi_d | 629.8 | MeV |
| Phi_s | 31.5 | MeV |
| X_c | -9.80 | MeV |
| F_{X_c} = lambda v^2/2 | 2.182e+10 | MeV^2 |

### Key Derivative: d^2(det M)/d(M^a_b)d(M^b_a) = -M_c

Verified explicitly from the Levi-Civita expansion of the determinant. The minus sign comes from the signature of the transposition permutation (ab).

### Off-Diagonal Mass-Squared: Case A

F_X = 0. No B-terms. All masses positive.

| Pair | m^2_F = |X_c/Lambda|^2 Phi_c^2 | m^2_soft = f_pi^2 | m^2_total |
|------|-----|---------|---------|
| u-d (spectator s) | 1.06 | 8464 | 8465 MeV^2 |
| u-s (spectator d) | 424 | 8464 | 8888 MeV^2 |
| d-s (spectator u) | 1980 | 8464 | 10444 MeV^2 |

**No tachyons. Vacuum is stable. No CKM mixing from meson VEVs.**

### Off-Diagonal Mass-Squared: Case B

F_{X_c} = lambda v^2/2 = 2.182e+10 MeV^2 (nonzero). B-terms from holomorphic mass:

B_{(M^a_b)(M^b_a)} = F_{X_c} * (-Phi_c / Lambda)

| Pair | m^2_diag | |B| | m^2_min = m^2_diag - |B| | Tachyonic? |
|------|----------|-----|---------|-----------|
| u-d (spectator s) | 8465 | 2.29e+09 | -2.29e+09 | YES |
| u-s (spectator d) | 8888 | 4.58e+10 | -4.58e+10 | YES |
| d-s (spectator u) | 10444 | 9.91e+10 | -9.91e+10 | YES |

The B-term exceeds the diagonal mass by 5-7 orders of magnitude. All sectors are deeply tachyonic.

**B-term hierarchy:**
- B(d-s) : B(u-s) : B(u-d) = Phi_u : Phi_d : Phi_s = M_u/M_s : M_d/M_s : 1 = 43.2 : 20.0 : 1

This is proportional to 1/m_quark (seesaw), NOT to sqrt(m_d/m_s). The tachyon does NOT reproduce the Oakes relation.

## Problem 5: Which Interpretation is Correct?

**X is a Lagrange multiplier. Case A is correct.**

Three independent arguments:

1. **Seiberg's result (1994):** For N_f = N_c, the low-energy description is in terms of gauge-invariant composites M, B, B~ subject to the quantum modified constraint det M - BB~ = Lambda^{2N_c}. The Lagrange multiplier X enforces this constraint. It has no kinetic term.

2. **No magnetic dual:** Seiberg duality maps N_f = N_c to N_c' = N_f - N_c = 0. There is no magnetic gauge group. The ISS mechanism requires N_f > N_c to produce a dynamical pseudo-modulus. For N_f = N_c, there is no magnetic theory and no dynamical X.

3. **Smooth quantum moduli space:** The constraint det M - BB~ = Lambda^6 defines a smooth variety (no singularities). The instanton-generated ADS superpotential W_dyn = Lambda^{2N_c+1}/det M enforces this at the quantum level.

**Consequence:** F_X = 0 at the vacuum. No B-terms. No tachyons. Case A is the physical answer.

## Problem 6: Soft-Term Tachyons?

At the shifted seesaw with F_X = 0:

m^2(Phi^a_b) = f_pi^2 + |X_c/Lambda|^2 Phi_c^2 > 0

Both terms are positive. The soft mass f_pi^2 = 8464 MeV^2 dominates over the F-term contribution (which is at most ~2000 MeV^2 for the d-s sector). No tachyonic directions exist from soft terms alone.

## Summary Assessment

### 1. Is the off-diagonal meson sector tachyonic?

**No.** For N_f = N_c = 3, X is a Lagrange multiplier. F_X = 0 at the vacuum. All off-diagonal meson masses are positive: m^2 = f_pi^2 + O(10^{0}) MeV^2 ~ 8464 MeV^2.

### 2. Source of CKM mixing (Case A)?

With X as Lagrange multiplier, there is NO tree-level CKM mixing from the meson vacuum. The vacuum is exactly diagonal. Possible sources of CKM mixing:

- Radiative corrections (one-loop Coleman-Weinberg potential)
- Higher-dimensional operators in the Kahler potential
- The Oakes relation as an INPUT from the quark mass matrix texture
- Non-perturbative effects (bion contributions)

### 3. Do Case B tachyons reproduce the Oakes relation?

No. The B-term hierarchy is B proportional to M_c = C/m_c (inverse quark mass). The ratios are m_s/m_u = 43.2 and m_s/m_d = 20.0, not sqrt(m_d/m_s) = 0.224. The tachyonic instability does not reproduce the Oakes relation.

However, qualitatively, the STRONGEST tachyon is in the d-s sector (where the spectator is the u quark, the lightest). This correctly identifies d-s as the sector with the largest flavor mixing, consistent with the large Cabibbo angle.

### 4. Status of the Cabibbo angle derivation

The Cabibbo angle does NOT emerge from the tree-level meson vacuum in either interpretation of X. The SQCD meson sector provides the mass hierarchy (seesaw M_i = C/m_i) and Koide structure (CW spectrum), but does not by itself generate CKM mixing.

The Oakes relation tan(theta_C) = sqrt(m_d/m_s) connects the Cabibbo angle to quark mass ratios, but in this framework the quark masses are inputs (m_i M^i_i in the superpotential). The connection between m_d/m_s and the CKM angle requires additional structure: a mass matrix texture, radiative corrections, or non-perturbative effects beyond the confined description.

**Bottom line:** The meson vacuum provides mass hierarchy and Koide. CKM mixing is a separate structural element.

### Critical constraint on the NMSSM coupling

For the shifted seesaw to have physical (positive) meson VEVs, the dimensionless coupling lambda_hat must satisfy:

lambda_hat < Lambda^2 / (v^2/2) = 2.97e-06

This is an extremely small coupling. If lambda_hat ~ O(1), the EWSB shift destroys the seesaw vacuum entirely (det M becomes deeply negative). This suggests the X-Higgs coupling is either:
- Absent (no direct X-Higgs coupling in the confined theory)
- Generated radiatively (suppressed by loop factors)
- Mediated through a separate singlet field (not X itself)
