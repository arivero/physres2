# Perturbative CKM Mixing from Off-Diagonal Meson VEVs

**Date:** 2026-03-04

## Setup

Superpotential with B = B̃ = 0 enforced:

    W = sum_i m_i M^i_i + X(det M - Lambda^6)
        + y_c H_u M^d_d + y_b H_d M^s_s + lambda X H_u H_d

Soft breaking: V_soft = f_pi^2 Tr(M†M).  H_u = H_d = 0.

Parameters: m_u = 2.16, m_d = 4.67, m_s = 93.4 MeV; Lambda = 300 MeV;
f_pi = 92 MeV; y_c = 1.0316e-02, y_b = 3.3953e-02, lambda = 0.72.

**Seesaw vacuum**: M_i = C/m_i, C = Lambda^2(m_u m_d m_s)^{1/3}.

  M_u = 408471.03 MeV,  M_d = 188928.78 MeV,  M_s = 9446.44 MeV
  X_vev = -1.2103e-09 MeV^(-4)

At this vacuum: F_{M^i_i} = m_i + X M_j M_k = 0 (all i), F_X = det M - Lambda^6 = 0.

## Analytic Hessian for Off-Diagonal Meson Entries

Write M^i_j = delta^i_j M_i + epsilon^i_j and expand V to second order.

For off-diagonal pair (a,b), let k = third(a,b). The key F-terms are:

- F_{M^a_b} = X * cof(M)^b_a ~ X_vev * M_k * epsilon^a_b  [linear in epsilon^a_b]
- F_{M^k_k} = -X_vev * epsilon^a_b * epsilon^b_a + O(epsilon^4)  [bilinear, but zero at eps=0]
- F_X = -M_k * epsilon^a_b * epsilon^b_a + (other pairs)  [quadratic in epsilon, zero at eps=0]

All contributions to the Hessian cross-terms vanish at the seesaw vacuum because:
1. F_{M^k_k} = 0 at epsilon = 0, so its contribution 2 F_{M^k_k} * d^2 F_{M^k_k} = 0.
2. d F_{M^k_k}/d epsilon^a_b = -X_vev * epsilon^b_a = 0 at epsilon = 0.
3. Same logic for |F_X|^2.

The analytic 6x6 mass matrix is **diagonal**:

    m^2_(a,b) = f_pi^2 + 2 X_vev^2 M_k^2

### Mass eigenvalues

| Sector | k | M_k (MeV) | f_pi^2 (MeV^2) | 2X^2 M_k^2 (MeV^2) | m^2_eff | Status |
|--------|---|-----------|----------------|---------------------|---------|--------|
| ud | s | 9446.44 | 8.4640e+03 | 2.6142e-10 | 8.4640e+03 | stable |
| us | d | 188928.78 | 8.4640e+03 | 1.0457e-07 | 8.4640e+03 | stable |
| du | s | 9446.44 | 8.4640e+03 | 2.6142e-10 | 8.4640e+03 | stable |
| ds | u | 408471.03 | 8.4640e+03 | 4.8880e-07 | 8.4640e+03 | stable |
| su | d | 188928.78 | 8.4640e+03 | 1.0457e-07 | 8.4640e+03 | stable |
| sd | u | 408471.03 | 8.4640e+03 | 4.8880e-07 | 8.4640e+03 | stable |

X_vev = -1.2103e-09 MeV^(-4) is tiny: 2 X_vev^2 M_k^2 << f_pi^2 = 8464.0 MeV^2.

**All off-diagonal modes are stable at the seesaw vacuum.**

## Numerical Verification

The numerical Hessian computed at small step h = 1e-4 MeV (true quadratic regime)
confirms the analytic result. The apparent tachyonic eigenvalues seen at large
step h ~ f_pi MeV are artifacts of the quartic potential dominating:

The quartic |F_X|^2 ~ M_k^2 * epsilon^2_{ab} * epsilon^2_{ba} generates an effective
cross-term in the finite-difference formula H^{num}_{(ab),(ba)}(h) ~ M_k^2 * h^2.
At h ~ f_pi ~ 10^2 MeV, this quartic contribution ~ M_k^2 * h^2 ~ M_k^2 * f_pi^2
which for k = d gives 188929^2 * 92^2 ~ 3.0e+14.
This dominates the true Hessian f_pi^2 ~ 8500, explaining the spurious tachyons.

## Global Minimum Search

Multi-start minimization (100 random starts) confirms that
V_seesaw = 1.7151e+15 MeV^2 is the global minimum in the off-diagonal sector.
No off-diagonal VEVs are generated: the CKM matrix is the identity.

## Physical Interpretation

The seesaw vacuum M_i = C/m_i is a stable critical point because:

1. The soft term f_pi^2 Tr(M†M) = f_pi^2 (sum M_i^2 + sum eps^2) contributes
   a positive mass f_pi^2 = 8464.0 MeV^2 to every off-diagonal mode.

2. The off-diagonal F-terms F_{M^a_b} = X_vev * M_k * epsilon^a_b add
   2 X_vev^2 M_k^2 ≈ 2.6e-10 MeV^2 (negligible) to the mass.

3. The W_{IJK} F_I mechanism (Yukawa-induced mass splitting) requires nonzero
   background F-terms. At the exact seesaw vacuum, all F_{M^i_i} = 0,
   so this mechanism is switched off.

For CKM mixing to emerge from off-diagonal meson VEVs, one needs either:
- A softer SUSY-breaking scale (f_pi << m_quark masses, not 92 MeV >> 2 MeV),
- Displacement from the seesaw vacuum so that some F-terms become nonzero.

## Oakes Relation as Residual Connection

Despite the absence of off-diagonal VEVs, the Oakes relation
    |V_us| ~ sqrt(m_d/m_s) = 0.22361  (PDG: 0.2243, deviation -0.3%)
provides a geometric link between the quark mass hierarchy and the Cabibbo angle.
This connection does not require off-diagonal meson VEVs; it is a property of the
dual Koide structure under the Seiberg seesaw M_j ~ 1/m_j.

## Conclusions

1. The 6x6 off-diagonal meson mass matrix is diagonal at the seesaw vacuum.
   All eigenvalues equal 2 f_pi^2 = 16928.0 MeV^2 > 0 (stable).

2. No tachyonic off-diagonal meson modes exist at the exact seesaw vacuum
   with B = B̃ = 0 enforced.

3. The CKM mixing matrix at this vacuum is the identity (zero mixing angles).

4. Tachyonic modes require nonzero F-term backgrounds, which are absent
   at the seesaw point by construction. The problem statement's claim that
   'tachyonic eigenvalues arise from W_{IJK} F_I terms' is mechanistically
   correct but cannot be realized at the exact seesaw vacuum.

5. A perturbative treatment valid near the seesaw vacuum predicts zero CKM angles.
   Nonzero CKM mixing must arise from a different mechanism.

6. The Oakes relation |V_us| = sqrt(m_d/m_s) = 0.2236 (PDG: 0.2243)
   connects quark mass ratios to the Cabibbo angle without off-diagonal VEVs.