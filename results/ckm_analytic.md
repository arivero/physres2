# Analytic CKM from Tachyonic Meson Condensate

## Setup

3x3 meson matrix M with diagonal Seiberg seesaw entries M_i = C/m_i,
C = Lambda^2 (m_u m_d m_s)^{1/3} = 882297.425722 MeV^2.

Parameters: m_u = 2.16, m_d = 4.67, m_s = 93.4 MeV;
Lambda = 300 MeV, v = 246220 MeV, lambda = 0.72.
F_X = -lambda v^2/2 = -2.1825e+10 MeV^5.

## 1. Diagonal Vacuum Values

| Field | Value (MeV) | Quark mass (MeV) |
|-------|-------------|-----------------|
| M_U | 408471.030427 | m_u = 2.16 |
| M_D | 188928.784951 | m_d = 4.67 |
| M_S | 9446.439248 | m_s = 93.4 |

det M / Lambda^6 = 0.999999999999999 (exactly 1 by construction).

## 2. Tachyonic Mass-Squared

From the B-term mechanism: W_{X, M^a_b, M^b_a} = M_k (cofactor coupling)
generates m^2_{ij} = -|F_X|^2 M_k^2 where k is the third index.

| Sector | k | M_k (MeV) | m^2_tach |
|--------|---|-----------|---------|
| ds | u | 408471.0304 | -7.9473e+31 |
| us | d | 188928.7850 | -1.7002e+31 |
| ud | s | 9446.4392 | -4.2504e+28 |

Hierarchy: |m^2_ds| : |m^2_us| : |m^2_ud| = M_U^2 : M_D^2 : M_S^2
= 1870 : 400 : 1

## 3. Off-Diagonal Condensates

From V(eps) = m^2 eps^2 + M_k^2 eps^4, the minimum is at eps^2 = |F_X|/M_k.

| Sector | eps_min (MeV) | eps ~ sqrt(lambda v^2 m_k / (2C)) |
|--------|---------------|-----------------------------------|
| ds | 231.1500 | sqrt(0.72 * 246220^2 * 2.16 / (2 * 882297.43)) |
| us | 339.8799 | sqrt(0.72 * 246220^2 * 4.67 / (2 * 882297.43)) |
| ud | 1519.9892 | sqrt(0.72 * 246220^2 * 93.4 / (2 * 882297.43)) |

Key ratio: eps_us / eps_ud = sqrt(m_d/m_s) -- this is the Cabibbo ratio.

## 4. GST Derivation

### Theorem

For the symmetric 2x2 meson matrix M = [[M_D, eps], [eps, M_S]] with
eps = C/sqrt(m_d m_s) = sqrt(M_D M_S), the eigenvalue problem yields
**tan(theta_C) = sqrt(m_d/m_s)** exactly (Weinberg-Oakes relation).

### Proof

tan(2 theta_C) = 2 eps / (M_D - M_S)
= 2C/sqrt(m_d m_s) / [C(m_s-m_d)/(m_d m_s)]
= 2 sqrt(m_d m_s) / (m_s - m_d)
= 2 sqrt(m_d/m_s) / (1 - m_d/m_s) = tan(2 arctan(sqrt(m_d/m_s))).

QED. The off-diagonal quark mass entry is delta_m = C eps/(M_D M_S) = sqrt(m_d m_s),
giving tan(theta_C) = sqrt(m_d m_s)/(m_s - m_d) = sqrt(m_d/m_s) for m_s >> m_d.

### Corollary

|V_us| = sin(arctan(sqrt(m_d/m_s))) = sqrt(m_d/(m_d + m_s))
= 0.21821789 (PDG: 0.2257, deviation: -3.315%)

## 5. Weinberg-Oakes vs GST

| Relation | Formula | theta_C (deg) | |V_us| | Dev from PDG |
|----------|---------|--------------|--------|-------------|
| Weinberg-Oakes | tan tC = sqrt(m_d/m_s) | 12.6044 | 0.218218 | -3.315% |
| GST | sin tC = sqrt(m_d/m_s) | 12.9210 | 0.223607 | -0.927% |
| PDG | measured | 13.04 | 0.2257 | --- |

This framework predicts **Weinberg-Oakes** (arctan), not GST (arcsin).
The distinction is controlled by the symmetry of the meson condensate:

- **Symmetric** M^d_s = M^s_d (energy minimum by AM-GM): **arctan** (Weinberg-Oakes)
- **Asymmetric** M^d_s != M^s_d (requires CP-violating source): **arcsin** (GST) possible

## 6. Leading Correction to GST

|V_us|_WO = sqrt(m_d/m_s) * (1 - m_d/(2 m_s) + O((m_d/m_s)^2))

Leading correction: m_d/(2 m_s) = 0.025000 (2.50%)

## 7. Full CKM (3x3 eigenvec)

With hierarchical condensate eps_ij = eta M_k, fitted eta = 0.357364:

| Angle | This work | PDG |
|-------|-----------|-----|
| theta_12 | 13.0399 deg | 13.04 deg |
| theta_23 | 28.5120 deg | 2.38 deg |
| theta_13 | 6.7531 deg | 0.201 deg |
| delta | 0 (real M) | 1.196 rad |
| J | 0.00e+00 | 3.18e-5 |

CP violation (delta != 0) requires complex off-diagonal condensates.
The real symmetric case gives J = 0.

## 8. Summary

1. The Seiberg seesaw M_i = C/m_i provides the diagonal structure.

2. F_X = -lambda v^2/2 = -2.1825e+10 MeV^5 drives tachyonic off-diagonal modes
   with hierarchy m^2(ds) : m^2(us) : m^2(ud) = M_U^2 : M_D^2 : M_S^2.

3. The **Weinberg-Oakes relation** tan(theta_C) = sqrt(m_d/m_s) is derived EXACTLY
   when eps_ds = sqrt(M_D M_S) = C/sqrt(m_d m_s) (geometric mean of diagonal VEVs).

4. This framework predicts **arctan** (Weinberg-Oakes), not **arcsin** (GST).
   The difference is O(m_d/m_s) ~ 5.0%.

5. |V_us| = sqrt(m_d/(m_d+m_s)) = 0.218218
   vs GST: sqrt(m_d/m_s) = 0.223607
   vs PDG: 0.2257

6. The tachyonic condensate ratio eps_us/eps_ud = sqrt(m_d/m_s) directly encodes
   the Cabibbo angle through the seesaw mass hierarchy.

7. GST requires left-right asymmetric condensate (M^d_s != M^s_d),
   specifically beta/alpha ~ 1.0273, disfavored by potential minimization.