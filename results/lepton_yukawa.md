# Lepton Yukawa from SQCD Meson Fermion Mass Matrix

## Setup

SU(3) SQCD with N_f = N_c = 3 flavors (u, d, s). Seiberg duality gives mesons
M^i_j and baryons B, Btilde. The superpotential is:

W = sum_i m_i M^i_i + X(det M - B Btilde - Lambda^6) + y_c Hu0 M^d_d + y_b Hd0 M^s_s

Parameters: m_u = 2.16, m_d = 4.67, m_s = 93.4 MeV, Lambda = 300 MeV,
v = 246220 MeV, y_c = 1.031598e-02, y_b = 3.395338e-02

Vacuum: M^i_i = C/m_i, X = -C/Lambda^6, Hu0 = Hd0 = v/sqrt(2)
where C = Lambda^2 (m_u m_d m_s)^(1/3) = 882297.425722 MeV^2

## 1. Central-Block Fermion Mass Matrix (6x6)

Fields: (M^u_u, M^d_d, M^s_s, X, Hu0, Hd0)

Nonzero entries:
- W[M^u_u, M^d_d] = X0 * M_s = -1.143288e-05
- W[M^u_u, M^s_s] = X0 * M_d = -2.286576e-04
- W[M^d_d, M^s_s] = X0 * M_u = -4.943662e-04
- W[M^u_u, X] = M_d * M_s = 1.784704e+09
- W[M^d_d, X] = M_u * M_s = 3.858597e+09
- W[M^s_s, X] = M_u * M_d = 7.717194e+10
- W[M^d_d, Hu0] = y_c = 1.031598e-02
- W[M^s_s, Hd0] = y_b = 3.395338e-02

## 2. Eigenvalues

| # | Eigenvalue (MeV) | |Eigenvalue| (MeV) |
|---|-----------------|---------------------|
| 1 | -7.69409336e-04 | 7.69409336e-04 |
| 2 | 7.79707136e-04 | 7.79707136e-04 |
| 3 | -1.04177239e-02 | 1.04177239e-02 |
| 4 | 1.04670166e-02 | 1.04670166e-02 |
| 5 | 7.72889485e+10 | 7.72889485e+10 |
| 6 | -7.72889485e+10 | 7.72889485e+10 |

## 3. Light Eigenstates

### Eigenvalue 1: -7.694093e-04 MeV

| Field | Component |
|-------|-----------|
| M^u_u | 0.70452138 |
| M^d_d | -0.00858109 |
| M^s_s | -0.01586395 |
| H_u^0 | 0.11505230 |
| H_d^0 | 0.70006235 |

### Eigenvalue 2: 7.797071e-04 MeV

| Field | Component |
|-------|-----------|
| M^u_u | -0.70921918 |
| M^d_d | 0.00874833 |
| M^s_s | 0.01596423 |
| H_u^0 | 0.11574511 |
| H_d^0 | 0.69518331 |

### Eigenvalue 3: -1.041772e-02 MeV

| Field | Component |
|-------|-----------|
| M^u_u | -0.00710744 |
| M^d_d | -0.70533982 |
| M^s_s | 0.03543136 |
| H_u^0 | 0.69846561 |
| H_d^0 | -0.11548005 |

### Comparison to lepton masses

| Eigenvalue (MeV) | Lepton mass (MeV) | Ratio |
|------------------|--------------------|-------|
| 7.694093e-04 | 0.510999 | 6.6414e+02 |
| 7.797071e-04 | 105.658376 | 1.3551e+05 |
| 1.041772e-02 | 1776.860000 | 1.7056e+05 |

The light eigenvalue ratios do NOT match lepton mass ratios.
A universal rescaling cannot relate them.

## 4. Koide Test

Q(light eigenvalues) = 0.4810142645
Deviation from 2/3: -27.8479%

## 5. NMSSM Extension

Adding W_NMSSM = lambda X Hu0 Hd0 with lambda = 0.72:

New entries: W[X, Hu0] = W[X, Hd0] = lambda v/sqrt(2) = 1.2535e+05,
W[Hu0, Hd0] = lambda X0 = -8.7142e-10

| # | Base (MeV) | NMSSM (MeV) | Change |
|---|-----------|-------------|--------|
| 1 | -7.694093e-04 | -7.694089e-04 | +0.0001% |
| 2 | 7.797071e-04 | 7.797067e-04 | -0.0001% |
| 3 | -1.041772e-02 | -1.041847e-02 | -0.0072% |
| 4 | 1.046702e-02 | 1.046634e-02 | -0.0064% |
| 5 | 7.728895e+10 | 7.728741e+10 | -0.0020% |
| 6 | -7.728895e+10 | -7.728741e+10 | +0.0020% |

NMSSM Koide Q = 0.4810220492 (-27.8467% from 2/3)

## 6. Effective Lepton Yukawas

If the three light mesinos are identified with (e, mu, tau), their effective Yukawas
would be y = sqrt(2) m_eigenvalue / v:

| Light eigenvalue (MeV) | y_eff | SM y_lepton |
|----------------------|-------|-------------|
| 7.694093e-04 | 4.419256e-09 | 2.935024e-06 (y_e) |
| 7.797071e-04 | 4.478403e-09 | 6.068699e-04 (y_mu) |
| 1.041772e-02 | 5.983627e-08 | 1.020575e-02 (y_tau) |

## 7. Off-Diagonal Meson Alternative

The off-diagonal meson pairs (M^i_j, M^j_i) have fermion mass eigenvalues
+/- |X0| M_k where k is the third flavor. These are proportional to 1/m_k.

| Pair | Eigenvalue (MeV) | Proportional to |
|------|-----------------|-----------------|
| (M^u_d, M^d_u) | 1.143288e-05 | 1/m_s |
| (M^u_s, M^s_u) | 2.286576e-04 | 1/m_d |
| (M^d_s, M^s_d) | 4.943662e-04 | 1/m_u |

Q(1/m_u, 1/m_d, 1/m_s) = 0.4425755922
Deviation from 2/3: -33.614%

For comparison: Q(1/m_d, 1/m_s, 1/m_b) = 0.6652098709 (0.22% from 2/3)

## Key Findings

1. The 6x6 central-block fermion mass matrix has a natural heavy/light split.
   The heavy eigenvalues are set by the meson VEVs (C/m_i ~ 10^3-10^5 MeV)
   while the light eigenvalues arise from the seesaw with the Yukawa couplings.

2. The Koide quotient for the three light eigenvalues is Q = 0.48101426.
   This is -27.8479% from 2/3.

3. The light eigenvalue mass ratios do NOT match lepton mass ratios.
   No universal rescaling can connect them to (e, mu, tau).

4. The NMSSM extension (lambda X Hu Hd) does not qualitatively change the spectrum
   for Lambda = 300 MeV, because the NMSSM shift to det M is tiny compared to Lambda^6.

5. The off-diagonal meson eigenvalues are proportional to 1/m_k. Their Koide quotient
   Q(1/m_u, 1/m_d, 1/m_s) = 0.442576 is far from 2/3.
   The dual Koide Q(1/m_d, 1/m_s, 1/m_b) = 0.665210 (0.22% from 2/3) uses
   the down-type quarks (d, s, b), not the light quarks (u, d, s).

