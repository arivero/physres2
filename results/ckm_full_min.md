# Full Numerical Minimization of V = sum|F_I|^2 + f_pi^2 Tr(M^dag M)

## Setup

Superpotential for N_f = N_c = 3 SQCD with Higgs + mu-term:

    W = sum_i m_i M^i_i + X(det M - BB~ - Lambda^6)
        + y_c H_u^0 M^d_d + y_b H_d^0 M^s_s
        + lambda X (H_u^0 H_d^0 - H_u^+ H_d^-)

Parameters: m_u = 2.16, m_d = 4.67, m_s = 93.4 MeV; Lambda = 300 MeV; v = 246220 MeV;
y_c = 2 x 1270/v = 1.032e-2, y_b = 2 x 4180/v = 3.395e-2, lambda = 0.72; f_pi = 92 MeV.

Variables: 20 real = 9 complex meson entries M^i_j + complex X.
Higgs fixed at H_u^0 = H_d^0 = v/sqrt(2) = 174104 MeV, B = B~ = 0.


## Starting Point: Seesaw Vacuum

The diagonal Seiberg seesaw gives M_i = C/m_i with C = Lambda^2 (m_u m_d m_s)^{1/3} = 882297 MeV^2:

| Field | Value (MeV) |
|-------|-------------|
| M_u | 408471 |
| M_d | 188929 |
| M_s | 9446 |

At this vacuum with Higgs fixed at v/sqrt(2):

- V = 4.76e20 MeV^2
- Dominated by |F_X|^2 = (lam v^2/2)^2 = 4.76e20 (>99.99%)
- F_X = lam * (v/sqrt(2))^2 = 2.18e10 MeV^5 is the dominant F-term
- Soft term f_pi^2 Tr(M^dag M) = 1.72e15 is subdominant
- F_{M^d_d} = m_d + y_c v/sqrt(2) = 1796 MeV, F_{M^s_s} = m_s + y_b v/sqrt(2) = 5911 MeV


## Optimization Strategy

1. **L-BFGS-B + Powell** from seesaw starting point
2. **40 random restarts** at perturbation scales 1%, 5%, 10%, 50% around seesaw
3. **Basin hopping** (25 iterations)
4. **Differential evolution** (popsize 40, 400 iterations)
5. **Refinement**: 3 rounds of L-BFGS-B -> Powell -> Nelder-Mead on best result
6. **Directional search**: displacement along 3 softest Hessian eigenvectors

Total optimization time: ~100 seconds.


## Results

### Best Minimum

V = 4.219e14 MeV^2  (reduction 1.13e6x from starting point)

Meson matrix |M^i_j| at minimum (MeV):

|       |    u     |    d      |    s     |
|-------|----------|-----------|----------|
| **u** |  76054   |  45070    |  68846   |
| **d** |  19842   | 182424    |  24510   |
| **s** |  54296   |   3622    |   7742   |

Key observation: off-diagonal entries are large (off-diag/diag ratio = 0.84),
comparable to diagonal entries. This is NOT a small perturbation around the seesaw.

### Potential Decomposition

| Contribution | Value (MeV^2) | Fraction |
|-------------|---------------|----------|
| sum |F_M|^2 | 5.04e7 | 0.00% |
| |F_X|^2 | 1.90e8 | 0.00% |
| sum |F_H|^2 | 3.61e6 | 0.00% |
| f_pi^2 Tr(M^dag M) | **4.22e14** | **100.0%** |

The soft term completely dominates at the minimum. The optimizer has driven |F_X|
from 2.18e10 down to ~1.38e4 by adjusting det M, but can't reduce the soft term
further because the large diagonal VEVs M_d ~ 182000 MeV are required to
approximately satisfy the Seiberg constraint det M ~ Lambda^6.

Other fields at minimum:
- X = -2.81e-8 + 3.63e-7 i (shifted from X_0 = -1.21e-9)
- F_X = 3134 - 13433i (reduced 1.6e6x from starting value)
- F_{Hu0} = 1882 MeV (comparable to starting value)
- F_{Hd0} = 256 MeV


### CKM Mixing Angles from SVD

SVD: M = U_L diag(sigma) V_R^dag

Singular values: 190455, 111330, 34380 MeV
Seesaw values:   408471, 188929, 9446 MeV

V_CKM = U_L^dag V_R:

|           | column 1 | column 2 | column 3 |
|-----------|----------|----------|----------|
| **row 1** | 0.934    | 0.320    | 0.158    |
| **row 2** | 0.352    | 0.866    | 0.354    |
| **row 3** | 0.055    | 0.383    | 0.922    |

Mixing angles:

| Angle | This work | PDG | Oakes | Ratio |
|-------|-----------|-----|-------|-------|
| theta_12 | 18.9 deg | 13.04 deg | 12.6 deg | 1.45 |
| theta_23 | 21.0 deg | 2.38 deg | --- | 8.8 |
| theta_13 | 9.08 deg | 0.201 deg | --- | 45.2 |
| J | 6.1e-3 | 3.18e-5 | --- | 192 |

The mixing angles are order-one and hierarchically wrong. The 1-2 mixing (Cabibbo)
is within a factor of 1.5 of PDG, but 2-3 and 1-3 are 10-50x too large.
The Jarlskog invariant is 200x the PDG value.


### Hessian Analysis

Physical Hessian eigenvalues at the best minimum:

| Mode | m^2 (MeV^2) | Status |
|------|-------------|--------|
| 1 | -7.33e21 | **tachyonic** (X direction) |
| 2 | -2.29e10 | **tachyonic** (M^s_s, M^s_u) |
| 3 | -1.31e10 | **tachyonic** (M^u_s, M^s_d) |
| 4 | -7.09e9 | **tachyonic** (M^u_u, M^s_d) |
| 5 | -6.15e9 | **tachyonic** (M^s_d, M^d_s) |
| 6 | -5.18e9 | **tachyonic** (M^u_u, M^d_s) |
| 7 | -1.69e9 | **tachyonic** (M^u_d, M^d_u) |
| 8 | -1.42e9 | **tachyonic** (M^d_u, M^u_d) |
| 9 | -1.06e9 | **tachyonic** (M^d_d, Im M^d_d) |
| 10-17 | +1.2e9 to +1.4e10 | positive |
| 18-20 | ~1.3e21 | positive (stiff X/Higgs) |

**9 tachyonic directions** -- this is a saddle point, not a true minimum.
The two stiffest eigenvalues (~10^21) correspond to the X direction coupled
to the determinant constraint det M = Lambda^6 - lam Hu Hd.

The most tachyonic physical mode (m^2 = -7.33e21) is along Im(X)/Re(X),
consistent with F_X driving the potential.

The specific tachyonic direction reported by agent 10D -- Im(0.87 M^d_d + 0.49 M^u_u) --
has positive curvature 3.17e19 at this minimum. That direction was tachyonic at
the original seesaw vacuum, not at this displaced minimum.


### Catalog of Distinct Minima

73 distinct stationary points were found. The 10 lowest:

| # | V (MeV^2) | theta_12 | theta_23 | theta_13 | off/diag |
|---|-----------|----------|----------|----------|----------|
| 1 | 5.00e14 | 19.5 | 18.0 | 6.3 | 0.84 |
| 2 | 5.04e14 | 18.6 | 52.0 | 8.7 | 0.80 |
| 3 | 5.21e14 | 20.5 | 43.5 | 7.8 | 0.81 |
| 4 | 5.35e14 | 20.1 | 46.9 | 8.6 | 0.83 |
| 5 | 5.36e14 | 20.2 | 46.1 | 8.5 | 0.83 |
| 6 | 5.38e14 | 20.0 | 46.1 | 8.7 | 0.82 |
| 14 | 6.16e14 | **13.0** | 32.4 | 12.2 | 0.69 |

Minimum #14 has theta_12 = 13.0 deg matching PDG (13.04 deg) exactly,
but theta_23 = 32.4 deg is still far off. No minimum achieves the full CKM hierarchy.

All minima have off-diagonal/diagonal ratios of 0.7-1.0, meaning the
off-diagonal VEVs are as large as the diagonal ones. This is NOT the
perturbative regime. The optimizer is exploring a flat landscape of saddle
points where the soft term dominates.


## Physical Interpretation

### Why the Minimization Cannot Produce PDG CKM

1. **F_X dominance sets the scale.** At the seesaw vacuum with Higgs on, F_X = lam v^2/2 = 2.18e10 MeV^5
   gives |F_X|^2 = 4.76e20, dwarfing the soft term (1.72e15) by 5 orders of magnitude.
   The optimizer first tries to satisfy F_X = 0 by adjusting det M, which requires
   det M = Lambda^6 - lam v^2/2. Since lam v^2/2 = 2.18e10 >> Lambda^6 = 7.29e14,
   this is impossible: det M would need to be negative and enormous.

2. **The potential is a compromise.** The optimizer finds a saddle point where X shifts
   to partially cancel F_X, the meson matrix distorts from diagonal to reduce det M,
   and the soft term pays the price for the large VEVs needed.

3. **The off-diagonal distortion is non-perturbative.** Because the potential landscape
   near the seesaw is very flat (all eigenvalues ~10^3-10^9 vs V ~ 10^14-10^20),
   the optimizer slides far from the diagonal, producing order-one mixing angles
   rather than the small hierarchical angles of the PDG CKM.

4. **The Weinberg-Oakes relation is not recovered.** The Oakes relation theta_C = arctan(sqrt(m_d/m_s)) = 12.6 deg
   arises from the ratio of quark mass parameters, but the minimization mixes all sectors
   simultaneously with no mechanism to separate the 1-2 mixing from 2-3 and 1-3.

### What This Means

The tree-level scalar potential V = sum|F|^2 + f_pi^2 Tr(M^dag M) with B = B~ = 0
and Higgs fixed at v/sqrt(2) does NOT have a minimum that reproduces the CKM matrix.
The reason is that the lam X Hu Hd coupling creates an enormous F_X that overwhelms
the structure of the meson sector.

Possible resolutions:
1. The CKM structure emerges at a different scale (e.g., from the Coleman-Weinberg
   potential in the ISS framework, not from the tree-level potential)
2. The Higgs should not be fixed at v/sqrt(2) at the QCD scale -- it should be
   integrated in as a dynamical field with its own minimization
3. B, B~ should be included (the baryonic escape route), where BB~ = -Lambda^6
   allows det M = 0 and the M = 0 vacuum as the true minimum (see offdiag_vacuum.md)
4. The mu-term coupling lam should be much smaller, or the X field should be
   treated separately from the meson sector


## Code

Script: `results/ckm_full_min.py`
