# One-loop effective Kahler potential in the O'Raifeartaigh model

## Model and conventions

The O'Raifeartaigh model has superpotential

W = f X + m phi phi_tilde + g X phi^2

with chiral superfields X, phi, phi_tilde and real parameters f, m, g > 0.  The tree-level scalar potential V = |dW/dX|^2 + |dW/dphi|^2 + |dW/dphi_tilde|^2 breaks SUSY via F_X = f.  The field X is a pseudo-modulus: classically flat at phi = phi_tilde = 0.

Dimensionless variables: v = g <X>/m (pseudo-modulus VEV), y = gf/m^2 (SUSY-breaking parameter).  All masses are given in units of m throughout.  The CW potential is given in units of m^4/(64 pi^2).

## 1. Mass matrices

**Fermion mass matrix** in the basis (psi_X, psi_phi, psi_phi_tilde) at phi = phi_tilde = 0:

    M_F / m = [[0, 0, 0],
               [0, 2v, 1],
               [0, 1, 0]]

The X-fermion (Goldstino) decouples with mass zero.  The 2x2 block has eigenvalues

    m_+/m = v + sqrt(v^2 + 1),    m_-/m = sqrt(v^2 + 1) - v

with product m_+ m_- = m^2 and sum m_+ + m_- = 2m sqrt(v^2+1).

**Scalar mass-squared matrix** decomposes into three decoupled blocks:

*X sector:* Both Re(X) and Im(X) have m^2 = 0 (flat direction at tree level).

*Re sector* (basis: Re phi, Re phi_tilde):

    M_Re^2 / m^2 = [[4v^2 + 1 + 2y,  2v],
                     [2v,              1 ]]

    trace = 4v^2 + 2 + 2y,  det = 1 + 2y

*Im sector* (basis: Im phi, Im phi_tilde):

    M_Im^2 / m^2 = [[4v^2 + 1 - 2y,  2v],
                     [2v,              1 ]]

    trace = 4v^2 + 2 - 2y,  det = 1 - 2y

The +/-2y splitting in the (1,1) entries arises from |F_X|^2 = |f + g phi^2|^2, which contributes +2gf to Re(phi)^2 and -2gf to Im(phi)^2.

**Tachyon condition:**  The Im sector determinant is 1 - 2y, which goes negative for y > 1/2.  One scalar becomes tachyonic at the origin when gf > m^2/2.


### Explicit mass eigenvalues

| v | y | Goldstino | mf_- | mf_+ | mb1^2 | mb2^2 | mb3^2 | mb4^2 |
|---|---|-----------|-------|-------|-------|-------|-------|-------|
| 0 | 0.1 | 0 | 1.000 | 1.000 | 0.800 | 1.000 | 1.000 | 1.200 |
| 0 | 0.5 | 0 | 1.000 | 1.000 | 0.000 | 1.000 | 1.000 | 2.000 |
| 0 | 1.0 | 0 | 1.000 | 1.000 | **-1.000** | 1.000 | 1.000 | 3.000 |
| 1 | 0.1 | 0 | 0.4142 | 2.4142 | 0.1414 | 0.2000 | 5.659 | 6.000 |
| 1 | 0.5 | 0 | 0.4142 | 2.4142 | 0.000 | 0.2984 | 5.000 | 6.702 |
| 1 | 1.0 | 0 | 0.4142 | 2.4142 | **-0.236** | 0.3945 | 4.236 | 7.606 |
| sqrt(3) | 0.1 | 0 | 0.2679 | 3.7321 | 0.0582 | 0.0850 | 13.74 | 14.11 |
| sqrt(3) | 0.5 | 0 | 0.2679 | 3.7321 | 0.000 | 0.1345 | 13.00 | 14.87 |
| sqrt(3) | 1.0 | 0 | 0.2679 | 3.7321 | **-0.0828** | 0.1898 | 12.08 | 15.81 |

Boldface entries indicate tachyonic modes.


## 2. Coleman-Weinberg effective potential

The one-loop CW potential is

    V_CW = (1/(64 pi^2)) STr[M^4 (ln |M^2|/mu^2 - 3/2)]

with the supertrace counting each real scalar d.o.f. with weight +1 and each Weyl fermion d.o.f. with weight -2.  We set mu = m.

For tachyonic scalars (m^2 < 0), we use the real part of the CW integrand: (m^2)^2 (ln|m^2| - 3/2), which equals |m^2|^2 (ln|m^2| - 3/2).

**Supertraces:**

- STr[M^2] = 0 for all v, y (consequence of canonical Kahler and renormalizable W).
- STr[M^4] = 8y^2, independent of v.

Since STr[M^4] is v-independent, the -3/2 constant contributes only to V_0(y), and the v-dependence comes entirely from STr[M^4 ln|M^2|].

The Taylor expansion of V_CW around v = 0 has only even powers (all mass eigenvalues are functions of v^2):

    V_CW(v) = V_0(y) + c_2(y) v^2 + c_4(y) v^4 + c_6(y) v^6 + ...

Coefficients computed via mpmath automatic differentiation at 50-digit precision:

| y | c_2 | c_4 | c_6 |
|---|-----|-----|-----|
| 0.01 | +0.002133 | -0.001920 | +0.001951 |
| 0.05 | +0.05339 | -0.04811 | +0.04897 |
| 0.10 | +0.2142 | -0.1939 | +0.1984 |
| 0.20 | +0.8677 | -0.7992 | +0.8369 |
| 0.30 | +1.9975 | -1.9015 | +2.0800 |
| 0.50 | +6.1807 | -7.0281 | +10.049 |
| 0.70 | +13.653 | -14.562 | +9.588 |
| **1.00** | **+23.550** | **-32/3** | -12.633 |
| 1.50 | +35.755 | +1.804 | -21.904 |
| 2.00 | +44.697 | +12.283 | -22.146 |

Note: at y = 1, c_4 = -32/3 exactly (verified to 30 decimal places).

The physical potential in terms of the field X:

    V_phys = V_0_phys + [c_2 g^2 m^2/(64 pi^2)] |X|^2
                       + [c_4 g^4/(64 pi^2)] |X|^4 + ...

The one-loop pseudo-modulus mass is m_X^2 = c_2 g^2 m^2/(64 pi^2), positive for all y > 0.


### Effective Kahler correction

The v^4 term in V_CW corresponds to a one-loop correction to the effective Kahler potential for X.  Writing V = [g^4/(16 pi^2)] C_4 |X|^4, the coefficient is C_4 = c_4/4.

In the weak-SUSY-breaking limit, the analytical structure is:

    c_4(y) = -(96/5) y^2 - (128/7) y^4 + O(y^6)

The leading term -96y^2/5 gives C_4 = -24y^2/5 in units of g^4/(16 pi^2).


## 3. Sign of the quartic coefficient

| y | c_4 | sign |
|---|-----|------|
| 0.01 | -0.00192 | negative |
| 0.10 | -0.194 | negative |
| 0.50 | -7.028 | negative |
| 1.00 | -32/3 = -10.667 | negative |
| 1.42 | ~0 | crossover |
| 1.50 | +1.804 | positive |
| 2.00 | +12.28 | positive |

The quartic coefficient c_4 is **negative** for y < y_cross = 1.4246 and positive above.

**Answer:** Yes, the one-loop correction generates a negative |X|^4 coefficient for moderate SUSY breaking (y < 1.42).  At y = 1 (gf = m^2), the value is

    c_4 = -32/3

In units of g^4/(16 pi^2 m^2), the quartic coefficient of (|X|/m)^4 is:

    C_4 = c_4/4 = -8/3

For the small-y regime: C_4 = -(24/5)y^2 in units of g^4/(16 pi^2).


## 4. Comparison with tree-level non-canonical Kahler

A tree-level Kahler K = |X|^2 + c|X|^4/m^2 with c = -g^2/12 produces the scalar potential

    V_tree = f^2 / (1 + 4c|X|^2/m^2) = f^2 / (1 - v^2/3)

which has a pole at v = sqrt(3) (the O'Raifeartaigh-Koide point).  Expanding:

    V_tree = f^2 [1 + v^2/3 + v^4/9 + ...]

The tree-level quartic coefficient of v^4 is f^2/9 = y^2 m^4/(9g^2) in physical units.

To compare with the one-loop result in the same units (m^4/(64 pi^2)):

| y | c_4 (1-loop) | tree quartic (CW units) | ratio |
|---|---|---|---|
| 0.1 | -0.194 | 0.702 | -0.276 |
| 0.5 | -7.028 | 17.55 | -0.401 |
| 1.0 | -10.667 | 70.18 | -0.152 |

The ratio is negative (opposite signs) and y-dependent.  **The one-loop c_4 is not equal to -g^2/12** in any normalization.

The tree-level Kahler with c = -g^2/12 and the one-loop CW potential produce quartic terms with different functional dependences on the model parameters:
- Tree-level: proportional to f^2/m^4 (i.e., to y^2/g^2)
- One-loop: proportional to g^4/(64 pi^2) times c_4(y)

The effective Kahler coefficient extracted from the quadratic CW term (c_eff = -c_2/(4y^2) in natural units) gives values of order 5-6, far from -1/12 = -0.083.

**Conclusion:** c = -g^2/12 does not arise from one-loop CW dynamics.  The pole structure at v = sqrt(3) is a tree-level geometric feature of the non-canonical Kahler, not a consequence of quantum corrections.


## 5. Numerical benchmarks at g/m = 0.1, 0.5, 1.0

### g/m = 0.1

| y | f/m | Scalar m^2 at v=0 | m_X^2/m^2 (1-loop) | c_2 | c_4 |
|---|-----|---|---|---|---|
| 0.1 | 1.0 | 0.80, 1.00, 1.00, 1.20 | 3.39e-06 | +0.214 | -0.194 |
| 0.5 | 5.0 | 0.00, 1.00, 1.00, 2.00 | 9.78e-05 | +6.181 | -7.028 |
| 1.0 | 10.0 | -1.00, 1.00, 1.00, 3.00 | 3.73e-04 | +23.55 | -10.67 |

### g/m = 0.5

| y | f/m | Scalar m^2 at v=0 | m_X^2/m^2 (1-loop) | c_2 | c_4 |
|---|-----|---|---|---|---|
| 0.1 | 0.2 | 0.80, 1.00, 1.00, 1.20 | 8.48e-05 | +0.214 | -0.194 |
| 0.5 | 1.0 | 0.00, 1.00, 1.00, 2.00 | 2.45e-03 | +6.181 | -7.028 |
| 1.0 | 2.0 | -1.00, 1.00, 1.00, 3.00 | 9.32e-03 | +23.55 | -10.67 |

### g/m = 1.0

| y | f/m | Scalar m^2 at v=0 | m_X^2/m^2 (1-loop) | c_2 | c_4 |
|---|-----|---|---|---|---|
| 0.1 | 0.1 | 0.80, 1.00, 1.00, 1.20 | 3.39e-04 | +0.214 | -0.194 |
| 0.5 | 0.5 | 0.00, 1.00, 1.00, 2.00 | 9.78e-03 | +6.181 | -7.028 |
| 1.0 | 1.0 | -1.00, 1.00, 1.00, 3.00 | 3.73e-02 | +23.55 | -10.67 |

Fermion masses at v = 0 are always 0 (Goldstino), 1, 1 in units of m, independent of y.  The c_2 and c_4 depend only on y, not on g/m separately.  The physical pseudo-modulus mass m_X^2 scales as (g/m)^2 at fixed y.


## Summary of answers

1. The fermion mass matrix is given above with eigenvalues 0, m(v +/- sqrt(v^2+1)).  The scalar mass-squared matrices decompose into Re and Im 2x2 blocks with determinants 1 + 2y and 1 - 2y.

2. The CW expansion gives c_2 and c_4 as tabulated.  At y = 1: c_2 = +23.55, c_4 = -32/3 (exactly).  In the small-y limit: c_4 = -(96/5)y^2 + O(y^4).

3. Yes, c_4 < 0 for y < 1.425.  At y = 1 the value is -32/3, giving C_4 = c_4/4 = -8/3 in units of g^4/(16 pi^2).

4. The one-loop c_4 is NOT equal to -g^2/12.  The tree-level Kahler pole at v = sqrt(3) is not reproduced by one-loop CW dynamics.

5. Mass eigenvalues and CW coefficients are given in the benchmark tables.  All are explicit functions of the single parameter y = gf/m^2.
