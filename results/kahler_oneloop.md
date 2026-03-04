# One-loop effective Kahler potential in the O'Raifeartaigh model

## Model and conventions

The O'Raifeartaigh model has superpotential

W = f X + m phi phi_tilde + g X phi^2

with chiral superfields X, phi, phi_tilde and real parameters f, m, g > 0.  The tree-level scalar potential V = |dW/dX|^2 + |dW/dphi|^2 + |dW/dphi_tilde|^2 breaks SUSY via F_X = f.  The field X is a pseudo-modulus: classically flat at phi = phi_tilde = 0.

Dimensionless variables: v = g <X>/m (pseudo-modulus VEV), y = gf/m^2 (SUSY-breaking parameter).  All masses are given in units of m throughout.

## 1. Mass matrices

The fermion mass matrix in the basis (psi_X, psi_phi, psi_phi_tilde) at phi = phi_tilde = 0 is

M_F / m = [[0, 0, 0], [0, 2v, 1], [0, 1, 0]]

The X-fermion (Goldstino) decouples with mass zero.  The 2x2 block has eigenvalues

m_+/m = v + sqrt(v^2 + 1),    m_-/m = sqrt(v^2 + 1) - v

with product m_+ m_- = m^2 and sum m_+ + m_- = 2m sqrt(v^2+1).

The scalar mass-squared matrix decomposes into three decoupled blocks:

**X sector:** Both Re(X) and Im(X) have m^2 = 0 (flat direction).

**Re sector** (basis: Re phi, Re phi_tilde):

M_Re^2 / m^2 = [[4v^2 + 1 + 2y,  2v], [2v, 1]]

trace = 4v^2 + 2 + 2y,  det = 1 + 2y.

**Im sector** (basis: Im phi, Im phi_tilde):

M_Im^2 / m^2 = [[4v^2 + 1 - 2y,  2v], [2v, 1]]

trace = 4v^2 + 2 - 2y,  det = 1 - 2y.

The +/-2y splitting in the (1,1) entries arises from expanding |F_X|^2 = |f + g phi^2|^2, which gives a mass contribution +2gf to Re(phi)^2 and -2gf to Im(phi)^2.

**Tachyon condition:**  The Im sector has det = 1 - 2y < 0 when y > 1/2, signaling that the origin phi = 0 becomes an unstable saddle point for large SUSY breaking.

### Explicit mass eigenvalues

| v | y | Goldstino | mf_- | mf_+ | mb1^2 | mb2^2 | mb3^2 | mb4^2 |
|---|---|-----------|-------|-------|-------|-------|-------|-------|
| 0 | 0.1 | 0 | 1.000 | 1.000 | 0.800 | 1.000 | 1.000 | 1.200 |
| 0 | 0.5 | 0 | 1.000 | 1.000 | 0.000 | 1.000 | 1.000 | 2.000 |
| 0 | 1.0 | 0 | 1.000 | 1.000 | -1.000 | 1.000 | 1.000 | 3.000 |
| 1 | 0.1 | 0 | 0.4142 | 2.4142 | 0.1414 | 0.2000 | 5.659 | 6.000 |
| 1 | 0.5 | 0 | 0.4142 | 2.4142 | 0.000 | 0.2984 | 5.000 | 6.702 |
| 1 | 1.0 | 0 | 0.4142 | 2.4142 | -0.236 | 0.3945 | 4.236 | 7.606 |
| sqrt(3) | 0.1 | 0 | 0.2679 | 3.7321 | 0.0582 | 0.0850 | 13.74 | 14.11 |
| sqrt(3) | 0.5 | 0 | 0.2679 | 3.7321 | 0.000 | 0.1345 | 13.00 | 14.87 |
| sqrt(3) | 1.0 | 0 | 0.2679 | 3.7321 | -0.0828 | 0.1898 | 12.08 | 15.81 |

Entries marked with negative mb^2 indicate tachyonic modes (at y = 1.0).


## 2. Coleman-Weinberg effective potential

The one-loop CW potential is

V_CW = (1/(64 pi^2)) STr[M^4 (ln M^2/mu^2 - 3/2)]

with the supertrace counting each real scalar d.o.f. with weight +1 and each Weyl fermion d.o.f. with weight -2.  We take the renormalization scale mu = m, so that ln(M^2/mu^2) = ln(M^2/m^2) is the log of the dimensionless mass-squared.

**Supertraces:**
- STr[M^2] = 0 for all v, y (verified numerically; consequence of canonical Kahler and renormalizable W).
- STr[M^4] = 8y^2, independent of v.  This means the -3/2 constant in the CW integrand contributes only to V_0, not to the v-dependence of V_CW.

The Taylor expansion of V_CW around v = 0 has only even powers of v (since all mass eigenvalues depend on v only through v^2):

V_CW(v) = V_0(y) + c_2(y) v^2 + c_4(y) v^4 + c_6(y) v^6 + ...

The coefficients were computed at 50-digit precision using mpmath's automatic differentiation:

| y | c_2 | c_4 | c_6 |
|---|-----|-----|-----|
| 0.01 | +0.00213342 | -0.00192018 | +0.00195080 |
| 0.05 | +0.05338682 | -0.04811473 | +0.04896612 |
| 0.10 | +0.21419659 | -0.19385759 | +0.19836618 |
| 0.20 | +0.86765617 | -0.79923491 | +0.83687500 |
| 0.50 | +6.18070978 | -7.02808622 | +10.0493796 |
| 1.00 | +19.5500424 | **-6.66666667** | -9.96658188 |
| 1.50 | +32.4818927 | -1.67069050 | -13.7157955 |
| 2.00 | +46.4718956 | +1.95130724 | -14.4243464 |

The physical effective potential in terms of the field X is:

V_phys = (m^4/(64 pi^2)) [V_0 + c_2 (g|X|/m)^2 + c_4 (g|X|/m)^4 + ...]

which gives:

V_phys = V_0_phys + [c_2 g^2 m^2/(64 pi^2)] |X|^2 + [c_4 g^4/(64 pi^2)] |X|^4 + ...

The one-loop mass squared for the pseudo-modulus X is:

m_X^2 = c_2 (g/m)^2 m^2 / (64 pi^2)

which is positive (c_2 > 0 for all y > 0), confirming that the origin is a one-loop minimum for X.


## 3. Sign of the quartic coefficient

The quartic coefficient c_4 is **negative** for y < y_cross = 1.7126 and positive for y > y_cross.

For moderate and weak SUSY breaking (gf/m^2 < 1.71), the one-loop correction generates a negative |X|^4 term.  This means the effective Kahler potential receives a negative quartic correction, consistent with the known result that the CW potential flattens the pseudo-modulus direction at large field values.

At y = 1 (gf = m^2), the quartic coefficient takes the exact value

c_4(y=1) = -20/3

verified to 30 decimal places.  This gives a physical quartic coefficient:

V_phys contains [g^4/(64 pi^2)] * (-20/3) * |X|^4 = -[5g^4/(48 pi^2)] |X|^4

Or equivalently, in units of g^4/(16 pi^2):

C_4 = c_4/4 = -5/3

In the small-y limit, c_4 approaches a linear behavior in y^2:

c_4 -> -19.2 y^2   as   y -> 0

The ratio c_4/y^2 is not exactly constant (it varies from -19.2 at small y to about -6.7 at y = 1), reflecting genuine higher-order dependence on the SUSY-breaking parameter.


## 4. Comparison with tree-level non-canonical Kahler

The tree-level Kahler potential K = |X|^2 + c|X|^4/m^2 with c < 0 gives

V_tree = f^2 / (1 + 4c|X|^2/m^2)

which has a pole at |X|_pole = m/(2 sqrt(|c|)).  In terms of v = g|X|/m, the pole is at v_pole = g/(2 sqrt(|c|)).  Requiring v_pole = sqrt(3) gives

c = -g^2/12

Expanding around the origin:

V_tree = f^2 [1 + v^2/3 + v^4/9 + ...]

The tree-level quartic coefficient is f^2/9 = y^2 m^4/(9g^2) (in physical units).

The one-loop quartic is c_4 m^4/(64 pi^2) (as a coefficient of v^4).  These are:

| y | c_4 (1-loop, CW units) | f^2/9 (tree, CW units) | ratio |
|---|---|---|---|
| 0.1 | -0.194 | 0.702 | -0.276 |
| 0.5 | -7.028 | 17.55 | -0.401 |
| 1.0 | -6.667 | 70.18 | -0.095 |

The ratio is negative (opposite signs) and varies with y.  The one-loop c_4 is not equal to -g^2/12 (or -1/12, or any y-independent number) in any normalization.  The tree-level Kahler coefficient c = -g^2/12 involves f^2 (through the potential V = f^2/K_{XX*}), making it y-dependent, while the one-loop CW potential has a different functional dependence on y.

The effective Kahler coefficient extracted from the quadratic term (matching c_2 to -4c_eff f^2/m^2) gives c_eff = -c_2/(4y^2), which equals -5.35, -6.18, -4.89 for y = 0.1, 0.5, 1.0 respectively -- all much larger in magnitude than -1/12 = -0.083 and y-dependent.

**Conclusion:** The one-loop Coleman-Weinberg computation does not reproduce c = -g^2/12.  The non-canonical Kahler pole at v = sqrt(3) is a tree-level statement about the cutoff structure, not a consequence of one-loop dynamics.


## 5. Benchmarks at g/m = 0.1, 0.5, 1.0

### g/m = 0.1

| y | f/m | m_X^2/m^2 (1-loop) | c_2 | c_4 | tachyon? |
|---|-----|-----|-----|-----|----------|
| 0.1 | 1.0 | 3.39e-06 | +0.214 | -0.194 | no |
| 0.5 | 5.0 | 9.78e-05 | +6.181 | -7.028 | no |
| 1.0 | 10.0 | 3.10e-04 | +19.55 | -6.667 | yes |

### g/m = 0.5

| y | f/m | m_X^2/m^2 (1-loop) | c_2 | c_4 | tachyon? |
|---|-----|-----|-----|-----|----------|
| 0.1 | 0.2 | 8.48e-05 | +0.214 | -0.194 | no |
| 0.5 | 1.0 | 2.45e-03 | +6.181 | -7.028 | no |
| 1.0 | 2.0 | 7.74e-03 | +19.55 | -6.667 | yes |

### g/m = 1.0

| y | f/m | m_X^2/m^2 (1-loop) | c_2 | c_4 | tachyon? |
|---|-----|-----|-----|-----|----------|
| 0.1 | 0.1 | 3.39e-04 | +0.214 | -0.194 | no |
| 0.5 | 0.5 | 9.78e-03 | +6.181 | -7.028 | no |
| 1.0 | 1.0 | 3.10e-02 | +19.55 | -6.667 | yes |

The c_2 and c_4 values depend on y alone (not on g/m separately).  The physical X mass scales as (g/m)^2 at fixed y.  Tachyonic scalars appear at the origin when y > 1/2 (Im-sector det = 1 - 2y < 0).


## Summary of answers

1. The fermion mass matrix is M_F/m = [[0,0,0],[0,2v,1],[0,1,0]] with eigenvalues 0, m(v +/- sqrt(v^2+1)).  The scalar mass-squared matrix decomposes into Re and Im 2x2 blocks with determinants 1+2y and 1-2y respectively.

2. The CW potential expanded as V_0 + c_2 v^2 + c_4 v^4 gives c_2 = +19.55 and c_4 = -20/3 at y = 1 (in units of m^4/(64pi^2)).  The Kahler coefficient C_4 = c_4/4 = -5/3 in units of g^4/(16pi^2).

3. Yes, c_4 is negative for y < 1.713.  At y = 1, the value is exactly c_4 = -20/3, giving C_4 = -5/3 in units of g^4/(16pi^2).

4. The one-loop c_4 is NOT equal to -g^2/12 in any normalization.  The tree-level Kahler c = -g^2/12 and the one-loop CW quartic have different functional dependences on the model parameters.

5. All mass eigenvalues are tabulated above for g/m = 0.1, 0.5, 1.0 at y = 0.1, 0.5, 1.0.
