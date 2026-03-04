# Canonical Rescaling of the Lagrange Multiplier X

## Setup

In SQCD with N_f = N_c = 3 (flavors u, d, s), the confined-phase superpotential is

    W = sum_i m_i M^i_i + X(det M - BB~ - Lambda^6) + y_c H_u M^d_d + y_b H_d M^s_s + lambda X (H_u.H_d)

The Lagrange multiplier X has [X] = [mass]^{-3}. The NMSSM coupling lambda X (H_u.H_d)
requires lambda to have dimension [mass]^4 for W to have correct dimension [mass]^3.
The problem statement's option (a) says [mass]^2 -- this is incorrect; the correct answer is [mass]^4.

## 1. Rescaled Field

Define S = X Lambda^4, so [S] = [mass]^{-3} [mass]^4 = [mass].

    <S> = X_0 Lambda^4 = (-1.210285e-09) * (8.100000e+09) = -9.803305 MeV

Since X_0 = -C/Lambda^6 and C = Lambda^2 (m_u m_d m_s)^{1/3} = 882297.4257 MeV^2:
    <S> = -C/Lambda^2 = -882297.4257/90000 = -9.803305 MeV

## 2. Superpotential in Terms of S

    W = sum_i m_i M^i_i
        + S (det M / Lambda^4 - BB~/Lambda^4 - Lambda^2)
        + y_c H_u M^d_d + y_b H_d M^s_s
        + 0.72 S (H_u.H_d)

New couplings:
- Constraint: S couples to det M / Lambda^4 (dim [mass]^{-4} = 1.2346e-10)
- Tadpole: S has a linear term -S Lambda^2 = -S * 90000 MeV^2
- NMSSM: lambda_new = 0.72 (dimensionless, since lambda_orig was [mass]^4)

## 3. Effective mu-Term

    mu_eff = lambda_new <S> = 0.72 * (-9.803305) = -7.058379 MeV
           = 0.007058 GeV

This is three orders of magnitude below the electroweak scale (~100 GeV).
The smallness is physical: <S> = -(m_u m_d m_s)^{1/3} / Lambda^2,
which is set by light quark masses over the confinement scale.

## 4. Higgs Quartic and Tree-Level Mass

At tan(beta) = 1, the NMSSM singlet contribution to the Higgs mass is:

    m_h^2 = lambda_new^2 v^2 sin^2(2 beta) / 2 = (0.72)^2 * (246.22 GeV)^2 / 2
    m_h = 125.3548 GeV

The rescaling does NOT change this. lambda_new = 0.72 in both bases
(the original lambda was already [mass]^4 = 0.72 * Lambda^4, giving 0.72 after
absorbing the Lambda^4 into S).

For m_h = 125.25 GeV, we need lambda = sqrt(2) * 125.25/246.22 = 0.7194.
The value 0.72 gives 125.4 GeV at tree level, -0.1 GeV below the observed mass.

## 5. Kahler Potential

Original: K = |X|^2 - |X|^4/(12 mu_K^2) + ...

In terms of S:
    K = |S|^2/Lambda^8 - |S|^4/(12 mu_K^2 Lambda^{16}) + ...

    K_{S Sbar} = 1/Lambda^8 - |S|^2/(3 mu_K^2 Lambda^{16})

The metric is NOT canonical (factor 1/Lambda^8). To get canonical kinetic terms,
one must rescale back to X = S/Lambda^4, recovering the original basis.

The pseudo-modulus pole is at:
    |S_pole| = sqrt(3) mu_K Lambda^4 = 4.2089e+12 MeV
    |X_pole| = sqrt(3) mu_K = 519.6152 MeV  (same physical location)

|<S>|/|S_pole| = 2.3292e-12 (VEV far from pole, consistent with perturbativity).

## 6. F-Terms and Scalar Potential

F_S = F_X / Lambda^4 (chain rule: dW/dS = (dW/dX)(dX/dS) = F_X / Lambda^4)

The physical scalar potential compensates via the inverse Kahler metric:
    V = K^{S Sbar} |F_S|^2 = Lambda^8 |F_X/Lambda^4|^2 = |F_X|^2

Numerical verification: V_X = 3.125132e+40, V_S = 3.125132e+40, ratio = 1.000000000000000

The scalar potential is EXACTLY invariant.

## 7. Summary

The rescaling X -> S = X Lambda^4 is a holomorphic field redefinition.
It changes NO physical observable:

| Quantity | Changed? | Reason |
|----------|----------|--------|
| Scalar potential V | No | K^{IJ} compensates F_I rescaling |
| Physical masses | No | Eigenvalues of V'' are basis-independent |
| Mixing angles | No | Diagonalization of mass matrix is basis-independent |
| F-terms | Yes (by Lambda^4) | But K^{IJ} compensates, leaving V unchanged |
| Kahler metric | Yes (1 -> 1/Lambda^8) | Non-canonical in S basis |
| Coupling lambda | No | 0.72 in both bases (correct dim. analysis) |
| mu_eff | No | 7.058379 MeV = 0.007058 GeV in both bases |
| m_h (tree) | No | 125.3548 GeV in both bases |

### The mu-Problem

The effective mu-term |mu_eff| = 7.058 MeV = 0.007058 GeV is far too small
for viable EWSB (requires ~100 GeV). This is a physical problem independent
of the field basis. It arises because the Seiberg seesaw VEV <X> (or <S>)
is set by light quark masses: <S> = -(m_u m_d m_s)^{1/3}/Lambda^2.

To get mu_eff ~ 100 GeV, one would need either:
- lambda >> 1 (non-perturbative)
- A separate source for the mu-term (Giudice-Masiero mechanism)
- The NMSSM singlet is not identified with the Seiberg Lagrange multiplier
- A different SQCD phase (N_f != N_c = 3)

### Numerical Summary Table

| Parameter | Value |
|-----------|-------|
| Lambda | 300.0 MeV |
| C = Lambda^2 (m_u m_d m_s)^{1/3} | 882297.4257 MeV^2 |
| <X> | -1.210285e-09 MeV^{-3} |
| <S> = <X> Lambda^4 | -9.803305 MeV |
| lambda_new | 0.72 (dimensionless) |
| mu_eff = lambda_new <S> | -7.058379 MeV |
| m_h (tree, NMSSM) | 125.3548 GeV |
| Kahler pole |S_pole| | 4.2089e+12 MeV |
| y_c | 0.0103159776 |
| y_b | 0.0339533750 |

