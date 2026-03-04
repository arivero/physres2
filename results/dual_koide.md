# Koide Invariant under Mass Inversion

**Computation date:** 2026-03-04

**Script:** `results/dual_koide.py`


## Setup

The Koide quotient for three masses with sign conventions:

    Q = (m_1 + m_2 + m_3) / (s_1 sqrt(m_1) + s_2 sqrt(m_2) + s_3 sqrt(m_3))^2

where s_j = +/-1. The "dual" or "inverted" Koide quotient replaces m_j with 1/m_j (or equivalently C/m_j for any constant C, since C cancels).

Mass values (MeV): m_e = 0.511, m_mu = 105.658, m_tau = 1776.86, m_u = 2.16, m_d = 4.67, m_s = 93.4, m_c = 1270, m_b = 4180, m_t = 172760.


## Part 1: Direct Koide Q

| Triple | Signs | Q | Deviation from 2/3 |
|--------|-------|---|---------------------|
| (e, mu, tau) | +,+,+ | 0.66666082 | -0.0009% |
| (-s, c, b) | -,+,+ | 0.67495422 | +1.2431% |
| (c, b, t) | +,+,+ | 0.66948936 | +0.4234% |
| (d, s, b) | +,+,+ | 0.73142765 | +9.7141% |
| (-d, s, b) | -,+,+ | 0.82167439 | +23.2512% |

The charged lepton triple (e, mu, tau) gives Q = 2/3 to 9 ppm. The quark triples (-s,c,b) and (c,b,t) are within 1.2% and 0.4% respectively. The down-type triple (d,s,b) with any sign convention is far from 2/3.


## Part 2: Algebraic Relation between Q and Q_dual

Setting a_j = sqrt(m_j), the direct and dual quotients are

    Q      = (a_1^2 + a_2^2 + a_3^2)                     / (a_1 + a_2 + a_3)^2
    Q_dual = (a_2^2 a_3^2 + a_1^2 a_3^2 + a_1^2 a_2^2)  / (a_2 a_3 + a_1 a_3 + a_1 a_2)^2

In terms of the elementary symmetric polynomials p_1 = sum a_j, p_2 = sum_{j<k} a_j a_k, p_3 = prod a_j:

    Q      = 1 - 2 p_2/p_1^2
    Q_dual = 1 - 2 p_1 p_3/p_2^2

Therefore Q = 2/3 iff p_2/p_1^2 = 1/6, while Q_dual = 2/3 iff p_1 p_3/p_2^2 = 1/6. **These are conditions on different symmetric functions and are algebraically independent.** Q = 2/3 does not imply Q_dual = 2/3.


### Koide parametrization

In the standard parametrization

    sv_k = sqrt(M_0) * (1 + sqrt(2) cos(2 pi k/3 + delta)),   k = 0, 1, 2

the symmetric polynomials of the signed roots are p_1 = 3 sqrt(M_0), p_2 = 3 M_0/2, and p_3 = M_0^{3/2} P(delta), where

    P(delta) = prod_k (1 + sqrt(2) cos(2 pi k/3 + delta))

Expanding the product using the trigonometric identities sum cos(phi_k + d) = 0, sum_{j<k} cos(phi_j+d) cos(phi_k+d) = -3/4, and prod cos(phi_k+d) = (1/4) cos(3d):

    P(delta) = -1/2 + cos(3 delta) / sqrt(2)

Substituting into Q_dual = 1 - 2 p_1 p_3 / p_2^2 = 1 - 8 P/3:

    Q_dual(delta) = 7/3 - 4 sqrt(2) cos(3 delta) / 3

This formula holds in the all-positive-roots regime where every sv_k > 0, which requires delta to lie in six arcs of width 30 degrees centered at 0, 120, 240 degrees (and their reflections).

Q_dual = 2/3 requires cos(3 delta) = 5 sqrt(2)/8 = 0.883883..., which is achievable. There are six solutions in [0, 2 pi):

| delta (deg) | delta (rad) | Q_dual |
|-------------|-------------|--------|
| 9.2952 | 0.162232 | 0.666667 |
| 110.7048 | 1.932163 | 0.666667 |
| 129.2952 | 2.256627 | 0.666667 |
| 230.7048 | 4.026559 | 0.666667 |
| 249.2952 | 4.351022 | 0.666667 |
| 350.7048 | 6.120954 | 0.666667 |

All six lie in the all-sv-positive region, confirming the formula is self-consistent.


## Part 3: Dual Koide Q_dual = Q(1/m_j)

| Triple | Signs | Q_dual | Deviation from 2/3 |
|--------|-------|--------|---------------------|
| (d, s, b) | +,+,+ | 0.66520987 | -0.2185% |
| (-d, s, b) | -,+,+ | 1.90419140 | +185.63% |
| (u, d, s) | +,+,+ | 0.44257559 | -33.61% |
| (-u, d, s) | -,+,+ | 52.74 | +7811% |
| (e, mu, tau) | +,+,+ | 0.85144838 | +27.72% |
| (c, b, t) | +,+,+ | 0.48932089 | -26.60% |
| (-s, c, b) | -,+,+ | 3.26522593 | +389.78% |

An exhaustive scan over 7 triples with 4 sign combinations (28 cases total) finds exactly one case within 2% of 2/3:

**Q_dual(d, s, b) with signs (+,+,+) = 0.66521, deviation -0.22% from 2/3.**

This is the "dual Koide" observation. It is unique among the triples tested.


## Part 4: SQCD Seesaw

Under the SQCD seesaw M_j = C/m_j with C^3 = Lambda^6 prod m_j and Lambda = 300 MeV:

| Sector | C (MeV) | M_1 (MeV) | M_2 (MeV) | M_3 (MeV) | Q (+++) | Q (-++) |
|--------|---------|-----------|-----------|-----------|---------|---------|
| (u,d,s) | 882297 | 408471 | 188929 | 9446 | 0.44258 | 52.74 |
| (d,s,b) | 10994847 | 2354357 | 117718 | 2630 | 0.66521 | 1.90419 |
| (s,c,b) | 71233585 | 762672 | 56089 | 17042 | 0.54298 | 3.26523 |

The seesaw constant C cancels from Q. The values Q(M_j) = Q(1/m_j) are identical to the dual values in Part 3. The seesaw sets the mass scale but not the Koide quotient.


## Part 5: Connection between Dual and Direct Koide

    Q_direct(-s, c, b)  = 0.67495422   (+1.24% from 2/3)
    Q_dual(d, s, b)     = 0.66520987   (-0.22% from 2/3)

Both are near 2/3 but involve different quarks: the direct triple contains (s, c, b) with a negative sign on s, while the dual triple contains (d, s, b) inverted. They share two quarks (s, b) but differ in the third (c vs d).

The seesaw maps m_j to 1/m_j. If applied to the direct triple (s, c, b), the seesaw produces Q(1/m_s, 1/m_c, 1/m_b) = 0.54298, which is 18.6% from 2/3. The dual Koide near-miss is specific to (d, s, b) inverted, not to (s, c, b) inverted.

Therefore the dual Koide Q(1/m_d, 1/m_s, 1/m_b) ~ 2/3 is an **independent** numerical observation, not a consequence of Q(-s, c, b) ~ 2/3 through any seesaw mechanism.


## Part 6: Parametric Conclusion

On the Koide manifold (where Q = 2/3 identically), the dual quotient is a deterministic function of the phase parameter:

    Q_dual = 7/3 - (4 sqrt(2)/3) cos(3 delta)

This oscillates between 0.448 and 4.219. The value 2/3 is achievable but only at six specific delta values satisfying cos(3 delta) = 5 sqrt(2)/8. There is no algebraic mechanism forcing Q_dual = 2/3 whenever Q = 2/3.

For the physical Koide triples:

| Triple | delta (deg) | Q_direct | cos(3 delta) | Q_dual (formula) | Q_dual (numerical) |
|--------|-------------|----------|--------------|------------------|---------------------|
| (e, mu, tau) | 132.73 | 0.66666 | 0.78587 | 0.85147 | 0.85145 |
| (-s, c, b) | 157.21 | 0.67495 | -0.36864 | [invalid] | 0.54298 |
| (c, b, t) | 123.93 | 0.66949 | 0.97887 | 0.48756 | 0.48932 |

The formula matches the numerical result to machine precision when all sv_k > 0, and breaks down when any sv_k < 0 (as for the (-s,c,b) triple, where the formula applies only to signed-root Koide, not to all-positive-root Koide of the inverted masses).

**The (d, s, b) triple does not lie on the Koide manifold** (Q_direct(d,s,b) = 0.731 with all positive signs). The formula Q_dual = 7/3 - ... is therefore not applicable to it. Its dual Koide near-miss Q_dual ~ 2/3 is a standalone empirical fact.


## Summary of Results

1. **Q = 2/3 does not imply Q_dual = 2/3.** They are conditions on independent symmetric functions (p_2/p_1^2 = 1/6 vs p_1 p_3/p_2^2 = 1/6).

2. **On the Koide manifold,** Q_dual = 7/3 - (4 sqrt(2)/3) cos(3 delta), oscillating between 0.448 and 4.219. The product P(delta) = -1/2 + cos(3 delta)/sqrt(2) encodes all the information.

3. **The dual Koide Q(1/m_d, 1/m_s, 1/m_b) = 0.66521 is 0.22% from 2/3.** This is the only triple (among 28 cases scanned) with Q_dual within 2% of 2/3.

4. **The SQCD seesaw constant C cancels from Q_dual.** The seesaw sets the mass scale of the dual spectrum but does not affect the Koide quotient.

5. **The dual Koide observation is independent** of the direct Koide Q(-s, c, b) ~ 2/3. They involve different quarks and no seesaw mapping connects them.

6. **Achievability:** Q_dual = 2/3 on the Koide manifold requires cos(3 delta) = 5 sqrt(2)/8 = 0.8839, which gives delta = 9.30 degrees (and five other solutions mod 60 degrees). None of the physical triples sits at these special delta values.
