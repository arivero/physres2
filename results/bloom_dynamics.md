# Bloom Dynamics: Instanton Landscape and Phase Selection

## Setup

The Koide parametrization writes the signed square roots as:

    z_k(delta) = 1 + sqrt(2) cos(delta + 2 pi k / 3),    k = 0, 1, 2

with masses m_k = M_0 z_k^2. The Koide ratio Q = sum(m) / (sum sqrt(m))^2 = 2/3
is an algebraic identity of this parametrization (for all delta where every z_k > 0).

The seed delta_0 = 3 pi/4 sits on the boundary z_0 = 0. The physical lepton value
is delta ~ 132.73 deg, displaced 2.27 deg from the seed.

The "bloom problem": what potential V(delta) selects the physical phase while
keeping the system on the Koide manifold?


## 1. The cos(9 delta) Landscape

A three-instanton superpotential W_3 = c_3 (det M)^3 / Lambda^18 generates a
potential V ~ -cos(9 delta + phi), since det M ~ e^{3i delta} and cubing gives e^{9i delta}.

This potential has 9 equally spaced minima at delta_n = 2 pi n / 9 (n = 0,...,8),
or equivalently 3 distinct minima within each Koide period [0, 2pi/3):

    delta = 0 deg,  40 deg,  80 deg

The physical lepton delta is at 132.73 deg. The nearest cos(9 delta) minimum
is at n=3, i.e., delta = 120 deg, at a distance of 12.73 deg.

The cos(9 delta) value at the physical point: cos(9 * 132.73 deg) = -0.416.
This is neither a minimum (+1) nor a maximum (-1) -- the physical value does NOT
sit at a cos(9 delta) extremum. The instanton potential alone does not select
the physical delta.


## 2. The Bion Contribution

The bion sum S_bloom(delta) = sum_k omega^k z_k, with omega = exp(2pi i/3), simplifies
analytically using the orthogonality of cube roots of unity:

    sum_k omega^k = 0   (the constant "1" in z_k drops out)

    S_bloom(delta) = sqrt(2) sum_k omega^k cos(delta + 2pi k/3)

Using cos theta = (e^{i theta} + e^{-i theta})/2:

    omega^k e^{i 2pi k/3} = e^{i 4pi k/3},    sum_k = 0
    omega^k e^{-i 2pi k/3} = 1,                sum_k = 3

Therefore:

    S_bloom(delta) = (3 sqrt(2)/2) e^{-i delta}

This is an exact result. The modulus |S_bloom|^2 = 9/2 = 4.5 for ALL delta.

The naive Kahler correction |S_bloom|^2 is a constant and generates no potential
for delta. The bion contribution is trivial for signed square roots.

When z_k is replaced by |z_k| (relevant when some z_k < 0, i.e., near or beyond
the seed), the modulus develops nontrivial delta-dependence with minima near 60 deg
where |S|^2 drops to 1.67. However, the physical leptons sit entirely in the z_k > 0
regime, so this absolute-value structure is irrelevant for the bloom.


## 3. The Combined Potential

Since |S_bloom| is constant, the only delta-dependent bion contribution comes from
the DIFFERENCE term:

    V_bion = A * |S_bloom(delta) - 2 S_bloom(delta_0)|^2
           = A * (9/2) * [5 - 4 cos(delta_0 - delta)]

This is a cosine potential centered at the seed, minimized at delta = delta_0 = 135 deg.

Combined with the instanton term:

    V(delta) = A * (9/2) [5 - 4 cos(delta_0 - delta)] + B cos(9 delta + phi)

The bion term pulls delta toward the seed (135 deg). The instanton term has its
nearest minimum at 120 deg (for phi = 0). The physical equilibrium at 132.7 deg
arises from competition between these two forces.

For the stationarity condition dV/ddelta = 0 at the physical delta, the ratio
B/A depends on the phase phi:

| phi (deg) | B/A     | d^2V/ddelta^2 |
|----------:|--------:|--------------:|
|         0 | -0.087  |          15.1 |
|        30 | -0.137  |           9.0 |
|        90 |  0.190  |          32.0 |
|       250 | -0.988  |          97.7 |
|       270 | -0.190  |          32.0 |

The solution exists for a continuous family of (B/A, phi). The maximum curvature
at the minimum (deepest well) occurs near phi = 250 deg with B/A ~ -1.

Physical picture:
- The bion term provides a restoring force toward the seed
- The instanton term provides 9-fold periodic kicks
- Equilibrium sits between seed (135 deg) and nearest instanton minimum (120 deg)
- The physical delta at 132.7 deg is 2.27 deg from the seed and 12.73 deg from 120 deg
- This implies the bion term dominates: it holds delta close to the seed,
  with a small perturbation from the instanton potential


## 4. Precision Delta Measurement

Using PDG masses: m_e = 0.51099895 MeV, m_mu = 105.6583755 MeV, m_tau = 1776.86 MeV.

Extraction:
- sqrt(M_0) = (sqrt(m_e) + sqrt(m_mu) + sqrt(m_tau)) / 3 = 17.7156 MeV^{1/2}
- M_0 = 313.841 MeV
- delta = 2.316625 rad = 132.7328 deg

Precision comparison:
- delta mod (2pi/3) = 0.222230 rad
- 2/9              = 0.222222 rad
- Residual = 7.41 x 10^{-6} rad

| Quantity | Value |
|---------:|------:|
| delta mod (2pi/3) | 0.2222296315 rad |
| 2/9 | 0.2222222222 rad |
| Residual | 7.41 x 10^{-6} rad |
| Residual / (2/9) | 33.3 ppm |
| Residual / (2pi/3) | 3.5 ppm |
| sigma(delta) from m_tau | 8.35 x 10^{-6} rad |
| Significance | 0.89 sigma |

The physical delta agrees with the rational value 2/9 to 33 ppm (relative to 2/9)
or 3.5 ppm (relative to the period 2pi/3). The deviation is 0.89 sigma given
the m_tau uncertainty of 0.12 MeV.

Note: the Koide Q for the physical lepton masses deviates from 2/3 by -9.2 ppm,
so the masses do not lie exactly on the parametrization manifold. The small
reconstruction errors (10^{-4} in sqrt(m)) reflect this.


## 5. Q = 2/3 Identically on the Koide Manifold

Algebraic proof:

    sum z_k = 3 + sqrt(2) * 0 = 3     (sum of equispaced cosines vanishes)
    sum z_k^2 = 3 + 0 + 2 * 3/2 = 6   (using sum cos^2 = 3/2)

    Q = sum(m) / (sum sqrt(m))^2 = M_0 * 6 / (sqrt(M_0) * 3)^2 = 6/9 = 2/3

This holds for all delta where every z_k > 0. At the boundary (z_k = 0, the seed),
Q is still 2/3. Therefore:

    dQ/ddelta = 0       (identically)
    d^2Q/ddelta^2 = 0   (identically)

A rotation of delta is automatically Q-preserving. The bloom problem is not about
protecting Q -- it is about what dynamics selects the VALUE of delta.

Numerical verification: max |Q - 2/3| = 7.8 x 10^{-16} over the z_k > 0 regime
(machine epsilon).

When some z_k < 0 (beyond the seed), the standard Koide identity breaks because
sqrt(m_k) = |z_k| is no longer equal to z_k. In that regime Q drops below 2/3.


## 6. Implications

The bloom problem decomposes into two distinct questions:

(a) **What holds the system on the Koide manifold?** Since Q = 2/3 is an exact
identity of the parametrization, any dynamics that keeps the system parametrized
by (M_0, delta) automatically preserves Q. The Koide condition is NOT a tuning
condition on the potential -- it is a geometric property of the mass manifold.

(b) **What selects the physical delta?** The three-instanton cos(9 delta) potential
provides a landscape of 9 possible phases, but the physical lepton delta does not
sit at any of the cos(9 delta) minima. A competing term (the bion / seed-centering
potential) is needed.

The bion term |S_bloom - 2 S_bloom(delta_0)|^2, despite the analytical simplification
that makes |S_bloom|^2 constant, provides a nontrivial RELATIVE potential when
anchored to the seed value. This creates a restoring force toward delta_0 = 3pi/4.

The equilibrium between the bion restoring force (toward 135 deg) and the instanton
landscape (nearest minimum at 120 deg) produces a physical delta at 132.7 deg.
The dominance of the bion term (2.27 deg displacement vs. 15 deg to the instanton
minimum) explains why delta sits so close to the seed.

The bloom rotation (seed to physical) is exactly:

    delta_seed mod (2pi/3) - delta_phys mod (2pi/3)
        = pi/12 - 2/9  (if 2/9 is exact)
        = 0.03958 rad = 2.268 deg

This is an algebraic number (difference of pi/12 and a rational). The seed phase
pi/12 is exact (from 3pi/4 mod 2pi/3); the physical phase 2/9 is empirical.

Whether the closeness of delta mod (2pi/3) to 2/9 (33 ppm, 0.89 sigma) is
accidental or has a deeper origin remains an open question.
