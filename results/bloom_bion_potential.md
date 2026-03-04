# Bion Potential on the Koide Manifold

## Setup

The Koide parametrization for a mass triple with Q = 2/3:

    z_k(delta) = v0 * (1 + sqrt(2) cos(delta + 2 pi k/3)),    k = 0, 1, 2
    m_k = z_k^2

Identities (algebraic, valid for all delta):
- sum z_k = 3 v0
- sum z_k^2 = 6 v0^2
- Q = sum(m_k)/(sum z_k)^2 = 6 v0^2 / 9 v0^2 = 2/3

The bion Kahler potential on the Koide manifold:

    V_bion(delta) = (zeta^2/Lambda^2) |sum_k s_k(delta) sqrt(m_k(delta))|^2

## Part 1: Masses in terms of (v0, delta)

    m_k(delta) = v0^2 * (1 + sqrt(2) cos(delta + 2 pi k/3))^2

At the seed delta_0 = 3 pi/4:
- z_0 = 0, z_1 = v0 * (1 + sqrt(2) cos(3pi/4 + 2pi/3)) = v0 * 0.634, z_2 = v0 * 2.366
- m_0 = 0 (one massless state)

## Part 2: Bion potential -- the sign assignment is critical

Since sqrt(m_k) = |z_k| and s_k is a sign, we have s_k sqrt(m_k) = s_k |z_k|.

**Case A (adiabatic signs, s_k = sign(z_k)):**

    s_k |z_k| = z_k  for all k
    sum s_k |z_k| = sum z_k = 3 v0
    V_bion = 9 v0^2 = CONSTANT for all delta

This is an exact algebraic identity. The adiabatic bion potential is **completely flat** on the Koide manifold. It cannot select delta.

Numerical verification: max|V - 9| = 1.4e-14 over [0, 2pi].

**Case B (frozen/topological signs, all s_k = +1):**

    V_frozen = (sum |z_k|)^2

When all z_k > 0: V = (sum z_k)^2 = 9 (same as constant).
When some z_k < 0: |z_k| != z_k, so V differs from 9.

    V_frozen range = [9.0, 14.657]

Nontrivial structure appears only when z_k crosses zero.

**Case C (one sign flipped, modeling the (-s,c,b) sector):**

    V_mixed = (-|z_0| + |z_1| + |z_2|)^2

This has zeros (V = 0) where -|z_0| + |z_1| + |z_2| = 0, and equals 9 in the all-positive regime.

## Part 3: Full map of V_bion(delta) -- local minima

Zero-crossings (z_k = 0) occur at 6 points in [0, 2pi]:

| k | delta (deg) |
|---|-------------|
| 0 | 135, 225    |
| 1 | 15, 105     |
| 2 | 255, 345    |

For 75% of the delta range, at least one z_k is negative.

V_frozen = 9 (flat plateau) in the region where all z_k > 0. This covers:
- [105 deg, 135 deg] (z_0 and z_1 both near zero, z_2 large)
- [225 deg, 255 deg] (z_0 and z_2 both near zero, z_1 large)
- [345 deg, 15 deg] (z_1 and z_2 both near zero, z_0 large)

Each plateau spans ~30 degrees. Between plateaus, V_frozen rises to max = 14.657 at delta = 60, 180, 300 deg (where one z_k hits its most negative value, -0.414).

The "minima" of V_frozen are not isolated points but **flat plateaus** at V = 9.

### Analytical structure near the seed

For delta < 3pi/4 (z_0 > 0): V_frozen = 9 (constant, flat)
For delta > 3pi/4 (z_0 < 0): V_frozen = (3 - 2z_0)^2 > 9, with z_0 = 1 + sqrt(2) cos(delta) < 0

The seed is a **corner minimum** -- flat on the left, rising on the right. The derivative is discontinuous at the seed.

    V(seed)            = 9.0000
    V(seed + 0.01 rad) = 9.1198
    V(seed - 0.01 rad) = 9.0000

## Part 4: Quark triple (-s, c, b)

Input masses (MeV): m_s = 93.4, m_c = 1270, m_b = 4180

Signed roots: z = (-9.664, 35.637, 64.653)

    v0(full) = (-sqrt(m_s) + sqrt(m_c) + sqrt(m_b))/3 = 30.209 MeV^{1/2}
    delta_phys = 157.21 deg (= 2.744 rad)
    Q(-s,c,b) = 0.6750  (deviates from 2/3 by 1.24%; quarks are near but not on the Koide manifold)

**v0-doubling confirmed:**

    v0(seed) = (sqrt(m_s) + sqrt(m_c))/3 = 15.100 MeV^{1/2}
    v0(full)/v0(seed) = 2.0005

**Bion potential at the physical quark delta:**

The physical delta = 157.21 deg is in the z_0 < 0 region (between 135 deg and 225 deg).

    Normalized z at physical delta: [-0.304, 1.178, 2.126]
    V_frozen(physical) = 13.015 > 9
    V_frozen is monotonically increasing past the seed
    => The bion potential DISFAVORS the bloom position
    => It acts as a RESTORING FORCE toward the seed

The physical quark delta is NOT at a minimum of V_bion. The frozen-sign bion potential pushes delta BACK toward the seed at 135 deg.

## Part 5: Lepton triple (e, mu, tau)

Input masses (MeV): m_e = 0.511, m_mu = 105.66, m_tau = 1776.86

Signed roots: z = (0.715, 10.279, 42.153)

    v0 = 17.716 MeV^{1/2}
    delta_phys = 132.73 deg (= 2.317 rad)

The physical lepton delta is at 132.73 deg, BEFORE the seed at 135 deg. All z_k > 0.

    V_frozen(lepton) = 9.0000 exactly

The bion potential is completely flat at the lepton position. It provides no delta-selection for leptons.

The lepton delta is notable for:

    delta mod (2pi/3) = 0.22223 rad
    2/9               = 0.22222 rad
    Residual          = 7.4e-6 rad  (33 ppm of 2/9)

## Part 6: What drives the bloom?

### Answer

The bion potential V_bion(delta) does **not** have a minimum at the physical quark bloom position.

**(a) Adiabatic signs:** V = 9 v0^2 = constant everywhere. No delta-selection at all. This is an exact algebraic identity of the Koide parametrization: the sum of signed square roots is fixed.

**(b) Frozen signs:** V has its minimum (a flat plateau at V = 9) at and before the seed (delta <= 135 deg), and rises monotonically at the physical quark delta (~157 deg) where V = 13.015. The bion potential acts as a restoring force toward the seed, **opposing** the bloom.

### The three-instanton alternative

The three-instanton potential V_3inst = -cos(9 delta) has its nearest minimum to the physical quark delta at:

    delta = 160 deg  (n=4 minimum)
    Distance from physical quark delta: 2.8 deg

This is strikingly close. The physical quark delta sits 2.8 degrees from a cos(9delta) minimum, but 22.2 degrees from the seed.

### Combined picture

V_total(delta) = A * (V_frozen - 9) + B * (-cos(9 delta))

where A is the bion coupling (restoring toward seed) and B is the three-instanton coupling (driving toward delta = 160 deg).

Required B/A for equilibrium at the physical quark delta:

    B/A = -2.07

The second derivative at the physical delta is positive (d^2V/ddelta^2 = 135.4), confirming a stable minimum.

Physical picture:
1. The three-instanton potential provides the attractive force that drives delta away from the seed (135 deg) toward the cos(9delta) minimum at 160 deg
2. The bion potential provides a restoring force keeping delta from overshooting
3. The equilibrium is the physical bloom position at ~157 deg

### Consistency with known results

- The bloom is **nonperturbative**, consistent with the earlier finding that CW masses are 5 orders of magnitude too small
- The driving force is multi-instanton (cos(9delta)), not semiclassical (bion)
- The bion potential does not select the bloom -- it **stabilizes** it
- For leptons (delta = 132.73 deg < seed), the bion potential is flat; the lepton delta is selected by the cos(9delta) minimum at 120 deg (distance 12.7 deg), with no bion restoring force acting (since all z_k > 0)

## Key Algebraic Result

The flatness of the adiabatic bion potential is the most important result. On the Koide manifold:

    sum_k sign(z_k) |z_k| = sum_k z_k = 3 v0

This holds because sign(z_k) |z_k| = z_k identically. The Koide constraint forces the sum of signed square roots to be constant, regardless of delta. Any bion-type potential built from sum s_k sqrt(m_k) with adiabatic sign-tracking is automatically constant on the Koide manifold.

The only escape is non-adiabatic sign dynamics (topological sector freezing), which produces a restoring force toward the seed, not a bloom-driving force.
