# Koide Bloom Delta Analysis

## Setup

The Koide parametrization writes three signed square roots as

    sigma_i = sqrt(M0) * (1 + sqrt(2) * cos(2*pi*(i-1)/3 + delta))

for i = 1, 2, 3, where sigma_i^2 = m_i and the sign convention is
per-triple (all positive for (e,mu,tau) and (c,b,t); negative first
entry for (-s,c,b)).

From the parametrization, since sum_i cos(phi_i + delta) = 0:

    sqrt(M0) = sum(sigma_i) / 3
    delta    = arctan2(-sum(sigma_i sin phi_i), sum(sigma_i cos phi_i))

These are exact when Q = 2/3; for Q != 2/3 they minimize the squared
residual between the parametrization and the data.

When delta = delta_seed = 3*pi/4, the first cosine factor is
cos(3*pi/4) = -1/sqrt(2), so sigma_1 = sqrt(M0) * (1 + sqrt(2) *
(-1/sqrt(2))) = 0. The seed always has its first mass zero.

Masses used (MeV): m_e = 0.51100, m_mu = 105.658, m_tau = 1776.86,
m_s = 93.4, m_c = 1270.0, m_b = 4180.0, m_t = 172760.0.


## Seed geometry

At delta_seed = 3*pi/4 = 2.356194 rad = 135.0 deg:

    c1 = 0  (exact)
    c2 = 1 + sqrt(2)*cos(2*pi/3 + 3*pi/4) = 0.63397460
    c3 = 1 + sqrt(2)*cos(4*pi/3 + 3*pi/4) = 2.36602540
    c3/c2 = (1 + sqrt(3)) / (sqrt(3) - 1) = 3.73205081  (exact)

So at the seed, the two non-zero masses are in the ratio c3^2/c2^2 = 13.928.


## Part 1: Full parametrization

### (e, mu, tau)

    M0_full    = 313.840920 MeV
    delta_full = 2.31662426 rad = 132.732793 deg
    delta_full / pi      = 0.73740440
    delta_full / (2pi/3) = 1.10610660
    Koide Q = 0.66666082  (deviation -0.0009% from 2/3)

### (-s, c, b)

    M0_full    = 912.555714 MeV
    delta_full = 2.74384326 rad = 157.210639 deg
    delta_full / pi      = 0.87339244
    delta_full / (2pi/3) = 1.31008865
    Koide Q = 0.67495422  (deviation +1.2431% from 2/3)

The (-s,c,b) triple does not satisfy Q = 2/3 at the PDG central values;
the parametrization values are the best-fit projection.

### (c, b, t)

    M0_full    = 29576.439062 MeV
    delta_full = 2.16303821 rad = 123.932960 deg
    delta_full / pi      = 0.68851645
    delta_full / (2pi/3) = 1.03277467
    Koide Q = 0.66948936  (deviation +0.4234% from 2/3)


## Part 2: Seed parameters

The Koide seed has sigma_1 = 0 at delta_seed = 3*pi/4. Given the two
non-zero signed roots of the full triple, the seed M0 is fit by
least squares: sqrt(M0_seed) = (sigma_2 * c2 + sigma_3 * c3) / (c2^2 + c3^2).

### (e, mu, tau)

    M0_seed = 313.592543 MeV
    Seed masses: (0, 126.04, 1755.52) MeV
    Actual non-zero masses: m_mu = 105.658, m_tau = 1776.86 MeV

Because m_e << m_mu, m_tau, the full triple is already close to the seed.
The seed masses differ from (m_mu, m_tau) by 19% and -1.2% respectively,
because the c3/c2 ratio at the seed (3.732) does not match m_tau/m_mu^(1/2)
at the actual lepton values.

### (-s, c, b)

    M0_seed = 856.181155 MeV
    Seed masses: (0, 344.12, 4792.97) MeV
    Actual non-zero masses: m_c = 1270.0, m_b = 4180.0 MeV

The Koide seed masses (344, 4793) are far from (m_c, m_b) because the
full triple has a large bloom rotation (22.2 deg, see Part 4) and Q
deviates 1.24% from 2/3.

### (c, b, t)

    M0_seed = 29150.597143 MeV
    Seed masses: (0, 11716.32, 163187.26) MeV
    Actual non-zero masses: m_b = 4180.0, m_t = 172760.0 MeV

Similarly, the seed masses differ substantially from (m_b, m_t) due to
the 11.1 deg bloom rotation.


## Part 3: v0-doubling

Here v0 = sqrt(M0) in the Koide parametrization, since sum(sigma_i)/3
= sqrt(M0) by construction.

| Triple      | M0_full    | M0_seed    | M0_full/M0_seed | v0_full/v0_seed |
|-------------|------------|------------|-----------------|-----------------|
| (e,mu,tau)  | 313.8409   | 313.5925   | 1.000792        | 1.000396        |
| (-s,c,b)    | 912.5557   | 856.1812   | 1.065844        | 1.032397        |
| (c,b,t)     | 29576.4391 | 29150.5971 | 1.014608        | 1.007278        |

None of these Koide-seed ratios is near 2. The M0_full/M0_seed ratios
are all close to 1, reflecting how much the lightest member contributes
to M0 relative to the seed (which uses only the two heavier members).

The v0-doubling factor of 2 from v0_doubling.md comes from a different
seed definition, detailed below.


## v0-doubling: bloom seed for (-s, c, b)

The Koide seed just computed (sigma_1 = 0) is not the seed relevant to
the v0-doubling observation. The relevant seed is the PRECEDING 2-mass
Koide triplet (s, c, 0), which has sigma_3 = 0 and both non-zero entries
with positive sign. This is the mass pair before the b quark bloomed in.

When sigma_3 = 0: cos(4*pi/3 + delta_seed3) = -1/sqrt(2), giving

    delta_seed3 = 3*pi/4 - 4*pi/3  (mod 2*pi)  =  4.450590 rad  =  1.41667 pi

At that seed, the v0 is

    v0_seed_bloom = (sqrt(m_s) + sqrt(m_c)) / 3  =  15.100476 MeV^(1/2)
    M0_seed_bloom = v0_seed_bloom^2               =  228.024364 MeV

The full (-s,c,b) triple has

    v0_full = (-sqrt(m_s) + sqrt(m_c) + sqrt(m_b)) / 3  =  30.208537 MeV^(1/2)
    M0_full = v0_full^2                                  =  912.555714 MeV

Ratios:

    v0_full / v0_seed_bloom = 2.000502  (near 2)
    M0_full / M0_seed_bloom = 4.002010  (near 4)

The bloom rotation relative to this seed is

    Delta_delta_bloom = delta_full - delta_seed3  =  -1.70675 rad  =  -97.789 deg

The v0-doubling implies the exact relation sqrt(m_b) = 3*sqrt(m_s) + sqrt(m_c),
which predicts

    m_b_predicted = (3*sqrt(93.4) + sqrt(1270))^2  =  4177.06 MeV

compared to PDG m_b = 4180 +/- 30 MeV (delta = -0.070%, well within uncertainty).


## Part 4: Bloom rotation

Bloom rotation Delta_delta = delta_full - delta_seed, where delta_seed = 3*pi/4
(sigma_1 = 0 convention), taken in the range (-pi, pi].

| Triple      | delta_full (rad) | delta_full (deg) | Delta_delta (rad) | Delta_delta (deg) | Ddelta/(2pi/3)  |
|-------------|------------------|------------------|-------------------|-------------------|-----------------|
| (e,mu,tau)  | 2.31662426       | 132.732793       | -0.03957023       | -2.267207         | -0.01889340     |
| (-s,c,b)    | 2.74384326       | 157.210639       | +0.38764877       | +22.210639        | +0.18508865     |
| (c,b,t)     | 2.16303821       | 123.932960       | -0.19315628       | -11.067040        | -0.09222533     |

The bloom rotation is not universal. The leptons rotate by only -2.3 deg
because m_e is tiny. The quark triples rotate in opposite directions.

The sign pattern is: (-s,c,b) rotates positive (delta increases from
seed); (e,mu,tau) and (c,b,t) rotate negative.


## Part 5: Summary table

| Triple      | M0_full    | delta_full | M0_seed    | delta_seed | M0r    | Ddelta (rad) | Ddelta (deg) |
|-------------|------------|------------|------------|------------|--------|--------------|--------------|
| (e,mu,tau)  | 313.8409   | 2.316624   | 313.5925   | 2.356194   | 1.0008 | -0.039570    | -2.26721     |
| (-s,c,b)    | 912.5557   | 2.743843   | 856.1812   | 2.356194   | 1.0658 | +0.387649    | +22.21064    |
| (c,b,t)     | 29576.4391 | 2.163038   | 29150.5971 | 2.356194   | 1.0146 | -0.193156    | -11.06704    |

| Triple      | v0r (Koide seed) | Ddelta/(2pi/3) | Ddelta/pi   | Q          |
|-------------|------------------|----------------|-------------|------------|
| (e,mu,tau)  | 1.000396         | -0.01889340    | -0.01259560 | 0.6666608  |
| (-s,c,b)    | 1.032397         | +0.18508865    | +0.12339244 | 0.6749542  |
| (c,b,t)     | 1.007278         | -0.09222533    | -0.06148355 | 0.6694894  |

(-s,c,b) bloom seed (preceding pair (s,c), sigma_3 = 0):

    v0_full / v0_seed_bloom = 2.000502   (near 2)
    M0_full / M0_seed_bloom = 4.002010   (near 4)


## Patterns observed

**Bloom rotations are not universal.** The three triples have different
Delta_delta values: -2.27, +22.21, -11.07 degrees. There is no single
"bloom rotation angle."

**Near-integer ratio between quark bloom rotations:**

    Dd(-s,c,b) / Dd(c,b,t) = -2.0069

This implies Dd(-s,c,b) + 2*Dd(c,b,t) = 0.001336 rad, which is
0.34% of Dd(-s,c,b). However, a 4% shift in m_s changes this combination
by 0.001 rad (comparable to the combination itself), so the near-2 ratio
is not a stable numerical fact given current quark mass uncertainties.

**v0-doubling is specific to the quark bloom.** The factor-of-2 in
v0_full/v0_seed applies only to (-s,c,b) with the bloom-seed definition
(preceding pair (s,c)). For leptons (factor 1.014) and (c,b,t) with its
natural bloom seed (b,t) (factor 1.074), there is no factor-of-2 doubling.

**Two distinct notions of "seed" give very different M0 ratios.** The
Koide seed (sigma_1 = 0) gives M0_full/M0_seed near 1 for all triples.
The bloom seed (preceding pair, sigma_3 = 0) gives M0_full/M0_seed = 4
for (-s,c,b). The distinction matters for interpreting v0-doubling.

**Lepton triple is a near-perfect Koide triple** (Q deviates only 0.0009%
from 2/3) with a tiny bloom rotation (-2.27 deg), consistent with m_e
being much smaller than m_mu and m_tau.

**Quark triples deviate more from Q = 2/3** at PDG central values: (+1.24%
for (-s,c,b), +0.42% for (c,b,t)), producing larger bloom rotations
relative to the seed angle.
