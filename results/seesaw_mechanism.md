# Seesaw Mechanism and Dual Koide Analysis

## Setup

In SU(3) SQCD with N_f = N_c = 3, the Seiberg duality maps quarks to
mesons with VEVs M_j = C'/m_j, where C' = Lambda^2 (m_d m_s m_b)^{1/3}.
This produces an inverted mass spectrum.

## 1. Seesaw Spectrum

With Lambda = 300 MeV:
- C' = 10994846.9670 MeV
- M_d = C'/m_d = 2354356.95 MeV
- M_s = C'/m_s = 117717.85 MeV
- M_b = C'/m_b = 2630.35 MeV

The hierarchy is inverted: M_d > M_s > M_b.

## 2. C' Cancellation

Q(M_d, M_s, M_b) = Q(1/m_d, 1/m_s, 1/m_b) exactly (C' cancels).
- Q(1/m_d, 1/m_s, 1/m_b) = 0.66520987  (2/3 deviation: -0.2185%)
- Q(m_d, m_s, m_b) = 0.73142765  (2/3 deviation: +9.7141%)

The seesaw transforms a non-Koide spectrum (Q = 0.731) into a near-Koide one (Q = 0.665).

## 3. Algebraic Condition

With symmetric polynomials of a_j = sqrt(m_j):
- Q(m) = 1 - 2 p2/p1^2
- Q(1/m) = 1 - 2 p1 p3/p2^2

These are independent conditions. Q_dual ~ 2/3 requires p1 p3/p2^2 ~ 1/6.
Actual: p1 p3/p2^2 = 0.16739506 vs 1/6 = 0.16666667 (+0.437%).
The leading-order approximation a_d/a_s = 0.2236 is 34% from 1/6;
subleading terms bring the full ratio to within 0.44%.

## 4. Parametric Analysis

### Charged leptons (e, mu, tau)
- M0 = 313.841127 MeV, delta = 132.7328 deg
- Q = 0.6666605115, Q_dual = 0.85144869

### Down-type quarks (d, s, b) [direct]
- M0 = 649.881243 MeV, delta = 126.3126 deg
- Q = 0.73142765 (NOT on Koide manifold)
- Q_dual = 0.66520987

### Inverse down-type quarks (1/m_d, 1/m_s, 1/m_b)
- M0_dual = 3.7595258732e-02 (1/MeV), delta_dual = 130.7021 deg
- This IS on the Koide manifold (Q = 0.66520987 ~ 2/3)
- cos(3 delta_dual) = 0.84706239
- Required for exact Q = 2/3: cos(3 delta) = 0.88388348

## 5. Delta Relations

| Triple | delta (deg) |
|--------|------------|
| (e, mu, tau) | 132.7328 |
| (-s, c, b) | 157.2106 |
| (c, b, t) | 123.9330 |
| (1/m_d, 1/m_s, 1/m_b) | 130.7021 |

No simple algebraic relation between delta_dual and the other
deltas has been identified. The dual Koide is an independent
numerical observation.

## 6. Mass Prediction from Exact Dual Koide

Setting Q(1/m_d, 1/m_s, 1/m_b) = 2/3 exactly with m_d = 4.67, m_s = 93.4:

- Predicted m_b = 4562.0 MeV
- PDG m_b = 4180 +/- 30 MeV
- Deviation: +9.139% (+12.73 sigma)

### Comparison of m_b predictions

| Method | m_b (MeV) | PDG deviation |
|--------|----------|--------------|
| PDG central | 4180 | -- |
| v0-doubling | 4177.1 | -0.070% |
| Dual Koide | 4562.0 | +9.139% |
| Overlap | 4159 | -0.502% |

## 7. Compatibility Test

v0-doubling condition: sqrt(m_b) = 3 sqrt(m_s) + sqrt(m_c)
Dual Koide condition: Q(1/m_d, 1/m_s, 1/m_b) = 2/3

v0-doubling gives m_b = 4177.1 MeV (uses m_s, m_c)
Dual Koide gives m_b = 4562.0 MeV (uses m_d, m_s)
Difference: 385.0 MeV (12.83 sigma)


The exact dual Koide condition predicts m_b = 4562 MeV (9% off PDG),
while v0-doubling predicts 4177 MeV (0.07% off). The 0.22% closeness
of Q_dual to 2/3 does NOT translate to a precise m_b prediction because
the Koide equation amplifies small Q deviations into large mass deviations.

However, both conditions CAN be made simultaneously compatible by
adjusting m_d to 4.604 MeV (0.14 sigma from PDG).

## Key Takeaways

1. The Seiberg seesaw maps the down-type quark spectrum to a near-Koide
   spectrum. This is not a consequence of Q(m) ~ 2/3 (which fails for d,s,b)
   but an independent algebraic property of the mass ratios.

2. The condition p1*p3/p2^2 = 1/6 holds to 0.44%. The leading-order
   approximation (sqrt(m_d)/sqrt(m_s) ~ 1/6) is 34% off;
   the sub-percent agreement requires all three masses.

3. IMPORTANT: The exact dual Koide condition Q = 2/3 predicts m_b = 4562 MeV
   (9% off PDG). The 0.22% closeness of Q to 2/3 does NOT translate to a
   precise mass prediction, because the Koide quadratic amplifies small Q
   deviations. The v0-doubling prediction (m_b = 4177, 0.07% off) is far superior.

4. Both conditions can be made simultaneously compatible with a small shift
   in m_d (0.14 sigma from PDG), suggesting they are not in fundamental tension.

5. The ISS seesaw mechanism provides a natural physical framework
   for the dual Koide: if SQCD confines with N_f = 3 down-type flavors,
   the meson VEVs automatically produce the inverted spectrum.
