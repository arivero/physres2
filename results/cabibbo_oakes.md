# Cabibbo Angle from Oakes Relation and Koide Constraints

**Computation date:** 2026-03-04

## Input Data

PDG quark masses (MeV):

| Quark | Mass (MeV) |
|-------|-----------|
| u | 2.16 |
| d | 4.67 |
| s | 93.4 |
| c | 1270.0 |
| b | 4180.0 |
| t | 172760.0 |

PDG Cabibbo angle: theta_C = 13.04 deg, |V_us| = 0.2257

## Part 1: Oakes Relation (Forward)

The Weinberg-Oakes relation: tan(theta_C) = sqrt(m_d / m_s)

- m_d / m_s = 0.050000
- sqrt(m_d / m_s) = 0.223607
- theta_C = arctan(sqrt(m_d/m_s)) = 12.6044 deg
- Deviation from PDG (13.04 deg): -0.4356 deg (-3.34%)

The implied |V_us| = sin(theta_C) = 0.218218, vs PDG 0.2257: deviation -3.32%

## Part 2: Inverse Problem

Given sin(theta_C) = 0.2257, what m_d/m_s does Oakes predict?

- From tan relation: m_d/m_s = tan^2(theta_C) = 0.053640
  Predicted m_d = 5.010 MeV (PDG: 4.67 MeV, deviation: +7.28%)
- From sin relation: m_d/m_s = sin^2/(1-sin^2) = 0.053640
  Predicted m_d = 5.010 MeV (PDG: 4.67 MeV, deviation: +7.28%)
- PDG m_d/m_s = 0.050000

## Part 3: Koide Q for All Light-Quark Triples

### Unsigned sqrt (all positive)

| Triple | Q | Q - 2/3 | % deviation |
|--------|---|---------|-------------|
| (m_u, m_c, m_t) | 0.849006 | +0.182340 | +27.351% |
| (m_u, m_d, m_s) | 0.567043 | -0.099624 | -14.944% |
| (m_d, m_s, m_b) | 0.731428 | +0.064761 | +9.714% |
| (m_u, m_s, m_b) | 0.744396 | +0.077729 | +11.659% |
| (m_u, m_c, m_b) | 0.526523 | -0.140144 | -21.022% |
| (m_u, m_d, m_c) | 0.828058 | +0.161391 | +24.209% |

### Signed sqrt (negative first entry)

| Triple | Q_signed | Q - 2/3 | % deviation |
|--------|----------|---------|-------------|
| (-sqrt(m_u), sqrt(m_c), sqrt(m_t)) | 0.860139 | +0.193472 | +29.021% |
| (-sqrt(m_u), sqrt(m_d), sqrt(m_s)) | 0.934629 | +0.267963 | +40.194% |
| (-sqrt(m_d), sqrt(m_s), sqrt(m_b)) | 0.821674 | +0.155008 | +23.251% |
| (-sqrt(m_u), sqrt(m_s), sqrt(m_c)) | 0.710777 | +0.044111 | +6.617% |
| (-sqrt(m_d), sqrt(m_c), sqrt(m_b)) | 0.566466 | -0.100200 | -15.030% |
| (-sqrt(m_u), sqrt(m_d), sqrt(m_c)) | 0.967478 | +0.300811 | +45.122% |

### Reference (established triples)

| Triple | Q | Q - 2/3 | % deviation |
|--------|---|---------|-------------|
| (m_e, m_mu, m_tau) | 0.666659 | -0.000008 | -0.001% |
| (-sqrt(m_s), sqrt(m_c), sqrt(m_b)) | 0.674954 | +0.008288 | +1.243% |
| (m_c, m_b, m_t) | 0.669489 | +0.002823 | +0.423% |

**Observation**: None of the light-quark triples (u,d,s) come close to Q = 2/3.
The best is Q(m_u, m_d, m_s) = 0.5670 (unsigned), -14.9% off.
For comparison, the established triples deviate by <1%.

## Part 4: Koide Predicts m_u from (u, c, t) Triple

Solving Q(x, m_c, m_t) = 2/3 for x (with all positive signs):

The quadratic p^2 - 4S*p + 3(b+c) - 2S^2 = 0 where p = sqrt(x), S = sqrt(m_c) + sqrt(m_t):

| Solution | sqrt(m_u) | m_u (MeV) | Q check | Ratio to PDG |
|----------|-----------|-----------|---------|-------------|
| 1 | 1739.125479 | 3024557.4304 | 0.6666666667 | 1400258.0696 |
| 2 | 65.999073 | 4355.8776 | 0.6666666667 | 2016.6100 |

With signed (-sqrt(x), +sqrt(m_c), +sqrt(m_t)):

| Solution | sqrt(m_u) | m_u (MeV) | Q check | Ratio to PDG |
|----------|-----------|-----------|---------|-------------|

## Part 5: Koide Predicts m_u or m_d from (u, d, s) Triple

### 5a: Solve for m_d given m_u and m_s (all positive)

| Solution | sqrt(m_d) | m_d (MeV) | Q check | Ratio to PDG |
|----------|-----------|-----------|---------|-------------|
| 1 | 43.648580 | 1905.1985 | 0.6666666667 | 407.9654 |
| 2 | 0.887666 | 0.7880 | 0.6666666667 | 0.1687 |

### 5b: Solve for m_u given m_d and m_s (all positive)

| Solution | sqrt(m_u) | m_u (MeV) | Q check | Ratio to PDG |
|----------|-----------|-----------|---------|-------------|
| 1 | 46.992333 | 2208.2794 | 0.6666666667 | 1022.3516 |
| 2 | 0.309210 | 0.0956 | 0.6666666667 | 0.0443 |

### 5c: Solve for m_u given m_d and m_s (signed: -sqrt(m_u))

| Solution | sqrt(m_u) | m_u (MeV) | Q check | Ratio to PDG |
|----------|-----------|-----------|---------|-------------|

## Part 6: Cross-Check — Oakes + Koide Determines Everything

Starting from theta_C = 13.04 deg and m_s = 93.4 MeV:

1. Oakes: m_d = m_s * tan^2(theta_C) = 5.010 MeV (PDG: 4.67, dev: +7.28%)

2. Koide Q(m_u, 5.0, 93.4) = 2/3 predicts m_u:

| Solution | m_u (MeV) | Q check | Dev from PDG |
|----------|-----------|---------|-------------|
| 1 (all +) | 2242.9483 | 0.6666666667 | +103740.20% |
| 2 (all +) | 0.0630 | 0.6666666667 | -97.09% |

### Parameter counting

| Constraint | Fixes | Source |
|-----------|-------|--------|
| Koide seed: m_c/m_s = (2+sqrt(3))^2 | m_c from m_s | sBootstrap |
| v0-doubling: sqrt(m_b) = 3*sqrt(m_s) + sqrt(m_c) | m_b from m_s, m_c | sBootstrap |
| (c,b,t) Koide Q = 2/3 | overdetermination test for m_t | sBootstrap |
| Top Yukawa at tan(beta) = 1 | m_t (external) | EWSB |
| Oakes: tan(theta_C) = sqrt(m_d/m_s) | m_d/m_s from CKM | Weinberg-Oakes |
| Koide Q(u,d,s) = 2/3 | m_u from m_d, m_s | hypothetical |

Total: 6 masses, 1 overall scale (m_s), 4+1 ratio constraints + 1 external = **all ratios fixed**

**But**: Q(m_u, m_d, m_s) = 2/3 is NOT satisfied by PDG masses.
With PDG masses, Q(m_u, m_d, m_s) = 0.5670 (unsigned), -14.9% from 2/3.
The Koide condition does not hold for the (u,d,s) triple at current precision.
So this last constraint is hypothetical, not empirical.

## Part 7: Alternative Oakes Relations

| Relation | Formula | Value | theta_C (deg) | |V_us| | Dev from PDG |
|----------|---------|-------|---------------|--------|-------------|
| Weinberg-Oakes | tan(tC) = sqrt(m_d/m_s) | 0.223607 | 12.6044 | 0.218218 | -3.315% |
| Gatto-Sartori-Tonin | sin(tC) = sqrt(m_d/m_s) | 0.223607 | 12.9210 | 0.223607 | -0.927% |
| Fritzsch | tan(tC) = sqrt(md/ms) - sqrt(mu/mc) | 0.182366 | 10.3352 | 0.179407 | -20.511% |
| PDG | measured | — | 13.0400 | 0.225700 | — |

**Notes:**
- Weinberg-Oakes and GST differ only in whether sqrt(m_d/m_s) is set equal to tan or sin;
  numerically they give theta_C = 12.604 vs 12.921 deg (difference: 0.317 deg).
- The Fritzsch correction subtracts sqrt(m_u/m_c) = 0.041241, reducing theta_C by about 2.269 deg.
- The Fritzsch correction *worsens* the match because sqrt(m_d/m_s) is already below |V_us|;
  subtracting sqrt(m_u/m_c) pushes it further away. With a complex phase, the Fritzsch
  texture allows |V_us| in [0.1824, 0.2648], which brackets the PDG value.
- The GST relation (sin = sqrt(m_d/m_s)) gives the best match at -0.9%.

## Summary Table

| # | Result | Value | Reference | Deviation |
|---|--------|-------|-----------|-----------|
| 1 | theta_C from Oakes (tan) | 12.604 deg | 13.04 deg | -3.34% |
| 2 | theta_C from GST (sin) | 12.921 deg | 13.04 deg | -0.91% |
| 3 | theta_C from Fritzsch | 10.335 deg | 13.04 deg | -20.74% |
| 4 | m_d/m_s from Oakes inverse | 0.0536 | 0.0500 (PDG) | +7.28% |
| 5 | Q(m_u, m_c, m_t) unsigned | 0.8490 | 0.6667 | +27.35% |
| 6 | Q(-su, sc, st) signed | 0.8601 | 0.6667 | +29.02% |
| 7 | Q(m_u, m_d, m_s) unsigned | 0.5670 | 0.6667 | -14.94% |
| 8 | Q(-su, sd, ss) signed | 0.9346 | 0.6667 | +40.19% |
| 9.1 | m_u from Q(u,c,t)=2/3 sol 1 | 3024557.43 MeV | 2.16 MeV | +140025707.0% |
| 9.2 | m_u from Q(u,c,t)=2/3 sol 2 | 4355.88 MeV | 2.16 MeV | +201561.0% |
| 10.1 | m_u from Q(u,d,s)=2/3 sol 1 | 2208.28 MeV | 2.16 MeV | +102135.2% |
| 10.2 | m_u from Q(u,d,s)=2/3 sol 2 | 0.10 MeV | 2.16 MeV | -95.6% |

## Conclusions

1. **Oakes relation works well.** tan(theta_C) = sqrt(m_d/m_s) gives theta_C = 12.60 deg,
   only -3.34% from the PDG value. This is a robust, long-established
   result (Weinberg 1977, Oakes 1969) and remains numerically accurate with current PDG masses.

2. **GST (sin = sqrt(m_d/m_s)) is the best variant**, at -0.9% from PDG.
   The Fritzsch correction (subtracting sqrt(m_u/m_c) = 0.0412) makes things worse
   because sqrt(m_d/m_s) is already below |V_us|. With a complex phase, Fritzsch brackets the PDG value.

3. **Koide does NOT hold for (u,d,s).** Q(m_u, m_d, m_s) = 0.5670 (unsigned),
   which is -14.9% from 2/3. This is qualitatively different from the ~0.5% deviations
   seen in the established Koide triples. The (u,d,s) triple does not satisfy the Koide condition.

4. **Koide on (u,c,t) also fails badly.** Q(m_u, m_c, m_t) = 0.8490 (unsigned),
   +27.4% off. The Koide condition constrains the established triples
   (-s,c,b), (c,b,t), (e,mu,tau), but not those containing m_u or m_d.

5. **Parameter counting remains at 2 free.** The Koide seed + v0-doubling + top Yukawa
   fix (m_c, m_b, m_t) in terms of m_s. The remaining free parameters are m_u and m_d.
   The Oakes relation connects m_d/m_s to theta_C, so if theta_C is taken as input,
   only m_u remains free. But no Koide condition currently fixes m_u.

6. **The free parameters (m_u, m_d) connect to CKM mixing**, not to further mass
   relations. This is consistent with the observation that ALL Koide triples near Q = 2/3
   are mixed (up+down type), suggesting the Koide structure is fundamentally intertwined
   with CKM mixing.