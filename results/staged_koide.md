# Staged Koide Analysis

Parametrization:  m_i = M0 * (1 + sqrt(2)*cos(2*pi*(i-1)/3 + delta))^2
This automatically gives Q = 2/3 for any (M0, delta).

Definitions:
  sigma_i  = s_i * sqrt(m_i)  (signed square root)
  v0       = (1/3) * sum_i sigma_i   (vacuum charge; equals sqrt(M0) in the parametrization)
  M0       = v0^2  (since sum_k f_k = 3 always)

Seed: mass slot k vanishes when f_k(delta) = 0, i.e.
  delta_seed = 3*pi/4 - 2*pi*k/3  (mod 2*pi)  [principal branch]
  delta_seed = 5*pi/4 - 2*pi*k/3  (mod 2*pi)  [second branch]

v0-doubling definition (from v0_doubling.md):
  v0(seed) = (1/3) * sum_{nonzero slots} sigma_i
  v0(full) = (1/3) * sum_i sigma_i  (all three)
  Ratio v0(full)/v0(seed) -- the v0-doubling ratio.

## Part 1: Fitted Parameters for Each Triple

### Triple (e, mu, tau)

Physical masses (MeV):  e = 0.51099895,  mu = 105.6583755,  tau = 1776.86
Signed sqrt masses (sigma):  (0.714842, 10.279026, 42.152817)  sqrt(MeV)
Koide Q = 0.66666051   (exact 2/3 = 0.66666667)
v0_full = (1/3)*sum(sigma) = 17.715562  sqrt(MeV)
M0 from v0^2 = 313.8411  MeV  (valid for exact Q=2/3)

Fitted (M0, delta) by least-squares on sigma values:
  M0_full    = 313.838229 MeV
  v0_full (from sqrt(M0_full)) = 17.715480  sqrt(MeV)
  delta_full = 4.41101984 rad  =  252.732820 deg  =  1.40407122 pi
  Index mapping: model slot (0,1,2) -> physical index (1, 2, 0)
  RMS fractional reconstruction error = 0.0005%

Reconstruction check (model slot order):
  slot 0: mu  model = 105.6560 MeV  physical = 105.6584 MeV  err = -0.0023%
  slot 1: tau  model = 1776.8626 MeV  physical = 1776.8600 MeV  err = +0.0001%
  slot 2: e  model = 0.5108 MeV  physical = 0.5110 MeV  err = -0.0448%

Seed parameters:
  delta_seed = 4.45058959 rad  =  255.000000 deg  =  1.41666667 pi
  Zero slot: model slot 2  ->  physical mass e
  Bloom rotation: Delta_delta = delta_full - delta_seed = -0.039570 rad  =  -2.267180 deg
  Delta_delta / (2*pi/3) = -0.018893

v0 values (physical definition: (1/3)*sum sigma):
  v0_full = (1/3)*(0.7148+10.2790+42.1528) = 17.715562  sqrt(MeV)
  v0_seed = (1/3)*(sigma_mu + sigma_tau) = (1/3)*(10.2790 + 42.1528) = 17.477281  sqrt(MeV)
  v0 ratio: v0_full / v0_seed = 1.013634   (ideal: 2.0000)

Seed triple masses (M0_seed = v0_seed^2 = 305.4554 MeV, delta_seed):
  mu_seed = 122.7698 MeV   physical = 105.6584 MeV   ratio = 1.161950
  tau_seed = 1709.9623 MeV   physical = 1776.8600 MeV   ratio = 0.962351
  e_seed = 3.765e-28 MeV  (zero slot)


### Triple (-s, c, b)

Physical masses (MeV):  s = 93.4,  c = 1270.0,  b = 4180.0
Signed sqrt masses (sigma):  (-9.664368, 35.637059, 64.652920)  sqrt(MeV)
Koide Q = 0.67495422   (exact 2/3 = 0.66666667)
v0_full = (1/3)*sum(sigma) = 30.208537  sqrt(MeV)
M0 from v0^2 = 912.5557  MeV  (valid for exact Q=2/3)

Fitted (M0, delta) by least-squares on sigma values:
  M0_full    = 923.865175 MeV
  v0_full (from sqrt(M0_full)) = 30.395151  sqrt(MeV)
  delta_full = 2.74384326 rad  =  157.210639 deg  =  0.87339244 pi
  Index mapping: model slot (0,1,2) -> physical index (0, 1, 2)
  RMS fractional reconstruction error = 0.6139%

Reconstruction check (model slot order):
  slot 0: s  model = 85.2750 MeV  physical = 93.4000 MeV  err = -8.6992%
  slot 1: c  model = 1280.9633 MeV  physical = 1270.0000 MeV  err = +0.8633%
  slot 2: b  model = 4176.9528 MeV  physical = 4180.0000 MeV  err = -0.0729%

Seed parameters:
  delta_seed = 2.35619449 rad  =  135.000000 deg  =  0.75000000 pi
  Zero slot: model slot 0  ->  physical mass s
  Bloom rotation: Delta_delta = delta_full - delta_seed = 0.387649 rad  =  22.210639 deg
  Delta_delta / (2*pi/3) = 0.185089

v0 values (physical definition: (1/3)*sum sigma):
  v0_full = (1/3)*(-9.6644+35.6371+64.6529) = 30.208537  sqrt(MeV)
  v0_seed = (1/3)*(sigma_c + sigma_b) = (1/3)*(35.6371 + 64.6529) = 33.429993  sqrt(MeV)
  v0 ratio: v0_full / v0_seed = 0.903636   (ideal: 2.0000)

Seed triple masses (M0_seed = v0_seed^2 = 1117.5644 MeV, delta_seed):
  s_seed = 0.000e+00 MeV  (zero slot)
  c_seed = 449.1757 MeV   physical = 1270.0000 MeV   ratio = 0.353682
  b_seed = 6256.2108 MeV   physical = 4180.0000 MeV   ratio = 1.496701


### Triple (c, b, t)

Physical masses (MeV):  c = 1270.0,  b = 4180.0,  t = 172760.0
Signed sqrt masses (sigma):  (35.637059, 64.652920, 415.644079)  sqrt(MeV)
Koide Q = 0.66948936   (exact 2/3 = 0.66666667)
v0_full = (1/3)*sum(sigma) = 171.978019  sqrt(MeV)
M0 from v0^2 = 29576.4391  MeV  (valid for exact Q=2/3)

Fitted (M0, delta) by least-squares on sigma values:
  M0_full    = 29701.534671 MeV
  v0_full (from sqrt(M0_full)) = 172.341332  sqrt(MeV)
  delta_full = 6.21454220 rad  =  356.067040 deg  =  1.97815022 pi
  Index mapping: model slot (0,1,2) -> physical index (2, 1, 0)
  RMS fractional reconstruction error = 0.2108%

Reconstruction check (model slot order):
  slot 0: t  model = 172635.9284 MeV  physical = 172760.0000 MeV  err = -0.0718%
  slot 1: b  model = 4256.5196 MeV  physical = 4180.0000 MeV  err = +1.8306%
  slot 2: c  model = 1316.7600 MeV  physical = 1270.0000 MeV  err = +3.6819%

Seed parameters:
  delta_seed = 6.02138592 rad  =  345.000000 deg  =  1.91666667 pi
  Zero slot: model slot 2  ->  physical mass c
  Bloom rotation: Delta_delta = delta_full - delta_seed = 0.193156 rad  =  11.067040 deg
  Delta_delta / (2*pi/3) = 0.092225

v0 values (physical definition: (1/3)*sum sigma):
  v0_full = (1/3)*(35.6371+64.6529+415.6441) = 171.978019  sqrt(MeV)
  v0_seed = (1/3)*(sigma_t + sigma_b) = (1/3)*(415.6441 + 64.6529) = 160.098999  sqrt(MeV)
  v0 ratio: v0_full / v0_seed = 1.074198   (ideal: 2.0000)

Seed triple masses (M0_seed = v0_seed^2 = 25631.6896 MeV, delta_seed):
  t_seed = 143488.1517 MeV   physical = 172760.0000 MeV   ratio = 0.830564
  b_seed = 10301.9858 MeV   physical = 4180.0000 MeV   ratio = 2.464590
  c_seed = 5.055e-27 MeV  (zero slot)


## Part 2: The Staged Picture

Hypothesis: quark mass spectrum arises in stages:
  Stage 1 (EW scale): (c,b,t) triple; top Yukawa dominates. Seed has c_seed = 0.
  Stage 2 (chiral):  (-s,c,b) triple; s acquires mass. Seed has s_seed = 0.
  Stage 3 (bloom):   s acquires its physical mass, rotating delta from delta_seed.

### Stage 1: (c, b, t) — EW-determined
Zero mass at seed: c
delta_seed = 6.021386 rad = 1.916667 pi = 345.0000 deg
delta_full = 6.214542 rad = 1.978150 pi = 356.0670 deg
Bloom Delta_delta = 0.193156 rad = 11.0670 deg

Seed triple masses:
  t_seed = 143488.1517 MeV   physical = 172760.0000 MeV   ratio = 0.830564
  b_seed = 10301.9858 MeV   physical = 4180.0000 MeV   ratio = 2.464590
  c_seed = 0  (zero slot; physical mass = 1270.0000 MeV)

### Stage 2: (-s, c, b) — chiral
Zero mass at seed: s
delta_seed = 2.356194 rad = 0.750000 pi = 135.0000 deg
delta_full = 2.743843 rad = 0.873392 pi = 157.2106 deg
Bloom Delta_delta = 0.387649 rad = 22.2106 deg

Seed triple masses:
  s_seed = 0  (zero slot; physical mass = 93.4000 MeV)
  c_seed = 449.1757 MeV   physical = 1270.0000 MeV   ratio = 0.353682
  b_seed = 6256.2108 MeV   physical = 4180.0000 MeV   ratio = 1.496701

### Stage 3: (e, mu, tau) — leptonic
Zero mass at seed: e
delta_seed = 4.450590 rad = 1.416667 pi = 255.0000 deg
delta_full = 4.411020 rad = 1.404071 pi = 252.7328 deg
Bloom Delta_delta = -0.039570 rad = -2.2672 deg

Seed triple masses:
  mu_seed = 122.7698 MeV   physical = 105.6584 MeV   ratio = 1.161950
  tau_seed = 1709.9623 MeV   physical = 1776.8600 MeV   ratio = 0.962351
  e_seed = 0  (zero slot; physical mass = 0.5110 MeV)


## Part 3: Overlap Prediction from Dual Q=2/3 Constraints

Given m_s = 93.4 MeV and m_t = 172760 MeV as inputs,
find (m_c, m_b) satisfying Q(-s,c,b) = 2/3 AND Q(c,b,t) = 2/3 simultaneously.

Solution 1:
  m_c (predicted) = 1369.1272 MeV   physical = 1270.0000 MeV   deviation = +7.805%
  m_b (predicted) = 4159.3740 MeV   physical = 4180.0000 MeV   deviation = -0.493%
  Verification: Q(-s,c,b) = 0.6666666667   Q(c,b,t) = 0.6666666667
                (both should equal 2/3 = 0.6666666667)


## Part 4: v0-Doubling

v0 = (1/3)*sum(sigma_i) where sigma_i are the signed square roots.
v0(seed) = (1/3)*sum of the two nonzero sigma values at the seed.
v0-doubling hypothesis: v0(full)/v0(seed) = 2.

### Triple (e, mu, tau)
  sigma values: e=0.714842, mu=10.279026, tau=42.152817  sqrt(MeV)
  v0_full = (1/3)*(0.7148+10.2790+42.1528) = 17.715562  sqrt(MeV)
  v0_seed = (1/3)*(sigma_mu+sigma_tau) = (1/3)*(10.2790+42.1528) = 17.477281  sqrt(MeV)
  v0_full/v0_seed = 1.013634   ideal: 2.000   deviation: 49.318%

### Triple (-s, c, b)
  sigma values: s=-9.664368, c=35.637059, b=64.652920  sqrt(MeV)
  v0_full = (1/3)*(-9.6644+35.6371+64.6529) = 30.208537  sqrt(MeV)
  v0_seed = (1/3)*(sigma_c+sigma_b) = (1/3)*(35.6371+64.6529) = 33.429993  sqrt(MeV)
  v0_full/v0_seed = 0.903636   ideal: 2.000   deviation: 54.818%

### Triple (c, b, t)
  sigma values: c=35.637059, b=64.652920, t=415.644079  sqrt(MeV)
  v0_full = (1/3)*(35.6371+64.6529+415.6441) = 171.978019  sqrt(MeV)
  v0_seed = (1/3)*(sigma_t+sigma_b) = (1/3)*(415.6441+64.6529) = 160.098999  sqrt(MeV)
  v0_full/v0_seed = 1.074198   ideal: 2.000   deviation: 46.290%

### Cross-check: explicit v0-doubling formula for quarks

From v0_doubling.md: the quark bloom has
  v0(seed) = (sqrt(m_s) + sqrt(m_c))/3
  v0(full) = (-sqrt(m_s) + sqrt(m_c) + sqrt(m_b))/3

  sqrt(m_s) = 9.664368, sqrt(m_c) = 35.637059, sqrt(m_b) = 64.652920  sqrt(MeV)
  v0(seed) = (9.664368 + 35.637059)/3 = 15.100476
  v0(full) = (-9.664368 + 35.637059 + 64.652920)/3 = 30.208537
  Ratio v0(full)/v0(seed) = 2.000502   (reported value: 2.0005)

v0-doubling prediction for m_b: set ratio = 2 exactly.
  2*v0(seed) = v0(full):
  2*(sqrt(m_s)+sqrt(m_c))/3 = (-sqrt(m_s)+sqrt(m_c)+sqrt(m_b))/3
  2*sqrt(m_s)+2*sqrt(m_c) = -sqrt(m_s)+sqrt(m_c)+sqrt(m_b)
  sqrt(m_b) = 3*sqrt(m_s) + sqrt(m_c)
  sqrt(m_b) = 3*9.664368 + 35.637059 = 64.630162
  m_b (predicted) = 4177.0578 MeV   physical = 4180.0000 MeV   deviation = -0.0704%

### Cross-consistency: (c,b,t) seed b-mass vs physical b

The (c,b,t) seed sets c=0. Seed masses are computed from M0_seed = v0_seed^2.
v0_seed(c,b,t) = 160.098999  sqrt(MeV)
M0_seed(c,b,t) = 25631.6896 MeV

  t_seed = 143488.1517 MeV   physical = 172760.0000 MeV   ratio = 0.830564
  b_seed = 10301.9858 MeV   physical = 4180.0000 MeV   ratio = 2.464590
  c_seed = 5.055e-27 MeV  (zero)


## Part 5: Summary Tables

### Table 1: Koide Q Values

| Triple | Signs | Q (direct) | Q - 2/3 |
|--------|-------|-----------|---------|
| (e, mu, tau) | (+,+,+) | 0.66666051 | -6.16e-06 |
| (-s, c, b) | (-,+,+) | 0.67495422 | +8.29e-03 |
| (c, b, t) | (+,+,+) | 0.66948936 | +2.82e-03 |

### Table 2: Full Triple Parameters

| Triple | M0_full (MeV) | delta_full / pi | v0_full (sqrt(MeV)) |
|--------|--------------|----------------|---------------------|
| (e, mu, tau) | 313.8382 | 1.404071 | 17.715562 |
| (-s, c, b) | 923.8652 | 0.873392 | 30.208537 |
| (c, b, t) | 29701.5347 | 1.978150 | 171.978019 |

### Table 3: Seed Parameters

| Triple | Zero mass | delta_seed / pi | v0_seed (sqrt(MeV)) | M0_seed (MeV) |
|--------|-----------|----------------|--------------------|--------------| 
| (e, mu, tau) | e | 1.416667 | 17.477281 | 305.4554 |
| (-s, c, b) | s | 0.750000 | 33.429993 | 1117.5644 |
| (c, b, t) | c | 1.916667 | 160.098999 | 25631.6896 |

### Table 4: Bloom Rotations

| Triple | Delta_delta (rad) | Delta_delta (deg) | Delta_delta / (2*pi/3) |
|--------|------------------|------------------|------------------------|
| (e, mu, tau) | -0.039570 | -2.2672 | -0.018893 |
| (-s, c, b) | 0.387649 | 22.2106 | 0.185089 |
| (c, b, t) | 0.193156 | 11.0670 | 0.092225 |

### Table 5: v0-Doubling Ratios

v0 = (1/3)*sum(sigma_i).  v0_seed = (1/3)*sum of nonzero sigma values at seed.

| Triple | Zero mass | v0_full | v0_seed | v0_full/v0_seed | Deviation from 2 |
|--------|-----------|--------|--------|----------------|-----------------|
| (e, mu, tau) | e | 17.71556 | 17.47728 | 1.013634 | 49.318% |
| (-s, c, b) | s | 30.20854 | 33.42999 | 0.903636 | 54.818% |
| (c, b, t) | c | 171.97802 | 160.09900 | 1.074198 | 46.290% |

### Table 6: Overlap Prediction

Inputs: m_s = 93.4 MeV,  m_t = 172760.0 MeV

| Mass | Predicted (MeV) | Physical (MeV) | Deviation |
|------|----------------|---------------|-----------|
| m_c | 1369.127 | 1270.000 | +7.805% |
| m_b | 4159.374 | 4180.000 | -0.493% |

### Table 7: Seed vs Physical Mass Comparison

Seed masses are computed from M0_seed = v0_seed^2 and delta_seed.
The nonzero seed masses are NOT equal to the physical masses because
f_k(delta_seed) differs from f_k(delta_full).

| Triple | Mass | Seed (MeV) | Physical (MeV) | seed/physical |
|--------|------|-----------|---------------|---------------|
| (e, mu, tau) | mu | 122.7698 | 105.6584 | 1.161950 |
| (e, mu, tau) | tau | 1709.9623 | 1776.8600 | 0.962351 |
| (e, mu, tau) | e | 0 | 0.5110 | — |
| (-s, c, b) | s | 0 | 93.4000 | — |
| (-s, c, b) | c | 449.1757 | 1270.0000 | 0.353682 |
| (-s, c, b) | b | 6256.2108 | 4180.0000 | 1.496701 |
| (c, b, t) | t | 143488.1517 | 172760.0000 | 0.830564 |
| (c, b, t) | b | 10301.9858 | 4180.0000 | 2.464590 |
| (c, b, t) | c | 0 | 1270.0000 | — |

---
All masses in MeV.  Exact Koide condition: Q = 2/3 = 0.666...


## Appendix: Seed Angle Reference

Principal seed angles delta_seed = 3*pi/4 - 2*pi*k/3 (mod 2*pi):
  k=0: delta_seed = 2.35619449 rad = 0.75000000 pi = 135.0000 deg,  f_k = 0.00e+00
  k=1: delta_seed = 0.26179939 rad = 0.08333333 pi = 15.0000 deg,  f_k = 0.00e+00
  k=2: delta_seed = 4.45058959 rad = 1.41666667 pi = 255.0000 deg,  f_k = 1.11e-15

Second branch delta_seed = 5*pi/4 - 2*pi*k/3 (mod 2*pi):
  k=0: delta_seed = 3.92699082 rad = 1.25000000 pi = 225.0000 deg,  f_k = -2.22e-16
  k=1: delta_seed = 1.83259571 rad = 0.58333333 pi = 105.0000 deg,  f_k = -2.22e-16
  k=2: delta_seed = 6.02138592 rad = 1.91666667 pi = 345.0000 deg,  f_k = 4.44e-16

Note: delta_seed/pi values are exact fractions:
  Principal branch: 3/4, 3/4-2/3=1/12, 3/4-4/3=-7/12 = 5/12 (mod 1) ... checking:
  k=0: 0.75 pi  (3/4 - 0*2/3 mod 1 = 0.7500000000)
  k=1: 0.08333333333333337 pi  (3/4 - 1*2/3 mod 1 = 0.0833333333)
  k=2: 0.41666666666666674 pi  (3/4 - 2*2/3 mod 1 = 0.4166666667)

  Second branch: 5/4, 5/4-2/3=7/12, 5/4-4/3=-1/12 = 11/12 ... checking:
  k=0: 0.25 pi  (5/4 - 0*2/3 mod 1 = 0.2500000000)
  k=1: 0.5833333333333334 pi  (5/4 - 1*2/3 mod 1 = 0.5833333333)
  k=2: 0.9166666666666667 pi  (5/4 - 2*2/3 mod 1 = 0.9166666667)
