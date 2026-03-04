# Exhaustive Koide Scan: SM Masses + EW Scales

**Formula**: Q = (m_1 + m_2 + m_3) / (sqrt(m_1) + sqrt(m_2) + sqrt(m_3))^2

Q = 2/3 is the Koide value. Q = 1/3 for equal masses. Q -> 1 for extreme hierarchy.

**Signed variant**: Q_s = (m_1 + m_2 + m_3) / (-sqrt(m_min) + sqrt(m_mid) + sqrt(m_max))^2

**Input masses** (GeV): m_e = 0.000511, m_mu = 0.10566, m_tau = 1.77686, m_u = 0.00216, m_d = 0.00467, m_s = 0.0934, m_c = 1.270, m_b = 4.180, m_t = 172.76, m_W = 80.369, m_Z = 91.188, m_H = 125.25, f_pi = 0.092, v_EW = 246.22

---

## 1. Top 30 Triples (Unsigned Q)

All C(14,3) = 364 unordered triples computed.

| Rank | Triple | Q | |Q - 2/3| | % dev |
|------|--------|---|----------|-------|
| 1 | **(e, mu, tau)** | **0.66666049** | **6.18e-06** | **0.0009%** |
| 2 | **(c, b, t)** | **0.66948936** | **2.82e-03** | **0.42%** |
| 3 | (b, W, f_pi) | 0.66137621 | 5.29e-03 | 0.79% |
| 4 | (mu, b, Z) | 0.67207606 | 5.41e-03 | 0.81% |
| 5 | (s, b, W) | 0.66111840 | 5.55e-03 | 0.83% |
| 6 | (tau, u, f_pi) | 0.66072962 | 5.94e-03 | 0.89% |
| 7 | (tau, u, s) | 0.65942090 | 7.25e-03 | 1.09% |
| 8 | (s, b, Z) | 0.67418691 | 7.52e-03 | 1.13% |
| 9 | (mu, b, W) | 0.65894828 | 7.72e-03 | 1.16% |
| 10 | (b, Z, f_pi) | 0.67443762 | 7.77e-03 | 1.17% |
| 11 | (e, tau, s) | 0.67791043 | 1.12e-02 | 1.69% |
| 12 | (tau, b, t) | 0.65475129 | 1.19e-02 | 1.79% |
| 13 | (e, tau, f_pi) | 0.67928207 | 1.26e-02 | 1.89% |
| 14 | (tau, c, Z) | 0.65340895 | 1.33e-02 | 1.99% |
| 15 | (e, u, f_pi) | 0.68266353 | 1.60e-02 | 2.40% |
| 16 | (e, u, s) | 0.68428338 | 1.76e-02 | 2.64% |
| 17 | (mu, tau, u) | 0.64868661 | 1.80e-02 | 2.70% |
| 18 | (e, c, f_pi) | 0.64549168 | 2.12e-02 | 3.18% |
| 19 | (tau, c, H) | 0.68842825 | 2.18e-02 | 3.26% |
| 20 | (tau, d, f_pi) | 0.64475470 | 2.19e-02 | 3.29% |

The gap between rank 1 and rank 2 is enormous: the (e, mu, tau) result is 450x closer to 2/3 than the next-best triple.

Rank 2 (c, b, t) is the only other triple within 1%.

Ranks 3-10 are a cluster of accidental near-hits at the 0.8%-1.2% level. These involve m_b paired with EW boson masses (W, Z) or QCD/light scales (f_pi), where the mass ratios happen to fall near the Koide ratio. They are NOT physically motivated combinations.

---

## 2. Best Signed-Variant Triples

The signed variant (flipping smallest sqrt) produces a different ranking:

| Rank | Triple | Q_unsigned | Q_signed | Best |Q - 2/3| | Which |
|------|--------|-----------|----------|---------------|-------|
| 1 | (e, mu, tau) | 0.66666049 | 0.70402878 | 6.18e-06 | unsigned |
| 2 | (c, b, t) | 0.66948936 | 0.90131389 | 2.82e-03 | unsigned |
| 7 | **(e, mu, c)** | 0.63288350 | **0.67355197** | **6.89e-03** | **signed** |
| 8 | **(c, b, f_pi)** | 0.45900376 | **0.67370237** | **7.04e-03** | **signed** |
| 13 | **(tau, c, f_pi)** | 0.41108588 | **0.67487905** | **8.21e-03** | **signed** |
| 14 | **(-s, c, b)** | 0.45851275 | **0.67495422** | **8.29e-03** | **signed** |
| 15 | **(tau, s, c)** | 0.41058571 | **0.67662196** | **9.96e-03** | **signed** |

The known (-s, c, b) Koide triple appears at signed rank 14, with Q_signed = 0.6750, |Q - 2/3| = 8.3e-3 (1.2%).

---

## 3. Known Koide Triples — Verification

| Triple | Q_unsigned | |Q - 2/3| | Q_signed | |Q_s - 2/3| |
|--------|-----------|----------|----------|-----------|
| (e, mu, tau) | 0.6666604894 | 6.18e-06 | 0.7040 | 3.74e-02 |
| (-s, c, b) | 0.4585 | 2.08e-01 | **0.6749542234** | **8.29e-03** |
| (c, b, t) | **0.6694893550** | **2.82e-03** | 0.9013 | 2.35e-01 |
| (d, s, b) | 0.7314 | 6.48e-02 | 0.8217 | 1.55e-01 |
| (u, c, t) | 0.8490 | 1.82e-01 | 0.8601 | 1.93e-01 |

Only (e, mu, tau) and (c, b, t) are close for the unsigned formula. Only (-s, c, b) is close in signed form. The down-type (d, s, b) and up-type (u, c, t) unsigned Koide ratios are far from 2/3 — they are NOT Koide triples in any meaningful sense.

---

## 4. Physically Motivated Triples

| Triple | Q | |Q - 2/3| | % dev | Q_signed | |Q_s - 2/3| |
|--------|---|----------|-------|----------|-----------|
| (W, Z, H) | 0.3364 | 3.30e-01 | 49.5% | 2.140 | 1.47 |
| (W, Z, t) | 0.3436 | 3.23e-01 | 48.5% | 1.827 | 1.16 |
| (t, H, v_EW) | 0.3397 | 3.27e-01 | 49.0% | 1.748 | 1.08 |
| (W, H, t) | 0.3412 | 3.25e-01 | 48.8% | 1.602 | 0.94 |
| (c, b, W) | 0.5827 | 8.40e-02 | 12.6% | 0.879 | 0.21 |
| (b, t, W) | 0.4411 | 2.26e-01 | 33.8% | **0.639** | **2.75e-02** |
| **(b, t, H)** | **0.4342** | **2.32e-01** | 34.9% | **0.608** | **5.85e-02** |
| (b, t, Z) | 0.4382 | 2.29e-01 | 34.3% | **0.629** | **3.78e-02** |
| (f_pi, tau, c) | 0.4111 | 2.56e-01 | 38.3% | **0.675** | **8.21e-03** |
| (f_pi, b, v_EW) | 0.7698 | 1.03e-01 | 15.5% | 0.824 | 0.16 |

**None of the pure EW boson triples are anywhere near Koide.** Q(W, Z, H) = 0.336, which is essentially 1/3 (the equal-mass limit, unsurprisingly since W, Z, H are within a factor of 1.6 of each other).

(b, t, H) unsigned Q = 0.434 is far from 2/3. In the first scan (which used the wrong formula Q = (sum m)^2/(3 sum m^2)), this appeared close to 2/3 — that was an artifact of a different quantity.

The signed variants of (b, t, W), (b, t, Z), and (b, t, H) are modestly interesting (Q_s = 0.61-0.64), but still 4-9% from 2/3. The signed (f_pi, tau, c) is the notable near-miss at Q_s = 0.675 (1.2% deviation).

---

## 5. Heavy-Sector Koide: Is (c, b, t) Alone?

**Reference**: Q(c, b, t) = 0.66949, |Q - 2/3| = 2.82e-3

All 35 triples from {c, b, t, W, Z, H, v_EW}, sorted by |Q - 2/3|:

| Triple | Q | |Q - 2/3| |
|--------|---|----------|
| **(c, b, t)** | **0.6695** | **2.82e-03** |
| (c, b, H) | 0.6336 | 3.31e-02 |
| (c, b, v_EW) | 0.7073 | 4.07e-02 |
| (c, b, Z) | 0.5972 | 6.95e-02 |
| (c, b, W) | 0.5827 | 8.40e-02 |
| all others | ... | > 0.17 |

**(c, b, t) is completely alone.** The next-closest heavy-sector triple is (c, b, H) at 5.0% deviation — an order of magnitude worse. No triple involving EW bosons comes within a factor of 10 of the (c, b, t) result.

The best signed-Q results from the heavy set:

| Triple | Q_signed | |Q_s - 2/3| |
|--------|----------|-----------|
| (b, W, Z) | 0.6479 | 1.88e-02 |
| (b, W, v_EW) | 0.6469 | 1.97e-02 |
| (b, W, H) | 0.6396 | 2.71e-02 |
| (b, t, W) | 0.6392 | 2.75e-02 |

The signed variants are interesting — the b quark paired with EW masses gives Q_s values around 0.64-0.65, closer to 2/3 than their unsigned values. But these are still 3-4% off, not competitive with the genuine Koide triples.

---

## 6. Koide Parametrization

For triples within 5% of Q = 2/3 (30 triples qualify), we extract M0 and delta from:

m_k = M0 * (1 + sqrt(2) * cos(delta + 2*pi*k/3))^2

where M0 = (m1 + m2 + m3)/6.

| Triple | Q | M0 (GeV) | delta/pi | Nearest frac | Error |
|--------|---|---------|----------|--------------|-------|
| **(e, mu, tau)** | **0.666660** | **0.3138** | **1.2626** | **14/11** | **0.0101** |
| **(c, b, t)** | **0.669489** | **29.7017** | **0.6895** | **7/10** | **0.0105** |
| (b, W, f_pi) | 0.661376 | 14.1068 | 1.2748 | 14/11 | 0.0020 |
| (s, b, W) | 0.661118 | 14.1071 | 1.2749 | 14/11 | 0.0022 |
| (s, b, Z) | 0.674187 | 15.9102 | 0.7265 | 8/11 | 0.0008 |
| (b, Z, f_pi) | 0.674438 | 15.9099 | 0.7266 | 8/11 | 0.0006 |
| (tau, c, Z) | 0.653409 | 15.7057 | 1.3311 | 12/9 | 0.0022 |
| (e, mu, d) | 0.640501 | 0.0185 | 1.2993 | 13/10 | 0.0007 |
| (e, mu, u) | 0.697368 | 0.0181 | 0.7001 | 7/10 | 0.0001 |
| (tau, d, f_pi) | 0.644755 | 0.3123 | 1.2869 | 9/7 | 0.0012 |
| (tau, d, s) | 0.643499 | 0.3125 | 1.2869 | 9/7 | 0.0011 |

The (e, mu, tau) delta = 1.263*pi. This is the well-known delta_0 value.

The (c, b, t) delta = 0.690*pi. Not close to any simple fraction to better than 1%.

The "accidental" near-Koide triples (ranks 3-10) cluster at delta/pi near 14/11 or 8/11. These are not obviously meaningful fractions.

The most notable parametrization result: **(e, mu, u) has delta/pi = 0.7001, within 0.01% of 7/10**. But this triple has Q = 0.697, which is 4.6% from 2/3, so it is not a strong Koide triple.

---

## 7. Pair Scan

Q_pair = (m_1 + m_2) / (sqrt(m_1) + sqrt(m_2))^2

Q_pair = 2/3 requires mass ratio r = (2 + sqrt(3))^2 = 13.928.

| Rank | Pair | Q_pair | |Q - 2/3| | Ratio |
|------|------|--------|----------|-------|
| 1 | **(c, f_pi)** | **0.66581** | **8.59e-04** | **13.804** |
| 2 | **(s, c)** | **0.66435** | **2.31e-03** | **13.597** |
| 3 | (mu, c) | 0.65250 | 1.42e-02 | 12.020 |
| 4 | (mu, tau) | 0.68478 | 1.81e-02 | 16.817 |
| 5 | (tau, s) | 0.69655 | 2.99e-02 | 19.024 |
| 6 | (b, W) | 0.69756 | 3.09e-02 | 19.227 |
| 7 | (tau, f_pi) | 0.69799 | 3.13e-02 | 19.314 |
| 8 | (d, f_pi) | 0.69987 | 3.32e-02 | 19.700 |
| 9 | (d, s) | 0.70130 | 3.46e-02 | 20.000 |
| 10 | (e, d) | 0.62644 | 4.02e-02 | 9.139 |

The pair (c, f_pi) has mass ratio 13.80, close to the Koide pair ratio 13.93. This means m_c / f_pi is close to (2 + sqrt(3))^2. Numerically: m_c = 1.270, f_pi = 0.092, ratio = 13.80 vs target 13.93 (0.9% off).

The pair (s, c) with ratio 13.60 is the second closest.

---

## 8. Mass-Squared Koide: A Different Formula

There is a second quantity that also equals 2/3 for specific mass triples:

**Q_msq = (m_1 + m_2 + m_3)^2 / (3 * (m_1^2 + m_2^2 + m_3^2))**

This is NOT the Koide formula (it has no square roots). Q_msq = 1 for equal masses, Q_msq -> 1/3 for extreme hierarchy. It equals 2/3 when the coefficient of variation of the masses is 1/sqrt(2), i.e., when (sum m)^2 = 2(sum m^2), equivalently sum(m_i * m_j)_{i<j} = (sum m_i^2)/2.

Top results for Q_msq:

| Rank | Triple | Q_msq | |Q - 2/3| | % dev |
|------|--------|-------|----------|-------|
| 1 | (e, mu, f_pi) | 0.66692 | 2.57e-04 | 0.04% |
| 2 | (e, mu, s) | 0.66755 | 8.86e-04 | 0.13% |
| 3 | **(b, t, H)** | **0.66825** | **1.59e-03** | **0.24%** |
| 4 | (mu, W, Z) | 0.66484 | 1.82e-03 | 0.27% |
| 5 | (s, W, Z) | 0.66475 | 1.92e-03 | 0.29% |
| 6 | (W, Z, f_pi) | 0.66474 | 1.93e-03 | 0.29% |
| 7-9 | (x, W, Z) | ~0.664 | ~2.6e-03 | ~0.4% |
| 11 | (b, W, H) | 0.66196 | 4.70e-03 | 0.71% |

In the Q_msq formula, **(b, t, H) is rank 3 at 0.24% deviation** — genuinely close. Note that (e, mu, tau) gives Q_msq = 0.373, far from 2/3 in this formula. **The two formulas pick out completely different triples.**

The algebraic content of Q_msq(b,t,H) = 2/3 is:
- (m_b + m_t + m_H)^2 = 2(m_b^2 + m_t^2 + m_H^2)
- Numerically: (302.19)^2 = 91,339 vs 2*(46,023) = 92,046. Off by 0.8%.
- This says the three EW-heaviest particles (b, t, H) have a specific spread relation.

Whether this is meaningful or accidental is unclear. The (W, Z) pair with ANY light particle (f_pi, mu, s, d, u, e) gives Q_msq near 2/3 because W and Z are close in mass (ratio 1.13), so Q_msq is dominated by (m_W + m_Z)^2 / (3(m_W^2 + m_Z^2)) = 0.6640. Adding any small third mass barely perturbs this.

---

## 9. Summary and Assessment

### The standard Koide formula Q = (sum m)/(sum sqrt m)^2

**Genuine Koide triples (Q within 1% of 2/3)**:
1. **(e, mu, tau)**: Q = 0.66666049, deviation 0.0009%. The original Koide relation. Stands alone at 450x precision above the next competitor.
2. **(c, b, t)**: Q = 0.66949, deviation 0.42%. The only other plausible Koide triple.

**Signed Koide**:
3. **(-s, c, b)**: Q_signed = 0.6750, deviation 1.2%. The known quark Koide chain.

**EW sector results for standard Koide**:
- No triple of EW bosons (W, Z, H) is anywhere near Koide. Q(W,Z,H) = 0.336 ~ 1/3.
- No mixed quark-boson triple reaches 1%.
- (b, t, H) is NOT a Koide triple in the standard formula: Q = 0.434 (35% off).
- The signed variant Q_s(b, t, W) = 0.639 is the best EW candidate, still 4% off.
- **(c, b, t) is completely isolated** in the heavy sector. Nothing involving EW bosons comes within a factor of 10.

### The mass-squared formula Q_msq = (sum m)^2/(3 sum m^2)

**(b, t, H) scores rank 3 at 0.24% deviation.** This IS a striking result, but in a different formula from the original Koide. It is not obvious whether Q_msq has any theoretical motivation comparable to Q_Koide.

The algebraic content is: (m_b + m_t + m_H)^2 ~ 2(m_b^2 + m_t^2 + m_H^2), which constrains the "spread" of the three heaviest particles coupling to the Higgs. Given that m_b << m_t, m_H, this approximately says (m_t + m_H)^2 ~ 2(m_t^2 + m_H^2), i.e., m_t/m_H + m_H/m_t ~ 2sqrt(2), giving m_t/m_H in the range [sqrt(2)-1, sqrt(2)+1]. The actual ratio m_t/m_H = 1.379 sits near sqrt(2) = 1.414. So (b,t,H) Q_msq ~ 2/3 is roughly saying **m_t/m_H ~ sqrt(2)**, which is a mildly interesting numerical relation.

### f_pi connections

Several triples involving f_pi = 0.092 GeV appear in the top 10 standard Koide list (ranks 3, 6, 10, 13), driven by the fact that f_pi creates favorable mass ratios with m_b and m_tau.

The pair (c, f_pi) has the pair-Koide ratio closest to 2/3, meaning **m_c/f_pi ~ (2+sqrt(3))^2 = 13.93** (actual: 13.80, 0.9% off). The pair (s, c) with m_c/m_s = 13.60 is second.

### Hierarchy of significance

For any formula, the question is: what fraction of random triples would score as well? Q_Koide ranges from 1/3 to 1, so a deviation of 6e-6 on a range of 2/3 represents a 9 ppm coincidence (for e, mu, tau). The (c, b, t) deviation of 2.8e-3 is 0.4% — notable but 450x less precise. The (b, t, H) result in Q_msq at 0.24% is comparable in precision to (c, b, t) in the standard formula, but in a formula with no established theoretical interpretation.