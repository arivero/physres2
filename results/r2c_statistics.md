# Statistical Significance — Results

## Task 1: Look-Elsewhere Corrected p-value for Koide Scan

SM best: |Q − 3/2| = 1.385 × 10⁻⁵ for (e, μ, τ)

Monte Carlo: 10,000 random sets of 9 masses (log-uniform), scan all 84 triplets × 8 signs.

| Mass range | Count ≤ SM best | p-value | σ |
|---|---|---|---|
| [0.1, 200000] MeV | 24/10000 | 0.24% | **2.8σ** |
| [0.5, 180000] MeV | 37/10000 | 0.37% | **2.7σ** |

**The Koide coincidence is ~2.8σ even after look-elsewhere correction.**

## Task 2: δ₀ mod 2π/3 ≈ 2/9

Observed: δ₀ mod 2π/3 = 0.222230, |residual from 2/9| = 7.78 × 10⁻⁶

| Test | p-value | σ |
|---|---|---|
| Near 2/9 specifically | 7.4 × 10⁻⁶ | **4.3σ** |
| Near ANY nice fraction p/q (q≤20) | 9.5 × 10⁻⁴ | **3.1σ** |
| Monte Carlo (near any p/q, q≤20) | 1.08 × 10⁻³ | ~3.1σ |

128 nice fractions p/q with q ≤ 20 exist in [0, 2π/3). Even after this generous look-elsewhere, the proximity to 2/9 is 3.1σ.

The δ₀ ≈ 2/9 coincidence is statistically significant at 3.1σ (look-elsewhere corrected) or 4.3σ (if 2/9 was predicted a priori).

### Note on "7 ppm vs 33 ppm"
The residual is 7.78 × 10⁻⁶.
- Relative to 2/9 = 0.222222: 7.78e-6 / 0.222222 = 35 ppm
- Relative to 2π/3 = 2.0944: 7.78e-6 / 2.0944 = 3.7 ppm
- The "7 ppm" claim in the initial computation used an intermediate denominator. **The fair statement is 35 ppm relative to 2/9.**

## Task 3: SU(5) Anomaly and Asymptotic Freedom

### Anomaly coefficients
| Rep | A(R) |
|-----|------|
| 24 (adjoint, real) | 0 |
| 15 (symmetric) | N+4 = 9 |
| 15̄ | −9 |
| **Total** | **0 ✓** |

24 ⊕ 15 ⊕ 15̄ is **anomaly-free**.

### Dynkin indices and beta function
| Rep | T(R) |
|-----|------|
| 24 | 5 |
| 15 | 7/2 |
| 15̄ | 7/2 |
| **Total** | **12** |

One-loop beta: b₀ = (11/3)×5 − (2/3)×12 = 55/3 − 8 = **31/3 ≈ 10.33**

**SU(5) with 24 ⊕ 15 ⊕ 15̄ as fermions is asymptotically free.** (b₀ > 0)
