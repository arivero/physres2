# Koide Phase Analysis

## All Koide Phases (δ in Koide parametrization √m_k = √M₀(1 + √2 cos(2πk/3 + δ)))

| Triple | Q (×3/2) | δ/π | δ mod 2π/3 | M₀ (MeV) |
|--------|----------|-----|------------|----------|
| (e, μ, τ) | 1.500014 | 0.7374 | 0.2222 | 313.84 |
| seed (0, s, c) | exact 3/2 | 0.7500 | π/12 = 0.2618 | — |
| (−s, c, b) | 1.4816 | 0.8734 | 0.6494 | 912.56 |
| (c, b, t) | 1.4937 | 0.6885 | 0.0686 | 29576 |

## Key Relations

### 2/9 coincidence
- δ_lep mod 2π/3 = 0.222230 vs 2/9 = 0.222222
- Residual: 7.4 × 10⁻⁶
- **33 ppm** (agent's "5 ppm" was incorrect — verified independently)
- 3.1σ look-elsewhere corrected, 4.3σ specific

### Bloom phase shift
- Seed δ = 3π/4 = 2.3562
- Full (−s,c,b) δ = 2.7438
- Shift = 0.3877 rad = 0.123π (no obvious simple fraction)

### M₀ scaling
- M₀(−s,c,b) / M₀(lep) = 2.908 (not exactly 3; 3× scaling is approximate)
- δ(−s,c,b) / δ(lep) = 1.184 (NOT 3; the 3× scaling refers to seed, not full)

### Cartan angle
- θ = π/6 − δ (mod 2π), verified numerically to ~10⁻⁷
- No Cartan angles are special multiples of π/6

### Z₉ discrete symmetry
- 2/9 = 2/(3×3) points to Z₉ in Koide modulus space
- Z₃ (mass permutation) is unique subgroup
- CW potential flat in δ at one loop
- Simplest Z₉ potential: V ⊃ −cos(9(δ − 2/9))
- 9 minima, 3 Z₃-equivalent at δ = 2/9 + 2kπ/3

### Coleman-Weinberg
- Minimal O'R CW depends only on |Φ₀|, not phase → δ flat direction
- ISS extension needed with Z₉-breaking terms
