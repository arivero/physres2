# Agent 3: Simultaneous 3×2 Grid — Results

## The Q = 3/2 Cone

Setting x_i = s_i √m_i, the condition Q = 3/2 becomes:

x₁² + x₂² + x₃² = 4(x₁x₂ + x₁x₃ + x₂x₃)

Matrix A = [[1,−2,−2],[−2,1,−2],[−2,−2,1]], eigenvalues **(3, 3, −3)** — Lorentzian signature (2,1).

Timelike direction: (1,1,1)/√3. The cone cuts the unit sphere as a circle:

u(θ) = (1,1,1)/√6 + cos(θ)(1,−1,0)/2 + sin(θ)(1,1,−2)/(2√3)

**Emergent O(2,1) symmetry** — the Koide condition defines a Lorentzian light cone.

## Solution Catalog

### TYPE 0: Fully Degenerate (aᵢ = bᵢ)

All 8 triplets identical. Any (m₁,m₂,m₃) satisfying Q=3/2 works.
**Family:** 2-parameter (θ on circle + overall scale M₀).

### TYPE 1: One Row Reflected (Vieta)

Start from (x₁,x₂,x₃) on cone. The cone equation is quadratic in each variable. By Vieta's formulas, the second root for row j is:

**yⱼ = 4(xₖ + xₗ) − xⱼ**

Key relations:
- xⱼ + yⱼ = 4(xₖ + xₗ)
- xⱼ · yⱼ = xₖ² + xₗ² − 4xₖxₗ

**Family:** 2-parameter (θ, M₀) × 3 row choices. Verified numerically.

Example (θ = 20°, row 1 reflected):

| 0.95419 | 0.00022 |
|---------|---------|
| 0.00138 | 0.00138 |
| 0.04443 | 0.04443 |

All 8 Q values = 1.500000 ✓

### TYPE 2: Two Rows Reflected (Discrete)

When two rows are reflected, the fourth triplet (y_i², y_j², x_k²) imposes an additional algebraic constraint:

**S/xₖ = (−9 + √133)/2 ≈ 1.2663**

These exist as a **discrete** (measure-zero) set — parameterized only by M₀ (no continuous angle freedom).

Verified example (rows 1,2 reflected, row 3 kept):

| 0.06897 | 14.0780 |
|---------|---------|
| 0.00001 | 25.4700 |
| 1.00000 | 1.00000 |

All 8 Q = 1.500000 ✓

### TYPE 3: All Three Rows Reflected — NO SOLUTIONS

Exhaustive numerical search (36000 angles + 500 random restarts) found zero Type 3 solutions. The three simultaneous cross-reflection constraints are overdetermined.

## Symmetries

1. **ℝ₊ scaling:** m → λm preserves Q
2. **S₃ row permutations:** reordering rows
3. **ℤ₂ column swap:** exchanging columns a,b
4. **ℤ₂ per-row swap:** exchanging aⱼ,bⱼ (Vieta roots)
5. **U(1) circle rotation:** continuous family for Types 0,1
6. **O(2,1) cone symmetry:** Lorentz group in (2+1)D on signed square roots

**Emergent result:** Vieta reflection yⱼ = 4(xₖ+xₗ) − xⱼ is an involution mapping the cone to itself. Composing two reflections is NOT generally an involution — explaining why Type 2 requires special conditions.

## Physical Scaling

For column-1 masses ≈ (0.122, 1.700, 3.640) GeV, closest exact Q=3/2 point gives M₀ ~ 5.462:

| Row reflected | a (GeV) | b (GeV) | ratio b/a |
|:---:|:---:|:---:|:---:|
| Row 1 (lightest) | 0.122 | 156.2 | 1280 |
| Row 2 (middle) | 1.700 | 59.7 | 35.1 |
| Row 3 (heaviest) | 3.640 | 22.1 | 6.08 |

Reflecting heaviest row gives b₃ ~ 22 GeV (in range of known particles).

## Hint Values k⁺, k⁻

k⁺ = 2+√3, k⁻ = 2−√3 satisfy k⁺k⁻ = 1, k⁺+k⁻ = 4.
However, (k⁻, 1, k⁺) does **not** satisfy Q = 3/2 with any sign choice (best Q ~ 1.17).

## Scripts Produced
- `sources/triplet_final.py` — full numerical verification
- `sources/triplet_catalog.py` — algebraic analysis of Type 2
- `sources/triplet3.py` — initial Vieta reflection derivation
