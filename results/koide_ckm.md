# Mixing Matrix from Overlapping Koide Mass Triples

**Computation date:** 2026-03-04

## Setup and Notation

The Koide parametrisation assigns a signed square-root value to each slot k:

    sv_k = √M₀ · (1 + √2 cos(2πk/3 + δ)),  k = 0, 1, 2

Algebraic identities (exact for any M₀ > 0, any δ, any permutation of k):

    Σ_k sv_k = 3√M₀         →  M₀ = (Σ sv_k)² / 9
    Σ_k sv_k² = 6M₀
    Q_signed ≡ Σsv² / (Σsv)² = 6M₀/9M₀ = 2/3

The invariant that equals 2/3 is **Q_signed** (using the signed sum Σsv).
For triples with a negative entry (like (−√m_s, √m_c, √m_b)),
Q_signed = 2/3 is the appropriate statement; Q_unsigned uses Σ|sv| and
equals 2/3 only when all entries are non-negative.

Each physical triple is described by two free parameters (M₀, δ),
with M₀ determined analytically from the signed sum. The best-fit δ
and the permutation mapping physical indices to Koide slots k are
found by minimising the fit residual over all 6 permutations.

Mass values (MeV): m_u = 2.16, m_d = 4.67, m_s = 93.4,
m_c = 1270, m_b = 4180, m_t = 172760.

## Part 1: Koide Parameters

### Triple A: (−√m_s, √m_c, √m_b)

Best-fit k-permutation (physical index → Koide slot): (0, 1, 2)

| Parameter | Value |
|-----------|-------|
| M₀ (MeV) | 912.555714 |
| √M₀ (MeV^1/2) | 30.208537 |
| δ (rad) | 2.74384327 |
| δ / (2π/3) | 1.310089 |
| δ (deg) | 157.2106 |
| Q_signed = Σsv²/(Σsv)² | 0.6666666667 |
| Q_unsigned = Σsv²/(Σ\|sv\|)² | 0.461008 |
| \|Q_signed − 2/3\| | 2.22e-16 |
| Max fit residual (MeV^1/2) | 4.866e-01 |

Fitted sv:  (-9.17775, 35.57081, 64.23255)
Target sv:  (-9.66437, 35.63706, 64.65292)
Residuals:  (-4.866e-01, 6.625e-02, 4.204e-01)

### Triple B: (√m_c, √m_b, √m_t)

Best-fit k-permutation (physical index → Koide slot): (0, 1, 2)

| Parameter | Value |
|-----------|-------|
| M₀ (MeV) | 29576.439062 |
| √M₀ (MeV^1/2) | 171.978019 |
| δ (rad) | 2.16303820 |
| δ / (2π/3) | 1.032775 |
| δ (deg) | 123.9330 |
| Q_signed = Σsv²/(Σsv)² | 0.6666666667 |
| Q_unsigned = Σsv²/(Σ\|sv\|)² | 0.666667 |
| \|Q_signed − 2/3\| | 1.11e-16 |
| Max fit residual (MeV^1/2) | 1.025e+00 |

Fitted sv:  (36.21069, 65.10447, 414.61890)
Target sv:  (35.63706, 64.65292, 415.64408)
Residuals:  (-5.736e-01, -4.516e-01, 1.025e+00)

Phase difference: δ_B − δ_A = -33.2777°

**Note on fit residuals.** The parametrisation has two free parameters (M₀, δ)
for three mass values, so an exact fit is not generically possible.
The residuals (~0.5–1 MeV^1/2) reflect how far the PDG masses lie from
the nearest Koide curve. Q_signed = 2/3 exactly by algebraic identity,
regardless of fit quality.

## Part 2: Rotation Matrix R

The 3D flavor space decomposes into orthogonal subspaces:

- **Democratic direction** e₀ = (1,1,1)/√3: carries the Koide scale M₀
  (invariant under δ-variation).
- **Cartan plane**: spanned by f₁ = (1,−1,0)/√2 (λ₃) and
  f₂ = (1,1,−2)/√6 (λ₈).

Within the Cartan plane the Koide direction vector rotates as δ changes.
The rotation from Triple A to Triple B is determined by:

    Δφ = φ_B − φ_A = 33.2777°   (mod 360°, in (−180°, +180°])

2D rotation matrix in Cartan plane (angle Δφ = 33.2777°):

```
R_2D = [[+0.836021  -0.548697]
        [+0.548697  +0.836021]]
```

Full 3D rotation (identity on e₀, R_2D on Cartan plane):

```
R_3D =
  [+0.890681  -0.262131  +0.371450]
  [+0.371450  +0.890681  -0.262131]
  [-0.262131  +0.371450  +0.890681]
```

det(R_3D) = 1.00000000  (proper rotation)

## Part 3: CKM-like Angles from R_3D

PDG parametrisation applied to R_3D:

| Angle | Value (deg) | sin |
|-------|-------------|-----|
| θ₁₂ | 16.3994 | 0.282331 |
| θ₂₃ | 16.3994 | 0.282331 |
| θ₁₃ | 21.8051 | 0.371450 |
| δ_CP | -0.0000 | — |

| Element | From R_3D | PDG | Ratio |
|---------|-----------|-----|-------|
| \|V_us\| | 0.262131 | 0.2243 | 1.1687 |
| \|V_cb\| | 0.262131 | 0.0422 | 6.2116 |
| \|V_ub\| | 0.371450 | 0.00394 | 94.2767 |

The Koide-basis rotation does not reproduce CKM magnitudes.
The Cartan rotation angle is Δφ = 33.3°,
which is 2.6× the Cabibbo angle.

## Part 4: Koide Direction Angles

| Quantity | Value | Cabibbo reference | Ratio |
|---------|-------|-------------------|-------|
| 3D angle (n̂_A · n̂_B) | 23.3630° | 13.00° | 1.7972 |
| Cartan \|Δφ\| | 33.2777° | 13.00° | 2.5598 |

The angle between the two Koide direction vectors is 23.36°,
approximately 1.8× the Cabibbo angle.

## Part 5: Oakes Relation

| Element | Formula | Computed | PDG | Deviation |
|---------|---------|----------|-----|-----------|
| \|V_us\| | √(m_d/m_s) | 0.223607 | 0.2243 | -0.3% |
| \|V_cb\| | √(m_s/m_b) | 0.149481 | 0.0422 | +254.2% |
| \|V_ub\| | √(m_d/m_b) | 0.033425 | 0.00394 | +748.3% |

- √(m_d/m_s) = 0.22361: reproduces |V_us| = 0.2243 at −0.3%. **Works well.**
- √(m_s/m_b) = 0.14948: vs |V_cb| = 0.0422, off by factor 3.5. **Fails.**
- √(m_d/m_b) = 0.03342: vs |V_ub| = 0.00394, off by factor 8.5. **Fails.**

## Part 6: Wolfenstein-like Expansion

Cabibbo parameter: λ = √(m_d/m_s) = 0.223607  (PDG |V_us| = 0.2243, ratio 0.9969)
λ² = 0.050000,  λ³ = 0.011180

| Check | Computed | λⁿ | Ratio |
|-------|----------|-----|-------|
| √(m_s/m_b) vs λ² | 0.149481 | 0.050000 | 2.9896 |
| √(m_d/m_b) vs λ³ | 0.033425 | 0.011180 | 2.9896 |

Both ratios are equal (algebraic identity: √(m_d/m_b)/λ³ = √(m_s/m_b)/λ²).
This common ratio is the Wolfenstein A parameter implied by masses.

| Wolfenstein parameter | From mass ratios | From PDG CKM | Ratio |
|----------------------|-----------------|--------------|-------|
| λ | 0.223607 | 0.224300 | 0.9969 |
| A | 2.9896 | 0.8388 | 3.5642 |
| \|ρ−iη\| | 1.0000 | 0.4163 | 2.4024 |

The λ parameter matches PDG at 0.3%. The Wolfenstein A from masses (2.990)
exceeds the PDG value (0.839) by a factor of 3.56.
|ρ−iη| from masses (1.0000) vs PDG (0.4163): ratio 2.4024.

## Summary Table

| Quantity | Computed | Reference | Ratio |
|---------|---------|-----------|-------|
| |V_us| Oakes | 0.22361 | 0.2243 [PDG] | 0.9969 |
| |V_cb| Oakes | 0.14948 | 0.0422 [PDG] | 3.5422 |
| |V_ub| Oakes | 0.033425 | 0.00394 [PDG] | 8.4835 |
| Wolfenstein A | 2.9896 | 0.83879 [PDG] | 3.5642 |
| Wolfenstein |ρ−iη| | 1 | 0.41625 [PDG] | 2.4024 |
| |V_us| from R_3D | 0.26213 | 0.2243 [PDG] | 1.1687 |
| |V_cb| from R_3D | 0.26213 | 0.0422 [PDG] | 6.2116 |
| |V_ub| from R_3D | 0.37145 | 0.00394 [PDG] | 94.2767 |
| 3D Koide angle (deg) | 23.363 | 13 [Cabibbo] | 1.7972 |
| Cartan |Δφ| (deg) | 33.278 | 13 [Cabibbo] | 2.5598 |

## Conclusions

1. **Koide parameters.** Both triples fitted successfully.
   M₀_A = 912.556 MeV, δ_A = 157.211°;
   M₀_B = 29576.439 MeV, δ_B = 123.933°;
   Δδ = -33.278°.
   Q_signed = 2/3 exactly by algebraic identity for both triples.
   Fit residuals (~0.49–1.03 MeV^1/2) measure how far PDG masses
   deviate from any Koide curve.

2. **Oakes relation.** √(m_d/m_s) = 0.2236 matches |V_us| = 0.2243
   at −0.3%. The generalisation fails for |V_cb| and |V_ub|,
   which are off by factors of 3.5 and 8.5.

3. **Wolfenstein power structure.** The hierarchy λ : λ² : λ³ is present
   in the mass ratios with λ = √(m_d/m_s) ≈ 0.2236 (PDG: 0.2243).
   The Wolfenstein A from masses (2.990) is 3.56× the PDG value (0.839).
   The mass-ratio ansatz captures the power-law hierarchy but not the scale.

4. **Koide-basis rotation R_3D.** Geometrically clean (det = +1),
   rotation angle Δφ = 33.3° ≈ 2.6× θ_C.
   CKM-like elements from R_3D are 1.2–94× too large.

5. **Overall.** The only precise mass-to-mixing connection found here
   is the Oakes/Cabibbo relation |V_us| ≈ √(m_d/m_s) (−0.3%).
   No construction based solely on Koide phase parameters reproduces
   |V_cb| or |V_ub| at the correct magnitude.