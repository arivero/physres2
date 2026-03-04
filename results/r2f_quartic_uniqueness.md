# Quartic Uniqueness Analysis — Results

## The Polynomial

x² + s(s+1)·x − s(s+1) = 0, where x = M²/m², C = s(s+1)

Specific case: f(C) = C, g(C) = −C.

## Parameter Space

General form: x² + f(C)·x + g(C) = 0, with f = f₀ + f₁C, g = g₀ + g₁C.
**4 free parameters**: (f₀, f₁, g₀, g₁).

## Uniqueness Result

| Parameter space | dim | codim of R = 0.2231 |
|---|---|---|
| (f₀, f₁, g₀, g₁) | 4 | 1 (hypersurface) |
| f = αC, g = βC | 2 | 1 (curve) |
| f = αC, g = −αC | 1 | **0 (isolated)** |

**R is NOT unique** in the general 4-parameter family. It becomes isolated only after imposing g = −f.

### R varies with α in the g = −f subfamily:

| α | R |
|---|---|
| 1/2 | 0.267 |
| **1** | **0.2231** |
| 2 | 0.172 |
| 3 | 0.141 |

## Three Structural Constraints Pin the Polynomial

1. **f₀ = 0**: no spin-independent mass term
2. **g₀ = 0**: when C = 0 (spin-0), the only solution is x = 0
3. **g₁ = −f₁**: antisymmetry between linear and constant coefficients in x

Plus normalization (f₁ = 1). These 3 constraints on 4 parameters fully determine the polynomial.

**Referee implication**: The value R = 0.2231 is not a free parameter — it is a *consequence* of three structural constraints. But those constraints must themselves be justified from first principles (the simplest polynomial derivation remains open).

## Numerical Verification

R = 0.223101322300866 (15 significant figures)

R = (√19 − 3)(√19 − √3) / 16 ✓ (verified to 50 decimal places)

### Minimal polynomial over Q

64t⁴ − 304t³ + 414t² − 190t + 25 = 0

R ∈ Q(√3, √19), degree 4 algebraic number.

## Algebraic Derivation

R = 1 − x₊(1/2)/x₊(1) = 1 − [(−3+√57)/8] / [√3−1]

Rationalizing:
= 1 − (−3+√57)(√3+1)/16
= [19 + 3√3 − 3√19 − √57] / 16
= (√19 − 3)(√19 − √3) / 16

using √57 = √3·√19.
