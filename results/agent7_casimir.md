# Agent 7: Poincaré Casimir Mass Splitting — Results

## Sign Convention Note

With C₂ = W_μ W^μ = −m²s(s+1), the polynomial M⁴ − M²C₂ + C₁C₂ = 0 becomes:

x² + s(s+1)x − s(s+1) = 0, where x = M²/m²

## Mass Eigenvalues

### s = 1/2

| Root | Exact | Numerical |
|------|-------|-----------|
| M²(1/2,+)/m² | (−3 + √57)/8 | 0.568729304374 |
| M²(1/2,−)/m² | (−3 − √57)/8 | −1.318729304374 |

### s = 1

| Root | Exact | Numerical |
|------|-------|-----------|
| M²(1,+)/m² | √3 − 1 | 0.732050808569 |
| M²(1,−)/m² | −√3 − 1 | −2.732050808569 |

## The Ratio R

$$R = 1 - \frac{M^2_{1/2,+}}{M^2_{1,+}} = \frac{(\sqrt{19}-3)(\sqrt{19}-\sqrt{3})}{16}$$

**R = 0.223101322301** (12 significant digits, verified by script execution)

Note: An earlier analytical estimate gave 0.223058 due to arithmetic error; the code-verified value is 0.223101322301.

### Comparison to Experiment

| Quantity | Value |
|----------|-------|
| R (algebraic) | 0.223101322301 |
| sin²θ_W (on-shell, from M_W/M_Z) | 0.22320 |
| sin²θ_W (PDG quoted) | 0.22306 ± 0.00033 |
| R − sin²θ_W (PDG) | +4.1 × 10⁻⁵ |
| Tension | **0.13σ** |

**Within the 1-sigma experimental error bar.**

## Mass Spectrum (setting M(1,+) = M_Z)

Scale parameter: m = M_Z / √(√3 − 1) = 106.577 GeV

| State | Mass (GeV) | Physical comparison | Agreement |
|-------|------------|-------------------|-----------|
| M(1,+) | 91.1876 | M_Z = 91.1876 GeV | exact (input) |
| M(1/2,+) | 80.374 | M_W = 80.3692 ± 0.0133 GeV | 0.4σ |
| |M(1,−)| | 176.16 | v/√2 = 174.10 GeV | ~1.2% |
| |M(1/2,−)| | 122.39 | — | — |

## Self-Check: Classical de Broglie Derivation

The relativistic standing-wave condition β²/√(1−β²) = √(j(j+1)) yields:

u² + j(j+1)u − j(j+1) = 0

This is **algebraically identical** to the Casimir eigenvalue equation with u = β² ↔ x = M²/m².

| Method | R value |
|--------|---------|
| Casimir (quantum) | 0.223101322301 |
| de Broglie (classical) | 0.223101322301 |
| Difference | **0** (exact, verified symbolically) |

The two derivations are not merely numerically close — they solve the same polynomial.
