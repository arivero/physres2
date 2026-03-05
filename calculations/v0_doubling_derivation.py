#!/usr/bin/env python3
"""
v₀-doubling derivation: homework problem in algebra.

Given the Koide parametrization √m_k = z₀[1 + √2 cos(δ + 2πk/3)], k=0,1,2:
(a) Prove Σ√m_k = 3z₀ for all δ.
(b) Compute seed masses at δ = 3π/4.
(c) Derive √m_b = 3√m_s + √m_c from z₀-doubling with signed square roots.
(d) Numerical verification.
(e) Origin of the factor 3.
"""

import numpy as np
import sympy as sp

print("=" * 72)
print("PART (a): Prove Σ√m_k = 3z₀ for all δ")
print("=" * 72)

delta = sp.Symbol('delta')
z0 = sp.Symbol('z_0', positive=True)

vals = [z0 * (1 + sp.sqrt(2) * sp.cos(delta + sp.Rational(2, 3) * sp.pi * k))
        for k in range(3)]

total = sp.simplify(sum(vals))
print(f"""
  Parametrization: √m_k = z₀[1 + √2 cos(δ + 2πk/3)],  k = 0, 1, 2

  Sum: Σ_k √m_k = Σ_k z₀[1 + √2 cos(δ + 2πk/3)]
                 = 3z₀ + z₀√2 · Σ_k cos(δ + 2πk/3)

  The cosine sum vanishes by the roots-of-unity identity:
    Σ_k cos(δ + 2πk/3) = Re[e^(iδ) · Σ_k ω^k]   where ω = e^(2πi/3)
                        = Re[e^(iδ) · (1 + ω + ω²)]
                        = Re[e^(iδ) · 0]
                        = 0

  Therefore: v₀ ≡ Σ_k √m_k = 3z₀   for ALL δ.  ∎

  SymPy verification: Σ√m_k = {total}
""")

print("=" * 72)
print("PART (b): Seed masses at δ = 3π/4")
print("=" * 72)

d_seed = sp.Rational(3, 4) * sp.pi
seed_coeffs = [sp.simplify(1 + sp.sqrt(2) * sp.cos(d_seed + sp.Rational(2, 3) * sp.pi * k))
               for k in range(3)]

print(f"\n  At δ = 3π/4, the three √m_k/z₀ coefficients are:")
labels = ['m₀', 'm₁', 'm₂']
for k in range(3):
    c = seed_coeffs[k]
    print(f"    √{labels[k]}/z₀ = {c} ≈ {float(c):.10f}")

mass_ratio = sp.simplify((seed_coeffs[2] / seed_coeffs[1]) ** 2)
print(f"""
  Exact values:
    √m₀/z₀ = 0                     [one mass vanishes!]
    √m₁/z₀ = (3 - √3)/2 ≈ 0.634
    √m₂/z₀ = (3 + √3)/2 ≈ 2.366

  Mass ratio: m₂/m₁ = [(3+√3)/(3-√3)]² = {mass_ratio} ≈ {float(mass_ratio):.6f}
  This is the O'Raifeartaigh ratio at gv/m = √3.

  Since √m₁ + √m₂ = 3z₀ (from part a, with √m₀ = 0):
    z₀(seed) = (√m₁ + √m₂)/3

  In terms of the O'R mass parameter m:
    √m₁ = (2 - √3)·m,   √m₂ = (2 + √3)·m
    Sum = 4m   ⟹   z₀(seed) = 4m/3
""")

print("=" * 72)
print("PART (c): Derive √m_b = 3√m_s + √m_c from z₀-doubling")
print("=" * 72)

print("""
  ┌─────────────────────────────────────────────────────────────────┐
  │  KEY INSIGHT: The (-s,c,b) Koide triple uses a SIGNED square   │
  │  root for the strange quark: (-√m_s, +√m_c, +√m_b).           │
  │  This sign is what produces the factor of 3.                   │
  └─────────────────────────────────────────────────────────────────┘

  SEED TRIPLE (at δ = 3π/4, the O'Raifeartaigh point):

    The Koide parametrization gives three signed square roots:
      σ₀√m₀ = z₀[1 + √2 cos(3π/4)]           = 0
      σ₁√m₁ = z₀[1 + √2 cos(3π/4 + 2π/3)]   = z₀(3-√3)/2 > 0
      σ₂√m₂ = z₀[1 + √2 cos(3π/4 + 4π/3)]   = z₀(3+√3)/2 > 0

    All signs are positive (σ_k = +1). The zero mass has σ₀ = +1 trivially.

    Identification with quarks:
      m₀ = 0   →  will bloom into m_s (the lightest in the triple)
      m₁       →  identified with m_c (charm)
      m₂       →  identified with m_b at the seed (NOT physical m_b)

    v₀(seed) = 0 + √m_s^(seed) + √m_c^(seed) = √m_s + √m_c
    z₀(seed) = (√m_s + √m_c)/3


  BLOOM TRIPLE (the physical (-s,c,b) Koide triple):

    When δ shifts away from 3π/4, the zero mass becomes nonzero.
    CRUCIALLY, the lightest mass in the (-s,c,b) triple has a
    NEGATIVE signed square root.  The parametrization gives:

      σ₀√m_s = z₀'[1 + √2 cos(δ')]        < 0  (negative branch)
      σ₁√m_c = z₀'[1 + √2 cos(δ'+2π/3)]   > 0
      σ₂√m_b = z₀'[1 + √2 cos(δ'+4π/3)]   > 0

    Wait — actually the Koide parametrization 1 + √2 cos(θ) is always
    non-negative (minimum = 1 - √2 ≈ -0.414), so individual terms CAN
    be negative. But the convention for the (-s,c,b) triple is that the
    signed sum is:

      v₀(bloom) = -√m_s + √m_c + √m_b = 3z₀'

    The MINUS sign on √m_s is the defining feature of the (-s,c,b) triple.
    It corresponds to the k=0 term having 1 + √2 cos(δ') < 0 for certain
    δ' values, or equivalently, it's the analytic continuation from the
    seed where this term was exactly zero.


  THE DERIVATION:

    z₀-doubling condition: z₀' = 2 · z₀(seed)

    ⟹  v₀(bloom)  =  3z₀'  =  6 · z₀(seed)  =  2 · v₀(seed)

    Substituting:
      v₀(bloom) = -√m_s + √m_c + √m_b
      v₀(seed)  =    0  + √m_s + √m_c     [where 0 is the vanishing mass]

    Note: the m_s and m_c in v₀(seed) are the SAME physical masses as in
    v₀(bloom) — the seed fixes the O'R mass ratio m_c/m_s = 7+4√3, and
    the bloom preserves these two masses while generating m_b.

    Therefore:
      -√m_s + √m_c + √m_b  =  2(√m_s + √m_c)

    Solving for √m_b:

    ┌──────────────────────────────────────────────┐
    │                                              │
    │   √m_b  =  3√m_s + √m_c                     │
    │                                              │
    │   The factor 3 = 2 + 1, where:              │
    │     2 comes from the doubling (z₀' = 2z₀)   │
    │     1 comes from the sign flip (+√m_s → -√m_s)│
    │                                              │
    └──────────────────────────────────────────────┘
""")

print("=" * 72)
print("PART (d): Numerical verification")
print("=" * 72)

m_s = 93.4     # MeV, MS-bar at 2 GeV (PDG 2024)
m_c = 1270.0   # MeV, MS-bar at m_c (PDG 2024)
m_b_pdg = 4180.0  # MeV, MS-bar at m_b (PDG 2024)

sq_s = np.sqrt(m_s)
sq_c = np.sqrt(m_c)
sq_b = np.sqrt(m_b_pdg)

print(f"""
  Input masses (PDG 2024, MS-bar):
    m_s = {m_s} MeV    →  √m_s = {sq_s:.6f}
    m_c = {m_c} MeV   →  √m_c = {sq_c:.6f}
    m_b = {m_b_pdg} MeV   →  √m_b = {sq_b:.6f}
""")

print(f"  O'Raifeartaigh mass ratio check:")
mc_ms_pred = 7 + 4 * np.sqrt(3)
print(f"    m_c/m_s = {m_c / m_s:.4f}")
print(f"    7 + 4√3 = {mc_ms_pred:.4f}")
print(f"    Deviation: {(m_c / m_s / mc_ms_pred - 1) * 100:+.2f}%")

print(f"\n  v₀-DOUBLING PREDICTION:")
mb_pred = (3 * sq_s + sq_c) ** 2
print(f"    √m_b(pred) = 3 × {sq_s:.4f} + {sq_c:.4f}")
print(f"               = {3 * sq_s:.4f} + {sq_c:.4f}")
print(f"               = {3 * sq_s + sq_c:.4f}")
print(f"    √m_b(PDG)  = {sq_b:.4f}")
print(f"    m_b(pred) = {mb_pred:.1f} MeV")
print(f"    m_b(PDG)  = {m_b_pdg:.1f} MeV")
print(f"    Deviation: {(mb_pred / m_b_pdg - 1) * 100:+.3f}%")

print(f"\n  v₀-DOUBLING RATIO CHECK:")
v0_seed = sq_s + sq_c                   # = 0 + √m_s + √m_c
v0_bloom = -sq_s + sq_c + sq_b          # = -√m_s + √m_c + √m_b
print(f"    v₀(seed)  =  0 + √m_s + √m_c        = {v0_seed:.6f}")
print(f"    v₀(bloom) = -√m_s + √m_c + √m_b     = {v0_bloom:.6f}")
print(f"    Ratio v₀(bloom)/v₀(seed)              = {v0_bloom / v0_seed:.6f}")
print(f"    Expected (exact doubling)              = 2.000000")
print(f"    Deviation from 2                       = {(v0_bloom / v0_seed - 2) * 1e4:.2f} × 10⁻⁴")

print(f"\n  KOIDE Q-PARAMETER for signed (-s, c, b) triple:")
Q_signed = (m_s + m_c + m_b_pdg) / v0_bloom ** 2
print(f"    Q = (m_s + m_c + m_b) / (-√m_s + √m_c + √m_b)²")
print(f"      = {m_s + m_c + m_b_pdg:.1f} / {v0_bloom ** 2:.1f}")
print(f"      = {Q_signed:.6f}")
print(f"    2/3 = {2/3:.6f}")
print(f"    Deviation from 2/3: {(Q_signed / (2/3) - 1) * 100:+.3f}%")

print(f"\n  DECOMPOSITION of the factor 3:")
print(f"    √m_b = 3√m_s + √m_c")
print(f"         = (2+1)√m_s + √m_c")
print(f"         = 2√m_s + (√m_s + √m_c)")
print(f"         ────────   ─────────────")
print(f"         from z₀'   = v₀(seed)")
print(f"         = 2z₀      (original two")
print(f"         doubling    nonzero masses)")
print(f"         the sign")
print(f"         flip +→-")
print(f"    ")
print(f"    More precisely:")
print(f"    v₀(bloom)  = -√m_s + √m_c + √m_b = 2 × (√m_s + √m_c) = 2 v₀(seed)")
print(f"    ⟹  √m_b  = 2√m_s + 2√m_c - √m_c + √m_s  ...no, just solve directly:")
print(f"    √m_b  = 2(√m_s + √m_c) - √m_c + √m_s = 3√m_s + √m_c")
print(f"                              └── +√m_s from moving -√m_s to RHS ──┘")

print("\n" + "=" * 72)
print("PART (e): Origin of the factor 3")
print("=" * 72)

print("""
  The integer 3 in √m_b = 3√m_s + √m_c has a CLEAN algebraic origin:

    3 = 2 + 1

  where:
    • The '2' comes from the z₀-DOUBLING condition (z₀' = 2z₀).
      This contributes +2√m_s to √m_b (doubling the seed).

    • The '1' comes from the SIGN FLIP in the (-s,c,b) Koide triple.
      At the seed, √m_s enters v₀(seed) with sign +1.
      At the bloom, √m_s enters v₀(bloom) with sign -1.
      Moving it to the other side of the equation adds another +√m_s.

  Explicitly:  v₀(bloom) = v₀(seed) × 2
              (-√m_s + √m_c + √m_b) = 2(+√m_s + √m_c)
               √m_b = 2√m_s + √m_c + √m_s
               √m_b = 3√m_s + √m_c

  The factor 3 is therefore NOT:
    ✗  N_c = 3 (number of colors)
    ✗  The 3 in v₀ = 3z₀ (sum of three terms)
    ✗  An independent numerical coincidence

  It IS:
    ✓  A consequence of z₀-doubling (factor 2) combined with the
       signed square root convention of the Koide parametrization
       (factor +1 from sign flip).

  The doubling itself (z₀' = 2z₀) still needs a physical explanation.
  The bion Kähler mechanism (δK ~ λ₂|S_bloom − 2·S_seed|²) has been
  proposed as the nonperturbative origin.
""")

print("=" * 72)
print("SUMMARY")
print("=" * 72)
print(f"""
  Inputs:    m_s = {m_s} MeV,  m_c = {m_c} MeV  (PDG 2024 MS-bar)
  Predicts:  m_b = {mb_pred:.1f} MeV
  Observed:  m_b = {m_b_pdg:.1f} MeV  (PDG 2024 MS-bar)
  Accuracy:  {abs(mb_pred / m_b_pdg - 1) * 100:.3f}%  ({abs(mb_pred - m_b_pdg):.1f} MeV)

  The derivation chain:
    O'Raifeartaigh at gv/m = √3
    → seed triple (0, m_s, m_c) with m_c/m_s = 7+4√3
    → z₀-doubling (z₀' = 2z₀, nonperturbative/bion)
    → signed Koide triple (-s, c, b)
    → √m_b = (2+1)√m_s + √m_c = 3√m_s + √m_c
    → m_b ≈ 4177 MeV  ✓
""")
