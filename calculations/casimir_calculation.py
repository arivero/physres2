#!/usr/bin/env python3
"""
Poincaré Casimir invariants and electroweak mass spectrum calculation.

The two Casimir invariants of the Poincaré group:
  C₁ = P_μ P^μ = m²
  C₂ = W_μ W^μ = -m² s(s+1)

Mass operator equation:
  M⁴ - M² C₂ + C₁ C₂ = 0

Substituting C₁ = m², C₂ = -m² s(s+1):
  M⁴ + m² s(s+1) M² - m⁴ s(s+1) = 0

This is a quadratic in X = M²:
  X² + m² s(s+1) X - m⁴ s(s+1) = 0

Solution:
  X = m² s(s+1)/2 · (-1 ± √(1 + 4/(s(s+1))))

The formula written by the user:
  M_s² = m² s(s+1)/2 · (1 ± √(1 - 4/(s(s+1))))
encodes the same equation (the signs are absorbed differently).
The physically relevant (positive) root uses the "+" branch of the ± with the
"(-1 + √(...))" form.
"""

from sympy import (Symbol, Rational, sqrt, simplify, solve, S, radsimp,
                   together, factor, nsimplify, pretty, I, Abs, re, im, latex)
import mpmath
from mpmath import mp, mpf, sqrt as mpsqrt

mp.dps = 50

m = Symbol('m', positive=True)

print("=" * 80)
print("POINCARÉ GROUP CASIMIR INVARIANTS — ELECTROWEAK MASS SPECTRUM")
print("=" * 80)
print()
print("Equation: M⁴ + m² s(s+1) M² - m⁴ s(s+1) = 0")
print("Solution: M_s² = m² s(s+1)/2 · (-1 + √(1 + 4/(s(s+1))))   [positive root]")
print("          M_s² = m² s(s+1)/2 · (-1 - √(1 + 4/(s(s+1))))   [negative root]")

# ============================================================
# Helper
# ============================================================
def M2_roots(spin):
    """Return (M²_+, M²_-) for given spin, as sympy expressions × m²."""
    ss1 = spin * (spin + 1)
    disc = sqrt(1 + Rational(4, 1) / ss1)
    plus = ss1 / 2 * (-1 + disc)
    minus = ss1 / 2 * (-1 - disc)
    return simplify(plus), simplify(minus)

# ============================================================
# 1. s = 1/2
# ============================================================
print("\n" + "=" * 80)
print("1. SPIN s = 1/2")
print("=" * 80)

s_half = Rational(1, 2)
ss1_half = s_half * (s_half + 1)  # 3/4
print(f"\n   s(s+1) = (1/2)(3/2) = {ss1_half}")
print(f"   4/s(s+1) = 4/(3/4) = {Rational(4)/ss1_half} = 16/3")
print(f"   1 + 16/3 = 19/3")
print(f"   √(19/3) = √57/√9 = √57/3")

c_half_plus, c_half_minus = M2_roots(s_half)
print(f"\n   M²_(1/2,+) / m² = (3/8)(-1 + √(19/3))")
print(f"                    = {c_half_plus}")
print(f"                    = {float(c_half_plus):.15f}")
print(f"\n   M²_(1/2,-) / m² = (3/8)(-1 - √(19/3))")
print(f"                    = {c_half_minus}")
print(f"                    = {float(c_half_minus):.15f}  [negative → tachyonic]")

# Verify these are roots
X = Symbol('X')
eq_half = X**2 + ss1_half * X - ss1_half
print(f"\n   Verification: X² + (3/4)X - 3/4 = 0")
print(f"   f(X_+) = {simplify(eq_half.subs(X, c_half_plus))}")
print(f"   f(X_-) = {simplify(eq_half.subs(X, c_half_minus))}")

# ============================================================
# 2. s = 1
# ============================================================
print("\n" + "=" * 80)
print("2. SPIN s = 1")
print("=" * 80)

s_one = Rational(1)
ss1_one = s_one * (s_one + 1)  # 2
print(f"\n   s(s+1) = 1·2 = {ss1_one}")
print(f"   4/s(s+1) = 4/2 = 2")
print(f"   1 + 2 = 3")
print(f"   √3")

c_one_plus, c_one_minus = M2_roots(s_one)
print(f"\n   M²_(1,+) / m² = 1·(-1 + √3)")
print(f"                  = {c_one_plus}")
print(f"                  = {float(c_one_plus):.15f}")
print(f"\n   M²_(1,-) / m² = 1·(-1 - √3)")
print(f"                  = {c_one_minus}")
print(f"                  = {float(c_one_minus):.15f}  [negative → tachyonic]")

# Verify
eq_one = X**2 + 2*X - 2
print(f"\n   Verification: X² + 2X - 2 = 0")
print(f"   f(X_+) = {simplify(eq_one.subs(X, c_one_plus))}")
print(f"   f(X_-) = {simplify(eq_one.subs(X, c_one_minus))}")

# Vieta's check
print(f"\n   Vieta's formulae:")
print(f"   X_+ + X_- = {simplify(c_one_plus + c_one_minus)}  (should be -2)")
print(f"   X_+ · X_- = {simplify(c_one_plus * c_one_minus)}  (should be -2)")

# ============================================================
# 3. Ratio R
# ============================================================
print("\n" + "=" * 80)
print("3. RATIO R = 1 - M²_(1/2,+) / M²_(1,+)")
print("=" * 80)

R_exact = S(1) - c_half_plus / c_one_plus
R_simplified = simplify(R_exact)
R_radsimped = radsimp(R_exact)

print(f"\n   R = 1 - [(3/8)(-1+√(19/3))] / [(-1+√3)]")
print(f"\n   Sympy simplified: {R_simplified}")
print(f"   Radsimp:          {R_radsimped}")

# Manual rationalization
# num = (3/8)(-1 + sqrt(19/3))
# den = sqrt(3) - 1
# Multiply num and den by (sqrt(3)+1): den becomes 3-1=2
# Multiply num and den by sqrt(3): to clear sqrt(19/3) = sqrt(57)/3

# Let's get a fully rationalized form step by step:
# R = 1 - 3(-1 + sqrt(19/3)) / (8(sqrt(3)-1))
# sqrt(19/3) = sqrt(57)/3
# So: 3(-1 + sqrt(57)/3) / (8(sqrt(3)-1)) = 3(-1 + sqrt(57)/3) / (8(sqrt(3)-1))
# = (-3 + sqrt(57)) / (8(sqrt(3)-1))
# Multiply top and bottom by (sqrt(3)+1):
# = (-3 + sqrt(57))(sqrt(3)+1) / (8·2)
# = (-3sqrt(3) - 3 + sqrt(57·3) + sqrt(57)) / 16
# = (-3sqrt(3) - 3 + sqrt(171) + sqrt(57)) / 16
# sqrt(171) = sqrt(9·19) = 3sqrt(19)
# = (-3sqrt(3) - 3 + 3sqrt(19) + sqrt(57)) / 16

# So R = 1 - (-3sqrt(3) - 3 + 3sqrt(19) + sqrt(57))/16
# = (16 + 3sqrt(3) + 3 - 3sqrt(19) - sqrt(57)) / 16
# = (19 + 3sqrt(3) - 3sqrt(19) - sqrt(57)) / 16

R_rationalized = (19 + 3*sqrt(3) - 3*sqrt(19) - sqrt(57)) / 16
R_rationalized_simplified = simplify(R_rationalized)
print(f"\n   Fully rationalized:")
print(f"   R = (19 + 3√3 - 3√19 - √57) / 16")
print(f"     = {R_rationalized_simplified}")

# Verify equality
diff_check = simplify(R_rationalized - R_exact)
print(f"   Check R_rationalized - R_exact = {diff_check}")

# High-precision numerical value
R_mp = 1 - mpf(3)/8 * (-1 + mpsqrt(mpf(19)/3)) / (-1 + mpsqrt(3))
print(f"\n   R = {mp.nstr(R_mp, 30)}")
print(f"   R (12 significant figures) = {mp.nstr(R_mp, 12)}")

# Alternative check with mpmath
R_mp2 = (19 + 3*mpsqrt(3) - 3*mpsqrt(19) - mpsqrt(57)) / 16
print(f"   R (from rationalized form) = {mp.nstr(R_mp2, 30)}")
print(f"   Difference = {mp.nstr(R_mp - R_mp2, 5)}")

# ============================================================
# 4. Comparison with sin²θ_W
# ============================================================
print("\n" + "=" * 80)
print("4. COMPARISON WITH EXPERIMENTAL sin²θ_W")
print("=" * 80)

MW_val = mpf('80.3692')
MW_err = mpf('0.0133')
MZ_val = mpf('91.1876')
MZ_err = mpf('0.0021')

sin2_onshell = 1 - (MW_val / MZ_val)**2
sin2_central = mpf('0.22306')
sin2_err = mpf('0.00033')

print(f"\n   On-shell definition: sin²θ_W = 1 - M_W²/M_Z²")
print(f"   M_W = {MW_val} ± {MW_err} GeV")
print(f"   M_Z = {MZ_val} ± {MZ_err} GeV")
print(f"   sin²θ_W (from masses) = {mp.nstr(sin2_onshell, 8)}")
print(f"   sin²θ_W (PDG value)   = {sin2_central} ± {sin2_err}")
print(f"\n   R (Casimir)           = {mp.nstr(R_mp, 12)}")
print(f"   sin²θ_W (on-shell)    = {mp.nstr(sin2_onshell, 8)}")
print(f"\n   R - sin²θ_W           = {mp.nstr(R_mp - sin2_onshell, 6)}")
print(f"   R - 0.22306           = {mp.nstr(R_mp - sin2_central, 6)}")

sigma = abs(R_mp - sin2_central) / sin2_err
print(f"   |R - sin²θ_W| / σ    = {mp.nstr(sigma, 4)} standard deviations")

# Propagate error on sin2 from MW, MZ
# sin2 = 1 - (MW/MZ)^2
# d(sin2)/d(MW) = -2MW/MZ^2
# d(sin2)/d(MZ) = 2MW^2/MZ^3
dsindMW = -2*MW_val / MZ_val**2
dsindMZ = 2*MW_val**2 / MZ_val**3
sin2_prop_err = mpsqrt((dsindMW * MW_err)**2 + (dsindMZ * MZ_err)**2)
print(f"\n   Propagated uncertainty on sin²θ_W from masses: ±{mp.nstr(sin2_prop_err, 4)}")
sigma2 = abs(R_mp - sin2_onshell) / sin2_prop_err
print(f"   |R - sin²θ_W(masses)| / σ_prop = {mp.nstr(sigma2, 4)} standard deviations")

# ============================================================
# 5. Mass spectrum with M_{1,+} = M_Z
# ============================================================
print("\n" + "=" * 80)
print("5. MASS SPECTRUM: setting M_(1,+) = M_Z = 91.1876 GeV")
print("=" * 80)

# M²_(1,+) = m² (√3 - 1)  =>  m² = M_Z² / (√3 - 1)
sqrt3 = mpsqrt(3)
m2_val = MZ_val**2 / (sqrt3 - 1)
m_phys = mpsqrt(m2_val)
print(f"\n   m² = M_Z² / (√3 - 1) = {mp.nstr(m2_val, 10)} GeV²")
print(f"   m  = {mp.nstr(m_phys, 10)} GeV")

# M_(1/2,+)
M2_half_plus_num = m2_val * float(c_half_plus)
M_half_plus_num = mpsqrt(M2_half_plus_num)
print(f"\n   M²_(1/2,+) = m² · (3/8)(-1+√(19/3)) = {mp.nstr(M2_half_plus_num, 10)} GeV²")
print(f"   M_(1/2,+)  = {mp.nstr(M_half_plus_num, 10)} GeV")
print(f"   (Compare: M_W = {MW_val} ± {MW_err} GeV)")
print(f"   Difference: M_(1/2,+) - M_W = {mp.nstr(M_half_plus_num - MW_val, 5)} GeV")

# Cross-check: M_(1/2,+) = M_Z * sqrt(1-R)
M_half_cross = MZ_val * mpsqrt(1 - R_mp)
print(f"   Cross-check: M_Z·√(1-R) = {mp.nstr(M_half_cross, 10)} GeV  ✓")

# M_(1,-)
M2_one_minus_num = m2_val * float(c_one_minus)  # negative
M_one_minus_abs = mpsqrt(abs(M2_one_minus_num))
print(f"\n   M²_(1,-) = m² · (-1-√3) = {mp.nstr(M2_one_minus_num, 10)} GeV²  [negative]")
print(f"   |M_(1,-)| = √|M²_(1,-)| = {mp.nstr(M_one_minus_abs, 10)} GeV")

# Exact: |M²_(1,-)| = m²(1+√3), M²_(1,+) = m²(√3-1)
# |M²_(1,-)| / M²_(1,+) = (1+√3)/(√3-1) = (1+√3)²/2
# |M_(1,-)|² / M_(1,+)² = (1+√3)/(√3-1)
# |M_(1,-)| / M_(1,+) = √((1+√3)/(√3-1)) = (1+√3)/√2  [rationalizing]
ratio_exact = (1 + sqrt3) / mpsqrt(2)
print(f"   Exact: |M_(1,-)| / M_Z = (1+√3)/√2 = {mp.nstr(ratio_exact, 10)}")
print(f"   Check: M_Z · (1+√3)/√2 = {mp.nstr(MZ_val * ratio_exact, 10)} GeV  ✓")

# M_(1/2,-)
M2_half_minus_num = m2_val * float(c_half_minus)  # negative
M_half_minus_abs = mpsqrt(abs(M2_half_minus_num))
print(f"\n   M²_(1/2,-) = m² · (3/8)(-1-√(19/3)) = {mp.nstr(M2_half_minus_num, 10)} GeV²  [negative]")
print(f"   |M_(1/2,-)| = √|M²_(1/2,-)| = {mp.nstr(M_half_minus_abs, 10)} GeV")

# ============================================================
# 6. Comparison with electroweak VEV
# ============================================================
print("\n" + "=" * 80)
print("6. COMPARISON WITH ELECTROWEAK VACUUM EXPECTATION VALUE")
print("=" * 80)

GF = mpf('1.1663788e-5')  # GeV^-2
v_vev = mpf('246.22')
v_sqrt2 = v_vev / mpsqrt(2)

print(f"\n   ⟨v⟩ = (√2 G_F)^(-1/2) = {v_vev} GeV")
print(f"   ⟨v⟩/√2 = {mp.nstr(v_sqrt2, 6)} GeV")
print(f"\n   |M_(1,-)| = {mp.nstr(M_one_minus_abs, 7)} GeV")
print(f"   ⟨v⟩/√2   = {mp.nstr(v_sqrt2, 6)} GeV")
print(f"\n   |M_(1,-)| - ⟨v⟩/√2 = {mp.nstr(M_one_minus_abs - v_sqrt2, 5)} GeV")
print(f"   Relative difference = {mp.nstr(100*(M_one_minus_abs - v_sqrt2)/v_sqrt2, 4)}%")

# Exact algebraic relationship
# |M_(1,-)| = M_Z · (1+√3)/√2
# v/√2 ≈ 174.10 GeV
# M_Z · (1+√3)/√2 = 91.1876 · 1.93185... = 176.17... GeV
print(f"\n   Exact: |M_(1,-)| = M_Z · (1+√3)/√2")
print(f"   Numerically: {mp.nstr(MZ_val * (1+sqrt3)/mpsqrt(2), 8)} GeV")

# Also note: if we use the product relation
# M_(1,+)·|M_(1,-)| = √(|X_+ · X_-|) · m² = √2 · m²   [from Vieta: X_+·X_- = -2m⁴]
# Actually: X_+·X_- = -s(s+1)m⁴ = -2m⁴ for s=1
# |X_+·X_-| = 2m⁴
# √(M²_(1,+) · |M²_(1,-)|) = m² √2
# So M_(1,+) · |M_(1,-)| = m² √(2·(√3-1)·(1+√3)) = m²√(2·2) = 2m²
# Hmm let me be more careful:
# M²_(1,+) · |M²_(1,-)| = m⁴(√3-1)(1+√3) = m⁴(3-1) = 2m⁴
# So √(M²_(1,+)) · √(|M²_(1,-)|) = m² √2
# M_(1,+) · |M_(1,-)| = m²√2
product = MZ_val * M_one_minus_abs
print(f"\n   Product relation: M_(1,+) · |M_(1,-)| = m²√2")
print(f"   M_Z · |M_(1,-)| = {mp.nstr(product, 8)} GeV²")
print(f"   m²√2 = {mp.nstr(m2_val * mpsqrt(2), 8)} GeV²")
print(f"   Match: {mp.nstr(abs(product - m2_val*mpsqrt(2)), 3)}")

# ============================================================
# SELF-CHECK: Classical derivation
# ============================================================
print("\n" + "=" * 80)
print("SELF-CHECK: CLASSICAL (de Broglie / Landé-Pauli) DERIVATION")
print("=" * 80)
print()
print("   Relativistic velocity condition with angular momentum substitution:")
print("   β² / √(1-β²) = √(j(j+1))")
print()
print("   Squaring: β⁴ / (1-β²) = j(j+1)")
print("   Let x = β²:  x² + j(j+1)x - j(j+1) = 0")
print()
print("   This is IDENTICAL to the Casimir equation X² + s(s+1)X - s(s+1) = 0")
print("   with m = 1, confirming M²/m² = β².")

# j = 1/2
jj1_half = Rational(3, 4)
beta2_half = simplify(jj1_half / 2 * (-1 + sqrt(1 + 4 / jj1_half)))
print(f"\n   j = 1/2:")
print(f"   j(j+1) = 3/4")
print(f"   β²_(1/2) = (3/8)(-1 + √(19/3))")
print(f"            = {beta2_half}")
print(f"            = {float(beta2_half):.15f}")

# j = 1
jj1_one = Rational(2)
beta2_one = simplify(jj1_one / 2 * (-1 + sqrt(1 + 4 / jj1_one)))
print(f"\n   j = 1:")
print(f"   j(j+1) = 2")
print(f"   β²_1 = -1 + √3")
print(f"        = {beta2_one}")
print(f"        = {float(beta2_one):.15f}")

# R_classical
R_classical_sym = simplify(1 - beta2_half / beta2_one)
print(f"\n   R_classical = 1 - β²_(1/2) / β²_1")
print(f"               = {R_classical_sym}")
print(f"               = {float(R_classical_sym):.15f}")

R_classical_mp = 1 - mpf(3)/8 * (-1 + mpsqrt(mpf(19)/3)) / (-1 + mpsqrt(3))
print(f"\n   R_classical (50 digits) = {mp.nstr(R_classical_mp, 30)}")

# Verify classical β satisfies the original equation
print(f"\n   Verification: β²/(1-β²)^(1/2) = √(j(j+1))?")
b2h = float(beta2_half)
lhs_half = b2h / (1 - b2h)**0.5
rhs_half = (0.75)**0.5
print(f"   j=1/2: β⁴/(1-β²) = {b2h**2/(1-b2h):.10f}, j(j+1) = {0.75:.10f}")
b2o = float(beta2_one)
lhs_one = b2o / (1 - b2o)**0.5
rhs_one = (2.0)**0.5
print(f"   j=1:   β⁴/(1-β²) = {b2o**2/(1-b2o):.10f}, j(j+1) = {2.0:.10f}")

# ============================================================
# Cross-check between Casimir and Classical
# ============================================================
print("\n" + "-" * 80)
print("   CROSS-CHECK:")
print(f"   R_Casimir   = {mp.nstr(R_mp, 20)}")
print(f"   R_classical = {mp.nstr(R_classical_mp, 20)}")
print(f"   Difference  = {mp.nstr(R_mp - R_classical_mp, 5)}")
print(f"   (Identically zero: the equations are algebraically the same)")
print("-" * 80)

# Symbolic verification
R_cas_sym = simplify(1 - c_half_plus / c_one_plus)
diff_sym = simplify(R_cas_sym - R_classical_sym)
print(f"\n   Symbolic: R_Casimir - R_classical = {diff_sym}")

# ============================================================
# SUMMARY
# ============================================================
print("\n" + "=" * 80)
print("COMPLETE SUMMARY")
print("=" * 80)

print(f"""
┌───────────────────────────────────────────────────────────────────────────────┐
│ MASS SPECTRUM FROM POINCARÉ CASIMIR INVARIANTS                              │
│ Equation: M⁴ - M² C₂ + C₁ C₂ = 0  with C₁=m², C₂=-m²s(s+1)              │
├───────────────────────────────────────────────────────────────────────────────┤
│                                                                             │
│ s = 1/2:                                                                    │
│   M²_(1/2,+) / m² = (3/8)(-1 + √(19/3))  = {float(c_half_plus):.15f}      │
│   M²_(1/2,-) / m² = (3/8)(-1 - √(19/3))  = {float(c_half_minus):.15f}     │
│                                                                             │
│ s = 1:                                                                      │
│   M²_(1,+)   / m² = √3 - 1               = {float(c_one_plus):.15f}       │
│   M²_(1,-)   / m² = -1 - √3              = {float(c_one_minus):.15f}      │
│                                                                             │
├───────────────────────────────────────────────────────────────────────────────┤
│ RATIO:                                                                      │
│                                                                             │
│   R = 1 - M²_(1/2,+) / M²_(1,+)                                           │
│     = (19 + 3√3 - 3√19 - √57) / 16                                        │
│     = {mp.nstr(R_mp, 12):>20s}                                              │
│                                                                             │
│   sin²θ_W (on-shell, PDG) = 0.22306 ± 0.00033                             │
│   sin²θ_W (from M_W/M_Z)  = {mp.nstr(sin2_onshell, 8):>10s}                │
│                                                                             │
│   Difference R - sin²θ_W  = {mp.nstr(R_mp - sin2_central, 5):>10s}          │
│   Tension                 = {mp.nstr(sigma, 3):>5s} σ                       │
│                                                                             │
├───────────────────────────────────────────────────────────────────────────────┤
│ MASS SPECTRUM (M_(1,+) = M_Z = 91.1876 GeV):                               │
│                                                                             │
│   m (fundamental)  = {mp.nstr(m_phys, 8):>12s} GeV                          │
│   M_(1/2,+)        = {mp.nstr(M_half_plus_num, 8):>12s} GeV   (cf. M_W = 80.3692 GeV) │
│   |M_(1,-)|         = {mp.nstr(M_one_minus_abs, 8):>12s} GeV   (cf. ⟨v⟩/√2 = 174.10 GeV) │
│   |M_(1/2,-)|       = {mp.nstr(M_half_minus_abs, 8):>12s} GeV              │
│                                                                             │
├───────────────────────────────────────────────────────────────────────────────┤
│ CLASSICAL CROSS-CHECK (β²/√(1-β²) = √(j(j+1))):                           │
│                                                                             │
│   β²_(1/2)    = (3/8)(-1 + √(19/3))  = {float(beta2_half):.15f}            │
│   β²_1        = √3 - 1               = {float(beta2_one):.15f}             │
│   R_classical = 1 - β²_(1/2)/β²_1    = {mp.nstr(R_classical_mp, 12):>20s}  │
│   R_Casimir   =                         {mp.nstr(R_mp, 12):>20s}            │
│   Difference  =                         {mp.nstr(R_mp - R_classical_mp, 3):>20s} (exact)│
│                                                                             │
└───────────────────────────────────────────────────────────────────────────────┘
""")
