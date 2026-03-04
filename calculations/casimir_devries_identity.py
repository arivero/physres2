#!/usr/bin/env python3
"""
Verify whether the Casimir ratio R = 1 - x+(1/2)/x+(1) is ALGEBRAICALLY
identical to the de Vries angle R_dV = (sqrt(19)-3)(sqrt(19)-sqrt(3))/16.
"""

from sympy import sqrt, Rational, simplify, nsimplify, expand, factor, radsimp, cancel
from sympy import Symbol, solve

print("=" * 72)
print("ALGEBRAIC IDENTITY CHECK: Casimir ratio vs de Vries angle")
print("=" * 72)

# Casimir ratio
# s=1/2: x^2 + (3/4)x - 3/4 = 0 → 4x^2 + 3x - 3 = 0 → x = (-3 + sqrt(9+48))/8 = (-3+sqrt(57))/8
# s=1:   x^2 + 2x - 2 = 0 → x = (-2 + sqrt(4+8))/2 = (-2+2sqrt(3))/2 = -1+sqrt(3)

x_half = (-3 + sqrt(57)) / 8
x_one = -1 + sqrt(3)

R_casimir = 1 - x_half / x_one

# de Vries angle
R_devries = (sqrt(19) - 3) * (sqrt(19) - sqrt(3)) / 16

print(f"\nR_casimir = 1 - (-3+sqrt(57))/(8*(-1+sqrt(3)))")
print(f"R_devries = (sqrt(19)-3)(sqrt(19)-sqrt(3))/16")

print(f"\nR_casimir numerical = {float(R_casimir):.18f}")
print(f"R_devries numerical = {float(R_devries):.18f}")

# Check if they are the same
diff = simplify(R_casimir - R_devries)
print(f"\nR_casimir - R_devries = {diff}")
print(f"simplify(diff) = {simplify(diff)}")
print(f"expand(diff) = {expand(diff)}")

if diff == 0:
    print("\n*** CONFIRMED: R_casimir = R_devries EXACTLY (algebraic identity) ***")
else:
    print(f"\n*** NOT identical. Difference = {diff} ***")
    print(f"*** Numerical difference = {float(diff):.2e} ***")

# Now prove the identity step by step
print("\n" + "=" * 72)
print("STEP-BY-STEP PROOF")
print("=" * 72)

print("""
R_casimir = 1 - (-3+sqrt(57)) / (8*(-1+sqrt(3)))

Rationalize the denominator by multiplying by (1+sqrt(3))/(1+sqrt(3)):
  denominator: 8*(-1+sqrt(3))*(1+sqrt(3))/(1+sqrt(3)) = 8*(3-1)/(1+sqrt(3)) = 16/(1+sqrt(3))
  Wait, let's do this carefully.
""")

# Step 1: Compute the ratio
ratio = x_half / x_one
print(f"ratio = x+(1/2)/x+(1) = (-3+sqrt(57))/(8*(sqrt(3)-1))")

# Rationalize denominator: multiply by (sqrt(3)+1)/(sqrt(3)+1)
# denom becomes 8*(3-1) = 16
ratio_rationalized = (-3 + sqrt(57)) * (sqrt(3) + 1) / (8 * (sqrt(3) - 1) * (sqrt(3) + 1))
ratio_rationalized_simplified = simplify(ratio_rationalized)
print(f"\nRationalized: (-3+sqrt(57))*(sqrt(3)+1) / (8*2)")
print(f"            = (-3+sqrt(57))*(sqrt(3)+1) / 16")

num = (-3 + sqrt(57)) * (sqrt(3) + 1)
num_expanded = expand(num)
print(f"\nNumerator expanded: {num_expanded}")
print(f"  = -3*sqrt(3) - 3 + sqrt(57)*sqrt(3) + sqrt(57)")
print(f"  = -3*sqrt(3) - 3 + sqrt(171) + sqrt(57)")

# Note: sqrt(171) = sqrt(9*19) = 3*sqrt(19)
# And sqrt(57) = sqrt(3*19) = sqrt(3)*sqrt(19)
print(f"  = -3*sqrt(3) - 3 + 3*sqrt(19) + sqrt(3)*sqrt(19)")

print(f"\nSo ratio = (-3 - 3*sqrt(3) + 3*sqrt(19) + sqrt(3)*sqrt(19)) / 16")
print(f"         = (3*(sqrt(19)-1) + sqrt(3)*(sqrt(19)-3)) / 16")
# Hmm, let me check the other way

print(f"\nAlternatively: ratio = (sqrt(19)*sqrt(3) + 3*sqrt(19) - 3*sqrt(3) - 3) / 16")
print(f"                     = (sqrt(19)*(sqrt(3)+3) - 3*(sqrt(3)+1)) / 16")

# Now R = 1 - ratio = (16 - numerator) / 16
R_num = 16 - num_expanded
R_num_simplified = simplify(R_num)
print(f"\nR = 1 - ratio = (16 - [{num_expanded}]) / 16")
print(f"  numerator of R = {R_num_simplified}")
print(f"  expanded = {expand(R_num_simplified)}")

# de Vries: (sqrt(19)-3)(sqrt(19)-sqrt(3))/16
dv_num = (sqrt(19) - 3) * (sqrt(19) - sqrt(3))
dv_num_expanded = expand(dv_num)
print(f"\nde Vries numerator: (sqrt(19)-3)(sqrt(19)-sqrt(3))")
print(f"  = {dv_num_expanded}")
print(f"  = 19 - sqrt(57) - 3*sqrt(19) + 3*sqrt(3)")

# Check equality of numerators
diff_num = simplify(R_num_simplified - dv_num_expanded)
print(f"\nDifference of numerators: {diff_num}")

if diff_num == 0:
    print("\n*** NUMERATORS ARE IDENTICAL ***")
    print("*** Therefore R_casimir = R_devries is an algebraic identity ***")

print("\n" + "=" * 72)
print("THE IDENTITY IN CLEAN FORM")
print("=" * 72)

print("""
THEOREM: The following identity holds exactly:

    1 - (-3+√57)/(8(√3-1)) = (√19-3)(√19-√3)/16

Proof:
  LHS = 1 - (-3+√57)/(8(√3-1))
      = (8(√3-1) - (-3+√57)) / (8(√3-1))
      = (8√3-8+3-√57) / (8(√3-1))

  Multiply num and denom by (√3+1):
      = (8√3-5-√57)(√3+1) / (8·2)
      = (8·3+8√3-5√3-5-√(57·3)-√57) / 16
      = (24+8√3-5√3-5-√171-√57) / 16
      = (19+3√3-3√19-√57) / 16

  Note √171 = 3√19 and √57 = √(3·19).

  RHS = (√19-3)(√19-√3)/16
      = (19 - √57 - 3√19 + 3√3) / 16

  LHS numerator = 19 + 3√3 - 3√19 - √57
  RHS numerator = 19 + 3√3 - 3√19 - √57

  These are identical. QED.
""")

# Additional observation: what does sqrt(57) have to do with anything?
print("=" * 72)
print("ALGEBRAIC STRUCTURE")
print("=" * 72)

print("""
Key observation: 57 = 3 × 19

The Casimir equation for s=1/2 introduces √57 = √(3·19).
The de Vries formula introduces √19 and √3 separately.
The identity works because √57 = √3 · √19, connecting the two.

The number 19 arises from the Casimir equation:
  4x^2 + 3x - 3 = 0  (s=1/2, cleared denominators)
  discriminant = 9 + 48 = 57 = 3 × 19

And for s=1:
  x^2 + 2x - 2 = 0
  discriminant = 4 + 8 = 12 = 4 × 3

So 3 appears in BOTH discriminants, while 19 = 57/3 is specific to s=1/2.

The de Vries formula (√19−3)(√19−√3)/16 was originally derived from
entirely different considerations (mass relations). This identity shows
it is equivalent to a ratio of SU(2) Casimir eigenvalue roots.
""")

# Final: compute to very high precision using mpmath
from mpmath import mp, mpf, sqrt as msqrt
mp.dps = 50

R_mp = 1 - (-3 + msqrt(57))/(8*(msqrt(3)-1))
R_dv_mp = (msqrt(19)-3)*(msqrt(19)-msqrt(3))/16
print(f"50-digit verification:")
print(f"R_casimir = {R_mp}")
print(f"R_devries = {R_dv_mp}")
print(f"Difference= {R_mp - R_dv_mp}")
