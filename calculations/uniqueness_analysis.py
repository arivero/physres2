#!/usr/bin/env python3
"""
Analysis of uniqueness for the polynomial eigenvalue equation
M^4 - M^2 * s(s+1) * m^2 + s(s+1) * m^4 = 0

Investigation of the ratio R = 1 - x_+(1/2) / x_+(1)
"""

from sympy import (
    symbols, sqrt, Rational, Integer, simplify, expand, radsimp,
    nsimplify, oo, solve, Symbol, factor, cancel, together, apart,
    collect, trigsimp, powsimp, sqrtdenest, numer, denom
)
from mpmath import mp, mpf, sqrt as mpsqrt, power

# ==========================================================================
# PART 1: General form and parameter count
# ==========================================================================
print("=" * 70)
print("PART 1: General form and parameter count")
print("=" * 70)
print()

print("""
SETUP: We seek polynomials P(M^2, m^2, s) = 0 that are:
  (a) At most degree 2 in M^2 (quartic in M)
  (b) Polynomial in C = s(s+1), the SO(3) Casimir
  (c) Homogeneous in (M, m)

Homogeneity means each monomial (M^2)^i * (m^2)^j has i + j = d (constant).
Since i <= 2, the monomials are:
  i=2: (M^2)^2 * (m^2)^{d-2} = M^4 * m^{2d-4}
  i=1: (M^2)^1 * (m^2)^{d-1} = M^2 * m^{2d-2}
  i=0: (M^2)^0 * (m^2)^d     = m^{2d}

Dividing by m^{2d}, with x = M^2/m^2:

  alpha(C) * x^2 + beta(C) * x + gamma(C) = 0

where alpha, beta, gamma are polynomials in C = s(s+1).

Normalizing (dividing by alpha), the MONIC form is:

  x^2 + f(C) * x + g(C) = 0

with f, g rational functions of C.

PARAMETER COUNT (for polynomials of bounded degree in C):

If alpha, beta, gamma are each polynomials of degree <= n in C:
  - Total coefficients: 3(n+1)
  - Minus 1 for overall scaling (we can always normalize the leading
    coefficient of x^2): 3(n+1) - 1 = 3n + 2 free parameters.

For the simplest case n = 1 (linear in C):
  alpha(C) = a0 + a1*C
  beta(C)  = b0 + b1*C
  gamma(C) = c0 + c1*C
  => 5 free parameters after normalization.

For the monic form with f, g each of degree <= 1 in C (rational degree 0):
  f(C) = f0 + f1*C
  g(C) = g0 + g1*C
  => 4 free parameters (f0, f1, g0, g1).

The SPECIFIC polynomial has:
  f(C) = C     => f0 = 0, f1 = 1
  g(C) = -C    => g0 = 0, g1 = -1
  i.e., 2 parameters are zero, 2 are +/-1.
""")

# ==========================================================================
# PART 2: Compute R for the specific polynomial
# ==========================================================================
print("=" * 70)
print("PART 2: Computing roots and R for the specific polynomial")
print("=" * 70)
print()

x = symbols('x')
C = symbols('C')

# The specific polynomial: x^2 + C*x - C = 0
poly_specific = x**2 + C*x - C

# General solution
roots = solve(poly_specific, x)
print("Roots of x^2 + C*x - C = 0:")
for i, r in enumerate(roots):
    print(f"  x_{i} = {r}")
print()

# s = 1/2: C = 3/4
C_half = Rational(3, 4)
# s = 1:   C = 2
C_one  = Rational(2)

x_half_roots = [r.subs(C, C_half) for r in roots]
x_one_roots  = [r.subs(C, C_one) for r in roots]

print(f"s = 1/2 (C = 3/4): roots = {x_half_roots}")
print(f"s = 1   (C = 2):   roots = {x_one_roots}")
print()

# Identify positive roots
x_half_pos = [r for r in x_half_roots if r.evalf() > 0][0]
x_one_pos  = [r for r in x_one_roots  if r.evalf() > 0][0]

print(f"Positive root for s=1/2: x_+(1/2) = {x_half_pos} = {simplify(x_half_pos)}")
print(f"Positive root for s=1:   x_+(1)   = {x_one_pos} = {simplify(x_one_pos)}")
print()

# Verify expected forms
x_half_expected = (-3 + sqrt(57)) / 8
x_one_expected  = -1 + sqrt(3)

print(f"Expected x_+(1/2) = (-3+sqrt(57))/8 = {x_half_expected}")
print(f"Expected x_+(1)   = -1+sqrt(3)      = {x_one_expected}")
print(f"Match s=1/2: {simplify(x_half_pos - x_half_expected) == 0}")
print(f"Match s=1:   {simplify(x_one_pos - x_one_expected) == 0}")
print()

# Compute R
R_sym = 1 - x_half_pos / x_one_pos
R_simplified = simplify(R_sym)
R_radsimp = radsimp(R_sym)
R_together = together(R_sym)

print(f"R = 1 - x_+(1/2)/x_+(1)")
print(f"  = 1 - [(-3+sqrt(57))/8] / [-1+sqrt(3)]")
print(f"  Simplified: {R_simplified}")
print(f"  Radsimp:    {R_radsimp}")
print(f"  Together:   {R_together}")
print()

# ==========================================================================
# PART 2b: Verify the proposed closed form R = (sqrt(19)-3)(sqrt(19)-sqrt(3))/16
# ==========================================================================
print("=" * 70)
print("PART 2b: Verifying R = (sqrt(19)-3)(sqrt(19)-sqrt(3))/16")
print("=" * 70)
print()

R_proposed = (sqrt(19) - 3) * (sqrt(19) - sqrt(3)) / 16

# Expand
R_prop_expanded = expand(R_proposed)
print(f"(sqrt(19)-3)(sqrt(19)-sqrt(3))/16 expanded:")
print(f"  = {R_prop_expanded}")
print()

# Numerical comparison
R_num = float(R_sym.evalf(50))
R_prop_num = float(R_proposed.evalf(50))
print(f"R (from roots):    {R_sym.evalf(50)}")
print(f"R (proposed form): {R_proposed.evalf(50)}")
print(f"Difference:        {(R_sym - R_proposed).evalf(50)}")
print()

# Symbolic verification
# R = 1 - (-3+sqrt(57)) / [8*(sqrt(3)-1)]
# Rationalize: multiply num and den by (sqrt(3)+1):
# = 1 - (-3+sqrt(57))*(sqrt(3)+1) / [8*2]
# = 1 - (-3+sqrt(57))*(sqrt(3)+1) / 16
#
# (-3+sqrt(57))*(sqrt(3)+1) = -3*sqrt(3) - 3 + sqrt(57)*sqrt(3) + sqrt(57)
#                            = -3*sqrt(3) - 3 + sqrt(171) + sqrt(57)
# Now 171 = 9*19, so sqrt(171) = 3*sqrt(19)
# And 57 = 3*19, so sqrt(57) = sqrt(3)*sqrt(19)
#
# = -3*sqrt(3) - 3 + 3*sqrt(19) + sqrt(3)*sqrt(19)
# = -3 - 3*sqrt(3) + 3*sqrt(19) + sqrt(3)*sqrt(19)
#
# So R = 1 - [-3 - 3*sqrt(3) + 3*sqrt(19) + sqrt(3)*sqrt(19)] / 16
#      = [16 + 3 + 3*sqrt(3) - 3*sqrt(19) - sqrt(3)*sqrt(19)] / 16
#      = [19 + 3*sqrt(3) - 3*sqrt(19) - sqrt(3)*sqrt(19)] / 16
#
# Now check (sqrt(19)-3)(sqrt(19)-sqrt(3)):
# = 19 - sqrt(3)*sqrt(19) - 3*sqrt(19) + 3*sqrt(3)
# = 19 + 3*sqrt(3) - 3*sqrt(19) - sqrt(57)
# where sqrt(57) = sqrt(3)*sqrt(19)
# = 19 + 3*sqrt(3) - 3*sqrt(19) - sqrt(3)*sqrt(19)
#
# MATCH! Both are 19 + 3*sqrt(3) - 3*sqrt(19) - sqrt(3)*sqrt(19).

print("ALGEBRAIC VERIFICATION:")
print()
print("Step 1: R = 1 - (-3+sqrt(57)) / [8*(sqrt(3)-1)]")
print()
print("Step 2: Rationalize denominator, multiply by (sqrt(3)+1)/(sqrt(3)+1):")
print("  R = 1 - (-3+sqrt(57))*(sqrt(3)+1) / [8*2]")
print("    = 1 - (-3+sqrt(57))*(sqrt(3)+1) / 16")
print()
print("Step 3: Expand numerator of fraction:")
print("  (-3+sqrt(57))*(sqrt(3)+1)")
print("  = -3*sqrt(3) - 3 + sqrt(171) + sqrt(57)")
print("  = -3*sqrt(3) - 3 + 3*sqrt(19) + sqrt(3)*sqrt(19)")
print("  [using sqrt(171) = sqrt(9*19) = 3*sqrt(19)]")
print("  [and sqrt(57) = sqrt(3*19) = sqrt(3)*sqrt(19)]")
print()
print("Step 4: R = [16 - (-3-3sqrt(3)+3sqrt(19)+sqrt(3)sqrt(19))] / 16")
print("         = [16 + 3 + 3sqrt(3) - 3sqrt(19) - sqrt(3)sqrt(19)] / 16")
print("         = [19 + 3sqrt(3) - 3sqrt(19) - sqrt(3)sqrt(19)] / 16")
print()
print("Step 5: Expand (sqrt(19)-3)(sqrt(19)-sqrt(3)):")
print("  = 19 - sqrt(3)sqrt(19) - 3sqrt(19) + 3sqrt(3)")
print("  = 19 + 3sqrt(3) - 3sqrt(19) - sqrt(3)sqrt(19)")
print()
print("Steps 4 and 5 match. QED: R = (sqrt(19)-3)(sqrt(19)-sqrt(3))/16")
print()

# Double-check symbolically
expr1 = 19 + 3*sqrt(3) - 3*sqrt(19) - sqrt(3)*sqrt(19)
expr2 = expand((sqrt(19)-3)*(sqrt(19)-sqrt(3)))
print(f"  19+3sqrt3-3sqrt19-sqrt3*sqrt19 = {expr1.evalf(30)}")
print(f"  (sqrt19-3)(sqrt19-sqrt3)       = {expr2.evalf(30)}")
print(f"  Difference: {simplify(expr1 - expr2)}")
print()

# ==========================================================================
# PART 3: Uniqueness analysis - codimension
# ==========================================================================
print("=" * 70)
print("PART 3: Uniqueness / codimension analysis")
print("=" * 70)
print()

print("""
Consider the family of monic quadratics in x:

  x^2 + f(C)*x + g(C) = 0

with f, g depending on C = s(s+1).

The positive root is:
  x_+(C) = [-f(C) + sqrt(f(C)^2 - 4*g(C))] / 2

The ratio is:
  R(f, g) = 1 - x_+(C_1) / x_+(C_2)

where C_1 = 3/4 (s=1/2) and C_2 = 2 (s=1).

This is ONE equation (R = R_target) constraining the function pair (f, g).

CASE A: f and g are both LINEAR in C.
  f(C) = f0 + f1*C,  g(C) = g0 + g1*C
  4 free parameters: (f0, f1, g0, g1)
  1 constraint: R = R_target
  => 3-dimensional family of solutions (codimension 1 in 4-space).

  The specific polynomial (f0=0,f1=1,g0=0,g1=-1) lives in this family.
  R = R_target is NOT isolated; it defines a codimension-1 surface.

CASE B: f(C) = alpha*C,  g(C) = beta*C  (proportional to C, no constant term).
  This is the NATURAL subclass where f(0)=g(0)=0 (massless limit consistency:
  when C=0, x^2=0 giving x=0, i.e., M=0).
  2 free parameters: (alpha, beta)
  1 constraint: R = R_target
  => 1-dimensional family (codimension 1 in 2-space, i.e., a curve).

  The specific polynomial has alpha=1, beta=-1.

CASE C: f(C) = alpha*C, g(C) = beta*C with the ADDITIONAL constraint
  that f and g are related by g = -f (i.e., beta = -alpha).
  1 free parameter: alpha
  But R depends on alpha! Let's check...
""")

# Case C analysis
alpha = symbols('alpha', positive=True)
f_C = alpha * C
g_C = -alpha * C

x_plus = (-f_C + sqrt(f_C**2 + 4*alpha*C)) / 2

# At C = 3/4
x_at_half = x_plus.subs(C, Rational(3,4))
x_at_half_simplified = simplify(x_at_half)

# At C = 2
x_at_one = x_plus.subs(C, Rational(2))
x_at_one_simplified = simplify(x_at_one)

R_alpha = 1 - x_at_half / x_at_one
R_alpha_simplified = simplify(R_alpha)

print(f"For f = alpha*C, g = -alpha*C:")
print(f"  x^2 + alpha*C*x - alpha*C = 0")
print(f"  x_+ = [-alpha*C + sqrt(alpha^2*C^2 + 4*alpha*C)] / 2")
print()
print(f"  x_+(3/4) = {x_at_half_simplified}")
print(f"  x_+(2)   = {x_at_one_simplified}")
print()

# Factor out alpha dependence
# x^2 + alpha*C*x - alpha*C = 0
# Let y = x/sqrt(alpha) if alpha > 0? No, that doesn't simplify.
# Actually: x = [-alpha*C + sqrt(alpha^2*C^2 + 4*alpha*C)] / 2
# Factor alpha from under sqrt: sqrt(alpha) * sqrt(alpha*C^2 + 4*C)
# = [-alpha*C + sqrt(alpha)*sqrt(alpha*C^2 + 4*C)] / 2

# The ratio x_+(3/4) / x_+(2):
# = [-alpha*(3/4) + sqrt(alpha)*sqrt(alpha*9/16 + 3)] / [-alpha*2 + sqrt(alpha)*sqrt(4*alpha + 8)]

# Let u = sqrt(alpha):
# = [-u^2*(3/4) + u*sqrt(u^2*9/16 + 3)] / [-2*u^2 + u*sqrt(4*u^2 + 8)]
# = u*[-(3/4)*u + sqrt(9*u^2/16 + 3)] / {u*[-2u + sqrt(4u^2+8)]}
# = [-(3/4)*u + sqrt(9u^2/16 + 3)] / [-2u + sqrt(4u^2+8)]

# So the ratio DOES depend on alpha (or equivalently, u = sqrt(alpha)).
# Therefore R depends on alpha, and R = R_target gives one equation for one unknown.
# This means: in the (alpha, beta=-alpha) one-parameter family, R_target picks
# out ISOLATED points.

print("Does R depend on alpha in the subfamily g = -f?")
print()

# Evaluate R at a few values of alpha
from sympy import N
for a_val in [Rational(1,2), 1, 2, 3]:
    R_val = R_alpha.subs(alpha, a_val)
    print(f"  alpha = {a_val}: R = {N(R_val, 15)}")
print()

print("YES, R depends on alpha. So within the 1-parameter family g=-f,")
print("R = R_target selects ISOLATED points (generically finitely many).")
print()

# Now the general Case B: f = alpha*C, g = beta*C
a_sym, b_sym = symbols('a b')
f_gen = a_sym * C
g_gen = b_sym * C

x_gen = (-f_gen + sqrt(f_gen**2 - 4*g_gen)) / 2
x_gen_half = x_gen.subs(C, Rational(3,4))
x_gen_one = x_gen.subs(C, Rational(2))

R_gen = 1 - x_gen_half / x_gen_one

print("CASE B detail: f = a*C, g = b*C (2 parameters)")
print(f"  x_+ at C=3/4: {simplify(x_gen_half)}")
print(f"  x_+ at C=2:   {simplify(x_gen_one)}")
print()
print("R = R_target gives 1 equation in 2 unknowns (a, b).")
print("=> Generically a 1-dimensional curve in (a,b)-space.")
print("=> The specific point (a=1, b=-1) is NOT isolated in Case B.")
print()

print("""
SUMMARY OF CODIMENSION ANALYSIS:

  Parameter space                        dim    #constraints   codim of R=R*
  ────────────────────────────────────   ────   ────────────   ─────────────
  (f0,f1,g0,g1) [linear in C]             4         1              1
  (alpha,beta) [f=αC, g=βC]               2         1              1
  alpha alone [f=αC, g=-αC]               1         1              0 (isolated)

  In the natural 2-parameter family f=αC, g=βC:
    - R = R_target defines a CURVE (1-dimensional family).
    - The point (1,-1) is one point on this curve.
    - CODIMENSION = 1 within this family.

  In the most constrained family f=αC, g=-αC:
    - R = R_target gives isolated points.
    - (α=1) is one such point.
    - CODIMENSION = 0 (point is isolated = codimension equals dimension = 1).
""")

# ==========================================================================
# PART 4: How constrained is f=C, g=-C?
# ==========================================================================
print("=" * 70)
print("PART 4: Constraints on f = C, g = -C among degree-2-in-x polynomials")
print("=" * 70)
print()

print("""
Among x^2 + f(C)*x + g(C) = 0 with f, g rational in C:

The choice f = C, g = -C is equivalent to specifying:
  (1) f and g are both linear in C with zero constant term: 2 constraints
      (f0 = 0, g0 = 0)
  (2) The coefficients are f1 = 1, g1 = -1: but one is fixed by normalization.

More precisely, in the linear family f = f0 + f1*C, g = g0 + g1*C:
  - Requiring f0 = 0: 1 constraint (removes constant in f)
  - Requiring g0 = 0: 1 constraint (removes constant in g)
  - Requiring g1 = -f1: 1 constraint (anti-symmetry between f and g)
  - Normalization fixes f1 = 1.

  Total: 3 constraints on 4 parameters = completely determined.

But the VALUE of R is determined by any point in the 4-parameter space.
The question is really: how large is the FIBER {(f0,f1,g0,g1) : R = R*}?

Answer: It's a 3-dimensional hypersurface in 4-space (codimension 1).
The point (0,1,0,-1) is just one point on this hypersurface.
The specific polynomial is determined by its structural properties
(linearity in C, no constant terms, antisymmetry g = -f), NOT solely by R.
""")

# ==========================================================================
# PART 5: High-precision numerical computation
# ==========================================================================
print("=" * 70)
print("PART 5: R to 15+ significant figures and verification")
print("=" * 70)
print()

mp.dps = 50  # 50 decimal places

sqrt57 = mpsqrt(57)
sqrt3  = mpsqrt(3)
sqrt19 = mpsqrt(19)

x_half_mp = (-3 + sqrt57) / 8
x_one_mp  = -1 + sqrt3

R_mp = 1 - x_half_mp / x_one_mp

R_proposed_mp = (sqrt19 - 3) * (sqrt19 - sqrt3) / 16

print(f"x_+(1/2) = (-3+sqrt(57))/8 = {x_half_mp}")
print(f"x_+(1)   = -1+sqrt(3)      = {x_one_mp}")
print()
print(f"R = 1 - x_+(1/2)/x_+(1)")
print(f"  = {R_mp}")
print()
print(f"(sqrt(19)-3)(sqrt(19)-sqrt(3))/16")
print(f"  = {R_proposed_mp}")
print()
print(f"Difference: {R_mp - R_proposed_mp}")
print()
print(f"R to 15 significant figures: {mp.nstr(R_mp, 15)}")
print(f"R to 20 significant figures: {mp.nstr(R_mp, 20)}")
print(f"R to 30 significant figures: {mp.nstr(R_mp, 30)}")
print()

# Additional verification: compute R^2 and see if it satisfies a nice polynomial
R_sq = R_mp ** 2
print(f"R^2 = {R_sq}")
print()

# What polynomial does R satisfy?
# R = (sqrt(19)-3)(sqrt(19)-sqrt(3))/16
# Let's find the minimal polynomial of R over Q.
# 16R = (sqrt(19)-3)(sqrt(19)-sqrt(3))
# 16R = 19 - sqrt(57) - 3*sqrt(19) + 3*sqrt(3)
# 16R - 19 = -sqrt(57) - 3*sqrt(19) + 3*sqrt(3)
#
# This involves sqrt(3) and sqrt(19), so R lies in Q(sqrt(3), sqrt(19)).
# Minimal polynomial over Q has degree 4.

R_sym_exact = (sqrt(19) - 3)*(sqrt(19) - sqrt(3)) / 16
minpoly_R = symbols('t')

# Use sympy's minpoly
from sympy import minpoly as compute_minpoly, QQ
t = symbols('t')
try:
    mp_result = compute_minpoly(R_sym_exact, t)
    print(f"Minimal polynomial of R over Q:")
    print(f"  {mp_result} = 0")
    print(f"  Degree: {mp_result.as_poly(t).degree()}")
except Exception as e:
    print(f"minpoly computation: {e}")
print()

# Verify by substitution
if 'mp_result' in dir():
    try:
        val = mp_result.subs(t, R_sym_exact)
        val_simplified = simplify(val)
        print(f"Substituting R into minimal polynomial: {val_simplified}")
    except:
        pass
print()

# ==========================================================================
# FINAL SUMMARY
# ==========================================================================
print("=" * 70)
print("FINAL SUMMARY")
print("=" * 70)
print(f"""
1. GENERAL FORM: After homogeneity, the equation reduces to
     x^2 + f(C)*x + g(C) = 0,  x = M^2/m^2, C = s(s+1)
   with f, g functions of C. For f, g linear in C: 4 free parameters.

2. R = 1 - x_+(1/2)/x_+(1) for the specific polynomial x^2 + Cx - C = 0:

     R = (sqrt(19) - 3)(sqrt(19) - sqrt(3)) / 16
       = {mp.nstr(R_mp, 15)}

3. UNIQUENESS: R = R_target is NOT an isolated point. It defines:
   - A codimension-1 hypersurface in the 4-parameter space (f0,f1,g0,g1)
   - A codimension-1 curve in the 2-parameter space (alpha, beta)
     for f=alpha*C, g=beta*C
   - An isolated point ONLY in the 1-parameter family f=alpha*C, g=-alpha*C

4. The choice f=C, g=-C is fully determined by THREE structural constraints:
   (i) f(0) = 0, (ii) g(0) = 0, (iii) g = -f, plus normalization.
   The value of R then FOLLOWS; it does not serve as an additional constraint.

5. Numerical value verified to 50 decimal places:
   R = {R_mp}
""")
