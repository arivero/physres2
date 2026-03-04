#!/usr/bin/env python3
"""
Complete analysis of Q=3/2 triplet problem.
Place 6 non-negative reals in a 3x2 grid, forming 8 triplets.
Find all configurations where Q = (s1*sqrt(m1)+s2*sqrt(m2)+s3*sqrt(m3))^2/(m1+m2+m3) = 3/2
for all 8 triplets (with appropriate sign choices).
"""

import numpy as np
from itertools import product

# =============================================================================
# PART 0: BASIC SETUP
# =============================================================================
print("="*80)
print("PART 0: BASIC SETUP")
print("="*80)

kp = 2 + np.sqrt(3)  # k+
km = 2 - np.sqrt(3)  # k-

print(f"k+ = 2 + sqrt(3) = {kp:.10f}")
print(f"k- = 2 - sqrt(3) = {km:.10f}")
print(f"k+ * k- = {kp*km:.10f}")
print(f"k+ + k- = {kp+km:.10f}")

# sqrt(k+) = (1+sqrt(3))/sqrt(2), sqrt(k-) = (sqrt(3)-1)/sqrt(2)
skp = np.sqrt(kp)
skm = np.sqrt(km)
print(f"\nsqrt(k+) = {skp:.10f} = (1+sqrt(3))/sqrt(2) = {(1+np.sqrt(3))/np.sqrt(2):.10f}")
print(f"sqrt(k-) = {skm:.10f} = (sqrt(3)-1)/sqrt(2) = {(np.sqrt(3)-1)/np.sqrt(2):.10f}")
print(f"sqrt(k+) + sqrt(k-) = {skp+skm:.10f}, sqrt(6) = {np.sqrt(6):.10f}")
print(f"sqrt(k+) - sqrt(k-) = {skp-skm:.10f}, sqrt(2) = {np.sqrt(2):.10f}")
print(f"sqrt(k+) * sqrt(k-) = {skp*skm:.10f} = 1")

def Q_val(m1, m2, m3, s1=1, s2=1, s3=1):
    """Compute Q = (s1*sqrt(m1) + s2*sqrt(m2) + s3*sqrt(m3))^2 / (m1+m2+m3)"""
    denom = m1 + m2 + m3
    if denom < 1e-30:
        return float('inf')
    num = (s1*np.sqrt(abs(m1)) + s2*np.sqrt(abs(m2)) + s3*np.sqrt(abs(m3)))**2
    return num / denom

print(f"\nQ(k-, 1, k+) all positive = {Q_val(km, 1, kp):.10f}")
print(f"Target: 3/2 = 1.5")

# Check all sign combos for the basic triplet (k-, 1, k+)
print("\nAll sign combos for (k-, 1, k+):")
for s1, s2, s3 in product([1,-1], repeat=3):
    q = Q_val(km, 1, kp, s1, s2, s3)
    print(f"  signs ({s1:+d},{s2:+d},{s3:+d}): Q = {q:.10f}",
          "  <-- Q=3/2" if abs(q-1.5)<1e-8 else "")

# =============================================================================
# PART 1: Characterize Q = 3/2 condition algebraically
# =============================================================================
print("\n" + "="*80)
print("PART 1: ALGEBRAIC ANALYSIS OF Q = 3/2")
print("="*80)

# Q = (s1*sqrt(a) + s2*sqrt(b) + s3*sqrt(c))^2 / (a+b+c) = 3/2
# => 2*(s1*sqrt(a) + s2*sqrt(b) + s3*sqrt(c))^2 = 3*(a+b+c)
# Expand LHS: 2*(a + b + c + 2*s1*s2*sqrt(ab) + 2*s1*s3*sqrt(ac) + 2*s2*s3*sqrt(bc))
# = 2(a+b+c) + 4(s1*s2*sqrt(ab) + s1*s3*sqrt(ac) + s2*s3*sqrt(bc))
# Set equal to 3(a+b+c):
# a + b + c = 4(s1*s2*sqrt(ab) + s1*s3*sqrt(ac) + s2*s3*sqrt(bc))
#
# With x=sqrt(a), y=sqrt(b), z=sqrt(c):
# x^2 + y^2 + z^2 = 4(s1*s2*xy + s1*s3*xz + s2*s3*yz)
#
# Or: x^2 + y^2 + z^2 - 4(e12*xy + e13*xz + e23*yz) = 0
# where e_ij = s_i * s_j

print("Q = 3/2 condition in terms of x=sqrt(a), y=sqrt(b), z=sqrt(c):")
print("  x^2 + y^2 + z^2 = 4*(e12*x*y + e13*x*z + e23*y*z)")
print("  where e_ij = s_i * s_j")
print()
print("Note: if all signs same (e12=e13=e23=+1):")
print("  x^2 + y^2 + z^2 - 4xy - 4xz - 4yz = 0")
print("  This is a quadratic in x: x^2 - 4(y+z)x + (y^2+z^2-4yz) = 0")
print("  x = 2(y+z) ± sqrt(4(y+z)^2 - (y^2+z^2-4yz))")
print("    = 2(y+z) ± sqrt(3y^2 + 12yz + 3z^2)")
print("    = 2(y+z) ± sqrt(3)*sqrt(y^2 + 4yz + z^2)")
print("    = 2(y+z) ± sqrt(3)*sqrt((y+z)^2 + 2yz)")

# Let's verify: for (k-, 1, k+), x=sqrt(k-)=(sqrt(3)-1)/sqrt(2), y=1, z=sqrt(k+)=(1+sqrt(3))/sqrt(2)
x0 = skm; y0 = 1.0; z0 = skp
lhs = x0**2 + y0**2 + z0**2
rhs = 4*(x0*y0 + x0*z0 + y0*z0)
print(f"\nVerify for (k-,1,k+): LHS = {lhs:.8f}, RHS = {rhs:.8f}")

# =============================================================================
# PART 2: FULLY DEGENERATE CASE (a_i = b_i for all i)
# =============================================================================
print("\n" + "="*80)
print("PART 2: FULLY DEGENERATE CASE (a_i = b_i)")
print("="*80)

# If a_i = b_i, then the two choices in each row give the same value.
# All 8 triplets collapse to a single triplet (a1, a2, a3).
# We need Q(a1,a2,a3) = 3/2 with some sign combination.
#
# In terms of ratios: set a2 = 1 (WLOG by scaling).
# x = sqrt(a1), y = 1, z = sqrt(a3)
# x^2 + 1 + z^2 = 4(e12*x + e13*xz + e23*z)
#
# For all-positive signs:
# x^2 + 1 + z^2 = 4x + 4xz + 4z
# x^2 - 4x(1+z) + (1 + z^2 - 4z) = 0
# x = [4(1+z) ± sqrt(16(1+z)^2 - 4(1+z^2-4z))] / 2
#   = 2(1+z) ± sqrt(4(1+z)^2 - (1+z^2-4z))
#   = 2(1+z) ± sqrt(3z^2 + 12z + 3)
#   = 2(1+z) ± sqrt(3)*sqrt(z^2 + 4z + 1)

print("For the all-positive-sign case with a2=1:")
print("  sqrt(a1) = 2(1+sqrt(a3)) ± sqrt(3)*sqrt(a3+4*sqrt(a3)+1)")
print()

# The basic solution (k-, 1, k+) should satisfy this
z = skp  # sqrt(k+) = sqrt(a3)
sol1 = 2*(1+z) - np.sqrt(3)*np.sqrt(z**2 + 4*z + 1)
sol2 = 2*(1+z) + np.sqrt(3)*np.sqrt(z**2 + 4*z + 1)
print(f"With z=sqrt(k+)={z:.8f}:")
print(f"  Solution 1: sqrt(a1) = {sol1:.8f}, a1 = {sol1**2:.8f}, k- = {km:.8f}")
print(f"  Solution 2: sqrt(a1) = {sol2:.8f}, a1 = {sol2**2:.8f}")

# So any degenerate solution is a permutation of (k-*M0, M0, k+*M0) times overall scale
# But we need to be more careful: different sign combinations give different constraints

# Let's think about this more carefully.
# With the degenerate grid | a1 a1 | all 8 triplets give same (a1,a2,a3)
#                          | a2 a2 | but with SAME values, the 8 sign combos
#                          | a3 a3 | correspond to different signs on sqrt
# So we need ALL 8 sign combinations of Q to equal 3/2 simultaneously.
# That's impossible unless some masses are zero!

# Wait, let me re-read the problem. "allowing sign flips on sqrt(m)"
# The problem says we choose ONE element from each row. If a_i = b_i,
# all 8 triplets give the same (a1, a2, a3). But we need Q=3/2 with
# "allowing sign flips". So for EACH triplet, we can choose signs independently?

# Re-reading: "Find all configurations of 6 non-negative reals such that
# all 8 triplets simultaneously satisfy Q = 3/2, allowing sign flips on sqrt(m)."

# I interpret this as: for each of the 8 triplets (m1,m2,m3), there EXISTS
# a choice of signs (s1,s2,s3) such that Q = 3/2.

print("\n" + "-"*40)
print("INTERPRETATION: For each triplet, there EXISTS signs s.t. Q=3/2")
print("-"*40)

# For degenerate case, all 8 triplets are the same (a1,a2,a3).
# We just need (a1,a2,a3) to satisfy Q=3/2 for SOME sign choice.
# A single sign choice works for all 8 (identical) triplets.
# So ANY (a1,a2,a3) with Q=3/2 for some signs gives a degenerate solution.

# But actually, since there are only 3 values in the degenerate case,
# the "8 triplets" are all identical. So the constraint is simply that
# the single triplet (a1,a2,a3) satisfies Q=3/2 with some signs.

# The fundamental solution family (up to scaling) comes from:
# x^2 + y^2 + z^2 = 4(±xy ± xz ± yz)  [with consistent signs e_ij = s_i*s_j]

# Note: e_ij = s_i*s_j means e12*e13*e23 = (s1*s2)(s1*s3)(s2*s3) = s1^2*s2^2*s3^2 = 1
# So we need e12*e13*e23 = +1. The valid sign patterns are:
# (+,+,+), (+,-,-), (-,+,-), (-,-,+) -- these are the 4 with product +1

print("\nSign patterns for e_ij (product must be +1):")
for e12, e13, e23 in product([1,-1], repeat=3):
    if e12*e13*e23 == 1:
        print(f"  e12={e12:+d}, e13={e13:+d}, e23={e23:+d}")

# For each valid sign pattern, solve x^2+y^2+z^2 = 4(e12*xy+e13*xz+e23*yz)
# with x,y,z >= 0

print("\nSolutions for each sign pattern (normalized with y=1):")

for e12, e13, e23 in product([1,-1], repeat=3):
    if e12*e13*e23 != 1:
        continue
    # x^2 - 4(e12 + e13*z)x + (1 + z^2 - 4*e23*z) = 0
    # Discriminant = 16(e12+e13*z)^2 - 4(1+z^2-4*e23*z)
    # = 16*e12^2 + 32*e12*e13*z + 16*e13^2*z^2 - 4 - 4z^2 + 16*e23*z
    # = 12 + (32*e12*e13+16*e23)*z + 12*z^2   [since e12^2=e13^2=1]
    # = 12(1 + z^2) + (32*e12*e13 + 16*e23)*z

    # For (+,+,+): disc = 12 + 12z^2 + 48z = 12(z^2+4z+1)
    # x = 2(1+z) ± sqrt(3(z^2+4z+1))

    # For (+,-,-): e12=+1,e13=-1,e23=-1
    # disc = 12 + 12z^2 + (-32-16)z = 12(1+z^2) - 48z = 12(z^2-4z+1)
    # x = 2(1-z) ± sqrt(3(z^2-4z+1))

    # For (-,+,-): e12=-1,e13=+1,e23=-1
    # disc = 12+12z^2+(-32-16)z wait let me redo
    # 32*e12*e13 + 16*e23 = 32*(-1)*(1) + 16*(-1) = -32-16 = -48
    # disc = 12(1+z^2) - 48z = 12(z^2-4z+1)
    # x = 2(-1+z) ± sqrt(3(z^2-4z+1))  [since 4(e12+e13*z) = 4(-1+z)]

    # For (-,-,+): e12=-1,e13=-1,e23=+1
    # 32*e12*e13+16*e23 = 32*(1)+16*(1) = 48
    # disc = 12(1+z^2)+48z = 12(z^2+4z+1)
    # x = 2(-1-z) ± sqrt(3(z^2+4z+1))  -- but x>=0, and -1-z < 0...
    # Need 2(-1-z) + sqrt(...) >= 0

    print(f"\n  Pattern e12={e12:+d}, e13={e13:+d}, e23={e23:+d}:")
    # Try z = sqrt(k+)
    z_vals = [skp, skm, 0.5, 1.0, 2.0]
    for z in [skp]:
        A = 1
        B = -4*(e12 + e13*z)
        C = 1 + z**2 - 4*e23*z
        disc = B**2 - 4*A*C
        if disc >= 0:
            x1 = (-B + np.sqrt(disc))/(2*A)
            x2 = (-B - np.sqrt(disc))/(2*A)
            for x in [x1, x2]:
                if x >= -1e-10:
                    x = max(x, 0)
                    a1, a2, a3 = x**2, 1.0, z**2
                    # verify
                    # find what s1,s2,s3 gives these e_ij
                    # e12=s1*s2, e13=s1*s3, e23=s2*s3
                    # s1=+1: s2=e12, s3=e13, check: s2*s3=e12*e13 should = e23
                    # e12*e13 = e23 (guaranteed by product=+1)
                    s1, s2, s3 = 1, e12, e13
                    q = Q_val(a1, a2, a3, s1, s2, s3)
                    print(f"    z=sqrt(k+): x={x:.6f}, (a1,a2,a3)=({a1:.6f},{a2:.6f},{a3:.6f}), Q={q:.8f}")

# =============================================================================
# PART 2b: Complete parametric family for degenerate case
# =============================================================================
print("\n" + "="*80)
print("PART 2b: PARAMETRIC FAMILY FOR ALL-POSITIVE SIGNS")
print("="*80)

# For e12=e13=e23=+1 (all signs same):
# x^2 + y^2 + z^2 = 4(xy + xz + yz)
# Set y=1: x = 2(1+z) ± sqrt(3(z^2+4z+1))
# For this to have real solutions: z^2+4z+1 >= 0, i.e. z >= -2+sqrt(3) or z <= -2-sqrt(3)
# Since z>=0, always real.

# The discriminant of z^2+4z+1=0 gives z = -2±sqrt(3)
# Since z >= 0, the inner square root is always positive.

# Let's parametrize: set z=t for t>=0
# x = 2(1+t) - sqrt(3)*sqrt(t^2+4t+1)  or  x = 2(1+t) + sqrt(3)*sqrt(t^2+4t+1)

# For the minus branch to give x>=0:
# 2(1+t) >= sqrt(3)*sqrt(t^2+4t+1)
# 4(1+t)^2 >= 3(t^2+4t+1)
# 4+8t+4t^2 >= 3t^2+12t+3
# t^2 - 4t + 1 >= 0
# t <= 2-sqrt(3) or t >= 2+sqrt(3)

# For 2-sqrt(3) < t < 2+sqrt(3), only the plus branch gives x>=0.
# At t = 2-sqrt(3) = k-: the minus branch gives x=0.
# At t = 2+sqrt(3) = k+: hmm let me check

print("Minus branch x = 2(1+t) - sqrt(3)*sqrt(t^2+4t+1):")
for t in [0, km, 1, kp, 5]:
    inner = t**2 + 4*t + 1
    x = 2*(1+t) - np.sqrt(3)*np.sqrt(inner)
    print(f"  t = {t:.6f}: x = {x:.6f}, x^2 = {x**2:.6f}")

print("\nPlus branch x = 2(1+t) + sqrt(3)*sqrt(t^2+4t+1):")
for t in [0, km, 1, kp, 5]:
    inner = t**2 + 4*t + 1
    x = 2*(1+t) + np.sqrt(3)*np.sqrt(inner)
    print(f"  t = {t:.6f}: x = {x:.6f}, x^2 = {x**2:.6f}")

# So for the (k-,1,k+) solution: z=sqrt(k+), the minus branch gives x=sqrt(k-)
# x = 2(1+sqrt(k+)) - sqrt(3)*sqrt(k+ + 4*sqrt(k+) + 1)
print(f"\nVerify: 2(1+sqrt(k+)) - sqrt(3)*sqrt(k++4sqrt(k+)+1)")
z = skp
val = 2*(1+z) - np.sqrt(3)*np.sqrt(z**2+4*z+1)
print(f"  = {val:.10f}, sqrt(k-) = {skm:.10f}")

# =============================================================================
# PART 3: NON-DEGENERATE SOLUTIONS - ONE PAIR DIFFERS
# =============================================================================
print("\n" + "="*80)
print("PART 3: NON-DEGENERATE CASE - SYSTEMATIC ANALYSIS")
print("="*80)

# Grid: | a1  b1 |
#        | a2  b2 |
#        | a3  b3 |
#
# 8 triplets from choosing left or right in each row.
# Triplet indexed by (c1,c2,c3) in {0,1}^3 where 0=left(a), 1=right(b).
# Triplet values: (row1[c1], row2[c2], row3[c3])
#
# For each triplet, need to find signs (s1,s2,s3) with
# (s1*sqrt(m1)+s2*sqrt(m2)+s3*sqrt(m3))^2 / (m1+m2+m3) = 3/2

# Let's use the parametrization: x_i = sqrt(a_i), y_i = sqrt(b_i)
# The 8 conditions become:
# For (c1,c2,c3): choose r_i = x_i if c_i=0, y_i if c_i=1
# Need: r_1^2 + r_2^2 + r_3^2 = 4(e12*r_1*r_2 + e13*r_1*r_3 + e23*r_2*r_3)
# for some valid sign pattern (e12,e13,e23) with e12*e13*e23 = 1.

# IMPORTANT INSIGHT: The sign pattern (e12,e13,e23) can differ between triplets!
# Each triplet chooses its own signs.

# Let me think about this differently using the structure.
#
# Consider row i: it provides either x_i^2 or y_i^2 as m_i.
# The constraint for each triplet involves sqrt of values from each row.
#
# Key idea: Let's define for each row: x_i, y_i >= 0 (square roots of a_i, b_i).
# For a triplet choosing (r_1, r_2, r_3) with r_i in {x_i, y_i}:
# Need signs (s_1, s_2, s_3) with s_i in {+1,-1} such that
# sum_i r_i^2 = 4 * sum_{i<j} s_i*s_j * r_i * r_j

# Let's define u_i = s_i * r_i. Then:
# sum_i u_i^2 = 4 * sum_{i<j} u_i * u_j
# sum_i u_i^2 - 4*sum_{i<j} u_i*u_j = 0
#
# Let S = u_1+u_2+u_3. Then S^2 = sum u_i^2 + 2*sum u_i*u_j.
# So sum u_i^2 = S^2 - 2*sum u_i*u_j.
# Substitute: S^2 - 2*sum u_i*u_j - 4*sum u_i*u_j = 0
# S^2 = 6 * sum_{i<j} u_i*u_j
# Also: sum_{i<j} u_i*u_j = (S^2 - sum u_i^2)/2
# So: S^2 = 6*(S^2 - sum u_i^2)/2 = 3*(S^2 - sum u_i^2)
# S^2 = 3*S^2 - 3*sum u_i^2
# 3*sum u_i^2 = 2*S^2
# sum u_i^2 = (2/3)*S^2
# This is exactly Q = S^2/sum r_i^2 = S^2/sum u_i^2 = S^2/((2/3)S^2) = 3/2. ✓

# So the condition is: sum u_i^2 = (2/3)*(u_1+u_2+u_3)^2
# where u_i = ±r_i = ±sqrt(row_i choice).

# Equivalently: u_1^2 + u_2^2 + u_3^2 - (2/3)*(u_1+u_2+u_3)^2 = 0
# This defines a quadratic cone in (u_1,u_2,u_3) space.

# Let T = u_1+u_2+u_3.
# (2/3)T^2 = sum u_i^2
# This means: variance of u_i = (1/3)*sum u_i^2 - (1/3*T)^2... let me think differently.

# Actually: sum u_i^2 = (2/3)(sum u_i)^2
# This is a homogeneous quadratic, so a cone.
# On the unit sphere sum u_i^2 = 1, it gives (sum u_i)^2 = 3/2.
# So u_1+u_2+u_3 = ±sqrt(3/2) on the unit sphere.
# This is the intersection of a plane with the unit sphere = a circle!

# So the solutions form a 1-parameter family (up to overall scale).
# Parametrize: u_1 = cos(θ)*v_1 + sin(θ)*w_1 (for some basis of the circle)

# Let's find the circle explicitly.
# Normal to the plane u_1+u_2+u_3 = sqrt(3/2) on unit sphere:
# n = (1,1,1)/sqrt(3), plane at distance sqrt(3/2)/sqrt(... hmm

# On unit sphere: (u1+u2+u3)^2 = 3/2
# The direction (1,1,1) has length sqrt(3).
# u.n_hat = (u1+u2+u3)/sqrt(3) = ±sqrt(3/2)/sqrt(3) = ±1/sqrt(2)
# So the circle has cos(α) = 1/sqrt(2), α = π/4.
# Radius of circle = sin(π/4) = 1/sqrt(2).

# Center of circle: (1/sqrt(2)) * (1,1,1)/sqrt(3) = (1,1,1)/sqrt(6)
# Two orthogonal vectors perpendicular to (1,1,1):
# v = (1,-1,0)/sqrt(2), w = (1,1,-2)/sqrt(6)

# Parametrize points on circle:
# u = (1,1,1)/sqrt(6) + (1/sqrt(2))*(cos θ * v + sin θ * w)
# = (1,1,1)/sqrt(6) + cos(θ)*(1,-1,0)/2 + sin(θ)*(1,1,-2)/(2*sqrt(3))

print("Circle parametrization on unit sphere:")
print("  u = (1,1,1)/sqrt(6) + cos(θ)*(1,-1,0)/2 + sin(θ)*(1,1,-2)/(2*sqrt(3))")
print()

def circle_point(theta):
    """Point on the Q=3/2 circle, unit sphere, positive sum branch"""
    u = np.array([1,1,1])/np.sqrt(6) + np.cos(theta)*np.array([1,-1,0])/2 + \
        np.sin(theta)*np.array([1,1,-2])/(2*np.sqrt(3))
    return u

# Verify
for theta in [0, np.pi/3, np.pi/2, np.pi]:
    u = circle_point(theta)
    print(f"  θ={theta:.4f}: u=({u[0]:.6f},{u[1]:.6f},{u[2]:.6f}), |u|^2={np.sum(u**2):.8f}, (sum u)^2={np.sum(u)**2:.8f}")

# Now, the actual masses are r_i = |u_i| with the sign encoded separately.
# A point on the circle gives (u_1,u_2,u_3) where u_i can be positive or negative.
# The mass in row i is u_i^2 (the sign is absorbed into the ±sqrt choice).

# For the grid problem: row i contributes either a_i or b_i.
# Setting u_i = s_i*sqrt(m_i), the 8 triplets correspond to 8 points on (possibly
# different) Q=3/2 circles, one for each choice of a_i vs b_i.

# KEY INSIGHT: Since u_i^2 = m_i (the mass from row i), and the circle lives on
# the unit sphere, we need to SCALE. A general solution has u = R * (point on circle)
# for some R > 0. Then m_i = R^2 * u_i^2 (on unit sphere), and sum m_i = R^2.

# For the grid: for a given row i, the two possible masses are a_i and b_i.
# On the circle, the component u_i determines sqrt(a_i) or sqrt(b_i).
# But the SAME point must be on a circle for each triplet choice!

# Let me reconsider. We have 8 triplets. Each must lie on the Q=3/2 cone.
# The cone is: sum u_i^2 = (2/3)(sum u_i)^2.
# This factors as: (u1+u2+u3)^2 = (3/2)(u1^2+u2^2+u3^2) ... no wait.
# sum u_i^2 = (2/3)(sum u_i)^2 <=> 3*sum u_i^2 = 2*(sum u_i)^2
# <=> 3(u1^2+u2^2+u3^2) - 2(u1+u2+u3)^2 = 0
# = u1^2+u2^2+u3^2 - 4(u1*u2+u1*u3+u2*u3) = 0 ... let me recompute
# 3(u1^2+u2^2+u3^2) - 2(u1^2+u2^2+u3^2+2u1u2+2u1u3+2u2u3) = 0
# u1^2+u2^2+u3^2 - 4(u1u2+u1u3+u2u3) = 0
# Yes, this matches the earlier equation with all e_ij = +1.

# But for the NEGATIVE sum branch (sum u = -sqrt(3/2)):
# Same cone, just other nappe. By u -> -u symmetry, this gives same masses.

# Now the question: we have 8 triplets, each lives on this cone.
# For triplet (c1,c2,c3), the point is (±sqrt(a_i or b_i depending on c_i)).
# The signs are chosen to put it on the cone.

# This is getting complex. Let me just do a numerical search.

print("\n" + "="*80)
print("PART 4: NUMERICAL SEARCH FOR SOLUTIONS")
print("="*80)

def check_grid(a, b, verbose=False):
    """
    Check if a 3x2 grid (a[0..2], b[0..2]) satisfies Q=3/2 for all 8 triplets.
    Returns (success, details).
    For each triplet, try all 8 sign combinations.
    """
    results = []
    all_ok = True
    for c1, c2, c3 in product([0,1], repeat=3):
        m = [a[0] if c1==0 else b[0],
             a[1] if c2==0 else b[1],
             a[2] if c3==0 else b[2]]
        found = False
        best_q = None
        best_signs = None
        for s1, s2, s3 in product([1,-1], repeat=3):
            q = Q_val(m[0], m[1], m[2], s1, s2, s3)
            if best_q is None or abs(q - 1.5) < abs(best_q - 1.5):
                best_q = q
                best_signs = (s1, s2, s3)
            if abs(q - 1.5) < 1e-8:
                found = True
                if not verbose:
                    break
        if not found:
            all_ok = False
        results.append({
            'choices': (c1,c2,c3),
            'masses': m,
            'found': found,
            'best_q': best_q,
            'best_signs': best_signs
        })
    if verbose:
        for r in results:
            status = "OK" if r['found'] else "FAIL"
            print(f"  Triplet {r['choices']}: m=({r['masses'][0]:.6f},{r['masses'][1]:.6f},{r['masses'][2]:.6f}) "
                  f"Q={r['best_q']:.8f} signs={r['best_signs']} [{status}]")
    return all_ok, results

# First test: fully degenerate with (k-, 1, k+)
print("\nTest 1: Fully degenerate (k-, 1, k+)")
a = [km, 1.0, kp]
b = [km, 1.0, kp]
ok, _ = check_grid(a, b, verbose=True)
print(f"All OK: {ok}")

# Test 2: permutation (1, k+, k-)
print("\nTest 2: Fully degenerate (1, k+, k-)")
a = [1.0, kp, km]
b = [1.0, kp, km]
ok, _ = check_grid(a, b, verbose=True)
print(f"All OK: {ok}")

# Test 3: (k+, k-, 1)
print("\nTest 3: Fully degenerate (k+, k-, 1)")
a = [kp, km, 1.0]
b = [kp, km, 1.0]
ok, _ = check_grid(a, b, verbose=True)
print(f"All OK: {ok}")
PYEOF