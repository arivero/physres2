#!/usr/bin/env python3
"""
Part 2: Correct analysis using the cone/circle parametrization.
Q = (s1*sqrt(m1)+s2*sqrt(m2)+s3*sqrt(m3))^2 / (m1+m2+m3) = 3/2
<=> sum u_i^2 = (2/3)(sum u_i)^2  where u_i = s_i*sqrt(m_i)
<=> the vector (u1,u2,u3) lies on the cone C: 3*sum(u_i^2) = 2*(sum u_i)^2
"""

import numpy as np
from itertools import product

def Q_val(m1, m2, m3, s1=1, s2=1, s3=1):
    denom = m1 + m2 + m3
    if denom < 1e-30:
        return float('inf')
    num = (s1*np.sqrt(abs(m1)) + s2*np.sqrt(abs(m2)) + s3*np.sqrt(abs(m3)))**2
    return num / denom

def check_grid(a, b, verbose=False):
    """Check if all 8 triplets satisfy Q=3/2 for some sign choice."""
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
            if abs(q - 1.5) < 1e-6:
                found = True
                best_q = q
                best_signs = (s1, s2, s3)
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
            print(f"  ({r['choices'][0]},{r['choices'][1]},{r['choices'][2]}): "
                  f"m=({r['masses'][0]:.8f},{r['masses'][1]:.8f},{r['masses'][2]:.8f}) "
                  f"Q={r['best_q']:.10f} signs={r['best_signs']} [{status}]")
    return all_ok, results

# =============================================================================
# The cone C: u1^2+u2^2+u3^2 - 4(u1*u2+u1*u3+u2*u3) = 0
# On the unit sphere, this gives (u1+u2+u3)^2 = 3/2
# i.e., u1+u2+u3 = ±sqrt(3/2) = ±sqrt(6)/2
#
# Circle parametrization (positive sum branch):
# u(θ) = (1/sqrt(6))*(1,1,1) + (1/sqrt(2))*[cosθ*(1,-1,0)/sqrt(2) + sinθ*(1,1,-2)/sqrt(6)]
# = (1/sqrt(6), 1/sqrt(6), 1/sqrt(6)) + cosθ*(1/2, -1/2, 0) + sinθ*(1/(2sqrt(3)), 1/(2sqrt(3)), -1/sqrt(3))
# =============================================================================

print("="*80)
print("CONE SOLUTIONS: Q=3/2 CIRCLE ON UNIT SPHERE")
print("="*80)

def circle_point(theta, branch=+1):
    """Point on Q=3/2 circle, unit sphere. branch=+1 or -1 for sum sign."""
    c = np.array([1,1,1]) / np.sqrt(6) * branch
    v1 = np.array([1,-1,0]) / 2
    v2 = np.array([1,1,-2]) / (2*np.sqrt(3))
    return c + np.cos(theta)*v1 + np.sin(theta)*v2

# Find the special angles where a component is zero
# u1 = 1/sqrt(6) + cosθ/2 + sinθ/(2sqrt(3)) = 0
# u2 = 1/sqrt(6) - cosθ/2 + sinθ/(2sqrt(3)) = 0
# u3 = 1/sqrt(6) - sinθ/sqrt(3) = 0

# u3 = 0: sinθ = sqrt(3)/sqrt(6) = 1/sqrt(2), θ = π/4 or 3π/4
# u2 = 0: 1/sqrt(6) + sinθ/(2sqrt(3)) = cosθ/2
#          2/(sqrt(6)*2) + sinθ/sqrt(3) = cosθ ... hmm let me redo
# u2 = 1/sqrt(6) - cosθ/2 + sinθ/(2sqrt(3))
# = 0 => cosθ/2 = 1/sqrt(6) + sinθ/(2sqrt(3))
# cosθ = 2/sqrt(6) + sinθ/sqrt(3)

# Let me just compute numerically
print("\nCircle points at various angles:")
print(f"{'theta':>8s} {'u1':>10s} {'u2':>10s} {'u3':>10s} {'sum':>10s} {'|u|^2':>10s} {'m1':>10s} {'m2':>10s} {'m3':>10s}")
for theta_deg in range(0, 360, 15):
    theta = np.radians(theta_deg)
    u = circle_point(theta)
    m = u**2
    print(f"{theta_deg:8d} {u[0]:10.6f} {u[1]:10.6f} {u[2]:10.6f} {sum(u):10.6f} {sum(u**2):10.6f} {m[0]:10.6f} {m[1]:10.6f} {m[2]:10.6f}")

# Find angles where components are zero
from scipy.optimize import brentq

def u1_val(theta): return circle_point(theta)[0]
def u2_val(theta): return circle_point(theta)[1]
def u3_val(theta): return circle_point(theta)[2]

print("\nZeros of components:")
for name, func in [("u1", u1_val), ("u2", u2_val), ("u3", u3_val)]:
    zeros = []
    for start in np.linspace(0, 2*np.pi, 100):
        try:
            z = brentq(func, start, start + 2*np.pi/100)
            # check not duplicate
            if not any(abs(z - zz) < 0.01 for zz in zeros):
                zeros.append(z)
        except:
            pass
    for z in sorted(zeros):
        u = circle_point(z)
        print(f"  {name}=0 at θ={np.degrees(z):.2f}°: u=({u[0]:.8f},{u[1]:.8f},{u[2]:.8f})")

print("\n" + "="*80)
print("KEY INSIGHT: Masses are u_i^2, so sign of u_i doesn't matter for mass")
print("="*80)

# The mass m_i = u_i^2. Two different points on the cone give the same
# masses if and only if |u_i| are the same (up to overall sign/scale).

# For the GRID problem: row i has two values a_i = x_i^2, b_i = y_i^2.
# For each triplet (c1,c2,c3), we choose r_i ∈ {x_i, y_i} and need to find
# signs s_i such that the point (s_1*r_i_choice) lies on the cone.

# Equivalently: there exist 8 points on the cone (one per triplet), where
# the i-th coordinate of each point has absolute value either x_i or y_i.

# So for row i, the absolute value of coordinate i takes one of two values: x_i or y_i.
# This is a system of 8 cone equations with 6 unknowns (x_i, y_i for i=1,2,3)
# plus 8*3 = 24 binary sign choices.

# Let me think about this using the circle parametrization.
# Scale out: if (u1,u2,u3) is on the cone, so is λ*(u1,u2,u3).
# So we work on the unit sphere circle.

# For 8 points on the cone, each has |u_i| ∈ {x_i, y_i}.
# On the unit sphere: u1^2+u2^2+u3^2 = 1, and each |u_i| is either x_i/R or y_i/R
# where R = sqrt(x_choice_1^2 + x_choice_2^2 + x_choice_3^2).
# But R depends on the choices! So working on the sphere isn't straightforward.

# Let me instead work directly on the cone.
# Cone equation: u1^2+u2^2+u3^2 = 4(u1*u2 + u1*u3 + u2*u3) ... NO!
# Let me re-derive.

# Q = 3/2 <=> (sum u_i)^2 / sum(u_i^2) = 3/2
# <=> 2*(sum u_i)^2 = 3*sum(u_i^2)
# <=> 2*(u1+u2+u3)^2 = 3*(u1^2+u2^2+u3^2)
# <=> 2*u1^2+2*u2^2+2*u3^2+4*u1*u2+4*u1*u3+4*u2*u3 = 3*u1^2+3*u2^2+3*u3^2
# <=> u1^2+u2^2+u3^2 - 4*(u1*u2+u1*u3+u2*u3) = 0     ... (*)

# WAIT. That's wrong. Let me redo:
# 2*(u1+u2+u3)^2 = 3*(u1^2+u2^2+u3^2)
# 2*u1^2+2*u2^2+2*u3^2 + 4*u1u2+4*u1u3+4*u2u3 = 3*u1^2+3*u2^2+3*u3^2
# 4*u1u2+4*u1u3+4*u2u3 = u1^2+u2^2+u3^2
# So: u1^2+u2^2+u3^2 = 4*(u1u2+u1u3+u2u3)   ... (**)

# This is correct! Note that (*) had the wrong sign. Let me verify:
# u1^2+u2^2+u3^2 - 4*(u1u2+u1u3+u2u3) = 0
# This IS the same as (**). OK.

# Now, this is a QUADRATIC form. Let's write it as u^T A u = 0 where
# A = I - 2*J + I = diag stuff... Let me just write it out:
# A = [[1, -2, -2], [-2, 1, -2], [-2, -2, 1]]
# Eigenvalues: trace = 3, the matrix is 3*I - 2*ones = I + 2*(I - ones)
# Hmm: A = I - 2*(ones_matrix). ones_matrix has eigenvalues 3 (eigvec (1,1,1))
# and 0,0 (eigvecs perp to (1,1,1)).
# So A has eigenvalues 1-2*3=-5 (for (1,1,1)) and 1-0=1,1 (for perps).
# Wait: ones_matrix = J/n? No. Let J = matrix of all 1s.
# A_ij = δ_ij - 2*(1 - δ_ij) = 3*δ_ij - 2
# Actually: A_ij = -2 for i≠j, 1 for i=j
# A = I + (-2)(J - I) = 3I - 2J where J = all-ones matrix
# Eigenvalues of J: 3 (for (1,1,1)/√3) and 0,0 (perps)
# Eigenvalues of A = 3I-2J: 3-6=-3 and 3,3

# So A has signature (+,+,-). The cone u^T A u = 0 is a real cone with one sheet.

# Hmm wait: the form is u1^2+u2^2+u3^2 - 4(u1u2+u1u3+u2u3)
# Let me factor it differently:
# = -(u1+u2+u3)^2 + (u1-u2)^2 + (u1-u3)^2 + (u2-u3)^2 ... no
# Actually:
# (u1-u2)^2 + (u1-u3)^2 + (u2-u3)^2 = 2*(u1^2+u2^2+u3^2) - 2*(u1u2+u1u3+u2u3)
# And our condition is u1^2+u2^2+u3^2 = 4*(u1u2+u1u3+u2u3)
# So sum of squares of differences = 2*4*(u1u2+u1u3+u2u3) - 2*(u1u2+u1u3+u2u3)
#   = 6*(u1u2+u1u3+u2u3)

# Let S = u1+u2+u3. Then:
# S^2 = u1^2+u2^2+u3^2 + 2*(u1u2+u1u3+u2u3) = 4P + 2P = 6P
# where P = u1u2+u1u3+u2u3
# So P = S^2/6.
# And u1^2+u2^2+u3^2 = 4P = 2S^2/3.
# On unit sphere: 1 = 2S^2/3, S^2 = 3/2. ✓

# Now for the grid problem:
# We need 8 points on cone (**), where point (c1,c2,c3) has
# |u_i| = x_i if c_i=0, |u_i| = y_i if c_i=1
# (with independent sign choices for each point).

# Let's denote the 8 points as u^(c1c2c3) = (σ1*r1(c1), σ2*r2(c2), σ3*r3(c3))
# where r_i(0)=x_i, r_i(1)=y_i, and σ_i ∈ {±1} (may depend on c).

# Cone condition for each point:
# r1(c1)^2 + r2(c2)^2 + r3(c3)^2 = 4*(σ1σ2*r1(c1)*r2(c2) + σ1σ3*r1(c1)*r3(c3) + σ2σ3*r2(c2)*r3(c3))

# Note: the sign products ε_ij = σ_i*σ_j determine everything about signs.
# And ε12*ε13*ε23 = 1 always.

# So for each triplet we need to find sign products ε_ij (with product 1) such that
# the cone equation holds. The 4 valid sign patterns correspond to choosing
# an overall "negative" coordinate among {none, 1, 2, 3}:
# - none flipped: ε_ij = +1 for all (i,j)
# - flip coord k: ε_ij = +1 if neither or both are k, -1 if exactly one is k

# Let me denote the sign pattern by which coordinate is "flipped": ∅, 1, 2, or 3.
# Pattern ∅: u1^2+u2^2+u3^2 = 4(u1u2+u1u3+u2u3)
# Pattern 1: u1^2+u2^2+u3^2 = 4(-u1u2-u1u3+u2u3)  (flip sign of terms with u1)
# Pattern 2: u1^2+u2^2+u3^2 = 4(-u1u2+u1u3-u2u3)
# Pattern 3: u1^2+u2^2+u3^2 = 4(u1u2-u1u3-u2u3)

# These can be written uniformly. For pattern k (k=1,2,3):
# Replace u_k -> -u_k in the all-positive equation.
# i.e., (-u_k)^2+... = 4*(-u_k*u_j+...) which is same as cone condition on (-u_k, u_j, u_l).

# So the cone equation u^T A u = 0 is invariant under u_k -> -u_k (it must be,
# since A has the symmetry). Wait, no:
# u^T A u = u1^2+u2^2+u3^2 - 4(u1u2+u1u3+u2u3)
# Under u1 -> -u1: u1^2+u2^2+u3^2 - 4(-u1u2-u1u3+u2u3) ≠ u^T A u in general.
# The cone is NOT invariant under coordinate reflections.

# So different sign patterns give different cones.
# Let me define: for sign pattern (ε12,ε13,ε23) with product 1:
# f(u1,u2,u3) = u1^2+u2^2+u3^2 - 4*(ε12*u1u2+ε13*u1u3+ε23*u2u3) = 0

# Now, for the grid problem, each of the 8 triplets must lie on one of these 4 cones
# (with the appropriate cone chosen based on the sign pattern).

# KEY: Instead of working with signs, let me use the u-variables directly.
# For each row i, there are two possible |u_i| values: x_i and y_i.
# The sign of u_i in each triplet is free.
# The constraint is: for each triplet, the vector (u1,u2,u3) lies on
# u^2+v^2+w^2 = 4(uv+uw+vw)  ... i.e. the SINGLE cone with ALL + signs.

# Wait, but the sign of u_i IS the sign choice s_i. So:
# u_i = s_i * |u_i| where |u_i| ∈ {x_i, y_i}
# The cone equation (with all the sign stuff absorbed):
# (s1*r1)^2+(s2*r2)^2+(s3*r3)^2 = 4*(s1*r1*s2*r2 + s1*r1*s3*r3 + s2*r2*s3*r3)
# = r1^2+r2^2+r3^2 = 4*(s1s2*r1r2 + s1s3*r1r3 + s2s3*r2r3)
# But we can choose s_i freely. So yes, each triplet chooses its own cone.

# Equivalently: define U_i = s_i * r_i (signed square root). Then:
# U1^2+U2^2+U3^2 = 4*(U1*U2+U1*U3+U2*U3)
# This is the SINGLE cone equation with real U_i (not necessarily positive).
# And |U_i| ∈ {x_i, y_i}.

# So the problem is: find x_i, y_i >= 0 such that for each of the 8 choices of
# |U_i| ∈ {x_i, y_i}, there exist signs for U_i making the cone equation hold.

# The cone: U1^2+U2^2+U3^2 = 4(U1*U2+U1*U3+U2*U3)
# = (U1+U2+U3)^2 - 2(U1^2+U2^2+U3^2) + 4(U1*U2+U1*U3+U2*U3) = 0 ... no
# Already established: this is equivalent to sum U_i^2 = (2/3)(sum U_i)^2
# i.e. 3*sum(U_i^2) = 2*(sum U_i)^2

print("\n" + "="*80)
print("REFORMULATION")
print("="*80)
print("""
The problem: find x_i, y_i >= 0 (i=1,2,3) such that for each of the 8
combinations (r_1,r_2,r_3) with r_i in {x_i, y_i}, there exist signs
ε_i in {+1,-1} such that:

  (ε_1*r_1 + ε_2*r_2 + ε_3*r_3)^2 = (3/2)*(r_1^2 + r_2^2 + r_3^2)

i.e., |ε_1*r_1 + ε_2*r_2 + ε_3*r_3| = sqrt(3/2) * sqrt(r_1^2+r_2^2+r_3^2)
""")

# This is a cleaner formulation. Let's explore it.

# For a SINGLE triplet (r_1,r_2,r_3) >= 0, the condition is:
# There exist ε_i in {±1} s.t. (ε_1*r_1+ε_2*r_2+ε_3*r_3)^2 = (3/2)*(r_1^2+r_2^2+r_3^2)
#
# The 8 possible values of |ε_1*r_1+ε_2*r_2+ε_3*r_3| come in pairs (by overall sign flip).
# So there are 4 distinct absolute values:
# |r_1+r_2+r_3|, |r_1+r_2-r_3|, |r_1-r_2+r_3|, |r_1-r_2-r_3|
# (with the convention r_i >= 0)
#
# We need at least one of these to equal sqrt(3/2)*||(r_1,r_2,r_3)||.

# Let p = r_1+r_2+r_3, and let R = sqrt(r_1^2+r_2^2+r_3^2).
# Need one of: |p|, |p-2r_3|, |p-2r_2|, |p-2r_1| = sqrt(3/2)*R
# (since r_1+r_2-r_3 = p-2r_3, etc.)

# For the grid: ALL 8 combinations of (r_1,r_2,r_3) must satisfy this.

# =============================================================================
# APPROACH: Think of it as 8 points on a sphere of varying radius
# =============================================================================

# Let's normalize: define t_i = r_i / R where R = ||(r_1,r_2,r_3)||.
# Then t_1^2+t_2^2+t_3^2 = 1 and we need:
# |ε_1*t_1+ε_2*t_2+ε_3*t_3| = sqrt(3/2) for some ε_i.

# The 4 candidates are:
# |t_1+t_2+t_3|, |t_1+t_2-t_3|, |t_1-t_2+t_3|, |t_1-t_2-t_3|
# All are <= sqrt(3) (by Cauchy-Schwarz). Need one to be sqrt(3/2).

# On the unit sphere with t_i >= 0 (first octant):
# t_1+t_2+t_3 is the largest of the 4, ranging from 1 to sqrt(3).
# It equals sqrt(3/2) on a "latitude circle" t_1+t_2+t_3 = sqrt(3/2).

# If t_1+t_2+t_3 < sqrt(3/2), we need one of the others to have absolute value sqrt(3/2).
# |t_1+t_2-t_3| = sqrt(3/2) means t_1+t_2-t_3 = ±sqrt(3/2).
# Since t_i >= 0 and t_1+t_2+t_3 <= sqrt(3), we have t_1+t_2-t_3 <= sqrt(3).
# t_1+t_2-t_3 = sqrt(3/2) means t_3 = t_1+t_2-sqrt(3/2).

# This is getting complex. Let me just do systematic numerics.

print("\n" + "="*80)
print("NUMERICAL SEARCH: Degenerate solutions (a_i = b_i)")
print("="*80)

# For degenerate: all 8 triplets are the same (a_1, a_2, a_3).
# Need: for the single triplet (r_1,r_2,r_3), there exist signs ε_i s.t.
# (ε_1 r_1+ε_2 r_2+ε_3 r_3)^2 = (3/2)(r_1^2+r_2^2+r_3^2)

# On unit sphere (r_1^2+r_2^2+r_3^2=1), need one of:
# |r_1+r_2+r_3|, |r_1+r_2-r_3|, |r_1-r_2+r_3|, |-r_1+r_2+r_3| = sqrt(3/2)

# For r_i >= 0, r_1+r_2+r_3 >= |any other combo|.
# If r_1+r_2+r_3 = sqrt(3/2), that's one family.
# If r_1+r_2-r_3 = sqrt(3/2), then r_3 = r_1+r_2-sqrt(3/2), small r_3.

# Let's parametrize on the unit sphere in first octant.
# r = (sinφ*cosλ, sinφ*sinλ, cosφ)

# Condition r_1+r_2+r_3 = sqrt(3/2):
print("Family A: r_1+r_2+r_3 = sqrt(3/2) with r_1^2+r_2^2+r_3^2 = 1")
print("This is the intersection of a plane and sphere = a circle.")
print()

# Let's find the (k-,1,k+) point on this circle
kp = 2 + np.sqrt(3)
km = 2 - np.sqrt(3)
r_test = np.array([np.sqrt(km), 1.0, np.sqrt(kp)])
R_test = np.linalg.norm(r_test)
t_test = r_test / R_test
print(f"(sqrt(k-), 1, sqrt(k+)) normalized: ({t_test[0]:.6f}, {t_test[1]:.6f}, {t_test[2]:.6f})")
print(f"  sum = {sum(t_test):.6f}, sqrt(3/2) = {np.sqrt(1.5):.6f}")
print(f"  sum of squares = {sum(t_test**2):.6f}")

# The sum exceeds sqrt(3/2)! Check if any other combo works:
sums = [
    t_test[0]+t_test[1]+t_test[2],
    t_test[0]+t_test[1]-t_test[2],
    t_test[0]-t_test[1]+t_test[2],
    -t_test[0]+t_test[1]+t_test[2],
]
print("  All 4 combo sums:", [f"{s:.6f}" for s in sums])
print("  Need one to be ±sqrt(3/2) = ±", f"{np.sqrt(1.5):.6f}")
# None of them match! So (k-,1,k+) does NOT satisfy Q=3/2 with any signs.
# As we saw earlier, the best Q is about 1.17.

# Let me find what DOES satisfy it.
# On the unit sphere first octant, r_1+r_2+r_3 = sqrt(3/2).
# Parametrize: center = (1,1,1)/sqrt(3) * sqrt(1/2) = (1,1,1)*sqrt(1/6)...
# Actually the center of the circle on the sphere at "latitude" sum=sqrt(3/2):
# The unit vector in (1,1,1) direction is n = (1,1,1)/sqrt(3).
# r·n = sqrt(3/2)/sqrt(3) = sqrt(1/2) = 1/sqrt(2)
# So the "height" along n is 1/sqrt(2), and the circle radius is sqrt(1 - 1/2) = 1/sqrt(2).

# Center of circle: (1/sqrt(2)) * (1,1,1)/sqrt(3) = (1,1,1)/sqrt(6)
# Two perpendicular directions on the sphere at this height:
# v1 = (1,-1,0)/sqrt(2), v2 = (1,1,-2)/sqrt(6)

# r(θ) = (1,1,1)/sqrt(6) + (1/sqrt(2))*[cosθ * v1 + sinθ * v2]
#       = (1,1,1)/sqrt(6) + cosθ*(1,-1,0)/2 + sinθ*(1,1,-2)/(2*sqrt(3))

def family_A_point(theta):
    return np.array([1,1,1])/np.sqrt(6) + np.cos(theta)*np.array([1,-1,0])/2 + \
           np.sin(theta)*np.array([1,1,-2])/(2*np.sqrt(3))

print("\n\nFamily A circle (unit sphere, positive sum = sqrt(3/2)):")
print(f"{'θ(°)':>6s} {'r1':>10s} {'r2':>10s} {'r3':>10s} {'r1+r2+r3':>10s} {'Σri²':>10s} {'m1':>10s} {'m2':>10s} {'m3':>10s} {'ratios':>30s}")

for theta_deg in range(0, 360, 30):
    theta = np.radians(theta_deg)
    r = family_A_point(theta)
    m = r**2  # masses (proportional to)
    s = sum(r)
    ss = sum(r**2)
    # Ratios of masses (sort first)
    ms = sorted(m)
    if ms[0] > 1e-10:
        ratios = f"{ms[1]/ms[0]:.4f}, {ms[2]/ms[0]:.4f}"
    else:
        ratios = "has zero"
    print(f"{theta_deg:6d} {r[0]:10.6f} {r[1]:10.6f} {r[2]:10.6f} {s:10.6f} {ss:10.6f} {m[0]:10.6f} {m[1]:10.6f} {m[2]:10.6f} {ratios:>30s}")

# The constraint for the degenerate case is simple: ANY point on this circle works.
# In terms of masses: m_i = r_i^2, and we need r_1+r_2+r_3 = sqrt(3/2) on unit sphere.
# (Or one of the other 3 sign combos equals sqrt(3/2), but on the first octant
# r_i >= 0, the all-positive sum is largest.)

# Actually, we need r_i >= 0 for real masses. Some points on the circle have r_i < 0.
# Those points correspond to choosing the NEGATIVE sign for sqrt(m_i).
# The mass is still m_i = r_i^2 >= 0, but the "natural" sign is negative.

# For the degenerate grid, ALL 8 triplets are the same masses. So we need the
# SAME set of masses to satisfy Q=3/2. Since the masses determine the 4 possible
# |signed sums|, and Q=3/2 needs one of them to equal sqrt(3/2)*R, the degenerate
# solutions are exactly the orbit of the circle under the 6 permutations of coordinates.

# But wait - for the degenerate grid, all triplets are the same. So we just need
# ONE configuration (m1,m2,m3) satisfying Q=3/2. This is a 2-parameter family
# (up to overall scale): the circle parametrized by θ.

# Now for the non-degenerate case, things get interesting.

print("\n" + "="*80)
print("NON-DEGENERATE SOLUTIONS: Systematic approach")
print("="*80)

# We need: for each of 8 combos (r1,r2,r3) with ri in {xi, yi}:
# ∃ signs εi s.t. (ε1*r1+ε2*r2+ε3*r3)^2 = (3/2)*(r1^2+r2^2+r3^2)
#
# With U_i = ε_i * r_i, this is 3*sum(U_i^2) = 2*(sum U_i)^2
# i.e., sum(U_i^2) = 4*sum_{i<j}(U_i*U_j)
#
# CRUCIAL: the sign choices ε_i can be DIFFERENT for different triplets.

# Strategy: use the constraint that PAIRS of triplets share coordinates.
# Triplets (0,0,0) and (1,0,0) differ only in row 1: x_1 vs y_1.
# Both must lie on the cone with some signs.

# Let me write the constraints explicitly for a few pairs.
# Triplet (0,0,0): U = (ε1*x1, ε2*x2, ε3*x3) on cone
# Triplet (1,0,0): V = (δ1*y1, δ2*x2, δ3*x3) on cone
# Both have same rows 2,3 values. But signs may differ!

# This makes it very hard to constrain. Let me try a direct numerical optimization.

from scipy.optimize import minimize, differential_evolution

def objective(params):
    """
    params = [x1, y1, x2, y2, x3, y3] (square roots of masses, all >= 0)
    For each of 8 triplets, compute min over sign choices of |Q - 3/2|.
    Return sum of squared deviations.
    """
    x = [params[0], params[2], params[4]]
    y = [params[1], params[3], params[5]]

    total_err = 0
    for c1, c2, c3 in product([0,1], repeat=3):
        r = [x[0] if c1==0 else y[0],
             x[1] if c2==0 else y[1],
             x[2] if c2==0 else y[2]]  # BUG: should be c3 for row 3
        # Fix:
        r = [x[0] if c1==0 else y[0],
             x[1] if c2==0 else y[1],
             x[2] if c3==0 else y[2]]

        R2 = r[0]**2 + r[1]**2 + r[2]**2
        if R2 < 1e-30:
            continue

        target = np.sqrt(1.5) * np.sqrt(R2)

        # 4 distinct |signed sums|
        sums = [
            abs(r[0]+r[1]+r[2]),
            abs(r[0]+r[1]-r[2]),
            abs(r[0]-r[1]+r[2]),
            abs(-r[0]+r[1]+r[2])
        ]

        min_dev = min((s - target)**2 for s in sums)
        total_err += min_dev

    return total_err

# Also try with the signed square roots (allowing negative = sign choice)
def objective2(params):
    """
    params = [u1a, u1b, u2a, u2b, u3a, u3b] (signed square roots)
    Masses: m_ia = u_ia^2, m_ib = u_ib^2
    """
    ua = [params[0], params[2], params[4]]
    ub = [params[1], params[3], params[5]]

    total_err = 0
    for c1, c2, c3 in product([0,1], repeat=3):
        # For each triplet, the absolute values are fixed, signs are free
        r = [abs(ua[0]) if c1==0 else abs(ub[0]),
             abs(ua[1]) if c2==0 else abs(ub[1]),
             abs(ua[2]) if c3==0 else abs(ub[2])]

        R2 = r[0]**2 + r[1]**2 + r[2]**2
        if R2 < 1e-30:
            continue

        target = np.sqrt(1.5) * np.sqrt(R2)

        sums = [
            abs(r[0]+r[1]+r[2]),
            abs(r[0]+r[1]-r[2]),
            abs(r[0]-r[1]+r[2]),
            abs(-r[0]+r[1]+r[2])
        ]

        min_dev = min((s - target)**2 for s in sums)
        total_err += min_dev

    return total_err

# Search for solutions
print("\nNumerical search using differential evolution...")
np.random.seed(42)

bounds = [(0, 5)] * 6
result = differential_evolution(objective, bounds, seed=42, maxiter=1000, tol=1e-15,
                                 popsize=50, mutation=(0.5,1.5))
print(f"Best objective: {result.fun:.2e}")
if result.fun < 1e-10:
    params = result.x
    x = [params[0], params[2], params[4]]
    y = [params[1], params[3], params[5]]
    print(f"Solution: x = ({x[0]:.8f}, {x[1]:.8f}, {x[2]:.8f})")
    print(f"          y = ({y[0]:.8f}, {y[1]:.8f}, {y[2]:.8f})")
    a = [xi**2 for xi in x]
    b = [yi**2 for yi in y]
    print(f"Masses a: ({a[0]:.8f}, {a[1]:.8f}, {a[2]:.8f})")
    print(f"Masses b: ({b[0]:.8f}, {b[1]:.8f}, {b[2]:.8f})")
    check_grid(a, b, verbose=True)
else:
    params = result.x
    x = [params[0], params[2], params[4]]
    y = [params[1], params[3], params[5]]
    print(f"Best found: x = ({x[0]:.6f}, {x[1]:.6f}, {x[2]:.6f})")
    print(f"            y = ({y[0]:.6f}, {y[1]:.6f}, {y[2]:.6f})")

# Try multiple random seeds
print("\nMultiple random restarts:")
best_overall = float('inf')
best_params = None
for seed in range(20):
    result = differential_evolution(objective, bounds, seed=seed, maxiter=500, tol=1e-15,
                                     popsize=30)
    if result.fun < best_overall:
        best_overall = result.fun
        best_params = result.x
    if result.fun < 1e-10:
        params = result.x
        x = [params[0], params[2], params[4]]
        y = [params[1], params[3], params[5]]
        a = [xi**2 for xi in x]
        b = [yi**2 for yi in y]
        print(f"  Seed {seed}: obj={result.fun:.2e}, a=({a[0]:.6f},{a[1]:.6f},{a[2]:.6f}), b=({b[0]:.6f},{b[1]:.6f},{b[2]:.6f})")

print(f"\nBest overall objective: {best_overall:.2e}")
if best_params is not None:
    params = best_params
    x = [params[0], params[2], params[4]]
    y = [params[1], params[3], params[5]]
    a = [xi**2 for xi in x]
    b = [yi**2 for yi in y]
    print(f"Best params: a=({a[0]:.8f},{a[1]:.8f},{a[2]:.8f})")
    print(f"             b=({b[0]:.8f},{b[1]:.8f},{b[2]:.8f})")
    ok, _ = check_grid(a, b, verbose=True)

print("\n\nDone with Part 2.")
