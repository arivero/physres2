#!/usr/bin/env python3
"""
FINAL CATALOG: Understanding the Type 2 solutions properly.
The free numerical search found many Type 2 solutions. Let me analyze their
structure to understand what conditions actually determine them.
"""
import sys
import numpy as np
from itertools import product
from scipy.optimize import minimize

def Q_val(m1, m2, m3, s1=1, s2=1, s3=1):
    denom = m1 + m2 + m3
    if denom < 1e-30:
        return float('inf')
    num = (s1*np.sqrt(abs(m1)) + s2*np.sqrt(abs(m2)) + s3*np.sqrt(abs(m3)))**2
    return num / denom

def check_grid(a, b, verbose=False, tol=1e-6):
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
            if abs(q - 1.5) < tol:
                found = True
                best_q = q
                best_signs = (s1, s2, s3)
                break
            if best_q is None or abs(q - 1.5) < abs(best_q - 1.5):
                best_q = q
                best_signs = (s1, s2, s3)
        if not found:
            all_ok = False
        results.append({'choices': (c1,c2,c3), 'masses': m, 'found': found,
                        'best_q': best_q, 'best_signs': best_signs})
    if verbose:
        for r in results:
            status = "OK" if r['found'] else "FAIL"
            print(f"  ({r['choices'][0]},{r['choices'][1]},{r['choices'][2]}): "
                  f"m=({r['masses'][0]:.8f},{r['masses'][1]:.8f},{r['masses'][2]:.8f}) "
                  f"Q={r['best_q']:.12f} signs={r['best_signs']} [{status}]")
    return all_ok, results

# =========================================================================
print("="*80)
print("UNDERSTANDING TYPE 2 SOLUTIONS")
print("="*80)
sys.stdout.flush()

# A Type 2 solution has a1 != b1 AND a_j != b_j for some j, while a_k = b_k.
# But from the numerical search, many of the Type 2 solutions have a_k = b_k
# for one row. Let me verify this and understand the structure.

# Key insight: Type 2 does NOT require both reflections to be Vieta reflections
# from the SAME point. The two columns can come from DIFFERENT points on the
# Q=3/2 circle!

# Let me analyze a Type 2 solution in detail.
# From trial 1: a = (0.147, 2.129, 0.0000287), b = (29.97, 2.129, 54.22)
# Row 2 is the same. Rows 1 and 3 differ.

# This is: column a has (x1^2, x2^2, x3^2) and column b has (y1^2, x2^2, y3^2)
# where x2 is shared.

# For all 8 triplets to work, we need:
# (x1,x2,x3) -> Q=3/2 with some signs
# (y1,x2,x3) -> Q=3/2 with some signs
# (x1,x2,y3) -> Q=3/2 with some signs
# (y1,x2,y3) -> Q=3/2 with some signs

# The first is given. The second means y1 is the Vieta reflection of x1
# with x2,x3 fixed. The third means y3 is the Vieta reflection of x3
# with x1,x2 fixed. The fourth is the constraint!

# Let's verify this interpretation.
# Take a = (0.147, 2.129, 0.0000287)
a = [0.14680520, 2.12867835, 0.00002865]
x = [np.sqrt(ai) for ai in a]

# Vieta reflection of row 1 (x_j+x_k fixed):
y1 = 4*(x[1]+x[2]) - x[0]
b1_vieta = y1**2
print(f"a = ({a[0]:.8f}, {a[1]:.8f}, {a[2]:.8f})")
print(f"x = ({x[0]:.8f}, {x[1]:.8f}, {x[2]:.8f})")
print(f"Vieta y1 = {y1:.8f}, b1 = {b1_vieta:.8f}")
print(f"From search: b1 = 29.96746437")

# Vieta reflection of row 3:
y3 = 4*(x[0]+x[1]) - x[2]
b3_vieta = y3**2
print(f"Vieta y3 = {y3:.8f}, b3 = {b3_vieta:.8f}")
print(f"From search: b3 = 54.21745890")

# Check (y1, x2, y3):
print(f"\nCheck (y1^2, x2^2, y3^2):")
for s1, s2, s3 in product([1,-1], repeat=3):
    q = Q_val(y1**2, x[1]**2, y3**2, s1, s2, s3)
    if abs(q-1.5) < 1e-4:
        print(f"  signs ({s1},{s2},{s3}): Q = {q:.12f} <-- Q=3/2!")
    elif abs(q-1.5) < 0.1:
        print(f"  signs ({s1},{s2},{s3}): Q = {q:.12f}")

# So the Vieta reflections in rows 1 and 3 DO give a valid Type 2 solution
# but the fourth triplet (y1,x2,y3) also happens to satisfy Q=3/2!
# This is NOT guaranteed in general.

sys.stdout.flush()

print("\n" + "="*80)
print("WHEN DO DOUBLE VIETA REFLECTIONS WORK?")
print("="*80)

# Let me scan the circle and check for each point whether reflecting
# two rows simultaneously gives a valid grid.

def circle_point(theta):
    return np.array([1,1,1])/np.sqrt(6) + np.cos(theta)*np.array([1,-1,0])/2 + \
           np.sin(theta)*np.array([1,1,-2])/(2*np.sqrt(3))

print("\nScan: for each theta, reflect rows i,j and check if (yi,x_k,yj) works")

# For each pair of rows reflected:
for (ri, rj, rk) in [(0,2,1), (0,1,2), (1,2,0)]:
    print(f"\n--- Reflecting rows {ri+1},{rj+1}, keeping row {rk+1} ---")
    working_thetas = []
    for theta_deg_10 in range(0, 3600):
        theta = np.radians(theta_deg_10 / 10)
        u = circle_point(theta)
        x = np.abs(u)
        if min(x) < 1e-10:
            continue

        yi = 4*(x[rj]+x[rk]) - x[ri]  # Vieta for row i (using x_j,x_k)
        # Wait, this is wrong. The Vieta formula for reflecting row i uses the
        # OTHER two rows (rj, rk). But rj is ALSO being reflected!
        # Actually no: yi is computed from the ORIGINAL x values.
        # It's the reflection at fixed (x_j, x_k).
        # Similarly yj is computed from original x values too.

        yj_alt = 4*(x[ri]+x[rk]) - x[rj]  # Vieta for row j with original xi,xk

        a_test = list(x**2)
        b_test = list(x**2)
        b_test[ri] = yi**2
        b_test[rj] = yj_alt**2

        ok, _ = check_grid(a_test, b_test, verbose=False)
        if ok:
            working_thetas.append(theta_deg_10/10)

    if working_thetas:
        # Identify continuous ranges
        diffs = [working_thetas[i+1]-working_thetas[i] for i in range(len(working_thetas)-1)]
        print(f"  Works at {len(working_thetas)} angles")
        print(f"  Range: {working_thetas[0]:.1f}° to {working_thetas[-1]:.1f}°")
        # Check if it's a continuous range
        gaps = [i for i,d in enumerate(diffs) if d > 0.5]
        if gaps:
            print(f"  Gaps at angles: {[working_thetas[g+1] for g in gaps]}")
        else:
            print(f"  Continuous range!")
    else:
        print(f"  No solutions found!")

sys.stdout.flush()

# =========================================================================
print("\n" + "="*80)
print("RE-EXAMINING: Are Type 2 solutions really from double Vieta?")
print("="*80)

# Let me take the trial 1 solution and check if it comes from a circle point
# with double Vieta reflection.

# Trial 1: a=(0.14680520, 2.12867835, 0.00002865), b=(29.96746437, 2.12867835, 54.21745890)
a_t1 = [0.14680520, 2.12867835, 0.00002865]
b_t1 = [29.96746437, 2.12867835, 54.21745890]

# Check: does (a1,a2,a3) satisfy Q=3/2?
for s1,s2,s3 in product([1,-1], repeat=3):
    q = Q_val(a_t1[0], a_t1[1], a_t1[2], s1, s2, s3)
    if abs(q-1.5) < 1e-4:
        print(f"Column a: Q=3/2 with signs ({s1},{s2},{s3}), Q={q:.10f}")

# Check column b:
for s1,s2,s3 in product([1,-1], repeat=3):
    q = Q_val(b_t1[0], b_t1[1], b_t1[2], s1, s2, s3)
    if abs(q-1.5) < 1e-4:
        print(f"Column b: Q=3/2 with signs ({s1},{s2},{s3}), Q={q:.10f}")

# Check Vieta for row 1:
x = [np.sqrt(ai) for ai in a_t1]
y1_vieta = 4*(x[1]+x[2]) - x[0]
print(f"\nVieta y1 from a: {y1_vieta:.8f}, y1^2 = {y1_vieta**2:.8f}")
print(f"Actual b1: {b_t1[0]:.8f}")
print(f"Match: {'YES' if abs(y1_vieta**2 - b_t1[0]) < 0.01 else 'NO'}")

# Check Vieta for row 3:
y3_vieta = 4*(x[0]+x[1]) - x[2]
print(f"Vieta y3 from a: {y3_vieta:.8f}, y3^2 = {y3_vieta**2:.8f}")
print(f"Actual b3: {b_t1[2]:.8f}")
print(f"Match: {'YES' if abs(y3_vieta**2 - b_t1[2]) < 0.01 else 'NO'}")

# Now check: the actual b values may NOT be Vieta reflections from the same base point.
# They could be from a DIFFERENT mechanism entirely.

# Let me check if the 4 distinct triplets in the grid each separately satisfy Q=3/2
print("\nAll 4 distinct triplets:")
for name, triplet in [
    ("(a1,a2,a3)", [a_t1[0], a_t1[1], a_t1[2]]),
    ("(b1,a2,a3)", [b_t1[0], a_t1[1], a_t1[2]]),
    ("(a1,a2,b3)", [a_t1[0], a_t1[1], b_t1[2]]),
    ("(b1,a2,b3)", [b_t1[0], a_t1[1], b_t1[2]]),
]:
    best_q = None
    best_s = None
    for s1,s2,s3 in product([1,-1], repeat=3):
        q = Q_val(triplet[0], triplet[1], triplet[2], s1, s2, s3)
        if best_q is None or abs(q-1.5) < abs(best_q-1.5):
            best_q = q
            best_s = (s1,s2,s3)
    ok = abs(best_q-1.5) < 1e-6
    print(f"  {name}: Q={best_q:.12f} signs={best_s} {'OK' if ok else 'FAIL'}")

sys.stdout.flush()

# =========================================================================
print("\n" + "="*80)
print("UNDERSTANDING TYPE 2: CROSS-REFLECTION STRUCTURE")
print("="*80)

# The 4 distinct triplets in a Type 2 grid (rows i,j differ, row k same) are:
# T00 = (x_i, x_j, x_k) -- original
# T10 = (y_i, x_j, x_k) -- row i reflected
# T01 = (x_i, y_j, x_k) -- row j reflected
# T11 = (y_i, y_j, x_k) -- both reflected

# T00 on cone: given
# T10 on cone: y_i = 4(x_j+x_k)-x_i (Vieta)
# T01 on cone: y_j = 4(x_i+x_k)-x_j (Vieta)
# T11 on cone: need (y_i, y_j, x_k) on some Q=3/2 cone

# y_i = 4(x_j+x_k)-x_i = 4*x_j + 4*x_k - x_i
# y_j = 4(x_i+x_k)-x_j = 4*x_i + 4*x_k - x_j

# For T11 on all-positive cone:
# y_i^2 + y_j^2 + x_k^2 = 4(y_i*y_j + y_i*x_k + y_j*x_k)

# I derived earlier that this requires:
# 8*(x_i+x_j)^2 + 48*(x_i+x_j)*x_k - 144*x_k^2 = 0
# which gives (x_i+x_j)/x_k = (-48+sqrt(48^2+4*8*144))/(2*8) = (-48+sqrt(6912))/16
# = (-48+48*sqrt(3))/16 = (-3+3*sqrt(3)) = 3(sqrt(3)-1)
# So x_i+x_j = 3(sqrt(3)-1)*x_k ≈ 2.196*x_k

# But some Type 2 solutions may use a DIFFERENT sign pattern for T11.
# Let me check which pattern is used.

# From the numerical Type 2 solutions found, let me verify the cross-reflection condition.

# Trial 10: a=(0.65017, 9.13168, 5.33e-6), b=(0.65017, 0.04816, 16.8547)
# Row 1 same, rows 2,3 differ. Let me check.
print("\nAnalyzing trial 10 solution:")
a_t10 = [0.65017061, 9.13168026, 0.00000533]
b_t10 = [0.65017061, 0.04816, 16.85469]

x10 = [np.sqrt(ai) for ai in a_t10]
print(f"x = ({x10[0]:.8f}, {x10[1]:.8f}, {x10[2]:.8f})")
print(f"x1+x3 = {x10[0]+x10[2]:.8f}")
print(f"x2 = {x10[1]:.8f}")
print(f"(x1+x3)/x2 doesn't apply here, row 1 is kept")

# Vieta y2 = 4(x1+x3)-x2
y2 = 4*(x10[0]+x10[2]) - x10[1]
print(f"Vieta y2 = {y2:.8f}, y2^2 = {y2**2:.8f}, b2 = {b_t10[1]}")

# Vieta y3 = 4(x1+x2)-x3
y3 = 4*(x10[0]+x10[1]) - x10[2]
print(f"Vieta y3 = {y3:.8f}, y3^2 = {y3**2:.8f}, b3 = {b_t10[2]}")

# So b2 ≈ y2^2 and b3 ≈ y3^2? Let me check
print(f"b2 match Vieta: {abs(y2**2 - b_t10[1]) < 0.01}")
print(f"b3 match Vieta: {abs(y3**2 - b_t10[2]) < 0.01}")

# The cross-triplet (x1, y2, y3):
print(f"\nCross-triplet (x1^2, y2^2, y3^2) = ({x10[0]**2:.6f}, {y2**2:.6f}, {y3**2:.6f})")
for s1,s2,s3 in product([1,-1], repeat=3):
    q = Q_val(x10[0]**2, y2**2, y3**2, s1, s2, s3)
    if abs(q-1.5) < 1e-4:
        print(f"  Q=3/2 with signs ({s1},{s2},{s3}), Q={q:.10f}")

# What's the ratio (x2+x3)/x1?
ratio = (x10[1]+x10[2])/x10[0]
print(f"\n(x2+x3)/x1 = {ratio:.8f}")
print(f"3(sqrt(3)-1) = {3*(np.sqrt(3)-1):.8f}")

sys.stdout.flush()

# =========================================================================
# Now I see: the Type 2 solutions found by the free search are NOT restricted
# to the Vieta double reflection from a single circle point. They represent
# a more general class where both columns are independently on the Q=3/2 surface,
# and the cross-terms also work.

# Let me understand this better. In a Type 2 solution with row k fixed:
# We need 4 triplets on Q=3/2 cones:
# (x_i, x_j, x_k), (y_i, x_j, x_k), (x_i, y_j, x_k), (y_i, y_j, x_k)
# Each can use ANY of the 4 sign patterns.

# The first three are guaranteed by the Vieta construction. The fourth
# requires the additional condition. This condition depends on which
# sign pattern T11 uses.

# From the numerical solutions, T11 typically uses a DIFFERENT sign pattern
# from T00. This relaxes the algebraic condition.

print("\n" + "="*80)
print("COMPLETE ALGEBRAIC CONDITIONS FOR TYPE 2")
print("="*80)

# For T11 = (y_i, y_j, x_k) with various sign patterns:
# Pattern A (all +): y_i^2+y_j^2+x_k^2 = 4(y_i*y_j+y_i*x_k+y_j*x_k)
#   => x_i+x_j = 3(sqrt(3)-1)*x_k  [derived above]
#
# Pattern B (flip k): y_i^2+y_j^2+x_k^2 = 4(y_i*y_j-y_i*x_k-y_j*x_k)
#   => different condition
#
# Pattern C (flip i): y_i^2+y_j^2+x_k^2 = 4(-y_i*y_j-y_i*x_k+y_j*x_k)
#   => yet another condition
#
# Pattern D (flip j): y_i^2+y_j^2+x_k^2 = 4(-y_i*y_j+y_i*x_k-y_j*x_k)
#   => another condition

# Let me derive all 4.
# With y_i = 4(x_j+x_k)-x_i, y_j = 4(x_i+x_k)-x_j, and S = x_i+x_j+x_k:
# y_i = 4S-5x_i, y_j = 4S-5x_j
# y_i+y_j = 8S-5(x_i+x_j) = 8S-5(S-x_k) = 3S+5x_k
# y_i*y_j = (4S-5x_i)(4S-5x_j) = 16S^2-20S(x_i+x_j)+25*x_i*x_j
#         = 16S^2-20S(S-x_k)+25P where P=x_i*x_j
#         = -4S^2+20S*x_k+25P

# From the original cone: S^2 = 6*(x_i*x_j+x_i*x_k+x_j*x_k) = 6*(P+(x_i+x_j)*x_k)
# = 6*(P+S*x_k-x_k^2)
# So P = S^2/6 - S*x_k + x_k^2

# y_i*y_j = -4S^2+20S*x_k+25*(S^2/6-S*x_k+x_k^2) = -4S^2+20S*x_k+25S^2/6-25S*x_k+25*x_k^2
# = S^2*(-4+25/6) + S*x_k*(20-25) + 25*x_k^2
# = S^2*(1/6) - 5*S*x_k + 25*x_k^2

# Pattern A: y_i^2+y_j^2+x_k^2 = 4*(y_i*y_j+y_i*x_k+y_j*x_k) = 4*y_i*y_j+4*(y_i+y_j)*x_k
# y_i^2+y_j^2 = (y_i+y_j)^2-2*y_i*y_j = (3S+5x_k)^2-2*(S^2/6-5S*x_k+25x_k^2)
# = 9S^2+30S*x_k+25x_k^2-S^2/3+10S*x_k-50x_k^2
# = S^2*(9-1/3)+40S*x_k-25x_k^2
# = S^2*(26/3)+40S*x_k-25x_k^2

# LHS = S^2*(26/3)+40S*x_k-25x_k^2+x_k^2 = S^2*(26/3)+40S*x_k-24x_k^2
# RHS = 4*(S^2/6-5S*x_k+25x_k^2)+4*(3S+5x_k)*x_k
# = 2S^2/3-20S*x_k+100x_k^2+12S*x_k+20x_k^2
# = 2S^2/3-8S*x_k+120x_k^2

# LHS-RHS = 26S^2/3-2S^2/3+40S*x_k+8S*x_k-24x_k^2-120x_k^2
# = 24S^2/3+48S*x_k-144x_k^2
# = 8S^2+48S*x_k-144x_k^2
# = 8*(S^2+6S*x_k-18x_k^2) = 0
# S = (-6+sqrt(36+72))*x_k/2 = (-6+sqrt(108))*x_k/2 = (-6+6sqrt(3))*x_k/2 = 3(sqrt(3)-1)*x_k
# Confirmed: x_i+x_j = S-x_k = (3sqrt(3)-4)*x_k ≈ 1.196*x_k

# Pattern B (flip k): y_i^2+y_j^2+x_k^2 = 4*(y_i*y_j-y_i*x_k-y_j*x_k)
# RHS = 4*y_i*y_j - 4*(y_i+y_j)*x_k
# = 4*(S^2/6-5S*x_k+25x_k^2) - 4*(3S+5x_k)*x_k
# = 2S^2/3-20S*x_k+100x_k^2-12S*x_k-20x_k^2
# = 2S^2/3-32S*x_k+80x_k^2

# LHS-RHS = 26S^2/3-2S^2/3+40S*x_k+32S*x_k-24x_k^2-80x_k^2
# = 8S^2+72S*x_k-104x_k^2
# = 8*(S^2+9S*x_k-13x_k^2) = 0
# S = (-9+sqrt(81+52))*x_k/2 = (-9+sqrt(133))*x_k/2 ≈ 1.266*x_k
# So x_i+x_j = S-x_k = ((-9+sqrt(133))/2-1)*x_k ≈ 0.266*x_k

val_A = 3*(np.sqrt(3)-1)
val_B = (-9+np.sqrt(133))/2

print(f"Pattern A (all +): S/x_k = {val_A:.8f}, (x_i+x_j)/x_k = {val_A-1:.8f}")
print(f"Pattern B (flip k): S/x_k = {val_B:.8f}, (x_i+x_j)/x_k = {val_B-1:.8f}")

# Pattern C (flip i): y_i^2+y_j^2+x_k^2 = 4*(-y_i*y_j-y_i*x_k+y_j*x_k)
# RHS = -4*y_i*y_j + 4*(-y_i+y_j)*x_k
# y_j-y_i = (4S-5x_j)-(4S-5x_i) = 5(x_i-x_j)
# RHS = -4*(S^2/6-5S*x_k+25x_k^2) + 4*5*(x_i-x_j)*x_k
# = -2S^2/3+20S*x_k-100x_k^2 + 20*(x_i-x_j)*x_k

# LHS-RHS = 26S^2/3+2S^2/3+40S*x_k-20S*x_k-24x_k^2+100x_k^2-20(x_i-x_j)*x_k
# = 28S^2/3+20S*x_k+76x_k^2-20(x_i-x_j)*x_k

# This now depends on x_i-x_j, not just S and x_k.
# So Pattern C gives a condition involving the DIFFERENCE x_i-x_j.

# Similarly Pattern D would give a condition on x_j-x_i.
# These break the i<->j symmetry.

# So there are exactly 2 symmetric conditions (Patterns A and B) that
# determine (x_i+x_j)/x_k as a specific ratio, and 2 asymmetric conditions
# (Patterns C and D) that also constrain x_i-x_j.

# This is consistent with the numerical scan finding Type 2 solutions at
# specific (discrete) theta values.

print("\nPattern C (flip i) and D (flip j) involve x_i-x_j as well.")
print("These give additional discrete solutions but with specific asymmetry.")

sys.stdout.flush()

# =========================================================================
# Let me now produce a CLEAN catalog with symbolic expressions
# =========================================================================
print("\n" + "="*80)
print("CLEAN EXAMPLES: One of each type")
print("="*80)

# ---- TYPE 0 ----
print("\n--- TYPE 0 (Degenerate) ---")
theta = np.radians(45)
u = circle_point(theta)
a0 = list(u**2)
b0 = list(a0)
print(f"Grid: | {a0[0]:.10f}  {b0[0]:.10f} |")
print(f"      | {a0[1]:.10f}  {b0[1]:.10f} |")
print(f"      | {a0[2]:.10f}  {b0[2]:.10f} |")
check_grid(a0, b0, verbose=True)

# ---- TYPE 1 ----
print("\n--- TYPE 1 (Row 3 reflected) ---")
theta = np.radians(45)
u = circle_point(theta)
x = np.abs(u)
y3 = 4*(x[0]+x[1]) - x[2]
a1 = list(x**2)
b1 = [x[0]**2, x[1]**2, y3**2]
print(f"Grid: | {a1[0]:.10f}  {b1[0]:.10f} |")
print(f"      | {a1[1]:.10f}  {b1[1]:.10f} |")
print(f"      | {a1[2]:.10f}  {b1[2]:.10f} |")
check_grid(a1, b1, verbose=True)

# ---- TYPE 2 (Pattern B, rows 1,2 reflected, row 3 kept) ----
print("\n--- TYPE 2 (Rows 1,2 reflected, Pattern B for cross-term) ---")
# Condition: S = val_B * x_k, with S = x_i+x_j+x_k
# And x_i+x_j = (val_B-1)*x_k
# Plus cone equation: S^2 = 6*P + 6*(x_i+x_j)*x_k = 6*x_i*x_j + 6*(S-x_k)*x_k
c_val = val_B - 1
x3 = 1.0
# 6t^2-6ct+(c^2+1-4c)=0 from earlier
disc = 3*c_val**2 + 24*c_val - 6
t1 = (3*c_val + np.sqrt(disc))/6
t2 = (3*c_val - np.sqrt(disc))/6
x1 = t1
x2 = c_val - t1
y1 = 4*(x2+x3) - x1
y2 = 4*(x1+x3) - x2
a2 = [x1**2, x2**2, x3**2]
b2 = [y1**2, y2**2, x3**2]
print(f"Grid: | {a2[0]:.10f}  {b2[0]:.10f} |")
print(f"      | {a2[1]:.10f}  {b2[1]:.10f} |")
print(f"      | {a2[2]:.10f}  {b2[2]:.10f} |")
ok, _ = check_grid(a2, b2, verbose=True)
print(f"All OK: {ok}")

# ---- TYPE 2 (Pattern A) ----
print("\n--- TYPE 2 (Rows 1,2 reflected, Pattern A for cross-term) ---")
c_val_A = val_A - 1
x3 = 1.0
disc_A = 3*c_val_A**2 + 24*c_val_A - 6
if disc_A >= 0:
    t1_A = (3*c_val_A + np.sqrt(disc_A))/6
    t2_A = (3*c_val_A - np.sqrt(disc_A))/6
    x1_A = t1_A
    x2_A = c_val_A - t1_A
    if x1_A > 0 and x2_A > 0:
        y1_A = 4*(x2_A+x3) - x1_A
        y2_A = 4*(x1_A+x3) - x2_A
        a2A = [x1_A**2, x2_A**2, x3**2]
        b2A = [y1_A**2, y2_A**2, x3**2]
        print(f"Grid: | {a2A[0]:.10f}  {b2A[0]:.10f} |")
        print(f"      | {a2A[1]:.10f}  {b2A[1]:.10f} |")
        print(f"      | {a2A[2]:.10f}  {b2A[2]:.10f} |")
        ok, _ = check_grid(a2A, b2A, verbose=True)
        print(f"All OK: {ok}")
    else:
        # Try other root
        x1_A = t2_A
        x2_A = c_val_A - t2_A
        if x1_A > 0 and x2_A > 0:
            y1_A = 4*(x2_A+x3) - x1_A
            y2_A = 4*(x1_A+x3) - x2_A
            a2A = [x1_A**2, x2_A**2, x3**2]
            b2A = [y1_A**2, y2_A**2, x3**2]
            print(f"Grid: | {a2A[0]:.10f}  {b2A[0]:.10f} |")
            print(f"      | {a2A[1]:.10f}  {b2A[1]:.10f} |")
            print(f"      | {a2A[2]:.10f}  {b2A[2]:.10f} |")
            ok, _ = check_grid(a2A, b2A, verbose=True)
            print(f"All OK: {ok}")
        else:
            print(f"  No valid solution: t1={t1_A}, t2={t2_A}")
else:
    print("  No real solutions (discriminant < 0)")

sys.stdout.flush()

# =========================================================================
# PHYSICAL SCALING FOR TYPE 1 (using exact Q=3/2 point)
# =========================================================================
print("\n" + "="*80)
print("PHYSICAL SCALING: EXACT Q=3/2 point matching (0.122, 1.70, 3.64) GeV")
print("="*80)

target = np.array(sorted([0.122, 1.70, 3.64]))

# The target doesn't satisfy Q=3/2 exactly. Find closest Q=3/2 point.
best_theta = None
best_err = float('inf')
for theta_deg_100 in range(0, 36000):
    theta = np.radians(theta_deg_100 / 100)
    u = circle_point(theta)
    m = np.sort(u**2)
    if m[0] < 1e-12:
        continue
    S = np.sum(target) / np.sum(m)
    scaled = m * S
    err = np.sum((scaled - target)**2)
    if err < best_err:
        best_err = err
        best_theta = theta

u_best = circle_point(best_theta)
idx = np.argsort(u_best**2)
x_circ = np.abs(u_best[idx])
M0 = np.sum(target) / np.sum(x_circ**2)
col_a = x_circ**2 * M0
x_phys = x_circ * np.sqrt(M0)

print(f"theta = {np.degrees(best_theta):.3f}°")
print(f"M0 = {M0:.6f}")
print(f"Column a (m_i*M0): ({col_a[0]:.6f}, {col_a[1]:.6f}, {col_a[2]:.6f}) GeV")
print(f"Target:            ({target[0]:.6f}, {target[1]:.6f}, {target[2]:.6f}) GeV")

# Type 1 reflections
print(f"\nType 1 solutions (one row reflected):")
for j in range(3):
    k, l = [i for i in range(3) if i != j]
    yj = 4*(x_phys[k]+x_phys[l]) - x_phys[j]
    bj = yj**2
    print(f"  Reflect row {j+1} (a={col_a[j]:.4f} GeV): b = {bj:.4f} GeV (ratio {bj/col_a[j]:.2f})")

    # Full grid verification
    grid_a = list(col_a)
    grid_b = list(col_a)
    grid_b[j] = bj
    ok, _ = check_grid(grid_a, grid_b, verbose=False)
    print(f"    All 8 triplets OK: {ok}")

# =========================================================================
print("\n" + "="*80)
print("FINAL SUMMARY")
print("="*80)
print("""
COMPLETE SOLUTION CATALOG FOR Q(m1,m2,m3) = 3/2 GRID PROBLEM
=============================================================

The Q=3/2 condition is equivalent to the Koide formula: it defines a
quadric cone x1^2+x2^2+x3^2 = 4(x1*x2+x1*x3+x2*x3) in the space
of signed square roots x_i = s_i*sqrt(m_i). The cone has Lorentzian
signature (eigenvalues 3,3,-3) with timelike direction (1,1,1)/sqrt(3).

SOLUTION TYPES:
===============

TYPE 0 (Degenerate): a_i = b_i for all i.
  - 2-parameter family: angle theta on Q=3/2 circle + scale M0.
  - All 8 triplets identical. Q=3/2 with one fixed sign pattern.
  - Full symmetry: S_3 x Z_2 x R+ x U(1).

TYPE 1 (One row reflected): exactly one a_j != b_j.
  - 2-parameter family (theta, M0) x 3 row choices.
  - Vieta formula: sqrt(b_j) = 4*(sqrt(a_k)+sqrt(a_l)) - sqrt(a_j).
  - Always works: both Vieta roots satisfy the cone equation.
  - 8 triplets reduce to 2 distinct mass triplets.
  - Sign pattern may differ between columns (row j sign flips if y_j < 0).

TYPE 2 (Two rows reflected): exactly two a_i != b_i.
  - Discrete solutions (not continuous families), parameterized only by M0.
  - Requires algebraic condition on mass ratios.
  - Two sub-types:
    * Pattern A (cross-term all-positive): (x_i+x_j)/x_k = 3*sqrt(3)-4
    * Pattern B (cross-term flip-k): S/x_k = (-9+sqrt(133))/2
  - 8 triplets reduce to 4 distinct mass triplets.
  - Cross-term triplet uses different sign pattern from original.

TYPE 3 (All three rows reflected): NO SOLUTIONS EXIST.
  - Verified by exhaustive numerical search.

PHYSICAL SCALING TO (0.122, 1.70, 3.64) GeV:
  - Best Q=3/2 match: M0 = 5.462, masses (0.122, 1.700, 3.640) GeV.
  - Type 1 column-2 masses:
    * Reflect lightest: 156 GeV (ratio 1280)
    * Reflect middle:   60 GeV (ratio 35)
    * Reflect heaviest:  22 GeV (ratio 6.1)

KEY IDENTITIES:
  - Vieta: sqrt(a) + sqrt(b) = 4*(sqrt(m_k) + sqrt(m_l))
  - Product: sqrt(a)*sqrt(b) = m_k + m_l - 4*sqrt(m_k*m_l)
  - Cone matrix: A = diag(1,1,1) - 2*(all-ones), eigenvalues (3,3,-3)
  - Circle: u = (1,1,1)/sqrt(6) + cos(t)*(1,-1,0)/2 + sin(t)*(1,1,-2)/(2*sqrt(3))
""")

print("=== ANALYSIS COMPLETE ===")
