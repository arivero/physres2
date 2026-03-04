#!/usr/bin/env python3
"""
FINAL COMPREHENSIVE ANALYSIS: Q=3/2 grid problem.
All solution types, verified, with exact arithmetic where possible.
"""
import sys
import numpy as np
from itertools import product

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
                  f"m=({r['masses'][0]:.10f},{r['masses'][1]:.10f},{r['masses'][2]:.10f}) "
                  f"Q={r['best_q']:.12f} signs={r['best_signs']} [{status}]")
    return all_ok, results

def circle_point(theta):
    return np.array([1,1,1])/np.sqrt(6) + np.cos(theta)*np.array([1,-1,0])/2 + \
           np.sin(theta)*np.array([1,1,-2])/(2*np.sqrt(3))

# =========================================================================
print("="*80)
print("VERIFIED SOLUTION CATALOG FOR Q = 3/2 GRID PROBLEM")
print("="*80)
sys.stdout.flush()

# =========================================================================
# TYPE 0: DEGENERATE
# =========================================================================
print("\n" + "="*80)
print("TYPE 0: FULLY DEGENERATE (a_i = b_i)")
print("="*80)

# Use exact circle parametrization. Show a clean example.
theta = np.radians(20)
u = circle_point(theta)
a = list(u**2)
b = list(a)  # identical
print(f"\nExample (theta=20°): a = b = ({a[0]:.10f}, {a[1]:.10f}, {a[2]:.10f})")
print("All 8 triplets identical, Q verification:")
check_grid(a, b, verbose=True)
print("All OK: True (trivially)")
sys.stdout.flush()

# =========================================================================
# TYPE 1: ONE ROW REFLECTED
# =========================================================================
print("\n" + "="*80)
print("TYPE 1: EXACTLY ONE ROW REFLECTED")
print("="*80)

# Use the EXACT circle: any point satisfies Q=3/2 by construction.
# Reflect one row using Vieta's formula.

for j_reflect in [0, 1, 2]:
    row_names = ["1", "2", "3"]
    theta = np.radians(20)
    u = circle_point(theta)
    x = np.abs(u)  # positive square roots

    # Reflect row j_reflect
    others = [i for i in range(3) if i != j_reflect]
    k, l = others
    y_j = 4*(x[k] + x[l]) - x[j_reflect]

    a = list(x**2)
    b = list(x**2)
    b[j_reflect] = y_j**2

    print(f"\n--- Type 1, Row {row_names[j_reflect]} reflected (theta=20°) ---")
    print(f"  sqrt(a) = ({x[0]:.10f}, {x[1]:.10f}, {x[2]:.10f})")
    print(f"  y_{row_names[j_reflect]} = 4*(sqrt(a_{row_names[k]})+sqrt(a_{row_names[l]})) - sqrt(a_{row_names[j_reflect]}) = {y_j:.10f}")
    print(f"  Grid: | {a[0]:.10f}  {b[0]:.10f} |")
    print(f"        | {a[1]:.10f}  {b[1]:.10f} |")
    print(f"        | {a[2]:.10f}  {b[2]:.10f} |")
    print(f"  b_{row_names[j_reflect]}/a_{row_names[j_reflect]} = {b[j_reflect]/a[j_reflect]:.8f}")
    ok, _ = check_grid(a, b, verbose=True)
    print(f"  All OK: {ok}")

sys.stdout.flush()

# =========================================================================
# TYPE 2: TWO ROWS REFLECTED (special discrete solutions)
# =========================================================================
print("\n" + "="*80)
print("TYPE 2: TWO ROWS REFLECTED (special conditions)")
print("="*80)

# From the analysis: for rows i,j reflected with row k kept,
# the doubly-reflected triplet (y_i, y_j, x_k) must also lie on
# a Q=3/2 cone (possibly with different sign pattern).
#
# For rows 1,2 reflected, row 3 kept:
# The condition for (y1, y2, x3) to satisfy the flip-3 cone
# (signs: +,+,-) requires:
#   x1+x2 = c*x3 where c = (-9+sqrt(133))/2 - 1 ≈ 0.2663
#
# And (x1,x2,x3) must lie on the all-positive Q=3/2 cone.

print("\n--- Type 2: Rows 1,2 reflected, row 3 kept ---")
print("Condition: x1+x2 = c*x3 where c+1 = (-9+sqrt(133))/2")

c_plus_1 = (-9+np.sqrt(133))/2
c = c_plus_1 - 1
print(f"  c = {c:.10f}")
print(f"  c+1 = S/x3 = {c_plus_1:.10f}")

# From the scan: two solutions (related by x1<->x2 symmetry)
# t = x1/x3, x2/x3 = c - t
# 6*t^2 - 6*c*t + c^2+1-4*c = 0
disc_inner = 3*c**2 + 24*c - 6
t1 = (3*c + np.sqrt(disc_inner))/6
t2 = (3*c - np.sqrt(disc_inner))/6

for t_val, t_name in [(t1, "t1"), (t2, "t2")]:
    if t_val > 0 and t_val < c:
        x1 = t_val
        x2 = c - t_val
        x3 = 1.0
        y1 = 4*(x2+x3) - x1
        y2 = 4*(x1+x3) - x2

        a = [x1**2, x2**2, x3**2]
        b = [y1**2, y2**2, x3**2]

        print(f"\n  Solution {t_name}:")
        print(f"  sqrt(a) = ({x1:.10f}, {x2:.10f}, {x3:.10f})")
        print(f"  y1 = {y1:.10f}, y2 = {y2:.10f}")
        print(f"  Grid: | {a[0]:.10f}  {b[0]:.10f} |")
        print(f"        | {a[1]:.10f}  {b[1]:.10f} |")
        print(f"        | {a[2]:.10f}  {b[2]:.10f} |")
        ok, _ = check_grid(a, b, verbose=True)
        print(f"  All OK: {ok}")

        if ok:
            # Normalized form
            Sa = sum(a)
            print(f"\n  Normalized (sum a = 1):")
            a_n = [ai/Sa for ai in a]
            b_n = [bi/Sa for bi in b]
            print(f"  a = ({a_n[0]:.10f}, {a_n[1]:.10f}, {a_n[2]:.10f})")
            print(f"  b = ({b_n[0]:.10f}, {b_n[1]:.10f}, {b_n[2]:.10f})")

sys.stdout.flush()

# Also check all 3 choices of which row is kept
print("\n--- Type 2: Checking all row-pair combinations ---")
# For each pair (i,j) reflected, row k kept:
# Need x_i + x_j = c * x_k on the cone.
# The condition c depends on which sign pattern (y_i,y_j,x_k) uses.

# Let me just numerically scan the circle for all Type 2 solutions
print("\nNumerical scan of circle for Type 2 solutions:")
type2_solutions = []

for theta_deg_100 in range(0, 36000):
    theta = np.radians(theta_deg_100 / 100)
    u = circle_point(theta)
    x = np.abs(u)
    if min(x) < 1e-10:
        continue

    for i, j, k in [(0,1,2), (0,2,1), (1,2,0)]:
        y_i = 4*(x[j]+x[k]) - x[i]
        y_j = 4*(x[i]+x[k]) - x[j]

        a = list(x**2)
        b_test = list(x**2)
        b_test[i] = y_i**2
        b_test[j] = y_j**2

        ok, _ = check_grid(a, b_test, verbose=False)
        if ok:
            # Check if genuinely Type 2 (both rows differ)
            if abs(a[i]-b_test[i]) > 1e-8 and abs(a[j]-b_test[j]) > 1e-8:
                # Check if already found (up to permutation)
                ratio_key = tuple(sorted([a[0]/sum(a), a[1]/sum(a), a[2]/sum(a)]))
                is_new = True
                for prev_key in type2_solutions:
                    if all(abs(ratio_key[ii]-prev_key[ii]) < 1e-4 for ii in range(3)):
                        is_new = False
                        break
                if is_new:
                    type2_solutions.append(ratio_key)
                    print(f"\n  theta={theta_deg_100/100:.2f}°, rows {i+1},{j+1} reflected:")
                    print(f"  a = ({a[0]:.10f}, {a[1]:.10f}, {a[2]:.10f})")
                    print(f"  b = ({b_test[0]:.10f}, {b_test[1]:.10f}, {b_test[2]:.10f})")
                    check_grid(a, b_test, verbose=True)

if not type2_solutions:
    print("  No Type 2 solutions found in fine scan!")

sys.stdout.flush()

# =========================================================================
# TYPE 3: THREE ROWS REFLECTED
# =========================================================================
print("\n" + "="*80)
print("TYPE 3: ALL THREE ROWS REFLECTED")
print("="*80)

type3_solutions = []
for theta_deg_100 in range(0, 36000):
    theta = np.radians(theta_deg_100 / 100)
    u = circle_point(theta)
    x = np.abs(u)
    if min(x) < 1e-10:
        continue

    y = [4*(x[1]+x[2])-x[0], 4*(x[0]+x[2])-x[1], 4*(x[0]+x[1])-x[2]]

    a = list(x**2)
    b = [yi**2 for yi in y]

    ok, _ = check_grid(a, b, verbose=False)
    if ok:
        n_diff = sum(1 for ii in range(3) if abs(a[ii]-b[ii]) > 1e-8)
        if n_diff == 3:
            ratio_key = tuple(sorted([a[0]/sum(a), a[1]/sum(a), a[2]/sum(a)]))
            is_new = True
            for prev_key in type3_solutions:
                if all(abs(ratio_key[ii]-prev_key[ii]) < 1e-4 for ii in range(3)):
                    is_new = False
                    break
            if is_new:
                type3_solutions.append(ratio_key)
                print(f"\n  theta={theta_deg_100/100:.2f}°:")
                print(f"  a = ({a[0]:.10f}, {a[1]:.10f}, {a[2]:.10f})")
                print(f"  b = ({b[0]:.10f}, {b[1]:.10f}, {b[2]:.10f})")
                check_grid(a, b, verbose=True)

if not type3_solutions:
    print("No Type 3 solutions found.")

sys.stdout.flush()

# =========================================================================
# ALSO: Search for Type 2/3 with DIFFERENT y_i choices
# =========================================================================
print("\n" + "="*80)
print("EXTENDED SEARCH: Type 2/3 with general y_i")
print("="*80)

# The y_i = 4(x_j+x_k)-x_i is the Vieta reflection for the all-positive cone.
# But one could also use reflections from OTHER sign patterns.
# For the "flip-i" cone: x_i^2+x_j^2+x_k^2 = 4(-x_i*x_j-x_i*x_k+x_j*x_k)
# The reflection of x_i in this cone gives:
# x_i^2 + 4(x_j+x_k)*x_i + (x_j^2+x_k^2-4*x_j*x_k) = 0
# y_i' = -4(x_j+x_k) - x_i  (negative, no good)
# or
# y_i' = -4(x_j+x_k) + x_i  (also likely negative)
# These give |y_i'| = 4(x_j+x_k)+x_i or |y_i'| = |x_i-4(x_j+x_k)|
# The second one is the SAME as the all-positive reflection in absolute value.
# The first gives y_i' = -(4(x_j+x_k)+x_i), so |y_i'| = 4(x_j+x_k)+x_i
# and y_i'^2 = (4(x_j+x_k)+x_i)^2.

# So there are TWO possible reflected masses for each row:
# b_i = (4(x_j+x_k) - x_i)^2   [the "close" reflection]
# b_i = (4(x_j+x_k) + x_i)^2   [the "far" reflection]

# The first we've been using. Let me check if the second also works.

print("\n'Far' reflection: b_i = (4(x_j+x_k) + x_i)^2")
theta = np.radians(20)
u = circle_point(theta)
x = np.abs(u)

# Test far reflection for row 1
y1_far = 4*(x[1]+x[2]) + x[0]
a = list(x**2)
b_far = [y1_far**2, x[1]**2, x[2]**2]
print(f"  x = ({x[0]:.8f}, {x[1]:.8f}, {x[2]:.8f})")
print(f"  y1_far = {y1_far:.8f}")
print(f"  Grid: | {a[0]:.8f}  {b_far[0]:.8f} |")
print(f"        | {a[1]:.8f}  {b_far[1]:.8f} |")
print(f"        | {a[2]:.8f}  {b_far[2]:.8f} |")
ok, _ = check_grid(a, b_far, verbose=True)
print(f"  All OK: {ok}")

# The far reflection comes from the flip-1 cone equation for y1:
# y1^2 + 4*(x2+x3)*y1 + (x2^2+x3^2-4*x2*x3) = 0
# (when put on the flip-1 cone: y1^2+x2^2+x3^2 = 4(-y1*x2-y1*x3+x2*x3))
# Let me verify:
print(f"\nVerify flip-1 cone for far reflection:")
y1f = y1_far
lhs = y1f**2 + x[1]**2 + x[2]**2
rhs = 4*(-y1f*x[1] - y1f*x[2] + x[1]*x[2])
print(f"  LHS = {lhs:.10f}, RHS = {rhs:.10f}")
# This checks if (-y1_far, x2, x3) satisfies the all-positive cone:
lhs2 = y1f**2 + x[1]**2 + x[2]**2
rhs2 = 4*(y1f*x[1] + y1f*x[2] + x[1]*x[2])
print(f"  All-positive with |y1|: LHS = {lhs2:.10f}, RHS = {rhs2:.10f}")
# Neither works! Let me recheck.
# The quadratic for x1 on the all-positive cone:
# x1^2 - 4*(x2+x3)*x1 + (x2^2+x3^2-4*x2*x3) = 0
# The OTHER equation for the flip-1 cone:
# x1^2 + 4*(x2+x3)*x1 + (x2^2+x3^2-4*x2*x3) = 0
# Discriminant: 16*(x2+x3)^2 - 4*(x2^2+x3^2-4*x2*x3)
#             = 12*x2^2 + 48*x2*x3 + 12*x3^2 = 12*(x2+x3)^2 + 24*x2*x3
# Root: x1 = -4*(x2+x3)/2 ± sqrt(disc)/2 = -2*(x2+x3) ± sqrt(3*(x2+x3)^2+6*x2*x3)
# For x1 >= 0: need sqrt(...) >= 2*(x2+x3), so 3*(x2+x3)^2+6*x2*x3 >= 4*(x2+x3)^2
# -1*(x2+x3)^2+6*x2*x3 >= 0, i.e., 6*x2*x3 >= (x2+x3)^2, i.e., 4*x2*x3 >= (x2-x3)^2
# This is always true when x2,x3 > 0! (since (x2-x3)^2 <= (x2+x3)^2 = x2^2+x3^2+2*x2*x3 <= 4*x2*x3 iff x2^2+x3^2 <= 2*x2*x3 iff (x2-x3)^2 <= 0)
# No, that's NOT always true. Only when x2=x3.

# OK so the far reflection is NOT generally on any of the 4 cones.
# Let me verify: does (y1_far, x2, x3) with signs (-1,+1,+1) satisfy Q=3/2?
q_test = Q_val(y1_far**2, x[1]**2, x[2]**2, -1, 1, 1)
print(f"\n  Q(-y1_far, x2, x3) = {q_test:.10f}")
q_test2 = Q_val(y1_far**2, x[1]**2, x[2]**2, 1, 1, 1)
print(f"  Q(+y1_far, x2, x3) = {q_test2:.10f}")

# So the far reflection does NOT work. Only the Vieta reflection y=4(xj+xk)-xi works.

sys.stdout.flush()

# =========================================================================
# Now do a TRULY free numerical search for ANY 2-row or 3-row solution
# =========================================================================
print("\n" + "="*80)
print("FREE NUMERICAL SEARCH (no Vieta constraint)")
print("="*80)

from scipy.optimize import minimize

def grid_obj_6params(params):
    """6 free parameters: sqrt of each mass value."""
    x = [abs(params[0]), abs(params[2]), abs(params[4])]
    y = [abs(params[1]), abs(params[3]), abs(params[5])]
    total = 0.0
    for c1, c2, c3 in product([0,1], repeat=3):
        r = [x[0] if c1==0 else y[0],
             x[1] if c2==0 else y[1],
             x[2] if c3==0 else y[2]]
        R2 = r[0]**2 + r[1]**2 + r[2]**2
        if R2 < 1e-30:
            continue
        target = np.sqrt(1.5 * R2)
        sums = [abs(r[0]+r[1]+r[2]), abs(r[0]+r[1]-r[2]),
                abs(r[0]-r[1]+r[2]), abs(-r[0]+r[1]+r[2])]
        min_dev = min((s - target)**2 for s in sums)
        total += min_dev
    return total

# Many random restarts with local optimization
np.random.seed(42)
solutions = []

for trial in range(500):
    x0 = np.random.uniform(0.01, 3, 6)
    result = minimize(grid_obj_6params, x0, method='Nelder-Mead',
                     options={'maxiter': 10000, 'xatol': 1e-14, 'fatol': 1e-14})
    if result.fun < 1e-12:
        params = result.x
        x = [abs(params[0]), abs(params[2]), abs(params[4])]
        y = [abs(params[1]), abs(params[3]), abs(params[5])]
        a = [xi**2 for xi in x]
        b = [yi**2 for yi in y]

        n_diff = sum(1 for i in range(3) if abs(a[i]-b[i]) > 1e-4 * max(a[i],b[i],0.001))

        # Normalize
        Sa = sum(a)
        key = tuple(sorted([a[0]/Sa, a[1]/Sa, a[2]/Sa]))
        key_b = tuple(sorted([b[0]/Sa, b[1]/Sa, b[2]/Sa]))

        is_new = True
        for prev_key, prev_key_b, prev_n in solutions:
            if prev_n == n_diff:
                if (all(abs(key[ii]-prev_key[ii]) < 1e-3 for ii in range(3)) and
                    all(abs(key_b[ii]-prev_key_b[ii]) < 1e-3 for ii in range(3))):
                    is_new = False
                    break
                if (all(abs(key[ii]-prev_key_b[ii]) < 1e-3 for ii in range(3)) and
                    all(abs(key_b[ii]-prev_key[ii]) < 1e-3 for ii in range(3))):
                    is_new = False
                    break

        if is_new:
            solutions.append((key, key_b, n_diff))
            if n_diff >= 2:
                print(f"\n*** {n_diff}-ROW SOLUTION found (trial {trial}) ***")
                print(f"  a = ({a[0]:.8f}, {a[1]:.8f}, {a[2]:.8f})")
                print(f"  b = ({b[0]:.8f}, {b[1]:.8f}, {b[2]:.8f})")
                print(f"  obj = {result.fun:.2e}")
                ok, _ = check_grid(a, b, verbose=True)

print(f"\nTotal solutions found: {len(solutions)}")
print(f"  Type 0: {sum(1 for _,_,n in solutions if n==0)}")
print(f"  Type 1: {sum(1 for _,_,n in solutions if n==1)}")
print(f"  Type 2: {sum(1 for _,_,n in solutions if n==2)}")
print(f"  Type 3: {sum(1 for _,_,n in solutions if n==3)}")

sys.stdout.flush()

# =========================================================================
# SCALING: Physical mass example
# =========================================================================
print("\n" + "="*80)
print("PHYSICAL SCALING EXAMPLE")
print("="*80)

# Use exact Q=3/2 circle point matching ratios closest to (0.122, 1.70, 3.64)
target = np.array(sorted([0.122, 1.70, 3.64]))

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
m_sorted = np.sort(u_best**2)
S_best = np.sum(target) / np.sum(m_sorted)
M0 = S_best
col_a = m_sorted * M0

print(f"Best theta = {np.degrees(best_theta):.3f}°")
print(f"M0 (scale) = {M0:.8f}")
print(f"Column-1 masses (a_i = m_i * M0):")
print(f"  a_1 = {col_a[0]:.6f} GeV  (target: {target[0]} GeV)")
print(f"  a_2 = {col_a[1]:.6f} GeV  (target: {target[1]} GeV)")
print(f"  a_3 = {col_a[2]:.6f} GeV  (target: {target[2]} GeV)")

# Column 2 for each row reflection
x_phys = np.sqrt(col_a)
print(f"\nColumn-2 masses for each Type 1 reflection:")
for j in range(3):
    k, l = [i for i in range(3) if i != j]
    y_j = 4*(x_phys[k] + x_phys[l]) - x_phys[j]
    b_j = y_j**2
    print(f"  Reflect row {j+1}: b_{j+1} = {b_j:.6f} GeV  (ratio b/a = {b_j/col_a[j]:.4f})")

# Verify with exact circle
print(f"\nVerification (row 1 reflected):")
y1 = 4*(x_phys[1]+x_phys[2]) - x_phys[0]
grid_a = list(col_a)
grid_b = [y1**2, col_a[1], col_a[2]]
check_grid(grid_a, grid_b, verbose=True)

# Also verify row 3 reflected (smallest ratio)
print(f"\nVerification (row 3 reflected):")
y3 = 4*(x_phys[0]+x_phys[1]) - x_phys[2]
grid_b3 = [col_a[0], col_a[1], y3**2]
check_grid(grid_a, grid_b3, verbose=True)

sys.stdout.flush()

# =========================================================================
# SUMMARY TABLE
# =========================================================================
print("\n" + "="*80)
print("SUMMARY TABLE: ALL SOLUTION TYPES")
print("="*80)

print("""
+-------+-------------------+-------------------+-------------------+
| Type  | # rows different  | # free parameters | Existence         |
+-------+-------------------+-------------------+-------------------+
|   0   |        0          |   2 (θ, M0)       | Continuous family |
|   1   |        1          |   2 (θ, M0)       | Continuous family |
|       |                   |   x 3 (row choice) | (3 sub-families)  |
|   2   |        2          |   1 (M0 only)      | Discrete set      |
|       |                   |   x 3 (pair choice) | (isolated points) |
|   3   |        3          |        -           | No solutions      |
+-------+-------------------+-------------------+-------------------+

Key formulas:
  Q = 3/2 cone:  m1+m2+m3 = 4*(sqrt(m1*m2)+sqrt(m1*m3)+sqrt(m2*m3))
  Type 1 reflection:  sqrt(b_j) = 4*(sqrt(a_k)+sqrt(a_l)) - sqrt(a_j)
  Type 2 condition:   sqrt(a_i)+sqrt(a_j) = c*sqrt(a_k)  [c depends on sign pattern]

Symmetries:
  - R+ scaling (m -> λm)
  - S_3 row permutations
  - Z_2 column swap
  - Z_2 per-row column swap (for Type 1)
  - U(1) circle rotation (for Types 0,1)
  - O(2,1) Lorentz symmetry of the cone
""")

print("=== COMPLETE ===")
