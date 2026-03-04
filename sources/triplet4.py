#!/usr/bin/env python3
"""
Part 4: Complete classification of Q=3/2 grid solutions.
Search for 2-row and 3-row non-degenerate solutions, including the special
condition S = 3(sqrt(3)-1)*x3 for 2-row solutions on the all-positive cone.
Also do a free numerical search for general solutions.
"""
import numpy as np
from itertools import product
from scipy.optimize import minimize, differential_evolution

def Q_val(m1, m2, m3, s1=1, s2=1, s3=1):
    denom = m1 + m2 + m3
    if denom < 1e-30:
        return float('inf')
    num = (s1*np.sqrt(abs(m1)) + s2*np.sqrt(abs(m2)) + s3*np.sqrt(abs(m3)))**2
    return num / denom

def check_grid(a, b, verbose=False):
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
            if abs(q - 1.5) < 1e-6:
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
                  f"Q={r['best_q']:.10f} signs={r['best_signs']} [{status}]")
    return all_ok, results

def circle_point(theta):
    return np.array([1,1,1])/np.sqrt(6) + np.cos(theta)*np.array([1,-1,0])/2 + \
           np.sin(theta)*np.array([1,1,-2])/(2*np.sqrt(3))

# =============================================================================
# PART A: 2-row solutions with special condition
# =============================================================================
print("="*80)
print("PART A: 2-ROW SOLUTIONS (rows 1,2 differ; row 3 same)")
print("="*80)

# Condition: x1+x2 = (3sqrt(3)-4)*x3 for all-positive cone compatibility
# of the reflected pair (y1,y2,x3).
# Equivalently: S = x1+x2+x3 = (3sqrt(3)-3)*x3

# Also need to check the other sign patterns.
# Let me scan the circle for points satisfying x1+x2 = c*x3 for various c values.

# For the Q=3/2 circle: u(θ) = (u1,u2,u3).
# x_i = |u_i| (we need all positive for the all-+ cone).
# Range of θ where all u_i > 0:
# u1 > 0 and u2 > 0 and u3 > 0.

# Find the ranges
thetas = np.linspace(0, 2*np.pi, 10000)
all_pos = []
for theta in thetas:
    u = circle_point(theta)
    if u[0] > 0 and u[1] > 0 and u[2] > 0:
        all_pos.append(theta)

if all_pos:
    print(f"All positive range: θ ∈ [{np.degrees(all_pos[0]):.1f}°, {np.degrees(all_pos[-1]):.1f}°]")

# In this range, x_i = u_i. The cone is automatically satisfied.
# Reflection: y_i = 4*(x_j+x_k) - x_i

# For rows 1,2 to differ: need (y1,y2,x3) on some cone.
# Check the condition for all-+ cone: x1+x2 = (3sqrt(3)-4)*x3

target_ratio_allplus = 3*np.sqrt(3) - 4
print(f"\nTarget x1+x2 / x3 for all-+ cone: {target_ratio_allplus:.8f}")

# Also for flip-3 cone: S² + 6*S*x3 - 18*x3² = 0 gives S = (3sqrt(3)-3)*x3
# Wait, I derived this for all-+ cone. Let me check each cone.

# For flip-3 cone: derived condition S²+9*S*x3 - 13*x3² = 0
# S = (-9+sqrt(133))/2 * x3
target_ratio_flip3 = (-9+np.sqrt(133))/2  # This is S/x3
target_x12_flip3 = target_ratio_flip3 - 1  # (x1+x2)/x3

print(f"Target x1+x2 / x3 for flip-3 cone: {target_x12_flip3:.8f}")
print(f"  (S/x3 = {target_ratio_flip3:.8f})")

# Let me also derive for flip-1 and flip-2 cones.
# Flip-1: ε12=-1, ε13=-1, ε23=+1
# (y1,y2,x3) on flip-1 cone:
# y1²+y2²+x3² = 4*(-y1*y2 - y1*x3 + y2*x3)
# With y1 = 4S-5x1, y2 = 4S-5x2, S=x1+x2+x3:
# This needs a separate derivation. Let me just check numerically.

# Scan the circle for points where the 2-row reflection works
print("\nScanning for 2-row reflection solutions...")

for theta_deg_10 in range(0, 3600):
    theta = np.radians(theta_deg_10 / 10)
    u = circle_point(theta)
    x1, x2, x3 = abs(u[0]), abs(u[1]), abs(u[2])
    if min(x1,x2,x3) < 1e-10:
        continue

    y1 = abs(4*(x2+x3) - x1)
    y2 = abs(4*(x1+x3) - x2)
    y3 = abs(4*(x1+x2) - x3)

    # Test all 3 pairs of rows
    for pair_name, b_test in [
        ("rows 1,2", [y1**2, y2**2, x3**2]),
        ("rows 1,3", [y1**2, x2**2, y3**2]),
        ("rows 2,3", [x1**2, y2**2, y3**2]),
    ]:
        a_test = [x1**2, x2**2, x3**2]
        ok, _ = check_grid(a_test, b_test, verbose=False)
        if ok:
            print(f"\n*** 2-ROW at θ={theta_deg_10/10:.1f}°, {pair_name} ***")
            print(f"  a = ({a_test[0]:.8f}, {a_test[1]:.8f}, {a_test[2]:.8f})")
            print(f"  b = ({b_test[0]:.8f}, {b_test[1]:.8f}, {b_test[2]:.8f})")
            ratio = (x1+x2)/x3
            print(f"  (x1+x2)/x3 = {ratio:.8f}")
            check_grid(a_test, b_test, verbose=True)

# Also try reflection with OTHER points on the circle for y_i
# The reflection y_i = 4(x_j+x_k) - x_i assumes the reflected point stays on
# the all-+ cone. But y_i could also be the OTHER root of the quadratic.
# Let's try: for each pair of rows, scan independently for compatible y values.

print("\n" + "="*80)
print("PART B: FREE NUMERICAL SEARCH FOR 2-ROW AND 3-ROW SOLUTIONS")
print("="*80)

def grid_objective(params):
    """Optimize over 6 square roots."""
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

print("\nDifferential evolution search (20 seeds)...")
all_solutions = []

for seed in range(30):
    bounds = [(0.01, 5)] * 6
    result = differential_evolution(grid_objective, bounds, seed=seed, maxiter=2000,
                                     tol=1e-16, popsize=40, mutation=(0.5,1.5))
    if result.fun < 1e-10:
        params = result.x
        x = [abs(params[0]), abs(params[2]), abs(params[4])]
        y = [abs(params[1]), abs(params[3]), abs(params[5])]
        a = [xi**2 for xi in x]
        b = [yi**2 for yi in y]

        # Check how many rows actually differ
        n_diff = sum(1 for i in range(3) if abs(a[i]-b[i]) > 1e-6 * max(a[i],b[i],1))

        # Normalize: divide by sum of column a
        Sa = sum(a)
        a_norm = [ai/Sa for ai in a]
        b_norm = [bi/Sa for bi in b]

        # Check if this is a new solution (not a permutation of existing)
        is_new = True
        for prev_a, prev_b, prev_nd in all_solutions:
            if prev_nd == n_diff:
                # Check if same up to row permutations and column swaps
                for perm in [(0,1,2),(0,2,1),(1,0,2),(1,2,0),(2,0,1),(2,1,0)]:
                    pa = [a_norm[p] for p in perm]
                    pb = [b_norm[p] for p in perm]
                    if (all(abs(pa[i]-prev_a[i]) < 1e-4 for i in range(3)) and
                        all(abs(pb[i]-prev_b[i]) < 1e-4 for i in range(3))):
                        is_new = False
                        break
                    # Also try column swap
                    if (all(abs(pb[i]-prev_a[i]) < 1e-4 for i in range(3)) and
                        all(abs(pa[i]-prev_b[i]) < 1e-4 for i in range(3))):
                        is_new = False
                        break

        if is_new:
            all_solutions.append((a_norm, b_norm, n_diff))
            print(f"\n  Seed {seed}: {n_diff} rows differ, obj={result.fun:.2e}")
            print(f"    a = ({a[0]:.6f}, {a[1]:.6f}, {a[2]:.6f})")
            print(f"    b = ({b[0]:.6f}, {b[1]:.6f}, {b[2]:.6f})")
            print(f"    a_norm = ({a_norm[0]:.6f}, {a_norm[1]:.6f}, {a_norm[2]:.6f})")
            print(f"    b_norm = ({b_norm[0]:.6f}, {b_norm[1]:.6f}, {b_norm[2]:.6f})")
            # Verify
            ok, _ = check_grid(a, b, verbose=True)

print(f"\n\nTotal distinct solutions found: {len(all_solutions)}")
print(f"  0-row (degenerate): {sum(1 for _,_,n in all_solutions if n==0)}")
print(f"  1-row: {sum(1 for _,_,n in all_solutions if n==1)}")
print(f"  2-row: {sum(1 for _,_,n in all_solutions if n==2)}")
print(f"  3-row: {sum(1 for _,_,n in all_solutions if n==3)}")

# =============================================================================
# PART C: Analytical understanding of 1-row solutions
# =============================================================================
print("\n" + "="*80)
print("PART C: COMPLETE ANALYTICAL DESCRIPTION")
print("="*80)

# Result so far: only 1-row solutions and 0-row (degenerate) solutions exist.
# (No 2-row or 3-row solutions found.)

# A 1-row solution has this structure (WLOG, row 1 differs):
# Column a: (x1², x2², x3²) where (x1,x2,x3) satisfies Q=3/2 cone
# Column b: ((4(x2+x3)-x1)², x2², x3²) same row 2,3

# The reflection formula: y1 = 4(x2+x3) - x1
# This comes from the quadratic x1² - 4(x2+x3)*x1 + (x2²+x3²-4x2x3) = 0
# which has roots x1 and y1 with x1 + y1 = 4(x2+x3), x1*y1 = x2²+x3²-4x2x3.

# So: a1*b1 = (x1*y1)² = (x2²+x3²-4x2x3)² = (a2+a3-4*sqrt(a2*a3))²
# And: sqrt(a1)+sqrt(b1) = x1+y1 = 4(x2+x3) = 4(sqrt(a2)+sqrt(a3))
#   OR sqrt(a1)+sqrt(b1) = x1+|y1| when y1>0
#   OR sqrt(a1)-sqrt(b1) = x1-|y1| = 4(x2+x3)-2y1... hmm

# Let me express in terms of the original masses.
# With a2=m2, a3=m3 fixed, the two column-1 values satisfy:
# sqrt(a1) + sqrt(b1) = 4*(sqrt(m2) + sqrt(m3))
# sqrt(a1) * sqrt(b1) = |m2 + m3 - 4*sqrt(m2*m3)|

# Note: m2+m3-4*sqrt(m2*m3) = (sqrt(m2)-sqrt(m3))^2 - 2*sqrt(m2*m3)*(2-1)
# Hmm, not so clean. Let me just verify the product formula.

print("\nReflection formulas (row j differs, rows k,l same):")
print("  sqrt(a_j) + sqrt(b_j) = 4*(sqrt(m_k) + sqrt(m_l))")
print("  sqrt(a_j) * sqrt(b_j) = |m_k + m_l - 4*sqrt(m_k*m_l)|")
print("  (where the product can be negative if sqrt(mk)+sqrt(ml) is too large)")
print()

# Verify with a specific example
theta = np.radians(20)
u = circle_point(theta)
x1, x2, x3 = u[0], u[1], u[2]
y1 = 4*(x2+x3) - x1
print(f"Example: θ=20°, x=({x1:.6f},{x2:.6f},{x3:.6f})")
print(f"  y1 = {y1:.6f}")
print(f"  x1+y1 = {x1+y1:.6f}, 4(x2+x3) = {4*(x2+x3):.6f}")
print(f"  x1*y1 = {x1*y1:.6f}, x2²+x3²-4x2x3 = {x2**2+x3**2-4*x2*x3:.6f}")
print(f"  a1 = {x1**2:.8f}, b1 = {y1**2:.8f}")
print(f"  a1*b1 = {x1**2*y1**2:.8f}, (x2²+x3²-4x2x3)² = {(x2**2+x3**2-4*x2*x3)**2:.8f}")

# Product: x1*y1 = x2² + x3² - 4*x2*x3 (from Vieta)
# = (x2-x3)² - 2*x2*x3
# For x2,x3 > 0 and close: this can be negative.
# When it's negative, y1 < 0 (since x1 > 0), meaning b1 = y1² but with flipped sign.

print(f"\n  x1*y1 = {x1*y1:.6f} ({'positive' if x1*y1>=0 else 'NEGATIVE'})")
print(f"  This means y1 is {'same sign as x1' if y1>=0 else 'OPPOSITE sign to x1'}")

# When does x1*y1 = x2²+x3²-4x2x3 < 0?
# x2²+x3²-4x2x3 < 0
# Let r = x2/x3: r² - 4r + 1 < 0 => r in (2-sqrt(3), 2+sqrt(3))
# i.e., k- < x2/x3 < k+

km = 2 - np.sqrt(3)
kp = 2 + np.sqrt(3)
print(f"\n  x1*y1 < 0  iff  k- < x2/x3 < k+  (i.e., {km:.6f} < x2/x3 < {kp:.6f})")
print(f"  In this range, the reflection y1 has opposite sign to x1.")
print(f"  This means the 'b' column triplet uses a different sign pattern than 'a' column.")

# =============================================================================
# PART D: Solution catalog and scaling
# =============================================================================
print("\n" + "="*80)
print("PART D: SOLUTION CATALOG")
print("="*80)

print("""
SOLUTION TYPE 0 (Fully Degenerate): a_i = b_i for all i
  Grid: | m1  m1 |    where (m1,m2,m3) satisfy Q=3/2 for some signs.
        | m2  m2 |
        | m3  m3 |
  This is a 2-parameter family (on the unit sphere): the Q=3/2 circle.
  Parametrized by angle θ on the circle:
    u(θ) = (1,1,1)/√6 + cosθ*(1,-1,0)/2 + sinθ*(1,1,-2)/(2√3)
  Masses: m_i = R² * u_i(θ)² for arbitrary scale R > 0.
  Since all 8 triplets are identical, any such configuration trivially works.

SOLUTION TYPE 1 (One row differs): exactly one pair (a_j, b_j) has a_j ≠ b_j.
  WLOG j=1. Start from (x1,x2,x3) on the Q=3/2 cone.
  Grid: | x1²               [4(x2+x3)-x1]²  |
        | x2²               x2²               |
        | x3²               x3²               |

  The reflection formula:
    y1 = 4*(x2+x3) - x1  (from Vieta's formulas for the cone quadratic in x1)

  Key identity: x1 + y1 = 4(x2+x3), x1*y1 = x2² + x3² - 4*x2*x3

  The 4 triplets with c1=0 give (x1², x2², x3²) -> Q=3/2 by construction.
  The 4 triplets with c1=1 give (y1², x2², x3²) -> Q=3/2 since y1 is the
  other root of the same quadratic.

  This is again a 2-parameter family (θ on circle + scale R), but with the
  ADDITIONAL choice of which row to reflect (3 choices: j=1,2,3).

  The b-column mass in the reflected row can be much larger than the a-column mass:
    b1/a1 = y1²/x1² = [4(x2+x3)/x1 - 1]²
  When x1 is small (near a zero of u1), y1 ≈ 4(x2+x3) is large, so b1 >> a1.
  When x1 is large, y1 ≈ -x1 + 4(x2+x3) can be small, so b1 << a1.

SOLUTION TYPE 2 (Two rows differ): NOT found in any search.
  The analytical condition for rows 1,2 to be simultaneously reflectable
  requires x1+x2 = (3√3-4)*x3 ≈ 1.196*x3 (for all-+ cone compatibility),
  which is a measure-zero set. But even this does not guarantee all 8 triplets work
  because the 4 mixed triplets need DIFFERENT sign patterns.

SOLUTION TYPE 3 (Three rows differ): NOT found in any search.
  Even more constrained; appears impossible.
""")

# =============================================================================
# PART E: Verify and display specific examples with physical scaling
# =============================================================================
print("="*80)
print("PART E: PHYSICAL SCALING AND SPECIFIC EXAMPLES")
print("="*80)

# Target: column-1 masses ≈ (0.122, 1.70, 3.64) GeV
# These look like they could be strange quark, charm, and some other particle masses

target_masses = [0.122, 1.70, 3.64]

# For a Q=3/2 solution, we need sqrt(m1)+sqrt(m2)+sqrt(m3) to satisfy the cone.
# Condition: sum(mi) = 4*sum_{i<j}(sqrt(mi*mj)) ... no wait.
# Cone: (sum sqrt(mi))^2 = (3/2)*sum(mi) for all-positive signs.
# Check:
m1, m2, m3 = target_masses
S = np.sqrt(m1) + np.sqrt(m2) + np.sqrt(m3)
R2 = m1 + m2 + m3
Q_target = S**2 / R2
print(f"Target masses: ({m1}, {m2}, {m3})")
print(f"Q with all positive signs: {Q_target:.6f}")
print(f"  (need 1.5, got {Q_target:.6f})")
print()

# These don't satisfy Q=3/2. We need to find masses that DO.
# From the circle parametrization, the masses are proportional to u_i(θ)².
# We need to find θ and R such that R²*u_i(θ)² = target_i.

# This means: u_i(θ)² ∝ target_i.
# On the circle: u_i(θ)² are specific functions of θ.
# We can try to fit, but the ratios are constrained.

# Let's see what ratios the circle gives:
print("Looking for circle point with mass ratios close to target...")
print(f"Target ratios: m2/m1 = {m2/m1:.4f}, m3/m1 = {m3/m1:.4f}")

best_match = None
best_err = float('inf')
for theta_deg_10 in range(0, 3600):
    theta = np.radians(theta_deg_10 / 10)
    u = circle_point(theta)
    masses = u**2
    # Sort to match ordering
    ms = sorted(masses)
    if ms[0] < 1e-12:
        continue
    r2 = ms[1]/ms[0]
    r3 = ms[2]/ms[0]
    err = (r2 - m2/m1)**2 + (r3 - m3/m1)**2
    if err < best_err:
        best_err = err
        best_match = (theta_deg_10/10, u, masses)

if best_match:
    theta_deg, u, masses = best_match
    ms = sorted(masses)
    scale = m1 / ms[0]
    scaled = [mi * scale for mi in sorted(masses)]
    print(f"\nBest circle match: θ = {theta_deg:.1f}°")
    print(f"  Circle masses (unsorted): ({masses[0]:.6f}, {masses[1]:.6f}, {masses[2]:.6f})")
    print(f"  Sorted: ({ms[0]:.6f}, {ms[1]:.6f}, {ms[2]:.6f})")
    print(f"  Ratios: {ms[1]/ms[0]:.4f}, {ms[2]/ms[0]:.4f}")
    print(f"  Scaled to match m1={m1}: ({scaled[0]:.4f}, {scaled[1]:.4f}, {scaled[2]:.4f})")
    print(f"  Target: ({m1}, {m2}, {m3})")
    print(f"  Residual error: {best_err:.6f}")

# The target ratios are m2/m1 = 13.93, m3/m1 = 29.84.
# From the circle, what's the max ratio range?
print("\nRange of mass ratios on Q=3/2 circle:")
ratios = []
for theta_deg_10 in range(0, 3600):
    theta = np.radians(theta_deg_10 / 10)
    u = circle_point(theta)
    masses = sorted(u**2, reverse=True)
    if masses[2] > 1e-10:
        ratios.append((masses[0]/masses[2], theta_deg_10/10))
max_ratio = max(ratios, key=lambda x: x[0])
print(f"Max ratio max(m)/min(m) = {max_ratio[0]:.2f} at θ={max_ratio[1]:.1f}°")

# The target ratio 3.64/0.122 = 29.8, while the max on the circle is much larger.
# But the INTERMEDIATE ratio might not match. Let me check more carefully.

# Actually, the user says "scale M0 so that column-1 masses approximate (0.122, 1.70, 3.64)".
# These look like down, strange, charm quark masses or similar.
# Let me search allowing negative u_i (which just means sign flip).

print("\n\nSearching with sign flips allowed...")
best_match2 = None
best_err2 = float('inf')
target_sorted = sorted(target_masses)

for theta_deg_10 in range(0, 3600):
    theta = np.radians(theta_deg_10 / 10)
    u = circle_point(theta)
    masses = sorted(u**2)  # smallest to largest
    if masses[0] < 1e-12:
        continue
    # Scale to match
    scale = target_sorted[0] / masses[0]
    scaled = [mi * scale for mi in masses]
    err = sum((scaled[i] - target_sorted[i])**2 for i in range(3))
    if err < best_err2:
        best_err2 = err
        best_match2 = (theta_deg_10/10, u, masses, scale, scaled)

if best_match2:
    theta_deg, u, masses, scale, scaled = best_match2
    print(f"Best match: θ = {theta_deg:.1f}°")
    print(f"  Circle u: ({u[0]:.6f}, {u[1]:.6f}, {u[2]:.6f})")
    print(f"  Circle masses (sorted): ({masses[0]:.6f}, {masses[1]:.6f}, {masses[2]:.6f})")
    print(f"  Scale factor: {scale:.4f}")
    print(f"  Scaled: ({scaled[0]:.4f}, {scaled[1]:.4f}, {scaled[2]:.4f})")
    print(f"  Target: ({target_sorted[0]:.4f}, {target_sorted[1]:.4f}, {target_sorted[2]:.4f})")
    print(f"  Error: {best_err2:.6f}")

    # Now compute column 2 by reflection (for each row)
    x = [np.sqrt(s) for s in scaled]
    for j in range(3):
        k, l = [i for i in range(3) if i != j]
        yj = 4*(x[k]+x[l]) - x[j]
        bj = yj**2
        print(f"\n  Reflecting row {j+1}: b_{j+1} = {bj:.6f}")
        print(f"    Ratio b_{j+1}/a_{j+1} = {bj/scaled[j]:.4f}")

# =============================================================================
# PART F: Try to use the actual problem hint with k+, k-
# =============================================================================
print("\n" + "="*80)
print("PART F: SOLUTIONS USING k+ AND k-")
print("="*80)

# The hint suggests triplet (k-, 1, k+) should satisfy Q=3/2.
# We showed it doesn't with the standard Q formula.
# But let me reconsider: maybe the problem uses a DIFFERENT sign convention.

# The problem says "Q(m1,m2,m3) = (√m1 + √m2 + √m3)² / (m1 + m2 + m3)"
# and "allowing sign flips on √m". So √m means we can choose ±√m for each.

# For (k-, 1, k+) with k- = 2-√3, k+ = 2+√3:
# √(k-) = (√3-1)/√2, √(k+) = (1+√3)/√2
# With all positive: (√(k-) + 1 + √(k+))² = (√6 + 1)² = 6 + 2√6 + 1 = 7 + 2√6
# Sum = k- + 1 + k+ = 5
# Q = (7+2√6)/5 ≈ 2.38 ≠ 3/2

# With sign flip on √(k-): (-√(k-) + 1 + √(k+))² = (√2 + 1)²  [since √(k+)-√(k-) = √2]
# Wait: -√(k-) + √(k+) = √2, so -√(k-) + 1 + √(k+) = 1 + √2
# Q = (1+√2)²/5 = (3+2√2)/5 ≈ 1.166 ≠ 3/2

# What about sign flip on 1?
# √(k-) - 1 + √(k+) = √6 - 1
# Q = (√6-1)²/5 = (7-2√6)/5 ≈ 0.420 ≠ 3/2

# Hmm, none work for (k-, 1, k+) with standard masses.

# WAIT. Let me re-read the problem: "the triplet (k⁻, 1, k⁺) satisfies Q = 3/2 with all positive roots"
# Let me check: maybe it's a different Q or different k values?

# Actually, let me reconsider the problem statement more carefully.
# "Q(m1,m2,m3) = (√m1 + √m2 + √m3)² / (m1 + m2 + m3) = 3/2"
# With m1=k-, m2=1, m3=k+:
# Q = (√(k-) + 1 + √(k+))² / (k- + 1 + k+) = (√6 + 1)² / 5 ≈ 2.38

# This is NOT 3/2. So the problem's assertion seems wrong... unless I'm computing wrong.
# Let me be very careful:
kp_exact = 2 + np.sqrt(3)
km_exact = 2 - np.sqrt(3)
skp = np.sqrt(kp_exact)
skm = np.sqrt(km_exact)
print(f"k+ = {kp_exact:.10f}")
print(f"k- = {km_exact:.10f}")
print(f"sqrt(k+) = {skp:.10f}")
print(f"sqrt(k-) = {skm:.10f}")
print(f"sqrt(k+) + sqrt(k-) = {skp+skm:.10f} = sqrt(6) = {np.sqrt(6):.10f}")
print(f"Sum of k's + 1 = {km_exact + 1 + kp_exact:.10f}")
print(f"(sqrt(k-)+1+sqrt(k+))^2 = {(skm+1+skp)**2:.10f}")
print(f"Q_all_pos = {(skm+1+skp)**2 / (km_exact+1+kp_exact):.10f}")

# Let me try k- = (2-sqrt(3))^2 = 7-4sqrt(3), k+ = (2+sqrt(3))^2 = 7+4sqrt(3)
kp2 = (2+np.sqrt(3))**2
km2 = (2-np.sqrt(3))**2
print(f"\nTrying k+ = (2+sqrt(3))^2 = {kp2:.6f}, k- = (2-sqrt(3))^2 = {km2:.6f}")
Q_test = (np.sqrt(km2) + 1 + np.sqrt(kp2))**2 / (km2 + 1 + kp2)
print(f"Q = {Q_test:.10f}")
# sqrt(k-) = 2-sqrt(3), sqrt(k+) = 2+sqrt(3)
# sum sqrt = (2-sqrt(3)) + 1 + (2+sqrt(3)) = 5
# sum m = (7-4sqrt(3)) + 1 + (7+4sqrt(3)) = 15
# Q = 25/15 = 5/3 ≠ 3/2

# Hmm. Let me try a GENERAL approach: given Q=3/2, what ratio m3/m1 with m2=1?
# (sqrt(m1) + 1 + sqrt(m3))^2 = (3/2)*(m1+1+m3)
# Let x = sqrt(m1), z = sqrt(m3):
# (x+1+z)^2 = (3/2)*(x^2+1+z^2)
# 2*(x+1+z)^2 = 3*(x^2+1+z^2)
# 2x^2+2+2z^2+4x+4z+4xz = 3x^2+3+3z^2
# x^2+z^2+1 = 4x+4z+4xz
# x^2 - 4x(1+z) + z^2+1-4z = 0
# x = 2(1+z) ± sqrt(4(1+z)^2 - z^2-1+4z) = 2(1+z) ± sqrt(3z^2+12z+3)
# = 2(1+z) ± sqrt(3(z^2+4z+1)) = 2(1+z) ± sqrt(3)*sqrt((z+2)^2-3)

# For z=1: x = 2*2 ± sqrt(3*6) = 4 ± 3sqrt(2)
# x = 4-3sqrt(2) ≈ -0.243 (negative, not valid with all + signs)
# x = 4+3sqrt(2) ≈ 8.243
# So m1 = (4+3sqrt(2))^2 ≈ 67.94 when m2=1, m3=1.
# That's very asymmetric!

# Let's try z = k- or z = k+:
for z_name, z_val in [("1", 1.0), ("k-", km_exact), ("k+", kp_exact)]:
    disc = 3*(z_val**2 + 4*z_val + 1)
    if disc >= 0:
        x1 = 2*(1+z_val) - np.sqrt(disc)
        x2 = 2*(1+z_val) + np.sqrt(disc)
        print(f"\nm2=1, m3={z_name}={z_val:.6f}:")
        print(f"  sqrt(m1) = {x1:.6f} or {x2:.6f}")
        print(f"  m1 = {x1**2:.6f} or {x2**2:.6f}")
        if x1 > 0:
            print(f"  Triplet: ({x1**2:.6f}, 1, {z_val:.6f})")
            Q_check = (x1+1+z_val)**2/(x1**2+1+z_val**2)
            print(f"  Q_check = {Q_check:.10f}")
        if x2 > 0:
            print(f"  Triplet: ({x2**2:.6f}, 1, {z_val:.6f})")
            Q_check = (x2+1+z_val)**2/(x2**2+1+z_val**2)
            print(f"  Q_check = {Q_check:.10f}")

# So with m2=1, m3=1: m1 = (4+3√2)² ≈ 67.94. Not a nice number.
# With m2=1, m3=k-: sqrt(m1) = 2(1+sqrt(k-))-sqrt(3*(k-+4*sqrt(k-)+1))

# Let me try to find solutions with m3/m1 = k+ or similar patterns
print("\n\nSearching for solutions with special ratios...")
for z_val in np.linspace(0.01, 5, 1000):
    disc = 3*(z_val**2 + 4*z_val + 1)
    if disc < 0:
        continue
    x1 = 2*(1+z_val) - np.sqrt(disc)
    if x1 > 1e-6:
        m1 = x1**2
        m3 = z_val**2
        # Check if m1, 1, m3 have nice ratios
        if abs(m3/m1 - kp_exact) < 0.01:
            print(f"  m3/m1 ≈ k+: z={z_val:.6f}, m1={m1:.6f}, m3={m3:.6f}, ratio={m3/m1:.6f}")
        if abs(m1 - km_exact) < 0.01:
            print(f"  m1 ≈ k-: z={z_val:.6f}, m1={m1:.6f}, m3={m3:.6f}")
        if abs(m3 - kp_exact) < 0.01:
            print(f"  m3 ≈ k+: z={z_val:.6f}, m1={m1:.6f}, m3={m3:.6f}")

print("\nDone with analysis.")
