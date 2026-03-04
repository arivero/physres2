#!/usr/bin/env python3
"""
Focused analysis of Q=3/2 grid problem.
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

def check_grid_full(a, b, verbose=False):
    """Check all 8 triplets, return detailed info."""
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
        results.append({
            'choices': (c1,c2,c3), 'masses': m, 'found': found,
            'best_q': best_q, 'best_signs': best_signs
        })
    if verbose:
        for r in results:
            status = "OK" if r['found'] else "FAIL"
            print(f"  ({r['choices'][0]},{r['choices'][1]},{r['choices'][2]}): "
                  f"m=({r['masses'][0]:.8f},{r['masses'][1]:.8f},{r['masses'][2]:.8f}) "
                  f"Q={r['best_q']:.10f} signs={r['best_signs']} [{status}]")
    return all_ok, results

def grid_objective(params):
    """Objective: sum of min squared deviations from Q=3/2 over all 8 triplets."""
    x1, y1, x2, y2, x3, y3 = params
    total = 0.0
    for c1, c2, c3 in product([0,1], repeat=3):
        r = [abs(x1) if c1==0 else abs(y1),
             abs(x2) if c2==0 else abs(y2),
             abs(x3) if c3==0 else abs(y3)]
        R2 = r[0]**2 + r[1]**2 + r[2]**2
        if R2 < 1e-30:
            continue
        target = np.sqrt(1.5 * R2)
        sums = [abs(r[0]+r[1]+r[2]), abs(r[0]+r[1]-r[2]),
                abs(r[0]-r[1]+r[2]), abs(-r[0]+r[1]+r[2])]
        min_dev = min((s - target)**2 for s in sums)
        total += min_dev
    return total

# ============================================================
# First, let's understand the constraint better analytically
# ============================================================
print("="*80)
print("ANALYTICAL APPROACH")
print("="*80)

# Consider triplets that differ in only one row.
# T1 = (x1, x2, x3) and T2 = (y1, x2, x3).
# Both must satisfy Q=3/2 (with possibly different signs).
#
# For T1: ε1*x1 + ε2*x2 + ε3*x3 = ±√(3/2) * √(x1²+x2²+x3²)
# For T2: δ1*y1 + δ2*x2 + δ3*x3 = ±√(3/2) * √(y1²+x2²+x3²)
#
# The cross terms with x2,x3 must work for BOTH.
# This is highly constraining.

# Let me try: what if we set x_i = y_i for i=2,3 (only row 1 differs)?
# Then T1 = (x1², x2², x3²) and T2 = (y1², x2², x3²).
# We need signs for each.

# For the all-positive-signs cone:
# x1² + x2² + x3² = 4(x1*x2 + x1*x3 + x2*x3)
# y1² + x2² + x3² = 4(y1*x2 + y1*x3 + x2*x3) [or different sign pattern]
#
# Subtracting (if same sign pattern):
# x1² - y1² = 4*(x1-y1)*(x2+x3)
# (x1-y1)(x1+y1) = 4*(x1-y1)*(x2+x3)
# If x1 ≠ y1: x1+y1 = 4*(x2+x3)

# And from the original: x1² + x2² + x3² = 4(x1*x2 + x1*x3 + x2*x3)
# = 4*x1*(x2+x3) + 4*x2*x3
# With S23 = x2+x3, P23 = x2*x3:
# x1² + S23² - 2*P23 = 4*x1*S23 + 4*P23
# x1² - 4*x1*S23 + S23² - 6*P23 = 0

# From x1+y1 = 4*S23: y1 = 4*S23 - x1
# And y1 must also satisfy the cone:
# y1² + S23² - 2*P23 = 4*y1*S23 + 4*P23
# y1² - 4*y1*S23 + S23² - 6*P23 = 0
# (4*S23-x1)² - 4*(4*S23-x1)*S23 + S23² - 6*P23 = 0
# 16*S23² - 8*x1*S23 + x1² - 16*S23² + 4*x1*S23 + S23² - 6*P23 = 0
# x1² - 4*x1*S23 + S23² - 6*P23 = 0
# Same equation! So it's automatically satisfied. Good.

# So: if (x1, x2, x3) is on the all-positive cone, then so is (4*(x2+x3)-x1, x2, x3),
# as long as 4*(x2+x3)-x1 >= 0.

print("\nKey identity: if (x1,x2,x3) on cone, then (4(x2+x3)-x1, x2, x3) also on cone")
print("  (provided 4(x2+x3)-x1 >= 0)")
print()

# Verify with the circle parametrization
def circle_point(theta):
    return np.array([1,1,1])/np.sqrt(6) + np.cos(theta)*np.array([1,-1,0])/2 + \
           np.sin(theta)*np.array([1,1,-2])/(2*np.sqrt(3))

# Pick a point on the circle
theta0 = 0.3
r = circle_point(theta0)
print(f"Point on circle at θ={theta0}: r = ({r[0]:.6f}, {r[1]:.6f}, {r[2]:.6f})")
y1_new = 4*(r[1]+r[2]) - r[0]
r_new = np.array([y1_new, r[1], r[2]])
Q_check = (r_new[0]+r_new[1]+r_new[2])**2 / (r_new[0]**2+r_new[1]**2+r_new[2]**2)
print(f"Reflected point: ({r_new[0]:.6f}, {r_new[1]:.6f}, {r_new[2]:.6f})")
print(f"  Q (all positive) = {Q_check:.8f}")

# Now, the masses are a1=x1², b1=y1²=(4(x2+x3)-x1)², and a_i=b_i for i=2,3.
# For 8 triplets: the 4 triplets with c1=0 use a1=x1², the 4 with c1=1 use b1=y1².
# But all use the same (x2²,x3²) for rows 2,3.
# Each of the 8 triplets needs Q=3/2 with SOME signs.

# Triplets with c1=0: (x1², x2², x3²) - needs Q=3/2 ✓ (by construction, all +)
# Triplets with c1=1: (y1², x2², x3²) - needs Q=3/2 ✓ (by construction, all +)
# But wait: c2 and c3 can be 0 or 1, and since a_i=b_i for i=2,3,
# ALL 8 triplets reduce to either (x1², x2², x3²) or (y1², x2², x3²).
# Both satisfy Q=3/2 with all positive signs. So this IS a valid solution!

print("\n*** SOLUTION TYPE 1: One pair differs, same sign pattern ***")
print("Given (x1, x2, x3) on the Q=3/2 cone (all positive signs),")
print("set y1 = 4(x2+x3) - x1, keep y2=x2, y3=x3.")
print("The grid is:")
print("  | x1²          (4(x2+x3)-x1)² |")
print("  | x2²          x2²             |")
print("  | x3²          x3²             |")
print()

# But we can also flip in each row. And we can use different sign patterns.
# Let me check: what if T1 uses all-positive but T2 uses a different sign pattern?

# For T1 = (x1,x2,x3) on all-positive cone:
# x1² + x2² + x3² = 4(x1*x2 + x1*x3 + x2*x3)  ... (A)
#
# For T2 = (y1,x2,x3) on "flip-1" cone (signs: ε12=-1, ε13=-1, ε23=+1):
# y1² + x2² + x3² = 4(-y1*x2 - y1*x3 + x2*x3)  ... (B)
#
# (A) - (B): x1²-y1² = 4*(x1+y1)*(x2+x3) - 8*x2*x3 + 4*y1*(x2+x3) - 4*x2*x3
# Hmm, let me redo:
# (A): x1² + x2² + x3² = 4*x1*(x2+x3) + 4*x2*x3
# (B): y1² + x2² + x3² = -4*y1*(x2+x3) + 4*x2*x3
# Subtract: x1² - y1² = 4*x1*(x2+x3) + 4*y1*(x2+x3) = 4*(x1+y1)*(x2+x3)
# So: (x1-y1)(x1+y1) = 4*(x1+y1)*(x2+x3)
# If x1+y1 ≠ 0: x1-y1 = 4*(x2+x3), i.e. y1 = x1 - 4*(x2+x3)

print("Mix of sign patterns:")
print("T1: all-positive, T2: flip-1")
print("  => y1 = x1 - 4*(x2+x3)")
print()

# Check if T2 actually satisfies (B):
# y1² + x2² + x3² = -4*y1*(x2+x3) + 4*x2*x3
# LHS = (x1-4S)² + x2²+x3² where S=x2+x3
# = x1²-8*x1*S+16S² + S²-2P where P=x2*x3
# = x1²-8*x1*S+17S²-2P
# RHS = -4*(x1-4S)*S + 4P = -4*x1*S+16S² + 4P
# LHS-RHS = x1²-8*x1*S+17S²-2P + 4*x1*S - 16S² - 4P
# = x1² - 4*x1*S + S² - 6P
# This should be 0 from (A): x1²+S²-2P = 4*x1*S+4P => x1²-4*x1*S+S²-6P = 0 ✓

print("Verified: if (x1,x2,x3) on all-+ cone, (x1-4(x2+x3),x2,x3) on flip-1 cone.")
print("  Mass: y1² = (x1-4(x2+x3))²")
print()

# So we have TWO reflection formulas for row 1:
# Same-sign: y1 = 4(x2+x3) - x1
# Cross-sign: y1 = x1 - 4(x2+x3) = -(4(x2+x3)-x1)
# These give y1² the same! y1² = (4(x2+x3)-x1)² in both cases.
# The mass b1 = y1² = (4(x2+x3)-x1)² regardless of sign convention.

print("BOTH reflections give the SAME mass: b1 = (4(x2+x3)-x1)²")
print("This is because flipping the sign of y1 changes the sign pattern but not the mass.")
print()

# Similarly, we can do the same for rows 2 and 3:
# b2_sqrt = |4(x1+x3) - x2|, so b2 = (4(x1+x3)-x2)²
# b3_sqrt = |4(x1+x2) - x3|, so b3 = (4(x1+x2)-x3)²

# But we need ALL 8 triplets to work. Let's see what happens when we allow
# multiple rows to differ.

print("="*80)
print("GENERAL REFLECTION: All three rows differ")
print("="*80)

# Start with (x1, x2, x3) on the all-positive cone.
# Define y_i = |4*(x_j+x_k) - x_i| for {i,j,k} = {1,2,3}.
# Then (y_i, x_j, x_k) is also on the cone for each i.
# But what about (y1, y2, x3)? Two rows differ.

# Check: is (y1, y2, x3) on ANY of the 4 cones?
# y1 = 4(x2+x3)-x1 (using the + branch)
# y2 = 4(x1+x3)-x2

# On the all-positive cone: need y1²+y2²+x3² = 4(y1*y2+y1*x3+y2*x3)?
theta0 = 0.3
r = circle_point(theta0)
x1, x2, x3 = abs(r[0]), abs(r[1]), abs(r[2])

# Make sure on cone
Q_orig = (x1+x2+x3)**2/(x1**2+x2**2+x3**2)
print(f"Original: ({x1:.6f},{x2:.6f},{x3:.6f}), Q={Q_orig:.8f}")

y1 = 4*(x2+x3) - x1
y2 = 4*(x1+x3) - x2
y3 = 4*(x1+x2) - x3

print(f"y1 = 4(x2+x3)-x1 = {y1:.6f}")
print(f"y2 = 4(x1+x3)-x2 = {y2:.6f}")
print(f"y3 = 4(x1+x2)-x3 = {y3:.6f}")

# Check (y1, x2, x3)
for s1, s2, s3 in product([1,-1], repeat=3):
    q = Q_val(y1**2, x2**2, x3**2, s1, s2, s3)
    if abs(q-1.5) < 1e-6:
        print(f"  (y1², x2², x3²): Q=3/2 with signs ({s1},{s2},{s3})")
        break

# Check (y1, y2, x3)
found = False
for s1, s2, s3 in product([1,-1], repeat=3):
    q = Q_val(y1**2, y2**2, x3**2, s1, s2, s3)
    if abs(q-1.5) < 1e-6:
        print(f"  (y1², y2², x3²): Q=3/2 with signs ({s1},{s2},{s3})")
        found = True
        break
if not found:
    print(f"  (y1², y2², x3²): NO sign combo gives Q=3/2!")
    for s1, s2, s3 in product([1,-1], repeat=3):
        q = Q_val(y1**2, y2**2, x3**2, s1, s2, s3)
        print(f"    signs ({s1},{s2},{s3}): Q = {q:.8f}")

# Check (y1, y2, y3)
found = False
for s1, s2, s3 in product([1,-1], repeat=3):
    q = Q_val(y1**2, y2**2, y3**2, s1, s2, s3)
    if abs(q-1.5) < 1e-6:
        print(f"  (y1², y2², y3²): Q=3/2 with signs ({s1},{s2},{s3})")
        found = True
        break
if not found:
    print(f"  (y1², y2², y3²): NO sign combo gives Q=3/2!")
    for s1, s2, s3 in product([1,-1], repeat=3):
        q = Q_val(y1**2, y2**2, y3**2, s1, s2, s3)
        print(f"    signs ({s1},{s2},{s3}): Q = {q:.8f}")

# So the simple reflection y_i = 4(x_j+x_k)-x_i does NOT automatically make
# ALL 8 triplets work. Only the ones where at most one row is reflected.

print("\n" + "="*80)
print("SYSTEMATIC: When can 2 rows differ?")
print("="*80)

# Need (y1, y2, x3) on some cone. Let's find what y1, y2 must satisfy.
# Using the all-positive cone (other cones just flip signs):
# y1² + y2² + x3² = 4(y1*y2 + y1*x3 + y2*x3)
#
# We have y1 = 4(x2+x3)-x1, y2 = 4(x1+x3)-x2.
# Let S = x1+x2+x3. Then:
# y1 = 4(S-x1)-x1 = 4S - 5x1
# y2 = 4(S-x2)-x2 = 4S - 5x2
# y1+y2 = 8S - 5(x1+x2) = 8S - 5(S-x3) = 3S + 5x3
# y1*y2 = (4S-5x1)(4S-5x2) = 16S² - 20S(x1+x2) + 25x1x2
#        = 16S² - 20S(S-x3) + 25x1x2 = -4S² + 20Sx3 + 25x1x2

# Cone cond: y1²+y2²+x3² = 4(y1y2+y1x3+y2x3)
# = 4(y1y2 + (y1+y2)x3)
# = 4(y1y2 + (3S+5x3)x3)
# = 4y1y2 + 4(3S+5x3)x3

# y1²+y2² = (y1+y2)² - 2y1y2 = (3S+5x3)² - 2y1y2

# LHS = (3S+5x3)² - 2y1y2 + x3²
# RHS = 4y1y2 + 12Sx3 + 20x3²

# LHS-RHS = (3S+5x3)² - 6y1y2 + x3² - 12Sx3 - 20x3²
# = 9S²+30Sx3+25x3² - 6y1y2 + x3² - 12Sx3 - 20x3²
# = 9S² + 18Sx3 + 6x3² - 6y1y2
# = 9S² + 18Sx3 + 6x3² - 6(-4S²+20Sx3+25x1x2)
# = 9S² + 18Sx3 + 6x3² + 24S² - 120Sx3 - 150x1x2
# = 33S² - 102Sx3 + 6x3² - 150x1x2

# This must be zero for (y1,y2,x3) to be on the all-positive cone.
# Using the original cone condition: S² = 6P where P = x1x2+x1x3+x2x3
# x1x2 = P - (x1+x2)x3 = P - (S-x3)x3 = P - Sx3 + x3²
# So: 150x1x2 = 150P - 150Sx3 + 150x3²

# 33S² - 102Sx3 + 6x3² - 150P + 150Sx3 - 150x3²
# = 33S² + 48Sx3 - 144x3² - 150P
# With S² = 6P, P = S²/6:
# = 33S² + 48Sx3 - 144x3² - 25S²
# = 8S² + 48Sx3 - 144x3²
# = 8(S² + 6Sx3 - 18x3²)
# = 8(S + (3+3√3)x3)(S + (3-3√3)x3)  ... let me check discriminant
# S² + 6Sx3 - 18x3² = 0 => S = (-6 ± √(36+72))x3/2 = (-6 ± √108)x3/2
# = (-6 ± 6√3)x3/2 = (-3 ± 3√3)x3
# S = (3√3-3)x3 or S = -(3√3+3)x3

# Since S,x3 > 0: S = (3√3-3)x3 = 3(√3-1)x3
# i.e., x1+x2+x3 = 3(√3-1)*x3
# i.e., x1+x2 = (3√3-4)*x3

# So (y1,y2,x3) is on the all-positive cone ONLY if x1+x2 = (3√3-4)*x3.
# This is a very specific condition!

val_3s3_4 = 3*np.sqrt(3) - 4
print(f"Condition for 2 rows to differ (all-+ cone): x1+x2 = {val_3s3_4:.6f} * x3")
print(f"  i.e., S = x1+x2+x3 = {3*(np.sqrt(3)-1):.6f} * x3")

# But we should also check the other 3 sign patterns for (y1,y2,x3).
# The condition changes for each pattern. Let me check them all.

# For the flip-3 cone: ε12=+1, ε13=-1, ε23=-1
# y1²+y2²+x3² = 4(y1y2 - y1x3 - y2x3)
# = 4y1y2 - 4(y1+y2)x3
# LHS-RHS = (3S+5x3)²-2y1y2+x3² - 4y1y2 + 4(3S+5x3)x3
# = 9S²+30Sx3+25x3²-6y1y2+x3²+12Sx3+20x3²
# = 9S²+42Sx3+46x3²-6y1y2
# = 9S²+42Sx3+46x3² - 6(-4S²+20Sx3+25x1x2)
# = 9S²+42Sx3+46x3²+24S²-120Sx3-150x1x2
# = 33S²-78Sx3+46x3²-150x1x2
# With x1x2 = P-Sx3+x3², P=S²/6:
# = 33S²-78Sx3+46x3²-150(S²/6-Sx3+x3²)
# = 33S²-78Sx3+46x3²-25S²+150Sx3-150x3²
# = 8S²+72Sx3-104x3²
# = 8(S²+9Sx3-13x3²)
# Disc: 81+52=133, √133≈11.53
# S = (-9±√133)x3/2
# S = (-9+√133)/2 * x3 ≈ 1.266*x3
# So x1+x2 ≈ 0.266*x3

val2 = (-9+np.sqrt(133))/2
print(f"\nFor flip-3 cone: x1+x2 = {val2-1:.6f} * x3, S = {val2:.6f} * x3")

# Hmm, (3√3-4) ≈ 1.196, so x1+x2 = 1.196*x3.
# And ((-9+√133)/2 - 1) ≈ 0.266, so x1+x2 = 0.266*x3.
# These are specific conditions.

# For flip-1: ε12=-1, ε13=-1, ε23=+1
# y1²+y2²+x3² = 4(-y1y2 - y1x3 + y2x3) = -4y1y2 + 4(-y1+y2)x3
# ... this gets complicated. Let me just do it numerically.

# Actually, wait. The sign pattern for (y1,y2,x3) doesn't have to be the same
# as for (x1,x2,x3). It can be any of the 4 valid patterns.
# And actually, we also need to check (x1,y2,x3), (x1,x2,y3), (y1,x2,y3), etc.

# Let me take a step back and think about this more carefully.
# The grid has 8 triplets. If all satisfy Q=3/2, then for any pair of triplets
# that differ in exactly one row, the two triplets are "compatible".

# Let me parameterize on the circle and look for solutions numerically.

print("\n" + "="*80)
print("DIRECT NUMERICAL SEARCH (small scale)")
print("="*80)

best_solutions = []

# Strategy: start from a point on the circle, define y_i by reflection,
# and see if the mixed triplets also satisfy Q=3/2.

# From the circle: u(θ) on unit sphere satisfies (u1+u2+u3)^2 = 3/2.
# Scale by R: (R*u1, R*u2, R*u3) are square roots. Masses = R^2 * u_i^2.
# The reflection for row i: y_i = 4*(x_j+x_k) - x_i where x_i = R*|u_i|, sign free.

for theta_deg in range(0, 360, 5):
    theta = np.radians(theta_deg)
    u = circle_point(theta)
    x1, x2, x3 = abs(u[0]), abs(u[1]), abs(u[2])

    # Check which points have x_i > 0 (non-zero mass)
    if min(x1,x2,x3) < 1e-10:
        continue

    # Reflections
    y1 = abs(4*(x2+x3) - x1)
    y2 = abs(4*(x1+x3) - x2)
    y3 = abs(4*(x1+x2) - x3)

    a = [x1**2, x2**2, x3**2]
    b = [y1**2, y2**2, y3**2]

    ok, results = check_grid_full(a, b, verbose=False)
    if ok:
        print(f"\n*** SOLUTION at θ={theta_deg}° ***")
        print(f"  a = ({a[0]:.8f}, {a[1]:.8f}, {a[2]:.8f})")
        print(f"  b = ({b[0]:.8f}, {b[1]:.8f}, {b[2]:.8f})")
        check_grid_full(a, b, verbose=True)
        best_solutions.append((theta_deg, a, b))

    # Also try reflecting only 1 row
    for row in [0, 1, 2]:
        b_test = list(a)  # copy
        if row == 0: b_test[0] = y1**2
        elif row == 1: b_test[1] = y2**2
        else: b_test[2] = y3**2
        ok, _ = check_grid_full(a, b_test, verbose=False)
        if ok and a != b_test:
            print(f"\n*** 1-ROW SOLUTION at θ={theta_deg}°, row {row+1} ***")
            print(f"  a = ({a[0]:.8f}, {a[1]:.8f}, {a[2]:.8f})")
            print(f"  b = ({b_test[0]:.8f}, {b_test[1]:.8f}, {b_test[2]:.8f})")
            check_grid_full(a, b_test, verbose=True)

    # Try reflecting 2 rows
    for r1, r2 in [(0,1), (0,2), (1,2)]:
        b_test = list(a)
        if r1==0: b_test[0] = y1**2
        elif r1==1: b_test[1] = y2**2
        else: b_test[2] = y3**2
        if r2==0: b_test[0] = y1**2
        elif r2==1: b_test[1] = y2**2
        else: b_test[2] = y3**2
        ok, _ = check_grid_full(a, b_test, verbose=False)
        if ok and a != b_test:
            print(f"\n*** 2-ROW SOLUTION at θ={theta_deg}°, rows {r1+1},{r2+1} ***")
            print(f"  a = ({a[0]:.8f}, {a[1]:.8f}, {a[2]:.8f})")
            print(f"  b = ({b_test[0]:.8f}, {b_test[1]:.8f}, {b_test[2]:.8f})")
            check_grid_full(a, b_test, verbose=True)

if not best_solutions:
    print("\nNo full 3-row solutions found with simple reflection.")

print("\n\nDone.")
