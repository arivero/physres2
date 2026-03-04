#!/usr/bin/env python3
"""
Solve systems of Diophantine equations over positive integers (N, r, s).

System A:  rs = 2N,  r(r+1)/2 = 2N,  r^2 + s^2 - 1 = 4N
System B:  rs = 2N,  r(r+1)/2 = 2N   (drop eq 3)
System C:  rs = 2N,  r(r+1)/2 = 2N,  r^2 + s^2 - 1 = 2N
"""

print("=" * 70)
print("SYSTEM A: rs = 2N, r(r+1)/2 = 2N, r^2 + s^2 - 1 = 4N")
print("=" * 70)

# --- METHOD (i): Algebraic elimination ---
print("\n--- Method (i): Algebraic elimination ---\n")

print("From (1): 2N = rs")
print("From (2): 2N = r(r+1)/2")
print()
print("Step 1: Eliminate N using (1) and (2):")
print("  rs = r(r+1)/2")
print("  Since r > 0, divide both sides by r:")
print("  s = (r+1)/2")
print()
print("For s to be a positive integer, r must be odd. Write r = 2k - 1 for k >= 1.")
print("Then s = k.")
print()
print("Step 2: Substitute s = (r+1)/2 into equation (3): r^2 + s^2 - 1 = 4N")
print("  From (1): 4N = 2rs = 2r * (r+1)/2 = r(r+1)")
print("  So equation (3) becomes:")
print("    r^2 + ((r+1)/2)^2 - 1 = r(r+1)")
print("    r^2 + (r+1)^2/4 - 1 = r^2 + r")
print("    (r+1)^2/4 - 1 = r")
print("    (r+1)^2 - 4 = 4r")
print("    r^2 + 2r + 1 - 4 = 4r")
print("    r^2 - 2r - 3 = 0")
print("    (r - 3)(r + 1) = 0")
print()
print("  So r = 3 or r = -1. Since r must be a positive integer, r = 3.")
print("  Then s = (3+1)/2 = 2, and 2N = rs = 6, so N = 3.")
print()
print("  ALGEBRAIC SOLUTION: (N, r, s) = (3, 3, 2)")
print()

# Verify
r_alg, s_alg = 3, 2
N_alg = (r_alg * s_alg) // 2
assert r_alg * s_alg == 2 * N_alg, "Eq (1) fails"
assert r_alg * (r_alg + 1) // 2 == 2 * N_alg, "Eq (2) fails"
assert r_alg**2 + s_alg**2 - 1 == 4 * N_alg, "Eq (3) fails"
print(f"  Verification: rs = {r_alg*s_alg} = 2*{N_alg} = {2*N_alg}  [OK]")
print(f"  Verification: r(r+1)/2 = {r_alg*(r_alg+1)//2} = 2*{N_alg} = {2*N_alg}  [OK]")
print(f"  Verification: r^2+s^2-1 = {r_alg**2+s_alg**2-1} = 4*{N_alg} = {4*N_alg}  [OK]")

# --- METHOD (ii): Brute-force enumeration ---
print("\n--- Method (ii): Brute-force enumeration (1 <= r, s <= 100) ---\n")

solutions_A_brute = []
for r in range(1, 101):
    for s in range(1, 101):
        val_rs = r * s
        if val_rs % 2 != 0:
            continue
        N = val_rs // 2
        if N <= 0:
            continue
        # Check eq (2): r(r+1)/2 = 2N
        if r * (r + 1) // 2 != 2 * N:
            continue
        # Check eq (3): r^2 + s^2 - 1 = 4N
        if r**2 + s**2 - 1 != 4 * N:
            continue
        solutions_A_brute.append((N, r, s))

print(f"  Found {len(solutions_A_brute)} solution(s):")
for sol in solutions_A_brute:
    print(f"    (N, r, s) = {sol}")

# Compare
print("\n--- Comparison of methods ---")
solutions_A_algebraic = [(3, 3, 2)]
if solutions_A_brute == solutions_A_algebraic:
    print("  AGREEMENT: Both methods find exactly (N, r, s) = (3, 3, 2).")
else:
    print("  DISAGREEMENT!")
    print(f"    Algebraic: {solutions_A_algebraic}")
    print(f"    Brute-force: {solutions_A_brute}")

# List all solutions with N <= 1000
print("\n--- System A solutions with N <= 1000 ---")
solutions_A_1000 = []
for r in range(1, 2001):
    # From eq (2): 2N = r(r+1)/2, so N = r(r+1)/4. Need r(r+1) divisible by 4.
    if (r * (r + 1)) % 4 != 0:
        continue
    N = r * (r + 1) // 4
    if N > 1000:
        break
    # From eq (1): s = 2N/r
    if (2 * N) % r != 0:
        continue
    s = (2 * N) // r
    if s <= 0:
        continue
    # Check eq (3)
    if r**2 + s**2 - 1 == 4 * N:
        solutions_A_1000.append((N, r, s))

print(f"  Found {len(solutions_A_1000)} solution(s) with N <= 1000:")
for sol in solutions_A_1000:
    print(f"    (N, r, s) = {sol}")

print("\n  UNIQUENESS PROOF:")
print("  The algebraic reduction shows r must satisfy r^2 - 2r - 3 = 0,")
print("  whose only positive root is r = 3. Since r determines s = (r+1)/2 = 2")
print("  and N = rs/2 = 3, the solution (N, r, s) = (3, 3, 2) is UNIQUE")
print("  over all positive integers (not just N <= 1000).")

print()
print("=" * 70)
print("SYSTEM B: rs = 2N, r(r+1)/2 = 2N  (only equations 1 and 2)")
print("=" * 70)

print("\n--- Step 1: Show 2N must be a triangular number ---\n")
print("  Equation (2) says 2N = r(r+1)/2 = T_r, the r-th triangular number.")
print("  So 2N is necessarily a triangular number.")
print("  Combined with (1): rs = r(r+1)/2, giving s = (r+1)/2 (for r > 0),")
print("  so r must be odd for s to be a positive integer.")

print("\n--- Step 2: All solutions with N <= 200 ---\n")

solutions_B = []
for r in range(1, 2000):
    # r must be odd
    if r % 2 == 0:
        continue
    s = (r + 1) // 2
    N_times_2 = r * s  # = r(r+1)/2
    if N_times_2 % 2 != 0:
        # N must be a positive integer
        # Actually N = rs/2. Check if rs is even.
        # rs = r * (r+1)/2. r is odd, (r+1)/2 = s. r*s: r odd, s = (r+1)/2.
        # If r = 1: s=1, rs=1, N=1/2 -> not integer!
        # Need rs even, i.e., r*s even. r is odd, so s must be even.
        # s = (r+1)/2 is even when r+1 divisible by 4, i.e., r = 3 mod 4.
        continue
    N = N_times_2 // 2
    if N <= 0 or N > 200:
        continue
    solutions_B.append((N, r, s))

# Also handle the case more carefully: N = r(r+1)/4
# Need r(r+1) divisible by 4. r odd means r+1 even.
# r+1 even: if r = 1 mod 4, then r+1 = 2 mod 4, so r(r+1)/4 = r*(r+1)/4.
#   r+1 = 2 mod 4 means (r+1)/2 is odd. r(r+1)/4 = r*((r+1)/2)/2.
#   r odd, (r+1)/2 odd, so r*(r+1)/2 is odd, not divisible by 2. N not integer.
# If r = 3 mod 4, then r+1 = 0 mod 4, so (r+1)/4 is integer, and N = r(r+1)/4 is integer.
# So r = 3 mod 4.

print(f"  {'N':>6s}  {'r':>6s}  {'s':>6s}  {'2N':>8s}  {'T_r':>8s}  {'r^2+s^2-1':>10s}  {'4N':>8s}  {'2N':>8s}  {'eq3(4N)?':>9s}  {'eqC(2N)?':>9s}")
print(f"  {'---':>6s}  {'---':>6s}  {'---':>6s}  {'---':>8s}  {'---':>8s}  {'---':>10s}  {'---':>8s}  {'---':>8s}  {'---':>9s}  {'---':>9s}")

for (N, r, s) in solutions_B:
    T_r = r * (r + 1) // 2
    lhs3 = r**2 + s**2 - 1
    fourN = 4 * N
    twoN = 2 * N
    eq3_sat = "YES" if lhs3 == fourN else "no"
    eqC_sat = "YES" if lhs3 == twoN else "no"
    print(f"  {N:>6d}  {r:>6d}  {s:>6d}  {twoN:>8d}  {T_r:>8d}  {lhs3:>10d}  {fourN:>8d}  {twoN:>8d}  {eq3_sat:>9s}  {eqC_sat:>9s}")

print()
print("=" * 70)
print("SYSTEM C: rs = 2N, r(r+1)/2 = 2N, r^2 + s^2 - 1 = 2N")
print("=" * 70)

print("\n--- Algebraic approach ---\n")
print("  From (1) and (2): s = (r+1)/2 (same as before), r odd.")
print("  Equation (C): r^2 + s^2 - 1 = 2N = rs")
print("  So: r^2 + ((r+1)/2)^2 - 1 = r*(r+1)/2")
print("    r^2 + (r+1)^2/4 - 1 = r(r+1)/2")
print("    4r^2 + (r+1)^2 - 4 = 2r(r+1)       [multiply by 4]")
print("    4r^2 + r^2 + 2r + 1 - 4 = 2r^2 + 2r")
print("    5r^2 + 2r - 3 = 2r^2 + 2r")
print("    3r^2 - 3 = 0")
print("    r^2 = 1")
print("    r = 1  (positive root)")
print()
print("  If r = 1: s = (1+1)/2 = 1, N = rs/2 = 1/2. NOT an integer!")
print("  So System C has NO positive integer solutions.")

# Brute-force verify
print("\n--- Brute-force verification (1 <= r, s <= 200) ---")
solutions_C = []
for r in range(1, 201):
    for s in range(1, 201):
        val_rs = r * s
        if val_rs % 2 != 0:
            continue
        N = val_rs // 2
        if N <= 0:
            continue
        if r * (r + 1) // 2 != 2 * N:
            continue
        if r**2 + s**2 - 1 != 2 * N:
            continue
        solutions_C.append((N, r, s))

if len(solutions_C) == 0:
    print("  Confirmed: NO solutions found. System C has no positive integer solutions.")
else:
    print(f"  Found {len(solutions_C)} solution(s): {solutions_C}")

print()
print("=" * 70)
print("BONUS: n_f = r + s for each System B solution, flag 2*n_f < 33")
print("=" * 70)
print()

# Recompute System B solutions with N <= 200
print(f"  {'N':>6s}  {'r':>6s}  {'s':>6s}  {'n_f=r+s':>8s}  {'2*n_f':>6s}  {'2*n_f<33?':>10s}")
print(f"  {'---':>6s}  {'---':>6s}  {'---':>6s}  {'---':>8s}  {'---':>6s}  {'---':>10s}")

for (N, r, s) in solutions_B:
    nf = r + s
    two_nf = 2 * nf
    flag = "YES" if two_nf < 33 else "no"
    print(f"  {N:>6d}  {r:>6d}  {s:>6d}  {nf:>8d}  {two_nf:>6d}  {flag:>10s}")

print()
print("=" * 70)
print("EXTENDED: System B solutions with N <= 1000")
print("=" * 70)
print()

solutions_B_1000 = []
for r in range(1, 10000):
    if r % 2 == 0:
        continue
    s = (r + 1) // 2
    rs = r * s
    if rs % 2 != 0:
        continue
    N = rs // 2
    if N <= 0 or N > 1000:
        continue
    solutions_B_1000.append((N, r, s))

print(f"  {'N':>6s}  {'r':>6s}  {'s':>6s}  {'n_f=r+s':>8s}  {'2*n_f':>6s}  {'2*n_f<33?':>10s}  {'r^2+s^2-1':>10s}  {'4N':>8s}  {'eq3?':>5s}")
print(f"  {'---':>6s}  {'---':>6s}  {'---':>6s}  {'---':>8s}  {'---':>6s}  {'---':>10s}  {'---':>10s}  {'---':>8s}  {'---':>5s}")

for (N, r, s) in solutions_B_1000:
    nf = r + s
    two_nf = 2 * nf
    flag = "YES" if two_nf < 33 else "no"
    lhs3 = r**2 + s**2 - 1
    fourN = 4 * N
    eq3 = "YES" if lhs3 == fourN else "no"
    print(f"  {N:>6d}  {r:>6d}  {s:>6d}  {nf:>8d}  {two_nf:>6d}  {flag:>10s}  {lhs3:>10d}  {fourN:>8d}  {eq3:>5s}")

print()
print("=" * 70)
print("SUMMARY")
print("=" * 70)
print()
print("System A: (N, r, s) = (3, 3, 2) is the UNIQUE solution.")
print("  - Algebraic and brute-force methods AGREE.")
print("  - Uniqueness is proven: the quadratic r^2 - 2r - 3 = 0 has")
print("    r = 3 as its only positive root.")
print()
print(f"System B: {len(solutions_B)} solutions with N <= 200.")
print("  - Only the solution (3, 3, 2) also satisfies eq (3): r^2+s^2-1 = 4N.")
print("  - None of the other solutions satisfy either r^2+s^2-1 = 4N or = 2N.")
print()
print("System C: NO positive integer solutions exist.")
print("  - The algebraic reduction gives r = 1, s = 1, N = 1/2 (not integer).")
