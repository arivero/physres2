#!/usr/bin/env python3
"""
Fit quark mass square roots to linear combinations of SU(5) Cartan charges.

Cartan generators of SU(5) in the fundamental:
  H_1 = diag(1,-1,0,0,0)
  H_2 = diag(0,1,-1,0,0)
  H_3 = diag(0,0,1,-1,0)
  H_4 = diag(0,0,0,1,-1)

Quark k (k=1..5) has charge Q_i(k) = (H_i)_{kk}.

Target: z_k = sum_i alpha_i * Q_i(k) + gamma
"""

import numpy as np
from itertools import product

# Quark masses in MeV (PDG central values)
m_quarks = np.array([2.2, 4.7, 93.4, 1270.0, 4180.0])
labels = ['u', 'd', 's', 'c', 'b']

# Square roots
z_pos = np.sqrt(m_quarks)
print("Square roots of quark masses (MeV^{1/2}):")
for i, (l, z) in enumerate(zip(labels, z_pos)):
    print(f"  z_{l} = {z:.6f}")
print()

# Cartan charges: H_i diagonal entries for fundamental 5
# H_1 = diag(1,-1,0,0,0)
# H_2 = diag(0,1,-1,0,0)
# H_3 = diag(0,0,1,-1,0)
# H_4 = diag(0,0,0,1,-1)
Q = np.array([
    [1, -1, 0, 0, 0],   # H_1
    [0,  1,-1, 0, 0],   # H_2
    [0,  0, 1,-1, 0],   # H_3
    [0,  0, 0, 1,-1],   # H_4
])

# ============================================================
# PART 1: Exact solution with 4 Cartan generators + constant
# ============================================================
print("=" * 70)
print("PART 1: Exact solution (4 Cartan + constant), all-positive z_k")
print("=" * 70)

# System: for each k=1..5:
#   z_k = alpha_1*Q_1(k) + alpha_2*Q_2(k) + alpha_3*Q_3(k) + alpha_4*Q_4(k) + gamma
# Matrix form: A * [alpha_1, alpha_2, alpha_3, alpha_4, gamma]^T = z

A_full = np.zeros((5, 5))
for k in range(5):
    for i in range(4):
        A_full[k, i] = Q[i, k]
    A_full[k, 4] = 1.0  # constant term

def solve_and_report(z_vals, sign_label=""):
    """Solve A_full * params = z_vals and report."""
    params = np.linalg.solve(A_full, z_vals)
    alpha = params[:4]
    gamma = params[4]

    # Verify
    z_fit = A_full @ params
    residual = np.max(np.abs(z_fit - z_vals))

    print(f"  alpha_1 = {alpha[0]:+.6f}")
    print(f"  alpha_2 = {alpha[1]:+.6f}")
    print(f"  alpha_3 = {alpha[2]:+.6f}")
    print(f"  alpha_4 = {alpha[3]:+.6f}")
    print(f"  gamma   = {gamma:+.6f}")
    print(f"  Max residual: {residual:.2e}")

    # Check ratios
    print(f"\n  Ratios:")
    for i in range(4):
        for j in range(i+1, 4):
            if abs(alpha[j]) > 1e-10:
                print(f"    alpha_{i+1}/alpha_{j+1} = {alpha[i]/alpha[j]:+.6f}")
    if abs(gamma) > 1e-10:
        for i in range(4):
            print(f"    alpha_{i+1}/gamma = {alpha[i]/gamma:+.6f}")

    # Check differences (since Cartan charges are differences of adjacent entries)
    # The actual z differences:
    print(f"\n  Successive differences of z_k:")
    for k in range(4):
        print(f"    z_{k+2} - z_{k+1} = {z_vals[k+1] - z_vals[k]:+.6f}")

    return params

print("\nAll-positive assignment:")
solve_and_report(z_pos)

# ============================================================
# PART 2: All 32 sign combinations
# ============================================================
print("\n" + "=" * 70)
print("PART 2: All 32 sign combinations")
print("=" * 70)

results = []

for signs in product([-1, 1], repeat=5):
    s = np.array(signs)
    z_signed = s * z_pos
    params = np.linalg.solve(A_full, z_signed)
    alpha = params[:4]
    gamma = params[4]

    # Compute a "niceness" score: how close are alphas to simple ratios?
    # Try: sum of |alpha_i| (smaller = nicer), or check rationality

    # Score 1: L2 norm of alphas (prefer small)
    norm = np.linalg.norm(alpha)

    # Score 2: check if alphas are approximately rational with small denominators
    # Try denominators 1..20
    def rationality_score(x, max_denom=20):
        best = float('inf')
        for d in range(1, max_denom+1):
            n = round(x * d)
            err = abs(x * d - n) / d
            if err < best:
                best = err
        return best

    rat_score = sum(rationality_score(a) for a in alpha) + rationality_score(gamma)

    sign_str = ''.join(['+' if s > 0 else '-' for s in signs])
    results.append((sign_str, params.copy(), norm, rat_score, s))

# Sort by rationality score
results.sort(key=lambda x: x[3])

print("\nTop 10 by rationality score (lower = more rational):")
print(f"{'Signs':>8s}  {'alpha_1':>10s} {'alpha_2':>10s} {'alpha_3':>10s} {'alpha_4':>10s} {'gamma':>10s}  {'|alpha|':>8s} {'rat_score':>10s}")
for sign_str, params, norm, rat_score, s in results[:10]:
    print(f"{sign_str:>8s}  {params[0]:+10.4f} {params[1]:+10.4f} {params[2]:+10.4f} {params[3]:+10.4f} {params[4]:+10.4f}  {norm:8.4f} {rat_score:10.6f}")

# Sort by norm
results.sort(key=lambda x: x[2])
print("\nTop 10 by smallest |alpha| norm:")
print(f"{'Signs':>8s}  {'alpha_1':>10s} {'alpha_2':>10s} {'alpha_3':>10s} {'alpha_4':>10s} {'gamma':>10s}  {'|alpha|':>8s}")
for sign_str, params, norm, rat_score, s in results[:10]:
    print(f"{sign_str:>8s}  {params[0]:+10.4f} {params[1]:+10.4f} {params[2]:+10.4f} {params[3]:+10.4f} {params[4]:+10.4f}  {norm:8.4f}")

# ============================================================
# PART 3: Detailed analysis of interesting sign combos
# ============================================================
print("\n" + "=" * 70)
print("PART 3: Detailed analysis of best sign combinations")
print("=" * 70)

# The Koide-relevant sign combo: (-s, c, b) or more generally signs related to sBootstrap
# Let's check specifically: all positive, alternating, and the (-u, d, s, c, b) etc.

interesting = [
    ('+', '+', '+', '+', '+'),
    ('-', '+', '+', '+', '+'),
    ('+', '-', '+', '+', '+'),
    ('+', '+', '-', '+', '+'),
    ('-', '-', '+', '+', '+'),
    ('-', '+', '-', '+', '+'),
    ('+', '-', '-', '+', '+'),
    ('-', '-', '-', '+', '+'),
    ('+', '+', '+', '-', '+'),
    ('+', '+', '-', '-', '+'),
]

for combo in interesting:
    s = np.array([1 if c == '+' else -1 for c in combo])
    z_signed = s * z_pos
    params = np.linalg.solve(A_full, z_signed)
    alpha = params[:4]
    gamma = params[4]

    sign_str = ''.join(combo)
    assigned = ', '.join([f"{c}{l}" for c, l in zip(combo, labels)])
    print(f"\nSigns: {sign_str}  ({assigned})")
    print(f"  alpha = [{alpha[0]:+.6f}, {alpha[1]:+.6f}, {alpha[2]:+.6f}, {alpha[3]:+.6f}]")
    print(f"  gamma = {gamma:+.6f}")

    # Express in terms of z values explicitly
    # Note: the Cartan structure means:
    # alpha_1 = z_1 - z_2 (from H_1 acting on positions 1,2)
    # Actually let's derive the exact relations
    # A_full * params = z  =>  params = A_full^{-1} * z
    # Let's compute A_full^{-1} once

A_inv = np.linalg.inv(A_full)
print("\n\nExplicit solution formulas (A_inv matrix):")
print("  params = A_inv @ z, where params = [alpha_1, alpha_2, alpha_3, alpha_4, gamma]")
print(f"\n  A_inv =")
for i in range(5):
    row = ' '.join([f"{A_inv[i,j]:+.4f}" for j in range(5)])
    print(f"    [{row}]")

print("\nSo:")
param_names = ['alpha_1', 'alpha_2', 'alpha_3', 'alpha_4', 'gamma']
for i in range(5):
    terms = []
    for j in range(5):
        if abs(A_inv[i,j]) > 1e-10:
            coeff = A_inv[i,j]
            if abs(coeff - round(coeff)) < 1e-10:
                coeff = int(round(coeff))
            terms.append(f"{coeff:+} * z_{labels[j]}")
    print(f"  {param_names[i]} = {' '.join(terms)}")

# ============================================================
# PART 4: Least-squares with 2 Cartan generators + constant
# ============================================================
print("\n" + "=" * 70)
print("PART 4: Least-squares with 2 Cartan generators + constant")
print("=" * 70)

# Try all pairs of Cartan generators
from itertools import combinations

for (i, j) in combinations(range(4), 2):
    A_red = np.zeros((5, 3))
    for k in range(5):
        A_red[k, 0] = Q[i, k]
        A_red[k, 1] = Q[j, k]
        A_red[k, 2] = 1.0

    # Least squares for all-positive z
    params_ls, residuals, rank, sv = np.linalg.lstsq(A_red, z_pos, rcond=None)
    z_fit = A_red @ params_ls
    rms_residual = np.sqrt(np.mean((z_fit - z_pos)**2))
    max_residual = np.max(np.abs(z_fit - z_pos))
    rel_errors = np.abs(z_fit - z_pos) / z_pos * 100

    print(f"\nH_{i+1}, H_{j+1} + constant:")
    print(f"  alpha_{i+1} = {params_ls[0]:+.4f}, alpha_{j+1} = {params_ls[1]:+.4f}, gamma = {params_ls[2]:+.4f}")
    print(f"  RMS residual: {rms_residual:.4f}")
    print(f"  Max residual: {max_residual:.4f}")
    print(f"  Fitted vs target:")
    for k in range(5):
        print(f"    z_{labels[k]}: fit={z_fit[k]:+.4f}  target={z_pos[k]:.4f}  err={rel_errors[k]:.1f}%")

# Also try all sign combos for best 2-Cartan fit
print("\n\nBest 2-Cartan fits across all sign combos:")
best_2cartan = []

for (i, j) in combinations(range(4), 2):
    A_red = np.zeros((5, 3))
    for k in range(5):
        A_red[k, 0] = Q[i, k]
        A_red[k, 1] = Q[j, k]
        A_red[k, 2] = 1.0

    for signs in product([-1, 1], repeat=5):
        s = np.array(signs)
        z_signed = s * z_pos
        params_ls, residuals, rank, sv = np.linalg.lstsq(A_red, z_signed, rcond=None)
        z_fit = A_red @ params_ls
        rms = np.sqrt(np.mean((z_fit - z_signed)**2))
        sign_str = ''.join(['+' if si > 0 else '-' for si in signs])
        best_2cartan.append((rms, i, j, sign_str, params_ls.copy(), z_fit.copy(), z_signed.copy()))

best_2cartan.sort(key=lambda x: x[0])
print(f"\n{'RMS':>8s} {'Cartans':>8s} {'Signs':>8s}  {'a_i':>10s} {'a_j':>10s} {'gamma':>10s}")
for rms, i, j, sign_str, params, z_fit, z_signed in best_2cartan[:15]:
    print(f"{rms:8.4f} H{i+1},H{j+1}    {sign_str:>8s}  {params[0]:+10.4f} {params[1]:+10.4f} {params[2]:+10.4f}")

# Show detail for the best one
rms, i, j, sign_str, params, z_fit, z_signed = best_2cartan[0]
print(f"\nBest 2-Cartan fit: H_{i+1}, H_{j+1}, signs={sign_str}")
print(f"  a_{i+1} = {params[0]:+.6f}, a_{j+1} = {params[1]:+.6f}, gamma = {params[2]:+.6f}")
for k in range(5):
    print(f"  z_{labels[k]}: fit={z_fit[k]:+.6f}  signed_target={z_signed[k]:+.6f}  diff={z_fit[k]-z_signed[k]:+.6f}")

# ============================================================
# PART 5: Pattern analysis
# ============================================================
print("\n" + "=" * 70)
print("PART 5: Pattern analysis of the exact solution")
print("=" * 70)

# For the all-positive case
params_exact = np.linalg.solve(A_full, z_pos)
alpha = params_exact[:4]
gamma = params_exact[4]

print("\nAll-positive exact solution:")
print(f"  alpha = [{', '.join([f'{a:+.6f}' for a in alpha])}]")
print(f"  gamma = {gamma:+.6f}")

print(f"\n  Successive alpha ratios:")
for i in range(3):
    print(f"    alpha_{i+2}/alpha_{i+1} = {alpha[i+1]/alpha[i]:+.6f}")

print(f"\n  alpha_i / gamma:")
for i in range(4):
    print(f"    alpha_{i+1}/gamma = {alpha[i]/gamma:+.6f}")

# The alphas are just differences of adjacent z values (from Cartan structure)
# Let's verify
print(f"\n  Note: from A_inv structure, the alphas are cumulative sums of z-differences:")
print(f"    alpha_1 = z_1 - z_2 = {z_pos[0] - z_pos[1]:+.6f}  vs alpha_1 = {alpha[0]:+.6f}")
# Actually let's just read off from A_inv
print(f"\n  Explicit alpha formulas from A_inv:")
for i in range(5):
    terms = []
    for j in range(5):
        c = A_inv[i, j]
        if abs(c) > 1e-10:
            if abs(c - round(c)) < 1e-10:
                c_str = f"{int(round(c)):+d}"
            else:
                c_str = f"{c:+.2f}"
            terms.append(f"{c_str}*z_{labels[j]}")
    print(f"    {param_names[i]} = {' '.join(terms)}")

# Check if the solution has structure related to Koide
print(f"\n  Koide-relevant quantities:")
v0 = np.sum(z_pos)
Q_koide = np.sum(m_quarks) / v0**2
print(f"    v_0 = sum(z_k) = {v0:.6f}")
print(f"    Q = sum(m_k) / v_0^2 = {Q_koide:.6f}  (2/3 = {2/3:.6f})")
print(f"    gamma = v_0/5 = {v0/5:.6f}  (vs actual gamma = {gamma:.6f})")
# gamma should be v_0/5 since sum of Cartan charges = 0 for each H_i
print(f"    gamma = mean(z_k) = {np.mean(z_pos):.6f}  ← must equal gamma (Cartan trace-free)")

# Check: do the alphas relate to any mass scale?
print(f"\n  Alpha magnitudes vs mass scales:")
for i in range(4):
    print(f"    |alpha_{i+1}| = {abs(alpha[i]):.4f} MeV^{{1/2}}, |alpha_{i+1}|^2 = {alpha[i]**2:.4f} MeV")

# The key structural insight: with trace-free Cartan generators,
# gamma = mean(z_k) is forced, and the alphas encode the DEVIATIONS from the mean
print(f"\n  Deviations from mean:")
z_dev = z_pos - gamma
for k in range(5):
    print(f"    z_{labels[k]} - gamma = {z_dev[k]:+.6f}")
print(f"    Sum of deviations = {np.sum(z_dev):.2e} (must be 0)")

# ============================================================
# PART 6: Alternative basis - Dynkin basis
# ============================================================
print("\n" + "=" * 70)
print("PART 6: Check with Gell-Mann-like diagonal generators")
print("=" * 70)

# Instead of simple roots, use the standard diagonal Gell-Mann matrices for SU(5)
# These are: lambda_3, lambda_8, lambda_15, lambda_24-like
# Normalized as Tr(T_a T_b) = 1/2 delta_ab

T3 = np.diag([1, -1, 0, 0, 0]) / 2
T8 = np.diag([1, 1, -2, 0, 0]) / (2*np.sqrt(3))
T15 = np.diag([1, 1, 1, -3, 0]) / (2*np.sqrt(6))
T24 = np.diag([1, 1, 1, 1, -4]) / (2*np.sqrt(10))

GM_charges = np.array([
    [T3[k,k] for k in range(5)],
    [T8[k,k] for k in range(5)],
    [T15[k,k] for k in range(5)],
    [T24[k,k] for k in range(5)],
])

A_GM = np.zeros((5, 5))
for k in range(5):
    for i in range(4):
        A_GM[k, i] = GM_charges[i, k]
    A_GM[k, 4] = 1.0

params_GM = np.linalg.solve(A_GM, z_pos)
beta = params_GM[:4]
gamma_GM = params_GM[4]

print("Gell-Mann basis (properly normalized):")
print(f"  beta_3  = {beta[0]:+.6f}")
print(f"  beta_8  = {beta[1]:+.6f}")
print(f"  beta_15 = {beta[2]:+.6f}")
print(f"  beta_24 = {beta[3]:+.6f}")
print(f"  gamma   = {gamma_GM:+.6f}")

print(f"\n  Ratios:")
for i in range(4):
    for j in range(i+1, 4):
        if abs(beta[j]) > 1e-10:
            print(f"    beta_{[3,8,15,24][i]}/beta_{[3,8,15,24][j]} = {beta[i]/beta[j]:+.6f}")

print(f"\n  beta magnitudes: {np.abs(beta)}")
print(f"  beta^2: {beta**2}")
print(f"  sum(beta^2) = {np.sum(beta**2):.6f}")

# Verify
z_check = A_GM @ params_GM
print(f"\n  Verification: max |z_fit - z_target| = {np.max(np.abs(z_check - z_pos)):.2e}")

# ============================================================
# PART 7: Weight space interpretation
# ============================================================
print("\n" + "=" * 70)
print("PART 7: Weight space — the 5 quarks as weight vectors")
print("=" * 70)

# Weights of fundamental 5 of SU(5) in Cartan basis (H_1,...,H_4)
weights = Q.T  # 5 x 4 matrix, each row is a weight vector
print("Weight vectors of the fundamental 5:")
for k in range(5):
    print(f"  w_{labels[k]} = {weights[k]}")

# The alpha vector in weight space
print(f"\n  alpha vector = {alpha}")
print(f"  |alpha| = {np.linalg.norm(alpha):.6f}")

# Dot products with simple roots
simple_roots = np.array([
    [1, -1, 0, 0],
    [0, 1, -1, 0],
    [0, 0, 1, -1],
])
# Actually the simple roots of SU(5) in the Cartan basis ARE the rows of Q (shifted)
# The weights of the fundamental are the rows of Q.T

# The alpha vector maps each weight to a z-deviation
# z_k - gamma = alpha . w_k
print(f"\n  Check: alpha . w_k = z_k - gamma?")
for k in range(5):
    dot = np.dot(alpha, weights[k])
    print(f"    alpha . w_{labels[k]} = {dot:+.6f}  vs  z_{labels[k]} - gamma = {z_pos[k] - gamma:+.6f}")

print("\n\nDone.")
