#!/usr/bin/env python3
"""Round 2-C: Statistical significance computations"""
import numpy as np
from itertools import combinations
import sys

# ============================================================
# Task 1: Look-elsewhere correction for Koide scan
# ============================================================
print("=" * 60)
print("TASK 1: Look-elsewhere corrected p-value for Koide scan")
print("=" * 60)

def Q_best(masses_9):
    """For 9 masses, find best |Q - 3/2| over all 84 triplets x 8 sign combos."""
    best = np.inf
    sqrts = np.sqrt(masses_9)
    for idx in combinations(range(9), 3):
        m = masses_9[list(idx)]
        sq = sqrts[list(idx)]
        denom = m.sum()
        for s1 in [1, -1]:
            for s2 in [1, -1]:
                for s3 in [1, -1]:
                    num = (s1*sq[0] + s2*sq[1] + s3*sq[2])**2
                    Q = num / denom
                    diff = abs(Q - 1.5)
                    if diff < best:
                        best = diff
    return best

# Actual SM value
sm_masses = np.array([0.51099895, 105.6583755, 1776.86, 2.16, 4.67, 93.4, 1270, 4180, 172690])
sm_best = Q_best(sm_masses)
print(f"SM best |Q - 3/2| = {sm_best:.6e}")

N_MC = 10000
np.random.seed(42)

# Log-uniform in [0.1, 200000]
count_wide = 0
for i in range(N_MC):
    if i % 1000 == 0:
        print(f"  Wide range: {i}/{N_MC}...", file=sys.stderr)
    log_masses = np.random.uniform(np.log(0.1), np.log(200000), 9)
    masses = np.exp(log_masses)
    if Q_best(masses) <= sm_best:
        count_wide += 1

p_wide = count_wide / N_MC
print(f"\nWide range [0.1, 200000] MeV:")
print(f"  p-value = {count_wide}/{N_MC} = {p_wide:.6f}")
print(f"  Equivalent sigma = {-1 if p_wide == 0 else np.abs(np.sqrt(2) * np.erfc(2*p_wide)):.2f}" if p_wide > 0 else "  p = 0 (> {:.1f} sigma)".format(np.sqrt(2) * 3.7))

# Log-uniform in [0.5, 180000] (tighter)
count_tight = 0
for i in range(N_MC):
    if i % 1000 == 0:
        print(f"  Tight range: {i}/{N_MC}...", file=sys.stderr)
    log_masses = np.random.uniform(np.log(0.5), np.log(180000), 9)
    masses = np.exp(log_masses)
    if Q_best(masses) <= sm_best:
        count_tight += 1

p_tight = count_tight / N_MC
print(f"\nTight range [0.5, 180000] MeV:")
print(f"  p-value = {count_tight}/{N_MC} = {p_tight:.6f}")

# Convert to sigma
from scipy.stats import norm
for label, p in [("Wide", p_wide), ("Tight", p_tight)]:
    if p > 0:
        sigma = norm.ppf(1 - p)
        print(f"  {label}: p = {p:.4e}, sigma = {sigma:.2f}")
    else:
        print(f"  {label}: p < 1/{N_MC} = {1/N_MC:.1e}, sigma > {norm.ppf(1-1/N_MC):.1f}")

# ============================================================
# Task 2: Statistical significance of δ₀ mod 2π/3 ≈ 2/9
# ============================================================
print("\n" + "=" * 60)
print("TASK 2: Is δ₀ mod 2π/3 ≈ 2/9 statistically significant?")
print("=" * 60)

# Nice fractions p/q with q <= 20
nice_fracs = set()
period = 2 * np.pi / 3
for q in range(1, 21):
    for p in range(0, q):
        val = p / q
        if val < period:  # mod period, so val in [0, period)
            nice_fracs.add((p, q, val))

# Actually δ₀ mod 2π/3 gives a value in [0, 2π/3). We compare to p/q as a raw number.
# The claim is δ₀ mod 2π/3 = 0.222230 ≈ 2/9 = 0.222222
# residual = 7.4e-6
residual_observed = abs(0.222230 - 2/9)
print(f"Observed residual: |0.222230 - 2/9| = {residual_observed:.2e}")

# For random δ₀ uniform in [0, 2π/3), compute δ₀ mod 2π/3 = δ₀ itself
# Check if it lands within residual_observed of ANY nice fraction p/q with q ≤ 20
nice_values = sorted(set(p/q for p, q, _ in nice_fracs if p/q < period))
nice_arr = np.array(nice_values)

N_MC2 = 100000
np.random.seed(123)
count_near = 0
for i in range(N_MC2):
    delta = np.random.uniform(0, period)
    # min distance to any nice fraction
    dists = np.abs(delta - nice_arr)
    if np.min(dists) <= residual_observed:
        count_near += 1

p_nice = count_near / N_MC2
print(f"\nFraction of random δ₀ within {residual_observed:.2e} of any p/q (q≤20): {count_near}/{N_MC2} = {p_nice:.6f}")
print(f"Number of nice fractions in [0, {period:.4f}): {len(nice_values)}")
print(f"Expected (analytic): 2 × {len(nice_values)} × {residual_observed:.2e} / {period:.4f} = {2*len(nice_values)*residual_observed/period:.6f}")

# Now specifically for 2/9
N_frac = len(nice_values)
p_specific = 2 * residual_observed / period  # probability of landing near ONE specific fraction
p_any = N_frac * p_specific  # union bound for any fraction
print(f"\np(near 2/9 specifically) = {p_specific:.2e}")
print(f"p(near any p/q, q≤20, union bound) = {p_any:.2e}")
if p_specific > 0:
    print(f"Sigma (specific to 2/9) = {norm.ppf(1-p_specific):.2f}")
    print(f"Sigma (any nice fraction) = {norm.ppf(1-p_any):.2f}")

# ============================================================
# Task 3: SU(5) anomaly coefficients
# ============================================================
print("\n" + "=" * 60)
print("TASK 3: SU(5) anomaly and Dynkin index computation")
print("=" * 60)

N = 5

# Anomaly coefficients A(R)
# For SU(N): A(fund) = 1, A(adj) = 0 (real rep), A(S²) = N+4, A(Λ²) = N-4
# A(R̄) = -A(R)
print("\nAnomaly coefficients A(R):")
print(f"  A(5)  = 1")
print(f"  A(5̄)  = -1")
print(f"  A(10) = A(Λ²(5)) = N-4 = {N-4}")
print(f"  A(10̄) = -(N-4) = {-(N-4)}")
print(f"  A(15) = A(S²(5)) = N+4 = {N+4}")
print(f"  A(15̄) = -(N+4) = {-(N+4)}")
print(f"  A(24) = A(adj) = 0  [adjoint is real for SU(N)]")
print(f"\n  A(24 ⊕ 15 ⊕ 15̄) = 0 + {N+4} + {-(N+4)} = 0  ✓ ANOMALY-FREE")

# Dynkin index T(R) = dim(R) * C₂(R) / dim(adj)
# Standard: T(fund) = 1/2, T(adj) = N, T(S²) = (N+2)/2, T(Λ²) = (N-2)/2
print(f"\nDynkin indices T(R):")
T_fund = 0.5
T_adj = N
T_S2 = (N + 2) / 2
T_A2 = (N - 2) / 2
print(f"  T(5)  = 1/2 = {T_fund}")
print(f"  T(10) = T(Λ²) = (N-2)/2 = {T_A2}")
print(f"  T(15) = T(S²) = (N+2)/2 = {T_S2}")
print(f"  T(24) = T(adj) = N = {T_adj}")

print(f"\n  T(24 ⊕ 15 ⊕ 15̄) = {T_adj} + {T_S2} + {T_S2} = {T_adj + 2*T_S2}")

# Beta function for SU(5) gauge theory
# b₀ = (11/3)C₂(G) - (2/3)T(matter) - (1/3)T(scalars)
# For pure gauge + matter in 24 ⊕ 15 ⊕ 15̄ (treating as Weyl fermions):
C2_G = N  # C₂(SU(N)) = N
T_matter = T_adj + 2 * T_S2  # 24 + 15 + 15̄ as fermions
b0_fermion = 11/3 * C2_G - 2/3 * T_matter
print(f"\nOne-loop beta function (SU(5) with 24⊕15⊕15̄ as fermions):")
print(f"  b₀ = (11/3)×{C2_G} - (2/3)×{T_matter}")
print(f"  b₀ = {11/3 * C2_G:.4f} - {2/3 * T_matter:.4f} = {b0_fermion:.4f}")
print(f"  Asymptotic freedom: {'YES' if b0_fermion > 0 else 'NO'} (b₀ {'>' if b0_fermion > 0 else '<'} 0)")

# For QCD: SU(3) with n_f flavors of fundamental fermions
# b₀ = 11 - 2n_f/3
# AF requires n_f < 16.5
print(f"\nFor comparison, QCD SU(3):")
print(f"  b₀ = 11 - 2n_f/3, AF requires n_f < 16.5")
print(f"  With n_f = 6: b₀ = {11 - 2*6/3:.1f}")

print("\n" + "=" * 60)
print("ALL TASKS COMPLETE")
print("=" * 60)
