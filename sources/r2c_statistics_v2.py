#!/usr/bin/env python3
"""Round 2-C: Statistical significance — no scipy dependency"""
import numpy as np
from itertools import combinations
import math, sys

def erfinv_approx(x):
    """Approximate inverse error function for p-value to sigma conversion."""
    # Winitzki approximation
    a = 0.147
    ln1mx2 = math.log(1 - x*x)
    s = np.sign(x)
    t = 2/(math.pi*a) + ln1mx2/2
    return s * math.sqrt(math.sqrt(t*t - ln1mx2/a) - t)

def p_to_sigma(p):
    """Convert one-sided p-value to sigma."""
    if p <= 0:
        return float('inf')
    if p >= 1:
        return 0
    # sigma such that P(X > sigma) = p for standard normal
    # = sqrt(2) * erfinv(1 - 2p)
    return math.sqrt(2) * erfinv_approx(1 - 2*p)

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

sm_masses = np.array([0.51099895, 105.6583755, 1776.86, 2.16, 4.67, 93.4, 1270, 4180, 172690])
sm_best = Q_best(sm_masses)
print(f"SM best |Q - 3/2| = {sm_best:.6e}")

N_MC = 10000
np.random.seed(42)

# Log-uniform in [0.1, 200000]
count_wide = 0
for i in range(N_MC):
    if i % 2000 == 0:
        print(f"  Wide range: {i}/{N_MC}...", file=sys.stderr)
    log_masses = np.random.uniform(np.log(0.1), np.log(200000), 9)
    masses = np.exp(log_masses)
    if Q_best(masses) <= sm_best:
        count_wide += 1

p_wide = count_wide / N_MC
sigma_wide = p_to_sigma(p_wide) if p_wide > 0 else "> 3.7"
print(f"\nWide range [0.1, 200000] MeV:")
print(f"  Count: {count_wide}/{N_MC}")
print(f"  p-value = {p_wide:.4e}")
print(f"  Equivalent sigma ~ {sigma_wide}")

# Log-uniform in [0.5, 180000]
count_tight = 0
np.random.seed(123)
for i in range(N_MC):
    if i % 2000 == 0:
        print(f"  Tight range: {i}/{N_MC}...", file=sys.stderr)
    log_masses = np.random.uniform(np.log(0.5), np.log(180000), 9)
    masses = np.exp(log_masses)
    if Q_best(masses) <= sm_best:
        count_tight += 1

p_tight = count_tight / N_MC
sigma_tight = p_to_sigma(p_tight) if p_tight > 0 else "> 3.7"
print(f"\nTight range [0.5, 180000] MeV:")
print(f"  Count: {count_tight}/{N_MC}")
print(f"  p-value = {p_tight:.4e}")
print(f"  Equivalent sigma ~ {sigma_tight}")

# ============================================================
# Task 2: Statistical significance of δ₀ mod 2π/3 ≈ 2/9
# ============================================================
print("\n" + "=" * 60)
print("TASK 2: Is δ₀ mod 2π/3 ≈ 2/9 statistically significant?")
print("=" * 60)

period = 2 * np.pi / 3

# Nice fractions p/q with q <= 20, values in [0, period)
nice_values = set()
for q in range(1, 21):
    for p in range(0, q):
        val = p / q
        if val < period:
            nice_values.add(val)
nice_arr = np.array(sorted(nice_values))
N_nice = len(nice_arr)

residual_observed = abs(0.222230 - 2/9)
print(f"Observed residual: |0.222230 - 2/9| = {residual_observed:.2e}")
print(f"Number of nice fractions p/q (q≤20) in [0, {period:.4f}): {N_nice}")

# Monte Carlo
N_MC2 = 100000
np.random.seed(456)
count_near = 0
for i in range(N_MC2):
    delta = np.random.uniform(0, period)
    dists = np.abs(delta - nice_arr)
    if np.min(dists) <= residual_observed:
        count_near += 1

p_nice = count_near / N_MC2
print(f"\nMonte Carlo ({N_MC2} trials):")
print(f"  Random δ₀ within {residual_observed:.2e} of ANY p/q (q≤20): {count_near}/{N_MC2}")
print(f"  p-value (any nice fraction) = {p_nice:.6f}")

# Analytic estimate
p_analytic_any = 2 * N_nice * residual_observed / period
p_analytic_specific = 2 * residual_observed / period
print(f"\nAnalytic estimates:")
print(f"  p(near 2/9 specifically) = 2×{residual_observed:.2e}/{period:.4f} = {p_analytic_specific:.2e}")
print(f"  sigma (specific to 2/9) ~ {p_to_sigma(p_analytic_specific):.1f}")
print(f"  p(near any p/q, q≤20) = {N_nice}× above = {p_analytic_any:.2e}")
print(f"  sigma (any nice fraction) ~ {p_to_sigma(p_analytic_any):.1f}")

# ============================================================
# Task 3: SU(5) anomaly coefficients
# ============================================================
print("\n" + "=" * 60)
print("TASK 3: SU(5) anomaly and Dynkin index computation")
print("=" * 60)

N = 5
print(f"\nAnomaly coefficients A(R) for SU({N}):")
print(f"  A(24) = 0  [adjoint is real]")
print(f"  A(15) = N+4 = {N+4}")
print(f"  A(15̄) = −{N+4}")
print(f"  Total: 0 + {N+4} + (−{N+4}) = 0  ✓ ANOMALY-FREE")

T_adj = N
T_S2 = (N + 2) / 2
print(f"\nDynkin indices T(R):")
print(f"  T(24) = {T_adj}, T(15) = T(15̄) = {T_S2}")
print(f"  T(total) = {T_adj + 2*T_S2}")

C2_G = N
T_matter = T_adj + 2 * T_S2
b0 = 11/3 * C2_G - 2/3 * T_matter
print(f"\nOne-loop beta (SU(5), 24⊕15⊕15̄ as fermions):")
print(f"  b₀ = (11/3)×{C2_G} − (2/3)×{T_matter} = {b0:.4f}")
print(f"  Asymptotic freedom: {'YES' if b0 > 0 else 'NO'}")

print("\n" + "=" * 60)
print("ALL TASKS COMPLETE")
print("=" * 60)
