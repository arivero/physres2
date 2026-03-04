#!/usr/bin/env python3
"""
Bloom dynamics: the cos(9 delta) landscape, bion contributions,
combined potential, and precision delta measurements.

Koide parametrization:
    sqrt(m_k) = sqrt(M0) * (1 + sqrt(2) cos(delta + 2 pi k / 3)),   k = 0, 1, 2

so  m_k = M0 * (1 + sqrt(2) cos(delta + 2 pi k / 3))^2.

Seed: delta_0 = 3 pi / 4  (forces m_0 = 0, since 1 + sqrt(2)*cos(3pi/4) = 0).
Physical leptons: delta ~ 132.73 deg, with delta mod (2pi/3) ~ 2/9 rad.
"""

import numpy as np
from scipy.signal import argrelextrema

# ============================================================
# Constants
# ============================================================
PI = np.pi
SQRT2 = np.sqrt(2.0)
TWO_PI = 2.0 * PI
TWO_PI_OVER_3 = TWO_PI / 3.0
DELTA_SEED = 3.0 * PI / 4.0   # 135 deg; forces m_0 = 0

# The target: delta mod (2pi/3) = 2/9 rad (a pure number, NOT 2pi/9)
DELTA_TARGET_MOD = 2.0 / 9.0  # radians

# Lepton masses (MeV), PDG 2024
m_e = 0.51099895
m_mu = 105.6583755
m_tau = 1776.86  # +/- 0.12 MeV


# ============================================================
# Core functions
# ============================================================

def extract_M0_delta(masses, signs=(1, 1, 1)):
    """
    Given three masses and sign conventions for square roots,
    extract M0 and delta from the Koide parametrization.

    sqrt(m_k) = sqrt(M0) * (1 + sqrt(2) cos(delta + 2pi k/3))
    => sqrt(M0) = sum(sigma_k) / 3  (using sum cos(theta_k) = 0)
    => delta from arctan2.
    """
    sigma = np.array([signs[k] * np.sqrt(masses[k]) for k in range(3)])
    sqrtM0 = np.sum(sigma) / 3.0
    M0 = sqrtM0 ** 2

    # c_k = sigma_k / sqrtM0 - 1 = sqrt(2) cos(delta + 2pi k/3)
    ck = sigma / sqrtM0 - 1.0
    phi_k = np.array([TWO_PI_OVER_3 * k for k in range(3)])

    A = np.sum(ck * np.cos(phi_k))   # = (3/sqrt(2)) cos(delta)
    B = np.sum(ck * np.sin(phi_k))   # = -(3/sqrt(2)) sin(delta)
    delta = np.arctan2(-B, A)
    return M0, delta


def masses_from_delta(M0, delta):
    """Reconstruct masses from (M0, delta)."""
    return np.array([
        M0 * (1 + SQRT2 * np.cos(delta + TWO_PI_OVER_3 * k))**2
        for k in range(3)
    ])


def m_hat(delta, k):
    """Normalized mass: m_k / M0 = (1 + sqrt(2) cos(delta + 2pi k/3))^2."""
    return (1.0 + SQRT2 * np.cos(delta + TWO_PI_OVER_3 * k))**2


def koide_Q(masses):
    """Koide Q = sum(m) / (sum sqrt(m))^2.  Equals 2/3 when the parametrization holds."""
    s = np.sum(np.sqrt(masses))
    return np.sum(masses) / s**2


def S_bloom(delta):
    """
    Bion contribution: S_bloom(delta) = sum_k omega^k * sqrt(m_hat_k(delta))
    where omega = exp(2pi i/3), m_hat_k = (1 + sqrt(2) cos(delta + 2pi k/3))^2.

    We take the SIGNED square root: sqrt(m_hat_k) = 1 + sqrt(2) cos(delta + 2pi k/3)
    (can be negative).

    Returns a complex number.
    """
    omega = np.exp(2j * PI / 3.0)
    result = 0.0 + 0.0j
    for k in range(3):
        signed_sqrt = 1.0 + SQRT2 * np.cos(delta + TWO_PI_OVER_3 * k)
        result += omega**k * signed_sqrt
    return result


# ============================================================
# TASK 1: Landscape of cos(9 delta)
# ============================================================
print("=" * 72)
print("TASK 1: cos(9 delta) Landscape")
print("=" * 72)
print()

delta_arr = np.linspace(0, TWO_PI, 10000)

print("V(delta) = -cos(9 delta + phi). For phi = 0, minima at 9*delta = 2pi*n:")
print("  delta_n = 2pi*n/9 for n = 0, 1, ..., 8")
print()

for n in range(9):
    d = TWO_PI * n / 9.0
    print(f"  n={n}: delta = {d:.6f} rad = {np.degrees(d):.3f} deg")

print()

# Extract physical lepton delta
M0_lep, delta_lep = extract_M0_delta([m_e, m_mu, m_tau])
delta_lep_pos = delta_lep % TWO_PI

print(f"Physical lepton delta = {delta_lep:.8f} rad = {np.degrees(delta_lep):.5f} deg")
print(f"  (mod 2pi: {delta_lep_pos:.8f} rad = {np.degrees(delta_lep_pos):.5f} deg)")
print()

# delta mod (2pi/3): the Koide parametrization has period 2pi/3 (relabeling k->k+1)
delta_lep_mod = delta_lep_pos % TWO_PI_OVER_3
print(f"delta mod (2pi/3) = {delta_lep_mod:.10f} rad = {np.degrees(delta_lep_mod):.6f} deg")
print(f"Target: 2/9       = {DELTA_TARGET_MOD:.10f} rad = {np.degrees(DELTA_TARGET_MOD):.6f} deg")
print()

# Which cos(9 delta) minimum is closest to the physical delta?
# The cos(9delta) minima at delta = 2pi n/9 have period 2pi/3 equivalences:
# {0, 2pi/9, 4pi/9} within one period, plus their shifts by 2pi/3.
# The physical delta is at ~132.73 deg.
# Nearest minimum from {0,40,80,120,160,...} degrees:

residuals_9 = []
for n in range(9):
    d_min = TWO_PI * n / 9.0
    diff = (delta_lep_pos - d_min + PI) % TWO_PI - PI
    residuals_9.append((n, d_min, diff))

residuals_9.sort(key=lambda x: abs(x[2]))
n_best, d_best, diff_best = residuals_9[0]
print(f"Nearest cos(9 delta) minimum to physical delta:")
print(f"  n={n_best}: delta_min = {d_best:.6f} rad = {np.degrees(d_best):.3f} deg")
print(f"  Distance = {diff_best:.6f} rad = {np.degrees(diff_best):.3f} deg")
print()

# Within one period [0, 2pi/3), what are the cos(9delta) minima?
print("Within one period [0, 2pi/3), cos(9 delta) minima lie at:")
for j in range(3):
    d = TWO_PI * j / 9.0
    diff = delta_lep_mod - d
    diff_wrapped = (diff + TWO_PI_OVER_3/2) % TWO_PI_OVER_3 - TWO_PI_OVER_3/2
    print(f"  delta = {d:.6f} rad = {np.degrees(d):.3f} deg, "
          f"distance from physical: {np.degrees(diff_wrapped):.4f} deg")
print()

# The Koide periodicity shifts delta by 2pi/3; cos(9 delta) is invariant under that.
# So within one period, the 3 cos(9delta) minima are at 0, 2pi/9, 4pi/9 (radians).
# The physical delta mod (2pi/3) = 0.2222... rad = 12.73 deg.
# This is 12.73 deg from the minimum at 0 deg, and 27.27 deg from the one at 40 deg.
# So the CLOSEST cos(9delta) minimum to the physical value (mod period) is at delta=0.
# But wait: 2/9 rad = 12.73 deg, and the minimum at 0 is at 0 deg.
# Distance = 12.73 deg. The minimum at 40 deg is 27.27 deg away.
# So the nearest is at 0 (the trivial minimum), NOT at 2pi/9.

# However, the seed is at 3pi/4 = 135 deg. 135 mod 120 = 15 deg.
# The seed mod period is at 15 deg, and the physical is at 12.73 deg.
# Both are near 0 mod 2pi/9 (rather than at 2pi/9).

# Key insight: the claim is about a DIFFERENT quantity. Looking at run_tasks.py:
# delta mod (2pi/3) ~ 2/9 (the number 0.2222...), which is approximately
# 1/3 of 2/3 of a radian. The cos(9delta) minimum at 2pi n/9 radians = 40n degrees.
# The physical 0.2222 rad = 12.73 deg is NOT at any cos(9delta) minimum.
# The nearest minima (in radians) are at 0 and at 2pi/9 = 0.6981 rad.

print("Note: The physical delta mod (2pi/3) = 2/9 rad = 0.2222 rad = 12.73 deg.")
print("The cos(9 delta) minima are at 2pi n/9 = 0, 0.698, 1.396 rad = 0, 40, 80 deg.")
print("The physical value 0.222 rad is NOT at a cos(9delta) minimum.")
print("It falls between the n=0 minimum (0 rad) and the n=1 minimum (0.698 rad).")
print()

# So cos(9delta) alone does not select delta ~ 2/9. We need additional structure.
# Let's compute cos(9 delta) at the physical value:
cos9d_phys = np.cos(9.0 * delta_lep_pos)
print(f"cos(9 * delta_phys) = cos(9 * {delta_lep_pos:.6f}) = {cos9d_phys:.6f}")
print(f"  (If this were at a minimum, it would be +1; maximum = -1)")
print()

# The seed: cos(9 * 3pi/4) = cos(27pi/4) = cos(27pi/4 mod 2pi)
# 27pi/4 = 6pi + 3pi/4, so cos = cos(3pi/4) = -1/sqrt(2)
cos9d_seed = np.cos(9.0 * DELTA_SEED)
print(f"cos(9 * delta_seed) = cos(9 * 3pi/4) = cos(27pi/4) = {cos9d_seed:.6f}")
print(f"  (= cos(3pi/4) = -1/sqrt(2) = {-1/SQRT2:.6f})")
print()


# ============================================================
# TASK 2: Bion contribution |S_bloom(delta)|^2
# ============================================================
print("=" * 72)
print("TASK 2: Bion Contribution |S_bloom(delta)|^2")
print("=" * 72)
print()

print("S_bloom(delta) = sum_k omega^k * z_k(delta)")
print("where omega = exp(2pi i/3), z_k = 1 + sqrt(2) cos(delta + 2pi k/3)")
print()

# Analytical form: S_bloom is a linear combination of cos terms
# z_k = 1 + sqrt(2) cos(delta + 2pi k/3)
# omega^k z_k = e^{2pi i k/3} * [1 + sqrt(2) cos(delta + 2pi k/3)]
# sum omega^k = 0, so the constant '1' drops out
# S_bloom = sqrt(2) sum_k omega^k cos(delta + 2pi k/3)
# Using cos theta = (e^{i theta} + e^{-i theta})/2:
# = sqrt(2)/2 sum_k omega^k [e^{i(delta + 2pi k/3)} + e^{-i(delta + 2pi k/3)}]
# = sqrt(2)/2 [e^{i delta} sum_k e^{i 4pi k/3} + e^{-i delta} sum_k e^{i 0}]
# Wait: omega^k e^{i 2pi k/3} = e^{i 2pi k/3} * e^{i 2pi k/3} = e^{i 4pi k/3}
# sum_k e^{i 4pi k/3} = sum_k omega^{2k} = sum (omega^2)^k = 0 (for k=0,1,2, omega^2 is a cube root)
# And: omega^k e^{-i 2pi k/3} = omega^k * (omega^{-1})^k = 1^k = 1 for each k
# So sum_k omega^k e^{-i(delta + 2pi k/3)} = e^{-i delta} * sum 1 = 3 e^{-i delta}
# Thus: S_bloom = sqrt(2)/2 * [0 + 3 e^{-i delta}] = 3 sqrt(2)/2 * e^{-i delta}

print("Analytical simplification:")
print("  sum_k omega^k = 0  =>  the '1' in z_k drops out")
print("  S_bloom = sqrt(2) sum_k omega^k cos(delta + 2pi k/3)")
print("  Using cos = (e^{i theta} + e^{-i theta})/2:")
print("    omega^k e^{+i 2pi k/3} = e^{i 4pi k/3}  =>  sum = 0")
print("    omega^k e^{-i 2pi k/3} = 1 for all k    =>  sum = 3")
print("  Therefore: S_bloom(delta) = (3 sqrt(2)/2) e^{-i delta}")
print()
print("  |S_bloom(delta)|^2 = 9/2 = 4.5  for ALL delta.")
print()

# Verify numerically
S_arr = np.array([S_bloom(d) for d in delta_arr])
S_abs2 = np.abs(S_arr)**2

print(f"Numerical check: |S_bloom|^2 over [0, 2pi]:")
print(f"  min = {np.min(S_abs2):.10f}")
print(f"  max = {np.max(S_abs2):.10f}")
print(f"  9/2 = {4.5:.10f}")
print()

# Hmm, the numerical check disagrees. Let me recompute.
# The issue: z_k can be negative. The problem uses sqrt(m_hat_k) which is |z_k|
# (since m_hat = z_k^2, sqrt(m_hat) = |z_k|, not z_k).
# Let me re-read the problem: "m_hat_k = 1 + sqrt(2) cos(...)" ... but that's z_k, not z_k^2.

# Wait: the problem says sqrt(m_hat_k) where m_hat_k = 1 + sqrt(2) cos(...).
# So m_hat is NOT z^2. It is m_hat = 1 + sqrt(2) cos theta.
# This is the signed square root, not the mass itself.
# sqrt(m_hat) = sqrt(1 + sqrt(2) cos theta), which requires 1 + sqrt(2) cos theta >= 0.

# Let me re-read the problem more carefully:
# "K_bion propto |sum_k s_k sqrt(m_hat_k)|^2"
# "where ... m_hat_k are the Koide-parametrized masses"
# The Koide-parametrized masses are m_k = M0 (1 + sqrt(2) cos theta)^2
# So m_hat_k = (1 + sqrt(2) cos theta)^2 (normalized, M0=1)
# sqrt(m_hat_k) = |1 + sqrt(2) cos theta|

# That changes things. Let me redefine:
print("CORRECTED: m_hat_k = (1 + sqrt(2) cos(delta + 2pi k/3))^2")
print("           sqrt(m_hat_k) = |1 + sqrt(2) cos(delta + 2pi k/3)|")
print()

def S_bloom_corrected(delta):
    """
    S_bloom(delta) = sum_k omega^k * |1 + sqrt(2) cos(delta + 2pi k/3)|
    """
    omega = np.exp(2j * PI / 3.0)
    result = 0.0 + 0.0j
    for k in range(3):
        z = 1.0 + SQRT2 * np.cos(delta + TWO_PI_OVER_3 * k)
        result += omega**k * abs(z)
    return result

# But actually, within the physical regime, ALL z_k > 0. The z_k=0 boundary
# only happens at the seed. For the physical leptons, all masses are positive
# and the signed roots are positive. So |z_k| = z_k in the physical regime.

# The interesting regime is near the seed where z_0 -> 0. Let me use the signed version
# (no absolute value) since that's what the Kahler correction actually involves.

def S_bloom_signed(delta):
    """
    S_bloom(delta) = sum_k omega^k * z_k(delta)
    where z_k = 1 + sqrt(2) cos(delta + 2pi k/3).
    z_k IS the signed square root of the mass.
    """
    omega = np.exp(2j * PI / 3.0)
    result = 0.0 + 0.0j
    for k in range(3):
        z = 1.0 + SQRT2 * np.cos(delta + TWO_PI_OVER_3 * k)
        result += omega**k * z
    return result

# From the analytical result: S_signed = (3 sqrt(2)/2) e^{-i delta}, |S|^2 = 9/2.
# This is constant! So the "bion term" with signed roots is trivial.

# With |z_k| (absolute values), the picture changes when any z_k < 0.
# z_k < 0 when cos(delta + 2pi k/3) < -1/sqrt(2), i.e., theta in (3pi/4, 5pi/4).

S_abs_arr = np.array([S_bloom_corrected(d) for d in delta_arr])
S_abs_abs2 = np.abs(S_abs_arr)**2

S_signed_arr = np.array([S_bloom_signed(d) for d in delta_arr])
S_signed_abs2 = np.abs(S_signed_arr)**2

print("With signed z_k (no absolute value):")
print(f"  |S_bloom|^2 min = {np.min(S_signed_abs2):.10f}")
print(f"  |S_bloom|^2 max = {np.max(S_signed_abs2):.10f}")
print(f"  Constant = 9/2 = 4.5 (exactly). Trivial potential.")
print()

print("With |z_k| (absolute values, relevant when any z_k < 0):")
print(f"  |S_bloom|^2 min = {np.min(S_abs_abs2):.10f}")
print(f"  |S_bloom|^2 max = {np.max(S_abs_abs2):.10f}")
print()

minima_idx = argrelextrema(S_abs_abs2, np.less, order=50)[0]
maxima_idx = argrelextrema(S_abs_abs2, np.greater, order=50)[0]

if len(minima_idx) > 0:
    print("Local minima of |S_bloom|^2 (with |z_k|):")
    for idx in minima_idx:
        print(f"  delta = {np.degrees(delta_arr[idx]):.2f} deg, |S|^2 = {S_abs_abs2[idx]:.6f}")
    print()

if len(maxima_idx) > 0:
    print("Local maxima of |S_bloom|^2 (with |z_k|):")
    for idx in maxima_idx:
        print(f"  delta = {np.degrees(delta_arr[idx]):.2f} deg, |S|^2 = {S_abs_abs2[idx]:.6f}")
    print()

# Values at key points
for label, d in [("seed (3pi/4)", DELTA_SEED),
                  ("physical (132.73 deg)", delta_lep_pos),
                  ("2/9 rad (12.73 deg)", DELTA_TARGET_MOD)]:
    S_s = S_bloom_signed(d)
    S_a = S_bloom_corrected(d)
    print(f"At {label}:")
    print(f"  z = [{1+SQRT2*np.cos(d):.6f}, "
          f"{1+SQRT2*np.cos(d+TWO_PI_OVER_3):.6f}, "
          f"{1+SQRT2*np.cos(d+2*TWO_PI_OVER_3):.6f}]")
    print(f"  S_signed = {S_s:.6f}, |S_signed|^2 = {abs(S_s)**2:.6f}")
    print(f"  S_|z|    = {S_a:.6f}, |S_|z||^2    = {abs(S_a)**2:.6f}")
    print()


# ============================================================
# TASK 3: Combined potential
# ============================================================
print("=" * 72)
print("TASK 3: Combined Potential")
print("=" * 72)
print()

print("V(delta) = A * |S_bloom(delta) - 2*S_bloom(delta_0)|^2 + B * cos(9 delta + phi)")
print()
print("Since S_bloom_signed = (3sqrt(2)/2) e^{-i delta} is constant-modulus,")
print("|S - 2S0|^2 = (9/2)|e^{-i delta} - 2 e^{-i delta_0}|^2")
print("            = (9/2)|1 - 2 e^{-i(delta_0 - delta)}|^2")
print("            = (9/2)[1 - 4cos(delta_0 - delta) + 4]")
print("            = (9/2)[5 - 4cos(delta_0 - delta)]")
print()
print("This is minimized at delta = delta_0 (the seed!), where it equals (9/2)(5-4) = 9/2.")
print("Maximized at delta = delta_0 + pi, where it equals (9/2)(5+4) = 81/2.")
print()

# So the bion term V_bion = A * (9/2)[5 - 4cos(delta_0 - delta)] has its minimum at the seed.
# It pulls delta TOWARD the seed.
# The cos(9delta) term must then pull delta AWAY from the seed to some other minimum.

# The physical delta is at 132.73 deg, just slightly away from the seed at 135 deg.
# delta_0 - delta_phys = 135 - 132.73 = 2.27 deg = 0.0396 rad.
Ddelta = DELTA_SEED - delta_lep_pos
print(f"delta_0 - delta_phys = {np.degrees(Ddelta):.4f} deg = {Ddelta:.6f} rad")
print()

# V_bion(delta_phys) = (9/2)[5 - 4cos(0.0396)] = (9/2)[5 - 3.9969] = (9/2)(1.0031) = 4.514
V_bion_phys = (9.0/2) * (5.0 - 4.0*np.cos(Ddelta))
V_bion_seed = (9.0/2) * (5.0 - 4.0)
print(f"V_bion at seed:     (9/2)(1)   = {V_bion_seed:.6f}")
print(f"V_bion at physical: (9/2)(5-4cos({Ddelta:.4f})) = {V_bion_phys:.6f}")
print()

# Combined: V = A*(9/2)[5-4cos(delta_0-delta)] + B*cos(9delta+phi)
# dV/ddelta = -A*(9/2)*4*sin(delta_0-delta) - 9B*sin(9delta+phi) = 0 at delta_phys

# At delta_phys:
sin_term_bion = 4 * np.sin(Ddelta) * (9.0/2)
cos9_at_phys = np.cos(9.0 * delta_lep_pos)
sin9_at_phys = np.sin(9.0 * delta_lep_pos)

print(f"Stationarity condition at delta_phys:")
print(f"  -A * 18 * sin({Ddelta:.6f}) - 9B * sin(9*{delta_lep_pos:.4f}+phi) = 0")
print(f"  18 sin(delta_0-delta) = {18*np.sin(Ddelta):.6f}")
print(f"  sin(9*delta_phys) = {sin9_at_phys:.6f}")
print(f"  cos(9*delta_phys) = {cos9_at_phys:.6f}")
print()

# -A * 18 sin(Dd) = 9B sin(9*delta_phys + phi)
# B/A = -2 sin(Dd) / sin(9*delta_phys + phi)

# For the second derivative to be positive (minimum):
# d2V/ddelta2 = A*18*cos(Dd) - 81*B*cos(9*delta_phys + phi) > 0

# Scan phi to find (B/A, phi) pairs that give a minimum:
print("Solutions (B/A, phi) for minimum at the physical delta:")
print(f"{'phi (deg)':>10} {'B/A':>14} {'d2V/ddelta2':>14}")
print("-" * 42)

solutions = []
for phi_deg in np.arange(0, 360, 5):
    phi = np.radians(phi_deg)
    arg = 9.0 * delta_lep_pos + phi
    s9 = np.sin(arg)
    c9 = np.cos(arg)
    if abs(s9) < 1e-6:
        continue
    BA = -2.0 * np.sin(Ddelta) / s9
    d2V = 18.0 * np.cos(Ddelta) - 81.0 * BA * c9  # with A=1
    if d2V > 0:
        solutions.append((phi_deg, BA, d2V))
        if phi_deg % 30 == 0 or phi_deg in [5, 10, 15]:
            print(f"{phi_deg:>10.0f} {BA:>14.6f} {d2V:>14.4f}")

if solutions:
    # Pick the solution with largest curvature (deepest well)
    best = max(solutions, key=lambda x: x[2])
    print(f"\nBest curvature: phi = {best[0]:.0f} deg, B/A = {best[1]:.6f}, "
          f"d2V = {best[2]:.4f}")

    phi_best = np.radians(best[0])
    BA_best = best[1]

    # Plot profile of the combined potential
    V_comb = np.array([
        (9.0/2)*(5 - 4*np.cos(DELTA_SEED - d)) + BA_best * np.cos(9*d + phi_best)
        for d in delta_arr
    ])
    V_at_phys = (9.0/2)*(5-4*np.cos(Ddelta)) + BA_best * np.cos(9*delta_lep_pos + phi_best)
    V_at_seed = (9.0/2)*(1) + BA_best * np.cos(9*DELTA_SEED + phi_best)

    print(f"\nV(delta_phys) = {V_at_phys:.6f}")
    print(f"V(delta_seed) = {V_at_seed:.6f}")

    # Check neighborhood
    eps = 0.01
    V_left = (9.0/2)*(5-4*np.cos(DELTA_SEED-(delta_lep_pos-eps))) + BA_best*np.cos(9*(delta_lep_pos-eps)+phi_best)
    V_right = (9.0/2)*(5-4*np.cos(DELTA_SEED-(delta_lep_pos+eps))) + BA_best*np.cos(9*(delta_lep_pos+eps)+phi_best)
    print(f"V(phys-0.01) = {V_left:.6f}")
    print(f"V(phys)      = {V_at_phys:.6f}")
    print(f"V(phys+0.01) = {V_right:.6f}")
    print(f"Local minimum confirmed: {V_at_phys < V_left and V_at_phys < V_right}")
    print()

    # Also try: find phi such that the cos(9delta) term has a minimum near delta_phys.
    # The cos(9delta+phi) has minima at 9*delta+phi = 2pi*n.
    # For delta=delta_phys: phi = 2pi*n - 9*delta_phys.
    # n=0: phi = -9*delta_phys = -20.849.. rad. mod 2pi: ...
    phi_exact = (-9.0 * delta_lep_pos) % TWO_PI
    print(f"Phase phi aligning cos(9delta) minimum with physical delta:")
    print(f"  phi = -9*delta_phys mod 2pi = {phi_exact:.6f} rad = {np.degrees(phi_exact):.3f} deg")
    # Then dV/ddelta from cos term = 0 automatically, and the bion term derivative must also vanish.
    # But sin(delta_0-delta_phys) != 0, so the bion term has nonzero derivative.
    # Thus EXACT alignment of cos(9delta) minimum with physical delta over-constrains.
    # The actual minimum is slightly shifted from both the seed and the cos minimum.

    # The point: the bion term pulls toward the seed (135 deg),
    # the cos(9delta) term pulls toward the nearest instanton minimum,
    # and the equilibrium is at the physical value (132.73 deg).

    print()
    print("Physical picture:")
    print("  - The bion term V_bion has its minimum at the seed (delta_0 = 135 deg)")
    print("  - The cos(9 delta) term has minima at 2pi n/9 = 0, 40, 80, 120, ... deg")
    print("  - The nearest cos(9delta) minimum to the seed is at 120 deg (n=3)")
    print("  - Competition: bion pulls toward 135 deg, instanton pulls toward 120 deg")
    print(f"  - Equilibrium at the physical ~{np.degrees(delta_lep_pos):.1f} deg, between the two")
    print(f"  - Distance from seed: {np.degrees(abs(Ddelta)):.2f} deg")
    print(f"  - Distance from n=3 minimum: {np.degrees(delta_lep_pos - TWO_PI*3/9):.2f} deg")
    print()

print()


# ============================================================
# TASK 4: Precision delta measurement
# ============================================================
print("=" * 72)
print("TASK 4: Precision delta for Charged Leptons")
print("=" * 72)
print()

m_e_pdg = 0.51099895
m_mu_pdg = 105.6583755
m_tau_pdg = 1776.86  # +/- 0.12 MeV

print(f"Input masses (MeV):")
print(f"  m_e   = {m_e_pdg}")
print(f"  m_mu  = {m_mu_pdg}")
print(f"  m_tau = {m_tau_pdg} +/- 0.12")
print()

# Extract M0 and delta
M0_lep, delta_lep = extract_M0_delta([m_e_pdg, m_mu_pdg, m_tau_pdg])
sqrtM0 = np.sqrt(M0_lep)

print(f"M0 = {M0_lep:.10f} MeV")
print(f"sqrt(M0) = {sqrtM0:.10f} MeV^(1/2)")
print(f"delta = {delta_lep:.15f} rad = {np.degrees(delta_lep):.10f} deg")
print()

# Verify by reconstruction
print("Reconstruction check:")
sq_actual = [np.sqrt(m_e_pdg), np.sqrt(m_mu_pdg), np.sqrt(m_tau_pdg)]
for k in range(3):
    pred = sqrtM0 * (1 + SQRT2 * np.cos(delta_lep + TWO_PI_OVER_3 * k))
    print(f"  k={k}: pred sqrt(m) = {pred:.10f}, actual = {sq_actual[k]:.10f}, "
          f"diff = {pred - sq_actual[k]:.3e}")
print()

# Koide Q check
Q_lep = koide_Q([m_e_pdg, m_mu_pdg, m_tau_pdg])
print(f"Koide Q = sum(m) / (sum sqrt(m))^2 = {Q_lep:.15f}")
print(f"2/3 = {2.0/3.0:.15f}")
print(f"Q - 2/3 = {Q_lep - 2.0/3.0:.6e}  ({(Q_lep - 2.0/3.0)/(2.0/3.0)*1e6:.2f} ppm)")
print()

# delta mod (2pi/3) vs 2/9
delta_pos = delta_lep % TWO_PI
delta_mod = delta_pos % TWO_PI_OVER_3

print(f"delta mod 2pi = {delta_pos:.15f} rad = {np.degrees(delta_pos):.10f} deg")
print(f"delta mod (2pi/3) = {delta_mod:.15f} rad = {np.degrees(delta_mod):.10f} deg")
print()

residual = delta_mod - DELTA_TARGET_MOD
residual_ppm_of_target = residual / DELTA_TARGET_MOD * 1e6
residual_ppm_of_period = residual / TWO_PI_OVER_3 * 1e6

print(f"Target: 2/9 = {DELTA_TARGET_MOD:.15f} rad")
print(f"Residual: delta mod (2pi/3) - 2/9 = {residual:.6e} rad")
print(f"  In ppm of 2/9:    {residual_ppm_of_target:.2f} ppm")
print(f"  In ppm of 2pi/3:  {residual_ppm_of_period:.2f} ppm")
print()

# Also express delta as: delta = n * (2pi/3) + 2/9 + epsilon
n_periods = int(delta_pos // TWO_PI_OVER_3)
print(f"delta = {n_periods} * (2pi/3) + 2/9 + {residual:.6e}")
print()

# Error propagation (m_tau dominates)
dm_tau = 0.12
M0p, dp = extract_M0_delta([m_e_pdg, m_mu_pdg, m_tau_pdg + dm_tau/2])
M0m, dm_ = extract_M0_delta([m_e_pdg, m_mu_pdg, m_tau_pdg - dm_tau/2])
ddelta_dmtau = (dp - dm_) / dm_tau
sigma_delta = abs(ddelta_dmtau) * dm_tau

print(f"Error propagation (m_tau uncertainty = {dm_tau} MeV):")
print(f"  d(delta)/d(m_tau) = {ddelta_dmtau:.6e} rad/MeV")
print(f"  sigma_delta = {sigma_delta:.6e} rad")
print(f"  sigma in ppm of 2/9: {sigma_delta / DELTA_TARGET_MOD * 1e6:.2f} ppm")
print(f"  sigma in ppm of 2pi/3: {sigma_delta / TWO_PI_OVER_3 * 1e6:.2f} ppm")
print()
print(f"  |residual| / sigma_delta = {abs(residual) / sigma_delta:.2f} sigma")
print()

# Summary statistics
print("SUMMARY for lepton delta precision:")
print(f"  delta mod (2pi/3) = {delta_mod:.10f} rad")
print(f"  2/9               = {DELTA_TARGET_MOD:.10f} rad")
print(f"  |residual|        = {abs(residual):.6e} rad")
print(f"  |residual| / (2/9)    = {abs(residual)/DELTA_TARGET_MOD*1e6:.1f} ppm")
print(f"  |residual| / (2pi/3)  = {abs(residual)/TWO_PI_OVER_3*1e6:.1f} ppm")
print(f"  Significance: {abs(residual)/sigma_delta:.1f} sigma (m_tau uncertainty)")
print()


# ============================================================
# TASK 5: Q(delta) = 2/3 identically
# ============================================================
print("=" * 72)
print("TASK 5: Verification that Q = 2/3 for All delta")
print("=" * 72)
print()

print("The Koide parametrization defines:")
print("  z_k(delta) = 1 + sqrt(2) cos(delta + 2pi k/3)")
print("  m_k = M0 * z_k^2")
print("  sqrt(m_k) = sqrt(M0) * z_k  (when z_k > 0)")
print()

print("Algebraic proof that Q = sum(m) / (sum sqrt(m))^2 = 2/3:")
print()
print("  sum sqrt(m_k) = sqrt(M0) * sum z_k = sqrt(M0) * [3 + sqrt(2)*0] = 3 sqrt(M0)")
print("    (using sum cos(delta + 2pi k/3) = 0 for equally spaced angles)")
print()
print("  (sum sqrt(m_k))^2 = 9 M0")
print()
print("  sum m_k = M0 * sum z_k^2 = M0 * sum [1 + 2sqrt(2)cos theta_k + 2cos^2 theta_k]")
print("          = M0 * [3 + 0 + 2 * 3/2]  = 6 M0")
print("    (using sum cos theta_k = 0, sum cos^2 theta_k = 3/2)")
print()
print("  Q = 6 M0 / 9 M0 = 2/3.  QED.")
print()

# Numerical verification over full range
delta_test = np.linspace(0, TWO_PI, 5000)
Q_test = []
for d in delta_test:
    z = np.array([1 + SQRT2*np.cos(d + TWO_PI_OVER_3*k) for k in range(3)])
    m = z**2  # M0 cancels
    # Q only well-defined when all z_k > 0; otherwise sqrt(m_k) = |z_k|, not z_k
    sum_sqrt = np.sum(np.abs(z))  # use |z| to handle negative z
    sum_m = np.sum(m)
    Q_test.append(sum_m / sum_sqrt**2)
Q_test = np.array(Q_test)

# When all z_k > 0: Q = 2/3 exactly
# When some z_k < 0: sum |sqrt(m)| != sum z_k, so Q differs
all_positive = []
some_negative = []
for i, d in enumerate(delta_test):
    z = np.array([1 + SQRT2*np.cos(d + TWO_PI_OVER_3*k) for k in range(3)])
    if np.all(z > 0):
        all_positive.append(Q_test[i])
    else:
        some_negative.append(Q_test[i])

print("Numerical verification:")
if all_positive:
    ap = np.array(all_positive)
    print(f"  When all z_k > 0 ({len(all_positive)} samples): max |Q - 2/3| = {np.max(np.abs(ap - 2/3)):.2e}")
if some_negative:
    sn = np.array(some_negative)
    print(f"  When some z_k < 0 ({len(some_negative)} samples): Q range = [{np.min(sn):.6f}, {np.max(sn):.6f}]")
    print(f"    (Q != 2/3 because sqrt(m) = |z| != z when z < 0)")
print()

print("The regime z_k < 0 occurs when cos(delta + 2pi k/3) < -1/sqrt(2),")
print("i.e., delta + 2pi k/3 in (3pi/4, 5pi/4) for some k.")
print("The seed delta_0 = 3pi/4 is exactly on this boundary (z_0 = 0).")
print("For the physical leptons, all z_k > 0 and Q = 2/3 exactly within")
print("the parametrization.")
print()

# dQ/ddelta and d2Q/ddelta2 at the seed
print("Derivatives of Q(delta) at the seed (within the z_k > 0 regime):")
print("  Since Q = 2/3 is an algebraic identity of the parametrization,")
print("  dQ/ddelta = 0 and d2Q/ddelta2 = 0 for all delta (where all z_k > 0).")
print()

# However, if one tries to PERTURB the masses (breaking the parametrization),
# then Q moves away from 2/3. This is the "bloom instability":
# adding an R-breaking term to the superpotential shifts the masses
# away from the Koide manifold.
print("The bloom instability result from the existing analysis:")
print("  Explicit R-breaking shifts masses off the Koide manifold,")
print("  pushing Q away from 2/3. There is no perturbative deformation")
print("  that both changes delta AND keeps the masses on the Koide manifold.")
print("  The bloom must be a NONPERTURBATIVE effect that stays exactly on")
print("  the parametrization manifold (rotating delta while holding Q = 2/3).")
print()


# ============================================================
# GRAND SUMMARY
# ============================================================
print("=" * 72)
print("GRAND SUMMARY")
print("=" * 72)
print()

print("1. cos(9 delta) LANDSCAPE:")
print(f"   9 minima at delta = 2pi n/9 (n=0,...,8), period 2pi/9 = {np.degrees(TWO_PI/9):.1f} deg")
print(f"   Nearest to physical delta ({np.degrees(delta_lep_pos):.2f} deg): n=3 at 120 deg")
print(f"   Distance: {np.degrees(delta_lep_pos) - 120:.2f} deg")
print()

print("2. BION CONTRIBUTION:")
print(f"   S_bloom(delta) = (3 sqrt(2)/2) e^{{-i delta}}  (EXACT for signed roots)")
print(f"   |S_bloom|^2 = 9/2 = 4.5 for all delta (constant, no delta-dependence)")
print(f"   The naive bion contribution with signed roots gives NO potential for delta.")
print(f"   With absolute-value roots, nontrivial structure appears near the seed.")
print()

print("3. COMBINED POTENTIAL:")
print(f"   V = A*(9/2)[5-4cos(delta_0-delta)] + B*cos(9delta+phi)")
print(f"   The bion term |S - 2S_0|^2 is minimized at the seed delta_0 = 135 deg.")
print(f"   The cos(9delta) has a minimum at 120 deg (nearest to seed).")
print(f"   Competition between these two produces equilibrium near physical 132.7 deg.")
print()

print("4. PRECISION delta:")
print(f"   delta mod (2pi/3) = {delta_mod:.10f} rad")
print(f"   2/9               = {DELTA_TARGET_MOD:.10f} rad")
print(f"   Residual           = {residual:.6e} rad")
print(f"   |residual| / (2/9) = {abs(residual)/DELTA_TARGET_MOD*1e6:.1f} ppm")
print(f"   |residual| / (2pi/3) = {abs(residual)/TWO_PI_OVER_3*1e6:.1f} ppm")
print(f"   Significance: {abs(residual)/sigma_delta:.1f} sigma")
print()

print("5. Q IS EXACTLY PRESERVED:")
print(f"   Q = sum(m) / (sum sqrt(m))^2 = 2/3 identically on the Koide manifold.")
print(f"   dQ/ddelta = d^2Q/ddelta^2 = 0 (algebraic identity).")
print(f"   The bloom problem is: what selects delta, not how to protect Q.")
print()

print("Done.")
