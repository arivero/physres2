#!/usr/bin/env python3
"""
Bion potential on the Koide manifold V_bion(delta).

Setup:
  Koide parametrization:  z_k = v0 * (1 + sqrt(2) cos(delta + 2 pi k/3))
  Mass: m_k = z_k^2,  with v0 = (sum z_k)/3
  Q = sum(m_k)/(sum z_k)^2 = 2/3  identically.

  The bion potential on the Koide manifold is:
    V_bion(delta) = (zeta^2/Lambda^2) |sum_k s_k(delta) sqrt(m_k(delta))|^2

  where s_k = sign(z_k) -- i.e. s_k = +1 if z_k > 0, s_k = -1 if z_k < 0.
  Since sqrt(m_k) = |z_k|, we have s_k sqrt(m_k) = z_k when z_k != 0.

  BUT the sign assignment flips when z_k passes through zero. The bion amplitude
  sees the PHYSICAL square root |z_k| with a sign that is locked to the
  monopole-instanton sector. When z_k crosses zero (at the seed), the sign
  s_k must be updated. This means:

    V_bion(delta) propto |sum_k s_k(delta) |z_k(delta)||^2

  The key question: is s_k(delta) = sign(z_k(delta)) at all times (adiabatic tracking),
  or does the sign lag behind a zero-crossing (non-adiabatic)?

  Case A (adiabatic): s_k = sign(z_k) always.
    Then sum s_k |z_k| = sum z_k = 3 v0 = constant on Koide manifold.
    V_bion = (zeta/Lambda)^2 * 9 v0^2 = CONSTANT.  No delta dependence. Trivial.

  Case B (non-adiabatic / sign-frozen): Signs are fixed in a given topological sector.
    E.g. for the (-s,c,b) triple: s = (+1, +1, +1) in the sector where all monopole
    charges are positive, regardless of the sign of z_k.
    Then sum s_k |z_k| = sum |z_k|, and V_bion = (zeta/Lambda)^2 * (sum |z_k|)^2.
    This DOES depend on delta.

  We compute both cases and also the mixed case where signs track z_k
  except at the specific crossing points.
"""

import numpy as np
from scipy.signal import argrelextrema
from scipy.optimize import minimize_scalar

# ============================================================
# Constants
# ============================================================
PI = np.pi
SQRT2 = np.sqrt(2.0)
TWO_PI = 2.0 * PI
TWO_PI_OVER_3 = TWO_PI / 3.0
DELTA_SEED = 3.0 * PI / 4.0  # 135 deg; forces z_0 = 0

# Physical masses (MeV)
m_e   = 0.51099895
m_mu  = 105.6583755
m_tau = 1776.86

m_s = 93.4
m_c = 1270.0
m_b = 4180.0
m_t = 172760.0


# ============================================================
# Core Koide functions
# ============================================================

def z_k(delta, k):
    """Signed square root: z_k = 1 + sqrt(2) cos(delta + 2pi k/3).
    Normalized (v0 = 1 effectively, i.e. we work with z/v0)."""
    return 1.0 + SQRT2 * np.cos(delta + TWO_PI_OVER_3 * k)

def z_all(delta):
    """All three z_k values."""
    return np.array([z_k(delta, k) for k in range(3)])

def m_k_normalized(delta, k):
    """Normalized mass m_k/v0^2 = z_k^2."""
    zk = z_k(delta, k)
    return zk**2

def extract_v0_delta(signed_roots):
    """Given signed roots (z_0, z_1, z_2), extract v0 and delta."""
    s = np.array(signed_roots, dtype=float)
    v0 = np.sum(s) / 3.0
    phi_k = np.array([TWO_PI_OVER_3 * k for k in range(3)])
    # z_k/v0 = 1 + sqrt(2) cos(delta + phi_k)
    # (z_k/v0 - 1) = sqrt(2) cos(delta + phi_k)
    c = s / v0 - 1.0
    A = np.sum(c * np.cos(phi_k))  # = (3/sqrt(2)) cos(delta) [from orthogonality]
    B = np.sum(c * np.sin(phi_k))  # = -(3/sqrt(2)) sin(delta)
    delta = np.arctan2(-B, A)
    return v0, delta

def koide_Q(masses, signs=None):
    """Koide quotient with signed roots.
    Q = sum(m_k) / (sum sigma_k)^2 where sigma_k = sign_k * sqrt(m_k)."""
    if signs is None:
        signs = np.ones(len(masses))
    sigma = np.array(signs) * np.sqrt(np.abs(masses))
    return np.sum(masses) / np.sum(sigma)**2


# ============================================================
# Bion potential functions
# ============================================================

def V_bion_adiabatic(delta):
    """Case A: signs track z_k adiabatically. s_k = sign(z_k).
    Then sum s_k |z_k| = sum z_k = 3 (constant for normalized v0=1).
    V = 9 (constant). No delta dependence."""
    z = z_all(delta)
    # sum s_k |z_k| = sum z_k when s_k = sign(z_k) and z_k != 0
    # At z_k = 0, the contribution vanishes
    val = np.sum(z)  # should be 3.0 always
    return val**2

def V_bion_frozen_positive(delta):
    """Case B: All signs frozen at s_k = +1 (the sector before any crossing).
    V = (sum |z_k|)^2."""
    z = z_all(delta)
    return np.sum(np.abs(z))**2

def V_bion_physical(delta):
    """The physical bion potential with sign assignment tracking the
    topological sector.

    In the (-s,c,b) convention: z_0 < 0 means s_0 should flip.
    Signs: s_k(delta) = sign(z_k(delta)).
    Then s_k * sqrt(m_k) = s_k * |z_k| = z_k  (by definition).
    So V = (sum z_k)^2 = 9 = constant.  Same as Case A.

    The ONLY way to get delta-dependence is if signs DON'T track.
    So let's compute with FIXED signs. We define:

    Given a reference delta_ref where all z_k > 0, freeze signs at +1.
    As delta evolves, z_k may become negative, but s_k stays +1.
    Then s_k sqrt(m_k) = |z_k| (not z_k), and the sum changes.
    """
    return V_bion_frozen_positive(delta)

def V_bion_mixed(delta, flip_index=0):
    """Mixed case: one sign is flipped (negative) while others are positive.
    This models the (-s,c,b) sector where z_0 = -sqrt(m_s) < 0.
    s_0 = -1, s_1 = s_2 = +1.
    V = (-|z_0| + |z_1| + |z_2|)^2.
    """
    z = z_all(delta)
    signs = np.ones(3)
    signs[flip_index] = -1.0
    val = np.sum(signs * np.abs(z))
    return val**2

def V_bion_sign_tracking(delta):
    """Signs determined by the Koide charge: s_k = +1 if z_k > 0, -1 if z_k < 0.
    V = |sum s_k sqrt(m_k)|^2 = |sum z_k|^2 but with correct handling at boundaries.

    Actually: s_k sqrt(m_k) = sign(z_k) * |z_k| = z_k.  So sum = sum z_k = 3.
    V = 9 always.

    This is the ADIABATIC case.  Only nontrivial if signs DON'T track.
    """
    z = z_all(delta)
    return np.sum(z)**2  # = 9 always


# ============================================================
# The three-instanton potential for comparison
# ============================================================

def V_three_instanton(delta, phi=0):
    """V = -cos(9 delta + phi)."""
    return -np.cos(9.0 * delta + phi)


# ============================================================
# COMPUTATION
# ============================================================

N = 10000
delta_arr = np.linspace(0, TWO_PI, N, endpoint=False)

print("=" * 72)
print("BION POTENTIAL ON THE KOIDE MANIFOLD")
print("=" * 72)
print()

# ============================================================
# PART 1: Express m_k(delta) in terms of v0 and delta
# ============================================================
print("=" * 72)
print("PART 1: Koide Mass Parametrization")
print("=" * 72)
print()

print("z_k(delta) = v0 * (1 + sqrt(2) cos(delta + 2 pi k/3))   for k = 0, 1, 2")
print("m_k(delta) = z_k(delta)^2 = v0^2 * (1 + sqrt(2) cos(delta + 2 pi k/3))^2")
print()
print("Properties (all algebraic identities):")
print("  sum z_k = 3 v0       [sum cos(delta + 2pi k/3) = 0]")
print("  sum z_k^2 = 6 v0^2   [sum cos^2 = 3/2]")
print("  Q = sum(m_k)/(sum z_k)^2 = 6 v0^2 / 9 v0^2 = 2/3")
print()

# Verify at specific deltas
for label, d in [("0", 0), ("pi/4", PI/4), ("3pi/4 (seed)", DELTA_SEED), ("pi", PI)]:
    z = z_all(d)
    print(f"  delta = {label}:  z = [{z[0]:.6f}, {z[1]:.6f}, {z[2]:.6f}]")
    print(f"    sum z = {np.sum(z):.10f} (should be 3)")
    print(f"    sum z^2 = {np.sum(z**2):.10f} (should be 6)")
    print(f"    m = [{z[0]**2:.6f}, {z[1]**2:.6f}, {z[2]**2:.6f}]")
    print()


# ============================================================
# PART 2: Bion potential V_bion(delta)
# ============================================================
print("=" * 72)
print("PART 2: Bion Potential V_bion(delta)")
print("=" * 72)
print()

print("V_bion(delta) = (zeta^2/Lambda^2) |sum_k s_k(delta) sqrt(m_k(delta))|^2")
print()
print("Key identity: sqrt(m_k) = |z_k|, and s_k sqrt(m_k) = s_k |z_k|")
print()
print("CASE A (adiabatic signs, s_k = sign(z_k)):")
print("  s_k |z_k| = z_k for all k")
print("  sum s_k |z_k| = sum z_k = 3 v0")
print("  V_bion = 9 v0^2 = CONSTANT.  No delta-selection.  TRIVIAL.")
print()

# Numerical verification of Case A
V_adiab = np.array([V_bion_adiabatic(d) for d in delta_arr])
print(f"  Numerical check: V_adiabatic range = [{np.min(V_adiab):.10f}, {np.max(V_adiab):.10f}]")
print(f"  Expected: 9.0 everywhere.  Confirmed: max deviation = {np.max(np.abs(V_adiab - 9.0)):.2e}")
print()

print("CASE B (frozen signs, all s_k = +1):")
print("  V_frozen = (sum |z_k|)^2")
print("  When all z_k > 0: sum |z_k| = sum z_k = 3, V = 9.")
print("  When some z_k < 0: |z_k| != z_k, so V differs from 9.")
print("  Nontrivial structure appears only when z_k crosses zero.")
print()

V_frozen = np.array([V_bion_frozen_positive(d) for d in delta_arr])
print(f"  V_frozen range = [{np.min(V_frozen):.6f}, {np.max(V_frozen):.6f}]")
print()

# Find where z_k < 0
# z_k = 0 when cos(delta + 2pi k/3) = -1/sqrt(2)
# => delta + 2pi k/3 = 3pi/4 or 5pi/4
# For k=0: delta = 3pi/4 or 5pi/4
# For k=1: delta = 3pi/4 - 2pi/3 = pi/12 or 5pi/4 - 2pi/3 = 7pi/12
# For k=2: delta = 3pi/4 - 4pi/3 = -7pi/12 = 17pi/12 or 5pi/4 - 4pi/3 = -pi/12 = 23pi/12
print("  Zero-crossings z_k = 0 occur at:")
crossings = []
for k in range(3):
    for branch in [3*PI/4, 5*PI/4]:
        d_cross = (branch - TWO_PI_OVER_3 * k) % TWO_PI
        crossings.append((k, d_cross))
        print(f"    k={k}: delta = {d_cross:.6f} rad = {np.degrees(d_cross):.2f} deg")
print()

# Identify regions where V_frozen != 9
print("  Regions where V_frozen differs from 9 (some z_k < 0):")
z0_neg = np.array([z_k(d, 0) < 0 for d in delta_arr])
z1_neg = np.array([z_k(d, 1) < 0 for d in delta_arr])
z2_neg = np.array([z_k(d, 2) < 0 for d in delta_arr])
any_neg = z0_neg | z1_neg | z2_neg
frac_neg = np.mean(any_neg)
print(f"    Fraction of delta range with some z_k < 0: {frac_neg:.4f}")
print(f"    (= {frac_neg * 360:.1f} degrees out of 360)")
print()

# Find minima and maxima of V_frozen
minima_idx = argrelextrema(V_frozen, np.less, order=50)[0]
maxima_idx = argrelextrema(V_frozen, np.greater, order=50)[0]

print("  Local minima of V_frozen:")
for idx in minima_idx:
    d = delta_arr[idx]
    z = z_all(d)
    print(f"    delta = {np.degrees(d):.3f} deg, V = {V_frozen[idx]:.6f}, "
          f"z = [{z[0]:.4f}, {z[1]:.4f}, {z[2]:.4f}]")
print()

print("  Local maxima of V_frozen:")
for idx in maxima_idx:
    d = delta_arr[idx]
    z = z_all(d)
    print(f"    delta = {np.degrees(d):.3f} deg, V = {V_frozen[idx]:.6f}, "
          f"z = [{z[0]:.4f}, {z[1]:.4f}, {z[2]:.4f}]")
print()


# ============================================================
# PART 2b: Case C - Mixed signs (one flipped)
# ============================================================
print("CASE C (one sign flipped, s_0 = -1, s_1 = s_2 = +1):")
print("  V_mixed = (-|z_0| + |z_1| + |z_2|)^2")
print("  This models a topological sector where the first monopole has opposite sign.")
print()

V_mixed0 = np.array([V_bion_mixed(d, flip_index=0) for d in delta_arr])
V_mixed1 = np.array([V_bion_mixed(d, flip_index=1) for d in delta_arr])
V_mixed2 = np.array([V_bion_mixed(d, flip_index=2) for d in delta_arr])

for flip, V_m, label in [(0, V_mixed0, "s_0=-1"), (1, V_mixed1, "s_1=-1"), (2, V_mixed2, "s_2=-1")]:
    print(f"  {label}: V range = [{np.min(V_m):.6f}, {np.max(V_m):.6f}]")
    minima_m = argrelextrema(V_m, np.less, order=50)[0]
    if len(minima_m) > 0:
        print(f"    Minima:")
        for idx in minima_m:
            d = delta_arr[idx]
            z = z_all(d)
            signs = np.ones(3)
            signs[flip] = -1
            val = np.sum(signs * np.abs(z))
            print(f"      delta = {np.degrees(d):.3f} deg, V = {V_m[idx]:.6f}, "
                  f"sum(s*|z|) = {val:.4f}")
    # Find zeros: V = 0 when sum s_k |z_k| = 0
    zeros = np.where(V_m < 0.01)[0]
    if len(zeros) > 0:
        # Report first and last
        print(f"    V ~ 0 near delta = {np.degrees(delta_arr[zeros[0]]):.2f} - "
              f"{np.degrees(delta_arr[zeros[-1]]):.2f} deg")
    print()


# ============================================================
# PART 3: Map V_bion(delta) for delta in [0, 2pi]
# ============================================================
print("=" * 72)
print("PART 3: Full Map of V_bion(delta)")
print("=" * 72)
print()

print("Computing all bion potential variants over [0, 2pi]...")
print()

# The physically relevant potential: with sign assignment from the problem statement
# s_k = +1 if z_k > 0, s_k = -1 if z_k < 0
# This gives V = (sum z_k)^2 = 9 = constant.
# So the "physical" bion potential defined in the problem is TRIVIAL.

# However, the problem statement says "s_k flips when m_k passes through zero".
# This implies that the signs DO track z_k, making V constant.
# But the physical setup of the bion Kahler potential has a different structure.

# Let me implement the problem's EXACT definition:
# V_bion(delta) = |s_0 sqrt(m_0) + s_1 sqrt(m_1) + s_2 sqrt(m_2)|^2
# where s_k = sign(z_k(delta)) and sqrt(m_k) = |z_k|.
# Then s_k sqrt(m_k) = sign(z_k) |z_k| = z_k, and sum = 3, V = 9.

# Unless... the sign assignment is TOPOLOGICAL, not local.
# In the bion picture, s_k labels the monopole-instanton species.
# The sign is determined by which root alpha_k the monopole sits on.
# It's NOT determined by whether z_k > 0 or < 0.

# The correct interpretation: the bion amplitude involves
# sum_k s_k sqrt(m_k) where s_k are FIXED (topological) signs,
# NOT tracking z_k. The three sectors correspond to:
# Sector (+++): s = (+1, +1, +1) -> V = (sum |z_k|)^2
# Sector (-++): s = (-1, +1, +1) -> V = (-|z_0| + |z_1| + |z_2|)^2
# Sector (+-+): etc.
# Sector (++-): etc.

# The total bion potential is a sum over all sectors:
# V_total = sum_{sectors} A_sector * |sum_k s_k sqrt(m_k)|^2

# For SU(3) with the extended Dynkin triangle, there are 3 adjacent
# monopole pairs (bion types). The signs come from the relative
# orientation of the monopole and anti-monopole.

# The simplest nontrivial case: V = (sum |z_k|)^2
# This is the frozen-sign (all +1) potential.

print("Key insight: The sign assignment s_k is TOPOLOGICAL (fixed per sector),")
print("not adiabatic (tracking z_k). So the nontrivial bion potential is:")
print()
print("  V_bion(delta) = (sum_k |z_k(delta)|)^2     [frozen signs, all s_k = +1]")
print()
print("This has nontrivial delta-dependence whenever some z_k < 0.")
print()

# Find all local minima with high resolution
print("Local minima of V_frozen = (sum |z_k|)^2:")
# Use finer grid
N_fine = 100000
d_fine = np.linspace(0, TWO_PI, N_fine, endpoint=False)
V_fine = np.array([V_bion_frozen_positive(d) for d in d_fine])
min_idx_fine = argrelextrema(V_fine, np.less, order=500)[0]

all_minima = []
for idx in min_idx_fine:
    d = d_fine[idx]
    # Refine with scipy
    res = minimize_scalar(V_bion_frozen_positive, bounds=(d - 0.01, d + 0.01), method='bounded')
    d_opt = res.x
    V_opt = res.fun
    z = z_all(d_opt)
    all_minima.append((d_opt, V_opt, z.copy()))
    print(f"  delta = {np.degrees(d_opt):.6f} deg = {d_opt:.8f} rad")
    print(f"    V = {V_opt:.8f}")
    print(f"    z = [{z[0]:.6f}, {z[1]:.6f}, {z[2]:.6f}]")
    print(f"    sum|z| = {np.sum(np.abs(z)):.8f}")
    print()

# Where V_frozen is NOT flat (= 9):
print("V_frozen profile at selected points:")
print(f"{'delta (deg)':>12} {'V_frozen':>12} {'z_0':>10} {'z_1':>10} {'z_2':>10} {'note':>15}")
print("-" * 75)
for d_deg in np.arange(0, 361, 15):
    d = np.radians(d_deg)
    z = z_all(d)
    V = V_bion_frozen_positive(d)
    neg = "neg" if any(z < 0) else ""
    print(f"{d_deg:>12.1f} {V:>12.6f} {z[0]:>10.4f} {z[1]:>10.4f} {z[2]:>10.4f} {neg:>15}")
print()


# ============================================================
# PART 4: Quark triple (-s, c, b)
# ============================================================
print("=" * 72)
print("PART 4: Quark Triple (-s, c, b)")
print("=" * 72)
print()

sr_quarks = [-np.sqrt(m_s), np.sqrt(m_c), np.sqrt(m_b)]
v0_quarks, delta_quarks = extract_v0_delta(sr_quarks)
delta_quarks_pos = delta_quarks % TWO_PI

print(f"Signed roots: z = [{sr_quarks[0]:.6f}, {sr_quarks[1]:.6f}, {sr_quarks[2]:.6f}]")
print(f"v0 = sum(z)/3 = {v0_quarks:.6f} MeV^(1/2)")
print(f"delta = {delta_quarks:.8f} rad = {np.degrees(delta_quarks):.6f} deg")
print(f"delta mod 2pi = {delta_quarks_pos:.8f} rad = {np.degrees(delta_quarks_pos):.6f} deg")
print()

# Verify reconstruction
print("Reconstruction check:")
for k in range(3):
    zk_pred = v0_quarks * (1 + SQRT2 * np.cos(delta_quarks + TWO_PI_OVER_3 * k))
    print(f"  k={k}: z_pred = {zk_pred:.6f}, z_actual = {sr_quarks[k]:.6f}, "
          f"diff = {zk_pred - sr_quarks[k]:.3e}")
print()

# Koide Q check
Q_quarks_unsigned = koide_Q([m_s, m_c, m_b])
Q_quarks_signed = koide_Q([m_s, m_c, m_b], signs=[-1, 1, 1])
print(f"Koide Q(-s,c,b) unsigned roots = {Q_quarks_unsigned:.8f}  (WRONG: uses |sqrt(m_s)| not -sqrt(m_s))")
print(f"Koide Q(-s,c,b) signed roots   = {Q_quarks_signed:.8f}  (2/3 = {2/3:.8f}, dev = {(Q_quarks_signed - 2/3)/(2/3)*100:.4f}%)")
print(f"  Q with signed roots: sum(m) = {m_s+m_c+m_b:.1f}, (sum sigma)^2 = {(-np.sqrt(m_s)+np.sqrt(m_c)+np.sqrt(m_b))**2:.1f}")
print()

# v0 doubling
sr_seed = [np.sqrt(m_s), np.sqrt(m_c)]
v0_seed = (np.sqrt(m_s) + np.sqrt(m_c)) / 3.0
v0_full = (-np.sqrt(m_s) + np.sqrt(m_c) + np.sqrt(m_b)) / 3.0
ratio_v0 = v0_full / v0_seed
print(f"v0(seed) = (sqrt(m_s) + sqrt(m_c))/3 = {v0_seed:.6f} MeV^(1/2)")
print(f"v0(full) = (-sqrt(m_s) + sqrt(m_c) + sqrt(m_b))/3 = {v0_full:.6f} MeV^(1/2)")
print(f"v0(full)/v0(seed) = {ratio_v0:.6f}  (target: 2.0005)")
print()

# Check where on the V_bion landscape the physical delta sits
# We need to use normalized delta (the delta from the Koide parametrization)
# The normalized z_k = (1 + sqrt(2) cos(delta + 2pi k/3))
# At the physical quark delta:
z_phys_q = z_all(delta_quarks_pos)
V_phys_q_frozen = V_bion_frozen_positive(delta_quarks_pos)
V_phys_q_mixed0 = V_bion_mixed(delta_quarks_pos, flip_index=0)

print(f"At physical quark delta = {np.degrees(delta_quarks_pos):.4f} deg:")
print(f"  Normalized z = [{z_phys_q[0]:.6f}, {z_phys_q[1]:.6f}, {z_phys_q[2]:.6f}]")
print(f"  All z > 0? {all(z_phys_q > 0)}")
print(f"  V_frozen(+++) = {V_phys_q_frozen:.8f}")
print(f"  V_mixed(s_0=-1) = {V_phys_q_mixed0:.8f}")
print()

# Check: is the physical delta at or near a minimum?
# For the frozen-sign potential:
eps = 1e-4
V_left = V_bion_frozen_positive(delta_quarks_pos - eps)
V_right = V_bion_frozen_positive(delta_quarks_pos + eps)
V_center = V_bion_frozen_positive(delta_quarks_pos)
is_min_frozen = V_center < V_left and V_center < V_right
print(f"  V_frozen at delta +/- epsilon:")
print(f"    V(delta-eps) = {V_left:.10f}")
print(f"    V(delta)     = {V_center:.10f}")
print(f"    V(delta+eps) = {V_right:.10f}")
print(f"    Local minimum? {is_min_frozen}")
print()

# For mixed sign potential:
V_left_m = V_bion_mixed(delta_quarks_pos - eps, 0)
V_right_m = V_bion_mixed(delta_quarks_pos + eps, 0)
V_center_m = V_bion_mixed(delta_quarks_pos, 0)
is_min_mixed = V_center_m < V_left_m and V_center_m < V_right_m
print(f"  V_mixed(s_0=-1) at delta +/- epsilon:")
print(f"    V(delta-eps) = {V_left_m:.10f}")
print(f"    V(delta)     = {V_center_m:.10f}")
print(f"    V(delta+eps) = {V_right_m:.10f}")
print(f"    Local minimum? {is_min_mixed}")
print()

# The quark delta is in the all-positive regime (check):
# For the normalized parametrization, all z_k > 0 when delta is not near any crossing.
# The crossings for z_0 = 0 are at delta = 3pi/4 and 5pi/4.
# The physical quark delta ~ 2.74 rad ~ 157 deg is between 3pi/4 = 135 deg and 5pi/4 = 225 deg.
# So z_0 < 0 at the physical quark delta!

print(f"Physical quark delta = {np.degrees(delta_quarks_pos):.2f} deg")
print(f"z_0 crossing at 3pi/4 = 135 deg and 5pi/4 = 225 deg")
print(f"z_0 < 0 in [135, 225] deg => z_0 at physical delta is {z_phys_q[0]:.6f}")
if z_phys_q[0] < 0:
    print(f"  CONFIRMED: z_0 < 0 at physical quark delta.")
    print(f"  This means the frozen-sign potential V = (sum |z_k|)^2 differs from 9.")
else:
    print(f"  z_0 > 0: frozen-sign potential equals 9 (flat).")
print()


# ============================================================
# PART 5: Lepton triple (e, mu, tau)
# ============================================================
print("=" * 72)
print("PART 5: Lepton Triple (e, mu, tau)")
print("=" * 72)
print()

sr_leptons = [np.sqrt(m_e), np.sqrt(m_mu), np.sqrt(m_tau)]
v0_leptons, delta_leptons = extract_v0_delta(sr_leptons)
delta_leptons_pos = delta_leptons % TWO_PI

print(f"Signed roots: z = [{sr_leptons[0]:.6f}, {sr_leptons[1]:.6f}, {sr_leptons[2]:.6f}]")
print(f"v0 = sum(z)/3 = {v0_leptons:.6f} MeV^(1/2)")
print(f"delta = {delta_leptons:.8f} rad = {np.degrees(delta_leptons):.6f} deg")
print(f"delta mod 2pi = {delta_leptons_pos:.8f} rad = {np.degrees(delta_leptons_pos):.6f} deg")
print()

# Check z_k signs at lepton delta
z_phys_l = z_all(delta_leptons_pos)
print(f"At physical lepton delta = {np.degrees(delta_leptons_pos):.4f} deg:")
print(f"  Normalized z = [{z_phys_l[0]:.6f}, {z_phys_l[1]:.6f}, {z_phys_l[2]:.6f}]")
print(f"  All z > 0? {all(z_phys_l > 0)}")
print()

V_phys_l_frozen = V_bion_frozen_positive(delta_leptons_pos)
print(f"  V_frozen(+++) = {V_phys_l_frozen:.8f}")

if all(z_phys_l > 0):
    print(f"  All z > 0 => V_frozen = (sum z_k)^2 = 9.  No information from bion potential.")
else:
    eps = 1e-4
    V_left_l = V_bion_frozen_positive(delta_leptons_pos - eps)
    V_right_l = V_bion_frozen_positive(delta_leptons_pos + eps)
    is_min_l = V_phys_l_frozen < V_left_l and V_phys_l_frozen < V_right_l
    print(f"  Local minimum? {is_min_l}")
print()

# delta mod (2pi/3) vs 2/9
delta_lep_mod = delta_leptons_pos % TWO_PI_OVER_3
print(f"delta mod (2pi/3) = {delta_lep_mod:.10f} rad")
print(f"2/9               = {2.0/9.0:.10f} rad")
print(f"Residual          = {delta_lep_mod - 2.0/9.0:.6e} rad")
print()


# ============================================================
# PART 6: Key Question - What drives the bloom?
# ============================================================
print("=" * 72)
print("PART 6: What Drives the Bloom?")
print("=" * 72)
print()

print("RESULT: The bion potential V_bion(delta) = |sum s_k sqrt(m_k)|^2")
print("on the Koide manifold is:")
print()
print("  (a) If signs track z_k adiabatically (s_k = sign z_k):")
print("      V = (sum z_k)^2 = 9 v0^2 = CONSTANT. No delta-selection.")
print()
print("  (b) If signs are frozen (topological, all s_k = +1):")
print("      V = (sum |z_k|)^2. Nontrivial ONLY in the region where")
print("      some z_k < 0 (i.e., near/past the seed at delta = 3pi/4).")
print()

# Analyze the frozen-sign potential near the seed
print("Behavior of V_frozen near the seed (delta = 3pi/4 = 135 deg):")
print()
d_near_seed = np.linspace(DELTA_SEED - 0.5, DELTA_SEED + 0.5, 1000)
V_near_seed = np.array([V_bion_frozen_positive(d) for d in d_near_seed])

print(f"  V_frozen at seed:              {V_bion_frozen_positive(DELTA_SEED):.8f}")
print(f"  V_frozen just before seed:     {V_bion_frozen_positive(DELTA_SEED - 0.01):.8f}")
print(f"  V_frozen just after seed:      {V_bion_frozen_positive(DELTA_SEED + 0.01):.8f}")
print()

# At the seed, z_0 = 0, z_1, z_2 > 0.  So sum |z_k| = sum z_k = 3, V = 9.
# Just past the seed (delta > 3pi/4), z_0 < 0. |z_0| = -z_0.
# sum |z_k| = |z_0| + z_1 + z_2 = -z_0 + z_1 + z_2 = (z_1+z_2-z_0) = 3 - 2z_0 > 3.
# So V > 9 just past the seed. The seed is a LOCAL MINIMUM of V_frozen.

# Analytical: V_frozen = (|z_0| + z_1 + z_2)^2 when z_0 < 0 and z_1,z_2 > 0
# = (-z_0 + z_1 + z_2)^2 = (3 - 2z_0)^2  [since z_0+z_1+z_2 = 3]
# = 9 - 12 z_0 + 4 z_0^2
# Now z_0 = 1 + sqrt(2) cos(delta) which is < 0 when delta > 3pi/4.
# dV/ddelta = (-12 + 8z_0) dz_0/ddelta = (-12 + 8z_0)(-sqrt(2) sin delta)
# At the seed, z_0 = 0: dV/ddelta = (-12)(-sqrt(2) sin(3pi/4)) = -12*(-1) = +12.
# Wait, sin(3pi/4) = 1/sqrt(2), so -sqrt(2)*sin(3pi/4) = -1.
# dV/ddelta = (-12)(-1) = 12 > 0.
# So V increases as delta increases past the seed. But BEFORE the seed (delta < 3pi/4),
# z_0 > 0 and V = 9 (flat). So the seed is a CORNER (not a smooth minimum).

print("Analytical behavior near the seed (delta = 3pi/4):")
print()
print("  For delta < 3pi/4 (z_0 > 0): V_frozen = 9 (constant)")
print("  For delta > 3pi/4 (z_0 < 0): V_frozen = (3 - 2z_0)^2 > 9")
print("      = 9 - 12 z_0 + 4 z_0^2  with z_0 = 1 + sqrt(2) cos(delta) < 0")
print()
print("  The seed delta = 3pi/4 is a CORNER minimum of V_frozen:")
print("  it's flat on the left (V=9) and rises on the right (V > 9).")
print("  NOT a smooth minimum -- the first derivative is discontinuous.")
print()

# Now check: does the physical quark delta sit in the rising part?
print(f"Physical quark delta = {np.degrees(delta_quarks_pos):.2f} deg > seed = 135 deg?")
if delta_quarks_pos > DELTA_SEED:
    print(f"  YES: delta > seed. V_frozen is ABOVE 9 at the physical quark point.")
    print(f"  V_frozen(physical) = {V_phys_q_frozen:.6f} > 9")
    print(f"  The bion potential PUSHES delta back toward the seed (135 deg).")
    print(f"  It does NOT favor the bloom direction.")
else:
    print(f"  NO: delta < seed. V_frozen = 9 (flat). No force from bion potential.")
print()

# Now consider the three-instanton potential
print("Three-instanton potential V_3inst = -cos(9 delta):")
print()
# Nearest cos(9delta) minimum to seed and to physical delta
cos9_at_seed = np.cos(9 * DELTA_SEED)
cos9_at_phys_q = np.cos(9 * delta_quarks_pos)
print(f"  cos(9 * delta_seed) = cos(27pi/4) = cos(3pi/4) = {cos9_at_seed:.6f}")
print(f"  cos(9 * delta_phys_quark) = {cos9_at_phys_q:.6f}")
print()

# Minima of -cos(9delta) at 9*delta = 2pi*n, i.e. delta = 2pi*n/9
print(f"  cos(9 delta) minima (V=-1) at delta = 2pi n/9:")
for n in range(9):
    d_min_3 = TWO_PI * n / 9
    dist_seed = np.degrees(d_min_3 - DELTA_SEED)
    dist_phys = np.degrees(d_min_3 - delta_quarks_pos)
    print(f"    n={n}: delta = {np.degrees(d_min_3):.2f} deg, "
          f"dist from seed = {dist_seed:+.2f} deg, "
          f"dist from phys = {dist_phys:+.2f} deg")
print()

# Combined analysis: which minimum of -cos(9delta) is nearest to the physical quark delta?
nearest_n_q = min(range(9), key=lambda n: abs(((TWO_PI*n/9 - delta_quarks_pos + PI) % TWO_PI) - PI))
d_nearest_q = TWO_PI * nearest_n_q / 9
dist_q = ((d_nearest_q - delta_quarks_pos + PI) % TWO_PI) - PI
print(f"Nearest -cos(9 delta) minimum to physical quark delta:")
print(f"  n={nearest_n_q}, delta = {np.degrees(d_nearest_q):.2f} deg, "
      f"distance = {np.degrees(dist_q):.2f} deg")
print()

# Nearest to lepton delta
nearest_n_l = min(range(9), key=lambda n: abs(((TWO_PI*n/9 - delta_leptons_pos + PI) % TWO_PI) - PI))
d_nearest_l = TWO_PI * nearest_n_l / 9
dist_l = ((d_nearest_l - delta_leptons_pos + PI) % TWO_PI) - PI
print(f"Nearest -cos(9 delta) minimum to physical lepton delta:")
print(f"  n={nearest_n_l}, delta = {np.degrees(d_nearest_l):.2f} deg, "
      f"distance = {np.degrees(dist_l):.2f} deg")
print()


# ============================================================
# PART 6b: Comprehensive potential landscape
# ============================================================
print("=" * 72)
print("PART 6b: Combined Bion + Three-Instanton Landscape")
print("=" * 72)
print()

# V_total(delta) = A * V_frozen(delta) + B * (-cos(9 delta))
# where A, B are relative couplings to be determined.

# The bion potential V_frozen = 9 in the z_k > 0 regime and > 9 otherwise.
# Define the excess: Delta_V = V_frozen - 9, which is >= 0 everywhere.
# V_total = 9A + A * Delta_V - B cos(9 delta)

# The competition is: the bion excess pushes toward the seed (where Delta_V = 0),
# while the instanton term pushes toward the nearest cos(9delta) minimum.

Delta_V = np.array([V_bion_frozen_positive(d) - 9.0 for d in delta_arr])
cos9 = np.array([-np.cos(9.0 * d) for d in delta_arr])

# Plot the two contributions
print("Profile along the bloom direction (delta around 135-160 deg):")
print(f"{'delta (deg)':>12} {'Delta_V':>12} {'-cos(9d)':>12}")
print("-" * 40)
for d_deg in np.arange(120, 175, 5):
    d = np.radians(d_deg)
    dv = V_bion_frozen_positive(d) - 9.0
    c9 = -np.cos(9.0 * d)
    print(f"{d_deg:>12.1f} {dv:>12.6f} {c9:>12.6f}")
print()

# For the combined potential: need A/B ratio.
# At the physical quark delta, dV/ddelta = 0 requires:
# A * dDelta_V/ddelta = B * 9 * sin(9*delta)
# At physical quark delta:
d_phys = delta_quarks_pos
eps = 1e-6
dDeltaV = (V_bion_frozen_positive(d_phys + eps) - V_bion_frozen_positive(d_phys - eps)) / (2*eps)
sin9_phys_q = 9.0 * np.sin(9.0 * d_phys)
if abs(sin9_phys_q) > 1e-10 and abs(dDeltaV) > 1e-10:
    BA_ratio = dDeltaV / sin9_phys_q
    print(f"Required B/A for minimum at physical quark delta:")
    print(f"  dDelta_V/ddelta = {dDeltaV:.6f}")
    print(f"  9 sin(9*delta) = {sin9_phys_q:.6f}")
    print(f"  B/A = {BA_ratio:.6f}")
    print()

    # Check second derivative (stability)
    d2DeltaV = (V_bion_frozen_positive(d_phys + eps) - 2*V_bion_frozen_positive(d_phys)
                + V_bion_frozen_positive(d_phys - eps)) / eps**2
    d2cos9 = -81.0 * np.cos(9.0 * d_phys)
    d2V_total = d2DeltaV + BA_ratio * d2cos9
    print(f"  d^2 V_bion / ddelta^2 = {d2DeltaV:.4f}")
    print(f"  d^2 (-cos9d) / ddelta^2 = {d2cos9:.4f}")
    print(f"  d^2 V_total / ddelta^2 = {d2DeltaV + BA_ratio * d2cos9:.4f}")
    print(f"  Stable minimum? {d2V_total > 0}")
    print()
else:
    print(f"  dDelta_V/ddelta = {dDeltaV:.6e}")
    print(f"  9 sin(9*delta) = {sin9_phys_q:.6e}")
    if abs(dDeltaV) < 1e-10:
        print(f"  Delta_V is flat at physical delta => already at bion minimum")
    print()


# ============================================================
# FINAL SUMMARY
# ============================================================
print("=" * 72)
print("FINAL SUMMARY AND ANSWER TO KEY QUESTION")
print("=" * 72)
print()

print("1. MASS PARAMETRIZATION:")
print(f"   m_k(delta) = v0^2 * (1 + sqrt(2) cos(delta + 2pi k/3))^2")
print(f"   with v0 = (sum z_k)/3, Q = 2/3 identically.")
print()

print("2. BION POTENTIAL ON KOIDE MANIFOLD:")
print(f"   V_bion = |sum s_k sqrt(m_k)|^2")
print(f"   With ADIABATIC signs (s_k = sign(z_k)): V = 9 v0^2 = CONSTANT.")
print(f"   The bion potential with adiabatic sign-tracking is FLAT on the Koide")
print(f"   manifold. It cannot select delta. THIS IS AN EXACT ALGEBRAIC RESULT.")
print()
print(f"   With FROZEN signs (all s_k = +1): V = (sum |z_k|)^2.")
print(f"   Equals 9 when all z_k > 0; rises when any z_k goes negative.")
print(f"   Has a corner minimum at the seed (delta = 3pi/4 = 135 deg).")
print()

print("3. LOCAL MINIMA OF V_bion (frozen signs):")
if len(all_minima) > 0:
    for d_opt, V_opt, z_opt in all_minima:
        print(f"   delta = {np.degrees(d_opt):.3f} deg, V = {V_opt:.6f}")
else:
    print("   The global minimum V = 9 is attained across the entire region where")
    print("   all z_k > 0. This region spans ~240 degrees out of 360 (two-thirds).")
    print("   The 'minimum' is a flat plateau, not a point.")
print()

print("4. QUARK TRIPLE (-s, c, b):")
print(f"   Physical delta = {np.degrees(delta_quarks_pos):.4f} deg")
print(f"   v0(seed) = {v0_seed:.6f} MeV^(1/2)")
print(f"   v0(full) = {v0_full:.6f} MeV^(1/2)")
print(f"   v0(full)/v0(seed) = {ratio_v0:.6f}  (confirming 2.0005)")
print()
if z_phys_q[0] < 0:
    print(f"   z_0 < 0 at physical delta => V_frozen = {V_phys_q_frozen:.4f} > 9.")
    print(f"   The bion potential DISFAVORS the bloom position!")
    print(f"   It pushes delta BACK toward the seed at 135 deg.")
else:
    print(f"   All z_k > 0 => V_frozen = 9. Bion potential gives no information.")
print()

print("5. LEPTON TRIPLE (e, mu, tau):")
print(f"   Physical delta = {np.degrees(delta_leptons_pos):.4f} deg")
if all(z_phys_l > 0):
    print(f"   All z_k > 0 => V_frozen = 9. Bion potential is FLAT here.")
    print(f"   No delta-selection from the bion potential for leptons.")
else:
    print(f"   Some z_k < 0 => V_frozen != 9.")
print()

print("6. KEY ANSWER:")
print()
print("   The bion potential V_bion(delta) does NOT have a minimum at the physical")
print("   quark bloom position. In fact:")
print()
print("   (a) With adiabatic signs: V = const everywhere. No delta-selection at all.")
print()
print("   (b) With frozen signs: V has its minimum at/before the seed (135 deg)")
print("       and RISES at the physical quark delta (~157 deg). The bion potential")
print("       acts as a RESTORING force toward the seed, opposing the bloom.")
print()
print("   The bloom must therefore be driven by a DIFFERENT mechanism. The")
print("   three-instanton potential V ~ -cos(9 delta) is the natural candidate:")
print(f"   its nearest minimum to the physical quark delta is at")
print(f"   delta = {np.degrees(d_nearest_q):.1f} deg (distance {abs(np.degrees(dist_q)):.1f} deg).")
print()
print("   The physical picture is: the three-instanton potential provides the")
print("   attractive force that drives delta away from the seed toward a cos(9delta)")
print("   minimum, while the bion potential provides a restoring force that keeps")
print("   delta from going too far. The equilibrium is the physical bloom position.")
print()
print("   This is consistent with the existing result that the bloom must be")
print("   NONPERTURBATIVE -- it is driven by multi-instanton effects (the")
print("   cos(9 delta) landscape) rather than by the semiclassical bion potential.")

print()
print("Done.")
