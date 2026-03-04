#!/usr/bin/env python3
"""
Scan all C(9,3)=84 triplets of SM fermion masses for near-Koide relations.

Q(m1,m2,m3) = (sqrt(m1) + sqrt(m2) + sqrt(m3))^2 / (m1 + m2 + m3)

Q = 3/2 is the Koide relation (theta = 45 degrees).
"""

import numpy as np
from itertools import combinations
from math import comb

# ── PDG 2024 masses, all in MeV ──────────────────────────────────────────
particles = {
    'e':  0.51099895,        # MeV
    'mu': 105.6583755,       # MeV
    'tau': 1776.86,          # MeV
    'u':  2.16,              # MeV
    'd':  4.67,              # MeV
    's':  93.4,              # MeV
    'c':  1270.0,            # MeV  (1.27 GeV)
    'b':  4180.0,            # MeV  (4.18 GeV)
    't':  172690.0,          # MeV  (172.69 GeV)
}

names = list(particles.keys())
masses = np.array([particles[n] for n in names])

# ── Sign combinations: 8 total (+++, --+, -+-, +--, -++, +-+, ++-, ---) ──
# Represent as tuples of signs for (m1, m2, m3)
sign_combos = [
    (+1, +1, +1),
    (-1, +1, +1),
    (+1, -1, +1),
    (+1, +1, -1),
    (-1, -1, +1),
    (-1, +1, -1),
    (+1, -1, -1),
    (-1, -1, -1),
]

def sign_label(signs):
    return ''.join('+' if s > 0 else '-' for s in signs)

def Q_direct(m1, m2, m3, s1=1, s2=1, s3=1):
    """Compute Q via the direct formula."""
    numerator = (s1*np.sqrt(m1) + s2*np.sqrt(m2) + s3*np.sqrt(m3))**2
    denominator = m1 + m2 + m3
    return numerator / denominator

def Q_angle(m1, m2, m3, s1=1, s2=1, s3=1):
    """Compute Q via the angle formula: Q = 3 cos^2(theta)."""
    S = m1 + m2 + m3
    cos_theta = (s1*np.sqrt(m1) + s2*np.sqrt(m2) + s3*np.sqrt(m3)) / np.sqrt(3*S)
    theta = np.arccos(np.clip(cos_theta, -1, 1))
    return 3 * cos_theta**2, np.degrees(theta)

# ── Verify C(9,3) = 84 ──────────────────────────────────────────────────
triplets = list(combinations(range(len(names)), 3))
assert len(triplets) == comb(9, 3) == 84, f"Expected 84 triplets, got {len(triplets)}"

# ── Main scan ────────────────────────────────────────────────────────────
results = []       # (deviation, triplet_names, signs_label, Q_val, theta_deg)
max_discrepancy = 0.0

for idx_triple in triplets:
    i, j, k = idx_triple
    m1, m2, m3 = masses[i], masses[j], masses[k]
    trip_name = f"({names[i]}, {names[j]}, {names[k]})"

    for signs in sign_combos:
        s1, s2, s3 = signs
        Q1 = Q_direct(m1, m2, m3, s1, s2, s3)
        Q2, theta = Q_angle(m1, m2, m3, s1, s2, s3)

        disc = abs(Q1 - Q2)
        if disc > max_discrepancy:
            max_discrepancy = disc

        if disc > 1e-10:
            print(f"WARNING: discrepancy {disc:.2e} for {trip_name} signs={sign_label(signs)}")

        dev = abs(Q1 - 1.5)
        results.append((dev, trip_name, sign_label(signs), Q1, theta))

# Sort by deviation
results.sort(key=lambda x: x[0])

# ── Output: ranked table of |Q - 3/2| < 0.05 ───────────────────────────
print("=" * 90)
print(f"{'Rank':>4}  {'Triplet':<22} {'Signs':>5}  {'Q':>14}  {'|Q - 3/2|':>14}  {'theta (deg)':>12}")
print("=" * 90)

rank = 0
for dev, trip_name, sl, Q_val, theta in results:
    if dev >= 0.05:
        break
    rank += 1
    print(f"{rank:4d}  {trip_name:<22} {sl:>5}  {Q_val:14.10f}  {dev:14.10e}  {theta:12.6f}")

print("=" * 90)
print(f"\nTotal entries with |Q - 3/2| < 0.05: {rank}")
print(f"Maximum discrepancy between direct and angle formulas: {max_discrepancy:.2e}")

# ── Charged lepton triplet to high precision ─────────────────────────────
print("\n" + "=" * 90)
print("CHARGED LEPTON TRIPLET (e, mu, tau) — HIGH PRECISION")
print("=" * 90)

me  = 0.51099895       # MeV
mmu = 105.6583755      # MeV
mtau = 1776.86         # MeV

Q_lep = Q_direct(me, mmu, mtau)
Q_lep_angle, theta_lep = Q_angle(me, mmu, mtau)

print(f"  Q (direct)      = {Q_lep:.10f}")
print(f"  Q (angle)       = {Q_lep_angle:.10f}")
print(f"  |Q - 3/2|       = {abs(Q_lep - 1.5):.10e}")
print(f"  theta            = {theta_lep:.8f} degrees")
print(f"  |theta - 45|     = {abs(theta_lep - 45.0):.8e} degrees")
print(f"  Discrepancy      = {abs(Q_lep - Q_lep_angle):.2e}")

# ── Extended precision with mpmath ───────────────────────────────────────
print("\n--- Using mpmath for extended precision ---")
try:
    from mpmath import mp, mpf, sqrt, acos, pi, cos
    mp.dps = 50  # 50 decimal places

    me_mp  = mpf('0.51099895')
    mmu_mp = mpf('105.6583755')
    mtau_mp = mpf('1776.86')

    S_mp = me_mp + mmu_mp + mtau_mp
    sqsum = sqrt(me_mp) + sqrt(mmu_mp) + sqrt(mtau_mp)
    Q_mp = sqsum**2 / S_mp
    cos_theta_mp = sqsum / sqrt(3 * S_mp)
    theta_mp = acos(cos_theta_mp) * 180 / pi

    print(f"  Q               = {mp.nstr(Q_mp, 15)}")
    print(f"  Q (10 sig figs) = {mp.nstr(Q_mp, 10)}")
    print(f"  |Q - 3/2|       = {mp.nstr(abs(Q_mp - mpf('1.5')), 10)}")
    print(f"  theta            = {mp.nstr(theta_mp, 12)} degrees")
    print(f"  theta (8 s.f.)   = {mp.nstr(theta_mp, 8)} degrees")
    print(f"  |theta - 45|     = {mp.nstr(abs(theta_mp - 45), 10)} degrees")
except ImportError:
    print("  mpmath not available; skipping extended precision.")

# ── Summary statistics ───────────────────────────────────────────────────
print("\n" + "=" * 90)
print("SUMMARY STATISTICS")
print("=" * 90)
total_combos = 84 * 8
print(f"Total triplet-sign combinations scanned: {total_combos}")
devs = [r[0] for r in results]
print(f"Closest to Q = 3/2: {results[0][1]} signs={results[0][2]}, "
      f"|Q-3/2| = {results[0][0]:.6e}")
print(f"Number with |Q - 3/2| < 0.001: {sum(1 for d in devs if d < 0.001)}")
print(f"Number with |Q - 3/2| < 0.01:  {sum(1 for d in devs if d < 0.01)}")
print(f"Number with |Q - 3/2| < 0.05:  {sum(1 for d in devs if d < 0.05)}")
