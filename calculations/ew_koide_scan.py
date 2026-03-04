#!/usr/bin/env python3
"""
Exhaustive Koide ratio scan over SM masses + EW scales.

The correct Koide formula:
  Q = (m1 + m2 + m3) / (sqrt(m1) + sqrt(m2) + sqrt(m3))^2
  Q = 2/3 for (e, mu, tau) to ~6e-6.

Signed variant (negative sqrt on smallest mass):
  Q_signed = (m1 + m2 + m3) / (-sqrt(m_min) + sqrt(m_mid) + sqrt(m_max))^2

Pair version:
  Q_pair = (m1 + m2) / (sqrt(m1) + sqrt(m2))^2
"""

import numpy as np
from itertools import combinations

# ============================================================
# Mass list (GeV)
# ============================================================
masses = {
    'm_e':    0.000510999,
    'm_mu':   0.1056584,
    'm_tau':  1.77686,
    'm_u':    0.00216,
    'm_d':    0.00467,
    'm_s':    0.0934,
    'm_c':    1.270,
    'm_b':    4.180,
    'm_t':    172.76,
    'm_W':    80.369,
    'm_Z':    91.1876,
    'm_H':    125.25,
    'f_pi':   0.092,
    'v_EW':   246.22,
}

names = list(masses.keys())
vals = [masses[n] for n in names]
N = len(names)

def Q_koide(m1, m2, m3):
    """Q = (m1+m2+m3) / (sqrt(m1)+sqrt(m2)+sqrt(m3))^2"""
    s = np.sqrt(m1) + np.sqrt(m2) + np.sqrt(m3)
    return (m1 + m2 + m3) / s**2

def Q_signed(m1, m2, m3):
    """Signed Koide: flip sign of smallest sqrt mass"""
    mm = sorted([m1, m2, m3])
    s = -np.sqrt(mm[0]) + np.sqrt(mm[1]) + np.sqrt(mm[2])
    if s == 0:
        return float('inf')
    return (mm[0] + mm[1] + mm[2]) / s**2

def Q_pair(m1, m2):
    """Q_pair = (m1+m2) / (sqrt(m1)+sqrt(m2))^2"""
    s = np.sqrt(m1) + np.sqrt(m2)
    return (m1 + m2) / s**2

def koide_parametrization(m1, m2, m3):
    """
    Koide parametrization: m_k = M0 (1 + sqrt(2) cos(delta + 2*pi*k/3))^2

    Given Q = 2/3 exactly:
      sum m_k = M0 * sum (1 + sqrt(2) cos(theta_k))^2
      = M0 * (3 + 2*sum cos^2(theta_k))  [since sum cos(theta_k) = 0]
      = M0 * (3 + 3) = 6*M0

    So M0 = (m1+m2+m3)/6.
    Then sqrt(m_k/M0) = 1 + sqrt(2) cos(delta + 2*pi*k/3)
    """
    sorted_m = sorted([m1, m2, m3])
    M0 = sum(sorted_m) / 6.0

    r = [np.sqrt(mi / M0) for mi in sorted_m]

    # r_k = 1 + sqrt(2) cos(delta + 2*pi*k/3)
    # cos(delta + 2*pi*k/3) = (r_k - 1) / sqrt(2)

    best_delta = None
    best_err = float('inf')

    for perm in [(0,1,2), (0,2,1), (1,0,2), (1,2,0), (2,0,1), (2,1,0)]:
        cos_vals = [(r[j] - 1) / np.sqrt(2) for j in range(3)]
        if any(abs(c) > 1.0 + 0.1 for c in cos_vals):
            continue

        for sign0 in [1, -1]:
            c0 = np.clip(cos_vals[0], -1, 1)
            angle0 = sign0 * np.arccos(c0)
            delta_try = angle0 - 2*np.pi*perm[0]/3

            err = 0
            for j in range(3):
                predicted = 1 + np.sqrt(2) * np.cos(delta_try + 2*np.pi*perm[j]/3)
                err += (predicted - r[j])**2

            if err < best_err:
                best_err = err
                best_delta = delta_try

    if best_delta is not None:
        best_delta = best_delta % (2 * np.pi)

    return M0, best_delta, best_err

# ============================================================
# TASK 1: All C(14,3) = 364 triples, sorted by |Q - 2/3|
# ============================================================

print("=" * 90)
print("TASK 1: ALL 364 TRIPLES — Q = (sum m) / (sum sqrt(m))^2")
print("Sorted by |Q - 2/3|")
print("=" * 90)

results = []
for i, j, k in combinations(range(N), 3):
    n1, n2, n3 = names[i], names[j], names[k]
    m1, m2, m3 = vals[i], vals[j], vals[k]
    Q = Q_koide(m1, m2, m3)
    Qs = Q_signed(m1, m2, m3)
    results.append((abs(Q - 2/3), Q, Qs, n1, n2, n3, m1, m2, m3))

results.sort(key=lambda x: x[0])

print(f"\nTop 30 closest to Q = 2/3 (unsigned):")
print(f"{'Rank':>4} {'Triple':>35} {'Q':>12} {'|Q-2/3|':>12} {'%dev':>8} {'Q_signed':>12} {'|Qs-2/3|':>12}")
print("-" * 105)

for rank, (dQ, Q, Qs, n1, n2, n3, m1, m2, m3) in enumerate(results[:30], 1):
    pct = dQ / (2/3) * 100
    label = f"({n1}, {n2}, {n3})"
    print(f"{rank:>4} {label:>35} {Q:>12.8f} {dQ:>12.2e} {pct:>7.4f}% {Qs:>12.8f} {abs(Qs-2/3):>12.2e}")

# ============================================================
# TASK 1b: Also sort by SIGNED variant
# ============================================================

print(f"\n\nTop 30 closest to Q = 2/3 (SIGNED variant, best of unsigned/signed):")
results_best = [(min(abs(Q-2/3), abs(Qs-2/3)), Q, Qs, n1, n2, n3, m1, m2, m3)
                 for _, Q, Qs, n1, n2, n3, m1, m2, m3 in results]
results_best.sort(key=lambda x: x[0])

print(f"{'Rank':>4} {'Triple':>35} {'Q':>12} {'Q_signed':>12} {'best|Q-2/3|':>14} {'which':>8}")
print("-" * 95)

for rank, (dQ, Q, Qs, n1, n2, n3, m1, m2, m3) in enumerate(results_best[:30], 1):
    label = f"({n1}, {n2}, {n3})"
    which = "signed" if abs(Qs-2/3) < abs(Q-2/3) else "unsign"
    print(f"{rank:>4} {label:>35} {Q:>12.8f} {Qs:>12.8f} {dQ:>14.2e} {which:>8}")

# ============================================================
# TASK 2: Detailed table for top 20 unsigned
# ============================================================

print("\n" + "=" * 90)
print("TASK 2: DETAILED VIEW — TOP 20 UNSIGNED + SIGNED VARIANTS")
print("=" * 90)

print(f"\n{'Rank':>4} {'Triple':>35} {'Q_unsign':>12} {'|Qu-2/3|':>12} {'Q_signed':>12} {'|Qs-2/3|':>12}")
print("-" * 95)

for rank, (dQ, Q, Qs, n1, n2, n3, m1, m2, m3) in enumerate(results[:20], 1):
    label = f"({n1}, {n2}, {n3})"
    print(f"{rank:>4} {label:>35} {Q:>12.8f} {abs(Q-2/3):>12.2e} {Qs:>12.8f} {abs(Qs-2/3):>12.2e}")

# ============================================================
# TASK 3: Physically motivated triples
# ============================================================

print("\n" + "=" * 90)
print("TASK 3: PHYSICALLY MOTIVATED TRIPLES")
print("=" * 90)

motivated = [
    ("(a) W, Z, H",       'm_W', 'm_Z', 'm_H'),
    ("(b) W, Z, t",       'm_W', 'm_Z', 'm_t'),
    ("(c) t, H, v_EW",    'm_t', 'm_H', 'v_EW'),
    ("(d) W, H, t",       'm_W', 'm_H', 'm_t'),
    ("(e) c, b, W",       'm_c', 'm_b', 'm_W'),
    ("(f) b, t, W",       'm_b', 'm_t', 'm_W'),
    ("(g) b, t, H",       'm_b', 'm_t', 'm_H'),
    ("(h) b, t, Z",       'm_b', 'm_t', 'm_Z'),
    ("(i) f_pi, tau, c",  'f_pi', 'm_tau', 'm_c'),
    ("(j) f_pi, b, v_EW", 'f_pi', 'm_b', 'v_EW'),
]

print(f"\n{'Label':>20} {'Q':>12} {'|Q-2/3|':>12} {'%dev':>8} {'Q_signed':>12} {'|Qs-2/3|':>12}")
print("-" * 90)

for label, n1, n2, n3 in motivated:
    m1, m2, m3 = masses[n1], masses[n2], masses[n3]
    Q = Q_koide(m1, m2, m3)
    Qs = Q_signed(m1, m2, m3)
    pct = abs(Q - 2/3) / (2/3) * 100
    print(f"{label:>20} {Q:>12.8f} {abs(Q-2/3):>12.2e} {pct:>7.3f}% {Qs:>12.8f} {abs(Qs-2/3):>12.2e}")

# ============================================================
# TASK 3b: Known Koide triples verification
# ============================================================

print("\n" + "=" * 90)
print("VERIFICATION: KNOWN KOIDE TRIPLES")
print("=" * 90)

known = [
    ("(e, mu, tau)",     'm_e', 'm_mu', 'm_tau'),
    ("(-s, c, b)",       'm_s', 'm_c', 'm_b'),
    ("(c, b, t)",        'm_c', 'm_b', 'm_t'),
    ("(d, s, b)",        'm_d', 'm_s', 'm_b'),
    ("(u, c, t)",        'm_u', 'm_c', 'm_t'),
]

print(f"\n{'Triple':>20} {'Q_unsigned':>14} {'|Qu-2/3|':>12} {'Q_signed':>14} {'|Qs-2/3|':>12}")
print("-" * 80)

for label, n1, n2, n3 in known:
    m1, m2, m3 = masses[n1], masses[n2], masses[n3]
    Q = Q_koide(m1, m2, m3)
    Qs = Q_signed(m1, m2, m3)
    print(f"{label:>20} {Q:>14.10f} {abs(Q-2/3):>12.2e} {Qs:>14.10f} {abs(Qs-2/3):>12.2e}")

# ============================================================
# TASK 4: Koide parametrization for triples within 5% of Q=2/3
# ============================================================

print("\n" + "=" * 90)
print("TASK 4: KOIDE PARAMETRIZATION m_k = M0*(1 + sqrt(2)*cos(delta + 2*pi*k/3))^2")
print("For triples within 5% of Q = 2/3 (unsigned)")
print("=" * 90)

within_5pct = [(dQ, Q, Qs, n1, n2, n3, m1, m2, m3)
               for dQ, Q, Qs, n1, n2, n3, m1, m2, m3 in results
               if abs(Q - 2/3) / (2/3) < 0.05]

print(f"\n{len(within_5pct)} triples within 5% of Q = 2/3")

# Simple fractions of pi to check
pi_fracs = []
for num in range(0, 25):
    for den in range(1, 13):
        fval = num/den
        if 0 <= fval <= 2:
            pi_fracs.append((num, den, fval))
pi_fracs = list(set(pi_fracs))
pi_fracs.sort(key=lambda x: x[2])

def find_nearest_frac(delta_over_pi):
    best = None
    best_err = float('inf')
    for num, den, fval in pi_fracs:
        err = abs(delta_over_pi - fval)
        if err < best_err:
            best_err = err
            best = f"{num}/{den}"
    return best, best_err

print(f"\n{'Triple':>35} {'Q':>10} {'M0(GeV)':>10} {'delta':>10} {'d/pi':>10} {'frac':>8} {'err':>8}")
print("-" * 105)

# First the known ones
for label, n1, n2, n3 in known:
    m1, m2, m3 = masses[n1], masses[n2], masses[n3]
    Q = Q_koide(m1, m2, m3)
    M0, delta, fit_err = koide_parametrization(m1, m2, m3)
    if delta is not None:
        d_pi = delta / np.pi
        frac, ferr = find_nearest_frac(d_pi)
        print(f"{'*'+label:>35} {Q:>10.6f} {M0:>10.4f} {delta:>10.6f} {d_pi:>10.6f} {frac:>8} {ferr:>8.4f}")

print("-" * 105)

for dQ, Q, Qs, n1, n2, n3, m1, m2, m3 in within_5pct:
    M0, delta, fit_err = koide_parametrization(m1, m2, m3)
    label = f"({n1}, {n2}, {n3})"

    if delta is not None:
        d_pi = delta / np.pi
        frac, ferr = find_nearest_frac(d_pi)
        marker = " <--" if ferr < 0.01 else ""
        print(f"{label:>35} {Q:>10.6f} {M0:>10.4f} {delta:>10.6f} {d_pi:>10.6f} {frac:>8} {ferr:>8.4f}{marker}")
    else:
        print(f"{label:>35} {Q:>10.6f} {M0:>10.4f if M0 else 'N/A':>10} {'N/A':>10} {'N/A':>10}")

# ============================================================
# TASK 5: Is (c,b,t) alone as heavy Koide? EW comparison
# ============================================================

print("\n" + "=" * 90)
print("TASK 5: HEAVY-SECTOR KOIDE — (c,b,t) vs EW TRIPLES")
print("=" * 90)

Q_cbt = Q_koide(masses['m_c'], masses['m_b'], masses['m_t'])
Q_cbt_s = Q_signed(masses['m_c'], masses['m_b'], masses['m_t'])
print(f"\nReference: Q(c,b,t) = {Q_cbt:.8f}, |Q-2/3| = {abs(Q_cbt-2/3):.2e}")
print(f"           Q_signed(c,b,t) = {Q_cbt_s:.8f}, |Qs-2/3| = {abs(Q_cbt_s-2/3):.2e}")

# Triples from {c,b,t,W,Z,H,v_EW}
heavy_all = {'m_c', 'm_b', 'm_t', 'm_W', 'm_Z', 'm_H', 'v_EW'}

print(f"\nAll triples from {{c,b,t,W,Z,H,v_EW}}, sorted by |Q-2/3|:")
print(f"{'Triple':>35} {'Q':>12} {'|Q-2/3|':>12} {'Q_signed':>12} {'|Qs-2/3|':>12}")
print("-" * 90)

heavy_results = []
for dQ, Q, Qs, n1, n2, n3, m1, m2, m3 in results:
    triple_set = {n1, n2, n3}
    if triple_set <= heavy_all:
        heavy_results.append((dQ, Q, Qs, n1, n2, n3))

heavy_results.sort(key=lambda x: x[0])
for dQ, Q, Qs, n1, n2, n3 in heavy_results:
    label = f"({n1}, {n2}, {n3})"
    marker = " ***" if abs(Q - 2/3) < 0.005 else ""
    print(f"{label:>35} {Q:>12.8f} {dQ:>12.2e} {Qs:>12.8f} {abs(Qs-2/3):>12.2e}{marker}")

# Triples with BEST signed Q
print(f"\nBest signed-Q triples from heavy set:")
heavy_signed = [(min(abs(Q-2/3), abs(Qs-2/3)), Q, Qs, n1, n2, n3) for _, Q, Qs, n1, n2, n3 in heavy_results]
heavy_signed.sort(key=lambda x: x[0])
print(f"{'Triple':>35} {'Q':>12} {'Q_signed':>12} {'best':>12}")
print("-" * 75)
for dQ, Q, Qs, n1, n2, n3 in heavy_signed[:15]:
    label = f"({n1}, {n2}, {n3})"
    print(f"{label:>35} {Q:>12.8f} {Qs:>12.8f} {dQ:>12.2e}")

# ============================================================
# TASK 5b: What about (b, t, H)?
# ============================================================

print("\n" + "=" * 90)
print("TASK 5b: SPOTLIGHT ON (b, t, H)")
print("=" * 90)

mb, mt, mH = masses['m_b'], masses['m_t'], masses['m_H']
Q_bth = Q_koide(mb, mt, mH)
Qs_bth = Q_signed(mb, mt, mH)
M0_bth, delta_bth, _ = koide_parametrization(mb, mt, mH)

print(f"Q(b, t, H) = {Q_bth:.10f}")
print(f"|Q - 2/3|  = {abs(Q_bth - 2/3):.2e}  ({abs(Q_bth-2/3)/(2/3)*100:.3f}%)")
print(f"Q_signed   = {Qs_bth:.10f}")
print(f"|Qs - 2/3| = {abs(Qs_bth - 2/3):.2e}")
if delta_bth is not None:
    print(f"M0 = {M0_bth:.4f} GeV")
    print(f"delta = {delta_bth:.6f} rad = {delta_bth/np.pi:.6f} pi")
    frac, ferr = find_nearest_frac(delta_bth/np.pi)
    print(f"Nearest fraction: {frac}*pi (error {ferr:.4f})")

    # Predict masses from M0, delta
    print(f"\nPrediction test (using M0, delta from fit):")
    for k in range(3):
        m_pred = M0_bth * (1 + np.sqrt(2) * np.cos(delta_bth + 2*np.pi*k/3))**2
        print(f"  k={k}: m_pred = {m_pred:.4f} GeV")
    print(f"  Actual: m_b={mb}, m_t={mt}, m_H={mH}")

# Compare with (e, mu, tau)
print(f"\nComparison:")
me, mmu, mtau = masses['m_e'], masses['m_mu'], masses['m_tau']
Q_emt = Q_koide(me, mmu, mtau)
M0_emt, delta_emt, _ = koide_parametrization(me, mmu, mtau)
print(f"  (e, mu, tau):  Q = {Q_emt:.10f}, |Q-2/3| = {abs(Q_emt-2/3):.2e}, delta/pi = {delta_emt/np.pi:.6f}")
print(f"  (b, t, H):     Q = {Q_bth:.10f}, |Q-2/3| = {abs(Q_bth-2/3):.2e}, delta/pi = {delta_bth/np.pi:.6f}")

# ============================================================
# BONUS: PAIR SCAN
# ============================================================

print("\n" + "=" * 90)
print("BONUS: PAIR SCAN — Q_pair = (a+b) / (sqrt(a)+sqrt(b))^2")
print("=" * 90)

print(f"\nQ_pair = 2/3 when mass ratio r = (2+sqrt(3))^2 = {(2+np.sqrt(3))**2:.4f}")
print(f"                          or r = (2-sqrt(3))^2 = {(2-np.sqrt(3))**2:.6f}")

pair_results = []
for i, j in combinations(range(N), 2):
    n1, n2 = names[i], names[j]
    m1, m2 = vals[i], vals[j]
    Qp = Q_pair(m1, m2)
    pair_results.append((abs(Qp - 2/3), Qp, n1, n2, m1, m2))

pair_results.sort(key=lambda x: x[0])

print(f"\nTop 25 pairs closest to Q_pair = 2/3:")
print(f"{'Rank':>4} {'Pair':>25} {'Q_pair':>12} {'|Q-2/3|':>12} {'ratio':>12}")
print("-" * 75)

for rank, (dQ, Qp, n1, n2, m1, m2) in enumerate(pair_results[:25], 1):
    ratio = max(m1, m2) / min(m1, m2)
    label = f"({n1}, {n2})"
    print(f"{rank:>4} {label:>25} {Qp:>12.8f} {dQ:>12.2e} {ratio:>12.4f}")

# ============================================================
# FULL TABLE: ALL 364 TRIPLES
# ============================================================

print("\n" + "=" * 90)
print("APPENDIX: ALL 364 TRIPLES (UNSIGNED Q)")
print("=" * 90)
print(f"\n{'Rank':>4} {'Triple':>35} {'Q':>12} {'|Q-2/3|':>12}")
print("-" * 70)

for rank, (dQ, Q, Qs, n1, n2, n3, m1, m2, m3) in enumerate(results, 1):
    label = f"({n1}, {n2}, {n3})"
    print(f"{rank:>4} {label:>35} {Q:>12.8f} {dQ:>12.2e}")
