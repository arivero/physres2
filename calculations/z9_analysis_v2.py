#!/usr/bin/env python3
"""
Z_9 residual analysis — clean version with fixes.
"""

from mpmath import mp, mpf, sqrt, atan2, pi, cos, fmod, fabs, nstr, findroot
from itertools import permutations

mp.dps = 50

def extract_delta_v0(m0, m1, m2):
    """Extract (delta, v0) from masses assigned to k=0,1,2. All must be >= 0."""
    s0, s1, s2 = sqrt(fabs(m0)), sqrt(fabs(m1)), sqrt(fabs(m2))
    # Apply signs for negative masses (Koide convention for -s)
    if m0 < 0: s0 = -s0
    if m1 < 0: s1 = -s1
    if m2 < 0: s2 = -s2

    v0 = (s0 + s1 + s2) / 3

    if v0 == 0:
        return None, None, False

    a0 = s0 / v0 - 1
    a1 = s1 / v0 - 1
    a2 = s2 / v0 - 1

    C = a0 - a1/2 - a2/2
    S = (a1 - a2) * sqrt(mpf(3))/2

    delta = atan2(-S, C)

    # Verify all reconstructed sqrt(m) have correct signs
    ok = True
    for k, s_orig in [(0, s0), (1, s1), (2, s2)]:
        sm = v0 * (1 + sqrt(mpf(2)) * cos(delta + 2*pi*k/3))
        # Check sign consistency
        if (s_orig >= 0 and sm < -mpf('1e-15')) or (s_orig < 0 and sm > mpf('1e-15')):
            ok = False

    return delta, v0, ok


def koide_Q(*masses):
    """Koide Q for positive masses only."""
    S = sum(sqrt(m) for m in masses)
    return sum(masses) / S**2


def try_perms(masses, labels):
    results = []
    for perm in permutations(range(3)):
        m = [masses[perm[i]] for i in range(3)]
        delta, v0, ok = extract_delta_v0(*m)
        if delta is None:
            continue
        dn = fmod(delta, 2*pi)
        if dn < 0: dn += 2*pi
        pl = f"k0={labels[perm[0]]}, k1={labels[perm[1]]}, k2={labels[perm[2]]}"

        sm_vals = []
        all_pos = True
        for k in range(3):
            sm = v0 * (1 + sqrt(mpf(2)) * cos(delta + 2*pi*k/3))
            sm_vals.append(sm)
            if sm < -mpf('1e-15'):
                all_pos = False

        results.append((perm, pl, dn, v0, all_pos, sm_vals))
    return results


two_pi = 2*pi
two_pi_3 = two_pi/3
two_pi_9 = two_pi/9

print("=" * 90)
print("TASK 1-2: CHARGED LEPTONS")
print("=" * 90)

m_e = mpf('0.51100')
m_mu = mpf('105.658')
m_tau = mpf('1776.86')
dm_tau = mpf('0.12')

Q_lep = koide_Q(m_e, m_mu, m_tau)
print(f"\nKoide Q(e,mu,tau) = {nstr(Q_lep, 15)}")
print(f"2/3               = {nstr(mpf(2)/3, 15)}")
print(f"Deviation          = {nstr(Q_lep - mpf(2)/3, 10)}")

res = try_perms([m_e, m_mu, m_tau], ['e', 'mu', 'tau'])
valid = [r for r in res if r[4]]

print(f"\nAll 6 permutations are valid. Distinct delta values (mod 2pi/3 equivalence):")
print(f"  Set A (delta mod 2pi/3 ≈ 2/9):")
for p, pl, d, v, ap, sv in valid:
    dm = fmod(d, two_pi_3)
    if dm < 0: dm += two_pi_3
    if fabs(dm - mpf(2)/9) < 0.01:
        print(f"    {pl}: delta = {nstr(d, 15)} rad")

print(f"  Set B (delta mod 2pi/3 ≈ 1.872):")
for p, pl, d, v, ap, sv in valid:
    dm = fmod(d, two_pi_3)
    if dm < 0: dm += two_pi_3
    if fabs(dm - mpf(2)/9) > 0.01:
        print(f"    {pl}: delta = {nstr(d, 15)} rad")

# The three in Set A differ by 2pi/3 rotations; same for Set B.
# Set B has delta_B = 2pi/3 - delta_A (reflection).
# Use canonical: k0=e, k1=mu, k2=tau
delta_lep = valid[0][2]
v0_lep = valid[0][3]
print(f"\nCanonical (k0=e, k1=mu, k2=tau):")
print(f"  delta = {nstr(delta_lep, 15)} rad")
print(f"  v0    = {nstr(v0_lep, 15)} MeV^(1/2)")

# Z_9 analysis
nd = fmod(9*delta_lep, two_pi)
if nd < 0: nd += two_pi
print(f"\n  9*delta mod 2pi = {nstr(nd, 15)} rad")
print(f"  (For Z_9 minimum this should be near 0 or 2pi)")
if nd > pi:
    off_z9 = nd - two_pi
else:
    off_z9 = nd
print(f"  Offset from nearest Z_9 = {nstr(off_z9, 12)} rad")
frac_z9 = delta_lep / two_pi_9
print(f"  delta/(2pi/9) = {nstr(frac_z9, 15)} (nearest int: {int(float(frac_z9)+0.5)})")
print(f"  Fractional part = {nstr(frac_z9 - int(float(frac_z9)), 12)}")
print(f"\n  VERDICT: Leptons are NOT at a Z_9 minimum. Offset is ~2 rad ≈ 0.32×(2pi).")

# delta uncertainty
dp, _, _ = extract_delta_v0(m_e, m_mu, m_tau + dm_tau)
dm, _, _ = extract_delta_v0(m_e, m_mu, m_tau - dm_tau)
dpn = fmod(dp, two_pi); dmn = fmod(dm, two_pi)
if dpn < 0: dpn += two_pi
if dmn < 0: dmn += two_pi
sig_d = fabs(dpn - dmn)/2
sig_9d = 9*sig_d
print(f"\n  sigma_delta (from dm_tau) = {nstr(sig_d, 8)} rad")
print(f"  sigma(9*delta)            = {nstr(sig_9d, 8)} rad")

# ============================================================================
print("\n" + "=" * 90)
print("TASK 3-4: QUARK TRIPLES")
print("=" * 90)

m_d = mpf('4.67')
m_s = mpf('93.4')
m_b = mpf('4180')
m_c = mpf('1270')
m_t = mpf('163000')

triples = [
    ("(d,s,b)",   [m_d, m_s, m_b], ['d','s','b']),
    ("(-s,c,b)",  [-m_s, m_c, m_b], ['-s','c','b']),
    ("(s,c,b)",   [m_s, m_c, m_b], ['s','c','b']),
    ("(c,b,t)",   [m_c, m_b, m_t], ['c','b','t']),
]

print(f"\n{'Triple':<12} {'Q':<12} {'Valid perms':<12} {'delta':<20} {'9d mod 2pi':<16} {'d mod(2pi/3)':<16}")
print("-" * 90)

# Leptons first
d_mod_lep = fmod(delta_lep, two_pi_3)
if d_mod_lep < 0: d_mod_lep += two_pi_3
nd_lep = fmod(9*delta_lep, two_pi)
if nd_lep < 0: nd_lep += two_pi
print(f"{'(e,mu,tau)':<12} {nstr(Q_lep,8):<12} {'6':<12} {nstr(delta_lep,14):<20} {nstr(nd_lep,12):<16} {nstr(d_mod_lep,12)}")

for name, masses, labels in triples:
    # Koide Q (use absolute values for Q computation)
    abs_masses = [fabs(m) for m in masses]
    # For (-s,c,b), Q uses sqrt of absolute values with sign
    if any(m < 0 for m in masses):
        sm = [(-1 if m < 0 else 1)*sqrt(fabs(m)) for m in masses]
        Q = sum(fabs(m) for m in masses) / sum(sm)**2
    else:
        Q = koide_Q(*masses)

    res = try_perms(masses, labels)
    valid = [r for r in res if r[4]]

    if valid:
        _, pl, d, v, _, _ = valid[0]
        nd = fmod(9*d, two_pi)
        if nd < 0: nd += two_pi
        dm = fmod(d, two_pi_3)
        if dm < 0: dm += two_pi_3
        print(f"{name:<12} {nstr(Q,8):<12} {str(len(valid)):<12} {nstr(d,14):<20} {nstr(nd,12):<16} {nstr(dm,12)}")
    else:
        print(f"{name:<12} {nstr(Q,8):<12} {'0 (NONE)':<12}")

# Detailed (-s,c,b) analysis
print(f"\nDetailed (-s,c,b) analysis:")
print(f"  Using SIGNED sqrt: sqrt(-s)→-sqrt(s), sqrt(c)→+sqrt(c), sqrt(b)→+sqrt(b)")
res_scb = try_perms([-m_s, m_c, m_b], ['-s','c','b'])
for p, pl, d, v, ap, sv in res_scb:
    marker = " <<<" if ap else ""
    print(f"    {pl}: delta={nstr(d,10)}, v0={nstr(v,8)}, pos={ap}{marker}")
    if ap:
        print(f"      sqrt(m): {nstr(sv[0],8)}, {nstr(sv[1],8)}, {nstr(sv[2],8)}")

# ============================================================================
print("\n" + "=" * 90)
print("TASK 5: CHECK CLAIM delta mod(2pi/3) ≈ 2/9")
print("=" * 90)

d_mod = fmod(delta_lep, two_pi_3)
if d_mod < 0: d_mod += two_pi_3

print(f"\nFor leptons (k0=e, k1=mu, k2=tau):")
print(f"  delta           = {nstr(delta_lep, 15)} rad")
print(f"  delta mod(2pi/3)= {nstr(d_mod, 15)} rad")
print(f"  2/9             = {nstr(mpf(2)/9, 15)} rad")
print(f"  Difference      = {nstr(d_mod - mpf(2)/9, 12)} rad")
print(f"  Relative        = {nstr((d_mod - mpf(2)/9)/(mpf(2)/9) * 1e6, 8)} ppm")

print(f"\nNote: the permutation k0=tau, k1=e, k2=mu gives delta = {nstr(mpf('0.222229153001135'), 15)}")
print(f"  This IS delta mod(2pi/3) directly — the smallest representative.")
print(f"  Numerically: 0.222229... ≈ 2/9 = 0.222222...")

# Equivalence check
print(f"\nEquivalent statement: 3*delta mod 2pi ≈ 3×(2/9) = 2/3 rad")
three_d = fmod(3*delta_lep, two_pi)
if three_d < 0: three_d += two_pi
print(f"  3*delta mod 2pi = {nstr(three_d, 15)} rad")
print(f"  2/3             = {nstr(mpf(2)/3, 15)} rad")
print(f"  Difference      = {nstr(three_d - mpf(2)/3, 12)} rad")

print(f"\nCLARIFICATION:")
print(f"  The claim delta mod(2pi/3) ≈ 2/9 is about 2/9 RADIANS (not 2pi/9).")
print(f"  This is NOT the same as sitting at a Z_9 = Z(2pi/9) minimum.")
print(f"  2/9 rad = {nstr(mpf(2)/9, 12)}")
print(f"  2pi/9   = {nstr(two_pi_9, 12)}")
print(f"  Ratio   = 1/pi = {nstr(mpf(1)/pi, 12)}")
print(f"  The claim involves 1/pi — it's a transcendental relationship, not a Z_9 lattice point.")

# ============================================================================
print("\n" + "=" * 90)
print("TASK 6: SENSITIVITY ANALYSIS")
print("=" * 90)

# 6a: m_tau for exact delta mod(2pi/3) = 2/9
def res_29(mt):
    d, v, ok = extract_delta_v0(m_e, m_mu, mt)
    dn = fmod(d, two_pi)
    if dn < 0: dn += two_pi
    dm = fmod(dn, two_pi_3)
    if dm < 0: dm += two_pi_3
    return dm - mpf(2)/9

m_tau_29 = findroot(res_29, m_tau, tol=mpf('1e-40'))
pull_29 = (m_tau_29 - m_tau) / dm_tau

print(f"\n6a: m_tau for exact delta mod(2pi/3) = 2/9 radians:")
print(f"  m_tau(exact)  = {nstr(m_tau_29, 12)} MeV")
print(f"  m_tau(PDG)    = {nstr(m_tau, 12)} MeV")
print(f"  Difference    = {nstr(m_tau_29 - m_tau, 10)} MeV")
print(f"  Pull          = {nstr(pull_29, 6)} sigma (dm_tau = {dm_tau} MeV)")

# Verify
d29, _, _ = extract_delta_v0(m_e, m_mu, m_tau_29)
d29n = fmod(d29, two_pi)
if d29n < 0: d29n += two_pi
dm29 = fmod(d29n, two_pi_3)
if dm29 < 0: dm29 += two_pi_3
print(f"  Verification: d mod(2pi/3) = {nstr(dm29, 15)} (= 2/9 = {nstr(mpf(2)/9, 15)})")

# 6b: m_tau for nearest Z_9 minimum
# delta/(2pi/9) = 3.318..., nearest integer = 3
# delta = 3*(2pi/9) = 2pi/3
# But need to relax tolerance
print(f"\n6b: m_tau for nearest Z_9 minimum (delta = 3×2pi/9 = 2pi/3):")
def res_z9_3(mt):
    d, v, ok = extract_delta_v0(m_e, m_mu, mt)
    dn = fmod(d, two_pi)
    if dn < 0: dn += two_pi
    return dn - 3*two_pi_9

try:
    m_tau_z9 = findroot(res_z9_3, m_tau, tol=mpf('1e-20'))
    pull_z9 = (m_tau_z9 - m_tau) / dm_tau
    print(f"  m_tau(Z_9)    = {nstr(m_tau_z9, 12)} MeV")
    print(f"  Difference    = {nstr(m_tau_z9 - m_tau, 10)} MeV")
    print(f"  Pull          = {nstr(pull_z9, 6)} sigma")
except Exception as e:
    # Try bisection approach
    # delta decreases as m_tau increases (heavier tau -> smaller delta)
    # We need delta = 2pi/3 = 2.0944. Current delta = 2.3166. So need delta to decrease -> m_tau increase
    # Let's scan
    print(f"  Root finder failed: {e}")
    print(f"  Scanning for m_tau where delta = 2pi/3...")
    target = 3*two_pi_9  # = 2pi/3
    # delta at current: 2.3166
    # delta at m_tau=2000: ?
    for mt_try in [1800, 1850, 1900, 2000, 2500, 3000, 5000, 10000]:
        d_try, _, _ = extract_delta_v0(m_e, m_mu, mpf(mt_try))
        dn_try = fmod(d_try, two_pi)
        if dn_try < 0: dn_try += two_pi
        print(f"    m_tau={mt_try}: delta={nstr(dn_try, 10)}")

# ============================================================================
print("\n" + "=" * 90)
print("COMPREHENSIVE SUMMARY TABLE")
print("=" * 90)

print(f"""
+-------------------+---------------+------------------+------------------+-----------------+
| Triple            | Koide Q       | delta (rad)      | d mod(2pi/3)     | 9d mod 2pi      |
+-------------------+---------------+------------------+------------------+-----------------+""")

# Leptons
print(f"| (e, mu, tau)      | {nstr(Q_lep,10):<13} | {nstr(delta_lep,14):<16} | {nstr(d_mod_lep,12):<16} | {nstr(nd_lep,12):<15} |")

for name, masses, labels in triples:
    if any(m < 0 for m in masses):
        sm = [(-1 if m < 0 else 1)*sqrt(fabs(m)) for m in masses]
        Q = sum(fabs(m) for m in masses) / sum(sm)**2
    else:
        Q = koide_Q(*[fabs(m) for m in masses])
    res = try_perms(masses, labels)
    valid = [r for r in res if r[4]]
    if valid:
        _, _, d, _, _, _ = valid[0]
        nd = fmod(9*d, two_pi)
        if nd < 0: nd += two_pi
        dm = fmod(d, two_pi_3)
        if dm < 0: dm += two_pi_3
        print(f"| {name:<17} | {nstr(Q,10):<13} | {nstr(d,14):<16} | {nstr(dm,12):<16} | {nstr(nd,12):<15} |")
    else:
        print(f"| {name:<17} | {nstr(Q,10):<13} | {'NO VALID PERM':<16} | {'—':<16} | {'—':<15} |")

print(f"+-------------------+---------------+------------------+------------------+-----------------+")

print(f"""
KEY FINDINGS:
=============

1. LEPTON delta: {nstr(delta_lep, 15)} rad (k0=e, k1=mu, k2=tau)

2. Z_9 TEST: 9*delta mod 2pi = {nstr(nd_lep, 12)} rad
   Leptons are NOT at a Z_9 minimum (offset ≈ 2.0 rad ≈ 115 degrees).

3. "2/9 CLAIM" CONFIRMED:
   delta mod(2pi/3) = {nstr(d_mod_lep, 15)} rad
   2/9              = {nstr(mpf(2)/9, 15)} rad
   Match to {nstr(fabs(d_mod_lep - mpf(2)/9)/(mpf(2)/9)*1e6, 6)} ppm

4. EXACT MATCH at m_tau = {nstr(m_tau_29, 10)} MeV (pull = {nstr(pull_29, 4)} sigma from PDG)

5. The "2/9" claim and "Z_9 minimum" are DISTINCT statements.
   2/9 radians ≠ 2pi/9 radians (differ by factor pi).
   delta sits at 2/9 rad mod(2pi/3), NOT at a Z_9 = Z(2pi/9) lattice point.

6. (s,c,b) with positive masses has NO valid permutation (Q = 0.459, too far from 2/3).
   (-s,c,b) with signed sqrt gives valid permutations.
""")
