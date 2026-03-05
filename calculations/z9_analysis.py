#!/usr/bin/env python3
"""
Z_9 residual analysis for Koide angular parametrization.

sqrt(m_k) = v0 [1 + sqrt(2) cos(delta + 2*pi*k/3)],  k=0,1,2

Given three masses, extract delta and check proximity to cos(9*delta) minima.
"""

from mpmath import mp, mpf, sqrt, atan2, pi, cos, sin, fmod, fabs, nstr
from itertools import permutations

mp.dps = 50

def extract_delta_and_v0(m0, m1, m2):
    """
    Given masses assigned as m_k for k=0,1,2, extract delta and v0.
    Returns (delta, v0, all_positive_flag).
    """
    if m0 <= 0 or m1 <= 0 or m2 <= 0:
        return None, None, False

    s0, s1, s2 = sqrt(m0), sqrt(m1), sqrt(m2)
    v0 = (s0 + s1 + s2) / 3

    a0 = s0 / v0 - 1
    a1 = s1 / v0 - 1
    a2 = s2 / v0 - 1

    # Projections onto cos/sin(2pi k/3)
    C = a0 + a1 * mpf('-0.5') + a2 * mpf('-0.5')
    S = a1 * sqrt(mpf(3))/2 + a2 * (-sqrt(mpf(3))/2)

    delta = atan2(-S, C)

    # Verify: all reconstructed sqrt(m_k) positive
    ok = True
    for k in range(3):
        sm_k = v0 * (1 + sqrt(mpf(2)) * cos(delta + 2*pi*k/3))
        if sm_k < mpf('-1e-20'):
            ok = False

    return delta, v0, ok


def compute_koide_Q(m1, m2, m3):
    return (m1 + m2 + m3) / (sqrt(m1) + sqrt(m2) + sqrt(m3))**2


def try_all_permutations(masses, labels):
    """Try all 6 permutations of mass-to-k assignment. Return list of valid ones."""
    results = []
    for perm in permutations(range(3)):
        m0, m1, m2 = masses[perm[0]], masses[perm[1]], masses[perm[2]]
        delta, v0, ok = extract_delta_and_v0(m0, m1, m2)
        if delta is not None:
            delta_norm = fmod(delta, 2*pi)
            if delta_norm < 0:
                delta_norm += 2*pi
            perm_label = f"k0={labels[perm[0]]}, k1={labels[perm[1]]}, k2={labels[perm[2]]}"

            # Check reconstructed sqrt values
            sm_vals = []
            all_pos = True
            for k in range(3):
                sm = v0 * (1 + sqrt(mpf(2)) * cos(delta + 2*pi*k/3))
                sm_vals.append(sm)
                if sm < mpf('-1e-20'):
                    all_pos = False

            results.append((perm, perm_label, delta_norm, v0, all_pos, sm_vals))
    return results


def print_sep():
    print("=" * 80)


# ============================================================================
# TASK 1: Charged leptons
# ============================================================================
print_sep()
print("TASK 1: Extract delta for charged leptons (e, mu, tau)")
print_sep()

m_e = mpf('0.51100')
m_mu = mpf('105.658')
m_tau = mpf('1776.86')

lepton_masses = [m_e, m_mu, m_tau]
lepton_labels = ['e', 'mu', 'tau']

Q_lep = compute_koide_Q(m_e, m_mu, m_tau)
print(f"\nKoide Q = {nstr(Q_lep, 15)}  (2/3 = {nstr(mpf(2)/3, 15)})")

results_lep = try_all_permutations(lepton_masses, lepton_labels)

for perm, perm_label, delta, v0, all_pos, sm_vals in results_lep:
    print(f"\n  {perm_label}: delta={nstr(delta,12)}, all_pos={all_pos}")
    print(f"    sqrt(m): {nstr(sm_vals[0],8)}, {nstr(sm_vals[1],8)}, {nstr(sm_vals[2],8)}")

# Find valid permutations
valid_lep = [(p, pl, d, v, ap, sv) for p, pl, d, v, ap, sv in results_lep if ap]
print(f"\nValid permutations (all sqrt positive): {len(valid_lep)}")

# Use first valid
perm, perm_label, delta_lep, v0_lep, _, _ = valid_lep[0]
print(f"\nUsing: {perm_label}")
print(f"  delta = {nstr(delta_lep, 15)} rad")
print(f"  v0    = {nstr(v0_lep, 15)} MeV^(1/2)")

# Show all valid deltas for reference
if len(valid_lep) > 1:
    print("\nAll valid deltas:")
    for _, pl, d, _, _, _ in valid_lep:
        print(f"  {pl}: delta = {nstr(d, 12)}")

# ============================================================================
# TASK 2: Z_9 residual for leptons
# ============================================================================
print()
print_sep()
print("TASK 2: Z_9 residual for leptons")
print_sep()

two_pi = 2 * pi
two_pi_9 = two_pi / 9

for _, perm_label, delta, _, _, _ in valid_lep:
    # 9*delta mod 2pi
    nd = fmod(9 * delta, two_pi)
    if nd < 0: nd += two_pi

    # Offset from nearest multiple of 2pi (i.e. nearest Z_9 minimum of cos(9d))
    if nd > pi:
        offset_9d = nd - two_pi
    else:
        offset_9d = nd

    # Also find nearest Z_9 minimum by checking delta/(2pi/9) fractional part
    n_nearest = int(float(delta / two_pi_9 + mpf('0.5')))
    offset_d = delta - n_nearest * two_pi_9

    print(f"\n  {perm_label}")
    print(f"    delta               = {nstr(delta, 15)} rad")
    print(f"    delta / (2pi/9)     = {nstr(delta / two_pi_9, 15)}")
    print(f"    Nearest Z_9 index   = {n_nearest}")
    print(f"    9*delta mod 2pi     = {nstr(nd, 15)} rad")
    print(f"    Offset (in 9d)      = {nstr(offset_9d, 12)} rad")
    print(f"    Offset (in d)       = {nstr(offset_d, 12)} rad")
    print(f"    Offset / (2pi/9)    = {nstr(offset_d / two_pi_9, 12)}")

# Uncertainty from dm_tau
dm_tau = mpf('0.12')
delta_plus, _, _ = extract_delta_and_v0(m_e, m_mu, m_tau + dm_tau)
delta_minus, _, _ = extract_delta_and_v0(m_e, m_mu, m_tau - dm_tau)
dp = fmod(delta_plus, two_pi);
if dp < 0: dp += two_pi
dm = fmod(delta_minus, two_pi);
if dm < 0: dm += two_pi

sigma_delta = fabs(dp - dm) / 2
sigma_9delta = 9 * sigma_delta

# Use the first valid delta for the offset
delta_use = valid_lep[0][2]
nd_use = fmod(9 * delta_use, two_pi)
if nd_use < 0: nd_use += two_pi
if nd_use > pi:
    offset_9d_use = nd_use - two_pi
else:
    offset_9d_use = nd_use

print(f"\nUncertainty (dm_tau = {dm_tau} MeV):")
print(f"  sigma_delta   = {nstr(sigma_delta, 10)} rad")
print(f"  sigma_9delta  = {nstr(sigma_9delta, 10)} rad")
print(f"  |offset_9d| / sigma_9d = {nstr(fabs(offset_9d_use)/sigma_9delta, 6)} sigma")

# ============================================================================
# TASK 3: Quark triples
# ============================================================================
print()
print_sep()
print("TASK 3: Extract delta for quark triples")
print_sep()

quark_triples = [
    ("(d, s, b)", [mpf('4.67'), mpf('93.4'), mpf('4180')], ['d', 's', 'b']),
    ("(s, c, b)", [mpf('93.4'), mpf('1270'), mpf('4180')], ['s', 'c', 'b']),
    ("(c, b, t)", [mpf('1270'), mpf('4180'), mpf('163000')], ['c', 'b', 't']),
]

quark_results = {}

for name, masses, labels in quark_triples:
    print(f"\n--- {name} ---")
    Q = compute_koide_Q(masses[0], masses[1], masses[2])
    print(f"Koide Q = {nstr(Q, 10)}  (2/3 = {nstr(mpf(2)/3, 10)})")

    results = try_all_permutations(masses, labels)
    valid = [(p, pl, d, v, ap, sv) for p, pl, d, v, ap, sv in results if ap]

    for _, pl, d, _, ap, sv in results:
        marker = " <<<" if ap else ""
        print(f"  {pl}: delta={nstr(d,10)}, pos={ap}{marker}")

    if valid:
        _, pl, d, v, _, _ = valid[0]
        quark_results[name] = (pl, d, v)
        print(f"\n  Selected: {pl}, delta = {nstr(d, 15)}")
    else:
        print(f"\n  NO valid permutation found!")
        quark_results[name] = None

# ============================================================================
# TASK 4: Z_9 residual for each triple
# ============================================================================
print()
print_sep()
print("TASK 4: Z_9 residuals for all triples")
print_sep()

print(f"\n{'Triple':<15} {'Q':<12} {'delta':<18} {'d/(2pi/9)':<18} {'9d mod 2pi':<16} {'Z9 offset (rad)'}")
print("-" * 95)

# Leptons
nd = fmod(9 * delta_use, two_pi)
if nd < 0: nd += two_pi
if nd > pi: off = nd - two_pi
else: off = nd
n_near = int(float(delta_use / two_pi_9 + 0.5))
off_d = delta_use - n_near * two_pi_9
print(f"{'(e,mu,tau)':<15} {nstr(Q_lep,8):<12} {nstr(delta_use,14):<18} {nstr(delta_use/two_pi_9,12):<18} {nstr(nd,12):<16} {nstr(off,10)}")

for name, masses, labels in quark_triples:
    Q = compute_koide_Q(masses[0], masses[1], masses[2])
    if quark_results[name]:
        pl, d, v = quark_results[name]
        nd = fmod(9 * d, two_pi)
        if nd < 0: nd += two_pi
        if nd > pi: off = nd - two_pi
        else: off = nd
        print(f"{name:<15} {nstr(Q,8):<12} {nstr(d,14):<18} {nstr(d/two_pi_9,12):<18} {nstr(nd,12):<16} {nstr(off,10)}")

print("\nAre leptons at a Z_9 minimum?")
nd_lep = fmod(9 * delta_use, two_pi)
if nd_lep < 0: nd_lep += two_pi
if nd_lep > pi: off_lep = nd_lep - two_pi
else: off_lep = nd_lep
print(f"  9*delta mod 2pi = {nstr(nd_lep, 12)}, offset = {nstr(off_lep, 10)} rad")
print(f"  This is {nstr(fabs(off_lep)/(two_pi), 6)} of a full cycle — NOT near a Z_9 minimum.")

# ============================================================================
# TASK 5: Check "delta_0 mod(2pi/3) ≈ 2/9"
# ============================================================================
print()
print_sep()
print("TASK 5: Check claim 'delta mod(2pi/3) ≈ 2/9'")
print_sep()

two_pi_3 = two_pi / 3
d_mod = fmod(delta_use, two_pi_3)
if d_mod < 0: d_mod += two_pi_3

print(f"\ndelta = {nstr(delta_use, 15)} rad")
print(f"delta mod (2pi/3) = {nstr(d_mod, 15)} rad")
print(f"2/9               = {nstr(mpf(2)/9, 15)} rad")
print(f"Difference         = {nstr(d_mod - mpf(2)/9, 10)} rad")
print(f"Relative diff      = {nstr((d_mod - mpf(2)/9) / (mpf(2)/9), 8)}")
print(f"  = {nstr((d_mod - mpf(2)/9) / (mpf(2)/9) * 1e6, 6)} ppm")

# Equivalently: 3*delta mod 2pi should be close to 3*(2/9) = 2/3 rad
three_d = fmod(3 * delta_use, two_pi)
if three_d < 0: three_d += two_pi
print(f"\n3*delta mod 2pi    = {nstr(three_d, 15)} rad")
print(f"3 × (2/9) = 2/3   = {nstr(mpf(2)/3, 15)} rad")
print(f"Difference         = {nstr(three_d - mpf(2)/3, 10)} rad")

# Is this the SAME as the Z_9 statement? Let's clarify.
# Z_9: delta = n*(2pi/9) for integer n. That means 9*delta = 2pi*n.
# The claim: delta mod(2pi/3) ≈ 2/9 (radians, not fractions!)
# delta ≈ 2/9 + k*(2pi/3) for some integer k
# delta ≈ 2/9 + 2*pi*k/3
# 9*delta ≈ 2 + 6*pi*k = 2 mod(6pi) -- NOT 2pi*n

# So this is a DIFFERENT claim from Z_9. Let's compute both.
print(f"\n--- Relationship between the two claims ---")
print(f"Z_9 claim: delta = n*(2pi/9) => 9*delta mod 2pi = 0")
print(f"  9*delta mod 2pi = {nstr(nd_lep, 12)} (FAR from 0)")
print(f"")
print(f"2/9 claim: delta mod(2pi/3) = 2/9 radians")
print(f"  delta mod(2pi/3) = {nstr(d_mod, 12)} rad")
print(f"  2/9              = {nstr(mpf(2)/9, 12)} rad")
print(f"  Match to {nstr(fabs(d_mod - mpf(2)/9)/(mpf(2)/9)*1e6, 4)} ppm")
print(f"")
print(f"Note: 2/9 radians ≠ 2pi/9 radians. The claim is delta mod(2pi/3) ≈ 2/9 (bare radians).")
print(f"  2/9  = {nstr(mpf(2)/9, 12)} rad")
print(f"  2pi/9 = {nstr(two_pi_9, 12)} rad")

# ============================================================================
# TASK 6: Sensitivity — m_tau for delta mod(2pi/3) = EXACTLY 2/9
# ============================================================================
print()
print_sep()
print("TASK 6: Sensitivity analysis")
print_sep()

from mpmath import findroot

# 6a: m_tau for exact delta mod(2pi/3) = 2/9
def residual_29(mt):
    d, v, ok = extract_delta_and_v0(m_e, m_mu, mt)
    d_n = fmod(d, two_pi)
    if d_n < 0: d_n += two_pi
    d_mod = fmod(d_n, two_pi_3)
    if d_mod < 0: d_mod += two_pi_3
    return d_mod - mpf(2)/9

print("6a: m_tau for exact delta mod(2pi/3) = 2/9")
try:
    m_tau_29 = findroot(residual_29, m_tau, tol=mpf('1e-30'))
    pull_29 = (m_tau_29 - m_tau) / dm_tau
    d_29, _, _ = extract_delta_and_v0(m_e, m_mu, m_tau_29)
    d_29n = fmod(d_29, two_pi)
    if d_29n < 0: d_29n += two_pi
    print(f"  m_tau = {nstr(m_tau_29, 12)} MeV")
    print(f"  PDG   = {nstr(m_tau, 12)} MeV")
    print(f"  Diff  = {nstr(m_tau_29 - m_tau, 10)} MeV")
    print(f"  Pull  = {nstr(pull_29, 6)} sigma")
    print(f"  delta = {nstr(d_29n, 15)} rad")
    dm29 = fmod(d_29n, two_pi_3)
    if dm29 < 0: dm29 += two_pi_3
    print(f"  d mod(2pi/3) = {nstr(dm29, 15)} (should be 2/9 = {nstr(mpf(2)/9,15)})")
except Exception as e:
    print(f"  Error: {e}")

# 6b: m_tau for exact Z_9 minimum (nearest)
# delta / (2pi/9) = 3.318... so nearest is n=3
# Want delta = 3*(2pi/9) = 6pi/9 = 2pi/3
print("\n6b: m_tau for exact Z_9 minimum (nearest: n=3, delta = 2pi/3)")
def residual_z9(mt):
    d, v, ok = extract_delta_and_v0(m_e, m_mu, mt)
    d_n = fmod(d, two_pi)
    if d_n < 0: d_n += two_pi
    return d_n - 3 * two_pi_9

try:
    m_tau_z9 = findroot(residual_z9, m_tau, tol=mpf('1e-30'))
    pull_z9 = (m_tau_z9 - m_tau) / dm_tau
    print(f"  m_tau = {nstr(m_tau_z9, 12)} MeV")
    print(f"  Diff  = {nstr(m_tau_z9 - m_tau, 10)} MeV")
    print(f"  Pull  = {nstr(pull_z9, 6)} sigma")
except Exception as e:
    print(f"  Error: {e}")

# Also try n=4: delta = 4*(2pi/9)
print("\n6c: m_tau for Z_9 at n=4, delta = 4*(2pi/9)")
def residual_z9_4(mt):
    d, v, ok = extract_delta_and_v0(m_e, m_mu, mt)
    d_n = fmod(d, two_pi)
    if d_n < 0: d_n += two_pi
    return d_n - 4 * two_pi_9

try:
    m_tau_z9_4 = findroot(residual_z9_4, mpf('1800'), tol=mpf('1e-30'))
    pull_z9_4 = (m_tau_z9_4 - m_tau) / dm_tau
    print(f"  m_tau = {nstr(m_tau_z9_4, 12)} MeV")
    print(f"  Diff  = {nstr(m_tau_z9_4 - m_tau, 10)} MeV")
    print(f"  Pull  = {nstr(pull_z9_4, 6)} sigma")
except Exception as e:
    print(f"  Error: {e}")

# ============================================================================
# Now redo with DIFFERENT convention: maybe the "standard" Koide convention
# has delta defined differently. Let's try: k=0->tau, k=1->mu, k=2->e
# (heaviest first) since in Koide's original, the phase is in 2nd quadrant
# ============================================================================
print()
print_sep()
print("CROSS-CHECK: delta for all valid lepton permutations — Z_9 and 2/9 tests")
print_sep()

for perm, perm_label, delta, v0, all_pos, sm_vals in results_lep:
    if not all_pos:
        continue

    # Z_9 test
    nd = fmod(9 * delta, two_pi)
    if nd < 0: nd += two_pi
    if nd > pi: off_z9 = nd - two_pi
    else: off_z9 = nd

    # 2/9 test
    dm = fmod(delta, two_pi_3)
    if dm < 0: dm += two_pi_3
    off_29 = dm - mpf(2)/9

    print(f"\n  {perm_label}")
    print(f"    delta = {nstr(delta, 15)} rad")
    print(f"    Z_9:  9d mod 2pi = {nstr(nd, 12)}, offset = {nstr(off_z9, 10)}")
    print(f"    2/9:  d mod(2pi/3) = {nstr(dm, 12)}, offset from 2/9 = {nstr(off_29, 10)} ({nstr(off_29/(mpf(2)/9)*1e6, 4)} ppm)")

# ============================================================================
# FINAL SUMMARY
# ============================================================================
print()
print_sep()
print("FINAL SUMMARY")
print_sep()

print(f"""
LEPTON TRIPLE (e, mu, tau):
  Koide Q = {nstr(Q_lep, 12)} (2/3 to {nstr(fabs(Q_lep - mpf(2)/3)/(mpf(2)/3)*100, 4)}%)
  delta = {nstr(delta_use, 15)} rad (permutation: {valid_lep[0][1]})
  v0 = {nstr(valid_lep[0][3], 12)} MeV^(1/2)

Z_9 MINIMUM TEST:
  9*delta mod 2pi = {nstr(nd_lep, 12)} rad
  Leptons are NOT at a Z_9 minimum.
  Offset = {nstr(fabs(off_lep), 10)} rad = {nstr(fabs(off_lep)/two_pi*100, 4)}% of full cycle.

"delta mod(2pi/3) = 2/9" CLAIM:
  delta mod(2pi/3) = {nstr(d_mod, 15)} rad
  2/9              = {nstr(mpf(2)/9, 15)} rad
  Match: {nstr(fabs(d_mod - mpf(2)/9)/(mpf(2)/9)*1e6, 4)} ppm
  This is a remarkable numerical coincidence but is NOT the same as a Z_9 minimum.
  2/9 radians ≠ 2pi/9 radians.
""")

# Recheck: maybe the claim is stated differently. Let's check delta_0 mod(2pi/3) = 2/9
# where delta_0 is defined differently. In the memory: "δ₀ mod 2π/3 ≈ 2/9"
# This literally means: take delta, reduce mod (2pi/3), and the result is ≈ 2/9 (radians).
# We confirmed: 0.222229... ≈ 0.222222... = 2/9. Match to 33 ppm.
# The memory also says "33 ppm at PDG 2024 central m_τ=1776.86"
print(f"Confirmed: delta mod(2pi/3) ≈ 2/9 to {nstr(fabs(d_mod - mpf(2)/9)/(mpf(2)/9)*1e6, 4)} ppm")
print(f"(Memory states 33 ppm — let's check: {nstr(fabs(d_mod - mpf(2)/9)/(mpf(2)/9)*1e6, 6)} ppm)")
