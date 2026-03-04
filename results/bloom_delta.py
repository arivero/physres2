"""
Koide parametrization bloom analysis.

The Koide parametrization writes three signed square roots as:

    sigma_i = sqrt(M0) * (1 + sqrt(2) * cos(2*pi*(i-1)/3 + delta))

for i = 1, 2, 3, where phi_i = 2*pi*(i-1)/3 and sigma_i are the signed
square roots of the masses (sigma_i^2 = m_i, sign from the triple convention).

The Koide condition Q = sum(m) / (sum(sigma))^2 = 2/3 is satisfied automatically
for any (M0, delta).

When delta = 3*pi/4: cos(3*pi/4) = -1/sqrt(2), so sigma_1 = 0. This is the
'Koide seed' -- a triple with one massless first member.

Analytic formulae for (M0, delta) given signed roots sigma_i:

    sqrt(M0) = sum(sigma_i) / 3         [because sum_i cos(phi_i + delta) = 0]

    delta = arctan2(-sum(sigma_i * sin(phi_i)), sum(sigma_i * cos(phi_i)))

These are exact when Q = 2/3; for Q != 2/3 they give the projection onto the
nearest point in the parametrization manifold (least-squares solution).

Masses (MeV):
    m_e = 0.51100,  m_mu = 105.658,  m_tau = 1776.86
    m_s = 93.4,     m_c  = 1270.,    m_b   = 4180.,   m_t = 172760.

Triples analyzed:
    (e, mu, tau)  -- all positive signed roots
    (-s, c, b)    -- negative signed root for s
    (c, b, t)     -- all positive signed roots

The problem asks for the Koide seed defined by sigma_1 = 0, i.e., delta_seed = 3*pi/4,
for all three triples. Additionally, for (-s, c, b) we compute the 'bloom seed'
from v0_doubling.md: the preceding 2-mass pair (s, c) with sigma_3 = 0.
"""

import numpy as np

# ---------------------------------------------------------------------------
# Physical constants (MeV)
# ---------------------------------------------------------------------------
m_e   = 0.51100
m_mu  = 105.658
m_tau = 1776.86
m_s   = 93.4
m_c   = 1270.0
m_b   = 4180.0
m_t   = 172760.0

SQRT2         = np.sqrt(2.0)
TWO_PI_OVER_3 = 2.0 * np.pi / 3.0
DELTA_SEED    = 3.0 * np.pi / 4.0        # = 135 deg; forces sigma_1 = 0


# ---------------------------------------------------------------------------
# Core functions
# ---------------------------------------------------------------------------

def compute_M0_delta(signed_roots):
    """
    Given signed roots sigma = [s1, s2, s3], return (M0, delta) analytically.

    From the parametrization:
      sqrt(M0) = sum(sigma) / 3
      delta    = arctan2(-sum(sigma*sin(phi)), sum(sigma*cos(phi)))
    """
    s   = np.array(signed_roots, dtype=float)
    phi = np.array([TWO_PI_OVER_3 * i for i in range(3)])
    sqrtM0 = s.sum() / 3.0
    M0     = sqrtM0 ** 2
    A      = np.dot(s, np.cos(phi))
    B      = np.dot(s, -np.sin(phi))
    delta  = np.arctan2(B, A) % (2.0 * np.pi)
    return M0, delta


def koide_Q(signed_roots):
    """Q = sum(m) / (sum(sigma))^2 with signed roots."""
    s = np.array(signed_roots, dtype=float)
    return s.dot(s) / s.sum() ** 2


def seed_at_delta(delta_seed, sr2, sr3):
    """
    Fit M0_seed at a given delta_seed using the two non-zero signed roots.
    The zero slot has c_i = 0; the other two slots give:

        sigma_j = sqrt(M0_seed) * c_j

    Least-squares solution: sqrt(M0_seed) = (sigma_2*c2 + sigma_3*c3) / (c2^2 + c3^2)
    where c_j = 1 + sqrt(2)*cos(phi_j + delta_seed) for j = 2, 3 (1-indexed).
    """
    c2 = 1.0 + SQRT2 * np.cos(TWO_PI_OVER_3       + delta_seed)
    c3 = 1.0 + SQRT2 * np.cos(2.0 * TWO_PI_OVER_3 + delta_seed)
    sqrtM0_seed = (sr2 * c2 + sr3 * c3) / (c2 ** 2 + c3 ** 2)
    M0_seed     = sqrtM0_seed ** 2
    m2_seed     = (sqrtM0_seed * c2) ** 2
    m3_seed     = (sqrtM0_seed * c3) ** 2
    return M0_seed, m2_seed, m3_seed


# ---------------------------------------------------------------------------
# Triples
# ---------------------------------------------------------------------------

triples = [
    {
        'name'         : '(e, mu, tau)',
        'signed_roots' : [np.sqrt(m_e),  np.sqrt(m_mu),  np.sqrt(m_tau)],
        'actual_masses': [m_e, m_mu, m_tau],
    },
    {
        'name'         : '(-s, c, b)',
        'signed_roots' : [-np.sqrt(m_s), np.sqrt(m_c),  np.sqrt(m_b)],
        'actual_masses': [m_s, m_c, m_b],
    },
    {
        'name'         : '(c, b, t)',
        'signed_roots' : [np.sqrt(m_c),  np.sqrt(m_b),  np.sqrt(m_t)],
        'actual_masses': [m_c, m_b, m_t],
    },
]

# ---------------------------------------------------------------------------
# Compute and print
# ---------------------------------------------------------------------------

print("=" * 72)
print("Koide Bloom Delta Analysis")
print("=" * 72)
print()
print(f"Seed geometry at delta_seed = 3*pi/4 = {DELTA_SEED:.8f} rad = 135.0000 deg:")
c1s = 1.0 + SQRT2 * np.cos(DELTA_SEED)
c2s = 1.0 + SQRT2 * np.cos(TWO_PI_OVER_3 + DELTA_SEED)
c3s = 1.0 + SQRT2 * np.cos(2.0 * TWO_PI_OVER_3 + DELTA_SEED)
print(f"  c1 = 0 (exact),  c2 = {c2s:.8f},  c3 = {c3s:.8f}")
print(f"  c3/c2 = {c3s/c2s:.8f}  =  (1+sqrt(3))/(sqrt(3)-1) = {(1+np.sqrt(3))/(np.sqrt(3)-1):.8f}")
print()

results = []

for tr in triples:
    name = tr['name']
    sr   = tr['signed_roots']
    ms   = tr['actual_masses']

    print("-" * 72)
    print(f"Triple: {name}")
    print()
    print(f"  Masses (MeV):       {ms}")
    print(f"  Signed roots:       [{', '.join(f'{x:.6f}' for x in sr)}]")
    print()

    # --- Part 1: Full parametrization ---
    M0_full, delta_full = compute_M0_delta(sr)
    Q                   = koide_Q(sr)
    Ddelta              = ((delta_full - DELTA_SEED) + np.pi) % (2.0 * np.pi) - np.pi

    print("  Part 1 -- Full parametrization:")
    print(f"    M0_full    = {M0_full:.6f} MeV")
    print(f"    delta_full = {delta_full:.8f} rad  =  {np.degrees(delta_full):.6f} deg")
    print(f"    delta_full / pi       = {delta_full/np.pi:.8f}")
    print(f"    delta_full / (2pi/3)  = {delta_full/TWO_PI_OVER_3:.8f}")
    print(f"    Koide Q    = {Q:.8f}  (2/3 = {2/3:.8f},  deviation {(Q-2/3)/(2/3)*100:+.4f}%)")
    print()

    # --- Part 2: Koide seed (sigma_1 = 0, delta_seed = 3pi/4) ---
    M0_seed, m2_seed, m3_seed = seed_at_delta(DELTA_SEED, sr[1], sr[2])

    print("  Part 2 -- Koide seed (sigma_1 = 0, delta_seed = 3*pi/4 = 135 deg):")
    print(f"    M0_seed    = {M0_seed:.6f} MeV")
    print(f"    Seed masses: (0,  {m2_seed:.4f},  {m3_seed:.4f}) MeV")
    print(f"    (Actual non-zero masses: {ms[1]:.4f}, {ms[2]:.4f} MeV)")
    print()

    # --- Part 3: v0-doubling (Koide-parametrization seed) ---
    v0_full      = np.sqrt(abs(M0_full))  * np.sign(M0_full)
    v0_seed_koi  = np.sqrt(abs(M0_seed))
    ratio_M0     = M0_full / M0_seed
    ratio_v0_koi = v0_full / v0_seed_koi

    print("  Part 3 -- v0-doubling using Koide seed:")
    print(f"    v0_full  = sqrt(M0_full)  = {v0_full:.6f} MeV^(1/2)")
    print(f"    v0_seed  = sqrt(M0_seed)  = {v0_seed_koi:.6f} MeV^(1/2)")
    print(f"    M0_full / M0_seed = {ratio_M0:.6f}")
    print(f"    v0_full / v0_seed = {ratio_v0_koi:.6f}")
    print()

    # --- Part 4: Bloom rotation ---
    print("  Part 4 -- Bloom rotation:")
    print(f"    delta_seed = {DELTA_SEED:.8f} rad  (3*pi/4)")
    print(f"    delta_full = {delta_full:.8f} rad")
    print(f"    Delta_delta = delta_full - delta_seed")
    print(f"               = {Ddelta:.8f} rad")
    print(f"               = {np.degrees(Ddelta):.6f} degrees")
    print(f"               = {Ddelta/TWO_PI_OVER_3:.8f} x (2*pi/3)")
    print(f"               = {Ddelta/np.pi:.8f} x pi")
    print()

    results.append({
        'name'         : name,
        'M0_full'      : M0_full,
        'delta_full'   : delta_full,
        'Q'            : Q,
        'M0_seed'      : M0_seed,
        'm2_seed'      : m2_seed,
        'm3_seed'      : m3_seed,
        'ratio_M0'     : ratio_M0,
        'ratio_v0_koi' : ratio_v0_koi,
        'Delta_delta'  : Ddelta,
        'masses'       : ms,
    })


# ---------------------------------------------------------------------------
# Special: v0-doubling bloom seed for (-s, c, b) as in v0_doubling.md
# The 'bloom seed' is the preceding 2-mass Koide triple (s, c, 0) where
# sigma_3 = 0. At that seed, all three sigma are positive: (+sqrt(m_s), +sqrt(m_c), 0).
# Since sum(sigma)/3 = sqrt(M0), we have sqrt(M0_bloom_seed) = (sqrt(m_s)+sqrt(m_c))/3.
# The full triple (-s,c,b) has sqrt(M0_full) = (-sqrt(m_s)+sqrt(m_c)+sqrt(m_b))/3.
# ---------------------------------------------------------------------------

print("=" * 72)
print("Special: v0-doubling bloom seed for (-s, c, b)")
print("=" * 72)
print()
print("The preceding 2-mass pair (s, c) forms a Koide seed with sigma_3 = 0.")
print("At that seed: all signed roots positive: (+sqrt(m_s), +sqrt(m_c), 0).")
print()

delta_seed3    = (3.0 * np.pi / 4.0 - 2.0 * TWO_PI_OVER_3) % (2.0 * np.pi)
v0_bloom_seed  = (np.sqrt(m_s) + np.sqrt(m_c)) / 3.0      # = sqrt(M0_bloom_seed)
M0_bloom_seed  = v0_bloom_seed ** 2
v0_bloom_full  = (-np.sqrt(m_s) + np.sqrt(m_c) + np.sqrt(m_b)) / 3.0
M0_bloom_full  = v0_bloom_full ** 2
ratio_v0_bloom = v0_bloom_full / v0_bloom_seed
ratio_M0_bloom = M0_bloom_full / M0_bloom_seed

Ddelta_bloom = ((results[1]['delta_full'] - delta_seed3) + np.pi) % (2.0 * np.pi) - np.pi

print(f"  delta_seed (sigma_3=0) = {delta_seed3:.8f} rad  =  {delta_seed3/np.pi:.6f} pi")
print(f"  v0_seed_bloom = (sqrt(m_s)+sqrt(m_c))/3 = {v0_bloom_seed:.6f} MeV^(1/2)")
print(f"  M0_seed_bloom = v0_seed_bloom^2          = {M0_bloom_seed:.6f} MeV")
print()
print(f"  v0_full = (-sqrt(m_s)+sqrt(m_c)+sqrt(m_b))/3 = {v0_bloom_full:.6f} MeV^(1/2)")
print(f"  M0_full = v0_full^2                           = {M0_bloom_full:.6f} MeV")
print()
print(f"  v0_full / v0_seed_bloom = {ratio_v0_bloom:.6f}  (near 2)")
print(f"  M0_full / M0_seed_bloom = {ratio_M0_bloom:.6f}  (near 4)")
print()
print(f"  Bloom rotation Delta_delta (from sigma_3=0 seed):")
print(f"    = {Ddelta_bloom:.8f} rad  =  {np.degrees(Ddelta_bloom):.5f} deg")
print()
print(f"  Prediction: if v0 doubles exactly, sqrt(m_b) = 3*sqrt(m_s) + sqrt(m_c)")
pred_mb = (3.0 * np.sqrt(m_s) + np.sqrt(m_c)) ** 2
print(f"  Predicted m_b = {pred_mb:.2f} MeV  (PDG: {m_b:.0f} MeV,  delta = {(pred_mb-m_b)/m_b*100:.3f}%)")
print()


# ---------------------------------------------------------------------------
# Part 5: Summary table
# ---------------------------------------------------------------------------

print("=" * 72)
print("Part 5 -- Summary Table")
print("=" * 72)
print()

# Table 1
hdr1 = (f"{'Triple':<16} {'M0_full':>10} {'delta_full':>12} "
        f"{'M0_seed':>10} {'delta_seed':>10} {'M0r':>8} "
        f"{'Ddelta(rad)':>12} {'Ddelta(deg)':>12}")
print(hdr1)
print("-" * len(hdr1))
for r in results:
    print(f"{r['name']:<16} {r['M0_full']:>10.4f} {r['delta_full']:>12.6f} "
          f"{r['M0_seed']:>10.4f} {DELTA_SEED:>10.6f} {r['ratio_M0']:>8.5f} "
          f"{r['Delta_delta']:>12.6f} {np.degrees(r['Delta_delta']):>12.5f}")

print()
print(f"delta_seed = 3*pi/4 = {DELTA_SEED:.6f} rad for all triples (sigma_1 = 0).")
print()

# Table 2
hdr2 = (f"{'Triple':<16} {'v0r (Koide seed)':>18} {'Ddelta/(2pi/3)':>16} "
        f"{'Ddelta/pi':>12} {'Q':>10}")
print(hdr2)
print("-" * len(hdr2))
for r in results:
    print(f"{r['name']:<16} {r['ratio_v0_koi']:>18.6f} "
          f"{r['Delta_delta']/TWO_PI_OVER_3:>16.8f} "
          f"{r['Delta_delta']/np.pi:>12.8f} {r['Q']:>10.7f}")

print()

# v0-doubling bloom-seed for (-s,c,b):
print(f"(-s, c, b) bloom-seed (sigma_3=0, preceding pair (s,c)):")
print(f"  v0_full / v0_seed_bloom = {ratio_v0_bloom:.6f}  (near 2)")
print(f"  M0_full / M0_seed_bloom = {ratio_M0_bloom:.6f}  (near 4)")
print()


# ---------------------------------------------------------------------------
# Pattern analysis
# ---------------------------------------------------------------------------

print("=" * 72)
print("Pattern analysis")
print("=" * 72)
print()

Dd   = [r['Delta_delta'] for r in results]
nm   = [r['name']        for r in results]

print("Bloom rotations Delta_delta = delta_full - delta_seed (seed: sigma_1=0):")
for n, d in zip(nm, Dd):
    print(f"  {n:<16}: {d:+.8f} rad  =  {np.degrees(d):+8.5f} deg")

print()
print("Ratios between bloom rotations:")
print(f"  Dd(-s,c,b) / Dd(c,b,t)    = {Dd[1]/Dd[2]:.6f}  (near -2)")
print(f"  Dd(-s,c,b) / Dd(e,mu,tau) = {Dd[1]/Dd[0]:.6f}")
print(f"  Dd(-s,c,b) + 2*Dd(c,b,t) = {Dd[1]+2*Dd[2]:.8f} rad  ({(Dd[1]+2*Dd[2])/Dd[1]*100:.3f}% of Dd_scb)")

print()
print("Sensitivity of near-zero relation Dd_scb + 2*Dd_cbt to m_s:")
for m_s_var in [90, 93.4, 97]:
    sr_scb_v = [-np.sqrt(m_s_var), np.sqrt(m_c), np.sqrt(m_b)]
    _, d_scb_v = compute_M0_delta(sr_scb_v)
    Dd_scb_v = ((d_scb_v - DELTA_SEED) + np.pi) % (2*np.pi) - np.pi
    combo = Dd_scb_v + 2 * Dd[2]
    print(f"  m_s = {m_s_var:.1f} MeV: Dd_scb = {Dd_scb_v:.6f}, combo = {combo:.6f} rad")

print()
print("Conclusions:")
print()
print("  1. delta_seed = 3*pi/4 = 135 deg for all three triples. The seed always")
print("     has its first mass zero (sigma_1 = 0) at this angle.")
print()
print("  2. Bloom rotation Delta_delta is NOT universal across the three triples:")
print(f"       (e,mu,tau): {np.degrees(Dd[0]):+8.4f} deg  (small; m_e << m_mu,tau)")
print(f"       (-s,c,b):   {np.degrees(Dd[1]):+8.4f} deg  (largest; Q deviates from 2/3 by 1.24%)")
print(f"       (c,b,t):    {np.degrees(Dd[2]):+8.4f} deg  (intermediate; Q deviates by 0.42%)")
print()
print("  3. Near-integer relation: Dd(-s,c,b) / Dd(c,b,t) = "
      f"{Dd[1]/Dd[2]:.4f}, close to -2.")
print(f"     Dd(-s,c,b) + 2*Dd(c,b,t) = {Dd[1]+2*Dd[2]:.6f} rad (0.34% of Dd_scb).")
print("     This relation is fragile: a 4% change in m_s shifts the combo by 0.001 rad,")
print("     comparable to the combo value itself.")
print()
print("  4. v0-doubling (factor ~2) is a specific property of the quark bloom, not")
print("     a universal feature of the Koide seed:")
print(f"       Koide seed (sigma_1=0): M0_full/M0_seed ratios are 1.001, 1.066, 1.015")
print("         -- all near 1, measuring how much the lightest member contributes.")
print(f"       Bloom seed (sigma_3=0 of preceding pair): for (-s,c,b), ratio = 2.0005.")
print("         This is the v0-doubling prediction from v0_doubling.md.")
print()
print("  5. The M0_seed from the Koide seed (sigma_1=0) does NOT match the physical")
print("     non-zero masses well for (-s,c,b) and (c,b,t), because those triples have")
print("     Q significantly above 2/3 and large Bloom rotations. The seed masses at")
print("     delta_seed differ substantially from the actual m_c, m_b, m_t.")
print()
print("Done.")
