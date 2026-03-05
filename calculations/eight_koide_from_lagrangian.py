"""
Theme 1 for the talk: do the 8 Koide equations from (t,s)x(b,u)x(c,d)
emerge from the Lagrangian?

Each quark pair contributes one member per triple.
Extended options: direct sqrt(m), magnetic 1/sqrt(m), signed -sqrt(m),
and limits u->0, t->inf.

For each triple close to Q=2/3, identify the Lagrangian origin.
"""

import math
import itertools

# PDG 2024 MSbar masses (MeV)
M = {'u': 2.16, 'd': 4.67, 's': 93.4, 'c': 1270.0, 'b': 4183.0, 't': 162500.0}

TARGET = 2/3

def Q(y1, y2, y3):
    s = y1 + y2 + y3
    if abs(s) < 1e-30: return float('nan')
    return (y1**2 + y2**2 + y3**2) / s**2

# The 8 base triples: one from each pair (t,s) x (u,b) x (c,d)
base_triples = [
    ('t', 'u', 'c'), ('t', 'u', 'd'), ('t', 'b', 'c'), ('t', 'b', 'd'),
    ('s', 'u', 'c'), ('s', 'u', 'd'), ('s', 'b', 'c'), ('s', 'b', 'd'),
]

# For each quark, the possible "charges" it can contribute:
# direct: sqrt(m), signed: -sqrt(m), magnetic: 1/sqrt(m), mag-signed: -1/sqrt(m)
# limits: u->0 gives charge 0; t->inf in magnetic gives charge 0

def charges(q):
    """Return list of (label, value) for quark q."""
    m = M[q]
    out = []
    out.append((f'+{q}', math.sqrt(m)))
    out.append((f'-{q}', -math.sqrt(m)))
    if m > 0:
        out.append((f'+1/{q}', 1/math.sqrt(m)))
        out.append((f'-1/{q}', -1/math.sqrt(m)))
    return out

# Also include limit versions
def charges_with_limits(q):
    """Like charges() but add the u->0 and t->inf limits."""
    out = charges(q)
    if q == 'u':
        out.append(('0(u)', 0.0))
    if q == 't':
        out.append(('0(1/t)', 0.0))  # 1/sqrt(inf) = 0
    return out

print("=" * 90)
print("ALL 8 BASE TRIPLES: (t,s) x (u,b) x (c,d)")
print("Each entry: direct, signed, magnetic, magnetic-signed, or limit")
print("Showing only |Q - 2/3| < 0.02 (1.5% relative)")
print("=" * 90)

results = []

for q1, q2, q3 in base_triples:
    for l1, v1 in charges_with_limits(q1):
        for l2, v2 in charges_with_limits(q2):
            for l3, v3 in charges_with_limits(q3):
                # Dimensional check: all nonzero finite entries should be same type
                types = []
                for lbl in [l1, l2, l3]:
                    if lbl.startswith('0('): continue
                    if '1/' in lbl: types.append('mag')
                    else: types.append('dir')
                dim_ok = len(set(types)) <= 1

                Qval = Q(v1, v2, v3)
                if math.isnan(Qval): continue
                dev = abs(Qval - TARGET)
                if dev < 0.02:
                    results.append({
                        'triple': f'({l1}, {l2}, {l3})',
                        'base': f'{q1}{q2}{q3}',
                        'Q': Qval,
                        'dev': dev,
                        'dim_ok': dim_ok,
                        'values': (v1, v2, v3),
                    })

results.sort(key=lambda r: r['dev'])

print(f"\n{'Rk':>3} {'Triple':<45} {'Base':>5} {'Q':>10} {'|Q-2/3|':>10} {'dim':>4}")
print("-" * 85)
for i, r in enumerate(results[:40], 1):
    d = 'ok' if r['dim_ok'] else ' *'
    print(f"{i:3d} {r['triple']:<45} {r['base']:>5} {r['Q']:10.6f} {r['dev']:10.6f} {d:>4}")

# Now classify by Lagrangian origin
print("\n" + "=" * 90)
print("CLASSIFICATION BY LAGRANGIAN ORIGIN")
print("=" * 90)

# Known Lagrangian mechanisms:
# A. O'Raifeartaigh seed: (0, m_s, m_c) => forced at gv/m = sqrt(3)
# B. Bion doubling: sqrt(m_b) = 3 sqrt(m_s) + sqrt(m_c)
# C. Yukawa eigenvalue: Q(c,b,t) = 2/3 at tan(beta)=1
# D. Seiberg seesaw: M_j = C/m_j => magnetic charges 1/sqrt(m_j)
# E. Lepton sector: SU(2) s-confinement with same O'R structure

for i, r in enumerate(results[:20], 1):
    trip = r['triple']
    base = r['base']
    Q_val = r['Q']
    dev_pct = r['dev'] / TARGET * 100

    # Try to identify mechanism
    mechanism = "?"
    # Check for seeds (one entry = 0)
    v1, v2, v3 = r['values']
    n_zero = sum(1 for v in [v1,v2,v3] if abs(v) < 1e-6)

    if n_zero >= 1:
        mechanism = "SEED (one mass = 0)"
    elif base == 'sbc' and r['dim_ok']:
        mechanism = "O'R + BION (electric quark chain)"
    elif base == 'tbc' and r['dim_ok']:
        mechanism = "YUKAWA eigenvalue (heavy triple)"
    elif '1/' in trip and r['dim_ok']:
        mechanism = "SEIBERG SEESAW (magnetic)"

    print(f"\n  {i}. {trip}")
    print(f"     Q = {Q_val:.6f} (dev {dev_pct:.2f}%), base = {base}")
    print(f"     Mechanism: {mechanism}")

# === Theme 1 specific: the t=1 vs t=inf question ===
print("\n" + "=" * 90)
print("THE t AND u LIMITS")
print("=" * 90)

print("""
The top quark enters the framework in two ways:
  1. DIRECT: sqrt(m_t) = 403.1 in triple (c, b, t)
  2. MAGNETIC LIMIT: 1/sqrt(m_t) = 0.00248 ≈ 0 (effectively a seed)

Similarly, the up quark:
  1. DIRECT: sqrt(m_u) = 1.47 ≈ 0 (effectively a seed)
  2. MAGNETIC: 1/sqrt(m_u) = 0.680 (comparable to 1/sqrt(m_d) = 0.463)
""")

# Direct (c,b,t) — the Yukawa eigenvalue triple
Q_cbt = Q(math.sqrt(M['c']), math.sqrt(M['b']), math.sqrt(M['t']))
print(f"Direct (c, b, t):     Q = {Q_cbt:.6f}")

# Magnetic (1/c, 1/b, 1/t) — t entry ≈ 0, so this is a seed
Q_mag_cbt = Q(1/math.sqrt(M['c']), 1/math.sqrt(M['b']), 1/math.sqrt(M['t']))
print(f"Magnetic (1/c, 1/b, 1/t): Q = {Q_mag_cbt:.6f}")

# Magnetic with t->0: seed (1/c, 1/b, 0)
Q_mag_cb0 = Q(1/math.sqrt(M['c']), 1/math.sqrt(M['b']), 0)
print(f"Seed (1/c, 1/b, 0):      Q = {Q_mag_cb0:.6f}")

# Direct with u->0: seed (0, c, d) for example
Q_seed_0cd = Q(0, math.sqrt(M['c']), math.sqrt(M['d']))
print(f"Seed (0, c, d):           Q = {Q_seed_0cd:.6f}")

# Magnetic u: 1/sqrt(m_u)
Q_mag_usd = Q(1/math.sqrt(M['u']), 1/math.sqrt(M['s']), 1/math.sqrt(M['d']))
print(f"Magnetic (1/u, 1/s, 1/d): Q = {Q_mag_usd:.6f}")

# The real question: which triples involving t or u are close to 2/3?
print(f"\nTriples involving t, close to Q=2/3:")
for r in results:
    if 't' in r['base']:
        print(f"  {r['triple']:<45} Q={r['Q']:.6f}  dev={r['dev']:.6f}  base={r['base']}")

print(f"\nTriples involving u, close to Q=2/3:")
for r in results:
    if 'u' in r['base']:
        print(f"  {r['triple']:<45} Q={r['Q']:.6f}  dev={r['dev']:.6f}  base={r['base']}")

# === The O'Raifeartaigh t parameter ===
print("\n" + "=" * 90)
print("O'RAIFEARTAIGH PARAMETER t = gv/m")
print("=" * 90)

print("""
The O'R fermion spectrum at general t:
  m_0 = 0 (goldstino)
  m_± = (sqrt(t²+1) ± t) m

Special values:
  t = 0:      m+ = m- = m (SUSY restored, degenerate)
  t = 1:      m+/m- = (sqrt(2)+1)/(sqrt(2)-1) = 3+2sqrt(2) ≈ 5.83
  t = sqrt(3): m+/m- = (2+sqrt(3))^2 ≈ 13.93 (Koide seed)
  t -> inf:   m+ -> 2tm, m- -> m/(2t) -> 0 (one mass decouples)

At t = sqrt(3), the UNIQUE property is Q = 2/3 exactly.
""")

# Verify: Q as function of t
print(f"  {'t':>6} {'m+/m':>10} {'m-/m':>10} {'m+/m-':>10} {'Q':>10}")
for t in [0, 0.5, 1.0, math.sqrt(2), math.sqrt(3), 2.0, 3.0, 5.0, 10.0, 100.0]:
    mp = (math.sqrt(t**2+1) + t)
    mm = (math.sqrt(t**2+1) - t)
    Qval = Q(0, math.sqrt(mm), math.sqrt(mp))
    ratio = mp/mm if mm > 1e-10 else float('inf')
    print(f"  {t:6.2f} {mp:10.4f} {mm:10.4f} {ratio:10.4f} {Qval:10.6f}")

print("""
Q = 2/3 is a SINGLE point: t = sqrt(3).
Not t = 1, not t = inf. The Lagrangian pins it.
As t -> inf, Q -> 1 (one mass dominates).
As t -> 0, Q -> 1/2 (degenerate pair + zero).
""")

# === Check: what if t = gv/m is NOT sqrt(3) for the top sector? ===
print("=" * 90)
print("WHAT IF THE TOP ENTERS VIA A DIFFERENT O'R VACUUM?")
print("=" * 90)

# The (c,b,t) triple has Q = 0.6627 (close but not exact 2/3)
# Does this correspond to a specific t parameter?
# If (0, m_c, m_b) is a seed with t_1 = sqrt(3), and then
# we need (m_c, m_b, m_t) to also satisfy Q = 2/3...
# These are different triples!

# The (c,b,t) triple does NOT have a zero entry.
# So it's NOT an O'R seed. It's a BLOOM.
# The bloom angle for (c,b,t) is:
z_cbt = [math.sqrt(M['c']), math.sqrt(M['b']), math.sqrt(M['t'])]
z0_cbt = sum(z_cbt) / 3
zi_cbt = [z - z0_cbt for z in z_cbt]
var_ratio = sum(x**2 for x in zi_cbt) / (3 * z0_cbt**2)
print(f"  (c, b, t) triple:")
print(f"    z0 = {z0_cbt:.2f}")
print(f"    <zi²>/z0² = {var_ratio:.6f}  (1.000 = exact Q=2/3)")
print(f"    Q = {Q(*z_cbt):.6f}")

# For (-s, c, b):
z_scb = [-math.sqrt(M['s']), math.sqrt(M['c']), math.sqrt(M['b'])]
z0_scb = sum(z_scb) / 3
zi_scb = [z - z0_scb for z in z_scb]
var_ratio_scb = sum(x**2 for x in zi_scb) / (3 * z0_scb**2)
print(f"\n  (-s, c, b) triple:")
print(f"    z0 = {z0_scb:.2f}")
print(f"    <zi²>/z0² = {var_ratio_scb:.6f}")
print(f"    Q = {Q(*z_scb):.6f}")

# For (1/d, 1/s, 1/b):
z_dsb = [1/math.sqrt(M['d']), 1/math.sqrt(M['s']), 1/math.sqrt(M['b'])]
z0_dsb = sum(z_dsb) / 3
zi_dsb = [z - z0_dsb for z in z_dsb]
var_ratio_dsb = sum(x**2 for x in zi_dsb) / (3 * z0_dsb**2)
print(f"\n  (1/d, 1/s, 1/b) triple:")
print(f"    z0 = {z0_dsb:.6f}")
print(f"    <zi²>/z0² = {var_ratio_dsb:.6f}")
print(f"    Q = {Q(*z_dsb):.6f}")

# === Relationship between the three triples ===
print("\n" + "=" * 90)
print("RELATIONSHIP: shared masses link the three triples")
print("=" * 90)

print(f"""
Triple 1 (electric seed+bloom):  (-s, c, b)     Q = {Q(*z_scb):.6f}
Triple 2 (heavy Yukawa):         (c, b, t)      Q = {Q(*z_cbt):.6f}
Triple 3 (magnetic seesaw):      (1/d, 1/s, 1/b) Q = {Q(*z_dsb):.6f}

Shared masses:
  Triples 1 & 2 share: m_c, m_b
  Triples 1 & 3 share: m_s, m_b
  All three share: m_b

The bion relation sqrt(m_b) = 3 sqrt(m_s) + sqrt(m_c) links triples 1 & 2.
The magnetic Koide (triple 3) then constrains m_d given m_s, m_b.

Free parameters after all three conditions:
  Input: m_s (scale), m_e, m_mu
  Predicted: m_c (O'R), m_b (bion), m_t (Yukawa), m_tau (lepton Koide)
  Predicted: m_d (magnetic Koide!)
  Remaining free: m_u, theta_23, theta_13, delta_CP

That's 7 free parameters, not 8.
The magnetic Koide eliminates m_d.
""")

# === Count the independent Koide equations ===
print("=" * 90)
print("HOW MANY INDEPENDENT KOIDE EQUATIONS?")
print("=" * 90)

print(f"""
From 6 quark masses, we have found:

  1. Q(-s, c, b) = 2/3     [electric, O'R + bion]
  2. Q(c, b, t) = 2/3      [heavy Yukawa at tan(beta)=1]
  3. Q(1/d, 1/s, 1/b) = 2/3  [magnetic seesaw]

Each equation eliminates one parameter:
  Eq 1: given m_s, fixes m_c/m_s ratio (→ m_c)
         Actually, Eq 1 + bion gives m_c AND m_b from m_s
  Eq 2: given m_c, m_b, fixes m_t
  Eq 3: given m_s, m_b, fixes m_d

Plus leptons:
  Eq 4: Q(e, mu, tau) = 2/3  [fixes m_tau from m_e, m_mu]

4 equations on 10 masses (6 quarks + 3 leptons + theta_C via GST).
Inputs: m_s, m_e, m_mu (3 parameters)
Predictions: m_c, m_b, m_t, m_d, m_tau, theta_C (6 predictions)
Remaining free: m_u, theta_23, theta_13, delta_CP (4 parameters)

Are there MORE Koide equations hiding in the 8 base triples?
""")

# Check ALL base triples systematically
print("\nSystematic check of all 8 base triples:")
print(f"  {'Triple':>12} {'Direct Q':>10} {'Best Q':>10} {'Best variant':>35}")
for q1, q2, q3 in base_triples:
    # Direct
    Qdir = Q(math.sqrt(M[q1]), math.sqrt(M[q2]), math.sqrt(M[q3]))

    # Find best variant (over signs and inversions)
    best_Q = None
    best_dev = 999
    best_label = ""
    for s1, s2, s3 in itertools.product([1, -1], repeat=3):
        for inv in [False, True]:  # all direct or all magnetic
            if inv:
                v1 = s1/math.sqrt(M[q1])
                v2 = s2/math.sqrt(M[q2])
                v3 = s3/math.sqrt(M[q3])
                lbl = f"({s1:+d}/√{q1},{s2:+d}/√{q2},{s3:+d}/√{q3})"
            else:
                v1 = s1*math.sqrt(M[q1])
                v2 = s2*math.sqrt(M[q2])
                v3 = s3*math.sqrt(M[q3])
                lbl = f"({s1:+d}√{q1},{s2:+d}√{q2},{s3:+d}√{q3})"
            Qval = Q(v1, v2, v3)
            if not math.isnan(Qval):
                dev = abs(Qval - TARGET)
                if dev < best_dev:
                    best_dev = dev
                    best_Q = Qval
                    best_label = lbl

    print(f"  {q1+q2+q3:>12} {Qdir:10.6f} {best_Q:10.6f} {best_label:>35}  (dev {best_dev:.4f})")

# Also check with limit substitutions
print(f"\n  With u->0 and t->inf limits:")
for q1, q2, q3 in base_triples:
    for (l1, v1_), (l2, v2_), (l3, v3_) in [(('t→∞', 0.0), (q2, math.sqrt(M[q2])), (q3, math.sqrt(M[q3]))),
                                               ((q1, math.sqrt(M[q1])), ('u→0', 0.0), (q3, math.sqrt(M[q3])))]:
        if q1 == 't' and l1 == 't→∞':
            Qval = Q(v1_, v2_, v3_)
            if not math.isnan(Qval):
                dev = abs(Qval - TARGET)
                if dev < 0.05:
                    print(f"    ({l1}, {l2}, {l3}): Q = {Qval:.6f}  (dev {dev:.4f})")
        if q2 == 'u' and l2 == 'u→0':
            Qval = Q(v1_, v2_, v3_)
            if not math.isnan(Qval):
                dev = abs(Qval - TARGET)
                if dev < 0.05:
                    print(f"    ({l1}, {l2}, {l3}): Q = {Qval:.6f}  (dev {dev:.4f})")
    # Magnetic limits
    if q1 == 't':
        # t->inf in magnetic = 0
        for s2, s3 in itertools.product([1, -1], repeat=2):
            v2m = s2/math.sqrt(M[q2])
            v3m = s3/math.sqrt(M[q3])
            Qval = Q(0, v2m, v3m)
            if not math.isnan(Qval):
                dev = abs(Qval - TARGET)
                if dev < 0.02:
                    print(f"    (0[1/t], {s2:+d}/√{q2}, {s3:+d}/√{q3}): Q = {Qval:.6f}  (dev {dev:.4f}) SEED")
