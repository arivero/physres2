"""
Koide invariant under mass inversion (dual/seesaw Koide).

The Koide quotient for three masses (m_1, m_2, m_3) with signs (s_1, s_2, s_3):

    Q = (m_1 + m_2 + m_3) / (s_1 sqrt(m_1) + s_2 sqrt(m_2) + s_3 sqrt(m_3))^2

Tasks:
  1. Compute Q for several standard and non-standard triples.
  2. Algebraic relation between Q and Q_dual under m_j -> C/m_j.
  3. Numerical Q_dual values for down-type quarks.
  4. SQCD seesaw M_j = C/m_j with Lambda = 300 MeV.
  5. Explore connection between dual Koide and direct Koide through seesaw.
  6. Check whether Q = 2/3 for (m_1, m_2, m_3) implies Q = 2/3 for (1/m_1, 1/m_2, 1/m_3).
"""

import numpy as np

# ---------------------------------------------------------------------------
# Masses (MeV)
# ---------------------------------------------------------------------------
m_e   = 0.51100
m_mu  = 105.658
m_tau = 1776.86

m_u   = 2.16
m_d   = 4.67
m_s   = 93.4
m_c   = 1270.0
m_b   = 4180.0
m_t   = 172760.0


def koide_Q(masses, signs=None):
    """
    Compute Q = sum(m_j) / (sum(s_j * sqrt(m_j)))^2.
    signs: list of +1/-1 for each mass. Default: all +1.
    """
    ms = np.array(masses, dtype=float)
    if signs is None:
        signs = np.ones(len(ms))
    ss = np.array(signs, dtype=float)
    numerator = np.sum(ms)
    denominator = np.sum(ss * np.sqrt(ms))**2
    return numerator / denominator


# ═══════════════════════════════════════════════════════════════════════════
# PART 1: Direct Koide Q for standard triples
# ═══════════════════════════════════════════════════════════════════════════
print("=" * 72)
print("PART 1: Direct Koide Q for several mass triples")
print("=" * 72)
print()

triples_part1 = [
    ("(e, mu, tau)  +,+,+",       [m_e, m_mu, m_tau],     [+1, +1, +1]),
    ("(-s, c, b)    -,+,+",       [m_s, m_c, m_b],        [-1, +1, +1]),
    ("(c, b, t)     +,+,+",       [m_c, m_b, m_t],        [+1, +1, +1]),
    ("(d, s, b)     +,+,+",       [m_d, m_s, m_b],        [+1, +1, +1]),
    ("(-d, s, b)    -,+,+",       [m_d, m_s, m_b],        [-1, +1, +1]),
]

for name, masses, signs in triples_part1:
    Q = koide_Q(masses, signs)
    dev_pct = (Q - 2/3) / (2/3) * 100
    sqrts = [s * np.sqrt(m) for s, m in zip(signs, masses)]
    print(f"  {name:<28}  Q = {Q:.8f}   (2/3 dev: {dev_pct:+.4f}%)")
    print(f"    masses: {masses}")
    print(f"    signed sqrt: [{', '.join(f'{x:.6f}' for x in sqrts)}]")
    print()


# ═══════════════════════════════════════════════════════════════════════════
# PART 2: Algebraic relation Q vs Q_dual under m_j -> C/m_j
# ═══════════════════════════════════════════════════════════════════════════
print("=" * 72)
print("PART 2: Algebraic analysis — Q vs Q_dual under mass inversion")
print("=" * 72)
print()

print("""Given masses (m1, m2, m3) with all signs positive:

  Q = (m1 + m2 + m3) / (sqrt(m1) + sqrt(m2) + sqrt(m3))^2

Under M_j = C/m_j (Seiberg seesaw):

  Q_dual = (C/m1 + C/m2 + C/m3) / (sqrt(C/m1) + sqrt(C/m2) + sqrt(C/m3))^2

Factor out C:

  Numerator:   C * (1/m1 + 1/m2 + 1/m3)
  Denominator: C * (1/sqrt(m1) + 1/sqrt(m2) + 1/sqrt(m3))^2

Therefore C cancels:

  Q_dual = (1/m1 + 1/m2 + 1/m3) / (1/sqrt(m1) + 1/sqrt(m2) + 1/sqrt(m3))^2

This is independent of C. Rewrite in terms of a_j = sqrt(m_j):

  Q      = (a1^2 + a2^2 + a3^2)                 / (a1 + a2 + a3)^2
  Q_dual = (1/a1^2 + 1/a2^2 + 1/a3^2)           / (1/a1 + 1/a2 + 1/a3)^2

Multiply num and denom of Q_dual by (a1*a2*a3)^2:

  Q_dual = (a2^2*a3^2 + a1^2*a3^2 + a1^2*a2^2) / (a2*a3 + a1*a3 + a1*a2)^2

Using symmetric polynomials of (a1, a2, a3):
  p1 = a1+a2+a3,   p2 = a1*a2+a2*a3+a1*a3,   p3 = a1*a2*a3

  Q      = (p1^2 - 2*p2) / p1^2    = 1 - 2*p2/p1^2
  Q_dual = (p2^2 - 2*p1*p3) / p2^2 = 1 - 2*p1*p3/p2^2

So Q = 2/3  iff  p2/p1^2 = 1/6
   Q_dual = 2/3  iff  p1*p3/p2^2 = 1/6

These are INDEPENDENT conditions. Q = 2/3 does NOT imply Q_dual = 2/3.
""")


# ═══════════════════════════════════════════════════════════════════════════
# Koide parametrization analysis of Q_dual
# ═══════════════════════════════════════════════════════════════════════════
print("--- Koide parametrization analysis ---")
print()
print("""For the Koide parametrization (all positive signs, all sv_k > 0):

  sv_k = sqrt(M0) * (1 + sqrt(2)*cos(2*pi*k/3 + delta)),  k=0,1,2
  m_k  = sv_k^2

This gives Q = 2/3 identically. The symmetric polynomials of the sv_k are:

  p1 = sum sv_k = 3*sqrt(M0)
  sum sv_k^2 = 6*M0
  p2 = (p1^2 - sum sv_k^2)/2 = (9*M0 - 6*M0)/2 = 3*M0/2
  p3 = prod sv_k = M0^(3/2) * P(delta)

where P(delta) = prod_k (1 + sqrt(2)*cos(2*pi*k/3 + delta)).

Q_dual = 1 - 2*p1*p3/p2^2
       = 1 - 2 * 3*sqrt(M0) * M0^(3/2) * P / (3*M0/2)^2
       = 1 - 2 * 3*M0^2*P / (9*M0^2/4)
       = 1 - 8*P/3

Now expand P(delta) analytically. Using:
  sum_k cos(phi_k + d) = 0  (for phi_k = 2*pi*k/3)
  sum_{j<k} cos(phi_j+d)*cos(phi_k+d) = -3/4
  prod_k cos(phi_k+d) = (1/4)*cos(3*d)  (triple product identity)

  P = prod_k (1 + sqrt(2)*cos(phi_k+d))
    = 1 + sqrt(2)*sum cos(phi_k+d) + 2*sum_{j<k} cos(phi_j+d)cos(phi_k+d)
      + 2*sqrt(2)*prod cos(phi_k+d)
    = 1 + 0 + 2*(-3/4) + 2*sqrt(2)*(1/4)*cos(3d)
    = -1/2 + cos(3*delta)/sqrt(2)

Therefore:

  Q_dual = 1 - 8*(-1/2 + cos(3d)/sqrt(2))/3
         = 1 + 4/3 - 8*cos(3d)/(3*sqrt(2))
         = 7/3 - 4*sqrt(2)*cos(3d)/3

This oscillates with period 2*pi/3 in delta.
Q_dual = 2/3 requires cos(3*delta) = 5*sqrt(2)/8 = 0.88388...

CAVEAT: This formula is valid only when all sv_k > 0, which requires
  1 + sqrt(2)*cos(phi_k + delta) > 0 for all k.
When any sv_k < 0 (which happens for certain delta ranges), the masses
are still positive (m_k = sv_k^2) but the "all positive signs" Q_dual
no longer matches the formula because the inversion 1/m_k does not
carry sign information.
""")

# Verify P(delta) formula
TWO_PI_OVER_3 = 2 * np.pi / 3

def P_func(delta):
    """P = prod_k (1 + sqrt(2)*cos(2*pi*k/3 + delta))"""
    P = 1.0
    for k in range(3):
        P *= (1 + np.sqrt(2) * np.cos(TWO_PI_OVER_3 * k + delta))
    return P

def P_analytic(delta):
    return -0.5 + np.cos(3*delta)/np.sqrt(2)

print("Verification of P = -1/2 + cos(3*delta)/sqrt(2):")
for d in [0, np.pi/6, np.pi/4, np.pi/3, np.pi/2, 2*np.pi/3, np.pi, 5.0]:
    P_n = P_func(d)
    P_a = P_analytic(d)
    print(f"  delta = {d:.4f}: P_numerical = {P_n:+.10f}, "
          f"P_analytic = {P_a:+.10f}, diff = {abs(P_n-P_a):.2e}")

print()

# Check which delta ranges have all sv_k > 0
def all_sv_positive(delta):
    for k in range(3):
        sv = 1 + np.sqrt(2) * np.cos(TWO_PI_OVER_3 * k + delta)
        if sv <= 0:
            return False
    return True

# Scan
n_scan = 10000
deltas = np.linspace(0, 2*np.pi, n_scan, endpoint=False)
valid_mask = np.array([all_sv_positive(d) for d in deltas])

# Find boundaries of valid regions
print("Regions of delta where all sv_k > 0 (formula valid):")
transitions = np.diff(valid_mask.astype(int))
starts = np.where(transitions == 1)[0] + 1
ends = np.where(transitions == -1)[0] + 1
# Handle wrapping
if valid_mask[0]:
    starts = np.concatenate([[0], starts])
if valid_mask[-1]:
    ends = np.concatenate([ends, [n_scan-1]])

for s, e in zip(starts, ends):
    print(f"  delta in [{np.degrees(deltas[s]):.1f}, {np.degrees(deltas[e]):.1f}] deg")

print()

# Verify Q_dual formula only in valid regions
print("Verification of Q_dual = 7/3 - 4*sqrt(2)*cos(3*delta)/3 (valid region only):")
print()

def koide_masses_from_delta(M0, delta):
    """Return (m0, m1, m2) from parametrization, all positive."""
    ms = []
    for k in range(3):
        sv = np.sqrt(M0) * (1 + np.sqrt(2) * np.cos(TWO_PI_OVER_3 * k + delta))
        ms.append(sv**2)
    return ms

M0_test = 1.0
for d in [0.0, 0.15, 0.30, 2.0, 4.0, 4.1, 6.0, 6.1]:
    if not all_sv_positive(d):
        continue
    ms = koide_masses_from_delta(M0_test, d)
    Q_dir = koide_Q(ms)
    Q_dual_num = koide_Q([1/m for m in ms])
    Q_dual_form = 7/3 - 4*np.sqrt(2)*np.cos(3*d)/3
    print(f"  delta = {d:.4f} ({np.degrees(d):7.2f} deg): "
          f"Q_direct = {Q_dir:.10f}, "
          f"Q_dual_num = {Q_dual_num:.10f}, "
          f"Q_dual_form = {Q_dual_form:.10f}, "
          f"diff = {abs(Q_dual_num - Q_dual_form):.2e}")

print()

# Now also verify in invalid regions to show the formula breaks
print("In invalid regions (some sv_k < 0), the formula does not match:")
for d in [np.pi/4, np.pi/2, np.pi, 3.0, 5.0]:
    if all_sv_positive(d):
        continue
    ms = koide_masses_from_delta(M0_test, d)
    Q_dir = koide_Q(ms)
    Q_dual_num = koide_Q([1/m for m in ms])
    Q_dual_form = 7/3 - 4*np.sqrt(2)*np.cos(3*d)/3
    print(f"  delta = {d:.4f} ({np.degrees(d):7.2f} deg): "
          f"Q_direct = {Q_dir:.8f} (!=2/3), "
          f"Q_dual_num = {Q_dual_num:.8f}, "
          f"Q_dual_form = {Q_dual_form:.8f}, "
          f"diff = {abs(Q_dual_num - Q_dual_form):.4f}")

print()
print("When sv_k < 0, the mass m_k = sv_k^2 is still positive, but Q_direct")
print("computed with all-positive signs gives Q != 2/3. The Koide identity")
print("Q = 2/3 holds for the SIGNED sum (sum of sv_k), not for sum of |sv_k|.")
print("The Q_dual formula assumes all sv_k > 0.")
print()

# Q_dual = 2/3 solutions
cos3d_crit = 5*np.sqrt(2)/8
print(f"Q_dual = 2/3 requires cos(3*delta) = 5*sqrt(2)/8 = {cos3d_crit:.10f}")
print()

base = np.arccos(cos3d_crit)
all_sols = []
for n in range(-3, 4):
    for sign in [1, -1]:
        sol = (sign * base + 2*np.pi*n) / 3
        if 0 <= sol < 2*np.pi:
            all_sols.append(sol)
all_sols = sorted(set([round(s, 12) for s in all_sols]))

print(f"All solutions delta in [0, 2*pi):")
for sol in all_sols:
    valid = all_sv_positive(sol)
    ms_sol = koide_masses_from_delta(1.0, sol)
    Q_d = koide_Q([1/m for m in ms_sol])
    Q_dir = koide_Q(ms_sol)
    print(f"  delta = {sol:.8f} rad = {np.degrees(sol):.4f} deg, "
          f"all_sv>0: {valid}, "
          f"Q_direct = {Q_dir:.8f}, "
          f"Q_dual = {Q_d:.10f}")

print()


# ═══════════════════════════════════════════════════════════════════════════
# PART 3: Dual Koide values for various triples
# ═══════════════════════════════════════════════════════════════════════════
print("=" * 72)
print("PART 3: Dual Koide Q_dual = Q(1/m_j) for various triples")
print("=" * 72)
print()

triples_part3 = [
    ("Q_dual(d,s,b) +,+,+",           [1/m_d, 1/m_s, 1/m_b],      [+1, +1, +1]),
    ("Q_dual(-d,s,b) -,+,+",          [1/m_d, 1/m_s, 1/m_b],      [-1, +1, +1]),
    ("Q_dual(u,d,s) +,+,+",           [1/m_u, 1/m_d, 1/m_s],      [+1, +1, +1]),
    ("Q_dual(-u,d,s) -,+,+",          [1/m_u, 1/m_d, 1/m_s],      [-1, +1, +1]),
    ("Q_dual(e,mu,tau) +,+,+",        [1/m_e, 1/m_mu, 1/m_tau],   [+1, +1, +1]),
    ("Q_dual(c,b,t) +,+,+",           [1/m_c, 1/m_b, 1/m_t],      [+1, +1, +1]),
    ("Q_dual(-s,c,b) -,+,+",          [1/m_s, 1/m_c, 1/m_b],      [-1, +1, +1]),
]

for name, inv_masses, signs in triples_part3:
    Q = koide_Q(inv_masses, signs)
    dev_pct = (Q - 2/3) / (2/3) * 100
    print(f"  {name:<32}  Q = {Q:.8f}   (2/3 dev: {dev_pct:+.4f}%)")

print()
print("Key result: Q_dual(d,s,b) with all +++ = 0.66521 is 0.22% from 2/3.")
print("This is the 'dual Koide' observation.")
print()


# ═══════════════════════════════════════════════════════════════════════════
# PART 4: SQCD seesaw
# ═══════════════════════════════════════════════════════════════════════════
print("=" * 72)
print("PART 4: SQCD seesaw M_j = C/m_j")
print("=" * 72)
print()

Lambda = 300.0  # MeV

# For uds sector: C^3 = Lambda^6 * m_u * m_d * m_s
C3_uds = Lambda**6 * m_u * m_d * m_s
C_uds = C3_uds**(1/3)
print(f"SQCD seesaw for (u, d, s):")
print(f"  Lambda = {Lambda} MeV")
print(f"  C^3 = Lambda^6 * m_u * m_d * m_s = {C3_uds:.6e} MeV^3")
print(f"  C   = {C_uds:.6f} MeV")
print()

M_u = C_uds / m_u
M_d = C_uds / m_d
M_s = C_uds / m_s

print(f"  M_u = C/m_u = {M_u:.4f} MeV")
print(f"  M_d = C/m_d = {M_d:.4f} MeV")
print(f"  M_s = C/m_s = {M_s:.4f} MeV")
print()

Q_uds_ppp = koide_Q([M_u, M_d, M_s], [+1, +1, +1])
Q_uds_mpp = koide_Q([M_u, M_d, M_s], [-1, +1, +1])

print(f"  Q(M_u, M_d, M_s) +,+,+ = {Q_uds_ppp:.8f}   (2/3 dev: {(Q_uds_ppp - 2/3)/(2/3)*100:+.4f}%)")
print(f"  Q(M_u, M_d, M_s) -,+,+ = {Q_uds_mpp:.8f}   (2/3 dev: {(Q_uds_mpp - 2/3)/(2/3)*100:+.4f}%)")
print()
print("  Note: Q(C/m_j) is independent of C, so these are the same as Q(1/m_j).")
print()

# For dsb sector
C3_dsb = Lambda**6 * m_d * m_s * m_b
C_dsb = C3_dsb**(1/3)

print(f"SQCD seesaw for (d, s, b):")
print(f"  C'^3 = Lambda^6 * m_d * m_s * m_b = {C3_dsb:.6e} MeV^3")
print(f"  C'   = {C_dsb:.4f} MeV")
print()

M_d2 = C_dsb / m_d
M_s2 = C_dsb / m_s
M_b2 = C_dsb / m_b

print(f"  M_d = C'/m_d = {M_d2:.4f} MeV")
print(f"  M_s = C'/m_s = {M_s2:.4f} MeV")
print(f"  M_b = C'/m_b = {M_b2:.4f} MeV")
print()

Q_dsb_ppp = koide_Q([M_d2, M_s2, M_b2], [+1, +1, +1])
Q_dsb_mpp = koide_Q([M_d2, M_s2, M_b2], [-1, +1, +1])

print(f"  Q(M_d, M_s, M_b) +,+,+ = {Q_dsb_ppp:.8f}   (2/3 dev: {(Q_dsb_ppp - 2/3)/(2/3)*100:+.4f}%)")
print(f"  Q(M_d, M_s, M_b) -,+,+ = {Q_dsb_mpp:.8f}   (2/3 dev: {(Q_dsb_mpp - 2/3)/(2/3)*100:+.4f}%)")
print()
print("  Q is C-independent; the seesaw constant only sets the mass scale, not Q.")
print()

# Also try the scb seesaw
C3_scb = Lambda**6 * m_s * m_c * m_b
C_scb = C3_scb**(1/3)
print(f"SQCD seesaw for (s, c, b):")
print(f"  C_scb^3 = Lambda^6 * m_s * m_c * m_b = {C3_scb:.6e} MeV^3")
print(f"  C_scb   = {C_scb:.4f} MeV")
print()

M_s3 = C_scb / m_s
M_c3 = C_scb / m_c
M_b3 = C_scb / m_b

print(f"  M_s = C_scb/m_s = {M_s3:.4f} MeV")
print(f"  M_c = C_scb/m_c = {M_c3:.4f} MeV")
print(f"  M_b = C_scb/m_b = {M_b3:.4f} MeV")
print()

Q_scb_ppp = koide_Q([M_s3, M_c3, M_b3], [+1, +1, +1])
Q_scb_mpp = koide_Q([M_s3, M_c3, M_b3], [-1, +1, +1])
print(f"  Q(M_s, M_c, M_b) +,+,+ = {Q_scb_ppp:.8f}   (2/3 dev: {(Q_scb_ppp - 2/3)/(2/3)*100:+.4f}%)")
print(f"  Q(-M_s, M_c, M_b) -,+,+ = {Q_scb_mpp:.8f}   (2/3 dev: {(Q_scb_mpp - 2/3)/(2/3)*100:+.4f}%)")
print()

# C' for down-type sector (all three down quarks)
print("What is C' if the seesaw applies to the full down-type sector?")
print()
print("  If C'^3 = Lambda^6 * m_d * m_s * m_b (SQCD with N_f=3 down quarks):")
print(f"    C' = {C_dsb:.4f} MeV  (= {C_dsb/1000:.4f} GeV)")
print(f"    This produces seesaw masses M_j = C'/m_j for j = d, s, b")
print(f"    at scales {M_d2:.0f}, {M_s2:.0f}, {M_b2:.0f} MeV.")
print()


# ═══════════════════════════════════════════════════════════════════════════
# PART 5: Dual Koide vs direct Koide connection
# ═══════════════════════════════════════════════════════════════════════════
print("=" * 72)
print("PART 5: Connection between dual Koide Q(1/m_d,1/m_s,1/m_b)")
print("        and direct Koide Q(-m_s,m_c,m_b)")
print("=" * 72)
print()

Q_direct_scb = koide_Q([m_s, m_c, m_b], [-1, +1, +1])
Q_dual_dsb   = koide_Q([1/m_d, 1/m_s, 1/m_b], [+1, +1, +1])

print(f"  Q_direct(-s, c, b) = {Q_direct_scb:.8f}   (dev from 2/3: {(Q_direct_scb-2/3)/(2/3)*100:+.4f}%)")
print(f"  Q_dual(d, s, b)    = {Q_dual_dsb:.8f}   (dev from 2/3: {(Q_dual_dsb-2/3)/(2/3)*100:+.4f}%)")
print()
print("  Both are near 2/3 but involve different quarks: direct uses (s,c,b),")
print("  dual uses (d,s,b). They are not algebraically identical.")
print()

# Exhaustive sign scan for dual triples
print("  Exhaustive sign scan for Q_dual across many triples:")
print()

dual_triples_all = [
    ("(d, s, b)",    [m_d, m_s, m_b]),
    ("(u, d, s)",    [m_u, m_d, m_s]),
    ("(s, c, b)",    [m_s, m_c, m_b]),
    ("(e, mu, tau)", [m_e, m_mu, m_tau]),
    ("(c, b, t)",    [m_c, m_b, m_t]),
    ("(u, c, t)",    [m_u, m_c, m_t]),
    ("(d, c, b)",    [m_d, m_c, m_b]),
]

sign_combos = [
    ([+1, +1, +1], "+,+,+"),
    ([-1, +1, +1], "-,+,+"),
    ([+1, -1, +1], "+,-,+"),
    ([+1, +1, -1], "+,+,-"),
]

best_duals = []

for tname, masses in dual_triples_all:
    inv_masses = [1/m for m in masses]
    for signs, slabel in sign_combos:
        Q = koide_Q(inv_masses, signs)
        dev = abs(Q - 2/3) / (2/3) * 100
        label = f"Q_dual{tname} {slabel}"
        if dev < 2.0:
            best_duals.append((label, Q, dev))
            marker = " <--"
        else:
            marker = ""
        print(f"    {label:<36}  Q = {Q:.8f}  (|dev|: {dev:.4f}%){marker}")
    print()

print("  Triples with Q_dual within 2% of 2/3:")
for label, Q, dev in sorted(best_duals, key=lambda x: x[2]):
    print(f"    {label:<36}  Q = {Q:.8f}  (|dev|: {dev:.4f}%)")

print()

# Explore seesaw mapping between sectors
print("  Is Q_dual(d,s,b) ~ 2/3 connected to Q(-s,c,b) ~ 2/3 through seesaw?")
print()
print("  The seesaw M_j = C/m_j maps m_j -> 1/m_j (up to scale).")
print("  Q_dual(d,s,b) IS the Koide quotient of seesaw-mapped (d,s,b) masses.")
print("  But Q_direct(-s,c,b) involves m_c (not m_d). The connection would")
print("  require a cross-sector seesaw mapping down-type to up-type quarks.")
print()
print("  Under such a cross-type seesaw of (s,c,b): Q(1/s, 1/c, 1/b).")

Q_dual_scb = koide_Q([1/m_s, 1/m_c, 1/m_b], [+1, +1, +1])
Q_dual_scb_m = koide_Q([1/m_s, 1/m_c, 1/m_b], [-1, +1, +1])
print(f"    Q_dual(s,c,b) +,+,+  = {Q_dual_scb:.8f}  (dev: {(Q_dual_scb-2/3)/(2/3)*100:+.4f}%)")
print(f"    Q_dual(-s,c,b) -,+,+ = {Q_dual_scb_m:.8f}  (dev: {(Q_dual_scb_m-2/3)/(2/3)*100:+.4f}%)")
print()
print("  Neither is close to 2/3. The dual Koide near-hit is specific to (d,s,b)+++.")
print()
print("  The dual Koide Q(1/m_d, 1/m_s, 1/m_b) = 2/3 is therefore an independent")
print("  observation, not a consequence of Q(-s, c, b) = 2/3 via seesaw.")
print()


# ═══════════════════════════════════════════════════════════════════════════
# PART 6: Full parametric analysis
# ═══════════════════════════════════════════════════════════════════════════
print("=" * 72)
print("PART 6: Parametric analysis — does Q=2/3 imply Q_dual=2/3?")
print("=" * 72)
print()

print("""On the Koide manifold (Q = 2/3 identically), we derived:

  Q_dual(delta) = 7/3 - 4*sqrt(2)*cos(3*delta)/3

This holds for delta values where all signed roots sv_k > 0.
""")

# Q_dual range
Q_dual_max = 7/3 + 4*np.sqrt(2)/3
Q_dual_min = 7/3 - 4*np.sqrt(2)/3
print(f"  Q_dual range: [{Q_dual_min:.6f}, {Q_dual_max:.6f}]")
print(f"                 = [7/3 - 4*sqrt(2)/3, 7/3 + 4*sqrt(2)/3]")
print(f"                 = [{7/3 - 4*np.sqrt(2)/3:.6f}, {7/3 + 4*np.sqrt(2)/3:.6f}]")
print()

# But the actually achievable range is limited by the all-sv-positive constraint
# sv_k = 1 + sqrt(2)*cos(phi_k + d) > 0  =>  cos(phi_k + d) > -1/sqrt(2)
# This means phi_k + d is NOT in the range (3pi/4, 5pi/4) mod 2pi
# The excluded delta regions are where any sv_k hits zero:
# 1 + sqrt(2)*cos(phi_k + d) = 0  =>  cos(phi_k + d) = -1/sqrt(2)
# => phi_k + d = 3pi/4 or 5pi/4 (mod 2pi)
# For k=0: d = 3pi/4 or 5pi/4
# For k=1: d = 3pi/4 - 2pi/3 or 5pi/4 - 2pi/3
# For k=2: d = 3pi/4 - 4pi/3 or 5pi/4 - 4pi/3

print("  Boundaries of the all-sv-positive region (sv_k = 0):")
boundaries = []
for k in range(3):
    for angle in [3*np.pi/4, 5*np.pi/4]:
        d = (angle - TWO_PI_OVER_3 * k) % (2*np.pi)
        boundaries.append((d, k, angle))
        print(f"    k={k}: delta = {d:.6f} rad = {np.degrees(d):.2f} deg "
              f"(phi_k+d = {np.degrees(angle):.0f} deg)")

boundaries.sort()
print()

# The formula Q_dual = 7/3 - 4sqrt2 cos(3d)/3 is valid between boundaries
# Let's compute the actual achievable Q_dual range in valid regions

valid_deltas = [d for d in deltas if all_sv_positive(d)]
if valid_deltas:
    Q_dual_valid = [7/3 - 4*np.sqrt(2)*np.cos(3*d)/3 for d in valid_deltas]
    print(f"  Achievable Q_dual in all-sv-positive regions:")
    print(f"    min = {min(Q_dual_valid):.8f}")
    print(f"    max = {max(Q_dual_valid):.8f}")
    print(f"    Q_dual = 2/3 is {'achievable' if min(Q_dual_valid) <= 2/3 <= max(Q_dual_valid) else 'NOT achievable'}.")
    print()

# Check physical triples' delta against the dual formula
print("  Physical triples and their Q_dual values:")
print()

def extract_delta(signed_roots):
    """Extract (M0, delta) from signed roots."""
    s = np.array(signed_roots, dtype=float)
    sqrtM0 = s.sum() / 3.0
    M0 = sqrtM0**2
    phi = np.array([TWO_PI_OVER_3 * k for k in range(3)])
    A = np.dot(s, np.cos(phi))
    B = np.dot(s, -np.sin(phi))
    delta = np.arctan2(B, A) % (2*np.pi)
    return M0, delta

phys_triples = [
    ("(e, mu, tau) +++",   [np.sqrt(m_e), np.sqrt(m_mu), np.sqrt(m_tau)],
     [m_e, m_mu, m_tau]),
    ("(-s, c, b) -++",     [-np.sqrt(m_s), np.sqrt(m_c), np.sqrt(m_b)],
     [m_s, m_c, m_b]),
    ("(c, b, t) +++",      [np.sqrt(m_c), np.sqrt(m_b), np.sqrt(m_t)],
     [m_c, m_b, m_t]),
]

for name, sr, masses in phys_triples:
    M0, delta = extract_delta(sr)
    Q_dir = koide_Q(masses, [np.sign(s) for s in sr])
    Q_dual_num = koide_Q([1/m for m in masses])
    cos3d = np.cos(3*delta)
    Q_dual_form = 7/3 - 4*np.sqrt(2)*cos3d/3
    valid = all_sv_positive(delta)

    print(f"  {name}:")
    print(f"    delta = {delta:.8f} rad = {np.degrees(delta):.4f} deg")
    print(f"    cos(3*delta) = {cos3d:.8f}")
    print(f"    Q_direct (signed) = {Q_dir:.8f}")
    print(f"    Q_dual (numerical, all +) = {Q_dual_num:.8f}  "
          f"(dev from 2/3: {(Q_dual_num-2/3)/(2/3)*100:+.4f}%)")
    print(f"    Q_dual (formula) = {Q_dual_form:.8f}  "
          f"(valid={valid}, formula assumes all sv>0 and Q=2/3)")
    print()


# ═══════════════════════════════════════════════════════════════════════════
# Summary
# ═══════════════════════════════════════════════════════════════════════════
print("=" * 72)
print("SUMMARY")
print("=" * 72)
print()

print("1. DIRECT KOIDE Q VALUES:")
print()
for name, masses, signs in triples_part1:
    Q = koide_Q(masses, signs)
    print(f"   {name:<28}  Q = {Q:.8f}")
print()

print("2. ALGEBRAIC RESULT:")
print("   Q = 2/3 does NOT imply Q_dual = 2/3.")
print("   In symmetric polynomial notation:")
print("     Q = 1 - 2*p2/p1^2         (p2/p1^2 = 1/6 for Q = 2/3)")
print("     Q_dual = 1 - 2*p1*p3/p2^2 (p1*p3/p2^2 = 1/6 for Q_dual = 2/3)")
print("   These are independent conditions on different symmetric functions.")
print()
print("   On the Koide manifold (Q = 2/3 exactly):")
print("     Q_dual = 7/3 - 4*sqrt(2)*cos(3*delta)/3")
print("     P(delta) = -1/2 + cos(3*delta)/sqrt(2)")
print(f"     Q_dual = 2/3 requires cos(3*delta) = 5*sqrt(2)/8 = {cos3d_crit:.8f}")
print()

print("3. DUAL KOIDE VALUES (Q of inverted masses, all positive signs):")
print()
for name, inv_masses, signs in triples_part3:
    Q = koide_Q(inv_masses, signs)
    dev = (Q - 2/3)/(2/3)*100
    print(f"   {name:<32}  Q = {Q:.8f}  ({dev:+.4f}%)")
print()

print("4. SQCD SEESAW:")
print("   Q(C/m_j) = Q(1/m_j) for any C > 0. The seesaw constant C only")
print("   sets the mass scale of the dual spectrum, not the Koide quotient.")
print(f"   For (u,d,s) with Lambda = {Lambda} MeV: C = {C_uds:.1f} MeV")
print(f"   For (d,s,b) with Lambda = {Lambda} MeV: C' = {C_dsb:.1f} MeV")
print()

print("5. DUAL vs DIRECT KOIDE:")
print(f"   Q(-s,c,b)       = {Q_direct_scb:.8f}  (direct,  {(Q_direct_scb-2/3)/(2/3)*100:+.4f}% from 2/3)")
print(f"   Q(1/d,1/s,1/b)  = {Q_dual_dsb:.8f}  (dual,    {(Q_dual_dsb-2/3)/(2/3)*100:+.4f}% from 2/3)")
print("   These involve DIFFERENT quarks and are algebraically independent.")
print("   Among all triples and sign combinations scanned, (d,s,b)+++ is the")
print("   ONLY dual triple within 2% of 2/3.")
print()

print("6. PARAMETRIZATION CONCLUSION:")
print("   On the Koide manifold, Q_dual oscillates as 7/3 - 4*sqrt(2)*cos(3d)/3.")
print(f"   Range: [{Q_dual_min:.4f}, {Q_dual_max:.4f}].")
print("   Q_dual = 2/3 is achievable but only at special delta values.")
print("   Q = 2/3 does NOT imply Q_dual = 2/3.")
print()
print("Done.")
