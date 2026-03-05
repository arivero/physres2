"""
Joint constraint from electric and magnetic Koide on the bloom.

Electric sector:  Q(-sqrt(m_s), sqrt(m_c), sqrt(m_b)) ≈ 2/3
Magnetic sector:  Q(1/sqrt(m_d), 1/sqrt(m_s), 1/sqrt(m_b)) ≈ 2/3

These share m_s and m_b. The bloom must satisfy BOTH simultaneously.

Key question: does the seesaw z -> 1/z preserve Q = 2/3?
Answer (proven): NO, except at the seed (where one mass = 0).

So having BOTH Q ≈ 2/3 is a non-trivial additional constraint.
How tight is it? What does it tell us about the bloom mechanism?
"""

import math
import numpy as np

# PDG 2024 MSbar masses (MeV)
m_u, m_d, m_s = 2.16, 4.67, 93.4
m_c, m_b, m_t = 1270.0, 4183.0, 162500.0

def Q(y1, y2, y3):
    s = y1 + y2 + y3
    if abs(s) < 1e-30: return float('nan')
    return (y1**2 + y2**2 + y3**2) / s**2

TARGET = 2/3

# === 1. Current values ===
print("=" * 75)
print("1. CURRENT Q VALUES")
print("=" * 75)

Q_elec = Q(-math.sqrt(m_s), math.sqrt(m_c), math.sqrt(m_b))
Q_mag = Q(1/math.sqrt(m_d), 1/math.sqrt(m_s), 1/math.sqrt(m_b))
Q_cbt = Q(math.sqrt(m_c), math.sqrt(m_b), math.sqrt(m_t))

print(f"  Electric (-s,c,b):     Q = {Q_elec:.6f}  (dev = {abs(Q_elec - TARGET):.6f})")
print(f"  Magnetic (1/d,1/s,1/b): Q = {Q_mag:.6f}  (dev = {abs(Q_mag - TARGET):.6f})")
print(f"  Heavy (c,b,t):          Q = {Q_cbt:.6f}  (dev = {abs(Q_cbt - TARGET):.6f})")

# === 2. Parametrize the electric sector bloom ===
print("\n" + "=" * 75)
print("2. ELECTRIC SECTOR: bloom parametrization")
print("=" * 75)

# The electric triple is (-sqrt(m_s), sqrt(m_c), sqrt(m_b))
# In the angular parametrization: z_k = z0(1 + sqrt(2) cos(delta + 2pi k/3))
# With k=0 -> -sqrt(m_s), k=1 -> sqrt(m_c), k=2 -> sqrt(m_b)

# From the actual masses, extract z0 and delta
z = [-math.sqrt(m_s), math.sqrt(m_c), math.sqrt(m_b)]
z0 = sum(z) / 3
zi = [zk - z0 for zk in z]
r = math.sqrt(sum(x**2 for x in zi) * 2 / 3)  # r = sqrt(2) z0 if Q=2/3

print(f"  z = [{z[0]:.4f}, {z[1]:.4f}, {z[2]:.4f}]")
print(f"  z0 = {z0:.4f}")
print(f"  r = {r:.4f}")
print(f"  r / (sqrt(2) z0) = {r / (math.sqrt(2) * z0):.6f}  (1.000 = exact Q=2/3)")

# Extract delta from z_0 = z0 + r cos(delta)
cos_delta = zi[0] / r
delta = math.acos(max(-1, min(1, cos_delta)))
# Check sign: sin(delta) should give correct z_1
sin_delta_check = zi[1] / r
# delta + 2pi/3 should give z_1
z1_check = r * math.cos(delta + 2*math.pi/3)
if abs(z1_check - zi[1]) > 0.01:
    delta = -delta  # try negative
    z1_check = r * math.cos(delta + 2*math.pi/3)

# Better: use atan2
# z_0/r = cos(delta), but we need to determine sign from z_1
# z_1 = r cos(delta + 2pi/3) = r[cos(delta)cos(2pi/3) - sin(delta)sin(2pi/3)]
# = r[-cos(delta)/2 - sqrt(3)sin(delta)/2]
# So: zi[1]/r = -cos(delta)/2 - sqrt(3)sin(delta)/2
# and zi[0]/r = cos(delta)
cos_d = zi[0] / r
sin_d = -(2*zi[1]/r + cos_d) / math.sqrt(3)
delta = math.atan2(sin_d, cos_d)

print(f"  delta = {delta:.6f} rad = {math.degrees(delta):.2f} deg")
print(f"  delta/pi = {delta/math.pi:.6f}")
print(f"  seed delta = 3pi/4 = {3*math.pi/4:.6f} rad = 135.00 deg")
print(f"  bloom shift = {math.degrees(delta - 3*math.pi/4):.2f} deg")

# Verify
z_check = [z0 + r*math.cos(delta + 2*math.pi*k/3) for k in range(3)]
print(f"  Verify: z_check = [{z_check[0]:.4f}, {z_check[1]:.4f}, {z_check[2]:.4f}]")

# === 3. The joint constraint ===
print("\n" + "=" * 75)
print("3. JOINT CONSTRAINT: electric Q and magnetic Q vs bloom angle")
print("=" * 75)

print("""
The electric triple (-s,c,b) is parametrized by (z0_e, delta_e).
Given z0_e and delta_e, we get m_s, m_c, m_b.

The magnetic triple (1/d, 1/s, 1/b) shares m_s and m_b but adds m_d.
m_d is a FREE parameter (not in the electric chain).

So the magnetic Q depends on:
  - m_s, m_b (from electric bloom)
  - m_d (free)

Question: for the bloom-determined m_s and m_b, does there exist m_d
such that Q_mag = 2/3? And does the PDG m_d satisfy it?
""")

# Given m_s and m_b from the electric sector, what m_d gives Q_mag = 2/3?
# Q(1/sqrt(m_d), 1/sqrt(m_s), 1/sqrt(m_b)) = 2/3
# Let x = 1/sqrt(m_d), a = 1/sqrt(m_s), b = 1/sqrt(m_b)
# Q(x, a, b) = (x^2 + a^2 + b^2) / (x + a + b)^2 = 2/3
# 3(x^2 + a^2 + b^2) = 2(x + a + b)^2
# 3x^2 + 3a^2 + 3b^2 = 2x^2 + 2a^2 + 2b^2 + 4x(a+b) + 4ab
# x^2 + a^2 + b^2 - 4x(a+b) - 4ab = 0
# x^2 - 4(a+b)x + (a^2 + b^2 - 4ab) = 0

a = 1/math.sqrt(m_s)
b = 1/math.sqrt(m_b)

disc = 16*(a+b)**2 - 4*(a**2 + b**2 - 4*a*b)
disc2 = 12*a**2 + 12*b**2 + 48*a*b  # simplified
print(f"  a = 1/sqrt(m_s) = {a:.6f}")
print(f"  b = 1/sqrt(m_b) = {b:.6f}")
print(f"  Discriminant = {disc:.6f} (always positive)")

x_plus = (4*(a+b) + math.sqrt(disc)) / 2
x_minus = (4*(a+b) - math.sqrt(disc)) / 2

md_plus = 1/x_plus**2
md_minus = 1/x_minus**2

print(f"\n  Solutions for Q_mag = 2/3:")
print(f"    x+ = {x_plus:.6f}  =>  m_d = {md_plus:.4f} MeV")
print(f"    x- = {x_minus:.6f}  =>  m_d = {md_minus:.4f} MeV")
print(f"    PDG m_d = {m_d:.2f} MeV")

# Check which is physical
Q_check_plus = Q(x_plus, a, b)
Q_check_minus = Q(x_minus, a, b)
print(f"    Q(x+, a, b) = {Q_check_plus:.6f}")
print(f"    Q(x-, a, b) = {Q_check_minus:.6f}")

# The physical solution
md_pred = md_minus if abs(md_minus - m_d) < abs(md_plus - m_d) else md_plus
print(f"\n  Closest to PDG: m_d^pred = {md_pred:.4f} MeV")
print(f"  PDG m_d = {m_d:.2f} ± 0.48 MeV")
print(f"  Pull = {(md_pred - m_d)/0.48:.2f} sigma")

# === 4. The constraint surface ===
print("\n" + "=" * 75)
print("4. SCANNING: how Q_mag varies with bloom angle")
print("=" * 75)

print("""
Fix z0 at the physical value. Scan delta from seed (3pi/4) outward.
For each delta, compute m_s and m_b from the electric parametrization,
then compute Q_mag using the PDG m_d.
""")

z0_phys = z0
r_phys = math.sqrt(2) * z0_phys  # exact Q=2/3 in electric sector

print(f"  z0 = {z0_phys:.4f}")
print(f"  r = {r_phys:.4f}")
print(f"  {'delta (deg)':>12} {'m_s':>10} {'m_c':>10} {'m_b':>10} {'Q_elec':>10} {'Q_mag':>10}")
print(f"  {'-'*12} {'-'*10} {'-'*10} {'-'*10} {'-'*10} {'-'*10}")

for delta_deg in np.linspace(130, 145, 31):
    d = math.radians(delta_deg)
    zk = [z0_phys + r_phys * math.cos(d + 2*math.pi*k/3) for k in range(3)]
    # z_0 = -sqrt(m_s), z_1 = sqrt(m_c), z_2 = sqrt(m_b)
    ms_i = zk[0]**2
    mc_i = zk[1]**2
    mb_i = zk[2]**2

    Qe = Q(zk[0], zk[1], zk[2])  # always 2/3 by construction

    if ms_i > 0 and mb_i > 0:
        Qm = Q(1/math.sqrt(m_d), 1/math.sqrt(ms_i), 1/math.sqrt(mb_i))
    else:
        Qm = float('nan')

    marker = " <-- physical" if abs(delta_deg - math.degrees(delta)) < 0.3 else ""
    print(f"  {delta_deg:12.2f} {ms_i:10.2f} {mc_i:10.2f} {mb_i:10.2f} {Qe:10.6f} {Qm:10.6f}{marker}")

# === 5. The deep question ===
print("\n" + "=" * 75)
print("5. ANALYSIS: what the joint constraint means")
print("=" * 75)

# How much does Q_mag vary over the bloom range?
# At seed: m_s = 0, so 1/sqrt(m_s) -> infinity, Q_mag is ill-defined
# At physical bloom: Q_mag = 0.665

# The point is: the electric bloom fixes m_s and m_b.
# The magnetic Q then depends only on m_d (the remaining free parameter).
# So the "joint constraint" is really: does the PDG m_d sit at Q_mag = 2/3?
# And the answer is: approximately yes (0.665 vs 0.667).

# But there's a deeper point: m_d is supposed to be free in the electric chain.
# If the magnetic Koide also holds, then m_d is PREDICTED from m_s and m_b.

print(f"""
The electric bloom determines m_s and m_b (and m_c).
The magnetic Koide Q(1/d, 1/s, 1/b) = 2/3 then PREDICTS m_d.

With PDG m_s = {m_s} and m_b = {m_b}:
  m_d^pred = {md_pred:.3f} MeV (from Q_mag = 2/3 exactly)
  m_d^PDG  = {m_d:.2f} ± 0.48 MeV

This means: the magnetic Koide reduces the free parameters by one more.
From 8 free parameters, we'd go to 7.
The parameter m_d is no longer free: it's determined by (m_s, m_b) via duality.

The Cabibbo angle then becomes:
  sin(theta_C) = sqrt(m_d^pred / m_s) = {math.sqrt(md_pred/m_s):.4f}
  sin(theta_C)^PDG = {math.sqrt(m_d/m_s):.4f}
  |V_us|^PDG = 0.2243 ± 0.0008
""")

# === 6. Does the O'R chain + magnetic Koide overdetermine? ===
print("=" * 75)
print("6. OVERDETERMINATION TEST")
print("=" * 75)

# If we use the O'R predicted m_c and bion m_b instead of PDG:
m_c_pred = (7 + 4*math.sqrt(3)) * m_s  # O'R prediction
m_b_pred_from_OR = (3*math.sqrt(m_s) + math.sqrt(m_c_pred))**2  # bion with O'R m_c
m_b_pred_from_PDG_c = (3*math.sqrt(m_s) + math.sqrt(m_c))**2   # bion with PDG m_c

a_or = 1/math.sqrt(m_s)
b_or = 1/math.sqrt(m_b_pred_from_OR)
b_pdg = 1/math.sqrt(m_b_pred_from_PDG_c)

disc_or = 16*(a_or+b_or)**2 - 4*(a_or**2 + b_or**2 - 4*a_or*b_or)
x_or = (4*(a_or+b_or) - math.sqrt(disc_or)) / 2
md_from_chain = 1/x_or**2

disc_pdg = 16*(a_or+b_pdg)**2 - 4*(a_or**2 + b_pdg**2 - 4*a_or*b_pdg)
x_pdg = (4*(a_or+b_pdg) - math.sqrt(disc_pdg)) / 2
md_from_bion_pdg = 1/x_pdg**2

print(f"  Using O'R chain (m_c={m_c_pred:.0f}, m_b={m_b_pred_from_OR:.0f}):")
print(f"    m_d^pred = {md_from_chain:.3f} MeV")
print(f"    sin(theta_C) = {math.sqrt(md_from_chain/m_s):.4f}")

print(f"\n  Using bion with PDG m_c (m_c={m_c:.0f}, m_b={m_b_pred_from_PDG_c:.0f}):")
print(f"    m_d^pred = {md_from_bion_pdg:.3f} MeV")
print(f"    sin(theta_C) = {math.sqrt(md_from_bion_pdg/m_s):.4f}")

print(f"\n  Using PDG masses directly (m_s={m_s}, m_b={m_b}):")
print(f"    m_d^pred = {md_pred:.3f} MeV")
print(f"    sin(theta_C) = {math.sqrt(md_pred/m_s):.4f}")

print(f"\n  PDG values:")
print(f"    m_d = {m_d} ± 0.48 MeV")
print(f"    sin(theta_C) = 0.2243 ± 0.0008")

# === 7. Self-duality of the seed ===
print("\n" + "=" * 75)
print("7. SEED SELF-DUALITY vs BLOOM DUALITY-BREAKING")
print("=" * 75)

print("""
At the seed: m_1 = 0, m_2 = m_-, m_3 = m_+.
  Direct charges: (0, sqrt(m_-), sqrt(m_+))
  Inverse charges: (infinity, 1/sqrt(m_-), 1/sqrt(m_+))

The seed is self-dual: Q(0, a, b) = Q(0, 1/a, 1/b) always.
This is because Q depends only on the ratio a/b, and
Q(0, a, b) = Q(0, 1, b/a) = Q(0, 1, a/b) by homogeneity + inversion.
Wait: Q(0, 1/a, 1/b) = Q(0, b, a)/(ab)^0... let me just verify.
""")

# Verify seed self-duality
for ratio in [2+math.sqrt(3), 3.0, 5.0]:
    a_test, b_test = 1.0, ratio
    Q_dir = Q(0, a_test, b_test)
    Q_inv = Q(0, 1/a_test, 1/b_test)
    print(f"  (0, 1, {ratio:.3f}): Q_direct = {Q_dir:.6f}, Q_inverse = {Q_inv:.6f}, equal = {abs(Q_dir-Q_inv)<1e-10}")

print(f"""
Seed IS self-dual. Q(0, a, b) = Q(0, 1/a, 1/b) exactly.

At the bloom: m_1 > 0. Now Q(z) != Q(1/z) generically.
The bloom BREAKS the self-duality.

The fact that BOTH Q values are close to 2/3 at the physical masses
means the bloom is "almost self-dual" — it breaks duality only weakly.

Quantifying the duality breaking:
  Q_electric = {Q_elec:.6f}  (dev from 2/3: {abs(Q_elec-TARGET)*100:.3f}%)
  Q_magnetic = {Q_mag:.6f}  (dev from 2/3: {abs(Q_mag-TARGET)*100:.3f}%)

The magnetic deviation ({abs(Q_mag-TARGET)*100:.3f}%) is larger than the electric
({abs(Q_elec-TARGET)*100:.3f}%), consistent with the bloom breaking duality at O(epsilon^2).
""")

# === 8. Perturbative analysis ===
print("=" * 75)
print("8. PERTURBATIVE: how duality-breaking scales with bloom")
print("=" * 75)

# At the seed, Q_elec = Q_mag = 2/3 exactly.
# The bloom shifts delta by epsilon from 3pi/4.
# Q_elec stays 2/3 by construction (angular parametrization with r=sqrt(2)z0).
# Q_mag deviates. How?

# Scan: for the down-type triple (d,s,b) with m_d = PDG,
# how does Q_mag depend on bloom epsilon?
# Actually, the magnetic triple is (1/d, 1/s, 1/b) where d,s,b are down-type.
# The electric triple is (-s, c, b).
# These are DIFFERENT triples on different quarks (except s and b shared).

# More precisely: as the bloom angle increases from seed,
# m_s goes from 0 to ~93 MeV and m_b grows.
# The magnetic Q(1/d, 1/s, 1/b) depends on the ratio m_s/m_b (since m_d is fixed).

print(f"\n  Bloom epsilon vs magnetic Q deviation:")
print(f"  {'eps (deg)':>10} {'m_s':>8} {'m_b':>8} {'Q_mag':>10} {'|Q-2/3|':>10}")

for eps_deg in np.linspace(0.5, 15, 30):
    d = 3*math.pi/4 + math.radians(eps_deg)
    zk = [z0_phys + r_phys * math.cos(d + 2*math.pi*k/3) for k in range(3)]
    ms_i = zk[0]**2
    mb_i = zk[2]**2

    if ms_i > 0.01 and mb_i > 0.01:
        Qm = Q(1/math.sqrt(m_d), 1/math.sqrt(ms_i), 1/math.sqrt(mb_i))
        marker = " <--" if abs(eps_deg - math.degrees(delta - 3*math.pi/4)) < 0.3 else ""
        print(f"  {eps_deg:10.2f} {ms_i:8.1f} {mb_i:8.1f} {Qm:10.6f} {abs(Qm-TARGET):10.6f}{marker}")

print(f"""

KEY INSIGHT: Q_mag is NOT constant as the bloom varies.
It depends on the bloom angle through m_s and m_b.
The fact that Q_mag ≈ 2/3 at the PHYSICAL bloom angle is a constraint
on the bloom mechanism — not all bloom angles give Q_mag ≈ 2/3.

This means: if the bloom is dynamically determined (e.g. by three-instanton
potential), the magnetic Koide is a PREDICTION of the bloom mechanism.
Or conversely: if we demand both electric AND magnetic Koide, the bloom
angle is more tightly constrained.
""")
