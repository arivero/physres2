#!/usr/bin/env python3
"""
Bloom correction, tau mass, and V_cb: are they connected?

The O'Raifeartaigh seed at t=sqrt(3) gives lepton mass ratio
m_+/m_- = (2+sqrt3)^2 = 13.93. Physical m_tau/m_mu = 16.82.
The bloom (R-breaking perturbation) rotates the Koide phase delta
by 2.3 degrees, creating the electron mass and shifting the
tau/muon ratio.

Key question: does the SAME bloom that fixes the tau mass
also determine V_cb through the bion eigenvalue mixing?
"""
import numpy as np
from scipy.optimize import minimize

# Physical masses (MeV)
m_e = 0.511
m_mu = 105.658
m_tau = 1776.86

m_u = 2.16
m_d = 4.70
m_s = 93.4
m_c = 1275.0
m_b = 4180.0
m_t = 172760.0

f_pi = 92.2
Lambda = 300.0

print("=" * 70)
print("BLOOM → TAU MASS → V_cb CONNECTION")
print("=" * 70)

# =====================================================================
# 1. Lepton bloom: from seed to physical
# =====================================================================
print("\n--- 1. Lepton Koide parametrisation ---\n")

# sqrt(m_k) = v0 * (1 + sqrt(2) * cos(delta + 2*pi*k/3))
v0_lep = (np.sqrt(m_e) + np.sqrt(m_mu) + np.sqrt(m_tau)) / 3
print(f"v0 = {v0_lep:.4f} MeV^{{1/2}}")

# Find physical delta
def koide_residual(delta, v0, targets):
    vals = [(v0 * (1 + np.sqrt(2)*np.cos(delta + 2*np.pi*k/3)))**2
            for k in range(3)]
    return sum((v - t)**2 for v, t in zip(sorted(vals), sorted(targets)))

best = min(np.linspace(0, 2*np.pi, 10000),
           key=lambda d: koide_residual(d, v0_lep, [m_e, m_mu, m_tau]))
res = minimize(koide_residual, best, args=(v0_lep, [m_e, m_mu, m_tau]),
               method='Nelder-Mead')
delta_lep = res.x[0] % (2*np.pi)

# Find nearest seed (where one eigenvalue = 0)
seeds = []
for k0 in range(3):
    for angle in [3*np.pi/4, 5*np.pi/4]:
        seeds.append((angle - 2*np.pi*k0/3) % (2*np.pi))

diffs = [((delta_lep - s + np.pi) % (2*np.pi)) - np.pi for s in seeds]
best_idx = min(range(len(diffs)), key=lambda i: abs(diffs[i]))
delta_seed = seeds[best_idx]
bloom_shift = diffs[best_idx]

print(f"delta_physical = {delta_lep:.4f} rad = {delta_lep*180/np.pi:.2f}°")
print(f"delta_seed     = {delta_seed:.4f} rad = {delta_seed*180/np.pi:.2f}°")
print(f"bloom shift    = {bloom_shift:.4f} rad = {bloom_shift*180/np.pi:.2f}°")

# Seed spectrum
vals_seed = [(v0_lep * (1 + np.sqrt(2)*np.cos(delta_seed + 2*np.pi*k/3)))**2
             for k in range(3)]
vals_seed_s = sorted(vals_seed)
print(f"\nSeed masses:     {vals_seed_s[0]:.3f}, {vals_seed_s[1]:.3f}, {vals_seed_s[2]:.3f} MeV")
print(f"Physical masses: {m_e:.3f}, {m_mu:.3f}, {m_tau:.3f} MeV")
if vals_seed_s[1] > 0:
    print(f"Seed ratio m+/m- = {vals_seed_s[2]/vals_seed_s[1]:.4f}")
print(f"Phys ratio m_tau/m_mu = {m_tau/m_mu:.4f}")
print(f"(2+√3)² = {(2+np.sqrt(3))**2:.4f}")

# =====================================================================
# 2. Quark bloom: same structure?
# =====================================================================
print("\n--- 2. Quark Koide parametrisation ---\n")

# The (-s, c, b) triple has Q ≈ 2/3
# sqrt(m_k) values with signs from the triple (-s, c, b):
z_s = -np.sqrt(m_s)  # negative for s
z_c = np.sqrt(m_c)
z_b = np.sqrt(m_b)

v0_q = (z_s + z_c + z_b) / 3
print(f"(-s, c, b) triple:")
print(f"  z_s = {z_s:.4f}, z_c = {z_c:.4f}, z_b = {z_b:.4f}")
print(f"  v0 = {v0_q:.4f} MeV^{{1/2}}")

Q_scb = (m_s + m_c + m_b) / (z_s + z_c + z_b)**2
print(f"  Q = {Q_scb:.6f} (should be 2/3 = {2/3:.6f})")
print(f"  Deviation from 2/3: {(Q_scb - 2/3)/(2/3)*100:.3f}%")

# Find delta for quarks
def koide_res_signed(delta, v0, targets_signed):
    """targets_signed = list of (sign, mass) pairs."""
    vals = [v0 * (1 + np.sqrt(2)*np.cos(delta + 2*np.pi*k/3))
            for k in range(3)]
    # Match by sorting absolute values
    vals_abs = sorted(enumerate(vals), key=lambda x: abs(x[1]))
    targets_abs = sorted(enumerate(targets_signed), key=lambda x: abs(x[1]))
    return sum((vals_abs[i][1] - targets_abs[i][1])**2 for i in range(3))

targets_q = [z_s, z_c, z_b]  # These are sqrt(m) with signs

best_q = min(np.linspace(0, 2*np.pi, 10000),
             key=lambda d: koide_res_signed(d, v0_q, targets_q))
res_q = minimize(koide_res_signed, best_q, args=(v0_q, targets_q),
                 method='Nelder-Mead')
delta_q = res_q.x[0] % (2*np.pi)

print(f"\n  delta_quark = {delta_q:.4f} rad = {delta_q*180/np.pi:.2f}°")

# Quark seed
seeds_q = []
for k0 in range(3):
    for angle in [3*np.pi/4, 5*np.pi/4]:
        seeds_q.append((angle - 2*np.pi*k0/3) % (2*np.pi))

diffs_q = [((delta_q - s + np.pi) % (2*np.pi)) - np.pi for s in seeds_q]
best_idx_q = min(range(len(diffs_q)), key=lambda i: abs(diffs_q[i]))
delta_seed_q = seeds_q[best_idx_q]
bloom_shift_q = diffs_q[best_idx_q]

print(f"  delta_seed = {delta_seed_q:.4f} rad = {delta_seed_q*180/np.pi:.2f}°")
print(f"  bloom shift = {bloom_shift_q:.4f} rad = {bloom_shift_q*180/np.pi:.2f}°")

# =====================================================================
# 3. Compare lepton and quark bloom shifts
# =====================================================================
print("\n--- 3. Bloom shift comparison ---\n")

print(f"Lepton bloom: Δδ = {bloom_shift:.4f} rad = {bloom_shift*180/np.pi:.2f}°")
print(f"Quark bloom:  Δδ = {bloom_shift_q:.4f} rad = {bloom_shift_q*180/np.pi:.2f}°")
print(f"Ratio: {bloom_shift/bloom_shift_q:.4f}")
print()
print("If both blooms come from the SAME mechanism (ε Φ₀Φ₂),")
print("their ratio should be related to the scale ratio between")
print("the lepton and quark O'R mass parameters.")

# =====================================================================
# 4. v0-doubling and V_cb
# =====================================================================
print("\n--- 4. v0-doubling structure and V_cb ---\n")

# The v0-doubling formula: sqrt(m_b) = 3*sqrt(m_s) + sqrt(m_c)
# In the Koide parametrisation for (-s, c, b):
# z_k = v0 * (1 + sqrt(2) cos(delta + 2pi*k/3))
# At the seed: z_0 = 0 (the zero eigenvalue, which blooms into s)
# After bloom: z_0 = -sqrt(m_s) (with the sign from the triple)
#
# The v0-doubling means: v0(bloom) = 2 * v0(seed)
# At the seed: v0_seed has only two nonzero z's (c and b at m_- and m_+)
# After doubling: v0_bloom = 2 * v0_seed, and the third z appears

# At the seed (before bloom), the two nonzero eigenvalues are:
# z_1 = v0_seed * (1 + sqrt(2)*cos(delta_seed + 2pi/3))
# z_2 = v0_seed * (1 + sqrt(2)*cos(delta_seed + 4pi/3))
# Their ratio should be (2+sqrt3)/(2-sqrt3) = (2+sqrt3)^2

# The v0-doubling coefficient "3" in sqrt(m_b) = 3*sqrt(m_s) + sqrt(m_c)
# means: the s contribution is 3× its "natural" weight.
# This "3" could come from the NUMBER OF COLORS N_c = 3.

coeff_s = 3 * np.sqrt(m_s)
coeff_c = np.sqrt(m_c)
sqrt_mb_pred = coeff_s + coeff_c
mb_pred = sqrt_mb_pred**2

print(f"v0-doubling: √m_b = 3√m_s + √m_c = {sqrt_mb_pred:.3f}")
print(f"  m_b(pred) = {mb_pred:.0f} MeV (PDG: {m_b:.0f}, dev: {(mb_pred-m_b)/m_b*100:.2f}%)")
print()

# The b quark in terms of the Koide eigenstates:
# At the seed, the eigenstates are:
#   |0> (massless, becomes s after bloom)
#   |-> (light: c or s depending on assignment)
#   |+> (heavy: b at seed, then modified by bloom)
#
# After bloom, the physical b is NOT purely |+>.
# The bloom mixes |0> into |+> and |->. The admixture is
# proportional to the bloom angle.
#
# If |b> = cos(alpha)|+> + sin(alpha)|->  (simplified 2-state mixing),
# then V_cb ~ sin(alpha).

# The bloom angle for quarks is bloom_shift_q ≈ some value.
# The mixing angle between the ISS eigenstates and the physical
# mass eigenstates is related to the bloom.

# But the CKM mixes UP-type with DOWN-type. The quark triple (-s,c,b)
# mixes s (down), c (up), b (down). The Koide rotation is WITHIN
# the (-s,c,b) system, mixing these three quarks.

# The CKM V_cb is the overlap between the state "paired with c by EW"
# and the physical b. In the Koide parametrisation, c and b are
# at different Koide eigenvalues (|-> and |+>). Their mixing comes
# from the bloom rotation.

# At leading order, the bloom angle gives a mixing:
# V_cb ~ |bloom_shift_q| * (geometrical factor)

# What geometrical factor? In a 3-state system with bloom rotation
# by angle dd, the mixing between the m_+ and m_- eigenstates is:
# V_{+-} ~ dd * d(z_+)/d(delta) / (z_+ - z_-)

# Let me compute this numerically.
dd = abs(bloom_shift_q)
# At the seed, the three eigenvalues are:
# z_k = v0 * (1 + sqrt(2)*cos(delta_seed + 2*pi*k/3))
# dz_k/d(delta) = -v0 * sqrt(2) * sin(delta_seed + 2*pi*k/3)

# First, identify which k corresponds to which quark at the physical point
vals_phys = [v0_q * (1 + np.sqrt(2)*np.cos(delta_q + 2*np.pi*k/3))
             for k in range(3)]
# Sort and match to (-sqrt(m_s), sqrt(m_c), sqrt(m_b))
order = sorted(range(3), key=lambda k: vals_phys[k])
# order[0] → most negative → -sqrt(m_s)
# order[1] → middle → sqrt(m_c) or similar
# order[2] → largest → sqrt(m_b)

print("Eigenvalue assignment at physical delta:")
names_sorted = ['-√m_s', '√m_c', '√m_b']
for i, name in enumerate(names_sorted):
    k = order[i]
    print(f"  k={k}: z = {vals_phys[k]:+.3f} → {name} = {targets_q[i]:+.3f}")

# Now compute the derivative at the seed
k_b = order[2]  # k index for the b quark
k_c = order[1]  # k index for the c quark
k_s = order[0]  # k index for the s quark

# At the seed:
z_seed = [v0_q * (1 + np.sqrt(2)*np.cos(delta_seed_q + 2*np.pi*k/3))
          for k in range(3)]
dz_dd = [-v0_q * np.sqrt(2) * np.sin(delta_seed_q + 2*np.pi*k/3)
         for k in range(3)]

print(f"\nAt the seed (delta = {delta_seed_q*180/np.pi:.1f}°):")
for k in range(3):
    print(f"  k={k}: z = {z_seed[k]:+.4f}, dz/dδ = {dz_dd[k]:+.4f}")

# The mixing between b and c eigenstates due to bloom:
# In perturbation theory, the off-diagonal matrix element is
# proportional to dd * dz/dd, and the mixing angle is:
# theta_mix ~ dd * |dz_c/dd| / |z_b - z_c| (at the seed)

z_b_seed = z_seed[k_b]
z_c_seed = z_seed[k_c]
dz_c_seed = dz_dd[k_c]
dz_b_seed = dz_dd[k_b]

print(f"\n  z_b(seed) = {z_b_seed:.4f}")
print(f"  z_c(seed) = {z_c_seed:.4f}")
print(f"  z_b - z_c = {z_b_seed - z_c_seed:.4f}")
print(f"  dz_c/dδ = {dz_c_seed:.4f}")
print(f"  dz_b/dδ = {dz_b_seed:.4f}")

# The bloom doesn't directly mix b and c (they're in different
# weak sectors). But it changes the Koide angle, which modifies
# how much "c character" the b state has.
# The relevant quantity is: how much does the b eigenstate change
# its overlap with the c direction when delta rotates by dd?

# More directly: the CKM V_cb is related to the MISMATCH between
# the up-type and down-type diagonalisation. In the O'R framework:
# - Down-type mass matrix diag(m_d, m_s, m_b) is diagonalised at delta_down
# - Up-type mass matrix diag(m_u, m_c, m_t) is diagonalised at delta_up
# - CKM = U_up† U_down

# For the (-s,c,b) triple, s is down-type, c is up-type, b is down-type.
# The triple MIXES up and down quarks — this is the cross-generation structure.
# The Koide rotation acts on all three simultaneously.
# V_cb comes from the angle between the c and b directions in this space.

# Let me try a simpler approach: the Fritzsch-like estimate
# V_cb ~ sqrt(m_s/m_b) - e^{i*phi} * sqrt(m_c/m_t)
# We already know this gives 0.064 (too large by 50%).
# But if the bloom modifies the effective mass ratios...

V_cb_fritzsch_min = abs(np.sqrt(m_s/m_b) - np.sqrt(m_c/m_t))
V_cb_PDG = 0.0422

print(f"\n  Fritzsch V_cb = |√(m_s/m_b) - √(m_c/m_t)|")
print(f"    = |{np.sqrt(m_s/m_b):.4f} - {np.sqrt(m_c/m_t):.4f}|")
print(f"    = {V_cb_fritzsch_min:.4f}")
print(f"  PDG: {V_cb_PDG:.4f}")
print(f"  Overshoot: {(V_cb_fritzsch_min - V_cb_PDG)/V_cb_PDG * 100:.0f}%")

# Now try using the SEED mass ratios instead of physical ones.
# At the seed, the (-s,c,b) triple has m_s_seed ≠ m_s_phys.
# The seed puts s at the zero eigenvalue, so m_s(seed) → 0!
# Then V_cb(seed) = sqrt(m_s_seed/m_b_seed) - sqrt(m_c/m_t)
#                  → 0 - sqrt(m_c/m_t) = sqrt(m_c/m_t) = 0.086
# That's WORSE. But at intermediate bloom...

print(f"\n  At the SEED (m_s → 0):")
print(f"    V_cb(seed) → √(m_c/m_t) = {np.sqrt(m_c/m_t):.4f}")
print(f"    (even worse — the seed has no s-quark mass to cancel)")

# The key insight: the Fritzsch formula has a PHASE e^{i*phi}.
# V_cb = |sqrt(m_s/m_b) - e^{i*phi} * sqrt(m_c/m_t)|
# For phi = 0: V_cb = 0.064 (real, minimum)
# For phi = pi: V_cb = 0.236 (maximum)
# To get V_cb = 0.042, we need a specific phase:
# |0.1496 - e^{i*phi} * 0.0859| = 0.0422
# Let x = 0.1496, y = 0.0859, V = 0.0422
# V^2 = x^2 + y^2 - 2xy cos(phi)
# cos(phi) = (x^2 + y^2 - V^2)/(2xy)

x = np.sqrt(m_s/m_b)
y = np.sqrt(m_c/m_t)
V = V_cb_PDG
cos_phi = (x**2 + y**2 - V**2) / (2*x*y)
if abs(cos_phi) <= 1:
    phi = np.arccos(cos_phi)
    print(f"\n  Fritzsch phase required for V_cb = {V:.4f}:")
    print(f"    cos(φ) = {cos_phi:.4f}")
    print(f"    φ = {phi:.4f} rad = {phi*180/np.pi:.1f}°")
    print(f"    This is a NON-TRIVIAL CP phase.")
else:
    print(f"\n  cos(φ) = {cos_phi:.4f} — outside [-1,1], no solution!")

# =====================================================================
# 5. The bloom angle as a V_cb proxy
# =====================================================================
print("\n--- 5. Bloom angle vs V_cb ---\n")

# Alternative approach: if V_cb measures how much the bloom
# has mixed the ISS eigenstates, then V_cb ~ sin(bloom_shift_q).
# Or V_cb ~ bloom_shift_q (for small angles).

V_cb_from_bloom = abs(bloom_shift_q)
print(f"Quark bloom shift: |Δδ| = {abs(bloom_shift_q):.4f} rad")
print(f"sin(Δδ) = {np.sin(abs(bloom_shift_q)):.4f}")
print(f"V_cb (PDG) = {V_cb_PDG:.4f}")
print(f"Ratio V_cb / |Δδ| = {V_cb_PDG / abs(bloom_shift_q):.4f}")
print(f"Ratio V_cb / sin(Δδ) = {V_cb_PDG / np.sin(abs(bloom_shift_q)):.4f}")

# Also check: V_us (Cabibbo) vs lepton bloom shift
V_us_PDG = 0.2250
print(f"\nLepton bloom shift: |Δδ| = {abs(bloom_shift):.4f} rad")
print(f"sin(Δδ) = {np.sin(abs(bloom_shift)):.4f}")
print(f"V_us (PDG) = {V_us_PDG:.4f}")
print(f"Ratio V_us / |Δδ| = {V_us_PDG / abs(bloom_shift):.4f}")

# =====================================================================
# 6. GST from bloom
# =====================================================================
print("\n--- 6. GST relation in Koide language ---\n")

# GST: sin(theta_C) = sqrt(m_d/m_s) = 0.2243
# Cabibbo angle = 0.2250 (PDG)
# In the Koide parametrisation, m_d/m_s depends on the bloom angle.
# At the seed, m_d = 0 (zero eigenvalue) → sin(theta_C) = 0.
# After bloom by dd, m_d appears, and sin(theta_C) ~ sqrt(m_d/m_s).
# The bloom CREATES the Cabibbo angle from zero.

print("The GST relation sin θ_C = √(m_d/m_s) connects the")
print("Cabibbo angle to the down-type bloom.")
print()

# For the (u,d,s) or (d,s) sector:
# At the isospin seed: m_d = 0 (all mass in m_s)
# After bloom: m_d appears, and θ_C = arcsin(√(m_d/m_s))
#
# The bloom shift in the down-type sector creates m_d from zero.
# This is the SAME mechanism as the electron mass creation:
# a rotation of delta away from 3pi/4.

# For the down-type sector (d, s) at the isospin seed:
# Only m_s is nonzero. m_d = 0.
# After bloom: m_d = (m_u + m_d) * something.
# GST says: θ_C = arcsin(√(m_d/m_s))

sin_C = np.sqrt(m_d/m_s)
theta_C = np.arcsin(sin_C)
print(f"sin θ_C = √(m_d/m_s) = √({m_d}/{m_s}) = {sin_C:.4f}")
print(f"θ_C = {theta_C:.4f} rad = {theta_C*180/np.pi:.2f}°")
print(f"PDG: sin θ_C = {V_us_PDG:.4f}")
print()

# Now the KEY question: is there a relation between the
# lepton bloom angle and the Cabibbo angle?
print("Is the lepton bloom related to the Cabibbo angle?")
print(f"  Lepton bloom: {abs(bloom_shift)*180/np.pi:.2f}°")
print(f"  Cabibbo angle: {theta_C*180/np.pi:.2f}°")
print(f"  Ratio: {abs(bloom_shift)/theta_C:.4f}")
print()

# And for the quark bloom vs V_cb:
print("Is the quark bloom related to V_cb?")
theta_cb = np.arcsin(V_cb_PDG)
print(f"  Quark bloom: {abs(bloom_shift_q)*180/np.pi:.2f}°")
print(f"  θ_cb = arcsin(V_cb) = {theta_cb*180/np.pi:.2f}°")
print(f"  Ratio: {abs(bloom_shift_q)/theta_cb:.4f}")

# =====================================================================
# 7. Summary
# =====================================================================
print("\n" + "=" * 70)
print("SUMMARY")
print("=" * 70)
print()
print("1. Lepton bloom: Δδ = 2.3° rotates the Koide phase,")
print("   creating m_e and shifting m_τ/m_μ from 13.93 to 16.82.")
print()
print(f"2. Quark (-s,c,b) bloom: Δδ = {abs(bloom_shift_q)*180/np.pi:.1f}°")
print()
if abs(cos_phi) <= 1:
    print(f"3. Fritzsch V_cb requires CP phase φ = {phi*180/np.pi:.0f}° to match PDG.")
else:
    print(f"3. Fritzsch V_cb CANNOT match PDG even with CP phase (cos φ = {cos_phi:.4f} > 1).")
print(f"   Without the phase: V_cb = {V_cb_fritzsch_min:.4f} (50% too large).")
print()
print(f"4. The bloom angle is NOT directly equal to V_cb or θ_C.")
print(f"   But the bloom CREATES m_d from zero, and sin θ_C = √(m_d/m_s).")
print(f"   So the Cabibbo angle is an INDIRECT consequence of the bloom.")
print()
print("5. The connection between bloom and CKM is through the")
print("   MASS RATIOS, not the bloom angle directly:")
print("   - Bloom creates m_d → gives θ_C via GST")
print("   - Bloom modifies m_s, m_c, m_b → changes Fritzsch V_cb")
print("   - The CP phase in V_cb is a SEPARATE parameter (not from bloom)")
