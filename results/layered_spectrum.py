#!/usr/bin/env python3
"""
Layered spectrum computation: parameter reduction via algebraic relations.
Pure mathematics — eigenvalue structure of a superpotential.
"""
import math

# ============================================================
# INPUTS
# ============================================================
# Free inputs (MeV)
m_s = 93.4       # strange quark mass
m_e = 0.511      # electron mass
m_mu = 105.658   # muon mass
m_d = 4.67       # down quark mass

# PDG comparison values
pdg = {
    'm_c': (1275, 25),       # MeV
    'm_b': (4180, 30),       # MeV
    'm_t': (172760, 300),    # MeV (pole mass)
    'm_tau': (1776.86, 0.12), # MeV
    'Vus': (0.2250, 0.0007),  # |V_us|
}

def pull(pred, central, sigma):
    return (pred - central) / sigma

print("=" * 70)
print("LAYERED SPECTRUM — PARAMETER REDUCTION")
print("=" * 70)

# ============================================================
# LEVEL 0: No relations — 13 free parameters
# ============================================================
print("\n--- LEVEL 0: No relations ---")
print("Free parameters: 13")
print("  Masses: m_u, m_d, m_s, m_c, m_b, m_t, m_e, m_mu, m_tau")
print("  Mixing: theta_12, theta_23, theta_13, delta_CP")

# ============================================================
# LEVEL 1: R1 — O'Raifeartaigh mass ratio
# m_c = (2 + sqrt(3))^2 * m_s
# ============================================================
print("\n--- LEVEL 1: R1 — O'Raifeartaigh mass ratio ---")
ratio_R1 = (2 + math.sqrt(3))**2
m_c_pred = ratio_R1 * m_s
print(f"  Ratio (2+√3)² = 7 + 4√3 = {ratio_R1:.10f}")
print(f"  m_c = {ratio_R1:.6f} × {m_s} = {m_c_pred:.2f} MeV")
p_c = pull(m_c_pred, *pdg['m_c'])
print(f"  PDG m_c = {pdg['m_c'][0]} ± {pdg['m_c'][1]} MeV")
print(f"  Pull: ({m_c_pred:.2f} - {pdg['m_c'][0]}) / {pdg['m_c'][1]} = {p_c:.2f}σ")
print(f"  Determined: m_c")
print(f"  Free parameters remaining: 12")

# ============================================================
# LEVEL 2: R2 — Bion mass relation
# sqrt(m_b) = 3*sqrt(m_s) + sqrt(m_c)
# ============================================================
print("\n--- LEVEL 2: R2 — Bion mass relation ---")
sqrt_mb = 3 * math.sqrt(m_s) + math.sqrt(m_c_pred)
m_b_pred = sqrt_mb**2
print(f"  √m_b = 3√m_s + √m_c = 3×{math.sqrt(m_s):.6f} + {math.sqrt(m_c_pred):.6f}")
print(f"       = {sqrt_mb:.6f}")
print(f"  m_b = {m_b_pred:.2f} MeV")
p_b = pull(m_b_pred, *pdg['m_b'])
print(f"  PDG m_b = {pdg['m_b'][0]} ± {pdg['m_b'][1]} MeV")
print(f"  Pull: {p_b:.2f}σ")
print(f"  Determined: m_b")
print(f"  Free parameters remaining: 11")

# ============================================================
# LEVEL 3: R3 — Yukawa eigenvalue constraint (Koide for quarks)
# (m_c + m_b + m_t) / (sqrt(m_c) + sqrt(m_b) + sqrt(m_t))^2 = 2/3
# Solve for m_t
# ============================================================
print("\n--- LEVEL 3: R3 — Yukawa eigenvalue constraint ---")
# Let x = sqrt(m_t). Then:
# m_c + m_b + x^2 = (2/3)(sqrt(m_c) + sqrt(m_b) + x)^2
# Let S = sqrt(m_c) + sqrt(m_b), M = m_c + m_b
# M + x^2 = (2/3)(S + x)^2 = (2/3)(S^2 + 2Sx + x^2)
# M + x^2 = (2/3)S^2 + (4/3)Sx + (2/3)x^2
# (1/3)x^2 - (4/3)Sx + M - (2/3)S^2 = 0
# x^2 - 4Sx + 3M - 2S^2 = 0
# x = (4S ± sqrt(16S^2 - 4(3M - 2S^2))) / 2
# x = 2S ± sqrt(4S^2 - 3M + 2S^2)
# x = 2S ± sqrt(6S^2 - 3M)
# x = 2S ± sqrt(3(2S^2 - M))

S_q = math.sqrt(m_c_pred) + math.sqrt(m_b_pred)
M_q = m_c_pred + m_b_pred
discriminant = 3 * (2 * S_q**2 - M_q)
print(f"  S = √m_c + √m_b = {S_q:.6f}")
print(f"  M = m_c + m_b = {M_q:.2f}")
print(f"  Discriminant = 3(2S² - M) = {discriminant:.4f}")
sqrt_disc = math.sqrt(discriminant)
x_plus = 2 * S_q + sqrt_disc
x_minus = 2 * S_q - sqrt_disc
m_t_plus = x_plus**2
m_t_minus = x_minus**2
print(f"  √m_t = 2S ± √(3(2S²-M))")
print(f"  Solution +: √m_t = {x_plus:.4f}, m_t = {m_t_plus:.2f} MeV")
print(f"  Solution -: √m_t = {x_minus:.4f}, m_t = {m_t_minus:.2f} MeV")
# Physical solution is the large one
m_t_pred = m_t_plus
p_t = pull(m_t_pred, *pdg['m_t'])
print(f"  Physical: m_t = {m_t_pred:.2f} MeV")
print(f"  PDG m_t = {pdg['m_t'][0]} ± {pdg['m_t'][1]} MeV")
print(f"  Pull: {p_t:.2f}σ")

# Verify Koide
Q_quark = (m_c_pred + m_b_pred + m_t_pred) / (math.sqrt(m_c_pred) + math.sqrt(m_b_pred) + math.sqrt(m_t_pred))**2
print(f"  Verification: Q(c,b,t) = {Q_quark:.15f} (should be 2/3 = {2/3:.15f})")
print(f"  Determined: m_t")
print(f"  Free parameters remaining: 10")

# ============================================================
# LEVEL 4: R4 — Lepton sector Koide
# (m_e + m_mu + m_tau) / (sqrt(m_e) + sqrt(m_mu) + sqrt(m_tau))^2 = 2/3
# ============================================================
print("\n--- LEVEL 4: R4 — Lepton sector Koide ---")
S_l = math.sqrt(m_e) + math.sqrt(m_mu)
M_l = m_e + m_mu
disc_l = 3 * (2 * S_l**2 - M_l)
sqrt_disc_l = math.sqrt(disc_l)
x_l_plus = 2 * S_l + sqrt_disc_l
x_l_minus = 2 * S_l - sqrt_disc_l
m_tau_plus = x_l_plus**2
m_tau_minus = x_l_minus**2
print(f"  S = √m_e + √m_μ = {S_l:.6f}")
print(f"  M = m_e + m_μ = {M_l:.6f}")
print(f"  Solution +: √m_τ = {x_l_plus:.6f}, m_τ = {m_tau_plus:.4f} MeV")
print(f"  Solution -: √m_τ = {x_l_minus:.6f}, m_τ = {m_tau_minus:.4f} MeV")
m_tau_pred = m_tau_plus
p_tau = pull(m_tau_pred, *pdg['m_tau'])
print(f"  Physical: m_τ = {m_tau_pred:.4f} MeV")
print(f"  PDG m_τ = {pdg['m_tau'][0]} ± {pdg['m_tau'][1]} MeV")
print(f"  Pull: {p_tau:.2f}σ")

Q_lepton = (m_e + m_mu + m_tau_pred) / (math.sqrt(m_e) + math.sqrt(m_mu) + math.sqrt(m_tau_pred))**2
print(f"  Verification: Q(e,μ,τ) = {Q_lepton:.15f}")
print(f"  Determined: m_τ")
print(f"  Free parameters remaining: 9")

# ============================================================
# LEVEL 5: R5 — Gatto-Sartori-Tonin
# sin(theta_12) = sqrt(m_d / m_s)
# ============================================================
print("\n--- LEVEL 5: R5 — Gatto-Sartori-Tonin ---")
sin_th12 = math.sqrt(m_d / m_s)
Vus_pred = sin_th12
print(f"  sin θ₁₂ = √(m_d/m_s) = √({m_d}/{m_s}) = {sin_th12:.6f}")
print(f"  |V_us| = {Vus_pred:.6f}")
p_Vus = pull(Vus_pred, *pdg['Vus'])
print(f"  PDG |V_us| = {pdg['Vus'][0]} ± {pdg['Vus'][1]}")
print(f"  Pull: {p_Vus:.2f}σ")
print(f"  Determined: θ₁₂")
print(f"  Free parameters remaining: 8")

# ============================================================
# CHARGE REPRESENTATION
# ============================================================
print("\n" + "=" * 70)
print("CHARGE REPRESENTATION")
print("=" * 70)

def compute_charges(label, sqrt_masses, check_koide=True):
    """
    Given sqrt(m_k), compute z_0, z_k, and check <z_k^2> = z_0^2.
    sqrt_masses can include negative entries (signed roots).
    """
    frak_z = list(sqrt_masses)  # these ARE the frak_z
    z0 = sum(frak_z) / 3
    z = [fz - z0 for fz in frak_z]
    
    # Masses
    masses = [fz**2 for fz in frak_z]
    
    # Check sum z_k = 0
    sum_z = sum(z)
    
    # <z_k^2> = (1/3) sum z_k^2
    mean_z2 = sum(zk**2 for zk in z) / 3
    
    # Koide condition: <z_k^2> = z_0^2
    koide_check = mean_z2 / z0**2 if z0 != 0 else float('inf')
    
    # Also compute Q directly
    sum_m = sum(masses)
    sum_sqrt = sum(frak_z)
    Q = sum_m / sum_sqrt**2 if sum_sqrt != 0 else float('inf')
    
    print(f"\n  {label}:")
    print(f"    √m = ({frak_z[0]:.6f}, {frak_z[1]:.6f}, {frak_z[2]:.6f})")
    print(f"    m   = ({masses[0]:.4f}, {masses[1]:.4f}, {masses[2]:.4f})")
    print(f"    z₀  = {z0:.6f}")
    print(f"    z   = ({z[0]:.6f}, {z[1]:.6f}, {z[2]:.6f})")
    print(f"    Σz  = {sum_z:.2e} (should be 0)")
    print(f"    ⟨z²⟩ = {mean_z2:.6f}")
    print(f"    z₀²  = {z0**2:.6f}")
    print(f"    ⟨z²⟩/z₀² = {koide_check:.10f} (should be 1 if Koide holds)")
    print(f"    Q = Σm/(Σ√m)² = {Q:.10f} (should be 2/3 = {2/3:.10f} if Koide)")
    
    return z0, z, masses, Q

# Level 1: Seed triple (0, m_s, m_c)
print("\n--- Level 1: Quark seed triple (0, m_s, m_c) ---")
sqrt_seed = [0, math.sqrt(m_s), math.sqrt(m_c_pred)]
z0_seed, z_seed, m_seed, Q_seed = compute_charges("Seed (0, m_s, m_c)", sqrt_seed)

# Verify the algebraic identity for seed
print(f"\n    Algebraic check: eigenvalues (0, m_-, m_+) with m_±=(2±√3)m")
m_ratio = m_c_pred / m_s
print(f"    m_c/m_s = {m_ratio:.6f} = (2+√3)²/(2-√3)² should be... wait")
print(f"    Actually: m_+ = (2+√3)² × m_s, m_- = m_s (seed is (0, m_s, m_c))")
print(f"    For (0, a, b): Q = (a+b)/(√a+√b)² ")
a_val = m_s
b_val = m_c_pred
Q_check = (a_val + b_val) / (math.sqrt(a_val) + math.sqrt(b_val))**2
print(f"    Q = ({a_val:.4f} + {b_val:.4f}) / ({math.sqrt(a_val):.4f} + {math.sqrt(b_val):.4f})²")
print(f"      = {a_val + b_val:.4f} / {(math.sqrt(a_val) + math.sqrt(b_val))**2:.4f} = {Q_check:.10f}")
# This should NOT be 2/3 for a 2-body seed
# The 2/3 for (0, m_-, m_+): Q = (0 + m_- + m_+)/(0 + √m_- + √m_+)²
# With m_± = (2±√3)²m: sum_m = ((2-√3)² + (2+√3)²)m = (7-4√3 + 7+4√3)m = 14m
# sum_sqrt = (2-√3 + 2+√3)√m = 4√m
# Q = 14m / 16m = 14/16 = 7/8... NO
# Wait: (0, m_-, m_+) with m_- = (2-√3)m, m_+ = (2+√3)m (not squared)
# Then sum_m = 0 + (2-√3)m + (2+√3)m = 4m
# sum_sqrt = 0 + √((2-√3)m) + √((2+√3)m) = (√(2-√3) + √(2+√3))√m
# Now √(2-√3) + √(2+√3): let's compute
v1 = math.sqrt(2 - math.sqrt(3))
v2 = math.sqrt(2 + math.sqrt(3))
print(f"\n    With linear eigenvalues (0, (2-√3)m, (2+√3)m):")
print(f"    √(2-√3) = {v1:.6f}, √(2+√3) = {v2:.6f}")
print(f"    Sum = {v1+v2:.6f}, Sum² = {(v1+v2)**2:.6f}")
print(f"    Q = 4m / ({(v1+v2)**2:.6f} × m) = {4/(v1+v2)**2:.10f}")
# That's 4/(v1+v2)^2. (v1+v2)^2 = 2-√3 + 2+√3 + 2√(4-3) = 4 + 2 = 6
print(f"    (v1+v2)² = (2-√3) + (2+√3) + 2√((2-√3)(2+√3)) = 4 + 2√(4-3) = 4+2 = 6")
print(f"    Q = 4/6 = 2/3 ✓")

# But in our actual setup: m_s is the (2-√3)² scaling of some base, not linear
# The MEMORY says: "eigenvalues (0, m_-, m_+) with m_± = (2±√3)m, Q = 4m/6m = 2/3"
# So the Koide identity holds for the triple (0, (2-√3)m, (2+√3)m) with LINEAR scaling
# In our case the relation is m_c = (2+√3)² × m_s, which means
# if m_s = (2-√3)² × m_base, then m_c = (2+√3)² × m_base ... but we didn't assume that.
# Let me just check: is m_s = (2-√3)²×m_base, m_c = (2+√3)²×m_base for some m_base?
# Then m_c/m_s = ((2+√3)/(2-√3))² = (2+√3)⁴ / 1 ... no
# Actually (2+√3)/(2-√3) = (2+√3)²/((2-√3)(2+√3)) = (2+√3)²/1 = (2+√3)²
# So m_c/m_s = (2+√3)⁴? No, we defined m_c = (2+√3)² m_s, so m_c/m_s = (2+√3)² = 7+4√3
# That matches (2+√3)² not (2+√3)⁴.
# The seed interpretation: m_base such that m_s = (2-√3)m_base, m_c = (2+√3)m_base
# Then m_c/m_s = (2+√3)/(2-√3) = (2+√3)² (rationalizing). Yes! That works.
# So m_s = (2-√3)m_base → m_base = m_s/(2-√3) = m_s(2+√3)
m_base = m_s * (2 + math.sqrt(3))
print(f"\n    Seed base mass: m_base = m_s/(2-√3) = m_s(2+√3) = {m_base:.4f} MeV")
print(f"    Check: (2-√3)×m_base = {(2-math.sqrt(3))*m_base:.4f} = m_s = {m_s}")
print(f"    Check: (2+√3)×m_base = {(2+math.sqrt(3))*m_base:.4f} = m_c = {m_c_pred:.4f}")

# So the proper seed triple for the Koide check is (0, (2-√3)m, (2+√3)m) = (0, m_s, m_c)
# with m = m_base. And Q = 2/3 exactly.
# Let's verify numerically with the actual masses:
Q_seed_check = (0 + m_s + m_c_pred) / (0 + math.sqrt(m_s) + math.sqrt(m_c_pred))**2
print(f"    Q(0, m_s, m_c) = {Q_seed_check:.10f}")
print(f"    2/3 = {2/3:.10f}")
print(f"    Difference: {abs(Q_seed_check - 2/3):.2e}")

# Level 2: Bloom triple (-√m_s, √m_c, √m_b)
print("\n--- Level 2: Quark bloom triple (-√m_s, √m_c, √m_b) ---")
sqrt_bloom = [-math.sqrt(m_s), math.sqrt(m_c_pred), math.sqrt(m_b_pred)]
z0_bloom, z_bloom, m_bloom, Q_bloom = compute_charges("Bloom (-√m_s, √m_c, √m_b)", sqrt_bloom)

# Level 3: Heavy triple (√m_c, √m_b, √m_t)
print("\n--- Level 3: Heavy quark triple (√m_c, √m_b, √m_t) ---")
sqrt_heavy = [math.sqrt(m_c_pred), math.sqrt(m_b_pred), math.sqrt(m_t_pred)]
z0_heavy, z_heavy, m_heavy, Q_heavy = compute_charges("Heavy (√m_c, √m_b, √m_t)", sqrt_heavy)

# Level 4: Lepton triple (√m_e, √m_μ, √m_τ)
print("\n--- Level 4: Lepton triple (√m_e, √m_μ, √m_τ) ---")
sqrt_lepton = [math.sqrt(m_e), math.sqrt(m_mu), math.sqrt(m_tau_pred)]
z0_lepton, z_lepton, m_lepton, Q_lepton_ch = compute_charges("Lepton (√m_e, √m_μ, √m_τ)", sqrt_lepton)

# ============================================================
# SUMMARY
# ============================================================
print("\n" + "=" * 70)
print("SUMMARY TABLE DATA")
print("=" * 70)

print(f"\nLevel | Relation | Predicts | Value (MeV) | PDG (MeV) | Pull")
print("-" * 70)
print(f"  0   | —        | —        | —           | —         | — | 13 free")
print(f"  1   | R1 (O'R) | m_c      | {m_c_pred:.1f}      | {pdg['m_c'][0]}±{pdg['m_c'][1]}   | {p_c:+.2f}σ | 12 free")
print(f"  2   | R2 (bion)| m_b      | {m_b_pred:.1f}      | {pdg['m_b'][0]}±{pdg['m_b'][1]}   | {p_b:+.2f}σ | 11 free")
print(f"  3   | R3 (Yuk) | m_t      | {m_t_pred:.0f}   | {pdg['m_t'][0]}±{pdg['m_t'][1]}  | {p_t:+.2f}σ | 10 free")
print(f"  4   | R4 (lep) | m_τ      | {m_tau_pred:.2f}  | {pdg['m_tau'][0]}±{pdg['m_tau'][1]}| {p_tau:+.2f}σ |  9 free")
print(f"  5   | R5 (GST) | |V_us|   | {Vus_pred:.4f}    | {pdg['Vus'][0]}±{pdg['Vus'][1]}| {p_Vus:+.2f}σ |  8 free")

# ============================================================
# CHARGE TABLE DATA
# ============================================================
print("\n" + "=" * 70)
print("CHARGE TABLE DATA")
print("=" * 70)

def print_charge_row(label, sqrt_masses):
    frak_z = list(sqrt_masses)
    z0 = sum(frak_z) / 3
    z = [fz - z0 for fz in frak_z]
    mean_z2 = sum(zk**2 for zk in z) / 3
    ratio = mean_z2 / z0**2 if z0 != 0 else float('inf')
    m_vals = [fz**2 for fz in frak_z]
    sum_m = sum(m_vals)
    sum_fz = sum(frak_z)
    Q = sum_m / sum_fz**2 if sum_fz != 0 else float('inf')
    print(f"  {label}")
    print(f"    𝔷 = ({frak_z[0]:+.4f}, {frak_z[1]:+.4f}, {frak_z[2]:+.4f})")
    print(f"    z₀ = {z0:.4f}")
    print(f"    z  = ({z[0]:+.4f}, {z[1]:+.4f}, {z[2]:+.4f})")
    print(f"    ⟨z²⟩ = {mean_z2:.4f},  z₀² = {z0**2:.4f},  ⟨z²⟩/z₀² = {ratio:.6f}")
    print(f"    Q = {Q:.6f}")
    return z0, z, ratio, Q

print("\nL1 — Seed:")
print_charge_row("(0, √m_s, √m_c)", [0, math.sqrt(m_s), math.sqrt(m_c_pred)])

print("\nL2 — Bloom:")
print_charge_row("(-√m_s, √m_c, √m_b)", [-math.sqrt(m_s), math.sqrt(m_c_pred), math.sqrt(m_b_pred)])

print("\nL3 — Heavy:")
print_charge_row("(√m_c, √m_b, √m_t)", [math.sqrt(m_c_pred), math.sqrt(m_b_pred), math.sqrt(m_t_pred)])

print("\nL4 — Lepton:")
print_charge_row("(√m_e, √m_μ, √m_τ)", [math.sqrt(m_e), math.sqrt(m_mu), math.sqrt(m_tau_pred)])

# ============================================================
# Detailed numerical verification
# ============================================================
print("\n" + "=" * 70)
print("DETAILED VERIFICATION")
print("=" * 70)

print(f"\nR1: (2+√3)² = {(2+math.sqrt(3))**2:.10f}")
print(f"    7 + 4√3  = {7 + 4*math.sqrt(3):.10f}")
print(f"    Match: {abs((2+math.sqrt(3))**2 - (7+4*math.sqrt(3))) < 1e-12}")

print(f"\nR2: √m_b = 3√m_s + √m_c")
print(f"    3×{math.sqrt(m_s):.6f} + {math.sqrt(m_c_pred):.6f} = {3*math.sqrt(m_s) + math.sqrt(m_c_pred):.6f}")
print(f"    (√m_b)² = {(3*math.sqrt(m_s) + math.sqrt(m_c_pred))**2:.4f}")

print(f"\nR3: Koide Q(c,b,t):")
print(f"    Numerator: {m_c_pred} + {m_b_pred:.4f} + {m_t_pred:.4f} = {m_c_pred + m_b_pred + m_t_pred:.4f}")
print(f"    Denominator: ({math.sqrt(m_c_pred):.4f} + {math.sqrt(m_b_pred):.4f} + {math.sqrt(m_t_pred):.4f})²")
denom_q = (math.sqrt(m_c_pred) + math.sqrt(m_b_pred) + math.sqrt(m_t_pred))**2
print(f"             = {denom_q:.4f}")
print(f"    Q = {(m_c_pred + m_b_pred + m_t_pred)/denom_q:.15f}")

print(f"\nR4: Koide Q(e,μ,τ):")
num_l = m_e + m_mu + m_tau_pred
denom_l = (math.sqrt(m_e) + math.sqrt(m_mu) + math.sqrt(m_tau_pred))**2
print(f"    Numerator: {num_l:.6f}")
print(f"    Denominator: {denom_l:.6f}")
print(f"    Q = {num_l/denom_l:.15f}")

print(f"\nR5: √(m_d/m_s) = √({m_d}/{m_s}) = {math.sqrt(m_d/m_s):.6f}")

# Additional: check the seed Koide is exact when using the O'Raifeartaigh structure
print(f"\n--- Exact algebraic verification of seed Q=2/3 ---")
print(f"For eigenvalues (0, (2-√3)m, (2+√3)m):")
print(f"  Sum of masses = 4m")
print(f"  Sum of √masses = (√(2-√3) + √(2+√3))√m")
sq1 = math.sqrt(2 - math.sqrt(3))
sq2 = math.sqrt(2 + math.sqrt(3))
print(f"  √(2-√3) = {sq1:.10f}")
print(f"  √(2+√3) = {sq2:.10f}")
print(f"  Sum = {sq1+sq2:.10f}")
print(f"  Sum² = {(sq1+sq2)**2:.10f} (should be 6)")
print(f"  Q = 4m / (6m) = 2/3 ✓")

# Also check: does the bloom triple satisfy Koide?
Q_bloom_check = (m_s + m_c_pred + m_b_pred) / (-math.sqrt(m_s) + math.sqrt(m_c_pred) + math.sqrt(m_b_pred))**2
print(f"\nBloom Koide check (using m=𝔷²):")
print(f"  Σm_k = m_s + m_c + m_b = {m_s + m_c_pred + m_b_pred:.4f}")
print(f"  (Σ𝔷_k)² = (-√m_s + √m_c + √m_b)² = {(-math.sqrt(m_s) + math.sqrt(m_c_pred) + math.sqrt(m_b_pred))**2:.4f}")
print(f"  Q = {Q_bloom_check:.10f}")
print(f"  Note: m_k = 𝔷_k² (the SQUARES of the signed roots)")
print(f"  So masses are still (m_s, m_c, m_b) but the Koide sum uses signed roots")

# The bloom triple: using the SIGNED roots, what is Q?
# Q = ((-√m_s)² + (√m_c)² + (√m_b)²) / ((-√m_s) + √m_c + √m_b)²
# = (m_s + m_c + m_b) / (-√m_s + √m_c + √m_b)²
# This is NOT the same as the standard Koide for (s,c,b).
# The question asks us to verify ⟨z²⟩ = z₀² for each triple.

# v0 doubling check
print(f"\n--- v₀ doubling ---")
v0_seed = (0 + math.sqrt(m_s) + math.sqrt(m_c_pred)) / 3
v0_bloom = (-math.sqrt(m_s) + math.sqrt(m_c_pred) + math.sqrt(m_b_pred)) / 3
v0_heavy = (math.sqrt(m_c_pred) + math.sqrt(m_b_pred) + math.sqrt(m_t_pred)) / 3
print(f"  v₀(seed)  = {v0_seed:.4f} √MeV")
print(f"  v₀(bloom) = {v0_bloom:.4f} √MeV")
print(f"  v₀(heavy) = {v0_heavy:.4f} √MeV")
print(f"  v₀(bloom)/v₀(seed) = {v0_bloom/v0_seed:.6f}")
print(f"  v₀(heavy)/v₀(bloom) = {v0_heavy/v0_bloom:.6f}")

