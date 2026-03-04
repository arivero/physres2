#!/usr/bin/env python3
"""
Fritzsch texture: complete computation for all 4 parts.
Clean version with correct sign conventions.
"""
import numpy as np
from numpy import sqrt, pi, cos, sin
from scipy.optimize import brentq

print("="*70)
print("PART 1: Fritzsch Texture — Eigenvalues and GST")
print("="*70)

# Characteristic polynomial of
#   M = (0  A  0)
#       (A* 0  B)
#       (0  B* C)
# is: λ³ - Cλ² - (|A|²+|B|²)λ + |A|²C = 0
# (taking A,B real and positive for now)
#
# Vieta: λ₁+λ₂+λ₃ = C, λ₁λ₂+λ₁λ₃+λ₂λ₃ = -(A²+B²), λ₁λ₂λ₃ = -A²C
#
# In the hierarchy |A|≪|B|≪C:
#   λ₃ ≈ C + B²/C                    (heaviest, positive)
#   λ₂ ≈ -B²/C                       (negative — middle mass)
#   λ₁ ≈ A²C/B² = A²/(B²/C)         (lightest, positive)
#
# Physical masses mᵢ = |λᵢ|:
#   m₃ ≈ C, m₂ ≈ B²/C, m₁ ≈ A²C/B² = A²/m₂
# So m₁/m₂ ≈ A²C²/B⁴ = (A/B)²(C/B)² — a double hierarchy.
# Equivalently: m₁m₂ ≈ A², m₂m₃ ≈ B² (to leading order).

# NUMERICAL VERIFICATION with good hierarchy:
for A, B, C in [(0.005, 0.3, 5.0), (0.01, 0.5, 10.0), (0.001, 0.1, 3.0)]:
    M = np.array([[0,A,0],[A,0,B],[0,B,C]])
    ev = np.linalg.eigvalsh(M)
    masses = sorted(np.abs(ev))
    print(f"\n  A={A}, B={B}, C={C}")
    print(f"  masses: {masses}")
    print(f"  m₁≈A²C/B² = {A**2*C/B**2:.6f} vs {masses[0]:.6f} ({abs(A**2*C/B**2 - masses[0])/masses[0]*100:.1f}%)")
    print(f"  m₂≈B²/C   = {B**2/C:.6f} vs {masses[1]:.6f} ({abs(B**2/C - masses[1])/masses[1]*100:.1f}%)")
    print(f"  m₃≈C      = {C:.6f} vs {masses[2]:.6f} ({abs(C - masses[2])/masses[2]*100:.1f}%)")

# GST RELATION: The 12-mixing from diagonalizing a SINGLE Fritzsch matrix.
# 2nd-order perturbation theory in A gives:
#   sin θ₁₂ ≈ A/|λ₂| × <2|V|1'> = A/(B²/C)·(B/√(B²+C²))
# But more directly from the 2×2 effective matrix for (1,2) after
# integrating out state 3:
#   M_eff²(1,2) ≈ (m₁  √m₁m₂)
#                  (√m₁m₂  m₂)
# Diag by rotation with tan θ = √(m₁/m₂).
# This is the GST relation: sin θ_C ≈ √(m₁/m₂).

print(f"\nGST verification:")
A,B,C = 0.005, 0.3, 5.0
M = np.array([[0,A,0],[A,0,B],[0,B,C]])
ev, U = np.linalg.eigh(M)
idx = np.argsort(np.abs(ev))
m = np.abs(ev[idx])
U_sorted = U[:,idx]
# Mixing: |U₁₂| = component of flavor-1 in mass eigenstate 2
# In the mass-ordered basis, the "12 mixing" = |U[0,1]| or |U[1,0]|
mix_12 = abs(U_sorted[1,0])  # how much of flavor-2 is in mass state 1
mix_21 = abs(U_sorted[0,1])  # how much of flavor-1 is in mass state 2
print(f"  |U₁₂| (flavor 1→mass 2) = {mix_21:.6f}")
print(f"  |U₂₁| (flavor 2→mass 1) = {mix_12:.6f}")
print(f"  √(m₁/m₂) = {sqrt(m[0]/m[1]):.6f}")
# These agree to 1st order.

print("\n" + "="*70)
print("PART 2: Fritzsch CKM relation — |V_us|")
print("="*70)

m_u, m_d = 2.16, 4.67          # MeV, MS-bar at 2 GeV
m_s, m_c = 93.4, 1270.0
m_b, m_t = 4180.0, 172760.0

r_d = sqrt(m_d/m_s)
r_u = sqrt(m_u/m_c)
V_us_PDG, V_us_err = 0.2243, 0.0008

print(f"\n  √(m_d/m_s) = {r_d:.6f}")
print(f"  √(m_u/m_c) = {r_u:.6f}")
print(f"\nFritzsch relation: |V_us| = |√(m_d/m_s) − e^{{iφ}} √(m_u/m_c)|")
print(f"  |V_us|² = m_d/m_s + m_u/m_c − 2√(m_d m_u/(m_s m_c)) cos φ")

# Key phase values
for phi_deg, phi in [(0,0), (90,pi/2), (180,pi)]:
    v = sqrt(r_d**2 + r_u**2 - 2*r_d*r_u*cos(phi))
    dev = (v - V_us_PDG)/V_us_err
    print(f"  φ = {phi_deg:>3d}°:  |V_us| = {v:.6f}  ({dev:+.1f}σ)")

# Optimal phase
cos_phi_opt = (r_d**2 + r_u**2 - V_us_PDG**2) / (2*r_d*r_u)
phi_opt = np.arccos(cos_phi_opt)
print(f"\n  Optimal φ = {np.degrees(phi_opt):.2f}°  (cos φ = {cos_phi_opt:.4f})")
print(f"  At this φ: |V_us| = {V_us_PDG} (by construction)")
print(f"  PDG: |V_us| = {V_us_PDG} ± {V_us_err}")

# Comment on the φ=π/2 case (Fritzsch's classic prediction):
V_us_classic = sqrt(r_d**2 + r_u**2)
print(f"\n  Classic Fritzsch (φ=π/2): |V_us| = {V_us_classic:.6f}")
print(f"  = √(m_d/m_s + m_u/m_c) = √({m_d/m_s:.6f} + {m_u/m_c:.6f})")
print(f"  θ_C = {np.degrees(np.arcsin(V_us_classic)):.2f}° vs PDG {np.degrees(np.arcsin(V_us_PDG)):.2f}°")
print(f"  Deviation: {(V_us_classic - V_us_PDG)/V_us_err:+.1f}σ")

print("\n" + "="*70)
print("PART 3: Oakes relation")
print("="*70)

m_pi = 139.57  # MeV (charged pion)
m_K  = 493.68  # MeV (charged kaon)

# GMOR relations:
#   m_π² = (m_u + m_d) B₀
#   m_K⁺² = (m_u + m_s) B₀
#   m_K⁰² = (m_d + m_s) B₀
#
# Using the K⁰ relation (appropriate for the d/s ratio):
#   m_d/m_s from K⁰: m_π²/m_K² = (m_u+m_d)/(m_d+m_s)
# In the chiral limit m_u → 0:
#   m_π²/m_K⁰² → m_d/(m_d+m_s) = 1/(1 + m_s/m_d)
# So m_s/m_d = m_K²/m_π² − 1 = (m_K²−m_π²)/m_π²
# → m_d/m_s = m_π²/(m_K²−m_π²) [at m_u = 0]
#
# With the (average K) Weinberg 1977 version:
#   m_π² = 2m̂ B₀  (m̂ = (m_u+m_d)/2),  m_K² = (m̂+m_s)B₀
#   m̂/m_s = m_π²/(2m_K²−m_π²)
#
# The problem states the Oakes form:
#   sin θ_C ≈ √(m_d/m_s) ≈ √((m_K²−m_π²)/(m_K²+m_π²))
#
# Let me check: GMOR with K⁺ (uses m_u+m_s):
#   m_π²=B₀(m_u+m_d), m_K²=B₀(m_u+m_s)
#   (m_K²−m_π²)/(m_K²+m_π²) = (m_s−m_d)/(m_s+m_d+2m_u)
# If m_u ≈ m_d/2: (m_s−m_d)/(m_s+2m_d) ≈ 1 − 3m_d/m_s + ...
# That gives ~0.85, not ~0.05.
#
# The formula √((m_K²−m_π²)/(m_K²+m_π²)) ≈ 0.923 is clearly cos θ_C, not sin.
# The correct Oakes relation for sin θ_C uses a DIFFERENT meson mass combination.
#
# The Weinberg 1977 formula gives the best GMOR estimate:
#   m_d/m_s ≈ m_π²/(2m_K² − m_π²)
#
# sin θ_C(Oakes/Weinberg) = √(m_π²/(2m_K²−m_π²)):

sC_GMOR = sqrt(m_pi**2 / (2*m_K**2 - m_pi**2))
print(f"\n  m_π = {m_pi} MeV, m_K = {m_K} MeV")
print(f"\n  Weinberg-GMOR: m̂/m_s = m_π²/(2m_K²−m_π²) = {m_pi**2/(2*m_K**2-m_pi**2):.6f}")
print(f"  sin θ_C(GMOR) = {sC_GMOR:.6f}")

# Alternative: m_u = 0 limit with K⁰
sC_Oakes_K0 = sqrt(m_pi**2 / (m_K**2 - m_pi**2))
print(f"\n  At m_u=0: m_d/m_s = m_π²/(m_K²−m_π²) = {m_pi**2/(m_K**2-m_pi**2):.6f}")
print(f"  sin θ_C(m_u=0) = {sC_Oakes_K0:.6f}")

# NOW: the formula as stated in the problem: √((m_K²−m_π²)/(m_K²+m_π²))
# This is cos θ_C in the Oakes relation, or perhaps the problem
# intends a different angle convention. Let me check:
# If tan θ_C = √(m_d/m_s), then:
#   sin θ_C = √(m_d/m_s)/√(1+m_d/m_s) = √(m_d/(m_d+m_s))
#   cos θ_C = 1/√(1+m_d/m_s) = √(m_s/(m_d+m_s))
# And from GMOR at m_u=0:
#   sin²θ_C = m_d/(m_d+m_s) = m_π²/m_K² (from m_K⁰²=(m_d+m_s)B₀)
#   cos²θ_C = m_s/(m_d+m_s) = (m_K²−m_π²)/m_K²
# That gives sin θ_C = m_π/m_K = 0.283 — too big.
#
# Or with the stated formula, if it's the Oakes (1969) TANGENT relation:
#   tan θ_C = √(m_d/m_s)
# Then sin θ_C = tan θ_C / √(1+tan²θ_C) = √(m_d/m_s)/√(1+m_d/m_s)
# At m_u=0: = √(m_π²/(m_K²-m_π²)) / √(1 + m_π²/(m_K²-m_π²))
#          = √(m_π²/(m_K²-m_π²)) / √(m_K²/(m_K²-m_π²))
#          = m_π/m_K = 0.283

# None of these give 0.224 from meson masses alone. The point is:
# the GMOR relation gives an APPROXIMATE value that depends on which 
# approximation is used for the quark mass ratio.

# The comparison requested:
print(f"\n--- Comparison table ---")
print(f"  {'Quantity':<50} {'Value':<10}")
print(f"  {'-'*60}")
print(f"  {'√(m_d/m_s) [GST, quark masses]':<50} {r_d:.6f}")
print(f"  {'m_π/√(2m_K²−m_π²) [Weinberg-GMOR]':<50} {sC_GMOR:.6f}")
print(f"  {'m_π/m_K [m_u=0 leading order]':<50} {m_pi/m_K:.6f}")
print(f"  {'|V_us|_PDG':<50} {V_us_PDG:.6f}")

# The Weinberg formula gives 0.204 which is 9% low — the m_u ≈ m_d isospin
# approximation is too crude. The exact GMOR result uses m_u ≠ m_d and
# the full set of meson masses.

# Using measured m_u/m_d ratio to improve:
# From GMOR: m_π² = (m_u+m_d)B₀, m_K⁺² = (m_u+m_s)B₀
# m_d/m_s = (m_K² - m_π² + m_π²·m_d/(m_u+m_d)) / (m_K²·(1 - m_u/(m_u+m_s)))
# This gets complicated. The point of the exercise is the comparison,
# not matching to arbitrary precision.

print("\n" + "="*70)
print("PART 4: Charge-language expression of GST")
print("="*70)

# Koide parametrization: √mₖ = z₀(1 + √2 cos(δ + 2πk/3))
# Constraints: Σzₖ = 0 (by construction), <zₖ²> = z₀² iff Q = 2/3 (identity).
#
# For a down-type Koide triple k=0→d, k=1→s, k=2→b:
#   sin²θ_C = m_d/m_s = [(1+√2 cos δ)/(1+√2 cos(δ+2π/3))]²
#
# At the SEED δ = 3π/4:
#   1 + √2 cos(3π/4) = 1 − 1 = 0 → m_d = 0, sin θ_C = 0
#   1 + √2 cos(17π/12) = (3−√3)/2    [proven analytically]
#   1 + √2 cos(π/12) = (3+√3)/2
#
# Small bloom δ = 3π/4 + ε:
#   Numerator ≈ −ε    [from Taylor expansion of cos(3π/4+ε)]
#   Denominator ≈ (3−√3)/2   [to zeroth order in ε]
#   sin θ_C ≈ 2|ε|/(3−√3) = |ε|√2/(√3−1) = |ε|(√3+1)/√2·... 
#   = 2|ε|/(3−√3) = 2|ε|(3+√3)/((3−√3)(3+√3)) = 2|ε|(3+√3)/6

# Simplify: 2/(3−√3) = 2(3+√3)/((9−3)) = (3+√3)/3

approx_coeff = 2/(3-sqrt(3))
approx_coeff_v2 = (3+sqrt(3))/3
print(f"\n  Small-bloom coefficient: 2/(3−√3) = {approx_coeff:.6f}")
print(f"  = (3+√3)/3 = {approx_coeff_v2:.6f}")

# Find ε for physical sin θ_C = √(m_d/m_s)
target = sqrt(m_d/m_s)

def sin_theta(delta):
    f0 = (1 + sqrt(2)*cos(delta))**2
    f1 = (1 + sqrt(2)*cos(delta + 2*pi/3))**2
    return sqrt(f0/f1) if f1 > 0 else float('inf')

# Root finding
delta_sol = brentq(lambda d: sin_theta(d) - target, 3*pi/4 + 0.01, 3*pi/4 + 0.5)
eps = delta_sol - 3*pi/4

print(f"\n  Target sin θ_C = {target:.6f}")
print(f"  Solution: ε = {eps:.6f} rad = {np.degrees(eps):.2f}°")
print(f"  δ = 3π/4 + ε = {delta_sol:.6f} rad = {np.degrees(delta_sol):.2f}°")
print(f"  sin θ_C(exact) = {sin_theta(delta_sol):.6f}")
print(f"  sin θ_C(approx) = (3+√3)ε/3 = {(3+sqrt(3))*eps/3:.6f}")
print(f"  Approx error: {abs((3+sqrt(3))*eps/3 - target)/target*100:.1f}%")

# The key physical content:
# (1) At the seed, m_d = 0 and θ_C = 0.
# (2) The Cabibbo angle is a MONOTONIC function of ε, the bloom distance.
# (3) sin θ_C = |ℨ_d/ℨ_s| = ratio of mass charges.

# What z₀ would be needed if (d,s,b) satisfied Koide?
f0 = (1 + sqrt(2)*cos(delta_sol))**2
f1 = (1 + sqrt(2)*cos(delta_sol + 2*pi/3))**2
f2 = (1 + sqrt(2)*cos(delta_sol + 4*pi/3))**2
z0_sq = m_s / f1
z0 = sqrt(z0_sq)
m_b_koide = z0_sq * f2
print(f"\n  If Koide held for (d,s,b):")
print(f"    z₀ = {z0:.4f} MeV^(1/2)")
print(f"    m_d(input) = {m_d} MeV ✓")
print(f"    m_s(input) = {m_s} MeV ✓")
print(f"    m_b(Koide) = {m_b_koide:.1f} MeV (PDG: {m_b})")
print(f"  Koide does NOT hold for (d,s,b): Q = {(m_d+m_s+m_b)/(sqrt(m_d)+sqrt(m_s)+sqrt(m_b))**2:.4f} vs 2/3")
print(f"  But the GST relation sin θ_C = √(m_d/m_s) = |ℨ_d/ℨ_s| is independent of m_b.")

# Constraint summary:
# Given Koide (Q=2/3) and the angle δ:
#   - z₀ sets the overall mass scale (e.g., z₀² = m_s/f₁)  
#   - δ encodes all mass ratios within the triple
#   - sin θ_C = √(m_d/m_s) is a function of δ alone
#   - The bloom ε = δ − 3π/4 is the single parameter controlling θ_C

# The mass-charge identity <zₖ²> = z₀² is an identity for ALL δ
# (it's equivalent to Q = 2/3 being the angular-average identity
#  Σcos²(δ+2πk/3) = 3/2 for all δ).
# So the constraint on (z₀, δ) from sin θ_C is simply: δ is fixed.
# z₀ remains free (it's the overall mass scale).

print(f"\n  Summary of constraints:")
print(f"  <zₖ²> = z₀² is AUTOMATIC (identity ∀δ, equivalent to Q=2/3)")
print(f"  sin θ_C = √(m_d/m_s) fixes δ = 3π/4 + ε with ε = {eps:.4f}")
print(f"  z₀ remains free (sets the mass scale)")
print(f"  sin θ_C ≈ (3+√3)ε/3 for small ε (bloom from seed)")

print("\n" + "="*70)
print("COMPLETE NUMERICAL SUMMARY")
print("="*70)

print(f"\n  Part 1 — Eigenvalues:")
print(f"    m₁ ≈ |A|²C/|B|²    (lightest)")
print(f"    m₂ ≈ |B|²/C        (middle)")
print(f"    m₃ ≈ C             (heaviest)")
print(f"    sin θ_C ≈ √(m₁/m₂)  [GST relation]")

print(f"\n  Part 2 — Fritzsch relation:")
print(f"    |V_us| = |√(m_d/m_s) − e^{{iφ}} √(m_u/m_c)|")
print(f"    √(m_d/m_s) = {r_d:.6f}")
print(f"    √(m_u/m_c) = {r_u:.6f}")
print(f"    φ = 0:   |V_us| = {abs(r_d-r_u):.6f}  ({(abs(r_d-r_u)-V_us_PDG)/V_us_err:+.1f}σ)")
print(f"    φ = π/2: |V_us| = {sqrt(r_d**2+r_u**2):.6f}  ({(sqrt(r_d**2+r_u**2)-V_us_PDG)/V_us_err:+.1f}σ)")
print(f"    φ_opt = {np.degrees(phi_opt):.1f}°: |V_us| = {V_us_PDG}")
print(f"    PDG:     |V_us| = {V_us_PDG} ± {V_us_err}")

print(f"\n  Part 3 — Oakes relation:")
print(f"    √(m_d/m_s)                = {r_d:.6f}  (−0.9σ from PDG)")
print(f"    m_π/√(2m_K²−m_π²) [GMOR]  = {sC_GMOR:.6f}  (leading-order chiral)")
print(f"    |V_us|_PDG                 = {V_us_PDG}")
theta_GST = np.degrees(np.arcsin(r_d))
theta_GMOR = np.degrees(np.arcsin(sC_GMOR))
theta_PDG = np.degrees(np.arcsin(V_us_PDG))
print(f"    θ_C(GST)  = {theta_GST:.2f}°")
print(f"    θ_C(GMOR) = {theta_GMOR:.2f}°")
print(f"    θ_C(PDG)  = {theta_PDG:.2f}°")

print(f"\n  Part 4 — Charge language:")
print(f"    sin θ_C = |ℨ_d/ℨ_s|")
print(f"    = |(1+√2 cos δ)/(1+√2 cos(δ+2π/3))|")
print(f"    Seed: δ₀ = 3π/4 → m_d = 0, θ_C = 0")
print(f"    Physical: ε = δ−δ₀ = {eps:.4f} rad = {np.degrees(eps):.1f}°")
print(f"    Small-bloom: sin θ_C ≈ (3+√3)ε/3")

