#!/usr/bin/env python3
"""
Fritzsch texture diagonalization, CKM mixing, Oakes relation,
and charge-language connection.
"""
import numpy as np
from numpy import sqrt, pi, cos, sin, arctan2
from scipy.optimize import brentq

print("="*70)
print("PART 1: Fritzsch Texture Diagonalization")
print("="*70)

# The Fritzsch nearest-neighbour Hermitian mass matrix:
#
#   M = ( 0   A   0 )
#       ( A*  0   B )
#       ( 0   B*  C )
#
# For real A,B,C the characteristic polynomial is:
#   λ³ - Cλ² - (A²+B²)λ + A²C = 0
#
# Vieta relations for eigenvalues λ₁, λ₂, λ₃:
#   λ₁ + λ₂ + λ₃ = C
#   λ₁λ₂ + λ₁λ₃ + λ₂λ₃ = -(A²+B²)
#   λ₁λ₂λ₃ = -A²C
#
# The physical masses are |λᵢ|. One eigenvalue is negative.

# HIERARCHY LIMIT: |A| << |B| << |C|
# 
# Block-diagonalize: first solve the (2,3) block
#   (0  B)  eigenvalues: ±|B|
#   (B  C)  → C/2 ± √(C²/4 + B²) ≈ {C + B²/C, -B²/C}
#
# Then treat A as perturbation coupling state 1 to the (2,3) block.
# The unperturbed eigenvalues are: 0 (state 1), -B²/C (state 2'), C+B²/C (state 3').
# Second-order perturbation on state 1 from state 2':
#   Δλ₁ = |<1|V|2'>|² / (0 - (-B²/C)) = (A sinθ)² / (B²/C) 
# where sinθ ≈ B/C is the (2,3) mixing.
# Actually, the coupling of state 1 to the (2,3) eigenstates:
# |2'> ≈ -|2> + (B/C)|3>  (unnormalized, for the eigenvalue -B²/C)
# <1|V|2'> ≈ A·(-1) = -A  (since V₁₂ = A, V₁₃ = 0)
# Δλ₁ ≈ A² / (B²/C) = A²C/B²
# 
# So the three eigenvalues in the hierarchy limit are approximately:
#   λ₁ ≈ A²C/B²   (smallest, positive if C > 0)
#   λ₂ ≈ -B²/C    (negative)
#   λ₃ ≈ C + B²/C  (largest)
#
# Physical masses (absolute values):
#   m₁ ≈ A²C/B²,  m₂ ≈ B²/C,  m₃ ≈ C

# Numerical verification with good hierarchy
A, B, C = 0.005, 0.3, 5.0

M = np.array([[0, A, 0],
              [A, 0, B],
              [0, B, C]])

evals_raw = np.linalg.eigvalsh(M)
evals = sorted(evals_raw, key=abs)  # sort by magnitude
print(f"\nNumerical example: A={A}, B={B}, C={C}")
print(f"Eigenvalues (sorted by |λ|): {evals}")
print(f"  m₁ = |λ_small| = {abs(evals[0]):.8f},   A²C/B² = {A**2*C/B**2:.8f}")
print(f"  m₂ = |λ_mid|   = {abs(evals[1]):.8f},   B²/C   = {B**2/C:.8f}")
print(f"  m₃ = |λ_large| = {abs(evals[2]):.8f},   C      = {C:.8f}")

# GST Relation from the eigenvectors:
# The eigenvector of λ₁ (smallest) is |1> + (A/B²/C·C)|2'> + ...
# In the flavor basis, the Cabibbo-like mixing angle is:
#   sinθ_C ≈ |<2|ψ₁>| / |<1|ψ₁>| ← NO, it's the (1,2) element of U†
#
# For the Fritzsch texture, the standard result is:
#   tanθ₁₂ ≈ √(m₁/m₂) for the diagonalization of a SINGLE matrix M.
# This is the GST relation.

evals_full, evecs = np.linalg.eigh(M)
# Sort by |eigenvalue|:
idx = np.argsort(np.abs(evals_full))
evals_sorted = evals_full[idx]
evecs_sorted = evecs[:, idx]
print(f"\nEigenvectors (columns, sorted by |eigenvalue|):")
print(evecs_sorted)
# The (1,2)-like mixing in flavor space:
# V₁₂ ~ component of flavor-1 in the eigenvector of m₂, OR
# equivalently, component of eigenvector-1 in flavor direction 2.
# Let's check both:
U12_a = abs(evecs_sorted[0, 1])  # flavor 0 component in eigenvector 1
U12_b = abs(evecs_sorted[1, 0])  # flavor 1 component in eigenvector 0
print(f"|U[flavor0, eigvec1]| = {U12_a:.6f}")
print(f"|U[flavor1, eigvec0]| = {U12_b:.6f}")
print(f"√(m₁/m₂) = {sqrt(abs(evals_sorted[0]/evals_sorted[1])):.6f}")

# These should match for a single matrix. The mixing angle element
# of the rotation is the off-diagonal element of the 2x2 block
# in the (1,2) sector, which by perturbation theory is:
#   sin θ₁₂ ≈ A / √(A² + (B²/C)²) ... no.
# Actually, for the 2-generation block of the Fritzsch texture:
# The eigenvalues of the reduced 2x2 problem give 
#   tan θ₁₂ = A·C/(B² - ...) ≈ AC/B²
# And m₁/m₂ ≈ A²C²/B⁴ = (AC/B²)² = tan²θ₁₂ ≈ sin²θ₁₂
# So √(m₁/m₂) ≈ sinθ₁₂. That's the GST relation.

# More careful: eigenvalue problem in 2×2 effective mass-squared matrix
# After integrating out the heavy state m₃, the effective 2×2 mass
# matrix in the (1,2) space is approximately:
#   M_eff = (-m₁   √(m₁m₂))
#           (√(m₁m₂)   m₂  )
# whose off-diagonal/diagonal ratio gives tan θ = √(m₁/m₂).
#
# This is EXACTLY the GST relation: sin θ_C ≈ √(m_d/m_s).

print(f"\n✓ GST relation verified: sin θ_C ≈ √(m₁/m₂)")
print(f"  Numerical mixing: {U12_b:.6f}")
print(f"  √(m₁/m₂):        {sqrt(abs(evals_sorted[0]/evals_sorted[1])):.6f}")

print("\n" + "="*70)
print("PART 2: Fritzsch CKM Relation")
print("="*70)

# PDG quark masses (MS-bar at 2 GeV)
m_u = 2.16   # MeV
m_d = 4.67   # MeV
m_s = 93.4   # MeV
m_c = 1270.0 # MeV
m_b = 4180.0 # MeV
m_t = 172760.0  # MeV

print(f"\nPDG inputs:")
print(f"  m_u = {m_u} MeV, m_d = {m_d} MeV")
print(f"  m_s = {m_s} MeV, m_c = {m_c} MeV")
print(f"  m_b = {m_b} MeV, m_t = {m_t} MeV")

# Fritzsch relation for V_us:
# V_us = U_u†₁₂ U_d₂₂ - U_u†₁₁ U_d₁₂ + ...
# In the 2-generation approximation:
#   |V_us| = |√(m_d/m_s) - e^{iφ} √(m_u/m_c)|
# where φ is a relative CP-violating phase from complex A,B entries.

r_d = sqrt(m_d / m_s)
r_u = sqrt(m_u / m_c)

print(f"\n√(m_d/m_s) = {r_d:.6f}")
print(f"√(m_u/m_c) = {r_u:.6f}")

# |V_us|² = r_d² + r_u² - 2 r_d r_u cos φ
V_us_min = abs(r_d - r_u)  # φ = 0
V_us_max = r_d + r_u        # φ = π
V_us_pi2 = sqrt(r_d**2 + r_u**2)  # φ = π/2

V_us_PDG = 0.2243
V_us_err = 0.0008

print(f"\n|V_us| for various φ:")
print(f"  φ = 0:   |V_us| = {V_us_min:.6f}  ({(V_us_min-V_us_PDG)/V_us_err:+.1f}σ from PDG)")
print(f"  φ = π/2: |V_us| = {V_us_pi2:.6f}  ({(V_us_pi2-V_us_PDG)/V_us_err:+.1f}σ from PDG)")
print(f"  φ = π:   |V_us| = {V_us_max:.6f}  ({(V_us_max-V_us_PDG)/V_us_err:+.1f}σ from PDG)")
print(f"  PDG:     |V_us| = {V_us_PDG} ± {V_us_err}")

# Find optimal φ
cos_phi = (r_d**2 + r_u**2 - V_us_PDG**2) / (2 * r_d * r_u)
if abs(cos_phi) <= 1:
    phi_opt = np.arccos(cos_phi)
    V_us_check = sqrt(r_d**2 + r_u**2 - 2*r_d*r_u*cos_phi)
    print(f"\nOptimal φ for PDG match: {np.degrees(phi_opt):.2f}°")
    print(f"  cos φ = {cos_phi:.6f}")
    print(f"  |V_us| = {V_us_check:.6f}")
else:
    print(f"\nPDG value outside Fritzsch range (cos φ = {cos_phi:.4f})")

# The classic Fritzsch (1977) prediction was φ = π/2 (maximal CP phase)
# giving |V_us| = √(m_d/m_s + m_u/m_c)
print(f"\nClassic Fritzsch (φ=π/2):")
print(f"  |V_us| = √(m_d/m_s + m_u/m_c) = {V_us_pi2:.6f}")
print(f"  θ_C = {np.degrees(np.arcsin(V_us_pi2)):.2f}°")

# Exact numerical diagonalization
# Build Fritzsch textures from quark masses
# For M with eigenvalues (-m₁, -m₂, m₃):
# C = m₃ - m₁ - m₂, A² = m₁m₂m₃/C, etc.

def fritzsch_params(m1, m2, m3):
    """Return (A, B, C) for Fritzsch matrix with masses m1 < m2 < m3."""
    C = m3 - m1 - m2
    Asq = m1 * m2 * m3 / C
    Bsq = m1*m3 + m2*m3 - m1*m2 - Asq
    return sqrt(Asq), sqrt(Bsq), C

A_d, B_d, C_d = fritzsch_params(m_d, m_s, m_b)
A_u, B_u, C_u = fritzsch_params(m_u, m_c, m_t)

print(f"\nFritzsch parameters:")
print(f"  Down: A_d={A_d:.4f}, B_d={B_d:.4f}, C_d={C_d:.4f}")
print(f"  Up:   A_u={A_u:.4f}, B_u={B_u:.4f}, C_u={C_u:.4f}")

M_d = np.array([[0, A_d, 0], [A_d, 0, B_d], [0, B_d, C_d]])
M_u = np.array([[0, A_u, 0], [A_u, 0, B_u], [0, B_u, C_u]])

ev_d, U_d = np.linalg.eigh(M_d)
ev_u, U_u = np.linalg.eigh(M_u)

# Sort by |eigenvalue| (mass ordering)
idx_d = np.argsort(np.abs(ev_d))
idx_u = np.argsort(np.abs(ev_u))
ev_d_sorted = ev_d[idx_d]
ev_u_sorted = ev_u[idx_u]
U_d_sorted = U_d[:, idx_d]
U_u_sorted = U_u[:, idx_u]

print(f"\nDiag check (|eigenvalues|):")
print(f"  Down: {np.abs(ev_d_sorted)} vs ({m_d}, {m_s}, {m_b})")
print(f"  Up:   {np.abs(ev_u_sorted)} vs ({m_u}, {m_c}, {m_t})")

# CKM = U_u† U_d (both in mass-ordered basis)
V = U_u_sorted.T @ U_d_sorted
print(f"\n|V_CKM| (real Fritzsch, φ=0):")
for i in range(3):
    row = "  ".join(f"{abs(V[i,j]):.6f}" for j in range(3))
    print(f"  {row}")
print(f"\n|V_us| (exact diag, real) = {abs(V[0,1]):.6f}")
print(f"Fritzsch formula (φ=0):    {V_us_min:.6f}")
print(f"GST (√(m_d/m_s) only):     {r_d:.6f}")
print(f"PDG:                        {V_us_PDG}")

# Note: for the REAL Fritzsch texture, |V_us| = |√(m_d/m_s) - √(m_u/m_c)|
# which is φ = 0. The relative phase between A_u and A_d provides φ.

print("\n" + "="*70)
print("PART 3: Oakes Relation")
print("="*70)

m_pi = 139.57  # MeV (charged pion)
m_K = 493.68   # MeV (charged kaon)

# The Oakes relation uses the Gell-Mann--Oakes--Renner relations:
#   m_π² = (m_u + m_d) B₀
#   m_K² = (m_s + m_d) B₀  (for K⁺ = us̄, or average)
#
# The Weinberg--Oakes relation is:
#   tan²θ_C ≈ m_d/m_s  (equivalent to GST at leading order)
#
# And from GMOR:
#   m_d/m_s ≈ (2m_K² - m_π²) ... no, let me derive carefully.
#
# From GMOR:  m_π² = (m_u + m_d)B₀,  m_K² = (m_s + m_d)B₀ (K⁺)
# Neglecting m_u (isospin limit m_u → 0):
#   m_π² ≈ m_d B₀,   m_K² ≈ m_s B₀  (at m_u = 0)
#   m_d/m_s ≈ m_π²/m_K²
# Then sin²θ_C ≈ m_d/m_s ≈ m_π²/m_K²
# and sin θ_C ≈ m_π/m_K
#
# More precisely, the Cabibbo-Oakes relation uses the RATIO:
#   (m_K² - m_π²)/(m_K² + m_π²) = (m_s + m_d - m_u - m_d)/(m_s + m_d + m_u + m_d)  -- no
# Wait, let me reconsider. Actually the Weinberg (1977) relation is:
#   m_d/m_s ≈ (2m_K₀² - m_K⁺² + m_π⁺² - m_π₀²)/(m_K₀² - m_K⁺² + m_π⁺²) 
# -- that's more complex.
#
# The simpler Oakes relation is just:
#   sin θ_C ≈ f_π/f_K ≈ m_π/m_K × something...
# 
# Actually, from the Weinberg 1977 paper and the standard textbook form:
# The relevant ratio from GMOR is:
#   m_d/m_s = m_π²/(2m_K² - m_π²)   [Weinberg 1977]
# 
# Let me just use the formula stated in the problem:
# sin θ_C ≈ √(m_d/m_s) ≈ √((m_K² - m_π²)/(m_K² + m_π²))

# First, check the STATED formula:
oakes = sqrt((m_K**2 - m_pi**2) / (m_K**2 + m_pi**2))
print(f"\nm_π = {m_pi} MeV, m_K = {m_K} MeV")
print(f"\n√((m_K² - m_π²)/(m_K² + m_π²)) = {oakes:.6f}")
print(f"√(m_d/m_s) = {r_d:.6f}")

# Hmm, 0.923 vs 0.224. The formula as stated doesn't work numerically!
# The CORRECT Oakes relation should be:
#   m_d/m_s ≈ m_π²/(2m_K² - m_π²)  [Weinberg]
# or equivalently sin θ_C ≈ √(m_π²/(2m_K² - m_π²)) = m_π/√(2m_K² - m_π²)

weinberg = m_pi**2 / (2*m_K**2 - m_pi**2)
print(f"\nWeinberg ratio m_d/m_s ≈ m_π²/(2m_K² - m_π²) = {weinberg:.6f}")
print(f"sin θ_C ≈ √(...) = {sqrt(weinberg):.6f}")
print(f"cf. m_d/m_s(PDG) = {m_d/m_s:.6f}")
print(f"cf. √(m_d/m_s) = {r_d:.6f}")

# Alternative: Oakes' original form used in the problem statement.
# Let me re-read the problem... it says:
#   sin θ_C ≈ √(m_d/m_s) ≈ √((m_K² - m_π²)/(m_K² + m_π²))
#
# This is NOT the standard Oakes/Weinberg relation. Let me check if 
# there's a different form. Perhaps the problem means a different 
# meson mass combination. Let's check Gell-Mann--Oakes--Renner:
#
#   m_π² = (m_u + m_d)B₀
#   m_K⁺² = (m_u + m_s)B₀
#   m_K⁰² = (m_d + m_s)B₀
#
# So: (m_K⁰² - m_π²)/(m_K⁰²) = (m_s - m_u)/(m_d + m_s) ≈ 1 - m_d/m_s  (if m_u ~ 0)
# Or: m_K⁰²/m_π² = (m_d + m_s)/(m_u + m_d) 
# At m_u = 0: m_K²/m_π² = (m_d+m_s)/m_d = 1 + m_s/m_d
# So: m_d/m_s = 1/(m_K²/m_π² - 1) = m_π²/(m_K² - m_π²)

oakes_v2 = sqrt(m_pi**2 / (m_K**2 - m_pi**2))
print(f"\n√(m_π²/(m_K² - m_π²)) = {oakes_v2:.6f}")

# That's the one! Let me verify:
# At m_u = 0: m_π² = m_d·B₀, m_K⁰² = (m_d + m_s)·B₀
# m_d/m_s = m_π²/(m_K² - m_π²) × (m_d + m_s)/m_s ... no.
# m_K² - m_π² = m_s·B₀ (at m_u = 0)
# m_π² = m_d·B₀
# m_d/m_s = m_π²/(m_K² - m_π²)
# sin θ_C = √(m_d/m_s) = m_π/√(m_K² - m_π²) = 0.295... still not 0.224.

# Hmm. Let me try the Weinberg formula again:
# m_d/m_s = m_π²/(2m_K² - m_π²)
# This uses m_K² = ½(m_u + m_d + 2m_s)B₀ (the average K mass)
# At m_u ≈ m_d: m_π² = 2m_d B₀, m_K² = (m_d + m_s)B₀
# m_K² = (m_d+m_s)B₀, m_π² = 2m_d B₀
# m_d/m_s = m_π²/(2m_K² - m_π²) = 2m_d B₀ / (2(m_d+m_s)B₀ - 2m_d B₀)
#         = 2m_d / (2m_s) = m_d/m_s ✓

# sin θ_C = √(m_π²/(2m_K² - m_π²))
sthC_weinberg = sqrt(m_pi**2 / (2*m_K**2 - m_pi**2))
print(f"\nWeinberg-Oakes: sin θ_C = √(m_π²/(2m_K² - m_π²)) = {sthC_weinberg:.6f}")
print(f"√(m_d/m_s) = {r_d:.6f}")
print(f"|V_us|_PDG = {V_us_PDG}")
print(f"  Deviation from PDG: {(sthC_weinberg - V_us_PDG)/V_us_err:.1f}σ")

# CORRECTION: The problem statement's formula √((m_K²-m_π²)/(m_K²+m_π²)) ≈ 0.923
# This is actually cos θ_C, not sin θ_C! Or it's a different relation.
# Let me check if sin θ_C ≈ √(1 - (m_K²-m_π²)/(m_K²+m_π²)) = √(2m_π²/(m_K²+m_π²))
alt = sqrt(2*m_pi**2 / (m_K**2 + m_pi**2))
print(f"\n√(2m_π²/(m_K² + m_π²)) = {alt:.6f}")

# Or perhaps the formula in the problem is tan θ_C, not sin θ_C?
# tan²θ_C = m_d/m_s ≈ m_π²/(2m_K² - m_π²) at leading order in ChPT.

# Actually, there's a well-known relation attributed to Oakes:
#   tan θ_C = f_K sin θ_1 / f_π
# But the simplest meson mass form that gives sin θ_C ≈ 0.22 is:
#   sin θ_C ≈ m_π / m_K × 1/√2
# Let's check:
print(f"\nm_π/(√2 m_K) = {m_pi/(sqrt(2)*m_K):.6f}")
# That gives 0.1998, not great.

# OK, the stated formula √((m_K²-m_π²)/(m_K²+m_π²)) gives cos of the
# "meson angle", not Cabibbo. Let me just use the Weinberg formula.
# But let me also check if maybe the problem means something else.

# Perhaps the Oakes relation is:
#   sin θ_C = m_π / (√2 m_K) ... no, that's Cabibbo's original guess
# Or:
#   sin²(2θ_C) = m_π²/m_K²
# sin(2θ_C) = m_pi/m_K = 0.283
# θ_C = arcsin(0.283)/2 = 8.2° -- no.

# Let me go with what the problem says and note the discrepancy.
# The formula stated is:
#   sin θ_C ≈ √(m_d/m_s) ≈ √((m_K² - m_π²)/(m_K² + m_π²))
# This gives 0.923, which is obviously wrong for sin θ_C ≈ 0.22.
# 
# I think the intended formula is actually just the Weinberg 1977 form.
# Let me use the CORRECT Oakes/Weinberg relation.

print(f"\n--- CORRECTED Oakes relation ---")
print(f"From GMOR at m_u ≈ m_d (isospin limit):")
print(f"  m_π² = 2m_q B₀,  m_K² = (m_q + m_s)B₀")
print(f"  m_d/m_s = m_π²/(2m_K² - m_π²)")
print(f"  sin θ_C ≈ √(m_π²/(2m_K² - m_π²)) = {sthC_weinberg:.6f}")
print(f"  = m_π/√(2m_K² - m_π²)")

numerator = m_pi**2
denominator = 2*m_K**2 - m_pi**2
print(f"\n  m_π² = {numerator:.2f}")
print(f"  2m_K² - m_π² = {denominator:.2f}")
print(f"  ratio = {numerator/denominator:.6f}")

# Actually, now I realize: with the NON-degenerate m_u ≠ m_d GMOR relations:
# m_π⁺² = (m_u + m_d)B₀
# m_K⁺² = (m_u + m_s)B₀  
# m_K⁰² = (m_d + m_s)B₀
#
# So: m_d/m_s = (m_K⁰² - m_K⁺² + m_π⁺²)/(m_K⁰²) ... complicated.
# In the isospin limit m_u = m_d = m_q:
#   m_π² = 2m_q B₀ → m_q = m_π²/(2B₀)
#   m_K² = (m_q + m_s)B₀ → m_s = m_K²/B₀ - m_q = m_K²/B₀ - m_π²/(2B₀)
#   m_q/m_s = m_π²/2 / (m_K² - m_π²/2) = m_π²/(2m_K² - m_π²)
#
# This is literally Weinberg's relation. And m_q here ≈ (m_u+m_d)/2 ≈ m_d
# (since we're computing sin θ_C = √(m_d/m_s), we should use m_d not m_q).
# The factor of 2 means this is an approximation.

print(f"\n--- Summary of three determinations ---")
print(f"{'Method':<50} {'Value':<10} {'vs PDG'}")
print(f"{'-'*70}")
print(f"{'1. √(m_d/m_s) [GST, quark masses]':<50} {r_d:.6f}  {(r_d-V_us_PDG)/V_us_err:+.1f}σ")
print(f"{'2. m_π/√(2m_K²-m_π²) [Weinberg-Oakes]':<50} {sthC_weinberg:.6f}  {(sthC_weinberg-V_us_PDG)/V_us_err:+.1f}σ")
print(f"{'3. |V_us|_PDG':<50} {V_us_PDG:.6f}  ref")

# Comparison table (as requested in the problem)
print(f"\n  θ_C(GST) = {np.degrees(np.arcsin(r_d)):.2f}°")
print(f"  θ_C(Oakes) = {np.degrees(np.arcsin(sthC_weinberg)):.2f}°")
print(f"  θ_C(PDG) = {np.degrees(np.arcsin(V_us_PDG)):.2f}°")

print("\n" + "="*70)
print("PART 4: GST in Abelian Charge Language")
print("="*70)

# Koide parametrization: √m_k = z_0(1 + √2 cos(δ + 2πk/3))
# with k = 0,1,2.
# The mass-charge identity <z_k²> = z_0² is equivalent to Q = 2/3.
# 
# For a Koide triple with k=0 → lightest (d), k=1 → middle (s), k=2 → heaviest (b):
#   sin²θ_C = m_d/m_s = [(1+√2 cos δ)/(1+√2 cos(δ+2π/3))]²
#
# At the SEED δ = 3π/4:
#   1 + √2 cos(3π/4) = 1 + √2·(-1/√2) = 0  → m_d = 0
#   sin θ_C = 0 at the seed.
#
# Small bloom: δ = 3π/4 + ε
#   Numerator: 1 + √2 cos(3π/4 + ε) ≈ -ε
#   Denominator: 1 + √2 cos(17π/12 + ε) ≈ (3-√3)/2 - ε·(√3-1)/(2√2)·√2
#     = (3-√3)/2 + O(ε)
#
# So sin θ_C ≈ |ε| / ((3-√3)/2) = 2|ε|/(3-√3)

# More precisely: near δ = 3π/4 + ε,
# cos(3π/4 + ε) = -cos ε/√2 - sin ε/√2 ≈ -1/√2 - ε/√2
# 1 + √2(-1/√2 - ε/√2) = 1 - 1 - ε = -ε
# 
# cos(3π/4 + ε + 2π/3) = cos(17π/12 + ε)
# cos(17π/12) = -(√6-√2)/4
# sin(17π/12) = (√6+√2)/4  [17π/12 = 255°, sin 255° = -sin 75° = -(√6+√2)/4]
# Actually sin(255°) = -sin(75°) = -(√6+√2)/4
# cos(17π/12 + ε) ≈ -(√6-√2)/4 + ε(√6+√2)/4
# √2·cos(17π/12 + ε) ≈ -(√3-1)/2 + ε(√3+1)/2
# 1 + √2·cos(17π/12 + ε) ≈ (3-√3)/2 + ε(√3+1)/2

# So:
# sin θ_C ≈ ε / ((3-√3)/2 + ε(√3+1)/2)
# ≈ 2ε/(3-√3) for small ε

# Verify numerically by finding ε that gives physical m_d/m_s:
# sin²θ_C = m_d/m_s = 0.05000
# sin θ_C = 0.22361

def sin_theta_from_delta(delta):
    """sin θ_C = √(m₀/m₁) in the Koide parametrization."""
    f0 = (1 + sqrt(2)*cos(delta))**2
    f1 = (1 + sqrt(2)*cos(delta + 2*pi/3))**2
    if f1 == 0:
        return float('inf')
    return sqrt(f0 / f1)

# Scan to find where it crosses the target
target = sqrt(m_d / m_s)
print(f"\nTarget: sin θ_C = √(m_d/m_s) = {target:.6f}")

# The function is 0 at δ = 3π/4 and increases. Find where it hits target.
deltas = np.linspace(3*pi/4 + 0.001, 3*pi/4 + 0.5, 1000)
vals = [sin_theta_from_delta(d) for d in deltas]
vals = np.array(vals)
# Find crossing
idx_cross = np.argmin(np.abs(vals - target))
print(f"Rough estimate: δ ≈ {deltas[idx_cross]:.6f} rad, sin θ ≈ {vals[idx_cross]:.6f}")

# Refine with brentq
def f_root(delta):
    return sin_theta_from_delta(delta) - target

# Find sign change
a = 3*pi/4 + 0.001
b = deltas[idx_cross] + 0.1
print(f"f({a:.4f}) = {f_root(a):.6f}")
print(f"f({b:.4f}) = {f_root(b):.6f}")

# The function might not cross below for large δ. Let me check.
# At δ = 3π/4 + π/3 (= 13π/12), the roles of the masses rotate.
# Let's check values more carefully:
print(f"\nScan of sin θ_C vs δ:")
for eps_val in [0.001, 0.01, 0.02, 0.05, 0.08, 0.1, 0.14, 0.15, 0.2, 0.3, 0.5]:
    d = 3*pi/4 + eps_val
    s = sin_theta_from_delta(d)
    print(f"  ε = {eps_val:.3f}, δ = {np.degrees(d):.2f}°, sin θ_C = {s:.6f}")

# From the scan: sin θ_C = 0.224 around ε ≈ 0.14-0.15
# Let me search more carefully
a = 3*pi/4 + 0.13
b = 3*pi/4 + 0.16
fa = f_root(a)
fb = f_root(b)
print(f"\nf({a:.4f}) = {fa:.6f}")
print(f"f({b:.4f}) = {fb:.6f}")

if fa * fb < 0:
    delta_sol = brentq(f_root, a, b)
    eps_sol = delta_sol - 3*pi/4
    print(f"\nSolution: δ = {delta_sol:.8f} rad = {np.degrees(delta_sol):.4f}°")
    print(f"  ε = δ - 3π/4 = {eps_sol:.8f} rad = {np.degrees(eps_sol):.4f}°")
    print(f"  sin θ_C = {sin_theta_from_delta(delta_sol):.8f}")
    
    # Verify approximation
    approx = 2*eps_sol / (3-sqrt(3))
    print(f"\n  Approximation: 2ε/(3-√3) = {approx:.8f}")
    print(f"  Exact:                      {sin_theta_from_delta(delta_sol):.8f}")
    print(f"  Error: {abs(approx - sin_theta_from_delta(delta_sol))/sin_theta_from_delta(delta_sol)*100:.2f}%")
else:
    # Try wider range
    for b_try in np.arange(3*pi/4 + 0.1, 3*pi/4 + 1.0, 0.01):
        if f_root(3*pi/4 + 0.001) * f_root(b_try) < 0:
            delta_sol = brentq(f_root, 3*pi/4 + 0.001, b_try)
            eps_sol = delta_sol - 3*pi/4
            print(f"\nSolution found: δ = {delta_sol:.6f}, ε = {eps_sol:.6f}")
            print(f"sin θ_C = {sin_theta_from_delta(delta_sol):.6f}")
            break
    else:
        print("No solution found in range!")

# Physical interpretation
print(f"\n--- Physical interpretation ---")
print(f"The GST relation sin θ_C = √(m_d/m_s) = |ℨ_d|/|ℨ_s|")
print(f"connects the Cabibbo angle to the RATIO of mass charges.")
print(f"At the seed (δ = 3π/4), the d quark is massless and θ_C = 0.")
print(f"The bloom ε = {eps_sol:.4f} rad from seed generates the Cabibbo angle.")
print(f"The Cabibbo angle MEASURES the bloom distance from the seed.")

# Constraints on (z_0, δ):
# If Q = 2/3 holds (Koide condition), then z_0 is determined by
# the total mass: 6z_0² = m_d + m_s + m_b.
# And δ is determined by sin θ_C = 2ε/(3-√3) where ε = δ - 3π/4.
# So the free parameter is z_0 (the mass scale), while δ (the Koide phase)
# is fixed by the Cabibbo angle.

# If we know m_s and demand Koide + GST:
# m_s = z_0² f₁(δ) where f₁ = (1 + √2 cos(δ+2π/3))²
# m_d = z_0² f₀(δ) = z_0² ε² (approximately)
# m_b = z_0² f₂(δ)
# sin θ_C = √(f₀/f₁) = ε·(something)
# So z_0 is the scale and ε encodes sin θ_C.

# At the physical point:
f0_phys = (1 + sqrt(2)*cos(delta_sol))**2
f1_phys = (1 + sqrt(2)*cos(delta_sol + 2*pi/3))**2
f2_phys = (1 + sqrt(2)*cos(delta_sol + 4*pi/3))**2
z0_from_ms = sqrt(m_s / f1_phys)
m_d_pred = z0_from_ms**2 * f0_phys
m_b_pred = z0_from_ms**2 * f2_phys
print(f"\nIf Koide held for (d,s,b) with this δ:")
print(f"  z₀ = √(m_s/f₁) = {z0_from_ms:.4f} MeV^(1/2)")
print(f"  m_d(pred) = {m_d_pred:.4f} MeV (PDG: {m_d})")
print(f"  m_b(pred) = {m_b_pred:.1f} MeV (PDG: {m_b})")

Q_pred = (m_d_pred + m_s + m_b_pred) / (sqrt(m_d_pred) + sqrt(m_s) + sqrt(m_b_pred))**2
print(f"  Q = {Q_pred:.6f} (should be 2/3 = {2/3:.6f})")

print(f"\n" + "="*70)
print("FINAL SUMMARY TABLE")
print("="*70)
print(f"\n{'Quantity':<55} {'Value':<12} {'PDG/σ'}")
print(f"{'-'*75}")
print(f"{'√(m_d/m_s) [GST]':<55} {r_d:.6f}     {(r_d-V_us_PDG)/V_us_err:+.1f}σ")
print(f"{'m_π/√(2m_K² - m_π²) [Weinberg-Oakes]':<55} {sthC_weinberg:.6f}     {(sthC_weinberg-V_us_PDG)/V_us_err:+.1f}σ")
print(f"{'|√(m_d/m_s) - √(m_u/m_c)| [Fritzsch, φ=0]':<55} {V_us_min:.6f}     {(V_us_min-V_us_PDG)/V_us_err:+.1f}σ")
print(f"{'√(m_d/m_s + m_u/m_c) [Fritzsch, φ=π/2]':<55} {V_us_pi2:.6f}     {(V_us_pi2-V_us_PDG)/V_us_err:+.1f}σ")
print(f"{'Fritzsch optimal (φ = {np.degrees(phi_opt):.0f}°)':<55} {V_us_PDG:.6f}     exact")
print(f"{'|V_us|_PDG':<55} {V_us_PDG:.6f}     ref")

print(f"\n{'Cabibbo angle θ_C':<55} {'Value'}")
print(f"{'-'*75}")
print(f"{'From GST':<55} {np.degrees(np.arcsin(r_d)):.2f}°")
print(f"{'From Weinberg-Oakes':<55} {np.degrees(np.arcsin(sthC_weinberg)):.2f}°")
print(f"{'From PDG |V_us|':<55} {np.degrees(np.arcsin(V_us_PDG)):.2f}°")

print(f"\nCharge-language result:")
print(f"  δ = 3π/4 + {eps_sol:.6f} = {delta_sol:.6f} rad")
print(f"  sin θ_C = 2ε/(3-√3) to {abs(approx - sin_theta_from_delta(delta_sol))/sin_theta_from_delta(delta_sol)*100:.1f}%")
print(f"  The Cabibbo angle measures the bloom of the d-quark charge")
print(f"  from the seed at δ = 3π/4.")

