#!/usr/bin/env python3
"""
CKM matrix from the sBootstrap Lagrangian.

We have:
  W = Tr(m̂ M) + X(det M^{(3)} - BB̃ - Λ^6) + c_3(det M^{(3)})^3/Λ^18
      + y_c H_u M^c_c + y_b H_d M^b_b + y_t H_u Q^t Q̄_t
      + λ_S S H_u·H_d + κ/3 S^3

Question: Does the STRUCTURE of this Lagrangian (seesaw + Seiberg constraint + det M)
generate CKM mixing beyond what is put in by hand through off-diagonal Yukawas?

We analyze three scenarios:
  A) Diagonal Yukawas only → CKM = identity
  B) Seesaw-induced one-loop off-diagonal meson VEVs from X det M coupling
  C) Off-diagonal Yukawa couplings Y^u_{cj}, Y^d_{bj} with j ≠ c,b
  D) Fritzsch-type texture as UV input (for comparison)
"""

import numpy as np
from numpy import sqrt, log, pi

# ============================================================
# INPUT PARAMETERS
# ============================================================
# Quark masses (MeV) - running masses at ~2 GeV for light, m_q for heavy
m_u = 2.16      # MeV
m_d = 4.67      # MeV
m_s = 93.4      # MeV
m_c = 1275.0    # MeV
m_b = 4180.0    # MeV
m_t = 172500.0  # MeV

# Electroweak parameters
v = 246000.0    # MeV (246 GeV)
v_u = v / sqrt(2)  # tan β = 1
v_d = v / sqrt(2)

# QCD parameters
Lambda = 1000.0  # MeV (Λ = 1 GeV)
alpha_s = 0.3

# PDG CKM values
V_us_pdg = 0.2243
V_cb_pdg = 0.0422
V_ub_pdg = 0.00394

print("=" * 80)
print("FULL CKM ANALYSIS FROM THE LAGRANGIAN")
print("=" * 80)

# ============================================================
# SCENARIO A: DIAGONAL YUKAWAS
# ============================================================
print("\n" + "=" * 80)
print("SCENARIO A: Strictly diagonal Yukawas")
print("=" * 80)

# ISS seesaw for confined sector (u, d, s)
# ⟨M^i_i⟩ = C / m_i where C = (Λ^6 m_u m_d m_s)^{1/3}
C_confined = (Lambda**6 * m_u * m_d * m_s)**(1.0/3.0)
print(f"\nConfined sector scale: C = (Λ^6 m_u m_d m_s)^{{1/3}} = {C_confined:.2f} MeV^2")
print(f"  = ({C_confined/Lambda**2:.4f}) Λ^2")

# Meson VEVs
M_uu = C_confined / m_u
M_dd = C_confined / m_d
M_ss = C_confined / m_s

print(f"\nDiagonal meson VEVs (confined sector):")
print(f"  ⟨M^u_u⟩ = C/m_u = {M_uu:.2f} MeV")
print(f"  ⟨M^d_d⟩ = C/m_d = {M_dd:.2f} MeV")
print(f"  ⟨M^s_s⟩ = C/m_s = {M_ss:.2f} MeV")

# For c, b: similar seesaw from Tr(m̂ M) term
# The mass term gives ⟨M^c_c⟩ = C'/m_c, ⟨M^b_b⟩ = C'/m_b
# where C' depends on the specific dynamics. For the diagonal Yukawa case,
# the physical masses are reproduced by construction:
# m_c^phys = y_c v_u ⟨M^c_c⟩ = y_c v_u C'/m_c
# This is just parameter fitting — y_c C' is chosen to give the right mass.

print(f"\nWith diagonal Yukawas:")
print(f"  M_U = diag(m_u, m_c, m_t) = diag({m_u}, {m_c}, {m_t}) MeV")
print(f"  M_D = diag(m_d, m_s, m_b) = diag({m_d}, {m_s}, {m_b}) MeV")
print(f"  CKM = I (identity matrix)")
print(f"\n  |V_us| = 0     (PDG: {V_us_pdg})")
print(f"  |V_cb| = 0     (PDG: {V_cb_pdg})")
print(f"  |V_ub| = 0     (PDG: {V_ub_pdg})")
print(f"\n  VERDICT: Diagonal Yukawas give ZERO CKM mixing. Trivially.")

# ============================================================
# SCENARIO B: ONE-LOOP SEESAW-INDUCED OFF-DIAGONAL VEVs
# ============================================================
print("\n" + "=" * 80)
print("SCENARIO B: One-loop off-diagonal meson VEVs from X det M^{(3)}")
print("=" * 80)

# The Seiberg constraint term: X det M^{(3)}
# det M^{(3)} = M^u_u M^d_d M^s_s - M^u_d M^d_u M^s_s - M^u_s M^s_u M^d_d
#             - M^d_s M^s_d M^u_u + M^u_d M^d_s M^s_u + M^u_s M^s_d M^d_u
#
# F-terms from W = Tr(m̂ M) + X det M^{(3)}:
#   F_{M^j_i} = m̂_i δ^j_i + X · cofactor(M)^j_i
#
# For diagonal M at tree level, the cofactor matrix is:
#   cof(M)^u_u = M^d_d M^s_s, etc.
#   cof(M)^u_d = -M^d_u M^s_s = 0 (since M^d_u = 0 at tree level)
#
# So at tree level, F_{M^j_i} = 0 for i ≠ j. The off-diagonal mesons are flat.

print("\nTree-level analysis:")
print("  det M^{(3)} = M^u_u M^d_d M^s_s (at diagonal VEV)")
print("  cofactor(M)^i_j for i≠j = 0 (at diagonal VEV)")
print("  → No tree-level off-diagonal meson VEVs")

# Now compute the mass matrix for off-diagonal mesons.
# The second derivatives of W give the mass terms for the off-diagonal mesons.
#
# ∂²W / ∂M^u_d ∂M^d_u = X · ∂²(det M)/∂M^u_d ∂M^d_u
#
# For det M = M^u_u M^d_d M^s_s - M^u_d M^d_u M^s_s - ...
# ∂²(det M)/∂M^d_u ∂M^u_d = -M^s_s  (from the -M^u_d M^d_u M^s_s term)
#
# Wait — careful with indices. det M^{(3)} for M^i_j:
# det = ε_{ijk} ε^{lmn} M^i_l M^j_m M^k_n / 6 (with i,j,k ∈ {u,d,s}, l,m,n ∈ {u,d,s})
#
# But in the ISS, M^i_j = q^i q̄_j, so ∂W/∂M^j_i gives the coupling.
#
# The off-diagonal meson mass matrix (for the pair M^u_d, M^d_u) comes from:
# W contains: X · det M^{(3)}
# The bilinear in off-diagonal mesons:
#   X · (-M^s_s) · M^u_d · M^d_u   (from det expansion)
# So the mass term is: m²_{M^u_d, M^d_u} = ⟨X⟩ · ⟨M^s_s⟩
#
# But what is ⟨X⟩? At the ISS vacuum:
# F_X = det M^{(3)} - BB̃ - Λ^6 = 0 (Seiberg constraint satisfied)
# The VEV of X is determined by the F-term conditions for M:
# F_{M^u_u} = m_u + X · ⟨M^d_d⟩⟨M^s_s⟩ = 0
# → ⟨X⟩ = -m_u / (⟨M^d_d⟩ · ⟨M^s_s⟩)

X_vev = -m_u / (M_dd * M_ss)
print(f"\nX VEV from F_{{M^u_u}} = 0:")
print(f"  ⟨X⟩ = -m_u / (⟨M^d_d⟩ · ⟨M^s_s⟩) = {X_vev:.6e} MeV^{{-1}}")

# Cross-check: F_{M^d_d} = m_d + X · ⟨M^u_u⟩ · ⟨M^s_s⟩ should also = 0
X_check_d = -m_d / (M_uu * M_ss)
X_check_s = -m_s / (M_uu * M_dd)
print(f"  Cross-check from F_{{M^d_d}}: ⟨X⟩ = {X_check_d:.6e} MeV^{{-1}}")
print(f"  Cross-check from F_{{M^s_s}}: ⟨X⟩ = {X_check_s:.6e} MeV^{{-1}}")

# These should all be equal. Let's check:
# ⟨X⟩ = -m_u/(M_dd · M_ss) = -m_u · m_d · m_s / (C² · 1) = -m_u m_d m_s / C²
# C² = (Λ^6 m_u m_d m_s)^{2/3}
# So ⟨X⟩ = -(m_u m_d m_s) / (Λ^6 m_u m_d m_s)^{2/3} = -(m_u m_d m_s)^{1/3} / Λ^4
X_analytic = -(m_u * m_d * m_s)**(1.0/3) / Lambda**4
print(f"\n  Analytic: ⟨X⟩ = -(m_u m_d m_s)^{{1/3}} / Λ^4 = {X_analytic:.6e} MeV^{{-1}}")
print(f"  This is the SAME for all three F-term equations (self-consistent).")

# Off-diagonal meson masses
# For the pair (M^u_d, M^d_u):
# W ⊃ X · (-⟨M^s_s⟩) · M^u_d · M^d_u + m_u M^u_u + m_d M^d_d + ...
# But wait: the mass term Tr(m̂ M) contributes m_u M^u_u, m_d M^d_d but NOT
# m_u M^u_d or m_d M^d_u (because m̂ is diagonal).
# The ONLY mass term for off-diagonal mesons comes from the X det M coupling.
#
# Mass of (M^u_d, M^d_u) pair: |⟨X⟩| · |⟨M^s_s⟩| = |X| · C/m_s
# Mass of (M^u_s, M^s_u) pair: |⟨X⟩| · |⟨M^d_d⟩| = |X| · C/m_d
# Mass of (M^d_s, M^s_d) pair: |⟨X⟩| · |⟨M^u_u⟩| = |X| · C/m_u

m2_ud = abs(X_vev) * M_ss
m2_us = abs(X_vev) * M_dd
m2_ds = abs(X_vev) * M_uu

print(f"\nOff-diagonal meson masses (from X det M coupling):")
print(f"  m(M^u_d, M^d_u) = |⟨X⟩| · ⟨M^s_s⟩ = {m2_ud:.6e} MeV")
print(f"  m(M^u_s, M^s_u) = |⟨X⟩| · ⟨M^d_d⟩ = {m2_us:.6e} MeV")
print(f"  m(M^d_s, M^s_d) = |⟨X⟩| · ⟨M^u_u⟩ = {m2_ds:.6e} MeV")

# Simplify: |⟨X⟩| · ⟨M^k_k⟩ = (m_u m_d m_s)^{1/3}/Λ^4 · C/m_k
#                                = (m_u m_d m_s)^{1/3}/Λ^4 · (Λ^6 m_u m_d m_s)^{1/3}/m_k
#                                = (m_u m_d m_s)^{2/3} Λ^2 / (Λ^4 m_k)
#                                = (m_u m_d m_s)^{2/3} / (Λ^2 m_k)
m2_ud_check = (m_u * m_d * m_s)**(2.0/3) / (Lambda**2 * m_s)
m2_us_check = (m_u * m_d * m_s)**(2.0/3) / (Lambda**2 * m_d)
m2_ds_check = (m_u * m_d * m_s)**(2.0/3) / (Lambda**2 * m_u)
print(f"\n  Simplified: m_{{ud}} = (m_u m_d m_s)^{{2/3}} / (Lam^2 m_s) = {m2_ud_check:.6e} MeV")
print(f"              m_{{us}} = (m_u m_d m_s)^{{2/3}} / (Lam^2 m_d) = {m2_us_check:.6e} MeV")
print(f"              m_{{ds}} = (m_u m_d m_s)^{{2/3}} / (Lam^2 m_u) = {m2_ds_check:.6e} MeV")

# These are VERY small (because they scale as m_quark/Λ).
# The off-diagonal mesons are NOT massless but very light.

# NOW: Does the Coleman-Weinberg effective potential generate off-diagonal VEVs?
# The CW potential at one loop:
#   V_CW = (1/64π²) STr[M⁴ log(M²/μ²)]
# where M is the mass matrix in the theory.
#
# The key question: is there a tadpole for off-diagonal mesons?
# At the diagonal VEV, the tadpole is:
#   ∂V_CW/∂M^u_d |_{M^u_d=0} = ?
#
# By the symmetry of the problem: the Lagrangian with diagonal m̂ and diagonal ⟨M⟩
# has a Z_2 symmetry M^i_j → -M^i_j for i ≠ j (off-diagonal parity).
# This means the tadpole VANISHES: ∂V_CW/∂M^i_j|_{diag} = 0 for i ≠ j.

print("\n" + "-" * 60)
print("One-loop Coleman-Weinberg analysis:")
print("-" * 60)
print("  The diagonal vacuum ⟨M⟩ = diag has a Z₂ symmetry:")
print("    M^i_j → -M^i_j for i ≠ j (off-diagonal sign flip)")
print("  This symmetry is preserved by:")
print("    - Tr(m̂ M) (m̂ diagonal)")
print("    - det M^{(3)} (determinant invariant under sign flip of off-diag)")
print("    - All F-terms")
print("  Therefore: ∂V_CW/∂M^i_j|_{diag} = 0 for i ≠ j")
print("  → NO off-diagonal meson VEVs generated at one loop")

# But wait — could two-loop effects or non-perturbative effects break this?
# Let's check if there's a SECOND-ORDER effect: the off-diagonal meson
# gets a VEV proportional to some small symmetry-breaking parameter.
#
# The only way to get off-diagonal VEVs is if m̂ is NOT diagonal, i.e.,
# if the quarks are not in the mass basis. But by definition, m̂ is the
# mass matrix IN the mass basis.

print("\n  The Z₂ symmetry is EXACT at all loop orders for diagonal m̂.")
print("  Off-diagonal meson VEVs require off-diagonal entries in m̂,")
print("  which would mean the quarks are not in the mass basis.")
print("\n  VERDICT: The seesaw + Seiberg constraint generates NO CKM mixing")
print("  at any order in perturbation theory, as long as the UV mass matrix")
print("  m̂ is diagonal.")

# ============================================================
# SCENARIO C: OFF-DIAGONAL YUKAWA COUPLINGS
# ============================================================
print("\n" + "=" * 80)
print("SCENARIO C: Off-diagonal Yukawa couplings")
print("=" * 80)

# The most general Yukawa consistent with gauge symmetry:
#   W_Y = Σ_j Y^u_{cj} H_u M^c_j + Σ_j Y^d_{bj} H_d M^b_j + y_t H_u Q^t Q̄_t
#
# After EWSB with ⟨H_u⟩ = v_u, ⟨H_d⟩ = v_d:
# The up-type mass matrix in the basis (u, c, t):
#   M_U = | m_u^{eff}     Y^u_{cu} v_u ⟨M^c_u⟩     0           |
#         | 0              Y^u_{cc} v_u ⟨M^c_c⟩     0           |
#         | 0              0                         y_t v_u      |
#
# Wait — this isn't right. Let me think more carefully.
#
# The light quarks (u, d, s) get their masses from the ISS seesaw.
# In the ISS, the physical quarks are NOT the original quarks q^i.
# They are linear combinations that diagonalize the seesaw mass matrix.
#
# Actually, in the ISS/Seiberg duality:
# - The ELECTRIC quarks q^i, q̄_j are confined
# - The MAGNETIC quarks χ^i, χ̄_j are the light degrees of freedom
# - The mesons M^i_j = q^i q̄_j are also light
# - The Yukawa coupling Y H M couples the Higgs to the meson
# - The meson VEV gives mass to the quarks through the ISS mechanism
#
# The FULL mass matrix for quarks comes from:
# W ⊃ h χ M χ̄ + μ χ χ̄  (ISS superpotential)
# After M gets a VEV, the quark mass is m_χ ~ h ⟨M⟩ ~ h C/m_q
#
# But the PHYSICAL SM quark mass comes from the Yukawa:
# W ⊃ Y^u_{ij} H_u M^i_j → after EWSB: Y^u_{ij} v_u ⟨M^i_j⟩
#
# For DIAGONAL Yukawas and diagonal ⟨M⟩:
#   m_u^{phys} ~ Y^u_{uu} v_u C/m_u
#   m_c^{phys} ~ Y^u_{cc} v_u C'/m_c
#   m_t^{phys} ~ y_t v_u
#
# The Yukawa couplings Y^u_{ij} are then CHOSEN to reproduce the masses.
# This is parameter fitting.
#
# For OFF-DIAGONAL Yukawas: Y^u_{cu} ≠ 0, Y^u_{cs} ≠ 0, etc.
# These would generate off-diagonal entries in the mass matrix.
# But these are FREE PARAMETERS — they're not predicted.

print("\nThe most general Yukawa:")
print("  W_Y = Σ_{j} Y^u_{cj} H_u M^c_j + Σ_{j} Y^d_{bj} H_d M^b_j + y_t H_u Q^t Q̄_t")
print("\nAfter EWSB with diagonal meson VEVs ⟨M^i_j⟩ = (C/m_i) δ^i_j:")
print("  Only the diagonal terms survive!")
print("  Y^u_{cu} v_u ⟨M^c_u⟩ = Y^u_{cu} v_u · 0 = 0")
print("  Y^u_{cc} v_u ⟨M^c_c⟩ = Y^u_{cc} v_u C'/m_c ≠ 0")
print("\n  Even with off-diagonal Yukawas, the diagonal meson VEVs")
print("  project out all off-diagonal contributions!")

print("\n  HOWEVER: if there were also off-diagonal terms like")
print("  Y^u_{uc} H_u M^u_c, these would contribute Y^u_{uc} v_u ⟨M^u_c⟩ = 0")
print("  (again zero because ⟨M^u_c⟩ = 0).")

print("\n  The ONLY way to get off-diagonal mass entries is:")
print("  1. Off-diagonal meson VEVs ⟨M^i_j⟩ ≠ 0 for i ≠ j (ruled out in Scenario B)")
print("  2. Direct Yukawa couplings like Y H_u q^c q̄_u (but these don't exist in")
print("     the magnetic dual — the elementary quarks are confined)")
print("  3. UV physics that makes m̂ non-diagonal")

print("\n  VERDICT: Off-diagonal Yukawas multiplied by diagonal meson VEVs = 0.")
print("  CKM mixing CANNOT come from the IR meson sector.")

# ============================================================
# SCENARIO B': WHAT IF m̂ HAS OFF-DIAGONAL ENTRIES?
# ============================================================
print("\n" + "=" * 80)
print("SCENARIO B': Non-diagonal mass matrix m̂ (UV input)")
print("=" * 80)

# If the UV mass matrix is:
#   m̂ = | m_uu  m_ud  m_us |
#       | m_du  m_dd  m_ds |
#       | m_su  m_sd  m_ss |
# Then the Seiberg constraint F_{M^j_i} = m̂_{ij} + X · cof(M)^j_i = 0
# gives ⟨M⟩ = -m̂^{-1} / X, i.e., M is proportional to the INVERSE of m̂.
#
# The inverse of a non-diagonal matrix IS non-diagonal.
# So off-diagonal ⟨M^i_j⟩ ≠ 0, and the Yukawa coupling Y H M generates mixing.

print("\nIf m̂ is non-diagonal (Fritzsch texture or similar):")
print("  F_{M^j_i} = m̂_{ij} + X · cof(M)^j_i = 0")
print("  → ⟨M⟩ ∝ m̂^{-1}  (matrix inverse)")
print("  → Off-diagonal meson VEVs generated!")
print("  → CKM mixing from the inverse mass matrix structure")

# Let's compute: if m̂ has a Fritzsch texture, what is ⟨M⟩?
# Fritzsch texture (6-zero):
#   m_u_matrix = | 0    A_u  0   |     m_d_matrix = | 0    A_d  0   |
#                | A_u* 0    B_u |                   | A_d* 0    B_d |
#                | 0    B_u* C_u |                   | 0    B_d* C_d |
#
# with |A_u| = √(m_u m_c), |B_u| = √(m_c m_t), |C_u| ≈ m_t (dominant)
# and  |A_d| = √(m_d m_s), |B_d| = √(m_s m_b), |C_d| ≈ m_b

# But wait — in our setup, only (u,d,s) are confined. c and b are separate.
# So the 3×3 ISS sector has flavors (u,d,s) only.
# The relevant mass matrix for the confined sector is:
#   m̂_confined = | m_uu  m_ud  m_us |
#                | m_du  m_dd  m_ds |
#                | m_su  m_sd  m_ss |
# for (u, d, s) flavors.

# In the Fritzsch texture for 3 generations, mixing comes from:
# WITHIN the up sector: u-c mixing (but c is NOT in the confined sector!)
# WITHIN the down sector: d-s mixing

# So the Fritzsch texture for the confined (u,d,s) sector would be:
# m̂_confined = diag(m_u, m_d, m_s) + off-diagonal terms
# But the off-diagonal terms mix u↔d, u↔s, d↔s — which are different-charge quarks!
# That makes no sense for a mass matrix. Mass matrices are WITHIN a charge sector.

print("\n  KEY INSIGHT: The confined sector (u,d,s) mixes DIFFERENT charge quarks!")
print("  u is charge +2/3, d and s are charge -1/3")
print("  The mass matrix m̂ in Tr(m̂ M) is flavor-diagonal BY CONSTRUCTION")
print("  because it comes from the original QCD mass terms m_i q̄^i q_i.")
print("  Off-diagonal terms m_{ud} q̄^u q_d would violate electric charge!")

print("\n  In the ISS, the meson M^i_j = q^i q̄_j carries charge Q_i - Q_j.")
print("  M^u_d has charge +2/3 - (-1/3) = +1 (charged meson → π⁺-like)")
print("  M^u_u has charge 0 (neutral)")
print("  The mass term Tr(m̂ M) must be charge-neutral → m̂ diagonal!")

print("\n  THEREFORE: m̂ cannot have off-diagonal entries between different-charge")
print("  quarks. The ISS mass matrix is necessarily diagonal in the mass basis.")
print("  CKM = identity from the confined sector alone.")

# ============================================================
# THE ACTUAL CKM MECHANISM: SEPARATE UP AND DOWN SECTORS
# ============================================================
print("\n" + "=" * 80)
print("THE ACTUAL STRUCTURE: Separate up and down mass matrices")
print("=" * 80)

print("""
In the Standard Model, the CKM matrix arises because the up-type and down-type
mass matrices are diagonalized by DIFFERENT unitary transformations:
  V_CKM = V_U† V_D

In our Lagrangian:
  UP sector (charge +2/3):   u, c, t
  DOWN sector (charge -1/3): d, s, b

The ISS confines (u, d, s) together, but the MASS TERMS are:
  m_u M^u_u + m_d M^d_d + m_s M^s_s  (all diagonal, charge-preserving)

The up-type 3×3 mass matrix has:
  Row 1 (u): mass from ISS seesaw → m_u
  Row 2 (c): mass from y_c H_u M^c_c → m_c  (M^c_c is NOT in the confined sector)
  Row 3 (t): mass from y_t H_u Q^t Q̄_t → m_t

The down-type 3×3 mass matrix has:
  Row 1 (d): mass from ISS seesaw → m_d
  Row 2 (s): mass from ISS seesaw → m_s
  Row 3 (b): mass from y_b H_d M^b_b → m_b

CKM mixing requires off-diagonal entries WITHIN the up-type matrix
or WITHIN the down-type matrix.

Can the ISS generate such entries?
""")

# The up-type mass matrix connects:
# u (confined, ISS seesaw) with c (composite meson M^c_c) and t (elementary)
# The ONLY way to get u-c mixing is if there's a coupling between the
# ISS sector (that produces u) and the meson M^c_c.
#
# In the superpotential, M^c_u = q^c q̄_u is a meson. If c is not in the
# confined sector but couples to it, then M^c_u could get a VEV.
# But M^c_u is a (anti)quark bilinear with charge +2/3 - (+2/3) = 0.
# So charge conservation allows it!

print("Charge analysis for up-type off-diagonal mesons:")
print("  M^c_u = c̄ u (or q^c q̄_u in sBootstrap): charge 0 ✓")
print("  This IS allowed by charge conservation!")
print()
print("The key question: does M^c_u get a VEV?")
print()
print("In the ISS, M^c_u = q^c q̄_u is a meson where c is NOT in the")
print("confined sector. So M^c_u is NOT one of the mesons in det M^{(3)}.")
print("The F-term for M^c_u is:")
print("  F_{M^u_c} = ∂W/∂M^u_c = m̂_{cu} = 0  (no off-diagonal mass term)")
print("  (because m̂ = diag(m_u, m_d, m_s, m_c, m_b) has no c-u entry)")
print()
print("The meson M^c_u has NO potential at tree level — it's a flat direction.")
print("Its VEV is zero by the minimum principle (no tadpole, positive mass²).")

# ============================================================
# ONE-LOOP CONTRIBUTION TO OFF-DIAGONAL MESON VEVs
# ============================================================
print("\n" + "=" * 80)
print("ONE-LOOP: Coleman-Weinberg contribution to off-diagonal VEVs")
print("=" * 80)

# The CW potential involves loops of the massive fields.
# For M^c_u (charge-neutral, up-type off-diagonal):
#
# The mass of M^c_u comes from the Kähler potential K ⊃ |M^c_u|²
# and any terms in W that couple to M^c_u.
#
# In the superpotential:
# W ⊃ m_c M^c_c + m_u M^u_u + ... (no M^c_u term)
# W ⊃ X det M^{(3)} (M^c_u is NOT in M^{(3)}, which only contains u,d,s)
#
# So M^c_u is MASSLESS at tree level (flat direction).
# At one loop, the CW potential generates a mass for M^c_u through
# diagrams involving the massive ISS modes and the meson M^c_c.

# The one-loop tadpole for M^c_u:
# For a tadpole to exist, we need a coupling of M^c_u to other fields
# that have VEVs. The coupling must be:
# W ⊃ something · M^c_u · (field with VEV)
#
# Looking at W: there is NO such coupling. The only terms involving M^c_j
# are m_c M^c_c (diagonal) and possibly Y^u_{cu} H_u M^c_u (if we allow it).
#
# With Y^u_{cu} H_u M^c_u: after EWSB, this gives a tadpole
# Y^u_{cu} v_u M^c_u, but this is a TREE-LEVEL mass term for M^c_u,
# not a VEV for M^c_u.

# Actually wait — Y^u_{cu} H_u M^c_u is a coupling between H_u and M^c_u.
# After EWSB with ⟨H_u⟩ = v_u, this gives a LINEAR term Y^u_{cu} v_u M^c_u
# in the effective potential. This would be a TADPOLE for M^c_u!
# If such a coupling exists with Y^u_{cu} ≠ 0, then M^c_u would get a VEV:
# ⟨M^c_u⟩ = -Y^u_{cu} v_u / m²_{M^c_u}

# But the mass of M^c_u is... what? From the Kähler potential, m² ~ Λ²
# or from the ISS, m² ~ h² μ² where h is the magnetic Yukawa and μ is the ISS scale.

# Actually, THIS is the key question. Let me be more careful.

print("""
The off-diagonal meson M^c_u has the following properties:
  - Charge: 0 (allowed)
  - Tree-level mass: comes from Kähler potential
  - Tree-level tadpole: ZERO (no off-diagonal mass term in m̂)
  - One-loop tadpole: ZERO (Z₂ symmetry M^c_u → -M^c_u is preserved)

The Z₂ symmetry argument:
  The Lagrangian is invariant under M^c_u → -M^c_u (and M^u_c → -M^u_c)
  because:
    - Tr(m̂ M) has only diagonal terms → invariant
    - det M^{(3)} doesn't involve M^c_u → invariant
    - y_c H_u M^c_c doesn't involve M^c_u → invariant
    - y_t H_u Q^t Q̄_t doesn't involve M^c_u → invariant

  This Z₂ is EXACT at all orders → ⟨M^c_u⟩ = 0 exactly.
""")

# ============================================================
# THE HONEST CONCLUSION
# ============================================================
print("=" * 80)
print("HONEST CONCLUSION: CKM = IDENTITY from this Lagrangian")
print("=" * 80)

print("""
The Lagrangian as written:
  W = Tr(m̂ M) + X(det M^{(3)} - BB̃ - Λ^6) + c₃(det M^{(3)})³/Λ¹⁸
      + y_c H_u M^c_c + y_b H_d M^b_b + y_t H_u Q^t Q̄_t
      + λ_S S H_u·H_d + κ/3 S³

produces CKM = I (identity matrix) because:

1. DIAGONAL m̂: The mass matrix m̂ = diag(m_u, m_d, m_s, m_c, m_b) is diagonal.
   This is not a choice — it's forced by charge conservation in the confined sector,
   which mixes u (charge +2/3) with d,s (charge -1/3). Off-diagonal entries
   would violate electric charge.

2. DIAGONAL MESON VEVs: The seesaw gives ⟨M^i_j⟩ = (C/m_i) δ^i_j.
   Off-diagonal VEVs are forbidden by the Z₂ parity M^i_j → -M^i_j (i ≠ j),
   which is exact at all orders.

3. DIAGONAL YUKAWAS: The Yukawa couplings y_c H_u M^c_c and y_b H_d M^b_b
   are diagonal by construction. Even if off-diagonal Yukawas Y^u_{cj} H_u M^c_j
   were added, they would multiply ⟨M^c_j⟩ = 0 for j ≠ c and give zero.

4. NO LOOP CORRECTIONS: The Z₂ parity prevents off-diagonal VEVs at all loop
   orders. The CW potential respects this symmetry.

5. THE ONLY WAY TO GET CKM: The mass matrix m̂ must have off-diagonal entries
   in the SAME CHARGE SECTOR. This means:
     - Up sector: m̂_up for (u, c, t) — but u is confined while c is not,
       and t is elementary. Different dynamics for each.
     - Down sector: m̂_down for (d, s, b) — d,s are confined while b is not.

   These off-diagonal entries are UV INPUT. They cannot be generated by the
   IR dynamics (ISS, seesaw, Seiberg constraint).

THEREFORE: CKM mixing in this framework is a UV input, not a prediction.
This is consistent with the memory note "CKM = UV input (Round 16D)".
""")

# ============================================================
# QUANTITATIVE: What UV input reproduces the observed CKM?
# ============================================================
print("=" * 80)
print("QUANTITATIVE: UV input needed to reproduce CKM")
print("=" * 80)

# If we ADD off-diagonal mass terms in the up and down sectors:
# M_U = | m_u     a_u    0   |     M_D = | m_d     a_d    0   |
#       | a_u'    m_c    b_u |           | a_d'    m_s    b_d |
#       | 0       b_u'   m_t |           | 0       b_d'   m_b |
# (Fritzsch-like texture)

# For the Fritzsch texture (Hermitian):
# |a_u| = √(m_u m_c), |b_u| = √(m_c m_t), diagonal ≈ m_t
# |a_d| = √(m_d m_s), |b_d| = √(m_s m_b), diagonal ≈ m_b

def fritzsch_texture(m1, m2, m3, phase_a=0, phase_b=0):
    """Construct Fritzsch Hermitian mass matrix."""
    a = np.sqrt(m1 * m2) * np.exp(1j * phase_a)
    b = np.sqrt(m2 * m3) * np.exp(1j * phase_b)
    M = np.array([
        [0,          a,           0],
        [np.conj(a), 0,           b],
        [0,          np.conj(b),  m3]
    ])
    return M

def diagonalize_hermitian(M):
    """Diagonalize M M†, return eigenvalues and unitary V."""
    MM = M @ M.conj().T
    eigenvalues, V = np.linalg.eigh(MM)
    # Sort by magnitude
    idx = np.argsort(eigenvalues)
    eigenvalues = eigenvalues[idx]
    V = V[:, idx]
    return np.sqrt(np.abs(eigenvalues)), V

# Standard Fritzsch (6-zero) texture
M_U_fritzsch = fritzsch_texture(m_u, m_c, m_t)
M_D_fritzsch = fritzsch_texture(m_d, m_s, m_b)

masses_u, V_U = diagonalize_hermitian(M_U_fritzsch)
masses_d, V_D = diagonalize_hermitian(M_D_fritzsch)

V_CKM = V_U.conj().T @ V_D

print(f"\nFritzsch 6-zero texture:")
print(f"  Up masses: {masses_u[0]:.2f}, {masses_u[1]:.2f}, {masses_u[2]:.2f} MeV")
print(f"  Down masses: {masses_d[0]:.2f}, {masses_d[1]:.2f}, {masses_d[2]:.2f} MeV")
print(f"\n  |V_CKM| =")
for i in range(3):
    print(f"    [{', '.join(f'{abs(V_CKM[i,j]):.6f}' for j in range(3))}]")

print(f"\n  |V_us| = {abs(V_CKM[0,1]):.4f}  (PDG: {V_us_pdg})")
print(f"  |V_cb| = {abs(V_CKM[1,2]):.4f}  (PDG: {V_cb_pdg})")
print(f"  |V_ub| = {abs(V_CKM[0,2]):.6f}  (PDG: {V_ub_pdg:.6f})")

# GST relation
V_us_GST = np.sqrt(m_d / m_s)
V_cb_GST = np.sqrt(m_s / m_b)
V_ub_GST = np.sqrt(m_u / m_t) * np.sqrt(m_d / m_s)  # product of two rotations

# The Fritzsch formula gives:
# |V_us| ≈ |√(m_d/m_s) - e^{iφ}√(m_u/m_c)|
# |V_cb| ≈ |√(m_s/m_b) - e^{iψ}√(m_c/m_t)|

V_us_fritzsch_approx = abs(np.sqrt(m_d/m_s) - np.sqrt(m_u/m_c))
V_cb_fritzsch_approx = abs(np.sqrt(m_s/m_b) - np.sqrt(m_c/m_t))

print(f"\n  Fritzsch approximate formulas:")
print(f"  |V_us| ≈ |√(m_d/m_s) - √(m_u/m_c)| = {V_us_fritzsch_approx:.4f}")
print(f"  |V_cb| ≈ |√(m_s/m_b) - √(m_c/m_t)| = {V_cb_fritzsch_approx:.4f}")
print(f"  GST: |V_us| ≈ √(m_d/m_s) = {V_us_GST:.4f}")

print(f"\n  Known problem: Fritzsch gives |V_cb| ≈ {V_cb_fritzsch_approx:.4f}")
print(f"  vs PDG {V_cb_pdg} — too large by ~40% (known since 1990s)")

# ============================================================
# THE BLOCK-DIAGONAL OBSTRUCTION (from memory)
# ============================================================
print("\n" + "=" * 80)
print("SUMMARY OF OBSTRUCTIONS TO CKM FROM IR DYNAMICS")
print("=" * 80)

print("""
Multiple independent arguments prove CKM cannot come from the meson sector:

1. Z₂ PARITY: Off-diagonal meson VEVs forbidden by M^i_j → -M^i_j symmetry.
   Exact at all loop orders. (This computation)

2. CHARGE CONSERVATION: The confined sector (u,d,s) mixes different charges.
   Off-diagonal mass terms m_{ud} would violate U(1)_em. (This computation)

3. BLOCK-DIAGONAL OBSTRUCTION (Round 15B): The 6×6 off-diagonal meson mass
   matrix decomposes into three independent 2×2 blocks. This structure is
   preserved by ALL polynomial Kähler invariants. CKM CANNOT come from
   meson sector.

4. FLAVOR-UNIVERSAL YUKAWA (Round 16A): y_j M_j = √2 C/v is the SAME for
   all confined quarks. Exact cancellation of m_j. No flavor discrimination
   → no mixing.

5. ALL 7 IR MECHANISMS KILLED (Round 16D): CKM comes from UV Fritzsch
   texture, transmitted by seesaw.
""")

# ============================================================
# SCENARIO D: GST + UV TEXTURE TRANSMITTED BY SEESAW
# ============================================================
print("=" * 80)
print("SCENARIO D: How the seesaw TRANSMITS a UV texture")
print("=" * 80)

print("""
If the UV mass matrix m̂ is non-diagonal (within each charge sector),
the seesaw INVERTS it:

  ⟨M^i_j⟩ ∝ (m̂^{-1})^i_j

The CKM matrix is then determined by the UV texture m̂.

For the up sector (u, c confined): if m̂_up has off-diagonal entry m̂_{uc},
then ⟨M⟩_up ∝ m̂_up^{-1} has off-diagonal entries that generate u-c mixing.

For the down sector (d, s confined): if m̂_down has off-diagonal entry m̂_{ds},
then ⟨M⟩_down ∝ m̂_down^{-1} generates d-s mixing.

The seesaw PRESERVES the mixing pattern but INVERTS the mass hierarchy.
""")

# Compute: if m̂_down has a Fritzsch structure in the (d,s) 2×2 block:
m_down_2x2 = np.array([
    [0, np.sqrt(m_d * m_s)],
    [np.sqrt(m_d * m_s), m_s]  # modified to get correct eigenvalues
])

# Actually, the standard Fritzsch 2×2 is:
# M = | 0   a |
#     | a   D |
# with eigenvalues λ± = (D ± √(D² + 4a²))/2
# For a = √(m_d m_s) and eigenvalues m_d, m_s:
# m_d + m_s = D, m_d * m_s = -a² + D*0 ...
# Actually for 2×2: eigenvalues of |0  a; a  D| are (D ± √(D²+4a²))/2
# We want them to be m_d and m_s (both positive).
# D = m_s - m_d (trace = sum of eigenvalues, but one eigenvalue is negative)
# No — let's be careful. det = -a², trace = D
# eigenvalues = (D ± √(D²+4a²))/2
# For a = √(m_d m_s): 4a² = 4 m_d m_s
# √(D²+4m_d m_s) = √((m_s-m_d)² + 4m_d m_s) if D = m_s - m_d
#                  = √(m_s²+m_d²+2m_d m_s) = m_s + m_d
# So eigenvalues = (m_s-m_d ± (m_s+m_d))/2 = m_s or -m_d
# Negative eigenvalue! The mass matrix gives -m_d.

# For M M†: eigenvalues are m_d² and m_s².
# The mixing angle is:
# tan θ = a / (m_s + m_d) ... no, need to compute properly.

# Diagonalization of Hermitian matrix M:
# M = U diag(λ₁, λ₂) U†
# For M = |0  a; a  D|:
# eigenvector for λ = m_s: (a, m_s) / √(a²+m_s²)
# eigenvector for λ = -m_d: (a, -m_d) / √(a²+m_d²)... no

# Let's just compute it numerically
print("2×2 Fritzsch texture for (d,s):")
a_ds = np.sqrt(m_d * m_s)
M_ds = np.array([[0, a_ds], [a_ds, m_s - m_d]])
evals, evecs = np.linalg.eigh(M_ds)
print(f"  M = |0     {a_ds:.2f}|")
print(f"      |{a_ds:.2f}  {m_s-m_d:.2f}|")
print(f"  Eigenvalues: {evals[0]:.4f}, {evals[1]:.4f} MeV")
print(f"  (should be close to -m_d={-m_d} and m_s={m_s})")

# The mixing angle
theta_ds = np.arctan2(evecs[0,1], evecs[1,1])
print(f"  Mixing angle theta_ds = {np.degrees(theta_ds):.2f} deg = {theta_ds:.4f} rad")
print(f"  sin(theta_ds) = {np.sin(theta_ds):.4f}")
print(f"  √(m_d/m_s) = {np.sqrt(m_d/m_s):.4f} (GST)")

# Now the SEESAW transmission:
# ⟨M⟩ ∝ m̂^{-1}
# If m̂ = U diag(m_d, m_s) U†, then m̂^{-1} = U diag(1/m_d, 1/m_s) U†
# The diagonalizing matrices are THE SAME.
# So the seesaw doesn't change the mixing angle — it just inverts the masses.

print(f"\n  SEESAW TRANSMISSION: m̂^{{-1}} has the SAME mixing angle as m̂.")
print(f"  The seesaw preserves V_{{ds}} exactly.")
print(f"  ⟨M⟩ = U · diag(C/m_d, C/m_s) · U† → same rotation U diagonalizes both.")

# ============================================================
# FULL 3×3 COMPUTATION WITH UV FRITZSCH TEXTURE
# ============================================================
print("\n" + "=" * 80)
print("FULL 3×3 COMPUTATION: Fritzsch texture as UV input")
print("=" * 80)

# Build up and down Fritzsch textures
# Up sector: (u, c, t) — but u is confined, c is composite meson, t is elementary
# Down sector: (d, s, b) — d,s are confined, b is composite meson

# For the Fritzsch texture to work here, we need:
# m̂_up = | 0      √(m_u m_c)    0          |
#         | √(m_u m_c)  0        √(m_c m_t)  |
#         | 0      √(m_c m_t)   m_t          |
#
# But u is in the confined sector (ISS) while c comes from a different mechanism.
# The off-diagonal entry √(m_u m_c) mixes an ISS quark with a composite meson quark.
# This is a UV coupling between the two sectors.

print("Up-type Fritzsch texture (UV input):")
a_u = np.sqrt(m_u * m_c)
b_u = np.sqrt(m_c * m_t)
print(f"  a_u = √(m_u · m_c) = {a_u:.2f} MeV")
print(f"  b_u = √(m_c · m_t) = {b_u:.2f} MeV")

print(f"\nDown-type Fritzsch texture (UV input):")
a_d = np.sqrt(m_d * m_s)
b_d = np.sqrt(m_s * m_b)
print(f"  a_d = √(m_d · m_s) = {a_d:.2f} MeV")
print(f"  b_d = √(m_s · m_b) = {b_d:.2f} MeV")

# Full CKM from Fritzsch (already computed above, repeat for clarity)
print(f"\nCKM from Fritzsch (6-zero, real entries, zero phases):")
print(f"  |V_us| = {abs(V_CKM[0,1]):.4f}   PDG: {V_us_pdg}   "
      f"deviation: {(abs(V_CKM[0,1])-V_us_pdg)/V_us_pdg*100:+.1f}%")
print(f"  |V_cb| = {abs(V_CKM[1,2]):.4f}   PDG: {V_cb_pdg}   "
      f"deviation: {(abs(V_CKM[1,2])-V_cb_pdg)/V_cb_pdg*100:+.1f}%")
print(f"  |V_ub| = {abs(V_CKM[0,2]):.6f} PDG: {V_ub_pdg:.6f} "
      f"deviation: {(abs(V_CKM[0,2])-V_ub_pdg)/V_ub_pdg*100:+.1f}%")

# With optimal phase
best_Vcb_dev = float('inf')
best_phase = 0
for phase in np.linspace(0, 2*np.pi, 1000):
    M_U_ph = fritzsch_texture(m_u, m_c, m_t, phase_a=0, phase_b=phase)
    M_D_ph = fritzsch_texture(m_d, m_s, m_b, phase_a=0, phase_b=0)
    _, VU = diagonalize_hermitian(M_U_ph)
    _, VD = diagonalize_hermitian(M_D_ph)
    V = VU.conj().T @ VD
    dev = abs(abs(V[1,2]) - V_cb_pdg)
    if dev < best_Vcb_dev:
        best_Vcb_dev = dev
        best_phase = phase
        best_V = V.copy()

print(f"\n  Optimizing phase for |V_cb|: best phase = {np.degrees(best_phase):.1f}°")
print(f"  |V_us| = {abs(best_V[0,1]):.4f}   PDG: {V_us_pdg}")
print(f"  |V_cb| = {abs(best_V[1,2]):.4f}   PDG: {V_cb_pdg}")
print(f"  |V_ub| = {abs(best_V[0,2]):.6f} PDG: {V_ub_pdg:.6f}")

# Also check with both phases free
print(f"\n  Scanning both phases for best overall fit...")
best_chi2 = float('inf')
for pa in np.linspace(0, 2*np.pi, 200):
    for pb in np.linspace(0, 2*np.pi, 200):
        M_U_ph = fritzsch_texture(m_u, m_c, m_t, phase_a=0, phase_b=pb)
        M_D_ph = fritzsch_texture(m_d, m_s, m_b, phase_a=0, phase_b=pa)
        _, VU = diagonalize_hermitian(M_U_ph)
        _, VD = diagonalize_hermitian(M_D_ph)
        V = VU.conj().T @ VD
        chi2 = ((abs(V[0,1]) - V_us_pdg)/0.0008)**2 + \
               ((abs(V[1,2]) - V_cb_pdg)/0.0008)**2 + \
               ((abs(V[0,2]) - V_ub_pdg)/0.00036)**2
        if chi2 < best_chi2:
            best_chi2 = chi2
            best_pa, best_pb = pa, pb
            best_V2 = V.copy()

print(f"  Best fit phases: φ_d = {np.degrees(best_pa):.1f}°, φ_u = {np.degrees(best_pb):.1f}°")
print(f"  |V_us| = {abs(best_V2[0,1]):.4f}   PDG: {V_us_pdg}")
print(f"  |V_cb| = {abs(best_V2[1,2]):.4f}   PDG: {V_cb_pdg}")
print(f"  |V_ub| = {abs(best_V2[0,2]):.6f} PDG: {V_ub_pdg:.6f}")
print(f"  χ² = {best_chi2:.1f}")

# ============================================================
# GST (Gatto-Sartori-Tonin) relation
# ============================================================
print("\n" + "=" * 80)
print("GST RELATION: sin θ_C = √(m_d/m_s)")
print("=" * 80)

V_us_GST = np.sqrt(m_d / m_s)
print(f"  √(m_d/m_s) = {V_us_GST:.4f}")
print(f"  PDG |V_us| = {V_us_pdg}")
print(f"  Deviation: {(V_us_GST - V_us_pdg)/V_us_pdg * 100:+.2f}%")
print(f"  σ = {(V_us_GST - V_us_pdg)/0.0008:.1f}σ")

# The GST relation IS the Fritzsch prediction for V_us when m_u/m_c corrections
# are small: |V_us| ≈ |√(m_d/m_s) - e^{iφ}√(m_u/m_c)|
# Without phase: |V_us| ≈ |√(m_d/m_s) - √(m_u/m_c)| or √(m_d/m_s) + √(m_u/m_c)
correction = np.sqrt(m_u / m_c)
print(f"\n  Fritzsch correction: √(m_u/m_c) = {correction:.4f}")
print(f"  |V_us| ≈ √(m_d/m_s) - √(m_u/m_c) = {V_us_GST - correction:.4f}")
print(f"  |V_us| ≈ √(m_d/m_s) + √(m_u/m_c) = {V_us_GST + correction:.4f}")

# ============================================================
# THE HONEST BOTTOM LINE
# ============================================================
print("\n" + "=" * 80)
print("THE BOTTOM LINE")
print("=" * 80)

print("""
QUESTION: Does the Lagrangian structure (ISS seesaw + Seiberg constraint + det M)
          provide CKM mixing beyond what's put in by hand?

ANSWER: NO.

The Lagrangian with diagonal Yukawas gives CKM = I (identity).
All loop corrections respect the Z₂ off-diagonal parity and cannot generate
off-diagonal meson VEVs. This is exact to all orders.

CKM mixing requires UV INPUT:
  - Off-diagonal entries in the quark mass matrix m̂ within each charge sector
  - These entries mix quarks of the same charge (u↔c, d↔s, etc.)
  - The seesaw TRANSMITS this UV texture faithfully (preserving mixing angles)
  - But it does NOT generate it

The ISS/seesaw framework is COMPATIBLE with CKM mixing (it can transmit a
Fritzsch texture or any other UV input) but it does not PREDICT it.

The free parameters of the CKM are:
  - 3 mixing angles (θ₁₂, θ₂₃, θ₁₃)
  - 1 CP phase (δ)
  All four are UV inputs.

STRUCTURAL INSIGHT: The reason CKM cannot come from the IR is that the meson
mass matrix is block-diagonal (Round 15B), and this block structure is protected
by the Z₂ parity of the off-diagonal mesons. No polynomial Kähler correction can
break it. Only UV physics (the mass matrix m̂) can introduce flavor mixing.
""")

# ============================================================
# WHAT WOULD NEED TO CHANGE FOR CKM TO BE PREDICTED?
# ============================================================
print("=" * 80)
print("WHAT WOULD NEED TO CHANGE FOR CKM TO EMERGE?")
print("=" * 80)

print("""
For CKM mixing to be a PREDICTION rather than an input, one would need:

1. A MECHANISM to generate off-diagonal mass terms m̂_{ij} (i ≠ j, same charge)
   from more fundamental parameters. Examples:
   - Higher-dimensional operators in a UV completion
   - Radiative generation from a GUT (e.g., SU(5) relates up and down sectors)
   - String theory compactification fixing the texture

2. ALTERNATIVELY: a modification of the Lagrangian that breaks the Z₂ parity.
   For example:
   - A term like ε M^u_d M^d_s M^s_u (cubic in off-diagonal mesons)
     → This would break Z₂ at tree level
     → But it would also break charge conservation (M^u_d has charge +1)
     → So this is forbidden by gauge symmetry

3. The ONLY remaining possibility is non-perturbative effects that break the
   Z₂ spontaneously. For example:
   - Instanton effects that generate ⟨M^u_d⟩ ≠ 0
   - But SU(3)_c instantons preserve charge → cannot do this
   - SU(5)_flavor instantons could, but SU(5)_flavor is a global symmetry
     and its instantons are suppressed (no gauge coupling)

CONCLUSION: Within the current framework, CKM is necessarily a UV input.
The seesaw transmits it from high energy to low energy, preserving the texture.
The Fritzsch ansatz (or a refinement of it) provides the simplest UV texture
consistent with observations.

The parameter budget for the full model is:
  - IR parameters: Λ, ⟨X⟩/μ = √3 (fixed by Koide), v₀-doubling coefficient = 3
  - UV inputs: m̂_up texture (3 real params + 1 phase), m̂_down texture (3 real + 1 phase)
  - Minus constraints: det m̂_up = m_u m_c m_t, det m̂_down = m_d m_s m_b (eigenvalues fixed)
  - Net CKM parameters: 3 angles + 1 phase (standard count)
""")
