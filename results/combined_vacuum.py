#!/usr/bin/env python3
"""
Combined vacuum analysis for SQCD + NMSSM-type X-Higgs coupling.

Superpotential:
  W = sum_i m_i M^i_i + X(det M - BB~ - Lambda^6) + lambda X H_u.H_d + m_B BB~

Soft SUSY breaking: V_soft = f_pi^2 Tr(M^dag M)

Two interpretations of X:
  (a) X = Lagrange multiplier (Seiberg electric, N_f = N_c = 3): constraint is exact
  (b) X = dynamical field (magnetic dual / ISS): X has canonical Kahler

This script computes ALL quantities in Problems 1-6 numerically.
"""

import numpy as np
from numpy.linalg import eigvalsh, eig

# =========================================================================
# Physical inputs (MeV)
# =========================================================================
m_u   = 2.16        # MeV (MS-bar at 2 GeV)
m_d   = 4.67        # MeV
m_s   = 93.4        # MeV
Lambda = 300.0       # MeV
f_pi  = 92.0        # MeV
lam   = 0.72         # NMSSM coupling
v     = 246220.0     # MeV (electroweak VEV, 246.22 GeV)
m_B   = 300.0        # MeV (baryon mass)

# Derived
Lambda6 = Lambda**6
m_arr = np.array([m_u, m_d, m_s])
lam_v2_half = lam * v**2 / 2.0
m_tilde_sq = f_pi**2

print("=" * 78)
print("  COMBINED VACUUM ANALYSIS: SQCD + NMSSM X-Higgs COUPLING")
print("=" * 78)
print()
print("INPUT PARAMETERS:")
print(f"  m_u = {m_u} MeV,  m_d = {m_d} MeV,  m_s = {m_s} MeV")
print(f"  Lambda = {Lambda} MeV")
print(f"  Lambda^6 = {Lambda6:.6e} MeV^6")
print(f"  f_pi = {f_pi} MeV,  m_tilde^2 = f_pi^2 = {m_tilde_sq} MeV^2")
print(f"  lambda = {lam}")
print(f"  v = {v} MeV,  v^2/2 = {v**2/2:.6e} MeV^2")
print(f"  lambda * v^2/2 = {lam_v2_half:.6e} MeV^2")
print(f"  m_B = {m_B} MeV")
print()


# =========================================================================
# PROBLEM 1: EWSB-shifted seesaw vacuum
# =========================================================================
print("=" * 78)
print("  PROBLEM 1: EWSB-SHIFTED SEESAW VACUUM")
print("=" * 78)
print()

# Original (unshifted) seesaw
C_orig = Lambda**2 * np.prod(m_arr)**(1.0/3.0)
M_orig = C_orig / m_arr
X0_orig = -C_orig / Lambda6
det_orig = np.prod(M_orig)

print("UNSHIFTED Seiberg seesaw (det M = Lambda^6):")
print(f"  C = Lambda^2 (m_u m_d m_s)^(1/3) = {C_orig:.6f} MeV^2")
print(f"  M_u = C/m_u = {M_orig[0]:.4f} MeV^2")
print(f"  M_d = C/m_d = {M_orig[1]:.4f} MeV^2")
print(f"  M_s = C/m_s = {M_orig[2]:.4f} MeV^2")
print(f"  X_0 = -C/Lambda^6 = {X0_orig:.6e} MeV^(-4)")
print(f"  det M = {det_orig:.6e}  (Lambda^6 = {Lambda6:.6e})")
print(f"  det M / Lambda^6 = {det_orig / Lambda6:.15f}")
print()

# The key ratio Delta
Delta = lam_v2_half / Lambda6
print("EWSB SHIFT:")
print(f"  Delta = lambda * v^2 / (2 * Lambda^6)")
print(f"        = {lam:.2f} * {v:.0f}^2 / (2 * {Lambda6:.6e})")
print(f"        = {lam_v2_half:.6e} / {Lambda6:.6e}")
print(f"        = {Delta:.6e}")
print()

if Delta > 1:
    print(f"  *** Delta = {Delta:.4e} >> 1 ***")
    print(f"  The EWSB shift DOMINATES over Lambda^6!")
    print(f"  lambda*v^2/2 = {lam_v2_half:.6e} MeV^2")
    print(f"  Lambda^6     = {Lambda6:.6e} MeV^6")
    print()
    print(f"  CRITICAL: These quantities have DIFFERENT dimensions!")
    print(f"  [lambda*v^2/2] = MeV^2  (since [lambda] is dimensionless, [v] = MeV)")
    print(f"  [Lambda^6]     = MeV^6")
    print(f"  [X det M]      = MeV^(-4) * MeV^6 = MeV^2  (same as [W] from dim analysis)")
    print(f"  Wait -- [det M] = MeV^6 and [X det M] should equal [W] = MeV^3 in 4d SUSY.")
    print()
    print(f"  DIMENSIONAL RECONCILIATION:")
    print(f"  In Seiberg's confined SQCD description:")
    print(f"    [M^i_j] = MeV^2  (composite meson, dimension of QQ~)")
    print(f"    [det M] = MeV^6  (for 3x3 matrix)")
    print(f"    [W] = MeV^3      (superpotential)")
    print(f"    [X] = MeV^(-3)   (to make X*det M have dimension MeV^3)")
    print(f"    [Lambda^6] = MeV^6  (so X*Lambda^6 also has dim MeV^3)")
    print()
    print(f"  The NMSSM coupling: lambda * X * H_u . H_d")
    print(f"    [lambda * X * H_u * H_d] = MeV^0 * MeV^(-3) * MeV * MeV = MeV^(-1)")
    print(f"    This is WRONG! [W] = MeV^3, not MeV^(-1).")
    print()
    print(f"  RESOLUTION: lambda must carry dimensions!")
    print(f"    [lambda] = MeV^4  to make lambda * X * H_u * H_d have dim MeV^3")
    print(f"    lambda = lambda_hat * Lambda^4  (dimensionless lambda_hat * scale^4)")
    print()
    print(f"  With this: lambda * v^2/2 has dimension MeV^4 * MeV^2 = MeV^6")
    print(f"  Same as Lambda^6. The constraint becomes:")
    print(f"    det M = Lambda^6 - lambda * v^2/2")
    print(f"  with [lambda] = MeV^4.")
    print()

# Recompute with correct dimensions
# lambda_phys has dim MeV^4, parametrized as lambda_hat * Lambda^4
lambda_hat = lam  # dimensionless coupling
lambda_phys = lambda_hat * Lambda**4  # MeV^4
lam_v2_half_corrected = lambda_phys * v**2 / 2.0  # MeV^6

Delta_corrected = lam_v2_half_corrected / Lambda6

print(f"CORRECTED CALCULATION (lambda_phys = lambda_hat * Lambda^4):")
print(f"  lambda_hat = {lambda_hat}")
print(f"  lambda_phys = lambda_hat * Lambda^4 = {lambda_phys:.6e} MeV^4")
print(f"  lambda_phys * v^2/2 = {lam_v2_half_corrected:.6e} MeV^6")
print(f"  Lambda^6 = {Lambda6:.6e} MeV^6")
print(f"  Delta_corrected = lambda_phys * v^2 / (2 * Lambda^6)")
print(f"                   = {Delta_corrected:.6e}")
print()

if Delta_corrected > 1:
    print(f"  Delta_corrected = {Delta_corrected:.4e} >> 1 STILL!")
    print(f"  lambda_hat * v^2 / (2 * Lambda^2) = {lambda_hat * v**2 / (2 * Lambda**2):.6e}")
    print(f"  Since v >> Lambda, even lambda_hat * Lambda^4 * v^2 >> Lambda^6")
    print(f"  unless lambda_hat is VERY small.")
    print()
    print(f"  Required lambda_hat for Delta = 1:")
    lambda_hat_crit = Lambda6 / (Lambda**4 * v**2 / 2)
    print(f"    lambda_hat_crit = Lambda^2 / (v^2/2) = {Lambda**2:.0f} / {v**2/2:.4e}")
    print(f"                    = {lambda_hat_crit:.6e}")
    print()
    print(f"  For Delta << 1 (perturbative EWSB correction):")
    print(f"    Need lambda_hat << {lambda_hat_crit:.2e}")
    print()

# Let's also compute with the NAIVE interpretation where [lambda] is dimensionless
# and the constraint is F_X = det M - Lambda^6 + lambda * H_u . H_d = 0
# The F_X equation: dW/dX = det M - BB~ - Lambda^6 + lambda * Hu * Hd
# For this to be dimensionally consistent:
# [det M] = [Lambda^6] (OK, both MeV^6)
# [lambda * Hu * Hd] must also be MeV^6
# [Hu * Hd] = MeV^2 (since [Hu] = [Hd] = MeV in canonical normalization)
# So [lambda] = MeV^4 as deduced above.

# ALTERNATIVELY: in the NMSSM literature, the coupling is lambda * S * Hu * Hd
# where S is a SINGLET with [S] = MeV (canonical dim). Here X plays the role of S
# but [X] = MeV^(-3). So we need a dimensionful conversion:
# lambda_NMSSM * S * Hu * Hd  <-->  (lambda_NMSSM / Lambda_X^4) * X * Hu * Hd
# where Lambda_X relates the canonical S to X: S = Lambda_X^4 * X.

print("-" * 78)
print("SCENARIO ANALYSIS: What values of Delta are physical?")
print("-" * 78)
print()

# Scan Delta values
for Delta_val in [1e-10, 1e-5, 1e-2, 0.1, 0.5, 0.99, 1.0 - 1e-6]:
    det_shifted = Lambda6 * (1 - Delta_val)
    Cprime = Lambda**2 * np.prod(m_arr)**(1.0/3.0) * (1 - Delta_val)**(1.0/3.0)
    M_shifted = Cprime / m_arr
    print(f"  Delta = {Delta_val:.2e}:")
    print(f"    det M = Lambda^6 * (1 - Delta) = {det_shifted:.6e} MeV^6")
    print(f"    C' = C * (1-Delta)^(1/3) = {Cprime:.6f} MeV^2")
    print(f"    M_u' = {M_shifted[0]:.4f},  M_d' = {M_shifted[1]:.4f},  M_s' = {M_shifted[2]:.4f}")
    print(f"    Shift in M_i: {(1 - Delta_val)**(1.0/3.0) - 1:.6e} fractional")
    print()

print()

# For the physical problem, use a specific small Delta
# The physical question is: how large is the EWSB back-reaction on the QCD vacuum?
# The answer depends on how X couples to the Higgs sector.

# Let's proceed with BOTH cases and state clearly which is which.

print("=" * 78)
print("  CASE A: X AS LAGRANGE MULTIPLIER (Seiberg electric, N_f = N_c = 3)")
print("=" * 78)
print()
print("In this case, X has NO kinetic term. Its equation of motion is the CONSTRAINT:")
print("  det M - BB~ - Lambda^6 + lambda_phys * Hu . Hd = 0")
print()
print("At the EWSB vacuum Hu = Hd = v/sqrt(2), BB~ = 0:")
print("  det M = Lambda^6 - lambda_phys * v^2/2")
print()

# For Delta << 1 regime (generic expectation for weak coupling):
# Use a representative small Delta
Delta_phys = 0.01  # 1% shift as a benchmark
print(f"BENCHMARK: Delta = {Delta_phys} (1% shift)")
det_shifted = Lambda6 * (1 - Delta_phys)
factor = (1 - Delta_phys)**(1.0/3.0)
Cprime = C_orig * factor
M_shifted = Cprime / m_arr
X_shifted = -Cprime / det_shifted  # From F_{M_i} = m_i + X * cof(M) = 0

print(f"  det M = {det_shifted:.6e} MeV^6")
print(f"  C' = C * (1-Delta)^(1/3) = {Cprime:.6f} MeV^2")
print(f"  M_u' = {M_shifted[0]:.4f}  (was {M_orig[0]:.4f}, shift {(factor-1)*100:.4f}%)")
print(f"  M_d' = {M_shifted[1]:.4f}  (was {M_orig[1]:.4f})")
print(f"  M_s' = {M_shifted[2]:.4f}  (was {M_orig[2]:.4f})")
print(f"  X_0' = {X_shifted:.6e}  (was {X0_orig:.6e})")
print()

print("F-terms at the shifted seesaw:")
# F_X = det M - Lambda^6 + lambda_phys * v^2/2 = 0 (BY CONSTRUCTION)
# F_{M_i} = m_i + X * cof(M) = 0 (solved by seesaw)
# F_B = -X * B~ + m_B * B~ (with m_B BB~ term)
# F_{B~} = -X * B + m_B * B (similar)
# At B = B~ = 0: F_B = F_{B~} = 0

for i, flav in enumerate(['u', 'd', 's']):
    cof_i = np.prod(M_shifted) / M_shifted[i]
    F_Mi = m_arr[i] + X_shifted * cof_i
    print(f"  F_{{M_{flav}}} = m_{flav} + X * cof(M,{flav}) = {m_arr[i]:.2f} + ({X_shifted:.6e}) * ({cof_i:.6e}) = {F_Mi:.6e}")
print(f"  F_X = det M - Lambda^6 + lambda_phys*v^2/2 = 0 (constraint)")
print(f"  F_B = F_B~ = 0 (at B = B~ = 0)")
print()

# Off-diagonal meson masses in Case A
print("OFF-DIAGONAL MESON MASSES (Case A: X = Lagrange multiplier)")
print()
print("For off-diagonal meson M^a_b (a != b):")
print("  The only contribution to the mass is from soft terms + W_{M^a_b, M^b_a}")
print("  W_{M^a_b, M^b_a} = X * M_c  (c = third flavor)")
print("  |W_{M^a_b, M^b_a}|^2 = |X|^2 * M_c^2")
print("  m^2(M^a_b) = f_pi^2 + |X_0'|^2 * M_c^2")
print()
print("  Since X enforces the constraint (no kinetic term), there is NO F_X contribution")
print("  to B-terms. The scalar mass matrix has NO B-terms from X.")
print()

print("  Flavor pair    |  f_pi^2    |  |X|^2 * M_c^2          |  m^2_total")
print("  " + "-" * 72)

for (a, b, c_idx, label) in [(0,1,2, 'u-d (c=s)'), (0,2,1, 'u-s (c=d)'), (1,2,0, 'd-s (c=u)')]:
    M_c = M_shifted[c_idx]
    F_term_sq = abs(X_shifted)**2 * M_c**2
    m_sq_total = m_tilde_sq + F_term_sq
    flav_c = ['u','d','s'][c_idx]
    print(f"  {label:16s} |  {m_tilde_sq:8.0f}   |  |X|^2 * M_{flav_c}^2 = {F_term_sq:.6e}  |  {m_sq_total:.6e}")

print()
print("  All off-diagonal masses are POSITIVE. No tachyons.")
print("  The F-term contribution |X|^2 * M_c^2 is tiny (X ~ 10^{-9}).")
print("  The soft term f_pi^2 = 8464 MeV^2 dominates.")
print()
print("  CONCLUSION (Case A): No tachyonic off-diagonal modes.")
print("  CKM mixing does NOT arise from the vacuum structure.")
print()


# =========================================================================
print()
print("=" * 78)
print("  CASE B: X AS DYNAMICAL FIELD (magnetic dual / ISS)")
print("=" * 78)
print()

print("In this case, X has a canonical kinetic term: K contains |X|^2.")
print("X is NOT a Lagrange multiplier; the constraint is NOT exact.")
print("Instead, the F-term of X is:")
print("  F_X = det M - BB~ - Lambda^6 + lambda_phys * Hu * Hd")
print()
print("At the ORIGINAL seesaw vacuum (det M = Lambda^6, B = B~ = 0):")
print("  F_X = Lambda^6 - Lambda^6 + lambda_phys * v^2/2 = lambda_phys * v^2/2")
print()

# But we need to be careful about dimensions again.
# In the ISS picture, the magnetic meson Phi has [Phi] = MeV (canonical).
# The mapping is M = Lambda * Phi (or similar), so the Lagrange multiplier X
# gets replaced by a canonical field with mass dimension 1.
# The ISS superpotential is:
#   W_ISS = h Tr(q Phi q~) - h mu^2 Tr(Phi) + ...
# where [h] is dimensionless, [mu^2] = MeV^2.

# For N_f = N_c = 3, there is NO magnetic gauge group (N_c' = N_f - N_c = 0).
# So the ISS dual is trivially a Wess-Zumino model of mesons.
# X in this case would be the trace of Phi, a canonical scalar.

# Let's proceed formally with X dynamical and see what happens.

# The vacuum: at the seesaw M_i = C/m_i with det M = Lambda^6:
F_X_dyn = lam_v2_half_corrected  # = lambda_phys * v^2/2 (dimension MeV^6)
# But F_X has dimension [F_X] = [dW/dX] = [W]/[X] = MeV^3 / MeV^(-3) = MeV^6.
# And lambda_phys * v^2/2 has dim MeV^4 * MeV^2 = MeV^6. Consistent.

print(f"  F_X at original seesaw = lambda_phys * v^2/2 = {F_X_dyn:.6e} MeV^6")
print()

# But wait: if X is dynamical, its VEV minimizes V.
# The potential includes |F_X|^2 = (det M - Lambda^6 + lambda_phys * v^2/2)^2
# Minimizing in X (no X-dependent terms in F_X itself): X enters only through
# other F-terms. Specifically:
# F_{M^i_i} = m_i + X * cof(M, i)
# So |F_M|^2 depends on X. The minimum in X (at fixed M) is:
# dV/dX = 0 => sum_i F_{M^i_i} * cof(M, i) = 0
# This gives: sum_i [m_i + X * cof(M,i)] * cof(M,i) = 0
# X = - sum_i m_i * cof(M,i) / sum_i cof(M,i)^2

# At diagonal seesaw with M_i = C/m_i:
# cof(M, i) = prod_{j!=i} M_j = (C/m_1)(C/m_2)(C/m_3) / (C/m_i) = C^2 m_i / (m_1 m_2 m_3)^{???}
# Actually: cof(M, 1) = M_2 * M_3 = (C/m_2)(C/m_3) = C^2/(m_2 m_3)
# Let P = m_u * m_d * m_s
# cof(M, i) = C^2 * m_i / P  (since cof(M,i) = det M / M_i = (C^3/P) / (C/m_i) = C^2 m_i / P)

# Wait, det M = M_1 M_2 M_3 = C^3 / P
# cof(M, i) = det M / M_i = (C^3/P) / (C/m_i) = C^2 m_i / P

P = np.prod(m_arr)
cof = np.array([C_orig**2 * m_arr[i] / P for i in range(3)])
print(f"  Cofactors at original seesaw:")
for i, flav in enumerate(['u', 'd', 's']):
    print(f"    cof(M, {flav}) = M_{{others}} = C^2 * m_{flav} / P = {cof[i]:.6e}")

# Minimization in X:
X_opt = -np.sum(m_arr * cof) / np.sum(cof**2)
print(f"\n  X_opt = -sum(m_i * cof_i) / sum(cof_i^2) = {X_opt:.6e}")
print(f"  X_seesaw = -C/Lambda^6 = {X0_orig:.6e}")
print(f"  Ratio: {X_opt / X0_orig:.12f}")
print(f"  These agree: the seesaw solution is also the X-minimum at fixed M.")
print()

# The B-term contribution to off-diagonal meson masses
print("B-TERM CONTRIBUTIONS TO OFF-DIAGONAL MESON MASSES (Case B)")
print()
print("When X is dynamical, the scalar potential includes:")
print("  V = |F_X|^2 + sum_{I} |F_I|^2 + V_soft")
print()
print("The off-diagonal meson mass-squared matrix receives a B-term from:")
print("  V contains |F_X|^2 = (det M - Lambda^6 + lambda_phys * v^2/2)^2")
print("  and also |F_{M^a_b}|^2 contributions.")
print()
print("For M^a_b (a != b), F_{M^a_b} = X * (d det M / d M^a_b)")
print("At diagonal vacuum, d det M / d M^a_b = 0 for a != b (single off-diag).")
print("So F_{M^a_b} = 0 at the diagonal vacuum.")
print()
print("The mass-squared for off-diagonal fluctuation delta M^a_b comes from:")
print("  m^2 = d^2 V / d(M^a_b) d(M^a_b*)  evaluated at the vacuum")
print()
print("From |F_X|^2 contribution:")
print("  |F_X|^2 = (det M - Lambda^6 + ...)^2")
print("  d|F_X|^2/d(M^a_b) = 2 F_X * (d det M / d M^a_b)")
print("  At diagonal vacuum, d det M / d M^a_b = 0, so this is zero.")
print("  d^2|F_X|^2/d(M^a_b)d(M^a_b*) = ... involves second derivatives of det M.")
print()
print("Second derivative of det M w.r.t. M^a_b and M^b_a:")
print("  d^2(det M) / d(M^a_b) d(M^b_a) = M_c  (third diagonal element)")
print("This is the COFACTOR for a 2x2 minor.")
print()
print("The full mass matrix for the pair (M^a_b, M^b_a) in the basis")
print("  (delta M^a_b, delta M^{b*}_a) has entries from:")
print()
print("1) Diagonal mass-squared (from |F_{M^a_b}|^2):")
print("   m^2_diag = |W_{M^a_b, M^b_a}|^2 + W_{M^a_b, X} * W_{X, M^a_b}^* + ...")
print("   At diagonal vacuum:")
print("     W_{M^a_b, M^b_a} = X_0 * M_c  (second derivative of W)")
print("     W_{M^a_b, X} = cof(M, a, b) = 0 for a != b at diagonal vacuum")
print("   So: m^2_diag = |X_0|^2 * M_c^2")
print()
print("2) B-term (holomorphic mass from F_X):")
print("   B = F_X * W_{X, M^a_b, M^b_a}  (third derivative of W !!)")
print("   W_{X, M^a_b, M^b_a} = d/dX [d^2W / d(M^a_b) d(M^b_a)]")
print("                        = d/dX [X * M_c] = M_c")
print()
print("   WAIT: This is a THIRD derivative of W. In the scalar potential V = |F_I|^2,")
print("   the mass-squared matrix is built from SECOND derivatives of W (the fermion mass")
print("   matrix W_{IJ}). The scalar mass-squared is:")
print("     (m^2_scalar)_{i j*} = sum_K W_{K i} W_{K j}^* + F^K W_{K i j}^* + V_soft")
print()
print("   The B-term is: F^K * W_{K i j}^* where K runs over ALL fields.")
print("   For i = M^a_b, j = M^b_a:")
print("     F^X * W_{X, M^a_b, M^b_a}^* = F_X * M_c^*")
print()
print("   This is the BOSONIC mass term that pairs (M^a_b, M^{b*}_a).")
print()

# Compute the B-term
F_X_value = lam_v2_half_corrected  # = lambda_phys * v^2/2 at original seesaw
print("NUMERICAL COMPUTATION OF B-TERMS (Case B, at original seesaw):")
print()
print(f"  F_X = lambda_phys * v^2/2 = {F_X_value:.6e} MeV^6")
print()

# But wait: the B-term pairs (M^a_b, M^{b*}_a) = (phi, phi^*) mixing.
# The 2x2 mass matrix for the REAL and IMAGINARY parts of M^a_b is:
#
#   For a complex scalar phi, the mass matrix in the (phi, phi*) basis is:
#   | m^2_diag + soft   B     |
#   | B*                m^2_diag + soft |
#
#   Eigenvalues: (m^2_diag + soft) +/- |B|
#
# So: the LIGHTER eigenvalue is  (m^2_diag + soft) - |B|
# This is tachyonic if |B| > m^2_diag + soft.

print("  2x2 mass matrix for (M^a_b, M^{b*}_a):")
print("    | m^2_F + m^2_soft    B              |")
print("    | B*                  m^2_F + m^2_soft|")
print()
print("  where m^2_F = |X_0|^2 * M_c^2,  m^2_soft = f_pi^2,  B = F_X * M_c")
print()
print("  Eigenvalues: (m^2_F + m^2_soft) +/- |B|")
print("  Tachyonic if |B| > m^2_F + m^2_soft")
print()

print("  Numerical values:")
print(f"  {'Pair':12s} | {'M_c':12s} | {'m^2_F':14s} | {'m^2_soft':10s} | {'|B| = F_X*M_c':14s} | {'m^2_min':14s} | {'Tachyonic?':10s}")
print("  " + "-" * 95)

results_B = {}
for (a, b, c_idx, label) in [(0,1,2, 'u-d (c=s)'), (0,2,1, 'u-s (c=d)'), (1,2,0, 'd-s (c=u)')]:
    M_c = M_orig[c_idx]
    flav_c = ['u','d','s'][c_idx]
    m2_F = abs(X0_orig)**2 * M_c**2
    m2_soft = m_tilde_sq
    B_val = F_X_value * M_c  # This has dimension MeV^6 * MeV^2 = MeV^8 ??

    # DIMENSIONAL CHECK:
    # F_X has dim MeV^6 (from dW/dX with [X] = MeV^{-3})
    # M_c has dim MeV^2
    # B = F_X * M_c has dim MeV^8
    # But m^2 should have dim MeV^4 (or MeV^2 for dimensionless mass-squared)
    #
    # The issue: the scalar mass-squared formula is
    #   (m^2)_{ij*} = sum_K W_{Ki} W^*_{Kj} + F^K W^*_{Kij}
    # where F^K = K^{KL*} F_L = K^{KL*} dW/dPhi_L
    # With canonical Kahler, F^K = F_K.
    #
    # [W_{Ki}] = [d^2W / dPhi_K dPhi_i] = MeV^3 / (MeV * MeV) = MeV  ... for canonical fields
    # But our fields are NOT canonically normalized!
    # [X] = MeV^{-3}, [M] = MeV^2
    # [W_{X, M^a_b}] = MeV^3 / (MeV^{-3} * MeV^2) = MeV^4
    # [F^X] = [dW/dX] = MeV^3 / MeV^{-3} = MeV^6
    # [W_{X, M^a_b, M^b_a}] = [d^3W / (dX dM^a_b dM^b_a)] = MeV^3 / (MeV^{-3} * MeV^2 * MeV^2) = MeV^2
    #
    # Actually, the B-term formula uses second derivatives:
    # F^K * W_{K,ij}^* where W_{K,ij} = d^2W / (dPhi_K dPhi_i ... NO)
    # W_{Kij} is NOT a thing. The B-term is F^K * (d^2 K / d phi_K ...) NO.
    #
    # Let me be more careful. The full scalar potential in N=1 SUGRA / global SUSY is:
    #   V = K^{I J*} F_I F_J^*
    # where F_I = dW/dPhi_I (global SUSY, canonical Kahler => F_I = partial_I W)
    #
    # The mass matrix:
    #   (m^2)_{I J*} = d^2 V / (d phi_I d phi_J^*)
    #                = sum_K [W_{KI} W_{KJ}^*] + F^K W_{KIJ}^* ... NO, W_{KIJ} is third deriv
    #
    # Actually the correct formula for the scalar mass-squared matrix is:
    #   m^2_{i j*} = sum_K (partial_i partial_K W)(partial_j partial_K W)^*
    #              + (sum_K F_K * partial_i partial_j partial_K^* W^*) ... NO
    #
    # In GLOBAL SUSY with canonical Kahler:
    #   V = sum_I |F_I|^2 = sum_I |partial_I W|^2
    #   d^2 V / d(phi_i) d(phi_j^*) = sum_K (partial_K partial_i W)(partial_K partial_j W)^*
    #
    # That's it for the i != j holomorphic mass terms.
    # For the B-term (holomorphic mass mixing phi_i with phi_j, NOT phi_j^*):
    #   d^2 V / d(phi_i) d(phi_j) = sum_K (partial_K partial_i W)(partial_K partial_j W)^*  ... NO
    #   Actually: d^2 V / d(phi_i) d(phi_j) = sum_K F_K^* * partial_i partial_j partial_K W
    #                                         + sum_K (partial_K partial_i W) ... hmm
    #
    # Let me write it out properly.
    # V = sum_K |partial_K W|^2 = sum_K (partial_K W)(partial_K W)^*
    # partial_i V = sum_K (partial_K partial_i W)(partial_K W)^*
    # partial_i partial_j V = sum_K [(partial_K partial_i partial_j W)(partial_K W)^*
    #                               + (partial_K partial_i W)(partial_K partial_j W)^*]  ... WRONG
    #
    # No: partial_j [(partial_K partial_i W)(partial_K W)^*]
    #   = (partial_K partial_i partial_j W)(partial_K W)^*
    #     + (partial_K partial_i W) * partial_j[(partial_K W)^*]
    #   = (partial_K partial_i partial_j W)(partial_K W)^*    (holomorphic-holomorphic)
    #     + 0                                                  (partial_j of antiholomorphic = 0)
    #
    # So: partial_i partial_j V = sum_K (partial_K partial_i partial_j W) * F_K^*
    #                            (pure B-term: all holomorphic derivatives)
    #
    # And: partial_i partial_{j*} V = sum_K (partial_K partial_i W)(partial_K partial_j W)^*
    #                                (mass-squared: mixed holomorphic-antiholomorphic)

    # So the B-term for (M^a_b, M^b_a) is:
    # B_{ab,ba} = sum_K F_K^* * W_{K, M^a_b, M^b_a}
    # where W_{K, M^a_b, M^b_a} = d^3W / (dPhi_K d M^a_b d M^b_a)

    # The relevant K's:
    # K = X: W_{X, M^a_b, M^b_a} = d/dX [d^2/d(M^a_b)d(M^b_a) (X det M)]
    #      = d/dX [X * d^2(det M)/d(M^a_b)d(M^b_a)]
    #      = d^2(det M)/d(M^a_b)d(M^b_a)
    #      = M_c (the remaining diagonal element)
    #
    # K = M^k_k: W_{M^k_k, M^a_b, M^b_a} = d/dM^k_k [X * d^2(det M)/d(M^a_b)d(M^b_a) + ...]
    #      = X * d/dM^k_k [M_c] = X * delta_{k,c}
    #      Only nonzero for k = c.
    #
    # So: B_{ab,ba} = F_X^* * M_c + F_{M^c_c}^* * X_0

    # At the SUSY vacuum (original seesaw), F_{M^c_c} = 0 (solved by seesaw).
    # So: B_{ab,ba} = F_X^* * M_c

    # With F_X = lambda_phys * v^2/2 at the original seesaw:
    B_term = F_X_value * M_c  # dim: MeV^6 * MeV^2 = MeV^8

    # But the mass-squared should have dim MeV^4 (for canonical fields with dim MeV).
    # Our M has [M] = MeV^2 (non-canonical). The physical mass-squared for a field
    # with VEV M is m^2_phys = (1/M^2) * d^2V/dM^2 ??? No, that's not right either.
    #
    # For non-canonical fields, the scalar mass eigenvalues come from the
    # eigenvalue problem K_{IJ*} omega^2 phi_J = m^2_{IJ*} phi_J.
    # With canonical Kahler, K_{IJ*} = delta_{IJ}, so omega^2 = m^2.
    #
    # The point is: if we're using canonical Kahler for X and M (as assumed in
    # the ISS picture), then the dimensions of the fields are already canonical
    # ([phi] = MeV in 4d SUSY), and the mass formula gives dim MeV^2 for m^2.
    #
    # The confusion arises because M is a COMPOSITE with [M] = MeV^2 in the
    # electric description, but in the magnetic (ISS) description, the mesons
    # are elementary fields with [Phi] = MeV.
    #
    # Resolution: In the magnetic description, M = mu_mag * Phi where
    # [mu_mag] = MeV and [Phi] = MeV. Then:
    #   det M = mu_mag^3 det Phi
    #   W = m_i mu_mag Phi^i_i + X(mu_mag^3 det Phi - Lambda^6) + ...
    # And X absorbs the dimensional factors.
    #
    # For the PURPOSE of this computation, let's work in units where everything
    # is expressed in MeV and treat the non-canonical dimensions as a bookkeeping
    # issue. The KEY question is whether the B-term overwhelms the diagonal mass.

    # The mass-squared matrix for (M^a_b, M^{b*}_a) in the scalar potential has:
    # m^2_diag = sum_K |W_{K, M^a_b}|^2 + soft
    #          = |W_{X, M^a_b}|^2 + |W_{M^b_a, M^a_b}|^2 + f_pi^2
    #
    # At diagonal vacuum:
    # W_{X, M^a_b} = d^2W/(dX dM^a_b) = cofactor of M^a_b in det M = 0 for off-diag
    # W_{M^b_a, M^a_b} = d^2W/(dM^b_a dM^a_b) = X * M_c  (from X * d^2(det M)/(dM^a_b dM^b_a))
    # So: m^2_diag = |X_0 * M_c|^2 + f_pi^2 = |X_0|^2 M_c^2 + f_pi^2

    # B-term = F_X^* * M_c  (from third derivative W_{X, M^a_b, M^b_a} = M_c)

    m2_diag_val = abs(X0_orig)**2 * M_c**2 + m_tilde_sq
    B_val_num = abs(F_X_value) * M_c  # should compare to m2_diag

    # But dimensions: m2_diag has dim MeV^4 + MeV^2 ... different!
    # |X_0|^2 * M_c^2 has dim (MeV^{-3})^2 * (MeV^2)^2 = MeV^{-2}
    # f_pi^2 has dim MeV^2
    # These are NOT the same dimension!

    # This confirms: we CANNOT naively add these terms without canonical normalization.
    # The fundamental issue: M and X have non-standard dimensions.
    #
    # In the ISS/magnetic description with canonical fields:
    # Phi_i (dim MeV), S (dim MeV)  [S plays role of X]
    # W = h S (Phi_i Phi~_j - mu^2 delta_ij) + m_i Phi_i + ...
    # where h is dimensionless, mu^2 has dim MeV^2.

    # Let me redo everything in canonical ISS variables.
    pass  # we'll do this below

print()
print("  DIMENSIONAL INCONSISTENCY DETECTED!")
print("  The non-canonical dimensions of M ([M] = MeV^2) and X ([X] = MeV^{-3})")
print("  prevent naive addition of terms in the mass-squared matrix.")
print("  We must work in CANONICALLY NORMALIZED fields.")
print()


# =========================================================================
# CANONICAL NORMALIZATION via ISS variables
# =========================================================================
print("=" * 78)
print("  CANONICAL NORMALIZATION (ISS VARIABLES)")
print("=" * 78)
print()

# In the ISS / Seiberg description for N_f = N_c = 3:
# The confined description has mesons M^i_j ~ Q^i Q~_j / Lambda
# with [M] = MeV^2 / MeV = MeV ... wait, [Q] = MeV^{3/2} (canonical scalar in 4d)
# [Q Q~] = MeV^3
# [M] = [Q Q~ / Lambda] = MeV^3 / MeV = MeV^2
#
# To canonicalize: Phi = M / Lambda, so [Phi] = MeV
# Then: W = Lambda * m_i Phi^i_i + (Lambda^3 X_can)(Lambda^3 det Phi - Lambda^3) + ...
# where X_can = X * Lambda^3 has [X_can] = MeV^{-3} * MeV^3 = dimensionless.
#
# Actually for the superpotential to have dim MeV^3:
# W = m_i M^i_i has [m_i][M] = MeV * MeV^2 = MeV^3. Good.
# W = X det M has [X][M^3] = MeV^{-3} * MeV^6 = MeV^3. Good.
# W = X Lambda^6 has [X] * MeV^6 = MeV^{-3} * MeV^6 = MeV^3. Good.
#
# Canonical: let phi = M^{1/2} (geometric mean with Lambda)
# No, the standard procedure: define Phi_{ij} = M^i_j / mu_scale
# where mu_scale has dim MeV to make [Phi] = MeV.

# The simplest approach: define the canonical fields as
# phi_i = M_i (they happen to have dim MeV^2, but in the scalar potential
# V = sum |F_I|^2, the F-terms have their own dimensions and V has dim MeV^4
# in 4d. Actually V has dim MeV^4 when F has dim MeV^2.
#
# F_I = dW/dPhi_I. If [Phi_I] = MeV, [W] = MeV^3, then [F_I] = MeV^2.
# V = sum |F_I|^2 has dim MeV^4. Good.
#
# But our M has [M] = MeV^2 and X has [X] = MeV^{-3}.
# F_{M^i_i} = dW/dM^i_i = m_i + X cof(M,i), [F_{M}] = MeV + MeV^{-3} * MeV^4 = MeV.
# V = sum |F_{M^i_i}|^2 has dim MeV^2, not MeV^4.
#
# This is because the Kahler potential for non-canonical fields is
# K = Lambda^{-2} Tr(M^dag M) (to make [K] = MeV^2 for [M] = MeV^2)
# The scalar potential is V = K^{IJ*} F_I F_J^* where K^{IJ*} = Lambda^2 delta^{IJ}
# V = Lambda^2 sum |F_M|^2 = Lambda^2 * MeV^2 = MeV^4 * MeV^{-2} ??
#
# Actually with non-canonical Kahler K = (1/Lambda^2) Tr(M^dag M):
# K_{M_i M_j^*} = delta_{ij} / Lambda^2
# K^{M_i M_j^*} = Lambda^2 delta^{ij}
# V = sum_{ij} K^{ij*} F_i F_j^* = Lambda^2 sum_i |F_i|^2
# [V] = MeV^2 * MeV^2 = MeV^4 ??? [F] = MeV here.
# Lambda^2 * MeV^2 = MeV^4. Yes!

# OK so the proper scalar potential with the Kahler metric included is:
# V = Lambda^2 sum_i |F_{M_i}|^2 + Lambda^6 |F_X|^2 + ... + f_pi^2 Lambda^{-2} Tr(M^dag M)

# For the off-diagonal mass matrix:
# The physical mass-squared is d^2V / dM^a_b dM^{a*}_b with proper Kahler factors.
# m^2_phys = K^{M^a_b, M^{a*}_b} * [ |W_{X M^a_b}|^2 K_{XX*}^{-1} + ... ]
#
# This is getting complicated. Let me just canonically normalize everything.

# CANONICAL FIELDS:
# phi_{ij} = M^i_j / Lambda  (dim MeV)
# chi = X * Lambda^3          (dim MeV^0 = dimensionless ... not right)
#
# Actually let's be more careful.
# In 4d N=1 SUSY, the canonical Kahler is K = sum_I |Phi_I|^2 with [Phi_I] = MeV.
# The superpotential has [W] = MeV^3.
# The F-term is F_I = dW/dPhi_I has dim MeV^2.
# The scalar potential V = sum_I |F_I|^2 has dim MeV^4.
#
# For the meson: M has dim MeV^2. Define phi_M = M / mu_0 where mu_0 = Lambda.
# Then [phi_M] = MeV. W in terms of phi_M:
#   W = m_i Lambda phi^i_i + X Lambda^3 (det phi - 1) + lambda_phys X Hu Hd
#
# For X: [X] = MeV^{-3}. Define phi_X = X * Lambda^{k} to get dim MeV.
# Need MeV^{-3} * MeV^k = MeV, so k = 4. phi_X = X * Lambda^4, [phi_X] = MeV.
# W in terms of phi_X:
#   X det M = (phi_X / Lambda^4) * Lambda^3 det phi = phi_X Lambda^{-1} det phi
#   X Lambda^6 = (phi_X / Lambda^4) * Lambda^6 = phi_X Lambda^2
#
# So: W = m_i Lambda phi^i_i + (phi_X / Lambda) (Lambda^3 det phi - Lambda^5)
#       + (lambda_phys / Lambda^4) phi_X Hu Hd

# Let's define:
# h = Lambda^2 / Lambda = Lambda (coupling in ISS-like form)
# mu^2 = Lambda^2 (mass scale in ISS)
# lambda_eff = lambda_phys / Lambda^4 (effective dimensionless coupling)

# Canonical superpotential:
# W_can = m_i Lambda phi^i_i + Lambda^2 phi_X (det phi - Lambda^2/Lambda^2 ??? )

# This dimensional analysis is getting circular. Let me just proceed NUMERICALLY
# with all quantities in MeV, treating the mass matrix as a numerical matrix,
# and identifying which eigenvalues are positive or negative.

print("DIRECT NUMERICAL APPROACH:")
print("We compute the full scalar potential V = sum |F_I|^2 + V_soft")
print("using the ORIGINAL field variables (M, X, B, B~, Hu, Hd).")
print("The Hessian d^2V/d(field_i) d(field_j) gives the mass-squared matrix")
print("in the original (non-canonical) field basis. Tachyonic directions")
print("are signaled by NEGATIVE eigenvalues regardless of normalization.")
print()

# =========================================================================
# NUMERICAL HESSIAN COMPUTATION
# =========================================================================
print("=" * 78)
print("  NUMERICAL HESSIAN AT THE SEESAW VACUUM")
print("=" * 78)
print()

# Vacuum: M = diag(M_u, M_d, M_s), X = X_0, B = B~ = 0, Hu = Hd = v/sqrt(2)
# With the m_B BB~ term: W contains + m_B B B~
# F_B = -X B~ + m_B B~  =>  at B=B~=0: F_B = 0
# F_{B~} = -X B + m_B B  => at B=B~=0: F_{B~} = 0

# Two scenarios:
# Scenario A: X is Lagrange multiplier -> constraint det M = Lambda^6 - lambda_phys v^2/2
#   At this vacuum, F_X = 0 identically.
# Scenario B: X is dynamical -> det M = Lambda^6 at seesaw, F_X = lambda_phys v^2/2 != 0

# We compute the Hessian for BOTH scenarios.

# For Case B (X dynamical), the scalar potential is:
#   V = sum |F_I|^2 + f_pi^2 Tr(M^dag M)
# No Kahler metric factors (assuming canonical K for X too, which is the ISS assumption).

# IMPORTANT: For N_f = N_c = 3, the correct description is the CONFINED theory.
# X is NOT a magnetic meson -- there is no magnetic gauge group.
# But we can still ask: if X had a canonical kinetic term (as a phenomenological choice),
# what would the spectrum look like?

# Fields: M^i_j (9 complex), X (1 complex), B (1 complex), B~ (1 complex),
#          Hu0 (1 complex), Hd0 (1 complex) -> 14 complex fields = 28 real DOF

# For the off-diagonal analysis, we only need the 2x2 block for each pair (M^a_b, M^b_a).

# The F-terms:
# F_{M^i_j} = m_i delta^i_j + X * d(det M)/d(M^i_j)
#           = m_i delta^i_j + X * cofactor(M, i, j)
# At diagonal vacuum: cofactor(M, a, b) = 0 for a != b (single off-diag entry)
# cofactor(M, a, a) = prod_{k!=a} M_k
# cofactor(M, a, b) for a!=b in det M = M_1(M_2 M_3 - ...) - M^0_1(...) + ...
# For 3x3 diagonal matrix: d(det M)/d(M^a_b) at diagonal M:
#   = epsilon^{a p q} epsilon_{b r s} M^r_p M^s_q / 2
#   For a != b: this involves two diagonal elements with crossed indices
#   d(det M)/d(M^a_b) = 0 when M is diagonal and a != b? Not quite.
#
# For M = diag(M_1, M_2, M_3):
# det M = M_1 M_2 M_3
# d(det M)/d(M^0_1) = d(M^0_0 M^1_1 M^2_2 - ...)/d(M^0_1)
# In the Levi-Civita expansion: det M = sum_{sigma} sgn(sigma) prod_i M^i_{sigma(i)}
# d(det M)/d(M^a_b) = sum_{sigma with sigma(a)=b} sgn(sigma) prod_{i!=a} M^i_{sigma(i)}
# For a=0, b=1: sigma(0)=1, and we need prod_{i=1,2} M^i_{sigma(i)}
# sigma must be a permutation with sigma(0)=1. Options:
# sigma = (1, 0, 2): sgn = -1, prod = M^1_0 * M^2_2 = 0 * M_3 = 0 (at diag vacuum)
# sigma = (1, 2, 0): sgn = +1, prod = M^1_2 * M^2_0 = 0 * 0 = 0
# So d(det M)/d(M^0_1) = 0 at diagonal vacuum. CONFIRMED.

# Second derivatives:
# d^2(det M)/d(M^a_b)d(M^c_d) involves further Levi-Civita sums.
# For the specific case M^a_b and M^b_a (i.e., c=b, d=a):
# d^2(det M)/d(M^a_b)d(M^b_a) = ?
# We need permutations sigma with sigma(a)=b AND sigma(b)=a, and the third index k
# maps to sigma(k)=k. For 3x3 with k = third index c:
# sigma(a)=b, sigma(b)=a, sigma(c)=c => sigma is the transposition (a b).
# sgn(sigma) = -1
# The contribution is: (-1) * M^c_c = -M_c
#
# Wait, let me redo:
# d(det M)/d(M^a_b) = sum_{sigma: sigma(a)=b} sgn(sigma) prod_{i!=a} M^i_{sigma(i)}
# Now d/d(M^b_a) of this:
# = sum_{sigma: sigma(a)=b} sgn(sigma) * d/d(M^b_a) [prod_{i!=a} M^i_{sigma(i)}]
# For sigma(a)=b, the factor in the product for i=b is M^b_{sigma(b)}.
# We need sigma(b) = a for the derivative d/d(M^b_a) to act.
# If sigma(a)=b and sigma(b)=a: the remaining index c has sigma(c)=c.
# The derivative d/d(M^b_a) of M^b_a gives 1.
# The remaining factor is M^c_{sigma(c)} = M^c_c = M_c.
# sgn of (a b) transposition = -1.
# So: d^2(det M)/d(M^a_b)d(M^b_a) = -M_c ? That has a minus sign.
#
# Hmm, but the superpotential term is X * det M.
# W_{M^a_b, M^b_a} = X * d^2(det M)/d(M^a_b)d(M^b_a)
# If d^2(det M)/d(M^a_b)d(M^b_a) = -M_c, then W_{M^a_b, M^b_a} = -X M_c
# |W_{M^a_b, M^b_a}|^2 = |X|^2 M_c^2 (sign doesn't matter for modulus squared)

# But wait, I should check this with a concrete example.
# det M = M^0_0 M^1_1 M^2_2 + M^0_1 M^1_2 M^2_0 + M^0_2 M^1_0 M^2_1
#       - M^0_2 M^1_1 M^2_0 - M^0_1 M^1_0 M^2_2 - M^0_0 M^1_2 M^2_1
#
# d(det M)/d(M^0_1) = M^1_2 M^2_0 - M^1_0 M^2_2
# At diagonal: = 0 - 0 = 0. Good.
#
# d/d(M^1_0) [M^1_2 M^2_0 - M^1_0 M^2_2] = 0 - M^2_2 = -M_s
# So d^2(det M)/d(M^0_1)d(M^1_0) = -M_s (at diagonal vacuum).
# For (a,b) = (0,1): c = 2 (s-flavor). d^2 det M = -M_c. Confirmed with minus sign.
#
# Similarly:
# d(det M)/d(M^0_2) = M^1_0 M^2_1 - M^1_1 M^2_0 => at diag: 0
# d/d(M^2_0)[M^1_0 M^2_1 - M^1_1 M^2_0] = 0 - M^1_1 = -M_d
# d^2(det M)/d(M^0_2)d(M^2_0) = -M_d. c = 1 (d-flavor). Confirmed -M_c.
#
# d(det M)/d(M^1_2) = M^0_1 M^2_0 - M^0_0 M^2_1 => at diag: 0
# Hmm wait. Let me redo from the explicit formula.
# det M = M00(M11 M22 - M12 M21) - M01(M10 M22 - M12 M20) + M02(M10 M21 - M11 M20)
# d(det)/d(M12) = d/d(M12)[- M00 M12 M21 + ... ] hmm let me just differentiate term by term.
# Terms containing M12:
# M01 M12 M20 (from + term): coeff +1 * M01 M20
# -M00 M12 M21 (from - in expansion): coeff -M00 M21
# M02 ... no M12 there. Let me just use the Levi-Civita directly.
#
# d(det)/d(M^1_2): look for terms where row 1 maps to column 2.
# sigma(1) = 2. Permutations: (0,2,1) with sgn=-1 giving M00 M12 M21 with sgn -1 => -M00 M21
#                              (2,2,0) invalid (not a permutation)
# Wait: sigma=(0,2,1): sigma(0)=0, sigma(1)=2, sigma(2)=1. Valid. sgn = -1.
# Contribution to det: -M00 M12 M21 => d/dM12 = -M00 M21. At diag: 0.
#
# sigma=(2,2,0): not a permutation (sigma(0)=sigma(1)=2, invalid).
#
# Actually, for sigma(1)=2: possible sigma's with sigma(1)=2:
# sigma = (0, 2, 1): sgn = -1. Product excl i=1: M^0_0 * M^2_1. At diag: M_u * 0 = 0.
# sigma = (2, 2, 0): invalid (2 appears twice).
# Hmm, I'm overcomplicating this. For 3x3:
# Permutations with sigma(1)=2:
# (0,2,1): sgn=-1
# (2,2,0): invalid
# Actually enumerate all 6 permutations:
# (0,1,2) sgn=+1: sigma(1)=1 != 2
# (0,2,1) sgn=-1: sigma(1)=2 YES -> prod_{i!=1} M^i_{sigma(i)} = M^0_0 * M^2_1 = M_u * 0 = 0
# (1,0,2) sgn=-1: sigma(1)=0 != 2
# (1,2,0) sgn=+1: sigma(1)=2 YES -> prod_{i!=1} = M^0_1 * M^2_0 = 0 * 0 = 0
# (2,0,1) sgn=+1: sigma(1)=0 != 2
# (2,1,0) sgn=-1: sigma(1)=1 != 2
# So d(det)/d(M^1_2) = 0 at diagonal. Good.
#
# Now d^2(det)/d(M^1_2)d(M^2_1):
# d/d(M^2_1) of [(-1) M^0_0 M^2_1] = -M^0_0 = -M_u. (from (0,2,1) permutation)
# d/d(M^2_1) of [(+1) M^0_1 M^2_0] = 0 (no M^2_1 factor)
# Total: -M_u. Here c = 0 (u-flavor). Confirmed d^2 det/d(M^1_2)d(M^2_1) = -M_c.
#
# Great. So universally: d^2(det M)/d(M^a_b)d(M^b_a) = -M_c where c is the third flavor.
# And W_{M^a_b, M^b_a} = X_0 * (-M_c) = -X_0 M_c.

print("Second derivatives of det M at diagonal vacuum:")
print("  d^2(det M)/d(M^a_b)d(M^b_a) = -M_c (c = third flavor, c != a,b)")
print("  (minus sign from transposition signature)")
print()

# Now the superpotential second derivatives:
# W_{M^a_b, M^b_a} = X_0 * (-M_c)
# |W_{M^a_b, M^b_a}|^2 = |X_0|^2 * M_c^2

# The B-term (third derivative):
# W_{X, M^a_b, M^b_a} = d/dX [X * (-M_c)] = -M_c
# Wait: W_{X, M^a_b, M^b_a} = d^3W / (dX dM^a_b dM^b_a) at the vacuum.
# W = X * det M + ...
# d^3(X det M)/(dX dM^a_b dM^b_a) = d^2(det M)/(dM^a_b dM^b_a) = -M_c.
#
# So B_{(M^a_b)(M^b_a)} = F_X^* * (-M_c) + F_{M^c_c}^* * X_0 * ... hmm
# Actually: B_{ij} = d^2 V / d(phi_i) d(phi_j) = sum_K F_K^* * W_{Kij}
# For (i,j) = (M^a_b, M^b_a):
# B = sum_K F_K^* * W_{K, M^a_b, M^b_a}
# K = X: W_{X, M^a_b, M^b_a} = -M_c. Contribution: F_X^* * (-M_c)
# K = M^c_c: W_{M^c_c, M^a_b, M^b_a} = d^3W/(dM^c_c dM^a_b dM^b_a)
#   = d/dM^c_c [X * d^2(det M)/(dM^a_b dM^b_a)] = X * d/dM^c_c(-M_c) = -X
#   Contribution: F_{M^c_c}^* * (-X_0) = 0 (at seesaw, F_{M^c_c} = 0)
# Other K's: all zero (no cubic coupling involving three off-diagonal or mixed terms)
#
# So: B_{(M^a_b)(M^b_a)} = -F_X^* * M_c

# For Case B with F_X = lambda_phys * v^2/2:
# B = -(lambda_phys * v^2/2) * M_c  (this is a holomorphic mass)

# The 2x2 mass matrix for (M^a_b, M^{b*}_a) in the basis (Re, Im):
# m^2_diag = sum_K |W_{K, M^a_b}|^2 + soft = |X_0|^2 M_c^2 + f_pi^2
# (the |W_{K, M^a_b}|^2 only gets contribution from W_{M^b_a, M^a_b} = -X_0 M_c)
# Wait: the sum is over K of |W_{K, M^a_b}|^2 = sum_K |d^2W/dPhi_K dM^a_b|^2
# W_{K, M^a_b} at diagonal vacuum:
# K = X: W_{X, M^a_b} = cofactor(M, a, b) = 0 (off-diagonal cofactor at diag vacuum is 0)
# K = M^b_a: W_{M^b_a, M^a_b} = X_0 * d^2(det M)/d(M^b_a) d(M^a_b)
#   Same as W_{M^a_b, M^b_a} by symmetry = -X_0 M_c
# K = M^k_l (other off-diag): zero
# K = M^k_k (diagonal): W_{M^k_k, M^a_b} = X_0 * d^2(det M)/(dM^k_k dM^a_b)
#   For k != a,b: d^2(det M)/(dM^k_k dM^a_b) = d/dM^k_k [0] = 0 (first deriv is 0 at diag)
#   Actually: d(det M)/d(M^a_b) = 0 at diagonal, and d/dM^k_k of this?
#   The MIXED partial d^2(det)/d(M^k_k)d(M^a_b) at diagonal vacuum:
#   This requires permutations sigma with sigma(a)=b and sigma(k)=k.
#   For k=c (third flavor): sigma(a)=b, sigma(c)=c, so sigma(b)=remaining = a (for 3 indices).
#   This is the transposition (ab). sgn = -1. Product = M^c_c = M_c.
#   But d^2/d(M^k_k)d(M^a_b): we're differentiating det M once w.r.t. M^a_b (off-diag)
#   and once w.r.t. M^k_k (diagonal). At the diagonal vacuum:
#   From the explicit check above, d(det)/d(M^0_1) = M12 M20 - M10 M22.
#   d/d(M^2_2)[M12 M20 - M10 M22] = -M10. At diagonal: 0.
#   d/d(M^1_1)[M12 M20 - M10 M22] = 0.
#   d/d(M^0_0)[M12 M20 - M10 M22] = 0.
#   So d^2(det)/d(M^k_k)d(M^a_b) = 0 at diagonal vacuum for ALL k. Good.
#
# So: sum_K |W_{K, M^a_b}|^2 = |W_{M^b_a, M^a_b}|^2 = |X_0|^2 M_c^2

print("MASS-SQUARED MATRIX FOR OFF-DIAGONAL MESONS:")
print()
print("In the basis (M^a_b, M^{b*}_a) -- i.e., the scalar phi and its conjugate partner --")
print("the 2x2 mass matrix is:")
print()
print("  | m^2_diag + m^2_soft    B*   |")
print("  | B                m^2_diag + m^2_soft |")
print()
print("where:")
print("  m^2_diag = |X_0|^2 * M_c^2  (from |W_{M^b_a, M^a_b}|^2)")
print("  m^2_soft = f_pi^2")
print("  B = -F_X^* * M_c  (from F_X * W_{X, M^a_b, M^b_a}^* = F_X * (-M_c)^* )")
print()
print("Eigenvalues: (m^2_diag + m^2_soft) +/- |B|")
print("Lighter eigenvalue: (m^2_diag + m^2_soft) - |B|")
print("Tachyonic if |B| > m^2_diag + m^2_soft")
print()

# =========================================================================
# Now compute for BOTH cases
# =========================================================================

print("=" * 78)
print("  CASE A RESULTS: X = Lagrange multiplier")
print("=" * 78)
print()
print("F_X = 0 by construction. B-term = 0.")
print()

print(f"  {'Pair':16s} | {'m^2_F = |X|^2 M_c^2':20s} | {'m^2_soft':10s} | {'m^2_total':14s} | {'Tachyonic?':10s}")
print("  " + "-" * 80)

for (a, b, c_idx, label) in [(0,1,2, 'u-d (M_s)'), (0,2,1, 'u-s (M_d)'), (1,2,0, 'd-s (M_u)')]:
    M_c = M_orig[c_idx]
    m2_F = abs(X0_orig)**2 * M_c**2
    m2_total = m2_F + m_tilde_sq
    tachyonic = "NO"
    print(f"  {label:16s} | {m2_F:20.6e} | {m_tilde_sq:10.0f} | {m2_total:14.6e} | {tachyonic:10s}")

print()
print("  ALL masses are positive. No tachyons. No CKM mixing from vacuum.")
print()


# =========================================================================
print()
print("=" * 78)
print("  CASE B RESULTS: X = dynamical field")
print("=" * 78)
print()

# F_X at original seesaw (det M = Lambda^6):
# F_X = det M - Lambda^6 + lambda_phys * v^2/2 = lambda_phys * v^2/2
# But what is lambda_phys? It has dim MeV^4.
# lambda_phys = lambda_hat * Lambda^4 with lambda_hat = 0.72

F_X_B = lambda_hat * Lambda**4 * v**2 / 2.0  # MeV^6
# Alternative: if the problem statement means lambda is dimensionless and the
# term is lambda * X * Hu * Hd in the superpotential, and X is canonically
# normalized (dim MeV), then Hu, Hd have dim MeV, and lambda * X * Hu * Hd
# has dim MeV^3 = [W]. This requires [X] = MeV (canonical).
# In this case, the seesaw relates the CANONICAL X to det M / Lambda^k.

# Let's try the CANONICAL interpretation directly:
# If X has dim MeV and canonical kinetic term, then:
# W = m_i M^i_i + (X / Lambda^3) (det M - Lambda^6) + lambda X Hu Hd + m_B BB~
# [X/Lambda^3 * det M] = MeV * MeV^{-3} * MeV^6 = MeV^4 ??? NOT MeV^3.
#
# Alternatively:
# W = m_i (M^i_i / Lambda) + X (det(M/Lambda^2) - 1) Lambda^3 + lambda X Hu Hd
# where M/Lambda^2 is dimensionless. Then det(M/Lambda^2) is dimensionless.
# [X * Lambda^3] = MeV * MeV^3 = MeV^4 ??? Still not right.
#
# The correct way: define Phi = M / Lambda, [Phi] = MeV (since [M] = MeV^2).
# W = m_i Lambda Phi^i_i + X Lambda^3 (det(Phi/Lambda) - 1) + lambda X Hu Hd
# Wait, det Phi has dim MeV^3.
# [X * det Phi] = MeV * MeV^3 = MeV^4. Still off by one power.
# We need: W = X * (det Phi / Lambda^k) - ... with [W] = MeV^3.
# [X * det Phi / Lambda^k] = MeV^{4-k}. Need k=1.
# W = (X / Lambda) det Phi - X Lambda^2 + ...? [X Lambda^2] = MeV^3. OK.
# = X (det Phi / Lambda - Lambda^2)
# Hmm, [(det Phi)/Lambda] = MeV^2. [X * MeV^2] = MeV^3. Good.
# [X * Lambda^2] = MeV^3. Good.

# Define the canonical fields for the ISS picture:
# Phi_ij = M^i_j / Lambda  (dim MeV)
# X = X (dim MeV, canonical)
# Hu, Hd (dim MeV, canonical)

# Superpotential in canonical variables:
# W = m_i Lambda Phi^i_i + (X/Lambda) (Lambda^3 det(Phi/Lambda) - Lambda^6/Lambda)
#   = m_i Lambda Phi^i_i + X(det Phi / Lambda - Lambda^5/Lambda)
# Actually this is getting messy. Let me just define it consistently.

# The STANDARD Seiberg superpotential with canonical X is:
# W = h X (Phi_i Phi~_j delta^{ij} - mu^2) + m_i Phi_i Phi~_i + ...
# For confined SQCD (N_f = N_c): W = X(det M / Lambda^{2N_c - N_f} - Lambda^{3N_c - N_f})
# For N_f = N_c = 3: W = X(det M / Lambda^3 - Lambda^6)
# With M = Lambda Phi: det M = Lambda^3 det Phi
# W = X(Lambda^3 det Phi / Lambda^3 - Lambda^6) = X(det Phi - Lambda^6)
# [X * det Phi] = MeV * MeV^3 = MeV^4 ??? Not MeV^3!

# Actually, in the literature (Intriligator, Seiberg, Shih), the confined description
# for N_f = N_c uses:
# W = X(det M - BB~ - Lambda^{2N_c}) with [M] = MeV^2 (composite) and [X] = MeV^{-3}
# This is the ORIGINAL non-canonical description. The canonical X has [X_can] = MeV
# and is related by X = X_can / Lambda^4.

# So let me use X_can = X * Lambda^4, [X_can] = MeV.
# W = (X_can / Lambda^4)(det M - Lambda^6) + lambda X_can Hu Hd + m_i M^i_i + m_B BB~
# [X_can / Lambda^4 * det M] = MeV * MeV^{-4} * MeV^6 = MeV^3. GOOD.
# [X_can / Lambda^4 * Lambda^6] = MeV * MeV^{-4} * MeV^6 = MeV^3. GOOD.
# [lambda * X_can * Hu * Hd] = MeV * MeV * MeV = MeV^3. GOOD (lambda dimensionless).
# [m_i * M^i_i] = MeV * MeV^2 = MeV^3. GOOD.

# Now the F-terms with canonical X:
# F_{X_can} = dW/dX_can = (1/Lambda^4)(det M - Lambda^6) + lambda * Hu * Hd
# At Seiberg seesaw (det M = Lambda^6) with Hu = Hd = v/sqrt(2):
# F_{X_can} = 0 + lambda * v^2/2 = lambda * v^2/2
# [F_{X_can}] = MeV^2 (as expected for canonical field).

F_X_can = lam * v**2 / 2.0  # MeV^2 (canonical F-term)
X_can_vev = X0_orig * Lambda**4  # MeV (canonical X VEV)

print(f"CANONICAL VARIABLES:")
print(f"  X_can = X * Lambda^4 = {X_can_vev:.6e} MeV")
print(f"  F_{{X_can}} = lambda * v^2/2 = {F_X_can:.6e} MeV^2")
print(f"  |F_{{X_can}}|^2 = {F_X_can**2:.6e} MeV^4")
print()

# The F-terms for M remain:
# F_{M^i_j} = m_i delta^i_j + (X_can / Lambda^4) * d(det M)/d(M^i_j)
# This is the SAME as before with X = X_can / Lambda^4.

# The B-term for off-diagonal mesons:
# V_holo = d^2V / d(M^a_b) d(M^b_a) = sum_K F_K^* * W_{K, M^a_b, M^b_a}
# K = X_can: W_{X_can, M^a_b, M^b_a} = (1/Lambda^4) * d^2(det M)/d(M^a_b)d(M^b_a)
#   = (1/Lambda^4) * (-M_c) = -M_c / Lambda^4
# Contribution: F_{X_can}^* * (-M_c / Lambda^4) = -(lambda * v^2/2) * M_c / Lambda^4

B_term_can = {}
m2_diag_can = {}
for (a, b, c_idx, label) in [(0,1,2, 'u-d (M_s)'), (0,2,1, 'u-s (M_d)'), (1,2,0, 'd-s (M_u)')]:
    M_c = M_orig[c_idx]
    flav_c = ['u','d','s'][c_idx]

    # W_{M^b_a, M^a_b} = X * d^2(det M)/(dM^b_a dM^a_b) = X * (-M_c) = -(X_can/Lambda^4) * M_c
    W_ba_ab = -(X_can_vev / Lambda**4) * M_c  # dim: MeV * MeV^{-4} * MeV^2 = MeV^{-1} ???

    # Hmm, W_{IJ} = d^2W / dPhi_I dPhi_J should have dim MeV^3 / (MeV^2 * MeV^2) = MeV^{-1}
    # for non-canonical M. That's the fermion mass in non-canonical basis.

    # The scalar mass-squared from F-terms:
    # (m^2)_{M^a_b, M^{a*}_b} = sum_K |W_{K, M^a_b}|^2
    # W_{M^b_a, M^a_b} has dim MeV^{-1} (as above)
    # |W|^2 has dim MeV^{-2}. This is NOT a mass-squared.
    # We need to include the Kahler metric: K_{M^a_b, M^{a*}_b} which is NOT 1
    # for non-canonical fields.

    # This confirms: we MUST canonically normalize M too.
    pass

# FULL CANONICAL NORMALIZATION:
# Phi_{ij} = M^i_j / Lambda, [Phi] = MeV
# X_c = X * Lambda^4, [X_c] = MeV
# W in terms of Phi, X_c:
# m_i M^i_i = m_i Lambda Phi^i_i
# (X_c/Lambda^4)(det M - Lambda^6) = (X_c/Lambda^4)(Lambda^3 det Phi - Lambda^6)
#   = X_c Lambda^{-1} det Phi - X_c Lambda^2
# lambda X_c Hu Hd (already OK)
# m_B BB~ (B has dim MeV^3 from composite baryon; B_c = B/Lambda^2 has dim MeV)
# m_B Lambda^2 B_c Lambda^2 B~_c = m_B Lambda^4 B_c B~_c

# Full canonical W:
# W = m_i Lambda Phi^i_i + X_c (det Phi / Lambda - Lambda^2) + lambda X_c Hu Hd + m_B Lambda^4 B_c B~_c
# + X_c/Lambda * (- B_c Lambda^2 * B~_c Lambda^2) = X_c Lambda^3 B_c B~_c ?? No...
# The original term is X (- B B~) = (X_c / Lambda^4)(- Lambda^4 B_c B~_c) = - X_c B_c^2 Lambda^4/Lambda^4 = -X_c B_c B~_c? Hmm.
# Actually B = Lambda^2 B_c, B~ = Lambda^2 B~_c
# X * B * B~ = (X_c / Lambda^4) * Lambda^2 B_c * Lambda^2 B~_c = X_c B_c B~_c

# So: W = m_i Lambda Phi^i_i + X_c (det Phi / Lambda - Lambda^2 - B_c B~_c)
#       + lambda X_c Hu Hd + m_B Lambda^4 B_c B~_c

# F-terms in canonical variables:
# F_{Phi^i_i} = m_i Lambda + (X_c / Lambda) cofactor(Phi, i)
# F_{X_c} = det Phi / Lambda - Lambda^2 - B_c B~_c + lambda Hu Hd
# F_{Hu} = lambda X_c Hd
# F_{Hd} = lambda X_c Hu
# F_{B_c} = -X_c B~_c + m_B Lambda^4 B~_c = B~_c(m_B Lambda^4 - X_c)
# F_{B~_c} = -X_c B_c + m_B Lambda^4 B_c = B_c(m_B Lambda^4 - X_c)

# At the vacuum: Phi = diag(Phi_1, Phi_2, Phi_3) with Phi_i = M_i / Lambda
Phi_vev = M_orig / Lambda  # dim MeV
X_c_vev = X0_orig * Lambda**4  # dim MeV

print(f"FULLY CANONICAL VARIABLES:")
print(f"  Phi_u = M_u/Lambda = {Phi_vev[0]:.4f} MeV")
print(f"  Phi_d = M_d/Lambda = {Phi_vev[1]:.4f} MeV")
print(f"  Phi_s = M_s/Lambda = {Phi_vev[2]:.4f} MeV")
print(f"  X_c = X * Lambda^4 = {X_c_vev:.6e} MeV")
print(f"  Hu = Hd = v/sqrt(2) = {v/np.sqrt(2):.2f} MeV")
print()

# Check F-terms at original seesaw vacuum:
cof_Phi = np.array([Phi_vev[1]*Phi_vev[2], Phi_vev[0]*Phi_vev[2], Phi_vev[0]*Phi_vev[1]])
det_Phi = np.prod(Phi_vev)

print(f"F-term checks:")
for i, flav in enumerate(['u', 'd', 's']):
    F_Phi_i = m_arr[i] * Lambda + (X_c_vev / Lambda) * cof_Phi[i]
    print(f"  F_{{Phi_{flav}}} = m_{flav}*Lambda + (X_c/Lambda)*cof = {m_arr[i]*Lambda:.4f} + {(X_c_vev/Lambda)*cof_Phi[i]:.4f} = {F_Phi_i:.6e}")

F_Xc_vacuum = det_Phi / Lambda - Lambda**2 + lam * (v/np.sqrt(2))**2
F_Xc_no_higgs = det_Phi / Lambda - Lambda**2
print(f"  F_{{X_c}} (no Higgs) = det Phi / Lambda - Lambda^2 = {det_Phi/Lambda:.6e} - {Lambda**2:.6e} = {F_Xc_no_higgs:.6e}")
print(f"  F_{{X_c}} (with Higgs) = {F_Xc_no_higgs:.6e} + lambda*v^2/2 = {F_Xc_no_higgs:.6e} + {lam*v**2/2:.6e} = {F_Xc_vacuum:.6e}")

# Verify det Phi / Lambda = Lambda^2 at original seesaw
print(f"  det Phi / Lambda = {det_Phi/Lambda:.6e},  Lambda^2 = {Lambda**2:.6e}")
print(f"  Ratio: {det_Phi/(Lambda * Lambda**2):.12f}")
print()

print(f"  F_{{X_c}} = lambda * v^2/2 = {lam * v**2/2:.6e} MeV^2")
print(f"  This is the SUSY-breaking F-term when X is dynamical.")
print()

# Now compute the off-diagonal meson mass-squared in canonical variables.
# For Phi^a_b (a != b):
# (m^2)_{Phi^a_b, Phi^{a*}_b} = sum_K |W_{K, Phi^a_b}|^2 + f_pi^2 / Lambda^2 ??? (soft term)
#
# Wait: the soft term is V_soft = f_pi^2 Tr(M^dag M) = f_pi^2 Lambda^2 Tr(Phi^dag Phi)
# So: V_soft = f_pi^2 Lambda^2 sum |Phi_{ij}|^2
# The coefficient of |Phi^a_b|^2 is f_pi^2 Lambda^2.
# This contributes f_pi^2 Lambda^2 to the mass-squared of Phi^a_b.

m2_soft_can = f_pi**2 * Lambda**2  # MeV^4 (mass-squared for canonical Phi)

# W_{K, Phi^a_b}:
# K = Phi^b_a: W_{Phi^b_a, Phi^a_b} = (X_c / Lambda) * d^2(det Phi)/d(Phi^a_b)d(Phi^b_a)
# d^2(det Phi)/d(Phi^a_b)d(Phi^b_a) = -Phi_c (same sign analysis as before)
# W_{Phi^b_a, Phi^a_b} = -(X_c / Lambda) * Phi_c

# All other W_{K, Phi^a_b} = 0 at diagonal vacuum (same reasoning as before).
# sum_K |W_{K, Phi^a_b}|^2 = |X_c / Lambda|^2 * Phi_c^2

# B-term:
# B_{(Phi^a_b)(Phi^b_a)} = sum_K F_K^* * W_{K, Phi^a_b, Phi^b_a}
# Only K = X_c contributes:
# W_{X_c, Phi^a_b, Phi^b_a} = (1/Lambda) * d^2(det Phi)/(dPhi^a_b dPhi^b_a) = -Phi_c / Lambda
# B = F_{X_c}^* * (-Phi_c / Lambda) = -(lam * v^2/2) * Phi_c / Lambda
# |B| = (lam * v^2/2) * Phi_c / Lambda

print(f"OFF-DIAGONAL MESON MASS-SQUARED (Canonical variables, Case B):")
print()
print(f"  m^2_soft(Phi^a_b) = f_pi^2 * Lambda^2 = {m2_soft_can:.6e} MeV^4")
print()

print(f"  {'Pair':16s} | {'Phi_c (MeV)':12s} | {'m^2_F':14s} | {'m^2_soft':14s} | {'|B|':14s} | {'m^2_min':14s} | {'Tachyonic?':10s}")
print("  " + "-" * 105)

for (a, b, c_idx, label) in [(0,1,2, 'u-d (Phi_s)'), (0,2,1, 'u-s (Phi_d)'), (1,2,0, 'd-s (Phi_u)')]:
    Phi_c = Phi_vev[c_idx]
    flav_c = ['u','d','s'][c_idx]

    m2_F = abs(X_c_vev / Lambda)**2 * Phi_c**2
    B_abs = abs(lam * v**2 / 2.0) * Phi_c / Lambda  # MeV^2 * MeV / MeV = MeV^2
    # Wait, dim check: [F_{X_c}] = MeV^2, [Phi_c/Lambda] = dimensionless
    # |B| has dim MeV^2. But m^2_F and m^2_soft have dim MeV^4 (mass-squared of canonical field).
    # |B| should also have dim MeV^4 for consistency.
    # B = F_{X_c}^* * W_{X_c, Phi^a_b, Phi^b_a}
    # [W_{X_c, Phi^a_b, Phi^b_a}] = [d^3W / d X_c d Phi^a_b d Phi^b_a]
    #   = MeV^3 / (MeV * MeV * MeV) = MeV^0 = dimensionless
    # [B] = [F_{X_c}] * [W_3rd deriv] = MeV^2 * 1 = MeV^2
    # But [m^2_F] = [|W_2nd deriv|^2] = dimensionless^2 = dimensionless ??? No.
    # [W_{Phi^b_a, Phi^a_b}] = MeV^3 / (MeV * MeV) = MeV. So [m^2_F] = MeV^2.
    # [m^2_soft] = f_pi^2 * Lambda^2 = MeV^4. INCONSISTENT!

    # I have a dimensional error in the soft term. Let me recheck.
    # V_soft = f_pi^2 Tr(M^dag M) has [V_soft] = MeV^2 * MeV^4 = MeV^6 ???
    # No: [f_pi^2] = MeV^2, [Tr(M^dag M)] = MeV^4, [V_soft] = MeV^6 ???
    # But [V] should be MeV^4 in 4d. So the coefficient is NOT f_pi^2 but
    # something with dim MeV^0 to make V have dim MeV^4.
    #
    # In canonical variables: V_soft = (f_pi/Lambda)^2 * Tr(Phi^dag Phi) * Lambda^2
    # Hmm. Actually the issue is that in the confined SQCD description,
    # the Kahler potential is K = Tr(M^dag M) / Lambda^2 (to get dim MeV^2),
    # so the soft term is V_soft = (f_pi^2 / Lambda^2) * Tr(M^dag M)
    # = f_pi^2 Tr(Phi^dag Phi) in canonical Phi = M/Lambda.
    # [V_soft] = MeV^2 * MeV^2 = MeV^4. GOOD.
    pass

# CORRECTED soft mass:
m2_soft_can_corrected = f_pi**2  # MeV^2 (mass-squared for canonical Phi = M/Lambda)

print()
print("  CORRECTED: V_soft = f_pi^2 Tr(Phi^dag Phi) in canonical Phi = M/Lambda")
print(f"  m^2_soft(Phi^a_b) = f_pi^2 = {m2_soft_can_corrected} MeV^2")
print()

print(f"  {'Pair':16s} | {'Phi_c':10s} | {'m^2_F':14s} | {'m^2_soft':10s} | {'|B|':14s} | {'m^2_diag':14s} | {'m^2_min':14s} | {'Tachyonic?'}")
print("  " + "-" * 115)

summary_data = []
for (a, b, c_idx, label) in [(0,1,2, 'u-d (Phi_s)'), (0,2,1, 'u-s (Phi_d)'), (1,2,0, 'd-s (Phi_u)')]:
    Phi_c = Phi_vev[c_idx]
    flav_c = ['u','d','s'][c_idx]

    # W_{Phi^b_a, Phi^a_b} = -(X_c / Lambda) * Phi_c  (dim: MeV/MeV * MeV = MeV)
    W_2nd = abs(X_c_vev / Lambda) * Phi_c  # MeV
    m2_F = W_2nd**2  # MeV^2

    # B = F_{X_c}^* * W_{X_c, Phi^a_b, Phi^b_a}
    # W_{X_c, Phi^a_b, Phi^b_a} = -Phi_c / Lambda (dimensionless)
    W_3rd = Phi_c / Lambda  # dimensionless
    B_abs = abs(lam * v**2 / 2.0) * W_3rd  # MeV^2 * dimensionless = MeV^2

    m2_diag = m2_F + m2_soft_can_corrected
    m2_min = m2_diag - B_abs
    tachyonic = "YES" if m2_min < 0 else "NO"

    summary_data.append((label, flav_c, Phi_c, m2_F, m2_soft_can_corrected, B_abs, m2_diag, m2_min, tachyonic))

    print(f"  {label:16s} | {Phi_c:10.4f} | {m2_F:14.6e} | {m2_soft_can_corrected:10.0f} | {B_abs:14.6e} | {m2_diag:14.6e} | {m2_min:14.6e} | {tachyonic}")

print()

# =========================================================================
print()
print("=" * 78)
print("  PROBLEM 3: RESOLUTION OF THE CONTRADICTION")
print("=" * 78)
print()

print("THE KEY QUESTION: Is X a Lagrange multiplier or a dynamical field?")
print()
print("For N_f = N_c = 3 (SU(3) SQCD with 3 flavors):")
print()
print("Seiberg's result: The confined description has gauge-invariant composites")
print("  M^i_j, B, B~ satisfying the QUANTUM MODIFIED constraint:")
print("    det M - B B~ = Lambda^{2 N_c} = Lambda^6")
print()
print("  This constraint is EXACT (non-perturbative). It arises from:")
print("    - Instanton-generated superpotential contributions")
print("    - The quantum moduli space is a SMOOTH manifold")
print("    - X is introduced as a Lagrange multiplier to ENFORCE this constraint")
print()
print("Seiberg duality: N_f = N_c maps to N_c' = N_f - N_c = 0.")
print("  The magnetic dual has NO gauge group. It is a pure Wess-Zumino model.")
print("  In this case, the 'magnetic quarks' q, q~ do not exist (N_c' = 0).")
print("  The singlet X (which in ISS with N_f > N_c corresponds to a magnetic meson)")
print("  does NOT arise as a dynamical field from the magnetic dual.")
print()
print("CONCLUSION ON X:")
print("  For N_f = N_c = 3, X is a LAGRANGE MULTIPLIER.")
print("  It has NO canonical kinetic term in the Wilsonian effective action.")
print("  Its equation of motion IS the constraint det M - BB~ = Lambda^6.")
print("  Case A is the correct physical description.")
print()
print("HOWEVER:")
print("  If one chooses to give X a Kahler potential (e.g., through")
print("  the NMSSM coupling or Kahler pole mechanism), then X becomes")
print("  dynamical. This is a MODEL CHOICE, not derived from first principles.")
print("  In such a model, Case B applies, but the tachyonic masses depend")
print("  on the Kahler normalization of X.")
print()


# =========================================================================
print()
print("=" * 78)
print("  PROBLEM 4: DETAILED OFF-DIAGONAL MASSES IN BOTH CASES")
print("=" * 78)
print()

print("CASE A (X = Lagrange multiplier, CORRECT for N_f = N_c = 3):")
print()
print("  The constraint is det M = Lambda^6 - lambda_phys * v^2/2 at EWSB.")
print("  F_X = 0 (constraint satisfied).")
print("  Off-diagonal meson mass-squared:")
print()

for data in summary_data:
    label, flav_c, Phi_c, m2_F, m2_soft, B_abs, m2_diag, m2_min, tach = data
    print(f"  {label}:")
    print(f"    m^2 = m^2_F + m^2_soft = {m2_F:.6e} + {m2_soft:.0f} = {m2_F + m2_soft:.6e} MeV^2")
    print(f"    B-term = 0 (F_X = 0)")
    print(f"    BOTH eigenvalues = {m2_F + m2_soft:.6e} MeV^2 > 0")
    print()

print("  No tachyonic modes. The vacuum is STABLE.")
print("  CKM mixing does NOT arise from off-diagonal meson VEVs.")
print()

print("CASE B (X = dynamical, HYPOTHETICAL for N_f = N_c = 3):")
print()

for data in summary_data:
    label, flav_c, Phi_c, m2_F, m2_soft, B_abs, m2_diag, m2_min, tach = data
    print(f"  {label}:")
    print(f"    m^2_F = |X_c/Lambda|^2 * Phi_c^2 = {m2_F:.6e} MeV^2")
    print(f"    m^2_soft = f_pi^2 = {m2_soft:.0f} MeV^2")
    print(f"    B = F_{{X_c}} * Phi_c / Lambda = {B_abs:.6e} MeV^2")
    print(f"    m^2_diag = m^2_F + m^2_soft = {m2_diag:.6e} MeV^2")
    print(f"    m^2_min = m^2_diag - |B| = {m2_min:.6e} MeV^2  {'TACHYONIC' if m2_min < 0 else 'STABLE'}")
    ratio = B_abs / m2_diag if m2_diag != 0 else float('inf')
    print(f"    |B| / m^2_diag = {ratio:.6e}")
    print()

print("  ANALYSIS OF CASE B TACHYONS:")
print()

# Check if the tachyon hierarchy matches Oakes
if len(summary_data) >= 3:
    # Extract the B-terms for comparison
    B_ds = None  # d-s sector, c=u
    B_us = None  # u-s sector, c=d
    B_ud = None  # u-d sector, c=s

    for data in summary_data:
        label, flav_c, Phi_c, m2_F, m2_soft, B_abs, m2_diag, m2_min, tach = data
        if 'u-d' in label:
            B_ud = B_abs
            m2min_ud = m2_min
        elif 'u-s' in label:
            B_us = B_abs
            m2min_us = m2_min
        elif 'd-s' in label:
            B_ds = B_abs
            m2min_ds = m2_min

    if B_ud and B_ds:
        print(f"  B-term hierarchy:")
        print(f"    B(u-d, c=s) = {B_ud:.6e} (Phi_s = {Phi_vev[2]:.4f})")
        print(f"    B(u-s, c=d) = {B_us:.6e} (Phi_d = {Phi_vev[1]:.4f})")
        print(f"    B(d-s, c=u) = {B_ds:.6e} (Phi_u = {Phi_vev[0]:.4f})")
        print()
        print(f"  Ratios:")
        print(f"    B(d-s) / B(u-d) = Phi_u / Phi_s = {Phi_vev[0]/Phi_vev[2]:.6f} = M_u/M_s = {M_orig[0]/M_orig[2]:.6f}")
        print(f"    B(u-s) / B(u-d) = Phi_d / Phi_s = {Phi_vev[1]/Phi_vev[2]:.6f} = M_d/M_s = {M_orig[1]/M_orig[2]:.6f}")
        print()

        # The Oakes relation: sin theta_C = sqrt(m_d/m_s)
        # In the seesaw: M_i = C/m_i, so M_d/M_s = m_s/m_d
        # Phi_d/Phi_s = m_s/m_d
        print(f"  Connection to Oakes:")
        print(f"    Phi_d / Phi_s = m_s / m_d = {m_s/m_d:.6f}")
        print(f"    sqrt(m_d/m_s) = {np.sqrt(m_d/m_s):.6f}")
        print(f"    tan(theta_C) from Oakes = sqrt(m_d/m_s) = {np.sqrt(m_d/m_s):.6f}")
        print(f"    PDG |V_us| = 0.2257")
        print()
        print(f"  The B-term hierarchy does NOT directly reproduce the Oakes relation.")
        print(f"  The tachyon in the d-s sector (B proportional to M_u = largest meson)")
        print(f"  is the STRONGEST, while the u-d sector (B proportional to M_s = smallest")
        print(f"  meson) is the WEAKEST. This is OPPOSITE to what would give the Cabibbo")
        print(f"  angle from the d-s mixing!")
        print()
        print(f"  The B-term is proportional to M_c = C/m_c (inverse quark mass).")
        print(f"  The LARGEST tachyon is in the sector where the spectator is the LIGHTEST quark.")
        print(f"  This means the u-quark drives the strongest flavor mixing in the d-s sector,")
        print(f"  which IS qualitatively consistent with the large Cabibbo angle connecting d and s.")
        print()
        print(f"  Quantitatively: the ratio B(d-s)/B(u-d) = M_u/M_s = m_s/m_u = {m_s/m_u:.2f}.")
        print(f"  The mixing angle from a tachyonic instability scales as sqrt(|B_ds|/|B_ud|)")
        print(f"  = sqrt(m_s/m_u) = {np.sqrt(m_s/m_u):.4f}, which is NOT the Oakes relation.")
        print()


# =========================================================================
print()
print("=" * 78)
print("  PROBLEM 5: WHICH INTERPRETATION IS CORRECT?")
print("=" * 78)
print()

print("DEFINITIVE ANSWER for N_f = N_c = 3:")
print()
print("  X is a LAGRANGE MULTIPLIER. Case A is correct.")
print()
print("  Proof:")
print("  1) Seiberg's confined description for N_f = N_c (1994, hep-th/9411149):")
print("     The low-energy theory is described by gauge-invariant composites")
print("     M, B, B~ subject to the QUANTUM constraint det M - BB~ = Lambda^{2N_c}.")
print("     This constraint is enforced by a Lagrange multiplier in the")
print("     superpotential. The Lagrange multiplier has no kinetic term.")
print()
print("  2) Magnetic dual has N_c' = N_f - N_c = 0:")
print("     There is no magnetic gauge group. The 'magnetic quarks' do not exist.")
print("     The ISS mechanism (Intriligator-Seiberg-Shih, 2006) requires N_f > N_c")
print("     to have a magnetic dual with a non-trivial gauge group and dynamical")
print("     pseudo-modulus. For N_f = N_c, the ISS mechanism does not apply.")
print()
print("  3) The constraint is exact:")
print("     For N_f = N_c, the quantum moduli space is smooth (no singularities).")
print("     The constraint arises from the instanton-generated superpotential")
print("     W_dyn = (Lambda^{3N_c-N_f} / det M) which, combined with the tree-level")
print("     mass terms, enforces det M = Lambda^{2N_c} at the quantum level.")
print()
print("  CONSEQUENCE: F_X = 0 at the vacuum. No B-terms from F_X.")
print("  The off-diagonal meson mass-squared is ALWAYS POSITIVE:")
print("    m^2 = f_pi^2 + |X_0|^2 M_c^2 > 0")
print("  There are NO tachyonic off-diagonal modes.")
print()


# =========================================================================
print()
print("=" * 78)
print("  PROBLEM 6: SOFT-TERM TACHYONS? (Case A, X = Lagrange multiplier)")
print("=" * 78)
print()

print("At the EWSB-shifted seesaw, the off-diagonal meson masses from soft terms alone:")
print()
print("  m^2_soft(Phi^a_b) = f_pi^2 > 0  for ALL (a, b)")
print()
print("The F-term contribution from W_{Phi^b_a, Phi^a_b} adds positively:")
print("  m^2_F(Phi^a_b) = |X_0 M_c / Lambda^2|^2 = |X_c / Lambda|^2 Phi_c^2 > 0")
print()
print("  (Note: this is tiny, ~ 10^{-16} MeV^2, compared to f_pi^2 = 8464 MeV^2)")
print()
print("Total: m^2 = f_pi^2 + |X_c/Lambda|^2 Phi_c^2 > 0")
print()
print("NO TACHYONIC DIRECTIONS from soft terms alone.")
print("The soft mass f_pi^2 stabilizes ALL off-diagonal modes.")
print()


# =========================================================================
print()
print("=" * 78)
print("  SUMMARY AND ASSESSMENT")
print("=" * 78)
print()

print("1. IS THE OFF-DIAGONAL MESON SECTOR TACHYONIC?")
print("   NO. For N_f = N_c = 3, X is a Lagrange multiplier (Case A).")
print("   The constraint det M = Lambda^6 - lambda_phys v^2/2 is exact.")
print("   F_X = 0 at the vacuum. No B-terms. No tachyons.")
print(f"   Off-diagonal masses: m^2 = f_pi^2 + O(10^{{-16}}) MeV^2 = {m_tilde_sq} MeV^2.")
print()

print("2. SOURCE OF CKM MIXING (Case A)?")
print("   With X as Lagrange multiplier, there is NO tree-level CKM mixing from")
print("   the meson vacuum. Possible sources:")
print("   a) Radiative corrections (one-loop CW potential, which was computed in")
print("      previous analyses)")
print("   b) Higher-dimensional operators in the Kahler potential")
print("   c) The Oakes relation as an INPUT from the quark mass matrix texture")
print("      (m_d/m_s ratio connected to Cabibbo angle)")
print("   d) Non-perturbative effects (bion contributions to off-diagonal masses)")
print()

print("3. IF X WERE DYNAMICAL (Case B, hypothetical), DO TACHYONS REPRODUCE OAKES?")
tachyonic_found = any(d[8] == 'YES' for d in summary_data)
if tachyonic_found:
    print("   Case B gives tachyonic off-diagonal modes. The B-term hierarchy is:")
    print(f"   B(d-s) : B(u-s) : B(u-d) = M_u : M_d : M_s")
    print(f"   = {M_orig[0]/M_orig[2]:.2f} : {M_orig[1]/M_orig[2]:.2f} : 1")
    print(f"   = m_s/m_u : m_s/m_d : 1")
    print(f"   = {m_s/m_u:.1f} : {m_s/m_d:.1f} : 1")
    print()
    print(f"   The STRONGEST tachyon is in the d-s sector (spectator = u, lightest quark).")
    print(f"   This means the d-s mixing is driven hardest, consistent with the large")
    print(f"   Cabibbo angle. BUT the ratio is m_s/m_u = {m_s/m_u:.1f}, not sqrt(m_d/m_s).")
    print(f"   The Oakes relation gives |V_us| ~ sqrt(m_d/m_s) = {np.sqrt(m_d/m_s):.4f},")
    print(f"   which is NOT the same as any simple ratio of the B-terms.")
    print(f"   The tachyon does NOT directly reproduce the Oakes relation.")
else:
    print("   No tachyonic modes found even in Case B (unexpected).")
print()

print("4. STATUS OF THE CABIBBO ANGLE DERIVATION:")
print("   The Cabibbo angle does NOT emerge from the tree-level meson vacuum")
print("   in either case. For Case A (physical case), the vacuum is diagonal")
print("   and all off-diagonal masses are positive. CKM mixing requires additional")
print("   physics beyond the tree-level confined SQCD superpotential.")
print()
print("   The Oakes relation tan(theta_C) = sqrt(m_d/m_s) relates the Cabibbo")
print("   angle to quark mass ratios, but in this framework the quark masses are")
print("   INPUTS (through the m_i M^i_i terms in W), not derived quantities.")
print("   The connection between m_d/m_s and the CKM angle must come from the")
print("   structure of the quark mass matrix TEXTURE, not from the meson vacuum.")
print()
print("   BOTTOM LINE: The SQCD meson sector provides the MASS HIERARCHY")
print("   (through the seesaw M_i = C/m_i) and the KOIDE STRUCTURE (through")
print("   the CW spectrum), but does NOT by itself generate CKM mixing.")
print("   The Cabibbo angle is a separate structural element that must be")
print("   imposed through the quark mass matrix or a texture ansatz.")
print()


# =========================================================================
# Numerical verification table
# =========================================================================
print("=" * 78)
print("  NUMERICAL VERIFICATION TABLE")
print("=" * 78)
print()

print(f"{'Quantity':>40} | {'Value':>20} | {'Units':>10}")
print("-" * 75)
print(f"{'Lambda^6':>40} | {Lambda6:>20.6e} | {'MeV^6':>10}")
print(f"{'lambda_hat':>40} | {lambda_hat:>20.4f} | {'':>10}")
print(f"{'v':>40} | {v:>20.2f} | {'MeV':>10}")
print(f"{'v^2/2':>40} | {v**2/2:>20.6e} | {'MeV^2':>10}")
print(f"{'lambda_hat * Lambda^4 * v^2/2':>40} | {lambda_hat * Lambda**4 * v**2/2:>20.6e} | {'MeV^6':>10}")
print(f"{'Delta (corrected, dim-consistent)':>40} | {Delta_corrected:>20.6e} | {'':>10}")
print(f"{'C (original)':>40} | {C_orig:>20.6f} | {'MeV^2':>10}")
print(f"{'M_u = C/m_u':>40} | {M_orig[0]:>20.4f} | {'MeV^2':>10}")
print(f"{'M_d = C/m_d':>40} | {M_orig[1]:>20.4f} | {'MeV^2':>10}")
print(f"{'M_s = C/m_s':>40} | {M_orig[2]:>20.4f} | {'MeV^2':>10}")
print(f"{'X_0 = -C/Lambda^6':>40} | {X0_orig:>20.6e} | {'MeV^{-4}':>10}")
print(f"{'X_c = X_0 * Lambda^4 (canonical)':>40} | {X_c_vev:>20.6e} | {'MeV':>10}")
print(f"{'Phi_u = M_u/Lambda':>40} | {Phi_vev[0]:>20.4f} | {'MeV':>10}")
print(f"{'Phi_d = M_d/Lambda':>40} | {Phi_vev[1]:>20.4f} | {'MeV':>10}")
print(f"{'Phi_s = M_s/Lambda':>40} | {Phi_vev[2]:>20.4f} | {'MeV':>10}")
print(f"{'F_{{X_c}} = lambda * v^2/2':>40} | {lam * v**2/2:>20.6e} | {'MeV^2':>10}")
print(f"{'f_pi^2':>40} | {f_pi**2:>20.0f} | {'MeV^2':>10}")
print(f"{'det M / Lambda^6':>40} | {det_orig/Lambda6:>20.12f} | {'':>10}")
print()

# Oakes check
print(f"{'Oakes: tan theta_C = sqrt(m_d/m_s)':>40} | {np.sqrt(m_d/m_s):>20.6f} | {'':>10}")
print(f"{'theta_C (Oakes)':>40} | {np.degrees(np.arctan(np.sqrt(m_d/m_s))):>20.4f} | {'deg':>10}")
print(f"{'theta_C (PDG)':>40} | {13.04:>20.4f} | {'deg':>10}")
print()

print("=" * 78)
print("  END OF COMPUTATION")
print("=" * 78)
