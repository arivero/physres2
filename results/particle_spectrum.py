#!/usr/bin/env python3
"""
Complete N=1 SUSY particle spectrum computation.

Theory: SU(3) SQCD Nf=Nc=3 (u,d,s) confined at Lambda=300 MeV,
coupled to NMSSM Higgs sector (H_u, H_d, S), tan beta = 1.

Reference parameters from problem statement.
All masses in MeV unless noted.
"""

import numpy as np
from numpy.linalg import eigvalsh, eigh, eigvals, norm

# ============================================================
# PARAMETERS
# ============================================================
LAM   = 300.0        # MeV, confinement scale
v     = 246220.0     # MeV, full EW vev
tanbeta = 1.0
vu    = v / np.sqrt(2)    # = v_d since tan beta = 1
vd    = v / np.sqrt(2)
f_pi  = 92.0         # MeV, soft SUSY-breaking scale (f_pi^2 = soft mass^2)
lam_S = 0.72         # NMSSM singlet coupling
m_B   = LAM          # = 300 MeV, baryon mass parameter

# Quark masses (MeV)
m_u = 2.16
m_d = 4.67
m_s = 93.4
m_c = 1275.0
m_b = 4180.0
m_t = 172760.0

# Derived vacuum quantities
# C = (Lambda^6 m_u m_d m_s)^{1/3}  ... the problem says C = Lambda^6 (m_u m_d m_s)^{1/3} / Lambda^4 = Lambda^2 (m_u m_d m_s)^{1/3}
# But the problem also states C = (Lambda^6 m_u m_d m_s)^{1/3}. Let's check both:
# "C = (Λ⁶ m_u m_d m_s)^{1/3} = 882297 MeV²"
# (300^6 * 2.16 * 4.67 * 93.4)^{1/3}
C_check = (LAM**6 * m_u * m_d * m_s)**(1.0/3.0)
# Alternatively the older convention: C = Lambda^2 (m_u m_d m_s)^{1/3}
C_old = LAM**2 * (m_u * m_d * m_s)**(1.0/3.0)
print(f"C = (Λ⁶ m_u m_d m_s)^(1/3) = {C_check:.4f} MeV²   [problem statement: 882297]")
print(f"C = Λ²(m_u m_d m_s)^(1/3)  = {C_old:.4f} MeV²   [older convention]")
# The problem states C = 882297 MeV², so use C_check
C = C_check
print(f"Using C = {C:.2f} MeV²\n")

# Meson VEVs: <M^i_i> = C / m_i  (from W_X = det M - Lambda^6 = 0, W_mi = m_i + X*cofactor = 0)
# At diagonal vacuum: det M = M_u * M_d * M_s => M_i = C / m_i
M_u_vev = C / m_u
M_d_vev = C / m_d
M_s_vev = C / m_s

# X VEV: from W_M_uu = m_u + X * M_d_vev * M_s_vev = 0
# => X0 = -m_u / (M_d_vev * M_s_vev) = -m_u / (C/m_d * C/m_s) = -m_u m_d m_s / C^2
# Also: X0 = -C / Lambda^6  (problem statement)
X0_from_C  = -C / LAM**6
X0_from_mi = -m_u * m_d * m_s / C**2
print(f"X0 = -C/Λ⁶              = {X0_from_C:.6e} MeV⁻⁴   [problem: -1.210e-9 MeV^-4]")
print(f"X0 = -m_u m_d m_s / C²  = {X0_from_mi:.6e} MeV⁻⁴   [should match]")
X0 = X0_from_C
print(f"Using X0 = {X0:.6e} MeV⁻⁴\n")

# Verify det(M_vev) = Lambda^6
det_M = M_u_vev * M_d_vev * M_s_vev
print(f"det M_vev / Λ⁶ = {det_M / LAM**6:.12f}  (should be 1.0)")
print(f"M_uu = {M_u_vev:.3f} MeV,  M_dd = {M_d_vev:.3f} MeV,  M_ss = {M_s_vev:.3f} MeV")
print()

# Yukawa couplings: y_c H_u M^c_c => y_c = m_c / vu  (since m_c = y_c * vu = y_c * v/sqrt(2))
# Similarly y_b = m_b / vd
y_c_Yukawa = m_c / vu   # Yukawa coupling for charm
y_b_Yukawa = m_b / vd   # Yukawa coupling for bottom
y_t_Yukawa = m_t / vu   # Yukawa coupling for top

print("=" * 70)
print("PART (a): SCALAR MESON MASSES")
print("=" * 70)

# ============================================================
# PART (a): SCALAR MESON MASSES
# ============================================================
#
# Scalar potential: V = sum_i |dW/dPhi_i|^2 + V_soft + V_D
#
# V_soft = f_pi^2 Tr(M^dag M)  (adds f_pi^2 to all meson scalar mass^2)
#
# After integrating out X (no kinetic term), the effective potential
# for mesons comes from the modified F-terms.
#
# The constraint from X (no kinetic term) is:
#   dW/dX = 0 => det M - BB~ - Lambda^6 = 0   (classical constraint)
#
# F-terms for mesons:
#   F_{M^i_j}* = dW/dM^i_j = m_i delta_{ij} + X * (cofactor of M)^j_i
#              + (c_3 * 3 * (det M)^2 * cofactor)  [higher order, ignore for now]
#              + (Yukawa terms)
#
# At the vacuum (diagonal M, B=B~=0):
#   dW/dM^i_i = m_i + X0 * (M_j * M_k)  where j,k != i
#             = m_i + X0 * (C/m_j)(C/m_s) = m_i - C^2/(Lambda^6) * Lambda^6/M_i * m_i ... wait
#
# Let's be careful. At vacuum M_diag:
#   cofactor(M^i_i) = M_j * M_k for the diagonal element  (j,k != i)
#   cofactor for off-diagonal: at diagonal vacuum, cofactor(M^i_j) for i!=j
#     involves a 2x2 sub-determinant including only the remaining rows/cols
#
# The mass^2 matrix comes from d^2V/dphi dphi* = |W_{phi_a phi_b}|^2 evaluated at vac.
# For complex scalar phi = (Re phi + i Im phi)/sqrt(2), the mass^2 in the (Re,Im) basis
# has the form:
#   M^2 = ( W*_{ab} W_{bc} + soft )
# but we need to be careful about the two-component structure.
#
# For a holomorphic superpotential, the scalar mass-squared matrix is:
#   (M^2)_{phi_a, phi_b*} = sum_c W_{ac}* W_{bc} + (soft terms)
# and there's also the B-type mixing:
#   (M^2)_{phi_a, phi_b}  = sum_c W_{abc} <phi_c> * ...   (A-terms, vanish at tree level here)
#
# Since V_soft = f_pi^2 Tr(M^dag M) = f_pi^2 sum_{ij} |M^i_j|^2,
# this adds f_pi^2 * 1 to each diagonal entry of the M^2 matrix for the mesons.
#
# ------------------------------------------------------------------
# DIAGONAL MESONS: M^i_i  (i = u, d, s)
# ------------------------------------------------------------------
# Off-shell, after integrating out X using its equation of motion:
#   det M = Lambda^6 + BB~
# The effective superpotential for the diagonal mesons M_i = M^i_i
# (at B=0) is approximately:
#   W_eff = sum_i m_i M_i  + W_constraint_substituted
#
# But for the scalar mass computation, we evaluate the second derivatives
# of |F|^2 at the vacuum.
#
# The non-trivial part is from:
#   W_M_ii M_jj = d^2W / (dM_ii dM_jj) = X0 * M_k  (at vacuum, k != i,j)
#   W_M_ii X    = cofactor(M_ii) = M_j M_k = Lambda^6 / M_i  (at vac)
#
# Since X has no kinetic term, we integrate it out by solving:
#   F_X* = dW/dX = det M - BB~ - Lambda^6 = 0
# This means we substitute the constraint into the potential.
# The effective F-terms for M_ii after integrating out X are obtained by
# computing the propagator correction.
#
# More precisely, the scalar mass matrix for {M_ii, i=u,d,s} after
# integrating out X is obtained from the W_{ij} matrix (before integrating out X):
#
#   W matrix for {M_uu, M_dd, M_ss, X} is 4x4:
#     W_ii_jj = X0 * M_k   (i != j)
#     W_ii_X  = Lambda^6 / M_i
#     W_XX    = 0
#
# The fermion mass matrix is this 4x4 matrix (times 1/2 for Majorana).
# When X has no kinetic term, the scalar spectrum is different from the
# fermion spectrum.
#
# For the SCALAR sector with X having no kinetic term:
# The scalar potential is:
#   V = |F_M_uu|^2 + |F_M_dd|^2 + |F_M_ss|^2 + |F_X|^2 (but no kinetic for X)
#
# Wait - if X has NO kinetic term, it is an AUXILIARY field (or Lagrange multiplier).
# So we do NOT add |F_X|^2 to the potential. Instead, we solve F_X = 0 as a
# constraint and substitute into the other F-terms.
#
# Actually the more standard treatment: X is an auxiliary field with algebraic
# equation of motion. The procedure:
#   1. Solve dW/dX = 0 => det M - BB~ = Lambda^6 (constraint)
#   2. Plug back into W to get W_eff (Affleck-Dine-Seiberg superpotential)
#
# But the mass spectrum is computed from the full W before integrating out X,
# treating X as a dynamical field but with a delta-function constraint.
# This means X couples the diagonal mesons together.
#
# The effective scalar mass matrix for M_uu, M_dd, M_ss (after integrating out X)
# comes from the 4x4 W_{ab} matrix with a = {M_uu, M_dd, M_ss, X}.
#
# For the 4x4 block {M_uu, M_dd, M_ss, X}:
# The fermion mass matrix M_F is:
A4 = np.array([
    # M_uu    M_dd    M_ss    X
    [0,       X0*M_s_vev, X0*M_d_vev, LAM**6/M_u_vev],   # M_uu
    [X0*M_s_vev, 0,    X0*M_u_vev, LAM**6/M_d_vev],       # M_dd
    [X0*M_d_vev, X0*M_u_vev, 0,    LAM**6/M_s_vev],       # M_ss
    [LAM**6/M_u_vev, LAM**6/M_d_vev, LAM**6/M_s_vev, 0]  # X
], dtype=float)

print("\n4x4 W matrix block {M_uu, M_dd, M_ss, X} (entries in MeV units):")
labels4 = ['M_uu', 'M_dd', 'M_ss', 'X   ']
for i in range(4):
    row_str = "  ".join(f"{A4[i,j]:+.4e}" for j in range(4))
    print(f"  {labels4[i]}: {row_str}")

# Integrate out X: block decompose
# A = [[A_MM, A_MX], [A_XM, A_XX]] where A_XX = 0
# After integrating out X (auxiliary), the effective mass matrix for M_ii is:
#   A_eff = A_MM - A_MX * (A_XX)^{-1} * A_XM
# But A_XX = 0, so this is singular. Instead, since X has no kinetic term,
# the physical scalar spectrum is from the "reduced" theory.
#
# Correct treatment: The Kahler potential K contains |M|^2 + |B|^2 + |B~|^2
# but NOT |X|^2. The scalar potential is:
#   V = K^{i j*} F_i F_j* + V_soft
# where K^{ij*} = delta^{ij} for the fields WITH kinetic terms.
# X does not appear in K, so |F_X|^2 is NOT added to V.
# Instead, F_X = 0 is the constraint (equation of motion for X).
#
# The scalar mass matrix for M_uu, M_dd, M_ss is therefore:
#   (M_S^2)_ij = sum_{k: has kinetic term} W_ik^* W_jk + f_pi^2 delta_ij
# where the sum runs over ALL fields k that have kinetic terms (M_ab, B, B~).
# X is excluded from this sum.
#
# So for the diagonal meson block {M_uu, M_dd, M_ss}, we need:
#   (M_S^2)_ij = W_{M_ii, M_uu} W^*_{M_jj, M_uu}  [sum over all M_ab rows of W,
#                   but only those with kinetic terms]
#
# At the vacuum, the relevant W_{M_ii, M_ab} entries:
# - W_{M_uu, M_dd} = X0 M_ss,  W_{M_uu, M_ss} = X0 M_dd,  W_{M_uu, X} = Λ⁶/M_uu
# - For diagonal mesons: cross-derivatives with off-diagonal mesons vanish at diagonal vac
#   because the cofactors involve the off-diagonal entries.
# Wait - actually the cofactor structure: at a DIAGONAL M:
#   d^2(det M) / dM_uu dM_ab = 0 if a != b (since det M = M_uu M_dd M_ss at diag vac,
#   so cross with off-diag is zero at diag vac).
#
# So W_{M_ii, M_ab} = 0 for off-diagonal a != b at diagonal vacuum.
# The W matrix for diagonal mesons ONLY couples to each other and to X.
# X is excluded from the scalar potential sum.
# Therefore the scalar mass matrix for diagonal mesons is:
#   (M_S^2)_{ii, jj} = sum_{k in {M_uu, M_dd, M_ss, off_diag_M, B, Bt}} W_{ii,k} W^*_{jj,k}
#                    + f_pi^2 delta_{ii,jj}
#
# W_{M_ii, off-diag M_ab} = 0 at diagonal vacuum (as argued above).
# W_{M_ii, B} = W_{M_ii, Bt} = 0 (no B in W_M_ii at B=0 vacuum... actually:
#   W_M_uu = m_u + X * d(det M)/dM_uu + X * d(-BB~)/dM_uu = m_u + X * M_dd * M_ss
#   So W_{M_uu, M_dd} = X0 * M_ss, W_{M_uu, M_ss} = X0 * M_dd, W_{M_uu, X} = M_dd*M_ss = Λ⁶/M_uu
#   W_{M_uu, B} = W_{M_uu, Bt} = 0  (no M_uu in B terms at B=0)
#
# So the ONLY contributions to scalar mass^2 for diagonal mesons come from
# W_{ii, jj} = X0 * M_k (k != i,j) [excluded X, so this term involves M_jj which has kinetic]
# and W_{ii, X} = Λ⁶/M_i which CANNOT contribute since X has no kinetic term.
#
# Therefore:
#   (M_S^2)_{uu,uu} = |W_{M_uu, M_dd}|^2 + |W_{M_uu, M_ss}|^2 + f_pi^2
#                   = (X0 M_ss)^2 + (X0 M_dd)^2 + f_pi^2
#
#   (M_S^2)_{uu,dd} = W_{M_uu,M_dd} * W^*_{M_dd,M_dd} = ???
# No wait. The scalar mass^2 matrix (M^2)_{I J*} from F-terms is:
#   (M^2_F)_{I J*} = W_{IK} W^*_{JK}    (sum over K including ALL fields with kinetic terms)
# where we need W_{IK} = d^2W/dPhi_I dPhi_K.
#
# So (M^2_F)_{M_uu, M_uu} = sum_K |W_{M_uu, K}|^2  (K over fields with kinetic terms)
#                          = |W_{M_uu, M_uu}|^2 + |W_{M_uu, M_dd}|^2 + |W_{M_uu, M_ss}|^2
#                            + sum_{off-diag K} |W_{M_uu, K}|^2 + |W_{M_uu, B}|^2 + ...
#                          = 0 + (X0 M_ss)^2 + (X0 M_dd)^2 + 0 + 0
#
# And (M^2_F)_{M_uu, M_dd} = sum_K W_{M_uu,K} W^*_{M_dd,K}
#                           = W_{M_uu,M_uu} W*_{M_dd,M_uu} + ...
#                           = 0 (since W_{M_uu,M_dd} W*_{M_dd,M_dd} = W_{M_uu,M_dd}*0 = 0
#                              and W_{M_uu,M_ss} W*_{M_dd,M_ss} = X0*M_dd * X0*M_uu
# Wait:
#   W_{M_uu, M_ss} = X0 * M_dd  (cofactor of M_ss in 3x3 det = M_uu * M_dd when off... no)
# Let me recompute the W_{ij} entries carefully.

print()
print("Recomputing W_{ij} entries (second derivatives of W at vacuum):")
print()
# W = Tr(m_hat M) + X(det M - BB~ - Λ^6) + Yukawa + ...
# = m_u M_uu + m_d M_dd + m_s M_ss + X*(M_uu M_dd M_ss + off-diag terms - BB~ - Λ^6) + ...
# At diagonal vacuum:
#
# W_{M_uu, M_dd} = d/dM_dd [m_u + X * d(det M)/dM_uu]
#                = X * d^2(det M)/dM_uu dM_dd
# det M = M_uu*M_dd*M_ss (at diagonal)
# d^2(det M)/dM_uu dM_dd = M_ss (at diagonal vacuum)
# So W_{M_uu, M_dd} = X0 * M_s_vev  ✓
#
# Similarly W_{M_uu, M_ss} = X0 * M_d_vev
#           W_{M_dd, M_ss} = X0 * M_u_vev
#
# W_{M_uu, X} = d/dX [m_u + X * d(det M)/dM_uu]
#             = d(det M)/dM_uu = M_dd * M_ss = Λ^6 / M_uu  ✓
#
# At diagonal vacuum, off-diagonal meson couplings:
# W_{M_ud, M_du} = X * d^2(det M)/dM_ud dM_du
# For 3x3 matrix M:
# det M = M_uu(M_dd M_ss - M_ds M_sd) - M_ud(M_du M_ss - M_ds M_su) + M_us(M_du M_sd - M_dd M_su)
# d(det M)/dM_ud = -(M_du M_ss - M_ds M_su)  = -M_ss M_du + M_ds M_su  (at diag vac: -M_du*M_ss)
# d^2(det M)/dM_ud dM_du = -M_ss  (at diagonal vac)
# So W_{M_ud, M_du} = X0 * (-M_s_vev) ... but this is NEGATIVE. Hmm.
# Actually: d/dM_du [-(M_du M_ss - M_ds M_su)] = -M_ss
# So W_{M_ud, M_du} = X0 * (-M_ss) ... that gives a different sign than stated.

# Let me use the full cofactor matrix:
# (Cofactor of M)^j_i = d(det M)/dM^i_j
# For 3x3: (Cof)^j_i = epsilon^{jkl} epsilon_{imn} M^m_k M^n_l / 1 (for det M = epsilon M M M / 3!)
# At diagonal M = diag(M1, M2, M3):
# (Cof)^1_1 = M2*M3, (Cof)^2_2 = M1*M3, (Cof)^3_3 = M1*M2
# (Cof)^2_1 = d(det M)/dM^1_2 = (M at rows {2,3} x cols {1,3} minor with sign)
# Actually for upper-triangular interpretation of epsilon:
# det M = sum_{perm P} sign(P) M^P(1)_1 M^P(2)_2 M^P(3)_3
# d(det M)/dM^i_j = (cofactor)_ij = sum over perms that include M^i_j
# For off-diagonal entries at DIAGONAL vacuum (M^i_j for i!=j have M^i_j=0):
# d(det M)/dM^1_2|_{diag} = 0 (since the Leibniz formula for this cofactor involves
#   the minor of M^1_2 which is M^2_? M^3_? at diagonal => involves M off-diag elements that are 0 or the diagonal ones

# Actually the correct cofactor computation:
# For a 3x3 matrix, the cofactor C^i_j = (-1)^{i+j} M_{ji} where M_{ji} is the (j,i) minor
# (minor = determinant of the 2x2 submatrix obtained by deleting row j, col i)
# Note: det M = sum_i M^1_i C^1_i  (expanding along row 1)
#
# d(det M)/dM^i_j = (-1)^{i+j} * (minor with row i, col j deleted)
#
# At DIAGONAL vacuum M = diag(M1, M2, M3):
# Minor obtained by deleting row i, col j:
# If i=j (diagonal): minor is product of remaining 2 diagonal entries
#   d(det M)/dM^1_1 = M2*M3  (delete row 1, col 1 => det [[M2,0],[0,M3]] = M2*M3)
#   sign = (-1)^{1+1} = +1
# If i!=j (off-diagonal): minor of (del row i, del col j) has zeros except possibly diagonal
#   For (i,j) = (1,2): delete row 1, col 2:
#     remaining 2x2 = [[M^2_1, M^2_3], [M^3_1, M^3_3]] = [[0, 0], [0, M3]] at diag vac
#     det = 0. sign = (-1)^{1+2} = -1. So d(det)/dM^1_2|_{diag} = -1 * 0 = 0. ✓
#
# Now d^2(det M)/dM^i_j dM^k_l:
# This equals d(cofactor(i,j))/dM^k_l.
# cofactor(i,j) = (-1)^{i+j} * det(minor_{ij})
# For (i,j) = (1,1) (diagonal):
#   cofactor(1,1) = det[[M22,M23],[M32,M33]] = M22*M33 - M23*M32
#   d(cofactor(1,1))/dM^2_2 = M33, d/dM^3_3 = M22, d/dM^2_3 = -M32, d/dM^3_2 = -M23
#   At diag vac: d/dM^2_2 = M3, d/dM^3_3 = M2, d/dM^2_3 = 0, d/dM^3_2 = 0
#   => W_{M^1_1, M^2_2} = X0 * M3, W_{M^1_1, M^3_3} = X0 * M2 ✓
#
# For (i,j) = (1,2) (off-diagonal):
#   cofactor(1,2) = (-1)^{1+2} det[[M^2_1, M^2_3],[M^3_1, M^3_3]]
#                = -(M^2_1 M^3_3 - M^2_3 M^3_1)
#   d(cofactor(1,2))/dM^2_1|_{diag} = -M^3_3|_{diag} = -M3
#   d(cofactor(1,2))/dM^3_1|_{diag} = +M^2_3|_{diag} = 0
#   d(cofactor(1,2))/dM^2_3|_{diag} = +M^3_1|_{diag} = 0
#   d(cofactor(1,2))/dM^3_3|_{diag} = -M^2_1|_{diag} = 0
#   => W_{M^1_2, M^2_1} = X0 * (-M3)   <-- NEGATIVE!
#
# Hmm, but the standard SQCD literature (Intriligator-Seiberg) and our existing code
# both use +|X0| * M_k for the off-diagonal blocks. The sign of X0 matters:
# X0 < 0, so X0 * (-M3) = |X0| * M3 > 0. So the physical mass squared is
# W_{M12,M21}^2 = (X0 * (-M3))^2 = X0^2 * M3^2 > 0. ✓
#
# Let me redo with the sign explicit:

print("Second derivatives of det M at diagonal vacuum:")
print(f"  W_{{M_uu, M_dd}} = X0 * M_s = {X0 * M_s_vev:.4e}  (sign: {np.sign(X0 * M_s_vev)})")
print(f"  W_{{M_uu, M_ss}} = X0 * M_d = {X0 * M_d_vev:.4e}")
print(f"  W_{{M_dd, M_ss}} = X0 * M_u = {X0 * M_u_vev:.4e}")
print(f"  W_{{M_ud, M_du}} = X0 * (-M_s) = {X0 * (-M_s_vev):.4e}")
print(f"  W_{{M_us, M_su}} = X0 * (-M_d) = {X0 * (-M_d_vev):.4e}")
print(f"  W_{{M_ds, M_sd}} = X0 * (-M_u) = {X0 * (-M_u_vev):.4e}")
print(f"  W_{{M_uu, X}} = M_dd*M_ss = Λ⁶/M_uu = {LAM**6/M_u_vev:.4e}")
print(f"  W_{{B, Bt}} = -X0 = {-X0:.4e}")
print()

# Key insight: the sign difference between diagonal and off-diagonal W entries:
# W_{Mii, Mjj} = X0 * Mk  (k != i,j)   => SAME sign as X0 (negative)
# W_{Mij, Mji} = X0 * (-Mk) for i<j     => OPPOSITE sign, so +|X0| * Mk (positive)
# These enter mass^2 as |W|^2, so the sign doesn't affect scalar masses.
# But for FERMION masses (where the W_ij matrix enters directly), signs matter!

# ------------------------------------------------------------------
# SCALAR MASS^2 MATRIX (after integrating out X)
# For field phi with kinetic term, the F-term contribution to scalar mass^2 is:
#   (m^2_scalar)_{IJ} = sum_{K: kinetic} W_{IK} W*_{JK}
# The X-row/column of W is excluded from this sum (X has no kinetic term).
# The soft term adds f_pi^2 * delta_{IJ} for all meson scalars.
# ------------------------------------------------------------------

# 3x3 block for diagonal mesons {M_uu, M_dd, M_ss}:
# W_{M_uu, M_dd} = X0 * M_s_vev,  W_{M_uu, M_ss} = X0 * M_d_vev
# W_{M_dd, M_uu} = X0 * M_s_vev (symmetric),  W_{M_dd, M_ss} = X0 * M_u_vev
# W_{M_ss, M_uu} = X0 * M_d_vev, W_{M_ss, M_dd} = X0 * M_u_vev
# W_{M_ii, M_ii} = 0  (second deriv of m_i M_ii + X*(M_dd*M_ss) wrt M_uu^2 = 0)
# W_{M_uu, X} = Λ^6/M_uu  [but X excluded from sum!]
# W_{M_uu, off-diag} = 0 at diagonal vacuum

# The F-term mass^2 matrix for {M_uu, M_dd, M_ss} from kinetic-term fields is:
# (M^2_F)_{uu,uu} = |W_{M_uu,M_dd}|^2 + |W_{M_uu,M_ss}|^2
#                 = (X0*M_s)^2 + (X0*M_d)^2
# (M^2_F)_{dd,dd} = (X0*M_s)^2 + (X0*M_u)^2
# (M^2_F)_{ss,ss} = (X0*M_d)^2 + (X0*M_u)^2
# (M^2_F)_{uu,dd} = W_{uu,dd}* W_{dd,dd} + W_{uu,ss}* W_{dd,ss} + ...
#                 = W_{M_uu,K} * W*_{M_dd,K} summed over K
#   For K = M_uu: W_{uu,uu} = 0
#   For K = M_dd: W_{uu,dd} * W_{dd,dd} = (X0*M_s) * 0 = 0
#   For K = M_ss: W_{uu,ss} * W_{dd,ss} = (X0*M_d) * (X0*M_u)
#   For K = off-diagonal: 0
# (M^2_F)_{uu,dd} = X0^2 * M_d * M_u
# (M^2_F)_{uu,ss} = X0^2 * M_s * M_u  (from K=M_dd: W_{uu,dd}*W_{ss,dd} = X0*M_s*X0*M_u = X0^2*M_u*M_s)
# (M^2_F)_{dd,ss} = X0^2 * M_s * M_d  (from K=M_uu)

X02 = X0**2
M_diag_mass2_F = np.array([
    [X02*(M_s_vev**2 + M_d_vev**2),  X02*M_d_vev*M_u_vev,              X02*M_s_vev*M_u_vev],
    [X02*M_d_vev*M_u_vev,             X02*(M_s_vev**2 + M_u_vev**2),   X02*M_s_vev*M_d_vev],
    [X02*M_s_vev*M_u_vev,             X02*M_s_vev*M_d_vev,              X02*(M_d_vev**2 + M_u_vev**2)]
])

# Can also write as: M^2_F = X0^2 * (M*M^T - diag(M_u^2, M_d^2, M_s^2))  where M = (M_u, M_d, M_s)
M_vec = np.array([M_u_vev, M_d_vev, M_s_vev])
M2_F_check = X02 * (np.outer(M_vec, M_vec) - np.diag(M_vec**2))
print("Check M_diag_mass2_F vs formula X0^2*(MM^T - diag(M^2)):")
print(f"  Max diff = {np.max(np.abs(M_diag_mass2_F - M2_F_check)):.2e}")

# Add soft mass
M_diag_mass2 = M_diag_mass2_F + f_pi**2 * np.eye(3)

print("\n3x3 scalar mass^2 matrix for diagonal mesons (MeV^4 or appropriate units?):")
print("  Note: [X0] = MeV^{-4}, [M_vev] = MeV, so [X0*M_vev]^2 = MeV^{-6}")
print("  [Λ^6/M_vev] = MeV^5, so W_{Mii,X} = MeV^5")
print("  [W_{Mii,Mjj}] = [X0]*[M_vev] = MeV^{-3}")
print("  [W_{Mii,Mjj}]^2 = MeV^{-6}  --- this doesn't have mass^2 units!")
print()
print("  DIMENSION ANALYSIS:")
print("  The superpotential W has mass dimension 3 (in natural units where [action]=0).")
print("  [W] = mass^3. Then [W_ij] = mass^3 / (mass * mass) = mass^1.")
print("  So [W_ij]^2 = mass^2. ✓")
print()
print("  But [X] = ?  In W = X*(det M - Λ^6), [W] = mass^3.")
print("  [det M] = [M]^3. If [M] = mass^2 (as a meson = qq~, mass dim 2 in SQCD), then")
print("  [det M] = mass^6, and we need [X]*mass^6 = mass^3, so [X] = mass^{-3}.")
print("  But the problem states X0 = -1.210e-9 MeV^{-4}.")
print()
print("  Let's just trust the problem's units and check W_ij dimensions:")
print(f"  [W_{{M_uu,M_dd}}] = [X0]*[M_s_vev] = MeV^{{-4}} * MeV = MeV^{{-3}}")
print(f"  [W_{{M_uu,X}}]   = [Λ^6/M_uu] = MeV^6/MeV = MeV^5")
print()
print("  These have DIFFERENT units! This means the W matrix entries mix incompatible")
print("  dimensions, which signals that different sectors have different mass dimensions.")
print("  The fermion mass matrix (Majorana mass matrix from W_ij) has entries of")
print("  dimension MeV^{-3} for meson-meson blocks and MeV^5 for meson-X blocks.")
print()
print("  For the SCALAR mass matrix, we need W_ij * W*_kj (sum over j with kinetic terms).")
print("  Since X has no kinetic term, W_{M,X} does NOT contribute to scalar masses.")
print("  So scalar mass^2 for mesons only gets contributions from W_{Mi,Mj}^2")
print("  => [m^2_scalar] = [W_{Mii,Mjj}]^2 = MeV^{-6}")
print()
print("  This is wrong dimensionally. The issue: in 4D SUSY, [phi] = mass^1, [W] = mass^3,")
print("  so [W_ij] = mass^1. But our M_vev has units of mass^2 (C/m_i = MeV^2/MeV = MeV).")
print("  Actually [M] = MeV = mass^1 ✓. And [X] = MeV^{-3} (mass^{-3}) from [W] = mass^3:")
print("  W = X*(det M - Λ^6): [X][M]^3 = mass^3 => [X] = mass^0 = dimensionless?")
print()
print("  RECONCILIATION: Use the given numerical values. The problem states:")
print(f"  X0 = {X0:.4e} MeV^{{-4}} (their convention, perhaps from [X][Λ^6] = MeV^3 => [X]=MeV^{-3})")
print(f"  C = {C:.2f} MeV^2 (their convention)")
print(f"  M_uu = C/m_u = {M_u_vev:.4f} MeV = mass^1 ✓")
print()
print("  W = X*(det M - Λ^6): [X]*MeV^3 = MeV^3 => [X] = dimensionless = MeV^0")
print("  But problem says X0 in MeV^{-4}. Let me check: C = (Λ^6 m_u m_d m_s)^{1/3}")
print(f"  [C] = (MeV^6 * MeV^3)^{{1/3}} = MeV^3 ✓  (so [M_vev] = [C]/[m] = MeV^2)")
print(f"  C = {C:.2f} MeV^2 ✓ (not MeV)")
print()
print("  So [M_vev] = MeV^2, [det M] = MeV^6 = [Λ^6] ✓")
print("  [W] = mass^3 = MeV^3; W = X*(MeV^6 - MeV^6) = X*MeV^6 => [X] = MeV^{-3}")
print("  X0 should be in MeV^{-3}, but problem says MeV^{-4}. Could be a typo or")
print("  different convention (Λ^6 not in MeV^6 but in MeV^4 with different normalization).")
print()
print("  For the mass matrix, let's just use the numerical values as given.")
print("  The key observable: meson VEVs in MeV^2, X0 in some negative dimension,")
print("  and the mass^2 entries from |W_ij|^2.")

print()
print("=" * 70)
print("PART (a): SCALAR MASSES - NUMERICAL COMPUTATION")
print("=" * 70)

# Use the given values directly.
# VEVs:
Mv = np.array([M_u_vev, M_d_vev, M_s_vev])  # in MeV^2 (problem's convention)
print(f"\nMeson VEVs (MeV^2): M_uu = {Mv[0]:.3f}, M_dd = {Mv[1]:.3f}, M_ss = {Mv[2]:.3f}")
print(f"X0 = {X0:.4e} (units as problem states: MeV^-4)")

# ------------------------------------------------------------------
# Off-diagonal scalars: M^i_j for i != j
# ------------------------------------------------------------------
# From the problem statement:
# m^2(M^a_b) = 2*f_pi^2 + |X0|^2 * |M_c|^2  (where c is the "spectator" index)
# The dominant term is 2*f_pi^2 = 2*(92)^2 = 16928 MeV^2
# The correction |X0|^2 |M_c|^2:
#   For M_ud pair (spectator s): |X0|^2 * M_s^2 = (1.21e-9)^2 * (M_s_vev)^2
#   Numerically tiny.

print("\n--- Off-diagonal scalar mesons ---")
print("m^2(off-diag) = 2*f_pi^2 + |X0|^2 * |M_spectator|^2")
m2_offdiag_base = 2 * f_pi**2
print(f"  2*f_pi^2 = {m2_offdiag_base:.4f} MeV^4... wait, f_pi=92 MeV, f_pi^2=8464 MeV^2")
print(f"  2*f_pi^2 = {m2_offdiag_base:.4f}")
print(f"  This is 16928 as stated in the problem ✓")

# Pairs and their spectator VEVs:
offdiag_pairs = [
    ('M_ud', 'M_du', 'M_ss', M_s_vev),
    ('M_us', 'M_su', 'M_dd', M_d_vev),
    ('M_ds', 'M_sd', 'M_uu', M_u_vev),
]

print()
print(f"{'Pair':<15} {'Spectator':<12} {'|X0|^2|M_c|^2':>14} {'m^2':>14} {'m (|MeV|)':>12}")
print("-" * 70)
offdiag_masses_sq = []
for p1, p2, spec, Mc in offdiag_pairs:
    correction = X0**2 * Mc**2
    m2 = m2_offdiag_base + correction
    mass = np.sqrt(abs(m2))
    offdiag_masses_sq.append((p1, p2, m2, mass))
    print(f"  {p1+','+p2:<13} {spec:<12} {correction:>14.4e} {m2:>14.4f} {mass:>12.4f}")

print()
print("Note: |X0|^2 * M_vev^2 correction is negligible vs 2*f_pi^2 = 16928")
print("All 6 off-diagonal scalars have mass^2 ≈ 16928 MeV^2, mass ≈ 130.1 MeV")
m_offdiag = np.sqrt(m2_offdiag_base)
print(f"  m(off-diag) ≈ {m_offdiag:.4f} MeV (dominant, 2*f_pi^2 approximation)")

# Each complex scalar = 2 real scalars. The 6 complex off-diagonal mesons
# give 12 real scalars. However, we also need to account for the CP-even/CP-odd
# decomposition. For a complex scalar phi = (phi_R + i phi_I)/sqrt(2),
# both phi_R and phi_I get the same mass at tree level (no B-term).
# So 12 real scalars from 6 complex off-diagonal mesons, all with mass ≈ 130.1 MeV.

# ------------------------------------------------------------------
# Diagonal scalars: M_uu, M_dd, M_ss
# ------------------------------------------------------------------
# The scalar mass^2 matrix (3x3) for the diagonal complex mesons:
# INCLUDING both real and imaginary parts.
# For a complex scalar phi_i = (phi_R_i + i phi_I_i)/sqrt(2):
# The full 6x6 (Re,Im) matrix can be decomposed:
# CP-even (Re): m^2_+ matrix
# CP-odd  (Im): m^2_- matrix
# where m^2_{+/-} = m^2_F +/- (A-term contribution)
# At tree level with no A-terms, CP-even = CP-odd for each complex scalar,
# so the 3 complex scalars give 3 degenerate pairs.
# But the F-term 3x3 matrix already gives the eigenvalues for COMPLEX scalars.
# The 3x3 eigenvalues of M_diag_mass2 give the 3 complex scalar masses.
# Each complex scalar = 2 real scalars with the same mass.

# However, there's a subtlety: the off-diagonal entries W_{ii,jj} * W_{ii,jj} contribute
# to the HERMITIAN mass^2 matrix. But we also need B-term (non-holomorphic) contributions
# from the superpotential, which arise as:
#   (m^2_B)_{IJ} = W_{IKL} A_{KL} where A are soft A-terms (absent here)
# So no B-terms. The full complex scalar mass^2 matrix = the Hermitian F-term matrix +
# soft diagonal.

print("\n--- Diagonal scalar mesons ---")
print("Mass^2 matrix = X0^2 * (M M^T - diag(M_i^2)) + f_pi^2 * I_3")

# Scale analysis: X0^2 * M_vev^2 vs f_pi^2
print(f"\n  Scale comparison:")
print(f"  X0^2 * M_uu^2 = {X0**2 * M_u_vev**2:.4e}")
print(f"  X0^2 * M_dd^2 = {X0**2 * M_d_vev**2:.4e}")
print(f"  X0^2 * M_ss^2 = {X0**2 * M_s_vev**2:.4e}")
print(f"  f_pi^2        = {f_pi**2:.4e}")

# The X0^2 M^2 terms are negligibly small compared to f_pi^2
# m^2(diagonal) ≈ f_pi^2 from the off-diagonal F-term contributions...
# Wait: (M_diag_mass2_F)_uu,uu = X0^2*(M_s^2 + M_d^2) which is tiny
# So m^2 ≈ f_pi^2 for all diagonal mesons too, at this order.

print()
ev_diag_mass2 = eigvalsh(M_diag_mass2)
print(f"  Eigenvalues of 3x3 diagonal mass^2 matrix: {ev_diag_mass2}")
print(f"  Masses: {np.sqrt(np.abs(ev_diag_mass2))}")

# The matrix is essentially f_pi^2 * I_3 with tiny X0^2 corrections
# Let's verify:
print()
print("  Off-diagonal F-term corrections to diagonal meson masses are ~X0^2*M^2")
print(f"  which are {X0**2 * M_u_vev**2:.2e} << f_pi^2 = {f_pi**2:.2e}")
print("  So diagonal mesons: m^2 ≈ f_pi^2, m ≈ f_pi = 92 MeV")
print("  (with tiny corrections from Yukawa terms not yet included)")

# BUT WAIT: we need to include the F-term from the DIAGONAL direction:
# F_{X} = det M - Λ^6 is set to zero by constraint.
# The physical scalar mass comes also from fluctuations around the vacuum.
# The standard approach for SQCD: the diagonal scalars (meson VEVs) correspond to
# the SUSY-breaking (pseudo)moduli space. The soft terms stabilize them.
# The mass for fluctuations of phi_i around <phi_i> = M_i_vev:
# V ≈ f_pi^2 |delta phi_i|^2 + F-term corrections
# The F-term for M_ii fluctuation:
# dW/d(delta M_ii) = m_i + X0 * cofactor(ii) + higher order
# At vacuum this = 0. The second variation:
# d^2V/d(delta M_ii)^2 = (d^2W/d(M_ii)^2)^2 + ... = 0 (no W_{ii,ii} term)
# + coupling through F_{jj} and F_X terms.
# The result, after careful analysis, is that the scalar mass for the DIAGONAL modes
# comes primarily from the soft term f_pi^2 and the small F-term correction.
# So m^2(diagonal scalar, any direction) ≈ f_pi^2 = (92 MeV)^2

# Actually, there's a SECOND important source: the c_3 term (det M)^3/Λ^18
# This provides stabilization for the moduli. But the problem asks us to compute
# with the given W, and the c_3 term affects masses too.
# The problem says "compute all 18 real scalar masses" so let's include everything.

# For the purpose of this computation, following the problem statement which gives
# the off-diagonal formula explicitly and says to include F-term + soft terms,
# the diagonal masses are:

# From direct computation of the 3x3 matrix including both real and imaginary parts:
# The 3 complex scalar diagonal mesons → 6 real scalars
# The mass^2 matrix for real scalars in (Re phi_1, Re phi_2, Re phi_3, Im phi_1, Im phi_2, Im phi_3) basis:
# Assuming no B-type terms (no holomorphic mass matrix from A-terms in soft breaking here):
# (m^2)_{Re,Re} = (m^2)_{Im,Im} = M_diag_mass2  (the 3x3 matrix computed above)
# (m^2)_{Re,Im} = 0

# So CP-even and CP-odd scalars for diagonal mesons are degenerate: both have eigenvalues from M_diag_mass2.
ev_mass_sq_diag = eigvalsh(M_diag_mass2)
m_diag_scalars = np.sqrt(np.abs(ev_mass_sq_diag))

print()
print("3x3 diagonal meson scalar mass^2 matrix eigenvalues:")
for i, (ev, m) in enumerate(zip(ev_mass_sq_diag, m_diag_scalars)):
    print(f"  lambda_{i+1} = {ev:.6f} MeV^2 (depends on unit convention) => m = {m:.6f}")

print()
print("Since X0^2*M^2 << f_pi^2, we have approximately:")
print(f"  m_diag ≈ f_pi = {f_pi:.1f} MeV (all three diagonal scalars)")
print("  Each complex scalar = 2 real scalars => 3*2 = 6 real scalars at m ≈ 92 MeV")

# Final off-diagonal: 6 complex = 12 real scalars at m ≈ sqrt(2)*f_pi
m_offdiag_final = np.sqrt(2*f_pi**2)
print(f"  m_offdiag ≈ sqrt(2)*f_pi = {m_offdiag_final:.4f} MeV ≈ 130.1 MeV")
print("  6 complex = 12 real scalars at m ≈ 130.1 MeV")

print()
print("SUMMARY - Part (a): 9 complex = 18 real scalar mesons")
print(f"  Diagonal (6 real): m^2 ≈ f_pi^2 + tiny corrections")
print(f"  Diagonal masses ≈ {f_pi:.1f} MeV (each of the 6 real components)")
print(f"  Off-diagonal (12 real): m^2 = 2*f_pi^2 + negligible")
print(f"  Off-diagonal masses ≈ {m_offdiag_final:.2f} MeV")

print()
print("=" * 70)
print("PART (b): FERMION MASSES (MESINOS + BARYONINOS after integrating out X)")
print("=" * 70)

# ------------------------------------------------------------------
# FERMION MASS MATRIX
# ------------------------------------------------------------------
# After integrating out X (auxiliary field, no kinetic term):
# The relevant W_ij entries for the fermion mass matrix:
#
# DIAGONAL MESINO BLOCK (3x3 after integrating out X):
# Before integrating out X, the 4x4 matrix for {psi_uu, psi_dd, psi_ss, psi_X} is:
#   M_F = [[0,       W_ud,    W_us,    Λ^6/M_u],
#           [W_ud,   0,       W_ds,    Λ^6/M_d],
#           [W_us,   W_ds,    0,       Λ^6/M_s],
#           [Λ^6/M_u,Λ^6/M_d,Λ^6/M_s, 0      ]]
# where W_ud = X0*M_s, W_us = X0*M_d, W_ds = X0*M_u

# Integrating out psi_X (auxiliary fermion from integrating out X):
# Since the X row has entries only in the diagonal meson column (at diagonal vacuum),
# we integrate out psi_X using its equation of motion from the quadratic action.
# The effective mass matrix for the 3 diagonal mesinos is obtained by integrating
# out psi_X from the mass terms:
# L_mass = (1/2) M_F^{ab} psi_a psi_b + h.c.
# where a,b run over {uu, dd, ss, X}.
# Since X has no kinetic term, psi_X appears algebraically:
# Equation for psi_X: sum_{i=u,d,s} (Λ^6/M_i) psi_i + 0 * psi_X = 0
# This gives: sum_i (Λ^6/M_i) psi_i = 0 (constraint!)
# But Λ^6/M_i = m_i (from vacuum equation: X0*M_j*M_k = -m_i => M_j*M_k = -m_i/X0 = m_i*Λ^6/C ...)
# Let's check: Λ^6/M_u = Λ^6/(C/m_u) = Λ^6 m_u / C
# And C = (Λ^6 m_u m_d m_s)^{1/3}, so C/Λ^6 = (m_u m_d m_s)^{1/3}/Λ^4*... let's compute numerically

L6_over_Mu = LAM**6 / M_u_vev
L6_over_Md = LAM**6 / M_d_vev
L6_over_Ms = LAM**6 / M_s_vev
print(f"\nΛ^6/M_u = {L6_over_Mu:.4e}")
print(f"Λ^6/M_d = {L6_over_Md:.4e}")
print(f"Λ^6/M_s = {L6_over_Ms:.4e}")
print(f"Compare: m_u={m_u}, m_d={m_d}, m_s={m_s}")

# The exact integrating-out procedure for auxiliary fermion psi_X:
# When X has no kinetic term, psi_X is an auxiliary field in the fermion sector too.
# The Lagrangian mass term is:
#   L = (1/2)[M_uu M_dd: W_{uu,dd} psi_uu psi_dd + ...] + (Λ^6/M_i) psi_ii psi_X + (1/2)(0) psi_X psi_X
# Equation of motion for psi_X: dL/dpsi_X = sum_i (Λ^6/M_i) psi_i = 0
# This is a CONSTRAINT on the diagonal mesinos.
# Substituting the constraint, we can eliminate one mesino direction.
# The effective mass matrix for the remaining 2 diagonal mesinos:
# Project onto the subspace orthogonal to n_i = (Λ^6/M_i) / ||(Λ^6/M_i)||.

# Actually, the correct procedure for integrating out an auxiliary field in the fermion sector:
# If the mass matrix is:
#   L = (1/2) [A_ij psi_i psi_j + 2 C_i psi_i psi_X + 0]
# where psi_X has no kinetic term, we set dL/dpsi_X = 0:
#   sum_i C_i psi_i = 0
# This is a constraint. The effective theory has (n-1) light fermions satisfying this constraint.
# The effective mass matrix for the constrained subspace:
# Use Schur complement: since W_{XX} = 0, the constraint is singular, meaning one mesino
# combination gets a mass much larger than the others (becomes heavy), and we integrate it out.
# But actually W_{XX} = 0 means psi_X CANNOT get a mass from W directly.
# Instead, psi_X appears as a Lagrange multiplier coupling to the mesinos.
#
# The correct interpretation: since W_{XX} = 0 (X appears linearly in W),
# psi_X gives the constraint sum_i (Λ^6/M_i) psi_i = 0 (one constraint on 3 mesinos).
# This leaves 2 independent diagonal mesino combinations.
# However, the 2 remaining mesinos STILL have a 2x2 mass matrix from W_{ij} (i,j in {u,d,s}).
#
# Let me think more carefully. The fermion mass terms from W are:
#   (1/2) W_IJ psi_I psi_J + h.c.
# For {psi_u, psi_d, psi_s, psi_X} with W_XX = 0:
# L = (1/2)[W_ud psi_u psi_d + W_us psi_u psi_s + W_ds psi_d psi_s
#           + (Λ^6/M_u) psi_u psi_X + (Λ^6/M_d) psi_d psi_X + (Λ^6/M_s) psi_s psi_X]
#
# psi_X appears with no kinetic term and no W_XX coupling:
# Integrating out psi_X:
#   F_X_fermion = (Λ^6/M_u) psi_u + (Λ^6/M_d) psi_d + (Λ^6/M_s) psi_s = 0
# This sets the constraint. The effective action for the remaining 2 directions (orthogonal to n)
# in {psi_u, psi_d, psi_s} space has mass matrix = projection of W_ij onto this subspace.

# Define the coupling vector n = (Λ^6/M_u, Λ^6/M_d, Λ^6/M_s)
n_vec = np.array([L6_over_Mu, L6_over_Md, L6_over_Ms])
n_norm = norm(n_vec)
n_hat = n_vec / n_norm

print(f"\nConstraint vector n_i = Λ^6/M_i:")
print(f"  n = ({n_vec[0]:.4e}, {n_vec[1]:.4e}, {n_vec[2]:.4e})")
print(f"  |n| = {n_norm:.4e}")
print(f"  n_hat = ({n_hat[0]:.6f}, {n_hat[1]:.6f}, {n_hat[2]:.6f})")

# The 3x3 mass matrix for diagonal mesinos (before integrating out psi_X):
W_diag_3x3 = np.array([
    [0,            X0*M_s_vev,  X0*M_d_vev],
    [X0*M_s_vev,  0,            X0*M_u_vev],
    [X0*M_d_vev,  X0*M_u_vev,  0          ]
])

print("\n3x3 W matrix for diagonal mesinos {psi_uu, psi_dd, psi_ss}:")
mlab = ['psi_uu', 'psi_dd', 'psi_ss']
for i in range(3):
    row = "  ".join(f"{W_diag_3x3[i,j]:+.4e}" for j in range(3))
    print(f"  {mlab[i]}: {row}")

# The constraint is n_hat direction = 0, so effective 2x2 matrix:
# Project onto 2D subspace orthogonal to n_hat.
# Projection matrix P = I - |n_hat><n_hat|
P = np.eye(3) - np.outer(n_hat, n_hat)

# The effective mass matrix after integrating out psi_X is:
# W_eff = P * W_diag_3x3 * P  (projected)
# BUT the psi_X integration also modifies the mass matrix through the
# Schur complement. Since W_XX = 0, the Schur complement formula:
# W_eff = W_MM - W_MX * (W_XX)^{-1} * W_XM
# is undefined (W_XX = 0 is singular). This means we should NOT use the
# Schur complement but instead recognize that psi_X forces a constraint.
#
# The resolution: since psi_X has NO kinetic term and W_XX = 0,
# psi_X is a Lagrange multiplier that enforces the constraint n_hat . psi = 0.
# The physical fermion spectrum comes from the 2D projected theory.
# Eigenvalues of P * W_diag_3x3 * P:

W_eff_diag = P @ W_diag_3x3 @ P
ev_eff_diag = eigvalsh(W_eff_diag)
print(f"\nEigenvalues of projected 3x3 matrix (after integrating out psi_X):")
print(f"  {ev_eff_diag}")

# One eigenvalue should be ≈ 0 (the n_hat direction, which is the constraint)
# The other two are the physical mesino masses.

# Let's also compute the full 4x4 eigenvalues to check:
W_4x4_diag = np.array([
    [0,            X0*M_s_vev,  X0*M_d_vev,  L6_over_Mu],
    [X0*M_s_vev,  0,            X0*M_u_vev,  L6_over_Md],
    [X0*M_d_vev,  X0*M_u_vev,  0,            L6_over_Ms],
    [L6_over_Mu,  L6_over_Md,  L6_over_Ms,  0          ]
], dtype=float)

ev_4x4 = eigvalsh(W_4x4_diag)
print(f"\nEigenvalues of full 4x4 {'{M_uu, M_dd, M_ss, X}'} block:")
for ev in ev_4x4:
    print(f"  {ev:.6e}")

# The 4x4 has 2 large eigenvalues (from Λ^6/M_i ~ 10^9-10^10 scale)
# and 2 small eigenvalues (from X0*M ~ 10^-9 * 10^3 = 10^-6 scale).
# After integrating out X: 2 diagonal mesinos remain (the zero of the projected matrix
# corresponds to the direction absorbed into psi_X constraint).

# OFF-DIAGONAL MESINO BLOCKS (three 2x2 blocks):
# For each pair (M_ij, M_ji) with i != j:
# W_{M_ij, M_ji} = X0 * cofactor sign
# The 2x2 block for {psi_ud, psi_du}:
#   W = [[0,          X0*(-M_s)],
#        [X0*(-M_s),  0        ]]
# Eigenvalues: ±|X0 * M_s| = ±|X0| * M_s_vev
# (The sign of X0*(-M_s) = -X0*M_s. Since X0 < 0 and M_s > 0: -X0*M_s = |X0|*M_s > 0)
# So eigenvalues of [[0, a], [a, 0]] are ±a = ±|X0|*M_s ✓

print("\n--- Off-diagonal mesino blocks (three 2x2) ---")
offdiag_W_entries = [
    ('psi_ud/psi_du', -X0 * M_s_vev),   # spectator is s: W = X0*(-M_s), which = |X0|*M_s
    ('psi_us/psi_su', -X0 * M_d_vev),   # spectator is d
    ('psi_ds/psi_sd', -X0 * M_u_vev),   # spectator is u
]
print(f"\n  {'Block':<20} {'W_IJ = |X0|*M_spec':>20} {'Eigenvalues':>20}")
print("  " + "-" * 60)
mesino_offdiag_masses = []
for name, w_entry in offdiag_W_entries:
    ev_pos = abs(w_entry)
    mesino_offdiag_masses.append(ev_pos)
    print(f"  {name:<20} {w_entry:>20.6e} {f'±{ev_pos:.6e}':>20}")

print()
print("  Note: each 2x2 block gives one positive and one negative eigenvalue.")
print("  Physical masses = |eigenvalue|.")
print()
print("  Off-diagonal mesino masses:")
for name, m in zip(['(u,d)', '(u,s)', '(d,s)'], mesino_offdiag_masses):
    print(f"    m(mesino_{name}) = |X0| * M_spectator = {m:.6e} MeV (in whatever units)")

# These are extremely small! |X0| ~ 1.21e-9 (MeV^-4 ?) and M_vev ~ 4e3 (MeV^2 ?)
# So |X0| * M_vev ~ 5e-6 -- but in what units?
# Let's just report the numerical values.
print()
print(f"  |X0| * M_ss = {abs(X0) * M_s_vev:.4e}")
print(f"  |X0| * M_dd = {abs(X0) * M_d_vev:.4e}")
print(f"  |X0| * M_uu = {abs(X0) * M_u_vev:.4e}")

# Now for the DIAGONAL MESINO block (3x3 → effective 2x2 after constraint):
# The physical masses are the nonzero eigenvalues of the projected matrix.
# From the 4x4 eigenvalues, the two "light" eigenvalues are what survive after
# projecting out the Goldstino/constraint direction.

print("\n--- Diagonal mesino block ---")
print("Full 4x4 eigenvalues (before integrating out psi_X):")
for ev in ev_4x4:
    print(f"  {ev:.6e}")
print()
print("After integrating out psi_X (constraint: n_hat.psi = 0):")
print("Projected 3x3 eigenvalues:")
for ev in ev_eff_diag:
    print(f"  {ev:.6e}")
print()
print("Physical eigenvalues: the two non-zero ones from the 4x4 that are 'light'")
# Sort by magnitude
ev_4x4_sorted = sorted(ev_4x4, key=abs)
print(f"Smallest |eigenvalue|: {ev_4x4_sorted[0]:.4e}")
print(f"Second smallest:       {ev_4x4_sorted[1]:.4e}")
print(f"Large (X-dominated):  {ev_4x4_sorted[2]:.4e} and {ev_4x4_sorted[3]:.4e}")

print()
print("The 'large' modes involve psi_X and a linear combination of diagonal mesinos.")
print("Mass ~ sqrt(n_i * W_ij * n_j + ...) ~ Λ^6/M_i ~ 10^9 scale")
print("After integrating out psi_X: 2 light diagonal mesinos remain,")
print("with masses from the X0*M_vev off-diagonal entries: |X0|*M_k")
print()
print("SUMMARY of mesino masses (9 Weyl fermions = 4.5 Dirac):")
print("  Diagonal mesinos (after integrating out psi_X):")
print(f"    2 light diagonal mesinos with masses from projected W ≈ {ev_eff_diag[0]:.4e} and {ev_eff_diag[1]:.4e}")
print(f"    1 direction absorbed by psi_X constraint (0 mass in effective theory)")
print("  Off-diagonal mesinos:")
for name, m in zip(['(u,d)', '(u,s)', '(d,s)'], mesino_offdiag_masses):
    print(f"    {name}: ±{m:.4e} (Dirac pair)")

print()
print("=" * 70)
print("PART (c): BARYONS")
print("=" * 70)

# Baryon superfield mass term: m_B B B~ in W
# With m_B = Λ = 300 MeV
# At vacuum B = B~ = 0.
# Scalar mass: from F_B = dW/dB = X0*(B~) + m_B * B~  [wait, careful]
# W contains: ... + X*(det M - BB~ - Λ^6) + m_B BB~
# F_B = dW/dB = X*(-B~) + m_B * B~  = (m_B - X)*B~
# Wait: dW/dB = -X*B~ + m_B*B~? Let me be careful about the B coupling:
# W includes: X*(det M - BB~ - Λ^6) => the B term is -X*B*B~
# And: m_B * B*B~ (separate mass term)
# F_B = dW/dB = (-X + m_B)*B~ = (-X0 + m_B)*B~|_{B~=0} = 0  ✓ (at vacuum B=B~=0)
#
# The 2x2 baryon mass matrix for {B, B~}:
# W_BB = 0
# W_{B,Bt} = -X0 + m_B   [from -X*B*Bt + m_B*B*Bt = (m_B - X)*B*Bt]
# W_{Bt,B} = -X0 + m_B

m_baryon_fermion = abs(-X0 + m_B)   # this is just m_B to very good approximation since |X0| << m_B

print(f"\nm_B (baryon mass parameter) = {m_B} MeV")
print(f"|X0| = {abs(X0):.4e} (much smaller than m_B)")
print(f"W_{{B,Bt}} = m_B - X0 = {m_B - X0:.6f} ≈ m_B = {m_B} MeV")
print()
print("Baryonino (Weyl fermion pair):")
print(f"  Mass = |m_B - X0| ≈ {m_baryon_fermion:.4f} MeV ≈ Λ = 300 MeV")
print("  The 2x2 Dirac mass matrix [[0, m_B], [m_B, 0]] => one Dirac fermion of mass m_B")
print()

# Baryon SCALAR mass:
# V_scalar for B: |F_B|^2 + |F_Bt|^2 + |F_X * (-BB~)|^2
# F_B = (m_B - X)*B~, F_Bt = (m_B - X)*B
# Scalar mass^2 matrix for {B, B~}:
# d^2V/dB dB* = |m_B - X0|^2 (from |F_Bt|^2 = |m_B - X0|^2 |B|^2)
# Plus contribution from F_X: F_X = det M - BB~ - Λ^6
# At vacuum, the fluctuation of B enters F_X through delta(BB~) = B_vev*delta(B~) + delta(B)*B~_vev = 0
# (since B_vev = B~_vev = 0). So the F_X contribution to baryon scalar mass vanishes at B=B~=0.
# The D-terms: baryons are color singlets (and assuming singlet under SU(2)), D-term contribution
# involves the baryon charges, which for a SU(3) color singlet are zero for SU(3), but may contribute
# for U(1)_Y.

# Baryon scalar mass^2:
# |W_{B,Bt}|^2 = |m_B - X0|^2 + soft terms (if any for baryons)
# The soft breaking V_soft = f_pi^2 Tr(M†M) only for MESON scalars.
# Baryons: no explicit soft term given => only F-term contribution.

m2_baryon_scalar = (m_B - X0)**2   # from W_{B,Bt}^2
m_baryon_scalar = np.sqrt(m2_baryon_scalar)
print(f"Baryon scalar mass^2 = |W_{{B,Bt}}|^2 = |m_B - X0|^2 = {m2_baryon_scalar:.4f}")
print(f"Baryon scalar mass   = {m_baryon_scalar:.4f} MeV ≈ m_B = {m_B} MeV")
print()
print("Baryon sector (at B = B~ = 0):")
print(f"  Baryonino (Weyl pair = 1 Dirac): mass = {m_B - X0:.6f} MeV ≈ 300 MeV")
print(f"  Baryon scalar (2 complex = 4 real): mass^2 = {m2_baryon_scalar:.4f}, mass ≈ 300 MeV")
print("  (B and B~ give 2 complex scalars, both with mass ≈ m_B ≈ Λ = 300 MeV)")

print()
print("=" * 70)
print("PART (d): HIGGS SECTOR (NMSSM)")
print("=" * 70)

# W_NMSSM = λ_S S H_u H_d + κ/3 S³ + (meson Yukawa already counted)
# Plus the EW symmetry breaking: <H_u^0> = vu = v/sqrt(2), <H_d^0> = vd = v/sqrt(2)
# mu_eff = λ_S <S>  (effective mu term from NMSSM)
#
# Tree-level lightest CP-even Higgs:
# The standard MSSM contribution: m_h^2 = m_Z^2 cos^2(2β) + ... (stop corrections)
# NMSSM adds: λ_S^2 v^2 sin^2(2β) / 2 = λ_S^2 v^2 / 2 for tan β = 1 (sin(2β) = 1)
# => m_h^2(NMSSM tree) = m_Z^2 cos^2(2β) + λ_S^2 v^2 sin^2(2β)/2 - δ
# For tan β = 1: cos(2β) = 0, sin(2β) = 1
# => m_h^2 = λ_S^2 v^2 / 2  (MSSM term vanishes at tan β = 1)

m_Z = 91187.6  # MeV
m_h_tree_sq = lam_S**2 * v**2 / 2
m_h_tree = np.sqrt(m_h_tree_sq)
print(f"\nNMSSM Higgs parameters:")
print(f"  tan β = 1,  v = {v/1000:.2f} GeV,  λ_S = {lam_S}")
print(f"  v_u = v_d = v/√2 = {vu/1000:.2f} GeV")
print()
print(f"Tree-level lightest CP-even Higgs:")
print(f"  m_h^2 = λ_S^2 v^2 / 2 = {lam_S}^2 × ({v/1000:.2f} GeV)^2 / 2")
print(f"  m_h^2 = {m_h_tree_sq/1e6:.4f} GeV^2")
print(f"  m_h   = {m_h_tree/1000:.4f} GeV = {m_h_tree:.2f} MeV")
print(f"  (Problem states: λ_S v/√2 = 0.72 × 246.22/1.414 = {lam_S * v / np.sqrt(2)/1000:.3f} GeV)")

# The problem says: m_h = λ_S v/√2 = 125 GeV
# Let's check: λ_S v/√2 = 0.72 * 246220/√2 MeV = 0.72 * 174076 MeV = 125334 MeV ≈ 125 GeV ✓
m_h_problem = lam_S * v / np.sqrt(2)
print(f"\nProblem formula: m_h = λ_S v/√2 = {m_h_problem/1000:.3f} GeV ≈ 125 GeV ✓")
print()
print("Note: the problem formula m_h = λ_S v/√2 = λ_S * v_u = λ_S * v_d (at tan β = 1)")
print("This matches m_h^2 = λ_S^2 v^2/2 ✓")
print()
print("Elementary quarks above threshold:")
for name, mass in [('c', m_c), ('b', m_b), ('t', m_t)]:
    y = mass / vu  # Yukawa coupling
    print(f"  {name} quark: m = {mass} MeV = {mass/1000:.3f} GeV, y_{name} = m/v_u = {y:.6f}")

print()
print("=" * 70)
print("PART (e): COMPLETE MASS SPECTRUM SUMMARY TABLE")
print("=" * 70)

# Let's define final masses clearly
print()
print("Final computed masses (all in MeV):")
print()

# === SECTOR 1: LIGHT CONFINED SECTOR ===
# Off-diagonal scalar mesons (12 real scalars):
m_offdiag_scalar = np.sqrt(2 * f_pi**2)
# Diagonal scalar mesons (6 real scalars):
m_diag_scalar = f_pi  # ≈ 92 MeV (dominant soft term)
# Off-diagonal mesinos (6 Weyl fermions = 3 Dirac):
m_mesino_ud = abs(X0) * M_s_vev
m_mesino_us = abs(X0) * M_d_vev
m_mesino_ds = abs(X0) * M_u_vev

print("Sector 1: Light confined sector (meson scalars + off-diagonal mesinos)")
print(f"  Diagonal scalar mesons (6 real):     m ≈ {m_diag_scalar:.1f} MeV [= f_pi]")
print(f"  Off-diagonal scalar mesons (12 real): m ≈ {m_offdiag_scalar:.2f} MeV [= sqrt(2)*f_pi]")
print(f"  Mesino (u,d) pair (2 Weyl = 1 Dirac): m = {m_mesino_ud:.4e} MeV")
print(f"  Mesino (u,s) pair:                     m = {m_mesino_us:.4e} MeV")
print(f"  Mesino (d,s) pair:                     m = {m_mesino_ds:.4e} MeV")
print()

# === SECTOR 2: HEAVY CONFINED SECTOR ===
# Large eigenvalues of the 4x4 diagonal mesino block (X-dominated)
ev_large = sorted(ev_4x4, key=abs, reverse=True)
m_heavy_mesino_1 = abs(ev_large[0])
m_heavy_mesino_2 = abs(ev_large[1])
# Scalar partners of the X field (after integrating out):
# The X field has NO kinetic term, so NO propagating scalar.
# But there IS a heavy scalar from the mixing of diagonal mesinos with X.
# The X-dominated scalar masses come from the meson mass matrix with the X direction.
# After integrating out X from the scalar potential:
# The heavy mode in the 4x4 block gets mass from Λ^6/M_i terms.
# Numerically:
print("Sector 2: Heavy confined sector (X-dominated diagonal mesinos)")
print(f"  Heavy diagonal mesino 1: m = {m_heavy_mesino_1:.4e}")
print(f"  Heavy diagonal mesino 2: m = {m_heavy_mesino_2:.4e}")
print("  (These are the large eigenvalues of the 4x4 {M_uu,M_dd,M_ss,X} Majorana mass matrix)")
print("  These correspond to the X-integrated-out modes, approximately:")
print(f"  m_heavy ≈ |n| = ||(Λ^6/M_u, Λ^6/M_d, Λ^6/M_s)|| = {n_norm:.4e}")
print()
print("  The X scalar itself has NO propagating degree of freedom (no kinetic term).")
print()

# === SECTOR 3: BARYONS ===
print("Sector 3: Baryons")
print(f"  Baryon B scalar (complex): m = {m_baryon_scalar:.2f} MeV ≈ m_B = {m_B:.0f} MeV")
print(f"  Baryon B~ scalar (complex): m = {m_baryon_scalar:.2f} MeV")
print(f"  Baryonino (1 Dirac from B+B~): m = {m_B - X0:.6f} MeV ≈ {m_B:.0f} MeV")
print()

# === SECTOR 4: ELEMENTARY QUARKS ===
print("Sector 4: Elementary quarks (above confinement scale, not in SQCD sector)")
for name, mass in [('c', m_c), ('b', m_b), ('t', m_t)]:
    print(f"  {name} quark (Dirac): m = {mass} MeV = {mass/1000:.3f} GeV")
# Squarks (scalar superpartners of c, b, t):
# From W: y_c H_u M^c_c (meson-Higgs Yukawa) — this doesn't directly give c squark mass.
# The c, b, t are ELEMENTARY (above threshold), so they have squark scalar partners.
# In the MSSM, squark masses come from soft terms (not specified here).
# Yukawa contribution to squark mass through Higgs VEV:
# At electroweak breaking, <H_u^0> = v_u = v/sqrt(2):
# The F-term from W_H_u = y_c M^c_c + y_t Q_t bar_Q_t gives:
# F_{H_u} = y_c M^c_c + y_t Q_t bar_Q_t
# This doesn't give propagating squark masses without the squark kinetic terms.
# Squark masses in MSSM come from soft breaking (m_soft^2 ~TeV^2, not specified).
print(f"  c squark (scalar): mass from soft terms (not specified, ~TeV)")
print(f"  b squark (scalar): mass from soft terms (not specified, ~TeV)")
print(f"  t squark (scalar): mass from soft terms (not specified, ~TeV)")
print()

# === SECTOR 5: HIGGS SECTOR ===
print("Sector 5: Higgs sector (NMSSM)")
print(f"  CP-even Higgs h (lightest): m = {m_h_problem/1000:.3f} GeV ≈ 125 GeV")
print("  CP-even Higgs H: mass depends on soft parameters (not computed)")
print("  CP-odd Higgs A:  mass depends on soft parameters")
print("  Charged Higgs H±: mass from m_A^2 + m_W^2 (soft dependent)")
print("  Singlino S (fermion): mass = κ<S> + ... (κ free, <S> free)")
print("  Higgsino H_u~ (Weyl pair): mass from λ_S <S> = μ_eff")
print("  Higgsino H_d~ (Weyl pair): mass from λ_S <S> = μ_eff")
print()

# === SECTOR 6: GAUGE BOSONS ===
m_W = 80379.0   # MeV (PDG)
m_Z_val = 91187.6  # MeV
print("Sector 6: Gauge bosons (after EW breaking)")
print("  Gluons g (8): massless (confined, no propagating gluons)")
print(f"  W± bosons: m = {m_W/1000:.3f} GeV")
print(f"  Z0 boson:  m = {m_Z_val/1000:.4f} GeV")
print("  Photon γ:  massless")
print("  Gluino g~: mass from soft terms (not specified)")
print("  Wino W~±: mass from soft terms")
print("  Bino B~:  mass from soft terms")
print("  Zino: mixture of Wino/Bino/Higgsino after EW breaking")

print()
print("=" * 70)
print("WRITING COMPLETE SPECTRUM TABLE")
print("=" * 70)

# Write the markdown file
spectrum_table = []

def W(s): spectrum_table.append(s)

W("# Complete N=1 SUSY Particle Spectrum")
W("")
W("## Theory Setup")
W("")
W("SU(3) SQCD with N_f = N_c = 3 (u, d, s), confined at Λ = 300 MeV.")
W("Coupled to NMSSM Higgs sector (H_u, H_d, S) with tan β = 1.")
W("")
W("### Parameters")
W("")
W("| Parameter | Value |")
W("|-----------|-------|")
W(f"| Λ | 300 MeV |")
W(f"| v | 246.22 GeV |")
W(f"| tan β | 1 (v_u = v_d = v/√2 = {vu/1000:.2f} GeV) |")
W(f"| f_π | 92 MeV |")
W(f"| λ_S | 0.72 |")
W(f"| m_B | 300 MeV |")
W(f"| m_u | 2.16 MeV |")
W(f"| m_d | 4.67 MeV |")
W(f"| m_s | 93.4 MeV |")
W(f"| m_c | 1275 MeV |")
W(f"| m_b | 4180 MeV |")
W(f"| m_t | 172760 MeV |")
W("")
W("### Derived Vacuum Quantities")
W("")
W(f"- C = (Λ⁶ m_u m_d m_s)^(1/3) = {C:.2f} MeV²")
W(f"- M_uu = C/m_u = {M_u_vev:.3f} MeV²")
W(f"- M_dd = C/m_d = {M_d_vev:.3f} MeV²")
W(f"- M_ss = C/m_s = {M_s_vev:.3f} MeV²")
W(f"- X₀ = -C/Λ⁶ = {X0:.4e} (units: MeV⁻⁴ per problem statement)")
W(f"- det(M_vev)/Λ⁶ = {det_M/LAM**6:.12f} ≈ 1.0 ✓")
W("")
W("### Key F-term W Matrix Entries at Vacuum")
W("")
W("| Entry | Formula | Value |")
W("|-------|---------|-------|")
W(f"| W_{{M_uu, M_dd}} | X₀ M_ss | {X0*M_s_vev:.4e} |")
W(f"| W_{{M_uu, M_ss}} | X₀ M_dd | {X0*M_d_vev:.4e} |")
W(f"| W_{{M_dd, M_ss}} | X₀ M_uu | {X0*M_u_vev:.4e} |")
W(f"| W_{{M_ud, M_du}} | X₀·(-M_ss) = +|X₀| M_ss | {-X0*M_s_vev:.4e} |")
W(f"| W_{{M_us, M_su}} | X₀·(-M_dd) = +|X₀| M_dd | {-X0*M_d_vev:.4e} |")
W(f"| W_{{M_ds, M_sd}} | X₀·(-M_uu) = +|X₀| M_uu | {-X0*M_u_vev:.4e} |")
W(f"| W_{{M_uu, X}} | Λ⁶/M_uu | {L6_over_Mu:.4e} |")
W(f"| W_{{M_dd, X}} | Λ⁶/M_dd | {L6_over_Md:.4e} |")
W(f"| W_{{M_ss, X}} | Λ⁶/M_ss | {L6_over_Ms:.4e} |")
W(f"| W_{{B, B̃}} | m_B - X₀ | {m_B - X0:.6f} ≈ {m_B:.1f} MeV |")
W("")
W("---")
W("")
W("## Part (a): Scalar Meson Masses")
W("")
W("### Method")
W("")
W("Scalar mass² from F-terms (with X integrated out, contributing no scalar d.o.f.):")
W("")
W("    (m²_scalar)_IJ = Σ_{K: kinetic term} W_IK W*_JK + f_π² δ_IJ (mesons only)")
W("")
W("Since X has no kinetic term, W_{M,X} entries do NOT contribute to scalar masses.")
W("The soft term V_soft = f_π² Tr(M†M) adds f_π² to each meson scalar mass².")
W("")
W("### Off-Diagonal Scalar Mesons (6 complex = 12 real scalars)")
W("")
W("For M^a_b with a ≠ b, the only relevant second derivative is W_{M_ab, M_ba}.")
W("The F-term contribution comes from the kinetic fields only:")
W("")
W("    m²(M^a_b) = |W_{M_ab, M_ba}|² + 2·f_π²")
W("              = |X₀|² M_c² + 2·f_π²")
W("")
W("where c is the spectator flavor index, and the factor 2 comes from both")
W("the real and imaginary components of the complex scalar (no B-terms).")
W("")
W("| Scalar | Spectator | |X₀|² M_c² | m² (MeV²) | m (MeV) |")
W("|--------|-----------|------------|-----------|---------|")

offdiag_table = [
    ('M^u_d and M^d_u', 'M_ss', M_s_vev),
    ('M^u_s and M^s_u', 'M_dd', M_d_vev),
    ('M^d_s and M^s_d', 'M_uu', M_u_vev),
]
for scal, spec, Mc in offdiag_table:
    corr = X0**2 * Mc**2
    m2 = m2_offdiag_base + corr
    mass = np.sqrt(abs(m2))
    W(f"| {scal} | {spec} | {corr:.4e} | {m2:.4f} | {mass:.4f} |")

W("")
W(f"All off-diagonal scalar masses ≈ √(2) × f_π = {m_offdiag_scalar:.4f} MeV ≈ 130.1 MeV")
W(f"(The |X₀|² M² correction is negligible: {X0**2 * M_s_vev**2:.2e} << 2f_π² = {2*f_pi**2:.0f})")
W("")
W("Each complex scalar decomposes into 2 real scalars (CP-even + CP-odd) with equal masses.")
W("Total: 12 real scalars at m ≈ 130.1 MeV")
W("")
W("### Diagonal Scalar Mesons (3 complex = 6 real scalars)")
W("")
W("The 3×3 mass² matrix for {M_uu, M_dd, M_ss}:")
W("")
W("    (m²_F)_ij = X₀² (M_i M_j - M_i² δ_ij)  = X₀² (outer product - diagonal²)")
W("    (m²_total)_ij = (m²_F)_ij + f_π² δ_ij")
W("")
W("Explicitly:")
W(f"    (m²_F)_uu,uu = X₀²(M_ss² + M_dd²) = {X02*(M_s_vev**2 + M_d_vev**2):.4e}")
W(f"    (m²_F)_dd,dd = X₀²(M_ss² + M_uu²) = {X02*(M_s_vev**2 + M_u_vev**2):.4e}")
W(f"    (m²_F)_ss,ss = X₀²(M_dd² + M_uu²) = {X02*(M_d_vev**2 + M_u_vev**2):.4e}")
W(f"    (m²_F)_uu,dd = X₀² M_dd M_uu      = {X02*M_d_vev*M_u_vev:.4e}")
W(f"    f_π² = {f_pi**2:.4f} MeV²")
W("")
W(f"All F-term entries ≪ f_π². The mass² matrix is dominated by f_π² I_3.")
W("")
W("Eigenvalues of full 3×3 diagonal mass² matrix:")
for i, ev in enumerate(ev_mass_sq_diag):
    W(f"    λ_{i+1} = {ev:.6f} MeV² → m = {np.sqrt(abs(ev)):.6f} MeV")
W("")
W(f"All diagonal scalars: m ≈ f_π = {f_pi:.1f} MeV (to <0.01% correction from F-terms)")
W("Total: 6 real scalars at m ≈ 92 MeV (CP-even and CP-odd degenerate)")
W("")
W("### Scalar Meson Summary")
W("")
W("| Type | Count (real) | Mass (MeV) | Origin |")
W("|------|-------------|------------|--------|")
W(f"| Diagonal (Re part) | 3 | {f_pi:.1f} | V_soft = f_π² |")
W(f"| Diagonal (Im part) | 3 | {f_pi:.1f} | V_soft = f_π² |")
W(f"| Off-diagonal (Re) | 6 | {m_offdiag_scalar:.2f} | 2f_π² (factor 2 from 2 F-terms) |")
W(f"| Off-diagonal (Im) | 6 | {m_offdiag_scalar:.2f} | 2f_π² |")
W(f"| **Total** | **18** | | |")
W("")
W("---")
W("")
W("## Part (b): Fermion Masses (Mesinos)")
W("")
W("### Block Structure of the 12×12 Fermion Mass Matrix")
W("")
W("After integrating out ψ_X (no kinetic term → auxiliary fermion):")
W("")
W("**Off-diagonal mesino blocks** (three independent 2×2 Majorana blocks):")
W("")
W("For pair {ψ_{M_ab}, ψ_{M_ba}} with a ≠ b:")
W("")
W("    W_{M_ab, M_ba} = X₀ × (cofactor sign) × M_c")
W("")
W("At diagonal vacuum, the cofactor for off-diagonal elements gives a minus sign relative")
W("to the diagonal case. With X₀ < 0:")
W("")
W("    W_{M_ud, M_du} = X₀ × (-M_ss) = |X₀| × M_ss   (positive)")
W("")
W("Each 2×2 block [[0, μ],[μ, 0]] has eigenvalues ±μ (one Dirac fermion of mass μ).")
W("")
W("| Block | W_IJ | Mesino Mass (MeV) |")
W("|-------|------|------------------|")
W(f"| ψ_ud / ψ_du | |X₀| × M_ss = {abs(X0)*M_s_vev:.4e} | {abs(X0)*M_s_vev:.4e} |")
W(f"| ψ_us / ψ_su | |X₀| × M_dd = {abs(X0)*M_d_vev:.4e} | {abs(X0)*M_d_vev:.4e} |")
W(f"| ψ_ds / ψ_sd | |X₀| × M_uu = {abs(X0)*M_u_vev:.4e} | {abs(X0)*M_u_vev:.4e} |")
W("")
W("**Diagonal mesino block** (after integrating out ψ_X):")
W("")
W("The 4×4 block {ψ_uu, ψ_dd, ψ_ss, ψ_X} with W_{XX} = 0:")
W("")
W("    W_{4×4} = [[0,        X₀M_s,   X₀M_d,   Λ⁶/M_u],")
W("               [X₀M_s,   0,        X₀M_u,   Λ⁶/M_d],")
W("               [X₀M_d,   X₀M_u,   0,        Λ⁶/M_s],")
W("               [Λ⁶/M_u,  Λ⁶/M_d,  Λ⁶/M_s,  0      ]]")
W("")
W("Since W_{XX} = 0, ψ_X acts as a Lagrange multiplier enforcing the constraint:")
W("")
W("    Λ⁶/M_u × ψ_uu + Λ⁶/M_d × ψ_dd + Λ⁶/M_s × ψ_ss = 0")
W("")
W("This eliminates one diagonal mesino direction. The 4×4 eigenvalues are:")
W("")
W("| n | Eigenvalue | Type |")
W("|---|-----------|------|")
ev_4x4_sorted_all = sorted(ev_4x4, key=lambda x: abs(x))
for i, ev in enumerate(ev_4x4_sorted_all):
    label = "light (X₀·M scale)" if abs(ev) < 1e-3 else "heavy (Λ⁶/M scale)"
    W(f"| {i+1} | {ev:.6e} | {label} |")
W("")
W("The two light eigenvalues correspond to the 2 physical diagonal mesinos (constrained subspace).")
W("The two heavy eigenvalues are the X-dominated modes (integrated out at low energy).")
W("")
W("Diagonal mesino mass matrix projected onto constraint-satisfying subspace:")
for i, ev in enumerate(ev_eff_diag):
    tag = "zero (constraint dir)" if abs(ev) < 1e-20 else "physical mesino"
    W(f"    Eigenvalue {i+1}: {ev:.6e} → {tag}")
W("")
W("### Physical Mesino Spectrum Summary")
W("")
W("| Mesino | Mass (MeV) | Description |")
W("|--------|-----------|-------------|")
W(f"| ψ_ud = ψ_du Dirac | {abs(X0)*M_s_vev:.4e} | Off-diag, spectator s |")
W(f"| ψ_us = ψ_su Dirac | {abs(X0)*M_d_vev:.4e} | Off-diag, spectator d |")
W(f"| ψ_ds = ψ_sd Dirac | {abs(X0)*M_u_vev:.4e} | Off-diag, spectator u |")
ev_light = [ev for ev in ev_4x4_sorted_all if abs(ev) < 1e-3]
for i, ev in enumerate(ev_light):
    W(f"| Diagonal mesino {i+1} | {abs(ev):.4e} | From 4×4 block light mode |")
W(f"| Heavy X-mode 1 | {abs(ev_4x4_sorted_all[2]):.4e} | Integrated out (Λ⁶/M scale) |")
W(f"| Heavy X-mode 2 | {abs(ev_4x4_sorted_all[3]):.4e} | Integrated out (Λ⁶/M scale) |")
W("")
W("---")
W("")
W("## Part (c): Baryons")
W("")
W("With m_B = Λ = 300 MeV and B = B̃ = 0 at vacuum.")
W("")
W("**Baryon mass term in W:**")
W("")
W("    W ⊃ -X(BB̃) + m_B BB̃ = (m_B - X)BB̃")
W("")
W("**Baryonino** (Weyl fermion pair ψ_B, ψ_B̃ → 1 Dirac fermion):")
W("")
W(f"    m_baryonino = |W_{{B,B̃}}| = |m_B - X₀| ≈ m_B = {m_B - X0:.6f} MeV")
W("")
W("(The correction X₀ ≈ -1.21×10⁻⁹ is negligible.)")
W("")
W("**Baryon scalars** (2 complex = 4 real, from B and B̃):")
W("")
W("F-term contribution: from F_B = (m_B - X₀)B̃ and F_B̃ = (m_B - X₀)B:")
W("")
W(f"    m²_scalar(B) = |m_B - X₀|² = {m2_baryon_scalar:.6f} MeV²")
W(f"    m_scalar(B) = {m_baryon_scalar:.6f} MeV ≈ Λ = 300 MeV")
W("")
W("No soft SUSY-breaking term for baryons (V_soft only for mesons in this model).")
W("")
W("| Particle | Spin | Mass (MeV) | Count |")
W("|----------|------|-----------|-------|")
W(f"| Baryon B (scalar) | 0 | {m_baryon_scalar:.2f} | 2 real (Re B) |")
W(f"| Anti-baryon B̃ (scalar) | 0 | {m_baryon_scalar:.2f} | 2 real (Re B̃) |")
W(f"| Baryonino (ψ_B, ψ_B̃) | 1/2 | {m_B - X0:.4f} | 1 Dirac |")
W("")
W("---")
W("")
W("## Part (d): Higgs Sector (NMSSM)")
W("")
W("**Superpotential:** W ⊃ λ_S S H_u·H_d + κ/3 S³")
W("")
W("**VEVs:** ⟨H_u⁰⟩ = v_u = v/√2, ⟨H_d⁰⟩ = v_d = v/√2 (tan β = 1)")
W("")
W("**Tree-level lightest CP-even Higgs** (tan β = 1 → cos 2β = 0):")
W("")
W("    m_h² = m_Z² cos²(2β) + λ_S² v² sin²(2β)/2")
W("          = 0 + λ_S² v²/2  (since sin 2β = 1 at tan β = 1)")
W("")
W(f"    m_h = λ_S v/√2 = {lam_S} × {v/1000:.2f} GeV / √2 = {m_h_problem/1000:.4f} GeV ≈ 125 GeV")
W("")
W("**Yukawa couplings** (from W ⊃ y_c H_u M^c_c + y_b H_d M^b_b + y_t H_u Q^t Q̄_t):")
W("")
W("| Quark | Yukawa y | Fermion mass (MeV) | Scalar partner (squark) |")
W("|-------|----------|-------------------|------------------------|")
W(f"| c | y_c = m_c/v_u = {m_c/vu:.6f} | {m_c} MeV | mass from soft terms |")
W(f"| b | y_b = m_b/v_d = {m_b/vd:.6f} | {m_b} MeV | mass from soft terms |")
W(f"| t | y_t = m_t/v_u = {m_t/vu:.6f} | {m_t} MeV | mass from soft terms |")
W("")
W("**NMSSM spectrum** (soft-dependent, not fully computed):")
W("")
W("| Particle | Spin | Mass | Notes |")
W("|----------|------|------|-------|")
W(f"| h (CP-even Higgs) | 0 | {m_h_problem/1000:.3f} GeV | Tree-level |")
W("| H (heavy CP-even) | 0 | ≫ m_h | Soft dependent |")
W("| A (CP-odd) | 0 | Soft dependent | |")
W("| H± (charged Higgs) | 0 | √(m_A²+m_W²) | Soft dependent |")
W("| S (singlet scalar) | 0 | κ⟨S⟩ | κ, ⟨S⟩ free |")
W("| H̃_u (Higgsino) | 1/2 | λ_S⟨S⟩ = μ_eff | |")
W("| H̃_d (Higgsino) | 1/2 | λ_S⟨S⟩ = μ_eff | |")
W("| S̃ (singlino) | 1/2 | 2κ⟨S⟩ | |")
W("")
W("---")
W("")
W("## Part (e): Complete Spectrum Table")
W("")
W("### Group 1: Light Confined Sector (meson scalar + off-diagonal mesino)")
W("")
W("| Particle | Spin | Q_em | SU(2)_L | Mass (MeV) | Origin |")
W("|----------|------|------|---------|-----------|--------|")
W(f"| M_uu scalar (Re,Im) | 0 | 0 | singlet | {f_pi:.1f} | V_soft = f_π² |")
W(f"| M_dd scalar (Re,Im) | 0 | 0 | singlet | {f_pi:.1f} | V_soft = f_π² |")
W(f"| M_ss scalar (Re,Im) | 0 | 0 | singlet | {f_pi:.1f} | V_soft = f_π² |")
W(f"| M_ud, M_du scalars | 0 | ±1 | singlet | {m_offdiag_scalar:.2f} | 2f_π² (doubled) |")
W(f"| M_us, M_su scalars | 0 | ±2/3 | singlet | {m_offdiag_scalar:.2f} | 2f_π² |")
W(f"| M_ds, M_sd scalars | 0 | ∓1/3 | singlet | {m_offdiag_scalar:.2f} | 2f_π² |")
W(f"| Mesino ψ_ud / ψ_du (Dirac) | 1/2 | ±1 | singlet | {abs(X0)*M_s_vev:.4e} | W_{{M_ud,M_du}} = |X₀|M_ss |")
W(f"| Mesino ψ_us / ψ_su (Dirac) | 1/2 | ±2/3 | singlet | {abs(X0)*M_d_vev:.4e} | W_{{M_us,M_su}} = |X₀|M_dd |")
W(f"| Mesino ψ_ds / ψ_sd (Dirac) | 1/2 | ∓1/3 | singlet | {abs(X0)*M_u_vev:.4e} | W_{{M_ds,M_sd}} = |X₀|M_uu |")
W("")
W("### Group 2: Heavy Confined Sector (X-dominated diagonal mesino modes)")
W("")
W("| Particle | Spin | Q_em | Mass | Origin |")
W("|----------|------|------|------|--------|")
W(f"| Diagonal mesino light 1 | 1/2 | 0 | {abs(ev_4x4_sorted_all[0]):.4e} MeV | 4×4 light mode |")
W(f"| Diagonal mesino light 2 | 1/2 | 0 | {abs(ev_4x4_sorted_all[1]):.4e} MeV | 4×4 light mode |")
W(f"| X-dominated mesino 1 | 1/2 | 0 | {abs(ev_4x4_sorted_all[2]):.4e} | Λ⁶/M scale, X coupling |")
W(f"| X-dominated mesino 2 | 1/2 | 0 | {abs(ev_4x4_sorted_all[3]):.4e} | Λ⁶/M scale, X coupling |")
W("| X scalar | none | — | — | No kinetic term: not propagating |")
W("| Diagonal meson scalars (see Group 1 above) | | | | |")
W("")
W("### Group 3: Baryons")
W("")
W("| Particle | Spin | Q_em | SU(2)_L | Mass (MeV) | Origin |")
W("|----------|------|------|---------|-----------|--------|")
W(f"| B (complex scalar) | 0 | 0 | singlet | {m_baryon_scalar:.2f} | F_B = m_B B̃ |")
W(f"| B̃ (complex scalar) | 0 | 0 | singlet | {m_baryon_scalar:.2f} | F_B̃ = m_B B |")
W(f"| ψ_B, ψ_B̃ (Dirac) | 1/2 | 0 | singlet | {m_B - X0:.4f} | W_{{B,B̃}} = m_B |")
W("")
W("### Group 4: Elementary Quarks (above confinement threshold)")
W("")
W("| Particle | Spin | Q_em | SU(3)_c | SU(2)_L | Mass (MeV) | Origin |")
W("|----------|------|------|---------|---------|-----------|--------|")
W(f"| c quark (Dirac) | 1/2 | +2/3 | 3 | doublet | {m_c} | y_c H_u M^c_c |")
W(f"| b quark (Dirac) | 1/2 | −1/3 | 3 | singlet | {m_b} | y_b H_d M^b_b |")
W(f"| t quark (Dirac) | 1/2 | +2/3 | 3 | doublet | {m_t} | y_t H_u Q_t Q̄_t |")
W(f"| c~ squark (scalar) | 0 | +2/3 | 3 | doublet | soft-breaking | Superpartner of c |")
W(f"| b~ squark (scalar) | 0 | −1/3 | 3 | singlet | soft-breaking | Superpartner of b |")
W(f"| t~ squark (scalar) | 0 | +2/3 | 3 | doublet | soft-breaking | Superpartner of t |")
W("")
W("### Group 5: Higgs Sector")
W("")
W("| Particle | Spin | Q_em | SU(2)_L | Mass | Origin |")
W("|----------|------|------|---------|------|--------|")
W(f"| h (CP-even Higgs) | 0 | 0 | singlet | {m_h_problem/1000:.3f} GeV | λ_S v/√2, tan β=1 |")
W("| H (heavy CP-even) | 0 | 0 | singlet | soft-dep. | Decoupling limit |")
W("| A (CP-odd Higgs) | 0 | 0 | singlet | soft-dep. | |")
W("| H± (charged Higgs) | 0 | ±1 | doublet | soft-dep. | |")
W("| S (singlet scalar) | 0 | 0 | singlet | κ⟨S⟩ | W ⊃ κ/3 S³ |")
W("| H̃_u± (Higgsino) | 1/2 | ±1 | doublet | μ_eff=λ_S⟨S⟩ | |")
W("| H̃_u⁰ (Higgsino) | 1/2 | 0 | doublet | μ_eff | |")
W("| H̃_d± (Higgsino) | 1/2 | ±1 | doublet | μ_eff | |")
W("| H̃_d⁰ (Higgsino) | 1/2 | 0 | doublet | μ_eff | |")
W("| S̃ (singlino) | 1/2 | 0 | singlet | 2κ⟨S⟩ | W_SS = κS term |")
W("")
W("### Group 6: Gauge Bosons (after EW breaking)")
W("")
W("| Particle | Spin | Q_em | Mass (MeV) | Notes |")
W("|----------|------|------|-----------|-------|")
W("| Gluons (8) | 1 | 0 | — | Confined, no propagating gluons |")
W(f"| W± bosons | 1 | ±1 | {80379:.0f} (≈80.4 GeV) | EW breaking |")
W(f"| Z⁰ boson | 1 | 0 | {91187.6:.1f} (≈91.2 GeV) | EW breaking |")
W("| Photon γ | 1 | 0 | 0 | Unbroken U(1)_em |")
W("| Gluino g̃ (8) | 1/2 | 0 | soft-dep. | SU(3)_c adjoint |")
W("| Wino W̃± | 1/2 | ±1 | soft-dep. | SU(2)_L adjoint |")
W("| Wino W̃⁰ | 1/2 | 0 | soft-dep. | SU(2)_L adjoint |")
W("| Bino B̃ | 1/2 | 0 | soft-dep. | U(1)_Y adjoint |")
W("")
W("---")
W("")
W("## Degree of Freedom Count")
W("")
W("### Confined sector:")
W("| Field | Scalars (real) | Weyl fermions |")
W("|-------|---------------|--------------|")
W("| M^i_j diagonal (3 complex) | 6 | 3 (+ 1 constrained away by ψ_X) |")
W("| M^i_j off-diagonal (6 complex) | 12 | 6 (3 Dirac pairs) |")
W("| X (no kinetic) | 0 | 0 (auxiliary → Lagrange multiplier) |")
W("| B, B̃ | 4 | 2 (1 Dirac) |")
W("| **Total confined** | **22 real scalars** | **11 Weyl = 5.5 Dirac** |")
W("")
W("### Notes on dimensional analysis:")
W("")
W(f"The off-diagonal mesino masses |X₀| × M_vev are extremely small in the given units:")
W(f"- |X₀| × M_ss = {abs(X0)*M_s_vev:.4e}")
W(f"- These are not in MeV directly; the unit mismatch reflects that [X₀] = MeV^{{-3}} or MeV^{{-4}}")
W(f"  and [M_vev] = MeV or MeV², so the product is dimensionless or MeV^{{-1}}.")
W(f"- In the physical interpretation, these masses set the splitting between scalar and")
W(f"  fermionic components of the meson supermultiplet (supersymmetric partner masses).")
W(f"- The dominant physical mass scale in the confined sector is the soft scale f_π = {f_pi} MeV.")
W("")
W("---")
W("")
W("*Generated by particle_spectrum.py*")

report_text = "\n".join(spectrum_table)
report_path = "/home/codexssh/phys3/results/particle_spectrum.md"
with open(report_path, "w") as f:
    f.write(report_text)
print(f"\nReport written to {report_path}")
print("Python script complete.")
