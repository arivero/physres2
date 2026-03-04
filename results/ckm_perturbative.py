#!/usr/bin/env python3
"""
Perturbative CKM mixing from off-diagonal meson VEVs.

Setup: N_f = N_c = 3 SQCD with baryons B = B~ = 0 enforced.
Superpotential:
    W = sum_i m_i M^i_i + X(det M - Lambda^6)
        + y_c H_u M^d_d + y_b H_d M^s_s + lambda X H_u H_d

with soft breaking V_soft = f_pi^2 Tr(M†M).

Diagonal seesaw vacuum: M^i_j = delta^i_j C/m_i,
  C = Lambda^2 (m_u m_d m_s)^{1/3}, F_X = det M - Lambda^6 = 0.
H_u = H_d = 0 before electroweak symmetry breaking.

Task: compute the 3x3 Hessian for off-diagonal meson entries at this
vacuum, identify tachyonic directions, find the nonzero VEVs by
including quartic terms, perform SVD, extract CKM angles.
"""

import numpy as np
from scipy.optimize import minimize_scalar, brentq
from scipy.linalg import svd
import sys

np.set_printoptions(precision=8, suppress=False, linewidth=120)

# ============================================================
# Physical parameters (MeV units throughout)
# ============================================================
m_u   =     2.16        # MeV
m_d   =     4.67        # MeV
m_s   =    93.4         # MeV
Lambda =   300.0        # MeV
v     = 246220.0        # MeV  (electroweak VEV)
m_c   =  1270.0         # MeV  (charm mass, for Yukawa)
m_b   =  4180.0         # MeV  (bottom mass, for Yukawa)
y_c   = 2.0 * m_c / v  # charm Yukawa
y_b   = 2.0 * m_b / v  # bottom Yukawa
lam   = 0.72            # lambda coupling in W term lambda X H_u H_d
f_pi  = 92.0            # MeV  (pion decay constant = soft mass)

masses = np.array([m_u, m_d, m_s])

# ============================================================
# Diagonal seesaw vacuum (B = B~ = 0 enforced)
# ============================================================
# M^i_j = delta^i_j M_i  with  M_i = C / m_i
# C = Lambda^2 (m_u m_d m_s)^{1/3}
# X satisfies  det M = Lambda^6  =>  X from F_{M^i_i} = 0
# F_{M^i_i} = m_i + X * (adj M)^i_i = 0
#
# For diagonal M: (adj M)^i_i = product of the OTHER two diagonal entries
# So X = -m_i / (M_j * M_k)  (same for all i by the seesaw construction)

C = Lambda**2 * np.prod(masses)**(1.0/3.0)
M_vev = C / masses  # [M_u, M_d, M_s]
Lambda6 = Lambda**6

# Verify det M = Lambda^6
det_check = np.prod(M_vev)
print("=" * 70)
print("DIAGONAL SEESAW VACUUM")
print("=" * 70)
print(f"  C = Lambda^2 (m_u m_d m_s)^(1/3) = {C:.6f} MeV^2")
print(f"  M_u = C/m_u = {M_vev[0]:.4f} MeV")
print(f"  M_d = C/m_d = {M_vev[1]:.4f} MeV")
print(f"  M_s = C/m_s = {M_vev[2]:.4f} MeV")
print(f"  det M = {det_check:.6e}  Lambda^6 = {Lambda6:.6e}  ratio = {det_check/Lambda6:.12f}")
print()

# X at the seesaw vacuum: from F_{M^u_u} = m_u + X * M_d * M_s = 0
X_vev = -m_u / (M_vev[1] * M_vev[2])
print(f"  X_vev = -m_u / (M_d M_s) = {X_vev:.6e} MeV^(-4)")

# Verify consistency via F_{M^d_d} and F_{M^s_s}:
Fdd_check = m_d + X_vev * M_vev[0] * M_vev[2]
Fss_check = m_s + X_vev * M_vev[0] * M_vev[1]
print(f"  F_{{M^d_d}} = m_d + X M_u M_s = {Fdd_check:.6e}  (should be 0 for pure seesaw)")
print(f"  F_{{M^s_s}} = m_s + X M_u M_d = {Fss_check:.6e}  (should be 0 for pure seesaw)")
print()

# Note: The exact seesaw condition X M_d M_s = -m_u AND X M_u M_s = -m_d AND X M_u M_d = -m_s
# is self-consistent iff m_u * M_u = m_d * M_d = m_s * M_s = C  (true by construction).
# So the diagonal vacuum is a consistent critical point.

# ============================================================
# Scalar potential at H_u = H_d = 0
# ============================================================
# V = sum_{i,j} |F_{M^i_j}|^2 + |F_X|^2 + |F_{H_u}|^2 + |F_{H_d}|^2
#     + f_pi^2 Tr(M†M)
#
# F_{M^i_j} = dW/dM^i_j = m_i delta^i_j + X * cof(M)^j_i
#             + y_c H_u delta^i_d delta^j_d + y_b H_d delta^i_s delta^j_s
# At H_u = H_d = 0:
#   F_{M^i_j} = m_i delta^i_j + X * cof(M)^j_i
#
# F_X = dW/dX = det M - Lambda^6  (+ lambda H_u H_d, but H = 0)
# F_{H_u} = dW/dH_u = y_c M^d_d  (+ lambda X H_d, but H_d = 0)
# F_{H_d} = dW/dH_d = y_b M^s_s  (+ lambda X H_u, but H_u = 0)
#
# At the diagonal vacuum with H_u = H_d = 0:
#   F_{M^i_i} = m_i + X * (M_j M_k) = 0  (by seesaw construction)
#   F_{M^i_j} = X * cof(M)^j_i  for i != j
#
# For a diagonal matrix: cof(M)^j_i = 0 for i != j
# So F_{M^i_j} = 0 for i != j at the diagonal vacuum.  Good.

print("=" * 70)
print("F-TERMS AT THE DIAGONAL VACUUM")
print("=" * 70)
# Build cofactor matrix for diagonal M
cof = np.zeros((3, 3))
cof[0, 0] = M_vev[1] * M_vev[2]
cof[1, 1] = M_vev[0] * M_vev[2]
cof[2, 2] = M_vev[0] * M_vev[1]
# Off-diagonal cofactors of a diagonal matrix are all zero.

F_M_diag = np.zeros((3, 3))
for i in range(3):
    F_M_diag[i, i] = masses[i] + X_vev * cof[i, i]

F_X_diag = det_check - Lambda6
F_Hu_diag = y_c * M_vev[1]  # y_c * M_d  (index 1 = d)
F_Hd_diag = y_b * M_vev[2]  # y_b * M_s  (index 2 = s)

labels = ['u', 'd', 's']
for i in range(3):
    for j in range(3):
        val = F_M_diag[i, j]
        print(f"  F_{{M^{labels[i]}_{labels[j]}}} = {val:.6e}  {'(should be 0)' if val == 0 else ''}")
print(f"  F_X     = {F_X_diag:.6e}  (should be 0)")
print(f"  F_{{H_u}} = y_c M_d = {F_Hu_diag:.6e}")
print(f"  F_{{H_d}} = y_b M_s = {F_Hd_diag:.6e}")
print()
print("  Note: F_{H_u} and F_{H_d} are nonzero at H=0: these drive H into a")
print("  nonzero VEV at EW breaking scale, separate from the QCD-scale analysis.")
print()

# ============================================================
# HESSIAN FOR OFF-DIAGONAL MESON ENTRIES
# ============================================================
#
# Write M^i_j = delta^i_j M_i + epsilon^i_j  (epsilon are the off-diagonal perturbations)
# Focus on the 6 off-diagonal entries: (0,1),(0,2),(1,0),(1,2),(2,0),(2,1)
# = (ud, us, du, ds, su, sd) in labels.
#
# To second order in epsilon:
#
#   V = V_0 + (1/2) sum_{i!=j, k!=l} H_{ij,kl} epsilon^i_j epsilon^k_l + O(epsilon^3)
#
# We compute H analytically.
#
# The relevant F-terms for off-diagonal perturbations:
#
# 1) F_{M^i_j} for i != j at off-diagonal epsilon^i_j:
#
#    F_{M^i_j} = X * cof(M)^j_i  evaluated at the perturbed M
#
#    To first order in epsilon:
#    cof(M)^j_i  for i != j in a mostly-diagonal matrix:
#    The (j,i) cofactor of M is obtained by deleting row j, column i.
#    For 3x3 with only epsilon^i_j nonzero (and all other off-diag = 0):
#
#    If we turn on only epsilon^{a}_{b} (a != b):
#    The cofactor cof^{b}_{a} = product of elements NOT in row b, column a
#                              = (diagonal entry at the remaining index) * epsilon (cross terms)
#                              + epsilon^2 term
#    Actually for the cofactor of a 3x3 matrix (deleting row j, col i),
#    we get a 2x2 determinant.
#
#    Let me be explicit. Take (i,j) = (0,1) [u-d off-diagonal]:
#    F_{M^0_1} = X * cof(M)^{1}_{0}
#    cof(M)^{1}_{0} = minor obtained by deleting row 1, col 0 of M
#    = det([[M[0,1], M[0,2]], [M[2,1], M[2,2]]])
#    = M[0,1]*M[2,2] - M[0,2]*M[2,1]
#    At the diagonal vacuum: M[0,1]=epsilon^0_1, M[0,2]=epsilon^0_2,
#                             M[2,1]=epsilon^2_1, M[2,2]=M_s
#    => cof^1_0 = epsilon^0_1 * M_s - epsilon^0_2 * epsilon^2_1
#    => F_{M^0_1} = X_vev * (epsilon^0_1 * M_s - epsilon^0_2 * epsilon^2_1) + O(eps^3)
#
#    First-order: F_{M^i_j}^{(1)} = X_vev * (delta M^i_j) * M_k  [k = the third index]
#
# Let me define a systematic expansion.
# Index notation: i,j,k are a permutation of {0,1,2}.
# For off-diagonal perturbation epsilon^a_b (a != b), define k = the remaining index (not a, not b).

print("=" * 70)
print("HESSIAN OF V FOR OFF-DIAGONAL MESON ENTRIES")
print("=" * 70)
print()
print("  Expanding V to second order in epsilon^i_j (i != j).")
print()

# The 6 off-diagonal entries, labeled as (row, col) pairs:
od_pairs = [(0,1),(0,2),(1,0),(1,2),(2,0),(2,1)]
# Abbreviation: pair (a,b) is labeled
od_labels = ['ud', 'us', 'du', 'ds', 'su', 'sd']

# For each off-diagonal pair (a,b), the "remaining" index is:
def third(a, b):
    return {0,1,2}.difference({a,b}).pop()

# ============================================================
# Analytic Hessian computation
# ============================================================
# Contributions to H^2 = mass^2 matrix:
#
# (A) Soft term: +f_pi^2 * |epsilon^a_b|^2   [diagonal in off-diag sector]
#     => contributes +f_pi^2 to H_{ab,ab}
#
# (B) |F_{M^a_b}|^2 terms (a != b):
#     F_{M^a_b} = X_vev * cof(M)^b_a
#     The leading-order contribution to cof(M)^b_a:
#       Delete row b, col a from M. The 2x2 remaining submatrix has:
#         - Two diagonal elements from M_vev at the other rows/cols
#         - Off-diagonal epsilon entries
#
#     For a 3x3 matrix, if we delete row b and column a (a != b):
#     The remaining 2x2 involves rows != b and columns != a.
#     Let rows_left = {r : r != b}, cols_left = {c : c != a}.
#     The minor is:
#       det([[M[r1,c1], M[r1,c2]], [M[r2,c1], M[r2,c2]]])
#     where r1 < r2 are the two rows != b, c1 < c2 are the two columns != a.
#
#     At the diagonal vacuum (only epsilon^a_b turned on):
#     The non-epsilon elements are M[i,i] = M_i.
#     The epsilon^a_b element appears at position (row=a, col=b) of M.
#     After deleting row b and col a:
#       - The retained entries from M are: the diagonal ones at rows != b, cols != a
#       - The retained epsilon is: epsilon^a_b appears at (row a, col b),
#         but col b is retained only if b != a (always true) -- wait:
#         We deleted COL a, not col b. So col b is retained.
#         So epsilon^a_b at position (a, b) IS in the 2x2 submatrix (row a retained since a != b,
#         col b retained since b != a).
#
#     Let me compute explicitly for each (a,b):

# I'll compute the full 6x6 Hessian numerically using the analytic formula,
# then verify with a numerical second-derivative check.

# Define the scalar potential function V(M) with diagonal enforced + off-diagonal perturbation
def V_full(epsilon_vec, M_diag_vev=M_vev, X=X_vev,
           f_pi_sq=f_pi**2, Lambda6=Lambda6,
           y_c=y_c, y_b=y_b):
    """
    Full scalar potential V at M = diag(M_vev) + off-diagonal epsilon,
    with H_u = H_d = 0 and B = B~ = 0.

    epsilon_vec: 6-vector for [M^0_1, M^0_2, M^1_0, M^1_2, M^2_0, M^2_1]
                              [ud,    us,    du,    ds,    su,    sd   ]
    """
    # Build 3x3 meson matrix
    M = np.diag(M_diag_vev.copy())
    M[0, 1] = epsilon_vec[0]  # ud
    M[0, 2] = epsilon_vec[1]  # us
    M[1, 0] = epsilon_vec[2]  # du
    M[1, 2] = epsilon_vec[3]  # ds
    M[2, 0] = epsilon_vec[4]  # su
    M[2, 1] = epsilon_vec[5]  # sd

    # Cofactors of M
    # cof[j, i] = cofactor of M at position (i, j)  [for F_{M^i_j} = X * cof^j_i]
    cof_M = np.zeros((3, 3))
    # Using direct formula
    cof_M[0, 0] = M[1,1]*M[2,2] - M[1,2]*M[2,1]
    cof_M[1, 0] = -(M[0,1]*M[2,2] - M[0,2]*M[2,1])  # del row 0, col 1; cofactor sign
    cof_M[2, 0] = M[0,1]*M[1,2] - M[0,2]*M[1,1]
    cof_M[0, 1] = -(M[1,0]*M[2,2] - M[1,2]*M[2,0])
    cof_M[1, 1] = M[0,0]*M[2,2] - M[0,2]*M[2,0]
    cof_M[2, 1] = -(M[0,0]*M[1,2] - M[0,2]*M[1,0])
    cof_M[0, 2] = M[1,0]*M[2,1] - M[1,1]*M[2,0]
    cof_M[1, 2] = -(M[0,0]*M[2,1] - M[0,1]*M[2,0])
    cof_M[2, 2] = M[0,0]*M[1,1] - M[0,1]*M[1,0]

    # F_{M^i_j} = m_i delta^i_j + X * cof^j_i
    # (at H_u = H_d = 0, the Yukawa terms contribute only to F_{M^d_d} and F_{M^s_s}
    #  via y_c H_u and y_b H_d, which vanish)
    F_M = np.zeros((3, 3))
    for i in range(3):
        for j in range(3):
            F_M[i, j] = X * cof_M[j, i]
            if i == j:
                F_M[i, j] += masses[i]

    # F_X = det M - Lambda^6
    detM = (M[0,0]*(M[1,1]*M[2,2]-M[1,2]*M[2,1])
           -M[0,1]*(M[1,0]*M[2,2]-M[1,2]*M[2,0])
           +M[0,2]*(M[1,0]*M[2,1]-M[1,1]*M[2,0]))
    F_X = detM - Lambda6

    # F_{H_u} = y_c M^d_d = y_c M[1,1],  F_{H_d} = y_b M^s_s = y_b M[2,2]
    F_Hu = y_c * M[1, 1]
    F_Hd = y_b * M[2, 2]

    # Scalar potential (all meson entries, including diagonal)
    V = (np.sum(F_M**2) + F_X**2 + F_Hu**2 + F_Hd**2
         + f_pi_sq * np.sum(M**2))

    return V


# Numerical Hessian using central differences
def numerical_hessian_6x6(h=1e-3):
    """6x6 Hessian for the off-diagonal sector at the diagonal vacuum."""
    x0 = np.zeros(6)
    V0 = V_full(x0)
    H = np.zeros((6, 6))
    for i in range(6):
        for j in range(i, 6):
            xpp = x0.copy(); xpp[i] += h; xpp[j] += h
            xpm = x0.copy(); xpm[i] += h; xpm[j] -= h
            xmp = x0.copy(); xmp[i] -= h; xmp[j] += h
            xmm = x0.copy(); xmm[i] -= h; xmm[j] -= h
            H[i, j] = (V_full(xpp) - V_full(xpm) - V_full(xmp) + V_full(xmm)) / (4*h*h)
            H[j, i] = H[i, j]
    return H, V0


# Use a step size comparable to the smallest natural scale
h_step = min(M_vev) * 1e-4  # fraction of smallest meson VEV
print(f"  Computing numerical Hessian (h = {h_step:.4e} MeV)...")
H6, V0 = numerical_hessian_6x6(h=h_step)
print()

print("  Off-diagonal Hessian (MeV^2) labeled as [ud, us, du, ds, su, sd]:")
header = f"{'':>6}" + "".join(f"{l:>14}" for l in od_labels)
print("  " + header)
for i in range(6):
    row = f"  {od_labels[i]:>6}" + "".join(f"{H6[i,j]:14.4e}" for j in range(6))
    print(row)
print()

# Eigenvalues of the 6x6 Hessian
evals_6, evecs_6 = np.linalg.eigh(H6)
print("  Eigenvalues of the 6x6 off-diagonal Hessian (MeV^2):")
for k, ev in enumerate(evals_6):
    status = "TACHYONIC" if ev < 0 else "stable"
    print(f"    lambda_{k+1} = {ev:14.6e}  [{status}]")
print()

# Identify tachyonic directions
tachyonic_mask = evals_6 < 0
n_tachyonic = np.sum(tachyonic_mask)
print(f"  Number of tachyonic directions: {n_tachyonic}")
print()

# ============================================================
# ANALYTIC HESSIAN — sector by sector
# ============================================================
print("=" * 70)
print("ANALYTIC HESSIAN: SECTOR-BY-SECTOR DECOMPOSITION")
print("=" * 70)
print()
print("  The 6x6 off-diagonal sector decomposes by flavor pair.")
print("  For each pair (a,b) with a != b, define k = third index.")
print()
print("  Key contributions to m^2_{(a,b)} (mass-squared of off-diagonal mode):")
print()
print("  (I)  Soft term:           +f_pi^2")
print("  (II) |F_{M^a_b}|^2:       +X_vev^2 * M_k^2   (from cofactor)")
print("         [note: X_vev < 0, X_vev^2 > 0]")
print("  (III)|F_{M^b_a}|^2:       +X_vev^2 * M_k^2   [same by symmetry if M is symmetric]")
print("  (IV) |F_X|^2 terms:       from det(M) off-diagonal expansion")
print("         det(M) to second order in epsilon^a_b:")
print("           d(det M)/d(epsilon^a_b)|_diag = 0  (off-diag of adj are zero at diag M)")
print("           d^2(det M)/(d epsilon^a_b d epsilon^a_b) = 0  [diagonal of Hessian of det]")
print("           d^2(det M)/(d epsilon^a_b d epsilon^b_a) = M_k  [cross term]")
print("         => F_X = det M - Lambda^6 = (sum of cross terms) * M_k * epsilon^a_b * epsilon^b_a")
print("         Contribution to H_{(a,b),(b,a)}: +2 * M_k^2 (from |F_X|^2 cross term)")
print("         Contribution to H_{(a,b),(a,b)}: 0  [F_X is zero at leading order]")
print("  (V)  Yukawa F-terms:      |F_{H_u}|^2 = y_c^2 M_d^2, |F_{H_d}|^2 = y_b^2 M_s^2")
print("         These do NOT depend on off-diagonal entries at first order.")
print()

# Compute the analytic 2x2 block Hessian for each sector (a,b) and (b,a)
def analytic_block_hessian(a, b):
    """
    Returns the 4x4 block of the Hessian for the (a,b), (b,a) pair,
    in the basis [epsilon^a_b, epsilon^b_a].
    Actually returns 2x2 for the (a,b)-(b,a) mixed block.
    """
    k = third(a, b)
    M_a = M_vev[a]
    M_b = M_vev[b]
    M_k = M_vev[k]

    # At the diagonal vacuum, the Hessian mixing (a,b) and (b,a):
    # d^2V / d(epsilon^a_b)^2:
    #   Soft: +f_pi^2
    #   |F_{M^a_b}|^2: F_{M^a_b} = X_vev * cof^b_a = X_vev * (M_k * epsilon^a_b - ...)
    #     cof^b_a = minor of M with row b, col a deleted
    #     = det([[M[a,b], M[a,k]], [M[k,b], M[k,k]]])     [since rows {not b} = {a,k}, cols {not a} = {b,k}]
    #     At diagonal + single epsilon^a_b:
    #       M[a,b] = epsilon^a_b, M[a,k] = 0 (diagonal for a=a col k, off-diag entry = 0)
    #       M[k,b] = 0 (off-diagonal entry epsilon^k_b = 0)
    #       M[k,k] = M_k
    #     => cof^b_a = epsilon^a_b * M_k - 0 = epsilon^a_b * M_k
    #   d(F_{M^a_b})/d(epsilon^a_b) = X_vev * M_k
    #   d^2|F_{M^a_b}|^2/d(epsilon^a_b)^2 = 2 * X_vev^2 * M_k^2
    #
    #   |F_{M^b_a}|^2: F_{M^b_a} = X_vev * cof^a_b
    #     cof^a_b = minor of M with row a, col b deleted
    #     Rows {not a} = {b,k}, cols {not b} = {a,k}
    #     = det([[M[b,a], M[b,k]], [M[k,a], M[k,k]]])
    #     M[b,a] = epsilon^b_a = 0 (only (a,b) is nonzero), M[b,k] = 0, M[k,a] = 0, M[k,k] = M_k
    #     => cof^a_b = 0
    #   d^2|F_{M^b_a}|^2/d(epsilon^a_b)^2 = 0
    #
    #   |F_X|^2: det M at second order in epsilon^a_b alone:
    #     d^2(det M)/d(epsilon^a_b)^2 = 0 [because det is multilinear, each element appears once]
    #     d|F_X|^2/d(epsilon^a_b)^2 = 2*(d(det M)/d(epsilon^a_b))^2 / (=0 at vac) = 0 at linear order
    #   Actually: det M = sum_perms sign * product of elements.
    #   The only permutation containing M[a,b] also contains either M[b,?] or M[?,a] terms
    #   from off-diagonal, which are all zero at the vacuum. So (d det M/d epsilon^a_b)|_{vac} = 0.
    #   Therefore (d F_X/d epsilon^a_b) = 0, and d^2|F_X|^2/d(epsilon^a_b)^2 = 0.
    #
    #   Summary: H_{(a,b),(a,b)} = f_pi^2 + 2 * X_vev^2 * M_k^2

    H_aa = f_pi**2 + 2.0 * X_vev**2 * M_k**2

    # d^2V / d(epsilon^a_b) d(epsilon^b_a):
    #   Soft: 0 (soft term is diagonal in field space)
    #   |F_{M^a_b}|^2 cross term:
    #     F_{M^a_b} = X_vev * cof^b_a, cof^b_a now from full expansion
    #     When we also have epsilon^b_a != 0:
    #     cof^b_a = det([[M[a,b], M[a,k]], [M[k,b], M[k,k]]])
    #             = epsilon^a_b * M_k - 0 = epsilon^a_b * M_k  [to leading order]
    #     No cross term from |F_{M^a_b}|^2 alone at O(epsilon^2).
    #     Wait — I need to be more careful. The cross derivative:
    #     d^2|F_{M^a_b}|^2 / (d eps^a_b d eps^b_a) = 2 * (d F^a_b/d eps^a_b)(d F^a_b/d eps^b_a)
    #       d F^a_b/d eps^b_a: does F^a_b depend on eps^b_a?
    #       cof^b_a involves row {a,k}, col {b,k}. The entry eps^b_a sits at (row b, col a) —
    #       row b is DELETED in cof^b_a (we delete row b). So eps^b_a does NOT appear in cof^b_a.
    #       => d F^a_b/d eps^b_a = 0. No cross term from |F^a_b|^2.
    #
    #   |F_{M^b_a}|^2 cross term:
    #     Similarly, d F^b_a/d eps^a_b: cof^a_b involves rows {b,k}, cols {a,k}.
    #     Entry eps^a_b sits at (row a, col b). Row a is deleted in cof^a_b. So no appearance.
    #     => d F^b_a/d eps^a_b = 0. No cross term from |F^b_a|^2.
    #
    #   |F_X|^2 cross term:
    #     det M to second order in epsilon^a_b AND epsilon^b_a:
    #     The permutation contributing M[a,b]*M[b,a] is the transposition (a b), which gives
    #     sign * M[a,b] * M[b,a] * M[k,k] = -(−1)^{...} * eps^a_b * eps^b_a * M_k
    #     For 3-cycle permutation: sign of transposition (a b) in S_3 is -1.
    #     det M = ... + (-1) * eps^a_b * eps^b_a * M_k + ...
    #     So d^2(det M)/(d eps^a_b d eps^b_a) = -M_k
    #     => d^2|F_X|^2/(d eps^a_b d eps^b_a) = 2 * (d det M/d eps^a_b)*(d det M/d eps^b_a)|_{eps=0}
    #        + 2 * F_X * (d^2 det M/d eps^a_b d eps^b_a)|_{eps=0}
    #     At the vacuum: F_X = 0, and (d det M/d eps^a_b)|_{eps=0} = 0 (as shown above).
    #     So the cross term also = 0 at the vacuum to this order.
    #
    #   Hmm, this suggests H_{(a,b),(b,a)} = 0. But the numerical Hessian will show if that's true.
    #   Actually, wait: there are also cross terms from F_{M^a_a} and F_{M^b_b} and F_{M^k_k}.
    #
    #   Reconsider: the DIAGONAL F-terms F_{M^i_i} also depend on off-diagonal entries!
    #   F_{M^a_a} = m_a + X * cof(M)^a_a
    #   cof(M)^a_a = det of 2x2 submatrix with row a, col a deleted = det[[M[b,b], M[b,k]], [M[k,b], M[k,k]]]
    #   When eps^a_b and eps^b_a are nonzero:
    #   M[b,b] = M_b, M[b,k] = 0 (eps^b_k = 0), M[k,b] = 0 (eps^k_b = 0), M[k,k] = M_k
    #   => cof(M)^a_a = M_b * M_k  [unchanged by eps^a_b or eps^b_a since those don't appear in the a-row or a-col of the 2x2]
    #   Wait: the 2x2 submatrix with row a, col a DELETED has entries:
    #     rows {b,k}, cols {b,k}  => M[b,b], M[b,k], M[k,b], M[k,k]
    #     None of these are eps^a_b (which is at row a) or eps^b_a (which is at col a but row b --
    #     actually row b is in the submatrix, col a is deleted).
    #     eps^b_a is at (row b, col a) — col a is deleted. So it's NOT in the submatrix.
    #   So cof(M)^a_a = M_b * M_k + 0, unchanged.
    #
    #   But what about cof(M)^b_b?
    #   cof(M)^b_b = det of submatrix with row b, col b deleted = rows {a,k}, cols {a,k}
    #   M[a,a] = M_a, M[a,k] = 0, M[k,a] = 0, M[k,k] = M_k
    #   Here eps^a_b is at (row a, col b) — col b is DELETED. Not in submatrix.
    #   eps^b_a is at (row b, col a) — row b is DELETED. Not in submatrix.
    #   => cof^b_b = M_a * M_k, unchanged.
    #
    #   So diagonal F-terms don't contribute cross terms. The cross term H_{(a,b),(b,a)} = 0
    #   analytically at the seesaw vacuum. This is confirmed by the block-diagonal structure.
    #
    #   But wait: there's also the coupling between M^a_b and M^b_a through the
    #   F_{M^a_b} itself at SECOND order. Let me be more careful.
    #   F_{M^a_b} at epsilon^a_b != 0 and epsilon^b_a != 0:
    #   cof(M)^b_a = det of submatrix (rows {a,k}, cols {b,k}):
    #     M[a,b] = eps^a_b, M[a,k] = 0, M[k,b] = 0, M[k,k] = M_k
    #     => cof^b_a = eps^a_b * M_k - 0 = eps^a_b * M_k
    #   Independent of eps^b_a. So F^a_b = X * eps^a_b * M_k.
    #
    #   For other off-diagonal couplings (e.g., eps^a_k and eps^k_a):
    #   F_{M^a_b} = X * cof^b_a where cof^b_a involves rows {a,k}, cols {b,k}.
    #   If eps^k_b != 0 (entry at row k, col b):
    #     M[k,b] = eps^k_b (NOW in the submatrix at position [1,0])
    #     => cof^b_a = eps^a_b * M_k - M[a,k] * eps^k_b = eps^a_b * M_k - 0 = eps^a_b * M_k
    #   Still no coupling. What if eps^a_k != 0?
    #     M[a,k] = eps^a_k (entry at row a, col k, in the submatrix at [0,1])
    #     => cof^b_a = eps^a_b * M_k - eps^a_k * M[k,b] = eps^a_b * M_k - eps^a_k * 0
    #   No coupling either (since M[k,b] = 0 at leading order).
    #
    # Summary: to leading order in epsilon, the 6 off-diagonal sectors
    # DECOUPLE from each other in the quadratic action! The 6x6 Hessian
    # should be (approximately) diagonal.
    #
    # Each diagonal entry:
    #   H_{(a,b),(a,b)} = f_pi^2 + 2 * X_vev^2 * M_k^2

    H_bb = f_pi**2 + 2.0 * X_vev**2 * M_k**2  # same as H_aa by symmetry
    H_ab = 0.0  # cross term is zero

    return H_aa, H_ab, H_bb


print("  Analytic diagonal mass-squared for each off-diagonal sector:")
print(f"  {'Sector':>6}  {'k':>4}  {'M_k (MeV)':>14}  {'f_pi^2':>14}  {'2X^2 M_k^2':>14}  {'m^2_eff':>14}  Status")
print("  " + "-" * 90)

analytic_m2 = {}
for idx, (a, b) in enumerate(od_pairs):
    k = third(a, b)
    M_k = M_vev[k]
    soft_term = f_pi**2
    X_term = 2.0 * X_vev**2 * M_k**2
    m2_eff = soft_term + X_term
    H_aa, H_ab, H_bb = analytic_block_hessian(a, b)
    analytic_m2[(a, b)] = m2_eff
    status = "TACHYONIC" if m2_eff < 0 else "stable"
    print(f"  {od_labels[idx]:>6}  {labels[k]:>4}  {M_k:14.4f}  {soft_term:14.4e}  "
          f"{X_term:14.4e}  {m2_eff:14.4e}  {status}")

print()
print(f"  Note: X_vev = {X_vev:.6e}, X_vev^2 = {X_vev**2:.6e}")
print()

# Cross check with numerical Hessian diagonal
print("  Cross-check: numerical vs analytic diagonal entries:")
print(f"  {'Sector':>6}  {'Numerical':>14}  {'Analytic':>14}  {'Ratio':>10}")
for idx, (a, b) in enumerate(od_pairs):
    num_val = H6[idx, idx]
    ana_val = analytic_m2[(a, b)]
    ratio = num_val / ana_val if ana_val != 0 else float('inf')
    print(f"  {od_labels[idx]:>6}  {num_val:14.4e}  {ana_val:14.4e}  {ratio:10.6f}")
print()

# ============================================================
# TACHYONIC DIRECTIONS: FIND NONZERO VEVS
# ============================================================
print("=" * 70)
print("TACHYONIC MODES: FINDING THE POTENTIAL MINIMUM")
print("=" * 70)
print()
print("  For each tachyonic off-diagonal sector, V(epsilon) looks like:")
print("  V(epsilon) = (1/2) m^2_eff * epsilon^2 + (1/4) lambda_4 * epsilon^4 + ...")
print()
print("  The quartic coefficient comes from higher-order terms in |F_M|^2 and |F_X|^2.")
print()
print("  For the off-diagonal mode epsilon^a_b, the quartic V_4:")
print("  The leading quartic from |F_{M^a_b}|^2:")
print("    F_{M^a_b} = X_vev * cof^b_a = X_vev * [eps^a_b * M_k + O(eps^3)]")
print("    |F_{M^a_b}|^2 = X_vev^2 * M_k^2 * eps^2 + O(eps^4)")
print("  No quartic from F_{M^a_b} at this order (cofactor is linear in eps).")
print()
print("  Quartic terms arise from |F_{M^i_j}|^2 where i=j (diagonal F-terms):")
print("    F_{M^a_a} = m_a + X_vev * cof^a_a")
print("    cof^a_a (rows/cols {b,k}): depends on eps^a_b only at O(eps^2)?")
print("  Actually at O(eps^2), the diagonal cofactors gain corrections.")
print()
print("  Compute V(epsilon) along the tachyonic direction numerically.")
print()

# For each tachyonic sector, scan V along the direction
# to find the minimum (quartic potential minimum)

vev_od = np.zeros((3, 3))  # off-diagonal VEVs
V_shifts = {}

for idx, (a, b) in enumerate(od_pairs):
    m2 = analytic_m2[(a, b)]
    if m2 >= 0:
        print(f"  Sector {od_labels[idx]}: m^2 = {m2:.4e} > 0, no tachyon, eps = 0.")
        continue

    print(f"  Sector {od_labels[idx]}: m^2 = {m2:.4e} < 0 (TACHYONIC)")

    # Scan V(epsilon) along this direction (all other off-diag = 0)
    def V_1d(eps_val):
        ev = np.zeros(6)
        ev[idx] = eps_val
        return V_full(ev)

    # Estimate the scale of the minimum:
    # At the minimum: dV/deps = m^2 eps + lambda_4 eps^3 = 0
    # => eps_min ~ sqrt(-m^2 / lambda_4)
    # Estimate lambda_4 numerically
    h4 = abs(m2)**0.5 * 0.01
    eps_grid = np.linspace(-5*h4, 5*h4, 200)
    V_grid = np.array([V_1d(e) for e in eps_grid])

    # Also scan larger range to find the true minimum
    eps_max = abs(M_vev[a]) * 0.5  # up to 50% of diagonal
    eps_wide = np.linspace(-eps_max, eps_max, 1000)
    V_wide = np.array([V_1d(e) for e in eps_wide])
    i_min_wide = np.argmin(V_wide)

    eps_guess = eps_wide[i_min_wide]

    # Refine with minimize_scalar
    bracket_lo = eps_wide[max(0, i_min_wide-10)]
    bracket_hi = eps_wide[min(999, i_min_wide+10)]
    if bracket_lo > bracket_hi:
        bracket_lo, bracket_hi = bracket_hi, bracket_lo

    # find where V'=0 more precisely
    from scipy.optimize import minimize_scalar as ms
    res = ms(V_1d, bounds=(bracket_lo - eps_max*0.1, bracket_hi + eps_max*0.1), method='bounded',
             options={'xatol': 1e-6})

    eps_min = res.x
    V_min = res.fun
    V_diag_here = V_1d(0.0)
    delta_V = V_min - V_diag_here

    print(f"    Estimated VEV: epsilon_min = {eps_min:.6e} MeV")
    print(f"    V(0)    = {V_diag_here:.6e},  V(eps_min) = {V_min:.6e}")
    print(f"    Delta_V = {delta_V:.6e} MeV^2")

    # Estimate quartic coupling from the 4th derivative
    delta_4 = 1e-2 * abs(eps_min) if abs(eps_min) > 1e-10 else 1.0
    V4_pts = [V_1d(e) for e in [-2*delta_4, -delta_4, 0, delta_4, 2*delta_4]]
    # 4th derivative: (V(2h) - 4V(h) + 6V(0) - 4V(-h) + V(-2h)) / h^4
    lam_4_num = (V4_pts[4] - 4*V4_pts[3] + 6*V4_pts[2] - 4*V4_pts[1] + V4_pts[0]) / (24 * delta_4**4)
    # lam_4_num is d^4V/deps^4 / 4! = quartic coefficient
    # Actually the quartic potential is V = (1/2) m^2 eps^2 + (1/24) d^4V/deps^4 * eps^4
    # The 4th derivative / 4! gives the coefficient in 4th-order expansion.
    # Minimum condition: m^2 + (1/6) d^4V/deps^4 * eps^2 = 0
    # => eps^2 = -6 m^2 / (d^4V/deps^4)
    d4V = (V4_pts[4] - 4*V4_pts[3] + 6*V4_pts[2] - 4*V4_pts[1] + V4_pts[0]) / (delta_4**4)
    if d4V > 0:
        eps_analytic = np.sqrt(-6 * m2 / d4V)
    else:
        eps_analytic = float('nan')

    print(f"    d^4V/deps^4 = {d4V:.6e}")
    if not np.isnan(eps_analytic):
        print(f"    Analytic estimate: |eps_min| ~ sqrt(-6 m^2 / d^4V) = {eps_analytic:.6e} MeV")

    # Store the VEV (use positive root)
    vev_od[a, b] = eps_min
    V_shifts[(a, b)] = delta_V
    print()

# ============================================================
# FULL MESON MATRIX WITH OFF-DIAGONAL VEVS
# ============================================================
print("=" * 70)
print("FULL MESON MATRIX WITH OFF-DIAGONAL VEVS")
print("=" * 70)
print()

M_full = np.diag(M_vev.copy())
for (a, b), vev in zip(od_pairs, vev_od[np.array(od_pairs)[:,0], np.array(od_pairs)[:,1]].tolist()):
    M_full[a, b] = vev

# Re-extract from vev_od
M_full = np.diag(M_vev.copy())
for (a, b) in od_pairs:
    M_full[a, b] = vev_od[a, b]

print("  Meson matrix M (MeV):")
print(f"  {'':>4}" + "".join(f"  {labels[j]:>14}" for j in range(3)))
for i in range(3):
    print(f"  {labels[i]:>4}" + "".join(f"  {M_full[i,j]:14.4e}" for j in range(3)))
print()

print("  Off-diagonal/diagonal ratios:")
for (a, b) in od_pairs:
    if M_full[a, a] != 0:
        ratio = abs(M_full[a, b]) / abs(M_full[a, a])
        print(f"  |M^{labels[a]}_{labels[b]}| / M^{labels[a]}_{labels[a]} = {ratio:.6e}")
print()

# ============================================================
# SVD AND CKM ANGLE EXTRACTION
# ============================================================
print("=" * 70)
print("SVD AND CKM MIXING ANGLES")
print("=" * 70)
print()

U_svd, sigma_svd, Vt_svd = svd(M_full)
V_svd = Vt_svd.T

print("  SVD: M = U @ diag(sigma) @ V^T")
print(f"  Singular values: {sigma_svd[0]:.6f},  {sigma_svd[1]:.6f},  {sigma_svd[2]:.6f} MeV")
print(f"  (Diagonal seesaw: {M_vev[0]:.6f},  {M_vev[1]:.6f},  {M_vev[2]:.6f}) MeV")
print()

print("  U (left singular vectors):")
for i in range(3):
    print(f"    [{U_svd[i,0]:+.8f}  {U_svd[i,1]:+.8f}  {U_svd[i,2]:+.8f}]")
print()
print("  V (right singular vectors):")
for i in range(3):
    print(f"    [{V_svd[i,0]:+.8f}  {V_svd[i,1]:+.8f}  {V_svd[i,2]:+.8f}]")
print()

# CKM-like mixing matrix V_CKM = U^† V (for down-type quarks)
# Convention: if M is the quark mass matrix diagonalized by U^† M V = diag(masses),
# the CKM matrix is V_CKM = U_up^† U_down.
# Here we have a single meson matrix M = M_{ij} coupling flavor i to j.
# The mixing matrix that appears in charged currents is V_CKM = U_L^dag V_R.
# For a Hermitian-like matrix, the SVD gives U_L = U, U_R = V.
# The CKM mixing is V_CKM = U^T V (for real matrices).

V_CKM = U_svd.T @ V_svd

# Ensure det > 0 (proper rotation convention)
if np.linalg.det(V_CKM) < 0:
    V_CKM[:, -1] *= -1

print("  CKM mixing matrix V_CKM = U^T V:")
for i in range(3):
    print(f"    [{V_CKM[i,0]:+.8f}  {V_CKM[i,1]:+.8f}  {V_CKM[i,2]:+.8f}]")
print()

# Extract angles using standard PDG parametrization
def extract_ckm_angles(V):
    """Extract theta_12, theta_23, theta_13 from a unitary matrix."""
    s13 = np.clip(abs(V[0, 2]), 0.0, 1.0)
    theta_13 = np.arcsin(s13)
    c13 = np.cos(theta_13)
    if c13 > 1e-10:
        s12 = np.clip(abs(V[0, 1]) / c13, 0.0, 1.0)
        s23 = np.clip(abs(V[1, 2]) / c13, 0.0, 1.0)
    else:
        s12 = 0.0
        s23 = 0.0
    theta_12 = np.arcsin(s12)
    theta_23 = np.arcsin(s23)
    return theta_12, theta_23, theta_13

theta_12, theta_23, theta_13 = extract_ckm_angles(V_CKM)

# PDG values
theta_12_PDG = np.radians(13.04)
theta_23_PDG = np.radians(2.38)
theta_13_PDG = np.radians(0.201)

print("  CKM angles (PDG standard parametrization):")
print(f"  {'Angle':>10}  {'Computed':>12}  {'PDG':>12}  {'Ratio':>10}")
print("  " + "-" * 50)
for name, val, pdg in [("theta_12", theta_12, theta_12_PDG),
                        ("theta_23", theta_23, theta_23_PDG),
                        ("theta_13", theta_13, theta_13_PDG)]:
    ratio = np.degrees(val) / np.degrees(pdg) if pdg != 0 else float('inf')
    print(f"  {name:>10}  {np.degrees(val):>10.4f}°  {np.degrees(pdg):>10.4f}°  {ratio:10.4f}")
print()

# Also report CKM matrix elements
Vus_comp = abs(V_CKM[0, 1])
Vcb_comp = abs(V_CKM[1, 2])
Vub_comp = abs(V_CKM[0, 2])
Vus_PDG = 0.2243
Vcb_PDG = 0.0422
Vub_PDG = 0.00394

print("  CKM matrix elements:")
print(f"  {'Element':>8}  {'Computed':>12}  {'PDG':>12}  {'Ratio':>10}")
print("  " + "-" * 48)
for name, val, pdg in [("Vus", Vus_comp, Vus_PDG),
                        ("Vcb", Vcb_comp, Vcb_PDG),
                        ("Vub", Vub_comp, Vub_PDG)]:
    ratio = val / pdg if pdg != 0 else float('inf')
    print(f"  {name:>8}  {val:>12.6f}  {pdg:>12.6f}  {ratio:10.4f}")
print()

# ============================================================
# ALTERNATIVE: PERTURBATIVE ANALYTIC ESTIMATE
# ============================================================
print("=" * 70)
print("PERTURBATIVE ANALYTIC ESTIMATE OF MIXING ANGLES")
print("=" * 70)
print()
print("  When the off-diagonal VEVs epsilon^a_b are small relative to the")
print("  diagonal entries, the mixing angle in the (a,b) sector is:")
print("    theta_{ab} ≈ |epsilon^a_b + epsilon^b_a| / |M_a - M_b|")
print()

for a, b in [(0, 1), (0, 2), (1, 2)]:
    eps_ab = vev_od[a, b]
    eps_ba = vev_od[b, a]
    dM = abs(M_vev[a] - M_vev[b])
    theta_approx = abs(eps_ab + eps_ba) / dM if dM > 0 else 0.0
    print(f"  theta_{labels[a]}{labels[b]}:")
    print(f"    eps^{labels[a]}_{labels[b]} = {eps_ab:.4e},  eps^{labels[b]}_{labels[a]} = {eps_ba:.4e}")
    print(f"    |M_{labels[a]} - M_{labels[b]}| = {dM:.4f}")
    print(f"    theta_approx = {np.degrees(theta_approx):.6f} deg")
    print()

# ============================================================
# POTENTIAL ENERGY BOOKKEEPING
# ============================================================
print("=" * 70)
print("POTENTIAL ENERGY AT THE OFF-DIAGONAL MINIMUM")
print("=" * 70)
print()

eps_vec_best = np.array([vev_od[a, b] for (a, b) in od_pairs])
V_od_min = V_full(eps_vec_best)
V_diag_here = V_full(np.zeros(6))

print(f"  V(diagonal seesaw vacuum) = {V_diag_here:.6e} MeV^2")
print(f"  V(with off-diagonal VEVs) = {V_od_min:.6e} MeV^2")
print(f"  Energy gain = {V_od_min - V_diag_here:.6e} MeV^2")
print()

# Decompose the potential at the minimum
M_check = np.diag(M_vev.copy())
for (a, b) in od_pairs:
    M_check[a, b] = vev_od[a, b]

def V_decomposed(M):
    """Return breakdown of V into F-term contributions."""
    cof_M = np.zeros((3, 3))
    cof_M[0,0] = M[1,1]*M[2,2] - M[1,2]*M[2,1]
    cof_M[1,0] = -(M[0,1]*M[2,2] - M[0,2]*M[2,1])
    cof_M[2,0] = M[0,1]*M[1,2] - M[0,2]*M[1,1]
    cof_M[0,1] = -(M[1,0]*M[2,2] - M[1,2]*M[2,0])
    cof_M[1,1] = M[0,0]*M[2,2] - M[0,2]*M[2,0]
    cof_M[2,1] = -(M[0,0]*M[1,2] - M[0,2]*M[1,0])
    cof_M[0,2] = M[1,0]*M[2,1] - M[1,1]*M[2,0]
    cof_M[1,2] = -(M[0,0]*M[2,1] - M[0,1]*M[2,0])
    cof_M[2,2] = M[0,0]*M[1,1] - M[0,1]*M[1,0]

    F_M = np.zeros((3,3))
    for i in range(3):
        for j in range(3):
            F_M[i,j] = X_vev * cof_M[j,i]
            if i == j:
                F_M[i,j] += masses[i]

    detM = (M[0,0]*(M[1,1]*M[2,2]-M[1,2]*M[2,1])
           -M[0,1]*(M[1,0]*M[2,2]-M[1,2]*M[2,0])
           +M[0,2]*(M[1,0]*M[2,1]-M[1,1]*M[2,0]))
    F_X = detM - Lambda6
    F_Hu = y_c * M[1,1]
    F_Hd = y_b * M[2,2]

    V_FM_diag = sum(F_M[i,i]**2 for i in range(3))
    V_FM_od = sum(F_M[i,j]**2 for i in range(3) for j in range(3) if i != j)
    V_FX = F_X**2
    V_FH = F_Hu**2 + F_Hd**2
    V_soft = f_pi**2 * np.sum(M**2)
    return {
        'F_M_diag': V_FM_diag,
        'F_M_od': V_FM_od,
        'F_X': V_FX,
        'F_H': V_FH,
        'soft': V_soft,
        'total': V_FM_diag + V_FM_od + V_FX + V_FH + V_soft,
        'F_M_diag_vals': [F_M[i,i] for i in range(3)],
        'F_X_val': F_X,
    }

breakdown = V_decomposed(M_check)
print("  V decomposition at off-diagonal minimum:")
print(f"    |F_M(diag)|^2  = {breakdown['F_M_diag']:.6e}")
print(f"    |F_M(offdiag)|^2 = {breakdown['F_M_od']:.6e}")
print(f"    |F_X|^2        = {breakdown['F_X']:.6e}")
print(f"    |F_H|^2        = {breakdown['F_H']:.6e}")
print(f"    f_pi^2 Tr(M^2) = {breakdown['soft']:.6e}")
print(f"    Total V        = {breakdown['total']:.6e}")
print()
print(f"  F_X = det(M) - Lambda^6 = {breakdown['F_X_val']:.6e}")
print(f"    (non-zero because M_ij != 0 changes det M)")
print()

# ============================================================
# STRUCTURE OF THE HESSIAN: WHY SOME SECTORS ARE TACHYONIC
# ============================================================
print("=" * 70)
print("STRUCTURE ANALYSIS: WHY CERTAIN SECTORS ARE TACHYONIC")
print("=" * 70)
print()
print("  m^2_eff = f_pi^2 + 2 X_vev^2 M_k^2")
print(f"  f_pi^2 = {f_pi**2:.4e} MeV^2")
print(f"  X_vev = {X_vev:.6e} MeV^(-4)")
print()
print("  The term 2 X_vev^2 M_k^2 is ALWAYS positive.")
print("  Therefore all sectors have m^2_eff > 0 at the seesaw vacuum!")
print()
print("  *** STABILITY RESULT: The diagonal seesaw vacuum is STABLE against ***")
print("  *** all off-diagonal meson perturbations when B = B~ = 0 enforced. ***")
print()
print("  This is in contrast to the statement in the problem that tachyonic")
print("  eigenvalues arise from W_{IJK} F_I terms. Let us investigate those:")
print()

# ============================================================
# W_{IJK} F_I TERMS (SUGRA-LIKE MASS MATRIX)
# ============================================================
print("  The full SUSY scalar mass matrix including D-terms and Kahler corrections")
print("  has the schematic form:")
print("  m^2_{IJ} = W_{IKL} W*^{JKL} - R_{IJ} + (D-term contributions)")
print()
print("  In global SUSY with Kahler K = phi_I phi_I^*:")
print("  (mass^2)_{IJ} = W_{IK} W*_{JK} [fermion mass] is not the scalar mass.")
print("  The scalar mass matrix is the Hessian of V = |F|^2 + f_pi^2 Tr(M^dag M).")
print()
print("  The W_{IJK} F_I term arises from expanding |F_I|^2 to second order:")
print("  V = |W_I|^2 = |W_I^(0) + W_{IJ} delta phi_J + (1/2) W_{IJK} delta phi_J delta phi_K + ...|^2")
print()
print("  The Hessian: H_{JK} = W_{IJ} W*_{IK} + W_I* W_{IJK} + ...")
print("  The W_I^* W_{IJK} term can contribute negative entries if W_I^* < 0.")
print()
print("  At the seesaw vacuum: F_{M^i_i} = m_i + X M_j M_k = 0  (all zero by construction).")
print("  So W_I^* = F_I^* = 0 for the diagonal entries => no W_{IJK} F_I contribution!")
print()
print("  => At the EXACT seesaw vacuum (where all F-terms vanish), there are")
print("     NO tachyonic directions. The Hessian is positive definite in the")
print("     off-diagonal meson sector.")
print()
print("  CAVEAT: The seesaw vacuum does NOT satisfy F_{H_u} = y_c M_d = 0")
print("  and F_{H_d} = y_b M_s = 0 at H = 0.")
print("  These F-terms source W_{IJK} F_I contributions through the coupling:")
print("  W = ... + y_c H_u M^d_d + y_b H_d M^s_s + lambda X H_u H_d")
print()
print("  Specifically: W_{H_u, M^d_d, X} = 0  [no 3-way coupling in simple W]")
print("  But W_{M^d_d, M^d_d, H_u} = 0, W_{X, M^d_d, H_u} = 0 [by linearity in H_u]")
print()
print("  The Yukawa couplings y_c and y_b contribute to F_Hu and F_Hd,")
print("  not to off-diagonal F-terms. So no tachyonic instability at the seesaw")
print("  vacuum (with B=B~=0 enforced and H_u=H_d=0).")
print()

# ============================================================
# WHAT IF WE INCLUDE THE YUKAWA COUPLING TO H?
# ============================================================
print("=" * 70)
print("CROSS-SECTOR INSTABILITY: HIGGS-MESON COUPLING")
print("=" * 70)
print()
print("  Including the H_u, H_d degrees of freedom in the Hessian,")
print("  there can be cross-sector mixing between M^d_d and H_u (via y_c)")
print("  and between M^s_s and H_d (via y_b).")
print()
print("  This does NOT directly create off-diagonal meson tachyons.")
print("  Rather, it causes M^d_d and M^s_s (the diagonal entries) to")
print("  relax from their seesaw values, breaking the seesaw condition.")
print()

# Compute optimal H_u, H_d given M = seesaw vacuum
# F_{H_u} = y_c M_d + lambda X H_d = 0  => H_u minimized when:
# d|F_{H_u}|^2/dH_u + d|F_{H_d}|^2/dH_u = 0
# F_{H_u} = y_c M_d (at H_d=0),  so |F_{H_u}|^2 = y_c^2 M_d^2
# This is independent of H_u! (F_Hu = y_c M_d depends on M_d, not H_u)
# Wait -- F_{H_u} = dW/dH_u = y_c M^d_d + lambda X H_d
# This does not depend on H_u at all -- it depends on M^d_d and H_d.
# So V contains terms |y_c M^d_d + lambda X H_d|^2 (from |F_{H_u}|^2).
# Minimizing over H_d: d/dH_d |lambda X H_d + y_c M^d_d|^2 = 0
#   => lambda X * 2(lambda X H_d + y_c M^d_d) = 0
#   => H_d = -y_c M^d_d / (lambda X)

H_d_opt = -y_c * M_vev[1] / (lam * X_vev)
H_u_opt = -y_b * M_vev[2] / (lam * X_vev)

print(f"  At the seesaw vacuum, minimizing V over H_u and H_d:")
print(f"    F_{{H_u}} = y_c M_d + lambda X H_d = 0  => H_d = -y_c M_d / (lambda X)")
print(f"    F_{{H_d}} = y_b M_s + lambda X H_u = 0  => H_u = -y_b M_s / (lambda X)")
print()
print(f"    lambda = {lam:.4f},  X_vev = {X_vev:.6e}")
print(f"    H_d (optimal) = {H_d_opt:.6e} MeV")
print(f"    H_u (optimal) = {H_u_opt:.6e} MeV")
print()

# These are the Higgs VEVs at the QCD scale. Compare to EW scale:
print(f"  Comparison with EW scale:")
print(f"    v/sqrt(2) = {v/np.sqrt(2):.1f} MeV = {v/np.sqrt(2)/1000:.1f} GeV")
print(f"    |H_d_opt| = {abs(H_d_opt):.4e} MeV  (ratio to v/sqrt(2): {abs(H_d_opt)/(v/np.sqrt(2)):.4e})")
print(f"    |H_u_opt| = {abs(H_u_opt):.4e} MeV  (ratio to v/sqrt(2): {abs(H_u_opt)/(v/np.sqrt(2)):.4e})")
print()

# ============================================================
# HESSIAN INCLUDING HIGGS DEGREES OF FREEDOM
# ============================================================
print("=" * 70)
print("FULL HESSIAN INCLUDING HIGGS DOF AT SEESAW VACUUM")
print("=" * 70)
print()

def V_full_with_H(params):
    """
    V(epsilon_od[6], H_u, H_d) at the seesaw vacuum.
    params = [eps_ud, eps_us, eps_du, eps_ds, eps_su, eps_sd, H_u, H_d]
    """
    epsilon_vec = params[:6]
    Hu = params[6]
    Hd = params[7]

    M = np.diag(M_vev.copy())
    M[0, 1] = epsilon_vec[0]
    M[0, 2] = epsilon_vec[1]
    M[1, 0] = epsilon_vec[2]
    M[1, 2] = epsilon_vec[3]
    M[2, 0] = epsilon_vec[4]
    M[2, 1] = epsilon_vec[5]

    cof_M = np.zeros((3, 3))
    cof_M[0,0] = M[1,1]*M[2,2] - M[1,2]*M[2,1]
    cof_M[1,0] = -(M[0,1]*M[2,2] - M[0,2]*M[2,1])
    cof_M[2,0] = M[0,1]*M[1,2] - M[0,2]*M[1,1]
    cof_M[0,1] = -(M[1,0]*M[2,2] - M[1,2]*M[2,0])
    cof_M[1,1] = M[0,0]*M[2,2] - M[0,2]*M[2,0]
    cof_M[2,1] = -(M[0,0]*M[1,2] - M[0,2]*M[1,0])
    cof_M[0,2] = M[1,0]*M[2,1] - M[1,1]*M[2,0]
    cof_M[1,2] = -(M[0,0]*M[2,1] - M[0,1]*M[2,0])
    cof_M[2,2] = M[0,0]*M[1,1] - M[0,1]*M[1,0]

    F_M = np.zeros((3, 3))
    for i in range(3):
        for j in range(3):
            F_M[i, j] = X_vev * cof_M[j, i]
            if i == j:
                F_M[i, j] += masses[i]
    F_M[1, 1] += y_c * Hu   # y_c H_u for M^d_d
    F_M[2, 2] += y_b * Hd   # y_b H_d for M^s_s

    detM = (M[0,0]*(M[1,1]*M[2,2]-M[1,2]*M[2,1])
           -M[0,1]*(M[1,0]*M[2,2]-M[1,2]*M[2,0])
           +M[0,2]*(M[1,0]*M[2,1]-M[1,1]*M[2,0]))
    F_X = detM - Lambda6 + lam * Hu * Hd  # W = X(det M - Lambda^6) + lambda X Hu Hd ?
    # Actually the lambda term is in W as: lambda X H_u H_d
    # => F_X = dW/dX = det M - Lambda^6 + lambda H_u H_d
    # F_{H_u} = dW/dH_u = y_c M^d_d + lambda X H_d
    # F_{H_d} = dW/dH_d = y_b M^s_s + lambda X H_u

    # Recompute F_X correctly
    F_X = detM - Lambda6 + lam * Hu * Hd
    F_Hu = y_c * M[1, 1] + lam * X_vev * Hd
    F_Hd = y_b * M[2, 2] + lam * X_vev * Hu

    V = (np.sum(F_M**2) + F_X**2 + F_Hu**2 + F_Hd**2
         + f_pi**2 * np.sum(M**2))
    return V


def numerical_hessian_8x8(h=1e-3):
    x0 = np.zeros(8)
    H = np.zeros((8, 8))
    for i in range(8):
        for j in range(i, 8):
            xpp = x0.copy(); xpp[i] += h; xpp[j] += h
            xpm = x0.copy(); xpm[i] += h; xpm[j] -= h
            xmp = x0.copy(); xmp[i] -= h; xmp[j] += h
            xmm = x0.copy(); xmm[i] -= h; xmm[j] -= h
            H[i, j] = (V_full_with_H(xpp) - V_full_with_H(xpm)
                      - V_full_with_H(xmp) + V_full_with_H(xmm)) / (4*h*h)
            H[j, i] = H[i, j]
    return H

h8 = min(M_vev) * 1e-4
print(f"  Computing 8x8 Hessian (h = {h8:.4e})...")
H8 = numerical_hessian_8x8(h=h8)

labels_8 = od_labels + ['H_u', 'H_d']
print()
print("  8x8 Hessian (off-diagonal mesons + Higgs):")
print(f"  {'':>6}" + "".join(f"{l:>12}" for l in labels_8))
for i in range(8):
    row = f"  {labels_8[i]:>6}" + "".join(f"{H8[i,j]:12.3e}" for j in range(8))
    print(row)
print()

evals_8 = np.linalg.eigvalsh(H8)
print("  Eigenvalues of 8x8 Hessian (MeV^2):")
for k, ev in enumerate(sorted(evals_8)):
    status = "TACHYONIC" if ev < 0 else "stable"
    print(f"    lambda_{k+1} = {ev:14.6e}  [{status}]")
print()

n_tachyonic_8 = np.sum(evals_8 < 0)
print(f"  Number of tachyonic directions: {n_tachyonic_8}")
print()

# ============================================================
# SUMMARY
# ============================================================
print("=" * 70)
print("SUMMARY")
print("=" * 70)
print()
print(f"  Physical parameters:")
print(f"    m_u = {m_u} MeV,  m_d = {m_d} MeV,  m_s = {m_s} MeV")
print(f"    Lambda = {Lambda} MeV,  f_pi = {f_pi} MeV")
print(f"    y_c = {y_c:.6e},  y_b = {y_b:.6e},  lambda = {lam}")
print()
print(f"  Seesaw vacuum:")
print(f"    M_u = {M_vev[0]:.4f},  M_d = {M_vev[1]:.4f},  M_s = {M_vev[2]:.4f} MeV")
print(f"    X_vev = {X_vev:.6e} MeV^-4")
print()
print(f"  Off-diagonal Hessian (B = B~ = 0, H_u = H_d = 0):")
print(f"    m^2 = f_pi^2 + 2 X_vev^2 M_k^2 = f_pi^2 + [positive] > 0")
for idx, (a, b) in enumerate(od_pairs):
    m2 = analytic_m2[(a, b)]
    status = "TACHYONIC" if m2 < 0 else "stable"
    print(f"    {od_labels[idx]}: m^2 = {m2:.4e} MeV^2  [{status}]")
print()
print(f"  All off-diagonal sectors are STABLE (m^2 > 0).")
print(f"  The tachyonic sectors are suppressed by f_pi^2 overwhelming the")
print(f"  destabilizing terms (which are also positive here).")
print()
print(f"  8x8 Hessian (including H_u and H_d):")
print(f"    Tachyonic directions: {n_tachyonic_8}")
print()
print(f"  Meson matrix M with off-diagonal VEVs (all zero since no tachyons):")
print(f"    M = diag({M_vev[0]:.4f}, {M_vev[1]:.4f}, {M_vev[2]:.4f}) MeV")
print()
print(f"  CKM angles at the diagonal vacuum (trivially):")
print(f"    theta_12 = 0.0000 deg  (PDG: 13.04 deg)")
print(f"    theta_23 = 0.0000 deg  (PDG: 2.38 deg)")
print(f"    theta_13 = 0.0000 deg  (PDG: 0.201 deg)")
print()
print(f"  Off-diagonal VEVs from tachyonic directions (if any existed):")
for (a, b) in od_pairs:
    if vev_od[a, b] != 0:
        print(f"    epsilon^{labels[a]}_{labels[b]} = {vev_od[a,b]:.6e} MeV")
print()
print(f"  SVD mixing angles from M_full:")
print(f"    theta_12 = {np.degrees(theta_12):.6f} deg  (PDG: {np.degrees(theta_12_PDG):.2f} deg)")
print(f"    theta_23 = {np.degrees(theta_23):.6f} deg  (PDG: {np.degrees(theta_23_PDG):.2f} deg)")
print(f"    theta_13 = {np.degrees(theta_13):.6f} deg  (PDG: {np.degrees(theta_13_PDG):.3f} deg)")
print()
print("  KEY FINDING: The seesaw vacuum is stable against off-diagonal")
print("  perturbations. The quadratic terms are ALL positive because:")
print("  (1) f_pi^2 > 0 (soft breaking)")
print("  (2) 2 X_vev^2 M_k^2 > 0 (F-term from cofactor, X_vev^2 > 0)")
print("  The problem statement's claim of tachyonic directions requires")
print("  either: (a) a different point in field space, or (b) the Yukawa")
print("  W_{IJK} F_I terms evaluated at a point where F_I != 0.")
print()
print("  Physical picture: at the seesaw vacuum, all diagonal F-terms vanish")
print("  by construction. The W_{IJK} F_I mechanism requires nonzero F-terms.")
print("  The Higgs F-terms F_{H_u} = y_c M_d != 0 and F_{H_d} = y_b M_s != 0")
print("  are the only nonzero F-terms, and they do not couple to off-diagonal")
print("  meson entries at leading order.")

sys.stdout.flush()

# ============================================================
# WRITE MARKDOWN REPORT
# ============================================================
lines = []
def md(s=""): lines.append(s)

md("# Perturbative CKM Mixing from Off-Diagonal Meson VEVs")
md()
md("**Date:** 2026-03-04")
md()
md("## Setup")
md()
md("The scalar potential is")
md()
md("    V = sum_{I} |F_I|^2 + f_pi^2 Tr(M†M)")
md()
md("with B = B̃ = 0 enforced, H_u = H_d = 0 at the QCD-scale analysis.")
md("The superpotential is")
md()
md("    W = sum_i m_i M^i_i + X(det M - Lambda^6)")
md("        + y_c H_u M^d_d + y_b H_d M^s_s + lambda X H_u H_d")
md()
md("The seesaw vacuum is M^i_j = delta^i_j C/m_i with")
md(f"C = Lambda^2 (m_u m_d m_s)^(1/3) = {C:.4f} MeV^2.")
md()
md(f"At this vacuum (with H_u = H_d = 0 and B = B̃ = 0), X satisfies")
md(f"X_vev = -m_u / (M_d M_s) = {X_vev:.6e} MeV^(-4),")
md(f"and all diagonal F_{{M^i_i}} = 0 by construction.")
md()
md("## Diagonal Seesaw Vacuum")
md()
md("| Field | Value (MeV) |")
md("|-------|-------------|")
md(f"| M_u = C/m_u | {M_vev[0]:.4f} |")
md(f"| M_d = C/m_d | {M_vev[1]:.4f} |")
md(f"| M_s = C/m_s | {M_vev[2]:.4f} |")
md()
md("## Hessian for Off-Diagonal Meson Entries")
md()
md("Write M^i_j = delta^i_j M_i + epsilon^i_j and expand V to second order.")
md("For each off-diagonal pair (a,b) with a != b, let k be the third index.")
md()
md("The mass-squared of the off-diagonal mode epsilon^a_b is:")
md()
md("    m^2_(a,b) = f_pi^2 + 2 X_vev^2 M_k^2")
md()
md("Origin of each term:")
md()
md("- **f_pi^2**: from the soft-breaking term f_pi^2 Tr(M†M), always positive.")
md("- **2 X_vev^2 M_k^2**: from |F_{M^a_b}|^2 = |X_vev * cof(M)^b_a|^2.")
md("  At the diagonal vacuum, cof(M)^b_a = epsilon^a_b * M_k (to leading order),")
md("  giving |F_{M^a_b}|^2 = X_vev^2 * M_k^2 * |epsilon^a_b|^2. The factor 2")
md("  comes from both |F_{M^a_b}|^2 and |F_{M^b_a}|^2 each contributing X_vev^2 M_k^2.")
md()
md("The |F_X|^2 and cross-sector terms vanish at the seesaw vacuum because:")
md("(1) d(det M)/d(epsilon^a_b) = 0 at the diagonal vacuum (off-diagonal cofactors vanish),")
md("(2) the 6 off-diagonal sectors decouple from each other at quadratic order.")
md()
md("### Results")
md()
md("| Sector | k | M_k (MeV) | f_pi^2 | 2X^2 M_k^2 | m^2_eff | Status |")
md("|--------|---|-----------|--------|------------|---------|--------|")
for idx, (a, b) in enumerate(od_pairs):
    k = third(a, b)
    M_k = M_vev[k]
    soft = f_pi**2
    Xterm = 2*X_vev**2*M_k**2
    m2 = analytic_m2[(a, b)]
    status = "TACHYONIC" if m2 < 0 else "stable"
    md(f"| {od_labels[idx]} | {labels[k]} | {M_k:.2f} | {soft:.3e} | {Xterm:.3e} | {m2:.3e} | {status} |")

md()
md("**All off-diagonal sectors are stable** (m^2_eff > 0).")
md()
md("## Physical Interpretation")
md()
md("The seesaw vacuum is stable against off-diagonal perturbations.")
md("Both contributions to m^2_eff are positive:")
md()
md("1. The soft term f_pi^2 = (92 MeV)^2 = 8464 MeV^2 provides a positive mass.")
md("2. The F-term contribution 2 X_vev^2 M_k^2 is positive (X_vev^2 > 0 always).")
md()
md("The problem statement's tachyonic directions require a nonzero W_{IJK} F_I")
md("contribution, where F_I refers to a background F-term. At the exact seesaw")
md("vacuum, all diagonal F-terms vanish by construction (F_{M^i_i} = 0).")
md("The only nonzero F-terms are:")
md()
md(f"    F_{{H_u}} = y_c M_d + lambda X H_d = {y_c * M_vev[1]:.4e} MeV  (at H_d = 0)")
md(f"    F_{{H_d}} = y_b M_s + lambda X H_u = {y_b * M_vev[2]:.4e} MeV  (at H_u = 0)")
md()
md("These couple M^d_d and M^s_s to the Higgs sector but do not create off-diagonal")
md("tachyons in the meson sector at leading order.")
md()
md("## Hessian Including Higgs Degrees of Freedom")
md()
md("Extending the analysis to the 8x8 Hessian (6 off-diagonal mesons + H_u + H_d):")
md()
md(f"- Number of tachyonic eigenvalues: {n_tachyonic_8}")
md()
evals_sorted = sorted(evals_8)
md("Eigenvalues (MeV^2):")
md()
md("| # | Eigenvalue (MeV^2) | Status |")
md("|---|-------------------|--------|")
for k, ev in enumerate(evals_sorted):
    status = "tachyonic" if ev < 0 else "stable"
    md(f"| {k+1} | {ev:.4e} | {status} |")
md()
md("## CKM Angles")
md()
md("Since no off-diagonal VEVs are generated, the meson matrix remains diagonal")
md("and the SVD yields no mixing angles. The CKM mixing is identically zero")
md("at the seesaw vacuum with B = B̃ = 0.")
md()
md("| Angle | Computed | PDG |")
md("|-------|---------|-----|")
md(f"| theta_12 (Cabibbo) | {np.degrees(theta_12):.4f} deg | 13.04 deg |")
md(f"| theta_23 | {np.degrees(theta_23):.4f} deg | 2.38 deg |")
md(f"| theta_13 | {np.degrees(theta_13):.4f} deg | 0.201 deg |")
md()
md("## Conclusions")
md()
md(f"1. The seesaw vacuum M_i = C/m_i is a stable critical point of V")
md(f"   against off-diagonal meson perturbations. The quadratic mass")
md(f"   m^2_eff = f_pi^2 + 2 X_vev^2 M_k^2 is positive for all 6 off-diagonal sectors.")
md()
md(f"2. No tachyonic off-diagonal VEVs are generated at the seesaw vacuum")
md(f"   with B = B̃ = 0 and H_u = H_d = 0.")
md()
md(f"3. The CKM mixing angles are zero at this vacuum. The meson matrix")
md(f"   remains diagonal and its SVD produces no mixing.")
md()
md(f"4. Tachyonic off-diagonal modes would require a point in field space")
md(f"   where F-terms are nonzero and couple to the off-diagonal sector.")
md(f"   The seesaw vacuum has all meson F-terms equal to zero by construction,")
md(f"   so the W_{{IJK}} F_I mechanism does not operate there.")
md()
md(f"5. The Higgs F-terms F_{{H_u}} = y_c M_d ≈ {y_c * M_vev[1]:.2e} and")
md(f"   F_{{H_d}} = y_b M_s ≈ {y_b * M_vev[2]:.2e} are nonzero but do not")
md(f"   generate off-diagonal meson tachyons at leading order.")
md()
md(f"6. If the problem's 'known' tachyonic sectors are to be realized,")
md(f"   the vacuum should be displaced from the seesaw point such that")
md(f"   some F-terms acquire nonzero backgrounds. The ds and us sectors")
md(f"   are the natural candidates because their third index is u (with the")
md(f"   smallest M_k = M_u = {M_vev[0]:.1f} MeV), giving the smallest positive")
md(f"   stabilizing term 2 X_vev^2 M_u^2 = {2*X_vev**2*M_vev[0]**2:.4e} MeV^2.")
md(f"   This is {2*X_vev**2*M_vev[0]**2/f_pi**2:.4f}x f_pi^2.")

report = "\n".join(lines)
with open("/home/codexssh/phys3/results/ckm_perturbative.md", "w") as f:
    f.write(report)

print()
print("=" * 70)
print(f"Markdown report written to results/ckm_perturbative.md")
print("=" * 70)
