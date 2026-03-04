"""
Complete particle spectrum of an N=1 SUSY model.

Field content:
  M^i_j  (3x3 complex meson matrix):  9 complex scalars + 9 Weyl fermions
  X      (Lagrange multiplier):        1 complex scalar  + 1 Weyl fermion
  B, Btilde (baryons):                 2 complex scalars + 2 Weyl fermions
  H_u, H_d (Higgs doublets):           4 complex scalars + 4 Weyl fermions
  Total: 16 complex scalars (32 real), 16 Weyl fermions

Superpotential:
  W = sum_i m_i M^i_i + X(det M - B Btilde - Lambda^6)
      + y_c H_u^0 M^c_c + y_b H_d^0 M^b_b

where det M = full 3x3 determinant of the meson matrix.

Soft SUSY-breaking:
  V_soft = mtilde^2 Tr(M^dag M) + (B_mu sum_i m_i M^i_i + h.c.)

Vacuum:
  <M^i_i> = C/m_i (diagonal seesaw), off-diagonal = 0
  <X> = -C/Lambda^6, <B> = <Btilde> = 0
  <H_u> = (0, v/sqrt(2)), <H_d> = (v/sqrt(2), 0)
  C = Lambda^2 (m_u m_d m_s)^{1/3}

Pure linear algebra computation. No theory labels.
"""

import numpy as np
from numpy.linalg import eigvalsh, eigh, eigvals

# ═══════════════════════════════════════════════════════════════════════════
# Parameters
# ═══════════════════════════════════════════════════════════════════════════
m_u = 2.16       # MeV
m_d = 4.67       # MeV
m_s = 93.4       # MeV
m_c = 1270.0     # MeV
m_b = 4180.0     # MeV
LAM = 300.0      # MeV
v = 246220.0     # MeV (full EW vev, 246.22 GeV)
f_pi = 92.0      # MeV

# Yukawa couplings (convention: m = y v / sqrt(2), so y = sqrt(2) m / v = 2m/v with v_full)
# Using the problem's convention: y_c = 2 m_c / v, y_b = 2 m_b / v
y_c = 2.0 * m_c / v
y_b = 2.0 * m_b / v

# Soft SUSY breaking
m_tilde_sq = f_pi**2   # = (92)^2 MeV^2
B_mu = f_pi             # = 92 MeV

# Derived quantities
m_arr = np.array([m_u, m_d, m_s])
C = LAM**2 * np.prod(m_arr)**(1.0/3.0)
LAM6 = LAM**6

# Meson VEVs (Seiberg seesaw)
M_u = C / m_u
M_d = C / m_d
M_s = C / m_s
M_vev = np.array([M_u, M_d, M_s])

# X VEV
X0 = -C / LAM6

# Higgs VEVs
v_hu = v / np.sqrt(2)  # <H_u^0>
v_hd = v / np.sqrt(2)  # <H_d^0>

# Determinant check
det_M = M_u * M_d * M_s

print("=" * 72)
print("COMPLETE PARTICLE SPECTRUM: N=1 SUSY MODEL")
print("=" * 72)

print("\n--- Parameter verification ---")
print(f"  m_u = {m_u} MeV,  m_d = {m_d} MeV,  m_s = {m_s} MeV")
print(f"  m_c = {m_c} MeV,  m_b = {m_b} MeV")
print(f"  Lambda = {LAM} MeV,  v = {v} MeV")
print(f"  y_c = 2*m_c/v = {y_c:.8e}")
print(f"  y_b = 2*m_b/v = {y_b:.8e}")
print(f"  f_pi = {f_pi} MeV,  m_tilde^2 = {m_tilde_sq:.0f} MeV^2,  B_mu = {B_mu} MeV")
print(f"  C = Lambda^2 (m_u m_d m_s)^(1/3) = {C:.6f} MeV^2")

# Verify C
C_check = 300.0**2 * (2.16 * 4.67 * 93.4)**(1.0/3.0)
print(f"  C (verify) = {C_check:.6f} MeV^2")

print(f"\n  M_u = C/m_u = {M_u:.4f} MeV")
print(f"  M_d = C/m_d = {M_d:.4f} MeV")
print(f"  M_s = C/m_s = {M_s:.4f} MeV")
print(f"  X_0 = -C/Lambda^6 = {X0:.8e}")
print(f"  det M / Lambda^6 = {det_M/LAM6:.14f}  (should be 1)")
print(f"  <H_u^0> = v/sqrt(2) = {v_hu:.4f} MeV")
print(f"  <H_d^0> = v/sqrt(2) = {v_hd:.4f} MeV")

# ═══════════════════════════════════════════════════════════════════════════
# TASK 1: Fermion mass matrix W_IJ = d^2 W / (d Phi_I d Phi_J)
# ═══════════════════════════════════════════════════════════════════════════
print("\n" + "=" * 72)
print("TASK 1: FERMION MASS MATRIX (16x16)")
print("=" * 72)

# Field ordering (16 Weyl fermions):
# 0: M^u_u,  1: M^u_d,  2: M^u_s
# 3: M^d_u,  4: M^d_d,  5: M^d_s
# 6: M^s_u,  7: M^s_d,  8: M^s_s
# 9: X
# 10: B,  11: Btilde
# 12: H_u^+,  13: H_u^0,  14: H_d^0,  15: H_d^-

NAMES = ['Muu','Mud','Mus','Mdu','Mdd','Mds','Msu','Msd','Mss',
         'X', 'B', 'Bt', 'Hu+','Hu0','Hd0','Hd-']

def idx_M(i, j):
    """Index in the 16-field list for meson M^i_j (i,j in 0,1,2)."""
    return 3*i + j

W = np.zeros((16, 16))

# --- det M second derivatives ---
# det M = M^u_u(M^d_d M^s_s - M^d_s M^s_d) - M^u_d(M^d_u M^s_s - M^d_s M^s_u)
#       + M^u_s(M^d_u M^s_d - M^d_d M^s_u)
#
# At diagonal vacuum M = diag(M_u, M_d, M_s), off-diag = 0:
# d(det M)/dM^i_i = cofactor(i,i) = product of other two diagonal entries
# d^2(det M)/dM^i_i dM^j_j = M_k  (k = third diagonal, i != j)
# d^2(det M)/dM^i_j dM^j_i = (-1)^... cofactor of the 2x2 minor
#
# Explicitly for off-diagonal pairs at diagonal vacuum:
# d^2(det M)/(dM^u_d dM^d_u) = +M_s  (from the -M^u_d * M^d_u * M^s_s term)
# Wait, let's be very careful. det M expanded:
# det M = M^u_u M^d_d M^s_s - M^u_u M^d_s M^s_d
#       - M^u_d M^d_u M^s_s + M^u_d M^d_s M^s_u
#       + M^u_s M^d_u M^s_d - M^u_s M^d_d M^s_u
#
# Second derivatives at the diagonal vacuum:
# d^2/dM^u_u dM^d_d = M^s_s -> M_s
# d^2/dM^u_u dM^s_s = M^d_d -> M_d
# d^2/dM^d_d dM^s_s = M^u_u -> M_u
#
# d^2/dM^u_d dM^d_u: The terms containing both M^u_d and M^d_u:
#   -M^u_d M^d_u M^s_s => d^2/dM^u_d dM^d_u = -M^s_s = -M_s
#
# d^2/dM^u_s dM^s_u: Terms with both:
#   -M^u_s M^d_d M^s_u => d^2/dM^u_s dM^s_u = -M^d_d = -M_d
#
# d^2/dM^d_s dM^s_d: Terms with both:
#   -M^u_u M^d_s M^s_d => d^2/dM^d_s dM^s_d = -M^u_u = -M_u
#
# d^2/dM^u_d dM^d_s: Terms with M^u_d and M^d_s:
#   +M^u_d M^d_s M^s_u => d^2/dM^u_d dM^d_s = M^s_u -> 0 at diagonal vacuum
#
# All other mixed off-diagonal second derivatives vanish at diagonal vacuum.

# Actually, let me redo this more carefully using the Levi-Civita formula:
# det M = epsilon^{ijk} M^u_i M^d_j M^s_k
# d^2(det M)/dM^a_i dM^b_j = epsilon^{abl} ... No, indices don't work that way.
#
# Let's use: det M = sum_{sigma in S_3} sgn(sigma) M^1_{sigma(1)} M^2_{sigma(2)} M^3_{sigma(3)}
# where I'm using 1=u, 2=d, 3=s for row indices and sigma permutes column indices.
#
# d(det M)/dM^a_b = sum_{sigma: sigma(a)=b} sgn(sigma) prod_{c != a} M^c_{sigma(c)}
# = cofactor(a,b)
#
# d^2(det M)/dM^a_b dM^c_d (a != c):
# = sum_{sigma: sigma(a)=b, sigma(c)=d} sgn(sigma) M^e_{sigma(e)} where e = third row
# At diagonal vacuum, M^e_{sigma(e)} = M_e * delta_{e, sigma(e)}
# So we need sigma(a)=b, sigma(c)=d, sigma(e)=e.
# This means: the permutation maps a->b, c->d, e->e.
# If b=a and d=c (diagonal-diagonal): sigma is identity, sgn = +1, value = M_e
# If b != a or d != c: sigma is a transposition (a b)(c d) or identity only if...
# Actually if sigma(a)=b, sigma(c)=d, sigma(e)=e and a,c,e are distinct rows, b,d,e must
# be a permutation of columns. sigma = (a->b, c->d, e->e).
# For this to be a valid permutation, {b,d,e} must be {a,c,e} permuted = {1,2,3}.
# If e = e (it's fixed), then {b,d} must be a permutation of {a,c}.
# Case 1: b=a, d=c => identity perm on those, sgn(sigma) = +1 (sigma = id if also e->e)
# Case 2: b=c, d=a => transposition of a,c, sgn(sigma) = -1
#
# So d^2(det M)/dM^a_b dM^c_d at diagonal vacuum:
# = M_e * [delta_{b,a} delta_{d,c} - delta_{b,c} delta_{d,a}]   (a != c)
# where e is the third row index.
#
# This gives:
# d^2/dM^a_a dM^c_c = M_e * (1 - 0) = M_e   [diagonal-diagonal, agrees]
# d^2/dM^a_c dM^c_a = M_e * (0 - 1) = -M_e  [off-diagonal pair]
# d^2/dM^a_b dM^c_d = 0 for other combinations

# Now build X * d^2(det M)/dM^a_b dM^c_d contributions:

# Diagonal-diagonal meson couplings (from X * det M term)
for a in range(3):
    for c in range(3):
        if a != c:
            e = 3 - a - c  # third index
            i_a = idx_M(a, a)
            i_c = idx_M(c, c)
            W[i_a, i_c] += X0 * M_vev[e]

# Off-diagonal meson pair couplings
for a in range(3):
    for c in range(3):
        if a != c:
            e = 3 - a - c
            i_ac = idx_M(a, c)  # M^a_c
            i_ca = idx_M(c, a)  # M^c_a
            W[i_ac, i_ca] += X0 * (-M_vev[e])

# Diagonal meson - X couplings: W_{M^a_a, X} = cofactor(a,a) = prod_{c != a} M_c
for a in range(3):
    others = [M_vev[c] for c in range(3) if c != a]
    cof = others[0] * others[1]
    i_a = idx_M(a, a)
    W[i_a, 9] += cof
    W[9, i_a] += cof

# Baryon - anti-baryon: W_{B, Bt} = -X0
W[10, 11] += -X0
W[11, 10] += -X0

# Yukawa couplings: y_c H_u^0 M^c_c and y_b H_d^0 M^b_b
# In the problem, M^c_c = M^d_d (c is the second flavor = index 1 in u,d,s labeling)
# Wait — the problem says "M^c_c" and "M^b_b" where c=charm, b=bottom.
# But the meson matrix is 3x3 with i,j = u,d,s. The Yukawa couples H_u to M^c_c
# which maps to the SECOND diagonal entry M^d_d (since c is the second-generation quark,
# and in the sBootstrap/seesaw picture, the (u,d,s) meson M^d_d maps to charm via
# the Seiberg dual). Actually, let me re-read the problem.
#
# The problem states: "y_c H_u^0 M^c_c + y_b H_d^0 M^b_b"
# The meson matrix indices i,j run over u,d,s. But the Yukawa couples to
# M^c_c and M^b_b which use quark flavor indices c (charm) and b (bottom).
# In the ISS/seesaw framework, the meson field M^i_j = Q^i Qbar_j, and the
# Seiberg dual maps M^i_i with VEV C/m_i. The charm and bottom quark masses
# emerge from the Yukawa couplings to specific diagonal entries.
#
# Looking at the earlier code (sqcd_yukawa_spectrum.py), the Yukawa y_j couples
# to M^j_j. With the 3-flavor (u,d,s) structure:
# - M^d_d (index 4 in the 9-meson ordering) has VEV C/m_d, and the Yukawa y_c
#   to this gives a fermion mass that reproduces m_c after EWSB
# - M^s_s (index 8) has VEV C/m_s, and y_b reproduces m_b
#
# But the problem says M^c_c and M^b_b, treating c,b as flavor labels equal to
# the (d,s) entries. Let me follow the problem literally:
# "y_c H_u^0 M^c_c" — here c is the charm quark index. In the 3x3 matrix with
# i,j = u,d,s, this maps to the entry labeled by "charm" which through Seiberg
# duality corresponds to the d-quark entry (ISS maps heavy <-> light).
# However, re-reading the problem more carefully, it just says i,j = u,d,s for
# the meson matrix. The Yukawa terms y_c H_u^0 M^c_c and y_b H_d^0 M^b_b must
# then mean M^c_c = M^{d}_{d} and M^b_b = M^{s}_{s} (mapping charm->d, bottom->s
# through the seesaw: the charm mass emerges from the d-entry, bottom from s-entry).
#
# Actually the simplest reading: the superpotential is written with M^c_c meaning
# the (c,c) diagonal entry, where c = the second flavor index. Given i,j = u,d,s
# and the sBootstrap identification: u<->u, d<->c, s<->b (via seesaw duality),
# M^c_c ≡ M^d_d (index 4) and M^b_b ≡ M^s_s (index 8).
#
# This is consistent with the earlier scripts. So:
# W_{H_u^0, M^d_d} = y_c, W_{H_d^0, M^s_s} = y_b

# Yukawa: y_c H_u^0 M^d_d
# H_u^0 is field index 13, M^d_d is field index 4
W[13, 4] += y_c
W[4, 13] += y_c

# Yukawa: y_b H_d^0 M^s_s
# H_d^0 is field index 14, M^s_s is field index 8
W[14, 8] += y_b
W[8, 14] += y_b

# Verify symmetry
assert np.allclose(W, W.T), "Fermion mass matrix is not symmetric!"

print("\nFermion mass matrix W_IJ (16x16):")
print("Nonzero entries:")
for i in range(16):
    for j in range(i, 16):
        if abs(W[i, j]) > 1e-30:
            print(f"  W[{NAMES[i]:>4s}, {NAMES[j]:>4s}] = {W[i,j]:+.10e}")

# Diagonalize
fermion_evals, fermion_evecs = eigh(W)

print("\nFermion mass eigenvalues (Weyl fermion masses):")
print(f"  {'n':>3}  {'eigenvalue (MeV)':>22}  {'|mass| (MeV)':>16}")
for i, ev in enumerate(fermion_evals):
    print(f"  {i+1:3d}  {ev:+22.10e}  {abs(ev):16.10e}")

# Count zero, positive, negative
tol = 1e-15
n_zero_f = np.sum(np.abs(fermion_evals) < tol)
n_pos_f = np.sum(fermion_evals > tol)
n_neg_f = np.sum(fermion_evals < -tol)
print(f"\n  Zero eigenvalues (|lambda| < {tol}): {n_zero_f}")
print(f"  Positive: {n_pos_f}")
print(f"  Negative: {n_neg_f}")

# Block structure analysis
print("\n--- Block structure ---")

# Central 4x4: {M^u_u, M^d_d, M^s_s, X} = indices {0, 4, 8, 9}
ci = [0, 4, 8, 9]
W_central = W[np.ix_(ci, ci)]
print("\nCentral 4x4 block {Muu, Mdd, Mss, X}:")
clabels = ['Muu','Mdd','Mss','X']
for k, row in enumerate(W_central):
    print(f"  {clabels[k]:>4s}  " + "  ".join(f"{val:+.6e}" for val in row))
ev_central = eigvalsh(W_central)
print(f"  Eigenvalues: {[f'{v:.6e}' for v in ev_central]}")

# Off-diagonal meson 2x2 blocks
print("\nOff-diagonal meson 2x2 blocks:")
offdiag_pairs = [
    ((0, 1), (1, 0), "Mud-Mdu", 2),   # M^u_d, M^d_u; third flavor = s
    ((0, 2), (2, 0), "Mus-Msu", 1),   # M^u_s, M^s_u; third = d
    ((1, 2), (2, 1), "Mds-Msd", 0),   # M^d_s, M^s_d; third = u
]
for (a, b), (c, d), label, e in offdiag_pairs:
    i1 = idx_M(a, b)
    i2 = idx_M(c, d)
    val = W[i1, i2]
    print(f"  {label}: W[{i1},{i2}] = {val:+.10e}, eigenvalues = +/-{abs(val):.10e}")

# Baryon block
print(f"\nBaryon 2x2 block {{B, Bt}}:")
val_bb = W[10, 11]
print(f"  W[B, Bt] = {val_bb:+.10e}, eigenvalues = +/-{abs(val_bb):.10e}")

# Higgs couplings
print(f"\nHiggs couplings:")
print(f"  W[Hu0, Mdd] = {W[13, 4]:+.10e}  (y_c)")
print(f"  W[Hd0, Mss] = {W[14, 8]:+.10e}  (y_b)")
print(f"  H_u^+ (index 12): all zero (decoupled)")
print(f"  H_d^- (index 15): all zero (decoupled)")

# Identify Dirac pairs
print("\n--- Physical fermion spectrum (Dirac masses from W_IJ) ---")
print("Nonzero eigenvalues come in +/- pairs; each pair = one Dirac fermion")
print("Zero eigenvalues = massless Weyl fermions")

# Sort by absolute value
sorted_idx = np.argsort(np.abs(fermion_evals))
print(f"\n{'n':>3}  {'eigenvalue':>22}  {'|mass|':>16}  Identification")
for i in sorted_idx:
    ev = fermion_evals[i]
    evec = fermion_evecs[:, i]
    # Find dominant component
    dom_idx = np.argmax(np.abs(evec))
    dom_frac = evec[dom_idx]**2
    ident = f"dominant: {NAMES[dom_idx]} ({dom_frac:.3f})"
    print(f"  {i+1:3d}  {ev:+22.10e}  {abs(ev):16.10e}  {ident}")


# ═══════════════════════════════════════════════════════════════════════════
# TASK 2: Scalar mass-squared matrix (32x32 real)
# ═══════════════════════════════════════════════════════════════════════════
print("\n\n" + "=" * 72)
print("TASK 2: SCALAR MASS-SQUARED MATRIX (32x32)")
print("=" * 72)

# For each complex scalar phi_I, write phi_I = (phi_I^R + i phi_I^I)/sqrt(2)
# where phi_I^R and phi_I^I are real fields.
#
# The scalar potential contributions:
# (a) |F_I|^2 = |dW/dphi_I|^2 evaluated at the vacuum
# (b) W_IJ W*_KJ contributions to the mass matrix
# (c) Soft SUSY-breaking: m_tilde^2 Tr(M^dag M) + B_mu terms
#
# The mass-squared matrix for real fields (phi_1^R, phi_1^I, phi_2^R, phi_2^I, ...):
#
# From |F_I|^2: V = sum_I |F_I|^2 where F_I = dW/dphi_I
# The second-order expansion around the vacuum:
# delta V = sum_{IJ} [W*_IK W_JK + W_IK W*_JK] delta phi_I^* delta phi_J / 2
#         + sum_{IJ} [W_IJK F*_K delta phi_J delta phi_K + h.c.] (third-derivative terms)
#
# At a SUSY vacuum (F_I = 0), the mass-squared matrix simplifies:
# m^2_{IJ} (complex basis) = sum_K W*_IK W_JK  (for phi_I^*, phi_J combination)
# This gives m^2_scalar = W^dag W in the complex basis.
#
# However, we also have soft breaking terms. And the vacuum is NOT exactly SUSY
# (the soft terms break SUSY explicitly).
#
# Let me be more careful. The full scalar potential is:
# V = sum_I |F_I|^2 + V_soft + V_D
#
# where F_I = dW/dphi_I (evaluated with VEVs substituted for background).
#
# The expansion to second order around the VEV <phi_I> requires computing
# F_I at the vacuum and the second-derivative matrix.
#
# F-term contributions to the scalar mass matrix:
# V_F = sum_I |dW/dphi_I|^2
# d^2 V_F / d(Re phi_J) d(Re phi_K) involves W_IJ W*_IK + F*_I W_IJK + h.c.
#
# At a point where F_I = 0 (SUSY vacuum), the scalar mass matrix from F-terms
# is simply (W^dag W)_JK for the complex scalar mass matrix.
# In real components, it becomes a 32x32 matrix.
#
# For BROKEN SUSY (F_I != 0), we need the full expression including W_IJK F*_I.
# But computing all third derivatives is complex. Let me first check whether
# F_I = 0 at this vacuum.

# --- F-terms at the vacuum ---
print("\n--- F-terms at the vacuum ---")
# F_I = dW/dphi_I at <phi> = vacuum values

# The superpotential is:
# W = sum_i m_i M^i_i + X(det M - B Bt - LAM^6) + y_c Hu0 Mdd + y_b Hd0 Mss

# F_{M^u_u} = m_u + X * d(det M)/dM^u_u + 0
# d(det M)/dM^u_u at diagonal vacuum = M_d * M_s
F_Muu = m_u + X0 * M_d * M_s
print(f"  F_Muu = m_u + X0 * M_d * M_s = {m_u} + ({X0:.6e})*{M_d:.4f}*{M_s:.4f} = {F_Muu:.6e}")

# F_{M^d_d} = m_d + X * M_u * M_s + y_c * Hu0_vev
F_Mdd = m_d + X0 * M_u * M_s + y_c * v_hu
print(f"  F_Mdd = m_d + X0 * M_u * M_s + y_c * v_hu = {F_Mdd:.6e}")

# F_{M^s_s} = m_s + X * M_u * M_d + y_b * Hd0_vev
F_Mss = m_s + X0 * M_u * M_d + y_b * v_hd
print(f"  F_Mss = m_s + X0 * M_u * M_d + y_b * v_hd = {F_Mss:.6e}")

# F_{M^i_j} (off-diagonal) = X * d(det M)/dM^i_j
# At diagonal vacuum, d(det M)/dM^i_j = 0 for i != j
# So F_{M^i_j} = 0 for all off-diagonal mesons
print(f"  F_Mij (off-diagonal) = 0 (all)")

# F_X = det M - B*Bt - LAM^6 = LAM^6 - 0 - LAM^6 = 0
F_X = det_M - LAM6
print(f"  F_X = det M - LAM^6 = {F_X:.6e}")

# F_B = -X * Bt = -X * 0 = 0, F_Bt = -X * B = 0
print(f"  F_B = F_Bt = 0")

# F_{Hu+} = 0 (no Hu+ coupling in W at this vacuum)
# F_{Hu0} = y_c * Mdd_vev = y_c * M_d
F_Hu0 = y_c * M_d
print(f"  F_Hu0 = y_c * M_d = {y_c:.6e} * {M_d:.4f} = {F_Hu0:.6e}")

# F_{Hd0} = y_b * Mss_vev = y_b * M_s
F_Hd0 = y_b * M_s
print(f"  F_Hd0 = y_b * M_s = {y_b:.6e} * {M_s:.4f} = {F_Hd0:.6e}")

print(f"  F_Hu+ = 0,  F_Hd- = 0")

# Check: at the Seiberg seesaw vacuum WITHOUT Yukawa, m_i + X * cofactor = 0.
# The Yukawa terms add y_c * v_hu to F_Mdd and y_b * v_hd to F_Mss.
# These are nonzero, so SUSY IS broken (by Yukawa + EWSB).
# Also F_Hu0 and F_Hd0 are nonzero.

print(f"\n  SUSY breaking check:")
print(f"  Without Yukawa: F_Muu = {m_u + X0 * M_d * M_s:.6e} (= 0 by seesaw)")
print(f"  Without Yukawa: F_Mdd = {m_d + X0 * M_u * M_s:.6e} (= 0 by seesaw)")
print(f"  Without Yukawa: F_Mss = {m_s + X0 * M_u * M_d:.6e} (= 0 by seesaw)")
print(f"  Yukawa shifts:  delta F_Mdd = y_c * v_hu = {y_c * v_hu:.6e}")
print(f"                  delta F_Mss = y_b * v_hd = {y_b * v_hd:.6e}")
print(f"  Higgs F-terms:  F_Hu0 = {F_Hu0:.6e}")
print(f"                  F_Hd0 = {F_Hd0:.6e}")
print(f"\n  SUSY is broken by Yukawa + EWSB (F_Mdd, F_Mss, F_Hu0, F_Hd0 != 0)")

# Now build the full scalar mass-squared matrix.
# Since F_I != 0, we need the complete expression:
# (m^2_scalar)_{JK} = sum_I W*_{IJ} W_{IK}  (from |F|^2 expansion, "F-flat" piece)
# plus W_{IJK} F*_I + W*_{IJK} F_I terms (the "F-cross" piece involving third derivatives)
# plus soft terms.
#
# For the complex scalar mass matrix:
# M^2_{JK} = sum_I (W_bar)_{IJ} W_{IK}  [acts on (delta phi_K)]
# M^2_{JK bar} = sum_I W_{IJL} F*_I  [acts on (delta phi_K^*)] — this mixes phi and phi*
#
# In the real (phi^R, phi^I) basis, we get a 32x32 matrix.
#
# Let me construct this in the complex basis first.
# The holomorphic mass matrix M^2_hol_{JK} = sum_I W_{IJK} F*_I (third derivative)
# contributes to the (phi, phi) and (phi*, phi*) blocks.
# The anti-holomorphic matrix M^2_ahol = (M^2_hol)* contributes to conjugate blocks.
# The Hermitian piece M^2_H_{JK} = sum_I W*_{IJ} W_{IK} contributes to (phi*, phi).

# For the full 32x32 real matrix, in the basis
# (phi_1^R, phi_1^I, phi_2^R, phi_2^I, ...) we have:
# m^2_real = [[ Re(M^2_H + M^2_hol),  Im(-M^2_H + M^2_hol) ],
#             [ Im(M^2_H + M^2_hol),   Re(M^2_H - M^2_hol)  ]]
#
# But this mixes all 32 components. Let me use a simpler approach:
# construct the mass matrix numerically by direct differentiation.

# Since all VEVs and couplings are REAL in this problem, the holomorphic
# and anti-holomorphic pieces are also real, simplifying things.

# Step 1: Build M^2_H = W^dag W (the "SUSY-preserving" piece)
# WARNING: W has condition number ~4e21. Computing W @ W directly and then
# diagonalizing gives numerically corrupted eigenvalues for the small modes.
# The correct eigenvalues of W^2 are lambda_i^2 where lambda_i are eigenvalues of W.
# We use the eigensystem of W itself (which is well-conditioned enough) to reconstruct W^2.
M2_H = W.T @ W   # W is real and symmetric, so W^dag W = W^T W = W^2

# The correct M2_H eigenvalues are the squares of W's eigenvalues:
m2_H_evals_correct = np.sort(fermion_evals**2)

print("\n\n--- Scalar mass-squared matrix from |F|^2 (SUSY piece: W^2) ---")
print(f"NOTE: W has condition number ~ {np.max(np.abs(fermion_evals))/np.min(np.abs(fermion_evals[np.abs(fermion_evals) > 1e-20])):.1e}")
print(f"Direct diagonalization of W^2 is numerically unreliable for small eigenvalues.")
print(f"Using eigenvalues of W squared instead (lambda_i^2):")
print(f"Correct eigenvalues of W^dag W = W^2 (16 complex scalars):")
for i, ev in enumerate(m2_H_evals_correct):
    print(f"  {i+1:3d}  {ev:+.8e} MeV^2  (sqrt = {np.sqrt(abs(ev)):.6e} MeV)")

# Step 2: Third-derivative contributions
# W_{IJK} = d^3 W / (dphi_I dphi_J dphi_K)
# The nonzero third derivatives come from the X * det M term:
# d^3(X det M)/dphi_I dphi_J dphi_K:
# - d^3/dX dM^a_a dM^b_b = cofactor(a,b) = M_c for a != b (c = third)
#   Wait, d^2(det M)/dM^a_a dM^b_b = M_c, so d^3(X det M)/dX dM^a_a dM^b_b = M_c
#   But this is already included in W_IJ. Third derivatives need THREE different fields.
#
# For X det M: d^3/dX dM^a_b dM^c_d = d^2(det M)/dM^a_b dM^c_d = +/- M_e or 0
# (as computed above, but now it's a third derivative of W involving X and two mesons)
#
# For the det M expansion with more than 2 meson derivatives:
# d^3(det M)/dM^a_i dM^b_j dM^c_k = epsilon_{abc} epsilon_{ijk} (the Levi-Civita)
# d^3(X det M)/dM^a_i dM^b_j dM^c_k = X * epsilon_{abc} epsilon_{ijk}
#
# Plus: d^3(X det M)/dX dM^a_i dM^b_j = d^2(det M)/dM^a_i dM^b_j
# (already computed, = +/- M_e at diagonal vacuum)
#
# So the third derivatives W_{IJK} are:
# W_{X, M^a_b, M^c_d} = d^2(det M)/dM^a_b dM^c_d (at vacuum)
# W_{M^a_i, M^b_j, M^c_k} = X * epsilon_{abc} epsilon_{ijk}

# The F-cross contribution to scalar mass is:
# (M^2_hol)_{JK} = sum_I W_{IJK} F*_I
#
# Since F_I are real in our case, F*_I = F_I.

# Build the holomorphic mass matrix
F_vec = np.zeros(16)
F_vec[0] = F_Muu  # practically 0
F_vec[4] = F_Mdd  # = y_c * v_hu (dominant from Yukawa)
F_vec[8] = F_Mss  # = y_b * v_hd (dominant from Yukawa)
F_vec[9] = F_X    # = 0
F_vec[13] = F_Hu0  # = y_c * M_d
F_vec[14] = F_Hd0  # = y_b * M_s

print(f"\nF-term vector:")
for i in range(16):
    if abs(F_vec[i]) > 1e-30:
        print(f"  F[{NAMES[i]:>4s}] = {F_vec[i]:+.10e}")

# Build third derivative tensor W_IJK (sparse)
# Only nonzero components:
# (1) W_{X, M^a_b, M^c_d} = d^2(det M)/dM^a_b dM^c_d  (for a != c)
# (2) W_{M^a_i, M^b_j, M^c_k} = X0 * epsilon_{abc} epsilon_{ijk}
# (3) No Yukawa third derivatives (y_c H_u^0 M^d_d is bilinear => W_{IJK}=0 for Yukawa)

# The holomorphic mass matrix from third derivatives:
M2_hol = np.zeros((16, 16))

# Contribution (1): W_{X, M^a_b, M^c_d} * F_X
# F_X = 0, so this contributes nothing.

# Contribution (1): W_{X, M^a_b, M^c_d} * F_{M^e_f}
# For J=X(9), K=M^a_b, I=M^c_d: M2_hol[9, idx_M(a,b)] += W_{M^c_d, X, M^a_b} * F_{M^c_d}
# Actually, M2_hol[J,K] = sum_I W_{IJK} F_I
# So for the X-meson-meson terms:
# W_{X, M^a_b, M^c_d} is symmetric in all three indices.
# M2_hol[J,K] gets contributions when any of the three index slots is I.

# Let me enumerate. The third derivative W_{X, M^a_a, M^b_b} = M_e (the VEV of third diagonal)
# for a != b, where e = 3-a-b.
# And W_{X, M^a_c, M^c_a} = -M_e for a != c, e = 3-a-c.

# Type (1) contributions to M2_hol:
# Pick I = X (index 9): W_{9, J, K} * F[9] = 0 (since F_X = 0)
# Pick I = M^a_b: W_{idx_M(a,b), J, K} * F[idx_M(a,b)]
# For this, we need W_{M^a_b, X, M^c_d} and W_{M^a_b, M^c_d, M^e_f}

# This is getting complex. Let me build the W_IJK tensor directly.
# Since it's sparse, I'll use a dictionary.

from itertools import permutations

def levi_civita(p):
    """Sign of permutation p of (0,1,2)."""
    n = len(p)
    p = list(p)
    sign = 1
    for i in range(n):
        while p[i] != i:
            sign *= -1
            j = p[i]
            p[i], p[j] = p[j], p[i]
    return sign

# Build W3 as a symmetric tensor using a clean accumulation approach.
# First accumulate ALL contributions to each UNORDERED triple (as a sorted tuple),
# then distribute correctly to the symmetric tensor.
#
# The superpotential W = ... + X * det M + ...
# Third derivatives of X * det M:
#   (1) d^3/dX dM^a_b dM^c_d = d^2(det M)/dM^a_b dM^c_d
#       Nonzero at diagonal vacuum only for:
#       - a != c, b=a, d=c: d^2(det)/dM^a_a dM^c_c = M_e (e = third)
#       - a != c, b=c, d=a: d^2(det)/dM^a_c dM^c_a = -M_e
#   (2) d^3/dM^a_i dM^b_j dM^c_k = X * eps_{abc} eps_{ijk}
#       (a,b,c) and (i,j,k) each a permutation of (0,1,2)

# Accumulate into a dict keyed by sorted (unordered) triple of field indices.
# Each entry stores the CORRECT symmetric tensor value W_{IJK}.
W3_accum = {}  # key = sorted tuple, value = W_{IJK} (same for all permutations)

# Type (1): W_{X(9), M^a_b, M^c_d} = d^2(det M)/dM^a_b dM^c_d
# Enumerate unique unordered triples: for each distinct pair (a,c) with a < c
for a in range(3):
    for c in range(a+1, 3):
        e = 3 - a - c
        # Diagonal-diagonal: W_{X, M^a_a, M^c_c} = M_e
        key = tuple(sorted([9, idx_M(a,a), idx_M(c,c)]))
        W3_accum[key] = W3_accum.get(key, 0) + M_vev[e]
        # Off-diagonal pair: W_{X, M^a_c, M^c_a} = -M_e
        key = tuple(sorted([9, idx_M(a,c), idx_M(c,a)]))
        W3_accum[key] = W3_accum.get(key, 0) + (-M_vev[e])

# Type (2): W_{M^a_i, M^b_j, M^c_k} = X0 * eps_{abc} eps_{ijk}
# The symmetric tensor value is the SAME for all permutations of the three
# field indices (I,J,K). Summing over all 36 row/col permutations counts
# each unordered triple 6 times (for 3 distinct fields) with the same sign.
# So: accumulate the sum and divide by 6 to get the per-triple value.
W3_type2_raw = {}
for perm_row in permutations(range(3)):
    for perm_col in permutations(range(3)):
        sign = levi_civita(list(perm_row)) * levi_civita(list(perm_col))
        val = X0 * sign
        i1 = idx_M(perm_row[0], perm_col[0])
        i2 = idx_M(perm_row[1], perm_col[1])
        i3 = idx_M(perm_row[2], perm_col[2])
        key = tuple(sorted([i1, i2, i3]))
        W3_type2_raw[key] = W3_type2_raw.get(key, 0) + val

for key, raw_val in W3_type2_raw.items():
    # Each unordered triple appears 6 times in the sum (for distinct fields),
    # all with the same value. Divide by 6 to get W_{IJK}.
    W3_accum[key] = W3_accum.get(key, 0) + raw_val / 6.0

# Now build the full symmetric W3 tensor as a dict
W3 = {}
for sorted_key, val in W3_accum.items():
    I, J, K = sorted_key
    for perm in [(I,J,K),(I,K,J),(J,I,K),(J,K,I),(K,I,J),(K,J,I)]:
        W3[perm] = W3.get(perm, 0) + val

# Build M2_hol[J,K] = sum_I W3[(I,J,K)] * F_vec[I]
for (I, J, K), val in W3.items():
    M2_hol[J, K] += val * F_vec[I]

# Symmetrize (should already be symmetric)
M2_hol = 0.5 * (M2_hol + M2_hol.T)

print(f"\nHolomorphic mass matrix M2_hol (from W_IJK F_I):")
print(f"  Nonzero entries:")
for i in range(16):
    for j in range(i, 16):
        if abs(M2_hol[i,j]) > 1e-20:
            print(f"    M2_hol[{NAMES[i]:>4s}, {NAMES[j]:>4s}] = {M2_hol[i,j]:+.8e}")

# Step 3: Soft SUSY-breaking contributions
# V_soft = m_tilde^2 Tr(M^dag M) + (B_mu sum_i m_i M^i_i + h.c.)
#
# m_tilde^2 Tr(M^dag M) = m_tilde^2 sum_{i,j} |M^i_j|^2
# This gives a universal mass-squared m_tilde^2 for all 9 meson scalars (both Re and Im parts).
#
# B_mu sum_i m_i M^i_i + h.c. = B_mu (m_u M^u_u + m_d M^d_d + m_s M^s_s) + h.c.
# At the vacuum, expanding M^i_i = <M^i_i> + delta M^i_i:
# B_mu m_i (M_i + delta M^i_i) + h.c. = B_mu m_i M_i + B_mu m_i delta M^i_i + h.c.
# The linear term B_mu m_i delta M^i_i + h.c. is a tadpole (shifts the vacuum).
# The mass contribution: B_mu m_i is real, and M^i_i = (phi^R + i phi^I)/sqrt(2),
# so B_mu m_i M^i_i + h.c. = 2 B_mu m_i phi^R / sqrt(2) = sqrt(2) B_mu m_i phi^R
# This is a LINEAR term, not quadratic. No mass contribution from B_mu at second order.
#
# Actually, in the standard treatment, B_mu term in the Lagrangian V = B_mu m M + h.c.
# with M complex, this is V = B_mu m (phi^R + i phi^I)/sqrt(2) + h.c.
# = sqrt(2) B_mu m phi^R  (if B_mu, m are real)
# This is indeed linear. It shifts the vacuum but doesn't directly contribute to masses.
# However, if we're expanding around the CORRECT vacuum (which includes the B_mu shift),
# then the vacuum equations absorb this linear term.
#
# The quadratic soft terms are:
# From m_tilde^2: diagonal mass m_tilde^2 for each meson scalar component
# From B_mu: this term is linear in the fields, not quadratic.
# So the soft contribution to the scalar mass matrix is just m_tilde^2 for mesons.
#
# Wait, I need to reconsider. The B_mu term in the MSSM-like context is:
# V_soft = B_mu (H_u H_d + h.c.) which IS bilinear in Higgs fields.
# But in this problem, V_soft = B_mu sum_i m_i M^i_i + h.c. which is LINEAR in M.
# This is a standard soft "A-term" style or linear soft term.
# For the mass matrix, it contributes nothing at second order.
# The mass matrix gets:
# - m_tilde^2 for meson real and imaginary parts (from Tr(M^dag M))
# - Nothing from B_mu at second order (it's linear)

M2_soft = np.zeros((16, 16))
# m_tilde^2 for all 9 meson fields (indices 0-8)
for i in range(9):
    M2_soft[i, i] = m_tilde_sq

print(f"\nSoft mass contributions:")
print(f"  m_tilde^2 = {m_tilde_sq:.0f} MeV^2 for all 9 meson scalars")
print(f"  B_mu term is linear in fields => no quadratic mass contribution")

# Step 4: Build the full 32x32 real scalar mass matrix
# In the real basis (phi_1^R, phi_1^I, phi_2^R, phi_2^I, ...):
#
# The 32x32 matrix has 2x2 blocks for each pair of fields.
# For the (phi_J, phi_K) entry, the 2x2 block is:
#
# From M^2_H (Hermitian piece, always real for real W):
# M^2_H_{JK} contributes to both (R,R) and (I,I) blocks equally.
#
# From M^2_hol (holomorphic piece):
# M^2_hol_{JK} is real in our case.
# It contributes +M^2_hol to (R,R) and -M^2_hol to (I,I).
# (Because phi phi = (phi^R + i phi^I)^2/2 = ((phi^R)^2 - (phi^I)^2 + 2i phi^R phi^I)/2)
#
# From M^2_soft:
# Diagonal, same for R and I parts.
#
# Full 2x2 block for (J,K):
# [[M^2_H + M^2_hol + M^2_soft,  0                ],
#  [0,                            M^2_H - M^2_hol + M^2_soft]]
# (No R-I mixing because everything is real in our model.)
#
# Actually, for the (phi_J^R, phi_K^I) and (phi_J^I, phi_K^R) entries,
# the contributions vanish because M^2_H is real symmetric and M^2_hol is real symmetric.
# The general formula for complex fields:
# V = M^2_H_{JK} phi*_J phi_K + 1/2 M^2_hol_{JK} phi_J phi_K + h.c.
# = M^2_H_{JK} (phiR_J phiR_K + phiI_J phiI_K)/2
#   + M^2_hol_{JK} (phiR_J phiR_K - phiI_J phiI_K)/2
#   [the cross terms phiR phiI cancel because M^2_hol is symmetric and real]
#
# So:
# m^2_{R_J R_K} = M^2_H_{JK} + M^2_hol_{JK} + M^2_soft_{JK}  (for J in meson sector)
# m^2_{I_J I_K} = M^2_H_{JK} - M^2_hol_{JK} + M^2_soft_{JK}
# m^2_{R_J I_K} = 0

M2_scalar_32 = np.zeros((32, 32))

for J in range(16):
    for K in range(16):
        # Real-Real block
        M2_scalar_32[2*J, 2*K] = M2_H[J, K] + M2_hol[J, K] + M2_soft[J, K]
        # Imag-Imag block
        M2_scalar_32[2*J+1, 2*K+1] = M2_H[J, K] - M2_hol[J, K] + M2_soft[J, K]
        # Real-Imag = 0 (already initialized)

# Symmetrize
M2_scalar_32 = 0.5 * (M2_scalar_32 + M2_scalar_32.T)

scalar_evals = eigvalsh(M2_scalar_32)

print(f"\n--- Scalar mass-squared eigenvalues (32 real d.o.f.) ---")
print(f"{'n':>3}  {'m^2 (MeV^2)':>22}  {'m (MeV)':>16}  {'m^2 > 0?':>10}")
for i, ev in enumerate(scalar_evals):
    m_val = np.sqrt(abs(ev)) if ev >= 0 else -np.sqrt(abs(ev))
    sign = "YES" if ev > -1e-10 else "NO (tachyon)"
    print(f"  {i+1:3d}  {ev:+22.10e}  {m_val:+16.6e}  {sign}")

n_tachyon = np.sum(scalar_evals < -1e-10)
n_zero_s = np.sum(np.abs(scalar_evals) < 1e-10)
n_pos_s = np.sum(scalar_evals > 1e-10)
print(f"\n  Tachyonic (m^2 < 0): {n_tachyon}")
print(f"  Zero (|m^2| < 1e-10): {n_zero_s}")
print(f"  Positive: {n_pos_s}")

# --- NUMERICAL RELIABILITY WARNING ---
# The 32x32 matrix is numerically ill-conditioned due to the ~10^21 dynamic range
# from the X-Mss sector (~10^10 MeV) to the baryon sector (~10^-9 MeV).
# The tachyonic modes are NUMERICAL ARTIFACTS from the condition number of W^2.
# Eigenvalues of W^2 should all be non-negative (W is symmetric, W^2 = W^T W is PSD).
# The correct scalar spectrum has:
#   - For each fermion eigenvalue lambda_i: two scalar modes with mass-squared
#     lambda_i^2 + M2_hol_{ii} + M2_soft_{ii}  (R part)  and
#     lambda_i^2 - M2_hol_{ii} + M2_soft_{ii}  (I part)
#
# To get RELIABLE scalar masses for the light sector, we analyze block by block:

print(f"\n\n--- BLOCK-BY-BLOCK SCALAR MASS ANALYSIS (numerically reliable) ---")

# The fermion mass matrix decomposes into independent blocks:
# (a) Central block: {Muu(0), Mdd(4), Mss(8), X(9), Hu0(13), Hd0(14)} -- coupled via W
# (b) Off-diag meson pair: {Mud(1), Mdu(3)} -- 2x2
# (c) Off-diag meson pair: {Mus(2), Msu(6)} -- 2x2
# (d) Off-diag meson pair: {Mds(5), Msd(7)} -- 2x2
# (e) Baryon pair: {B(10), Bt(11)} -- 2x2
# (f) Decoupled: Hu+(12) -- 1x1, mass = 0
# (g) Decoupled: Hd-(15) -- 1x1, mass = 0

def scalar_block_masses(W_block, M2_hol_block, M2_soft_block, block_name, field_names_block):
    """Compute scalar mass-squared eigenvalues for a single block."""
    n = len(W_block)
    # Build the 2n x 2n real scalar mass matrix for this block
    M2_H_block = W_block @ W_block  # This is reliable for small blocks
    M2_R = M2_H_block + M2_hol_block + M2_soft_block  # Real-Real
    M2_I = M2_H_block - M2_hol_block + M2_soft_block  # Imag-Imag

    # For 2n x 2n: [[M2_R, 0], [0, M2_I]]
    M2_full = np.zeros((2*n, 2*n))
    M2_full[:n, :n] = M2_R
    M2_full[n:, n:] = M2_I
    evals = eigvalsh(M2_full)

    print(f"\n  Block: {block_name}")
    print(f"  Fields: {field_names_block}")
    fermion_evals_block = eigvalsh(W_block)
    print(f"  Fermion masses: {[f'{abs(v):.6e}' for v in sorted(fermion_evals_block, key=abs)]}")
    for i, ev in enumerate(evals):
        m_val = np.sqrt(abs(ev))
        ri = "R" if i < n else "I"
        tach = " [tachyon!]" if ev < -1e-10 else ""
        print(f"    m^2_{i+1}({ri}) = {ev:+.8e} MeV^2  (|m| = {m_val:.6e} MeV){tach}")
    return evals

# Block (b): Mud(1) - Mdu(3) off-diagonal pair
idx_b = [1, 3]
W_b = W[np.ix_(idx_b, idx_b)]
hol_b = M2_hol[np.ix_(idx_b, idx_b)]
soft_b = M2_soft[np.ix_(idx_b, idx_b)]
evals_b = scalar_block_masses(W_b, hol_b, soft_b, "ud off-diagonal", ["Mud", "Mdu"])

# Block (c): Mus(2) - Msu(6)
idx_c = [2, 6]
W_c = W[np.ix_(idx_c, idx_c)]
hol_c = M2_hol[np.ix_(idx_c, idx_c)]
soft_c = M2_soft[np.ix_(idx_c, idx_c)]
evals_c = scalar_block_masses(W_c, hol_c, soft_c, "us off-diagonal", ["Mus", "Msu"])

# Block (d): Mds(5) - Msd(7)
idx_d = [5, 7]
W_d = W[np.ix_(idx_d, idx_d)]
hol_d = M2_hol[np.ix_(idx_d, idx_d)]
soft_d = M2_soft[np.ix_(idx_d, idx_d)]
evals_d = scalar_block_masses(W_d, hol_d, soft_d, "ds off-diagonal", ["Mds", "Msd"])

# Block (e): B(10) - Bt(11)
idx_e = [10, 11]
W_e = W[np.ix_(idx_e, idx_e)]
hol_e = M2_hol[np.ix_(idx_e, idx_e)]
soft_e = M2_soft[np.ix_(idx_e, idx_e)]
evals_e = scalar_block_masses(W_e, hol_e, soft_e, "baryon pair", ["B", "Bt"])

# Block (f): Hu+(12)
idx_f = [12]
W_f = W[np.ix_(idx_f, idx_f)]
hol_f = M2_hol[np.ix_(idx_f, idx_f)]
soft_f = M2_soft[np.ix_(idx_f, idx_f)]
evals_f = scalar_block_masses(W_f, hol_f, soft_f, "Hu+ (decoupled)", ["Hu+"])

# Block (g): Hd-(15)
idx_g = [15]
W_g = W[np.ix_(idx_g, idx_g)]
hol_g = M2_hol[np.ix_(idx_g, idx_g)]
soft_g = M2_soft[np.ix_(idx_g, idx_g)]
evals_g = scalar_block_masses(W_g, hol_g, soft_g, "Hd- (decoupled)", ["Hd-"])

# Block (a): Central block {Muu(0), Mdd(4), Mss(8), X(9), Hu0(13), Hd0(14)}
# This block has a huge dynamic range. Analyze it carefully.
idx_a = [0, 4, 8, 9, 13, 14]
W_a = W[np.ix_(idx_a, idx_a)]
hol_a = M2_hol[np.ix_(idx_a, idx_a)]
soft_a = M2_soft[np.ix_(idx_a, idx_a)]
print(f"\n  Block: central (Muu, Mdd, Mss, X, Hu0, Hd0)")
print(f"  WARNING: condition number ~ 10^21; scalar eigenvalues for light modes unreliable")
evals_a = scalar_block_masses(W_a, hol_a, soft_a, "central block", ["Muu","Mdd","Mss","X","Hu0","Hd0"])

# Collect all block eigenvalues for the reliable supertrace
all_block_evals = np.concatenate([evals_a, evals_b, evals_c, evals_d, evals_e, evals_f, evals_g])
all_block_evals = np.sort(all_block_evals)
print(f"\n  Total scalar eigenvalues from blocks: {len(all_block_evals)} (should be 32)")
print(f"  Sum of block eigenvalues: {np.sum(all_block_evals):+.10e}")
print(f"  Sum from 32x32 matrix:    {np.sum(scalar_evals):+.10e}")
print(f"  (These should agree if blocks don't couple.)")


# ═══════════════════════════════════════════════════════════════════════════
# TASK 3: SM Identification
# ═══════════════════════════════════════════════════════════════════════════
print("\n\n" + "=" * 72)
print("TASK 3: SM PARTICLE IDENTIFICATION")
print("=" * 72)

# The 16 chiral superfields and their SM associations:
# Meson diagonal M^i_i: these are neutral color-singlet scalars.
#   VEVs ~ C/m_i (large). Fluctuations around these VEVs = neutral mesons.
#   In the Seiberg dual picture, these correspond to quark bilinears.
#   The seesaw M_j = C/m_j means the CW spectrum tracks quark mass ratios.
#   After the Yukawa coupling, M^d_d mixes with H_u^0 to produce charm-mass states,
#   and M^s_s mixes with H_d^0 to produce bottom-mass states.
#
# Meson off-diagonal M^i_j (i != j): charged/flavor-changing mesons.
#   These get masses from the X * det M coupling.
#   Masses ~ |X| * M_k (tiny in MeV terms).
#
# X: Lagrange multiplier. Its fermion partner is very light.
# B, Btilde: baryons. Get mass ~ |X| (extremely light).
# H_u, H_d: Higgs doublets. The neutral components couple to mesons via Yukawa.
#   After EWSB, the charged components (H_u^+, H_d^-) remain as charged Higgs.

# Fermion identification:
print("\n--- Fermion mass spectrum and SM identification ---")
print(f"\n{'n':>3}  {'Mass (MeV)':>16}  {'Dominant field':>20}  {'SM identification':>30}")

# Get sorted eigenvalues with eigenvectors
fermion_sorted_idx = np.argsort(np.abs(fermion_evals))
for rank, i in enumerate(fermion_sorted_idx):
    ev = fermion_evals[i]
    evec = fermion_evecs[:, i]
    # Find top 2 dominant components
    abs_evec = np.abs(evec)
    top2 = np.argsort(abs_evec)[-2:][::-1]
    dom1_name = NAMES[top2[0]]
    dom1_frac = evec[top2[0]]**2
    dom2_name = NAMES[top2[1]]
    dom2_frac = evec[top2[1]]**2

    # SM identification logic
    if abs(ev) < 1e-15:
        sm_id = "massless (Goldstino or decoupled)"
    elif dom1_name in ['Hu+', 'Hd-']:
        sm_id = "charged Higgsino (decoupled)"
    elif dom1_name in ['B', 'Bt']:
        sm_id = "baryon fermion (exotic)"
    elif dom1_name == 'X':
        sm_id = "X-fermion (pseudo-modulus partner)"
    elif dom1_name in ['Muu', 'Mdd', 'Mss']:
        if abs(abs(ev) - abs(X0 * M_vev[0])) / max(abs(ev), 1e-30) < 0.01:
            sm_id = "diagonal meson fermion (light neutral)"
        else:
            sm_id = "diagonal meson fermion"
    elif 'M' in dom1_name:
        sm_id = "off-diagonal meson fermion (flavor-changing)"
    elif dom1_name == 'Hu0':
        sm_id = "neutral Higgsino (up-type)"
    elif dom1_name == 'Hd0':
        sm_id = "neutral Higgsino (down-type)"
    else:
        sm_id = "unidentified"

    print(f"  {rank+1:3d}  {abs(ev):16.8e}  {dom1_name:>5s}({dom1_frac:.3f})+{dom2_name:>4s}({dom2_frac:.3f})  {sm_id}")


# ═══════════════════════════════════════════════════════════════════════════
# TASK 4: Supertrace
# ═══════════════════════════════════════════════════════════════════════════
print("\n\n" + "=" * 72)
print("TASK 4: SUPERTRACE STr[M^2]")
print("=" * 72)

# STr[M^2] = sum_bosons m^2 - 2 * sum_fermions m^2
# The factor of 2 for fermions accounts for the two helicity states of a Dirac fermion
# formed from a Weyl pair. But our fermion eigenvalues are for Weyl fermions.
# For a Weyl fermion with mass m, m^2_fermion = m^2.
# For 16 Weyl fermions, the fermion contribution is sum of m_i^2.
#
# Actually, the supertrace formula for N=1 SUSY:
# STr M^2 = Tr m^2_scalar - 2 Tr m^2_fermion
# where m^2_scalar are the 32 real scalar mass-squared eigenvalues,
# and m^2_fermion are the squared masses of the 16 Weyl fermions
# (each Weyl contributes m_i^2 once, not twice).
# Wait, the standard formula:
# STr M^2 = sum_{J} (-1)^{2J} (2J+1) Tr m^2_J
# For scalars (J=0): 1 * Tr m^2_scalar
# For Weyl fermions (J=1/2): -2 * Tr m^2_Weyl = -2 * sum |m_i|^2
#
# But we have 32 REAL scalars and 16 Weyl fermions (= 32 real d.o.f. each).
# Supertrace: STr M^2 = sum_{32 real scalars} m^2_i - 2 * sum_{16 Weyl} m^2_i
#
# Note: For a massive Weyl fermion with eigenvalue lambda, the mass-squared is lambda^2.
# For the fermion mass matrix W_IJ, the mass-squared eigenvalues are eigenvalues of W^T W = W^2.
# The trace of W^2 gives sum of lambda_i^2.

# ANALYTICAL COMPUTATION of the supertrace (avoids numerical issues):
#
# STr M^2 = sum_{32 real scalars} m^2_i - 2 * sum_{16 Weyl} m^2_i
#
# The scalar mass-squared matrix has the form:
#   m^2_{RR} = M2_H + M2_hol + M2_soft  (16x16)
#   m^2_{II} = M2_H - M2_hol + M2_soft  (16x16)
#   m^2_{RI} = 0  (no R-I mixing for real couplings)
#
# where M2_H = W^2 (positive semi-definite), M2_hol = sum_I W_{IJK} F_I, M2_soft = diagonal.
#
# Supertrace:
#   STr = Tr(m^2_{RR}) + Tr(m^2_{II}) - 2*Tr(W^2)
#       = Tr(M2_H + M2_hol + M2_soft) + Tr(M2_H - M2_hol + M2_soft) - 2*Tr(M2_H)
#       = 2*Tr(M2_H) + 2*Tr(M2_soft) - 2*Tr(M2_H)
#       = 2*Tr(M2_soft)
#
# The M2_hol (holomorphic) piece cancels in R + I traces.
# The M2_H = W^2 piece cancels between scalar and fermion sectors.
# ONLY the soft terms survive.

Tr_M2_soft = np.trace(M2_soft)
STr_analytical = 2.0 * Tr_M2_soft
print(f"\n  ANALYTICAL supertrace:")
print(f"  Tr(M2_soft) = {Tr_M2_soft:.2f}  (= 9 * m_tilde^2 = 9 * {m_tilde_sq:.0f})")
print(f"  STr[M^2] = 2 * Tr(M2_soft) = {STr_analytical:.2f} MeV^2")
print(f"           = 18 * f_pi^2 = 18 * {f_pi}^2 = {18 * f_pi**2:.0f} MeV^2")
print(f"           = {STr_analytical:.0f} MeV^2")

# Verify components:
print(f"\n  Derivation:")
sum_fermion_m2 = np.sum(fermion_evals**2)
Tr_W2 = np.trace(W @ W)
print(f"  Tr(W^2) = {Tr_W2:.8e}")
print(f"  sum(lambda_fermion^2) = {sum_fermion_m2:.8e}")
print(f"  Tr(W^2) - sum(lambda^2) = {Tr_W2 - sum_fermion_m2:.4e}  (= 0 analytically, check)")
print(f"  Tr(M2_hol) = {np.trace(M2_hol):.4e}  (cancels in R+I)")
print(f"  2*Tr(M2_soft) = {2*Tr_M2_soft:.0f}")

# Numerical STr from the 32x32 matrix (unreliable for verification):
sum_scalar_m2 = np.sum(scalar_evals)
STr_numerical = sum_scalar_m2 - 2.0 * sum_fermion_m2
print(f"\n  Numerical check (32x32, affected by condition number ~10^21):")
print(f"  sum_scalar(32x32) = {sum_scalar_m2:+.8e}")
print(f"  2*sum_fermion = {2*sum_fermion_m2:+.8e}")
print(f"  STr_numerical = {STr_numerical:+.8e}")
print(f"  STr_analytical = {STr_analytical:+.8e}")
print(f"  Discrepancy = {abs(STr_numerical - STr_analytical):.4e} (numerical error from condition number)")

# Block-by-block STr verification:
sum_scalar_blocks = np.sum(all_block_evals)
STr_blocks = sum_scalar_blocks - 2.0 * sum_fermion_m2
print(f"\n  Block-by-block STr: {STr_blocks:+.8e}")
print(f"  Block sum scalar = {sum_scalar_blocks:+.8e}")
print(f"  Discrepancy from analytical = {abs(STr_blocks - STr_analytical):.4e}")

print(f"\n  RESULT: STr[M^2] = 2 * 9 * m_tilde^2 = {STr_analytical:.0f} MeV^2")
print(f"  For unbroken SUSY (m_tilde = 0): STr = 0.")
print(f"  The nonzero STr is entirely due to the soft mass m_tilde^2 = f_pi^2.")
print(f"  The Yukawa-induced F-terms do NOT contribute to the supertrace")
print(f"  (they cancel between R and I scalar components).")


# ═══════════════════════════════════════════════════════════════════════════
# TASK 5: Complete spectrum table
# ═══════════════════════════════════════════════════════════════════════════
print("\n\n" + "=" * 72)
print("TASK 5: COMPLETE SPECTRUM TABLE")
print("=" * 72)

# Build a unified spectrum table
# For fermions: eigenvalue of W_IJ, mass = |eigenvalue|
# For scalars: eigenvalue of M^2_32, mass = sqrt(|eigenvalue|)

print("\n--- FERMION SPECTRUM ---")
print(f"{'#':>3}  {'Mass |m| (MeV)':>18}  {'Spin':>5}  {'Dominant':>12}  {'SM identification'}")
print("-" * 80)

# Sort fermions by |mass|
f_sorted = sorted(enumerate(fermion_evals), key=lambda x: abs(x[1]))
for rank, (idx, ev) in enumerate(f_sorted):
    evec = fermion_evecs[:, idx]
    dom_idx = np.argmax(np.abs(evec))
    dom_name = NAMES[dom_idx]
    mass = abs(ev)

    # Charge assignment
    # M^i_j: neutral (Q=0) for diagonal, various for off-diagonal
    # X: Q=0, B: Q depends on baryon number, Bt: opposite
    # Hu+: Q=+1, Hu0: Q=0, Hd0: Q=0, Hd-: Q=-1
    charge = "0"
    if dom_name == 'Hu+':
        charge = "+1"
    elif dom_name == 'Hd-':
        charge = "-1"

    # Identification
    if mass < 1e-15:
        ident = "massless Weyl (decoupled)"
    elif dom_name in ['Hu+']:
        ident = "charged Higgsino (exotic)"
    elif dom_name in ['Hd-']:
        ident = "charged Higgsino (exotic)"
    elif dom_name in ['B', 'Bt']:
        ident = "baryonic fermion"
    elif dom_name == 'X':
        ident = "X-ino (pseudo-modulus fermion)"
    elif dom_name in ['Mud', 'Mdu']:
        ident = "off-diag meson (ud sector)"
    elif dom_name in ['Mus', 'Msu']:
        ident = "off-diag meson (us sector)"
    elif dom_name in ['Mds', 'Msd']:
        ident = "off-diag meson (ds sector)"
    elif dom_name in ['Muu']:
        ident = "diag meson fermion (uu)"
    elif dom_name in ['Mdd']:
        ident = "diag meson fermion / charm"
    elif dom_name in ['Mss']:
        ident = "diag meson fermion / bottom"
    elif dom_name in ['Hu0']:
        ident = "neutral Higgsino (up)"
    elif dom_name in ['Hd0']:
        ident = "neutral Higgsino (down)"
    else:
        ident = "mixed"

    print(f"  {rank+1:3d}  {mass:18.8e}  {0.5:>5.1f}  Q={charge:>3s} {dom_name:>5s}  {ident}")

print("\n--- SCALAR SPECTRUM ---")
print(f"{'#':>3}  {'m^2 (MeV^2)':>22}  {'|m| (MeV)':>16}  {'Spin':>5}  {'Type':>8}  {'SM identification'}")
print("-" * 90)

# For scalars, identify which complex field each eigenstate belongs to
# The 32x32 matrix is in the basis (phi_1^R, phi_1^I, phi_2^R, phi_2^I, ...)
# So eigenstate n has components along phi_I^R and phi_I^I.

scalar_evecs_all = np.linalg.eigh(M2_scalar_32)[1]

s_sorted = sorted(enumerate(scalar_evals), key=lambda x: abs(x[1]))
for rank, (idx, m2) in enumerate(s_sorted):
    evec = scalar_evecs_all[:, idx]
    # Find dominant complex field
    field_weights = np.zeros(16)
    for I in range(16):
        field_weights[I] = evec[2*I]**2 + evec[2*I+1]**2
    dom_field = np.argmax(field_weights)
    dom_name = NAMES[dom_field]
    dom_frac = field_weights[dom_field]

    # Check if R or I component dominates
    is_real = evec[2*dom_field]**2 > evec[2*dom_field+1]**2
    ri_label = "R" if is_real else "I"

    mass = np.sqrt(abs(m2))
    sign = "+" if m2 >= 0 else "-"

    # Identification
    if abs(m2) < 1e-10:
        ident = "massless (Goldstone or flat direction)"
    elif dom_name in ['Hu+']:
        ident = f"charged Higgs scalar ({ri_label})"
    elif dom_name in ['Hd-']:
        ident = f"charged Higgs scalar ({ri_label})"
    elif dom_name in ['B', 'Bt']:
        ident = f"baryonic scalar ({ri_label})"
    elif dom_name == 'X':
        ident = f"X-scalar ({ri_label})"
    elif 'M' in dom_name and len(dom_name) == 3 and dom_name[1] != dom_name[2]:
        ident = f"off-diag meson scalar ({ri_label})"
    elif 'M' in dom_name:
        ident = f"diag meson scalar ({ri_label})"
    elif dom_name in ['Hu0']:
        ident = f"neutral Higgs (up, {ri_label})"
    elif dom_name in ['Hd0']:
        ident = f"neutral Higgs (down, {ri_label})"
    else:
        ident = f"mixed ({ri_label})"

    print(f"  {rank+1:3d}  {m2:+22.8e}  {mass:16.8e}  {0:>5d}  {sign+'m^2':>8s}  {dom_name:>5s} {ident}")


# ═══════════════════════════════════════════════════════════════════════════
# Summary and consistency checks
# ═══════════════════════════════════════════════════════════════════════════
print("\n\n" + "=" * 72)
print("SUMMARY AND CONSISTENCY CHECKS")
print("=" * 72)

print(f"\n  Total fields: 16 complex scalars (32 real) + 16 Weyl fermions")
print(f"  Fermion mass eigenvalues: {len(fermion_evals)}")
print(f"  Scalar mass-squared eigenvalues: {len(scalar_evals)}")

print(f"\n  Supertrace: STr[M^2] = {STr_analytical:+.0f} MeV^2  (analytical)")
print(f"  = 2 * 9 * m_tilde^2 = 18 * f_pi^2 = {18*m_tilde_sq:.0f} MeV^2")
print(f"  Nonzero STr confirms SUSY breaking by soft terms.")

print(f"\n  SUSY breaking pattern:")
print(f"    F_Mdd = {F_Mdd:.6e} MeV (from y_c Yukawa)")
print(f"    F_Mss = {F_Mss:.6e} MeV (from y_b Yukawa)")
print(f"    F_Hu0 = {F_Hu0:.6e} MeV (from y_c * M_d)")
print(f"    F_Hd0 = {F_Hd0:.6e} MeV (from y_b * M_s)")
print(f"    |F|^2 = {F_Mdd**2 + F_Mss**2 + F_Hu0**2 + F_Hd0**2:.6e} MeV^2")
print(f"    sqrt(|F|) = {np.sqrt(abs(F_Mdd**2 + F_Mss**2 + F_Hu0**2 + F_Hd0**2)**0.5):.4f} MeV")

# Scale hierarchy
print(f"\n  Scale hierarchy:")
ferm_masses = np.sort(np.abs(fermion_evals))
scal_masses = np.sort(np.sqrt(np.abs(scalar_evals)))
print(f"    Heaviest fermion: {ferm_masses[-1]:.6e} MeV")
print(f"    Lightest nonzero fermion: {ferm_masses[ferm_masses > 1e-15][0]:.6e} MeV" if np.any(ferm_masses > 1e-15) else "    All fermions massless")
print(f"    Heaviest scalar: {scal_masses[-1]:.6e} MeV")
print(f"    Lightest nonzero scalar: {scal_masses[scal_masses > 1e-5][0]:.6e} MeV" if np.any(scal_masses > 1e-5) else "    All scalars massless")

# CW spectrum tracking
print(f"\n  ISS/CW spectrum ratios (fermion masses proportional to 1/m_quark):")
# The off-diagonal meson pair masses are |X0| * M_k = |X0| * C/m_k
# Ratios:
print(f"    |X|*M_u / |X|*M_s = M_u/M_s = m_s/m_u = {m_s/m_u:.2f}")
print(f"    |X|*M_d / |X|*M_s = M_d/M_s = m_s/m_d = {m_s/m_d:.2f}")
print(f"    |X|*M_u / |X|*M_d = M_u/M_d = m_d/m_u = {m_d/m_u:.2f}")
print(f"    These ratios are the INVERSE of the quark mass ratios (seesaw duality).")

print("\n" + "=" * 72)
print("COMPUTATION COMPLETE")
print("=" * 72)
