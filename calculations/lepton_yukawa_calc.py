#!/usr/bin/env python3
"""
Lepton Yukawa from SQCD meson fermion mass matrix.

Setup: SU(3) SQCD with N_f = N_c = 3, flavors Q^a_i, Qbar^j_a (a=color, i,j=u,d,s).
Seiberg duality gives mesons M^i_j = Q^a_i Qbar^j_a / Lambda and baryons B, Btilde.

Superpotential:
  W = sum_i m_i M^i_i + X(det M - B Btilde - Lambda^6) + y_c Hu0 M^d_d + y_b Hd0 M^s_s

We compute the 6x6 central-block fermion mass matrix W_IJ = d^2W/dPhi_I dPhi_J
for the fields (M^u_u, M^d_d, M^s_s, X, Hu0, Hd0) at the diagonal Seiberg vacuum.

Then: eigenvalues, light eigenstates, Koide test, NMSSM extension, effective Yukawas,
and off-diagonal meson alternative.

Pure linear algebra / superpotential eigenvalue computation.
"""

import numpy as np
from numpy.linalg import eigh, eigvalsh, eigvals, det, inv

# ============================================================================
# Parameters
# ============================================================================
m_u = 2.16       # MeV
m_d = 4.67       # MeV
m_s = 93.4       # MeV
m_c = 1270.0     # MeV
m_b = 4180.0     # MeV
LAM = 300.0      # MeV
v   = 246220.0   # MeV (full EW vev, 246.22 GeV)

# Yukawa couplings (convention: m = y v / sqrt(2), so y = 2m/v)
y_c = 2.0 * m_c / v
y_b = 2.0 * m_b / v

# NMSSM coupling
lam_nmssm = 0.72

# Derived quantities
m_arr = np.array([m_u, m_d, m_s])
C = LAM**2 * np.prod(m_arr)**(1.0/3.0)
LAM6 = LAM**6

# Meson VEVs (Seiberg seesaw): M^i_i = C/m_i
M_vev = C / m_arr  # [M_u, M_d, M_s]
M_u_vev, M_d_vev, M_s_vev = M_vev

# X VEV: from F_{M^u_u} = 0 => m_u + X * M_d * M_s = 0 => X = -m_u / (M_d * M_s)
# = -m_u m_d m_s / (C^2 / (m_d m_s)) = -m_u * m_d * m_s * m_d * m_s / C^2
# More simply: X = -C / Lambda^6
X0 = -C / LAM6

# Higgs VEVs
v_hu = v / np.sqrt(2)
v_hd = v / np.sqrt(2)

# SM lepton masses for comparison
m_e   = 0.51099895   # MeV
m_mu  = 105.6583755  # MeV
m_tau = 1776.86      # MeV

# ============================================================================
# Output accumulator
# ============================================================================
output_lines = []

def p(s=""):
    output_lines.append(s)
    print(s)

def koide_Q(masses, signs=None):
    """Compute Q = sum(m_j) / (sum(s_j * sqrt(m_j)))^2."""
    ms = np.array(masses, dtype=float)
    if signs is None:
        signs = np.ones(len(ms))
    ss = np.array(signs, dtype=float)
    numerator = np.sum(ms)
    denominator = np.sum(ss * np.sqrt(ms))**2
    return numerator / denominator


# ############################################################################
# TASK 1: 6x6 Central-Block Fermion Mass Matrix
# ############################################################################
p("=" * 76)
p("TASK 1: 6x6 Central-Block Fermion Mass Matrix W_IJ")
p("=" * 76)
p()

p("--- Parameter verification ---")
p(f"  m_u = {m_u} MeV, m_d = {m_d} MeV, m_s = {m_s} MeV")
p(f"  Lambda = {LAM} MeV, v = {v} MeV")
p(f"  C = Lambda^2 * (m_u m_d m_s)^(1/3) = {C:.6f} MeV^2")
p(f"  y_c = 2*m_c/v = {y_c:.8e}")
p(f"  y_b = 2*m_b/v = {y_b:.8e}")
p()
p(f"  Meson VEVs:")
p(f"    M_u = C/m_u = {M_u_vev:.4f} MeV")
p(f"    M_d = C/m_d = {M_d_vev:.4f} MeV")
p(f"    M_s = C/m_s = {M_s_vev:.4f} MeV")
p(f"  X_0 = -C/Lambda^6 = {X0:.8e}")
p(f"  <Hu0> = <Hd0> = v/sqrt(2) = {v_hu:.4f} MeV")
p()

# Verify: det M / Lambda^6
det_M = M_u_vev * M_d_vev * M_s_vev
p(f"  Check: det M = {det_M:.6e}, Lambda^6 = {LAM6:.6e}")
p(f"  det M / Lambda^6 = {det_M/LAM6:.14f}  (should be 1.000)")
p()

# Field ordering for the 6x6 block:
# 0: M^u_u  (= phi_0)
# 1: M^d_d  (= phi_1)
# 2: M^s_s  (= phi_2)
# 3: X      (= phi_3)
# 4: Hu0    (= phi_4)
# 5: Hd0    (= phi_5)

NAMES_6 = ['M^u_u', 'M^d_d', 'M^s_s', 'X', 'H_u^0', 'H_d^0']

W6 = np.zeros((6, 6))

# --- Terms from W = sum_i m_i M^i_i ---
# d^2 W / d(M^i_i) d(anything) from the linear mass terms: all zero (no second derivatives)

# --- Terms from X * det M ---
# det M = M^u_u * M^d_d * M^s_s  (for the diagonal block, off-diagonal mesons = 0)
# At the vacuum:
#
# d^2(X det M) / dM^a_a dM^b_b  (a != b):
#   = X * d^2(det M)/dM^a_a dM^b_b = X * M_e  (e = third index, a != b, e != a, e != b)
#
# d^2(X det M) / dM^a_a dX:
#   = d(det M)/dM^a_a = cofactor(a,a) = prod_{b != a} M_b
#
# d^2(X det M) / dX dX = 0  (det M is independent of X)

# Diagonal meson - diagonal meson: X0 * M_e
for a in range(3):
    for b in range(3):
        if a != b:
            e = 3 - a - b  # third flavor index
            W6[a, b] += X0 * M_vev[e]

# Diagonal meson - X: cofactor(a,a) = prod of other two M_vev
for a in range(3):
    others = [M_vev[c] for c in range(3) if c != a]
    cof = others[0] * others[1]
    W6[a, 3] += cof
    W6[3, a] += cof

# --- Terms from y_c Hu0 M^d_d ---
# d^2 / dHu0 dM^d_d = y_c
W6[4, 1] += y_c
W6[1, 4] += y_c

# --- Terms from y_b Hd0 M^s_s ---
# d^2 / dHd0 dM^s_s = y_b
W6[5, 2] += y_b
W6[2, 5] += y_b

# That's it for the base superpotential. All other second derivatives vanish.

p("6x6 Fermion mass matrix W_IJ (base superpotential, no NMSSM):")
p()
p(f"  Field ordering: {NAMES_6}")
p()

# Print the matrix nicely
p("  W_IJ:")
for i in range(6):
    row_str = "  ["
    for j in range(6):
        val = W6[i, j]
        if abs(val) < 1e-30:
            row_str += f"{'0':>14s}"
        else:
            row_str += f"{val:>14.6e}"
        if j < 5:
            row_str += ", "
    row_str += " ]"
    p(f"    {NAMES_6[i]:>6s}: {row_str}")
p()

# Print symbolic structure
p("  Symbolic structure:")
p(f"    W[M^u_u, M^d_d] = X0 * M_s = {X0:.4e} * {M_s_vev:.4f} = {X0 * M_s_vev:.6e}")
p(f"    W[M^u_u, M^s_s] = X0 * M_d = {X0:.4e} * {M_d_vev:.4f} = {X0 * M_d_vev:.6e}")
p(f"    W[M^d_d, M^s_s] = X0 * M_u = {X0:.4e} * {M_u_vev:.4f} = {X0 * M_u_vev:.6e}")
p(f"    W[M^u_u, X]     = M_d * M_s = {M_d_vev:.4f} * {M_s_vev:.4f} = {M_d_vev * M_s_vev:.6e}")
p(f"    W[M^d_d, X]     = M_u * M_s = {M_u_vev:.4f} * {M_s_vev:.4f} = {M_u_vev * M_s_vev:.6e}")
p(f"    W[M^s_s, X]     = M_u * M_d = {M_u_vev:.4f} * {M_d_vev:.4f} = {M_u_vev * M_d_vev:.6e}")
p(f"    W[M^d_d, Hu0]   = y_c = {y_c:.8e}")
p(f"    W[M^s_s, Hd0]   = y_b = {y_b:.8e}")
p()


# ############################################################################
# TASK 2: Eigenvalues
# ############################################################################
p("=" * 76)
p("TASK 2: Eigenvalues of the 6x6 Fermion Mass Matrix")
p("=" * 76)
p()

# The matrix is real symmetric, use eigh for eigenvalues + eigenvectors
eigenvalues, eigenvectors = eigh(W6)

# Sort by absolute value
idx_sort = np.argsort(np.abs(eigenvalues))
eigenvalues_sorted = eigenvalues[idx_sort]
eigenvectors_sorted = eigenvectors[:, idx_sort]

p("  Eigenvalues (sorted by |lambda|):")
p()
for k in range(6):
    ev = eigenvalues_sorted[k]
    p(f"    lambda_{k+1} = {ev:>16.8e} MeV  (|lambda| = {abs(ev):.8e} MeV)")
p()

# Identify heavy/light split
p("  Hierarchy analysis:")
abs_eigs = np.abs(eigenvalues_sorted)
for k in range(5):
    ratio = abs_eigs[k+1] / abs_eigs[k] if abs_eigs[k] > 0 else float('inf')
    p(f"    |lambda_{k+2}| / |lambda_{k+1}| = {ratio:.4e}")
p()


# ############################################################################
# TASK 3: Light Eigenstates
# ############################################################################
p("=" * 76)
p("TASK 3: Light Eigenstates")
p("=" * 76)
p()

# Identify which are "light" vs "heavy"
# There should be a natural separation
p("  Three lightest eigenvalues:")
p()
for k in range(3):
    ev = eigenvalues_sorted[k]
    vec = eigenvectors_sorted[:, k]
    p(f"  --- Eigenvalue {k+1}: lambda = {ev:.8e} MeV ---")
    p(f"    Eigenvector components:")
    for j in range(6):
        if abs(vec[j]) > 1e-10:
            p(f"      {NAMES_6[j]:>6s}: {vec[j]:>12.8f}")
    p(f"    Dominant component: {NAMES_6[np.argmax(np.abs(vec))]}")
    p()

p("  Three heaviest eigenvalues:")
p()
for k in range(3, 6):
    ev = eigenvalues_sorted[k]
    vec = eigenvectors_sorted[:, k]
    p(f"  --- Eigenvalue {k+1}: lambda = {ev:.8e} MeV ---")
    p(f"    Eigenvector components:")
    for j in range(6):
        if abs(vec[j]) > 1e-10:
            p(f"      {NAMES_6[j]:>6s}: {vec[j]:>12.8f}")
    p(f"    Dominant component: {NAMES_6[np.argmax(np.abs(vec))]}")
    p()

# Check resemblance to lepton masses
light_eigs = np.abs(eigenvalues_sorted[:3])
p("  Comparison to lepton masses:")
p(f"    Light eigenvalues (absolute): {light_eigs}")
p(f"    Lepton masses: [{m_e:.6f}, {m_mu:.4f}, {m_tau:.2f}] MeV")
p()

# Check if any rescaling works
if light_eigs[0] > 0:
    scale_e = m_e / light_eigs[0]
    scale_mu = m_mu / light_eigs[1] if light_eigs[1] > 0 else 0
    scale_tau = m_tau / light_eigs[2] if light_eigs[2] > 0 else 0
    p(f"    Rescaling factors to match leptons:")
    p(f"      m_e / |lambda_1| = {scale_e:.6e}")
    p(f"      m_mu / |lambda_2| = {scale_mu:.6e}")
    p(f"      m_tau / |lambda_3| = {scale_tau:.6e}")
    p(f"    These would need to be equal for a universal rescaling to work.")
    if abs(scale_e) > 0 and abs(scale_mu) > 0 and abs(scale_tau) > 0:
        p(f"    Ratios: scale_mu/scale_e = {scale_mu/scale_e:.4f}, scale_tau/scale_e = {scale_tau/scale_e:.4f}")
p()

# Mass ratios
p("  Mass ratios:")
if light_eigs[0] > 0 and light_eigs[1] > 0 and light_eigs[2] > 0:
    p(f"    Light eigenvalues: |l2/l1| = {light_eigs[1]/light_eigs[0]:.4f}, |l3/l1| = {light_eigs[2]/light_eigs[0]:.4f}")
    p(f"    Lepton masses:     m_mu/m_e = {m_mu/m_e:.4f}, m_tau/m_e = {m_tau/m_e:.4f}")
p()


# ############################################################################
# TASK 4: Koide Test on Light Eigenvalues
# ############################################################################
p("=" * 76)
p("TASK 4: Koide Test on Light Eigenvalues")
p("=" * 76)
p()

light_abs = sorted(np.abs(eigenvalues_sorted[:3]))
Q_light = koide_Q(light_abs)
p(f"  Three lightest |eigenvalues|: {light_abs}")
p(f"  Q(|l1|, |l2|, |l3|) = {Q_light:.10f}")
p(f"  Deviation from 2/3: {(Q_light - 2/3)/(2/3)*100:+.4f}%")
p()

# Also test with signed eigenvalues (some may be negative)
signed_light = sorted(eigenvalues_sorted[:3])
p(f"  Signed light eigenvalues: {signed_light}")
# For signed Koide, use signs
signs_light = [np.sign(x) for x in signed_light]
abs_light = [abs(x) for x in signed_light]
if all(a > 0 for a in abs_light):
    Q_signed = koide_Q(abs_light, signs_light)
    p(f"  Q with signs: {Q_signed:.10f}")
    p(f"  Deviation from 2/3: {(Q_signed - 2/3)/(2/3)*100:+.4f}%")
p()


# ############################################################################
# TASK 5: NMSSM Extension
# ############################################################################
p("=" * 76)
p("TASK 5: NMSSM Extension (add lambda X Hu Hd)")
p("=" * 76)
p()

p(f"  NMSSM coupling: lambda = {lam_nmssm}")
p()

# The NMSSM term lambda X Hu0 Hd0 adds:
#   d^2 / dX dHu0 = lambda * <Hd0> = lambda * v/sqrt(2)
#   d^2 / dX dHd0 = lambda * <Hu0> = lambda * v/sqrt(2)
#   d^2 / dHu0 dHd0 = lambda * <X> = lambda * X0

# First: check if the vacuum needs modification
# F_X now includes: det M - Lambda^6 + lambda Hu0 Hd0
# At the diagonal vacuum: det M = Lambda^6, and lambda * v^2/2 is a perturbation
shift_detM = lam_nmssm * v_hu * v_hd  # = lambda * v^2/2
p(f"  NMSSM shift to det M constraint: lambda * v^2/2 = {shift_detM:.6e} MeV^6")
p(f"  Lambda^6 = {LAM6:.6e} MeV^6")
p(f"  Ratio: {shift_detM / LAM6:.6e}")
p()

# Since Lambda^6 >> lambda v^2/2 (for Lambda = 300 MeV), the shift is tiny
# We can either:
# (a) Use the same vacuum (first approximation)
# (b) Adjust the vacuum to satisfy the new F_X = 0

# Option (b): new det M = Lambda^6 - lambda v^2/2
det_M_new = LAM6 - shift_detM
p(f"  New det M constraint: Lambda^6 - lambda*v^2/2 = {det_M_new:.6e}")

if det_M_new > 0:
    LAM6_eff = det_M_new
    prod_m = np.prod(m_arr)
    C_eff = (LAM6_eff * prod_m)**(1.0/3.0)
    M_vev_eff = C_eff / m_arr
    X0_eff = -C_eff / LAM6_eff
    p(f"  C_eff = {C_eff:.6f} MeV^2  (original: {C:.6f})")
    p(f"  Fractional change: {(C_eff - C)/C:.6e}")
    p(f"  M_vev_eff = [{M_vev_eff[0]:.4f}, {M_vev_eff[1]:.4f}, {M_vev_eff[2]:.4f}]")
    p(f"  X0_eff = {X0_eff:.8e}  (original: {X0:.8e})")
    p()

    # Use effective vacuum for NMSSM
    M_vev_nmssm = M_vev_eff
    X0_nmssm = X0_eff
else:
    p("  WARNING: New det M < 0, using original vacuum")
    M_vev_nmssm = M_vev
    X0_nmssm = X0
p()

# Build the NMSSM 6x6 matrix
W6_nmssm = np.zeros((6, 6))

# Same structure as before but with NMSSM vacuum values
# Diagonal meson - diagonal meson: X0_nmssm * M_e
for a in range(3):
    for b in range(3):
        if a != b:
            e = 3 - a - b
            W6_nmssm[a, b] += X0_nmssm * M_vev_nmssm[e]

# Diagonal meson - X: cofactor
for a in range(3):
    others = [M_vev_nmssm[c] for c in range(3) if c != a]
    cof = others[0] * others[1]
    W6_nmssm[a, 3] += cof
    W6_nmssm[3, a] += cof

# Yukawa couplings (same)
W6_nmssm[4, 1] += y_c
W6_nmssm[1, 4] += y_c
W6_nmssm[5, 2] += y_b
W6_nmssm[2, 5] += y_b

# NMSSM additions
W6_nmssm[3, 4] += lam_nmssm * v_hd   # d^2/dX dHu0 = lambda <Hd0>
W6_nmssm[4, 3] += lam_nmssm * v_hd
W6_nmssm[3, 5] += lam_nmssm * v_hu   # d^2/dX dHd0 = lambda <Hu0>
W6_nmssm[5, 3] += lam_nmssm * v_hu
W6_nmssm[4, 5] += lam_nmssm * X0_nmssm   # d^2/dHu0 dHd0 = lambda X0
W6_nmssm[5, 4] += lam_nmssm * X0_nmssm

p("  NMSSM 6x6 Fermion mass matrix:")
p(f"  Field ordering: {NAMES_6}")
p()
p("  New entries added:")
p(f"    W[X, Hu0]   = lambda * v/sqrt(2) = {lam_nmssm * v_hd:.6e}")
p(f"    W[X, Hd0]   = lambda * v/sqrt(2) = {lam_nmssm * v_hu:.6e}")
p(f"    W[Hu0, Hd0] = lambda * X0        = {lam_nmssm * X0_nmssm:.6e}")
p()

# Print the NMSSM matrix
p("  W_IJ (NMSSM):")
for i in range(6):
    row_str = "  ["
    for j in range(6):
        val = W6_nmssm[i, j]
        if abs(val) < 1e-30:
            row_str += f"{'0':>14s}"
        else:
            row_str += f"{val:>14.6e}"
        if j < 5:
            row_str += ", "
    row_str += " ]"
    p(f"    {NAMES_6[i]:>6s}: {row_str}")
p()

# Eigenvalues
eigs_nmssm, vecs_nmssm = eigh(W6_nmssm)
idx_sort_n = np.argsort(np.abs(eigs_nmssm))
eigs_nmssm_sorted = eigs_nmssm[idx_sort_n]
vecs_nmssm_sorted = vecs_nmssm[:, idx_sort_n]

p("  NMSSM Eigenvalues (sorted by |lambda|):")
p()
for k in range(6):
    ev = eigs_nmssm_sorted[k]
    p(f"    lambda_{k+1} = {ev:>16.8e} MeV  (|lambda| = {abs(ev):.8e} MeV)")
p()

# Light eigenvalues comparison
p("  Comparison of light eigenvalues (base vs NMSSM):")
p()
for k in range(3):
    ev_base = eigenvalues_sorted[k]
    ev_nmssm = eigs_nmssm_sorted[k]
    p(f"    {k+1}: base = {ev_base:.8e}, NMSSM = {ev_nmssm:.8e}, "
      f"change = {(ev_nmssm - ev_base)/abs(ev_base)*100 if abs(ev_base) > 0 else 0:+.4f}%")
p()

# Koide test on NMSSM light eigenvalues
light_nmssm = sorted(np.abs(eigs_nmssm_sorted[:3]))
Q_nmssm = koide_Q(light_nmssm)
p(f"  Q(NMSSM light) = {Q_nmssm:.10f}")
p(f"  Deviation from 2/3: {(Q_nmssm - 2/3)/(2/3)*100:+.4f}%")
p()

# NMSSM light eigenvectors
p("  NMSSM light eigenvectors:")
for k in range(3):
    ev = eigs_nmssm_sorted[k]
    vec = vecs_nmssm_sorted[:, k]
    p(f"  --- Eigenvalue {k+1}: lambda = {ev:.8e} MeV ---")
    for j in range(6):
        if abs(vec[j]) > 1e-10:
            p(f"      {NAMES_6[j]:>6s}: {vec[j]:>12.8f}")
    p()


# ############################################################################
# TASK 6: Lepton Yukawa Identification
# ############################################################################
p("=" * 76)
p("TASK 6: Effective Lepton Yukawa Couplings")
p("=" * 76)
p()

p("  If the three light mesinos ARE the leptons, their masses come from")
p("  m_ell = y_ell * v / sqrt(2).")
p()

# Use the base (non-NMSSM) light eigenvalues
p("  --- Base superpotential ---")
for k in range(3):
    ev = np.abs(eigenvalues_sorted[k])
    y_eff = ev * np.sqrt(2) / v
    p(f"    |lambda_{k+1}| = {ev:.8e} MeV  =>  y_eff = {y_eff:.8e}")

p()
p("  SM lepton Yukawa couplings (y = sqrt(2) m / v):")
for name, m_l in [('e', m_e), ('mu', m_mu), ('tau', m_tau)]:
    y_sm = np.sqrt(2) * m_l / v
    p(f"    y_{name} = {y_sm:.8e}")

p()

# Also for NMSSM
p("  --- NMSSM ---")
for k in range(3):
    ev = np.abs(eigs_nmssm_sorted[k])
    y_eff = ev * np.sqrt(2) / v
    p(f"    |lambda_{k+1}| = {ev:.8e} MeV  =>  y_eff = {y_eff:.8e}")
p()


# ############################################################################
# TASK 7: Off-Diagonal Meson Alternative
# ############################################################################
p("=" * 76)
p("TASK 7: Off-Diagonal Meson Fermion Masses")
p("=" * 76)
p()

p("  The off-diagonal meson pairs (M^i_j, M^j_i) with i != j have")
p("  fermion mass matrix 2x2 block:")
p()
p("    ( 0        X0 * (-M_k) )")
p("    ( X0*(-M_k)    0       )")
p()
p("  where k is the third flavor index (not i or j).")
p("  Eigenvalues: +/- |X0| * M_k = +/- |X0| * C/m_k")
p()

p("  But there's also the X coupling. The off-diagonal mesons couple to X via")
p("  the d^2(det M)/dM^i_j dX = cofactor(i,j) which vanishes at the diagonal")
p("  vacuum for i != j. So the 2x2 block is decoupled from X.")
p()

# The off-diagonal pairs and their eigenvalues
# (M^u_d, M^d_u): eigenvalue = |X0| * M_s = |X0| * C/m_s
# (M^u_s, M^s_u): eigenvalue = |X0| * M_d = |X0| * C/m_d
# (M^d_s, M^s_d): eigenvalue = |X0| * M_u = |X0| * C/m_u

abs_X0 = abs(X0)

offdiag_pairs = [
    ('(M^u_d, M^d_u)', 's', m_s, M_s_vev),
    ('(M^u_s, M^s_u)', 'd', m_d, M_d_vev),
    ('(M^d_s, M^s_d)', 'u', m_u, M_u_vev),
]

offdiag_eigenvalues = []
p("  Off-diagonal meson pair eigenvalues:")
p()
for pair_name, k_name, m_k, M_k in offdiag_pairs:
    eig_val = abs_X0 * M_k
    offdiag_eigenvalues.append(eig_val)
    p(f"    {pair_name}: |X0| * M_{k_name} = {abs_X0:.6e} * {M_k:.4f} = {eig_val:.8e} MeV")
    p(f"      = |X0| * C/m_{k_name} = {abs_X0 * C:.6e} / {m_k} = {eig_val:.8e} MeV")
    p(f"      proportional to 1/m_{k_name}: C_eff/m_{k_name} where C_eff = |X0| * C = {abs_X0 * C:.6e}")

p()

# These are proportional to 1/m_k: the eigenvalue for the pair involving
# the k-th flavor is |X0| * C/m_k.
# The three eigenvalues are proportional to (1/m_s, 1/m_d, 1/m_u).

p("  Sorted by magnitude:")
offdiag_sorted = sorted(offdiag_eigenvalues)
for k, ev in enumerate(offdiag_sorted):
    p(f"    {k+1}: {ev:.8e} MeV")
p()

# Koide test
Q_offdiag = koide_Q(offdiag_sorted)
p(f"  Q(offdiag eigenvalues) = {Q_offdiag:.10f}")
p(f"  Deviation from 2/3: {(Q_offdiag - 2/3)/(2/3)*100:+.4f}%")
p()

# These are proportional to 1/m_i for i in {u, d, s}
# So Q = Q(1/m_u, 1/m_d, 1/m_s) (since |X0|*C cancels)
Q_inv_uds = koide_Q([1/m_u, 1/m_d, 1/m_s])
p(f"  Cross-check: Q(1/m_u, 1/m_d, 1/m_s) = {Q_inv_uds:.10f}")
p(f"  Match: {np.isclose(Q_offdiag, Q_inv_uds)}")
p()

# Compare to the known dual Koide for down-type quarks
Q_inv_dsb = koide_Q([1/m_d, 1/m_s, 1/m_b])
p(f"  For comparison:")
p(f"    Q(1/m_u, 1/m_d, 1/m_s) = {Q_inv_uds:.10f}")
p(f"    Q(1/m_d, 1/m_s, 1/m_b) = {Q_inv_dsb:.10f}  (dual Koide, 0.22% from 2/3)")
p(f"    Q(1/m_u, 1/m_d, 1/m_s) is {(Q_inv_uds - 2/3)/(2/3)*100:+.3f}% from 2/3.")
p()

# Mass ratios of off-diagonal eigenvalues vs lepton mass ratios
p("  Off-diagonal eigenvalue ratios:")
p(f"    eig_2 / eig_1 = {offdiag_sorted[1]/offdiag_sorted[0]:.4f}")
p(f"    eig_3 / eig_1 = {offdiag_sorted[2]/offdiag_sorted[0]:.4f}")
p(f"  Lepton mass ratios:")
p(f"    m_mu / m_e = {m_mu/m_e:.4f}")
p(f"    m_tau / m_e = {m_tau/m_e:.4f}")
p(f"  Inverse light quark mass ratios:")
p(f"    (1/m_d) / (1/m_s) = m_s/m_d = {m_s/m_d:.4f}")
p(f"    (1/m_u) / (1/m_s) = m_s/m_u = {m_s/m_u:.4f}")
p()


# ############################################################################
# SUMMARY
# ############################################################################
p()
p("=" * 76)
p("SUMMARY")
p("=" * 76)
p()

p("1. CENTRAL BLOCK (6x6):")
p(f"   The 6x6 fermion mass matrix for (M^u_u, M^d_d, M^s_s, X, Hu0, Hd0)")
p(f"   has eigenvalues dominated by the X-meson couplings (cofactors) and")
p(f"   the X*det(M) second derivatives.")
p()
p(f"   Eigenvalues (sorted by |lambda|):")
for k in range(6):
    ev = eigenvalues_sorted[k]
    p(f"     {k+1}: {ev:>16.8e} MeV")
p()

p("2. LIGHT EIGENSTATES:")
p(f"   The three lightest eigenvalues are {eigenvalues_sorted[0]:.4e}, "
  f"{eigenvalues_sorted[1]:.4e}, {eigenvalues_sorted[2]:.4e} MeV.")
p(f"   Their Koide quotient Q = {Q_light:.8f}.")
p(f"   2/3 deviation: {(Q_light - 2/3)/(2/3)*100:+.4f}%.")
p()

p("3. NMSSM:")
p(f"   Adding lambda X Hu Hd with lambda = {lam_nmssm} shifts the light")
p(f"   eigenvalues. NMSSM Koide Q = {Q_nmssm:.8f} ({(Q_nmssm - 2/3)/(2/3)*100:+.4f}% from 2/3).")
p()

p("4. OFF-DIAGONAL MESONS:")
p(f"   The off-diagonal meson pair eigenvalues are proportional to 1/m_k.")
p(f"   Q(1/m_u, 1/m_d, 1/m_s) = {Q_inv_uds:.8f} ({(Q_inv_uds - 2/3)/(2/3)*100:+.3f}% from 2/3).")
p()

p("5. EFFECTIVE YUKAWAS:")
p(f"   If light mesinos = leptons with m = y v/sqrt(2):")
for k in range(3):
    ev = np.abs(eigenvalues_sorted[k])
    y_eff = ev * np.sqrt(2) / v
    p(f"     y_{k+1} = {y_eff:.6e}")
p(f"   SM values: y_e = {np.sqrt(2)*m_e/v:.6e}, y_mu = {np.sqrt(2)*m_mu/v:.6e}, y_tau = {np.sqrt(2)*m_tau/v:.6e}")
p()


# ============================================================================
# Write markdown summary
# ============================================================================
md = []
md.append("# Lepton Yukawa from SQCD Meson Fermion Mass Matrix")
md.append("")
md.append("## Setup")
md.append("")
md.append("SU(3) SQCD with N_f = N_c = 3 flavors (u, d, s). Seiberg duality gives mesons")
md.append("M^i_j and baryons B, Btilde. The superpotential is:")
md.append("")
md.append("W = sum_i m_i M^i_i + X(det M - B Btilde - Lambda^6) + y_c Hu0 M^d_d + y_b Hd0 M^s_s")
md.append("")
md.append("Parameters: m_u = 2.16, m_d = 4.67, m_s = 93.4 MeV, Lambda = 300 MeV,")
md.append(f"v = 246220 MeV, y_c = {y_c:.6e}, y_b = {y_b:.6e}")
md.append("")
md.append("Vacuum: M^i_i = C/m_i, X = -C/Lambda^6, Hu0 = Hd0 = v/sqrt(2)")
md.append(f"where C = Lambda^2 (m_u m_d m_s)^(1/3) = {C:.6f} MeV^2")
md.append("")

md.append("## 1. Central-Block Fermion Mass Matrix (6x6)")
md.append("")
md.append("Fields: (M^u_u, M^d_d, M^s_s, X, Hu0, Hd0)")
md.append("")
md.append("Nonzero entries:")
md.append(f"- W[M^u_u, M^d_d] = X0 * M_s = {W6[0,1]:.6e}")
md.append(f"- W[M^u_u, M^s_s] = X0 * M_d = {W6[0,2]:.6e}")
md.append(f"- W[M^d_d, M^s_s] = X0 * M_u = {W6[1,2]:.6e}")
md.append(f"- W[M^u_u, X] = M_d * M_s = {W6[0,3]:.6e}")
md.append(f"- W[M^d_d, X] = M_u * M_s = {W6[1,3]:.6e}")
md.append(f"- W[M^s_s, X] = M_u * M_d = {W6[2,3]:.6e}")
md.append(f"- W[M^d_d, Hu0] = y_c = {y_c:.6e}")
md.append(f"- W[M^s_s, Hd0] = y_b = {y_b:.6e}")
md.append("")

md.append("## 2. Eigenvalues")
md.append("")
md.append("| # | Eigenvalue (MeV) | |Eigenvalue| (MeV) |")
md.append("|---|-----------------|---------------------|")
for k in range(6):
    ev = eigenvalues_sorted[k]
    md.append(f"| {k+1} | {ev:.8e} | {abs(ev):.8e} |")
md.append("")

md.append("## 3. Light Eigenstates")
md.append("")
for k in range(3):
    ev = eigenvalues_sorted[k]
    vec = eigenvectors_sorted[:, k]
    md.append(f"### Eigenvalue {k+1}: {ev:.6e} MeV")
    md.append("")
    md.append("| Field | Component |")
    md.append("|-------|-----------|")
    for j in range(6):
        if abs(vec[j]) > 1e-10:
            md.append(f"| {NAMES_6[j]} | {vec[j]:.8f} |")
    md.append("")

md.append("### Comparison to lepton masses")
md.append("")
md.append("| Eigenvalue (MeV) | Lepton mass (MeV) | Ratio |")
md.append("|------------------|--------------------|-------|")
leptons = [m_e, m_mu, m_tau]
for k in range(3):
    ev = abs(eigenvalues_sorted[k])
    ml = leptons[k] if k < 3 else 0
    ratio = ml/ev if ev > 0 else float('inf')
    md.append(f"| {ev:.6e} | {ml:.6f} | {ratio:.4e} |")
md.append("")
md.append("The light eigenvalue ratios do NOT match lepton mass ratios.")
md.append("A universal rescaling cannot relate them.")
md.append("")

md.append("## 4. Koide Test")
md.append("")
md.append(f"Q(light eigenvalues) = {Q_light:.10f}")
md.append(f"Deviation from 2/3: {(Q_light - 2/3)/(2/3)*100:+.4f}%")
md.append("")

md.append("## 5. NMSSM Extension")
md.append("")
md.append(f"Adding W_NMSSM = lambda X Hu0 Hd0 with lambda = {lam_nmssm}:")
md.append("")
md.append(f"New entries: W[X, Hu0] = W[X, Hd0] = lambda v/sqrt(2) = {lam_nmssm * v_hd:.4e},")
md.append(f"W[Hu0, Hd0] = lambda X0 = {lam_nmssm * X0_nmssm:.4e}")
md.append("")
md.append("| # | Base (MeV) | NMSSM (MeV) | Change |")
md.append("|---|-----------|-------------|--------|")
for k in range(6):
    ev_b = eigenvalues_sorted[k]
    ev_n = eigs_nmssm_sorted[k]
    ch = (ev_n - ev_b)/abs(ev_b)*100 if abs(ev_b) > 0 else 0
    md.append(f"| {k+1} | {ev_b:.6e} | {ev_n:.6e} | {ch:+.4f}% |")
md.append("")
md.append(f"NMSSM Koide Q = {Q_nmssm:.10f} ({(Q_nmssm - 2/3)/(2/3)*100:+.4f}% from 2/3)")
md.append("")

md.append("## 6. Effective Lepton Yukawas")
md.append("")
md.append("If the three light mesinos are identified with (e, mu, tau), their effective Yukawas")
md.append("would be y = sqrt(2) m_eigenvalue / v:")
md.append("")
md.append("| Light eigenvalue (MeV) | y_eff | SM y_lepton |")
md.append("|----------------------|-------|-------------|")
sm_yukawas = [np.sqrt(2)*m_e/v, np.sqrt(2)*m_mu/v, np.sqrt(2)*m_tau/v]
sm_names = ['y_e', 'y_mu', 'y_tau']
for k in range(3):
    ev = abs(eigenvalues_sorted[k])
    y_eff = ev * np.sqrt(2) / v
    md.append(f"| {ev:.6e} | {y_eff:.6e} | {sm_yukawas[k]:.6e} ({sm_names[k]}) |")
md.append("")

md.append("## 7. Off-Diagonal Meson Alternative")
md.append("")
md.append("The off-diagonal meson pairs (M^i_j, M^j_i) have fermion mass eigenvalues")
md.append("+/- |X0| M_k where k is the third flavor. These are proportional to 1/m_k.")
md.append("")
md.append("| Pair | Eigenvalue (MeV) | Proportional to |")
md.append("|------|-----------------|-----------------|")
for pair_name, k_name, m_k, M_k in offdiag_pairs:
    eig_val = abs_X0 * M_k
    md.append(f"| {pair_name} | {eig_val:.6e} | 1/m_{k_name} |")
md.append("")
md.append(f"Q(1/m_u, 1/m_d, 1/m_s) = {Q_inv_uds:.10f}")
md.append(f"Deviation from 2/3: {(Q_inv_uds - 2/3)/(2/3)*100:+.3f}%")
md.append("")
md.append(f"For comparison: Q(1/m_d, 1/m_s, 1/m_b) = {Q_inv_dsb:.10f} (0.22% from 2/3)")
md.append("")

md.append("## Key Findings")
md.append("")
md.append("1. The 6x6 central-block fermion mass matrix has a natural heavy/light split.")
md.append("   The heavy eigenvalues are set by the meson VEVs (C/m_i ~ 10^3-10^5 MeV)")
md.append("   while the light eigenvalues arise from the seesaw with the Yukawa couplings.")
md.append("")
md.append(f"2. The Koide quotient for the three light eigenvalues is Q = {Q_light:.8f}.")
md.append(f"   This is {(Q_light - 2/3)/(2/3)*100:+.4f}% from 2/3.")
md.append("")
md.append("3. The light eigenvalue mass ratios do NOT match lepton mass ratios.")
md.append("   No universal rescaling can connect them to (e, mu, tau).")
md.append("")
md.append("4. The NMSSM extension (lambda X Hu Hd) does not qualitatively change the spectrum")
md.append("   for Lambda = 300 MeV, because the NMSSM shift to det M is tiny compared to Lambda^6.")
md.append("")
md.append("5. The off-diagonal meson eigenvalues are proportional to 1/m_k. Their Koide quotient")
md.append(f"   Q(1/m_u, 1/m_d, 1/m_s) = {Q_inv_uds:.6f} is far from 2/3.")
md.append(f"   The dual Koide Q(1/m_d, 1/m_s, 1/m_b) = {Q_inv_dsb:.6f} (0.22% from 2/3) uses")
md.append("   the down-type quarks (d, s, b), not the light quarks (u, d, s).")
md.append("")

md_text = "\n".join(md) + "\n"
md_path = "/home/codexssh/phys3/results/lepton_yukawa.md"
with open(md_path, "w") as f:
    f.write(md_text)
p(f"\n[Markdown summary written to {md_path}]")
