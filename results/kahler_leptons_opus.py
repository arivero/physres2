"""
Physical fermion masses from Kahler corrections at the Seiberg vacuum.

In N_f = N_c = 3 SQCD, the superpotential fermion mass matrix W_IJ
has eigenvalues proportional to C/m_i ~ 1/m_quark. These are far too
small and too degenerate to match lepton masses.

The PHYSICAL fermion mass matrix is not W_IJ but rather:

    m_phys = e^{K/2} (K^{-1/2}) W_IJ (K^{-1/2})

where K^{IJ-bar} is the Kahler metric and K^{-1/2} factors come from
canonical normalization of the kinetic terms.

At the Seiberg vacuum M_i = C/m_i, the meson VEVs are INVERSELY
proportional to the quark masses:
    M_u >> M_d >> M_s

A Kahler metric K_II that depends on the field value will be VERY
different for different meson fields, and this field-dependent
normalization can dramatically reshape the physical mass spectrum.

This script investigates whether bion-induced or generic Kahler
corrections can produce a lepton-like spectrum with Koide Q ~ 2/3.
"""

import numpy as np
from scipy.optimize import brentq, minimize
from itertools import combinations

# ============================================================================
# Physical parameters
# ============================================================================

m_u = 2.16       # MeV
m_d = 4.67       # MeV
m_s = 93.4       # MeV
m_c = 1270.0     # MeV
m_b = 4180.0     # MeV
m_t = 172760.0   # MeV

m_e   = 0.51099895  # MeV
m_mu  = 105.6583755 # MeV
m_tau = 1776.86     # MeV

Lambda = 300.0    # MeV (confinement scale)
v_ew  = 246220.0  # MeV (electroweak VEV)

# Yukawa couplings for charm and bottom via Higgs
y_c = 2 * m_c / v_ew
y_b = 2 * m_b / v_ew

# Seiberg vacuum scale
C = Lambda**2 * (m_u * m_d * m_s)**(1.0/3.0)

# Quark masses used for the SQCD sector
m_quarks = np.array([m_u, m_d, m_s])
quark_labels = ['u', 'd', 's']

# Meson VEVs at the Seiberg vacuum: M_i = C/m_i
M_vev = C / m_quarks

# Lagrange multiplier VEV
X0 = -C / Lambda**6

# Lepton masses (target)
lepton_masses = np.array([m_e, m_mu, m_tau])

# ============================================================================
# Koide function
# ============================================================================

def koide_Q(masses, signs=None):
    """Compute Q = sum(m) / (sum(s*sqrt(m)))^2."""
    ms = np.array(masses, dtype=float)
    if signs is None:
        signs = np.ones(len(ms))
    ss = np.array(signs, dtype=float)
    num = np.sum(ms)
    den = np.sum(ss * np.sqrt(ms))**2
    if den < 1e-30:
        return np.nan
    return num / den


# ============================================================================
# Print header
# ============================================================================

print("=" * 76)
print("PHYSICAL FERMION MASSES FROM KAHLER CORRECTIONS AT THE SEIBERG VACUUM")
print("=" * 76)
print()

# ============================================================================
# PART 0: Baseline — superpotential fermion masses (no Kahler correction)
# ============================================================================

print("=" * 76)
print("PART 0: Baseline — W_IJ eigenvalues at the Seiberg vacuum")
print("=" * 76)
print()

print("Parameters:")
print(f"  m_u = {m_u} MeV,  m_d = {m_d} MeV,  m_s = {m_s} MeV")
print(f"  Lambda = {Lambda} MeV")
print(f"  C = Lambda^2 (m_u m_d m_s)^(1/3) = {C:.6f} MeV^2")
print()

print("Meson VEVs M_i = C/m_i:")
for i, (label, mi, Mi) in enumerate(zip(quark_labels, m_quarks, M_vev)):
    print(f"  M_{label} = {C:.2f}/{mi} = {Mi:.4f} MeV")
print()

print(f"  X = -C/Lambda^6 = {X0:.6e}")
print()

# Off-diagonal block eigenvalues: +-|X| * M_k for each pair
# These are the fermion masses from the 2x2 blocks
offdiag_masses = np.abs(X0) * M_vev
print("Off-diagonal block fermion masses |X|*M_k (seesaw masses):")
for label, mi, mass in zip(quark_labels, m_quarks, offdiag_masses):
    print(f"  |X|*M_{label} = {mass:.6e} MeV  (proportional to 1/m_{label})")
print()

# These masses are proportional to 1/m_quark = C/(m_quark * something)
# The hierarchy is INVERTED: heaviest for lightest quark
print("Mass hierarchy (INVERTED):")
print(f"  |X|*M_u / |X|*M_s = M_u/M_s = m_s/m_u = {m_s/m_u:.2f}")
print(f"  |X|*M_d / |X|*M_s = M_d/M_s = m_s/m_d = {m_s/m_d:.2f}")
print()

Q_seesaw = koide_Q(offdiag_masses)
print(f"Koide Q(|X|*M_u, |X|*M_d, |X|*M_s) = {Q_seesaw:.8f}")
print(f"  = Q(1/m_u, 1/m_d, 1/m_s) = {koide_Q(1.0/m_quarks):.8f}")
print(f"  Deviation from 2/3: {(Q_seesaw - 2/3)/(2/3)*100:+.4f}%")
print()

# Compare with lepton Koide
Q_lep = koide_Q(lepton_masses)
print(f"Lepton Koide Q(e, mu, tau) = {Q_lep:.10f}")
print(f"  Deviation from 2/3: {(Q_lep - 2/3)/(2/3)*100:+.6f}%")
print()

print("PROBLEM: The seesaw masses are ~10^{-4} to ~10^{-9} MeV,")
print("far too small for leptons. And Q = {:.4f}, far from 2/3.".format(Q_seesaw))
print("We need Kahler corrections to fix BOTH the scale and the hierarchy.")
print()


# ============================================================================
# PART 1: Generic Kahler correction — field-dependent metric
# ============================================================================

print("=" * 76)
print("PART 1: Kahler correction with K = Tr(M+M)(1 + c*Tr(M+M)/Lambda^2)")
print("=" * 76)
print()

print("With K = sum_i |M_i|^2 * (1 + c * sum_j |M_j|^2 / Lambda^2),")
print("the Kahler metric for the diagonal meson M_i is:")
print()
print("  K_{ii} = d^2 K / d M_i d M_i* = 1 + c*(sum_j |M_j|^2)/Lambda^2 + c*|M_i|^2/Lambda^2")
print("         = 1 + c*(Tr M+M + |M_i|^2)/Lambda^2")
print()

# At the Seiberg vacuum:
TrMM = np.sum(M_vev**2)
print(f"At the vacuum: Tr(M+M) = sum_i M_i^2 = {TrMM:.4f} MeV^2")
print(f"  M_u^2 = {M_vev[0]**2:.4f}")
print(f"  M_d^2 = {M_vev[1]**2:.4f}")
print(f"  M_s^2 = {M_vev[2]**2:.4f}")
print()

# Kahler metric for each flavor
def kahler_metric_flavor_universal(c, p=1):
    """
    K = Tr(M+M) * (1 + c * (Tr M+M)^p / Lambda^{2p})
    K_ii = 1 + c * (p+1) * (Tr M+M)^p / Lambda^{2p} + ...

    More precisely for K = sum_i |M_i|^2 + c * (sum_i |M_i|^2)^2 / Lambda^2:
    K_ii = 1 + 2c * sum_j |M_j|^2 / Lambda^2  (for i != j, partial only wrt M_i)

    Actually, let's be more careful.
    K = sum_i |M_i|^2 + c/Lambda^2 * (sum_i |M_i|^2)^2
    dK/dM_i* = M_i + 2c/Lambda^2 * (sum_j |M_j|^2) * M_i
    d^2K/(dM_i dM_i*) = 1 + 2c/Lambda^2 * (sum_j |M_j|^2) + 2c/Lambda^2 * |M_i|^2
                       = 1 + 2c/Lambda^2 * (TrMM + |M_i|^2)
    """
    K_ii = np.zeros(3)
    for i in range(3):
        K_ii[i] = 1.0 + 2*c/Lambda**2 * (TrMM + M_vev[i]**2)
    return K_ii


print("K_ii = 1 + 2c/Lambda^2 * (Tr M+M + |M_i|^2):")
print()
for c_test in [0.001, 0.01, 0.1, 1.0]:
    K_ii = kahler_metric_flavor_universal(c_test)
    print(f"  c = {c_test}:")
    for i, label in enumerate(quark_labels):
        print(f"    K_{label}{label} = {K_ii[i]:.6f}")
    print()

# Physical fermion masses: W_ij / sqrt(K_ii * K_jj)
# For the diagonal off-diagonal blocks, the fermion mass for the (ij)-(ji) pair
# is W_{ij,ji} / sqrt(K_{ij,ij} * K_{ji,ji}).
# But for diagonal mesons, the relevant mass matrix involves the central 4x4 block.
# Let's focus on the off-diagonal pair masses which give the simplest structure.

print("For the off-diagonal meson pairs (M_ij, M_ji) with i!=j:")
print("  Tree-level fermion mass = |X| * M_k  (k = third flavor)")
print("  Physical mass = |X| * M_k / sqrt(K_{ij,ij} * K_{ji,ji})")
print()
print("Since M_ij are off-diagonal mesons, their Kahler metric is just:")
print("  K_{ij,ij} = 1 + 2c/Lambda^2 * TrMM  (no extra |M_ij|^2 term at vacuum)")
print("  because <M_ij> = 0 at the diagonal vacuum.")
print()

def physical_offdiag_masses(c):
    """Physical masses of the three off-diagonal pairs."""
    # At the diagonal vacuum, off-diagonal mesons have zero VEV
    # So their Kahler metric is just 1 + 2c*TrMM/Lambda^2
    K_offdiag = 1.0 + 2*c*TrMM/Lambda**2
    if K_offdiag <= 0:
        return None
    # The tree-level mass is |X|*M_k, normalized by sqrt(K_offdiag * K_offdiag)
    return offdiag_masses / K_offdiag

print("FLAVOR-UNIVERSAL Kahler correction (same K for all off-diagonal mesons):")
print("  Physical masses = |X|*M_k / K_offdiag")
print("  This only rescales all masses by a common factor.")
print("  Q is UNCHANGED by universal rescaling.")
print()

for c_test in [0.001, 0.01, 0.1, 1.0, 10.0]:
    m_phys = physical_offdiag_masses(c_test)
    if m_phys is not None:
        Q = koide_Q(m_phys)
        print(f"  c = {c_test:8.3f}: K_offdiag = {1 + 2*c_test*TrMM/Lambda**2:.4f}, "
              f"Q = {Q:.8f} (unchanged)")

print()
print("CONCLUSION: Flavor-universal Kahler corrections cannot change the")
print("Koide ratio. We NEED flavor-dependent corrections.")
print()


# ============================================================================
# PART 2: Flavor-dependent Kahler — the key mechanism
# ============================================================================

print("=" * 76)
print("PART 2: Flavor-dependent Kahler corrections")
print("=" * 76)
print()

print("The crucial insight: for DIAGONAL meson fields M_i = M^i_i,")
print("a Kahler potential like")
print()
print("  K = sum_i |M_i|^2 * (1 + c_i * |M_i|^{2(p-1)} / Lambda^{2(p-1)})")
print()
print("gives a flavor-DEPENDENT Kahler metric:")
print()
print("  K_{ii} = 1 + p * c_i * |M_i|^{2(p-1)} / Lambda^{2(p-1)}")
print()
print("At the Seiberg vacuum M_i = C/m_i, this becomes:")
print()
print("  K_{ii} = 1 + p * c_i * (C/m_i)^{2(p-1)} / Lambda^{2(p-1)}")
print()
print("Since M_u >> M_d >> M_s, the correction is LARGEST for the LIGHTEST quark.")
print("This can dramatically reshape the physical mass spectrum.")
print()

# ── Mechanism: K_ii ~ |M_i|^{2p} ──

print("-" * 76)
print("Case: K_ii ~ 1 + c * (M_i/Lambda)^{2(p-1)}")
print()
print("Physical fermion mass for the central block diagonal entry:")
print("  The W_IJ matrix for diagonal mesons involves the 4x4 central block.")
print("  But at leading order, the dominant effect of Kahler normalization")
print("  is simply to rescale each field by 1/sqrt(K_ii).")
print()
print("  For the CW mass from the off-diagonal blocks:")
print("  m_phys(k) = |X|*M_k / sqrt(K_{ij,ij} * K_{ji,ji})")
print()
print("  For the DIAGONAL meson masses (central 4x4 block),")
print("  the entries W[M_ii, X] = Lambda^6/M_i get rescaled by")
print("  1/sqrt(K_{ii} * K_{XX}).")
print()

# Let's be precise. The off-diagonal blocks have <M_ij> = 0.
# The correction to K_{ij,ij} comes from cross-terms in K.
# For K = sum_a sum_b |M^a_b|^2 * f(|M^a_b|^2), at the diagonal vacuum:
#
# If K = sum_a |M_a|^2 + c * sum_a |M_a|^{2p}/Lambda^{2p-2} + ...
# Then for diagonal fields: K_{ii,ii} = 1 + c*p*(p)*|M_i|^{2(p-1)}/Lambda^{2(p-1)} + ...
#   (more precisely, for K_i = |M_i|^2 + c|M_i|^{2p}/Lambda^{2p-2}:
#    K_{ii} = d^2/dM_i dM_i* [|M_i|^2 + c|M_i|^{2p}/Lambda^{2p-2}]
#           = 1 + c*p^2*|M_i|^{2(p-1)}/Lambda^{2(p-1)})
#
# For off-diagonal fields M_ij with <M_ij>=0:
#   K_{ij,ij} = 1  (the correction vanishes at zero VEV for p > 1)
#   UNLESS there are cross-coupling terms in K.

# So the off-diagonal pair masses are NOT affected by single-field Kahler corrections.
# Only the DIAGONAL meson sector (the central 4x4 block) is modified.

print("IMPORTANT: For single-field Kahler K_i = |M_i|^2 + c|M_i|^{2p}/Lambda^{2p-2},")
print("the correction vanishes for off-diagonal mesons at the vacuum (since <M_ij>=0).")
print("Only DIAGONAL meson fields M_ii have nonzero VEVs and get modified K_ii.")
print()
print("So we must work with the FULL central 4x4 block (M_11, M_22, M_33, X)")
print("and include Kahler normalization for the diagonal mesons.")
print()


# ============================================================================
# PART 3: Physical mass matrix from the central block with Kahler
# ============================================================================

print("=" * 76)
print("PART 3: Central 4x4 block with Kahler normalization")
print("=" * 76)
print()

# The central 4x4 superpotential mass matrix W_IJ at the Seiberg vacuum
# (indices: M_11, M_22, M_33, X):
#
# W[M_ii, M_jj] = X * M_k  (i != j, k = third)
# W[M_ii, X] = Lambda^6/M_i
# W[X, X] = 0

def build_central_4x4():
    """Build the W_IJ matrix for (M_11, M_22, M_33, X)."""
    W = np.zeros((4, 4))
    L6 = Lambda**6

    # Off-diagonal between diagonal mesons
    W[0, 1] = W[1, 0] = X0 * M_vev[2]   # M11-M22: third = M33
    W[0, 2] = W[2, 0] = X0 * M_vev[1]   # M11-M33: third = M22
    W[1, 2] = W[2, 1] = X0 * M_vev[0]   # M22-M33: third = M11

    # Diagonal meson × X
    W[0, 3] = W[3, 0] = L6 / M_vev[0]   # M11-X
    W[1, 3] = W[3, 1] = L6 / M_vev[1]   # M22-X
    W[2, 3] = W[3, 2] = L6 / M_vev[2]   # M33-X

    return W

W_central = build_central_4x4()

print("W_IJ for the central 4x4 block (M_11, M_22, M_33, X):")
labels = ['M_uu', 'M_dd', 'M_ss', 'X']
for i in range(4):
    row = "  " + labels[i] + ":  "
    row += "  ".join(f"{W_central[i,j]:+.4e}" for j in range(4))
    print(row)
print()

evals_tree = np.linalg.eigvalsh(W_central)
print("Tree-level eigenvalues (no Kahler correction):")
for i, ev in enumerate(sorted(evals_tree)):
    print(f"  lambda_{i+1} = {ev:+.8e}")
print()

# Now apply Kahler normalization.
# The physical mass matrix is:
#   m_phys_IJ = W_IJ / sqrt(K_II * K_JJ)
# where K_II is the Kahler metric for field I.

def physical_mass_matrix(c_arr, p=2):
    """
    Compute the physical mass matrix from the central 4x4 block.

    c_arr: array of [c_u, c_d, c_s, c_X] — Kahler coefficients.
    p: power in K_i = |M_i|^2 + c_i * |M_i|^{2p}/Lambda^{2p-2}

    K_ii = 1 + c_i * p^2 * |M_i|^{2(p-1)} / Lambda^{2(p-1)}
    K_XX = 1 + c_X * p^2 * |X|^{2(p-1)} / Lambda^{2(p-1)}

    Actually, for K = |phi|^2 + c|phi|^{2p}/Lambda^{2p-2}:
    d/d(phi*) K = phi + c*p*|phi|^{2(p-1)}*phi / Lambda^{2(p-1)}
    d^2/(dphi dphi*) K = 1 + c*p^2*|phi|^{2(p-1)} / Lambda^{2(p-1)}
    """
    K_diag = np.zeros(4)
    for i in range(3):
        K_diag[i] = 1.0 + c_arr[i] * p**2 * M_vev[i]**(2*(p-1)) / Lambda**(2*(p-1))

    # For X (which has tiny VEV):
    K_diag[3] = 1.0 + c_arr[3] * p**2 * abs(X0)**(2*(p-1)) / Lambda**(2*(p-1))

    # Check positivity
    if np.any(K_diag <= 0):
        return None, K_diag

    # Physical mass matrix
    W = build_central_4x4()
    m_phys = np.zeros((4, 4))
    for i in range(4):
        for j in range(4):
            m_phys[i, j] = W[i, j] / np.sqrt(K_diag[i] * K_diag[j])

    return m_phys, K_diag


print("Structure of the Kahler metric at the Seiberg vacuum:")
print()
print(f"  For K_i = |M_i|^2 + c * |M_i|^{{2p}} / Lambda^{{2p-2}}:")
print()
print(f"  K_{{uu}} = 1 + c * p^2 * (M_u/Lambda)^{{2(p-1)}}  where M_u/Lambda = {M_vev[0]/Lambda:.4f}")
print(f"  K_{{dd}} = 1 + c * p^2 * (M_d/Lambda)^{{2(p-1)}}  where M_d/Lambda = {M_vev[1]/Lambda:.4f}")
print(f"  K_{{ss}} = 1 + c * p^2 * (M_s/Lambda)^{{2(p-1)}}  where M_s/Lambda = {M_vev[2]/Lambda:.4f}")
print()
print("The ratios are:")
print(f"  M_u/Lambda = {M_vev[0]/Lambda:.6f}  (>> 1)")
print(f"  M_d/Lambda = {M_vev[1]/Lambda:.6f}  (>> 1)")
print(f"  M_s/Lambda = {M_vev[2]/Lambda:.6f}  (~ 1)")
print()
print("For p >= 2, the Kahler correction for M_u is ENORMOUSLY larger than for M_s.")
print("This is the key: canonical normalization suppresses the M_u mass eigenvalue")
print("relative to the M_s one, potentially creating a realistic hierarchy.")
print()


# ============================================================================
# PART 4: The key insight — what happens when K_ii ~ M_i^2
# ============================================================================

print("=" * 76)
print("PART 4: The key scaling argument")
print("=" * 76)
print()

print("If the Kahler metric is dominated by the correction (c >> 1 or large M_i):")
print("  K_ii ~ c * (M_i/Lambda)^{2(p-1)}")
print()
print("Then the physical mass from the off-diagonal entries W[M_ii, M_jj] = X*M_k:")
print("  m_phys ~ X*M_k / sqrt(K_ii * K_jj)")
print("         ~ X*M_k / [c * (M_i*M_j)^{(p-1)} / Lambda^{2(p-1)}]")
print()
print("For the M_ii-X entries: W[M_ii, X] = Lambda^6/M_i:")
print("  m_phys ~ (Lambda^6/M_i) / sqrt(K_ii * K_XX)")
print("         ~ Lambda^6 / (M_i * sqrt(K_ii)) * 1/sqrt(K_XX)")
print()
print("In the large-K limit, the physical eigenvalues of the central block")
print("are controlled by the RATIOS of Kahler metric entries.")
print()
print("Most important: if K_ii dominates and K_ii ~ M_i^{2(p-1)},")
print("then the normalization factor 1/sqrt(K_ii) ~ 1/M_i^{p-1} ~ m_i^{p-1}/C^{p-1}.")
print()
print("For p=2: 1/sqrt(K_ii) ~ m_i/C. The W_IJ entry Lambda^6/M_i = Lambda^6*m_i/C")
print("  gets rescaled to ~ Lambda^6 * m_i^2/C^2 -- proportional to m_i^2!")
print()
print("This means the PHYSICAL diagonal coupling goes as m_i^2, producing the")
print("CORRECT ordering (m_u < m_d < m_s) and a much larger hierarchy than 1/m.")
print()

# Let's compute this explicitly for p=2

print("-" * 76)
print("Explicit computation for p = 2 (quartic Kahler)")
print("-" * 76)
print()
print("K = sum_i [|M_i|^2 + c * |M_i|^4 / Lambda^2]")
print("K_ii = 1 + 4c * M_i^2/Lambda^2")
print()

def scan_p2(c_val):
    """Compute physical masses for p=2, universal c."""
    c_arr = np.array([c_val, c_val, c_val, c_val])
    m_phys_mat, K_diag = physical_mass_matrix(c_arr, p=2)
    if m_phys_mat is None:
        return None, K_diag, None
    evals = np.sort(np.abs(np.linalg.eigvalsh(m_phys_mat)))
    return evals, K_diag, m_phys_mat


print(f"{'c':>12} | {'K_uu':>14} {'K_dd':>14} {'K_ss':>14} | {'|lam_1|':>14} {'|lam_2|':>14} {'|lam_3|':>14} {'|lam_4|':>14} | {'Q(1,2,3)':>10}")
print("-" * 145)

c_values = [1e-8, 1e-6, 1e-4, 1e-3, 0.01, 0.1, 0.5, 1.0, 5.0, 10.0, 100.0, 1000.0]

for c_val in c_values:
    evals, K_diag, _ = scan_p2(c_val)
    if evals is not None:
        Q3 = koide_Q(evals[:3]) if len(evals) >= 3 and np.all(evals[:3] > 0) else np.nan
        print(f"{c_val:12.2e} | {K_diag[0]:14.4e} {K_diag[1]:14.4e} {K_diag[2]:14.4e} | "
              f"{evals[0]:14.6e} {evals[1]:14.6e} {evals[2]:14.6e} {evals[3]:14.6e} | {Q3:10.6f}")

print()


# ============================================================================
# PART 5: General power p scan
# ============================================================================

print("=" * 76)
print("PART 5: General power p scan — Q(three lightest) as function of c and p")
print("=" * 76)
print()

p_values = [1, 1.5, 2, 3, 4, 5]
c_scan = np.logspace(-6, 4, 500)

best_results = []  # (p, c, Q, evals)

print(f"{'p':>4} {'c':>12} | {'Q(3 lightest)':>14} | {'|lam_1|':>14} {'|lam_2|':>14} {'|lam_3|':>14} | {'ratio 3/1':>10} {'ratio 2/1':>10}")
print("-" * 120)

for p_val in p_values:
    best_Q_dev = 1e10
    best_for_p = None

    for c_val in c_scan:
        c_arr = np.array([c_val]*4)
        m_phys_mat, K_diag = physical_mass_matrix(c_arr, p=p_val)
        if m_phys_mat is None:
            continue
        if np.any(K_diag <= 0):
            continue

        evals = np.sort(np.abs(np.linalg.eigvalsh(m_phys_mat)))
        if len(evals) < 3 or np.any(evals[:3] <= 0):
            continue

        Q = koide_Q(evals[:3])
        if np.isnan(Q):
            continue

        dev = abs(Q - 2./3.)
        if dev < best_Q_dev:
            best_Q_dev = dev
            best_for_p = (p_val, c_val, Q, evals.copy(), K_diag.copy())

    if best_for_p is not None:
        p_val, c_val, Q, evals, K_diag = best_for_p
        r31 = evals[2]/evals[0] if evals[0] > 0 else np.inf
        r21 = evals[1]/evals[0] if evals[0] > 0 else np.inf
        print(f"{p_val:4.1f} {c_val:12.4e} | {Q:14.8f} | {evals[0]:14.6e} {evals[1]:14.6e} {evals[2]:14.6e} | {r31:10.2f} {r21:10.2f}")
        best_results.append(best_for_p)

print()

# Check if any Q is close to 2/3
print("Summary of best Q values for each p:")
for p_val, c_val, Q, evals, K_diag in best_results:
    print(f"  p = {p_val:.1f}: best Q = {Q:.8f} at c = {c_val:.4e} "
          f"(deviation from 2/3: {(Q-2/3)/(2/3)*100:+.4f}%)")
print()


# ============================================================================
# PART 6: Flavor-dependent c_i — the full exploration
# ============================================================================

print("=" * 76)
print("PART 6: Flavor-dependent Kahler coefficients")
print("=" * 76)
print()

print("Allow different c_i for different meson fields.")
print("K_ii = 1 + c_i * p^2 * (M_i/Lambda)^{2(p-1)}")
print()
print("This is the most general case. We scan over (c_u, c_d, c_s) for fixed p.")
print()

# For p=2, scan over c_u, c_d, c_s independently
# The parameter space is 3D. We'll do a smart scan.

def compute_spectrum_flavor_dep(c_u, c_d, c_s, c_X=0.0, p=2):
    """Compute physical eigenvalues with flavor-dependent Kahler."""
    c_arr = np.array([c_u, c_d, c_s, c_X])
    m_phys_mat, K_diag = physical_mass_matrix(c_arr, p=p)
    if m_phys_mat is None:
        return None, None
    evals = np.sort(np.abs(np.linalg.eigvalsh(m_phys_mat)))
    return evals, K_diag


# Strategy: scan c_u/c_s and c_d/c_s ratios at fixed c_s, since Q depends on ratios
print("Strategy: fix c_s = 1, vary c_u and c_d (since Q depends on ratios)")
print("p = 2")
print()

c_s_fixed = 1.0
p_fixed = 2

# Dense scan
n_scan = 200
log_cu_range = np.linspace(-4, 4, n_scan)  # c_u/c_s from 10^-4 to 10^4
log_cd_range = np.linspace(-4, 4, n_scan)

best_Q_dev_fd = 1e10
best_fd = None

Q_grid = np.full((n_scan, n_scan), np.nan)

for i, log_cu in enumerate(log_cu_range):
    for j, log_cd in enumerate(log_cd_range):
        c_u_val = c_s_fixed * 10**log_cu
        c_d_val = c_s_fixed * 10**log_cd
        evals, K_diag = compute_spectrum_flavor_dep(c_u_val, c_d_val, c_s_fixed, p=p_fixed)
        if evals is None or len(evals) < 3:
            continue
        if np.any(evals[:3] <= 0):
            continue
        Q = koide_Q(evals[:3])
        if np.isnan(Q):
            continue
        Q_grid[i, j] = Q
        dev = abs(Q - 2./3.)
        if dev < best_Q_dev_fd:
            best_Q_dev_fd = dev
            best_fd = (c_u_val, c_d_val, c_s_fixed, Q, evals.copy(), K_diag.copy())

if best_fd is not None:
    c_u_b, c_d_b, c_s_b, Q_b, evals_b, K_b = best_fd
    print(f"Best Q = {Q_b:.8f} at:")
    print(f"  c_u/c_s = {c_u_b/c_s_b:.4e}")
    print(f"  c_d/c_s = {c_d_b/c_s_b:.4e}")
    print(f"  Eigenvalues: {evals_b}")
    if evals_b[0] > 0:
        print(f"  Mass ratios: {evals_b[1]/evals_b[0]:.4f}, {evals_b[2]/evals_b[0]:.4f}")
    print(f"  Lepton ratios: {m_mu/m_e:.4f}, {m_tau/m_e:.4f}")
    print(f"  K_diag: {K_b}")
    print()

# Now do a finer scan around the best point
if best_fd is not None:
    log_cu_best = np.log10(c_u_b/c_s_b)
    log_cd_best = np.log10(c_d_b/c_s_b)

    log_cu_fine = np.linspace(log_cu_best - 1, log_cu_best + 1, 500)
    log_cd_fine = np.linspace(log_cd_best - 1, log_cd_best + 1, 500)

    for log_cu in log_cu_fine:
        for log_cd in log_cd_fine:
            c_u_val = c_s_fixed * 10**log_cu
            c_d_val = c_s_fixed * 10**log_cd
            evals, K_diag = compute_spectrum_flavor_dep(c_u_val, c_d_val, c_s_fixed, p=p_fixed)
            if evals is None or len(evals) < 3 or np.any(evals[:3] <= 0):
                continue
            Q = koide_Q(evals[:3])
            if np.isnan(Q):
                continue
            dev = abs(Q - 2./3.)
            if dev < best_Q_dev_fd:
                best_Q_dev_fd = dev
                best_fd = (c_u_val, c_d_val, c_s_fixed, Q, evals.copy(), K_diag.copy())

    c_u_b, c_d_b, c_s_b, Q_b, evals_b, K_b = best_fd
    print(f"After fine scan, best Q = {Q_b:.10f}")
    print(f"  c_u/c_s = {c_u_b/c_s_b:.6e}")
    print(f"  c_d/c_s = {c_d_b/c_s_b:.6e}")
    print(f"  Eigenvalues: {['%.6e' % e for e in evals_b]}")
    print(f"  Deviation from 2/3: {(Q_b - 2/3)/(2/3)*100:+.6f}%")
    print()


# ============================================================================
# PART 7: Direct optimization — find c_u, c_d, c_s that give Q = 2/3
# ============================================================================

print("=" * 76)
print("PART 7: Direct optimization for Q = 2/3")
print("=" * 76)
print()

def objective_Q23(log_c, p=2):
    """Objective: minimize |Q - 2/3|^2."""
    c_u = 10**log_c[0]
    c_d = 10**log_c[1]
    c_s = 10**log_c[2]
    c_X = 0.0  # X has tiny VEV, correction negligible
    evals, K_diag = compute_spectrum_flavor_dep(c_u, c_d, c_s, c_X, p=p)
    if evals is None or len(evals) < 3 or np.any(evals[:3] <= 0):
        return 1e10
    Q = koide_Q(evals[:3])
    if np.isnan(Q):
        return 1e10
    return (Q - 2./3.)**2


# Try multiple starting points
best_opt = None
best_obj = 1e10

starting_points = [
    [0, 0, 0],
    [-2, -1, 0],
    [2, 1, 0],
    [-4, -2, 0],
    [4, 2, 0],
    [1, 0.5, 0],
    [-1, -0.5, 0],
    [3, 1.5, 0],
]

if best_fd is not None:
    starting_points.append([np.log10(c_u_b), np.log10(c_d_b), np.log10(c_s_b)])

for p_try in [2, 3, 4]:
    for x0 in starting_points:
        try:
            result = minimize(lambda x: objective_Q23(x, p=p_try), x0,
                            method='Nelder-Mead',
                            options={'xatol': 1e-12, 'fatol': 1e-16, 'maxiter': 50000})
            if result.fun < best_obj:
                best_obj = result.fun
                best_opt = (result.x, p_try, result.fun)
        except:
            pass

if best_opt is not None:
    log_c_opt, p_opt, obj_opt = best_opt
    c_u_opt = 10**log_c_opt[0]
    c_d_opt = 10**log_c_opt[1]
    c_s_opt = 10**log_c_opt[2]

    evals_opt, K_opt = compute_spectrum_flavor_dep(c_u_opt, c_d_opt, c_s_opt, 0.0, p=p_opt)
    Q_opt = koide_Q(evals_opt[:3])

    print(f"Best optimization result (p = {p_opt}):")
    print(f"  c_u = {c_u_opt:.6e}")
    print(f"  c_d = {c_d_opt:.6e}")
    print(f"  c_s = {c_s_opt:.6e}")
    print(f"  Q = {Q_opt:.12f}")
    print(f"  Deviation from 2/3: {(Q_opt - 2/3)/(2/3)*100:+.8f}%")
    print(f"  Eigenvalues: {['%.6e' % e for e in evals_opt]}")
    print(f"  K_diag: {['%.6e' % k for k in K_opt]}")
    print()

    if evals_opt[0] > 0:
        print(f"  Mass ratios (physical):")
        print(f"    m_2/m_1 = {evals_opt[1]/evals_opt[0]:.4f}  (target m_mu/m_e = {m_mu/m_e:.4f})")
        print(f"    m_3/m_1 = {evals_opt[2]/evals_opt[0]:.4f}  (target m_tau/m_e = {m_tau/m_e:.4f})")
        print(f"    m_3/m_2 = {evals_opt[2]/evals_opt[1]:.4f}  (target m_tau/m_mu = {m_tau/m_mu:.4f})")
    print()

print()


# ============================================================================
# PART 8: The simplest scenario — large c limit
# ============================================================================

print("=" * 76)
print("PART 8: Large-c limit (Kahler dominates) — analytical")
print("=" * 76)
print()

print("When K_ii >> 1 for all diagonal mesons (large c), the normalization")
print("factor is 1/sqrt(K_ii) ~ 1/(sqrt(c)*p*(M_i/Lambda)^{p-1}).")
print()
print("The physical mass matrix becomes (for p=2, c universal):")
print()
print("  m_phys[i,j] = W[i,j] / (2c * M_i * M_j / Lambda^2)")
print("  m_phys[i,3] = W[i,X] / (2c^{1/2} * M_i / Lambda)  (since K_XX ~ 1)")
print()
print("W[M_ii, M_jj] = X*M_k, so:")
print("  m_phys[i,j] ~ X*M_k*Lambda^2 / (2c*M_i*M_j) = X*Lambda^2*C/(2c*m_k*M_i*M_j)")
print("  Since M_i = C/m_i: = X*Lambda^2*C/(2c*(C/m_i)*(C/m_j)) * C/m_k")
print("  = X*Lambda^2*m_i*m_j / (2c*C*m_k)")
print()
print("W[M_ii, X] = Lambda^6/M_i = Lambda^6*m_i/C, so:")
print("  m_phys[i,X] ~ Lambda^6*m_i / (C * 2*sqrt(c)*M_i/Lambda)")
print("  = Lambda^7*m_i^2 / (2*sqrt(c)*C^2)")
print()
print("The key entry: m_phys[i,X] ~ m_i^2 / C^2 (up to common factors)")
print()
print("In the large-c limit, the physical M_ii-X coupling goes as m_i^2.")
print("This is a DRAMATIC reshaping compared to the tree-level 1/m_i.")
print()

# Let's compute the physical matrix explicitly in the large-c limit
print("Physical 4x4 matrix in the large-c limit (p=2, c >> 1, K_XX ~ 1):")
print()

c_large = 1e6
c_arr_large = np.array([c_large]*3 + [0])  # no correction for X (tiny VEV)
m_phys_large, K_large = physical_mass_matrix(c_arr_large, p=2)

if m_phys_large is not None:
    print("  Kahler metric:")
    for i, label in enumerate(labels):
        print(f"    K_{label} = {K_large[i]:.6e}")
    print()

    print("  Physical mass matrix:")
    for i in range(4):
        row = "  " + labels[i] + ":  "
        row += "  ".join(f"{m_phys_large[i,j]:+.4e}" for j in range(4))
        print(row)
    print()

    evals_large = np.linalg.eigvalsh(m_phys_large)
    print("  Eigenvalues:")
    for i, ev in enumerate(sorted(evals_large)):
        print(f"    lambda_{i+1} = {ev:+.8e}")
    print()

    abs_evals = np.sort(np.abs(evals_large))
    print("  |Eigenvalues| (sorted):")
    for i, ev in enumerate(abs_evals):
        print(f"    |lambda_{i+1}| = {ev:.8e}")
    print()

    Q_large = koide_Q(abs_evals[:3])
    print(f"  Q(three lightest) = {Q_large:.8f}")
    print(f"  Deviation from 2/3: {(Q_large - 2/3)/(2/3)*100:+.4f}%")
    print()

    if abs_evals[0] > 0:
        print(f"  Mass ratios:")
        print(f"    m_2/m_1 = {abs_evals[1]/abs_evals[0]:.4f}")
        print(f"    m_3/m_1 = {abs_evals[2]/abs_evals[0]:.4f}")
        print()

print()


# ============================================================================
# PART 9: The m_i^2 spectrum — direct Koide check
# ============================================================================

print("=" * 76)
print("PART 9: Koide ratio of m_quark^2 spectrum")
print("=" * 76)
print()

print("In the large-c limit (p=2), the physical M_ii-X entries go as m_i^2.")
print("If these dominate the eigenvalue structure, the physical masses")
print("would be proportional to m_i^2.")
print()
print("Let's check Q for various powers of the quark masses:")
print()

for power in [0.5, 1, 1.5, 2, 3, 4, -0.5, -1, -1.5, -2]:
    m_powered = m_quarks**power
    Q_p = koide_Q(m_powered)
    print(f"  Q(m^{power:+.1f}) = {Q_p:.8f}  "
          f"(m_u^p={m_powered[0]:.4e}, m_d^p={m_powered[1]:.4e}, m_s^p={m_powered[2]:.4e})"
          f"  dev from 2/3: {(Q_p-2/3)/(2/3)*100:+.3f}%")

print()

# Also check with signed square roots (the (-s,c,b) convention)
print("Interesting: Q(m^2) for (u,d,s) quarks =", koide_Q(m_quarks**2))
print("Compare: Q(m) for leptons =", koide_Q(lepton_masses))
print()

# Check if there's a power alpha such that Q(m_quark^alpha) = 2/3
print("Scanning for alpha such that Q(m_quark^alpha) = 2/3:")

alpha_scan = np.linspace(-3, 5, 100000)
Q_alpha = np.array([koide_Q(m_quarks**alpha) for alpha in alpha_scan])
valid = np.isfinite(Q_alpha)

# Find crossings of Q = 2/3
crossings = []
for i in range(len(alpha_scan)-1):
    if valid[i] and valid[i+1]:
        if (Q_alpha[i] - 2./3.) * (Q_alpha[i+1] - 2./3.) < 0:
            alpha_cross = brentq(lambda a: koide_Q(m_quarks**a) - 2./3.,
                                alpha_scan[i], alpha_scan[i+1])
            crossings.append(alpha_cross)

for ac in crossings:
    m_pow = m_quarks**ac
    Q_check = koide_Q(m_pow)
    print(f"  alpha = {ac:.8f}: Q = {Q_check:.10f}")
    print(f"    m_u^alpha = {m_pow[0]:.6e}")
    print(f"    m_d^alpha = {m_pow[1]:.6e}")
    print(f"    m_s^alpha = {m_pow[2]:.6e}")
    print(f"    Ratios: m_d^a/m_u^a = {m_pow[1]/m_pow[0]:.4f}, m_s^a/m_u^a = {m_pow[2]/m_pow[0]:.4f}")
    print(f"    Compare lepton: m_mu/m_e = {m_mu/m_e:.4f}, m_tau/m_e = {m_tau/m_e:.4f}")
    print()

print()


# ============================================================================
# PART 10: Numerical optimization — match lepton mass ratios
# ============================================================================

print("=" * 76)
print("PART 10: Can Kahler corrections match lepton mass ratios?")
print("=" * 76)
print()

print("Target: physical mass ratios m_mu/m_e = {:.4f}, m_tau/m_e = {:.4f}".format(m_mu/m_e, m_tau/m_e))
print("  with Q(m_e, m_mu, m_tau) = {:.10f}".format(koide_Q(lepton_masses)))
print()

# Optimize c_u, c_d, c_s to match lepton mass ratios
def objective_leptons(log_c, p=2):
    """Minimize (Q - 2/3)^2 + penalty for wrong ratios."""
    c_u = 10**log_c[0]
    c_d = 10**log_c[1]
    c_s = 10**log_c[2]
    evals, K_diag = compute_spectrum_flavor_dep(c_u, c_d, c_s, 0.0, p=p)
    if evals is None or len(evals) < 3 or np.any(evals[:3] <= 0):
        return 1e10

    Q = koide_Q(evals[:3])
    if np.isnan(Q):
        return 1e10

    # Match ratios
    r21 = evals[1]/evals[0]
    r31 = evals[2]/evals[0]
    r21_target = m_mu/m_e
    r31_target = m_tau/m_e

    return (Q - 2./3.)**2 + 0.01*(np.log(r21/r21_target))**2 + 0.01*(np.log(r31/r31_target))**2


best_lep = None
best_lep_obj = 1e10

for p_try in [2, 3, 4, 5]:
    for x0 in starting_points:
        try:
            result = minimize(lambda x: objective_leptons(x, p=p_try), x0,
                            method='Nelder-Mead',
                            options={'xatol': 1e-12, 'fatol': 1e-16, 'maxiter': 100000})
            if result.fun < best_lep_obj:
                best_lep_obj = result.fun
                best_lep = (result.x, p_try, result.fun)
        except:
            pass

if best_lep is not None:
    log_c_lep, p_lep, obj_lep = best_lep
    c_u_lep = 10**log_c_lep[0]
    c_d_lep = 10**log_c_lep[1]
    c_s_lep = 10**log_c_lep[2]

    evals_lep, K_lep = compute_spectrum_flavor_dep(c_u_lep, c_d_lep, c_s_lep, 0.0, p=p_lep)
    Q_lep_result = koide_Q(evals_lep[:3])

    print(f"Best match to lepton ratios + Q = 2/3 (p = {p_lep}):")
    print(f"  c_u = {c_u_lep:.6e}")
    print(f"  c_d = {c_d_lep:.6e}")
    print(f"  c_s = {c_s_lep:.6e}")
    print(f"  Q = {Q_lep_result:.10f}")
    print(f"  Deviation from 2/3: {(Q_lep_result - 2/3)/(2/3)*100:+.8f}%")
    print(f"  Eigenvalues: {['%.6e' % e for e in evals_lep]}")
    print()

    if evals_lep[0] > 0:
        r21 = evals_lep[1]/evals_lep[0]
        r31 = evals_lep[2]/evals_lep[0]
        print(f"  Mass ratio m_2/m_1 = {r21:.4f}  (target: {m_mu/m_e:.4f})")
        print(f"  Mass ratio m_3/m_1 = {r31:.4f}  (target: {m_tau/m_e:.4f})")
        print(f"  Mass ratio m_3/m_2 = {evals_lep[2]/evals_lep[1]:.4f}  (target: {m_tau/m_mu:.4f})")
    print()

print()


# ============================================================================
# PART 11: Analytical structure — what controls Q?
# ============================================================================

print("=" * 76)
print("PART 11: Analytical structure — how Kahler reshapes the spectrum")
print("=" * 76)
print()

print("Consider the physical mass matrix in the large-c limit (p=2, universal c).")
print("Define the normalized field phi_i = M_i / sqrt(K_ii) where K_ii = 4c*M_i^2/Lambda^2.")
print("Then phi_i = Lambda / (2*sqrt(c)).")
print()
print("ALL normalized fields have the SAME magnitude in the large-c limit!")
print("This means the diagonal mesons become DEGENERATE after Kahler normalization.")
print()
print("But the off-diagonal entries W[M_ii, M_jj] = X*M_k get rescaled differently:")
print("  m_phys[i,j] = X*M_k / sqrt(K_ii*K_jj) = X*M_k*Lambda^2/(4c*M_i*M_j)")
print("             = X*Lambda^2/(4c) * M_k/(M_i*M_j)")
print("             = X*Lambda^2/(4c) * (C/m_k)/((C/m_i)*(C/m_j))")
print("             = X*Lambda^2*m_i*m_j / (4c*C*m_k)")
print()
print("And the M_ii-X entries:")
print("  m_phys[i,X] = Lambda^6/(M_i*sqrt(K_ii)) = Lambda^6*m_i/(C*2*sqrt(c)*M_i/Lambda)")
print("              = Lambda^7*m_i^2 / (2*sqrt(c)*C^2)")
print()
print("So in the large-c, p=2 limit:")
print()

# Compute the explicit large-c physical matrix
print("  Let A = X*Lambda^2/(4c*C), B = Lambda^7/(2*sqrt(c)*C^2)")
print()
print("  Physical matrix:")
print("  [   0        A*m_u*m_d/m_s  A*m_u*m_s/m_d  B*m_u^2 ]")
print("  [ A*m_u*m_d/m_s  0          A*m_d*m_s/m_u  B*m_d^2 ]")
print("  [ A*m_u*m_s/m_d  A*m_d*m_s/m_u  0          B*m_s^2 ]")
print("  [ B*m_u^2        B*m_d^2        B*m_s^2    0       ]")
print()

# In this limit, B >> A (because Lambda^7/C^2 >> X*Lambda^2/C)
# X = -C/Lambda^6, so A = -Lambda^2/(4c*Lambda^6) * Lambda^2 = -1/(4c*Lambda^4) ... let me compute
A_val = X0 * Lambda**2 / (4*c_large*C)
B_val = Lambda**7 / (2*np.sqrt(c_large)*C**2)

print(f"  For c = {c_large:.0e}:")
print(f"    A = {A_val:.6e}")
print(f"    B = {B_val:.6e}")
print(f"    B/A = {B_val/A_val:.4e}")
print()
print(f"  The M_ii-X entries (B*m_i^2) dominate over the M_ii-M_jj entries (A*m_i*m_j/m_k)")
print(f"  when B*m_i^2 >> A*m_i*m_j/m_k, i.e., B/A >> m_j/(m_i*m_k).")
print(f"  Since B/A = {B_val/abs(A_val):.4e}, this is satisfied for ALL entries.")
print()
print("  In this regime, the 4x4 matrix is dominated by the M_ii-X couplings.")
print("  The three eigenvalues of the 3+1 system with the column vector")
print("  (B*m_u^2, B*m_d^2, B*m_s^2) coupling to X are determined by the secular equation.")
print()

# The dominant 3+1 block: X couples to the three M_ii with strengths B*m_i^2
# This is a rank-1 perturbation on the zero block.
# Eigenvalues: one large (~ B * sqrt(m_u^4 + m_d^4 + m_s^4))
# and two zero (in the B >> A limit)
# So the three lightest come from the A-corrections and are much smaller.

print("  In the extreme large-c limit, the X coupling produces ONE large eigenvalue")
print("  and leaves TWO exactly zero (from the rank-1 structure).")
print("  The off-diagonal A entries then split these two into nonzero values.")
print()
print("  This means the three lightest eigenvalues are NOT simply proportional")
print("  to m_i^2. The structure is more complex.")
print()

# Let's verify by computing at several c values
print("Evolution of eigenvalues with increasing c (p=2, universal):")
print()
print(f"  {'c':>12} | {'|lam_1|':>14} {'|lam_2|':>14} {'|lam_3|':>14} {'|lam_4|':>14} | "
      f"{'Q(1,2,3)':>10} | {'r21':>8} {'r31':>8}")
print("-" * 100)

for c_test in [1e-2, 1e-1, 1, 10, 100, 1e3, 1e4, 1e5, 1e6, 1e8, 1e10]:
    evals_t, _, _ = scan_p2(c_test)
    if evals_t is not None and np.all(evals_t[:3] > 0):
        Q_t = koide_Q(evals_t[:3])
        r21 = evals_t[1]/evals_t[0] if evals_t[0] > 0 else 0
        r31 = evals_t[2]/evals_t[0] if evals_t[0] > 0 else 0
        print(f"  {c_test:12.2e} | {evals_t[0]:14.6e} {evals_t[1]:14.6e} "
              f"{evals_t[2]:14.6e} {evals_t[3]:14.6e} | {Q_t:10.6f} | {r21:8.2f} {r31:8.2f}")

print()


# ============================================================================
# PART 12: Non-universal Kahler — independent c_i optimization
# ============================================================================

print("=" * 76)
print("PART 12: Systematic search — flavor-dependent c_i, multiple p")
print("=" * 76)
print()

# Optimize ONLY for Q = 2/3 (the mass ratios are a separate issue)
print("Search for Q = 2/3 with flavor-dependent coefficients:")
print()

results_table = []

for p_try in [1.5, 2, 2.5, 3, 4, 5]:
    best_for_p = None
    best_dev_p = 1e10

    # Systematic grid + optimization
    for log_cu in np.linspace(-6, 6, 50):
        for log_cd in np.linspace(-6, 6, 50):
            c_u_try = 10**log_cu
            c_d_try = 10**log_cd
            c_s_try = 1.0

            evals, K_diag = compute_spectrum_flavor_dep(c_u_try, c_d_try, c_s_try, 0.0, p=p_try)
            if evals is None or len(evals) < 3 or np.any(evals[:3] <= 0):
                continue
            Q = koide_Q(evals[:3])
            if np.isnan(Q):
                continue
            dev = abs(Q - 2./3.)
            if dev < best_dev_p:
                best_dev_p = dev
                best_for_p = (c_u_try, c_d_try, c_s_try, Q, evals.copy(), K_diag.copy())

    # Refine with optimization from best grid point
    if best_for_p is not None:
        x0_refine = [np.log10(best_for_p[0]), np.log10(best_for_p[1]), np.log10(best_for_p[2])]
        try:
            result = minimize(lambda x: objective_Q23(x, p=p_try), x0_refine,
                            method='Nelder-Mead',
                            options={'xatol': 1e-14, 'fatol': 1e-18, 'maxiter': 100000})
            c_u_r = 10**result.x[0]
            c_d_r = 10**result.x[1]
            c_s_r = 10**result.x[2]
            evals_r, K_r = compute_spectrum_flavor_dep(c_u_r, c_d_r, c_s_r, 0.0, p=p_try)
            if evals_r is not None and len(evals_r) >= 3 and np.all(evals_r[:3] > 0):
                Q_r = koide_Q(evals_r[:3])
                if not np.isnan(Q_r) and abs(Q_r - 2./3.) < best_dev_p:
                    best_for_p = (c_u_r, c_d_r, c_s_r, Q_r, evals_r.copy(), K_r.copy())
        except:
            pass

    if best_for_p is not None:
        c_u_f, c_d_f, c_s_f, Q_f, evals_f, K_f = best_for_p
        r21_f = evals_f[1]/evals_f[0] if evals_f[0] > 0 else 0
        r31_f = evals_f[2]/evals_f[0] if evals_f[0] > 0 else 0
        results_table.append((p_try, c_u_f, c_d_f, c_s_f, Q_f, evals_f, r21_f, r31_f))

print(f"{'p':>4} | {'c_u':>12} {'c_d':>12} {'c_s':>12} | {'Q':>14} | {'dev %':>10} | {'r21':>8} {'r31':>8}")
print("-" * 100)

for p, cu, cd, cs, Q, ev, r21, r31 in results_table:
    dev_pct = (Q - 2./3.)/(2./3.)*100
    print(f"{p:4.1f} | {cu:12.4e} {cd:12.4e} {cs:12.4e} | {Q:14.10f} | {dev_pct:+10.6f} | {r21:8.2f} {r31:8.2f}")

print()
print("Comparison with lepton ratios:")
print(f"  m_mu/m_e = {m_mu/m_e:.4f}")
print(f"  m_tau/m_e = {m_tau/m_e:.4f}")
print()


# ============================================================================
# PART 13: The deep insight — K_ii ~ M_i^{2n} gives m_phys ~ m_quark^{2n}
# ============================================================================

print("=" * 76)
print("PART 13: Analytical formula — physical mass as function of Kahler power")
print("=" * 76)
print()

print("In the limit where the Kahler correction dominates (K_ii >> 1),")
print("with K_ii ~ (M_i/Lambda)^{2(p-1)}:")
print()
print("  Normalization factor: 1/sqrt(K_ii) ~ (Lambda/M_i)^{p-1} = (m_i*Lambda/C)^{p-1}")
print()
print("  The W[M_ii, X] = Lambda^6/M_i = Lambda^6*m_i/C entry becomes:")
print("  m_phys[i,X] ~ (Lambda^6*m_i/C) * (m_i*Lambda/C)^{p-1} * (1/sqrt(K_XX))")
print("              ~ Lambda^{5+p} * m_i^p / C^p")
print("              ~ m_i^p  (up to common factor)")
print()
print("So the physical M_ii-X coupling goes as m_i^p.")
print()
print("For what p does Q(m_u^p, m_d^p, m_s^p) = 2/3?")
print()

# We already scanned this. Let's be more precise.
for ac in crossings:
    m_pow = m_quarks**ac
    Q_check = koide_Q(m_pow)
    print(f"  alpha = {ac:.10f}: Q(m^alpha) = {Q_check:.12f}")
    # Compare with lepton ratios
    if m_pow[0] > 0:
        print(f"    Ratios: {m_pow[1]/m_pow[0]:.4f}, {m_pow[2]/m_pow[0]:.4f}")
        print(f"    Lepton: {m_mu/m_e:.4f}, {m_tau/m_e:.4f}")
    print()

print("The values of alpha where Q(m_u^alpha, m_d^alpha, m_s^alpha) = 2/3 are found above.")
print("These correspond to the Kahler power p = alpha.")
print()
print("NOTE: These are the powers in the large-K limit. At finite c, the")
print("actual eigenvalues of the 4x4 matrix differ from the simple m_i^p scaling")
print("because the off-diagonal M_ii-M_jj entries contribute and because the")
print("4x4 structure introduces mixing between the M_ii modes and X.")
print()


# ============================================================================
# PART 14: Off-diagonal block with diagonal-meson Kahler
# ============================================================================

print("=" * 76)
print("PART 14: Full 13x13 physical spectrum with Kahler")
print("=" * 76)
print()

print("Including off-diagonal mesons, baryons, and the Higgs field.")
print("Off-diagonal mesons M_ij at the vacuum have <M_ij> = 0,")
print("so their Kahler metric K_{ij,ij} = 1 (no correction).")
print()
print("Their physical masses = |X|*M_k (unchanged by Kahler).")
print("These are the seesaw masses, proportional to 1/m_k.")
print()
print("The FULL physical spectrum has:")
print("  - 3 pairs from off-diagonal blocks: masses ~ 1/m_k (UNCHANGED)")
print("  - 1 pair from baryon block: mass ~ |X| (UNCHANGED)")
print("  - 4 from the central block (MODIFIED by Kahler)")
print("  - 1 zero from H (or modified by Yukawa)")
print()
print("Only the 4 central-block masses are affected by diagonal-meson Kahler.")
print("The lepton-like spectrum must come from the central block.")
print()

# Build the full 13x13 physical mass matrix
def build_full_physical_13x13(c_diag, p=2, y_arr=None):
    """
    Build the full 13x13 physical fermion mass matrix with Kahler normalization.

    Fields: M11, M12, M13, M21, M22, M23, M31, M32, M33, X, B, Bbar, H

    Kahler metric:
    - Diagonal mesons M_ii: K_ii = 1 + c * p^2 * M_i^{2(p-1)} / Lambda^{2(p-1)}
    - Off-diagonal mesons M_ij: K_{ij} = 1 (zero VEV)
    - X: K_XX = 1 (tiny VEV, negligible correction)
    - B, Bbar: K = 1
    - H: K = 1
    """
    if y_arr is None:
        y_arr = np.zeros(3)

    # Build the W_IJ matrix
    W = np.zeros((13, 13))
    L6 = Lambda**6
    M1, M2, M3 = M_vev

    # Diagonal-diagonal meson couplings
    W[0, 4] = W[4, 0] = X0 * M3
    W[0, 8] = W[8, 0] = X0 * M2
    W[4, 8] = W[8, 4] = X0 * M1

    # Off-diagonal meson pairs
    W[1, 3] = W[3, 1] = X0 * M3
    W[2, 6] = W[6, 2] = X0 * M2
    W[5, 7] = W[7, 5] = X0 * M1

    # Diagonal meson x X
    W[0, 9] = W[9, 0] = L6 / M1
    W[4, 9] = W[9, 4] = L6 / M2
    W[8, 9] = W[9, 8] = L6 / M3

    # Baryons
    W[10, 11] = W[11, 10] = -X0

    # Higgs
    W[12, 0] = W[0, 12] = y_arr[0]
    W[12, 4] = W[4, 12] = y_arr[1]
    W[12, 8] = W[8, 12] = y_arr[2]

    # Kahler metric (diagonal)
    K = np.ones(13)
    for i_diag, field_idx in enumerate([0, 4, 8]):  # M11, M22, M33
        K[field_idx] = 1.0 + c_diag * p**2 * M_vev[i_diag]**(2*(p-1)) / Lambda**(2*(p-1))

    if np.any(K <= 0):
        return None, None, W, K

    # Physical mass matrix
    m_phys = np.zeros((13, 13))
    for i in range(13):
        for j in range(13):
            m_phys[i, j] = W[i, j] / np.sqrt(K[i] * K[j])

    evals = np.linalg.eigvalsh(m_phys)
    return m_phys, evals, W, K


# Compute for several c values
print("Full 13x13 physical spectrum (p=2, universal c on diagonal mesons):")
print()
print(f"  {'c':>10} | " + " ".join(f"{'|lam_'+str(i+1)+'|':>14}" for i in range(6)) + " | Q(1,2,3)")
print("-" * 120)

for c_test in [0, 1e-4, 1e-2, 1, 100, 1e4, 1e6, 1e8]:
    _, evals_full, _, _ = build_full_physical_13x13(c_test, p=2)
    if evals_full is not None:
        abs_evals = np.sort(np.abs(evals_full))
        # Skip exact zeros
        nonzero = abs_evals[abs_evals > 1e-15]
        Q_full = koide_Q(nonzero[:3]) if len(nonzero) >= 3 else np.nan
        print(f"  {c_test:10.2e} | " +
              " ".join(f"{nonzero[i]:14.6e}" if i < len(nonzero) else f"{'---':>14}"
                       for i in range(6)) + f" | {Q_full:8.6f}")

print()


# ============================================================================
# PART 15: Summary and assessment
# ============================================================================

print("=" * 76)
print("SUMMARY AND ASSESSMENT")
print("=" * 76)
print()

print("1. BASELINE PROBLEM:")
print(f"   At the Seiberg vacuum M_i = C/m_i with C = {C:.2f} MeV^2,")
print(f"   the superpotential fermion masses are proportional to 1/m_quark.")
print(f"   The eigenvalues span {offdiag_masses[2]:.2e} to {offdiag_masses[0]:.2e} MeV.")
print(f"   Q(seesaw masses) = {Q_seesaw:.6f}, far from 2/3.")
print()

print("2. FLAVOR-UNIVERSAL KAHLER:")
print("   A universal correction K = Tr(M+M)(1 + c*Tr(M+M)/Lambda^2)")
print("   rescales all masses by a common factor. Q is UNCHANGED.")
print("   This cannot work.")
print()

print("3. FLAVOR-DEPENDENT KAHLER — THE KEY MECHANISM:")
print("   When K_ii depends on the FIELD VALUE M_i at the vacuum,")
print("   the canonical normalization is different for each flavor.")
print("   Since M_u >> M_d >> M_s at the Seiberg vacuum,")
print("   the correction is largest for the lightest quark (u).")
print()

print("4. LARGE-c SCALING (p=2):")
print("   Physical M_ii-X coupling goes as m_quark^2.")
print("   This INVERTS the 1/m_quark hierarchy to the correct m_quark^2 ordering.")
print("   However, Q(m_u^2, m_d^2, m_s^2) must be checked.")
print()

Q_m2 = koide_Q(m_quarks**2)
print(f"5. Q(m_u^2, m_d^2, m_s^2) = {Q_m2:.8f}")
print(f"   Deviation from 2/3: {(Q_m2 - 2/3)/(2/3)*100:+.3f}%")
print()

if crossings:
    print("6. CRITICAL KAHLER POWERS:")
    for ac in crossings:
        print(f"   Q(m^alpha) = 2/3 at alpha = {ac:.8f}")
        m_pow = m_quarks**ac
        print(f"     Ratios: m_d^a/m_u^a = {m_pow[1]/m_pow[0]:.4f}, m_s^a/m_u^a = {m_pow[2]/m_pow[0]:.4f}")
        print(f"     Compare lepton: {m_mu/m_e:.4f}, {m_tau/m_e:.4f}")
    print()

print("7. FULL 4x4 BLOCK STRUCTURE:")
print("   The actual eigenvalues of the 4x4 central block are not simply m_i^p")
print("   because of mixing between M_ii and X modes and off-diagonal M_ii-M_jj")
print("   entries. Numerical optimization finds:")
if results_table:
    best_res = min(results_table, key=lambda x: abs(x[4] - 2./3.))
    print(f"   Best Q = {best_res[4]:.10f} at p = {best_res[0]}")
    print(f"   with c_u = {best_res[1]:.4e}, c_d = {best_res[2]:.4e}, c_s = {best_res[3]:.4e}")
    print(f"   Mass ratio m_2/m_1 = {best_res[6]:.2f}, m_3/m_1 = {best_res[7]:.2f}")
print()

print("8. ASSESSMENT:")
print("   Flavor-dependent Kahler corrections CAN produce Q = 2/3 for the")
print("   physical fermion masses at the Seiberg vacuum. The mechanism is:")
print("   (a) The seesaw vacuum M_i = C/m_i creates a huge hierarchy in field VEVs.")
print("   (b) A Kahler potential K ~ |M|^{2p} gives K_ii ~ (C/m_i)^{2(p-1)},")
print("       which is LARGEST for the LIGHTEST quark.")
print("   (c) Canonical normalization 1/sqrt(K_ii) suppresses the M_u mode most,")
print("       converting the 1/m_i seesaw spectrum to an m_i^p spectrum.")
print("   (d) For the right power p, Q(m^p) = 2/3.")
print()
print("   The mass RATIOS, however, do not match (e, mu, tau) for (u, d, s) inputs.")
print("   The quark mass hierarchy (m_s/m_u ~ 43) is too compressed compared to")
print("   the lepton hierarchy (m_tau/m_e ~ 3477).")
print()
print("   For the ratios to match, one would need either:")
print("   (i)  A very large power p ~ {:.1f}".format(np.log(m_tau/m_e)/np.log(m_s/m_u)) + " to amplify the hierarchy,")
print("   (ii) Additional structure in the superpotential (Yukawa couplings), or")
print("   (iii) A different choice of quark flavors for the SQCD sector.")
print()

print("=" * 76)
print("DONE")
print("=" * 76)
