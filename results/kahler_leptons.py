"""
Kähler potential corrections to fermion masses in SQCD with N_f = N_c = 3.

Superpotential:
    W = sum_i m_i M^i_i + X(det M - Lambda^6) + y_c H_u M^d_d + y_b H_d M^s_s

At the Seiberg vacuum: M_i^vac = Lambda^2 / m_i (electric frame).

The physical fermion mass matrix is:
    m_phys = e^{K/2} K^{-1/2}_{IJ} W_{IJ} K^{-1/2}_{IJ}

This script investigates whether bion-induced Kähler corrections can
generate a lepton-like spectrum from the tree-level W_IJ matrix.

Numerical inputs (MeV):
    m_u = 2.16, m_d = 4.67, m_s = 93.4
    Lambda = 300
    v = 246220
    y_c = 2*1270/v, y_b = 2*4180/v
"""

import numpy as np
from scipy.linalg import eig, eigh, sqrtm
from itertools import product as iproduct

# ============================================================
# Physical parameters
# ============================================================

m_u = 2.16       # MeV
m_d = 4.67       # MeV
m_s = 93.4       # MeV
m_c = 1270.0     # MeV  (charm quark mass)
m_b = 4180.0     # MeV  (bottom quark mass)
m_t = 173000.0   # MeV  (top quark mass, for reference)
Lambda = 300.0   # MeV  (dynamical scale)
v = 246220.0     # MeV  (Higgs VEV)

y_c = 2.0 * m_c / v   # ~= 2*1270/246220
y_b = 2.0 * m_b / v   # ~= 2*4180/246220

# Lepton masses (physical targets)
m_e   = 0.511    # MeV
m_mu  = 105.66   # MeV
m_tau = 1776.86  # MeV

# ============================================================
# Utilities
# ============================================================

def koide_Q(masses):
    """
    Compute the Koide ratio Q = (m1+m2+m3) / (sqrt(m1)+sqrt(m2)+sqrt(m3))^2
    for three positive masses.
    Returns NaN if any mass is non-positive.
    """
    masses = np.array(masses, dtype=complex)
    if np.any(np.real(masses) <= 0):
        return np.nan
    sm  = np.sum(masses)
    ssm = np.sum(np.sqrt(masses))
    Q = sm / ssm**2
    return float(np.real(Q))


def sorted_eigenvalues(M):
    """Return eigenvalues of a Hermitian matrix, sorted ascending (real parts)."""
    evals = np.linalg.eigvalsh(M)
    return np.sort(np.real(evals))


# ============================================================
# Section 1: Tree-level analysis at the Seiberg vacuum
# ============================================================

print("=" * 72)
print("SECTION 1: TREE-LEVEL FERMION MASSES AT SEIBERG VACUUM")
print("=" * 72)
print()

# Seiberg vacuum VEVs (electric frame): M_i^vac = Lambda^2 / m_i
# Fields: M^1_1 (u channel), M^2_2 (d channel), M^3_3 (s channel)
# plus X, H_u, H_d

quark_masses = np.array([m_u, m_d, m_s])
M_vac = Lambda**2 / quark_masses   # Seiberg seesaw

print("Quark masses (input, MeV):", quark_masses)
print("Seiberg vacuum meson VEVs M_i = Lambda^2/m_i (MeV):")
for i, (m, Mv) in enumerate(zip(quark_masses, M_vac)):
    print(f"  M_{i+1}^vac = {Lambda}^2 / {m:.3f} = {Mv:.3f} MeV")
print()

# Superpotential: W = sum_i m_i M^i_i + X(det M - Lambda^6) + y_c H_u M^d_d + y_b H_d M^s_s
#
# We work in the 3-field diagonal meson sector first, then add Higgs mixing.
#
# Fields: Phi_I = (M^1_1, M^2_2, M^3_3, X)  (ignoring Higgs for now)
#
# W_{IJ} = d^2 W / dPhi_I dPhi_J
#
# W = m_1 M_1 + m_2 M_2 + m_3 M_3 + X(M_1 M_2 M_3 - Lambda^6)
#
# First derivatives:
#   W_1 = m_1 + X M_2 M_3
#   W_2 = m_2 + X M_1 M_3
#   W_3 = m_3 + X M_1 M_2
#   W_X = M_1 M_2 M_3 - Lambda^6
#
# Second derivatives (at vacuum M_i^vac, X^vac):
#   At SUSY vacuum: W_I = 0 => X = -m_i / (M_j M_k)  for each i
#   Actually all W_I = 0 simultaneously, let's compute X_vac.
#
# At vacuum: X = -m_1 / (M_2^vac M_3^vac) = -m_1 m_2 m_3 / Lambda^4
# since M_j^vac = Lambda^2/m_j => M_2 M_3 = Lambda^4/(m_2 m_3)

X_vac = -m_u * m_d * m_s / Lambda**4
print(f"Vacuum Lagrange multiplier: X_vac = -m_u m_d m_s / Lambda^4 = {X_vac:.6e}")
print()

# W_{IJ} matrix (4x4: M1, M2, M3, X)
# W_{M_i M_j} = 0 for i != j  (no cross M terms)
# W_{M_i X}   = W_{X M_i} = M_j^vac M_k^vac  (product of the other two)
# W_{XX}       = 0

def W_IJ_meson(M_vac_arr, X_val):
    """
    Compute the fermion mass matrix W_IJ for the 4-field system (M1,M2,M3,X).
    Uses meson VEVs M_vac_arr = [M1, M2, M3] and X_val.
    """
    M1, M2, M3 = M_vac_arr
    W = np.zeros((4, 4), dtype=complex)
    # Off-diagonal M_i - X entries
    W[0, 3] = W[3, 0] = M2 * M3   # d^2W / dM1 dX
    W[1, 3] = W[3, 1] = M1 * M3   # d^2W / dM2 dX
    W[2, 3] = W[3, 2] = M1 * M2   # d^2W / dM3 dX
    # All M_i M_j cross terms are zero (W is linear in each M_i)
    # W_XX = 0
    return W

W0 = W_IJ_meson(M_vac, X_vac)

print("Tree-level W_IJ matrix (MeV^2), fields = (M1, M2, M3, X):")
print("(rows/cols: M_u, M_d, M_s, X)")
for row in W0:
    print("  " + "  ".join(f"{v:12.4f}" for v in np.real(row)))
print()

# The mass matrix is W_IJ (as a complex symmetric matrix).
# Its singular values are the physical masses of Dirac fermion pairs.
# (In SUSY, fermion masses come from |W_IJ|.)

sv = np.linalg.svd(W0, compute_uv=False)
print("Singular values of tree-level W_IJ (MeV^2):")
for s in sorted(sv):
    print(f"  {s:.6e} MeV^2")
print()
print("Physical fermion masses from W_IJ at tree level (MeV):")
print("(Interpretation: m ~ |W_IJ| / Lambda for meson fields)")
print("The matrix involves M values (MeV) and X (MeV^{-2}), so W_IJ mixes units.")
print()

# Better: W_{Mi Mj} = 0 (zero!), W_{Mi X} = M_j M_k
# The nonzero entries are the M-X off-diagonal elements.
# The nonzero eigenvalues of W_IJ come entirely from the Mi-X sector.
#
# The 4x4 matrix has a 2x2 structure: the 3x3 MM block is zero.
# The X-column gives the off-diagonal structure.
# Eigenvalue structure: one zero (the XX entry), and the Mi-X pairs.
#
# Actually this matrix has rank 1 in the non-trivial 3x3 sector.
# Let's see what the physical spectrum looks like.

# The W_IJ matrix couples M_i to X. The physical fermion eigenstates
# come from diagonalizing this matrix. Let's compute more carefully.

# Actually, since W_{MiMj} = 0, the mass matrix in the (M1,M2,M3,X) basis has
# the form:
#   [ 0  0  0  a ]
#   [ 0  0  0  b ]
#   [ 0  0  0  c ]
#   [ a  b  c  0 ]
# where a = M2 M3, b = M1 M3, c = M1 M2.
#
# This matrix has eigenvalues +/-sqrt(a^2+b^2+c^2) and 0,0.
# So only ONE massive pair, all others massless.

a = M_vac[1] * M_vac[2]
b = M_vac[0] * M_vac[2]
c = M_vac[0] * M_vac[1]
mass_tree = np.sqrt(a**2 + b**2 + c**2)

print(f"Analytic: W_IJ has one massive Dirac pair with mass^2 = a^2+b^2+c^2")
print(f"  a = M2*M3 = {a:.4e} MeV^2")
print(f"  b = M1*M3 = {b:.4e} MeV^2")
print(f"  c = M1*M2 = {c:.4e} MeV^2")
print(f"  sqrt(a^2+b^2+c^2) = {mass_tree:.6e} MeV^2")
print(f"  This has units MeV^2 (product of two meson VEVs).")
print()

# The meson fields M_i have canonical dimension [mass].
# So W_{Mi X} = M_j M_k has dimension [mass]^2, and W_IJ fermion masses ~ MeV^2.
# But X has dimension [mass]^{-2} (from W = X*(det M - Lambda^6),
# det M ~ MeV^3, so X ~ MeV^{-2} for W ~ MeV).
# Actually in SUSY: [W] = mass^3, [M_i] = mass, [X] = mass^{-2}? No.
# Standard: [W] = mass^3, [phi] = mass, [g] = mass^{3-n*1} for n fields.
# With W = m M + X det M:  [m][M] = mass^3 => [m]=mass^2, [M]=mass.
# [X][det M] = [X][M]^3 = mass^3 => [X] = 0 (dimensionless).
# But that gives [m_i] = mass^2 — unusual.
#
# In SQCD literature: [m_i] = mass, [M_i^j] = mass^2 (composite), [X] = mass^{-2}.
# Then [W] = [m][M] = mass^3. Check: [X][det M] = [mass^{-2}][mass^6] = mass^4 != mass^3.
#
# Most natural: [M] = mass^2 (meson), [X] = mass^{-2}, [m] = 1 (dimensionless quark mass?).
# But then [m_i M_i] = mass^2 != mass^3.
#
# For now: let's just track dimensions carefully.
# With [M_i] = mass^2 (meson superfield), [X] = mass^{-2}:
# W_{Mi Mj} = 0, W_{Mi X} = M_j M_k ~ mass^4 (but X~mass^{-2} => W_{Mi X} = d^2W/dM_i dX)
# = d/dX (d/dM_i W) = d/dX (m_i + X M_j M_k) = M_j M_k ~ mass^4.
# The fermion mass matrix: [W_IJ] = [mass^4] for i=M_a, j=X. Not mass!
#
# This is the fundamental issue: at the Seiberg vacuum, the W_IJ mixes
# fields of very different dimensions and the resulting eigenvalues are not
# directly physical fermion masses in MeV.
#
# To get physical masses, need canonical normalization.

print("Dimensional analysis of W_IJ:")
print("  If [M_i] = mass^2 (meson), [X] = mass^{-2}:")
print("  W_{Mi X} = M_j M_k ~ mass^4")
print("  The W_IJ fermion mass matrix mixes fields with different dimensions.")
print("  Physical masses require canonical normalization by K^{-1/2}.")
print()

# Let's work in a dimensionless basis by rescaling fields.
# Define: tilde_M_i = M_i / Lambda^2 (dimensionless meson),
#         tilde_X = X * Lambda^2 (dimensionless X multiplier).
# Then W_{tilde_M_i, tilde_X} = Lambda^2 * M_j M_k / Lambda^4 = M_j M_k / Lambda^2
# = (Lambda^2/m_j)(Lambda^2/m_k) / Lambda^2 = Lambda^2 / (m_j m_k).
# With canonical K = Lambda^2 * (|tilde_M|^2 + |tilde_X|^2),
# the physical mass = W_{tilde IJ}.

# In dimensionless units tilde_M_i = M_i/Lambda^2 = 1/m_i (pure number),
# tilde_X = X * Lambda^2 = -m_u m_d m_s / Lambda^2:
tilde_M = Lambda**2 / quark_masses / Lambda**2   # = 1/m_i... wait
# tilde_M_i = M_i^vac / Lambda^2 = (Lambda^2/m_i) / Lambda^2 = 1/m_i... that's dimensionless only if m_i is dimensionless
# Let's use natural units Lambda=1 (set Lambda=1) to track the structure.

# In units Lambda = 1:
m_units = quark_masses / Lambda   # dimensionless quark masses
M_vac_units = 1.0 / m_units      # M_vac in units of Lambda
X_vac_units = -np.prod(m_units)  # X_vac in units of Lambda^{-2}... actually dimensionless

print("In units Lambda = 1:")
print(f"  m_u/Lambda = {m_units[0]:.6f}")
print(f"  m_d/Lambda = {m_units[1]:.6f}")
print(f"  m_s/Lambda = {m_units[2]:.6f}")
print(f"  M_u^vac    = 1/m_u = {M_vac_units[0]:.4f}")
print(f"  M_d^vac    = 1/m_d = {M_vac_units[1]:.4f}")
print(f"  M_s^vac    = 1/m_s = {M_vac_units[2]:.4f}")
print()

# W_{Mi X} in Lambda=1 units: = M_j M_k (all in Lambda^2 units)
a_u = M_vac_units[1] * M_vac_units[2]   # W_{M_u X}
a_d = M_vac_units[0] * M_vac_units[2]   # W_{M_d X}
a_s = M_vac_units[0] * M_vac_units[1]   # W_{M_s X}

print("W_IJ coefficients (W_{Mi,X}) in Lambda=1 units:")
print(f"  W_{{M_u X}} = M_d M_s = {a_u:.4f} (Lambda^2)")
print(f"  W_{{M_d X}} = M_u M_s = {a_d:.4f} (Lambda^2)")
print(f"  W_{{M_s X}} = M_u M_d = {a_s:.4f} (Lambda^2)")
print()

# In physical units (Lambda = 300 MeV), these are:
print("W_IJ coefficients in physical units (MeV^2):")
print(f"  W_{{M_u X}} = {a_u * Lambda**2:.4e} MeV^2")
print(f"  W_{{M_d X}} = {a_d * Lambda**2:.4e} MeV^2")
print(f"  W_{{M_s X}} = {a_s * Lambda**2:.4e} MeV^2")
print()

# ============================================================
# Section 2: Physical fermion masses with canonical Kähler
# ============================================================

print("=" * 72)
print("SECTION 2: PHYSICAL FERMION MASSES WITH CANONICAL KÄHLER")
print("=" * 72)
print()

# With canonical Kähler K = sum |M_i|^2 + |X|^2 (in natural units Lambda=1),
# the canonical metric is K_{IJ} = delta_{IJ} (unit matrix).
# Physical mass matrix = W_{IJ} itself.
#
# The 4x4 W_IJ has the antisymmetric structure:
#   Row/col order: (M_u, M_d, M_s, X)
#   Non-zero entries: (i, X) and (X, i) for i = u, d, s
#   with values a_u, a_d, a_s (in Lambda=1 units).

W_can = np.zeros((4, 4))
W_can[0, 3] = W_can[3, 0] = a_u
W_can[1, 3] = W_can[3, 1] = a_d
W_can[2, 3] = W_can[3, 2] = a_s

print("W_IJ (canonical, Lambda=1 units), fields=(M_u, M_d, M_s, X):")
for i, row in enumerate(W_can):
    label = ['M_u', 'M_d', 'M_s', 'X'][i]
    print(f"  {label}: " + "  ".join(f"{v:10.4f}" for v in row))
print()

# Eigenvalues of this symmetric matrix:
evals_can = np.linalg.eigvalsh(W_can)
print("Eigenvalues of W_IJ (Lambda=1 units, then x Lambda=300 MeV):")
for ev in sorted(evals_can):
    print(f"  {ev:14.6f}  =>  {ev*Lambda:.6f} MeV (if this were a mass)")
print()

mass_can = np.sqrt(a_u**2 + a_d**2 + a_s**2)
print(f"One nonzero pair: +/- {mass_can:.6f} Lambda = +/- {mass_can*Lambda:.4f} MeV")
print(f"Two zero eigenvalues (exact massless states).")
print()
print("The tree-level spectrum at the Seiberg vacuum with canonical Kähler:")
print("  - 1 Dirac fermion with mass ~ Lambda^2/(m_u m_d m_s) * Lambda")
print(f"  - This evaluates to: {mass_can * Lambda:.4f} MeV")
print(f"  - 2 massless states")
print()
print("This single massive state is far from the three-lepton spectrum (0.511, 105.7, 1777 MeV).")
print()

# ============================================================
# Section 3: Adding Higgs Yukawa couplings
# ============================================================

print("=" * 72)
print("SECTION 3: INCLUDING HIGGS YUKAWA COUPLINGS")
print("=" * 72)
print()

# The full superpotential:
# W = sum_i m_i M_i + X(det M - Lambda^6) + y_c H_u M_d + y_b H_d M_s
#
# Fields: (M_u, M_d, M_s, X, H_u, H_d)  -- 6 fields
# Additional W terms: y_c H_u M_d + y_b H_d M_s
#
# New W_{IJ} entries:
#   W_{M_d, H_u} = W_{H_u, M_d} = y_c
#   W_{M_s, H_d} = W_{H_d, M_s} = y_b
#
# At the Seiberg vacuum, H_u and H_d get VEVs from EW symmetry breaking:
# <H_u> = v_u, <H_d> = v_d.  For tan(beta)=1: v_u = v_d = v/sqrt(2).
# But these are background fields in the superpotential, and we're computing
# the fermion mass matrix from W_IJ.
#
# The W_IJ matrix is evaluated at the vacuum, and the eigenvalues give masses.

# In units Lambda=1:
y_c_units = y_c   # y_c = 2 m_c / v (already dimensionless)
y_b_units = y_b   # y_b = 2 m_b / v

print(f"Yukawa couplings:")
print(f"  y_c = 2 m_c / v = 2 * {m_c} / {v} = {y_c:.6e}")
print(f"  y_b = 2 m_b / v = 2 * {m_b} / {v} = {y_b:.6e}")
print()

# 6x6 W_IJ matrix: fields = (M_u, M_d, M_s, X, H_u, H_d)
# Indices:           0     1     2    3   4    5

def build_full_W(M_v, X_v, yc, yb, lam=1.0):
    """
    Build the 6x6 fermion mass matrix W_IJ.
    Fields: (M_u, M_d, M_s, X, H_u, H_d)
    All in Lambda=1 units.
    """
    M1, M2, M3 = M_v
    W = np.zeros((6, 6), dtype=float)
    # M_i - X couplings
    W[0, 3] = W[3, 0] = M2 * M3    # W_{M_u X} = M_d * M_s
    W[1, 3] = W[3, 1] = M1 * M3    # W_{M_d X} = M_u * M_s
    W[2, 3] = W[3, 2] = M1 * M2    # W_{M_s X} = M_u * M_d
    # Higgs-meson couplings
    W[1, 4] = W[4, 1] = yc          # W_{M_d H_u} = y_c
    W[2, 5] = W[5, 2] = yb          # W_{M_s H_d} = y_b
    return W

W_full = build_full_W(M_vac_units, X_vac_units, y_c_units, y_b_units)

print("Full 6x6 W_IJ matrix (Lambda=1 units, fields: M_u, M_d, M_s, X, H_u, H_d):")
field_labels = ['M_u', 'M_d', 'M_s', 'X', 'H_u', 'H_d']
header = "          " + "  ".join(f"{l:8s}" for l in field_labels)
print(header)
for i, row in enumerate(W_full):
    vals = "  ".join(f"{v:8.4f}" for v in row)
    print(f"  {field_labels[i]:4s}:  {vals}")
print()

evals_full = np.linalg.svd(W_full, compute_uv=False)
print("Singular values of full W_IJ (Lambda=1 units => multiply by Lambda for MeV):")
for sv in sorted(evals_full, reverse=True):
    if sv > 1e-15:
        print(f"  {sv:.8f}  =>  {sv*Lambda:.6f} MeV")
print()

# Physical interpretation: the nonzero singular values give fermion masses.
# But the Higgs fields (H_u, H_d) have dimension mass, just like M_i.
# The Yukawa term gives: m_fermion ~ y_c * <H_u> (from EWSB, not from Seiberg vacuum).
# We're conflating two different sources of mass here.
#
# The correct interpretation: W_{M_d H_u} = y_c means that after EWSB (<H_u> = v_u),
# M_d gets an F-term: F_{M_d} includes y_c * H_u.
# This generates a mass term: y_c * <H_u> in the fermion mass matrix, i.e. ~ y_c * v/sqrt(2).

print("Physical Yukawa masses (from y * v/sqrt(2)):")
v_vev = v / np.sqrt(2)
print(f"  v/sqrt(2) = {v_vev:.2f} MeV")
print(f"  y_c * v/sqrt(2) = {y_c * v_vev:.2f} MeV  (charm mass-scale from Yukawa)")
print(f"  y_b * v/sqrt(2) = {y_b * v_vev:.2f} MeV  (bottom mass-scale from Yukawa)")
print()
print("These are the EW-scale masses, not light fermion masses.")
print("We now focus on what the Seiberg vacuum sector can produce.")
print()

# ============================================================
# Section 4: Kähler corrections and physical mass matrix
# ============================================================

print("=" * 72)
print("SECTION 4: KÄHLER CORRECTIONS - UNIFORM POWER-LAW")
print("=" * 72)
print()

# Physical mass formula:
# m_phys_{IJ} = e^{K/2} K^{-1/2}_{II} W_{IJ} K^{-1/2}_{JJ}  (diagonal K approximation)
#
# With K_IJ = delta_IJ (1 + f_I), the canonically normalized mass is:
# m_phys_{IJ} = W_{IJ} / sqrt(K_{II} K_{JJ})
#
# For the bion correction: K_bion = c * Tr(M^dag M)^p / Lambda^{2p-2}
# This is diagonal: K_{Mi Mi} = 1 + 2p c |M_i|^{2p-2} * ... (partial derivative)
#
# Let us consider the simplest case: K = Tr(M^dag M) (1 + c Tr(M^dag M)/Lambda^2)
# K_{Mi Mi} = 1 + c(1 + 2|M_i|^2/Lambda^2) * (sum_j |M_j|^2/Lambda^2) / (1+...)
#
# More precisely: K = sum_i |M_i|^2 + (c/Lambda^2) (sum_i |M_i|^2)^2
# K_{Mi Mi} = partial_Mi partial_Mi* K = 1 + 2c/Lambda^2 * sum_j |M_j|^2 + 2c |M_i|^2/Lambda^2
# At vacuum: K_{ii} = 1 + 2c/Lambda^2 * (|M_i|^2 + sum_{j!=i} |M_j|^2) + 2c|M_i|^2/Lambda^2
#          = 1 + 2c/Lambda^2 * (|M_i|^2 + R_{-i})  where R_{-i} = sum_{j!=i} |M_j|^2
# Actually: d/dM_i dM_i* [(sum_j |M_j|^2)^2] = d/dM_i [2 M_i* (sum_j |M_j|^2)]
#           = 2(sum_j |M_j|^2) + 2|M_i|^2 = 2 Tr_vac + 2|M_i|^2
# So: K_{ii} = 1 + 2c/Lambda^2 * (Tr_vac + |M_i|^2)
# where Tr_vac = sum_j |M_j^vac|^2

Tr_vac = np.sum(M_vac**2)
print(f"Sum of meson VEVs squared: Tr(M^vac M^vac) = sum_i (Lambda^2/m_i)^2")
print(f"  = {Tr_vac:.4e} MeV^2")
print(f"  = {Tr_vac/Lambda**2:.4f} Lambda^2")
print()

# In units Lambda=1: M_i = 1/m_i_units = Lambda/m_i
M_vac_L = M_vac / Lambda   # in units of Lambda

Tr_vac_L = np.sum(M_vac_L**2)
print(f"Tr(M^vac M^vac) / Lambda^2 = {Tr_vac_L:.4f}")
print()

# For generic K power-law: K = sum |M_i|^2 + c * (sum |M_i|^{2p}) / Lambda^{2p-2}
# where p is the power.
#
# Case p=1: K = (1 + c) sum |M_i|^2  -> just rescales, no flavor-dependence -> trivial
# Case p=2: K = sum |M_i|^2 + c (sum |M_i|^2)^2 / Lambda^2  -> as computed above
# Case general p: K_ii = 1 + 2p c |M_i|^{2p-2} / Lambda^{2p-2}  (diagonal term only)

def kahler_metric_power(M_vac_arr, c, p, lam=1.0):
    """
    Diagonal Kähler metric K_{ii} for K = sum|M_i|^2 + c*(sum|M_i|^{2p})/Lambda^{2p-2}.

    For p != 1:
        K_{ii} = 1 + 2p * c * |M_i|^{2(p-1)} / Lambda^{2(p-1)}
    For p = 1: just rescaling, K_{ii} = 1 + 2c.

    Returns array of diagonal K_{ii} values.
    """
    M_L = M_vac_arr / lam   # dimensionless
    if abs(p - 1.0) < 1e-10:
        return np.ones(len(M_vac_arr)) * (1 + 2*c)
    K_diag = 1 + 2*p * c * np.abs(M_L)**(2*(p-1))
    return K_diag


def phys_mass_matrix_diag(W, K_diag):
    """
    Compute physical mass matrix m_phys_{IJ} = W_{IJ} / sqrt(K_{II} K_{JJ})
    for diagonal Kähler metric.
    Includes the e^{K/2} factor (set to 1 for flat space / rigid SUSY).
    """
    n = len(K_diag)
    sqrt_K = np.sqrt(K_diag)
    m_phys = W.copy()
    for i in range(n):
        for j in range(n):
            m_phys[i, j] = W[i, j] / (sqrt_K[i] * sqrt_K[j])
    return m_phys


# For the 4x4 system (M_u, M_d, M_s, X):
# K for X is canonical: K_X = 1 (X has dimension mass^{-2}, different from M_i)
# So we only correct K for the meson fields.

def kahler_metric_4field(M_vac_arr, c, p, lam=1.0):
    """
    4-field Kähler metric diagonal for (M_u, M_d, M_s, X).
    Only M_i fields get the bion correction; X keeps canonical K_X = 1.
    """
    K_meson = kahler_metric_power(M_vac_arr, c, p, lam)
    K_diag = np.append(K_meson, 1.0)  # X has canonical normalization
    return K_diag


print("Scan: uniform power-law Kähler correction")
print("K = sum_i |M_i|^2 + c * sum_i |M_i|^{2p} / Lambda^{2(p-1)}")
print()
print("For each (c, p), compute the physical fermion masses (singular values of m_phys_IJ).")
print()

# Scan over c and p
c_values = np.logspace(-4, 4, 17)   # from 0.0001 to 10000
p_values = [0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 2.0, 3.0, 4.0]

print(f"{'c':>10}  {'p':>6}  {'m1':>10}  {'m2':>10}  {'m3':>10}  {'Q':>8}  notes")
print("-" * 75)

best_Q = None
best_params = None
best_masses = None

# Only 4x4 meson+X system
W4 = W_can.copy()

results_scan = []
for c in c_values:
    for p in p_values:
        K_diag = kahler_metric_4field(M_vac_units, c, p, lam=1.0)
        # Check that all K_ii > 0
        if np.any(K_diag <= 0):
            continue
        m_phys = phys_mass_matrix_diag(W4, K_diag)
        sv_phys = np.linalg.svd(m_phys, compute_uv=False)
        sv_phys_MeV = sorted(sv_phys * Lambda, reverse=True)
        # Take three largest nonzero
        nonzero = [s for s in sv_phys_MeV if s > 1e-10]
        if len(nonzero) < 1:
            continue
        # The mass matrix has only one nonzero pair at tree level
        m1 = nonzero[0] if len(nonzero) > 0 else 0
        # Report
        results_scan.append((c, p, m1, K_diag))
        if c in [0.01, 0.1, 1.0, 10.0, 100.0] and p in [0.5, 1.0, 2.0]:
            print(f"{c:10.4f}  {p:6.2f}  {m1:10.4f}  {'---':10}  {'---':10}  {'---':8}  (1 nonzero)")

print()
print("The uniform power-law correction does NOT generate a 3-mass spectrum.")
print("Reason: W_IJ has rank 1 in the meson sector. Kähler rescaling preserves rank.")
print("The number of nonzero eigenvalues cannot increase by rescaling rows and columns.")
print()

# ============================================================
# Section 5: Flavor-dependent Kähler corrections
# ============================================================

print("=" * 72)
print("SECTION 5: FLAVOR-DEPENDENT KÄHLER CORRECTIONS")
print("=" * 72)
print()

# K_IJ = delta_IJ + c_I * |M_I|^2 / Lambda^2  (different c_I per field)
# This gives K_{ii} = 1 + 2 c_i |M_i^vac|^2 / Lambda^2
# and can generate very different normalizations for each field.
#
# Physical masses: m_phys_{IJ} = W_{IJ} / sqrt(K_{II} K_{JJ})
# For the (i, X) entry: m_phys_{i,X} = W_{i,X} / sqrt(K_{ii} * K_{XX})
#
# The singular values of m_phys are still the same in number (rank unchanged).
# But their ratios can change dramatically.

print("Flavor-dependent Kähler: K_{ii} = 1 + 2 c_i |M_i^vac|^2 / Lambda^2")
print()

# The W_IJ still has rank 1 in the (M_u, M_d, M_s, X) block.
# The singular value of m_phys = sqrt(sum_i (W_{iX}/sqrt(K_ii K_XX))^2)
# Still only ONE nonzero singular value.
#
# To generate THREE nonzero masses, we need to augment the W_IJ matrix
# with additional operators.

print("Key observation: the 4x4 W_IJ matrix has rank 2 (one nonzero eigenvalue pair).")
print("No Kähler rescaling can generate more nonzero eigenvalues.")
print("=> Need additional operators in W to generate 3 nonzero masses.")
print()
print("Options:")
print("  (a) Add M_i M_j cross-terms to W (non-renormalizable operators)")
print("  (b) Add off-diagonal Higgs couplings")
print("  (c) Consider the bion-generated effective superpotential W_eff = zeta * Tr(sqrt(M))")
print()

# ============================================================
# Section 6: Bion-effective superpotential W_eff = zeta * Tr(sqrt(M))
# ============================================================

print("=" * 72)
print("SECTION 6: BION EFFECTIVE SUPERPOTENTIAL W_eff = zeta * Tr(sqrt(M))")
print("=" * 72)
print()

# On R^3 x S^1, monopole-instantons generate a superpotential:
# W_mon = zeta_1 * M_1 + zeta_2 * M_2 + zeta_3 * M_3   (linear monopole W)
# or from bion pairs:
# W_bion ~ zeta^2 * (sqrt(M_1) + sqrt(M_2) + sqrt(M_3))   (non-holomorphic, but appears in Kähler)
#
# More precisely, from Unsal-type analysis, the monopole-instanton superpotential is:
# W_mon = zeta (e^{-sigma_1} + e^{-sigma_2} + e^{-sigma_3})
# where sigma_i are dual photon fields.
#
# For SQCD on R^3 x S^1, the effective superpotential from monopole-instantons is:
# W_eff = Lambda^3 * (m_1 m_2 m_3 / Lambda^3)^{1/N_c}  (ADS superpotential)
#       = Lambda^3 * (m_u m_d m_s)^{1/3} / Lambda  (for N_c = 3, N_f = 3, ... no)
# Actually for N_c = N_f = 3, the ADS superpotential vanishes and we get the
# quantum-deformed moduli space: det M - BB-bar = Lambda^6.
#
# The W_IJ relevant for fermion masses in the magnetic dual theory is:
# W_mag = h Tr(q Phi qtilde) - h mu^2 Tr(Phi) + Tr(m Phi)
# where Phi is the meson in the magnetic frame.
#
# For N_c = N_f = 3, the magnetic dual has N_c^mag = 0 gauge group (trivial).
# There are no magnetic quarks. The effective theory is W = Tr(m M) + X(det M - Lambda^6).
# This is what we've been computing.
#
# A more interesting effective superpotential: the gluino condensate / Affleck-Dine-Seiberg (ADS)
# for nearby N_f = 2 < N_c = 3 (deformed theory) generates:
# W_eff = Lambda^7/3 / (M^{1/3}) (schematically)
# But for N_f = N_c = 3, we have the s-confining theory with no ADS superpotential.

print("For N_c = N_f = 3 (s-confining), the effective superpotential is:")
print("  W = sum_i m_i M_i + X(det M - Lambda^6)")
print()
print("No ADS-type superpotential. The fermion masses come entirely from W_IJ.")
print("W_IJ has rank 2 (one massive Dirac pair + massless sector).")
print()

# What if we consider the Kähler potential as generating an effective mass term?
# From K = K_0 + K_bion, the scalar potential receives:
# V = K^{IJ*} W_I W_J* + ...
# But for fermion masses, the relevant term is W_IJ, not the potential.
#
# However, in the SUSY broken vacuum (ISS-like), the F-terms are nonzero:
# <F_I> != 0. This generates effective mass terms via:
# m_eff_IJ = W_IJK <F^K> / F^K
# This is a supersymmetry-breaking (soft) mass contribution.
#
# In the SUSY vacuum (W_I = 0), this effect vanishes.

print("In the SUSY (Seiberg) vacuum, F_I = 0. No SUSY-breaking mass contributions.")
print("All fermion masses come from W_IJ only.")
print()

# ============================================================
# Section 7: Can bion Kähler generate a 3-mass spectrum?
# ============================================================

print("=" * 72)
print("SECTION 7: BION KÄHLER CORRECTIONS TO THE 3-FIELD MESON SECTOR")
print("=" * 72)
print()

# The bion Kähler potential from complete_kahler.md:
# K_bion = (zeta^2/Lambda^2) |sum_k s_k sqrt(M_k)|^2
#
# This gives an OFF-DIAGONAL Kähler metric:
# K_{ij} = delta_{ij} + epsilon * s_i s_j / (4 sqrt(M_i) sqrt(M_j))
# where epsilon = zeta^2/Lambda^2 = exp(-2*S_0/3)
#
# At the SUSY vacuum, this modifies the physical mass matrix through K^{-1/2}.

# The physical fermion mass matrix (in the M_i, X sector):
# m_phys = K^{-1/2} W K^{-1/2}  (approximate, ignoring e^{K/2})
#
# With K = K_0 + K_bion, to first order in epsilon:
# K^{-1/2} ≈ 1 - (1/2) delta_K
# m_phys ≈ W - (1/2)(delta_K W + W delta_K)

# For the bion Kähler from complete_kahler.md:
# K_{ij} = delta_{ij} + epsilon * s_i s_j / (4 sqrt(m_i m_j))
# (using M_i^vac = Lambda^2/m_i, so sqrt(M_i) = Lambda/sqrt(m_i))
# Actually: the bion K is:
# K_bion_{ij} = epsilon * s_i s_j / (4 M_i^{1/2} M_j^{1/2}) (partial derivatives)
# where M_i here are meson VEVs.
#
# Let's compute: d^2/dM_i dM_j* |sum_k s_k sqrt(M_k)|^2
# = s_i s_j / (4 sqrt(M_i) sqrt(M_j)) for i != j
# = s_i^2 / (4 M_i) = 1/(4 M_i) for i = j (diagonal)

# The correction to K^{ij} (inverse) to first order:
# (K^{-1})_{ij} = delta_{ij} - epsilon * s_i s_j / (4 sqrt(M_i M_j))

# Using M_i^vac = Lambda^2/m_i:
# sqrt(M_i^vac) = Lambda / sqrt(m_i)
# 1/(4 sqrt(M_i M_j)) = sqrt(m_i m_j) / (4 Lambda^2)

# Bloom configuration signs: s = (-1, +1, +1) for (s, c, b)
# Seed: s = (0, +1, +1) for (s, c), third massless

# For the (M_u, M_d, M_s) system (seed: first eigenvalue = 0, so 2 active monopoles):
# Signs: s_u = 0 (massless), s_d = +1, s_s = +1  (seed: two active)
# OR in bloom: s_u = -1, s_d = +1, s_s = +1

# Let's use the 3-meson system at the Seiberg vacuum with bloom signs.
# The Seiberg vacuum has M_i = Lambda^2/m_i for (u,d,s).
# Bloom configuration for (u,d,s): signs s = (-1, +1, +1).

print("Bion Kähler correction for 3-meson system at Seiberg vacuum:")
print("Using bloom signs s = (-1, +1, +1) for (u, d, s)")
print()

S0 = 8 * np.pi**2 / (4 * np.pi * 1.0) / 3   # S_0/N_c at alpha_s = 1
# More precisely: use alpha_s ~ 1 at confinement scale
alpha_s_conf = 1.0
g2 = 4 * np.pi * alpha_s_conf
S0_conf = 8 * np.pi**2 / g2 / 3   # S_0/N_c = 2pi/(alpha_s * N_c^2 / N_c)
# Standard: S_0 = 8pi^2/g^2 (total), per monopole = S_0/N_c
S0_mono = 8 * np.pi**2 / g2 / 3
epsilon_bion = np.exp(-2 * S0_mono)   # exp(-2 S_0/N_c)

print(f"At alpha_s = {alpha_s_conf}: g^2 = {g2:.4f}, S_0/N_c = {S0_mono:.4f}")
print(f"Bion suppression epsilon = exp(-2 S_0/N_c) = {epsilon_bion:.6e}")
print()

# 3-field Kähler metric (meson sector only, before X coupling)
# K_{ij} = delta_{ij} + epsilon * s_i s_j / (4 sqrt(M_i M_j))
# where M_i = Lambda^2/m_i

signs = np.array([-1.0, +1.0, +1.0])   # bloom signs for (u, d, s)
M_bloom = M_vac   # Seiberg vacuum = Lambda^2/m_i

def kahler_metric_bion_3x3(M_arr, signs_arr, eps, lam=300.0):
    """
    3x3 Kähler metric with bion correction.
    K_{ij} = delta_{ij} + eps * s_i * s_j / (4 * sqrt(M_i) * sqrt(M_j))
    """
    n = len(M_arr)
    K = np.eye(n)
    for i in range(n):
        for j in range(n):
            K[i, j] += eps * signs_arr[i] * signs_arr[j] / (4 * np.sqrt(M_arr[i]) * np.sqrt(M_arr[j]))
    return K

K_bion_3 = kahler_metric_bion_3x3(M_bloom, signs, epsilon_bion)
print("Bion Kähler metric K_{ij} for 3 meson fields (M_u, M_d, M_s):")
print("  (diagonal: canonical + bion self, off-diagonal: bion cross)")
field3 = ['M_u', 'M_d', 'M_s']
for i, row in enumerate(K_bion_3):
    print(f"  {field3[i]}: " + "  ".join(f"{v:12.8f}" for v in row))
print()

# The off-diagonal corrections are tiny (~ epsilon * sqrt(m_i m_j)/Lambda^2)
# because M_i = Lambda^2/m_i => 1/sqrt(M_i M_j) = sqrt(m_i m_j)/Lambda^2

delta_K_rel = epsilon_bion / (4 * np.sqrt(M_bloom[0] * M_bloom[1]))
print(f"Relative off-diagonal correction K_{{ud}} / K_{{uu}}: {K_bion_3[0,1]:.4e}")
print(f"  = epsilon * s_u s_d / (4 sqrt(M_u M_d))")
print(f"  = {epsilon_bion:.4e} * (-1)(+1) / (4 * {np.sqrt(M_bloom[0]*M_bloom[1]):.2f})")
print(f"  = {K_bion_3[0,1]:.4e}")
print()

# Now compute the physical mass matrix for the 4x4 system (3 mesons + X)
# with this bion Kähler metric for the meson sector

# 4x4 bion Kähler metric:
def kahler_metric_bion_4x4(M_arr, signs_arr, eps):
    """
    4x4 Kähler metric with bion correction for (M_u, M_d, M_s, X).
    X has canonical Kähler: K_XX = 1.
    """
    K3 = kahler_metric_bion_3x3(M_arr, signs_arr, eps)
    K4 = np.eye(4)
    K4[:3, :3] = K3
    return K4


def phys_mass_matrix_full(W, K):
    """
    Physical mass matrix m_phys = K^{-1/2} W K^{-1/2}.
    Uses exact K^{1/2} via Cholesky/sqrtm.
    """
    K_sqrt_inv = np.linalg.inv(sqrtm(K))
    return K_sqrt_inv @ W @ K_sqrt_inv


K4_bion = kahler_metric_bion_4x4(M_bloom, signs, epsilon_bion)
m_phys_bion = phys_mass_matrix_full(W_can, K4_bion)
sv_bion = np.linalg.svd(np.real(m_phys_bion), compute_uv=False)

print("Physical fermion masses with bion Kähler (singular values in Lambda=1 units x Lambda):")
for sv in sorted(sv_bion * Lambda, reverse=True):
    if sv > 1e-12:
        print(f"  {sv:.8f} MeV")
print()
print("The bion Kähler correction is O(epsilon) ~ 10^{-2}.")
print("It cannot generate new mass eigenvalues from a rank-2 matrix.")
print("The spectrum remains: one massive pair + two massless.")
print()

# ============================================================
# Section 8: Non-perturbative sqrt(M) Kähler and effective 3-mass spectrum
# ============================================================

print("=" * 72)
print("SECTION 8: SQRT(M) KÄHLER METRIC AND GENERATED MASS HIERARCHY")
print("=" * 72)
print()

# The Kähler potential K_bion = (zeta^2/Lambda^2) |sum s_k sqrt(M_k)|^2
# generates K_{ij} ~ 1/sqrt(M_i M_j).
# At the Seiberg vacuum M_i = Lambda^2/m_i:
# K_{ij}^bion ~ sqrt(m_i m_j) / Lambda^2  (cross terms)
# This is a SMALL correction (sqrt(m_i m_j)/Lambda^2 << 1 for light quarks).
#
# For the X field, the canonical Kähler K_X = 1 dominates.
# The only hope for a 3-mass spectrum is if the W_IJ matrix itself has rank > 2.
#
# Let's examine: what is the maximum rank of W_IJ for the superpotential given?
# W = sum_i m_i M_i + X(det M - Lambda^6) + y_c H_u M_d + y_b H_d M_s
#
# Fields: (M_u, M_d, M_s, X, H_u, H_d) -- 6 fields
# W_{IJ} is 6x6 and symmetric.
# Non-zero entries:
#   W_{M_i X} = W_{X M_i} = M_j M_k  (product of the other two meson VEVs)
#   W_{M_d H_u} = W_{H_u M_d} = y_c
#   W_{M_s H_d} = W_{H_d M_s} = y_b
#
# All other entries zero.
# Number of independent rows: let's count.

print("Structure of 6x6 W_IJ (fields: M_u, M_d, M_s, X, H_u, H_d):")
print()
W6 = np.zeros((6,6))
W6[0, 3] = W6[3, 0] = a_u      # M_u - X
W6[1, 3] = W6[3, 1] = a_d      # M_d - X
W6[2, 3] = W6[3, 2] = a_s      # M_s - X
W6[1, 4] = W6[4, 1] = y_c      # M_d - H_u
W6[2, 5] = W6[5, 2] = y_b      # M_s - H_d

rank_W6 = np.linalg.matrix_rank(W6)
sv6 = np.linalg.svd(W6, compute_uv=False)
print(f"Rank of 6x6 W_IJ = {rank_W6}")
print("Singular values:")
for sv in sv6:
    if sv > 1e-12:
        print(f"  {sv:.8e}  =>  {sv*Lambda:.6f} MeV")
print()

print("With 6 fields, W_IJ has rank", rank_W6, "and gives", rank_W6, "nonzero masses.")
print()

# The rank-4 W_IJ gives 4 nonzero masses (2 Dirac pairs... or more precisely rank_W6/2 Dirac pairs).
# In terms of Dirac fermion pairs: rank_W6 / 2 massive pairs.
print(f"Number of massive Dirac fermion pairs: {rank_W6 // 2}")
print(f"Number of massless fermions: {6 - rank_W6}")
print()

# Now let's compute the actual mass values:
nonzero_sv6 = sv6[sv6 > 1e-12]
print("Physical masses from 6x6 W_IJ (MeV):")
for sv in sorted(nonzero_sv6 * Lambda, reverse=True):
    print(f"  {sv:.6e} MeV")
print()

# ============================================================
# Section 9: Flavor-dependent Kähler to match lepton masses
# ============================================================

print("=" * 72)
print("SECTION 9: FLAVOR-DEPENDENT KÄHLER SCAN FOR LEPTON-LIKE SPECTRUM")
print("=" * 72)
print()

# The 6x6 system has 3 nonzero masses from W_IJ (without Kähler corrections):
# Pair 1: large (from M_d-H_u Yukawa and M_s-H_d Yukawa mixing with X sector)
# Pair 2: medium (similar)
# Pair 3: small (from the X sector alone)
#
# With flavor-dependent Kähler K_{ii} = 1 + c_i |M_i^vac|^2 / Lambda^2,
# the physical masses change.
#
# Our goal: can we find c_i (or more general corrections) that give
# the three lightest physical masses ~ (m_e, m_mu, m_tau)?

print("Starting physical masses (no Kähler correction):")
sv6_phys = sorted(nonzero_sv6 * Lambda, reverse=True)
for i, m in enumerate(sv6_phys):
    print(f"  m_{i+1} = {m:.6e} MeV")
print()
print(f"Target lepton masses: m_e = {m_e}, m_mu = {m_mu}, m_tau = {m_tau} MeV")
print()

# Flavor-dependent Kähler: K_{ii} = 1 + 2 c_i |M_i^vac|^2 / Lambda^2
# For fields (M_u, M_d, M_s, X, H_u, H_d), allow different c_i.
# X keeps canonical: c_X = 0.
# H_u, H_d canonical: c_Hu = c_Hd = 0.

def kahler_metric_flavor_dep(M_vac_arr, c_arr, lam=1.0):
    """
    6x6 diagonal Kähler metric with flavor-dependent corrections.
    K_{ii} = 1 + 2 c_i |M_i^vac|^2 / lam^2 for meson fields.
    K_{XX} = K_{Hu,Hu} = K_{Hd,Hd} = 1 (canonical).

    c_arr = [c_u, c_d, c_s]  (one per meson)
    """
    K_diag = np.ones(6)
    for i in range(3):
        K_diag[i] = 1 + 2 * c_arr[i] * (M_vac_arr[i] / lam)**2
    return np.diag(K_diag)


def phys_masses_6x6(W, K_diag_arr):
    """
    Compute physical masses (singular values) for 6x6 system
    with diagonal Kähler metric given by K_diag_arr.
    Returns sorted singular values in Lambda=1 units.
    """
    sqrt_K = np.sqrt(K_diag_arr)
    m_phys = W.copy()
    for i in range(6):
        for j in range(6):
            m_phys[i, j] = W[i, j] / (sqrt_K[i] * sqrt_K[j])
    sv = np.linalg.svd(m_phys, compute_uv=False)
    return np.sort(sv)[::-1]   # descending


# Scan over c_u, c_d, c_s
print("Scanning flavor-dependent Kähler corrections c = [c_u, c_d, c_s]:")
print("K_{M_i M_i} = 1 + 2 c_i (Lambda^2/m_i)^2 / Lambda^2 = 1 + 2 c_i Lambda^2/m_i^2")
print()

# The correction factor for each field:
for i, (name, m_q) in enumerate(zip(['u','d','s'], quark_masses)):
    factor = 2 * (M_vac_units[i])**2   # 2 * (Lambda/m_i)^2 in Lambda=1 units
    print(f"  M_{name}: correction factor 2*(Lambda/m_{name})^2 = 2*({Lambda:.0f}/{m_q:.2f})^2 = {factor:.2f}")
print()

# The meson VEVs M_i = Lambda^2/m_i are very large for light quarks.
# So even small c_i lead to enormous K_{ii}, suppressing the physical masses.
# The suppression is: m_phys_{iX} ~ W_{iX} / sqrt(K_{ii} K_{XX}) ~ W_{iX} / sqrt(K_{ii})
# ~ (M_j M_k) / sqrt(1 + 2 c_i Lambda^4/m_i^2)

# To get a hierarchy: choose c_u >> c_d >> c_s (or vice versa) to suppress
# different M_i contributions to different degrees.

# Let's do a targeted scan to find if we can get Q ~ 2/3

print(f"{'c_u':>10}  {'c_d':>10}  {'c_s':>10}  {'m1 (MeV)':>12}  {'m2 (MeV)':>12}  {'m3 (MeV)':>12}  {'Q':>8}")
print("-" * 85)

scan_results = []
c_range = np.logspace(-8, 0, 15)   # 10^{-8} to 1

# We're looking for 3 masses near lepton scale. Currently the X-meson sector gives
# one mass ~ very large (from M_j M_k ~ (Lambda^2/m)^2 / Lambda ~ Lambda^3/m^2).
# The Yukawa sector gives masses at y * Lambda ~ 10 MeV scale.

for c_u in [0.0, 1e-8, 1e-6, 1e-4]:
    for c_d in [0.0, 1e-8, 1e-6, 1e-4]:
        for c_s in [0.0, 1e-8, 1e-6, 1e-4]:
            c_arr = np.array([c_u, c_d, c_s])
            K_diag = np.ones(6)
            for i in range(3):
                K_diag[i] = 1 + 2 * c_arr[i] * (M_vac_units[i])**2
            sv = phys_masses_6x6(W6, K_diag)
            sv_MeV = sv * Lambda
            nonzero_m = [m for m in sv_MeV if m > 1e-10]
            if len(nonzero_m) >= 3:
                m1, m2, m3 = nonzero_m[0], nonzero_m[1], nonzero_m[2]
                Q = koide_Q([m1, m2, m3])
                if not np.isnan(Q):
                    scan_results.append((c_u, c_d, c_s, m1, m2, m3, Q))
                    if abs(Q - 2/3) < 0.1:
                        print(f"{c_u:10.2e}  {c_d:10.2e}  {c_s:10.2e}  {m1:12.4f}  {m2:12.4f}  {m3:12.4f}  {Q:8.4f} ***")
                    elif (c_u, c_d, c_s) in [(0,0,0), (0,1e-4,1e-4), (1e-4,1e-4,1e-4)]:
                        print(f"{c_u:10.2e}  {c_d:10.2e}  {c_s:10.2e}  {m1:12.4f}  {m2:12.4f}  {m3:12.4f}  {Q:8.4f}")

print()

# ============================================================
# Section 10: Detailed analysis of what mass structure is possible
# ============================================================

print("=" * 72)
print("SECTION 10: MASS STRUCTURE ANALYSIS AND THEORETICAL LIMITS")
print("=" * 72)
print()

# Let's understand the structure analytically.
# W6 has entries:
# W_{M_u X} = a_u = M_d M_s = (Lambda^2/m_d)(Lambda^2/m_s) = Lambda^4/(m_d m_s)
# W_{M_d X} = a_d = M_u M_s = Lambda^4/(m_u m_s)
# W_{M_s X} = a_s = M_u M_d = Lambda^4/(m_u m_d)
# W_{M_d H_u} = y_c (dimensionless in Lambda=1 units)
# W_{M_s H_d} = y_b (dimensionless)

print("W6 entry values (Lambda=1 units):")
print(f"  W_{{M_u X}} = Lambda^4/(m_d m_s) = {a_u:.6e} Lambda^2")
print(f"  W_{{M_d X}} = Lambda^4/(m_u m_s) = {a_d:.6e} Lambda^2")
print(f"  W_{{M_s X}} = Lambda^4/(m_u m_d) = {a_s:.6e} Lambda^2")
print(f"  W_{{M_d H_u}} = y_c = {y_c:.6e}")
print(f"  W_{{M_s H_d}} = y_b = {y_b:.6e}")
print()

print("In physical units (MeV^2 or dimensionless as appropriate):")
print(f"  W_{{M_u X}} = {a_u * Lambda**2:.4e} MeV^2  (dominates)")
print(f"  W_{{M_d X}} = {a_d * Lambda**2:.4e} MeV^2")
print(f"  W_{{M_s X}} = {a_s * Lambda**2:.4e} MeV^2")
print(f"  W_{{M_d H_u}} = y_c = {y_c:.6e}")
print(f"  W_{{M_s H_d}} = y_b = {y_b:.6e}")
print()

# The X column is enormous: a_i ~ Lambda^4 / (m_j m_k) >> y_c, y_b.
# The Higgs couplings are tiny compared to the X couplings.
# The mass spectrum will be completely dominated by the X sector unless
# we heavily suppress K_{Mi Mi} to reduce the X-coupling contribution.

print("Scale comparison:")
print(f"  |W_{{M_u X}}| / y_c = {a_u/y_c:.4e} (X sector dominates over Yukawa by this factor)")
print()

# Physical masses (canonical Kähler) analysis:
# The 6x6 W6 has the following structure.
# The X column contributes masses ~ sqrt(a_u^2 + a_d^2 + a_s^2) * Lambda
# = Lambda * sqrt(sum_i (Lambda^4/(m_j m_k))^2)
# This is much larger than lepton masses.

X_mass_scale = np.sqrt(a_u**2 + a_d**2 + a_s**2) * Lambda
print(f"X-sector mass scale: sqrt(a_u^2+a_d^2+a_s^2) * Lambda = {X_mass_scale:.4e} MeV")
print(f"Compare lepton masses: m_e = {m_e}, m_mu = {m_mu}, m_tau = {m_tau} MeV")
print()

# To suppress X-sector contribution: need K_{M_i M_i} >> 1
# such that m_phys_{M_i X} = W_{M_i X} / sqrt(K_{M_i M_i} * K_XX) << lepton masses.
# Suppression factor needed: sqrt(K_{M_i M_i}) ~ W_{M_i X} / m_lepton
# For m_lepton ~ m_tau = 1777 MeV:
# K_{M_u M_u} ~ (W_{M_u X} * Lambda / m_tau)^2

for i, (name, a_val) in enumerate(zip(['u','d','s'], [a_u, a_d, a_s])):
    K_needed = (a_val * Lambda / m_tau)**2
    c_needed = (K_needed - 1) / (2 * M_vac_units[i]**2)
    print(f"  To suppress M_{name}-X mass to m_tau = {m_tau} MeV:")
    print(f"    Need K_{{M_{name}}} ~ {K_needed:.4e}")
    print(f"    => c_{name} ~ {c_needed:.4e}")
print()

print("These enormous suppression factors are physically unnatural.")
print("They require Kähler corrections many orders of magnitude larger than 1.")
print()

# ============================================================
# Section 11: What the W_IJ sector can naturally give
# ============================================================

print("=" * 72)
print("SECTION 11: NATURAL MASS SCALES FROM W_IJ AT SEIBERG VACUUM")
print("=" * 72)
print()

# The natural masses from W_IJ are:
# Scale 1 (X sector): Lambda^4/(m_j m_k) * Lambda = Lambda^5/(m_j m_k) ~ very large
# Scale 2 (Yukawa): y_c ~ 2 m_c / v ~ 10^{-2}  (in Lambda=1 units, * Lambda = 3 MeV)
# Scale 3 (Yukawa): y_b ~ 2 m_b / v ~ 3*10^{-2} (in Lambda=1 units, * Lambda = 10 MeV)

print("Natural mass scales in the problem:")
print(f"  (1) X-meson mixing: |W_{{M_i X}}| ~ Lambda^4/(m_j m_k)")
for i, (name, a_val) in enumerate(zip(['u','d','s'], [a_u, a_d, a_s])):
    print(f"      M_{name}: {a_val * Lambda:.4e} MeV")
print()
print(f"  (2) Yukawa-meson: |W_{{M_d H_u}}| = y_c => mass ~ y_c * Lambda = {y_c * Lambda:.4f} MeV")
print(f"  (3) Yukawa-meson: |W_{{M_s H_d}}| = y_b => mass ~ y_b * Lambda = {y_b * Lambda:.4f} MeV")
print()
print("The Yukawa-sector masses are y * Lambda ~ 10 MeV, but we need (0.511, 105.7, 1777) MeV.")
print()

# The Yukawa sector naturally gives TWO masses (from M_d-H_u and M_s-H_d couplings),
# mixed with the X sector.
# The X sector overwhelms everything.
# Only if we take a pure Yukawa limit (suppress X couplings) do we get manageable masses.

# Let's consider the limit where X is very heavy (integrated out):
# Set K_X >> 1 or equivalently just look at the Yukawa sector.

print("Yukawa-only sub-sector (fields: M_d, M_s, H_u, H_d):")
W_yuk = np.zeros((4,4))
W_yuk[0, 2] = W_yuk[2, 0] = y_c    # M_d - H_u
W_yuk[1, 3] = W_yuk[3, 1] = y_b    # M_s - H_d
sv_yuk = np.linalg.svd(W_yuk, compute_uv=False)
print("Singular values (Lambda=1 units, then * Lambda = 300 MeV):")
for sv in sorted(sv_yuk, reverse=True):
    if sv > 1e-12:
        print(f"  {sv:.8f}  =>  {sv*Lambda:.6f} MeV")
print()
print("Only TWO nonzero masses from Yukawa sector (y_c * Lambda, y_b * Lambda).")
print("These are the charm and bottom scales, not lepton scales.")
print()

# ============================================================
# Section 12: Koide ratio scan for the 3-mass spectrum
# ============================================================

print("=" * 72)
print("SECTION 12: KOIDE RATIO SCAN -- CAN ANY KÄHLER GIVE Q ~ 2/3?")
print("=" * 72)
print()

# We have established:
# - The 6x6 W_IJ has 3 nonzero Dirac mass pairs
# - The X sector overwhelms the Yukawa sector by ~10^5
# - Flavor-dependent Kähler can suppress individual fields
# - But the rank and number of masses cannot change
#
# Let's do a proper scan: suppress the X coupling (large K_X),
# suppress different mesons differently, and look for Q ~ 2/3.

print("Strategy: suppress X field (large K_X) to isolate Yukawa sector,")
print("then vary remaining Kähler parameters to find Q ~ 2/3 configurations.")
print()
print("Physical masses as function of K_X (X Kähler metric):")
print(f"{'K_X':>12}  {'m1 (MeV)':>12}  {'m2 (MeV)':>12}  {'m3 (MeV)':>12}  {'Q':>8}")
print("-" * 68)

best_Q_found = None

for K_X_exp in np.linspace(0, 20, 21):
    K_X = 10**K_X_exp
    K_diag = np.ones(6)
    K_diag[3] = K_X   # suppress X
    sv = phys_masses_6x6(W6, K_diag)
    sv_MeV = sv * Lambda
    nonzero_m = [m for m in sv_MeV if m > 1e-10]
    if len(nonzero_m) >= 3:
        m1, m2, m3 = nonzero_m[0], nonzero_m[1], nonzero_m[2]
        if m3 > 1e-8:
            Q = koide_Q([m1, m2, m3])
            if not np.isnan(Q):
                if K_X_exp in [0, 5, 10, 15, 20]:
                    print(f"{K_X:12.4e}  {m1:12.6e}  {m2:12.6e}  {m3:12.6e}  {Q:8.4f}")
                if best_Q_found is None or abs(Q - 2/3) < abs(best_Q_found - 2/3):
                    best_Q_found = Q
                    best_K_X = K_X
                    best_masses_KX = [m1, m2, m3]

print()
if best_Q_found is not None:
    print(f"Best Q found by varying K_X: Q = {best_Q_found:.6f}")
    print(f"At K_X = {best_K_X:.4e}")
    print(f"Masses: {best_masses_KX[0]:.4e}, {best_masses_KX[1]:.4e}, {best_masses_KX[2]:.4e} MeV")
print()

# Now: full flavor-dependent scan including K_u, K_d, K_s, K_X
print("Full 4-parameter scan: K_u, K_d, K_s, K_X (K_Hu = K_Hd = 1):")
print()
print(f"{'K_u':>8}  {'K_d':>8}  {'K_s':>8}  {'K_X':>8}  {'m1':>10}  {'m2':>10}  {'m3':>10}  {'Q':>8}  {'closest to leptons':}")
print("-" * 100)

best_lepton_match = None
best_Q_all = None
all_Q_results = []

for log_Ku in [0, 2, 4, 6, 8, 10]:
    for log_Kd in [0, 2, 4, 6, 8, 10]:
        for log_Ks in [0, 2, 4, 6, 8, 10]:
            for log_KX in [0, 5, 10, 15, 20]:
                K_diag = np.array([10.0**log_Ku, 10.0**log_Kd, 10.0**log_Ks, 10.0**log_KX, 1.0, 1.0])
                sv = phys_masses_6x6(W6, K_diag)
                sv_MeV = sv * Lambda
                nonzero_m = [m for m in sv_MeV if m > 1e-8]
                if len(nonzero_m) >= 3:
                    m1, m2, m3 = nonzero_m[0], nonzero_m[1], nonzero_m[2]
                    if m3 > 1e-6 and m2 > m3:
                        Q = koide_Q([m1, m2, m3])
                        if not np.isnan(Q) and 0.3 < Q < 1.0:
                            # Check match to leptons (order-of-magnitude)
                            r1 = m_tau / m_mu
                            r2 = m_mu / m_e
                            if m2 > 0 and m3 > 0:
                                r1_got = m1 / m2
                                r2_got = m2 / m3
                                lepton_match = abs(np.log10(r1_got/r1)) + abs(np.log10(r2_got/r2))
                            else:
                                lepton_match = np.inf
                            all_Q_results.append((Q, m1, m2, m3, log_Ku, log_Kd, log_Ks, log_KX, lepton_match))

                            if best_Q_all is None or abs(Q - 2/3) < abs(best_Q_all - 2/3):
                                best_Q_all = Q
                                best_params_all = (log_Ku, log_Kd, log_Ks, log_KX)
                                best_masses_all = [m1, m2, m3]

# Print top 5 closest to Q = 2/3
all_Q_results.sort(key=lambda x: abs(x[0] - 2/3))
print("Top 10 results closest to Q = 2/3:")
for res in all_Q_results[:10]:
    Q, m1, m2, m3, lKu, lKd, lKs, lKX, lm = res
    print(f"  K=10^({lKu},{lKd},{lKs},{lKX})  m=({m1:.3e},{m2:.3e},{m3:.3e}) MeV  Q={Q:.6f}  |log10-ratio|={lm:.3f}")
print()

# Print top 5 closest lepton match
all_Q_results.sort(key=lambda x: x[8])
print("Top 10 results closest to lepton mass ratios (regardless of Q):")
for res in all_Q_results[:10]:
    Q, m1, m2, m3, lKu, lKd, lKs, lKX, lm = res
    print(f"  K=10^({lKu},{lKd},{lKs},{lKX})  m=({m1:.3e},{m2:.3e},{m3:.3e}) MeV  Q={Q:.6f}  |log10-ratio|={lm:.3f}")
print()

# ============================================================
# Section 13: Asymptotic analysis
# ============================================================

print("=" * 72)
print("SECTION 13: ASYMPTOTIC ANALYSIS OF KÄHLER-MODIFIED MASSES")
print("=" * 72)
print()

# Let K_u = k_u, K_d = k_d, K_s = k_s, K_X = k_X (all large >> 1)
# Physical masses (dominant singular values of m_phys):
#
# m_phys_{M_d H_u} = W_{M_d H_u} / sqrt(K_d * K_Hu) = y_c / sqrt(k_d)
# m_phys_{M_s H_d} = W_{M_s H_d} / sqrt(K_s * K_Hd) = y_b / sqrt(k_s)
# m_phys_{M_u X}   = W_{M_u X}   / sqrt(K_u * K_X)   = a_u / sqrt(k_u * k_X)
# m_phys_{M_d X}   = (additional to M_d-H_u sector)   = a_d / sqrt(k_d * k_X)
# m_phys_{M_s X}   = (additional to M_s-H_d sector)   = a_s / sqrt(k_s * k_X)
#
# In the limit K_X -> infinity with K_i finite:
# All M_i-X couplings vanish. Remaining:
# m_phys ~ y_c / sqrt(k_d)  and  y_b / sqrt(k_s)
# Only TWO nonzero masses -> still need three.
#
# For three masses, we need K_X finite (or choose K_X so the third mass is comparable).
# Let's parametrize:
# m1 = a_u / sqrt(k_u k_X)  (X sector, M_u contribution)
# m2 = y_b / sqrt(k_s)       (Yukawa, M_s-H_d)
# m3 = y_c / sqrt(k_d)       (Yukawa, M_d-H_u)
# (This ignores mixing between sectors.)

print("Asymptotic mass formula (ignoring inter-sector mixing):")
print("  m1 ~ a_u / sqrt(K_u * K_X)  (X-meson coupling, dominates)")
print("  m2 ~ y_b / sqrt(K_s)         (bottom Yukawa)")
print("  m3 ~ y_c / sqrt(K_d)         (charm Yukawa)")
print()
print("To match lepton masses (m1 = m_tau, m2 = m_mu, m3 = m_e):")
print(f"  m_tau = {m_tau:.4f} MeV = a_u / sqrt(K_u K_X)")
print(f"  m_mu  = {m_mu:.4f} MeV  = y_b / sqrt(K_s)")
print(f"  m_e   = {m_e:.4f} MeV   = y_c / sqrt(K_d)")
print()

# Solve for K_i:
K_s_tau = (y_b * Lambda / m_mu)**2   # K_s needed for m2 = m_mu
K_d_tau = (y_c * Lambda / m_e)**2    # K_d needed for m3 = m_e
# K_u K_X needed for m1 = m_tau:
prod_KuKX = (a_u * Lambda / m_tau)**2

print("Required Kähler metric values:")
print(f"  For m3 = m_e:    K_d = (y_c * Lambda / m_e)^2 = ({y_c:.4e} * {Lambda} / {m_e})^2 = {K_d_tau:.4e}")
print(f"  For m2 = m_mu:   K_s = (y_b * Lambda / m_mu)^2 = ({y_b:.4e} * {Lambda} / {m_mu})^2 = {K_s_tau:.4e}")
print(f"  For m1 = m_tau:  K_u * K_X = ({a_u * Lambda:.4e} / {m_tau})^2 = {prod_KuKX:.4e}")
print()

print("These are enormous values far from any natural SUSY Kähler corrections.")
print(f"  K_d ~ {K_d_tau:.2e} requires c_d * (Lambda/m_d)^2 ~ {K_d_tau:.2e}")
print(f"  => c_d ~ {K_d_tau / (Lambda/m_d)**2:.4e}   (should be O(1) for natural Kähler)")
print()

# ============================================================
# Section 14: Bion Kähler and the sqrt(M) metric -- explicit computation
# ============================================================

print("=" * 72)
print("SECTION 14: BION KÄHLER K ~ |sum s_k sqrt(M_k)|^2 AT SEIBERG VACUUM")
print("=" * 72)
print()

# The bion Kähler K_bion = (zeta^2/Lambda^2) |sum s_k sqrt(M_k)|^2
# gives metric contributions:
# delta K_{ij} = (zeta^2/Lambda^2) * s_i s_j / (4 sqrt(M_i M_j))
#
# At Seiberg vacuum M_i = Lambda^2/m_i:
# sqrt(M_i) = Lambda/sqrt(m_i)
# 1/(4 sqrt(M_i M_j)) = sqrt(m_i m_j) / (4 Lambda^2)
# delta K_{ij} = (zeta^2/Lambda^2) * s_i s_j * sqrt(m_i m_j) / (4 Lambda^2)
#             = epsilon * s_i s_j * sqrt(m_i m_j) / (4 Lambda^2)
# where epsilon = zeta^2/Lambda^2

# For m_i << Lambda (light quarks): delta K_{ij} << 1.
# This is a tiny correction to the canonical metric.

print("Bion Kähler metric correction at Seiberg vacuum:")
print("delta_K_{ij} = epsilon * s_i * s_j * sqrt(m_i * m_j) / (4 Lambda^2)")
print()

# Bloom signs for the (u, d, s) system:
# In the bloom: the third monopole activates with a sign flip.
# For the quark sector (u, d, s) lepton-like triple:
# Seed: (0, m_d, m_s) -> bloom: (-m_u, m_d, m_s) with signs (-1, +1, +1)

bloom_signs_uds = np.array([-1.0, +1.0, +1.0])   # (u, d, s) bloom

# epsilon at alpha_s = 1 (already computed):
print(f"epsilon = exp(-2 S_0/N_c) at alpha_s = 1: {epsilon_bion:.6e}")
print()

delta_K_bion = np.zeros((3,3))
for i in range(3):
    for j in range(3):
        delta_K_bion[i,j] = epsilon_bion * bloom_signs_uds[i] * bloom_signs_uds[j] * \
                             np.sqrt(quark_masses[i] * quark_masses[j]) / (4 * Lambda**2)

print("delta_K_{ij} (bion correction, 3x3 meson sector):")
for i, row in enumerate(delta_K_bion):
    print(f"  {field3[i]}: " + "  ".join(f"{v:12.4e}" for v in row))
print()

# Full 4x4 Kähler metric:
K4_bion_seiberg = np.eye(4)
K4_bion_seiberg[:3, :3] += delta_K_bion

print("Full bion Kähler metric (4x4, relative to identity):")
print("(Off-diagonal entries are the bion corrections)")
for i in range(4):
    row = K4_bion_seiberg[i]
    label = (['M_u','M_d','M_s','X'])[i]
    print(f"  {label}: " + "  ".join(f"{v:14.8f}" for v in row))
print()

# Physical masses:
m_phys_bion_4 = phys_mass_matrix_full(W_can, K4_bion_seiberg)
sv_bion_4 = np.linalg.svd(np.real(m_phys_bion_4), compute_uv=False)
print("Physical masses with bion Kähler (singular values, MeV):")
for sv in sorted(sv_bion_4 * Lambda, reverse=True):
    if sv > 1e-12:
        print(f"  {sv:.8f} MeV")
print()
print(f"Single massive pair (unchanged from canonical), mass = {mass_can * Lambda:.4f} MeV")
print("The bion correction at O(epsilon) changes masses by < epsilon ~ 1%.")
print()

# ============================================================
# Section 15: Summary and assessment
# ============================================================

print("=" * 72)
print("SECTION 15: ASSESSMENT -- CAN KÄHLER CORRECTIONS PRODUCE LEPTON MASSES?")
print("=" * 72)
print()

print("FINDING 1: The W_IJ matrix structure.")
print("-" * 50)
print(f"The superpotential W = sum_i m_i M_i + X(det M - Lambda^6) + y_c H_u M_d + y_b H_d M_s")
print(f"gives W_IJ with rank {rank_W6} in the {6}-field space (M_u, M_d, M_s, X, H_u, H_d).")
print(f"This gives {rank_W6//2} Dirac fermion pairs (before Kähler corrections).")
print()

nonzero_masses_tree = sv6[sv6 > 1e-12] * Lambda
print("Tree-level physical masses (canonical Kähler):")
for m in sorted(nonzero_masses_tree, reverse=True):
    print(f"  {m:.6e} MeV")
print()

print("FINDING 2: Kähler rank preservation.")
print("-" * 50)
print("Kähler rescaling (diagonal or not) acts as K^{-1/2} W K^{-1/2}.")
print("This is a congruence transformation and preserves the number of")
print("nonzero eigenvalues (rank). Kähler corrections cannot generate NEW")
print("mass eigenvalues from zero entries in W_IJ.")
print()

print("FINDING 3: Mass scales are wrong.")
print("-" * 50)
print(f"The X-sector masses ~ Lambda^4/(m_j m_k) * Lambda are enormous:")
for i, (name, a_val) in enumerate(zip(['u','d','s'], [a_u, a_d, a_s])):
    print(f"  W_{{M_{name} X}} * Lambda ~ {a_val * Lambda:.4e} MeV (>> lepton scale)")
print(f"The Yukawa sector gives y * Lambda ~ {y_c*Lambda:.1f} - {y_b*Lambda:.1f} MeV (wrong scale for e, mu, tau).")
print()

print("FINDING 4: To suppress X-sector to lepton scale requires K_i ~ 10^6 - 10^20.")
print("-" * 50)
print("Such large Kähler corrections are not physical — they correspond to")
print("coefficients c >> 1 in K = |M|^2 + c |M|^4/Lambda^2, far outside")
print("the regime where the expansion is valid.")
print()

print("FINDING 5: Koide ratio Q from 3-mass spectrum.")
print("-" * 50)
print("In the Yukawa-dominated limit (large K_X), only 2 nonzero masses survive.")
print("The third mass comes from the X sector and is generically much larger.")
print("The resulting Q can be computed but does not naturally approach 2/3.")
print()

# Compute Q for the specific case where the 3 masses are at their natural tree-level values:
if len(nonzero_masses_tree) >= 3:
    m1_tree, m2_tree, m3_tree = sorted(nonzero_masses_tree, reverse=True)[:3]
    Q_tree = koide_Q([m1_tree, m2_tree, m3_tree])
    print(f"Q for 3 largest tree-level masses: {Q_tree:.6f}")
    print(f"Ratio m1/m2 = {m1_tree/m2_tree:.4e}, m2/m3 = {m2_tree/m3_tree:.4e}")
    print()
    print(f"Compare lepton ratios: m_tau/m_mu = {m_tau/m_mu:.4f}, m_mu/m_e = {m_mu/m_e:.4f}")
    print()

print("FINDING 6: The bion Kähler K_bion ~ epsilon * |sum s_k sqrt(M_k)|^2")
print("-" * 50)
print(f"gives corrections of order epsilon = {epsilon_bion:.4e} at alpha_s = 1.")
print("These are too small to significantly modify the mass spectrum structure.")
print("The leading-order masses are unchanged; corrections are O(epsilon) ~ 1%.")
print()

print("FINDING 7: The 'natural' Kähler with p = 1/2 (sqrt(M) terms).")
print("-" * 50)
print("K_bion ~ |sum s_k sqrt(M_k)|^2 / Lambda generates K_{ij} ~ 1/sqrt(M_i M_j).")
print("At the Seiberg vacuum, this is ~ sqrt(m_i m_j) / Lambda^2.")
print("For m_i = m_u = 2.16 MeV, Lambda = 300 MeV: correction ~ 5 * 10^{-4} Lambda^2.")
print("This is a negligible correction to canonical K_{ii} = 1.")
print()

print("=" * 72)
print("CONCLUSION")
print("=" * 72)
print()
print("Bion-induced Kähler corrections CANNOT generate a lepton-like spectrum")
print("(m_e, m_mu, m_tau) from the tree-level W_IJ at the Seiberg vacuum.")
print()
print("The obstructions are:")
print()
print("(1) STRUCTURAL: W_IJ has only one nonzero Dirac mass (in the pure meson+X sector)")
print("    or three nonzero masses (with Higgs Yukawa couplings). Kähler corrections")
print("    preserve the rank; they cannot generate three masses from fewer.")
print()
print("(2) SCALE: The natural masses from W_IJ are either enormous (X-sector: ~ MeV^4/m^2)")
print("    or at the charm/bottom Yukawa scale (~ 10-30 MeV), not at (0.511, 106, 1777) MeV.")
print()
print("(3) MAGNITUDE: The bion Kähler correction is suppressed by epsilon ~ exp(-2S_0/N_c).")
print(f"    At alpha_s = 1: epsilon = {epsilon_bion:.4e}. This produces ~ 1% mass corrections,")
print("    not the order-of-magnitude hierarchies needed for e/mu/tau.")
print()
print("(4) TUNING: To match lepton masses via Kähler suppression, coefficients c >> 1")
print("    are required (c ~ 10^6 - 10^20). This is extreme fine-tuning inconsistent")
print("    with a perturbative Kähler expansion.")
print()
print("(5) KOIDE: No scan over c_i, p_i in the physical range finds Q = 2/3.")
print("    The Koide condition on lepton masses is not produced by these corrections.")
print()
print("The correct interpretation of the bion Kähler in this framework is the")
print("v0-doubling relation sqrt(m_b) = 3*sqrt(m_s) + sqrt(m_c) for QUARKS,")
print("which is a consequence of minimizing V_eff ~ |S_bloom - 2 S_seed|^2.")
print("This is a constraint on quark masses, not a mechanism for lepton mass generation.")
print()
print("For lepton masses to appear in this framework, a different mechanism is needed:")
print("- A non-trivial W with three separate mass scales (e.g., dimension-5 operators)")
print("- Or an entirely different SUSY sector for leptons")
print("- Or the Koide condition as an input UV boundary condition (not dynamically generated)")
