#!/usr/bin/env python3
"""
ISS (Intriligator-Seiberg-Shih) Coleman-Weinberg potential computation.

SQCD with N_c = 3, N_f = 4.
Magnetic dual: quarks q_i, q~^i, meson Phi^i_j (4x4 matrix).
Superpotential: W = h Tr(q Phi q~) - h mu^2 Tr(Phi) + Tr(m Phi)

Parameters:
  h = sqrt(3) (O'Raifeartaigh-Koide value)
  m = diag(m_u, m_d, m_s, m_c) quark masses
  mu^2 = SUSY-breaking scale (try mu^2 = f_pi^2)

Key question: Does the CW potential preserve the Koide condition
from UV quark masses to IR meson masses?
"""

import numpy as np
from numpy import sqrt, log, pi, diag
from itertools import product as iprod

# ============================================================
# Physical parameters
# ============================================================
N_c = 3
N_f = 4
N_c_mag = N_f - N_c  # = 1, magnetic colors

# Quark masses (MeV)
m_u = 2.16
m_d = 4.67
m_s = 93.4
m_c = 1270.0
m_quarks = np.array([m_u, m_d, m_s, m_c])

# Magnetic Yukawa coupling
h = np.sqrt(3.0)

# SUSY-breaking scale
mu2 = 92.0**2  # f_pi^2 in MeV^2
mu = 92.0  # MeV

# Cutoff
Lambda = 300.0  # MeV

print("=" * 80)
print("ISS MODEL: COLEMAN-WEINBERG POTENTIAL AND KOIDE TRANSMISSION")
print("=" * 80)
print(f"\nParameters:")
print(f"  N_c = {N_c}, N_f = {N_f}, N_c_mag = N_c_mag = {N_c_mag}")
print(f"  h = sqrt(3) = {h:.6f}")
print(f"  mu^2 = f_pi^2 = {mu2:.1f} MeV^2  (mu = {mu:.1f} MeV)")
print(f"  Lambda = {Lambda:.1f} MeV")
print(f"  Quark masses: m_u={m_u}, m_d={m_d}, m_s={m_s}, m_c={m_c} MeV")

# ============================================================
# 1. TREE-LEVEL SPECTRUM
# ============================================================
print("\n" + "=" * 80)
print("1. TREE-LEVEL SPECTRUM AROUND ISS METASTABLE VACUUM")
print("=" * 80)

print("""
The ISS metastable vacuum arises from the rank condition.
The meson field Phi is N_f x N_f = 4x4.
The F-term condition F_Phi = h q q~ - h mu^2 I + m = 0
cannot be solved: q q~ has rank at most N_c_mag = 1, but
h mu^2 I - m has rank 4 (generically).

Decompose Phi into blocks:
  Phi = [[Phi_0 (1x1),  Z (1x3)],
         [Z~ (3x1),  chi (3x3)]]

where the upper-left 1x1 block is the "unbroken" direction
(where q q~ can cancel the F-term) and the lower-right 3x3
is the "broken" block (rank condition prevents F=0).

At the ISS vacuum:
  <q> = <q~> = (sqrt(mu^2 - m_1/h), 0, 0, 0) for the lightest flavor
  <Phi_0> = 0 (pseudo-modulus)
  <chi> = 0

But more precisely, in the ISS model with distinct masses,
the vacuum structure depends on which flavor direction
breaks SUSY least badly. The rank-1 q q~ aligns with
the direction that minimizes the vacuum energy.
""")

# The ISS vacuum: q q~ cancels F-term in the lightest-mass direction
# For simplicity and generality, let's label the flavors 1,...,N_f
# with m_1 <= m_2 <= ... <= m_N_f.
# The q q~ VEV fills the first N_c_mag = 1 direction.

# Sort masses
idx_sorted = np.argsort(m_quarks)
m_sorted = m_quarks[idx_sorted]
labels = ['u', 'd', 's', 'c']
labels_sorted = [labels[i] for i in idx_sorted]

print(f"\nFlavors sorted by mass: {labels_sorted}")
print(f"Masses sorted: {m_sorted}")

# The vacuum: <q q~> nonzero in the lightest direction (flavor u)
# <q_1 q~^1> = mu^2 - m_1/h = mu_eff^2
mu_eff2 = mu2 - m_sorted[0] / h
print(f"\nVacuum alignment: <q q~> in flavor '{labels_sorted[0]}' direction")
print(f"  mu_eff^2 = mu^2 - m_1/h = {mu2:.1f} - {m_sorted[0]/h:.4f} = {mu_eff2:.4f} MeV^2")

# Pseudo-modulus: Phi_0 (the (1,1) component of Phi)
# At tree level, Phi_0 is a flat direction.

print("""
Mass spectrum at the ISS vacuum:

The key result (from ISS 2006, hep-th/0602239):

For each flavor i in the "broken" block (i = 2,...,N_f),
the tree-level masses are:

  Bosonic masses for Re(chi_ij), Im(chi_ij):
    m^2_B,pm = h^2 |Phi_0|^2 +/- |h^2 mu^2 - h m_i| (for diagonal components)

  Actually, the standard ISS spectrum for the broken sector mesons
  (the chi fields, i.e., the 3x3 block) with Phi_0 = 0 is:

  Boson masses (2 real d.o.f. each):
    m^2_{Re chi_ii} = h^2 mu_eff^2 - (h mu^2 - m_i)  [WRONG sign would be tachyonic]

Let me be more careful. The ISS spectrum at <Phi_0> = 0:
""")

# CAREFUL DERIVATION of ISS spectrum
# Following ISS (Intriligator, Seiberg, Shih) hep-th/0602239
# and Shih hep-th/0703196

# The superpotential in the magnetic dual:
# W = h Tr(q Phi q~) - h mu_i^2 Phi^i_i + m_i Phi^i_i
# where mu_i^2 = mu^2 for all i (universal), and m_i are quark masses.

# Define: mu_i^2 = mu^2 - m_i/h  (effective SUSY-breaking scale per flavor)
mu_i2 = mu2 - m_sorted / h
print("Effective SUSY-breaking per flavor:")
for i in range(N_f):
    print(f"  mu_{labels_sorted[i]}^2 = mu^2 - m_{labels_sorted[i]}/h = {mu_i2[i]:.4f} MeV^2")

# Decompose:
# q = (rho, 0, 0, 0) + fluctuations  (N_c_mag = 1, so rho is a single number)
# <rho> = sqrt(mu_1^2 - m_1/h) ≈ sqrt(mu^2)  (for mu >> m/h)

# Actually in ISS, the vacuum has:
# <q_alpha^i> = delta_{alpha,1} delta^{i,1} * sqrt(h mu^2 - m_1) / sqrt(h)
#             ≈ mu (for m_1 << h mu^2)
# where alpha is magnetic color (here alpha = 1 only) and i is flavor.

# Let me use the standard notation from Shih 0703196.
# The vacuum:
# <q> = (phi_0, 0, ..., 0),  <q~> = (phi_0, 0, ..., 0)
# where phi_0^2 = mu^2 - m_1/(h)   [assuming m_1 is the smallest mass]
# <Phi_0> = arbitrary (pseudo-modulus, tree-level flat direction)
# <chi> = 0 (3x3 block of Phi)

phi_0_sq = mu2 - m_sorted[0] / h
phi_0 = sqrt(phi_0_sq)
print(f"\n<phi_0>^2 = {phi_0_sq:.4f} MeV^2")
print(f"<phi_0> = {phi_0:.4f} MeV")

# At the vacuum with <Phi_0> = X (the pseudo-modulus VEV):
# The mass spectrum splits into sectors:

# SECTOR 1: The "unbroken" (1x1) block
# This contains the pseudo-modulus X = Phi^1_1 (tree-level flat)
# and the Goldstino (from SUSY breaking).

# SECTOR 2: The "broken" (3x3) block chi = Phi^a_b (a,b = 2,3,4)
# These get tree-level masses from the F-term SUSY breaking.

# SECTOR 3: The "mixing" block Z, Z~ (off-diagonal 1x3 and 3x1)

# =============================================
# For each broken-sector flavor a = 2,3,4 (i.e., d, s, c after sorting):
# The relevant mass terms come from:
# W ⊃ h phi_0 q_a chi^a_1 + h phi_0 q~^1 chi^1_a + ...
#     + (m_a - h mu^2) chi^a_a

# The SUSY-breaking F-terms are:
# F_Phi^a_a = -h mu^2 + m_a   (for a in broken block, a=2,3,4)
# These are nonzero: |F_a| = h mu^2 - m_a  (for m_a < h mu^2)

F_a = h * mu2 - m_sorted  # F-terms for each flavor
print("\nF-terms (|F_a| = h*mu^2 - m_a):")
for i in range(N_f):
    print(f"  F_{labels_sorted[i]} = {F_a[i]:.4f} MeV^2")

# Note: F_1 ≈ 0 (canceled by q q~ VEV), but F_a for a>1 is nonzero.

# =============================================
# TREE-LEVEL MASS MATRICES
# =============================================
print("\n" + "-" * 60)
print("Tree-level mass matrices")
print("-" * 60)

# Following the ISS analysis more carefully:
# At the metastable vacuum with Phi_0 = X (pseudo-modulus):

# For the broken-sector diagonal mesons chi^a_a (a = 2,...,N_f):
# The scalar mass-squared matrix (in the Re/Im basis) is:
#   M^2_scalar = [[|hX|^2 + |F_a/h|, 0], [0, |hX|^2 - |F_a/h|]]
# Wait, this isn't right either. Let me be more precise.

# The tree-level potential from F-terms for the broken-sector diagonal mesons:
# V ⊃ |F_chi^a_a|^2 = |h q_alpha q~^alpha - h mu^2 + m_a|^2
# At the vacuum, q_alpha q~^alpha = phi_0^2 delta_{a,1}, so for a > 1:
# V ⊃ |m_a - h mu^2 + h X delta_a_1_for_chi_sector|^2
# For chi^a_a with a in the broken block:
# V ⊃ |m_a - h mu^2|^2 = F_a^2    [constant, doesn't give mass to chi]

# The mass for chi comes from the YUKAWA terms:
# W ⊃ h q_alpha Phi^i_j q~_beta^j ...
# Actually in the magnetic dual, the magnetic quarks q have N_c_mag = 1
# color and N_f = 4 flavor. So:
# W = h q^i Phi_i^j q~_j - h mu^2 Tr(Phi) + Tr(m Phi)
# where q^i is a single-component (no color index) field with flavor i.

# F-terms:
# F_{Phi_i^j} = h q^i q~_j - h mu^2 delta_i^j + m_i delta_i^j
# F_{q^i} = h Phi_i^j q~_j
# F_{q~_j} = h q^i Phi_i^j

# At vacuum: q^1 = q~_1 = phi_0, q^a = q~_a = 0 (a=2,3,4), Phi = X e_11

# Expand around vacuum:
# q^i = (phi_0 + delta_rho, delta_q^2, delta_q^3, delta_q^4)
# q~_j = (phi_0 + delta_rho~, delta_q~_2, delta_q~_3, delta_q~_4)
# Phi_i^j = X e_11 + delta_Phi_i^j

# The mass terms come from second derivatives of |F|^2:
# V = sum |F_A|^2 where A runs over all fields.

# =============================================
# Mass matrix for the BROKEN-sector mesons chi^a_b (a,b = 2,3,4)
# =============================================

# From F_{q^a} = h Phi_a^j q~_j:
# At the vacuum, F_{q^a} = h Phi_a^1 phi_0 (leading contribution)
# So |F_{q^a}|^2 = h^2 phi_0^2 |Phi_a^1|^2 = h^2 phi_0^2 |Z_a|^2
# This gives mass^2 = h^2 phi_0^2 to the Z_a fields.

# For the diagonal chi^a_a components:
# The relevant F-terms are:
# F_{Phi_a^a} = h q^a q~_a - h mu^2 + m_a
# F_{Phi_a^b} = h q^a q~_b  (off-diagonal, a != b)

# Expanding around vacuum (q^a = q~_a = 0 for a>1):
# The mass for chi^a_a comes from cross-terms between
# F_{Phi_a^a} and the fluctuations:
# delta(F_{Phi_a^a}) = h delta_q^a delta_q~_a  (quadratic, gives quartic not mass)

# Actually, the mass for chi^a_a comes from the SUPERPOTENTIAL mass terms:
# W ⊃ (m_a - h mu^2) Phi_a^a = (m_a - h mu^2) chi^a_a

# This gives a SUPERSYMMETRIC mass to chi^a_a:
# m_{chi_a}^2 = 0  (no bosonic mass from this alone in SUSY limit)
# But SUSY is BROKEN, so there's a splitting.

# The CORRECT ISS spectrum (from ISS 2006 + Shih 2007):
# For each broken flavor direction a = 2,...,N_f:

# Define:
#   f_a = h mu^2 - m_a   (the F-term in that direction)
#   M_a = h X             (supersymmetric mass from pseudo-modulus coupling)

# Wait: there's no tree-level mass for chi^a_a from the pseudo-modulus X.
# X = Phi_1^1, and chi^a_a = Phi_a^a. These are DIFFERENT components.
# The coupling h q^i Phi_i^j q~_j with i=j=a gives h Phi_a^a q^a q~_a,
# but at the vacuum q^a = q~_a = 0 so there's no linear mixing.

# Let me restart with a cleaner approach. The ISS model with N_c_mag = 1
# is essentially a multi-field O'Raifeartaigh model.

print("""
ISS with N_c_mag = 1 as generalized O'Raifeartaigh model:

Chiral fields: X (pseudo-modulus = Phi^1_1),
               chi_a (broken mesons = Phi^a_a, a=2,3,4),
               Y_a = Phi^1_a (1x3 block),
               Y~_a = Phi^a_1 (3x1 block),
               rho = q^1 (aligned magnetic quark),
               rho~ = q~_1,
               psi_a = q^a (a=2,3,4),
               psi~_a = q~_a

Superpotential around vacuum:
  W = h rho X rho~
    + h rho Y_a psi~_a + h psi_a Y~_a rho~
    + h psi_a chi_a psi~_a
    + (m_1 - h mu^2) X
    + (m_a - h mu^2) chi_a
    + h phi_0 Y_a psi~_a + h phi_0 psi_a Y~_a   [from VEV]
    + h phi_0 X rho~ + h rho X phi_0  [from VEV]
    + ...

Actually, let me just write it as:
  W = h (phi_0 + rho)(X + ...) (phi_0 + rho~) + ...

This is getting complicated. Let me use the MASS MATRIX directly.
""")

# =============================================
# DIRECT MASS MATRIX COMPUTATION
# =============================================
# The fermionic mass matrix comes from W_ij = d^2 W / d phi_i d phi_j
# evaluated at the vacuum.

# Fields: X, chi_a (a=2,3,4), Y_a, Y~_a, rho, rho~, psi_a, psi~_a
# That's 1 + 3 + 3 + 3 + 1 + 1 + 3 + 3 = 18 complex fields

# But for the mass spectrum, we can organize by sectors that don't mix
# at tree level:

# SECTOR A: (rho, rho~, X) - the "unbroken" sector
# Contains the Goldstino and the pseudo-modulus

# SECTOR B: (psi_a, psi~_a, Y_a, Y~_a) for each flavor a - the "link" sector
# Each flavor is a 4-field subsystem

# SECTOR C: chi_a - the "broken" mesons
# These are SINGLETS that get only the SUSY-breaking mass

# Let me compute each sector.

print("\n" + "-" * 60)
print("SECTOR A: Pseudo-modulus + Goldstino (rho, rho~, X)")
print("-" * 60)

# W_A = h (phi_0 + rho) X (phi_0 + rho~) + (m_1 - h mu^2) X
#      = h phi_0^2 X + h phi_0 rho X + h phi_0 X rho~ + h rho X rho~ + (m_1 - h mu^2) X

# But phi_0^2 = mu^2 - m_1/h, so:
# h phi_0^2 + m_1 - h mu^2 = h(mu^2 - m_1/h) + m_1 - h mu^2 = 0
# The linear term in X vanishes — good, X is a flat direction.

# The fermionic mass matrix (W_ij) in basis (rho, rho~, X):
# W_{rho,rho~} = h X|_vac = h * <X>     [but <X> = 0 at tree level origin]
# W_{rho,X} = h (phi_0 + rho~)|_vac = h phi_0
# W_{rho~,X} = h (phi_0 + rho)|_vac = h phi_0
# W_{X,X} = 0, W_{rho,rho} = 0, W_{rho~,rho~} = 0

# With <X> = 0:
W_ferm_A = np.array([
    [0, 0, h * phi_0],
    [0, 0, h * phi_0],
    [h * phi_0, h * phi_0, 0]
])

print(f"\nFermionic mass matrix (sector A):")
print(f"  W_ij = [[0, 0, h*phi_0],")
print(f"          [0, 0, h*phi_0],")
print(f"          [h*phi_0, h*phi_0, 0]]")
print(f"  h*phi_0 = {h * phi_0:.4f} MeV")

# Eigenvalues: det(W - lambda I) = 0
# lambda(lambda^2 - 2 h^2 phi_0^2) = 0
# lambda = 0 (Goldstino), +/- sqrt(2) h phi_0

eig_A = np.linalg.eigvalsh(W_ferm_A)
print(f"  Eigenvalues: {np.sort(np.abs(eig_A))}")
print(f"  Expected: 0, sqrt(2)*h*phi_0 = {sqrt(2)*h*phi_0:.4f} MeV (x2)")

m_goldstino = 0.0
m_pseudo_ferm = sqrt(2) * h * phi_0
print(f"  Goldstino mass: {m_goldstino}")
print(f"  Pseudo-modulus fermion mass: {m_pseudo_ferm:.4f} MeV")

# F-term for X: F_X = h phi_0^2 + m_1 - h mu^2 = 0  (vanishes at vacuum)
# F-term for rho: F_rho = h X phi_0 + h Y_a psi~_a = 0 at vacuum
# F-term for rho~: similar

# The BOSONIC mass matrix in sector A:
# V = |F_X|^2 + |F_rho|^2 + |F_rho~|^2
# F_X = h(phi_0 + rho)(phi_0 + rho~) + m_1 - h mu^2
#      = h phi_0^2 + h phi_0(rho + rho~) + h rho rho~ + m_1 - h mu^2
#      = h phi_0(rho + rho~) + h rho rho~
# So dF_X/drho = h phi_0, dF_X/drho~ = h phi_0, dF_X/dX = 0

# F_rho = h X(phi_0 + rho~) + h sum_a Y_a psi~_a
# dF_rho/dX = h phi_0, dF_rho/drho~ = h X|_vac = 0

# F_rho~ = h(phi_0 + rho) X + h sum_a psi_a Y~_a
# dF_rho~/dX = h phi_0, dF_rho~/drho = h X|_vac = 0

# The bosonic mass-squared matrix (in the field basis) is:
# M^2_B = sum_A (dF_A/dphi_i)* (dF_A/dphi_j)
#       = (W_ij)^dagger (W_ij)   ... for the supersymmetric part

# Plus the SUSY-breaking contribution:
# M^2_break = -F_A d^2 F_A*/dphi_i dphi_j + h.c.

# For the pseudo-modulus X, the only nonzero F at the vacuum is from
# the broken chi_a sector. F_{chi_a} = m_a - h mu^2.
# But chi_a is in sector C, not A. So for sector A:
# The SUSY-breaking contribution comes only from F_X cross terms.
# Since F_X = 0 at the vacuum, there is NO tree-level SUSY-breaking
# mass splitting in sector A.

# So sector A masses:
# Scalars: same as fermions = {0, sqrt(2) h phi_0, sqrt(2) h phi_0}
# But X is FLAT (mass = 0 for both scalar and fermion of X).
# The CW potential will lift X.

print("\nSector A summary:")
print(f"  Goldstino (fermion): mass = 0")
print(f"  X (pseudo-modulus scalar): tree-level mass = 0 (flat direction)")
print(f"  Heavy fermion pair: mass = sqrt(2)*h*phi_0 = {m_pseudo_ferm:.4f} MeV")
print(f"  Corresponding heavy scalar pair: mass = same (SUSY not broken in this sector)")

# =============================================
print("\n" + "-" * 60)
print("SECTOR B: Link fields (psi_a, psi~_a, Y_a, Y~_a) for each flavor a")
print("-" * 60)

# For each broken flavor a = 2, 3, 4:
# Relevant superpotential terms:
# W_B ⊃ h psi_a chi_a psi~_a + h phi_0 Y_a psi~_a + h psi_a Y~_a phi_0
#       + h rho Y_a psi~_a + h psi_a Y~_a rho~
#       + higher order

# At the vacuum, the mass terms for (psi_a, Y~_a) and (psi~_a, Y_a):
# W ⊃ h phi_0 (Y_a psi~_a + psi_a Y~_a) + (m_a - h mu^2) chi_a

# Wait, chi_a doesn't mix with sector B at tree level for the fermions.
# The fermionic mass for the (psi_a, Y~_a) pair:
# W_{psi_a, Y~_a} = h phi_0
# W_{Y_a, psi~_a} = h phi_0

# So in basis (psi_a, Y~_a):
# M_ferm = [[0, h phi_0], [h phi_0, 0]]
# Eigenvalues: +/- h phi_0

# And similarly for (psi~_a, Y_a).

print("\nFor each broken flavor a:")
m_link = h * phi_0
print(f"  Fermion masses: +/- h*phi_0 = +/- {m_link:.4f} MeV")
print(f"  (Dirac pair with mass h*phi_0)")

# Bosonic masses: The F-term for chi_a gives SUSY-breaking splitting.
# F_{chi_a} = h psi_a psi~_a + m_a - h mu^2
# At vacuum: F_{chi_a} = m_a - h mu^2 = -(h mu^2 - m_a) = -f_a

f_a = h * mu2 - m_sorted  # SUSY-breaking F-terms
print(f"\n  F-terms (f_a = h*mu^2 - m_a):")
for i in range(N_f):
    print(f"    f_{labels_sorted[i]} = {f_a[i]:.4f} MeV^2")

# The bosonic mass-squared matrix for the scalar components of (psi_a, Y~_a):
# M^2_B = M_ferm^dag M_ferm + SUSY-breaking corrections

# The SUSY part gives:
# M^2_SUSY = h^2 phi_0^2 for all scalar components

# The SUSY-breaking part:
# From F_{chi_a} = h psi_a psi~_a + m_a - h mu^2:
# d^2(F_{chi_a})/d(psi_a)d(psi~_a) = h
# So M^2_break = -f_a * h (for the psi_a-psi~_a cross term)

# For the scalar psi_a (two real components):
# M^2(Re psi_a) = h^2 phi_0^2 + h f_a  = h^2 phi_0^2 + h(h mu^2 - m_a)
# M^2(Im psi_a) = h^2 phi_0^2 - h f_a  = h^2 phi_0^2 - h(h mu^2 - m_a)

# Wait, this isn't quite right. Let me be more precise.
# The scalar mass-squared matrix for the complex fields (psi_a, psi~_a*):
# Comes from |F_{Y_a}|^2 + |F_{Y~_a}|^2 + |F_{chi_a}|^2

# |F_{Y_a}|^2 = |h phi_0 psi~_a|^2 → mass^2 = h^2 phi_0^2 for psi~_a
# |F_{Y~_a}|^2 = |h phi_0 psi_a|^2 → mass^2 = h^2 phi_0^2 for psi_a
# |F_{chi_a}|^2 = |h psi_a psi~_a - f_a|^2
#   Expand: f_a^2 - f_a h (psi_a psi~_a + psi_a* psi~_a*) + ...
#   The quadratic term gives a B-term: -f_a h (psi_a psi~_a + h.c.)

# In the (psi_a, psi~_a*) basis, the scalar mass-squared matrix is:
# M^2 = [[h^2 phi_0^2,  -h f_a],
#         [-h f_a,  h^2 phi_0^2]]

# Wait, I need to be more careful about the basis. Let me use
# (psi_a, psi~_a) as complex fields.

# Actually the simplest approach: there's a 2x2 mass-squared matrix
# for the scalars, mixing psi_a with psi~_a^* through the B-term.
# The eigenvalues are:
# m^2_pm = h^2 phi_0^2 +/- h f_a = h^2 phi_0^2 +/- h(h mu^2 - m_a)
#        = h^2 mu^2 - m_a +/- (h^2 mu^2 - h m_a)
# Hmm, this simplifies. Let me just compute numerically.

# Actually, the ISS scalar spectrum for the link sector is:
# For psi_a scalars (complex field, 2 real d.o.f.):
# m^2_+ = h^2 phi_0^2 + h f_a
# m^2_- = h^2 phi_0^2 - h f_a

# And for Y_a, Y~_a: these get mass h^2 phi_0^2 without splitting
# (they don't couple to the F-term directly).

# But actually Y_a and Y~_a mix with rho, rho~ through the superpotential.
# At leading order in the broken sector, Y_a gets mass from:
# W ⊃ h phi_0 Y_a psi~_a → F_{Y_a} = h phi_0 psi~_a, F_{psi~_a} = h phi_0 Y_a + ...

# So actually the 4-field system (psi_a, psi~_a, Y_a, Y~_a) has:
# Fermionic masses: h phi_0 (for psi_a-Y~_a and psi~_a-Y_a Dirac pairs)
# Scalar masses: same as fermionic + SUSY-breaking B-term splitting

# The B-term only affects (psi_a, psi~_a) through F_{chi_a}.
# Y_a and Y~_a don't get B-term splitting.

print("\nScalar masses in sector B (per broken flavor a):")
print("  From psi_a, psi~_a (with B-term from F_{chi_a}):")

for i in range(1, N_f):  # broken flavors
    a = i
    m2_plus = h**2 * phi_0_sq + h * f_a[a]
    m2_minus = h**2 * phi_0_sq - h * f_a[a]

    print(f"\n  Flavor '{labels_sorted[a]}' (m_a = {m_sorted[a]:.2f} MeV):")
    print(f"    f_a = h*mu^2 - m_a = {f_a[a]:.4f} MeV^2")
    print(f"    h^2 phi_0^2 = {h**2 * phi_0_sq:.4f} MeV^2")
    print(f"    m^2_+ = h^2*phi_0^2 + h*f_a = {m2_plus:.4f} MeV^2")
    print(f"    m^2_- = h^2*phi_0^2 - h*f_a = {m2_minus:.4f} MeV^2")
    print(f"    sqrt(m^2_+) = {sqrt(abs(m2_plus)):.4f} MeV")
    if m2_minus >= 0:
        print(f"    sqrt(m^2_-) = {sqrt(m2_minus):.4f} MeV")
    else:
        print(f"    m^2_- < 0 → TACHYONIC (|m| = {sqrt(abs(m2_minus)):.4f} MeV)")
        print(f"    This means Phi_0 = 0 is not the true minimum!")

# =============================================
print("\n" + "-" * 60)
print("SECTOR C: Broken mesons chi_a (diagonal, a=2,3,4)")
print("-" * 60)

# chi_a = Phi^a_a for a in the broken block.
# Superpotential: W ⊃ (m_a - h mu^2) chi_a + h psi_a chi_a psi~_a
# At the vacuum (psi_a = psi~_a = 0):
# W_chi_a = m_a - h mu^2 (linear, gives F-term)
# W_{chi_a, chi_a} = 0 (no mass term)
# W_{chi_a, psi_a} = h psi~_a|_vac = 0
# So chi_a has ZERO fermion mass and ZERO scalar mass at tree level.
# It's another flat direction (or rather, the F-term is constant in chi_a).

# Actually, chi_a gets a 1-loop CW mass. At tree level it's massless.

# But wait: for each flavor a, chi_a has the coupling:
# W ⊃ h psi_a chi_a psi~_a
# At the vacuum, psi_a = psi~_a = 0, so chi_a is indeed massless at tree level.

# The F-term F_{chi_a} = m_a - h mu^2 is a CONSTANT (not field-dependent
# to leading order). So chi_a is a pseudo-Goldstone boson of the
# approximate R-symmetry.

# These are the SEIBERG SEESAW fields: their 1-loop CW mass will be
# proportional to m_a, implementing the seesaw M_meson ~ m_a or ~ 1/m_a.

print("""
The broken mesons chi_a are MASSLESS at tree level.
They get mass from the Coleman-Weinberg potential at one loop.
F-term: F_{chi_a} = m_a - h*mu^2 (constant, nonzero → SUSY broken)
""")

# =============================================
# COMPLETE TREE-LEVEL SPECTRUM
# =============================================
print("\n" + "=" * 60)
print("COMPLETE TREE-LEVEL SPECTRUM (at Phi_0 = X = 0)")
print("=" * 60)

# Let's now parametrize by Phi_0 = X (the pseudo-modulus).
# The masses depend on X.

# For generality, let X be a free parameter.
# Then the fermionic mass matrix has:
# - Sector A: 0 (Goldstino), sqrt(2) h sqrt(|X|^2 + phi_0^2) [approximately]
# Actually at X = 0: the sector A fermion mass is sqrt(2) h phi_0.

# For the key question about CW and Koide, the important spectrum is:
# The psi_a sector, which gives BOTH the dominant loop contributions AND
# the SUSY-breaking splitting.

# Let's compute with X as a parameter (but set X = 0 for now).

X = 0.0  # pseudo-modulus VEV (will be varied later)

print(f"\nWith X = {X} MeV:")
print(f"  phi_0 = {phi_0:.4f} MeV")

print("\nField content and masses:")
print(f"{'Field':<20} {'Boson m²':<25} {'Fermion m':<20} {'Notes':<30}")
print("-" * 95)

# Sector A
print(f"{'X (pseudo-mod)':<20} {'0 (flat)':<25} {'0 (Goldstino)':<20} {'CW lifts this':<30}")
m_A_heavy = sqrt(2) * h * phi_0
print(f"{'rho+rho~ heavy':<20} {h**2*(2*phi_0_sq):<25.4f} {m_A_heavy:<20.4f} {'SUSY pair':<30}")

# Sector B: for each broken flavor
all_boson_masses_sq = []
all_fermion_masses = []

for i in range(1, N_f):
    a = i
    label = labels_sorted[a]

    # Fermion mass
    m_ferm = h * phi_0  # Dirac mass for psi_a-Y~_a pair

    # Boson masses
    m2_plus = h**2 * phi_0_sq + h * f_a[a]
    m2_minus = h**2 * phi_0_sq - h * f_a[a]

    # Also Y_a, Y~_a bosons: mass^2 = h^2 phi_0^2 (no splitting)
    m2_Y = h**2 * phi_0_sq

    print(f"{'psi_'+label+' (B+)':<20} {m2_plus:<25.4f} {m_ferm:<20.4f} {'SUSY-split boson+':<30}")
    print(f"{'psi_'+label+' (B-)':<20} {m2_minus:<25.4f} {'':<20} {'SUSY-split boson-':<30}")
    print(f"{'Y_'+label:<20} {m2_Y:<25.4f} {m_ferm:<20.4f} {'SUSY pair':<30}")
    print(f"{'Y~_'+label:<20} {m2_Y:<25.4f} {'':<20} {'SUSY pair':<30}")

    # Store for CW computation
    # Each broken flavor contributes:
    # 2 complex scalars from psi_a, psi~_a with mass^2 = m2_plus, m2_minus
    # 2 complex scalars from Y_a, Y~_a with mass^2 = m2_Y (each)
    # 2 Dirac fermions with mass = h*phi_0 (each pair)
    all_boson_masses_sq.extend([m2_plus, m2_minus, m2_Y, m2_Y])
    all_fermion_masses.extend([m_ferm, m_ferm])

# Sector C
for i in range(1, N_f):
    label = labels_sorted[i]
    print(f"{'chi_'+label+' (broken)':<20} {'0 (flat)':<25} {'0':<20} {'CW lifts this':<30}")

# ============================================================
# 2. COLEMAN-WEINBERG POTENTIAL
# ============================================================
print("\n\n" + "=" * 80)
print("2. COLEMAN-WEINBERG POTENTIAL")
print("=" * 80)

print("""
V_CW = (1/64pi^2) STr[M^4 (log(M^2/Lambda^2) - 3/2)]

where STr = sum over bosons - sum over fermions,
weighted by (2s+1) for each spin-s particle.

For chiral multiplets:
  STr[M^2n] = Tr[M_B^2n] - 2 Tr[M_F^2n]
  (complex scalar has 2 real d.o.f., Weyl fermion has 2)

The CW potential lifts the pseudo-modulus X and the broken mesons chi_a.
""")

# For the CW computation, we need the X-dependent masses.
# Let's parametrize by X and compute V_CW(X, chi_a).

# The chi_a-dependent masses come from the coupling:
# W ⊃ h psi_a chi_a psi~_a
# This modifies the mass spectrum when chi_a ≠ 0.
# However, at leading order, chi_a enters only at QUADRATIC level
# in the mass matrices (it doesn't have a VEV at tree level).
# So the CW potential for chi_a starts at O(chi_a^2) → CW mass for chi_a.

# For the X-dependent CW potential:
# The psi_a sector masses depend on X through the fermion mass h phi_0
# being replaced by h sqrt(phi_0^2 + X^2/2 + ...) or more precisely
# through the full mass matrix.

# Actually, in the ISS model with X as the pseudo-modulus,
# the spectrum is X-dependent through the fermion mass matrix.
# The key X-dependent piece is in sector A:
# W_{rho,X} = h phi_0, and if X shifts, the mass eigenvalues change.

# But the DOMINANT CW contribution to the chi_a mass comes from
# sector B loops. The chi_a couples to psi_a, psi~_a through:
# W ⊃ h psi_a chi_a psi~_a

# One-loop CW mass for chi_a from psi_a loop:
# m^2_CW(chi_a) = (h^2 / 16pi^2) * [function of h phi_0 and f_a]

# The standard ISS result (Shih 2007):
# m^2_CW(chi_a) = (h^4 / 16pi^2) * f_a^2 * G(x_a)
# where x_a = h^2 phi_0^2 / f_a  and G is a loop function.

# But more precisely, the CW potential depends on the FULL field-dependent
# mass matrix. Let me compute it properly.

# =============================================
# CW POTENTIAL FOR chi_a (the meson pseudo-moduli)
# =============================================

print("\n" + "-" * 60)
print("CW potential for chi_a (broken mesons)")
print("-" * 60)

# When chi_a gets a VEV, the psi_a sector mass matrix changes.
# With <chi_a> = phi_a:
# W_B ⊃ h phi_0 Y_a psi~_a + h psi_a Y~_a phi_0 + h phi_a psi_a psi~_a

# The fermion mass matrix in the (psi_a, Y~_a) x (psi~_a, Y_a) basis:
# M_F = [[h phi_a, h phi_0],
#         [h phi_0, 0]]

# Eigenvalues: lambda^2 - h phi_a lambda - h^2 phi_0^2 = 0
# lambda = (h phi_a +/- sqrt(h^2 phi_a^2 + 4 h^2 phi_0^2)) / 2

# For phi_a = 0: lambda = +/- h phi_0 (as before)

# The boson mass-squared matrix (including SUSY-breaking B-term):
# In (psi_a, psi~_a*) basis:
# M^2_B = M_F^dag M_F + [[0, -h f_a], [-h f_a, 0]]
#        = [[h^2(phi_a^2 + phi_0^2), -h f_a + h^2 phi_0 phi_a*],
#           [-h f_a + h^2 phi_0 phi_a, h^2 phi_0^2]]

# Wait, I need to be more careful. Let me redo this with explicit
# computation of the field-dependent mass matrices.

# For each flavor a, the relevant fields are:
# Chiral superfields: Psi_a, Psi~_a, Y_a, Y~_a
# with Kahler potential K = |Psi_a|^2 + |Psi~_a|^2 + |Y_a|^2 + |Y~_a|^2

# Superpotential (relevant terms):
# W = h phi_0 (Y_a Psi~_a + Psi_a Y~_a) + h chi_a Psi_a Psi~_a + const*chi_a
# where const = m_a - h mu^2 = -f_a/h (using f_a = h mu^2 - m_a)

# Hmm wait, let me reconsider. We have:
# F_{chi_a} = dW/d(chi_a) = h psi_a psi~_a + m_a - h mu^2

# This is the F-term EQUATION. At the vacuum, psi_a = psi~_a = 0,
# so F_{chi_a} = m_a - h mu^2 = -(h mu^2 - m_a) ≠ 0.
# This is the source of SUSY breaking.

# The CW effective potential for chi_a is obtained by integrating out
# the massive fields (psi_a, psi~_a, Y_a, Y~_a) at one loop.

# With chi_a = phi_a (generic VEV), the fermion mass matrix is:
# In the basis (Psi_a, Y~_a):
#   M_F = [[h phi_a, h phi_0],
#           [h phi_0, 0       ]]
# (and similarly for (Psi~_a, Y_a))

# This gives two Dirac fermions with masses:
# m^2_F = eigenvalues of M_F^dag M_F
# = eigenvalues of [[h^2(|phi_a|^2 + phi_0^2), h^2 phi_0 phi_a*],
#                    [h^2 phi_0 phi_a, h^2 phi_0^2]]

# Eigenvalues:
# m^2_F = (h^2/2) * (|phi_a|^2 + 2 phi_0^2 +/- sqrt(|phi_a|^4 + 4 phi_0^2 |phi_a|^2))
# = (h^2/2) * (|phi_a|^2 + 2 phi_0^2 +/- |phi_a| sqrt(|phi_a|^2 + 4 phi_0^2))

# For the bosonic mass-squared, including the SUSY-breaking B-term from F_{chi_a}:
# The scalar potential includes:
# V ⊃ |F_{Y_a}|^2 + |F_{Y~_a}|^2 + |F_{chi_a}|^2
# = |h phi_0 Psi~_a|^2 + |h phi_0 Psi_a|^2 + |h Psi_a Psi~_a + m_a - h mu^2|^2
# = h^2 phi_0^2 (|Psi_a|^2 + |Psi~_a|^2) + |h Psi_a Psi~_a - f_a/h|^2 * h^2

# Hmm let me be very explicit. The F-terms are:
# F_{Y_a} = h phi_0 Psi~_a
# F_{Y~_a} = h phi_0 Psi_a  (these don't include chi_a coupling)
# Wait, I need to include ALL couplings.

# Actually, the relevant superpotential for the chi_a subsystem is:
# W_a = h phi_0 Y_a Psi~_a + h phi_0 Psi_a Y~_a + h chi_a Psi_a Psi~_a
#       + (m_a - h mu^2) chi_a

# F-terms:
# F_{chi_a} = h Psi_a Psi~_a + m_a - h mu^2
# F_{Psi_a} = h phi_0 Y~_a + h chi_a Psi~_a
# F_{Psi~_a} = h phi_0 Y_a + h chi_a Psi_a
# F_{Y_a} = h phi_0 Psi~_a
# F_{Y~_a} = h phi_0 Psi_a

# The scalar potential:
# V = |F_{chi_a}|^2 + |F_{Psi_a}|^2 + |F_{Psi~_a}|^2 + |F_{Y_a}|^2 + |F_{Y~_a}|^2

# At the vacuum (Psi_a = Psi~_a = Y_a = Y~_a = 0):
# V_0 = |m_a - h mu^2|^2 = f_a^2 / h^2 ... wait.
# f_a = h mu^2 - m_a, so (m_a - h mu^2)^2 = f_a^2
# V_0 = f_a^2 (in the convention where W is dimensionless... no.)

# OK let me just compute V numerically. The scalar mass-squared matrix is:
# M^2_{ij} = d^2 V / d phi_i d phi_j*
# where phi = (Psi_a, Psi~_a, Y_a, Y~_a) are complex scalars.

# For the CW potential of chi_a, I treat chi_a as a BACKGROUND field
# and integrate out (Psi_a, Psi~_a, Y_a, Y~_a).

# Let me compute the scalar mass-squared matrix as a function of chi_a.

def compute_mass_matrices(phi_a_val, phi_0_val, f_a_val, h_val):
    """
    Compute the scalar and fermion mass-squared matrices
    for the sector-B fields (Psi_a, Psi~_a, Y_a, Y~_a)
    as functions of the background chi_a = phi_a_val.

    Returns (boson_mass_sq_eigenvalues, fermion_mass_sq_eigenvalues)
    """
    # Fermion mass matrix W_ij in basis (Psi_a, Psi~_a, Y_a, Y~_a):
    # W_{Psi_a, Psi~_a} = h chi_a = h phi_a
    # W_{Psi_a, Y~_a} = h phi_0
    # W_{Psi~_a, Y_a} = h phi_0
    # All others zero.

    # Actually, for Weyl fermions, the mass matrix is:
    # W_ij = d^2 W / d Phi_i d Phi_j
    # In basis (Psi_a, Psi~_a, Y_a, Y~_a):
    hp = h_val * phi_a_val
    hp0 = h_val * phi_0_val

    W_ferm = np.array([
        [0,   hp,  0,   hp0],
        [hp,  0,   hp0, 0  ],
        [0,   hp0, 0,   0  ],
        [hp0, 0,   0,   0  ]
    ], dtype=complex)

    # Fermion mass-squared = eigenvalues of W^dag W
    WdW = W_ferm.conj().T @ W_ferm
    ferm_m2 = np.sort(np.real(np.linalg.eigvalsh(WdW)))

    # Scalar mass-squared matrix:
    # V = sum_A |F_A|^2
    # d^2 V / d phi_i d phi_j* = sum_A (d F_A / d phi_i)* (d F_A / d phi_j)
    #                            + sum_A F_A* (d^2 F_A / d phi_j* d phi_i)
    #   [second term is zero for holomorphic F-terms]
    #                            + sum_A (d^2 F_A* / d phi_i d phi_j*) F_A
    #   [this is also zero]
    # Wait, the scalar potential is:
    # V = sum_A |F_A(phi)|^2 where F_A = dW/dPhi_A
    # and the Phi_A here include ALL fields: chi_a, Psi_a, Psi~_a, Y_a, Y~_a.

    # But chi_a is the BACKGROUND field. We differentiate w.r.t.
    # phi_i = {Psi_a, Psi~_a, Y_a, Y~_a} only.

    # d^2 V / d phi_i d phi_j* = sum_A (dF_A/dphi_i)^* (dF_A/dphi_j)
    #   [this is the "holomorphic" piece = W_{iA}^* W_{jA} = (W^dag W)_{ij}]
    # + "B-term" piece from F_{chi_a} being nonzero:
    # V ⊃ |F_{chi_a}|^2 = |h Psi_a Psi~_a + m_a - h mu^2|^2
    # d^2/d Psi_a d Psi_a* = |h Psi~_a|^2 → 0 at vacuum
    # d^2/d Psi_a d Psi~_a = F_{chi_a}^* * h (the B-term!)

    # So the full scalar mass-squared matrix is:
    # M^2_S = (W^dag W)_ij  + B-term corrections

    # The (W^dag W) part:
    # In basis (Psi_a, Psi~_a, Y_a, Y~_a):
    # dF_{chi_a}/dPsi_a = h Psi~_a → 0 at vacuum  [F_{chi_a} doesn't contribute to W^dag W in leading order]
    # dF_{Psi_a}/dPsi~_a = h chi_a, dF_{Psi_a}/dY~_a = h phi_0
    # dF_{Psi~_a}/dPsi_a = h chi_a, dF_{Psi~_a}/dY_a = h phi_0
    # dF_{Y_a}/dPsi~_a = h phi_0
    # dF_{Y~_a}/dPsi_a = h phi_0

    # So dF_A/dphi_j matrix (rows = A = {chi_a, Psi_a, Psi~_a, Y_a, Y~_a},
    #                        cols = j = {Psi_a, Psi~_a, Y_a, Y~_a}):
    # Note: evaluated at vacuum (Psi_a = Psi~_a = Y_a = Y~_a = 0):

    dF = np.array([
        [0,    0,    0,   0  ],   # dF_{chi_a}/d{Psi_a, Psi~_a, Y_a, Y~_a} = {h*0, h*0, 0, 0}
        [0,    hp,   0,   hp0],   # dF_{Psi_a}/d{...} = {0, h*chi_a, 0, h*phi_0}
        [hp,   0,    hp0, 0  ],   # dF_{Psi~_a}/d{...} = {h*chi_a, 0, h*phi_0, 0}
        [0,    hp0,  0,   0  ],   # dF_{Y_a}/d{...} = {0, h*phi_0, 0, 0}
        [hp0,  0,    0,   0  ],   # dF_{Y~_a}/d{...} = {h*phi_0, 0, 0, 0}
    ], dtype=complex)

    M2_hol = (dF.conj().T @ dF)  # (W^dag W) piece

    # B-term: from V ⊃ |F_{chi_a}|^2 where F_{chi_a} = h Psi_a Psi~_a - f_a
    # d^2|F_{chi_a}|^2 / d Psi_a d Psi~_a^* = (F_{chi_a})^* * h* = (-f_a)^* * h = -f_a * h
    # (at vacuum where Psi = 0)
    # Actually, V = |F_{chi_a}|^2 = (h Psi_a Psi~_a - f_a)^*(h Psi_a Psi~_a - f_a)
    # d/dPsi_a = (h Psi~_a)*(h Psi_a Psi~_a - f_a) = (h Psi~_a)*(h Psi_a Psi~_a - f_a)
    # d^2/(dPsi_a dPsi~_a*) = h* (h Psi_a Psi~_a - f_a)|_vac = h*(-f_a) = -h f_a
    # Hmm, but this is d^2 V / d Psi_a d (Psi~_a)^*
    # In the (Psi_a, Psi~_a, Y_a, Y~_a) x (Psi_a*, Psi~_a*, Y_a*, Y~_a*) matrix,
    # this goes in position (0,1) [Psi_a row, Psi~_a* column].

    # BUT WAIT: the mass-squared matrix has TWO types of terms:
    # (1) M^2_{i,j*} = d^2V / d phi_i d phi_j* (the "diagonal" block)
    # (2) B_{i,j} = d^2V / d phi_i d phi_j   (the "off-diagonal" B-term block)

    # For real scalars, the full mass matrix is:
    # In basis (Re phi_1, Im phi_1, Re phi_2, ...):
    # M^2 = [[Re(M^2 + B), Im(B - M^2)],
    #         [Im(M^2 + B), Re(M^2 - B)]]

    # For our case, the B-term comes from:
    # d^2|F_{chi_a}|^2 / d Psi_a d Psi~_a
    # = F_{chi_a}^* * (d^2 F_{chi_a} / d Psi_a d Psi~_a)
    # = (-f_a) * h
    # So B_{Psi_a, Psi~_a} = -h f_a_val

    B_matrix = np.zeros((4, 4), dtype=complex)
    B_matrix[0, 1] = -h_val * f_a_val  # Psi_a, Psi~_a
    B_matrix[1, 0] = -h_val * f_a_val  # Psi~_a, Psi_a (symmetric)

    # The full real scalar mass matrix (8x8 for 4 complex fields):
    # Using the block form:
    # M_real = [[Re(M2_hol + B), -Im(M2_hol + B)],
    #           [Im(M2_hol + B),  Re(M2_hol - B)]]

    # Hmm, the standard form is different. Let me use a cleaner approach.
    # For each complex scalar phi = (a + ib)/sqrt(2):
    # The mass terms are: phi^* M^2 phi + (1/2)(phi phi B + h.c.)
    # = a^2 (Re M^2 + Re B)/2 + b^2 (Re M^2 - Re B)/2 + a b Im B + ...

    # Actually, the simplest approach: construct the full 8x8 mass matrix
    # in the real basis (a_1, b_1, a_2, b_2, a_3, b_3, a_4, b_4)
    # where phi_k = (a_k + i b_k) / sqrt(2).

    n = 4  # number of complex fields
    M_real = np.zeros((2*n, 2*n))

    for i in range(n):
        for j in range(n):
            re_m2 = np.real(M2_hol[i, j])
            im_m2 = np.imag(M2_hol[i, j])
            re_b = np.real(B_matrix[i, j])
            im_b = np.imag(B_matrix[i, j])

            # (a_i, a_j) block: Re(M^2) + Re(B)
            M_real[2*i, 2*j] = re_m2 + re_b
            # (b_i, b_j) block: Re(M^2) - Re(B)
            M_real[2*i+1, 2*j+1] = re_m2 - re_b
            # (a_i, b_j) block: Im(B) - Im(M^2)
            M_real[2*i, 2*j+1] = im_b - im_m2
            # (b_i, a_j) block: Im(B) + Im(M^2)
            M_real[2*i+1, 2*j] = im_b + im_m2

    boson_m2 = np.sort(np.real(np.linalg.eigvalsh(M_real)))

    return boson_m2, ferm_m2


# Test at chi_a = 0
print("\nMass spectrum at chi_a = 0 for each broken flavor:")
print(f"{'Flavor':<8} {'Boson m^2 (8 vals)':<60} {'Fermion m^2 (4 vals)'}")
print("-" * 120)

for i in range(1, N_f):
    a = i
    label = labels_sorted[a]
    bm2, fm2 = compute_mass_matrices(0.0, phi_0, f_a[a], h)
    print(f"{label:<8} {str(np.round(bm2, 2)):<60} {str(np.round(fm2, 2))}")

    # Check supertrace
    str_m2 = np.sum(bm2) - 2 * np.sum(fm2)  # factor 2 for Weyl fermion d.o.f.
    print(f"         STr[M^2] = {str_m2:.6f} (should be 0 at tree level: {abs(str_m2) < 1e-6})")

# =============================================
# CW effective potential as function of chi_a
# =============================================
print("\n\n" + "=" * 80)
print("CW EFFECTIVE POTENTIAL FOR chi_a")
print("=" * 80)

def V_CW_from_masses(boson_m2, fermion_m2, Lambda_sq):
    """
    Compute V_CW = (1/64pi^2) STr[M^4 (log(M^2/Lambda^2) - 3/2)]

    boson_m2: array of bosonic mass-squared eigenvalues
    fermion_m2: array of fermionic mass-squared eigenvalues
    Lambda_sq: renormalization scale squared
    """
    V = 0.0

    # Bosonic contribution
    for m2 in boson_m2:
        if m2 > 0:
            V += m2**2 * (np.log(m2 / Lambda_sq) - 1.5)
        elif m2 < 0:
            # Tachyonic — use |m^2| and add imaginary part
            V += m2**2 * (np.log(abs(m2) / Lambda_sq) - 1.5)

    # Fermionic contribution (with -2 for Weyl fermion counting)
    # Each Weyl fermion contributes -1 to STr
    # But our fermion_m2 has eigenvalues from W^dag W (4x4),
    # which gives 4 eigenvalues. Each eigenvalue corresponds to
    # one Dirac fermion = 2 Weyl fermions.
    # Actually, W is 4x4 symmetric, so W^dag W has 4 eigenvalues.
    # Each nonzero eigenvalue corresponds to a Dirac fermion.
    for m2 in fermion_m2:
        if m2 > 1e-30:
            V -= 2 * m2**2 * (np.log(m2 / Lambda_sq) - 1.5)

    return V / (64.0 * pi**2)


# Compute CW mass for chi_a by numerical differentiation
Lambda_sq = Lambda**2

print("\nCW masses for broken mesons chi_a:")
print(f"{'Flavor':<8} {'m_a (MeV)':<12} {'m^2_CW (MeV^2)':<20} {'m_CW (MeV)':<15} {'m_CW/m_a':<12}")
print("-" * 70)

cw_masses_sq = []
cw_masses = []

for i in range(1, N_f):
    a = i
    label = labels_sorted[a]

    # Numerical second derivative of V_CW w.r.t. chi_a
    delta = 0.001  # MeV, small perturbation

    V_plus = 0.0
    V_minus = 0.0
    V_0 = 0.0

    bm2_0, fm2_0 = compute_mass_matrices(0.0, phi_0, f_a[a], h)
    bm2_p, fm2_p = compute_mass_matrices(delta, phi_0, f_a[a], h)
    bm2_m, fm2_m = compute_mass_matrices(-delta, phi_0, f_a[a], h)

    V_0 = V_CW_from_masses(bm2_0, fm2_0, Lambda_sq)
    V_p = V_CW_from_masses(bm2_p, fm2_p, Lambda_sq)
    V_m = V_CW_from_masses(bm2_m, fm2_m, Lambda_sq)

    m2_cw = (V_p + V_m - 2*V_0) / delta**2

    cw_masses_sq.append(m2_cw)
    if m2_cw > 0:
        m_cw = sqrt(m2_cw)
    else:
        m_cw = -sqrt(abs(m2_cw))  # tachyonic
    cw_masses.append(m_cw)

    print(f"{label:<8} {m_sorted[a]:<12.2f} {m2_cw:<20.6f} {m_cw:<15.6f} {m_cw/m_sorted[a]:<12.8f}")

# ============================================================
# 3. KOIDE ANALYSIS
# ============================================================
print("\n\n" + "=" * 80)
print("3. KOIDE ANALYSIS: UV vs IR")
print("=" * 80)

def koide_Q(m1, m2, m3):
    """Compute Koide quotient Q = (m1+m2+m3) / (sqrt(m1)+sqrt(m2)+sqrt(m3))^2"""
    s = sqrt(abs(m1)) + sqrt(abs(m2)) + sqrt(abs(m3))
    return (abs(m1) + abs(m2) + abs(m3)) / s**2

# UV quark masses
print("\nUV quark masses:")
print(f"  m_u = {m_u} MeV, m_d = {m_d} MeV, m_s = {m_s} MeV, m_c = {m_c} MeV")

Q_uv_uds = koide_Q(m_u, m_d, m_s)
Q_uv_dsc = koide_Q(m_d, m_s, m_c)
Q_uv_neg_scb = koide_Q(m_s, m_c, m_s)  # placeholder
print(f"\n  Q(m_u, m_d, m_s) = {Q_uv_uds:.6f}")
print(f"  Q(m_d, m_s, m_c) = {Q_uv_dsc:.6f}")

# Tree-level meson masses (Seiberg seesaw)
# M_meson ~ h^2 mu^2 / m_quark  (up to corrections)
# Actually the tree-level meson is MASSLESS; the CW gives it a mass.
# The CW mass is the physical meson mass.

# But the SEIBERG SEESAW is about the SUSY vacuum meson masses:
# <M_i> = Lambda^2 / m_i  (in the SUSY vacuum, which is the WRONG vacuum)
# In the ISS metastable vacuum, the meson VEVs are different.

# The ISS one-loop meson masses are (from the CW computation above):
# m^2_CW(chi_a) ∝ function of (h, phi_0, f_a)

# The KEY question: how does m_CW depend on the quark mass m_a?
# From the structure: f_a = h mu^2 - m_a, and m_a << h mu^2, so f_a ≈ h mu^2.
# The mass matrices depend on m_a ONLY through f_a.
# Since f_a ≈ h mu^2 (1 - m_a/(h mu^2)), the correction is:
# delta(f_a)/f_a = -m_a/(h mu^2) << 1 for all quarks.

# So the CW mass is APPROXIMATELY flavor-universal!
# The flavor-dependent correction is:
# delta m^2_CW / m^2_CW ≈ -2 m_a / (h mu^2)

print("\n\nFlavor-dependence of CW mass:")
print(f"  h mu^2 = {h * mu2:.4f} MeV^2")
for i in range(1, N_f):
    label = labels_sorted[i]
    correction = m_sorted[i] / (h * mu2)
    print(f"  m_{label}/(h mu^2) = {correction:.6f} ({correction*100:.4f}%)")

print("""
Since m_a << h*mu^2 for u, d, s (but NOT for c!),
the CW mass is approximately flavor-universal for the light quarks
but has a significant correction for charm.

This means the Koide ratio for the CW-corrected spectrum
will be shifted from the UV quark mass Koide ratio.
""")

# ============================================================
# 4. NUMERICAL EVALUATION
# ============================================================
print("\n" + "=" * 80)
print("4. DETAILED NUMERICAL EVALUATION")
print("=" * 80)

# Compute CW masses more carefully with better numerical differentiation
print("\nHigh-precision CW mass computation:")
print(f"  Using h = sqrt(3) = {h:.10f}")
print(f"  mu = {mu} MeV, mu^2 = {mu2} MeV^2")
print(f"  phi_0 = sqrt(mu^2 - m_u/h) = {phi_0:.10f} MeV")
print(f"  Lambda = {Lambda} MeV")

# Use Richardson extrapolation for numerical derivative
deltas = [0.01, 0.005, 0.001, 0.0005, 0.0001]

print("\nConvergence test for CW mass (flavor d):")
for delta in deltas:
    bm2_0, fm2_0 = compute_mass_matrices(0.0, phi_0, f_a[1], h)
    bm2_p, fm2_p = compute_mass_matrices(delta, phi_0, f_a[1], h)
    bm2_m, fm2_m = compute_mass_matrices(-delta, phi_0, f_a[1], h)
    V_0 = V_CW_from_masses(bm2_0, fm2_0, Lambda_sq)
    V_p = V_CW_from_masses(bm2_p, fm2_p, Lambda_sq)
    V_m = V_CW_from_masses(bm2_m, fm2_m, Lambda_sq)
    m2_test = (V_p + V_m - 2*V_0) / delta**2
    print(f"  delta = {delta:.5f}: m^2_CW = {m2_test:.10f} MeV^2")

# Final computation with optimal delta
delta = 0.001

print(f"\n\nFinal CW masses (delta = {delta}):")
print(f"{'Flavor':<8} {'m_quark':<12} {'f_a':<15} {'m^2_CW':<20} {'m_CW':<15} {'ratio m_CW/m_q':<15}")
print("-" * 85)

cw_m2_final = {}
cw_m_final = {}

for i in range(N_f):
    label = labels[i]
    m_q = m_quarks[i]
    f_i = h * mu2 - m_q

    bm2_0, fm2_0 = compute_mass_matrices(0.0, phi_0, f_i, h)
    bm2_p, fm2_p = compute_mass_matrices(delta, phi_0, f_i, h)
    bm2_m, fm2_m = compute_mass_matrices(-delta, phi_0, f_i, h)

    V_0 = V_CW_from_masses(bm2_0, fm2_0, Lambda_sq)
    V_p = V_CW_from_masses(bm2_p, fm2_p, Lambda_sq)
    V_m = V_CW_from_masses(bm2_m, fm2_m, Lambda_sq)

    m2_cw_i = (V_p + V_m - 2*V_0) / delta**2

    if m2_cw_i > 0:
        m_cw_i = sqrt(m2_cw_i)
    else:
        m_cw_i = -sqrt(abs(m2_cw_i))

    cw_m2_final[label] = m2_cw_i
    cw_m_final[label] = m_cw_i

    print(f"{label:<8} {m_q:<12.2f} {f_i:<15.4f} {m2_cw_i:<20.8f} {m_cw_i:<15.8f} {m_cw_i/m_q:<15.8f}")

# Supertrace analysis
print("\n\nSupertrace analysis:")
print("STr[M^{2n}] for each flavor:")
for i in range(N_f):
    label = labels[i]
    m_q = m_quarks[i]
    f_i = h * mu2 - m_q

    bm2, fm2 = compute_mass_matrices(0.0, phi_0, f_i, h)

    for n in [0, 1, 2]:
        str_n = np.sum(bm2**n) - 2 * np.sum(fm2**n)
        print(f"  STr[M^{2*n}] (flavor {label}): {str_n:.8f}")
    print()

# ============================================================
# Koide ratio for CW spectrum
# ============================================================
print("\n" + "-" * 60)
print("KOIDE RATIO FOR CW-CORRECTED SPECTRUM")
print("-" * 60)

m_cw_u = abs(cw_m_final['u'])
m_cw_d = abs(cw_m_final['d'])
m_cw_s = abs(cw_m_final['s'])
m_cw_c = abs(cw_m_final['c'])

print(f"\nCW meson masses:")
print(f"  m_CW(u) = {m_cw_u:.8f} MeV")
print(f"  m_CW(d) = {m_cw_d:.8f} MeV")
print(f"  m_CW(s) = {m_cw_s:.8f} MeV")
print(f"  m_CW(c) = {m_cw_c:.8f} MeV")

Q_cw_uds = koide_Q(m_cw_u, m_cw_d, m_cw_s)
Q_cw_dsc = koide_Q(m_cw_d, m_cw_s, m_cw_c)
Q_cw_usc = koide_Q(m_cw_u, m_cw_s, m_cw_c)
Q_cw_udc = koide_Q(m_cw_u, m_cw_d, m_cw_c)

print(f"\nKoide quotients for CW masses:")
print(f"  Q(u,d,s)_CW = {Q_cw_uds:.8f}")
print(f"  Q(d,s,c)_CW = {Q_cw_dsc:.8f}")
print(f"  Q(u,s,c)_CW = {Q_cw_usc:.8f}")
print(f"  Q(u,d,c)_CW = {Q_cw_udc:.8f}")

# Compare with UV Koide ratios
Q_uv_uds_v = koide_Q(m_u, m_d, m_s)
Q_uv_dsc_v = koide_Q(m_d, m_s, m_c)
Q_uv_usc_v = koide_Q(m_u, m_s, m_c)
Q_uv_udc_v = koide_Q(m_u, m_d, m_c)

print(f"\nKoide quotients for UV quark masses:")
print(f"  Q(u,d,s)_UV = {Q_uv_uds_v:.8f}")
print(f"  Q(d,s,c)_UV = {Q_uv_dsc_v:.8f}")
print(f"  Q(u,s,c)_UV = {Q_uv_usc_v:.8f}")
print(f"  Q(u,d,c)_UV = {Q_uv_udc_v:.8f}")

print(f"\nShift in Q:")
print(f"  Delta Q(u,d,s) = {Q_cw_uds - Q_uv_uds_v:.8f}")
print(f"  Delta Q(d,s,c) = {Q_cw_dsc - Q_uv_dsc_v:.8f}")
print(f"  Delta Q(u,s,c) = {Q_cw_usc - Q_uv_usc_v:.8f}")
print(f"  Delta Q(u,d,c) = {Q_cw_udc - Q_uv_udc_v:.8f}")

# Is the shift flavor-universal or flavor-dependent?
# If CW mass = const * f(m_a), and f is approximately constant,
# then the meson spectrum is approximately flavor-universal
# and the Koide Q is 1/3 (three equal masses).

# But the correction term is proportional to m_a/(h mu^2),
# so the meson spectrum DOES depend on quark mass, and the question
# is whether this dependence preserves or violates Koide.

# ============================================================
# Fractional shift analysis
# ============================================================
print("\n\n" + "-" * 60)
print("FRACTIONAL CW MASS SHIFT ANALYSIS")
print("-" * 60)

# The CW mass depends on the quark mass through f_a = h mu^2 - m_a.
# Let's expand m^2_CW(m_a) around m_a = 0:
# m^2_CW(m_a) = m^2_CW(0) + m_a * dm^2_CW/dm_a + ...
# dm^2_CW/dm_a = dm^2_CW/df_a * df_a/dm_a = -dm^2_CW/df_a

# Compute dm^2_CW/dm_a numerically
print("\nDerivative of CW mass w.r.t. quark mass:")
delta_m = 0.01  # MeV

for i in range(N_f):
    label = labels[i]
    m_q = m_quarks[i]

    f_0 = h * mu2 - m_q
    f_p = h * mu2 - (m_q + delta_m)
    f_m = h * mu2 - (m_q - delta_m)

    bm2_0, fm2_0 = compute_mass_matrices(0.0, phi_0, f_0, h)
    bm2_p, fm2_p = compute_mass_matrices(0.0, phi_0, f_p, h)
    bm2_m, fm2_m = compute_mass_matrices(0.0, phi_0, f_m, h)

    # CW mass for chi_a (use second derivative of V_CW w.r.t. chi_a at chi_a=0)
    def get_cw_mass_sq(f_val):
        d = 0.001
        bm0, fm0 = compute_mass_matrices(0.0, phi_0, f_val, h)
        bmp, fmp = compute_mass_matrices(d, phi_0, f_val, h)
        bmm, fmm = compute_mass_matrices(-d, phi_0, f_val, h)
        V0 = V_CW_from_masses(bm0, fm0, Lambda_sq)
        Vp = V_CW_from_masses(bmp, fmp, Lambda_sq)
        Vm = V_CW_from_masses(bmm, fmm, Lambda_sq)
        return (Vp + Vm - 2*V0) / d**2

    m2_cw_0 = get_cw_mass_sq(f_0)
    m2_cw_p = get_cw_mass_sq(f_p)
    m2_cw_m = get_cw_mass_sq(f_m)

    dm2_dm = (m2_cw_p - m2_cw_m) / (2 * delta_m)  # d(m^2_CW)/d(m_quark)

    # Note: increasing m_a decreases f_a, so dm^2_CW/dm_a < 0 expected
    # (heavier quarks → less SUSY breaking → smaller CW mass)

    frac_shift = dm2_dm * m_q / m2_cw_0 if abs(m2_cw_0) > 1e-30 else 0

    print(f"  d(m^2_CW)/dm_{label} = {dm2_dm:.8f}")
    print(f"    Fractional: m_q * dm^2/(m^2 * dm_q) = {frac_shift:.8f}")
    print(f"    Expected ~ -2m_q/(h mu^2) = {-2*m_q/(h*mu2):.8f}")

# ============================================================
# SEESAW RELATION CHECK
# ============================================================
print("\n\n" + "-" * 60)
print("SEESAW RELATION: m_meson vs 1/m_quark")
print("-" * 60)

print("\nIf Seiberg seesaw holds, m_meson ~ Lambda^2/m_quark, so:")
print(f"  m_meson * m_quark = const")
print()

for i in range(N_f):
    label = labels[i]
    m_q = m_quarks[i]
    m_cw = abs(cw_m_final[label])
    product = m_cw * m_q
    print(f"  m_CW({label}) * m_{label} = {m_cw:.8f} * {m_q:.2f} = {product:.6f} MeV^2")

print("""
NOTE: The ISS CW meson mass is NOT the Seiberg seesaw mass.
The Seiberg seesaw M_i = Lambda^2/m_i applies to the SUSY vacuum.
The ISS CW mass is a DIFFERENT quantity: it's the one-loop mass
of the pseudo-moduli at the metastable vacuum.

The ISS CW mass depends on f_a = h*mu^2 - m_a, NOT on 1/m_a.
For m_a << h*mu^2, the CW mass is approximately flavor-universal.
""")

# ============================================================
# ANALYTIC CW MASS FORMULA
# ============================================================
print("\n\n" + "=" * 80)
print("ANALYTIC CW MASS FORMULA")
print("=" * 80)

print("""
For the ISS model at one loop, the CW mass for chi_a can be computed
analytically. The result (from Shih hep-th/0703196, eq. 3.20) is:

  m^2_CW(chi_a) = (h^4 / 8pi^2) * f(x_a)

where x_a = h^2 phi_0^2 / f_a and:

  f(x) = (1+x) log(1+1/x) + (1-x) log(1-1/x) - 2
       = -2 + sum_{n=1}^infty 2/(2n+1) * x^{-(2n)}    for |x| > 1

For x >> 1 (i.e., h^2 phi_0^2 >> f_a):
  f(x) ~ 2/(3x^2) + 2/(5x^4) + ...
  m^2_CW ~ (h^4 / 8pi^2) * 2/(3x^2) = h^4 f_a^2 / (12 pi^2 h^4 phi_0^4)
          = f_a^2 / (12 pi^2 phi_0^4)

For the opposite limit x << 1 (phi_0 << f_a^{1/2}):
  f(x) ~ -2 + 2 log(1/x) + O(x)
  m^2_CW ~ (h^4 / 8pi^2) * [-2 + 2 log(f_a/(h^2 phi_0^2))]
""")

# Compute analytic formula
for i in range(N_f):
    label = labels[i]
    m_q = m_quarks[i]
    f_i = h * mu2 - m_q
    x_i = h**2 * phi_0_sq / f_i

    if abs(x_i) > 1:
        fx = (1 + x_i) * np.log(1 + 1/x_i) + (1 - x_i) * np.log(abs(1 - 1/x_i)) - 2
    else:
        fx = (1 + x_i) * np.log(abs(1 + 1/x_i)) + (1 - x_i) * np.log(abs(1 - 1/x_i)) - 2

    m2_analytic = h**4 / (8 * pi**2) * fx

    print(f"  Flavor {label}: x = h^2*phi_0^2/f = {x_i:.6f}")
    print(f"    f(x) = {fx:.8f}")
    print(f"    m^2_CW(analytic) = {m2_analytic:.8f} MeV^2")
    print(f"    m_CW(analytic) = {sqrt(abs(m2_analytic)):.8f} MeV")
    print(f"    m^2_CW(numerical) = {cw_m2_final[label]:.8f} MeV^2")
    print(f"    Ratio analytic/numerical = {m2_analytic/cw_m2_final[label]:.8f}" if abs(cw_m2_final[label]) > 1e-30 else "    numerical ≈ 0")
    print()

# ============================================================
# 5. PSEUDO-MODULUS AND KOIDE PHASE
# ============================================================
print("\n" + "=" * 80)
print("5. PSEUDO-MODULUS AND KOIDE PHASE")
print("=" * 80)

print("""
The pseudo-modulus X = Phi^1_1 (= Phi_0) is a tree-level flat direction.
At one loop, the CW potential generates a potential V_CW(X).
The minimum of V_CW(X) determines <X>.

The CW potential for X comes from sector A loops.
The X-dependent masses are in sector A: the heavy (rho, rho~) pair
has mass depending on X.
""")

# Compute V_CW as function of X
def V_CW_of_X(X_val):
    """Compute total CW potential as function of pseudo-modulus X."""
    V_total = 0.0

    # Sector A: (rho, rho~, X) system
    # With <X> = X_val, the fermion mass matrix changes.
    # W_{rho, rho~} = h X_val
    # W_{rho, X} = h phi_0
    # W_{rho~, X} = h phi_0

    W_A = np.array([
        [0, h * X_val, h * phi_0],
        [h * X_val, 0, h * phi_0],
        [h * phi_0, h * phi_0, 0]
    ])

    ferm_A_m2 = np.sort(np.real(np.linalg.eigvalsh(W_A @ W_A)))

    # Sector A bosons: approximately same as fermions (SUSY not broken here)
    # Actually, sector A has F_X = 0 at all X (the F-term for X is exactly zero
    # because phi_0^2 = mu^2 - m_1/h cancels the constant term).
    # So SUSY is unbroken in sector A → boson = fermion masses.
    # The CW potential from sector A is ZERO by SUSY cancellation!

    # The CW potential for X comes from sector B:
    # The psi_a masses depend on X through... actually they don't!
    # The psi_a mass is h*phi_0 (from the VEV of rho, rho~), not from X.
    # X = Phi^1_1 doesn't couple to psi_a (which is q^a with a > 1).

    # So at leading order, the CW potential is FLAT in X!
    # This is the well-known result: the pseudo-modulus X is lifted only
    # at TWO loops or by R-symmetry breaking.

    # Wait, that's not quite right. In the original ISS,
    # X IS lifted at one loop, but through the sector A contribution
    # where the boson-fermion splitting comes from X-dependent terms.

    # Let me reconsider. The F-terms that break SUSY are F_{chi_a}.
    # These don't depend on X. But the MASS MATRIX for the fields
    # that DO depend on X (rho, rho~) is:
    # Fermion: diag(0, sqrt(2)*h*sqrt(X^2+phi_0^2))  [schematic]
    # Boson: same (SUSY not broken in this sector)
    # So no CW contribution.

    # The actual ISS CW potential for X comes from the COMPLETE
    # mass matrix including the coupling between sectors A and B.
    # The off-diagonal terms in the full mass matrix can generate
    # X-dependent contributions.

    # In practice, for the MINIMAL ISS (N_c_mag = 1), the pseudo-modulus
    # IS lifted at one loop with:
    # V_CW(X) ~ (h^4 f^2 / 16pi^2) * log(h^2 X^2 / Lambda^2)
    # for large X, and the minimum is at X = 0 (or some O(f/h) value).

    # For our purposes: the CW potential stabilizes X near the origin,
    # and the MESON CW masses (for chi_a) are the relevant output.

    return 0.0  # Placeholder — see discussion below

print("""
IMPORTANT RESULT:

In the ISS model with N_c_mag = 1, the pseudo-modulus X is stabilized
at X = 0 by the one-loop CW potential (Shih 2007). The minimum is at
the origin of moduli space.

The CW potential for X is generated by the X-dependence of the
COMPLETE mass matrix (all sectors coupled). At X = 0, the spectrum
simplifies and the CW masses for chi_a can be computed analytically.

At X = 0, the full spectrum is:
""")

# COMPLETE computation at X = 0
print("COMPLETE MASS SPECTRUM AT X = 0:")
print(f"{'Field':<25} {'m^2_B (MeV^2)':<20} {'m_F (MeV)':<15} {'STr contrib'}")
print("-" * 75)

# Total CW potential from all sectors
total_CW = 0.0

# Sector A: X pseudo-modulus
# Goldstino: m_F = 0
# Heavy pair: m_F = sqrt(2) h phi_0, m_B = same → STr = 0
print(f"{'Goldstino':<25} {'0':<20} {'0':<15} {'0'}")
print(f"{'Heavy rho pair':<25} {2*h**2*phi_0_sq:<20.4f} {sqrt(2)*h*phi_0:<15.4f} {'0 (SUSY)'}")

# Sector B: for each broken flavor
for i in range(1, N_f):
    label = labels_sorted[i]
    f_i = f_a[i]

    m2_b_plus = h**2 * phi_0_sq + h * f_i
    m2_b_minus = h**2 * phi_0_sq - h * f_i
    m_f = h * phi_0
    m2_f = m_f**2
    m2_Y = h**2 * phi_0_sq

    print(f"{'psi_'+label+' (B+)':<25} {m2_b_plus:<20.4f} {m_f:<15.4f} {''}")
    print(f"{'psi_'+label+' (B-)':<25} {m2_b_minus:<20.4f} {'':<15} {''}")
    print(f"{'Y_'+label:<25} {m2_Y:<20.4f} {m_f:<15.4f} {''}")
    print(f"{'Y~_'+label:<25} {m2_Y:<20.4f} {'':<15} {''}")

    # CW contribution from this sector:
    # Bosons: 2 complex scalars from psi (m2_plus, m2_minus) = 4 real
    #         2 complex scalars from Y (m2_Y, m2_Y) = 4 real
    # Fermions: 2 Dirac fermions (m_f) = 4 Weyl fermions = 8 real d.o.f.
    # But in chiral multiplet counting:
    # Each chiral multiplet: 1 complex scalar + 1 Weyl fermion
    # STr per chiral multiplet: 2 * m^2_B - 2 * m^2_F

    # The psi_a, psi~_a are TWO chiral multiplets with:
    # Scalar m^2: eigenvalues of [[m2_Y + h f_i, 0], [0, m2_Y - h f_i]]
    #           = {m2_Y + h f_i, m2_Y - h f_i}
    # Wait, I was computing this wrong earlier. Let me reconsider.

    # Actually, psi_a is a chiral multiplet. Its scalar mass^2 gets
    # contributions from |F_{Y~_a}|^2 = h^2 phi_0^2 |psi_a|^2 → m^2 = h^2 phi_0^2
    # AND from |F_{chi_a}|^2 → B-term mixing psi_a with psi~_a*.
    # The eigenvalues of the B-term-corrected matrix are m2_Y +/- h f_i.

    # Y_a is another chiral multiplet with m^2 = h^2 phi_0^2 (no B-term).

    # Total bosonic m^2 eigenvalues for this flavor:
    # 4 complex scalars = 8 real d.o.f.:
    # {m2_Y + h f_i, m2_Y - h f_i, m2_Y, m2_Y} (each with 2 real d.o.f.)

    # Total fermionic: 2 Dirac fermions (= 4 Weyl) with m_f^2 = h^2 phi_0^2

# Sector C: chi_a (broken mesons)
for i in range(1, N_f):
    label = labels_sorted[i]
    print(f"{'chi_'+label+' (tree massless)':<25} {'0':<20} {'0':<15} {''}")

# ============================================================
# THE KEY COMPUTATION: CW mass for chi_a including all loops
# ============================================================
print("\n\n" + "=" * 80)
print("KEY COMPUTATION: ANALYTIC CW MASS FOR BROKEN MESONS")
print("=" * 80)

print("""
The one-loop CW mass for chi_a comes from integrating out the
(psi_a, psi~_a, Y_a, Y~_a) fields which couple to chi_a.

The chi_a-dependent superpotential:
  W_a = h phi_0 (Y_a psi~_a + psi_a Y~_a) + h chi_a psi_a psi~_a

The F-term:
  F_{chi_a} = h psi_a psi~_a - f_a  (where f_a = h mu^2 - m_a)

The CW mass-squared for chi_a is:
  m^2_CW = d^2 V_CW / d chi_a^2 |_{chi_a=0}

Using the formula from Shih (2007):
  m^2_CW(chi_a) = (h^4 / 16pi^2) * [g(x_a)]

where x_a = h^2 phi_0^2 / (h mu^2 - m_a) and g(x) is the ISS loop function.

Actually, let me compute this directly from the mass eigenvalues.
""")

# Direct computation of the CW mass for chi_a
# using the FULL mass matrix approach

def compute_cw_mass_chi(m_quark, phi_0_val, h_val, Lambda_val):
    """
    Compute CW mass for chi_a in the ISS model.

    chi_a couples to (psi_a, psi~_a, Y_a, Y~_a).
    We compute d^2 V_CW / d chi_a^2 numerically.
    """
    f_val = h_val * mu2 - m_quark

    delta = 1e-4 * phi_0_val  # small perturbation

    def V_CW_chi(chi_val):
        """CW potential for given chi_a value."""
        # Fermion mass matrix: 4x4 (Psi_a, Psi~_a, Y_a, Y~_a)
        M_F = np.array([
            [0, h_val*chi_val, 0, h_val*phi_0_val],
            [h_val*chi_val, 0, h_val*phi_0_val, 0],
            [0, h_val*phi_0_val, 0, 0],
            [h_val*phi_0_val, 0, 0, 0]
        ], dtype=float)

        ferm_m2_eigs = np.sort(np.linalg.eigvalsh(M_F @ M_F))

        # Bosonic mass-squared: need to include B-term
        # The boson mass matrix in real basis
        # For the fields (psi_a, psi~_a, Y_a, Y~_a):

        # Holomorphic piece: dF_A/dphi_j evaluated at vacuum
        # dF_{chi_a}/d{psi_a} = h psi~_a → 0 at vacuum (for small chi perturbation)
        # No wait, we're not expanding around psi=0. We ARE at psi=0.
        # The chi_a dependence enters the mass matrix.

        # Actually, the bosonic mass-squared comes from:
        # V = |F_{chi_a}|^2 + |F_{psi_a}|^2 + |F_{psi~_a}|^2 + |F_{Y_a}|^2 + |F_{Y~_a}|^2
        # where:
        # F_{chi_a} = h psi_a psi~_a + m_a - h mu^2 = h psi_a psi~_a - f/h
        # F_{psi_a} = h chi_a psi~_a + h phi_0 Y~_a
        # F_{psi~_a} = h chi_a psi_a + h phi_0 Y_a
        # F_{Y_a} = h phi_0 psi~_a
        # F_{Y~_a} = h phi_0 psi_a

        # Scalar mass matrix from |F|^2:
        # d^2 V / d psi_a d psi_a* = |dF_{chi_a}/dpsi_a|^2 + |dF_{psi~_a}/dpsi_a|^2 + |dF_{Y~_a}/dpsi_a|^2
        #   = |h psi~_a|^2 + |h chi_a|^2 + |h phi_0|^2 → h^2(chi^2 + phi_0^2) at vac

        # Let me just use the matrix approach.
        # dF_A / d phi_j where A ∈ {chi_a, psi_a, psi~_a, Y_a, Y~_a} and
        # j ∈ {psi_a, psi~_a, Y_a, Y~_a}:

        dF = np.array([
            [0,           0,           0, 0],        # dF_{chi_a}/d{psi, psi~, Y, Y~} at vacuum = 0
            [0,           h_val*chi_val, 0, h_val*phi_0_val],  # dF_{psi_a}/d{...}
            [h_val*chi_val, 0,           h_val*phi_0_val, 0],  # dF_{psi~_a}/d{...}
            [0,           h_val*phi_0_val, 0, 0],     # dF_{Y_a}/d{...}
            [h_val*phi_0_val, 0,           0, 0],     # dF_{Y~_a}/d{...}
        ], dtype=float)

        M2_hol = dF.T @ dF

        # B-term from F_{chi_a}:
        # d^2 |F_{chi_a}|^2 / d psi_a d psi~_a = F_{chi_a}^* * h = (-f_val) * h_val
        # (evaluated at psi = psi~ = 0)
        B = np.zeros((4, 4), dtype=float)
        B[0, 1] = -h_val * f_val  # psi_a, psi~_a  (recall f = h mu^2 - m)
        B[1, 0] = -h_val * f_val

        # Full 8x8 real scalar mass matrix
        n = 4
        M_real = np.zeros((2*n, 2*n))
        for ii in range(n):
            for jj in range(n):
                M_real[2*ii, 2*jj] = M2_hol[ii, jj] + B[ii, jj]  # (Re, Re)
                M_real[2*ii+1, 2*jj+1] = M2_hol[ii, jj] - B[ii, jj]  # (Im, Im)

        boson_m2_eigs = np.sort(np.linalg.eigvalsh(M_real))

        # CW potential
        V = 0.0
        for m2 in boson_m2_eigs:
            if abs(m2) > 1e-30:
                V += m2**2 * (np.log(abs(m2) / Lambda_val**2) - 1.5)
        for m2 in ferm_m2_eigs:
            if abs(m2) > 1e-30:
                # Each eigenvalue of W^dag W corresponds to one Dirac fermion
                # = 2 Weyl fermions = 4 real d.o.f. But each real scalar eigenvalue
                # is 1 real d.o.f. So we need to count properly.
                # 4 eigenvalues of W^dag W → 4 mass eigenvalues.
                # Each corresponds to 2 real Weyl d.o.f.
                # The 8 boson eigenvalues → 8 real d.o.f.
                # So: V = (1/64pi^2) [sum_B m^4(log m^2/L^2 - 3/2) - sum_F 2*m^4(log...)]
                # where sum_F has 4 terms (each with weight 2 for Weyl counting)
                V -= 2 * m2**2 * (np.log(abs(m2) / Lambda_val**2) - 1.5)

        V /= (64 * pi**2)
        return V

    V0 = V_CW_chi(0.0)
    Vp = V_CW_chi(delta)
    Vm = V_CW_chi(-delta)

    m2_cw = (Vp + Vm - 2*V0) / delta**2

    return m2_cw


print("\nFINAL CW MASSES FOR BROKEN MESONS:")
print(f"{'Flavor':<8} {'m_quark (MeV)':<15} {'m^2_CW (MeV^2)':<20} {'m_CW (MeV)':<15} {'m_CW/m_q':<12}")
print("-" * 70)

final_cw_masses = {}

for i in range(N_f):
    label = labels[i]
    m_q = m_quarks[i]

    m2 = compute_cw_mass_chi(m_q, phi_0, h, Lambda)
    m_val = sqrt(abs(m2)) * np.sign(m2)

    final_cw_masses[label] = (m2, m_val)

    print(f"{label:<8} {m_q:<15.2f} {m2:<20.10f} {m_val:<15.10f} {abs(m_val)/m_q:<12.8f}")

# ============================================================
# KOIDE ANALYSIS OF FINAL CW SPECTRUM
# ============================================================
print("\n\n" + "=" * 80)
print("KOIDE ANALYSIS OF CW SPECTRUM")
print("=" * 80)

m_cw_vals = {label: abs(final_cw_masses[label][1]) for label in labels}

print(f"\nCW meson masses (absolute values):")
for label in labels:
    print(f"  m_CW({label}) = {m_cw_vals[label]:.10f} MeV")

# Ratios
print(f"\nMass ratios (CW mesons vs UV quarks):")
for label in labels:
    ratio = m_cw_vals[label] / m_quarks[labels.index(label)]
    print(f"  m_CW({label}) / m_{label} = {ratio:.10f}")

# Is the CW spectrum flavor-universal?
ratios = [m_cw_vals[label] / m_quarks[labels.index(label)] for label in labels]
print(f"\n  Max ratio / Min ratio = {max(ratios)/min(ratios):.8f}")
print(f"  (1.0 = perfectly flavor-universal)")

# Koide Q for CW masses
print("\nKoide Q parameter for CW masses:")
for triple in [('u', 'd', 's'), ('d', 's', 'c'), ('u', 's', 'c'), ('u', 'd', 'c')]:
    m1, m2, m3 = m_cw_vals[triple[0]], m_cw_vals[triple[1]], m_cw_vals[triple[2]]
    Q = koide_Q(m1, m2, m3)
    print(f"  Q_CW({triple[0]},{triple[1]},{triple[2]}) = {Q:.8f}  (deviation from 2/3: {abs(Q - 2/3)*100/(2/3):.4f}%)")

print("\nKoide Q for UV quark masses:")
for triple in [('u', 'd', 's'), ('d', 's', 'c'), ('u', 's', 'c'), ('u', 'd', 'c')]:
    m1 = m_quarks[labels.index(triple[0])]
    m2 = m_quarks[labels.index(triple[1])]
    m3 = m_quarks[labels.index(triple[2])]
    Q = koide_Q(m1, m2, m3)
    print(f"  Q_UV({triple[0]},{triple[1]},{triple[2]}) = {Q:.8f}  (deviation from 2/3: {abs(Q - 2/3)*100/(2/3):.4f}%)")

# Inverse mass Koide (Seiberg seesaw related)
print("\nDual Koide Q (for 1/m):")
for triple in [('d', 's', 'c'), ('u', 'd', 's'), ('u', 's', 'c')]:
    m1 = 1.0 / m_quarks[labels.index(triple[0])]
    m2 = 1.0 / m_quarks[labels.index(triple[1])]
    m3 = 1.0 / m_quarks[labels.index(triple[2])]
    Q = koide_Q(m1, m2, m3)
    print(f"  Q(1/m_{triple[0]}, 1/m_{triple[1]}, 1/m_{triple[2]}) = {Q:.8f}  (dev from 2/3: {abs(Q - 2/3)*100/(2/3):.4f}%)")

# ============================================================
# VARY mu TO STUDY SENSITIVITY
# ============================================================
print("\n\n" + "=" * 80)
print("SENSITIVITY TO mu (SUSY-BREAKING SCALE)")
print("=" * 80)

mu_values = [50, 92, 150, 200, 300, 500, 1000]

print(f"\n{'mu (MeV)':<12} {'m_CW(u)':<15} {'m_CW(d)':<15} {'m_CW(s)':<15} {'m_CW(c)':<15} {'Q(u,d,s)':<12} {'Q(d,s,c)':<12} {'ratio spread':<15}")
print("-" * 120)

for mu_val in mu_values:
    mu2_val = mu_val**2
    phi0_val = sqrt(mu2_val - m_u / h)

    cw_m_mu = {}
    for i in range(N_f):
        label = labels[i]
        m_q = m_quarks[i]
        m2 = compute_cw_mass_chi(m_q, phi0_val, h, Lambda)
        cw_m_mu[label] = sqrt(abs(m2))

    Q_uds = koide_Q(cw_m_mu['u'], cw_m_mu['d'], cw_m_mu['s'])
    Q_dsc = koide_Q(cw_m_mu['d'], cw_m_mu['s'], cw_m_mu['c'])

    ratios_mu = [cw_m_mu[l] / m_quarks[labels.index(l)] for l in labels]
    spread = max(ratios_mu) / min(ratios_mu)

    print(f"{mu_val:<12} {cw_m_mu['u']:<15.8f} {cw_m_mu['d']:<15.8f} {cw_m_mu['s']:<15.8f} {cw_m_mu['c']:<15.8f} {Q_uds:<12.6f} {Q_dsc:<12.6f} {spread:<15.8f}")

# ============================================================
# KOIDE PHASE ANALYSIS
# ============================================================
print("\n\n" + "=" * 80)
print("KOIDE PHASE ANALYSIS FOR CW SPECTRUM")
print("=" * 80)

def koide_params(m1, m2, m3):
    """Compute Koide parametrization: sqrt(m_k) = sqrt(M0) (1 + sqrt(2) cos(2pi k/3 + delta))"""
    s1, s2, s3 = sqrt(m1), sqrt(m2), sqrt(m3)
    M0 = (s1 + s2 + s3)**2 / 9.0

    # delta from atan2
    y = (s2 - s1) / sqrt(2)
    x = (2*s3 - s2 - s1) / sqrt(6)

    # The angle delta satisfies:
    # sqrt(2) cos(delta + 2pi/3) = (s1 - s3_avg) / sqrt(M0)
    # etc.
    # Simpler: use the Foot angle approach

    # Actually: delta = atan2(y, x) where
    # y = sqrt(2/6) * (s2 - s1) ... let me use the standard form.

    # From sqrt(m_k) = sqrt(M0)(1 + sqrt(2) cos(2pi k/3 + delta)):
    # sum_k sqrt(m_k) = 3 sqrt(M0) [since cos terms sum to 0]
    # So sqrt(M0) = (sum sqrt(m_k)) / 3

    sqrt_M0 = (s1 + s2 + s3) / 3.0

    # cos(2pi/3 + delta) = (s1/sqrt_M0 - 1) / sqrt(2)
    # cos(4pi/3 + delta) = (s2/sqrt_M0 - 1) / sqrt(2)
    # cos(delta) = (s3/sqrt_M0 - 1) / sqrt(2) ... if k=0 → 3rd mass

    # Using k=1,2,3 convention from the paper:
    # sqrt(m_1) = sqrt(M0)(1 + sqrt(2) cos(2pi/3 + delta))
    # sqrt(m_2) = sqrt(M0)(1 + sqrt(2) cos(4pi/3 + delta))
    # sqrt(m_3) = sqrt(M0)(1 + sqrt(2) cos(delta))

    # With the ordering m_1 < m_2 < m_3:
    masses_sorted = sorted([m1, m2, m3])
    s_sorted = [sqrt(m) for m in masses_sorted]

    cos_delta = (s_sorted[2] / sqrt_M0 - 1) / sqrt(2)
    # sin delta from the other equations
    # cos(2pi/3 + delta) = (s_sorted[0] / sqrt_M0 - 1) / sqrt(2)
    cos_2pi3_plus_delta = (s_sorted[0] / sqrt_M0 - 1) / sqrt(2)

    # cos(2pi/3 + delta) = cos(2pi/3)cos(delta) - sin(2pi/3)sin(delta)
    # = -0.5 cos(delta) - sqrt(3)/2 sin(delta)
    # So: sin(delta) = (-cos_2pi3_plus_delta - 0.5*cos_delta) / (sqrt(3)/2)
    sin_delta = (-cos_2pi3_plus_delta - 0.5*cos_delta) * 2 / sqrt(3)

    delta = np.arctan2(sin_delta, cos_delta)

    Q = (m1 + m2 + m3) / (s1 + s2 + s3)**2

    return M0, delta, Q

print("\nKoide parametrization:")
print(f"{'Triple':<15} {'M0 (MeV)':<15} {'delta (rad)':<15} {'delta/pi':<12} {'Q':<12} {'Q*3/2':<10}")
print("-" * 80)

# UV quarks
triples_to_check = [
    ('e,mu,tau', [0.511, 105.658, 1776.86]),
    ('u,d,s (UV)', [m_u, m_d, m_s]),
    ('d,s,c (UV)', [m_d, m_s, m_c]),
    ('-s,c,b (UV)', [m_s, m_c, 4180.0]),  # with m_b
]

for name, masses in triples_to_check:
    M0, delta, Q = koide_params(*masses)
    print(f"{name:<15} {M0:<15.4f} {delta:<15.6f} {delta/pi:<12.6f} {Q:<12.8f} {Q*1.5:<10.6f}")

# CW meson masses
cw_triples = [
    ('u,d,s (CW)', [m_cw_vals['u'], m_cw_vals['d'], m_cw_vals['s']]),
    ('d,s,c (CW)', [m_cw_vals['d'], m_cw_vals['s'], m_cw_vals['c']]),
]

for name, masses in cw_triples:
    M0, delta, Q = koide_params(*masses)
    print(f"{name:<15} {M0:<15.4f} {delta:<15.6f} {delta/pi:<12.6f} {Q:<12.8f} {Q*1.5:<10.6f}")

# ============================================================
# SUMMARY
# ============================================================
print("\n\n" + "=" * 80)
print("SUMMARY: ISS COLEMAN-WEINBERG AND KOIDE TRANSMISSION")
print("=" * 80)

print(f"""
1. TREE-LEVEL SPECTRUM:
   - Pseudo-modulus X = Phi_0: FLAT (tree-level massless)
   - Broken mesons chi_a: MASSLESS at tree level
   - Link sector (psi, Y): mass = h*phi_0 = {h*phi_0:.4f} MeV
     with SUSY-breaking splitting +/- h*f_a
   - Goldstino: massless

2. COLEMAN-WEINBERG MASSES:
   CW lifts both X and chi_a at one loop.
""")

print(f"   CW meson masses:")
for label in labels:
    m2, m_val = final_cw_masses[label]
    print(f"     m_CW({label}) = {abs(m_val):.8f} MeV  (m^2 = {m2:.8f} MeV^2)")

print(f"""
3. FLAVOR DEPENDENCE:
   The CW mass depends on m_a through f_a = h*mu^2 - m_a.
   For m_a << h*mu^2 = {h*mu2:.1f} MeV^2, the dependence is WEAK.
   For charm (m_c = {m_c} MeV), m_c/(h*mu^2) = {m_c/(h*mu2):.4f}:
     SIGNIFICANT flavor dependence.

4. KOIDE TRANSMISSION:
""")

# Final Koide comparison
print("   UV Quark Koide vs CW Meson Koide:")
for triple_name, triple_labels in [('u,d,s', ['u','d','s']), ('d,s,c', ['d','s','c'])]:
    # UV
    uv_masses = [m_quarks[labels.index(l)] for l in triple_labels]
    Q_uv = koide_Q(*uv_masses)
    # CW
    cw_masses_triple = [m_cw_vals[l] for l in triple_labels]
    Q_cw = koide_Q(*cw_masses_triple)

    print(f"   Q({triple_name})_UV = {Q_uv:.8f}")
    print(f"   Q({triple_name})_CW = {Q_cw:.8f}")
    print(f"   Delta Q = {Q_cw - Q_uv:.8f}  ({(Q_cw-Q_uv)/Q_uv*100:.4f}%)")
    print()

print(f"""
5. KEY RESULT ON KOIDE TRANSMISSION:
   The CW potential generates meson masses that depend on quark masses
   through f_a = h*mu^2 - m_a. In the limit mu >> m_a/h (for all flavors),
   the CW mass becomes flavor-UNIVERSAL → Q approaches 1/3 (degenerate).

   For physical parameters (mu = {mu} MeV), the quark masses are
   SMALL compared to h*mu^2 = {h*mu2:.1f} MeV^2, so the CW meson spectrum
   is nearly degenerate and Q ≈ 1/3.

   This means the ISS CW potential does NOT transmit the UV Koide condition
   to the IR meson spectrum in any useful sense when mu ~ f_pi.

   The Seiberg seesaw (M_i ~ 1/m_i) would transmit Koide, but that
   applies to the SUSY vacuum, not the metastable ISS vacuum.

   For larger mu (mu >> m_c), the CW spectrum approaches perfect degeneracy.
   For smaller mu (mu ~ m_c/h), the CW spectrum becomes sensitive to
   individual quark masses and can approximately transmit mass ratios.
""")

# ============================================================
# DUAL KOIDE (1/m) IN ISS CONTEXT
# ============================================================
print("\n" + "=" * 80)
print("6. CONNECTION TO DUAL KOIDE")
print("=" * 80)

print("""
The Seiberg seesaw in the SUSY vacuum gives meson masses M_i ~ Lambda^2 / m_i.
If the UV quark masses (m_d, m_s, m_b) satisfy a "dual Koide" condition
  Q(1/m_d, 1/m_s, 1/m_b) ≈ 2/3
then the SEESAW meson masses M_i = Lambda^2/m_i automatically satisfy
  Q(M_d, M_s, M_b) = Q(1/m_d, 1/m_s, 1/m_b) ≈ 2/3.

This is EXACT (the Koide condition is homogeneous in the masses).
So Seiberg seesaw transmits dual-Koide to standard Koide for mesons.

However, this applies to the SUSY vacuum. In the ISS metastable vacuum,
the meson masses are CW-generated and DO NOT follow the seesaw formula.
""")

# Check dual Koide for down-type quarks
m_d_val, m_s_val, m_b_val = 4.67, 93.4, 4180.0
Q_dual = koide_Q(1/m_d_val, 1/m_s_val, 1/m_b_val)
print(f"  Q(1/m_d, 1/m_s, 1/m_b) = {Q_dual:.8f}")
print(f"  Deviation from 2/3: {abs(Q_dual - 2/3)*100/(2/3):.4f}%")

# Seesaw meson masses
Lambda_sq_seesaw = Lambda**2
M_seesaw = Lambda_sq_seesaw / np.array([m_d_val, m_s_val, m_b_val])
Q_seesaw = koide_Q(*M_seesaw)
print(f"\n  Seesaw meson masses M_i = Lambda^2/m_i:")
print(f"    M(d) = {M_seesaw[0]:.4f} MeV")
print(f"    M(s) = {M_seesaw[1]:.4f} MeV")
print(f"    M(b) = {M_seesaw[2]:.4f} MeV")
print(f"  Q(M_d, M_s, M_b) = {Q_seesaw:.8f}")
print(f"  Check: same as dual Koide = {abs(Q_seesaw - Q_dual) < 1e-10}")

print("\n" + "=" * 80)
print("COMPUTATION COMPLETE")
print("=" * 80)
