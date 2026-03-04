"""
Seiberg effective theory spectrum for SU(3)_c SQCD.

Part 1: N_f = N_c = 3 (confining, quantum-modified moduli space)
Part 2: ISS regime N_f = 4, N_c = 3

Computes:
- Meson VEVs M_j (Seiberg seesaw)
- Fermion mass matrix d^2W/dPhi_I dPhi_J
- Scalar mass-squared matrix
- Supertrace STr[M^2]
- Koide parameters Q for various triples
"""

import numpy as np
from itertools import combinations

# ============================================================
# Physical inputs
# ============================================================
Lambda = 300.0  # MeV
m_u = 2.16      # MeV
m_d = 4.67      # MeV
m_s = 93.4      # MeV
m_c = 1270.0    # MeV
m_b = 4180.0    # MeV

masses_all = {'u': m_u, 'd': m_d, 's': m_s, 'c': m_c, 'b': m_b}

def koide_Q(m1, m2, m3):
    """Koide ratio Q = (m1+m2+m3)^2 / (3*(m1^2+m2^2+m3^2)) using sqrt masses"""
    # Actually standard Koide: Q = (sum sqrt(m))^2 / (3 * sum m)
    # But let's use the mass version: Q = (sum m)^2 / (3 sum m^2)
    # No wait -- the standard Koide formula is:
    # Q = (sqrt(m1) + sqrt(m2) + sqrt(m3))^2 / (3(m1 + m2 + m3))
    s1 = np.sqrt(m1) + np.sqrt(m2) + np.sqrt(m3)
    s2 = m1 + m2 + m3
    return s1**2 / (3.0 * s2)

def koide_Q_direct(m1, m2, m3):
    """Koide Q using the masses directly (not sqrt)."""
    return koide_Q(m1, m2, m3)

print("=" * 80)
print("PART 1: N_f = N_c = 3 -- CONFINING PHASE (QUANTUM-MODIFIED MODULI SPACE)")
print("=" * 80)

# ============================================================
# 1A. Seiberg seesaw for (u, d, s)
# ============================================================
print("\n--- 1A. Seiberg Seesaw: (u, d, s) ---\n")

flavors_uds = ['u', 'd', 's']
m_uds = np.array([m_u, m_d, m_s])

# M_j = Lambda^2 * (prod m_k)^{1/3} / m_j
# C = Lambda^2 * (m_1 m_2 m_3)^{1/3}
C_uds = Lambda**2 * np.prod(m_uds)**(1.0/3.0)
M_uds = C_uds / m_uds

print(f"C = Lambda^2 * (m_u * m_d * m_s)^(1/3) = {C_uds:.4f} MeV^2")
print(f"  = {Lambda}^2 * ({np.prod(m_uds):.4f})^(1/3) = {Lambda**2:.0f} * {np.prod(m_uds)**(1./3.):.4f}")
print()

# X (Lagrange multiplier)
X_uds = np.prod(m_uds)**(1.0/3.0) / Lambda**2
print(f"X = (m_u m_d m_s)^(1/3) / Lambda^2 = {X_uds:.6e} MeV^(-1)")
print()

for i, f in enumerate(flavors_uds):
    print(f"M_{f} = {M_uds[i]:.4f} MeV")

print(f"\nVerification: det M = M_u * M_d * M_s = {np.prod(M_uds):.4f}")
print(f"Lambda^6 = {Lambda**6:.4f}")
print(f"Match: {np.isclose(np.prod(M_uds), Lambda**6)}")

print(f"\nMeson VEV ratios:")
print(f"M_u / M_d = {M_uds[0]/M_uds[1]:.4f} = m_d/m_u = {m_d/m_u:.4f}")
print(f"M_u / M_s = {M_uds[0]/M_uds[2]:.4f} = m_s/m_u = {m_s/m_u:.4f}")
print(f"M_d / M_s = {M_uds[1]/M_uds[2]:.4f} = m_s/m_d = {m_s/m_d:.4f}")

print(f"\nm_j * M_j products (should all be equal = C):")
for i, f in enumerate(flavors_uds):
    print(f"  m_{f} * M_{f} = {m_uds[i] * M_uds[i]:.4f} MeV^2")

# Koide checks
print(f"\nKoide Q(m_u, m_d, m_s) = {koide_Q(m_u, m_d, m_s):.6f}")
print(f"Koide Q(M_u, M_d, M_s) = {koide_Q(M_uds[0], M_uds[1], M_uds[2]):.6f}")
print(f"Koide Q(1/m_u, 1/m_d, 1/m_s) = {koide_Q(1/m_u, 1/m_d, 1/m_s):.6f}")

# Check: Q(M_j) should equal Q(1/m_j) since M_j = C/m_j and Koide is scale-invariant
print(f"\nQ(M) = Q(1/m) [scale invariance]: {np.isclose(koide_Q(M_uds[0], M_uds[1], M_uds[2]), koide_Q(1/m_u, 1/m_d, 1/m_s))}")

# ============================================================
# 1B. Seiberg seesaw for (s, c, b)
# ============================================================
print("\n--- 1B. Seiberg Seesaw: (s, c, b) ---\n")

flavors_scb = ['s', 'c', 'b']
m_scb = np.array([m_s, m_c, m_b])

C_scb = Lambda**2 * np.prod(m_scb)**(1.0/3.0)
M_scb = C_scb / m_scb

print(f"C = Lambda^2 * (m_s * m_c * m_b)^(1/3) = {C_scb:.4f} MeV^2")

for i, f in enumerate(flavors_scb):
    print(f"M_{f} = {M_scb[i]:.4f} MeV")

print(f"\nVerification: det M = {np.prod(M_scb):.4f}, Lambda^6 = {Lambda**6:.4f}")
print(f"Match: {np.isclose(np.prod(M_scb), Lambda**6)}")

print(f"\nm_j * M_j products:")
for i, f in enumerate(flavors_scb):
    print(f"  m_{f} * M_{f} = {m_scb[i] * M_scb[i]:.4f} MeV^2")

print(f"\nKoide Q(m_s, m_c, m_b) = {koide_Q(m_s, m_c, m_b):.6f}")
print(f"Koide Q(M_s, M_c, M_b) = {koide_Q(M_scb[0], M_scb[1], M_scb[2]):.6f}")
print(f"Koide Q(1/m_s, 1/m_c, 1/m_b) = {koide_Q(1/m_s, 1/m_c, 1/m_b):.6f}")

# ============================================================
# 1C. Seiberg seesaw for (d, s, b) -- the dual Koide triple
# ============================================================
print("\n--- 1C. Seiberg Seesaw: (d, s, b) -- the dual Koide triple ---\n")

flavors_dsb = ['d', 's', 'b']
m_dsb = np.array([m_d, m_s, m_b])

C_dsb = Lambda**2 * np.prod(m_dsb)**(1.0/3.0)
M_dsb = C_dsb / m_dsb

print(f"C = Lambda^2 * (m_d * m_s * m_b)^(1/3) = {C_dsb:.4f} MeV^2")

for i, f in enumerate(flavors_dsb):
    print(f"M_{f} = {M_dsb[i]:.4f} MeV")

print(f"\nKoide Q(m_d, m_s, m_b) = {koide_Q(m_d, m_s, m_b):.6f}")
print(f"Koide Q(M_d, M_s, M_b) = {koide_Q(M_dsb[0], M_dsb[1], M_dsb[2]):.6f}")
print(f"Koide Q(1/m_d, 1/m_s, 1/m_b) = {koide_Q(1/m_d, 1/m_s, 1/m_b):.6f}")
print(f"\n** Q(1/m_d, 1/m_s, 1/m_b) = {koide_Q(1/m_d, 1/m_s, 1/m_b):.6f}")
print(f"   Deviation from 2/3: {abs(koide_Q(1/m_d, 1/m_s, 1/m_b) - 2./3.):.6f}")
print(f"   Fractional deviation: {abs(koide_Q(1/m_d, 1/m_s, 1/m_b) - 2./3.)/(2./3.)*100:.4f}%")

# ============================================================
# PART 2: FERMION MASS MATRIX
# ============================================================
print("\n" + "=" * 80)
print("PART 2: FERMION MASS MATRIX (N_f = N_c = 3)")
print("=" * 80)

# The fields are: M^i_j (9 components for N_f=3) and X (Lagrange multiplier)
# Total: 10 chiral superfields
#
# W = sum_j m_j M^j_j + X (det M - Lambda^6)
#
# At the vacuum: M diagonal with M_j = C/m_j, X = -(product m)^{1/3}/Lambda^2
# (negative sign from the F-term equation m_j + X * cofactor = 0)
#
# The fermion mass matrix is F_{IJ} = d^2W / dPhi_I dPhi_J
# evaluated at the vacuum.
#
# Fields: Phi_I = {M^1_1, M^2_2, M^3_3, M^1_2, M^2_1, M^1_3, M^3_1, M^2_3, M^3_2, X}
#
# First, we need derivatives of det M.
# For diagonal M: det M = M^1_1 * M^2_2 * M^3_3
#
# d(det M)/d(M^i_j) = cofactor(M)^j_i  (with proper index structure)
# For diagonal vacuum:
# d(det M)/d(M^j_j) = product of other diagonal M's = det M / M^j_j
# d(det M)/d(M^i_j) for i != j: only nonzero for epsilon-type terms
#
# Second derivatives d^2(det M)/d(M^a_b)d(M^c_d):
# For N_f = 3: det M = (1/6) epsilon^{ijk} epsilon_{lmn} M^l_i M^m_j M^n_k
# So d^2(det M)/d(M^a_b)d(M^c_d) = sum of terms with one remaining M

# Let me be very explicit. For 3x3 matrix:
# det M = M^1_1(M^2_2 M^3_3 - M^2_3 M^3_2) - M^1_2(M^2_1 M^3_3 - M^2_3 M^3_1) + M^1_3(M^2_1 M^3_2 - M^2_2 M^3_1)

# Let's define the fields in a flat index: M_{(i,j)} where (i,j) runs over all 9 entries
# and X as the 10th field.

# I'll work with numerical second derivatives of W.

Nf = 3
Nfields = Nf**2 + 1  # 9 meson + 1 multiplier = 10

# Use (u, d, s) for the concrete case
m_vec = np.array([m_u, m_d, m_s])

def det3(M):
    """Determinant of 3x3 matrix."""
    return np.linalg.det(M)

def superpotential_W(M_flat, X_val, m_vec, Lambda):
    """
    W = Tr(m M) + X(det M - Lambda^6)
    M_flat is the 9-component flattened M^i_j (row-major: M[i,j])
    """
    M = M_flat.reshape(3, 3)
    # Tr(m M) = sum_j m_j * M^j_j (only diagonal of m contributes)
    trm = np.sum(m_vec * np.diag(M))
    d = det3(M)
    return trm + X_val * (d - Lambda**6)

# Vacuum values
M_vac = np.diag(M_uds)  # diagonal meson VEV
X_vac = -m_vec[0] / (det3(M_vac) / M_uds[0])  # from F_{M^1_1} = 0
# Check: m_1 + X * (M_2 * M_3) = 0 => X = -m_1/(M_2*M_3)
X_check = -m_vec[0] / (M_uds[1] * M_uds[2])
print(f"\nX_vac = {X_vac:.6e}")
print(f"X_check = {X_check:.6e}")
print(f"X from formula = -(prod m)^(1/3)/Lambda^2 = {-np.prod(m_vec)**(1./3.)/Lambda**2:.6e}")

# Verify F-terms vanish
for j in range(3):
    cofM = det3(M_vac) / M_uds[j]  # cofactor for diagonal
    Fj = m_vec[j] + X_vac * cofM
    print(f"F_{{M^{j+1}_{j+1}}} = {Fj:.6e}")

FX = det3(M_vac) - Lambda**6
print(f"F_X = det M - Lambda^6 = {FX:.4e}")

# ============================================================
# Fermion mass matrix: d^2W/dPhi_I dPhi_J
# ============================================================
print("\n--- Fermion Mass Matrix ---\n")

# Fields indexed as:
# 0-8: M^i_j in row-major order: M^1_1, M^1_2, M^1_3, M^2_1, M^2_2, M^2_3, M^3_1, M^3_2, M^3_3
# 9: X

# Numerical second derivatives via finite differences
eps = 1e-4  # step size in MeV

def field_to_state(M_flat, X_val):
    """Pack into single state vector."""
    return np.concatenate([M_flat, [X_val]])

M_vac_flat = M_vac.flatten()  # row-major
state_vac = field_to_state(M_vac_flat, X_vac)

def W_from_state(state):
    """Evaluate W from state vector."""
    M_flat = state[:9]
    X_val = state[9]
    return superpotential_W(M_flat, X_val, m_vec, Lambda)

# Check W at vacuum
W_vac = W_from_state(state_vac)
print(f"W at vacuum = {W_vac:.6f} MeV^2")

# Compute d^2W/dPhi_I dPhi_J numerically
fermion_mass = np.zeros((Nfields, Nfields))

for I in range(Nfields):
    for J in range(I, Nfields):
        # d^2W / dPhi_I dPhi_J by central differences
        state_pp = state_vac.copy()
        state_pm = state_vac.copy()
        state_mp = state_vac.copy()
        state_mm = state_vac.copy()

        state_pp[I] += eps; state_pp[J] += eps
        state_pm[I] += eps; state_pm[J] -= eps
        state_mp[I] -= eps; state_mp[J] += eps
        state_mm[I] -= eps; state_mm[J] -= eps

        d2W = (W_from_state(state_pp) - W_from_state(state_pm)
               - W_from_state(state_mp) + W_from_state(state_mm)) / (4 * eps**2)

        fermion_mass[I, J] = d2W
        fermion_mass[J, I] = d2W

print("Fermion mass matrix (10x10):")
print("Field labels: M^1_1, M^1_2, M^1_3, M^2_1, M^2_2, M^2_3, M^3_1, M^3_2, M^3_3, X")
print()

# Print in a readable format
labels = ['M11', 'M12', 'M13', 'M21', 'M22', 'M23', 'M31', 'M32', 'M33', 'X']
print(f"{'':>6s}", end='')
for l in labels:
    print(f"{l:>12s}", end='')
print()

for i in range(Nfields):
    print(f"{labels[i]:>6s}", end='')
    for j in range(Nfields):
        val = fermion_mass[i, j]
        if abs(val) < 1e-6:
            print(f"{'0':>12s}", end='')
        else:
            print(f"{val:>12.4f}", end='')
    print()

# Eigenvalues
fermion_evals = np.linalg.eigvalsh(fermion_mass)
print(f"\nFermion mass matrix eigenvalues (MeV):")
for i, ev in enumerate(sorted(fermion_evals)):
    print(f"  lambda_{i+1} = {ev:>16.6f} MeV")

# Physical fermion masses = |eigenvalues|
print(f"\nPhysical fermion masses:")
for i, ev in enumerate(sorted(np.abs(fermion_evals))):
    print(f"  |m_f|_{i+1} = {ev:>16.6f} MeV")

# ============================================================
# PART 3: SCALAR MASS-SQUARED MATRIX
# ============================================================
print("\n" + "=" * 80)
print("PART 3: SCALAR MASS-SQUARED MATRIX")
print("=" * 80)

# In N=1 SUSY with canonical Kahler potential K = sum |Phi_I|^2,
# the scalar potential is V = sum_I |F_I|^2
# where F_I = dW/dPhi_I (complex conjugate of dW*/dPhi*_I)
#
# The scalar mass-squared matrix for the REAL scalar fields
# (phi_I^R, phi_I^I for each complex Phi_I) is:
#
# M^2_scalar = ( F*_{IK} F_{KJ} + F*_I F_{IJ}    , ...  )
#              ( ...                                , ...  )
#
# In the SUSY vacuum where F_I = 0 for all I:
# M^2_scalar_IJ = sum_K F*_{IK} F_{KJ} = (F^dag F)_{IJ}
#
# and the mass-squared eigenvalues of the complex scalars are
# eigenvalues of F^dag F, which are |lambda_fermion|^2.
#
# So in an EXACT SUSY vacuum: m^2_scalar_n = m^2_fermion_n
# and STr[M^2] = 0 automatically.

print("\n--- Scalar masses in SUSY vacuum ---\n")
print("Since F_I = 0 at the SUSY vacuum (all F-terms vanish by construction),")
print("the scalar mass-squared matrix is M^2 = F^dag F.")
print()

# F here is the fermion mass matrix (which is real symmetric in our case)
scalar_mass_sq = fermion_mass @ fermion_mass  # F^T F = F^2 for real symmetric

scalar_evals_sq = np.linalg.eigvalsh(scalar_mass_sq)
print("Scalar mass-squared eigenvalues (MeV^2):")
for i, ev in enumerate(sorted(scalar_evals_sq)):
    print(f"  m^2_{i+1} = {ev:>16.4f} MeV^2  =>  m = {np.sqrt(abs(ev)):>12.4f} MeV")

# Check: scalar masses should equal |fermion masses|
print("\nVerification: scalar masses vs |fermion masses|:")
scalar_masses = np.sqrt(np.abs(sorted(scalar_evals_sq)))
fermion_masses_sorted = sorted(np.abs(fermion_evals))
for i in range(Nfields):
    print(f"  scalar: {scalar_masses[i]:>12.4f} MeV  vs  fermion: {fermion_masses_sorted[i]:>12.4f} MeV")

# ============================================================
# PART 4: SUPERTRACE
# ============================================================
print("\n" + "=" * 80)
print("PART 4: SUPERTRACE CHECK")
print("=" * 80)

# STr[M^2] = sum_bosons m^2_b - sum_fermions m^2_f
# For complex scalars: each complex scalar has 2 real d.o.f.
# Actually, in N=1 SUSY: each chiral multiplet has 1 complex scalar (2 real) + 1 Weyl fermion (2 real)
# The supertrace is: STr[M^2] = sum_n (2*0+1)*m^2_{scalar,n} - sum_n (2*1/2+1)*m^2_{fermion,n}
# Wait -- for the standard supertrace over helicity states:
# STr[M^2] = sum (-1)^{2J} (2J+1) m^2_J
# For scalars (J=0): +1 * m^2  (but we have complex scalars = 2 real scalars)
# For fermions (J=1/2): -2 * m^2  (Weyl fermion = 2 helicity states)
#
# Each chiral multiplet contributes: 2*m^2_scalar - 2*m^2_fermion
# In SUSY vacuum: m^2_scalar = m^2_fermion, so STr = 0.

# But let's compute it explicitly using the mass matrices

# For complex scalars, the mass-squared matrix eigenvalues are already computed
# Each eigenvalue of the complex scalar mass-squared matrix corresponds to one complex scalar
# = 2 real scalar d.o.f.
#
# For Weyl fermions, each eigenvalue of the fermion mass matrix gives mass for one Weyl fermion
# = 2 real fermionic d.o.f.

# STr[M^2] = sum_n 2 * m^2_scalar_n - sum_n 2 * m^2_fermion_n
# = 2 * [Tr(F^dag F) - Tr(F^dag F)] = 0

# But let me compute explicitly:
str_m2_scalars = np.sum(scalar_evals_sq)
str_m2_fermions = np.sum(fermion_evals**2)

print(f"\nSum of scalar m^2: {str_m2_scalars:.6f} MeV^2")
print(f"Sum of fermion m^2: {str_m2_fermions:.6f} MeV^2")
print(f"STr[M^2] = 2 * (scalar - fermion) = {2*(str_m2_scalars - str_m2_fermions):.6e} MeV^2")

# Also compute: do it per d.o.f.
print(f"\nAlternatively:")
print(f"  Tr(M^2_scalar) = Tr(F^dag F) = {np.trace(scalar_mass_sq):.6f}")
print(f"  Tr(M^2_fermion) = Tr(F^2) = {np.trace(fermion_mass @ fermion_mass):.6f}")
print(f"  These are identical because F is real symmetric: F^dag F = F^T F = F^2")
print(f"  Therefore STr[M^2] = 0 identically.")

# ============================================================
# PART 5: ANALYTICAL STRUCTURE OF THE FERMION MASS MATRIX
# ============================================================
print("\n" + "=" * 80)
print("PART 5: ANALYTICAL STRUCTURE")
print("=" * 80)

# Let's understand the fermion mass matrix analytically.
#
# W = m_1 M^1_1 + m_2 M^2_2 + m_3 M^3_3 + X(det M - Lambda^6)
#
# d^2W / dM^a_b dM^c_d = X * d^2(det M) / dM^a_b dM^c_d
# d^2W / dM^a_b dX = d(det M) / dM^a_b = cofactor(M^a_b) evaluated at vacuum
# d^2W / dX dX = 0
#
# The only linear terms in W (giving first derivatives) are m_j M^j_j.
# These don't contribute to second derivatives.
#
# So the fermion mass matrix has two types of entries:
# (1) M-M block: X * (second derivative of det M)
# (2) M-X entries: cofactor of M at vacuum
# (3) X-X: 0
#
# For det M of a 3x3 matrix, the second derivative is:
# d^2(det M)/dM^a_b dM^c_d = epsilon^{acp} epsilon_{bdq} M^q_p
# (Levi-Civita contraction leaving one M factor)
#
# At the diagonal vacuum M = diag(M_1, M_2, M_3):
# This is nonzero only when (a,b) and (c,d) pick different rows/cols and
# the remaining index p,q picks a diagonal element.

print("\n--- Analytical fermion mass matrix entries ---\n")

# Let me compute this analytically
# d^2(det M)/dM^a_b dM^c_d at diagonal M:
# = epsilon_{ace} epsilon_{bdf} M^f_e
# Only nonzero when {a,c,e} is a permutation of {1,2,3} and {b,d,f} is too.
# At diagonal vacuum, M^f_e is nonzero only for f=e.
# So we need e=f, and {a,c,e} permutation, {b,d,e} permutation.
# This means a!=c!=e, b!=d!=e, and e is the remaining index.

# For diagonal entries (a=b, c=d):
# d^2(det M)/dM^a_a dM^c_c = epsilon_{ace} epsilon_{ace} M_e = M_e
# where e is the remaining index (if a!=c).
# If a=c: epsilon_{aae}=0, so zero.

# For off-diagonal entries (a!=b):
# d^2(det M)/dM^a_b dM^c_d = epsilon_{ace} epsilon_{bde} M_e
# Need a,c,e all different and b,d,e all different.
# e is fixed by a,c: e = 6-a-c (using 1-indexing where indices are 1,2,3)
# Similarly need e = 6-b-d.
# So 6-a-c = 6-b-d => a+c = b+d.
# And all distinct.

# Let me just build it explicitly.
print("Building analytical fermion mass matrix...\n")

# Using 0-indexed (i,j) for M^{i+1}_{j+1}
# Field index: I = 3*i + j for M^{i+1}_{j+1}, I=9 for X

# Levi-Civita
def levi_civita(i, j, k):
    """3D Levi-Civita symbol (0-indexed)."""
    if (i, j, k) in [(0,1,2), (1,2,0), (2,0,1)]:
        return 1
    elif (i, j, k) in [(0,2,1), (2,1,0), (1,0,2)]:
        return -1
    return 0

# Cofactor at diagonal vacuum
# cof(M)^j_i = (1/2) epsilon^{jac} epsilon_{ibd} M^b_a M^d_c
# At diagonal vacuum: M^b_a = M_a * delta^b_a
# cof(M)^j_i = (1/2) epsilon^{jac} epsilon_{iac} M_a M_c  for diagonal
# = delta^j_i * (product of other two M's)

# Actually let me just compute cofactor numerically
def cofactor_matrix(M):
    """Cofactor matrix of 3x3."""
    cof = np.zeros((3,3))
    for i in range(3):
        for j in range(3):
            minor = np.delete(np.delete(M, i, axis=0), j, axis=1)
            cof[i,j] = (-1)**(i+j) * np.linalg.det(minor)
    return cof

cof_vac = cofactor_matrix(M_vac)
print("Cofactor matrix at vacuum:")
for i in range(3):
    for j in range(3):
        if abs(cof_vac[i,j]) > 1e-6:
            print(f"  cof(M)^{i+1}_{j+1} = {cof_vac[i,j]:.4f}")

# Build analytical fermion mass matrix
F_analytic = np.zeros((10, 10))

# M-M block: F_{(a,b),(c,d)} = X * d^2(det M)/dM^a_b dM^c_d
# d^2(det M)/dM^a_b dM^c_d = sum_e,f epsilon_{ace} epsilon_{bdf} M^f_e
# At diagonal vacuum M^f_e = M_e * delta_{fe}
# = sum_e epsilon_{ace} epsilon_{bde} M_e

for a in range(3):
    for b in range(3):
        I = 3*a + b
        for c in range(3):
            for d in range(3):
                J = 3*c + d
                val = 0.0
                for e in range(3):
                    val += levi_civita(a, c, e) * levi_civita(b, d, e) * M_uds[e]
                F_analytic[I, J] = X_vac * val

# M-X block: F_{(a,b), X} = d(det M)/dM^a_b at vacuum = cof(M)^b_a (transpose of cofactor)
# Wait: d(det M)/dM^i_j = cof(M)^j_i
# Actually for a matrix with upper-lower indices M^i_j:
# det M = (1/6) eps^{ijk} eps_{lmn} M^l_i M^m_j M^n_k
# d(det M)/dM^a_b = (1/2) eps^{abc} eps_{bef} M^e_? ...
# Let me just use the numerical result.

# d(det M)/dM^a_b at diagonal vacuum:
# For b=a: d(det M)/dM^a_a = product of M_k for k != a
# For b!=a: zero at diagonal vacuum? Let's check.

for a in range(3):
    for b in range(3):
        I = 3*a + b
        # Compute d(det M)/dM^a_b numerically
        M_p = M_vac.copy()
        M_m = M_vac.copy()
        M_p[a, b] += eps
        M_m[a, b] -= eps
        ddetdM = (det3(M_p) - det3(M_m)) / (2*eps)
        F_analytic[I, 9] = ddetdM
        F_analytic[9, I] = ddetdM

# X-X: 0 (already zero)

print("\nAnalytical fermion mass matrix (nonzero entries):")
for i in range(10):
    for j in range(i, 10):
        if abs(F_analytic[i, j]) > 1e-6:
            print(f"  F[{labels[i]}, {labels[j]}] = {F_analytic[i,j]:.4f}")

# Compare with numerical
print(f"\nMax difference analytical vs numerical: {np.max(np.abs(F_analytic - fermion_mass)):.6e}")

# Eigenvalues of analytical matrix
evals_analytic = np.linalg.eigvalsh(F_analytic)
print(f"\nFermion mass eigenvalues (analytical):")
for i, ev in enumerate(sorted(evals_analytic)):
    print(f"  lambda_{i+1} = {ev:>16.6f} MeV")

# ============================================================
# PART 6: ISS REGIME (N_f = 4, N_c = 3)
# ============================================================
print("\n" + "=" * 80)
print("PART 6: ISS REGIME (N_f = 4, N_c = 3) -- (u, d, s, c)")
print("=" * 80)

# ISS setup: N_f = 4 > N_c = 3, magnetic dual has N_c^mag = N_f - N_c = 1
#
# Magnetic superpotential:
# W = h q^i Phi_i^j qtilde_j - h mu^2 Tr(Phi) + Tr(m Phi)
#
# where Phi is the 4x4 meson matrix, q^i and qtilde_j are magnetic quarks
# with a single magnetic color.
#
# mu^2 = Lambda^2 (matching to electric theory)
# h is the magnetic Yukawa (we'll take h = sqrt(3) as in iss_cw_koide.md)

print("\n--- ISS Tree-Level Spectrum ---\n")

h_mag = np.sqrt(3.0)
mu_ISS = 92.0  # MeV (from matching, same as iss_cw_koide.md)
# Actually, let's use mu^2 = Lambda_mag^2 matching. Different conventions exist.
# Following the previous computation in iss_cw_koide.md:
mu_sq = mu_ISS**2  # 8464 MeV^2

m_ISS = np.array([m_u, m_d, m_s, m_c])  # 4 flavors
flav_ISS = ['u', 'd', 's', 'c']

# ISS vacuum: q^1 = qtilde_1 = phi_0
phi_0 = np.sqrt(mu_sq - m_ISS[0]/h_mag)
print(f"phi_0 = sqrt(mu^2 - m_u/h) = {phi_0:.4f} MeV")

# F-terms for broken directions (a = 2,3,4 = d,s,c)
F_ISS = np.array([h_mag * mu_sq - m_ISS[a] for a in range(1, 4)])
print(f"\nF-terms (SUSY breaking):")
for i, f in enumerate(flav_ISS[1:]):
    print(f"  F_{f} = h*mu^2 - m_{f} = {F_ISS[i]:.4f} MeV^2")

# Tree-level masses:
print(f"\n--- Link field sector (Sector B) ---")
print(f"Fermion mass (universal): m_F = h * phi_0 = {h_mag * phi_0:.4f} MeV")
print()

print(f"{'Flavor':>8s} {'m_q (MeV)':>12s} {'m^2_B+ (MeV^2)':>16s} {'m^2_B- (MeV^2)':>16s} {'m_F (MeV)':>12s}")
for i, f in enumerate(flav_ISS[1:]):
    m2_plus = h_mag**2 * phi_0**2 + h_mag * F_ISS[i]
    m2_minus = h_mag**2 * phi_0**2 - h_mag * F_ISS[i]
    m_F = h_mag * phi_0
    print(f"{f:>8s} {m_ISS[i+1]:>12.2f} {m2_plus:>16.2f} {m2_minus:>16.2f} {m_F:>12.2f}")

# Goldstino
print(f"\nGoldstino: massless (1 state)")
print(f"Heavy fermion pair: m = sqrt(2) * h * phi_0 = {np.sqrt(2) * h_mag * phi_0:.4f} MeV")

# Pseudo-modulus
print(f"\nPseudo-modulus X: tree-level FLAT, stabilized by CW at origin")

# CW masses for broken mesons
print(f"\n--- Coleman-Weinberg masses for broken mesons chi_a ---")
print(f"(From iss_cw_koide.md, verified previously)")
print(f"m_CW ~ 24-27 MeV for all flavors (nearly degenerate)")

# Supertrace in ISS
print(f"\n--- Supertrace in ISS ---")
print(f"\nFor each link-field sector:")
for i, f in enumerate(flav_ISS[1:]):
    m2_plus = h_mag**2 * phi_0**2 + h_mag * F_ISS[i]
    m2_minus = h_mag**2 * phi_0**2 - h_mag * F_ISS[i]
    m2_Y = h_mag**2 * phi_0**2  # Y field (no splitting)
    m2_F = (h_mag * phi_0)**2

    # 2 complex scalars with m^2_plus, 2 with m^2_minus, 2 Y-scalars with m^2_Y
    # Wait, let me be more careful about counting.
    # Each broken flavor a has:
    # - q^a, qtilde_a: 2 complex scalars + 2 Weyl fermions
    # - Phi^1_a, Phi^a_1: 2 complex scalars + 2 Weyl fermions
    # Total: 4 complex scalars, 4 Weyl fermions per flavor
    # But mixing occurs...

    # From the iss_cw_koide.md computation:
    # Boson masses: m^2_{B+}, m^2_{B-} (from q-qtilde mixing), m^2_Y (from Phi off-diag)
    # Fermion mass: m_F (universal)

    # STr for this sector:
    # Scalars: 1*(m^2_plus) + 1*(m^2_minus) + 2*(m^2_Y) = m^2_plus + m^2_minus + 2*m^2_Y
    # = 2*h^2*phi_0^2 + 2*h^2*phi_0^2 = 4*h^2*phi_0^2
    # Fermions: 4 * m^2_F = 4*h^2*phi_0^2
    # STr = 4*h^2*phi_0^2 - 4*h^2*phi_0^2 = 0 ... but that's tree-level in the unbroken SUSY part.

    # Actually in the ISS with SUSY breaking, STr[M^2] = 0 still holds at tree level
    # because the SUSY breaking is F-term type and we have canonical Kahler.

    str_boson = m2_plus + m2_minus + 2*m2_Y
    str_fermion = 4 * m2_F

    # Hmm, need to be more careful with the counting. Let me just note:
    str_sector = (m2_plus + m2_minus) - 2 * m2_F
    # = (h^2 phi_0^2 + h F) + (h^2 phi_0^2 - h F) - 2 h^2 phi_0^2 = 0
    print(f"  STr[M^2] (link sector {f}) = {str_sector:.6e} MeV^2  (= 0 exactly)")

print(f"\nSTr[M^2] = 0 at tree level in ISS (as expected: F-term breaking with canonical Kahler)")

# ============================================================
# COMPREHENSIVE KOIDE TABLE
# ============================================================
print("\n" + "=" * 80)
print("COMPREHENSIVE KOIDE TABLE")
print("=" * 80)

# All triples from {u,d,s,c,b}
all_flavors = ['u', 'd', 's', 'c', 'b']
all_masses = [m_u, m_d, m_s, m_c, m_b]

print(f"\n{'Triple':>12s} {'Q(m)':>10s} {'Q(1/m)':>10s} {'Q(1/m)-2/3':>12s} {'%dev':>8s}")
print("-" * 60)

for combo in combinations(range(5), 3):
    names = ','.join([all_flavors[i] for i in combo])
    ms = [all_masses[i] for i in combo]
    inv_ms = [1.0/all_masses[i] for i in combo]

    q_m = koide_Q(*ms)
    q_inv = koide_Q(*inv_ms)
    dev = q_inv - 2./3.
    pct = abs(dev) / (2./3.) * 100

    print(f"({names:>8s})  {q_m:>10.6f} {q_inv:>10.6f} {dev:>+12.6f} {pct:>8.4f}")

# Special: leptons
m_e = 0.51099895  # MeV
m_mu = 105.6583755
m_tau = 1776.86
print(f"\n{'(e,mu,tau)':>12s}  {koide_Q(m_e, m_mu, m_tau):>10.6f} {koide_Q(1/m_e, 1/m_mu, 1/m_tau):>10.6f}")

# Seesaw meson Koide for various triples
print(f"\n\n{'Triple':>12s} {'C (MeV^2)':>14s} {'M_1':>12s} {'M_2':>12s} {'M_3':>12s} {'Q(M)':>10s}")
print("-" * 80)

for combo in combinations(range(5), 3):
    names = ','.join([all_flavors[i] for i in combo])
    ms = np.array([all_masses[i] for i in combo])
    C = Lambda**2 * np.prod(ms)**(1./3.)
    Ms = C / ms
    q = koide_Q(*Ms)
    print(f"({names:>8s})  {C:>14.4f} {Ms[0]:>12.4f} {Ms[1]:>12.4f} {Ms[2]:>12.4f} {q:>10.6f}")

print("\n\nDone.")
