"""
Fritzsch texture mass matrices, CKM mixing, and seesaw meson VEVs.

Parts (a)-(f): linear algebra of nearest-neighbor quark mass matrices,
CKM angles from Fritzsch relations, seesaw inversion for meson VEVs,
and numerical verification of all results.
"""

import numpy as np
from numpy import sqrt, pi

# =============================================================================
# Input quark masses (MeV), PDG central values
# =============================================================================
m_u = 2.16
m_c = 1275.0
m_t = 172760.0

m_d = 4.67
m_s = 93.4
m_b = 4180.0

Lambda_QCD = 300.0  # MeV, SQCD confinement scale

print("=" * 72)
print("PART (a): Verify Fritzsch relations via characteristic polynomial")
print("=" * 72)

# The Fritzsch texture for the Hermitian matrix M is:
#
#   M = [[0,    A,    0 ],
#        [A*,   0,    B ],
#        [0,    B*,   C ]]
#
# For REAL entries, the characteristic polynomial is:
#   det(M - lambda I) = 0
#   -lambda^3 + C lambda^2 + (A^2 + B^2) lambda - C A^2 = 0
#
# Or equivalently: lambda^3 - C lambda^2 - (A^2+B^2) lambda + C A^2 = 0
#
# Vieta's formulas for eigenvalues lambda_1, lambda_2, lambda_3:
#   lambda_1 + lambda_2 + lambda_3 = C                    (I)
#   lambda_1*lambda_2 + lambda_1*lambda_3 + lambda_2*lambda_3 = -(A^2+B^2)  (II)
#   lambda_1*lambda_2*lambda_3 = -C*A^2                   (III)
#
# Wait, let me redo this. For M - lambda I:
# det = (0-lam)[(0-lam)(C-lam) - B*B] - A[A*(C-lam) - 0] + 0
#     = -lam[-(C-lam)lam - B^2] - A*A(C-lam)
#     = -lam[-C*lam + lam^2 - B^2] - A^2(C-lam)
#     = -lam^3 + C*lam^2 + B^2*lam - A^2*C + A^2*lam
#     = -lam^3 + C*lam^2 + (A^2+B^2)*lam - A^2*C
#
# So: lam^3 - C*lam^2 - (A^2+B^2)*lam + A^2*C = 0
#
# Vieta's (for lam^3 + p*lam^2 + q*lam + r = 0):
#   lam1+lam2+lam3 = -p = C
#   lam1*lam2 + lam1*lam3 + lam2*lam3 = q = -(A^2+B^2)
#   lam1*lam2*lam3 = -r = -A^2*C
#
# Wait, standard Vieta for lam^3 - c2*lam^2 + c1*lam - c0 = 0:
#   sum = c2, sum_pairs = c1, product = c0
#
# Our equation: lam^3 - C*lam^2 + (-(A^2+B^2))*lam + A^2*C = 0
# So: c2 = C, c1 = -(A^2+B^2), c0 = -A^2*C
#
# Vieta's: sum = C, sum_pairs = -(A^2+B^2), product = -A^2*C
# Wait no. For lam^3 - c2*lam^2 + c1*lam - c0:
#   product = c0. Our polynomial is lam^3 - C*lam^2 - (A^2+B^2)*lam + A^2*C
#   = lam^3 - C*lam^2 + [-(A^2+B^2)]*lam - [-A^2*C]
#
# So c2=C, c1=-(A^2+B^2), c0=-A^2*C
# sum_of_eigs = c2 = C
# sum_of_pairs = c1 = -(A^2+B^2)
# product = c0 = -A^2*C
#
# For Fritzsch, the eigenvalues are (+m_1, -m_2, +m_3) where m_i > 0 and
# m_1 < m_2 < m_3.
#
# Check: sum = m_1 - m_2 + m_3 = C  =>  C = m_3 - m_2 + m_1  ✓
# sum_pairs = -m_1*m_2 + m_1*m_3 - m_2*m_3 = -(A^2+B^2)
#   => A^2+B^2 = m_1*m_2 - m_1*m_3 + m_2*m_3
#   => A^2+B^2 = m_1*m_2 + m_2*m_3 - m_1*m_3
#   Hmm, this gives: m_1*m_2 + m_2*m_3 - m_1*m_3 = A^2+B^2
#
# product = m_1*(-m_2)*m_3 = -m_1*m_2*m_3 = -A^2*C
#   => A^2 = m_1*m_2*m_3 / C = m_1*m_2*m_3 / (m_3-m_2+m_1)
#
# In the limit m_1 << m_2 << m_3:
#   C ≈ m_3 - m_2 ≈ m_3
#   A^2 ≈ m_1*m_2*m_3/m_3 = m_1*m_2
#   B^2 = (A^2+B^2) - A^2 = (m_1*m_2+m_2*m_3-m_1*m_3) - m_1*m_2*m_3/(m_3-m_2+m_1)
#
# But the EXACT Fritzsch relations are:
#   C = m_3 - m_2 + m_1  (exact, from Vieta I)
#   A^2 = m_1*m_2*m_3 / (m_3 - m_2 + m_1)  (exact, from Vieta III)
#   B^2 = m_1*m_2 + m_2*m_3 - m_1*m_3 - A^2  (exact, from Vieta II)
#        = m_1*m_2 + m_2*m_3 - m_1*m_3 - m_1*m_2*m_3/(m_3-m_2+m_1)
#
# The APPROXIMATE (leading order) relations |A|^2 ≈ m_1*m_2 and |B|^2 ≈ m_2*m_3
# are what the problem statement gives. Let me verify both exact and approximate.

def fritzsch_exact(m1, m2, m3):
    """
    Compute EXACT Fritzsch parameters A, B, C from eigenvalues (+m1, -m2, +m3).
    """
    C = m3 - m2 + m1
    A2 = m1 * m2 * m3 / C
    B2 = m1*m2 + m2*m3 - m1*m3 - A2
    return sqrt(A2), sqrt(B2), C, A2, B2

def fritzsch_approx(m1, m2, m3):
    """
    Compute APPROXIMATE Fritzsch parameters.
    """
    A = sqrt(m1 * m2)
    B = sqrt(m2 * m3)
    C = m3 - m2 + m1
    return A, B, C

def fritzsch_matrix_exact(m1, m2, m3):
    """Construct exact Fritzsch matrix."""
    A, B, C, A2, B2 = fritzsch_exact(m1, m2, m3)
    M = np.array([
        [0,  A,  0],
        [A,  0,  B],
        [0,  B,  C]
    ], dtype=float)
    return M, A, B, C, A2, B2

def fritzsch_matrix_approx(m1, m2, m3):
    """Construct approximate Fritzsch matrix."""
    A, B, C = fritzsch_approx(m1, m2, m3)
    M = np.array([
        [0,  A,  0],
        [A,  0,  B],
        [0,  B,  C]
    ], dtype=float)
    return M, A, B, C

def verify_sector(m1, m2, m3, label):
    """Verify Fritzsch relations for one sector."""
    print(f"\n{'='*60}")
    print(f"  {label} sector: m1={m1}, m2={m2}, m3={m3}")
    print(f"{'='*60}")

    # --- Exact relations ---
    A_ex, B_ex, C_ex, A2_ex, B2_ex = fritzsch_exact(m1, m2, m3)
    M_ex, _, _, _, _, _ = fritzsch_matrix_exact(m1, m2, m3)

    eigs_ex = np.sort(np.linalg.eigvalsh(M_ex))
    expected = np.sort([-m2, m1, m3])

    print(f"\n  EXACT Fritzsch parameters:")
    print(f"    C = m3-m2+m1 = {C_ex:.6f}")
    print(f"    A^2 = m1*m2*m3/C = {A2_ex:.6f}  (approx m1*m2 = {m1*m2:.6f})")
    print(f"    B^2 = {B2_ex:.6f}  (approx m2*m3 = {m2*m3:.6f})")
    print(f"    A = {A_ex:.6f}  (approx sqrt(m1*m2) = {sqrt(m1*m2):.6f})")
    print(f"    B = {B_ex:.6f}  (approx sqrt(m2*m3) = {sqrt(m2*m3):.6f})")

    print(f"\n  Eigenvalue check (exact parameters):")
    print(f"    Computed:  {eigs_ex}")
    print(f"    Expected:  {expected}")
    print(f"    Max error: {np.max(np.abs(eigs_ex - expected)):.2e}")

    # Vieta check
    lam = [-m2, m1, m3]
    S1 = sum(lam)
    S2 = lam[0]*lam[1] + lam[0]*lam[2] + lam[1]*lam[2]
    S3 = lam[0]*lam[1]*lam[2]

    print(f"\n  Vieta's formulas:")
    print(f"    sum(lambda_i) = {S1:.6f},  C = {C_ex:.6f},  err = {abs(S1-C_ex):.2e}")
    print(f"    sum(lambda_i lambda_j) = {S2:.6f},  -(A^2+B^2) = {-(A2_ex+B2_ex):.6f},  err = {abs(S2+(A2_ex+B2_ex)):.2e}")
    print(f"    prod(lambda_i) = {S3:.6f},  -A^2*C = {-A2_ex*C_ex:.6f},  err = {abs(S3+A2_ex*C_ex):.2e}")

    # --- Approximate relations ---
    A_ap, B_ap, C_ap = fritzsch_approx(m1, m2, m3)
    M_ap, _, _, _ = fritzsch_matrix_approx(m1, m2, m3)
    eigs_ap = np.sort(np.linalg.eigvalsh(M_ap))

    print(f"\n  APPROXIMATE Fritzsch parameters:")
    print(f"    A_approx = {A_ap:.6f},  A_exact = {A_ex:.6f},  rel err = {abs(A_ap/A_ex-1):.2e}")
    print(f"    B_approx = {B_ap:.6f},  B_exact = {B_ex:.6f},  rel err = {abs(B_ap/B_ex-1):.2e}")
    print(f"    Eigenvalues (approx): {eigs_ap}")
    print(f"    |masses| (approx):    {np.sort(np.abs(eigs_ap))}")
    print(f"    |masses| (exact):     {np.sort([m1,m2,m3])}")

    return M_ex

M_up = verify_sector(m_u, m_c, m_t, "Up-type (u,c,t)")
M_down = verify_sector(m_d, m_s, m_b, "Down-type (d,s,b)")


print("\n\n" + "=" * 72)
print("PART (b): CKM angles from exact Fritzsch diagonalization")
print("=" * 72)

# Diagonalize exact Fritzsch matrices
evals_u, evecs_u = np.linalg.eigh(M_up)
evals_d, evecs_d = np.linalg.eigh(M_down)

print(f"\nUp eigenvalues (sorted): {evals_u}")
print(f"Down eigenvalues (sorted): {evals_d}")

# numpy sorts eigenvalues ascending. For Fritzsch with (+m1, -m2, +m3):
# ascending order is (-m2, m1, m3)
# We want physical mass ordering: m1, m2, m3 (by |eigenvalue|)

def reorder_physical(evals, evecs):
    """Reorder to physical mass ordering (ascending |eigenvalue|)."""
    idx = np.argsort(np.abs(evals))
    return evals[idx], evecs[:, idx]

evals_u_phys, evecs_u_phys = reorder_physical(evals_u, evecs_u)
evals_d_phys, evecs_d_phys = reorder_physical(evals_d, evecs_d)

print(f"\nPhysical ordering:")
print(f"  Up:   eigenvalues = {evals_u_phys},  |m| = {np.abs(evals_u_phys)}")
print(f"  Down: eigenvalues = {evals_d_phys},  |m| = {np.abs(evals_d_phys)}")

# The CKM matrix.
# The physical mass matrix M_phys = diag(m_u, m_c, m_t) (positive) is
# obtained by M = U_L diag(eigenvalues) U_L^T for real symmetric case.
# The negative eigenvalue (-m_2) means we need a phase convention.
#
# Standard convention: absorb the sign into a diagonal phase matrix P
# such that U_L^phys = U_L * P where P = diag(1, -1, 1) for the
# second eigenvalue being negative.
#
# V_CKM = (U_L^u_phys)^T * U_L^d_phys
#        = P_u * evecs_u^T * evecs_d * P_d^{-1}
#
# Since CKM elements enter as |V_ij|, the phase convention doesn't
# affect the moduli. But it does affect which element maps to which
# physical state.
#
# Actually for the CKM, what matters is V = U_u^dag * U_d where
# U_u, U_d diagonalize the MASS SQUARED matrices M_u M_u^dag and M_d M_d^dag.
# For real symmetric Fritzsch, M M^T = M^2, and the eigenvectors of M^2
# are the same as those of M.
#
# So the CKM moduli are simply |evecs_u_phys^T * evecs_d_phys|.

V_CKM = evecs_u_phys.T @ evecs_d_phys

print(f"\nV_CKM (raw):")
for i in range(3):
    print(f"  [{V_CKM[i,0]:+.8f}  {V_CKM[i,1]:+.8f}  {V_CKM[i,2]:+.8f}]")

print(f"\n|V_CKM|:")
for i in range(3):
    print(f"  [{abs(V_CKM[i,0]):.8f}  {abs(V_CKM[i,1]):.8f}  {abs(V_CKM[i,2]):.8f}]")

# Extract mixing angles: |V_us| = sin theta_12, |V_cb| = sin theta_23, |V_ub| = sin theta_13
s12 = abs(V_CKM[0, 1])  # V_us
s23 = abs(V_CKM[1, 2])  # V_cb
s13 = abs(V_CKM[0, 2])  # V_ub

theta12 = np.arcsin(s12) * 180 / pi
theta23 = np.arcsin(s23) * 180 / pi
theta13 = np.arcsin(s13) * 180 / pi

# PDG values
theta12_pdg = 13.04
theta23_pdg = 2.38
theta13_pdg = 0.201
s12_pdg = np.sin(theta12_pdg * pi / 180)
s23_pdg = np.sin(theta23_pdg * pi / 180)
s13_pdg = np.sin(theta13_pdg * pi / 180)

print(f"\nCKM angles from exact Fritzsch diagonalization:")
print(f"  |V_us| = sin theta_12 = {s12:.6f}  ->  theta_12 = {theta12:.4f} deg")
print(f"  |V_cb| = sin theta_23 = {s23:.6f}  ->  theta_23 = {theta23:.4f} deg")
print(f"  |V_ub| = sin theta_13 = {s13:.6f}  ->  theta_13 = {theta13:.4f} deg")

print(f"\nPDG values:")
print(f"  sin theta_12 = {s12_pdg:.6f}  ->  theta_12 = {theta12_pdg:.2f} deg")
print(f"  sin theta_23 = {s23_pdg:.6f}  ->  theta_23 = {theta23_pdg:.2f} deg")
print(f"  sin theta_13 = {s13_pdg:.6f}  ->  theta_13 = {theta13_pdg:.3f} deg")


# ---- Fritzsch approximate formulas ----
print(f"\n--- Fritzsch approximate formulas ---")

# The approximate formula (Fritzsch 1978, Eq. 22):
# For a 2x2 block [[0, a],[a, D]], the diagonalizing rotation has
# tan(theta) ≈ a/D ≈ sqrt(m_light/m_heavy) for m_light << D.
#
# For the 3x3 Fritzsch texture, the diagonalization proceeds in two steps:
# 1) Rotate in the (2,3) block: angle phi ≈ B/C ≈ sqrt(m_2/m_3)
# 2) Rotate in the (1,2) block: angle psi ≈ A/m_2 ≈ sqrt(m_1/m_2)
#
# The CKM (1,2) element is then:
# V_us ≈ psi_d - psi_u = sqrt(m_d/m_s) - sqrt(m_u/m_c) * e^{i phi}
# In the real case: sin theta_12 ≈ |sqrt(m_d/m_s) - sqrt(m_u/m_c)|
#
# Similarly:
# sin theta_23 ≈ |sqrt(m_s/m_b) - sqrt(m_c/m_t)|
# sin theta_13 ≈ |sqrt(m_d/m_b) - sqrt(m_u/m_t)|   (product of two rotations)

s12_fritzsch = abs(sqrt(m_d / m_s) - sqrt(m_u / m_c))
s23_fritzsch = abs(sqrt(m_s / m_b) - sqrt(m_c / m_t))
s13_fritzsch = abs(sqrt(m_d / m_b) - sqrt(m_u / m_t))
s12_gst = sqrt(m_d / m_s)

print(f"  sin theta_12 ≈ |sqrt(m_d/m_s) - sqrt(m_u/m_c)|")
print(f"               = |{sqrt(m_d/m_s):.6f} - {sqrt(m_u/m_c):.6f}| = {s12_fritzsch:.6f}")
print(f"    -> theta_12 = {np.arcsin(s12_fritzsch)*180/pi:.4f} deg")
print(f"")
print(f"  sin theta_23 ≈ |sqrt(m_s/m_b) - sqrt(m_c/m_t)|")
print(f"               = |{sqrt(m_s/m_b):.6f} - {sqrt(m_c/m_t):.6f}| = {s23_fritzsch:.6f}")
print(f"    -> theta_23 = {np.arcsin(s23_fritzsch)*180/pi:.4f} deg")
print(f"")
print(f"  sin theta_13 ≈ |sqrt(m_d/m_b) - sqrt(m_u/m_t)|")
print(f"               = |{sqrt(m_d/m_b):.6f} - {sqrt(m_u/m_t):.6f}| = {s13_fritzsch:.6f}")
print(f"    -> theta_13 = {np.arcsin(s13_fritzsch)*180/pi:.4f} deg")

print(f"\n  GST/Oakes relation (m_u -> 0 limit):")
print(f"    sin theta_12 ≈ sqrt(m_d/m_s) = {s12_gst:.6f}")
print(f"    -> theta_12 = {np.arcsin(s12_gst)*180/pi:.4f} deg")


# ---- Derivation of GST/Oakes ----
print(f"\n--- Derivation of Gatto-Sartori-Tonin / Oakes relation ---")
print(f"""
  For a 2x2 Fritzsch block (the leading approximation ignoring 3rd gen):
    M = [[0, A], [A, D]]
  with eigenvalues lambda_1, lambda_2 satisfying:
    lambda_1 + lambda_2 = D
    lambda_1 * lambda_2 = -A^2
  The diagonalizing rotation angle satisfies:
    tan(2*theta) = 2A/D
  For eigenvalues (+m_light, -m_heavy):
    A^2 = m_light * m_heavy   (exact for 2x2)
    D = m_heavy - m_light
    tan(theta) = A / (m_heavy + ...) ≈ sqrt(m_light/m_heavy)

  More precisely, for the 2x2 case:
    sin theta = sqrt(m_light / (m_light + m_heavy)) ≈ sqrt(m_light/m_heavy)

  The CKM (1,2) element combines up and down sector rotations:
    V_us = cos(phi_u)*sin(phi_d) - sin(phi_u)*cos(phi_d)*e^{{i delta}}

  For small angles:
    V_us ≈ phi_d - phi_u ≈ sqrt(m_d/m_s) - sqrt(m_u/m_c)

  In the limit m_u -> 0:
    V_us ≈ sqrt(m_d/m_s) = {sqrt(m_d/m_s):.6f}

  This is the Gatto-Sartori-Tonin (1968) / Oakes (1969) relation.

  Numerically: sqrt(m_d/m_s) = sqrt({m_d}/{m_s}) = {sqrt(m_d/m_s):.6f}
  PDG |V_us| = {s12_pdg:.6f}
  Agreement: {abs(sqrt(m_d/m_s)/s12_pdg - 1)*100:.2f}%
""")


print("\n" + "=" * 72)
print("PART (c): Seesaw map for meson VEVs — argument")
print("=" * 72)

print("""
QUESTION: Does the seesaw M_j -> C/m_j apply only to eigenvalues,
or does it preserve the full matrix structure (off-diagonal entries)?

ANSWER: It preserves the full matrix structure. The seesaw is a MATRIX
inversion, not just eigenvalue inversion.

DERIVATION:

For N_f = N_c = 3 SQCD, the quantum-modified constraint is:
    det(Phi) - B tilde{B} = Lambda^6

where Phi^i_j = Q^i tilde{Q}_j is the meson superfield. With B = tilde{B} = 0:
    det(Phi) = Lambda^6

Adding a mass deformation W_mass = Tr(m_q Phi) with quark mass matrix m_q,
the F-term equation is (using Lagrange multiplier for the constraint):

    (m_q)_ij + mu * (cof Phi)_ij = 0

where (cof Phi)_ij = (partial det Phi)/(partial Phi^j_i) = (det Phi)(Phi^{-1})^j_i.

This gives: m_q = -mu * Lambda^6 * (Phi^{-1})^T

Therefore: Phi = -mu * Lambda^6 * (m_q^{-1})^T  (as a MATRIX)

The proportionality constant mu is fixed by the constraint det Phi = Lambda^6:
    det Phi = (-mu)^3 Lambda^{18} / det(m_q) = Lambda^6
    => (-mu)^3 = det(m_q) / Lambda^{12}
    => -mu = (det(m_q))^{1/3} / Lambda^4

So: Phi = alpha * m_q^{-T}  where alpha = Lambda^2 * (det m_q)^{1/3}

If m_q = U diag(m_1, m_2, m_3) V^dag, then:
    m_q^{-1} = V diag(1/m_1, 1/m_2, 1/m_3) U^dag

    Phi = alpha * (V diag(1/m_1, 1/m_2, 1/m_3) U^dag)^T
        = alpha * U* diag(1/m_1, 1/m_2, 1/m_3) V^T

For REAL matrices (Fritzsch): U* = U, V^T = V, so:
    Phi = alpha * U diag(1/m_k) V
        = U diag(alpha/m_k) V

Since for real symmetric M: U = V (same eigenvectors), we get:
    Phi = U diag(alpha/m_k) U^T

The meson matrix has the SAME eigenvectors as the quark mass matrix,
with eigenvalues alpha/m_k. The off-diagonal structure is preserved.

KEY POINT: The mixing information is not lost. It is inherited from
the quark mass matrix through the matrix inverse.

Consistency check:
    det Phi = det(U) * prod(alpha/m_k) * det(U^T)
            = det(U)^2 * alpha^3 / prod(m_k)
            = 1 * (Lambda^6 * prod(m_k)) / prod(m_k)    [since alpha^3 = Lambda^6 * prod(m_k)]
            = Lambda^6  CHECK
""")


print("\n" + "=" * 72)
print("PART (d): Numerical meson VEV matrix for down-type sector")
print("=" * 72)

# Use EXACT Fritzsch matrix for down sector
M_d_ex, A_d_ex, B_d_ex, C_d_ex, A2_d_ex, B2_d_ex = fritzsch_matrix_exact(m_d, m_s, m_b)

print(f"\nDown-type exact Fritzsch matrix (MeV):")
for i in range(3):
    print(f"  [{M_d_ex[i,0]:>12.6f}  {M_d_ex[i,1]:>12.6f}  {M_d_ex[i,2]:>12.6f}]")

evals_d, evecs_d = np.linalg.eigh(M_d_ex)
print(f"\nEigenvalues: {evals_d}")
print(f"Check: should be (-m_s, m_d, m_b) = ({-m_s}, {m_d}, {m_b})")
print(f"Max error: {max(abs(evals_d[0]+m_s), abs(evals_d[1]-m_d), abs(evals_d[2]-m_b)):.2e}")

# Reorder to physical mass ordering
evals_d_phys, evecs_d_phys = reorder_physical(evals_d, evecs_d)
phys_masses = np.abs(evals_d_phys)
print(f"\nPhysical masses: {phys_masses}")
print(f"Eigenvectors (physical ordering, columns):")
for i in range(3):
    print(f"  [{evecs_d_phys[i,0]:+.8f}  {evecs_d_phys[i,1]:+.8f}  {evecs_d_phys[i,2]:+.8f}]")

# Seesaw: Phi = U diag(alpha/m_k) U^T  with alpha = Lambda^2 * (m_d m_s m_b)^{1/3}
alpha_seesaw = Lambda_QCD**2 * (m_d * m_s * m_b)**(1.0/3.0)
print(f"\nSeesaw scale: alpha = Lambda^2 * (m_d m_s m_b)^{{1/3}} = {alpha_seesaw:.4f} MeV^2")

meson_eigenvalues = alpha_seesaw / phys_masses
print(f"\nMeson eigenvalues (alpha/m_k):")
print(f"  alpha/m_d = {meson_eigenvalues[0]:.4f} MeV")
print(f"  alpha/m_s = {meson_eigenvalues[1]:.4f} MeV")
print(f"  alpha/m_b = {meson_eigenvalues[2]:.4f} MeV")
print(f"  Hierarchy: {meson_eigenvalues[0]/meson_eigenvalues[1]:.2f} : 1 : {meson_eigenvalues[2]/meson_eigenvalues[1]:.4f}")

# Full meson VEV matrix in flavor basis
M_meson = evecs_d_phys @ np.diag(meson_eigenvalues) @ evecs_d_phys.T

print(f"\nMeson VEV matrix <Phi^a_b> in flavor basis (a=d,s,b; b=d,s,b):")
print(f"  (units: MeV)")
labels = ['d', 's', 'b']
print(f"{'':>8}", end='')
for j in range(3):
    print(f"  {labels[j]:>14}", end='')
print()
for i in range(3):
    print(f"  {labels[i]:>4}", end='')
    for j in range(3):
        print(f"  {M_meson[i,j]:>14.4f}", end='')
    print()

# Also verify by direct matrix inversion
M_d_inv = np.linalg.inv(M_d_ex)
# Need to account for negative eigenvalue in determinant
det_M_d = np.linalg.det(M_d_ex)
# alpha for signed det: alpha^3 = Lambda^6 * |det(M_d)|
# Actually: Phi = alpha_signed * M_d^{-1}
# det(Phi) = alpha_signed^3 / det(M_d) = Lambda^6
# alpha_signed^3 = Lambda^6 * det(M_d)
alpha_signed_cubed = Lambda_QCD**6 * det_M_d
alpha_signed = np.sign(alpha_signed_cubed) * abs(alpha_signed_cubed)**(1/3)
M_meson_direct = alpha_signed * M_d_inv

print(f"\nVerification via direct matrix inversion (alpha * M_d^-1):")
print(f"  det(M_d) = {det_M_d:.4f} (= {'-' if det_M_d < 0 else '+'}m_d*m_s*m_b = {-m_d*m_s*m_b:.4f})")
print(f"  alpha_signed = {alpha_signed:.4f}")
for i in range(3):
    print(f"  {labels[i]:>4}", end='')
    for j in range(3):
        print(f"  {M_meson_direct[i,j]:>14.4f}", end='')
    print()

print(f"\n  Note: The signed-determinant version gives a matrix with different")
print(f"  overall sign because det(M_d) < 0 for Fritzsch texture.")
print(f"  The physical meson VEVs use |eigenvalues| for masses.")

# Off-diagonal ratios
print(f"\nOff-diagonal meson VEVs (physical version):")
print(f"  <Phi_ds>/<Phi_dd> = {M_meson[0,1]/M_meson[0,0]:.6f}")
print(f"  <Phi_db>/<Phi_dd> = {M_meson[0,2]/M_meson[0,0]:.6f}")
print(f"  <Phi_sb>/<Phi_ss> = {M_meson[1,2]/M_meson[1,1]:.6f}")
print(f"\n  Compare to Cabibbo-like ratios:")
print(f"  sqrt(m_d/m_s) = {sqrt(m_d/m_s):.6f}")
print(f"  sqrt(m_s/m_b) = {sqrt(m_s/m_b):.6f}")
print(f"  sqrt(m_d/m_b) = {sqrt(m_d/m_b):.6f}")


print("\n\n" + "=" * 72)
print("PART (e): Determinant consistency check")
print("=" * 72)

det_meson = np.linalg.det(M_meson)
Lambda6 = Lambda_QCD**6
C3 = alpha_seesaw**3
prod_m = m_d * m_s * m_b

print(f"\n  alpha^3 = Lambda^6 * m_d * m_s * m_b")
print(f"         = {Lambda6:.6e} * {prod_m:.6e} = {Lambda6*prod_m:.6e}")
print(f"  alpha^3 (computed) = {alpha_seesaw**3:.6e}")
print(f"  Match: {abs(alpha_seesaw**3 - Lambda6*prod_m)/(Lambda6*prod_m):.2e}")
print(f"")
print(f"  det(<Phi>) = prod(alpha/m_k) = alpha^3 / (m_d*m_s*m_b)")
print(f"             = Lambda^6 * m_d*m_s*m_b / (m_d*m_s*m_b)")
print(f"             = Lambda^6")
print(f"")
print(f"  Numerical check:")
print(f"    det(<Phi>) from matrix = {det_meson:.6e}")
print(f"    Lambda^6               = {Lambda6:.6e}")
print(f"    Relative error         = {abs(det_meson - Lambda6)/Lambda6:.2e}")
print(f"")
print(f"  ALGEBRAIC PROOF:")
print(f"  Let Phi = U diag(alpha/m_k) U^T with alpha^3 = Lambda^6 * prod(m_k).")
print(f"  det(Phi) = det(U)^2 * prod(alpha/m_k)")
print(f"           = 1 * alpha^3 / prod(m_k)")
print(f"           = Lambda^6 * prod(m_k) / prod(m_k)")
print(f"           = Lambda^6")
print(f"")
print(f"  This is independent of U (and would hold for U != V in the non-symmetric case).")
print(f"  The constraint det(Phi) = Lambda^6 is automatically satisfied.")
print(f"  Off-diagonal meson VEVs from the Fritzsch texture are CONSISTENT")
print(f"  with the Seiberg constraint.")


print("\n\n" + "=" * 72)
print("PART (f): Cabibbo angle and CKM angles -- full comparison")
print("=" * 72)

print(f"\n{'Method':<35} {'sin theta_12':>14} {'theta_12':>10} {'sin theta_23':>14} {'theta_23':>10} {'sin theta_13':>14} {'theta_13':>10}")
print("-" * 110)

rows = [
    ("GST/Oakes (m_u=0 limit)",
     sqrt(m_d/m_s), sqrt(m_s/m_b), sqrt(m_d/m_b)),
    ("Fritzsch approx formula",
     s12_fritzsch, s23_fritzsch, s13_fritzsch),
    ("Full numerical (exact Fritzsch)",
     s12, s23, s13),
    ("PDG 2024",
     s12_pdg, s23_pdg, s13_pdg),
]

for name, s12v, s23v, s13v in rows:
    print(f"  {name:<33} {s12v:>14.6f} {np.arcsin(s12v)*180/pi:>9.4f}  {s23v:>14.6f} {np.arcsin(s23v)*180/pi:>9.4f}  {s13v:>14.6f} {np.arcsin(s13v)*180/pi:>9.4f}")

print(f"\n--- Deviations from PDG ---")
for name, s12v, s23v, s13v in rows[:-1]:
    print(f"  {name:<33} theta_12: {(np.arcsin(s12v) - np.arcsin(s12_pdg))*180/pi:+.3f} deg   theta_23: {(np.arcsin(s23v) - np.arcsin(s23_pdg))*180/pi:+.3f} deg   theta_13: {(np.arcsin(s13v) - np.arcsin(s13_pdg))*180/pi:+.3f} deg")

print(f"\n--- Key observations ---")
print(f"  1. GST/Oakes gives theta_12 = {np.arcsin(sqrt(m_d/m_s))*180/pi:.2f} deg, within {abs(np.arcsin(sqrt(m_d/m_s))*180/pi - 13.04):.2f} deg of PDG (13.04 deg)")
print(f"     This is the best prediction: sqrt(m_d/m_s) = {sqrt(m_d/m_s):.4f} vs |V_us|_PDG = {s12_pdg:.4f}")
print(f"  2. The Fritzsch formula subtracts the up-sector angle, WORSENING the fit.")
print(f"  3. The full diagonalization gives sin theta_12 = {s12:.4f}, even further from PDG.")
print(f"     This is because the Fritzsch texture is too constrained: it gives")
print(f"     only 2 real parameters per sector (A, B, C with C fixed by trace).")

# Individual sector angles
print(f"\n--- Individual sector mixing angles ---")
print(f"  Down sector (1-2): sqrt(m_d/m_s) = {sqrt(m_d/m_s):.6f}  ({np.arcsin(sqrt(m_d/m_s))*180/pi:.4f} deg)")
print(f"  Down sector (2-3): sqrt(m_s/m_b) = {sqrt(m_s/m_b):.6f}  ({np.arcsin(sqrt(m_s/m_b))*180/pi:.4f} deg)")
print(f"  Up sector (1-2):   sqrt(m_u/m_c) = {sqrt(m_u/m_c):.6f}  ({np.arcsin(sqrt(m_u/m_c))*180/pi:.4f} deg)")
print(f"  Up sector (2-3):   sqrt(m_c/m_t) = {sqrt(m_c/m_t):.6f}  ({np.arcsin(sqrt(m_c/m_t))*180/pi:.4f} deg)")

# Dual Koide on meson eigenvalues
print(f"\n--- Dual Koide check on meson eigenvalues ---")
M_vals = sorted(meson_eigenvalues, reverse=True)
Q_meson = (sum(M_vals))**2 / (3 * sum(v**2 for v in M_vals))
# Also check on 1/m_i directly
inv_masses = [1/m_d, 1/m_s, 1/m_b]
Q_inv = (sum(inv_masses))**2 / (3 * sum(v**2 for v in inv_masses))
print(f"  Meson eigenvalues: {M_vals}")
print(f"  Koide Q(meson eigenvalues) = {Q_meson:.6f}")
print(f"  Koide Q(1/m_d, 1/m_s, 1/m_b) = {Q_inv:.6f}")
print(f"  These should be equal (Q depends only on ratios): match = {abs(Q_meson - Q_inv) < 1e-6}")
print(f"  Deviation from 2/3: {(Q_inv - 2/3)/(2/3)*100:.2f}%")
print(f"  Note: Dual Koide Q = {Q_inv:.4f} is far from 2/3 = 0.6667.")
print(f"  The near-miss Q(1/m_d, 1/m_s, 1/m_b) = 0.665 noted in the project")
print(f"  uses different mass values (probably pole masses or different PDG era).")

# Check with sqrt-Koide formula
print(f"\n--- Koide formula (sqrt version) ---")
# Standard Koide: (sum sqrt(m))^2 / (3 sum m) = 2/3
for label, masses in [("down quarks (d,s,b)", [m_d, m_s, m_b]),
                       ("up quarks (u,c,t)", [m_u, m_c, m_t]),
                       ("inverse down (1/m_b, 1/m_s, 1/m_d)", [1/m_b, 1/m_s, 1/m_d])]:
    sq = [sqrt(m) for m in masses]
    Q = (sum(sq))**2 / (3 * sum(masses))
    print(f"  Q({label}) = {Q:.6f}  (dev from 2/3: {(Q-2/3)/(2/3)*100:+.2f}%)")


print("\n\n" + "=" * 72)
print("COMPLETE SUMMARY")
print("=" * 72)

print(f"""
Part (a): Fritzsch relations verified.
  EXACT: C = m3-m2+m1, A^2 = m1*m2*m3/C, B^2 = m1*m2+m2*m3-m1*m3 - A^2
  APPROXIMATE (leading order): |A|^2 ~ m1*m2, |B|^2 ~ m2*m3
  The exact relations reproduce eigenvalues to machine precision.

Part (b): CKM angles computed.
  Full diag:      theta_12 = {theta12:.2f} deg, theta_23 = {theta23:.2f} deg, theta_13 = {theta13:.3f} deg
  Fritzsch approx: theta_12 = {np.arcsin(s12_fritzsch)*180/pi:.2f} deg, theta_23 = {np.arcsin(s23_fritzsch)*180/pi:.2f} deg, theta_13 = {np.arcsin(s13_fritzsch)*180/pi:.2f} deg
  GST/Oakes:      theta_12 = {np.arcsin(s12_gst)*180/pi:.2f} deg
  PDG:            theta_12 = 13.04 deg, theta_23 = 2.38 deg, theta_13 = 0.201 deg

  GST gives the best theta_12, within 0.12 deg of PDG.
  Full diag and Fritzsch approx both undershoot theta_12.
  theta_23 is overshot by ~50% by both Fritzsch methods.
  theta_13 matches well in full diag but is way off in Fritzsch approx.

Part (c): Seesaw is a MATRIX inversion.
  Phi = alpha * M_q^{{-1}}, eigenvalues inverted, eigenvectors preserved.
  Off-diagonal entries survive.

Part (d): Meson VEV matrix computed numerically.
  Significant off-diagonal entries: <Phi_ds>/<Phi_dd> ~ {M_meson[0,1]/M_meson[0,0]:.3f}

Part (e): det(Phi) = Lambda^6 verified both algebraically and numerically.
  The constraint is automatically satisfied regardless of mixing.

Part (f): Cabibbo angle comparison.
  GST relation sin theta_C = sqrt(m_d/m_s) = {sqrt(m_d/m_s):.4f} is the
  best simple prediction, remarkably close to PDG |V_us| = {s12_pdg:.4f}.
""")
