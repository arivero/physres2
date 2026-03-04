#!/usr/bin/env python3
"""
Parts 4-5: Near-pole analysis + summary + final plots.
Fix: use % formatting to avoid {0bar0} conflict with .format()
Also: investigate STr[M²] = 0 vs expected 2y.
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize_scalar
from scipy.integrate import quad

# Core functions
def fermion_masses_sq(t):
    s = np.sqrt(t**2 + 1)
    return np.array([0.0, (s - t)**2, (s + t)**2])

def scalar_masses_sq(t, y):
    Mplus = np.array([[4*t**2 + 1 - 2*y, 2*t],
                       [2*t,               1  ]])
    Mminus = np.array([[4*t**2 + 1 + 2*y, 2*t],
                        [2*t,               1  ]])
    eig_p = np.linalg.eigvalsh(Mplus)
    eig_m = np.linalg.eigvalsh(Mminus)
    return np.concatenate([[0.0, 0.0], eig_p, eig_m])

def CW_contribution(msq_arr, n_dof_arr, mu_sq=1.0):
    V = 0.0
    for msq, ndof in zip(msq_arr, n_dof_arr):
        if msq > 1e-30:
            V += ndof * msq**2 * (np.log(msq / mu_sq) - 1.5)
    return V / (64 * np.pi**2)

def compute_V_CW(t, y, mu_sq=1.0):
    smsq = scalar_masses_sq(t, y)
    fmsq = fermion_masses_sq(t)
    all_msq = np.concatenate([smsq, fmsq])
    all_ndof = np.concatenate([np.ones(6), -2*np.ones(3)])
    return CW_contribution(all_msq, all_ndof, mu_sq)

def K00(t, c=-1.0/12):
    return 1 + 4*c*t**2

def V_tree(t, y, c=-1.0/12):
    k = K00(t, c)
    if k <= 0:
        return np.inf
    return y**2 / k

t_pole = np.sqrt(3)
y = 0.1

# ============================================================
# STr[M²] check
# ============================================================
print("STr[M^2] investigation:")
print("  For O'Raifeartaigh: STr[M^2] = Tr[m^2_scalar] - 2*Tr[m^2_fermion]")
for tt in [0.0, 0.5, 1.0]:
    smsq = scalar_masses_sq(tt, y)
    fmsq = fermion_masses_sq(tt)
    str2 = np.sum(smsq) - 2*np.sum(fmsq)
    print(f"  t={tt}: sum(scalar m^2) = {np.sum(smsq):.6f}, 2*sum(fermion m^2) = {2*np.sum(fmsq):.6f}, STr = {str2:.6f}")
    print(f"    scalar eigenvalues: {np.sort(smsq)}")
    print(f"    fermion eigenvalues: {np.sort(fmsq)}")

# STr[M^2] = 0 because the Phi_0 sector contributes zero to both.
# The (Phi_1, Phi_2) sector: 4 real scalars with eigenvalues from M+,M-
# 2 Weyl fermions with eigenvalues (s-t)^2, (s+t)^2
# 
# Tr[scalar m^2 in (1,2) sector] = Tr(M+) + Tr(M-) = 2*Tr(M_F^2) = 2*(4t^2+2) 
# 2*Tr[fermion m^2] = 2*((s-t)^2 + (s+t)^2) = 2*(2(t^2+1)) = 4(t^2+1) = 4t^2+4
# But Tr(M_F^2) = 4t^2+1+1 = 4t^2+2, so Tr_scalar = 2*(4t^2+2) = 8t^2+4
# 2*Tr_fermion = 4t^2+4
# STr = 8t^2+4 - (4t^2+4) = 4t^2 ≠ 0 ??

# Let me recheck
print("\n  Manual check at t=1, y=0.1:")
t_test = 1.0
# M_F^2 eigenvalues in (1,2) sector: those of [[5,2],[2,1]]
MF2 = np.array([[5, 2],[2, 1]])
eig_f = np.linalg.eigvalsh(MF2)
print(f"  M_F^2 eigenvalues: {eig_f}")
print(f"  Fermion m^2 (Weyl): {fermion_masses_sq(t_test)}")
# Note: fermion masses are eigenvalues of W_ij, not W^dag W!
# The fermion mass-squared matrix is |W_ij|^2 = W^dag W
# But the eigenvalues of W_ij are the fermion masses (not mass-squared)
# Eigenvalues of W_ij: 0, and roots of λ^2 - 2t λ - 1 = 0 → λ = t ± sqrt(t^2+1)
# Mass-squared: (t + sqrt(t^2+1))^2, (t - sqrt(t^2+1))^2 = (sqrt(t^2+1) - t)^2

# The scalar mass matrix eigenvalues:
# For real parts: eigenvalues of M_F^2 + B = [[5-0.2, 2],[2, 1]]
M_plus = np.array([[5-0.2, 2],[2, 1]])
M_minus = np.array([[5+0.2, 2],[2, 1]])
eig_p = np.linalg.eigvalsh(M_plus)
eig_m = np.linalg.eigvalsh(M_minus)
print(f"  Scalar M+ eigenvalues: {eig_p}")
print(f"  Scalar M- eigenvalues: {eig_m}")

# Wait — the scalar mass matrix in the (Phi_1, Phi_2) sector is NOT M_F^2 ± B.
# It's more nuanced. The 4x4 real scalar mass matrix is:
# In the basis (phi_1, phi_1*, phi_2, phi_2*) or equivalently (Re phi_1, Im phi_1, Re phi_2, Im phi_2):
#
# The mass-squared matrix for complex scalar fields is:
# ( (W^dag W)_ij + B_ij     0              )     (for the "holomorphic" mass)
# ( 0                       (W^dag W)_ij - B_ij )
# 
# Wait, I need to be more careful. Let me use the standard result.
# 
# For the O'Raifeartaigh model, the scalar mass-squared eigenvalues are well known:
# In the (Phi_1, Phi_2) sector, 4 real scalars with masses:
# m^2 = (t^2 + 1) ± sqrt(4t^2 + y^2)   (2 eigenvalues, each doubly degenerate? No.)
#
# Actually the standard reference (e.g., Shih 0703196) gives:
# The (Phi_1, Phi_2) boson mass eigenvalues from the 4x4 matrix are:
# m^2_{1,2} = t^2 + 1/2 ± sqrt(t^2 + 1/4 + y)    [with appropriate scaling]
# Hmm, let me just compute the full 4x4 matrix.

# The scalar potential at the vacuum Phi_1 = Phi_2 = 0:
# V = |W_0|^2 + |W_1|^2 + |W_2|^2
# = |f + g Phi_1^2|^2 + |m Phi_2 + 2g v Phi_1|^2 + |m Phi_1|^2
#
# Expanding to quadratic order in (Phi_1, Phi_2):
# |f + g Phi_1^2|^2 ≈ f^2 + f g (Phi_1^2 + Phi_1*^2) + ...
# The f g Phi_1^2 term contributes to the (Phi_1, Phi_1) and (Phi_1*, Phi_1*) entries
# = f g (Phi_1^2 + c.c.) → this is the B-term

# |m Phi_2 + 2gv Phi_1|^2 = m^2|Phi_2|^2 + 4g^2v^2|Phi_1|^2 + 2mgv(Phi_1*Phi_2 + c.c.)
# |m Phi_1|^2 = m^2|Phi_1|^2
# 
# In m=1 units with t=gv:
# Diagonal: |Phi_1|^2 has coeff (4t^2 + 1), |Phi_2|^2 has coeff 1
# Off-diagonal: Phi_1* Phi_2 has coeff 2t
# B-term: Phi_1^2 has coeff y (= gf = y in m=1 units, since f = y/g and we set g=1... wait)
# 
# If g=1, m=1: t = v, y = f, so gf = f = y. And W_{110} = 2g = 2.
# B_{11} = W_{110} F_0* = 2·(-f) = -2y
# For the REAL scalar fields, this gives a splitting.
# The Phi_0 sector has F^0 = f (nonzero), but the Phi_0 mass matrix 
# at tree level vanishes (W_00 = W_01 = W_02 = 0 → no mass for Phi_0).

# OK let me just compute the full mass matrix numerically as a check.
# Using the 4x4 hermitian matrix approach for (Phi_1, Phi_2, Phi_1*, Phi_2*):

# The mass matrix M^2 in the (Phi, Phi*) basis:
# M^2 = ( A   B  )
#        ( B*  A* )
# where A_ij = sum_k (W_ki)* W_kj and B_ij = sum_k W_ijk F^{k*}
# 
# Here A = W^dag W restricted to (1,2):
# A = ( 4t^2+1,  2t )
#     ( 2t,      1  )
#
# B = ( -2y,  0 )
#     ( 0,    0 )
# (since B_11 = W_110 * F^{0*} = 2*(-y) = -2y, all others 0)
#
# The 4x4:
# ( 4t^2+1,  2t,    -2y,  0 )
# ( 2t,      1,      0,   0 )
# ( -2y,     0,  4t^2+1,  2t )  [conjugate block]
# ( 0,       0,      2t,  1  )

t_test = 1.0
y_test = 0.1
M4 = np.array([
    [4*t_test**2+1,  2*t_test,     -2*y_test,  0],
    [2*t_test,       1,             0,          0],
    [-2*y_test,      0,     4*t_test**2+1,  2*t_test],
    [0,              0,             2*t_test,   1]
])
eig4 = np.sort(np.linalg.eigvalsh(M4))
print(f"\n  4x4 matrix eigenvalues: {eig4}")

# Compare with my scalar_masses_sq function (excluding Phi_0 zeros):
smsq = scalar_masses_sq(t_test, y_test)
smsq_12 = np.sort(smsq[2:])  # skip the two Phi_0 zeros
print(f"  My function (1,2 sector): {smsq_12}")

# They should match!
# The 2x2 block approach: eigenvalues of (A ± B)
# A+B = ( 4t^2+1-2y, 2t )   and   A-B = ( 4t^2+1+2y, 2t )
#        ( 2t,        1  )            ( 2t,          1 )
# These give eigenvalues for the "holomorphic" and "anti-holomorphic" parts.
# But the 4x4 matrix mixes them...

# Actually for REAL fields, we diagonalize in (Re, Im) basis:
# Phi = (Re + i Im)/sqrt(2)
# The unitary transformation U = (1  1; -i  i)/sqrt(2) maps (Phi, Phi*) → (Re, Im)/sqrt(2)
# Under this, A+B → real part mass matrix, A-B → imaginary part mass matrix? Not quite.
# 
# For the 4x4 above, the eigenvalues of the FULL matrix give the 4 physical masses.
# Let me check if they agree with my A±B approach.

Amat = np.array([[4*t_test**2+1, 2*t_test],[2*t_test, 1]])
Bmat = np.array([[-2*y_test, 0],[0, 0]])
eig_ApB = np.sort(np.linalg.eigvalsh(Amat + Bmat))
eig_AmB = np.sort(np.linalg.eigvalsh(Amat - Bmat))
print(f"  A+B eigenvalues: {eig_ApB}")
print(f"  A-B eigenvalues: {eig_AmB}")
print(f"  Combined sorted: {np.sort(np.concatenate([eig_ApB, eig_AmB]))}")

# The 4x4 eigenvalues should equal the union of A+B and A-B eigenvalues
# because B is real symmetric → the 4x4 block-diagonalizes.
# Let me verify:
print(f"  4x4 sorted:      {eig4}")
print(f"  Match: {np.allclose(np.sort(np.concatenate([eig_ApB, eig_AmB])), eig4)}")

# Great. Now check STr[M^2]:
print(f"\n  STr[M^2] with correct counting:")
# 4 real scalars from (Phi_1,Phi_2), 2 from Phi_0 (zero)
# 3 Weyl fermions (0, and two from (1,2) sector)
# STr = sum(4 scalar m^2) - 2*sum(2 nonzero fermion m^2)
tr_scal = np.sum(eig4)
tr_ferm = 2 * np.sum(fermion_masses_sq(t_test)[1:])
print(f"  Tr[scalar m^2] = {tr_scal:.6f}")  
print(f"  2*Tr[fermion m^2] = {tr_ferm:.6f}")
print(f"  STr[M^2] = {tr_scal - tr_ferm:.6f}")
# Expected: STr[M^2] = 2*|F_0|^2 * (n_scalars_0 - n_fermions_0 - 1)... 
# Actually for O'Raifeartaigh: STr[M^2] = 0 in the minimal model
# because of the tree-level flat direction.
# The mass sum rule says STr[M^2] = 0 when SUSY is broken only by F-terms
# with renormalizable interactions and canonical Kahler.
# This is because the F-term contributions cancel exactly:
# sum_scalars m^2_i = Tr(A+B) + Tr(A-B) = 2 Tr(A)
# 2*sum_fermions m^2_i = 2 Tr(A)  [same matrix!]
# So STr = 0. ✓

print("\n  STr[M^2] = 0 confirmed. This is the standard mass sum rule")
print("  for F-term SUSY breaking with canonical Kahler potential.")

# Now: STr[M^4]
str4 = np.sum(eig4**2) - 2*np.sum(fermion_masses_sq(t_test)**2)
print(f"  STr[M^4] = {str4:.6f}")
# This should be 8y^2 from the B-term:
# STr[M^4] = Tr[(A+B)^2] + Tr[(A-B)^2] - 2*Tr[A^2]
# = 2Tr[A^2] + 2Tr[B^2] - 2Tr[A^2] = 2Tr[B^2]
# B = diag(-2y, 0) → B^2 = diag(4y^2, 0) → Tr[B^2] = 4y^2
# STr[M^4] = 8y^2
print(f"  Expected 8y^2 = {8*y_test**2:.6f}")

# ============================================================
# Now the MAIN conclusion plot
# ============================================================
print("\n" + "=" * 60)
print("MAIN RESULT: Why the CW potential increases monotonically")
print("=" * 60)

# V_CW = (1/64pi^2) STr[M^4(log M^2/mu^2 - 3/2)]
# Since STr[M^2] = 0 and STr[M^4] = 8y^2 (constant in t!),
# the CW potential's t-dependence comes from STr[M^4 log M^2]:
# V_CW(t) = (1/64pi^2) [STr[M^4 log M^2] - (3/2) STr[M^4]]
# = (1/64pi^2) [STr[M^4 log M^2] - 12y^2]
# 
# Since STr[M^4] is t-independent, the shape comes entirely from 
# STr[M^4 log M^2].

print("  The CW potential shape comes from STr[M^4 log M^2].")
print("  At large t: all masses grow ~ t^2, so log M^2 ~ 2 log t.")
print("  STr[M^4 log M^2] ~ 8y^2 * 2 log t → increasing.")
print("  This is why V_CW increases monotonically for all y.")

# Shih's result: the CW potential of the canonical O'Raifeartaigh model
# stabilizes the pseudo-modulus at v = 0. This is CORRECT.
# Reference: Shih hep-ph/0703196

# ============================================================
# FINAL COMPREHENSIVE PLOT
# ============================================================

fig, axes = plt.subplots(1, 3, figsize=(18, 5.5))

# Panel 1: CW potential alone
ax = axes[0]
t_arr = np.linspace(0.01, 3.5, 500)
for yy in [0.05, 0.1, 0.5, 1.0]:
    V = np.array([compute_V_CW(t, yy) for t in t_arr])
    V -= V[0]
    ax.plot(t_arr, V / yy**2, linewidth=1.5, label=f'y = {yy}')
ax.axvline(t_pole, color='k', linestyle='--', alpha=0.5, label='pole')
ax.set_xlabel('t = gv/m', fontsize=13)
ax.set_ylabel(r'$(V_{CW}(t) - V_{CW}(0))/y^2$', fontsize=13)
ax.set_title('CW potential (rescaled by y^2)', fontsize=14)
ax.legend(fontsize=10)
ax.grid(True, alpha=0.3)

# Panel 2: V_tree + V_CW 
ax = axes[1]
t_arr2 = np.linspace(0.01, t_pole - 0.02, 500)
for yy in [0.05, 0.1, 0.3]:
    Veff = np.array([V_tree(t, yy) + compute_V_CW(t, yy) for t in t_arr2])
    Veff -= Veff[0]
    ax.plot(t_arr2, Veff, linewidth=1.5, label=f'y = {yy}')
ax.axvline(t_pole, color='k', linestyle='--', alpha=0.5, label='pole')
ax.set_xlabel('t = gv/m', fontsize=13)
ax.set_ylabel(r'$V_{eff}(t) - V_{eff}(0)$', fontsize=13)
ax.set_title('Total potential with Kahler pole', fontsize=14)
ax.legend(fontsize=10)
ax.grid(True, alpha=0.3)

# Panel 3: Log-scale near pole showing divergences
ax = axes[2]
yy = 0.1
t_near = np.linspace(0.5, t_pole - 0.02, 500)
Vtr = np.array([V_tree(t, yy) for t in t_near])
Vcw = np.array([compute_V_CW(t, yy) for t in t_near])
ax.semilogy(t_near, Vtr, 'r-', linewidth=2, label=r'$V_{tree}$')
ax.semilogy(t_near, np.abs(Vcw), 'b-', linewidth=2, label=r'$|V_{CW}|$')
ax.semilogy(t_near, Vtr + Vcw, 'k-', linewidth=2, label=r'$V_{eff}$')
ax.axvline(t_pole, color='k', linestyle='--', alpha=0.5)
ax.set_xlabel('t = gv/m', fontsize=13)
ax.set_ylabel('V (log scale)', fontsize=13)
ax.set_title(f'Near-pole behavior (y = {yy})', fontsize=14)
ax.legend(fontsize=10)
ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('/home/codexssh/phys3/results/kahler_pole_potential_physics.png', dpi=150)
plt.close()
print("\nFinal plot saved: results/kahler_pole_potential_physics.png")

# ============================================================
# SUMMARY
# ============================================================
print("\n" + "=" * 70)
print("PART 5: COMPLETE SUMMARY")
print("=" * 70)

t_crit = 1.549193
sigma_max = np.sqrt(3)*np.pi/4

print(f"""
SUMMARY: Kahler pole stabilization of the O'Raifeartaigh pseudo-modulus
=========================================================================

Setup:
  W = f Phi_0 + m Phi_1 Phi_2 + g Phi_0 Phi_1^2
  K = |Phi_0|^2 - (1/12)|Phi_0|^4 / (m/g)^2 + |Phi_1|^2 + |Phi_2|^2
  
  Kahler metric: K_00bar = 1 - t^2/3,  where t = gv/m
  Pole at t = sqrt(3) = {t_pole:.6f}
  Proper distance to pole: sigma_max = sqrt(3) pi/4 = {sigma_max:.6f} (FINITE)

Part 1 results:
  The Coleman-Weinberg potential V_CW(t) is MONOTONICALLY INCREASING 
  from t = 0 for ALL values of y = gf/m^2.
  
  There is NO minimum at t ~ 0.49 or anywhere else.
  The pseudo-modulus is stabilized at v = 0 (the origin).
  This agrees with the standard result (Shih 2007).
  
  Analytic reason: STr[M^2] = 0 (mass sum rule), STr[M^4] = 8y^2 
  (t-independent). The t-dependence is entirely in STr[M^4 log M^2],
  which increases monotonically because all masses grow with t.

Part 2 results:
  V_eff = V_tree + V_CW = y^2/(1 - t^2/3) + V_CW(t)
  
  V_tree increases monotonically and diverges at the pole.
  V_CW increases monotonically.
  Both contributions push toward t = 0.
  
  The minimum is at t = 0 for ALL values of y. There is NO value 
  of y that moves the minimum to t = sqrt(3).
  The Kahler pole acts as a REPULSIVE WALL, not an attractive well.

Part 3 results:
  Kahler dressing replaces y -> y_eff = y / K_00bar in the B-term.
  Near the pole, y_eff diverges, which:
  (a) Makes V_CW_dressed diverge faster than V_tree (as 1/eps^4 vs 1/eps)
  (b) Creates a TACHYON at t = {t_crit:.4f} (K = 0.20, y_eff = 0.50)
      Distance from pole: {t_pole - t_crit:.4f}
  
  The model breaks down before the field reaches the pole.
  V_CW_dressed increases even faster than undressed → still no minimum.

Part 4 results:
  The pole is at FINITE proper distance sigma_max = {sigma_max:.4f}.
  V_tree -> infinity as 1/eps at the boundary.
  V_CW is finite and increasing at the boundary.
  The boundary is an impenetrable wall in the potential.
  
  The minimum AT the pole is impossible because V -> infinity.
  The minimum cannot be pushed arbitrarily close to the pole 
  because both V_tree and V_CW increase monotonically.

CONCLUSION:
  The Kahler pole stabilization FAILS because the pole creates a WALL,
  not a WELL. The tree-level Kahler pole V_tree = f^2/(1 - t^2/3) 
  diverges positively at the boundary, and the Coleman-Weinberg 
  potential increases monotonically in the same direction. There is 
  no competition between a rising V_tree and a falling V_CW — both 
  increase together.

  For stabilization at t = sqrt(3), one needs a mechanism that creates 
  a NEGATIVE contribution to the potential that grows near the pole. 
  Candidates: nonperturbative superpotential (instantons, strong dynamics), 
  or a different origin for the c = -1/12 coefficient that modifies the 
  mass spectrum rather than just the Kahler metric.

  The c = -1/12 tree-level Kahler pole is a BOUNDARY of field space, 
  not a vacuum. The pseudo-modulus sits at the origin v = 0, far from 
  the pole.
""")

