#!/usr/bin/env python3
"""
Supermultiplet restoration: track fermion-scalar pairs from SUSY limit
through SUSY breaking to the physical spectrum.

In the O'Raifeartaigh superpotential W = fX + m Phi1 Phi3 + g X Phi1^2,
the SUSY limit (f -> 0) makes each fermion eigenvalue degenerate with
a scalar eigenvalue. Turning f back on splits them.

This script:
1. Computes fermion and scalar mass spectra at general f
2. Identifies the degenerate pairs at f=0 (SUSY restored)
3. Tracks the splitting as f increases
4. Checks m_pi^2 - m_mu^2 = f_pi^2 against the B-term splitting
5. Identifies which meson pairs with which lepton
"""
import numpy as np
from numpy.linalg import eigvalsh

print("=" * 70)
print("SUPERMULTIPLET RESTORATION")
print("=" * 70)

# =====================================================================
# 1. O'Raifeartaigh mass matrices
# =====================================================================
# W = f X + m Phi1 Phi3 + g X Phi1^2
# Vacuum: <X> = v, <Phi_i> = 0, F_X = f
# t = gv/m, at Koide point t = sqrt(3)

m_param = 1.0  # normalized
g = 1.0

print("\n--- 1. Mass matrices at t = sqrt(3) ---\n")

t = np.sqrt(3)
v = t * m_param / g

# Fermion mass matrix (2x2 block in Phi1, Phi3 subspace)
# M_F = [[2gv, m], [m, 0]]
# Eigenvalues: m_pm = sqrt(g^2 v^2 + m^2) +/- gv
m_plus = np.sqrt(g**2 * v**2 + m_param**2) + g * v
m_minus = np.sqrt(g**2 * v**2 + m_param**2) - g * v

print(f"Fermion eigenvalues (normalized):")
print(f"  psi_X:  0  (goldstino)")
print(f"  psi_2:  0  (pseudo-modulus)")
print(f"  psi_-:  {m_minus:.6f}  = (2-sqrt3) = {2 - np.sqrt(3):.6f}")
print(f"  psi_+:  {m_plus:.6f}  = (2+sqrt3) = {2 + np.sqrt(3):.6f}")
print(f"  Ratio m+/m- = {m_plus/m_minus:.4f} = (2+sqrt3)^2 = {(2+np.sqrt(3))**2:.4f}")
print(f"  Product m+ * m- = {m_plus * m_minus:.6f} = m^2 = {m_param**2:.6f}")

# Scalar mass matrix
# The scalar mass^2 matrix includes:
# - F-term: |M_F|^2 (same eigenvalues squared)
# - B-term: F_X * d^2(dW/dX)/dPhi_i dPhi_j = f * 2g * delta_i1 delta_j1
# For the real components sigma (Re) and pi (Im):
# m^2_sigma = |M_F|^2 + B = |M_F|^2 + 2gf  (for Phi1 direction)
# m^2_pi = |M_F|^2 - B = |M_F|^2 - 2gf  (for Phi1 direction)
# Other directions: no B-term splitting

print("\n--- 2. Scalar mass^2 spectrum ---\n")
print("The B-term splits ONLY the Phi1 direction:")
print("  m^2_sigma(Phi1) = |M_F(Phi1)|^2 + 2gf")
print("  m^2_pi(Phi1)    = |M_F(Phi1)|^2 - 2gf")
print("  m^2(Phi2)       = 0  (pseudo-modulus, both sigma and pi)")
print("  m^2(Phi3)       = m^2  (no B-term)")
print("  m^2(X)          = 0 + 2gf for sigma_X, 0 - 2gf for pi_X")
print()
print("But the mass eigenstates MIX Phi1 and Phi3!")
print("Need to diagonalize the full scalar mass^2 matrix.")

# Full scalar mass^2 matrix
# In the basis (X, Phi1, Phi2, Phi3), the fermion mass matrix is:
# M_F = [[0, 0, 0, 0],
#         [0, 2gv, 0, m],
#         [0, 0, 0, 0],
#         [0, m, 0, 0]]
# (ignoring the off-diagonal X-Phi1 terms which vanish at <Phi1>=0)

# F-term scalar mass^2 = M_F^dag M_F (for holomorphic masses)
# In the (Phi1, Phi3) subspace:
MF_block = np.array([[2*g*v, m_param],
                      [m_param, 0]])
MFsq = MF_block.T @ MF_block  # = M_F^T M_F for real M_F

print("\n--- 3. Full scalar mass^2 in (Phi1, Phi3) subspace ---\n")
print(f"M_F^T M_F = ")
print(f"  [{MFsq[0,0]:.4f}, {MFsq[0,1]:.4f}]")
print(f"  [{MFsq[1,0]:.4f}, {MFsq[1,1]:.4f}]")

# B-term matrix (only in Phi1-Phi1 entry)
# B = d^2 W / dX dPhi_i = 2g delta_{i,1}
# So B-term contribution to scalar mass^2 = f * 2g in (1,1) position

# For sigma (real parts): m^2 = M_F^T M_F + f * B_matrix
# For pi (imaginary parts): m^2 = M_F^T M_F - f * B_matrix

print("\n--- 4. Spectrum as function of f (SUSY breaking parameter) ---\n")
print(f"{'f/m':>8s} | {'sigma1':>10s} {'sigma3':>10s} {'pi1':>10s} {'pi3':>10s} | "
      f"{'ferm_-':>10s} {'ferm_+':>10s} | {'split_low':>10s} {'split_high':>10s}")
print("-" * 110)

B_matrix = np.array([[2*g, 0],
                      [0, 0]])

for f_ratio in [0, 0.01, 0.05, 0.1, 0.2, 0.5, 1.0]:
    f_val = f_ratio * m_param

    # Sigma (real) mass^2
    Msq_sigma = MFsq + f_val * B_matrix
    eig_sigma = np.sort(eigvalsh(Msq_sigma))

    # Pi (imaginary) mass^2
    Msq_pi = MFsq - f_val * B_matrix
    eig_pi = np.sort(eigvalsh(Msq_pi))

    # Fermion masses (unchanged by f)
    ferm = np.sort([m_minus, m_plus])

    # Splitting: scalar mass - fermion mass for each pair
    # At f=0: sigma = pi = |fermion mass|, so mass^2 = fermion^2
    # Track which scalar eigenvalue is closest to which fermion

    sigma_masses = np.sqrt(np.maximum(eig_sigma, 0))
    pi_masses = np.sqrt(np.maximum(eig_pi, 0))

    split_low = eig_sigma[0] - ferm[0]**2  # m^2 splitting for lower pair
    split_high = eig_sigma[1] - ferm[1]**2

    print(f"{f_ratio:8.3f} | {sigma_masses[0]:10.6f} {sigma_masses[1]:10.6f} "
          f"{pi_masses[0]:10.6f} {pi_masses[1]:10.6f} | "
          f"{ferm[0]:10.6f} {ferm[1]:10.6f} | "
          f"{split_low:10.6f} {split_high:10.6f}")

# =====================================================================
# 5. SUSY limit: identify the pairs
# =====================================================================
print("\n--- 5. SUSY limit (f=0): degenerate pairs ---\n")

f_val = 0
Msq_sigma_susy = MFsq
eig_sigma_susy = np.sort(eigvalsh(Msq_sigma_susy))
sigma_susy = np.sqrt(np.maximum(eig_sigma_susy, 0))

Msq_pi_susy = MFsq  # same as sigma when f=0
eig_pi_susy = np.sort(eigvalsh(Msq_pi_susy))
pi_susy = np.sqrt(np.maximum(eig_pi_susy, 0))

print("At f=0:")
print(f"  Fermions: psi_- = {m_minus:.6f}, psi_+ = {m_plus:.6f}")
print(f"  Scalars (sigma): {sigma_susy[0]:.6f}, {sigma_susy[1]:.6f}")
print(f"  Scalars (pi):    {pi_susy[0]:.6f}, {pi_susy[1]:.6f}")
print()
print("Supermultiplet identification at SUSY point:")
print(f"  Multiplet 1: fermion {m_minus:.4f} <-> sigma {sigma_susy[0]:.4f} + pi {pi_susy[0]:.4f}")
print(f"  Multiplet 2: fermion {m_plus:.4f} <-> sigma {sigma_susy[1]:.4f} + pi {pi_susy[1]:.4f}")
print(f"  Multiplet X: goldstino <-> pseudo-modulus X (both massless)")
print(f"  Multiplet 2: psi_2 <-> Phi_2 (both massless, CW mass)")

# =====================================================================
# 6. Eigenvector analysis: which Phi is in which multiplet?
# =====================================================================
print("\n--- 6. Eigenvector decomposition ---\n")
print("The fermion mass eigenstates mix Phi1 and Phi3:")

# Diagonalize M_F in (Phi1, Phi3) subspace
evals_F, evecs_F = np.linalg.eigh(MF_block.T @ MF_block)
# But we want the eigenvectors of M_F itself
# M_F |v> = lambda |v>
# Since M_F is real symmetric... no it's not symmetric!
# M_F = [[2gv, m], [m, 0]]
# Eigenvalues of M_F: lambda = gv +/- sqrt(g^2v^2 + m^2)
# = t*m +/- sqrt(t^2+1)*m
# At t=sqrt(3): lambda_+ = (sqrt(3) + 2)m, lambda_- = (sqrt(3) - 2)m (negative!)
# So m_- = |lambda_-| = (2 - sqrt(3))m

# The eigenvectors of M_F:
# (2gv - lambda) v1 + m v3 = 0
# For lambda_+:
lam_plus = g*v + np.sqrt(g**2*v**2 + m_param**2)
lam_minus = g*v - np.sqrt(g**2*v**2 + m_param**2)
# Note: lam_minus is NEGATIVE, |lam_minus| = m_minus

# Eigenvector for lambda_+: (2gv - lam_+) v1 + m v3 = 0
# v3/v1 = -(2gv - lam_+)/m = (lam_+ - 2gv)/m
ratio_plus = (lam_plus - 2*g*v) / m_param
norm_plus = np.sqrt(1 + ratio_plus**2)
v1_plus = 1/norm_plus
v3_plus = ratio_plus/norm_plus

ratio_minus = (lam_minus - 2*g*v) / m_param
norm_minus = np.sqrt(1 + ratio_minus**2)
v1_minus = 1/norm_minus
v3_minus = ratio_minus/norm_minus

print(f"  psi_+ (mass {m_plus:.4f}): {v1_plus:.4f} * Phi1 + {v3_plus:.4f} * Phi3")
print(f"  psi_- (mass {m_minus:.4f}): {v1_minus:.4f} * Phi1 + {v3_minus:.4f} * Phi3")
print()

# Verify orthogonality
dot = v1_plus * v1_minus + v3_plus * v3_minus
print(f"  Orthogonality check: v_+ . v_- = {dot:.2e}")

# =====================================================================
# 7. Physical identification with mesons and leptons
# =====================================================================
print("\n--- 7. Physical spectrum (lepton sector) ---\n")

m_OR_lep = 470.76  # MeV, lepton sector O'R scale
m_e = 0.511
m_mu = 105.658
m_tau = 1776.86
f_pi = 92.2  # MeV

m_minus_phys = (2 - np.sqrt(3)) * m_OR_lep
m_plus_phys = (2 + np.sqrt(3)) * m_OR_lep

print(f"O'R scale: m_OR = {m_OR_lep:.2f} MeV")
print(f"O'R eigenvalues:")
print(f"  m_- = (2-sqrt3) * {m_OR_lep:.2f} = {m_minus_phys:.2f} MeV")
print(f"  m_+ = (2+sqrt3) * {m_OR_lep:.2f} = {m_plus_phys:.2f} MeV")
print(f"  m_+/m_- = {m_plus_phys/m_minus_phys:.4f}")
print()
print("Actual lepton masses:")
print(f"  m_e = {m_e:.3f} MeV")
print(f"  m_mu = {m_mu:.3f} MeV")
print(f"  m_tau = {m_tau:.2f} MeV")
print(f"  (m_e + m_mu + m_tau)/4 = {(m_e + m_mu + m_tau)/4:.2f} MeV = m_OR ✓")
print()

print("B-term splitting: m^2_scalar - m^2_fermion = +/- 2gf")
print(f"  If this equals f_pi^2 = {f_pi**2:.0f} MeV^2:")
print(f"  2gf = f_pi^2 → f = {f_pi**2/(2*g):.0f} MeV^2 (for g=1)")
print()
print("Check founding relation:")
m_pi = 139.57
print(f"  m_pi^2 - m_mu^2 = {m_pi**2 - m_mu**2:.0f} MeV^2")
print(f"  f_pi^2 = {f_pi**2:.0f} MeV^2")
print(f"  Ratio: {(m_pi**2 - m_mu**2)/f_pi**2:.4f}")
print()

# =====================================================================
# 8. Meson candidates for each supermultiplet
# =====================================================================
print("\n--- 8. Meson partner candidates ---\n")

print("In SUSY limit, each supermultiplet has:")
print("  fermion mass = scalar mass (degenerate)")
print()
print("When SUSY breaks (f > 0):")
print("  m^2_sigma = m^2_fermion + 2gf  (sigma = real part)")
print("  m^2_pi    = m^2_fermion - 2gf  (pi = imaginary part)")
print()
print("The pseudoscalar meson is the 'pi' component.")
print("For the founding relation: m^2_pi(meson) = m^2_fermion(lepton) + Delta")
print("where Delta depends on the B-term direction.")
print()

# If m_pi^2 - m_mu^2 = f_pi^2, then the pion is the sigma (heavier)
# component and the muon is the fermion, with splitting = f_pi^2.
# OR: the pion is the pi component with m^2_pi = m^2_mu + f_pi^2
# (the scalar gets HEAVIER, which means m^2_scalar > m^2_fermion,
# which is the sigma channel with +2gf)

print("Founding relation: m_pi^2 = m_mu^2 + f_pi^2")
print(f"  {m_pi**2:.0f} = {m_mu**2:.0f} + {f_pi**2:.0f}")
print(f"  {m_pi**2:.0f} ≈ {m_mu**2 + f_pi**2:.0f} (match: {abs(m_pi**2 - m_mu**2 - f_pi**2)/(m_pi**2)*100:.1f}%)")
print()
print("This means the pion is the SCALAR partner of the muon,")
print("split upward by the SUSY-breaking B-term.")
print()

# Cross-generation pairings (ts), (ub), (cd)
print("Cross-generation supermultiplet structure:")
print("  The quark pairings (t,s), (u,b), (c,d) mean that")
print("  within each supermultiplet, one up-type and one down-type")
print("  quark from DIFFERENT standard generations are paired.")
print()
print("  Charged mesons with Q = -1 (lepton partners):")
print("  From (t,s): s t-bar → Q = -1/3 - 2/3 = -1")
print("  From (u,b): b u-bar → Q = -1/3 - 2/3 = -1  (= B- meson)")
print("  From (c,d): d c-bar → Q = -1/3 - 2/3 = -1  (= D- meson)")
print()

m_Bmeson = 5279.34  # MeV
m_Dmeson = 1869.66  # MeV

print("If m^2_meson = m^2_lepton + f_pi^2:")
for name, m_mes, m_lep_name in [("pi (ud-bar)", m_pi, "mu"),
                                   ("D- (dc-bar)", m_Dmeson, "?"),
                                   ("B- (bu-bar)", m_Bmeson, "?")]:
    m_lep_pred = np.sqrt(max(m_mes**2 - f_pi**2, 0))
    print(f"  {name}: m_meson = {m_mes:.1f} MeV → m_lepton = {m_lep_pred:.1f} MeV ({m_lep_name})")

print()
print("If m^2_meson = m^2_lepton + C (non-universal splitting):")
print("  pi → mu: C = m_pi^2 - m_mu^2 = {:.0f} MeV^2 = f_pi^2".format(m_pi**2 - m_mu**2))
print("  D → tau: C = m_D^2 - m_tau^2 = {:.0f} MeV^2".format(m_Dmeson**2 - m_tau**2))
print("  B → ???: C = m_B^2 - m_???^2")
print()
print(f"  D^2 - tau^2 = {m_Dmeson**2 - m_tau**2:.0f} MeV^2")
print(f"  sqrt(D^2 - tau^2) = {np.sqrt(m_Dmeson**2 - m_tau**2):.1f} MeV")
print(f"  Compare: f_pi = {f_pi} MeV, so this is NOT f_pi^2!")
print()

# =====================================================================
# 9. The real question: what IS the SUSY-breaking pattern?
# =====================================================================
print("\n--- 9. Non-universal SUSY breaking ---\n")
print("The B-term is NOT universal: it acts only on the Phi1 direction.")
print("The mass eigenstates are mixtures of Phi1 and Phi3.")
print("So the splitting depends on the Phi1 component of each eigenstate.")
print()
print("Splitting of eigenstate |psi_k> = a_k |Phi1> + b_k |Phi3>:")
print("  Delta m^2_k = 2gf * |a_k|^2")
print()
print(f"  psi_+ has Phi1 component |a_+|^2 = {v1_plus**2:.6f}")
print(f"  psi_- has Phi1 component |a_-|^2 = {v1_minus**2:.6f}")
print()
print("So the heavier eigenstate (tau) gets a SMALLER B-term splitting")
print("and the lighter eigenstate (muon) gets a LARGER B-term splitting.")
print()

# At t = sqrt(3), the eigenvectors are:
# psi_+ ∝ (1, ratio_+) where ratio_+ = (lam_+ - 2gv)/m
# The Phi1 content of each eigenstate:
phi1_content_plus = v1_plus**2
phi1_content_minus = v1_minus**2
print(f"  B-term splitting for psi_+ (tau pair): 2gf * {phi1_content_plus:.4f}")
print(f"  B-term splitting for psi_- (mu pair):  2gf * {phi1_content_minus:.4f}")
print(f"  Ratio of splittings: {phi1_content_minus/phi1_content_plus:.4f}")
print()

# Physical numbers
# If 2gf * |a_-|^2 = f_pi^2 = 8501 (for the muon-pion pair):
two_gf = f_pi**2 / phi1_content_minus
print(f"  From mu-pi pair: 2gf = {two_gf:.0f} MeV^2")
print(f"  Then tau-meson splitting = 2gf * {phi1_content_plus:.4f} = {two_gf * phi1_content_plus:.0f} MeV^2")
m_tau_partner = np.sqrt(m_tau**2 + two_gf * phi1_content_plus)
print(f"  m(tau partner) = sqrt(m_tau^2 + {two_gf * phi1_content_plus:.0f}) = {m_tau_partner:.1f} MeV")
print()

# What meson has this mass?
print(f"  Predicted tau meson partner mass: {m_tau_partner:.1f} MeV")
print(f"  D meson: {m_Dmeson:.1f} MeV")
print(f"  D_s meson: 1968.3 MeV")
print(f"  Deviation from D: {abs(m_tau_partner - m_Dmeson)/m_Dmeson*100:.1f}%")
print(f"  Deviation from D_s: {abs(m_tau_partner - 1968.3)/1968.3*100:.1f}%")

# =====================================================================
# 10. Summary
# =====================================================================
print("\n" + "=" * 70)
print("SUMMARY")
print("=" * 70)
print()
print("In the O'Raifeartaigh model at t = sqrt(3):")
print()
print("SUSY LIMIT (f=0):")
print("  Multiplet A: psi_- (mass m_-) ↔ scalar_- (mass m_-) [degenerate]")
print("  Multiplet B: psi_+ (mass m_+) ↔ scalar_+ (mass m_+) [degenerate]")
print("  Multiplet X: goldstino ↔ pseudo-modulus X [both massless]")
print("  Multiplet 2: psi_2 ↔ Phi_2 [both massless, CW mass]")
print()
print("BROKEN SUSY (f > 0):")
print("  The B-term splits scalars from fermions NON-UNIFORMLY:")
print(f"  Multiplet A (mu-pion): Delta m^2 = 2gf × {phi1_content_minus:.4f}")
print(f"  Multiplet B (tau-???): Delta m^2 = 2gf × {phi1_content_plus:.4f}")
print()
print(f"  Setting Delta_A = f_pi^2 (founding relation):")
print(f"  → Tau meson partner: {m_tau_partner:.0f} MeV")
print()
print("OPEN QUESTIONS:")
print("  1. Which specific meson is the tau partner?")
print("  2. What determines the electron mass (bloom of zero eigenvalue)?")
print("  3. How do the cross-generation pairings (ts),(ub),(cd) emerge?")
print("  4. Is the Phi_2 pseudo-modulus = neutrino?")

# =====================================================================
# 11. Bloom correction: electron mass from zero eigenvalue
# =====================================================================
print("\n" + "=" * 70)
print("BLOOM CORRECTION: ELECTRON FROM ZERO EIGENVALUE")
print("=" * 70)
print()

# In the O'R model, the goldstino (psi_X) and pseudo-modulus (psi_2)
# are both massless at tree level. The bloom (R-symmetry breaking
# perturbation epsilon*Phi_0*Phi_2) gives mass to psi_X.
#
# At leading order in epsilon, the mass is:
#   m_electron ~ 3^{5/4} |epsilon/m| * m  (from bloom direction calculation)
#
# We can use the electron mass to fix epsilon:
m_e = 0.511  # MeV
m_mu = 105.658  # MeV
m_tau = 1776.86  # MeV

# The O'R scale m is related to the strange quark mass:
# m_s = z0^2 * m_- = z0^2 * (2 - sqrt3) * m
# With z0-doubling: z0(bloom) = 2 * z0(seed)
# From the one-parameter spectrum: m_u + m_d = S = 7.05 MeV
# The isospin seed gives m_s = (2+sqrt3)^2 * S = 98.2 MeV
# Then m_- * z0^2 = m_s gives the O'R scale

# But for the lepton sector, the key quantity is the RATIO:
# m_e / m_mu (or m_e / m_tau)
# The bloom gives m_e from the zero eigenvalue, so:
# m_e/m_mu = epsilon_eff / m_-

# From the mass chain: m_e = 0.511 MeV
# The O'R eigenvalues at t=sqrt(3): (0, 2-sqrt3, 2+sqrt3)
# After bloom: (epsilon_bloom, 2-sqrt3, 2+sqrt3)
# The electron comes from epsilon_bloom

print("The O'R eigenvalues before bloom: (0, m_-, m_+)")
print(f"  m_-/m_+ = {m_minus/m_plus:.6f} = (2-sqrt3)^2 = {(2-np.sqrt(3))**2:.6f}")
print()

# What epsilon_bloom/m gives m_e/m_mu?
# In the Koide parametrization: z_k = z0 * sqrt(1 + sqrt(2) * cos(delta + 2pi*k/3))
# The seed has delta = 3pi/4 where cos(delta) = -1/sqrt(2), forcing z_0 = 0
# Bloom rotates delta away from 3pi/4 by a small angle dd
# At leading order: z_0^2 ~ 2*dd * z_0_coeff
# This gives m_e = z_0^2 proportional to dd

# For leptons, seed at delta_seed = 3*pi/4:
delta_seed = 3 * np.pi / 4
z0 = 1.0  # normalized

# Physical lepton Koide relation:
# Σm_k / (Σ√m_k)² = 2/3
sum_sqrt = np.sqrt(m_e) + np.sqrt(m_mu) + np.sqrt(m_tau)
sum_m = m_e + m_mu + m_tau
Q_lep = sum_m / sum_sqrt**2
print(f"Lepton Koide Q = {Q_lep:.6f} (should be ~2/3 = {2/3:.6f})")

# The v0 and delta for leptons:
# Parametrization: m_k = v0^2 (1 + sqrt(2) cos(delta + 2pi*k/3))
# where k labels need NOT be k=0→e, k=1→mu, k=2→tau.
# v0 = sum(sqrt(m)) / 3
v0_lep = sum_sqrt / 3
print(f"v0_lepton = {v0_lep:.4f} MeV^{{1/2}}")
print(f"v0^2 = {v0_lep**2:.4f} MeV")
print()

# Extract delta using the Cartan angle formula.
# The Cartan vector r = sum sqrt(m_k) w_k where w_k are SU(3) weights.
# Convention: w_0 = (0,1), w_1 = (-sqrt3/2, -1/2), w_2 = (sqrt3/2, -1/2)
# The key relation: theta = pi/6 - delta (exact identity)
#
# But the ASSIGNMENT of (e, mu, tau) to (k=0, k=1, k=2) matters.
# At the seed delta=3pi/4, the zero eigenvalue is at k=0.
# After bloom, the smallest eigenvalue STAYS at k=0 (continuous deformation).
# So: k=0 → e (lightest), k=1 → ?, k=2 → ?
#
# From cos(delta + 2pi*k/3) we need the three values to be ordered correctly.
# Let's find delta by numerical fit.

from scipy.optimize import minimize_scalar

masses_lep = np.array([m_e, m_mu, m_tau])

# IMPORTANT: The Koide parametrization is:
#   sqrt(m_k) = v_0 * (1 + sqrt(2) * cos(delta + 2*pi*k/3))
# NOT m_k = v_0^2 * (1 + sqrt(2) * cos(...)).
# The squaring is OUTSIDE: m_k = [v_0 * (1 + sqrt(2)*cos(...))]^2

def koide_residual(delta):
    """Find delta that reproduces lepton masses (up to permutation of k)."""
    vals = [(v0_lep * (1 + np.sqrt(2) * np.cos(delta + 2*np.pi*k/3)))**2 for k in range(3)]
    vals_sorted = sorted(vals)
    masses_sorted = sorted(masses_lep)
    return sum((v - m)**2 for v, m in zip(vals_sorted, masses_sorted))

# Scan for minimum
best_delta = None
best_res = 1e30
for trial in np.linspace(0, 2*np.pi, 1000):
    r = koide_residual(trial)
    if r < best_res:
        best_res = r
        best_delta = trial

# Refine
from scipy.optimize import minimize
result = minimize(koide_residual, best_delta, method='Nelder-Mead')
delta_lep = result.x[0] % (2*np.pi)

print(f"Lepton Koide phase delta = {delta_lep:.6f} rad = {delta_lep*180/np.pi:.2f}°")
print(f"  Residual = {result.fun:.2e}")

# Show the three eigenvalues at this delta
vals_at_delta = []
for k in range(3):
    v = (v0_lep * (1 + np.sqrt(2) * np.cos(delta_lep + 2*np.pi*k/3)))**2
    vals_at_delta.append(v)
vals_sorted = sorted(vals_at_delta)
print(f"  Eigenvalues (sorted): {vals_sorted[0]:.4f}, {vals_sorted[1]:.4f}, {vals_sorted[2]:.4f}")
print(f"  Masses (sorted):      {m_e:.4f}, {m_mu:.4f}, {m_tau:.4f}")

# Seed comparison
print()
print(f"  Seed delta = 3pi/4 = {3*np.pi/4:.4f} rad = 135°")
# delta mod 2pi could be near 3pi/4 or offset by 2pi/3
# Find the shift relative to nearest seed equivalent
seed_candidates = [3*np.pi/4, 3*np.pi/4 + 2*np.pi/3, 3*np.pi/4 + 4*np.pi/3,
                   3*np.pi/4 - 2*np.pi/3, 3*np.pi/4 - 4*np.pi/3]
shifts = [(delta_lep - s) % (2*np.pi) for s in seed_candidates]
# Make shifts in [-pi, pi]
shifts = [(s if s < np.pi else s - 2*np.pi) for s in shifts]
best_shift = min(shifts, key=abs)
best_seed = seed_candidates[shifts.index(best_shift)]
print(f"  Nearest seed equivalent: {best_seed:.4f} rad = {best_seed*180/np.pi:.1f}°")
print(f"  Bloom shift dd = {best_shift:.4f} rad = {best_shift*180/np.pi:.2f}°")
print(f"  |dd|/(2pi/3) = {abs(best_shift)/(2*np.pi/3):.4f}")
if abs(best_shift) < 0.5:
    print(f"  → Small perturbation on the seed ✓")
else:
    print(f"  → NOT a small perturbation! Bloom is large for leptons.")

# =====================================================================
# 12. Cross-generation supermultiplet pairings
# =====================================================================
print("\n" + "=" * 70)
print("CROSS-GENERATION SUPERMULTIPLET PAIRINGS")
print("=" * 70)
print()

# The sBootstrap supermultiplets pair quarks ACROSS standard generations:
# (t,s), (u,b), (c,d)
# Each pair (q_up, q_down) sits in a chiral multiplet with a lepton partner.
# Under the Seiberg seesaw, M^i_j = Q^i Q_j / Lambda, so:
# M^c_s = Q^c Q_s / Lambda → the cs meson is a scalar partner of...?
#
# The cross-generation pairings arise because the O'R eigenvalues
# are m_+ = (2+sqrt3)*m, m_- = (2-sqrt3)*m, and the Seiberg seesaw
# inverts masses: M_j ~ 1/m_j.
#
# In the quark sector:
# O'R eigenvalue m_+ → quark mass proportional to (2+sqrt3)^2 → charm (or top)
# O'R eigenvalue m_- → quark mass proportional to (2-sqrt3)^2 → strange (or down)
# Zero eigenvalue → (nearly) massless → up (or electron)
#
# But the seesaw M_j ~ Lambda^2/m_j means:
# Light quark → heavy meson, heavy quark → light meson.
# So strange (light) → heavy M_s, charm (heavy) → light M_c.
# The meson masses are INVERTED relative to quark masses.

print("Seiberg seesaw inversion:")
print("  Quark masses: m_u < m_d < m_s < m_c < m_b < m_t")
print("  Meson masses: M ~ Lambda^2/m → inverted")
print()

# The key connection: leptons = Miyazawa partners of MESONS (not quarks).
# If the supermultiplet contains (quark_up, quark_down, lepton, meson):
# then the lepton mass should be related to the MESON mass, not quark mass.
#
# At t = sqrt(3), the three O'R eigenstates pair as:
# Fermion psi_X (mass 0) ↔ Scalar X (mass 0) → electron ↔ pion
# Fermion psi_- (mass m_-) ↔ Scalar_- (mass m_-) → muon ↔ kaon?
# Fermion psi_+ (mass m_+) ↔ Scalar_+ (mass m_+) → tau ↔ D meson?
#
# After SUSY breaking by B-term:
# The scalar gets shifted: m_scalar^2 = m_fermion^2 + delta_B
# where delta_B depends on Phi_1 content

print("SUSY-limit supermultiplet identification:")
print()
print("  Eigenstate  | Fermion (lepton) | Scalar (meson) | B-term shift")
print("  " + "-" * 65)
print(f"  psi_X (m=0) |     electron     |     pion       |  minimal (~0)")
print(f"  psi_- (m_-) |      muon        |     kaon?      |  2gf×sin²(π/12) = f_π²")
print(f"  psi_+ (m_+) |      tau         |     D/D_s?     |  2gf×cos²(π/12) = 14×f_π²")
print()

# Check the kaon assignment for the muon partner:
m_K = 493.677  # MeV (K±)
m_pi = 139.570  # MeV (pi±)

print("Testing kaon as muon partner:")
print(f"  m_K^2 - m_mu^2 = {m_K**2 - m_mu**2:.0f} MeV^2")
print(f"  m_pi^2 - m_mu^2 = {m_pi**2 - m_mu**2:.0f} MeV^2 = {(m_pi**2 - m_mu**2)/f_pi**2:.4f} × f_π²")
print(f"  f_π^2 = {f_pi**2:.0f} MeV^2")
print()
print("  The founding relation says pi-mu: m_pi^2 - m_mu^2 ≈ f_pi^2")
print("  If kaon were the muon partner instead of pion:")
print(f"    m_K^2 - m_mu^2 = {m_K**2 - m_mu**2:.0f} MeV^2 = {(m_K**2 - m_mu**2)/f_pi**2:.1f} × f_π²")
print("  This is 27× too large — kaon is NOT the muon's supermultiplet partner.")
print()

# The pion assignment works:
print("Pion as muon partner (founding relation):")
print(f"  m_pi^2 - m_mu^2 = {m_pi**2 - m_mu**2:.0f} MeV^2")
print(f"  f_pi^2 = {f_pi**2:.0f} MeV^2")
print(f"  Ratio = {(m_pi**2 - m_mu**2)/f_pi**2:.4f} (should be ~1)")
print()

# And the predicted tau partner:
print("From non-universal B-term splitting:")
print(f"  Muon pair splitting = {f_pi**2:.0f} MeV^2 (= f_π²)")
print(f"  Tau pair splitting = {f_pi**2 * phi1_content_plus/phi1_content_minus:.0f} MeV^2 = {phi1_content_plus/phi1_content_minus:.2f} × f_π²")
tau_partner_mass = np.sqrt(m_tau**2 + f_pi**2 * phi1_content_plus / phi1_content_minus)
print(f"  Tau partner mass = sqrt(m_tau^2 + {f_pi**2 * phi1_content_plus/phi1_content_minus:.0f}) = {tau_partner_mass:.1f} MeV")
print()

# What about the electron partner?
# Electron is from the zero eigenvalue (psi_X = goldstino).
# Its scalar partner is X (pseudo-modulus), also massless in SUSY limit.
# After SUSY breaking, X gets CW mass ~ gf/(4pi), but psi_X stays massless
# (absorbed as goldstino, or gets tiny bloom mass).
# The electron partner is NOT a standard meson — it's the pseudo-modulus X.
# In the sBootstrap, X maps to... what?
# If X is the pion: the pion is the pseudo-Goldstone boson, which is
# the lightest pseudoscalar. But we already used pion for muon!

print("The electron partner problem:")
print("  psi_X is the goldstino (mass 0 before bloom)")
print("  Its scalar partner X is the pseudo-modulus (mass 0 + CW)")
print("  After bloom: m_e = 0.511 MeV, so X gets mass ~0.511 MeV too?")
print("  But X is the pseudo-modulus = det(M)/Lambda^3 in ISS")
print("  This is a SINGLET, not a meson in the adjoint.")
print()
print("  Possible resolution: the X field maps to a pion-like state")
print("  (Goldstone of approximate chiral symmetry).")
print("  If m_X ~ m_e ~ 0.511 MeV, this is much lighter than m_pi.")
print("  The physical pion may be a MIXTURE of X and meson states.")

# =====================================================================
# 13. The (2+sqrt3)^2 structural constant
# =====================================================================
print("\n" + "=" * 70)
print("THE (2+√3)² STRUCTURAL CONSTANT")
print("=" * 70)
print()

r = (2 + np.sqrt(3))**2
print(f"(2+√3)² = {r:.6f}")
print()
print("This ratio appears in THREE independent ways:")
print(f"  1. Fermion mass ratio: m_+/m_- = {m_plus/m_minus:.4f}")
print(f"  2. Quark mass ratio: m_c/m_s = {r:.2f} (predicted: 1369/98.2 = {1369/98.2:.1f})")
print(f"  3. B-term splitting ratio: cos²(π/12)/sin²(π/12) = {np.cos(np.pi/12)**2/np.sin(np.pi/12)**2:.4f}")
print()
print("All three are consequences of t = gv/m = √3:")
print(f"  - The fermion mass matrix [[2t, 1],[1, 0]] has eigenvalue ratio (t+√(t²+1))/(−t+√(t²+1))")
print(f"  - At t=√3: ratio = (√3+2)/(−√3+2) = (2+√3)/(2−√3) = (2+√3)² = {r:.4f}")
print(f"  - The eigenvector mixing angle: tan(2θ) = 1/t = 1/√3 → θ = π/12")
print(f"  - cos²(π/12)/sin²(π/12) = (same ratio) = {r:.4f}")
print()
print("The structural constant (2+√3)² ≈ 13.93 is thus a CONSEQUENCE")
print("of the Koide condition gv/m = √3. It controls both the mass hierarchy")
print("and the SUSY-breaking splitting pattern simultaneously.")

# =====================================================================
# 11. Consistency: does the non-universal splitting reproduce
#     the meson spectrum at higher order?
# =====================================================================
print("\n\n" + "=" * 70)
print("EXTENDED: MESON-LEPTON PAIRING AT PHYSICAL MASSES")
print("=" * 70)

print("\n--- 11. Checking known meson-lepton mass relations ---\n")

# Known pseudoscalar mesons with Q = -1 or Q = 0 (neutral)
mesons = {
    'pi±': (139.57, 'ud-bar'),
    'K±':  (493.68, 'us-bar'),
    'D±':  (1869.66, 'cd-bar'),
    'D_s±': (1968.35, 'cs-bar'),
    'B±':  (5279.34, 'ub-bar'),
    'B_c±': (6274.47, 'cb-bar'),
}

leptons_phys = {'e': m_e, 'mu': m_mu, 'tau': m_tau}

print("m^2_meson - m^2_lepton for all candidate pairs:")
print(f"{'Meson':>8s} {'m_M':>8s} | {'Lepton':>6s} {'m_L':>8s} | {'m_M²-m_L²':>12s} {'sqrt':>8s} {'ratio/f_pi²':>12s}")
print("-" * 80)

for mname, (m_M, content) in sorted(mesons.items(), key=lambda x: x[1][0]):
    for lname, m_L in sorted(leptons_phys.items(), key=lambda x: x[1]):
        diff = m_M**2 - m_L**2
        if diff > 0:
            sqrt_diff = np.sqrt(diff)
            ratio = diff / f_pi**2
            marker = " ← founding" if mname == 'pi±' and lname == 'mu' else ""
            marker = marker or (" ← tau partner?" if abs(sqrt_diff - 581) < 50 and lname == 'tau' else "")
            print(f"{mname:>8s} {m_M:8.1f} | {lname:>6s} {m_L:8.3f} | {diff:12.0f} {sqrt_diff:8.1f} {ratio:12.4f}{marker}")

print()
print("--- Key mass-squared differences ---")
print(f"  pi - mu:  {m_pi**2 - m_mu**2:.0f} MeV² = {(m_pi**2 - m_mu**2)/f_pi**2:.4f} × f_pi²  ← founding relation")
print(f"  K - mu:   {493.68**2 - m_mu**2:.0f} MeV² = {(493.68**2 - m_mu**2)/f_pi**2:.2f} × f_pi²")
print(f"  K - tau:  {493.68**2 - m_tau**2:.0f} MeV² (NEGATIVE — K lighter than tau)")
print(f"  D - tau:  {1869.66**2 - m_tau**2:.0f} MeV² = {(1869.66**2 - m_tau**2)/f_pi**2:.2f} × f_pi²")
print(f"  D_s - tau: {1968.35**2 - m_tau**2:.0f} MeV² = {(1968.35**2 - m_tau**2)/f_pi**2:.2f} × f_pi²")
print()

# The non-universal prediction: splitting ratio = |a_+|²/|a_-|² = 0.933/0.067 = 13.93
# That is: (m²_tau_partner - m²_tau) / (m²_pi - m²_mu) = 13.93
ratio_predicted = phi1_content_plus / phi1_content_minus
actual_pi_mu = m_pi**2 - m_mu**2
predicted_tau_split = actual_pi_mu * ratio_predicted
predicted_tau_partner = np.sqrt(m_tau**2 + predicted_tau_split)

print("--- Non-universal B-term prediction ---")
print(f"  Splitting ratio |a_+|²/|a_-|² = {ratio_predicted:.4f}")
print(f"  = (2+√3)²/(2-√3)² ... wait, is this (m+/m-)²?")
print(f"    (m+/m-)² = {(m_plus/m_minus)**2:.4f}")
print(f"    Ratio = {ratio_predicted:.4f}")
print(f"    These are NOT the same.")
print()
print(f"  From pi-mu: Delta² = {actual_pi_mu:.0f} MeV²")
print(f"  Predicted tau splitting: {ratio_predicted:.2f} × {actual_pi_mu:.0f} = {predicted_tau_split:.0f} MeV²")
print(f"  Predicted tau partner mass: sqrt({m_tau**2:.0f} + {predicted_tau_split:.0f}) = {predicted_tau_partner:.1f} MeV")
print()

# Check: what is the ratio 0.933/0.067 exactly?
# At t = sqrt(3), the eigenvectors are:
# psi_+ = cos(alpha) Phi1 + sin(alpha) Phi3
# psi_- = -sin(alpha) Phi1 + cos(alpha) Phi3
# where tan(2alpha) = 2m/(2gv) = 1/t = 1/sqrt(3) → 2alpha = pi/6 → alpha = pi/12
alpha = np.arctan2(m_param, 2*g*v) / 1  # actually: tan(alpha) = m/(2gv - lam) for each eigenvector

# More carefully: the mixing angle
# M_F = [[2gv, m], [m, 0]]
# For a 2x2 matrix [[a,b],[b,0]], the eigenvectors have angle
# tan(2theta) = 2b/(a-0) = 2m/(2gv) = 1/t
theta_mix = 0.5 * np.arctan2(2*m_param, 2*g*v)
print(f"  Mixing angle theta = {theta_mix:.6f} rad = {np.degrees(theta_mix):.2f}°")
print(f"  cos²(theta) = {np.cos(theta_mix)**2:.6f} (= Phi1 content of psi_+)")
print(f"  sin²(theta) = {np.sin(theta_mix)**2:.6f} (= Phi1 content of psi_-)")
print(f"  tan(2*theta) = {np.tan(2*theta_mix):.6f} = 1/sqrt(3) = {1/np.sqrt(3):.6f}")
print(f"  So 2*theta = pi/6 → theta = pi/12 = 15°")
print()
print(f"  The Phi1 content ratio is:")
print(f"    cos²(pi/12) / sin²(pi/12) = {np.cos(np.pi/12)**2/np.sin(np.pi/12)**2:.4f}")
print(f"    = (2+sqrt3)/(2-sqrt3) = {(2+np.sqrt(3))/(2-np.sqrt(3)):.4f}")
print(f"    = (2+sqrt3)² = {(2+np.sqrt(3))**2:.4f}")
print()
print("  REMARKABLE: the B-term splitting ratio equals the O'R mass ratio!")
print("  The same (2+sqrt3)² ≈ 13.93 that gives m_c/m_s also gives")
print("  the ratio of SUSY-breaking splittings between tau and muon pairs.")
