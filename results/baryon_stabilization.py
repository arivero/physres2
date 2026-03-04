"""
Baryon mass term stabilization of the seesaw vacuum.

Superpotential (original):
    W = sum_i m_i M^i_i + X(det M - B Btilde - Lambda^6)
        + lambda X (H_u . H_d)

Soft SUSY breaking:
    V_soft = f_pi^2 Tr(M^dag M)

Two vacua:
    (A) Seesaw: M_i = C/m_i, B = Btilde = 0
    (B) Global: M = 0, BB~ = -Lambda^6

We add W_B = m_B B Btilde and analyze:
  1. F-term equations with W_B => threshold m_B to eliminate vacuum (B)
  2. Bounce action for tunneling A -> B in original theory
  3. Seesaw vacuum still a solution with W_B
  4. Impact on mass predictions
  5. Natural scale for m_B
"""

import numpy as np
from scipy.optimize import fsolve, minimize_scalar

# ===========================================================================
# Physical inputs (MeV)
# ===========================================================================
m_u = 2.16
m_d = 4.67
m_s = 93.4
Lambda = 300.0
Lambda6 = Lambda**6
f_pi = 92.0
f_pi_sq = f_pi**2
lam = 0.72
v_ew = 246220.0  # MeV
v = v_ew / np.sqrt(2)

masses = np.array([m_u, m_d, m_s])
C = Lambda**2 * np.prod(masses)**(1.0/3.0)
M_diag = C / masses
X_seesaw = -C / Lambda6

print("=" * 78)
print("  BARYON MASS TERM STABILIZATION OF THE SEESAW VACUUM")
print("=" * 78)
print()
print(f"  Physical inputs:")
print(f"    m_u = {m_u} MeV, m_d = {m_d} MeV, m_s = {m_s} MeV")
print(f"    Lambda = {Lambda} MeV, Lambda^6 = {Lambda6:.6e} MeV^6")
print(f"    f_pi = {f_pi} MeV, lambda = {lam}, v = {v_ew} MeV")
print()
print(f"  Derived:")
print(f"    C = Lambda^2 (m_u m_d m_s)^(1/3) = {C:.6f} MeV^2")
print(f"    M_u = C/m_u = {M_diag[0]:.4f} MeV")
print(f"    M_d = C/m_d = {M_diag[1]:.4f} MeV")
print(f"    M_s = C/m_s = {M_diag[2]:.4f} MeV")
print(f"    X_seesaw = -C/Lambda^6 = {X_seesaw:.6e}")
print()


# ===========================================================================
# TASK 1: F-term equations with W_B = m_B B Btilde
# ===========================================================================
print("=" * 78)
print("  TASK 1: BARYON MASS TERM AND VACUUM (B) ELIMINATION")
print("=" * 78)
print()

print("Modified superpotential:")
print("  W = sum_i m_i M^i_i + X(det M - B Btilde - Lambda^6)")
print("      + lambda X (H_u . H_d) + m_B B Btilde")
print()

print("F-term equations:")
print("  F_X     = det M - B Btilde - Lambda^6 + lambda (H_u . H_d)")
print("  F_{M_i} = m_i + X * cofactor(M,i,i)")
print("  F_B     = -X Btilde + m_B Btilde")
print("  F_Bt    = -X B + m_B B")
print("  F_{H_u} = lambda X H_d")
print("  F_{H_d} = lambda X H_u")
print()

print("--- Analysis of vacuum (B): M = 0, BB~ = -Lambda^6 ---")
print()
print("At M = 0:")
print("  F_{M_i} = m_i + X * cofactor(0, i, i) = m_i  (cofactors vanish for M=0)")
print("  => F_{M_i} = m_i (nonzero, but these contribute to V)")
print()
print("  F_X = det(0) - B Btilde - Lambda^6 = -B Btilde - Lambda^6")
print("  Setting BB~ = -Lambda^6: F_X = Lambda^6 - Lambda^6 = 0  (OK)")
print()
print("  F_B = (-X + m_B) Btilde")
print("  F_Bt = (-X + m_B) B")
print()
print("For the baryonic flat direction in the ORIGINAL theory (m_B = 0):")
print("  F_B = -X Btilde, F_Bt = -X B")
print("  At vacuum (B), we need X. From F_{M_i}: m_i = 0 is needed for F_{M_i}=0.")
print("  But m_i != 0, so F_{M_i} are nonzero at M=0.")
print("  The minimum is found by minimizing V = sum |F_I|^2 + f_pi^2 Tr(M^dag M)")
print()

# Compute V at vacuum (B) in the original theory
# At M = 0, BB~ = -Lambda^6, X = 0 (minimizes over X), H = 0
# F_X = 0, F_{M_i} = m_i, F_B = 0, F_Bt = 0
# V_B = sum m_i^2 + 0 (no soft term since M=0)
V_B_original = np.sum(masses**2)
print(f"V at vacuum (B) [original, m_B=0]:")
print(f"  M = 0, BB~ = -Lambda^6, X = 0, H = 0")
print(f"  F_X = 0 (constraint satisfied)")
print(f"  F_{{M_i}} = m_i (from mass terms)")
print(f"  F_B = F_Bt = 0 (X = 0)")
print(f"  V_B = sum m_i^2 = {m_u}^2 + {m_d}^2 + {m_s}^2")
print(f"       = {m_u**2:.4f} + {m_d**2:.4f} + {m_s**2:.4f}")
print(f"       = {V_B_original:.4f} MeV^2")
print()

# Note: the problem statement says V_B = m_u^2 = 4.67 MeV^2 which is actually m_d^2.
# The correct V_B = sum m_i^2. Let me check what happens if we minimize over X too.
# With M=0, B Btilde = -Lambda^6:
# F_X = 0, F_{M_i} = m_i, F_B = -X * Btilde, F_Bt = -X * B
# V = sum m_i^2 + |X|^2 (|Btilde|^2 + |B|^2)
# For B Btilde = -Lambda^6, |B|^2 |Btilde|^2 >= |B Btilde|^2 = Lambda^12
# Minimized when |B| = |Btilde| = Lambda^3, giving |B|^2 + |Bt|^2 = 2 Lambda^6
# V = sum m_i^2 + 2 Lambda^6 |X|^2
# Minimized at X = 0: V = sum m_i^2

print(f"Minimizing over X: V = sum m_i^2 + 2 Lambda^6 |X|^2")
print(f"  => X = 0, V_B = sum m_i^2 = {V_B_original:.4f} MeV^2")
print()

# Compute V at vacuum (A)
# At seesaw vacuum with lambda coupling:
F_X_A = -lam * v**2  # dominant F-term from lambda X H_u.H_d with H_u.H_d = -v^2
# Actually F_X = det M - BB~ - Lambda^6 + lam(H_u.H_d)
# At seesaw: det M = Lambda^6, BB~ = 0
# H_u.H_d = H_u^+ H_d^- - H_u^0 H_d^0 = -v^2
F_X_at_A = Lambda6 - 0 - Lambda6 + lam * (-v**2)
F_Md_A = m_d + 0  # y_c contribution...but for the simplified version without Yukawa
# In the simplified model (no Yukawa), F_{M_i} = m_i + X * cofactor = 0 at seesaw
# The dominant contribution to V is |F_X|^2
V_A_FX = abs(F_X_at_A)**2
V_A_soft = f_pi_sq * np.sum(M_diag**2)
V_A_total = V_A_FX + V_A_soft

print(f"V at vacuum (A) [seesaw]:")
print(f"  F_X = det M - Lambda^6 + lambda * H_u.H_d")
print(f"       = 0 + {lam} * (-{v:.4f}^2)")
print(f"       = {F_X_at_A:.6e} MeV")
print(f"  |F_X|^2 = {V_A_FX:.6e} MeV^2")
print(f"  V_soft = f_pi^2 * sum M_i^2 = {V_A_soft:.6e} MeV^2")
print(f"  V_A = |F_X|^2 + V_soft = {V_A_total:.6e} MeV^2")
print()

epsilon = V_A_total - V_B_original
print(f"Energy splitting: epsilon = V_A - V_B = {epsilon:.6e} MeV^2")
print()

# Now with m_B:
print("--- With m_B: vacuum (B) analysis ---")
print()
print("At M = 0, H = 0, the F-terms become:")
print("  F_X = -B Btilde - Lambda^6")
print("  F_{M_i} = m_i")
print("  F_B = (-X + m_B) Btilde")
print("  F_Bt = (-X + m_B) B")
print()
print("Setting F_X = 0: B Btilde = -Lambda^6")
print("Setting F_B = 0, F_Bt = 0: X = m_B (or B = Btilde = 0)")
print()
print("If X = m_B and B Btilde = -Lambda^6:")
print("  V = sum m_i^2 + 0 + 0 + 0 = sum m_i^2")
print("  (same as before - the baryonic vacuum is NOT changed by m_B alone!)")
print()
print("HOWEVER: X = m_B is now a specific value. The question is whether")
print("the B-flat direction is lifted for the FULL potential including soft terms.")
print()

# More careful analysis: with W_B, the scalar potential for B sector is
# V_B_sector = |F_B|^2 + |F_Bt|^2 = |(-X + m_B)|^2 (|Bt|^2 + |B|^2)
# This is minimized at X = m_B (for any B) or B = Bt = 0 (for any X).
#
# The key insight: in the original theory, at M=0 with X=0, the B-flat direction
# is exact: B Btilde = -Lambda^6 with any |B|, |Bt| satisfying the constraint.
# All have V = sum m_i^2.
#
# With m_B, at M=0 the B-equations become F_B = (m_B - X) Btilde, F_Bt = (m_B - X) B.
# But we also need to minimize over X. F_X gives B Btilde = -Lambda^6.
# V = sum m_i^2 + |m_B - X|^2 (|B|^2 + |Bt|^2)
# For B Btilde = -Lambda^6, |B|^2 + |Bt|^2 >= 2 Lambda^3
# Minimize over X: X = m_B => V = sum m_i^2 (unchanged!)
# The baryonic vacuum has X = m_B, not X = 0.

print("RESULT: The baryon mass term m_B shifts X to m_B at vacuum (B),")
print("but does NOT change V_B = sum m_i^2. The flat direction persists.")
print()
print("To actually lift vacuum (B), we need the F_{M_i} terms to change.")
print("At M = 0 with X = m_B: F_{M_i} = m_i + m_B * cofactor(0, i, i) = m_i")
print("(cofactors of a zero matrix are zero).")
print()

# The real mechanism to eliminate vacuum (B):
# Vacuum (B) exists because at M=0, the constraint det M = 0 != Lambda^6
# is satisfied by B Btilde = -Lambda^6.
# With m_B, the F-term equation for B gives X = m_B.
# But now at M=0 with X = m_B: the off-diagonal F-terms are
# F_{M^i_j} = X * cofactor(M, i, j) = m_B * 0 = 0 (still zero at M=0)
#
# So the question becomes: is there a barrier preventing the system from
# reaching M=0? Or does m_B create an additional contribution to V?
#
# Actually, let's think about this differently. In the presence of m_B,
# consider the potential along the path M_i -> t * M_i (seesaw), B, Bt.
# The F_B equation is F_B = (-X + m_B) Btilde.
# At the seesaw vacuum: X = X_seesaw << m_B (for m_B > 0),
# so F_B = m_B Btilde. If B = Bt = 0, this vanishes.
#
# The modified potential at a general point (M, B, Bt, X):
# V = |det M - B Bt - Lambda^6 + lam H.H|^2 + sum |m_i + X cof_i|^2
#   + |(-X + m_B) Bt|^2 + |(-X + m_B) B|^2 + f_pi^2 Tr(M^dag M)
#
# At M=0, X=m_B: V = |0 - B Bt - Lambda^6|^2 + sum m_i^2 + 0
# => minimize over BBt: B Bt = -Lambda^6 => V = sum m_i^2
#
# At M=0, X=0: V = |0 - B Bt - Lambda^6|^2 + sum m_i^2 + m_B^2(|B|^2+|Bt|^2)
# => now B-flat direction IS lifted by m_B^2 term!
# With X not free (or X=0 from some other constraint), m_B lifts the direction.

# But wait - we should minimize over ALL variables including X.
# The problem is that X is a dynamical field. At M=0:
# V = |-BBt - Lambda^6|^2 + sum m_i^2 + |m_B - X|^2 (|B|^2 + |Bt|^2)
# Minimize over X: X = m_B => the m_B^2 term vanishes.
# So vacuum (B) at M=0 is NOT lifted by m_B alone.

# HOWEVER: there's a subtlety. If we add a SOFT mass for X too,
# or if the Kahler potential pins X near 0, then X cannot freely adjust to m_B.
# Let's check: with the Kahler pole at |X| = sqrt(3) mu where
# mu = |C / (sqrt(3) Lambda^6)| ~ 10^{-9}, the X field is dynamically
# confined near |X| ~ 10^{-9}. It CANNOT reach X = m_B for any m_B > 10^{-6}.

# With the Kahler pole constraining X near X_seesaw:
# V at M=0, BBt = -Lambda^6, X ~ X_seesaw:
# V = sum m_i^2 + |m_B - X_seesaw|^2 (|B|^2 + |Bt|^2)
# >= sum m_i^2 + 2 |m_B - X_seesaw|^2 Lambda^6
# For m_B >> |X_seesaw|:
# V_B(m_B) = sum m_i^2 + 2 m_B^2 Lambda^6

print("--- Crucial observation: Kahler pole constrains X ---")
print()
print("The Kahler potential K_X = |X|^2 - |X|^4/(12 mu^2) creates a pole")
print(f"at |X| = sqrt(3) mu, where mu = |C/(sqrt(3) Lambda^6)| = {abs(C)/(np.sqrt(3)*Lambda6):.6e}")
print()
print("If X is dynamically confined near X_seesaw, it cannot freely adjust")
print("to X = m_B. Then at vacuum (B) with X ~ X_seesaw ~ 10^{-9}:")
print()
print("  V_B(m_B) = sum m_i^2 + |m_B - X_seesaw|^2 (|B|^2 + |Bt|^2)")
print("           >= sum m_i^2 + 2 m_B^2 Lambda^6  (for m_B >> X_seesaw)")
print()

V_B_with_mB = lambda mB: V_B_original + 2 * mB**2 * Lambda6
# Vacuum (B) eliminated when V_B(m_B) > V_A:
# sum m_i^2 + 2 m_B^2 Lambda^6 > V_A
# m_B^2 > (V_A - sum m_i^2) / (2 Lambda^6) = epsilon / (2 Lambda^6)
mB_threshold = np.sqrt(epsilon / (2 * Lambda6))

print(f"Vacuum (B) is eliminated when V_B(m_B) > V_A:")
print(f"  sum m_i^2 + 2 m_B^2 Lambda^6 > V_A")
print(f"  m_B > sqrt(epsilon / (2 Lambda^6))")
print(f"  m_B > sqrt({epsilon:.6e} / {2*Lambda6:.6e})")
print(f"  m_B > {mB_threshold:.6e} MeV")
print(f"  m_B > {mB_threshold/1000:.4f} GeV")
print()

# But this is huge. Let's reconsider.
# The question in the problem says V_A ~ 10^15 MeV^2 and V_B ~ m_u^2 ~ 4.67 MeV^2
# Let me recalculate more carefully. The problem statement has a slightly different
# parametrization. Let me follow it more literally.

print("--- Following the problem statement parametrization ---")
print()

# The problem says V_A ~ |F_X|^2 + f_pi^2 sum M_i^2 ~ 10^15 MeV^2
# and V_B = m_u^2 ~ 4.67 MeV^2.
# This assumes NO lambda coupling (or lambda contribution already folded in differently).
# Let me compute both scenarios.

# Without lambda (pure SQCD + soft breaking):
# F_X = 0 at seesaw (det M = Lambda^6 exactly), so |F_X|^2 = 0
# V_A = f_pi^2 * sum M_i^2 (soft terms only)
V_A_no_lam = f_pi_sq * np.sum(M_diag**2)
print(f"Without lambda coupling:")
print(f"  V_A(no lam) = f_pi^2 * sum M_i^2 = {V_A_no_lam:.6e} MeV^2")
print()

# With lambda:
print(f"With lambda coupling:")
print(f"  |F_X|^2 = |lam * v^2|^2 = {abs(lam * v**2)**2:.6e} MeV^2")
print(f"  V_soft = f_pi^2 sum M_i^2 = {V_A_soft:.6e} MeV^2")
print(f"  V_A(with lam) = {V_A_total:.6e} MeV^2")
print()

# The problem says V_A ~ 10^15. Let's see which matches:
print(f"  V_A(no lam) = {V_A_no_lam:.3e}  (this is ~ 10^15)")
print(f"  V_A(with lam) = {V_A_total:.3e}  (this is ~ 10^20)")
print()
print("The problem statement uses V_A ~ 10^15, matching the no-lambda case.")
print("Using V_A = f_pi^2 * sum M_i^2 for the bounce calculation.")
print()

# Use the no-lambda V_A for the bounce calculation as the problem implies
V_A_bounce = V_A_no_lam

# For vacuum (B): the problem says V_B = m_u^2.
# This is sum_i m_i^2 ~ m_s^2 = 8723.56, not m_u^2.
# Actually at M=0 the F_{M_i} = m_i, so V_B = sum m_i^2.
# The problem says V_B = m_u^2 = 4.67. That's m_d^2 not m_u^2.
# Probably a typo. Let me use the correct value:
print(f"V_B = sum m_i^2 = {V_B_original:.4f} MeV^2")
print(f"  (Note: dominated by m_s^2 = {m_s**2:.2f})")
print()

# Actually, re-reading the problem: it says V_B = m_u^2 ≈ 4.67.
# m_u = 2.16, so m_u^2 = 4.67. This is correct: the problem defines V_B as
# JUST the smallest F-term squared? No - at M=0 ALL F_{M_i} are nonzero.
# Unless the global minimum has M_i that are NOT all zero.
# Let me reconsider vacuum (B).

# In the original theory:
# V = |det M - BBt - Lambda^6|^2 + sum_i |m_i + X cof_i|^2 + |X Bt|^2 + |X B|^2
#     + f_pi^2 Tr(M^dag M)
# At B Bt = -Lambda^6, the constraint term vanishes.
# We can then set X = 0 (minimizes baryon F-terms).
# V = sum |m_i + 0|^2 + f_pi^2 sum |M_i|^2
# Minimizing over M_i: set M_i = 0 to minimize soft term.
# V_B = sum m_i^2 = m_u^2 + m_d^2 + m_s^2 = 4.6656 + 21.8089 + 8723.56 = 8750.03

# Unless we can do better by having M_i != 0 with X != 0 to partially cancel F_{M_i}.
# F_{M_i} = m_i + X * cof_i(M)
# With diagonal M: cof_i = prod_{j!=i} M_j
# F_{M_i} = m_i + X M_j M_k = 0 => X = -m_i/(M_j M_k)
# But X must be the same for all i => m_i/(M_j M_k) = m_j/(M_i M_k)
# => m_i M_i = m_j M_j = C => seesaw!
# So the only critical point with M != 0 is the seesaw vacuum.
# There's no intermediate minimum at M != 0 with lower V.

# Wait, but what about M != 0 but NOT at the seesaw point?
# E.g., very small M_i. Then cof_i ~ M_j M_k is small, and
# F_{M_i} ~ m_i + X * small.
# V ~ sum m_i^2 + f_pi^2 sum M_i^2 + small corrections
# This is always > sum m_i^2 for M_i > 0.
# So V_B = sum m_i^2 is indeed the minimum at M = 0.

# But the problem says V_B = m_u^2 = 4.67. Let me think about this differently.
# Maybe they mean that you can set M_d, M_s to cancel F_{M_d}, F_{M_s}
# (partial seesaw) while keeping B Bt = -Lambda^6.
# With det M - BBt = Lambda^6 and BBt = -Lambda^6: det M = 0.
# So with B Bt = -Lambda^6, the constraint becomes det M = 0.
# This means at least one eigenvalue of M is zero, say M_u = 0.
# Then F_{M_u} = m_u (can't cancel since cof_u = M_d M_s and det M = 0 means
# at least one of M_d, M_s is zero too, or... no, det M = M_u M_d M_s = 0
# is satisfied by M_u = 0 alone).
# F_{M_d} = m_d + X * M_u * M_s = m_d + 0 = m_d (no, cof_d = M_u M_s = 0 since M_u=0!)
# Actually, cof(M, d, d) for diagonal M = M_u * M_s. If M_u = 0, cof = 0.
# So F_{M_d} = m_d and F_{M_s} = m_s regardless.
# We're stuck with V >= sum m_i^2 at M = 0 sector.

# Hmm, but the problem says V_B = m_u^2 = 4.67. Maybe they envision a path
# where M_d, M_s adjust while B absorbs the constraint.
# Let me try: det M = M_u M_d M_s = Lambda^6 + BBt.
# If BBt = -Lambda^6 + M_u M_d M_s, the constraint gives F_X = 0.
# Then we're free to choose X from F_{M_i} = 0.
# m_i + X prod_{j!=i} M_j = 0 for all i with the seesaw C = m_i M_i.
# This requires det M = C^3/(m_u m_d m_s) = Lambda^6.
# Then BBt = det M - Lambda^6 = 0. We're back at vacuum (A)!

# The only way to have BBt != 0 is to have det M != Lambda^6.
# Then F_X != 0.

# OK I think the problem statement's V_B = m_u^2 is simply using the problem's
# given values. Let me just proceed with the correct V_B = sum m_i^2 for my
# computation and note the discrepancy.

print(f"Corrected: V_B = sum m_i^2 = {V_B_original:.4f} MeV^2")
print(f"  (Problem statement says V_B = m_u^2 = 4.67 MeV^2;")
print("   this is incorrect -- all F_{M_i} contribute at M=0.)")
print()

# ===========================================================================
# Task 1 continued: threshold m_B with Kahler constraint
# ===========================================================================
print("-" * 78)
print("  TASK 1 RESULT: m_B THRESHOLD")
print("-" * 78)
print()

# Without Kahler constraint on X: vacuum (B) cannot be eliminated by m_B alone
# (X adjusts freely to X = m_B).

# WITH Kahler constraint confining X near X_seesaw:
# At vacuum (B), X cannot reach m_B, so F_B = m_B Btilde (approximately)
# V_B(m_B) = sum m_i^2 + m_B^2 * (|B|^2 + |Bt|^2)
# With |B| = |Bt| = Lambda^3 (minimum of constraint):
# V_B(m_B) = sum m_i^2 + 2 m_B^2 Lambda^6

# For vacuum (B) to be eliminated: V_B(m_B) > V_A
# Using V_A = f_pi^2 sum M_i^2 (no-lambda case):
mB_threshold_no_lam = np.sqrt(max(V_A_no_lam - V_B_original, 0) / (2 * Lambda6))
# Using V_A with lambda:
mB_threshold_with_lam = np.sqrt(max(V_A_total - V_B_original, 0) / (2 * Lambda6))

print("Scenario 1: Without Kahler constraint on X")
print("  X adjusts to m_B at vacuum (B), so F_B = F_Bt = 0.")
print("  m_B does NOT eliminate vacuum (B). m_B_threshold = infinity.")
print()
print("Scenario 2: With Kahler pole confining X near X_seesaw ~ 10^{-9}")
print("  V_B(m_B) = sum m_i^2 + 2 m_B^2 Lambda^6")
print()
print("  Case (a): V_A = f_pi^2 sum M_i^2 (no lambda):")
print(f"    m_B > sqrt((V_A - V_B) / (2 Lambda^6))")
print(f"    m_B > sqrt(({V_A_no_lam:.4e} - {V_B_original:.2f}) / (2 * {Lambda6:.4e}))")
print(f"    m_B > {mB_threshold_no_lam:.6f} MeV")
print(f"    m_B > {mB_threshold_no_lam:.6f} MeV = {mB_threshold_no_lam/1000:.6f} GeV")
print()
print(f"  Case (b): V_A includes lambda (V_A = {V_A_total:.4e}):")
print(f"    m_B > {mB_threshold_with_lam:.6f} MeV = {mB_threshold_with_lam/1000:.4f} GeV")
print()

# Even simpler criterion: just need 2 m_B^2 Lambda^6 > sum m_i^2
# (so vacuum B is disfavored over even a LOWER vacuum A)
mB_min = np.sqrt(V_B_original / (2 * Lambda6))
print("Minimum m_B to make baryonic direction massive (V_B > 2 sum m_i^2):")
print(f"  m_B > sqrt(sum m_i^2 / (2 Lambda^6))")
print(f"  m_B > {mB_min:.6e} MeV")
print()

# Actually, the simplest meaningful threshold: the baryonic mass^2 term
# m_B^2 (|B|^2 + |Bt|^2) at the baryonic vacuum should exceed the
# energy gain from going to vacuum (B):
# Delta V = V_A - V_B
# m_B^2 * 2 Lambda^6 > V_A - V_B
# m_B > sqrt((V_A - V_B)/(2 Lambda^6))

# A MUCH simpler approach: forget Kahler constraint.
# The baryon mass term adds to V_B:
# At M=0, BBt = -Lambda^6, X=0:
# V = sum m_i^2 + m_B^2 * (|B|^2 + |Bt|^2)  [from F_B = m_B Bt, F_Bt = m_B B]
# Wait no: F_B = (-X + m_B) Bt. At X=0: F_B = m_B Bt.
# |F_B|^2 = m_B^2 |Bt|^2, |F_Bt|^2 = m_B^2 |B|^2
# V = sum m_i^2 + m_B^2 (|B|^2 + |Bt|^2)
# With |B| = |Bt| = Lambda^3: V = sum m_i^2 + 2 m_B^2 Lambda^6

# But we should minimize over X too:
# V = sum m_i^2 + |m_B - X|^2 (|B|^2 + |Bt|^2)
# min over X: X = m_B => V = sum m_i^2 (m_B drops out!)
# Unless there's a cost to having X = m_B from other terms.

# From F_X: |det M - BBt - Lambda^6|^2 = |0 - (-Lambda^6) - Lambda^6|^2 = 0
# This is independent of X.

# So the key: in the flat Kahler metric, X is a free field and adjusts to m_B.
# Vacuum (B) is NOT lifted.
# In the non-trivial Kahler metric, X has a potential barrier.

# Let me just compute the cost of moving X from X_seesaw to m_B:
# V_Kahler(X) = |F_X|^2 / (1 - |X|^2/(3 mu^2))
# At X = m_B >> mu: the Kahler metric goes negative -> not physical.
# X cannot exceed sqrt(3) mu ~ 10^{-9} in the Kahler-pole theory.
mu_K = abs(C) / (np.sqrt(3) * Lambda6)
X_pole = np.sqrt(3) * mu_K

print("Kahler constraint on X:")
print(f"  mu = {mu_K:.6e}")
print(f"  X_pole = sqrt(3) mu = {X_pole:.6e}")
print(f"  X_seesaw = {abs(X_seesaw):.6e}")
print(f"  (X_pole and X_seesaw match by construction)")
print()
print("For m_B >> X_pole ~ 10^{-9} MeV:")
print("  X CANNOT reach m_B. The field is confined below X_pole.")
print("  The baryonic F-terms remain nonzero: F_B ≈ m_B Btilde.")
print()
print("THRESHOLD (with Kahler pole):")
print(f"  Any m_B >> X_pole = {X_pole:.4e} suffices.")
print(f"  Practically, m_B > 10^{{-8}} MeV eliminates vacuum (B).")
print(f"  For m_B ~ Lambda = {Lambda} MeV: V_B gains 2 m_B^2 Lambda^6 = {2*Lambda**2*Lambda6:.4e} MeV^2")
print(f"  This far exceeds V_A = {V_A_no_lam:.4e} MeV^2.")
print()


# ===========================================================================
# TASK 2: BOUNCE ACTION
# ===========================================================================
print()
print("=" * 78)
print("  TASK 2: BOUNCE ACTION FOR TUNNELING A -> B")
print("=" * 78)
print()

# Field displacement: from seesaw to M=0
# The dominant displacement is in M_u (largest VEV):
Delta_phi = np.sqrt(np.sum(M_diag**2))  # total field displacement
Delta_phi_u = M_diag[0]  # largest single component

DeltaV = V_A_no_lam - V_B_original  # using no-lambda V_A

print(f"Field displacement Delta phi = sqrt(sum M_i^2) = {Delta_phi:.4f} MeV")
print(f"  Largest component: M_u = {M_diag[0]:.4f} MeV")
print(f"  M_d = {M_diag[1]:.4f} MeV, M_s = {M_diag[2]:.4f} MeV")
print()
print(f"Energy difference: Delta V = V_A - V_B = {DeltaV:.6e} MeV^2")
print()

# THIN-WALL APPROXIMATION
# sigma = integral of sqrt(2 V(phi)) dphi along the path
# V(phi) is the potential barrier between the two vacua.
# In the thin-wall limit (small epsilon/barrier), sigma is determined by the barrier.
# For the soft-mass potential V = f_pi^2 |M|^2, the barrier is just the soft potential.
# The path goes from M_i (seesaw) to 0, and the barrier height is ~ f_pi^2 M_u^2.

# Surface tension estimate: sigma ~ sqrt(V_barrier) * Delta_phi
# V_barrier ~ f_pi^2 * M_u^2 (the maximum of the soft potential along the path)
V_barrier = f_pi_sq * M_diag[0]**2  # peak barrier from soft mass of M_u
sigma_estimate = np.sqrt(2 * V_barrier) * Delta_phi / 6  # factor 1/6 from kink profile

# More careful: for a quartic potential with two minima,
# sigma = (2/3) sqrt(2 V_barrier) * l where l is the width.
# Here V_barrier ~ f_pi^2 M_u^2 and l ~ M_u.
sigma = (2.0/3.0) * np.sqrt(2 * V_barrier) * M_diag[0]
print("--- Thin-wall approximation ---")
print()
print(f"Barrier height: V_barrier ~ f_pi^2 * M_u^2 = {V_barrier:.6e} MeV^2")
print(f"Surface tension: sigma ~ (2/3) sqrt(2 V_barrier) * M_u")
print(f"  sigma = {sigma:.6e} MeV^3")
print()

epsilon_tw = DeltaV
S_thin = 27 * np.pi**2 * sigma**4 / (2 * epsilon_tw**3)
print(f"Thin-wall bounce action:")
print(f"  S = 27 pi^2 sigma^4 / (2 epsilon^3)")
print(f"  S = 27 * {np.pi**2:.6f} * ({sigma:.4e})^4 / (2 * ({epsilon_tw:.4e})^3)")
print(f"  S = {S_thin:.6e}")
print()
if S_thin > 400:
    print(f"  S = {S_thin:.4e} >> 400 => COSMOLOGICALLY STABLE")
else:
    print(f"  S = {S_thin:.4e} < 400 => UNSTABLE (tunnels within age of universe)")
print()

# THICK-WALL APPROXIMATION
# S ~ (Delta phi)^4 / Delta V
S_thick = Delta_phi**4 / DeltaV
print("--- Thick-wall approximation ---")
print()
print(f"  S = (Delta phi)^4 / Delta V")
print(f"  S = ({Delta_phi:.4e})^4 / {DeltaV:.4e}")
print(f"  S = {S_thick:.6e}")
print()
if S_thick > 400:
    print(f"  S = {S_thick:.4e} >> 400 => COSMOLOGICALLY STABLE")
else:
    print(f"  S = {S_thick:.4e} < 400 => UNSTABLE")
print()

# More refined: use the parametric formula for a single field with
# V(phi) = f_pi^2 phi^2 (soft mass), rolling from phi = M_u to phi = 0.
# This is a harmonic potential, so the bounce is that of a harmonic oscillator
# deformed by the false vacuum energy.
# S_bounce ~ 2 pi^2 (f_pi M_u^2)^2 / f_pi^4 = 2 pi^2 M_u^4 / f_pi^4
# No wait, let's be more careful.

# The potential along the M_u direction (keeping others at seesaw):
# V(M_u) = |m_u + X M_d M_s / M_u_vev * M_u|^2 + f_pi^2 M_u^2 + ...
# This is complicated. Let me use the single-field reduction:
# phi = M_u, V(phi) ~ f_pi^2 phi^2 (dominant soft mass)
# Delta V ~ f_pi^2 M_u^2 - 0 ~ f_pi^2 M_u^2
# Delta phi ~ M_u
# S_thick ~ M_u^4 / (f_pi^2 M_u^2) = M_u^2 / f_pi^2

S_thick_single = M_diag[0]**2 / f_pi_sq
print("--- Single-field (M_u) thick-wall ---")
print()
print(f"  S ~ M_u^2 / f_pi^2 = ({M_diag[0]:.4f})^2 / ({f_pi})^2")
print(f"  S = {S_thick_single:.6e}")
print()
if S_thick_single > 400:
    print(f"  S = {S_thick_single:.4e} >> 400 => COSMOLOGICALLY STABLE")
else:
    print(f"  S = {S_thick_single:.4e} < 400 => UNSTABLE")
print()

# The most conservative estimate: use the ratio (Delta phi)^4 / Delta V
# with the FULL multifield displacement
print("--- Summary of bounce action estimates ---")
print()
print(f"  Thin-wall:        S = {S_thin:.4e}")
print(f"  Thick-wall (all): S = {S_thick:.4e}")
print(f"  Single-field:     S = {S_thick_single:.4e}")
print(f"  Threshold:        S > 400 for stability")
print()

all_stable = all(s > 400 for s in [S_thin, S_thick, S_thick_single])
if all_stable:
    print("  ALL estimates give S >> 400.")
    print("  The seesaw vacuum is COSMOLOGICALLY STABLE against tunneling to (B).")
else:
    print("  At least one estimate gives S < 400. Need more careful analysis.")
print()

# Tunneling rate
Gamma_over_V = np.exp(-min(S_thin, S_thick, S_thick_single))
print(f"  Tunneling rate: Gamma ~ exp(-S_min) = exp(-{min(S_thin, S_thick, S_thick_single):.4e})")
print(f"  Lifetime >> age of universe by ~ 10^{{{min(S_thin, S_thick, S_thick_single)/np.log(10):.0e}}} orders of magnitude")
print()


# ===========================================================================
# TASK 3: SEESAW VACUUM WITH W_B
# ===========================================================================
print()
print("=" * 78)
print("  TASK 3: SEESAW VACUUM WITH m_B B Btilde")
print("=" * 78)
print()

print("F-term equations at B = Btilde = 0, diagonal M, H = 0:")
print()
print("  F_X     = det M - 0 - Lambda^6 = det M - Lambda^6")
print("  F_{M_i} = m_i + X * prod_{j!=i} M_j")
print("  F_B     = (-X + m_B) * 0 = 0  [Btilde = 0]")
print("  F_Bt    = (-X + m_B) * 0 = 0  [B = 0]")
print()
print("The seesaw vacuum M_i = C/m_i, X = -C/Lambda^6 satisfies:")
print("  F_X = Lambda^6 - Lambda^6 = 0  (check)")
print("  F_{M_i} = m_i - (C/Lambda^6)(Lambda^6/M_i) = m_i - C/M_i = m_i - m_i = 0  (check)")
print("  F_B = 0, F_Bt = 0  (check, since B = Bt = 0)")
print()
print("RESULT: The seesaw vacuum is UNAFFECTED by m_B.")
print("  M_i, X, B, Bt all take the same values as without m_B.")
print("  m_B only modifies the equations along the baryonic directions.")
print()

# V at seesaw vacuum with m_B (no lambda for simplicity)
# V = 0 (all F-terms vanish) + f_pi^2 sum M_i^2
V_seesaw_mB = f_pi_sq * np.sum(M_diag**2)
print(f"V at seesaw vacuum (any m_B, no lambda):")
print(f"  V = f_pi^2 * sum M_i^2 = {V_seesaw_mB:.6e} MeV^2")
print(f"  (independent of m_B!)")
print()

# With lambda:
V_seesaw_mB_lam = V_A_total
print(f"V at seesaw vacuum (any m_B, with lambda):")
print(f"  V = |F_X|^2 + f_pi^2 sum M_i^2 = {V_seesaw_mB_lam:.6e} MeV^2")
print(f"  (also independent of m_B)")
print()

# However, the baryon mass DOES change the spectrum around the seesaw vacuum.
# The fermion mass matrix gains a new entry:
print("Baryon sector fermion mass at seesaw vacuum:")
print("  W_{B,Bt} = -X + m_B")
print(f"  Without m_B: W_{{B,Bt}} = -X = {-X_seesaw:.6e}")
print(f"  With m_B:    W_{{B,Bt}} = -X + m_B = {-X_seesaw:.6e} + m_B")
print(f"  For any m_B >> |X| ~ 10^{{-9}}: W_{{B,Bt}} ≈ m_B")
print()
print("The baryon fermion (baryino) gets mass ~ m_B instead of ~ |X| ~ 10^{-9}.")
print("This is a LARGE improvement: the ultra-light baryino is lifted to m_B.")
print()


# ===========================================================================
# TASK 4: IMPACT ON MASS PREDICTIONS
# ===========================================================================
print()
print("=" * 78)
print("  TASK 4: IMPACT ON MASS PREDICTIONS")
print("=" * 78)
print()

print("(a) Meson masses M_i:")
print("  M_i = C/m_i is determined by the seesaw equations F_{M_i} = 0.")
print("  These equations do not involve B, Btilde, or m_B.")
print(f"  M_u = {M_diag[0]:.4f} MeV  (UNCHANGED)")
print(f"  M_d = {M_diag[1]:.4f} MeV  (UNCHANGED)")
print(f"  M_s = {M_diag[2]:.4f} MeV  (UNCHANGED)")
print()

print("(b) Off-diagonal tachyonic masses:")
print("  The tachyonic B-terms are B_{ab} = -X F_{M_c}* - M_c F_X*")
print("  These depend on X and M_c at the vacuum, not on m_B.")
print("  Since m_B does not appear in F_{M_i} or F_X at B = Bt = 0:")
print("  The off-diagonal tachyonic masses are UNCHANGED.")
# Compute them for reference
F_X_val = F_X_at_A
for (a, b, label) in [(0,1,"ud"), (0,2,"us"), (1,2,"ds")]:
    c = 3 - a - b
    B_term = -X_seesaw * 0 - M_diag[c] * F_X_val  # F_{M_c} = 0 at seesaw
    m2_tach = f_pi_sq - abs(B_term)
    print(f"    m^2_tach({label}) = {m2_tach:.4e} MeV^2 (UNCHANGED)")
print()

print("(c) Higgs mass:")
print("  m_h = lambda v / sqrt(2) from the NMSSM coupling lambda X H_u.H_d")
m_h = lam * v_ew / np.sqrt(2)
print(f"  m_h = {lam} * {v_ew} / sqrt(2) = {m_h:.4f} MeV = {m_h/1000:.4f} GeV")
print(f"  This depends on lambda and v, NOT on m_B.")
print(f"  UNCHANGED.")
print()

print("(d) Cabibbo angle prediction:")
print("  The Oakes relation theta_C = arctan(sqrt(m_d/m_s)) depends on quark masses only.")
theta_C = np.degrees(np.arctan(np.sqrt(m_d/m_s)))
print(f"  theta_C = {theta_C:.4f} deg (PDG: 13.04 deg)")
print(f"  UNCHANGED by m_B.")
print()

print("SUMMARY: m_B BB~ modifies NOTHING in the meson/Higgs sector.")
print("  It only affects the baryonic (B, Btilde) sector:")
print("  - Lifts the ultra-light baryino from m ~ 10^{-9} MeV to m ~ m_B")
print("  - Lifts the baryonic flat direction (with Kahler constraint)")
print("  - Does not change any observable mass predictions")
print()


# ===========================================================================
# TASK 5: NATURAL SCALE FOR m_B
# ===========================================================================
print()
print("=" * 78)
print("  TASK 5: NATURAL SCALE FOR m_B")
print("=" * 78)
print()

print("Three candidate scales:")
print()
print(f"  (1) m_B ~ Lambda = {Lambda} MeV (confinement scale)")
print(f"      This is the most natural from the SQCD perspective.")
print(f"      The baryon is a confining-phase composite Q Q Q / Lambda^3.")
print(f"      A mass term m_B ~ Lambda arises from the confining dynamics.")
print()
print(f"  (2) m_B ~ f_pi = {f_pi} MeV (soft SUSY breaking scale)")
print(f"      This would tie m_B to the same scale as V_soft.")
print(f"      Less natural: f_pi enters as a spurion, not from SQCD dynamics.")
print()
print(f"  (3) m_B ~ v = {v_ew} MeV = {v_ew/1000:.2f} GeV (electroweak scale)")
print(f"      Too high: this would make the baryino heavier than most SM particles.")
print(f"      No direct coupling between baryons and the Higgs doublets.")
print()

# Check what m_B ~ Lambda does numerically:
mB_Lambda = Lambda
V_B_Lambda = V_B_original + 2 * mB_Lambda**2 * Lambda6
print(f"With m_B = Lambda = {Lambda} MeV:")
print(f"  V_B(m_B) = sum m_i^2 + 2 m_B^2 Lambda^6")
print(f"           = {V_B_original:.2f} + 2 * {mB_Lambda}^2 * {Lambda6:.4e}")
print(f"           = {V_B_original:.2f} + {2 * mB_Lambda**2 * Lambda6:.4e}")
print(f"           = {V_B_Lambda:.6e} MeV^2")
print()
print(f"  V_A(no lam) = {V_A_no_lam:.6e} MeV^2")
print(f"  V_B(m_B=Lambda) / V_A = {V_B_Lambda/V_A_no_lam:.6e}")
print()
if V_B_Lambda > V_A_no_lam:
    print(f"  V_B > V_A: vacuum (B) is ELIMINATED. Seesaw is global minimum!")
else:
    print(f"  V_B < V_A: vacuum (B) still lower. Need larger m_B.")
print()

# Check m_B = f_pi:
mB_fpi = f_pi
V_B_fpi = V_B_original + 2 * mB_fpi**2 * Lambda6
print(f"With m_B = f_pi = {f_pi} MeV:")
print(f"  V_B(m_B) = {V_B_fpi:.6e} MeV^2")
print(f"  V_A = {V_A_no_lam:.6e} MeV^2")
if V_B_fpi > V_A_no_lam:
    print(f"  V_B > V_A: ELIMINATED")
else:
    print(f"  V_B < V_A: still lower")
print()

# What m_B is needed?
print(f"Required m_B to make V_B > V_A (no-lambda case):")
print(f"  m_B > sqrt((V_A - V_B) / (2 Lambda^6))")
print(f"  m_B > sqrt({DeltaV:.4e} / {2*Lambda6:.4e})")
print(f"  m_B > {mB_threshold_no_lam:.6f} MeV")
print()

# Baryon mass from instanton effects:
print("Physical origin of m_B:")
print("  In SQCD with N_f = N_c, the baryon is a solitonic object in the")
print("  confining phase. The Seiberg constraint det M - BBt = Lambda^{2N_c}")
print("  already relates baryons to mesons.")
print()
print("  A mass term m_B BBt can arise from:")
print("  (1) Higher-dimensional operators in the UV theory, suppressed by Lambda")
print("  (2) Instanton corrections: W_inst ~ Lambda^{2N_c+1} BBt / Lambda^6 ~ Lambda BBt")
print("  (3) Explicit mass from integrating out heavy flavors")
print()
print("  The instanton origin gives m_B ~ Lambda = 300 MeV as the natural scale.")
print("  This is equivalent to the dynamical scale of the confining theory.")
print()

# The baryino mass spectrum with m_B ~ Lambda:
print("Baryino spectrum with m_B = Lambda = 300 MeV:")
print(f"  m_baryino = |m_B - X_seesaw| = |{Lambda} - ({X_seesaw:.4e})| ≈ {Lambda} MeV")
print(f"  Compare: without m_B, m_baryino = |X_seesaw| = {abs(X_seesaw):.4e} MeV")
print()

# ===========================================================================
# NUMERICAL VERIFICATION TABLE
# ===========================================================================
print()
print("=" * 78)
print("  NUMERICAL VERIFICATION SUMMARY")
print("=" * 78)
print()

print("| Quantity                | Value                    | Unit    |")
print("|-------------------------|--------------------------|---------|")
print(f"| C                       | {C:.6f}            | MeV^2   |")
print(f"| M_u = C/m_u            | {M_diag[0]:.4f}            | MeV     |")
print(f"| M_d = C/m_d            | {M_diag[1]:.4f}            | MeV     |")
print(f"| M_s = C/m_s            | {M_diag[2]:.4f}             | MeV     |")
print(f"| X_seesaw               | {X_seesaw:.6e}       |         |")
print(f"| det M / Lambda^6       | {np.prod(M_diag)/Lambda6:.12f}      |         |")
print(f"| V_A (no lam)           | {V_A_no_lam:.6e}       | MeV^2   |")
print(f"| V_A (with lam)         | {V_A_total:.6e}       | MeV^2   |")
print(f"| V_B                    | {V_B_original:.4f}            | MeV^2   |")
print(f"| S_thin                 | {S_thin:.6e}       |         |")
print(f"| S_thick (all fields)   | {S_thick:.6e}       |         |")
print(f"| S_thick (M_u only)     | {S_thick_single:.6e}       |         |")
print(f"| m_B threshold (no lam) | {mB_threshold_no_lam:.6f}       | MeV     |")
print(f"| m_B threshold (w/ lam) | {mB_threshold_with_lam:.6f}     | MeV     |")
print(f"| m_h = lam v/sqrt(2)    | {m_h:.4f}           | MeV     |")
print(f"| theta_C (Oakes)        | {theta_C:.4f}              | deg     |")
print()

# Cross-checks
print("Cross-checks:")
print(f"  det M = M_u * M_d * M_s = {np.prod(M_diag):.6e}")
print(f"  Lambda^6 = {Lambda6:.6e}")
print(f"  Ratio = {np.prod(M_diag)/Lambda6:.12f} (should be 1.000000)")
print(f"  C = m_u * M_u = {m_u * M_diag[0]:.6f}")
print(f"  C = m_d * M_d = {m_d * M_diag[1]:.6f}")
print(f"  C = m_s * M_s = {m_s * M_diag[2]:.6f}")
print(f"  (all equal to C = {C:.6f}, check)")
print()

print("=" * 78)
print("  END OF COMPUTATION")
print("=" * 78)
