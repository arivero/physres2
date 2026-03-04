#!/usr/bin/env python3
"""
Canonical rescaling of the Lagrange multiplier X in Seiberg SQCD + NMSSM.

The Lagrange multiplier X enforcing det M - BB~ - Lambda^6 = 0 has
non-canonical dimension [X] = [mass]^{-3}. The rescaling S = X * Lambda^4
defines a field with [S] = [mass], enabling a standard NMSSM coupling
lambda_new S H_u.H_d with dimensionless lambda_new.

This script computes:
  1. The rescaled VEV <S> = X_0 * Lambda^4
  2. The full superpotential in terms of S with explicit couplings
  3. The effective mu-term mu_eff = lambda_new <S>
  4. The Higgs quartic and tree-level m_h at tan(beta) = 1
  5. The Kahler potential and canonical metric in terms of S
  6. F-terms in both bases and the physical scalar potential
  7. Summary of physical (in)dependence on the field redefinition
"""

import numpy as np

# =========================================================================
# Parameters
# =========================================================================
m_u = 2.16       # MeV (MS-bar at 2 GeV)
m_d = 4.67       # MeV
m_s = 93.4       # MeV
m_c = 1270.0     # MeV
m_b = 4180.0     # MeV
LAM = 300.0      # MeV (confinement scale)
v   = 246220.0   # MeV (full EW VEV)

# Yukawa couplings at tan(beta) = 1
y_c = 2.0 * m_c / v
y_b = 2.0 * m_b / v

# NMSSM coupling (original, dimensionful)
lam_orig = 0.72  # [mass]^2 for W = lam_orig * X * (H_u . H_d) to have dim [mass]^3

# Derived
m_arr = np.array([m_u, m_d, m_s])
prod_m = np.prod(m_arr)
prod_m_cbrt = prod_m**(1.0/3.0)
C = LAM**2 * prod_m_cbrt        # MeV^2
LAM4 = LAM**4                   # MeV^4
LAM6 = LAM**6                   # MeV^6
LAM8 = LAM**8                   # MeV^8

# Meson VEVs (Seiberg seesaw)
M_vev = C / m_arr               # [M_u, M_d, M_s] in MeV^2
M_u, M_d, M_s = M_vev

# det <M> = product of diagonal VEVs
det_M = M_u * M_d * M_s         # MeV^6

# X VEV from F_{M_i} = 0: m_i + X * cofactor_i = 0
# -> X = -m_i / (M_j * M_k) for any i (all equivalent at the seesaw VEV)
X0 = -m_u / (M_d * M_s)         # MeV^{-3} -- note: negative
# Cross-check
X0_check1 = -m_d / (M_u * M_s)
X0_check2 = -m_s / (M_u * M_d)
X0_check3 = -C / LAM6           # = -(prod m)^{1/3} / Lambda^4

# Higgs VEVs at tan(beta) = 1
v_hu = v / np.sqrt(2)
v_hd = v / np.sqrt(2)

# =========================================================================
# Output lines
# =========================================================================
out = []

def p(s=""):
    out.append(s)
    print(s)


# =========================================================================
# TASK 0: Verify parameters and X VEV
# =========================================================================
p("=" * 76)
p("TASK 0: Parameter verification")
p("=" * 76)
p()
p(f"  Quark masses:  m_u = {m_u} MeV,  m_d = {m_d} MeV,  m_s = {m_s} MeV")
p(f"  Yukawas:       y_c = 2*m_c/v = {y_c:.10f}")
p(f"                 y_b = 2*m_b/v = {y_b:.10f}")
p(f"  Lambda = {LAM} MeV")
p(f"  Lambda^4 = {LAM4:.6e} MeV^4")
p(f"  Lambda^6 = {LAM6:.6e} MeV^6")
p(f"  C = Lambda^2 * (m_u m_d m_s)^{{1/3}} = {C:.6f} MeV^2")
p(f"  (m_u m_d m_s)^{{1/3}} = {prod_m_cbrt:.6f} MeV")
p()
p(f"  Meson VEVs (seesaw):")
p(f"    <M_u> = C/m_u = {M_u:.6f} MeV^2")
p(f"    <M_d> = C/m_d = {M_d:.6f} MeV^2")
p(f"    <M_s> = C/m_s = {M_s:.6f} MeV^2")
p()
p(f"  det<M> = M_u * M_d * M_s = {det_M:.6e} MeV^6")
p(f"  Lambda^6 = {LAM6:.6e} MeV^6")
p(f"  Constraint: det<M> - Lambda^6 = {det_M - LAM6:.6e} MeV^6")
p(f"    (Should vanish for SUSY vacuum with B = Btilde = 0.)")
p()

p(f"  X VEV (three equivalent computations):")
p(f"    X0 = -m_u / (M_d * M_s) = {X0:.10e} MeV^{{-3}}")
p(f"    X0 = -m_d / (M_u * M_s) = {X0_check1:.10e} MeV^{{-3}}")
p(f"    X0 = -m_s / (M_u * M_d) = {X0_check2:.10e} MeV^{{-3}}")
p(f"    X0 = -C / Lambda^6      = {X0_check3:.10e} MeV^{{-3}}")
p(f"    Consistency: max deviation = {max(abs(X0 - X0_check1), abs(X0 - X0_check2), abs(X0 - X0_check3)):.2e}")
p()


# =========================================================================
# TASK 1: Rescaled field S = X * Lambda^4
# =========================================================================
p("=" * 76)
p("TASK 1: Rescaled field S = X * Lambda^4")
p("=" * 76)
p()

S0 = X0 * LAM4

p("  Dimensional analysis:")
p("    [X] = [mass]^{-3}  (from [X * det M] = [mass]^3, [det M] = [mass]^6)")
p("    [S] = [X] * [Lambda^4] = [mass]^{-3} * [mass]^4 = [mass]")
p()
p(f"  <S> = X0 * Lambda^4 = ({X0:.10e}) * ({LAM4:.6e})")
p(f"      = {S0:.10e} MeV")
p(f"      = {S0:.6f} MeV")
p(f"      = {S0 / 1e3:.10f} GeV")
p()

# Check: S0 = -C / Lambda^2
S0_check = -C / LAM**2
p(f"  Cross-check: <S> = -C / Lambda^2 = {S0_check:.10e} MeV")
p(f"    (Since X0 = -C/Lambda^6, and S = X*Lambda^4 => S0 = -C/Lambda^2)")
p(f"    Match: {abs(S0 - S0_check) < 1e-10}")
p()

# Magnitude
p(f"  |<S>| = {abs(S0):.6f} MeV = {abs(S0)/1e3:.10f} GeV")
p()


# =========================================================================
# TASK 2: Full superpotential in terms of S
# =========================================================================
p("=" * 76)
p("TASK 2: Superpotential in terms of S")
p("=" * 76)
p()

p("  Original superpotential (in terms of X):")
p("    W = sum_i m_i M^i_i + X (det M - B Btilde - Lambda^6)")
p("        + y_c H_u^0 M^d_d + y_b H_d^0 M^s_s")
p("        + lambda X (H_u . H_d)")
p()
p("  The coupling lambda X (H_u.H_d) has dimensions:")
p("    [lambda] * [X] * [H_u] * [H_d] = [lambda] * [mass]^{-3} * [mass] * [mass]")
p("    = [lambda] * [mass]^{-1}")
p("  For [W] = [mass]^3, we need [lambda] = [mass]^4, i.e. lambda has dim [mass]^4.")
p("  With numerical value lambda_orig = 0.72 (stated as dimensionless in problem),")
p("  this means the problem statement's lambda implicitly absorbs Lambda^4.")
p()
p("  Substituting X = S / Lambda^4:")
p()
p("  W = sum_i m_i M^i_i + (S / Lambda^4)(det M - B Btilde - Lambda^6)")
p("      + y_c H_u^0 M^d_d + y_b H_d^0 M^s_s")
p("      + lambda_orig * (S / Lambda^4) * (H_u . H_d)")
p()
p("  Rewrite grouping by Lambda^4:")
p()
p("  W = sum_i m_i M^i_i")
p("      + S * (det M / Lambda^4 - B Btilde / Lambda^4 - Lambda^2)")
p("      + y_c H_u^0 M^d_d + y_b H_d^0 M^s_s")
p("      + lambda_new * S * (H_u . H_d)")
p()

# New coupling constants
lam_new = lam_orig / LAM4

p("  New coupling constants:")
p(f"    lambda_new = lambda_orig / Lambda^4 = {lam_orig} / {LAM4:.6e}")
p(f"              = {lam_new:.10e} MeV^{{-4}}")
p()
p("  WAIT -- this gives lambda_new with dimensions [mass]^{-4}.")
p("  Let us recheck: the term S * (H_u . H_d) has dimensions")
p("    [mass] * [mass] * [mass] = [mass]^3")
p("  So lambda_new must be DIMENSIONLESS for [W] = [mass]^3.")
p()
p("  The resolution: the original lambda = 0.72 already HAS dimension [mass]^4.")
p("  The problem statement says lambda = 0.72, but for W to work, this is")
p("  lambda_dim = 0.72 * Lambda^4 when Lambda is set to 1 in natural units.")
p("  Equivalently, if we interpret the problem as working in units of Lambda,")
p("  then lambda = 0.72 is dimensionless and Lambda = 1.")
p()
p("  INTERPRETATION A: lambda is dimensionless, Lambda is the unit.")
p("  Then lambda_new = lambda / Lambda^4 has dim [mass]^{-4}.")
p("  But S (H_u.H_d) has dim [mass]^3, so lambda_new * S * (H_u.H_d) has dim [mass]^{-1}.")
p("  This FAILS. So lambda CANNOT be dimensionless in the X basis.")
p()
p("  INTERPRETATION B: lambda has dim [mass]^2 (stated in the problem).")
p("  Then [lambda * X * H_u * H_d] = [mass]^2 * [mass]^{-3} * [mass]^2 = [mass]^1.")
p("  This also fails: we need [mass]^3.")
p()
p("  INTERPRETATION C (CORRECT): lambda has dim [mass]^4.")
p("  [lambda * X * H_u * H_d] = [mass]^4 * [mass]^{-3} * [mass]^2 = [mass]^3. Correct!")
p("  The problem's option (a) says [mass]^2 but the correct answer is [mass]^4.")
p()
p("  With lambda = lam_phys * Lambda^4 where lam_phys = 0.72:")
p("  W_NMSSM = lam_phys * Lambda^4 * X * (H_u.H_d)")
p("          = lam_phys * (X * Lambda^4) * (H_u.H_d)")
p("          = lam_phys * S * (H_u.H_d)")
p()
p("  So lambda_new = lam_phys = 0.72 (dimensionless).")
p()

# Set the correct interpretation
lam_new_correct = lam_orig  # dimensionless

p(f"  RESULT: lambda_new = {lam_new_correct}")
p()

p("  Similarly, the constraint term:")
p("    X (det M - BB~ - Lambda^6) = (S/Lambda^4)(det M - BB~ - Lambda^6)")
p("    = S (det M / Lambda^4) - S (BB~ / Lambda^4) - S Lambda^2")
p()
p("  The new couplings in the constraint sector are:")

kappa_det = 1.0 / LAM4     # coupling of S * det M
kappa_BB = 1.0 / LAM4      # coupling of S * B Btilde
kappa_tad = LAM**2          # linear tadpole coefficient, S * Lambda^2

p(f"    kappa_det  = 1/Lambda^4 = {kappa_det:.10e} MeV^{{-4}}  (S * det M coupling)")
p(f"    kappa_BB   = 1/Lambda^4 = {kappa_BB:.10e} MeV^{{-4}}  (S * B Btilde coupling)")
p(f"    kappa_tad  = Lambda^2   = {kappa_tad:.6e} MeV^2    (S tadpole = -S * Lambda^2)")
p()
p("  Complete superpotential:")
p("    W = sum_i m_i M^i_i")
p("        + S (det M / Lambda^4 - BB~/Lambda^4 - Lambda^2)")
p("        + y_c H_u^0 M^d_d + y_b H_d^0 M^s_s")
p("        + 0.72 * S * (H_u . H_d)")
p()

# Verify VEVs are consistent
# F_S = det<M>/Lambda^4 - Lambda^2 should vanish
FS_vev = det_M / LAM4 - LAM**2
p("  VEV consistency checks:")
p(f"    F_S = det<M>/Lambda^4 - Lambda^2 = {det_M/LAM4:.6e} - {LAM**2:.6e}")
p(f"        = {FS_vev:.6e} MeV^2  (should vanish)")
p()

# F_{M_i} = m_i + S * d(det M / Lambda^4)/dM_i = m_i + (S/Lambda^4) * cofactor_i
# cofactor_i = product of other two M_j, M_k
# At the seesaw VEV with S = S0:
for i, (name, mi, cof) in enumerate([
    ("u", m_u, M_d * M_s),
    ("d", m_d, M_u * M_s),
    ("s", m_s, M_u * M_d)
]):
    FM_i = mi + S0 * cof / LAM4
    p(f"    F_{{M_{name}}} = m_{name} + <S> * cofactor_{name}/Lambda^4")
    p(f"              = {mi} + ({S0:.6e}) * ({cof:.6e}) / ({LAM4:.6e})")
    p(f"              = {mi} + ({S0 * cof / LAM4:.10e})")
    p(f"              = {FM_i:.10e} MeV  (should vanish)")
    p()


# =========================================================================
# TASK 3: Effective mu-term
# =========================================================================
p("=" * 76)
p("TASK 3: Effective mu-term")
p("=" * 76)
p()

mu_eff = lam_new_correct * S0

p("  mu_eff = lambda_new * <S>")
p(f"         = {lam_new_correct} * ({S0:.10e} MeV)")
p(f"         = {mu_eff:.10e} MeV")
p(f"         = {mu_eff:.6f} MeV")
p(f"         = {mu_eff / 1e3:.10f} GeV")
p()

p(f"  |mu_eff| = {abs(mu_eff):.6f} MeV = {abs(mu_eff)/1e3:.10f} GeV")
p()

p("  Required for viable EWSB: |mu| ~ 100-1000 GeV")
p(f"  Actual: |mu_eff| = {abs(mu_eff)/1e3:.6f} GeV")
p()

if abs(mu_eff) / 1e3 < 1.0:
    p("  RESULT: |mu_eff| << 100 GeV. The mu-term is FAR TOO SMALL.")
    p()
    p("  This is because <S> = -C/Lambda^2, and C = Lambda^2 * (m_u m_d m_s)^{1/3}")
    p("  is set by light quark masses:")
    p(f"    C = {C:.4f} MeV^2")
    p(f"    <S> = -C/Lambda^2 = -{C:.4f}/{LAM**2:.0f} = {S0:.6f} MeV")
    p(f"    |<S>| = {abs(S0):.6f} MeV << m_Z = 91188 MeV")
    p()
    p("  The Seiberg seesaw VEV for X is set by the RATIO of the quark mass")
    p("  geometric mean to Lambda^2. Since (m_u m_d m_s)^{1/3} ~ 10 MeV")
    p("  and Lambda ~ 300 MeV, we get <S> ~ 10 MeV.")
    p()
    p("  To get mu_eff ~ 100 GeV, we would need either:")
    p("    (a) lambda_new ~ 100 GeV / |<S>| = {:.1f}  (non-perturbative)".format(
        100e3 / abs(S0)))
    p("    (b) A different source for the mu-term (Giudice-Masiero, etc.)")
    p("    (c) N_f != N_c = 3 (different Seiberg phase)")
    p("    (d) The NMSSM singlet S is NOT identified with the Seiberg X")
else:
    p("  RESULT: |mu_eff| is in the right ballpark for EWSB.")

p()


# =========================================================================
# TASK 4: Higgs quartic and tree-level m_h
# =========================================================================
p("=" * 76)
p("TASK 4: Higgs quartic and tree-level Higgs mass")
p("=" * 76)
p()

p("  In the NMSSM, the singlet coupling lambda S (H_u.H_d) generates")
p("  a quartic for the Higgs:")
p("    |F_S|^2 contains |lambda_new (H_u.H_d)|^2 = lambda_new^2 |H_u.H_d|^2")
p()
p("  At tan(beta) = 1 (v_u = v_d = v/sqrt(2)):")
p("    The D-term quartic vanishes along the D-flat direction.")
p("    The NMSSM singlet contributes:")
p("      Delta m_h^2 = lambda_new^2 * v^2 * sin^2(2 beta) / 2")
p("                  = lambda_new^2 * v^2 / 2   [since sin(2 beta) = 1]")
p()

mh2_NMSSM = lam_new_correct**2 * v**2 / 2
mh_NMSSM = np.sqrt(mh2_NMSSM)

p(f"  lambda_new = {lam_new_correct}")
p(f"  v = {v} MeV")
p(f"  Delta m_h^2 = {lam_new_correct}^2 * {v}^2 / 2 = {mh2_NMSSM:.6e} MeV^2")
p(f"  Delta m_h   = sqrt({mh2_NMSSM:.6e}) = {mh_NMSSM:.2f} MeV")
p(f"              = {mh_NMSSM/1e3:.4f} GeV")
p()

p("  Observed m_h = 125.25 GeV")
p(f"  Tree-level prediction: m_h = {mh_NMSSM/1e3:.4f} GeV")
p(f"  Ratio to observed: {mh_NMSSM / 125.25e3:.4f}")
p()

p("  KEY POINT: The Higgs mass formula depends only on lambda_new and v.")
p("  The rescaling X -> S does NOT change lambda_new (it is the same 0.72")
p("  in both bases once dimensional analysis is done correctly).")
p("  Therefore the Higgs quartic and m_h are UNCHANGED by the rescaling.")
p()

# What lambda would give 125 GeV?
lam_needed = np.sqrt(2) * 125.25e3 / v
p(f"  For m_h = 125.25 GeV: lambda_needed = sqrt(2) * m_h / v = {lam_needed:.6f}")
p(f"  Actual lambda_new = {lam_new_correct}")
p(f"  The tree-level Higgs mass is {mh_NMSSM/1e3:.2f} GeV, which is")
if mh_NMSSM / 1e3 > 125.25:
    p(f"  ABOVE the observed value. Radiative corrections would push it further up.")
else:
    deficit = 125.25e3 - mh_NMSSM
    p(f"  {deficit/1e3:.1f} GeV below the observed value.")
    p(f"  Stop loops or other radiative corrections must supply the remainder.")

p()


# =========================================================================
# TASK 5: Kahler potential in terms of S
# =========================================================================
p("=" * 76)
p("TASK 5: Kahler potential in terms of S")
p("=" * 76)
p()

p("  Original Kahler potential (for X):")
p("    K = |X|^2 - |X|^4 / (12 mu_K^2) + ...")
p()
p("  Substituting X = S / Lambda^4:")
p("    |X|^2 = |S|^2 / Lambda^8")
p("    |X|^4 = |S|^4 / Lambda^16")
p()
p("    K = |S|^2 / Lambda^8 - |S|^4 / (12 mu_K^2 Lambda^16) + ...")
p()
p("  The canonical Kahler metric is:")
p("    K_{S Sbar} = d^2 K / (dS dSbar)")
p("              = 1/Lambda^8 - |S|^2 / (3 mu_K^2 Lambda^16) + ...")
p()
p("  At the vacuum <S> = S0:")

# For the mu_K parameter, use a representative value
mu_K = LAM  # representative
mu_K2 = mu_K**2

K_SS_leading = 1.0 / LAM8
K_SS_correction = abs(S0)**2 / (3 * mu_K2 * LAM**16)
K_SS_total = K_SS_leading - K_SS_correction

p(f"    (Using mu_K = Lambda = {mu_K} MeV as representative)")
p()
p(f"    K_{{S Sbar}} = 1/Lambda^8 - |<S>|^2 / (3 mu_K^2 Lambda^16)")
p(f"               = {K_SS_leading:.10e}")
p(f"                 - {K_SS_correction:.10e}")
p(f"               = {K_SS_total:.10e} MeV^{{-8}}")
p()
p(f"    Ratio correction/leading = {K_SS_correction / K_SS_leading:.10e}")
p(f"    The correction is negligible: |<S>|^2 / (3 mu_K^2 Lambda^8) ~ {K_SS_correction / K_SS_leading:.2e}")
p()

p("  CRITICAL OBSERVATION: The Kahler metric for S is NOT canonical.")
p("  K_{S Sbar} = 1/Lambda^8, not 1.")
p("  The kinetic term is (1/Lambda^8) |dS|^2, not |dS|^2.")
p()
p("  To get a canonical kinetic term, define the TRULY canonical field:")
p("    S_hat = S / Lambda^4 = X")
p("  which takes us back to X! The field X, despite having non-standard mass")
p("  dimension, has a canonical Kahler potential K = |X|^2.")
p()
p("  ALTERNATIVELY, if we keep S, we must carry the factor 1/Lambda^8")
p("  through all physical computations. The F-term scalar potential becomes:")
p("    V = K^{S Sbar} |F_S|^2 + ... = Lambda^8 |F_S|^2 + ...")
p()
p("  Is the pseudo-modulus pole preserved?")
p("  The quartic term -|S|^4/(12 mu_K^2 Lambda^16) in the Kahler potential")
p("  gives a contribution to K_{S Sbar} proportional to -|S|^2.")
p("  At the zero of K_{S Sbar}:")
p("    1/Lambda^8 = |S_pole|^2 / (3 mu_K^2 Lambda^16)")
p("    |S_pole|^2 = 3 mu_K^2 Lambda^8")
p("    |S_pole| = sqrt(3) mu_K Lambda^4")
p()

S_pole = np.sqrt(3) * mu_K * LAM4
X_pole = np.sqrt(3) * mu_K  # the original X pole

p(f"    |S_pole| = sqrt(3) * mu_K * Lambda^4")
p(f"            = sqrt(3) * {mu_K} * {LAM4:.6e}")
p(f"            = {S_pole:.6e} MeV")
p()
p(f"    In terms of X: |X_pole| = |S_pole|/Lambda^4 = sqrt(3) * mu_K = {X_pole:.4f} MeV")
p(f"    (This is the same pole location as in the X formulation.)")
p()
p(f"    |<S>| / |S_pole| = {abs(S0) / S_pole:.10e}")
p(f"    |<X>| / |X_pole| = {abs(X0) / (np.sqrt(3) * mu_K):.10e}")
p(f"    The VEV is far from the Kahler pole. The pseudo-modulus structure is preserved.")
p()


# =========================================================================
# TASK 6: F-terms and physical scalar potential
# =========================================================================
p("=" * 76)
p("TASK 6: F-terms in both bases")
p("=" * 76)
p()

p("  In N=1 SUGRA / global SUSY, the scalar potential is:")
p("    V = K^{I Jbar} F_I Fbar_J")
p("  where F_I = dW/dPhi^I (the holomorphic F-term, NO Kahler metric factor).")
p()
p("  In the X basis:")
p("    F_X = dW/dX = det M - BB~ - Lambda^6 + lambda_dim (H_u.H_d)")
p("  At the vacuum (B = Btilde = 0, <M> = seesaw, <H> ~ v):")
p()

# F_X at vacuum (ignoring Higgs contribution for the constraint piece)
FX_constraint = det_M - LAM6
FX_Higgs = lam_orig * LAM4 * v_hu * v_hd  # lambda_dim = lam_orig * Lambda^4

p(f"    F_X (constraint) = det<M> - Lambda^6 = {FX_constraint:.6e} MeV^6")
p(f"    F_X (Higgs)      = lambda_dim * <H_u.H_d> = {lam_orig}*{LAM4:.2e}*{v_hu:.0f}*{v_hd:.0f}")
p(f"                     = {FX_Higgs:.6e} MeV^6")
p(f"    F_X (total)      = {FX_constraint + FX_Higgs:.6e} MeV^6")
p()

# In the Kahler potential K = |X|^2, the metric is K_{X Xbar} = 1
# So K^{X Xbar} = 1
# V_X = K^{X Xbar} |F_X|^2 = |F_X|^2

V_X = abs(FX_constraint + FX_Higgs)**2

p(f"    K_{{X Xbar}} = 1 (canonical)")
p(f"    V_X = K^{{X Xbar}} |F_X|^2 = |F_X|^2 = {V_X:.6e} MeV^12")
p()

p("  In the S basis:")
p("    F_S = dW/dS = det M / Lambda^4 - BB~/Lambda^4 - Lambda^2 + lambda_new (H_u.H_d)")
p()

FS_constraint = det_M / LAM4 - LAM**2
FS_Higgs = lam_new_correct * v_hu * v_hd

p(f"    F_S (constraint) = det<M>/Lambda^4 - Lambda^2 = {FS_constraint:.6e} MeV^2")
p(f"    F_S (Higgs)      = lambda_new * <H_u.H_d> = {lam_new_correct}*{v_hu:.0f}*{v_hd:.0f}")
p(f"                     = {FS_Higgs:.6e} MeV^2")
p(f"    F_S (total)      = {FS_constraint + FS_Higgs:.6e} MeV^2")
p()

# The Kahler metric for S is K_{S Sbar} = 1/Lambda^8
# So K^{S Sbar} = Lambda^8
# V_S = K^{S Sbar} |F_S|^2 = Lambda^8 |F_S|^2

V_S = LAM8 * abs(FS_constraint + FS_Higgs)**2

p(f"    K_{{S Sbar}} = 1/Lambda^8 = {1/LAM8:.10e} MeV^{{-8}}")
p(f"    K^{{S Sbar}} = Lambda^8 = {LAM8:.6e} MeV^8")
p(f"    V_S = K^{{S Sbar}} |F_S|^2 = Lambda^8 * |F_S|^2 = {V_S:.6e} MeV^12")
p()

# Verify equality
p("  VERIFICATION: V_X = V_S ?")
p(f"    V_X = {V_X:.10e} MeV^12")
p(f"    V_S = {V_S:.10e} MeV^12")
p(f"    Ratio V_X / V_S = {V_X / V_S if V_S != 0 else float('inf'):.15f}")
p()

# The key identity: F_S = F_X / Lambda^4
# K^{S Sbar} = Lambda^8 = (Lambda^4)^2
# So V_S = Lambda^8 * |F_X/Lambda^4|^2 = Lambda^8 * |F_X|^2 / Lambda^8 = |F_X|^2 = V_X

p("  Algebraic proof:")
p("    F_S = dW/dS = d(W(X=S/Lambda^4))/dS = (dW/dX) * (dX/dS) = F_X / Lambda^4")
p("    K_{S Sbar} = (dX/dS)(dXbar/dSbar) * K_{X Xbar} = K_{X Xbar} / Lambda^8 = 1/Lambda^8")
p("    K^{S Sbar} = Lambda^8")
p("    V_S = Lambda^8 * |F_X / Lambda^4|^2 = Lambda^8 * |F_X|^2 / Lambda^8 = |F_X|^2 = V_X")
p()
p("  CONFIRMED: The scalar potential is INVARIANT under the field redefinition.")
p()


# Also check F_M terms
p("  F-term cross-check for mesons:")
p("    In X basis: F_{M_i} = m_i + X * cofactor_i")
p("    In S basis: F_{M_i} = m_i + (S/Lambda^4) * cofactor_i")
p("    These are IDENTICAL (S/Lambda^4 = X by definition).")
p("    The meson F-terms do not change at all.")
p()


# =========================================================================
# TASK 7: Summary
# =========================================================================
p("=" * 76)
p("TASK 7: Summary — does X -> S change any physical observable?")
p("=" * 76)
p()

p("  The rescaling S = X * Lambda^4 is a holomorphic field redefinition.")
p("  Under such a redefinition:")
p()
p("  1. SUPERPOTENTIAL: Changes form but not physics. W(X) -> W(S/Lambda^4).")
p("     The equations of motion dW/dPhi = 0 are equivalent in both bases.")
p()
p("  2. VACUUM: <S> = Lambda^4 <X>. The location changes but the physics doesn't.")
p(f"     <X> = {X0:.6e} MeV^{{-3}},  <S> = {S0:.6f} MeV")
p()
p("  3. COUPLING CONSTANTS: lambda_new = lambda_orig / Lambda^4 in the naive reading.")
p("     But lambda_orig was ALREADY dimensionful ([mass]^4), so lambda_new = 0.72")
p("     (dimensionless). No physical coupling changes.")
p()
p("  4. MU-TERM: mu_eff = lambda_new <S> = 0.72 * ({:.6f} MeV) = {:.6f} MeV".format(
    S0, mu_eff))
p(f"     |mu_eff| = {abs(mu_eff)/1e3:.6f} GeV")
p("     This is the SAME as lambda_dim <X> = (0.72 * Lambda^4) * <X>")
val = lam_orig * LAM4 * X0
p(f"     = {lam_orig} * {LAM4:.6e} * ({X0:.6e}) = {val:.6f} MeV")
p(f"     Confirmed equal: |difference| = {abs(mu_eff - val):.2e} MeV")
p()
p("  5. HIGGS MASS: m_h^2 = lambda_new^2 v^2 / 2 depends on lambda_new and v.")
p("     Since lambda_new = 0.72 in both bases, m_h is UNCHANGED.")
p(f"     m_h (tree) = {mh_NMSSM/1e3:.4f} GeV")
p()
p("  6. F-TERMS: F_S = F_X / Lambda^4. The PHYSICAL potential V = K^{IJ} F_I Fbar_J")
p("     compensates via K^{S Sbar} = Lambda^8. V is exactly invariant.")
p()
p("  7. KAHLER METRIC: K_{S Sbar} = 1/Lambda^8 (non-canonical). The pseudo-modulus")
p("     pole at |X_pole| = sqrt(3) mu_K maps to |S_pole| = sqrt(3) mu_K Lambda^4.")
p("     The pole structure is preserved; only the coordinate label changes.")
p()
p("  CONCLUSION: The rescaling X -> S = X Lambda^4 changes NO physical observable.")
p("  It is a pure holomorphic field redefinition. Masses, mixing angles, F-terms,")
p("  and the scalar potential are all invariant.")
p()
p("  The only advantage of the S basis is NOTATIONAL: it avoids writing a field")
p("  with non-standard mass dimension, making the NMSSM analogy more transparent.")
p("  But the price is a non-canonical Kahler metric K_{S Sbar} = 1/Lambda^8.")
p()
p("  The mu-problem remains: |mu_eff| = {:.3f} MeV = {:.6f} GeV,".format(
    abs(mu_eff), abs(mu_eff)/1e3))
p("  which is three orders of magnitude below the electroweak scale.")
p("  This is a PHYSICAL fact, independent of the field basis.")
p()

# =========================================================================
# Numerical summary table
# =========================================================================
p("=" * 76)
p("NUMERICAL SUMMARY")
p("=" * 76)
p()
p("  Quantity                    X basis                    S basis")
p("  " + "-" * 72)
p(f"  Field dimension            [mass]^{{-3}}               [mass]")
p(f"  VEV                        {X0:.6e} MeV^{{-3}}    {S0:.6f} MeV")
p(f"  Kahler metric              1 (canonical)              1/Lambda^8 = {1/LAM8:.4e} MeV^{{-8}}")
p(f"  Inverse metric             1                          Lambda^8 = {LAM8:.4e} MeV^8")
p(f"  NMSSM coupling             0.72 * Lambda^4            0.72 (dimensionless)")
p(f"  F-term                     {FX_constraint + FX_Higgs:.4e} MeV^6        {FS_constraint + FS_Higgs:.4e} MeV^2")
p(f"  Scalar potential V         {V_X:.4e} MeV^12       {V_S:.4e} MeV^12")
p(f"  mu_eff                     {mu_eff:.6f} MeV           {mu_eff:.6f} MeV")
p(f"  m_h (tree, NMSSM)          {mh_NMSSM/1e3:.4f} GeV            {mh_NMSSM/1e3:.4f} GeV")
p(f"  Kahler pole at             |X| = {X_pole:.4f} MeV            |S| = {S_pole:.4e} MeV")
p()

p("Done.")


# =========================================================================
# Write markdown
# =========================================================================
md = []

md.append("# Canonical Rescaling of the Lagrange Multiplier X")
md.append("")
md.append("## Setup")
md.append("")
md.append("In SQCD with N_f = N_c = 3 (flavors u, d, s), the confined-phase superpotential is")
md.append("")
md.append("    W = sum_i m_i M^i_i + X(det M - BB~ - Lambda^6) + y_c H_u M^d_d + y_b H_d M^s_s + lambda X (H_u.H_d)")
md.append("")
md.append("The Lagrange multiplier X has [X] = [mass]^{-3}. The NMSSM coupling lambda X (H_u.H_d)")
md.append("requires lambda to have dimension [mass]^4 for W to have correct dimension [mass]^3.")
md.append("The problem statement's option (a) says [mass]^2 -- this is incorrect; the correct answer is [mass]^4.")
md.append("")
md.append("## 1. Rescaled Field")
md.append("")
md.append("Define S = X Lambda^4, so [S] = [mass]^{-3} [mass]^4 = [mass].")
md.append("")
md.append(f"    <S> = X_0 Lambda^4 = ({X0:.6e}) * ({LAM4:.6e}) = {S0:.6f} MeV")
md.append("")
md.append(f"Since X_0 = -C/Lambda^6 and C = Lambda^2 (m_u m_d m_s)^{{1/3}} = {C:.4f} MeV^2:")
md.append(f"    <S> = -C/Lambda^2 = -{C:.4f}/{LAM**2:.0f} = {S0:.6f} MeV")
md.append("")
md.append("## 2. Superpotential in Terms of S")
md.append("")
md.append("    W = sum_i m_i M^i_i")
md.append("        + S (det M / Lambda^4 - BB~/Lambda^4 - Lambda^2)")
md.append("        + y_c H_u M^d_d + y_b H_d M^s_s")
md.append("        + 0.72 S (H_u.H_d)")
md.append("")
md.append("New couplings:")
md.append(f"- Constraint: S couples to det M / Lambda^4 (dim [mass]^{{-4}} = {kappa_det:.4e})")
md.append(f"- Tadpole: S has a linear term -S Lambda^2 = -S * {kappa_tad:.0f} MeV^2")
md.append("- NMSSM: lambda_new = 0.72 (dimensionless, since lambda_orig was [mass]^4)")
md.append("")
md.append("## 3. Effective mu-Term")
md.append("")
md.append(f"    mu_eff = lambda_new <S> = 0.72 * ({S0:.6f}) = {mu_eff:.6f} MeV")
md.append(f"           = {abs(mu_eff)/1e3:.6f} GeV")
md.append("")
md.append("This is three orders of magnitude below the electroweak scale (~100 GeV).")
md.append("The smallness is physical: <S> = -(m_u m_d m_s)^{1/3} / Lambda^2,")
md.append("which is set by light quark masses over the confinement scale.")
md.append("")
md.append("## 4. Higgs Quartic and Tree-Level Mass")
md.append("")
md.append("At tan(beta) = 1, the NMSSM singlet contribution to the Higgs mass is:")
md.append("")
md.append(f"    m_h^2 = lambda_new^2 v^2 sin^2(2 beta) / 2 = (0.72)^2 * ({v/1e3:.2f} GeV)^2 / 2")
md.append(f"    m_h = {mh_NMSSM/1e3:.4f} GeV")
md.append("")
md.append("The rescaling does NOT change this. lambda_new = 0.72 in both bases")
md.append("(the original lambda was already [mass]^4 = 0.72 * Lambda^4, giving 0.72 after")
md.append("absorbing the Lambda^4 into S).")
md.append("")
md.append("For m_h = 125.25 GeV, we need lambda = sqrt(2) * 125.25/246.22 = {:.4f}.".format(lam_needed))
md.append("The value 0.72 gives {:.1f} GeV at tree level, {:.1f} GeV below the observed mass.".format(
    mh_NMSSM/1e3, (125.25e3 - mh_NMSSM)/1e3))
md.append("")
md.append("## 5. Kahler Potential")
md.append("")
md.append("Original: K = |X|^2 - |X|^4/(12 mu_K^2) + ...")
md.append("")
md.append("In terms of S:")
md.append("    K = |S|^2/Lambda^8 - |S|^4/(12 mu_K^2 Lambda^{16}) + ...")
md.append("")
md.append("    K_{S Sbar} = 1/Lambda^8 - |S|^2/(3 mu_K^2 Lambda^{16})")
md.append("")
md.append("The metric is NOT canonical (factor 1/Lambda^8). To get canonical kinetic terms,")
md.append("one must rescale back to X = S/Lambda^4, recovering the original basis.")
md.append("")
md.append("The pseudo-modulus pole is at:")
md.append(f"    |S_pole| = sqrt(3) mu_K Lambda^4 = {S_pole:.4e} MeV")
md.append(f"    |X_pole| = sqrt(3) mu_K = {X_pole:.4f} MeV  (same physical location)")
md.append("")
md.append(f"|<S>|/|S_pole| = {abs(S0)/S_pole:.4e} (VEV far from pole, consistent with perturbativity).")
md.append("")
md.append("## 6. F-Terms and Scalar Potential")
md.append("")
md.append("F_S = F_X / Lambda^4 (chain rule: dW/dS = (dW/dX)(dX/dS) = F_X / Lambda^4)")
md.append("")
md.append("The physical scalar potential compensates via the inverse Kahler metric:")
md.append("    V = K^{S Sbar} |F_S|^2 = Lambda^8 |F_X/Lambda^4|^2 = |F_X|^2")
md.append("")
md.append(f"Numerical verification: V_X = {V_X:.6e}, V_S = {V_S:.6e}, ratio = {V_X/V_S if V_S != 0 else 0:.15f}")
md.append("")
md.append("The scalar potential is EXACTLY invariant.")
md.append("")
md.append("## 7. Summary")
md.append("")
md.append("The rescaling X -> S = X Lambda^4 is a holomorphic field redefinition.")
md.append("It changes NO physical observable:")
md.append("")
md.append("| Quantity | Changed? | Reason |")
md.append("|----------|----------|--------|")
md.append("| Scalar potential V | No | K^{IJ} compensates F_I rescaling |")
md.append("| Physical masses | No | Eigenvalues of V'' are basis-independent |")
md.append("| Mixing angles | No | Diagonalization of mass matrix is basis-independent |")
md.append("| F-terms | Yes (by Lambda^4) | But K^{IJ} compensates, leaving V unchanged |")
md.append("| Kahler metric | Yes (1 -> 1/Lambda^8) | Non-canonical in S basis |")
md.append("| Coupling lambda | No | 0.72 in both bases (correct dim. analysis) |")
md.append("| mu_eff | No | {:.6f} MeV = {:.6f} GeV in both bases |".format(abs(mu_eff), abs(mu_eff)/1e3))
md.append("| m_h (tree) | No | {:.4f} GeV in both bases |".format(mh_NMSSM/1e3))
md.append("")
md.append("### The mu-Problem")
md.append("")
md.append("The effective mu-term |mu_eff| = {:.3f} MeV = {:.6f} GeV is far too small".format(
    abs(mu_eff), abs(mu_eff)/1e3))
md.append("for viable EWSB (requires ~100 GeV). This is a physical problem independent")
md.append("of the field basis. It arises because the Seiberg seesaw VEV <X> (or <S>)")
md.append("is set by light quark masses: <S> = -(m_u m_d m_s)^{1/3}/Lambda^2.")
md.append("")
md.append("To get mu_eff ~ 100 GeV, one would need either:")
md.append("- lambda >> 1 (non-perturbative)")
md.append("- A separate source for the mu-term (Giudice-Masiero mechanism)")
md.append("- The NMSSM singlet is not identified with the Seiberg Lagrange multiplier")
md.append("- A different SQCD phase (N_f != N_c = 3)")
md.append("")
md.append("### Numerical Summary Table")
md.append("")
md.append("| Parameter | Value |")
md.append("|-----------|-------|")
md.append(f"| Lambda | {LAM} MeV |")
md.append(f"| C = Lambda^2 (m_u m_d m_s)^{{1/3}} | {C:.4f} MeV^2 |")
md.append(f"| <X> | {X0:.6e} MeV^{{-3}} |")
md.append(f"| <S> = <X> Lambda^4 | {S0:.6f} MeV |")
md.append(f"| lambda_new | 0.72 (dimensionless) |")
md.append(f"| mu_eff = lambda_new <S> | {mu_eff:.6f} MeV |")
md.append(f"| m_h (tree, NMSSM) | {mh_NMSSM/1e3:.4f} GeV |")
md.append(f"| Kahler pole |S_pole| | {S_pole:.4e} MeV |")
md.append(f"| y_c | {y_c:.10f} |")
md.append(f"| y_b | {y_b:.10f} |")
md.append("")

md_text = "\n".join(md) + "\n"
md_path = "/home/codexssh/phys3/results/mu_term_analysis.md"
with open(md_path, "w") as f:
    f.write(md_text)

print(f"\n[Markdown written to {md_path}]")
