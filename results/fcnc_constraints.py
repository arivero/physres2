"""
FCNC Constraints on Off-Diagonal Meson Scalars in Confined SQCD
================================================================

Model: N=1 SUSY with SU(3) SQCD (N_f = N_c = 3), confined phase.
Low-energy spectrum: 9 meson scalars M^i_j and mesino partners psi^i_j.

The off-diagonal mesons (i != j) carry flavour quantum numbers and can
mediate FCNC at tree level (scalar exchange) or one loop (mesino box).

This script estimates:
  (a) Tree-level K-Kbar mixing from off-diagonal meson VEV
  (b) Mesino box-diagram contribution to K-Kbar mixing
  (c) Effective 4-quark Wilson coefficient vs UTfit bound
  (d) D-Dbar and B-Bbar constraints
  (e) Required M_FCNC scale
  (f) Mass at the true (stabilized) vacuum

All numbers in MeV unless otherwise stated.
"""

import numpy as np

# ═══════════════════════════════════════════════════════════════════════════
# Physical inputs
# ═══════════════════════════════════════════════════════════════════════════

# Quark masses (MS-bar at 2 GeV)
m_u = 2.16       # MeV
m_d = 4.67       # MeV
m_s = 93.4       # MeV
m_c = 1270.0     # MeV
m_b = 4180.0     # MeV
m_t = 172760.0   # MeV

# QCD / SQCD scales
Lambda_QCD = 300.0   # MeV
Lambda6 = Lambda_QCD**6
f_pi = 92.0          # MeV (pion decay constant)
m_tilde_sq = f_pi**2 # 8464 MeV^2  -- universal soft scalar mass^2

# Electroweak
v_ew = 246220.0  # MeV (full electroweak VEV, 246.22 GeV)
G_F = 1.1663788e-5 * 1e-6  # MeV^-2 (Fermi constant: 1.166e-5 GeV^-2)

# Seiberg seesaw VEVs
m_arr = np.array([m_u, m_d, m_s])
C = Lambda_QCD**2 * np.prod(m_arr)**(1.0/3.0)
M_vev = C / m_arr  # M_i = C/m_i
M_uu, M_dd, M_ss = M_vev

# X VEV
X0 = -C / Lambda6

# Yukawa couplings
y_c = 2.0 * m_c / v_ew
y_b = 2.0 * m_b / v_ew

# ═══════════════════════════════════════════════════════════════════════════
# Meson parameters
# ═══════════════════════════════════════════════════════════════════════════

# Kaon system
m_K = 497.611     # MeV (K0 mass)
f_K = 155.7       # MeV (kaon decay constant)
Dm_K_exp = 3.484e-12  # MeV (K_L - K_S mass difference)
eps_K_exp = 2.228e-3   # (indirect CP violation parameter)
# Bag parameter
B_K = 0.7625      # (lattice, RBC/UKQCD)

# D meson system
m_D = 1864.84     # MeV (D0 mass)
f_D = 212.0       # MeV (D decay constant)
x_D_exp = 0.0039  # (D0-D0bar mixing parameter x = Dm/Gamma)
Gamma_D = 1.0 / (410.1e-15 * 3e23 / 1e13)  # rough: tau_D = 410 fs => Gamma
# More careful: tau_D ~ 4.1e-13 s, convert to MeV^-1: hbar = 6.582e-22 MeV*s
hbar = 6.582119569e-22  # MeV*s
tau_D = 410.1e-15  # seconds
Gamma_D_mev = hbar / tau_D  # MeV
Dm_D_exp = x_D_exp * Gamma_D_mev  # MeV

# B_d system
m_Bd = 5279.66    # MeV
f_Bd = 190.0      # MeV
Dm_Bd_exp = 3.337e-10  # MeV (from Delta m_d = 0.5065 ps^-1)

# B_s system
m_Bs = 5366.88    # MeV
f_Bs = 230.3      # MeV
Dm_Bs_exp = 1.169e-8   # MeV (from Delta m_s = 17.749 ps^-1)

# UTfit Wilson coefficient bounds at 2 GeV
# |C1| bounds for Delta F = 2 operators
# From UTfit collaboration (arXiv:0707.0636, updated)
# C1 for K-Kbar: |Re C1| < 9.0e-13 GeV^-2 = 9.0e-7 MeV^-2
# C1 for D-Dbar: |C1| < 4.3e-13 GeV^-2
# C1 for B_d-Bbar_d: |C1| < 2.3e-11 GeV^-2
# C1 for B_s-Bbar_s: |C1| < 1.0e-9 GeV^-2

C1_bound_K  = 9.0e-13 * 1e6   # GeV^-2 -> MeV^-2 = 9.0e-7
C1_bound_D  = 4.3e-13 * 1e6   # MeV^-2
C1_bound_Bd = 2.3e-11 * 1e6   # MeV^-2
C1_bound_Bs = 1.0e-9  * 1e6   # MeV^-2

# ═══════════════════════════════════════════════════════════════════════════
# Output accumulator
# ═══════════════════════════════════════════════════════════════════════════
output_lines = []
md_lines = []

def p(s=""):
    output_lines.append(s)
    print(s)

def M(s=""):
    md_lines.append(s)

# ═══════════════════════════════════════════════════════════════════════════
# PART (a): Tree-level K-Kbar mixing from off-diagonal meson exchange
# ═══════════════════════════════════════════════════════════════════════════

p("=" * 76)
p("PART (a): Tree-level K-Kbar mixing from off-diagonal meson exchange")
p("=" * 76)
p()

# The off-diagonal meson M^d_s is a composite Q_tilde^d Q_s / Lambda.
# It couples to d and s quarks with strength proportional to 1/Lambda
# (or more precisely, through the superpotential coupling).
#
# At the diagonal seesaw vacuum, the effective 4-quark coupling from
# tree-level exchange of a scalar M^d_s with mass m_scalar is:
#
#   C_1^(tree) = g_eff^2 / (2 * m_scalar^2)
#
# where g_eff is the coupling at the d-s-scalar vertex.
#
# In the SQCD context, the meson M^d_s = Q_tilde^d Q_s / Lambda
# has a Yukawa coupling to quarks of order:
#   g_eff ~ Lambda / f_composite ~ 1
# (strongly coupled, g ~ 4pi from NDA)
#
# More precisely, from the Kahler potential in the magnetic theory:
#   K = M^dag M / Lambda^2
# which gives canonical normalization M_can = M / Lambda.
# The coupling is then g ~ 1 (from the superpotential linear in M).
#
# For the F-term contribution: the superpotential W = m_i M^i_i + ...
# gives a mass-insertion coupling:
#   vertex d_L s_R: ~ m_s (from F_{M^d_s} backreaction)
#
# However, the dominant coupling at the confinement scale is from
# the strong dynamics. Using NDA (naive dimensional analysis):
#   g_NDA = 4*pi / sqrt(N_c) ~ 4*pi / sqrt(3) ~ 7.3
#
# For a conservative estimate, use g_eff = 1 (perturbative)
# and g_eff = 4*pi (NDA) as upper bound.

# Scalar mass at the seesaw vacuum:
# From the Hessian, the off-diagonal meson masses^2 at the seesaw vacuum
# come from two sources:
#   1. The superpotential second derivative: W_{M^d_s, M^s_d} = X * M_u
#      This gives mass^2 = |X*M_u|^2 at the SUSY level
#   2. The soft term: m_tilde^2

# At the M=0 vacuum (true minimum from offdiag_vacuum analysis):
#   All meson scalars have mass^2 = 2*m_tilde^2 = 2*f_pi^2
#   (from the Hessian eigenvalues, see offdiag_vacuum.md)
m_scalar_M0 = np.sqrt(2 * m_tilde_sq)  # ~ 130 MeV

# At the seesaw vacuum:
#   The off-diagonal mesons have mass^2 from the F-term Hessian
#   plus the soft mass. The F-term mass is |X*M_k|^2 where k
#   is the third flavor.
#   For M^d_s and M^s_d: the third flavor is u, so
#   m_F^2 = (X0 * M_uu)^2 for the fermion mass,
#   For scalars: m_scalar^2 ~ |X0*M_uu| + m_tilde^2 (approximate)

# But at the seesaw vacuum the off-diagonal scalar masses come from
# the full scalar potential Hessian. The fermion mass for the (d,s)
# pair from W is |X0 * M_uu|:
m_mesino_ds = abs(X0) * M_uu
m_mesino_us = abs(X0) * M_dd
m_mesino_ud = abs(X0) * M_ss

p("Mesino masses (fermion partners of off-diagonal mesons) at seesaw vacuum:")
p(f"  m_mesino(d,s pair) = |X|*M_u = {m_mesino_ds:.6e} MeV")
p(f"  m_mesino(u,s pair) = |X|*M_d = {m_mesino_us:.6e} MeV")
p(f"  m_mesino(u,d pair) = |X|*M_s = {m_mesino_ud:.6e} MeV")
p()

p("Off-diagonal scalar masses at different vacua:")
p(f"  At M=0 vacuum:    m_scalar = sqrt(2)*f_pi = {m_scalar_M0:.2f} MeV")
p(f"  Soft mass alone:  m_scalar = f_pi = {f_pi:.0f} MeV")
p()

# The coupling of the off-diagonal meson to quarks:
# In the magnetic dual theory, the meson M^i_j = Q_tilde^i Q_j / Lambda
# couples to the elementary quarks through the quark-meson matching.
# The effective coupling at low energies (below Lambda) is:
#
# L_eff = (g / Lambda) * d_bar * s * M^d_s + h.c.
#
# where g ~ 4*pi (NDA) or g ~ 1 (perturbative matching).
# The canonical normalization of M gives:
# L_eff = g_eff * d_bar * s * M^d_s_can
# with g_eff ~ g * M_vev / Lambda (from expanding around the VEV)
# or simply g_eff ~ 1 from the Kahler term normalization.
#
# The simplest estimate: the 4-fermion operator from tree-level exchange
# of M^d_s with mass m_S and coupling g at each vertex:
#
#   C_1^(tree) = g^2 / (2 * m_S^2)
#
# This generates the effective Lagrangian:
#   L_eff = C_1 (d_bar_L gamma_mu s_L)^2

# Case 1: g_eff = 1 (perturbative, minimal coupling)
g_pert = 1.0

# Case 2: g_eff = 4*pi / sqrt(N_c) (NDA)
g_NDA = 4.0 * np.pi / np.sqrt(3.0)

# Case 3: g_eff from superpotential structure
# The superpotential coupling m_i M^i_i means the vertex is
# effectively g ~ m_q / M_vev (from expanding W around the VEV)
# For the d-s coupling: g ~ sqrt(m_d * m_s) / sqrt(M_dd * M_ss) (geometric mean)
g_supo = np.sqrt(m_d * m_s) / np.sqrt(M_dd * M_ss)

p("Effective coupling estimates for d-s-M^d_s vertex:")
p(f"  g_perturbative = {g_pert:.4f}")
p(f"  g_NDA = 4*pi/sqrt(3) = {g_NDA:.4f}")
p(f"  g_superpotential = sqrt(m_d*m_s)/sqrt(M_d*M_s) = {g_supo:.6e}")
p()

# The Wilson coefficient C1 from tree-level scalar exchange:
# C1 = g^2 / (2 * m_S^2)
# where the factor of 2 comes from the (LL)(LL) operator normalization.

# We compute C1 for various (g, m_S) combinations:

p("Tree-level Wilson coefficient C1 = g^2 / (2 * m_S^2):")
p()
p(f"{'g_eff':<20} {'m_S (MeV)':<16} {'C1 (MeV^-2)':<18} {'C1 (GeV^-2)':<18} {'C1/C1_bound':<14}")
p("-" * 86)

scenarios = [
    ("g=1, m=f_pi",           g_pert,  f_pi),
    ("g=1, m=sqrt(2)*f_pi",   g_pert,  m_scalar_M0),
    ("g=1, m=300",             g_pert,  300.0),
    ("g=1, m=1000",            g_pert,  1000.0),
    ("g=1, m=10000",           g_pert,  10000.0),
    ("g=NDA, m=f_pi",          g_NDA,   f_pi),
    ("g=NDA, m=sqrt(2)*f_pi",  g_NDA,   m_scalar_M0),
    ("g=NDA, m=300",            g_NDA,   300.0),
    ("g=NDA, m=1000",           g_NDA,   1000.0),
    ("g=NDA, m=10000",          g_NDA,   10000.0),
    ("g=supo, m=f_pi",          g_supo,  f_pi),
    ("g=supo, m=sqrt(2)*f_pi",  g_supo,  m_scalar_M0),
]

results_C1 = []
for label, g, mS in scenarios:
    C1 = g**2 / (2.0 * mS**2)
    C1_gev = C1 * 1e-6  # MeV^-2 -> GeV^-2
    ratio = C1 / C1_bound_K
    results_C1.append((label, g, mS, C1, C1_gev, ratio))
    p(f"  {label:<20} {mS:<16.2f} {C1:<18.6e} {C1_gev:<18.6e} {ratio:<14.2e}")

p()
p(f"UTfit bound on |Re(C1)| for K-Kbar: {C1_bound_K:.2e} MeV^-2 = {C1_bound_K*1e-6:.2e} GeV^-2")
p()

# K-Kbar mass difference contribution
# Dm_K^(NP) = 2 * Re(<Kbar| H_eff |K>) / (2 * m_K)
# where <Kbar| H_eff |K> = C1 * <Kbar| O1 |K>
# and <Kbar| O1 |K> = (8/3) * f_K^2 * m_K^2 * B_K
# (using the vacuum insertion approximation with bag factor)
#
# So: Dm_K^(NP) = C1 * (8/3) * f_K^2 * m_K * B_K

p("K-Kbar mass difference from new physics:")
p(f"  Dm_K^(NP) = C1 * (8/3) * f_K^2 * m_K * B_K")
p(f"  with f_K = {f_K} MeV, m_K = {m_K} MeV, B_K = {B_K}")
p()

hadronic_K = (8.0/3.0) * f_K**2 * m_K * B_K  # MeV^3
p(f"  Hadronic matrix element factor = {hadronic_K:.6e} MeV^3")
p()

p(f"{'Scenario':<24} {'Dm_K^NP (MeV)':<18} {'Dm_K^NP/Dm_K^exp':<20}")
p("-" * 62)
for label, g, mS, C1, C1_gev, ratio in results_C1:
    Dm_K_NP = C1 * hadronic_K
    ratio_DmK = Dm_K_NP / Dm_K_exp
    p(f"  {label:<24} {Dm_K_NP:<18.6e} {ratio_DmK:<20.2e}")

p()

# ═══════════════════════════════════════════════════════════════════════════
# PART (b): Mesino box-diagram contribution
# ═══════════════════════════════════════════════════════════════════════════

p("=" * 76)
p("PART (b): Mesino box-diagram contribution to K-Kbar mixing")
p("=" * 76)
p()

# Box diagram with mesino propagators:
#
#   d ---[g]--- mesino_1 ---[g]--- s
#   |                               |
#   s ---[g]--- mesino_2 ---[g]--- d
#
# The amplitude is:
#   C1^(box) = g^4 / (16*pi^2) * F(m1/M, m2/M)
# where m1, m2 are the mesino masses circulating in the box,
# M is the external momentum scale (~ m_K for K-Kbar),
# and F is the Inami-Lim function.
#
# For mesinos much heavier than m_K:
#   F(x1, x2) ~ 1 / (m_mesino^2) * [log terms]
#
# The mesino masses at the seesaw vacuum are:
#   m_mesino(d,s) = |X| * M_u ~ 4.94e-4 MeV  (extremely light!)
#
# At the M=0 vacuum: mesino masses come from the F-term mass matrix.
# With M=0, all W_{ij} involving off-diag mesons vanish (no cofactor terms).
# So mesino masses at M=0 come purely from the soft-breaking gaugino mass
# (not included in the superpotential) or are exactly zero.
#
# This is a PROBLEM: at the seesaw vacuum, mesinos are O(1e-4) MeV,
# and at M=0 they are massless. Both are extremely light.

p("Mesino masses at the seesaw vacuum:")
p(f"  m_mesino(d-s pair): {m_mesino_ds:.6e} MeV  (= |X|*M_u)")
p(f"  m_mesino(u-s pair): {m_mesino_us:.6e} MeV  (= |X|*M_d)")
p(f"  m_mesino(u-d pair): {m_mesino_ud:.6e} MeV  (= |X|*M_s)")
p()

# Inami-Lim box function for equal masses:
# S(x) = x * [1/(1-x)^2 + 3/(4(1-x)^3) * ln(x)] for x = m^2/M_W^2
# But here we are not in the SM; the box is between mesinos.
#
# For the general box with two internal masses m1, m2 and external quarks:
# When m_loop << m_external (which is the case here! mesino << m_K):
#   The box integral is dominated by the UV and gives:
#   C1^(box) ~ g^4 / (16*pi^2 * m_mesino^2)
#   (with log corrections from running)
#
# When m_loop >> m_external:
#   C1^(box) ~ g^4 / (16*pi^2 * m_loop^2)

# For the (d,s) sector at the seesaw vacuum:
# The mesino mass is m_mesino_ds ~ 5e-4 MeV << m_K
# So the box integral is ~ g^4/(16*pi^2 * m_mesino^2)
# This is HUGE because the mesino is so light.

# However, at the M=0 vacuum, there are no meson VEVs, so the mesino
# coupling to quarks also vanishes. The coupling goes as g ~ M_vev/Lambda
# at the seesaw vacuum, or g ~ VEV/Lambda.

# Let's compute with g = 1 and g = NDA, using the mesino mass from the
# seesaw vacuum:

p("Box-diagram Wilson coefficient C1^(box):")
p(f"  C1^(box) = g^4 / (16*pi^2 * m_mesino^2)")
p()

box_factor = 1.0 / (16.0 * np.pi**2)

for g_label, g_val in [("g=1", g_pert), ("g=NDA", g_NDA)]:
    for m_label, m_val in [("m_mesino(ds)", m_mesino_ds),
                           ("m_mesino(us)", m_mesino_us),
                           ("m_mesino(ud)", m_mesino_ud)]:
        if m_val > 0:
            C1_box = g_val**4 * box_factor / m_val**2
            ratio_box = C1_box / C1_bound_K
            p(f"  {g_label}, {m_label} = {m_val:.4e} MeV:")
            p(f"    C1^(box) = {C1_box:.4e} MeV^-2 = {C1_box*1e-6:.4e} GeV^-2")
            p(f"    C1^(box)/C1_bound = {ratio_box:.4e}")
p()

p("CRITICAL: The mesino box diagrams give C1 values that exceed the UTfit")
p("bound by many orders of magnitude if the mesinos are as light as the")
p("seesaw vacuum predicts (~5e-4 MeV). However:")
p("  1. At the M=0 vacuum, the effective coupling vanishes with the VEVs.")
p("  2. The mesino masses themselves depend on which vacuum is realized.")
p("  3. The relevant mass scale for the box is max(m_mesino, Lambda_QCD).")
p()

# ═══════════════════════════════════════════════════════════════════════════
# PART (c): Wilson coefficient comparison with UTfit bounds
# ═══════════════════════════════════════════════════════════════════════════

p("=" * 76)
p("PART (c): Tree-level Wilson coefficient vs UTfit bound")
p("=" * 76)
p()

# For off-diagonal meson scalars with mass^2 = f_pi^2 (= 8464 MeV^2):
m_S = f_pi  # 92 MeV
p(f"Off-diagonal meson scalar mass: m_S = f_pi = {m_S} MeV")
p(f"  m_S^2 = {m_S**2} MeV^2")
p()

# Effective 4-quark coupling from tree-level exchange at the kaon scale:
# C1 = g^2 / (2 * m_S^2)

for g_label, g_val in [("g=1", 1.0), ("g=4pi/sqrt(3) (NDA)", g_NDA),
                        ("g=0.1", 0.1), ("g=0.01", 0.01)]:
    C1_tree = g_val**2 / (2.0 * m_S**2)
    ratio_utfit = C1_tree / C1_bound_K
    p(f"  {g_label}:")
    p(f"    C1 = {g_val}^2 / (2 * {m_S}^2) = {C1_tree:.6e} MeV^-2")
    p(f"    C1 = {C1_tree*1e-6:.6e} GeV^-2")
    p(f"    UTfit bound: {C1_bound_K:.2e} MeV^-2 = {C1_bound_K*1e-6:.2e} GeV^-2")
    p(f"    Ratio C1/C1_bound = {ratio_utfit:.4e}")
    p(f"    {'EXCLUDED' if ratio_utfit > 1 else 'ALLOWED'}")
    p()

# Required coupling for m_S = f_pi to satisfy the bound:
g_max_K = np.sqrt(2.0 * m_S**2 * C1_bound_K)
p(f"Maximum coupling for m_S = {m_S} MeV to satisfy K-Kbar bound:")
p(f"  g_max = sqrt(2 * m_S^2 * C1_bound) = {g_max_K:.6e}")
p(f"  This is extremely small: {g_max_K:.2e}")
p()

# ═══════════════════════════════════════════════════════════════════════════
# PART (d): D-Dbar and B-Bbar constraints
# ═══════════════════════════════════════════════════════════════════════════

p("=" * 76)
p("PART (d): D-Dbar and B-Bbar constraints")
p("=" * 76)
p()

systems = [
    ("K-Kbar (Delta S=2)",  "d-s", m_K, f_K, B_K, Dm_K_exp,  C1_bound_K,  "Re(C1)"),
    ("D-Dbar (Delta C=2)",  "u-c", m_D, f_D, 0.76, Dm_D_exp,  C1_bound_D,  "|C1|"),
    ("B_d-Bbar (Delta B=2)","d-b", m_Bd, f_Bd, 1.26, Dm_Bd_exp, C1_bound_Bd, "|C1|"),
    ("B_s-Bbar (Delta B=2)","s-b", m_Bs, f_Bs, 1.26, Dm_Bs_exp, C1_bound_Bs, "|C1|"),
]

p(f"{'System':<24} {'Bound on C1':<18} {'C1(g=1,m=92)':<18} {'Ratio':<12} {'Status'}")
p("-" * 90)

for sys_name, flavs, m_meson, f_meson, B_bag, Dm_exp, C1_bound, bound_type in systems:
    C1_g1 = g_pert**2 / (2.0 * f_pi**2)
    ratio = C1_g1 / C1_bound
    status = "EXCLUDED" if ratio > 1 else "ALLOWED"
    p(f"  {sys_name:<24} {C1_bound:<18.4e} {C1_g1:<18.6e} {ratio:<12.2e} {status}")

p()

# Detailed analysis for each system
for sys_name, flavs, m_meson, f_meson, B_bag, Dm_exp, C1_bound, bound_type in systems:
    p(f"--- {sys_name} ---")
    p(f"  Meson mass: {m_meson:.2f} MeV")
    p(f"  Decay constant: {f_meson:.1f} MeV")
    p(f"  Bag parameter: {B_bag}")
    p(f"  Experimental Dm: {Dm_exp:.4e} MeV")
    p(f"  UTfit bound {bound_type}: {C1_bound:.4e} MeV^-2 = {C1_bound*1e-6:.4e} GeV^-2")
    p()

    # Hadronic matrix element
    H_meson = (8.0/3.0) * f_meson**2 * m_meson * B_bag

    # New physics contribution
    for g_label, g_val in [("g=1", 1.0), ("g=NDA", g_NDA)]:
        C1_calc = g_val**2 / (2.0 * f_pi**2)
        Dm_NP = C1_calc * H_meson
        ratio_Dm = Dm_NP / Dm_exp if Dm_exp > 0 else float('inf')
        ratio_C1 = C1_calc / C1_bound
        p(f"  {g_label}: C1 = {C1_calc:.4e}, Dm^NP = {Dm_NP:.4e} MeV, "
          f"Dm^NP/Dm^exp = {ratio_Dm:.2e}, C1/bound = {ratio_C1:.2e}")

    # Required scalar mass for g = 1:
    m_S_req = np.sqrt(g_pert**2 / (2.0 * C1_bound))
    p(f"  Required m_S for g=1 to satisfy bound: {m_S_req:.2f} MeV = {m_S_req/1e3:.4f} GeV")

    # Required scalar mass for g = NDA:
    m_S_req_NDA = np.sqrt(g_NDA**2 / (2.0 * C1_bound))
    p(f"  Required m_S for g=NDA to satisfy bound: {m_S_req_NDA:.2f} MeV = {m_S_req_NDA/1e3:.2f} GeV")
    p()

# ═══════════════════════════════════════════════════════════════════════════
# PART (e): Critical question — M_FCNC scale
# ═══════════════════════════════════════════════════════════════════════════

p("=" * 76)
p("PART (e): Required FCNC scale M_FCNC")
p("=" * 76)
p()

p("For tree-level scalar exchange, the Wilson coefficient is:")
p("  C1 = g^2 / (2 * M_FCNC^2)")
p()
p("Solving for M_FCNC from the bound C1 < C1_bound:")
p("  M_FCNC > g / sqrt(2 * C1_bound)")
p()

p(f"{'System':<24} {'C1_bound (MeV^-2)':<20} {'M_FCNC(g=1) (MeV)':<22} {'M_FCNC(g=1) (GeV)':<20} {'M_FCNC(g=NDA) (GeV)':<22}")
p("-" * 108)

M_FCNC_results = {}
for sys_name, flavs, m_meson, f_meson, B_bag, Dm_exp, C1_bound, bound_type in systems:
    M_FCNC_g1 = np.sqrt(g_pert**2 / (2.0 * C1_bound))
    M_FCNC_NDA = np.sqrt(g_NDA**2 / (2.0 * C1_bound))
    M_FCNC_results[sys_name] = (M_FCNC_g1, M_FCNC_NDA)
    p(f"  {sys_name:<24} {C1_bound:<20.4e} {M_FCNC_g1:<22.2f} {M_FCNC_g1/1e3:<20.4f} {M_FCNC_NDA/1e3:<22.2f}")

p()

# The model has m_scalar = f_pi = 92 MeV at the soft mass level.
# Compare with required M_FCNC:
p("Comparison with model's scalar mass m_S = f_pi = 92 MeV:")
p()

for sys_name, (M_g1, M_NDA) in M_FCNC_results.items():
    ratio_g1 = M_g1 / f_pi
    ratio_NDA = M_NDA / f_pi
    p(f"  {sys_name}:")
    p(f"    M_FCNC(g=1)/f_pi   = {ratio_g1:.2f}  (need scalar {ratio_g1:.0f}x heavier)")
    p(f"    M_FCNC(g=NDA)/f_pi = {ratio_NDA:.2f}  (need scalar {ratio_NDA:.0f}x heavier)")

p()
p("CONCLUSION: The K-Kbar constraint is by far the most stringent.")

M_K_g1, M_K_NDA = M_FCNC_results["K-Kbar (Delta S=2)"]
p(f"  For g=1:   M_FCNC > {M_K_g1:.0f} MeV = {M_K_g1/1e3:.1f} GeV")
p(f"  For g=NDA: M_FCNC > {M_K_NDA:.0f} MeV = {M_K_NDA/1e3:.0f} GeV")
p(f"  Model gives: m_S = {f_pi} MeV")
p()

# Tension factor
tension_g1 = M_K_g1 / f_pi
tension_NDA = M_K_NDA / f_pi
p(f"  Tension (g=1):   factor of {tension_g1:.1f}")
p(f"  Tension (g=NDA): factor of {tension_NDA:.0f}")
p()

p("Is m_tilde^2 = f_pi^2 compatible?")
p(f"  NO for g=1: need m_S > {M_K_g1:.0f} MeV, have {f_pi} MeV.")
p(f"  NO for g=NDA: need m_S > {M_K_NDA:.0f} MeV, have {f_pi} MeV.")
p()

# What soft mass is needed?
p("Required soft mass m_tilde for the off-diagonal mesons:")
p(f"  For g=1:   m_tilde > {M_K_g1:.0f} MeV (soft mass^2 > {M_K_g1**2:.2e} MeV^2)")
p(f"  For g=NDA: m_tilde > {M_K_NDA:.0f} MeV (soft mass^2 > {M_K_NDA**2:.2e} MeV^2)")
p()

# Alternative: what coupling is allowed at f_pi?
g_allowed_K = np.sqrt(2.0 * f_pi**2 * C1_bound_K)
p(f"Alternatively, if m_S = f_pi = {f_pi} MeV is fixed:")
p(f"  Maximum coupling from K-Kbar: g < {g_allowed_K:.6e}")
p(f"  This is far below both g=1 and g=NDA.")
p()

# ═══════════════════════════════════════════════════════════════════════════
# PART (f): Mass from tachyonic stabilization
# ═══════════════════════════════════════════════════════════════════════════

p("=" * 76)
p("PART (f): Off-diagonal scalar mass at the true vacuum (tachyonic stabilization)")
p("=" * 76)
p()

# At the M=0 vacuum (true global minimum from offdiag_vacuum analysis):
# The Hessian gives d^2V/d(M^i_j)^2 = 2*m_tilde^2 for ALL meson directions.
# This means the off-diagonal meson scalar mass^2 = 2*m_tilde^2 = 2*f_pi^2.

p("At the M=0 vacuum (true global minimum):")
p(f"  Scalar mass^2 = 2*m_tilde^2 = 2*f_pi^2 = {2*m_tilde_sq:.0f} MeV^2")
p(f"  Scalar mass = sqrt(2)*f_pi = {np.sqrt(2)*f_pi:.2f} MeV")
p()

# At the seesaw vacuum (metastable):
# The off-diagonal scalar masses get contributions from:
# 1. Soft mass: m_tilde^2 = f_pi^2
# 2. SUSY F-term mass squared (from second derivative of |F|^2)
# For the (M^d_s, M^s_d) pair at the seesaw vacuum:
#   The 2x2 SUSY mass matrix has eigenvalues ±|X*M_u|
#   (from W_{M^d_s, M^s_d} = X*M_u = X0*M_uu)
#   The scalar mass^2 matrix (including soft terms) is:
#   m_scalar^2 = |W_{ij}|^2 + m_tilde^2
#              = (X0*M_uu)^2 + f_pi^2
#   But at the seesaw vacuum, X0*M_uu is extremely small:

SUSY_mass_ds = abs(X0 * M_uu)
scalar_mass_sq_seesaw = SUSY_mass_ds**2 + m_tilde_sq

p("At the seesaw vacuum (metastable):")
p(f"  SUSY contribution: |W_(M^d_s, M^s_d)|^2 = |X*M_u|^2 = ({SUSY_mass_ds:.4e})^2 = {SUSY_mass_ds**2:.4e} MeV^2")
p(f"  Soft contribution: m_tilde^2 = {m_tilde_sq:.0f} MeV^2")
p(f"  Total: {scalar_mass_sq_seesaw:.4f} MeV^2")
p(f"  Scalar mass = {np.sqrt(scalar_mass_sq_seesaw):.4f} MeV")
p(f"  Note: SUSY contribution is negligible (ratio: {SUSY_mass_ds**2/m_tilde_sq:.2e})")
p()

# How heavy must the off-diagonal scalars be?
p("Required off-diagonal scalar masses from FCNC constraints:")
p()

for sys_name, (M_g1, M_NDA) in M_FCNC_results.items():
    p(f"  {sys_name}:")
    p(f"    g=1:   m_scalar > {M_g1:.0f} MeV")
    p(f"    g=NDA: m_scalar > {M_NDA:.0f} MeV")

p()

# Can the off-diagonal scalars get a large mass from the tachyonic stabilization?
# At the M=0 vacuum, the mass is sqrt(2)*f_pi ~ 130 MeV.
# This is NOT enough for any scenario.
#
# The key question is whether there exists a vacuum intermediate between
# M=0 and the seesaw, where off-diagonal scalars have large masses.
# The answer from the offdiag_vacuum analysis is NO: the global minimum
# is at M=0.
#
# Possible mechanisms to raise the off-diagonal scalar masses:
# 1. Non-universal soft masses: m_tilde^2(offdiag) >> m_tilde^2(diag)
# 2. D-term contributions from gauging the flavor symmetry
# 3. Additional superpotential terms beyond the minimal model
# 4. Running effects (RG enhancement above Lambda)

p("Possible mechanisms to raise off-diagonal scalar masses:")
p()
p("  1. Non-universal soft masses:")
p(f"     Need m_tilde^2(offdiag) > ({M_K_g1:.0f} MeV)^2 = {M_K_g1**2:.2e} MeV^2")
p(f"     Current: m_tilde^2 = f_pi^2 = {m_tilde_sq:.0f} MeV^2")
p(f"     Enhancement factor: {M_K_g1**2/m_tilde_sq:.0f} (for g=1)")
p()

p("  2. Alignment / decoupling:")
p("     If the off-diagonal mesons are decoupled from the low-energy")
p("     spectrum (e.g., they are integrated out at scale Lambda),")
p("     then M_FCNC ~ Lambda = 300 MeV.")
ratio_Lambda_g1 = M_K_g1 / Lambda_QCD
ratio_Lambda_NDA = M_K_NDA / Lambda_QCD
p(f"     Lambda = {Lambda_QCD} MeV vs M_FCNC(g=1) = {M_K_g1:.0f} MeV")
p(f"     Ratio: {ratio_Lambda_g1:.1f} (g=1), {ratio_Lambda_NDA:.0f} (g=NDA)")
p(f"     {'MARGINAL for g=1' if ratio_Lambda_g1 < 10 else 'INSUFFICIENT even at Lambda for g=1'}")
p()

p("  3. Suppressed coupling via alignment:")
p("     In supersymmetric models with aligned squark masses,")
p("     the FCNC coupling is suppressed by mixing angles.")
p("     If the meson-quark coupling carries a CKM-like suppression:")
p(f"     g_eff ~ V_ds ~ sin(theta_C) = 0.225")
g_Cabibbo = 0.225
M_FCNC_Cab = np.sqrt(g_Cabibbo**2 / (2.0 * C1_bound_K))
p(f"     M_FCNC(g=V_ds) = {M_FCNC_Cab:.0f} MeV = {M_FCNC_Cab/1e3:.2f} GeV")
p(f"     This helps but still requires m_scalar > {M_FCNC_Cab:.0f} MeV >> f_pi")
p()

p("  4. GIM cancellation:")
p("     In the SM, K-Kbar mixing is GIM-suppressed by (m_c^2 - m_u^2)/M_W^2.")
p("     In the meson model, an analogous mechanism could arise if")
p("     the off-diagonal meson exchange conserves flavour in the")
p("     massless limit. This requires the coupling matrix to be")
p("     proportional to the mass matrix (alignment).")
p()

# ═══════════════════════════════════════════════════════════════════════════
# SUMMARY
# ═══════════════════════════════════════════════════════════════════════════

p("=" * 76)
p("SUMMARY")
p("=" * 76)
p()

p("1. TREE-LEVEL FCNC from off-diagonal meson scalar exchange:")
p(f"   - Scalar mass at M=0 vacuum: sqrt(2)*f_pi = {m_scalar_M0:.1f} MeV")
p(f"   - Scalar mass at seesaw vacuum: ~f_pi = {f_pi} MeV")
p(f"   - Wilson coefficient C1(g=1, m=92): {g_pert**2/(2*f_pi**2):.4e} MeV^-2")
p(f"   - UTfit K-Kbar bound: {C1_bound_K:.2e} MeV^-2")
p(f"   - VIOLATION FACTOR: {g_pert**2/(2*f_pi**2)/C1_bound_K:.0e}")
p()

p("2. MESINO BOX DIAGRAMS:")
p(f"   - Mesino mass at seesaw vacuum: {m_mesino_ds:.2e} MeV (tiny)")
p(f"   - Box gives even larger FCNC contribution")
p(f"   - At M=0 vacuum: mesinos are massless (worse)")
p()

p("3. REQUIRED M_FCNC SCALE (from K-Kbar, most stringent):")
p(f"   - For g=1:   M_FCNC > {M_K_g1:.0f} MeV = {M_K_g1/1e3:.1f} GeV")
p(f"   - For g=NDA: M_FCNC > {M_K_NDA:.0f} MeV = {M_K_NDA/1e3:.0f} GeV")
p(f"   - Model value: f_pi = {f_pi} MeV")
p(f"   - TENSION: factor of {tension_g1:.0f} (g=1) to {tension_NDA:.0f} (g=NDA)")
p()

p("4. COMPATIBILITY of m_tilde^2 = f_pi^2:")
p(f"   INCOMPATIBLE with K-Kbar constraints by {tension_g1:.0f}x to {tension_NDA:.0f}x.")
p(f"   The off-diagonal mesons MUST acquire masses >> f_pi from some")
p(f"   mechanism not included in the minimal model.")
p()

p("5. POSSIBLE RESOLUTIONS:")
p("   (a) Non-universal soft masses: m_tilde^2(offdiag) >> f_pi^2")
p(f"       Need enhancement: ~{M_K_g1**2/m_tilde_sq:.0f}x for g=1, ~{M_K_NDA**2/m_tilde_sq:.0f}x for g=NDA")
p("   (b) Off-diagonal mesons decouple at Lambda (integrated out):")
p(f"       Then M_FCNC ~ Lambda = {Lambda_QCD} MeV, still marginal for g=1")
p("   (c) Flavour alignment: coupling g_eff ~ V_ds * m_q/Lambda")
p(f"       Reduces M_FCNC to {M_FCNC_Cab:.0f} MeV with Cabibbo suppression")
p("   (d) GIM-like mechanism from the meson matrix structure")
p("   (e) Off-diagonal scalars confined above Lambda (not in the")
p("       low-energy spectrum at all)")
p()

p("6. KEY NUMBERS TABLE:")
p()
p(f"   f_pi               = {f_pi:.0f} MeV")
p(f"   m_scalar(model)     = {f_pi:.0f} - {m_scalar_M0:.0f} MeV")
p(f"   M_FCNC (K, g=1)    = {M_K_g1:.0f} MeV = {M_K_g1/1e3:.1f} GeV")
p(f"   M_FCNC (K, g=NDA)  = {M_K_NDA:.0f} MeV = {M_K_NDA/1e3:.0f} GeV")

M_D_g1, M_D_NDA = M_FCNC_results["D-Dbar (Delta C=2)"]
M_Bd_g1, M_Bd_NDA = M_FCNC_results["B_d-Bbar (Delta B=2)"]
M_Bs_g1, M_Bs_NDA = M_FCNC_results["B_s-Bbar (Delta B=2)"]

p(f"   M_FCNC (D, g=1)    = {M_D_g1:.0f} MeV = {M_D_g1/1e3:.1f} GeV")
p(f"   M_FCNC (Bd, g=1)   = {M_Bd_g1:.0f} MeV = {M_Bd_g1/1e3:.2f} GeV")
p(f"   M_FCNC (Bs, g=1)   = {M_Bs_g1:.0f} MeV = {M_Bs_g1/1e3:.4f} GeV")
p()

p("Done.")

# ═══════════════════════════════════════════════════════════════════════════
# Write markdown report
# ═══════════════════════════════════════════════════════════════════════════

M("# FCNC Constraints on Off-Diagonal Meson Scalars in Confined SQCD")
M()
M("## Model Setup")
M()
M("N=1 SUSY with confined SU(3) SQCD (N_f = N_c = 3). The low-energy spectrum")
M("contains 9 complex meson scalars M^i_j and their fermionic superpartners")
M("(mesinos). The meson scalars have universal soft mass m_tilde^2 = f_pi^2 = (92 MeV)^2.")
M()
M("### Diagonal meson VEVs (Seiberg seesaw)")
M()
M(f"| Field | VEV (MeV) |")
M(f"|-------|-----------|")
M(f"| M^u_u | {M_uu:.0f} |")
M(f"| M^d_d | {M_dd:.0f} |")
M(f"| M^s_s | {M_ss:.0f} |")
M(f"| X | {X0:.4e} MeV^{{-4}} |")
M()
M("### Mesino masses at the seesaw vacuum")
M()
M(f"| Pair | Mass (MeV) | Expression |")
M(f"|------|-----------|------------|")
M(f"| d-s  | {m_mesino_ds:.4e} | |X| M_u |")
M(f"| u-s  | {m_mesino_us:.4e} | |X| M_d |")
M(f"| u-d  | {m_mesino_ud:.4e} | |X| M_s |")
M()
M("The mesino masses are extremely small (O(10^-4) MeV), far below the kaon mass.")
M()

M("---")
M()
M("## (a) Tree-Level K-Kbar Mixing")
M()
M("The off-diagonal meson M^d_s mediates Delta S = 2 transitions at tree level.")
M("The Wilson coefficient for the (LL)(LL) operator is:")
M()
M("C_1 = g^2 / (2 m_S^2)")
M()
M("where g is the effective d-s-scalar coupling and m_S is the scalar mass.")
M()
M("### Results for m_S = f_pi = 92 MeV")
M()
M("| Coupling | C_1 (GeV^-2) | UTfit bound (GeV^-2) | Ratio | Status |")
M("|----------|-------------|---------------------|-------|--------|")
for g_label, g_val in [("g = 1", 1.0), ("g = NDA (7.3)", g_NDA), ("g = V_ds (0.225)", g_Cabibbo)]:
    C1_val = g_val**2 / (2.0 * f_pi**2)
    C1_gev = C1_val * 1e-6
    ratio = C1_val / C1_bound_K
    status = "EXCLUDED" if ratio > 1 else "allowed"
    M(f"| {g_label} | {C1_gev:.4e} | 9.0e-13 | {ratio:.2e} | {status} |")
M()
M(f"The K-Kbar mass difference contribution Dm_K = C_1 * (8/3) f_K^2 m_K B_K exceeds")
M(f"the experimental value by a factor of {g_pert**2/(2*f_pi**2) * hadronic_K / Dm_K_exp:.0e} for g = 1.")
M()

M("---")
M()
M("## (b) Mesino Box Diagrams")
M()
M("Box diagrams with mesino propagators give:")
M()
M("C_1^(box) = g^4 / (16 pi^2 m_mesino^2)")
M()
M(f"With m_mesino(d-s) = {m_mesino_ds:.4e} MeV, these contributions are")
M(f"even more severely excluded than tree-level exchange:")
M()
C1_box_g1 = g_pert**4 * box_factor / m_mesino_ds**2
M(f"- g = 1: C_1^(box) = {C1_box_g1:.4e} MeV^-2, ratio to bound = {C1_box_g1/C1_bound_K:.2e}")
C1_box_NDA = g_NDA**4 * box_factor / m_mesino_ds**2
M(f"- g = NDA: C_1^(box) = {C1_box_NDA:.4e} MeV^-2, ratio to bound = {C1_box_NDA/C1_bound_K:.2e}")
M()
M("At the M = 0 vacuum, mesinos are massless (all F-term couplings vanish),")
M("but the effective coupling to quarks also vanishes.")
M()

M("---")
M()
M("## (c) Wilson Coefficient Comparison")
M()
M("### Effective 4-quark coupling at the kaon scale")
M()
M(f"For off-diagonal mesons with m_S^2 = f_pi^2 = {m_tilde_sq:.0f} MeV^2:")
M()
M(f"- C_1(g=1) = 1/(2 * {f_pi}^2) = {1.0/(2*f_pi**2):.6e} MeV^-2 = {1.0/(2*f_pi**2)*1e-6:.6e} GeV^-2")
M(f"- UTfit bound: |Re C_1| < 9.0e-13 GeV^-2")
M(f"- **Ratio: {1.0/(2*f_pi**2) / C1_bound_K:.0e}**")
M()
M(f"Maximum coupling allowed at m_S = {f_pi} MeV: g < {g_allowed_K:.2e}")
M()

M("---")
M()
M("## (d) D-Dbar and B-Bbar Constraints")
M()
M("| System | M_FCNC (g=1) | M_FCNC (g=NDA) | Model m_S | Tension factor |")
M("|--------|-------------|---------------|-----------|---------------|")
for sys_name, (M_g1, M_NDA) in M_FCNC_results.items():
    t_g1 = M_g1 / f_pi
    M(f"| {sys_name} | {M_g1:.0f} MeV ({M_g1/1e3:.2f} GeV) | {M_NDA:.0f} MeV ({M_NDA/1e3:.2f} GeV) | {f_pi:.0f} MeV | {t_g1:.1f}x (g=1) |")
M()
M("K-Kbar is the most stringent constraint by far.")
M()

M("---")
M()
M("## (e) Required M_FCNC Scale")
M()
M("The critical question: is m_tilde^2 = f_pi^2 compatible with FCNC bounds?")
M()
M("**No.** The K-Kbar constraint requires:")
M()
M(f"- For g = 1: M_FCNC > {M_K_g1/1e3:.1f} GeV (f_pi = 0.092 GeV)")
M(f"- For g = NDA: M_FCNC > {M_K_NDA/1e3:.0f} GeV")
M()
M(f"The off-diagonal meson scalars at m_S = f_pi are **{tension_g1:.0f}x too light** (g = 1)")
M(f"or **{tension_NDA:.0f}x too light** (g = NDA) to satisfy K-Kbar mixing constraints.")
M()
M(f"Required soft mass squared: m_tilde^2 > {M_K_g1**2:.2e} MeV^2 (g = 1), much larger than f_pi^2 = {m_tilde_sq:.0f} MeV^2.")
M()

M("---")
M()
M("## (f) Mass at the Stabilized Vacuum")
M()
M("At the M = 0 vacuum (true global minimum of the scalar potential):")
M(f"- Off-diagonal scalar mass^2 = 2 m_tilde^2 = 2 f_pi^2 = {2*m_tilde_sq:.0f} MeV^2")
M(f"- Scalar mass = sqrt(2) f_pi = {m_scalar_M0:.1f} MeV")
M()
M("At the seesaw vacuum (metastable):")
M(f"- SUSY contribution: |X M_u|^2 = {SUSY_mass_ds**2:.4e} MeV^2 (negligible)")
M(f"- Soft contribution: f_pi^2 = {m_tilde_sq:.0f} MeV^2 (dominant)")
M(f"- Total: m_scalar ~ {np.sqrt(scalar_mass_sq_seesaw):.1f} MeV")
M()
M("In **neither** vacuum do the off-diagonal scalars acquire masses large enough")
M("to satisfy FCNC bounds. The tachyonic stabilization at M = 0 gives")
M("m_scalar ~ 130 MeV, which is only marginally larger than f_pi.")
M()

M("---")
M()
M("## Possible Resolutions")
M()
M("1. **Non-universal soft masses**: If SUSY-breaking generates m_tilde^2(offdiag) >> f_pi^2")
M(f"   while keeping m_tilde^2(diag) ~ f_pi^2, the off-diagonal scalars can be pushed above")
M(f"   M_FCNC. Required enhancement: ~{M_K_g1**2/m_tilde_sq:.0f}x (for g = 1).")
M()
M("2. **Off-diagonal mesons decouple at Lambda**: If these states are confined above the")
M(f"   QCD scale and not present in the low-energy spectrum, M_FCNC ~ Lambda = {Lambda_QCD} MeV,")
M(f"   which is marginal for g = 1 (requires further suppression).")
M()
M("3. **Flavour alignment / GIM mechanism**: If the meson-quark coupling matrix is aligned")
M("   with the quark mass matrix, FCNC couplings are CKM-suppressed. With g_eff ~ V_ds ~ 0.225,")
M(f"   M_FCNC drops to {M_FCNC_Cab/1e3:.1f} GeV, still requiring scalars well above f_pi.")
M()
M("4. **Gauging flavour SU(3)**: D-term contributions from a gauged flavour symmetry could")
M("   generate large masses for off-diagonal scalars proportional to the D-term VEV.")
M()
M("5. **RG running**: Above the confinement scale, the Kahler potential receives large anomalous")
M("   dimension corrections that could enhance off-diagonal scalar masses.")
M()

M("---")
M()
M("## Key Numbers Summary")
M()
M("| Quantity | Value |")
M("|----------|-------|")
M(f"| f_pi | {f_pi} MeV |")
M(f"| m_tilde^2 = f_pi^2 | {m_tilde_sq:.0f} MeV^2 |")
M(f"| m_scalar (M = 0 vacuum) | {m_scalar_M0:.1f} MeV |")
M(f"| m_scalar (seesaw vacuum) | {np.sqrt(scalar_mass_sq_seesaw):.1f} MeV |")
M(f"| m_mesino(d-s) at seesaw | {m_mesino_ds:.4e} MeV |")
M(f"| M_FCNC (K-Kbar, g=1) | {M_K_g1:.0f} MeV = {M_K_g1/1e3:.1f} GeV |")
M(f"| M_FCNC (K-Kbar, g=NDA) | {M_K_NDA:.0f} MeV = {M_K_NDA/1e3:.0f} GeV |")
M(f"| M_FCNC (D-Dbar, g=1) | {M_D_g1:.0f} MeV = {M_D_g1/1e3:.1f} GeV |")
M(f"| M_FCNC (Bd-Bbar, g=1) | {M_Bd_g1:.0f} MeV = {M_Bd_g1/1e3:.2f} GeV |")
M(f"| Tension factor (K, g=1) | {tension_g1:.0f}x |")
M(f"| Tension factor (K, g=NDA) | {tension_NDA:.0f}x |")
M(f"| Max coupling at m_S = f_pi | {g_allowed_K:.2e} |")
M()
M("---")
M()
M("*Generated by fcnc_constraints.py*")

md_text = "\n".join(md_lines) + "\n"
md_path = "/home/codexssh/phys3/results/fcnc_constraints.md"
with open(md_path, "w") as f:
    f.write(md_text)

print(f"\n[Markdown report written to {md_path}]")
