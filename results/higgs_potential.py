"""
Two-Higgs-doublet scalar potential at tan(beta) = 1
=====================================================

N=1 SUSY 2HDM with meson singlet Yukawa couplings from an ISS seesaw.
Computes: EWSB conditions, Higgs mass matrices, tree-level spectrum,
          soft-mass implications from F-term SUSY breaking.
"""

import numpy as np

# =================================================================
# Numerical inputs
# =================================================================
v_ew    = 246.22e3       # MeV  (electroweak VEV)
m_Z     = 91.1876e3      # MeV
g       = 0.653          # SU(2)_L coupling
gp      = 0.350          # U(1)_Y coupling
m_c     = 1270.0         # MeV  (charm quark mass, MS-bar)
m_b     = 4180.0         # MeV  (bottom quark mass, MS-bar)
m_t     = 172.76e3       # MeV  (top quark pole mass)
m_h_obs = 125.25e3       # MeV  (observed Higgs mass)
C       = 882297.0       # MeV^2  (SQCD seesaw scale)

# ISS meson VEVs: <M_j> = C / m_j
M_c_vev = C / m_c
M_b_vev = C / m_b

# tan(beta) = 1  =>  v_u = v_d = v/sqrt(2)
tan_beta = 1.0
cos2beta = np.cos(2 * np.arctan(tan_beta))  # = 0
sin2beta = np.sin(2 * np.arctan(tan_beta))  # = 1
v_u = v_ew / np.sqrt(2)
v_d = v_ew / np.sqrt(2)

# Yukawa couplings at tan(beta) = 1:
#   m_c = y_c * v_u / sqrt(2) = y_c * v / 2  =>  y_c = 2 m_c / v
y_c = 2 * m_c / v_ew
y_b = 2 * m_b / v_ew
y_t = np.sqrt(2) * m_t / v_ew   # top: m_t = y_t v sin(beta) / sqrt(2)

# Gauge couplings
gz2 = g**2 + gp**2
m_W_tree = g * v_ew / 2
m_Z_tree = np.sqrt(gz2) * v_ew / 2

# SUSY-breaking F-term scale
F_susy = C / v_ew

# =================================================================
# Print inputs
# =================================================================
print("=" * 70)
print("  TWO-HIGGS-DOUBLET SCALAR POTENTIAL AT tan(beta) = 1")
print("=" * 70)

print("\n  Inputs:")
print(f"    v       = {v_ew/1e3:.4f} GeV")
print(f"    m_Z     = {m_Z/1e3:.6f} GeV")
print(f"    g       = {g:.4f},  g' = {gp:.4f}")
print(f"    m_c     = {m_c:.1f} MeV,  m_b = {m_b:.1f} MeV,  m_t = {m_t/1e3:.2f} GeV")
print(f"    C       = {C:.0f} MeV^2")

print(f"\n  Derived:")
print(f"    v_u = v_d = v/sqrt(2) = {v_u/1e3:.4f} GeV")
print(f"    y_c = 2 m_c / v = {y_c:.8f}")
print(f"    y_b = 2 m_b / v = {y_b:.8f}")
print(f"    y_t = sqrt(2) m_t / v = {y_t:.6f}")
print(f"    <M_c> = C/m_c = {M_c_vev:.2f} MeV")
print(f"    <M_b> = C/m_b = {M_b_vev:.2f} MeV")
print(f"    cos(2 beta) = {cos2beta:.10f}")

# =================================================================
# 1. F-term structure
# =================================================================
print("\n" + "=" * 70)
print("  1. F-TERM STRUCTURE")
print("=" * 70)

# The superpotential W_H = mu_H H_u.H_d + y_c H_u M_c + y_b H_d M_b
# gives F-terms:
#   F_{H_u}^* = mu_H H_d^0 + y_c <M_c>
#   F_{H_d}^* = mu_H H_u^0 + y_b <M_b>
#   F_{M_c}^* = y_c H_u^0
#   F_{M_b}^* = y_b H_d^0

# The key ISS identity: y_j <M_j> = (2 m_j/v)(C/m_j) = 2C/v
yc_Mc = y_c * M_c_vev
yb_Mb = y_b * M_b_vev
two_C_over_v = 2 * C / v_ew

print(f"\n  y_c <M_c> = {yc_Mc:.6f} MeV")
print(f"  y_b <M_b> = {yb_Mb:.6f} MeV")
print(f"  2C/v       = {two_C_over_v:.6f} MeV")
print(f"  => y_c <M_c> = y_b <M_b> = 2C/v  (ISS seesaw identity)")

# Meson F-terms at the EWSB vacuum:
F_Mc = y_c * v_u / np.sqrt(2)  # = y_c v/2 = m_c
F_Mb = y_b * v_d / np.sqrt(2)  # = y_b v/2 = m_b

print(f"\n  F_{{M_c}} = y_c v_u/sqrt(2) = {F_Mc:.2f} MeV  (= m_c)")
print(f"  F_{{M_b}} = y_b v_d/sqrt(2) = {F_Mb:.2f} MeV  (= m_b)")
print(f"  The meson F-terms equal the quark masses: F_{{M_j}} = m_j.")

print(f"\n  SUSY-breaking scale: F = C/v = {F_susy:.4f} MeV = {F_susy/1e3:.6f} GeV")

# =================================================================
# 2. Scalar potential and EWSB
# =================================================================
print("\n" + "=" * 70)
print("  2. SCALAR POTENTIAL AND EWSB CONDITIONS")
print("=" * 70)

# The full scalar potential for neutral components (setting mu_H = 0):
#
# V = y_c^2 |H_u^0|^2 + y_b^2 |H_d^0|^2                  (F-terms from mesons)
#   + m_Hu^2 |H_u^0|^2 + m_Hd^2 |H_d^0|^2                (soft masses)
#   - b (H_u^0 H_d^0 + h.c.)                               (soft bilinear, b > 0)
#   + (g^2 + g'^2)/8 (|H_u^0|^2 - |H_d^0|^2)^2            (D-terms)
#
# Define effective masses:
#   m_1^2 = m_Hd^2 + y_b^2
#   m_2^2 = m_Hu^2 + y_c^2
#
# V = m_2^2 |H_u^0|^2 + m_1^2 |H_d^0|^2 - b(H_u^0 H_d^0 + h.c.) + V_D
#
# Minimization at v_u = v_d = v/sqrt(2) (where D-term vanishes):
#   m_2^2 v_u - b v_d = 0  =>  m_2^2 = b
#   m_1^2 v_d - b v_u = 0  =>  m_1^2 = b
#
# Both give:  m_1^2 = m_2^2 = b

print(f"""
  With mu_H = 0 and effective masses m_1^2 = m_Hd^2 + y_b^2,  m_2^2 = m_Hu^2 + y_c^2:

  Minimization at tan(beta) = 1 (D-term vanishes):
    m_1^2 = m_2^2 = b

  Constraint 1 (difference):
    m_Hu^2 - m_Hd^2 = y_b^2 - y_c^2 = {y_b**2 - y_c**2:.6e}

  Constraint 2 (each equation):
    m_Hu^2 = b - y_c^2
    m_Hd^2 = b - y_b^2

  The soft bilinear b remains as a free parameter.
  It determines the pseudoscalar mass: m_A^2 = m_1^2 + m_2^2 = 2b.""")

# =================================================================
# 3. CP-even Higgs mass matrix
# =================================================================
print("\n" + "=" * 70)
print("  3. CP-EVEN HIGGS MASS MATRIX")
print("=" * 70)

# Standard MSSM result in (h_d, h_u) basis:
#   M^2 = [[m_A^2 sin^2 beta + m_Z^2 cos^2 beta,  -(m_A^2 + m_Z^2) sin beta cos beta],
#          [-(m_A^2 + m_Z^2) sin beta cos beta,     m_A^2 cos^2 beta + m_Z^2 sin^2 beta]]
#
# At tan beta = 1: sin beta = cos beta = 1/sqrt(2), sin^2 = cos^2 = 1/2
#
# M^2 = [[(m_A^2 + m_Z^2)/2,  -(m_A^2 + m_Z^2)/2],
#        [-(m_A^2 + m_Z^2)/2,   (m_A^2 + m_Z^2)/2]]
#
# Eigenvalues:
#   m_h^2 = 0
#   m_H^2 = m_A^2 + m_Z^2
#
# Eigenvectors: h = (h_u + h_d)/sqrt(2)  (massless),
#               H = (h_u - h_d)/sqrt(2)  (massive)

print(f"""
  In the (h_d, h_u) basis, the MSSM CP-even mass matrix at tan(beta) = 1:

    M^2 = (m_A^2 + m_Z^2)/2 * [[1, -1], [-1, 1]]

  Eigenvalues:
    m_h^2 = 0                   (D-flat direction, massless at tree level)
    m_H^2 = m_A^2 + m_Z^2      (D-term direction)

  The light eigenstate h = (h_u + h_d)/sqrt(2) is the D-flat direction.
  The only quartic coupling is from D-terms, which vanishes along h_u = h_d.

  MSSM tree-level upper bound:
    m_h <= m_Z |cos(2 beta)| = {m_Z/1e3:.4f} * |{cos2beta:.4f}| = 0 GeV""")

# =================================================================
# 4. Full Higgs spectrum at tree level
# =================================================================
print("\n" + "=" * 70)
print("  4. TREE-LEVEL HIGGS SPECTRUM")
print("=" * 70)

print(f"""
  CP-even:   m_h = 0,           m_H = sqrt(m_A^2 + m_Z^2)
  CP-odd:    m_A = sqrt(2b),    G^0 = 0 (eaten by Z)
  Charged:   m_{{H^pm}} = sqrt(m_A^2 + m_W^2),  G^pm = 0 (eaten by W)

  m_W (tree) = g v/2 = {m_W_tree/1e3:.4f} GeV

  Spectrum for various m_A:
  {'m_A (GeV)':>12}  {'m_H (GeV)':>12}  {'m_H^pm (GeV)':>14}  {'b (GeV^2)':>12}  {'m_h':>6}""")

for m_A_GeV in [50, 100, 200, 300, 500, 1000]:
    m_A_MeV = m_A_GeV * 1e3
    m_H_GeV = np.sqrt(m_A_MeV**2 + m_Z**2) / 1e3
    m_Hpm_GeV = np.sqrt(m_A_MeV**2 + m_W_tree**2) / 1e3
    b_GeV2 = (m_A_MeV / 1e3)**2 / 2
    print(f"  {m_A_GeV:12.0f}  {m_H_GeV:12.2f}  {m_Hpm_GeV:14.2f}  {b_GeV2:12.0f}  {'0':>6}")

# =================================================================
# 5. Soft masses from F-terms
# =================================================================
print("\n" + "=" * 70)
print("  5. SOFT MASSES FROM F-TERM SUSY BREAKING")
print("=" * 70)

# The Yukawa F-terms contribute effective mass-squared terms:
#   delta m_Hu^2 = y_c^2     (from |F_{M_c}|^2 = y_c^2 |H_u|^2)
#   delta m_Hd^2 = y_b^2     (from |F_{M_b}|^2 = y_b^2 |H_d|^2)
#
# From the EWSB conditions:
#   m_Hu^2 = b - y_c^2
#   m_Hd^2 = b - y_b^2

print(f"""
  F-term contributions to scalar masses (absorbed into m_1^2, m_2^2):
    y_c^2 = {y_c**2:.6e}   (from |F_{{M_c}}|^2)
    y_b^2 = {y_b**2:.6e}   (from |F_{{M_b}}|^2)

  Soft masses from EWSB (in terms of free parameter b):
    m_Hu^2 = b - y_c^2 = b - {y_c**2:.6e}
    m_Hd^2 = b - y_b^2 = b - {y_b**2:.6e}
    m_Hu^2 - m_Hd^2 = y_b^2 - y_c^2 = {y_b**2 - y_c**2:.6e}

  The Higgs F-terms are flavor-universal:
    F_{{H_u}} = y_c <M_c> = 2C/v = {yc_Mc:.4f} MeV = {yc_Mc/1e3:.6f} GeV
    F_{{H_d}} = y_b <M_b> = 2C/v = {yb_Mb:.4f} MeV = {yb_Mb/1e3:.6f} GeV

  If b is generated by F-term mediation at a scale M_med:
    b ~ F_{{H_u}} F_{{H_d}} / M_med = (2C/v)^2 / M_med
    (2C/v)^2 = {(2*C/v_ew)**2:.4f} MeV^2 = {(2*C/v_ew)**2/1e6:.6e} GeV^2""")

# =================================================================
# 6. Radiative corrections needed for m_h = 125 GeV
# =================================================================
print("\n" + "=" * 70)
print("  6. RADIATIVE CORRECTIONS")
print("=" * 70)

# Leading-log one-loop stop correction at tan(beta) = 1:
#   Delta m_h^2 = (3 y_t^4 v^2 sin^4 beta) / (8 pi^2) * ln(m_stop^2/m_t^2)
#   sin^4(beta) = 1/4 at tan(beta) = 1
#   Delta m_h^2 = (3 y_t^4 v^2) / (32 pi^2) * ln(m_stop^2/m_t^2)
#
# To get Delta m_h^2 = m_h_obs^2:
log_needed = 32 * np.pi**2 * m_h_obs**2 / (3 * y_t**4 * v_ew**2)
m_stop = m_t * np.exp(log_needed / 2)

# NMSSM singlet: extra quartic lambda^2 |H_u.H_d|^2
#   Delta m_h^2 = lambda^2 v^2 sin^2(2 beta) / 2 = lambda^2 v^2 / 2  at sin(2beta)=1
lambda_NMSSM = np.sqrt(2) * m_h_obs / v_ew

print(f"""
  At tan(beta) = 1, the tree-level m_h = 0 requires large radiative corrections.

  (a) Stop loops (leading-log, MSSM):
      Delta m_h^2 = 3 y_t^4 v^2 sin^4(beta) / (8 pi^2) * ln(m_stop^2/m_t^2)
      sin^4(beta) = {np.sin(np.arctan(1))**4:.4f}  (= 1/4)
      For m_h = {m_h_obs/1e3:.2f} GeV:
        ln(m_stop^2/m_t^2) = {log_needed:.2f}
        m_stop = {m_stop/1e3:.0f} GeV
      This is unphysically large. Pure MSSM at tan(beta) = 1 is disfavored.

  (b) NMSSM singlet coupling (lambda S H_u H_d):
      Delta m_h^2 = lambda^2 v^2 / 2
      For m_h = {m_h_obs/1e3:.2f} GeV:
        lambda = sqrt(2) m_h / v = {lambda_NMSSM:.4f}
      This is perturbative (lambda < 1). The NMSSM resolves the tan(beta)=1 problem.""")

# =================================================================
# 7. ISS identity: flavor-universal Higgs F-term
# =================================================================
print("\n" + "=" * 70)
print("  7. ISS SEESAW: FLAVOR-UNIVERSAL F-TERM")
print("=" * 70)

print(f"""
  The ISS seesaw VEV <M_j> = C/m_j combined with y_j = 2m_j/v gives:

    y_j <M_j> = (2 m_j / v)(C / m_j) = 2C/v

  This cancellation is exact and flavor-independent. It means:
  - The F-term contribution to the Higgs potential from each flavor is identical
  - The Higgs sector does not distinguish charm from bottom
  - SUSY breaking in the Higgs sector is flavor-blind despite different Yukawas

  Numerical verification:
    y_c <M_c> = {y_c:.8f} * {M_c_vev:.4f} = {y_c * M_c_vev:.6f} MeV
    y_b <M_b> = {y_b:.8f} * {M_b_vev:.4f} = {y_b * M_b_vev:.6f} MeV
    2C/v       = 2 * {C:.0f} / {v_ew:.0f}       = {2*C/v_ew:.6f} MeV
    Match: {abs(y_c * M_c_vev - 2*C/v_ew) < 1e-6} (to machine precision)""")

# =================================================================
# Summary
# =================================================================
print("\n" + "=" * 70)
print("  SUMMARY")
print("=" * 70)
print(f"""
  At tan(beta) = 1 in the SUSY 2HDM with ISS meson Yukawas:

  1. The D-term quartic vanishes along the D-flat direction h_u = h_d.
  2. The meson Yukawa F-terms (y_c^2, y_b^2) contribute mass terms, not quartics.
  3. The lightest CP-even Higgs is massless at tree level: m_h = 0.
  4. The MSSM bound m_h <= m_Z|cos 2beta| = 0 is saturated.
  5. The ISS seesaw makes the Higgs F-terms flavor-universal: y_j<M_j> = 2C/v.
  6. The F-term of each meson equals the quark mass: F_{{M_j}} = m_j.
  7. Achieving m_h = 125 GeV requires either:
     - Absurdly heavy stops (m_stop ~ {m_stop/1e3:.0e} GeV), or
     - An NMSSM-like singlet with lambda = {lambda_NMSSM:.4f}.
""")

print("Done.")
