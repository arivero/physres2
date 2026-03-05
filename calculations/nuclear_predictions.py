"""
Nuclear-scale predictions from the sBootstrap Lagrangian.

Compute:
1. f_pi^{SUSY} = sqrt(m_pi^2 - m_mu^2) vs f_pi^{PDG}
2. Goldberger-Treiman: g_A prediction with both f_pi values
3. g_piNN from GT relation
4. One-pion-exchange nuclear potential parameters
5. Sigma meson (diagonal mode) at f_pi = 92 MeV: range comparison
6. Deuteron binding: effective range theory estimate
7. Muon mass prediction from SUSY: m_mu = sqrt(m_pi^2 - f_pi^2)
8. Muonic atom parameters
"""

import math

# === Constants (PDG 2024) ===
m_pi_pm = 139.57039    # MeV, charged pion
m_pi_0  = 134.9768     # MeV, neutral pion
m_mu    = 105.6584     # MeV, muon
m_e     = 0.51100      # MeV, electron
m_N     = 938.272       # MeV, nucleon (average p,n)
m_p     = 938.272       # MeV, proton
m_n     = 939.565       # MeV, neutron
f_pi    = 92.07         # MeV, pion decay constant (particle physics convention)
g_A     = 1.2756        # axial coupling (neutron beta decay)
g_piNN  = 13.17         # pion-nucleon coupling (Nijmegen)
alpha   = 1.0/137.036   # fine structure constant
hbarc   = 197.3269804   # MeV fm

# Quark masses MSbar (PDG 2024)
m_u = 2.16   # MeV
m_d = 4.67   # MeV
m_s = 93.4   # MeV
Lambda_QCD = 300.0  # MeV (approximate)

print("=" * 75)
print("1. SUSY-PREDICTED f_pi")
print("=" * 75)

f_pi_susy = math.sqrt(m_pi_pm**2 - m_mu**2)
print(f"  m_pi± = {m_pi_pm} MeV")
print(f"  m_mu  = {m_mu} MeV")
print(f"  f_pi^SUSY = sqrt(m_pi^2 - m_mu^2) = {f_pi_susy:.4f} MeV")
print(f"  f_pi^PDG  = {f_pi:.2f} MeV")
print(f"  Ratio: f_pi^SUSY / f_pi^PDG = {f_pi_susy/f_pi:.6f}")
print(f"  Discrepancy: {abs(f_pi_susy - f_pi)/f_pi * 100:.2f}%")

ratio_sq = (m_pi_pm**2 - m_mu**2) / f_pi**2
print(f"\n  (m_pi^2 - m_mu^2)/f_pi^2 = {ratio_sq:.4f}")
print(f"  Deviation from 1: {abs(ratio_sq - 1)*100:.2f}%")

print("\n" + "=" * 75)
print("2. GOLDBERGER-TREIMAN RELATION")
print("=" * 75)

# GT: g_piNN * f_pi = g_A * m_N
# => g_A^GT = g_piNN * f_pi / m_N
g_A_pdg = g_piNN * f_pi / m_N
g_A_susy = g_piNN * f_pi_susy / m_N

print(f"  g_A^exp = {g_A}")
print(f"  g_piNN = {g_piNN}")
print(f"  m_N = {m_N} MeV")
print(f"\n  With f_pi^PDG = {f_pi} MeV:")
print(f"    g_A^GT = g_piNN * f_pi / m_N = {g_A_pdg:.4f}")
print(f"    GT discrepancy: |g_A^GT - g_A^exp|/g_A^exp = {abs(g_A_pdg - g_A)/g_A*100:.3f}%")
print(f"    ({abs(g_A_pdg - g_A)/g_A*100:.1f}% or {(g_A_pdg - g_A)/0.0013:.1f}sigma)")

print(f"\n  With f_pi^SUSY = {f_pi_susy:.2f} MeV:")
print(f"    g_A^GT = g_piNN * f_pi / m_N = {g_A_susy:.4f}")
print(f"    GT discrepancy: |g_A^GT - g_A^exp|/g_A^exp = {abs(g_A_susy - g_A)/g_A*100:.3f}%")
print(f"    ({abs(g_A_susy - g_A)/g_A*100:.1f}% or {(g_A_susy - g_A)/0.0013:.1f}sigma)")

ratio_improve = abs(g_A_pdg - g_A) / abs(g_A_susy - g_A)
print(f"\n  ==> SUSY f_pi reduces GT discrepancy by factor {ratio_improve:.1f}")

print("\n  Alternatively: predict g_piNN from (g_A, f_pi, m_N):")
g_piNN_pdg = g_A * m_N / f_pi
g_piNN_susy = g_A * m_N / f_pi_susy
print(f"    g_piNN from f_pi^PDG:  {g_piNN_pdg:.3f}  (exp: {g_piNN})")
print(f"    g_piNN from f_pi^SUSY: {g_piNN_susy:.3f}  (exp: {g_piNN})")

print("\n" + "=" * 75)
print("3. ONE-PION-EXCHANGE NUCLEAR POTENTIAL")
print("=" * 75)

# Range of OPEP
r_pi = hbarc / m_pi_pm
r_pi_susy = hbarc / (math.sqrt(2) * f_pi_susy)  # sBootstrap: m_pi = sqrt(2) f_pi
r_pi_pdg = hbarc / m_pi_pm

print(f"  Pion Compton wavelength:")
print(f"    1/m_pi (PDG) = {r_pi_pdg:.4f} fm")
print(f"    1/m_pi (sBootstrap: sqrt(2)*f_pi) = {hbarc/(math.sqrt(2)*f_pi_susy):.4f} fm")
print(f"    1/m_pi (sBootstrap: sqrt(2)*f_pi^PDG) = {hbarc/(math.sqrt(2)*f_pi):.4f} fm")

# Sigma meson (diagonal mode) at f_pi
m_sigma_sb = f_pi_susy  # sBootstrap diagonal mode
r_sigma_sb = hbarc / m_sigma_sb
m_sigma_phys = 450.0  # approximate f_0(500)
r_sigma_phys = hbarc / m_sigma_phys

print(f"\n  Sigma meson range:")
print(f"    sBootstrap: m_sigma = f_pi = {m_sigma_sb:.1f} MeV => range = {r_sigma_sb:.3f} fm")
print(f"    Physical f_0(500): m ~ {m_sigma_phys} MeV => range = {r_sigma_phys:.3f} fm")
print(f"    Ratio: sBootstrap sigma range / pion range = {r_sigma_sb / r_pi_pdg:.2f}")

# OPEP strength at 1 fm
V_opep_1fm = -(g_piNN**2 / (4 * math.pi)) * (m_pi_pm / (3 * m_N)) * math.exp(-m_pi_pm / hbarc) * hbarc
# Actually: V_central = -(g^2/4pi)(m_pi^3/(12 m_N^2))(e^{-x}/x) where x = m_pi r
# At r = 1 fm: x = m_pi * 1 fm / hbarc = 139.57 * 1 / 197.33 = 0.707
x_1fm = m_pi_pm * 1.0 / hbarc
V_central_1fm = -(g_piNN**2 / (4*math.pi)) * (m_pi_pm**3 / (12 * m_N**2)) * math.exp(-x_1fm) / x_1fm
print(f"\n  OPEP central part at r = 1 fm:")
print(f"    x = m_pi*r/hbarc = {x_1fm:.4f}")
print(f"    V_C(1 fm) = {V_central_1fm:.2f} MeV")

# With sigma exchange at sBootstrap mass
x_sigma_1fm = m_sigma_sb * 1.0 / hbarc
V_sigma_1fm_generic = math.exp(-x_sigma_1fm) / x_sigma_1fm  # in units of coupling
print(f"\n  Sigma Yukawa factor at r = 1 fm:")
print(f"    sBootstrap (92 MeV): e^(-x)/x = {math.exp(-x_sigma_1fm)/x_sigma_1fm:.4f} (x = {x_sigma_1fm:.3f})")
print(f"    Physical (450 MeV):  e^(-x)/x = {math.exp(-m_sigma_phys/hbarc)/(m_sigma_phys/hbarc):.4f} (x = {m_sigma_phys/hbarc:.3f})")
print(f"    Ratio: {(math.exp(-x_sigma_1fm)/x_sigma_1fm)/(math.exp(-m_sigma_phys/hbarc)/(m_sigma_phys/hbarc)):.2f}")

print("\n" + "=" * 75)
print("4. DEUTERON BINDING ENERGY ESTIMATE")
print("=" * 75)

# Effective range theory: the deuteron is a weakly bound state
# E_d = -hbar^2 / (m_N * a_s^2) where a_s is the triplet scattering length
# Experimental: a_t = 5.424 fm, E_d = 2.224 MeV
a_t_exp = 5.424  # fm, triplet scattering length
r_t_exp = 1.759  # fm, effective range
E_d_exp = 2.2246 # MeV

# From effective range theory: 1/a + sqrt(m_N E_d) / hbar = r_eff * m_N * E_d / (2*hbar^2) + ...
kappa = math.sqrt(m_N * E_d_exp / 2) / hbarc  # in fm^{-1}
print(f"  Deuteron kappa = sqrt(m_N E_d / 2) / hbarc = {kappa:.4f} fm^-1")
print(f"  Deuteron size: 1/kappa = {1/kappa:.2f} fm")
print(f"  Experimental a_t = {a_t_exp} fm, r_t = {r_t_exp} fm")
print(f"  E_d = {E_d_exp} MeV")

# The pion range is r_pi = 1.414 fm. The deuteron size is 4.3 fm.
# The deuteron is held together by the pion exchange tail.
print(f"\n  Ratio deuteron size / pion range = {1/kappa / r_pi_pdg:.2f}")
print(f"  Ratio deuteron size / sigma range (sBootstrap) = {1/kappa / r_sigma_sb:.2f}")
print(f"  The deuteron lives in the PION EXCHANGE regime, not the sigma regime.")

# Simple square-well estimate: V_0 * R^2 = pi^2 hbar^2 / (4 m_red)
# For barely bound state: V_0 = pi^2 hbar^2 / (4 m_red R^2)
m_red = m_N / 2  # reduced mass
R_eff = r_pi_pdg  # range ~ pion wavelength
V_0_sq = math.pi**2 * hbarc**2 / (4 * m_red * R_eff**2)
print(f"\n  Square-well estimate (depth for just-bound state):")
print(f"    R = r_pi = {R_eff:.3f} fm")
print(f"    V_0 = pi^2 hbarc^2/(4 m_red R^2) = {V_0_sq:.1f} MeV")
print(f"    Typical nuclear potential depth: ~50 MeV")

print("\n" + "=" * 75)
print("5. MUON MASS FROM SUSY AND MUONIC ATOMS")
print("=" * 75)

m_mu_susy = math.sqrt(m_pi_pm**2 - f_pi**2)
print(f"  m_mu^SUSY = sqrt(m_pi^2 - f_pi^2) = sqrt({m_pi_pm**2:.1f} - {f_pi**2:.1f})")
print(f"           = sqrt({m_pi_pm**2 - f_pi**2:.1f}) = {m_mu_susy:.4f} MeV")
print(f"  m_mu^exp = {m_mu} MeV")
print(f"  Discrepancy: {abs(m_mu_susy - m_mu)/m_mu * 100:.3f}%")
print(f"  ({abs(m_mu_susy - m_mu):.3f} MeV, {abs(m_mu_susy - m_mu)/0.0000035:.0f} sigma in m_mu measurement)")

# Muonic hydrogen Bohr radius
a_mu = hbarc / (alpha * m_mu)  # fm
a_mu_susy = hbarc / (alpha * m_mu_susy)
print(f"\n  Muonic hydrogen Bohr radius:")
print(f"    a_mu (exp m_mu) = {a_mu:.1f} fm")
print(f"    a_mu (SUSY m_mu) = {a_mu_susy:.1f} fm")
print(f"    Shift: {abs(a_mu_susy - a_mu):.1f} fm ({abs(a_mu_susy - a_mu)/a_mu*100:.2f}%)")

# Muonic hydrogen binding energy
E_bind_mu = alpha**2 * m_mu / 2  # keV
E_bind_mu_susy = alpha**2 * m_mu_susy / 2
print(f"\n  Muonic hydrogen 1s binding energy:")
print(f"    E_1s (exp m_mu) = {E_bind_mu*1000:.3f} eV = {E_bind_mu:.5f} keV")
print(f"    E_1s (SUSY m_mu) = {E_bind_mu_susy*1000:.3f} eV = {E_bind_mu_susy:.5f} keV")
print(f"    Shift: {abs(E_bind_mu_susy - E_bind_mu)*1000:.1f} eV")

# Muon-catalyzed fusion: ddmu formation
print(f"\n  Muon-catalyzed fusion:")
print(f"    dtmu molecular ion binding ~ 0.66 eV (J=1 state)")
print(f"    The SUSY shift in E_1s ({abs(E_bind_mu_susy - E_bind_mu)*1000:.0f} eV) is much LARGER")
print(f"    than the molecular binding (0.66 eV)")
print(f"    => The 0.5% SUSY prediction for m_mu is not precise enough for uCF rates")
print(f"    => But it EXPLAINS m_mu as derived quantity, not free parameter")

print("\n" + "=" * 75)
print("6. PION-NUCLEON SIGMA TERM")
print("=" * 75)

sigma_piN_exp = 59.0  # MeV (lattice + dispersive)
sigma_piN_err = 7.0

# GOR: m_pi^2 f_pi^2 = -(m_u + m_d) <qbar q>
# sigma_piN = m_hat * d/dm_hat (m_N) where m_hat = (m_u + m_d)/2
# sigma_piN / f_pi^2 = m_pi^2 * (fraction of m_N from quark mass)
m_hat = (m_u + m_d) / 2  # MeV
print(f"  m_hat = (m_u + m_d)/2 = {m_hat:.3f} MeV")
print(f"  sigma_piN^exp = {sigma_piN_exp} +/- {sigma_piN_err} MeV")

# From SUSY: m_pi^2 = 2 f_pi^2 (off-diagonal soft mass)
# Then: -(m_u+m_d)<qbar q> = m_pi^2 f_pi^2 = 2 f_pi^4
# => |<qbar q>| = 2 f_pi^4 / (m_u + m_d)
qqbar = 2 * f_pi**4 / (m_u + m_d)
qqbar_cube = qqbar**(1/3)
print(f"\n  From SUSY m_pi^2 = 2 f_pi^2:")
print(f"    |<qbar q>| = 2 f_pi^4 / (m_u + m_d)")
print(f"              = 2 * {f_pi:.1f}^4 / {m_u+m_d:.2f}")
print(f"              = {qqbar:.0f} MeV^3")
print(f"    |<qbar q>|^(1/3) = {qqbar_cube:.1f} MeV")
print(f"    Standard: |<qbar q>|^(1/3) ~ 260 MeV")
print(f"    Discrepancy: {abs(qqbar_cube - 260)/260*100:.1f}%")

# sigma_piN ~ m_hat * m_N * <sigma_piN_content>
# Chiral perturbation: sigma_piN = m_hat * <N|uu+dd|N>
# Alternatively: sigma_piN ~ m_pi^2 * (dlog m_N / dlog m_pi^2) * m_N
# In sBootstrap: m_N ~ Lambda_QCD (mostly gluonic), so sigma_piN ~ m_hat/Lambda * m_N?
sigma_naive = m_hat * m_N / Lambda_QCD
print(f"\n  Naive estimate: sigma_piN ~ m_hat * m_N / Lambda")
print(f"    = {m_hat:.2f} * {m_N:.0f} / {Lambda_QCD:.0f}")
print(f"    = {sigma_naive:.1f} MeV (vs exp {sigma_piN_exp} MeV)")

print("\n" + "=" * 75)
print("7. SUMMARY TABLE")
print("=" * 75)

print("""
  Quantity                   sBootstrap           Experiment         Status
  ─────────────────────────────────────────────────────────────────────────
  f_pi^SUSY                  {:.2f} MeV          {:.2f} MeV       {:.1f}% off
  m_pi (off-diag, sqrt(2)*f) {:.1f} MeV           {:.2f} MeV       {:.1f}% off
  GT discrepancy (f_pi^PDG)  {:.2f}%              0%                standard
  GT discrepancy (f_pi^SUSY) {:.2f}%              0%                4x improved
  g_piNN (from GT+SUSY f_pi) {:.2f}               {:.2f}            {:.1f}% off
  m_mu (from SUSY relation)  {:.2f} MeV          {:.4f} MeV      {:.2f}% off
  m_sigma (diagonal mode)    {:.1f} MeV           ~450 MeV         different!
  |<qbar q>|^(1/3)           {:.1f} MeV           ~260 MeV         {:.1f}% off
""".format(
    f_pi_susy, f_pi, abs(f_pi_susy - f_pi)/f_pi*100,
    math.sqrt(2)*f_pi_susy, m_pi_pm, abs(math.sqrt(2)*f_pi_susy - m_pi_pm)/m_pi_pm*100,
    abs(g_A_pdg - g_A)/g_A*100,
    abs(g_A_susy - g_A)/g_A*100,
    g_piNN_susy, g_piNN, abs(g_piNN_susy - g_piNN)/g_piNN*100,
    m_mu_susy, m_mu, abs(m_mu_susy-m_mu)/m_mu*100,
    m_sigma_sb,
    qqbar_cube, abs(qqbar_cube-260)/260*100,
))
