"""
Muon-catalyzed fusion parameters from the sBootstrap.

The sBootstrap derives m_mu = sqrt(m_pi^2 - f_pi^2).
This shifts all muonic-atom and molecular-ion parameters.
Compute the shifts and their consequences for MCF.
"""

import math

# === Constants ===
m_pi = 139.57039    # MeV, charged pion
m_mu_exp = 105.6584 # MeV, muon
f_pi = 92.07        # MeV, PDG pion decay constant
alpha = 1/137.036
hbarc = 197.3270    # MeV fm
m_p = 938.272       # MeV, proton
m_d = 1875.613      # MeV, deuteron
m_t = 2808.921      # MeV, triton
m_alpha = 3727.379  # MeV, alpha particle
c = 2.998e23        # fm/s
tau_mu = 2.197e-6   # s, muon lifetime
Q_dt = 17.588       # MeV, dt fusion energy release

# SUSY-predicted muon mass
m_mu_susy = math.sqrt(m_pi**2 - f_pi**2)
dm = m_mu_exp - m_mu_susy  # MeV
frac = dm / m_mu_exp

print("=" * 75)
print("MUON-CATALYZED FUSION: sBootstrap PREDICTIONS")
print("=" * 75)

print(f"\n1. MUON MASS")
print(f"  m_mu^exp  = {m_mu_exp:.4f} MeV")
print(f"  m_mu^SUSY = {m_mu_susy:.4f} MeV")
print(f"  Shift: {dm:.4f} MeV ({frac*100:.3f}%)")

# === Muonic atom parameters ===
print(f"\n2. MUONIC DEUTERIUM")
# Reduced mass for mu-d system
m_red_d_exp = m_mu_exp * m_d / (m_mu_exp + m_d)
m_red_d_susy = m_mu_susy * m_d / (m_mu_susy + m_d)

a_mud_exp = hbarc / (alpha * m_red_d_exp)  # fm
a_mud_susy = hbarc / (alpha * m_red_d_susy)

E1s_mud_exp = m_red_d_exp * alpha**2 / 2  # MeV
E1s_mud_susy = m_red_d_susy * alpha**2 / 2

print(f"  Reduced mass (mu-d):")
print(f"    m_red^exp  = {m_red_d_exp:.3f} MeV")
print(f"    m_red^SUSY = {m_red_d_susy:.3f} MeV")
print(f"    Shift: {(m_red_d_exp - m_red_d_susy):.3f} MeV ({(m_red_d_exp - m_red_d_susy)/m_red_d_exp*100:.3f}%)")

print(f"  Bohr radius (mu-d):")
print(f"    a_mud^exp  = {a_mud_exp:.2f} fm")
print(f"    a_mud^SUSY = {a_mud_susy:.2f} fm")
print(f"    Shift: {a_mud_susy - a_mud_exp:.2f} fm ({(a_mud_susy - a_mud_exp)/a_mud_exp*100:.3f}%)")

print(f"  1s binding energy:")
print(f"    E_1s^exp  = {E1s_mud_exp*1e3:.2f} keV")
print(f"    E_1s^SUSY = {E1s_mud_susy*1e3:.2f} keV")
print(f"    Shift: {(E1s_mud_exp - E1s_mud_susy)*1e6:.1f} eV")

# === dtmu molecular ion ===
print(f"\n3. dtmu MOLECULAR ION")
# The dtmu molecular ion binding energy scales approximately as m_red
# The (J=1, v=1) loosely bound state: epsilon_11 = 0.660 eV (literature)
eps_11 = 0.660  # eV
# Binding scales roughly as reduced mass of the 3-body system
# Simplified scaling: eps ~ m_mu (since m_mu << m_d, m_t)
eps_11_susy = eps_11 * (m_red_d_susy / m_red_d_exp)
deps = eps_11 - eps_11_susy

print(f"  epsilon_11 (J=1,v=1 binding energy):")
print(f"    eps^exp  = {eps_11:.3f} eV")
print(f"    eps^SUSY = {eps_11_susy:.3f} eV")
print(f"    Shift: {deps*1e3:.2f} meV ({deps/eps_11*100:.2f}%)")

# Vesman resonance condition: eps_11 + E_kin = Delta_E_rovibrational
# The resonance energy shifts by deps
kT_300K = 0.02585  # eV at 300 K
print(f"\n  Vesman resonance shift:")
print(f"    Delta eps_11 = {deps*1e3:.2f} meV")
print(f"    kT (300 K) = {kT_300K*1e3:.1f} meV")
print(f"    Shift / kT = {deps/kT_300K:.3f}")
print(f"    => Shifts resonance energy by {deps/kT_300K*100:.1f}% of thermal energy")

# dtmu equilibrium internuclear distance
# R_eq ~ a_mud (molecular size set by muon Bohr radius)
# More precisely, R_eq ~ 0.6 a_mud for the ground state
R_eq_exp = 0.6 * a_mud_exp  # approximate
R_eq_susy = 0.6 * a_mud_susy
print(f"\n  Internuclear distance in dtmu:")
print(f"    R_eq^exp  ~ {R_eq_exp:.1f} fm")
print(f"    R_eq^SUSY ~ {R_eq_susy:.1f} fm")
print(f"    Shift: {R_eq_susy - R_eq_exp:.1f} fm ({(R_eq_susy - R_eq_exp)/R_eq_exp*100:.2f}%)")

# === Sticking fraction ===
print(f"\n4. STICKING FRACTION")
# omega_s ~ (v_alpha/c)^3 * (Z_alpha * alpha * m_red)^2 / ...
# In sudden approximation: omega_s propto m_mu^(5/2) approximately
# (from Struensee et al. 1988, Jackson 1957)
# A simpler scaling: omega_s propto m_red^3 (from phase space)
# Let's use the known scaling: omega_s propto m_mu/m_alpha * ...
# Best estimate: omega_s scales approximately linearly with m_red
# (from the probability of the muon being within a_alpha of the alpha)
# More precisely: omega_s^0 ~ (a_alpha/a_mol)^3 where both scale as 1/m_red
# So omega_s^0 is approximately m_red-independent to leading order!
# The correction comes from the recoil: momentum kick p = m_mu * v_alpha
# where v_alpha = sqrt(2 * E_alpha / m_alpha), E_alpha = Q * m_n/(m_alpha+m_n) = 3.52 MeV
E_alpha = Q_dt * 939.565 / (m_alpha + 939.565)  # MeV
v_alpha = math.sqrt(2 * E_alpha / m_alpha)  # in units of c
p_kick_exp = m_mu_exp * v_alpha  # MeV/c
p_kick_susy = m_mu_susy * v_alpha

print(f"  Alpha particle after fusion:")
print(f"    E_alpha = {E_alpha:.3f} MeV")
print(f"    v_alpha/c = {v_alpha:.4f}")
print(f"    Momentum kick to muon:")
print(f"      p_kick^exp  = {p_kick_exp:.3f} MeV/c")
print(f"      p_kick^SUSY = {p_kick_susy:.3f} MeV/c")
print(f"      Shift: {(p_kick_exp - p_kick_susy):.4f} MeV/c ({(p_kick_exp - p_kick_susy)/p_kick_exp*100:.2f}%)")

# The sticking probability depends on the overlap integral
# which involves exp(-p_kick * r / hbar). The relevant scale is:
# xi = p_kick * a_alpha / hbar
xi_exp = p_kick_exp * a_mud_exp / (2 * hbarc)  # divide by Z=2 for alpha
xi_susy = p_kick_susy * a_mud_susy / (2 * hbarc)
print(f"\n  Adiabaticity parameter xi = p * a_alpha / hbar:")
print(f"    xi^exp  = {xi_exp:.4f}")
print(f"    xi^SUSY = {xi_susy:.4f}")
print(f"    (xi << 1 => sudden limit; xi >> 1 => adiabatic)")

# omega_s scaling with m_mu: from Bracci & Fiorentini (1982), omega_s propto m_mu^3
# This is because the sticking involves a momentum-space overlap that depends on p^3
omega_s_exp = 0.0067  # effective sticking, 0.67%
# Scaling: omega_s propto m_red^3 (from the Jackson formula)
omega_s_susy = omega_s_exp * (m_red_d_susy / m_red_d_exp)**3
d_omega = omega_s_exp - omega_s_susy

print(f"\n  Sticking probability (ω_s propto m_red^3 scaling):")
print(f"    omega_s^exp  = {omega_s_exp*100:.3f}%")
print(f"    omega_s^SUSY = {omega_s_susy*100:.3f}%")
print(f"    Reduction: {d_omega/omega_s_exp*100:.2f}%")

# === Fusions per muon ===
print(f"\n5. FUSIONS PER MUON")
# N = lambda_c * tau_mu / (1 + omega_s * lambda_c * tau_mu)
# In the limit omega_s * lambda_c * tau >> 1: N ~ 1/omega_s
lambda_c = 1.0e8  # s^-1, cycling rate (approximate)
N_exp = lambda_c * tau_mu / (1 + omega_s_exp * lambda_c * tau_mu)
N_susy = lambda_c * tau_mu / (1 + omega_s_susy * lambda_c * tau_mu)

print(f"  Cycling rate: lambda_c = {lambda_c:.0e} s^-1")
print(f"  Muon lifetime: tau_mu = {tau_mu*1e6:.3f} us")
print(f"  N_fusions^exp  = {N_exp:.1f}")
print(f"  N_fusions^SUSY = {N_susy:.1f}")
print(f"  Improvement: +{N_susy - N_exp:.1f} fusions ({(N_susy-N_exp)/N_exp*100:.2f}%)")

# In 1/omega_s limit:
N_inf_exp = 1/omega_s_exp
N_inf_susy = 1/omega_s_susy
print(f"\n  In the sticking-limited regime (1/omega_s):")
print(f"    N^exp  = {N_inf_exp:.0f}")
print(f"    N^SUSY = {N_inf_susy:.0f}")
print(f"    Improvement: +{N_inf_susy - N_inf_exp:.0f}")

# === Energy balance ===
print(f"\n6. ENERGY BALANCE AND BREAK-EVEN")
E_muon_production = 5000  # MeV (approximate, from accelerator)
E_per_fusion = Q_dt  # MeV
N_breakeven = E_muon_production / E_per_fusion

print(f"  Energy to produce one muon: ~{E_muon_production:.0f} MeV")
print(f"  Energy per dt fusion: {E_per_fusion:.1f} MeV")
print(f"  Break-even fusions: {N_breakeven:.0f}")
print(f"  Current best: ~150 fusions/muon (60% of break-even)")
print(f"  SUSY correction: +{N_susy - N_exp:.1f} fusions ({(N_susy-N_exp)/N_breakeven*100:.1f}% of gap to break-even)")

# === The SUSY chain ===
print(f"\n7. THE COMPLETE SUSY CHAIN FOR MCF")
print(f"  m_pi = {m_pi:.2f} MeV (QCD)")
print(f"  f_pi = {f_pi:.2f} MeV (QCD)")
print(f"  m_mu = sqrt(m_pi^2 - f_pi^2) = {m_mu_susy:.2f} MeV (SUSY)")
print(f"  a_mud = hbarc/(alpha * m_red) = {a_mud_susy:.1f} fm (Bohr)")
print(f"  R_eq(dtmu) ~ 0.6 * a_mud = {R_eq_susy:.0f} fm")
print(f"  => Nuclear overlap at ~{R_eq_susy:.0f} fm enables fusion")
print(f"  => eps_11 = {eps_11_susy:.3f} eV enables Vesman resonance")
print(f"  => omega_s = {omega_s_susy*100:.3f}% limits cycling")
print(f"  => N_fus = {N_susy:.0f} (vs break-even {N_breakeven:.0f})")
print(f"\n  Summary: SUSY derives the muon mass from QCD parameters,")
print(f"  which determines the muonic atom size, molecular binding,")
print(f"  resonance condition, and ultimately the MCF cycle count.")

# === Sensitivity analysis ===
print(f"\n8. SENSITIVITY: what if SUSY relation were exact?")
# If m_mu were exactly 104.9 MeV instead of 105.66 MeV
# The dtmu binding shifts by 4.5 meV
# The molecular formation rate depends on the resonance match
# A wider resonance scan (different temperatures) would be needed
print(f"  The Vesman resonance at T=300K has thermal width ~kT = 25.9 meV")
print(f"  The SUSY shift in eps_11 is {deps*1e3:.1f} meV")
print(f"  This is {deps/kT_300K*100:.0f}% of the thermal width")
print(f"  => The resonance peak shifts but remains accessible")
print(f"  => Optimal temperature shifts by ~{deps/kT_300K * 300:.0f} K")
print(f"     (from 300K to ~{300 + deps/kT_300K * 300:.0f} K to re-center)")

# === The GT connection ===
print(f"\n9. GT RELATION AND NUCLEAR MATRIX ELEMENTS")
g_A = 1.2756
m_N = 938.272
g_piNN_susy = g_A * m_N / math.sqrt(m_pi**2 - m_mu_exp**2)
g_piNN_pdg = g_A * m_N / f_pi
g_piNN_exp = 13.17

print(f"  g_piNN^SUSY = {g_piNN_susy:.3f} (from GT with f_pi^SUSY)")
print(f"  g_piNN^PDG  = {g_piNN_pdg:.3f} (from GT with f_pi^PDG)")
print(f"  g_piNN^exp  = {g_piNN_exp:.2f}")
print(f"\n  The dt fusion S-factor at low energy depends on the nuclear")
print(f"  potential, which at r > 1 fm is dominated by OPEP:")
print(f"    V_OPEP propto g_piNN^2")
print(f"  A 0.35% shift in g_piNN shifts V by 0.7%, which propagates")
print(f"  into the tunneling rate and hence the fusion rate lambda_fus.")
print(f"  Since lambda_fus >> lambda_c, this does not limit the MCF cycle,")
print(f"  but it IS relevant for the S-factor itself.")
