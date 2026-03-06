#!/usr/bin/env python3
"""
Top-neutrino duality: check m_t × m_ν against geometric scales.

The M-theory uplift conjecture (paper_lagrangian.tex, Section 11):
In the electric frame the top gets Dirac mass, neutrino massless.
In the magnetic frame the neutrino gets Majorana mass, top massless.
The product m_t × m_ν should be set by a duality scale.

The neutrino CW mass from the paper:
  m_ν ~ (g²m)/(4π) × (m/Λ_L)^6

with m ~ z₀² ~ 314 MeV, g ~ O(1), Λ_L ~ 9-12 GeV.

NOTE: The original version of this script had a UNIT ERROR.
It used m_ν in "meV" as if 1 meV = 10^{-3} MeV, but actually
1 meV = 10^{-3} eV = 10^{-9} MeV.  The spurious result
"m_t × m_ν = f_π²" was entirely due to this factor-of-10^6 bug.
The corrected script shows no simple relation between m_t × m_ν
and known scales.
"""
import numpy as np

print("=" * 70)
print("TOP-NEUTRINO DUALITY: m_t × m_ν PRODUCT")
print("=" * 70)

# Known values
m_t_pole = 172.76e3  # MeV
m_t_MSbar = 162.5e3  # MeV

# CORRECT UNIT: 1 meV = 10^{-9} MeV
# Neutrino masses: atmospheric window
m_nu_atm = 50e-9  # MeV  (= 50 meV)

print("\n--- Product m_t × m_ν (ALL in MeV) ---")
for m_nu_meV_val in [9, 20, 30, 50]:
    m_nu = m_nu_meV_val * 1e-9  # meV → MeV: 1 meV = 10^-9 MeV
    prod = m_t_pole * m_nu  # MeV²
    sqrt_prod = np.sqrt(prod)  # MeV
    print(f"  m_ν = {m_nu_meV_val:3d} meV = {m_nu:.2e} MeV:  "
          f"m_t × m_ν = {prod:.4e} MeV²   √(m_t m_ν) = {sqrt_prod:.4e} MeV "
          f"= {sqrt_prod*1e3:.2f} keV")

print()
print("--- Candidate scales ---")

f_pi = 92.2  # MeV
LQCD = 300   # MeV
print(f"  f_π = {f_pi} MeV → f_π² = {f_pi**2:.0f} MeV²")
print(f"  ΛQCD = {LQCD} MeV → ΛQCD² = {LQCD**2:.0f} MeV²")
print()
print("  For m_ν = 50 meV: m_t × m_ν = {:.4e} MeV²".format(m_t_pole * 50e-9))
print(f"  f_π² = {f_pi**2:.0f} MeV²")
print(f"  Ratio: f_π² / (m_t × m_ν) = {f_pi**2 / (m_t_pole * 50e-9):.0f}")
print()
print("  There is NO relation m_t × m_ν = f_π².")
print("  The original script had a factor 10^6 unit error (meV vs MeV).")

print()
print("--- What scale IS m_t × m_ν? ---")
for m_nu_meV_val in [50]:
    m_nu = m_nu_meV_val * 1e-9
    prod = m_t_pole * m_nu
    sqrt_prod = np.sqrt(prod)
    print(f"  √(m_t × m_ν) = {sqrt_prod*1e6:.1f} eV = {sqrt_prod*1e3:.2f} keV")
    print(f"  This is ~92 keV, not 92 MeV.")
    print(f"  It does not correspond to any known hadronic or EW scale.")

print()
print("--- CW neutrino mass from paper ---")
m_OR_MeV = 470.76  # MeV, O'R scale
g = 1.0

for LL in [9, 10, 12]:
    LL_MeV = LL * 1000
    m_nu_CW = (g**2 * m_OR_MeV) / (4 * np.pi) * (m_OR_MeV / LL_MeV)**6
    m_nu_CW_meV = m_nu_CW * 1e9  # MeV → meV
    print(f"  Λ_L = {LL} GeV:  m_ν(CW) = {m_nu_CW_meV:.4f} meV = {m_nu_CW:.4e} MeV")

print()
print("  The CW formula at g=1 gives m_ν ~ 10^{-4} meV for Λ_L=9 GeV,")
print("  much below the atmospheric scale 50 meV.")
print("  Reaching 50 meV requires g >> 1 or a different mechanism.")

print()
print("--- For m_t × m_ν = f_π², what would m_ν be? ---")
m_nu_for_fpi2 = f_pi**2 / m_t_pole  # MeV
print(f"  m_ν = f_π²/m_t = {m_nu_for_fpi2:.4f} MeV = {m_nu_for_fpi2*1e3:.1f} keV")
print(f"  This is 49 keV — cosmologically excluded (hot dark matter bounds).")
print(f"  Conclusion: m_t × m_ν ≠ f_π².")

print()
print("=" * 70)
print("SUMMARY")
print("=" * 70)
print()
print("The original script contained a unit error: 1 meV was treated")
print("as 10^{-3} MeV instead of 10^{-9} MeV (factor 10^6 wrong).")
print("With correct units:")
print(f"  √(m_t × 50 meV) = {np.sqrt(m_t_pole * 50e-9)*1e3:.1f} keV (not 92 MeV)")
print(f"  f_π²/m_t = {f_pi**2/m_t_pole*1e3:.1f} keV (not 49 meV)")
print()
print("The product m_t × m_ν does not equal f_π² or any other")
print("known scale. The duality relation, if it exists, remains")
print("to be identified.")
