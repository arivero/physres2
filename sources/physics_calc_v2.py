#!/usr/bin/env python3
"""
Additional checks and cleaner presentation.
"""
import numpy as np

# ============================================================
# Verify key results with full precision
# ============================================================

# Leptons
m_e = 0.51099895
m_mu = 105.6583755
m_tau = 1776.86

# Mesons
m_pi0 = 134.9768
m_pipm = 139.57039
m_K0 = 497.611
m_Kpm = 493.677
m_eta = 547.862
m_etap = 957.78

m_eta8 = np.sqrt((4*m_K0**2 - m_pi0**2)/3)

print("="*70)
print("DETAILED RESULTS (all values in MeV unless stated)")
print("="*70)

print("\n--- Derived quantity ---")
print(f"  m_eta8 = sqrt((4*m_K0^2 - m_pi0^2)/3) = {m_eta8:.6f} MeV")
print(f"  Compare: m_eta = {m_eta:.6f} MeV")
print(f"  eta-eta8 mass difference = {m_eta - m_eta8:.6f} MeV")

print("\n--- SET A (6 sig figs) ---")
A1 = np.sqrt(m_e) * np.sqrt(m_tau - m_mu)
A2 = np.sqrt(m_mu) * np.sqrt(m_tau - m_e)
A3 = np.sqrt(m_tau) * np.sqrt(m_mu - m_e)
print(f"  sqrt(m_e)   * sqrt(m_tau - m_mu) = {A1:.6g} MeV")
print(f"  sqrt(m_mu)  * sqrt(m_tau - m_e)  = {A2:.6g} MeV")
print(f"  sqrt(m_tau) * sqrt(m_mu - m_e)   = {A3:.6g} MeV")

print("\n--- SET B (6 sig figs) ---")
B1 = m_pi0 - m_mu
B2 = m_pipm - m_mu
B3 = m_eta8 - m_pi0
B4 = m_eta - m_pi0
print(f"  m_pi0  - m_mu  = {B1:.6g} MeV")
print(f"  m_pipm - m_mu  = {B2:.6g} MeV")
print(f"  m_eta8 - m_pi0 = {B3:.6g} MeV")
print(f"  m_eta  - m_pi0 = {B4:.6g} MeV")

print("\n--- SET C: de Vries' isospin relation ---")
m_Z = 91187.6
lhs = ((m_pipm / m_pi0) - 1)**2
rhs = m_mu / m_Z
print(f"  [(m_pi+- / m_pi0) - 1]^2 = {lhs:.6g}")
print(f"  m_mu / m_Z                = {rhs:.6g}")
print(f"  Discrepancy               = {abs(lhs-rhs)/rhs*100:.4f}%")

print("\n--- SET D: EM decay width scaling ---")
print("  Testing whether (Gamma_EM / M^3) is approximately constant")
print("  i.e. (Gamma_EM / M^3)^(-1/2) should give similar values if scaling holds\n")

particles = [
    ("pi0",   134.9768,   7.73e-6,    "Gamma(gamma gamma) = 7.73 eV"),
    ("eta",   547.862,    0.516e-3,   "Gamma(gamma gamma) = 0.516 keV"),
    ("omega", 782.66,     0.632e-3,   "Gamma(e+e-) = 0.632 keV"),
    ("J/psi", 3096.900,   5.53e-3,    "Gamma(e+e-) = 5.53 keV"),
    ("Z0",    91187.6,    83.984,     "Gamma(e+e-) = 83.984 MeV"),
]

print(f"  {'Particle':<8} {'M (MeV)':<12} {'Gamma_EM':<24} {'Gamma/M^3 (MeV^-2)':<22} {'(Gamma/M^3)^(-1/2) (GeV)'}")
print(f"  {'-'*8} {'-'*12} {'-'*24} {'-'*22} {'-'*25}")
for name, M, G, desc in particles:
    r = G / M**3
    inv = 1.0 / np.sqrt(r) / 1000  # GeV
    print(f"  {name:<8} {M:<12.4f} {desc:<24} {r:<22.6g} {inv:<25.6g}")

# Note: pi0 and eta have very similar (Gamma/M^3)^(-1/2) ~ 564 GeV
# This is because Gamma(P -> gamma gamma) ~ alpha^2 * M^3 / (64 pi^3 f_pi^2)
# and f_pi ~ 92 MeV. Let's check:
f_pi = 92.07  # MeV, PDG pion decay constant
alpha = 1/137.036
Gamma_theory_pi0 = alpha**2 * m_pi0**3 / (64 * np.pi**3 * f_pi**2)
print(f"\n  Theory check (chiral anomaly):")
print(f"  Gamma(pi0->gg) = alpha^2 * m_pi0^3 / (64*pi^3*f_pi^2)")
print(f"                 = {Gamma_theory_pi0*1e6:.3f} eV  (vs PDG 7.73 eV)")
print(f"  The 1/(64 pi^3 f_pi^2) factor gives (Gamma/M^3)^(-1/2):")
fac = alpha**2 / (64 * np.pi**3 * f_pi**2)
print(f"  = 1/sqrt({fac:.6g}) MeV = {1/np.sqrt(fac)/1000:.4f} GeV")

print("\n" + "="*70)
print("TABLE OF COINCIDENCES (sorted by precision)")
print("="*70)

results = [
    ("sqrt(m_mu)*sqrt(m_tau - m_e)", A2, "sqrt(m_tau)*sqrt(m_mu - m_e)", A3),
    ("sqrt(m_mu)*sqrt(m_tau - m_e)", A2, "m_eta8 - m_pi0", B3),
    ("sqrt(m_e)*sqrt(m_tau - m_mu)", A1, "m_pi0 - m_mu", B1),
    ("sqrt(m_tau)*sqrt(m_mu - m_e)", A3, "m_eta8 - m_pi0", B3),
]

print(f"\n{'#':<3} {'Quantity 1':<35} {'Value (MeV)':<14} {'Quantity 2':<35} {'Value (MeV)':<14} {'Discrepancy'}")
print(f"{'-'*3} {'-'*35} {'-'*14} {'-'*35} {'-'*14} {'-'*11}")
for idx, (n1, v1, n2, v2) in enumerate(results, 1):
    disc = abs(v1 - v2) / max(v1, v2) * 100
    print(f"{idx:<3} {n1:<35} {v1:<14.6g} {n2:<35} {v2:<14.6g} {disc:.4f}%")

print("\nNote: Coincidence #1 is purely leptonic (Koide-like).")
print("      Coincidence #3 links the electron mass to the pion-muon mass gap.")
print("      Coincidences #2 and #4 link lepton geometric means to")
print("      the Gell-Mann-Okubo eta8 - pion mass gap.")
