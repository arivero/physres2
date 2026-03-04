#!/usr/bin/env python3
"""
Computational physics analysis of lepton and meson mass coincidences.
Uses PDG 2024 values.
"""
import numpy as np
from itertools import combinations
from math import comb

# ============================================================
# PDG 2024 masses (MeV)
# ============================================================
# Leptons
m_e = 0.51099895
m_mu = 105.6583755
m_tau = 1776.86

# Pseudoscalar mesons
m_pi0 = 134.9768
m_pipm = 139.57039
m_K0 = 497.611
m_Kpm = 493.677
m_eta = 547.862
m_etap = 957.78

# Gell-Mann-Okubo octet eta
m_eta8 = np.sqrt((4*m_K0**2 - m_pi0**2)/3)
print(f"m_eta8 = {m_eta8:.6f} MeV")

# Z boson
m_Z = 91187.6  # MeV

# ============================================================
# SET A: Mixed lepton products
# sqrt(m_a) * sqrt(m_b - m_c) for all {a,b,c} in {e,mu,tau}
# with m_b > m_c
# ============================================================
leptons = {'e': m_e, 'mu': m_mu, 'tau': m_tau}
lepton_names = ['e', 'mu', 'tau']

set_A = {}
print("\n" + "="*60)
print("SET A: sqrt(m_a) * sqrt(m_b - m_c)")
print("="*60)

for trio in combinations(lepton_names, 3):
    masses = [(n, leptons[n]) for n in trio]
    for i, (na, ma) in enumerate(masses):
        others = [masses[j] for j in range(3) if j != i]
        others.sort(key=lambda x: x[1], reverse=True)
        nb, mb = others[0]
        nc, mc = others[1]
        val = np.sqrt(ma) * np.sqrt(mb - mc)
        label = f"sqrt(m_{na}) * sqrt(m_{nb} - m_{nc})"
        set_A[label] = val
        print(f"  {label} = {val:.6g} MeV")

# ============================================================
# SET B: Meson-lepton differences
# ============================================================
set_B = {}
print("\n" + "="*60)
print("SET B: Meson-lepton differences")
print("="*60)

quantities_B = [
    ("m_pi0 - m_mu", m_pi0 - m_mu),
    ("m_pipm - m_mu", m_pipm - m_mu),
    ("m_eta8 - m_pi0", m_eta8 - m_pi0),
    ("m_eta - m_pi0", m_eta - m_pi0),
]
for label, val in quantities_B:
    set_B[label] = val
    print(f"  {label} = {val:.6g} MeV")

# ============================================================
# SET C: de Vries' isospin relation
# ============================================================
print("\n" + "="*60)
print("SET C: de Vries' isospin relation")
print("="*60)
lhs = ((m_pipm / m_pi0) - 1)**2
rhs = m_mu / m_Z
print(f"  |(m_pipm / m_pi0) - 1|^2 = {lhs:.6g}")
print(f"  m_mu / m_Z               = {rhs:.6g}")
print(f"  Ratio LHS/RHS            = {lhs/rhs:.6g}")
print(f"  Discrepancy              = {abs(lhs - rhs)/max(lhs, rhs)*100:.4f}%")

# ============================================================
# SET D: Electromagnetic decay width scaling
# ============================================================
print("\n" + "="*60)
print("SET D: (Gamma_EM / M^3)^(-1/2) in GeV")
print("="*60)

particles_D = [
    ("pi0",    134.9768,   7.73e-6),      # Gamma(gamma gamma) = 7.73 eV
    ("eta",    547.862,    0.516e-3),      # Gamma(gamma gamma) ~ 0.516 keV
    ("omega",  782.66,     0.632e-3),      # Gamma(e+e-) ~ 0.632 keV
    ("J/psi",  3096.900,   5.53e-3),       # Gamma(e+e-) = 5.53 keV
    ("Z0",     91187.6,    83.984),        # Gamma(e+e-) = 83.984 MeV
]

print(f"  {'Particle':<10} {'M (MeV)':<12} {'Gamma_EM (MeV)':<18} {'(Gamma/M^3)^(-1/2) (GeV)'}")
print(f"  {'-'*10} {'-'*12} {'-'*18} {'-'*25}")
for name, M, Gamma_EM in particles_D:
    ratio = Gamma_EM / M**3   # MeV^(-2)
    inv_sqrt = 1.0 / np.sqrt(ratio)  # MeV
    inv_sqrt_GeV = inv_sqrt / 1000.0  # GeV
    print(f"  {name:<10} {M:<12.4f} {Gamma_EM:<18.6g} {inv_sqrt_GeV:<25.6g}")

# ============================================================
# SCAN: All pairs from Sets A and B, find |x-y|/max(x,y) < 1%
# ============================================================
print("\n" + "="*60)
print("COINCIDENCE SCAN: Sets A and B pairs with < 1% agreement")
print("="*60)

all_quantities = {}
all_quantities.update(set_A)
all_quantities.update(set_B)

names = list(all_quantities.keys())
values = [all_quantities[n] for n in names]
N = len(names)
n_pairs = comb(N, 2)

print(f"\nTotal quantities N = {N}")
print(f"Total pairwise comparisons C(N,2) = {n_pairs}")

coincidences = []
for i in range(N):
    for j in range(i+1, N):
        x, y = values[i], values[j]
        if max(abs(x), abs(y)) == 0:
            continue
        disc = abs(x - y) / max(abs(x), abs(y))
        if disc < 0.01:
            coincidences.append((names[i], x, names[j], y, disc*100))

coincidences.sort(key=lambda t: t[4])
k_found = len(coincidences)

print(f"\nCoincidences found at <1% level: {k_found}")
print()
print(f"{'#':<3} {'Expression 1':<40} {'Value 1 (MeV)':<16} {'Expression 2':<25} {'Value 2 (MeV)':<16} {'Discrepancy'}")
print(f"{'='*3} {'='*40} {'='*16} {'='*25} {'='*16} {'='*11}")
for idx, (n1, v1, n2, v2, disc) in enumerate(coincidences, 1):
    print(f"{idx:<3} {n1:<40} {v1:<16.6g} {n2:<25} {v2:<16.6g} {disc:.4f}%")

# ============================================================
# LOOK-ELSEWHERE CORRECTION: Monte Carlo
# ============================================================
print("\n" + "="*60)
print("LOOK-ELSEWHERE CORRECTION (Monte Carlo)")
print("="*60)

np.random.seed(42)
n_trials = 10000
threshold = 0.01

def simulate_trial():
    """Draw random masses and compute the same functional forms."""
    # 3 lepton masses, log-uniform in [0.1, 2000]
    log_leps = np.random.uniform(np.log(0.1), np.log(2000), 3)
    leps = np.sort(np.exp(log_leps))  # sorted ascending

    # Set A: 3 quantities
    setA_sim = [
        np.sqrt(leps[0]) * np.sqrt(leps[2] - leps[1]),
        np.sqrt(leps[1]) * np.sqrt(leps[2] - leps[0]),
        np.sqrt(leps[2]) * np.sqrt(leps[1] - leps[0]),
    ]

    # For Set B: draw 5 underlying masses to form 4 differences
    log_mesons = np.random.uniform(np.log(0.1), np.log(2000), 5)
    mesons = np.exp(log_mesons)
    setB_sim = [
        abs(mesons[0] - mesons[1]),
        abs(mesons[2] - mesons[1]),
        abs(mesons[3] - mesons[0]),
        abs(mesons[4] - mesons[0]),
    ]

    all_sim = setA_sim + setB_sim

    # Count coincidences at 1% level
    n_coinc = 0
    for i in range(len(all_sim)):
        for j in range(i+1, len(all_sim)):
            x, y = all_sim[i], all_sim[j]
            if max(abs(x), abs(y)) > 0:
                if abs(x - y) / max(abs(x), abs(y)) < threshold:
                    n_coinc += 1
    return n_coinc

# Run Monte Carlo
coinc_counts = np.array([simulate_trial() for _ in range(n_trials)])

print(f"\nMonte Carlo: {n_trials} trials")
print(f"N = {N} quantities, C(N,2) = {n_pairs} pairwise comparisons")
print(f"Actual coincidences found: k = {k_found}")
print(f"\nDistribution of coincidences under null hypothesis:")
max_k = max(coinc_counts.max(), k_found) + 2
for k_val in range(max_k):
    count = np.sum(coinc_counts == k_val)
    pct = count / n_trials * 100
    marker = " <-- observed" if k_val == k_found else ""
    if count > 0 or k_val <= k_found + 1:
        print(f"  k = {k_val}: {count:>5} trials ({pct:6.2f}%){marker}")

p_value = np.mean(coinc_counts >= k_found)
print(f"\nP(k >= {k_found}) = {p_value:.4f}")
if p_value > 0:
    from scipy.stats import norm
    n_sigma = norm.ppf(1 - p_value)
    print(f"  Equivalent significance: {n_sigma:.2f} sigma")
else:
    print(f"  No trials reached k >= {k_found} in {n_trials} trials")
    print(f"  P < 1/{n_trials} = {1/n_trials:.1e}")

# Also compute expected number of coincidences analytically
# For log-uniform masses, the probability of a single pair agreeing to 1%
# is approximately 2 * 0.01 = 0.02 (for ratio near 1)
# But this is rough; the MC is more reliable
p_single_approx = 2 * threshold  # rough estimate
expected_approx = n_pairs * p_single_approx
print(f"\nRough analytic estimate:")
print(f"  P(single pair agrees to 1%) ~ {p_single_approx:.3f}")
print(f"  Expected coincidences ~ {n_pairs} * {p_single_approx:.3f} = {expected_approx:.2f}")
print(f"  MC mean: {coinc_counts.mean():.2f} +/- {coinc_counts.std():.2f}")

print("\n" + "="*60)
print("FINAL SUMMARY")
print("="*60)
print(f"\n{k_found} coincidence(s) found at the <1% level among {n_pairs} pairs.")
print(f"Look-elsewhere corrected p-value: {p_value:.4f}")
if p_value > 0:
    print(f"Significance: {n_sigma:.2f} sigma")
