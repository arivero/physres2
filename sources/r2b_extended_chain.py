#!/usr/bin/env python3
"""Round 2-B: Extended mass chain computations"""
import numpy as np
from itertools import combinations

# Masses in MeV (PDG 2023)
masses = {
    'e': 0.51099895, 'mu': 105.6583755, 'tau': 1776.86,
    'u': 2.16, 'd': 4.67, 's': 93.4, 'c': 1270, 'b': 4180, 't': 172690
}
names = ['e', 'mu', 'tau', 'u', 'd', 's', 'c', 'b', 't']

def Q_val(m1, m2, m3, s1, s2, s3):
    """Compute Q for given masses and signs."""
    sq = s1*np.sqrt(m1) + s2*np.sqrt(m2) + s3*np.sqrt(m3)
    return sq**2 / (m1 + m2 + m3)

# ============================================================
# Task 1: Full ranked table of all 84 triplets
# ============================================================
print("=" * 80)
print("TASK 1: Full ranked table of all 84 triplets")
print("=" * 80)

results = []
for idx in combinations(range(9), 3):
    n = tuple(names[i] for i in idx)
    m = tuple(masses[names[i]] for i in idx)
    best_diff = float('inf')
    best_Q = 0
    best_signs = (1,1,1)
    for s1 in [1, -1]:
        for s2 in [1, -1]:
            for s3 in [1, -1]:
                Q = Q_val(m[0], m[1], m[2], s1, s2, s3)
                diff = abs(Q - 1.5)
                if diff < best_diff:
                    best_diff = diff
                    best_Q = Q
                    best_signs = (s1, s2, s3)
    results.append((best_diff, n, best_signs, best_Q))

results.sort()
print(f"{'Rank':>4} {'Triplet':<16} {'Signs':<8} {'Q':>14} {'|Q-3/2|':>12}")
print("-" * 60)
for rank, (diff, n, signs, Q) in enumerate(results, 1):
    sign_str = ''.join('+' if s > 0 else '-' for s in signs)
    print(f"{rank:>4} ({n[0]:>3},{n[1]:>3},{n[2]:>3}) {sign_str:<8} {Q:>14.10f} {diff:>12.4e}")

# ============================================================
# Task 2: Second scaling step
# ============================================================
print("\n" + "=" * 80)
print("TASK 2: Second scaling step M₂ = 3M₁, δ₂ = 3δ₁")
print("=" * 80)

# Lepton fit parameters
M0 = 313.8316 # MeV (will refine)

# Fit M0 and delta0 from leptons
# sqrt(m_k) = sqrt(M0) * (1 + sqrt(2) * cos(2*pi*k/3 + delta0))
# k = 0, 1, 2 for e, mu, tau
from scipy.optimize import minimize

def lepton_residual(params):
    M0, delta0 = params
    sqM0 = np.sqrt(M0)
    predicted = []
    for k in range(3):
        sqm = sqM0 * (1 + np.sqrt(2) * np.cos(2*np.pi*k/3 + delta0))
        predicted.append(sqm**2)
    actual = [masses['e'], masses['mu'], masses['tau']]
    return sum((p - a)**2 for p, a in zip(predicted, actual))

# Try without scipy first
# We know M0 ≈ 313.84, delta0 ≈ 2.31662
# Simple grid refinement
best_res = float('inf')
best_M0, best_d0 = 313.84, 2.31662
for M0_try in np.linspace(313.80, 313.90, 100):
    for d0_try in np.linspace(2.3160, 2.3175, 100):
        r = lepton_residual([M0_try, d0_try])
        if r < best_res:
            best_res = r
            best_M0, best_d0 = M0_try, d0_try

# Finer refinement
for M0_try in np.linspace(best_M0 - 0.01, best_M0 + 0.01, 1000):
    for d0_try in np.linspace(best_d0 - 0.001, best_d0 + 0.001, 1000):
        r = lepton_residual([M0_try, d0_try])
        if r < best_res:
            best_res = r
            best_M0, best_d0 = M0_try, d0_try

M0 = best_M0
delta0 = best_d0

print(f"Lepton fit: M₀ = {M0:.6f} MeV, δ₀ = {delta0:.8f} rad")
print(f"δ₀ mod 2π/3 = {delta0 % (2*np.pi/3):.8f} (cf. 2/9 = {2/9:.8f})")

# Check lepton predictions
print("\nLepton predictions:")
sqM0 = np.sqrt(M0)
for k, name in enumerate(['e', 'mu', 'tau']):
    sqm = sqM0 * (1 + np.sqrt(2) * np.cos(2*np.pi*k/3 + delta0))
    m_pred = sqm**2
    m_pdg = masses[name]
    print(f"  {name:>3}: predicted = {m_pred:.6f} MeV, PDG = {m_pdg:.6f} MeV, dev = {abs(m_pred-m_pdg)/m_pdg*100:.4f}%")

# First scaling: M1 = 3*M0, delta1 = 3*delta0
M1 = 3 * M0
delta1 = 3 * delta0
print(f"\nFirst scaling: M₁ = 3M₀ = {M1:.4f} MeV, δ₁ = 3δ₀ = {delta1:.6f} rad")
sqM1 = np.sqrt(M1)
print("Quark predictions (M₁, δ₁):")
quark_names_1 = ['s', 'c', 'b']
for k, name in enumerate(quark_names_1):
    sqm = sqM1 * (1 + np.sqrt(2) * np.cos(2*np.pi*k/3 + delta1))
    m_pred = sqm**2
    m_pdg = masses[name]
    print(f"  {name:>3}: predicted = {m_pred:.4f} MeV, PDG = {m_pdg:.4f} MeV, dev = {abs(m_pred-m_pdg)/m_pdg*100:.2f}%")

# Second scaling: M2 = 3*M1 = 9*M0, delta2 = 3*delta1 = 9*delta0
M2 = 3 * M1
delta2 = 3 * delta1
print(f"\nSecond scaling: M₂ = 9M₀ = {M2:.4f} MeV, δ₂ = 9δ₀ = {delta2:.6f} rad")
sqM2 = np.sqrt(M2)
print("Predicted masses (M₂, δ₂):")
for k in range(3):
    sqm = sqM2 * (1 + np.sqrt(2) * np.cos(2*np.pi*k/3 + delta2))
    m_pred = sqm**2
    print(f"  k={k}: m = {m_pred:.4f} MeV = {m_pred/1000:.4f} GeV")
    # Check Q for the triplet
print()

# Compute Q for second scaling triplet
m_second = []
for k in range(3):
    sqm = sqM2 * (1 + np.sqrt(2) * np.cos(2*np.pi*k/3 + delta2))
    m_second.append(sqm**2)
Q_second = Q_val(m_second[0], m_second[1], m_second[2], 1, 1, 1)
print(f"Q for second scaling triplet (all +): {Q_second:.15f}")
# Try with sign flip
for s1 in [1, -1]:
    for s2 in [1, -1]:
        for s3 in [1, -1]:
            Q = Q_val(m_second[0], m_second[1], m_second[2], s1, s2, s3)
            if abs(Q - 1.5) < 0.01:
                sign_str = ''.join('+' if s > 0 else '-' for s in [s1,s2,s3])
                print(f"  Signs {sign_str}: Q = {Q:.15f}, |Q-3/2| = {abs(Q-1.5):.2e}")

# Check if any match known particles
print("\nComparison to known particles:")
known = {'u': 2.16, 'd': 4.67, 's': 93.4, 'c': 1270, 'b': 4180, 't': 172690}
for k, m in enumerate(m_second):
    print(f"  k={k}: {m:.4f} MeV = {m/1000:.4f} GeV")
    for name, mpdg in known.items():
        ratio = m / mpdg
        if 0.5 < ratio < 2.0:
            print(f"    ~ {name} (PDG: {mpdg} MeV, ratio = {ratio:.3f})")

# ============================================================
# Task 3: Complete 4-branch catalog
# ============================================================
print("\n" + "=" * 80)
print("TASK 3: Complete 4-branch catalog for descending chain")
print("=" * 80)

def solve_middle_mass(m1, m3, outer_sign, inner_sign):
    """Solve Q(m1, m2, m3) = 3/2 for m2.

    (s1*sqrt(m1) + s2*sqrt(m2) + s3*sqrt(m3))^2 = 3/2 * (m1 + m2 + m3)

    With outer_sign for s1=s3 and inner_sign for s2:
    Let S = outer_sign*(sqrt(m1) + sqrt(m3)), and x = inner_sign*sqrt(m2)
    Then (S + x)^2 = 3/2 * (m1 + x^2 + m3)
    S^2 + 2Sx + x^2 = 3/2*m1 + 3/2*x^2 + 3/2*m3
    -1/2*x^2 + 2Sx + S^2 - 3/2*(m1+m3) = 0
    x^2 - 4Sx - 2S^2 + 3*(m1+m3) = 0

    Actually let me be more careful. The signs s1, s2, s3 are independent.
    Let me parameterize differently.

    For Q = 3/2: (a + b)^2 = 3/2 * (m1 + m2 + m3)
    where a = s1*sqrt(m1) + s3*sqrt(m3) and b = s2*sqrt(m2)

    a^2 + 2ab + b^2 = 3/2*m1 + 3/2*m2 + 3/2*m3
    Note b^2 = m2.
    a^2 + 2ab + m2 = 3/2*m1 + 3/2*m2 + 3/2*m3
    2ab = 3/2*m1 + 1/2*m2 + 3/2*m3 - a^2

    Let u = sqrt(m2), so m2 = u^2, b = s2*u.
    2a*s2*u = 3/2*m1 + 1/2*u^2 + 3/2*m3 - a^2
    1/2*u^2 - 2a*s2*u + (3/2*(m1+m3) - a^2) = 0
    u^2 - 4a*s2*u + (3*(m1+m3) - 2a^2) = 0

    u = [4a*s2 ± sqrt(16a^2 - 4*(3*(m1+m3) - 2a^2))] / 2
    u = 2a*s2 ± sqrt(4a^2 - 3*(m1+m3) + 2a^2)
    u = 2a*s2 ± sqrt(6a^2 - 3*(m1+m3))
    """
    results = []
    for s1 in [1, -1]:
        for s3 in [1, -1]:
            a = s1*np.sqrt(m1) + s3*np.sqrt(m3)
            disc = 6*a**2 - 3*(m1 + m3)
            if disc < 0:
                continue
            sqrt_disc = np.sqrt(disc)
            for s2 in [1, -1]:
                for pm in [1, -1]:
                    u = 2*a*s2 + pm*sqrt_disc
                    if u < 0:
                        continue  # sqrt(m2) must be non-negative
                    m2 = u**2
                    # Verify
                    Q = Q_val(m1, m2, m3, s1, s2, s3)
                    if abs(Q - 1.5) < 1e-8:
                        results.append({
                            'signs': (s1, s2, s3),
                            'm2': m2,
                            'Q_residual': abs(Q - 1.5),
                            'u': u
                        })
    # Remove duplicates (same m2 to 6 sig figs)
    unique = []
    seen = set()
    for r in results:
        key = round(r['m2'], 3)
        if key not in seen:
            seen.add(key)
            unique.append(r)
    return unique

mt = masses['t']  # 172690 MeV
mb = masses['b']  # 4180 MeV

print("\n--- Step 1: Q(m_b, m_?, m_t) = 3/2 ---")
step1 = solve_middle_mass(mb, mt, 1, 1)
step1.sort(key=lambda x: x['m2'])
for i, sol in enumerate(step1):
    s = ''.join('+' if x > 0 else '-' for x in sol['signs'])
    print(f"  Branch {i+1}: signs={s}, m = {sol['m2']:.6f} MeV = {sol['m2']/1000:.6f} GeV, |Q-3/2| = {sol['Q_residual']:.2e}")

print(f"\n  Found {len(step1)} distinct solutions")

print("\n--- Step 2: For each Step 1 result, solve Q(m_step1, m_?, m_b) = 3/2 ---")
all_step2 = []
for i, s1_sol in enumerate(step1):
    m_s1 = s1_sol['m2']
    step2 = solve_middle_mass(m_s1, mb, 1, 1)
    step2.sort(key=lambda x: x['m2'])
    print(f"\n  From Step 1 branch {i+1} (m = {m_s1:.4f} MeV):")
    for j, sol in enumerate(step2):
        s = ''.join('+' if x > 0 else '-' for x in sol['signs'])
        print(f"    Branch {j+1}: signs={s}, m = {sol['m2']:.6f} MeV, |Q-3/2| = {sol['Q_residual']:.2e}")
        all_step2.append({'parent': i, 'm_parent': m_s1, 'm2': sol['m2'], 'signs': sol['signs'], 'Q_res': sol['Q_residual']})

print(f"\n  Total Step 2 solutions: {len(all_step2)}")

print("\n--- Step 3: For interesting Step 2 results, solve Q(m_step2, m_?, m_step1) = 3/2 ---")
# Focus on the chain that gives physical masses
for s2_sol in all_step2:
    m_s2 = s2_sol['m2']
    m_s1 = s2_sol['m_parent']
    if m_s2 < 0.001 or m_s2 > 1e6:
        continue  # skip extreme values
    step3 = solve_middle_mass(m_s2, m_s1, 1, 1)
    step3.sort(key=lambda x: x['m2'])
    if len(step3) > 0:
        print(f"\n  Chain: m_t -> {m_s1:.2f} -> {m_s2:.4f} MeV:")
        for k, sol in enumerate(step3):
            s = ''.join('+' if x > 0 else '-' for x in sol['signs'])
            m = sol['m2']
            print(f"    Branch {k+1}: signs={s}, m = {m:.6f} MeV, |Q-3/2| = {sol['Q_residual']:.2e}")
            # Flag echoes
            for name, mpdg in masses.items():
                if abs(m - mpdg) / mpdg < 0.15:
                    print(f"      ** echoes {name} (PDG: {mpdg} MeV, dev = {abs(m-mpdg)/mpdg*100:.1f}%) **")

print("\n" + "=" * 80)
print("ALL TASKS COMPLETE")
print("=" * 80)
