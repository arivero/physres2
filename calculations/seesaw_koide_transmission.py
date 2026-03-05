#!/usr/bin/env python3
"""
Seesaw-Koide transmission analysis.

Pure algebra: when does the Koide ratio Q(a,b,c) = (a^2+b^2+c^2)/(a+b+c)^2
transmit through an inversion map z_k -> 1/z_k?

Parts (a)-(d) as specified.
"""

import numpy as np
from scipy.optimize import brentq

# ==============================================================================
# Utility
# ==============================================================================

def Q(a, b, c):
    """Koide ratio for three positive reals."""
    num = a**2 + b**2 + c**2
    den = (a + b + c)**2
    return num / den

def koide_param(delta, z0=1.0):
    """Standard Koide parametrization: z_k = z0(1 + sqrt(2)*cos(delta + 2*pi*k/3))."""
    return [z0 * (1 + np.sqrt(2) * np.cos(delta + 2*np.pi*k/3)) for k in range(3)]

# ==============================================================================
# Part (a): Leading correction Q(1/z) - 2/3 when Q(z) = 2/3
# ==============================================================================

def part_a():
    print("=" * 72)
    print("PART (a): Q(1/z) - 2/3 as a function of mass ratios")
    print("=" * 72)

    # Analytic derivation:
    # Q(z) = (z1^2 + z2^2 + z3^2) / (z1 + z2 + z3)^2
    # Q(1/z) = (1/z1^2 + 1/z2^2 + 1/z3^2) / (1/z1 + 1/z2 + 1/z3)^2
    #
    # Let S_n = sum z_k^n. Then:
    #   Q(z) = S_2 / S_1^2
    #   Q(1/z) = S_{-2} / S_{-1}^2
    #
    # Q(z) = 2/3 means S_2 = (2/3) S_1^2, i.e., 3*S_2 = 2*S_1^2.
    #
    # For Q(1/z):
    #   S_{-2} = sum 1/z_k^2 = (z2*z3)^2 + (z1*z3)^2 + (z1*z2)^2 / (z1*z2*z3)^2
    #          = e_2^2 - 2*e_1*e_3) / e_3^2    [where e_k are elementary symmetric polynomials]
    #
    # Actually, let's use Newton's identities more carefully.
    # Let e1 = z1+z2+z3, e2 = z1*z2+z1*z3+z2*z3, e3 = z1*z2*z3.
    # Then:
    #   S_{-1} = 1/z1 + 1/z2 + 1/z3 = e2/e3
    #   S_{-2} = 1/z1^2 + 1/z2^2 + 1/z3^2 = (e2^2 - 2*e1*e3) / e3^2
    #
    # So Q(1/z) = S_{-2}/S_{-1}^2 = (e2^2 - 2*e1*e3)/e3^2 / (e2/e3)^2
    #           = (e2^2 - 2*e1*e3) / e2^2
    #           = 1 - 2*e1*e3/e2^2
    #
    # Therefore:
    #   Q(1/z) = 1 - 2*e1*e3/e2^2
    #
    # And the correction:
    #   Q(1/z) - 2/3 = 1/3 - 2*e1*e3/e2^2
    #                 = (e2^2 - 6*e1*e3) / (3*e2^2)

    print()
    print("EXACT FORMULA (no approximation needed):")
    print()
    print("  Q(1/z) = 1 - 2*e1*e3 / e2^2")
    print()
    print("  where e1 = z1+z2+z3, e2 = z1*z2+z1*z3+z2*z3, e3 = z1*z2*z3")
    print()
    print("  Correction: Q(1/z) - 2/3 = (e2^2 - 6*e1*e3) / (3*e2^2)")
    print()
    print("  This vanishes iff e2^2 = 6*e1*e3.")
    print()

    # Verify numerically with Koide parametrization
    print("Numerical verification with Koide parametrization at various delta:")
    print(f"  {'delta':>10s}  {'Q(z)':>10s}  {'Q(1/z)':>10s}  {'Q(1/z)-2/3':>12s}  {'formula':>12s}")

    for delta in [0.1, 0.5, 1.0, 1.5, 2.0, np.pi, 3*np.pi/4]:
        zz = koide_param(delta)
        z1, z2, z3 = zz
        if all(z > 0 for z in zz):
            Qd = Q(z1, z2, z3)
            Qi = Q(1/z1, 1/z2, 1/z3)
            e1 = z1 + z2 + z3
            e2 = z1*z2 + z1*z3 + z2*z3
            e3 = z1*z2*z3
            formula = (e2**2 - 6*e1*e3) / (3*e2**2)
            print(f"  {delta:10.4f}  {Qd:10.6f}  {Qi:10.6f}  {Qi-2/3:12.6e}  {formula:12.6e}")

    print()

    # Condition for Q(1/z) close to 2/3:
    # Need e2^2 ≈ 6*e1*e3
    # In terms of mass ratios r = m2/m1, s = m3/m1 (with z_k = sqrt(m_k)):
    # e1 = sqrt(m1)(1 + sqrt(r) + sqrt(s))
    # e2 = m1(sqrt(r) + sqrt(s) + sqrt(rs))
    # e3 = m1^{3/2} * sqrt(rs)
    # e2^2 / (e1*e3) = (sqrt(r)+sqrt(s)+sqrt(rs))^2 / ((1+sqrt(r)+sqrt(s))*sqrt(rs))
    # This should equal 6 for Q(1/z) = 2/3.

    print("Condition for Q(1/z) = 2/3 in terms of mass ratios r = m2/m1, s = m3/m1:")
    print("  (sqrt(r) + sqrt(s) + sqrt(rs))^2 = 6*(1 + sqrt(r) + sqrt(s))*sqrt(rs)")
    print()

    # When one mass vanishes (z1 -> 0): e3 -> 0, so Q(1/z) -> 1, which is NOT 2/3.
    # But the seed case is z1 = 0, so we must be careful.
    # Actually for the seed z1 = 0, we can't compute 1/z1. The statement in the problem
    # says Q(0,a,b) = Q(0,1/a,1/b) -- the zero stays zero.

    print("Special case z1 = 0 (seed): The inversion is applied only to nonzero entries.")
    print("  Q(0, a, b) = (a^2+b^2)/(a+b)^2 = Q(0, 1/a, 1/b) = (1/a^2+1/b^2)/(1/a+1/b)^2")
    print("  Both equal 1 - 2ab/(a^2+b^2) ... wait, let me verify.")

    # Q(0,a,b) = (a^2+b^2)/(a+b)^2
    # Q(0,1/a,1/b) = (1/a^2+1/b^2)/(1/a+1/b)^2 = ((a^2+b^2)/(ab)^2) / ((a+b)/(ab))^2
    #              = (a^2+b^2)/(a+b)^2
    # YES, identical.

    a_test, b_test = 3.0, 7.0
    print(f"  Check: Q(0,{a_test},{b_test}) = {Q(0, a_test, b_test):.10f}")
    print(f"         Q(0,1/{a_test},1/{b_test}) = {Q(0, 1/a_test, 1/b_test):.10f}")
    print()

# ==============================================================================
# Part (b): Q(1/z) as function of delta in Koide parametrization
# ==============================================================================

def part_b():
    print("=" * 72)
    print("PART (b): Q(1/z) in Koide parametrization as function of delta")
    print("=" * 72)
    print()

    # z_k = z0(1 + sqrt(2)*cos(delta + 2*pi*k/3)), k=0,1,2
    # z0 cancels from Q(z) but NOT from Q(1/z).
    # Actually z0 also cancels from Q(1/z):
    #   Q(1/z) = sum(1/z_k^2) / (sum(1/z_k))^2
    #   1/z_k = 1/(z0 * f_k) where f_k = 1 + sqrt(2)*cos(delta + 2*pi*k/3)
    #   sum(1/z_k^2) = (1/z0^2) * sum(1/f_k^2)
    #   (sum(1/z_k))^2 = (1/z0)^2 * (sum(1/f_k))^2
    # So z0 cancels. Q(1/z) depends only on delta.

    print("z0 cancels from Q(1/z), so it depends only on delta.")
    print()

    # Compute Q(1/z) on a fine grid
    deltas = np.linspace(0.01, np.pi - 0.01, 10000)
    Q_inv = np.zeros_like(deltas)
    valid = np.ones_like(deltas, dtype=bool)

    for i, d in enumerate(deltas):
        zz = koide_param(d)
        if all(z > 0 for z in zz):
            Q_inv[i] = Q(1/zz[0], 1/zz[1], 1/zz[2])
        else:
            valid[i] = False

    # Find where Q(1/z) = 2/3
    # First, let's check at delta = 3*pi/4 (seed)
    d_seed = 3*np.pi/4
    zz_seed = koide_param(d_seed)
    print(f"At delta = 3*pi/4 (seed):")
    print(f"  z = ({zz_seed[0]:.6f}, {zz_seed[1]:.6f}, {zz_seed[2]:.6f})")
    if zz_seed[0] < 1e-10:
        print(f"  z_0 ~ 0 (seed has one vanishing mass)")
        # For seed, Q(1/z) is computed on the two nonzero entries
        if abs(zz_seed[0]) < 1e-10:
            Qi_seed = Q(0, 1/zz_seed[1], 1/zz_seed[2])
        else:
            Qi_seed = Q(1/zz_seed[0], 1/zz_seed[1], 1/zz_seed[2])
    else:
        Qi_seed = Q(1/zz_seed[0], 1/zz_seed[1], 1/zz_seed[2])
    print(f"  Q(z) = {Q(*zz_seed):.10f}")
    print(f"  Q(1/z) with zero entry preserved = {Q(0, 1/zz_seed[1], 1/zz_seed[2]):.10f}")
    print()

    # Now find zeros of Q(1/z) - 2/3 in the valid range
    # Need all z_k > 0. This requires |cos(...)| < 1/sqrt(2) for all k,
    # i.e., delta in (arccos(1/sqrt(2)), pi - arccos(1/sqrt(2))) roughly
    # Actually the positivity condition is more nuanced.

    # Let's find the valid range first
    def all_positive(d):
        zz = koide_param(d)
        return all(z > 0 for z in zz)

    # The three z_k are positive when 1 + sqrt(2)*cos(d + 2*pi*k/3) > 0 for all k
    # i.e., cos(d + 2*pi*k/3) > -1/sqrt(2) for all k
    # cos(x) > -1/sqrt(2) iff x not in (3*pi/4, 5*pi/4) mod 2*pi
    # So d + 2*pi*k/3 not in (3*pi/4, 5*pi/4) mod 2*pi for k=0,1,2

    print("Finding delta values where all three z_k > 0:")

    # Scan for valid range boundaries
    test_deltas = np.linspace(0, 2*np.pi, 100000)
    valid_mask = np.array([all_positive(d) for d in test_deltas])
    transitions = np.where(np.diff(valid_mask.astype(int)))[0]

    print("  Valid regions (all z_k > 0):")
    if len(transitions) > 0:
        # Print transition points
        for t in transitions:
            d_trans = (test_deltas[t] + test_deltas[t+1]) / 2
            entering = valid_mask[t+1]
            print(f"    delta = {d_trans:.6f} ({d_trans/np.pi:.4f}*pi) : {'entering' if entering else 'leaving'} valid region")
    print()

    # Now find Q(1/z) = 2/3 crossings
    def Q_inv_minus_twothirds(d):
        zz = koide_param(d)
        if any(z <= 0 for z in zz):
            return np.nan
        return Q(1/zz[0], 1/zz[1], 1/zz[2]) - 2/3

    # Scan for sign changes
    print("Scanning for delta where Q(1/z) = 2/3:")
    solutions = []

    # Use fine grid within valid regions
    fine_deltas = np.linspace(0.001, 2*np.pi - 0.001, 100000)
    prev_val = None
    prev_d = None

    for d in fine_deltas:
        zz = koide_param(d)
        if all(z > 0 for z in zz):
            val = Q_inv_minus_twothirds(d)
            if not np.isnan(val) and prev_val is not None and not np.isnan(prev_val):
                if prev_val * val < 0:
                    # Sign change - find root
                    try:
                        root = brentq(Q_inv_minus_twothirds, prev_d, d)
                        solutions.append(root)
                    except:
                        pass
            prev_val = val
            prev_d = d
        else:
            prev_val = None
            prev_d = None

    if solutions:
        print(f"  Found {len(solutions)} solution(s):")
        for s in solutions:
            zz = koide_param(s)
            print(f"    delta = {s:.10f} ({s/np.pi:.10f}*pi)")
            print(f"    z = ({zz[0]:.6f}, {zz[1]:.6f}, {zz[2]:.6f})")
            print(f"    Q(z) = {Q(*zz):.10f}")
            print(f"    Q(1/z) = {Q(1/zz[0], 1/zz[1], 1/zz[2]):.10f}")
            print()
    else:
        print("  NO solution found in the interior of the valid region.")
        print("  Q(1/z) = 2/3 only at the boundary (seed) where one z_k -> 0.")
        print()

    # Print Q(1/z) at a selection of delta values
    print("Table of Q(1/z) vs delta:")
    print(f"  {'delta':>10s} {'delta/pi':>10s}  {'Q(1/z)':>12s}  {'Q(1/z)-2/3':>14s}  {'valid':>6s}")
    for d_frac in np.arange(0.05, 1.0, 0.05):
        d = d_frac * np.pi
        zz = koide_param(d)
        if all(z > 0 for z in zz):
            Qi = Q(1/zz[0], 1/zz[1], 1/zz[2])
            print(f"  {d:10.6f} {d_frac:10.4f}  {Qi:12.8f}  {Qi-2/3:14.8e}  {'yes':>6s}")
        else:
            print(f"  {d:10.6f} {d_frac:10.4f}  {'---':>12s}  {'---':>14s}  {'no':>6s}")
    print()

    # Analysis near the boundary
    print("Behavior near boundary delta -> 3*pi/4 (seed):")
    for eps in [0.01, 0.001, 0.0001, 0.00001]:
        d = 3*np.pi/4 - eps
        zz = koide_param(d)
        if all(z > 0 for z in zz):
            Qi = Q(1/zz[0], 1/zz[1], 1/zz[2])
            print(f"  delta = 3pi/4 - {eps:.0e}:  z_min = {min(zz):.6e}  Q(1/z) = {Qi:.10f}  Q(1/z)-2/3 = {Qi-2/3:.6e}")
    print()

    # Key result: does Q(1/z) approach 2/3 as delta -> 3*pi/4?
    print("RESULT: As delta -> 3*pi/4 (one z_k -> 0), Q(1/z) -> ?")
    # When z_0 -> 0, 1/z_0 -> inf. The Q ratio is dominated by the largest entry:
    # Q(1/z_0, 1/z_1, 1/z_2) -> (1/z_0)^2 / (1/z_0)^2 = 1/3... no.
    # Q(A, b, c) with A >> b,c -> A^2 / A^2 = 1/3? No: (A^2+b^2+c^2)/(A+b+c)^2 -> 1/3?
    # No: -> A^2/A^2 * 1/(1+...)^2 -> 1.
    # Wait: (A^2)/(A)^2 = 1. And (A^2+b^2+c^2)/(A+b+c)^2 -> 1 as A -> inf. So Q -> 1/3.
    # Hmm, let me think again.
    # Q(A,b,c) = (A^2+b^2+c^2)/(A+b+c)^2 -> A^2/A^2 = 1 as A -> inf? No!
    # (A+b+c)^2 = A^2 + 2A(b+c) + (b+c)^2
    # So Q -> A^2/(A^2 + 2A(b+c) + ...) -> 1/(1 + 2(b+c)/A + ...) -> 1 as A -> inf.
    # No wait. That's Q -> 1/3 only if ALL three are equal. For one dominating: Q -> 1.
    # Hmm no. Let A -> inf:
    # Q = (A^2 + small)/(A + small)^2 ≈ A^2/A^2 = 1.
    # But we want Q = 1/3 to be the equal case. Q = (3a^2)/(3a)^2 = 3a^2/9a^2 = 1/3. Yes.
    # So for one dominant: Q -> 1, not 1/3.
    # But the "seed" with one ZERO gives Q(0,a,b) which is between 1/3 and 1.
    # And Q(0,1/a,1/b) = Q(0,a,b) (proven in part a).
    # The issue is that the limit as z_0 -> 0 of Q(1/z_0, 1/z_1, 1/z_2) is 1 (not 2/3).
    # But Q(0, 1/z_1, 1/z_2) = Q(0, z_1, z_2) = 2/3 at the seed.
    # So the "seed self-duality" holds only when we keep the zero as zero (don't invert it).

    print("  Q(1/z_0, 1/z_1, 1/z_2) -> 1 as z_0 -> 0 (one entry diverges).")
    print("  The seed self-duality Q(0,a,b) = Q(0,1/a,1/b) is a DIFFERENT operation:")
    print("  it keeps the zero fixed rather than inverting it.")
    print()

    return solutions

# ==============================================================================
# Part (c): Physical down-type quarks and seesaw consistency
# ==============================================================================

def part_c():
    print("=" * 72)
    print("PART (c): Seesaw consistency for down-type quarks")
    print("=" * 72)
    print()

    # Masses in MeV
    m_d = 4.67
    m_s = 93.4
    m_b = 4183.0

    # "Inverse" charges: xi_k = 1/sqrt(m_k)
    xi_d = 1/np.sqrt(m_d)
    xi_s = 1/np.sqrt(m_s)
    xi_b = 1/np.sqrt(m_b)

    Q_inv_phys = Q(xi_d, xi_s, xi_b)
    Q_dir_phys = Q(np.sqrt(m_d), np.sqrt(m_s), np.sqrt(m_b))

    print(f"Down-type quark masses: m_d = {m_d}, m_s = {m_s}, m_b = {m_b} MeV")
    print()
    print(f"Direct Koide Q(sqrt(m_d), sqrt(m_s), sqrt(m_b)) = {Q_dir_phys:.6f}")
    print(f"  Deviation from 2/3: {Q_dir_phys - 2/3:.6e}")
    print()
    print(f"Inverse Koide Q(1/sqrt(m_d), 1/sqrt(m_s), 1/sqrt(m_b)) = {Q_inv_phys:.6f}")
    print(f"  Deviation from 2/3: {Q_inv_phys - 2/3:.6e}")
    print()

    # Question: find M_k = C^2/m_k such that Q(sqrt(M_1), sqrt(M_2), sqrt(M_3)) = 2/3
    # sqrt(M_k) = C/sqrt(m_k) = C * xi_k
    # Q(C*xi_d, C*xi_s, C*xi_b) = Q(xi_d, xi_s, xi_b) since C cancels from Q.
    # So Q(sqrt(M)) = Q(1/sqrt(m)) = 0.6652.
    # This is NOT 2/3.

    print("Seesaw masses M_k = C^2/m_k:")
    print("  sqrt(M_k) = C/sqrt(m_k) = C * xi_k")
    print("  Q(sqrt(M)) = Q(C*xi_d, C*xi_s, C*xi_b) = Q(xi_d, xi_s, xi_b)  [C cancels]")
    print(f"  = {Q_inv_phys:.6f}")
    print()
    print(f"  This is NOT 2/3. Deviation = {Q_inv_phys - 2/3:.4e}")
    print()
    print("  CONCLUSION: There is NO value of C that makes Q(sqrt(M)) = 2/3 with")
    print("  M_k = C^2/m_k. The seesaw does not transmit Koide from direct to inverse")
    print("  masses for down-type quarks.")
    print()

    # However, let's check what C would be if we just asked for seesaw relation
    # For Seiberg seesaw: M_j = Lambda^2 * (m_1*m_2*m_3)^{1/3} / m_j
    # Let's compute what C would need to be
    geo_mean = (m_d * m_s * m_b)**(1/3)
    print(f"Geometric mean (m_d*m_s*m_b)^(1/3) = {geo_mean:.4f} MeV")
    print()

    # For each quark: M_j = Lambda^2 * geo_mean / m_j, so C^2 = Lambda^2 * geo_mean
    # is NOT the same as M_j = C^2/m_j. There's an extra geo_mean factor.
    # Actually M_j = (Lambda^2 * geo_mean) / m_j, so C^2 = Lambda^2 * geo_mean.
    # C = Lambda * geo_mean^{1/2}

    print("Seiberg seesaw form: M_j = Lambda^2 * (m_d*m_s*m_b)^{1/3} / m_j")
    print("  This IS of the form M_j = C^2/m_j with C^2 = Lambda^2 * (m_d*m_s*m_b)^{1/3}")
    print(f"  So C = Lambda * {np.sqrt(geo_mean):.4f} MeV^(1/2)")
    print()

    # Compute the seesaw masses for Lambda = 1 GeV (ISS scale)
    Lambda = 1000  # MeV
    C2 = Lambda**2 * geo_mean
    M_d = C2 / m_d
    M_s = C2 / m_s
    M_b = C2 / m_b

    print(f"For Lambda = {Lambda} MeV:")
    print(f"  C^2 = {C2:.2f} MeV^2")
    print(f"  C = {np.sqrt(C2):.2f} MeV")
    print(f"  M_d = {M_d:.2f} MeV,  M_s = {M_s:.2f} MeV,  M_b = {M_b:.2f} MeV")
    print(f"  Q(sqrt(M_d), sqrt(M_s), sqrt(M_b)) = {Q(np.sqrt(M_d), np.sqrt(M_s), np.sqrt(M_b)):.6f}")
    print()

    # Now let's ask the inverse question: what masses would HAVE Q(1/sqrt(m)) = 2/3?
    # Using Koide parametrization: 1/sqrt(m_k) = z0(1 + sqrt(2)*cos(delta + 2*pi*k/3))
    # Then m_k = 1/[z0(1 + sqrt(2)*cos(delta + 2*pi*k/3))]^2
    # We need to find delta, z0 to fit the down-type quarks.

    # Fit z_k = 1/sqrt(m_k) to Koide parametrization
    xi = [xi_d, xi_s, xi_b]  # ordered by size: xi_b < xi_s < xi_d
    xi_sorted = sorted(xi, reverse=True)  # [xi_d, xi_s, xi_b]

    z0_inv = sum(xi) / 3
    print(f"Inverse charges: xi_d = {xi_d:.6f}, xi_s = {xi_s:.6f}, xi_b = {xi_b:.6f}")
    print(f"z0 = mean(xi) = {z0_inv:.6f}")
    print(f"(xi_k - z0)/z0 = ({(xi_d-z0_inv)/z0_inv:.6f}, {(xi_s-z0_inv)/z0_inv:.6f}, {(xi_b-z0_inv)/z0_inv:.6f})")
    deviation_ratios = [(x - z0_inv)/z0_inv for x in xi]
    # For Koide: (xi_k - z0)/z0 = sqrt(2)*cos(delta + 2*pi*k/3)
    # Sum of cos(...) = 0 (check), sum of cos^2 = 3/2 (Koide condition)
    sum_dev_sq = sum(d**2 for d in deviation_ratios)
    print(f"Sum of squared deviations / z0^2 = {sum_dev_sq:.6f} (should be 3/2 = 1.5 for Koide)")
    print(f"  Corresponds to Q = {(1 + 2*sum_dev_sq/3) / (1 + 2*sum_dev_sq/3 + (1 - sum_dev_sq/1.5 + 1)):.6f}")
    # Actually Q = (sum z_k^2)/(sum z_k)^2 = (3*z0^2 + sum dev^2) / (3*z0)^2
    # = (3*z0^2*(1 + sum_dev_sq/3)) / (9*z0^2) -- wait no.
    # sum z_k^2 = sum (z0 + z0*d_k)^2 = z0^2 * sum(1 + d_k)^2 = z0^2 * (3 + 2*sum(d_k) + sum(d_k^2))
    # = z0^2 * (3 + 0 + sum_dev_sq) since sum(d_k) = 0
    # (sum z_k)^2 = (3*z0)^2 = 9*z0^2
    # Q = z0^2*(3 + sum_dev_sq) / (9*z0^2) = (3 + sum_dev_sq)/9
    Q_check = (3 + sum_dev_sq) / 9
    print(f"  Q = (3 + {sum_dev_sq:.6f})/9 = {Q_check:.6f}")
    print(f"  Direct computation: Q = {Q_inv_phys:.6f}")
    print()

    # How close is Q(1/sqrt(m)) to 2/3?
    print(f"SUMMARY: Q(1/sqrt(m_d), 1/sqrt(m_s), 1/sqrt(m_b)) = {Q_inv_phys:.6f}")
    print(f"  Deviation from 2/3: {Q_inv_phys - 2/3:.4e} ({(Q_inv_phys - 2/3)*100/(2/3):.2f}%)")
    print(f"  This is remarkably close: only {(Q_inv_phys - 2/3)/(2/3)*100:.2f}% deviation.")
    print()

# ==============================================================================
# Part (d): Seiberg seesaw transmits Q exactly
# ==============================================================================

def part_d():
    print("=" * 72)
    print("PART (d): Seiberg seesaw transmission")
    print("=" * 72)
    print()

    # Seiberg seesaw for N_f = N_c = 3:
    # M_j = Lambda^2 * (m_1*m_2*m_3)^{1/3} / m_j
    #
    # sqrt(M_j) = Lambda * (m_1*m_2*m_3)^{1/6} / sqrt(m_j)
    #           = Lambda * (m_1*m_2*m_3)^{1/6} * (1/sqrt(m_j))
    #
    # So sqrt(M_j) = A / sqrt(m_j) where A = Lambda * (m_1*m_2*m_3)^{1/6}
    # is a CONSTANT (independent of j).
    #
    # Therefore:
    # Q(sqrt(M_1), sqrt(M_2), sqrt(M_3)) = Q(A/sqrt(m_1), A/sqrt(m_2), A/sqrt(m_3))
    #                                     = Q(1/sqrt(m_1), 1/sqrt(m_2), 1/sqrt(m_3))
    # since the constant A cancels from the ratio.
    #
    # QED: Q(sqrt(M)) = Q(1/sqrt(m)) exactly.

    print("PROOF that Q(sqrt(M)) = Q(1/sqrt(m)):")
    print()
    print("  Seiberg seesaw: M_j = Lambda^2 * (m_1*m_2*m_3)^{1/3} / m_j")
    print()
    print("  sqrt(M_j) = Lambda * (m_1*m_2*m_3)^{1/6} / sqrt(m_j)")
    print("            = A / sqrt(m_j)")
    print()
    print("  where A = Lambda * (m_1*m_2*m_3)^{1/6} is j-independent.")
    print()
    print("  Q is homogeneous of degree 0:")
    print("    Q(A*x, A*y, A*z) = Q(x, y, z) for any A > 0.")
    print()
    print("  Therefore: Q(sqrt(M_1), sqrt(M_2), sqrt(M_3))")
    print("           = Q(A/sqrt(m_1), A/sqrt(m_2), A/sqrt(m_3))")
    print("           = Q(1/sqrt(m_1), 1/sqrt(m_2), 1/sqrt(m_3))    [QED]")
    print()

    # Now prove Q(1/sqrt(m)) != Q(sqrt(m)) unless one mass vanishes
    print("PROOF that Q(1/sqrt(m)) != Q(sqrt(m)) generically:")
    print()
    print("  From part (a): Q(1/z) = 1 - 2*e1*e3/e2^2")
    print("  And:            Q(z) = S_2/S_1^2 where S_n = sum z_k^n")
    print()
    print("  Q(z) = 2/3 means 3*S_2 = 2*S_1^2.")
    print("  Q(1/z) = 2/3 means 3*e2^2 = 2*(3*e2^2 + ... ) -- let's be precise.")
    print()
    print("  Q(1/z) = 2/3 iff 1 - 2*e1*e3/e2^2 = 2/3 iff e1*e3/e2^2 = 1/6")
    print("  i.e., e2^2 = 6*e1*e3.")
    print()
    print("  Both Q(z)=2/3 AND Q(1/z)=2/3 require:")
    print("    3*(z1^2+z2^2+z3^2) = 2*(z1+z2+z3)^2   ... (I)")
    print("    (z1*z2+z1*z3+z2*z3)^2 = 6*(z1+z2+z3)*(z1*z2*z3)  ... (II)")
    print()

    # Let's check this system. Use z_k = z0(1 + sqrt(2)*cos(delta + 2*pi*k/3))
    # which satisfies (I) for all delta.
    # Condition (II) becomes a condition on delta.

    print("Using Koide parametrization z_k = z0(1 + sqrt(2)*cos(delta + 2*pi*k/3)):")
    print("  (I) is automatically satisfied.")
    print("  (II) becomes a transcendental equation in delta.")
    print()

    # Compute e1, e2, e3 symbolically in terms of delta
    # e1 = 3*z0 (always)
    # e2 = sum_{j<k} z0^2 * (1+sqrt(2)*c_j)(1+sqrt(2)*c_k)
    #    = z0^2 * [3 + 2*sqrt(2)*sum(c_j) + 2*sum_{j<k} c_j*c_k]
    # sum c_j = sum cos(delta + 2*pi*k/3) = 0
    # sum_{j<k} c_j*c_k = [( sum c_j)^2 - sum c_j^2]/2 = [0 - 3/2]/2 = -3/4
    # So e2 = z0^2 * [3 + 0 + 2*(-3/4)] = z0^2 * [3 - 3/2] = 3*z0^2/2

    # e3 = z0^3 * prod(1 + sqrt(2)*c_k)
    # prod(1 + sqrt(2)*c_k) = 1 + sqrt(2)*sum(c_k) + 2*sum_{j<k}c_j*c_k + 2*sqrt(2)*c_0*c_1*c_2
    #    = 1 + 0 + 2*(-3/4) + 2*sqrt(2)*c_0*c_1*c_2
    #    = 1 - 3/2 + 2*sqrt(2)*c_0*c_1*c_2
    #    = -1/2 + 2*sqrt(2)*c_0*c_1*c_2
    #
    # c_0*c_1*c_2 = prod cos(delta + 2*pi*k/3) = (1/4)*cos(3*delta)
    # [standard identity: prod cos(x+2*pi*k/3) = (1/4)*cos(3x)]
    #
    # So e3 = z0^3 * [-1/2 + 2*sqrt(2)*(1/4)*cos(3*delta)]
    #       = z0^3 * [-1/2 + (sqrt(2)/2)*cos(3*delta)]
    #       = z0^3 * [-1/2 + cos(3*delta)/sqrt(2)]

    print("  e1 = 3*z0")
    print("  e2 = (3/2)*z0^2")
    print("  e3 = z0^3 * [-1/2 + cos(3*delta)/sqrt(2)]")
    print()
    print("  Condition (II): e2^2 = 6*e1*e3")
    print("    (9/4)*z0^4 = 6 * 3*z0 * z0^3 * [-1/2 + cos(3*delta)/sqrt(2)]")
    print("    9/4 = 18 * [-1/2 + cos(3*delta)/sqrt(2)]")
    print("    [-1/2 + cos(3*delta)/sqrt(2)] = 9/72 = 1/8")
    print("    cos(3*delta)/sqrt(2) = 1/8 + 1/2 = 5/8")
    print("    cos(3*delta) = 5*sqrt(2)/8")
    print()

    val = 5*np.sqrt(2)/8
    print(f"  cos(3*delta) = 5*sqrt(2)/8 = {val:.10f}")
    print(f"  But |cos(3*delta)| <= 1, and 5*sqrt(2)/8 = {val:.6f} < 1. So solutions exist!")
    print()

    delta_sol = np.arccos(val) / 3
    print(f"  3*delta = arccos(5*sqrt(2)/8) = {np.arccos(val):.10f} rad")
    print(f"  delta = {delta_sol:.10f} rad = {delta_sol/np.pi:.10f}*pi")
    print()

    # Verify
    zz = koide_param(delta_sol)
    print("  Verification at this delta:")
    print(f"    z = ({zz[0]:.10f}, {zz[1]:.10f}, {zz[2]:.10f})")
    if all(z > 0 for z in zz):
        Qd = Q(*zz)
        Qi = Q(1/zz[0], 1/zz[1], 1/zz[2])
        print(f"    Q(z) = {Qd:.12f}")
        print(f"    Q(1/z) = {Qi:.12f}")
        print(f"    Q(1/z) - 2/3 = {Qi - 2/3:.6e}")
    else:
        print(f"    Some z_k < 0! Not physical.")
    print()

    # All solutions: 3*delta = +/- arccos(5*sqrt(2)/8) + 2*pi*n
    print("  All solutions: delta = [+/- arccos(5*sqrt(2)/8) + 2*pi*n] / 3")
    print()
    alpha = np.arccos(val)
    print("  Enumerating solutions in [0, 2*pi):")
    all_sols = []
    for sign in [1, -1]:
        for n in range(-3, 4):
            d = (sign * alpha + 2*np.pi*n) / 3
            d_mod = d % (2*np.pi)
            if 0 <= d_mod < 2*np.pi:
                all_sols.append(d_mod)
    all_sols = sorted(set([round(d, 10) for d in all_sols]))

    print(f"  {'delta':>12s}  {'delta/pi':>12s}  {'all z>0':>8s}  {'Q(z)':>10s}  {'Q(1/z)':>10s}")
    for d in all_sols:
        zz = koide_param(d)
        pos = all(z > 0 for z in zz)
        if pos:
            Qd = Q(*zz)
            Qi = Q(1/zz[0], 1/zz[1], 1/zz[2])
            print(f"  {d:12.8f}  {d/np.pi:12.8f}  {'yes':>8s}  {Qd:10.6f}  {Qi:10.6f}")
        else:
            print(f"  {d:12.8f}  {d/np.pi:12.8f}  {'NO':>8s}  {'---':>10s}  {'---':>10s}")
    print()

    # Now verify that Q(z) != Q(1/z) unless one mass vanishes
    print("Proving Q(z) != Q(1/z) for generic z (all nonzero):")
    print()
    print("  We showed Q(z) = 2/3 for all delta (Koide parametrization).")
    print("  Q(1/z) = 2/3 requires cos(3*delta) = 5*sqrt(2)/8.")
    print("  At the seed delta = 3*pi/4: cos(3*3*pi/4) = cos(9*pi/4) = cos(pi/4) = 1/sqrt(2)")
    print(f"    = {1/np.sqrt(2):.10f}")
    print(f"    vs 5*sqrt(2)/8 = {val:.10f}")
    print(f"    NOT EQUAL. So the seed does NOT satisfy condition (II).")
    print()
    print("  But at the seed, z_0 = 0 (one mass vanishes), so 1/z_0 is undefined.")
    print("  The 'seed self-duality' Q(0,a,b) = Q(0,1/a,1/b) is a different statement")
    print("  that doesn't require condition (II).")
    print()

    # Check: when one z_k = 0, e3 = 0, and condition (II) becomes e2^2 = 0, so e2 = 0.
    # But e2 = sum z_i*z_j > 0 for positive z. So condition (II) is NEVER satisfied
    # when one z = 0 (with the other two positive).
    # This confirms that the "all three" inversion Q(1/z_1,1/z_2,1/z_3) = 2/3
    # is a DIFFERENT condition from the seed self-duality.

    print("  Key insight: Q(sqrt(M)) = Q(1/sqrt(m)) [seesaw transmits the INVERSE Koide].")
    print("  This is generically NOT the same as Q(sqrt(m)) [the DIRECT Koide].")
    print("  They coincide only at the special delta values where cos(3*delta) = 5*sqrt(2)/8.")
    print()

    # Numerical check with physical masses
    print("NUMERICAL CHECK with physical masses:")
    print()

    # Leptons
    m_e, m_mu, m_tau = 0.511, 105.658, 1776.86  # MeV
    Q_lep = Q(np.sqrt(m_e), np.sqrt(m_mu), np.sqrt(m_tau))
    Q_lep_inv = Q(1/np.sqrt(m_e), 1/np.sqrt(m_mu), 1/np.sqrt(m_tau))
    print(f"Leptons (e, mu, tau): {m_e}, {m_mu}, {m_tau} MeV")
    print(f"  Q(sqrt(m)) = {Q_lep:.8f}  (deviation from 2/3: {Q_lep-2/3:.4e})")
    print(f"  Q(1/sqrt(m)) = {Q_lep_inv:.8f}  (deviation from 2/3: {Q_lep_inv-2/3:.4e})")
    print()

    # Down quarks
    m_d, m_s, m_b = 4.67, 93.4, 4183.0
    Q_down = Q(np.sqrt(m_d), np.sqrt(m_s), np.sqrt(m_b))
    Q_down_inv = Q(1/np.sqrt(m_d), 1/np.sqrt(m_s), 1/np.sqrt(m_b))
    print(f"Down quarks (d, s, b): {m_d}, {m_s}, {m_b} MeV")
    print(f"  Q(sqrt(m)) = {Q_down:.8f}  (deviation from 2/3: {Q_down-2/3:.4e})")
    print(f"  Q(1/sqrt(m)) = {Q_down_inv:.8f}  (deviation from 2/3: {Q_down_inv-2/3:.4e})")
    print()

    # Koide triple (c, b, t)
    m_c, m_b2, m_t = 1280.0, 4183.0, 173000.0
    Q_cbt = Q(np.sqrt(m_c), np.sqrt(m_b2), np.sqrt(m_t))
    Q_cbt_inv = Q(1/np.sqrt(m_c), 1/np.sqrt(m_b2), 1/np.sqrt(m_t))
    print(f"(c, b, t): {m_c}, {m_b2}, {m_t} MeV")
    print(f"  Q(sqrt(m)) = {Q_cbt:.8f}  (deviation from 2/3: {Q_cbt-2/3:.4e})")
    print(f"  Q(1/sqrt(m)) = {Q_cbt_inv:.8f}  (deviation from 2/3: {Q_cbt_inv-2/3:.4e})")
    print()

    # Koide triple (-s, c, b)
    m_s2 = 93.4
    Q_scb = Q(np.sqrt(m_s2), np.sqrt(m_c), np.sqrt(m_b2))  # all positive; "-s" means the phase
    Q_scb_inv = Q(1/np.sqrt(m_s2), 1/np.sqrt(m_c), 1/np.sqrt(m_b2))
    print(f"(s, c, b): {m_s2}, {m_c}, {m_b2} MeV")
    print(f"  Q(sqrt(m)) = {Q_scb:.8f}  (deviation from 2/3: {Q_scb-2/3:.4e})")
    print(f"  Q(1/sqrt(m)) = {Q_scb_inv:.8f}  (deviation from 2/3: {Q_scb_inv-2/3:.4e})")
    print()

    # Seiberg meson masses
    print("Seiberg seesaw meson masses M_j = Lambda^2*(m_d*m_s*m_b)^{1/3}/m_j:")
    Lambda = 1000  # arbitrary, cancels from Q
    geo = (m_d * m_s * m_b)**(1.0/3)
    M = [Lambda**2 * geo / m for m in [m_d, m_s, m_b]]
    Q_meson = Q(np.sqrt(M[0]), np.sqrt(M[1]), np.sqrt(M[2]))
    print(f"  Lambda = {Lambda} MeV (arbitrary, cancels from Q)")
    print(f"  M_d = {M[0]:.2f}, M_s = {M[1]:.2f}, M_b = {M[2]:.2f} MeV")
    print(f"  Q(sqrt(M)) = {Q_meson:.8f}")
    print(f"  Q(1/sqrt(m)) = {Q_down_inv:.8f}")
    print(f"  Difference: {abs(Q_meson - Q_down_inv):.2e}  (should be exactly 0)")
    print()

    print("=" * 72)
    print("MASTER SUMMARY")
    print("=" * 72)
    print()
    print("(a) Q(1/z) = 1 - 2*e1*e3/e2^2 (exact, closed form).")
    print("    Q(1/z) = 2/3 iff e2^2 = 6*e1*e3.")
    print("    Seed self-duality is special: Q(0,a,b) = Q(0,1/a,1/b) identically.")
    print()
    print("(b) In Koide parametrization, Q(1/z) depends only on delta.")
    print("    Q(1/z) = 2/3 iff cos(3*delta) = 5*sqrt(2)/8 = 0.88388...")
    print(f"    Solution: delta = {delta_sol:.8f} = {delta_sol/np.pi:.8f}*pi")
    print("    The seed (delta = 3*pi/4) does NOT satisfy this (different mechanism).")
    print()
    print("(c) For physical down quarks:")
    print(f"    Q(1/sqrt(m)) = {Q_down_inv:.6f} (0.22% from 2/3).")
    print("    Seesaw M_k = C^2/m_k gives Q(sqrt(M)) = Q(1/sqrt(m)), so no new info.")
    print("    C is not determined by the Koide condition alone (cancels).")
    print()
    print("(d) Seiberg seesaw M_j = Lambda^2*(prod m)^{1/3}/m_j gives")
    print("    Q(sqrt(M)) = Q(1/sqrt(m)) EXACTLY (proven: overall scale cancels).")
    print("    This equals Q(sqrt(m)) only if cos(3*delta) = 5*sqrt(2)/8,")
    print("    which generically fails. The seesaw transmits the INVERSE Koide,")
    print("    not the direct Koide.")
    print()

# ==============================================================================
# MAIN
# ==============================================================================

if __name__ == "__main__":
    part_a()
    print()
    solutions_b = part_b()
    print()
    part_c()
    print()
    part_d()
