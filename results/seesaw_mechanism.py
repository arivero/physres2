"""
Seesaw mechanism and dual Koide analysis.

In SU(3) SQCD with N_f = N_c = 3, the Seiberg duality maps quarks to
mesons with VEVs M^i_j = delta^i_j C/m_i, where C = Lambda^2 (prod m_i)^{1/3}.
This produces an inverted mass spectrum: seesaw masses proportional to 1/m_quark.

The dual Koide observation: Q(1/m_d, 1/m_s, 1/m_b) = 0.66521 (0.22% from 2/3).
This is the ONLY triple out of 28 scanned whose inverse masses satisfy Q ~ 2/3.

This script investigates:
  1. Numerical seesaw spectrum for down-type quarks
  2. Cancellation of C in the Koide quotient
  3. Algebraic condition for Q(1/m) ~ 2/3 while Q(m) != 2/3
  4. Koide parametrization fits for (e,mu,tau), (d,s,b), and (1/m_d, 1/m_s, 1/m_b)
  5. Relation between delta values across triples
  6. Mass prediction from dual Koide = 2/3 exactly
  7. Compatibility test: dual Koide vs v0-doubling predictions for m_b
"""

import numpy as np
from scipy.optimize import brentq, minimize_scalar

# ============================================================================
# Physical masses (MeV), PDG values
# ============================================================================

m_e   = 0.51099895
m_mu  = 105.6583755
m_tau = 1776.86

m_u   = 2.16
m_d   = 4.67
m_s   = 93.4
m_c   = 1270.0
m_b   = 4180.0
m_t   = 172760.0

SQRT2         = np.sqrt(2.0)
TWO_PI_OVER_3 = 2.0 * np.pi / 3.0


# ============================================================================
# Core functions
# ============================================================================

def koide_Q(masses, signs=None):
    """
    Compute Q = sum(m_j) / (sum(s_j * sqrt(m_j)))^2.
    signs: list of +1/-1 for each mass. Default: all +1.
    """
    ms = np.array(masses, dtype=float)
    if signs is None:
        signs = np.ones(len(ms))
    ss = np.array(signs, dtype=float)
    numerator = np.sum(ms)
    denominator = np.sum(ss * np.sqrt(ms))**2
    return numerator / denominator


def f_vec(delta):
    """Return f = (f0, f1, f2) where f_k = 1 + sqrt(2)*cos(2*pi*k/3 + delta)."""
    return np.array([
        1 + SQRT2 * np.cos(delta),
        1 + SQRT2 * np.cos(delta + TWO_PI_OVER_3),
        1 + SQRT2 * np.cos(delta + 2 * TWO_PI_OVER_3),
    ])


def fit_koide_parametrization(masses, signs=None):
    """
    Fit Koide parametrization sqrt(m_k) = sqrt(M0) * (1 + sqrt(2)*cos(2pi*k/3 + delta))
    to given masses.

    Returns (M0, delta) by analytical extraction from signed roots.
    For Q ~ 2/3 triples this is essentially exact.

    When signs are given, signed_roots = signs * sqrt(masses).
    The parametrization assigns slot 0 to the lightest/most-negative signed root.
    """
    ms = np.array(masses, dtype=float)
    if signs is None:
        signs = np.ones(len(ms))
    ss = np.array(signs, dtype=float)

    # Signed square roots
    sigma = ss * np.sqrt(ms)

    # Sort sigma to assign to slots: smallest signed root -> slot 0
    idx = np.argsort(sigma)
    sigma_sorted = sigma[idx]

    # Analytical extraction
    # sum sigma_k = 3 * sqrt(M0)
    sqrtM0 = sigma_sorted.sum() / 3.0
    M0 = sqrtM0**2

    # delta from the two Fourier components
    phi = np.array([TWO_PI_OVER_3 * k for k in range(3)])
    A = np.dot(sigma_sorted, np.cos(phi))
    B = np.dot(sigma_sorted, -np.sin(phi))
    # A = sqrt(2) * sqrt(M0) * cos(delta), B = sqrt(2) * sqrt(M0) * sin(delta)
    delta = np.arctan2(B, A) % (2.0 * np.pi)

    return M0, delta, idx


def koide_masses_from_params(M0, delta):
    """Return three masses from the Koide parametrization."""
    f = f_vec(delta)
    return M0 * f**2


def all_sv_positive(delta):
    """Check if all signed roots are positive at given delta."""
    for k in range(3):
        sv = 1 + SQRT2 * np.cos(TWO_PI_OVER_3 * k + delta)
        if sv <= 0:
            return False
    return True


# ============================================================================
# Output accumulator
# ============================================================================

output_lines = []

def p(s=""):
    output_lines.append(s)
    print(s)


# ############################################################################
# PART 1: Down-type seesaw spectrum
# ############################################################################

p("=" * 76)
p("PART 1: Down-type seesaw spectrum")
p("=" * 76)
p()
p("For N_f = 3 SQCD with flavors (d, s, b):")
p("  M_j = C'/m_j  where C' = Lambda^2 * (m_d * m_s * m_b)^{1/3}")
p()

Lambda = 300.0  # MeV, representative SQCD scale

C_prime_cubed = Lambda**6 * m_d * m_s * m_b
C_prime = C_prime_cubed**(1.0/3.0)

p(f"  Lambda = {Lambda:.0f} MeV")
p(f"  C'^3 = Lambda^6 * m_d * m_s * m_b = {C_prime_cubed:.6e} MeV^3")
p(f"  C'   = {C_prime:.6f} MeV")
p()

M_D = C_prime / m_d
M_S = C_prime / m_s
M_B = C_prime / m_b

p("  Seesaw masses:")
p(f"    M_d = C'/m_d = {M_D:.4f} MeV")
p(f"    M_s = C'/m_s = {M_S:.4f} MeV")
p(f"    M_b = C'/m_b = {M_B:.4f} MeV")
p()
p(f"  Mass hierarchy is INVERTED: M_d > M_s > M_b")
p(f"  Ratio M_d/M_b = {M_D/M_B:.2f}  (= m_b/m_d = {m_b/m_d:.2f})")
p()


# ############################################################################
# PART 2: Q on seesaw spectrum = Q on inverse masses
# ############################################################################

p("=" * 76)
p("PART 2: Verify Q(M_d, M_s, M_b) = Q(1/m_d, 1/m_s, 1/m_b)")
p("=" * 76)
p()

Q_seesaw = koide_Q([M_D, M_S, M_B])
Q_inverse = koide_Q([1.0/m_d, 1.0/m_s, 1.0/m_b])
Q_direct_dsb = koide_Q([m_d, m_s, m_b])

p(f"  Q(M_d, M_s, M_b)          = {Q_seesaw:.10f}")
p(f"  Q(1/m_d, 1/m_s, 1/m_b)    = {Q_inverse:.10f}")
p(f"  Difference                 = {abs(Q_seesaw - Q_inverse):.2e}  (machine epsilon)")
p()
p(f"  C' cancels: confirmed.")
p()
p(f"  For comparison:")
p(f"  Q(m_d, m_s, m_b)  [direct] = {Q_direct_dsb:.10f}  (2/3 dev: {(Q_direct_dsb - 2/3)/(2/3)*100:+.4f}%)")
p(f"  Q(1/m_d, 1/m_s, 1/m_b)    = {Q_inverse:.10f}  (2/3 dev: {(Q_inverse - 2/3)/(2/3)*100:+.4f}%)")
p()
p(f"  Direct: Q = 0.731 (9.7% off). Inverse: Q = 0.665 (0.22% off).")
p(f"  The seesaw TRANSFORMS a non-Koide spectrum into a near-Koide one.")
p()


# ############################################################################
# PART 3: Algebraic condition for Q(1/m) ~ 2/3 while Q(m) != 2/3
# ############################################################################

p("=" * 76)
p("PART 3: Algebraic condition — when does Q(1/m) ~ 2/3 while Q(m) != 2/3?")
p("=" * 76)
p()

p("""Let a_j = sqrt(m_j). Define symmetric polynomials:
  p1 = a1 + a2 + a3
  p2 = a1*a2 + a2*a3 + a1*a3
  p3 = a1*a2*a3

Then:
  Q(m)   = (a1^2 + a2^2 + a3^2) / (a1 + a2 + a3)^2
         = 1 - 2*p2/p1^2

  Q(1/m) = (1/a1^2 + 1/a2^2 + 1/a3^2) / (1/a1 + 1/a2 + 1/a3)^2
         = (p2^2 - 2*p1*p3) / p2^2
         = 1 - 2*p1*p3/p2^2

So:
  Q(m)   = 2/3  iff  p2/p1^2 = 1/6
  Q(1/m) = 2/3  iff  p1*p3/p2^2 = 1/6

These are INDEPENDENT conditions on different combinations of symmetric polynomials.
""")

# Verify numerically
a_d = np.sqrt(m_d)
a_s = np.sqrt(m_s)
a_b = np.sqrt(m_b)

p1 = a_d + a_s + a_b
p2 = a_d*a_s + a_s*a_b + a_d*a_b
p3 = a_d*a_s*a_b

ratio_direct = p2 / p1**2
ratio_dual = p1 * p3 / p2**2

p("Numerical verification for (d, s, b):")
p(f"  a_d = {a_d:.6f}, a_s = {a_s:.6f}, a_b = {a_b:.6f}")
p(f"  p1 = {p1:.6f}")
p(f"  p2 = {p2:.6f}")
p(f"  p3 = {p3:.6f}")
p()
p(f"  p2/p1^2    = {ratio_direct:.8f}  (1/6 = {1/6:.8f} for Q = 2/3)")
p(f"  p1*p3/p2^2 = {ratio_dual:.8f}  (1/6 = {1/6:.8f} for Q_dual = 2/3)")
p()
p(f"  Q(m)   = 1 - 2*{ratio_direct:.8f} = {1 - 2*ratio_direct:.8f}  (actual: {Q_direct_dsb:.8f})")
p(f"  Q(1/m) = 1 - 2*{ratio_dual:.8f} = {1 - 2*ratio_dual:.8f}  (actual: {Q_inverse:.8f})")
p()

# What makes the dual condition nearly satisfied?
p("The condition p1*p3/p2^2 = 1/6 can be rewritten:")
p("  6*p1*p3 = p2^2")
p("  6*(a1+a2+a3)*(a1*a2*a3) = (a1*a2 + a2*a3 + a1*a3)^2")
p()
six_p1_p3 = 6 * p1 * p3
p2_sq = p2**2
p(f"  6*p1*p3 = {six_p1_p3:.4f}")
p(f"  p2^2    = {p2_sq:.4f}")
p(f"  Ratio   = {six_p1_p3/p2_sq:.8f}  (1.000 for exact dual Koide)")
p()

# The key asymmetry: large mass hierarchy
p("Physical intuition: with mass hierarchy m_d << m_s << m_b,")
p("  p1 ~ a_b (dominated by heaviest)")
p("  p2 ~ a_s*a_b (dominated by products with heaviest)")
p("  p3 = a_d*a_s*a_b (all three)")
p()
p("  p2/p1^2 ~ (a_s*a_b)/(a_b)^2 = a_s/a_b << 1/6  (too small: Q > 2/3)")
p("  p1*p3/p2^2 ~ (a_b * a_d*a_s*a_b) / (a_s*a_b)^2 = a_d/a_s")
p(f"  Actual a_d/a_s = {a_d/a_s:.4f}  vs 1/6 = {1/6:.4f}")
p()
p("  The leading-order approximation a_d/a_s = 0.224 is NOT close to 1/6 = 0.167.")
p(f"  The deviation: a_d/a_s - 1/6 = {a_d/a_s - 1/6:.5f} ({(a_d/a_s - 1/6)/(1/6)*100:.2f}%)")
p("  The 0.22% closeness of Q_dual to 2/3 comes from the FULL expression")
p("  p1*p3/p2^2 = 0.1674, not from the crude approximation.")
p("  Subleading terms involving all three masses conspire to bring p1*p3/p2^2")
p("  much closer to 1/6 than the leading a_d/a_s term alone.")
p()


# ############################################################################
# PART 4: Parametric analysis
# ############################################################################

p("=" * 76)
p("PART 4: Koide parametrization fits")
p("=" * 76)
p()

# --- 4a: Charged leptons ---

p("--- 4a: Charged leptons (e, mu, tau) ---")
p()

M0_lep, delta_lep, idx_lep = fit_koide_parametrization([m_e, m_mu, m_tau], [+1, +1, +1])
Q_lep = koide_Q([m_e, m_mu, m_tau])
Q_dual_lep = koide_Q([1/m_e, 1/m_mu, 1/m_tau])

# Compute Q_dual via the formula for the Koide manifold
cos3d_lep = np.cos(3 * delta_lep)
Q_dual_formula_lep = 7/3 - 4*SQRT2*cos3d_lep/3

p(f"  M0    = {M0_lep:.6f} MeV")
p(f"  delta = {delta_lep:.8f} rad = {np.degrees(delta_lep):.4f} deg")
p(f"  Q(e,mu,tau)        = {Q_lep:.10f}  (2/3 dev: {(Q_lep - 2/3)/(2/3)*100:+.6f}%)")
p(f"  Q_dual(e,mu,tau)   = {Q_dual_lep:.10f}  (2/3 dev: {(Q_dual_lep - 2/3)/(2/3)*100:+.4f}%)")
p(f"  cos(3*delta)       = {cos3d_lep:.10f}")
p(f"  Q_dual (formula)   = {Q_dual_formula_lep:.10f}  (valid={all_sv_positive(delta_lep)})")
p()

# Reconstruct and verify
m_recon_lep = koide_masses_from_params(M0_lep, delta_lep)
sorted_actual_lep = np.array([m_e, m_mu, m_tau])[idx_lep]
p("  Reconstruction check (sorted by signed root):")
for k in range(3):
    p(f"    slot {k}: model = {m_recon_lep[k]:.6f} MeV, actual = {sorted_actual_lep[k]:.6f} MeV, "
      f"err = {(m_recon_lep[k] - sorted_actual_lep[k])/sorted_actual_lep[k]*100:+.6f}%")
p()


# --- 4b: Down-type quarks (d, s, b) ---

p("--- 4b: Down-type quarks (d, s, b) ---")
p()

M0_dsb, delta_dsb, idx_dsb = fit_koide_parametrization([m_d, m_s, m_b], [+1, +1, +1])
Q_dsb = koide_Q([m_d, m_s, m_b])
Q_dual_dsb_num = koide_Q([1/m_d, 1/m_s, 1/m_b])

cos3d_dsb = np.cos(3 * delta_dsb)
Q_dual_formula_dsb = 7/3 - 4*SQRT2*cos3d_dsb/3

p(f"  M0    = {M0_dsb:.6f} MeV")
p(f"  delta = {delta_dsb:.8f} rad = {np.degrees(delta_dsb):.4f} deg")
p(f"  Q(d,s,b)           = {Q_dsb:.10f}  (2/3 dev: {(Q_dsb - 2/3)/(2/3)*100:+.4f}%)")
p(f"  Q_dual(d,s,b)      = {Q_dual_dsb_num:.10f}  (2/3 dev: {(Q_dual_dsb_num - 2/3)/(2/3)*100:+.4f}%)")
p(f"  cos(3*delta)       = {cos3d_dsb:.10f}")
p(f"  Q_dual (formula)   = {Q_dual_formula_dsb:.10f}  (valid={all_sv_positive(delta_dsb)})")
p()

# Note: since Q(d,s,b) != 2/3, the formula Q_dual = 7/3 - 4sqrt2 cos(3d)/3
# only applies on the Koide manifold and won't match here.
p("  NOTE: Q(d,s,b) = 0.731 != 2/3, so (d,s,b) is NOT on the Koide manifold.")
p("  The formula Q_dual = 7/3 - 4*sqrt(2)*cos(3*delta)/3 does NOT apply.")
p("  The dual Koide observation is a DIRECT numerical fact about the masses,")
p("  not a consequence of the parametrization formula.")
p()


# --- 4c: INVERSE down-type quarks (1/m_d, 1/m_s, 1/m_b) ---

p("--- 4c: Inverse down-type quarks (1/m_d, 1/m_s, 1/m_b) ---")
p()
p("  Since Q(1/m_d, 1/m_s, 1/m_b) ~ 2/3, these inverse masses ARE")
p("  near the Koide manifold. We fit the parametrization to them.")
p()

inv_masses_dsb = [1.0/m_d, 1.0/m_s, 1.0/m_b]

M0_dual, delta_dual, idx_dual = fit_koide_parametrization(inv_masses_dsb, [+1, +1, +1])

p(f"  Inverse masses (1/MeV):")
p(f"    1/m_d = {1/m_d:.8f}")
p(f"    1/m_s = {1/m_s:.8f}")
p(f"    1/m_b = {1/m_b:.8f}")
p()
p(f"  M0_dual    = {M0_dual:.10e} (1/MeV)")
p(f"  delta_dual = {delta_dual:.8f} rad = {np.degrees(delta_dual):.4f} deg")
p()

# Check what delta_dual is close to
p(f"  delta_dual / pi       = {delta_dual/np.pi:.8f}")
p(f"  delta_dual / (pi/6)   = {delta_dual/(np.pi/6):.6f}")
p(f"  delta_dual / (2pi/3)  = {delta_dual/TWO_PI_OVER_3:.6f}")
p()

# Critical value for dual Koide = 2/3: cos(3*delta) = 5*sqrt(2)/8
cos3d_crit = 5*SQRT2/8
delta_crit = np.arccos(cos3d_crit) / 3.0  # principal value
p(f"  For Q_dual = 2/3 on the Koide manifold: cos(3*delta) = 5*sqrt(2)/8 = {cos3d_crit:.8f}")
p(f"  This gives delta = arccos(5*sqrt(2)/8)/3 = {delta_crit:.8f} rad = {np.degrees(delta_crit):.4f} deg")
p()

cos3d_dual = np.cos(3 * delta_dual)
Q_dual_formula_dual = 7/3 - 4*SQRT2*cos3d_dual/3

p(f"  Actual cos(3*delta_dual)     = {cos3d_dual:.10f}")
p(f"  Required cos(3*delta)        = {cos3d_crit:.10f}")
p(f"  Deviation:                     {(cos3d_dual - cos3d_crit)/cos3d_crit*100:+.4f}%")
p()
p(f"  Q_dual formula check = {Q_dual_formula_dual:.10f}  (valid={all_sv_positive(delta_dual)})")
p(f"  Q_dual numerical     = {Q_dual_dsb_num:.10f}")
p()

# Reconstruction of inverse masses
m_recon_dual = koide_masses_from_params(M0_dual, delta_dual)
sorted_inv = np.array(inv_masses_dsb)[idx_dual]
p("  Reconstruction check:")
for k in range(3):
    p(f"    slot {k}: model = {m_recon_dual[k]:.10e}, actual = {sorted_inv[k]:.10e}, "
      f"err = {(m_recon_dual[k] - sorted_inv[k])/sorted_inv[k]*100:+.6f}%")
p()


# ############################################################################
# PART 5: Connecting the two Koide phases
# ############################################################################

p("=" * 76)
p("PART 5: Connecting the two Koide phases — delta relations")
p("=" * 76)
p()

# Fit all the known near-Koide triples
triples_info = {
    '(e,mu,tau) +++':    {'masses': [m_e, m_mu, m_tau], 'signs': [+1, +1, +1]},
    '(-s,c,b) -++':      {'masses': [m_s, m_c, m_b],    'signs': [-1, +1, +1]},
    '(c,b,t) +++':       {'masses': [m_c, m_b, m_t],    'signs': [+1, +1, +1]},
}

p("Delta values for all near-Koide triples:")
p()

delta_values = {}
for name, info in triples_info.items():
    M0_t, delta_t, idx_t = fit_koide_parametrization(info['masses'], info['signs'])
    Q_t = koide_Q(info['masses'], info['signs'])
    delta_values[name] = delta_t
    p(f"  {name:<20}: delta = {delta_t:.8f} rad = {np.degrees(delta_t):8.4f} deg,  "
      f"Q = {Q_t:.8f},  M0 = {M0_t:.4f} MeV")

p()
p(f"  Inverse (d,s,b)     : delta_dual = {delta_dual:.8f} rad = {np.degrees(delta_dual):8.4f} deg")
p()

# Q_dual = 2/3 requires cos(3*delta) = 5*sqrt(2)/8
# For the dsb direct triple, compute delta and check relation to delta_dual
p("Algebraic relations between delta values:")
p()

# The direct (d,s,b) delta
delta_direct_dsb = delta_dsb
p(f"  delta(d,s,b direct)  = {np.degrees(delta_direct_dsb):.4f} deg")
p(f"  delta(d,s,b inverse) = {np.degrees(delta_dual):.4f} deg")
p(f"  delta(e,mu,tau)      = {np.degrees(delta_values['(e,mu,tau) +++']):.4f} deg")
p(f"  delta(-s,c,b)        = {np.degrees(delta_values['(-s,c,b) -++']):.4f} deg")
p(f"  delta(c,b,t)         = {np.degrees(delta_values['(c,b,t) +++']):.4f} deg")
p()

# Check simple relations
delta_scb = delta_values['(-s,c,b) -++']
delta_cbt = delta_values['(c,b,t) +++']
delta_emt = delta_values['(e,mu,tau) +++']

p("Simple ratios:")
p(f"  delta_dual / delta(e,mu,tau)  = {delta_dual / delta_emt:.6f}")
p(f"  delta_dual / delta(-s,c,b)   = {delta_dual / delta_scb:.6f}")
p(f"  delta_dual / delta(c,b,t)    = {delta_dual / delta_cbt:.6f}")
p()

# Check sums/differences
p("Differences (mod 2pi/3):")
for name, dv in delta_values.items():
    diff = (delta_dual - dv) % TWO_PI_OVER_3
    p(f"  delta_dual - delta{name} = {np.degrees(diff):.4f} deg (mod 120)")

p()

# Check whether delta_dual has a special value
p("Is delta_dual close to a simple fraction of pi?")
for num in range(1, 25):
    for den in range(1, 25):
        test = np.pi * num / den
        if abs(test - delta_dual) < 0.01 or abs(test - delta_dual + 2*np.pi) < 0.01:
            p(f"  delta_dual ~ {num}/{den} * pi = {np.degrees(test):.4f} deg  "
              f"(actual: {np.degrees(delta_dual):.4f} deg, "
              f"diff: {np.degrees(abs(test - delta_dual)):.4f} deg)")

p()

# The special value: delta for Q_dual = 2/3 on Koide manifold
p("The critical delta where Q_dual = 2/3 (on the Koide manifold):")
p(f"  cos(3*delta_crit) = 5*sqrt(2)/8 = {cos3d_crit:.8f}")

# All solutions in [0, 2pi)
all_crits = []
base = np.arccos(cos3d_crit)
for n in range(-3, 4):
    for sign in [1, -1]:
        sol = (sign * base + 2*np.pi*n) / 3
        sol_mod = sol % (2*np.pi)
        if 0 <= sol_mod < 2*np.pi:
            all_crits.append(sol_mod)
all_crits = sorted(set([round(s, 10) for s in all_crits]))

p(f"  All critical delta in [0, 2pi):")
for dc in all_crits:
    valid = all_sv_positive(dc)
    p(f"    delta = {dc:.8f} rad = {np.degrees(dc):.4f} deg  (all sv>0: {valid})")

p()

# Closest critical delta to delta_dual
dists = [(abs(delta_dual - dc), dc) for dc in all_crits]
dists.sort()
closest_crit = dists[0][1]
p(f"  Closest to delta_dual ({np.degrees(delta_dual):.4f} deg):")
p(f"    delta_crit = {np.degrees(closest_crit):.4f} deg, difference = {np.degrees(abs(delta_dual - closest_crit)):.4f} deg")
p()


# ############################################################################
# PART 6: Mass prediction from dual Koide = 2/3 exactly
# ############################################################################

p("=" * 76)
p("PART 6: Mass prediction from exact dual Koide")
p("=" * 76)
p()

p("If Q(1/m_d, 1/m_s, 1/m_b) = 2/3 exactly, this is equivalent to:")
p()
p("  3*(1/m_d + 1/m_s + 1/m_b) = 2*(1/sqrt(m_d) + 1/sqrt(m_s) + 1/sqrt(m_b))^2")
p()
p("Given m_d and m_s, solve for m_b.")
p()

def dual_koide_residual(mb, md, ms):
    """Q(1/md, 1/ms, 1/mb) - 2/3"""
    inv_m = [1.0/md, 1.0/ms, 1.0/mb]
    return koide_Q(inv_m) - 2.0/3.0


# Solve numerically
# First check the sign at boundaries
p(f"  Inputs: m_d = {m_d} MeV, m_s = {m_s} MeV")
p()

# Scan for the root
mb_scan = np.linspace(100, 50000, 500000)
residuals = np.array([dual_koide_residual(mb, m_d, m_s) for mb in mb_scan])

# Find sign changes
sign_changes = []
for i in range(len(residuals)-1):
    if np.isfinite(residuals[i]) and np.isfinite(residuals[i+1]):
        if residuals[i] * residuals[i+1] < 0:
            sign_changes.append(i)

p(f"  Found {len(sign_changes)} sign change(s) in residual over [100, 50000] MeV")

mb_predictions_dual = []
for i in sign_changes:
    mb_sol = brentq(dual_koide_residual, mb_scan[i], mb_scan[i+1],
                    args=(m_d, m_s), xtol=1e-10)
    Q_check = koide_Q([1/m_d, 1/m_s, 1/mb_sol])
    mb_predictions_dual.append(mb_sol)
    dev_pdg = (mb_sol - m_b) / m_b * 100
    p(f"    m_b = {mb_sol:.4f} MeV  (Q_check = {Q_check:.12f})")
    p(f"    PDG m_b = {m_b:.0f} MeV,  deviation = {dev_pdg:+.4f}%")
    p(f"    v0-doubling prediction m_b = 4177 MeV,  deviation from dual = {(mb_sol - 4177)/4177*100:+.4f}%")

p()

# Also solve analytically: Q(1/m) = 2/3 means
# 3*(1/m_d + 1/m_s + 1/m_b) = 2*(1/sqrt(m_d) + 1/sqrt(m_s) + 1/sqrt(m_b))^2
# Let x = 1/sqrt(m_b)
# 3*(1/m_d + 1/m_s + x^2) = 2*(1/sqrt(m_d) + 1/sqrt(m_s) + x)^2
# Let A = 1/sqrt(m_d) + 1/sqrt(m_s), B = 1/m_d + 1/m_s
# 3*B + 3*x^2 = 2*(A + x)^2 = 2*A^2 + 4*A*x + 2*x^2
# x^2 - 4*A*x + 3*B - 2*A^2 = 0

p("Analytical solution:")
p()
A = 1/np.sqrt(m_d) + 1/np.sqrt(m_s)
B = 1/m_d + 1/m_s
disc = 16*A**2 - 4*(3*B - 2*A**2)
# = 16*A^2 - 12*B + 8*A^2 = 24*A^2 - 12*B

p(f"  A = 1/sqrt(m_d) + 1/sqrt(m_s) = {A:.10f}")
p(f"  B = 1/m_d + 1/m_s = {B:.10f}")
p(f"  Discriminant = 24*A^2 - 12*B = {disc:.10f}")
p()

if disc >= 0:
    sq_disc = np.sqrt(disc)
    x1 = (4*A + sq_disc) / 2
    x2 = (4*A - sq_disc) / 2
    p(f"  x = 1/sqrt(m_b):")
    p(f"    x1 = {x1:.10f}  -> m_b = {1/x1**2:.4f} MeV")
    p(f"    x2 = {x2:.10f}  -> m_b = {1/x2**2:.4f} MeV")
    p()

    for xi, label in [(x1, "x1"), (x2, "x2")]:
        if xi > 0:
            mb_anal = 1/xi**2
            Q_check = koide_Q([1/m_d, 1/m_s, 1/mb_anal])
            p(f"  Solution {label}: m_b = {mb_anal:.4f} MeV")
            p(f"    Q(1/m_d, 1/m_s, 1/m_b) = {Q_check:.12f}")
            p(f"    PDG deviation: {(mb_anal - m_b)/m_b*100:+.4f}%")
            p(f"    Absolute difference: {mb_anal - m_b:+.2f} MeV")

p()

# Summary of m_b predictions
p("Summary of m_b predictions:")
p(f"  PDG:            m_b = {m_b:.0f} +/- 30 MeV")
p(f"  v0-doubling:    m_b = 4177 MeV  (0.07%, 0.1 sigma)")
if mb_predictions_dual:
    mb_pred = mb_predictions_dual[0]  # take the physical solution
    p(f"  Dual Koide:     m_b = {mb_pred:.1f} MeV  ({(mb_pred - m_b)/m_b*100:+.3f}%)")
    p(f"  Overlap:        m_b = 4159 MeV  (0.5%)")
p()


# ############################################################################
# PART 7: Compatibility test — dual Koide vs v0-doubling
# ############################################################################

p("=" * 76)
p("PART 7: Compatibility test — dual Koide vs v0-doubling")
p("=" * 76)
p()

p(f"Inputs: m_s = {m_s} MeV, m_c = {m_c} MeV, m_d = {m_d} MeV")
p()

# v0-doubling prediction
sq_s = np.sqrt(m_s)
sq_c = np.sqrt(m_c)
mb_v0 = (3*sq_s + sq_c)**2
p(f"v0-doubling: sqrt(m_b) = 3*sqrt(m_s) + sqrt(m_c)")
p(f"  = 3 * {sq_s:.6f} + {sq_c:.6f} = {3*sq_s + sq_c:.6f}")
p(f"  m_b(v0) = {mb_v0:.4f} MeV")
p()

# Dual Koide prediction (using m_d)
# Already computed above; use the physical root
# Need to identify which root is the physical one (near 4180)
physical_mb_dual = None
for mb_sol in mb_predictions_dual:
    if abs(mb_sol - m_b) < 1000:
        physical_mb_dual = mb_sol
        break

if physical_mb_dual is None and disc >= 0:
    # Use analytical solutions
    for xi in [x1, x2]:
        if xi > 0:
            mb_test = 1/xi**2
            if abs(mb_test - m_b) < 1000:
                physical_mb_dual = mb_test
                break

if physical_mb_dual is not None:
    p(f"Dual Koide: Q(1/m_d, 1/m_s, 1/m_b) = 2/3 exactly")
    p(f"  m_b(dual) = {physical_mb_dual:.4f} MeV")
    p()
    p(f"Comparison:")
    p(f"  m_b(v0-doubling) = {mb_v0:.4f} MeV")
    p(f"  m_b(dual Koide)  = {physical_mb_dual:.4f} MeV")
    p(f"  Difference:        {mb_v0 - physical_mb_dual:+.4f} MeV ({(mb_v0 - physical_mb_dual)/physical_mb_dual*100:+.4f}%)")
    p(f"  PDG m_b:           {m_b:.0f} +/- 30 MeV")
    p()
    p(f"  v0-doubling deviation from PDG: {(mb_v0 - m_b)/m_b*100:+.4f}% = {(mb_v0 - m_b)/30:+.2f} sigma")
    p(f"  Dual Koide deviation from PDG:  {(physical_mb_dual - m_b)/m_b*100:+.4f}% = {(physical_mb_dual - m_b)/30:+.2f} sigma")
    p()

    # Are they simultaneously satisfiable?
    p("Can both conditions hold simultaneously?")
    p()
    p("  v0-doubling condition: sqrt(m_b) = 3*sqrt(m_s) + sqrt(m_c)")
    p("  Dual Koide condition:  Q(1/m_d, 1/m_s, 1/m_b) = 2/3")
    p()
    p("  These involve different input sets:")
    p("    v0-doubling uses (m_s, m_c)")
    p("    Dual Koide uses (m_d, m_s)")
    p()
    p("  They predict different m_b values, so they are NOT exactly compatible")
    p(f"  with the current PDG masses. The tension is {abs(mb_v0 - physical_mb_dual):.2f} MeV")
    p(f"  ({abs(mb_v0 - physical_mb_dual)/30:.2f} sigma given PDG uncertainty).")
    p()

    # What m_s would make them compatible?
    p("  If we demand BOTH conditions hold with m_d = 4.67, m_c = 1270:")
    p("  Then m_b = (3*sqrt(m_s) + sqrt(m_c))^2 AND Q(1/m_d, 1/m_s, 1/m_b) = 2/3")
    p("  This constrains m_s.")
    p()

    def combined_residual(ms_test):
        sq_s_t = np.sqrt(ms_test)
        sq_c_t = np.sqrt(m_c)
        mb_test = (3*sq_s_t + sq_c_t)**2
        return koide_Q([1/m_d, 1/ms_test, 1/mb_test]) - 2/3

    # Scan
    ms_scan = np.linspace(50, 150, 100000)
    resid_ms = np.array([combined_residual(ms) for ms in ms_scan])

    sign_changes_ms = []
    for i in range(len(resid_ms)-1):
        if np.isfinite(resid_ms[i]) and np.isfinite(resid_ms[i+1]):
            if resid_ms[i] * resid_ms[i+1] < 0:
                sign_changes_ms.append(i)

    if sign_changes_ms:
        for i in sign_changes_ms:
            ms_sol = brentq(combined_residual, ms_scan[i], ms_scan[i+1], xtol=1e-12)
            mb_compat = (3*np.sqrt(ms_sol) + np.sqrt(m_c))**2
            Q_check = koide_Q([1/m_d, 1/ms_sol, 1/mb_compat])
            p(f"  Solution: m_s = {ms_sol:.4f} MeV  (PDG: {m_s} +/- 8 MeV)")
            p(f"            m_b = {mb_compat:.4f} MeV  (PDG: {m_b} +/- 30 MeV)")
            p(f"            Q_check = {Q_check:.12f}")
            p(f"            m_s deviation from PDG: {(ms_sol - m_s)/m_s*100:+.3f}%")
            p(f"            m_b deviation from PDG: {(mb_compat - m_b)/m_b*100:+.3f}%")
    else:
        p("  No simultaneous solution found in the scanned range.")

p()

# Also check: what if we use m_d as free parameter instead?
p("  Alternatively, fix m_s = 93.4 and m_c = 1270, and find m_d that makes")
p("  dual Koide exact with m_b from v0-doubling:")
p()

mb_v0_fixed = mb_v0  # m_b from v0-doubling with m_s=93.4, m_c=1270

def md_residual(md_test):
    return koide_Q([1/md_test, 1/m_s, 1/mb_v0_fixed]) - 2/3

md_scan = np.linspace(1, 20, 100000)
resid_md = np.array([md_residual(md) for md in md_scan])

sign_changes_md = []
for i in range(len(resid_md)-1):
    if np.isfinite(resid_md[i]) and np.isfinite(resid_md[i+1]):
        if resid_md[i] * resid_md[i+1] < 0:
            sign_changes_md.append(i)

if sign_changes_md:
    for i in sign_changes_md:
        md_sol = brentq(md_residual, md_scan[i], md_scan[i+1], xtol=1e-12)
        Q_check = koide_Q([1/md_sol, 1/m_s, 1/mb_v0_fixed])
        p(f"  Solution: m_d = {md_sol:.4f} MeV  (PDG: {m_d} +/- 0.48 MeV)")
        p(f"            m_b = {mb_v0_fixed:.4f} MeV  (from v0-doubling)")
        p(f"            Q_check = {Q_check:.12f}")
        p(f"            m_d deviation from PDG: {(md_sol - m_d)/m_d*100:+.3f}%")
        p(f"            m_d deviation in sigma: {(md_sol - m_d)/0.48:+.2f}")


# ############################################################################
# FINAL SUMMARY
# ############################################################################

p()
p()
p("=" * 76)
p("FINAL SUMMARY")
p("=" * 76)
p()

p("1. SEESAW SPECTRUM: The Seiberg seesaw M_j = C'/m_j inverts the down-type")
p("   quark hierarchy. The scale C' depends on Lambda but cancels in Q.")
p()

p("2. DUAL KOIDE: Q(1/m_d, 1/m_s, 1/m_b) = {:.8f}, which is {:.4f}% from 2/3.".format(
    Q_inverse, (Q_inverse - 2/3)/(2/3)*100))
p("   The seesaw maps a non-Koide spectrum (Q = 0.731) to a near-Koide one.")
p()

p("3. ALGEBRAIC ORIGIN: Q_dual ~ 2/3 requires p1*p3/p2^2 ~ 1/6.")
p(f"   Actual: p1*p3/p2^2 = {ratio_dual:.8f}, 1/6 = {1/6:.8f} ({(ratio_dual - 1/6)/(1/6)*100:+.3f}%).")
p("   The leading-order approximation (a_d/a_s) is 34% off from 1/6;")
p("   subleading terms from all three masses bring the full ratio to within 0.44%.")
p()

p("4. PARAMETRIZATION: The inverse masses fit the Koide parametrization with")
p(f"   delta_dual = {np.degrees(delta_dual):.4f} deg.")
if physical_mb_dual:
    p(f"   The dual Koide condition Q = 2/3 predicts m_b = {physical_mb_dual:.1f} MeV.")
p()

p("5. MASS PREDICTIONS for m_b:")
p(f"   PDG:            {m_b:.0f} +/- 30 MeV")
p(f"   v0-doubling:    {mb_v0:.1f} MeV  ({(mb_v0 - m_b)/m_b*100:+.3f}%)")
if physical_mb_dual:
    p(f"   Dual Koide:     {physical_mb_dual:.1f} MeV  ({(physical_mb_dual - m_b)/m_b*100:+.3f}%)")
p(f"   Overlap:        4159 MeV  (-0.502%)")
p()

if physical_mb_dual:
    p("6. COMPATIBILITY: The dual Koide condition Q(1/m_d, 1/m_s, 1/m_b) = 2/3")
    p(f"   predicts m_b = {physical_mb_dual:.0f} MeV, which is {(physical_mb_dual - m_b)/30:.1f} sigma from PDG.")
    p("   This is a POOR prediction. The v0-doubling prediction (m_b = 4177) is far superior.")
    p()
    p("   Key insight: Q_dual = 0.665 being close to 2/3 does not mean the exact")
    p("   condition Q_dual = 2/3 predicts m_b well. The Koide equation is a quadratic")
    p("   in 1/sqrt(m_b), and the 0.22% deviation in Q translates to a 9% deviation in m_b.")
    p()
    p("   However, both conditions CAN be made simultaneously compatible by adjusting")
    p("   m_d within 0.14 sigma of PDG (m_d = 4.604 vs 4.67 +/- 0.48).")

p()
p("Done.")


# ============================================================================
# Write markdown summary
# ============================================================================

md_lines = []

md_lines.append("# Seesaw Mechanism and Dual Koide Analysis")
md_lines.append("")
md_lines.append("## Setup")
md_lines.append("")
md_lines.append("In SU(3) SQCD with N_f = N_c = 3, the Seiberg duality maps quarks to")
md_lines.append("mesons with VEVs M_j = C'/m_j, where C' = Lambda^2 (m_d m_s m_b)^{1/3}.")
md_lines.append("This produces an inverted mass spectrum.")
md_lines.append("")
md_lines.append("## 1. Seesaw Spectrum")
md_lines.append("")
md_lines.append(f"With Lambda = {Lambda:.0f} MeV:")
md_lines.append(f"- C' = {C_prime:.4f} MeV")
md_lines.append(f"- M_d = C'/m_d = {M_D:.2f} MeV")
md_lines.append(f"- M_s = C'/m_s = {M_S:.2f} MeV")
md_lines.append(f"- M_b = C'/m_b = {M_B:.2f} MeV")
md_lines.append("")
md_lines.append("The hierarchy is inverted: M_d > M_s > M_b.")
md_lines.append("")
md_lines.append("## 2. C' Cancellation")
md_lines.append("")
md_lines.append("Q(M_d, M_s, M_b) = Q(1/m_d, 1/m_s, 1/m_b) exactly (C' cancels).")
md_lines.append(f"- Q(1/m_d, 1/m_s, 1/m_b) = {Q_inverse:.8f}  (2/3 deviation: {(Q_inverse - 2/3)/(2/3)*100:+.4f}%)")
md_lines.append(f"- Q(m_d, m_s, m_b) = {Q_direct_dsb:.8f}  (2/3 deviation: {(Q_direct_dsb - 2/3)/(2/3)*100:+.4f}%)")
md_lines.append("")
md_lines.append("The seesaw transforms a non-Koide spectrum (Q = 0.731) into a near-Koide one (Q = 0.665).")
md_lines.append("")
md_lines.append("## 3. Algebraic Condition")
md_lines.append("")
md_lines.append("With symmetric polynomials of a_j = sqrt(m_j):")
md_lines.append("- Q(m) = 1 - 2 p2/p1^2")
md_lines.append("- Q(1/m) = 1 - 2 p1 p3/p2^2")
md_lines.append("")
md_lines.append("These are independent conditions. Q_dual ~ 2/3 requires p1 p3/p2^2 ~ 1/6.")
md_lines.append(f"Actual: p1 p3/p2^2 = {ratio_dual:.8f} vs 1/6 = {1/6:.8f} ({(ratio_dual - 1/6)/(1/6)*100:+.3f}%).")
md_lines.append(f"The leading-order approximation a_d/a_s = {a_d/a_s:.4f} is 34% from 1/6;")
md_lines.append(f"subleading terms bring the full ratio to within 0.44%.")
md_lines.append("")
md_lines.append("## 4. Parametric Analysis")
md_lines.append("")
md_lines.append("### Charged leptons (e, mu, tau)")
md_lines.append(f"- M0 = {M0_lep:.6f} MeV, delta = {np.degrees(delta_lep):.4f} deg")
md_lines.append(f"- Q = {Q_lep:.10f}, Q_dual = {Q_dual_lep:.8f}")
md_lines.append("")
md_lines.append("### Down-type quarks (d, s, b) [direct]")
md_lines.append(f"- M0 = {M0_dsb:.6f} MeV, delta = {np.degrees(delta_dsb):.4f} deg")
md_lines.append(f"- Q = {Q_dsb:.8f} (NOT on Koide manifold)")
md_lines.append(f"- Q_dual = {Q_dual_dsb_num:.8f}")
md_lines.append("")
md_lines.append("### Inverse down-type quarks (1/m_d, 1/m_s, 1/m_b)")
md_lines.append(f"- M0_dual = {M0_dual:.10e} (1/MeV), delta_dual = {np.degrees(delta_dual):.4f} deg")
md_lines.append(f"- This IS on the Koide manifold (Q = {Q_inverse:.8f} ~ 2/3)")
md_lines.append(f"- cos(3 delta_dual) = {cos3d_dual:.8f}")
md_lines.append(f"- Required for exact Q = 2/3: cos(3 delta) = {cos3d_crit:.8f}")
md_lines.append("")
md_lines.append("## 5. Delta Relations")
md_lines.append("")
md_lines.append("| Triple | delta (deg) |")
md_lines.append("|--------|------------|")
md_lines.append(f"| (e, mu, tau) | {np.degrees(delta_emt):.4f} |")
md_lines.append(f"| (-s, c, b) | {np.degrees(delta_scb):.4f} |")
md_lines.append(f"| (c, b, t) | {np.degrees(delta_cbt):.4f} |")
md_lines.append(f"| (1/m_d, 1/m_s, 1/m_b) | {np.degrees(delta_dual):.4f} |")
md_lines.append("")
md_lines.append(f"No simple algebraic relation between delta_dual and the other")
md_lines.append(f"deltas has been identified. The dual Koide is an independent")
md_lines.append(f"numerical observation.")
md_lines.append("")
md_lines.append("## 6. Mass Prediction from Exact Dual Koide")
md_lines.append("")
md_lines.append("Setting Q(1/m_d, 1/m_s, 1/m_b) = 2/3 exactly with m_d = 4.67, m_s = 93.4:")
md_lines.append("")

if physical_mb_dual:
    dev_pct = (physical_mb_dual - m_b)/m_b*100
    dev_sig = (physical_mb_dual - m_b)/30
    md_lines.append(f"- Predicted m_b = {physical_mb_dual:.1f} MeV")
    md_lines.append(f"- PDG m_b = {m_b:.0f} +/- 30 MeV")
    md_lines.append(f"- Deviation: {dev_pct:+.3f}% ({dev_sig:+.2f} sigma)")

md_lines.append("")
md_lines.append("### Comparison of m_b predictions")
md_lines.append("")
md_lines.append("| Method | m_b (MeV) | PDG deviation |")
md_lines.append("|--------|----------|--------------|")
md_lines.append(f"| PDG central | {m_b:.0f} | -- |")
md_lines.append(f"| v0-doubling | {mb_v0:.1f} | {(mb_v0 - m_b)/m_b*100:+.3f}% |")
if physical_mb_dual:
    md_lines.append(f"| Dual Koide | {physical_mb_dual:.1f} | {(physical_mb_dual - m_b)/m_b*100:+.3f}% |")
md_lines.append(f"| Overlap | 4159 | -0.502% |")
md_lines.append("")

md_lines.append("## 7. Compatibility Test")
md_lines.append("")
md_lines.append("v0-doubling condition: sqrt(m_b) = 3 sqrt(m_s) + sqrt(m_c)")
md_lines.append("Dual Koide condition: Q(1/m_d, 1/m_s, 1/m_b) = 2/3")
md_lines.append("")
md_lines.append(f"v0-doubling gives m_b = {mb_v0:.1f} MeV (uses m_s, m_c)")
if physical_mb_dual:
    md_lines.append(f"Dual Koide gives m_b = {physical_mb_dual:.1f} MeV (uses m_d, m_s)")
    md_lines.append(f"Difference: {abs(mb_v0 - physical_mb_dual):.1f} MeV ({abs(mb_v0 - physical_mb_dual)/30:.2f} sigma)")
    md_lines.append("")
    md_lines.append("")
    md_lines.append("The exact dual Koide condition predicts m_b = 4562 MeV (9% off PDG),")
    md_lines.append("while v0-doubling predicts 4177 MeV (0.07% off). The 0.22% closeness")
    md_lines.append("of Q_dual to 2/3 does NOT translate to a precise m_b prediction because")
    md_lines.append("the Koide equation amplifies small Q deviations into large mass deviations.")
    md_lines.append("")
    md_lines.append("However, both conditions CAN be made simultaneously compatible by")
    md_lines.append("adjusting m_d to 4.604 MeV (0.14 sigma from PDG).")

md_lines.append("")
md_lines.append("## Key Takeaways")
md_lines.append("")
md_lines.append("1. The Seiberg seesaw maps the down-type quark spectrum to a near-Koide")
md_lines.append("   spectrum. This is not a consequence of Q(m) ~ 2/3 (which fails for d,s,b)")
md_lines.append("   but an independent algebraic property of the mass ratios.")
md_lines.append("")
md_lines.append("2. The condition p1*p3/p2^2 = 1/6 holds to 0.44%. The leading-order")
md_lines.append("   approximation (sqrt(m_d)/sqrt(m_s) ~ 1/6) is 34% off;")
md_lines.append("   the sub-percent agreement requires all three masses.")
md_lines.append("")
md_lines.append("3. IMPORTANT: The exact dual Koide condition Q = 2/3 predicts m_b = 4562 MeV")
md_lines.append("   (9% off PDG). The 0.22% closeness of Q to 2/3 does NOT translate to a")
md_lines.append("   precise mass prediction, because the Koide quadratic amplifies small Q")
md_lines.append("   deviations. The v0-doubling prediction (m_b = 4177, 0.07% off) is far superior.")
md_lines.append("")
md_lines.append("4. Both conditions can be made simultaneously compatible with a small shift")
md_lines.append("   in m_d (0.14 sigma from PDG), suggesting they are not in fundamental tension.")
md_lines.append("")
md_lines.append("5. The ISS seesaw mechanism provides a natural physical framework")
md_lines.append("   for the dual Koide: if SQCD confines with N_f = 3 down-type flavors,")
md_lines.append("   the meson VEVs automatically produce the inverted spectrum.")

md_text = "\n".join(md_lines) + "\n"
md_path = "/home/codexssh/phys3/results/seesaw_mechanism.md"
with open(md_path, "w") as f:
    f.write(md_text)

print(f"\n[Markdown summary written to {md_path}]")
