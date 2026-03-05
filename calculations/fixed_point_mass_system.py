"""
Coupled algebraic mass system: exact solution and numerical analysis.

Four relations among positive reals (m0, m1, m2, m3, m4, m5):
  1. m2 = R * (m0 + m1)         [ratio relation]
  2. m3 = R * m2                 [chain relation]
  3. sqrt(m4) = 3*sqrt(m2) + sqrt(m3)  [doubling relation]
  4. Q(1/sqrt(m1), 1/sqrt(m2), 1/sqrt(m4)) = 2/3  [energy balance on inverses]

where R = (2 + sqrt(3))^2 = 7 + 4*sqrt(3).

Part (d): Q(sqrt(m2), sqrt(m3), sqrt(m5)) = 2/3 determines m5.
"""

import numpy as np
from scipy.optimize import fsolve, brentq

# ─── Constants ───────────────────────────────────────────────────────────
R = (2 + np.sqrt(3))**2
print(f"R = (2 + sqrt(3))^2 = {R:.10f}")
print(f"  = 7 + 4*sqrt(3)  = {7 + 4*np.sqrt(3):.10f}")
print()

# ─── Experimental targets (MeV) ─────────────────────────────────────────
targets = {
    'm0': (2.16, 0.49),
    'm1': (4.67, 0.48),
    'm2': (93.4, 4.8),
    'm3': (1270.0, 25.0),
    'm4': (4183.0, 30.0),
    'm5': (162500.0, None),  # target for part (d)
}

# ─── Koide Q function ───────────────────────────────────────────────────
def Q(x, y, z):
    """Q = (x^2 + y^2 + z^2) / (x + y + z)^2"""
    return (x**2 + y**2 + z**2) / (x + y + z)**2


# =====================================================================
# PART (a): Express all masses as functions of m0
# =====================================================================
print("=" * 70)
print("PART (a): All masses as functions of m0")
print("=" * 70)
print()

print("From equations 1-3, we need m1 to express m2, m3, m4.")
print("Equation 4 (energy balance) constrains the system further.")
print()
print("Strategy: Use eq. 1 to write m2 = R*(m0+m1), eq. 2 gives m3 = R^2*(m0+m1),")
print("eq. 3 gives sqrt(m4) = 3*sqrt(R*(m0+m1)) + sqrt(R^2*(m0+m1)).")
print()
print("Let S = m0 + m1. Then:")
print("  m2 = R*S")
print("  m3 = R^2*S")
print("  sqrt(m4) = sqrt(S) * (3*sqrt(R) + R)")
print("  m4 = S * (3*sqrt(R) + R)^2")
print()

# Compute (3*sqrt(R) + R)^2
coeff_m4 = (3*np.sqrt(R) + R)**2
print(f"  (3*sqrt(R) + R)^2 = {coeff_m4:.10f}")
print(f"  3*sqrt(R) = {3*np.sqrt(R):.10f}")
print(f"  3*sqrt(R) + R = {3*np.sqrt(R) + R:.10f}")
print()

print("So all masses depend on S = m0 + m1:")
print("  m1 = S - m0")
print("  m2 = R * S")
print("  m3 = R^2 * S")
print("  m4 = (3*sqrt(R) + R)^2 * S")
print()

print("Now check whether eq. 4 is automatic or constrains S (equivalently m1).")
print()
print("Substituting into Q(1/sqrt(m1), 1/sqrt(m2), 1/sqrt(m4)) = 2/3:")
print("  x = 1/sqrt(S - m0)")
print("  y = 1/sqrt(R*S)")
print("  z = 1/sqrt((3*sqrt(R)+R)^2 * S)")
print()
print("This depends on BOTH m0 and S, so eq. 4 provides a nontrivial")
print("relation between m0 and S (i.e., between m0 and m1).")
print("The system has exactly ONE free parameter (m0).")
print()

# ─── Solve for m1 given m0 ──────────────────────────────────────────────
def solve_system(m0):
    """Given m0, solve eq. 4 for m1, then compute all masses."""
    coeff4 = (3*np.sqrt(R) + R)**2

    def equation4(m1):
        if m1 <= 0:
            return 1e10
        S = m0 + m1
        m2 = R * S
        m4 = coeff4 * S
        x = 1.0 / np.sqrt(m1)
        y = 1.0 / np.sqrt(m2)
        z = 1.0 / np.sqrt(m4)
        return Q(x, y, z) - 2.0/3.0

    # Bracket search
    m1_sol = brentq(equation4, 0.01, 100.0)
    S = m0 + m1_sol
    m2 = R * S
    m3 = R**2 * S
    m4 = coeff4 * S
    return m1_sol, m2, m3, m4

# Test at target m0
m0_test = 2.16
m1_sol, m2_sol, m3_sol, m4_sol = solve_system(m0_test)
print(f"At m0 = {m0_test}:")
print(f"  m1 = {m1_sol:.6f}  (target: {targets['m1'][0]} +/- {targets['m1'][1]})")
print(f"  m2 = {m2_sol:.4f}  (target: {targets['m2'][0]} +/- {targets['m2'][1]})")
print(f"  m3 = {m3_sol:.2f}  (target: {targets['m3'][0]} +/- {targets['m3'][1]})")
print(f"  m4 = {m4_sol:.2f}  (target: {targets['m4'][0]} +/- {targets['m4'][1]})")
print()

# Verify Q = 2/3
x = 1.0/np.sqrt(m1_sol)
y = 1.0/np.sqrt(m2_sol)
z = 1.0/np.sqrt(m4_sol)
Q_check = Q(x, y, z)
print(f"  Verification: Q(1/sqrt(m1), 1/sqrt(m2), 1/sqrt(m4)) = {Q_check:.15f}")
print(f"  Deviation from 2/3: {Q_check - 2/3:.2e}")
print()

# ─── Check: is eq. 4 automatically satisfied? ───────────────────────────
print("--- Self-consistency check ---")
print("If eq. 4 were automatic (redundant), Q would equal 2/3 for ANY m1.")
print("Testing several m1 values with m0 fixed:")
print()
m0_fixed = 2.16
for m1_trial in [2.0, 4.0, 6.0, 8.0, 10.0, 20.0]:
    S = m0_fixed + m1_trial
    m2_t = R * S
    m4_t = (3*np.sqrt(R) + R)**2 * S
    x = 1.0/np.sqrt(m1_trial)
    y = 1.0/np.sqrt(m2_t)
    z = 1.0/np.sqrt(m4_t)
    Qval = Q(x, y, z)
    print(f"  m1 = {m1_trial:6.2f}  => Q = {Qval:.10f}  (deviation from 2/3: {Qval - 2/3:+.6e})")

print()
print("CONCLUSION: Q depends on m1, so eq. 4 is NOT automatic.")
print("Eq. 4 provides a genuine constraint, reducing the system to 1 free parameter (m0).")
print()

# =====================================================================
# PART (b): Find optimal m0
# =====================================================================
print("=" * 70)
print("PART (b): Optimal m0 to match experimental targets")
print("=" * 70)
print()

# Chi-squared fit
def chi2(m0):
    try:
        m1, m2, m3, m4 = solve_system(m0)
    except:
        return 1e10
    chi2_val = 0.0
    chi2_val += ((m0 - targets['m0'][0]) / targets['m0'][1])**2
    chi2_val += ((m1 - targets['m1'][0]) / targets['m1'][1])**2
    chi2_val += ((m2 - targets['m2'][0]) / targets['m2'][1])**2
    chi2_val += ((m3 - targets['m3'][0]) / targets['m3'][1])**2
    chi2_val += ((m4 - targets['m4'][0]) / targets['m4'][1])**2
    return chi2_val

# Scan m0
print("Scanning m0 from 0.5 to 5.0:")
m0_scan = np.linspace(0.5, 5.0, 1000)
chi2_scan = []
for m0_val in m0_scan:
    try:
        chi2_scan.append(chi2(m0_val))
    except:
        chi2_scan.append(1e10)
chi2_scan = np.array(chi2_scan)
best_idx = np.argmin(chi2_scan)
m0_best_approx = m0_scan[best_idx]
print(f"  Approximate minimum at m0 = {m0_best_approx:.4f}, chi2 = {chi2_scan[best_idx]:.4f}")
print()

# Refine with Brent
from scipy.optimize import minimize_scalar
result = minimize_scalar(chi2, bounds=(0.5, 5.0), method='bounded')
m0_best = result.x
print(f"  Refined minimum at m0 = {m0_best:.8f}, chi2 = {result.fun:.8f}")
print()

m1_best, m2_best, m3_best, m4_best = solve_system(m0_best)
print(f"  Optimal masses:")
print(f"    m0 = {m0_best:.6f}   (target: {targets['m0'][0]} +/- {targets['m0'][1]}, pull: {(m0_best - targets['m0'][0])/targets['m0'][1]:+.2f} sigma)")
print(f"    m1 = {m1_best:.6f}   (target: {targets['m1'][0]} +/- {targets['m1'][1]}, pull: {(m1_best - targets['m1'][0])/targets['m1'][1]:+.2f} sigma)")
print(f"    m2 = {m2_best:.4f}   (target: {targets['m2'][0]} +/- {targets['m2'][1]}, pull: {(m2_best - targets['m2'][0])/targets['m2'][1]:+.2f} sigma)")
print(f"    m3 = {m3_best:.2f}   (target: {targets['m3'][0]} +/- {targets['m3'][1]}, pull: {(m3_best - targets['m3'][0])/targets['m3'][1]:+.2f} sigma)")
print(f"    m4 = {m4_best:.2f}   (target: {targets['m4'][0]} +/- {targets['m4'][1]}, pull: {(m4_best - targets['m4'][0])/targets['m4'][1]:+.2f} sigma)")
print(f"    chi2/dof = {result.fun:.4f} / 4 = {result.fun/4:.4f}")
print()

# Also show individual chi2 contributions
print("  Individual chi2 contributions:")
pulls = [
    ('m0', m0_best, targets['m0']),
    ('m1', m1_best, targets['m1']),
    ('m2', m2_best, targets['m2']),
    ('m3', m3_best, targets['m3']),
    ('m4', m4_best, targets['m4']),
]
for name, val, (tgt, err) in pulls:
    p = (val - tgt) / err
    print(f"    {name}: pull = {p:+.4f} sigma, chi2_contrib = {p**2:.4f}")
print()

# ─── Also try: what m0 makes each mass hit its central value? ─────────
print("--- What m0 is required by each mass individually? ---")
for name, m_target, m_err in [
    ('m2', targets['m2'][0], targets['m2'][1]),
    ('m3', targets['m3'][0], targets['m3'][1]),
    ('m4', targets['m4'][0], targets['m4'][1]),
]:
    # m2 = R*(m0+m1), m3 = R^2*(m0+m1), m4 = coeff*(m0+m1)
    # Each determines S = m0+m1, then eq. 4 determines m1(m0), so m0 is fixed
    if name == 'm2':
        S_needed = m_target / R
    elif name == 'm3':
        S_needed = m_target / R**2
    elif name == 'm4':
        S_needed = m_target / (3*np.sqrt(R) + R)**2
    print(f"  {name} = {m_target} requires S = m0+m1 = {S_needed:.4f}")

print()

# =====================================================================
# PART (c): Sensitivity dm_i/dm_0
# =====================================================================
print("=" * 70)
print("PART (c): Sensitivity analysis: dm_i / dm_0")
print("=" * 70)
print()

# Numerical differentiation
dm0 = 1e-6
m1_p, m2_p, m3_p, m4_p = solve_system(m0_best + dm0)
m1_m, m2_m, m3_m, m4_m = solve_system(m0_best - dm0)

dm1_dm0 = (m1_p - m1_m) / (2 * dm0)
dm2_dm0 = (m2_p - m2_m) / (2 * dm0)
dm3_dm0 = (m3_p - m3_m) / (2 * dm0)
dm4_dm0 = (m4_p - m4_m) / (2 * dm0)

print(f"At m0 = {m0_best:.6f}:")
print(f"  dm1/dm0 = {dm1_dm0:.6f}")
print(f"  dm2/dm0 = {dm2_dm0:.4f}")
print(f"  dm3/dm0 = {dm3_dm0:.2f}")
print(f"  dm4/dm0 = {dm4_dm0:.2f}")
print()
print(f"  Relative sensitivities (elasticity d ln m_i / d ln m_0):")
print(f"  d ln m1 / d ln m0 = {dm1_dm0 * m0_best / m1_best:.6f}")
print(f"  d ln m2 / d ln m0 = {dm2_dm0 * m0_best / m2_best:.6f}")
print(f"  d ln m3 / d ln m0 = {dm3_dm0 * m0_best / m3_best:.6f}")
print(f"  d ln m4 / d ln m0 = {dm4_dm0 * m0_best / m4_best:.6f}")
print()

# Propagated uncertainty from m0 alone
sigma_m0 = targets['m0'][1]
print(f"  Propagated uncertainties from sigma(m0) = {sigma_m0}:")
print(f"  sigma(m1) = {abs(dm1_dm0) * sigma_m0:.4f}  (expt: {targets['m1'][1]})")
print(f"  sigma(m2) = {abs(dm2_dm0) * sigma_m0:.2f}  (expt: {targets['m2'][1]})")
print(f"  sigma(m3) = {abs(dm3_dm0) * sigma_m0:.1f}  (expt: {targets['m3'][1]})")
print(f"  sigma(m4) = {abs(dm4_dm0) * sigma_m0:.1f}  (expt: {targets['m4'][1]})")
print()

# =====================================================================
# PART (d): Energy balance on direct masses -> m5
# =====================================================================
print("=" * 70)
print("PART (d): Determine m5 from Q(sqrt(m2), sqrt(m3), sqrt(m5)) = 2/3")
print("=" * 70)
print()

# Use the best-fit m2, m3
m2_use = m2_best
m3_use = m3_best

print(f"Using m2 = {m2_use:.4f}, m3 = {m3_use:.4f}")
print()

# Q(sqrt(m2), sqrt(m3), sqrt(m5)) = 2/3
# Let a = sqrt(m2), b = sqrt(m3), c = sqrt(m5)
# (a^2 + b^2 + c^2) / (a + b + c)^2 = 2/3
# 3(a^2 + b^2 + c^2) = 2(a + b + c)^2
# 3a^2 + 3b^2 + 3c^2 = 2a^2 + 2b^2 + 2c^2 + 4ab + 4ac + 4bc
# a^2 + b^2 + c^2 - 4ab - 4ac - 4bc = 0
# c^2 - 4(a+b)c + (a^2 + b^2 - 4ab) = 0

a = np.sqrt(m2_use)
b = np.sqrt(m3_use)

print(f"a = sqrt(m2) = {a:.6f}")
print(f"b = sqrt(m3) = {b:.6f}")
print()

# Quadratic in c = sqrt(m5)
A_coeff = 1.0
B_coeff = -4.0 * (a + b)
C_coeff = a**2 + b**2 - 4*a*b

discriminant = B_coeff**2 - 4*A_coeff*C_coeff
print(f"Quadratic: c^2 - 4(a+b)c + (a^2 + b^2 - 4ab) = 0")
print(f"  Coefficients: A={A_coeff}, B={B_coeff:.6f}, C={C_coeff:.6f}")
print(f"  Discriminant = {discriminant:.6f}")
print()

c_plus = (-B_coeff + np.sqrt(discriminant)) / (2*A_coeff)
c_minus = (-B_coeff - np.sqrt(discriminant)) / (2*A_coeff)

m5_plus = c_plus**2
m5_minus = c_minus**2

print(f"Solutions:")
print(f"  c+ = {c_plus:.6f}  => m5 = {m5_plus:.2f}")
print(f"  c- = {c_minus:.6f}  => m5 = {m5_minus:.4f}")
print()

# The physical solution (larger one, matching target ~162500)
m5_phys = m5_plus
print(f"Physical solution: m5 = {m5_phys:.2f}")
print(f"  Target: {targets['m5'][0]}")
print(f"  Ratio m5/target = {m5_phys/targets['m5'][0]:.6f}")
print(f"  Percentage deviation = {100*(m5_phys - targets['m5'][0])/targets['m5'][0]:.2f}%")
print()

# Verify
Q_check_d = Q(np.sqrt(m2_use), np.sqrt(m3_use), np.sqrt(m5_phys))
print(f"  Verification: Q(sqrt(m2), sqrt(m3), sqrt(m5)) = {Q_check_d:.15f}")
print(f"  Deviation from 2/3: {Q_check_d - 2/3:.2e}")
print()

# ─── Also compute m5 using exact target masses ──────────────────────────
print("--- m5 from exact target m2, m3 ---")
a_t = np.sqrt(targets['m2'][0])
b_t = np.sqrt(targets['m3'][0])
B_t = -4.0 * (a_t + b_t)
C_t = a_t**2 + b_t**2 - 4*a_t*b_t
disc_t = B_t**2 - 4*C_t
c_t_plus = (-B_t + np.sqrt(disc_t)) / 2
m5_from_targets = c_t_plus**2
print(f"  Using m2={targets['m2'][0]}, m3={targets['m3'][0]}:")
print(f"  m5 = {m5_from_targets:.2f}")
print(f"  Target: {targets['m5'][0]}")
print(f"  Deviation = {100*(m5_from_targets - targets['m5'][0])/targets['m5'][0]:.2f}%")
print()

# =====================================================================
# PART (d) continued: the other solution
# =====================================================================
print("--- The smaller solution ---")
print(f"  m5_small = {m5_minus:.6f}")
print(f"  This corresponds to sqrt(m5) = {c_minus:.6f}")
print()

# =====================================================================
# SUMMARY TABLE
# =====================================================================
print("=" * 70)
print("SUMMARY")
print("=" * 70)
print()
print(f"Free parameter: m0 = {m0_best:.6f}")
print(f"R = (2+sqrt(3))^2 = {R:.10f}")
print(f"S = m0 + m1 = {m0_best + m1_best:.6f}")
print()
print(f"{'Mass':>6s}  {'Predicted':>12s}  {'Target':>10s}  {'Unc':>6s}  {'Pull (sigma)':>12s}")
print("-" * 55)
for name, val, (tgt, err) in [
    ('m0', m0_best, targets['m0']),
    ('m1', m1_best, targets['m1']),
    ('m2', m2_best, targets['m2']),
    ('m3', m3_best, targets['m3']),
    ('m4', m4_best, targets['m4']),
]:
    pull = (val - tgt) / err
    print(f"{name:>6s}  {val:12.4f}  {tgt:10.2f}  {err:6.2f}  {pull:+12.4f}")
print()
print(f"m5 = {m5_phys:.2f}  (target: {targets['m5'][0]:.0f})")
print()

# ─── Analytic expressions ────────────────────────────────────────────────
print("=" * 70)
print("ANALYTIC STRUCTURE")
print("=" * 70)
print()
print("Let S = m0 + m1. Then:")
print(f"  m2 = R * S = {R:.6f} * S")
print(f"  m3 = R^2 * S = {R**2:.4f} * S")
print(f"  m4 = (3*sqrt(R) + R)^2 * S = {(3*np.sqrt(R) + R)**2:.4f} * S")
print()
print("Eq. 4 determines m1 = m1(m0). Defining t = m0/m1:")
print()

t_best = m0_best / m1_best
print(f"  t = m0/m1 = {t_best:.8f}")
print(f"  m1 = S/(1+t), m0 = tS/(1+t)")
print()

# Write out the Q = 2/3 equation more explicitly
print("The Q = 2/3 constraint becomes (after substituting):")
print("  Let u = 1/sqrt(m1), v = 1/sqrt(R*S), w = 1/sqrt(C*S)")
print(f"  where C = (3*sqrt(R)+R)^2 = {(3*np.sqrt(R)+R)**2:.6f}")
print()
print("  Q = (u^2 + v^2 + w^2) / (u + v + w)^2 = 2/3")
print()
print("  Since v = 1/sqrt(R*S) and w = 1/sqrt(C*S),")
print("  and u = 1/sqrt(S-m0) = 1/sqrt(m1),")
print("  this mixes m0 and S nontrivially through m1 = S - m0.")
print()

# ─── Ratio check ─────────────────────────────────────────────────────────
print("=" * 70)
print("RATIO CHECKS")
print("=" * 70)
print()
print(f"  m2/m1 = {m2_best/m1_best:.4f}  (R = {R:.4f})")
print(f"  m3/m2 = {m3_best/m2_best:.4f}  (R = {R:.4f})")
print(f"  m3/m1 = {m3_best/m1_best:.4f}  (R^2 = {R**2:.4f})")
print(f"  m4/m2 = {m4_best/m2_best:.4f}  ((3+sqrt(R))^2... no, see formula)")
print(f"  sqrt(m4)/sqrt(m2) = {np.sqrt(m4_best)/np.sqrt(m2_best):.6f}")
print(f"  3 + sqrt(m3/m2) = {3 + np.sqrt(m3_best/m2_best):.6f}")
print(f"  (should be equal by eq. 3: sqrt(m4) = 3*sqrt(m2) + sqrt(m3))")
print(f"  => sqrt(m4)/sqrt(m2) = 3 + sqrt(R) = {3 + np.sqrt(R):.6f}")
print()

# ─── Symbolic exact solution attempt ─────────────────────────────────────
print("=" * 70)
print("EXACT SYMBOLIC ANALYSIS")
print("=" * 70)
print()
print("Can we solve Q = 2/3 analytically for the ratio t = m0/m1?")
print()
print("Let r = m0/S = t/(1+t). Then m1 = S(1-r), m2 = RS, m4 = CS.")
print("  u^2 = 1/(S(1-r))")
print("  v^2 = 1/(RS)")
print("  w^2 = 1/(CS)")
print()
print("  Q = [1/((1-r)S) + 1/(RS) + 1/(CS)] / [1/sqrt((1-r)S) + 1/sqrt(RS) + 1/sqrt(CS)]^2")
print()
print("Factor out 1/S from numerator and 1/sqrt(S) from denominator:")
print("  Q = (1/S)[1/(1-r) + 1/R + 1/C] / (1/S)[1/sqrt(1-r) + 1/sqrt(R) + 1/sqrt(C)]^2")
print("  Q = [1/(1-r) + 1/R + 1/C] / [1/sqrt(1-r) + 1/sqrt(R) + 1/sqrt(C)]^2")
print()
print("S cancels completely! Q depends only on r = m0/(m0+m1), not on the overall scale.")
print("This means: given the STRUCTURE (eqs 1-3), eq. 4 fixes the RATIO m0/m1,")
print("and the overall scale (e.g., m0 itself) remains free.")
print()

# Solve for r analytically
def Q_of_r(r):
    """Q as function of r = m0/S, with R and C fixed."""
    C = (3*np.sqrt(R) + R)**2
    if r <= 0 or r >= 1:
        return 1.0  # invalid
    a = 1.0/(1-r)
    b = 1.0/R
    c = 1.0/C
    num = a + b + c
    sa = 1.0/np.sqrt(1-r)
    sb = 1.0/np.sqrt(R)
    sc = 1.0/np.sqrt(C)
    den = (sa + sb + sc)**2
    return num / den

# Find r such that Q(r) = 2/3
r_sol = brentq(lambda r: Q_of_r(r) - 2/3, 0.001, 0.999)
t_sol = r_sol / (1 - r_sol)
print(f"EXACT ratio (numerical):")
print(f"  r = m0/(m0+m1) = {r_sol:.12f}")
print(f"  t = m0/m1 = {t_sol:.12f}")
print(f"  m1/m0 = {1/t_sol:.12f}")
print()
print(f"  Check: Q(r_sol) = {Q_of_r(r_sol):.15f}")
print()

# So the system is: fix ratio m0/m1, then m0 is free parameter
print(f"All masses as functions of m0:")
print(f"  m1 = m0 / t = m0 / {t_sol:.10f} = {1/t_sol:.10f} * m0")
print(f"  S = m0 + m1 = m0 * (1 + 1/t) = {1 + 1/t_sol:.10f} * m0")
S_over_m0 = 1 + 1/t_sol
print(f"  m2 = R * S = {R * S_over_m0:.8f} * m0")
print(f"  m3 = R^2 * S = {R**2 * S_over_m0:.6f} * m0")
C = (3*np.sqrt(R) + R)**2
print(f"  m4 = C * S = {C * S_over_m0:.4f} * m0")
print()

# ─── What m0 gives the best fit? ────────────────────────────────────────
# With fixed ratio, each mass is proportional to m0.
# m_i = k_i * m0 for known k_i
k1 = 1.0 / t_sol
k2 = R * S_over_m0
k3 = R**2 * S_over_m0
k4 = C * S_over_m0

print("Proportionality constants k_i (where m_i = k_i * m_0):")
print(f"  k1 = {k1:.10f}")
print(f"  k2 = {k2:.8f}")
print(f"  k3 = {k3:.6f}")
print(f"  k4 = {k4:.4f}")
print()

# Weighted least squares for m0
# Minimize sum_i ((k_i * m0 - target_i) / sigma_i)^2
# Include m0 itself: target_0 = 2.16, k_0 = 1
all_k = [1.0, k1, k2, k3, k4]
all_tgt = [targets['m0'][0], targets['m1'][0], targets['m2'][0], targets['m3'][0], targets['m4'][0]]
all_sig = [targets['m0'][1], targets['m1'][1], targets['m2'][1], targets['m3'][1], targets['m4'][1]]

# Analytic solution: m0 = sum(k_i * target_i / sigma_i^2) / sum(k_i^2 / sigma_i^2)
num_wls = sum(k * t / s**2 for k, t, s in zip(all_k, all_tgt, all_sig))
den_wls = sum(k**2 / s**2 for k, t, s in zip(all_k, all_tgt, all_sig))
m0_wls = num_wls / den_wls

print(f"Weighted least squares optimal m0 = {m0_wls:.8f}")
print()
print(f"Final predicted masses:")
for name, k, (tgt, sig) in zip(['m0','m1','m2','m3','m4'], all_k,
                                 [targets['m0'], targets['m1'], targets['m2'], targets['m3'], targets['m4']]):
    val = k * m0_wls
    pull = (val - tgt) / sig
    print(f"  {name} = {val:12.4f}  (target: {tgt:10.2f} +/- {sig}, pull: {pull:+.4f} sigma)")

chi2_final = sum(((k*m0_wls - tgt)/sig)**2 for k, tgt, sig in zip(all_k, all_tgt, all_sig))
ndof = len(all_k) - 1  # 5 measurements, 1 parameter
print(f"\n  chi2 = {chi2_final:.4f}, ndof = {ndof}, chi2/ndof = {chi2_final/ndof:.4f}")
print()

# ─── m5 at optimal ───────────────────────────────────────────────────────
m2_final = k2 * m0_wls
m3_final = k3 * m0_wls
a_f = np.sqrt(m2_final)
b_f = np.sqrt(m3_final)
B_f = -4*(a_f + b_f)
C_f = a_f**2 + b_f**2 - 4*a_f*b_f
disc_f = B_f**2 - 4*C_f
c_f = (-B_f + np.sqrt(disc_f)) / 2
m5_final = c_f**2

print(f"m5 from Q(sqrt(m2), sqrt(m3), sqrt(m5)) = 2/3:")
print(f"  m5 = {m5_final:.2f}  (target: {targets['m5'][0]:.0f})")
print(f"  Deviation: {100*(m5_final - targets['m5'][0])/targets['m5'][0]:.2f}%")
print()

# ─── Can we express m5/m0 analytically? ──────────────────────────────────
print("=" * 70)
print("m5 AS A FUNCTION OF m0")
print("=" * 70)
print()
print(f"m5 = k5 * m0 where k5 is determined by the quadratic.")
print(f"Since m2 = k2*m0, m3 = k3*m0:")
print(f"  a = sqrt(k2*m0), b = sqrt(k3*m0)")
print(f"  c = 2(a+b) + sqrt(3(a+b)^2 + 4ab)  ... from quadratic")
print(f"  c = sqrt(m0) * [2(sqrt(k2)+sqrt(k3)) + sqrt(3(sqrt(k2)+sqrt(k3))^2 + 4*sqrt(k2*k3))]")
print(f"  m5 = c^2 = m0 * [...]^2")
print()

sk2 = np.sqrt(k2)
sk3 = np.sqrt(k3)
c_coeff = 2*(sk2 + sk3) + np.sqrt(3*(sk2+sk3)**2 + 4*sk2*sk3)
# Actually let me redo this properly from the quadratic
# c^2 - 4(a+b)c + (a^2+b^2-4ab) = 0
# a = sqrt(k2*m0) = sqrt(k2)*sqrt(m0)
# Let a' = sqrt(k2), b' = sqrt(k3)
# c = sqrt(m0) * c', where c'^2 - 4(a'+b')c' + (a'^2 + b'^2 - 4a'b') = 0
aprime = np.sqrt(k2)
bprime = np.sqrt(k3)
Bp = -4*(aprime + bprime)
Cp = aprime**2 + bprime**2 - 4*aprime*bprime
discp = Bp**2 - 4*Cp
cprime_plus = (-Bp + np.sqrt(discp)) / 2
k5 = cprime_plus**2

print(f"k5 = m5/m0 = {k5:.4f}")
print(f"m5/m0 = {m5_final/m0_wls:.4f}")
print(f"m5 = {k5:.4f} * m0 = {k5 * m0_wls:.2f}")
print()

# =====================================================================
# PART (d) VARIANT: Q on (m3, m4, m5) instead of (m2, m3, m5)
# =====================================================================
print("=" * 70)
print("PART (d) VARIANT: Q(sqrt(m3), sqrt(m4), sqrt(m5)) = 2/3")
print("  (using the next triple in the chain)")
print("=" * 70)
print()

# Using best-fit m3, m4
m3_v = m3_best
m4_v = m4_best
av = np.sqrt(m3_v)
bv = np.sqrt(m4_v)
Bv = -4*(av + bv)
Cv = m3_v + m4_v - 4*np.sqrt(m3_v*m4_v)
discv = Bv**2 - 4*Cv
cv_plus = (-Bv + np.sqrt(discv))/2
cv_minus = (-Bv - np.sqrt(discv))/2
m5_var_plus = cv_plus**2
m5_var_minus = cv_minus**2

print(f"Using m3 = {m3_v:.4f}, m4 = {m4_v:.4f}")
print(f"Quadratic: c^2 - 4(sqrt(m3)+sqrt(m4))c + (m3+m4-4*sqrt(m3*m4)) = 0")
print(f"  B = {Bv:.6f}")
print(f"  C = {Cv:.6f}")
print(f"  discriminant = {discv:.6f}")
print()
print(f"Solutions:")
print(f"  c+ = {cv_plus:.4f}  => m5 = {m5_var_plus:.2f}")
print(f"  c- = {cv_minus:.4f}  => m5 = {m5_var_minus:.4f}")
print(f"Target: {targets['m5'][0]:.0f}")
print(f"Deviation of m5+: {100*(m5_var_plus - targets['m5'][0])/targets['m5'][0]:.2f}%")
print()

# Verify
Q_var = Q(np.sqrt(m3_v), np.sqrt(m4_v), np.sqrt(m5_var_plus))
print(f"Verification: Q(sqrt(m3), sqrt(m4), sqrt(m5)) = {Q_var:.15f}")
print()

# Also with PDG central values
av2 = np.sqrt(1270.0)
bv2 = np.sqrt(4183.0)
Bv2 = -4*(av2+bv2)
Cv2 = 1270.0 + 4183.0 - 4*np.sqrt(1270.0*4183.0)
discv2 = Bv2**2 - 4*Cv2
cv2_plus = (-Bv2 + np.sqrt(discv2))/2
m5_pdg = cv2_plus**2
print(f"With PDG targets m3=1270, m4=4183:")
print(f"  m5 = {m5_pdg:.2f}")
print(f"  sqrt(m5) = {cv2_plus:.2f}")
print(f"  Deviation from 162500: {100*(m5_pdg - 162500)/162500:.2f}%")
print()

# =====================================================================
# FINAL COMPREHENSIVE SUMMARY
# =====================================================================
print("=" * 70)
print("FINAL COMPREHENSIVE SUMMARY")
print("=" * 70)
print()
print("The system of 4 equations in 6 unknowns has the following structure:")
print()
print("1. Equations 1-3 express m2, m3, m4 in terms of S = m0 + m1.")
print("2. Equation 4 (energy balance) fixes the RATIO m0/m1 = t.")
print("   Crucially, S cancels from Q, so eq. 4 does NOT fix the scale.")
print("3. Therefore the system has exactly ONE free parameter: m0 (or equivalently S).")
print()
print(f"The universal ratio is: t = m0/m1 = {t_sol:.10f}")
print(f"                        m0/m1 = {t_sol:.10f}")
print()
print(f"All masses scale linearly with m0:")
print(f"  m0 = m0")
print(f"  m1 = {k1:.8f} * m0")
print(f"  m2 = {k2:.6f} * m0")
print(f"  m3 = {k3:.4f} * m0")
print(f"  m4 = {k4:.2f} * m0")
print(f"  m5 = {k5:.2f} * m0")
print()
print(f"Best-fit m0 = {m0_wls:.6f} (weighted least squares)")
print(f"chi2 = {chi2_final:.4f} for {ndof} dof (p-value: {1 - __import__('scipy').stats.chi2.cdf(chi2_final, ndof):.4f})")
print()
print(f"Sensitivities (all masses proportional to m0, so dm_i/dm_0 = k_i):")
print(f"  dm1/dm0 = {k1:.6f}")
print(f"  dm2/dm0 = {k2:.4f}")
print(f"  dm3/dm0 = {k3:.2f}")
print(f"  dm4/dm0 = {k4:.1f}")
print(f"  dm5/dm0 = {k5:.1f}")
