"""
Bion-induced effective potential for a pseudo-modulus in SU(3) on R^3 x S^1.

We compute the bion potential arising from correlated monopole-instanton /
anti-monopole-instanton pairs in SU(3) gauge theory with N_f massive
fundamental flavors on R^3 x S^1.  The calculation uses the extended
Dynkin diagram for SU(3), enumerates all magnetic bion types, determines
their mass dependence through fermionic zero modes, expands the bion
potential to quartic order in a pseudo-modulus X, extracts the quartic
Kahler correction c_bion, and checks the condition c_bion = -1/12.

Numerical inputs:
    Lambda_QCD = 300 MeV
    alpha_s(Lambda_QCD) = 0.3
    N_c = 3
    L = 1 / Lambda_QCD
"""

import numpy as np
from fractions import Fraction
from scipy.optimize import brentq

# ---------------------------------------------------------------------------
# Physical parameters
# ---------------------------------------------------------------------------
Lambda_QCD = 300.0       # MeV
alpha_s = 0.3            # at confinement scale (rough)
Nc = 3
g2 = 4 * np.pi * alpha_s  # g^2 = 4 pi alpha_s
L = 1.0 / Lambda_QCD     # natural compactification, MeV^{-1}

# Monopole-instanton action
S0 = 8 * np.pi**2 / (g2 * Nc)           # = 2 pi / (alpha_s N_c)
S0_alt = 2 * np.pi / (alpha_s * Nc)     # same formula, check
S_bion = 2 * S0 / Nc                    # bion action

print("=" * 72)
print("BION-INDUCED KAHLER POTENTIAL FOR SU(3) ON R^3 x S^1")
print("=" * 72)
print()
print(f"Parameters:")
print(f"  Lambda_QCD = {Lambda_QCD} MeV")
print(f"  alpha_s    = {alpha_s}")
print(f"  g^2        = 4 pi alpha_s = {g2:.6f}")
print(f"  N_c        = {Nc}")
print(f"  L          = 1/Lambda   = {L:.6e} MeV^(-1)")
print(f"  S_0        = 8 pi^2/(g^2 N_c) = {S0:.6f}")
print(f"  S_0 (alt)  = 2 pi/(alpha_s N_c) = {S0_alt:.6f}")
print(f"  S_bion     = 2 S_0/N_c   = {S_bion:.6f}")
print(f"  exp(-S_bion) = {np.exp(-S_bion):.6e}")
print()

# ===================================================================
# TASK 1: Enumerate magnetic bion types for SU(3)
# ===================================================================

print("=" * 72)
print("TASK 1: Magnetic Bion Types for SU(3)")
print("=" * 72)
print()

# SU(3) simple roots in the weight basis
# alpha_1 = (1, -1, 0), alpha_2 = (0, 1, -1)
# Affine root: alpha_0 = -(alpha_1 + alpha_2) = (-1, 0, 1)
#
# Extended Dynkin diagram for SU(3):
#
#       alpha_0
#      /       \
#  alpha_1 --- alpha_2
#
# It is a triangle (Z_3 symmetric).  Every pair of roots is adjacent.

roots = {
    'alpha_1': np.array([1, -1, 0]),
    'alpha_2': np.array([0, 1, -1]),
    'alpha_0': np.array([-1, 0, 1]),
}

adjacent_pairs = [
    ('alpha_1', 'alpha_2'),
    ('alpha_2', 'alpha_0'),
    ('alpha_0', 'alpha_1'),
]

print("SU(3) extended Dynkin diagram: triangle (Z_3 symmetric)")
print()
print("Simple roots (in weight space R^3):")
for name, vec in roots.items():
    print(f"  {name} = {vec}")
print()
print(f"Affine root: alpha_0 = -(alpha_1 + alpha_2) = {roots['alpha_0']}")
print()

print("Adjacent pairs in the extended Dynkin diagram:")
for a, b in adjacent_pairs:
    print(f"  ({a}, {b})")
print()

# Magnetic bions: for each adjacent pair (a, b), the bion is M_a M_b-bar.
# The conjugate M_b M_a-bar is also a distinct object.
# Total including orientation: 6. Distinct up to conjugation: 3.

bion_types = []
for a, b in adjacent_pairs:
    bion_types.append((a, b))
    bion_types.append((b, a))

print("Magnetic bion types (monopole_a - anti-monopole_b):")
for i, (a, b) in enumerate(bion_types):
    diff = roots[a] - roots[b]
    print(f"  B_{i+1}: ({a}, {b}-bar)  charge = {diff}")
print()
print(f"Total distinct bion types (including orientation): {len(bion_types)}")
print(f"Distinct up to conjugation: {len(adjacent_pairs)}")
print()
print("Each adjacent pair contributes a bion and its conjugate (anti-bion).")
print("Together they produce a cosine potential in the dual photon.")
print()

# ===================================================================
# TASK 2: Mass dependence from fermionic zero modes
# ===================================================================

print("=" * 72)
print("TASK 2: Quark Mass Dependence of Bion Amplitudes")
print("=" * 72)
print()

# In the center-symmetric vacuum, each monopole-instanton on any root
# (alpha_0, alpha_1, alpha_2) has exactly one fermionic zero mode per
# fundamental flavor.
#
# For a bion on (alpha_a, alpha_b):
#   A_ab propto Product_f (m_f L)^{n_{a,f} + n_{b,f}}
#              = Product_f (m_f L)^{1 + 1}
#              = Product_f (m_f L)^2
#              = [Product_f m_f]^2 L^{2 N_f}

print("Zero-mode counting:")
print("  Each root (alpha_0, alpha_1, alpha_2) has n_{a,f} = 1")
print("  for every fundamental flavor f = 1, ..., N_f.")
print()
print("Mass dependence for bion on (alpha_a, alpha_b):")
print("  A_ab propto Product_f (m_f L)^{n_{a,f} + n_{b,f}}")
print("             = Product_f (m_f L)^2")
print("             = [Product_f (m_f L)]^2")
print()
print("All three bion types have IDENTICAL mass dependence.")
print("(Z_3 symmetry of extended Dynkin diagram, democratic zero modes.)")
print()

m_quarks = {'u': 2.16, 'd': 4.67, 's': 93.4, 'c': 1270.0, 'b': 4180.0}

print("Numerical example with N_f = 5 light quarks (u, d, s, c, b):")
print()
for name, m in m_quarks.items():
    print(f"  m_{name} L = m_{name}/Lambda = {m/Lambda_QCD:.6f}")

mass_product = np.prod(list(m_quarks.values()))
mass_factor = (mass_product * L**5)**2
print()
print(f"  Product(m_f) = {mass_product:.4e} MeV^5")
print(f"  [Product(m_f L)]^2 = {mass_factor:.4e} (dimensionless)")
print()

# ===================================================================
# TASK 3: Bion potential as function of pseudo-modulus X
# ===================================================================

print("=" * 72)
print("TASK 3: Bion Potential V_bion(X) and Quartic Expansion")
print("=" * 72)
print()

# Holonomy charges for the three roots in center-symmetric SU(3) vacuum:
#   q_0 = 0,  q_1 = 1,  q_2 = -1   (mod 3)

q = {'alpha_0': 0, 'alpha_1': 1, 'alpha_2': -1}

print("Holonomy charges of roots (center-symmetric vacuum):")
for name, charge in q.items():
    print(f"  q_{name[-1]} = {charge}")
print()

# For each adjacent pair, the cosine argument is 2 pi (q_a - q_b) X / Lambda
print("Charge differences for adjacent pairs:")

charge_diffs = []
for a, b in adjacent_pairs:
    dq = q[a] - q[b]
    charge_diffs.append(dq)
    print(f"  ({a}, {b}): q_{a[-1]} - q_{b[-1]} = {q[a]} - ({q[b]}) = {dq}")

print()

# The three terms:
# (alpha_1, alpha_2): dq = 2   -> cos(4 pi X/Lambda)
# (alpha_2, alpha_0): dq = -1  -> cos(-2 pi X/Lambda) = cos(2 pi X/Lambda)
# (alpha_0, alpha_1): dq = -1  -> cos(-2 pi X/Lambda) = cos(2 pi X/Lambda)
#
# Sum: cos(4 pi X/Lambda) + 2 cos(2 pi X/Lambda)

print("V_bion(X) = A exp(-S_bion) [cos(4 pi X/Lambda) + 2 cos(2 pi X/Lambda)]")
print()

# Taylor expansion.  Let phi = 2 pi X / Lambda.
# cos(2 phi) = 1 - 2 phi^2 + (2/3) phi^4 - ...
# 2 cos(phi) = 2 - phi^2 + (1/12) phi^4 - ...
# Sum = 3 - 3 phi^2 + (2/3 + 1/12) phi^4 - ...
#     = 3 - 3 phi^2 + (3/4) phi^4 - ...

print("Taylor expansion (phi = 2 pi X / Lambda):")
print()

c2_cos2phi = Fraction(-2, 1)
c4_cos2phi = Fraction(2, 3)
c2_2cosphi = Fraction(-1, 1)
c4_2cosphi = Fraction(1, 12)

c0_total = Fraction(3, 1)
c2_total = c2_cos2phi + c2_2cosphi
c4_total = c4_cos2phi + c4_2cosphi

print(f"  cos(2 phi): phi^2 coeff = {c2_cos2phi}, phi^4 coeff = {c4_cos2phi}")
print(f"  2 cos(phi): phi^2 coeff = {c2_2cosphi}, phi^4 coeff = {c4_2cosphi}")
print()
print(f"  Sum:  constant = {c0_total}")
print(f"        phi^2 coeff = {c2_total}")
print(f"        phi^4 coeff = {c4_total} = {float(c4_total):.6f}")
print()

# In terms of X:
# phi^2 = (2 pi)^2 (X/Lambda)^2
# phi^4 = (2 pi)^4 (X/Lambda)^4
c2_num = float(c2_total) * (2*np.pi)**2
c4_num = float(c4_total) * (2*np.pi)**4

print("In terms of X/Lambda:")
print(f"  V_bion = A exp(-S_bion) [3  {c2_total}*(2pi)^2 (X/Lambda)^2")
print(f"                              + {c4_total}*(2pi)^4 (X/Lambda)^4 + ...]")
print()
print(f"        = A exp(-S_bion) [3  - 12 pi^2 (X/Lambda)^2")
print(f"                              + 12 pi^4 (X/Lambda)^4 + ...]")
print()
print(f"Numerical coefficients:")
print(f"  (X/Lambda)^0: 3")
print(f"  (X/Lambda)^2: {c2_num:.6f}")
print(f"  (X/Lambda)^4: {c4_num:.6f}")
print()

# Quartic coefficient explicitly
print(f"Quartic coefficient of |X|^4/Lambda^4:")
print(f"  = (3/4)(2 pi)^4 = 12 pi^4 = {(3.0/4)*(2*np.pi)**4:.4f}")
print()

# Numerical verification of the expansion
print("Numerical verification (comparing exact vs Taylor at small X/Lambda):")
phi_test = np.array([0.01, 0.05, 0.1, 0.2, 0.3])
for phi in phi_test:
    exact = np.cos(2*phi) + 2*np.cos(phi)
    taylor = 3 - 3*phi**2 + 0.75*phi**4
    print(f"  phi = {phi:.2f}: exact = {exact:.10f}, Taylor = {taylor:.10f}, "
          f"error = {abs(exact-taylor):.2e}")
print()

# ===================================================================
# TASK 4: When does the quartic coefficient equal -1/12?
# ===================================================================

print("=" * 72)
print("TASK 4: Condition for Quartic Coefficient = -1/12")
print("=" * 72)
print()

# The bion amplitude (overall normalization).
# Semiclassical computation gives the bion potential as:
#   V_bion = - (Lambda^4 / g^4) exp(-S_bion) M^2 [cos(4pi X/Lambda) + 2 cos(2pi X/Lambda)]
#
# The minus sign is for MAGNETIC bions (center-stabilizing).
# M^2 = [Product_f (m_f L)]^2 is the mass factor.
#
# The quartic term (coefficient of |X|^4) is:
#   V^(4) = -(Lambda^4/g^4) exp(-S_bion) M^2 * (3/4)(2pi)^4 / Lambda^4 * |X|^4
#         = -(3/4)(2pi)^4 / g^4 * exp(-S_bion) * M^2 * |X|^4
#
# We want to match this to a Kahler-type quartic:
#   V^(4) = c_bion * Lambda^2 * |X|^4 / Lambda^4  (so c_bion is dimensionless)
# i.e.  c_bion / Lambda^2 * |X|^4 = -(3/4)(2pi)^4 / g^4 * exp(-S_bion) * M^2 * |X|^4
#
# => c_bion = -(3/4)(2pi)^4 Lambda^2 / g^4 * exp(-S_bion) * M^2
#
# BUT this has wrong dimensions (Lambda^2 is dimensionful).
#
# The resolution: c_bion is defined as the DIMENSIONLESS coefficient
# in K_bion = c_bion |X|^4/Lambda^2.  This is already dimensionless
# if |X| has dimension [mass] and Lambda has dimension [mass].
#
# The bion potential is V_bion with dimensions [mass]^4.
# To extract c_bion, we write:
#   V_bion^(4) = c_bion * |X|^4 / Lambda^2  (dimensions: [mass]^4 / [mass]^2 * [mass]^4 / [mass]^4... )
#
# No -- let me be more careful.
# K = |X|^2 + c_bion |X|^4/Lambda^2
# K_{XX*} = 1 + 4 c_bion |X|^2/Lambda^2
# V = |F_X|^2 / K_{XX*} + V_bion(X)
#
# The bion potential V_bion contributes DIRECTLY to the potential.
# The Kahler correction c_bion modifies the kinetic term.
# These are different physical effects.
#
# For the purpose of this problem, we identify the bion contribution
# to the DIRECT potential, not the Kahler correction.
# Then c_bion is defined by matching:
#   V_bion^(4) = c_bion * |X|^4 / Lambda^2
#
# This gives c_bion = V_bion^(4) * Lambda^2 / |X|^4
#                   = -(3/4)(2pi)^4 / g^4 * exp(-S_bion) * M^2 * Lambda^2
#
# Hmm, still dimensionful.  The issue is that V_bion^(4) has an overall
# mass scale.  Let me define things more carefully.
#
# Actually, the bion potential per unit 3-volume has dimension [mass]^4.
# The semiclassical result is:
#   V_bion = -C_0 * (1/(g^2 L))^3 * exp(-S_bion) * M^2 * [cos(...) + 2cos(...)]
#
# where 1/(g^2 L) ~ Lambda^3/g^2... This gets complicated.
# Let me just define c_bion as the COEFFICIENT in the potential expansion:
#
#   V_bion = V_0 + V_2 |X/Lambda|^2 + V_4 |X/Lambda|^4 + ...
#
# Then c_bion = V_4 / Lambda^4  (dimensionless ratio).
#
# From the expansion:
#   V_4 = A_0 * (3/4)(2pi)^4
# where A_0 = -(C_0/g^4) * Lambda^4 * exp(-S_bion) * M^2 is the overall amplitude.
#
# So c_bion = -(C_0/g^4) * (3/4)(2pi)^4 * exp(-S_bion) * M^2

# Substituting g^4 = 16 pi^2 alpha_s^2:
# c_bion = -(3/4)(2pi)^4 / (16 pi^2 alpha_s^2) * C_0 * exp(-S_bion) * M^2
#        = -3 pi^2 / alpha_s^2 * C_0 * exp(-S_bion) * M^2

def c_bion_func(a_s, Nf_mass_factor, C0=1.0, Nc_val=3):
    """Compute c_bion for given alpha_s and mass factor M^2."""
    s0 = 2*np.pi/(a_s * Nc_val)
    sb = 2*s0/Nc_val
    return -3*np.pi**2/a_s**2 * C0 * np.exp(-sb) * Nf_mass_factor

# With M^2 = 1 (N_f = 1, m = Lambda, or all masses equal Lambda):
print("c_bion = -3 pi^2 C_0 exp(-S_bion) M^2 / alpha_s^2")
print()
print("where S_bion = 4 pi / (alpha_s N_c^2) and M^2 = [Product_f(m_f/Lambda)]^2")
print()

# Let's define f(alpha_s) = 3 pi^2 / alpha_s^2 * exp(-4 pi/(9 alpha_s))
# for N_c = 3 with M^2 = 1.  We want f(alpha_s) = 1/12.

def f_Nf1(a_s):
    """Magnitude of c_bion for M^2 = 1, C_0 = 1."""
    sb = 4*np.pi/(9*a_s)
    return 3*np.pi**2/a_s**2 * np.exp(-sb)

# Analyze the function f(alpha_s)
print("Behavior of |c_bion| = 3 pi^2/alpha_s^2 * exp(-4pi/(9 alpha_s))")
print("as a function of alpha_s (with M^2 = 1, C_0 = 1):")
print()

a_s_range = np.logspace(-1.5, 2, 10000)
f_vals = np.array([f_Nf1(a) for a in a_s_range])
max_idx = np.argmax(f_vals)
max_as = a_s_range[max_idx]
max_f = f_vals[max_idx]

print(f"  Maximum: |c_bion|_max = {max_f:.6f} at alpha_s = {max_as:.4f}")
print()

# The maximum of f(a) = 3pi^2/a^2 exp(-4pi/(9a))
# df/da = 3pi^2 [-2/a^3 + 4pi/(9a^4)] exp(-4pi/(9a)) = 0
# => -2/a + 4pi/(9a^2) = 0
# => a = 2pi/9 = 0.6981

a_s_peak = 2*np.pi/9
f_peak = f_Nf1(a_s_peak)
print(f"  Analytical peak: alpha_s* = 2 pi/9 = {a_s_peak:.6f}")
print(f"  |c_bion| at peak = {f_peak:.6f}")
print()

print(f"  Since the max |c_bion| = {f_peak:.4f} is much larger than 1/12 = {1/12:.4f},")
print(f"  the equation |c_bion| = 1/12 has TWO solutions (one on each side of the peak).")
print()

# Find both roots
# Left root: small alpha_s, exponential suppression dominates
# Right root: large alpha_s, 1/alpha_s^2 suppression dominates

# Left root
try:
    a_s_left = brentq(lambda a: f_Nf1(a) - 1.0/12, 0.05, a_s_peak - 0.01)
    print(f"  Left solution:  alpha_s = {a_s_left:.6f}")
    sb_left = 4*np.pi/(9*a_s_left)
    print(f"    S_bion = {sb_left:.4f}, exp(-S_bion) = {np.exp(-sb_left):.6e}")
except Exception as e:
    a_s_left = None
    print(f"  Left solution: not found ({e})")

# Right root
try:
    a_s_right = brentq(lambda a: f_Nf1(a) - 1.0/12, a_s_peak + 0.01, 100.0)
    print(f"  Right solution: alpha_s = {a_s_right:.6f}")
    sb_right = 4*np.pi/(9*a_s_right)
    print(f"    S_bion = {sb_right:.4f}, exp(-S_bion) = {np.exp(-sb_right):.6e}")
except Exception as e:
    a_s_right = None
    print(f"  Right solution: not found ({e})")

print()

# Detailed scan
print(f"{'alpha_s':>10} | {'|c_bion|':>14} | {'ratio to 1/12':>14}")
print("-" * 44)
for a_s_scan in [0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5, a_s_peak, 1.0, 1.5, 2.0, 3.0, 5.0, 10.0]:
    cb = f_Nf1(a_s_scan)
    print(f"{a_s_scan:10.4f} | {cb:14.6f} | {cb/(1/12):14.6f}")
print()

# ===================================================================
# TASK 5: c_bion formula
# ===================================================================

print("=" * 72)
print("TASK 5: c_bion in Terms of Physical Parameters")
print("=" * 72)
print()

print("General formula:")
print()
print("  c_bion = - (3 pi^2 C_0 / alpha_s^2) exp(-4 pi/(N_c^2 alpha_s))")
print("           x [Product_f (m_f / Lambda)]^2")
print()
print("where:")
print("  S_0 = 2 pi / (alpha_s N_c)")
print("  S_bion = 2 S_0 / N_c = 4 pi / (alpha_s N_c^2)")
print("  C_0 = O(1) moduli-space prefactor (set to 1)")
print("  M^2 = [Product_f (m_f / Lambda)]^2 = mass zero-mode factor")
print()

# Numerical evaluation with physical masses
M2_dimless = np.prod([m/Lambda_QCD for m in m_quarks.values()])**2
c_bion_at_03 = c_bion_func(alpha_s, M2_dimless)

print(f"Numerical evaluation at alpha_s = {alpha_s}, N_f = 5:")
print()
print(f"  S_bion = {S_bion:.4f}")
print(f"  exp(-S_bion) = {np.exp(-S_bion):.6e}")
print(f"  3 pi^2 / alpha_s^2 = {3*np.pi**2/alpha_s**2:.4f}")
print(f"  M^2 = {M2_dimless:.6e}")
print(f"  c_bion = {c_bion_at_03:.6e}")
print(f"  Compare -1/12 = {-1.0/12:.6f}")
print()

# With M^2 = 1
c_bion_M1 = c_bion_func(alpha_s, 1.0)
print(f"With M^2 = 1 (all masses at Lambda scale):")
print(f"  c_bion = {c_bion_M1:.6f}")
print()

# ===================================================================
# TASK 6: Natural combination for c_bion = -1/12
# ===================================================================

print("=" * 72)
print("TASK 6: Is There a Natural c_bion = -1/12?")
print("=" * 72)
print()

print("The condition |c_bion| = 1/12 with M^2 = 1 has two solutions:")
print()

if a_s_left is not None:
    print(f"  (a) alpha_s = {a_s_left:.6f}  (weak-coupling side)")
    print(f"      S_bion = {4*np.pi/(9*a_s_left):.4f}")
    print(f"      This is semiclassically controlled: S_bion >> 1.")
    print()

if a_s_right is not None:
    print(f"  (b) alpha_s = {a_s_right:.6f}  (strong-coupling side)")
    print(f"      S_bion = {4*np.pi/(9*a_s_right):.4f}")
    print(f"      This is NOT semiclassically controlled.")
    print()

# Alternative: can we match with M^2 adjusted?
# Need: 3 pi^2 / alpha_s^2 * exp(-S_bion) * M^2 = 1/12
# At alpha_s = 0.3: 3 pi^2/0.09 * exp(-4.654) * M^2 = 1/12
# => M^2 = 1/(12 * 328.99 * 0.00952) = 1/(12 * 3.132) = 1/37.59 = 0.0266

M2_needed_03 = 1.0 / (12 * 3*np.pi**2/alpha_s**2 * np.exp(-S_bion))
print(f"At the given alpha_s = {alpha_s}:")
print(f"  M^2 needed for c_bion = -1/12: {M2_needed_03:.6f}")
print(f"  Actual M^2 (5 flavors): {M2_dimless:.6e}")
print(f"  Ratio: {M2_dimless/M2_needed_03:.6e}")
print()

# For N_f = 1 with mass m:
# M^2 = (m/Lambda)^2
# m/Lambda needed:
m_needed = np.sqrt(M2_needed_03)
print(f"  For N_f = 1: m/Lambda needed = {m_needed:.6f}")
print(f"    m = {m_needed * Lambda_QCD:.1f} MeV")
print()

# Check if this is near any physical quark mass
print("  Compare with quark masses:")
for name, m in m_quarks.items():
    ratio = m / (m_needed * Lambda_QCD)
    print(f"    m_{name} = {m:.1f} MeV, ratio = {ratio:.4f}")
print()

# The interesting case: alpha_s at the peak
print(f"At the peak alpha_s = 2 pi/9 = {a_s_peak:.4f}:")
M2_needed_peak = 1.0 / (12 * f_peak)
print(f"  |c_bion|_max = {f_peak:.4f}")
print(f"  M^2 needed = {M2_needed_peak:.6e}")
print(f"  For N_f = 1: m/Lambda = {np.sqrt(M2_needed_peak):.6f}")
print()

# KEY: Is there a NATURAL value of alpha_s where c_bion = -1/12
# without fine-tuning?
print("=" * 72)
print("KEY RESULTS")
print("=" * 72)
print()

print("The function |c_bion(alpha_s)| = 3 pi^2/alpha_s^2 * exp(-4 pi/(9 alpha_s))")
print("has a single maximum at alpha_s = 2 pi/9.")
print()
print(f"  Peak value: {f_peak:.4f}")
print(f"  Target:     {1.0/12:.4f}")
print()

if f_peak > 1.0/12:
    print(f"  The peak exceeds 1/12 by a factor of {f_peak * 12:.2f}.")
    print(f"  Therefore c_bion = -1/12 IS achievable for two values of alpha_s")
    print(f"  (one on each side of the peak).")
    print()

    if a_s_left is not None:
        g2_left = 4*np.pi*a_s_left
        S0_left = 2*np.pi/(a_s_left*Nc)
        print(f"  Solution (a): alpha_s = {a_s_left:.6f}")
        print(f"    g^2 = {g2_left:.4f}")
        print(f"    S_0 = {S0_left:.4f}")
        print(f"    S_bion = {4*np.pi/(9*a_s_left):.4f}")
        print(f"    Semiclassical reliability: {'YES' if 4*np.pi/(9*a_s_left) > 2 else 'marginal'}")
        print()

    if a_s_right is not None:
        g2_right = 4*np.pi*a_s_right
        S0_right = 2*np.pi/(a_s_right*Nc)
        print(f"  Solution (b): alpha_s = {a_s_right:.6f}")
        print(f"    g^2 = {g2_right:.4f}")
        print(f"    S_0 = {S0_right:.4f}")
        print(f"    S_bion = {4*np.pi/(9*a_s_right):.4f}")
        print(f"    Semiclassical reliability: {'YES' if 4*np.pi/(9*a_s_right) > 2 else 'NO (strong coupling)'}")
        print()

# With physical quark masses (N_f = 5):
print("With N_f = 5 physical quark masses:")
f_vals_5 = np.array([f_Nf1(a) * M2_dimless for a in a_s_range])
max_idx_5 = np.argmax(f_vals_5)
max_f_5 = f_vals_5[max_idx_5]
max_as_5 = a_s_range[max_idx_5]

print(f"  Max |c_bion| = {max_f_5:.6e} at alpha_s = {max_as_5:.4f}")
print(f"  This is {max_f_5/(1.0/12):.2e} times 1/12.")

if max_f_5 < 1.0/12:
    print(f"  With physical masses, |c_bion| NEVER reaches 1/12.")
    print(f"  The light-quark mass suppression (m_u m_d << Lambda) kills it.")
else:
    print(f"  Solution exists (barely).")
print()

# ===================================================================
# Physical interpretation
# ===================================================================

print("=" * 72)
print("PHYSICAL INTERPRETATION")
print("=" * 72)
print()
print("1. For M^2 = 1 (all quark masses ~ Lambda), the bion Kahler correction")
print("   c_bion = -1/12 is achieved at two values of alpha_s.")
if a_s_left is not None:
    print(f"   The semiclassically valid solution is alpha_s = {a_s_left:.4f}.")
print()
print("2. With 5 physical flavors, the light-quark mass factor M^2 ~ 4e-6")
print("   suppresses c_bion well below 1/12 at any coupling.")
print()
print("3. In the ISS model with N_f = N_c + 1 = 4 heavy flavors")
print("   (all at the Lambda scale), the bion correction IS large enough")
print("   to produce c_bion = -1/12.")
print()
print("4. The condition c_bion = -1/12 translates to:")
print("   exp(-4 pi/(9 alpha_s)) = alpha_s^2 / (36 pi^2)")
print("   This is a transcendental equation with solutions at")
if a_s_left is not None:
    print(f"   alpha_s = {a_s_left:.4f} and alpha_s = {a_s_right:.4f}.")
print()

# ===================================================================
# COMPLETE SUMMARY
# ===================================================================

print("=" * 72)
print("COMPLETE SUMMARY TABLE")
print("=" * 72)
print()
print("Task 1: 3 bion types (6 with orientation) from SU(3) extended Dynkin triangle")
print()
print("Task 2: A_ab propto [Product_f(m_f L)]^2, same for all bion types")
print()
print("Task 3: V_bion = A e^{-S_b} [cos(4pi X/L) + 2 cos(2pi X/L)]")
print(f"        Quartic: (3/4)(2pi)^4 = {(3./4)*(2*np.pi)**4:.2f}")
print()
print("Task 4: c_4 = -1/12 requires 3pi^2 e^{-S_b}/alpha_s^2 = 1/12")
if a_s_left:
    print(f"        Solutions: alpha_s = {a_s_left:.4f} or {a_s_right:.4f} (for M^2 = 1)")
print()
print("Task 5: c_bion = -3 pi^2 C_0 exp(-4pi/(9 alpha_s)) M^2 / alpha_s^2")
print()
print("Task 6: With M^2 = 1, c_bion = -1/12 at alpha_s =", end=" ")
if a_s_left:
    print(f"{a_s_left:.4f} (semiclassical)")
else:
    print("(no semiclassical solution)")
print(f"        With 5 physical flavors: NOT achievable (mass suppression)")
print(f"        With N_f = 4, m ~ Lambda (ISS): achievable")
