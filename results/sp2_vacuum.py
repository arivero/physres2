"""
Sp(2) SQCD with N_f = 3: Complete vacuum analysis.

Convention: Sp(2) = rank-2 symplectic group, i.e., Sp(4) as a matrix group
(4x4 antisymmetric symplectic form). Fundamentals are 4-dimensional.
N_f = 3 means 3 chiral superfields Q^i (i = 1,...,3) in the fundamental (4-dim).
The 2*N_f = 6 flavor indices arise from the pseudo-real structure:
the combined index runs over i = 1,...,6 (= 2*N_f) in the antisymmetric meson.

This is the s-confining case: N_f = N_c + 1 for Sp(N_c) with N_c = 2.

Mesons: M^{ij} = Q^i J Q^j (antisymmetric 6x6, J = symplectic form).
Pfaffian: Pf(M) = product of Pfaffian eigenvalues = lambda_1 lambda_2 lambda_3.
"""

import numpy as np
from scipy.optimize import fsolve

print("=" * 78)
print("  Sp(2) SQCD with N_f = 3: COMPLETE ANALYSIS")
print("  Convention: Sp(2) = rank-2 symplectic, fundamentals are 4-dim")
print("=" * 78)


# ============================================================================
# PART 1: CONFINED SPECTRUM
# ============================================================================

print("\n" + "=" * 78)
print("PART 1: CONFINED SPECTRUM")
print("=" * 78)

N_c = 2       # rank of Sp(N_c), gauge group = Sp(2*N_c) = Sp(4)
N_f = 3       # fundamental flavors
dim_fund = 2 * N_c   # fundamental of Sp(2N_c) is 2N_c-dimensional = 4
n_half = 2 * N_f     # 2N_f = 6 flavor indices for the antisymmetric meson

n_mesons = n_half * (n_half - 1) // 2
print(f"""
Gauge group: Sp(2*{N_c}) = Sp({2*N_c}) (rank {N_c})
Fundamental representation: dimension {dim_fund}
Number of flavors: N_f = {N_f}
Number of flavor half-indices: 2*N_f = {n_half}

Mesons: M^{{ij}} = Q^i_a J^{{ab}} Q^j_b  (contracted over color a,b = 1,...,{dim_fund})

Since J is antisymmetric, M^{{ij}} = -M^{{ji}}: the meson matrix is ANTISYMMETRIC.
The 6x6 antisymmetric matrix M has:
  Independent components = C(6,2) = {n_mesons}

One-loop beta function coefficient:
  b_0 = 3(N_c + 1) - N_f = 3*{N_c+1} - {N_f} = {3*(N_c+1) - N_f}
  Asymptotically free: YES

Confinement: N_f = {N_f} = N_c + 1 = {N_c + 1}
  -> This is the s-CONFINING case (Intriligator-Pouliot 1995).

Quantum-modified constraint:
  Pf(M) = Lambda^{{2(N_c+1)}} = Lambda^{2*(N_c+1)} = Lambda^6

  where Pf(M) is the Pfaffian of the 6x6 antisymmetric matrix M.
  For a 2n x 2n antisymmetric matrix: Pf(M)^2 = det(M).
  Here 2n = 6, so Pf(M) is a degree-3 polynomial in the M^{{ij}}.

Confined superpotential (with Lagrange multiplier X):
  W_conf = X * (Pf(M) - Lambda^6)

  This is the exact low-energy description. X enforces the quantum
  constraint Pf(M) = Lambda^6.
""")


# ============================================================================
# PART 2: O'RAIFEARTAIGH VACUUM
# ============================================================================

print("=" * 78)
print("PART 2: O'RAIFEARTAIGH VACUUM")
print("=" * 78)

print("""
Tree-level superpotential deformation:
  W_tree = g Tr(M * A) + m Tr(M)

where A is a fixed antisymmetric 6x6 background (spurion) matrix and m is a
mass parameter. The total superpotential is:

  W_total = X (Pf(M) - Lambda^6) + g Tr(M A) + m Tr(M J)

where J is the symplectic form (used as a "trace" for antisymmetric matrices:
Tr(M) = sum_{i<j} M^{ij} J_{ij}).

F-TERM EQUATIONS:

  F_X = Pf(M) - Lambda^6 = 0
  F_{M^{ij}} = X * (d Pf(M) / d M^{ij}) + g A_{ij} + m J_{ij} = 0

For a 6x6 antisymmetric matrix in Pfaffian-eigenvalue form
  M = diag(lambda_1 J_2, lambda_2 J_2, lambda_3 J_2)

where J_2 = [[0,1],[-1,0]], the Pfaffian is:
  Pf(M) = lambda_1 * lambda_2 * lambda_3

and its derivatives give:
  d Pf(M) / d lambda_k = Pf(M) / lambda_k = lambda_1 lambda_2 lambda_3 / lambda_k

With A = diag(a_1 J_2, a_2 J_2, a_3 J_2), the F-term equations reduce to:

  (I)   lambda_1 lambda_2 lambda_3 = Lambda^6
  (II)  X Lambda^6 / lambda_k + g a_k + m = 0,   k = 1,2,3

From (II):
  lambda_k = -X Lambda^6 / (g a_k + m)

Substituting into (I):
  (-X Lambda^6)^3 / [(g a_1 + m)(g a_2 + m)(g a_3 + m)] = Lambda^6
  => X^3 = -Product_k (g a_k + m) / Lambda^{12}

This ALWAYS has a solution (cubic equation with at least one real root).
=> The naive construction preserves SUSY.

SUSY BREAKING requires additional structure:
the O'Raifeartaigh mechanism with incompatible F-terms.
""")

print("""
--- THE WORKING O'RAIFEARTAIGH CONSTRUCTION ---

Embed a standard O'Raifeartaigh model into the meson sector:

  W = f X + m M_12 M_34 + g X M_12^2 + (eps/Lambda^3) Pf(M)

where M_12 and M_34 are specific meson components (a "complementary pair"),
X is a singlet, f is the SUSY-breaking F-term, and the Pfaffian term
encodes the nonperturbative dynamics.

At eps = 0 (ignoring the Pfaffian), this is the standard O'Raifeartaigh model:
  F_X    = f + g M_12^2    => wants M_12 != 0
  F_{M34} = m M_12          => wants M_12 = 0
  CONFLICT => SUSY BROKEN, F_X = f != 0.

At eps != 0, a distant SUSY-preserving vacuum appears at M ~ Lambda^3/eps.
The metastable SUSY-breaking vacuum near the origin persists, with
tunneling rate Gamma ~ exp(-(Lambda^3/eps)^4).

The CONDITION for O'Raifeartaigh SUSY breaking:
  1. f, g, m all nonzero
  2. The cubic and linear terms in M_12 are incompatible
  3. The pseudo-modulus X has a flat direction at tree level (lifted by CW)
""")


# ============================================================================
# PART 3: PFAFFIAN EIGENVALUES
# ============================================================================

print("=" * 78)
print("PART 3: PFAFFIAN EIGENVALUES")
print("=" * 78)

print("""
A 6x6 antisymmetric matrix M can be brought to canonical (Darboux/normal) form
by a congruence transformation U^T M U (where U is in SO(6) or Sp(6)):

  M_canonical = diag(lambda_1 J_2, lambda_2 J_2, lambda_3 J_2)

       = [[  0,       lambda_1,  0,        0,        0,       0       ],
          [-lambda_1, 0,         0,        0,        0,       0       ],
          [  0,       0,         0,        lambda_2, 0,       0       ],
          [  0,       0,        -lambda_2, 0,        0,       0       ],
          [  0,       0,         0,        0,        0,       lambda_3],
          [  0,       0,         0,        0,       -lambda_3, 0      ]]

The lambda_k >= 0 are the "Pfaffian eigenvalues" (or Williamson normal form
eigenvalues). They are related to the ordinary eigenvalues of M by:
  eigenvalues of M = {+i lambda_1, -i lambda_1, +i lambda_2, -i lambda_2,
                      +i lambda_3, -i lambda_3}

PFAFFIAN:
  Pf(M) = lambda_1 * lambda_2 * lambda_3

QUANTUM CONSTRAINT:
  lambda_1 * lambda_2 * lambda_3 = Lambda^6

VACUUM in terms of Pfaffian eigenvalues:

At the O'Raifeartaigh vacuum (in the sector where M is block-diagonal and
the background A is also block-diagonal):

  From Part 2: lambda_k = -X Lambda^6 / (g a_k + m)

  With the constraint: X^3 = -Product_k(g a_k + m) / Lambda^{12}

  If a_k = 0 for all k (no background splitting):
    All three eigenvalues are equal: lambda_1 = lambda_2 = lambda_3 = Lambda^2
    (since Lambda^6 = Lambda^2 * Lambda^2 * Lambda^2)

  With nonzero a_k, the eigenvalues split according to the background.
""")


# ============================================================================
# PART 4: KOIDE CONDITION
# ============================================================================

print("=" * 78)
print("PART 4: KOIDE CONDITION FROM gv/m = sqrt(3)")
print("=" * 78)

print("""
At the metastable O'Raifeartaigh vacuum (M = 0, X = v = pseudo-modulus VEV):

The fermion mass matrix W_{IJ} in the (M_12, M_34) sector is:

  W_F = [[2gv,  m],
         [m,    0]]

Eigenvalues: m_pm = gv +/- sqrt(g^2 v^2 + m^2)

Setting t = gv/m:
  m_+ = m(t + sqrt(t^2 + 1))
  m_- = m(sqrt(t^2 + 1) - t)

Including the Goldstino (X fermion, mass = 0) and all flat directions,
the Koide triple is (0, m_-, m_+).
""")

# Compute Q for gv/m = sqrt(3)
t_koide = np.sqrt(3.0)
m_plus  = t_koide + np.sqrt(t_koide**2 + 1.0)    # in units of m
m_minus = np.sqrt(t_koide**2 + 1.0) - t_koide     # in units of m

def koide_Q_triple(m1, m2, m3):
    """Koide quality factor Q = sum(m_i) / (sum(sqrt(m_i)))^2."""
    masses = np.array([m1, m2, m3], dtype=float)
    sqrts  = np.sqrt(masses)
    return np.sum(masses) / np.sum(sqrts)**2

Q_seed = koide_Q_triple(0.0, m_minus, m_plus)

print(f"For t = gv/m = sqrt(3) = {t_koide:.12f}:")
print(f"  sqrt(t^2 + 1) = sqrt(4) = {np.sqrt(t_koide**2 + 1):.12f}")
print(f"  m_+ = (2 + sqrt(3)) m = {m_plus:.12f} m")
print(f"  m_- = (2 - sqrt(3)) m = {m_minus:.12f} m")
print(f"  m_0 = 0 (Goldstino)")
print(f"  m_+ * m_- = {m_plus * m_minus:.15f} (should be 1)")
print()
print(f"Koide quality factor Q(0, m_-, m_+):")
print(f"  Sum(m_i) = 0 + {m_minus:.10f} + {m_plus:.10f} = {m_minus + m_plus:.10f}")

a = np.sqrt(m_minus)
b = np.sqrt(m_plus)
print(f"  sqrt(m_-) = {a:.10f}")
print(f"  sqrt(m_+) = {b:.10f}")
print(f"  sqrt(m_-) * sqrt(m_+) = {a*b:.15f}  (= sqrt(m_+ m_-) = 1)")
print(f"  (0 + sqrt(m_-) + sqrt(m_+))^2 = {(a + b)**2:.10f}")
print(f"  = m_- + m_+ + 2*sqrt(m_- m_+) = {m_minus + m_plus:.6f} + 2 = {m_minus + m_plus + 2:.6f}")
print(f"  Q = {m_minus + m_plus:.10f} / {(a + b)**2:.10f} = {Q_seed:.15f}")
print(f"  2/3 = {2.0/3.0:.15f}")
print(f"  |Q - 2/3| = {abs(Q_seed - 2.0/3.0):.2e}")
print()
print("ANALYTIC PROOF:")
print("  Sum(m_i) = (2-sqrt3) + (2+sqrt3) = 4")
print("  (Sum sqrt(m_i))^2 = (sqrt(2-sqrt3) + sqrt(2+sqrt3))^2")
print("                    = (2-sqrt3) + (2+sqrt3) + 2*sqrt((2-sqrt3)(2+sqrt3))")
print("                    = 4 + 2*sqrt(4-3) = 4 + 2 = 6")
print("  Q = 4/6 = 2/3 EXACTLY.  QED.")


# ============================================================================
# PART 5: BLOOM ROTATION
# ============================================================================

print("\n\n" + "=" * 78)
print("PART 5: BLOOM ROTATION")
print("=" * 78)

SQRT2 = np.sqrt(2.0)
TWO_PI_OVER_3 = 2.0 * np.pi / 3.0
DELTA_SEED = 3.0 * np.pi / 4.0   # = 135 deg; forces z_0 = 0

def koide_parametrization(M0, delta):
    """Compute signed roots z_k and masses m_k from Koide parametrization."""
    z = np.array([1.0 + SQRT2 * np.cos(delta + TWO_PI_OVER_3 * k) for k in range(3)])
    masses = M0 * z**2
    return masses, z

def fit_M0_delta(signed_roots):
    """Given signed roots z_k, extract (M0, delta)."""
    s = np.array(signed_roots, dtype=float)
    sqrtM0 = s.sum() / 3.0
    M0 = sqrtM0**2
    phi = np.array([TWO_PI_OVER_3 * k for k in range(3)])
    A = np.dot(s, np.cos(phi))
    B = np.dot(s, -np.sin(phi))
    delta = np.arctan2(B, A) % (2.0 * np.pi)
    return M0, delta

print("""
--- Part 5a: What rotation preserves Q = 2/3? ---

The Koide parametrization writes three signed square roots as:
  z_k(delta) = 1 + sqrt(2) cos(delta + 2*pi*k/3),   k = 0, 1, 2

with m_k = M_0 * z_k^2. Then:

  Sum(z_k) = 3 + sqrt(2) * Sum_k cos(delta + 2*pi*k/3) = 3 + 0 = 3
  Sum(z_k^2) = 3 + 2*Sum_k cos^2(delta + 2*pi*k/3) + 2*sqrt(2)*Sum_k cos(...)
             = 3 + 2*(3/2) + 0 = 6

  Q = M_0 * Sum(z_k^2) / (M_0 * (Sum z_k)^2)  [only valid when all z_k > 0]
    = 6/9 ... NO WAIT.

  Q = Sum(m_k) / (Sum sqrt(m_k))^2
    = M_0 Sum(z_k^2) / (sqrt(M_0) Sum(z_k))^2   [if all z_k > 0]
    = Sum(z_k^2) / (Sum(z_k))^2
    = 6 / 9 ... that gives 2/3!

RESULT: Q = 2/3 is an ALGEBRAIC IDENTITY of the Koide parametrization.
  ANY pure delta-rotation at fixed M_0 preserves Q = 2/3 exactly.
  (Valid when all z_k > 0, or when signed roots are used consistently.)

  This is the fundamental property: the bloom is a Q-preserving rotation
  of delta at fixed M_0.
""")

# Numerical verification
print("Numerical verification:")
for delta_test in [0.5, 1.0, 1.5, 2.0, 2.3, np.pi/3, np.pi/2]:
    ms, zs = koide_parametrization(1.0, delta_test)
    if all(zs > 0):
        Q = koide_Q_triple(ms[0], ms[1], ms[2])
        print(f"  delta = {np.degrees(delta_test):7.2f} deg: z = ({zs[0]:.4f}, {zs[1]:.4f}, {zs[2]:.4f}), Q = {Q:.15f}")
    else:
        print(f"  delta = {np.degrees(delta_test):7.2f} deg: z = ({zs[0]:.4f}, {zs[1]:.4f}, {zs[2]:.4f})  [boundary/negative z]")

print(f"""
--- Part 5b: Delta for (m_e, m_mu, m_tau) ---

The seed spectrum is at delta_seed = 3*pi/4 = 135 deg, where z_0 = 0.
The physical lepton masses are: m_e = 0.511 MeV, m_mu = 105.658 MeV, m_tau = 1776.86 MeV.
""")

m_e   = 0.51100
m_mu  = 105.658
m_tau = 1776.86

# Fit Koide parameters
sr_lepton = [np.sqrt(m_e), np.sqrt(m_mu), np.sqrt(m_tau)]
M0_lep, delta_lep = fit_M0_delta(sr_lepton)

# Bloom rotation from seed
delta_bloom_lep = delta_lep - DELTA_SEED
delta_bloom_lep = ((delta_bloom_lep + np.pi) % (2.0 * np.pi)) - np.pi

print(f"Signed roots: ({sr_lepton[0]:.6f}, {sr_lepton[1]:.6f}, {sr_lepton[2]:.6f})")
print(f"Fitted M_0 = {M0_lep:.6f} MeV")
print(f"Fitted delta = {delta_lep:.8f} rad = {np.degrees(delta_lep):.6f} deg")
print(f"Seed delta = {np.degrees(DELTA_SEED):.4f} deg")
print(f"Bloom rotation: Delta_delta = {delta_bloom_lep:.8f} rad = {np.degrees(delta_bloom_lep):.6f} deg")

# Verify reconstruction
ms_recon, zs_recon = koide_parametrization(M0_lep, delta_lep)
print(f"\nReconstruction check:")
print(f"  Reconstructed: ({ms_recon[0]:.4f}, {ms_recon[1]:.4f}, {ms_recon[2]:.4f}) MeV")
print(f"  Physical:      ({m_e:.4f}, {m_mu:.4f}, {m_tau:.4f}) MeV")
Q_lep = koide_Q_triple(m_e, m_mu, m_tau)
print(f"  Q_lepton = {Q_lep:.10f}  (2/3 = {2/3:.10f}, deviation = {(Q_lep - 2/3)/(2/3)*100:.4f}%)")

print(f"""
--- Part 5c: M_0 for the lepton sector ---

From the fit:
  M_0 = (sum(sqrt(m_k))/3)^2 = {M0_lep:.6f} MeV
  sqrt(M_0) = {np.sqrt(M0_lep):.6f} MeV^(1/2)

The seed spectrum at delta = 3*pi/4 has mass ratios (0, 2-sqrt(3), 2+sqrt(3)):
  Seed masses: (0, {M0_lep * (2 - np.sqrt(3)):.4f}, {M0_lep * (2 + np.sqrt(3)):.4f}) MeV

The bloom rotates from delta = 135 deg to delta = {np.degrees(delta_lep):.4f} deg,
redistributing mass while preserving Q = 2/3 and M_0.
""")

# Also do quark sectors
m_s = 93.4
m_c = 1270.0
m_b = 4180.0
m_t = 172760.0

sr_scb = [-np.sqrt(m_s), np.sqrt(m_c), np.sqrt(m_b)]
M0_scb, delta_scb = fit_M0_delta(sr_scb)
bloom_scb = ((delta_scb - DELTA_SEED + np.pi) % (2*np.pi)) - np.pi

sr_cbt = [np.sqrt(m_c), np.sqrt(m_b), np.sqrt(m_t)]
M0_cbt, delta_cbt = fit_M0_delta(sr_cbt)
bloom_cbt = ((delta_cbt - DELTA_SEED + np.pi) % (2*np.pi)) - np.pi

print("Bloom parameters for all triples:")
print(f"{'Triple':<16} {'M_0 (MeV)':>12} {'delta (deg)':>12} {'bloom (deg)':>12} {'Q':>10}")
print("-" * 66)
for name, M0, delta, bloom, sr in [
    ("(e, mu, tau)", M0_lep, delta_lep, delta_bloom_lep, sr_lepton),
    ("(-s, c, b)", M0_scb, delta_scb, bloom_scb, sr_scb),
    ("(c, b, t)", M0_cbt, delta_cbt, bloom_cbt, sr_cbt),
]:
    m_arr = np.array(sr)**2
    Q_val = koide_Q_triple(*m_arr)
    print(f"{name:<16} {M0:>12.4f} {np.degrees(delta):>12.4f} {np.degrees(bloom):>12.4f} {Q_val:>10.7f}")

# seed M_0 for the lepton sector
print(f"\nLepton sector: M_0 = {M0_lep:.4f} MeV")
print(f"  This is the common mass scale. The seed has masses")
print(f"  (0, (2-sqrt3)*M_0, (2+sqrt3)*M_0) = (0, {M0_lep*(2-np.sqrt(3)):.4f}, {M0_lep*(2+np.sqrt(3)):.4f}) MeV")


# ============================================================================
# PART 6: BION CORRECTIONS — Sp(2) vs SU(3)
# ============================================================================

print("\n\n" + "=" * 78)
print("PART 6: BION CORRECTIONS — Sp(2) vs SU(3)")
print("=" * 78)

print("""
ROOT SYSTEM AND BION STRUCTURE:

SU(3) (A_2):
  Extended Dynkin diagram: triangle (all 3 pairs adjacent)
  Center: Z_3
  Bion types: 3 (up to conjugation), related by Z_3 symmetry
  Bion potential: V ~ cos(2phi) + 2*cos(phi)  [Z_3 invariant]

Sp(2) = USp(4) (C_2):
  Root system: C_2 with simple roots alpha_1 (short), alpha_2 (long)
  Extended Dynkin diagram: chain  alpha_0 =>= alpha_1 =<= alpha_2
    (NOT a cycle — alpha_0 and alpha_2 are NOT adjacent)
  Center: Z_2
  Bion types: 2 (up to conjugation)
    B_{01}: monopoles on (alpha_0, alpha_1)
    B_{12}: monopoles on (alpha_1, alpha_2)
  Bion potential: V ~ A_1 cos(p_1 phi) + A_2 cos(p_2 phi)  [Z_2 invariant]

  Dual Coxeter number: h^v = 3 (comarks: a_0^v = a_1^v = a_2^v = 1)
  All monopole-instanton actions are equal: S_k = 2*pi/(3*alpha_s)

KEY STRUCTURAL DIFFERENCE:
  SU(3): 3 bion types, cyclic Z_3 symmetry, triangle Dynkin diagram
  Sp(2): 2 bion types, Z_2 symmetry, chain Dynkin diagram

  The number of bion types affects the angular structure of the potential
  in the Koide delta-space:
    SU(3): V ~ cos(3*delta)  [from Z_3]  -> 3 minima per period
    Sp(2): V ~ cos(2*delta)  [from Z_2]  -> 2 minima per period
""")

print("""
v_0-DOUBLING RATIO:

The v_0-doubling is defined as:
  v_0(full) / v_0(seed)

where:
  v_0(seed) = (sqrt(m_s) + sqrt(m_c)) / 3   [seed triple (s, c, 0)]
  v_0(full) = (-sqrt(m_s) + sqrt(m_c) + sqrt(m_b)) / 3  [full triple (-s, c, b)]
""")

# Compute v_0 doubling for quarks
v0_seed = (np.sqrt(m_s) + np.sqrt(m_c)) / 3.0
v0_full = (-np.sqrt(m_s) + np.sqrt(m_c) + np.sqrt(m_b)) / 3.0
ratio_v0 = v0_full / v0_seed

print(f"Numerical values:")
print(f"  v_0(seed) = ({np.sqrt(m_s):.4f} + {np.sqrt(m_c):.4f}) / 3 = {v0_seed:.6f} MeV^(1/2)")
print(f"  v_0(full) = (-{np.sqrt(m_s):.4f} + {np.sqrt(m_c):.4f} + {np.sqrt(m_b):.4f}) / 3 = {v0_full:.6f} MeV^(1/2)")
print(f"  Ratio = {ratio_v0:.6f}")
print()
print(f"  The prediction: sqrt(m_b) = 3*sqrt(m_s) + sqrt(m_c)")
pred_mb = (3*np.sqrt(m_s) + np.sqrt(m_c))**2
print(f"  Predicted m_b = {pred_mb:.2f} MeV (PDG: 4180 +/- 30 MeV)")
print(f"  Deviation: {abs(pred_mb - 4180)/4180 * 100:.3f}%")

print(f"""
ANSWER: The v_0 ratio = {ratio_v0:.4f} is a KINEMATIC result.

The factor of 2 arises purely from the Koide algebra:
  v_0(full) = (sum of signed roots in full triple) / 3
  v_0(seed) = (sum of signed roots in seed pair) / 3

When the bloom converts (s, c, 0) -> (-s, c, b) by:
  1. Flipping the sign of sqrt(m_s)  (changes sum by -2*sqrt(m_s))
  2. Replacing 0 by sqrt(m_b)        (changes sum by +sqrt(m_b))

For the sum to double: sqrt(m_b) - 2*sqrt(m_s) = sqrt(m_s) + sqrt(m_c)
=> sqrt(m_b) = 3*sqrt(m_s) + sqrt(m_c)  [exactly the v_0-doubling prediction]

This is INDEPENDENT of the gauge group. It depends only on the mass values.

Therefore: v_0(full)/v_0(seed) = {ratio_v0:.4f} for BOTH SU(3) and Sp(2).

The gauge group affects HOW the bloom happens dynamically:
  SU(3): cos(3*delta) potential from Z_3 center -> 3 bion types
  Sp(2): cos(2*delta) potential from Z_2 center -> 2 bion types

  The Sp(2) potential has minima at delta = 0 and delta = pi.
  The seed at delta = 3*pi/4 = 135 deg is:
    - 45 deg from the Sp(2) minimum at pi (= 180 deg)
    - 75 deg from the SU(3) minimum at 2*pi/3 (= 120 deg)

  So the Sp(2) bloom rotation is SHORTER (45 deg vs 75 deg),
  but the kinematic v_0 ratio is the same: {ratio_v0:.4f}.
""")

# Additional: lepton v0 ratio for comparison
v0_seed_lep = (np.sqrt(m_mu) + np.sqrt(m_tau)) / 3.0
v0_full_lep = (np.sqrt(m_e) + np.sqrt(m_mu) + np.sqrt(m_tau)) / 3.0
ratio_lep = v0_full_lep / v0_seed_lep
print(f"For comparison, the lepton v_0 ratio (seed (mu, tau, 0) -> full (e, mu, tau)):")
print(f"  v_0(full)/v_0(seed) = {ratio_lep:.6f}")
print(f"  This is NOT 2; the doubling is specific to the quark sector where")
print(f"  the sign flip (-s) contributes decisively.")


# ============================================================================
# SUMMARY
# ============================================================================

print("\n\n" + "=" * 78)
print("COMPLETE SUMMARY OF RESULTS")
print("=" * 78)

print(f"""
1. CONFINED SPECTRUM:
   Sp(2) = USp(4) with N_f = 3 flavors.
   Mesons M^{{ij}} (antisymmetric 6x6): {n_mesons} independent components.
   Quantum constraint: Pf(M) = Lambda^6.
   Confined superpotential: W = X(Pf(M) - Lambda^6).

2. O'RAIFEARTAIGH VACUUM:
   W = f X + m M_12 M_34 + g X M_12^2 + (eps/Lambda^3) Pf(M)
   F-terms: F_X = f (SUSY broken), F_M34 = m M_12 (conflicts with F_X).
   Metastable vacuum at origin; SUSY vacuum at M ~ Lambda^3/eps.

3. PFAFFIAN EIGENVALUES:
   M -> diag(lambda_1 J_2, lambda_2 J_2, lambda_3 J_2)
   Pf(M) = lambda_1 lambda_2 lambda_3 = Lambda^6
   At O'R vacuum: lambda_k = -X Lambda^6 / (g a_k + m)

4. KOIDE CONDITION:
   At the metastable vacuum with gv/m = sqrt(3):
   Fermion spectrum: (0, (2-sqrt(3))m, (2+sqrt(3))m)
   Q = 4/6 = 2/3 EXACTLY.
   This is the unique value of t = gv/m giving Q = 2/3.

5. BLOOM ROTATION:
   a) ANY delta-rotation at fixed M_0 preserves Q = 2/3 (algebraic identity).
   b) Lepton delta = {np.degrees(delta_lep):.4f} deg, bloom from seed = {np.degrees(delta_bloom_lep):.4f} deg.
   c) Lepton M_0 = {M0_lep:.4f} MeV.

6. BION CORRECTIONS (Sp(2) vs SU(3)):
   SU(3): 3 bion types (Z_3), cos(3*delta) potential, triangle Dynkin diagram.
   Sp(2): 2 bion types (Z_2), cos(2*delta) potential, chain Dynkin diagram.
   v_0(full)/v_0(seed) = {ratio_v0:.4f} for BOTH groups (kinematic, not dynamic).
   The gauge group changes the bloom mechanism but not the v_0 ratio.
""")

print("Done.")
