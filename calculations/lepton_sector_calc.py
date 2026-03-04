#!/usr/bin/env python3
"""
SU(2) SQCD lepton sector: O'Raifeartaigh model for s-confining theory.

SU(2) with N_f = 3 fundamentals Psi^i_a (i=1,2,3 flavor, a=1,2 color).
This is the s-confining case N_f = N_c + 1 = 3.

Confined degrees of freedom: antisymmetric mesons L^{ij} = eps^{ab} Psi^i_a Psi^j_b
  Three independent: L^{12}, L^{13}, L^{23}

Superpotential:
  W = f X_L + m L^{12} L^{23} + g X_L (L^{12})^2 + eps L^{12} L^{13} L^{23} / Lambda_L^5

Parts (a)-(f): F-terms, fermion mass matrix, spectrum, Koide verification,
mediation, confinement scale, goldstino identification.
"""

import numpy as np
from numpy import sqrt, array, pi

# Physical constants
m_e   = 0.51100    # MeV
m_mu  = 105.658    # MeV
m_tau = 1776.86    # MeV
M0_lepton = m_e + m_mu + m_tau  # total lepton mass

print("=" * 72)
print("PART (a): F-TERM EQUATIONS")
print("=" * 72)

print("""
Superpotential:
  W = f X_L + m L^{12} L^{23} + g X_L (L^{12})^2 + eps L^{12} L^{13} L^{23} / Lambda_L^5

Fields: X_L, L^{12}, L^{13}, L^{23}
Abbreviations: L1 = L^{12}, L2 = L^{13}, L3 = L^{23}

F-term equations (dW/dPhi_I = 0):

(1) F_{X_L} = dW/dX_L = f + g (L1)^2

(2) F_{L1}  = dW/dL1  = m L3 + 2g X_L L1 + eps (L2 L3) / Lambda_L^5

(3) F_{L2}  = dW/dL2  = eps L1 L3 / Lambda_L^5

(4) F_{L3}  = dW/dL3  = m L1 + eps L1 L2 / Lambda_L^5

SUSY breaking analysis:
  From (1): F_{X_L} = f + g (L1)^2.
  Setting all F-terms to zero simultaneously:

  From (3): eps L1 L3 / Lambda^5 = 0.
    If eps != 0: either L1 = 0 or L3 = 0.

  Case A: L1 = 0 (metastable vacuum at small eps):
    - (1) becomes: F_{X_L} = f != 0.  SUSY IS BROKEN.
    - (2) becomes: m L3 = 0 => L3 = 0.
    - (3) becomes: 0 = 0 (automatic).
    - (4) becomes: 0 = 0 (automatic).

    Metastable vacuum: <L1> = <L3> = 0, X_L and L2 are pseudo-moduli.
    F_{X_L} = f != 0 => SUSY broken by rank condition.

  Case B: L3 = 0 (from eq (3)):
    - (4) becomes: m L1 + eps L1 L2 / Lambda^5 = 0 => L1 (m + eps L2/Lambda^5) = 0
    - Sub-case B1: L1 = 0 => same as Case A.
    - Sub-case B2: L1 != 0, L2 = -m Lambda^5 / eps (far away for small eps)
      - (1): f + g L1^2 = 0 => L1^2 = -f/g  (requires f/g < 0)
      - (2): 0 + 2g X_L L1 + eps L2 * 0 / Lambda^5 = 0 => X_L = 0 (since L3=0)
      - This is the DISTANT SUSY-PRESERVING vacuum (Pfaffian-restored).
      - Parametrically far: L2 ~ m Lambda^5 / eps >> Lambda for eps << 1.

CONCLUSION: At the metastable vacuum (Case A),
  <L1> = <L3> = 0
  <X_L>, <L2> = pseudo-moduli (lifted at one loop)
  F_{X_L} = f != 0 => SUSY BROKEN
""")


print("=" * 72)
print("PART (b): FERMION MASS MATRIX W_{IJ}")
print("=" * 72)

print("""
Second derivatives W_{IJ} = d^2 W / dPhi_I dPhi_J:

Basis: (X_L, L1, L2, L3) where L1=L^{12}, L2=L^{13}, L3=L^{23}.

From W = f X_L + m L1 L3 + g X_L L1^2 + eps L1 L2 L3 / Lambda^5:

  W_{X_L, X_L} = 0
  W_{X_L, L1}  = 2g L1
  W_{X_L, L2}  = 0
  W_{X_L, L3}  = 0

  W_{L1, L1}   = 2g X_L + eps*0/Lambda^5  [d^2(eps L1 L2 L3)/dL1^2 = 0]
                  Wait: d(eps L1 L2 L3/Lambda^5)/dL1 = eps L2 L3/Lambda^5
                  d^2/dL1^2 = 0.
                  So W_{L1,L1} = 2g X_L.

  W_{L1, L2}   = eps L3 / Lambda^5
  W_{L1, L3}   = m + eps L2 / Lambda^5

  W_{L2, L2}   = 0
  W_{L2, L3}   = eps L1 / Lambda^5

  W_{L3, L3}   = 0

Full 4x4 matrix:

         X_L        L1                    L2                  L3
X_L [    0,      2g L1,                  0,                  0              ]
L1  [  2g L1,    2g X_L,          eps L3/Lam^5,     m + eps L2/Lam^5       ]
L2  [    0,    eps L3/Lam^5,             0,          eps L1/Lam^5          ]
L3  [    0,   m + eps L2/Lam^5,   eps L1/Lam^5,           0               ]
""")


print("=" * 72)
print("PART (c): MASS MATRIX AT METASTABLE VACUUM")
print("=" * 72)

print("""
At the metastable vacuum: <L1> = <L3> = 0, <X_L> = v (pseudo-modulus), <L2> arbitrary.

Setting eps -> 0 (metastable regime), the 4x4 reduces to:

         X_L    L1     L2    L3
X_L [    0,     0,     0,    0   ]
L1  [    0,    2gv,    0,    m   ]
L2  [    0,     0,     0,    0   ]
L3  [    0,     m,     0,    0   ]

X_L decouples (zero row/column) => goldstino direction.
L2 decouples (zero row/column) => another massless mode (pseudo-modulus partner).

The 2x2 massive block for (L1, L3):

  M_2x2 = [[2gv,  m],
            [m,    0]]

Characteristic polynomial: lambda^2 - 2gv * lambda - m^2 = 0

  lambda = gv +/- sqrt((gv)^2 + m^2)

Physical masses (absolute values):
  |lambda_+| = gv + sqrt((gv)^2 + m^2)
  |lambda_-| = |gv - sqrt((gv)^2 + m^2)| = sqrt((gv)^2 + m^2) - gv

Including the two massless modes, the full spectrum is:

  (0, 0, sqrt((gv)^2 + m^2) - gv, sqrt((gv)^2 + m^2) + gv)

The FERMION MASS SPECTRUM (L-sector only, 3x3 block excluding X_L):

  M_3x3 = [[2gv,  0,  m],
            [0,    0,  0],
            [m,    0,  0]]

Eigenvalues: 0, and lambda^2 - 2gv*lambda - m^2 = 0.

At t = gv/m = sqrt(3):

  lambda_+ = m(sqrt(3) + sqrt(3+1)) = m(sqrt(3) + 2) = m(2 + sqrt(3))
  lambda_- = m(sqrt(3) - 2) => |lambda_-| = m(2 - sqrt(3))

Spectrum: (0, (2-sqrt(3))m, (2+sqrt(3))m)
""")

# Numerical verification
m_unit = 1.0  # mass unit
gv = sqrt(3.0) * m_unit
M3 = array([
    [2*gv,   0,       m_unit],
    [0,      0,       0     ],
    [m_unit, 0,       0     ]
])

evals = np.linalg.eigvalsh(M3)
phys_masses = np.sort(np.abs(evals))

print(f"Numerical verification (m=1, gv/m = sqrt(3)):")
print(f"  M_3x3 eigenvalues: {evals}")
print(f"  Physical masses:   {phys_masses}")
print(f"  Analytic: (0, 2-sqrt(3), 2+sqrt(3)) = (0, {2-sqrt(3):.10f}, {2+sqrt(3):.10f})")
print(f"  Match lambda_0: {abs(phys_masses[0]) < 1e-14}")
print(f"  Match lambda_-: {abs(phys_masses[1] - (2-sqrt(3))) < 1e-14}")
print(f"  Match lambda_+: {abs(phys_masses[2] - (2+sqrt(3))) < 1e-14}")
print()

# Koide Q verification
m0, mm, mp = 0.0, (2-sqrt(3))*m_unit, (2+sqrt(3))*m_unit
num = m0 + mm + mp
denom = (sqrt(m0) + sqrt(mm) + sqrt(mp))**2
Q = num / denom
print(f"Koide Q = ({m0} + {mm:.10f} + {mp:.10f}) / (sqrt({m0}) + sqrt({mm:.10f}) + sqrt({mp:.10f}))^2")
print(f"        = {num:.10f} / {denom:.10f}")
print(f"        = {Q:.15f}")
print(f"  2/3   = {2/3:.15f}")
print(f"  |Q - 2/3| = {abs(Q - 2/3):.2e}")
print()

# Algebraic proof
print("Algebraic proof Q = 2/3:")
print(f"  Numerator = (2-sqrt(3)) + (2+sqrt(3)) = 4m")
print(f"  A = sqrt(2-sqrt(3)), B = sqrt(2+sqrt(3))")
print(f"  A*B = sqrt((2-sqrt(3))(2+sqrt(3))) = sqrt(4-3) = 1")
print(f"  A^2 + B^2 = (2-sqrt(3)) + (2+sqrt(3)) = 4")
print(f"  (A+B)^2 = A^2 + B^2 + 2AB = 4 + 2 = 6")
print(f"  Denominator = m*(A+B)^2 = 6m")
print(f"  Q = 4m / 6m = 2/3.  QED.")
print()


print("=" * 72)
print("PART (d): MEDIATION AND LEPTON MASS SCALE")
print("=" * 72)

print("""
The O'Raifeartaigh seed spectrum is (0, (2-sqrt(3))m, (2+sqrt(3))m).

In the Koide parametrization, masses are:
  m_k = M0 (1 + sqrt(2) cos(delta + 2*pi*k/3))^2

At the seed delta = 3*pi/4:
  m_1 = 0  (the cos = -1/sqrt(2) makes the factor vanish)
  m_2 = M0 (1 + sqrt(2) cos(3pi/4 + 2pi/3))^2
  m_3 = M0 (1 + sqrt(2) cos(3pi/4 + 4pi/3))^2

The seed total mass:
  M0_seed = m_1 + m_2 + m_3 = (2-sqrt(3))m_OR + (2+sqrt(3))m_OR = 4*m_OR

where m_OR is the O'Raifeartaigh mass parameter.

The physical lepton masses have:
  M0_phys = m_e + m_mu + m_tau = 1882.67 MeV

From the Koide parametrization:
  v0 = sqrt(M0/6)   [where M0 = sum of masses and v0^2 = M0/6 from
  sum_k m_k = 3M0 with M0 = (sum sigma_k / 3)^2, so sum m_k = 6*v0^2
  Wait, let me be more careful.]
""")

# Koide parametrization: sigma_k = v0 (1 + sqrt(2) cos(delta + 2pi*k/3))
# m_k = sigma_k^2 = v0^2 (1 + sqrt(2) cos(...))^2
# sum sigma_k = 3*v0 (since sum cos(... + 2pi*k/3) = 0)
# sum m_k = v0^2 * sum_k (1 + sqrt(2) cos(...))^2
#         = v0^2 * [3 + 2*sqrt(2)*sum cos(...) + 2*sum cos^2(...)]
#         = v0^2 * [3 + 0 + 2*(3/2)]  (since sum cos^2 = 3/2)
#         = v0^2 * [3 + 3] = 6*v0^2

v0_phys_sq = M0_lepton / 6.0
v0_phys = sqrt(v0_phys_sq)
print(f"Physical lepton Koide parameters:")
print(f"  M0 = m_e + m_mu + m_tau = {M0_lepton:.4f} MeV")
print(f"  sum m_k = 6 v0^2  =>  v0^2 = M0/6 = {v0_phys_sq:.6f} MeV")
print(f"  v0 = {v0_phys:.6f} MeV^(1/2)")
print()

# Physical delta for leptons
sigma_e   = sqrt(m_e)
sigma_mu  = sqrt(m_mu)
sigma_tau = sqrt(m_tau)
v0_check  = (sigma_e + sigma_mu + sigma_tau) / 3.0
print(f"  v0 from signed roots: (sqrt(m_e) + sqrt(m_mu) + sqrt(m_tau))/3 = {v0_check:.6f}")
print(f"  v0^2 = {v0_check**2:.6f} MeV")
print(f"  6*v0^2 = {6*v0_check**2:.4f} MeV  (should match M0 = {M0_lepton:.4f})")
print()

# Compute physical delta for leptons
phi = array([0, 2*pi/3, 4*pi/3])
sigmas = array([sigma_e, sigma_mu, sigma_tau])
A_cos = np.dot(sigmas, np.cos(phi))
B_sin = np.dot(sigmas, -np.sin(phi))
delta_phys = np.arctan2(B_sin, A_cos) % (2*pi)
print(f"  Physical delta = {delta_phys:.8f} rad = {np.degrees(delta_phys):.4f} deg")
print(f"  Seed delta = 3pi/4 = {3*pi/4:.8f} rad = 135.0000 deg")
print(f"  Bloom rotation = {delta_phys - 3*pi/4:.8f} rad = {np.degrees(delta_phys - 3*pi/4):.4f} deg")
print()

# Koide Q for leptons
Q_lepton = M0_lepton / (sigma_e + sigma_mu + sigma_tau)**2
print(f"  Koide Q(e,mu,tau) = {Q_lepton:.10f}")
print(f"  2/3 = {2/3:.10f}")
print(f"  deviation = {(Q_lepton - 2/3)/(2/3) * 100:.4f}%")
print()

# Seed scale
print("Seed scale from O'Raifeartaigh model:")
print(f"  Seed total mass = 4 * m_OR")
print(f"  Seed v0_seed^2 = (4*m_OR)/6 = 2*m_OR/3")
print()

# v0-doubling for leptons
# From bloom_mechanism.md: v0-doubling does NOT hold for leptons (ratio 1.01)
# But let's check what m_OR would be with and without v0-doubling.

# Without v0-doubling (seed = physical):
m_OR_direct = M0_lepton / 4.0
print(f"Case 1: No v0-doubling (seed IS the physical spectrum)")
print(f"  4*m_OR = M0 = {M0_lepton:.4f} MeV")
print(f"  m_OR = M0/4 = {m_OR_direct:.4f} MeV")
print(f"  v0_seed = v0_phys = {v0_phys:.6f} MeV^(1/2)")
print()

# With v0-doubling (v0_phys = 2 * v0_seed):
print(f"Case 2: With v0-doubling (v0_phys = 2 * v0_seed)")
v0_seed_doubled = v0_check / 2.0
M0_seed_doubled = 6 * v0_seed_doubled**2
m_OR_doubled = M0_seed_doubled / 4.0
print(f"  v0_seed = v0_phys/2 = {v0_seed_doubled:.6f} MeV^(1/2)")
print(f"  M0_seed = 6*v0_seed^2 = {M0_seed_doubled:.4f} MeV")
print(f"  m_OR = M0_seed/4 = {m_OR_doubled:.4f} MeV")
print()

# Actually, from the existing analysis, v0-doubling does NOT hold for leptons.
# The bloom_mechanism.md says ratio = 1.01 for leptons.
# So we should use Case 1.

# Let's compute the Koide-parametrization seed for leptons (sigma_1=0 at delta=3pi/4)
delta_seed = 3*pi/4
c1 = 1 + sqrt(2)*np.cos(delta_seed)
c2 = 1 + sqrt(2)*np.cos(2*pi/3 + delta_seed)
c3 = 1 + sqrt(2)*np.cos(4*pi/3 + delta_seed)
print(f"Seed coefficients at delta = 3pi/4:")
print(f"  c1 = {c1:.15f}  (should be ~0)")
print(f"  c2 = {c2:.10f}")
print(f"  c3 = {c3:.10f}")
print()

# Fit v0_seed to the two non-zero physical roots
# sigma_2 = v0_seed * c2, sigma_3 = v0_seed * c3
# Least-squares: v0_seed = (sigma_mu * c2 + sigma_tau * c3) / (c2^2 + c3^2)
v0_seed_fit = (sigma_mu * c2 + sigma_tau * c3) / (c2**2 + c3**2)
M0_seed_fit = v0_seed_fit**2
m2_seed = (v0_seed_fit * c2)**2
m3_seed = (v0_seed_fit * c3)**2
m_OR_from_seed = (m2_seed + m3_seed) / 4.0  # since seed total = 4*m_OR

print(f"Koide seed fitted to (e,mu,tau) physical roots:")
print(f"  v0_seed = {v0_seed_fit:.6f} MeV^(1/2)")
print(f"  M0_seed = {M0_seed_fit:.6f} MeV")
print(f"  Seed masses: (0, {m2_seed:.4f}, {m3_seed:.4f}) MeV")
print(f"  Physical:    ({m_e:.4f}, {m_mu:.4f}, {m_tau:.4f}) MeV")
print(f"  Seed total = {m2_seed + m3_seed:.4f} MeV")
print(f"  m_OR (from seed total/4) = {m_OR_from_seed:.4f} MeV")
print(f"  v0_phys / v0_seed = {v0_check / v0_seed_fit:.6f}")
print()

# Now map to OR model: seed spectrum = (0, (2-sqrt(3))*m_OR, (2+sqrt(3))*m_OR)
# Check that the seed mass ratio matches
ratio_seed_masses = m3_seed / m2_seed if m2_seed > 0 else float('inf')
ratio_OR = (2+sqrt(3)) / (2-sqrt(3))
print(f"Seed mass ratio m3/m2 = {ratio_seed_masses:.6f}")
print(f"O'Raifeartaigh ratio (2+sqrt(3))/(2-sqrt(3)) = {ratio_OR:.6f} = (2+sqrt(3))^2 = {(2+sqrt(3))**2:.6f}")
print(f"Match: {abs(ratio_seed_masses - ratio_OR) / ratio_OR * 100:.4f}% deviation")
print()

# The seed at delta=3*pi/4 should automatically give the OR ratio
# because the Koide parametrization at delta=3*pi/4 with Q=2/3 is
# equivalent to the OR spectrum.
# Let's verify: at delta=3pi/4, the masses are:
#   m_k = v0^2 * c_k^2
# So m_2/m_3 = (c_2/c_3)^2

ratio_c = (c2/c3)**2
print(f"c2^2/c3^2 = {ratio_c:.10f}")
print(f"(2-sqrt(3))/(2+sqrt(3)) = {(2-sqrt(3))/(2+sqrt(3)):.10f}")
# These should be the same if the parametrization is consistent
print(f"Match: {abs(ratio_c - (2-sqrt(3))/(2+sqrt(3))) < 1e-10}")
print()

# Now compute m_OR for physical leptons
# Method: identify the total seed mass with 4*m_OR
# For leptons (approximately Koide-exact), the seed is close to the physical:
# v0 ratio ~ 1.01 means the bloom is tiny.
# Direct: use the seed masses to get m_OR.

m_OR_lepton = (m2_seed + m3_seed) / 4.0
print(f"O'Raifeartaigh mass parameter for lepton sector:")
print(f"  m_OR = (seed total)/4 = {m_OR_lepton:.4f} MeV")
print(f"  Physical masses from OR spectrum after bloom:")
print(f"    m_OR = {m_OR_lepton:.4f} MeV")
print(f"    (2-sqrt(3))*m_OR = {(2-sqrt(3))*m_OR_lepton:.4f} MeV")
print(f"    (2+sqrt(3))*m_OR = {(2+sqrt(3))*m_OR_lepton:.4f} MeV")
print()


print("=" * 72)
print("PART (e): CONFINEMENT SCALE Lambda_L")
print("=" * 72)

print("""
For SU(2) SQCD with N_f = 3 (s-confining), the confinement scale Lambda_L
sets the overall dimensionful scale.

The dynamical superpotential is W_dyn = Pf(L) / Lambda_L^5 = L^{12} L^{13} L^{23} / Lambda_L^5.

The meson fields L^{ij} have mass dimension 2 (composites of two fundamental fields).
The O'Raifeartaigh deformation parameters have dimensions:
  [f] = mass^3    (from f X_L, where [X_L] = mass)
  [m] = mass^0 = dimensionless?  No.

Wait, let me be more careful about dimensions.

In 4d N=1 SUSY, the superpotential has mass dimension 3.
  [W] = 3
  [X_L] = 1  (elementary field)
  [L^{ij}] = 2  (composite of two fundamentals, each dim 1 at free-field level)

Actually, in the s-confining description, the L^{ij} are the confined
degrees of freedom and should be treated as having dimension 1 (they are
the low-energy effective fields). The relation L ~ Psi Psi introduces a
factor of Lambda_L to make dimensions work:
  L_eff = Psi Psi / Lambda_L  => [L_eff] = 1.

Let me use the canonical normalization where all fields have dimension 1.
Then:
  [f X_L] = 3  =>  [f] = 2
  [m L1 L3] = 3  =>  [m] = 1
  [g X_L L1^2] = 3  =>  [g] = 0 (dimensionless)
  [eps L1 L2 L3 / Lambda^{2N_c - 1}] = [eps * L^3 / Lambda^{2*2-1}] = 3
  For SU(2): 2*N_c - 1 = 3, so Lambda^3.
  [eps L1 L2 L3 / Lambda^3] = [eps] + 3 - 3 = 3 => [eps] = 3? No.
  [eps] * 3 / 3 should give 3. So [eps] = 3... that's wrong.

Let me reconsider. For SU(N_c) with N_f = N_c + 1 (s-confining):
  The dynamical scale appears as Lambda^{3N_c - N_f} = Lambda^{3*2 - 3} = Lambda^3.
  Wait, the instanton factor for SU(N_c) SQCD is Lambda^{b_0} where
  b_0 = 3N_c - N_f = 6 - 3 = 3.

  The s-confining superpotential for SU(2) with N_f = 3 is:
  W_dyn = Pf(M) / Lambda^{b_0} where M^{ij} ~ Psi^i Psi^j (dimension 2).
  [Pf(M)] for 3x3 antisymmetric (actually the Pfaffian of 2N_f x 2N_f...
  no, for SU(2), the mesons are L^{ij} = eps^{ab} Psi^i_a Psi^j_b).

  With [L^{ij}] = 2 and Pf(L) = L^{12} L^{13} L^{23}, [Pf] = 6.
  W_dyn = Pf(L) / Lambda^{b_0} has [W_dyn] = 6 - b_0 = 3 => b_0 = 3. Checks out.

  Now with canonically normalized L_eff^{ij} = L^{ij} / Lambda_L:
  [L_eff] = 1.
  W_dyn = Lambda_L^3 L_eff^{12} L_eff^{13} L_eff^{23} / Lambda_L^3
        = L_eff^{12} L_eff^{13} L_eff^{23}
  That's dimensionless, not dim 3. Wrong.

  Correct: L^{ij} has dim 2, there are 3 of them in the Pfaffian (dim 6),
  divided by Lambda^3 gives dim 3. But L_eff = L/Lambda has dim 1.
  Pf(L_eff) = L_eff^{12} L_eff^{13} L_eff^{23} has dim 3.
  And L = Lambda * L_eff, so Pf(L) = Lambda^3 Pf(L_eff).
  W_dyn = Lambda^3 Pf(L_eff) / Lambda^3 = Pf(L_eff). Dim 3. OK.

  So in terms of canonical fields:
    W_dyn = L_eff^{12} L_eff^{13} L_eff^{23}   (coefficient 1, no Lambda)

  And the full superpotential is:
    W = f X_L + m L1 L3 + g X_L L1^2 + eps L1 L2 L3

  where all fields have dim 1 and:
  [f] = 2, [m] = 1, [g] = 0, [eps] = 0.

  The "eps" is just a dimensionless coupling for the Pfaffian term in canonical fields.
""")

print(f"With canonical normalization:")
print(f"  L_eff = L / Lambda_L  ([L_eff] = 1)")
print(f"  W = f X_L + m L1_eff L3_eff + g X_L (L1_eff)^2 + eps L1_eff L2_eff L3_eff")
print(f"  [f] = mass^2, [m] = mass, [g] = dimensionless, [eps] = dimensionless")
print()

print(f"The O'Raifeartaigh mass parameter m sets the fermion mass scale:")
print(f"  Spectrum: (0, (2-sqrt(3))*m, (2+sqrt(3))*m)")
print(f"  Total mass = 4*m")
print()

# For the lepton sector, the seed total mass should be ~ M0_lepton
# (since v0 doubling ratio is ~1.01, not 2)
print(f"Lepton sector scale setting:")
print(f"  Seed total = 4*m_OR = {4*m_OR_lepton:.4f} MeV")
print(f"  Physical total M0 = {M0_lepton:.4f} MeV")
print(f"  Ratio M0_phys/M0_seed = {M0_lepton/(4*m_OR_lepton):.6f}  (near 1 since bloom is small)")
print()

# The confinement scale Lambda_L
# In the s-confining theory, the physical meson mass ~ m_OR and the confinement
# scale Lambda_L are related. The O'R parameters are:
#   m ~ h * mu  (ISS mapping: Yukawa * VEV)
#   f ~ h * mu^2  (SUSY-breaking scale)
# where mu ~ Lambda_L (the meson VEV is of order Lambda_L in the confined description).
#
# More precisely, the ISS-like mapping gives:
#   m = h * Lambda_L  (h = magnetic Yukawa, Lambda_L = confinement scale)
#   f = h * Lambda_L^2
#
# If h ~ O(1), then m ~ Lambda_L.

print(f"Confinement scale estimate:")
print(f"  If m_OR ~ Lambda_L (with h ~ O(1)):")
m_OR_val = m_OR_lepton  # the fitted value
print(f"    Lambda_L ~ m_OR = {m_OR_val:.2f} MeV")
print()
print(f"  For comparison, from MEMORY.md: 'Λ_L ~ 300 MeV gives M₀ = 313.8 MeV'")
print(f"  But our fitted m_OR = {m_OR_val:.2f} MeV gives M0_seed = {4*m_OR_val:.2f} MeV")
print(f"  And our v0^2 = M0/6 = {M0_lepton/6:.2f} MeV gives M0 = {M0_lepton:.2f} MeV")
print()

# Let's check: if Lambda_L = 300 MeV and m_OR = Lambda_L, what M0 would we get?
Lambda_test = 300.0  # MeV
m_OR_test = Lambda_test
M0_test = 4 * m_OR_test
print(f"  If Lambda_L = 300 MeV, m_OR = Lambda_L = 300 MeV:")
print(f"    M0_seed = 4*m_OR = {M0_test:.0f} MeV")
print(f"    But M0_phys = {M0_lepton:.2f} MeV, so need bloom factor = {M0_lepton/M0_test:.4f}")
print()

# More likely: m_OR is not equal to Lambda_L. The relation depends on the Yukawa h.
# m_OR = h * Lambda_L. With the Koide condition gv/m = sqrt(3) and g = h,
# the condition is <X>/Lambda_L = sqrt(3) (same as ISS).
#
# If Lambda_L sets the scale and h is the free coupling:
# m_OR = h * Lambda_L, so 4*h*Lambda_L = M0_seed.

# Let's instead just directly state what m_OR must be.
print(f"Direct determination of m_OR from lepton masses:")
print(f"  The O'R seed at gv/m = sqrt(3) gives spectrum (0, (2-sqrt(3))*m, (2+sqrt(3))*m)")
print(f"  with total mass 4m.")
print()
print(f"  The physical lepton Koide seed has total mass approximately equal to M0_phys")
print(f"  (since the bloom for leptons is very small: v0 ratio = 1.01).")
print()
print(f"  Therefore: m_OR = M0_phys / 4 = {M0_lepton/4:.4f} MeV")
print(f"  Seed spectrum: (0, {(2-sqrt(3))*M0_lepton/4:.4f}, {(2+sqrt(3))*M0_lepton/4:.4f}) MeV")
print()

# Compare with actual lepton masses at the seed
m_OR_final = M0_lepton / 4.0
seed_minus = (2-sqrt(3)) * m_OR_final
seed_plus  = (2+sqrt(3)) * m_OR_final
print(f"  O'R seed:  (0, {seed_minus:.4f}, {seed_plus:.4f}) MeV")
print(f"  After bloom to delta_phys = {np.degrees(delta_phys):.4f} deg:")
# Reconstruct masses from Koide parametrization at physical delta
v0_val = v0_check
for k in range(3):
    c_k = 1 + sqrt(2)*np.cos(delta_phys + 2*pi*k/3)
    m_k = v0_val**2 * c_k**2
    print(f"    m_{k+1} = {m_k:.4f} MeV  (physical: {[m_e, m_mu, m_tau][k]:.4f} MeV)")

print()


print("=" * 72)
print("PART (f): THE GOLDSTINO")
print("=" * 72)

print("""
At the SUSY-breaking metastable vacuum:
  F_{X_L} = f != 0
  All other F-terms vanish (at leading order in eps).

The goldstino is the fermionic component of the superfield whose F-term
breaks SUSY. In a general O'Raifeartaigh model, the goldstino direction
in field space is:

  psi_goldstino = sum_I F_I psi_I / sqrt(sum_I |F_I|^2)

where psi_I is the fermionic component of superfield Phi_I and F_I = <dW/dPhi_I>.

At the metastable vacuum:
  F_{X_L} = f    (nonzero)
  F_{L1}  = 0
  F_{L2}  = 0
  F_{L3}  = 0

Therefore:
  psi_goldstino = psi_{X_L}

The goldstino is PURELY the fermionic component of X_L.

This is consistent with the mass matrix structure: the X_L row and column
of W_{IJ} are entirely zero at the vacuum (since <L1> = 0), giving
psi_{X_L} a zero mass eigenvalue. It is the Nambu-Goldstone fermion
of broken SUSY.

If SUGRA is included:
  The gravitino absorbs the goldstino (super-Higgs mechanism).
  Gravitino mass: m_{3/2} = f / (sqrt(3) M_Pl)
  where f is the SUSY-breaking F-term and M_Pl is the reduced Planck mass.

  With f = m_OR * Lambda_L ~ (470 MeV)^2 (if Lambda_L ~ m_OR ~ 470 MeV):
  m_{3/2} ~ (470 MeV)^2 / (sqrt(3) * 2.4 x 10^{21} MeV)
          ~ 5.3 x 10^{-17} MeV
          ~ 5.3 x 10^{-11} eV

  This is an extraordinarily small gravitino mass — far below any direct
  detection threshold, but potentially relevant for cosmology.
""")

# Gravitino mass estimate
M_Pl_reduced = 2.435e21  # MeV (reduced Planck mass)
Lambda_L_est = m_OR_final  # using m_OR as estimate for Lambda_L
f_est = m_OR_final * Lambda_L_est  # f ~ m * Lambda in natural units
m_32 = f_est / (sqrt(3) * M_Pl_reduced)
print(f"Gravitino mass estimate:")
print(f"  f ~ m_OR * Lambda_L ~ ({m_OR_final:.1f} MeV)^2 = {f_est:.1f} MeV^2")
print(f"  m_3/2 = f / (sqrt(3) M_Pl) = {m_32:.3e} MeV = {m_32*1e6:.3e} eV")
print()

# But the actual F-term dimension is [f] = mass^2 in the superpotential W = f*X_L.
# With [W] = mass^3, [X_L] = mass, [f] = mass^2.
# The SUSY-breaking order parameter is F_{X_L} = f (dimension mass^2).
# Gravitino mass: m_{3/2} = <F> / (sqrt(3) M_Pl)

# If f = h * mu^2 with h ~ O(1) and mu ~ Lambda_L:
f_susy_breaking = Lambda_L_est**2  # f ~ Lambda_L^2
m_32_v2 = f_susy_breaking / (sqrt(3) * M_Pl_reduced)
print(f"Alternative (f = Lambda_L^2):")
print(f"  f ~ Lambda_L^2 = ({Lambda_L_est:.1f} MeV)^2 = {f_susy_breaking:.1f} MeV^2")
print(f"  m_3/2 = {m_32_v2:.3e} MeV = {m_32_v2*1e6:.3e} eV")
print()


print("=" * 72)
print("COMPLETE SUMMARY")
print("=" * 72)

print(f"""
MODEL: SU(2) SQCD with N_f = 3 (s-confining)
FIELDS: X_L, L^{{12}}, L^{{13}}, L^{{23}} (Lagrange multiplier + 3 mesons)

SUPERPOTENTIAL:
  W = f X_L + m L^{{12}} L^{{23}} + g X_L (L^{{12}})^2 + eps Pf(L)/Lambda_L^3

METASTABLE VACUUM:
  <L^{{12}}> = <L^{{23}}> = 0
  <X_L>, <L^{{13}}> = pseudo-moduli (lifted at one loop)
  F_{{X_L}} = f != 0 => SUSY BROKEN

FERMION MASS SPECTRUM (at gv/m = sqrt(3)):
  Goldstino:  m_0 = 0   (psi_{{X_L}})
  Light:      m_- = (2 - sqrt(3)) * m_OR = {(2-sqrt(3))*m_OR_final:.4f} MeV
  Heavy:      m_+ = (2 + sqrt(3)) * m_OR = {(2+sqrt(3))*m_OR_final:.4f} MeV

KOIDE VERIFICATION:
  Q(0, m_-, m_+) = 4m / 6m = 2/3 EXACTLY
  This is the Koide seed with delta = 3*pi/4.

SCALE:
  m_OR = M0_lepton / 4 = {m_OR_final:.4f} MeV
  Lambda_L ~ m_OR ~ {m_OR_final:.0f} MeV  (if magnetic Yukawa h ~ O(1))

GOLDSTINO:
  psi_goldstino = psi_{{X_L}}  (purely the X_L fermion)

BLOOM:
  Seed delta = 135.00 deg
  Physical delta = {np.degrees(delta_phys):.4f} deg
  Bloom rotation = {np.degrees(delta_phys - 3*pi/4):.4f} deg  (very small for leptons)
  v0_phys / v0_seed = {v0_check/v0_seed_fit:.4f}  (near 1; no v0-doubling for leptons)
""")
