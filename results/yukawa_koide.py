"""
Yukawa coupling matrices from the Koide parametrization with tan beta = 1.

The Koide parametrization writes three signed square roots as:

    sqrt(m_k) = sqrt(M0) * (1 + sqrt(2) * cos(2*pi*k/3 + delta)),  k = 0,1,2

In a two-Higgs-doublet model with tan beta = v_u/v_d = 1, both vevs equal
v/sqrt(2) and the Yukawa couplings are y_k = 2*m_k/v for all fermion types.

The Koide parametrization transfers directly to Yukawa couplings:
    sqrt(y_k) = sqrt(y0) * (1 + sqrt(2) * cos(2*pi*k/3 + delta))
with y0 = 2*M0/v.

This script computes:
1. Koide parameters (M0, delta) and Yukawa couplings for each triple
2. Cartan subalgebra decomposition of sqrt(Y) matrices
3. Seed Yukawa matrices at delta = 3*pi/4 and bloom rotations
4. Analytic bloom rotation in the Cartan plane
5. Numerical verification of the Koide condition on Yukawa couplings
"""

import numpy as np

# ---------------------------------------------------------------------------
# Physical constants
# ---------------------------------------------------------------------------
v_higgs = 246.22  # GeV, electroweak vev

# Masses in MeV
m_e   = 0.51100
m_mu  = 105.658
m_tau = 1776.86
m_s   = 93.4
m_c   = 1270.0
m_b   = 4180.0
m_t   = 172760.0

# Convert to GeV for Yukawa computation
GeV = 1.0e3  # 1 GeV = 1000 MeV

SQRT2         = np.sqrt(2.0)
SQRT3         = np.sqrt(3.0)
SQRT6         = np.sqrt(6.0)
TWO_PI_OVER_3 = 2.0 * np.pi / 3.0
DELTA_SEED    = 3.0 * np.pi / 4.0  # 135 deg, forces sigma_0 = 0


# ---------------------------------------------------------------------------
# Gell-Mann diagonal matrices (3x3)
# ---------------------------------------------------------------------------
lambda3 = np.diag([1.0, -1.0, 0.0])
lambda8 = np.diag([1.0, 1.0, -2.0]) / SQRT3


# ---------------------------------------------------------------------------
# Core functions
# ---------------------------------------------------------------------------

def compute_M0_delta(signed_roots):
    """
    Given signed roots sigma = [s0, s1, s2], return (M0, delta).
    sqrt(M0) = sum(sigma) / 3
    delta = arctan2(-sum(sigma*sin(phi)), sum(sigma*cos(phi)))
    """
    s = np.array(signed_roots, dtype=float)
    phi = np.array([TWO_PI_OVER_3 * k for k in range(3)])
    sqrtM0 = s.sum() / 3.0
    M0 = sqrtM0 ** 2
    A = np.dot(s, np.cos(phi))
    B = np.dot(s, -np.sin(phi))
    delta = np.arctan2(B, A)
    return M0, delta


def koide_Q_signed(signed_roots):
    """Q_signed = sum(m) / (sum(sigma))^2 with signed roots."""
    s = np.array(signed_roots, dtype=float)
    return s.dot(s) / s.sum() ** 2


def masses_from_M0_delta(M0, delta):
    """Compute three masses from Koide parameters (M0, delta)."""
    sqrtM0 = np.sqrt(M0)
    sv = np.array([sqrtM0 * (1.0 + SQRT2 * np.cos(TWO_PI_OVER_3 * k + delta))
                    for k in range(3)])
    return sv ** 2, sv


def yukawa_from_mass(m_MeV):
    """y = 2*m/v with m in MeV, v in GeV. Returns dimensionless y."""
    return 2.0 * (m_MeV / GeV) / v_higgs


def cartan_decomposition(sqrt_y):
    """
    Decompose sqrt(Y) = diag(sqrt(y_0), sqrt(y_1), sqrt(y_2)) as:
        sqrt(Y) = sqrt(y0) * (I + (sqrt(6)/2) * (n1*lambda3 + n2*lambda8))
    where sqrt(y0) = Tr(sqrt(Y))/3 = (sum of sqrt(y_k))/3.

    Returns sqrt_y0, n_hat = (n1, n2), theta (angle in Cartan plane).
    """
    sy = np.array(sqrt_y, dtype=float)
    sqrt_y0 = sy.sum() / 3.0

    # Traceless part: T = sqrt(Y) - sqrt(y0)*I
    T = np.diag(sy) - sqrt_y0 * np.eye(3)

    # Project onto lambda3 and lambda8:
    # Tr(lambda3 * T) = T[0,0] - T[1,1]
    # Tr(lambda8 * T) = (T[0,0] + T[1,1] - 2*T[2,2]) / sqrt(3)
    # From sqrt(Y) = sqrt(y0)*(I + sqrt(6)/2 * (n1*lambda3 + n2*lambda8)):
    #   T = sqrt(y0) * sqrt(6)/2 * (n1*lambda3 + n2*lambda8)
    # So: Tr(lambda3 * T) = sqrt(y0)*sqrt(6)/2 * n1 * Tr(lambda3^2)
    #                      = sqrt(y0)*sqrt(6)/2 * n1 * 2
    # Similarly for lambda8: Tr(lambda8 * T) = sqrt(y0)*sqrt(6)/2 * n2 * Tr(lambda8^2)
    #                                        = sqrt(y0)*sqrt(6)/2 * n2 * 2

    # Tr(lambda3^2) = 2, Tr(lambda8^2) = 2
    tr_l3_T = T[0, 0] - T[1, 1]
    tr_l8_T = (T[0, 0] + T[1, 1] - 2.0 * T[2, 2]) / SQRT3

    coeff = sqrt_y0 * SQRT6  # sqrt(y0) * sqrt(6)/2 * 2 = sqrt(y0)*sqrt(6)
    if abs(coeff) < 1e-15:
        return sqrt_y0, np.array([0.0, 0.0]), 0.0

    n1 = tr_l3_T / coeff
    n2 = tr_l8_T / coeff

    # n_hat should be a unit vector
    n_norm = np.sqrt(n1 ** 2 + n2 ** 2)
    n_hat = np.array([n1 / n_norm, n2 / n_norm]) if n_norm > 1e-15 else np.array([1.0, 0.0])

    theta = np.arctan2(n_hat[1], n_hat[0])

    return sqrt_y0, n_hat, theta, n_norm


def seed_masses_at_delta(delta_s, M0_s):
    """Compute three masses at given (M0, delta), returning masses and signed roots."""
    sqM0 = np.sqrt(M0_s)
    sv = np.array([sqM0 * (1.0 + SQRT2 * np.cos(TWO_PI_OVER_3 * k + delta_s))
                    for k in range(3)])
    return sv ** 2, sv


# ---------------------------------------------------------------------------
# Define the three triples
# ---------------------------------------------------------------------------

triples = [
    {
        'name': '(e, mu, tau)',
        'signed_roots': [np.sqrt(m_e), np.sqrt(m_mu), np.sqrt(m_tau)],
        'masses_MeV': [m_e, m_mu, m_tau],
        'labels': ['e', 'mu', 'tau'],
    },
    {
        'name': '(-s, c, b)',
        'signed_roots': [-np.sqrt(m_s), np.sqrt(m_c), np.sqrt(m_b)],
        'masses_MeV': [m_s, m_c, m_b],
        'labels': ['s', 'c', 'b'],
    },
    {
        'name': '(c, b, t)',
        'signed_roots': [np.sqrt(m_c), np.sqrt(m_b), np.sqrt(m_t)],
        'masses_MeV': [m_c, m_b, m_t],
        'labels': ['c', 'b', 't'],
    },
]

# ═══════════════════════════════════════════════════════════════════════════
# Part 1: Koide parameters and Yukawa couplings
# ═══════════════════════════════════════════════════════════════════════════

print("=" * 74)
print("PART 1: Koide parameters and Yukawa couplings (tan beta = 1)")
print("=" * 74)
print(f"\nElectroweak vev: v = {v_higgs} GeV")
print(f"tan beta = 1: v_u = v_d = v/sqrt(2) = {v_higgs/SQRT2:.4f} GeV")
print(f"Yukawa coupling: y_k = 2 m_k / v  (for all fermion types at tan beta = 1)")
print()

results = []

for tr in triples:
    name = tr['name']
    sr = tr['signed_roots']
    masses = tr['masses_MeV']
    labels = tr['labels']

    M0, delta = compute_M0_delta(sr)
    Q = koide_Q_signed(sr)

    # Yukawa couplings
    y = [yukawa_from_mass(m) for m in masses]
    sqrt_y = [np.sqrt(yi) for yi in y]

    # Signed sqrt(y): same sign pattern as signed roots of masses
    signs = [np.sign(s) for s in sr]
    signed_sqrt_y = [sgn * sy for sgn, sy in zip(signs, sqrt_y)]

    # Koide Yukawa parametrization: y0 = 2*M0/v
    y0 = 2.0 * (M0 / GeV) / v_higgs
    sqrt_y0 = np.sqrt(y0)

    # Verify: sqrt_y0 should equal sum(signed_sqrt_y)/3
    sqrt_y0_check = sum(signed_sqrt_y) / 3.0

    # Verify y0 = (sum(signed_sqrt_y))^2 / 9
    y0_check = sum(signed_sqrt_y) ** 2 / 9.0

    # Q on Yukawa couplings (signed)
    Q_y = sum(y) / sum(signed_sqrt_y) ** 2

    print("-" * 74)
    print(f"Triple: {name}")
    print(f"  Masses (MeV):  {masses}")
    print(f"  Signed roots:  [{', '.join(f'{s:.6f}' for s in sr)}]")
    print()
    print(f"  Koide parameters:")
    print(f"    M0     = {M0:.6f} MeV")
    print(f"    sqrt(M0) = {np.sqrt(M0):.6f} MeV^(1/2)")
    print(f"    delta  = {delta:.8f} rad = {np.degrees(delta):.4f} deg")
    print(f"    Q      = {Q:.10f}  (2/3 = {2/3:.10f})")
    print()
    print(f"  Yukawa couplings (y_k = 2*m_k/v):")
    for i, (lbl, m, yi) in enumerate(zip(labels, masses, y)):
        print(f"    y_{lbl} = 2 * {m/GeV:.6e} GeV / {v_higgs} GeV = {yi:.6e}")
    print()
    print(f"  Koide Yukawa parametrization:")
    print(f"    y0 = 2*M0/v = 2 * {M0/GeV:.6e} / {v_higgs} = {y0:.6e}")
    print(f"    sqrt(y0) = {sqrt_y0:.6e}")
    print(f"    Check: (sum signed sqrt(y))^2 / 9 = {y0_check:.6e}  (ratio to y0: {y0_check/y0:.10f})")
    print(f"    Check: sum(signed sqrt(y)) / 3    = {sqrt_y0_check:.6e}  (ratio to sqrt(y0): {sqrt_y0_check/sqrt_y0:.10f})")
    print()
    print(f"  Reconstructed sqrt(y_k) from Koide parametrization:")
    for k in range(3):
        c_k = 1.0 + SQRT2 * np.cos(TWO_PI_OVER_3 * k + delta)
        reconstructed = sqrt_y0 * c_k
        print(f"    k={k}: sqrt(y0)*(1+sqrt(2)*cos(2*pi*{k}/3+delta)) = {reconstructed:.6e}  "
              f"(signed sqrt(y_{labels[k]}): {signed_sqrt_y[k]:+.6e}, ratio: {reconstructed/signed_sqrt_y[k]:.10f})")
    print()
    print(f"  Q on Yukawa couplings (signed): {Q_y:.10f}  (2/3 = {2/3:.10f})")
    print()

    results.append({
        'name': name,
        'labels': labels,
        'masses': masses,
        'signed_roots': sr,
        'M0': M0,
        'delta': delta,
        'Q': Q,
        'y': y,
        'sqrt_y': sqrt_y,
        'signed_sqrt_y': signed_sqrt_y,
        'y0': y0,
        'sqrt_y0': sqrt_y0,
        'Q_y': Q_y,
        'signs': signs,
    })


# ═══════════════════════════════════════════════════════════════════════════
# Part 2: Cartan subalgebra decomposition
# ═══════════════════════════════════════════════════════════════════════════

print("\n" + "=" * 74)
print("PART 2: Cartan subalgebra decomposition of sqrt(Y)")
print("=" * 74)
print()
print("Decomposition: sqrt(Y) = sqrt(y0) * (I + (sqrt(6)/2) * n_norm * (n1*lambda3 + n2*lambda8))")
print("where n_hat = (n1, n2) = (cos theta, sin theta) is a unit vector in the Cartan plane.")
print()

for r in results:
    name = r['name']
    signed_sqrt_y = r['signed_sqrt_y']
    sqrt_y0 = r['sqrt_y0']
    delta = r['delta']

    # Use absolute sqrt(y) for the diagonal matrix
    # (the Cartan decomposition acts on the physical sqrt(Y) = diag(sqrt(y_k)))
    # But for the signed version (with negative entries), we use signed_sqrt_y
    sqrt_y0_c, n_hat, theta, n_norm = cartan_decomposition(signed_sqrt_y)

    print("-" * 74)
    print(f"Triple: {name}")
    print()
    print(f"  sqrt(Y) = diag({', '.join(f'{s:+.6e}' for s in signed_sqrt_y)})")
    print(f"  Tr(sqrt(Y))/3 = sqrt(y0) = {sqrt_y0_c:.6e}")
    print()
    print(f"  Traceless part T = sqrt(Y) - sqrt(y0)*I:")
    T_diag = [s - sqrt_y0_c for s in signed_sqrt_y]
    print(f"    T = diag({', '.join(f'{t:+.6e}' for t in T_diag)})")
    print()
    print(f"  Cartan coefficients:")
    print(f"    n_norm (radius in Cartan plane) = {n_norm:.8f}")
    print(f"    n_hat = (cos theta, sin theta) = ({n_hat[0]:+.8f}, {n_hat[1]:+.8f})")
    print(f"    theta = {theta:.8f} rad = {np.degrees(theta):.4f} deg")
    print()

    # Verify: n_norm should equal 1 for exact Koide (since r = sqrt(2)*v0)
    print(f"  Verification:")
    print(f"    n_norm should be 1 for exact Koide: n_norm = {n_norm:.8f}")
    print()

    # Relation between theta and delta
    # The Koide parametrization gives:
    #   signed_sqrt_y_k = sqrt(y0) * (1 + sqrt(2)*cos(2*pi*k/3 + delta))
    # The Cartan plane coordinates are:
    #   n1 = (s0 - s1) / (sqrt(2) * sqrt(y0) * sqrt(3))  via lambda3
    #   n2 = (s0 + s1 - 2*s2) / (sqrt(6) * sqrt(y0))      via lambda8
    # We can compute these analytically from delta.

    # Analytic theta from delta:
    # c_k = 1 + sqrt(2)*cos(2*pi*k/3 + delta)
    # n1 propto (c_0 - c_1) = sqrt(2)*(cos(delta) - cos(2*pi/3 + delta))
    #   = sqrt(2) * 2 * sin(pi/3) * sin(pi/3 + delta)
    #   = 2*sqrt(2) * (sqrt(3)/2) * sin(pi/3 + delta)
    #   = sqrt(6) * sin(pi/3 + delta)
    # n2 propto (c_0 + c_1 - 2*c_2)/sqrt(3)
    #   cos(delta) + cos(2pi/3+delta) - 2*cos(4pi/3+delta)
    #   Using sum-to-product: ... = 3*cos(delta + pi/3 + pi/3 - pi) ...
    # Actually, let's just compute numerically.
    phi_arr = np.array([0, TWO_PI_OVER_3, 2*TWO_PI_OVER_3])
    c_arr = SQRT2 * np.cos(phi_arr + delta)
    # lambda3 component: sum c_k * lambda3_kk = c_0 - c_1
    # lambda8 component: sum c_k * lambda8_kk = (c_0 + c_1 - 2*c_2)/sqrt(3)
    n1_analytic = c_arr[0] - c_arr[1]
    n2_analytic = (c_arr[0] + c_arr[1] - 2*c_arr[2]) / SQRT3
    theta_analytic = np.arctan2(n2_analytic, n1_analytic)

    print(f"    theta (from Cartan decomp) = {np.degrees(theta):.6f} deg")
    print(f"    theta (analytic from delta) = {np.degrees(theta_analytic):.6f} deg")
    print(f"    Difference: {np.degrees(theta - theta_analytic):.2e} deg")
    print()
    print(f"    Koide delta = {np.degrees(delta):.4f} deg")
    print(f"    theta - delta = {np.degrees(theta - delta):.4f} deg")

    # Analytic relation: theta = delta + constant?
    # From the components:
    #   n1 = sqrt(2)*(cos(delta) - cos(2pi/3+delta))
    #      = sqrt(2)*2*sin(pi/3)*sin(pi/3+delta) = sqrt(6)*sin(pi/3+delta)
    #   n2 = sqrt(2)*(cos(delta) + cos(2pi/3+delta) - 2*cos(4pi/3+delta))/sqrt(3)
    # Let's compute cos(delta) + cos(2pi/3+delta) - 2*cos(4pi/3+delta):
    #   = cos d + cos d cos(2p/3) - sin d sin(2p/3) - 2*(cos d cos(4p/3) - sin d sin(4p/3))
    #   = cos d + cos d*(-1/2) - sin d*(sqrt(3)/2) - 2*(cos d*(-1/2) - sin d*(-sqrt(3)/2))
    #   = cos d - cos d/2 - sqrt(3)*sin d/2 + cos d + sqrt(3)*sin d
    #   = 2*cos d - cos d/2 + sqrt(3)*sin d/2
    #   Hmm, let me just use: cos a + cos b = 2cos((a+b)/2)cos((a-b)/2)
    # Actually the cleanest way:
    #   cos(delta) + cos(2pi/3+delta) = 2*cos(pi/3+delta)*cos(pi/3) = cos(pi/3+delta)
    #   So n2_num = cos(pi/3+delta) - 2*cos(4pi/3+delta)
    #   cos(4pi/3+delta) = cos(delta+pi/3+pi) = -cos(delta+pi/3)
    #   n2_num = cos(pi/3+delta) + 2*cos(pi/3+delta) = 3*cos(pi/3+delta)
    #   n2 = sqrt(2)*3*cos(pi/3+delta)/sqrt(3) = sqrt(6)*cos(pi/3+delta)
    # And n1 = sqrt(6)*sin(pi/3+delta)
    # Therefore: theta = arctan2(n2, n1) = arctan2(cos(pi/3+delta), sin(pi/3+delta))
    #          = pi/2 - (pi/3 + delta)
    #          = pi/6 - delta
    # Check: theta = pi/6 - delta

    theta_predicted = np.pi/6 - delta
    # Normalize to same branch
    theta_mod = theta % (2*np.pi)
    theta_pred_mod = theta_predicted % (2*np.pi)
    print(f"    Analytic: theta = pi/6 - delta = {np.degrees(theta_predicted):.4f} deg")
    print(f"    theta (mod 2pi) = {np.degrees(theta_mod):.4f} deg")
    print(f"    pi/6 - delta (mod 2pi) = {np.degrees(theta_pred_mod):.4f} deg")
    diff = ((theta - theta_predicted + np.pi) % (2*np.pi)) - np.pi
    print(f"    Residual (theta - (pi/6-delta)) mod 2pi: {np.degrees(diff):.2e} deg")
    print()

    r['theta'] = theta
    r['n_hat'] = n_hat
    r['n_norm'] = n_norm


# ═══════════════════════════════════════════════════════════════════════════
# Part 3: Seed Yukawa matrices and bloom rotations
# ═══════════════════════════════════════════════════════════════════════════

print("\n" + "=" * 74)
print("PART 3: Seed Yukawa matrices at delta = 3*pi/4 and bloom rotations")
print("=" * 74)
print()
print(f"Seed angle: delta_seed = 3*pi/4 = {np.degrees(DELTA_SEED):.1f} deg")
print(f"At this angle: 1 + sqrt(2)*cos(delta_seed) = 1 + sqrt(2)*(-1/sqrt(2)) = 0")
print(f"so the k=0 mass vanishes: m_0 = 0.")
print()

for r in results:
    name = r['name']
    M0 = r['M0']
    delta = r['delta']
    y0 = r['y0']

    Ddelta = ((delta - DELTA_SEED) + np.pi) % (2.0 * np.pi) - np.pi

    # Seed: at delta_seed, with the same M0
    seed_m, seed_sv = seed_masses_at_delta(DELTA_SEED, M0)
    seed_y = [yukawa_from_mass(m) for m in seed_m]

    # Seed Cartan angle
    theta_seed = np.pi/6 - DELTA_SEED  # from the relation theta = pi/6 - delta

    print("-" * 74)
    print(f"Triple: {name}")
    print()
    print(f"  Physical delta = {delta:.8f} rad = {np.degrees(delta):.4f} deg")
    print(f"  Seed delta     = {DELTA_SEED:.8f} rad = {np.degrees(DELTA_SEED):.4f} deg")
    print(f"  Bloom rotation = Delta_delta = {Ddelta:.8f} rad = {np.degrees(Ddelta):.4f} deg")
    print()
    print(f"  Seed masses (MeV): [{seed_m[0]:.6e}, {seed_m[1]:.6f}, {seed_m[2]:.6f}]")
    print(f"  Seed Yukawa:       [{seed_y[0]:.6e}, {seed_y[1]:.6e}, {seed_y[2]:.6e}]")
    print()
    print(f"  Seed sqrt(Y) = diag(0, sqrt(y_1^seed), sqrt(y_2^seed))")
    print(f"               = diag(0, {np.sqrt(seed_y[1]):.6e}, {np.sqrt(seed_y[2]):.6e})")
    print()
    print(f"  Seed Cartan angle: theta_seed = pi/6 - 3pi/4 = -7pi/12 = {np.degrees(theta_seed):.2f} deg")
    print(f"  Phys Cartan angle: theta_phys = pi/6 - delta = {np.degrees(np.pi/6 - delta):.4f} deg")
    print(f"  Bloom in Cartan plane: Delta_theta = theta_phys - theta_seed = -Delta_delta")
    print(f"    = {np.degrees(-(Ddelta)):.4f} deg  (sign flip: theta rotates opposite to delta)")
    print()

    # Yukawa matrix
    print(f"  Diagonal Yukawa matrix Y = diag(y_0, y_1, y_2):")
    for k, (lbl, yi) in enumerate(zip(r['labels'], r['y'])):
        print(f"    Y[{k},{k}] = y_{lbl} = {yi:.6e}")
    print()

    r['Ddelta'] = Ddelta
    r['seed_m'] = seed_m
    r['seed_y'] = seed_y
    r['theta_seed'] = theta_seed


# ═══════════════════════════════════════════════════════════════════════════
# Part 4: Bloom rotation as Cartan plane rotation
# ═══════════════════════════════════════════════════════════════════════════

print("\n" + "=" * 74)
print("PART 4: Bloom rotation = Cartan plane rotation")
print("=" * 74)
print()
print("The Koide parametrization at fixed M0:")
print("  sqrt(y_k) = sqrt(y0) * (1 + sqrt(2)*cos(2*pi*k/3 + delta))")
print()
print("The Cartan plane angle is theta = pi/6 - delta.")
print("A bloom delta -> delta + Delta shifts theta -> theta - Delta.")
print("This is a rigid rotation in the (lambda3, lambda8) plane by angle -Delta.")
print()
print("Writing sqrt(Y) = sqrt(y0) * (I + sqrt(6)/2 * (n1*lambda3 + n2*lambda8)),")
print("the bloom acts as:")
print("  (n1, n2) -> R(-Delta) * (n1, n2)")
print("where R(alpha) is the 2D rotation matrix by angle alpha.")
print()

for r in results:
    name = r['name']
    Dd = r['Ddelta']
    delta = r['delta']
    M0 = r['M0']
    y0 = r['y0']

    # Rotation matrix in Cartan plane
    alpha = -Dd  # theta shifts by -Ddelta
    R = np.array([[np.cos(alpha), -np.sin(alpha)],
                  [np.sin(alpha),  np.cos(alpha)]])

    print("-" * 74)
    print(f"Triple: {name}")
    print(f"  Bloom angle Delta_delta = {np.degrees(Dd):.4f} deg")
    print(f"  Cartan rotation angle   = {np.degrees(alpha):.4f} deg")
    print()
    print(f"  R(-Delta_delta) = [[{R[0,0]:+.8f}  {R[0,1]:+.8f}]")
    print(f"                     [{R[1,0]:+.8f}  {R[1,1]:+.8f}]]")
    print()

    # Verify: apply to seed n_hat, should get physical n_hat
    theta_seed = r['theta_seed']
    n_seed = np.array([np.cos(theta_seed), np.sin(theta_seed)])
    n_phys_predicted = R @ n_seed
    theta_phys_predicted = np.arctan2(n_phys_predicted[1], n_phys_predicted[0])
    theta_phys_actual = r['theta']

    print(f"  Seed Cartan direction:  n_seed = ({n_seed[0]:+.8f}, {n_seed[1]:+.8f})")
    print(f"  R * n_seed = ({n_phys_predicted[0]:+.8f}, {n_phys_predicted[1]:+.8f})")
    print(f"  Physical n_hat = ({r['n_hat'][0]:+.8f}, {r['n_hat'][1]:+.8f})")
    print(f"  theta predicted: {np.degrees(theta_phys_predicted):.4f} deg")
    print(f"  theta actual:    {np.degrees(theta_phys_actual):.4f} deg")
    diff_theta = ((theta_phys_predicted - theta_phys_actual + np.pi) % (2*np.pi)) - np.pi
    print(f"  Residual: {np.degrees(diff_theta):.2e} deg")
    print()

    # Masses as functions of (M0, Delta_delta)
    print(f"  Masses from bloom: m_k(M0, Delta_delta) = M0*(1 + sqrt(2)*cos(2*pi*k/3 + 3*pi/4 + Delta))^2")
    print(f"  With M0 = {M0:.4f} MeV, Delta = {Dd:.6f} rad:")
    for k in range(3):
        c_k = 1.0 + SQRT2 * np.cos(TWO_PI_OVER_3 * k + DELTA_SEED + Dd)
        m_k = M0 * c_k ** 2
        print(f"    k={k}: c_k = {c_k:.8f},  m_k = {m_k:.4f} MeV  (actual: {r['masses'][k]:.4f} MeV)")
    print()


# ═══════════════════════════════════════════════════════════════════════════
# Part 5: Numerical verification
# ═══════════════════════════════════════════════════════════════════════════

print("\n" + "=" * 74)
print("PART 5: Numerical verification of Koide condition on Yukawa couplings")
print("=" * 74)
print()
print(f"All masses (MeV):")
print(f"  m_e = {m_e},  m_mu = {m_mu},  m_tau = {m_tau}")
print(f"  m_s = {m_s},  m_c = {m_c},  m_b = {m_b},  m_t = {m_t}")
print(f"  v = {v_higgs} GeV")
print()

for r in results:
    name = r['name']
    y = r['y']
    signed_sqrt_y = r['signed_sqrt_y']
    Q_y = r['Q_y']

    sum_y = sum(y)
    sum_signed_sqrt_y = sum(signed_sqrt_y)

    print(f"Triple: {name}")
    print(f"  Yukawa couplings:       [{', '.join(f'{yi:.6e}' for yi in y)}]")
    print(f"  Signed sqrt(y):         [{', '.join(f'{s:+.6e}' for s in signed_sqrt_y)}]")
    print(f"  Sum y_k:                {sum_y:.6e}")
    print(f"  (Sum signed sqrt(y))^2: {sum_signed_sqrt_y**2:.6e}")
    print(f"  Q = Sum(y) / (Sum signed sqrt(y))^2 = {Q_y:.10f}")
    print(f"  2/3 = {2/3:.10f}")
    print(f"  |Q - 2/3| = {abs(Q_y - 2/3):.2e}")
    print()
    print(f"  Verification: Q_mass = {r['Q']:.10f}  (should equal Q_yukawa = {Q_y:.10f})")
    print(f"  They are identical because y_k = (2/v)*m_k scales out:")
    print(f"    Q = Sum(y) / (Sum sqrt(y))^2 = Sum((2/v)*m) / (Sum sqrt(2/v)*sqrt(m))^2")
    print(f"      = (2/v)*Sum(m) / ((2/v)*(Sum sqrt(m))^2) = Sum(m)/(Sum sqrt(m))^2")
    print()

# Final summary table
print("=" * 74)
print("SUMMARY TABLE")
print("=" * 74)
print()

hdr = f"{'Triple':<14} {'M0(MeV)':>10} {'y0':>12} {'delta(deg)':>12} {'Ddelta(deg)':>12} {'theta(deg)':>12} {'Q':>12}"
print(hdr)
print("-" * len(hdr))
for r in results:
    print(f"{r['name']:<14} {r['M0']:>10.4f} {r['y0']:>12.4e} "
          f"{np.degrees(r['delta']):>12.4f} {np.degrees(r['Ddelta']):>12.4f} "
          f"{np.degrees(r['theta']):>12.4f} {r['Q']:>12.8f}")

print()
print("Key identity: y0 = 2*M0/v (exact)")
print("Key identity: theta = pi/6 - delta (exact)")
print("Key identity: Q(y) = Q(m) = 2/3 (signed, exact for Koide parametrization)")
print("Bloom rotation Delta_delta acts as rigid rotation by -Delta_delta in Cartan plane.")
print()

# Yukawa matrix table
print("Diagonal Yukawa matrices Y = diag(y_0, y_1, y_2):")
print()
for r in results:
    name = r['name']
    y = r['y']
    labels = r['labels']
    print(f"  {name}:")
    print(f"    Y = diag({', '.join(f'{yi:.4e}' for yi in y)})")
    print(f"      = diag(y_{labels[0]}, y_{labels[1]}, y_{labels[2]})")
    print(f"    y_t/y_b = {172760.0/4180.0:.2f}" if 't' in labels else "", end="")
    print()

print()
print("Done.")


# ═══════════════════════════════════════════════════════════════════════════
# Generate Markdown report
# ═══════════════════════════════════════════════════════════════════════════

lines = []
def md(s=""): lines.append(s)

md("# Yukawa Coupling Matrices from the Koide Parametrization (tan beta = 1)")
md()
md("## Setup")
md()
md("The Koide parametrization writes three signed square roots as:")
md()
md("    sqrt(m_k) = sqrt(M0) * (1 + sqrt(2) * cos(2*pi*k/3 + delta)),  k = 0, 1, 2")
md()
md("where M0 = (sum signed sqrt(m_k))^2 / 9 and delta is the Koide phase.")
md()
md("In a two-Higgs-doublet model with tan beta = v_u/v_d = 1,")
md("both vevs equal v/sqrt(2) and the Yukawa couplings are")
md()
md("    y_k = 2 m_k / v")
md()
md("for all fermion types (up, down, charged lepton). Here v = 246.22 GeV.")
md()
md("The Koide parametrization transfers to Yukawa couplings:")
md()
md("    sqrt(y_k) = sqrt(y0) * (1 + sqrt(2) * cos(2*pi*k/3 + delta))")
md()
md("with y0 = 2*M0/v. The Koide phase delta is identical for masses and Yukawa couplings.")
md()
md("PDG masses used (MeV): m_e = 0.51100, m_mu = 105.658, m_tau = 1776.86,")
md("m_s = 93.4, m_c = 1270.0, m_b = 4180.0, m_t = 172760.0.")
md()

# Part 1
md("## Koide Parameters and Yukawa Couplings")
md()

for r in results:
    name = r['name']
    md(f"### {name}")
    md()
    md("| Parameter | Value |")
    md("|-----------|-------|")
    md(f"| M0 (MeV) | {r['M0']:.6f} |")
    md(f"| delta (rad) | {r['delta']:.8f} |")
    md(f"| delta (deg) | {np.degrees(r['delta']):.4f} |")
    md(f"| Q (signed) | {r['Q']:.10f} |")
    md(f"| y0 = 2*M0/v | {r['y0']:.6e} |")
    md()
    md("| Fermion | m (MeV) | y = 2m/v | sqrt(y) (signed) |")
    md("|---------|---------|----------|-------------------|")
    for lbl, m, yi, ssy in zip(r['labels'], r['masses'], r['y'], r['signed_sqrt_y']):
        md(f"| {lbl} | {m:.4f} | {yi:.6e} | {ssy:+.6e} |")
    md()
    md(f"Verification: (sum signed sqrt(y))^2 / 9 = {sum(r['signed_sqrt_y'])**2/9:.6e} = y0 = {r['y0']:.6e}")
    md()

# Part 2
md("## Cartan Subalgebra Decomposition")
md()
md("The diagonal matrix sqrt(Y) = diag(sqrt(y_0), sqrt(y_1), sqrt(y_2)) decomposes as:")
md()
md("    sqrt(Y) = sqrt(y0) * (I + (sqrt(6)/2) * (n1 * lambda3 + n2 * lambda8))")
md()
md("where lambda3 = diag(1, -1, 0) and lambda8 = diag(1, 1, -2)/sqrt(3) are")
md("the diagonal Gell-Mann matrices, and n_hat = (n1, n2) = (cos theta, sin theta)")
md("is a unit vector in the Cartan plane.")
md()
md("The Cartan angle theta is related to the Koide phase by an exact identity:")
md()
md("    theta = pi/6 - delta")
md()
md("This follows from the trigonometric identities:")
md()
md("    n1 propto cos(delta) - cos(2*pi/3 + delta) = sqrt(3) * sin(pi/3 + delta)")
md("    n2 propto cos(delta) + cos(2*pi/3 + delta) - 2*cos(4*pi/3 + delta) = 3*cos(pi/3 + delta)")
md()
md("so theta = arctan2(cos(pi/3 + delta), sin(pi/3 + delta)) = pi/2 - pi/3 - delta = pi/6 - delta.")
md()

for r in results:
    name = r['name']
    md(f"### {name}")
    md()
    md("| Cartan parameter | Value |")
    md("|------------------|-------|")
    md(f"| n_hat | ({r['n_hat'][0]:+.6f}, {r['n_hat'][1]:+.6f}) |")
    md(f"| n_norm | {r['n_norm']:.8f} |")
    md(f"| theta (deg) | {np.degrees(r['theta']):.4f} |")
    md(f"| pi/6 - delta (deg) | {np.degrees(np.pi/6 - r['delta']):.4f} |")
    md()

# Part 3
md("## Seed Yukawa Matrices and Bloom Rotations")
md()
md("At delta_seed = 3*pi/4 = 135 deg, the k=0 coefficient vanishes:")
md("1 + sqrt(2)*cos(3*pi/4) = 1 - 1 = 0, so m_0 = y_0 = 0.")
md()
md("The bloom rotation Delta_delta = delta_phys - 3*pi/4 parametrizes how far")
md("the physical triple has rotated from the seed configuration.")
md()

for r in results:
    name = r['name']
    md(f"### {name}")
    md()
    md("| Bloom parameter | Value |")
    md("|-----------------|-------|")
    md(f"| delta_phys (deg) | {np.degrees(r['delta']):.4f} |")
    md(f"| delta_seed (deg) | 135.0000 |")
    md(f"| Delta_delta (deg) | {np.degrees(r['Ddelta']):.4f} |")
    md(f"| Delta_delta (rad) | {r['Ddelta']:.8f} |")
    md()
    md(f"Seed Yukawa couplings at delta = 3*pi/4 (with same M0 = {r['M0']:.4f} MeV):")
    md()
    md("| k | Seed m (MeV) | Seed y | Physical m (MeV) | Physical y |")
    md("|---|--------------|--------|------------------|------------|")
    for k in range(3):
        md(f"| {k} | {r['seed_m'][k]:.4e} | {r['seed_y'][k]:.4e} | {r['masses'][k]:.4f} | {r['y'][k]:.4e} |")
    md()

# Part 4
md("## Bloom as Cartan Plane Rotation")
md()
md("Since theta = pi/6 - delta, shifting delta by Delta_delta shifts theta by -Delta_delta.")
md("In the Cartan plane, the bloom acts as a rigid rotation:")
md()
md("    (n1, n2) -> R(-Delta_delta) * (n1, n2)")
md()
md("where R(alpha) = [[cos alpha, -sin alpha], [sin alpha, cos alpha]].")
md()
md("Equivalently, the traceless part of sqrt(Y) rotates in the (lambda3, lambda8) plane")
md("without changing its norm (which is fixed at sqrt(6)/2 for Q = 2/3) or the trace")
md("(which is 3*sqrt(y0) and set by M0).")
md()
md("The three masses after a bloom of Delta_delta from the seed are:")
md()
md("    m_k(M0, Delta) = M0 * (1 + sqrt(2) * cos(2*pi*k/3 + 3*pi/4 + Delta))^2")
md()
md("At Delta = 0 (seed): m_0 = 0, m_1 = M0*(1+sqrt(2)*cos(2*pi/3+3*pi/4))^2,")
md("m_2 = M0*(1+sqrt(2)*cos(4*pi/3+3*pi/4))^2.")
md()
md("Summary of bloom rotations:")
md()
md("| Triple | Delta_delta (deg) | Delta_delta (rad) | Cartan rotation (deg) |")
md("|--------|-------------------|-------------------|-----------------------|")
for r in results:
    md(f"| {r['name']} | {np.degrees(r['Ddelta']):.4f} | {r['Ddelta']:.6f} | {np.degrees(-r['Ddelta']):.4f} |")
md()

# Part 5
md("## Numerical Verification")
md()
md("The Koide condition Q = sum(y_k) / (sum signed sqrt(y_k))^2 = 2/3")
md("is algebraically identical to Q = sum(m_k) / (sum signed sqrt(m_k))^2 = 2/3")
md("because y_k = (2/v)*m_k and the factor (2/v) cancels between numerator and denominator:")
md()
md("    Q = sum((2/v)*m) / (sum(sqrt(2/v)*sqrt(m)))^2 = (2/v)*sum(m) / ((2/v)*(sum sqrt(m))^2) = sum(m)/(sum sqrt(m))^2")
md()

md("| Triple | Q (signed) | |Q - 2/3| |")
md("|--------|------------|-----------|")
for r in results:
    md(f"| {r['name']} | {r['Q_y']:.10f} | {abs(r['Q_y'] - 2/3):.2e} |")
md()
md("For (e, mu, tau), Q is 2/3 to 6e-6 (0.001%), confirming the near-exact Koide relation.")
md("For (-s, c, b) and (c, b, t), the deviations (1.2% and 0.4%) reflect how far the")
md("PDG masses lie from the nearest Koide curve. The algebraic identity Q_signed = 2/3")
md("holds exactly for the *fitted* parametrization values, not for the raw PDG masses.")
md()

# Summary
md("## Summary")
md()
md("| Triple | M0 (MeV) | y0 | delta (deg) | Delta_delta (deg) | theta (deg) |")
md("|--------|----------|-----|-------------|-------------------|-------------|")
for r in results:
    md(f"| {r['name']} | {r['M0']:.4f} | {r['y0']:.4e} | {np.degrees(r['delta']):.4f} | {np.degrees(r['Ddelta']):.4f} | {np.degrees(r['theta']):.4f} |")
md()

md("Key results:")
md()
md("1. The Koide parametrization transfers directly from masses to Yukawa couplings")
md("   via y0 = 2*M0/v. The phase delta is identical.")
md()
md("2. The Cartan subalgebra decomposition of sqrt(Y) has an exact relationship")
md("   theta = pi/6 - delta between the Cartan plane angle and the Koide phase.")
md()
md("3. The bloom from the seed (delta = 3*pi/4, one mass zero) to the physical")
md("   triple is a rigid rotation in the Cartan plane by angle -Delta_delta.")
md()
md("4. The Q = 2/3 condition is invariant under the rescaling m -> y = 2m/v,")
md("   as it must be since Q is a ratio of homogeneous-degree-1 expressions.")
md()
md("5. At tan beta = 1 the top Yukawa coupling is y_t = 2*m_t/v = 1.403,")
md(f"   which is O(1) as expected. The full range spans y_e = {results[0]['y'][0]:.2e}")
md(f"   to y_t = {results[2]['y'][2]:.4f}, a hierarchy of {results[2]['y'][2]/results[0]['y'][0]:.0e}.")

report = "\n".join(lines)
with open("/home/codexssh/phys3/results/yukawa_koide.md", "w") as f:
    f.write(report)

print("\n\nMarkdown report written to /home/codexssh/phys3/results/yukawa_koide.md")
