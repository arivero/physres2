"""
Pseudomodulus stabilization via Coleman-Weinberg potential
in the O'Raifeartaigh model: W = f*X + m*phi*phit + g*X*phi^2

Computes tree-level spectrum, one-loop CW potential, and minimum.
"""

import numpy as np
from scipy.optimize import minimize_scalar
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import warnings
warnings.filterwarnings('ignore')

# ---------------------------------------------------------------------------
# Tree-level mass spectrum at phi = phit = 0
# ---------------------------------------------------------------------------
# Fields: X (pseudo-modulus), phi, phit (chiral superfields)
# Dimensionless variable: x = gX/m
# SUSY-breaking parameter: y = gf/m^2
#
# Superpotential second derivatives at minimum (phi=phit=0):
#   W_XX = 0
#   W_Xphi = 0
#   W_Xphit = 0
#   W_phiphi = 2gX  =>  2m*x  (in units of m)
#   W_phiphit = m
#   W_phitphit = 0
#
# Fermionic mass matrix in (phi, phit) block (X is Goldstino, mass=0):
#   M_ferm_2x2 = [[2gX, m], [m, 0]]  = m * [[2x, 1], [1, 0]]
#
# Fermionic masses: eigenvalues of M_ferm_2x2
#   lambda^2 - 2gX*lambda - m^2 = 0
#   lambda = gX +/- sqrt(g^2 X^2 + m^2) = m*(x +/- sqrt(x^2+1))
#   masses: |m*(x + sqrt(x^2+1))| and |m*(x - sqrt(x^2+1))|
#   Since sqrt(x^2+1) > |x|, both are real. The second is negative for x>0:
#   mass_f1 = m*(x + sqrt(x^2+1))
#   mass_f2 = m*(sqrt(x^2+1) - x)   [always positive]
#
# Scalar mass-squared matrix (4x4 in Re/Im phi, Re/Im phit basis):
# The complex scalar kinetic/mass structure comes from:
#   V = |f + g phi^2|^2 + |m phit + 2gX phi|^2 + |m phi|^2
#
# At phi=phit=0, expanding to quadratic order in fluctuations:
# Writing phi = a + ib, phit = c + id (all real), the mass-squared matrix is:
#
# From |F_X|^2 = |f + g phi^2|^2:
#   Expanding: f^2 + 2f*g*(a^2-b^2) + ... => contributes +2fg to (a,a) and -2fg to (b,b)
#   i.e. +2gf to Re(phi)^2, -2gf to Im(phi)^2
#
# From |F_phi|^2 = |m phit + 2gX phi|^2:
#   = m^2(c^2+d^2) + 4g^2X^2(a^2+b^2) + 4gXm(ac+bd)
#   => 4g^2X^2 for (a,a) and (b,b); m^2 for (c,c) and (d,d); 4gXm for (a,c) and (b,d) cross
#
# From |F_phit|^2 = |m phi|^2 = m^2(a^2+b^2):
#   => m^2 for (a,a) and (b,b)
#
# Full 4x4 matrix in basis (Re phi, Re phit, Im phi, Im phit):
# (Re phi and Re phit mix; Im phi and Im phit mix separately)
#
# Re sector (2x2):
#   [[4g^2X^2 + m^2 + 2gf,  2gXm ],
#    [2gXm,                  m^2  ]]
#
# Im sector (2x2):
#   [[4g^2X^2 + m^2 - 2gf,  2gXm ],
#    [2gXm,                  m^2  ]]

def fermion_masses(x):
    """
    Dimensionless fermion masses (in units of m).
    x = gX/m.
    Returns (mf1, mf2) where mf1 >= mf2 > 0. Goldstino is massless.
    """
    sq = np.sqrt(x**2 + 1.0)
    mf1 = x + sq   # always positive
    mf2 = sq - x   # always positive (sq > x for x > 0)
    # For x < 0, mf2 = sq - x > sq > 0 still, mf1 = x + sq could be small
    # Both always positive since sq = sqrt(x^2+1) > |x|
    return np.abs(mf1), np.abs(mf2)


def boson_masses_sq(x, y):
    """
    Dimensionless boson mass-squared values (in units of m^2).
    x = gX/m, y = gf/m^2.
    Returns array of 4 values: [mb1_sq, mb2_sq, mb3_sq, mb4_sq].
    Also includes mass-squared of X (=0, it's flat at tree level).

    The (phi, phit) sector gives 4 real scalars:
    Re-sector eigenvalues (2x2 matrix):
      A = [[4x^2+1+2y, 2x], [2x, 1]]
    Im-sector eigenvalues (2x2 matrix):
      B = [[4x^2+1-2y, 2x], [2x, 1]]

    Plus the X scalar which is flat: m^2_X = 0.
    """
    # Re-sector
    A = np.array([[4*x**2 + 1 + 2*y, 2*x],
                  [2*x,               1.0]])
    eig_re = np.linalg.eigvalsh(A)

    # Im-sector
    B = np.array([[4*x**2 + 1 - 2*y, 2*x],
                  [2*x,               1.0]])
    eig_im = np.linalg.eigvalsh(B)

    return np.concatenate([eig_re, eig_im])


def supertrace(x, y):
    """
    STr[M^2] = sum(-1)^{2J}(2J+1) m^2
    = (bosons m^2) - 2*(fermions m^2)
    In units of m^2.
    Goldstino and X contribute: Goldstino has m=0, X scalar has m=0.
    """
    bsq = boson_masses_sq(x, y)  # 4 values
    mf1, mf2 = fermion_masses(x)

    str_m2 = np.sum(bsq) - 2*(mf1**2 + mf2**2)
    # X boson contributes 0, Goldstino contributes 0
    return str_m2


def supertrace_m4(x, y):
    """
    STr[M^4] = sum(-1)^{2J}(2J+1) m^4
    """
    bsq = boson_masses_sq(x, y)
    mf1, mf2 = fermion_masses(x)

    str_m4 = np.sum(bsq**2) - 2*(mf1**4 + mf2**4)
    return str_m4


def cw_potential(x, y, include_X_scalar=True):
    """
    One-loop Coleman-Weinberg potential in units of m^4/(64 pi^2).

    V_CW = (1/64pi^2) * STr[M^4 (ln(M^2/mu^2) - 3/2)]

    Here we use a renormalization scale choice mu^2 = m^2 (so ln(M^2/mu^2) = ln(m^2_dimensionless)).
    We work in units where V is in units of m^4/(64pi^2).

    So V_CW (in units of m^4/(64pi^2)) =
        sum_bosons m^4_b * (ln(m^2_b) - 3/2)
        - 2 * sum_fermions m^4_f * (ln(m^2_f) - 3/2)

    where masses are dimensionless (in units of m).
    """
    bsq = boson_masses_sq(x, y)
    mf1, mf2 = fermion_masses(x)

    V = 0.0
    # Scalar contributions (factor 1 each for complex -> 2 real, but we have 4 real scalars)
    for msq in bsq:
        if msq > 1e-30:
            V += msq**2 * (np.log(msq) - 1.5)
        # If msq <= 0 or zero: massless scalar contributes 0 (or signals tachyon)

    # If including X scalar (massless at tree level), it contributes 0
    # (it's the pseudo-modulus whose mass we're computing)

    # Fermionic contributions (factor of -2 for Dirac, but these are Weyl)
    # Each Weyl fermion contributes with factor -1 in the CW potential
    # The two non-zero Weyl fermions with masses mf1, mf2:
    for mf in [mf1, mf2]:
        if mf > 1e-15:
            mfsq = mf**2
            V -= mfsq**2 * (np.log(mfsq) - 1.5)

    # Goldstino: massless, contributes 0

    return V


def find_cw_minimum(y, x_range=(-5, 5), n_points=1000):
    """
    Find the minimum of V_CW(x) for given y.
    Returns (x_min, V_min, success).
    """
    x_vals = np.linspace(x_range[0], x_range[1], n_points)
    V_vals = np.array([cw_potential(xi, y) for xi in x_vals])

    # Find approximate minimum
    # Check for tachyons (negative mass-squared)
    valid = np.array([all(boson_masses_sq(xi, y) > -1e-10) for xi in x_vals])

    # The CW potential can have minimum at x=0 or x>0
    # We look for the global minimum ignoring tachyonic regions
    idx_min = np.argmin(V_vals)
    x_approx = x_vals[idx_min]

    # Refine with scipy
    try:
        result = minimize_scalar(
            lambda xi: cw_potential(xi, y),
            bounds=(x_range[0], x_range[1]),
            method='bounded',
            options={'xatol': 1e-10}
        )
        x_min = result.x
        V_min = result.fun
        return x_min, V_min, result.success, V_vals, x_vals
    except Exception as e:
        return x_approx, V_vals[idx_min], False, V_vals, x_vals


def koide_Q(masses):
    """
    Koide ratio Q = (sum m) / (sum sqrt(m))^2
    masses: array of 3 positive values
    """
    m = np.array(masses)
    if np.any(m <= 0):
        return np.nan
    Q = np.sum(m) / (np.sum(np.sqrt(m)))**2
    return Q


# ---------------------------------------------------------------------------
# Main computation
# ---------------------------------------------------------------------------

print("=" * 70)
print("O'Raifeartaigh Model: Pseudomodulus Stabilization via CW Potential")
print("=" * 70)
print()

# Part 1: Tree-level spectrum
print("PART 1: Tree-level mass spectrum at phi = phit = 0")
print("-" * 70)
print()
print("Dimensionless fermion masses (in units of m), as functions of x = gX/m:")
print("  mf1 = |x + sqrt(x^2 + 1)|")
print("  mf2 = |sqrt(x^2 + 1) - x|")
print()
print(f"{'x':>8} | {'mf1':>12} | {'mf2':>12} | {'mf1*mf2':>12} | {'mf1+mf2':>12}")
print("-" * 60)
for x in [-2, -1, -0.5, 0, 0.5, 1, np.sqrt(3), 2, 3]:
    mf1, mf2 = fermion_masses(x)
    print(f"{x:8.4f} | {mf1:12.6f} | {mf2:12.6f} | {mf1*mf2:12.6f} | {mf1+mf2:12.6f}")
print()
print("Note: mf1 * mf2 = 1 (product of eigenvalues = det = -m^2, so |det| = m^2 = 1 in units)")
print()

# Part 2: Boson spectrum for a specific y
y_test = 1.0
print(f"Boson mass-squared (in units of m^2) for y = gf/m^2 = {y_test}:")
print(f"{'x':>8} | {'mb1^2':>10} | {'mb2^2':>10} | {'mb3^2':>10} | {'mb4^2':>10}")
print("-" * 55)
for x in [0, 0.5, 1.0, np.sqrt(3), 2.0]:
    bsq = boson_masses_sq(x, y_test)
    bsq_sorted = sorted(bsq)
    print(f"{x:8.4f} | {bsq_sorted[0]:10.5f} | {bsq_sorted[1]:10.5f} | {bsq_sorted[2]:10.5f} | {bsq_sorted[3]:10.5f}")
print()

# Part 5: Supertrace check
print("PART 5: Supertrace verification")
print("-" * 70)
print()
print("STr[M^2] should = 0 for O'Raifeartaigh model (canonical Kahler)")
print()
y_vals = [0.1, 0.5, 1.0, 2.0]
x_check = [0, 1.0, np.sqrt(3)]
print(f"{'y':>6} | {'x':>8} | {'STr[M^2]':>14} | {'STr[M^4]':>14}")
print("-" * 50)
for y in y_vals:
    for x in x_check:
        str2 = supertrace(x, y)
        str4 = supertrace_m4(x, y)
        print(f"{y:6.2f} | {x:8.4f} | {str2:14.6f} | {str4:14.6f}")
print()

# Parts 2, 3, 4: CW potential and minimum
print("PARTS 2-4: CW potential and minimum")
print("-" * 70)
print()

results = {}
fig, axes = plt.subplots(2, 2, figsize=(12, 10))
axes = axes.flatten()

for idx, y in enumerate(y_vals):
    x_min, V_min, success, V_vals, x_vals = find_cw_minimum(y, x_range=(0, 6))
    results[y] = {
        'x_min': x_min,
        'V_min': V_min,
        'success': success,
        'V_vals': V_vals,
        'x_vals': x_vals
    }

    # Also check x=0
    V_at_0 = cw_potential(0, y)

    ax = axes[idx]
    ax.plot(x_vals, V_vals, 'b-', linewidth=2, label='V_CW(x)')
    ax.axvline(x=x_min, color='r', linestyle='--', alpha=0.7, label=f'x_min = {x_min:.4f}')
    ax.axvline(x=np.sqrt(3), color='g', linestyle=':', alpha=0.7, label=f'√3 = {np.sqrt(3):.4f}')
    ax.scatter([x_min], [V_min], color='red', s=100, zorder=5)
    ax.set_xlabel('x = gX/m', fontsize=12)
    ax.set_ylabel('V_CW [units of m⁴/(64π²)]', fontsize=11)
    ax.set_title(f'y = gf/m² = {y}', fontsize=13)
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3)
    ax.set_xlim(0, 5)

plt.suptitle("Coleman-Weinberg Potential for O'Raifeartaigh Pseudomodulus", fontsize=14)
plt.tight_layout()
plt.savefig('/home/codexssh/phys3/results/pseudomodulus_cw_potential.png', dpi=150, bbox_inches='tight')
plt.close()
print("Saved: pseudomodulus_cw_potential.png")
print()

# Summary table
print(f"{'y':>6} | {'x_min':>10} | {'x_min/√3':>10} | {'V_CW(x_min)':>14} | {'V_CW(0)':>14}")
print("-" * 62)
sqrt3 = np.sqrt(3)
for y in y_vals:
    r = results[y]
    V0 = cw_potential(0, y)
    print(f"{y:6.2f} | {r['x_min']:10.6f} | {r['x_min']/sqrt3:10.6f} | {r['V_min']:14.6f} | {V0:14.6f}")
print()

# More detailed scan over y
print("Fine scan over y:")
y_scan = np.linspace(0.05, 3.0, 60)
x_mins = []
x_min_over_sqrt3 = []
for y in y_scan:
    xm, vm, ok, _, _ = find_cw_minimum(y, x_range=(0, 8))
    x_mins.append(xm)
    x_min_over_sqrt3.append(xm / sqrt3)

# Plot x_min vs y
fig2, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

ax1.plot(y_scan, x_mins, 'b-', linewidth=2)
ax1.axhline(y=sqrt3, color='r', linestyle='--', label=f'√3 ≈ {sqrt3:.4f}')
ax1.set_xlabel('y = gf/m²', fontsize=13)
ax1.set_ylabel('x_min = gX_min/m', fontsize=13)
ax1.set_title('CW minimum vs SUSY-breaking parameter', fontsize=13)
ax1.legend(fontsize=11)
ax1.grid(True, alpha=0.3)

ax2.plot(y_scan, x_min_over_sqrt3, 'b-', linewidth=2)
ax2.axhline(y=1.0, color='r', linestyle='--', label='x_min/√3 = 1')
ax2.set_xlabel('y = gf/m²', fontsize=13)
ax2.set_ylabel('x_min / √3', fontsize=13)
ax2.set_title('CW minimum in units of √3', fontsize=13)
ax2.legend(fontsize=11)
ax2.grid(True, alpha=0.3)

plt.suptitle("Location of CW Minimum", fontsize=14)
plt.tight_layout()
plt.savefig('/home/codexssh/phys3/results/pseudomodulus_xmin_vs_y.png', dpi=150, bbox_inches='tight')
plt.close()
print("Saved: pseudomodulus_xmin_vs_y.png")
print()

# Part 6: Koide ratio at the CW minimum
print("PART 6: Koide ratio at CW minimum")
print("-" * 70)
print()
print("At the CW minimum, compute boson masses and Koide Q = (Σm)/(Σ√m)²")
print()
print(f"{'y':>6} | {'x_min':>8} | {'mb1':>10} | {'mb2':>10} | {'mb3':>10} | {'mb4':>10} | {'Q(123)':>10} | {'Q(234)':>10}")
print("-" * 90)

Q_vals = []
for y in y_vals:
    xm = results[y]['x_min']
    bsq = boson_masses_sq(xm, y)
    bsq_sorted = sorted(bsq)

    # Compute boson masses (positive square roots, skip any zero/negative)
    b_masses = []
    for bsq_val in bsq_sorted:
        if bsq_val > 1e-12:
            b_masses.append(np.sqrt(bsq_val))
        else:
            b_masses.append(0.0)

    if len(b_masses) >= 3:
        # Koide on 3 heaviest
        m_nonzero = [m for m in b_masses if m > 1e-10]
        if len(m_nonzero) >= 3:
            Q123 = koide_Q(m_nonzero[:3])
            Q234 = koide_Q(m_nonzero[1:4]) if len(m_nonzero) >= 4 else np.nan
        else:
            Q123 = np.nan
            Q234 = np.nan
    else:
        Q123 = np.nan
        Q234 = np.nan

    Q_vals.append(Q123)

    bm_str = " | ".join([f"{m:10.5f}" for m in b_masses])
    print(f"{y:6.2f} | {xm:8.4f} | {bm_str} | {Q123:10.6f} | {Q234:10.6f}")

print()
print(f"Reference: Koide Q = 2/3 = {2/3:.6f}")
print()

# Scan Q over y
print("Q-ratio scan over y (at CW minimum):")
y_scan2 = np.linspace(0.05, 3.0, 60)
Q_scan_123 = []
Q_scan_234 = []
x_scan2 = []

for y in y_scan2:
    xm, vm, ok, _, _ = find_cw_minimum(y, x_range=(0, 8))
    x_scan2.append(xm)
    bsq = boson_masses_sq(xm, y)
    bsq_sorted = sorted(bsq)
    b_masses = [np.sqrt(max(v, 0)) for v in bsq_sorted]
    m_nonzero = [m for m in b_masses if m > 1e-10]
    if len(m_nonzero) >= 3:
        Q_scan_123.append(koide_Q(m_nonzero[:3]))
        if len(m_nonzero) >= 4:
            Q_scan_234.append(koide_Q(m_nonzero[1:4]))
        else:
            Q_scan_234.append(np.nan)
    else:
        Q_scan_123.append(np.nan)
        Q_scan_234.append(np.nan)

fig3, ax = plt.subplots(figsize=(10, 6))
ax.plot(y_scan2, Q_scan_123, 'b-', linewidth=2, label='Q(lightest 3 bosons)')
ax.plot(y_scan2, Q_scan_234, 'g--', linewidth=2, label='Q(2nd–4th bosons)')
ax.axhline(y=2/3, color='r', linestyle='--', linewidth=2, label='Q = 2/3')
ax.set_xlabel('y = gf/m²', fontsize=13)
ax.set_ylabel('Koide Q ratio', fontsize=13)
ax.set_title('Koide Q at CW Minimum vs SUSY-Breaking Parameter', fontsize=13)
ax.legend(fontsize=11)
ax.grid(True, alpha=0.3)
ax.set_ylim(0.3, 1.0)
plt.tight_layout()
plt.savefig('/home/codexssh/phys3/results/pseudomodulus_koide.png', dpi=150, bbox_inches='tight')
plt.close()
print("Saved: pseudomodulus_koide.png")
print()

# Detailed analysis at the O'Raifeartaigh-Koide point x = sqrt(3)
print("DETAILED ANALYSIS at x = √3 (the O'Raifeartaigh-Koide point)")
print("-" * 70)
print()
for y in y_vals:
    x = sqrt3
    mf1, mf2 = fermion_masses(x)
    bsq = boson_masses_sq(x, y)
    bsq_sorted = sorted(bsq)
    b_masses = [np.sqrt(max(v, 0)) for v in bsq_sorted]
    V = cw_potential(x, y)
    Vmin = cw_potential(results[y]['x_min'], y)

    m_nonzero = [m for m in b_masses if m > 1e-10]
    Q = koide_Q(m_nonzero[:3]) if len(m_nonzero) >= 3 else np.nan

    print(f"y = {y}:")
    print(f"  Fermion masses: mf1 = {mf1:.6f} m, mf2 = {mf2:.6f} m")
    print(f"  Boson masses:   mb1 = {b_masses[0]:.6f} m, mb2 = {b_masses[1]:.6f} m")
    print(f"                  mb3 = {b_masses[2]:.6f} m, mb4 = {b_masses[3]:.6f} m")
    print(f"  V_CW(√3) = {V:.8f}   V_CW(x_min) = {Vmin:.8f}")
    print(f"  Koide Q(3 lightest bosons) = {Q:.8f}   (2/3 = {2/3:.8f})")
    str2 = supertrace(x, y)
    str4 = supertrace_m4(x, y)
    print(f"  STr[M^2] = {str2:.6f}   STr[M^4] = {str4:.6f}")
    print()

# Check: at x=sqrt(3), what are the exact fermion mass eigenvalues?
print("Fermion masses at x = √3 exactly:")
x = sqrt3
mf1, mf2 = fermion_masses(x)
print(f"  mf1 = √3 + 2 = {sqrt3 + 2:.6f}  (computed: {mf1:.6f})")
print(f"  mf2 = 2 - √3 = {2 - sqrt3:.6f}  (computed: {mf2:.6f})")
print(f"  mf1 * mf2 = {mf1*mf2:.8f}  (should be 1)")
print(f"  mf1 + mf2 = {mf1+mf2:.8f}  (should be 4)")
print()

# Analytical check of boson masses at x=sqrt(3)
print("Boson mass-squared eigenvalues at x = √3, for various y:")
print("  Re-sector matrix: [[4*3+1+2y, 2√3], [2√3, 1]] = [[13+2y, 2√3], [2√3, 1]]")
print("  Im-sector matrix: [[13-2y, 2√3], [2√3, 1]]")
print()
for y in [0.5, 1.0, 2.0]:
    # Re sector
    A = np.array([[13 + 2*y, 2*sqrt3],
                  [2*sqrt3,  1.0]])
    eig_re = np.linalg.eigvalsh(A)
    # Im sector
    B = np.array([[13 - 2*y, 2*sqrt3],
                  [2*sqrt3,  1.0]])
    eig_im = np.linalg.eigvalsh(B)
    all_bsq = np.sort(np.concatenate([eig_re, eig_im]))
    bm = np.sqrt(np.maximum(all_bsq, 0))
    print(f"  y={y}: m^2 = [{all_bsq[0]:.4f}, {all_bsq[1]:.4f}, {all_bsq[2]:.4f}, {all_bsq[3]:.4f}]")
    print(f"         m   = [{bm[0]:.4f}, {bm[1]:.4f}, {bm[2]:.4f}, {bm[3]:.4f}]")
    m3 = bm[bm > 1e-10]
    if len(m3) >= 3:
        print(f"         Koide Q(3 lightest) = {koide_Q(m3[:3]):.8f}")
    print()

# Check whether minimum is at x=sqrt(3) for some y (gv/m = sqrt(3))
print("CHECKING: Is x_min = √3 for any y?")
print("-" * 70)
# The derivative of V_CW at x=sqrt(3) must be zero
# Let's compute dV_CW/dx numerically at x=sqrt(3) for each y
print(f"{'y':>8} | {'x_min':>10} | {'x_min-√3':>12} | {'dV/dx at √3':>14}")
print("-" * 50)
for y in np.linspace(0.1, 3.0, 20):
    xm, vm, ok, _, _ = find_cw_minimum(y, x_range=(0, 8))
    dx = 1e-7
    dV = (cw_potential(sqrt3 + dx, y) - cw_potential(sqrt3 - dx, y)) / (2*dx)
    print(f"{y:8.3f} | {xm:10.6f} | {xm - sqrt3:12.6f} | {dV:14.6f}")

print()
print("=" * 70)
print("Summary: See pseudomodulus_cw.md for full analysis")
print("=" * 70)
