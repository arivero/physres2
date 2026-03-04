"""
Seiberg effective theory spectrum for SU(3)_c SQCD.
N_f = N_c = 3 confining phase with quantum-modified moduli space.

Careful analytical + numerical treatment.
"""

import numpy as np
from itertools import combinations

# ============================================================
# Physical inputs (MS-bar at 2 GeV unless noted)
# ============================================================
Lambda = 300.0  # MeV
m_u = 2.16      # MeV
m_d = 4.67      # MeV
m_s = 93.4      # MeV
m_c = 1270.0    # MeV
m_b = 4180.0    # MeV

masses_all = {'u': m_u, 'd': m_d, 's': m_s, 'c': m_c, 'b': m_b}

def koide_Q(m1, m2, m3):
    """Standard Koide: Q = (sum sqrt(m))^2 / (3 sum m)"""
    s = np.sqrt(m1) + np.sqrt(m2) + np.sqrt(m3)
    return s**2 / (3.0 * (m1 + m2 + m3))

# ============================================================
# PART 1: SEIBERG SEESAW — N_f = N_c = 3
# ============================================================

def seiberg_seesaw(m_vec, Lambda, labels):
    """
    For N_f = N_c = 3 with quantum-modified constraint det M = Lambda^6.

    Superpotential: W = Tr(m M) + X(det M - B Bbar - Lambda^6)
    At B = Bbar = 0 vacuum:

    F-term for M^j_j: m_j + X * (det M / M_j) = 0
    F-term for X: det M = Lambda^6

    Solution: m_j * M_j = const = C for all j
    C = Lambda^2 * (prod m_j)^{1/3}
    M_j = C / m_j
    X = -C / Lambda^6 (from m_j + X * det M / M_j = 0)
      = -(prod m_j)^{1/3} / Lambda^4  ... wait let me redo.

    From m_j + X * (Lambda^6 / M_j) = 0:
    X = -m_j * M_j / Lambda^6 = -C / Lambda^6
    """
    Nf = len(m_vec)
    assert Nf == 3

    prod_m = np.prod(m_vec)
    C = Lambda**2 * prod_m**(1.0/3.0)
    M_vev = C / m_vec
    X_vev = -C / Lambda**6

    # Verify
    assert np.isclose(np.prod(M_vev), Lambda**6, rtol=1e-10), f"det M = {np.prod(M_vev)} != {Lambda**6}"
    for j in range(3):
        cofactor_j = np.prod(M_vev) / M_vev[j]
        Fj = m_vec[j] + X_vev * cofactor_j
        assert abs(Fj) < 1e-6 * abs(m_vec[j]), f"F_{j} = {Fj} nonzero"

    print(f"\n  C = Lambda^2 * (prod m)^(1/3) = {C:.6f} MeV^2")
    print(f"  X = -C / Lambda^6 = {X_vev:.6e} MeV^(-4)")
    print()
    for j in range(3):
        print(f"  M_{labels[j]} = {M_vev[j]:.4f} MeV   (m_{labels[j]}*M_{labels[j]} = {m_vec[j]*M_vev[j]:.4f})")

    print(f"\n  det M = {np.prod(M_vev):.4f},  Lambda^6 = {Lambda**6:.4f}")

    # Koide
    q_m = koide_Q(*m_vec)
    q_M = koide_Q(*M_vev)
    q_inv = koide_Q(*(1.0/m_vec))
    print(f"\n  Q(m) = {q_m:.6f}")
    print(f"  Q(M) = Q(1/m) = {q_M:.6f}")
    print(f"  Verify Q(1/m) = {q_inv:.6f}")
    print(f"  |Q(1/m) - 2/3| = {abs(q_inv - 2./3.):.6f} ({abs(q_inv - 2./3.)/(2./3.)*100:.4f}%)")

    return M_vev, X_vev, C

print("=" * 80)
print("PART 1: SEIBERG SEESAW MESON VEVS (N_f = N_c = 3)")
print("=" * 80)

for trip_labels, trip_masses in [
    (['u','d','s'], [m_u, m_d, m_s]),
    (['d','s','b'], [m_d, m_s, m_b]),
    (['s','c','b'], [m_s, m_c, m_b]),
    (['u','d','c'], [m_u, m_d, m_c]),
    (['u','d','b'], [m_u, m_d, m_b]),
]:
    print(f"\n--- Triple: ({','.join(trip_labels)}) ---")
    seiberg_seesaw(np.array(trip_masses), Lambda, trip_labels)


# ============================================================
# PART 2: FULL SPECTRUM FOR (u,d,s) — ANALYTICAL
# ============================================================
print("\n\n" + "=" * 80)
print("PART 2: FERMION MASS MATRIX — ANALYTICAL (u,d,s)")
print("=" * 80)

m_vec = np.array([m_u, m_d, m_s])
labels_3 = ['u', 'd', 's']

prod_m = np.prod(m_vec)
C = Lambda**2 * prod_m**(1.0/3.0)
M = C / m_vec  # M_1, M_2, M_3
X = -C / Lambda**6

print(f"\nVacuum: M = ({M[0]:.4f}, {M[1]:.4f}, {M[2]:.4f}) MeV")
print(f"        X = {X:.6e} MeV^(-4)")

# Fields: 9 mesons M^i_j (i=row, j=col in 0-indexed) + X
# Index map: field I = 3i+j for M^{i+1}_{j+1}, field 9 = X
#
# W = sum_j m_j M^j_j + X(det M - Lambda^6)
#
# Fermion mass matrix W_{IJ} = d^2W/dPhi_I dPhi_J evaluated at vacuum.
#
# There are three types of second derivatives:
#
# (a) W_{(ij)(kl)} = X * d^2(det M) / dM^i_j dM^k_l
#     For 3x3: d^2(det)/dM^i_j dM^k_l = epsilon_{ikp} epsilon_{jlq} M^p_q
#     At diagonal vacuum (M^p_q = M_q delta_pq):
#     = sum_p epsilon_{ikp} epsilon_{jlp} M_p  (since q=p required)
#
# (b) W_{(ij),X} = d(det M)/dM^i_j = cofactor(j,i)
#     At diagonal vacuum: (det M / M_j) * delta_{ij} = (M_k M_l) * delta_{ij}
#     where {j,k,l} = {0,1,2}
#
# (c) W_{XX} = 0

print("\n--- Building the 10x10 fermion mass matrix analytically ---")

def levi3(i,j,k):
    """3D Levi-Civita (0-indexed)"""
    if (i,j,k) in [(0,1,2),(1,2,0),(2,0,1)]: return 1
    if (i,j,k) in [(0,2,1),(2,1,0),(1,0,2)]: return -1
    return 0

W2 = np.zeros((10, 10))

# (a) M-M block
for i in range(3):
    for j in range(3):
        I = 3*i + j
        for k in range(3):
            for l in range(3):
                J = 3*k + l
                val = 0.0
                for p in range(3):
                    val += levi3(i,k,p) * levi3(j,l,p) * M[p]
                W2[I,J] = X * val

# (b) M-X cross terms
for i in range(3):
    for j in range(3):
        I = 3*i + j
        # d(det M)/dM^i_j at diagonal vacuum
        # = epsilon_{ikl} epsilon_{jmn} M^k_m M^l_n / 2
        # At diagonal: only m=k, n=l survive
        cof = 0.0
        for k in range(3):
            for l in range(3):
                cof += levi3(i,k,l) * levi3(j,k,l) * M[k] * M[l]
        cof /= 2.0  # factor from the epsilon product
        # Actually: d(det M)/dM^i_j = cofactor_{ji}
        # For diagonal M: cof_{ji} = 0 if i!=j
        # For i=j: cof_{ii} = product of other two diagonal entries
        W2[I, 9] = cof
        W2[9, I] = cof

# (c) X-X = 0 (already)

# Print nonzero entries
print("\nNonzero entries of W_{IJ}:")
field_labels = [f'M{i}{j}' for i in range(3) for j in range(3)] + ['X']
for I in range(10):
    for J in range(I, 10):
        if abs(W2[I,J]) > 1e-15:
            print(f"  W[{field_labels[I]:>3s},{field_labels[J]:>3s}] = {W2[I,J]:>20.8e}")

# Eigenvalues
evals = np.linalg.eigvalsh(W2)
print(f"\nEigenvalues of W_{{IJ}} (sorted by magnitude):")
for i, ev in enumerate(sorted(evals, key=abs)):
    print(f"  |lambda_{i+1}| = {abs(ev):>16.8e} MeV   (sign: {'+'if ev>=0 else '-'})")

# ============================================================
# ANALYTICAL BLOCK STRUCTURE
# ============================================================
print("\n--- Block Structure Analysis ---")
print()
print("The 10 fields decompose into blocks under the unbroken flavor symmetry")
print("(which is nothing, since m_u != m_d != m_s, but still useful to organize):")
print()
print("Diagonal mesons (M_11, M_22, M_33) + X form a 4x4 block.")
print("Off-diagonal pairs (M_ij, M_ji) for i<j form 2x2 blocks.")
print()

# Extract the diagonal + X block (fields 0, 4, 8, 9)
diag_idx = [0, 4, 8, 9]  # M00, M11, M22, X
W2_diag = W2[np.ix_(diag_idx, diag_idx)]

print("4x4 block (M_11, M_22, M_33, X):")
print()
dlabels = ['M_uu', 'M_dd', 'M_ss', 'X']
print(f"{'':>8s}", end='')
for l in dlabels:
    print(f"{l:>16s}", end='')
print()
for i in range(4):
    print(f"{dlabels[i]:>8s}", end='')
    for j in range(4):
        print(f"{W2_diag[i,j]:>16.6e}", end='')
    print()

evals_diag = np.linalg.eigvalsh(W2_diag)
print(f"\nEigenvalues of 4x4 diagonal block:")
for ev in sorted(evals_diag, key=abs):
    print(f"  {ev:>16.8e} MeV")

# Off-diagonal 2x2 blocks: (M_ij, M_ji) for i<j
print()
for i,j in [(0,1), (0,2), (1,2)]:
    I1 = 3*i + j  # M^i_j
    I2 = 3*j + i  # M^j_i
    block = W2[np.ix_([I1,I2], [I1,I2])]
    ev_block = np.linalg.eigvalsh(block)
    fi, fj = labels_3[i], labels_3[j]
    print(f"2x2 block (M_{fi}{fj}, M_{fj}{fi}):")
    print(f"  [{block[0,0]:>14.6e}, {block[0,1]:>14.6e}]")
    print(f"  [{block[1,0]:>14.6e}, {block[1,1]:>14.6e}]")
    print(f"  Eigenvalues: {ev_block[0]:>14.6e}, {ev_block[1]:>14.6e}")
    print()

# ============================================================
# UNDERSTAND THE OFF-DIAGONAL BLOCKS
# ============================================================
print("--- Off-diagonal block analysis ---")
print()
print("For the (M^i_j, M^j_i) block where i != j:")
print("W[M^i_j, M^j_i] = X * epsilon_{ijk} epsilon_{jik} M_k")
print("                  = X * (-1) * M_k  (since epsilon_{ijk}*epsilon_{jik} = -1)")
print("  where k is the remaining index.")
print()
print("The 2x2 block is:")
print("  [  0      , -X*M_k ]")
print("  [ -X*M_k  ,   0    ]")
print()
print("Eigenvalues: +/- |X| * M_k")
print()

for i,j in [(0,1), (0,2), (1,2)]:
    k = 3 - i - j
    fi, fj, fk = labels_3[i], labels_3[j], labels_3[k]
    val = abs(X) * M[k]
    print(f"  ({fi}{fj},{fj}{fi}): mass = |X|*M_{fk} = {abs(X):.6e} * {M[k]:.4f} = {val:.8e} MeV")

# ============================================================
# UNDERSTAND THE DIAGONAL + X BLOCK
# ============================================================
print("\n--- Diagonal + X block analysis ---")
print()
print("W[M_ii, M_jj] = X * M_k  (where k is the third index, i!=j!=k)")
print("W[M_ii, X] = M_j * M_k   (product of the OTHER two diagonal VEVs)")
print("W[X, X] = 0")
print()

# Build it analytically
W2_diag_analytic = np.zeros((4, 4))
for a in range(3):
    for b in range(3):
        if a == b:
            W2_diag_analytic[a, a] = 0  # W[M_ii, M_ii] = 0 (no repeated index in epsilon)
        else:
            k = 3 - a - b
            W2_diag_analytic[a, b] = X * M[k]
    # M-X cross term
    others = [M[k] for k in range(3) if k != a]
    W2_diag_analytic[a, 3] = others[0] * others[1]
    W2_diag_analytic[3, a] = others[0] * others[1]
W2_diag_analytic[3, 3] = 0

print("Analytical 4x4 block:")
for i in range(4):
    print(f"  {dlabels[i]:>8s}", end='')
    for j in range(4):
        print(f"{W2_diag_analytic[i,j]:>16.6e}", end='')
    print()

# Check agreement
print(f"\nMax diff numerical vs analytical: {np.max(np.abs(W2_diag - W2_diag_analytic)):.6e}")

evals_4x4 = np.linalg.eigvalsh(W2_diag_analytic)
print(f"\nEigenvalues of analytical 4x4 block:")
for ev in sorted(evals_4x4, key=abs):
    print(f"  {ev:>16.8e} MeV")

# ============================================================
# PART 3: SCALAR MASS-SQUARED = |W_IJ|^2 (in SUSY vacuum)
# ============================================================
print("\n\n" + "=" * 80)
print("PART 3: SCALAR MASS-SQUARED MATRIX (N_f = N_c = 3)")
print("=" * 80)

print("\nIn the SUSY vacuum, all F-terms vanish: F_I = dW/dPhi_I = 0.")
print("The scalar potential is V = sum_I |F_I|^2 = 0 at the vacuum.")
print()
print("The scalar mass-squared matrix (for complex scalars) is:")
print("  (M^2_s)_{IJ} = sum_K W*_{IK} W_{KJ}  =  (W^dag W)_{IJ}")
print()
print("Since W_{IJ} is real and symmetric in our case:")
print("  M^2_s = W^T W = W^2")
print()
print("The eigenvalues of M^2_s are the squares of the eigenvalues of W_{IJ}.")
print("Thus: m^2_{scalar,n} = (m_{fermion,n})^2 for each n.")
print("The scalar and fermion spectra are EXACTLY degenerate, as required by SUSY.")

scalar_mass_sq = evals**2
print(f"\nScalar mass-squared eigenvalues (from fermion eigenvalues squared):")
for i, (mf, ms2) in enumerate(sorted(zip(np.abs(evals), evals**2), key=lambda x: x[0])):
    print(f"  m_f = {mf:>16.8e} MeV  =>  m^2_s = {ms2:>16.8e} MeV^2")

# ============================================================
# PART 4: SUPERTRACE
# ============================================================
print("\n\n" + "=" * 80)
print("PART 4: SUPERTRACE")
print("=" * 80)

print("""
STr[M^2] = sum_n (-1)^{2J_n} (2J_n + 1) m_n^2

For each chiral multiplet (1 complex scalar + 1 Weyl fermion):
  Contribution = (+1)(m^2_scalar) + (-1)(2)(m^2_fermion)

Wait -- more carefully:
  Complex scalar: 2 real d.o.f., J=0, contributes 2 * m^2_s
  Weyl fermion: 2 helicity states, J=1/2, contributes -2 * m^2_f

  STr[M^2] = sum_n [2*m^2_{scalar,n} - 2*m^2_{fermion,n}]

In SUSY vacuum: m^2_{scalar,n} = m^2_{fermion,n} for all n.
Therefore STr[M^2] = 0 IDENTICALLY.

This is not a check of the computation -- it's guaranteed by the SUSY algebra.
The SUSY vacuum preserves supersymmetry, so the supertrace must vanish.
""")

str_m2 = 2 * np.sum(scalar_mass_sq) - 2 * np.sum(evals**2)
print(f"Explicit check: STr[M^2] = 2*sum(m^2_s) - 2*sum(m^2_f) = {str_m2:.6e} MeV^2")
print(f"(Zero to machine precision, as guaranteed by SUSY.)")

# Higher supertraces
str_m4 = 2 * np.sum(scalar_mass_sq**2) - 2 * np.sum(evals**4)
print(f"\nSTr[M^4] = {str_m4:.6e} MeV^4  (also zero in SUSY vacuum)")

# ============================================================
# PART 5: PHYSICAL INTERPRETATION OF THE SPECTRUM
# ============================================================
print("\n\n" + "=" * 80)
print("PART 5: PHYSICAL INTERPRETATION OF THE SPECTRUM")
print("=" * 80)

print(f"""
The 10 chiral fields at the vacuum of the Seiberg N_f=N_c=3 theory are:
  - 3 diagonal mesons: M_uu, M_dd, M_ss (diagonal M^j_j)
  - 6 off-diagonal mesons: M_ud, M_du, M_us, M_su, M_ds, M_sd
  - 1 Lagrange multiplier: X (enforcing det M = Lambda^6)

The spectrum organizes into:

1. THREE 2x2 off-diagonal blocks (M^i_j, M^j_i):
   Each has eigenvalues +/- |X| * M_k where k is the third flavor.
""")

for i,j in [(0,1), (0,2), (1,2)]:
    k = 3 - i - j
    fi, fj, fk = labels_3[i], labels_3[j], labels_3[k]
    mass = abs(X) * M[k]
    print(f"   ({fi}{fj},{fj}{fi}): m = |X|*M_{fk} = {mass:.8e} MeV = {mass:.4e} MeV")

print(f"""
   These masses are TINY because |X| = {abs(X):.6e} MeV^(-4) is extremely small
   (it scales as m^(1/3) / Lambda^4).

2. ONE 4x4 block (M_uu, M_dd, M_ss, X):
   This is the "physical" sector containing the diagonal mesons and the constraint.
""")

# Let's compute the 4x4 eigenvalues more carefully
print("   4x4 block eigenvalues:")
for ev in sorted(evals_4x4, key=abs):
    print(f"     {ev:>16.8e} MeV")

# Characteristic equation
# The 4x4 matrix has a specific structure. Let me work it out.
# Denote alpha_k = X * M_k and beta_k = prod_{j!=k} M_j
# Then W2_diag has off-diagonal (a,b) entries alpha_k (k = third index)
# and (a,3) entries beta_a.

print(f"\n   Structure parameters:")
for k in range(3):
    alpha_k = X * M[k]
    others = [j for j in range(3) if j != k]
    beta_k = M[others[0]] * M[others[1]]
    print(f"     alpha_{labels_3[k]} = X*M_{labels_3[k]} = {alpha_k:.8e}")
    print(f"     beta_{labels_3[k]} = M_{{{''.join(labels_3[j] for j in others)}}} = {beta_k:.4f}")

# ============================================================
# PART 6: COMPREHENSIVE KOIDE TABLE
# ============================================================
print("\n\n" + "=" * 80)
print("PART 6: COMPREHENSIVE KOIDE TABLE FOR SEESAW SPECTRA")
print("=" * 80)

all_f = ['u', 'd', 's', 'c', 'b']
all_m = [m_u, m_d, m_s, m_c, m_b]

print(f"\n{'Triple':>10s} | {'Q(m)':>8s} | {'Q(1/m)':>8s} | {'Q(1/m)-2/3':>11s} | {'% dev':>7s} | Note")
print("-" * 70)

for combo in combinations(range(5), 3):
    names = ','.join([all_f[i] for i in combo])
    ms = [all_m[i] for i in combo]
    q_m = koide_Q(*ms)
    q_inv = koide_Q(*(1.0/m for m in ms))
    dev = q_inv - 2./3.
    pct = abs(dev) / (2./3.) * 100
    note = ''
    if pct < 2.0:
        note = '<-- near 2/3'
    if abs(q_m - 2./3.) < 0.02:
        note += ' Q(m)~2/3'
    print(f"({names:>8s}) | {q_m:>8.6f} | {q_inv:>8.6f} | {dev:>+11.6f} | {pct:>7.3f} | {note}")

# Special triples from the sBootstrap program
print(f"\n\nSpecial triples from sBootstrap program:")
print(f"  Charged leptons: Q(e,mu,tau) = {koide_Q(0.511, 105.66, 1776.86):.6f}")
print(f"  (-s,c,b) [overlap]: Q = {koide_Q(m_s, m_c, m_b):.6f}")
print(f"  (c,b,t) with m_t=173000: Q = {koide_Q(m_c, m_b, 173000.):.6f}")

# ============================================================
# PART 7: ISS REGIME (N_f = 4)
# ============================================================
print("\n\n" + "=" * 80)
print("PART 7: ISS REGIME (N_f = 4, N_c = 3)")
print("=" * 80)

h = np.sqrt(3.0)
mu_sq = 92.0**2  # MeV^2 (from matching)
m_ISS = np.array([m_u, m_d, m_s, m_c])
flav_ISS = ['u', 'd', 's', 'c']

phi_0 = np.sqrt(mu_sq - m_ISS[0] / h)
print(f"\nphi_0 = {phi_0:.4f} MeV")
print(f"h = sqrt(3) = {h:.6f}")
print(f"mu^2 = {mu_sq:.2f} MeV^2")
print(f"h*mu^2 = {h*mu_sq:.2f} MeV^2")

# F-terms
print(f"\nF-terms:")
for a in range(4):
    if a == 0:
        F_a = 0  # u direction is the aligned direction (B=0 to leading order)
        print(f"  F_{flav_ISS[a]} = 0 (aligned direction)")
    else:
        F_a = h * mu_sq - m_ISS[a]
        print(f"  F_{flav_ISS[a]} = h*mu^2 - m_{flav_ISS[a]} = {F_a:.4f} MeV^2")

print(f"\n--- Full ISS tree-level spectrum ---")
print()

# Sector A: pseudo-modulus
print("SECTOR A (pseudo-modulus sector):")
m_heavy = np.sqrt(2) * h * phi_0
print(f"  Goldstino: m = 0")
print(f"  Heavy fermion: m = sqrt(2)*h*phi_0 = {m_heavy:.4f} MeV")
print(f"  Pseudo-modulus (scalar): m^2 = 0 (flat, lifted by CW)")
print()

# Sector B: link fields
print("SECTOR B (link fields, one per broken flavor):")
print(f"{'Flav':>6s} {'m_q':>8s} {'F_a':>12s} {'m^2_B+':>14s} {'m^2_B-':>14s} {'m_F':>10s} {'STr_B':>12s}")
for a in range(1, 4):
    F_a = h * mu_sq - m_ISS[a]
    m2_plus = h**2 * phi_0**2 + h * F_a
    m2_minus = h**2 * phi_0**2 - h * F_a
    m_F = h * phi_0
    str_B = m2_plus + m2_minus - 2 * m_F**2
    print(f"{flav_ISS[a]:>6s} {m_ISS[a]:>8.2f} {F_a:>12.2f} {m2_plus:>14.2f} {m2_minus:>14.2f} {m_F:>10.2f} {str_B:>12.6f}")

print()

# Sector C: broken mesons
print("SECTOR C (broken mesons chi_a):")
print("  Tree-level: FLAT (massless)")
print("  CW-lifted masses (from iss_cw_koide.md):")
cw_masses = {'u': 27.30, 'd': 27.29, 's': 27.08, 'c': 24.48}
for f in ['u', 'd', 's', 'c']:
    print(f"    m_CW({f}) = {cw_masses[f]:.2f} MeV")

print(f"\n--- ISS Supertrace ---")
print(f"  STr[M^2] = 0 at tree level (exact, by construction)")
print(f"  STr[M^4] != 0 (drives CW potential)")

# ============================================================
# PART 8: SEIBERG SEESAW WITH ALL 5 FLAVORS
# ============================================================
print("\n\n" + "=" * 80)
print("PART 8: SEESAW WITH 5 FLAVORS")
print("=" * 80)

print("""
For N_f = 5, N_c = 3: The Seiberg effective theory is different.
N_f > N_c + 1 = 4, so we are in the free magnetic phase (N_f > 3N_c/2 = 4.5 is FALSE
for N_f=5, N_c=3 since 5 > 4.5). Actually N_f = 5 > 3N_c/2 = 4.5, so:

  N_f = 5 is in the conformal window if 3N_c/2 < N_f < 3N_c, i.e., 4.5 < 5 < 9.
  So N_f = 5 is in the conformal window (IR free magnetic dual).

The magnetic dual has N_c^mag = N_f - N_c = 2 colors.
This is NOT the confining N_f = N_c case, so the quantum-modified constraint
det M - B Bbar = Lambda^6 does not apply directly.

For the seesaw to work with all 5 flavors simultaneously, one would need
to either:
  (a) Take N_c = 5 (not QCD), or
  (b) Treat the 5 flavors in overlapping N_c = 3 blocks

The sBootstrap approach uses option (b): overlapping triples like
(u,d,s), (s,c,b), etc.
""")

# Compute seesaw for all possible triples
print("Seesaw meson VEVs M_j = Lambda^2 * (prod m)^(1/3) / m_j  for all triples:")
print()
print(f"{'Triple':>10s} | {'M_1 (MeV)':>14s} | {'M_2 (MeV)':>14s} | {'M_3 (MeV)':>14s} | {'Q(M)':>8s}")
print("-" * 70)

for combo in combinations(range(5), 3):
    names = [all_f[i] for i in combo]
    ms = np.array([all_m[i] for i in combo])
    C_trip = Lambda**2 * np.prod(ms)**(1./3.)
    M_trip = C_trip / ms
    q = koide_Q(*M_trip)
    label = ','.join(names)
    print(f"({label:>8s}) | {M_trip[0]:>14.2f} | {M_trip[1]:>14.2f} | {M_trip[2]:>14.2f} | {q:>8.6f}")


# ============================================================
# FINAL SUMMARY
# ============================================================
print("\n\n" + "=" * 80)
print("SUMMARY OF KEY RESULTS")
print("=" * 80)

print("""
1. SEIBERG SEESAW (N_f = N_c = 3):
   M_j = Lambda^2 * (prod m_k)^(1/3) / m_j  =>  M_j proportional to 1/m_j

   This maps Koide Q(m_1,m_2,m_3) to Q(1/m_1,1/m_2,1/m_3) exactly
   (scale invariance of the Koide ratio).

2. DUAL KOIDE (inverse mass Koide):
   The closest triple to Q=2/3 in the inverse-mass sector is:
""")

best_combo = None
best_dev = 1.0
for combo in combinations(range(5), 3):
    ms = [all_m[i] for i in combo]
    q = koide_Q(*(1.0/m for m in ms))
    dev = abs(q - 2./3.)
    if dev < best_dev:
        best_dev = dev
        best_combo = combo

names = [all_f[i] for i in best_combo]
ms_best = [all_m[i] for i in best_combo]
q_best = koide_Q(*(1.0/m for m in ms_best))
print(f"   ({','.join(names)}): Q(1/m) = {q_best:.6f}, |Q-2/3| = {best_dev:.6f} ({best_dev/(2./3.)*100:.3f}%)")

print(f"""
   NOTE: The MEMORY file claims Q(1/m_d,1/m_s,1/m_b) = 0.665.
   The correct value is {koide_Q(1/m_d, 1/m_s, 1/m_b):.6f}.
   This appears to be an error in a previous agent's computation (iss_cw_koide.md line 151).

3. FERMION SPECTRUM (confining N_f=N_c=3, (u,d,s)):
   The spectrum has 3 distinct mass scales:
""")

unique_masses = sorted(set([f"{abs(ev):.4e}" for ev in evals]))
for m in sorted(set([abs(ev) for ev in evals])):
    if m > 1e-15:
        count = np.sum(np.isclose(np.abs(evals), m, rtol=1e-3))
        print(f"     m = {m:.6e} MeV  (multiplicity {int(count)})")

heavy_mass = abs(evals[np.argmax(np.abs(evals))])
print(f"""
   The very large mass (~{heavy_mass:.4e} MeV) comes from the
   X-M_diag mixing: the Lagrange multiplier X and diagonal mesons form
   a heavy supermultiplet at the scale Lambda^6 / (M_i * M_j).

   The tiny masses (~1e-4 MeV) come from the off-diagonal meson pairs,
   with masses |X|*M_k.

   The zero eigenvalues correspond to fields with no quadratic coupling.

4. SCALAR SPECTRUM:
   In the SUSY vacuum, m^2(scalar) = m^2(fermion) exactly (by SUSY).
   No SUSY-breaking splittings.

5. SUPERTRACE:
   STr[M^2] = 0 (guaranteed by unbroken SUSY at the confining vacuum).
   STr[M^4] = 0 (same reason).

6. ISS SPECTRUM (N_f=4, metastable):
   SUSY IS broken (rank condition), but STr[M^2] = 0 at tree level
   (canonical Kahler + F-term breaking).

   The CW potential gives m_CW ~ 24-27 MeV for broken mesons (nearly degenerate).
   CW does NOT transmit Koide (pushes Q -> 1/3 = degenerate limit).
   The seesaw DOES transmit Koide (exactly, by scale invariance).
""")
