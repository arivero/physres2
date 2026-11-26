"""
Chiral limit m_u = m_d = 0 of SU(3) SQCD with N_f = N_c = 3 (issue #2).

Language: Python 3 (sympy for exact algebra, numpy for numerics).

Checks, in order:
  1. Rank-one mass theory W = m_s M33 + X(det M - B Bt - L^6):
     an explicit ideal-membership certificate m_s L^6 in (F-terms),
     so no finite common zero of the F-terms exists.
  2. The trajectory M = diag(t, t, L^6/t^2), B = Bt = 0, X = -m_s/t^2:
     every F-term, including F_X, F_B, F_Bt and all off-diagonal F_Mij.
  3. Scalar potential along the trajectory with the weakly Higgsed
     (D-flat squark) metric K = 2 Tr sqrt(M^dag M) on the light block,
     both in the constrained three-flavour description and in the
     two-flavour ADS description W = L2^7/det M_light, L2^7 = m_s L^6.
  4. Massive control: m_u, m_d, m_s all nonzero; the three branches
     C_k^3 = L^6 m_u m_d m_s and their equal-light-mass limit m -> 0.
  5. A stabilizing deformation (squark soft masses) in the two-flavour
     weakly Higgsed theory: stationarity, full scalar Hessian with
     SU(3) D-terms, and the fermion mass matrix with gauginos.

Run:  python3 calculations/sqcd_chiral_limit_runaway.py
"""

import numpy as np
import sympy as sp

np.set_printoptions(linewidth=120)

# ----------------------------------------------------------------------
# 1. No finite F-zero: explicit certificate
# ----------------------------------------------------------------------
print("=" * 72)
print("1. Rank-one mass theory: certificate that no finite F-zero exists")
print("=" * 72)

ms, L = sp.symbols("m_s Lambda", nonzero=True)
X, B, Bt = sp.symbols("X B Bt")
Msym = sp.Matrix(3, 3, lambda i, j: sp.Symbol(f"M{i+1}{j+1}"))
detM = Msym.det()
W = ms * Msym[2, 2] + X * (detM - B * Bt - L**6)

F = {}
for i in range(3):
    for j in range(3):
        F[(i, j)] = sp.diff(W, Msym[i, j])
F_B, F_Bt, F_X = sp.diff(W, B), sp.diff(W, Bt), sp.diff(W, X)

# Row-1 Laplace expansion: sum_j M_1j F_1j = X det M  (no m_s term in row 1)
row1 = sum(Msym[0, j] * F[(0, j)] for j in range(3))
assert sp.expand(row1 - X * detM) == 0
# L^6 X = sum_j M_1j F_1j + B F_B - X F_X
lamX = row1 + B * F_B - X * F_X
assert sp.expand(lamX - L**6 * X) == 0
# m_s L^6 = L^6 F_33 - cof_33 * (L^6 X)
cof33 = Msym[0, 0] * Msym[1, 1] - Msym[0, 1] * Msym[1, 0]
cert = L**6 * F[(2, 2)] - cof33 * lamX
assert sp.expand(cert - ms * L**6) == 0
print("m_s*Lambda^6 = Lambda^6 F_M33 - cof_33 (sum_j M_1j F_M1j + B F_B - X F_X)")
print("  -> identity verified; m_s*Lambda^6 lies in the F-term ideal,")
print("     so the F-terms have no common zero at finite fields.\n")

# ----------------------------------------------------------------------
# 2. All F-terms on the runaway trajectory
# ----------------------------------------------------------------------
print("=" * 72)
print("2. F-terms on M = diag(t,t,L^6/t^2), B = Bt = 0, X = -m_s/t^2")
print("=" * 72)
t = sp.Symbol("t", nonzero=True)
traj = {Msym[i, j]: 0 for i in range(3) for j in range(3) if i != j}
traj.update({Msym[0, 0]: t, Msym[1, 1]: t, Msym[2, 2]: L**6 / t**2,
             B: 0, Bt: 0, X: -ms / t**2})
for (i, j), f in F.items():
    val = sp.simplify(f.subs(traj))
    if val != 0:
        print(f"  F_M{i+1}{j+1} = {val}")
for name, f in (("F_B", F_B), ("F_Bt", F_Bt), ("F_X", F_X)):
    print(f"  {name} = {sp.simplify(f.subs(traj))}")
print("  all other F_Mij = 0\n")

# ----------------------------------------------------------------------
# 3. Potential along the trajectory
# ----------------------------------------------------------------------
print("=" * 72)
print("3. Scalar potential along the runaway")
print("=" * 72)
# (a) constrained description with induced metric: V = |a|^2 - |<a,n>|^2/|n|^2
#     a = dW_tree = m_s dM33, n = d(det M); inverse metric diag(g11, g11, g33)
tp, g33 = sp.symbols("t_abs g33", positive=True)
msa, La = sp.symbols("ms_abs L_abs", positive=True)
g11 = 2 * tp                       # from K = 2 Tr sqrt(M^dag M), K_{t tbar} = 1/(2|t|)
n1 = La**6 / tp                    # |d det / dM11| = |d det / dM22|
n3 = tp**2                         # |d det / dM33|
a2 = g33 * msa**2
an = g33 * msa * n3
nn = 2 * g11 * n1**2 + g33 * n3**2
V_constr = sp.simplify(a2 - an**2 / nn)
V_lead = 4 * msa**2 * La**12 / tp**5
print("  constrained, induced metric:  V =", sp.factor(V_constr))
print("  ratio to 4|m_s|^2|L|^12/|t|^5:", sp.simplify(V_constr / V_lead))
for label, g in (("g33 = 1 (canonical meson)", 1),
                 ("g33 = 2|M33| = 2L^6/t^2 (quark-coordinate estimate)", 2 * La**6 / tp**2)):
    r = sp.simplify((V_constr / V_lead).subs(g33, g))
    print(f"    {label}: ratio = {r}  -> {sp.limit(r, tp, sp.oo)} as |t| -> oo")

# (b) two-flavour ADS description in D-flat squark coordinates Q = Qt^T = v diag(1,1)
v, L7 = sp.symbols("v L7", positive=True)  # L7 = |Lambda_2^7| = |m_s| |Lambda_3|^6
W_ads_grad = L7 / v**5                     # |dW/dQ| for each of the 4 nonzero squarks
V_ads = 4 * W_ads_grad**2
print("  two-flavour ADS, squark coordinates: V =", V_ads,
      "= 4|L2^7|^2/|t|^5 with |t| = v^2")
print("  -> both descriptions give V = 4|m_s|^2 |L3|^12/|t|^5 -> 0 (runaway)\n")

# (c) holomorphic matching: integrate out M33 and X on B = 0
M11, M22, M12, M21 = sp.symbols("M11 M22 M12 M21")
detl = M11 * M22 - M12 * M21
M33_sol = L**6 / detl                        # constraint with M_{3a} = M_{a3} = 0
print("  W_eff = m_s M33 |_(det M = L^6) =", sp.simplify(ms * M33_sol),
      " (= Lambda_2^7/det M_light)\n")

# ----------------------------------------------------------------------
# 4. Massive control: three branches and the equal-light-mass limit
# ----------------------------------------------------------------------
print("=" * 72)
print("4. All masses nonzero: three SUSY branches, m_u = m_d = m -> 0")
print("=" * 72)
Lam = 300.0      # MeV
m_s_val = 93.4   # MeV
omega = np.exp(2j * np.pi / 3)


def branch(mu, md, msv, k):
    C = (Lam**6 * mu * md * msv) ** (1 / 3) * omega**k
    Mv = np.diag([C / mu, C / md, C / msv])
    Xv = -C / Lam**6
    return C, Mv, Xv


def fterm_residual(mvec, Mv, Xv):
    cof = np.linalg.det(Mv) * np.linalg.inv(Mv).T
    FM = np.diag(mvec) + Xv * cof
    FX = np.linalg.det(Mv) - Lam**6
    return np.max(np.abs(FM)) / np.max(np.abs(mvec)), abs(FX) / Lam**6


print(f"  Lambda_3 = {Lam} MeV, m_s = {m_s_val} MeV")
print("  m (MeV)    k   |M_uu|=|M_dd| (MeV^2)   |M_ss| (MeV^2)   |X| (MeV^-4)   "
      "M_ss*M_uu^2/L^6   resid")
masses = [4.67, 1.0, 1e-1, 1e-2, 1e-3, 1e-4, 1e-6]
rows = []
for m in masses:
    for k in range(3):
        C, Mv, Xv = branch(m, m, m_s_val, k)
        rF, rX = fterm_residual(np.array([m, m, m_s_val]), Mv, Xv)
        prod = Mv[2, 2] * Mv[0, 0] ** 2 / Lam**6
        print(f"  {m:9.2e}  {k}   {abs(Mv[0,0]):14.6e}        {abs(Mv[2,2]):12.6e}   "
              f"{abs(Xv):11.4e}   {prod.real:+.6f}{prod.imag:+.1e}i   {max(rF, rX):.1e}")
        if k == 0:
            rows.append((m, abs(Mv[0, 0]), abs(Mv[2, 2]), abs(Xv)))
rows = np.array(rows)
sl_M = np.polyfit(np.log(rows[:, 0]), np.log(rows[:, 1]), 1)[0]
sl_S = np.polyfit(np.log(rows[:, 0]), np.log(rows[:, 2]), 1)[0]
sl_X = np.polyfit(np.log(rows[:, 0]), np.log(rows[:, 3]), 1)[0]
print(f"  fitted exponents in m: M_uu = M_dd ~ m^{sl_M:.4f}, "
      f"M_ss ~ m^{sl_S:.4f}, X ~ m^{sl_X:.4f}  (exact: -1/3, 2/3, 2/3)")
# Each branch lies on the runaway trajectory with t_k = Lam^2 (m_s/m)^{1/3} omega^k
m = 1e-3
for k in range(3):
    C, Mv, Xv = branch(m, m, m_s_val, k)
    tk = Mv[0, 0]
    print(f"  k={k}: t_k = {tk:.6e};  M_ss - L^6/t_k^2 = {abs(Mv[2,2]-Lam**6/tk**2):.1e};"
          f"  X + m_s/t_k^2 = {abs(Xv + m_s_val/tk**2):.1e}")
# Two-flavour ADS with light mass m: F = m - L2^7/(t^3) = 0
for k in range(3):
    tk = (m_s_val * Lam**6 / m) ** (1 / 3) * omega**k
    print(f"  ADS check k={k}: |m - m_s L^6/t^3| = {abs(m - m_s_val*Lam**6/tk**3):.1e}")
print()

# ----------------------------------------------------------------------
# 5. Soft-mass stabilization in the weakly Higgsed two-flavour theory
# ----------------------------------------------------------------------
print("=" * 72)
print("5. Squark soft mass mt^2 (|Q|^2+|Qt|^2): stationarity and fluctuations")
print("=" * 72)
# Units: Lambda_3 = 1, m_s = 1, so Lambda_2^7 = 1.  mt << m_s for control.
gcoup = 1.0
mt = 0.01
L7v = 1.0
v0 = (5 * L7v**2 / mt**2) ** (1 / 12)
print(f"  units Lambda_3 = m_s = 1, g = {gcoup}, mt = {mt}")
print(f"  predicted vev v0 = (5 |L2^7|^2/mt^2)^(1/12) = {v0:.6f}  (v0/Lambda_3 > 1 required)")

# fields: Q (3 colour x 2 flavour), Qt (2 flavour x 3 colour); 12 complex scalars
q = sp.symbols("q0:12")
Qs = sp.Matrix(3, 2, q[:6])
Qts = sp.Matrix(2, 3, q[6:])
Ws = L7v / (Qts * Qs).det()
grad = [sp.diff(Ws, x) for x in q]
hess = [[sp.diff(gi, x) for x in q] for gi in grad]
third = [[[sp.diff(hess[i][j], x) for x in q] for j in range(12)] for i in range(12)]
f_grad = sp.lambdify(q, grad, "numpy")
f_hess = sp.lambdify(q, hess, "numpy")
f_third = sp.lambdify(q, third, "numpy")

# SU(3) generators T^a = lambda^a/2 acting on the 12-component field vector
lam = np.zeros((8, 3, 3), dtype=complex)
lam[0][0, 1] = lam[0][1, 0] = 1
lam[1][0, 1], lam[1][1, 0] = -1j, 1j
lam[2][0, 0], lam[2][1, 1] = 1, -1
lam[3][0, 2] = lam[3][2, 0] = 1
lam[4][0, 2], lam[4][2, 0] = -1j, 1j
lam[5][1, 2] = lam[5][2, 1] = 1
lam[6][1, 2], lam[6][2, 1] = -1j, 1j
lam[7] = np.diag([1, 1, -2]) / np.sqrt(3)
T = lam / 2


def gen12(Ta):
    """Generator on the 12-vector (Q[c,f] row-major, then Qt[f,c] row-major)."""
    G = np.zeros((12, 12), dtype=complex)
    for f in range(2):
        for c in range(3):
            for c2 in range(3):
                G[c * 2 + f, c2 * 2 + f] += Ta[c, c2]                 # Q -> T Q
                G[6 + f * 3 + c, 6 + f * 3 + c2] += -Ta[c2, c]        # Qt -> -Qt T
    return G


G = np.array([gen12(T[a]) for a in range(8)])


def potential(phi):
    Wi = np.array(f_grad(*phi), dtype=complex)
    D = np.array([np.real(np.conj(phi) @ G[a] @ phi) for a in range(8)])
    return np.sum(np.abs(Wi) ** 2) + 0.5 * gcoup**2 * np.sum(D**2) + mt**2 * np.sum(np.abs(phi) ** 2)


phi0 = np.zeros(12, dtype=complex)
phi0[0 * 2 + 0] = v0      # Q[colour 1, flavour 1]
phi0[1 * 2 + 1] = v0      # Q[colour 2, flavour 2]
phi0[6 + 0 * 3 + 0] = v0  # Qt[flavour 1, colour 1]
phi0[6 + 1 * 3 + 1] = v0  # Qt[flavour 2, colour 2]

Wi = np.array(f_grad(*phi0), dtype=complex)
Wij = np.array(f_hess(*phi0), dtype=complex)
Wijk = np.array(f_third(*phi0), dtype=complex)
D = np.array([np.real(np.conj(phi0) @ G[a] @ phi0) for a in range(8)])
dV = Wij @ np.conj(Wi) + mt**2 * np.conj(phi0) \
    + gcoup**2 * sum(D[a] * (np.conj(phi0) @ G[a]) for a in range(8))
print(f"  max |D^a| at v0 = {np.max(np.abs(D)):.1e} (D-flat)")
print(f"  max |dV/dphi| at v0 = {np.max(np.abs(dV)):.1e}  (stationary)")

# scalar mass matrix in (dphi, dphi^*) basis
# V2 = dphi^dag A dphi + 1/2 (dphi^T Bm dphi + c.c.)
A = np.conj(Wij) @ Wij               # A[j,i] = sum_k Wbar_jk W_ki
Bm = np.einsum("ijk,k->ij", Wijk, np.conj(Wi))
for a in range(8):
    Gphi = G[a] @ phi0               # (T^a phi)_i
    phiG = np.conj(phi0) @ G[a]      # (phibar T^a)_i
    A = A + gcoup**2 * (np.outer(Gphi, phiG) + D[a] * G[a])
    Bm = Bm + gcoup**2 * np.outer(phiG, phiG)
A = A + mt**2 * np.eye(12)
Mfull = np.block([[A, np.conj(Bm)], [Bm, np.conj(A)]])
assert np.allclose(Mfull, Mfull.conj().T, atol=1e-12)
ev = np.sort(np.linalg.eigvalsh(Mfull))
tol = 1e-9 * np.max(np.abs(ev))
nzero = int(np.sum(np.abs(ev) < tol))
print(f"  scalar mass^2 eigenvalues (24 real modes), units of mt^2:")
print("   ", np.array2string(ev / mt**2, precision=4, suppress_small=True))
print(f"  zero modes: {nzero} (expected 8 eaten + 3 pions + 1 R-axion = 12);"
      f" negative: {int(np.sum(ev < -tol))}")
# finite-difference cross-check of the Hessian along random real directions
rng = np.random.default_rng(1)
h = 1e-4
for _ in range(3):
    dvec = rng.normal(size=12) + 1j * rng.normal(size=12)
    dvec /= np.linalg.norm(dvec)
    fd = (potential(phi0 + h * dvec) - 2 * potential(phi0) + potential(phi0 - h * dvec)) / h**2
    an = 2 * np.real(np.conj(dvec) @ A @ dvec) + 2 * np.real(dvec @ Bm @ dvec)
    print(f"  Hessian check along random direction: finite-diff {fd:.6e}, analytic {an:.6e}")

# fermion mass matrix: 12 chiral fermions + 8 gauginos
MF = np.zeros((20, 20), dtype=complex)
MF[:12, :12] = Wij
for a in range(8):
    col = np.sqrt(2) * gcoup * (np.conj(phi0) @ G[a])
    MF[:12, 12 + a] = col
    MF[12 + a, :12] = col
sv = np.sort(np.linalg.svd(MF, compute_uv=False))
print("  fermion masses (singular values, units of mt):")
print("   ", np.array2string(sv / mt, precision=4, suppress_small=True))
print(f"  smallest fermion mass / mt = {sv[0]/mt:.4f}: no massless fermion (no Goldstino)")
print("  radial D-flat mode predicted analytically: m^2 = 12 mt^2")

# Physical-scale estimate of the soft minimum
for mt_phys in (92.0, 9.2, 0.92):
    r = (5 * m_s_val**2 / mt_phys**2) ** (1 / 12)
    print(f"  m_s = {m_s_val} MeV, mt = {mt_phys:6.2f} MeV -> v0/Lambda_3 = {r:.3f}")
