#!/usr/bin/env python3
"""
Local metastable branch of the R-breaking deformation (issue #1).

Language: Python 3 (sympy for exact algebra, mpmath for high-precision
numerics). Derivation: results/rbreaking_metastable_branch.tex.

Model (rigid N=1, canonical Kahler):
    W = f X + m A B + g X A^2 + eps X B,   K = |X|^2 + |A|^2 + |B|^2,
    f, m, g > 0, eps real;  S = m^2 + eps^2,  eta = 2 g f/m^2,  z = 2 g f/S.

Checks:
  A. exact (sympy): stationary point P_eps, F-vector, W_ij F^j = 0, fermion
     singular masses, Goldstino, full 6x6 canonical scalar Hessian, both
     2x2 blocks, determinants, supertrace, Rayleigh quotient, small-eps
     eigenvalues;
  B. one-loop pseudomodulus curvature at eps = 0 with weight -2 per Weyl
     fermion versus the analytic mu_CW^2 = g^2 m^2 B(eta)/(8 pi^2);
  C. full one-loop effective potential (tree + heavy-mode Coleman-Weinberg)
     near P_eps: corrected stationary point, complete Hessian, light
     curvatures for the family eta = 1/2, eps^2 = mu_CW^2/7 at several g,
     and the location of the stability boundary;
  D. rigid Ward identity at the corrected stationary point in the one-loop
     Kahler effective theory (heavy modes integrated out): the covariant
     mass operator annihilates the F-vector, the bare W'' does not.

Run:  python3 results/rbreaking_metastable_branch.py
"""

import sympy as sp
import mpmath as mp

# ======================================================================
# A. Exact tree-level analysis
# ======================================================================
print("=" * 74)
print("A. Exact tree-level analysis at P_eps")
print("=" * 74)
f, m, g = sp.symbols("f m g", positive=True)
e = sp.Symbol("epsilon", real=True, nonzero=True)
S = m**2 + e**2
z = 2 * g * f / S
xr = sp.symbols("xX xA xB", real=True)
yi = sp.symbols("yX yA yB", real=True)
phi = [(xr[k] + sp.I * yi[k]) / sp.sqrt(2) for k in range(3)]
Xs, As, Bs = sp.symbols("X A B")
Wsym = f * Xs + m * As * Bs + g * Xs * As**2 + e * Xs * Bs
flds = [Xs, As, Bs]
W1 = [sp.diff(Wsym, v) for v in flds]
W2 = sp.Matrix(3, 3, lambda i, j: sp.diff(Wsym, flds[i], flds[j]))

B0 = -e * f / S
pt = {Xs: 0, As: 0, Bs: B0}
F = [sp.simplify(w.subs(pt)) for w in W1]
print("  F-vector (W_X, W_A, W_B) =", [sp.factor(x) for x in F])
assert sp.simplify(F[0] - f * m**2 / S) == 0 and sp.simplify(F[1] + f * m * e / S) == 0
assert F[2] == 0
WF = sp.simplify(W2.subs(pt) * sp.Matrix(F))
print("  W_ij F_j =", list(WF), " (stationarity, real parameters)")
assert all(sp.simplify(c) == 0 for c in WF)
V0 = sp.simplify(sum(x**2 for x in F))
print("  V(P_eps) =", sp.factor(V0))
MF = W2.subs(pt)
ev = sp.simplify((MF.T * MF).eigenvals())
print("  fermion mass^2 (singular values^2):", ev)
gold = sp.Matrix([m, -e, 0]) / sp.sqrt(S)
assert sp.simplify(MF * gold) == sp.zeros(3, 1)
print("  Goldstino psi_G = (m psi_X - eps psi_A)/sqrt(S): M_F psi_G = 0")

# canonical real Hessian of V = sum |W_i|^2
Vexpr = sum(sp.expand(w.subs(dict(zip(flds, phi)))) *
            sp.expand(sp.conjugate(w.subs(dict(zip(flds, phi))))) for w in W1)
reals = list(xr) + list(yi)
at = {xr[0]: 0, xr[1]: 0, xr[2]: sp.sqrt(2) * B0, yi[0]: 0, yi[1]: 0, yi[2]: 0}
H = sp.Matrix(6, 6, lambda i, j: sp.simplify(sp.diff(Vexpr, reals[i], reals[j]).subs(at)))
HR = H.extract([0, 1], [0, 1])
HI = H.extract([3, 4], [3, 4])
HR_issue = sp.Matrix([[e**2, e * m * (1 - z)], [e * m * (1 - z), m**2 * (1 + z)]])
HI_issue = sp.Matrix([[e**2, e * m * (1 + z)], [e * m * (1 + z), m**2 * (1 - z)]])
assert sp.simplify(HR - HR_issue) == sp.zeros(2, 2)
assert sp.simplify(HI - HI_issue) == sp.zeros(2, 2)
assert sp.simplify(H[2, 2] - S) == 0 and sp.simplify(H[5, 5] - S) == 0
Hexp = sp.zeros(6, 6)
Hexp[0:2, 0:2] = HR_issue
Hexp[3:5, 3:5] = HI_issue
Hexp[2, 2] = Hexp[5, 5] = S
assert sp.simplify(H - Hexp) == sp.zeros(6, 6)
print("  H_R, H_I equal the issue's blocks; both B components have mass^2 = S;")
print("  all Re-Im and (X,A)-B cross entries vanish")
dR, dI = sp.factor(HR.det()), sp.factor(HI.det())
print("  det H_R =", dR, "   det H_I =", dI)
assert sp.simplify(dR - e**2 * m**2 * z * (3 - z)) == 0
assert sp.simplify(dI + e**2 * m**2 * z * (3 + z)) == 0
print("  STr M^2 = tr H - 2 sum m_F^2 =", sp.simplify(H.trace() - 2 * 2 * S))
assert sp.simplify(H.trace() - 4 * S) == 0
n = sp.Matrix([m, -e]) / sp.sqrt(S)
rq = sp.simplify((n.T * HI * n)[0])
print("  Rayleigh quotient n^T H_I n =", sp.factor(rq))
assert sp.simplify(rq + 6 * g * f * m**2 * e**2 / S**2) == 0
# small-eps light eigenvalues via Schur complement
eta = sp.Symbol("eta", positive=True)
t = sp.Symbol("t", positive=True)  # t = eps/m
subs_eta = {f: eta * m**2 / (2 * g)}
for name, blk in (("R", HR), ("I", HI)):
    schur = sp.simplify((blk[0, 0] - blk[0, 1]**2 / blk[1, 1]).subs(subs_eta).subs(e, t * m))
    print(f"  light m_{name}^2 / eps^2 = ",
          sp.simplify(sp.series(schur / (t**2 * m**2), t, 0, 2).removeO()))


# ======================================================================
# Numerical machinery (mpmath)
# ======================================================================
mp.mp.dps = 40


class Model:
    def __init__(self, f, m, g, eps, lam2=None):
        self.f, self.m, self.g, self.eps = mp.mpf(f), mp.mpf(m), mp.mpf(g), mp.mpf(eps)
        self.lam2 = self.m**2 if lam2 is None else mp.mpf(lam2)

    def W1(self, X, A, B):
        return [self.f + self.g * A**2 + self.eps * B,
                self.m * B + 2 * self.g * X * A,
                self.m * A + self.eps * X]

    def W2(self, X, A, B):
        g, m, e = self.g, self.m, self.eps
        return mp.matrix([[0, 2 * g * A, e], [2 * g * A, 2 * g * X, m], [e, m, 0]])

    def Vtree(self, X, A, B):
        return sum(abs(w)**2 for w in self.W1(X, A, B))

    def scalar_m2(self, X, A, B):
        """Six real-field mass^2 from [[Wbar W, conj(C)], [C, conj(Wbar W)]]."""
        W2 = self.W2(X, A, B)
        W1 = self.W1(X, A, B)
        Am = W2.H * W2
        Am = mp.matrix([[Am[i, j] for j in range(3)] for i in range(3)])
        Cm = mp.matrix(3, 3)
        Cm[0, 1] = Cm[1, 0] = 2 * self.g * mp.conj(W1[1])   # W_XAA conj(W_A)
        Cm[1, 1] = 2 * self.g * mp.conj(W1[0])              # W_AAX conj(W_X)
        M = mp.matrix(6, 6)
        for i in range(3):
            for j in range(3):
                M[i, j] = Am[i, j]
                M[i + 3, j + 3] = mp.conj(Am[i, j])
                M[i, j + 3] = mp.conj(Cm[i, j])
                M[i + 3, j] = Cm[i, j]
        return sorted(mp.eighe(M, eigvals_only=True), key=lambda x: mp.re(x))

    def fermion_m2(self, X, A, B):
        W2 = self.W2(X, A, B)
        return sorted(mp.eighe(W2.H * W2, eigvals_only=True), key=lambda x: mp.re(x))

    def Vcw_heavy(self, X, A, B):
        """(1/64 pi^2)[sum_heavy scalars - 2 sum_heavy Weyl] m^4 (log m^2/Lam^2 - 3/2)."""
        sc = self.scalar_m2(X, A, B)
        sc_heavy = sorted(sc, key=lambda x: abs(x))[2:]
        fe_heavy = self.fermion_m2(X, A, B)[1:]
        h = lambda x: mp.re(x)**2 * (mp.log(mp.re(x) / self.lam2) - mp.mpf(3) / 2)
        return (sum(h(x) for x in sc_heavy) - 2 * sum(h(x) for x in fe_heavy)) / (64 * mp.pi**2)

    def Veff(self, *r):
        X, A, B = [(r[k] + 1j * r[k + 3]) / mp.sqrt(2) for k in range(3)]
        return self.Vtree(X, A, B) + self.Vcw_heavy(X, A, B)


def B_eta(eta):
    eta = mp.mpf(eta)
    return ((1 + eta)**2 * mp.log(1 + eta) - (1 - eta)**2 * mp.log(1 - eta)) / eta - 2


def mu2_cw(g, m, eta):
    return mp.mpf(g)**2 * mp.mpf(m)**2 * B_eta(eta) / (8 * mp.pi**2)


# ======================================================================
# B. Pseudomodulus curvature at eps = 0
# ======================================================================
print()
print("=" * 74)
print("B. One-loop pseudomodulus curvature (weight -2 per Weyl fermion)")
print("=" * 74)
for etav in (mp.mpf(1) / 4, mp.mpf(1) / 2, mp.mpf(3) / 4):
    gv, mv = mp.mpf("0.1"), mp.mpf(1)
    md = Model(etav * mv**2 / (2 * gv), mv, gv, 0)
    Vx = lambda xr_: md.Vcw_heavy(xr_ / mp.sqrt(2), 0, 0)
    curv0 = mp.diff(Vx, 0, 2)
    print(f"  eta = {float(etav):.2f}: d^2V/dX_R^2 at X=0: {mp.nstr(curv0, 12)}; "
          f"analytic mu_CW^2 = {mp.nstr(mu2_cw(gv, mv, etav), 12)}; "
          f"ratio = {mp.nstr(curv0 / mu2_cw(gv, mv, etav), 10)}")
print("  positive series: B(eta) = sum_n 4 eta^(2n)/((2n-1)(2n)(2n+1));",
      "B(1/2) =", mp.nstr(B_eta(0.5), 12), "=",
      mp.nstr(mp.nsum(lambda k: 4 * mp.mpf(0.5)**(2 * k) / ((2 * k - 1) * (2 * k) * (2 * k + 1)),
                      [1, mp.inf]), 12))


# ======================================================================
# C. Corrected stationary branch and light curvatures
# ======================================================================
print()
print("=" * 74)
print("C. Tree + heavy-mode CW near P_eps: family eta = 1/2, eps^2 = mu_CW^2/7")
print("=" * 74)


def branch(gv, ratio, mv=1, etav=mp.mpf(1) / 2, full=True):
    gv, mv = mp.mpf(gv), mp.mpf(mv)
    mu2 = mu2_cw(gv, mv, etav)
    eps = mp.sqrt(ratio * mu2)
    md = Model(etav * mv**2 / (2 * gv), mv, gv, eps)
    S_ = mv**2 + eps**2
    b0 = -eps * md.f / S_
    dVb = lambda xb: mp.diff(lambda u: md.Veff(0, 0, u, 0, 0, 0), xb)
    xb = mp.findroot(dVb, mp.sqrt(2) * b0)
    p = [0, 0, xb, 0, 0, 0]
    grad = [mp.diff(lambda u, k=k: md.Veff(*[u if j == k else p[j] for j in range(6)]), p[k])
            for k in range(6)]
    if not full:
        HI = mp.matrix(2, 2)
        idx = [3, 4]
    else:
        idx = list(range(6))
    Hm = mp.matrix(len(idx), len(idx))
    for a, i in enumerate(idx):
        for b, j in enumerate(idx):
            if b < a:
                Hm[a, b] = Hm[b, a]
                continue
            orders = [0] * 6
            orders[i] += 1
            orders[j] += 1
            Hm[a, b] = mp.diff(md.Veff, tuple(p), tuple(orders))
    evs = mp.eigsy(Hm, eigvals_only=True)
    return md, mu2, xb / mp.sqrt(2), b0, grad, Hm, sorted(evs)


ell = lambda gv: mp.mpf(gv)**2 / (8 * mp.pi**2)
for gv in ("0.3", "0.1", "0.03"):
    md, mu2, Bst, b0, grad, Hm, evs = branch(gv, mp.mpf(1) / 7)
    HR = mp.matrix([[Hm[0, 0], Hm[0, 1]], [Hm[1, 0], Hm[1, 1]]])
    HI = mp.matrix([[Hm[3, 3], Hm[3, 4]], [Hm[4, 3], Hm[4, 4]]])
    mR2 = min(mp.eigsy(HR, eigvals_only=True))
    mI2 = min(mp.eigsy(HI, eigvals_only=True))
    blk = [0, 0, 1, 2, 2, 3]
    off = max(abs(Hm[i, j]) for i in range(6) for j in range(6) if blk[i] != blk[j])
    print(f"  g = {gv:>4} (ell = {mp.nstr(ell(gv), 3)}): eps/m = {mp.nstr(md.eps, 4)}, "
          f"B/B_tree - 1 = {mp.nstr(Bst / b0 - 1, 3)}, max|grad V| = {mp.nstr(max(abs(x) for x in grad), 2)}")
    print(f"      m_I^2/mu^2 = {mp.nstr(mI2 / mu2, 8)} (leading 1/2),  "
          f"m_R^2/mu^2 = {mp.nstr(mR2 / mu2, 8)} (leading 47/42 = {mp.nstr(mp.mpf(47)/42, 8)})")
    print(f"      heavy eigenvalues/m^2 = {[mp.nstr(x, 6) for x in evs[2:]]};  "
          f"max cross-block entry = {mp.nstr(off, 2)}")
    print(f"      deviations/ell: I {mp.nstr((mI2 / mu2 - mp.mpf(1)/2) / ell(gv), 4)}, "
          f"R {mp.nstr((mR2 / mu2 - mp.mpf(47)/42) / ell(gv), 4)}")

# stability boundary: eps^2/mu^2 where m_I^2 = 0 (leading prediction 2/7 at eta = 1/2)
print("  stability boundary r* = eps^2/mu_CW^2 with m_I^2 = 0 (leading (1-eta)/(eta(3+eta)) = 2/7):")
for gv in ("0.1", "0.03"):
    def mI2_of(r, gv=gv):
        md, mu2, Bst, b0, grad, Hm, evs = branch(gv, r, full=False)
        return min(mp.eigsy(Hm, eigvals_only=True)) / mu2
    rstar = mp.findroot(mI2_of, mp.mpf(2) / 7)
    print(f"    g = {gv}: r* = {mp.nstr(rstar, 8)}  (2/7 = {mp.nstr(mp.mpf(2)/7, 8)}), "
          f"(r* - 2/7)/ell = {mp.nstr((rstar - mp.mpf(2)/7) / ell(gv), 4)}")


# ======================================================================
# D. Rigid Ward identity at the corrected stationary point (one loop)
# ======================================================================
print()
print("=" * 74)
print("D. Ward identity at the corrected stationary point (one loop)")
print("=" * 74)
# Goldstino theorem for the 1PI action at zero momentum:
#   (W_ij + Sigma_ij) <Wbar^j> = d_i V_eff,   <W_j> = W_j(phi) + dW_j,
# with dW_j = (1/2) W_jkl <phi^k phi^l> the auxiliary-field tadpole.
# Here W_Bkl = 0, so Sigma_Bj = 0 at one loop, and at X = A = 0 the X and A
# rows vanish by the bosonic Z2 (X, A, B) -> (-X, -A, B). The nontrivial
# B row is  W_Bj (Wbar^j(phi*) + dWbar^j) = 0. The bare W'' alone misses
# W_Bj Wbar^j(phi*) = S dB != 0. The propagator below is the heavy-mode
# part of (1/16 pi^2) M^2 (log M^2/Lam^2 - 1) (MS-bar, the same truncation
# as V_CW), with <phi^k phi^l> = G[k, 3+l] in the (dphi, dphibar) basis.
mp.mp.dps = 40
for gv in ("0.1", "0.03"):
    md, mu2, Bst, b0, grad, Hm, evs = branch(gv, mp.mpf(1) / 7, full=False)
    X, A, B = mp.mpf(0), mp.mpf(0), Bst
    W1 = md.W1(X, A, B)
    W2 = md.W2(X, A, B)
    Wbar = [mp.conj(w) for w in W1]
    tree_rows = [sum(W2[i, j] * Wbar[j] for j in range(3)) for i in range(3)]
    # scalar mass matrix and heavy-mode propagator
    Am = W2.H * W2
    Cm = mp.matrix(3, 3)
    Cm[0, 1] = Cm[1, 0] = 2 * md.g * mp.conj(W1[1])
    Cm[1, 1] = 2 * md.g * mp.conj(W1[0])
    M = mp.matrix(6, 6)
    for i in range(3):
        for j in range(3):
            M[i, j] = Am[i, j]
            M[i + 3, j + 3] = mp.conj(Am[i, j])
            M[i, j + 3] = mp.conj(Cm[i, j])
            M[i + 3, j] = Cm[i, j]
    E, Q = mp.eighe(M)
    order = sorted(range(6), key=lambda k: abs(E[k]))
    G = mp.matrix(6, 6)
    for k in order[2:]:
        lam = mp.re(E[k])
        w = lam * (mp.log(lam / md.lam2) - 1) / (16 * mp.pi**2)
        for a in range(6):
            for b in range(6):
                G[a, b] += w * Q[a, k] * mp.conj(Q[b, k])
    pp = lambda k, l: G[k, 3 + l]          # <phi^k phi^l>
    dW = [md.g * pp(1, 1),                 # (1/2) W_XAA <AA> * 2/2 -> g <AA>
          2 * md.g * pp(0, 1),             # (1/2)(W_AXA + W_AAX) <XA> -> 2g <XA>
          mp.mpf(0)]                       # W_Bkl = 0
    dWbar = [mp.conj(x) for x in dW]
    corr_rows = [sum(W2[i, j] * (Wbar[j] + dWbar[j]) for j in range(3)) for i in range(3)]
    dVB = (mp.diff(lambda u: md.Vcw_heavy(0, 0, (mp.sqrt(2) * B + u) / mp.sqrt(2)), 0)
           * mp.sqrt(2) / 2)               # holomorphic d/dB of V_CW (real direction)
    S_ = md.m**2 + md.eps**2
    print(f"  g = {gv}: B*/B_tree - 1 = {mp.nstr(Bst / b0 - 1, 4)}")
    print(f"     bare rows W_ij Wbar^j(phi*)          = {[mp.nstr(x, 4) for x in tree_rows]}")
    print(f"       (B row equals S dB = {mp.nstr(S_ * (Bst - b0), 4)})")
    print(f"     auxiliary tadpole dW_X, dW_A          = {mp.nstr(dW[0], 4)}, {mp.nstr(dW[1], 4)}")
    print(f"     W_Bj dWbar^j = {mp.nstr(corr_rows[2] - tree_rows[2], 6)},  dV_CW/dB = {mp.nstr(dVB, 6)}")
    print(f"     corrected rows W_ij (Wbar^j + dWbar^j) = {[mp.nstr(x, 3) for x in corr_rows]}")
    print(f"     ratio |corrected B row| / |bare B row| = "
          f"{mp.nstr(abs(corr_rows[2]) / abs(tree_rows[2]), 3)}")
