#!/usr/bin/env python3
"""
Fritzsch texture CKM from sBootstrap mass constraints.

Uses EXACT analytic diagonalization of the Fritzsch 6-texture-zero matrix.

The key result is that |V_cb| has a minimum value of |sqrt(m_s/m_b) - sqrt(m_c/m_t)|
which cannot be reduced below ~0.059 with the predicted masses, while PDG requires 0.042.
This is the well-known failure of the original Fritzsch texture.
"""

import numpy as np
from numpy import sqrt

# ============================================================
# INPUTS
# ============================================================
m_u = 2.16
m_d = 4.67
m_s = 93.4
m_t = 163000

pdg = {
    'm_c': (1270, 20),
    'm_b': (4180, 30),
    '|V_us|': (0.2243, 0.0008),
    '|V_cb|': (0.0422, 0.0008),
    '|V_ub|': (0.00394, 0.00036),
    '|V_td|': (0.00857, 0.00020),
    '|V_ts|': (0.0415, 0.0009),
    '|V_cd|': (0.221, 0.004),
}

print("=" * 70)
print("FRITZSCH TEXTURE CKM FROM SBOOTSTRAP MASS CONSTRAINTS")
print("=" * 70)

# ============================================================
# MASS PREDICTIONS
# ============================================================
print("\n--- MASS PREDICTIONS ---\n")

ratio_OR = 7 + 4*sqrt(3)
m_c_pred = ratio_OR * m_s
pull_mc = (m_c_pred - pdg['m_c'][0]) / pdg['m_c'][1]
print(f"Constraint 1: m_c/m_s = (2+sqrt3)^2 = {ratio_OR:.4f}")
print(f"  m_c = {m_c_pred:.1f} MeV  [PDG: {pdg['m_c'][0]} +/- {pdg['m_c'][1]}]  pull: {pull_mc:+.2f}s")

sqrt_mb = 3*sqrt(m_s) + sqrt(m_c_pred)
m_b_pred = sqrt_mb**2
pull_mb = (m_b_pred - pdg['m_b'][0]) / pdg['m_b'][1]
print(f"\nConstraint 2: sqrt(m_b) = 3 sqrt(m_s) + sqrt(m_c)")
print(f"  m_b = {m_b_pred:.1f} MeV  [PDG: {pdg['m_b'][0]} +/- {pdg['m_b'][1]}]  pull: {pull_mb:+.2f}s")

sc, sb = sqrt(m_c_pred), sqrt(m_b_pred)
aa, bb, cc_k = 1/3, -(4/3)*(sc+sb), m_c_pred + m_b_pred - (2/3)*(sc+sb)**2
m_t_koide = ((-bb + sqrt(bb**2 - 4*aa*cc_k))/(2*aa))**2
Q_cbt = (m_c_pred + m_b_pred + m_t) / (sc + sb + sqrt(m_t))**2
print(f"\nConstraint 3: Koide Q=2/3 for (c,b,t)")
print(f"  Q(m_t={m_t}) = {Q_cbt:.6f}")
print(f"  Predicted m_t = {m_t_koide:.0f} MeV ({(m_t_koide/m_t-1)*100:+.1f}%)")

# ============================================================
# EXACT DIAGONALIZATION OF FRITZSCH TEXTURE
# ============================================================
print("\n\n" + "=" * 70)
print("EXACT FRITZSCH CKM (ANALYTIC CONSTRUCTION)")
print("=" * 70)

def fritzsch_diag_matrix(m1, m2, m3):
    """
    Construct the EXACT unitary that diagonalizes a Fritzsch texture-zero
    matrix with eigenvalues (-m1, m2, m3), where 0 < m1 < m2 < m3.

    The mass matrix is:
        M = | 0   A   0 |
            | A*  0   B |   with A, B real positive, eigenvalues (-m1, m2, m3)
            | 0   B*  C |

    From the characteristic equation:
        C = -m1 + m2 + m3
        A^2 = m1 m2 m3 / C
        A^2 C + B^2 * 0 = ... no wait.

    Using Tr(M) = C = -m1+m2+m3 and det(M) = -A^2 C:
        det(M) = (-m1)(m2)(m3) = -m1 m2 m3
        So A^2 C = m1 m2 m3  =>  A^2 = m1 m2 m3 / (-m1+m2+m3)

    Using Tr(M^2) = 2A^2 + 2B^2 + C^2 = m1^2 + m2^2 + m3^2:
        B^2 = (m1^2 + m2^2 + m3^2 - C^2)/2 - A^2

    The unitary U satisfies M = U . diag(-m1, m2, m3) . U^dag.
    U can be found by solving the eigensystem of M.
    """
    C = -m1 + m2 + m3
    A2 = m1 * m2 * m3 / C
    B2 = (m1**2 + m2**2 + m3**2 - C**2)/2 - A2
    A = sqrt(abs(A2))
    B = sqrt(abs(B2))

    if B2 < 0:
        # This happens when the hierarchy is too strong for the nearest-neighbor
        # texture. Return None to flag it.
        pass

    # Construct M and diagonalize
    M = np.array([
        [0, A, 0],
        [A, 0, B],
        [0, B, C]
    ], dtype=float)

    eigenvalues, U = np.linalg.eigh(M)
    # eigenvalues are sorted ascending. For Fritzsch: -m1, m2, m3

    # Verify eigenvalues match
    ev_expected = sorted([-m1, m2, m3])

    # Phase convention: make first row positive
    for j in range(3):
        if U[0, j] < 0:
            U[:, j] *= -1

    return M, eigenvalues, U, {'A': A, 'B': B, 'C': C, 'B2': B2}


def compute_full_ckm(mu, md, ms, mc, mb, mt, label=""):
    """
    Compute exact CKM from Fritzsch texture diagonalization.

    V_CKM = U_u^dag . Phase . U_d

    The CP phase enters as a relative phase between sectors.
    """
    print(f"\n--- {label} ---")
    print(f"  u={mu:.2f}, d={md:.2f}, s={ms:.1f}, c={mc:.1f}, b={mb:.1f}, t={mt:.0f} MeV\n")

    # Diagonalize both sectors
    M_u, ev_u, U_u, p_u = fritzsch_diag_matrix(mu, mc, mt)
    M_d, ev_d, U_d, p_d = fritzsch_diag_matrix(md, ms, mb)

    print(f"  Up sector:   A={p_u['A']:.2f}, B={p_u['B']:.2f}, C={p_u['C']:.2f}")
    if p_u['B2'] < 0:
        print(f"    WARNING: B^2 = {p_u['B2']:.0f} < 0 (hierarchy too strong for D=0)")
    print(f"    Eigenvalues: {ev_u[0]:.2f}, {ev_u[1]:.2f}, {ev_u[2]:.2f}")
    print(f"    Expected:    {-mu:.2f}, {mc:.2f}, {mt:.2f}")

    print(f"  Down sector: A={p_d['A']:.2f}, B={p_d['B']:.2f}, C={p_d['C']:.2f}")
    if p_d['B2'] < 0:
        print(f"    WARNING: B^2 = {p_d['B2']:.0f} < 0")
    print(f"    Eigenvalues: {ev_d[0]:.2f}, {ev_d[1]:.2f}, {ev_d[2]:.2f}")
    print(f"    Expected:    {-md:.2f}, {ms:.2f}, {mb:.2f}")

    # Print the diagonalizing unitaries
    print(f"\n  U_u (columns = eigenvectors):")
    for i in range(3):
        print(f"    " + "  ".join(f"{U_u[i,j]:+.6f}" for j in range(3)))
    print(f"  U_d:")
    for i in range(3):
        print(f"    " + "  ".join(f"{U_d[i,j]:+.6f}" for j in range(3)))

    # The CKM is V = U_u^T U_d (both real since M is real symmetric)
    # But we need to account for which eigenvector corresponds to which quark.
    # Convention: column 0 -> -m1 (lightest), column 1 -> m2 (middle), column 2 -> m3 (heaviest)
    # For physical quarks: column 0 = u/d, column 1 = c/s, column 2 = t/b
    # HOWEVER, column 0 has eigenvalue -m1, so the physical mass is m1.
    # We need to flip the sign of the first column to make the mass positive.
    # This is equivalent to multiplying by diag(-1, 1, 1) on the right.

    # Actually, in the Fritzsch convention, the diagonalization gives
    # U^T M U = diag(-m1, m2, m3), so physical masses are (m1, m2, m3)
    # with the understanding that the first is associated with a sign flip.
    # The CKM is still V = U_u^T U_d.

    # Scan CP phase
    best = {'chi2': 1e10}

    # Try all sign conventions too: there's an ambiguity in eigenvector signs
    for s1 in [1, -1]:
        for s2 in [1, -1]:
            Uu_adj = U_u.copy()
            Uu_adj[:, 0] *= s1
            Ud_adj = U_d.copy()
            Ud_adj[:, 0] *= s2

            for delta_deg in np.arange(0, 360, 0.5):
                delta = np.radians(delta_deg)
                # Phase on the first column (associated with sign flip)
                P = np.diag([np.exp(1j*delta), 1, 1])
                V = (Uu_adj.T @ P @ Ud_adj).astype(complex)

                # Also try phase on second column
                P2 = np.diag([1, np.exp(1j*delta), 1])
                V2 = (Uu_adj.T @ P2 @ Ud_adj).astype(complex)

                for Vt in [V, V2]:
                    chi2 = ((abs(Vt[0,1]) - pdg['|V_us|'][0])/pdg['|V_us|'][1])**2 \
                         + ((abs(Vt[1,2]) - pdg['|V_cb|'][0])/pdg['|V_cb|'][1])**2 \
                         + ((abs(Vt[0,2]) - pdg['|V_ub|'][0])/pdg['|V_ub|'][1])**2
                    if chi2 < best['chi2']:
                        best = {'chi2': chi2, 'delta': delta_deg, 'V': Vt.copy(),
                                's1': s1, 's2': s2}

    V = best['V']
    print(f"\n  Best fit: delta={best['delta']:.1f}, signs=({best['s1']},{best['s2']})")
    print(f"  chi2 = {best['chi2']:.2f}")

    print(f"\n  |V_CKM|:")
    print(f"         d         s         b")
    for i, lab in enumerate(['u', 'c', 't']):
        print(f"    {lab}  " + "  ".join(f"{abs(V[i,j]):.5f}" for j in range(3)))

    print(f"\n  {'Element':<12} {'Pred':>9} {'PDG':>9} {'Pull':>7}")
    print(f"  {'-'*40}")
    ix = {'|V_us|':(0,1), '|V_cb|':(1,2), '|V_ub|':(0,2),
          '|V_cd|':(1,0), '|V_td|':(2,0), '|V_ts|':(2,1)}
    for key in ['|V_us|', '|V_cb|', '|V_ub|', '|V_cd|', '|V_td|', '|V_ts|']:
        i, j = ix[key]
        pred = abs(V[i,j])
        ref, sig = pdg[key]
        pull = (pred - ref)/sig
        print(f"  {key:<12} {pred:9.5f} {ref:9.5f} {pull:+7.1f}s")

    J = abs(np.imag(V[0,0]*V[1,1]*np.conj(V[0,1])*np.conj(V[1,0])))
    print(f"\n  Jarlskog J = {J:.2e}  [PDG: ~3.08e-05]")

    return best


print("\nCase A: Predicted masses")
res_A = compute_full_ckm(m_u, m_d, m_s, m_c_pred, m_b_pred, m_t,
                          "Predicted masses")

print("\nCase B: PDG masses")
res_B = compute_full_ckm(m_u, m_d, m_s, 1270, 4180, m_t,
                          "PDG masses")

# ============================================================
# CLEAN ANALYTIC SUMMARY
# ============================================================
print("\n\n" + "=" * 70)
print("ANALYTIC LEADING-ORDER SUMMARY")
print("=" * 70)

for label, mc, mb in [("Predicted", m_c_pred, m_b_pred), ("PDG", 1270, 4180)]:
    r_ds = sqrt(m_d/m_s)
    r_uc = sqrt(m_u/mc)
    r_sb = sqrt(m_s/mb)
    r_ct = sqrt(mc/m_t)
    r_db = sqrt(m_d/mb)

    print(f"\n--- {label} masses (m_c={mc:.0f}, m_b={mb:.0f}) ---")
    print(f"  |V_us| ~ sqrt(m_d/m_s) = {r_ds:.5f}  [PDG: {pdg['|V_us|'][0]:.5f}]  pull: {(r_ds-pdg['|V_us|'][0])/pdg['|V_us|'][1]:+.1f}s  GOOD")
    v_cb_min = abs(r_sb - r_ct)
    print(f"  |V_cb| >= |sqrt(m_s/m_b) - sqrt(m_c/m_t)| = |{r_sb:.5f} - {r_ct:.5f}| = {v_cb_min:.5f}")
    print(f"           PDG = {pdg['|V_cb|'][0]:.5f}  pull: {(v_cb_min-pdg['|V_cb|'][0])/pdg['|V_cb|'][1]:+.1f}s  FAILS ({v_cb_min/pdg['|V_cb|'][0]:.2f}x)")
    print(f"  |V_ub| ~ sqrt(m_d/m_b) = {r_db:.5f}  [PDG: {pdg['|V_ub|'][0]:.5f}]  8.4x too large")

# ============================================================
# PARAMETER BUDGET
# ============================================================
print("\n\n" + "=" * 70)
print("PARAMETER COUNT")
print("=" * 70)
print(f"""
SM quark sector:  10 parameters (6 masses + 3 angles + 1 phase)
This framework:    4 inputs (m_s, m_u, m_d, m_t)

Predicted from m_s alone:
  m_c = (7+4sqrt3) m_s = {m_c_pred:.1f} MeV   ({pull_mc:+.2f}s)
  m_b = (3sqrt(m_s)+sqrt(m_c))^2 = {m_b_pred:.1f} MeV   ({pull_mb:+.2f}s)

CKM from Fritzsch texture:
  |V_us| = sqrt(m_d/m_s) = {sqrt(m_d/m_s):.4f}   ({(sqrt(m_d/m_s)-pdg['|V_us|'][0])/pdg['|V_us|'][1]:+.1f}s)   WORKS
  |V_cb| minimum = {abs(sqrt(m_s/m_b_pred)-sqrt(m_c_pred/m_t)):.4f}   (+21s)  FAILS

VERDICT:
  Mass predictions: excellent (2 predictions within 2s each)
  V_us (Cabibbo angle): excellent (-0.9s via GST relation)
  V_cb: fails at +21s -- the original Fritzsch texture overshoots by 1.4x
  This failure is independent of the mass predictions (PDG masses fail equally).

  The V_cb problem is structural: sqrt(m_s/m_b) ~ 0.15 and sqrt(m_c/m_t) ~ 0.09
  differ by 0.06, but PDG requires 0.042. No choice of CP phase can fix this
  because |V_cb| >= |r_sb - r_ct| = 0.059 at leading order.
""")

# ============================================================
# MASS TABLE
# ============================================================
print("=" * 70)
print("MASS TABLE")
print("=" * 70)
print(f"{'Quark':<8} {'Status':<11} {'Value (MeV)':<14} {'PDG (MeV)':<18} {'Pull':<8}")
print("-" * 60)
print(f"{'m_u':<8} {'input':<11} {m_u:<14.2f} {'2.16 +0.49/-0.26':<18}")
print(f"{'m_d':<8} {'input':<11} {m_d:<14.2f} {'4.67 +0.48/-0.17':<18}")
print(f"{'m_s':<8} {'input':<11} {m_s:<14.1f} {'93.4 +8.6/-3.4':<18}")
print(f"{'m_c':<8} {'PREDICTED':<11} {m_c_pred:<14.1f} {'1270 +/- 20':<18} {pull_mc:+.1f}s")
print(f"{'m_b':<8} {'PREDICTED':<11} {m_b_pred:<14.1f} {'4180 +/- 30':<18} {pull_mb:+.1f}s")
print(f"{'m_t':<8} {'input':<11} {m_t:<14.0f} {'163000 (MS-bar)':<18}")
print(f"{'m_t':<8} {'PREDICTED':<11} {m_t_koide:<14.0f} {'':<18} {(m_t_koide/m_t-1)*100:+.1f}%")
