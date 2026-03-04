"""
One-loop effective Kahler potential in the O'Raifeartaigh model.

Superpotential: W = f X + m phi phi_tilde + g X phi^2
Chiral superfields: X (pseudo-modulus), phi, phi_tilde
Parameters: f, m, g real and positive

Dimensionless variables:
  v = g <X> / m       (pseudo-modulus VEV)
  y = g f / m^2       (SUSY-breaking parameter)

All masses in units of m.  CW potential in units of m^4/(64 pi^2).
"""

import numpy as np
from mpmath import mp, mpf, sqrt as msqrt, log as mlog, diff, taylor, matrix, eig

# Use high precision for reliable Taylor coefficients
mp.dps = 50  # 50 decimal places


# ======================================================================
#  Mass eigenvalues (mpmath for arbitrary precision)
# ======================================================================

def fermion_masses_sq_mp(v):
    """Squared fermion masses (mpmath). Returns [0, mf_-^2, mf_+^2]."""
    v = mpf(v)
    sq = msqrt(v**2 + 1)
    return [mpf(0), (sq - v)**2, (v + sq)**2]


def boson_masses_sq_mp(v, y):
    """Four boson mass^2 from (phi, phi_tilde) sector (mpmath)."""
    v, y = mpf(v), mpf(y)
    # Re sector: eigenvalues of [[4v^2+1+2y, 2v],[2v, 1]]
    tr_Re = 4*v**2 + 2 + 2*y
    det_Re = 1 + 2*y
    disc_Re = tr_Re**2 - 4*det_Re
    lb_Re_p = (tr_Re + msqrt(disc_Re)) / 2
    lb_Re_m = (tr_Re - msqrt(disc_Re)) / 2

    # Im sector: eigenvalues of [[4v^2+1-2y, 2v],[2v, 1]]
    tr_Im = 4*v**2 + 2 - 2*y
    det_Im = 1 - 2*y
    disc_Im = tr_Im**2 - 4*det_Im
    lb_Im_p = (tr_Im + msqrt(disc_Im)) / 2
    lb_Im_m = (tr_Im - msqrt(disc_Im)) / 2

    return sorted([lb_Re_p, lb_Re_m, lb_Im_p, lb_Im_m])


def h_cw_mp(msq):
    """CW integrand m^4(ln|m^2| - 3/2), mpmath.
    For tachyonic modes (m^2 < 0), we use the real part of the CW
    integral: (m^2)^2 (ln|m^2| - 3/2).  Since (m^2)^2 = |m^2|^2
    regardless of sign, this is |m^2|^2 (ln|m^2| - 3/2)."""
    amsq = abs(msq)
    if amsq < mpf('1e-30'):
        return mpf(0)
    return amsq**2 * (mlog(amsq) - mpf(3)/2)


def V_CW_mp(v, y):
    """CW potential in units of m^4/(64 pi^2), mpmath precision."""
    v, y = mpf(v), mpf(y)
    bsq = boson_masses_sq_mp(v, y)
    fsq = fermion_masses_sq_mp(v)
    V = sum(h_cw_mp(b) for b in bsq)
    V -= sum(2 * h_cw_mp(f) for f in fsq)
    return V


# ======================================================================
#  Taylor coefficients via mpmath.taylor()
# ======================================================================

def get_taylor_coefficients(y_val, order=4):
    """
    Use mpmath.taylor to get exact coefficients of V_CW(v) around v=0.

    V_CW is a function of v^2 (proven analytically: all mass eigenvalues
    depend only on v^2).  So V_CW(v) = sum c_{2n} v^{2n}.

    mpmath.taylor gives coefficients of the Taylor series of f(v).
    Since V is even in v, odd coefficients vanish.
    """
    y_val = mpf(y_val)
    f = lambda v: V_CW_mp(v, y_val)
    # Get Taylor coefficients of f(v) around v=0 up to degree 2*order
    # taylor(ctx, f, x, n) returns [f(x), f'(x), f''(x)/2!, ...]
    # i.e., coefficients a_k such that f(x+h) = sum a_k h^k
    coeffs = taylor(f, mpf(0), 2*order)
    # coeffs[k] = (1/k!) * d^k f/dv^k at v=0
    # So the coefficient of v^k is coeffs[k].
    # For even function: coeffs[1] = coeffs[3] = ... = 0
    # c_0 = coeffs[0], c_2 = coeffs[2], c_4 = coeffs[4], c_6 = coeffs[6]
    result = {}
    result['V_0'] = coeffs[0]
    result['c_2'] = coeffs[2]  # coeff of v^2
    result['c_4'] = coeffs[4]  # coeff of v^4
    if 2*order >= 6:
        result['c_6'] = coeffs[6]
    if 2*order >= 8:
        result['c_8'] = coeffs[8]
    return result


# ======================================================================
#  Numpy wrappers for tables
# ======================================================================

def fermion_masses_np(v):
    sq = np.sqrt(v**2 + 1)
    return np.array([0.0, sq - v, v + sq])

def boson_masses_sq_np(v, y):
    M_Re = np.array([[4*v**2 + 1 + 2*y, 2*v], [2*v, 1.0]])
    M_Im = np.array([[4*v**2 + 1 - 2*y, 2*v], [2*v, 1.0]])
    return np.sort(np.concatenate([np.linalg.eigvalsh(M_Re),
                                   np.linalg.eigvalsh(M_Im)]))


# ======================================================================
#  Main computation
# ======================================================================

def main():
    print("=" * 75)
    print("ONE-LOOP EFFECTIVE KAHLER POTENTIAL: O'RAIFEARTAIGH MODEL")
    print("W = f X + m phi phi_tilde + g X phi^2")
    print("=" * 75)

    # ==================================================================
    # PART 1: Mass matrices
    # ==================================================================
    print("\n" + "=" * 75)
    print("PART 1: MASS MATRICES")
    print("=" * 75)

    print("""
Fermion mass matrix M_F (basis: psi_X, psi_phi, psi_phi_tilde):

  M_F / m = [[0,  0,  0],
             [0, 2v,  1],
             [0,  1,  0]]

  v = g<X>/m (dimensionless pseudo-modulus VEV).

  Eigenvalues:  0  (Goldstino),
                m_+ = m (v + sqrt(v^2+1)),
                m_- = m (sqrt(v^2+1) - v).
  Product:  m_+ * m_- = m^2.
  Sum:      m_+ + m_- = 2m sqrt(v^2+1).

Scalar mass-squared matrix (phi, phi_tilde sector; X is flat):

  Re sector:  M_Re^2/m^2 = [[4v^2+1+2y,  2v],  [2v, 1]]
              trace = 4v^2 + 2 + 2y,   det = 1 + 2y

  Im sector:  M_Im^2/m^2 = [[4v^2+1-2y,  2v],  [2v, 1]]
              trace = 4v^2 + 2 - 2y,   det = 1 - 2y

  The +/-2y splitting in the (1,1) entry arises from |F_X|^2 = |f + g phi^2|^2,
  which gives +2gf to Re(phi)^2 and -2gf to Im(phi)^2.

  Tachyon:  Im sector det < 0 when y > 1/2. One scalar becomes tachyonic at
  the origin for gf > m^2/2 (the R-symmetry-breaking vacuum is preferred).

  X sector: m^2 = 0 (tree-level flat direction, lifted at one loop).
""")

    # Numerical table
    print("Explicit mass eigenvalues (units of m):\n")
    hdr = f"  {'v':>6} {'y':>5} | {'m0(G)':>8} {'mf_-':>10} {'mf_+':>10} | {'mb1^2':>9} {'mb2^2':>9} {'mb3^2':>9} {'mb4^2':>9}"
    print(hdr)
    print("  " + "-" * (len(hdr) - 2))

    for y in [0.1, 0.5, 1.0]:
        for v in [0.0, 0.5, 1.0, float(np.sqrt(3))]:
            fm = fermion_masses_np(v)
            bsq = boson_masses_sq_np(v, y)
            tach = " *" if np.any(bsq < -1e-10) else ""
            print(f"  {v:6.3f} {y:5.2f} | {fm[0]:8.4f} {fm[1]:10.6f} {fm[2]:10.6f} |"
                  f" {bsq[0]:9.5f} {bsq[1]:9.5f} {bsq[2]:9.5f} {bsq[3]:9.5f}{tach}")
        print()

    # ==================================================================
    # PART 2: CW potential and Taylor expansion
    # ==================================================================
    print("=" * 75)
    print("PART 2: COLEMAN-WEINBERG POTENTIAL")
    print("=" * 75)

    print("""
CW potential:  V_CW = (1/64pi^2) STr[M^4 (ln M^2/mu^2 - 3/2)]

Counting (per real d.o.f.): +1 for scalars, -1 for fermions.
  4 real scalars from (phi, phi_tilde): weight +1 each.
  2 massive Weyl fermions: weight -2 each (2 d.o.f. per Weyl).
  Goldstino (massless Weyl): contributes 0.
  X scalars (massless): contributes 0.

Supertraces:
  STr[M^2] = 0  (canonical Kahler, verified numerically).
  STr[M^4] = 8y^2  (v-independent).
  Since STr[M^4] is constant, the -3/2 piece does not contribute to
  the v-dependent part of V_CW.
""")

    # Supertrace check
    for y in [0.1, 0.5, 1.0]:
        bsq = boson_masses_sq_np(0, y)
        fsq = np.array([0.0, 1.0, 1.0])  # fermion m^2 at v=0
        str2 = np.sum(bsq) - 2*np.sum(fsq)
        str4 = np.sum(bsq**2) - 2*np.sum(fsq**2)
        print(f"  y={y}: STr[M^2]={str2:.1e}, STr[M^4]={str4:.6f} (expect {8*y**2:.6f})")
    print()

    # V_CW table
    print("V_CW(v) - V_CW(0)  (units of m^4/(64pi^2)):\n")
    for y in [0.1, 0.5, 1.0]:
        print(f"  y = {y}:")
        V0 = float(V_CW_mp(0, y))
        for v in [0.0, 0.2, 0.5, 1.0, 1.5, 2.0]:
            Vv = float(V_CW_mp(v, y))
            print(f"    v={v:4.1f}:  V_CW = {Vv:14.8f},  DeltaV = {Vv - V0:14.10f}")
        print()

    # ==================================================================
    # Taylor coefficients
    # ==================================================================
    print("=" * 75)
    print("TAYLOR EXPANSION: V_CW(v) = V_0 + c_2 v^2 + c_4 v^4 + c_6 v^6 + ...")
    print("=" * 75)
    print()
    print("(Computed via mpmath.taylor at 50-digit precision)")
    print()

    y_values = [0.01, 0.05, 0.1, 0.2, 0.3, 0.5, 0.7, 1.0, 1.5, 2.0, 3.0]
    c2_dict = {}
    c4_dict = {}

    print(f"  {'y':>6} | {'c_2':>18} | {'c_4':>18} | {'c_6':>18}")
    print("  " + "-" * 65)

    for y in y_values:
        tc = get_taylor_coefficients(y, order=4)
        c2 = float(tc['c_2'])
        c4 = float(tc['c_4'])
        c6 = float(tc['c_6'])
        c2_dict[y] = c2
        c4_dict[y] = c4
        print(f"  {y:6.3f} | {c2:+18.10f} | {c4:+18.10f} | {c6:+18.10f}")

    print()

    # Cross-check with direct V evaluation
    print("Cross-check c_2: compare c_2 v^2 with V_CW(v)-V_CW(0) at small v:")
    for y in [0.1, 0.5, 1.0]:
        v_test = 0.01
        V0 = float(V_CW_mp(0, y))
        Vv = float(V_CW_mp(v_test, y))
        dV = Vv - V0
        c2v2 = c2_dict[y] * v_test**2
        c4v4 = c4_dict[y] * v_test**4
        predicted = c2v2 + c4v4
        print(f"  y={y}: DeltaV(v={v_test}) = {dV:.14e}, c_2*v^2+c_4*v^4 = {predicted:.14e}, ratio = {dV/predicted:.10f}")
    print()

    # ==================================================================
    # PART 3: Sign of c_4
    # ==================================================================
    print("=" * 75)
    print("PART 3: SIGN AND MAGNITUDE OF c_4")
    print("=" * 75)
    print()

    print(f"  {'y':>6} | {'c_4':>18} | {'sign':>6} | {'c_4/y^2':>18}")
    print("  " + "-" * 56)
    for y in y_values:
        c4 = c4_dict[y]
        sign = "< 0" if c4 < 0 else "> 0"
        ratio = c4/y**2
        print(f"  {y:6.3f} | {c4:+18.10f} | {sign:>6} | {ratio:+18.8f}")

    print()

    # Find crossover
    neg_y = [y for y in y_values if c4_dict[y] < 0]
    pos_y = [y for y in y_values if c4_dict[y] > 0]
    if neg_y and pos_y:
        # Binary search for crossover
        y_lo = max(y for y in y_values if c4_dict[y] < 0)
        y_hi = min(y for y in y_values if y > y_lo and c4_dict[y] > 0)
        for _ in range(50):
            y_mid = (y_lo + y_hi) / 2
            tc = get_taylor_coefficients(y_mid, order=2)
            c4_mid = float(tc['c_4'])
            if c4_mid < 0:
                y_lo = y_mid
            else:
                y_hi = y_mid
        y_cross = (y_lo + y_hi) / 2
        print(f"  Crossover c_4 = 0 at y = {y_cross:.8f}")
    else:
        y_cross = None
        print("  No sign change detected in the scanned range.")
    print()

    # Value at y = 1 -- seems special
    c4_y1 = c4_dict[1.0]
    print(f"  Notable: at y = 1.0, c_4 = {c4_y1:.12f}")
    print(f"           -20/3 = {-20/3:.12f}")
    print(f"           Difference: {c4_y1 - (-20/3):.4e}")
    print()

    # Is c_4(y=1) exactly -20/3?  Check at higher precision
    mp.dps = 80
    tc_hp = get_taylor_coefficients(1.0, order=3)
    c4_hp = tc_hp['c_4']
    print(f"  High-precision (80 digits): c_4(y=1) = {mp.nstr(c4_hp, 30)}")
    print(f"  -20/3 = {mp.nstr(-mpf(20)/3, 30)}")
    print(f"  Difference: {mp.nstr(c4_hp + mpf(20)/3, 15)}")
    mp.dps = 50
    print()

    # Physical coefficient
    print("Physical quartic coefficient of |X|^4:")
    print()
    print("  V_phys = (m^4/(64pi^2)) * c_4 * (g|X|/m)^4 = [g^4/(64pi^2)] * c_4 * |X|^4")
    print("         = [g^4/(16pi^2)] * (c_4/4) * |X|^4")
    print()
    print(f"  {'y':>6} | {'c_4':>14} | {'C_4 = c_4/4':>14}")
    print("  " + "-" * 38)
    for y in [0.1, 0.5, 1.0]:
        c4 = c4_dict[y]
        print(f"  {y:6.3f} | {c4:+14.8f} | {c4/4:+14.8f}")
    print()

    # In units of g^4/(16pi^2 m^2):
    # The coefficient of |X|^4 is [g^4/(16pi^2)] * C_4 (dimensionless).
    # This has units of [mass^{-4}] * [mass^4] = dimensionless ... no.
    # V has dimension mass^4, |X|^4 has dimension mass^4, so the coefficient
    # of |X|^4 in V is dimensionless.  g^4/(16pi^2) is dimensionless (g is dimensionless
    # in the natural unit convention where m is the unit).
    # So C_4 = c_4/4 IS the answer in units of g^4/(16pi^2).
    #
    # But the question says "in units of g^4/(16pi^2 m^2)". Since V/|X|^4 is
    # dimensionless and g^4/(16pi^2 m^2) has dimension mass^{-2}, there's a
    # dimensional mismatch. The intended normalization is probably:
    # V = [g^4 m^2/(16pi^2)] * alpha * |X/m|^4  where alpha is dimensionless.
    # Then the coeff of |X|^4 is g^4/(16pi^2 m^2) * alpha.
    # alpha = c_4/4.

    print(f"  Answer: the quartic coefficient in units of g^4/(16pi^2 m^2) is")
    print(f"  alpha = c_4/4, where V ⊃ [g^4 m^2/(16pi^2)] * alpha * (|X|/m)^4.")
    print(f"  At y = 1: alpha = {c4_dict[1.0]/4:.8f}")
    print(f"  At y = 0.5: alpha = {c4_dict[0.5]/4:.8f}")
    print(f"  At y = 0.1: alpha = {c4_dict[0.1]/4:.8f}")
    print()

    # ==================================================================
    # PART 4: Comparison with tree-level Kahler c = -g^2/12
    # ==================================================================
    print("=" * 75)
    print("PART 4: COMPARISON WITH NON-CANONICAL KAHLER c = -g^2/12")
    print("=" * 75)

    print("""
Tree-level with K = |X|^2 + c|X|^4/Lambda^2  (c < 0, Lambda = m):

  K_{XX*} = 1 + 4c v^2/g^2 ... but we need to be careful with variables.

  In terms of v = g|X|/m:
  |X|^2 = v^2 m^2/g^2, so |X|^4/m^2 = v^4 m^2/g^4.
  K = v^2 m^2/g^2 + c v^4 m^2/g^4
  K_{XX*} = d^2K/dX dX* = 1 + 4c v^2 m^2/(g^2 m^2)  ... hmm.

  Let me work directly with X.
  K = |X|^2 + c|X|^4/m^2.  (c dimensionless, Lambda = m.)
  K_{XX*} = 1 + 4c|X|^2/m^2 = 1 + 4c v^2/g^2   (since |X| = vm/g).
  V_tree = f^2/K_{XX*} = f^2/(1 + 4cv^2/g^2).

  Pole: 4|c|v^2/g^2 = 1 => v_pole = g/(2 sqrt(|c|)).
  For v_pole = sqrt(3): sqrt(3) = g/(2sqrt(|c|)) => |c| = g^2/12, c = -g^2/12.

  Expansion: V_tree = f^2 * sum_n (-4c/g^2)^n v^{2n}
  = f^2 * sum_n (v^2/3)^n   (for c = -g^2/12)
  Quartic term: f^2 v^4/9 = f^2/9 * v^4.
  In m^4 units: f^2 = y^2 m^4/g^2, so quartic = y^2 m^4/(9g^2) * v^4.

  One-loop quartic: c_4 * v^4 * m^4/(64pi^2).

  Ratio (tree/loop) = y^2 m^4/(9g^2) / [c_4 m^4/(64pi^2)]
                    = 64pi^2 y^2/(9 g^2 c_4).

  For g = m = 1: ratio = 64 pi^2 y^2 / (9 c_4).
""")

    print("Numerical comparison (g/m = 1):")
    print(f"  {'y':>6} | {'c_4(loop)':>14} | {'f^2/9 (tree)':>14} | {'ratio':>12}")
    print("  " + "-" * 52)
    for y in [0.1, 0.5, 1.0]:
        c4 = c4_dict[y]
        # Tree quartic in the SAME units (m^4/(64pi^2)):
        # V_tree quartic = f^2/9 v^4 = y^2/(9*g^2) m^4 v^4.
        # In m^4/(64pi^2) units: divide by m^4/(64pi^2) => 64pi^2 y^2/9  (with g=1)
        tree_q = 64 * np.pi**2 * y**2 / 9
        ratio = c4 / tree_q if abs(tree_q) > 1e-20 else float('nan')
        print(f"  {y:6.3f} | {c4:+14.8f} | {tree_q:14.8f} | {ratio:+12.8f}")

    print()
    print("  The ratio is negative (one-loop quartic has opposite sign to tree")
    print("  quartic from c = -g^2/12) and is y-dependent.  The one-loop effective")
    print("  c_4 is NOT equal to -g^2/12 (or -1/12) in any normalization.")
    print()
    print("  Conclusion: c = -g^2/12 does not arise from one-loop CW computation.")
    print("  It would require a specific relation between y and g that is not")
    print("  naturally present in the model.")
    print()

    # ==================================================================
    # PART 5: Benchmarks at g/m = 0.1, 0.5, 1.0
    # ==================================================================
    print("=" * 75)
    print("PART 5: NUMERICAL BENCHMARKS")
    print("=" * 75)
    print()

    for gm in [0.1, 0.5, 1.0]:
        print(f"--- g/m = {gm} ---\n")

        for y in [0.1, 0.5, 1.0]:
            fm_ratio = y / gm  # f/m
            print(f"  y = {y} (f/m = {fm_ratio:.4f}):")

            # Tree-level mass spectrum at v = 0
            fsq_0 = np.array([0.0, 1.0, 1.0])
            bsq_0 = boson_masses_sq_np(0, y)
            fm_0 = np.sqrt(np.maximum(fsq_0, 0))
            bm_0 = np.sqrt(np.maximum(bsq_0, 0))

            print(f"    Tree-level spectrum at v=0 (units of m):")
            print(f"      Fermions: 0 (Goldstino), {fm_0[1]:.6f}, {fm_0[2]:.6f}")
            print(f"      Scalars (m^2): {bsq_0[0]:.6f}, {bsq_0[1]:.6f}, {bsq_0[2]:.6f}, {bsq_0[3]:.6f}")
            print(f"      Scalars (m):   {bm_0[0]:.6f}, {bm_0[1]:.6f}, {bm_0[2]:.6f}, {bm_0[3]:.6f}")
            if np.any(bsq_0 < -1e-10):
                print("      WARNING: tachyonic scalar at origin")

            c2 = c2_dict[y]
            c4 = c4_dict[y]

            # One-loop mass for X
            # V_phys = [m^4/(64pi^2)] c_2 v^2 + ...
            # = [m^4/(64pi^2)] c_2 (g|X|/m)^2 + ...
            # = [c_2 g^2 m^2/(64pi^2)] |X|^2 + ...
            # m_X^2 = c_2 g^2 m^2/(64pi^2)  [from d^2V/d|X|^2 for complex X]
            # In units of m^2: m_X^2/m^2 = c_2 (g/m)^2 / (64pi^2)
            mx_sq = c2 * gm**2 / (64 * np.pi**2)

            print(f"    One-loop X mass: m_X^2/m^2 = {mx_sq:.6e}")
            print(f"    Coefficients: c_2 = {c2:+.8f}, c_4 = {c4:+.8f}")
            print()

        print()

    # ==================================================================
    # Summary
    # ==================================================================
    print("=" * 75)
    print("SUMMARY OF RESULTS")
    print("=" * 75)
    print(f"""
1. FERMION MASS MATRIX (3x3, basis X, phi, phi_tilde):
   M_F/m = [[0,0,0],[0,2v,1],[0,1,0]]
   Eigenvalues: 0, m(v + sqrt(v^2+1)), m(sqrt(v^2+1) - v)
   Product = m^2, Sum = 2m*sqrt(v^2+1).

2. SCALAR MASS-SQUARED MATRIX (phi, phi_tilde sector, 4 real d.o.f.):
   Re: [[4v^2+1+2y, 2v],[2v, 1]], det = 1+2y, tr = 4v^2+2+2y
   Im: [[4v^2+1-2y, 2v],[2v, 1]], det = 1-2y, tr = 4v^2+2-2y

3. COLEMAN-WEINBERG POTENTIAL (units of m^4/(64pi^2)):
   V_CW(v) = V_0 + c_2 v^2 + c_4 v^4 + ...

   y = 0.1: c_2 = {c2_dict[0.1]:+.8f}, c_4 = {c4_dict[0.1]:+.8f}
   y = 0.5: c_2 = {c2_dict[0.5]:+.8f}, c_4 = {c4_dict[0.5]:+.8f}
   y = 1.0: c_2 = {c2_dict[1.0]:+.8f}, c_4 = {c4_dict[1.0]:+.8f}

4. SIGN OF c_4:
   c_4 < 0 for y < {y_cross:.4f} (confirmed for y = 0.01 through {max(neg_y):.2f}).
   c_4 > 0 for y > {y_cross:.4f}.
   The one-loop correction DOES generate a negative |X|^4 coefficient
   in the regime of moderate SUSY breaking.

   At y = 1: c_4 = {c4_dict[1.0]:.10f} (compare -20/3 = {-20/3:.10f}).

5. COMPARISON WITH c = -g^2/12:
   The one-loop c_4 is NOT equal to -g^2/12 (or any y-independent constant).
   c_4 depends on y = gf/m^2. The tree-level non-canonical Kahler with
   c = -g^2/12 gives a quartic proportional to f^2/9, which is a different
   functional dependence.
""")


if __name__ == "__main__":
    main()
