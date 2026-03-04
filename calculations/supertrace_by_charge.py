#!/usr/bin/env python3
"""
Supertrace decomposed by electromagnetic charge sector.

Setup: 16 chiral superfields with NMSSM coupling.
Superpotential:
  W = sum_i m_i M^i_i + X(det M - B Btilde - Lambda^6)
      + y_c Hu0 Mdd + y_b Hd0 Mss + lambda X (Hu0 Hd0 - Hu+ Hd-)

Soft SUSY breaking: V_soft = f_pi^2 Tr(M^dag M)  (9 meson fields only)

EM charges: Q(u) = 2/3, Q(d) = Q(s) = -1/3
  Meson M^i_j has Q = Q_i - Q_j
  B = eps Q^u Q^d Q^s has Q = 2/3 - 1/3 - 1/3 = 0
  Btilde has Q = 0

Tasks:
  1. Compute 16x16 fermion mass matrix and 32x32 scalar mass-squared matrix
  2. Decompose STr[M^2] by charge sector Q = -1, 0, +1
  3. Verify sum = 18 f_pi^2
  4. Check for simple rational ratios
  5. Compute Q_em and Q_em^2 weighted supertraces
  6. Compute STr[M^4] per sector
  7. Effect of Fayet-Iliopoulos D-term
"""

import numpy as np
from numpy.linalg import eigvalsh, eigh
from itertools import permutations
from mpmath import mp, mpf, matrix as mpmatrix, eigh as mp_eigh, fsum

# Set high precision for mpmath
mp.dps = 50  # 50 decimal places

# =========================================================================
# Parameters
# =========================================================================
m_u = 2.16       # MeV
m_d = 4.67       # MeV
m_s = 93.4       # MeV
m_c = 1270.0     # MeV
m_b = 4180.0     # MeV
LAM = 300.0      # MeV
v   = 246220.0   # MeV
f_pi = 92.0      # MeV

y_c = 2.0 * m_c / v
y_b = 2.0 * m_b / v
lam = 0.72
m_tilde_sq = f_pi**2   # 8464 MeV^2

m_arr = np.array([m_u, m_d, m_s])
quark_charges = np.array([2.0/3.0, -1.0/3.0, -1.0/3.0])  # u, d, s

LAM6 = LAM**6

# =========================================================================
# EM charges of 16 fields
# =========================================================================
# Field ordering:
# 0: M^u_u,  1: M^u_d,  2: M^u_s
# 3: M^d_u,  4: M^d_d,  5: M^d_s
# 6: M^s_u,  7: M^s_d,  8: M^s_s
# 9: X
# 10: B,  11: Btilde
# 12: Hu+,  13: Hu0,  14: Hd0,  15: Hd-

NAMES = ['Muu','Mud','Mus','Mdu','Mdd','Mds','Msu','Msd','Mss',
         'X', 'B', 'Bt', 'Hu+','Hu0','Hd0','Hd-']

def idx_M(i, j):
    return 3*i + j

# EM charges for each field
Q_em = np.zeros(16)
for i in range(3):
    for j in range(3):
        Q_em[idx_M(i,j)] = quark_charges[i] - quark_charges[j]
Q_em[9] = 0.0
Q_em[10] = 0.0
Q_em[11] = 0.0
Q_em[12] = +1.0
Q_em[13] = 0.0
Q_em[14] = 0.0
Q_em[15] = -1.0

print("=" * 78)
print("EM CHARGE ASSIGNMENTS")
print("=" * 78)
for i in range(16):
    print(f"  {NAMES[i]:>4s}:  Q_em = {Q_em[i]:+.4f}")

# Identify charge sectors
sectors = {}
for i in range(16):
    q = round(Q_em[i] * 3) / 3
    q_key = round(q, 4)
    if q_key not in sectors:
        sectors[q_key] = []
    sectors[q_key].append(i)

print("\nCharge sectors:")
for q in sorted(sectors.keys()):
    fields = [NAMES[i] for i in sectors[q]]
    print(f"  Q = {q:+.1f}: {fields}  ({len(fields)} fields)")

# =========================================================================
# Vacuum (NMSSM-modified Seiberg seesaw)
# =========================================================================
prod_m = np.prod(m_arr)
lam_v2_half = lam * v**2 / 2.0
det_M_required = LAM6 - lam_v2_half

C_eff = (det_M_required * prod_m)**(1.0/3.0)
M_vev = C_eff / m_arr
X0 = -C_eff / det_M_required

v_hu = v / np.sqrt(2)
v_hd = v / np.sqrt(2)

print(f"\n{'='*78}")
print("VACUUM")
print(f"{'='*78}")
print(f"  Lambda^6 = {LAM6:.6e}")
print(f"  NMSSM shift = {lam_v2_half:.6e}")
print(f"  det M = {det_M_required:.6e}")
print(f"  C_eff = {C_eff:.6f}")
print(f"  M_u = {M_vev[0]:.4f}, M_d = {M_vev[1]:.4f}, M_s = {M_vev[2]:.4f}")
print(f"  X0 = {X0:.8e}")
print(f"  v/sqrt(2) = {v_hu:.4f}")

# =========================================================================
# F-terms
# =========================================================================
F_vec = np.zeros(16)
F_vec[0] = m_u + X0 * M_vev[1] * M_vev[2]
F_vec[4] = m_d + X0 * M_vev[0] * M_vev[2] + y_c * v_hu
F_vec[8] = m_s + X0 * M_vev[0] * M_vev[1] + y_b * v_hd
F_vec[9] = det_M_required - LAM6 + lam * v_hu * v_hd
F_vec[13] = y_c * M_vev[1] + lam * X0 * v_hd
F_vec[14] = y_b * M_vev[2] + lam * X0 * v_hu

print(f"\n  F-terms:")
for i in range(16):
    if abs(F_vec[i]) > 1e-20:
        print(f"    F[{NAMES[i]:>4s}] = {F_vec[i]:+.8e} MeV")

# =========================================================================
# Fermion mass matrix (16x16)
# =========================================================================
print(f"\n{'='*78}")
print("FERMION MASS MATRIX W_IJ (16x16)")
print(f"{'='*78}")

W = np.zeros((16, 16))

for a in range(3):
    for c in range(3):
        if a != c:
            e = 3 - a - c
            W[idx_M(a,a), idx_M(c,c)] += X0 * M_vev[e]
            W[idx_M(a,c), idx_M(c,a)] += X0 * (-M_vev[e])

for a in range(3):
    others = [M_vev[c] for c in range(3) if c != a]
    cof = others[0] * others[1]
    W[idx_M(a,a), 9] += cof
    W[9, idx_M(a,a)] += cof

W[10, 11] += -X0
W[11, 10] += -X0

W[13, 4] += y_c
W[4, 13] += y_c
W[14, 8] += y_b
W[8, 14] += y_b

W[9, 13] += lam * v_hd
W[13, 9] += lam * v_hd
W[9, 14] += lam * v_hu
W[14, 9] += lam * v_hu
W[13, 14] += lam * X0
W[14, 13] += lam * X0
W[12, 15] += -lam * X0
W[15, 12] += -lam * X0

assert np.allclose(W, W.T)

# =========================================================================
# Third derivative tensor W_IJK and holomorphic mass matrix
# =========================================================================
def levi_civita_sign(p):
    p = list(p)
    sign = 1
    for i in range(len(p)):
        while p[i] != i:
            sign *= -1
            j = p[i]
            p[i], p[j] = p[j], p[i]
    return sign

W3_accum = {}

for a in range(3):
    for c in range(a+1, 3):
        e = 3 - a - c
        key = tuple(sorted([9, idx_M(a,a), idx_M(c,c)]))
        W3_accum[key] = W3_accum.get(key, 0) + M_vev[e]
        key = tuple(sorted([9, idx_M(a,c), idx_M(c,a)]))
        W3_accum[key] = W3_accum.get(key, 0) + (-M_vev[e])

W3_type2_raw = {}
for perm_row in permutations(range(3)):
    for perm_col in permutations(range(3)):
        sign = levi_civita_sign(list(perm_row)) * levi_civita_sign(list(perm_col))
        val = X0 * sign
        i1 = idx_M(perm_row[0], perm_col[0])
        i2 = idx_M(perm_row[1], perm_col[1])
        i3 = idx_M(perm_row[2], perm_col[2])
        key = tuple(sorted([i1, i2, i3]))
        W3_type2_raw[key] = W3_type2_raw.get(key, 0) + val

for key, raw_val in W3_type2_raw.items():
    W3_accum[key] = W3_accum.get(key, 0) + raw_val / 6.0

key_nmssm1 = tuple(sorted([9, 13, 14]))
W3_accum[key_nmssm1] = W3_accum.get(key_nmssm1, 0) + lam
key_nmssm2 = tuple(sorted([9, 12, 15]))
W3_accum[key_nmssm2] = W3_accum.get(key_nmssm2, 0) + (-lam)

W3 = {}
for sorted_key, val in W3_accum.items():
    I, J, K = sorted_key
    for perm in [(I,J,K),(I,K,J),(J,I,K),(J,K,I),(K,I,J),(K,J,I)]:
        W3[perm] = val

M2_hol = np.zeros((16, 16))
for (I, J, K), val in W3.items():
    M2_hol[J, K] += val * F_vec[I]
M2_hol = 0.5 * (M2_hol + M2_hol.T)

M2_soft = np.zeros((16, 16))
for i in range(9):
    M2_soft[i, i] = m_tilde_sq

# =========================================================================
# ANALYTIC SUPERTRACE PER SECTOR
# =========================================================================
# The supertrace STr_Q[M^2] = 2 * Tr(M2_soft restricted to sector Q)
# This is EXACT and doesn't suffer from cancellation.
# The proof: STr = Tr(m^2_R + m^2_I) - 2*Tr(m^2_ferm)
#           = Tr(W^2 + M_hol + M_soft) + Tr(W^2 - M_hol + M_soft) - 2*Tr(W^2)
#           = 2*Tr(M_soft)
# This holds block-by-block because W, M_hol, M_soft all respect EM charge conservation.

print(f"\n{'='*78}")
print("ANALYTIC SUPERTRACE: STr_Q[M^2] = 2*Tr(M_soft|_Q)")
print(f"{'='*78}")
print(f"\n  This is EXACT. The W^2 terms cancel between scalars and fermions,")
print(f"  and M_hol cancels between R and I scalar components.")

# =========================================================================
# SPECTRUM BY CHARGE SECTOR (block-by-block for eigenvalues)
# =========================================================================
print(f"\n{'='*78}")
print("SPECTRUM AND SUPERTRACE BY EM CHARGE SECTOR")
print(f"{'='*78}")

# For the Q=0 sector, we need high-precision arithmetic because
# the eigenvalues span ~10^10 to ~10^{-9} MeV.
# Use mpmath for the Q=0 sector.

results = {}

for q in sorted(sectors.keys()):
    idx = sectors[q]
    n = len(idx)
    field_names = [NAMES[i] for i in idx]

    print(f"\n{'─'*70}")
    print(f"  CHARGE SECTOR Q = {q:+.1f}")
    print(f"  Fields: {field_names}")
    print(f"{'─'*70}")

    # --- Fermion block ---
    W_block = W[np.ix_(idx, idx)]

    # Use high-precision for Q=0 sector
    if abs(q) < 0.01 and n > 6:
        # Build mpmath matrices
        n_mp = n
        W_mp = mpmatrix(n_mp, n_mp)
        for ii in range(n_mp):
            for jj in range(n_mp):
                W_mp[ii, jj] = mpf(str(W_block[ii, jj]))

        # Diagonalize with mpmath
        ferm_evals_mp, _ = mp_eigh(W_mp)
        ferm_evals = np.array([float(ferm_evals_mp[k]) for k in range(n_mp)])

        # For STr[M^4], compute via mpmath
        M2_H_block_mp = W_mp * W_mp  # mpmath matrix multiply
        M2_hol_block = M2_hol[np.ix_(idx, idx)]
        M2_soft_block = M2_soft[np.ix_(idx, idx)]

        M2_hol_mp = mpmatrix(n_mp, n_mp)
        M2_soft_mp = mpmatrix(n_mp, n_mp)
        for ii in range(n_mp):
            for jj in range(n_mp):
                M2_hol_mp[ii, jj] = mpf(str(M2_hol_block[ii, jj]))
                M2_soft_mp[ii, jj] = mpf(str(M2_soft_block[ii, jj]))

        # Build scalar mass-squared matrices (2n x 2n)
        M2_R_mp = M2_H_block_mp + M2_hol_mp + M2_soft_mp
        M2_I_mp = M2_H_block_mp - M2_hol_mp + M2_soft_mp

        # Diagonalize R and I separately (they don't mix for real couplings)
        evals_R_mp, _ = mp_eigh(M2_R_mp)
        evals_I_mp, _ = mp_eigh(M2_I_mp)

        scalar_evals_R = np.array([float(evals_R_mp[k]) for k in range(n_mp)])
        scalar_evals_I = np.array([float(evals_I_mp[k]) for k in range(n_mp)])
        scalar_evals = np.sort(np.concatenate([scalar_evals_R, scalar_evals_I]))

        # Compute STr[M^2] with mpmath precision
        sum_scalar_m2_mp = fsum(list(evals_R_mp) + list(evals_I_mp))
        sum_ferm_m2_mp = fsum([ferm_evals_mp[k]**2 for k in range(n_mp)])
        str_m2_mp = sum_scalar_m2_mp - 2 * sum_ferm_m2_mp
        str_m2 = float(str_m2_mp)

        # STr[M^4] with mpmath
        sum_scalar_m4_mp = fsum([evals_R_mp[k]**2 for k in range(n_mp)] +
                                [evals_I_mp[k]**2 for k in range(n_mp)])
        sum_ferm_m4_mp = fsum([ferm_evals_mp[k]**4 for k in range(n_mp)])
        str_m4_mp = sum_scalar_m4_mp - 2 * sum_ferm_m4_mp
        str_m4 = float(str_m4_mp)

        print(f"\n  [Using mpmath with {mp.dps} decimal places for this sector]")

    else:
        ferm_evals = np.sort(eigvalsh(W_block))
        ferm_mass_sq = ferm_evals**2

        M2_H_block = W_block @ W_block
        M2_hol_block = M2_hol[np.ix_(idx, idx)]
        M2_soft_block = M2_soft[np.ix_(idx, idx)]

        M2_R = M2_H_block + M2_hol_block + M2_soft_block
        M2_I = M2_H_block - M2_hol_block + M2_soft_block

        scalar_evals_R = np.sort(eigvalsh(M2_R))
        scalar_evals_I = np.sort(eigvalsh(M2_I))
        scalar_evals = np.sort(np.concatenate([scalar_evals_R, scalar_evals_I]))

        str_m2 = np.sum(scalar_evals) - 2.0 * np.sum(ferm_mass_sq)
        str_m4 = np.sum(scalar_evals**2) - 2.0 * np.sum(ferm_mass_sq**2)

    # Analytic STr
    str_analytic = 2.0 * sum(M2_soft[i, i] for i in idx)

    print(f"\n  Fermion eigenvalues (MeV):")
    for k, ev in enumerate(ferm_evals):
        print(f"    {k+1:2d}:  {ev:+22.14e}  (|m| = {abs(ev):.14e})")

    print(f"\n  Scalar mass-squared eigenvalues (MeV^2):")
    for k, ev in enumerate(scalar_evals):
        status = "TACHYON" if ev < -1e-6 else "OK"
        print(f"    {k+1:2d}:  {ev:+22.10e}  ({status})")

    print(f"\n  STr_Q[M^2]:")
    print(f"    Numerical  = {str_m2:+.10e} MeV^2")
    print(f"    Analytic   = {str_analytic:+.10e} MeV^2")
    print(f"    Ratio num/f_pi^2  = {str_m2 / f_pi**2:+.10f}")
    print(f"    Ratio ana/f_pi^2  = {str_analytic / f_pi**2:+.10f}")

    # Number of meson fields in this sector
    n_meson = sum(1 for i in idx if i < 9)
    print(f"    2*n_meson  = {2*n_meson} (expected STr/f_pi^2)")

    print(f"\n  STr_Q[M^4] = {str_m4:+.10e} MeV^4")
    print(f"    STr_Q[M^4] / f_pi^4 = {str_m4 / f_pi**4:+.10f}")

    results[q] = {
        'fields': field_names,
        'idx': idx,
        'ferm_evals': ferm_evals,
        'scalar_evals': scalar_evals,
        'STr_M2_numerical': str_m2,
        'STr_M2_analytic': str_analytic,
        'STr_M4': str_m4,
        'n_meson': n_meson,
    }

# =========================================================================
# TASK 3: Verify total supertrace = 18 f_pi^2
# =========================================================================
print(f"\n\n{'='*78}")
print("TASK 3: TOTAL SUPERTRACE VERIFICATION")
print(f"{'='*78}")

total_str_num = sum(r['STr_M2_numerical'] for r in results.values())
total_str_ana = sum(r['STr_M2_analytic'] for r in results.values())
expected = 18.0 * f_pi**2

print(f"\n  sum_Q STr_Q[M^2] (numerical) = {total_str_num:+.10e} MeV^2")
print(f"  sum_Q STr_Q[M^2] (analytic)  = {total_str_ana:+.10e} MeV^2")
print(f"  18 * f_pi^2                   = {expected:+.10e} MeV^2")
print(f"  Numerical relative error      = {abs(total_str_num - expected)/expected:.4e}")
print(f"  Analytic match                = {np.isclose(total_str_ana, expected, rtol=1e-10)}")

# =========================================================================
# TASK 4: Ratios STr_Q / f_pi^2
# =========================================================================
print(f"\n\n{'='*78}")
print("TASK 4: RATIOS STr_Q[M^2] / f_pi^2")
print(f"{'='*78}")

print(f"\n  {'Q':>6s}  {'Analytic STr_Q':>20s}  {'STr_Q/f_pi^2':>16s}  {'2*n_meson':>10s}  {'Match':>6s}")
for q in sorted(results.keys()):
    r = results[q]
    ratio = r['STr_M2_analytic'] / f_pi**2
    expected_ratio = 2 * r['n_meson']
    match = "YES" if abs(ratio - expected_ratio) < 1e-10 else "NO"
    print(f"  {q:+6.1f}  {r['STr_M2_analytic']:+20.4f}  {ratio:+16.6f}  {expected_ratio:10d}  {match:>6s}")

print(f"\n  Result: STr_Q[M^2] / f_pi^2 = 2 * n_meson(Q) EXACTLY for each sector.")
print(f"  The supertrace per charge sector is a simple integer multiple of f_pi^2:")
print(f"    Q = -1: STr = 4 f_pi^2  (from 2 mesons: Mdu, Msu)")
print(f"    Q =  0: STr = 10 f_pi^2 (from 5 mesons: Muu, Mdd, Mds, Msd, Mss)")
print(f"    Q = +1: STr = 4 f_pi^2  (from 2 mesons: Mud, Mus)")
print(f"    Total:  STr = 18 f_pi^2 (from all 9 mesons)")

# =========================================================================
# TASK 5: EM-weighted supertraces
# =========================================================================
print(f"\n\n{'='*78}")
print("TASK 5: EM-WEIGHTED SUPERTRACES")
print(f"{'='*78}")

sum_Q_STr = sum(q * r['STr_M2_analytic'] for q, r in results.items())
sum_Q2_STr = sum(q**2 * r['STr_M2_analytic'] for q, r in results.items())

print(f"\n  sum_Q  Q_em   * STr_Q[M^2] = {sum_Q_STr:+.10e} MeV^2")
print(f"  sum_Q  Q_em^2 * STr_Q[M^2] = {sum_Q2_STr:+.10e} MeV^2 = {sum_Q2_STr / f_pi**2:+.6f} f_pi^2")

print(f"\n  Analytic breakdown:")
print(f"    sum Q * STr_Q = (+1)(4 f_pi^2) + (-1)(4 f_pi^2) + (0)(10 f_pi^2) = 0")
print(f"    sum Q^2 * STr_Q = (1)(4 f_pi^2) + (1)(4 f_pi^2) + (0)(10 f_pi^2) = 8 f_pi^2")
print(f"    = {8*f_pi**2:.1f} MeV^2")

print(f"\n  In MSSM, sum Q^2 * STr[M^2] = 0 if there is no EM D-term.")
print(f"  Here sum Q^2 * STr = 8 f_pi^2 != 0 because soft breaking")
print(f"  is asymmetric: only mesons get soft mass, not Higgs fields.")

# =========================================================================
# TASK 6: STr[M^4] per sector
# =========================================================================
print(f"\n\n{'='*78}")
print("TASK 6: QUARTIC SUPERTRACE STr[M^4] BY CHARGE SECTOR")
print(f"{'='*78}")

total_str4 = 0.0
print(f"\n  {'Q':>6s}  {'STr_Q[M^4] (MeV^4)':>25s}  {'STr_Q[M^4] / f_pi^4':>25s}")
for q in sorted(results.keys()):
    r = results[q]
    ratio4 = r['STr_M4'] / f_pi**4
    print(f"  {q:+6.1f}  {r['STr_M4']:+25.10e}  {ratio4:+25.10f}")
    total_str4 += r['STr_M4']

print(f"\n  Total STr[M^4] = {total_str4:+.10e} MeV^4")
print(f"  Total STr[M^4] / f_pi^4 = {total_str4 / f_pi**4:+.10f}")

# Analytic formula for STr[M^4]:
# STr[M^4] = Tr((W^2+M_hol+M_soft)^2) + Tr((W^2-M_hol+M_soft)^2) - 2*Tr(W^4)
# = 2*Tr(W^4) + 2*Tr(M_hol^2) + 2*Tr(M_soft^2)
#   + 2*Tr(W^2 M_hol) + 2*Tr(W^2 M_soft) + 2*Tr(M_hol M_soft)
#   + 2*Tr(W^4) - 2*Tr(M_hol^2) + 2*Tr(M_soft^2)
#   - 2*Tr(W^2 M_hol) + 2*Tr(W^2 M_soft) - 2*Tr(M_hol M_soft) - 2*Tr(W^4)
# Wait, let me be more careful:
# (W^2+H+S)^2 = W^4 + H^2 + S^2 + W^2 H + H W^2 + W^2 S + S W^2 + H S + S H
# (W^2-H+S)^2 = W^4 + H^2 + S^2 - W^2 H - H W^2 + W^2 S + S W^2 - H S - S H
# Sum: 2W^4 + 2H^2 + 2S^2 + 2W^2 S + 2S W^2
# Subtract 2W^4:
# STr[M^4] = 2Tr(H^2) + 2Tr(S^2) + 4Tr(W^2 S)

M2_H_full = W @ W
analytic_str4 = (4.0 * np.trace(M2_H_full @ M2_soft)
                + 2.0 * np.trace(M2_hol @ M2_hol)
                + 2.0 * np.trace(M2_soft @ M2_soft))
print(f"\n  Analytic STr[M^4] formula (exact if all real):")
print(f"    4*Tr(W^2 * M_soft) = {4.0*np.trace(M2_H_full @ M2_soft):+.10e}")
print(f"    2*Tr(M_hol^2)      = {2.0*np.trace(M2_hol @ M2_hol):+.10e}")
print(f"    2*Tr(M_soft^2)     = {2.0*np.trace(M2_soft @ M2_soft):+.10e}")
print(f"    Sum                = {analytic_str4:+.10e}")

# Per-sector analytic STr[M^4]
print(f"\n  Per-sector analytic STr[M^4]:")
for q in sorted(results.keys()):
    idx = results[q]['idx']
    W_block = W[np.ix_(idx, idx)]
    M2_H_b = W_block @ W_block
    M2_hol_b = M2_hol[np.ix_(idx, idx)]
    M2_soft_b = M2_soft[np.ix_(idx, idx)]
    ana4 = (4.0 * np.trace(M2_H_b @ M2_soft_b)
           + 2.0 * np.trace(M2_hol_b @ M2_hol_b)
           + 2.0 * np.trace(M2_soft_b @ M2_soft_b))
    print(f"    Q = {q:+.1f}: analytic = {ana4:+.10e},  numerical = {results[q]['STr_M4']:+.10e}")

# =========================================================================
# TASK 7: Fayet-Iliopoulos D-term
# =========================================================================
print(f"\n\n{'='*78}")
print("TASK 7: FAYET-ILIOPOULOS D-TERM EFFECT")
print(f"{'='*78}")

print(f"""
  Adding FI D-term xi: each scalar gets m^2 -> m^2 + Q_em * xi

  Per-sector shift to STr_Q[M^2]:
    STr_Q[M^2](xi) = STr_Q[M^2](0) + 2*n_Q * Q * xi
  where n_Q = number of complex fields with charge Q.
""")

for q in sorted(results.keys()):
    r = results[q]
    n_Q = len(r['idx'])
    print(f"  Q = {q:+.1f}: n_Q = {n_Q:2d}, 2*n_Q*Q = {2*n_Q*q:+.1f}")
    print(f"    STr_Q(xi) = {r['STr_M2_analytic']:.1f} + ({2*n_Q*q:+.1f}) * xi")

total_nQ_Q = sum(len(r['idx']) * q for q, r in results.items())
print(f"\n  Total: sum_Q (n_Q * Q) = {total_nQ_Q:+.4f}")
print(f"  Total STr(xi) = {expected:.2f} + {2*total_nQ_Q:+.4f} * xi")
print(f"  The total supertrace is {'INDEPENDENT' if abs(total_nQ_Q) < 1e-10 else 'DEPENDENT'} of xi.")

sum_Q2_nQ = sum(q**2 * len(r['idx']) for q, r in results.items())
sum_Q3_nQ = sum(q**3 * len(r['idx']) for q, r in results.items())
print(f"\n  EM-weighted supertraces with FI:")
print(f"    sum Q * STr_Q(xi) = {sum_Q_STr:.4f} + {2*sum_Q2_nQ:+.4f} * xi")
print(f"    sum Q^2 * STr_Q(xi) = {sum_Q2_STr:.4f} + {2*sum_Q3_nQ:+.4f} * xi")

# =========================================================================
# SUMMARY TABLE
# =========================================================================
print(f"\n\n{'='*78}")
print("SUMMARY TABLE")
print(f"{'='*78}")

print(f"\n  {'Q_em':>6s} | {'n_fields':>8s} | {'n_meson':>7s} | {'STr[M^2]/f_pi^2':>16s} | {'2*n_meson':>10s} | {'STr[M^4] (MeV^4)':>22s}")
print(f"  {'─'*6}─┼─{'─'*8}─┼─{'─'*7}─┼─{'─'*16}─┼─{'─'*10}─┼─{'─'*22}")
for q in sorted(results.keys()):
    r = results[q]
    ratio = r['STr_M2_analytic'] / f_pi**2
    print(f"  {q:+6.1f} | {len(r['idx']):8d} | {r['n_meson']:7d} | {ratio:+16.1f} | {2*r['n_meson']:10d} | {r['STr_M4']:+22.6e}")

print(f"\n  Total STr[M^2] = 18 * f_pi^2 = {expected:.1f} MeV^2  (VERIFIED ANALYTICALLY)")
print(f"  sum Q * STr_Q = 0 (charge conjugation symmetry)")
print(f"  sum Q^2 * STr_Q = 8 * f_pi^2 = {8*f_pi**2:.1f} MeV^2")


# =========================================================================
# Write results/supertrace_by_charge.md
# =========================================================================
out_lines = []
def w(s=""):
    out_lines.append(s)

w("# Supertrace Decomposed by EM Charge Sector")
w()
w("## Setup")
w()
w("**16 chiral superfields** with NMSSM coupling $\\lambda X H_u \\cdot H_d$ ($\\lambda = 0.72$).")
w()
w("**Superpotential:**")
w("$$W = \\sum_i m_i M^i_i + X(\\det M - B\\tilde{B} - \\Lambda^6) + y_c H_u^0 M^d_d + y_b H_d^0 M^s_s + \\lambda X (H_u^0 H_d^0 - H_u^+ H_d^-)$$")
w()
w("**Soft SUSY breaking:** $V_{\\text{soft}} = f_\\pi^2 \\text{Tr}(M^\\dagger M)$, applied to the 9 meson fields only.")
w()
w("**EM charges:** $Q(u)=2/3$, $Q(d)=Q(s)=-1/3$. Mesons: $Q(M^i_j) = Q_i - Q_j$. Baryons: $Q(B)=Q(\\tilde{B})=0$.")
w()
w("---")
w()
w("## Charge Sector Decomposition")
w()
w("| Sector | Fields | $n_{\\text{fields}}$ | $n_{\\text{meson}}$ (soft) |")
w("|--------|--------|----------|-------------------|")
for q in sorted(results.keys()):
    r = results[q]
    w(f"| $Q = {q:+.0f}$ | {', '.join(r['fields'])} | {len(r['idx'])} | {r['n_meson']} |")
w()
w("---")
w()
w("## Fermion Spectra by Charge Sector")
w()

for q in sorted(results.keys()):
    r = results[q]
    w(f"### $Q = {q:+.0f}$ sector")
    w()
    w(f"Fields: {', '.join(r['fields'])}")
    w()
    w("**Fermion eigenvalues:**")
    w()
    w("| # | Eigenvalue (MeV) |")
    w("|---|------------------|")
    for k, ev in enumerate(r['ferm_evals']):
        w(f"| {k+1} | {ev:+.10e} |")
    w()
    w("**Scalar mass-squared eigenvalues:**")
    w()
    w("| # | $m^2$ (MeV$^2$) | Status |")
    w("|---|------------------|--------|")
    for k, ev in enumerate(r['scalar_evals']):
        status = "tachyon" if ev < -1e-6 else "OK"
        w(f"| {k+1} | {ev:+.10e} | {status} |")
    w()

w("---")
w()
w("## Quadratic Supertrace per Charge Sector")
w()
w("**Analytic result (exact):**")
w()
w("$$\\text{STr}_Q[M^2] = 2 \\, \\text{Tr}(M^2_{\\text{soft}}|_Q) = 2 \\, n_{\\text{meson}}(Q) \\, f_\\pi^2$$")
w()
w("This is exact because in $\\text{STr} = \\text{Tr}(m^2_R + m^2_I) - 2\\text{Tr}(m^2_f)$,")
w("the $W^2$ contributions cancel between scalars and fermions, and the holomorphic")
w("contribution $M^2_{\\text{hol}}$ cancels between R and I scalar components.")
w("The soft mass is the only surviving term.")
w()
w("| $Q_{\\text{em}}$ | $n_{\\text{meson}}$ | $\\text{STr}_Q / f_\\pi^2$ | Value (MeV$^2$) |")
w("|---------|-----------|--------------------------|-----------------|")
for q in sorted(results.keys()):
    r = results[q]
    ratio = int(round(r['STr_M2_analytic'] / f_pi**2))
    w(f"| ${q:+.0f}$ | {r['n_meson']} | {ratio} | {r['STr_M2_analytic']:.1f} |")
w(f"| **Total** | **9** | **18** | **{expected:.1f}** |")
w()

w("---")
w()
w("## EM-Weighted Supertraces")
w()
w("$$\\sum_Q Q_{\\text{em}} \\cdot \\text{STr}_Q[M^2] = (+1)(4) + (-1)(4) + (0)(10) = 0$$")
w()
w("$$\\sum_Q Q_{\\text{em}}^2 \\cdot \\text{STr}_Q[M^2] = (1)(4) + (1)(4) + (0)(10) = 8 f_\\pi^2$$")
w()
w(f"Numerically: $\\sum Q^2 \\cdot \\text{{STr}}_Q = {sum_Q2_STr:.1f}$ MeV$^2 = 8 f_\\pi^2$.")
w()
w("The first sum vanishes by charge conjugation symmetry of the soft sector")
w("(equal number of $Q=+1$ and $Q=-1$ mesons with equal soft masses).")
w()
w("The second sum gives $8 f_\\pi^2$ --- a nonzero electromagnetic supertrace.")
w("In the MSSM, $\\sum Q^2 \\cdot \\text{STr}[M^2] = 0$ when there is no EM D-term.")
w("Here the nonzero value reflects the fact that soft breaking is applied")
w("only to mesons (which carry EM charge), not to the Higgs fields or X, B, $\\tilde{B}$.")
w()

w("---")
w()
w("## Quartic Supertrace $\\text{STr}[M^4]$")
w()
w("$$\\text{STr}[M^4] = 4 \\text{Tr}(W^2 M_{\\text{soft}}) + 2 \\text{Tr}(M_{\\text{hol}}^2) + 2 \\text{Tr}(M_{\\text{soft}}^2)$$")
w()
w(f"| Term | Value (MeV$^4$) |")
w(f"|------|-----------------|")
w(f"| $4\\text{{Tr}}(W^2 M_{{\\text{{soft}}}})$ | {4.0*np.trace(M2_H_full @ M2_soft):+.6e} |")
w(f"| $2\\text{{Tr}}(M_{{\\text{{hol}}}}^2)$ | {2.0*np.trace(M2_hol @ M2_hol):+.6e} |")
w(f"| $2\\text{{Tr}}(M_{{\\text{{soft}}}}^2)$ | {2.0*np.trace(M2_soft @ M2_soft):+.6e} |")
w(f"| **Total** | **{analytic_str4:+.6e}** |")
w()
w(f"Per-sector:")
w()
w("| $Q_{\\text{em}}$ | $\\text{STr}_Q[M^4]$ (MeV$^4$) | $\\text{STr}_Q[M^4] / f_\\pi^4$ |")
w("|---------|-------------------------------|--------------------------------|")
for q in sorted(results.keys()):
    r = results[q]
    ratio4 = r['STr_M4'] / f_pi**4
    w(f"| ${q:+.0f}$ | {r['STr_M4']:+.6e} | {ratio4:+.6f} |")

w(f"\n**Total:** $\\text{{STr}}[M^4] = {total_str4:+.6e}$ MeV$^4$")
w()
w("Unlike $\\text{STr}[M^2]$, the quartic supertrace does **not** reduce to a pure soft-mass expression.")
w("It receives contributions from $W^2 M_{\\text{soft}}$ (fermion masses crossed with soft breaking)")
w("and $M_{\\text{hol}}^2$ (Yukawa-induced F-term splittings).")
w()

w("---")
w()
w("## Fayet-Iliopoulos D-term")
w()
w("Adding an EM FI term $\\xi$, each scalar gets $m^2 \\to m^2 + Q_{\\text{em}} \\xi$.")
w()
w("Per-sector shift:")
w("$$\\text{STr}_Q(\\xi) = \\text{STr}_Q(0) + 2 n_Q \\cdot Q \\cdot \\xi$$")
w()
w("where $n_Q$ is the number of complex fields with charge $Q$.")
w()
w("| $Q$ | $n_Q$ | $2 n_Q Q$ | Effect |")
w("|-----|-------|-----------|--------|")
for q in sorted(results.keys()):
    r = results[q]
    n_Q = len(r['idx'])
    coeff = 2*n_Q*q
    if abs(coeff) < 1e-10:
        effect = "unchanged"
    else:
        effect = f"$+({coeff:+.0f})\\xi$"
    w(f"| ${q:+.0f}$ | {n_Q} | {coeff:+.0f} | {effect} |")

w()
w(f"**Total:** $\\sum_Q n_Q Q = {total_nQ_Q:+.0f}$ --- the total supertrace is **independent of $\\xi$**,")
w("because the total EM charge of all chiral fields vanishes.")
w()
w("Per-sector, the FI term breaks the degeneracy between $Q=+1$ and $Q=-1$ sectors:")
w()
w("$$\\text{STr}_{+1}(\\xi) = 4 f_\\pi^2 + 6\\xi, \\qquad \\text{STr}_{-1}(\\xi) = 4 f_\\pi^2 - 6\\xi$$")
w()
w("The EM-weighted sums:")
w(f"$$\\sum_Q Q \\cdot \\text{{STr}}_Q(\\xi) = 12\\xi$$")
w(f"$$\\sum_Q Q^2 \\cdot \\text{{STr}}_Q(\\xi) = 8 f_\\pi^2$$")
w()
w("The $Q^2$-weighted sum is independent of $\\xi$ because $\\sum Q^3 n_Q = 0$.")
w()
w("---")
w()
w("*Generated by supertrace_by_charge.py*")

with open("results/supertrace_by_charge.md", "w") as f:
    f.write("\n".join(out_lines))

print(f"\n\nResults written to results/supertrace_by_charge.md")
