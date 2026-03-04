#!/usr/bin/env python3
"""
Full numerical minimization of the SQCD + Higgs scalar potential.

Superpotential:
    W = sum_i m_i M^i_i + X(det M - B B~ - Lambda^6)
        + y_c H_u^0 M^d_d + y_b H_d^0 M^s_s
        + lambda X (H_u^0 H_d^0 - H_u^+ H_d^-)

Scalar potential:
    V = sum_I |F_I|^2 + f_pi^2 Tr(M^dag M)

Variables: 20 real = 9 complex meson entries + complex X
           (Higgs fixed at v/sqrt(2), B = B~ = 0)
"""

import numpy as np
from scipy.optimize import minimize, differential_evolution, basinhopping
from scipy.linalg import svd
import sys
import time
import warnings
warnings.filterwarnings('ignore')

# =====================================================================
# Physical inputs (MeV)
# =====================================================================
m_u = 2.16
m_d = 4.67
m_s = 93.4
m_c = 1270.0
m_b = 4180.0
Lambda = 300.0
Lambda6 = Lambda**6
v_ew = 246220.0

y_c = 2.0 * m_c / v_ew
y_b = 2.0 * m_b / v_ew
lam = 0.72

f_pi = 92.0
m_tilde_sq = f_pi**2

Hu0_fixed = v_ew / np.sqrt(2)
Hd0_fixed = v_ew / np.sqrt(2)

# Seesaw vacuum
masses_q = np.array([m_u, m_d, m_s])
C = Lambda**2 * np.prod(masses_q)**(1.0/3.0)
M_diag = C / masses_q
X_diag = -m_u / (M_diag[1] * M_diag[2])
det_diag = np.prod(M_diag)
FX_seesaw = lam * Hu0_fixed * Hd0_fixed

print("=" * 78)
print("  FULL CKM MINIMIZATION: SQCD + HIGGS + MU-TERM")
print("=" * 78)
print()
print("Physical inputs:")
print(f"  m_u = {m_u}, m_d = {m_d}, m_s = {m_s} MeV")
print(f"  Lambda = {Lambda} MeV, Lambda^6 = {Lambda6:.6e}")
print(f"  v/sqrt(2) = {v_ew/np.sqrt(2):.4f} MeV")
print(f"  y_c = {y_c:.8e}, y_b = {y_b:.8e}, lambda = {lam}")
print(f"  f_pi = {f_pi} MeV, m_tilde^2 = {m_tilde_sq}")
print()
print("Seesaw vacuum:")
for i, fl in enumerate(['u', 'd', 's']):
    print(f"  M_{fl} = {M_diag[i]:.4f}")
print(f"  X_0 = {X_diag:.6e}")
print(f"  F_X(seesaw, Higgs on) = lam*v^2/2 = {FX_seesaw:.6e}")
print()
sys.stdout.flush()

# =====================================================================
# Variable packing: 20 real = M(3x3 complex) + X(complex)
# =====================================================================
X_scale = abs(X_diag) if abs(X_diag) > 1e-20 else 1e-9
M_scale = np.zeros((3, 3))
for i in range(3):
    for j in range(3):
        M_scale[i, j] = max(M_diag[i] if i == j else np.sqrt(M_diag[i]*M_diag[j]), 1.0)

SCALES = np.zeros(20)
for i in range(3):
    for j in range(3):
        idx = 2*(3*i+j)
        SCALES[idx] = M_scale[i, j]
        SCALES[idx+1] = M_scale[i, j]
SCALES[18] = X_scale
SCALES[19] = X_scale


def unpack(params):
    p = params * SCALES
    M = np.zeros((3, 3), dtype=complex)
    for i in range(3):
        for j in range(3):
            idx = 2*(3*i+j)
            M[i, j] = p[idx] + 1j*p[idx+1]
    X = p[18] + 1j*p[19]
    return M, X


def cofactor_matrix(M):
    cof = np.empty((3, 3), dtype=complex)
    cof[0, 0] = M[1, 1]*M[2, 2] - M[1, 2]*M[2, 1]
    cof[1, 1] = M[0, 0]*M[2, 2] - M[0, 2]*M[2, 0]
    cof[2, 2] = M[0, 0]*M[1, 1] - M[0, 1]*M[1, 0]
    cof[1, 0] = -(M[0, 1]*M[2, 2] - M[0, 2]*M[2, 1])
    cof[0, 1] = -(M[1, 0]*M[2, 2] - M[1, 2]*M[2, 0])
    cof[2, 0] = M[0, 1]*M[1, 2] - M[0, 2]*M[1, 1]
    cof[0, 2] = M[1, 0]*M[2, 1] - M[1, 1]*M[2, 0]
    cof[2, 1] = -(M[0, 0]*M[1, 2] - M[0, 2]*M[1, 0])
    cof[1, 2] = -(M[0, 0]*M[2, 1] - M[0, 1]*M[2, 0])
    return cof


def det3(M):
    return (M[0, 0]*(M[1, 1]*M[2, 2] - M[1, 2]*M[2, 1])
          - M[0, 1]*(M[1, 0]*M[2, 2] - M[1, 2]*M[2, 0])
          + M[0, 2]*(M[1, 0]*M[2, 1] - M[1, 1]*M[2, 0]))


def scalar_potential(params):
    M, X = unpack(params)
    cof = cofactor_matrix(M)
    detM = det3(M)
    mass_vec = np.array([m_u, m_d, m_s])

    V = 0.0
    for i in range(3):
        for j in range(3):
            Fij = X * cof[j, i]
            if i == j:
                Fij += mass_vec[i]
            if i == 1 and j == 1:
                Fij += y_c * Hu0_fixed
            if i == 2 and j == 2:
                Fij += y_b * Hd0_fixed
            V += abs(Fij)**2

    FX = detM - Lambda6 + lam * Hu0_fixed * Hd0_fixed
    V += abs(FX)**2

    F_Hu0 = y_c * M[1, 1] + lam * X * Hd0_fixed
    F_Hd0 = y_b * M[2, 2] + lam * X * Hu0_fixed
    V += abs(F_Hu0)**2 + abs(F_Hd0)**2

    V += m_tilde_sq * np.sum(np.abs(M)**2).real
    return V.real


def scalar_potential_detailed(params):
    M, X = unpack(params)
    cof = cofactor_matrix(M)
    detM = det3(M)
    mass_vec = np.array([m_u, m_d, m_s])

    FM = np.zeros((3, 3), dtype=complex)
    for i in range(3):
        for j in range(3):
            FM[i, j] = X * cof[j, i]
            if i == j:
                FM[i, j] += mass_vec[i]
            if i == 1 and j == 1:
                FM[i, j] += y_c * Hu0_fixed
            if i == 2 and j == 2:
                FM[i, j] += y_b * Hd0_fixed

    FX = detM - Lambda6 + lam * Hu0_fixed * Hd0_fixed
    F_Hu0 = y_c * M[1, 1] + lam * X * Hd0_fixed
    F_Hd0 = y_b * M[2, 2] + lam * X * Hu0_fixed

    V = np.sum(np.abs(FM)**2)
    V += abs(FX)**2
    V += abs(F_Hu0)**2 + abs(F_Hd0)**2
    V += m_tilde_sq * np.sum(np.abs(M)**2)

    return V.real, {
        'FM': FM, 'FX': FX, 'F_Hu0': F_Hu0, 'F_Hd0': F_Hd0,
        'M': M.copy(), 'X': X, 'detM': detM
    }


# =====================================================================
# (a) Starting point
# =====================================================================
print("=" * 78)
print("  (a) STARTING POINT")
print("=" * 78)
print()

p0 = np.zeros(20)
for i in range(3):
    idx = 2*(3*i+i)
    p0[idx] = M_diag[i] / SCALES[idx]
p0[18] = X_diag / SCALES[18]

V_start, info_start = scalar_potential_detailed(p0)
print(f"V(seesaw + Higgs) = {V_start:.6e} MeV^2")

FM0 = info_start['FM']
V_FM = np.sum(np.abs(FM0)**2).real
V_FX = abs(info_start['FX'])**2
V_FH = abs(info_start['F_Hu0'])**2 + abs(info_start['F_Hd0'])**2
V_soft = m_tilde_sq * np.sum(np.abs(info_start['M'])**2).real
print(f"  |F_M|^2 = {V_FM:.4e}, |F_X|^2 = {V_FX:.4e}, |F_H|^2 = {V_FH:.4e}, soft = {V_soft:.4e}")
print(f"  F_X dominates: lam*v^2/2 = {FX_seesaw:.4e}")
print()
sys.stdout.flush()


# =====================================================================
# (b) Multi-strategy minimization
# =====================================================================
print("=" * 78)
print("  (b) MINIMIZATION (20 real variables, 40 restarts + global methods)")
print("=" * 78)
print()
sys.stdout.flush()

t0 = time.time()


def minimize_fast(func, x0):
    """L-BFGS-B then Powell."""
    best_x, best_v = x0.copy(), func(x0)
    try:
        r = minimize(func, best_x, method='L-BFGS-B',
                     options={'maxiter': 30000, 'maxfun': 200000, 'ftol': 1e-22})
        if r.fun < best_v: best_x, best_v = r.x.copy(), r.fun
    except Exception: pass
    try:
        r = minimize(func, best_x, method='Powell',
                     options={'maxiter': 50000, 'maxfev': 300000, 'ftol': 1e-22, 'xtol': 1e-13})
        if r.fun < best_v: best_x, best_v = r.x.copy(), r.fun
    except Exception: pass
    return best_v, best_x


def minimize_full(func, x0):
    """L-BFGS-B -> Powell -> Nelder-Mead."""
    best_v, best_x = minimize_fast(func, x0)
    try:
        r = minimize(func, best_x, method='Nelder-Mead',
                     options={'maxiter': 200000, 'maxfev': 500000,
                              'xatol': 1e-14, 'fatol': 1e-22, 'adaptive': True})
        if r.fun < best_v: best_x, best_v = r.x.copy(), r.fun
    except Exception: pass
    return best_v, best_x


# From seesaw
V_seesaw, p_seesaw = minimize_fast(scalar_potential, p0)
print(f"From seesaw: V = {V_seesaw:.6e}")
all_results = [(V_seesaw, p_seesaw.copy(), 'seesaw')]

# Random restarts
np.random.seed(42)
scales = [0.01, 0.05, 0.10, 0.50]
for trial in range(40):
    scale = scales[trial % 4]
    p_trial = p0.copy()

    for i in range(3):
        for j in range(3):
            idx = 2*(3*i+j)
            if i != j:
                ref = np.sqrt(M_diag[i]*M_diag[j]) / SCALES[idx]
                p_trial[idx] = np.random.randn() * scale * ref
                p_trial[idx+1] = np.random.randn() * scale * ref * 0.3
            else:
                p_trial[idx] *= (1 + np.random.randn() * scale)
                p_trial[idx+1] = np.random.randn() * scale * 0.1 * abs(p_trial[idx])

    p_trial[18] *= (1 + np.random.randn() * scale)
    p_trial[19] = np.random.randn() * scale * 0.1 * abs(p_trial[18])

    V_t, p_t = minimize_fast(scalar_potential, p_trial)
    all_results.append((V_t, p_t.copy(), f'restart_{trial}_s{scale}'))

    if (trial+1) % 10 == 0:
        print(f"  {trial+1}/40 done, best: {min(r[0] for r in all_results):.4e}")
        sys.stdout.flush()

# Basin hopping
print("Basin hopping...")
sys.stdout.flush()
try:
    class StepTaker:
        def __init__(self, s=0.2): self.s = s
        def __call__(self, x):
            x_new = x.copy()
            for i in range(3):
                for j in range(3):
                    idx = 2*(3*i+j)
                    if i != j:
                        x_new[idx] += np.random.randn() * self.s
                        x_new[idx+1] += np.random.randn() * self.s * 0.3
                    else:
                        x_new[idx] *= (1 + np.random.randn() * self.s * 0.1)
            x_new[18] *= (1 + np.random.randn() * self.s * 0.1)
            return x_new

    bh = basinhopping(scalar_potential, p0, niter=25,
                       minimizer_kwargs={'method': 'L-BFGS-B',
                                         'options': {'maxiter': 20000, 'ftol': 1e-20}},
                       take_step=StepTaker(0.3), seed=123)
    all_results.append((bh.fun, bh.x.copy(), 'basin_hop'))
    print(f"  V = {bh.fun:.4e}")
except Exception as e:
    print(f"  Failed: {e}")

# Differential evolution
print("Differential evolution...")
sys.stdout.flush()
try:
    bounds = []
    for i in range(3):
        for j in range(3):
            idx = 2*(3*i+j)
            if i == j:
                c = M_diag[i] / SCALES[idx]
                w = abs(c) * 3.0
                bounds.append((c - w, c + w))
                bounds.append((-w*0.3, w*0.3))
            else:
                ref = np.sqrt(M_diag[i]*M_diag[j]) / SCALES[idx]
                bounds.append((-ref*1.5, ref*1.5))
                bounds.append((-ref*0.5, ref*0.5))
    cX = X_diag / SCALES[18]
    wX = abs(cX) * 5.0
    bounds.append((cX - wX, cX + wX))
    bounds.append((-wX*0.3, wX*0.3))

    de = differential_evolution(scalar_potential, bounds,
                                 maxiter=400, popsize=40, tol=1e-18,
                                 seed=456, polish=True)
    all_results.append((de.fun, de.x.copy(), 'diff_evol'))
    print(f"  V = {de.fun:.4e}")
except Exception as e:
    print(f"  Failed: {e}")

all_results.sort(key=lambda x: x[0])
print(f"\nBest 10 minima:")
for i, (V, _, lab) in enumerate(all_results[:10]):
    print(f"  {i+1}. V = {V:.10e}  ({lab})")
print()

# Refine best
print("Refining best (3 rounds of full chain)...")
sys.stdout.flush()
p_best = all_results[0][1].copy()
for rnd in range(3):
    V_best, p_best = minimize_full(scalar_potential, p_best)
    print(f"  Round {rnd+1}: V = {V_best:.10e}")
    sys.stdout.flush()

V_best_final, info_best = scalar_potential_detailed(p_best)
t_opt = time.time() - t0
print(f"\nBest: V = {V_best_final:.10e}  ({t_opt:.0f}s)")
print()
sys.stdout.flush()


# =====================================================================
# (c-d) CKM angles from SVD + PDG comparison
# =====================================================================
print("=" * 78)
print("  (c-d) CKM MIXING ANGLES")
print("=" * 78)
print()

labels_f = ['u', 'd', 's']
M_min = info_best['M']

print("Meson matrix at minimum:")
print(f"{'':>5}", end="")
for j in range(3):
    print(f"  {labels_f[j]+' Re':>14} {labels_f[j]+' Im':>14}", end="")
print()
for i in range(3):
    print(f"  {labels_f[i]:>3}", end="")
    for j in range(3):
        print(f"  {M_min[i,j].real:14.4f} {M_min[i,j].imag:14.4f}", end="")
    print()
print()

print("|M^i_j| (MeV):")
for i in range(3):
    row = f"  {labels_f[i]:>3}"
    for j in range(3):
        row += f"  {abs(M_min[i,j]):14.4f}"
    print(row)
print()

U_L, sigma_vals, Vt_R = svd(M_min)
V_R = Vt_R.conj().T

print(f"Singular values: {sigma_vals[0]:.4f}, {sigma_vals[1]:.4f}, {sigma_vals[2]:.4f}")
print(f"Seesaw values:   {M_diag[0]:.4f}, {M_diag[1]:.4f}, {M_diag[2]:.4f}")
for i in range(3):
    print(f"  sigma_{i+1}/M_{labels_f[i]} = {sigma_vals[i]/M_diag[i]:.6f}")
print()

V_CKM = U_L.conj().T @ V_R
if np.linalg.det(V_CKM).real < 0:
    V_CKM[:, -1] *= -1

print("|V_CKM|:")
for i in range(3):
    print("  [" + "  ".join(f"{abs(V_CKM[i,j]):10.6f}" for j in range(3)) + "  ]")
print()

s13 = abs(V_CKM[0, 2])
theta_13 = np.arcsin(np.clip(s13, 0, 1))
c13 = np.cos(theta_13)
s12 = abs(V_CKM[0, 1]) / c13 if c13 > 1e-10 else 0
s23 = abs(V_CKM[1, 2]) / c13 if c13 > 1e-10 else 0
theta_12 = np.arcsin(np.clip(s12, 0, 1))
theta_23 = np.arcsin(np.clip(s23, 0, 1))

theta_C_oakes = np.degrees(np.arctan(np.sqrt(m_d/m_s)))

J_CKM = abs(np.imag(V_CKM[0,0]*V_CKM[1,1]*np.conj(V_CKM[0,1])*np.conj(V_CKM[1,0])))

print("CKM angles (degrees):")
print(f"  {'':>12} {'This':>10} {'PDG':>10} {'Oakes':>10} {'ratio':>10}")
print(f"  {'theta_12':>12} {np.degrees(theta_12):10.4f} {13.04:10.2f} {theta_C_oakes:10.2f} {np.degrees(theta_12)/13.04:10.4f}")
print(f"  {'theta_23':>12} {np.degrees(theta_23):10.4f} {2.38:10.2f} {'---':>10} {np.degrees(theta_23)/2.38:10.4f}")
print(f"  {'theta_13':>12} {np.degrees(theta_13):10.4f} {0.201:10.3f} {'---':>10} {np.degrees(theta_13)/0.201:10.4f}")
print(f"  {'J':>12} {J_CKM:10.4e} {3.18e-5:10.2e}")
print()

# Quark masses
print("Seesaw masses from SVD (m_q = C/sigma_q):")
for i in range(3):
    m_q = C / sigma_vals[i] if sigma_vals[i] > 0 else float('inf')
    print(f"  sigma_{i+1} = {sigma_vals[i]:.2f} -> m = {m_q:.4f} MeV (input {masses_q[i]:.2f})")
print()

print("Yukawa masses:")
print(f"  y_c * |M^d_d| = {y_c*abs(M_min[1,1]):.2f} MeV (target m_c = {m_c})")
print(f"  y_b * |M^s_s| = {y_b*abs(M_min[2,2]):.2f} MeV (target m_b = {m_b})")
print()
sys.stdout.flush()


# =====================================================================
# (e) Hessian analysis
# =====================================================================
print("=" * 78)
print("  (e) HESSIAN AT MINIMUM")
print("=" * 78)
print()

def numerical_hessian(func, x0, h=1e-5):
    n = len(x0)
    H = np.zeros((n, n))
    for i in range(n):
        for j in range(i, n):
            xpp = x0.copy(); xpm = x0.copy(); xmp = x0.copy(); xmm = x0.copy()
            hi = h * max(abs(x0[i]), 1.0)
            hj = h * max(abs(x0[j]), 1.0)
            xpp[i] += hi; xpp[j] += hj
            xpm[i] += hi; xpm[j] -= hj
            xmp[i] -= hi; xmp[j] += hj
            xmm[i] -= hi; xmm[j] -= hj
            H[i, j] = (func(xpp) - func(xpm) - func(xmp) + func(xmm)) / (4*hi*hj)
            H[j, i] = H[i, j]
    return H

print("Computing Hessian...")
sys.stdout.flush()
H_scaled = numerical_hessian(scalar_potential, p_best, h=1e-5)
eigvals_scaled = np.sort(np.linalg.eigvalsh(H_scaled))

# Physical Hessian
H_phys = H_scaled / np.outer(SCALES, SCALES)
eigvals_phys, eigvecs_phys = np.linalg.eigh(H_phys)
idx_sort = np.argsort(eigvals_phys)
eigvals_phys = eigvals_phys[idx_sort]
eigvecs_phys = eigvecs_phys[:, idx_sort]

# Use a more reasonable threshold for stability
# The physical eigenvalues range from ~10^6 to ~10^21
# Use absolute threshold: if |eigenvalue| < 1.0 it's effectively zero
abs_threshold = 1.0  # 1 MeV^{-2} or similar mixed units

n_neg_phys = sum(1 for ev in eigvals_phys if ev < -abs_threshold)
n_zero_phys = sum(1 for ev in eigvals_phys if abs(ev) <= abs_threshold)
n_pos_phys = 20 - n_neg_phys - n_zero_phys

print("Physical Hessian eigenvalues (d^2V/dphi_i dphi_j):")
variable_names = []
for i in range(3):
    for j in range(3):
        variable_names.append(f"Re(M^{labels_f[i]}_{labels_f[j]})")
        variable_names.append(f"Im(M^{labels_f[i]}_{labels_f[j]})")
variable_names.append("Re(X)")
variable_names.append("Im(X)")

for i in range(20):
    ev = eigvals_phys[i]
    if ev < -abs_threshold:
        status = "TACHYONIC"
    elif abs(ev) <= abs_threshold:
        status = "~ZERO"
    else:
        status = "positive"
    print(f"  m^2_{i+1:2d} = {ev:14.6e}  {status}")

print()
print(f"Stability: {n_neg_phys} tachyonic, {n_zero_phys} flat, {n_pos_phys} positive")
if n_neg_phys == 0:
    print("=> TRUE LOCAL MINIMUM")
else:
    print(f"=> SADDLE POINT with {n_neg_phys} tachyonic direction(s)")
print()
sys.stdout.flush()


# =====================================================================
# (f) Tachyonic mode analysis
# =====================================================================
print("=" * 78)
print("  (f) TACHYONIC MODES")
print("=" * 78)
print()

tach_count = 0
tach_modes = []
for i in range(20):
    ev = eigvals_phys[i]
    if ev < -abs_threshold:
        tach_count += 1
        vec = eigvecs_phys[:, i]
        tach_modes.append((ev, vec, i))
        print(f"Mode {tach_count}: m^2 = {ev:.6e}")
        sorted_k = np.argsort(-np.abs(vec))
        for k in sorted_k[:5]:
            if abs(vec[k]) > 0.05:
                print(f"    {variable_names[k]:>25}: {vec[k]:+.4f}")
        print()

if tach_count == 0:
    print("No tachyonic modes. Showing softest 5 modes:")
    for i in range(min(5, 20)):
        ev = eigvals_phys[i]
        vec = eigvecs_phys[:, i]
        print(f"  Mode {i+1}: m^2 = {ev:.6e}")
        sorted_k = np.argsort(-np.abs(vec))
        for k in sorted_k[:3]:
            print(f"      {variable_names[k]:>25}: {vec[k]:+.4f}")
        print()

# Check specific direction: Im(0.87 M^d_d + 0.49 M^u_u)
dir_test = np.zeros(20)
dir_test[2*(3*1+1)+1] = 0.87  # Im(M^d_d)
dir_test[2*(3*0+0)+1] = 0.49  # Im(M^u_u)
dir_test /= np.linalg.norm(dir_test)
curv = dir_test @ H_phys @ dir_test
print(f"Curvature along Im(0.87 M^d_d + 0.49 M^u_u): {curv:.6e}")
print()
sys.stdout.flush()


# =====================================================================
# Explore along tachyonic/soft directions (limited search)
# =====================================================================
print("=" * 78)
print("  EXPLORING ALONG SOFTEST EIGENVECTORS")
print("=" * 78)
print()

# Use the 3 softest scaled eigenvectors with a few displacement magnitudes
_, evecs_scaled = np.linalg.eigh(H_scaled)
displacements = [0.01, 0.1, 0.5, 1.0, 2.0]

results_disp = []
n_explore = min(3, 20)

for mode in range(n_explore):
    evec = evecs_scaled[:, mode]
    for d in displacements:
        for sign in [+1, -1]:
            p_d = p_best + sign * d * evec
            V_d, p_d_min = minimize_fast(scalar_potential, p_d)
            results_disp.append((V_d, p_d_min.copy(), f'mode{mode}_d{sign*d:.2f}'))

results_disp.sort(key=lambda x: x[0])
print(f"Best 5 from directional search:")
for i, (V, _, lab) in enumerate(results_disp[:5]):
    marker = " ***" if V < V_best_final * 0.999 else ""
    print(f"  {i+1}. V = {V:.10e}  ({lab}){marker}")
print()

if results_disp[0][0] < V_best_final * (1 - 1e-10):
    print(f"IMPROVED: {V_best_final:.6e} -> {results_disp[0][0]:.6e}")
    p_best2 = results_disp[0][1].copy()
    for _ in range(3):
        V_best2, p_best2 = minimize_full(scalar_potential, p_best2)

    V_best2_final, info_best2 = scalar_potential_detailed(p_best2)
    print(f"Refined: V = {V_best2_final:.10e}")

    # CKM at improved point
    M_min2 = info_best2['M']
    U2, sig2, Vt2 = svd(M_min2)
    V2 = U2.conj().T @ Vt2.conj().T
    if np.linalg.det(V2).real < 0:
        V2[:, -1] *= -1
    s13_2 = abs(V2[0, 2])
    t13_2 = np.arcsin(np.clip(s13_2, 0, 1))
    c13_2 = np.cos(t13_2)
    s12_2 = abs(V2[0, 1]) / c13_2 if c13_2 > 1e-10 else 0
    s23_2 = abs(V2[1, 2]) / c13_2 if c13_2 > 1e-10 else 0
    t12_2 = np.arcsin(np.clip(s12_2, 0, 1))
    t23_2 = np.arcsin(np.clip(s23_2, 0, 1))

    print(f"|V_CKM| at improved minimum:")
    for i in range(3):
        print("  [" + "  ".join(f"{abs(V2[i,j]):10.6f}" for j in range(3)) + "  ]")
    print(f"  theta_12 = {np.degrees(t12_2):.4f}, theta_23 = {np.degrees(t23_2):.4f}, theta_13 = {np.degrees(t13_2):.4f}")
    print()

    # Update
    V_best_final, info_best, p_best = V_best2_final, info_best2, p_best2
    M_min, sigma_vals, V_CKM = M_min2, sig2, V2
    theta_12, theta_23, theta_13 = t12_2, t23_2, t13_2

    # Recompute Hessian
    print("Recomputing Hessian at improved minimum...")
    H_scaled = numerical_hessian(scalar_potential, p_best, h=1e-5)
    H_phys = H_scaled / np.outer(SCALES, SCALES)
    eigvals_phys = np.sort(np.linalg.eigvalsh(H_phys))
    n_neg_phys = sum(1 for ev in eigvals_phys if ev < -abs_threshold)
    n_zero_phys = sum(1 for ev in eigvals_phys if abs(ev) <= abs_threshold)
    print(f"  {n_neg_phys} tachyonic, {n_zero_phys} flat, {20-n_neg_phys-n_zero_phys} positive")
    print()
else:
    print("No improvement found.")
    print()

sys.stdout.flush()


# =====================================================================
# Collect distinct minima and their CKM angles
# =====================================================================
print("=" * 78)
print("  DISTINCT MINIMA CATALOG")
print("=" * 78)
print()

all_combined = all_results + results_disp
all_combined.sort(key=lambda x: x[0])

distinct = []
for V_val, p_val, label in all_combined:
    is_new = True
    for V_d, _, _ in distinct:
        if abs(V_val - V_d) < 1e-6 * max(abs(V_d), 1.0):
            is_new = False
            break
    if is_new:
        distinct.append((V_val, p_val, label))

print(f"{len(distinct)} distinct stationary points:")
print(f"  {'#':>3} {'V (MeV^2)':>16} {'theta_12':>10} {'theta_23':>10} {'theta_13':>10} {'offdiag/diag':>14} {'label'}")
for i, (V_val, p_val, label) in enumerate(distinct[:20]):
    V_c, ic = scalar_potential_detailed(p_val)
    Mc = ic['M']
    Uc, sc, Vtc = svd(Mc)
    Vc = Uc.conj().T @ Vtc.conj().T
    if np.linalg.det(Vc).real < 0:
        Vc[:, -1] *= -1
    s13c = abs(Vc[0, 2])
    t13c = np.degrees(np.arcsin(np.clip(s13c, 0, 1)))
    c13c = np.cos(np.radians(t13c))
    s12c = abs(Vc[0, 1]) / c13c if c13c > 1e-10 else 0
    s23c = abs(Vc[1, 2]) / c13c if c13c > 1e-10 else 0
    t12c = np.degrees(np.arcsin(np.clip(s12c, 0, 1)))
    t23c = np.degrees(np.arcsin(np.clip(s23c, 0, 1)))
    od = sum(abs(Mc[a, b]) for a in range(3) for b in range(3) if a != b)
    dg = sum(abs(Mc[a, a]) for a in range(3))
    print(f"  {i+1:3d} {V_c:16.6e} {t12c:10.2f} {t23c:10.2f} {t13c:10.3f} {od/dg:14.4f}  {label}")
print()
sys.stdout.flush()


# =====================================================================
# GRAND SUMMARY
# =====================================================================
print("=" * 78)
print("  GRAND SUMMARY")
print("=" * 78)
print()

M_final = info_best['M']
FM_f = info_best['FM']

print(f"Starting V = {V_start:.6e},  Best V = {V_best_final:.10e}")
print(f"Reduction: {V_start/V_best_final:.2f}x")
print()

print("|M^i_j| at minimum (MeV):")
for i in range(3):
    print(f"  {labels_f[i]:>3}" + "".join(f"  {abs(M_final[i,j]):14.4f}" for j in range(3)))
print()

print("Potential decomposition:")
V_FM_f = np.sum(np.abs(FM_f)**2).real
V_FX_f = abs(info_best['FX'])**2
V_FH_f = abs(info_best['F_Hu0'])**2 + abs(info_best['F_Hd0'])**2
V_soft_f = m_tilde_sq * np.sum(np.abs(M_final)**2).real
print(f"  |F_M|^2  = {V_FM_f:.6e}  ({V_FM_f/V_best_final*100:.1f}%)")
print(f"  |F_X|^2  = {V_FX_f:.6e}  ({V_FX_f/V_best_final*100:.1f}%)")
print(f"  |F_H|^2  = {V_FH_f:.6e}  ({V_FH_f/V_best_final*100:.1f}%)")
print(f"  soft      = {V_soft_f:.6e}  ({V_soft_f/V_best_final*100:.1f}%)")
print()

print(f"X = {info_best['X']:.6e}")
print(f"F_X = {info_best['FX']:.6e}")
print(f"F_Hu0 = {info_best['F_Hu0']:.6e}")
print(f"F_Hd0 = {info_best['F_Hd0']:.6e}")
print()

print("CKM mixing angles (degrees):")
print(f"  {'':>12} {'This':>10} {'PDG':>10} {'Oakes':>10}")
print(f"  {'theta_12':>12} {np.degrees(theta_12):10.4f} {13.04:10.2f} {theta_C_oakes:10.2f}")
print(f"  {'theta_23':>12} {np.degrees(theta_23):10.4f} {2.38:10.2f}")
print(f"  {'theta_13':>12} {np.degrees(theta_13):10.4f} {0.201:10.3f}")
print()

J_f = abs(np.imag(V_CKM[0,0]*V_CKM[1,1]*np.conj(V_CKM[0,1])*np.conj(V_CKM[1,0])))
print(f"J = {J_f:.6e} (PDG: 3.18e-5)")
print()

print(f"Singular values: {sigma_vals[0]:.2f}, {sigma_vals[1]:.2f}, {sigma_vals[2]:.2f}")
print(f"Seesaw values:   {M_diag[0]:.2f}, {M_diag[1]:.2f}, {M_diag[2]:.2f}")
print()

print(f"Hessian: {n_neg_phys} tachyonic, {n_zero_phys} flat, {20-n_neg_phys-n_zero_phys} positive")
if n_neg_phys == 0:
    print("=> TRUE LOCAL MINIMUM")
else:
    print(f"=> SADDLE POINT ({n_neg_phys} tachyonic)")
print()

print(f"Oakes relation: theta_C = arctan(sqrt(m_d/m_s)) = {theta_C_oakes:.4f} deg (PDG: 13.04 deg)")
print()
print("=" * 78)
print("  END")
print("=" * 78)
