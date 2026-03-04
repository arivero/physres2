#!/usr/bin/env python3
"""
Full numerical minimization of the SQCD + Higgs scalar potential.

Superpotential:
    W = sum_i m_i M^i_i + X(det M - B B~ - Lambda^6)
        + y_c H_u^0 M^d_d + y_b H_d^0 M^s_s
        + lambda X (H_u^0 H_d^0 - H_u^+ H_d^-)

Scalar potential:
    V = sum_I |F_I|^2 + f_pi^2 Tr(M^dag M)

Approach:
  (a) 20 real variables: 9 complex meson entries + complex X, Higgs fixed at v/sqrt(2), B = B~ = 0
  (b) Multiple optimizers (L-BFGS-B, Powell, Nelder-Mead) + basin hopping + differential evolution
      with 40 random restarts at perturbation scales 1%, 5%, 10%, 50%
  (c) SVD of meson matrix -> CKM angles
  (d) Compare to PDG
  (e) Hessian eigenvalue check
  (f) Include diagonal imaginary modes explicitly; check tachyonic direction
"""

import numpy as np
from scipy.optimize import minimize, differential_evolution, dual_annealing
from scipy.linalg import svd
import sys
import warnings
import time
warnings.filterwarnings('ignore')

# =====================================================================
# Physical inputs (MeV)
# =====================================================================
m_u = 2.16          # MeV
m_d = 4.67          # MeV
m_s = 93.4          # MeV
m_c = 1270.0        # MeV
m_b = 4180.0        # MeV
Lambda = 300.0       # MeV
Lambda6 = Lambda**6
v_ew = 246220.0      # MeV (246.22 GeV)

# Yukawa couplings: m = y v / 2
y_c = 2.0 * m_c / v_ew
y_b = 2.0 * m_b / v_ew
lam = 0.72           # mu-term coupling

# Soft term
f_pi = 92.0         # MeV
m_tilde_sq = f_pi**2  # 8464 MeV^2

# Higgs VEVs (fixed)
Hu0_fixed = v_ew / np.sqrt(2)
Hd0_fixed = v_ew / np.sqrt(2)

# =====================================================================
# Diagonal (Seiberg seesaw) vacuum
# =====================================================================
masses_q = np.array([m_u, m_d, m_s])
C = Lambda**2 * np.prod(masses_q)**(1.0/3.0)
M_diag = C / masses_q
X_diag = -m_u / (M_diag[1] * M_diag[2])
det_diag = np.prod(M_diag)

# F_X at the seesaw vacuum with Higgs on:
# F_X = det M - Lambda^6 + lam * Hu0 * Hd0 = 0 + lam * v^2/2
FX_seesaw = lam * Hu0_fixed * Hd0_fixed  # = lam * v^2/2 = 2.18e10

print("=" * 78)
print("  FULL CKM MINIMIZATION: SQCD + HIGGS + MU-TERM")
print("=" * 78)
print()
print("Physical inputs:")
print(f"  m_u = {m_u} MeV,  m_d = {m_d} MeV,  m_s = {m_s} MeV")
print(f"  m_c = {m_c} MeV,  m_b = {m_b} MeV")
print(f"  Lambda = {Lambda} MeV,  Lambda^6 = {Lambda6:.6e} MeV^6")
print(f"  v_ew = {v_ew} MeV,  v/sqrt(2) = {v_ew/np.sqrt(2):.4f} MeV")
print(f"  y_c = {y_c:.8e},  y_b = {y_b:.8e},  lambda = {lam}")
print(f"  f_pi = {f_pi} MeV,  m_tilde^2 = {m_tilde_sq} MeV^2")
print()
print("Diagonal vacuum (Seiberg seesaw):")
print(f"  C = Lambda^2 (m_u m_d m_s)^(1/3) = {C:.6f} MeV^2")
for i, fl in enumerate(['u', 'd', 's']):
    print(f"  M_{fl} = C/m_{fl} = {M_diag[i]:.6f} MeV")
print(f"  X_0 = -m_u/(M_d M_s) = {X_diag:.6e}")
print(f"  det M / Lambda^6 = {det_diag / Lambda6:.15f}")
print(f"  F_X(seesaw, Higgs on) = lam*v^2/2 = {FX_seesaw:.6e} MeV^5")
print()
sys.stdout.flush()


# =====================================================================
# Variable packing: 20 real variables
# =====================================================================
# Layout (B = B~ = 0, Higgs fixed):
#   [0..17] : M^i_j complex (9 entries), Re/Im interleaved:
#             M[i,j] -> p[2*(3*i+j)], p[2*(3*i+j)+1]
#   [18,19] : X (Re, Im)

# Scaling: we use physical units directly (no rescaling) for clarity,
# but we do normalize X by a reference scale to improve conditioning.
X_scale = abs(X_diag) if abs(X_diag) > 1e-20 else 1e-9
M_scale = np.zeros((3, 3))
for i in range(3):
    for j in range(3):
        if i == j:
            M_scale[i, j] = max(M_diag[i], 1.0)
        else:
            M_scale[i, j] = max(np.sqrt(M_diag[i] * M_diag[j]), 1.0)

SCALES = np.zeros(20)
for i in range(3):
    for j in range(3):
        idx = 2 * (3*i + j)
        SCALES[idx] = M_scale[i, j]
        SCALES[idx+1] = M_scale[i, j]
SCALES[18] = X_scale
SCALES[19] = X_scale


def unpack(params):
    """Unpack 20 scaled real params -> complex M(3x3), complex X."""
    p = params * SCALES
    M = np.zeros((3, 3), dtype=complex)
    for i in range(3):
        for j in range(3):
            idx = 2 * (3*i + j)
            M[i, j] = p[idx] + 1j * p[idx+1]
    X = p[18] + 1j * p[19]
    return M, X


def cofactor_matrix(M):
    """Cofactor matrix C where C[j,i] = cofactor of M[i,j]."""
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
    """3x3 determinant."""
    return (M[0, 0]*(M[1, 1]*M[2, 2] - M[1, 2]*M[2, 1])
          - M[0, 1]*(M[1, 0]*M[2, 2] - M[1, 2]*M[2, 0])
          + M[0, 2]*(M[1, 0]*M[2, 1] - M[1, 1]*M[2, 0]))


def scalar_potential(params):
    """V = sum |F_I|^2 + f_pi^2 Tr(M^dag M), with B=B~=0, Higgs fixed."""
    M, X = unpack(params)

    cof = cofactor_matrix(M)
    detM = det3(M)
    mass_vec = np.array([m_u, m_d, m_s])

    V = 0.0

    # F_{M^i_j} = m_i delta_{ij} + X * cof[j,i]
    #           + y_c * Hu0 * delta_{i,1} delta_{j,1}   (d-d component, i=1 means d)
    #           + y_b * Hd0 * delta_{i,2} delta_{j,2}   (s-s component, i=2 means s)
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

    # F_X = det M - Lambda^6 + lam * (Hu0 * Hd0)   [B=B~=0, H+=H-=0]
    FX = detM - Lambda6 + lam * Hu0_fixed * Hd0_fixed
    V += abs(FX)**2

    # F_{Hu0} = y_c * M^d_d + lam * X * Hd0
    F_Hu0 = y_c * M[1, 1] + lam * X * Hd0_fixed
    V += abs(F_Hu0)**2

    # F_{Hd0} = y_b * M^s_s + lam * X * Hu0
    F_Hd0 = y_b * M[2, 2] + lam * X * Hu0_fixed
    V += abs(F_Hd0)**2

    # Soft term
    V += m_tilde_sq * np.sum(np.abs(M)**2).real

    return V.real


def scalar_potential_detailed(params):
    """Return V and dictionary of all F-terms and fields."""
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
# (a) Starting point: diagonal seesaw with Higgs on
# =====================================================================
print("=" * 78)
print("  (a) STARTING POINT: SEESAW VACUUM WITH HIGGS FIXED AT v/sqrt(2)")
print("=" * 78)
print()

p0 = np.zeros(20)
# Diagonal mesons
for i in range(3):
    idx = 2 * (3*i + i)
    p0[idx] = M_diag[i] / SCALES[idx]
# X
p0[18] = X_diag / SCALES[18]

V_start, info_start = scalar_potential_detailed(p0)

print(f"V at starting point = {V_start:.6e} MeV^2")
print()

# Decompose
FM0 = info_start['FM']
V_FM = np.sum(np.abs(FM0)**2).real
V_FX = abs(info_start['FX'])**2
V_FH = abs(info_start['F_Hu0'])**2 + abs(info_start['F_Hd0'])**2
V_soft = m_tilde_sq * np.sum(np.abs(info_start['M'])**2).real

print("Potential decomposition:")
print(f"  sum |F_M|^2       = {V_FM:.6e}")
print(f"  |F_X|^2           = {V_FX:.6e}  (dominant: F_X = lam*v^2/2)")
print(f"  sum |F_H|^2       = {V_FH:.6e}")
print(f"  f_pi^2 Tr(M^dag M)= {V_soft:.6e}")
print(f"  Total              = {V_FM + V_FX + V_FH + V_soft:.6e}")
print()

# Individual F-terms
labels_f = ['u', 'd', 's']
print("F-terms at starting point:")
for i in range(3):
    for j in range(3):
        val = FM0[i, j]
        if abs(val) > 1e-6:
            print(f"  F_{{M^{labels_f[i]}_{labels_f[j]}}} = {val:.6e}")
print(f"  F_X   = {info_start['FX']:.6e}")
print(f"  F_Hu0 = {info_start['F_Hu0']:.6e}")
print(f"  F_Hd0 = {info_start['F_Hd0']:.6e}")
print()
sys.stdout.flush()


# =====================================================================
# (b) Multi-strategy minimization
# =====================================================================
print("=" * 78)
print("  (b) MULTI-STRATEGY MINIMIZATION (20 real variables)")
print("=" * 78)
print()
sys.stdout.flush()

t_start = time.time()


def minimize_chain(func, x0, label=""):
    """Sequential: L-BFGS-B -> Powell -> Nelder-Mead."""
    best_x = x0.copy()
    best_v = func(x0)

    # L-BFGS-B
    try:
        r = minimize(func, best_x, method='L-BFGS-B',
                     options={'maxiter': 50000, 'maxfun': 300000, 'ftol': 1e-25})
        if r.fun < best_v:
            best_x, best_v = r.x.copy(), r.fun
    except Exception:
        pass

    # Powell
    try:
        r = minimize(func, best_x, method='Powell',
                     options={'maxiter': 100000, 'maxfev': 500000,
                              'ftol': 1e-25, 'xtol': 1e-14})
        if r.fun < best_v:
            best_x, best_v = r.x.copy(), r.fun
    except Exception:
        pass

    # Nelder-Mead
    try:
        r = minimize(func, best_x, method='Nelder-Mead',
                     options={'maxiter': 300000, 'maxfev': 600000,
                              'xatol': 1e-15, 'fatol': 1e-25, 'adaptive': True})
        if r.fun < best_v:
            best_x, best_v = r.x.copy(), r.fun
    except Exception:
        pass

    return best_v, best_x


# --- From seesaw starting point ---
print("Minimizing from seesaw starting point...")
sys.stdout.flush()
V_seesaw, p_seesaw = minimize_chain(scalar_potential, p0, "seesaw")
print(f"  V = {V_seesaw:.10e}")
print()

all_results = [(V_seesaw, p_seesaw.copy(), 'seesaw')]

# --- Random restarts ---
np.random.seed(42)
n_restarts = 40
perturbation_scales = [0.01, 0.05, 0.10, 0.50]  # 1%, 5%, 10%, 50%

print(f"Running {n_restarts} random restarts with scales {perturbation_scales}...")
sys.stdout.flush()

for trial in range(n_restarts):
    scale = perturbation_scales[trial % len(perturbation_scales)]
    p_trial = p0.copy()

    # Perturb off-diagonal mesons (these are zero at the seesaw, so use absolute perturbation)
    for i in range(3):
        for j in range(3):
            if i != j:
                idx = 2 * (3*i + j)
                # Scale relative to geometric mean of neighboring diagonal entries
                ref = np.sqrt(M_diag[i] * M_diag[j]) / SCALES[idx]
                p_trial[idx] = np.random.randn() * scale * ref
                p_trial[idx+1] = np.random.randn() * scale * ref * 0.3

    # Perturb diagonal mesons (relative)
    for i in range(3):
        idx = 2 * (3*i + i)
        p_trial[idx] *= (1 + np.random.randn() * scale)
        p_trial[idx+1] = np.random.randn() * scale * 0.1 * abs(p_trial[idx])

    # Perturb X (relative)
    p_trial[18] *= (1 + np.random.randn() * scale)
    p_trial[19] = np.random.randn() * scale * 0.1 * abs(p_trial[18])

    V_trial, p_trial_min = minimize_chain(scalar_potential, p_trial, f'restart_{trial}')
    all_results.append((V_trial, p_trial_min.copy(), f'restart_{trial}_s{scale}'))

    if (trial + 1) % 10 == 0:
        print(f"  {trial+1}/{n_restarts} done, best so far: {min(r[0] for r in all_results):.6e}")
        sys.stdout.flush()

# --- Basin hopping ---
print("\nRunning basin hopping (20 iterations)...")
sys.stdout.flush()
try:
    from scipy.optimize import basinhopping

    class StepTaker:
        def __init__(self, stepsize=0.1):
            self.stepsize = stepsize
        def __call__(self, x):
            x_new = x.copy()
            # Larger steps for off-diagonal, smaller for diagonal
            for i in range(3):
                for j in range(3):
                    idx = 2*(3*i+j)
                    if i != j:
                        x_new[idx] += np.random.randn() * self.stepsize
                        x_new[idx+1] += np.random.randn() * self.stepsize * 0.3
                    else:
                        x_new[idx] *= (1 + np.random.randn() * self.stepsize * 0.1)
                        x_new[idx+1] += np.random.randn() * self.stepsize * 0.01 * abs(x_new[idx])
            x_new[18] *= (1 + np.random.randn() * self.stepsize * 0.1)
            x_new[19] += np.random.randn() * self.stepsize * 0.01 * abs(x_new[18])
            return x_new

    bh_kwargs = {'method': 'L-BFGS-B', 'options': {'maxiter': 20000, 'ftol': 1e-22}}
    bh_result = basinhopping(scalar_potential, p0, niter=20,
                              minimizer_kwargs=bh_kwargs,
                              take_step=StepTaker(0.2),
                              seed=123)
    all_results.append((bh_result.fun, bh_result.x.copy(), 'basin_hopping'))
    print(f"  Basin hopping: V = {bh_result.fun:.6e}")
except Exception as e:
    print(f"  Basin hopping failed: {e}")

# --- Differential evolution (global) ---
print("\nRunning differential evolution (bounded search)...")
sys.stdout.flush()
try:
    # Set bounds: diagonal mesons within factor 10 of seesaw, off-diag within diagonal scale, X within factor 10
    bounds_de = []
    for i in range(3):
        for j in range(3):
            idx = 2*(3*i+j)
            if i == j:
                center = M_diag[i] / SCALES[idx]
                width = abs(center) * 2.0
                bounds_de.append((center - width, center + width))  # Re
                bounds_de.append((-width * 0.3, width * 0.3))       # Im
            else:
                ref = np.sqrt(M_diag[i]*M_diag[j]) / SCALES[idx]
                bounds_de.append((-ref, ref))     # Re
                bounds_de.append((-ref*0.3, ref*0.3))  # Im
    center_X = X_diag / SCALES[18]
    width_X = abs(center_X) * 5.0
    bounds_de.append((center_X - width_X, center_X + width_X))
    bounds_de.append((-width_X * 0.3, width_X * 0.3))

    de_result = differential_evolution(scalar_potential, bounds_de,
                                        maxiter=300, popsize=30, tol=1e-18,
                                        seed=456, polish=True)
    all_results.append((de_result.fun, de_result.x.copy(), 'diff_evolution'))
    print(f"  Differential evolution: V = {de_result.fun:.6e}")
except Exception as e:
    print(f"  Differential evolution failed: {e}")

# Sort all results
all_results.sort(key=lambda x: x[0])

print(f"\nBest 10 results:")
for i, (V_val, _, label) in enumerate(all_results[:10]):
    print(f"  {i+1}. V = {V_val:.10e}  ({label})")
print()

# Refine the best result
print("Refining best result (5 rounds)...")
sys.stdout.flush()
p_best = all_results[0][1].copy()
for _ in range(5):
    V_best, p_best = minimize_chain(scalar_potential, p_best, 'refine')
    # Also try with very tight tolerances
    try:
        r = minimize(scalar_potential, p_best, method='Powell',
                     options={'maxiter': 500000, 'maxfev': 1000000,
                              'ftol': 1e-30, 'xtol': 1e-16})
        if r.fun < V_best:
            p_best, V_best = r.x.copy(), r.fun
    except Exception:
        pass

V_best_final, info_best = scalar_potential_detailed(p_best)
t_opt = time.time() - t_start
print(f"  V (refined) = {V_best_final:.10e} MeV^2")
print(f"  Optimization time: {t_opt:.1f} s")
print()
sys.stdout.flush()


# =====================================================================
# (c) CKM mixing angles from SVD
# =====================================================================
print("=" * 78)
print("  (c) CKM MIXING ANGLES FROM MESON MATRIX SVD")
print("=" * 78)
print()

M_min = info_best['M']

print("Meson matrix M^i_j at minimum:")
print(f"  {'':>6}", end="")
for j in range(3):
    print(f"  {labels_f[j]+' Re':>14} {labels_f[j]+' Im':>14}", end="")
print()
for i in range(3):
    print(f"  {labels_f[i]:>6}", end="")
    for j in range(3):
        print(f"  {M_min[i,j].real:14.6f} {M_min[i,j].imag:14.6f}", end="")
    print()
print()

# Magnitudes
print("Meson matrix |M^i_j| (MeV):")
for i in range(3):
    row_str = f"  {labels_f[i]:>3}  "
    for j in range(3):
        row_str += f"  {abs(M_min[i,j]):14.6f}"
    print(row_str)
print()

# SVD: M = U_L @ diag(sigma) @ V_R^dag
U_L, sigma_vals, Vt_R = svd(M_min)
V_R = Vt_R.conj().T

print("SVD: M = U_L diag(sigma) V_R^dag")
print(f"  Singular values: {sigma_vals[0]:.6f}, {sigma_vals[1]:.6f}, {sigma_vals[2]:.6f}")
print(f"  Seesaw values:   {M_diag[0]:.6f}, {M_diag[1]:.6f}, {M_diag[2]:.6f}")
for i in range(3):
    ratio = sigma_vals[i] / M_diag[i]
    print(f"    sigma_{i+1}/M_diag_{labels_f[i]} = {ratio:.8f}")
print()

# CKM-like matrix: V_CKM = U_L^dag V_R
V_CKM = U_L.conj().T @ V_R

# Ensure det ~ +1
det_V = np.linalg.det(V_CKM)
if det_V.real < 0:
    V_CKM[:, -1] *= -1
    det_V = np.linalg.det(V_CKM)

print(f"|V_CKM| (absolute values):")
for i in range(3):
    row_str = "  ["
    for j in range(3):
        row_str += f"  {abs(V_CKM[i,j]):10.8f}"
    row_str += "  ]"
    print(row_str)
print(f"  det(V_CKM) = {det_V:.6f}")
print()

# Standard CKM parametrization angles
s13 = abs(V_CKM[0, 2])
theta_13 = np.arcsin(np.clip(s13, 0, 1))
c13 = np.cos(theta_13)

if c13 > 1e-10:
    s12 = abs(V_CKM[0, 1]) / c13
    s23 = abs(V_CKM[1, 2]) / c13
else:
    s12 = 0
    s23 = 0

theta_12 = np.arcsin(np.clip(s12, 0, 1))
theta_23 = np.arcsin(np.clip(s23, 0, 1))

# (d) PDG comparison
theta_12_pdg = 13.04
theta_23_pdg = 2.38
theta_13_pdg = 0.201
theta_C_oakes = np.degrees(np.arctan(np.sqrt(m_d / m_s)))

J_CKM = abs(np.imag(V_CKM[0,0]*V_CKM[1,1]*np.conj(V_CKM[0,1])*np.conj(V_CKM[1,0])))
J_pdg = 3.18e-5

print("CKM-like mixing angles (standard parametrization):")
print(f"  {'':>12} {'This work':>12} {'PDG':>12} {'Oakes':>12}")
print(f"  {'theta_12':>12} {np.degrees(theta_12):12.4f} {theta_12_pdg:12.2f} {theta_C_oakes:12.2f}")
print(f"  {'theta_23':>12} {np.degrees(theta_23):12.4f} {theta_23_pdg:12.2f} {'---':>12}")
print(f"  {'theta_13':>12} {np.degrees(theta_13):12.4f} {theta_13_pdg:12.3f} {'---':>12}")
print(f"  {'J':>12} {J_CKM:12.4e} {J_pdg:12.2e} {'---':>12}")
print()

# Quark masses from seesaw inversion
print("Quark masses from singular values (seesaw: m_q = C/sigma_q):")
for i in range(3):
    m_seesaw = C / sigma_vals[i] if sigma_vals[i] > 1e-10 else float('inf')
    print(f"  sigma_{i+1} = {sigma_vals[i]:.4f} -> m = {m_seesaw:.4f} MeV  (input m_{labels_f[i]} = {masses_q[i]:.2f})")
print()

# Yukawa masses
m_c_Yukawa = y_c * abs(M_min[1, 1])
m_b_Yukawa = y_b * abs(M_min[2, 2])
print(f"Yukawa masses at minimum:")
print(f"  y_c * |M^d_d| = {m_c_Yukawa:.4f} MeV  (input m_c = {m_c} MeV)")
print(f"  y_b * |M^s_s| = {m_b_Yukawa:.4f} MeV  (input m_b = {m_b} MeV)")
print()
sys.stdout.flush()


# =====================================================================
# (e) Hessian eigenvalues
# =====================================================================
print("=" * 78)
print("  (e) HESSIAN EIGENVALUES AT THE MINIMUM")
print("=" * 78)
print()
sys.stdout.flush()


def numerical_hessian(func, x0, h=1e-5):
    """Central-difference Hessian with adaptive step size."""
    n = len(x0)
    H = np.zeros((n, n))
    f0 = func(x0)
    for i in range(n):
        for j in range(i, n):
            xpp = x0.copy(); xpm = x0.copy()
            xmp = x0.copy(); xmm = x0.copy()
            hi = h * max(abs(x0[i]), 1.0)
            hj = h * max(abs(x0[j]), 1.0)
            xpp[i] += hi; xpp[j] += hj
            xpm[i] += hi; xpm[j] -= hj
            xmp[i] -= hi; xmp[j] += hj
            xmm[i] -= hi; xmm[j] -= hj
            H[i, j] = (func(xpp) - func(xpm) - func(xmp) + func(xmm)) / (4*hi*hj)
            H[j, i] = H[i, j]
    return H


print("Computing 20x20 Hessian (central differences)...")
sys.stdout.flush()
H_scaled = numerical_hessian(scalar_potential, p_best, h=1e-5)
eigvals_scaled = np.linalg.eigvalsh(H_scaled)
eigvals_sorted = np.sort(eigvals_scaled)

# Physical Hessian: convert from scaled to physical variables
H_phys = H_scaled / np.outer(SCALES, SCALES)
eigvals_phys, eigvecs_phys = np.linalg.eigh(H_phys)
idx_sort = np.argsort(eigvals_phys)
eigvals_phys = eigvals_phys[idx_sort]
eigvecs_phys = eigvecs_phys[:, idx_sort]

print()
print("Scaled Hessian eigenvalues:")
threshold = 1e-4 * abs(eigvals_sorted[-1]) if abs(eigvals_sorted[-1]) > 0 else 1e-6
n_neg = 0
n_zero = 0
for i, ev in enumerate(eigvals_sorted):
    if ev < -threshold:
        status = "TACHYONIC"
        n_neg += 1
    elif abs(ev) < threshold:
        status = "~ZERO (flat)"
        n_zero += 1
    else:
        status = "OK"
    print(f"  lambda_{i+1:2d} = {ev:14.6e}  {status}")

print()
print(f"Summary: {n_neg} negative, {n_zero} near-zero, {20 - n_neg - n_zero} positive")
if n_neg == 0:
    print("=> TRUE LOCAL MINIMUM (no tachyonic directions)")
else:
    print(f"=> SADDLE POINT with {n_neg} tachyonic direction(s)")
print()

# Physical eigenvalues (mass^2 in appropriate units)
print("Physical Hessian eigenvalues (d^2V/dphi_i dphi_j):")
for i in range(20):
    ev = eigvals_phys[i]
    if ev < -threshold / SCALES[0]**2:
        status = "TACHYONIC"
    elif abs(ev) < threshold / SCALES[0]**2:
        status = "~ZERO"
    else:
        status = "OK"
    print(f"  m^2_{i+1:2d} = {ev:14.6e}  {status}")
print()
sys.stdout.flush()


# =====================================================================
# (f) Tachyonic direction analysis
# =====================================================================
print("=" * 78)
print("  (f) TACHYONIC DIRECTION ANALYSIS")
print("=" * 78)
print()

# Analyze the tachyonic eigenvectors
variable_names = []
for i in range(3):
    for j in range(3):
        variable_names.append(f"Re(M^{labels_f[i]}_{labels_f[j]})")
        variable_names.append(f"Im(M^{labels_f[i]}_{labels_f[j]})")
variable_names.append("Re(X)")
variable_names.append("Im(X)")

tachyonic_count = 0
for i in range(20):
    ev = eigvals_phys[i]
    if ev < -threshold / SCALES[0]**2:
        tachyonic_count += 1
        vec = eigvecs_phys[:, i]
        print(f"Tachyonic mode {tachyonic_count}: m^2 = {ev:.6e}")
        print(f"  Eigenvector components (|c| > 0.1):")
        for k in range(20):
            if abs(vec[k]) > 0.1:
                print(f"    {variable_names[k]:>20}: {vec[k]:+.6f}")
        print()

if tachyonic_count == 0:
    # Check for near-zero modes
    print("No tachyonic modes found. Checking near-zero modes:")
    for i in range(min(5, 20)):
        ev = eigvals_phys[i]
        vec = eigvecs_phys[:, i]
        print(f"  Mode {i+1}: m^2 = {ev:.6e}")
        print(f"    Dominant components:")
        sorted_idx = np.argsort(-np.abs(vec))
        for k in sorted_idx[:4]:
            print(f"      {variable_names[k]:>20}: {vec[k]:+.6f}")
        print()

# Check the specific tachyonic direction reported: Im(0.87 M^d_d + 0.49 M^u_u)
print("Checking reported tachyonic direction Im(0.87 M^d_d + 0.49 M^u_u):")
# Build direction vector in scaled space
direction = np.zeros(20)
# Im(M^d_d) is at index 2*(3*1+1)+1 = 7
# Im(M^u_u) is at index 2*(3*0+0)+1 = 1
direction[7] = 0.87 * SCALES[7]   # Scale to physical
direction[1] = 0.49 * SCALES[1]
direction /= np.linalg.norm(direction)

# Curvature along this direction in scaled space
m2_diag_im = direction @ H_scaled @ direction
print(f"  Curvature along this direction (scaled): {m2_diag_im:.6e}")
# In physical space
direction_phys = np.zeros(20)
direction_phys[7] = 0.87
direction_phys[1] = 0.49
direction_phys /= np.linalg.norm(direction_phys)
m2_diag_im_phys = direction_phys @ H_phys @ direction_phys
print(f"  Curvature along this direction (physical): {m2_diag_im_phys:.6e}")
print()
sys.stdout.flush()


# =====================================================================
# Additional: explore along tachyonic directions
# =====================================================================
print("=" * 78)
print("  EXPLORING TACHYONIC DIRECTIONS")
print("=" * 78)
print()

# If there are tachyonic directions, try minimizing along them
if n_neg > 0 or True:
    # Try displacing along the most negative scaled eigenvector and reminimizing
    _, eigvecs_scaled = np.linalg.eigh(H_scaled)

    displacements_to_try = [0.001, 0.01, 0.1, 0.5, 1.0, 2.0]

    results_displaced = []

    # For each of the most negative eigenvalue directions
    n_try = min(5, 20)
    for mode_idx in range(n_try):
        evec = eigvecs_scaled[:, mode_idx]
        ev = eigvals_sorted[mode_idx]

        for disp in displacements_to_try:
            p_displaced = p_best + disp * evec
            V_disp_start = scalar_potential(p_displaced)
            V_disp, p_disp = minimize_chain(scalar_potential, p_displaced, f'mode{mode_idx}_d{disp}')
            results_displaced.append((V_disp, p_disp.copy(), f'mode_{mode_idx}_disp_{disp}', ev))

            # Also try negative displacement
            p_displaced_neg = p_best - disp * evec
            V_disp_neg, p_disp_neg = minimize_chain(scalar_potential, p_displaced_neg, f'mode{mode_idx}_d-{disp}')
            results_displaced.append((V_disp_neg, p_disp_neg.copy(), f'mode_{mode_idx}_disp_{-disp}', ev))

    results_displaced.sort(key=lambda x: x[0])

    print(f"Best 10 results from tachyonic displacement search:")
    for i, (V_val, _, label, ev) in enumerate(results_displaced[:10]):
        marker = " ***" if V_val < V_best_final * 0.999 else ""
        print(f"  {i+1}. V = {V_val:.10e}  ({label}, eigenvalue={ev:.2e}){marker}")
    print()

    # If we found something better, update
    if results_displaced[0][0] < V_best_final * (1 - 1e-10):
        print(f"IMPROVED: V went from {V_best_final:.10e} to {results_displaced[0][0]:.10e}")
        V_best_final_new = results_displaced[0][0]
        p_best_new = results_displaced[0][1].copy()

        # Refine again
        for _ in range(3):
            V_best_final_new, p_best_new = minimize_chain(scalar_potential, p_best_new, 'refine2')

        V_best_final2, info_best2 = scalar_potential_detailed(p_best_new)

        # Recompute CKM
        M_min2 = info_best2['M']
        U_L2, sigma2, Vt_R2 = svd(M_min2)
        V_R2 = Vt_R2.conj().T
        V_CKM2 = U_L2.conj().T @ V_R2
        if np.linalg.det(V_CKM2).real < 0:
            V_CKM2[:, -1] *= -1

        s13_2 = abs(V_CKM2[0, 2])
        t13_2 = np.arcsin(np.clip(s13_2, 0, 1))
        c13_2 = np.cos(t13_2)
        s12_2 = abs(V_CKM2[0, 1]) / c13_2 if c13_2 > 1e-10 else 0
        s23_2 = abs(V_CKM2[1, 2]) / c13_2 if c13_2 > 1e-10 else 0
        t12_2 = np.arcsin(np.clip(s12_2, 0, 1))
        t23_2 = np.arcsin(np.clip(s23_2, 0, 1))

        print(f"\nImproved minimum: V = {V_best_final2:.10e}")
        print(f"\n|V_CKM| at improved minimum:")
        for i in range(3):
            row_str = "  ["
            for j in range(3):
                row_str += f"  {abs(V_CKM2[i,j]):10.8f}"
            row_str += "  ]"
            print(row_str)

        print(f"\nCKM angles at improved minimum:")
        print(f"  theta_12 = {np.degrees(t12_2):.4f} deg  (PDG: {theta_12_pdg})")
        print(f"  theta_23 = {np.degrees(t23_2):.4f} deg  (PDG: {theta_23_pdg})")
        print(f"  theta_13 = {np.degrees(t13_2):.4f} deg  (PDG: {theta_13_pdg})")

        # Update best
        V_best_final = V_best_final2
        p_best = p_best_new
        info_best = info_best2
        M_min = M_min2
        theta_12 = t12_2
        theta_23 = t23_2
        theta_13 = t13_2
        sigma_vals = sigma2
        V_CKM = V_CKM2

        # Recompute Hessian
        print("\nRecomputing Hessian at improved minimum...")
        H_scaled = numerical_hessian(scalar_potential, p_best, h=1e-5)
        eigvals_sorted = np.sort(np.linalg.eigvalsh(H_scaled))
        H_phys = H_scaled / np.outer(SCALES, SCALES)
        eigvals_phys, eigvecs_phys = np.linalg.eigh(H_phys)
        idx_sort = np.argsort(eigvals_phys)
        eigvals_phys = eigvals_phys[idx_sort]
        eigvecs_phys = eigvecs_phys[:, idx_sort]

        n_neg = sum(1 for ev in eigvals_sorted if ev < -threshold)
        n_zero = sum(1 for ev in eigvals_sorted if abs(ev) < threshold)

        print("Hessian eigenvalues at improved minimum:")
        for i, ev in enumerate(eigvals_sorted):
            if ev < -threshold:
                status = "TACHYONIC"
            elif abs(ev) < threshold:
                status = "~ZERO"
            else:
                status = "OK"
            print(f"  lambda_{i+1:2d} = {ev:14.6e}  {status}")
        print(f"Summary: {n_neg} tachyonic, {n_zero} flat, {20-n_neg-n_zero} positive")
        print()
    else:
        print("No improvement found from tachyonic displacement search.")
        print()

sys.stdout.flush()


# =====================================================================
# Collect all distinct minima
# =====================================================================
print("=" * 78)
print("  COLLECTING DISTINCT MINIMA")
print("=" * 78)
print()

# Merge all results
all_combined = all_results + results_displaced
all_combined.sort(key=lambda x: x[0])

# Cluster by V value (consider same if within 0.1%)
distinct_minima = []
for V_val, p_val, label, *_ in all_combined:
    is_new = True
    for V_d, _, _ in distinct_minima:
        if abs(V_val - V_d) < 1e-6 * max(abs(V_d), 1.0):
            is_new = False
            break
    if is_new:
        distinct_minima.append((V_val, p_val, label))

print(f"Found {len(distinct_minima)} distinct stationary points:")
for i, (V_val, p_val, label) in enumerate(distinct_minima[:15]):
    V_check, info_check = scalar_potential_detailed(p_val)
    M_check = info_check['M']
    U_c, sig_c, Vt_c = svd(M_check)
    V_ckm_c = U_c.conj().T @ Vt_c.conj().T
    if np.linalg.det(V_ckm_c).real < 0:
        V_ckm_c[:, -1] *= -1
    s13_c = abs(V_ckm_c[0, 2])
    t13_c = np.degrees(np.arcsin(np.clip(s13_c, 0, 1)))
    c13_c = np.cos(np.radians(t13_c))
    s12_c = abs(V_ckm_c[0, 1]) / c13_c if c13_c > 1e-10 else 0
    s23_c = abs(V_ckm_c[1, 2]) / c13_c if c13_c > 1e-10 else 0
    t12_c = np.degrees(np.arcsin(np.clip(s12_c, 0, 1)))
    t23_c = np.degrees(np.arcsin(np.clip(s23_c, 0, 1)))

    off_diag = sum(abs(M_check[i,j]) for i in range(3) for j in range(3) if i != j)
    diag_sum = sum(abs(M_check[i,i]) for i in range(3))

    print(f"  {i+1}. V = {V_check:.6e}  theta12={t12_c:.2f}  theta23={t23_c:.2f}  "
          f"theta13={t13_c:.3f}  offdiag/diag={off_diag/diag_sum:.4f}  ({label})")
print()
sys.stdout.flush()


# =====================================================================
# GRAND SUMMARY
# =====================================================================
print("=" * 78)
print("  GRAND SUMMARY")
print("=" * 78)
print()

print(f"Starting point: V = {V_start:.6e} MeV^2")
print(f"Best minimum:   V = {V_best_final:.10e} MeV^2")
print(f"Reduction factor: {V_start / V_best_final:.4f}")
print()

M_final = info_best['M']
print("Meson matrix at best minimum (|M^i_j|, MeV):")
for i in range(3):
    row_str = f"  {labels_f[i]:>3}  "
    for j in range(3):
        row_str += f"  {abs(M_final[i,j]):14.6f}"
    print(row_str)
print()

print(f"Other fields at minimum:")
print(f"  X    = {info_best['X']:.6e}")
print(f"  detM = {info_best['detM']:.6e}")
print(f"  F_X  = {info_best['FX']:.6e}")
print()

# F-term decomposition
FM_final = info_best['FM']
V_FM_final = np.sum(np.abs(FM_final)**2).real
V_FX_final = abs(info_best['FX'])**2
V_FH_final = abs(info_best['F_Hu0'])**2 + abs(info_best['F_Hd0'])**2
V_soft_final = m_tilde_sq * np.sum(np.abs(M_final)**2).real

print("Potential decomposition at minimum:")
print(f"  sum |F_M|^2       = {V_FM_final:.6e}  ({V_FM_final/V_best_final*100:.1f}%)")
print(f"  |F_X|^2           = {V_FX_final:.6e}  ({V_FX_final/V_best_final*100:.1f}%)")
print(f"  sum |F_H|^2       = {V_FH_final:.6e}  ({V_FH_final/V_best_final*100:.1f}%)")
print(f"  f_pi^2 Tr(M^dag M)= {V_soft_final:.6e}  ({V_soft_final/V_best_final*100:.1f}%)")
print()

print("CKM mixing angles:")
print(f"  {'':>12} {'This work':>12} {'PDG':>12} {'Oakes':>12} {'Ratio':>12}")
print(f"  {'theta_12':>12} {np.degrees(theta_12):12.4f} {theta_12_pdg:12.2f} {theta_C_oakes:12.2f} {np.degrees(theta_12)/theta_12_pdg:12.4f}")
print(f"  {'theta_23':>12} {np.degrees(theta_23):12.4f} {theta_23_pdg:12.2f} {'---':>12} {np.degrees(theta_23)/theta_23_pdg:12.4f}")
print(f"  {'theta_13':>12} {np.degrees(theta_13):12.4f} {theta_13_pdg:12.3f} {'---':>12} {np.degrees(theta_13)/theta_13_pdg:12.4f}")
print()

J_final = abs(np.imag(V_CKM[0,0]*V_CKM[1,1]*np.conj(V_CKM[0,1])*np.conj(V_CKM[1,0])))
print(f"Jarlskog invariant: J = {J_final:.6e}  (PDG: {J_pdg:.2e})")
print()

print(f"Singular values: {sigma_vals[0]:.4f}, {sigma_vals[1]:.4f}, {sigma_vals[2]:.4f}")
print(f"Seesaw values:   {M_diag[0]:.4f}, {M_diag[1]:.4f}, {M_diag[2]:.4f}")
print()

print(f"Hessian: {n_neg} tachyonic, {n_zero} flat, {20-n_neg-n_zero} positive")
if n_neg == 0:
    print(f"  => TRUE LOCAL MINIMUM")
else:
    print(f"  => SADDLE POINT")
print()

# Weinberg-Oakes comparison
print(f"Weinberg-Oakes: theta_C = arctan(sqrt(m_d/m_s)) = {theta_C_oakes:.4f} deg")
print(f"  PDG Cabibbo angle = {theta_12_pdg} deg")
print(f"  This computation  = {np.degrees(theta_12):.4f} deg")
print()

print("=" * 78)
print("  END OF COMPUTATION")
print("=" * 78)
