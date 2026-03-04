#!/usr/bin/env python3
"""
O'Raifeartaigh model with Kähler pole: stabilization analysis.

Superpotential: W = f Φ₀ + m Φ₁ Φ₂ + g Φ₀ Φ₁²
Kähler: K = |Φ₀|² + c|Φ₀|⁴/Λ² + |Φ₁|² + |Φ₂|²

Units: m = 1. Parameters: t = gv/m, y = gf/m².
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize_scalar, brentq
from scipy.integrate import quad

# ============================================================
# Mass spectrum of the O'Raifeartaigh model
# ============================================================
# At vacuum: Φ₁ = Φ₂ = 0, Φ₀ = v (pseudo-modulus)
# t = gv/m, y = gf/m²
#
# Superpotential derivatives:
# W_0 = f + g Φ₁² → f at vacuum
# W_1 = m Φ₂ + 2g Φ₀ Φ₁ → 0
# W_2 = m Φ₁ → 0
#
# F_0 = f (nonzero → SUSY breaking)
# F_1 = F_2 = 0
#
# Fermion mass matrix W_ij:
# W_00 = 0, W_01 = W_10 = 2gΦ₁ = 0, W_02 = W_20 = 0
# W_11 = 2gΦ₀ = 2gv = 2tm, W_12 = W_21 = m, W_22 = 0
#
# Fermion masses: eigenvalues of W_ij
# det(W_ij - λδ_ij) = -λ[(2tm - λ)(-λ) - m²] = -λ[λ² - 2tmλ - m²]
# = -λ(λ - tm - m√(t²+1))(λ - tm + m√(t²+1))
#
# So: λ₀ = 0, λ± = tm ± m√(t²+1) → |λ±| = m(√(t²+1) ± t)
# Fermion masses squared: 0, m²(√(t²+1) - t)², m²(√(t²+1) + t)²
# In m=1 units: 0, (√(t²+1) - t)², (√(t²+1) + t)²
#
# Scalar masses: eigenvalues of the 6×6 matrix
# M²_scalar has eigenvalues from (M²_F ± B):
#
# In the (Φ₁, Φ₂) sector:
# M²_F = W†W restricted:
# = ( |W_11|² + |W_01|² ,  W_11* W_12 + W_01* W_02 )
#   ( W_12* W_11 + W_02* W_01,  |W_12|² + |W_02|² )
# = ( 4t² + 0,  2t·1 )  [in m=1]
#   ( 2t,        1    )
# Wait, W†W:
# (W†)_ij = W_ji*
# (W†W)_ij = Σ_k W_ki* W_kj
# For i,j in {1,2}:
# (W†W)_11 = |W_01|² + |W_11|² + |W_21|² = 0 + 4t² + 1
# (W†W)_12 = W_01*W_02 + W_11*W_12 + W_21*W_22 = 0 + 2t·1 + 1·0 = 2t
# (W†W)_21 = 2t
# (W†W)_22 = |W_02|² + |W_12|² + |W_22|² = 0 + 1 + 0 = 1
#
# B_ij = Σ_k W_ijk F_k*  (only k=0 contributes, F_0* = f)
# W_ij0 = ∂W_ij/∂Φ_0:
# W_110 = ∂(2gΦ₀)/∂Φ₀ = 2g → B_11 = 2g·f = 2y (in m=1 units, y=gf)
# All others zero → B = diag(2y, 0) in the (1,2) sector
#
# Wait, but F_0 = -∂W*/∂Φ₀* which at vacuum gives -f* (assuming f real, F_0 = -f)
# Convention: F_i = -W_i*, so F_0* = -W_0 = -(f) → F_0* = -f
# B_ij = W_ij0 · F_0* = W_ij0 · (-f)
# B_11 = 2g · (-f) = -2gf = -2y·m² = -2y [in m=1 units]
# 
# Scalar mass matrices:
# M²_+ = M²_F + B (real parts)
# M²_- = M²_F - B (imaginary parts)
#
# M²_+ = ( 4t² + 1 - 2y,  2t )
#         ( 2t,             1  )
# M²_- = ( 4t² + 1 + 2y,  2t )
#         ( 2t,             1  )
#
# For Φ₀: tree-level mass = 0 (flat direction, both real components)

def fermion_masses_sq(t):
    """Returns array of 3 fermion mass-squared values (Weyl fermions)."""
    s = np.sqrt(t**2 + 1)
    return np.array([0.0, (s - t)**2, (s + t)**2])

def scalar_masses_sq(t, y):
    """Returns array of 6 real scalar mass-squared values."""
    # Φ₁-Φ₂ sector
    Mplus = np.array([[4*t**2 + 1 - 2*y, 2*t],
                       [2*t,               1  ]])
    Mminus = np.array([[4*t**2 + 1 + 2*y, 2*t],
                        [2*t,               1  ]])
    
    eig_p = np.linalg.eigvalsh(Mplus)
    eig_m = np.linalg.eigvalsh(Mminus)
    
    # Φ₀: m² = 0, 0 at tree level
    return np.concatenate([[0.0, 0.0], eig_p, eig_m])

def CW_contribution(msq_arr, n_dof_arr, mu_sq=1.0):
    """CW potential from an array of masses with given dof (sign: +1 boson, -1 fermion)."""
    V = 0.0
    for msq, ndof in zip(msq_arr, n_dof_arr):
        if msq > 1e-30:
            V += ndof * msq**2 * (np.log(msq / mu_sq) - 1.5)
    return V / (64 * np.pi**2)

def compute_V_CW(t, y, mu_sq=1.0):
    """Full supertrace CW potential."""
    smsq = scalar_masses_sq(t, y)
    fmsq = fermion_masses_sq(t)
    
    # Scalars: each real scalar contributes +1
    # Fermions: each Weyl fermion contributes -2 (2 real d.o.f., fermionic)
    
    all_msq = np.concatenate([smsq, fmsq])
    all_ndof = np.concatenate([np.ones(6), -2*np.ones(3)])
    
    return CW_contribution(all_msq, all_ndof, mu_sq)

# ============================================================
# PART 1: Plot V_CW and find minimum
# ============================================================
print("=" * 60)
print("PART 1: Coleman-Weinberg potential")
print("=" * 60)

y = 0.1
t_arr = np.linspace(0.01, 4.0, 2000)
V_cw = np.array([compute_V_CW(t, y) for t in t_arr])

# Find minimum
res = minimize_scalar(lambda t: compute_V_CW(t, y), bounds=(0.01, 3.0), method='bounded')
t_min_cw = res.x
V_min_cw = res.fun
print(f"  y = {y}")
print(f"  CW minimum at t = {t_min_cw:.6f}")
print(f"  V_CW(t_min) = {V_min_cw:.8e}")
print(f"  V_CW(0.01) = {compute_V_CW(0.01, y):.8e}")
print(f"  V_CW(√3) = {compute_V_CW(np.sqrt(3), y):.8e}")

# Check if minimum is at origin
# V_CW should decrease from origin → the pseudomodulus is pushed to v=0
# or some finite value
print(f"\n  V_CW values:")
for tt in [0.01, 0.1, 0.3, 0.49, 0.5, 1.0, np.sqrt(3), 2.0, 3.0]:
    print(f"    t={tt:.3f}: V_CW = {compute_V_CW(tt, y):.8e}")

# Try different y values
print(f"\n  CW minima for various y:")
for yy in [0.01, 0.05, 0.1, 0.3, 0.5, 0.8, 1.0, 1.5, 2.0]:
    res2 = minimize_scalar(lambda t: compute_V_CW(t, yy), bounds=(0.01, 5.0), method='bounded')
    print(f"    y={yy:.2f}: t_min = {res2.x:.4f}")

fig, axes = plt.subplots(1, 2, figsize=(14, 5))

# Left: V_CW
ax = axes[0]
V_shifted = V_cw - V_cw[0]
ax.plot(t_arr, V_shifted, 'b-', linewidth=2)
ax.axvline(t_min_cw, color='r', linestyle='--', alpha=0.7, label=f'min at t={t_min_cw:.3f}')
ax.axvline(np.sqrt(3), color='g', linestyle='--', alpha=0.7, label=f't = sqrt(3) = {np.sqrt(3):.3f}')
ax.set_xlabel('t = gv/m', fontsize=13)
ax.set_ylabel('V_CW(t) - V_CW(0)', fontsize=13)
ax.set_title(f'CW potential (y = {y})', fontsize=14)
ax.legend(fontsize=11)
ax.grid(True, alpha=0.3)

# Right: V_CW for several y values
ax = axes[1]
for yy in [0.05, 0.1, 0.3, 0.5, 1.0]:
    V_y = np.array([compute_V_CW(t, yy) for t in t_arr])
    V_y -= V_y[0]
    ax.plot(t_arr, V_y, linewidth=1.5, label=f'y={yy}')
ax.axvline(np.sqrt(3), color='g', linestyle='--', alpha=0.5)
ax.set_xlabel('t = gv/m', fontsize=13)
ax.set_ylabel('V_CW(t) - V_CW(0)', fontsize=13)
ax.set_title('CW potential for various y', fontsize=14)
ax.legend(fontsize=11)
ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('/home/codexssh/phys3/results/kahler_pole_CW_part1.png', dpi=150)
plt.close()
print("  Plot saved: results/kahler_pole_CW_part1.png")

# ============================================================
# PART 2: Total effective potential = tree + CW
# ============================================================
print("\n" + "=" * 60)
print("PART 2: Total effective potential with Kähler pole")
print("=" * 60)

def K00(t, c=-1.0/12):
    """Kähler metric K_{0bar0} = 1 + 4c t²  (with Λ = m/g, v = t m/g)
    K = |Φ₀|² + c|Φ₀|⁴/Λ²
    K_{0bar0} = 1 + 4c|Φ₀|²/Λ² = 1 + 4c (gv/m)²·(g/m)²·(m/g)² 
    Wait: |Φ₀|²/Λ² = v²/(m/g)² = (gv/m)² = t²
    So K_{0bar0} = 1 + 4ct² = 1 - t²/3 for c=-1/12
    """
    return 1 + 4*c*t**2

def V_tree(t, y, c=-1.0/12):
    """Tree-level potential V = f²/K_{0bar0} = y²/(g² K_{0bar0})
    In m=1 units: f = y/g, f² = y²/g²
    But we want V in units of m⁴... V = |F₀|² / K_{0bar0}
    |F₀|² = |f|² = (y·m²/g)² / m⁴... 
    
    Actually let's be careful. In m=1 units:
    V_tree = f² / K_{0bar0}  [the F-term potential with non-canonical K]
    f = y·m²/g = y/g (m=1)
    But y = gf/m² = gf, so f = y/g.
    f² = y²/g²
    
    But we haven't fixed g. The CW potential also depends on m=1 only.
    The tree potential is V_tree = f² = y²/g² at the canonical point.
    
    Hmm, this is an issue. The tree-level potential is f², and the CW 
    correction is O(m⁴/(16π²)). For the two to compete, we need f² ~ m⁴/(16π²),
    i.e., y/g ~ m²/(4π) → y ~ g/(4π).
    
    Let me write everything in m=1 units with f = y/g:
    V_tree = f²/K_{0bar0} = (y/g)² / K_{0bar0}
    V_CW ~ m⁴/(64π²) · h(t,y) = h(t,y)/(64π²)
    
    For them to compete: (y/g)² ~ 1/(64π²)
    
    Actually for the pseudomodulus stabilization problem, we typically 
    work with the CW potential alone (the tree-level f² is constant 
    if K is canonical). The tree level potential WITH the Kähler pole
    adds the 1/K_{0bar0} factor.
    
    Let me just set g=1 for simplicity. Then f = y, and:
    V_tree = y² / K_{0bar0}(t)
    """
    k = K00(t, c)
    if k <= 0:
        return np.inf
    return y**2 / k

def V_eff(t, y, c=-1.0/12, mu_sq=1.0):
    """Total effective potential."""
    return V_tree(t, y, c) + compute_V_CW(t, y, mu_sq)

# Find minima for various y
print("\n  Total potential minima (c = -1/12, Λ = m/g = 1 with g=1):")
print(f"  Pole at t = √3 = {np.sqrt(3):.6f}")

t_pole = np.sqrt(3)
t_search = np.linspace(0.01, t_pole - 0.01, 2000)

fig, axes = plt.subplots(2, 2, figsize=(14, 10))

y_values = [0.01, 0.05, 0.1, 0.3]
for idx, yy in enumerate(y_values):
    ax = axes[idx // 2][idx % 2]
    
    V_total = np.array([V_eff(t, yy) for t in t_search])
    V_tree_arr = np.array([V_tree(t, yy) for t in t_search])
    V_cw_arr = np.array([compute_V_CW(t, yy) for t in t_search])
    
    # Shift for visibility
    V_ref = V_total[0]
    
    ax.plot(t_search, V_total - V_ref, 'b-', linewidth=2, label='V_total - V(0)')
    ax.plot(t_search, V_tree_arr - V_tree_arr[0], 'r--', linewidth=1.5, label='V_tree - V_tree(0)')
    ax.plot(t_search, V_cw_arr - V_cw_arr[0], 'g--', linewidth=1.5, label='V_CW - V_CW(0)')
    
    # Find minimum of total
    res3 = minimize_scalar(lambda t: V_eff(t, yy), bounds=(0.01, t_pole - 0.01), method='bounded')
    t_min_eff = res3.x
    
    ax.axvline(t_min_eff, color='b', linestyle=':', alpha=0.7, label=f'min at t={t_min_eff:.3f}')
    ax.axvline(t_pole, color='k', linestyle='--', alpha=0.5, label=f'pole t={t_pole:.3f}')
    
    ax.set_xlabel('t = gv/m', fontsize=12)
    ax.set_ylabel('V - V(0)', fontsize=12)
    ax.set_title(f'y = {yy}', fontsize=13)
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)
    
    print(f"  y = {yy:.3f}: t_min = {t_min_eff:.4f}, ratio t_min/√3 = {t_min_eff/t_pole:.4f}")

plt.suptitle('Total effective potential V_tree + V_CW', fontsize=15, y=1.01)
plt.tight_layout()
plt.savefig('/home/codexssh/phys3/results/kahler_pole_CW_part2.png', dpi=150)
plt.close()
print("  Plot saved: results/kahler_pole_CW_part2.png")

# Scan over y to find if minimum ever reaches √3
print("\n  Scanning y to find when minimum → √3:")
y_scan = np.logspace(-3, 1, 200)
t_mins = []
for yy in y_scan:
    res4 = minimize_scalar(lambda t: V_eff(t, yy), bounds=(0.02, t_pole - 0.005), method='bounded')
    t_mins.append(res4.x)
t_mins = np.array(t_mins)

fig, ax = plt.subplots(figsize=(8, 5))
ax.semilogx(y_scan, t_mins, 'b-', linewidth=2)
ax.axhline(t_pole, color='r', linestyle='--', label=f't = sqrt(3)')
ax.set_xlabel('y = gf/m²', fontsize=13)
ax.set_ylabel('t_min', fontsize=13)
ax.set_title('Position of V_eff minimum vs y', fontsize=14)
ax.legend(fontsize=12)
ax.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig('/home/codexssh/phys3/results/kahler_pole_CW_tmin_vs_y.png', dpi=150)
plt.close()
print("  Plot saved: results/kahler_pole_CW_tmin_vs_y.png")

# ============================================================
# PART 3: Kähler-dressed CW potential
# ============================================================
print("\n" + "=" * 60)
print("PART 3: Kähler-dressed mass matrix in CW potential")
print("=" * 60)

def compute_V_CW_dressed(t, y, c=-1.0/12, mu_sq=1.0):
    """CW potential with Kähler-dressed masses.
    For the Φ₀ sector: masses get divided by K_{0bar0}.
    For Φ₁, Φ₂: K_{1bar1} = K_{2bar2} = 1, no change.
    
    The dressing affects the physical masses that enter loops.
    For Φ₀ (tree-level m² = 0): stays 0.
    
    But the mixing between Φ₀ and Φ₁,Φ₂ in the mass matrix 
    gets dressed by K^{-1/2}. For Φ₀ with K_{0bar0} ≠ 1, the 
    canonical normalization changes.
    
    Actually, the full mass matrix includes Φ₀. Let me construct it properly.
    
    The fermion mass matrix in the canonical basis:
    Ŵ_ij = (K_{i\bar{i}})^{-1/2} W_ij (K_{j\bar{j}})^{-1/2}
    
    W_ij:
    W_00 = 0, W_01 = 0, W_02 = 0
    W_11 = 2t, W_12 = 1, W_22 = 0
    
    Dressed:
    Ŵ_00 = 0 / K_{0bar0} = 0
    Ŵ_01 = 0 / sqrt(K_{0bar0}) = 0  
    Ŵ_11 = 2t  (K_{1bar1} = 1)
    Ŵ_12 = 1
    Ŵ_22 = 0
    
    So fermion masses are UNCHANGED (because W has no dependence on 
    Φ₀ in the mass matrix rows/columns — Φ₀ is massless).
    
    For scalars: the F-term contribution has the dressing.
    The scalar mass matrix involves F_0 = f, and the K-metric affects
    the physical F-term: F̂_0 = F_0 / sqrt(K_{0bar0})
    
    Actually the tree-level scalar potential is:
    V = K^{i\bar{j}} |W_i|² = |W_0|²/K_{0bar0} + |W_1|² + |W_2|²
    
    The second derivative gives the scalar mass matrix.
    For the B-term: B_ij = W_ij0 F̂^{0*} where F̂^0 = F^0/sqrt(K_{0bar0})
    Wait, more carefully:
    F̂^0 = K^{0\bar0} F_0 = F_0/K_{0bar0}  (no, that's the equation of motion)
    
    The physical (canonical) F-term is F̂_0 = F_0/sqrt(K_{0bar0}).
    And the cubic coupling W_{110} in canonical normalization:
    Ŵ_{110} = W_{110}/sqrt(K_{0bar0}) = 2g/sqrt(K_{0bar0})
    
    So the effective B-term in canonical normalization:
    B̂_11 = Ŵ_{110} · F̂^{0*} = (2g/sqrt(K_{0bar0})) · (f/sqrt(K_{0bar0}))
          = 2gf/K_{0bar0} = 2y/K_{0bar0}
    
    So the ONLY effect of the Kähler dressing on the mass matrix is:
    y → y/K_{0bar0} in the B-term splitting.
    The fermion masses are unchanged.
    """
    k = K00(t, c)
    if k <= 0:
        return np.inf
    
    y_eff = y / k  # dressed y
    
    smsq = scalar_masses_sq(t, y_eff)
    fmsq = fermion_masses_sq(t)
    
    all_msq = np.concatenate([smsq, fmsq])
    all_ndof = np.concatenate([np.ones(6), -2*np.ones(3)])
    
    return CW_contribution(all_msq, all_ndof, mu_sq)

# Compare dressed vs undressed
print("\n  Comparing dressed vs undressed CW:")
for yy in [0.1, 0.5, 1.0]:
    for tt in [0.5, 1.0, 1.5]:
        v_und = compute_V_CW(tt, yy)
        v_dre = compute_V_CW_dressed(tt, yy)
        k = K00(tt)
        print(f"    y={yy}, t={tt}: K={k:.4f}, V_CW={v_und:.6e}, V_CW_dressed={v_dre:.6e}, ratio={v_dre/v_und:.4f}")

# Total potential with dressed CW
def V_eff_dressed(t, y, c=-1.0/12, mu_sq=1.0):
    return V_tree(t, y, c) + compute_V_CW_dressed(t, y, c, mu_sq)

print("\n  Total potential minima with dressed CW:")
for yy in [0.01, 0.05, 0.1, 0.3, 0.5, 1.0]:
    res5 = minimize_scalar(lambda t: V_eff_dressed(t, yy), bounds=(0.02, t_pole - 0.005), method='bounded')
    res5b = minimize_scalar(lambda t: V_eff(t, yy), bounds=(0.02, t_pole - 0.005), method='bounded')
    print(f"  y={yy:.3f}: undressed t_min={res5b.x:.4f}, dressed t_min={res5.x:.4f}")

# Plot comparison
fig, axes = plt.subplots(1, 2, figsize=(14, 5))
for idx, yy in enumerate([0.1, 0.5]):
    ax = axes[idx]
    V_und = np.array([V_eff(t, yy) for t in t_search])
    V_dre = np.array([V_eff_dressed(t, yy) for t in t_search])
    
    ax.plot(t_search, V_und - V_und[0], 'b-', linewidth=2, label='undressed CW')
    ax.plot(t_search, V_dre - V_dre[0], 'r-', linewidth=2, label='dressed CW')
    ax.axvline(t_pole, color='k', linestyle='--', alpha=0.5, label='pole')
    ax.set_xlabel('t', fontsize=13)
    ax.set_ylabel('V_eff - V_eff(0)', fontsize=13)
    ax.set_title(f'y = {yy}', fontsize=14)
    ax.legend(fontsize=11)
    ax.grid(True, alpha=0.3)

plt.suptitle('Dressed vs undressed CW in total V_eff', fontsize=14)
plt.tight_layout()
plt.savefig('/home/codexssh/phys3/results/kahler_pole_CW_part3.png', dpi=150)
plt.close()
print("  Plot saved: results/kahler_pole_CW_part3.png")

# ============================================================
# PART 4: Analysis near the pole
# ============================================================
print("\n" + "=" * 60)
print("PART 4: Behavior near the pole")
print("=" * 60)

# Near the pole t → √3:
# K_{0bar0} = 1 - t²/3 → 0
# V_tree = y²/(1 - t²/3) → ∞
# V_CW(t) is smooth and finite at t = √3

# Proper distance σ:
# dσ/dt = sqrt(K_{0bar0}) = sqrt(1 - t²/3)
# σ(t) = ∫₀ᵗ sqrt(1 - t'²/3) dt'
# This is σ = (t/2)sqrt(1-t²/3) + (√3/2)arcsin(t/√3)
# σ_max = σ(√3) = (√3/2)·(π/2) = √3π/4

def sigma_of_t(t):
    """Proper distance from origin to t."""
    val, _ = quad(lambda tp: np.sqrt(max(1 - tp**2/3, 0)), 0, t)
    return val

sigma_max = sigma_of_t(t_pole - 1e-10)
print(f"  σ_max = {sigma_max:.6f}")
print(f"  √3·π/4 = {np.sqrt(3)*np.pi/4:.6f}")

# Near the pole: 1 - t²/3 = ε → t = √3·√(1-ε) ≈ √3(1 - ε/2)
# σ_max - σ ≈ ∫_{t}^{√3} √(1-t'²/3) dt' ≈ ∫ √ε dε' ... 
# V_tree ~ y²/ε → ∞ as ε → 0
# V_CW(√3) = finite constant

print(f"\n  V_CW at t→√3: {compute_V_CW(t_pole - 0.001, y):.8e}")
print(f"  V_tree at t→√3-0.001: {V_tree(t_pole - 0.001, y):.8e}")
print(f"  V_tree at t→√3-0.01: {V_tree(t_pole - 0.01, y):.8e}")
print(f"  V_tree at t→√3-0.1: {V_tree(t_pole - 0.1, y):.8e}")

# The tree potential diverges at the pole → V_eff → ∞ at the pole.
# The CW potential is finite. So V_eff has a wall at the pole.
# The minimum of V_eff is somewhere BEFORE the pole.

# The question: can the minimum be pushed arbitrarily close to the pole?
# As y → 0: V_tree → 0, V_CW dominates → minimum at CW minimum
# As y → ∞: V_tree dominates, V_tree increases monotonically → minimum at t=0
# For intermediate y: minimum is between 0 and √3

# Can V_eff have minimum AT the pole? No — V_tree → ∞ there.
# The pole is ALWAYS a wall, never a minimum.

# But the question asks about the limit. Let's look at V_eff'(t) = 0:
# V_tree'(t) = 2y²t/(3(1-t²/3)²) > 0 always → tree pushes toward t=0
# V_CW'(t) can be negative (pushes toward larger t)
# At the minimum: V_tree'(t_min) + V_CW'(t_min) = 0

# Near the pole: V_tree' ~ 1/(1-t²/3)² → ∞, V_CW' finite
# So V_eff' → +∞ near pole → minimum cannot be at the pole.

print("\n  Conclusion: V_tree diverges at the pole while V_CW is finite.")
print("  The pole is always a wall. The minimum is always strictly before the pole.")
print("  The tree-level Kähler pole REPELS, not attracts.")

# Check: what if we include the dressed CW which also diverges?
print(f"\n  Dressed V_CW near pole:")
for eps in [0.1, 0.01, 0.001]:
    tt = t_pole - eps
    k = K00(tt)
    y_eff = y / k
    print(f"    t = √3 - {eps}: K = {k:.6f}, y_eff = {y_eff:.4f}, V_CW_dressed = {compute_V_CW_dressed(tt, y):.6e}")

# The dressed CW also diverges! Because y_eff → ∞ as K → 0
# Let's check the behavior
t_near_pole = np.linspace(t_pole - 0.5, t_pole - 0.001, 500)
V_tree_near = np.array([V_tree(t, y) for t in t_near_pole])
V_cw_undressed_near = np.array([compute_V_CW(t, y) for t in t_near_pole])
V_cw_dressed_near = np.array([compute_V_CW_dressed(t, y) for t in t_near_pole])

fig, axes = plt.subplots(1, 2, figsize=(14, 5))

ax = axes[0]
ax.semilogy(t_near_pole, V_tree_near, 'r-', linewidth=2, label='V_tree')
ax.semilogy(t_near_pole, np.abs(V_cw_undressed_near), 'b-', linewidth=2, label='|V_CW| (undressed)')
ax.semilogy(t_near_pole, np.abs(V_cw_dressed_near), 'g-', linewidth=2, label='|V_CW| (dressed)')
ax.axvline(t_pole, color='k', linestyle='--', alpha=0.5)
ax.set_xlabel('t', fontsize=13)
ax.set_ylabel('|V|', fontsize=13)
ax.set_title(f'Potentials near pole (y={y})', fontsize=14)
ax.legend(fontsize=11)
ax.grid(True, alpha=0.3)

# Total with dressed CW
V_total_dressed_near = np.array([V_eff_dressed(t, y) for t in t_near_pole])
V_total_undressed_near = np.array([V_eff(t, y) for t in t_near_pole])

ax = axes[1]
ax.plot(t_near_pole, V_total_undressed_near, 'b-', linewidth=2, label='V_tree + V_CW (undressed)')
ax.plot(t_near_pole, V_total_dressed_near, 'r-', linewidth=2, label='V_tree + V_CW (dressed)')
ax.axvline(t_pole, color='k', linestyle='--', alpha=0.5)
ax.set_xlabel('t', fontsize=13)
ax.set_ylabel('V_eff', fontsize=13)
ax.set_title(f'Total potential near pole (y={y})', fontsize=14)
ax.legend(fontsize=11)
ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('/home/codexssh/phys3/results/kahler_pole_CW_part4.png', dpi=150)
plt.close()
print("  Plot saved: results/kahler_pole_CW_part4.png")

# ============================================================
# Check with the STANDARD CW result for O'Raifeartaigh
# ============================================================
print("\n" + "=" * 60)
print("VERIFICATION: Standard CW pseudo-modulus stabilization")
print("=" * 60)

# The standard result (Intriligator, Seiberg, Shih 0602239) is that
# the CW potential stabilizes the pseudo-modulus at the ORIGIN v=0.
# Let me verify this.

# For y << 1, the CW potential should be (ISS eq 3.10):
# V_CW ∝ y² log(1 + 4t²) + ... which has minimum at t=0.

# Let me check the sign of V_CW'(0+):
dt = 0.001
dV = compute_V_CW(dt, y) - compute_V_CW(0.01, y)
print(f"  V_CW(dt={dt}) - V_CW(0.01) = {dV:.8e}  (sign: {'positive, min at 0' if dV > 0 else 'negative, min away from 0'})")

# More careful: compute V_CW at many small t values
t_small = np.linspace(0.001, 0.5, 100)
V_small = np.array([compute_V_CW(t, y) for t in t_small])
print(f"  V_CW monotonically increasing from 0? {np.all(np.diff(V_small) > 0)}")
print(f"  V_CW monotonically decreasing from 0? {np.all(np.diff(V_small) < 0)}")

# If V_CW increases from 0, the minimum IS at the origin.
# The CW potential pushes the pseudo-modulus to v = 0.
# Then V_tree (with Kähler) = y²/(1 - t²/3) has minimum at t = 0 too.
# Both want t = 0 → minimum IS at origin, not at the pole!

# Let me recheck: maybe the standard result says the minimum IS at origin
# The "minimum at t ≈ 0.49" claim needs checking.

# Actually let me look at this more carefully for y = 1 (strong SUSY breaking)
print("\n  CW potential shape for y = 1.0:")
y_test = 1.0
t_test = np.linspace(0.01, 3.0, 300)
V_test = np.array([compute_V_CW(t, y_test) for t in t_test])
V_test -= V_test[0]

# Find if there's a minimum
idx_min = np.argmin(V_test)
print(f"  Minimum of V_CW at t = {t_test[idx_min]:.4f}")
print(f"  V_CW values at t = 0.01, 0.5, 1, 2, 3:")
for tt in [0.01, 0.5, 1.0, 2.0, 3.0]:
    print(f"    t={tt}: V_CW = {compute_V_CW(tt, y_test):.8e}")

# The issue is that for y > some value, V_CW may develop a minimum
# away from origin. Let me check the SUSY breaking parameter
# The relevant expansion parameter is y²
print("\n  Checking V_CW slope at origin for various y:")
for yy in [0.01, 0.1, 0.5, 1.0, 2.0, 5.0, 10.0]:
    slope = (compute_V_CW(0.02, yy) - compute_V_CW(0.01, yy)) / 0.01
    print(f"    y={yy:5.2f}: dV_CW/dt|_0 ≈ {slope:.6e}")

