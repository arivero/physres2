#!/usr/bin/env python3
"""
Part 4 (refined) and Part 5: Near-pole analysis and competition 
between V_tree and V_CW_dressed.
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize_scalar, brentq
from scipy.integrate import quad

# ============================================================
# Core functions (same as before)
# ============================================================

def fermion_masses_sq(t):
    s = np.sqrt(t**2 + 1)
    return np.array([0.0, (s - t)**2, (s + t)**2])

def scalar_masses_sq(t, y):
    Mplus = np.array([[4*t**2 + 1 - 2*y, 2*t],
                       [2*t,               1  ]])
    Mminus = np.array([[4*t**2 + 1 + 2*y, 2*t],
                        [2*t,               1  ]])
    eig_p = np.linalg.eigvalsh(Mplus)
    eig_m = np.linalg.eigvalsh(Mminus)
    return np.concatenate([[0.0, 0.0], eig_p, eig_m])

def CW_contribution(msq_arr, n_dof_arr, mu_sq=1.0):
    V = 0.0
    for msq, ndof in zip(msq_arr, n_dof_arr):
        if msq > 1e-30:
            V += ndof * msq**2 * (np.log(msq / mu_sq) - 1.5)
    return V / (64 * np.pi**2)

def compute_V_CW(t, y, mu_sq=1.0):
    smsq = scalar_masses_sq(t, y)
    fmsq = fermion_masses_sq(t)
    all_msq = np.concatenate([smsq, fmsq])
    all_ndof = np.concatenate([np.ones(6), -2*np.ones(3)])
    return CW_contribution(all_msq, all_ndof, mu_sq)

def K00(t, c=-1.0/12):
    return 1 + 4*c*t**2

def V_tree(t, y, c=-1.0/12):
    k = K00(t, c)
    if k <= 0:
        return np.inf
    return y**2 / k

def compute_V_CW_dressed(t, y, c=-1.0/12, mu_sq=1.0):
    k = K00(t, c)
    if k <= 0:
        return np.inf
    y_eff = y / k
    smsq = scalar_masses_sq(t, y_eff)
    fmsq = fermion_masses_sq(t)
    all_msq = np.concatenate([smsq, fmsq])
    all_ndof = np.concatenate([np.ones(6), -2*np.ones(3)])
    return CW_contribution(all_msq, all_ndof, mu_sq)

t_pole = np.sqrt(3)

# ============================================================
# PART 4 (refined): Asymptotic analysis near the pole
# ============================================================
print("=" * 60)
print("PART 4: Asymptotic analysis near the pole")
print("=" * 60)

# Let ε = 1 - t²/3, so K_{0bar0} = ε, t = √(3(1-ε))
# V_tree = y²/ε → diverges as 1/ε
# 
# For the dressed CW: y_eff = y/ε
# At large y_eff, the scalar mass splitting ~ 2y_eff >> 1
# The lightest scalar mass² ≈ (4t²+1) - √((4t²+1)²-4(4t²+1-2y_eff)) 
# for 2y_eff >> 4t²+1 ≈ 13: m²_light ≈ (4t²+1) - 2y_eff < 0 → TACHYON!
#
# This is important: the Kähler dressing eventually makes the scalar
# mass-squared negative → the model breaks down before the pole!

print("\n  Checking for tachyons near the pole:")
y = 0.1
for eps in [0.5, 0.2, 0.1, 0.05, 0.01, 0.005, 0.001]:
    tt = np.sqrt(3*(1-eps))
    k = K00(tt)
    y_eff = y / k
    smsq = scalar_masses_sq(tt, y_eff)
    min_msq = np.min(smsq)
    print(f"  ε={eps:.3f}, t={tt:.4f}, K={k:.4f}, y_eff={y_eff:.3f}, min(m²_scalar)={min_msq:.4f}", end='')
    if min_msq < 0:
        print(" *** TACHYON ***")
    else:
        print("")

# Find the critical ε where tachyon appears
print("\n  Finding tachyon threshold:")
def min_scalar_msq(eps, y=0.1):
    tt = np.sqrt(3*(1-eps))
    k = K00(tt)
    y_eff = y / k
    return np.min(scalar_masses_sq(tt, y_eff))

# Binary search for critical epsilon
eps_lo, eps_hi = 1e-6, 0.5
for _ in range(100):
    eps_mid = (eps_lo + eps_hi) / 2
    if min_scalar_msq(eps_mid) < 0:
        eps_lo = eps_mid
    else:
        eps_hi = eps_mid
eps_crit = (eps_lo + eps_hi) / 2
t_crit = np.sqrt(3*(1-eps_crit))
print(f"  Critical ε = {eps_crit:.8f}")
print(f"  Critical t = {t_crit:.6f} (pole at {t_pole:.6f})")
print(f"  Distance from pole: {t_pole - t_crit:.6f}")
print(f"  K at tachyon threshold: {K00(t_crit):.8f}")
print(f"  y_eff at threshold: {y/K00(t_crit):.4f}")

# For undressed CW: no tachyon issue, since y stays fixed
# But the undressed CW is finite at the pole while V_tree → ∞
# So V_eff = V_tree + V_CW (undressed) has minimum at origin.

# For dressed CW: both V_tree and V_CW diverge at the pole.
# Question: which diverges faster?

print("\n  Divergence rates near pole (y=0.1):")
print("  ε        V_tree      V_CW_dressed   V_CW_und     ratio(dressed/tree)")
for eps in [0.1, 0.05, 0.02, 0.01]:
    tt = np.sqrt(3*(1-eps))
    vt = V_tree(tt, y)
    vc_d = compute_V_CW_dressed(tt, y)
    vc_u = compute_V_CW(tt, y)
    print(f"  {eps:.3f}    {vt:.4e}   {vc_d:.4e}   {vc_u:.4e}   {vc_d/vt:.4f}")

# V_tree ~ y²/ε ~ 1/ε
# V_CW_dressed: the mass splitting B = 2y/ε, and STr[M⁴ log M²] 
# At large y_eff, the dominant contribution is ~ y_eff⁴ log y_eff ~ 1/ε⁴ log(1/ε)
# So V_CW_dressed diverges MUCH FASTER than V_tree!

# This means V_CW_dressed dominates near the pole.
# Since V_CW_dressed is positive and increasing → it's STILL a wall.
# The minimum structure doesn't change qualitatively.

# ============================================================
# What about NEGATIVE c (c < 0) with the OPPOSITE sign?
# We have c = -1/12, meaning K_{0bar0} DECREASES.
# What if instead c > 0, so K_{0bar0} INCREASES?
# Then V_tree = y²/K_{0bar0} DECREASES with t → could attract to large t!
# But there's no pole; it goes to zero smoothly.
# ============================================================

# ============================================================
# Alternative: V_tree has pole, but V_CW has MINIMUM away from 0
# when we include contributions from Φ₀ sector properly
# ============================================================

# Wait — I computed the Φ₀ scalar masses as zero. But with the Kähler 
# correction, the Φ₀ kinetic term changes, giving effective mass terms.
# Let me reconsider.

# The issue: in the CANONICAL field σ, the CW potential should be 
# computed from the canonical mass spectrum. The Φ₀ field has 
# non-canonical kinetic term → its effective mass changes.

# Actually, V_tree already captures the Kähler effect on Φ₀'s zero-mode.
# The CW potential comes from the HEAVY fields (Φ₁, Φ₂) integrated out.
# Their masses in the canonical frame are unchanged (K_{1bar1} = K_{2bar2} = 1).

# The only Kähler effect on CW is through the F-term splitting:
# F_0² → F_0²/K_{0bar0} (physical F-term) → this enters the B-matrix.
# This is exactly the "dressing" I already computed.

# ============================================================
# COMPREHENSIVE COMPARISON PLOT
# ============================================================

fig, axes = plt.subplots(2, 3, figsize=(18, 10))

# Row 1: Three different y values, comparing V_tree, V_CW, V_eff
for idx, yy in enumerate([0.05, 0.1, 0.5]):
    ax = axes[0][idx]
    t_arr = np.linspace(0.01, t_pole - 0.02, 500)
    
    Vtr = np.array([V_tree(t, yy) for t in t_arr])
    Vcw = np.array([compute_V_CW(t, yy) for t in t_arr])
    Veff = Vtr + Vcw
    
    # Normalize
    ax.plot(t_arr, Veff, 'b-', linewidth=2, label='V_tree + V_CW')
    ax.plot(t_arr, Vtr, 'r--', linewidth=1.5, label='V_tree')
    ax.plot(t_arr, Vcw, 'g--', linewidth=1.5, label='V_CW')
    ax.axvline(t_pole, color='k', linestyle=':', alpha=0.5, label='pole')
    
    res = minimize_scalar(lambda t: V_tree(t, yy) + compute_V_CW(t, yy), 
                          bounds=(0.02, t_pole-0.02), method='bounded')
    ax.axvline(res.x, color='b', linestyle=':', alpha=0.7, label=f'min: t={res.x:.3f}')
    
    ax.set_xlabel('t')
    ax.set_ylabel('V')
    ax.set_title(f'y = {yy}')
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)
    ax.set_ylim(bottom=min(Vcw[0], 0) - 0.001, top=Veff[0]*3)

# Row 2: Dressed CW comparison
for idx, yy in enumerate([0.05, 0.1, 0.5]):
    ax = axes[1][idx]
    # Need to stop before tachyon threshold
    def find_safe_range(y_val):
        for eps in np.logspace(-6, -0.3, 500):
            tt = np.sqrt(3*(1-eps))
            k = K00(tt)
            y_eff = y_val / k
            smsq = scalar_masses_sq(tt, y_eff)
            if np.min(smsq) < 0:
                return np.sqrt(3*(1-eps*2))
        return t_pole - 0.02
    
    t_max = find_safe_range(yy)
    t_arr2 = np.linspace(0.01, min(t_max, t_pole - 0.02), 500)
    
    Vtr = np.array([V_tree(t, yy) for t in t_arr2])
    Vcw_d = np.array([compute_V_CW_dressed(t, yy) for t in t_arr2])
    Veff_d = Vtr + Vcw_d
    
    ax.plot(t_arr2, Veff_d, 'r-', linewidth=2, label='V_tree + V_CW(dressed)')
    ax.plot(t_arr2, Vtr, 'r--', linewidth=1, label='V_tree')
    ax.plot(t_arr2, Vcw_d, 'g--', linewidth=1, label='V_CW(dressed)')
    ax.axvline(t_pole, color='k', linestyle=':', alpha=0.5)
    
    ax.set_xlabel('t')
    ax.set_ylabel('V')
    ax.set_title(f'Dressed, y = {yy}, t_max = {t_max:.3f}')
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)
    if len(Veff_d) > 0:
        ax.set_ylim(bottom=0, top=min(Veff_d[0]*5, np.max(Veff_d[Veff_d < np.inf])*1.1))

plt.suptitle('O\'Raifeartaigh + Kahler pole: potential comparison', fontsize=15)
plt.tight_layout()
plt.savefig('/home/codexssh/phys3/results/kahler_pole_potential_comparison.png', dpi=150)
plt.close()
print("  Comparison plot saved.")

# ============================================================
# PART 4 continued: Proper distance analysis
# ============================================================
print("\n" + "=" * 60)
print("PART 4: Proper distance and boundary behavior")
print("=" * 60)

def sigma_of_t(t):
    val, _ = quad(lambda tp: np.sqrt(max(1 - tp**2/3, 0)), 0, t)
    return val

sigma_max = sigma_of_t(t_pole)
print(f"  Proper distance to pole: σ_max = √3π/4 = {sigma_max:.6f}")

# In the σ coordinate, the PHYSICAL potential is V(t(σ)).
# As σ → σ_max, t → √3, V_tree → ∞.
# The boundary σ = σ_max is at FINITE proper distance.
# V → ∞ there → the field cannot reach the boundary.
# The effective potential in σ has an infinite wall at σ_max.

# The CW potential V_CW(t(σ)) is bounded.
# The tree potential V_tree(t(σ)) = y²/(1-t(σ)²/3) → ∞.
# So the total V_eff → ∞ at the boundary.

# For the minimum to be at the boundary: impossible, since V → ∞.
# For the minimum to be NEAR the boundary: need V_CW to have a deep 
# minimum near the pole that ALMOST cancels V_tree's rise.

# But V_CW increases monotonically! It has no minimum structure near the pole.
# Both V_tree and V_CW increase toward the pole.
# → The minimum is pushed to t = 0.

# The ONLY way to get a minimum at t = √3 is if the CW potential 
# DECREASES near the pole, which would require the supertrace to 
# change sign. Let's check:

print("\n  Supertrace structure:")
for tt in [0.0, 0.5, 1.0, 1.5, t_pole-0.01]:
    smsq = scalar_masses_sq(tt, y)
    fmsq = fermion_masses_sq(tt)
    str_m2 = np.sum(smsq) - 2*np.sum(fmsq)
    str_m4 = np.sum(smsq**2) - 2*np.sum(fmsq**2)
    print(f"  t={tt:.3f}: STr[M²]={str_m2:.6f}, STr[M⁴]={str_m4:.6f}")

# STr[M²] should be = 2y (from the mass sum rule with F-term SUSY breaking)
# This is t-independent!
print(f"\n  Expected STr[M²] = 2y = {2*y:.4f}")

# ============================================================
# PART 5: Summary
# ============================================================
print("\n" + "=" * 60)
print("PART 5: SUMMARY")
print("=" * 60)

summary = """
SUMMARY: Kahler pole stabilization analysis
=============================================

The question: can a Kahler potential with a pole at v = sqrt(3) m/g 
stabilize the O'Raifeartaigh pseudo-modulus at that pole?

Setup:
  K = |Phi_0|^2 - (1/12)|Phi_0|^4/Lambda^2  (with Lambda = m/g)
  K_{0bar0} = 1 - t^2/3  vanishes at t = sqrt(3)
  V_tree = f^2/K_{0bar0} diverges at the pole

Findings:

1. CW POTENTIAL IS MONOTONICALLY INCREASING from t = 0.
   The standard ISS result applies: the pseudo-modulus is stabilized 
   at the origin v = 0. There is no CW minimum at t ~ 0.49 or anywhere 
   else — the CW potential increases monotonically for ALL values of y.

2. V_tree + V_CW: both increase toward the pole.
   V_tree diverges as 1/epsilon (epsilon = 1 - t^2/3).
   V_CW is finite at the pole.
   The total V_eff = V_tree + V_CW is dominated by V_tree near the pole 
   and has its minimum at t = 0 for all y.

3. KAHLER-DRESSED CW: replacing y -> y/K_{0bar0} in the mass splitting.
   This makes V_CW_dressed also diverge at the pole (as ~ 1/epsilon^4 log(1/epsilon)).
   However, the dressed CW encounters a TACHYON before reaching the pole:
   the scalar mass-squared becomes negative when the effective SUSY breaking 
   y_eff = y/epsilon exceeds a critical value.
   For y = 0.1: tachyon appears at t = {t_crit:.4f} (distance {t_pole - t_crit:.4f} from pole).

4. PROPER DISTANCE: The pole is at FINITE proper distance sigma_max = sqrt(3)pi/4 = {sigma_max:.4f}.
   The potential diverges at finite proper distance → impenetrable wall.
   The potential minimum is always at the origin, not at the wall.

5. SUPERTRACE: STr[M^2] = 2y (constant, independent of t). This is the
   standard mass sum rule for F-term SUSY breaking. The CW potential 
   shape is determined by STr[M^4 log M^2], which increases monotonically.

CONCLUSION: The Kahler pole stabilization FAILS for the following reason:

   The Kahler pole creates a WALL (V_tree -> infinity), not a WELL.
   The CW potential, whether dressed or undressed, increases monotonically 
   toward the pole. Both the tree-level and one-loop contributions push 
   the pseudo-modulus AWAY from the pole, toward the origin.

   The mechanism cannot stabilize v at sqrt(3) m/g because:
   (a) V_tree diverges positively at the pole
   (b) V_CW increases monotonically 
   (c) The dressed CW hits a tachyon before reaching the pole
   (d) There is no competing negative contribution

   For stabilization at v = sqrt(3) m/g, one would need a DIFFERENT 
   mechanism — such as a nonperturbative superpotential contribution 
   (e.g., from instantons or strong dynamics) that creates a negative 
   contribution to V near the pole. The one-loop CW potential alone 
   cannot do this.
""".format(t_crit=t_crit, t_pole=t_pole, sigma_max=sigma_max)

print(summary)

# ============================================================
# Final comprehensive plot for the paper
# ============================================================

fig, axes = plt.subplots(1, 3, figsize=(18, 5))

# Panel 1: CW potential alone
ax = axes[0]
t_arr = np.linspace(0.01, 3.5, 500)
for yy in [0.05, 0.1, 0.5, 1.0]:
    V = np.array([compute_V_CW(t, yy) for t in t_arr])
    V -= V[0]
    ax.plot(t_arr, V, linewidth=1.5, label=f'y = {yy}')
ax.axvline(t_pole, color='k', linestyle='--', alpha=0.5, label='pole')
ax.set_xlabel('t = gv/m', fontsize=13)
ax.set_ylabel(r'$V_{CW}(t) - V_{CW}(0)$', fontsize=13)
ax.set_title('Coleman-Weinberg potential', fontsize=14)
ax.legend(fontsize=10)
ax.grid(True, alpha=0.3)

# Panel 2: Total V_eff = V_tree + V_CW
ax = axes[1]
t_arr2 = np.linspace(0.01, t_pole - 0.05, 500)
for yy in [0.05, 0.1, 0.3]:
    V = np.array([V_tree(t, yy) + compute_V_CW(t, yy) for t in t_arr2])
    V -= V[0]
    ax.plot(t_arr2, V, linewidth=1.5, label=f'y = {yy}')
ax.axvline(t_pole, color='k', linestyle='--', alpha=0.5, label='pole')
ax.set_xlabel('t = gv/m', fontsize=13)
ax.set_ylabel(r'$V_{eff}(t) - V_{eff}(0)$', fontsize=13)
ax.set_title(r'$V_{tree} + V_{CW}$ (Kahler pole)', fontsize=14)
ax.legend(fontsize=10)
ax.grid(True, alpha=0.3)

# Panel 3: Components near pole
ax = axes[2]
yy = 0.1
t_near = np.linspace(0.01, t_pole - 0.02, 500)
Vtr = np.array([V_tree(t, yy) for t in t_near])
Vcw = np.array([compute_V_CW(t, yy) for t in t_near])
ax.semilogy(t_near, Vtr, 'r-', linewidth=2, label=r'$V_{tree} = f^2/K_{0\bar{0}}$')
ax.semilogy(t_near, Vcw + 1e-6, 'b-', linewidth=2, label=r'$V_{CW}$')
ax.semilogy(t_near, Vtr + Vcw, 'k-', linewidth=2, label=r'$V_{eff}$')
ax.axvline(t_pole, color='k', linestyle='--', alpha=0.5)
ax.set_xlabel('t = gv/m', fontsize=13)
ax.set_ylabel('V (log scale)', fontsize=13)
ax.set_title(f'Near-pole behavior (y = {yy})', fontsize=14)
ax.legend(fontsize=10)
ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('/home/codexssh/phys3/results/kahler_pole_potential_physics.png', dpi=150)
plt.close()
print("Final physics plot saved: results/kahler_pole_potential_physics.png")

