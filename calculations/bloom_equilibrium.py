#!/usr/bin/env python3
"""
Bloom equilibrium: competition between three-instanton and bion potentials.

SETUP:
- Koide parametrization: z_k = z_0 * (1 + sqrt(2)*cos(delta + 2*pi*k/3)), k=0,1,2
- z_0 = (sum of signed sqrt-masses)/3
- Physical triple (-s, c, b): z_s = -sqrt(m_s), z_c = sqrt(m_c), z_b = sqrt(m_b)
- Seed at delta = 3*pi/4 (one mass zero, Q = 2/3)

POTENTIALS:
  V_3inst(delta) proportional to |f(delta)|^{2p} where f(delta) = prod_k (1 + sqrt(2)*cos(delta + 2pi*k/3))
                = -1/2 + cos(3*delta)/sqrt(2)

  V_bion(delta) = [S_unsigned(delta) - 3]^2
  where S_unsigned = sum_k |1 + sqrt(2)*cos(delta + 2*pi*k/3)|
  (equals 3 when all z_k >= 0; exceeds 3 when any z_k < 0)

The three-instanton DRIVES bloom (wants large |f|, away from seed where f=0).
The bion OPPOSES bloom (wants S=3, penalizes sign flips).

We explore multiple powers p for V_3inst since the exact F-term structure
gives different effective powers depending on the field content.
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.optimize import minimize_scalar, brentq, minimize

# ============================================================
# Core functions
# ============================================================

def f_delta(delta):
    """Product identity: prod_k(1+sqrt(2)*cos(delta+2pi*k/3)) = -1/2 + cos(3*delta)/sqrt(2)"""
    return -0.5 + np.cos(3.0 * delta) / np.sqrt(2.0)

def S_unsigned(delta):
    """sum_k |1 + sqrt(2)*cos(delta + 2*pi*k/3)|"""
    total = 0.0
    for k in range(3):
        val = 1.0 + np.sqrt(2.0) * np.cos(delta + 2.0 * np.pi * k / 3.0)
        total += abs(val)
    return total

# Vectorized versions for plotting
def f_delta_v(delta):
    return -0.5 + np.cos(3.0 * delta) / np.sqrt(2.0)

def S_unsigned_v(delta):
    total = np.zeros_like(delta)
    for k in range(3):
        val = 1.0 + np.sqrt(2.0) * np.cos(delta + 2.0 * np.pi * k / 3.0)
        total += np.abs(val)
    return total

def V_3inst(delta, p=1):
    """Three-instanton: -|f(delta)|^{2p}. NEGATIVE because it drives bloom
    (the instanton LOWERS energy away from f=0)."""
    return -np.abs(f_delta(delta))**(2*p)

def V_bion(delta):
    """Bion: [S_unsigned - 3]^2. Penalizes sign flips."""
    s = S_unsigned(delta)
    return (s - 3.0)**2

def V_eff(delta, r, p=1):
    """V_eff = V_3inst + r * V_bion = -|f|^{2p} + r*(S-3)^2"""
    return V_3inst(delta, p) + r * V_bion(delta)

# ============================================================
# 1. PRODUCT IDENTITY VERIFICATION
# ============================================================
print("=" * 70)
print("1. PRODUCT IDENTITY VERIFICATION")
print("=" * 70)
for d in [0.0, 0.7, np.pi/4, np.pi/2, np.pi, 3*np.pi/4]:
    p = 1.0
    for k in range(3):
        p *= (1 + np.sqrt(2) * np.cos(d + 2*np.pi*k/3))
    f = f_delta(d)
    print(f"  delta = {d:.4f}: product = {p:.12f}, f(delta) = {f:.12f}, match = {np.isclose(p, f)}")

# ============================================================
# 2. SEED PROPERTIES
# ============================================================
print("\n" + "=" * 70)
print("2. SEED at delta = 3*pi/4 = 135 deg")
print("=" * 70)

d_seed = 3*np.pi/4
print(f"  f(seed) = {f_delta(d_seed):.2e}  (zero: V_3inst = 0)")
print(f"  S(seed) = {S_unsigned(d_seed):.10f}  (= 3: V_bion = 0)")

print("\n  Mass spectrum at seed:")
for k in range(3):
    zk = 1 + np.sqrt(2)*np.cos(d_seed + 2*np.pi*k/3)
    print(f"    k={k}: z_k/z_0 = {zk:+.8f}, m_k/z_0^2 = {zk**2:.8f}")

# The O'R masses are (0, 2-sqrt(3), 2+sqrt(3)) ~ (0, 0.268, 3.732)
# In the parametrization these are (z/z_0)^2 = (0, 0.402, 5.598)
# The ratio m_2/m_1 = 5.598/0.402 = 13.93
# Note: (2+sqrt(3))/(2-sqrt(3)) = (3.732/0.268) = 13.93 as a mass ratio
# But sqrt(m_2)/sqrt(m_1) = 2.366/0.634 = 3.73

# ============================================================
# 3. PHYSICAL ANGLE for (-s, c, b)
# ============================================================
print("\n" + "=" * 70)
print("3. PHYSICAL BLOOM ANGLE: (-s, c, b) triple")
print("=" * 70)

m_s, m_c, m_b = 93.4, 1275.0, 4180.0  # MeV (PDG MSbar)

zs, zc, zb = -np.sqrt(m_s), np.sqrt(m_c), np.sqrt(m_b)
z0 = (zs + zc + zb) / 3.0

print(f"  m_s = {m_s}, m_c = {m_c}, m_b = {m_b} MeV")
print(f"  z_s = {zs:.6f}, z_c = {zc:.6f}, z_b = {zb:.6f}")
print(f"  z_0 = {z0:.6f}")

# Koide Q = sum(m_k) / (sum z_k)^2
Q = (m_s + m_c + m_b) / (zs + zc + zb)**2
print(f"  Koide Q = {Q:.8f}  (2/3 = 0.66666667, dev = {abs(Q-2/3):.6f} = {abs(Q-2/3)*100/Q:.2f}%)")

# Extract delta from k=0 entry: z_s/z_0 = 1 + sqrt(2)*cos(delta)
cos_d = (zs/z0 - 1) / np.sqrt(2)
print(f"  cos(delta) = {cos_d:.10f}")

# delta in [0, pi] since cos_d < 0
d_phys = np.arccos(np.clip(cos_d, -1, 1))

# Verify: also try 2pi - d_phys
for d_try in [d_phys, 2*np.pi - d_phys]:
    zk_recon = [z0*(1+np.sqrt(2)*np.cos(d_try + 2*np.pi*k/3)) for k in range(3)]
    match = np.isclose(zk_recon[0], zs, rtol=1e-6)
    if match:
        d_phys = d_try
        break

print(f"  delta_phys = {d_phys:.10f} rad = {np.degrees(d_phys):.6f} deg = {d_phys/np.pi:.8f} * pi")

print("\n  Reconstruction (exact for k=0, approximate for k=1,2 since Q != 2/3):")
for k in range(3):
    zk = z0*(1+np.sqrt(2)*np.cos(d_phys + 2*np.pi*k/3))
    mk = zk**2
    names = ['s (signed)', 'c', 'b']
    targets = [m_s, m_c, m_b]
    print(f"    k={k}: z_k = {zk:+.4f}, m_k = {mk:.2f} MeV  (target {names[k]} = {targets[k]:.2f}, err = {abs(mk-targets[k])/targets[k]*100:.2f}%)")

# Potential values at physical angle
print(f"\n  f(d_phys) = {f_delta(d_phys):.10f}")
print(f"  S(d_phys) = {S_unsigned(d_phys):.10f}")
print(f"  V_bion(d_phys) = {V_bion(d_phys):.10f}")

# ============================================================
# 4. BOUNDARY ANALYSIS
# ============================================================
print("\n" + "=" * 70)
print("4. REGIONS WHERE V_BION = 0 (all z_k >= 0)")
print("=" * 70)

# z_k/z_0 = 0 when cos(delta + 2pi*k/3) = -1/sqrt(2)
# i.e., delta + 2pi*k/3 = 3pi/4 or 5pi/4 (mod 2pi)
boundaries = sorted(set(
    (base - 2*np.pi*k/3) % (2*np.pi)
    for k in range(3)
    for base in [3*np.pi/4, 5*np.pi/4]
))
print("  Boundary angles (one z_k = 0):")
for d in boundaries:
    print(f"    {np.degrees(d):7.2f} deg")

print("\n  Positive-region intervals (V_bion = 0):")
for i in range(len(boundaries)):
    d1 = boundaries[i]
    d2 = boundaries[(i+1) % len(boundaries)]
    if d2 <= d1:
        d2 += 2*np.pi
    d_mid = (d1 + d2) / 2
    d_mid_mod = d_mid % (2*np.pi)
    all_pos = all(1+np.sqrt(2)*np.cos(d_mid_mod+2*np.pi*k/3) >= 0 for k in range(3))
    if all_pos:
        print(f"    [{np.degrees(d1):.1f}, {np.degrees(d2 % (2*np.pi)):.1f}] deg  (width = {np.degrees(d2-d1):.1f} deg)")

print(f"\n  delta_phys = {np.degrees(d_phys):.2f} deg is {'INSIDE' if V_bion(d_phys) < 1e-10 else 'OUTSIDE'} the positive region")
print(f"  (V_bion = {V_bion(d_phys):.6f})")

# ============================================================
# 5. ZEROS AND EXTREMA OF V_3inst (= -|f|^{2p})
# ============================================================
print("\n" + "=" * 70)
print("5. STRUCTURE OF V_3inst = -|f(delta)|^{2p}")
print("=" * 70)

# Zeros of f: cos(3d) = sqrt(2)/2 → 3d = pi/4, 7pi/4 + 2n*pi
# i.e., d = pi/12, 7pi/12, 3pi/4 (=pi/12 + 2pi/3), ...
print("  Zeros of f (V_3inst = 0, local MAXIMA of V_3inst since V is negative):")
f_zeros = sorted(set(
    d % (2*np.pi)
    for n in range(6)
    for d in [(np.pi/4 + 2*n*np.pi)/3, (7*np.pi/4 + 2*n*np.pi)/3]
    if 0 <= d % (2*np.pi) < 2*np.pi
))
for d in f_zeros:
    in_pos = all(1+np.sqrt(2)*np.cos(d+2*np.pi*k/3) >= -1e-10 for k in range(3))
    print(f"    delta = {np.degrees(d):7.2f} deg, f = {f_delta(d):+.2e}, in_positive_region = {in_pos}")

# Extrema of f: df/ddelta = 0 → -3*sin(3d)/sqrt(2) = 0 → sin(3d) = 0 → 3d = n*pi
print("\n  Extrema of |f| (minima of V_3inst, i.e., deepest wells):")
f_extrema = [(n*np.pi/3, f_delta(n*np.pi/3)) for n in range(7) if 0 <= n*np.pi/3 <= 2*np.pi]
for d, fv in f_extrema:
    print(f"    delta = {np.degrees(d):7.1f} deg, f = {fv:+.8f}, |f| = {abs(fv):.8f}")

# Maximum |f| values:
f_max = -0.5 + 1/np.sqrt(2)   # at cos(3d) = 1, i.e., d = 0, 2pi/3, 4pi/3
f_min = -0.5 - 1/np.sqrt(2)   # at cos(3d) = -1, i.e., d = pi/3, pi, 5pi/3
print(f"\n  |f|_max_type_1 = {abs(f_max):.8f} at d = 0, 120, 240 deg")
print(f"  |f|_max_type_2 = {abs(f_min):.8f} at d = 60, 180, 300 deg")
print(f"  V_3inst deepest wells at d = 60, 180, 300 deg (|f| = {abs(f_min):.4f})")

# ============================================================
# 6. EQUILIBRIUM: V_eff = -|f|^{2p} + r*(S-3)^2
# ============================================================
print("\n" + "=" * 70)
print("6. EQUILIBRIUM ANALYSIS")
print("=" * 70)

# The three-instanton drives toward large |f| (d = 60, 180, 300 deg are deepest).
# The bion penalizes leaving the positive region.
# The physical angle (159 deg) is just OUTSIDE the positive region boundary at 135 deg.
# The nearest V_3inst minimum is at d = 180 deg (|f| = 1.207).

print(f"\n  Physical angle:    {np.degrees(d_phys):.2f} deg")
print(f"  Nearest f-zero:   135.00 deg (= seed)")
print(f"  Nearest |f|-max:  180.00 deg")
print(f"  Distance to seed: {np.degrees(d_phys - d_seed):.2f} deg")
print(f"  Distance to 180:  {np.degrees(np.pi - d_phys):.2f} deg")

# The equilibrium should be between the seed (135 deg) and the f-maximum (180 deg),
# pulled toward 180 by the instanton and pulled back by the bion.
# delta_phys = 159 deg is at fractional position (159-135)/(180-135) = 24/45 = 0.53
# So it's slightly past halfway toward the instanton well.

frac = (np.degrees(d_phys) - 135) / (180 - 135)
print(f"  Fractional position in [135, 180] interval: {frac:.4f}")

# Detailed scan in the [135, 180] range
print(f"\n  Detailed V_eff in [120, 200] deg range (p=1, selected r values):")
print(f"  {'delta (deg)':>12s}  {'V_3inst':>12s}  {'V_bion':>12s}  {'r=0.1':>12s}  {'r=1':>12s}  {'r=10':>12s}")

for d_deg in range(120, 201, 5):
    d = np.radians(d_deg)
    v3 = V_3inst(d, p=1)
    vb = V_bion(d)
    print(f"  {d_deg:12d}  {v3:12.6f}  {vb:12.6f}  {v3+0.1*vb:12.6f}  {v3+1.0*vb:12.6f}  {v3+10*vb:12.6f}")

# ============================================================
# 7. Find minima for various r and p
# ============================================================
print("\n" + "=" * 70)
print("7. MINIMA OF V_EFF(delta) = -|f|^{2p} + r*(S-3)^2")
print("=" * 70)

def find_global_min(r, p=1, lo=0.01, hi=2*np.pi-0.01, n_scan=5000):
    """Find global minimum of V_eff in [lo, hi]."""
    deltas = np.linspace(lo, hi, n_scan)
    veffs = np.array([V_eff(d, r, p) for d in deltas])

    # Find all local minima
    mins = []
    for i in range(1, len(veffs)-1):
        if veffs[i] < veffs[i-1] and veffs[i] < veffs[i+1]:
            try:
                result = minimize_scalar(lambda d: V_eff(d, r, p),
                                        bounds=(deltas[max(0,i-3)], deltas[min(n_scan-1,i+3)]),
                                        method='bounded')
                if result.success:
                    mins.append((result.x, result.fun))
            except:
                mins.append((deltas[i], veffs[i]))

    if not mins:
        idx = np.argmin(veffs)
        return deltas[idx], veffs[idx]

    mins.sort(key=lambda x: x[1])
    return mins[0]

# Focus on the [100, 200] range near the physical angle
print("\n  Focused search in [100, 200] deg:")
print(f"  {'p':>4s}  {'r':>10s}  {'d_min (deg)':>12s}  {'V_eff':>12s}  {'|d-d_phys|':>12s}")
print(f"  {'----':>4s}  {'----------':>10s}  {'----------':>12s}  {'----------':>12s}  {'----------':>12s}")

for p in [1, 2, 3, 6]:
    for r in [0.01, 0.1, 0.5, 1.0, 2.0, 5.0, 10.0, 50.0, 100.0]:
        d_min, v_min = find_global_min(r, p, lo=np.radians(100), hi=np.radians(200))
        gap_deg = abs(np.degrees(d_min) - np.degrees(d_phys))
        marker = "  <---" if gap_deg < 2 else ""
        print(f"  {p:4d}  {r:10.2f}  {np.degrees(d_min):12.4f}  {v_min:12.6f}  {gap_deg:12.4f}{marker}")
    print()

# ============================================================
# 8. For each p, find r* where d_min = d_phys
# ============================================================
print("\n" + "=" * 70)
print("8. FINDING r* WHERE d_min = d_phys FOR EACH p")
print("=" * 70)

for p in [1, 2, 3, 6]:
    print(f"\n  p = {p}:")

    # The gradient at d_phys: dV/ddelta = 0 gives the equilibrium condition
    # dV_3inst/dd + r * dV_bion/dd = 0 at d = d_phys
    # → r = -dV_3inst/dd / dV_bion/dd at d = d_phys

    eps = 1e-8
    dv3 = (V_3inst(d_phys+eps, p) - V_3inst(d_phys-eps, p)) / (2*eps)
    dvb = (V_bion(d_phys+eps) - V_bion(d_phys-eps)) / (2*eps)

    if abs(dvb) > 1e-15:
        r_analytic = -dv3 / dvb
        print(f"    Analytic estimate (gradient=0): r* = {r_analytic:.6f}")

        # Verify: is this a minimum (not maximum)?
        d2v3 = (V_3inst(d_phys+eps, p) - 2*V_3inst(d_phys, p) + V_3inst(d_phys-eps, p)) / eps**2
        d2vb = (V_bion(d_phys+eps) - 2*V_bion(d_phys) + V_bion(d_phys-eps)) / eps**2
        d2v = d2v3 + r_analytic * d2vb
        print(f"    d2V/dd2 at d_phys with r*: {d2v:.6f} ({'MINIMUM' if d2v > 0 else 'MAXIMUM'})")

        if r_analytic > 0 and d2v > 0:
            # Verify numerically
            d_check, v_check = find_global_min(r_analytic, p, lo=np.radians(100), hi=np.radians(200))
            print(f"    Numerical check: d_min = {np.degrees(d_check):.4f} deg (target {np.degrees(d_phys):.4f} deg)")
            print(f"    Residual: {abs(np.degrees(d_check) - np.degrees(d_phys)):.4f} deg")
            print(f"    r* = {r_analytic:.6f}")
            print(f"    Is r* of order 1? {'YES' if 0.1 < r_analytic < 10 else 'NO'}")

            # Also find global minimum (not just local in [100,200])
            d_global, v_global = find_global_min(r_analytic, p)
            print(f"    GLOBAL min at r*: d = {np.degrees(d_global):.4f} deg, V = {v_global:.6f}")
            print(f"    Is [100,200] min the global min? {'YES' if np.isclose(v_check, v_global, rtol=0.01) else 'NO (global minimum elsewhere!)'}")
    else:
        print(f"    dV_bion/dd = 0 at d_phys (critical point of bion)")

# ============================================================
# 9. Angular distances for r = 0.1, 1, 10
# ============================================================
print("\n" + "=" * 70)
print("9. ANGULAR DISTANCES (focused on [100,200] deg)")
print("=" * 70)

for p in [1, 2, 3, 6]:
    print(f"\n  p = {p}:")
    for r in [0.1, 1.0, 10.0]:
        d_min, v_min = find_global_min(r, p, lo=np.radians(100), hi=np.radians(200))
        gap = abs(np.degrees(d_min) - np.degrees(d_phys))
        print(f"    r = {r:6.1f}: d_min = {np.degrees(d_min):.4f} deg, |gap| = {gap:.4f} deg")

# ============================================================
# 10. v_0-doubling angle
# ============================================================
print("\n" + "=" * 70)
print("10. v_0-DOUBLING ANALYSIS")
print("=" * 70)

# V_bion = 0 when S_unsigned = 3 (all z_k positive, the "positive region")
# The boundary of the positive region is where one z_k = 0.
#
# The v_0-doubling relation (sqrt(m_b) = 3*sqrt(m_s) + sqrt(m_c)) is about
# z_0 doubling, not delta. At the angular level, the bion constrains delta
# to be NEAR the positive region boundary.
#
# Physical delta = 159 deg is 24 deg past the boundary at 135 deg.
# S_unsigned at delta_phys:

print(f"  S_unsigned(d_phys) = {S_unsigned(d_phys):.8f}")
print(f"  S_unsigned(d_seed) = {S_unsigned(d_seed):.8f} (= 3)")
print(f"  Excess: S - 3 = {S_unsigned(d_phys) - 3:.8f}")
print(f"  V_bion(d_phys) = {V_bion(d_phys):.8f}")

# The excess S - 3 is the "cost" of the sign flip. It equals 2*|z_s/z_0|:
z_s_over_z0 = 1 + np.sqrt(2)*np.cos(d_phys)
print(f"\n  z_s/z_0 = {z_s_over_z0:.8f} (negative: sign flip)")
print(f"  2*|z_s/z_0| = {2*abs(z_s_over_z0):.8f}")
print(f"  S - 3 = {S_unsigned(d_phys) - 3:.8f}")
print(f"  Match: {np.isclose(S_unsigned(d_phys) - 3, 2*abs(z_s_over_z0))}")

# The v_0-doubling angle is where the SEED evolves from 3pi/4 to d_phys,
# with z_0 doubling. The constraint is that the bion energy is "affordable"
# given the instanton driving force.

# ============================================================
# 11. PLOTS
# ============================================================
print("\n" + "=" * 70)
print("11. GENERATING PLOTS")
print("=" * 70)

delta_arr = np.linspace(0, 2*np.pi, 3000)
f_arr = f_delta_v(delta_arr)
S_arr = S_unsigned_v(delta_arr)

fig, axes = plt.subplots(2, 2, figsize=(16, 12))

# Plot 1: f(delta) and S(delta)
ax = axes[0, 0]
ax.plot(np.degrees(delta_arr), f_arr, 'b-', linewidth=2, label=r'$f(\delta)$')
ax.plot(np.degrees(delta_arr), S_arr - 3, 'r-', linewidth=2, label=r'$S(\delta) - 3$')
ax.axhline(0, color='gray', linewidth=0.5)
ax.axvline(np.degrees(d_phys), color='green', linestyle='--', linewidth=2, label=f'phys ({np.degrees(d_phys):.0f} deg)')
ax.axvline(135, color='gray', linestyle=':', linewidth=2, label='seed (135 deg)')
ax.set_xlabel('delta (deg)', fontsize=13)
ax.set_ylabel('Value', fontsize=13)
ax.set_title('Angular functions f(delta) and S(delta)-3', fontsize=13)
ax.legend(fontsize=10)
ax.set_xlim(0, 360)
ax.grid(True, alpha=0.3)

# Plot 2: Individual potentials in [100, 220] deg range
ax = axes[0, 1]
d_zoom = np.linspace(np.radians(100), np.radians(220), 500)
v3_zoom = np.array([-abs(f_delta(d))**2 for d in d_zoom])
vb_zoom = np.array([V_bion(d) for d in d_zoom])

ax.plot(np.degrees(d_zoom), v3_zoom, 'b-', linewidth=2, label=r'$-|f|^2$ (instanton)')
ax.plot(np.degrees(d_zoom), vb_zoom, 'r-', linewidth=2, label=r'$(S-3)^2$ (bion)')
ax.axvline(np.degrees(d_phys), color='green', linestyle='--', linewidth=2)
ax.axvline(135, color='gray', linestyle=':', linewidth=2)
ax.axvline(180, color='purple', linestyle=':', linewidth=1, alpha=0.5)
ax.set_xlabel('delta (deg)', fontsize=13)
ax.set_ylabel('Potential', fontsize=13)
ax.set_title('Potentials near physical angle (p=1)', fontsize=13)
ax.legend(fontsize=10)
ax.grid(True, alpha=0.3)

# Plot 3: V_eff for various r (p=1), zoomed
ax = axes[1, 0]
for r, color in [(0.1, 'blue'), (0.5, 'orange'), (1.0, 'red'), (5.0, 'green'), (20.0, 'purple')]:
    veff = np.array([V_eff(d, r, p=1) for d in d_zoom])
    ax.plot(np.degrees(d_zoom), veff, color=color, linewidth=1.5, label=f'r = {r}')
ax.axvline(np.degrees(d_phys), color='green', linestyle='--', linewidth=2, alpha=0.5)
ax.axvline(135, color='gray', linestyle=':', linewidth=2, alpha=0.5)
ax.set_xlabel('delta (deg)', fontsize=13)
ax.set_ylabel(r'$V_{eff}$', fontsize=13)
ax.set_title(r'$V_{eff} = -|f|^2 + r(S-3)^2$, p=1', fontsize=13)
ax.legend(fontsize=10)
ax.grid(True, alpha=0.3)

# Plot 4: d_min vs r for multiple p
ax = axes[1, 1]
r_scan = np.logspace(-2, 3, 300)
for p, color, ls in [(1, 'blue', '-'), (2, 'red', '-'), (3, 'green', '--'), (6, 'purple', ':')]:
    d_mins_p = []
    for r in r_scan:
        dm, _ = find_global_min(r, p, lo=np.radians(100), hi=np.radians(200), n_scan=1000)
        d_mins_p.append(np.degrees(dm))
    ax.semilogx(r_scan, d_mins_p, color=color, linewidth=2, linestyle=ls, label=f'p = {p}')

ax.axhline(np.degrees(d_phys), color='green', linestyle='--', linewidth=2, alpha=0.7, label=f'phys ({np.degrees(d_phys):.1f} deg)')
ax.axhline(135, color='gray', linestyle=':', linewidth=2, alpha=0.5, label='seed (135 deg)')
ax.axhline(180, color='purple', linestyle=':', linewidth=1, alpha=0.3)
ax.set_xlabel('r = B/A', fontsize=13)
ax.set_ylabel('d_min (deg)', fontsize=13)
ax.set_title('Equilibrium angle vs r', fontsize=13)
ax.legend(fontsize=9, loc='best')
ax.set_ylim(130, 185)
ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('/home/codexssh/phys3/calculations/bloom_equilibrium.png', dpi=150, bbox_inches='tight')
print("  Saved: bloom_equilibrium.png")

# ============================================================
# 12. FINAL SUMMARY
# ============================================================
print("\n" + "=" * 70)
print("FINAL SUMMARY")
print("=" * 70)

print(f"""
  FORMULATION:
    V_eff(delta) = -|f(delta)|^{{2p}} + r * [S(delta) - 3]^2
    f(delta) = -1/2 + cos(3*delta)/sqrt(2)
    S(delta) = sum_k |1 + sqrt(2)*cos(delta + 2*pi*k/3)|

  KEY ANGLES:
    Seed:      delta = 3*pi/4 = 135.00 deg
    Physical:  delta = {np.degrees(d_phys):.4f} deg
    Nearest instanton well: 180.00 deg

  PHYSICAL ANGLE is {np.degrees(d_phys - d_seed):.2f} deg past seed,
  {np.degrees(np.pi - d_phys):.2f} deg before the instanton well at 180 deg.
  Fractional position in [135, 180]: {(np.degrees(d_phys)-135)/45:.4f}

  POTENTIAL VALUES at delta_phys:
    f(d_phys) = {f_delta(d_phys):.6f}
    S(d_phys) = {S_unsigned(d_phys):.6f}
    V_bion    = {V_bion(d_phys):.6f}

  ANALYTIC r* (gradient condition dV/ddelta = 0 at d_phys):""")

for p in [1, 2, 3, 6]:
    eps = 1e-8
    dv3 = (V_3inst(d_phys+eps, p) - V_3inst(d_phys-eps, p)) / (2*eps)
    dvb = (V_bion(d_phys+eps) - V_bion(d_phys-eps)) / (2*eps)
    d2v3 = (V_3inst(d_phys+eps, p) - 2*V_3inst(d_phys, p) + V_3inst(d_phys-eps, p)) / eps**2
    d2vb = (V_bion(d_phys+eps) - 2*V_bion(d_phys) + V_bion(d_phys-eps)) / eps**2
    if abs(dvb) > 1e-15:
        r_star = -dv3/dvb
        d2v = d2v3 + r_star * d2vb
        is_min = d2v > 0
        print(f"    p = {p}: r* = {r_star:.6f}, d2V/dd2 = {d2v:.4f} ({'MIN' if is_min else 'MAX'}), order-1? {'YES' if 0.1 < r_star < 10 else 'NO'}")

print(f"""
  INTERPRETATION:
    The three-instanton potential drives delta from the seed (135 deg)
    toward the nearest well at 180 deg where |f| is maximal.
    The bion Kahler potential resists, penalizing the sign flip of z_s.
    The equilibrium r* gives the relative strength needed to park
    the angle at {np.degrees(d_phys):.1f} deg (the physical value).
""")
