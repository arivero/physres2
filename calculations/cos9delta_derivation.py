"""
Derivation: How W_3 = c_3 (det M)^3 / Lambda^18 generates cos(9 delta) potential.

Mass parametrization:
  sqrt(M_k) = z_0 [1 + sqrt(2) cos(delta + 2 pi k / 3)],  k = 1,2,3
  M_k = z_0^2 [1 + sqrt(2) cos(delta + 2 pi k/3)]^2
"""

import numpy as np
from numpy import pi, sqrt, cos, sin
from numpy.fft import fft
import sympy as sp
from sympy import (symbols, cos as scos, sin as ssin, sqrt as ssqrt,
                   expand, trigsimp, simplify, Rational, pi as spi,
                   Poly, binomial, integrate)

# Power-reduction table: cos^n(w) = sum_m c_m cos(m w)
COSN_TABLE = {
    0: {0: Rational(1)},
    1: {1: Rational(1)},
    2: {0: Rational(1, 2), 2: Rational(1, 2)},
    3: {1: Rational(3, 4), 3: Rational(1, 4)},
    4: {0: Rational(3, 8), 2: Rational(4, 8), 4: Rational(1, 8)},
    5: {1: Rational(10, 16), 3: Rational(5, 16), 5: Rational(1, 16)},
    6: {0: Rational(10, 32), 2: Rational(15, 32), 4: Rational(6, 32), 6: Rational(1, 32)},
}

delta = symbols('delta', real=True)
z0 = symbols('z_0', positive=True)

# ============================================================
# Part (a): Analytic computation of f(delta)
# ============================================================
print("=" * 70)
print("PART (a): Analytic computation of f(delta)")
print("=" * 70)

phases = [Rational(2, 3) * spi * k for k in range(1, 4)]
factors = [1 + ssqrt(2) * scos(delta + phi) for phi in phases]

print("\nThe three factors:")
for k, fk in enumerate(factors, 1):
    print(f"  Factor {k}: {fk}")

f_delta = factors[0] * factors[1] * factors[2]
f_trig = trigsimp(expand(f_delta))
print(f"\nf(delta) = {f_trig}")

# Verify against analytic formula
f_analytic = Rational(-1, 2) + scos(3 * delta) / ssqrt(2)
diff = trigsimp(f_trig - f_analytic)
print(f"Verification: f_trig - [-1/2 + cos(3d)/sqrt(2)] = {diff}")

print()
print("--- Derivation using symmetric function identities ---")
print()
print("Let c_k = cos(delta + 2 pi k/3), k = 1,2,3")
print()
print("  sum c_k = 0")
print("  sum_{i<j} c_i c_j = -3/4")
print("  prod c_k = cos(3 delta) / 4")
print()
print("prod(1 + sqrt(2) c_k)")
print("  = 1 + sqrt(2)*0 + 2*(-3/4) + 2*sqrt(2)*cos(3d)/4")
print("  = -1/2 + cos(3d)/sqrt(2)")
print()
print("RESULT:  f(delta) = -1/2 + cos(3 delta) / sqrt(2)")

# ============================================================
# Part (b): [f(delta)]^6 Fourier decomposition
# ============================================================
print()
print("=" * 70)
print("PART (b): [f(delta)]^6 and its Fourier decomposition")
print("=" * 70)

print("\nf^6 = [-1/2 + cos(3d)/sqrt(2)]^6")
print("\nLet w = 3 delta, u = cos(w).  Then f = -1/2 + u/sqrt(2).")
print()

# Binomial expansion: f^6 = sum C(6,n) (-1/2)^{6-n} (u/sqrt(2))^n
u = symbols('u')
poly_expr = sum(binomial(6, n) * Rational(-1, 2)**(6 - n) * (u / ssqrt(2))**n
                for n in range(7))
poly_expr = expand(poly_expr)
print(f"f^6 as polynomial in u = cos(w):")
print(f"  {poly_expr}")
print()

# Extract coefficients
p = Poly(poly_expr, u)
u_coeffs = p.all_coeffs()[::-1]  # u^0, u^1, ..., u^6
print("Coefficients of u^n:")
for i, c in enumerate(u_coeffs):
    print(f"  u^{i}: {c} = {float(c):.10f}")
print()

# Convert to Fourier series using power-reduction
fourier_coeffs = {}  # harmonic m (in w) -> exact coefficient
for n in range(7):
    c_n = u_coeffs[n]
    if c_n == 0:
        continue
    for m, cm in COSN_TABLE[n].items():
        if m not in fourier_coeffs:
            fourier_coeffs[m] = Rational(0)
        fourier_coeffs[m] += c_n * cm

print("Exact Fourier coefficients:  f^6 = sum_m a_{3m} cos(3m delta)")
print()
print(f"{'m(in w)':>8} {'3m(in d)':>10} {'Exact':>30} {'Decimal':>18}")
print("-" * 70)
for m in sorted(fourier_coeffs.keys()):
    c = simplify(fourier_coeffs[m])
    if c != 0:
        print(f"{m:>8d} {3*m:>10d} {str(c):>30} {float(c):>18.10f}")

print()

# Numerical verification via FFT
N_fft = 4096
d_vals = np.linspace(0, 2 * pi, N_fft, endpoint=False)
def f_num(d):
    return -0.5 + np.cos(3 * d) / np.sqrt(2)

f6_vals = f_num(d_vals)**6
fft_coeffs = fft(f6_vals) / N_fft

print("Numerical FFT verification:")
print(f"{'Harmonic':>10} {'FFT cosine coeff':>20}")
print("-" * 35)
for n in range(20):
    a = fft_coeffs[0].real if n == 0 else 2 * fft_coeffs[n].real
    if abs(a) > 1e-12:
        print(f"{n:>10d} {a:>20.10f}")

print()

# ============================================================
# Part (c): F-term derivatives
# ============================================================
print("=" * 70)
print("PART (c): F-term derivatives  dW_3/dM_k")
print("=" * 70)
print()
print("W_3 = c_3 (M_1 M_2 M_3)^3 / Lambda^18")
print()
print("dW_3/dM_k = 3 c_3 (det M)^3 / (M_k Lambda^18)")
print("          = 3 c_3 z_0^{16} [f(d)]^6 / {Lambda^18 [g_k(d)]^2}")
print()
print("where g_k(d) = 1 + sqrt(2) cos(d + 2 pi k/3),  M_k = z_0^2 g_k^2")
print()
print("The F-term scalar potential is V = sum_k |dW/dM_k|^2.")
print("From ISS: dW_ISS/dM_k = m_k  (quark masses, d-independent)")
print()
print("V = sum_k |m_k + 3 c_3 (det M)^3 / (M_k Lambda^18)|^2")
print()
print("Cross-term (linear in c_3, dominant for small c_3):")
print("  V_cross = 6 c_3 z_0^{14} / Lambda^18 * [f(d)]^6 * sum_k m_k/g_k(d)^2")
print()

# ============================================================
# Part (d): Effective potential V_eff(delta) — harmonic content
# ============================================================
print("=" * 70)
print("PART (d): Effective potential harmonics")
print("=" * 70)
print()

# Pure W_3^2 term: V_3 = 9 c_3^2 z_0^32 [f]^12 / Lambda^36 * sum 1/g_k^4
# Cross term:      V_x = 6 c_3 z_0^14 [f]^6 / Lambda^18 * sum m_k/g_k^2

# For the cross-term (degenerate m_k = m):
# V_x ~ [f(d)]^6 * sum_k 1/g_k(d)^2
# This product is well-defined (f^6 has 6th-order zeros cancelling 2nd-order poles)

# Compute h_x(d) = f(d)^6 * sum_k 1/g_k(d)^2  numerically
h_x = np.zeros(N_fft)
for k in range(1, 4):
    gk = 1 + np.sqrt(2) * np.cos(d_vals + 2 * pi * k / 3)
    h_x += f_num(d_vals)**6 / gk**2

fft_hx = fft(h_x) / N_fft
print("Cross-term shape  h_x(d) = [f]^6 * sum 1/g_k^2:")
print(f"{'Harmonic':>10} {'Cosine coeff':>18}")
print("-" * 32)
for n in range(25):
    a = fft_hx[0].real if n == 0 else 2 * fft_hx[n].real
    if abs(a) > 1e-10:
        print(f"{n:>10d} {a:>18.10f}")
print()

# Pure |W_3|^2 term: V_3 ~ [f(d)]^12 * sum 1/g_k^4
h_3 = np.zeros(N_fft)
for k in range(1, 4):
    gk = 1 + np.sqrt(2) * np.cos(d_vals + 2 * pi * k / 3)
    h_3 += f_num(d_vals)**12 / gk**4

fft_h3 = fft(h_3) / N_fft
print("Pure W_3^2 shape  h_3(d) = [f]^12 * sum 1/g_k^4:")
print(f"{'Harmonic':>10} {'Cosine coeff':>18}")
print("-" * 32)
for n in range(37):
    a = fft_h3[0].real if n == 0 else 2 * fft_h3[n].real
    if abs(a) > 1e-10:
        print(f"{n:>10d} {a:>18.10f}")
print()

# ============================================================
# Part (e): Relative amplitudes and dominance analysis
# ============================================================
print("=" * 70)
print("PART (e): Is cos(9 delta) dominant?")
print("=" * 70)
print()

# Analyze [f]^6 harmonics
a0 = float(fourier_coeffs[0])
print("[f(delta)]^6 harmonic amplitudes (oscillating terms only):")
print()
oscillating = {}
for m in sorted(fourier_coeffs.keys()):
    if m > 0:
        c = float(fourier_coeffs[m])
        if abs(c) > 1e-15:
            oscillating[3 * m] = c

max_amp = max(abs(v) for v in oscillating.values())
for h in sorted(oscillating.keys()):
    amp = oscillating[h]
    print(f"  cos({h:2d} d): amplitude = {amp:>12.8f}, "
          f"relative to max = {amp/max_amp:>8.4f}")

print()
print(f"cos(9 delta) / cos(3 delta) = {oscillating[9]/oscillating[3]:.6f}")
print(f"|cos(9 delta)| / |cos(3 delta)| = {abs(oscillating[9]/oscillating[3]):.6f}")
print()

# ============================================================
# Part (f): The Z_9 structure analysis
# ============================================================
print("=" * 70)
print("PART (f): Why cos(9 delta) is physically singled out")
print("=" * 70)
print()
print("FACT: [f(delta)]^6 contains ALL harmonics cos(3n delta), n=0..6.")
print("      cos(9 delta) is the 3rd harmonic, with ~34% of the amplitude")
print("      of the dominant cos(3 delta) term.")
print()
print("The statement 'W_3 generates cos(9 delta)' is imprecise.")
print("More precisely:")
print()
print("(1) W_3 = c_3 (det M)^3/Lambda^18 is the LOWEST Z_18-invariant")
print("    instanton operator. det M ~ Lambda^6 e^{i alpha}, where alpha")
print("    is the U(1)_A phase. (det M)^3 ~ e^{3i alpha}. In the PHASE")
print("    variable alpha of det M, this gives cos(3 alpha).")
print()
print("(2) The Koide parametrization writes M_k in terms of (z_0, delta).")
print("    The phase of det M = prod M_k involves a DIFFERENT angular")
print("    variable than delta. det M is REAL on the Koide manifold")
print("    (all M_k real and positive for physical masses).")
print()
print("(3) The cos(9 delta) arises because f(delta) = -1/2 + cos(3d)/sqrt(2)")
print("    maps the Koide angle to a tripled frequency, and then the 6th")
print("    power generates harmonics up to cos(18 delta). The cos(9 delta)")
print("    term IS present but is not dominant.")
print()
print("(4) The Z_9 structure (9 minima) requires the POTENTIAL to have")
print("    9 equally-spaced minima. With all harmonics cos(3n delta) present,")
print("    the actual number of minima depends on the relative coefficients.")
print()

# Count actual minima of [f]^6
from scipy.signal import argrelmin
f6_fine = f_num(np.linspace(0, 2*pi, 100000, endpoint=False))**6
minima = argrelmin(f6_fine, order=5, mode='wrap')[0]
print(f"Number of local minima of [f]^6 in [0, 2 pi): {len(minima)}")

maxima = argrelmin(-f6_fine, order=5, mode='wrap')[0]
print(f"Number of local maxima of [f]^6 in [0, 2 pi): {len(maxima)}")
print()

# What about V = |dW_3/dM_k|^2?
# Count minima of h_3
h3_fine = np.zeros(100000)
d_fine = np.linspace(0, 2*pi, 100000, endpoint=False)
for k in range(1, 4):
    gk = 1 + np.sqrt(2) * np.cos(d_fine + 2 * pi * k / 3)
    # Avoid division by zero
    mask = np.abs(gk) > 1e-10
    h3_fine[mask] += f_num(d_fine[mask])**12 / gk[mask]**4
    h3_fine[~mask] = 0  # f^12/g^4 = g^8 P^12 -> 0

minima_V = argrelmin(h3_fine, order=5, mode='wrap')[0]
print(f"Number of local minima of V_3 = [f]^12 * sum 1/g_k^4 in [0, 2 pi): {len(minima_V)}")
print()

# ============================================================
# SUMMARY
# ============================================================
print("=" * 70)
print("COMPLETE SUMMARY")
print("=" * 70)
print()
print("1. EXACT RESULT:")
print("   f(delta) = prod_{k=1}^3 [1 + sqrt(2) cos(delta + 2 pi k/3)]")
print("            = -1/2 + cos(3 delta)/sqrt(2)")
print()
print("2. (det M)^3 = z_0^18 [f(delta)]^6,  with Fourier expansion:")
print()
for m in sorted(fourier_coeffs.keys()):
    c = simplify(fourier_coeffs[m])
    if c != 0:
        print(f"   a_{3*m:2d} = {str(c):>20s} = {float(c):>14.10f}")
print()
print("3. F-term derivatives:")
print("   dW_3/dM_k = 3 c_3 (det M)^3 / (M_k Lambda^18)")
print()
print("4. The scalar potential V_eff(delta) contains ALL harmonics cos(3n delta)")
print("   for n = 0, 1, 2, ..., 6.  cos(9 delta) is one of them but NOT dominant.")
print()
print("5. The cos(3 delta) term has amplitude 2.9x larger than cos(9 delta).")
print("   For pure Z_9 vacuum structure, lower harmonics must be cancelled.")
print()

# ============================================================
# PLOT
# ============================================================
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

fig, axes = plt.subplots(2, 2, figsize=(14, 10))

d = np.linspace(0, 2 * pi, 1000)

# Panel 1: f(delta)
ax = axes[0, 0]
ax.plot(d / pi, f_num(d), 'b-', linewidth=2)
ax.axhline(y=0, color='gray', linestyle='--', alpha=0.5)
ax.set_xlabel(r'$\delta / \pi$')
ax.set_ylabel(r'$f(\delta)$')
ax.set_title(r'$f(\delta) = -\frac{1}{2} + \frac{\cos 3\delta}{\sqrt{2}}$')
ax.grid(True, alpha=0.3)

# Panel 2: f^6
ax = axes[0, 1]
ax.plot(d / pi, f_num(d)**6, 'r-', linewidth=2)
ax.set_xlabel(r'$\delta / \pi$')
ax.set_ylabel(r'$[f(\delta)]^6$')
ax.set_title(r'$[f(\delta)]^6 \propto (\det M)^3$')
ax.grid(True, alpha=0.3)

# Panel 3: Fourier harmonics
ax = axes[1, 0]
for m_w in [1, 2, 3, 4]:
    amp = float(fourier_coeffs.get(m_w, 0))
    if abs(amp) > 1e-15:
        ax.plot(d / pi, amp * np.cos(3 * m_w * d), linewidth=1.5,
                label=fr'$a_{{{3*m_w}}} \cos({3*m_w}\delta)$, amp={amp:.4f}')
ax.set_xlabel(r'$\delta / \pi$')
ax.set_ylabel('Amplitude')
ax.set_title(r'Leading Fourier harmonics of $[f]^6$')
ax.legend(fontsize=8)
ax.grid(True, alpha=0.3)

# Panel 4: Bar chart of amplitudes
ax = axes[1, 1]
harmonics = [3*m for m in sorted(fourier_coeffs.keys()) if m > 0]
amplitudes = [float(fourier_coeffs[m]) for m in sorted(fourier_coeffs.keys()) if m > 0]
colors = ['steelblue' if h != 9 else 'crimson' for h in harmonics]
ax.bar([str(h) for h in harmonics], [abs(a) for a in amplitudes], color=colors)
ax.set_xlabel('Harmonic $n$ in $\\cos(n\\delta)$')
ax.set_ylabel('$|a_n|$')
ax.set_title(r'Harmonic amplitudes of $[f(\delta)]^6$ (red = $\cos 9\delta$)')
ax.grid(True, alpha=0.3, axis='y')

plt.tight_layout()
plt.savefig('/home/codexssh/phys3/results/cos9delta_decomposition.png', dpi=150, bbox_inches='tight')
print("\nPlot saved to results/cos9delta_decomposition.png")
