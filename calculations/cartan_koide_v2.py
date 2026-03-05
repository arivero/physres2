#!/usr/bin/env python3
"""
Cartan-Koide identity: independent verification.

Key subtlety: δ must be extracted INDEPENDENTLY from θ to make the
verification non-tautological. We extract δ from the mass ratios
(solving the cosine parametrization) and θ from the Cartan vector.
"""

import numpy as np

# Weight vectors
w = np.array([
    [ 0.5,  1/(2*np.sqrt(3))],
    [-0.5,  1/(2*np.sqrt(3))],
    [ 0.0, -1/np.sqrt(3)    ],
])

triples = {
    "Charged leptons (e, μ, τ)": [0.51099895, 105.6583755, 1776.86],
    "Down-type quarks (d, s, b)": [4.67, 93.4, 4180.0],
    "Heavy quarks (c, b, t)":     [1270.0, 4180.0, 163000.0],
}

print("="*75)
print("INDEPENDENT VERIFICATION: δ from masses, θ from Cartan vector")
print("="*75)

for name, masses in triples.items():
    m = np.array(masses)
    sq = np.sqrt(m)
    z0 = np.sum(sq) / 3.0
    Q = np.sum(m) / np.sum(sq)**2

    # ── Extract δ from masses ONLY (no Cartan vector) ──
    # c_k = (√m_k / z₀ - 1) / √2 should equal cos(δ + 2πk/3)
    c = [(sq[k]/z0 - 1) / np.sqrt(2) for k in range(3)]

    # From c_0 = cos δ and c_1 = cos(δ+2π/3):
    # c_1 = -c_0/2 - √3 sin δ / 2
    # → sin δ = -(2c_1 + c_0)/√3
    cos_d = c[0]
    sin_d = -(2*c[1] + c[0]) / np.sqrt(3)
    delta_mass = np.arctan2(sin_d, cos_d)

    # Verify: reconstruct all three
    recon = [z0*(1 + np.sqrt(2)*np.cos(delta_mass + 2*np.pi*k/3)) for k in range(3)]
    recon_err = max(abs(recon[k] - sq[k]) for k in range(3))

    # ── Extract θ from Cartan vector ONLY ──
    r = sum(sq[k] * w[k] for k in range(3))
    R = np.linalg.norm(r)
    theta_cartan = np.arctan2(r[1], r[0])

    # ── Compare ──
    predicted_theta = np.pi/6 - delta_mass
    diff = theta_cartan - predicted_theta

    # Normalize diff to [-π, π]
    diff = (diff + np.pi) % (2*np.pi) - np.pi

    print(f"\n{'─'*65}")
    print(f"  {name}")
    print(f"{'─'*65}")
    print(f"  Q (Koide)          = {Q:.6f}")
    print(f"  z₀                 = {z0:.6f}")
    print(f"  c_k (cosines)      = [{c[0]:.6f}, {c[1]:.6f}, {c[2]:.6f}]")
    print(f"  Σ c_k              = {sum(c):.2e}  (should be 0)")
    print(f"  Σ c_k²             = {sum(x**2 for x in c):.6f}  (3/2 if Q=2/3)")
    print(f"  √m recon error     = {recon_err:.2e}")
    print(f"  δ (from masses)    = {delta_mass:.10f} rad = {np.degrees(delta_mass):.6f}°")
    print(f"  θ (from Cartan)    = {theta_cartan:.10f} rad = {np.degrees(theta_cartan):.6f}°")
    print(f"  π/6 − δ            = {predicted_theta:.10f} rad")
    print(f"  θ − (π/6 − δ)      = {diff:.2e} rad")

    # Check Σc_k² vs 3/2 (this IS the Koide test)
    sum_c2 = sum(x**2 for x in c)
    print(f"\n  NOTE: Σc_k² = {sum_c2:.6f}. Deviation from 3/2: {sum_c2 - 1.5:.4e}")
    print(f"  This measures departure from Q=2/3.")
    print(f"  But θ = π/6 − δ holds REGARDLESS of Q value: |Δθ| = {abs(diff):.2e}")

print("\n" + "="*75)
print("CONCLUSION")
print("="*75)
print("""
The relation θ = π/6 − δ is verified to machine precision for all three
triples. Crucially, the down-type quarks have Q = 0.731 (far from 2/3),
yet the angular relation still holds exactly. This is because:

  θ = π/6 − δ is a KINEMATIC identity about how weight vectors project
  cosines onto the Cartan plane. It does NOT require Q = 2/3.

What DOES require Q = 2/3:
  - The parametrization √m_k = z₀[1 + √2 cos(δ+2πk/3)] having z₀ = Σ√m/3
  - The radius formula R = z₀√6/2
  - The cos² identity Σc_k² = 3/2

The angular relation θ(δ) is purely geometric; the Koide constraint
fixes the radial part R(z₀).
""")
