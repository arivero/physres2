#!/usr/bin/env python3
"""
Cartan-Koide identity verification.

Tasks 4-6: Numerical verification that θ = π/6 − δ for the Cartan plane
representation of Koide-parametrized mass triples, plus the cos² identity.
"""

import numpy as np

# ── Weight vectors of SU(3) fundamental ──
w = np.array([
    [ 0.5,  1/(2*np.sqrt(3))],   # w_0
    [-0.5,  1/(2*np.sqrt(3))],   # w_1
    [ 0.0, -1/np.sqrt(3)    ],   # w_2
])

print("="*70)
print("CHECK: Sum of weight vectors =", w.sum(axis=0), "(should be [0,0])")
print("="*70)

# ── Mass triples ──
triples = {
    "Charged leptons (e, μ, τ)": [0.51099895, 105.6583755, 1776.86],
    "Down-type quarks (d, s, b)": [4.67, 93.4, 4180.0],
    "Heavy quarks (c, b, t)":     [1270.0, 4180.0, 163000.0],
}

def extract_koide_params(masses):
    """Extract z₀ and δ from a mass triple using the angular parametrization.

    √m_k = z₀[1 + √2 cos(δ + 2πk/3)],  k = 0, 1, 2

    From this:
      Σ √m_k = 3 z₀   (since Σ cos(δ+2πk/3) = 0)
      z₀ = (Σ √m_k) / 3

    Then δ is extracted from the Cartan vector angle.
    """
    sq = np.sqrt(masses)
    z0 = np.sum(sq) / 3.0

    # Extract δ from the parametrization directly
    # √m_0 = z₀[1 + √2 cos δ]  →  cos δ = (√m_0/z₀ - 1)/√2
    # √m_1 = z₀[1 + √2 cos(δ + 2π/3)]
    # Use atan2 on the Cartan vector components for robustness

    # Build ratios c_k = (√m_k / z₀ - 1) / √2 = cos(δ + 2πk/3)
    c = [(sq[k]/z0 - 1)/np.sqrt(2) for k in range(3)]

    # Verify these are consistent cosines
    # From c_0 = cos δ, c_1 = cos(δ+2π/3):
    # Use the Cartan vector method: r_x, r_y computed directly
    r = sum(sq[k] * w[k] for k in range(3))

    R = np.sqrt(r[0]**2 + r[1]**2)
    theta = np.arctan2(r[1], r[0])

    # From analytics: θ = π/6 - δ  →  δ = π/6 - θ
    delta = np.pi/6 - theta

    # Also extract δ independently from the cosine formula
    cos_delta = (sq[0]/z0 - 1) / np.sqrt(2)
    sin_delta_from_1 = (sq[1]/z0 - 1) / np.sqrt(2)  # this is cos(δ+2π/3)
    # cos(δ+2π/3) = cos δ cos(2π/3) - sin δ sin(2π/3) = -cos δ/2 - √3 sin δ/2
    # → sin δ = -(2·sin_delta_from_1 + cos_delta) / √3
    sin_delta = -(2*sin_delta_from_1 + cos_delta) / np.sqrt(3)
    delta_direct = np.arctan2(sin_delta, cos_delta)

    return z0, delta, delta_direct, R, theta, r

def koide_Q(masses):
    """Compute Koide quotient Q = (Σm_k) / (Σ√m_k)²"""
    sq = np.sqrt(masses)
    return np.sum(masses) / np.sum(sq)**2

print("\n" + "="*70)
print("TASK 4: Numerical verification of θ = π/6 − δ")
print("="*70)

for name, masses in triples.items():
    masses = np.array(masses)
    z0, delta_from_theta, delta_direct, R, theta, r = extract_koide_params(masses)
    Q = koide_Q(masses)

    # Predicted R from formula
    R_predicted = z0 * np.sqrt(6) / 2

    # Check: reconstruct masses from delta_direct
    sq_recon = [z0 * (1 + np.sqrt(2)*np.cos(delta_direct + 2*np.pi*k/3)) for k in range(3)]

    print(f"\n{'─'*60}")
    print(f"  {name}")
    print(f"{'─'*60}")
    print(f"  masses      = {masses}")
    print(f"  √masses     = {np.sqrt(masses)}")
    print(f"  z₀          = {z0:.6f}")
    print(f"  Q (Koide)   = {Q:.6f}  (2/3 = {2/3:.6f})")
    print(f"  Cartan vec  = ({r[0]:.6f}, {r[1]:.6f})")
    print(f"  R (numeric) = {R:.6f}")
    print(f"  R (= z₀√6/2)= {R_predicted:.6f}")
    print(f"  R ratio     = {R/R_predicted:.10f}  (should be 1.0)")
    print(f"  θ (Cartan)  = {theta:.8f} rad = {np.degrees(theta):.4f}°")
    print(f"  δ (direct)  = {delta_direct:.8f} rad = {np.degrees(delta_direct):.4f}°")
    print(f"  δ (from θ)  = {delta_from_theta:.8f} rad")
    print(f"  π/6 − δ     = {np.pi/6 - delta_direct:.8f} rad")
    print(f"  θ − (π/6−δ) = {theta - (np.pi/6 - delta_direct):.2e} rad  ← should be 0")
    print(f"  δ_direct − δ_from_θ = {delta_direct - delta_from_theta:.2e}")

# ── Task 5: Geometric interpretation ──
print("\n" + "="*70)
print("TASK 5: Geometric interpretation")
print("="*70)
print("""
1. BLOOM = RIGID ROTATION in Cartan plane:
   Since θ = π/6 − δ, changing δ by Δδ rotates θ by −Δδ.
   The bloom (varying δ at fixed z₀) traces a CIRCLE of radius R = z₀√6/2
   in the Cartan (λ₃, λ₈) plane. The rotation is clockwise (negative sense).

2. TOPOLOGY: The constraint surface is S¹ × R₊:
   - S¹ parametrized by δ (or equivalently θ) — the angular degree of freedom
   - R₊ parametrized by z₀ — the overall mass scale
   At fixed z₀, the bloom traces a circle. Varying z₀ scales the radius.

3. THE KOIDE CONSTRAINT ⟨z²⟩ = z₀² AND R:
   The Koide condition Q = 2/3 is equivalent to ⟨z²⟩ = z₀²,
   where z_k = √m_k. This is:
     (1/3)Σ m_k = ((1/3)Σ √m_k)²

   From the parametrization: Σ m_k = z₀²[3 + 2·Σcos²(δ+2πk/3)]
   Since Σcos²(δ+2πk/3) = 3/2 (proven below), we get Σm_k = 3z₀²(1+1) = 6z₀²...

   Actually: Σm_k = z₀² Σ[1+√2 cos(δ+2πk/3)]²
           = z₀²[3 + 2√2·0 + 2·(3/2)]  = z₀²·6
   And (Σ√m_k)² = (3z₀)² = 9z₀²
   So Q = 6z₀²/(9z₀²) = 2/3.  ✓

   The Koide condition is AUTOMATIC in this parametrization — it's built in.
   The radius R = z₀√6/2 encodes the scale; Q=2/3 is the constraint that
   ALLOWS the angular parametrization to exist.
""")

# ── Task 6: The cos² identity ──
print("="*70)
print("TASK 6: Σ cos²(δ + 2πk/3) = 3/2  (identity proof + verification)")
print("="*70)

print("""
PROOF:
  cos²α = (1 + cos 2α)/2

  Σ_{k=0}^{2} cos²(δ + 2πk/3) = Σ [1 + cos(2δ + 4πk/3)] / 2
    = 3/2 + (1/2) Σ cos(2δ + 4πk/3)

  The sum Σ cos(2δ + 4πk/3) is a sum of three unit vectors at 120° separation
  (with initial phase 2δ), which vanishes by the same identity as Σ cos(δ+2πk/3)=0.

  Therefore: Σ cos²(δ + 2πk/3) = 3/2.  □

  CONNECTION TO R:
  R² = r_x² + r_y² = 6z₀²/4 · [sin²(δ+π/3) + cos²(δ+π/3)] = 6z₀²/4 = 3z₀²/2

  Meanwhile, Σ m_k = z₀²[3 + 2·(3/2)] = 6z₀² = 4R²
  So R² = (Σ m_k)/4 = (3/4)·⟨z²⟩·3 ...

  More directly: R = z₀√6/2, so R² = 3z₀²/2.
  The identity Σcos² = 3/2 ensures R is δ-INDEPENDENT.
  If the identity failed, different Koide phases would give different R at the same z₀.
  The cos² identity IS the reason bloom = pure rotation (no radial component).
""")

# Numerical verification of cos² identity
print("Numerical verification of Σ cos²(δ + 2πk/3) = 3/2:")
print("-" * 50)
test_deltas = np.linspace(0, 2*np.pi, 13)  # 13 points including endpoints
for d in test_deltas:
    s = sum(np.cos(d + 2*np.pi*k/3)**2 for k in range(3))
    print(f"  δ = {d:6.3f} rad ({np.degrees(d):7.2f}°):  Σcos² = {s:.15f}")

# ── Summary table ──
print("\n" + "="*70)
print("SUMMARY TABLE")
print("="*70)
print(f"{'Triple':<30} {'δ (rad)':>10} {'θ (rad)':>10} {'π/6−δ':>10} {'|Δ|':>12} {'R/R_pred':>10}")
print("-"*85)
for name, masses in triples.items():
    masses = np.array(masses)
    z0, _, delta, R, theta, _ = extract_koide_params(masses)
    R_pred = z0 * np.sqrt(6) / 2
    diff = abs(theta - (np.pi/6 - delta))
    short = name.split("(")[0].strip()
    print(f"{short:<30} {delta:>10.6f} {theta:>10.6f} {np.pi/6-delta:>10.6f} {diff:>12.2e} {R/R_pred:>10.8f}")

print("\nAll |Δ| values should be at machine epsilon (~10⁻¹⁶).")
print("All R/R_pred values should be 1.0 to machine precision.")
print("The identity θ = π/6 − δ is EXACT (algebraic, not numerical).")
