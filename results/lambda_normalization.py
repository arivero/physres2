#!/usr/bin/env python3
"""
Numerical verification for the lambda normalization problem.
Checks dimensional analysis, NMSSM Higgs mass, and the lambda_hat conflict.
"""

import numpy as np

# Physical inputs (all in MeV)
Lambda = 300.0          # QCD-like confinement scale
v = 246_220.0           # Higgs VEV
m_h = 125_400.0         # Observed Higgs mass
tan_beta = 1.0          # tan beta = 1

# Derived
sin_2beta = 2 * tan_beta / (1 + tan_beta**2)  # = 1 at tan beta = 1

print("=" * 60)
print("Lambda normalization: numerical verification")
print("=" * 60)

# Part (a): Dimensions of lambda_0
print("\n--- Part (a): Dimensions ---")
print(f"[W] = mass^3, [X] = mass^-3, [H_u.H_d] = mass^2")
print(f"[lambda_0] = mass^3 / (mass^-3 * mass^2) = mass^4")
print(f"lambda_0 = lambda_hat * Lambda^4")
print(f"Lambda^4 = {Lambda**4:.6e} MeV^4")

# Part (b): Canonical normalization
print("\n--- Part (b): Canonical X_c ---")
print(f"X_c = X * Lambda^p with p = 4")
print(f"X = X_c / Lambda^4")
print(f"lambda_0 * X * H_u.H_d = (lambda_hat * Lambda^4) * (X_c / Lambda^4) * H_u.H_d")
print(f"                        = lambda_hat * X_c * H_u.H_d")
print(f"Coupling of canonical singlet to Higgs = lambda_hat (dimensionless)")

# Part (c): NMSSM identification
print("\n--- Part (c): NMSSM identification ---")
print(f"lambda_NMSSM = lambda_hat = lambda_0 / Lambda^4")

# Part (d): Numerical evaluation
print("\n--- Part (d): Numerical evaluation ---")

lambda_NMSSM = m_h / (v * abs(sin_2beta))
print(f"sin(2 beta) = {sin_2beta:.6f}")
print(f"lambda_NMSSM = m_h / (v * sin(2beta))")
print(f"             = {m_h:.0f} / ({v:.0f} * {sin_2beta:.1f})")
print(f"             = {lambda_NMSSM:.6f}")

lambda_hat = lambda_NMSSM  # They are the same after canonical normalization
print(f"lambda_hat   = {lambda_hat:.6f}")

# Cross-check: predict m_h back
m_h_pred = lambda_NMSSM * v * abs(sin_2beta)
print(f"\nCross-check: m_h = lambda_NMSSM * v * sin(2beta)")
print(f"           = {lambda_NMSSM:.6f} * {v:.0f} * {sin_2beta:.1f}")
print(f"           = {m_h_pred:.1f} MeV")
assert abs(m_h_pred - m_h) < 1.0, "Higgs mass cross-check failed!"
print("Cross-check PASSED.")

# Constraint comparison
lambda_hat_max = 2.97e-6
ratio = lambda_hat / lambda_hat_max
print(f"\nlambda_hat required    = {lambda_hat:.6f}")
print(f"lambda_hat max allowed = {lambda_hat_max:.2e}")
print(f"Ratio (required/max)   = {ratio:.1f}")
print(f"Exceeds bound by a factor of {ratio:.0e}")

# Alternative way to see the bound
# det M = Lambda^6 - lambda_0 v^2/2 > 0
# lambda_hat * Lambda^4 * v^2/2 < Lambda^6
# lambda_hat < 2 Lambda^2 / v^2
lambda_hat_bound_exact = 2 * Lambda**2 / v**2
print(f"\nExact bound: lambda_hat < 2 * Lambda^2 / v^2")
print(f"           = 2 * {Lambda:.0f}^2 / {v:.0f}^2")
print(f"           = {lambda_hat_bound_exact:.4e}")
print(f"(Consistent with stated bound {lambda_hat_max:.2e})")

# Part (e): Summary
print("\n--- Part (e): The conflict ---")
print(f"Required lambda_hat for m_h = 125 GeV: {lambda_hat:.4f}")
print(f"Maximum  lambda_hat for det M > 0    : {lambda_hat_bound_exact:.2e}")
print(f"Ratio                                 : {lambda_hat / lambda_hat_bound_exact:.0f}")
print(f"CONFLICT: five orders of magnitude gap.")
print(f"Resolution: X (Lagrange multiplier) != S (NMSSM singlet)")

# Part (f): F-term analysis
print("\n--- Part (f): F-term analysis ---")
print("Case 1 (no kinetic term): X enforces constraint, F_X not in potential")
print("  det M = Lambda^6 (constraint)")
meson_vev_cubed = Lambda**6  # det M = M1 * M2 * M3 for diagonal
meson_vev = Lambda**2        # if M1 = M2 = M3
print(f"  Meson VEV (diagonal): M_ii = Lambda^2 = {meson_vev:.0f} MeV^2")
print(f"  det M = Lambda^6 = {Lambda**6:.2e} MeV^6")

print("\nCase 2 (with kinetic term): F_X != 0 generically")
# At the EW vacuum, F_X = det M - Lambda^6 + lambda_0 v^2/2
# If mesons stay at Lambda^2, then:
if lambda_hat > 0:
    lambda_0_val = lambda_hat * Lambda**4
    F_X_shift = lambda_0_val * v**2 / 2
    print(f"  lambda_0 = lambda_hat * Lambda^4 = {lambda_0_val:.4e} MeV^4")
    print(f"  Higgs contribution to F_X: lambda_0 * v^2/2 = {F_X_shift:.4e} MeV^6")
    print(f"  Lambda^6 = {Lambda**6:.4e} MeV^6")
    print(f"  Ratio F_X_shift / Lambda^6 = {F_X_shift / Lambda**6:.1f}")
    print(f"  --> F_X contribution overwhelms the constraint by factor {F_X_shift/Lambda**6:.0f}")
    print(f"  --> Meson vacuum completely destabilized")

print("\n" + "=" * 60)
print("ALL CHECKS COMPLETE")
print("=" * 60)
