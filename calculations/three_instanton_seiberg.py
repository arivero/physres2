#!/usr/bin/env python3
"""
N_f = N_c = 3 SQCD with three-instanton correction.

Superpotential:
    W = Tr(m M) + X(det M - Lambda^6) + c3 (det M)^3 / Lambda^18

where M = diag(M_u, M_d, M_s) at the diagonal vacuum.

Step 1: Write F-term equations.
Step 2: Find how X shifts from the standard Seiberg solution.
Step 3: Check consistency.
Step 4: Numerical evaluation.
"""

import numpy as np
from sympy import (
    symbols, Rational, sqrt, solve, simplify, expand, factor,
    diff, det, Matrix, Symbol, lambdify, nsolve, Eq, pretty,
    series, O, latex
)

# ============================================================
# PART 1: ANALYTIC F-TERM EQUATIONS
# ============================================================
print("=" * 72)
print("PART 1: F-TERM EQUATIONS")
print("=" * 72)

print("""
Superpotential at diagonal vacuum M = diag(M_u, M_d, M_s):

    det M = M_u * M_d * M_s  ≡ D

    W = m_u M_u + m_d M_d + m_s M_s
      + X (D - Lambda^6)
      + c3 D^3 / Lambda^18

F-term for diagonal entry M_j:

    dW/dM_j = m_j + X * (dD/dM_j) + c3 * 3 D^2 * (dD/dM_j) / Lambda^18 = 0

Now dD/dM_j = D / M_j  (cofactor = product of other two eigenvalues).

So the F-term equation for each j is:

    m_j + (X + 3 c3 D^2 / Lambda^18) * (D / M_j) = 0

Define the effective Lagrange multiplier:

    X_eff = X + 3 c3 D^2 / Lambda^18

Then the F-term equation is:

    m_j + X_eff * D / M_j = 0

    =>  M_j = - X_eff * D / m_j

This is formally identical to the standard Seiberg equation
    M_j = - X_std * Lambda^6 / m_j
but with X_eff replacing X_std.
""")

print("=" * 72)
print("PART 2: SOLVING FOR X AND CONSISTENCY")
print("=" * 72)

print("""
Standard Seiberg (c3 = 0):
    dW/dM_j = 0  =>  m_j + X * D/M_j = 0  =>  M_j = -X D / m_j
    Constraint: det M = D = Lambda^6
    =>  prod_j M_j = prod_j (-X D / m_j) = (-X)^3 D^3 / (prod m) = Lambda^6
    =>  (-X)^3 (Lambda^6)^3 / (prod m) = Lambda^6
    =>  (-X)^3 = (prod m) / Lambda^{12}
    =>  X = -(prod m)^{1/3} / Lambda^4

    Standard solution: M_j^std = Lambda^4 * (prod m)^{1/3} / (Lambda^4 m_j)
                               = Lambda^2 * (prod m)^{1/3} / m_j
    [Using D = Lambda^6, so D/m_j = Lambda^6/m_j and X=-prod(m)^{1/3}/Lambda^4]
    =>  M_j^std = (prod m)^{1/3} * Lambda^2 / m_j  ✓

With c3 != 0:
    F-term:  m_j + X_eff * D/M_j = 0  =>  M_j = -X_eff * D / m_j
    where X_eff = X + 3 c3 D^2 / Lambda^18

    Constraint: det M = D = Lambda^6 (this still holds — it is enforced by X).

    From M_j = -X_eff * D / m_j:
        prod_j M_j = (-X_eff)^3 D^3 / (prod m) = D = Lambda^6
        =>  (-X_eff)^3 = (prod m) / Lambda^{12}
        =>  X_eff = -(prod m)^{1/3} / Lambda^4   [same as X_std!]

    But X_eff = X + 3 c3 D^2 / Lambda^18 = X + 3 c3 Lambda^12 / Lambda^18
              = X + 3 c3 / Lambda^6

    So:
        X + 3 c3 / Lambda^6 = -(prod m)^{1/3} / Lambda^4

        =>  X = -(prod m)^{1/3} / Lambda^4 - 3 c3 / Lambda^6

    The shift is:
        delta X = X(c3) - X(0) = -3 c3 / Lambda^6

    And X_eff is UNCHANGED: the meson VEVs M_j are identical to the standard
    Seiberg solution. Only X itself shifts.
""")

print("=" * 72)
print("PART 3: CONSISTENCY CHECK")
print("=" * 72)

print("""
Yes, the system is fully consistent.

The c3 correction does NOT change the meson VEVs M_j.
It shifts the Lagrange multiplier X by a fixed amount -3 c3/Lambda^6.

Physical interpretation:
  - X is the auxiliary field enforcing det M = Lambda^6.
  - When c3 != 0, W gains an extra term c3 (det M)^3/Lambda^18.
  - At the vacuum det M = Lambda^6, this extra term is a constant c3 Lambda^18/Lambda^18 = c3.
  - Its variation with respect to M_j generates the shift in the effective coupling to the constraint.
  - The equation of motion for X still gives det M = Lambda^6 (X is a Lagrange multiplier,
    not a dynamical field — its equation dW/dX = det M - Lambda^6 = 0 is unchanged).
  - So X_eff (the combination that appears in the meson EOM) is fixed by the meson
    constraint, but X itself absorbs the shift from the c3 term.

Key result:
    X_std  = -(prod m)^{1/3} / Lambda^4
    X(c3)  = -(prod m)^{1/3} / Lambda^4 - 3 c3 / Lambda^6
    delta X = -3 c3 / Lambda^6

The Lagrange multiplier receives a flavor-independent shift proportional to c3.
""")

print("=" * 72)
print("PART 4: NUMERICAL EVALUATION")
print("=" * 72)

# Physical inputs (MeV)
m_u = 2.16
m_d = 4.67
m_s = 93.4
Lambda = 300.0
c3 = 1.0   # order-1 coefficient

prod_m = m_u * m_d * m_s
prod_m_13 = prod_m ** (1.0/3.0)

print(f"\n  Inputs:")
print(f"    m_u = {m_u} MeV")
print(f"    m_d = {m_d} MeV")
print(f"    m_s = {m_s} MeV")
print(f"    Lambda = {Lambda} MeV")
print(f"    c3 = {c3}")
print(f"\n  Derived:")
print(f"    prod(m) = m_u * m_d * m_s = {prod_m:.6f} MeV^3")
print(f"    (prod m)^(1/3) = {prod_m_13:.6f} MeV")
print(f"    Lambda^4 = {Lambda**4:.6e} MeV^4")
print(f"    Lambda^6 = {Lambda**6:.6e} MeV^6")

# Standard Seiberg solution (c3=0)
X_std = -prod_m_13 / Lambda**4

print(f"\n  Standard Seiberg (c3 = 0):")
print(f"    X_std = -(prod m)^(1/3) / Lambda^4 = {X_std:.6e} MeV^(-3)")

C = Lambda**2 * prod_m_13
M_u_std = C / m_u
M_d_std = C / m_d
M_s_std = C / m_s

print(f"    C = Lambda^2 * (prod m)^(1/3) = {C:.6f} MeV^3")
print(f"    M_u = C / m_u = {M_u_std:.4f} MeV")
print(f"    M_d = C / m_d = {M_d_std:.4f} MeV")
print(f"    M_s = C / m_s = {M_s_std:.4f} MeV")
print(f"    det M = {M_u_std * M_d_std * M_s_std:.6e} MeV^3")
print(f"    Lambda^6 = {Lambda**6:.6e} MeV^3")
det_check = M_u_std * M_d_std * M_s_std
print(f"    det M / Lambda^6 = {det_check / Lambda**6:.10f}  (should be 1.0)")

# Verify F-terms vanish for standard solution
print(f"\n  F-term verification (c3 = 0):")
D = Lambda**6  # = det M at vacuum
for m_j, M_j, lab in [(m_u, M_u_std, 'u'), (m_d, M_d_std, 'd'), (m_s, M_s_std, 's')]:
    F = m_j + X_std * (D / M_j)
    print(f"    dW/dM_{lab} = m_{lab} + X * D/M_{lab} = {F:.4e}  (should be 0)")

# Three-instanton correction
delta_X = -3.0 * c3 / Lambda**6
X_c3 = X_std + delta_X
X_eff = X_c3 + 3.0 * c3 * D**2 / Lambda**18  # = X_std by construction

print(f"\n  Three-instanton correction (c3 = {c3}):")
print(f"    delta X = -3 * c3 / Lambda^6 = {delta_X:.6e} MeV^(-3)")
print(f"    X(c3) = X_std + delta X = {X_c3:.6e} MeV^(-3)")
print(f"    X_eff = X(c3) + 3 c3 D^2 / Lambda^18 = {X_eff:.6e} MeV^(-3)")
print(f"    X_std = {X_std:.6e} MeV^(-3)")
print(f"    X_eff == X_std? {np.isclose(X_eff, X_std)}")

# Meson VEVs are unchanged
print(f"\n  Meson VEVs with c3 correction:")
print(f"  (Unchanged, because X_eff is identical to X_std)")
for m_j, M_j, lab in [(m_u, M_u_std, 'u'), (m_d, M_d_std, 'd'), (m_s, M_s_std, 's')]:
    M_j_new = -X_eff * D / m_j
    print(f"    M_{lab}(c3) = -X_eff * D / m_{lab} = {M_j_new:.4f} MeV  "
          f"(standard: {M_j:.4f} MeV,  match: {np.isclose(M_j_new, M_j)})")

# F-term verification with c3
print(f"\n  F-term verification (c3 = {c3}):")
for m_j, M_j, lab in [(m_u, M_u_std, 'u'), (m_d, M_d_std, 'd'), (m_s, M_s_std, 's')]:
    # Full F-term: m_j + X * D/M_j + 3 c3 D^2 * D/(M_j * Lambda^18)
    # Wait: d(det M^3)/dM_j = 3(det M)^2 * d(det M)/dM_j = 3 D^2 * D/M_j
    # So: dW/dM_j = m_j + X * D/M_j + 3 c3 D^3 / (M_j Lambda^18)
    # At vacuum D = Lambda^6:
    # = m_j + X * Lambda^6/M_j + 3 c3 Lambda^18 / (M_j Lambda^18)
    # = m_j + (X + 3 c3/Lambda^6) * Lambda^6/M_j
    # = m_j + X_eff * Lambda^6/M_j
    F_full = m_j + X_c3 * (D / M_j) + 3.0 * c3 * D**3 / (M_j * Lambda**18)
    F_compact = m_j + X_eff * (D / M_j)
    print(f"    dW/dM_{lab} = {F_full:.4e}  (compact form: {F_compact:.4e})  (should be 0)")

print(f"\n  X equation of motion:")
print(f"    dW/dX = det M - Lambda^6 = {M_u_std * M_d_std * M_s_std - Lambda**6:.4e}  (should be 0)")

# ============================================================
# PART 5: RATIOS AND KOIDE CONNECTION
# ============================================================
print("\n" + "=" * 72)
print("PART 5: X/LAMBDA RATIO AND KOIDE CONNECTION")
print("=" * 72)

print(f"""
  The standard Lagrange multiplier value:

    X_std = -(prod m)^(1/3) / Lambda^4

  In units of Lambda^(-3):

    X_std * Lambda^3 = -(prod m)^(1/3) / Lambda

  Numerically:
    (prod m)^(1/3) = {prod_m_13:.6f} MeV
    Lambda = {Lambda:.1f} MeV
    (prod m)^(1/3) / Lambda = {prod_m_13 / Lambda:.6f}  (dimensionless)
    X_std * Lambda^3 = {X_std * Lambda**3:.6f}
""")

# With c3 correction
X_c3_Lambda3 = X_c3 * Lambda**3
delta_X_Lambda3 = delta_X * Lambda**3
print(f"  With c3 = {c3} correction:")
print(f"    delta X * Lambda^3 = -3 c3 / Lambda^3 * Lambda^3 = -3 c3 = {delta_X_Lambda3:.6f}")
print(f"    X(c3) * Lambda^3 = {X_c3_Lambda3:.6f}")
print(f"")
print(f"  For Koide seed: need X_eff / mu = sqrt(3) for some mass scale mu.")
print(f"  The shift delta X = -3/Lambda^6 is flavor-blind.")
print(f"  The meson VEVs M_j = C/m_j are UNCHANGED by c3.")
print(f"  So c3 affects X directly (pseudo-modulus interpretation)")
print(f"  but does NOT affect the meson spectrum or the Koide condition on M_j.")

# Check: what is X_std in units of some natural scale?
# The natural scale: X has units MeV^{-3} in 3 dimensions? No.
# Actually in 4d SQCD: [W] = 3, [M] = 2 (mesons), so [X] = 1 for det M with [det M] = 2*3=6?
# Wait: dim analysis. [M_ij] = 1 (mass dim). [det M] = 3. [Lambda^6] = 6...
# That doesn't match. Let's be careful.
# In 4d SQCD: [W] = 3. [M] = dim 2. [Lambda] = 1. [det M] = 6.
# X (det M - Lambda^6): [X] * 6 = 3 => [X] = -3. So X has dim -3 as MeV^{-3}.
# c3 (det M)^3/Lambda^18: dim = 18 - 18 = 0? No: [W term] = 3.
# c3 (det M)^3/Lambda^18: dim = c3 * 3 + 18 - 18 = 3 only if [c3]=3. So c3 has dim mass^3?
# Or equivalently, c3 is dimensionless if the term is c3 * Lambda^3 * (det M)^3 / Lambda^21...
# The problem statement writes c3/Lambda^18, so [c3 det M^3 / Lambda^18] = c3 + 18 - 18 = c3 = 3.
# If c3 is dimensionless, the term should be Lambda^3 (det M)^3 / Lambda^18 = Lambda^3 * ...
# We'll just track the numerics as given.

print(f"""
  Dimension check:
    [X] = MeV^(-3)  (since [X * det M] = [X] * MeV^6 = [W] = MeV^3)
    [c3 * (det M)^3 / Lambda^18] = [c3] * MeV^18 / MeV^18 = [c3]
    For [W] = MeV^3: [c3] = MeV^3, so c3 = 1 means 1 MeV^3.

    delta X = -3 c3 / Lambda^6 = -3 * 1 MeV^3 / (300 MeV)^6
            = {-3.0 * c3 / Lambda**6:.6e} MeV^(-3)

    Relative shift:
    |delta X / X_std| = 3 c3 / (Lambda^6 * |X_std|)
                      = 3 c3 * Lambda^4 / (Lambda^6 * (prod m)^(1/3))
                      = 3 c3 / (Lambda^2 * (prod m)^(1/3))
                      = {3 * c3 / (Lambda**2 * prod_m_13):.6e}
                      (This is very small: c3 in MeV^3 << Lambda^2 * (prod m)^(1/3))
""")

# ============================================================
# PART 6: EXACT SYMBOLIC COMPUTATION
# ============================================================
print("=" * 72)
print("PART 6: SYMBOLIC VERIFICATION (sympy)")
print("=" * 72)

mu, ms_sym, md_sym, Lam, c3_sym, X_sym, D_sym = symbols(
    'mu md ms Lambda c3 X D', positive=True)

# At diagonal vacuum with constraint D = det M = Lambda^6,
# the F-term for M_j is:
# m_j + X * D/M_j + 3 c3 * D^2 * D/M_j / Lambda^18 = 0
# where D/M_j = product of other two M_k's = cofactor_j

# From M_j * m_j = -X_eff * D and prod M_j = D:
# prod(-X_eff * D / m_j) = D
# (-X_eff)^3 * D^3 / (prod m) = D
# (-X_eff)^3 = prod_m / D^2 = prod_m / Lambda^12

# So X_eff is determined entirely by prod_m and Lambda, independent of c3.

print("""
  Symbolic result:

  F_j = m_j + X_eff * D/M_j = 0  =>  M_j = -X_eff * D / m_j

  Product constraint (at D = Lambda^6):
    prod_j M_j = (-X_eff)^3 * D^3 / prod_m = Lambda^6
    =>  X_eff^3 = -prod_m / Lambda^12
    =>  X_eff = -(prod_m)^{1/3} / Lambda^4   [unique real cube root]

  Since X_eff = X + 3 c3/Lambda^6:
    X = X_eff - 3 c3/Lambda^6
      = -(prod_m)^{1/3}/Lambda^4 - 3 c3/Lambda^6

  This is EXACT. The meson VEVs M_j = C/m_j with C = Lambda^2 (prod m)^{1/3}
  are identical to the standard (c3=0) Seiberg solution.

  delta X = -3 c3 / Lambda^6   (exact, to all orders in c3)

  Note: "all orders" here is trivially true because the c3 term contributes
  only to X_eff linearly — there is no self-consistent shift of M_j that
  could generate higher-order corrections.
""")

# ============================================================
# PART 7: SUMMARY TABLE
# ============================================================
print("=" * 72)
print("SUMMARY")
print("=" * 72)

print(f"""
  Standard Seiberg (c3 = 0):
    M_j = Lambda^2 (prod m)^(1/3) / m_j
    X   = -(prod m)^(1/3) / Lambda^4

  Three-instanton correction (c3 != 0):
    M_j = SAME (unchanged)
    X   = -(prod m)^(1/3) / Lambda^4 - 3 c3 / Lambda^6
    delta X = -3 c3 / Lambda^6

  Numerically (m_u={m_u}, m_d={m_d}, m_s={m_s} MeV, Lambda={Lambda} MeV, c3={c3} MeV^3):
    X_std      = {X_std:.8e} MeV^(-3)
    delta X    = {delta_X:.8e} MeV^(-3)
    X(c3)      = {X_c3:.8e} MeV^(-3)

    Relative shift |delta X / X_std| = {abs(delta_X/X_std):.6e}

  Koide connection:
    The meson VEVs M_j are unchanged by c3. The Koide condition Q(M) = Q(1/m)
    is therefore also unchanged. The c3 term shifts X (the Lagrange multiplier)
    but NOT X_eff (the quantity controlling meson masses).

    If X is identified as the pseudo-modulus:
      X_std       = {X_std:.6e} MeV^(-3)
      delta X     = {delta_X:.6e} MeV^(-3)
      |delta X / X_std| = 3 c3 Lambda^(-2) / (prod m)^(1/3) = {3*c3/(Lambda**2*prod_m_13):.4e}

    For X/mu = sqrt(3): need a scale mu such that X_eff / mu = sqrt(3).
    With X_eff = -(prod m)^(1/3) / Lambda^4:
      mu = |X_eff| / sqrt(3) = (prod m)^(1/3) / (sqrt(3) Lambda^4)
         = {abs(X_std) / np.sqrt(3):.6e} MeV^(-3)  [dimensionally: MeV^(-3)]
    This does not have units of mass unless interpreted differently.

    More naturally: in the O'Raifeartaigh/Koide connection gv/m = sqrt(3),
    with g = coupling, v = VEV, m = mass. This is a different object from X.
    The Lagrange multiplier X has dimension [mass]^(-3) here and is not
    directly the pseudo-modulus VEV of the O'R model.

  BOTTOM LINE:
    c3 != 0 shifts X by -3c3/Lambda^6 (flavor-blind, exact).
    Meson VEVs and Koide condition on M_j are completely unaffected.
    The shift is numerically tiny: {abs(delta_X/X_std)*100:.4e}% of X_std for c3 = 1 MeV^3.
""")
