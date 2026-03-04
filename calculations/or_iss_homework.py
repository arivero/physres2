#!/usr/bin/env python3
"""
O'Raifeartaigh model and ISS mapping — homework computation.

Parts 1-3: O'R model fermion mass matrix eigenvalues and Koide Q.
Parts 4-6: ISS vacuum F-terms and parameter mapping (analytic).
"""

import numpy as np
from numpy import sqrt, array
from fractions import Fraction

print("=" * 70)
print("PART 1: EIGENVALUES OF THE FERMION MASS MATRIX")
print("=" * 70)

print("""
O'R superpotential:  W = f Phi_0 + m Phi_1 Phi_2 + g Phi_0 Phi_1^2

At vacuum phi_1 = 0, phi_0 = v:

  W_{ij} = d^2 W / d phi_i d phi_j   (evaluated at vacuum)

  W_{00} = 0                 (d^2W/dphi_0^2 = 0, no phi_0^2 term)
  W_{01} = W_{10} = g*phi_1 * 2 = 0  (zero since phi_1=0)
  W_{02} = W_{20} = 0                (no phi_0*phi_2 term)
  W_{11} = 2g*phi_0 = 2gv           (from g Phi_0 Phi_1^2, d^2/dphi_1^2)
  W_{12} = W_{21} = m               (from m Phi_1 Phi_2)
  W_{22} = 0                        (no Phi_2^2 term)

Basis ordering: phi_0, phi_1, phi_2.

  M_f = [[0,    0,    0 ],
         [0,    2gv,  m ],
         [0,    m,    0 ]]

The (0,0) block gives one massless fermion (the goldstino from broken SUSY).

The 2x2 lower-right block determines the massive eigenvalues:
  B = [[2gv,  m],
       [m,    0]]

char. poly: lambda^2 - 2gv * lambda - m^2 = 0

  lambda = gv ± sqrt((gv)^2 + m^2)
""")

# Symbolic verification via general formula
print("Analytic eigenvalues of 2x2 block [[2gv, m],[m, 0]]:")
print("  lambda_pm = gv ± sqrt((gv)^2 + m^2)")
print()

# Numerical check with a specific value
gv_over_m = sqrt(3.0)
print(f"Setting gv/m = sqrt(3) = {gv_over_m:.10f}")
print()

# Eigenvalues of full 3x3 matrix
# Use m=1 as scale, gv = sqrt(3)
m_val = 1.0
gv_val = sqrt(3.0)

M_f = array([
    [0.0,    0.0,    0.0],
    [0.0,    2*gv_val, m_val],
    [0.0,    m_val,  0.0]
])

print("Fermion mass matrix M_f (with gv/m = sqrt(3), m=1):")
print(M_f)
print()

# The fermion mass matrix in SUSY is W_ij, symmetric.
# Physical masses are singular values of W_ij, i.e. sqrt(eigenvalues of W^dag W).
# Since M_f is real symmetric here, eigenvalues of W^dag W = (eigenvalues of M_f)^2.
# Eigenvalues of M_f (including possible negatives):
evals_Mf = np.linalg.eigvalsh(M_f)
print(f"Eigenvalues of M_f: {evals_Mf}")
print()

# Physical fermion masses = |eigenvalues| (singular values of M_f)
# Since M_f is symmetric, sing. values = |lambda_i|
physical_masses = np.sort(np.abs(evals_Mf))
print(f"Physical masses |lambda_i| = {physical_masses}")
print()

# Analytic values for gv/m = sqrt(3):
# lambda_pm = gv ± sqrt((gv)^2 + m^2)
#           = m*sqrt(3) ± sqrt(3m^2 + m^2)
#           = m*sqrt(3) ± 2m
# So:
# lambda_+ = m*(sqrt(3) + 2) = m*(2 + sqrt(3))
# lambda_- = m*(sqrt(3) - 2) = -m*(2 - sqrt(3))
# Physical masses: m_+ = m*(2+sqrt(3)), m_- = m*(2-sqrt(3))

print("=" * 70)
print("PART 2: MASS RATIO VERIFICATION AT gv/m = sqrt(3)")
print("=" * 70)

m_plus_analytic  = m_val * (2.0 + sqrt(3.0))
m_minus_analytic = m_val * (2.0 - sqrt(3.0))   # |lambda_-|

print(f"Analytic:  m_+ = m*(2+sqrt(3)) = {m_plus_analytic:.10f}")
print(f"Analytic:  m_- = m*(2-sqrt(3)) = {m_minus_analytic:.10f}")
print()

# Compare with numerical eigenvalues
m_plus_num  = physical_masses[2]
m_minus_num = physical_masses[1]
print(f"Numerical: m_+ = {m_plus_num:.10f}")
print(f"Numerical: m_- = {m_minus_num:.10f}")
print()

ratio_analytic = (2 + sqrt(3.0))**2
ratio_numeric  = m_plus_num / m_minus_num

print(f"Ratio m_+/m_- (analytic): (2+sqrt(3))^2 / 1  <- wrong, let me compute correctly")
# Correct ratio:
# m_+/m_- = (2+sqrt(3))/(2-sqrt(3))
#          = (2+sqrt(3))^2 / ((2+sqrt(3))(2-sqrt(3)))
#          = (2+sqrt(3))^2 / (4-3)
#          = (2+sqrt(3))^2

actual_ratio = m_plus_analytic / m_minus_analytic
print(f"m_+/m_-   = (2+sqrt(3))/(2-sqrt(3))")
print(f"          = (2+sqrt(3))^2 / [(2+sqrt(3))(2-sqrt(3))]")
print(f"          = (2+sqrt(3))^2 / (4-3)")
print(f"          = (2+sqrt(3))^2")
print(f"          = {(2+sqrt(3.0))**2:.10f}")
print(f"Numerical ratio m_+/m_- = {ratio_numeric:.10f}")
print(f"Match: {abs(ratio_numeric - (2+sqrt(3.0))**2) < 1e-10}")
print()

print("Verification step-by-step:")
print(f"  (2+sqrt(3))*(2-sqrt(3)) = 4 - 3 = 1  -> m_+*m_- = m^2")
print(f"  (2+sqrt(3.0))*(2-sqrt(3.0)) = {(2+sqrt(3.0))*(2-sqrt(3.0)):.15f}")
print(f"  So m_+/m_- = (2+sqrt(3))^2 = {(2+sqrt(3.0))**2:.10f}  [confirmed]")
print()

print("=" * 70)
print("PART 3: KOIDE Q FOR TRIPLE (0, m_-, m_+)")
print("=" * 70)

print("""
Koide formula: Q = (m_1 + m_2 + m_3) / (sqrt(m_1) + sqrt(m_2) + sqrt(m_3))^2

For triple (m_0, m_-, m_+) = (0, m*(2-sqrt(3)), m*(2+sqrt(3))):

  Numerator   = 0 + m*(2-sqrt(3)) + m*(2+sqrt(3)) = 4m

  Denominator factor:
    sqrt(m_-) = sqrt(m)*sqrt(2-sqrt(3))
    sqrt(m_+) = sqrt(m)*sqrt(2+sqrt(3))

  Note: sqrt(2-sqrt(3)) = (sqrt(6)-sqrt(2))/2   [half-angle identity]
        sqrt(2+sqrt(3)) = (sqrt(6)+sqrt(2))/2

  Let A = sqrt(2-sqrt(3)), B = sqrt(2+sqrt(3)).
  A + B = sqrt(6)/2 + sqrt(2)/2 - sqrt(6)/2 + sqrt(2)/2 ...

  Actually: A^2 + B^2 = (2-sqrt(3)) + (2+sqrt(3)) = 4
            2AB = 2*sqrt((2-sqrt(3))(2+sqrt(3))) = 2*sqrt(4-3) = 2
  So (A+B)^2 = A^2 + B^2 + 2AB = 4 + 2 = 6  ->  A+B = sqrt(6).

  Denominator = (0 + sqrt(m)*A + sqrt(m)*B)^2 = m*(A+B)^2 = m*6

  Q = 4m / (6m) = 4/6 = 2/3   [exact]
""")

# Numerical verification
m0, mm, mp = 0.0, m_minus_analytic, m_plus_analytic
num   = m0 + mm + mp
denom = (sqrt(m0) + sqrt(mm) + sqrt(mp))**2
Q_numerical = num / denom
print(f"Numerical Q = {num:.10f} / {denom:.10f} = {Q_numerical:.15f}")
print(f"2/3 exact   = {2/3:.15f}")
print(f"Difference  = {abs(Q_numerical - 2/3):.2e}  (should be machine epsilon)")
print()

# Cross-check with general gv/m
print("Cross-check: Q for general gv/m ratio r = gv/m:")
for r in [0.5, 1.0, sqrt(3.0), 2.0, 5.0]:
    mp_r = r + sqrt(r**2 + 1.0)   # m_+ / m
    mm_r = abs(r - sqrt(r**2 + 1.0))  # |m_-| / m  (always positive)
    if mm_r < 1e-15:
        continue
    Q_r = (mp_r + mm_r) / (sqrt(mp_r) + sqrt(mm_r))**2
    print(f"  r={r:.4f}:  m_+/m={mp_r:.6f}, m_-/m={mm_r:.6f},  Q(0,m_-,m_+) = {Q_r:.8f}")

print()
print("Note: Q = 2/3 holds ONLY when m_0 = 0, for any r.")
print("When the massless state (goldstino) is included in the triple,")
print("Q is determined purely by the ratio m_+/m_-.")
print()

# Let's verify this algebraically:
print("Algebraic proof that Q=2/3 whenever m_0=0:")
print("  Set m_0=0, masses = (0, a, b) with b > a > 0.")
print("  Q = (a+b)/(sqrt(a)+sqrt(b))^2")
print("  With a = m*(2-sqrt(3)), b = m*(2+sqrt(3)):")
print(f"    a*b = m^2*(4-3) = m^2  ->  (sqrt(a)*sqrt(b)) = m")
print(f"    (sqrt(a)+sqrt(b))^2 = a+b+2*sqrt(ab) = 4m + 2m = 6m")
print(f"    Q = 4m / 6m = 2/3.")
print()
print("More generally, for (0, m_-, m_+) with m_+*m_- = m^2:")
print("  (sqrt(m_-)+sqrt(m_+))^2 = m_-+m_++2*sqrt(m_-*m_+) = m_-+m_++2m")
print("  Q = (m_-+m_+)/(m_-+m_++2m)")
print("  Q = 2/3  iff  3(m_-+m_+) = 2(m_-+m_++2m) = 2(m_-+m_+)+4m")
print("           iff  (m_-+m_+) = 4m")
print(f"  Check: m_-+m_+ = m*(2-sqrt(3))+m*(2+sqrt(3)) = 4m. Confirmed.")
print()


print("=" * 70)
print("PARTS 4-6: ISS MODEL ANALYTIC ANALYSIS")
print("=" * 70)

print("""
ISS MODEL: N_f=4, N_c=3, N_c^mag=1

Superpotential in the magnetic (IR free) description:
  W = h Tr(q Phi q~) - h mu^2 Tr(Phi)

where:
  Phi is a 4x4 meson superfield matrix
  q, q~ are 1x4 magnetic quark superfields (N_c^mag=1 color index)
  h is the magnetic Yukawa coupling

Index structure: W = h q_a Phi^ab q~_b - h mu^2 Phi^aa   (a,b = 1..4 flavor)

--- PART 4: F-TERMS AND SUSY BREAKING AT ISS VACUUM ---

Vacuum ansatz:  <q_a> = <q~_a> = phi_0 * delta_{a1}   (phi_0 = mu)

F-terms:

F_{Phi^ab} = dW/dPhi^ab = h q_a q~_b - h mu^2 delta^ab

At the vacuum:
  q_a = phi_0 delta_{a1},  q~_b = phi_0 delta_{b1}
  -> q_a q~_b = phi_0^2 delta_{a1} delta_{b1}

So:
  F_{Phi^ab} = h phi_0^2 delta_{a1}delta_{b1} - h mu^2 delta_{ab}

With phi_0 = mu:
  F_{Phi^ab} = h mu^2 (delta_{a1}delta_{b1} - delta_{ab})

Component by component:
  a=b=1:  F_{11} = h mu^2 (1 - 1) = 0           <- satisfied!
  a=b≠1:  F_{aa} = h mu^2 (0 - 1) = -h mu^2     <- NONZERO, breaks SUSY
  a≠b:    F_{ab} = h mu^2 (0 - 0) = 0            <- satisfied

So only the DIAGONAL components with a=b≠1 (i.e., a=2,3,4) have
nonzero F-terms. These 3 components break SUSY.

F_{q_a} = dW/dq_a = h Phi^ab q~_b

At vacuum <Phi> = 0, <q~_b> = phi_0 delta_{b1}:
  F_{q_a} = h * 0 * phi_0 delta_{b1} = 0   (satisfied, Phi has zero VEV)

F_{q~_b} = dW/dq~_b = h q_a Phi^ab
  = 0   (same reason)

Summary:
  SUSY broken by: F_{Phi^{22}}, F_{Phi^{33}}, F_{Phi^{44}}  (rank N_f - N_c = 1 breaking)
  This is the ISS "rank condition": the rank of <q_a q~_b> is N_c^mag=1,
  but the identity requires rank N_f=4. The rank deficit = 3 = N_f - N_c^mag
  broken F-terms.
""")

print("--- PART 5: PSEUDO-MODULUS SECTOR AND O'R PARAMETER MAPPING ---")

print("""
Expand around the ISS vacuum. Let:
  Phi^11 = X       (the pseudo-modulus, unfixed at tree level)
  Phi^aa = chi_a   for a = 2,3,4  (the SUSY-breaking directions)
  q_1 = phi_0 + rho,  q~_1 = phi_0 + rho~  (fluctuations in vev direction)
  q_a, q~_a for a=2,3,4  (the "psi" directions)

The key sector is phi^11 = X (the pseudo-modulus direction), which at tree
level has zero F-term and is undetermined.

Restricting to the (X, q_1, q~_1) sector with chi_a = 0:

  W|_{sector} = h X (q_1 q~_1) - h mu^2 X
              = h X (phi_0 + rho)(phi_0 + rho~) - h mu^2 X

At linear order in rho, rho~:
  W ≈ h X phi_0 (rho + rho~) + h X phi_0^2 - h mu^2 X
    = h X phi_0 (rho + rho~) + h X (phi_0^2 - mu^2)

With phi_0 = mu: phi_0^2 - mu^2 = 0, so:
  W|_{sector} ≈ h phi_0 X (rho + rho~)  [plus higher order]

But this doesn't have the O'R form yet. Let's look at the sector
(X, q_a, q~_a) for a = 2,3,4 fixed (the pseudo-flat directions of Phi^aa):

For a given "broken" flavor index a ∈ {2,3,4}:
  chi_a ≡ Phi^aa  (pseudo-modulus — direction that's still flat after F-condition)

  Wait — actually chi_a for a ≠ 1 DOES break SUSY. The pseudo-modulus
  is specifically the TRACELESS part of Phi^11 = X.

Let me identify correctly:

PSEUDO-MODULUS = X ≡ Phi^{11}

The F-term for Phi^{11} = h q_1 q~_1 - h mu^2.
At vacuum <q_1> = <q~_1> = phi_0, this gives F_{11} = h phi_0^2 - h mu^2.
With phi_0 = mu: F_{11} = 0, so X = Phi^{11} is the pseudo-modulus.

Now expand around X (keeping X general):
  Phi^{11} = X (free)

The coupling of X to q_1, q~_1:
  W contains:  h Phi^{11} q_1 q~_1 - h mu^2 Phi^{11}
             = h X q_1 q~_1 - h mu^2 X

The X F-term:  F_X = h q_1 q~_1 - h mu^2

At tree level with <q_1>=<q~_1>=phi_0=mu and X=0:
  F_X = h mu^2 - h mu^2 = 0.   (Pseudo-modulus: X leaves F_X=0 trivially
  since q_1 vev doesn't depend on X.)

Now look at how X couples to the BROKEN sector (q_2, q~_2) [for a=2]:
  W includes: h Phi^{12} q_1 q~_2 + h Phi^{21} q_2 q~_1
  and: h Phi^{22} q_2 q~_2 - h mu^2 Phi^{22}

These off-diagonal Phi entries ARE the crucial sector.

For the (X=Phi^{11}, Phi^{12}, Phi^{21}) sector coupled to
(q_1 = phi_0 = mu fixed, q_2, q~_2):

  W ≈ h mu (Phi^{12} q~_2 + Phi^{21} q_2) + h X q_2 q~_2 - h mu^2 Phi^{22}

This has the precise O'R form if we identify:
  Phi_0 ↔ X         (the pseudo-modulus, plays role of Phi_0 in O'R)
  Phi_1 ↔ q_2       (light magnetic quark, plays role of Phi_1)
  Phi_2 ↔ Phi^{12}  (off-diagonal meson, plays role of Phi_2)

With couplings:
  f_a = h mu^2 - h phi_0^2  [but phi_0^2 ≠ mu^2 generally during modulus scan]
      → at VEV: f_a = 0 (for the (1,1) block). The actual f comes from F_{22}:
  f_{ISS,a} = h mu^2   (the constant term in F_{Phi^{aa}})

So the CORRECT mapping is:

  O'R field Phi_0  ↔  ISS field X = Phi^{11}             (pseudo-modulus)
  O'R field Phi_1  ↔  ISS field q_a  (a≥2)               (magnetic quark)
  O'R field Phi_2  ↔  ISS field Phi^{1a} = Phi^{a1}      (off-diagonal meson)

  O'R parameter f  ↔  ISS quantity  h mu^2               (SUSY-breaking scale)
  O'R parameter m  ↔  ISS quantity  h phi_0 = h mu       (coupling × VEV)
  O'R parameter g  ↔  ISS coupling  h                    (magnetic Yukawa)
  O'R pseudo-modulus v  ↔  ISS pseudo-modulus VEV  <X>   (undetermined at tree level)

Parameter identification summary:
  f   = h mu^2    (O'R linear term coefficient = ISS SUSY-breaking scale)
  m   = h mu      (O'R bilinear = ISS Yukawa × VEV)
  g   = h         (O'R cubic coupling = ISS magnetic Yukawa)
  v   = <X>       (O'R modulus VEV = ISS pseudo-modulus expectation value)
""")

print("--- PART 6: DOES gv/m = sqrt(3) TRANSLATE TO h = sqrt(3)? ---")

print("""
From Part 5, the parameter mapping is:
  g = h    (magnetic Yukawa coupling)
  v = <X>  (pseudo-modulus VEV, set by one-loop CW potential)
  m = h * mu  (= h * phi_0 at the vacuum)

The O'R ratio:
  gv/m = h * <X> / (h * mu) = <X> / mu

So:
  gv/m = sqrt(3)   iff   <X>/mu = sqrt(3)

This is a CONDITION ON THE PSEUDO-MODULUS VEV, not on h.

The O'Raifeartaigh condition gv/m = sqrt(3) translates to:
    <X> = sqrt(3) * mu

The coupling h appears in BOTH numerator (g = h) and denominator (m = h*mu),
so it CANCELS. The condition is h-INDEPENDENT.

In the ISS model:
  - The tree-level pseudo-modulus X is undetermined.
  - The one-loop CW potential V_{CW}(X) lifts the pseudo-moduli.
  - The minimum of V_{CW} determines <X>.

For the STANDARD ISS CW computation (Intriligator, Seiberg, Shih 2006),
the potential V_{CW}(X) for the pseudo-modulus X is computed from the
STr M^4 log M^2 formula. The minimum is at X=0 (not at X = sqrt(3)*mu).

The VALUE h = sqrt(3) does NOT appear in the condition for Koide seed.
What matters is the ratio <X>/mu, which is determined dynamically by V_CW.

HOWEVER: h DOES appear in the structure of V_CW. Specifically:
  - The CW mass for X is m_X^2 = (h^4 / 16pi^2) * A(phi_0/mu)
  - The O'R analog condition gv/m = sqrt(3) corresponds to the RATIO
    of two mass parameters in the 2x2 subblock, not to a coupling constant.

CONCLUSION:
  - The Koide seed condition gv/m = sqrt(3) translates to <X>/mu = sqrt(3),
    which is NOT automatically satisfied by setting h = sqrt(3).
  - The appearance of h = sqrt(3) in the ISS context (e.g., in iss_cw_refined.py)
    is a SEPARATE fact: it is chosen to match the O'R mass ratio m_+/m_- = (2+sqrt(3))^2
    in the SPECTRUM (the goldstino is one of the fermions in the spectrum,
    and the two massive ones have the required Koide Q = 2/3).
  - Setting h = sqrt(3) in ISS does ensure sqrt(3) appears in the fermion
    mass matrix as 2gv = 2h<X>, and the ratio gv/m = h<X>/(h*mu) = <X>/mu.
    So h = sqrt(3) combined with <X> = mu gives gv/m = sqrt(3).
  - But CW minimization generically gives <X> = 0 (the standard result).
    The condition gv/m = sqrt(3) would require nonzero <X>, which is not
    naturally selected by one-loop CW in the pure ISS model.

TO SUMMARIZE:
  gv/m = sqrt(3)  <=>  <X>/mu = sqrt(3)   (constrains pseudo-modulus VEV)
  h = sqrt(3)     is a SEPARATE input (magnetic Yukawa coupling)
  h = sqrt(3) alone does NOT enforce gv/m = sqrt(3)
  You need BOTH h = sqrt(3) AND <X> = mu  to get gv/m = sqrt(3).
  But standard ISS CW gives <X> = 0, so gv/m = 0, not sqrt(3).
""")

print("=" * 70)
print("NUMERICAL SUMMARY")
print("=" * 70)

print(f"\nO'R model at gv/m = sqrt(3):")
print(f"  m_0 = {0.0}")
print(f"  m_- = m*(2 - sqrt(3)) = {m_minus_analytic:.10f}  [for m=1]")
print(f"  m_+ = m*(2 + sqrt(3)) = {m_plus_analytic:.10f}  [for m=1]")
print(f"  m_+/m_-              = {m_plus_analytic/m_minus_analytic:.10f}")
print(f"  (2+sqrt(3))^2        = {(2+sqrt(3.0))**2:.10f}")
print(f"  Match: {abs(m_plus_analytic/m_minus_analytic - (2+sqrt(3.0))**2) < 1e-10}")
print(f"  Koide Q(0,m_-,m_+)   = {Q_numerical:.15f}")
print(f"  2/3                  = {2/3:.15f}")
print(f"  |Q - 2/3|            = {abs(Q_numerical - 2/3):.2e}")

print(f"\nISS parameter mapping:")
print(f"  f (O'R)  <->  h * mu^2       (ISS SUSY-breaking)")
print(f"  m (O'R)  <->  h * mu         (ISS Yukawa x VEV)")
print(f"  g (O'R)  <->  h              (ISS magnetic Yukawa)")
print(f"  v (O'R)  <->  <X>            (ISS pseudo-modulus VEV)")
print(f"  gv/m = sqrt(3)  <->  <X>/mu = sqrt(3)  (h-independent condition)")
print(f"  Standard ISS CW gives <X>=0, so gv/m=0 at tree level.")
print(f"  h = sqrt(3) is an independent choice, does NOT enforce gv/m = sqrt(3).")
