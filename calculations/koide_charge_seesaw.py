"""
Analysis of the z_0 + z_i charge decomposition under the Seiberg seesaw.

Direct Koide: z_k = sqrt(m_k), energy balance <z_i^2> = z_0^2  =>  Q = 2/3
Inverse Koide: xi_k = 1/sqrt(m_k), energy balance <xi_i^2> = xi_0^2  =>  Q = 2/3

Questions:
  1. Does the seesaw z -> 1/z preserve the Koide parametrization?
  2. How do the two z_0+z_i decompositions (direct and inverse) relate?
  3. What about the seed self-duality?
"""

import math
import numpy as np

# PDG 2024 MSbar masses in MeV
masses = {'u': 2.16, 'd': 4.67, 's': 93.4, 'c': 1270.0, 'b': 4183.0, 't': 162500.0}

def koide_Q(y):
    """Q from three charge values y_k (mass = y_k^2)."""
    s = sum(y)
    if abs(s) < 1e-30: return float('nan')
    return sum(yi**2 for yi in y) / s**2

def decompose(y):
    """Decompose y_k = y0 + yi, return y0, perturbations, energy balance ratio."""
    y0 = sum(y) / 3
    yi = [yk - y0 for yk in y]
    var = sum(x**2 for x in yi) / 3
    return y0, yi, var / y0**2 if abs(y0) > 1e-30 else float('nan')

print("=" * 75)
print("1. CHARGE DECOMPOSITIONS FOR THE THREE INDEPENDENT KOIDE TRIPLES")
print("=" * 75)

triples = [
    ("(-s, c, b)", [-math.sqrt(masses['s']), math.sqrt(masses['c']), math.sqrt(masses['b'])]),
    ("(c, b, t)",  [math.sqrt(masses['c']), math.sqrt(masses['b']), math.sqrt(masses['t'])]),
    ("(1/d, 1/s, 1/b)", [1/math.sqrt(masses['d']), 1/math.sqrt(masses['s']), 1/math.sqrt(masses['b'])]),
]

for name, y in triples:
    Q = koide_Q(y)
    y0, yi, ratio = decompose(y)
    print(f"\n{name}:")
    print(f"  charges:  [{y[0]:.5f}, {y[1]:.5f}, {y[2]:.5f}]")
    print(f"  z0 = {y0:.5f}")
    print(f"  z_i = [{yi[0]:.5f}, {yi[1]:.5f}, {yi[2]:.5f}]")
    print(f"  <z_i^2> = {sum(x**2 for x in yi)/3:.5f}")
    print(f"  z0^2    = {y0**2:.5f}")
    print(f"  <z_i^2>/z0^2 = {ratio:.6f}  (1.000 = exact energy balance)")
    print(f"  Q = {Q:.6f}  (= 1/3 + ratio/3 = {1/3+ratio/3:.6f})")

print("\n" + "=" * 75)
print("2. SEESAW MAP z -> 1/z ON THE DIRECT TRIPLE (-s, c, b)")
print("=" * 75)

z_direct = [-math.sqrt(masses['s']), math.sqrt(masses['c']), math.sqrt(masses['b'])]
z_inv = [1/z for z in z_direct]

print(f"\nDirect triple (-s,c,b):  z = [{z_direct[0]:.4f}, {z_direct[1]:.4f}, {z_direct[2]:.4f}]")
print(f"Q_direct = {koide_Q(z_direct):.6f}")

print(f"\nInverted 1/z:  xi = [{z_inv[0]:.6f}, {z_inv[1]:.6f}, {z_inv[2]:.6f}]")
print(f"Q(1/z) = {koide_Q(z_inv):.6f}")
print(f"Note: Q(-s,c,b) = {koide_Q(z_direct):.6f} != Q(1/z) = {koide_Q(z_inv):.6f}")
print("=> Seesaw BREAKS the Koide condition for non-seed triples.")

print(f"\nAttempt Koide parametrization of the inverted triple:")
xi0, xi_i, xi_ratio = decompose(z_inv)
print(f"  xi0 = {xi0:.6f}")
print(f"  xi_i = [{xi_i[0]:.6f}, {xi_i[1]:.6f}, {xi_i[2]:.6f}]")
print(f"  <xi_i^2>/xi0^2 = {xi_ratio:.6f}")
print(f"  This corresponds to Q = {1/3+xi_ratio/3:.6f}")

print("\n" + "=" * 75)
print("3. SEED SELF-DUALITY: Q(0,a,b) = Q(0,1/a,1/b)")
print("=" * 75)

ms, mc = math.sqrt(masses['s']), math.sqrt(masses['c'])
seed = [0, ms, mc]
seed_inv = [0, 1/ms, 1/mc]  # note: 1/0 -> 0 (already zero)

print(f"\nSeed (0, sqrt(s), sqrt(c)) = [0, {ms:.4f}, {mc:.4f}]")
print(f"Q_seed = {koide_Q(seed):.6f}")
print(f"\nInverted seed (0, 1/sqrt(s), 1/sqrt(c)) = [0, {1/ms:.6f}, {1/mc:.6f}]")
print(f"Q_seed_inv = {koide_Q(seed_inv):.6f}")
print(f"\nIdentical: Q(0,a,b) = Q(0,1/a,1/b) by algebraic identity Q(a,b) = Q(1/a,1/b).")
print(f"The seed is a FIXED POINT of the seesaw map.")

print("\n" + "=" * 75)
print("4. PERTURBATIVE ANALYSIS: seesaw preserves energy balance when z_i << z_0")
print("=" * 75)

print("""
For z_k = z_0 + z_i with z_i << z_0:
  1/z_k = 1/(z_0 + z_i) ≈ (1/z_0)(1 - z_i/z_0 + z_i^2/z_0^2 - ...)

Dual charge: xi_k = 1/z_k ≈ xi_0 - z_i/z_0^2 + ...
  where xi_0 = 1/z_0

Dual perturbation: xi_i = -z_i/z_0^2 = -z_i * xi_0^2 / (1/z_0) = -z_i * xi_0/z_0

So: <xi_i^2> ≈ (xi_0/z_0)^2 * <z_i^2> = xi_0^2 * (<z_i^2>/z_0^2)

If <z_i^2> = z_0^2 (exact energy balance), then <xi_i^2> ≈ xi_0^2.
CONCLUSION: Seesaw preserves energy balance perturbatively (z_i << z_0).
            But real triples have z_i ~ z_0, so corrections matter.
""")

# Check: how "perturbative" are the actual triples?
for name, y in triples:
    y0, yi, ratio = decompose(y)
    max_pert = max(abs(x) for x in yi) / abs(y0) if abs(y0) > 0 else float('inf')
    print(f"  {name:20s}: max|z_i|/z_0 = {max_pert:.3f}  {'(perturbative)' if max_pert < 0.5 else '(NON-perturbative)'}")

print("\n" + "=" * 75)
print("5. RIGOROUS PROOF: Koide parametrization not preserved under z -> 1/z")
print("=" * 75)

print("""
Claim: If z_k = z_0(1 + sqrt(2) cos(delta + 2*pi*k/3)), then 1/z_k CANNOT
       have the form xi_0(1 + sqrt(2) cos(delta' + 2*pi*k/3)) for any xi_0, delta',
       UNLESS one z_k = 0.

Proof: z_k * xi_k = 1 for all k.  Product:
  z_0 * xi_0 * (1 + sqrt(2) cos theta_k)(1 + sqrt(2) cos phi_k) = 1

Let psi_k = (delta+delta')/2 + 2*pi*k/3 and A = 2*sqrt(2)*cos((delta-delta')/2).
The k-dependent part becomes:  2*cos^2(psi_k) + A*cos(psi_k) - 1 = 0.

This quadratic in cos(psi_k) can have at most 2 roots.
But psi_0, psi_1, psi_2 differ by 2*pi/3, giving 3 distinct values of cos(psi_k)
(generically). So the equation cannot hold for all three k.

Exception: when cos(psi_1) = cos(psi_2) = -1/2, cos(psi_0) = 1.
This requires psi_0 = 0 mod 2*pi, hence delta' = -delta.
Then A = 2*sqrt(2)*cos(delta), and the two conditions give A = -1,
so cos(delta) = -1/(2*sqrt(2)).  But this implies z_k = 0 for one k
(the seed angle). QED: non-seed triples are not self-dual.
""")

# Verify numerically: check if any non-trivial angle works
print("Numerical verification:")
for delta in np.linspace(0.1, 2*np.pi-0.1, 100):
    z = [1 + np.sqrt(2)*np.cos(delta + 2*np.pi*k/3) for k in range(3)]
    if min(abs(zk) for zk in z) < 0.01:  # skip near-zero (seed)
        continue
    xi = [1/zk for zk in z]
    Q_z = koide_Q(z)
    Q_xi = koide_Q(xi)
    if abs(Q_z - 2/3) < 0.001 and abs(Q_xi - 2/3) < 0.01:
        print(f"  delta={delta:.3f}: Q(z)={Q_z:.6f}, Q(1/z)={Q_xi:.6f}, diff={abs(Q_z-Q_xi):.6f}")
        break
else:
    print("  No non-seed delta found where BOTH Q(z) and Q(1/z) are near 2/3.")
    print("  (Q(z)=2/3 always, but Q(1/z) generically != 2/3)")

print("\n  Sample: delta=2.5 (near physical (-s,c,b)):")
delta = 2.5
z = [1 + np.sqrt(2)*np.cos(delta + 2*np.pi*k/3) for k in range(3)]
xi = [1/zk for zk in z]
print(f"    z = [{z[0]:.4f}, {z[1]:.4f}, {z[2]:.4f}]")
print(f"    Q(z) = {koide_Q(z):.6f}")
print(f"    1/z = [{xi[0]:.4f}, {xi[1]:.4f}, {xi[2]:.4f}]")
print(f"    Q(1/z) = {koide_Q(xi):.6f}")

print("\n" + "=" * 75)
print("6. THE TWO INDEPENDENT ENERGY-BALANCE SECTORS")
print("=" * 75)

print("""
ELECTRIC sector (quark masses, z = sqrt(m)):
  (-s, c, b):    <z_i^2>/z0^2 = {r1:.4f}   (Q = {q1:.4f})
  (c, b, t):     <z_i^2>/z0^2 = {r2:.4f}   (Q = {q2:.4f})

MAGNETIC sector (meson VEVs ∝ 1/m, xi = 1/sqrt(m)):
  (1/d, 1/s, 1/b): <xi_i^2>/xi0^2 = {r3:.4f}   (Q = {q3:.4f})

These are INDEPENDENT energy-balance conditions constraining different quarks.
The electric sector constrains (s,c,b,t). The magnetic sector constrains (d,s,b).
Together: 5 quarks constrained by 3 Koide conditions + v0-doubling = 4 constraints.
Free parameters remaining: m_u, m_d (connected to Cabibbo angle), plus m_s as scale.
""".format(
    r1=decompose(triples[0][1])[2], q1=koide_Q(triples[0][1]),
    r2=decompose(triples[1][1])[2], q2=koide_Q(triples[1][1]),
    r3=decompose(triples[2][1])[2], q3=koide_Q(triples[2][1]),
))

print("=" * 75)
print("7. ALTERNATIVE: DUAL CHARGE AS MESON VEV")
print("=" * 75)

print("""
In the ISS magnetic dual, meson VEVs scale as <M^j_j> = Lambda^2 / m_j.
The "dual charge" is:
  xi_j = sqrt(<M^j_j>) / Lambda = 1/sqrt(m_j)  (in natural units Lambda=1)

So the dual Koide Q(1/d, 1/s, 1/b) = 0.6652 is literally:
  Q(sqrt(<M_d>), sqrt(<M_s>), sqrt(<M_b>)) = 0.6652

= a Koide condition on the meson condensate charges.

The ELECTRIC theory has energy balance on quark masses.
The MAGNETIC theory has energy balance on meson VEVs.
The seed is self-dual: the same condition in both descriptions.
The bloom breaks duality: different Q values in the two descriptions.
""")
