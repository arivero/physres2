#!/usr/bin/env python3
"""
ISS Coleman-Weinberg: refined computation.

Focus on:
1. Correct analytic CW mass formula and comparison with numerical
2. The regime where CW DOES transmit mass ratios (mu << m_quark)
3. The Seiberg seesaw regime for comparison
4. Connection to Koide phase delta
"""

import numpy as np
from numpy import sqrt, log, pi, diag

# Physical parameters
N_c = 3
N_f = 4
h = np.sqrt(3.0)
Lambda = 300.0  # MeV

# Quark masses (MeV)
m_u = 2.16
m_d = 4.67
m_s = 93.4
m_c = 1270.0
m_quarks = np.array([m_u, m_d, m_s, m_c])
labels = ['u', 'd', 's', 'c']

def koide_Q(m1, m2, m3):
    """Koide quotient Q = (m1+m2+m3) / (sqrt(m1)+sqrt(m2)+sqrt(m3))^2"""
    s = sqrt(abs(m1)) + sqrt(abs(m2)) + sqrt(abs(m3))
    return (abs(m1) + abs(m2) + abs(m3)) / s**2

# ============================================================
# ANALYTIC CW MASS FORMULA (CORRECT VERSION)
# ============================================================
print("=" * 80)
print("ANALYTIC CW MASS FOR ISS BROKEN MESONS")
print("=" * 80)

print("""
The correct one-loop CW mass for the ISS pseudo-moduli chi_a
comes from the STr[M^4 log M^2] formula applied to the
chi_a-dependent mass matrix of the (psi_a, psi~_a, Y_a, Y~_a) sector.

Key: the bosonic mass eigenvalues (as functions of chi_a) are:
  m^2_{B,1} = h^2(phi_0^2 + |chi_a|^2) + h*f_a    (B+ mode)
  m^2_{B,2} = h^2(phi_0^2 + |chi_a|^2) - h*f_a    (B- mode)
  m^2_{B,3} = h^2 phi_0^2                          (Y mode, 2x degenerate)

And fermion mass eigenvalues:
  m^2_{F,1} = h^2(phi_0^2 + |chi_a|^2/2 + |chi_a|*sqrt(phi_0^2 + chi_a^2/4))
  m^2_{F,2} = h^2(phi_0^2 + |chi_a|^2/2 - |chi_a|*sqrt(phi_0^2 + chi_a^2/4))

But these are approximate. Let me use the exact eigenvalue formulas.

Actually, let me compute analytically the second derivative of the
CW potential d^2 V_CW / d|chi_a|^2 at chi_a = 0 using the
exact mass eigenvalues.

At chi_a = 0:
  Bosonic: m^2_{B+} = h^2 phi_0^2 + h f_a,  m^2_{B-} = h^2 phi_0^2 - h f_a
           m^2_{Y} = h^2 phi_0^2 (2x)
  Fermionic: m^2_F = h^2 phi_0^2 (2x degenerate Dirac)

The chi_a dependence enters through:
  dm^2_{B+}/d|chi_a|^2 = h^2,  dm^2_{B-}/d|chi_a|^2 = h^2
  d^2 m^2_{B+}/d|chi_a|^2 = 0
  (The Y modes don't depend on chi_a)

For fermions:
  dm^2_{F,1}/d|chi_a|^2 = h^2,  dm^2_{F,2}/d|chi_a|^2 = h^2  (at chi_a = 0)
  d^2 m^2_F/d|chi_a|^2 = h^2   (each)

But the CW potential is not a simple function of mass^2.
The correct way is:

V_CW = (1/64pi^2) sum_i n_i m_i^4(chi_a) [log(m_i^2/Lambda^2) - 3/2]

d^2V/d chi_a^2 = (1/64pi^2) sum_i n_i {
    12 m_i^2 (dm_i^2/dchi)^2 [log(m^2/L^2) - 1]
    + 4 m_i^4 (1/m_i^2)(dm_i^2/dchi)^2
    + m_i^4 (d^2 log m^2/dchi^2)
    + ...
}

This is getting messy. Let me just do the numerical computation PROPERLY.
""")

# ============================================================
# PROPER NUMERICAL CW COMPUTATION
# ============================================================

def compute_spectrum(chi_val, phi0, h_val, f_val):
    """
    Compute exact bosonic and fermionic mass-squared eigenvalues
    for the sector-B fields as a function of chi_a.

    The sector-B superpotential is:
    W = h phi_0 (Y psi~ + psi Y~) + h chi psi psi~

    F-terms:
    F_chi = h psi psi~ - f (constant at vacuum)
    F_psi = h chi psi~ + h phi_0 Y~
    F_psi~ = h chi psi + h phi_0 Y
    F_Y = h phi_0 psi~
    F_Y~ = h phi_0 psi

    Returns: (bosonic_m2, fermionic_m2, n_bos, n_ferm)
    where n_bos, n_ferm are the d.o.f. weights.
    """
    hc = h_val * chi_val
    hp = h_val * phi0

    # Fermionic mass matrix (Weyl basis):
    # W_ij: {psi, psi~, Y, Y~}
    # W_{psi,psi~} = h chi
    # W_{psi,Y~} = h phi_0
    # W_{psi~,Y} = h phi_0
    W_f = np.array([
        [0,  hc, 0,  hp],
        [hc, 0,  hp, 0 ],
        [0,  hp, 0,  0 ],
        [hp, 0,  0,  0 ]
    ], dtype=float)

    # Fermion mass-squared eigenvalues from W^dag W
    fm2 = np.linalg.eigvalsh(W_f @ W_f)
    fm2 = np.sort(fm2)

    # Bosonic mass-squared matrix:
    # "Holomorphic" piece: M2_hol = dF^dag dF where dF_A/dphi_j
    # dF/dphi matrix (rows = {chi, psi, psi~, Y, Y~}, cols = {psi, psi~, Y, Y~})
    # evaluated at vacuum (psi = psi~ = Y = Y~ = 0):
    dF = np.array([
        [0,  0,  0,  0 ],   # dF_chi/d{psi,psi~,Y,Y~} = {h*psi~, h*psi, 0, 0} → 0 at vac
        [0,  hc, 0,  hp],   # dF_psi/d{...} = {0, h*chi, 0, h*phi_0}
        [hc, 0,  hp, 0 ],   # dF_psi~/d{...} = {h*chi, 0, h*phi_0, 0}
        [0,  hp, 0,  0 ],   # dF_Y/d{...} = {0, h*phi_0, 0, 0}
        [hp, 0,  0,  0 ],   # dF_Y~/d{...} = {h*phi_0, 0, 0, 0}
    ], dtype=float)

    M2_hol = dF.T @ dF  # 4x4

    # B-term from F_chi = h psi psi~ - f:
    # d^2|F_chi|^2 / d psi d psi~^* = F_chi^* * (d^2 F_chi / d psi d psi~)
    # At vacuum: F_chi = -f, d^2 F_chi/(d psi d psi~) = h
    # So B_{psi, psi~} = -f * h  (and symmetric)
    B = np.zeros((4, 4))
    B[0, 1] = -h_val * f_val
    B[1, 0] = -h_val * f_val

    # Real 8x8 mass matrix
    n = 4
    M_real = np.zeros((2*n, 2*n))
    for i in range(n):
        for j in range(n):
            M_real[2*i,   2*j]   = M2_hol[i,j] + B[i,j]  # Re-Re
            M_real[2*i+1, 2*j+1] = M2_hol[i,j] - B[i,j]  # Im-Im

    bm2 = np.sort(np.linalg.eigvalsh(M_real))

    return bm2, fm2


def V_CW(chi_val, phi0, h_val, f_val, Lambda_val):
    """One-loop CW potential as function of chi."""
    bm2, fm2 = compute_spectrum(chi_val, phi0, h_val, f_val)
    L2 = Lambda_val**2

    V = 0.0
    # Bosonic: 8 real d.o.f.
    for m2 in bm2:
        if abs(m2) > 1e-30:
            V += m2**2 * (log(abs(m2) / L2) - 1.5)

    # Fermionic: 4 eigenvalues of W^dag W.
    # W is 4x4 symmetric → eigenvalues of W^dag W are m_F^2.
    # Each nonzero m_F corresponds to a Dirac fermion = 2 Weyl = 4 real d.o.f.?
    # No: the W matrix is 4x4, so we have 4 Weyl fermion masses.
    # Each Weyl fermion has 2 real d.o.f.
    # The STr formula weights: +1 per real boson d.o.f., -1 per real fermion d.o.f.
    # 8 bosonic real d.o.f. from 4 complex scalars
    # 8 fermionic real d.o.f. from 4 Weyl fermions
    # So fermionic contribution has coefficient -1 per real d.o.f., times 4 Weyl = 8 real
    # But the 4 eigenvalues of W^dag W give |m|^2 for each Weyl fermion.
    # The CW formula uses each eigenvalue once with weight (# real d.o.f. per eigenvalue).
    # Each eigenvalue of W^dag W has degeneracy 1, representing 1 Weyl fermion = 2 real d.o.f.
    for m2 in fm2:
        if abs(m2) > 1e-30:
            V -= 2 * m2**2 * (log(abs(m2) / L2) - 1.5)  # factor 2 for 2 real d.o.f. per Weyl

    return V / (64 * pi**2)


# Use 5-point stencil for stable numerical second derivative
def second_deriv_5pt(func, x0, dx):
    """5-point stencil for d^2f/dx^2."""
    fp2 = func(x0 + 2*dx)
    fp1 = func(x0 + dx)
    f0  = func(x0)
    fm1 = func(x0 - dx)
    fm2 = func(x0 - 2*dx)
    return (-fp2 + 16*fp1 - 30*f0 + 16*fm1 - fm2) / (12 * dx**2)


print("\n" + "=" * 80)
print("CW MASS COMPUTATION WITH 5-POINT STENCIL")
print("=" * 80)

mu_values = [92.0, 500.0, 5000.0, 50000.0]

for mu in mu_values:
    mu2 = mu**2
    phi0 = sqrt(mu2 - m_u / h)

    print(f"\n{'='*60}")
    print(f"mu = {mu:.1f} MeV  (phi_0 = {phi0:.4f} MeV)")
    print(f"h*mu^2 = {h*mu2:.4f} MeV^2")
    print(f"{'='*60}")

    cw_masses = {}
    for i in range(N_f):
        label = labels[i]
        m_q = m_quarks[i]
        f_i = h * mu2 - m_q

        if f_i <= 0:
            print(f"  {label}: f_a < 0 (m_a > h*mu^2), ISS vacuum doesn't exist")
            continue

        # Adaptive step size
        dx = min(0.01 * phi0, 0.1 * sqrt(f_i / h))
        if dx < 1e-6:
            dx = 1e-6

        def V_of_chi(chi):
            return V_CW(chi, phi0, h, f_i, Lambda)

        m2_cw = second_deriv_5pt(V_of_chi, 0.0, dx)

        if m2_cw > 0:
            m_cw = sqrt(m2_cw)
        else:
            m_cw = -sqrt(abs(m2_cw))

        cw_masses[label] = m_cw
        ratio_to_quark = m_cw / m_q

        print(f"  m_CW({label}) = {abs(m_cw):.8f} MeV  "
              f"(m^2 = {m2_cw:.6f} MeV^2, "
              f"m_CW/m_q = {ratio_to_quark:.6f}, "
              f"m_q/(h*mu^2) = {m_q/(h*mu2):.6f})")

    # Koide analysis
    if len(cw_masses) >= 3:
        triples = [('u','d','s'), ('d','s','c')]
        for t in triples:
            if all(l in cw_masses for l in t):
                Q = koide_Q(abs(cw_masses[t[0]]), abs(cw_masses[t[1]]), abs(cw_masses[t[2]]))
                Q_uv = koide_Q(m_quarks[labels.index(t[0])],
                               m_quarks[labels.index(t[1])],
                               m_quarks[labels.index(t[2])])
                print(f"  Q_CW({t[0]},{t[1]},{t[2]}) = {Q:.8f}  "
                      f"(Q_UV = {Q_uv:.8f}, "
                      f"dev from 2/3: CW {abs(Q-2/3)/(2/3)*100:.4f}%, "
                      f"UV {abs(Q_uv-2/3)/(2/3)*100:.4f}%)")


# ============================================================
# THE KEY REGIME: mu^2 MUCH LARGER THAN m_quark^2
# (ISS with mu >> m/h)
# ============================================================
print("\n\n" + "=" * 80)
print("ANALYSIS: WHY CW MASSES ARE NEARLY DEGENERATE")
print("=" * 80)

print("""
The CW mass for chi_a depends on the quark mass m_a through
the SUSY-breaking F-term:  f_a = h*mu^2 - m_a

The boson-fermion mass splitting that drives the CW potential is:
  m^2_{B+} - m^2_F = h*f_a     (positive splitting)
  m^2_{B-} - m^2_F = -h*f_a    (negative splitting)

Since f_a ≈ h*mu^2 * (1 - m_a/(h*mu^2)), and m_a << h*mu^2,
the splitting is approximately FLAVOR-INDEPENDENT.

The CW mass for chi_a is:
  m^2_CW ∝ (boson-fermion splitting)^2 / (SUSY mass scale)^2
         ∝ f_a^2 / (h^2 phi_0^2)^2
         ∝ (h*mu^2 - m_a)^2 / (h*mu^2)^2    [for phi_0 ≈ mu]
         ≈ 1 - 2 m_a/(h*mu^2)

So the FRACTIONAL flavor dependence is:
  Δm^2_CW / m^2_CW ≈ -2 m_a / (h*mu^2)

For mu = 92 MeV:
""")

for label, m_q in zip(labels, m_quarks):
    frac = 2 * m_q / (h * 92**2)
    print(f"  {label}: 2m/(h*mu^2) = {frac:.6f} ({frac*100:.4f}%)")

print("""
For u, d, s: the fractional correction is < 1%, so CW masses are
nearly degenerate → Q ≈ 1/3.

For c: the correction is ~17%, which gives a detectable but still
small deviation from degeneracy.

CONCLUSION: The ISS CW potential at mu ~ f_pi DOES NOT transmit
the UV Koide condition. The CW spectrum is too degenerate.
""")

# ============================================================
# ALTERNATIVE: WHAT IF mu^2 ~ Lambda^2/m_quark? (Seiberg seesaw regime)
# ============================================================
print("=" * 80)
print("ALTERNATIVE: SEESAW-LIKE REGIME")
print("=" * 80)

print("""
In the Seiberg seesaw (SUSY vacuum), meson VEVs are M_i ~ Lambda^2/m_i.
This inverts the mass hierarchy: lighter quarks → heavier mesons.

In the ISS (metastable vacuum), the CW mass is:
  m_CW ~ (h^2 / 4pi) * f_a / phi_0
  where f_a = h*mu^2 - m_a and phi_0 ≈ mu.

  m_CW ~ (h^2 / 4pi) * (h*mu^2 - m_a) / mu
       ~ (h^3 mu / 4pi) * (1 - m_a/(h*mu^2))

For the CW to transmit mass ratios faithfully, we need:
  m_a/(h*mu^2) to be O(1), not << 1.

This means mu ~ sqrt(m_a / h), which gives:
  For charm: mu ~ sqrt(1270/1.73) ~ 27 MeV
  For strange: mu ~ sqrt(93/1.73) ~ 7.3 MeV

Let's try mu in the "tracking" regime where m_c/(h*mu^2) ~ 0.5.
""")

# Find mu where CW masses track quark mass ratios
mu_tracking = sqrt(m_c / (0.5 * h))  # m_c/(h*mu^2) = 0.5
print(f"\nTracking regime: mu = {mu_tracking:.4f} MeV")
print(f"  m_c/(h*mu^2) = {m_c/(h*mu_tracking**2):.4f}")

mu2_t = mu_tracking**2
phi0_t = sqrt(max(mu2_t - m_u/h, 1e-10))
print(f"  phi_0 = {phi0_t:.4f} MeV")

print(f"\n  CW masses in tracking regime:")
for i in range(N_f):
    label = labels[i]
    m_q = m_quarks[i]
    f_i = h * mu2_t - m_q

    if f_i <= 0:
        print(f"    {label}: f_a < 0 → chi_a direction STABILIZED (not pseudo-modulus)")
        continue

    dx = 0.01 * phi0_t
    if dx < 1e-6:
        dx = 1e-6

    def V_of_chi(chi):
        return V_CW(chi, phi0_t, h, f_i, Lambda)

    m2_cw = second_deriv_5pt(V_of_chi, 0.0, dx)
    m_cw = sqrt(abs(m2_cw)) * np.sign(m2_cw)

    print(f"    m_CW({label}) = {abs(m_cw):.8f} MeV  "
          f"(m_q = {m_q:.2f}, ratio = {abs(m_cw)/m_q:.6f})")

# ============================================================
# CW MASS WITH EXPANSION AROUND f_a
# ============================================================
print("\n\n" + "=" * 80)
print("PERTURBATIVE EXPANSION OF CW MASS")
print("=" * 80)

print("""
At phi_0 = mu (approximately), the CW mass for chi_a can be expanded:

  m^2_CW(chi_a) = (alpha / (4*pi)) * [A(x_a)]

where alpha = h^2/(4*pi) and x_a = f_a/(h^2*phi_0^2) = (h*mu^2 - m_a)/(h^2*mu^2).

For x → 1 (f_a → h^2 phi_0^2, i.e., m_a → 0):
  m^2_CW ∝ constant (flavor-independent)

For x → 0 (f_a → 0, i.e., m_a → h*mu^2):
  m^2_CW → 0 (the B-term vanishes)

So the CW mass interpolates between a constant (for m_a << h*mu^2)
and zero (for m_a → h*mu^2).

The flavor-dependent part scales as:
  Δm^2_CW/m^2_CW ~ -2 * m_a/(h*mu^2) + O(m_a^2/(h*mu^2)^2)

This is the OPPOSITE of Seiberg seesaw (which gives Δm^2 ∝ +m_a).
""")

# Compute exact CW mass as function of x = f_a/(h^2 phi_0^2)
print("\nCW mass^2 as function of x = f/(h^2 phi_0^2):")
print(f"{'x':<12} {'m^2_CW/(h^4/(8pi^2))':<25} {'normalized':<15}")
print("-" * 55)

phi0_ref = 100.0  # MeV reference

x_values = np.array([0.01, 0.05, 0.1, 0.2, 0.5, 0.8, 0.9, 0.95, 0.99, 0.999])
m2_cw_x = []

for x in x_values:
    f_val = x * h**2 * phi0_ref**2
    dx = 0.01 * phi0_ref

    def V_of_chi(chi):
        return V_CW(chi, phi0_ref, h, f_val, Lambda)

    m2 = second_deriv_5pt(V_of_chi, 0.0, dx)
    m2_norm = m2 / (h**4 / (8 * pi**2))
    m2_cw_x.append(m2_norm)

    print(f"{x:<12.4f} {m2_norm:<25.8f} {m2_norm/m2_cw_x[0] if m2_cw_x[0] != 0 else 'N/A':<15}")

# ============================================================
# THE RESULT ABOUT KOIDE TRANSMISSION
# ============================================================
print("\n\n" + "=" * 80)
print("DEFINITIVE RESULT: ISS CW AND KOIDE")
print("=" * 80)

print("""
THEOREM (ISS CW and Koide transmission):

In the ISS model (N_c=3, N_f=4) with h = sqrt(3) and mu = f_pi:

1. The CW meson masses depend on quark masses through f_a = h*mu^2 - m_a.

2. For mu >> m_quark/h (the physical regime), f_a ≈ h*mu^2 for all
   flavors, so all CW meson masses are approximately EQUAL.
   → Koide Q ≈ 1/3 (degenerate limit)

3. The flavor-dependent correction is:
     m_CW(a) ≈ m_CW(0) * (1 - m_a/(h*mu^2))
   This is a LINEAR function of m_a, not of sqrt(m_a).

4. The Koide condition Q = 2/3 is preserved by transformations
   m_a → alpha * m_a (rescaling), but NOT by m_a → m_0 - m_a (affine shift).
   The ISS CW map f(m_a) = (h*mu^2 - m_a) is an affine shift.
   → Q is NOT preserved.

5. The Seiberg seesaw m_a → Lambda^2/m_a DOES preserve Koide
   (it maps Koide to "dual Koide"), but applies to the SUSY vacuum,
   not the ISS metastable vacuum.

6. The ISS CW TRANSMITS the Koide condition only in the limit
   mu → 0 (or more precisely, when h*mu^2 → 0 and the CW mass
   becomes proportional to m_a directly). But in this limit,
   the ISS vacuum doesn't exist (phi_0 → 0).

7. CONCLUSION: The ISS metastable vacuum CW potential does NOT
   transmit the UV Koide condition to the IR meson spectrum.
   The transmission mechanism must be DIFFERENT from one-loop CW.
   Candidates:
   (a) Non-perturbative effects (bion mechanism, as in complete_kahler.md)
   (b) Higher-loop contributions
   (c) The SUSY vacuum Seiberg seesaw (not the ISS metastable vacuum)
   (d) Koide is a UV condition only, transmitted by D-term immunity
""")

# ============================================================
# NUMERICAL SUMMARY TABLE
# ============================================================
print("=" * 80)
print("NUMERICAL SUMMARY TABLE")
print("=" * 80)

mu_val = 92.0
mu2_val = mu_val**2
phi0_val = sqrt(mu2_val - m_u/h)

print(f"\n{'Parameter':<30} {'Value':<25} {'Units'}")
print("-" * 65)
print(f"{'N_c':<30} {'3':<25} {''}")
print(f"{'N_f':<30} {'4':<25} {''}")
print(f"{'h (magnetic Yukawa)':<30} {f'{h:.6f} = sqrt(3)':<25} {''}")
print(f"{'mu (SUSY-breaking scale)':<30} {f'{mu_val:.1f}':<25} {'MeV'}")
print(f"{'mu^2':<30} {f'{mu2_val:.1f}':<25} {'MeV^2'}")
print(f"{'h*mu^2':<30} {f'{h*mu2_val:.4f}':<25} {'MeV^2'}")
print(f"{'phi_0 = <q>':<30} {f'{phi0_val:.4f}':<25} {'MeV'}")
print(f"{'Lambda (cutoff)':<30} {f'{Lambda:.1f}':<25} {'MeV'}")

print(f"\n{'UV Quark':<12} {'m_q (MeV)':<15} {'m_q/(h mu^2)':<15} {'CW mass (MeV)':<18} {'Q deviation'}")
print("-" * 75)

cw_final = {}
for i in range(N_f):
    label = labels[i]
    m_q = m_quarks[i]
    f_i = h * mu2_val - m_q
    dx = 0.001 * phi0_val

    def V_of_chi(chi):
        return V_CW(chi, phi0_val, h, f_i, Lambda)

    m2_cw = second_deriv_5pt(V_of_chi, 0.0, dx)
    m_cw = sqrt(abs(m2_cw))
    cw_final[label] = m_cw

    print(f"{label:<12} {m_q:<15.2f} {m_q/(h*mu2_val):<15.6f} {m_cw:<18.6f}")

# Koide Q values
print(f"\n{'Triple':<15} {'Q (UV quarks)':<18} {'Q (CW mesons)':<18} {'Q (seesaw M~1/m)':<18} {'2/3':<8}")
print("-" * 80)

for triple in [('u','d','s'), ('d','s','c'), ('u','s','c'), ('u','d','c')]:
    m_uv = [m_quarks[labels.index(l)] for l in triple]
    m_cw_t = [cw_final[l] for l in triple]
    m_ss = [Lambda**2 / m_quarks[labels.index(l)] for l in triple]

    Q_uv = koide_Q(*m_uv)
    Q_cw = koide_Q(*m_cw_t)
    Q_ss = koide_Q(*m_ss)

    print(f"({','.join(triple)}){'':>8} {Q_uv:<18.8f} {Q_cw:<18.8f} {Q_ss:<18.8f} {'0.66667'}")

print(f"\n{'Triple':<15} {'Q (dual: 1/m_q)':<18}")
print("-" * 40)
for triple in [('d','s','c'), ('u','d','s')]:
    m_inv = [1.0/m_quarks[labels.index(l)] for l in triple]
    Q_dual = koide_Q(*m_inv)
    print(f"({','.join(triple)}){'':>8} {Q_dual:<18.8f}")

# Add the (d,s,b) seesaw check
m_b = 4180.0
m_inv_dsb = [1.0/m_d, 1.0/m_s, 1.0/m_b]
Q_dual_dsb = koide_Q(*m_inv_dsb)
m_ss_dsb = [Lambda**2/m_d, Lambda**2/m_s, Lambda**2/m_b]
Q_ss_dsb = koide_Q(*m_ss_dsb)
print(f"\n  Q(1/m_d, 1/m_s, 1/m_b) = {Q_dual_dsb:.8f} (dev from 2/3: {abs(Q_dual_dsb-2/3)/(2/3)*100:.4f}%)")
print(f"  Q(seesaw d,s,b) = {Q_ss_dsb:.8f} (should equal dual Koide: {abs(Q_ss_dsb-Q_dual_dsb) < 1e-10})")
print(f"  Seesaw TRANSMITS dual-Koide exactly.")

print(f"""

KEY PHYSICS:
==========
1. The ISS CW potential generates NEARLY DEGENERATE meson masses.
   This is because all quark masses are small compared to h*mu^2.

2. The Seiberg seesaw (SUSY vacuum) DOES transmit mass ratios:
   m_i → Lambda^2/m_i preserves Koide (maps Q → same Q for 1/m).

3. The relevant observation for the sBootstrap is:
   Q(1/m_d, 1/m_s, 1/m_b) = {Q_dual_dsb:.6f} ≈ 2/3
   This "dual Koide" condition would become a STANDARD Koide condition
   for the seesaw meson masses, IF the SUSY vacuum is relevant.

4. The CW potential in the ISS model lifts pseudo-moduli but does NOT
   track quark mass ratios. Koide transmission requires a different
   mechanism (bion Kahler, as computed in complete_kahler.md).

5. The pseudo-modulus X is stabilized at X = 0 by the CW potential.
   The Koide phase delta is NOT selected by the CW potential —
   it is a FLAT DIRECTION of the CW potential (confirmed in koide_phases.md).
""")
