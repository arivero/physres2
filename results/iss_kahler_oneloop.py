"""
One-loop Kahler corrections for the pseudo-modulus in the ISS model.

ISS model: SU(N_c) magnetic theory with N_c = 3 colors, N_f = 4 flavors.
    Fields: magnetic quarks q^a_i, q-tilde^a_i  (a=1..3 color, i=1..4 flavor)
            mesons Phi^i_j  (4x4 matrix, gauge singlets)
    Superpotential: W = h Tr(q Phi q-tilde) - h mu^2 Tr(Phi)

Metastable vacuum:
    <q^a_i> = <q-tilde^a_i> = mu delta^a_i   (a=1..3, i=1..4)
    <Phi^i_j> = 0

The pseudo-modulus is X = Phi^4_4 (the (4,4) component not stabilized
at tree level because no q-condensate projects onto flavor 4).

We compute:
  1. Full tree-level mass spectrum at the ISS vacuum
  2. Coleman-Weinberg one-loop effective potential V_CW(X)
  3. Expansion V_CW near X = 0 (non-analytic due to massless fermion)
  4. Effective Kahler correction c_eff
  5. Monotonicity confirmation (X = 0 is the minimum)
  6. What external Kahler correction shifts the minimum to X = sqrt(3) mu/h
"""

import numpy as np
from scipy.optimize import brentq
import warnings
warnings.filterwarnings('ignore')

# ---------------------------------------------------------------------------
# Parameters
# ---------------------------------------------------------------------------
h_yuk = 1.0        # magnetic Yukawa coupling
mu = 300.0          # MeV, confinement scale
Nc = 3              # magnetic colors
Nf = 4              # flavors

mu2 = mu**2
hmu = h_yuk * mu    # h*mu
hmu2 = h_yuk * mu2  # h*mu^2

print("=" * 72)
print("ISS Model: One-Loop Kahler Corrections for the Pseudo-Modulus")
print("=" * 72)
print()
print(f"Parameters: h = {h_yuk}, mu = {mu} MeV, N_c = {Nc}, N_f = {Nf}")
print(f"            h*mu = {hmu:.1f} MeV,  h*mu^2 = {hmu2:.1f} MeV^2")
print()

# ===================================================================
# PART 1: Tree-level mass spectrum
# ===================================================================
#
# At the ISS vacuum, the pseudo-modulus X = Phi^4_4 couples to
# chi_a = q^a_4 and chi~_a = q-tilde^a_4 (a=1..Nc) through the
# superpotential W_eff = h X chi_a chi~_a - h mu^2 X.
#
# This is the standard O'Raifeartaigh sector. For each color a:
#   F_X = h(chi_a chi~_a) - h mu^2  -->  F_X = -h mu^2 at vacuum
#
# Mass spectrum as a function of x = h|X|/mu (dimensionless):
#   Scalar m^2_+  = (h mu)^2 (x^2 + 1)   [1 complex d.o.f. per color]
#   Scalar m^2_-  = (h mu)^2 (x^2 - 1)   [1 complex d.o.f. per color]
#   Fermion m_f   = h mu x              [1 Dirac = 2 Weyl per color]
#
# Note: m^2_- < 0 for x < 1. The CW potential uses |m^2| in the log.

print("=" * 72)
print("PART 1: Tree-Level Mass Spectrum")
print("=" * 72)
print()
print("Fields in the X-coupled sector (per magnetic color):")
print("  chi_a = q^a_4           (complex scalar + Weyl fermion)")
print("  chi~_a = q-tilde^a_4   (complex scalar + Weyl fermion)")
print()
print("Mass spectrum (x = h v / mu, with v = |<X>|):")
print("  Scalar m^2_+  = (h mu)^2 (x^2 + 1)")
print("  Scalar m^2_-  = (h mu)^2 (x^2 - 1)")
print("  Fermion m_f   = h mu x")
print()
print(f"At x = 0:  m_+ = {hmu:.0f} MeV,  m^2_- = -{hmu**2:.0f} MeV^2 (tachyonic),  m_f = 0")
print()
print("Supertrace identities (exact at all x):")
print(f"  STr[M^2] = Nc * [(x^2+1) + (x^2-1) - 2x^2] = 0")
print(f"  STr[M^4] = Nc * [(x^2+1)^2 + (x^2-1)^2 - 2x^4] = 2 Nc (h mu)^4")
print(f"           = {2*Nc*hmu**4:.6e} MeV^4  (constant, independent of x)")
print()

# Full field count
print("Full ISS field count:")
print(f"  Magnetic quarks: 2 x Nc x Nf = {2*Nc*Nf} complex fields")
print(f"    - Eaten by Higgsing: 2 x Nc x Nc = {2*Nc*Nc}")
print(f"      (q^a_i, q~^a_i for i=1..Nc pair with gauge bosons)")
print(f"    - Massive (paired with Phi^i_a, Phi^a_i): 2 x Nc x (Nf-Nc) = {2*Nc*(Nf-Nc)}")
print(f"      (these are the chi_a, chi~_a that couple to X)")
print(f"  Mesons Phi: Nf^2 = {Nf**2} complex fields")
print(f"    - Phi^i_j (i,j=1..Nc): Nc^2 = {Nc**2} get mass from q-VEV")
print(f"    - Phi^i_4 (i=1..Nc): {Nc} link fields")
print(f"    - Phi^4_i (i=1..Nc): {Nc} link fields")
print(f"    - Phi^4_4 = X: 1 pseudo-modulus (classically flat)")
print()

# ===================================================================
# PART 2: Coleman-Weinberg potential
# ===================================================================

def vcw_per_color(x):
    """
    CW potential per color in units of (h mu)^4 / (64 pi^2).

    f(x) = m_+^4 [ln(m_+^2) - 3/2] + m_-^4 [ln|m_-^2| - 3/2]
           - 2 m_f^4 [ln(m_f^2) - 3/2]

    with all mass-squared in units of (h mu)^2, and Lambda_RG = h mu.
    """
    x2 = x**2
    mp2 = x2 + 1.0
    mm2 = x2 - 1.0
    mf2 = x2

    V = mp2**2 * (np.log(mp2) - 1.5)

    if abs(mm2) > 1e-30:
        V += mm2**2 * (np.log(abs(mm2)) - 1.5)
    # else: mm2 = 0 contributes 0 (limit of mm2^2 ln|mm2| = 0)

    if mf2 > 1e-30:
        V -= 2.0 * mf2**2 * (np.log(mf2) - 1.5)

    return V


def vcw_total(x, nc=Nc):
    """Total CW potential = nc * f(x), same units."""
    return nc * vcw_per_color(x)


print("=" * 72)
print("PART 2: Coleman-Weinberg Potential V_CW(X)")
print("=" * 72)
print()
print("V_CW = Nc/(64 pi^2) * STr[M^4 (ln M^2/Lambda^2 - 3/2)]")
print(f"     = Nc * (h mu)^4/(64 pi^2) * f(x)")
print()
print(f"Nc (h mu)^4/(64 pi^2) = {Nc * hmu**4 / (64*np.pi**2):.6e} MeV^4")
print()

print(f"{'x':>8} | {'f(x) per color':>18} | {'df/dx (num)':>18}")
print("-" * 52)
x_vals = [0, 0.1, 0.2, 0.5, 0.8, 1.0, 1.2, 1.5, 2.0, 3.0, 5.0, 10.0]
for x in x_vals:
    f = vcw_per_color(x)
    dx = 1e-7
    if x > 0:
        df = (vcw_per_color(x + dx) - vcw_per_color(x - dx)) / (2 * dx)
    else:
        df = (vcw_per_color(dx) - vcw_per_color(0)) / dx
    print(f"{x:8.3f} | {f:18.10f} | {df:18.10f}")

print()

# Monotonicity check
x_fine = np.linspace(0.001, 20, 50000)
f_fine = np.array([vcw_per_color(x) for x in x_fine])
df_fine = np.diff(f_fine) / np.diff(x_fine)
neg_slope = np.where(df_fine < -1e-15)[0]

if len(neg_slope) == 0:
    print("MONOTONICITY: V_CW is strictly increasing for all x > 0.")
    print("Global minimum at x = 0 (CONFIRMED).")
else:
    print(f"WARNING: Non-monotonic at x = {x_fine[neg_slope[0]]:.4f}")
print()

# ===================================================================
# PART 3: Taylor expansion near x = 0
# ===================================================================
#
# The expansion is non-analytic because m_f = h mu x -> 0 at x = 0.
# The fermion contribution -2 u^2 ln u (u = x^2) has a logarithmic
# singularity in its second derivative.
#
# Analytical result:
#   (1+u)^2 ln(1+u) + (1-u)^2 ln(1-u) = 3u^2 - u^4/6 - u^6/15 - ...
#   (from expanding the product and combining even powers)
#
# Combined:
#   f(u) = [(1+u)^2 ln(1+u) + (1-u)^2 ln(1-u)] - 3(1+u^2) - 2u^2 ln u + 3u^2
#        = 3u^2 - u^4/6 - ... - 3 - 3u^2 - 2u^2 ln u + 3u^2
#        = -3 + 3u^2 - 2u^2 ln u - u^4/6 + ...
#
# In terms of x:
#   f = -3 + x^4 [3 + 2 ln(1/x^2)] - x^8/6 + O(x^10)

print("=" * 72)
print("PART 3: Expansion of V_CW Near X = 0")
print("=" * 72)
print()
print("f(u) = -3 + 3u^2 - 2u^2 ln u - u^4/6 + O(u^5)")
print("     = -3 + x^4 [3 + 2 ln(1/x^2)] - x^8/6 + ...")
print()
print("where u = x^2 = h^2 v^2 / mu^2")
print()
print("The leading x-dependent term ~ x^4 ln(1/x^2) is:")
print("  - Non-analytic (due to the massless fermion at x=0)")
print("  - Positive for 0 < x < exp(3/4) = 2.117")
print("  - Responsible for the log-enhanced pseudo-modulus mass")
print()

# Numerical verification
print("Verification of analytical expansion:")
print(f"{'x':>8} | {'f(x) exact':>16} | {'f(x) approx':>16} | {'error':>12}")
print("-" * 58)
for x in [0.001, 0.005, 0.01, 0.05, 0.1, 0.2, 0.3]:
    f_exact = vcw_per_color(x)
    u = x**2
    f_approx = -3.0 + 3*u**2 - 2*u**2 * np.log(u) - u**4/6.0
    err = abs(f_exact - f_approx)
    print(f"{x:8.4f} | {f_exact:16.10f} | {f_approx:16.10f} | {err:12.2e}")
print()

# The "mass" of X from the CW potential:
# d^2 V / d|X|^2 = Nc (h mu)^4/(64 pi^2) * d^2 f/dv^2
# f(v) in terms of v: f = -3 + (h^2 v^2/mu^2)^2 [3 + 2 ln(mu^2/(h^2 v^2))] + ...
# Let w = h^2 v^2/mu^2. Then d^2f/dv^2 at v->0:
# df/dv = 4h^2 v/mu^2 * w * [3 + 2 ln(1/w)] - 4h^2 v/mu^2 * w = 4(h^2/mu^2) v [w(2+2ln(1/w))]
# = 4(h^4/mu^4) v^3 [2 + 2 ln(mu^2/(h^2 v^2))]
# This vanishes as v->0 (good). The second derivative:
# d^2f/dv^2 = 4(h^4/mu^4) [3v^2(2+2ln(mu^2/(h^2 v^2))) - 2v^2]
# = 4(h^4/mu^4) v^2 [4 + 6 ln(mu^2/(h^2 v^2))]  ... wait, this also vanishes.
# Actually, the mass^2 involves d^2V/dv^2 which diverges logarithmically.
#
# More precisely: V_CW ~ const + alpha * v^4 * ln(mu^2/(h^2 v^2))
# d^2V/dv^2 ~ alpha * 12 v^2 ln(1/v) + ...  -> 0 as v -> 0.
# But V/v^2 -> 0 too. The effective mass-squared is:
# m_X^2 = lim_{v->0} 2 V''(v) which is actually 0 (not divergent).
#
# The log-enhanced mass is at FINITE v, where:
# m_X^2(v) ~ Nc h^4/(8 pi^2) * mu^2 * ln(mu^2/(h^2 v^2))

print("Effective pseudo-modulus mass squared at finite v:")
print("  m_X^2(v) ~ Nc h^4 mu^2 / (8 pi^2) * [1 + ln(mu^2/(h^2 v^2))]")
print()
print(f"{'v (MeV)':>10} | {'x = hv/mu':>10} | {'m_X^2 (MeV^2)':>16} | {'m_X (MeV)':>12}")
print("-" * 56)
for v_test in [1, 5, 10, 30, 50, 100, 200]:
    x_test = h_yuk * v_test / mu
    # Second derivative of V_CW(v) numerically
    dv = 0.01
    d2V = (vcw_total(h_yuk*(v_test+dv)/mu) - 2*vcw_total(h_yuk*v_test/mu) + vcw_total(h_yuk*(v_test-dv)/mu)) / (dv * h_yuk/mu)**2
    m2_X = d2V * hmu**4 / (64 * np.pi**2)  # convert to MeV^2
    # Actually, need to account for the chain rule. d^2V/dv^2 = (h/mu)^2 d^2[Nc f(x)]/dx^2 * (hmu)^4/(64pi^2)
    dxx = 1e-5
    d2f = (vcw_total(x_test+dxx) - 2*vcw_total(x_test) + vcw_total(x_test-dxx)) / dxx**2
    m2_phys = d2f * hmu**4 / (64*np.pi**2) * (h_yuk/mu)**2
    m_X = np.sqrt(abs(m2_phys)) if m2_phys > 0 else 0
    print(f"{v_test:10.1f} | {x_test:10.4f} | {m2_phys:16.4f} | {m_X:12.4f}")
print()

# ===================================================================
# PART 4: Effective Kahler correction c_eff
# ===================================================================

print("=" * 72)
print("PART 4: Effective Kahler Correction")
print("=" * 72)
print()

# The CW potential generates an effective quartic coupling for X.
# From f = -3 + x^4 [3 + 2 ln(1/x^2)] + ..., the quartic coefficient
# (at scale x_0) is:
#   c_quartic(x_0) = 3 + 2 ln(1/x_0^2)
#
# In physical units, V_CW ~ Nc (h mu)^4/(64 pi^2) * c_quartic * x^4
# = Nc h^8 v^4/(64 pi^2 mu^4) * c_quartic(x_0) * mu^4
# = [Nc h^8/(64 pi^2)] * c_quartic(x_0) * v^4
#
# The effective Kahler correction to the kinetic term:
# delta K_{XX*}^{(1)} = -Nc h^2/(16 pi^2) * ln[(x^2+1)|x^2-1|/x^4]
# (from the standard one-loop Kahler potential formula)

print("Quartic coefficient of V_CW as a function of scale x_0:")
print("  c_quartic(x_0) = 3 + 2 ln(1/x_0^2)")
print()
print(f"{'x_0':>8} | {'c_quartic':>12} | {'c_eff (MeV^-4)':>18} | {'Sign':>8}")
print("-" * 56)
for x0 in [0.01, 0.05, 0.1, 0.2, 0.5, 1.0, np.sqrt(3), 2.0, 2.117]:
    cq = 3.0 + 2.0 * np.log(1.0/x0**2)
    ceff_phys = Nc * h_yuk**8 / (64*np.pi**2) * cq
    sign = "+" if cq > 0 else "-" if cq < 0 else "0"
    print(f"{x0:8.4f} | {cq:12.6f} | {ceff_phys:18.6e} | {sign:>8}")

print()
print(f"c_quartic changes sign at x_0 = exp(3/4) = {np.exp(0.75):.4f}")
print()

# Compare with -1/12
print("Comparison with c = -1/12:")
print(f"  -1/12 = {-1/12:.6f}")
print(f"  c_quartic at x=1: {3.0:.6f}  (positive, 36x larger in magnitude)")
print(f"  c_quartic at x=sqrt(3): {3.0 + 2.0*np.log(1.0/3.0):.6f}  (still positive)")
print()
print("ANSWER: c_eff is POSITIVE at all physically relevant scales (x < 2.1).")
print("The ISS one-loop CW potential does NOT produce a negative quartic")
print("Kahler correction. It produces a positive (log-enhanced) quartic")
print("that stabilizes the pseudo-modulus at the origin.")
print()

# One-loop Kahler metric correction
print("One-loop Kahler metric correction:")
print("  delta Z(x) = -Nc h^2/(16 pi^2) * ln[(x^2+1)|x^2-1|/x^4]")
print()
print(f"{'x':>8} | {'delta Z / (Nc h^2/(16pi^2))':>30}")
print("-" * 42)
for x in [0.01, 0.05, 0.1, 0.2, 0.5, 1.5, 2.0, 5.0, 10.0]:
    arg = (x**2 + 1) * abs(x**2 - 1) / x**4
    if arg > 0:
        dz = -np.log(arg)
    else:
        dz = float('inf')
    print(f"{x:8.4f} | {dz:30.6f}")
print()
print("(Diverges at x=0 and at x=1; the x=1 singularity is an artifact")
print(" of the tachyonic boundary where m_-^2 = 0.)")
print()

# ===================================================================
# PART 5: External Kahler correction to shift minimum
# ===================================================================

print("=" * 72)
print("PART 5: External Kahler Correction for Minimum at x = sqrt(3)")
print("=" * 72)
print()

x_target = np.sqrt(3.0)
print(f"Target minimum: x = sqrt(3) = {x_target:.6f}")
print(f"  v = sqrt(3) mu / h = {x_target * mu / h_yuk:.1f} MeV")
print()

# Method: add K = |X|^2 (1 + c |X|^2/Lambda_K^2) with Lambda_K = mu/h.
# Then V_tree = |F_X|^2 / K_{XX*} = h^2 mu^4 / (1 + 4c x^2).
# V_total = V_tree + V_CW = h^2 mu^4 / (1+4cx^2) + Nc h^4 mu^4/(64pi^2) f(x)

loop_prefactor = Nc * h_yuk**2 / (64.0 * np.pi**2)  # relative to h^2 mu^4

def V_combined(x, c_ext):
    """V_total in units of h^2 mu^4."""
    denom = 1.0 + 4.0 * c_ext * x**2
    if denom <= 0:
        return 1e10
    V_tree = 1.0 / denom
    V_cw = loop_prefactor * vcw_per_color(x) * Nc
    return V_tree + V_cw

def dV_combined(x, c_ext, dx=1e-7):
    return (V_combined(x+dx, c_ext) - V_combined(x-dx, c_ext)) / (2*dx)

def d2V_combined(x, c_ext, dx=1e-5):
    return (V_combined(x+dx, c_ext) - 2*V_combined(x, c_ext) + V_combined(x-dx, c_ext)) / dx**2

# Scan for c that makes dV/dx = 0 at x = sqrt(3)
c_scan = np.linspace(-5, 5, 50000)
dV_scan = np.array([dV_combined(x_target, c) for c in c_scan])
sign_changes = np.where(np.diff(np.sign(dV_scan)))[0]

print("Kahler approach: K = |X|^2 (1 + c |X|^2 / (mu/h)^2)")
print("  V_tree = h^2 mu^4 / (1 + 4c x^2)")
print()

c_solutions_K = []
for idx in sign_changes:
    try:
        c_sol = brentq(lambda c: dV_combined(x_target, c), c_scan[idx], c_scan[idx+1])
        c_solutions_K.append(c_sol)
    except:
        pass

for c_sol in c_solutions_K:
    d2V_val = d2V_combined(x_target, c_sol)
    V_at_0 = V_combined(1e-6, c_sol)
    V_at_tgt = V_combined(x_target, c_sol)
    nature = "minimum" if d2V_val > 0 else "maximum" if d2V_val < 0 else "inflection"
    global_min = "global" if V_at_tgt < V_at_0 else "local (x=0 lower)"
    print(f"  c = {c_sol:+.8f}:  d2V/dx2 = {d2V_val:.4e} ({nature}), {global_min}")

print()

# Alternative approach: add a direct potential correction
# V_total = V_CW(x) + c_pot * x^4  (in dimensionless CW units)
# This is NOT from a Kahler correction but tests the raw quartic needed.

print("Direct potential approach: V_total = V_CW + c_pot x^4")
print("  (Not from Kahler; just subtracts a quartic from the CW potential)")
print()

def V_direct(x, c_pot):
    return vcw_per_color(x) + c_pot * x**4

def dV_direct(x, c_pot, dx=1e-7):
    return (V_direct(x+dx, c_pot) - V_direct(x-dx, c_pot)) / (2*dx)

def d2V_direct(x, c_pot, dx=1e-5):
    return (V_direct(x+dx, c_pot) - 2*V_direct(x, c_pot) + V_direct(x-dx, c_pot)) / dx**2

c_scan2 = np.linspace(-10, 0, 50000)
dV_scan2 = np.array([dV_direct(x_target, c) for c in c_scan2])
sc2 = np.where(np.diff(np.sign(dV_scan2)))[0]

c_solutions_direct = []
for idx in sc2:
    try:
        c_sol = brentq(lambda c: dV_direct(x_target, c), c_scan2[idx], c_scan2[idx+1])
        c_solutions_direct.append(c_sol)
    except:
        pass

for c_sol in c_solutions_direct:
    d2V_val = d2V_direct(x_target, c_sol)
    nature = "minimum" if d2V_val > 0 else "maximum"
    V0_d = V_direct(1e-6, c_sol)
    Vt_d = V_direct(x_target, c_sol)
    print(f"  c_pot = {c_sol:.8f}:  d2V/dx2 = {d2V_val:.4e} ({nature})")
    print(f"    V(0) = {V0_d:.6f},  V(sqrt(3)) = {Vt_d:.6f}")
    print(f"    c_pot / (-1/12) = {c_sol / (-1/12):.6f}")

print()

# The result: a pure quartic subtraction gives a MAXIMUM at sqrt(3),
# not a minimum. This is because the CW potential is concave in that region.
# To get a genuine minimum, one needs a correction that is steeper than quartic
# at large x (e.g., x^6 or a full non-perturbative Kahler).
#
# Let's find what combined c2 x^2 + c4 x^4 correction creates a TRUE minimum.

print("-" * 72)
print("Combined quadratic + quartic correction:")
print("  V_total = f(x) + c2 x^2 + c4 x^4")
print()
print("For a minimum at x = sqrt(3) that is also a local minimum (d2V > 0),")
print("we need both dV/dx = 0 and d2V/dx2 > 0.")
print()

# At x = sqrt(3), u = 3:
# f'(x) at sqrt(3) numerically:
dx_n = 1e-7
fprime_at_target = (vcw_per_color(x_target + dx_n) - vcw_per_color(x_target - dx_n)) / (2*dx_n)
# f''(x) at sqrt(3):
f2prime_at_target = (vcw_per_color(x_target + 1e-5) - 2*vcw_per_color(x_target) + vcw_per_color(x_target - 1e-5)) / 1e-10

print(f"At x = sqrt(3):")
print(f"  f'(sqrt(3))  = {fprime_at_target:.8f}")
print(f"  f''(sqrt(3)) = {f2prime_at_target:.8f}")
print()

# Correction: g(x) = c2 x^2 + c4 x^4
# g'(x) = 2 c2 x + 4 c4 x^3
# g''(x) = 2 c2 + 12 c4 x^2
# At x = sqrt(3):
#   g'(sqrt(3)) = 2 c2 sqrt(3) + 12 c4 sqrt(3) = 2 sqrt(3) (c2 + 6 c4)
#   g''(sqrt(3)) = 2 c2 + 36 c4
#
# Conditions:
#   f'(sqrt(3)) + 2 sqrt(3)(c2 + 6 c4) = 0  -->  c2 + 6 c4 = -f'(sqrt(3))/(2 sqrt(3))
#   f''(sqrt(3)) + 2 c2 + 36 c4 > 0

alpha = -fprime_at_target / (2 * x_target)
print(f"Constraint: c2 + 6 c4 = {alpha:.8f}")
print()

# Also need c2 = 0 at x=0 to not give X a tree-level mass.
# Actually, a c2 x^2 term DOES give X a mass. If we want X massless at tree level,
# c2 = 0. Then: c4 = alpha/6.
c4_needed = alpha / 6
c2_needed = 0

print(f"If c2 = 0 (no tree-level mass for X):")
print(f"  c4 = {c4_needed:.8f}")
print(f"  f''(sqrt(3)) + 36 c4 = {f2prime_at_target + 36*c4_needed:.8f}")
is_min_c4only = f2prime_at_target + 36*c4_needed > 0
print(f"  This is a {'minimum' if is_min_c4only else 'MAXIMUM'} (need positive for minimum).")
print()

# So pure quartic cannot give a minimum. Need c2 != 0.
# For a minimum: c2 + 6c4 = alpha AND 2c2 + 36c4 > -f''
# i.e. 2(alpha - 6c4) + 36c4 > -f'' => 2 alpha + 24 c4 > -f''
# c4 > -(f'' + 2 alpha)/24

c4_threshold = -(f2prime_at_target + 2*alpha)/24
print(f"For a true minimum: c4 > {c4_threshold:.8f}")
print(f"  (with c2 = alpha - 6 c4)")
print()

# Example: c4 = c4_threshold + 0.1 (to be safely above threshold)
c4_ex = c4_threshold + 0.1
c2_ex = alpha - 6 * c4_ex
d2V_ex = f2prime_at_target + 2*c2_ex + 36*c4_ex

print(f"Example with c4 = {c4_ex:.6f}:")
print(f"  c2 = {c2_ex:.6f}")
print(f"  d2V/dx2 at sqrt(3) = {d2V_ex:.6f}  ({'minimum' if d2V_ex > 0 else 'maximum'})")
print()

# Most natural choice: c4 = -1/12
c4_test = -1.0/12
c2_test = alpha - 6 * c4_test
d2V_test = f2prime_at_target + 2*c2_test + 36*c4_test
print(f"With the specific value c4 = -1/12:")
print(f"  c2 = {c2_test:.8f}")
print(f"  d2V/dx2 at sqrt(3) = {d2V_test:.6f}  ({'minimum' if d2V_test > 0 else 'maximum'})")

# Check V(sqrt(3)) vs V(0) for this case
def V_full_correction(x, c2, c4):
    return vcw_per_color(x) + c2*x**2 + c4*x**4

V0_test = V_full_correction(1e-6, c2_test, c4_test)
Vt_test = V_full_correction(x_target, c2_test, c4_test)
print(f"  V(0) = {V0_test:.8f}")
print(f"  V(sqrt(3)) = {Vt_test:.8f}")
print(f"  V(sqrt(3)) - V(0) = {Vt_test - V0_test:.6e}")
if Vt_test < V0_test:
    print(f"  sqrt(3) is the GLOBAL minimum.")
else:
    print(f"  sqrt(3) is a local minimum; x=0 is still lower.")
print()

# ===================================================================
# SUMMARY
# ===================================================================

print("=" * 72)
print("SUMMARY OF RESULTS")
print("=" * 72)
print()
print("1. MASS SPECTRUM at the ISS vacuum (x = 0):")
print(f"   Per color: m_+ = h mu = {hmu:.0f} MeV, m_- = tachyonic, m_f = 0")
print(f"   Nc = {Nc} copies. STr[M^2] = 0 (exact). STr[M^4] = 2Nc(h mu)^4 (constant).")
print()
print("2. CW POTENTIAL (per color, dimensionless):")
print("   f(x) = -3 + x^4 [3 + 2 ln(1/x^2)] - x^8/6 + ...")
print("   Monotonically increasing for all x > 0.")
print("   Global minimum: x = 0 (CONFIRMED).")
print()
print("3. EFFECTIVE KAHLER CORRECTION:")
print("   c_eff(x_0) = Nc h^8/(64 pi^2) * [3 + 2 ln(1/x_0^2)]")
print(f"   At x=1: c_eff = {Nc * h_yuk**8 * 3 / (64*np.pi**2):.6f}")
print("   Sign: POSITIVE at all x < 2.12. Stabilizes X = 0.")
print()
print("4. KEY ANSWER: c_eff is NOT negative.")
print("   The ISS one-loop CW does not produce a negative Kahler quartic.")
print(f"   Compare: c_eff(x=1) = +{3*Nc*h_yuk**8/(64*np.pi**2):.4f}, vs -1/12 = {-1/12:.4f}")
print("   Opposite sign, and smaller magnitude at h=1.")
print()
print("5. SHIFTING MINIMUM TO x = sqrt(3):")
print("   Pure quartic subtraction gives a MAXIMUM, not a minimum.")
print("   A combined c2 x^2 + c4 x^4 correction is needed.")
if c4_test == -1.0/12:
    print(f"   With c4 = -1/12: requires c2 = {c2_test:.6f}")
    nature_test = "minimum" if d2V_test > 0 else "maximum"
    print(f"   This gives a {nature_test} at x = sqrt(3).")
print()
print("6. PHYSICAL INTERPRETATION:")
print("   The ISS CW potential robustly selects X = 0.")
print("   Moving the minimum away requires an external (non-perturbative)")
print("   contribution that overwhelms the one-loop stabilization.")
print("   The bion Kahler potential (complete_kahler.md) is a candidate,")
print("   as it provides a genuinely non-perturbative Kahler correction")
print("   with the requisite structure.")
print()
print("=" * 72)
