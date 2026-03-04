"""
Mixing matrix from overlapping Koide mass triples.
Pure computation — no theory labels.
"""

import numpy as np
from scipy.optimize import minimize_scalar
import itertools

# ── masses (MeV) ──────────────────────────────────────────────────────────────
m_u  =    2.16
m_d  =    4.67
m_s  =   93.4
m_c  = 1270.0
m_b  = 4180.0
m_t  = 172760.0

# ── Koide parametrisation ────────────────────────────────────────────────────
#
# The parametrisation:
#   sv_k = √M₀ · (1 + √2 cos(2πk/3 + δ)),  k = 0,1,2
#
# Key algebraic identities:
#   Σ_k (1 + √2 cos(2πk/3+δ)) = 3        → Σ sv_k = 3√M₀
#   Σ_k (1 + √2 cos(2πk/3+δ))² = 6       → Σ sv_k² = 6M₀
#
# Therefore:   Q_signed ≡ (Σ sv_k²) / (Σ sv_k)² = 6M₀ / 9M₀ = 2/3
#
# Note: this Q uses the SIGNED sum (Σ sv_k), not Σ|sv_k|.
# The physical Koide Q is usually stated with unsigned √m_k, but for a triple
# with a negative sign entry the invariant that equals 2/3 is Q_signed.
# When all entries are positive, Q_signed = Q_unsigned = 2/3.
#
# M₀ is fixed by Σ sv = 3√M₀, so √M₀ = (Σ target_sv)/3.
# δ is then the only free parameter.

def koide_sv(M0, delta, perm=(0, 1, 2)):
    """sv[i] = √M₀ · (1 + √2 cos(2π·perm[i]/3 + δ))"""
    k = np.array(perm, dtype=float)
    return np.sqrt(M0) * (1.0 + np.sqrt(2) * np.cos(2*np.pi*k/3 + delta))


def koide_Q_signed(sv):
    """Q_signed = Σsv² / (Σsv)²  — equals 2/3 for ANY M₀,δ by identity."""
    sv = np.array(sv, dtype=float)
    return sv.dot(sv) / sv.sum()**2


def koide_Q_unsigned(sv):
    """Q_unsigned = Σsv² / (Σ|sv|)² — equals 2/3 only when all sv_k ≥ 0."""
    sv = np.array(sv, dtype=float)
    return sv.dot(sv) / np.abs(sv).sum()**2


def extract_koide(target_sv, try_permutations=True):
    """
    Fit M₀ and δ (and optionally the best k-permutation) to target_sv.

    Returns (M0, delta, perm, sv_fit, max_residual).
    """
    target = np.array(target_sv, dtype=float)
    # M₀ from sum constraint — exact regardless of δ and perm:
    sqrt_M0 = target.sum() / 3.0
    if sqrt_M0 <= 0:
        raise ValueError("Sum of target sv must be positive (√M₀ > 0).")
    M0 = sqrt_M0 ** 2

    perms = list(itertools.permutations([0, 1, 2])) if try_permutations else [(0, 1, 2)]

    best = None
    for perm in perms:
        def residual(delta, p=perm):
            return np.sum((koide_sv(M0, delta, p) - target)**2)

        # Coarse grid to find global minimum
        ds = np.linspace(0, 2*np.pi, 3600, endpoint=False)
        rv = np.array([residual(d) for d in ds])
        d0 = ds[np.argmin(rv)]

        res = minimize_scalar(residual, bounds=(d0-0.05, d0+0.05),
                              method='bounded', options={'xatol': 1e-14})
        delta_opt = res.x
        sv_fit = koide_sv(M0, delta_opt, perm)
        max_resid = np.max(np.abs(sv_fit - target))

        if best is None or max_resid < best[4]:
            best = (M0, delta_opt, perm, sv_fit, max_resid)

    return best   # (M0, delta, perm, sv_fit, max_residual)


# ══════════════════════════════════════════════════════════════════════════════
# Part 1: Koide parameters
# ══════════════════════════════════════════════════════════════════════════════
print("=" * 70)
print("PART 1: Koide parameters for the two triples")
print("=" * 70)

# Triple A: signed roots (-√m_s, √m_c, √m_b)
target_A = np.array([-np.sqrt(m_s), np.sqrt(m_c), np.sqrt(m_b)])
M0_A, delta_A, perm_A, sv_A, mr_A = extract_koide(target_A)

# Triple B: signed roots (√m_c, √m_b, √m_t)
target_B = np.array([np.sqrt(m_c), np.sqrt(m_b), np.sqrt(m_t)])
M0_B, delta_B, perm_B, sv_B, mr_B = extract_koide(target_B)

for lbl, target, M0, delta, perm, sv, mr in [
    ("A: (-√m_s, √m_c, √m_b)", target_A, M0_A, delta_A, perm_A, sv_A, mr_A),
    ("B: (√m_c, √m_b, √m_t)",  target_B, M0_B, delta_B, perm_B, sv_B, mr_B),
]:
    resid = target - sv
    Q_s  = koide_Q_signed(sv)
    Q_u  = koide_Q_unsigned(sv)
    print(f"\nTriple {lbl}")
    print(f"  k-perm (physical → Koide slot): {perm}")
    print(f"  Input sv:         {target}")
    print(f"  M₀  = {M0:.6f} MeV    √M₀ = {np.sqrt(M0):.6f} MeV^(1/2)")
    print(f"  δ   = {delta:.8f} rad  =  {np.degrees(delta):.4f}°  =  {delta/(2*np.pi/3):.6f}×(2π/3)")
    print(f"  Fitted sv:        {sv}")
    print(f"  Residuals:        {resid}")
    print(f"  Max |residual|:   {mr:.3e}  MeV^(1/2)")
    print(f"  Q_signed  (Σsv²/(Σsv)²):    {Q_s:.10f}  [should be 2/3 = {2/3:.10f}]")
    print(f"  Q_unsigned (Σsv²/(Σ|sv|)²): {Q_u:.10f}  [2/3 only if all sv>0]")

print(f"\nδ_B − δ_A = {delta_B - delta_A:.8f} rad = {np.degrees(delta_B - delta_A):.4f}°")

# Sanity: Q of PDG masses (unsigned)
Q_pdg_A = np.array([m_s, m_c, m_b])
Q_pdg_B = np.array([m_c, m_b, m_t])
print(f"\nDiagnostic — Q_unsigned of raw PDG mass triples:")
print(f"  Q(m_s, m_c, m_b) = {np.sum(Q_pdg_A)/(np.sum(np.sqrt(Q_pdg_A)))**2:.6f}")
print(f"  Q(m_c, m_b, m_t) = {np.sum(Q_pdg_B)/(np.sum(np.sqrt(Q_pdg_B)))**2:.6f}")
print(f"  Neither triple satisfies Q=2/3 exactly with PDG masses.")
print(f"  The fit residuals reflect how far the PDG masses are from the Koide curve.")

# ══════════════════════════════════════════════════════════════════════════════
# Part 2: Rotation matrix
# ══════════════════════════════════════════════════════════════════════════════
print("\n" + "=" * 70)
print("PART 2: Rotation matrix between the Koide eigenbases")
print("=" * 70)

# Orthonormal basis for 3D flavor space
e0 = np.array([1, 1, 1], dtype=float) / np.sqrt(3)   # democratic
f1 = np.array([1, -1, 0], dtype=float) / np.sqrt(2)  # λ₃
f2 = np.array([1, 1, -2], dtype=float) / np.sqrt(6)  # λ₈

# Koide 2D direction in Cartan plane for a given (delta, perm)
def cartan_direction(delta, perm=(0, 1, 2)):
    """Unit 2D vector in Cartan plane (f₁, f₂) for Koide phase δ."""
    # Cartan component of sv (per unit √M₀): √2 cos(2π·perm[j]/3 + δ), j=0,1,2
    v = np.array([np.sqrt(2)*np.cos(2*np.pi*perm[j]/3 + delta) for j in range(3)])
    x = np.dot(v, f1)
    y = np.dot(v, f2)
    r = np.hypot(x, y)
    return np.array([x, y]) / r if r > 1e-15 else np.array([1.0, 0.0])

nA_2D = cartan_direction(delta_A, perm_A)
nB_2D = cartan_direction(delta_B, perm_B)
phi_A = np.arctan2(nA_2D[1], nA_2D[0])
phi_B = np.arctan2(nB_2D[1], nB_2D[0])
dphi  = (phi_B - phi_A + np.pi) % (2*np.pi) - np.pi   # in (-π, π]

print(f"\nCartan-plane directions:")
print(f"  Triple A: φ_A = {np.degrees(phi_A):.4f}°,  n̂_A = ({nA_2D[0]:.6f}, {nA_2D[1]:.6f})")
print(f"  Triple B: φ_B = {np.degrees(phi_B):.4f}°,  n̂_B = ({nB_2D[0]:.6f}, {nB_2D[1]:.6f})")
print(f"  Δφ = {np.degrees(dphi):.4f}°")

R2D = np.array([[np.cos(dphi), -np.sin(dphi)],
                [np.sin(dphi),  np.cos(dphi)]])
print(f"\n2D rotation R_2D (Cartan plane, Δφ = {np.degrees(dphi):.4f}°):")
print(f"  [[{R2D[0,0]:+.6f}  {R2D[0,1]:+.6f}]")
print(f"   [{R2D[1,0]:+.6f}  {R2D[1,1]:+.6f}]]")

# Full 3D: identity on e0, R_2D on (f1,f2)
B = np.column_stack([e0, f1, f2])
Rblock = np.eye(3)
Rblock[1:, 1:] = R2D
R_3D = B @ Rblock @ B.T

print(f"\nFull 3D rotation R_3D:")
print(np.array2string(R_3D, precision=6, suppress_small=True))
print(f"det(R_3D) = {np.linalg.det(R_3D):.8f}")
print("R_3D^T R_3D:")
print(np.array2string(R_3D.T @ R_3D, precision=6, suppress_small=True))

# ══════════════════════════════════════════════════════════════════════════════
# Part 3: CKM-like angles from R_3D
# ══════════════════════════════════════════════════════════════════════════════
print("\n" + "=" * 70)
print("PART 3: CKM-like angles from R_3D vs PDG")
print("=" * 70)

def ckm_angles(V):
    s13 = float(np.clip(np.abs(V[0, 2]), 0, 1))
    t13 = np.arcsin(s13);  c13 = np.cos(t13)
    if c13 < 1e-10:
        return np.arctan2(abs(V[1,0]), abs(V[1,1])), 0.0, t13, 0.0
    t23 = np.arcsin(float(np.clip(abs(V[1,2])/c13, 0, 1)))
    t12 = np.arcsin(float(np.clip(abs(V[0,1])/c13, 0, 1)))
    dCP = float(-np.angle(V[0,2])) if abs(V[0,2]) > 1e-10 else 0.0
    return t12, t23, t13, dCP

t12, t23, t13, dCP = ckm_angles(R_3D)
s12, c12 = np.sin(t12), np.cos(t12)
s23, c23 = np.sin(t23), np.cos(t23)
s13v, c13v = np.sin(t13), np.cos(t13)

Vus_R = s12*c13v;  Vcb_R = s23*c13v;  Vub_R = s13v

Vus_pdg = 0.2243;  Vcb_pdg = 0.0422;  Vub_pdg = 0.00394

print(f"\nFrom R_3D (PDG parametrisation):")
print(f"  θ₁₂ = {np.degrees(t12):.4f}°  sin = {s12:.6f}")
print(f"  θ₂₃ = {np.degrees(t23):.4f}°  sin = {s23:.6f}")
print(f"  θ₁₃ = {np.degrees(t13):.4f}°  sin = {s13v:.6f}")
print(f"  δ_CP = {np.degrees(dCP):.4f}°")
print(f"\n  |V_us| from R_3D = {Vus_R:.6f}  (PDG {Vus_pdg}, ratio {Vus_R/Vus_pdg:.4f})")
print(f"  |V_cb| from R_3D = {Vcb_R:.6f}  (PDG {Vcb_pdg}, ratio {Vcb_R/Vcb_pdg:.4f})")
print(f"  |V_ub| from R_3D = {Vub_R:.6f}  (PDG {Vub_pdg}, ratio {Vub_R/Vub_pdg:.4f})")

# ══════════════════════════════════════════════════════════════════════════════
# Part 4: Koide direction angles
# ══════════════════════════════════════════════════════════════════════════════
print("\n" + "=" * 70)
print("PART 4: Koide direction angles in 3D flavor space and Cartan plane")
print("=" * 70)

def koide_3D_direction(delta, perm=(0, 1, 2)):
    v = np.array([1.0 + np.sqrt(2)*np.cos(2*np.pi*perm[j]/3 + delta) for j in range(3)])
    return v / np.linalg.norm(v)

n3D_A = koide_3D_direction(delta_A, perm_A)
n3D_B = koide_3D_direction(delta_B, perm_B)
angle_3D = np.degrees(np.arccos(np.clip(np.dot(n3D_A, n3D_B), -1, 1)))

nAf1, nAf2 = np.dot(n3D_A, f1), np.dot(n3D_A, f2)
nBf1, nBf2 = np.dot(n3D_B, f1), np.dot(n3D_B, f2)
phi_A_c = np.degrees(np.arctan2(nAf2, nAf1))
phi_B_c = np.degrees(np.arctan2(nBf2, nBf1))
dphi_c  = phi_B_c - phi_A_c

print(f"\n3D direction vectors:")
print(f"  n̂_A = {n3D_A}")
print(f"  n̂_B = {n3D_B}")
print(f"\n  3D angle(n̂_A, n̂_B) = {angle_3D:.4f}°  (Cabibbo = 13.00°, ratio {angle_3D/13.:.4f})")
print(f"\nCartan-plane projections:")
print(f"  φ_A = {phi_A_c:.4f}°,  φ_B = {phi_B_c:.4f}°,  |Δφ| = {abs(dphi_c):.4f}°")
print(f"  (Cabibbo = 13.00°, ratio {abs(dphi_c)/13.:.4f})")

# ══════════════════════════════════════════════════════════════════════════════
# Part 5: Oakes relation
# ══════════════════════════════════════════════════════════════════════════════
print("\n" + "=" * 70)
print("PART 5: Oakes relation")
print("=" * 70)

oakes_us = np.sqrt(m_d / m_s)
oakes_cb = np.sqrt(m_s / m_b)
oakes_ub = np.sqrt(m_d / m_b)

for lbl, val, pdg in [
    ("|V_us| ~ √(m_d/m_s)", oakes_us, Vus_pdg),
    ("|V_cb| ~ √(m_s/m_b)", oakes_cb, Vcb_pdg),
    ("|V_ub| ~ √(m_d/m_b)", oakes_ub, Vub_pdg),
]:
    print(f"\n  {lbl}: {val:.6f}  PDG: {pdg:.6f}  ratio: {val/pdg:.4f}  dev: {(val-pdg)/pdg*100:+.1f}%")

# ══════════════════════════════════════════════════════════════════════════════
# Part 6: Wolfenstein
# ══════════════════════════════════════════════════════════════════════════════
print("\n" + "=" * 70)
print("PART 6: Wolfenstein-like expansion")
print("=" * 70)

lam = np.sqrt(m_d / m_s)
A_mass = oakes_cb / lam**2
A_pdg  = Vcb_pdg  / Vus_pdg**2
rho_mass = oakes_ub / (A_mass * lam**3)
rho_pdg  = Vub_pdg  / (A_pdg  * Vus_pdg**3)

print(f"\nλ = √(m_d/m_s) = {lam:.6f}  (PDG {Vus_pdg}, ratio {lam/Vus_pdg:.4f})")
print(f"λ² = {lam**2:.6f},  λ³ = {lam**3:.6f}")
print(f"\n√(m_s/m_b) / λ² = {oakes_cb:.6f} / {lam**2:.6f} = {A_mass:.4f}  [Wolfenstein A from masses]")
print(f"√(m_d/m_b) / λ³ = {oakes_ub:.6f} / {lam**3:.6f} = {oakes_ub/lam**3:.4f}  [= A·|ρ−iη| from masses]")
print(f"(note: √(m_d/m_b)/λ³ = √(m_s/m_b)/λ² algebraically, so A·|ρ−iη| = A → |ρ−iη|=1 from masses)")
print(f"\nWolfenstein A:     masses={A_mass:.4f},  PDG={A_pdg:.4f},  ratio={A_mass/A_pdg:.4f}")
print(f"Wolfenstein |ρ−iη|: masses={rho_mass:.4f}  PDG={rho_pdg:.4f}   ratio={rho_mass/rho_pdg:.4f}")

# ══════════════════════════════════════════════════════════════════════════════
# Summary
# ══════════════════════════════════════════════════════════════════════════════
print("\n" + "=" * 70)
print("SUMMARY TABLE")
print("=" * 70)
print(f"\n{'Quantity':<40} {'Computed':<13} {'PDG/ref':<12} {'Ratio'}")
print("-" * 72)
for lbl, val, ref in [
    ("|V_us| Oakes √(m_d/m_s)",    oakes_us,  Vus_pdg),
    ("|V_cb| Oakes √(m_s/m_b)",    oakes_cb,  Vcb_pdg),
    ("|V_ub| Oakes √(m_d/m_b)",    oakes_ub,  Vub_pdg),
    ("Wolfenstein A (masses)",      A_mass,    A_pdg),
    ("Wolfenstein |ρ−iη| (masses)", rho_mass,  rho_pdg),
    ("|V_us| from R_3D",            Vus_R,     Vus_pdg),
    ("|V_cb| from R_3D",            Vcb_R,     Vcb_pdg),
    ("|V_ub| from R_3D",            Vub_R,     Vub_pdg),
    ("3D Koide angle (deg)",        angle_3D,  13.00),
    ("Cartan |Δφ| (deg)",           abs(dphi_c), 13.00),
]:
    print(f"  {lbl:<38} {val:<13.5g} {ref:<12.5g} {val/ref:.4f}")

# ══════════════════════════════════════════════════════════════════════════════
# Markdown report
# ══════════════════════════════════════════════════════════════════════════════
lines = []
def md(s=""): lines.append(s)

md("# Mixing Matrix from Overlapping Koide Mass Triples")
md()
md("**Computation date:** 2026-03-04")
md()
md("## Setup and Notation")
md()
md("The Koide parametrisation assigns a signed square-root value to each slot k:")
md()
md("    sv_k = √M₀ · (1 + √2 cos(2πk/3 + δ)),  k = 0, 1, 2")
md()
md("Algebraic identities (exact for any M₀ > 0, any δ, any permutation of k):")
md()
md("    Σ_k sv_k = 3√M₀         →  M₀ = (Σ sv_k)² / 9")
md("    Σ_k sv_k² = 6M₀")
md("    Q_signed ≡ Σsv² / (Σsv)² = 6M₀/9M₀ = 2/3")
md()
md("The invariant that equals 2/3 is **Q_signed** (using the signed sum Σsv).")
md("For triples with a negative entry (like (−√m_s, √m_c, √m_b)),")
md("Q_signed = 2/3 is the appropriate statement; Q_unsigned uses Σ|sv| and")
md("equals 2/3 only when all entries are non-negative.")
md()
md("Each physical triple is described by two free parameters (M₀, δ),")
md("with M₀ determined analytically from the signed sum. The best-fit δ")
md("and the permutation mapping physical indices to Koide slots k are")
md("found by minimising the fit residual over all 6 permutations.")
md()
md("Mass values (MeV): m_u = 2.16, m_d = 4.67, m_s = 93.4,")
md("m_c = 1270, m_b = 4180, m_t = 172760.")
md()
md("## Part 1: Koide Parameters")
md()

for lbl, M0, delta, perm, target, sv, mr in [
    ("Triple A: (−√m_s, √m_c, √m_b)",
     M0_A, delta_A, perm_A, target_A, sv_A, mr_A),
    ("Triple B: (√m_c, √m_b, √m_t)",
     M0_B, delta_B, perm_B, target_B, sv_B, mr_B),
]:
    resid = target - sv
    Q_s = koide_Q_signed(sv)
    Q_u = koide_Q_unsigned(sv)
    md(f"### {lbl}")
    md()
    md(f"Best-fit k-permutation (physical index → Koide slot): {perm}")
    md()
    md("| Parameter | Value |")
    md("|-----------|-------|")
    md(f"| M₀ (MeV) | {M0:.6f} |")
    md(f"| √M₀ (MeV^1/2) | {np.sqrt(M0):.6f} |")
    md(f"| δ (rad) | {delta:.8f} |")
    md(f"| δ / (2π/3) | {delta/(2*np.pi/3):.6f} |")
    md(f"| δ (deg) | {np.degrees(delta):.4f} |")
    md(f"| Q_signed = Σsv²/(Σsv)² | {Q_s:.10f} |")
    md(f"| Q_unsigned = Σsv²/(Σ\\|sv\\|)² | {Q_u:.6f} |")
    md(f"| \\|Q_signed − 2/3\\| | {abs(Q_s - 2/3):.2e} |")
    md(f"| Max fit residual (MeV^1/2) | {mr:.3e} |")
    md()
    md(f"Fitted sv:  ({sv[0]:.5f}, {sv[1]:.5f}, {sv[2]:.5f})")
    md(f"Target sv:  ({target[0]:.5f}, {target[1]:.5f}, {target[2]:.5f})")
    md(f"Residuals:  ({resid[0]:.3e}, {resid[1]:.3e}, {resid[2]:.3e})")
    md()

md(f"Phase difference: δ_B − δ_A = {np.degrees(delta_B-delta_A):.4f}°")
md()
md("**Note on fit residuals.** The parametrisation has two free parameters (M₀, δ)")
md("for three mass values, so an exact fit is not generically possible.")
md("The residuals (~0.5–1 MeV^1/2) reflect how far the PDG masses lie from")
md("the nearest Koide curve. Q_signed = 2/3 exactly by algebraic identity,")
md("regardless of fit quality.")
md()
md("## Part 2: Rotation Matrix R")
md()
md("The 3D flavor space decomposes into orthogonal subspaces:")
md()
md("- **Democratic direction** e₀ = (1,1,1)/√3: carries the Koide scale M₀")
md("  (invariant under δ-variation).")
md("- **Cartan plane**: spanned by f₁ = (1,−1,0)/√2 (λ₃) and")
md("  f₂ = (1,1,−2)/√6 (λ₈).")
md()
md("Within the Cartan plane the Koide direction vector rotates as δ changes.")
md("The rotation from Triple A to Triple B is determined by:")
md()
md(f"    Δφ = φ_B − φ_A = {np.degrees(dphi):.4f}°   (mod 360°, in (−180°, +180°])")
md()
md("2D rotation matrix in Cartan plane (angle Δφ = "
   + f"{np.degrees(dphi):.4f}°):")
md()
md("```")
md(f"R_2D = [[{R2D[0,0]:+.6f}  {R2D[0,1]:+.6f}]")
md(f"        [{R2D[1,0]:+.6f}  {R2D[1,1]:+.6f}]]")
md("```")
md()
md("Full 3D rotation (identity on e₀, R_2D on Cartan plane):")
md()
md("```")
md("R_3D =")
for row in R_3D:
    md(f"  [{row[0]:+.6f}  {row[1]:+.6f}  {row[2]:+.6f}]")
md("```")
md()
md(f"det(R_3D) = {np.linalg.det(R_3D):.8f}  (proper rotation)")
md()
md("## Part 3: CKM-like Angles from R_3D")
md()
md("PDG parametrisation applied to R_3D:")
md()
md("| Angle | Value (deg) | sin |")
md("|-------|-------------|-----|")
md(f"| θ₁₂ | {np.degrees(t12):.4f} | {s12:.6f} |")
md(f"| θ₂₃ | {np.degrees(t23):.4f} | {s23:.6f} |")
md(f"| θ₁₃ | {np.degrees(t13):.4f} | {s13v:.6f} |")
md(f"| δ_CP | {np.degrees(dCP):.4f} | — |")
md()
md("| Element | From R_3D | PDG | Ratio |")
md("|---------|-----------|-----|-------|")
md(f"| \\|V_us\\| | {Vus_R:.6f} | {Vus_pdg} | {Vus_R/Vus_pdg:.4f} |")
md(f"| \\|V_cb\\| | {Vcb_R:.6f} | {Vcb_pdg} | {Vcb_R/Vcb_pdg:.4f} |")
md(f"| \\|V_ub\\| | {Vub_R:.6f} | {Vub_pdg} | {Vub_R/Vub_pdg:.4f} |")
md()
md("The Koide-basis rotation does not reproduce CKM magnitudes.")
md(f"The Cartan rotation angle is Δφ = {np.degrees(abs(dphi)):.1f}°,")
md(f"which is {np.degrees(abs(dphi))/13.:.1f}× the Cabibbo angle.")
md()
md("## Part 4: Koide Direction Angles")
md()
md("| Quantity | Value | Cabibbo reference | Ratio |")
md("|---------|-------|-------------------|-------|")
md(f"| 3D angle (n̂_A · n̂_B) | {angle_3D:.4f}° | 13.00° | {angle_3D/13.:.4f} |")
md(f"| Cartan \\|Δφ\\| | {abs(dphi_c):.4f}° | 13.00° | {abs(dphi_c)/13.:.4f} |")
md()
md(f"The angle between the two Koide direction vectors is {angle_3D:.2f}°,")
md(f"approximately {angle_3D/13.:.1f}× the Cabibbo angle.")
md()
md("## Part 5: Oakes Relation")
md()
md("| Element | Formula | Computed | PDG | Deviation |")
md("|---------|---------|----------|-----|-----------|")
for lbl, formula, val, pdg in [
    ("\\|V_us\\|", "√(m_d/m_s)", oakes_us, Vus_pdg),
    ("\\|V_cb\\|", "√(m_s/m_b)", oakes_cb, Vcb_pdg),
    ("\\|V_ub\\|", "√(m_d/m_b)", oakes_ub, Vub_pdg),
]:
    md(f"| {lbl} | {formula} | {val:.6f} | {pdg} | {(val-pdg)/pdg*100:+.1f}% |")

md()
md(f"- √(m_d/m_s) = {oakes_us:.5f}: reproduces |V_us| = {Vus_pdg} at −0.3%. **Works well.**")
md(f"- √(m_s/m_b) = {oakes_cb:.5f}: vs |V_cb| = {Vcb_pdg}, off by factor {oakes_cb/Vcb_pdg:.1f}. **Fails.**")
md(f"- √(m_d/m_b) = {oakes_ub:.5f}: vs |V_ub| = {Vub_pdg}, off by factor {oakes_ub/Vub_pdg:.1f}. **Fails.**")
md()
md("## Part 6: Wolfenstein-like Expansion")
md()
md(f"Cabibbo parameter: λ = √(m_d/m_s) = {lam:.6f}  (PDG |V_us| = {Vus_pdg}, ratio {lam/Vus_pdg:.4f})")
md(f"λ² = {lam**2:.6f},  λ³ = {lam**3:.6f}")
md()
md("| Check | Computed | λⁿ | Ratio |")
md("|-------|----------|-----|-------|")
md(f"| √(m_s/m_b) vs λ² | {oakes_cb:.6f} | {lam**2:.6f} | {A_mass:.4f} |")
md(f"| √(m_d/m_b) vs λ³ | {oakes_ub:.6f} | {lam**3:.6f} | {oakes_ub/lam**3:.4f} |")
md()
md("Both ratios are equal (algebraic identity: √(m_d/m_b)/λ³ = √(m_s/m_b)/λ²).")
md("This common ratio is the Wolfenstein A parameter implied by masses.")
md()
md("| Wolfenstein parameter | From mass ratios | From PDG CKM | Ratio |")
md("|----------------------|-----------------|--------------|-------|")
md(f"| λ | {lam:.6f} | {Vus_pdg:.6f} | {lam/Vus_pdg:.4f} |")
md(f"| A | {A_mass:.4f} | {A_pdg:.4f} | {A_mass/A_pdg:.4f} |")
md(f"| \\|ρ−iη\\| | {rho_mass:.4f} | {rho_pdg:.4f} | {rho_mass/rho_pdg:.4f} |")
md()
md(f"The λ parameter matches PDG at 0.3%. The Wolfenstein A from masses ({A_mass:.3f})")
md(f"exceeds the PDG value ({A_pdg:.3f}) by a factor of {A_mass/A_pdg:.2f}.")
md(f"|ρ−iη| from masses ({rho_mass:.4f}) vs PDG ({rho_pdg:.4f}): ratio {rho_mass/rho_pdg:.4f}.")
md()
md("## Summary Table")
md()
md("| Quantity | Computed | Reference | Ratio |")
md("|---------|---------|-----------|-------|")
for lbl, val, ref, src in [
    ("|V_us| Oakes",               oakes_us,    Vus_pdg,   "PDG"),
    ("|V_cb| Oakes",               oakes_cb,    Vcb_pdg,   "PDG"),
    ("|V_ub| Oakes",               oakes_ub,    Vub_pdg,   "PDG"),
    ("Wolfenstein A",              A_mass,      A_pdg,     "PDG"),
    ("Wolfenstein |ρ−iη|",         rho_mass,    rho_pdg,   "PDG"),
    ("|V_us| from R_3D",           Vus_R,       Vus_pdg,   "PDG"),
    ("|V_cb| from R_3D",           Vcb_R,       Vcb_pdg,   "PDG"),
    ("|V_ub| from R_3D",           Vub_R,       Vub_pdg,   "PDG"),
    ("3D Koide angle (deg)",       angle_3D,    13.00,     "Cabibbo"),
    ("Cartan |Δφ| (deg)",          abs(dphi_c), 13.00,     "Cabibbo"),
]:
    md(f"| {lbl} | {val:.5g} | {ref:.5g} [{src}] | {val/ref:.4f} |")

md()
md("## Conclusions")
md()
md(f"1. **Koide parameters.** Both triples fitted successfully.")
md(f"   M₀_A = {M0_A:.3f} MeV, δ_A = {np.degrees(delta_A):.3f}°;")
md(f"   M₀_B = {M0_B:.3f} MeV, δ_B = {np.degrees(delta_B):.3f}°;")
md(f"   Δδ = {np.degrees(delta_B-delta_A):.3f}°.")
md(f"   Q_signed = 2/3 exactly by algebraic identity for both triples.")
md(f"   Fit residuals (~{mr_A:.2f}–{mr_B:.2f} MeV^1/2) measure how far PDG masses")
md(f"   deviate from any Koide curve.")
md()
md(f"2. **Oakes relation.** √(m_d/m_s) = {oakes_us:.4f} matches |V_us| = {Vus_pdg}")
md(f"   at −0.3%. The generalisation fails for |V_cb| and |V_ub|,")
md(f"   which are off by factors of {oakes_cb/Vcb_pdg:.1f} and {oakes_ub/Vub_pdg:.1f}.")
md()
md(f"3. **Wolfenstein power structure.** The hierarchy λ : λ² : λ³ is present")
md(f"   in the mass ratios with λ = √(m_d/m_s) ≈ {lam:.4f} (PDG: {Vus_pdg}).")
md(f"   The Wolfenstein A from masses ({A_mass:.3f}) is {A_mass/A_pdg:.2f}× the PDG value ({A_pdg:.3f}).")
md(f"   The mass-ratio ansatz captures the power-law hierarchy but not the scale.")
md()
md(f"4. **Koide-basis rotation R_3D.** Geometrically clean (det = +1),")
md(f"   rotation angle Δφ = {np.degrees(abs(dphi)):.1f}° ≈ {np.degrees(abs(dphi))/13.:.1f}× θ_C.")
md(f"   CKM-like elements from R_3D are {Vus_R/Vus_pdg:.1f}–{Vub_R/Vub_pdg:.0f}× too large.")
md()
md("5. **Overall.** The only precise mass-to-mixing connection found here")
md(f"   is the Oakes/Cabibbo relation |V_us| ≈ √(m_d/m_s) (−0.3%).")
md("   No construction based solely on Koide phase parameters reproduces")
md("   |V_cb| or |V_ub| at the correct magnitude.")

report = "\n".join(lines)
with open("/home/codexssh/phys3/results/koide_ckm.md", "w") as f:
    f.write(report)

print("\n\nMarkdown report written to /home/codexssh/phys3/results/koide_ckm.md")
