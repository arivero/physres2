#!/usr/bin/env python3
"""
Tau meson partner prediction.

The non-universal B-term splitting at t = sqrt(3) gives:
  muon-pion:  m_pi^2 - m_mu^2 = f_pi^2   (founding relation, 2% accuracy)
  tau-???:    m_???^2 - m_tau^2 = (2+sqrt3)^2 * f_pi^2

This predicts the tau's meson partner at 1810 MeV.
Which meson is it? What corrections shift this number?
"""
import numpy as np

print("=" * 70)
print("TAU MESON PARTNER PREDICTION")
print("=" * 70)

# Physical constants
m_e = 0.511       # MeV
m_mu = 105.658    # MeV
m_tau = 1776.86   # MeV
m_pi = 139.570    # MeV (pi+-)
f_pi = 92.2       # MeV

# The founding relation
founding = m_pi**2 - m_mu**2
print(f"\n--- Founding relation ---")
print(f"  m_pi^2 - m_mu^2 = {founding:.0f} MeV^2")
print(f"  f_pi^2           = {f_pi**2:.0f} MeV^2")
print(f"  Ratio: {founding/f_pi**2:.4f} (should be 1)")
print(f"  Deviation: {(founding/f_pi**2 - 1)*100:.2f}%")

# The B-term splitting ratio
ratio = (2 + np.sqrt(3))**2
print(f"\n--- Non-universal B-term splitting ---")
print(f"  cos^2(pi/12) / sin^2(pi/12) = (2+sqrt3)^2 = {ratio:.4f}")
print(f"  Tau splitting = {ratio:.2f} * (m_pi^2 - m_mu^2)")
print(f"               = {ratio * founding:.0f} MeV^2")

# Prediction
tau_split = ratio * founding
m_partner_sq = m_tau**2 + tau_split
m_partner = np.sqrt(m_partner_sq)
print(f"\n  m(tau partner)^2 = m_tau^2 + {tau_split:.0f}")
print(f"                   = {m_tau**2:.0f} + {tau_split:.0f}")
print(f"                   = {m_partner_sq:.0f}")
print(f"  m(tau partner)   = {m_partner:.1f} MeV")

# =====================================================================
# 1. Census of mesons near 1810 MeV
# =====================================================================
print("\n" + "=" * 70)
print("1. MESON CENSUS NEAR 1810 MeV")
print("=" * 70)

# All well-established charged mesons from PDG 2024
# Format: (name, mass_MeV, J^PC, quark_content, notes)
mesons = [
    ("pi+-",      139.570, "0-", "ud",   "pseudo-Goldstone"),
    ("K+-",       493.677, "0-", "us",   ""),
    ("D+-",      1869.66,  "0-", "cd",   ""),
    ("D0",       1864.84,  "0-", "cu",   "neutral"),
    ("D*+-",     2010.26,  "1-", "cd",   "vector"),
    ("D_s+-",    1968.35,  "0-", "cs",   ""),
    ("D*_s+-",   2112.2,   "1-", "cs",   "vector"),
    ("D_0*(2300)",2343.,   "0+", "cd",   "scalar, broad"),
    ("D_s0*(2317)",2317.8, "0+", "cs",   "scalar, narrow"),
    # Some less common states near 1800
    ("D(2S)+-",  2621.,    "0-", "cd",   "radial excitation"),
    ("B+-",      5279.34,  "0-", "ub",   ""),
    ("B_c+-",    6274.47,  "0-", "cb",   ""),
    # eta mesons (neutral but relevant for spectrum)
    ("eta_c(1S)", 2983.9,  "0-", "cc",   "neutral"),
    # f0 mesons near 1800
    ("f0(1710)",  1720.,   "0+", "ss/gg","neutral, scalar"),
    ("f0(1770)",  1784.,   "0+", "??",   "neutral, needs confirmation"),
    # a0 mesons
    ("a0(1450)",  1474.,   "0+", "ud",   "scalar"),
    # K mesons
    ("K*(1680)",  1718.,   "1-", "us",   "vector"),
    ("K_0*(1430)",1425.,   "0+", "us",   "scalar"),
    ("K(1830)",   1874.,   "0-", "us",   "needs confirmation"),
]

print(f"\n  {'Meson':<18s} {'mass':>8s} {'J^P':>4s} {'qq':>4s} {'Δm':>8s} {'m^2-m_τ^2':>12s} {'ratio':>8s}")
print("  " + "-" * 72)

for name, mass, jp, qq, note in sorted(mesons, key=lambda x: x[1]):
    if 1500 < mass < 2500:
        dm = mass - m_partner
        dsq = mass**2 - m_tau**2
        r = dsq / founding if founding != 0 else 0
        flag = ""
        if abs(dm) < 80:
            flag = " ← CANDIDATE"
        print(f"  {name:<18s} {mass:>8.1f} {jp:>4s} {qq:>4s} {dm:>+8.1f} {dsq:>12.0f} {r:>8.2f}{flag}")

# =====================================================================
# 2. Detailed analysis of D+- candidate
# =====================================================================
print("\n" + "=" * 70)
print("2. D+- AS TAU PARTNER")
print("=" * 70)

m_D = 1869.66
print(f"\n  Predicted tau partner: {m_partner:.1f} MeV")
print(f"  D+- mass:             {m_D:.1f} MeV")
print(f"  Discrepancy:          {m_D - m_partner:.1f} MeV ({(m_D - m_partner)/m_D*100:.1f}%)")
print()

# If D is the tau partner, what splitting ratio does it imply?
D_split = m_D**2 - m_tau**2
ratio_D = D_split / founding
print(f"  m_D^2 - m_tau^2 = {D_split:.0f} MeV^2")
print(f"  Ratio to founding: {ratio_D:.2f}")
print(f"  Predicted ratio:   {ratio:.2f}")
print(f"  Discrepancy:       {(ratio_D - ratio)/ratio*100:.1f}%")
print()

# Quantum number check
print("  Quantum numbers:")
print(f"    tau:  charge=-1, spin=1/2, lepton number=1, color singlet")
print(f"    D+-:  charge=+1, spin=0,   lepton number=0, color singlet")
print(f"    D-:   charge=-1, spin=0,   lepton number=0, color singlet")
print()
print("  The D- (cd-bar) has the SAME charge as the tau-.")
print("  In the SUSY limit, the tau- (fermion) and D- (scalar)")
print("  would form a chiral multiplet with mass m_+.")
print()
print("  Quark content of D-: c-bar d (or equivalently d c-bar)")
print("  This is an off-diagonal meson M^d_c in ISS language.")
print("  The ISS meson M^d_c carries flavor quantum numbers of")
print("  a d-quark and a c-antiquark.")

# =====================================================================
# 3. What corrections could close the 60 MeV gap?
# =====================================================================
print("\n" + "=" * 70)
print("3. CORRECTIONS TO THE 1810 MeV PREDICTION")
print("=" * 70)

gap = m_D - m_partner
print(f"\n  Gap: m_D - m(pred) = {gap:.1f} MeV ({gap/m_D*100:.1f}%)")
print()

# Correction 1: CW radiative correction
# The CW potential adds positive mass^2 to scalars.
# The mesino CW mass^2 at the Seiberg vacuum is tiny (~10^-9 of tree).
# But the MESON CW mass^2 could be larger.
print("  a) Coleman-Weinberg correction:")
print("     CW adds positive m^2 to scalars, pushing them UP.")
print("     The mesino CW correction is tiny (~10^-9 of tree).")
print("     But the meson CW correction from heavy quark loops")
print("     could be O(alpha_s/pi) * m_c^2 ~ 50-100 MeV^2.")
m_cw_correction = 0.1 * 1275**2  # rough: alpha_s/pi * m_c^2
delta_m_from_cw = m_cw_correction / (2 * m_partner)  # delta_m ≈ delta_m^2 / (2m)
print(f"     Estimate: δm^2 ~ (α_s/π)m_c^2 ~ {m_cw_correction:.0f} MeV^2")
print(f"     → δm ~ {delta_m_from_cw:.0f} MeV")
print()

# Correction 2: Bloom correction to the tau mass
# The O'R seed gives m_tau(seed) = 1472 MeV, physical = 1777 MeV.
# The bloom shifts the tau mass by 21%. If the meson partner
# shifts by a DIFFERENT amount, this modifies the splitting.
print("  b) Bloom correction:")
print("     The bloom shifts the lepton (fermion) and meson (scalar)")
print("     differently if the bloom enters the Kahler potential.")
print("     A 3% shift in the splitting would close the gap.")
print(f"     Need: δ(splitting) = {gap**2 + 2*m_partner*gap:.0f} MeV^2")
print(f"     = {(gap**2 + 2*m_partner*gap)/tau_split*100:.1f}% of tree splitting")
print()

# Correction 3: The founding relation itself is 2% off
# m_pi^2 - m_mu^2 = 8316, f_pi^2 = 8501 → 2.2% deficit
# If we use the EXACT founding relation instead of f_pi^2:
print("  c) Founding relation accuracy:")
print(f"     m_pi^2 - m_mu^2 = {founding:.0f} MeV^2")
print(f"     f_pi^2 = {f_pi**2:.0f} MeV^2")
print(f"     Using the exact founding value (not f_pi^2):")
tau_split_exact = ratio * founding  # already using founding
print(f"     Tau splitting = {ratio:.4f} × {founding:.0f} = {tau_split_exact:.0f} MeV^2")
print(f"     m(partner) = {np.sqrt(m_tau**2 + tau_split_exact):.1f} MeV (same, since we already used founding)")
print()

# Correction 4: Higher-order in the B-term expansion
# The formula uses first-order perturbation theory: delta_m^2 = B * |a_k|^2
# At second order, there are cross-terms between the B-term and
# the fermion mass matrix.
print("  d) Second-order B-term:")
# At first order: m_scalar^2 = m_fermion^2 + B * Phi1_content
# At second order: corrections of order B^2 / (m_+^2 - m_-^2)
y_phys = founding / (2 * np.sin(np.pi/12)**2 * (m_mu/(2-np.sqrt(3)))**2)
B_param = 2 * y_phys  # in O'R units
dm2_gap = m_D**2 - m_partner**2  # the gap in m^2
print(f"     y = {y_phys:.4f}")
print(f"     Second-order correction ~ y^2 * m^2 = {y_phys**2:.6f} × m^2")
m_OR = m_mu / (2 - np.sqrt(3))
second_order = y_phys**2 * m_OR**2
print(f"     ~ {second_order:.0f} MeV^2")
print(f"     → δm ~ {second_order/(2*m_partner):.0f} MeV")
print(f"     Gap to fill: {gap:.0f} MeV → need {dm2_gap:.0f} MeV^2")
print(f"     Second-order gives {second_order:.0f} MeV^2 — {'sufficient' if second_order > dm2_gap else 'insufficient'}")
print()

# =====================================================================
# 4. Alternative: D_s as tau partner
# =====================================================================
print("=" * 70)
print("4. ALTERNATIVE: D_s AS TAU PARTNER")
print("=" * 70)

m_Ds = 1968.35
Ds_split = m_Ds**2 - m_tau**2
ratio_Ds = Ds_split / founding
print(f"\n  m_Ds = {m_Ds:.2f} MeV")
print(f"  m_Ds^2 - m_tau^2 = {Ds_split:.0f} MeV^2")
print(f"  Ratio to founding: {ratio_Ds:.2f}")
print(f"  Predicted ratio:   {ratio:.2f}")
print(f"  Discrepancy:       {(ratio_Ds - ratio)/ratio*100:.0f}%")
print()
print("  The D_s is 158 MeV ABOVE the prediction (8.7%).")
print("  Much worse than D+- (60 MeV, 3.2%).")
print()
print("  Quark content: D_s- = c-bar s")
print("  This is meson M^s_c — strangeness-charm sector.")
print("  The D- (M^d_c) is closer to the prediction.")

# =====================================================================
# 5. Alternative: f0(1770) or K(1830) as tau partner
# =====================================================================
print("\n" + "=" * 70)
print("5. NEUTRAL/STRANGE CANDIDATES NEAR 1810 MeV")
print("=" * 70)
print()
print("  f0(1770): mass = 1784 MeV, J^P = 0+, neutral")
print(f"    m^2 - m_tau^2 = {1784**2 - m_tau**2:.0f} MeV^2")
print(f"    Ratio: {(1784**2 - m_tau**2)/founding:.2f}")
print(f"    Discrepancy from pred: {1784 - m_partner:.0f} MeV ({(1784-m_partner)/m_partner*100:.1f}%)")
print("    But: neutral, so cannot pair with charged tau")
print()
print("  K(1830): mass = 1874 MeV, J^P = 0-, charged (us-bar)")
print(f"    m^2 - m_tau^2 = {1874**2 - m_tau**2:.0f} MeV^2")
print(f"    Ratio: {(1874**2 - m_tau**2)/founding:.2f}")
print(f"    Discrepancy from pred: {1874 - m_partner:.0f} MeV ({(1874-m_partner)/m_partner*100:.1f}%)")
print("    But: 'needs confirmation' (PDG), broad state")

# =====================================================================
# 6. The electron partner
# =====================================================================
print("\n" + "=" * 70)
print("6. THE ELECTRON PARTNER")
print("=" * 70)
print()

# The electron is the goldstino multiplet: mass 0 at tree level,
# bloom gives 0.511 MeV. Its scalar partner is the pseudo-modulus X.
# X gets CW mass ~ 13 MeV.
# In the physical spectrum, what meson is at ~13 MeV?
# NOTHING. The lightest meson is pi0 at 135 MeV.
#
# But the prediction isn't m(partner) ~ m_e ~ 0.5 MeV.
# It's m(partner)^2 = m_e^2 + B_goldstino where B_goldstino ~ 0
# (the goldstino direction has no Phi1 content).
# So the electron partner mass ~ m_e ~ 0.5 MeV (barely shifted).
#
# Wait — the goldstino is psi_X, which has NO Phi1 content (it's the
# X direction, not Phi1 or Phi3). So its B-term splitting is ZERO.
# Its partner scalar has the same mass: m_X_scalar ~ 0 + CW.
#
# There is no 0.5 MeV meson. The electron partner is NOT a standard meson.
# It's the pseudo-modulus X, which in the ISS is det(M)/Lambda^3.
# This is a SINGLET under SU(5)_f, not an adjoint meson.

print("The electron sits in the goldstino multiplet (psi_X).")
print("Its scalar partner is the pseudo-modulus X.")
print("B-term splitting for this multiplet: ZERO (no Phi1 content).")
print()
print("The pseudo-modulus gets a CW mass ~ 13 MeV (one-loop).")
print("The electron gets a bloom mass ~ 0.5 MeV.")
print("These are DIFFERENT: the multiplet is strongly broken.")
print()
print("X = det(M)/Lambda^3 in the ISS — a SINGLET, not a meson.")
print("There is no 0.5 MeV charged meson. The electron partner")
print("is not a state in the observed meson spectrum.")
print()
print("This is consistent: the goldstino multiplet is the MOST")
print("broken, absorbing most of the SUSY breaking. Its members")
print("are the farthest from degeneracy.")

# =====================================================================
# 7. Prediction table
# =====================================================================
print("\n" + "=" * 70)
print("7. PREDICTION TABLE")
print("=" * 70)
print()

# The three O'R multiplets:
# Multiplet X (goldstino): electron ↔ pseudo-modulus
# Multiplet - (light): muon ↔ pion (CONFIRMED by founding relation)
# Multiplet + (heavy): tau ↔ D meson (PREDICTED)

print(f"  {'Multiplet':<12s} {'Fermion':<12s} {'m_f (MeV)':<12s} {'Scalar':<12s} {'m_s (MeV)':<12s} {'m_s^2-m_f^2':>14s} {'Status':<12s}")
print("  " + "-" * 86)
print(f"  {'X (0)':<12s} {'electron':<12s} {'0.511':<12s} {'X (singlet)':<12s} {'~13 (CW)':<12s} {'~170':>14s} {'broken':<12s}")
print(f"  {'- (m_-)':<12s} {'muon':<12s} {m_mu:<12.3f} {'pi+-':<12s} {m_pi:<12.3f} {m_pi**2-m_mu**2:>14.0f} {'CONFIRMED':<12s}")
print(f"  {'+ (m_+)':<12s} {'tau':<12s} {m_tau:<12.3f} {'D+- (?)':<12s} {m_D:<12.3f} {m_D**2-m_tau**2:>14.0f} {'predicted':<12s}")
print()

# Ratios
print(f"  Splitting ratios (normalised to muon-pion):")
print(f"    (m_pi^2 - m_mu^2) / (m_pi^2 - m_mu^2) = 1.000  [reference]")
print(f"    (m_D^2 - m_tau^2) / (m_pi^2 - m_mu^2) = {(m_D**2-m_tau**2)/founding:.2f}")
print(f"    Predicted ratio: (2+√3)^2 = {ratio:.2f}")
print(f"    Match: {abs((m_D**2-m_tau**2)/founding - ratio)/ratio*100:.1f}% off")
print()

# What if we use D0 instead of D+-?
m_D0 = 1864.84
print(f"  With D0 (neutral, {m_D0:.2f} MeV):")
print(f"    m_D0^2 - m_tau^2 = {m_D0**2 - m_tau**2:.0f} MeV^2")
print(f"    Ratio: {(m_D0**2 - m_tau**2)/founding:.2f}")
print(f"    Match: {abs((m_D0**2-m_tau**2)/founding - ratio)/ratio*100:.1f}% off")
print()

# =====================================================================
# 8. The 3% discrepancy: what it means
# =====================================================================
print("=" * 70)
print("8. THE 3.2% DISCREPANCY")
print("=" * 70)
print()
print(f"  Predicted: {m_partner:.1f} MeV")
print(f"  D+- mass:  {m_D:.2f} MeV")
print(f"  Gap:       {m_D - m_partner:.1f} MeV")
print()
print("  Possible interpretations:")
print("  A) The D meson IS the tau partner, and the 3% comes from")
print("     radiative corrections (CW, strong coupling, etc.)")
print("  B) The D meson is NOT the tau partner, and the true partner")
print("     is an as-yet undiscovered state at 1810 MeV.")
print("  C) The founding relation m_pi^2 - m_mu^2 = f_pi^2 is not")
print("     exact, and the actual B-term ratio differs from (2+sqrt3)^2.")
print()
print("  In favor of A:")
print("    - 3% is within typical radiative correction range")
print("    - D+- has the right charge (-1 for D- paired with tau-)")
print("    - D+- is J^P = 0- (pseudoscalar, matching the O'R scalar)")
print(f"    - Second-order B-term contributes ~{second_order:.0f} MeV^2 → {second_order/(2*m_partner):.0f} MeV")
print(f"    - CW correction could add ~{delta_m_from_cw:.0f} MeV")
print()
print("  In favor of B:")
print("    - A scalar meson at 1810 MeV with J^P = 0+ could exist")
print("      but be hard to detect (broad, overlapping with D)")
print()

# Final check: is the gap CLOSED by second-order + CW?
total_correction_m2 = second_order + m_cw_correction
corrected_partner = np.sqrt(m_tau**2 + tau_split + total_correction_m2)
print(f"  Including corrections:")
print(f"    Second-order B-term:  +{second_order:.0f} MeV^2")
print(f"    CW (rough estimate):  +{m_cw_correction:.0f} MeV^2")
print(f"    Corrected prediction: {corrected_partner:.1f} MeV")
print(f"    New gap to D+-:       {m_D - corrected_partner:.1f} MeV ({abs(m_D - corrected_partner)/m_D*100:.1f}%)")

# =====================================================================
# 9. REALITY CHECK: the m^2 splitting is NOT close
# =====================================================================
print("\n" + "=" * 70)
print("9. REALITY CHECK: m^2 vs mass")
print("=" * 70)
print()
print("WARNING: The 3% mass match is MISLEADING.")
print("The m^2 splittings are very different:")
print(f"  Predicted:  m^2(partner) - m^2(tau) = {tau_split:.0f} MeV^2")
print(f"  If D+-:     m^2(D) - m^2(tau)       = {m_D**2-m_tau**2:.0f} MeV^2")
print(f"  Ratio: {(m_D**2-m_tau**2)/tau_split:.2f} (should be 1.00)")
print()
print("The masses are close (3%) only because m_tau^2 DOMINATES:")
print(f"  m_tau^2 = {m_tau**2:.0f} MeV^2 (97% of m_D^2)")
print(f"  The splitting is only {tau_split/m_tau**2*100:.1f}% of m_tau^2.")
print(f"  ANY meson within ~200 MeV of m_tau gives a mass 'match' <10%.")
print()
print("For comparison, the founding relation is a STRONG test:")
print(f"  m_mu^2 = {m_mu**2:.0f} MeV^2")
print(f"  m_pi^2 - m_mu^2 = {founding:.0f} MeV^2 = {founding/m_mu**2*100:.0f}% of m_mu^2")
print(f"  The splitting is comparable to the mass — a genuine prediction.")
print()

# What the B-term ACTUALLY predicts for the tau sector:
# The tau partner should be at m_tau + delta, where delta ~ 33 MeV
# (from tau_split / (2*m_tau) ≈ 115829 / 3554 = 32.6 MeV)
delta_m_pred = tau_split / (2 * m_tau)
print(f"The B-term predicts m(partner) = m_tau + {delta_m_pred:.1f} MeV = {m_tau + delta_m_pred:.1f} MeV")
print(f"The D+- is at m_tau + {m_D - m_tau:.1f} MeV.")
print(f"The prediction is that the partner is {delta_m_pred:.0f} MeV above the tau.")
print(f"The D is {m_D - m_tau:.0f} MeV above the tau.")
print(f"Ratio: {(m_D - m_tau)/delta_m_pred:.2f}")
print()
print("The D is 2.8× farther from the tau than predicted.")
print("This is NOT a 3% discrepancy — it's a factor-3 failure")
print("in the SPLITTING, masked by the dominance of m_tau^2.")
print()
print("HONEST ASSESSMENT:")
print("  The (2+√3)^2 ratio predicts a tau partner 33 MeV above m_tau.")
print(f"  The D+- is 93 MeV above m_tau.")
print(f"  The meson f0(1770) at 1784 MeV is only 7 MeV above m_tau")
print(f"  — much closer to the predicted 33 MeV shift.")
print(f"  But f0(1770) is neutral and J^P = 0+, not 0-.")
