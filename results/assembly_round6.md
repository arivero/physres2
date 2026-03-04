# Assembly Round 6: Kähler Origin, Vacuum Self-Consistency, Yukawa-Koide

## Agent A: One-Loop CW Kähler (kahler_oneloop.md)

**Task**: Compute the Coleman-Weinberg effective potential for the O'Raifeartaigh pseudo-modulus and extract the quartic Kähler coefficient. Does one-loop CW generate c = -1/12?

**Key results**:

1. **c₄(y=1) = -20/3 exactly** (verified to 30 decimal places via mpmath). In units of g⁴/(16π²), the physical quartic coefficient C₄ = -5/3.

2. **c₄ is negative for y < y_cross = 1.7126.** The CW potential does produce a negative quartic in the relevant regime. However:

3. **c₄ ≠ -g²/12 in any normalization.** The tree-level Kähler coefficient c = -g²/12 involves f² (through V = f²/K_{XX̄}), giving a different parametric dependence than the one-loop CW. The effective Kähler coefficient extracted from c₂ gives c_eff ≈ -5 to -6, far from -1/12 = -0.083.

4. **STr[M⁴] = 8y² (constant in v)**, confirming the -3/2 constant doesn't affect v-dependence.

**Conclusion**: The Kähler pole at v = √3 is a **tree-level statement** about the Kähler cutoff structure, not a one-loop consequence. The CW potential stabilizes the origin (c₂ > 0), and while c₄ is negative, it's orders of magnitude too large to match c = -1/12.

---

## Agent B: Vacuum Structure (vacuum_structure.md)

**Task**: Derive F-term equations for the full Lagrangian with three Kähler terms. Check self-consistency of the Seiberg vacuum.

**Key results**:

1. **F_X = 0 at the Seiberg vacuum** (det M₃ = Λ⁶ identically). SUSY is unbroken in the confined (u,d,s) sector.

2. **Yukawa F-terms break SUSY**: F_{H_u} = y_c M_cc ≠ 0, F_{H_d} = y_b M_bb ≠ 0. The key insight: y_j M_j = (m_j/v)(C/m_j) = C/v, so **F-terms are flavor-universal**: F = C/v_EW = 289.6 MeV.

3. **Kähler pole condition determines μ**: |X| = √3 μ combined with the seesaw X = -C/Λ⁶ gives μ = |C/(√3 Λ⁶)|. This is a consistency relation, not a tuning.

4. **Bion corrections are perturbatively small**: O(ζ²) = O(1.5 × 10⁻²) for the diagonal metric, O(10⁻⁵) for strange meson.

5. **Parameter budget**: The vacuum is fully determined by **(m_u, m_d, m_s, Λ)** with m_c and m_b predicted.

6. **Goldstino direction** (with Yukawa SUSY breaking): ψ_G ∝ (C/v)(ψ_{H_u} + ψ_{H_d}) — equal-weight Higgsino combination.

**Conclusion**: Self-consistent vacuum. Three Kähler terms play distinct roles: canonical kinetics, X stabilization via pole, bion v₀-doubling.

---

## Agent C: Yukawa-Koide at tan β = 1 (yukawa_koide.md)

**Task**: Compute Yukawa coupling matrices from the Koide parametrization with tan β = 1 (v_u = v_d = v/√2).

**Key results**:

1. **Koide phase transfers identically to Yukawas**: y₀ = 2M₀/v, δ unchanged. Q is invariant under the rescaling m → y = 2m/v.

2. **Exact Cartan identity**: θ = π/6 - δ, where θ is the angle of the Koide vector in the (λ₃, λ₈) Cartan plane.

3. **Bloom = rigid Cartan rotation**: Shifting δ by Δδ rotates (n₁, n₂) by -Δδ in the Cartan plane. The norm ||n|| is fixed at √6/2 for Q = 2/3, and the trace (3√y₀) is fixed by M₀.

4. **Top Yukawa y_t = 1.403** at tan β = 1 — O(1) as expected for near-infrared fixed point.

5. **Bloom rotations confirmed**: -2.27° (leptons), +22.21° (-s,c,b), -11.07° (c,b,t).

**Conclusion**: The Koide parametrization maps cleanly to the Yukawa sector. Bloom is geometrically a Cartan plane rotation. The exact identity θ = π/6 - δ connects the Koide phase directly to flavor-space geometry.

---

## Agent D: ISS One-Loop Kähler (iss_kahler_oneloop.md)

**Task**: Compute one-loop CW corrections from ISS magnetic quarks. Does the ISS mechanism generate a negative quartic Kähler correction?

**Key results**:

1. **ISS CW is monotonically increasing** for x > 0, confirming X = 0 as the global CW minimum.

2. **c_eff is POSITIVE**: +0.014 at x = 1, vs. the target -1/12 = -0.083. The ISS perturbative loop produces the **wrong sign**.

3. **Quartic subtraction cannot create minimum at √3**: The CW potential is concave at x = √3, so any quartic correction gives a maximum, not a minimum.

4. **Kähler metric route**: Moving the minimum to x = √3 through the Kähler metric requires c ≈ +2.7 — unnaturally large.

5. **Supertrace identities**: STr[M²] = 0, STr[M⁴] = 2N_c(hμ)⁴ = const. Standard results confirmed.

**Conclusion**: The ISS perturbative CW mechanism robustly stabilizes X = 0 and produces a positive quartic correction. An external **non-perturbative** contribution (e.g., bion Kähler potential from center-symmetric compactification) is needed to shift the vacuum away from the origin.

---

## Synthesis: Physical Origin of c = -1/12

The two Kähler agents (A and D) converge on the same conclusion through different routes:

| Approach | Quartic sign | Magnitude | c = -1/12? |
|----------|-------------|-----------|------------|
| O'Raifeartaigh CW | negative | c₄ = -20/3 (too large) | No |
| ISS CW | positive | c_eff = +0.014 | No (wrong sign) |
| Tree-level Kähler pole | negative | c = -1/12 (exact) | By construction |

**The Kähler pole mechanism is tree-level.** The coefficient c = -1/12 must come from the UV completion or non-perturbative effects — not from integrating out loops. This sharpens the open question: what generates the tree-level Kähler structure K = |X|² - |X|⁴/(12μ²)?

Possible origins:
- **Bion-induced Kähler**: Center-symmetric compactification generates non-perturbative Kähler corrections of the form K ~ exp(-S₀)|X|⁴/Λ². The sign depends on the topological sector.
- **UV embedding**: The Kähler metric singularity at |X| = √3 μ may signal the breakdown of the low-energy description, analogous to the Landau pole. The UV theory (electric SQCD?) would resolve the singularity.
- **Composite structure**: If X = Φ⁴₄ is a confined meson, its Kähler potential is determined by the confining dynamics, not by canonical normalization.

## Key Numbers for the Paper

| Quantity | Value | Status |
|----------|-------|--------|
| c = -1/12 | exact tree-level Kähler | origin: non-perturbative/UV |
| F = C/v_EW | 289.6 MeV | SUSY breaking scale (universal) |
| y_t (tan β = 1) | 1.403 | O(1), near-infrared FP |
| θ = π/6 - δ | exact identity | Cartan-Koide connection |
| c₄(y=1) | -20/3 exactly | CW quartic (O'Raifeartaigh) |
| c_eff(ISS, x=1) | +0.014 | positive, wrong sign |
| Parameters | (m_u, m_d, m_s, Λ) | 4 free, m_c and m_b predicted |

## Status of Open Problems After Round 6

1. **Physical origin of c = -1/12**: SHARPENED — not radiative. Needs non-perturbative mechanism or UV completion. ★
2. **CKM mixing**: Still open from Round 5. First Oakes relation works; higher orders need new mechanism.
3. **15 vs 10̄ Pauli conflict**: Unchanged. SU(6) embedding separates them into different reps.
4. **Bloom mechanism**: Geometrically clarified (Cartan plane rotation), but dynamical origin still open.
5. **Vacuum structure**: RESOLVED — self-consistent with flavor-universal SUSY breaking F = C/v_EW.
6. **Yukawa-Koide transfer**: RESOLVED — exact identity θ = π/6 - δ connects Koide phase to Cartan geometry.
