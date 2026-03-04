# Assembly Round 11: μ-term, FCNC, CKM Quantitative, Lepton Alternatives

## Strategy

After Round 10, the most critical unresolved issues are:

1. **μ-term / X normalization** — X has [mass]^{-3}, making λXH_u·H_d dimensionally inconsistent as written. Need: rescale X̂ = X·Λ⁴ → canonical NMSSM singlet, rederive all couplings.
2. **FCNC constraints** — Mesinos at QCD scale mediate tree-level FCNC. Bounds from K-K̄ mixing, D-D̄ mixing, B-B̄ mixing could rule out the framework.
3. **CKM quantitative** — Agent 10D found θ₁₂ = 2.6° (vs 13.04° PDG) from saddle point. Need full 32-variable minimization or analytic treatment.
4. **Lepton sector alternatives** — Kähler approach exhausted. What else? Separate sector? Non-renormalizable W? UV boundary condition?

All agents use opus (odd round rule).

## Agents

- **11A** (opus): μ-term dimensional analysis — COMPLETE ✅
- **11B** (opus): FCNC constraints — COMPLETE ✅
- **11C** (opus): CKM full minimization — RUNNING
- **11D** (opus): Lepton mass alternatives — COMPLETE ✅

---

## Agent Results

### 11A: μ-term Dimensional Analysis ⭐

The rescaling X → Ŝ = X·Λ⁴ is a **pure field redefinition** that changes no physical observable. The scalar potential, masses, mixing angles are all invariant.

Key results:
- ⟨Ŝ⟩ = -(m_u m_d m_s)^{1/3} = -9.803 MeV
- μ_eff = λ⟨Ŝ⟩ = 0.72 × (-9.803) = **-7.058 MeV** (3 orders below ~100 GeV needed for EWSB)
- λ_new = 0.72 (dimensionless) in both bases
- m_h = λv/√2 = 125.35 GeV — unchanged
- F_X/Λ⁴ = F_Ŝ, but K^{ŜŜ̄} = Λ⁸ compensates → V is exactly invariant

**The μ-problem is physical**: ⟨S⟩ is set by the geometric mean of light quark masses. Model needs separate μ-source (Giudice-Masiero, radiative generation, or X ≠ NMSSM singlet).

### 11B: FCNC Constraints

**Reframed after user observation**: Off-diagonal mesons M^d_s ARE the SM kaons — they're not new d.o.f. Their K-K̄ contribution is standard QCD, not new physics. The question is whether SUSY couplings (λXH·H, Yukawa-meson) generate extra FCNC.

Original agent analysis (treating as new scalars):
- K-K̄: m_S = f_π too light by 8× for g = 1
- D-D̄: too light by 12×
- B-B̄: marginal

But this overcounts because mesons ARE the SM mesons. The real FCNC analysis requires:
- Computing the Yukawa-Higgs mediated flavour violation
- Checking for new box diagrams from the NMSSM coupling
- Neither introduces elementary scalar exchange at tree level (the mesons are composite)

### 11D: Lepton Mass Alternatives ⭐

Six mechanisms surveyed. Results:

| Mechanism | Q = 2/3? | Hierarchy? | Status |
|-----------|----------|------------|--------|
| Sp(2) N_f=3 sector | Yes (O'R-Koide) | Yes (bloom) | **Most promising** |
| UV boundary at Λ~300 MeV | Preserved (QED) | Deferred | Complementary |
| Non-renormalizable ops | Tuning only | Poor | Ruled out |
| Radiative | No | Qualitative | Unlikely |
| Lepton-meson duality | Wrong quarks | No | Ruled out |
| Duality cascade | No (Q~1/3) | No | Ruled out |

**Sp(2) proposal**:
- Sp(2) SQCD with N_f = 3 antisymmetric flavors → 3 Pfaffian eigenvalues
- O'Raifeartaigh-Koide mechanism (gv/m = √3) generates Q = 2/3 at seed
- Bloom rotation handles hierarchy; Λ_L ~ 300 MeV → M_0 = 313.8 MeV
- QED running preserves Q to O(α²/π²) ~ 10⁻⁵
- v₀(lepton)/v₀(seed) = 1.014 (NOT 2.0; different bion dynamics for Sp(2))
- Next steps: verify Sp(2) vacuum structure, compute bion potential, specify mediation

### 11C: CKM Full Minimization & Cabibbo-Oakes Analysis

Script ran 32-variable minimization (9 complex meson entries + X + B + B̃ + H_u + H_d). Partial completion (10/30 random restarts before timeout).

From the diagonal seesaw starting point:
- V_start = 4.76 × 10²⁰ MeV² (dominated by |F_X|² = 4.76 × 10²⁰)
- After minimization from seesaw: V = 1.71 × 10¹⁵ MeV² (soft-term dominated)

**Cabibbo-Oakes analysis** (cabibbo_oakes.md):
- Weinberg-Oakes: tan(θ_C) = √(m_d/m_s) → θ_C = 12.60° (PDG: 13.04°, -3.3%)
- GST: sin(θ_C) = √(m_d/m_s) → θ_C = 12.92° (PDG: 13.04°, **-0.9%**)
- Fritzsch correction worsens match (subtracts √(m_u/m_c) = 0.041)
- **GST is the best variant at -0.9%**
- Koide does NOT hold for (u,d,s): Q = 0.567, -14.9% from 2/3
- Parameter counting: 6 masses, 4+1 Koide constraints + 1 external (y_t) → **2 free parameters (m_u, m_d)**
- Oakes connects m_d/m_s to θ_C, so with θ_C as input, only m_u remains free

---

## Synthesis

### Resolved in Round 11:
1. ✅ **μ-problem is physical** — μ_eff = -7 MeV, needs separate source
2. ✅ **FCNC reframed** — mesons are SM particles, not extra d.o.f.
3. ✅ **Lepton sector path identified** — Sp(2) O'R-Koide is viable

### Round 12 priorities (even round = random sonnet/opus):
1. **Sp(2) vacuum verification** — Does N_f = 3 Sp(2) SQCD produce O'Raifeartaigh vacuum?
2. **CKM analytic** — Derive Oakes relation analytically from the tachyonic condensate structure
3. **Paper polishing** — Clean up open problems section, add Sp(2) details, prepare for submission
4. **Brainstorm** — What remains before the North Star is reached?
