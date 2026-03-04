# Assembly Round 18: Paper Split Preparation

## Strategy

Round 18 (even round) prepared for the paper split:
1. **O'R mass relations in Lagrangian language** — vocabulary for Paper A (18A, opus)
2. **AMSB soft spectrum** — can AMSB lift mesinos? (18B, opus)
3. **Paper section audit** — classify every section for the split (18C, opus)
4. **Complete lepton Lagrangian** — explicit SU(2) superpotential (18D, opus)

## Agents

- **18A** (opus): O'R mass relations — COMPLETE ✅
- **18B** (opus): AMSB spectrum — COMPLETE ✅
- **18C** (opus): Paper section audit — COMPLETE ✅
- **18D** (opus): Lepton Lagrangian — COMPLETE ✅

---

## Agent Results

### 18A: O'Raifeartaigh Mass Relations ⭐⭐⭐

**Result: Complete Lagrangian-language derivation of all five mass predictions.**

The document (or_mass_relations.md) presents:
1. **O'Raifeartaigh mass ratio**: m_+/m_- = (2+√3)² ≈ 13.93 (structural constant of the vacuum)
2. **Q = 2/3 structural invariant**: algebraic identity, proved from AB = 1, A²+B² = 4, (A+B)² = 6
3. **Bion mass relation**: √m_b = 3√m_s + √m_c. Preserves Q ≈ 2/3 (Q(-s,c,b) = 0.6727, 0.90% from 2/3)
4. **Yukawa eigenvalue constraint**: Solving Q(c,b,t) = 2/3 for m_t. From predicted inputs: m_t = 171,250 MeV (−0.9%)
5. **The m_s chain**: m_s → m_c (+2.0%) → m_b (+1.3%) → m_t (−0.9%)
6. **Lepton sector**: Q(e,μ,τ) = 2/3 predicts m_τ = 1776.97 MeV (+0.006%)

**Key for Paper A**: All mass predictions are presented as consequences of the superpotential structure, not as empirical formulas. The vocabulary is "O'Raifeartaigh mass ratio", "bion mass relation", "Yukawa eigenvalue constraint".

### 18B: AMSB Soft Spectrum ⭐⭐ (NEGATIVE RESULT)

**Result: AMSB scalar masses are TACHYONIC. Mesino fermion masses are NOT lifted.**

Key findings:

1. **Scalar meson soft masses**: m²_AMSB = −(1/4) γ̇_M m²_{3/2}. At α_s = 1, m_{3/2} = 50 TeV: m² = −(15,915 GeV)². **Tachyonic** (the known slepton problem).

2. **Fermion masses protected by holomorphy**: AMSB generates soft terms (scalar m², B-terms, A-terms) but does NOT generate fermion mass terms. The mesino masses |X₀|M_j remain at 11–494 eV. The claim δm_mesino ~ γ_M × m_{3/2} was WRONG — it confused B-term with fermion mass.

3. **Vacuum destabilization**: The tachyonic m² overwhelms f_π² by factor 10¹⁰, destroying the seesaw vacuum. At the new vacuum (incalculable), all masses reorganize at the AMSB scale.

4. **Implication**: AMSB cannot simply "lift" mesinos to LEP safety. It destabilizes the vacuum entirely. The mesino mass problem remains open (requires non-perturbative Kähler or new mechanism).

### 18C: Paper Section Audit ⭐⭐⭐ CRITICAL

**Result: Complete classification of all 45 sections/subsections. 11 sections must be split.**

Paper A (Lagrangian) gets:
- Composite scalar counting + Diophantine uniqueness (Sec. 8)
- Diquark sector + Pauli theorem (Sec. 9)
- SO(32) embedding (Sec. 10 Part A)
- Nearly all of SUSY breaking section (Sec. 11): O'R mechanism, Kähler, bion, Lagrangian, predictions
- Most Open Problems (diquark, bloom, lepton, baryon, FCNC, Higgs)

Paper B (Koide evidence) gets:
- Koide formula presentation + statistical significance (Sec. 2)
- Empirical seed triples (Sec. 3.2–3.5)
- Meson Koide triple (Sec. 5)
- Mass relations and iterative chains (Sec. 7)
- Look-elsewhere corrections, Monte Carlo
- δ₀ mod 2π/3 = 2/9 observation

**Most pervasive modification**: Replace "Koide" → "energy-balance condition/O'Raifeartaigh mass ratio" throughout Paper A. 12 specific paragraphs identified.

**Content gaps for Paper A**:
1. Self-contained forward derivation of Q = 2/3 from the superpotential (currently scattered)
2. Clean linear presentation of m_s → m_c → m_b → m_t chain (currently fragmented)
3. Spectrum table as centerpiece (currently buried at line 2978)
4. tan β = 1 prediction needs more prominence

### 18D: Complete Lepton Lagrangian ⭐⭐⭐

**Result: Full SU(2) SQCD Lagrangian with all couplings specified.**

Key findings:

1. **Superpotential**: W = f X_L + m L₁L₃ + g X_L L₁² + ε L₁L₂L₃ with ε = 1/Λ_L³

2. **CW does NOT select gv/m = √3**: The Coleman-Weinberg potential has minimum at v_min ~ 0.3–0.5, far below √3. The Koide condition requires a tree-level non-canonical Kähler: K = −|X_L|⁴/(12μ_K²). Self-consistency: μ_K = √2 m/g.

3. **SM quantum numbers**: UV preons with Y = −1/2 (SU(2)_W singlets) → composites with Y = −1 (right-handed charged leptons). SM coupling to doublets and Higgs is an open question.

4. **Parameter count**: 5 total (m, g, f, ε, δ_bloom), but only 2 determine physical masses: m_OR (scale) and δ_bloom (electron mass).

5. **Goldstino = ψ_{X_L}** (purely). ψ_{L₂} is the pseudo-modulus partner, massless at tree level.

---

## Synthesis

### The picture after Round 18

Round 18 delivers the vocabulary and structure for the paper split:

| Component | Status | Paper |
|-----------|--------|-------|
| O'R mass ratio language | Ready | A |
| Bion mass relation language | Ready | A |
| Yukawa constraint language | Ready | A |
| Paper section classification | Complete (45 sections) | Both |
| Lepton Lagrangian | Complete | A |
| AMSB soft spectrum | Computed — tachyonic (negative) | Open problem |
| Kähler stabilization of gv/m = √3 | Requires non-canonical Kähler | Open problem |

### New negative results

| Result | Source | Impact |
|--------|--------|--------|
| AMSB scalar masses tachyonic | 18B | Cannot simply lift mesinos; vacuum destabilizes |
| Mesino masses protected by holomorphy | 18B | No SUSY-breaking correction to fermion masses |
| CW does NOT select gv/m = √3 | 18D | Kähler origin needed (assumed, not derived) |

### Paper A structure (from audit)

1. Introduction (rewrite: Lagrangian motivation)
2. Bootstrap uniqueness (N=3, r=3, s=2) — Sec. 8
3. SU(5) flavor and diquarks — Sec. 9
4. SO(32) embedding — Sec. 10A
5. SUSY breaking: O'Raifeartaigh + seesaw — Sec. 11 (rewrite)
6. The mass prediction chain: m_s → m_c → m_b → m_t (NEW centerpiece)
7. Lepton sector: SU(2) O'Raifeartaigh — from 18D
8. Complete particle spectrum — from 17B
9. Open problems and FCNC
10. Conclusions

### Paper B structure (from audit)

1. Introduction (the Koide formula as empirical observation)
2. The energy-balance condition Q = 2/3
3. Seed triples: (0,s,c), (0,π,D_s), leptons
4. The quark triple (−s,c,b) and meson triple (−π,D_s,B)
5. The heavy triple (c,b,t)
6. Mass relations and iterative chains
7. Statistical assessment and look-elsewhere
8. The de Vries angle and Casimir equation
9. Connection to the Lagrangian (summary of Paper A results)
10. Conclusions

### Priorities for Round 19

1. **Draft Paper A abstract and introduction** — in Lagrangian language
2. **Build the prediction chain section** — the centerpiece of Paper A
3. **Draft Paper B introduction** — empirical framing
4. **Resolve the open Kähler problem** — can we say anything about K = −|X|⁴/(12μ²)?

---

*Generated: 2026-03-04*
*Based on: or_mass_relations.md, amsb_spectrum.md, paper_audit.md, lepton_lagrangian.md*
