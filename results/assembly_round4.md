# Assembly Round 4: Agent Synthesis

## Agent Results

### A. SQCD + Yukawa Spectrum (sqcd_yukawa_spectrum.md)

Full 13×13 fermion mass matrix at the Seiberg vacuum (N_f = N_c = 3, u,d,s).

**Key findings:**
1. At y=0, spectrum reproduces earlier sqcd_spectrum.md results exactly
2. Adding Yukawa coupling y₃ to heaviest flavor lifts the H zero mode at O(y²)
3. Scale hierarchy:
   - Heavy pair: ±7.73 × 10¹⁰ MeV (from cofactor M_i M_j terms)
   - Light modes: ~10⁻⁵ to 10⁻⁴ MeV (from |X| M_k terms)
   - Baryon: ±1.21 × 10⁻⁹ MeV (very small)
4. Characteristic central scale: √C = √(Λ²(m_u m_d m_s)^{1/3}) ≈ 939.3 MeV
   - Compare proton mass 938.3 MeV — coincidence from Λ ≈ Λ_QCD and dimensional transmutation
5. No Koide triples found in the eigenvalue spectrum (best Q ≈ 0.40)
   - Expected: this is the (u,d,s) sector with Q(m_u,m_d,m_s) = 0.567

**For the paper:** The Yukawa perturbation analysis shows how the Higgs mode emerges from the confining vacuum. The spectrum structure (seesaw + perturbative Yukawa) is consistent with the staged picture.

### B. Koide-CKM Mixing (koide_ckm.md)

Rotation matrix from overlapping Koide eigenbases for (-s,c,b) and (c,b,t).

**Key findings:**
1. Oakes first relation confirmed: |V_us| ≈ √(m_d/m_s) = 0.2236 vs 0.2243 (-0.3%)
2. Higher Oakes relations FAIL:
   - |V_cb| ≈ √(m_s/m_b) = 0.1495 vs 0.0422 (ratio 3.54)
   - |V_ub| ≈ √(m_d/m_b) = 0.0334 vs 0.00394 (ratio 8.48)
3. Direct Koide eigenbasis overlap gives CIRCULANT mixing matrix (all angles ~90°)
   - This is because the Koide parametrization produces democratic/circulant eigenvectors
   - The overlap R is not unitary (det R = 0.844)
4. Koide direction angle between triples: 23.36° (not Cabibbo)
5. Cartan-plane projected angle: 33.28° (not Cabibbo either)

**NEGATIVE RESULT:** CKM mixing does NOT emerge from direct overlap of Koide eigenbases. The Koide parametrization's circulant structure makes the eigenvectors too symmetric. The CKM hierarchy (small off-diagonal elements) requires an additional mechanism beyond the Koide condition.

**For the paper:** This sharpens the CKM open problem. The first Oakes relation (|V_us| ≈ √(m_d/m_s)) works at 0.3% — this is the free-parameter connection noted in the memory. Higher-order Oakes relations need additional suppression factors. The paper should state this honestly as an open problem.

### C. Pseudo-modulus (incomplete)

Agent produced pseudomodulus_vev.py but no writeup. Needs retry.

## Gap Status Update

Incorporating these results into the gap analysis:

- **CKM mixing**: Now definitively an open problem. First Oakes relation confirmed. Higher relations need new mechanism. Koide eigenbasis overlap doesn't give CKM.
- **SQCD+Yukawa spectrum**: Complete. Shows how Higgs mode emerges perturbatively from confining vacuum.
- **Pseudo-modulus stabilization**: Still open. Needs retry.
- **Staged symmetry breaking (Unfolding)**: New conceptual framework from user. Needs computational verification.

### D. SU(6) Branching (su6_branching.md)

Full branching rules for SU(6) representations under SU(5)×U(1) and SU(3)×SU(3)×U(1).

**Key structural finding for the sBootstrap:**

The sBootstrap uses 24 + 15 + 15̄ under SU(5). Under SU(6) → SU(5) × U(1):

- 24 ⊂ 35 (adjoint of SU(6))
- 15 ⊂ 21 = Sym²(6) (symmetric diquarks)
- 15̄ ⊂ 21̄

So the natural SU(6) composite content is **35 + 21 + 21̄ = 77 states**:

    35 → 24(0) + 1(0) + 5(+6) + 5̄(−6)
    21 → 15(+2) + 5(−4) + 1(−10)
    21̄ → 15̄(−2) + 5̄(+4) + 1(+10)

Total: 24 + 15 + 15̄ + 3×1 + 2×5 + 2×5̄ = 54 + 23 extra

The 23 extra states from the SU(6) → SU(5) breaking are:
- **3 singlets**: U(1) generator, tt̄ diquark scalar, conjugate
- **2×5**: tq̄ mesonic composites (top with each light quark)
- **2×5̄**: conjugates

**Physical interpretation**: At Stage 0 (Unfolding), all 6 flavors are active with SU(6) symmetry. The top = "elephant" is the 6th flavor: 6 → 5(+1) + 1(−5). Integrating out the top at Stage 1 (EW breaking) reduces SU(6) → SU(5) and removes the 23 top-channel composites. The surviving 54 = 24 + 15 + 15̄ are the sBootstrap content.

**The 15 vs 10̄ tension reframed**: Under SU(6), symmetric diquarks are in Sym²(6) = 21, antisymmetric diquarks in ∧²(6) = 15. The Pauli-preferred antisymmetric diquarks (10 of SU(5)) live in the SU(6) 15, while the sBootstrap's symmetric diquarks (15 of SU(5)) live in the SU(6) 21. The SU(6) embedding separates them cleanly into different representations.

## Next Agent Tasks (Round 5)

Following Manhattan methodology — pure math tasks, no theory context:

1. **Bloom δ-rotation**: Compute (M₀, δ) for all known Koide triples. Verify seeds at δ=3π/4. Check v₀-doubling.
2. **SU(6) branching**: SU(6) → SU(5) × U(1) for fundamental, adjoint, symmetric, antisymmetric reps.
3. **Pseudo-modulus CW**: Retry one-loop Coleman-Weinberg computation.
4. **Complete superpotential classification**: Most general renormalizable W for the field content.
