# Assembly Round 12: Sp(2) Vacuum, CKM Analytic, Brainstorm

## Strategy

Round 12 (even round) addresses:
1. **Sp(2) vacuum verification** — Does the lepton sector proposal from 11D actually work?
2. **CKM analytic** — Derive Cabibbo angle from tachyonic condensate structure
3. **Brainstorm** — 10 actionable items for paper completion

## Agents

- **12A** (opus): Sp(2) vacuum verification — COMPLETE ✅
- **12B** (opus): CKM analytic derivation — RUNNING
- **12C** (sonnet): Brainstorm — COMPLETE ✅

---

## Agent Results

### 12A: Sp(2) Vacuum Verification ⭐⭐

**Major result**: The O'Raifeartaigh-Koide mechanism is verified for the s-confining SU(2) = Sp(1) theory.

Key findings:
1. **ISS fails for s-confining theories**: SU(2) with N_f = 3 is s-confining (N_f = N_c + 1). Adding singlets X cannot create rank-condition SUSY breaking because all mesons are singlets. Polynomial F-term equations always have solutions.

2. **O'Raifeartaigh + Pfaffian works**: The superpotential
   W = fX + m L^{12}L^{23} + gX(L^{12})² + (ε/Λ⁵) L^{12}L^{13}L^{23}
   has a metastable SUSY-breaking vacuum at the origin with F_X = f ≠ 0, and a distant SUSY vacuum at L^{13} ~ Λ⁵/ε.

3. **Fermion spectrum**: At the metastable vacuum with t = gx₀/m:
   - Spectrum = (0 [Goldstino], 0 [flat direction], m₋, m₊)
   - m₊ = m(t + √(t²+1)), m₋ = m(√(t²+1) - t)
   - m₊ · m₋ = m² (exact)

4. **Koide Q = 2/3 exactly at t = √3**:
   Q = √(t²+1) / (√(t²+1) + 1) → Q = 2/3 ⟺ t = √3
   Seed spectrum: (0, (2-√3)m, (2+√3)m) — **verified to machine epsilon**

5. **Stability**: Vacuum is tachyon-free for gf/m² < 1/2. STr[M²] = 0 (exact).

6. **Pseudo-modulus**: CW potential gives t_min ~ 0.3-0.5 (NOT √3). Kähler pole K_X = |X|² - |X|⁴/(12μ²) pins ⟨X⟩ = √3·μ exactly, with μ = m/g.

7. **Universality**: Q = 2/3 is a property of the O'Raifeartaigh structure W = fX + m φφ̃ + gXφ², independent of gauge group. The dynamical superpotential (Pfaffian) provides metastability, not mass structure.

### 12C: Brainstorm ⭐

Priority-ordered actionable items:
1. **Baryon stabilization** (easy, HIGH): Add W_B = m_B BB̃ to lift baryonic flat direction. Makes seesaw the true vacuum.
2. **FCNC resolution** (medium, HIGH): Explicitly identify M^i_j = SM mesons, compute NMSSM-loop FCNC.
3. **CKM analytic** (medium, HIGH): Perturbative derivation of Cabibbo angle from tachyonic condensate.
4. **Vacuum tunneling** (medium, HIGH): Bounce action from seesaw to M=0 minimum.
5. **Yukawa flavor assignment** (easy, medium-high): Derive which meson couples to which Yukawa.
6. **Notation fixes** (easy, medium): Q overloading, v₀² rounding.
7. **Statistical framing** (easy-medium, medium): Prediction vs discovery distinction.
8. **μ-problem** (hard, HIGH): Accept as open, be explicit.
9. **Lepton sector** (hard/easy, medium): Sp(2) or demote to input.
10. **V_soft mediation** (medium, medium): Identify mediation mechanism or state as assumption.

**Single most important calculation**: Perturbative CKM from seesaw vacuum with B = B̃ = 0 enforced (combines Items 1 + 3).

### 12B: CKM Analytic Derivation ⭐⭐

**The B-term mechanism**: The NMSSM coupling generates F_X = λv²/2 = 2.18 × 10¹⁰ MeV⁵ at the seesaw vacuum. The three-point vertex W_{X, M^a_b, M^b_a} = M_k (where k is the third index) contracted with ⟨F_X⟩ gives holomorphic B-terms:

m²_eff(ab) = 2f_π² + 2|X₀|²M_k² − 2M_k F_X

The B-term dominates by ~10¹² over f_π²:
- |m²_ds| = 1.78 × 10¹⁶ MeV² (M_U = 408471)
- |m²_us| = 8.25 × 10¹⁵ MeV² (M_D = 188929)
- |m²_ud| = 4.12 × 10¹⁴ MeV² (M_S = 9446)

**Off-diagonal VEVs**: ε_ab = √(F_X/M_k) = √(λv²m_k/(2C))
- ε_ds = 231 MeV, ε_us = 340 MeV, ε_ud = 1520 MeV

**Key result — the Oakes ratio is EXACT**:
ε_us/ε_ud = √(M_S/M_D) = **√(m_d/m_s) = 0.22361 exactly** = tan θ_C (Oakes)

The Weinberg-Oakes relation is algebraically encoded in the tachyon hierarchy.

**But**: The absolute condensate is 183× too small for the Cabibbo angle. The physical normalization comes from PCAC/GOR (f_π² m_K² ~ m_s ⟨q̄q⟩), not from the tachyon minimum. The tachyon gives the RATIO (which IS the Cabibbo angle), not the MAGNITUDE.

---

## Synthesis

### Major advance in Round 12:
- ✅ **Sp(2) lepton mechanism verified** — O'Raifeartaigh + Pfaffian produces exact Q = 2/3 seed, with metastability from Pfaffian and pseudo-modulus pinned by Kähler pole
- ✅ **Universality established** — Same O'Raifeartaigh structure works for quarks (SU(3) SQCD) and leptons (SU(2)/Sp(1)), with Q = 2/3 as a universal consequence of t = √3

### Round 13 priorities (odd round = opus only):
1. **Baryon stabilization + CKM** — Add W_B, derive Cabibbo angle analytically
2. **Vacuum tunneling estimate** — Bounce action at seesaw
3. **Yukawa flavor assignment** — Half-page derivation
4. **Brainstorm** — Review remaining gaps
